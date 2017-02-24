import os
import tempfile
import subprocess
import shutil
from collections import Counter

import nexthiv
from nexthiv.align import get_alignment
from nexthiv.utils import get_data_directory

from Bio import AlignIO
from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment

import numpy as np
import pandas as pd
import networkx as nx

from BioExt.tn93 import tn93

def tn93_matrix(query,matchMode,minOverlap):
    nnodes=len(query)
    names=[q.id for q in query]
    m = np.zeros((nnodes,nnodes),dtype=np.float64)
    L = len(str(query[0].seq))
    for i in range(nnodes-1):
        q = query[i]
        for j in range(i+1,nnodes):
            r = query[j]
            d = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if d==-1:
                d=np.nan
            m[i,j]=d
            m[j,i]=d
    return(names,m)

def tn93_distance(query,matchMode,minOverlap,threshold=1.0,names=["ID1","ID2","Distance"]):
    nnodes=len(query)
    L = len(str(query[0].seq))
    result=[]
    for i in range(nnodes-1):
        q = query[i]
        for j in range(i+1,nnodes):
            r = query[j]
            d = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if d<=threshold:
                result.append({names[0]:q.id,names[1]:r.id,names[2]:d})
    return(result)

def tn93_distance_ref(query,ref,matchMode,minOverlap,threshold=0.05,names=["ID1","ID2","Distance"]):
    L = len(str(query[0].seq))
    result=[]
    for i in range(len(query)):
        q = query[i]
        for j in range(len(ref)):
            r = ref[j]
            d = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if d<=threshold:
                result.append({names[0]:q.id,names[1]:r.id,names[2]:d})
    return(result)

def tn93_closest(query,matchMode,minOverlap,piddict=None,names=["id","ALTER","MINDST","PID"]):
    nnodes=len(query)
    mind={q.id:[None,1.0] for q in query}
    L = len(str(query[0].seq))
    for i in range(nnodes-1):
        q = query[i]
        for j in range(i+1,nnodes):
            r = query[j]
            if piddict==None or piddict[q.id]!=piddict[r.id]:
                newd = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
                if newd < mind[q.id][1]:
                    mind[q.id] = [r.id,newd]
                    if piddict!=None:
                        mind[q.id].append(piddict[r.id])
                if newd < mind[r.id][1]:
                    mind[r.id] = [q.id,newd]
                    if piddict!=None:
                        mind[r.id].append(piddict[q.id])
    if piddict==None:
        result = [{names[0]:q.id,names[1]:mind[q.id][0],names[2]:mind[q.id][1]} for q in query]
    else:
        result = [{names[0]:q.id,names[1]:mind[q.id][0],names[2]:mind[q.id][1],names[3]:mind[q.id][2]} for q in query]
    return(result)

def tn93_closest_ref(query,refs,matchMode,minOverlap,names=["id","REFALTER","REFMINDST"]):
    result=[]
    L = len(str(query[0].seq))
    for q in query:
        d = 1.0
        rid = None
        for r in refs:
            newd = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if newd < d:
                d = newd
                rid=r.id
        result.append({names[0]:q.id,names[1]:rid,names[2]:d})
    return(result)

def cluster(cfg,baseline=False,names=["ID1","ID2","Distance"]):
    NEXTHIV_DB=cfg["db"]["name"]
    tn93=cfg["clustering"]
    AMBIGUITIES=int(tn93["ambiguities"])
    MIN_OVERLAP=int(tn93["min_overlap"])
    THRESHOLD=int(tn93["threshold"])
    msa=get_alignment(cfg,baseline)
    dst=tn93_distance(msa,AMBIGUITIES,MIN_OVERLAP,THRESHOLD,names)
    # Convert to dataframe
    return(dst)

# need to update; not using BAM right now, also need to use BioExt
def cluster_refs(cfg,baseline=False):
    NEXTHIV_DB=cfg["db"]["name"]
    tn93=cfg["programs"]["tn93"]
    tmp_path = tempfile.mkdtemp(prefix='nexthiv-')
    msa=get_alignment(cfg,baseline)
    # Write reference to temporary directory
    dd=get_data_directory()
    bamfile=os.path.join(dd,"hiv_refs_prrt_trim.bam")
    # Convert into FASTA
    try:
        for a in msa:
            basename = "nexthiv"
            OUTPUT_FASTA_FN  = os.path.join(tmp_path, basename+'.fas')
            JSON_TN93_FN  = os.path.join(tmp_path, basename+'_user.tn93output.json')
            TN93DIST=tn93["cmd"]
            OUTPUT_TN93_FN = os.path.join(tmp_path, basename+'_user.tn93output.csv')
            THRESHOLD=str(tn93["threshold"])
            AMBIGUITIES=tn93["ambiguities"]
            MIN_OVERLAP=str(tn93["min_overlap"])
            FRACTION=str(tn93["fraction"])
            output_format="csv"
            AlignIO.write(a,OUTPUT_FASTA_FN,format="fasta")
            with open(JSON_TN93_FN, 'w') as tn93_fh:
                tn93_process = [TN93DIST, '-q', '-o', OUTPUT_TN93_FN, '-t',
                    THRESHOLD, '-a', AMBIGUITIES, '-l',
                    MIN_OVERLAP, '-g', FRACTION if AMBIGUITIES == 'resolve' else '1.0',
                    '-f', output_format, OUTPUT_FASTA_FN]
                subprocess.check_call(tn93_process,stdout=tn93_fh,stderr=tn93_fh)
            dst=pd.read_csv(OUTPUT_TN93_FN)
            # Extract minimum from dst
    finally:
        shutil.rmtree(tmp_path)
    return(dst)

def insert_distances(cfg,baseline=False):
    # Retrieve sequence names from sequences table
    NEXTHIV_DB=cfg["db"]["name"]
    tbl=cfg["clustering"]["distances_table"]
    db=nexthiv.db.db_setup(name="rethinkdb")
    if not db.table_exists(cfg,tbl):
        stop("Table "+tbl+" does not exist.")
    data=cluster(cfg,baseline,names=["id","ALTER","DST"])
    db.db_insert(cfg,tbl,data)

def insert_clustering(cfg,baseline=False):
    # Retrieve sequence names from sequences table
    NEXTHIV_DB=cfg["db"]["name"]
    tbl=cfg["clustering"]["table"]
    db=nexthiv.db.db_setup(name="rethinkdb")
    if not db.table_exists(cfg,tbl):
        stop("Table "+tbl+" does not exist.")
    dst=cluster(cfg,baseline)
    iddict=db.get_dict(cfg,cfg["sequence"]["table"],"id",cfg["sequence"]["pid"])
    pid1=[iddict[x] for x in dst["ID1"]]
    pid2=[iddict[x] for x in dst["ID2"]]
    idx=[x[0]!=x[1] for x in zip(pid1,pid2)]
    dst=dst[idx]
    dst=dst.sort_values(by=["ID1","Distance"],ascending=True)
    mindst=dst.drop_duplicates(subset="ID1",keep="first")
    data=[]
    for row in mindst.iterrows():
        data.append({"id":row[1]["ID1"], "ALTER":row[1]["ID2"], "MINDST":row[1]["Distance"]})
    db.db_insert(cfg,tbl,data)
