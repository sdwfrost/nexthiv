import os
import tempfile
import subprocess
import shutil

import nexthiv
from nexthiv.align import get_alignment

from Bio import AlignIO
from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment

import pandas as pd

def cluster(cfg,baseline=False):
    NEXTHIV_DB=cfg["db"]["name"]
    tn93=cfg["programs"]["tn93"]
    tmp_path = tempfile.mkdtemp(prefix='nexthiv-')
    try:
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
        msa=get_alignment(cfg,baseline)
        AlignIO.write(msa,OUTPUT_FASTA_FN,format="fasta")
        with open(JSON_TN93_FN, 'w') as tn93_fh:
            tn93_process = [TN93DIST, '-q', '-o', OUTPUT_TN93_FN, '-t',
                    THRESHOLD, '-a', AMBIGUITIES, '-l',
                    MIN_OVERLAP, '-g', FRACTION if AMBIGUITIES == 'resolve' else '1.0',
                    '-f', output_format, OUTPUT_FASTA_FN]
            subprocess.check_call(tn93_process,stdout=tn93_fh,stderr=tn93_fh)
        dst=pd.read_csv(OUTPUT_TN93_FN)
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
    dst=cluster(cfg,baseline)
    dst=dst.sort_values(by=["ID1","Distance"],ascending=True)
    data=[]
    for row in dst.iterrows():
        data.append({"id":row[1]["ID1"], "ALTER":row[1]["ID2"], "DST":row[1]["Distance"]})
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
