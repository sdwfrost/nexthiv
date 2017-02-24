import os
import inspect
import numpy as np
import nexthiv
import pysam
from copy import copy
from BioExt.io._SamBamIO import _to_seqrecord
from BioExt.misc import gapful

def get_data_directory():
    dd=os.path.join(os.path.dirname(inspect.getfile(nexthiv)),"data")
    return(dd)

def get_seqids(cfg):
    db=nexthiv.db.db_setup(name=cfg["db"]["backend"])
    sq=cfg["sequence"]
    tbl=sq["table"]
    seqname="id" # update
    sids=db.get_ids(cfg,tbl,idcol=seqname)
    return(sids)

def get_baseline_ids(cfg):
    db=nexthiv.db.db_setup(name=cfg["db"]["backend"])
    sq=cfg["sequence"]
    tbl=sq["table"]
    seqname="id" # update
    pid=sq["pid"]
    ordering=sq["ordering"]
    df=db.get_dataframe(cfg,tbl,[seqname,pid,ordering])
    df=df.sort_values(by=[pid,ordering],ascending=True)
    df=df.drop_duplicates(subset=pid,keep="first")
    id_list=df["id"].tolist()
    return(id_list)

def insert_baseline_into_table(cfg,tbl,name="BASELINE"):
    sids=get_seqids(cfg)
    bids=set(get_baseline_ids(cfg))
    data=[{"id":id ,name: id in bids} for id in sids]
    db=nexthiv.db.db_setup(name=cfg["db"]["backend"])
    db.db_update_by_id(cfg,tbl,data)

def bam2records(bam_file, start=None, end=None):
    samfile = None
    try:
        # Index bam file in order to samfile.fetch
        #pysam.index(bam_file)
        samfile = pysam.Samfile(bam_file, 'rb')
        length = samfile.header['SQ'][0]['LN']
        fetch_args = []
        if start is not None or end is not None:
            fetch_args.extend([samfile.getrname(0), start, end])
        seqs = [gapful(_to_seqrecord(samfile, record), insertions=False) for record in samfile.fetch(*fetch_args)]
        result = [(seq + ('-' * (length - len(seq))))[start:end] for seq in seqs]
    finally:
        if samfile is not None:
            samfile.close()
    return(result)

def extract_annotations(s,column,sep):
    x=copy(s)
    x=x.split(sep)
    return(x[column])

def invert_dict(d):
    inv_d = {}
    for k, v in d.items():
        inv_d[v] = inv_d.get(v, [])
        inv_d[v].append(k)
    return(inv_d)

def np_to_phylip(nm,m,fn):
    f=open(fn,'w')
    n=len(nm)
    f.write(str(n)+"\n")
    for i in range(n):
        f.write(nm[i])
        for j in range(n):
            f.write("\t"+str(m[i,j]))
        f.write("\n")
    f.close()

def get_list(cfg,tbl,col):
    db=nexthiv.db.db_setup(name=cfg["db"]["backend"])
    l=db.get_ids(cfg,tbl,col)
    return(l)

def get_set(cfg,tbl,col):
    l=get_list(cfg,tbl,col)
    s=set(l)
    return(s)
