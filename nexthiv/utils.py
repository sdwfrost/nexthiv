import os
import inspect
import nexthiv

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
