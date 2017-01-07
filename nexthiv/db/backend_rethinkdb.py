import os
import csv
import json
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError
import pandas as pd

# Register TSV as a dialect
csv.register_dialect('tsv',delimiter='\t',quoting=csv.QUOTE_NONE)

def get_connection(cfg):
    HOST=cfg["db"]["host"]
    PORT=cfg["db"]["port"]
    NEXTHIV_DB=cfg["db"]["name"]
    connection = r.connect(host=HOST, port=PORT)
    return(connection)

def db_init(cfg):
    NEXTHIV_DB=cfg["db"]["name"]
    connection = get_connection(cfg)
    tbls=cfg["tables"]
    try:
        r.db_create(NEXTHIV_DB).run(connection)
    except RqlRuntimeError:
        stop('nexthiv database already exists.')
    for tbl in tbls:
        try:
            r.db(NEXTHIV_DB).table_create(tbl).run(connection)
        except RqlRuntimeError:
            stop(tbl+" table already exists.")
    connection.close()

def tsv2json(f):
    csv_reader=csv.DictReader(f,dialect='tsv')
    data=json.dumps([row for row in csv_reader])
    return(json.loads(data))

def db_insert(cfg,tbl,data,chunksize=100000):
    NEXTHIV_DB=cfg["db"]["name"]
    connection=get_connection(cfg)
    chunked_data=[data[i:i+chunksize] for i in range(0,len(data),chunksize)]
    if not r.db(NEXTHIV_DB).table_list().contains(tbl):
        stop("Table "+tbl+" does not exist.")
    try:
        for d in chunked_data:
            r.db(NEXTHIV_DB).table(tbl).insert(d).run(connection)
    except RqlRuntimeError:
        stop("Error importing "+tbl+".")
    finally:
        connection.close()

def db_update(cfg,tbl,data,chunksize=100000):
    NEXTHIV_DB=cfg["db"]["name"]
    connection=get_connection(cfg)
    chunked_data=[data[i:i+chunksize] for i in range(0,len(data),chunksize)]
    if not r.db(NEXTHIV_DB).table_list().contains(tbl):
        stop("Table "+tbl+" does not exist.")
    try:
        for d in chunked_data:
            r.db(NEXTHIV_DB).table(tbl).update(d).run(connection)
    except RqlRuntimeError:
        stop("Error importing "+tbl+".")
    finally:
        connection.close()

def db_update_by_id(cfg,tbl,data):
    NEXTHIV_DB=cfg["db"]["name"]
    connection=get_connection(cfg)
    if not r.db(NEXTHIV_DB).table_list().contains(tbl):
        stop("Table "+tbl+" does not exist.")
    try:
        for d in data:
            id=d["id"]
            d2={k: v for k, v in d.items() if k!="id"}
            r.db(NEXTHIV_DB).table(tbl).get_all(id).update(d2).run(connection)
    except RqlRuntimeError:
        stop("Error importing "+tbl+".")
    finally:
        connection.close()

def get_ids(cfg,tbl,idcol="id"):
    NEXTHIV_DB=cfg["db"]["name"]
    connection = get_connection(cfg)
    try:
        curs=r.db(NEXTHIV_DB).table(tbl).pluck(idcol).run(connection)
    except RqlRuntimeError:
        print('Error!')
    finally:
        ids=[x[idcol] for x in curs]
        connection.close()
    return(ids)

def get_seqs(cfg):
    NEXTHIV_DB=cfg["db"]["name"]
    SEQUENCE=cfg["sequence"]["column"]
    TBL=cfg["sequence"]["table"]
    connection = get_connection(cfg)
    try:
        curs=r.db(NEXTHIV_DB).table(TBL).pluck("id",SEQUENCE).run(connection)
        # print('Sequence extraction completed.')
    except RqlRuntimeError:
        print('Error!')
    finally:
        s=[x for x in curs]
        connection.close()
    return(s)

def get_aligned_seqs(cfg,ids=None):
    NEXTHIV_DB=cfg["db"]["name"]
    SEQUENCE=cfg["sequence"]["processed_name"]
    TBL=cfg["sequence"]["processed_table"]
    connection = get_connection(cfg)
    try:
        if ids is None:
            curs=r.db(NEXTHIV_DB).table(TBL).pluck("id",SEQUENCE).run(connection)
            s=[x for x in curs]
        else:
            s=[]
            for id in ids:
                curs=r.db(NEXTHIV_DB).table(TBL).get_all(id).pluck("id",SEQUENCE).run(connection)
                s.append(list(curs)[0])
        # print('Sequence extraction completed.')
    except RqlRuntimeError:
        print('Error!')
    finally:

        connection.close()
    return(s)

def insert_seqs(cfg,records,tbl,colname):
    NEXTHIV_DB=cfg["db"]["name"]
    connection = get_connection(cfg)
    if not r.db(NEXTHIV_DB).table_list().contains(tbl):
        stop("Table "+tbl+" does not exist.")
    try:
        for record in records:
            seqname=record.name
            d={'id':seqname,cfg['sequence']['processed_name']:str(record.seq)}
            r.db(NEXTHIV_DB).table(tbl).insert(d).run(connection)
            #r.db(NEXTHIV_DB).table(tbl).filter({"id":seqname}).update({colname:str(record.seq)}).run(connection)
        print('Inserted sequences into '+tbl+' table completed.')
    except RqlRuntimeError:
        print('Error!')
    finally:
        connection.close()

def get_dataframe(cfg,tbl,cols=None):
    NEXTHIV_DB=cfg["db"]["name"]
    connection = get_connection(cfg)
    try:
        if cols is None:
            curs=r.db(NEXTHIV_DB).table(tbl).run(connection)
        else:
            curs=r.db(NEXTHIV_DB).table(tbl).pluck(cols).run(connection)
    except RqlRuntimeError:
        print('Error!')
    finally:
        s=[x for x in curs]
        connection.close()
    df=pd.DataFrame.from_records(s)
    return(df)

def get_dict(cfg,tbl,ky,val):
    NEXTHIV_DB=cfg["db"]["name"]
    connection = get_connection(cfg)
    try:
        curs=r.db(NEXTHIV_DB).table(tbl).pluck(ky,val).run(connection)
    except RqlRuntimeError:
        print('Error!')
    finally:
        s=[x for x in curs]
        connection.close()
    d={x[ky]:x[val] for x in s}
    return(d)

def table_exists(cfg,tbl):
    NEXTHIV_DB=cfg["db"]["name"]
    return(r.db(NEXTHIV_DB).table_list().contains(tbl))
