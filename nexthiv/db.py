import os
import csv
import json
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError
import pandas as pd

# Register TSV as a dialect
csv.register_dialect('tsv',delimiter='\t',quoting=csv.QUOTE_NONE)

def get_connection(cfg):
    RDB_HOST=cfg["rethinkdb"]["host"]
    RDB_PORT=cfg["rethinkdb"]["port"]
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    connection = r.connect(host=RDB_HOST, port=RDB_PORT)
    return(connection)

def db_setup(cfg):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
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
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
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

def get_seqs(cfg):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    SEQUENCE=cfg["sequence"]["name"]
    connection = get_connection(cfg)
    try:
        curs=r.db(NEXTHIV_DB).table('sequences').pluck("id",SEQUENCE).run(connection)
        # print('Sequence extraction completed.')
    except RqlRuntimeError:
        print('Error!')
    finally:
        s=[x for x in curs]
        connection.close()
    return(s)

def insert_seqs(cfg,records,tbl,colname):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
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

def get_dataframe(cfg,tbl):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    connection = get_connection(cfg)
    try:
        curs=r.db(NEXTHIV_DB).table(tbl).run(connection)
    except RqlRuntimeError:
        print('Error!')
    finally:
        s=[x for x in curs]
        connection.close()
    df=pd.DataFrame.from_records(s)
    return(df)

def get_dict(cfg,tbl,ky,val):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
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
