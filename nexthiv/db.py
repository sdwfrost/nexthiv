import os
import csv
import json
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError

# Register TSV as a dialect
csv.register_dialect('tsv',delimiter='\t',quoting=csv.QUOTE_NONE)

def get_connection(cfg):
    RDB_HOST=cfg["rethinkdb"]["host"]
    RDB_PORT=cfg["rethinkdb"]["port"]
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    connection = r.connect(host=RDB_HOST, port=RDB_PORT)
    return(connection)

def dbSetup(cfg):
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

def dbInsert(cfg,tbl,data,chunksize=100000):
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
