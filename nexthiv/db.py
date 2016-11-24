import os
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError

def dbSetup(cfg):
    RDB_HOST=cfg["rethinkdb"]["host"]
    RDB_PORT=cfg["rethinkdb"]["port"]
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    connection = r.connect(host=RDB_HOST, port=RDB_PORT)
    tbls=cfg["tables"]
    try:
        r.db_create(NEXTHIV_DB).run(connection)
    except RqlRuntimeError:
        print('nexthiv database already exists.')
    for tbl in tbls:
        try:
            r.db(NEXTHIV_DB).table_create(tbl).run(connection)
        except RqlRunTimeError:
    connection.close()
