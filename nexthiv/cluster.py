import os
import tempfile
import subprocess

from nexthiv.db import get_connection
from nexthiv.align import get_alignment

import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError

from Bio import AlignIO
from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment

import pandas as pd

def cluster(cfg):
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    tn93=cfg["programs"]["tn93"]
    tmp_path = tempfile.mkdtemp(prefix='nexthiv-')
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
    msa=get_alignment(cfg)
    AlignIO.write(msa,OUTPUT_FASTA_FN,format="fasta")
    with open(JSON_TN93_FN, 'w') as tn93_fh:
        tn93_process = [TN93DIST, '-q', '-o', OUTPUT_TN93_FN, '-t',
                    THRESHOLD, '-a', AMBIGUITIES, '-l',
                    MIN_OVERLAP, '-g', FRACTION if AMBIGUITIES == 'resolve' else '1.0',
                    '-f', output_format, OUTPUT_FASTA_FN]
        subprocess.check_call(tn93_process,stdout=tn93_fh,stderr=tn93_fh)
    dst=pd.read_csv(OUTPUT_TN93_FN)
    return(dst)

def insert_clustering(cfg):
    # Retrieve sequence names from sequences table
    NEXTHIV_DB=cfg["rethinkdb"]["db"]
    connection = get_connection(cfg)
    if not r.db(NEXTHIV_DB).table_list().contains(cfg):
        stop("Table "+tbl+" does not exist.")
    dst=cluster(cfg)
    dst=dst.sort_values(by=["ID1","Distance"],ascending=True)
    mindst=dst.drop_duplicates(subset="ID1",keep="first")
    try:
        tbl=cfg["clustering"]["table"]
        for row in mindst.iterrows():
            r.db(NEXTHIV_DB).table(tbl).insert({"id":row[1]["ID1"], "ALTER":row[1]["ID2"], "MINDST":row[1]["Distance"]}).run(connection)
    except RqlRuntimeError:
        print('Error!')
    try:
        tbl=cfg["clustering"]["distances_table"]
        for row in dst.iterrows():
            r.db(NEXTHIV_DB).table(tbl).insert({"id":row[1]["ID1"], "ALTER":row[1]["ID2"], "DST":row[1]["Distance"]}).run(connection)
    except RqlRuntimeError:
        print('Error!')
    finally:
        connection.close()
