import nexthiv
from nexthiv.refseqs import HXB2_POL
from nexthiv.utils import get_baseline_ids

from copy import copy, deepcopy

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

from BioExt.args import (
    add_alphabet,
    add_reference,
    add_scorematrix
    )
from BioExt.io import BamIO
from BioExt.misc import compute_cigar, gapless
from BioExt.scorematrices import (
    DNAScoreMatrix,
    FrequenciesError,
    ProteinScoreMatrix,
    HIV_BETWEEN_F
    )
from BioExt.uds import _align_par
from BioExt.misc import compute_cigar, gapful, gapless

import pandas as pd

def dict2seqrec(cfg,s):
    SEQUENCE=cfg["sequence"]["column"]
    seqrec=[]
    for row in s:
        rec=SeqRecord(Seq(row[SEQUENCE],IUPAC.ambiguous_dna),
                      id=row['id'],
                      name=row['id'],
                      description=row['id'])
        seqrec.append(rec)
    return(seqrec)

def dict2aln(cfg,s):
    SEQUENCE_ALIGNED=cfg['sequence']['processed_name']
    seqrec=[]
    for row in s:
        if row[SEQUENCE_ALIGNED].count("-")<len(row[SEQUENCE_ALIGNED]):
            rec=SeqRecord(Seq(row[SEQUENCE_ALIGNED],Gapped(IUPAC.ambiguous_dna)),
                      id=row['id'],
                      name=row['id'],
                      description=row['id'])
            seqrec.append(rec)
    msa=MultipleSeqAlignment(seqrec)
    return(msa)

def align_seqs(cfg,reference=HXB2_POL,do_codon=True,score_matrix=HIV_BETWEEN_F.load(),reverse_complement=False,expected_identity=None,quiet=True):
    db=nexthiv.db.db_setup(name="rethinkdb")
    seqdict=db.get_seqs(cfg)
    records=dict2seqrec(cfg,seqdict)
    # Set up discards
    discards=[]
    def discard(record):
        discards.append(record)
    # Set up output alignment
    alignment = MultipleSeqAlignment([])
    alignment_length = len(reference)
    def suffix_pad(record):
        deficit = alignment_length - len(record)
        if deficit > 0:
            return SeqRecord(
               Seq(''.join((str(record.seq), '-' * deficit)), record.seq.alphabet),
               id=record.id,
               name=record.name,
               dbxrefs=copy(record.dbxrefs),
               description=record.description,
               annotations=copy(record.annotations),
               )
        return record
    def output(records):
        for record in records:
            alignment.append(suffix_pad(gapful(gapless(record), insertions=False)))
    # Align in parallel
    _align_par(
            reference,
            records,
            score_matrix,
            do_codon,
            reverse_complement,
            expected_identity,
            discard,output,quiet
            )
    # Return
    return(alignment)

def trim_alignment(cfg,a):
    startpos=cfg["sequence"]["startpos"]
    endpos=cfg["sequence"]["endpos"]
    trimal = deepcopy(a)
    for row in trimal:
        row.seq=row.seq[startpos:endpos]
    return(trimal)

def get_alignment(cfg,baseline=False):
    db=nexthiv.db.db_setup(name=cfg["db"]["backend"])
    if baseline:
        id_list=get_baseline_ids(cfg)
        s=db.get_aligned_seqs(cfg,id_list)
    else:
        s=db.get_aligned_seqs(cfg)
    msa=dict2aln(cfg,s)
    return(msa)
