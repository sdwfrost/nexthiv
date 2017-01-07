import os
import tempfile
import subprocess
import shutil
import random

from nexthiv.align import get_alignment

from Bio import AlignIO, SeqIO

def make_fasttree(cfg,alnFile,logfile,outfile):
    cmd=cfg["programs"]["fasttree"]["cmd"]
    model=cfg["programs"]["fasttree"]["model"]
    threads=cfg["programs"]["fasttree"]["threads"] # need to include as environment variables
    with open(logfile, 'w') as lh:
        cl = [cmd, '-gtr', '-gamma', '-nt', '-out', outfile, alnFile]
        try:
            subprocess.check_call(cl,stdout=lh,stderr=lh)
        except subprocess.CalledProcessError as e:
            print(e)

def make_iqtree(cfg,alnFile,logfile):
    cmd=cfg["programs"]["iqtree"]["cmd"]
    model=cfg["programs"]["iqtree"]["model"]
    threads=cfg["programs"]["iqtree"]["threads"]
    bootstrap=cfg["programs"]["iqtree"]["bootstrap_samples"]
    with open(logfile, 'w') as lh:
        cl = [cmd, '-s', alnFile, '-m', model, '-nt', str(threads), '-alrt', str(bootstrap),'-bb',str(bootstrap)]
        try:
            subprocess.check_call(cl,stdout=lh,stderr=lh)
        except subprocess.CalledProcessError as e:
            print(e)

def phylogeny(cfg,baseline=False,subsample=0):
    tmp_path = tempfile.mkdtemp(prefix='nexthiv-')
    if subsample<0:
        stop("Sample size must be positive")
    try:
        basename = "nexthiv"
        INPUT_FASTA=os.path.join(tmp_path, basename+'.fas')
        LOGFILE=os.path.join(tmp_path, basename+'.log')
        msa=get_alignment(cfg,baseline) # update to use masked alignment
        if subsample==0:
            AlignIO.write(msa,INPUT_FASTA,format="fasta")
        else:
            msal = len(msa)
            idx = random.sample(range(msal),subsample)
            idx.sort()
            submsa=[msa[i,:] for i in idx]
            SeqIO.write(submsa,INPUT_FASTA,format="fasta")
        if cfg["phylogeny"]["program"]=="iqtree":
            make_iqtree(cfg,INPUT_FASTA,LOGFILE)
        if cfg["phylogeny"]["program"]=="fasttree":
            OUTFILE=os.path.join(tmp_path, basename+'.nwk')
            make_fasttree(cfg,INPUT_FASTA,LOGFILE,OUTFILE)
    finally:
        print("Phylogeny finished")
        # shutil.rmtree(tmp_path)
