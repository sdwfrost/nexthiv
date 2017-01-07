import os
import re
import itertools

import nexthiv
from nexthiv.align import get_alignment
from nexthiv.utils import get_data_directory

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement as rc
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Data.CodonTable import TranslationError

from BioExt.misc import translate_ambiguous

def get_positions(cfg):
    records=get_alignment(cfg)
    prlen=cfg["resistance"]["protease_length"]
    rtlen=cfg["resistance"]["rt_length"]
    prnames = ["PR"+str(x) for x in range(1,prlen+1)]
    rtnames = ["RT"+str(x) for x in range(1,rtlen+1)]
    positions=[]
    for row in records:
        pos={}
        pos["id"]=row.name # update
        idx=0
        names=prnames
        try:
            aa=translate_ambiguous(row,gap_char="-",trim_gaps=False)
        except:
            mutseq=row.seq.tomutable()
            for i in range(0,len(row),3):
                codon=row[i:i+3]
                codonstr=str(codon.seq)
                gapcount=codonstr.count('-')
                if gapcount==1 or gapcount==2:
                    for j in range(3):
                        if codonstr[j]=="-":
                            mutseq[i+j]="N"
            newseq=mutseq.toseq()
            aa=translate_ambiguous(newseq,gap_char="-",trim_gaps=False)
        for a in aa:
            if idx==prlen and names==prnames:
                idx=0
                names=rtnames
            pos[names[idx]]=a
            idx+=1
        positions.append(pos)
    return(positions)

# function to calculate PI resistance score
def get_resistance_PI(resI,sRule,cRule):
    # read in the score rule files
    piS = [line.strip() for line in open(sRule,'r')]
    piCS = [line.strip() for line in open(cRule,'r')]
    #len(piS)

    # create a list of drug scores intialised with 0
    #dScore = [0 for x in range(nD)]
    #print(dScore)
    pre = 'PR'
    # get the name of the drugs from header line of piS
    words = piS[0].split()
    drugs = words[2:]
    #print(drugs)

    # get positions with resiatance mutations
    mutPos = set()

    # create a dictionary of scores for each drug
    dScore = {}
    for i in range(len(drugs)):
        dScore[words[i+2]] = 0

    #print(dScore)

    # calculate based on single score
    for j in range(1,len(piS)):
        words = piS[j].split()
        pos = pre + words[0]
        aa = words[1]
        #print(pos,aa)
        if aa in resI[pos]:
            #print(pos,aa)
            mutPos.add(pos)
            scores = words[2:]
            for z in range(len(drugs)):
                dScore[drugs[z]] += int(scores[z])

    # calculate based on combination score
    for j in range(1,len(piCS)): # len(piCS)
        words = piCS[j].split()

        if len(words) == 0:
            continue

        mFlag = 0 # to check how many rules matched

        # split the rule into individual ones
        rules = words[1].split('+')
        #print(rules)

        # get positions for the combined mutations
        tPos = set()

        for r in rules:
            #print(r)

            m  = [x.start() for x in re.finditer('\d', r)]
            if m:
                pos = pre + r[m[0]:(m[-1]+1)]
                aas = list(r[(m[-1]+1):])
                #print(pos,aas)
                com = list(set(aas) & set(resI[pos]))
                #print(aas,res[0][pos])
                if com:
                    mFlag += 1
                    tPos.add(pos)

        if mFlag == len(rules):
            #print(j,pos,com)
            dScore[words[0]] += int(words[2])
            mutPos = mutPos | tPos

    return dScore, mutPos
#*************************************************************

# function to get NRTI scores
def get_resistance_NRTI(resI,sRule,cRule):
    # read in the score rule files
    piS = [line.strip() for line in open(sRule,'r')]
    piCS = [line.strip() for line in open(cRule,'r')]
    #len(piS)
    pre = 'RT'
    # create a list of drug scores intialised with 0
    #dScore = [0 for x in range(nD)]
    #print(dScore)

    # get the name of the drugs from header line of piS
    words = piS[0].split()
    drugs = words[2:]
    #print(drugs)

    # get positions with resiatance mutations
    mutPos = set()

    # create a dictionary of scores for each drug
    dScore = {}
    for i in range(len(drugs)):
        dScore[words[i+2]] = 0

    #print(dScore)

    # calculate based on single score
    for j in range(1,len(piS)):
        words = piS[j].split()
        pos = pre + words[0]
        aas = list(words[1])
        com = list(set(aas) & set(resI[pos]))
        #print(pos,aa)
        if com:
            #print(pos,com)
            mutPos.add(pos)
            scores = words[2:]
            for z in range(len(drugs)):
                dScore[drugs[z]] += int(scores[z])


    # calculate based on combination score
    for j in range(1,len(piCS)): # len(piCS)
        words = piCS[j].split()

        if len(words) == 0:
            continue

        mFlag = 0 # to check how many rules matched

        # split the rule into individual ones
        rules = words[1].split('+')
        #print(rules)

        # get positions for the combined mutations
        tPos = set()

        for r in rules:
            #print(r)

            m  = [x.start() for x in re.finditer('\d', r)]
            if m:
                pos = pre + r[m[0]:(m[-1]+1)]
                aas = list(r[(m[-1]+1):])
                #print(pos,aas)
                com = list(set(aas) & set(resI[pos]))
                #print(aas,res[0][pos])
                if com:
                    mFlag += 1
                    tPos.add(pos)

        if mFlag == len(rules):
            #print(j,pos,com)
            # get all the drugs name
            drs = words[0].split(',')
            for d in drs:
                dScore[d] += int(words[2])
            mutPos = mutPos | tPos

    return dScore, mutPos


# function to get NNRTI scores
def get_resistance_NNRTI(resI,sRule,cRule):
    # read in the score rule files
    piS = [line.strip() for line in open(sRule,'r')]
    piCS = [line.strip() for line in open(cRule,'r')]
    #len(piS)

    # create a list of drug scores intialised with 0
    #dScore = [0 for x in range(nD)]
    #print(dScore)
    pre = 'RT'
    # get the name of the drugs from header line of piS
    words = piS[0].split()
    drugs = words[2:]
    #print(drugs)

    # get positions with resiatance mutations
    mutPos = set()


    # create a dictionary of scores for each drug
    dScore = {}
    for i in range(len(drugs)):
        dScore[words[i+2]] = 0

    #print(dScore)

    # calculate based on single score
    for j in range(1,len(piS)):
        words = piS[j].split()
        pos = pre + words[0]
        aa = words[1]
        #print(pos,aa)
        if aa in resI[pos]:
            mutPos.add(pos)
            scores = words[2:]
            for z in range(len(drugs)):
                dScore[drugs[z]] += int(scores[z])

    # calculate based on combination score
    for j in range(1,len(piCS)): # len(piCS)
        words = piCS[j].split()

        if len(words) == 0:
            continue

        mFlag = 0 # to check how many rules matched

        # split the rule into individual ones
        rules = words[1].split('+')
        #print(rules)

        # get positions for combined mutations
        tPos = set()

        for r in rules:
            #print(r)

            m  = [x.start() for x in re.finditer('\d', r)]
            if m:
                pos = pre + r[m[0]:(m[-1]+1)]
                aas = list(r[(m[-1]+1):])
                #print(pos,aas)
                com = list(set(aas) & set(resI[pos]))
                #print(aas,res[0][pos])
                if com:
                    mFlag += 1
                    tPos.add(pos)

        if mFlag == len(rules):
            dScore[words[0]] += int(words[2])
            mutPos = mutPos | tPos

    return dScore, mutPos

# function to calculate PI resistance score
def get_sdrm(resI,sRule):
    # read in the score rule files
    piS = [line.strip() for line in open(sRule,'r')]

    # get SDRM positions
    mutPos = list()

    for j in range(1,len(piS)):
        words = piS[j].split('\t')
        pre = words[2]
        wt = words[3]
        pos = words[4]
        prepos = pre + pos
        aa = words[5]
        cls = words[6]
        if cls!="INI":
            if aa in resI[prepos]:
                mutPos.append([pre,pos,wt,aa,cls])

    d={}
    d["SDRM_COLUMN"]=[]
    d["PI_MUTATIONS"]=[]
    d["NRTI_MUTATIONS"]=[]
    d["NNRTI_MUTATIONS"]=[]
    d["PI_MUTATION_COUNT"]=0
    d["NRTI_MUTATION_COUNT"]=0
    d["NNRTI_MUTATION_COUNT"]=0
    for p in mutPos:
        d["SDRM_COLUMN"].append(p[0]+p[1])
        d[p[4]+"_MUTATIONS"].append(p[2]+p[1]+p[3])
        d[p[4]+"_MUTATION_COUNT"]+=1

    return d

def get_scores(cfg):
    pos=get_positions(cfg)
    dd=get_data_directory()
    data=[]
    for p in pos:
        pid={"id":p["id"]} # need to update
        pi_scores, pi_positions=get_resistance_PI(p,os.path.join(dd,"Scores_PI.txt"),os.path.join(dd,"combinationScores_PI.txt"))
        max_pi_score={"MAX_PI_SCORE":max(pi_scores.values())}
        #pi_positions={"PI_MUTATIONS":list(pi_positions)}
        nrti_scores, nrti_positions=get_resistance_NRTI(p,os.path.join(dd,"Scores_NRTI.txt"),os.path.join(dd,"combinationScores_NRTI.txt"))
        max_nrti_score={"MAX_NRTI_SCORE":max(nrti_scores.values())}
        #nrti_positions={"NRTI_MUTATIONS":list(nrti_positions)}
        nnrti_scores, nnrti_positions=get_resistance_NNRTI(p,os.path.join(dd,"Scores_NNRTI.txt"),os.path.join(dd,"combinationScores_NNRTI.txt"))
        max_nnrti_score={"MAX_NNRTI_SCORE":max(nnrti_scores.values())}
        #nnrti_positions={"NNRTI_MUTATIONS":list(nnrti_positions)}
        sdrm=get_sdrm(p,os.path.join(dd,"SDRM_2009.txt"))
        d=dict(itertools.chain(pid.items(),pi_scores.items(),nrti_scores.items(),nnrti_scores.items()))
        d=dict(itertools.chain(d.items(),max_pi_score.items(),max_nrti_score.items(),max_nnrti_score.items()))
        d=dict(itertools.chain(d.items(),sdrm.items()))
        data.append(d)
    return(data)
