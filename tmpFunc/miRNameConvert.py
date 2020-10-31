#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/7/14
Description     :

"""
import sys
import pdb
import BioinfoComm, Comm
from Comm import print0

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

def convertJiangData():
    # input
    inputFn = r'C:\Users\Administrator\Desktop\mirna_regulatory_score_all.txt'
    miRBaseFn = r'D:\Dropbox\project\MotifFinding\sources\exo_miRNA_raw_data\miRBase_pre_mature_IDs_conversion.txt'

    # output
    outputFn = r'C:\Users\Administrator\Desktop\mirna_regulatory_score_all_formatted.txt'

    miRNANm2Seq, prevNm2seq, prevNm2matNm = loadMiRNABase(miRBaseFn)
    with open(inputFn) as inputFileobj, open(outputFn, 'w') as outputFileobj:
        titleLine = inputFileobj.readline()
        outputFileobj.write('%s' % titleLine)
        for line in inputFileobj:
            line = line.strip()
            items = line.split('\t')
            gene = items[0]
            prevMiRNm = 'hsa-%s' % items[1]
            scores = items[2:]
            # matMiRNm = prevNm2matNm.get(prevMiRNm, prevMiRNm)
            matMiRNm = formatMiRName(prevMiRNm, prevNm2matNm)
            outputLine = '%s\t%s\t%s' % (gene, matMiRNm, '\t'.join(scores))
            outputFileobj.write('%s\n' % outputLine)
#end_func

def formatMiRName(prevMiRNm, prevNm2matNm):
    prevMiRNm = prevMiRNm.upper()
    # if '*' in prevMiRNm:
        # print 1
    matNm = prevNm2matNm.get(prevMiRNm, prevMiRNm)
    matNm = matNm.lower()
    matNm = matNm.replace('mir', 'miR')
    return matNm
#end_func

def loadMiRNABase(filenm):
    miRNANm2Seq = {} # mature_name 1&2
    prevNm2seq = {} # mature_previous_IDs 1&2
    prevNm2matNm = {}

    # load sequences - map from miRNA name to sequences
    with open(filenm) as f:
        f.readline()
        for lineIdx, line in enumerate(f):
            items = line.strip().split('\t')
            mat1Nm = items[3].upper()
            mat1PrevNm = items[5].upper()
            seq1 = items[6].upper()
            miRNANm2Seq[mat1Nm] = seq1
            prevNmLis = mat1PrevNm.split(';')
            for prevNm in prevNmLis:
                prevNm2seq[prevNm] = seq1
                prevNm2matNm[prevNm] = mat1Nm
            if len(items) == 7: continue
            mat2Nm = items[7].upper()
            mat2PrevNm = items[9].upper()
            seq2 = items[10].upper()
            miRNANm2Seq[mat2Nm] = seq2
            prevNmLis = mat2PrevNm.split(';')
            for prevNm in prevNmLis:
                prevNm2seq[prevNm] = seq2
                prevNm2matNm[prevNm] = mat1Nm

    return miRNANm2Seq, prevNm2seq, prevNm2matNm
#end_func

def test():
    convertJiangData()
#end_test

def main():
    if PRINT_INTO_LOG:
        sys.stdout = open(LOG_FN, 'w')
        test()
        sys.stdout.close()
    else:
        test()
#end_main


if __name__ == "__main__":
    if ENABLE_CPROFILE:
        import cProfile
        import pstats
        cp = cProfile.Profile()
        cp.runcall(main)
        ps = pstats.Stats(cp)
        sorted_stats = ps.sort_stats('time')
        sorted_stats.print_stats(10)
    else:
       main()
#end_if
