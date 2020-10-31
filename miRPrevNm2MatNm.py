#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/8/16
Description     :

"""
import sys

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

def ExoName2seq():
    sys.path.append(r'D:\project\MotifFinding\kmerGraph')
    import seperateSample
    inputFn = r"D:\Dropbox\Motif analysis\2017-6-19 cow milk RNA rlt\sonal's data\exoName.txt"
    miRBaseFn = r'D:\project\MotifFinding\sources\exo_miRNA_raw_data\miRBase_pre_mature_IDs_conversion.txt'
    outputFn = r"D:\Dropbox\Motif analysis\2017-6-19 cow milk RNA rlt\sonal's data\userInput.fa"
    miRNANm2Seq, prevNm2seq = seperateSample.loadMiRNABase(miRBaseFn)
    with open(inputFn) as inputObj, open(outputFn, 'w') as outputObj:
        for line in inputObj:
            miRNm = line.strip().upper()
            if miRNm in miRNANm2Seq:
                seq = miRNANm2Seq[miRNm]
                outputLine = '>%s\n%s\n' % (miRNm, seq)
                outputObj.write(outputLine)
            elif miRNm in prevNm2seq:
                seq = prevNm2seq[miRNm]
                outputLine = '>%s\n%s\n' % (miRNm, seq)
                outputObj.write(outputLine)
            else:
                miRNm1 = miRNm + '-3P'
                miRNm2 = miRNm + '-5P'
                seq1 = miRNANm2Seq.get(miRNm1, prevNm2seq.get(miRNm1, ''))
                seq2 = miRNANm2Seq.get(miRNm2, prevNm2seq.get(miRNm2, ''))
                outputLine1 = '>%s\n%s\n' % (miRNm, seq1)
                outputLine2 = '>%s\n%s\n' % (miRNm, seq2)
                if seq1 and not seq2: outputObj.write(outputLine1)
                if seq2 and not seq1: outputObj.write(outputLine2)
                if not seq1 and not seq2:
                    print 'no match: %s' % miRNm
                elif seq1 and seq2:
                    print 'two match: %s' % miRNm
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

def Fetch3DBSeqCnt():
    # input
    fn = r'C:\Users\Administrator\Desktop\new 3.txt'
    # miRBaseFn = r'D:\project\MotifFinding\sources\exo_miRNA_raw_data\miRBase_pre_mature_IDs_conversion.txt'
    miRBaseFn = r'D:\Dropbox\project\MotifFinding\sources\exo_miRNA_raw_data\miRBase_pre_mature_IDs_conversion.txt'

    # output
    speciesMap = {'Homo sapiens':'hsa', 'Mus musculus': 'mmu', 'Bos taurus': 'bta'}

    miRNANm2Seq, prevNm2seq, prevNm2matNm = loadMiRNABase(miRBaseFn)
    with open(fn) as f:
        for line in f:
            line = line.strip()
            items = line.split('\t')
            prevMiRNm = items[0]
            spe = items[1]
            spe = speciesMap[spe]
            prevMiRNm = '%s-%s' % (spe, prevMiRNm)
            prevMiRNm = prevMiRNm.upper()
            matMiRNm = prevNm2matNm.get(prevMiRNm, prevMiRNm)
            print matMiRNm
#end_func

def formatMiRNm():
    inputFn = r'C:\Users\Administrator\Desktop\2\reg_names - Copy (2).txt'
    miRBaseFn = r'D:\Dropbox\project\MotifFinding\sources\exo_miRNA_raw_data\miRBase_pre_mature_IDs_conversion.txt'
    outputFn = r'C:\Users\Administrator\Desktop\2\reg_names_mature_name.txt'

    miRNANm2Seq, prevNm2seq, prevNm2matNm = loadMiRNABase(miRBaseFn)
    matNmLis = []
    with open(inputFn) as inputFileobj, open(outputFn, 'w') as outputFileobj:
        for line in inputFileobj:
            items = line.strip().split('\t')
            for prevNm in items:
                prevMiRNm = prevNm.upper()
                matMiRNm = prevNm2matNm.get(prevMiRNm, prevMiRNm).lower().replace('mir', 'miR')
                matNmLis.append(matMiRNm)
            outputLine = '\t'.join(matNmLis)
            outputFileobj.write('%s\n' % outputLine)
#end_func

def test():
    formatMiRNm()
#end_test

def main():
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
