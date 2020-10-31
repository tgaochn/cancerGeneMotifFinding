#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName:

Author:
    Tian Gao
Email:
    tgaochn@gmail.com
CreationDate:
    12/10/2017
Description:

"""

import BioinfoComm

def func1():
    fn1 = r'D:\Dropbox\Motif analysis\final result\2. result part\2.2 sw620 rlt\sw620_seq112\userInput.fa'
    fn2 = r'D:\Dropbox\tmp\Arid3a_3875.1_v1_deBruijn.txt'
    outputFn = r'D:\Dropbox\tmp\1_output.txt'

    seqLis, seqCnt, minLen, maxLen = BioinfoComm.loadSinglelineSeq(fn1)

    with open(outputFn, 'w') as outputFileobj, open(fn2) as fileobj2:
        lines = fileobj2.readlines()
        for seqId, seq in enumerate(seqLis):
            scoreLine = lines[seqId]
            score = scoreLine.split('\t')[0]
            outputFileobj.write('%s\t%s\n' % (score, seq))
#end_func

def func2():
    fn = 'Comm.py'
    with open(fn) as f:
        line = f.readline()
        print line
        for lineIdx, line in enumerate(f):
            if lineIdx > 10:break
        print line
        print f.readline()
        print f.readline()
#end_func

def func3():
    # fn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\toolComparison\MERCI - done\1.txt'
    # fn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\toolComparison\MERCI - done\2.txt'
    fn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\toolComparison\MERCI - done\3.txt'

    motifDic = {}
    with open(fn) as f:
        for line in f:
            charLis = line.strip().split(' ')
            kmer = ''.join(charLis)
            motifLen = len(kmer)
            motifDic.setdefault(motifLen, [])
            motifDic[motifLen].append(kmer)

    motifMatDic = {}
    for motifLen, kmerLis in motifDic.iteritems():
        motifMatDic.setdefault(motifLen, [])
        motif = BioinfoComm.SegmentMerge(kmerLis)
        print motifLen
        print kmerLis
        print motif
        print '=='
        for i in range(motifLen):
            tmpDic = {'A':0, 'C':0, 'G':0, 'T':0}
            for kmer in kmerLis:
                curChar = kmer[i]
                tmpDic[curChar] += 1
            motifMatDic[motifLen].append(tmpDic)

        print 'PO\tA\tC\tG\tU'
        for idx, curDic in enumerate(motifMatDic[motifLen]):
            print '%s\t%s\t%s\t%s\t%s' % (idx + 1, curDic['A'], curDic['C'], curDic['G'], curDic['T'])

        print '=' * 20


def main():
    # func1()
    # func2()
    func3()
#end_main

if __name__ == "__main__":
    main()
#end_if
