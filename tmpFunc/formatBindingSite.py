#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/8/3
Description     :

"""
import sys
import pdb
import BioinfoComm, Comm
from Comm import print0

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

def handleBlock(curLines, outputFileobj, seqID2Class):
    if 'MIMAT0000092' not in curLines[0]: return
    seqID, _, _ = curLines[0].partition('_MIMA')
    miRLen = int(curLines[2].split('\t')[-1])
    flagStr = curLines[4].split('\t')[0][:miRLen]
    formattedFlagStr = flagStr.replace('.', '0').replace('(', '1')
    classID = seqID2Class[seqID]
    outputLine = '%s\t%s\t%s' % (seqID, formattedFlagStr, classID)
    outputFileobj.write('%s\n' % outputLine)
#end_func

def fun1():
    seqID2ClassFn = r'C:\Users\Administrator\Desktop\bindingSite\mmc1.txt'
    bindingInfoFn = r'C:\Users\Administrator\Desktop\bindingSite\mmc2.txt'
    outputFn = r'C:\Users\Administrator\Desktop\bindingSite\bindingFlag_miR92a.txt'

    with open(seqID2ClassFn) as seqID2ClassFileobj, open(bindingInfoFn) as bindingInfoFileobj, open(outputFn, 'w') as outputFileobj:
        seqID2Class = {}
        for line in seqID2ClassFileobj:
            line = line.strip()
            if not line or line[0] == '#': continue
            items = line.split('\t')
            if items[0] == 'seq_ID': continue
            seqID = items[0]
            foldingClass = items[21]
            seqID2Class[seqID] = foldingClass

        curLines = []
        for line in bindingInfoFileobj:
            line = line.strip()
            if not line:
                handleBlock(curLines, outputFileobj, seqID2Class)
                curLines = []
            else:
                curLines.append(line)
        if curLines:
            handleBlock(curLines, outputFileobj, seqID2Class)
#end_func

def test():
    fun1()
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
