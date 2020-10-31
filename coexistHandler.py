#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/7/19
Description     :

"""
import copy
import sys
import math
import re
import BioinfoComm, Comm
import webServer
from Comm import print0
import os
import patternSimi
from scipy.special import comb
import FetchInfoGivenPattern
from fun import func2

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

def FetchCopvalueAndDistri(inputDir, offsetTuple):
    outputFn = os.path.join(inputDir, 'patternPair2info.txt')
    pvalueFn = os.path.join(inputDir, 'coexistPvalue.txt')
    miFn = os.path.join(inputDir, 'patternPair2mutualInfo.txt')
    minOffset, maxOffset = offsetTuple

    patternPair2mi = {}
    with open(pvalueFn) as pvalueFileobj, open(miFn) as miFileobj, open(outputFn, 'w') as outputFileobj:
        # load MI
        for line in miFileobj:
            items = line.strip().split('\t')
            p1 = items[0]
            p2 = items[2]
            mi = items[5]
            patternPair2mi[(p1, p2)] = mi

        titleLis = ['pattern1', 'pattern2', 'offset', 'mergedPattern', 'coverage', 'seqCnt', 'covPercentage', 'pvalue', 'mutualInfo']
        titleLine = '\t'.join(titleLis)
        outputFileobj.write('%s\n' % titleLine)

        for line in pvalueFileobj:
            items = line.strip().split('\t')
            patternPair = tuple(items[0].split('\w{%s,%s}' % (minOffset, maxOffset)))
            covStr = items[1]
            cov = int(covStr.partition('(')[0])
            seqCnt = int(covStr.partition(' ')[2][:-1])
            covPrec = '%s%%' % int(cov * 100 / seqCnt)
            pvalue = float(items[2])
            mi = patternPair2mi[(patternPair[0], patternPair[1])]
            offsetStr = '[%s, %s]' % (minOffset, maxOffset)
            outputLine = '\t'.join([patternPair[0], patternPair[1], offsetStr, items[0], str(cov), str(seqCnt), str(covPrec)[:5], str(pvalue), mi])
            outputFileobj.write('%s\n' % outputLine)

    # outputDistri(miDistriLis)
    # outputDistri(covPrecDistriLis)
    # outputNm(miDistriLis)
    # outputNm(covPrecDistriLis)
#end_func

def outputNm(distriLis):
    maxLen = max(map(lambda x: len(x), distriLis))
    for i in xrange(maxLen):
        curItems = []
        for curDatalis in distriLis:
            if len(curDatalis) <= i:
                curStr = '0'
            else:
                curStr = '%2f' % curDatalis[i]
            curItems.append(curStr)
        curLine = '\t'.join(curItems)
        print curLine
# end_func

def outputDistri(distriLis):
    rltLis = []
    partitionCnt = 25
    for curDatalis in distriLis:
        boundryLis = range(partitionCnt)
        boundryLis = map(lambda x: x * 1.0 / partitionCnt, boundryLis)
        countLis, boundaryLis = Comm.CountAmountInRanges(curDatalis, boundryLis)
        rltLis.append(countLis)

    for i in xrange(partitionCnt):
        curLine = []
        for curDatalis in rltLis:
            if len(curDatalis) <= i:
                curStr = ''
            else:
                curStr = '%2f' % curDatalis[i]
            curLine.append(curStr)
        curLine = '\t'.join(curLine)
        print curLine
# end_func

def FetchAllClusterCoexist():
    """
    fetch co-exist pattern in each cluster of a certain RNA
    input: result data
    output: raw coexist data
    """
    baseinputDir = r'C:\Users\Administrator\Desktop\coexist\result\mir16'
    # baseinputDir = r'C:\Users\Administrator\Desktop\coexist1\input'
    # baseoutputDir = r'C:\Users\Administrator\Desktop\coexist'
    baseoutputDir = r'C:\Users\Administrator\Desktop\coexist\raw coexist\mir16'

    for localDir in os.listdir(baseinputDir):
        curInputDir = os.path.join(baseinputDir, localDir)
        curOutputDir = os.path.join(baseoutputDir, localDir)
        if not os.path.exists(curOutputDir): os.makedirs(curOutputDir)
        FetchCoexist(curInputDir, curOutputDir)
#end_func

def FetchOneDatasetCoexist(curInputDir, curOutputDir, offsetTuple):
    """
    fetch coexist motif in a certain dataset
    :return:
    """
    # curInputDir = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster1'
    # curInputDir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\result\normal\cluster3'
    # curInputDir = r'/home/tgao/GraphBasedMotifFinding/output/coexistInput/bindingSiteCluster2'
    # curOutputDir = r'C:\Users\Administrator\Desktop\444\bindingSite_cluster3'
    # curOutputDir = r'/home/tgao/GraphBasedMotifFinding/output/finalRlt/coexist/bindingSiteClass2'
    # offsetTuple = (7, 13)
    FetchCoexist(curInputDir, curOutputDir, offsetTuple)
#end_func

def FetchCoexist(inputDir, outputDir, offsetTuple=(2, 5)):
    """
    fetch co-exist pattern from input
    output pattern pair coexist count and mutual information
    """
    inputFn = os.path.join(inputDir, 'userInput.fa')
    rltFn = os.path.join(inputDir, 'finalRlt.txt')
    logFn = os.path.join(inputDir, 'log.txt')
    seqLis, seqCnt, titleLis = BioinfoComm.loadMultilineSeq(inputFn)
    lengthLis = [2, 3, 4, 5, 6]
    length2patterDic = {}
    pvalueOutputFn = os.path.join(outputDir, 'patternPair2locOffsetInfo.txt') # this pvalue is calculated by formula, which doesn't work
    miOutputFn = os.path.join(outputDir, 'patternPair2mutualInfo.txt')
    if not os.path.exists(outputDir): os.makedirs(outputDir)

    # handle with 2mer
    with open(logFn) as logFileobj:
        line = logFileobj.readline()
        while 'begin to calculate on layer motifLength = 2' not in line:
            line = logFileobj.readline()
        logFileobj.readline()
        logFileobj.readline()
        line = logFileobj.readline()
        dimerSet = eval(line.split(':')[-1])
        dimerSet.remove('__')
        allPattern2cov, _ = BioinfoComm.FetchPatternCovInSeqLis(dimerSet, seqLis)
        length2patterDic[2] = allPattern2cov

    # handle with others
    for motifLength in lengthLis[1:]:
        pattern2cov = patternSimi.loadRlt(rltFn, motifLength)
        length2patterDic[motifLength] = pattern2cov
        allPattern2cov = dict(allPattern2cov, **pattern2cov)

    patternPair2locOffsetInfo = []
    patternPair2mi = []
    # for longLength in lengthLis[1:]:
    for longLength in lengthLis:
        longPatternLis = length2patterDic[longLength].keys()
        for shortLength in lengthLis:
            # if longLength <= shortLength: continue
            if longLength < shortLength: continue
            shortPatternLis = length2patterDic[shortLength].keys()

            # fetch location
            longPattern2locLis = findLoc(seqLis, longPatternLis)
            shortPattern2locLis = findLoc(seqLis, shortPatternLis)

            # fetch cov seqID
            allPattern2locLis = dict(longPattern2locLis, **shortPattern2locLis)
            pattern2covSeqIdSet = {}
            for pattern, locLis in allPattern2locLis.iteritems():
                pattern2covSeqIdSet.setdefault(pattern, set())
                for seqId in xrange(seqCnt):
                    if locLis[seqId]:
                        pattern2covSeqIdSet[pattern].add(seqId)

            # fetch offset
            patternPair2locOffsetLis = {}
            for longPattern, longLocLis in longPattern2locLis.iteritems():
                for shortPattern, shortLocLis in shortPattern2locLis.iteritems():
                    if longPattern == shortPattern: continue
                    patternPair2locOffsetLis.setdefault((longPattern, shortPattern), [])
                    patternPair2locOffsetLis.setdefault((shortPattern, longPattern), [])
                    for seqId in xrange(len(longLocLis)):
                        curLongLocLis = longLocLis[seqId]
                        curShortLocLis = shortLocLis[seqId]
                        minOffset = Comm.FetchLisMinOffset(curLongLocLis, curShortLocLis)
                        patternPair2locOffsetLis[(longPattern, shortPattern)].append(minOffset)
                        if longLength == shortLength: continue
                        minOffset = Comm.FetchLisMinOffset(curShortLocLis, curLongLocLis)
                        patternPair2locOffsetLis[(shortPattern, longPattern)].append(minOffset)

            # format offset
            inRangeCnt = 0
            outRangeCnt = 0
            for (pattern1, pattern2), locOffsetLis in patternPair2locOffsetLis.iteritems():
                # print (pattern1, pattern2), len(locOffsetLis)
                motifLength = BioinfoComm.getPatternLength(pattern1)
                # curLocOffsetDic = {seqId + 1: offset - motifLength for seqId, offset in enumerate(locOffsetLis) if offset != -1 and 2 + motifLength <= offset <= 5 + motifLength}
                curLocOffsetDic = {seqId + 1: offset - motifLength for seqId, offset in enumerate(locOffsetLis) if offset != -1 and offsetTuple[0] + motifLength <= offset <= offsetTuple[1] + motifLength}
                coexistCnt = len(curLocOffsetDic)
                cov1 = allPattern2cov[pattern1]
                cov2 = allPattern2cov[pattern2]
                patternStr = '%s-%s' % (pattern1, pattern2)

                # use formula to get pvalue for co-exist pattern - this formula for pvalue doesn't work
                coexistSeqIdSet = set(map(lambda x:x-1, curLocOffsetDic.keys()))
                mi = FetchMutualInfomation(pattern1, pattern2, coexistCnt, coexistSeqIdSet, cov1, pattern2covSeqIdSet[pattern1], cov2, pattern2covSeqIdSet[pattern2], seqCnt, set(range(seqCnt)))
                rawPvalue = fetchCoexistRawPvalue(coexistCnt, cov1, cov2, seqCnt, patternStr)
                # print pattern1, pattern2, mi, 'j=%s, n1=%s, n2=%s' % (coexistCnt, cov1, cov2)
                if rawPvalue >= 0:
                    inRangeCnt += 1
                    adjPvalue = rawPvalue * len(patternPair2locOffsetLis)
                else:
                    outRangeCnt += 1
                    adjPvalue = rawPvalue
                patternPair2mi.append((pattern1, cov1, pattern2, cov2, coexistCnt, mi))
                # if adjPvalue > 0.05: continue
                # print pattern1, cov1, pattern2, cov2, coexistCnt, adjPvalue
                patternPair2locOffsetInfo.append((pattern1, cov1, pattern2, cov2, coexistCnt, adjPvalue, curLocOffsetDic))
            # print inRangeCnt, outRangeCnt

    # output
    patternPair2locOffsetInfo.sort(key=lambda x:x[2], reverse=True)
    patternPair2mi.sort(key=lambda x:x[5], reverse=True)
    with open(pvalueOutputFn, 'w') as pvalueOutputFileobj, open(miOutputFn, 'w') as miOutputFileobj:# this pvalue is calculated by formula, which doesn't work
        for pattern1, cov1, pattern2, cov2, coexistCnt, adjPvalue, curLocOffsetDic in patternPair2locOffsetInfo:
            cov1 = '%s(in %s)' % (cov1, seqCnt)
            cov2 = '%s(in %s)' % (cov2, seqCnt)
            coexistCnt = '%s(%s%%)' % (coexistCnt, coexistCnt * 100 / seqCnt)
            outputLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (pattern1, cov1, pattern2, cov2, coexistCnt, adjPvalue, curLocOffsetDic)
            pvalueOutputFileobj.write('%s\n' % outputLine)
        for pattern1, cov1, pattern2, cov2, coexistCnt, mi in patternPair2mi:
            cov1 = '%s(in %s)' % (cov1, seqCnt)
            cov2 = '%s(in %s)' % (cov2, seqCnt)
            coexistCnt = '%s(%s%%)' % (coexistCnt, coexistCnt * 100 / seqCnt)
            outputLine = '%s\t%s\t%s\t%s\t%s\t%s' % (pattern1, cov1, pattern2, cov2, coexistCnt, mi)
            miOutputFileobj.write('%s\n' % outputLine)
        print0()
#end_func

def findLoc(seqLis, patternLis):
    """
    get a dict including location information
    {
        pattern: [(seqTitle, [matchIdx1, matchIdx2])]
    }
    """
    pattern2locLis = {}
    for pattern in patternLis:
        pattern2locLis.setdefault(pattern, [])
        patternLen = webServer.GetPatternLength(pattern)
        for seq in seqLis:
            matchSegSet = set(re.findall(pattern, seq))
            idxLis = []
            for matchSeg in matchSegSet:
                tmpSeq = copy.deepcopy(seq)
                lastMatchIdx = 0
                while matchSeg in tmpSeq:
                    matchIdx = tmpSeq.index(matchSeg)
                    idxLis.append(matchIdx + lastMatchIdx)
                    lastMatchIdx = matchIdx + lastMatchIdx + patternLen
                    tmpSeq = tmpSeq[matchIdx + patternLen:]
            pattern2locLis[pattern].append(idxLis)
    return pattern2locLis
#end_func


def FetchMutualInfomation(pattern1, pattern2, coexistCnt, coexistSeqIdSet, cov1, p1CovSeqId, cov2, p2CovSeqId, seqCnt, allSeqId):
    # version 1: separate px, py, pxy
    # px1 = cov1 * 1.0 / seqCnt
    # px0 = 1 - cov1 * 1.0 / seqCnt
    # py1 = cov2 * 1.0 / seqCnt
    # py0 = 1 - cov2 * 1.0 / seqCnt
    #
    # pxy11 = coexistCnt * 1.0 / seqCnt
    # pxy01 = py1 - pxy11
    # pxy10 = px1 - pxy11
    # pxy00 = py0 - pxy10
    # if pxy00 + pxy01 != px0: print 'error'

    if not coexistSeqIdSet:
        return 0.0

    # version 2: only consider pxy
    pxy10IdSet = p1CovSeqId - p2CovSeqId
    pxy01IdSet = p2CovSeqId - p1CovSeqId
    pxy00IdSet = allSeqId - p1CovSeqId - p2CovSeqId | (p1CovSeqId.intersection(p2CovSeqId) - coexistSeqIdSet)

    # sudo-count
    if not pxy10IdSet:
        pxy10IdSet.add(pxy00IdSet.pop())
    if not pxy01IdSet:
        pxy01IdSet.add(pxy00IdSet.pop())

    pxy11 = coexistCnt * 1.0 / seqCnt
    pxy10 = len(pxy10IdSet) * 1.0 / seqCnt
    pxy01 = len(pxy01IdSet) * 1.0 / seqCnt
    pxy00 = len(pxy00IdSet) * 1.0 / seqCnt
    px1 = pxy11 + pxy10
    px0 = pxy01 + pxy00
    py1 = pxy11 + pxy01
    py0 = pxy10 + pxy00

    Hx = GetEntrop(px1, px0)
    Hy = GetEntrop(py1, py0)

    mi00 = GetSingleMi(pxy00, px0, py0) if min(pxy00, px0, py0) > 0 else 0
    mi01 = GetSingleMi(pxy01, px0, py1) if min(pxy01, px0, py1) > 0 else 0
    mi10 = GetSingleMi(pxy10, px1, py0) if min(pxy10, px1, py0) > 0 else 0
    mi11 = GetSingleMi(pxy11, px1, py1) if min(pxy11, px1, py1) > 0 else 0
    Ixy = mi00 + mi01 + mi10 + mi11

    nmi = 2 * Ixy * 1.0 / (Hx + Hy)

    return nmi
    # return Ixy
#end_func

def GetSingleMi(pxy, px, py):
    return pxy * math.log(pxy / (px * py), 2)
#end_func

def GetEntrop(px1, px0):
    return - px1 * math.log(px1, 2) - px0 * math.log(px0, 2)
#end_func

def fetchCoexistRawPvalue(j, n1, n2, n, patternStr):
    """
    :param j: co-exist cnt
    :param n1: pattern1 cov
    :param n2: pattern2 cov
    :param n: total seq cnt
    """
    rawPvalue = 0
    # formula 1 & 2
    jLis = range(n + 1)[j:]
    # jLis = filter(lambda x:max(0, n1 + n2 - n) <= x <= min(n1, n2), jLis) # for formula 1
    if not jLis:
        print 'not in range', patternStr, max(0, n1 + n2 - n), min(n1, n2), n - j + 1, j
        return 1
    for j in jLis:
        # formula 1 from paper
        # comb_n_n2 = comb(n - n2, n1 - j) if comb(n - n2, n1 - j) > 0 else 1
        # pj = 1.0 * comb(n, j) / comb(n, n2) * comb(n - j, n2 - j) / comb(n, n1) * comb_n_n2
        # rawPvalue += pj

        # formula 2 from paper
        comb_n_n1 = comb(n - n1, n2 - j) if comb(n - n1, n2 - j) > 0 else 1
        pj = 1.0 * comb(n1, j) * comb_n_n1 / comb(n, n2)
        rawPvalue += pj

    # print 'in range', patternStr, rawPvalue
    return rawPvalue
#end_func

def tmpFormat():
    covFn = r'D:\Dropbox\Motif analysis\final result\3. supporting materials\Figure_code\data\covPlot.txt'
    miFn = r'D:\Dropbox\Motif analysis\final result\3. supporting materials\Figure_code\data\miPlot.txt'
    combFn = r'D:\Dropbox\Motif analysis\final result\3. supporting materials\Figure_code\data\mi_cov.txt'
    classMap = {0: 'I', 1: 'II', 2: 'III', 3: 'IV', 4: 'V'}
    with open(covFn) as covFileobj, open(miFn) as miFileobj, open(combFn, 'w') as combFileobj:
        covFileobj.readline()
        miFileobj.readline()
        combFileobj.write('mi\tcov\tclass\n')
        id2covLis = {}
        for lineIdx, line in enumerate(covFileobj):
            items = line.strip().split('\t')
            covLis = items[1:]
            id2covLis[lineIdx] = covLis

        for lineIdx, line in enumerate(miFileobj):
            items = line.strip().split('\t')
            miLis = items[1:]
            covLis = id2covLis.get(lineIdx, ['0'] * 5)

            for classId, mi in enumerate(miLis):
                cov = covLis[classId]
                classNm = classMap[classId]
                if eval(mi) == 0 and eval(cov) == 0: continue
                outputLine = '\t'.join([mi, cov, classNm])
                combFileobj.write('%s\n' % outputLine)

#end_func

def formatScatterPlot():
    """
    this function generates data for scatter plot
    """
    coexistDir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\co-exist motif\merged length\bindingSite\coexistWithDiffOffsetAll(4nd version)\dataForScatterPlot'
    outputFn = os.path.join(coexistDir, 'mi_cov_significant.txt')
    classDic = {}
    for localDir in os.listdir(coexistDir):
        curInputDir = os.path.join(coexistDir, localDir)
        if not os.path.isdir(curInputDir): continue
        coexistFn = os.path.join(curInputDir, 'patternPair2info.txt')
        classDic[localDir] = []
        with open(coexistFn) as coexistFileobj:
            for line in coexistFileobj:
                items = line.strip().split('\t')
                cov = items[6]
                cov = eval(cov[:-1]) * 0.01
                pvalue = eval(items[7])
                if pvalue > 0.05: continue
                mi = items[8]
                classDic[localDir].append((mi, cov))

    with open(outputFn, 'w') as outputFileobj:
        outputFileobj.write('mi\tcov\tclass\n')
        for classNm, dataLis in classDic.iteritems():
            for mi, cov in dataLis:
                outputLine = '%s\t%s\t%s' % (mi, cov, classNm)
                outputFileobj.write('%s\n' % outputLine)
#end_func

def func1():
    if len(sys.argv) > 1:
        _, inputDir, outputDir, alphabet, samplingFreq, minOffset, maxOffset = sys.argv
        samplingFreq = int(samplingFreq)
        offsetTuple = (int(minOffset), int(maxOffset))
    else:
        clusterIdx = 4

        # inputDir = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster%s' % clusterIdx
        # inputDir = r'coexist/input/cluster1'
        inputDir = r'coexist/input/cluster%s' % clusterIdx
        # outputDir = r'cooccurance/cluster%s' % clusterIdx
        # outputDir = r'coexist/output/cluster1'
        outputDir = r'coexist/output/cluster%s' % clusterIdx
        offsetTuple1 = (2, 5)
        offsetTuple2 = (4, 7)
        offsetTuple3 = (7, 13)
        offsetTuple4 = (2, 7)
        offsetTuple5 = (2, 5)

        offsetTuple = [offsetTuple1, offsetTuple2, offsetTuple3, offsetTuple4, offsetTuple5][clusterIdx - 1]

        alphabet = 'ACGU'
        # samplingFreq = 100
        samplingFreq = 100000

    print inputDir
    print outputDir
    print offsetTuple
    print alphabet
    print samplingFreq


    userInputFn = os.path.join(inputDir, 'userInput.fa')
    pattern2offsetFn = os.path.join(outputDir, 'patternPair2locOffsetInfo.txt')
    pvalueRltFn = os.path.join(outputDir, 'coexistPvalue.txt')

    FetchOneDatasetCoexist(inputDir, outputDir, offsetTuple)
    targetPatternSet = FetchInfoGivenPattern.FetchpatternWithRE(pattern2offsetFn, offsetTuple[0], offsetTuple[1])

    # fetch IC and P-value of patterns
    g = FetchInfoGivenPattern.FetchInfoWithPattern(userInputFn, alphabet, samplingFreq, targetPatternSet, pvalueRltFn)
    g.pipelineWithRE()
    FetchCopvalueAndDistri(outputDir, offsetTuple)


    # tmpFormat()
    # formatScatterPlot()
#end_test

def func2():
    offsetTuple1 = (2, 5)
    offsetTuple2 = (4, 7)
    offsetTuple3 = (7, 13)
    offsetTuple4 = (2, 7)
    offsetTuple5 = (2, 5)

    for clusterIdx in range(6)[1:]:
        outputDir = r'coexist/output/cluster%s' % clusterIdx
        offsetTuple = [offsetTuple1, offsetTuple2, offsetTuple3, offsetTuple4, offsetTuple5][clusterIdx - 1]
        FetchCopvalueAndDistri(outputDir, offsetTuple)
#end_func


def main():
    # func1()
    func2()
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
