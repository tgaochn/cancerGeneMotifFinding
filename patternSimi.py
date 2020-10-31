#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/18
Description     :

"""
import BioinfoComm, Comm
import os
import math
import patternClustering
import DimerGraph
import logging

INFO = logging.getLogger('root').info
ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

SUDO_COUNT = 0.000001

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

def loadRlt(fn, motifLength):
    """
    load pattern info(IC, Pvalue...) from finalRlt.txt
    """
    pattern2cov = {}
    flag = False
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line[0] == 'p':
                if int(line[-1]) == motifLength:
                    flag = True
                else:
                    flag = False
            elif flag:
                items = line.split('\t')
                pattern = items[0]
                cov = int(items[1])
                pattern2cov[pattern] = cov
    return pattern2cov
#end_func

def loadRawRlt(fn):
    """
    load pattern info(IC, Pvalue...) from finalRlt.txt
    """
    pattern2cov = {}
    continueFlag = 1
    with open(fn) as f:
        for line in f:
            if line.strip() == '=====motif length:3=====': continueFlag = 0
            if continueFlag: continue
            if line.strip() == '=====motif length:4=====': break

            items = line.strip().split('\t')

            if items[0] == 'pattern':continue
            if items[0] == '=====motif length:3=====':continue

            pattern = items[0]
            cov = items[2]
            pattern2cov[pattern] = cov
    return pattern2cov
#end_func

def loadSinglePWM(pwmFn):
    """
    load a single PWM as a matrix
    """
    pwmKmerMatrix = []
    with open(pwmFn) as pwmF:
        for line in pwmF:
            line = line.strip()
            if not line or line[0] == '>': continue
            pwmKmerMatrix.append(line)
    return pwmKmerMatrix
#end_func

def loadSinglePPM(ppmFn):
    """
    load a single PWM as a matrix
    """
    ppmMatrix = []
    with open(ppmFn) as ppmFileobj:
        _, char1, char2, char3, char4 = ppmFileobj.readline().strip().split('\t')
        for line in ppmFileobj:
            items = line.strip().split('\t')
            items = map(lambda x:eval(x), items[1:])
            ppmMatrix.append(items)
    return ppmMatrix
#end_func

def loadAllPwm(pattern2cov, pwmDir, alphabet='ACGT'):
    """
    load all the PWM file in the dir
    colCnt[0] = [0, 312, 0, 138] means that in pos 0, there are 0 A, 312 C, 0 G, 138 T
    """
    pattern2pwm = {}
    pattern2colCnt = {}
    for pattern, cov in pattern2cov.iteritems():
        pwmFn = os.path.join(pwmDir, 'weighted_kmer-%s.fa' % pattern)
        pwmKmerMatrix = loadSinglePWM(pwmFn)
        pattern2colCnt[pattern] = []
        for colId in range(len(pwmKmerMatrix[0])):
            strInCol = map(lambda x:x[colId], pwmKmerMatrix)
            colCntLis = [strInCol.count(curChar) for curChar in alphabet]
            pattern2colCnt[pattern].append(colCntLis)
        pattern2pwm[pattern] = pwmKmerMatrix
    return pattern2pwm, pattern2colCnt
#end_func

def loadAllPPM(pattern2cov, ppmDir, alphabet='ACGT'):
    """
    load all the PWM file in the dir
    colCnt[0] = [0, 312, 0, 138] means that in pos 0, there are 0 A, 312 C, 0 G, 138 T
    """
    pattern2colCnt = {}
    for pattern, cov in pattern2cov.iteritems():
        ppmFn = os.path.join(ppmDir, 'PPM-%s.txt' % pattern)
        colCntLis = loadSinglePPM(ppmFn)
        pattern2colCnt[pattern] = colCntLis
    return pattern2colCnt
#end_func

def CalcSimi(pattern2colCntMatrix):
    """
    get the similarity of each pattern pair
    """
    patternSimiLis = []
    for patternX, colCntMatrixX in pattern2colCntMatrix.iteritems():
        for patternY, colCntMatrixY in pattern2colCntMatrix.iteritems():
            if patternX <= patternY:continue
            # simi = PearsonChiSquare(colCntMatrixX, colCntMatrixY)
            simi = PCC(colCntMatrixX, colCntMatrixY)
            patternSimiLis.append((patternX, patternY, simi))
    return patternSimiLis
#end_func

def PCC(colCntMatrixX, colCntMatrixY):
    """
    Pearson correlation coefficient
    """
    colPCCLis = []

    for colId, colCntX in enumerate(colCntMatrixX):
        colCntY = colCntMatrixY[colId]

        colProbX = map(lambda x:x * 1.0 / sum(colCntX), colCntX)
        colProbY = map(lambda x:x * 1.0 / sum(colCntY), colCntY)

        avgX = 1.0 / len(colProbX)
        avgY = 1.0 / len(colProbY)

        sdXLis = [(Xa - avgX) ** 2 for Xa in colProbX]
        sdX = math.sqrt(sum(sdXLis))
        sdYLis = [(Ya - avgY) ** 2 for Ya in colProbY]
        sdY = math.sqrt(sum(sdYLis))

        denominator = sdX * sdY
        numerator = sum([(Xa - avgX) * (colProbY[idx] - avgY) for idx, Xa in enumerate(colProbX)])

        colPCC = numerator * 1.0 / denominator
        colPCCLis.append(colPCC)
    #end_func

    simi = sum(colPCCLis)
    return simi
#end_func

def PearsonChiSquare(colCntMatrixX, colCntMatrixY):
    """
    Pearson Chi Square method
    """
    # TODO: this method has some problems
    colChiSquare = []
    chiSquare = 0

    for colId, colCntX in enumerate(colCntMatrixX):
        colCntY = colCntMatrixY[colId]
        Nx = sum(colCntX)
        Ny = sum(colCntY)
        N = Nx + Ny
        sumA_T = [cntX + colCntY[curChar] for curChar, cntX in enumerate(colCntX)]

        for curChar, nxj_0 in enumerate(colCntX):
            Nj = sumA_T[curChar]
            nxj_e = Nx * Nj * 1.0 / N
            if nxj_0 - nxj_e != 0:
                chiSquare += (nxj_0 - nxj_e) ** 2 / nxj_e

        for curChar, nyj_0 in enumerate(colCntY):
            Nj = sumA_T[curChar]
            nyj_e = Ny * Nj * 1.0 / N
            if nyj_0 - nyj_e != 0:
                chiSquare += (nyj_0 - nyj_e) ** 2 / nyj_e

        colChiSquare.append(chiSquare)
        chiSquare = 0
    #end_for

    simi = SUDO_COUNT + reduce(lambda x, y: x * y, colChiSquare)
    return simi
#end_func

def WriteCov(pattern2cov, pattern2covFn):
    """
    write dict(pattern2cov) into file 
    """
    with open(pattern2covFn, 'w') as fileObj:
        for pattern, cov in pattern2cov.iteritems():
            outputLine = '%s\t%s' % (pattern, cov)
            fileObj.write('%s\n' % outputLine)
#end_func

def WriteSimi(patternSimiLis, pattern2simiFn):
    """
    write list(patternSimiLis) into file 
    """
    with open(pattern2simiFn, 'w') as pattern2simiObj:
        for pattern1, pattern2, simi in patternSimiLis:
            if simi == 0: simi += SUDO_COUNT
            if simi == 1: simi -= SUDO_COUNT
            pattern2simiOutputLine = '%s\t%s\t%s' % (pattern1, pattern2, simi)
            pattern2simiObj.write('%s\n' % pattern2simiOutputLine)
#end_func

def fetchPwmSimiOneCell():
    """
    fetch patterns' similarity of one cell line using PWM
    """
    # input file
    patternRltFn = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\raw result\sw620-seq39\finalRlt - 4mer.txt'
    pwmDir = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\raw result\sw620-seq39\weightedKmer\4'

    # output file
    pattern2simiFn = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\separate simi 39 112\pattern2simi - sw620 - seq39.txt'
    pattern2covFn = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\separate simi 39 112\pattern2cov - sw620 - seq39.txt'

    # load file
    pattern2cov = loadRlt(patternRltFn)
    pattern2pwm, pattern2colCntMatrix = loadAllPwm(pattern2cov, pwmDir)
    # pattern2cov = {'[TG][TG]TG':1, '[TG]G[CG]A':1}

    # calculation
    patternSimiLis = CalcSimi(pattern2colCntMatrix)
    # normalization
    simiLis = map(lambda x:x[2], patternSimiLis)
    minSimi = min(simiLis)
    maxSimi = max(simiLis)
    # patternSimiLis = [(pattern1, pattern2, (simi - minSimi) * 1.0 / (maxSimi - minSimi)) for pattern1, pattern2, simi in patternSimiLis]

    # output
    WriteCov(pattern2cov, pattern2covFn)
    WriteSimi(patternSimiLis, pattern2simiFn)
#end_test


def mergeDict(pattern2cov1, patternTitle1, pattern2cov2, patternTitle2):
    mergedDict = {}
    for pattern, cov in pattern2cov1.iteritems():
        newPattern = '%s-%s' % (pattern, patternTitle1)
        mergedDict[newPattern] = cov
    for pattern, cov in pattern2cov2.iteritems():
        newPattern = '%s-%s' % (pattern, patternTitle2)
        mergedDict[newPattern] = cov
    return mergedDict
#end_func

def fetchPwmSimiTwoCell():
    """
    fetch patterns' similarity of two cell lines using PWM
    eg: sw620, seq39 vs. seq112
    """
    # # input file
    # patternRltFn1 = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\separate result\sw620-seq112\finalRlt - 4mer.txt'
    # patternRltFn2 = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\separate result\sw620-seq39\finalRlt - 4mer.txt'
    # pwmDir1 = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\separate result\sw620-seq112\weightedKmer\4'
    # pwmDir2 = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\separate result\sw620-seq39\weightedKmer\4'
    # patternTitle1 = '112'
    # patternTitle2 = '39'
    #
    # # output file
    # pattern2simiFn = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\pattern2simi - sw620 - seq39_112.txt'
    # pattern2covFn = r'D:\Dropbox\Motif analysis\2017-5-22 all cell line & merge39_112\merge 39_112\pattern2cov - sw620 - seq39_112.txt'

    # input file
    patternRltFn1 = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\merge_simi\3mer\miR_3mer.txt'
    patternRltFn2 = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\merge_simi\3mer\RNA_3mer.txt'
    pwmDir1 = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\miR\PPM\3'
    pwmDir2 = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\RNA\PPM\3'
    patternTitle1 = 'miR'
    patternTitle2 = 'RNA'

    # output file
    pattern2simiFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\merge_simi\3mer\pattern2simi.txt'
    pattern2covFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk\merge_simi\3mer\pattern2cov.txt'

    # load file
    pattern2cov1 = loadRlt(patternRltFn1, 3)
    pattern2cov2 = loadRlt(patternRltFn2, 3)
    pattern2cov = mergeDict(pattern2cov1, patternTitle1, pattern2cov2, patternTitle2)

    # pattern2pwm1, pattern2colCntMatrix1 = loadAllPwm(pattern2cov1, pwmDir1)
    pattern2colCntMatrix1 = loadAllPPM(pattern2cov1, pwmDir1)
    # pattern2pwm2, pattern2colCntMatrix2 = loadAllPwm(pattern2cov2, pwmDir2)
    pattern2colCntMatrix2 = loadAllPPM(pattern2cov2, pwmDir2)
    pattern2colCntMatrix = mergeDict(pattern2colCntMatrix1, patternTitle1, pattern2colCntMatrix2, patternTitle2)

    # calculation
    patternSimiLis = CalcSimi(pattern2colCntMatrix)
    # normalization
    simiLis = map(lambda x:x[2], patternSimiLis)
    minSimi = min(simiLis)
    maxSimi = max(simiLis)
    patternSimiLis = [(pattern1, pattern2, (simi - minSimi) * 1.0 / (maxSimi - minSimi)) for pattern1, pattern2, simi in patternSimiLis]

    # output
    WriteCov(pattern2cov, pattern2covFn)
    WriteSimi(patternSimiLis, pattern2simiFn)

#end_func


def fetchClusterCore(workDir, motifLength, alphabet='ACGT'):
    """
    fetch patterns' similarity of one cell line using PWM
    """
    # input file
    patternRltFn = os.path.join(workDir, 'finalRlt.txt')
    pwmDir = os.path.join(workDir, 'weightedKmer', '%s' % motifLength)
    userDataFn = os.path.join(workDir, 'userInput.fa')

    # load file
    pattern2cov = loadRawRlt(patternRltFn)
    pattern2pwm, pattern2colCntMatrix = loadAllPwm(pattern2cov, pwmDir, alphabet=alphabet)
    clusterFn = os.path.join(workDir, 'cluster.txt')

    # calculation
    patternSimiLis = CalcSimi(pattern2colCntMatrix)

    # build pattern similarity graph and do clustering
    pattern2pwmSimi = {(patter1, pattern2): simi for patter1, pattern2, simi in patternSimiLis}
    clusterLis = patternClustering.patternClustering(pattern2cov, pattern2pwmSimi)
    for cluster in clusterLis:
        outputLine = '\t'.join(cluster)
        INFO(outputLine)

    Comm.print9(isLogging=True)

    # output file
    pattern2simiFn = os.path.join(workDir, 'pattern2simi.txt')
    pattern2covFn = os.path.join(workDir, 'pattern2cov.txt')
    # output
    WriteCov(pattern2cov, pattern2covFn)
    WriteSimi(patternSimiLis, pattern2simiFn)

    # TODO: have modified here
    return

    # for each cluster, find some cores
    userSeqLis, seqCnt, _, _ = BioinfoComm.loadSinglelineSeq(userDataFn)
    kmerMatrix = DimerGraph.DivideSeq(userSeqLis)
    dimerGraph = DimerGraph.BuildGraph(kmerMatrix, alphabet=alphabet)
    K2Paths = DimerGraph.SearchPaths(dimerGraph, 8)

    userKmer2Cov, kmer2visDetail = DimerGraph.PathDic2PathCov(K2Paths, 8)
    overallUserKmer2VisDetail = {motifLength: kmer2visDetail}
    userKmer2seqIdSet = {kmer: set(map(lambda x: x[0], visLis)) for kmer, visLis in overallUserKmer2VisDetail[motifLength].iteritems()}
    for pattern, cov in pattern2cov.iteritems():
        curCov, curSeq, userKmer2seqIdSet = BioinfoComm.FetchPatternCov(pattern, userSeqLis, userKmer2seqIdSet)


    clusterCoreLis = patternClustering.findCore(clusterLis, userKmer2seqIdSet, seqCnt)

    motifSet = set()
    motifCandSet = set()
    for clusterCores in clusterCoreLis:

        for idx, (core, totalSeq, restSeq) in enumerate(clusterCores):
            if idx == 0:
                motifSet.add(core)
            else:
                motifCandSet.add(core)
            outputLis = [core, str(pattern2cov[core]), str(totalSeq), str(restSeq)]
            outputLine = '\t'.join(outputLis)
            INFO(outputLine)
        INFO('==')
    INFO('========')

    uncoveredSeqIdSet = set(range(seqCnt)) - reduce(lambda x, y: x | y, map(lambda x: userKmer2seqIdSet[x], motifSet))
    # the top pattern in each cluster already cover all sequences
    while motifCandSet and uncoveredSeqIdSet:
        bestMotifCand = motifCandSet.pop()
        motifCandSet.add(bestMotifCand)
        bestCov = uncoveredSeqIdSet - userKmer2seqIdSet[bestMotifCand]
        for curMotifCand in motifCandSet:
            curCov = uncoveredSeqIdSet - userKmer2seqIdSet[curMotifCand]

            if len(curCov) < len(bestCov):  # cover more sequences in the rest
                bestMotifCand = curMotifCand
                bestCov = curCov
            elif len(curCov) == len(bestCov) and pattern2cov[curMotifCand] > pattern2cov[bestMotifCand]:
                # cover equal sequences in the rest but have better overall coverage
                bestMotifCand = curMotifCand
        newUncoveredSeqIdSet = uncoveredSeqIdSet - userKmer2seqIdSet[bestMotifCand]
        if len(newUncoveredSeqIdSet) == len(uncoveredSeqIdSet): break
        uncoveredSeqIdSet = newUncoveredSeqIdSet
        motifSet.add(bestMotifCand)
        motifCandSet.remove(bestMotifCand)

    with open(clusterFn, 'w') as clusterFile:
        for motif in motifSet:
            outputLis = [motif, str(pattern2cov[motif])]
            outputLine = '\t'.join(outputLis)
            clusterFile.write('%s\n' % outputLine)
            INFO(outputLine)
    INFO('================')

    # for clusterCores in clusterCoreLis:
    #     for core, totalSeq, coveredSeq in clusterCores:
    #         outputLis = [core, str(pattern2cov[core]), str(totalSeq), str(coveredSeq)]
    #         outputLine = '\t'.join(outputLis)
    #         print outputLine
    #     print '==' * 10
#end_test

def main():
    # if PRINT_INTO_LOG:
    #     sys.stdout = open(LOG_FN, 'w')
        
    # fetchPwmSimiOneCell()
    # fetchPwmSimiTwoCell()
    # workDir = r'C:\Users\Administrator\Desktop\DKO'
    # motifLength = 4
    # fetchClusterCore(workDir, motifLength, alphabet='ACGT')
    # fetchClusterCore(workDir, motifLength, alphabet='ACGU')
    a = loadRlt(r'D:\Dropbox\Motif analysis\final result\2. result part\2.5 cow milk(unfinished)\RNA\finalRlt.txt', 6)
    print a





    # if PRINT_INTO_LOG:
    #     sys.stdout.close()
#end_main

def RunAllCluter():
    initInputDir = r'C:\Users\Administrator\Desktop\2017-5-17 different starting node'
    folderLis = os.listdir(initInputDir)
    for folder in folderLis:
        folder = os.path.join(initInputDir, folder)
        print folder
        fetchClusterCore(folder, 3, alphabet='ACGT')

#end_func

def GetAllSimi():
    initInputDir = r'C:\Users\Administrator\Desktop\all cell line\finished'
    pattern2simiFn = r'C:\Users\Administrator\Desktop\all cell line\pattern2simi - allCell.txt'
    pattern2covFn = r'C:\Users\Administrator\Desktop\all cell line\pattern2cov - allCell.txt'

    folderLis = os.listdir(initInputDir)
    pattern2cov = {}
    pattern2colCntMatrix = {}
    # folderLis = ['evpdia_hsa_Caput_epithelial_cell_input', 'evpdia_hsa_Caput_luminal_fluid_input']
    for folder in folderLis:
        fullDir = os.path.join(initInputDir, folder)
        clusterFn = os.path.join(fullDir, 'cluster.txt')
        curPattern2cov = {}

        for line in open(clusterFn):
            items = line.strip().split('\t')
            pattern = items[0]
            cov = items[1]
            curPattern2cov[pattern] = cov

        pwmDir = os.path.join(fullDir, 'weightedKmer', '4')
        curPattern2pwm, curPattern2colCntMatrix = loadAllPwm(curPattern2cov, pwmDir)
        curPattern2colCntMatrix = {'%s_%s' % (k, folder): v for k, v in curPattern2colCntMatrix.iteritems()}
        curPattern2cov = {'%s_%s' % (k, folder): v for k, v in curPattern2cov.iteritems()}

        pattern2colCntMatrix = dict(pattern2colCntMatrix, **curPattern2colCntMatrix)
        pattern2cov = dict(pattern2cov, **curPattern2cov)

    patternSimiLis = CalcSimi(pattern2colCntMatrix)
    WriteCov(pattern2cov, pattern2covFn)
    WriteSimi(patternSimiLis, pattern2simiFn)



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
        # main()
        # RunAllCluter()
        # GetAllSimi()
        fetchPwmSimiTwoCell()
#end_if
