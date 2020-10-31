#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/14
Description     :
    cmd: python GraphBasedMotifFinding.py XX 6 ACGT 100000 False
    cmd: python GraphBasedMotifFinding.py XX 6 ACGT 100000 True

"""
import sys
import os

import BioinfoComm, Comm
from Comm import print0, print9, PrintWithTime
import SignificaceEvaluation
import Backpack
import outputRlt

SINGLE_LENGTH_MODE = True
TARGET_LENGTH = 5

ENABLE_CPROFILE = False
ENABLE_PARALLEL = False
CORE_NUM = 16

ENABLE_MIR_SAMPLING = False

MIR_REF_FN = './data/allHumanMiRNA(HSA).fa'
RNA_REF_FN = './data/human_CDS.fasta'

class FetchInfoWithPattern:
    def __init__(self, userInputFn, alphabet, samplingFreq, targetPatternLis, outputFn):
        self.alphabet = alphabet
        self.samplingFreq = samplingFreq
        self.enableMiRSampling = False
        self.targetPatternLis = targetPatternLis
        self.outputFn = outputFn

        # input file
        self.userDataFn = userInputFn

        # global variables
        self.overallSigKmer = [{}, {}] if self.enableMiRSampling else [{}]
        self.wholeTree = [[], []] if self.enableMiRSampling else [[]]
        self.allInTreeKmerSet = [set(), set()] if self.enableMiRSampling else [set()]
        self.overallKmer2Cov = {}
        self.overallRefKmer2CovThreshold = [{}, {}] if self.enableMiRSampling else [{}]
        self.overallRefKmer2covIdxMatrix = [{}, {}] if self.enableMiRSampling else [{}]
        self.K2kmers = {}
        self.pattern2kmerSet = {}

        # initialization functions
        self.ShowBasicInfo()
    #end_init

    def ShowBasicInfo(self):
        PrintWithTime("start...", isLogging=True)
        print("current data file: %s" % self.userDataFn)
        print("miRNA ref data file: %s" % MIR_REF_FN)
        print("RNA ref data file: %s" % RNA_REF_FN)
        print("sampling count: %s" % self.samplingFreq)
        print("enable MiR Sampling: %s" % self.enableMiRSampling)
        print("enable single length mode: %s" % SINGLE_LENGTH_MODE)
        if SINGLE_LENGTH_MODE: print("single length: %s" % TARGET_LENGTH)
        print9(isLogging=True)
    #end_func

    # @profile
    def pipelineWithRE(self):
        # load user data
        userSeqLis, seqCnt, minSeqLen, maxSeqLen = BioinfoComm.loadSinglelineSeq(self.userDataFn)

        # load ref data
        sampledSeqMatrix = SignificaceEvaluation.SampleRNARef(RNA_REF_FN, seqCnt, self.samplingFreq, minSeqLen, maxSeqLen, self.alphabet)

        PrintWithTime("finish loading data", isLogging=True)
        print("number of sequences: %s" % seqCnt)
        print("min length: %s" % minSeqLen)
        print("max length: %s" % maxSeqLen)
        print9(isLogging=True)

        userPattern2cov = {}
        pattern2Pvalue = {}
        for pattern in self.targetPatternLis:
            # get user data pattern coverage
            userCov = BioinfoComm.FetchCovWithRE(pattern, userSeqLis)
            userPattern2cov[pattern] = userCov

            # get reference data pattern coverage
            covLis = []
            for refSeqLis in sampledSeqMatrix:
                refCov = BioinfoComm.FetchCovWithRE(pattern, refSeqLis)
                covLis.append(refCov)
            distriLis = outputRlt.FetchCovDistri(covLis, seqCnt)
            # rawPvalue = BioinfoComm.FetchPvalueFromBG(distriLis, userCov)
            rawPvalue = BioinfoComm.FetchPvalueFromBG_FindNonZero(distriLis, userCov) # new pvalue method
            adjustedPvalue = min(rawPvalue * len(self.targetPatternLis), 1.0)  # adjust P-value by multiplying the length
            pattern2Pvalue[pattern] = adjustedPvalue
            print('p-value of pattern %s is %s' % (pattern, adjustedPvalue))

        # output Pvalue
        print('outputting Pvalue')
        with open(self.outputFn, 'w') as outputFileobj:
            for pattern, pvalue in pattern2Pvalue.iteritems():
                cov = userPattern2cov[pattern]
                outputLine = '%s\t%s(in %s)\t%s' % (pattern, cov, seqCnt, pvalue)
                outputFileobj.write('%s\n' % outputLine)
    #end_func

    def pipelineNoRE(self):
        # load user data
        userSeqLis, seqCnt, minSeqLen, maxSeqLen = BioinfoComm.loadSinglelineSeq(self.userDataFn)

        # load ref data
        sampledSeqMatrix = []
        RNASampledSeqMatrix = SignificaceEvaluation.SampleRNARef(RNA_REF_FN, seqCnt, self.samplingFreq, minSeqLen, maxSeqLen, self.alphabet)
        sampledSeqMatrix.append(RNASampledSeqMatrix)

        PrintWithTime("finish loading data", isLogging=True)
        print("number of sequences: %s" % seqCnt)
        print("min length: %s" % minSeqLen)
        print("max length: %s" % maxSeqLen)
        print9(isLogging=True)

        refKmer2covIdxMatrix = {}
        userPattern2cov = {}
        kmer2seqIdSet = {}
        userKmer2Cov = {}
        userKmer2seqIdInt = {}
        pattern2IC = {}
        for pattern in self.targetPatternLis:
            print("now fetching IC and kmer coverage of pattern: %s" % pattern)
            # get user data pattern coverage
            patternLen = BioinfoComm.getPatternLength(pattern)
            kmers = BioinfoComm.FetchAllKmerFromPattern(pattern)
            kmer2visit = BioinfoComm.FetchSeqVis(userSeqLis, kmers)
            for kmer, visLis in kmer2visit.iteritems():
                kmer2seqIdSet.setdefault(kmer, set(map(lambda x: x[0], visLis)))
                userKmer2Cov.setdefault(kmer, len(kmer2seqIdSet[kmer]))
                if kmer not in userKmer2seqIdInt:
                    seqIdxSet = kmer2seqIdSet[kmer]
                    covIdxBitArray = BioinfoComm.IdxLis2bin(seqIdxSet, seqCnt)
                    covIdxBitInt = Comm.bitarray2int(covIdxBitArray)
                    userKmer2seqIdInt[kmer] = covIdxBitInt
            curUserPattern2cov, _ = outputRlt.FetchPatternSetCov(self.targetPatternLis, userSeqLis, userKmer2seqIdInt, seqCnt)
            userPattern2cov = dict(userPattern2cov, **curUserPattern2cov)

            # get reference data kmer coverage
            for kmer in kmers:
                covIdxMatrix, covLis = SignificaceEvaluation.FetchCovSeqDetail(sampledSeqMatrix[0], kmer, self.samplingFreq)
                refKmer2covIdxMatrix[kmer] = covIdxMatrix
            self.overallKmer2Cov.setdefault(patternLen, userKmer2Cov)

            # get reference data pattern IC
            IC, _ = Backpack.FetchPatternWeighedIC(pattern, {}, self.overallKmer2Cov)
            pattern2IC[pattern] = IC

        # fetch pattern's p-value and print
        print('calculating Pvalue of patterns')
        pattern2pvalue, __ = outputRlt.FetchPatternPvalue(0, seqCnt, pattern2IC, userPattern2cov, sampledSeqMatrix[0], refKmer2covIdxMatrix, self.samplingFreq, enableFilter=False)

        # output Pvalue
        print('outputting Pvalue')
        with open(self.outputFn, 'w') as outputFileobj:
            for pattern, pvalue in pattern2pvalue.iteritems():
                IC = pattern2IC[pattern]
                cov = userPattern2cov[pattern]
                outputLine = '%s\t%s(in %s)\t%s\t%s' % (pattern, cov, seqCnt, pvalue, IC)
                outputFileobj.write('%s\n' % outputLine)
    #end_func
#end_class

def filterOffset(offsetDic):
    """
    filter the patterns with a certain offset for just once
    """
    offsetFreqDic = {}
    for seqId, offset in offsetDic.iteritems():
        offsetFreqDic.setdefault(offset, 0)
        offsetFreqDic[offset] += 1
    offsetDic = {seqId: offset for seqId, offset in offsetDic.iteritems() if offsetFreqDic[offset] > 1}
    return offsetDic
#end_func

def FetchpatternNoRE(pattern2offsetFn, alphabet):
    # load offset file and userinput file to fetch patterns
    targetPatternSet = set()
    with open(pattern2offsetFn) as pattern2MIFileobj:
        for line in pattern2MIFileobj:
            items = line.strip().split('\t')
            oriP1 = items[0]
            oriP2 = items[2]
            offsetDic = eval(items[6])
            offsetDic = filterOffset(offsetDic)
            offsetLis = offsetDic.values()
            for offset in offsetLis:
                pattern = '%s%s%s' % (oriP1, ('[%s]' % alphabet) * offset, oriP2)
                targetPatternSet.add(pattern)
    return targetPatternSet
#end_func

def FetchpatternWithRE(pattern2offsetFn, minOffset, maxOffset):
    # load offset file and userinput file to fetch patterns
    targetPatternSet = set()
    with open(pattern2offsetFn) as pattern2MIFileobj:
        for line in pattern2MIFileobj:
            items = line.strip().split('\t')
            oriP1 = items[0]
            oriP2 = items[2]
            offsetDic = eval(items[6])
            if not offsetDic: continue
            pattern = '%s\w{%s,%s}%s' % (oriP1, minOffset, maxOffset, oriP2)
            targetPatternSet.add(pattern)
    return targetPatternSet
#end_func

def main():
    if len(sys.argv) > 1:
        _, dataDir, alphabet, samplingFreq, minOffset, maxOffset = sys.argv
        userInputFn = os.path.join(dataDir, 'userInput.fa')
        pattern2offsetFn = os.path.join(dataDir, 'patternPair2locOffsetInfo.txt')
        outputFn = os.path.join(dataDir, 'coexistPvalue.txt')
        samplingFreq = int(samplingFreq)
        minOffset = int(minOffset)
        maxOffset = int(maxOffset)
    else:
        # userInputFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\input\bindingSiteCluster1.fa'
        userInputFn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster3\userInput.fa'
        # pattern2offsetFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\co-exist motif\merged length\bindingSite\normalized mutual info\bindingSiteCluster1\patternPair2locOffsetInfo.txt'
        pattern2offsetFn = r'C:\Users\Administrator\Desktop\444\bindingSite_cluster3\patternPair2locOffsetInfo.txt'
        outputFn = r'C:\Users\Administrator\Desktop\444\bindingSite_cluster3\1.txt'
        alphabet = 'ACGU'
        samplingFreq = 100
        minOffset = 2
        maxOffset = 5

    # targetPatternSet = FetchpatternNoRE(pattern2offsetFn, alphabet)
    targetPatternSet = FetchpatternWithRE(pattern2offsetFn, minOffset, maxOffset)

    # fetch IC and P-value of patterns
    g = FetchInfoWithPattern(userInputFn, alphabet, samplingFreq, targetPatternSet, outputFn)
    g.pipelineWithRE()
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
