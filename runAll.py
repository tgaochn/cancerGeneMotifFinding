#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/22
Description     :

        qsub -o /home/tgao/GraphBasedMotifFinding/output/runJob -e /home/tgao/GraphBasedMotifFinding/output/runJob/err.log -wd /home/tgao/GraphBasedMotifFinding run.sh

"""
import sys
import os
import time

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

WORKING_DIR = r'/home/tgao/GraphBasedMotifFinding'
SCRIPT_FN = r'/home/tgao/GraphBasedMotifFinding/GraphBasedMotifFinding.py'

def RunAllCellLine():
    # input dir
    initInputDir = '/home/tgao/GraphBasedMotifFinding/allCellLine/input/unfinished'

    # output dir
    initOutputDir = '/home/tgao/GraphBasedMotifFinding/allCellLine/output/unfinished'

    fileLis = os.listdir(initInputDir)
    fileLis = filter(lambda x:x[-2:] == 'fa', fileLis)
    # fileLis = ['u2t_hnRNPA2B1_exomiRNA_v2_formatted.fa', 'cell_2016_103seq_GGCU.fa', 'elife_mirna_seq163.fa'] + fileLis
    fileLis = fileLis

    for inputFn in fileLis:
        baseNm = os.path.splitext(inputFn)[0]
        outputDir = os.path.join(initOutputDir, baseNm)
        inputFn = os.path.join(initInputDir, inputFn)
        if not os.path.exists(outputDir): os.makedirs(outputDir)

        cmd = 'cp \"%s\" \"%s/userInput.fa\"' % (inputFn, outputDir)
        # print cmd
        os.system(cmd)

        cmd = 'python GraphBasedMotifFinding.py \"%s\" 8 ACGT 100000' % outputDir
        # print cmd
        os.system(cmd)
#end_test

def RunAllRNA():
    # input dir
    initInputDirLis = [
        (r'/home/tgao/GraphBasedMotifFinding/paperValidation/hnRNPA2B1_v2_seq30_GGAG/', 'ACGU', True),
        (r'/home/tgao/GraphBasedMotifFinding/paperValidation/cell2016_seq103_GGCU/', 'ACGU', True),
    ]

    RunJob(initInputDirLis)
#end_func

def RunJob(inputDirLis, sampling=100000):
    for inputDir, alphabet, enableRNASampling in inputDirLis:
        cmd = 'python GraphBasedMotifFinding.py \"%s\" 6 %s %s %s' % (inputDir, alphabet, sampling, enableRNASampling)
        os.system(cmd)

        jobDoneFn = os.path.join(inputDir, 'job Done!!')
        while not os.path.exists(jobDoneFn):
            time.sleep(10)

#end_func

def main():
    if PRINT_INTO_LOG:
        sys.stdout = open(LOG_FN, 'w')
        
    # RunAllCellLine()
    RunAllRNA()

    if PRINT_INTO_LOG:
        sys.stdout.close()
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
