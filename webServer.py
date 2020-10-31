#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     :
Author          : Tian Gao
CreationDate    : 2017/5/7
Description     :
    command line : python webServer.py -port=8000


"""
import tornado.ioloop
import tornado.web
import tornado.httpserver
import tornado.options
import tornado.httpserver
import tornado.web
import tornado.gen
import tornado.httpclient
import tornado.concurrent
from tornado.options import define, options
import os
import uuid
import shutil
import re
import pickle
import time
import BioinfoComm
import smtplib
from email.mime.text import MIMEText
from email.header import Header


REFRESH_INTERVAL = 10 # refresh waiting page for XX seconds
SAMPLING_FREQ = 100
ALPHABET = 'ACGU'
PORT = 7001

SHOWING_URL = 'http://sbbi-panda.unl.edu'
BASE_URL = 'http://cse-jcui-08.unl.edu'
# WORKING_DIR = '/home/tgao/GraphBasedMotifFinding'
WORKING_DIR = '/home/tgao/GraphBasedMotifFinding_master'
OUTPUT_BASE_DIR = '/home/tgao/GraphBasedMotifFinding/output/http'
HTML_BASE_DIR = '/var/www/html/MDS2/tmpRlt'
# HTML_BASE_DIR = '/var/www/html/GraphBasedMotifFinding09/tmpRlt'
HTTP_BASE_DIR = 'http://sbbi-panda.unl.edu/MDS2/tmpRlt'
# HTTP_BASE_DIR = 'http://sbbi-panda.unl.edu/GraphBasedMotifFinding09/tmpRlt'
EXISTING_RLT_DIR = '/var/www/html/MDS2/existingRlt'
# EXISTING_RLT_DIR = '/var/www/html/GraphBasedMotifFinding09/existingRlt'

TEMPLATES_DIR = 'http/templates'
STATIC_DIR = 'http/static'

# all the global variables are stored in this dict
OVERALL_DICT = {}

define("port", default=PORT, help="run on the given port", type=int)

# mail config
MAIL_HOST = "cse-mail.unl.edu"
MAIL_USER = "***"
MAIL_PWD = "***"
SENDER = 'bioinfo@cse.unl.edu'

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.redirect("/input")
    #end_func
#end_class

class InputHandler(tornado.web.RequestHandler):
    def get(self):
        self.jobNm = str(uuid.uuid4())
        self.render("input.html", )
    #end_func
#end_class

class RequestRltHandler(tornado.web.RequestHandler):
    def get(self):
        self.render("request.html", )
    #end_func
#end_class

class SendingDownloadHandler(tornado.web.RequestHandler):
    """
    send an email after the user request software
    """
    def post(self):
        # get arguments
        name = self.get_argument('nm')
        pos = self.get_argument('pos', None)
        piNm = self.get_argument('PI_Nm')
        organ = self.get_argument('organ')
        email = self.get_argument('email')

        print name, pos, piNm, organ, email

        # send email
        userInfoFn = os.path.join(OUTPUT_BASE_DIR, 'downloadUserInfo.txt')
        with open(userInfoFn, 'a') as userInfoFileobj:
            outputLine = '\t'.join([name, pos, piNm, organ, email])
            userInfoFileobj.write('%s\n' % outputLine)


        if email:
            receivers = [email]
            content = 'Please visit the following link to download the software:\n%s/MDS2/MDS2.tar' % SHOWING_URL
            subject = 'MDS2 downloading link'
            message = MIMEText(content, 'plain', 'utf-8')
            message['Subject'] = Header(subject, 'utf-8')

            try:
                smtpObj = smtplib.SMTP()
                smtpObj.connect(MAIL_HOST, 25)  # 25 is PORT
                smtpObj.login(MAIL_USER, MAIL_PWD)
                smtpObj.sendmail(SENDER, receivers, message.as_string())
                self.write("email is sent successfully! Please check your email.<br>")
                self.write('<a href="%s:%s/input">back to MDS2</a>' % (BASE_URL, PORT))
            except smtplib.SMTPException:
                print "fail to send the email."



#end_func

class AllCellLineHandler(tornado.web.RequestHandler):
    """
    give a webpage to access to all known the cell line
    """
    def get(self):
        self.existingRlt = {'vesiclepedia': [('hsa_Breast_cancer_cells', 'Human breast cancer cells'), ('hsa_B_cells', 'Human B cells'), ('hsa_Lung_cancer_cells', 'Human lung cancer cells'), ('hsa_Serum', 'Human serum cells'), ('hsa_T_cells', 'Human T cells'), ('hsa_Urine', 'Human urine')],
                            'evpedia': [('hsa_Caput_epithelial_cell_input', 'Human caput epithelial cells'), ('hsa_Caput_luminal_fluid_input', 'Human caput luminal fluid'), ('hsa_Cauda_epithelial_cell_input', 'Human cauda epithelial cells'), ('hsa_Cauda_luminal_fluid_input', 'Human cauda luminal fluid'), ('hsa_Colorectal_cancer_cell_SW480__input', 'Human colorectal cancer cell - SW480'), ('hsa_Lung_cancer_cell_DMS563__input', 'Human lung cancer cell - DMS563'), ('hsa_Lung_cancer_cell_NCI_H69__input', 'Human lung cancer cell - NCI H69'), ('hsa_Seminal_plasma_input', 'Human seminal plasma')],
                            'exocarta': [('hsa_B_cells', 'Human B cells'), ('hsa_Colon_carcinoma_cells', 'Human colon carcinoma cells'), ('hsa_Colorectal_cancer_cells', 'Human colorectal cancer cells'), ('hsa_Endothelial_cells', 'Human endothelial cells'), ('hsa_Plasma', 'Human plasma'), ('hsa_Serum', 'Human plasma'), ('hsa_T_cells', 'Human T cells')]
                            }
        self.website = {
            'vesiclepedia': "http://www.microvesicles.org/",
            'evpedia': "http://student4.postech.ac.kr/evpedia2_xe/xe/",
            'exocarta': "http://www.exocarta.org/",
        }
        self.render("allCellLine.html", )
    #end_func
#end_class

class inputErrHandler(tornado.web.RequestHandler):
    """
    give a webpage to show that the input format is incorrect
    """
    def get(self):
        self.render("inputErr.html", )
    #end_func
#end_class

class PreprocessHandler(tornado.web.RequestHandler):
    """
    generate necessary files and submit a job after click 'submit' button in the input web page
    """
    def post(self, jobNm):
        outputDir = os.path.join(OUTPUT_BASE_DIR, jobNm)
        errLog = os.path.join(outputDir, 'err.log')
        inputFn = os.path.join(outputDir, 'userInput.fa')
        negSeqFn = os.path.join(outputDir, 'negSeq.fa')
        inputInfoFn = os.path.join(outputDir, 'inputInfo.txt')
        emailFn = os.path.join(outputDir, 'email.txt')

        # get all the arguments in input page
        userEmail = self.get_argument('userEmail', None)
        enableNegSeq = True if self.get_arguments('negSeqCheck') else False
        enableMiRBG = True if self.get_arguments('BGmiRCheck') else False
        covPrecThres = self.get_argument('covThres', None)

        print enableNegSeq
        print enableMiRBG
        print covPrecThres

        # create folder
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        print outputDir

        # get input and write into a local file
        targetFn = ''
        inputType = self.get_argument('inputType', None)
        with open(inputFn, 'w') as inputObj:
            if inputType == 'typeIn': # type in sequences
                content = self.get_argument('typeInSeq', None).replace('\r\n', '\n')
                inputObj.write(content.upper())
                targetFn = 'typeInSeq.fa'
            elif inputType == 'upload': # upload sequences
                if self.request.files.get('uploadSeq', None):
                    uploadFile = self.request.files['uploadSeq'][0]
                    inputObj.write(uploadFile['body'])
                    targetFn = uploadFile['filename']

        # get negative sequences input and write into a local file
        if enableNegSeq:
            inputType = self.get_argument('inputTypeNeg', None)
            with open(negSeqFn, 'w') as negSeqFileobj:
                if inputType == 'typeIn': # type in sequences
                    content = self.get_argument('typeInSeqNeg', None).replace('\r\n', '\n')
                    negSeqFileobj.write(content.upper())
                elif inputType == 'upload': # upload sequences
                    if self.request.files.get('uploadSeqNeg', None):
                        uploadFile = self.request.files['uploadSeqNeg'][0]
                        negSeqFileobj.write(uploadFile['body'])

        seqLis, seqCnt, titleLis, minSeqLen, maxSeqLen = BioinfoComm.loadOnlineInput(inputFn)
        seqLis = [item for item in seqLis if item]
        titleLis = [item for item in titleLis if item]
        inputCharSet = reduce(lambda x, y: set(x) | set(y), seqLis)
        if not IsInputValid(seqLis, seqCnt, titleLis, minSeqLen, maxSeqLen, inputCharSet):self.redirect('/inputErr')

        # write the info about inputted seq into a file
        with open(inputInfoFn, 'w') as inputInfoFileobj:
            inputInfoFileobj.write('%s\n%s\n' % (targetFn, '\t'.join(map(lambda x:str(x), [seqCnt, minSeqLen, maxSeqLen]))))

        # write user email
        with open(emailFn, 'w') as emailFileobj:
            if userEmail:emailFileobj.write(userEmail)

        # generate shell script
        runFn = os.path.join(outputDir, 'run.sh')
        alphabet = ''.join(inputCharSet)
        shCmd = '#!/usr/bin/env bash\npython GraphBasedMotifFinding.py %s 6 %s %s %s %s 2 %s' % (outputDir, alphabet, SAMPLING_FREQ, enableMiRBG, enableNegSeq, covPrecThres)
        print shCmd
        with open(runFn, 'w') as runFile:
            runFile.write(shCmd)

        RunJob(outputDir, jobNm, errLog, WORKING_DIR, runFn)
        self.redirect('/processing/%s' % jobNm)
    #end_func
#end_class

class ProcessingHandler(tornado.web.RequestHandler):
    def get(self, jobNm):
        self.refreshInterval = REFRESH_INTERVAL
        self.jobNm = jobNm
        self.rltUrl = '%s:%s/result/%s' % (BASE_URL, PORT, jobNm)

        curOutputDir = os.path.join(OUTPUT_BASE_DIR, jobNm)
        inputInfoFn = os.path.join(curOutputDir, 'inputInfo.txt')
        userInput = os.path.join(curOutputDir, 'userInput.fa')

        with open(inputInfoFn) as inputInfoFileobj:
            inputFnInfo = [inputInfoFileobj.readline().strip()]
            inputFnInfo.append(inputInfoFileobj.readline().strip().split('\t'))

        self.inputFnInfoStr = '%s sequences with length from %s to %s in file %s' % (inputFnInfo[1][0], inputFnInfo[1][1], inputFnInfo[1][2], inputFnInfo[0])
        self.submitTime = time.ctime(os.stat(userInput).st_ctime)
        self.expireTime = time.ctime(os.stat(userInput).st_ctime + 259200) # 3 days to expire

        if not IsJobDone(jobNm):
            self.render("processing.html")
        else:

            self.redirect("/postprocessing/%s" % jobNm)
#end_class

class PostprocessHandler(tornado.web.RequestHandler):
    # copy necessary files from output dir to html dir
    def get(self, jobNm):
        # input
        outputDir = os.path.join(OUTPUT_BASE_DIR, jobNm)
        userInput = os.path.join(outputDir, 'userInput.fa')
        rltFn = os.path.join(outputDir, 'finalRlt.txt')
        inputInfoFn = os.path.join(outputDir, 'inputInfo.txt')
        pattern2kmerSetFn = os.path.join(outputDir, 'pattern2kmerSet.txt')
        rltHtml = os.path.join(WORKING_DIR, 'http', 'templates', 'result.html')
        emailFn = os.path.join(outputDir, 'email.txt')
        # visHtml = os.path.join(WORKING_DIR, 'http', 'templates', 'selectAndVis.html')

        # output
        curHtmlDir = os.path.join(HTML_BASE_DIR, jobNm)
        locFn = os.path.join(curHtmlDir, 'motifLocation.txt')

        # copy files
        if not os.path.exists(curHtmlDir): os.makedirs(curHtmlDir)
        shutil.copy(userInput, curHtmlDir)
        shutil.copy(rltFn, curHtmlDir)
        shutil.copy(inputInfoFn, curHtmlDir)
        shutil.copy(pattern2kmerSetFn, curHtmlDir)
        shutil.copy(rltHtml, curHtmlDir)
        shutil.copy(emailFn, curHtmlDir)
        # shutil.copy(visHtml, curHtmlDir)

        outputPngDir = os.path.join(outputDir, 'png')
        displayPngDir = os.path.join(curHtmlDir, 'png')

        # if not os.path.exists(displayPngDir): os.makedirs(displayPngDir)
        shutil.copytree(outputPngDir, displayPngDir)

        # finding location
        seqLis = loadInputSeq(userInput)
        patternInfoLis = loadRlt(rltFn, jobNm)
        pattern2locInfo, maxSeqLen = findLoc(seqLis, patternInfoLis)

        with open(locFn, 'w') as locFile:
            pickle.dump((pattern2locInfo, maxSeqLen), locFile)

        # send email
        emailFn = os.path.join(curHtmlDir, 'email.txt')
        with open(emailFn) as emailFileobj:
            email = emailFileobj.readline().strip()

        if email:
            receivers = [email]
            content = 'Please visit the following link to see the result:\n%s:%s/result/%s' % (BASE_URL, PORT, jobNm)
            subject = 'Your job %s is done' % jobNm
            message = MIMEText(content, 'plain', 'utf-8')
            message['Subject'] = Header(subject, 'utf-8')

            try:
                smtpObj = smtplib.SMTP()
                smtpObj.connect(MAIL_HOST, 25)  # 25 is PORT
                smtpObj.login(MAIL_USER, MAIL_PWD)
                smtpObj.sendmail(SENDER, receivers, message.as_string())
                print "email is sent successfully!"
            except smtplib.SMTPException:
                print "fail to send the email."


        self.redirect("/result/%s" % jobNm)
    #end_func
#end_class

class ResultHandler(tornado.web.RequestHandler):
    def get(self, jobNm):
        curHtmlDir = os.path.join(HTML_BASE_DIR, jobNm)
        rltFn = os.path.join(curHtmlDir, 'finalRlt.txt')
        rltHtml = os.path.join(curHtmlDir, 'result.html')
        locFn = os.path.join(curHtmlDir, 'motifLocation.txt')
        inputInfoFn = os.path.join(curHtmlDir, 'inputInfo.txt')
        pattern2kmerSetFn = os.path.join(curHtmlDir, 'pattern2kmerSet.txt')

        # load match location
        locFile = open(locFn)
        self.locDict, self.maxSeqLen = pickle.load(locFile)
        locFile.close()

        # format sequence name
        for (length, pattern), locInfo in self.locDict.iteritems():
            for i in xrange(len(locInfo)):
                rawSeqNm = locInfo[i][1]
                formattedNm = rawSeqNm.lower().replace('mir', 'miR')
                locInfo[i][1] = formattedNm

        # load pattern info
        self.jobNm = jobNm
        self.patternInfoLis = loadRlt(rltFn, jobNm) # this one is empty if there is no significant motif
        self.motifAllLenLis = map(lambda x:str(x[0]), self.patternInfoLis) if self.patternInfoLis else []
        self.motifDefaultLen = self.motifAllLenLis[0] if self.patternInfoLis else 0
        self.motifAllLenStr = ','.join(self.motifAllLenLis) if self.patternInfoLis else ''
        self.pattern2Idx = {}
        for (patternLen, pattern), values in self.locDict.iteritems():
            self.pattern2Idx[pattern] = str(len(self.pattern2Idx))

        # load input file info
        with open(inputInfoFn) as inputInfoFileobj:
            self.inputFnInfo = [inputInfoFileobj.readline().strip()]
            # self.inputFnInfo = ['typeInSeq.fa']
            self.inputFnInfo.append(inputInfoFileobj.readline().strip().split('\t'))
        if not self.inputFnInfo[1][0]:
            self.inputFnInfo = ['typeInSeq.fa', self.inputFnInfo[0]]

        # load pattern2kmerSet
        pattern2kmerSetFileobj = open(pattern2kmerSetFn)
        self.pattern2kmerSet = pickle.load(pattern2kmerSetFileobj)
        pattern2kmerSetFileobj.close()

        self.render(rltHtml)
    #end_func
#end_class

class existingRltHandler(tornado.web.RequestHandler):
    # copy necessary files from existing cell line dir to html dir
    def get(self, jobNm):
        # input
        outputDir = os.path.join(EXISTING_RLT_DIR, jobNm)
        userInput = os.path.join(outputDir, 'userInput.fa')
        rltFn = os.path.join(outputDir, 'finalRlt.txt')
        inputInfoFn = os.path.join(outputDir, 'inputInfo.txt')
        pattern2kmerSetFn = os.path.join(outputDir, 'pattern2kmerSet.txt')
        rltHtml = os.path.join(WORKING_DIR, 'http', 'templates', 'result.html')

        # generate input info
        seqLis, seqCnt, titleLis, minSeqLen, maxSeqLen = BioinfoComm.loadOnlineInput(userInput)
        with open(inputInfoFn, 'w') as inputInfoFileobj:
            inputInfoFileobj.write('%s\n%s\n' % ('typeInSeq.fa', '\t'.join(map(lambda x:str(x), [seqCnt, minSeqLen, maxSeqLen]))))

        # output
        curHtmlDir = os.path.join(HTML_BASE_DIR, jobNm)
        locFn = os.path.join(curHtmlDir, 'motifLocation.txt')

        # copy files
        if not os.path.exists(curHtmlDir): os.makedirs(curHtmlDir)
        shutil.copy(userInput, curHtmlDir)
        shutil.copy(rltFn, curHtmlDir)
        shutil.copy(inputInfoFn, curHtmlDir)
        shutil.copy(pattern2kmerSetFn, curHtmlDir)
        shutil.copy(rltHtml, curHtmlDir)

        outputPngDir = os.path.join(outputDir, 'png')
        displayPngDir = os.path.join(curHtmlDir, 'png')

        if not os.path.exists(outputPngDir): os.makedirs(displayPngDir)
        if not os.path.exists(displayPngDir): shutil.copytree(outputPngDir, displayPngDir)



        # finding location
        seqLis = loadInputSeq(userInput)
        patternInfoLis = loadRlt(rltFn, jobNm)
        pattern2locInfo, maxSeqLen = findLoc(seqLis, patternInfoLis)

        with open(locFn, 'w') as locFile:
            pickle.dump((pattern2locInfo, maxSeqLen), locFile)

        self.redirect("/result/%s" % jobNm)
    #end_func
#end_class

def findLoc(seqLis, patternInfoLis):
    """
    get a dict including location information
    {
        pattern: [(seqTitle, [matchIdx1, matchIdx2])]
    }
    """
    pattern2locInfo = {}
    maxSeqLen = 0
    for motifLength, items in patternInfoLis:
        for pattern, logoFn, coverage, IC, pvalueRNA, pvalueMiR in items:
            patternLen = GetPatternLength(pattern)
            pattern2locInfo[(patternLen, pattern)] = []
            for seqID, (seqTitle, seq) in enumerate(seqLis):
                seqLen = len(seq)
                if seqLen > maxSeqLen: maxSeqLen = seqLen
                matchSegSet = set(re.findall(pattern, seq))
                idxLis = []
                for matchSeg in matchSegSet:
                    while matchSeg in seq:
                        matchIdx = seq.index(matchSeg)
                        idxLis.append(matchIdx)
                        seq = seq[matchIdx + len(matchSeg):]
                pattern2locInfo[(patternLen, pattern)].append([seqID + 1, seqTitle, seqLen, idxLis])
    return pattern2locInfo, maxSeqLen
#end_func

def loadRlt(rltFn, jobNm):
    patternInfoLis = []

    with open(rltFn) as rltFile:
        for line in rltFile:
            line = line.strip()
            if line[:7] == 'pattern':  # for different length
                length = int(line.split(':')[1])
                patternInfoLis.append((length, []))
                continue
            items = line.split('\t')
            pattern = items[0]
            coverage = items[1]
            IC = items[2]
            pvalueRNA = items[3]
            pvalueMiR = items[4]
            logoFn = os.path.join(HTTP_BASE_DIR, jobNm, 'png', str(length), 'logo-%s.png' % pattern)
            patternInfoLis[-1][1].append((pattern, logoFn, coverage, IC, pvalueRNA, pvalueMiR))
    return patternInfoLis
#end_func

def loadInputSeq(fastaFn):
    seqLis = []
    seqTitle = ''

    with open(fastaFn) as fastaF:
        for line in fastaF:
            line = line.strip()
            if line[0] == '>':
                seqTitle = line[1:]
            elif line:
                seqLis.append((seqTitle, line))
    return seqLis
# end_func

def IsJobDone(jobNm):
    """
    a job is done when it doesn't exist
    """
    jobDoneFn = os.path.join(OUTPUT_BASE_DIR, jobNm, 'job Done!!')
    return os.path.exists(jobDoneFn)
#end_func

def IsInputValid(seqLis, seqCnt, titleLis, minSeqLen, maxSeqLen, inputCharSet):
    """
    check whether the input file is validated
    """
    if len(seqLis) == len(titleLis) and 0 < seqCnt <= 300 and 8 <= minSeqLen <= 40 and 8 <= maxSeqLen <= 40 and (inputCharSet == {'A', 'C', 'G', 'U'} or inputCharSet == {'A', 'C', 'G', 'T'}):
        return True
    else:
        return False
#end_func

def RunJob(outputDir, jobNm, errLog, WORKING_DIR, runFn):
    cmd = 'sh %s &' % runFn
    os.system(cmd)
    return

    cmd = 'qsub -o %s -N %s -e %s -wd %s %s' % (outputDir, jobNm, errLog, WORKING_DIR, runFn)
    print cmd
    submitInfo = os.popen(cmd).readlines()
    # TODO: handle with error in submitting jobs
    if 'has been submitted' not in submitInfo:
        print 'error'
#end_func

def GetPatternLength(patternStr):
    strLis = []
    while patternStr.find('[') != -1:
        head, _, rest = patternStr.partition('[')
        segPattern, _, patternStr = rest.partition(']')
        strLis.extend(list(head))
        strLis.append(segPattern)
    strLis.extend(list(patternStr))
    return len(strLis)
#end_func

def GetApplication():
    handlers = [
        (r"/", MainHandler),
        (r"/input", InputHandler),
        (r"/allCellLine", AllCellLineHandler),
        (r"/inputErr", inputErrHandler),
        (r"/preprocess/([\S]*)", PreprocessHandler),
        (r"/processing/([\S]*)", ProcessingHandler),
        (r"/postprocessing/([\S]*)", PostprocessHandler),
        (r"/result/([\S]*)", ResultHandler),
        (r"/existingRlt/([\S]*)", existingRltHandler),
        (r"/request", RequestRltHandler),
        (r"/sendingDownload", SendingDownloadHandler),
    ]

    settings = dict(
        template_path=os.path.join(os.path.dirname(__file__), TEMPLATES_DIR),
        static_path=os.path.join(os.path.dirname(__file__), STATIC_DIR),
        debug=True,
    )

    return tornado.web.Application(handlers, **settings)
#end_func

def main():
    tornado.options.parse_command_line()
    app = GetApplication()
    http_server = tornado.httpserver.HTTPServer(app)
    http_server.listen(options.port)
    tornado.ioloop.IOLoop.instance().start()
#end_main

if __name__ == "__main__":
    main()
#end_if
