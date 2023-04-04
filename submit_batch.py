#! /usr/bin/env python
#
# example: python submit_batch.py -c -N -1 -n 1 -p testtnp myfiles.txt
# -c = just create and do not submit (remove it to submit)
# this is to write 1 job per file in the dataset (-N = run on all the events in a file; -n = 1job/file)
#
# to copy the output on EOS: python submit_batch.py -c -N -1 -n 1 -p testtnp --eos=ok myfiles.txt
#
# if files are accessed via Grid, before you need:
# voms-proxy-init --voms cms --valid 168:00
# mv /tmp/x509up_u7291 ~/X509_USER_PROXY
# setenv X509_USER_PROXY ~/X509_USER_PROXY
# echo $X509_USER_PROXY 
#
# myfiles.txt   contains the nanoAODs

import os
import sys
import re
import time
import commands
import optparse
import datetime

def makeCondorFile(jobdir, srcFiles, options):
    dummy_exec = open(jobdir+'/dummy_exec.sh','w')
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('bash $*\n')
    dummy_exec.close()
     
    condor_file_name = jobdir+'/condor_submit.condor'
    condor_file = open(condor_file_name,'w')
    condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = true
Log        = {jd}/$(ProcId).log
Output     = {jd}/$(ProcId).out
Error      = {jd}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
request_memory = 1000
+MaxRuntime = {rt}\n
'''.format(de=os.path.abspath(dummy_exec.name), jd=os.path.abspath(jobdir), rt=int(options.runtime*3600), here=os.environ['PWD'] ) )
    for sf in srcFiles:
        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(sf)))
    condor_file.close()
    return condor_file_name

def main():
#######################################
### usage  
#######################################
    usage = '''usage: %prog [opts] dataset'''
    parser = optparse.OptionParser(usage=usage)
    now = datetime.datetime.now()
    defaultoutputdir='job_'+str(now.year)+str(now.month)+str(now.day)+"_"+str(now.hour)+str(now.minute)+str(now.second)
        
    parser.add_option('-q', '--queue',       action='store',     dest='queue',       help='run in batch in queue specified as option (default -q 8nh)', default='8nh')
    parser.add_option('-n', '--nfileperjob', action='store',     dest='nfileperjob', help='split the jobs with n files read/batch job'                , default=1,   type='int')
    parser.add_option('-p', '--prefix',      action='store',     dest='prefix',      help='the prefix to be added to the output'                      , default=defaultoutputdir)
    parser.add_option('-a', '--application', action='store',     dest='application', help='the executable to be run'                                  , default='/afs/cern.ch/work/c/crovelli/BPhys/BParkAnalysisToolsRun2Pisa/BParkApp')
    parser.add_option('-c', '--create',      action='store_true',dest='create',      help='create only the jobs, do not submit them'                  , default=False)
    parser.add_option('-t', '--testnjobs',   action='store',     dest='testnjobs',   help='submit only the first n jobs'                              , default=1000000, type='int')
    parser.add_option('-N', '--neventsjob', action='store',      dest='neventsjob',  help='split the jobs with n events  / batch job'                 , default=200,   type='int')
    parser.add_option('-T', '--eventsperfile', action='store',   dest='eventsperfile',  help='number of events per input file'                        , default=-1,   type='int')
    parser.add_option('-r', '--runtime',     action='store',     dest='runtime',     help='New runtime for condor resubmission in hours. default None: will take the original one.', default=4        , type=int);
    parser.add_option('--eos',               action='store',     dest='eos',         help='copy the output in the specified EOS path'                 , default='')
    parser.add_option('--scheduler',         action='store',     dest='scheduler',   help='select the batch scheduler (lsf,condor). Default=condor'   , default='condor')
    (opt, args) = parser.parse_args()
    print args

    if len(args) != 1:
        print len(args)
        print usage
        sys.exit(1)
    inputlist = args[0]

    output = os.path.splitext(os.path.basename(inputlist))[0]
    print "output = "+output

    print "the outputs will be in the directory: "+opt.prefix

    diskoutputdir = ''
    diskoutputmain = diskoutputdir+"/"+opt.prefix+"/"+output

    jobdir = opt.prefix+"/"+output
    logdir = jobdir+"/log/"
    os.system("mkdir -p "+jobdir)
    os.system("mkdir -p "+logdir)
    os.system("mkdir -p "+jobdir+"/src/")
    os.system("mkdir -p "+jobdir+"/cfg/")

    outputroot = diskoutputmain+"/root/"

    #look for the current directory
    #######################################
    pwd = os.environ['PWD']
    scramarch = os.environ['SCRAM_ARCH']
    #######################################
    inputListfile=open(inputlist)
    inputfiles = inputListfile.readlines()
    ijob=0
    jobdir = opt.prefix+"/"+output

    srcfiles = []
    while (len(inputfiles) > 0):
        L = []
        for line in range(min(opt.nfileperjob,len(inputfiles))):
            ntpfile = inputfiles.pop()
            ntpfile = ntpfile.rstrip('\n')
            ####ntpfile = re.sub(r'/eos/cms','',ntpfile.rstrip())
            if ntpfile != '':
                L.append("\'"+ntpfile+"\',\n")
            print ntpfile    

            # prepare the txt with root files
            icfgfilename = pwd+"/"+opt.prefix+"/"+output+"/cfg/tnp_"+str(ijob)+".txt"
            icfgfile = open(icfgfilename,'w')
            icfgfile.write(ntpfile)
            icfgfile.close()

            # prepare the script to run
            #rootoutputfile = output+'_'+str(ijob)+'.root'
            rootoutputfile = output+'_'+str(ijob)
            print "rootoutputfile = "+rootoutputfile
            outputname = jobdir+"/src/submit_"+str(ijob)+".src"
            outputfile = open(outputname,'w')
            outputfile.write('#!/bin/bash\n')
            outputfile.write('cd '+pwd+'\n')
            if opt.scheduler=='condor':
                rootoutputfile = '/tmp/'+rootoutputfile
            outputfile.write('echo $PWD\n')
            #outputfile.write(opt.application+' '+icfgfilename+' \n')
            outputfile.write(opt.application+' '+icfgfilename+' '+rootoutputfile+' \n')
            if(opt.eos!=''): 
                outputfile.write('cp '+rootoutputfile+'.root /eos/cms/store/user/crovelli/ParkingForPisa/'+output+'_'+str(ijob)+'.root\n')
                outputfile.write('rm '+rootoutputfile+'.root')
            outputfile.close()
            logfile = logdir+output+"_"+str(ijob)+".log"
            scriptfile = pwd+"/"+outputname
            if opt.scheduler=='condor':
                srcfiles.append(outputname)
            else:
                print "ERROR. Only Condor scheduler available"
                sys.exit(1)
            ijob = ijob+1
            if(ijob==opt.testnjobs): break
            if (opt.eventsperfile == -1): break

    if opt.scheduler=='condor':
        cf = makeCondorFile(jobdir,srcfiles,opt)
        subcmd = 'condor_submit {rf} '.format(rf = cf)
        if opt.create:
            print 'running dry, printing the commands...'
            print subcmd
        else:
            print 'submitting for real...'
            os.system(subcmd)

if __name__ == "__main__":
        main()

