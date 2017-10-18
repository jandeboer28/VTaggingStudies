#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re

#locatin /nfs/dust/cms/user/clseitz/DarkMatterMC/LHE_Grid_Scalar_Jul25/DMScalar_ttbar01j_Mphi100_Mchi20_g1_44965/Events/run_01/
def createJobs(f , jobs, i, EventsPerJob, outfolder):
    cmd = 'cmsRun config_doSubstr.py ' + f + ' '+ str(EventsPerJob) + ' ' + str(i) + ' ' + outfolder + '\n'
    print cmd
    jobs.write(cmd)
    return 1

def submitJobs(jobList, nchunks, outfolder, batchSystem):
    print 'Reading joblist'
    jobListName = jobList
    print jobList
#    subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobListName)
    subCmd = 'qsub -t 1-%s -o %s/logs/ %s %s' %(nchunks,outfolder,batchSystem,jobListName)
    print 'Going to submit', nchunks, 'jobs with', subCmd
    os.system(subCmd)

    return 1

def getNEvents(f):
	nEvents = int(sys.argv[3])
	print nEvents
	return nEvents
if __name__ == "__main__":

    outfolder = "Output"

    #should probably rewrite with option parser
    batchSystem = 'psibatch_runner.sh'
    if len(sys.argv) > 1:
        if sys.argv[1].find("S")!=-1: pattern = "WprimeToWZToWhadZhad"
        if sys.argv[1].find("B")!=-1: pattern = "QCD_Pt-15to7000"
        print 'Location of input files', pattern
    else:
        print "No location given, give folder with files"
        exit(0)

    if len(sys.argv) > 2:
        outfolder = sys.argv[2]
        print 'Output goes here: ', outfolder                                                                                                                                    
    else: 
        print "Using default output folder: ", outfolder                                                                                                                            

    try: os.stat(outfolder) 
    except: os.mkdir(outfolder)

    try: os.stat(outfolder+'/logs/') 
    except: os.mkdir(outfolder+'/logs/')
    
#    pattern = "datacardsABCD_2p1bins_fullscan2"
    filelist = glob.glob('/pnfs/psi.ch/cms/trivcat/store/mc/RunIISpring16MiniAODv2/'+pattern+'*/MINIAODSIM/*/*/*.root')
    print filelist
    jobList = 'joblist%s.txt'%pattern
    jobs = open(jobList, 'w')
    nChunks = 0
    for i,f in enumerate(filelist):
        print f
        nEvents = getNEvents(f)
        EventsPerJob = min(nEvents,10000)
        nJobs = nEvents / EventsPerJob
        createJobs(f,jobs,i,EventsPerJob, outfolder)
        nChunks = nChunks+1

    jobs.close()
    submit = raw_input("Do you also want to submit the jobs to the batch system? [y/n] ")
    if submit == 'y' or submit=='Y':
        submitJobs(jobList,nChunks, outfolder, batchSystem)
    else:
        print "Not submitting jobs"


    
