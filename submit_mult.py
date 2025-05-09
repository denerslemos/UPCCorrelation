#!/usr/bin/env python

'''Add important stuff'''
import os.path
import optparse

outputfolder = "/eos/user/d/ddesouza/UPChistos"
cmsswrepo = "/afs/cern.ch/work/d/ddesouza/UIC/SPRACE/CMSSW_13_0_5/src/"
pwdrepo = "/afs/cern.ch/work/d/ddesouza/UIC/SPRACE/CMSSW_13_0_5/src/UPCCorrelation/"
systematics = 0
inPut = 'PbPb_UPC_syst_'+str(systematics)

os.system("mkdir -p cond")
os.system("mkdir -p "+str(outputfolder)+"/"+str(inPut))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/PF_Run23_UPC_nonemptyfiles -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_syst"+str(systematics)+"_job_ -f tomorrow -c 2 -n 29 -s PbPb_UPC_syst"+str(systematics)+" -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
