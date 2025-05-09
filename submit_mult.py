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
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job01.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job01_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job01 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job02.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job02_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job02 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job03.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job03_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job03 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job04.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job04_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job04 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job05.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job05_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job05 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job06.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job06_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job06 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job07.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job07_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job07 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job08.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job08_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job08 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job09.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job09_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job09 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job10.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job10_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job10 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job11.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job11_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job11 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job12.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job12_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job12 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job13.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job13_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job13 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job14.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job14_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job14 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job15.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job15_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job15 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job16.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job16_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job16 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job17.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job17_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job17 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job18.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job18_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job18 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job19.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job19_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job19 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job20.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job20_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job20 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job21.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job21_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job21 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job22.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job22_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job22 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job23.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job23_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job23 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job24.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job24_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job24 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job25.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job25_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job25 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job26.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job26_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job26 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job27.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job27_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job27 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job28.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job28_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job28 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
os.system("python3 HTCondor_submit_multiplicity.py -i inputdataset/forcondor/listfiles_job29.list -o "+str(outputfolder)+"/"+str(inPut)+"/PbPb_UPC_job29_ -f tomorrow -c 2 -n 1 -s PbPb_UPC_job29 -u "+str(systematics)+" -w "+str(cmsswrepo)+" -p "+str(pwdrepo))
