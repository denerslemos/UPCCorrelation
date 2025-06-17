# HBT correlation code using UPC's
This code is going to be used to measure the HBT correlations in 1 and 3D for charged hadrons in pO, OO and NeNe collisions with 2025 dataset.
## Instructions
The code was updated to work on LXPLUS9 in order to use more resources from the cluster in HTCondor. To login in LXPLUS9 use:
```
ssh username@lxplus9.cern.ch
```
Once logged into LXPLUS machines setup the CMSSW version (to be updated) as follows:
```
cmsrel CMSSW_13_0_5
cd CMSSW_13_0_5/src/
cmsenv
```
and work inside of the ```src``` folder.
The code was created to generate histograms as function of multiplicity.
- To download the code use:
```
git clone https://github.com/denerslemos/UPCCorrelation.git
```
You can run a test using (no mixing, no 3D):
```
root -l -b -q "correlation_UPC.C(\"inputdataset/inputtestx.txt\", \"output\", 0, 1, 10, 5, 1.0, 1, 0, 0)" &> outlog.txt &
```
Adding mixing (no 3D):
```
root -l -b -q "correlation_UPC.C(\"inputdataset/inputtestx.txt\", \"output\", 0, 0, 10, 5, 1.0, 1, 0, 0)" &> outlog.txt &
```
Adding 3D (no mixing):
```
root -l -b -q "correlation_UPC.C(\"inputdataset/inputtestx.txt\", \"output\", 0, 1, 10, 5, 1.0, 0, 0, 0)" &> outlog.txt &
```
Adding all:
```
root -l -b -q "correlation_UPC.C(\"inputdataset/inputtestx.txt\", \"output\", 0, 0, 10, 5, 1.0, 0, 0, 0)" &> outlog.txt &
```
Test will run over ~27M events. Note that track efficiency correction is not added. Also track information is not available yet as well as other track informations. Also event filters must be checked/studies.

If you wanna submit a condor job for the entire dataset (please edit [lines 7 to 9](https://github.com/denerslemos/UPCCorrelation/blob/main/submit_mult.py#L7-L9) ) in order , use
```
python3 submit_mult.py
```
This will submit 29 jobs for a single systematics. It can be changed later (loop over systematics, be carefull reading same file multiple times). Also does not include 3D studies.
