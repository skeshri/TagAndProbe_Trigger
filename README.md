# TagAndProbe_Trigger

The ntuple making part works for 2021 and 2023 samples.
Although the macros works for all years.

cmsrel CMSSW_10_6_20  
cd CMSSW_9_4_13/src  
cmsenv  
git clone git@github.com:arunhep/TagAndProbe_Trigger.git
git checkout EGMHLT_run3
scramv1 b

## For Test Run (check the input file and triggers)
cd $CMSSW_BASE/src/TagAndProbe_Trigger/NtupleProducer/test   
<br>  
cmsRun runNtupler.py  
<br>  

## For crab job submission (check the paths etc in the file)

python multicrab_Ele.py

### Macros to produce the plots
cd TagAndProbe_Trigger/TagAndProbeMacros
./run.sh 
