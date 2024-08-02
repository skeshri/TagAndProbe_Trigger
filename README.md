# TagAndProbe_Trigger 
### Part 1: Ntupler

```
cmsrel CMSSW_14_0_12  
cd CMSSW_14_0_12/src  
cmsenv  
git cms-init  
git cms-addpkg RecoEgamma/ElectronIdentification
git clone https://github.com/skeshri/TagAndProbe_Trigger.git  
scram b -j 4   
```
*For test run*
```
cd $CMSSW_BASE/src/TagAndProbe_Trigger/NtupleProducer/test     
cmsRun runNtupler.py  
```
### More instructions to come for Efficiency calculations

