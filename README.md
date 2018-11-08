# TagAndProbe_Trigger

#https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_92X_samples_R

cmsrel CMSSW_9_4_0
cd CMSSW_9_4_0/src
cmsenv
git cms-init
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
scram b -j 9


# Add the area containing the MVA weights (from cms-data, to appear in “external”).
# Note: the “external” area appears after “scram build” is run at least once, as above
#
cd $CMSSW_BASE/external
# below, you may have a different architecture, this is just one example from lxplus
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
# Go back to the src/
cd $CMSSW_BASE/src


