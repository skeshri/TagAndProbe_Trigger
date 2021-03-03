name = 'DYEE_Run3_v1'

dataset = {
   #'DY2021' : '/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v2/MINIAODSIM',
   'DY2021_incl' : '/DYJets_incl_MLL-50_TuneCP5_14TeV-madgraphMLM-pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v1/MINIAODSIM',
  # 'DY2023' : '/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v2/MINIAODSIM', 
   'DY2023_incl' : '/DYJets_incl_MLL-50_TuneCP5_14TeV-madgraphMLM-pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v1/MINIAODSIM', 
   }


#nevents = -1 
#lumisPerJob = {
#   'Run2017B':        100,
#   'Run2017C':        100,
#   'Run2017D':        100,
#   'Run2017E':        100,
#   'Run2017F':        100,
#   }

listOfSamples = [
   #'DY2021',        
   'DY2021_incl',        
  # 'DY2023',        
   'DY2023_incl',        
   ]


if __name__ == '__main__':

   from CRABClient.UserUtilities import config
   config = config()

   from CRABAPI.RawCommand import crabCommand
   from multiprocessing import Process

   def submit(config):
       res = crabCommand('submit', config = config)

   config.General.workArea = 'crab_'+name
   config.General.transferLogs = False
   config.General.transferOutputs = True
   config.JobType.allowUndistributedCMSSW = True
   config.JobType.pluginName = 'Analysis'
   config.JobType.psetName = 'runNtupler.py'
   config.JobType.outputFiles = ['TnP_ntuple.root']

   config.Data.inputDBS = 'global'
   config.Data.splitting = 'Automatic'
#   config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
   config.Data.publication = False
   config.Data.totalUnits = -1
   config.Data.outLFNDirBase = '/store/group/phys_egamma/arun/TriggerEff_Run3/' + name

   config.Site.storageSite = 'T2_CH_CERN'
 #  config.Site.blacklist = ['T2_BR_SPRACE', 'T2_US_Wisconsin', 'T1_RU_JINR', 'T2_RU_JINR', 'T2_EE_Estonia']

   listOfSamples.reverse()
   for sample in listOfSamples:

      config.General.requestName = sample
      config.Data.inputDataset = dataset[sample]
#      config.Data.unitsPerJob = lumisPerJob[sample]
      config.Data.outputDatasetTag = sample
      p = Process(target=submit, args=(config,))
      p.start()
      p.join()
