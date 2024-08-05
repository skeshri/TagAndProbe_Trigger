
name = 'Scouting_Run3'

dataset = {
   #'Run2017B' : '/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD',
   #'Run2017C' : '/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD', 
   #'Run2017D' : '/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD', 
   #'Run2017E' : '/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD', 
   #'Run2023D' : '/ScoutingPFMonitor/Run2023D-PromptReco-v2/MINIAOD',
   'Run2023A' : '/ScoutingPFMonitor/Run2023A-PromptReco-v2/MINIAOD',
   'Run2023B': '/ScoutingPFMonitor/Run2023B-PromptReco-v1/MINIAOD',
   'Run2023Cv1': '/ScoutingPFMonitor/Run2023C-PromptReco-v1/MINIAOD',
   'Run2023Cv2': '/ScoutingPFMonitor/Run2023C-PromptReco-v2/MINIAOD',
   'Run2023Cv3': '/ScoutingPFMonitor/Run2023C-PromptReco-v3/MINIAOD',
   'Run2023Cv4': '/ScoutingPFMonitor/Run2023C-PromptReco-v4/MINIAOD',
   'Run2023Dv1' : '/ScoutingPFMonitor/Run2023D-PromptReco-v1/MINIAOD', 
   }


#nevents = -1 
lumisPerJob = {
   #'Run2017B':        100,
   #'Run2017C':        100,
   #'Run2017D':        100,
   #'Run2017E':        100,
   'Run2023A':        200,
   'Run2023B':        200,
   'Run2023Cv1':        200,
   'Run2023Cv2':        200,
   'Run2023Cv3':        200,
   'Run2023Cv4':        200,
   'Run2023Dv1':        200,
   }

listOfSamples = [
   #'Run2017B',        
   #'Run2017C',        
   #'Run2017D',        
   #'Run2017E',        
   'Run2023A',        
   'Run2023B',        
   'Run2023Cv1',        
   'Run2023Cv2',        
   'Run2023Cv3',        
   'Run2023Cv4',        
   'Run2023Dv1',        
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

   config.JobType.pluginName = 'Analysis'
   config.JobType.psetName = 'runNtupler.py'
   config.JobType.outputFiles = ['TnP_ntuple.root']

   config.Data.inputDBS = 'global'
   config.Data.splitting = 'LumiBased'
   #config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
   config.Data.publication = False
   config.Data.totalUnits = -1
   config.Data.outLFNDirBase = '/store/user/skeshri/HLT_Egamma_ntuple/' + name

   config.Site.storageSite = 'T3_CH_CERNBOX'
 #  config.Site.blacklist = ['T2_BR_SPRACE', 'T2_US_Wisconsin', 'T1_RU_JINR', 'T2_RU_JINR', 'T2_EE_Estonia']

   listOfSamples.reverse()
   for sample in listOfSamples:

      config.General.requestName = sample
      config.Data.inputDataset = dataset[sample]
      config.Data.unitsPerJob = lumisPerJob[sample]
      config.Data.outputDatasetTag = sample
      p = Process(target=submit, args=(config,))
      p.start()
      p.join()
