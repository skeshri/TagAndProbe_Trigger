import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 10

# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '') 
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc','') 

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



#import FWCore.Utilities.FileUtils as FileUtils
# 
#inputFilesAOD = cms.untracked.vstring()
#
#
#inputFilesAOD = cms.untracked.vstring( FileUtils.loadListFromFile (os.environ['CMSSW_BASE']+'/src/Efficiency/Analyzer/test/'+'inputFiles.txt') )
#
inputFilesAOD = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/data/Run2017C/SingleElectron/AOD/12Sep2017-v1/70000/80EF56E5-69A6-E711-AB37-48FD8E2824D7.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/data/Run2022C/EGamma/MINIAOD/PromptReco-v1/000/355/862/00000/55aa652c-6cd1-48d9-8020-4f4821e01b24.root'
'root://cms-xrd-global.cern.ch//store/data/Run2023D/ScoutingPFMonitor/MINIAOD/PromptReco-v2/000/370/717/00000/9d9eba7d-91e6-4a91-a6d0-b4eab50bd743.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2022D/EGamma/MINIAOD/PromptReco-v1/000/357/542/00000/0656fded-83e1-4ea9-8023-f3c606c5df26.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_94X_mc2017_realistic_v11_ext1-v1/920000/76190FD4-93A5-E811-AA90-A4BF0107E164.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2017B/MuonEG/MINIAOD/12Sep2017-v1/100000/0423E958-32A3-E711-94AF-1CC1DE1D0AD4.root'
    )
#
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False 
if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    pileupProductName = "addPileupInfo"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "TnP_ntuple.root"
    pileupProductName = "slimmedAddPileupInfo"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff',
		 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_Iso_V1_cff',
		 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff',
]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
#
# Configure the ntupler module
#
process.ntupler = cms.EDAnalyzer('Ntupler',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 pileup   = cms.InputTag( pileupProductName ),
                                 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 genEventInfoProduct = cms.InputTag('generator'),
                                 #
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 muons        = cms.InputTag("muons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 vertices     = cms.InputTag("offlinePrimaryVertices"),
                                 conversions  = cms.InputTag('allConversions'),
                                 triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
                                 triggerSummaryTag = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                 l1EGTag      = cms.InputTag("caloStage2Digis","EGamma","RECO"),
                                 l1MuonTag    = cms.InputTag("gmtStage2Digis","Muon","RECO"),

				 #pathsToSave  = cms.vstring( "HLT_Ele35_WPTight_Gsf_v",
                 #                                            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                 #                                            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
			#				                                 "HLT_IsoMu27_v",
             #                                                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
              #                                               "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
               #                                              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				#				),
				 pathsToSave  = cms.vstring( "DST_HLTMuon_Run3_PFScoutingPixelTracking_v",
                                             "DST_Run3_EG16_EG12_PFScoutingPixelTracking_v",
                                             "DST_Run3_EG30_PFScoutingPixelTracking_v"
                                             "DST_Run3_JetHT_PFScoutingPixelTracking_v",
								),
                            #     filterToMatch= cms.vstring(
                            #      			  "hltEle35noerWPTightGsfTrackIsoFilter",
				  			#  "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
				  			#  "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
				 			#  "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07",
							#  "hltL3fL1DoubleMu155fPreFiltered8",
							#  "hltL3fL1DoubleMu155fFiltered17",
							#  "hltDiMuon178RelTrkIsoFiltered0p4",
							#  "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12",
							#  "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
							#  "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",
							#  "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"
						#				),
                        filterToMatch = cms.vstring(
                            "hltPreDSTHLTMuonRun3PFScoutingPixelTracking",
                            "hltDoubleEG16EG12CaloIdLHEFilter",
                            "hltSingleEG30CaloIdLHEFilter",
                            "hltPreDSTRun3JetHTPFScoutingPixelTracking",
                            ),
				HLTprocess = cms.string("HLT"),

				#
                                 # Objects specific to MiniAOD format
                                #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
                                 muonsMiniAOD = cms.InputTag("slimmedMuons"),
                                 trigger     = cms.InputTag("TriggerResults", "", "HLT"),
			         prescale = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag('slimmedPatTrigger'),
                                 # Effective areas for computing PU correction for isolations
				#effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),
                                 #effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
                                 #effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),
                                 effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Run3_Winter22/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_122X.txt"),
                                 # ID decisions (common to all formats)
                                 #
                                 # all IDs listed below are available given the content of "my_id_modules" defined above.
                                 # only one is exercised for this example.
                                 #
                             
                                  #eleIdMapLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
                                  #eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
                                  #eleIdMapTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
                                  eleIdMapLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-loose"),
                                  eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-medium"),
                                  eleIdMapTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-RunIIIWinter22-V1-tight"),
				  #eleMVA90noIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-RunIIIWinter22-noIso-V1-wp90'),
				  #eleMVA80noIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-RunIIIWinter22-noIso-V1-wp80'),
				  #eleMVALoosenoIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose'),
				  #eleMVA90Iso    =  cms.InputTag('egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp90'),
				  #eleMVA80Iso    =  cms.InputTag('egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp80'),
				  #eleMVALooseIso =  cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose'),
				  #eleMVAValuesMapTokenIso = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values'),		
				  #eleMVAValuesMapTokenNoIso = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values'),	
                                  muInputTag = cms.InputTag("gmtStage2Digis","Muon","RECO"),
                                 isMC = cms.bool(False),
                                 doMuon = cms.bool(False),
                                 doEle = cms.bool(True)

                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )


process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
