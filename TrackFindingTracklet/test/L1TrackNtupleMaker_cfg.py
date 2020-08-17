############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import os
process = cms.Process("L1TrackNtuple")

############################################################
# edit options here
############################################################

GEOMETRY = "D49"
L1TRKALGO = 'HYBRID'  # L1 tracking algorithm: 'HYBRID' (baseline, 4par fit) or 'HYBRID_DISPLACED' (extended, 5par fit)

WRITE_DATA = False

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.categories.append('Tracklet')
process.MessageLogger.categories.append('L1track')
process.MessageLogger.Tracklet = cms.untracked.PSet(limit = cms.untracked.int32(-1))

if GEOMETRY == "D49": 
    print "using geometry " + GEOMETRY + " (tilted)"
    process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
else:
    print "this is not a valid geometry!!!"

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

# Get list of MC datasets from repo, or specify yourself.

def getTxtFile(txtFileName): 
  return FileUtils.loadListFromFile(os.environ['CMSSW_BASE']+'/src/'+txtFileName)

if GEOMETRY == "D49":
    #inputMC = ["/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/F7BF4AED-51F1-9D47-B86D-6C3DDA134AB9.root"]
    inputMC = ["/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/2D153BE5-611B-3E40-B3D1-E1D5A35A53E6.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/2C3F5156-3A3C-074D-9F89-462EC5E6E058.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/2B390F84-395E-D44B-8E4E-4ACB9F526817.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/26321D7B-7A96-E442-B644-C4BED07AA807.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/25D90FC3-49EF-7448-89F6-935F8D4BB8AC.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/20A54E94-A2D2-FA4A-B5F7-DFF768ACC6F5.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/1F59B52D-BD93-CA4A-BB59-64E00A436EEF.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/1D2905A5-5CB9-BB41-8C59-64427CEE7082.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/19678708-B88E-A547-884E-CD1D924B2D3C.root","/store/relval/CMSSW_11_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200-v1/20000/1818154D-11CB-9841-932A-A48845A430B6.root"]
else:
    print "this is not a valid geometry!!!"
    
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputMC))

process.TFileService = cms.Service("TFileService", fileName = cms.string('TTbar_PU200_'+GEOMETRY+'_test.root'), closeFileFast = cms.untracked.bool(True))
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))


############################################################
# L1 tracking: remake stubs?
############################################################

#process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
#from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
#process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")

#from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
#TTClusterAssociatorFromPixelDigis.digiSimLinks = cms.InputTag("simSiPixelDigis","Tracker")

#process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)
#process.TTClusterStubTruth = cms.Path(process.TrackTriggerAssociatorClustersStubs) 

############################################################
# L1 tracking
############################################################

# Load Quality params to change algorithm and model location
process.load("L1Trigger.TrackTrigger.L1TrackClassifier_cfi")
process.TrackQualityParams.Quality_Algorithm = cms.string("MVA")
#process.TrackQualityParams.ONNXmodel = cms.string("../../TrackTrigger/ML_data/FakeIDGBDT/GBDT_model.onnx")
#process.TrackQualityParams.ONNXmodel = cms.string("../../TrackTrigger/ML_data/GBDT_sklearnxgb_model.onnx")
#process.TrackQualityParams.ONNXInputName = cms.string("feature_input") #for GBDT ONNX model
#process.TrackQualityParams.ONNXInputName = cms.string("float_input") #for sklearn model
#process.TrackQualityParams.in_features = cms.vstring(["eta","bendchi2","nstubs","nlaymiss_interior","chi2rphi","chi2rz"]) #for skl GBDT model
#process.TrackQualityParams.in_features = cms.vstring(["log_chi2","log_bendchi2","log_chi2rphi","log_chi2rz",
#             "nstubs","lay1_hits","lay2_hits","lay3_hits","lay4_hits",
#             "lay5_hits","lay6_hits","disk1_hits","disk2_hits",
#             "disk3_hits","disk4_hits","disk5_hits","rinv","tanl","z0","dtot","ltot"])

process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")

# HYBRID: prompt tracking
if (L1TRKALGO == 'HYBRID'):
    process.TTTracksEmulation = cms.Path(process.L1HybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1HybridTracksWithAssociators)
    NHELIXPAR = 4
    L1TRK_NAME  = "TTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

# HYBRID: extended tracking
elif (L1TRKALGO == 'HYBRID_DISPLACED'):
    process.TTTracksEmulation = cms.Path(process.L1ExtendedHybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1ExtendedHybridTracksWithAssociators)
    NHELIXPAR = 5
    L1TRK_NAME  = "TTTracksFromExtendedTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigisExtended"
    
# LEGACY ALGORITHM (EXPERTS ONLY): TRACKLET  
elif (L1TRKALGO == 'TRACKLET'):
    print "\n WARNING - this is not a recommended algorithm! Please use HYBRID (HYBRID_DISPLACED)!"
    print "\n To run the tracklet-only algorithm, please ensure you have commented out #define USEHYBRID in interface/Settings.h + recompiled! \n"
    process.TTTracksEmulation = cms.Path(process.L1HybridTracks)
    process.TTTracksEmulationWithTruth = cms.Path(process.L1HybridTracksWithAssociators)
    NHELIXPAR = 4
    L1TRK_NAME  = "TTTracksFromTrackletEmulation"
    L1TRK_LABEL = "Level1TTTracks"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"

# LEGACY ALGORITHM (EXPERTS ONLY): TMTT  
elif (L1TRKALGO == 'TMTT'):
    print "\n WARNING - this is not a recommended algorithm! Please use HYBRID (HYBRID_DISPLACED)! \n"
    process.load("L1Trigger.TrackFindingTMTT.TMTrackProducer_Ultimate_cff")
    L1TRK_PROC  =  process.TMTrackProducer
    L1TRK_NAME  = "TMTrackProducer"
    L1TRK_LABEL = "TML1TracksKF4ParamsComb"
    L1TRUTH_NAME = "TTTrackAssociatorFromPixelDigis"
    NHELIXPAR = 4
    L1TRK_PROC.EnableMCtruth = cms.bool(False) # Reduce CPU use by disabling internal histos.
    L1TRK_PROC.EnableHistos  = cms.bool(False)
    process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
    process.load("SimTracker.TrackTriggerAssociation.TrackTriggerAssociator_cff")
    process.TTTrackAssociatorFromPixelDigis.TTTracks = cms.VInputTag( cms.InputTag(L1TRK_NAME, L1TRK_LABEL) )
    process.TTTracksEmulation = cms.Path(process.offlineBeamSpot*L1TRK_PROC)
    process.TTTracksEmulationWithTruth = cms.Path(process.offlineBeamSpot*L1TRK_PROC*process.TrackTriggerAssociatorTracks)

else:
    print "ERROR: Unknown L1TRKALGO option"
    exit(1)



############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(NHELIXPAR), # use 4 or 5-parameter L1 tracking?
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag(L1TRK_NAME, L1TRK_LABEL),         # TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag(L1TRUTH_NAME, L1TRK_LABEL),  # MCTruth input 
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       # tracking in jets (--> requires AK4 genjet collection present!)
                                       TrackingInJets = cms.bool(False),
                                       GenJetInputTag = cms.InputTag("ak4GenJets", ""),
                                       # Was track quality performed?
                                       TrackQuality = cms.bool(True)
                                       )

process.ana = cms.Path(process.L1TrackNtuple)

# use this if you want to re-run the stub making
# process.schedule = cms.Schedule(process.TTClusterStub,process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.ana)

# use this if cluster/stub associators not available 
# process.schedule = cms.Schedule(process.TTClusterStubTruth,process.TTTracksEmulationWithTruth,process.ana)

# use this to only run tracking + track associator
process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth,process.ana)


############################################################
# write output dataset?
############################################################

if (WRITE_DATA):
  process.writeDataset = cms.OutputModule("PoolOutputModule",
      splitLevel = cms.untracked.int32(0),
      eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
      outputCommands = process.RAWSIMEventContent.outputCommands,
      fileName = cms.untracked.string('output_dataset.root'), ## ADAPT IT ##
      dataset = cms.untracked.PSet(
          filterName = cms.untracked.string(''),
          dataTier = cms.untracked.string('GEN-SIM')
      )
  )
  process.writeDataset.outputCommands.append('keep  *TTTrack*_*_*_*')
  process.writeDataset.outputCommands.append('keep  *TTStub*_*_*_*')

  process.pd = cms.EndPath(process.writeDataset)
  process.schedule.append(process.pd)

