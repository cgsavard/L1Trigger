import FWCore.ParameterSet.Config as cms
from L1Trigger.TrackTrigger.TrackQualityParams_cfi import *

TTTracksFromTrackletEmulation = cms.EDProducer("L1FPGATrackProducer",
                                               TTStubSource = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                               readMoreMcTruth = cms.bool(True),
                                               MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                               MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                               TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                               TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                               BeamSpotSource = cms.InputTag("offlineBeamSpot"),
                                               asciiFileName = cms.untracked.string(""),
                                               # (if running on CRAB use "../../fitpattern.txt" etc instead)
                                               Extended=cms.bool(False),
                                               Hnpar=cms.uint32(4),
                                               fitPatternFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/fitpattern.txt'),
                                               memoryModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/memorymodules_hourglass.dat'),
                                               processingModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/processingmodules_hourglass.dat'),
                                               wiresFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/wires_hourglass.dat'),
                                               DTCLinkFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/calcNumDTCLinks.txt'),
                                               DTCLinkLayerDiskFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/dtclinklayerdisk.dat'),
                                               moduleCablingFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/modules_T5v3_27SP_nonant_tracklet.dat'),
                                               # Quality Flag and Quality params
                                               Quality =cms.bool(True),
                                               TrackQualityPSet = cms.PSet(TrackQualityParams)

    )

TTTracksFromExtendedTrackletEmulation = TTTracksFromTrackletEmulation.clone(
                                               Extended=cms.bool(True),
                                               Hnpar=cms.uint32(5),
                                               memoryModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/memorymodules_hourglassExtended.dat'),
                                               processingModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/processingmodules_hourglassExtended.dat'),
                                               wiresFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/wires_hourglassExtended.dat'),
                                               # specifying where the TrackletEngineDisplaced(TED)/TripletEngine(TRE) tables are located
                                               tableTEDFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/table_TED/table_TED_D1PHIA1_D2PHIA1.txt'),
                                               tableTREFile = cms.FileInPath('L1Trigger/TrackFindingTracklet/data/table_TRE/table_TRE_D1AD2A_1.txt'),
                                               Quality =cms.bool(False),
                                               TrackQualityPSet = cms.PSet(TrackQualityParams)
    )


