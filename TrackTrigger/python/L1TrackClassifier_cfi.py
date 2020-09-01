import FWCore.ParameterSet.Config as cms

TrackQualityParams = cms.PSet(Quality_Algorithm = cms.string("NN"), #None, Cut, NN, GBDT
                              ONNXmodel = cms.string("../../TrackTrigger/ML_data/FakeIDNN/NN_model.onnx"),
                              ONNXInputName = cms.string("input_1"),
                              #Vector of strings of training features, in the order that the model was trained with
                              in_features = cms.vstring(["log_chi2","log_bendchi2","log_chi2rphi","log_chi2rz",
                                                                            "nstubs","lay1_hits","lay2_hits","lay3_hits","lay4_hits",
                                                                            "lay5_hits","lay6_hits","disk1_hits","disk2_hits",
                                                                            "disk3_hits","disk4_hits","disk5_hits","rinv","tanl",
                                                                            "z0","dtot","ltot"]),

