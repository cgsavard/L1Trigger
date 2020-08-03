#ifndef L1Trigger_TrackTrigger_interface_Quality_h
#define L1Trigger_TrackTrigger_interface_Quality_h

#include <iostream>
#include <set>
#include <vector>
#include <memory>
#include <string>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"


class Quality{
  public:
    Quality();
    Quality(std::string Algorithm,
            std::string ONNXmodel,
            std::string ONNXInputName,
            std::string ONNXOutputName,
            std::vector<std::string> in_features
            );

    Quality(std::string Algorithm,
            float maxZ0,
            float maxEta, 
            float chi2dofMax,
            float bendchi2Max,
            float minPt,
            int nStubsmin
            );

    Quality(edm::ParameterSet Params);

    ~Quality() = default;

    std::vector<float> Feature_Transform(TTTrack <Ref_Phase2TrackerDigi_> aTrack, 
					 std::vector<std::string> in_features);
    
    void Prediction(TTTrack <Ref_Phase2TrackerDigi_> &aTrack);

    void Set_Cut_Parameters(std::string Algorithm,float maxZ0, float maxEta, float chi2dofMax, float bendchi2Max, 
			    float minPt, int nStubmin);
    void Set_ONNX_Model(std::string Algorithm, std::string ONNXmodel, std::string ONNXInputName,
			std::string ONNXOutputName, std::vector<std::string> in_features);



  private:
    std::string Algorithm_ = "None";
    std::string ONNXmodel_;
    std::string ONNXInputName_;
    std::string ONNXOutputName_;
    std::vector<std::string> in_features_;
    float maxZ0_ = 15;
    float maxEta_ = 2.4; 
    float chi2dofMax_ = 40;
    float bendchi2Max_ = 2.4;
    float minPt_ = 2.0;
    int nStubsmin_ = 4;
    
  };
#endif
