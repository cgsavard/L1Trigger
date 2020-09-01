/*
Track Quality Header file

C.Brown 28/07/20
*/

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
    //Default Constructor
    Quality() {};

    Quality(edm::ParameterSet Params);

    //Default Destructor
    ~Quality() = default;
    
    // Controls the conversion between TTTrack features and ML model training features
    std::vector<float> Feature_Transform(TTTrack <Ref_Phase2TrackerDigi_> aTrack, 
					                               std::vector<std::string> in_features);
    
    // Passed by reference a track without MVA filled, fills the track's MVA field
    void Prediction(TTTrack <Ref_Phase2TrackerDigi_> &aTrack);

    void Set_ONNX_Model(std::string Algorithm, std::string ONNXmodel, std::string ONNXInputName,
			std::vector<std::string> in_features);


  private:
    // Private Memember Data
    std::string Algorithm_ = "None";
    std::string ONNXmodel_ = "None";
    std::string ONNXInputName_;
    std::vector<std::string> in_features_;
    cms::Ort::ONNXRuntime Runtime(model_path=ONNXmodel_);
    
  };
#endif
