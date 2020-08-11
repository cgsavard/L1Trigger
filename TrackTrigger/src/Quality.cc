/*
Track Quality Body file

C.Brown & C.Savard 07/2020
*/

#include "L1Trigger/TrackTrigger/interface/Quality.h"

using namespace std;

//Constructors

Quality::Quality() {
      Algorithms_.clear();
      ONNXmodels_.clear();
      ONNXInputNames_.clear();
      ONNXOutputNames_.clear();
      in_features_.clear();


    };


Quality::Quality(edm::ParameterSet Params){
    std::cout << "parameter set constructor" << endl;
    vector<string> Algorithms = Params.getParameter<vector<string>>("Quality_Algorithm");
    // Unpacks EDM parameter set itself to save unecessary processing within TrackProducers
    std::cout << "parameter set constructor 2" << endl;
    vector<string> Models=  Params.getParameter<vector<<string>>("ONNXmodel")
    vector<string> Inputs = Params.getParameter<vector<<string>>("ONNXInputName")
    vector<string> Outputs= Params.getParameter<vector<<string>>("ONNXOutputName")
    vector<string> in_features = Params.getParameter<vector<string>>("in_features")


    for (int i = 0; i<Algorithms.size(); ++i){
        std::cout<< "loop through algorithm constructors:" << Algorithms[i] << endl;
    

    if (Algorithms[i] == "Cut"){
        Set_Cut_Parameters(Algorithms[i],
                           (float)Params.getParameter<double>("maxZ0"),
                           (float)Params.getParameter<double>("maxEta"),
                           (float)Params.getParameter<double>("chi2dofMax"),
                           (float)Params.getParameter<double>("bendchi2Max"),
                           (float)Params.getParameter<double>("minPt"),
                           Params.getParameter<int>("nStubsmin")); 
    }

    else {

        Set_ONNX_Model(Algorithms[i],
                       Models[i],
                       Inputs[i],
                       Outputs[i],
                       in_features
                       ); 
    }

    }

}

 
vector<float> Quality::Feature_Transform(TTTrack < Ref_Phase2TrackerDigi_ > aTrack, std::vector<std::string> in_features) {
    std::cout<< "transform features" << endl;
    // List input features for MVA in proper order below, the features options are 
    // {"log_chi2","log_chi2rphi","log_chi2rz","log_bendchi2","nstubs","lay1_hits","lay2_hits",
    // "lay3_hits","lay4_hits","lay5_hits","lay6_hits","disk1_hits","disk2_hits","disk3_hits",
    // "disk4_hits","disk5_hits","rinv","tanl","z0","dtot","ltot","chi2","chi2rz","chi2rphi",
    // "bendchi2","pt","eta","nlaymiss_interior"}
    
    std::vector<float> transformed_features;

    // The following converts the 7 bit hitmask in the TTTrackword to an expected
    // 11 bit hitmask based on the eta of the track
    std::vector<int> hitpattern_binary = {0,0,0,0,0,0,0};
    std::vector<int> hitpattern_expanded_binary = {0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<float> eta_bins = {0.0,0.2,0.41,0.62,0.9,1.26,1.68,2.08,2.4};


    // Expected hitmap table, each row corresponds to an eta bin, each value corresponds to 
    // the expected layer in the expanded hit pattern. The expanded hit pattern should be
    // 11 bits but contains a 12th element so this hitmap table is symmetric, the 12th element 
    // is removed.
    int hitmap[8][7] = {{0, 1,  2,  3,  4,  5,  11},
                        {0, 1,  2,  3,  4,  5,  11},
                        {0, 1,  2,  3,  4,  5,  11},
                        {0, 1,  2,  3,  4,  5,  11},
                        {0, 1,  2,  3,  4,  5,  11},
                        {0, 1,  2,  6,  7,  8,  9 },
                        {0, 1,  7,  8,  9, 10,  11},
                        {0, 6,  7,  8,  9, 10,  11}};

    // iterate through bits of the hitpattern and compare to 1 filling the hitpattern binary vector
    int tmp_trk_hitpattern = aTrack.hitPattern();
    for (int i = 6; i >=0; i--){
      int k = tmp_trk_hitpattern >> i;
            if (k&1)
                hitpattern_binary[i] = 1;
        }

    // calculate number of missed interior layers from hitpattern
    int nbits = floor(log2(tmp_trk_hitpattern)) + 1;
    int lay_i = 0;
    int tmp_trk_nlaymiss_interior = 0;
    bool seq = 0;
    for (int i=0; i<nbits; i++){
      lay_i = ((1<<i)&tmp_trk_hitpattern)>>i; //0 or 1 in ith bit (right to left)
    
      if (lay_i && !seq) seq = 1; //sequence starts when first 1 found
      if (!lay_i && seq) tmp_trk_nlaymiss_interior++;
    }


    float eta = abs(aTrack.eta());
    int eta_size = static_cast<int>(eta_bins.size());
    // First iterate through eta bins

    for (int j=0; j<eta_size; j++)
        {
          if (eta >= eta_bins[j] && eta < eta_bins[j+1]) // if track in eta bin
          {
              // Iterate through hitpattern binary
              for (int k=0; k<=6; k++)
                  // Fill expanded binary entries using the expected hitmap table positions 
                  hitpattern_expanded_binary[hitmap[j][k]] = hitpattern_binary[k];   
          }
        }

    //hitpattern_expanded_binary.pop_back(); //remove final unused bit
    int tmp_trk_ltot = 0;
    //calculate number of layer hits
    for (int i=0; i<6; ++i)
    {
      tmp_trk_ltot += hitpattern_expanded_binary[i];
    }
    

    int tmp_trk_dtot = 0;
    //calculate number of disk hits
    for (int i=6; i<11; ++i)
    {
      tmp_trk_dtot += hitpattern_expanded_binary[i];
    }

    
    // While not strictly necessary to define these parameters,
    // it is included so each variable is named to avoid confusion
    float tmp_trk_big_invr   = 500*abs(aTrack.rInv());
    float tmp_trk_tanl  = abs(aTrack.tanL());
    float tmp_trk_z0   = abs( aTrack.z0() );
    float tmp_trk_pt = aTrack.momentum().perp();
    float tmp_trk_eta = aTrack.eta();
    float tmp_trk_chi2 = aTrack.chi2();
    float tmp_trk_chi2rphi = aTrack.chi2XY();
    float tmp_trk_chi2rz = aTrack.chi2Z();
    float tmp_trk_bendchi2 = aTrack.stubPtConsistency();
    float tmp_trk_log_chi2 = log(tmp_trk_chi2);
    float tmp_trk_log_chi2rphi = log(tmp_trk_chi2rphi);
    float tmp_trk_log_chi2rz = log(tmp_trk_chi2rz);
    float tmp_trk_log_bendchi2 = log(tmp_trk_bendchi2);

    // fill feature map
    std::map<std::string,float> feat_map; 
    feat_map["log_chi2"] = tmp_trk_log_chi2;
    feat_map["log_chi2rphi"] = tmp_trk_log_chi2rphi;
    feat_map["log_chi2rz"] = tmp_trk_log_chi2rz;
    feat_map["log_bendchi2"] = tmp_trk_log_bendchi2;
    feat_map["chi2"] = tmp_trk_chi2;
    feat_map["chi2rphi"] = tmp_trk_chi2rphi;
    feat_map["chi2rz"] = tmp_trk_chi2rz;
    feat_map["bendchi2"] = tmp_trk_bendchi2;
    feat_map["nstubs"] = float(tmp_trk_dtot+tmp_trk_ltot);
    feat_map["lay1_hits"] = float(hitpattern_expanded_binary[0]);
    feat_map["lay2_hits"] = float(hitpattern_expanded_binary[1]);
    feat_map["lay3_hits"] = float(hitpattern_expanded_binary[2]);
    feat_map["lay4_hits"] = float(hitpattern_expanded_binary[3]);
    feat_map["lay5_hits"] = float(hitpattern_expanded_binary[4]);
    feat_map["lay6_hits"] = float(hitpattern_expanded_binary[5]);
    feat_map["disk1_hits"] = float(hitpattern_expanded_binary[6]);
    feat_map["disk2_hits"] = float(hitpattern_expanded_binary[7]);
    feat_map["disk3_hits"] = float(hitpattern_expanded_binary[8]);
    feat_map["disk4_hits"] = float(hitpattern_expanded_binary[9]);
    feat_map["disk5_hits"] = float(hitpattern_expanded_binary[10]);
    feat_map["rinv"] = tmp_trk_big_invr;
    feat_map["tanl"] = tmp_trk_tanl;
    feat_map["z0"] = tmp_trk_z0;
    feat_map["dtot"] = float(tmp_trk_dtot);
    feat_map["ltot"] = float(tmp_trk_ltot);
    feat_map["pt"] = tmp_trk_pt;
    feat_map["eta"] = tmp_trk_eta;
    feat_map["nlaymiss_interior"] = float(tmp_trk_nlaymiss_interior);

    // fill tensor with track params
    for (std::string feat : in_features) 
      transformed_features.push_back(feat_map[feat]);

    return transformed_features;

    }
    
void Quality::Prediction(TTTrack < Ref_Phase2TrackerDigi_ > &aTrack) {
    
    int onnx_algorithm_i = 0;
    for (string Algorithm : this->Algorithms_){
        std::cout<< "loopin through predictions" << Algorithm << endl;
        
        if (Algorithm == "Cut"){
            // Get Track parameters
            float trk_pt = aTrack.momentum().perp();
            float trk_bend_chi2 = aTrack.stubPtConsistency();
            float trk_z0 = aTrack.z0();
            float trk_eta = aTrack.momentum().eta();
            float trk_chi2 = aTrack.chi2();
            const auto& stubRefs = aTrack.getStubRefs();
            int nStubs = stubRefs.size();

            float classification = 0.0; // Default classification is 0

            if (trk_pt >= this->minPt_ && 
                abs(trk_z0) < this->maxZ0_ && 
                abs(trk_eta) < this->maxEta_ && 
                trk_chi2 < this->chi2dofMax_ && 
                trk_bend_chi2 < this->bendchi2Max_ && 
                nStubs >= this->nStubsmin_) classification = 1.0;
                // Classification updated to 1 if conditions are met

            //aTrack.settrkMVA1(classification);
            aTrack.settrkMVAvector(classification);
            std::cout<< "cut predicted" << Algorithm << endl;

        }

        else {
                // Setup ONNX input and output names and arrays
                vector<string> ortinput_names;
                vector<string> ortoutput_names;
                cms::Ort::FloatArrays ortinput;
                cms::Ort::FloatArrays ortoutputs;

                ortinput_names.push_back(this->ONNXInputNames_.at(onnx_algorithm_i));
                ortoutput_names.push_back(this->ONNXOutputNames_.at(onnx_algorithm_i));

                vector<float> Transformed_features = Feature_Transform(aTrack,this->in_features_);
                cms::Ort::ONNXRuntime Runtime(this->ONNXmodels_.at(onnx_algorithm_i)); //Setup ONNX runtime


                //ONNX runtime recieves a vector of vectors of floats so push back the input
                // vector of float to create a 1,1,21 ortinput
                ortinput.push_back(Transformed_features);

                // batch_size 1 as only one set of transformed features is being processed
                int batch_size = 1;
                // Run classification on a batch of 1
                ortoutputs = Runtime.run(ortinput_names,ortinput,ortoutput_names,batch_size); 
                // access first value of nested vector
                if (Algorithm == "NN"){
                    aTrack.settrkMVAvector(ortoutputs[0][0]);
                    std::cout<< "nn predicted" << Algorithm << endl;
                    onnx_algorithm_i++;

                }

                // The ortoutput_names vector for the GBDT is left blank due to issues returning the correct
                // output, instead the GBDT will fill the ortoutputs with both the class prediciton and the class 
                // probabilities. 
                //ortoutputs[0][0] = class prediction based on a 0.5 threshold
                //ortoutputs[1][0] = negative class probability
                //ortoutputs[1][1] = positive class probability
                
                if (Algorithm == "GBDT"){
                    aTrack.settrkMVAvector(ortoutputs[1][1]);
                    std::cout<< "gbdt predicted" << Algorithm << endl;
                    onnx_algorithm_i++;
                }

        }
        }
}

void Quality::Set_Cut_Parameters(string Algorithm,float maxZ0, float maxEta, float chi2dofMax,float bendchi2Max, float minPt, int nStubmin) {

    Algorithms_.push_back(Algorithm);
    maxZ0_ = maxZ0;
    maxEta_ = maxEta; 
    chi2dofMax_ = chi2dofMax;
    bendchi2Max_ = bendchi2Max;
    minPt_ = minPt;
    nStubsmin_ = nStubmin;

}

void Quality::Set_ONNX_Model(string Algorithm,string ONNXmodel,string ONNXInputName,string ONNXOutputName, vector<string> in_features) {

    Algorithms_.push_back(Algorithm);
    ONNXmodels_.push_back(ONNXmodel);
    ONNXInputNames_.push_back(ONNXInputName);
    ONNXOutputNames_.push_back(ONNXOutputName);
    in_features_ = in_features;

}
