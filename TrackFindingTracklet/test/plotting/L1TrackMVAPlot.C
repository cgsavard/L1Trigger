// ----------------------------------------------------------------------------------------------------------------
// Basic example ROOT script for making tracking performance plots using the ntuples produced by L1TrackNtupleMaker.cc
//
// e.g. in ROOT do: .L L1TrackMVAPlot.C++, L1TrackMVAPlot("TTbar_PU200_D49")
//
// By Claire Savard, July 2020
// Based off of L1TrackNtuplePlot.C
// ----------------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TGraph.h"
#include <TError.h>
#include "tdrstyle.C"
#include "CMS_lumi.h"



#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void SetPlotStyle();

// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------

void L1TrackMVAPlot(TString type,
                    TString type_dir = "",
		    TString treeName = ""){
  // type:              this is the name of the input file you want to process (minus ".root" extension)
  // type_dir:          this is the directory containing the input file you want to process. Note that this must end with a "/", as in "EventSets/"

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  
  //gROOT->LoadMacro("tdrstyle.C");
  
  SetPlotStyle();
  // ----------------------------------------------------------------------------------------------------------------
  // define input options

  // these are the LOOSE cuts, baseline scenario for efficiency and rate plots ==> configure as appropriate
  TCanvas c;

  //TString DIR = "./";

  int L1Tk_minNstub = 4;
  float L1Tk_maxChi2 = 999999;
  float L1Tk_maxChi2dof = 999999.;

  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TrackNtuple" + treeName + "/eventTree");
  tree->Add(type_dir + type + ".root");

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..."
         << endl;  //cout's kept in this file as it is an example standalone plotting script, not running in central CMSSW
    return;
  }

  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches
  // all L1 tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_chi2;
  vector<float>* trk_chi2rphi;
  vector<float>* trk_chi2rz;
  vector<int>* trk_nstub;
  vector<int>* trk_lhits;
  vector<int>* trk_dhits;
  vector<int>* trk_seed;
  vector<int>* trk_hitpattern;
  vector<unsigned int>* trk_phiSector;
  vector<int>* trk_fake;
  vector<int>* trk_genuine;
  vector<int>* trk_loose;
  vector<float>* trk_MVA1;

  TBranch* b_trk_pt;
  TBranch* b_trk_eta;
  TBranch* b_trk_phi;
  TBranch* b_trk_chi2;
  TBranch* b_trk_chi2rphi;
  TBranch* b_trk_chi2rz;
  TBranch* b_trk_nstub;
  TBranch* b_trk_lhits;
  TBranch* b_trk_dhits;
  TBranch* b_trk_phiSector;
  TBranch* b_trk_seed;
  TBranch* b_trk_hitpattern;
  TBranch* b_trk_fake;
  TBranch* b_trk_genuine;
  TBranch* b_trk_loose;
  TBranch* b_trk_MVA1;

  trk_pt = 0;
  trk_eta = 0;
  trk_phi = 0;
  trk_chi2 = 0;
  trk_chi2rphi = 0;
  trk_chi2rz = 0;
  trk_nstub = 0;
  trk_lhits = 0;
  trk_dhits = 0;
  trk_phiSector = 0;
  trk_seed = 0;
  trk_hitpattern = 0;
  trk_fake = 0;
  trk_genuine = 0;
  trk_loose = 0;
  trk_MVA1 = 0;

  tree->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  tree->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
  tree->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
  tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
  tree->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
  tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
  tree->SetBranchAddress("trk_lhits", &trk_lhits, &b_trk_lhits);
  tree->SetBranchAddress("trk_dhits", &trk_dhits, &b_trk_dhits);
  tree->SetBranchAddress("trk_phiSector", &trk_phiSector, &b_trk_phiSector);
  tree->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
  tree->SetBranchAddress("trk_hitpattern", &trk_hitpattern, &b_trk_hitpattern);
  tree->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
  tree->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
  tree->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
  tree->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);

  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h_trk_MVA1 = new TH1F("trk_MVA1", "; MVA1; L1 tracks", 50, 0, 1);

  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------

  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;

  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  vector<float> MVA1s;
  vector<float> fakes;
  for (int i = 0; i < nevt; i++) {
    tree->GetEntry(i, 0);

    for (int it = 0; it < (int)trk_pt->size(); it++) {
      // ----------------------------------------------------------------------------------------------------------------
      // track properties

      float MVA1 = trk_MVA1->at(it);
      float fake = trk_fake->at(it);

      MVA1s.push_back(MVA1);
      fakes.push_back(fake);

      h_trk_MVA1->Fill(MVA1);
    }
  }

  // -------------------------------------------------------------------------------------------
  // create ROC curve
  // -------------------------------------------------------------------------------------------

  vector<float> TPR;
  vector<float> FPR;
  int n = 30; //num of entries on ROC curve
  for (int i=0; i<n; i++){
    float dt = (float)i/n;
    int TP = 0;
    int FP = 0;
    int P = 0;
    int N = 0;
    for (int k=0; k<MVA1s.size(); k++){
      if (fakes.at(k)){
	P++;
	if (MVA1s.at(k)>dt) TP++;
      }else{
	N++;
	if (MVA1s.at(k)>dt) FP++;
      }
    }
    TPR.push_back((float)TP/P);
    FPR.push_back((float)FP/N);
  }

  TGraph* ROC = new TGraph(n, &FPR[0], &TPR[0]);
  ROC->SetName("ROC_MVA1");
  ROC->SetTitle("ROC curve; FPR; TPR");

  // -------------------------------------------------------------------------------------------
  // output file for histograms
  // -------------------------------------------------------------------------------------------

  TFile* fout;
  fout = new TFile(type_dir + "MVAoutput_" + type + treeName + ".root", "recreate");

  // -------------------------------------------------------------------------------------------
  // draw and save plots
  // -------------------------------------------------------------------------------------------



  h_trk_MVA1->SetLineColor(kOrange+7);
  h_trk_MVA1->SetFillColor(kOrange-2);
  h_trk_MVA1->Draw();
  c.SaveAs("MVAhist.pdf");


  h_trk_MVA1->SetMarkerColor(kOrange+7);
  ROC->SetMarkerStyle(20);
  ROC->SetMarkerSize(0.5);
  ROC->Draw("alp");
  c.SaveAs("ROC.pdf");

  h_trk_MVA1->Write();
  ROC->Write();
  

  fout->Close();
}


void SetPlotStyle(){
  setTDRStyle();

  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42, "x");
  gStyle->SetTitleFont(42, "x");
  gStyle->SetLabelFont(42, "y");
  gStyle->SetTitleFont(42, "y");
  gStyle->SetLabelFont(42, "z");
  gStyle->SetTitleFont(42, "z");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetLabelSize(0.05, "z");
  gStyle->SetTitleSize(0.05, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2, "[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);



}
