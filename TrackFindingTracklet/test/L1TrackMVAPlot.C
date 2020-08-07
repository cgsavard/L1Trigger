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

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

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

  SetPlotStyle();
  TCanvas c;

  // ----------------------------------------------------------------------------------------------------------------
  // define input options

  // these are the LOOSE cuts, baseline scenario for efficiency and rate plots ==> configure as appropriate
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

  TH1F* h_trk_MVA1_real = new TH1F("trk_MVA1_real",";MVA1; L1 tracks", 50, 0, 1);
  h_trk_MVA1_real->SetLineColor(1);
  TH1F* h_trk_MVA1_fake = new TH1F("trk_MVA1_fake",";MVA1; L1 tracks", 50, 0, 1);
  h_trk_MVA1_fake->SetLineColor(2);

  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------

  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;

  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  vector<float> MVA1s;
  vector<float> fakes;
  vector<float> etas;
  vector<float> pts;
  for (int i = 0; i < nevt; i++) {
    tree->GetEntry(i, 0);

    for (int it = 0; it < (int)trk_pt->size(); it++) {
      // ----------------------------------------------------------------------------------------------------------------
      // track properties

      float MVA1 = trk_MVA1->at(it);
      float fake = trk_fake->at(it);
      float eta = trk_eta->at(it);
      float pt = trk_pt->at(it);

      MVA1s.push_back(MVA1);
      fakes.push_back(fake);
      etas.push_back(eta);
      pts.push_back(pt);

      h_trk_MVA1->Fill(MVA1);
      if (fake==1.) h_trk_MVA1_real->Fill(MVA1);
      else if (fake==0.) h_trk_MVA1_fake->Fill(MVA1);
    }
  }

  // -------------------------------------------------------------------------------------------
  // create ROC curve
  // -------------------------------------------------------------------------------------------

  vector<float> TPR;
  vector<float> FPR;
  int n = 100; //num of entries on ROC curve
  for (int i=0; i<n; i++){
    float dt = (float)i/(n-1); //make sure it starts at (0,0) and ends at (1,1)
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

  // calculate AUC
  float AUC = 0.0;
  for (int i=0; i<n-1; i++){
    AUC += (TPR[i]+TPR[i+1])/2*(FPR[i]-FPR[i+1]);
  }

  TGraph* ROC = new TGraph(n, FPR.data(), TPR.data());
  ROC->SetName("ROC_MVA1");
  ROC->SetTitle(("ROC curve (AUC = "+to_string(AUC)+"); FPR; TPR").c_str());  

  // -------------------------------------------------------------------------------------------
  // create TPR vs. eta and FPR vs. eta
  // -------------------------------------------------------------------------------------------

  vector<float> TPR_eta, TPR_eta_err;
  vector<float> FPR_eta, FPR_eta_err;
  vector<float> eta_range_TPR, eta_range_TPR_err;
  vector<float> eta_range_FPR, eta_range_FPR_err;
  n = 20;
  float eta_low = -2.4;
  float eta_high = 2.4;
  float eta_temp = eta_low;
  float eta_step = (eta_high-eta_low)/n;
  float dt = .5;
  for (int ct=0; ct<n; ct++){
    int TP = 0;
    int FP = 0;
    int P = 0;
    int N = 0;
    for (int k=0; k<etas.size(); k++){
      if (etas.at(k)>eta_temp && etas.at(k)<=(eta_temp+eta_step)){
	if (fakes.at(k)){
	  P++;
	  if (MVA1s.at(k)>dt) TP++;
	}else{
	  N++;
	  if (MVA1s.at(k)>dt) FP++;
	}
      }
    }
    if (P>0){
      TPR_eta.push_back((float)TP/P);
      TPR_eta_err.push_back(sqrt(TP*(P-TP)/pow(P,3)));
      eta_range_TPR.push_back(eta_temp+eta_step/2); //halfway in bin
      eta_range_TPR_err.push_back(eta_step/2);
    }if (N>0){
      FPR_eta.push_back((float)FP/N);
      FPR_eta_err.push_back(sqrt(FP*(N-FP)/pow(N,3)));
      eta_range_FPR.push_back(eta_temp+eta_step/2); //halfway in bin
      eta_range_FPR_err.push_back(eta_step/2);
    }
    eta_temp += eta_step;
  }

  TGraphErrors* TPR_vs_eta = new TGraphErrors(TPR_eta.size(), eta_range_TPR.data(), TPR_eta.data(),
  					      eta_range_TPR_err.data(), TPR_eta_err.data());
  TPR_vs_eta->SetName("TPR_vs_eta_MVA1");
  TPR_vs_eta->SetTitle("TPR vs. #eta; #eta; TPR");

  TGraphErrors* FPR_vs_eta = new TGraphErrors(FPR_eta.size(), eta_range_FPR.data(), FPR_eta.data(),
					      eta_range_FPR_err.data(), FPR_eta_err.data());
  FPR_vs_eta->SetName("FPR_vs_eta_MVA1");
  FPR_vs_eta->SetTitle("FPR vs. #eta; #eta; FPR");

  // -------------------------------------------------------------------------------------------
  // create TPR vs. pt and FPR vs. pt
  // -------------------------------------------------------------------------------------------

  vector<float> TPR_pt, TPR_pt_err;
  vector<float> FPR_pt, FPR_pt_err;
  vector<float> pt_range_TPR, pt_range_TPR_err;
  vector<float> pt_range_FPR, pt_range_FPR_err;
  n = 10;
  float logpt_low = log10(2); //set low pt in log
  float logpt_high = log10(100); //set high pt in log
  float logpt_temp = logpt_low;
  float logpt_step = (logpt_high-logpt_low)/n;
  dt = .5;
  for (int ct=0; ct<n; ct++){
    int TP = 0;
    int FP = 0;
    int P = 0;
    int N = 0;
    for (int k=0; k<pts.size(); k++){
      if (pts.at(k)>pow(10,logpt_temp) && pts.at(k)<=(pow(10,logpt_temp+logpt_step))){
	if (fakes.at(k)){
	  P++;
	  if (MVA1s.at(k)>dt) TP++;
	}else{
	  N++;
	  if (MVA1s.at(k)>dt) FP++;
	}
      }
      }
    if (P>0){
      TPR_pt.push_back((float)TP/P);
      TPR_pt_err.push_back(sqrt(TP*(P-TP)/pow(P,3)));
      pt_range_TPR.push_back((pow(10,logpt_temp)+pow(10,logpt_temp+logpt_step))/2); //halfway in bin
      pt_range_TPR_err.push_back((pow(10,logpt_temp+logpt_step)-pow(10,logpt_temp))/2);
    }if (N>0){
      FPR_pt.push_back((float)FP/N);
      FPR_pt_err.push_back(sqrt(FP*(N-FP)/pow(N,3)));
      pt_range_FPR.push_back((pow(10,logpt_temp)+pow(10,logpt_temp+logpt_step))/2); //halfway in bin
      pt_range_FPR_err.push_back((pow(10,logpt_temp+logpt_step)-pow(10,logpt_temp))/2);
    }
    logpt_temp += logpt_step;
  }

  TGraphErrors* TPR_vs_pt = new TGraphErrors(TPR_pt.size(), pt_range_TPR.data(), TPR_pt.data(),
					     pt_range_TPR_err.data(), TPR_pt_err.data());
  TPR_vs_pt->SetName("TPR_vs_pt_MVA1");
  TPR_vs_pt->SetTitle("TPR vs. p_{T}; p_{T}; TPR");

  TGraphErrors* FPR_vs_pt = new TGraphErrors(FPR_pt.size(), pt_range_FPR.data(), FPR_pt.data(),
					     pt_range_FPR_err.data(), FPR_pt_err.data());
  FPR_vs_pt->SetName("FPR_vs_pt_MVA1");
  FPR_vs_pt->SetTitle("FPR vs. p_{T}; p_{T}; FPR");

  // -------------------------------------------------------------------------------------------
  // output file for histograms and graphs
  // -------------------------------------------------------------------------------------------

  TFile* fout;
  fout = new TFile(type_dir + "MVAoutput_" + type + treeName + ".root", "recreate");

  // -------------------------------------------------------------------------------------------
  // draw and save plots
  // -------------------------------------------------------------------------------------------

  h_trk_MVA1->Draw();
  h_trk_MVA1->Write();

  h_trk_MVA1_real->Draw();
  h_trk_MVA1_fake->Draw("same");
  TLegend* leg = new TLegend();
  leg->AddEntry(h_trk_MVA1_real,"real","l");
  leg->AddEntry(h_trk_MVA1_fake,"fake","l");
  leg->Draw("same");
  c.Write("trk_MVA1_rf");

  ROC->Draw();
  ROC->Write();
  c.SaveAs("ROC.pdf");

  TPR_vs_eta->Draw("ap");
  TPR_vs_eta->Write();
  c.SaveAs("TPR_vs_eta.pdf");

  FPR_vs_eta->Draw("ap");
  FPR_vs_eta->Write();
  c.SaveAs("FPR_vs_eta.pdf");

  TPR_vs_pt->Draw("ap");
  TPR_vs_pt->Write();
  c.SaveAs("TPR_vs_eta.pdf");

  FPR_vs_pt->Draw("ap");
  FPR_vs_pt->Write();
  c.SaveAs("FPR_vs_eta.pdf");

  fout->Close();
}

void SetPlotStyle() {
  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
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
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

