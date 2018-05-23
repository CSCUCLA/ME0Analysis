#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

void go() {

  TString trees[] = {"ttbar","wjets","singlet",""};



  auto proc = [](TString pre){
    TFile * fI = new TFile(pre +"_tree.root","Read");
    TTree *tI =0;
    fI->GetObject("Events",tI);
    tI->SetBranchStatus("*",1);
    // turn off other jets
	tI->SetBranchStatus("fj1p5_pt"        ,0);
	tI->SetBranchStatus("fj1p5_eta"       ,0);
	tI->SetBranchStatus("fj1p5_phi"       ,0);
	tI->SetBranchStatus("fj1p5_mass"      ,0);
	tI->SetBranchStatus("fj1p5_sd_massFFJ",0);
	tI->SetBranchStatus("fj1p5_csv"       ,0);
	tI->SetBranchStatus("fj1p5_bbcsv"     ,0);
	tI->SetBranchStatus("fj1p5_tau1"      ,0);
	tI->SetBranchStatus("fj1p5_tau2"      ,0);
	tI->SetBranchStatus("fj1p5_tau3"      ,0);
	tI->SetBranchStatus("fj1p5_sd_pt"     ,0);
	tI->SetBranchStatus("fj1p5_sd_eta"    ,0);
	tI->SetBranchStatus("fj1p5_sd_phi"    ,0);
	tI->SetBranchStatus("fj1p5_sd_mass"   ,0);
	tI->SetBranchStatus("fj1p5_sj1_pt"    ,0);
	tI->SetBranchStatus("fj1p5_sj1_eta"   ,0);
	tI->SetBranchStatus("fj1p5_sj1_phi"   ,0);
	tI->SetBranchStatus("fj1p5_sj1_mass"  ,0);
	tI->SetBranchStatus("fj1p5_sj1_csv"   ,0);
	tI->SetBranchStatus("fj1p5_sj2_pt"    ,0);
	tI->SetBranchStatus("fj1p5_sj2_eta"   ,0);
	tI->SetBranchStatus("fj1p5_sj2_phi"   ,0);
	tI->SetBranchStatus("fj1p5_sj2_mass"  ,0);
	tI->SetBranchStatus("fj1p5_sj2_csv"   ,0);
	tI->SetBranchStatus("jet_pt"          ,0);
	tI->SetBranchStatus("jet_eta"         ,0);
	tI->SetBranchStatus("jet_phi"         ,0);
	tI->SetBranchStatus("jet_mass"        ,0);
	tI->SetBranchStatus("jet_csv"         ,0);

    TFile * fO = new TFile(pre +"_lep0p8FJSkim.root","recreate");
    fO->cd();
//    TTree * tO = tI->CopyTree("lepton_pt@.size()==1 && (fj_pt@.size() > 1 || fj1p5_pt@.size() > 1)");
    TTree * tO = tI->CopyTree("lepton_pt@.size()==1");




    tO->Write();
    fO->Close();
  };

  for(unsigned int iT = 0; trees[iT][0]; ++iT)proc(trees[iT]);
}

#endif

void MakeCutSkims(){
  go();
}
