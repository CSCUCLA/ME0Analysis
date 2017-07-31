#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "vector"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/Types.h"

using namespace std;
using ASTypes::CylLorentzVectorF;
//void go() {
//
//  TString trees[] = {"bkg_tree", ""};
//
//
//
//  auto proc = [](TString pre){
//    TFile * fI =   new TFile(pre +".root","Read");
//    TTree *tI =0;
//    fI->GetObject("Events",tI);
//    tI->SetBranchStatus("*",1);
//
//    TFile * fO = new TFile(pre +"_lepFJSkim.root","recreate");
//    fO->cd();
//    TTree * tO = tI->CopyTree("lepton_pt@.size()==1 && fj_pt@.size() > 0");
//    tO->Write();
//    fO->Close();
//  };
//
//  for(unsigned int iT = 0; trees[iT][0]; ++iT)proc(trees[iT]);
//}



void go()  {
	TFile* file = new TFile("tree_tau3mu_gen.root","Read");
	TTreeReader reader("demo/data", file);
	TTreeReaderArray<double> mu_pt(reader, "mu_pt");
	TTreeReaderArray<double> mu_eta(reader, "mu_eta");
	TTreeReaderArray<double> mu_phi(reader, "mu_phi");
	TTreeReaderArray<double> mu_rec(reader, "mu_rec");

	HistGetter hists;

	while (reader.Next()) {
		std::vector<CylLorentzVectorF> muons;
		for(unsigned int iM = 0; iM < 3; ++iM){
			muons.push_back(CylLorentzVectorF(mu_pt[iM],mu_eta[iM],mu_phi[iM],.105));
		}

		int nCentral = 0;
		int nCSC = 0;
		int nME0 = 0;

		for(const auto& muon : muons ){
			if(std::fabs(muon.eta()) < 2.4 && muon.pt() > 3) nCentral++;
			else if(std::fabs(muon.eta()) > 1.6 && std::fabs(muon.eta()) < 2.4 && muon.P() > 5) nCSC++;
			else if(std::fabs(muon.eta()) > 2.4 && std::fabs(muon.eta()) < 2.8 && muon.P() > 8) nME0++;
		}

		hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(0);
		if(nCentral >= 2) hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(1);
		if(nCentral >= 3) hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(2);
		if(nCentral >= 2 && nCSC > 0 ) hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(3);
		if(nCentral >= 2 && nME0 > 0 ) hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(4);
		if(nCentral >= 2 && (nME0 +nCSC) > 0 ) hists.getOrMake1D("categories",";Incl/2 central/3 central/2 central + 1 CSC/2 central + 1 ME0/2 central + 1 ME0 or CSC",6,-0.5,5.5)->Fill(5);

		int nInclCentral = 0;
		int nPT2Central = 0;
		int nInclME0 = 0;
		for(const auto& muon : muons ){

			if(std::fabs(muon.eta()) < 2.4) nInclCentral++;
			if(std::fabs(muon.eta()) < 2.4 && muon.pt() > 2) nPT2Central++;
			else if(std::fabs(muon.eta()) > 2.4 && std::fabs(muon.eta()) < 2.8 && muon.pt() > 1) nInclME0++;
		}

		hists.getOrMake2D("TwoDCat",";N ME0 muons; nCentral muons",4,-0.5,3.5,4,-0.5,3.5)->Fill(nInclME0,nInclCentral);
		hists.getOrMake2D("TwoD2PTCat",";N ME0 muons; nCentral muons",4,-0.5,3.5,4,-0.5,3.5)->Fill(nInclME0,nPT2Central);



	}

	hists.write("outplots.root");

}

#endif

void LookAtTauThreemuTree(){
  go();
}
