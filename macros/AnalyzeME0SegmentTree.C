
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "HistoPlotting/include/Plotter.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TGraph.h"
#include "TFitter.h"
#include "TProfile.h"
#include "TRandom3.h"

using namespace std;
using namespace PhysicsUtilities;
using ASTypes::CylLorentzVectorF;
typedef std::vector<std::tuple<double,double,double>>  MeasureList;//pt DPhi eta



class Analyzer : public BaseTupleAnalyzer{
public:
	TString glbPrefix = "";

	Analyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){
		setBranchAddress("simMuon_pt"               ,&simMuon_pt               );
		setBranchAddress("simMuon_eta"              ,&simMuon_eta              );
		setBranchAddress("simMuon_phi"              ,&simMuon_phi              );
		setBranchAddress("simMuon_q"                ,&simMuon_q                );
		setBranchAddress("simMuon_gen_eta"          ,&simMuon_gen_eta          );
		setBranchAddress("simMuon_gen_phi"          ,&simMuon_gen_phi          );
		setBranchAddress("simMuon_gen_dphi"         ,&simMuon_gen_dphi         );
		setBranchAddress("simMuon_gen_deta"         ,&simMuon_gen_deta         );
		setBranchAddress("simMuon_gen_x"            ,&simMuon_gen_x            );
		setBranchAddress("simMuon_gen_y"            ,&simMuon_gen_y            );
		setBranchAddress("simMuon_gen_dx"           ,&simMuon_gen_dx           );
		setBranchAddress("simMuon_gen_dy"           ,&simMuon_gen_dy           );
		setBranchAddress("simMuon_segmentIDX"       ,&simMuon_segmentIDX       );
		setBranchAddress("simMuon_segment_quality"  ,&simMuon_segment_quality  );
		setBranchAddress("simMuon_segment_nGoodHits",&simMuon_segment_nGoodHits);
		setBranchAddress("simMuon_segment_nBadHits" ,&simMuon_segment_nBadHits );
		setBranchAddress("segment_eta"              ,&segment_eta              );
		setBranchAddress("segment_phi"              ,&segment_phi              );
		setBranchAddress("segment_dphi"             ,&segment_dphi             );
		setBranchAddress("segment_deta"             ,&segment_deta             );
		setBranchAddress("segment_x"                ,&segment_x                );
		setBranchAddress("segment_y"                ,&segment_y                );
		setBranchAddress("segment_dx"               ,&segment_dx               );
		setBranchAddress("segment_dy"               ,&segment_dy               );
	}


	std::vector<float>* simMuon_pt               = new vector<float>;
	std::vector<float>* simMuon_eta              = new vector<float>;
	std::vector<float>* simMuon_phi              = new vector<float>;
	std::vector<int>  * simMuon_q                = new vector<int>  ;

	std::vector<float>* simMuon_gen_eta          = new vector<float>;
	std::vector<float>* simMuon_gen_phi          = new vector<float>;
	std::vector<float>* simMuon_gen_dphi         = new vector<float>;
	std::vector<float>* simMuon_gen_deta         = new vector<float>;
	std::vector<float>* simMuon_gen_x            = new vector<float>;
	std::vector<float>* simMuon_gen_y            = new vector<float>;
	std::vector<float>* simMuon_gen_dx           = new vector<float>;
	std::vector<float>* simMuon_gen_dy           = new vector<float>;

	std::vector<int>  * simMuon_segmentIDX       = new vector<int>  ;
	std::vector<int>  * simMuon_segment_quality  = new vector<int>  ;
	std::vector<int>  * simMuon_segment_nGoodHits= new vector<int>  ;
	std::vector<int>  * simMuon_segment_nBadHits = new vector<int>  ;
	std::vector<float>* segment_eta              = new vector<float>;
	std::vector<float>* segment_phi              = new vector<float>;
	std::vector<float>* segment_dphi             = new vector<float>;
	std::vector<float>* segment_deta             = new vector<float>;
	std::vector<float>* segment_x                = new vector<float>;
	std::vector<float>* segment_y                = new vector<float>;
	std::vector<float>* segment_dx               = new vector<float>;
	std::vector<float>* segment_dy               = new vector<float>;

	void write(TString fileName){ plotter.write(fileName);}
	HistGetter plotter;

	virtual void runAEvent() {
	}
};

class FindParams : public Analyzer {
public:

	MeasureList genEtaList;
	MeasureList genThetaList;
	MeasureList genSinThetaList;
	MeasureList recoEtaList;
	MeasureList recoThetaList;
	MeasureList recoSinThetaList;
	FindParams(std::string fileName, std::string treeName) : Analyzer(fileName,treeName) {}

	virtual void runAEvent() {
		for(unsigned int iS = 0; iS < simMuon_pt->size(); ++iS){
			CylLorentzVectorF muon(simMuon_pt->at(iS),simMuon_eta->at(iS),simMuon_phi->at(iS),0.105);
			if(muon.pt() < 2) continue;
			double genTheta = muon.theta() < TMath::PiOver2() ?  muon.theta() : TMath::Pi() - muon.theta() ;
			genEtaList     .emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),std::fabs(muon.eta()) );
			genThetaList   .emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),genTheta );
			genSinThetaList.emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),std::sin(genTheta) );
			if(simMuon_segment_quality->at(iS) == 3) continue;
			const int idx = simMuon_segmentIDX->at(iS);
			const double theta = PhysicsUtilities::etaToTheta(std::fabs(segment_eta->at(idx)));
			recoEtaList     .emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),std::fabs(segment_eta->at(idx)) );
			recoThetaList   .emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),theta );
			recoSinThetaList.emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),std::sin(theta) );

			auto solPT = [](double dPhi, double eta, double a, double b) -> double { return (a + b*eta)/dPhi; };

			double genSolPT = solPT(std::fabs(simMuon_gen_dphi->at(iS)),genTheta,-0.004497,0.1373 );
			double recoSolPT = solPT(std::fabs(segment_dphi->at(idx)),theta,-0.004497,0.1373   );

			double genSolPT2 = solPT(std::fabs(simMuon_gen_dphi->at(iS)),std::fabs(muon.eta()),0.07292,-0.02116    );
			double recoSolPT2 =  solPT(std::fabs(segment_dphi->at(idx)),std::fabs(segment_eta->at(idx)),0.07292,-0.02116  );

			plotter.getOrMake1D("gen_eta_PToPT",";reco pt / gen pt"       ,100,0,10)->Fill(  genSolPT2/muon.pt()  );
			plotter.getOrMake1D("gen_theta_PToPT",";reco pt / gen pt"     ,100,0,10)->Fill(  genSolPT/muon.pt()  );
			plotter.getOrMake1D("reco_eta_PToPT",";reco pt / gen pt"      ,100,0,10)->Fill(  recoSolPT2/muon.pt()  );
			plotter.getOrMake1D("reco_theta_PToPT",";reco pt / gen pt",100,0,10)->Fill(  recoSolPT/muon.pt()  );



			TString selString = "";
			if(muon.pt() < 5) selString = "pt1to5_";
			else if(muon.pt() < 15) selString = "pt5to15_";
			else selString = "pt15to30_";
//			if(std::fabs(muon.eta()) < 2.2) selString += "eta2to2p2";
//			else if(std::fabs(muon.eta()) < 2.4) selString += "eta2p2to2p4";
//			else if(std::fabs(muon.eta()) < 2.6) selString += "eta2p2to2p6";
//			else selString += "eta2p6to2p8";

			if(std::fabs(muon.eta()) < 2.4) selString += "eta2to2p2";
			else selString += "eta2p6to2p8";


 			plotter.getOrMake1D(TString::Format("gen_%s_PToPT",selString.Data()),";segment p_{T} / gen p_{T}",100,0,4)->Fill(  genSolPT/muon.pt()  );
 			plotter.getOrMake1D(TString::Format("reco_%s_PToPT",selString.Data()),";segment p_{T} / gen p_{T}",100,0,4)->Fill(  recoSolPT/muon.pt()  );




		}
	}

	void solve() {

	    TFile * file = TFile::Open("profiles.root","RECREATE");
		auto makeSolution = [&](const MeasureList& list, TString name, int nBins, double min, double max){
			TProfile * n = new TProfile(TString::Format("%s_N",name.Data()),";|eta|",nBins,min,max);
			for(const auto& val : list ){
				if(std::get<1>(val) > .1) continue;
				if(std::get<1>(val) < 0.0008) continue;
				if(std::get<0>(val) > 10) continue;
				n->Fill(std::get<2>(val) ,  1/(std::get<0>(val)*std::get<1>(val)));
			}

			file->cd();
			n->Write();
		};
		makeSolution(genEtaList,"genEtaList",32,2,2.8);
		makeSolution(genThetaList,"genThetaList",32,.121,.269);
		makeSolution(genSinThetaList,"genSinThetaList",32,.121,.269);
		makeSolution(recoEtaList,"recoEtaList",32,2,2.8);
		makeSolution(recoThetaList,"recoThetaList",32,.121,.269);
		makeSolution(recoSinThetaList,"recoSinThetaList",32,.121,.269);
		file->Close();

	}

};

class TestParams : public Analyzer {
public:
	TRandom3 * rand = new TRandom3(12345);
	double getMeasuredPT(double dPhi, double theta, double a = -0.004497, double b = 0.1373 )const { return (a + b*theta)/dPhi;}
	CylLorentzVectorF getMuon(unsigned int idx) const { return CylLorentzVectorF(simMuon_pt->at(idx),simMuon_eta->at(idx),simMuon_phi->at(idx),0.105);}
	CylLorentzVectorF getMeasuredGen(unsigned int idx) const {
		auto muon = getMuon(idx);
		double genTheta = muon.theta() < TMath::PiOver2() ?  muon.theta() : TMath::Pi() - muon.theta() ;
		return CylLorentzVectorF(getMeasuredPT(std::fabs(simMuon_gen_dphi->at(idx)), genTheta  ),simMuon_eta->at(idx),simMuon_phi->at(idx),0.105);
	}
	CylLorentzVectorF getMeasuredReco(unsigned int idx) const {
		const double theta = PhysicsUtilities::etaToTheta(std::fabs(segment_eta->at(idx)))   ;
		return CylLorentzVectorF(getMeasuredPT(std::fabs(segment_dphi->at(idx)), theta  ),segment_eta->at(idx),segment_phi->at(idx),0.105 );
	}
	CylLorentzVectorF getRandomTrack(const double& pt) { return CylLorentzVectorF(pt,rand->Uniform(-2.4,2.4),rand->Uniform(-TMath::Pi(),TMath::Pi()),.105);}

	TestParams(std::string fileName, std::string treeName) : Analyzer(fileName,treeName) {}

	void fakeStudies() {


		int nGMuons = 0;
		for(const auto& qual : *simMuon_segment_quality ) if(qual < 3) nGMuons++;
		if(nGMuons != 2) return;
			plotter.getOrMake1D(TString::Format("%s_fs_nEvents",glbPrefix.Data()),";nEvents",1,0,2)->Fill(  1  );

			auto pt5SingleMu  = getRandomTrack(5);
			auto pt20SingleMu  = getRandomTrack(20);


		for(unsigned int idx = 0; idx < segment_eta->size(); ++idx){
			bool good = false; for(const auto& trueIDX : *simMuon_segmentIDX ) if(trueIDX == idx) good = true;
			if(good) continue;
			auto measuredMuon = getMeasuredReco(idx);
 			plotter.getOrMake1D(TString::Format("%s_fs_incl_pt",glbPrefix.Data()),";bkg. segment p_{T}",60,0,30)->Fill(  measuredMuon.pt() );
 			if(std::fabs(measuredMuon.eta()) < 2.4) continue;
 			plotter.getOrMake1D(TString::Format("%s_fs_pt",glbPrefix.Data()),";bkg. segment p_{T}",60,0,30)->Fill(  measuredMuon.pt() );

 			if(measuredMuon.pt() > 2) {
 				plotter.getOrMake1D(TString::Format("%s_fs_singleMu5_me0Mugt2_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt5SingleMu).mass() );
 				plotter.getOrMake1D(TString::Format("%s_fs_singleMu20_me0Mugt2_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt20SingleMu).mass() );

 			}
 			if(measuredMuon.pt() > 5) {
 				plotter.getOrMake1D(TString::Format("%s_fs_singleMu5_me0Mugt5_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt5SingleMu).mass() );
 				plotter.getOrMake1D(TString::Format("%s_fs_singleMu20_me0Mugt5_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt20SingleMu).mass() );

 			}




		}

	}

	virtual void runAEvent() {
		fakeStudies();
	}



};

#endif

void AnalyzeME0SegmentTree(std::string fileName, std::string prefix, std::string outName){

//	FindParams a (fileName,"Events");
//	a.glbPrefix = prefix;
//	a.analyze();
//	a.solve();
//	a.write(outName);

	TestParams a (fileName,"Events");
	a.glbPrefix = prefix;
	a.analyze();
	a.write(outName);
}

