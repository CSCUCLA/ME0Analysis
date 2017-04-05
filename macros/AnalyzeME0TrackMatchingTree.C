
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TVector2.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
using namespace std;
double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
double deltaR2(const float  eta1, float phi1,  const float eta2, float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2)*TVector2::Phi_mpi_pi(phi1 - phi2)  + (eta1-eta2)*(eta1-eta2);}

class Analyzer : public BaseTupleAnalyzer{
public:
	  TString glbPrefix = "";
	  unsigned int nP = 0;
	  unsigned int nS = 0;
	  double maxPhi2SC=0;
	  double phiC=0;
	  double dPhiC=0;
	  double etaC= 0;

  Analyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){
	  simMuon_pt             = new std::vector<float> ;
	  simMuon_eta            = new std::vector<float> ;
	  simMuon_phi            = new std::vector<float> ;
	  simMuon_q              = new std::vector<int>   ;
	  simMuon_segmentQuality = new std::vector<int>   ;
	  simMuon_segmentIDX     = new std::vector<int>   ;
		 simMuon_trackIDX    =new std::vector<int>   ;
		 simMuon_track_pt    =new std::vector<float> ;
		 simMuon_track_eta   =new std::vector<float> ;
		 simMuon_track_phi   =new std::vector<float> ;
		 me0Muon_p           =new std::vector<float> ;

	  me0Muon_pt             = new std::vector<float> ;
	  me0Muon_eta            = new std::vector<float> ;
	  me0Muon_phi            = new std::vector<float> ;
	  me0Muon_q              = new std::vector<int>   ;
	  me0Muon_trackIDX       = new std::vector<int>   ;
	  me0Muon_segIDX         = new std::vector<int>   ;
	  me0Muon_truthType      = new std::vector<int>   ;
	  me0Muon_track_eta      = new std::vector<float> ;
	  me0Muon_track_phi      = new std::vector<float> ;
	  me0Muon_track_dphi     = new std::vector<float> ;
	  me0Muon_track_deta     = new std::vector<float> ;
	  me0Muon_track_x        = new std::vector<float> ;
	  me0Muon_track_y        = new std::vector<float> ;
	  me0Muon_track_dx       = new std::vector<float> ;
	  me0Muon_track_dy       = new std::vector<float> ;
	  me0Muon_track_sigx     = new std::vector<float> ;
	  me0Muon_track_sigy     = new std::vector<float> ;
	  me0Muon_track_sigdx    = new std::vector<float> ;
	  me0Muon_track_sigdy    = new std::vector<float> ;
	  me0Muon_segment_eta    = new std::vector<float> ;
	  me0Muon_segment_phi    = new std::vector<float> ;
	  me0Muon_segment_dphi   = new std::vector<float> ;
	  me0Muon_segment_deta   = new std::vector<float> ;
	  me0Muon_segment_x      = new std::vector<float> ;
	  me0Muon_segment_y      = new std::vector<float> ;
	  me0Muon_segment_dx     = new std::vector<float> ;
	  me0Muon_segment_dy     = new std::vector<float> ;
	  me0Muon_segment_sigx   = new std::vector<float> ;
	  me0Muon_segment_sigy   = new std::vector<float> ;
	  me0Muon_segment_sigdx  = new std::vector<float> ;
	  me0Muon_segment_sigdy  = new std::vector<float> ;

	  setBranchAddress(     "simMuon_pt"              ,&simMuon_pt             );
	  setBranchAddress(     "simMuon_eta"             ,&simMuon_eta            );
	  setBranchAddress(     "simMuon_phi"             ,&simMuon_phi            );
	  setBranchAddress(     "simMuon_q"               ,&simMuon_q              );
	  setBranchAddress(     "simMuon_segmentQuality"  ,&simMuon_segmentQuality );
	  setBranchAddress(     "simMuon_segmentIDX"      ,&simMuon_segmentIDX     );
	  setBranchAddress(     "simMuon_trackIDX"        ,&simMuon_trackIDX );
	  setBranchAddress(     "simMuon_track_pt"        ,&simMuon_track_pt );
	  setBranchAddress(     "simMuon_track_eta"       ,&simMuon_track_eta);
	  setBranchAddress(     "simMuon_track_phi"       ,&simMuon_track_phi);

	  setBranchAddress(     "me0Muon_p"               ,&me0Muon_p              );
	  setBranchAddress(     "me0Muon_pt"              ,&me0Muon_pt             );
	  setBranchAddress(     "me0Muon_eta"             ,&me0Muon_eta            );
	  setBranchAddress(     "me0Muon_phi"             ,&me0Muon_phi            );
	  setBranchAddress(     "me0Muon_q"               ,&me0Muon_q              );
	  setBranchAddress(     "me0Muon_trackIDX"        ,&me0Muon_trackIDX       );
	  setBranchAddress(     "me0Muon_segIDX"          ,&me0Muon_segIDX         );
	  setBranchAddress(     "me0Muon_truthType"       ,&me0Muon_truthType      );
	  setBranchAddress(     "me0Muon_track_eta"       ,&me0Muon_track_eta      );
	  setBranchAddress(     "me0Muon_track_phi"       ,&me0Muon_track_phi      );
	  setBranchAddress(     "me0Muon_track_dphi"      ,&me0Muon_track_dphi     );
	  setBranchAddress(     "me0Muon_track_deta"      ,&me0Muon_track_deta     );
	  setBranchAddress(     "me0Muon_track_x"         ,&me0Muon_track_x        );
	  setBranchAddress(     "me0Muon_track_y"         ,&me0Muon_track_y        );
	  setBranchAddress(     "me0Muon_track_dx"        ,&me0Muon_track_dx       );
	  setBranchAddress(     "me0Muon_track_dy"        ,&me0Muon_track_dy       );
	  setBranchAddress(     "me0Muon_track_sigx"      ,&me0Muon_track_sigx     );
	  setBranchAddress(     "me0Muon_track_sigy"      ,&me0Muon_track_sigy     );
	  setBranchAddress(     "me0Muon_track_sigdx"     ,&me0Muon_track_sigdx    );
	  setBranchAddress(     "me0Muon_track_sigdy"     ,&me0Muon_track_sigdy    );
	  setBranchAddress(     "me0Muon_segment_eta"     ,&me0Muon_segment_eta    );
	  setBranchAddress(     "me0Muon_segment_phi"     ,&me0Muon_segment_phi    );
	  setBranchAddress(     "me0Muon_segment_dphi"    ,&me0Muon_segment_dphi   );
	  setBranchAddress(     "me0Muon_segment_deta"    ,&me0Muon_segment_deta   );
	  setBranchAddress(     "me0Muon_segment_x"       ,&me0Muon_segment_x      );
	  setBranchAddress(     "me0Muon_segment_y"       ,&me0Muon_segment_y      );
	  setBranchAddress(     "me0Muon_segment_dx"      ,&me0Muon_segment_dx     );
	  setBranchAddress(     "me0Muon_segment_dy"      ,&me0Muon_segment_dy     );
	  setBranchAddress(     "me0Muon_segment_sigx"    ,&me0Muon_segment_sigx   );
	  setBranchAddress(     "me0Muon_segment_sigy"    ,&me0Muon_segment_sigy   );
	  setBranchAddress(     "me0Muon_segment_sigdx"   ,&me0Muon_segment_sigdx  );
	  setBranchAddress(     "me0Muon_segment_sigdy"   ,&me0Muon_segment_sigdy  );






  }

  std::vector<float>  * simMuon_pt             ;
  std::vector<float>  * simMuon_eta            ;
  std::vector<float>  * simMuon_phi            ;
  std::vector<int>    * simMuon_q              ;
  std::vector<int>    * simMuon_segmentQuality ;
  std::vector<int>    * simMuon_segmentIDX     ;
  std::vector<int>    * simMuon_trackIDX       ;
  std::vector<float>  * simMuon_track_pt       ;
  std::vector<float>  * simMuon_track_eta      ;
  std::vector<float>  * simMuon_track_phi      ;
  std::vector<float>  * me0Muon_p              ;
  std::vector<float>  * me0Muon_pt             ;
  std::vector<float>  * me0Muon_eta            ;
  std::vector<float>  * me0Muon_phi            ;
  std::vector<int>    * me0Muon_q              ;
  std::vector<int>    * me0Muon_trackIDX       ;
  std::vector<int>    * me0Muon_segIDX         ;
  std::vector<int>    * me0Muon_truthType      ;
  std::vector<float>  * me0Muon_track_eta      ;
  std::vector<float>  * me0Muon_track_phi      ;
  std::vector<float>  * me0Muon_track_dphi     ;
  std::vector<float>  * me0Muon_track_deta     ;
  std::vector<float>  * me0Muon_track_x        ;
  std::vector<float>  * me0Muon_track_y        ;
  std::vector<float>  * me0Muon_track_dx       ;
  std::vector<float>  * me0Muon_track_dy       ;
  std::vector<float>  * me0Muon_track_sigx     ;
  std::vector<float>  * me0Muon_track_sigy     ;
  std::vector<float>  * me0Muon_track_sigdx    ;
  std::vector<float>  * me0Muon_track_sigdy    ;
  std::vector<float>  * me0Muon_segment_eta    ;
  std::vector<float>  * me0Muon_segment_phi    ;
  std::vector<float>  * me0Muon_segment_dphi   ;
  std::vector<float>  * me0Muon_segment_deta   ;
  std::vector<float>  * me0Muon_segment_x      ;
  std::vector<float>  * me0Muon_segment_y      ;
  std::vector<float>  * me0Muon_segment_dx     ;
  std::vector<float>  * me0Muon_segment_dy     ;
  std::vector<float>  * me0Muon_segment_sigx   ;
  std::vector<float>  * me0Muon_segment_sigy   ;
  std::vector<float>  * me0Muon_segment_sigdx  ;
  std::vector<float>  * me0Muon_segment_sigdy  ;

  std::vector<int>    goodMatches;
  int nTrueGoodMuons = 0;

  double matchQuantity(int idx)  {
//				  return deltaR2(me0Muon_track_eta->at(idx),me0Muon_track_phi->at(idx),me0Muon_segment_eta->at(idx),me0Muon_segment_phi->at(idx) );
	  const float dPhi =  std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
	  const float phi =  std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
	  const float eta =  std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));

	  const float xres =  std::fabs(me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx)) /
			  std::sqrt( me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)  + me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)   );

	  const float yres =  std::fabs(me0Muon_track_y->at(idx) - me0Muon_segment_y->at(idx)) /
			  std::sqrt( me0Muon_track_sigy->at(idx)*me0Muon_track_sigy->at(idx)  + me0Muon_segment_sigy->at(idx)*me0Muon_segment_sigy->at(idx)   );

	  const float dxres =  std::fabs(me0Muon_track_dx->at(idx) - me0Muon_segment_dx->at(idx)) /
			  std::sqrt( me0Muon_track_sigdx->at(idx)*me0Muon_track_sigdx->at(idx)  + me0Muon_segment_sigdx->at(idx)*me0Muon_segment_sigdx->at(idx)   );

	  const float dyres =  std::fabs(me0Muon_track_dy->at(idx) - me0Muon_segment_dy->at(idx)) /
			  std::sqrt( me0Muon_track_sigdy->at(idx)*me0Muon_track_sigdy->at(idx)  + me0Muon_segment_sigdy->at(idx)*me0Muon_segment_sigdy->at(idx)   );

	  return (10*dPhi) + phi + eta;
  }

void doMatching(bool applyMatching = false) {
	  goodMatches.clear();
	  nTrueGoodMuons =0;
  	std::vector<pair<float,int>> allMatches;
  	std::vector<int> usedSegments;
  	std::vector<int> usedTracks;

  	for(unsigned int iM = 0; iM < me0Muon_p->size(); ++iM ){
  		if(applyMatching && !isAMatch(iM)) continue;
//  		if(std::fabs(me0Muon_segment_dphi->at(iM)) > 0.013) continue;
  		allMatches.emplace_back(matchQuantity(iM),iM);
  	}
	std::sort(allMatches.begin(),allMatches.end(),
			[](const std::pair<float,int> & a,
					const std::pair<float,int> & b) { return a.first < b.first;});

	for(unsigned int iM =0; iM < allMatches.size(); ++iM){
		int matchIDX = allMatches[iM].second;
		bool isUnique = true;
		for(unsigned int iS = 0; iS < usedSegments.size(); ++iS){
			if(usedSegments[iS] == me0Muon_segIDX->at(matchIDX) ) {
				isUnique = false;
				break;
			}
		}
		for(unsigned int iT = 0; iT < usedTracks.size(); ++iT){
			if(usedTracks[iT] == me0Muon_trackIDX->at(matchIDX) ) {
				isUnique = false;
				break;
			}
		}
		if(!isUnique) continue;
		goodMatches.push_back(matchIDX);
		usedSegments.push_back(me0Muon_segIDX->at(matchIDX));
		usedTracks.push_back(me0Muon_trackIDX->at(matchIDX));
		if(me0Muon_truthType->at(matchIDX) == 0 ){
			bool found = false;
			bool goodPT = false;
			for(unsigned int iSM = 0; iSM < simMuon_pt->size(); ++iSM){
					if(simMuon_segmentIDX->at(iSM) != me0Muon_segIDX->at(matchIDX)) continue;
					if(simMuon_trackIDX->at(iSM) != me0Muon_trackIDX->at(matchIDX)) continue;
					found = true;
					if(std::fabs(simMuon_track_pt->at(iSM)/simMuon_pt->at(iSM) - 1) <= 0.2){
						goodPT = true;
					}
					break;
			}
			if(!found) cout << endl << "SOMETHING IS WRONG!"<< endl;
			if(goodPT) nTrueGoodMuons++;
		}
	}
  };

bool isAMatch(int idx, bool doExtra = false){
	  if(me0Muon_pt->at(idx) < 1) return false;
	  const double dPhi = std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
	  const double phi  = std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
	  const double eta  = std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));



	  const double dPhiS2C  = std::min(std::max(.3/me0Muon_p->at(idx),maxPhi2SC),dPhiC);
	  const double phiS2C  = std::min(std::max(1.2/me0Muon_p->at(idx),maxPhi2SC),phiC);


	  if(doExtra)
		  return eta < etaC && phi < phiC && dPhi < dPhiC && dPhi < dPhiS2C && phi < phiS2C ;
	  else
		  return eta < etaC && phi < phiC && dPhi < dPhiC;
}
  void makeOneDCutPlots(const int idx, TString prefix ){

	  auto makePlotGrp= [&] (TString name,TString title,TString unit, int nBins,float varMax,  float value, float error) {
		    plotter.getOrMake1D(TString::Format("%s_resid_%s",prefix.Data(),name.Data())    ,TString::Format(";track - segment %s %s;a.u."      , title.Data(), unit.Data()),nBins,-1*varMax ,varMax)->Fill(value);
		    plotter.getOrMake1D(TString::Format("%s_abs_resid_%s",prefix.Data(),name.Data()),TString::Format(";|track - segment %s| %s;a.u."    , title.Data(), unit.Data()),nBins,0   ,varMax)->Fill(std::fabs(value));
		    if(error < 0) return;
		    plotter.getOrMake1D(TString::Format("%s_pull_%s",prefix.Data(),name.Data())     ,TString::Format(";track - segment %s / #sigma;a.u."  , title.Data())           ,nBins,-10 ,10)->Fill(value/error);
		    plotter.getOrMake1D(TString::Format("%s_abs_pull_%s",prefix.Data(),name.Data()) ,TString::Format(";|track - segment %s / #sigma|;a.u.", title.Data())           ,nBins,0   ,10)->Fill(std::fabs(value/error ));
	  };
	  makePlotGrp("x","x","[cm]",1000,15,me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx), std::sqrt(me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)+me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)) );

	  makePlotGrp("xs","xs","[cm]",1000,30,me0Muon_pt->at(idx)*(me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx)), std::sqrt(me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)+me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)) );

	  makePlotGrp("y","y","[cm]",1000,15,me0Muon_track_y->at(idx) - me0Muon_segment_y->at(idx), std::sqrt(me0Muon_track_sigy->at(idx)*me0Muon_track_sigy->at(idx)+me0Muon_segment_sigy->at(idx)*me0Muon_segment_sigy->at(idx)) );
	  makePlotGrp("dx","#Deltax","",1000,5,me0Muon_track_dx->at(idx) - me0Muon_segment_dx->at(idx), std::sqrt(me0Muon_track_sigdx->at(idx)*me0Muon_track_sigdx->at(idx)+me0Muon_segment_sigdx->at(idx)*me0Muon_segment_sigdx->at(idx)) );
	  makePlotGrp("dy","#Deltay","",1000,5,me0Muon_track_dy->at(idx) - me0Muon_segment_dy->at(idx), std::sqrt(me0Muon_track_sigdy->at(idx)*me0Muon_track_sigdy->at(idx)+me0Muon_segment_sigdy->at(idx)*me0Muon_segment_sigdy->at(idx)) );

	  makePlotGrp("eta","#eta","",1000,.5,me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx), -1 );
	  makePlotGrp("phi","#phi","",1000,.5,me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx), -1 );
	  makePlotGrp("deta","#Delta#eta","",1000,.5,me0Muon_track_deta->at(idx) - me0Muon_segment_deta->at(idx), -1 );
	  makePlotGrp("dphi","#Delta#phi","",1000,.5,me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx), -1 );

	  makePlotGrp("seg_dphi","(seg) #Delta#phi","",1000,.05,me0Muon_segment_dphi->at(idx), -1 );

	  double alpha = .019;
	  double dPhiS = std::max(0.007, alpha/ me0Muon_pt->at(idx) ) ;

	  double valueDPS1 = (me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx))* me0Muon_pt->at(idx);
	  if(std::fabs((me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx))) > 0.0065) valueDPS1 = 2.0;
	  if(std::fabs((me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx))) < .0014 ) valueDPS1 = 0.0;
	  makePlotGrp("dphiS1","#Delta#phi S1","",1000,1.0,  valueDPS1, -1 );

	  double valueDPS2 = (me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx))* me0Muon_p->at(idx);
	  if(std::fabs((me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx))) > 0.0065) valueDPS2 = 2.0;
	  makePlotGrp("dphiS2","#Delta#phi S2","",1000,2.0,  valueDPS2, -1 );

	  double valuePS1 = (me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))* me0Muon_pt->at(idx);
	  if(std::fabs((me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))) > 0.05) valuePS1 = 2.0;
	  if(std::fabs((me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))) < .0014 ) valuePS1 = 0.0;
	  makePlotGrp("phiS1","#phi S1","",1000,1.0,  valuePS1, -1 );

	  double valuePS2 = (me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))* me0Muon_p->at(idx);
	  if(std::fabs((me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))) > 0.05) valuePS2 = 2.0;
	  if(std::fabs((me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))) < .0014 ) valuePS2 = 0.0;


	  makePlotGrp("phiS2","#phi S2","",1000,2.0,  valuePS2, -1 );




	  makePlotGrp("phis","#phi s","",1000,1,(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx))* me0Muon_pt->at(idx), -1 );

	  plotter.getOrMake2D(TString::Format("%s_phixdphi",prefix.Data()) ,";phi; dphi",100,-.1,.1,100,-.01,.01)->
			  					  Fill(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx),me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));



  }

  void oneDCutPlots(TString prefix){
	  if(nTrueGoodMuons != 2) return;
	  //eff
	  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM ){
		  const int me0MIDX = goodMatches[iRM];
		  TString muonStr = prefix + "_";

		  muonStr += me0Muon_truthType->at(me0MIDX) == 0 ? "signal" : "fake";
		  muonStr += "_";
		  float pt = me0Muon_pt->at(me0MIDX);
		  if(pt < 2){
			  muonStr += "pt_eq0p5to2";
		  } else if (pt < 3){
			  muonStr += "pt_eq2to3";
		  } else if (pt < 5){
			  muonStr += "pt_eq3to5";
		  } else if (pt < 10){
			  muonStr += "pt_eq5to10";
		  } else if (pt < 20){
			  muonStr += "pt_eq10to20";
		  } else muonStr += "pt_geq20";
		  makeOneDCutPlots(me0MIDX,muonStr);
		  if(std::fabs(me0Muon_track_dphi->at(me0MIDX) - me0Muon_segment_dphi->at(me0MIDX)) < 0.01){
			  makeOneDCutPlots(me0MIDX,muonStr+"_absDPhilt0p1");
			  if(std::fabs(me0Muon_track_eta->at(me0MIDX) - me0Muon_segment_eta->at(me0MIDX)) < 0.1 && std::fabs(me0Muon_track_phi->at(me0MIDX) - me0Muon_segment_phi->at(me0MIDX)) < .1){
				  makeOneDCutPlots(me0MIDX,muonStr+"_absDPhilt0p1_absPhiEtaLt0p1");
			  }
		  }
	  }

  }

  void doEffList(const int idx, const float pt, TString prefix) {
	  double dPhi = std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
	  double phi  = std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
	  double eta  = std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));

	  double dPhiS2C  = std::min(std::max(.3/me0Muon_p->at(idx),maxPhi2SC),dPhiC);
	  double phiS2C  = std::min(std::max(1.2/me0Muon_p->at(idx),maxPhi2SC),phiC);


	  auto doInt = [&] (TString name, float pt){
		  	if(pt > 2)  plotter.getOrMake1D(TString::Format("%s_%s_Int",prefix.Data(),name.Data()),";inner track p_{T}; a.u.",4,-.5,3.5)->Fill(0);
		  	if(pt > 3)  plotter.getOrMake1D(TString::Format("%s_%s_Int",prefix.Data(),name.Data()),";inner track p_{T}; a.u.",4,-.5,3.5)->Fill(1);
		  	if(pt > 5)  plotter.getOrMake1D(TString::Format("%s_%s_Int",prefix.Data(),name.Data()),";inner track p_{T}; a.u.",4,-.5,3.5)->Fill(2);
		  	if(pt > 10) plotter.getOrMake1D(TString::Format("%s_%s_Int",prefix.Data(),name.Data()),";inner track p_{T}; a.u.",4,-.5,3.5)->Fill(3);
	  };

	    plotter.getOrMake1D(TString::Format("%s_noCuts_pt",prefix.Data()),";inner track p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    if(eta < etaC){
	    	plotter.getOrMake1D(TString::Format("%s_dEta_lt0p08_pt",prefix.Data()),";inner track p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    	doInt("eta",pt);
	    }
	    if(eta < etaC && phi < phiC){
	    	plotter.getOrMake1D(TString::Format("%s_eta_lt0p08_phi_lt0p05_pt",prefix.Data()),";p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    	doInt("eta_phi",pt);
	    }
	    if(eta < etaC && phi < phiC && dPhi < dPhiC ){
	    	plotter.getOrMake1D(TString::Format("%s_eta_lt0p08_phi_lt0p05_dphi_lt0p0065_pt",prefix.Data()),";inner track p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    	doInt("eta_phi_dPhi",pt);
	    }
	    if(eta < etaC && phi < phiC && dPhi < dPhiC && phi < phiS2C ){
	    	plotter.getOrMake1D(TString::Format("%s_eta_lt0p08_phi_lt0p05_dphi_lt0p0065_tP_pt",prefix.Data()),";inner track p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    	doInt("eta_phi_dPhi_S1",pt);
	    }
	    if(eta < etaC && phi < phiC && dPhi < dPhiC  && phi < phiS2C && dPhi < dPhiS2C ){
	    	plotter.getOrMake1D(TString::Format("%s_eta_lt0p08_phi_lt0p05_dphi_lt0p0065_tPtDP_pt",prefix.Data()),";inner track p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
	    	doInt("eta_phi_dPhi_S2",pt);
	    }
  }

  void effCuts(TString prefix) {
	  //eff
	  int nSM = 0;
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
		  float pt = simMuon_pt->at(iM);
		    plotter.getOrMake1D(TString::Format("%s_real_muon_incl_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
		    bool goodTrackPT = goodTrack && std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) <= 0.2;
		    int  goodNearSegment = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearSegment = iRM;
		    	break;

		    }
		    bool passCut = goodNearSegment >= 0 && isAMatch(goodNearSegment);
		    bool goodMatch = false;
		    if(goodNearSegment >= 0)
		    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
		    	if(goodMatches[iRM] != goodNearSegment) continue;
		    	goodMatch = true;
		    	break;
		    }
		    if(goodTrackPT) plotter.getOrMake1D(TString::Format("%s_real_muon_goodTrackPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrackPT && goodSegment) plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegmentAndGoodPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrackPT && passCut)     plotter.getOrMake1D(TString::Format("%s_real_muon_passCutsAndGoodPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrackPT && passCut && goodMatch)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodME0Muon_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);

//		    if(goodTrackPT && passCut && !goodMatch && pt > 10){
//		    					  cout <<"M :(" << simMuon_pt->at(iM) <<","<< simMuon_eta->at(iM) <<","<< simMuon_phi->at(iM)<<
//		    							  simMuon_segmentIDX->at(iM) <<" "<< simMuon_trackIDX->at(iM) <<" "<<simMuon_segmentQuality->at(iM) <<endl;
//
//	    						  cout <<"C "<< goodNearSegment <<" : ("<<me0Muon_pt->at(goodNearSegment) <<","<<me0Muon_eta->at(goodNearSegment) <<","<<me0Muon_phi->at(goodNearSegment) <<") -> ("
//	    								   <<me0Muon_track_eta->at(goodNearSegment) <<","<<me0Muon_track_phi->at(goodNearSegment) <<") (" <<me0Muon_segment_eta->at(goodNearSegment) <<","<<me0Muon_segment_phi->at(goodNearSegment) <<") -> "
//	    								   << me0Muon_segIDX->at(goodNearSegment) <<endl;
//
//		    					  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
//		    					    	if((me0Muon_segIDX->at(goodMatches[iRM]) != simMuon_segmentIDX->at(iM)) &&
//		    					    	(me0Muon_trackIDX->at(goodMatches[iRM]) != simMuon_trackIDX->at(iM))) continue;
//
//		    						  cout <<"G "<< goodMatches[iRM] <<" : ("<<me0Muon_pt->at(goodMatches[iRM]) <<","<<me0Muon_eta->at(goodMatches[iRM]) <<","<<me0Muon_phi->at(goodMatches[iRM]) <<") -> ("
//		    								   <<me0Muon_track_eta->at(goodMatches[iRM]) <<","<<me0Muon_track_phi->at(goodMatches[iRM]) <<") (" <<me0Muon_segment_eta->at(goodMatches[iRM]) <<","<<me0Muon_segment_phi->at(goodMatches[iRM]) <<") -> "
//		    								   << me0Muon_segIDX->at(goodMatches[iRM]) <<endl;
//
//		    					  }
//		    }
//

		    if(goodSegment) plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegment_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrackPT && goodMatch) doEffList(goodNearSegment,me0Muon_pt->at(goodNearSegment),prefix + "_real_muon");
		    if(goodTrackPT && goodMatch) nSM++;


	  }
if(true){
//		  if(nSM==2){
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForFakes",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);
		  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
			  if(me0Muon_truthType->at(goodMatches[iRM]) == 0) continue;
			  doEffList(goodMatches[iRM],me0Muon_pt->at(goodMatches[iRM]),prefix + "_fake_muon");
		  }

	  }


  }

  void makeAnnaPlot(TString prefix) {
	  //eff
	  int nSM = 0;
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){

		    bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
		    bool goodTrackPT = goodTrack && std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) <= 0.2;
		    int  goodNearSegment = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearSegment = iRM;
		    	break;

		    }
		    bool passCut = goodNearSegment >= 0 && isAMatch(goodNearSegment);
		    bool goodMatch = false;
		    if(goodNearSegment >= 0)
		    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
		    	if(goodMatches[iRM] != goodNearSegment) continue;
		    	goodMatch = true;
		    	break;
		    }
		    if(!goodTrackPT) continue;
		    float pt = simMuon_track_pt->at(iM);
		    if(pt < 2) continue;
		    plotter.getOrMake1D(TString::Format("%s_real_muon_passTrack_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    if(goodSegment) plotter.getOrMake1D(TString::Format("%s_real_muon_passSegment_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    if(passCut && goodMatch)   plotter.getOrMake1D(TString::Format("%s_real_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    if(passCut && goodMatch && isAMatch(goodNearSegment,true) )   plotter.getOrMake1D(TString::Format("%s_real_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    if(goodTrackPT && goodMatch) nSM++;


	  }
if(true){
//		  if(nSM==2){
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForFakesA",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);
		  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
			  if((me0Muon_pt->at(goodMatches[iRM]) < 2)) continue;
			  if(me0Muon_truthType->at(goodMatches[iRM]) == 0) continue;
			  if(!isAMatch(goodMatches[iRM])) continue;
			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(me0Muon_pt->at(goodMatches[iRM]));
			  if(isAMatch(goodMatches[iRM], true)) plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(me0Muon_pt->at(goodMatches[iRM]));
		  }

	  }
  }

//  void checkUniqueness(TString prefix) {
//	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
//		  float pt = simMuon_pt->at(iM);
//		    if(pt < 2) continue;
//		    bool goodTrack   = simMuon_foundTrack->at(iM);
//		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
//		    bool goodME0Muon = simMuon_ME0IDX->at(iM) >= 0 && me0Muon_truthType->at(simMuon_ME0IDX->at(iM)) == 0;
//
//		    if(simMuon_ME0IDX->at(iM) >= 0){
//				  plotter.getOrMake1D(TString::Format("%s_real_muon_ptMatch",prefix.Data()),";pt track - sim; a.u.",1000,-100,100)->Fill( me0Muon_pt->at(simMuon_ME0IDX->at(iM)) - simMuon_pt->at(iM) ) ;
//				  plotter.getOrMake2D(TString::Format("%s_real_muon_eta_phi_match",prefix.Data()),";eta track - sim;phi track - sim",1000,-2,2,1000,-2,2)->Fill( me0Muon_eta->at(simMuon_ME0IDX->at(iM)) - simMuon_eta->at(iM), me0Muon_phi->at(simMuon_ME0IDX->at(iM)) - simMuon_phi->at(iM) ) ;
//
//		    }
//
//		    if(!goodME0Muon) continue;
//
//
//
//
//			  TString muonStr = "";
//			  if (pt < 3){
//				  muonStr += "pt_eq2to3";
//			  } else if (pt < 5){
//				  muonStr += "pt_eq3to5";
//			  } else if (pt < 20){
//				  muonStr += "pt_eq5to20";
//			  } else muonStr += "pt_geq20";
//
//			  int nMatches = 0;
//			  bool isClosest = true;
//
//			  float closestDPHi = -1;
//			  float closestPhi = -1;
//			  float closestETA = -1;
//
//			  auto matchDR2 = [&] (int idx) -> double {
////				  return deltaR2(me0Muon_track_eta->at(idx),me0Muon_track_phi->at(idx),me0Muon_segment_eta->at(idx),me0Muon_segment_phi->at(idx) );
//				  const float dPhi =  std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
//				  const float phi =  std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
//				  const float eta =  std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));
//
//				  const float xres =  std::fabs(me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx)) /
//						  std::sqrt( me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)  + me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)   );
//
//				  const float yres =  std::fabs(me0Muon_track_y->at(idx) - me0Muon_segment_y->at(idx)) /
//						  std::sqrt( me0Muon_track_sigy->at(idx)*me0Muon_track_sigy->at(idx)  + me0Muon_segment_sigy->at(idx)*me0Muon_segment_sigy->at(idx)   );
//
//				  const float dxres =  std::fabs(me0Muon_track_dx->at(idx) - me0Muon_segment_dx->at(idx)) /
//						  std::sqrt( me0Muon_track_sigdx->at(idx)*me0Muon_track_sigdx->at(idx)  + me0Muon_segment_sigdx->at(idx)*me0Muon_segment_sigdx->at(idx)   );
//
//				  const float dyres =  std::fabs(me0Muon_track_dy->at(idx) - me0Muon_segment_dy->at(idx)) /
//						  std::sqrt( me0Muon_track_sigdy->at(idx)*me0Muon_track_sigdy->at(idx)  + me0Muon_segment_sigdy->at(idx)*me0Muon_segment_sigdy->at(idx)   );
//
////				  return (10*dPhi)*(10*dPhi) + phi*phi + eta*eta/4;
//				  return (10*dPhi) + phi + eta;
////				  return phi;
////				  return xres*xres + dxres*dxres;
//			  };
//
//
//			  int simME0MIDX = simMuon_ME0IDX->at(iM);
//			  double signalDR = matchDR2(simME0MIDX);
//			  const float sigDPhi =  std::fabs(me0Muon_track_dphi->at(simME0MIDX) - me0Muon_segment_dphi->at(simME0MIDX));
//			  const float sigPhi =  std::fabs(me0Muon_track_phi->at(simME0MIDX) - me0Muon_segment_phi->at(simME0MIDX));
//			  const float sigEta =  std::fabs(me0Muon_track_eta->at(simME0MIDX) - me0Muon_segment_eta->at(simME0MIDX));
//
//			  if(sigDPhi > 0.0065) continue;
//			  if(sigPhi > 0.05) continue;
//			  if(sigEta > 0.08) continue;
//
//			  if( me0Muon_pt->at(simME0MIDX)/simMuon_pt->at(iM) > 1.2 || me0Muon_pt->at(simME0MIDX)/simMuon_pt->at(iM) < 0.8  ) continue;
//
//			  for(unsigned int iR = 0; iR < me0Muon_pt->size(); ++iR ){
//				  if(iR == simME0MIDX) continue;
//				  if(simMuon_segmentIDX->at(iM) != me0Muon_segIDX->at(iR)) continue;
//				  const float bkgDPhi =  std::fabs(me0Muon_track_dphi->at(iR) - me0Muon_segment_dphi->at(iR));
//				  const float bkgPhi =  std::fabs(me0Muon_track_phi->at(iR) - me0Muon_segment_phi->at(iR));
//				  const float bkgEta =  std::fabs(me0Muon_track_eta->at(iR) - me0Muon_segment_eta->at(iR));
//				  if(bkgDPhi > 0.0065) continue;
//				  if(bkgPhi > 0.05) continue;
//				  if(bkgEta > 0.08) continue;
//
//				  plotter.getOrMake2D(TString::Format("%s_%s_real_muon_bkg_phixdphi",prefix.Data(), muonStr.Data()) ,";bkg phi; bkg dphi",100,-.1,.1,100,-.01,.01)->
//				  					  Fill(me0Muon_track_phi->at(iR) - me0Muon_segment_phi->at(iR),me0Muon_track_dphi->at(iR) - me0Muon_segment_dphi->at(iR));
//
//				  nMatches++;
//				  if(closestDPHi <= 0 ||bkgDPhi < closestDPHi  ) closestDPHi = bkgDPhi;
//				  if(closestPhi <= 0 ||bkgPhi < closestPhi  ) closestPhi = bkgPhi;
//				  if(closestETA <= 0 ||bkgEta < closestETA  ) closestETA = bkgEta;
//				  if(matchDR2(iR) < signalDR) isClosest = false;
//
//			  }
//
//			  if(!isClosest){
//				  cout <<"M :(" << simMuon_pt->at(iM) <<","<< simMuon_eta->at(iM) <<","<< simMuon_phi->at(iM)<<") -> (" << sigDPhi<<":"<<sigPhi<<":"<<sigEta<<" ) -> "<<
//						  simMuon_segmentIDX->at(iM) <<" "<< simMuon_ME0IDX->at(iM) <<" "<<simMuon_segmentQuality->at(iM) <<endl;
//				  for(unsigned int iR = 0; iR < me0Muon_pt->size(); ++iR ){
//					  if(simMuon_segmentIDX->at(iM) != me0Muon_segIDX->at(iR)) continue;
//					  const float bkgDPhi =  std::fabs(me0Muon_track_dphi->at(iR) - me0Muon_segment_dphi->at(iR));
//					  const float bkgPhi =  std::fabs(me0Muon_track_phi->at(iR) - me0Muon_segment_phi->at(iR));
//					  const float bkgEta =  std::fabs(me0Muon_track_eta->at(iR) - me0Muon_segment_eta->at(iR));
//					  if(bkgDPhi > 0.0065) continue;
//					  if(bkgPhi > 0.05) continue;
//					  if(bkgEta > 0.08) continue;
//
//					  cout <<"T "<< iR <<" : ("<<me0Muon_pt->at(iR) <<","<<me0Muon_eta->at(iR) <<","<<me0Muon_phi->at(iR) <<") -> ("
//							   <<me0Muon_track_eta->at(iR) <<","<<me0Muon_track_phi->at(iR) <<") (" <<me0Muon_segment_eta->at(iR) <<","<<me0Muon_segment_phi->at(iR) <<") -> "
//							   << me0Muon_segIDX->at(iR) <<" ("<<  bkgDPhi<<":"<<bkgPhi<<":"<<bkgEta<<" )"<<endl;
//
//				  }
//			  }
//
//
//			  plotter.getOrMake1D(TString::Format("%s_%s_real_muon_numberOfMatches",prefix.Data(), muonStr.Data()),";# of matches; a.u.",100,-0.5,99.5)->Fill(nMatches);
//			  plotter.getOrMake1D(TString::Format("%s_%s_real_muon_isClosest",prefix.Data(), muonStr.Data()),";isClosest; a.u.",2,-0.5,1.5)->Fill(isClosest);
//
//			  plotter.getOrMake1D(TString::Format("%s_%s_real_muon_ptRatio",prefix.Data(), muonStr.Data()),";PT track/ pT muon; a.u.",600,-3,3)->Fill(me0Muon_pt->at(simME0MIDX)/simMuon_pt->at(iM));
//			 if(!isClosest) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_badMatch_ptRatio",prefix.Data(), muonStr.Data()),";PT track/ pT muon; a.u.",600,-3,3)->Fill(me0Muon_pt->at(simME0MIDX)/simMuon_pt->at(iM));
//			 if(isClosest) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_goodMatch_ptRatio",prefix.Data(), muonStr.Data()),";PT track/ pT muon; a.u.",600,-3,3)->Fill(me0Muon_pt->at(simME0MIDX)/simMuon_pt->at(iM));
//
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_bkg_dphi",prefix.Data(), muonStr.Data()),";bkg dphi; a.u.",100,0,.1)->Fill(closestDPHi);
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_bkg_phi",prefix.Data(), muonStr.Data()),";bkg phi; a.u.",100,0,.1)->Fill(closestPhi);
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_bkg_eta",prefix.Data(), muonStr.Data()),";bkg eta; a.u.",100,0,.1)->Fill(closestETA);
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_sig_dphi",prefix.Data(), muonStr.Data()),";sig dphi; a.u.",100,0,.1)->Fill(sigDPhi);
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_sig_phi",prefix.Data(), muonStr.Data()) ,";sig phi; a.u.",100,0,.1)->Fill(sigPhi);
//			  if(nMatches) plotter.getOrMake1D(TString::Format("%s_%s_real_muon_sig_eta",prefix.Data(), muonStr.Data()) ,";sig eta; a.u.",100,0,.1)->Fill(sigEta);
//
//			  if(nMatches) plotter.getOrMake2D(TString::Format("%s_%s_real_muon_sig_phixdphi",prefix.Data(), muonStr.Data()) ,";sig phi; sig dphi",100,-.1,.1,100,-.01,.01)->
//					  Fill(me0Muon_track_phi->at(simME0MIDX) - me0Muon_segment_phi->at(simME0MIDX),me0Muon_track_dphi->at(simME0MIDX) - me0Muon_segment_dphi->at(simME0MIDX));
//	  }
//  }


  void quickTest(TString prefix) {
	  int nGMuons = 0;
	  for(unsigned int iM = 0; iM < simMuon_eta->size(); ++iM){
		  if(simMuon_pt->at(iM) < 2) continue;
		  if(simMuon_trackIDX->at(iM) < 0 ) continue;
		  if(std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) > 0.2) continue;
		  nGMuons++;
		  plotter.getOrMake1D(TString::Format("%s_real_muon_incl_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(simMuon_pt->at(iM));
		  bool found = false;
		  for(int iRM : goodMatches ) {
			  if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM) ) continue;
			  if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM) ) continue;
			  found = true;
			  break;
		  }
		  if(found) plotter.getOrMake1D(TString::Format("%s_real_muon_good_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(simMuon_pt->at(iM));
	  }

	  for(int iRM : goodMatches ) {
		  if(me0Muon_pt->at(iRM) < 2 ) continue;
		  if(me0Muon_truthType->at(iRM) == 0 ) continue;
		  plotter.getOrMake1D(TString::Format("%s_fake_muon_pt",prefix.Data()),";fake muon p_{T} [GeV]; a.u.",2,2,4)->Fill(std::min(me0Muon_pt->at(iRM),float(3.9)));
		  if(nGMuons == 2) plotter.getOrMake1D(TString::Format("%s_nG2_fake_muon_pt",prefix.Data()),";fake muon p_{T} [GeV]; a.u.",2,2,4)->Fill(std::min(me0Muon_pt->at(iRM),float(3.9)) );
	  }

	  plotter.getOrMake1D(TString::Format("%s_nEvent",prefix.Data()),";nEvents",2,0,2)->Fill(.5);
	  if(nGMuons == 2) plotter.getOrMake1D(TString::Format("%s_nEvent",prefix.Data()),";nEvents",2,0,2)->Fill(1.5);





  }

  virtual void runAEvent() {
	  doMatching(false);
	  oneDCutPlots(glbPrefix);
	  effCuts(glbPrefix);
	  makeAnnaPlot(glbPrefix);
//	  checkUniqueness(glbPrefix);
//	  quickTest(glbPrefix);
  }

  void write(TString fileName){ plotter.write(fileName);}
  HistGetter plotter;
};

#endif

void AnalyzeME0TrackMatchingTree(std::string fileName,unsigned int nP, unsigned int nS, std::string outFileName){
	Analyzer a(fileName,"Events");
	a.nP = nP;
	a.nS = nS;
	  a.maxPhi2SC = 2*.35/float(nS);
	  a.phiC = 0.05;
	  a.dPhiC = 0.0065;
	  a.etaC = nP >= 8 ? 0.06 : 0.08;

  a.glbPrefix = TString::Format("p%us%u",nP,nS);
  a.analyze();
  a.write(outFileName);
}
