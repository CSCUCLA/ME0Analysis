

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
using namespace std;
double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
double deltaR2(const float  eta1, float phi1,  const float eta2, float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2)*TVector2::Phi_mpi_pi(phi1 - phi2)  + (eta1-eta2)*(eta1-eta2);}
float getP(float pt, float eta, float phi){
	TVector3 mom;mom.SetXYZ(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*sinh(eta));
	return  mom.Mag();
};


class Analyzer : public BaseTupleAnalyzer{
public:
	  TString glbPrefix = "";
	  unsigned int nP = 0;
	  unsigned int nS = 0;
	  double alphaPhi2SC = 0;
	  double maxPhi2SC=0;
	  double alphaDPhi2SC = 0;
	  double maxDPhi2SC=0;
	  double phiC=0;
	  double dPhiC=0;
	  double etaC= 0;
	  bool isPureBKG = false;

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

     simMuon_gen_nLays  = new std::vector<int>   ;
     simMuon_gen_eta    = new std::vector<float> ;
     simMuon_gen_phi    = new std::vector<float> ;
     simMuon_gen_dphi   = new std::vector<float> ;
     simMuon_gen_deta   = new std::vector<float> ;
     simMuon_gen_x      = new std::vector<float> ;
     simMuon_gen_y      = new std::vector<float> ;
     simMuon_gen_dx     = new std::vector<float> ;
     simMuon_gen_dy     = new std::vector<float> ;


		 me0Muon_p           =new std::vector<float> ;

	  me0Muon_pt             = new std::vector<float> ;
	  me0Muon_eta            = new std::vector<float> ;
	  me0Muon_phi            = new std::vector<float> ;
	  me0Muon_q              = new std::vector<int>   ;
	  me0Muon_trackIDX       = new std::vector<int>   ;
	  me0Muon_segIDX         = new std::vector<int>   ;
	  me0Muon_trackTruthType   = new std::vector<int>   ;
	  me0Muon_segTrackTruthType= new std::vector<int>   ;
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

	  setBranchAddress("simMuon_gen_nLays", &simMuon_gen_nLays );
	  setBranchAddress("simMuon_gen_eta"  , &simMuon_gen_eta   );
	  setBranchAddress("simMuon_gen_phi"  , &simMuon_gen_phi   );
	  setBranchAddress("simMuon_gen_dphi" , &simMuon_gen_dphi  );
	  setBranchAddress("simMuon_gen_deta" , &simMuon_gen_deta  );
	  setBranchAddress("simMuon_gen_x"    , &simMuon_gen_x     );
	  setBranchAddress("simMuon_gen_y"    , &simMuon_gen_y     );
	  setBranchAddress("simMuon_gen_dx"   , &simMuon_gen_dx    );
	  setBranchAddress("simMuon_gen_dy"   , &simMuon_gen_dy    );

	  setBranchAddress(     "me0Muon_p"               ,&me0Muon_p              );
	  setBranchAddress(     "me0Muon_pt"              ,&me0Muon_pt             );
	  setBranchAddress(     "me0Muon_eta"             ,&me0Muon_eta            );
	  setBranchAddress(     "me0Muon_phi"             ,&me0Muon_phi            );
	  setBranchAddress(     "me0Muon_q"               ,&me0Muon_q              );
	  setBranchAddress(     "me0Muon_trackTruthType",&me0Muon_trackTruthType   );
	  setBranchAddress(     "me0Muon_segTrackTruthType",&me0Muon_segTrackTruthType);
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
  std::vector<int>    * simMuon_gen_nLays       ;
  std::vector<float>  * simMuon_gen_eta         ;
  std::vector<float>  * simMuon_gen_phi         ;
  std::vector<float>  * simMuon_gen_dphi        ;
  std::vector<float>  * simMuon_gen_deta        ;
  std::vector<float>  * simMuon_gen_x           ;
  std::vector<float>  * simMuon_gen_y           ;
  std::vector<float>  * simMuon_gen_dx          ;
  std::vector<float>  * simMuon_gen_dy          ;
  std::vector<float>  * me0Muon_p              ;
  std::vector<float>  * me0Muon_pt             ;
  std::vector<float>  * me0Muon_eta            ;
  std::vector<float>  * me0Muon_phi            ;
  std::vector<int>    * me0Muon_q              ;
  std::vector<int>    * me0Muon_trackTruthType    ;
  std::vector<int>    * me0Muon_segTrackTruthType ;
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
//			if(goodPT) nTrueGoodMuons++;
			nTrueGoodMuons++;
		}
	}
  };

bool isAMatch(int idx, bool doExtra = false){
	  if(me0Muon_pt->at(idx) < 1) return false;
	  const double dPhi = std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
	  const double phi  = std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
	  const double eta  = std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));

	  const double dPhiS2C  = std::min(std::max(alphaDPhi2SC/me0Muon_p->at(idx),maxDPhi2SC),dPhiC);
	  const double phiS2C  = std::min(std::max(alphaPhi2SC/me0Muon_p->at(idx),maxPhi2SC),phiC);


	  if(doExtra)
		  return eta < etaC && phi < phiC && dPhi < dPhiC && dPhi < dPhiS2C && phi < phiS2C ;
	  else
		  return eta < etaC && phi < phiC && dPhi < dPhiC;

//	  double rotDPhi =  diag.transA(dPhi,phi);
//	  double rotPhi =  diag.transB(dPhi,phi);
//	  double decPhi = std::sqrt(rotDPhi*rotDPhi+rotPhi*rotPhi);
//	  	  const double decphiS2C  = std::min(std::max(100.0/me0Muon_p->at(idx),1.5),3.8);
////	  	if(doExtra)
////	  	  return eta < etaC && decPhi < decphiS2C;
////	  	else
////	  		return eta < etaC && decPhi < 3.8;
//
//		  	if(doExtra)
//		  		return eta < etaC && decPhi < decphiS2C;
//		  	else
//		  		return eta < etaC && phi < phiC && dPhi < dPhiC && dPhi < dPhiS2C && phi < phiS2C;



}
  void makeOneDCutPlots(const int idx, TString prefix ){
	  double xBins[] = {2,3,5,10,15,20,30};
	  int nBinsX = 6;

	  auto makePlotGrp= [&] (TString name,TString title,TString unit, int nBins,float varMax,  float value, float error) {
		    plotter.getOrMake1D(TString::Format("%s_resid_%s",prefix.Data(),name.Data())    ,TString::Format(";track - segment %s %s;a.u."      , title.Data(), unit.Data()),nBins,-1*varMax ,varMax)->Fill(value);
		    plotter.getOrMake1D(TString::Format("%s_abs_resid_%s",prefix.Data(),name.Data()),TString::Format(";|track - segment %s| %s;a.u."    , title.Data(),unit.Data()),nBins,0   ,varMax)->Fill(std::fabs(value));
		    plotter.getOrMake2D(TString::Format("%s_resid_%s_byPT",prefix.Data(),name.Data())    ,TString::Format(";track - segment %s %s; track p_{T}"      , title.Data(), unit.Data()),nBins,-1*varMax ,varMax, nBinsX,xBins)->Fill(value,me0Muon_pt->at(idx));

		    if(error < 0) return;
		    plotter.getOrMake1D(TString::Format("%s_pull_%s",prefix.Data(),name.Data())     ,TString::Format(";track - segment %s / #sigma;a.u."  , title.Data())           ,nBins,-10 ,10)->Fill(value/error);
		    plotter.getOrMake1D(TString::Format("%s_abs_pull_%s",prefix.Data(),name.Data()) ,TString::Format(";|track - segment %s / #sigma|;a.u.", title.Data())           ,nBins,0   ,10)->Fill(std::fabs(value/error ));
		    plotter.getOrMake2D(TString::Format("%s_pull_%s_byPT",prefix.Data(),name.Data())    ,TString::Format(";|track - segment %s / #sigma; track p_{T}"      , title.Data()),nBins,-4 ,4,  nBinsX,xBins)->Fill(value/error,me0Muon_pt->at(idx));

	  };
	  makePlotGrp("x","x","[cm]",1000,15,me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx), std::sqrt(me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)+me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)) );

	  makePlotGrp("xs","xs","[cm]",1000,30,me0Muon_pt->at(idx)*(me0Muon_track_x->at(idx) - me0Muon_segment_x->at(idx)), std::sqrt(me0Muon_track_sigx->at(idx)*me0Muon_track_sigx->at(idx)+me0Muon_segment_sigx->at(idx)*me0Muon_segment_sigx->at(idx)) );

	  makePlotGrp("y","y","[cm]",1000,15,me0Muon_track_y->at(idx) - me0Muon_segment_y->at(idx), std::sqrt(me0Muon_track_sigy->at(idx)*me0Muon_track_sigy->at(idx)+me0Muon_segment_sigy->at(idx)*me0Muon_segment_sigy->at(idx)) );
	  makePlotGrp("dx","#Deltax","",1000,0.7,me0Muon_track_dx->at(idx) - me0Muon_segment_dx->at(idx), std::sqrt(me0Muon_track_sigdx->at(idx)*me0Muon_track_sigdx->at(idx)+me0Muon_segment_sigdx->at(idx)*me0Muon_segment_sigdx->at(idx)) );
	  makePlotGrp("dy","#Deltay","",1000,4,me0Muon_track_dy->at(idx) - me0Muon_segment_dy->at(idx), std::sqrt(me0Muon_track_sigdy->at(idx)*me0Muon_track_sigdy->at(idx)+me0Muon_segment_sigdy->at(idx)*me0Muon_segment_sigdy->at(idx)) );

	  makePlotGrp("eta","#eta","",1000,.3,me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx), -1 );
	  makePlotGrp("phi","#phi","",1000,.3,me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx), -1 );
	  makePlotGrp("deta","#Delta#eta","",1000,.5,me0Muon_track_deta->at(idx) - me0Muon_segment_deta->at(idx), -1 );
	  makePlotGrp("dphi","#Delta#phi","",1000,.03,me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx), -1 );

	  const double dPhi = (me0Muon_track_dphi->at(idx) -me0Muon_segment_dphi->at(idx));
	  const double  phi = (me0Muon_track_phi->at(idx) -me0Muon_segment_phi->at(idx));
	  double rotDPhi =  diag.transA(dPhi,phi);
	  double rotPhi =  diag.transB(dPhi,phi);
	  makePlotGrp("decPhi","dec#phi","",1000,10,std::sqrt(rotDPhi*rotDPhi+rotPhi*rotPhi), -1 );
	  double valuedPS2 = std::sqrt(rotDPhi*rotDPhi+rotPhi*rotPhi)* me0Muon_p->at(idx);
	  if(std::sqrt(rotDPhi*rotDPhi+rotPhi*rotPhi) > 3.8 ) valuedPS2 = 500;
	  if(std::sqrt(rotDPhi*rotDPhi+rotPhi*rotPhi) < 1.2 ) valuedPS2 = 0.0;
	  makePlotGrp("decPhiS2","dec #phi S2","",1000,300.0,  valuedPS2, -1 );

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
		  TString muonStrPT = muonStr;
		  muonStrPT += "_";

		  float pt = me0Muon_pt->at(me0MIDX);
		  if(pt < 2){
			  muonStrPT += "pt_eq0p5to2";
		  } else if (pt < 3){
			  muonStrPT += "pt_eq2to3";
		  } else if (pt < 5){
			  muonStrPT += "pt_eq3to5";
		  } else if (pt < 10){
			  muonStrPT += "pt_eq5to10";
		  } else if (pt < 20){
			  muonStrPT += "pt_eq10to20";
		  } else muonStrPT += "pt_geq20";
//		  makeOneDCutPlots(me0MIDX,muonStr);
		  makeOneDCutPlots(me0MIDX,muonStrPT);
//		  if(std::fabs(me0Muon_track_dphi->at(me0MIDX) - me0Muon_segment_dphi->at(me0MIDX)) < 0.01){
//			  makeOneDCutPlots(me0MIDX,muonStr+"_absDPhilt0p1");
//			  if(std::fabs(me0Muon_track_eta->at(me0MIDX) - me0Muon_segment_eta->at(me0MIDX)) < 0.1 && std::fabs(me0Muon_track_phi->at(me0MIDX) - me0Muon_segment_phi->at(me0MIDX)) < .1){
//				  makeOneDCutPlots(me0MIDX,muonStr+"_absDPhilt0p1_absPhiEtaLt0p1");
//			  }
//		  }
	  }

  }

  void doEffList(const int idx, const float pt, TString prefix) {
	  double dPhi = std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
	  double phi  = std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
	  double eta  = std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));

	  double dPhiS2C  = std::min(std::max(alphaDPhi2SC/me0Muon_p->at(idx),maxDPhi2SC),dPhiC);
	  double phiS2C  = std::min(std::max(alphaPhi2SC/me0Muon_p->at(idx),maxPhi2SC),phiC);


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

		    auto getP = [&] (float pt, float eta, float phi)->float {
		    	TVector3 mom;mom.SetXYZ(pt*TMath::Cos(phi), pt*TMath::Sin(phi), pt*sinh(eta));
		    	return  mom.Mag();
		    };

		  float pt = simMuon_pt->at(iM);
		    plotter.getOrMake1D(TString::Format("%s_real_muon_incl_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    plotter.getOrMake1D(TString::Format("%s_real_muon_incl_p",prefix.Data()),";true muon |p| [GeV]; a.u.",30,0,150)->
		    		Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)));

		    plotter.getOrMake2D(TString::Format("%s_real_muon_incl_p_v_eta",prefix.Data()),";true muon |p| [GeV];true muon |#eta|",150,0,150,4,1.8,2.8)->Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)),std::fabs(simMuon_eta->at(iM)));
		    plotter.getOrMake2D(TString::Format("%s_real_muon_incl_pt_v_eta",prefix.Data()),";true muon p_{T} [GeV];true muon |#eta|",60,0,30,4,1.8,2.8)->Fill(simMuon_pt->at(iM),std::fabs(simMuon_eta->at(iM)));

		    if(pt >= 3) plotter.getOrMake1D(TString::Format("%s_real_muon_incl_eta",prefix.Data()),";true muon |#eta|; a.u.",12,1.8,3.0)->Fill(std::fabs(simMuon_eta->at(iM)));
		    bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
//		    bool goodTrackPT = goodTrack && std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) <= 0.2;
		    int  goodNearSegment = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearSegment = iRM;
		    	break;

		    }
		    bool passCut = goodNearSegment >= 0 && isAMatch(goodNearSegment,true);
		    bool goodMatch = false;
		    if(goodNearSegment >= 0)
		    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
		    	if(goodMatches[iRM] != goodNearSegment) continue;
		    	goodMatch = true;
		    	break;
		    }

		    int goodNearAnySegment = -1;
		    int goodNearAnySegment_PC = -1;
		    int goodNearAnySegment_GM = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearAnySegment = iRM;
		    	if(!isAMatch(iRM,true)) continue;
		    	goodNearAnySegment_PC = iRM;

			    for(unsigned int iRM2 = 0; iRM2 < goodMatches.size(); ++iRM2){
				    	if(goodMatches[iRM2] != iRM) continue;
				    	goodNearAnySegment_GM = iRM2;
				    	break;
				    }

		    }

		    if(goodNearSegment >= 0){
				double ptsolved = (0.073359 -0.02116*std::fabs(me0Muon_segment_eta->at(goodNearSegment)))/std::fabs(me0Muon_segment_dphi->at(goodNearSegment));
			    plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegmentAndTrack_solvedPToPT",prefix.Data()),";segment p_{T} / true muon p_{T}; a.u.",100,0,10)->Fill(ptsolved/pt);

		    }
		    if(goodSegment)   plotter.getOrMake2D(TString::Format("%s_real_muon_goodSegment_pt_v_eta",prefix.Data()),";true muon p_{T} [GeV];true muon |#eta|",60,0,30,4,1.8,2.8)->Fill(simMuon_pt->at(iM),std::fabs(simMuon_eta->at(iM)));
		    if(goodSegment)   plotter.getOrMake2D(TString::Format("%s_real_muon_goodSegment_p_v_eta",prefix.Data()),";true muon |p| [GeV];true muon |#eta|",150,0,150,4,1.8,2.8)->Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)),std::fabs(simMuon_eta->at(iM)));
		    if(goodSegment)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegment_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodSegment)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegment_p",prefix.Data()),";true muon |p| [GeV]; a.u.",30,0,150)->
		    		Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)));

		    if(goodSegment && pt >= 3) plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegment_eta",prefix.Data()),";true muon |#eta|; a.u.",12,1.8,3.0)->Fill(std::fabs(simMuon_eta->at(iM)));







		    if(goodTrack)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodTrack_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
//		    if(goodTrackPT) plotter.getOrMake1D(TString::Format("%s_real_muon_goodTrackPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrack && goodSegment) plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegmentAndGoodPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrack && passCut)     plotter.getOrMake1D(TString::Format("%s_real_muon_goodMatch_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrack && passCut && goodMatch)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodME0Muon_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);




		    if(goodTrack )   plotter.getOrMake1D(TString::Format("%s_real_muon_goodTrack_p",prefix.Data()),";true muon |p| [GeV]; a.u.",30,0,150)->
		    		Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)));
		    if(goodTrack && pt >= 3) plotter.getOrMake1D(TString::Format("%s_real_muon_goodTrack_eta",prefix.Data()),";true muon |#eta|; a.u.",12,1.8,3.0)->Fill(std::fabs(simMuon_eta->at(iM)));



			    if(goodTrack && goodSegment)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegmentAndGoodPT_p",prefix.Data()),";true muon |p| [GeV]; a.u.",30,0,150)->
			    		Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)));
			    if(goodTrack && goodSegment && pt >= 3) plotter.getOrMake1D(TString::Format("%s_real_muon_goodSegmentAndGoodPT_eta",prefix.Data()),";true muon |#eta|; a.u.",12,1.8,3.0)->Fill(std::fabs(simMuon_eta->at(iM)));


		    if(goodTrack && passCut)   plotter.getOrMake1D(TString::Format("%s_real_muon_goodMatch_p",prefix.Data()),";true muon |p| [GeV]; a.u.",30,0,150)->
		    		Fill(getP(simMuon_pt->at(iM),simMuon_eta->at(iM),simMuon_phi->at(iM)));
		    if(goodTrack && passCut && pt >= 3) plotter.getOrMake1D(TString::Format("%s_real_muon_goodMatch_eta",prefix.Data()),";true muon |#eta|; a.u.",12,1.8,3.0)->Fill(std::fabs(simMuon_eta->at(iM)));




		    if(goodTrack && goodNearAnySegment >= 0) plotter.getOrMake1D(TString::Format("%s_real_muon_ns_goodSegmentAndGoodPT_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrack && goodNearAnySegment_PC >= 0)     plotter.getOrMake1D(TString::Format("%s_real_muon_ns_goodMatch_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);
		    if(goodTrack && goodNearAnySegment_GM >= 0)   plotter.getOrMake1D(TString::Format("%s_real_muon_ns_goodME0Muon_pt",prefix.Data()),";true muon p_{T} [GeV]; a.u.",60,0,30)->Fill(pt);



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

		    if(goodTrack && goodMatch) doEffList(goodNearSegment,me0Muon_pt->at(goodNearSegment),prefix + "_real_muon");
		    if(goodTrack && goodMatch) nSM++;


	  }
		  if(isPureBKG || nSM==2){
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForFakes",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);
		  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
			  if(me0Muon_truthType->at(goodMatches[iRM]) == 0) continue;
			  doEffList(goodMatches[iRM],me0Muon_pt->at(goodMatches[iRM]),prefix + "_fake_muon");
		  }
		  for(unsigned int iRM = 0; iRM < me0Muon_truthType->size(); ++iRM){
			  if(me0Muon_truthType->at(iRM) == 0) continue;
			  doEffList(iRM,me0Muon_pt->at(iRM),prefix + "_nodisamb_fake_muon");
		  }

	  }


  }

  void getNPossibleMatches(TString prefix){
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
		  const double pt = simMuon_pt->at(iM);
//		  if(pt < 3 || pt > 5) continue;
		  if(simMuon_segmentQuality->at(iM) == 3) continue;
		  if(simMuon_trackIDX->at(iM) < 0) continue;

		    int  goodNearSegment = -1;
		    int nOtherMatches = 0;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;

		    	if(me0Muon_trackIDX->at(iRM) == simMuon_trackIDX->at(iM)) goodNearSegment = iRM;
		    	else  if(isAMatch(iRM,true)) nOtherMatches++;
		    }
		    if(goodNearSegment < 0 || !isAMatch(goodNearSegment,true)) continue;

		    plotter.getOrMake1D(TString::Format("%s_real_muon_otherTracks",prefix.Data()),";# of bkg. tracks matched to muon; a.u.",100,-0.5,99.5)->Fill(nOtherMatches);
	  }

	  std::map<int,int> numberOfTracksPerSeg;
	  std::map<int,int> numberOfSegPerTrack;
	  for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		  if(!isAMatch(iRM,true)) continue;
		  numberOfTracksPerSeg[me0Muon_segIDX->at(iRM)]++;
		  numberOfSegPerTrack[me0Muon_trackIDX->at(iRM)]++;
	  }
	  for(auto& iS : numberOfTracksPerSeg){
		    plotter.getOrMake1D(TString::Format("%s_all_nTracksPerSeg",prefix.Data()),";# of tracks matched to a single segment; a.u.",100,-0.5,99.5)->Fill(iS.second);

	  }
	  for(auto& iS : numberOfSegPerTrack){
		    plotter.getOrMake1D(TString::Format("%s_all_nSegPerTrack",prefix.Data()),";# of segments matched to a single track; a.u.",100,-0.5,99.5)->Fill(iS.second);
	  }

  }
  void makeAnnaPlot(TString prefix) {
	  //eff
	  int nSM = 0;
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){

		    bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
//		    bool goodTrackPT = goodTrack && std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) <= 0.2;
		    int  goodNearSegment = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearSegment = iRM;
		    	break;

		    }
		    bool passCut = goodNearSegment >= 0 && isAMatch(goodNearSegment,true);
		    bool goodMatch = false;
		    if(goodNearSegment >= 0)
		    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
		    	if(goodMatches[iRM] != goodNearSegment) continue;
		    	goodMatch = true;
		    	break;
		    }
		    if(simMuon_trackIDX->at(iM) < 0) continue;
		    float pt = simMuon_track_pt->at(iM);
		    float eta = simMuon_track_eta->at(iM);
		    if(pt < 2) continue;
		    plotter.getOrMake1D(TString::Format("%s_real_muon_passTrack_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    plotter.getOrMake1D(TString::Format("%s_real_muon_passTrack_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_real_muon_pt3to5_passTrack_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));

		    if(goodSegment){
		    	if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_real_muon_pt3to5_passSegment_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passSegment_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passSegment_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    }
		    if(passCut && goodMatch){
		    	if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_real_muon_pt3to5_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    }
		    if(passCut && goodMatch && isAMatch(goodNearSegment,true) ){
		    	if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_real_muon_pt3to5_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
		    	plotter.getOrMake1D(TString::Format("%s_real_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
		    }
//		    if(goodTrackPT && goodMatch) nSM++;
		    if(goodMatch && passCut) nSM++;


	  }

		  if(isPureBKG || nSM==2){
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForFakesA",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);
//		  for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
//			  if((me0Muon_pt->at(goodMatches[iRM]) < 2)) continue;
//			  if(me0Muon_truthType->at(goodMatches[iRM]) == 0) continue;
//			  if(!isAMatch(goodMatches[iRM])) continue;
//
//			  const float pt = me0Muon_pt->at(goodMatches[iRM]);
//			  const float eta = me0Muon_eta->at(goodMatches[iRM]);
//
//
//			  if(pt >= 2) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq2_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(pt >= 3) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq3_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(pt >= 10) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq10_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(pt >= 20) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq20_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//
//
//			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(me0Muon_pt->at(goodMatches[iRM]));
//			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_fake_muon_pt3to5_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(isAMatch(goodMatches[iRM], true)){
//				  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(me0Muon_pt->at(goodMatches[iRM]));
//				  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//				  if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_fake_muon_pt3to5_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  }
//		  }
		  for(unsigned int idx = 0; idx < me0Muon_p->size(); ++idx){
			  if((me0Muon_pt->at(idx) < 2)) continue;
			  if(me0Muon_truthType->at(idx) == 0) continue;
			  if(!isAMatch(idx,true)) continue;

//			  if(me0Muon_pt->at(idx) > 10 ){
//				  cout << "---------------------------------------------"<<endl;
//				  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
//					  cout << simMuon_pt->at(iM) <<" "<< simMuon_eta->at(iM)<< " "<<simMuon_phi->at(iM) <<" "<< simMuon_segmentQuality->at(iM)
//						   << " " << simMuon_segmentIDX->at(iM) <<" "<< simMuon_trackIDX->at(iM)<<endl;
//
//			  }
//				  cout <<"-"<<endl;
//				  for(unsigned int iM = 0; iM < me0Muon_eta->size(); ++iM ){
//					  cout << me0Muon_pt->at(iM) <<" "<< me0Muon_eta->at(iM)<< " "<<me0Muon_phi->at(iM) <<" "<< me0Muon_truthType->at(iM)
//						   << " " << me0Muon_segIDX->at(iM) <<" "<< me0Muon_trackIDX->at(iM) <<" "<< isAMatch(iM,true) <<endl;
//
//			  }
//			  }


			    bool goodMatch = false;
			    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
			    	if(goodMatches[iRM] != idx) continue;
			    	goodMatch = true;
			    	break;
			    }

			  const float pt = me0Muon_pt->at(idx);
			  const float eta = me0Muon_eta->at(idx);
			  const float totMom = me0Muon_p->at(idx);


			  if(pt >= 2) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq2_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
			  if(pt >= 3) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq3_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
			  if(pt >= 5 ) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq5_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
			  if(pt >= 10) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq10_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
			  if(pt >= 20) plotter.getOrMake1D(TString::Format("%s_fake_muon_ptgeq20_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));


			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_p",prefix.Data()),";pixel track |p| [GeV]; ",30,0,150)->Fill(totMom);
			  if(goodMatch)plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_GM_p",prefix.Data()),";pixel track |p| [GeV]; ",30,0,150)->Fill(totMom);
			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
			  if(goodMatch)plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_GM_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);

			  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));

			  if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_fake_muon_pt3to5_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(isAMatch(goodMatches[iRM], true)){
//				  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(me0Muon_pt->at(goodMatches[iRM]));
//				  plotter.getOrMake1D(TString::Format("%s_fake_muon_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//				  if(pt >= 3 && pt < 5) plotter.getOrMake1D(TString::Format("%s_fake_muon_pt3to5_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  }
		  }


//		  for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
//			  if((me0Muon_pt->at(iRM) < 2)) continue;
//			  if(me0Muon_truthType->at(iRM) == 0) continue;
//			  if(!isAMatch(iRM)) continue;
//
//			  const float pt = me0Muon_pt->at(iRM);
//			  const float eta = me0Muon_eta->at(iRM);
//
//			  plotter.getOrMake1D(TString::Format("%s_nodisamb_fake_muon_passME0Muon_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
//			  plotter.getOrMake1D(TString::Format("%s_nodisamb_fake_muon_passME0Muon_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  if(isAMatch(iRM, true)){
//				  plotter.getOrMake1D(TString::Format("%s_nodisamb_fake_muon_passME0MuonExtra_pt",prefix.Data()),";pixel track p_{T} [GeV]; ",9,2,11)->Fill(pt);
//				  plotter.getOrMake1D(TString::Format("%s_nodisamb_fake_muon_passME0MuonExtra_eta",prefix.Data()),";pixel track |#eta|; ",12,1.8,3.0)->Fill(std::fabs(eta));
//			  }
//		  }



	  }
  }

  void getBKGComp(TString prefix) {
	  //eff
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForBComp",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);

		  for(unsigned int idx = 0; idx < me0Muon_p->size(); ++idx){
			  if((me0Muon_pt->at(idx) < 2)) continue;
			  if(!isAMatch(idx,true)) continue;
			    bool goodMatch = false;
			    for(unsigned int iRM = 0; iRM < goodMatches.size(); ++iRM){
			    	if(goodMatches[iRM] != idx) continue;
			    	goodMatch = true;
			    	break;
			    }

			  const float pt = me0Muon_pt->at(idx);
			  const float eta = me0Muon_eta->at(idx);
			  int type = 5;
			  if(me0Muon_trackTruthType->at(idx) == 0 && me0Muon_segTrackTruthType->at(idx) == 0 ) type = 0;
			  else if(me0Muon_trackTruthType->at(idx) == 0 && me0Muon_segTrackTruthType->at(idx) != 0 ) type = 1;
			  else if(me0Muon_segTrackTruthType->at(idx) == 0 ) type = 2;
			  else if(me0Muon_segTrackTruthType->at(idx) == 1 ) type = 3;
			  else if(me0Muon_segTrackTruthType->at(idx) == 2 ) type = 4;

			  plotter.getOrMake1D(TString::Format("%s_fake_muon_backgroundComp",prefix.Data()),";type; ",10,-0.5,9.5)->Fill(type);
		  }
  }

  void getVariableCompPlots(TString prefix) {




  }

  class Diagonalizer {
  public:
	  void add(double a, double b) { n++; sA+=a; sB +=b; sAA+=a*a; sAB+=a*b; sBB+=b*b;}
	  void calcCov() {
		  sigAA = sAA/n  - (sA/n)*(sA/n);
		  sigBB = sBB/n  - (sB/n)*(sB/n);
		  sigAB = sAB/n  - (sA/n)*(sB/n);
	  }
	  void calcEigen() {
		  double D = std::sqrt((sigAA-sigBB)*(sigAA-sigBB) + 4*sigAB*sigAB);
		  lA = (sigAA+sigBB)/2 +D/2;
		  lB = (sigAA+sigBB)/2 -D/2;
		  vA1 = 2*sigAB;
		  vA2 = (sigBB - sigAA + D);
		  double ANorm = std::sqrt(vA1*vA1 + vA2*vA2);
		  vA1 /= ANorm;
		  vA2 /= ANorm;
		  vB1 = 2*sigAB;
		  vB2 = (sigBB - sigAA - D);
		  double BNorm = std::sqrt(vB1*vB1 + vB2*vB2);
		  vB1 /= BNorm;
		  vB2 /= BNorm;
	  }
	  void calcTrans () {
		  double isLA = 1./std::sqrt(lA);
		  double isLB = 1./std::sqrt(lB);
		  W11 = isLA*vA1;
		  W12 = isLB*vB1;
		  W21 = isLA*vA2;
		  W22 = isLB*vB2;
	  }
	  void fillTrans(double w11, double w12, double w21, double w22) {
		  W11 = w11;
		  W12 = w12;
		  W21 = w21;
		  W22 = w22;
	  }
	  double transA(double a, double b) const {return  W11*a + W21*b;}
	  double transB(double a, double b) const {return  W12*a + W22*b;}
	  int n=0;
	  double sA=0;
	  double sB=0;
	  double sAA=0;
	  double sBB=0;
	  double sAB=0;
	  double sigAA=0;
	  double sigAB=0;
	  double sigBB=0;

	  double lA = 0;
	  double lB = 0;
	  double vA1 = 0;
	  double vA2 = 0;
	  double vB1 = 0;
	  double vB2 = 0;

	  double W11=0;
	  double W12=0;
	  double W21=0;
	  double W22=0;
  };

  Diagonalizer diags[6];
  Diagonalizer diag;

  void makeResPlot(TString prefix) {
	  //eff
	  int nSM = 0;
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
		  const float absETA = std::fabs(simMuon_eta->at(iM));
		  if(absETA < 2.0 || absETA > 2.8) continue;
//		  if(simMuon_q->at(iM) < 0) continue;

		    bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
		    bool goodSegment = simMuon_segmentQuality->at(iM) != 3;
//		    bool goodTrackPT = goodTrack && std::fabs(simMuon_track_pt->at(iM)/simMuon_pt->at(iM) - 1) <= 0.2;
		    int  goodNearSegment = -1;
		    for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
		    	if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
		    	if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
		    	goodNearSegment = iRM;
		    	break;
		    }

//		    if(!goodTrackPT) continue;
		    if(goodNearSegment < 0 ) continue;
		    const  float pt = simMuon_track_pt->at(iM);
		    const  float eta = simMuon_track_eta->at(iM);
		    if(pt < 2) continue;

		    TString pTSTRING = "all";
		    if( pt < 3) pTSTRING = "pteq2to3";
		    else if( pt < 5) pTSTRING = "pteq3to5";
		    else if( pt < 10) pTSTRING = "pteq5to10";
		    else if( pt < 15) pTSTRING = "pteq10to15";
		    else if( pt < 20) pTSTRING = "pteq15to20";
		    else pTSTRING = "ptgeq20";

		    int iP = 0;
		    if( pt < 3) iP = 0;
		    else if( pt < 5)  iP = 1;
		    else if( pt < 10) iP = 2;
		    else if( pt < 15) iP = 3;
		    else if( pt < 20) iP = 4;
		    else iP = 5;


		    auto fillRes =  [&](TString recoName, TString name, TString title, int nBinsX, float minX, float maxX, float genQuant, float recoQuant ) {
			    plotter.getOrMake1D(TString::Format("%s_%s_real_muon_%s_genmreco_%s",prefix.Data(), pTSTRING.Data(),recoName.Data(), name.Data()),TString::Format(";sim %s - reco %s",title.Data(),title.Data()),nBinsX,minX,maxX)->Fill(genQuant - recoQuant );
			    plotter.getOrMake1D(TString::Format("%s_%s_real_muon_%s_value_%s",prefix.Data(), pTSTRING.Data(),recoName.Data(), name.Data()),TString::Format(";sim %s",title.Data()),nBinsX,minX,maxX)->Fill(recoQuant);
			    plotter.getOrMake1D(TString::Format("%s_%s_real_muon_%s_res_%s",prefix.Data(), pTSTRING.Data(),recoName.Data(), name.Data()),TString::Format(";(sim %s - reco %s) / gen. %s",title.Data(),title.Data(),title.Data()),1000,-5,5)->Fill( genQuant == 0 ? 5.0 :  (genQuant - recoQuant)/genQuant );
			    if(recoName == "track") {
			    	plotter.getOrMake1D(TString::Format("%s_%s_real_muon_gen_value_%s",prefix.Data(), pTSTRING.Data(), name.Data()),TString::Format(";sim %s",title.Data()),nBinsX,minX,maxX)->Fill(genQuant);
			    }

		    };
		    fillRes("track","eta","#eta"         ,200,-.05,0.05  ,simMuon_gen_eta->at(iM) ,me0Muon_track_eta ->at(goodNearSegment)  );
		    fillRes("segment","eta","#eta"       ,200,-.05,0.05  ,simMuon_gen_eta->at(iM) ,me0Muon_segment_eta ->at(goodNearSegment)  );
		    fillRes("track","phi","#phi"         ,200,-.05,0.05  ,simMuon_gen_phi->at(iM) ,me0Muon_track_phi ->at(goodNearSegment)  );
		    fillRes("segment","phi","#phi"       ,200,-.05,0.05  ,simMuon_gen_phi->at(iM) ,me0Muon_segment_phi ->at(goodNearSegment)  );
		    fillRes("track","deta","#Delta#eta"  ,200,-.2,0.2  ,simMuon_gen_deta->at(iM),me0Muon_track_deta->at(goodNearSegment) );
		    fillRes("segment","deta","#Delta#eta",200,-.2,0.2  ,simMuon_gen_deta->at(iM),me0Muon_segment_deta->at(goodNearSegment) );
		    fillRes("track","dphi","#Delta#phi"  ,200,-.006,0.006  ,simMuon_gen_dphi->at(iM),me0Muon_track_dphi->at(goodNearSegment) );
		    fillRes("segment","dphi","#Delta#phi",200,-.006,0.006  ,simMuon_gen_dphi->at(iM),me0Muon_segment_dphi->at(goodNearSegment) );

	    	plotter.getOrMake1D(TString::Format("%s_%s_real_muon_track_recomreco_phi",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #phi - reco #phi"),200,-.05,.05)->Fill(me0Muon_track_phi->at(goodNearSegment) -me0Muon_segment_phi->at(goodNearSegment));
	    	plotter.getOrMake1D(TString::Format("%s_%s_real_muon_track_recomreco_dphi",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi"),200,-.015,.015)->Fill(me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment));

	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibyPhi",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #phi - reco #phi"),50,-.015,.015,50,-.05,.05)->Fill(me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment),me0Muon_track_phi->at(goodNearSegment) -me0Muon_segment_phi->at(goodNearSegment));
	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibyDphi",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #Delta#phi"),50,-.015,.015,50,-.015,.015)->Fill(me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment),me0Muon_track_dphi->at(goodNearSegment));
	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibyeta",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #eta"),50,-.015,.015,50,1.8,3.0)->Fill(me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment),std::fabs(me0Muon_track_eta->at(goodNearSegment)));
	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibydeta",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #eta - reco #eta "),50,-.015,.015,50,-.05,0.05)->Fill(me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment),me0Muon_track_eta->at(goodNearSegment) -me0Muon_segment_eta->at(goodNearSegment));

	    	const double dPhi = (me0Muon_track_dphi->at(goodNearSegment) -me0Muon_segment_dphi->at(goodNearSegment));
	    	const double  phi = (me0Muon_track_phi->at(goodNearSegment) -me0Muon_segment_phi->at(goodNearSegment));
	    	diags[iP].add(dPhi,phi);
//	    	double rotDPhi =  dPhi*std::cos(TMath::PiOver4()) - phi*std::sin(TMath::PiOver4());
//	    	double rotPhi =  dPhi*std::sin(TMath::PiOver4()) + phi*std::cos(TMath::PiOver4());
//	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibyPhit10",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #phi - reco #phi"),50,-.05,.05,50,-.05,.05)->Fill(rotDPhi,rotPhi);
	    	double rotDPhi =  diag.transA(dPhi,phi);
	    	double rotPhi =  diag.transB(dPhi,phi);

	    	plotter.getOrMake2D(TString::Format("%s_%s_real_muon_track_recomreco_dphibyPhit10",prefix.Data(), pTSTRING.Data()),TString::Format(";gen. #Delta#phi - reco #Delta#phi ;gen. #phi - reco #phi"),500,-10,10,500,-10,10)->Fill(rotDPhi,rotPhi);
	  }

  }

  void runCalc() {
	  for(unsigned int iD = 0; iD < 6; ++iD){
		  cout << endl << iD << endl;
		  auto& diag = diags[iD];
		  diag.calcCov();
		  cout << diag.sigAA <<" "<< diag.sigAB <<endl;
		  cout << diag.sigAB <<" "<< diag.sigBB <<endl;
		  cout <<"----"<<endl;
		  diag.calcEigen();
		  cout << "("<<diag.lA <<","<<diag.lB<<")"<<endl;
		  cout << "("<<diag.vA1 <<","<<diag.vA2<<") ("<<diag.vB1 <<","<<diag.vB2<<")"<<endl;
		  cout <<"----"<<endl;
		  diag.calcTrans();
		  cout << diag.W11 <<" "<< diag.W21 <<endl;
		  cout << diag.W12 <<" "<< diag.W22 <<endl;
		  cout <<"----"<<endl;
		  cout << "("<<diag.transA(1,0) <<","<<diag.transB(1,0)<<") ("<<diag.transA(0,1) <<","<<diag.transB(0,1)<<") ("<<diag.transA(.0065,.05) <<","<<diag.transB(.0065,.05)<<")"<<endl;
		  cout << endl;

	  }
  }

  void testWorkingPoints(TString prefix) {
	  static float pBins[] = {0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50};
	  static const int nPBins = 15;

//	    static float pBins[] = {0,1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,75,100,125,150,200,300};
//	      static const int nPBins = 21;

//	  static float ptBins[] = {0,1,2,3,4,5,6,7,8,9,10,12,15,20};
//	  static const int nPtBins = 13;

	    static float ptBins[] = {0,1,2,3,5,10,15,20,25,35,40};
	      static const int nPtBins = 10;

//	    static float ptBins[] = {0,1,2,3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,60,80,100};
//	      static const int nPtBins = 22;

	  static float etaBins[] = {1.8,2.0,2.2,2.4,2.6,2.8,3.0};
	  static const int nEtaBins = 6;

	  auto passMatch = [&](int idx, double phiC, double etaC,double dPhiC)-> bool{
		  if(me0Muon_pt->at(idx) < 1) return false;
		  const double dPhi = std::fabs(me0Muon_track_dphi->at(idx) - me0Muon_segment_dphi->at(idx));
		  const double phi  = std::fabs(me0Muon_track_phi->at(idx) - me0Muon_segment_phi->at(idx));
		  const double eta  = std::fabs(me0Muon_track_eta->at(idx) - me0Muon_segment_eta->at(idx));
		  const double dPhiS2C  = std::min(std::max(alphaDPhi2SC/me0Muon_p->at(idx),maxDPhi2SC),dPhiC);
		  const double phiS2C  = std::min(std::max(alphaPhi2SC/me0Muon_p->at(idx),maxPhi2SC),phiC);
		  return eta < etaC && phi < phiC && dPhi < dPhiC && dPhi < dPhiS2C && phi < phiS2C ;
	  };

//	  auto passMed = [&](int idx) -> bool { return passMatch(idx,0.032202,0.0483636,0.00412121);};
//	  auto passLoose = [&](int idx) -> bool { return passMatch(idx,0.0241212,0.042303,0.00381818);};

	  auto passMed = [&](int idx) -> bool { return passMatch(idx,0.0564444,0.0766465,0.00957576);};
	  auto passTight = [&](int idx) -> bool { return passMatch(idx,0.032202,0.0483636,0.00412121);};

	  auto fillSigPlots = [&] (TString type, TString titleType, float pt,float mom,float absEta) {
		    plotter.getOrMake1D(TString::Format("%s_%s_pt",prefix.Data(),type.Data()),TString::Format(";%s p_{T} [GeV]; muon ID efficiency",titleType.Data()),nPtBins,ptBins)->Fill(pt);
		    plotter.getOrMake1D(TString::Format("%s_%s_p",prefix.Data(),type.Data()),TString::Format(";%s |p| [GeV]; muon ID efficiency",titleType.Data()),nPBins,pBins)->Fill(mom);
		    plotter.getOrMake2D(TString::Format("%s_%s_pt_v_eta",prefix.Data(),type.Data()),TString::Format(";%s p_{T} [GeV]; %s |#eta|",titleType.Data(),titleType.Data()),nPtBins,ptBins,nEtaBins,etaBins)->Fill(pt,absEta);
		    plotter.getOrMake2D(TString::Format("%s_%s_p_v_eta",prefix.Data(),type.Data()),TString::Format(";%s |p| [GeV]; %s |#eta|",titleType.Data(),titleType.Data()),nPBins,pBins,nEtaBins,etaBins)->Fill(mom,absEta);

	  };


	  //eff
	  int nSM = 0;
	  int nMOOAcc = 0;
	  int nMIAcc = 0;
	  for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
		  const float pt  = simMuon_pt->at(iM);
		  const float eta = simMuon_eta->at(iM);
		  const float absEta = std::fabs(eta);
		  const float mom = getP(pt,eta,simMuon_phi->at(iM));
		  if(absEta < 1.4 || absEta < 4.0) nMOOAcc ++;
		  else nMIAcc++;
		  if(absEta < 2.0 || absEta >2.8) continue;

		  int  matchIDX = -1;
		  for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
			  if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
			  if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
			  matchIDX = iRM;
			  break;

		  }

		  fillSigPlots("signal_incl","true muon",pt,mom,absEta);
		  if(matchIDX < 0) continue;
		  if(isAMatch(matchIDX,true)) fillSigPlots("signal_wp_std","true muon",pt,mom,absEta);
		  if(passMed(matchIDX)) fillSigPlots("signal_wp_95","true muon",pt,mom,absEta);
		  if(passTight(matchIDX)) fillSigPlots("signal_wp_90","true muon",pt,mom,absEta);
	  }

	  //bkg
	  if(isPureBKG || (nMOOAcc == 2 && nMIAcc == 0)){
		  plotter.getOrMake1D(TString::Format("%s_nEvtsForTestWP",prefix.Data()),";nEvtsForFakes; a.u.",1,0,2)->Fill(1);
		  for(unsigned int idx = 0; idx < me0Muon_p->size(); ++idx){
              const float pt  = me0Muon_pt->at(idx);
              const float eta = me0Muon_eta->at(idx);
              const float absEta = std::fabs(eta);
              if((me0Muon_pt->at(idx) < 2)) continue;
			  if(me0Muon_truthType->at(idx) == 0) continue;

//			  if(passTight(idx) && me0Muon_pt->at(idx) > 20){
//			      getTree()->Scan("*","","",1,eventNumber);
//			  }

			  const float mom = me0Muon_p->at(idx);
			  if(isAMatch(idx,true)) fillSigPlots("bkg_wp_std","pixel track",pt,mom,absEta);
			  if(passMed(idx)) fillSigPlots("bkg_wp_95","pixel track",pt,mom,absEta);
			  if(passTight(idx)) fillSigPlots("bkg_wp_90","pixel track",pt,mom,absEta);
		  }

	  }
  }



  virtual void runAEvent() {
	  doMatching(false);
//	  oneDCutPlots(glbPrefix);
	  effCuts(glbPrefix);
//	  makeAnnaPlot(glbPrefix);
//	  makeResPlot(glbPrefix);
//	  getNPossibleMatches(glbPrefix);
//	  getBKGComp(glbPrefix);
	  testWorkingPoints(glbPrefix);
  }

  void write(TString fileName){ plotter.write(fileName);}
  HistGetter plotter;
};

#endif

void AnalyzeME0TrackMatchingTree(std::string fileName,unsigned int nP, unsigned int nS, std::string outFileName, bool isPureBKG){
	Analyzer a(fileName,"Events");
	a.isPureBKG = isPureBKG;
	a.nP = nP;
	a.nS = nS;
//	  a.maxPhi2SC = 2*.35/float(nS);
	a.alphaPhi2SC = 1.2;
	a.maxPhi2SC = a.alphaPhi2SC/100;
	a.alphaDPhi2SC = 0.2;
	a.maxDPhi2SC = a.alphaDPhi2SC/100;

	  a.phiC = 0.05;
	  a.dPhiC = 0.0065;
	  a.etaC = nP >= 8 ? 0.06 : 0.08;

//    a.phiC = 0.032202;
//    a.dPhiC = 0.00412121;
//    a.etaC = nP >= 8 ? 0.0483636 : 0.0483636;

  a.glbPrefix = TString::Format("p%us%u",nP,nS);
  a.diag.fillTrans(10.0985,844.29 ,84.9847,-100.325);
  a.analyze();
  a.write(outFileName);
//  a.runCalc();
}
