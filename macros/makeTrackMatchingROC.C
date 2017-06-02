
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "HistoPlotting/include/Plotter.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TGraph.h"
using namespace std;



class ROCMaker {
public:
	class BinnedVar {
	public:
		unsigned int nS;
		float min;
		float max;

		BinnedVar(int nS, float min,float max) :nS(nS),min(min),max(max){}
		float get(int unsigned iS) const {return min + float(iS)*(max-min)/float(nS-1);}
		int find(float val) const {return (val - min) < 0 ? 0 : std::ceil((val - min)*float(nS-1)/(max-min)  );}// returns first int that is greater than val

	};

	class Totals {
	public:
		const unsigned int nS1;
		const unsigned int nS2;
		const unsigned int nS3;
		const unsigned int nTot;
		float * values;
		Totals(unsigned int nS1, unsigned int nS2, unsigned int nS3) : nS1(nS1),nS2(nS2),nS3(nS3),nTot(nS1*nS2*nS3){
			values = new float[nTot];
			for(unsigned int iF = 0; iF < nTot; ++iF) values[iF] = 0;
		}
		~Totals(){delete values;}
		float& at(unsigned int iS1, unsigned int iS2, unsigned int iS3){ return values[iS1 + nS1*(iS2 + nS2*iS3)];}
	};



	ROCMaker() :
		phiCut(100,.012,.212),
		etaCut(100,.012,.212),
		dPhiCut(100,.002,.032),
		signalValues(100,100,100),
		bkgValues(100,100,100)
	{

	}

	void add(float phi, float eta, float dPhi, Totals* vals ){
		const int minIDXPhi = phiCut.find(phi);
		const int minIDXEta = etaCut.find(eta);
		const int minIDXDPhi = dPhiCut.find(dPhi);
//		cout << phi <<" "<<eta<<" "<< dPhi<<" "<< minIDXPhi<<" "<<minIDXEta<<" "<<minIDXDPhi<<endl;
		if(minIDXPhi >= phiCut.nS || minIDXEta >= etaCut.nS || minIDXDPhi >= dPhiCut.nS) return;

		for(unsigned int iS1 = minIDXPhi; iS1 < phiCut.nS; ++iS1 )
			for(unsigned int iS2 = minIDXEta; iS2 < etaCut.nS; ++iS2 )
				for(unsigned int iS3 = minIDXDPhi; iS3 < dPhiCut.nS; ++iS3 ){
//					cout << vals->at(iS1,iS2,iS3) ;
					vals->at(iS1,iS2,iS3) += 1.0;
//					cout << " "<<vals->at(iS1,iS2,iS3) <<endl;
				}

	}
	void addSignal(float phi, float eta, float dPhi) { nTotSig++; add(phi,eta,dPhi,&signalValues);}
	void addEmptySignal() { nTotSig++;}
	void addBKG(float phi, float eta, float dPhi) { add(phi,eta,dPhi,&bkgValues);}
	void addBKGEvent(){nBKGEvents++;};

	void compute() {
		int nP = 0;
		TGraph * graph = new TGraph();
//		vector<float> effs = {0.99,0.985,.98,.979,.978,.977,.976,.975,.97,.965,.96,.955,.95,.945,.94,.93,.92,.91,.90,.89,.88,.87,.86,.85,.80,.75,.7,.6,.5,.4,.3};
		vector<float> effs = {0.99,0.985,.98,.979,.978,.977,.976,.975,.97,.965,.96,.955,.95,.945,.94,.93,.92,.91,.90,.89,.88,.87,.86,.85,.84,.83,.82,.81,.80,.79,.78,.77,.76,.75};
		for(unsigned int iE = 0; iE < effs.size(); ++iE){
			const float minSignal = effs[iE]*nTotSig;
			float minBKG =-1;
			int minBKGS1 =-1;int minBKGS2 =-1;int minBKGS3 =-1;

			for(unsigned int iS1 = 0; iS1 < phiCut.nS; ++iS1)
				for(unsigned int iS2 = 0; iS2 < etaCut.nS; ++iS2)
					for(unsigned int iS3 = 0; iS3 < dPhiCut.nS; ++iS3){
//						cout << "("<<iS1 <<","<<iS2 <<","<<iS3<<","<<signalValues.at(iS1,iS2,iS3)<<") ";
						if(signalValues.at(iS1,iS2,iS3) < minSignal ) continue;
						if(minBKG < 0 || bkgValues.at(iS1,iS2,iS3) < minBKG){
							minBKG =bkgValues.at(iS1,iS2,iS3);
							minBKGS1 =iS1; minBKGS2 = iS2; minBKGS3 = iS3;
						}
						break; //no need to cut looser
					}

			if(minBKG >= 0){
				cout << effs[iE] <<"\t"<< signalValues.at(minBKGS1,minBKGS2,minBKGS3)/nTotSig <<"\t"<<minBKG/nBKGEvents<<"\t"<< phiCut.get(minBKGS1)<<"\t"<< etaCut.get(minBKGS2)<<"\t"<< dPhiCut.get(minBKGS3)<<endl;
				graph->SetPoint(nP,minBKG/nBKGEvents,signalValues.at(minBKGS1,minBKGS2,minBKGS3)/nTotSig);
				nP++;
			}



		}

		graph->SetLineColor  (1);
		graph->SetLineWidth  (4);
		graph->SetLineStyle  (1);
		graph->SetMarkerStyle(20);
		graph->SetMarkerColor(1);
		graph->SetMarkerSize (1);
		Drawing::Drawable1D drawableF("A P L","ROC",Drawing::GRAPH,graph,false);
		drawableF.graphAxisHist = new TH1F("axisHist",";N. of p_{T} > 3 GeV background matches per BX;p_{T} 3-5 GeV signal muon efficiency",100,0,1);
		Plotter  * p = new Plotter;
		p->addDrawable(drawableF);

		p->draw(false,"ROC");
		p->setXTitle("N. of p_{T} > 3 GeV background matches per BX");
		p->setYTitle("p_{T} 3-5 GeV signal muon efficiency");



	}

	BinnedVar phiCut;
	BinnedVar etaCut;
	BinnedVar dPhiCut;
	Totals signalValues;
	Totals bkgValues;
	float nTotSig= 0;
	float nBKGEvents=0;


};



double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
double deltaR2(const float  eta1, float phi1,  const float eta2, float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2)*TVector2::Phi_mpi_pi(phi1 - phi2)  + (eta1-eta2)*(eta1-eta2);}

class TrackMatchingAnalyzer : public BaseTupleAnalyzer{
public:
	TString glbPrefix = "";

	TrackMatchingAnalyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){
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

	virtual void runAEvent() {
	}
};


class SignalAnalyzer : public TrackMatchingAnalyzer {
public:
	ROCMaker* rocMaker;
	SignalAnalyzer(ROCMaker* rocMaker, std::string fileName, std::string treeName) :rocMaker(rocMaker),TrackMatchingAnalyzer(fileName,treeName){}
	virtual void runAEvent() {
		for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
			const float pt = simMuon_pt->at(iM);
			if(pt < 2) continue;
			if(pt > 4) continue;

			bool goodTrack   = simMuon_trackIDX->at(iM) >= 0;
			int  goodNearSegment = -1;
			for(unsigned int iRM = 0; iRM < me0Muon_pt->size(); ++iRM){
				if(me0Muon_segIDX->at(iRM) != simMuon_segmentIDX->at(iM)) continue;
				if(me0Muon_trackIDX->at(iRM) != simMuon_trackIDX->at(iM)) continue;
				goodNearSegment = iRM;
				break;
			}

			const double dPhi = goodNearSegment >= 0 ? std::fabs(me0Muon_track_dphi->at(goodNearSegment) - me0Muon_segment_dphi->at(goodNearSegment)) : -1;
			const double phi  = goodNearSegment >= 0 ? std::fabs(me0Muon_track_phi->at(goodNearSegment) - me0Muon_segment_phi->at(goodNearSegment)) : -1;
			const double eta  = goodNearSegment >= 0 ? std::fabs(me0Muon_track_eta->at(goodNearSegment) - me0Muon_segment_eta->at(goodNearSegment)) : -1;

			bool passPDepCuts = goodNearSegment >= 0 && phi < std::max(1.2/me0Muon_p->at(goodNearSegment),1.2/100)
			&& dPhi < std::max(0.2/me0Muon_p->at(goodNearSegment),0.2/100) ;
//			cout << goodTrack <<" "<< goodNearSegment<<" "<<passPDepCuts<<endl;
			if(passPDepCuts )
				rocMaker->addSignal(phi,eta,dPhi);
			else
				rocMaker->addEmptySignal();
		}

	}
};

class BKGAnalyzer :public TrackMatchingAnalyzer {
public:
	ROCMaker* rocMaker;
	BKGAnalyzer(ROCMaker* rocMaker, std::string fileName, std::string treeName) :rocMaker(rocMaker),TrackMatchingAnalyzer(fileName,treeName){}
	virtual void runAEvent() {
		rocMaker->addBKGEvent();

		for(unsigned int iRM = 0; iRM < me0Muon_truthType->size(); ++iRM){
			if(me0Muon_truthType->at(iRM) == 0) continue;
			const double dPhi = std::fabs(me0Muon_track_dphi->at(iRM) - me0Muon_segment_dphi->at(iRM));
			const double phi  = std::fabs(me0Muon_track_phi->at(iRM) - me0Muon_segment_phi->at(iRM));
			const double eta  = std::fabs(me0Muon_track_eta->at(iRM) - me0Muon_segment_eta->at(iRM));
			if(me0Muon_pt->at(iRM) < 2) continue;
			bool passPDepCuts = phi < std::max(1.2/me0Muon_p->at(iRM),1.2/100)
			&& dPhi < std::max(0.2/me0Muon_p->at(iRM),0.2/100) ;
			if(passPDepCuts){
				rocMaker->addBKG(phi,eta,dPhi);
			}
		}
	}
};



#endif

void makeTrackMatchingROC(){
	ROCMaker maker;
	cout << "Load Signal!"<<endl;
	SignalAnalyzer sigAnalyzer(&maker,"/Users/nmccoll/Dropbox/Work/Projects/ME0/3_13_17_trackMuonMatching/TDRVersion/trackMatchingTreeForPOGHP_p8s384.root","Events");
	sigAnalyzer.analyze();
	cout << "Load BKG!"<<endl;
	BKGAnalyzer bkgAnalyzer(&maker,"/Users/nmccoll/Dropbox/Work/Projects/ME0/3_13_17_trackMuonMatching/TDRVersion/trackMatchingTreeForPOGHP_NU_p8s384.root","Events");
	bkgAnalyzer.analyze();
	cout <<"Compute!"<<endl;
	maker.compute();
}
