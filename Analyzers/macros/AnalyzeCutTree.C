
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "HistoPlotting/include/Plotter.h"
#include "BaseCutTree.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TGraph.h"
#include "TFitter.h"

class HiggsSolver {
 public:
  TFitter *minimizer;

  HiggsSolver();
  ~HiggsSolver();


  // operations
  static double higgsSolverFunction(double nz_,
			       double lx_, double ly_, double lz_,  double le_,
			       double jx_, double jy_, double jz_, double je_,
			       double metx_, double mety_
			       );

  static void minuitFunctionWrapper(int& nDim, double* gout, double& result, double *par, int flg);
  double massMinimization   (const CylLorentzVectorF& lep_, const CylLorentzVectorF& jet_,  const CylLorentzVectorF& met_, CylLorentzVectorF& higgsOutput);


};


HiggsSolver::HiggsSolver() : minimizer(new TFitter(4))
  {
  // less print-outs
    cout.precision(11);
    double p1 = -1;
    minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);

    // tell minimizer about the function to be minimized
    minimizer->SetFCN(minuitFunctionWrapper);
  }

HiggsSolver::~HiggsSolver(){
  delete minimizer;
}

double HiggsSolver::higgsSolverFunction(double nz_,
	       double lx_, double ly_, double lz_,  double le_,
	       double jx_, double jy_, double jz_, double je_,
	       double metx_, double mety_
			       ) {

  cout.precision(11);

  const double nE = std::sqrt(metx_*metx_+mety_*mety_+nz_*nz_);
  const double px = jx_+lx_+metx_;
  const double py = jy_+ly_+mety_;
  const double pz = jz_+lz_+nz_;
  const double W1_mass2 =  (nE + le_)*(nE + le_) - (px+lx_)*(px+lx_) - (py+ly_)*(py+ly_) - (pz+lz_)*(pz+lz_);
  const double W2_mass2 =  je_*je_ - jx_*jx_ -jy_*jy_ -jz_*jz_;
  const double minWMass = std::sqrt(std::min(W1_mass2,W2_mass2));
  const double maxWMass = std::sqrt(std::max(W1_mass2,W2_mass2));

//  const double metRE = 0.20;
//  const double jetRE = 0.05;
//  const double lepRE = 0.05;
//
//  const double Emsq = je_*je_*jetRE*jetRE + le_*le_*lepRE*lepRE + metRE*(metx_*metx_+mety_*mety_);

  const double massDiff = 125 - std::sqrt(  (nE + le_ + je_)*(nE + le_ + je_) - px*px - py*py - pz*pz      );
  return massDiff*massDiff + (maxWMass-80)*(maxWMass-80);// + (minWMass-30)*(minWMass-30)/25;

} // ~ end of Topness function



void HiggsSolver::minuitFunctionWrapper(int& nDim, double* gout, double& result, double* par, int flg) {


  result = higgsSolverFunction(par[0],par[1],par[2],par[3],
			   par[4],par[5],par[6],par[7],
			   par[8],par[9],par[10]
			   );

} // ~end of minuit function



double HiggsSolver::massMinimization(const CylLorentzVectorF& lep_, const CylLorentzVectorF& jet_,  const CylLorentzVectorF& met_, CylLorentzVectorF& higgsOutput) {
  cout.precision(11);


  double lx_ = lep_.px();
  double ly_ = lep_.py();
  double lz_ = lep_.pz();
  double le_ = lep_.energy();

  double jx_ = jet_.px();
  double jy_ = jet_.py();
  double jz_ = jet_.pz();
  double je_ = jet_.energy();

  double metx_  = met_.px();
  double mety_  = met_.py();


  // Define parameters [param number, param name, init val, estimated distance to min, bla, bla] // 300,3000,-3000,3000
  minimizer->SetParameter(0,"nz_",0,500,-3000,3000);

  // fixed parameters
  minimizer->SetParameter(1,"lx_",lx_,0,lx_-0.001,lx_+0.001);
  minimizer->SetParameter(2,"ly_",ly_,0,ly_-0.001,ly_+0.001);
  minimizer->SetParameter(3,"lz_",lz_,0,lz_-0.001,lz_+0.001);
  minimizer->SetParameter(4,"le_",le_,0,le_-0.001,le_+0.001);
  minimizer->SetParameter(5,"jx_",jx_,0,jx_-0.001,jx_+0.001);
  minimizer->SetParameter(6,"jy_",jy_,0,jy_-0.001,jy_+0.001);
  minimizer->SetParameter(7,"jz_",jz_,0,jz_-0.001,jz_+0.001);
  minimizer->SetParameter(8,"je_",je_,0,je_-0.001,je_+0.001);
  minimizer->SetParameter(9,"metx_",metx_,0,metx_-0.001,metx_+0.001);
  minimizer->SetParameter(10,"mety_",mety_,0,mety_-0.001,mety_+0.001);

  minimizer->FixParameter(1);
  minimizer->FixParameter(2);
  minimizer->FixParameter(3);
  minimizer->FixParameter(4);
  minimizer->FixParameter(5);
  minimizer->FixParameter(6);
  minimizer->FixParameter(7);
  minimizer->FixParameter(8);
  minimizer->FixParameter(9);
  minimizer->FixParameter(10);

  // Run the simplex minimizer to get close to the minimum [no good precision]
  minimizer->ExecuteCommand("SIMPLEX",0,0);

  // Run the migrad minimizer to precisely estimate the minimum
  //    minimizer->ExecuteCommand("MIGRAD",0,0);


  //Get the best fit values
  double nz_fit = minimizer->GetParameter(0);


  // get the function value at best fit

  higgsOutput = ASTypes::CartLorentzVector(lx_+jx_+metx_,ly_+jy_+mety_,lz_+jz_+nz_fit, std::sqrt(metx_*metx_+mety_*mety_+nz_fit*nz_fit) + le_ + je_);
  return  higgsSolverFunction(nz_fit,lx_,ly_,lz_,
		  le_,jx_,jy_,jz_,
		  je_,metx_,mety_);


} // ~ end of Topness Minimization()

class Analyzer : public BaseCutAnalyzer {
public:
    HiggsSolver solver;


	Analyzer(std::string fileName, std::string treeName) : BaseCutAnalyzer(fileName,treeName){}
	void makeLeptonPlots(TString prefix){
		plotter.getOrMake1D(TString::Format("%s_ht",prefix.Data()),";H_{T} [GeV]; arbitrary units",200,0,2000)->Fill(ht,nWeight);
		plotter.getOrMake1D(TString::Format("%s_lepton_pt",prefix.Data()),";lepton p_{T} [GeV]; arbitrary units",40,0,200)->Fill(lepton_pt->size() ? lepton_pt->at(0) : float(0.0),nWeight);
		plotter.getOrMake1D(TString::Format("%s_nLeptons",prefix.Data()),";N. of leptons; arbitrary units",5,-0.5,4.5)->Fill(lepton_pt->size(),nWeight);
		if(lepton_pt->size() && lepton_pt->at(0) >= 20 )plotter.getOrMake1D(TString::Format("%s_goodLep_nLeptons",prefix.Data()),";N. of leptons; arbitrary units",5,-0.5,4.5)->Fill(lepton_pt->size(),nWeight);
	}

	void makeHbbPlots(TString prefix){
		if(process == 10){
			CylLorentzVectorF genLEPMET =CylLorentzVectorF(w1d1_pt,w1d1_eta,w1d1_phi,w1d1_mass)+CylLorentzVectorF(w1d2_pt,w1d2_eta,w1d2_phi,w1d2_mass);
			CylLorentzVectorF genHBB(hbb_pt,hbb_eta,hbb_phi,hbb_mass);
			float genDPHI =  std::fabs(deltaPhi(genLEPMET.phi(),genHBB.phi()));
			plotter.getOrMake1D(TString::Format("%s_gen_dPHI_lepWHbb",prefix.Data()),";gen #Delta#phi(lep W, Hbb); arbitrary units",320,0,3.2)->Fill(genDPHI,nWeight);
		}

		plotter.getOrMake1D(TString::Format("%s_nFJ",prefix.Data()),";N. of fat jets; arbitrary units",10,-0.5,9.5)->Fill(fj_pt->size(),nWeight);
		if(!fj_pt->size()) return;

		CylLorentzVectorF lepton = getLep();
		CylLorentzVectorF metV = getMET();
		CylLorentzVectorF metPLep = lepton+metV;
		float maxDPHI = 0;
		for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
			CylLorentzVectorF jet = getFJ(iJ);
			float dPhi =  std::fabs(deltaPhi(metPLep.phi(),jet.phi()));
			if(dPhi > maxDPHI) maxDPHI = dPhi;
		}

		plotter.getOrMake1D(TString::Format("%s_maxFHDPhi",prefix.Data()),";maximum |#Delta#phi(reco W, fj)|; arbitrary units",320,0,3.2)->Fill(maxDPHI,nWeight);


		int nHighDPHI = 0;
		int sigClosestIDX = -1;
		int nPass = 0;

		if(isSignal()){
			CylLorentzVectorF hbb(hbb_pt,hbb_eta,hbb_phi,hbb_mass);
			for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
			CylLorentzVectorF jet = getFJ(iJ);
			float dPhi =  std::fabs(deltaPhi(metPLep.phi(),jet.phi()));
			if(dPhi < TMath::PiOver2()) continue;
			if(sigClosestIDX < 0 || deltaR2(hbb,jet) < deltaR2(hbb,getFJ(sigClosestIDX))  ) sigClosestIDX = iJ;
		}}


		for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){

			CylLorentzVectorF jet = getFJ(iJ);
			float dPhi =  std::fabs(deltaPhi(metPLep.phi(),jet.phi()));
			if(dPhi < TMath::PiOver2()) continue;

			CylLorentzVectorF sdJet = getSDFJ(iJ);
			bool passMass = (sdJet.mass() > 100 && sdJet.mass() < 150);
			bool passTau = (fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ)) < 0.7;
			bool passCSV = (std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)) > 0.6);

			nHighDPHI++;
			if(passMass && passTau && passCSV) nPass++;

			if(isSignal() && sigClosestIDX != iJ ) continue;

			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_pt",prefix.Data()),";fj p_{T} [GeV]; arbitrary units",200,0,2000)->Fill(jet.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_mass",prefix.Data()),";fj mass [GeV]; arbitrary units",250,0,500)->Fill(jet.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_sd_pt",prefix.Data()),";fj softdrop p_{T} [GeV]; arbitrary units",200,0,2000)->Fill(sdJet.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_sd_mass",prefix.Data()),";fj  softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(sdJet.mass(),nWeight);

			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_tau2otau1",prefix.Data()),";fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_csv",prefix.Data()),";fj csv; arbitrary units",100,0,1)->Fill(fj_csv->at(iJ),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_minSJPT",prefix.Data()),";fj min subjet p_{T} [GeV]; arbitrary units",100,0,500)->Fill(std::min(fj_sj1_pt->at(iJ),fj_sj2_pt->at(iJ)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_minSJCSV",prefix.Data()),";fj min subjet csv; arbitrary units",100,0,1)->Fill(std::min(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_maxSJCSV",prefix.Data()),";fj max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);

			if(passCSV && passTau){
				plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_oneM_sd_mass",prefix.Data()),";fj  softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(sdJet.mass(),nWeight);
				plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_oneM_mass",prefix.Data()),";fj mass [GeV]; arbitrary units",250,0,500)->Fill(jet.mass(),nWeight);
			}
			if(passCSV && passMass) plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_oneM_tau2otau1",prefix.Data()),";fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ),nWeight);
			if(passTau && passMass){
				plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_oneM_maxSJCSV",prefix.Data()),";fj max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);
				if(passCSV) plotter.getOrMake1D(TString::Format("%s_highDPHI_fj_oneM_minSJCSV",prefix.Data()),";fj min subjet csv; arbitrary units",100,0,1)->Fill(std::min(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);
			}
		}

		plotter.getOrMake1D(TString::Format("%s_nHighDPhi",prefix.Data()),";N. of fat jets (away); arbitrary units",10,-0.5,9.5)->Fill(nHighDPHI,nWeight);
		plotter.getOrMake1D(TString::Format("%s_nHighDPhi_pass",prefix.Data()),";N. of fat jets pass hbb selection; arbitrary units",10,-0.5,9.5)->Fill(nPass,nWeight);
	}

	void makeWjjPlots(TString prefix){
		if(isSignal()){
			CylLorentzVectorF genLep(w1d1_pt,w1d1_eta,w1d1_phi,w1d1_mass);
			CylLorentzVectorF genLEPMET =genLep+CylLorentzVectorF(w1d2_pt,w1d2_eta,w1d2_phi,w1d2_mass);
			CylLorentzVectorF genWjj =CylLorentzVectorF(w2d1_pt,w2d1_eta,w2d1_phi,w2d1_mass)+CylLorentzVectorF(w2d2_pt,w2d2_eta,w2d2_phi,w2d2_mass);

			plotter.getOrMake1D(TString::Format("%s_gen_dPHI_lepWWjj",prefix.Data()),";gen #Delta#phi(lep W, Wjj); arbitrary units",320,0,3.2)->Fill(std::fabs(deltaPhi(genLEPMET.phi(),genWjj.phi())),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_dR_lepWWjj",prefix.Data()),";gen #Delta#R(lep W, Wjj); arbitrary units",320,0,3.2)->Fill(std::sqrt(deltaR2(genLEPMET,genWjj)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_dPHI_lepWjj",prefix.Data()),";gen #Delta#phi(lep, Wjj); arbitrary units",320,0,3.2)->Fill(std::fabs(deltaPhi(genLep.phi(),genWjj.phi())),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_dR_lepWjj",prefix.Data()),";gen #Delta#R(lep, Wjj); arbitrary units",320,0,3.2)->Fill(std::sqrt(deltaR2(genLep,genWjj)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_WWjj_pt",prefix.Data()),";gen Wjj p_{T}; arbitrary units",200,0,2000)->Fill(genWjj.pt(),nWeight);

			if(genWjj.pt() > 200){
				plotter.getOrMake1D(TString::Format("%s_gen_Wjj_ptgeq200_dPHI_lepWWjj",prefix.Data()),";gen #Delta#phi(lep W, Wjj); arbitrary units",320,0,3.2)->Fill(std::fabs(deltaPhi(genLEPMET.phi(),genWjj.phi())),nWeight);
				plotter.getOrMake1D(TString::Format("%s_gen_Wjj_ptgeq200_dR_lepWWjj",prefix.Data()),";gen #Delta#R(lep W, Wjj); arbitrary units",320,0,3.2)->Fill(std::sqrt(deltaR2(genLEPMET,genWjj)),nWeight);
				plotter.getOrMake1D(TString::Format("%s_gen_Wjj_ptgeq200_dPHI_lepWjj",prefix.Data()),";gen #Delta#phi(lep, Wjj); arbitrary units",320,0,3.2)->Fill(std::fabs(deltaPhi(genLep.phi(),genWjj.phi())),nWeight);
				plotter.getOrMake1D(TString::Format("%s_gen_Wjj_ptgeq200_dR_lepWjj",prefix.Data()),";gen #Delta#R(lep, Wjj); arbitrary units",320,0,3.2)->Fill(std::sqrt(deltaR2(genLep,genWjj)),nWeight);

			}

		}

		CylLorentzVectorF lepton = getLep();
		CylLorentzVectorF metV = getMET();
		CylLorentzVectorF metPLep = lepton+metV;

		int wJJIDX = -1;
		double minDR2 = -1;
		for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
			auto jet = getFJ(iJ);
			const double absDPhi = absDeltaPhi(jet,metPLep);
			if(absDPhi > 0.6) continue;
			const double dr2 = deltaR2(jet,lepton);
			if(wJJIDX < 0 || dr2 < minDR2){
				minDR2 =dr2;
				wJJIDX = iJ;
			}
		}



		auto makeJetPlot = [&](const int iJ,const CylLorentzVectorF& jet, const CylLorentzVectorF& sdJet, TString jetName ){
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_pt"       ,prefix.Data(),jetName.Data()),";fj p_{T} [GeV]; arbitrary units",200,0,2000)->Fill(jet.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_mass"     ,prefix.Data(),jetName.Data()),";fj mass [GeV]; arbitrary units",250,0,500)->Fill(jet.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_sd_pt"    ,prefix.Data(),jetName.Data()),";fj softdrop p_{T} [GeV]; arbitrary units",200,0,2000)->Fill(sdJet.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_sd_mass"  ,prefix.Data(),jetName.Data()),";fj  softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(sdJet.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_tau2otau1",prefix.Data(),jetName.Data()),";fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_csv"      ,prefix.Data(),jetName.Data()),";fj csv; arbitrary units",100,0,1)->Fill(fj_csv->at(iJ),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj%sfj_maxSJCSV" ,prefix.Data(),jetName.Data()),";fj max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);

			bool passMass = sdJet.mass() > 10 && sdJet.mass() < 120;
			if(passMass){
				plotter.getOrMake1D(TString::Format("%s_wjj%sfj_oneM_csv"      ,prefix.Data(),jetName.Data()),";fj csv; arbitrary units",100,0,1)->Fill(fj_csv->at(iJ),nWeight);
				plotter.getOrMake1D(TString::Format("%s_wjj%sfj_oneM_maxSJCSV" ,prefix.Data(),jetName.Data()),";fj max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)),nWeight);
				plotter.getOrMake1D(TString::Format("%s_wjj%sfj_oneM_tau2otau1",prefix.Data(),jetName.Data()),";fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ),nWeight);
			}

		};


		if(wJJIDX < 0){
			plotter.getOrMake1D(TString::Format("%s_wJJ_cand",prefix.Data()),";fat jet near lepton; arbitrary units",3,-0.5,2.5)->Fill(0.0,nWeight);
			return;
		}

		auto jet = getFJ(wJJIDX);
		CylLorentzVectorF sdJet = getSDFJ(wJJIDX);

		makeJetPlot(wJJIDX,jet,sdJet,"_");

		CylLorentzVectorF sdJetLepSub = sdJet;
		CylLorentzVectorF jetLepSub = jet;
		const double deltaRLep = deltaR(jet,lepton);
		if(deltaRLep < 0.8){
			jetLepSub = lepton.pt() > jetLepSub.pt() ? CylLorentzVectorF() : jetLepSub - lepton;
			sdJetLepSub = lepton.pt() > sdJetLepSub.pt() ? CylLorentzVectorF() : sdJetLepSub - lepton;
		}
		makeJetPlot(wJJIDX,jetLepSub,sdJetLepSub,"_lepSub_");

		if(sdJetLepSub.pt() > 0 && sdJetLepSub.mass() > 10 && sdJetLepSub.mass() < 120 && std::max(fj_sj1_csv->at(wJJIDX),fj_sj2_csv->at(wJJIDX)) < 0.9) plotter.getOrMake1D(TString::Format("%s_wJJ_cand",prefix.Data()),";fat jet near lepton; arbitrary units",3,-0.5,2.5)->Fill(2.0,nWeight);
		else plotter.getOrMake1D(TString::Format("%s_wJJ_cand",prefix.Data()),";fat jet near lepton; arbitrary units",3,-0.5,2.5)->Fill(1.0,nWeight);


	}


	void makeAK4WjjPlots(TString prefix){
		if(isSignal()){
			CylLorentzVectorF genLep(w1d1_pt,w1d1_eta,w1d1_phi,w1d1_mass);
			CylLorentzVectorF genLEPMET =genLep+CylLorentzVectorF(w1d2_pt,w1d2_eta,w1d2_phi,w1d2_mass);
			CylLorentzVectorF genWj1(w2d1_pt,w2d1_eta,w2d1_phi,w2d1_mass);
			CylLorentzVectorF genWj2(w2d2_pt,w2d2_eta,w2d2_phi,w2d2_mass);

			plotter.getOrMake1D(TString::Format("%s_gen_min_Wj",prefix.Data()),";gen min Wj p_{T} [GeV]; arbitrary units",200,0,1000)->Fill(std::min(w2d1_pt,w2d2_pt),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_max_Wj",prefix.Data()),";gen max Wj p_{T} [GeV]; arbitrary units",200,0,1000)->Fill(std::max(w2d1_pt,w2d2_pt),nWeight);

			plotter.getOrMake1D(TString::Format("%s_gen_dR_Wj1j2"       ,prefix.Data()),";gen #DeltaR(Wj1, Wj2); arbitrary units",320,0,3.2)->Fill(std::sqrt(deltaR2(genWj1,genWj2)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_min_dR_lepWWj"  ,prefix.Data()),";gen min #DeltaR(lepW, Wj); arbitrary units",320,0,3.2)->Fill(std::min(std::sqrt(deltaR2(genWj2,genLEPMET)),std::sqrt(deltaR2(genWj1,genLEPMET))),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_max_dR_lepWWj"  ,prefix.Data()),";gen min #DeltaR(lepW, Wj); arbitrary units",320,0,3.2)->Fill(std::max(std::sqrt(deltaR2(genWj2,genLEPMET)),std::sqrt(deltaR2(genWj1,genLEPMET))),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_min_dR_lepWj"   ,prefix.Data()),";gen min #DeltaR(lep, Wj); arbitrary units",320,0,3.2)->Fill(std::min(std::sqrt(deltaR2(genWj2,genLep)),std::sqrt(deltaR2(genWj1,genLep))),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_max_dR_lepWj"   ,prefix.Data()),";gen min #DeltaR(lep, Wj); arbitrary units",320,0,3.2)->Fill(std::max(std::sqrt(deltaR2(genWj2,genLep)),std::sqrt(deltaR2(genWj1,genLep))),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_min_dPhi_lepWWj",prefix.Data()),";gen min #Delta#phi(lepW, Wj); arbitrary units",320,0,3.2)->Fill(std::min(absDeltaPhi(genWj2,genLEPMET),absDeltaPhi(genWj1,genLEPMET)),nWeight);
			plotter.getOrMake1D(TString::Format("%s_gen_max_dPhi_lepWWj",prefix.Data()),";gen min #Delta#phi(lepW, Wj); arbitrary units",320,0,3.2)->Fill(std::max(absDeltaPhi(genWj2,genLEPMET),absDeltaPhi(genWj1,genLEPMET)),nWeight);

		}
		auto lepton = getLep();
		CylLorentzVectorF metV = getMET();
		CylLorentzVectorF metPLep = lepton+metV;
		vector<CylLorentzVectorF> jets;
		for(unsigned int iJ = 0; iJ < jet_pt->size(); ++iJ) jets.push_back(getJet(iJ));
		double distToLep = 0;
		int jetOverlapIDX = findNearestDR(lepton,jets,distToLep,0.4);
		if(jetOverlapIDX >= 0) {
			CylLorentzVectorF newJetPT = jets[jetOverlapIDX];
			if(lepton.pt() > newJetPT.pt() || (newJetPT - lepton).pt() < 20 ) newJetPT = CylLorentzVectorF();
			else newJetPT = newJetPT - lepton;
			jets[jetOverlapIDX] =newJetPT;
		}
		vector<std::pair<CylLorentzVectorF,float> > jetsNearW;
		for(unsigned int iJ= 0; iJ < jets.size(); ++iJ){
			if(jets[iJ].pt() < 20 ) continue;
			if(absDeltaPhi(jets[iJ],metPLep) > 0.8) continue;
			jetsNearW.emplace_back(jets[iJ],jet_csv->at(iJ) );
		}
		int nPairs = 0;
		int bestPairIJ1 = -1;
		int bestPairIJ2 = -1;
		for(unsigned int iJ1= 0; iJ1 < jetsNearW.size(); ++iJ1){
			for(unsigned int iJ2= iJ1+1; iJ2 < jetsNearW.size(); ++iJ2){
				if(deltaR2(jetsNearW[iJ1].first,jetsNearW[iJ2].first) > 1.0*1.0) continue;
				nPairs++;
				if(bestPairIJ1 || deltaR2(jetsNearW[iJ1].first+jetsNearW[iJ2].first,lepton) < deltaR2(jetsNearW[bestPairIJ1].first+jetsNearW[bestPairIJ2].first,lepton)  ){
					bestPairIJ1 = iJ1; bestPairIJ2 = iJ2;
				}
		}
		}
		plotter.getOrMake1D(TString::Format("%s_nak4_nearW",prefix.Data()),";N. ak4 nearW; arbitrary units",20,-0.5,19.5)->Fill(jetsNearW.size(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_nak4_pairs_nearW",prefix.Data()),";N. ak4 pairs nearW; arbitrary units",20,-0.5,19.5)->Fill(nPairs,nWeight);

		if(!jetsNearW.size())plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(0.0,nWeight);
		else if(jetsNearW.size() == 1){
			if(jetsNearW[0].first.mass() > 30 && jetsNearW[0].second < CSV_MEDIUM) plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(2.0,nWeight);
			else plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(1.0,nWeight);
		} else{
			if(!nPairs)plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(3.0,nWeight);
			else if(jetsNearW[bestPairIJ1].second < CSV_MEDIUM && jetsNearW[bestPairIJ2].second < CSV_MEDIUM) plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(5.0,nWeight);
			else plotter.getOrMake1D(TString::Format("%s_ak4_categories",prefix.Data()),";None/1j fail/1j pass/+1j fail/ 2j fail/2j pass; arbitrary units",6,-0.5,5.5)->Fill(4.0,nWeight);
		}





		if(!jetsNearW.size()) return;

		if(jetsNearW.size() == 1){
			plotter.getOrMake1D(TString::Format("%s_ak4_1j_mass"     ,prefix.Data()),";1j ak4 mass [GeV]; arbitrary units",250,0,500)->Fill(jetsNearW[0].first.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_ak4_1j_csv"     ,prefix.Data()),";1j ak4 csv; arbitrary units",100,0,1)->Fill(jetsNearW[0].second,nWeight);
		}
		if(!nPairs) return;

		for(unsigned int iJ1= 0; iJ1 < jetsNearW.size(); ++iJ1){
			for(unsigned int iJ2= iJ1+1; iJ2 < jetsNearW.size(); ++iJ2){
				if(deltaR2(jetsNearW[iJ1].first,jetsNearW[iJ2].first) > 1.0*1.0) continue;

				plotter.getOrMake1D(TString::Format("%s_ak4_2j_mass"     ,prefix.Data()),";2j ak4 mass [GeV]; arbitrary units",250,0,500)->Fill((jetsNearW[iJ1].first+jetsNearW[iJ2].first).mass(),nWeight);
				plotter.getOrMake1D(TString::Format("%s_ak4_2j_mincsv"     ,prefix.Data()),";2j ak4 min csv; arbitrary units",100,0,1)->Fill(std::min(jetsNearW[iJ1].second,jetsNearW[iJ2].second),nWeight);
				plotter.getOrMake1D(TString::Format("%s_ak4_2j_maxcsv"     ,prefix.Data()),";2j ak4 max csv; arbitrary units",100,0,1)->Fill(std::max(jetsNearW[iJ1].second,jetsNearW[iJ2].second),nWeight);
				if(iJ1 == bestPairIJ1 && iJ2 == bestPairIJ2 ){
					plotter.getOrMake1D(TString::Format("%s_best_ak4_2j_mass"     ,prefix.Data()),";2j ak4 mass [GeV]; arbitrary units",250,0,500)->Fill((jetsNearW[iJ1].first+jetsNearW[iJ2].first).mass(),nWeight);
					plotter.getOrMake1D(TString::Format("%s_best_ak4_2j_mincsv"     ,prefix.Data()),";2j ak4 min csv; arbitrary units",100,0,1)->Fill(std::min(jetsNearW[iJ1].second,jetsNearW[iJ2].second),nWeight);
					plotter.getOrMake1D(TString::Format("%s_best_ak4_2j_maxcsv"     ,prefix.Data()),";2j ak4 max csv; arbitrary units",100,0,1)->Fill(std::max(jetsNearW[iJ1].second,jetsNearW[iJ2].second),nWeight);
				}

			}
		}




	}

	void makeAK4METPlots(TString prefix, const CylLorentzVectorF& met, const CylLorentzVectorF& lepton, const unsigned int ihbb, const CylLorentzVectorF& hbb){
		CylLorentzVectorF metPLep = lepton+met;
		vector<CylLorentzVectorF> jets;
		for(unsigned int iJ = 0; iJ < jet_pt->size(); ++iJ) jets.push_back(getJet(iJ));
		double distToLep = 0;
		int jetOverlapIDX = findNearestDR(lepton,jets,distToLep,0.4);
		if(jetOverlapIDX >= 0) {
			CylLorentzVectorF newJetPT = jets[jetOverlapIDX];
			if(lepton.pt() > newJetPT.pt() || (newJetPT - lepton).pt() < 20 ) newJetPT = CylLorentzVectorF();
			else newJetPT = newJetPT - lepton;
			jets[jetOverlapIDX] =newJetPT;
		}
		vector<std::pair<CylLorentzVectorF,float> > jetsNearW;
		for(unsigned int iJ= 0; iJ < jets.size(); ++iJ){
			if(jets[iJ].pt() < 20 ) continue;
			if(absDeltaPhi(jets[iJ],metPLep) > 0.8) continue;
			jetsNearW.emplace_back(jets[iJ],jet_csv->at(iJ) );
		}
		int nPairs = 0;
		int bestPairIJ1 = -1;
		int bestPairIJ2 = -1;
		for(unsigned int iJ1= 0; iJ1 < jetsNearW.size(); ++iJ1){
			for(unsigned int iJ2= iJ1+1; iJ2 < jetsNearW.size(); ++iJ2){
				if(deltaR2(jetsNearW[iJ1].first,jetsNearW[iJ2].first) > 1.0*1.0) continue;
				nPairs++;
				if(bestPairIJ1 || deltaR2(jetsNearW[iJ1].first+jetsNearW[iJ2].first,lepton) < deltaR2(jetsNearW[bestPairIJ1].first+jetsNearW[bestPairIJ2].first,lepton)  ){
					bestPairIJ1 = iJ1; bestPairIJ2 = iJ2;
				}
		}
		}

		TString exPref = prefix;
		CylLorentzVectorF vis;
		if(!jetsNearW.size()) return;
		else if(jetsNearW.size() == 1){
			if(jetsNearW[0].first.mass() > 30 && jetsNearW[0].second < CSV_MEDIUM){ exPref += "_ak41jpass"; vis = jetsNearW[0].first;}
			else return;
		} else{
			if(!nPairs) return;
			else if(jetsNearW[bestPairIJ1].second < CSV_MEDIUM && jetsNearW[bestPairIJ2].second < CSV_MEDIUM){ exPref += "_ak42jpass"; vis = jetsNearW[bestPairIJ1].first + jetsNearW[bestPairIJ2].first;}
			else return;
		}


		//find corrected hbb
		double hbbAlpha = 125/hbb.mass();
		auto corrHbb = hbb * hbbAlpha;
		auto corrMET = met +hbb - hbb;
		auto corrSolH = getnz(corrMET,vis) + vis;
		CylLorentzVectorF corrSolvedHH = corrSolH+corrHbb;

		plotter.getOrMake1D(TString::Format("%s_corr_hWW_mass",exPref.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(corrSolH.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_corr_hh_mass",exPref.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrSolvedHH.mass(),nWeight);



	}


	void makeMETPlots(TString prefix, const CylLorentzVectorF& met, const CylLorentzVectorF& lepton, const unsigned int ihbb, const CylLorentzVectorF& hbb, const unsigned int iwjj,const CylLorentzVectorF& wjj){
		CylLorentzVectorF vis = lepton + wjj;
		double mt = transverseMass(lepton + wjj,met);


		plotter.getOrMake1D(TString::Format("%s_met",prefix.Data()),";#slash{E}_{T} [GeV]; arbitrary units",200,0,1000)->Fill(met.pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_dPhi_met_lep",prefix.Data()),";#Delta#phi(#slash{E}_{T},lepton); arbitrary units",320,0,3.2)->Fill(absDeltaPhi(lepton,met),nWeight);
		plotter.getOrMake1D(TString::Format("%s_mt_met_lep",prefix.Data()),";m_{T}(lepton) [GeV]; arbitrary units",200,0,1000)->Fill(transverseMass(lepton,met),nWeight);
		plotter.getOrMake1D(TString::Format("%s_mt_met_lepWjj",prefix.Data()),";m_{T}(lepton+wJJ) [GeV]; arbitrary units",200,0,1000)->Fill(mt,nWeight);

//
//		CylLorentzVectorF hWW;
//		double chisq = solver.massMinimization(lepton,wjj,met,hWW);
//		plotter.getOrMake1D(TString::Format("%s_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(hWW.mass(),nWeight);
//		plotter.getOrMake1D(TString::Format("%s_hWW_pt",prefix.Data()),";hWW p_{T} [GeV]; arbitrary units",200,0,1000)->Fill(hWW.pt(),nWeight);
//		CylLorentzVectorF hh = hWW + hbb;
//		plotter.getOrMake1D(TString::Format("%s_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(hh.mass(),nWeight);
//		plotter.getOrMake1D(TString::Format("%s_hh_pt",prefix.Data()),";hh p_{T} [GeV]; arbitrary units",200,0,1000)->Fill(hh.pt(),nWeight);
//		plotter.getOrMake1D(TString::Format("%s_hWW_chi2",prefix.Data()),";chi2; arbitrary units",2000,0,100000)->Fill(chisq,nWeight);
//
//		plotter.getOrMake1D(TString::Format("%s_hh_dPhi",prefix.Data()),";#Delta#phi(hbb,hWW); arbitrary units",320,0,3.2)->Fill(absDeltaPhi(hWW,hbb),nWeight);
//		plotter.getOrMake1D(TString::Format("%s_hh_dEta",prefix.Data()),";#Delta#eta(hbb,hWW); arbitrary units",500,0,5.0)->Fill(absDeltaEta(hWW,hbb),nWeight);




		auto solH = getnz(met,vis) + vis;
		CylLorentzVectorF solvedHH = solH+hbb;


		plotter.getOrMake1D(TString::Format("%s_solved_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_solved_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(solvedHH.mass(),nWeight);

		plotter.getOrMake1D(TString::Format("%s_minW_mass",prefix.Data()),";min W mass [GeV]; arbitrary units",200,0,1000)->Fill(std::min(solH.mass(),wjj.mass()),nWeight);
		plotter.getOrMake1D(TString::Format("%s_maxW_mass",prefix.Data()),";max W mass [GeV]; arbitrary units",200,0,1000)->Fill(std::max(solH.mass(),wjj.mass()),nWeight);


		TString extra = "hh_mass_lt1200";
		if(solvedHH.mass() > 1400) extra= "hh_mass_geq1400";
		plotter.getOrMake1D(TString::Format("%s_%s_met",prefix.Data(),extra.Data()),";#slash{E}_{T} [GeV]; arbitrary units",200,0,1000)->Fill(met.pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_mt_met_lep",prefix.Data(),extra.Data()),";m_{T}(lepton) [GeV]; arbitrary units",200,0,1000)->Fill(transverseMass(lepton,met),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_lepW_pt",prefix.Data(),extra.Data()),";lep. W p_{T} [GeV]; arbitrary units",200,0,1000)->Fill((lepton+met).pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_lepW_pt_o_hWW_pt",prefix.Data(),extra.Data()),";lep. W p_{T} / hWW p{T}; arbitrary units",200,0,2)->Fill((lepton+met).pt()/solH.pt(),nWeight);

		plotter.getOrMake1D(TString::Format("%s_%s_wjj_fj_sd_mass"  ,prefix.Data(),extra.Data()),";wjj  softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(wjj.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_wjj_fj_tau2otau1",prefix.Data(),extra.Data()),";wjj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(iwjj) == 0 ? 99.0 : fj_tau2->at(iwjj)/fj_tau1->at(iwjj),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_wjj_fj_maxSJCSV" ,prefix.Data(),extra.Data()),";wjj max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(iwjj),fj_sj2_csv->at(iwjj)),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_wjj_fj_minSJCSV" ,prefix.Data(),extra.Data()),";wjj min subjet csv; arbitrary units",100,0,1)->Fill(std::min(fj_sj1_csv->at(iwjj),fj_sj2_csv->at(iwjj)),nWeight);

		plotter.getOrMake1D(TString::Format("%s_%s_hbb_fj_mass"  ,prefix.Data(),extra.Data()),";hbb mass [GeV]; arbitrary units",250,0,500)->Fill(fj_mass->at(ihbb),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_hbb_fj_sd_mass"  ,prefix.Data(),extra.Data()),";hbb  softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(hbb.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_hbb_fj_tau2otau1",prefix.Data(),extra.Data()),";hbb #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(fj_tau1->at(ihbb) == 0 ? 99.0 : fj_tau2->at(ihbb)/fj_tau1->at(ihbb),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_hbb_fj_maxSJCSV" ,prefix.Data(),extra.Data()),";hbb max subjet csv; arbitrary units",100,0,1)->Fill(std::max(fj_sj1_csv->at(ihbb),fj_sj2_csv->at(ihbb)),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_hbb_fj_minSJCSV" ,prefix.Data(),extra.Data()),";hbb min subjet csv; arbitrary units",100,0,1)->Fill(std::min(fj_sj1_csv->at(ihbb),fj_sj2_csv->at(ihbb)),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_nFJ"  ,prefix.Data(),extra.Data()),";N. of fat jets; arbitrary units",15,-0.5,14.5)->Fill(fj_pt->size(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_nAK4"  ,prefix.Data(),extra.Data()),";N. of ak4 jets; arbitrary units",15,-0.5,14.5)->Fill(jet_pt->size(),nWeight);


		float maxCSV = 0;
		int nExtraAK4 = 0;
		for(unsigned int iJ = 0; iJ < jet_pt->size(); ++iJ){
			if(deltaR2(jet_eta->at(iJ),jet_phi->at(iJ),hbb.eta(),hbb.phi() ) < 0.8*0.8 ) continue;
			if(deltaR2(jet_eta->at(iJ),jet_phi->at(iJ),wjj.eta(),wjj.phi() ) < 0.8*0.8 ) continue;
			if(deltaR2(jet_eta->at(iJ),jet_phi->at(iJ),lepton.eta(),lepton.phi() ) < 0.4*0.4 ) continue;
			nExtraAK4++;
			maxCSV = std::max(jet_csv->at(iJ),maxCSV);
		}

		plotter.getOrMake1D(TString::Format("%s_%s_nExtraAK4"  ,prefix.Data(),extra.Data()),";N. of extra ak4 jets; arbitrary units",15,-0.5,14.5)->Fill(nExtraAK4,nWeight);
		plotter.getOrMake1D(TString::Format("%s_%s_max_extra_ak4" ,prefix.Data(),extra.Data()),";extra ak4 max csv; arbitrary units",100,0,1)->Fill(maxCSV,nWeight);


		double minCSV = std::min(fj_sj1_csv->at(ihbb),fj_sj2_csv->at(ihbb));
		if(minCSV < CSV_MEDIUM)
			plotter.getOrMake1D(TString::Format("%s_loose_solved_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(solvedHH.mass(),nWeight);
		else
			plotter.getOrMake1D(TString::Format("%s_tight_solved_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(solvedHH.mass(),nWeight);



		//find corrected hbb
		double hbbAlpha = 125/hbb.mass();
		auto corrHbb = hbb * hbbAlpha;
		auto corrMET = met +hbb - corrHbb;
		auto corrSolH = getnz(corrMET,vis) + vis;
		CylLorentzVectorF corrSolvedHH = corrSolH+corrHbb;
		CylLorentzVectorF corrPartSolvedHH = solH+corrHbb;

		plotter.getOrMake1D(TString::Format("%s_corrPart_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrPartSolvedHH.mass(),nWeight);

		plotter.getOrMake1D(TString::Format("%s_corr_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(corrSolH.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_corr_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrSolvedHH.mass(),nWeight);
		if(minCSV < CSV_MEDIUM){
			plotter.getOrMake1D(TString::Format("%s_loose_corrPart_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrPartSolvedHH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_loose_corr_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrSolvedHH.mass(),nWeight);
			if(isSignal()) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",3,-0.5,2.5)->Fill(1.0,nWeight);

		}else{
			plotter.getOrMake1D(TString::Format("%s_tight_corrPart_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrPartSolvedHH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_tight_corr_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrSolvedHH.mass(),nWeight);
			if(isSignal()) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",3,-0.5,2.5)->Fill(2.0,nWeight);

		}

		CylLorentzVectorF b1(fj_sj1_pt->at(ihbb),fj_sd_eta->at(ihbb),fj_sd_phi->at(ihbb),0);
		CylLorentzVectorF b2(fj_sj2_pt->at(ihbb),fj_sd_eta->at(ihbb),fj_sd_phi->at(ihbb),0);
		double minLepMT = std::min(transverseMass(b1+lepton,met),transverseMass(b2+lepton,met));
		double minHadMass = std::min((wjj+b1).mass(),(wjj+b2).mass());
		double minLepMass = std::min((getnz(corrMET,vis)+lepton+b1).mass(),(getnz(corrMET,vis)+lepton+b2).mass());
		plotter.getOrMake1D(TString::Format("%s_minLepMT"  ,prefix.Data()),";minLepMT [GeV]; arbitrary units",500,0,1000)->Fill(minLepMT,nWeight);
		plotter.getOrMake1D(TString::Format("%s_minHadMass"  ,prefix.Data()),";minHadMass[GeV]; arbitrary units",500,0,1000)->Fill(minHadMass,nWeight);
		plotter.getOrMake1D(TString::Format("%s_minLepMass"  ,prefix.Data()),";minLepMass [GeV]; arbitrary units",500,0,1000)->Fill(minLepMass,nWeight);


	}



	virtual void runAEvent() {
		TString prefix = "";
		switch ( process ){
		case 2:
			prefix = "ttbar";
			break;
		case 3:
			prefix = "wlnu";
			break;
		case 6:
			prefix = "t";
			break;
		default:
			prefix = "signal_" +  	glbPrefix;
		}

		nWeight = lumi*weight;
		if(isSignal()){
			if(glbPrefix == "1000"){
				nWeight = lumi * 0.8241887906 * 20 / 1710;
			} else {
				nWeight = lumi * 0.8241887906 * 5 / 7940;
			}
		}
		if(isSignal()) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",3,-0.5,2.5)->Fill(0.0,nWeight);

		if(process == 10 && decayType != 5 ) return;
		if(ht < 400) return;



//		makeLeptonPlots(prefix);
		if(!lepton_pt->size() || lepton_pt->at(0) < 20) return;
//		makeHbbPlots(prefix);

		CylLorentzVectorF lepton = getLep();
		CylLorentzVectorF metV = getMET();
		CylLorentzVectorF metPLep = lepton+metV;

		int hbbIDX = -1;
		for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
			CylLorentzVectorF jet = getFJ(iJ);
			float dPhi =  std::fabs(deltaPhi(metPLep.phi(),jet.phi()));
			if(dPhi < TMath::PiOver2()) continue;

			CylLorentzVectorF sdJet = getSDFJ(iJ);
			bool passMass = (sdJet.mass() > 100 && sdJet.mass() < 150);
			bool passTau = (fj_tau1->at(iJ) == 0 ? 99.0 : fj_tau2->at(iJ)/fj_tau1->at(iJ)) < 0.7;
			bool passCSV = (std::max(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)) > CSV_MEDIUM && std::min(fj_sj1_csv->at(iJ),fj_sj2_csv->at(iJ)) > CSV_LOOSE );
			if(!passMass || !passTau || !passCSV) continue;
			if(hbbIDX < 0 || jet.pt() > fj_pt->at(hbbIDX)) hbbIDX = iJ;
		}

		if(hbbIDX < 0) return;
		makeWjjPlots(prefix);
//		makeAK4WjjPlots(prefix);
//		makeAK4METPlots(prefix,  metV, lepton, hbbIDX, getSDFJ(hbbIDX));

		int wJJIDX = -1;
		double minDR2 = -1;
		for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
			auto jet = getFJ(iJ);
			const double absDPhi = absDeltaPhi(jet,metPLep);
			if(absDPhi > 0.6) continue;
			const double dr2 = deltaR2(jet,lepton);
			if(wJJIDX < 0 || dr2 < minDR2){
				minDR2 =dr2;
				wJJIDX = iJ;
			}
		}
		if(wJJIDX < 0) return;

		CylLorentzVectorF wjjLepSub = getFJ(wJJIDX);
		CylLorentzVectorF sdWjjLepSub = getSDFJ(wJJIDX);

		const double deltaRLep = deltaR(wjjLepSub,lepton);
		if(deltaRLep < 0.8){
			wjjLepSub = lepton.pt() > wjjLepSub.pt() ? CylLorentzVectorF() : wjjLepSub - lepton;
			sdWjjLepSub = lepton.pt() > sdWjjLepSub.pt() ? CylLorentzVectorF() : sdWjjLepSub - lepton;
		}

		bool goodWJJ = (sdWjjLepSub.pt() > 0 && sdWjjLepSub.mass() > 10 && sdWjjLepSub.mass() < 120 && std::max(fj_sj1_csv->at(wJJIDX),fj_sj2_csv->at(wJJIDX)) < CSV_TIGHT);
		if(!goodWJJ) return;

		makeMETPlots(prefix,metV,lepton,hbbIDX,getSDFJ(hbbIDX),wJJIDX,sdWjjLepSub);


	}
};



#endif

void AnalyzeCutTree(std::string fileName, std::string prefix, std::string outName){

	Analyzer a (fileName,"Events");
	a.glbPrefix = prefix;
	a.analyze();
	a.write(outName);
}
