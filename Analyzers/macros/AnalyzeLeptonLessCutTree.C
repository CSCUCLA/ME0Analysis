
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


class Analyzer : public BaseCutAnalyzer {
public:


	Analyzer(std::string fileName, std::string treeName) : BaseCutAnalyzer(fileName,treeName){}


	void fillGenHbbPlots(TString prefix){
		CylLorentzVectorF hbb_b1(b1_pt,b1_eta,b1_phi,b1_mass);
		CylLorentzVectorF hbb_b2(b2_pt,b2_eta,b2_phi,b2_mass);

		CylLorentzVectorF wjj_genLep(w1d1_pt,w1d1_eta,w1d1_phi,w1d1_mass);
		CylLorentzVectorF wjj_d1(w2d1_pt,w2d1_eta,w2d1_phi,w2d1_mass);
		CylLorentzVectorF wjj_d2(w2d2_pt,w2d2_eta,w2d2_phi,w2d2_mass);


		int cat = 4;
		if(std::fabs(genHBB.eta()) > 2.4) cat = 0;
		else if (genHBB.pt() < 100) cat = 1;
		else if (deltaR2(hbb_b1,hbb_b2) > 0.8*0.8) cat = 2;
		else if(deltaR2(wjj_genLep,genHBB) <= deltaR2(wjj_genLep,genWjj)) cat=3;

		if(cat > 0)plotter.getOrMake2D(TString::Format("%s_gen_Hbb_pt_v_dr",prefix.Data()),";gen hbb p_{T}; dR(bb)",100,0,1000,100,0,3.2)->Fill(genHBB.pt(),deltaR(hbb_b1,hbb_b2),nWeight);
		plotter.getOrMake1D(TString::Format("%s_gen_Hbb_gen",prefix.Data()),";eta>2.4/pt<100/dr(bb)>0.8/dr(genlep,genHbb)<dr(genlep,genWjj)/good; arbitrary units",5,-0.5,4.5)->Fill(cat,nWeight);

		int sigClosestIDX = -1;
		for(unsigned int iJ = 0; iJ < fatjets.size(); ++iJ){
		if(sigClosestIDX < 0 || deltaR2(genHBB,fatjets[iJ]) < deltaR2(genHBB,fatjets[sigClosestIDX])  ) sigClosestIDX = iJ; }

		std::vector<std::pair<float,int>> rankedInclJets;
		for(unsigned int iJ = 0; iJ < fatjets.size(); ++iJ){ rankedInclJets.emplace_back(fatjets[iJ].pt(),iJ);}
		std::sort(rankedInclJets.begin(), rankedInclJets.end(),PhysicsUtilities::greaterFirst<float,int>());

		int inclRankOfGoodJet = -1;
		for(unsigned int iJ = 0; iJ < rankedInclJets.size(); ++iJ){
			if(rankedInclJets[iJ].second != sigClosestIDX ) continue;
			inclRankOfGoodJet = iJ;
			break;
		}

		int recoCat = 4;
		if(cat != 4) recoCat = 0;
		else if (sigClosestIDX < 0 ) recoCat = 1;
		else if (inclRankOfGoodJet < 0 || inclRankOfGoodJet > 1 ) recoCat = 2;
		else if(fatjets.size() > 1) {
			const auto& genHbbFJ =  fatjets[sigClosestIDX];
			const auto otherJet = fatjets[ rankedInclJets[inclRankOfGoodJet == 0 ? 1 : 0].second];
			if(deltaR2(lepton,genHbbFJ) <= deltaR2(lepton,otherJet)) recoCat = 3;
		}

		plotter.getOrMake1D(TString::Format("%s_gen_Hbb_reco",prefix.Data()),";fail gen/NoJet/rank < 1/closer to lep/good; arbitrary units",5,-0.5,4.5)->Fill(recoCat,nWeight);
		plotter.getOrMake1D(TString::Format("%s_genHbbfj_inclrank",prefix.Data()),";gen hbb fj rank; arbitrary units",11,-1.5,9.5)->Fill(inclRankOfGoodJet,nWeight);
		if(sigClosestIDX >= 0) plotter.getOrMake1D(TString::Format("%s_gen_dR_genHbb_hbb",prefix.Data()),";Min[#DeltaR(gen hbb,fj)]; arbitrary units",320,0,3.2)->Fill( deltaR(genHBB, fatjets[sigClosestIDX]),nWeight);

		int wjjcat = 4;
		if(std::fabs(genWjj.eta()) > 2.4) wjjcat = 0;
		else if (genWjj.pt() < 100) wjjcat = 1;
		else if (deltaR2(wjj_d1,wjj_d2) > 0.8*0.8) wjjcat = 2;
		else if(deltaR2(wjj_genLep,genWjj) >= deltaR2(wjj_genLep,genHBB)) wjjcat=3;


		int wjjClosestIDX = -1;
		for(unsigned int iJ = 0; iJ < fatjets.size(); ++iJ){
		if(wjjClosestIDX < 0 || deltaR2(genWjj,fatjets[iJ]) < deltaR2(genWjj,fatjets[wjjClosestIDX])  ) wjjClosestIDX = iJ; }
		int rankOfWJJJet = -1;
		for(unsigned int iJ = 0; iJ < rankedInclJets.size(); ++iJ){
			if(rankedInclJets[iJ].second != wjjClosestIDX ) continue;
			rankOfWJJJet = iJ;
			break;
		}

		int recoWjjCat = 4;
		if(recoCat != 4) recoWjjCat =0;
		else if(wjjcat != 4) recoWjjCat = 1;
		else if (wjjClosestIDX < 0 ) recoWjjCat = 2;
		else if (rankOfWJJJet < 0 || rankOfWJJJet > 1 ) recoWjjCat = 3;

		if(wjjcat > 0)plotter.getOrMake2D(TString::Format("%s_gen_wjj_pt_v_dr",prefix.Data()),";gen wjj p_{T}; dR(jj)",100,0,1000,100,0,3.2)->Fill(genWjj.pt(),deltaR(wjj_d1,wjj_d2),nWeight);
		plotter.getOrMake1D(TString::Format("%s_gen_wjj_gen",prefix.Data()),";eta>2.4/pt<100/dr(jj)>0.8/dr(genlep,genHbb)<dr(genlep,genWjj)/good; arbitrary units",5,-0.5,4.5)->Fill(wjjcat,nWeight);
		plotter.getOrMake1D(TString::Format("%s_gen_wjj_reco",prefix.Data()),";fail hbb /fail gen/NoJet/rank < 1/goodd; arbitrary units",5,-0.5,4.5)->Fill(recoWjjCat,nWeight);
		if(wjjcat >= 2)
			plotter.getOrMake1D(TString::Format("%s_genWjjfj_rank",prefix.Data()),";gen wjj fj rank; arbitrary units",11,-1.5,9.5)->Fill(rankOfWJJJet,nWeight);

		if(true){
			float minDR2 = 0;
			int wJJIDX = -1;
			bool passHbb = false;
			bool passWjj = false;
			for(unsigned int iJ = 0; iJ < fj_pt->size(); ++iJ){
				CylLorentzVectorF jet = getFJ(iJ);
				float dPhi =  std::fabs(deltaPhi(lepMET.phi(),jet.phi()));
				if(dPhi >= TMath::PiOver2() && iJ == sigClosestIDX) passHbb = true;
				if(dPhi > 0.6) continue;
				const double dr2 = deltaR2(jet,lepton);
				if(wJJIDX < 0 || dr2 < minDR2){
					minDR2 =dr2;
					wJJIDX = iJ;
				}
			}
			if(wJJIDX > 0 && wJJIDX == wjjClosestIDX) passWjj = true;
			if(cat > 2 && wjjcat > 2 && passWjj && passHbb)
				plotter.getOrMake1D(TString::Format("%s_gen_Hbb_recoComp",prefix.Data()),";pass old/pass new",4,-0.5,3.5)->Fill(0.0,nWeight);
			if(recoWjjCat == 4 && recoCat == 4)
				plotter.getOrMake1D(TString::Format("%s_gen_Hbb_recoComp",prefix.Data()),";pass old/pass new",4,-0.5,3.5)->Fill(1.0,nWeight);
		}




	}

	void makeHbbPlots(TString prefix){
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_mass",prefix.Data()),";hbb fj mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_sd_mass",prefix.Data()),";hbb fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->sdMom.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_sdfj_mass",prefix.Data()),";hbb fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->sd_massFFJ,nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_tau2otau1",prefix.Data()),";hbb fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(hbbFJ->tau2otau1(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_csv",prefix.Data()),";hbb fj csv; arbitrary units",100,0,1)->Fill(hbbFJ->csv,nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_bbcsv",prefix.Data()),";hbb fj bb csv; arbitrary units",100,0,1)->Fill(hbbFJ->bbcsv,nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->minSJCSV(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_maxsdcsv",prefix.Data()),";hbb fj bb max sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_maxMed_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV() >= CSV_MEDIUM ? hbbFJ->minSJCSV() :0.0,nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_fj_maxTight_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV() >= CSV_TIGHT? hbbFJ->minSJCSV() :0.0,nWeight);


		bool passMass = (hbbFJ->sd_massFFJ > 90 && hbbFJ->sd_massFFJ < 140);
		bool passTau = hbbFJ->tau2otau1() < 0.5;
		bool passCSV = (hbbFJ->maxSJCSV() > CSV_MEDIUM && hbbFJ->minSJCSV() > CSV_LOOSE );

		if(passMass && passTau){
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_csv",prefix.Data()),";hbb fj csv; arbitrary units",100,0,1)->Fill(hbbFJ->csv,nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_bbcsv",prefix.Data()),";hbb fj bb csv; arbitrary units",100,0,1)->Fill(hbbFJ->bbcsv,nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->minSJCSV(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_maxsdcsv",prefix.Data()),";hbb fj bb max sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_maxMed_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV() >= CSV_MEDIUM ? hbbFJ->minSJCSV() :0.0,nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_maxTight_minsdcsv",prefix.Data()),";hbb fj bb min sj csv; arbitrary units",100,0,1)->Fill(hbbFJ->maxSJCSV() >= CSV_TIGHT? hbbFJ->minSJCSV() :0.0,nWeight);
		}
		if(passMass&&passCSV){
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_tau2otau1",prefix.Data()),";hbb fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(hbbFJ->tau2otau1(),nWeight);
		}
		if(passTau&&passCSV){
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_mass",prefix.Data()),";hbb fj mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_sd_mass",prefix.Data()),";hbb fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->sdMom.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_sdfj_mass",prefix.Data()),";hbb fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->sd_massFFJ,nWeight);
		}
		if(passTau&&passCSV&&passMass){
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_rawmass",prefix.Data()),";hbb fj mass [GeV]; arbitrary units",250,0,500)->Fill(hbbFJ->mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_tau3otau2",prefix.Data()),";hbb fj #tau_{3}/#tau_{2}; arbitrary units",100,0,1)->Fill(hbbFJ->tau3otau2(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hbb_fj_oM_tau3otau1",prefix.Data()),";hbb fj #tau_{3}/#tau_{1}; arbitrary units",100,0,1)->Fill(hbbFJ->tau3otau1(),nWeight);

		}
	}

	void makeWjjPlots(TString prefix){
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_mass",prefix.Data()),";wjj fj mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_sd_mass",prefix.Data()),";wjj fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->sdMom.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_sdfj_mass",prefix.Data()),";wjj fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->sd_massFFJ,nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_tau2otau1",prefix.Data()),";wjj fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(wjjFJ->tau2otau1(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_csv",prefix.Data()),";wjj fj csv; arbitrary units",100,0,1)->Fill(wjjFJ->csv,nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_minsdcsv",prefix.Data()),";wjj fj min sj csv; arbitrary units",100,0,1)->Fill(wjjFJ->minSJCSV(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wjj_fj_maxsdcsv",prefix.Data()),";wjj fj max sj csv; arbitrary units",100,0,1)->Fill(wjjFJ->maxSJCSV(),nWeight);


		bool passMass = (wjjFJ->sd_massFFJ > 10);
		bool passTau = wjjFJ->tau2otau1() < 0.65;
		bool passCSV = (wjjFJ->maxSJCSV()  < CSV_TIGHT );

		if(passMass && passTau){
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_csv",prefix.Data()),";wjj fj csv; arbitrary units",100,0,1)->Fill(wjjFJ->csv,nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_minsdcsv",prefix.Data()),";wjj fj min sj csv; arbitrary units",100,0,1)->Fill(wjjFJ->minSJCSV(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_maxsdcsv",prefix.Data()),";wjj fj max sj csv; arbitrary units",100,0,1)->Fill(wjjFJ->maxSJCSV(),nWeight);
		}
		if(passMass&&passCSV){
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_tau2otau1",prefix.Data()),";wjj fj #tau_{2}/#tau_{1}; arbitrary units",100,0,1)->Fill(wjjFJ->tau2otau1(),nWeight);
		}
		if(passTau&&passCSV){
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_mass",prefix.Data()),";wjj fj mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_sd_mass",prefix.Data()),";wjj fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->sdMom.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_wjj_fj_oM_sdfj_mass",prefix.Data()),";wjj fj softdrop mass [GeV]; arbitrary units",250,0,500)->Fill(wjjFJ->sd_massFFJ,nWeight);
		}
	}
	void makeHHPlots(TString prefix){
		auto lepW = neutrino +lepton;

		bool tightHbb = hbbFJ->minSJCSV() > CSV_MEDIUM;
		bool deltaRCut = deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5;
		bool onShellWjj = wjjFJ->sd_massFFJ > 60;

		plotter.getOrMake1D(TString::Format("%s_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);

		if(tightHbb){
			plotter.getOrMake1D(TString::Format("%s_tightHbb_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_tightHbb_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		} else {
			plotter.getOrMake1D(TString::Format("%s_looseHbb_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_looseHbb_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		}
		if(deltaRCut){
		if(tightHbb){
			plotter.getOrMake1D(TString::Format("%s_DR_tightHbb_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_DR_tightHbb_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_DR_tightHbb_hbb_mass",prefix.Data()),";hbb mass [GeV]; arbitrary units",200,0,1000)->Fill(hbbFJ->sd_massFFJ,nWeight);
		} else {
			plotter.getOrMake1D(TString::Format("%s_DR_looseHbb_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_DR_looseHbb_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		}
		}

		if(onShellWjj){
			plotter.getOrMake1D(TString::Format("%s_onShellWjj_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_onShellWjj_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		} else {
			plotter.getOrMake1D(TString::Format("%s_offShellWjj_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_offShellWjj_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		}

		if(onShellWjj && tightHbb){
			plotter.getOrMake1D(TString::Format("%s_tightH_onShell_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_tightH_onShell_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		} else if(onShellWjj && !tightHbb){
			plotter.getOrMake1D(TString::Format("%s_looseH_onShell_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_looseH_onShell_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		} else if(!onShellWjj && tightHbb){
			plotter.getOrMake1D(TString::Format("%s_tightH_offShell_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_tightH_offShell_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		} else{
			plotter.getOrMake1D(TString::Format("%s_looseH_offShell_hWW_mass",prefix.Data()),";hWW mass [GeV]; arbitrary units",200,0,1000)->Fill(solH.mass(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_looseH_offShell_hh_mass",prefix.Data()),";hh mass [GeV]; arbitrary units",600,0,3000)->Fill(corrHH.mass(),nWeight);
		}


		plotter.getOrMake1D(TString::Format("%s_neutrino_pt",prefix.Data()),";neutrino p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( neutrino.pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_neutrino_eta",prefix.Data()),";neutrino |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs(neutrino.eta()),nWeight);
		plotter.getOrMake1D(TString::Format("%s_neutrino_lepton_dEta",prefix.Data()),";|#Delta#eta(neutrino,lepton)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta(neutrino,lepton),nWeight);
		plotter.getOrMake1D(TString::Format("%s_neutrino_lepton_dR",prefix.Data()),";|#DeltaR(neutrino,lepton)|; arbitrary units",50,0,5.0)->Fill( deltaR(neutrino,lepton),nWeight);

		plotter.getOrMake1D(TString::Format("%s_wlnu_pt",prefix.Data()),";wlnu p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( (neutrino+lepton).pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wlnu_eta",prefix.Data()),";wlnu |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs((neutrino+lepton).eta()),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wlnu_mass",prefix.Data()),";wlnu |#eta|; arbitrary units",100,0,1000)->Fill( std::fabs((neutrino+lepton).mass()),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wlnu_wjj_dEta",prefix.Data()),";|#Delta#eta(wlnu,wjj)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta((neutrino+lepton),wjjFJ->sdMom),nWeight);
		plotter.getOrMake1D(TString::Format("%s_wlnu_wjj_dR",prefix.Data()),";|#DeltaR(wlnu,wjj)|; arbitrary units",50,0,5.0)->Fill( deltaR((neutrino+lepton),wjjFJ->sdMom),nWeight);

		plotter.getOrMake1D(TString::Format("%s_wlnu_wjj_rPTo2",prefix.Data()),";|#DeltaR(wlnu,wjj)*p_{T}/2|; arbitrary units",100,0,1000)->Fill( deltaR(lepW,wjjFJ->sdMom) *solH.pt()/2 ,nWeight);

		plotter.getOrMake1D(TString::Format("%s_hWW_pt",prefix.Data()),";hWW p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( solH.pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hWW_eta",prefix.Data()),";hWW |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs(solH.eta()),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hWW_hbb_dEta",prefix.Data()),";|#Delta#eta(hWW,hbb)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta(solH,hbbFJ->sdMom),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hWW_hbb_dR",prefix.Data()),";|#DeltaR(hWW,hbb)|; arbitrary units",50,0,5.0)->Fill( deltaR(solH,hbbFJ->sdMom),nWeight);

		plotter.getOrMake1D(TString::Format("%s_wJJ_pt",prefix.Data()),";wjj p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( wjjFJ->pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hbb_pt",prefix.Data()),";hbb p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( hbbFJ->pt(),nWeight);
		plotter.getOrMake1D(TString::Format("%s_dR_wJJ_lepton",prefix.Data()),";#DeltaR(wjj,lep)]; arbitrary units",640,0,6.4)->Fill( deltaR(*wjjFJ, lepton),nWeight);
		plotter.getOrMake1D(TString::Format("%s_dR_hbb_lepton",prefix.Data()),";#DeltaR(hbb,lep)]; arbitrary units",640,0,6.4)->Fill( deltaR(*hbbFJ, lepton),nWeight);
		plotter.getOrMake1D(TString::Format("%s_hm_dPhi_hbb_lepton",prefix.Data()),";#Delta#phi(hbb,wJJ)]; arbitrary units",320,0,3.2)->Fill( absDeltaPhi(*hbbFJ, lepton),nWeight);

		if(corrHH.mass() > 1400 && corrHH.mass() < 1600){
			plotter.getOrMake1D(TString::Format("%s_hm_wJJ_pt",prefix.Data()),";wjj p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( wjjFJ->pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_hbb_pt",prefix.Data()),";hbb p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( hbbFJ->pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_dR_wJJ_lepton",prefix.Data()),";#DeltaR(wjj,lep)]; arbitrary units",640,0,6.4)->Fill( deltaR(*wjjFJ, lepton),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_dR_hbb_lepton",prefix.Data()),";#DeltaR(hbb,lep)]; arbitrary units",640,0,6.4)->Fill( deltaR(*hbbFJ, lepton),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_dPhi_hbb_lepton",prefix.Data()),";#Delta#phi(hbb,wJJ)]; arbitrary units",320,0,3.2)->Fill( absDeltaPhi(*hbbFJ, lepton),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_neutrino_pt",prefix.Data()),";neutrino p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( neutrino.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_neutrino_eta",prefix.Data()),";neutrino |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs(neutrino.eta()),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_neutrino_lepton_dEta",prefix.Data()),";|#Delta#eta(neutrino,lepton)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta(neutrino,lepton),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_neutrino_lepton_dR",prefix.Data()),";|#DeltaR(neutrino,lepton)|; arbitrary units",50,0,5.0)->Fill( deltaR(neutrino,lepton),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_pt",prefix.Data()),";wlnu p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( (neutrino+lepton).pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_eta",prefix.Data()),";wlnu |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs((neutrino+lepton).eta()),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_mass",prefix.Data()),";wlnu |#eta|; arbitrary units",100,0,1000)->Fill( std::fabs((neutrino+lepton).mass()),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_wjj_dEta",prefix.Data()),";|#Delta#eta(wlnu,wjj)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta((neutrino+lepton),wjjFJ->sdMom),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_wjj_dR",prefix.Data()),";|#DeltaR(wlnu,wjj)|; arbitrary units",50,0,5.0)->Fill( deltaR((neutrino+lepton),wjjFJ->sdMom),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_wlnu_wjj_rPTo2",prefix.Data()),";|#DeltaR(wlnu,wjj)*p_{T}/2|; arbitrary units",100,0,1000)->Fill( deltaR(lepW,wjjFJ->sdMom) *solH.pt()/2 ,nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_hWW_pt",prefix.Data()),";hWW p_{T} [GeV]; arbitrary units",200,0,2000)->Fill( solH.pt(),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_hWW_eta",prefix.Data()),";hWW |#eta|; arbitrary units",50,0,5.0)->Fill( std::fabs(solH.eta()),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_hWW_hbb_dEta",prefix.Data()),";|#Delta#eta(hWW,hbb)|; arbitrary units",50,0,5.0)->Fill( absDeltaEta(solH,hbbFJ->sdMom),nWeight);
			plotter.getOrMake1D(TString::Format("%s_hm_hWW_hbb_dR",prefix.Data()),";|#DeltaR(hWW,hbb)|; arbitrary units",50,0,5.0)->Fill( deltaR(solH,hbbFJ->sdMom),nWeight);
		}

		if(process == 2){
			auto tops = getTTBar();
			plotter.getOrMake2D(TString::Format("%s_tttmass_v_recoMass",prefix.Data()),"; gen tt mass [GeV];reco hh mass",200,0,2000,200,0,2000)->Fill( (tops.first.top + tops.second.top).mass(),corrHH.mass(),nWeight);

		}



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
				nWeight = lumi * 0.8241887906 * 20 / 8000;
			} else {
				nWeight = lumi * 0.8241887906 * 5 / 7940;
			}
		}
		plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5)->Fill(0.0,nWeight);
		plotter.getOrMake2D(TString::Format("%s_eff_byDecay",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5,6,-0.5,5.5)->Fill(0.0,decayType,nWeight);

//		if(process == 10 && decayType < 4 ) return;
		if(ht < 400) return;
		if(lepton_pt->size() != 1 || lepton_pt->at(0) < 20) return;

		plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5)->Fill(1.0,nWeight);
		plotter.getOrMake2D(TString::Format("%s_eff_byDecay",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5,6,-0.5,5.5)->Fill(1.0,decayType,nWeight);

		fillFJJets(fatjets);
		fillFJ1p5Jets(fatjets1p5);
		if(isSignal()){
			genLEPMET =CylLorentzVectorF(w1d1_pt,w1d1_eta,w1d1_phi,w1d1_mass)+CylLorentzVectorF(w1d2_pt,w1d2_eta,w1d2_phi,w1d2_mass);
			genHBB= CylLorentzVectorF(hbb_pt,hbb_eta,hbb_phi,hbb_mass);
			genWjj= CylLorentzVectorF(w2d1_pt,w2d1_eta,w2d1_phi,w2d1_mass)+CylLorentzVectorF(w2d2_pt,w2d2_eta,w2d2_phi,w2d2_mass);
		}

		lepton = getLep();
		metV = getMET();
		lepMET = lepton+metV;

		passHBB = false;
		passWJJ = false;
		hbbFJ = 0;
		wjjFJ = 0;
		if(fatjets.size() >= 2) {
			int iJHbb = -1;
			int iJWjj = -1;
			for(unsigned int iJ = 0; iJ < fatjets.size(); ++iJ){
				double dr2 = deltaR2(fatjets[iJ],lepton);
				if(dr2 > 2*2){
					if(iJHbb < 0 ||fatjets[iJ].pt() > fatjets[iJHbb].pt()  )
						iJHbb = iJ;
				} else if(dr2 < 0.8*0.8){
					if(iJWjj < 0 || dr2 <deltaR2(fatjets[iJWjj],lepton) )
						iJWjj = iJ;
				}
			}
			if(iJWjj>= 0) wjjFJ = &fatjets[iJWjj];
			if(iJHbb>= 0) hbbFJ = &fatjets[iJHbb];
//			std::sort(fatjets.begin(),fatjets.end(),PhysicsUtilities::greaterPT<FatJet>());
//			if(deltaR2(fatjets[0],lepton) < deltaR2(fatjets[1],lepton)){
//				wjjFJ = &fatjets[0];
//				hbbFJ = &fatjets[1];
//			} else {
//				wjjFJ = &fatjets[1];
//				hbbFJ = &fatjets[0];
//			}
			if(wjjFJ && wjjFJ->sd_massFFJ > 10 && wjjFJ->tau2otau1() < 0.65 && wjjFJ->maxSJCSV()  < CSV_TIGHT )
				passWJJ = true;
			if(hbbFJ && hbbFJ->sd_massFFJ > 90 && hbbFJ->sd_massFFJ < 140 && hbbFJ->tau2otau1() < 0.5 && hbbFJ->maxSJCSV()  > CSV_MEDIUM  && hbbFJ->minSJCSV() > CSV_LOOSE  )
				passHBB = true;
		}


		corrHbb  = CylLorentzVectorF();
		neutrino = CylLorentzVectorF();
		solH     = CylLorentzVectorF();
		corrHH   = CylLorentzVectorF();
		if(hbbFJ && wjjFJ && lepton.pt() > 0){
			double hbbAlpha = 125/hbbFJ->sdMom.mass();
			corrHbb= hbbFJ->sdMom * hbbAlpha;
			CylLorentzVectorF visHWW = lepton +wjjFJ->sdMom;
			neutrino =getnz(metV,visHWW );
			solH   = neutrino + visHWW;
			corrHH = solH + hbbFJ->sdMom;
		}

		if(isSignal() && decayType == 5) fillGenHbbPlots(prefix);
		if(hbbFJ) makeHbbPlots(prefix);
		if(passHBB && wjjFJ) makeWjjPlots(prefix);
		if(passWJJ && passHBB && lepton.pt() > 0) makeHHPlots(prefix);

		//round two
		if(passWJJ && hbbFJ && lepton.pt() > 0 ){
			if(corrHH.mass() > 900 && corrHH.mass() < 1100) makeHbbPlots(prefix + "_mass900to1100");
			if(corrHH.mass() > 1400 && corrHH.mass() < 1800) makeHbbPlots(prefix + "_mass1400to1800");
		}
		if(wjjFJ && passHBB && hbbFJ->minSJCSV() > CSV_MEDIUM && lepton.pt() > 0 ){
			if(corrHH.mass() > 900 && corrHH.mass() < 1100) makeWjjPlots(prefix + "_mass900to1100");
			if(corrHH.mass() > 1400 && corrHH.mass() < 1800) makeWjjPlots(prefix + "_mass1400to1800");
		}


		bool deltaRCut = wjjFJ && lepton.pt() > 0 && hbbFJ &&deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5;
		if(passHBB) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5)->Fill(2.0,nWeight);
		if(passHBB&&passWJJ&& deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5)->Fill(3.0,nWeight);
		if(passHBB&&passWJJ&& hbbFJ->minSJCSV() > CSV_MEDIUM&& deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5) plotter.getOrMake1D(TString::Format("%s_eff",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5)->Fill(4.0,nWeight);



		if(passHBB) plotter.getOrMake2D(TString::Format("%s_eff_byDecay",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5,6,-0.5,5.5)->Fill(2.0,decayType,nWeight);
		if(passHBB&&passWJJ&& deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5) plotter.getOrMake2D(TString::Format("%s_eff_byDecay",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5,6,-0.5,5.5)->Fill(3.0,decayType,nWeight);
		if(passHBB&&passWJJ&& hbbFJ->minSJCSV() > CSV_MEDIUM&& deltaR((neutrino+lepton),wjjFJ->sdMom) < 0.5) plotter.getOrMake2D(TString::Format("%s_eff_byDecay",prefix.Data()),";eff; arbitrary units",5,-0.5,4.5,6,-0.5,5.5)->Fill(4.0,decayType,nWeight);


	}

	std::vector<FatJet> fatjets;
	std::vector<FatJet> fatjets1p5;
	const FatJet * hbbFJ = 0;
	const FatJet * wjjFJ = 0;
	bool passHBB = false;
	bool passWJJ = false;

	CylLorentzVectorF corrHbb;
	CylLorentzVectorF neutrino;
	CylLorentzVectorF solH;
	CylLorentzVectorF corrHH;

	CylLorentzVectorF genLEPMET;
	CylLorentzVectorF genHBB;
	CylLorentzVectorF genWjj;
	CylLorentzVectorF lepton;
	CylLorentzVectorF  metV;
	CylLorentzVectorF lepMET;

};



#endif

void AnalyzeLeptonLessCutTree(std::string fileName, std::string prefix, std::string outName){

	Analyzer a (fileName,"Events");
	a.glbPrefix = prefix;
	a.analyze();
	a.write(outName);
}
