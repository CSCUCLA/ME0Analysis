// user include files

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/GEMRecHit/interface/ME0Segment.h"


#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../interface/ME0Helper.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class ME0DigiAnalyzer : public edm::EDAnalyzer {
public:
	explicit ME0DigiAnalyzer(const edm::ParameterSet&);
	~ME0DigiAnalyzer();


private:
	virtual void beginJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	void makeGeneralDigiPlots(const ME0Geometry* mgeom, const ME0DigiPreRecoCollection& digis, TString name);
	void compDigiPlots(const ME0Geometry* mgeom, const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname);
	void getDigiCuts(const ME0Geometry* mgeom,edm::Handle<std::vector<SimTrack>>&  simTrackH, edm::Handle<std::vector<PSimHit> >&  simHitH,
			const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis, const ME0DigiPreRecoMap& digiMap, TString sname);
	void fitSegment(const ME0Geometry* mgeom,  edm::Handle<std::vector<SimTrack>>&  simTrackH, edm::Handle<std::vector<PSimHit> >&  simHitH,
			const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname);
	virtual void endJob() {};


private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoMap>                 oldToNewMapToken_  ;
    edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
    edm::EDGetTokenT<std::vector<SimTrack>> track_token;



	TString outFileName;
	TString runName;
	bool runNewDigis;
	HistGetter hists;
};


ME0DigiAnalyzer::ME0DigiAnalyzer(const edm::ParameterSet& iConfig)
: outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
  runName(iConfig.getUntrackedParameter<std::string>("runName"))
{
	dToken_        = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
	runNewDigis        = iConfig.getParameter<std::string>("newDigiCollection") != "";
	newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
	oldToNewMapToken_  = consumes<ME0DigiPreRecoMap>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
	  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
	  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );

}


ME0DigiAnalyzer::~ME0DigiAnalyzer() {
	hists.write(outFileName);
}

void ME0DigiAnalyzer::makeGeneralDigiPlots(const ME0Geometry* mgeom, const ME0DigiPreRecoCollection& digis, TString sname){
	for(const auto& d: digis ){
		const auto* epart = mgeom->etaPartition(d.first);

		for (ME0DigiPreRecoCollection::const_iterator idigi = d.second.first;
				idigi != d.second.second;idigi++) {

			hists.getOrMake1D(TString::Format("%sdigiNumber",sname.Data()),";Total, Prompt, NonPrompt",3,-.5,2.5)->Fill(0);
			if(idigi->prompt())hists.getOrMake1D(TString::Format("%sdigiNumber",sname.Data()),";Total, Prompt, NonPrompt",3,-.5,2.5)->Fill(1);
			else hists.getOrMake1D(TString::Format("%sdigiNumber",sname.Data()),";Total, Prompt, NonPrompt",3,-.5,2.5)->Fill(2);

			hists.getOrMake1D(TString::Format("%stof",sname.Data()),";time of flight [ns]",400,-200,200)->Fill(idigi->tof());
			hists.getOrMake2D(TString::Format("%shitLoc",sname.Data()),";local x position [cm] ;local y position [cm]",4000,-30,30,100,-50,50)->Fill(idigi->x(),idigi->y());

			TString loc;
			if(idigi->x() < -15.0 )loc = "x_leqm15";
			else if (idigi->x() < 15.0 )loc = "x_m15to15";
			else if (idigi->x() < 15.0 )loc = "x_m15to15";
			else loc = "x_geq15";
			if (idigi->y() < -6 ) loc += "y_leqm6";
			else if (idigi->y() < 6 ) loc += "y_m6to6";
			else loc += "y_geq6";

			hists.getOrMake1D(TString::Format("%s%s_xError",sname.Data(),loc.Data()),";local x position error [cm]",2000,0,5)->Fill(idigi->ex());
			hists.getOrMake1D(TString::Format("%s%s_yError",sname.Data(),loc.Data()),";local y position error [cm]",100,0,20)->Fill(idigi->ey());




			int partID = TMath::Abs(idigi->pdgid());
			hists.getOrMake1D(TString::Format("%spartTypes",sname.Data()),";PDGID",2201,-.5,2200.5)->Fill(partID < 2199 ? partID : 2200);

			TString prefix = "";
			if(partID == 11) prefix = "ele";
			else if(partID == 2112) prefix = "neutron";
			else if(partID == 22) prefix = "photon";
			else prefix = "other";
			hists.getOrMake1D(TString::Format("%s%s_tof",sname.Data(),prefix.Data()),";tof",400,-200,200)->Fill(idigi->tof());
			hists.getOrMake2D(TString::Format("%s%s_byLay_tof",sname.Data(),prefix.Data()),";tof ;layer",400,-200,200,8,-0.5,7.5)->Fill(idigi->tof(),d.first.layer());
			hists.getOrMake1D(TString::Format("%s%s_hits",sname.Data(),prefix.Data()),";layer",8,-0.5,7.5)->Fill(d.first.layer());
			;
			hists.getOrMake2D(TString::Format("%s%s_byLay_hitsByETA",sname.Data(),prefix.Data()),";|eta| ;layer",500,0,5,8,-0.5,7.5)->Fill(TMath::Abs(epart->toGlobal(LocalPoint(idigi->x(),idigi->y(),0)).eta()),d.first.layer());
			if(d.first.layer() == 1)
				hists.getOrMake1D(TString::Format("%s%s_hitsByRadius",sname.Data(),prefix.Data()),";radius",1000,0,200)->Fill(epart->toGlobal(LocalPoint(idigi->x(),idigi->y(),0)).perp());
			if(d.first.layer() == 1)
				hists.getOrMake1D(TString::Format("%s%s_hitsByETA",sname.Data(),prefix.Data()),";|eta|",500,0,5)->Fill(TMath::Abs(epart->toGlobal(LocalPoint(idigi->x(),idigi->y(),0)).eta()));

		}

	}


}

void ME0DigiAnalyzer::compDigiPlots(const ME0Geometry* mgeom, const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname){
	for(const auto& d: oldDigis ){
		const auto* epart = mgeom->etaPartition(d.first);
		auto newDigiL = newDigis.get(d.first);
		auto oldToNewMap = digiMap.get(d.first);
		for (ME0DigiPreRecoCollection::const_iterator idigi = d.second.first;
				idigi != d.second.second;idigi++) {
			int newIDX =  oldToNewMap.first[int(idigi - d.second.first)];
			if(newIDX >= 0){
				const auto * newDigi = &newDigiL.first[newIDX];
				hists.getOrMake2D(TString::Format("%sodtonewd_tof",sname.Data()),";original time of flight [ns]; new time of flight [ns]",400,-200,200,400,-200,200)->Fill(idigi->tof(),newDigi->tof());
				hists.getOrMake1D(TString::Format("%sodtonewd_xdiff",sname.Data()),";local x position (original - new) [cm]",1000,-5,5)->Fill(idigi->x() - newDigi->x());
				hists.getOrMake1D(TString::Format("%sodtonewd_ydiff",sname.Data()),";local y position (original - new) [cm]",1000,-20,20)->Fill(idigi->y() - newDigi->y());

				auto oldG = epart->toGlobal(LocalPoint(idigi->x(),idigi->y(),0));
				auto newG = epart->toGlobal(LocalPoint(newDigi->x(),newDigi->y(),0));

				hists.getOrMake1D(TString::Format("%sodtonewd_phidiff",sname.Data()),";#phi position (original - new)",1000,-.003,.003)->Fill(TVector2::Phi_mpi_pi(oldG.phi() - newG.phi()) );
				hists.getOrMake1D(TString::Format("%sodtonewd_etadiff",sname.Data()),";#eta position (original - new)",1000,-1,1)->Fill(oldG.eta() - newG.eta());

				hists.getOrMake1D(TString::Format("%sodtonewd_xresid",sname.Data()),";local x position (original - new)/(new error)",1000,-5,5)->Fill((idigi->x() - newDigi->x())/newDigi->ex());
				hists.getOrMake1D(TString::Format("%sodtonewd_yresid",sname.Data()),";local y position (original - new)/(new error)",1000,-5,5)->Fill((idigi->y() - newDigi->y())/newDigi->ey());
			} else {
				if(TMath::Abs(idigi->tof()) > 30 )
					hists.getOrMake1D(TString::Format("%sodtonewd_whyFail",sname.Data()),";TOF,Neutron, Other",4,-.5,3.5)->Fill(0);
				else if(!idigi->prompt()) {
					hists.getOrMake1D(TString::Format("%sodtonewd_whyFail",sname.Data()),";TOF,Neutron, Other",4,-.5,3.5)->Fill(1);
				} else{
					hists.getOrMake1D(TString::Format("%sodtonewd_whyFail",sname.Data()),";TOF,Neutron, Other",4,-.5,3.5)->Fill(2);
					hists.getOrMake2D(TString::Format("%sodtonewd_otherFailLoc",sname.Data()),";local x position [cm] ;local y position [cm]",4000,-30,30,100,-50,50)->Fill(idigi->x(),idigi->y());
				}
			}



		}

	}

}

void ME0DigiAnalyzer::getDigiCuts(const ME0Geometry* mgeom,edm::Handle<std::vector<SimTrack>>&  simTrackH, edm::Handle<std::vector<PSimHit> >&  simHitH,
		const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis, const ME0DigiPreRecoMap& digiMap, TString sname) {

	for (const auto& simTrack : *simTrackH) {
		if(TMath::Abs(simTrack.type()) != 13) continue;
		const float pt = simTrack.momentum().pt();
		const float absETA = TMath::Abs(simTrack.momentum().eta());

		if(pt < 1) continue;
		if(absETA > 2.8 || absETA < 2.0 ) continue;

		TString ptName = sname;
		if(pt < 3 ) ptName += "ptleq3_";
		else if(pt < 5 ) ptName += "pteq3to5_";
		else if(pt < 20 ) ptName += "pteq5to20_";
		else ptName += "ptgeq20_";

		//collect simhits
		std::vector<const PSimHit* > simHits = ME0Helper::getMatchedSimHits(&simTrack,*simHitH);

		//Now get simhit segment info
		ME0Helper::ME0DigiList origMatchedDigis;
		std::vector<int> origDigiNewDigiIDX;
		ME0Helper::ME0DigiList newMatchedDigis = ME0Helper::getMatchedDigis(simHits,oldDigis,&newDigis,&digiMap,&origMatchedDigis,&origDigiNewDigiIDX);
		const int nHits = origMatchedDigis.size();
		bool inOneChamber = true;
		bool withinTimeBins = true;

		for(unsigned int iD = 0; iD < origMatchedDigis.size(); ++iD){
			if(iD && origMatchedDigis[iD].first.chamberId() != origMatchedDigis[0].first.chamberId() ) inOneChamber = false;
			if(origDigiNewDigiIDX[iD] >= 0) continue;
			withinTimeBins = false;
			hists.getOrMake1D(TString::Format("%s_missedDigi_tof",sname.Data()),";time of flight [ns]",400,-200,200)->Fill(origMatchedDigis[iD].second->tof());
		}
		int nPlotHits = nHits+2;
		if(!inOneChamber) nPlotHits = 0;
		else if(!withinTimeBins) nPlotHits = 1;
		else if(nHits > 6) nPlotHits = 9;

		if(!inOneChamber){
			std::map<unsigned int, unsigned int> chs;
			for(unsigned int iD = 0; iD < origMatchedDigis.size(); ++iD){
				chs[origMatchedDigis[iD].first.chamberId()]++;
			}
			unsigned int max = 0;
			for(auto it = chs.begin(); it != chs.end(); ++it){
				if(it->second > max) max = it->second;
			}
			hists.getOrMake1D(TString::Format("%s_nMaxHits",sname.Data()),";nHits",10,-0.5,9.5)->Fill(max);

		}
		hists.getOrMake1D(TString::Format("%s_muonClassifer",sname.Data()),";>1Chamber/Outside of time bin/ nHits (+2)",10,-0.5,9.5)->Fill(nPlotHits);
		hists.getOrMake1D(TString::Format("%s_muonClassifer",ptName.Data()),";>1Chamber/Outside of time bin/ nHits (+2)",10,-0.5,9.5)->Fill(nPlotHits);
		if(!inOneChamber || !withinTimeBins) continue;
		if(nHits < 3) continue;

		TString hitName = ptName;
		if(nHits == 3 ) hitName += "nHits3_";
		else if(nHits == 4 || nHits == 5 ) hitName += "nHits4to5_";
		else hitName += "nHits6_";

		std::vector<float> tofs(nHits,0);
		std::vector<float> etas(nHits,0);
		std::vector<float> phis(nHits,0);
		std::vector<float> nphis(nHits,0);
		std::vector<float> zs(nHits,0);
		std::vector<float> lys(nHits,0);
		std::vector<float> rs(nHits,0);
		std::vector<float> netas(nHits,0);
		int iSMin =-1;
		int iSMax =-1;
		float avgTOF = 0;
		for(unsigned int iD = 0; iD < origMatchedDigis.size(); ++iD){
			const ME0DetId refID=origMatchedDigis[iD].first;
			const ME0EtaPartition * part = mgeom->etaPartition(refID);
			const GlobalPoint gp = part->toGlobal(LocalPoint(origMatchedDigis[iD].second->x(),origMatchedDigis[iD].second->y(),0));
			tofs[iD] = newMatchedDigis[iD].second->tof();
			avgTOF += newMatchedDigis[iD].second->tof();
			etas[iD] = gp.eta();
			phis[iD] = TVector2::Phi_mpi_pi(gp.phi());
			zs[iD]   = gp.z();
			lys[iD] = origMatchedDigis[iD].second->y();
			rs[iD] = gp.perp();
//			std::cout <<refID<<" "<< gp.z()<<std::endl;
			if(iSMin == -1 || TMath::Abs(gp.z()) < TMath::Abs(zs[iSMin]) ) {
				iSMin = iD;
//				std::cout << "Min" <<std::endl;
			}
			if(iSMax == -1 || TMath::Abs(gp.z()) > TMath::Abs(zs[iSMax]) ) {
				iSMax = iD;
//				std::cout << "Max" <<std::endl;
			}
			const auto *nd =  newMatchedDigis[iD].second;
			netas[iD] = part->toGlobal(LocalPoint(nd->x(),nd->y(),0)).eta();
			nphis[iD] = part->toGlobal(LocalPoint(nd->x(),nd->y(),0)).phi();
		}
//		std::cout << etas[iSMin] <<" "<< etas[iSMax] <<" "<<etas[iSMin]-etas[iSMax]<<std::endl;
		avgTOF/=nHits;
		std::vector<float> sorted_tofs=  tofs; std::sort(sorted_tofs.begin(),sorted_tofs.end());
		std::vector<float> sorted_etas=  etas; std::sort(sorted_etas.begin(),sorted_etas.end());
		std::vector<float> sorted_phis=  phis; std::sort(sorted_phis.begin(),sorted_phis.end());
		std::vector<float> sorted_lys=  lys; std::sort(sorted_lys.begin(),sorted_lys.end());
		std::vector<float> sorted_rs=  rs; std::sort(sorted_rs.begin(),sorted_rs.end());
		std::vector<float> sorted_netas=  netas; std::sort(sorted_netas.begin(),sorted_netas.end());

		float avgETA = (netas[iSMax] + netas[iSMin])/2;
		float maxETAAway = -1;
		for( int iS = 0; iS < nHits; ++iS){
			if(iS == iSMax) continue;
			if(iS == iSMin) continue;
			float etaAway = TMath::Abs(netas[iS]) -TMath::Abs(avgETA);
			if(maxETAAway < 0 || TMath::Abs(etaAway) > TMath::Abs(maxETAAway)) maxETAAway = etaAway;
		}

		auto getMaxPhiAway = [&](const ME0Helper::ME0DigiList& digis) -> float{
			float maxPhiAway = -1;
			float m = (digis[iSMin].second->x() - digis[iSMax].second->x())/(zs[iSMin] - zs[iSMax]);
			float b = digis[iSMax].second->x() - zs[iSMax]*m;
			float m2 = (digis[iSMin].second->y() - digis[iSMax].second->y())/(zs[iSMin] - zs[iSMax]);
			float b2 = digis[iSMax].second->y() - zs[iSMax]*m2;
			for(int iS = 0; iS < nHits; ++iS){
				if(iS == iSMax) continue;
				if(iS == iSMin) continue;
				float extrapx = m*zs[iS] + b;
				float extrapy = m2*zs[iS] + b2;
				const ME0EtaPartition * part = mgeom->etaPartition(digis[iS].first);
				const GlobalPoint gp = part->toGlobal(LocalPoint(digis[iS].second->x(),digis[iS].second->y(),0));
				const GlobalPoint extgp = part->toGlobal(LocalPoint(extrapx,extrapy,0));
				float diffphi = TVector2::Phi_mpi_pi(extgp.phi() - gp.phi());
				if(maxPhiAway < 0 || std::fabs(diffphi) > std::fabs(maxPhiAway)) maxPhiAway = diffphi;
			}

			return maxPhiAway;
		};

		float maxOrigPhiAway = getMaxPhiAway(origMatchedDigis);
		float maxNewPhiAway = getMaxPhiAway(newMatchedDigis);

		auto getProjectedETAValues = [&](const ME0Helper::ME0DigiList& digis, float& seedDiff, float& maxCenDiff){
			const ME0EtaPartition * partMin = mgeom->etaPartition(digis[iSMin].first);
			const GlobalPoint gpMin = partMin->toGlobal(LocalPoint(digis[iSMin].second->x(),digis[iSMin].second->y(),0));
			const ME0EtaPartition * partMax = mgeom->etaPartition(digis[iSMax].first);
			const GlobalPoint gpMax = partMax->toGlobal(LocalPoint(digis[iSMax].second->x(),digis[iSMax].second->y(),0));

			double extPerp = gpMin.perp() * gpMax.z()/gpMin.z();
			double extX = extPerp*TMath::Cos(double(gpMin.phi()));
			double extY = extPerp*TMath::Sin(double(gpMin.phi()));
			GlobalPoint extGP(extX,extY,gpMax.z());
			seedDiff = TMath::Abs(extGP.eta()) - TMath::Abs(gpMax.eta());

			float m = (digis[iSMin].second->x() - digis[iSMax].second->x())/(zs[iSMin] - zs[iSMax]);
			float b = digis[iSMax].second->x() - zs[iSMax]*m;
			float m2 = (digis[iSMin].second->y() - digis[iSMax].second->y())/(zs[iSMin] - zs[iSMax]);
			float b2 = digis[iSMax].second->y() - zs[iSMax]*m2;
			maxCenDiff =-1;
			for(int iS = 0; iS < nHits; ++iS){
				if(iS == iSMax) continue;
				if(iS == iSMin) continue;
				float extrapx = m*zs[iS] + b;
				float extrapy = m2*zs[iS] + b2;
				const ME0EtaPartition * part = mgeom->etaPartition(digis[iS].first);
				const GlobalPoint gp = part->toGlobal(LocalPoint(digis[iS].second->x(),digis[iS].second->y(),0));
				const GlobalPoint extgp = part->toGlobal(LocalPoint(extrapx,extrapy,0));
				float diffeta = TMath::Abs(extgp.eta()) - TMath::Abs(gp.eta()) ;
				if(maxCenDiff < 0 || std::fabs(diffeta) > std::fabs(maxCenDiff)) maxCenDiff = diffeta;
			}


		};
		float seedETAFITDIFF=-1;
		float cenETAFITDIFF=-1;
		getProjectedETAValues(newMatchedDigis,seedETAFITDIFF,cenETAFITDIFF);

		//TOF plots
		hists.getOrMake1D(TString::Format("%s_trueProperties_avgTOF",hitName.Data()),";<tof> [ns]",400,-200,200)->Fill(avgTOF);
		hists.getOrMake1D(TString::Format("%s_trueProperties_seedTOFDiff",hitName.Data()),";|tof_{seed 1}-tof_{seed 2}| [ns]",200,0,200)->Fill(TMath::Abs(tofs[iSMin]-tofs[iSMax]));
		hists.getOrMake1D(TString::Format("%s_trueProperties_TOFDiff",hitName.Data()),";max[|tof_{i}-tof_{j}|] [ns]",200,0,200)->Fill(TMath::Abs(sorted_tofs.front()-sorted_tofs.back()));
		//eta plots
		hists.getOrMake1D(TString::Format("%s_trueProperties_seednETADiff",hitName.Data()),";|#eta_{seed 1}-#eta_{seed 2}| [ns]",400,-0.4,0.4)->Fill(    TMath::Abs(netas[iSMax])-TMath::Abs(netas[iSMin])     );
		hists.getOrMake1D(TString::Format("%s_trueProperties_nETADiff",hitName.Data()),";max[|#eta_{i}-#eta_{j}|] [ns]",400,-0.4,0.4)->Fill(sorted_netas.back() - sorted_netas.front());
		hists.getOrMake1D(TString::Format("%s_trueProperties_nETACen",hitName.Data()),";cen[|#eta_{i}-#eta_{j}|] [ns]",400,-0.4,0.4)->Fill(maxETAAway);

		hists.getOrMake1D(TString::Format("%s_trueProperties_ETASeedFitDiff",hitName.Data()),";seed fit[|#eta_{i}-#eta_{j}|] [ns]",400,-0.4,0.4)->Fill(seedETAFITDIFF);
		hists.getOrMake1D(TString::Format("%s_trueProperties_ETACenFitDiff",hitName.Data()),";cen fit[|#eta_{i}-#eta_{j}|] [ns]",400,-0.4,0.4)->Fill(cenETAFITDIFF);

		hists.getOrMake1D(TString::Format("%s_trueProperties_seedETADiff",hitName.Data()),";|#eta_{seed 1}-#eta_{seed 2}| [ns]",400,-0.4,0.4)->Fill(    etas[iSMax]-etas[iSMin]     );
		hists.getOrMake1D(TString::Format("%s_trueProperties_ETADiff",hitName.Data()),";max[|#eta_{i}-#eta_{j}|] [ns]",400,-0.4,0.4)->Fill(sorted_etas.back() - sorted_etas.front());
		hists.getOrMake1D(TString::Format("%s_trueProperties_seedYDiff",hitName.Data()),";|y_{seed 1}-y_{seed 2}| [ns]",500,-25,25)->Fill(    lys[iSMax]-lys[iSMin]     );
		hists.getOrMake1D(TString::Format("%s_trueProperties_YDiff",hitName.Data()),";max[|y_{i}-y_{j}|] [ns]",500,-25,25)->Fill(sorted_lys.back() - sorted_lys.front());
		hists.getOrMake1D(TString::Format("%s_trueProperties_seedRDiff",hitName.Data()),";|R_{seed 1}-R_{seed 2}| [ns]",500,-25,25)->Fill(    rs[iSMax]-rs[iSMin]     );
		hists.getOrMake1D(TString::Format("%s_trueProperties_RDiff",hitName.Data()),";max[|R_{i}-R_{j}|] [ns]",500,-25,25)->Fill(sorted_rs.back() - sorted_rs.front());
		//eta plots
		hists.getOrMake1D(TString::Format("%s_trueProperties_seedPHINDiff",hitName.Data()),";|#phi_{seed 1}-#phi_{seed 2}| [ns]",3500,0,.35)->Fill(TMath::Abs(TVector2::Phi_mpi_pi(nphis[iSMin]-nphis[iSMax])) );
		hists.getOrMake1D(TString::Format("%s_trueProperties_seedPHIDiff",hitName.Data()),";|#phi_{seed 1}-#phi_{seed 2}| [ns]",3500,0,.35)->Fill(TMath::Abs(TVector2::Phi_mpi_pi(phis[iSMin]-phis[iSMax])) );
		hists.getOrMake1D(TString::Format("%s_trueProperties_PHIDiff",hitName.Data()),";max[|#phi_{i}-#phi_{j}|] [ns]",3500,0,0.35)->Fill(TMath::Abs(TVector2::Phi_mpi_pi(sorted_phis.front()-sorted_phis.back())));

		hists.getOrMake1D(TString::Format("%s_trueProperties_PHICenDiff",hitName.Data()),";orig max[|#phi_{i}-#phi_{j}|] [ns]",3500,0,0.35)->Fill(std::fabs(maxOrigPhiAway));
		hists.getOrMake1D(TString::Format("%s_trueProperties_PHICenNDiff",hitName.Data()),";new max[|#phi_{i}-#phi_{j}|] [ns]",3500,0,0.35)->Fill(std::fabs(maxNewPhiAway));
	}


}

void ME0DigiAnalyzer::fitSegment(const ME0Geometry* mgeom,  edm::Handle<std::vector<SimTrack>>&  simTrackH, edm::Handle<std::vector<PSimHit> >&  simHitH,
		const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname){
	for (const auto& simTrack : *simTrackH) {
		if(TMath::Abs(simTrack.type()) != 13) continue;

		//collect simhits
		std::vector<const PSimHit* > simHits = ME0Helper::getMatchedSimHits(&simTrack,*simHitH);

		//Now get simhit segment info
		ME0Helper::ME0DigiList origMatchedDigis;
		ME0Helper::ME0DigiList newMatchedDigis = ME0Helper::getMatchedDigis(simHits,oldDigis,&newDigis,&digiMap,&origMatchedDigis);

		//Get the simtrack info...w/o any fitting or digitizing
		auto origSimHitProp = ME0Helper::getSimTrackProperties(mgeom,simHits);
		//basic filtering
		if(origSimHitProp.nLaysHit != 6 || !origSimHitProp.oneChamber) continue;

		//Fit the digis
		auto * origSegment = ME0Helper::buildSegment(mgeom,origMatchedDigis);
		auto * newSegment  = ME0Helper::buildSegment(mgeom,newMatchedDigis);
		//And get properties
		ME0Helper::SegmentProperties origSegmentProp =origSegment  ? ME0Helper::getSegmentProperties(mgeom->etaPartition(origSegment->recHits()[0]->rawId()),origSegment) : ME0Helper::SegmentProperties();
		ME0Helper::SegmentProperties newSegmentProp  =newSegment   ? ME0Helper::getSegmentProperties(mgeom->etaPartition(newSegment->recHits()[0]->rawId()),newSegment)  : ME0Helper::SegmentProperties();
//		ME0Helper::SegmentProperties origSegmentProp =ME0Helper::SegmentProperties();
//		ME0Helper::SegmentProperties newSegmentProp  =ME0Helper::SegmentProperties();

		//pt regions
		TString name = sname;
		if(simTrack.momentum().pt() < 3 ) name += "ptleq3_";
		else if(simTrack.momentum().pt() < 5 ) name += "pteq3to5_";
		else if(simTrack.momentum().pt() < 20 ) name += "pteq5to20_";
		else name += "ptgeq20_";
		if(origSegment == 0){
			std::cout <<"GGGG "<< origSimHitProp.nLaysHit <<" "<< origMatchedDigis.size() << std::endl;
		}
		if(origSegment)
		hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_chi2",name.Data()),";chi2",1000,0,100)
				->Fill(origSegment->chi2());
		if(origSegment)
		hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_chi2Prob",name.Data()),";chi2",1000,0,1)
				->Fill(TMath::Prob(origSegment->chi2(),2*6-4));

		if(origSegment && TMath::Abs(TVector2::Phi_mpi_pi(origSimHitProp.cenPhi - origSegmentProp.cenPhi)) > 0.01){
			hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_fail_chi2",name.Data()),";chi2",1000,0,100)
					->Fill(origSegment->chi2());
			hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_fail_chi2Prob",name.Data()),";chi2",1000,0,1)
					->Fill(TMath::Prob(origSegment->chi2(),2*6-4));

			cout << endl<< "("<<origSegmentProp.beginEta<<","<<origSegmentProp.beginPhi<<") ("<<origSegmentProp.cenEta<<","<<origSegmentProp.cenPhi
					<<") ("<<origSegmentProp.endEta<<","<<origSegmentProp.endPhi<<") "
					<<origSegmentProp.segAtCenter.point.x()<<","<<origSegmentProp.segAtCenter.point.y()<<" :: "<<origSegmentProp.initialPoint.x()<<","<<origSegmentProp.initialPoint.y()<<endl;
			for(unsigned int iH = 0; iH < origMatchedDigis.size(); ++iH){
//				cout <<"("<<origMatchedDigis[iH].first.layer()<<","<<mgeom->etaPartition(origMatchedDigis[iH].first)->surface().position()<<","<< mgeom->etaPartition(origMatchedDigis[iH].first)->surface().normalVector() <<","<<origMatchedDigis[iH].second->x()<<","<<origMatchedDigis[iH].second->y()<<") ";
				cout <<"("<<origMatchedDigis[iH].first.layer()<<","<<origMatchedDigis[iH].second->x()<<","<<origMatchedDigis[iH].second->y()<<") ";

			}


		}

		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_dphi",name.Data()),";#Delta#phi (SimHit)",1000,-1,1)
				->Fill(origSimHitProp.dPhi);

		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_origSeg_phidiff",name.Data()),";#phi position (SimHit - orig. seg)",100000,-.3,.3)
				->Fill(origSegment ? TVector2::Phi_mpi_pi(origSimHitProp.cenPhi - origSegmentProp.cenPhi) : 99.0 );
		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_origSeg_etadiff",name.Data()),";#eta position (SimHit - orig. seg)",1000,-1,1)
				->Fill(origSegment ? origSimHitProp.cenEta - origSegmentProp.cenEta : 99.0 );
		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_origSeg_dphidiff",name.Data()),";#Delta#phi (SimHit - orig. seg)",1000,-.003,.003)
				->Fill(origSegment ? origSimHitProp.dPhi - origSegmentProp.dPhi : 99.0 );
		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_origSeg_detadiff",name.Data()),";#Delta#eta  (SimHit - orig. seg)",1000,-1,1)
				->Fill(origSegment ? origSimHitProp.dEta - origSegmentProp.dEta : 99.0 );

		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_newSeg_phidiff",name.Data()),";#phi position (SimHit - new. seg)",1000,-.003,.003)
				->Fill(newSegment ? TVector2::Phi_mpi_pi(origSimHitProp.cenPhi - newSegmentProp.cenPhi) : 99.0 );
		hists.getOrMake1D(TString::Format("%ssegmentcheck_simHit_newSeg_dphidiff",name.Data()),";#Delta#phi (SimHit - new. seg)",1000,-.003,.003)
				->Fill(newSegment ? origSimHitProp.dPhi - newSegmentProp.dPhi : 99.0 );

		delete origSegment;
		delete newSegment;
	}

}
void
ME0DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	edm::Handle<ME0DigiPreRecoCollection> digisH;
	iEvent.getByToken(dToken_,digisH);

	edm::Handle<ME0DigiPreRecoCollection> digisNH;
	iEvent.getByToken(newDToken_,digisNH);

	edm::Handle<ME0DigiPreRecoMap> digisM;
	iEvent.getByToken(oldToNewMapToken_,digisM);

		edm::Handle <std::vector<SimTrack> > tracks;
	  iEvent.getByToken(track_token,tracks);

	  edm::Handle<std::vector<PSimHit> >  simHitH ;
	  iEvent.getByToken(shToken_,simHitH);

	edm::ESHandle<ME0Geometry> me0g;
	iSetup.get<MuonGeometryRecord>().get(me0g);
	const ME0Geometry* mgeom = &*me0g;

	hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);


	//First general info plots
	makeGeneralDigiPlots(mgeom,*digisH,TString::Format("%sstd_",runName.Data()));
	if(runNewDigis) makeGeneralDigiPlots(mgeom,*digisNH,TString::Format("%snew_",runName.Data()));
//	if(runNewDigis) compDigiPlots(mgeom,*digisH,*digisNH,*digisM,TString::Format("%s",runName.Data()));
//	if(runNewDigis) fitSegment(mgeom,tracks,simHitH,*digisH,*digisNH,*digisM,runName.Data());
	if(runNewDigis) getDigiCuts(mgeom,tracks,simHitH,*digisH,*digisNH,*digisM,runName.Data());



}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0DigiAnalyzer);
