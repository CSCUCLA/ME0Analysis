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


#include "../../AnalysisSupport/interface/HistGetter.h"
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

void ME0DigiAnalyzer::fitSegment(const ME0Geometry* mgeom,  edm::Handle<std::vector<SimTrack>>&  simTrackH, edm::Handle<std::vector<PSimHit> >&  simHitH,
		const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname){
	for (const auto& simTrack : *simTrackH) {
		if(TMath::Abs(simTrack.type()) != 13) continue;

		//collect simhits
		std::vector<const PSimHit* > simHits = ME0Helper::getMatchedSimHits(&simTrack,simHitH);

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
		ME0Helper::SegmentProperties origSegmentProp =origSegment  ? ME0Helper::getSegmentProperties(mgeom,origSegment) : ME0Helper::SegmentProperties();
		ME0Helper::SegmentProperties newSegmentProp  =newSegment   ? ME0Helper::getSegmentProperties(mgeom,newSegment)  : ME0Helper::SegmentProperties();

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

//		if(origSegment && TMath::Abs(TVector2::Phi_mpi_pi(origSimHitProp.cenPhi - origSegmentProp.cenPhi)) > 0.01){
//			hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_fail_chi2",name.Data()),";chi2",1000,0,100)
//					->Fill(origSegment->chi2());
//			hists.getOrMake1D(TString::Format("%ssegmentcheck_origSeg_fail_chi2Prob",name.Data()),";chi2",1000,0,1)
//					->Fill(TMath::Prob(origSegment->chi2(),2*6-4));
//
//			cout << endl<< "("<<origSegmentProp.beginEta<<","<<origSegmentProp.beginPhi<<") ("<<origSegmentProp.cenEta<<","<<origSegmentProp.cenPhi
//					<<") ("<<origSegmentProp.endEta<<","<<origSegmentProp.endPhi<<") "
//					<<origSegmentProp.segAtCenter.point.x()<<","<<origSegmentProp.segAtCenter.point.y()<<" :: "<<origSegmentProp.initialPoint.x()<<","<<origSegmentProp.initialPoint.y()<<endl;
//			for(unsigned int iH = 0; iH < origMatchedDigis.size(); ++iH){
//				cout <<"("<<origMatchedDigis[iH].first.layer()<<","<<origMatchedDigis[iH].second->x()<<","<<origMatchedDigis[iH].second->y()<<") ";
//			}
//
//
//		}

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
	if(runNewDigis) compDigiPlots(mgeom,*digisH,*digisNH,*digisM,TString::Format("%s",runName.Data()));
	if(runNewDigis) fitSegment(mgeom,tracks,simHitH,*digisH,*digisNH,*digisM,runName.Data());





}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0DigiAnalyzer);
