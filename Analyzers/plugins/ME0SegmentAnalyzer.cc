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
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class ME0SegmentAnalyzer : public edm::EDAnalyzer {
public:
	explicit ME0SegmentAnalyzer(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
	  runName(iConfig.getUntrackedParameter<std::string>("runName"))
	{
		dToken_        = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
		  newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  oldToNewMapToken_  = consumes<ME0DigiPreRecoMap>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
		  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
		  segToken_    = consumes<ME0SegmentCollection>( edm::InputTag(iConfig.getParameter<std::string>("segmentCollection")) );
		  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag(iConfig.getParameter<std::string>("recHitCollection")) );

	}

	~ME0SegmentAnalyzer() {
		hists.write(outFileName);
	}



private:
	virtual void beginJob() {};
	virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	void makeSegmentPlots(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, TString name){
		hists.getOrMake1D(TString::Format("%snSegments",name.Data()),";# of segments per event",200,-0.5,190.5)->Fill(segments.size());
		for(const auto& seg : segments ){
			const auto* ch = mgeom->chamber(seg.me0DetId().chamberId());
			GlobalPoint seggp = ch->toGlobal(seg.localPosition());
			LocalPoint seglp = seg.localPosition();
			const auto params = seg.parameters();
			hists.getOrMake1D(TString::Format("%s_segment_eta",name.Data()),";segment #eta",500,-5,5)->Fill(seggp.eta());
			hists.getOrMake1D(TString::Format("%s_segment_phi",name.Data()),";segment #phi",500,-3.2,3.2)->Fill(TVector2::Phi_mpi_pi(seggp.phi()));
			hists.getOrMake1D(TString::Format("%s_segment_nHits",name.Data()),";segment # of hits",10,-0.5,9.5)->Fill(seg.recHits().size());
			hists.getOrMake1D(TString::Format("%s_segment_chi2",name.Data()),";segment chi2",500,0,200)->Fill(seg.chi2());
			hists.getOrMake1D(TString::Format("%s_segment_chi2oDeg",name.Data()),";segment chi2/# of degrees of freedom",500,0,200)->Fill(seg.chi2()/seg.degreesOfFreedom());
			hists.getOrMake1D(TString::Format("%s_segment_chi2Prob",name.Data()),";segment chi2 prob",500,0,1)->Fill(TMath::Prob(seg.chi2(),seg.degreesOfFreedom()));
			hists.getOrMake1D(TString::Format("%s_segment_locx",name.Data()),";segment local x",500,-40,40)->Fill(seglp.x());
			hists.getOrMake1D(TString::Format("%s_segment_locy",name.Data()),";segment local y",500,-50,50)->Fill(seglp.y());
			hists.getOrMake1D(TString::Format("%s_segment_locz",name.Data()),";segment local z",500,-50,50)->Fill(seglp.z());
			hists.getOrMake1D(TString::Format("%s_segment_dxdz",name.Data()),";segment dx/dz",500,-5,5)->Fill(params[0]);
			hists.getOrMake1D(TString::Format("%s_segment_dydz",name.Data()),";segment dx/dz",500,-5,5)->Fill(params[1]);
			hists.getOrMake1D(TString::Format("%s_segment_tof",name.Data()),";segment tof",400,-200,200)->Fill(seg.time());
		}

	}


	std::string getDigiString(const ME0Helper::DigiInfo& digi) {
		std::string out;
		if(digi.muonIDXs.size() ) out = "M(" + std::to_string(digi.muonIDXs[0]) + ")";
		else if(digi.simHits.size() ) out = "S(" + std::to_string(digi.simHits[0]->particleType())+")";
		else {
			bool nonPrompt = true;
			for(unsigned int iD = 0; iD < digi.origDigis.size(); ++iD){
				if(digi.origDigis[iD]->prompt()){
					nonPrompt = false;
					out = "P("+std::to_string(digi.origDigis[iD]->pdgid())+")";
					break;
				}
			}
			if(nonPrompt) out = "N("+std::to_string(digi.digi->pdgid())+")";
		}
		return out;
	}
	void printSegmentInfo(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, const ME0RecHitCollection& recHits, const std::vector<SimTrack>&  simTrackH, const std::vector<PSimHit>&  simHitH,
			const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap){
		auto simMuons = ME0Helper::fillSimMuons(simTrackH,simHitH);
		ME0Helper::DigiInfoMap digiInfo;
		ME0Helper::fillDigiInfoMap(newDigis,digiMap,oldDigis,simHitH,digiInfo);
		ME0Helper::associateSimMuons(simMuons,digiInfo);
		ME0Helper::SegmentCatMap segmentCats;
		ME0Helper::fillSegmentCategories(segments,recHits,digiInfo,simMuons,segmentCats);

		std::cout << "<<<<<<<<<<<< ------  NEW EVENT  ------ >>>>>>>>>>>>>>" <<std::endl;
//		for(auto it : oldDigis){
//			std::cout << it.first <<" ";
//			for(auto id = it.second.first; id != it.second.second; ++id) cout << "or("<<id->x() <<","<<id->y()<<") "<< id->tof() <<" ";
//			std::cout<<std::endl;
//		}


		for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
			const auto* muon = simMuons[iM].track;
			std::cout << "SimMuon "<< iM <<" : ("<< muon->momentum().pt() <<","<<muon->momentum().eta()<<","<<muon->momentum().phi()<<") "
					<< simMuons[iM].simHits.size() <<" simHits "<<  simMuons[iM].digiInfos.size() <<" digis ";
			for(unsigned int iMS = 0; iMS < simMuons[iM].segments.size(); ++iMS)
				std::cout <<" S "<<simMuons[iM].segments[iMS].first.first <<" "<< simMuons[iM].segments[iMS].first.second <<" "<<simMuons[iM].segments[iMS].second <<" ";
			std::cout << std::endl;;
			std::map<ME0DetId,std::pair<std::vector<const PSimHit*>,std::vector<const ME0Helper::DigiInfo*>>> chs;
			for(const auto&  as : simMuons[iM].digiInfos ) { chs[as.first].second.push_back(&digiInfo[as.first][as.second]);  }
			for(const auto*  as : simMuons[iM].simHits ) { chs[as->detUnitId()].first.push_back(as);  }

			for(const auto& as : chs) {std::cout << "("<< as.first <<") : ";
			for(unsigned int iS = 0; iS < as.second.first.size(); ++iS) std::cout << "sh"<<as.second.first[iS]->entryPoint()<<" "<<as.second.first[iS]->tof() <<" " ;
			for(unsigned int iS = 0; iS < as.second.second.size(); ++iS) std::cout << "d("<<as.second.second[iS]->digi->x()<<","<<as.second.second[iS]->digi->y()<<") ";
//			auto od = oldDigis.get(as.first);
//			for(auto id = od.first; id != od.second; ++id) if(id->prompt()) cout << "or("<<id->x() <<","<<id->y()<<") "<< id->tof() <<" ";

//			const auto& dc = digiInfo[as.first];
//			for(const auto& id : dc) {
//				std::cout <<"NND("<<id.digi->x()<<","<<id.digi->y() <<","<<id.digi->tof() <<") " << id.muonIDXs.size() <<" ";
//				for(const auto* od : id.origDigis) std::cout <<" o("<<od->x()<<","<<od->y()<<","<< od->tof() <<") " << od->pdgid() <<" ";
//				for(const auto* od : id.simHits) std::cout <<" s("<<od->entryPoint()<<","<< od->tof() <<") " << od->particleType() <<" ";
//			}


			std::cout<< std::endl;

			}
		}
		std::cout << segments.size() <<" segments in "<< segments.id_size() <<"vchambers"<< std::endl;
		for(auto iC = segments.id_begin(); iC != segments.id_end(); ++iC){
			auto ch_segs = segments.get(*iC);
			std::cout << *iC <<" with "<< int(ch_segs.second - ch_segs.first) <<" segments"<< std::endl;
			for(auto iS = ch_segs.first; iS != ch_segs.second; ++iS){
				auto params = iS->parameters();
				auto glb = mgeom->chamber(*iC)->toGlobal(iS->localPosition());
				auto glbdr = mgeom->chamber(*iC)->toGlobal(iS->localDirection());
				std:: cout << "pos: "<<  iS->localPosition() <<" err: "<< iS->localPositionError() << "(dx/dz,dy/dz): ("<< params[0] <<","<< params[1]<<") "
						<< "glb (eta,phi): ("<< glb.eta()<<","<<glb.phi()<<")" << " glbdir (eta,phi): ("<< glbdr.eta()<<","<<glbdr.phi()<<") "
						<< "chi2: "<< iS->chi2() <<" chi2/DOF " << iS->chi2()/iS->degreesOfFreedom() << " cat: "<< ME0Helper::SegmentGenCatNames[segmentCats[*iC][iS - ch_segs.first].first]
						<<" "<<segmentCats[*iC][iS - ch_segs.first].second<<std::endl;




				  const auto& rh = iS->specificRecHits();
				  for(unsigned int iR = 0; iR < rh.size(); ++iR){
					  int idx = ME0Helper::getRecHitIndex(recHits,rh[iR]);
					  std::string type = idx >= 0 ? getDigiString(digiInfo[rh[iR].me0Id()][idx] ) : "!(-1)";
//					  if(idx >= 0)
//						  std::cout << rh[iR].me0Id() <<" "<< rh[iR].localPosition()<<" ("<<digiInfo[rh[iR].me0Id()][idx].digi->x() <<","<<digiInfo[rh[iR].me0Id()][idx].digi->y()<<") "<< type << std::endl;
					  std::cout << rh[iR].me0Id() <<" "<< rh[iR].localPosition()<< " "<< rh[iR].tof()<<" " << type << std::endl;

//						std::cout <<"GRRR"<<std::endl;
//						auto digis = newDigis.get(rh[iR].me0Id() );
//						auto rhs = recHits.get(rh[iR].me0Id() );
//						auto& di = digiInfo[rh[iR].me0Id() ];
//						for(unsigned int iH = 0; iH < di.size(); ++iH){
//							cout << (rhs.first + iH)->localPosition() <<" ("<<(digis.first + iH)->x() <<","<<(digis.first + iH)->y()<<") (" <<di[iH].digi->x() <<","<<di[iH].digi->y()<<std::endl;
//						}
//						std::cout <<"GRRRENNNDDDD"<<std::endl;


				  }

			}

		}
	}

	void segmentEfficiencyPlots(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, const ME0RecHitCollection& recHits, const std::vector<SimTrack>&  simTrackH, const std::vector<PSimHit>&  simHitH,
			const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap, TString sname){
		auto simMuons = ME0Helper::fillSimMuons(simTrackH,simHitH);
		ME0Helper::DigiInfoMap digiInfo;
		ME0Helper::fillDigiInfoMap(newDigis,digiMap,oldDigis,simHitH,digiInfo);
		ME0Helper::associateSimMuons(simMuons,digiInfo);
		ME0Helper::SegmentCatMap segmentCats;
		ME0Helper::fillSegmentCategories(segments,recHits,digiInfo,simMuons,segmentCats);
		int nGM = 0;


		for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
			const auto& muon = simMuons[iM];
			const float pt = muon.track->momentum().pt();
			const float absETA = TMath::Abs(muon.track->momentum().eta());

			if(pt < 1) continue;
			if(absETA > 2.8 || absETA < 2.0 ) continue;
			std::map<ME0DetId,int> chMap;
			for(const auto& di :muon.simHits ) chMap[ME0DetId(di->detUnitId()).chamberId()]++;
			if(chMap.size() > 1) continue;
			nGM++;
			std::string ptstr = "pteq1to3_";
			if(pt >= 3 && pt < 5 ) ptstr  = "pteq3to5_";
			else if(pt >= 5 && pt < 20 ) ptstr = "pteq5to20_";
			else  if(pt >= 20) ptstr = "ptgeq20_";

//			std::string catstr = "MISSING_SEGMENT";
//			if(muon.segments.size()) catstr = ME0Helper::SegmentGenCatNames[segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first];
			int typeC = 3; //Perfect/MissOne/MissMultiple/Lost
			if(muon.segments.size()) {
				int nRD = (segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second)->nRecHits();
				int nSD = muon.segments[0].second;
				if(segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first == ME0Helper::MUON_COMP_PURE) typeC = 0;
				else if(nRD - nSD <= 1 ) typeC = 1;
				else if( nSD > (nRD - nSD)) typeC = 2;
			}
			std::string catstr = "lost";
			std::string shortCatstr = "lost";
			std::string compstr = muon.segments.size() ?  ME0Helper::SegmentGenCatNames[segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first] : "MISSING_SEGMENT";
			switch(typeC){
			case 0 :
				catstr = "perfect";
				shortCatstr = "perfect";
				break;
			case 1:
				catstr = "oneFakeHit";
				shortCatstr = "dirty";
				break;
			case 2:
				catstr = "twoFakeHits";
				shortCatstr = "dirty";
				break;
			default:
				catstr = "lost";
				shortCatstr = "lost";
				break;
			}


//
//			if(catstr != "perfect"){
//				int nRD = muon.segments.size() ? (segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second)->nRecHits() : 0;
//				int nSD = muon.segments.size() ? muon.segments[0].second : 0;
//				std::cout << nRD <<" "<< nSD <<" "<< catstr <<std::endl;
//				printSegmentInfo(mgeom,segments,recHits,simTrackH,simHitH,oldDigis,newDigis,digiMap);
//			}


			hists.getOrMake1D(TString::Format("%sall_muon_pt",sname.Data()),";muon p_{T}",60,0,30)->Fill(pt);
			hists.getOrMake1D(TString::Format("%s%s_muon_pt",sname.Data(),catstr.c_str()),";muon p_{T}",60,0,30)->Fill(pt);


			auto doShortTypePlots = [&] (TString prefix){
				const ME0Segment * segment =  muon.segments.size() ?  &*(segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second) : 0;
				int nHits = segment ? segment->nRecHits() : 0 ;
				hists.getOrMake1D(TString::Format("%s_nHits",prefix.Data()),";# of hits",7,-0.5,6.5)->Fill(nHits);

				if(!segment) return;
				float chi2 = segment->chi2();
				int ndof = segment->degreesOfFreedom();
				hists.getOrMake1D(TString::Format("%s_chi2",prefix.Data()),";#chi^{2}",200,0,100)->Fill(chi2);
				hists.getOrMake1D(TString::Format("%s_chi2ondof",prefix.Data()),";#chi^{2}/dof",200,0,100)->Fill(chi2/ndof);
				auto prop = ME0Helper::getSegmentProperties(mgeom->chamber(muon.segments[0].first.first),&*segment);
				hists.getOrMake1D(TString::Format("%s_dPhi",prefix.Data()),";#Delta#phi",200,0,.1)->Fill(std::fabs(prop.dPhi));

				auto stp = ME0Helper::getSimTrackProperties(mgeom,muon.simHits);
				hists.getOrMake1D(TString::Format("%s_simMreco_phi",prefix.Data()),";#phi_{sim} - #phi_{seg}",200,-0.01,0.01)->Fill(stp.cenPhi - prop.cenPhi);
				hists.getOrMake1D(TString::Format("%s_simMreco_eta",prefix.Data()),";#eta_{sim} - #eta_{seg}",200,-0.25,0.25)->Fill(stp.cenEta - prop.cenEta);
				hists.getOrMake1D(TString::Format("%s_simMreco_dphi",prefix.Data()),";#Delta#phi_{sim} - #Delta#phi_{seg}",200,-0.01,0.01)->Fill(stp.dPhi - prop.dPhi);
				hists.getOrMake1D(TString::Format("%s_simMreco_deta",prefix.Data()),";#Delta#eta_{sim} - #Delta#eta_{seg}",200,-0.25,0.25)->Fill(stp.dEta - prop.dEta);






			};
			doShortTypePlots(TString::Format("%sall_muon_seg",sname.Data()));
			doShortTypePlots(TString::Format("%s%s_muon_seg",sname.Data(),shortCatstr.c_str()));
			doShortTypePlots(TString::Format("%s%s_muon_seg",sname.Data(),(ptstr + shortCatstr).c_str()));

			doShortTypePlots(TString::Format("%s%s_muon_seg",sname.Data(),compstr.c_str()));
			doShortTypePlots(TString::Format("%s%s_muon_seg",sname.Data(),(ptstr + compstr).c_str()));


		}

		if(nGM == 2){

			int nFakes = 0;
			for(auto iC = segments.id_begin(); iC != segments.id_end(); ++iC){
				auto ch_segs = segments.get(*iC);
				auto& ch_cats = segmentCats[*iC];
				for(auto iS = ch_segs.first; iS != ch_segs.second; ++iS){
					auto cat =ch_cats[iS-ch_segs.first].first;
					if(cat <= ME0Helper::MUON_MISS_DIRTY_NEUT) continue;
					nFakes++;
					int nHits = iS->nRecHits();
					float chi2 = iS->chi2();
					int ndof = iS->degreesOfFreedom();
					std::string catstr = ME0Helper::SegmentGenCatNames[cat];
					hists.getOrMake1D(TString::Format("%sall_fake_seg_nHits",sname.Data()),";# of hits",7,-0.5,6.5)->Fill(nHits);
					hists.getOrMake1D(TString::Format("%s%s_fake_seg_nHits",sname.Data(),catstr.c_str()),";# of hits",7,-0.5,6.5)->Fill(nHits);

					hists.getOrMake1D(TString::Format("%sall_fake_seg_chi2",sname.Data()),";#chi^{2}",200,0,100)->Fill(chi2);
					hists.getOrMake1D(TString::Format("%s%s_fake_seg_chi2",sname.Data(),catstr.c_str()),";#chi^{2}",200,0,100)->Fill(chi2);

					hists.getOrMake1D(TString::Format("%sall_fake_seg_chi2ondof",sname.Data()),";#chi^{2}",200,0,100)->Fill(chi2/ndof);
					hists.getOrMake1D(TString::Format("%s%s_fake_seg_chi2ondof",sname.Data(),catstr.c_str()),";#chi^{2}",200,0,100)->Fill(chi2/ndof);

					auto prop = ME0Helper::getSegmentProperties(mgeom->chamber(*iC),&*iS);
					hists.getOrMake1D(TString::Format("%s%s_fake_seg_dPhi",sname.Data(),(catstr).c_str()),";#Delta#phi",50,0,.025)->Fill(std::fabs(prop.dPhi));
					hists.getOrMake1D(TString::Format("%sall_fake_seg_dPhi",sname.Data()),";#Delta#phi",50,0,.025)->Fill(std::fabs(prop.dPhi));

				}
			}
			hists.getOrMake1D(TString::Format("%sfake_nSegs",sname.Data()),";# of segments",100,-0.5,99.5)->Fill(nFakes);

		}
	}



private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoMap>                 oldToNewMapToken_  ;
    edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
    edm::EDGetTokenT<ME0SegmentCollection>          segToken_ ;
    edm::EDGetTokenT<std::vector<SimTrack>> track_token;
    edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;


	TString outFileName;
	TString runName;
	HistGetter hists;
};




void
ME0SegmentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

	  edm::Handle<ME0SegmentCollection >  segH ;
	  iEvent.getByToken(segToken_,segH);

	  edm::Handle<ME0RecHitCollection> rechitsH;
	  iEvent.getByToken(rhToken_,rechitsH);

	edm::ESHandle<ME0Geometry> me0g;
	iSetup.get<MuonGeometryRecord>().get(me0g);
	const ME0Geometry* mgeom = &*me0g;

	hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);


	//First general info plots
	makeSegmentPlots(mgeom,*segH,TString::Format("%sgen_",runName.Data()));
//	printSegmentInfo(mgeom,*segH,*rechitsH,*tracks,*simHitH,*digisH,*digisNH,*digisM);
	segmentEfficiencyPlots(mgeom,*segH,*rechitsH,*tracks,*simHitH,*digisH,*digisNH,*digisM,TString::Format("%s_",runName.Data()));

}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0SegmentAnalyzer);
