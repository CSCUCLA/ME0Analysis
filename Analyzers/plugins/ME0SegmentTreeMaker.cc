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
#include "AnalysisSupport/TreeInterface/interface/TreeWriter.h"

#include "AnalysisSupport/Utilities/interface/Types.h"
#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../interface/ME0Helper.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include <TVector2.h>
#include <TMath.h>
#include <TRandom3.h>

using namespace std;
using ASTypes::size;

class ME0SegmentTreeMaker : public edm::EDAnalyzer {
public:
	explicit ME0SegmentTreeMaker(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
	  runName(iConfig.getUntrackedParameter<std::string>("runName")), 	  tree (outFileName,"Events","")

	{
		dToken_        = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
		  newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  oldToNewMapToken_  = consumes<ME0DigiPreRecoMap>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
		  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
		  segToken_    = consumes<ME0SegmentCollection>( edm::InputTag(iConfig.getParameter<std::string>("segmentCollection")) );
		  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag(iConfig.getParameter<std::string>("recHitCollection")) );
		  tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );

		  simMuon_pt                 = tree.addMulti<float>  (     "simMuon_pt"                      ,0);
		  simMuon_eta                = tree.addMulti<float>  (     "simMuon_eta"                     ,0);
		  simMuon_phi                = tree.addMulti<float>  (     "simMuon_phi"                     ,0);
		  simMuon_q                  = tree.addMulti<int>    (     "simMuon_q"                       ,0);
		  simMuon_gen_eta            = tree.addMulti<float>  (     "simMuon_gen_eta"             ,0);
		  simMuon_gen_phi            = tree.addMulti<float>  (     "simMuon_gen_phi"             ,0);
		  simMuon_gen_dphi           = tree.addMulti<float>  (     "simMuon_gen_dphi"            ,0);
		  simMuon_gen_deta           = tree.addMulti<float>  (     "simMuon_gen_deta"            ,0);
		  simMuon_gen_x              = tree.addMulti<float>  (     "simMuon_gen_x"               ,0);
		  simMuon_gen_y              = tree.addMulti<float>  (     "simMuon_gen_y"               ,0);
		  simMuon_gen_dx             = tree.addMulti<float>  (     "simMuon_gen_dx"              ,0);
		  simMuon_gen_dy             = tree.addMulti<float>  (     "simMuon_gen_dy"              ,0);
		  simMuon_segmentIDX         = tree.addMulti<int>    (     "simMuon_segmentIDX"              ,0);
		  simMuon_segment_quality    = tree.addMulti<int>    (     "simMuon_segment_quality"              ,0);
		  simMuon_segment_nGoodHits  = tree.addMulti<int>    (     "simMuon_segment_nGoodHits"              ,0);
		  simMuon_segment_nBadHits   = tree.addMulti<int>    (     "simMuon_segment_nBadHits"              ,0);
		  segment_eta                = tree.addMulti<float>  (     "segment_eta"              ,0);
		  segment_phi                = tree.addMulti<float>  (     "segment_phi"              ,0);
		  segment_dphi               = tree.addMulti<float>  (     "segment_dphi"             ,0);
		  segment_deta               = tree.addMulti<float>  (     "segment_deta"             ,0);
		  segment_x                  = tree.addMulti<float>  (     "segment_x"                ,0);
		  segment_y                  = tree.addMulti<float>  (     "segment_y"                ,0);
		  segment_dx                 = tree.addMulti<float>  (     "segment_dx"               ,0);
		  segment_dy                 = tree.addMulti<float>  (     "segment_dy"               ,0);


		  tree.book();


	}


	~ME0SegmentTreeMaker() {
		  tree.write();

	}



private:
	virtual void beginJob() {};
	virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

	void fillTree(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, const ME0RecHitCollection& recHits, const std::vector<SimTrack>&  simTrackH, const std::vector<PSimHit>&  simHitH,
			const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap,
			const edm::Handle<TrackingParticleCollection>& trParticles){

		auto simMuons = ME0Helper::fillSimMuonsByTP(trParticles,simTrackH,simHitH);
		ME0Helper::DigiInfoMap digiInfo;
		ME0Helper::fillDigiInfoMap(newDigis,digiMap,oldDigis,simHitH,digiInfo);
		ME0Helper::associateSimMuons(simMuons,digiInfo);
		ME0Helper::SegmentCatMap segmentCats;
		ME0Helper::fillSegmentCategories(segments,recHits,digiInfo,simMuons,segmentCats);

		tree.reset();

		  std::map<int,int> indxMap;
		  int idx = 0;
		  int cleanedidx = 0;

		  for(auto iC = segments.id_begin(); iC != segments.id_end(); ++iC){
			  auto ch_segs = segments.get(*iC);
			  const auto * chamber = mgeom->chamber(*iC);
			  for(auto iS = ch_segs.first; iS != ch_segs.second; ++iS){
				  const bool passTime = std::fabs(iS->time()) <11.0;
				  if(!passTime){
					  indxMap[idx] = -1; idx++; continue;
				  }
				  indxMap[idx] = cleanedidx;
				  idx++; cleanedidx++;

				  GlobalPoint segGlobPos(chamber->toGlobal(iS->localPosition()));
				  tree.fillMulti(segment_eta   ,float(segGlobPos.eta()));
				  tree.fillMulti(segment_phi   ,float(segGlobPos.phi()));

				  float segDPhi, segDEta;
				  auto extrap = [&] (const LocalPoint& point, const LocalVector& dir, double extZ) -> LocalPoint {
					  double extX = point.x()+extZ*dir.x()/dir.z();
					  double extY = point.y()+extZ*dir.y()/dir.z();
					  return LocalPoint(extX,extY,extZ);
				  };
				  auto getDPhiDEta= [&](const LocalPoint& point, const LocalVector& dir, float& dPhi, float& dEta)  {
					  float chamberZ = chamber->position().z();
					  LocalPoint projHigh = extrap(point,dir, ME0Helper::endOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
					  LocalPoint projLow = extrap(point,dir, ME0Helper::beginOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
					  auto globLow  = chamber->toGlobal(projLow );
					  auto globHigh = chamber->toGlobal(projHigh);
					  dPhi =  TVector2::Phi_mpi_pi(globHigh.phi() - globLow.phi());
					  dEta =  globHigh.eta() - globLow.eta();
				  };
				  getDPhiDEta(iS->localPosition(),iS->localDirection(),segDPhi,segDEta);
				  tree.fillMulti(segment_dphi  ,float(segDPhi));
				  tree.fillMulti(segment_deta  ,float(segDEta));
				  tree.fillMulti(segment_x     ,float(iS->localPosition().x()));
				  tree.fillMulti(segment_y     ,float(iS->localPosition().y()));
				  tree.fillMulti(segment_dx    ,float(iS->localDirection().x()/iS->localDirection().z()));
				  tree.fillMulti(segment_dy    ,float(iS->localDirection().y()/iS->localDirection().z()));


			  }
		  }


		  for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
				const auto& muon = simMuons[iM];
				const float absETA = TMath::Abs(muon.trPart->eta());
				if(absETA > 2.8 || absETA < 2.0 ) continue;

			    tree.fillMulti(simMuon_pt             ,float(muon.trPart->pt()));
			    tree.fillMulti(simMuon_eta            ,float(muon.trPart->eta()));
			    tree.fillMulti(simMuon_phi            ,float(muon.trPart->phi()));
			    tree.fillMulti(simMuon_q              ,int(muon.trPart->charge()));

				int typeC = 3; //Perfect/MissOne/MissMultiple/Lost
				int nGoodHits = 0;
				int nTotalHits = 0;
				if(muon.segments.size()) {
					nTotalHits = (segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second)->nRecHits();
					nGoodHits = muon.segments[0].second;
					if(segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first == ME0Helper::MUON_COMP_PURE) typeC = 0;
					else if(nTotalHits - nGoodHits <= 1 ) typeC = 1;
					else if( nGoodHits > (nTotalHits - nGoodHits)) typeC = 2;
					//Time cut
					const ME0Segment * segment = &*(segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second);
					if(std::fabs(segment->time()) >= 11.0) typeC = 3;
				}


			    tree.fillMulti(simMuon_segment_quality ,  typeC);//Perfect / MissOne/ MissSome /Lost
			    tree.fillMulti(simMuon_segment_nGoodHits ,  nGoodHits);
			    tree.fillMulti(simMuon_segment_nBadHits ,  nTotalHits - nGoodHits);
			    int segIDX = -1;

			    if(typeC <= 2){
			    	segIDX = int(segments.get(muon.segments[0].first.first).first - segments.begin()) +  muon.segments[0].first.second;
			    }
			    tree.fillMulti(simMuon_segmentIDX     ,  segIDX >= 0 ? indxMap[segIDX] : -1);

			    auto genProp = ME0Helper::getSimTrackProperties(mgeom,muon.simHits);
		    	tree.fillMulti(simMuon_gen_eta       ,float(genProp.cenEta));
		    	tree.fillMulti(simMuon_gen_phi       ,float(genProp.cenPhi));
		    	tree.fillMulti(simMuon_gen_dphi      ,float(genProp.dPhi));
		    	tree.fillMulti(simMuon_gen_deta      ,float(genProp.dEta));
		    	tree.fillMulti(simMuon_gen_x         ,float(genProp.theOrigin.x()));
		    	tree.fillMulti(simMuon_gen_y         ,float(genProp.theOrigin.y()));
		    	tree.fillMulti(simMuon_gen_dx        ,genProp.nLaysHit > 1 ? float(genProp.theLocalDirection.x() / genProp.theLocalDirection.z() )  : float(0.0));
		    	tree.fillMulti(simMuon_gen_dy        ,genProp.nLaysHit > 1 ? float(genProp.theLocalDirection.y() / genProp.theLocalDirection.z() )  : float(0.0));


		  }




		  tree.fillTree();
	}

private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoMap>                 oldToNewMapToken_  ;
    edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
    edm::EDGetTokenT<ME0SegmentCollection>          segToken_ ;
    edm::EDGetTokenT<std::vector<SimTrack>> track_token;
    edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;
    edm::EDGetTokenT<TrackingParticleCollection>          tpToken_ ;


	TString outFileName;
	TString runName;


    TreeWriter tree;
    size simMuon_pt                    =-1;
    size simMuon_eta                   =-1;
    size simMuon_phi                   =-1;
    size simMuon_q                     =-1;
    size simMuon_segmentIDX            =-1;

    size simMuon_gen_eta                   =-1;
    size simMuon_gen_phi                   =-1;
    size simMuon_gen_dphi                  =-1;
    size simMuon_gen_deta                  =-1;
    size simMuon_gen_x                     =-1;
    size simMuon_gen_y                     =-1;
    size simMuon_gen_dx                    =-1;
    size simMuon_gen_dy                    =-1;

    size simMuon_segment_quality               =-1;
    size simMuon_segment_nGoodHits             =-1;
    size simMuon_segment_nBadHits              =-1;
    size segment_eta                   =-1;
    size segment_phi                   =-1;
    size segment_dphi                  =-1;
    size segment_deta                  =-1;
    size segment_x                     =-1;
    size segment_y                     =-1;
    size segment_dx                    =-1;
    size segment_dy                    =-1;


};




void
ME0SegmentTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

	  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
	  iEvent.getByToken(tpToken_,TPCollectionH);


	  fillTree(mgeom,*segH,*rechitsH,*tracks,*simHitH,*digisH,*digisNH,*digisM,TPCollectionH);
}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0SegmentTreeMaker);
