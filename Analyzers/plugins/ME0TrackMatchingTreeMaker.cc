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
#include "AnalysisSupport/Utilities/interface/Types.h"
#include "AnalysisSupport/TreeInterface/interface/TreeWriter.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../interface/ME0Helper.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include <DataFormats/MuonReco/interface/ME0Muon.h>

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class ME0TrackMatchingTreeMaker : public edm::EDAnalyzer {
public:
	explicit ME0TrackMatchingTreeMaker(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
	  runName(iConfig.getUntrackedParameter<std::string>("runName")),
	  tree (outFileName,"Events","")
	{
		dToken_        = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
		  newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  oldToNewMapToken_  = consumes<ME0DigiPreRecoMap>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
		  vert_token = consumes<std::vector<SimVertex>>( edm::InputTag("g4SimHits","","SIM") );
		  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
		  segToken_    = consumes<ME0SegmentCollection>( edm::InputTag(iConfig.getParameter<std::string>("segmentCollection")) );
		  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag(iConfig.getParameter<std::string>("recHitCollection")) );


		  strToken_   = consumes<reco::SimToRecoCollection>( edm::InputTag("trackingParticleRecoTrackAsssociation") );
		  tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );
		  retrack_token_ = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );

		  muonsToken = consumes<std::vector<reco::ME0Muon>>(iConfig.getParameter<edm::InputTag>("muonsTag"));


		  simMuon_pt             = tree.addMulti<float>  (     "simMuon_pt"                      ,0);
		  simMuon_eta            = tree.addMulti<float>  (     "simMuon_eta"                     ,0);
		  simMuon_phi            = tree.addMulti<float>  (     "simMuon_phi"                     ,0);
		  simMuon_q              = tree.addMulti<int>    (     "simMuon_q"                       ,0);
		  simMuon_segmentQuality = tree.addMulti<int>    (     "simMuon_segmentQuality"          ,0);
		  simMuon_segmentIDX     = tree.addMulti<int>    (     "simMuon_segmentIDX"              ,0);
		 simMuon_trackIDX        = tree.addMulti<int>    (     "simMuon_trackIDX"                ,0);
		 simMuon_track_pt        = tree.addMulti<float>  (     "simMuon_track_pt"                ,0);
		 simMuon_track_eta       = tree.addMulti<float>  (     "simMuon_track_eta"               ,0);
		 simMuon_track_phi       = tree.addMulti<float>  (     "simMuon_track_phi"               ,0);

		 simMuon_gen_nLays       = tree.addMulti<int>    (     "simMuon_gen_nLays"                ,0);
		 simMuon_gen_eta         = tree.addMulti<float>  (     "simMuon_gen_eta"                ,0);
		 simMuon_gen_phi         = tree.addMulti<float>  (     "simMuon_gen_phi"                ,0);
		 simMuon_gen_dphi        = tree.addMulti<float>  (     "simMuon_gen_dphi"               ,0);
		 simMuon_gen_deta        = tree.addMulti<float>  (     "simMuon_gen_deta"               ,0);
		 simMuon_gen_x           = tree.addMulti<float>  (     "simMuon_gen_x"                  ,0);
		 simMuon_gen_y           = tree.addMulti<float>  (     "simMuon_gen_y"                  ,0);
		 simMuon_gen_dx          = tree.addMulti<float>  (     "simMuon_gen_dx"                 ,0);
		 simMuon_gen_dy          = tree.addMulti<float>  (     "simMuon_gen_dy"                 ,0);

		 me0Muon_p              = tree.addMulti<float>  (     "me0Muon_p"                      ,0);
		  me0Muon_pt             = tree.addMulti<float>  (     "me0Muon_pt"                      ,0);
		  me0Muon_eta            = tree.addMulti<float>  (     "me0Muon_eta"                     ,0);
		  me0Muon_phi            = tree.addMulti<float>  (     "me0Muon_phi"                     ,0);
		  me0Muon_q              = tree.addMulti<int>    (     "me0Muon_q"                       ,0);
		  me0Muon_trackTruthType    = tree.addMulti<int>    (     "me0Muon_trackTruthType"                       ,0);
		  me0Muon_segTrackTruthType = tree.addMulti<int>    (     "me0Muon_segTrackTruthType"                       ,0);
		 me0Muon_trackIDX         = tree.addMulti<int>    (     "me0Muon_trackIDX"                  ,0);
		  me0Muon_segIDX         = tree.addMulti<int>    (     "me0Muon_segIDX"                  ,0);
		  me0Muon_truthType      = tree.addMulti<int>    (     "me0Muon_truthType"               ,0);
		  me0Muon_track_eta      = tree.addMulti<float>  (     "me0Muon_track_eta"               ,0);
		  me0Muon_track_phi      = tree.addMulti<float>  (     "me0Muon_track_phi"               ,0);
		  me0Muon_track_dphi     = tree.addMulti<float>  (     "me0Muon_track_dphi"              ,0);
		  me0Muon_track_deta     = tree.addMulti<float>  (     "me0Muon_track_deta"              ,0);
		  me0Muon_track_x        = tree.addMulti<float>  (     "me0Muon_track_x"                 ,0);
		  me0Muon_track_y        = tree.addMulti<float>  (     "me0Muon_track_y"                 ,0);
		  me0Muon_track_dx       = tree.addMulti<float>  (     "me0Muon_track_dx"                ,0);
		  me0Muon_track_dy       = tree.addMulti<float>  (     "me0Muon_track_dy"                ,0);
		  me0Muon_track_sigx     = tree.addMulti<float>  (     "me0Muon_track_sigx"              ,0);
		  me0Muon_track_sigy     = tree.addMulti<float>  (     "me0Muon_track_sigy"              ,0);
		  me0Muon_track_sigdx    = tree.addMulti<float>  (     "me0Muon_track_sigdx"             ,0);
		  me0Muon_track_sigdy    = tree.addMulti<float>  (     "me0Muon_track_sigdy"             ,0);
		  me0Muon_segment_eta    = tree.addMulti<float>  (     "me0Muon_segment_eta"             ,0);
		  me0Muon_segment_phi    = tree.addMulti<float>  (     "me0Muon_segment_phi"             ,0);
		  me0Muon_segment_dphi   = tree.addMulti<float>  (     "me0Muon_segment_dphi"            ,0);
		  me0Muon_segment_deta   = tree.addMulti<float>  (     "me0Muon_segment_deta"            ,0);
		  me0Muon_segment_x      = tree.addMulti<float>  (     "me0Muon_segment_x"               ,0);
		  me0Muon_segment_y      = tree.addMulti<float>  (     "me0Muon_segment_y"               ,0);
		  me0Muon_segment_dx     = tree.addMulti<float>  (     "me0Muon_segment_dx"              ,0);
		  me0Muon_segment_dy     = tree.addMulti<float>  (     "me0Muon_segment_dy"              ,0);
		  me0Muon_segment_sigx   = tree.addMulti<float>  (     "me0Muon_segment_sigx"            ,0);
		  me0Muon_segment_sigy   = tree.addMulti<float>  (     "me0Muon_segment_sigy"            ,0);
		  me0Muon_segment_sigdx  = tree.addMulti<float>  (     "me0Muon_segment_sigdx"           ,0);
		  me0Muon_segment_sigdy  = tree.addMulti<float>  (     "me0Muon_segment_sigdy"           ,0);

		  tree.book();




	}

	~ME0TrackMatchingTreeMaker() {
//		hists.write(outFileName);
		  tree.write();

	}



private:
	virtual void beginJob() {};
	virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);




void fillTree(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, const ME0RecHitCollection& recHits, const std::vector<SimTrack>&  simTrackH, const std::vector<PSimHit>&  simHitH,
		const std::vector<SimVertex>& verts, const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap,
		const edm::Handle<TrackingParticleCollection>& trParticles, const reco::SimToRecoCollection& simToReco , const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp,
		const std::vector<reco::ME0Muon> & me0Muons, const reco::TrackCollection& retracks ){

	auto simMuons = ME0Helper::fillSimMuons(simTrackH,simHitH);
	ME0Helper::DigiInfoMap digiInfo;
	ME0Helper::fillDigiInfoMap(newDigis,digiMap,oldDigis,simHitH,digiInfo);
	ME0Helper::associateSimMuons(simMuons,digiInfo);
	ME0Helper::SegmentCatMap segmentCats;
	ME0Helper::fillSegmentCategories(segments,recHits,digiInfo,simMuons,segmentCats);
	ME0Helper::associateSimMuonsToTracks(simMuons, trParticles, simToReco);
//	auto trackTruths = ME0Helper::assignTrackTruth(trParticles,verts,simToReco);

	  tree.reset();

	  std::vector<int> goodSeg  (segments.size(),0);
	  std::vector<int> goodTrack(retracks.size(),0);


	  for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
			const auto& muon = simMuons[iM];
			const float absETA = TMath::Abs(muon.track->momentum().eta());
			if(absETA > 2.8 || absETA < 2.0 ) continue;

		    tree.fillMulti(simMuon_pt             ,float(muon.track->momentum().pt()));
		    tree.fillMulti(simMuon_eta            ,float(muon.track->momentum().eta()));
		    tree.fillMulti(simMuon_phi            ,float(muon.track->momentum().phi()));
		    tree.fillMulti(simMuon_q              ,  int(muon.track->charge()));

		    bool foundTrack = !muon.recoTrack.isNull();

			int typeC = 3; //Perfect/MissOne/MissMultiple/Lost
			if(muon.segments.size()) {
				int nRD = (segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second)->nRecHits();
				int nSD = muon.segments[0].second;
				if(segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first == ME0Helper::MUON_COMP_PURE) typeC = 0;
				else if(nRD - nSD <= 1 ) typeC = 1;
				else if( nSD > (nRD - nSD)) typeC = 2;
				//Time cut
				const ME0Segment * segment = &*(segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second);
				if(std::fabs(segment->time()) >= 11.0) typeC = 3;
			}
		    tree.fillMulti(simMuon_segmentQuality ,  typeC);//Perfect / MissOne/ MissSome /Lost
		    int segIDX = -1;

		    if(typeC <= 2){
		    	segIDX = int(segments.get(muon.segments[0].first.first).first - segments.begin()) +  muon.segments[0].first.second;
		    	goodSeg[   segIDX  ] =1;
		    }
		    tree.fillMulti(simMuon_segmentIDX     ,  segIDX);

		    auto genProp = ME0Helper::getSimTrackProperties(mgeom,muon.simHits);
		    tree.fillMulti(simMuon_gen_nLays     ,genProp.nLaysHit);
	    	tree.fillMulti(simMuon_gen_eta       ,float(genProp.cenEta));
	    	tree.fillMulti(simMuon_gen_phi       ,float(genProp.cenPhi));
	    	tree.fillMulti(simMuon_gen_dphi      ,float(genProp.dPhi));
	    	tree.fillMulti(simMuon_gen_deta      ,float(genProp.dEta));
	    	tree.fillMulti(simMuon_gen_x         ,float(genProp.theOrigin.x()));
	    	tree.fillMulti(simMuon_gen_y         ,float(genProp.theOrigin.y()));
	    	tree.fillMulti(simMuon_gen_dx        ,genProp.nLaysHit > 1 ? float(genProp.theLocalDirection.x() / genProp.theLocalDirection.z() )  : float(0.0));
	    	tree.fillMulti(simMuon_gen_dy        ,genProp.nLaysHit > 1 ? float(genProp.theLocalDirection.y() / genProp.theLocalDirection.z() )  : float(0.0));


//		    int me0MuonIDX = -1;
//		    if(foundTrack)
//		    for(unsigned int iM = 0; iM < me0Muons.size(); ++iM){
//		    	if(me0Muons[iM].innerTrack().key() != muon.recoTrack.key() ) continue;
//		    	me0MuonIDX = iM;
//		    	goodTrack[iM] = 1;
//		    }
		    if(foundTrack) goodTrack[muon.recoTrack.key()] = 1;
		    tree.fillMulti(simMuon_trackIDX         , foundTrack ? int(muon.recoTrack.key()) : int(-1));
		    tree.fillMulti(simMuon_track_pt         , foundTrack ? float(muon.recoTrack->pt() ) : float(0));
		    tree.fillMulti(simMuon_track_eta        , foundTrack ? float(muon.recoTrack->eta()) : float(0));
		    tree.fillMulti(simMuon_track_phi        , foundTrack ? float(muon.recoTrack->phi()) : float(0));
	  }

	    for(unsigned int iM = 0; iM < me0Muons.size(); ++iM){
	    	const auto& me0Muon = me0Muons[iM];
			tree.fillMulti(me0Muon_p             ,float(me0Muon.innerTrack()->p()));
	    	tree.fillMulti(me0Muon_pt            ,float(me0Muon.innerTrack()->pt()));
	    	tree.fillMulti(me0Muon_eta           ,float(me0Muon.innerTrack()->eta()));
	    	tree.fillMulti(me0Muon_phi           ,float(me0Muon.innerTrack()->phi()));
	    	tree.fillMulti(me0Muon_q             ,int(me0Muon.innerTrack()->charge()));

//	    	auto trkTruthInfo = trackTruths.find(me0Muon.innerTrack().key());
//	    	auto trkTruth = trkTruthInfo == trackTruths.end() ? ME0Helper::OTHER : trkTruthInfo->second.second;
//	    	tree.fillMulti(me0Muon_trackTruthType     ,int(trkTruth	));
	    	tree.fillMulti(me0Muon_trackTruthType     ,int(0	));

	    	ME0Helper::SegmentTrackTruthType segTrkTruth = ME0Helper::NEUTRON_SEG;
//	    	if(trkTruthInfo != trackTruths.end() && trkTruthInfo->second.first >= 0 )
//	    		segTrkTruth = getSegmentMatch(digiInfo,recHits,me0Muon.me0segment(),TrackingParticleRef(trParticles, trkTruthInfo->second.first));
	    	tree.fillMulti(me0Muon_segTrackTruthType  ,int(segTrkTruth));


			 tree.fillMulti(me0Muon_trackIDX, int(me0Muon.innerTrack().key()));
	    	tree.fillMulti(me0Muon_segIDX        ,int(me0Muon.me0segid()));

	    	int truthType  = 3;
	    	if(goodTrack[me0Muon.innerTrack().key()] && goodSeg[   me0Muon.me0segid()  ]  ) truthType = 0;
	    	else if(goodTrack[iM] && !goodSeg[   me0Muon.me0segid()  ]  ) truthType = 1;
	    	else if(!goodTrack[iM] && goodSeg[   me0Muon.me0segid()  ]  ) truthType = 2;
	    	tree.fillMulti(me0Muon_truthType     ,int(truthType)); //RealMu  / Fake Seg / Fake Track / Fake Fake

	    	tree.fillMulti(me0Muon_track_eta     ,float(me0Muon.globalTrackPosAtSurface().eta()));
	    	tree.fillMulti(me0Muon_track_phi     ,float(me0Muon.globalTrackPosAtSurface().phi()));

	    	LocalTrajectoryParameters ltp(me0Muon.localTrackPosAtSurface() ,me0Muon.localTrackMomAtSurface(),me0Muon.trackCharge());
	    	auto locTrackP = ltp.vector();

	    	const auto * chamber = mgeom->chamber(me0Muon.me0segment().me0DetId());
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
	    	float trackDPhi, trackDEta;
	    	getDPhiDEta(LocalPoint(locTrackP[3],locTrackP[4],0),LocalVector(locTrackP[1],locTrackP[2],1),trackDPhi,trackDEta);
	    	tree.fillMulti(me0Muon_track_dphi    ,float(trackDPhi));
	    	tree.fillMulti(me0Muon_track_deta    ,float(trackDEta));
	    	tree.fillMulti(me0Muon_track_x       ,float(locTrackP[3]));
	    	tree.fillMulti(me0Muon_track_y       ,float(locTrackP[4]));
	    	tree.fillMulti(me0Muon_track_dx      ,float(locTrackP[1]));
	    	tree.fillMulti(me0Muon_track_dy      ,float(locTrackP[2]));
	    	tree.fillMulti(me0Muon_track_sigx    ,float(std::sqrt(me0Muon.localTrackCov()[3][3])));
	    	tree.fillMulti(me0Muon_track_sigy    ,float(std::sqrt(me0Muon.localTrackCov()[4][4])));
	    	tree.fillMulti(me0Muon_track_sigdx   ,float(std::sqrt(me0Muon.localTrackCov()[1][1])));
	    	tree.fillMulti(me0Muon_track_sigdy   ,float(std::sqrt(me0Muon.localTrackCov()[2][2])));

	    	GlobalPoint segGlobPos(chamber->toGlobal(me0Muon.me0segment().localPosition()));
	    	tree.fillMulti(me0Muon_segment_eta   ,float(segGlobPos.eta()));
	    	tree.fillMulti(me0Muon_segment_phi   ,float(segGlobPos.phi()));

	    	float segDPhi, segDEta;
	    	getDPhiDEta(me0Muon.me0segment().localPosition(),me0Muon.me0segment().localDirection(),segDPhi,segDEta);
	    	tree.fillMulti(me0Muon_segment_dphi  ,float(segDPhi));
	    	tree.fillMulti(me0Muon_segment_deta  ,float(segDEta));
	    	tree.fillMulti(me0Muon_segment_x     ,float(me0Muon.me0segment().localPosition().x()));
	    	tree.fillMulti(me0Muon_segment_y     ,float(me0Muon.me0segment().localPosition().y()));
	    	tree.fillMulti(me0Muon_segment_dx    ,float(me0Muon.me0segment().localDirection().x()/me0Muon.me0segment().localDirection().z()));
	    	tree.fillMulti(me0Muon_segment_dy    ,float(me0Muon.me0segment().localDirection().y()/me0Muon.me0segment().localDirection().z()));
	    	tree.fillMulti(me0Muon_segment_sigx  ,float(std::sqrt(me0Muon.me0segment().parametersError()[2][2])));
	    	tree.fillMulti(me0Muon_segment_sigy  ,float(std::sqrt(me0Muon.me0segment().parametersError()[3][3])));
	    	tree.fillMulti(me0Muon_segment_sigdx ,float(std::sqrt(me0Muon.me0segment().parametersError()[0][0])));
	    	tree.fillMulti(me0Muon_segment_sigdy ,float(std::sqrt(me0Muon.me0segment().parametersError()[1][1])));
	    }



	  tree.fillTree();


}

private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoMap>                 oldToNewMapToken_  ;
    edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
    edm::EDGetTokenT<std::vector<SimVertex>> vert_token;

    edm::EDGetTokenT<ME0SegmentCollection>          segToken_ ;
    edm::EDGetTokenT<std::vector<SimTrack>> track_token;
    edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;

    edm::EDGetTokenT<reco::TrackCollection>           retrack_token_ ;
    edm::EDGetTokenT<reco::SimToRecoCollection>          strToken_ ;
    edm::EDGetTokenT<TrackingParticleCollection>          tpToken_ ;

    edm::EDGetTokenT<std::vector<reco::ME0Muon>>          muonsToken ;


	TString outFileName;
	TString runName;

    TreeWriter tree;
    ASTypes::size simMuon_pt                    =-1;
    ASTypes::size simMuon_eta                   =-1;
    ASTypes::size simMuon_phi                   =-1;
    ASTypes::size simMuon_q                     =-1;
    ASTypes::size simMuon_segmentQuality        =-1; //Perfect / MissOne/ MissSome /Lost
    ASTypes::size simMuon_segmentIDX            =-1;
    ASTypes::size simMuon_trackIDX              =-1;
    ASTypes::size simMuon_track_pt              =-1;
    ASTypes::size simMuon_track_eta             =-1;
    ASTypes::size simMuon_track_phi             =-1;
    ASTypes::size simMuon_gen_nLays             =-1;
    ASTypes::size simMuon_gen_eta               =-1;
    ASTypes::size simMuon_gen_phi               =-1;
    ASTypes::size simMuon_gen_dphi              =-1;
    ASTypes::size simMuon_gen_deta              =-1;
    ASTypes::size simMuon_gen_x                 =-1;
    ASTypes::size simMuon_gen_y                 =-1;
    ASTypes::size simMuon_gen_dx                =-1;
    ASTypes::size simMuon_gen_dy                =-1;
    ASTypes::size me0Muon_p                     =-1;
    ASTypes::size me0Muon_pt                    =-1;
    ASTypes::size me0Muon_eta                   =-1;
    ASTypes::size me0Muon_phi                   =-1;
    ASTypes::size me0Muon_q                     =-1;
    ASTypes::size me0Muon_trackTruthType        =-1;
    ASTypes::size me0Muon_segTrackTruthType     =-1;
    ASTypes::size me0Muon_trackIDX              =-1;
    ASTypes::size me0Muon_segIDX                =-1;
    ASTypes::size me0Muon_truthType             =-1; //RealMu  / Fake Seg / Fake Track / Fake Fake
    ASTypes::size me0Muon_track_eta             =-1;
    ASTypes::size me0Muon_track_phi             =-1;
    ASTypes::size me0Muon_track_dphi            =-1;
    ASTypes::size me0Muon_track_deta            =-1;
    ASTypes::size me0Muon_track_x               =-1;
    ASTypes::size me0Muon_track_y               =-1;
    ASTypes::size me0Muon_track_dx              =-1;
    ASTypes::size me0Muon_track_dy              =-1;
    ASTypes::size me0Muon_track_sigx            =-1;
    ASTypes::size me0Muon_track_sigy            =-1;
    ASTypes::size me0Muon_track_sigdx           =-1;
    ASTypes::size me0Muon_track_sigdy           =-1;
    ASTypes::size me0Muon_segment_eta           =-1;
    ASTypes::size me0Muon_segment_phi           =-1;
    ASTypes::size me0Muon_segment_dphi          =-1;
    ASTypes::size me0Muon_segment_deta          =-1;
    ASTypes::size me0Muon_segment_x             =-1;
    ASTypes::size me0Muon_segment_y             =-1;
    ASTypes::size me0Muon_segment_dx            =-1;
    ASTypes::size me0Muon_segment_dy            =-1;
    ASTypes::size me0Muon_segment_sigx          =-1;
    ASTypes::size me0Muon_segment_sigy          =-1;
    ASTypes::size me0Muon_segment_sigdx         =-1;
    ASTypes::size me0Muon_segment_sigdy         =-1;

};




void
ME0TrackMatchingTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

	  edm::Handle <std::vector<SimVertex> > verts;
	  iEvent.getByToken(vert_token,verts);

	  edm::Handle<ME0SegmentCollection >  segH ;
	  iEvent.getByToken(segToken_,segH);

	  edm::Handle<ME0RecHitCollection> rechitsH;
	  iEvent.getByToken(rhToken_,rechitsH);

	edm::ESHandle<ME0Geometry> me0g;
	iSetup.get<MuonGeometryRecord>().get(me0g);
	const ME0Geometry* mgeom = &*me0g;


	  edm::ESHandle<MagneticField> bField;
	  iSetup.get<IdealMagneticFieldRecord>().get(bField);

	  edm::ESHandle<Propagator> ThisshProp;
	  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", ThisshProp);

	  edm::Handle <reco::TrackCollection > retracks;
	  iEvent.getByToken(retrack_token_,retracks);

	  edm::Handle<reco::SimToRecoCollection> simRecH;
	  iEvent.getByToken(strToken_,simRecH);

	  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
	  iEvent.getByToken(tpToken_,TPCollectionH);


	  edm::Handle<std::vector<reco::ME0Muon>> muonCollectionH;
	  iEvent.getByToken(muonsToken, muonCollectionH);

	  fillTree(mgeom,*segH,*rechitsH,*tracks,*simHitH,*verts,*digisH,*digisNH,*digisM,TPCollectionH, *simRecH,bField,ThisshProp, *muonCollectionH, *retracks);


}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0TrackMatchingTreeMaker);
