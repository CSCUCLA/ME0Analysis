#ifndef ME0Helper_h
#define ME0Helper_h

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "DataFormats/GEMRecHit/interface/ME0Segment.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include <TVector2.h>
#include <TMath.h>
#include "../interface/MuonSegFit.h"


namespace ME0Helper{
const double beginOfDet  = 527;
const double centerOfDet = 539.5;
const double endOfDet    = 552;

std::vector<const PSimHit*> getMatchedSimHits(const SimTrack * simTrack, const std::vector<PSimHit>&  simHitH);
struct SimMuon {
	SimMuon(const SimTrack * track, std::vector<const PSimHit*> hits) : track(track), simHits(hits){}
	const SimTrack * track = 0;
	std::vector<const PSimHit*> simHits;
	std::vector<std::pair<ME0DetId,unsigned int >> digiInfos;
	std::vector< std::pair<std::pair<ME0DetId,int>, int >  > segments;
};
typedef std::vector<SimMuon> SimMuons;
SimMuons fillSimMuons(const std::vector<SimTrack>& simTracks, const std::vector<PSimHit>&  simHits);

struct DigiInfo {
	const ME0DigiPreReco* digi = 0;
	std::vector<const ME0DigiPreReco*> origDigis;
	std::vector<const PSimHit*> simHits;
	//Filled later
	std::vector<const SimTrack*> simTracks;
	std::vector<int> muonIDXs;
};

typedef std::map<ME0DetId,std::vector<DigiInfo> > DigiInfoMap;
void fillDigiInfoMap(const  ME0DigiPreRecoCollection& digis, const ME0DigiPreRecoMap& digiMap,
		const ME0DigiPreRecoCollection& oldDigis, const std::vector<PSimHit>&  simHitH,  DigiInfoMap& digiInfoMap);

void associateSimMuons(SimMuons& simMuons, DigiInfoMap& digiInfoMap );

int getRecHitIndex(const ME0RecHitCollection& recHits, const ME0RecHit& rh);

enum SegmentGenCat { MUON_COMP_PURE, MUON_COMP_DIRTY_TRACK,MUON_COMP_DIRTY_NEUT, MUON_MISS_PURE, MUON_MISS_DIRTY_TRACK,MUON_MISS_DIRTY_NEUT,
					 FAKE_TRACK_PURE, FAKE_NEUT_PURE, FAKE_MUON_MIX, FAKE_OTHER, SEGMENT_ERROR };
const std::string SegmentGenCatNames[] = {"MUON_COMP_PURE", "MUON_COMP_DIRTY_TRACK","MUON_COMP_DIRTY_NEUT", "MUON_MISS_PURE", "MUON_MISS_DIRTY_TRACK","MUON_MISS_DIRTY_NEUT",
		 "FAKE_TRACK_PURE", "FAKE_NEUT_PURE", "FAKE_MUON_MIX", "FAKE_OTHER", "SEGMENT_ERROR" };

typedef std::map<ME0DetId, std::vector<std::pair<SegmentGenCat,int> > > SegmentCatMap;
void fillSegmentCategories(const ME0SegmentCollection& segments,
		const ME0RecHitCollection& recHits,const DigiInfoMap& digiInfo,
		SimMuons& muons, SegmentCatMap& segCatMap);



typedef std::vector<std::pair<ME0DetId,const ME0DigiPreReco*> > ME0DigiList;
ME0DigiList getMatchedDigis(const std::vector<const PSimHit*>& simHits,
		const ME0DigiPreRecoCollection& oldDigis,const  ME0DigiPreRecoCollection * newDigis = 0,const ME0DigiPreRecoMap*  digiMap =0 ,
		ME0DigiList * origDigiList = 0, std::vector<int> * origDigiNewDigiIDX = 0
);


struct SimHitProperties{
	double cenEta = -1;
	double cenPhi = -1;
	double lay1Eta = -1;
	double lay1Phi = -1;
	double dPhi = -1;
	double dEta = -1;
	int nLaysHit = 0;
	bool oneChamber = true;
	bool inFirst= false;
};
SimHitProperties getSimTrackProperties( const ME0Geometry* mgeom, const std::vector<const PSimHit* >& simHits);

ME0Segment * buildSegment(const ME0Geometry* mgeom, const std::vector<std::pair<ME0DetId,const ME0DigiPreReco*> >& hits);
struct SegmentExtrapPoint {
	LocalPoint point;
	LocalError pointError;
	double dxdz = -1;
	double dydz = -1;
	double cov_dxdz_dxdz = -1;
	double cov_dydz_dydz = -1;
};
struct SegmentProperties{
	//Local center info
	SegmentExtrapPoint segAtCenter;
	LocalPoint initialPoint;
	double cenEta   = -1;
	double cenPhi   = -1;
	double beginEta = -1;
	double beginPhi = -1;
	double endEta   = -1;
	double endPhi   = -1;
	double dPhi     = -1;
	double dEta     = -1;
};
SegmentProperties getSegmentProperties(const GeomDet* loc, const ME0Segment * segment);


}


#endif
