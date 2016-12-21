#ifndef ME0Helper_h
#define ME0Helper_h

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include <TVector2.h>
#include <TMath.h>
#include "../interface/MuonSegFit.h"


namespace ME0Helper{
const double beginOfDet  = 527;
const double centerOfDet = 539.5;
const double endOfDet    = 552;


std::vector<const PSimHit*> getMatchedSimHits(const SimTrack * simTrack, edm::Handle<std::vector<PSimHit> >&  simHitH){
	std::vector<const PSimHit*> hitList;
    for (const auto& simHit : (*simHitH)) {
        if(simHit.eventId() != simTrack->eventId()) continue;
        if(simHit.trackId() != simTrack->trackId()) continue;
        hitList.push_back(&simHit);
    }
    return hitList;
}

typedef std::vector<std::pair<ME0DetId,const ME0DigiPreReco*> > ME0DigiList;

ME0DigiList getMatchedDigis(const std::vector<const PSimHit*>& simHits,
		const ME0DigiPreRecoCollection& oldDigis,const  ME0DigiPreRecoCollection * newDigis = 0,const ME0DigiPreRecoMap*  digiMap =0 , ME0DigiList * origDigiList = 0
) {
	if(origDigiList) origDigiList->clear();
	ME0DigiList hitList;
	for(const auto* simHit : simHits){
        const ME0DetId detID = simHit->detUnitId();
        const ME0DigiPreRecoCollection::Range& range = oldDigis.get(detID);
        unsigned int iD = 0;
        double minDR2 = -1;
        const ME0DigiPreReco * digi = 0;
        for (ME0DigiPreRecoCollection::const_iterator idigi = range.first;
            idigi != range.second;idigi++) {
          double dr2 = (simHit->entryPoint().x() - idigi->x())*(simHit->entryPoint().x() - idigi->x()) + (simHit->entryPoint().y() - idigi->y())*(simHit->entryPoint().y() - idigi->y());
          if(minDR2 < 0 || dr2 < minDR2){
            digi = &(*idigi);
            minDR2 = dr2;
            iD = int(idigi-range.first);
          }
        }
        if(minDR2 < 0 ) {std::cout << "Could not find matching digi!!!" <<std::endl; continue;}
        if(newDigis && digiMap){
            const ME0DigiPreRecoCollection::Range& newRange = newDigis->get(detID);
            const ME0DigiPreRecoMap::Range& oldToNewMap = digiMap->get(detID);
        	int newIDX =  oldToNewMap.first[iD];
        	if(newIDX >= 0) hitList.emplace_back(detID,&newRange.first[newIDX]);
        	if(origDigiList) origDigiList->emplace_back(detID,digi);
        }
        else  {
        	hitList.emplace_back(detID,digi);
        }
	}
    return hitList;
}


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
SimHitProperties getSimTrackProperties( const ME0Geometry* mgeom, const std::vector<const PSimHit* >& simHits) {

	SimHitProperties prop;

    std::vector<std::vector<const PSimHit*> > laysHit(6);
    for(const auto& h : simHits){
      laysHit[ME0DetId(h->detUnitId()).layer() - 1].push_back(h);
    }
    std::vector<const PSimHit *> prunedHits;
    for(auto& l : laysHit){
      if(l.size() >0) {
    	  prop.nLaysHit++;
        prunedHits.push_back(l.front());
      }
    }

    for(unsigned int iH = 0; iH < prunedHits.size(); ++iH){
    	if(!iH) continue;
    	if(ME0DetId(prunedHits[iH]->detUnitId()).chamber() !=  ME0DetId(prunedHits[0]->detUnitId()).chamber())
    		prop.oneChamber = false;
    }


    auto getGlobal = [&](const PSimHit * sh) -> GlobalPoint {
    	return mgeom->etaPartition(sh->detUnitId())->toGlobal(sh->entryPoint());
    };

    if(laysHit[0].size()){
    	prop.inFirst = true;
    	prop.lay1Eta = getGlobal(laysHit[0].front()).eta();
    	prop.lay1Phi = getGlobal(laysHit[0].front()).phi();
    }

    if(prop.nLaysHit > 1){
        auto getProj = [](double zValue, const GlobalPoint& up, const GlobalPoint& down) -> GlobalPoint {

          auto getCenter = [](double v1,double v2,double z1,double z2,double zc)->double{
            if(z2 == z1) return -99;
            double m = (v2 - v1)/(z2 - z1);
            double b =  (z2*v1 - z1*v2)/(z2 - z1);
            return m*zc+b;
          };

          const double zc = up.z() > 0 ? zValue : -1*zValue;
          double xc = getCenter(up.x(),down.x(),up.z(),down.z(),zc);
          double yc = getCenter(up.y(),down.y(),up.z(),down.z(),zc);
          return GlobalPoint(xc,yc,zc);
        };

        const int upInd = floor(prunedHits.size()/2);
        GlobalPoint centerUp = getGlobal(prunedHits[upInd]);
        GlobalPoint centerDown = getGlobal(prunedHits[upInd-1]);
        GlobalPoint centerPt = getProj(centerOfDet, centerUp, centerDown);
        GlobalPoint upPt = getProj(endOfDet, getGlobal(prunedHits[prunedHits.size() - 1]), getGlobal(prunedHits[prunedHits.size() - 2]));
        GlobalPoint downPt = getProj(beginOfDet, getGlobal(prunedHits[1]), getGlobal(prunedHits[0]));

        prop.cenEta = centerPt.eta();
        prop.cenPhi = centerPt.phi();
        prop.dPhi = TVector2::Phi_mpi_pi(upPt.phi() - downPt.phi());
        prop.dEta = upPt.eta() - downPt.eta();

    }
    return prop;
}

ME0Segment * buildSegment(const ME0Geometry* mgeom, const std::vector<std::pair<ME0DetId,const ME0DigiPreReco*> >& hits){
	if(hits.size() < 2) return 0;

	MuonSegFit::MuonRecHitContainer muonRecHits;
	std::vector<const ME0RecHit*> muonRH;

	const ME0DetId refID=hits[0].first;
	const ME0EtaPartition * refPart = mgeom->etaPartition(refID);
	for(const auto& hp : hits ){
	    const ME0EtaPartition * thePartition   =   mgeom->etaPartition(hp.first);

	    GlobalPoint gp = thePartition->toGlobal(LocalPoint(hp.second->x(),hp.second->y(),0));
	    const LocalPoint lp = refPart->toLocal(gp);

	    double errorXX = hp.second->ex()*hp.second->ex();
	    double errorYY = hp.second->ex()*hp.second->ex();
	    double errorXY = hp.second->corr()*hp.second->ex()*hp.second->ey();
	    LocalError le(std::max(errorXX,1.6E-5), std::max(errorXY,1.6E-6), std::max(errorYY,1.6E-5));

	    ME0RecHit* recHit = new ME0RecHit(refID,hp.second->tof(),lp,le);
	    muonRH.push_back(recHit);
	    MuonSegFit::MuonRecHitPtr trkRecHit(recHit);
	    muonRecHits.push_back(trkRecHit);
	}
	  MuonSegFit  * sfit_ = new MuonSegFit(muonRecHits);
	  sfit_->fit();

	  // obtain all information necessary to make the segment:
	  LocalPoint protoIntercept      = sfit_->intercept();
	  LocalVector protoDirection     = sfit_->localdir();
	  AlgebraicSymMatrix protoErrors = sfit_->covarianceMatrix();
	  double protoChi2               = sfit_->chi2();

	  ME0Segment * segment =new ME0Segment(muonRH, protoIntercept, protoDirection, protoErrors, protoChi2, 0.0, 5.0);

	  for (auto rh:muonRecHits) rh.reset();
	  delete sfit_;

	  return segment;
}
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
SegmentProperties getSegmentProperties(const ME0Geometry* mgeom, const ME0Segment * segment){

	const ME0EtaPartition * refPart = mgeom->etaPartition(segment->me0DetId());
	auto ld = segment->localDirection();
	auto lp = segment->localPosition();
	AlgebraicSymMatrix covMatrix = segment->parametersError();

	auto extrap = [&] (double extZ) -> SegmentExtrapPoint {
	    //covmatrix wors like 0: dx/dz dy/dz, x, y
	    double extX = lp.x()+extZ*ld.x()/ld.z();
	    double extY = lp.y()+extZ*ld.y()/ld.z();

	    double exx   =  extZ*extZ *covMatrix [0][0] + 2*extZ*covMatrix[0][2] + covMatrix[2][2];
	    double eyy   =  extZ*extZ *covMatrix [1][1] + 2*extZ*covMatrix[1][3] + covMatrix[3][3];
	    double exy   =  extZ*extZ *covMatrix [0][1] +   extZ*(covMatrix[1][2]+covMatrix[0][3]) + covMatrix[2][3];
	    double edxdx = covMatrix [0][0];
	    double edydy = covMatrix [1][1];
	    SegmentExtrapPoint extP;
	    extP.point =LocalPoint(extX,extY,extZ);
	    extP.pointError =LocalError(exx,exy,eyy);
	    extP.dxdz = ld.x()/ld.z();
	    extP.dydz = ld.y()/ld.z();
	    extP.cov_dxdz_dxdz = edxdx;
	    extP.cov_dydz_dydz = edydy;
	    return extP;
	  };

	  auto layCen =  refPart->toGlobal(LocalPoint(0,0,0));
	  auto locLow  = extrap( beginOfDet  - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());
	  auto lochigh = extrap( endOfDet    - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());
	  auto loccen  = extrap( centerOfDet - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());

	  auto globlow  = refPart->toGlobal(locLow .point);
	  auto globhigh = refPart->toGlobal(lochigh.point);
	  auto globcen  = refPart->toGlobal(loccen .point);

	  SegmentProperties properties;

	  properties.segAtCenter = loccen;
	  properties.initialPoint = lp;
	  properties.cenEta   =  globcen.eta();
	  properties.cenPhi   =  globcen.phi();
	  properties.beginEta =  globlow.eta();
	  properties.beginPhi =  globlow.phi();
	  properties.endEta   =  globhigh.eta();
	  properties.endPhi   =  globhigh.phi();
	  properties.dPhi     =  TVector2::Phi_mpi_pi(globhigh.phi() - globlow.phi());
	  properties.dEta     =  globhigh.eta() - globlow.eta();
	  return properties;
}


}


#endif
