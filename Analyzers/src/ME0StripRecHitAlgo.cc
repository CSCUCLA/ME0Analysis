///** \file
// *
// *  $Date: 2013/04/24 17:16:35 $
// *  $Revision: 1.1 $
// *  \author M. Maggi -- INFN Bari
//*/
//
//#include "../interface/ME0StripRecHitAlgo.h"
//
//
//
//
//#include "DataFormats/MuonDetId/interface/ME0DetId.h"
//#include "FWCore/Utilities/interface/Exception.h"
//
//#include "Geometry/Records/interface/MuonGeometryRecord.h"
//
//
//ME0StripRecHitAlgo::ME0StripRecHitAlgo(const edm::ParameterSet& config) : nP(config.getParameter<int>("nEtaParts")),nStrips(config.getParameter<int>("nStrips")) {
//}
//
//ME0StripRecHitAlgo::~ME0StripRecHitAlgo()
//{
//}
//
//void ME0StripRecHitAlgo::setES(const edm::EventSetup& setup)
//{
//  // Get the ME0 Geometry
//  setup.get<MuonGeometryRecord>().get(me0Geom);
//}
//
//// First Step
//bool ME0StripRecHitAlgo::compute(const ME0DigiPreReco& digi,
//            LocalPoint& Point,
//            LocalError& error)  const
//{
//  LocalPoint loctemp2(digi.x(),digi.y(),0.);
//  Point = loctemp2;
//  LocalError loerr2(digi.ex()*digi.ex(),digi.corr()*digi.ex()*digi.ey(),digi.ey()*digi.ey());
//  error = loerr2;
//  return true;
//}
//
//
//#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
//#include "TMath.h"
//#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
//
//TrapezoidalStripTopology * ME0StripRecHitAlgo::buildTopo(std::vector<float>& _p) const {
//    float b = _p[0];
//    float B = _p[1];
//    float h = _p[2];
//    float r0 = h*(B + b)/(B - b);
//    float striplength = h*2;
//    float strips = _p[3];
//    float pitch = (b + B)/strips;
//    int nstrip =static_cast<int>(strips);
//    return new TrapezoidalStripTopology(nstrip, pitch, striplength, r0);
//}
//
//
//
//
//// Build all hits in the range associated to the layerId, at the 1st step.
//edm::OwnVector<ME0RecHit> ME0StripRecHitAlgo::reconstruct(const ME0DetId& me0Id,
//  const ME0DigiPreRecoCollection::Range& digiRange){
//  edm::OwnVector<ME0RecHit> result;
//
//
//    const auto * mainEtaPartition =  me0Geom->etaPartition(me0Id);
//
//    LocalPoint localTop(0,mainEtaPartition->specificTopology().stripLength()/2);
//    LocalPoint localBottom(0,-1*mainEtaPartition->specificTopology().stripLength()/2);
//
//    GlobalPoint globalTop = mainEtaPartition->toGlobal(localTop);
//    GlobalPoint globalBottom = mainEtaPartition->toGlobal(localBottom);
//
//    double bottomDistanceFromBeam = TMath::Sqrt(globalBottom.x()*globalBottom.x() + globalBottom.y()*globalBottom.y() );
//    double topDistanceFromBeam = TMath::Sqrt(globalTop.x()*globalTop.x() + globalTop.y()*globalTop.y() );
//    double middleDistanceFromBeam = (topDistanceFromBeam + bottomDistanceFromBeam)/2;
//
//    std::vector<float> origPars = mainEtaPartition->specs()->parameters();
//
//
//    double etaTop = globalTop.eta();
//    double etaBottom = globalBottom.eta();
//    double zBottom = globalBottom.z();
//
//    std::vector<double> rollTops; rollTops.reserve(nP);
//    for(unsigned int iE = 0; iE < nP; ++iE){
//      double eta = (etaTop -etaBottom)*double(iE + 1)/double(nP) + etaBottom;
//      double distFromBeam = TMath::Abs(2 * zBottom /(TMath::Exp(-1*eta) - TMath::Exp(eta)  ));
//      rollTops.push_back(distFromBeam - middleDistanceFromBeam);
//    }
//
//
//    std::vector<TrapezoidalStripTopology * > stripTopos;stripTopos.reserve(nP);
//    for(unsigned int iE = 0; iE < nP; ++iE){
//      std::vector<float> params(4,nStrips);
//
//      auto getWidth = [&] ( float locY ) -> float { return (origPars[2]*(origPars[1]+origPars[0]) +locY*(origPars[1] - origPars[0]) )/(2*origPars[2]);};
//
//      params[0] = iE == 0 ?  origPars[0] : getWidth(rollTops[iE -1]);
//      params[1] = iE +1 == nP ?  origPars[1] : getWidth(rollTops[iE]);
//      params[2] = ((iE + 1 == nP ? localTop.y() : rollTops[iE] ) - (iE  == 0 ? localBottom.y() : rollTops[iE-1] ) )/2;
//
//      stripTopos.emplace_back(buildTopo(params));
//    }
//
//    auto getEtaPartition = [&] (double locY) -> TrapezoidalStripTopology * {
//      unsigned int etaPart = nP -1;
//      for(unsigned int iE = 0; iE < nP; ++iE ){
//        if(locY <  rollTops[iE]) {etaPart = iE; break;}
//      }
//      return stripTopos[etaPart];
//    };
//
//
//  for (ME0DigiPreRecoCollection::const_iterator digi = digiRange.first;
//       digi != digiRange.second;digi++) {
//
//    TrapezoidalStripTopology * etaPartition = getEtaPartition(digi->y());
//
//    LocalPoint origPoint(digi->x(), digi->y() - (etaPartition->radius() - middleDistanceFromBeam  )  ,0.);
//    int strip = etaPartition->channel(origPoint);
//    float stripF = float(strip) -0.5;
//
//    LocalPoint  newPoint = etaPartition->localPosition(stripF);
//    LocalError  newError = etaPartition->localError(stripF, 1./sqrt(12.));
//    LocalPoint  translatedToChamberPoint(newPoint.x(),newPoint.y() + (etaPartition->radius() - middleDistanceFromBeam  ));
//
//    ME0RecHit* recHit = new ME0RecHit(me0Id,digi->tof(),translatedToChamberPoint,newError);
//    result.push_back(recHit);
//  }
//
//  for(unsigned int iE = 0; iE < stripTopos.size(); ++iE) delete stripTopos[iE];
//  return result;
//}
//
