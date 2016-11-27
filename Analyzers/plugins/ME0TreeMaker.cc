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

#include "../../ME0Analysis/plugins/HistGetter.h"
#include "../interface/TreeWriter.h"
#include "../interface/MuonSegFit.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "DataFormats/GEMRecHit/interface/ME0Segment.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"


#include <TVector2.h>
#include <TMath.h>

double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return deta*deta + dphi*dphi;
}




template<typename Field, typename Object>
struct lesserFirst : public std::binary_function<const std::pair<Field,Object>&, const std::pair<Field,Object>&, Bool_t> {
  Bool_t operator()(const std::pair<Field,Object>& x, const std::pair<Field,Object>& y) const {
    return x.first < y.first;
  }
};

struct HitCollection {

  void add(const ME0DetId& id, const PSimHit * sh, const ME0DigiPreReco * d, const ME0RecHit * rh){
    detIDs.push_back(id);
    simHits.push_back(sh);
    digis.push_back(d);
    recHits.push_back(rh);
  }

  void sort() {
    unsigned int nL = detIDs.size();

    //Find majority chamber
    std::vector<int> chambers; chambers.reserve(nL);
    for(unsigned int i = 0; i < nL ; ++i) chambers.push_back(detIDs[i].chamber());
    std::sort(chambers.begin(), chambers.end());
    int nCM = 0;
    int chamber = -1;
    int tCM =0;
    int tchamber =chambers.front();
    for(unsigned int i = 0; i < chambers.size(); ++i){
      if(chambers[i] == tchamber ) {
        tCM++;
        continue;
      }
      if(tCM > nCM) {
        chamber = chambers[i-1];
        nCM = tCM;
      }
      tchamber = chambers[i];
      tCM = 1;
    }
    if(tCM > nCM) {
      chamber = chambers.back();
    }

    //sort by layer
    std::vector<std::pair<int,int> > sortables; sortables.reserve(nL);
    for(unsigned int i = 0; i < nL ; ++i) {
      if(detIDs[i].chamber() != chamber) continue;
      sortables.emplace_back(detIDs[i].layer(),i);
    }
    std::sort(sortables.begin(), sortables.end(), lesserFirst<int,int>());
    std::vector<      ME0DetId>         _detIDs   ;_detIDs   .reserve(nL);
    std::vector<const PSimHit*>        _simHits   ;_simHits  .reserve(nL);
    std::vector<const ME0DigiPreReco*> _digis     ;_digis    .reserve(nL);
    std::vector<const ME0RecHit*>      _recHits   ;_recHits  .reserve(nL);


    for(unsigned int i = 0; i < sortables.size(); ++i){
      if(i && detIDs[sortables[i-1].second].layer() == detIDs   [sortables[i].second].layer()) continue;
      _detIDs   .push_back(detIDs   [sortables[i].second]);
      _simHits  .push_back(simHits  [sortables[i].second]);
      _digis    .push_back(digis    [sortables[i].second]);
      _recHits  .push_back(recHits  [sortables[i].second]);
    }
    detIDs    =_detIDs   ;
    simHits   =_simHits  ;
    digis     =_digis    ;
    recHits   =_recHits  ;
  }
  void clear(){
    detIDs    .clear();
    simHits   .clear();
    digis     .clear();
    recHits   .clear();
  }



  std::vector<ME0DetId>              detIDs;
  std::vector<const PSimHit*>        simHits;
  std::vector<const ME0DigiPreReco*> digis;
  std::vector<const ME0RecHit*>      recHits;

};



class ME0TreeMaker : public edm::EDAnalyzer {
    public:
        explicit ME0TreeMaker(const edm::ParameterSet&);
        ~ME0TreeMaker();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);

        void getHits(edm::Handle<std::vector<PSimHit> >&  simHitH,edm::Handle<ME0DigiPreRecoCollection>& digis,edm::Handle<ME0RecHitCollection>& rechits,
            const SimTrack * simTrack, HitCollection& hits
        );
        void getSimInfo(const ME0Geometry* mgeom, const HitCollection& hits, float& gen_eta, float& gen_phi, float& gen_delta_eta, float& gen_delta_phi);
        std::pair<ME0DetId,ME0Segment*> buildSegment(const ME0Geometry* mgeom, const HitCollection& hits) const;
        void getSegmentProperties(const ME0Segment * seg, const ME0EtaPartition* partition,
            float& cen_eta, float& cen_phi, float& delta_eta, float& delta_phi) const;

        FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
                 int , const AlgebraicSymMatrix66& ,
                 const MagneticField* );

        FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& ,
                 int , const AlgebraicSymMatrix55& ,
                 const MagneticField* );

        void getFromFTS(const FreeTrajectoryState& ,
            GlobalVector& , GlobalVector& ,
            int& , AlgebraicSymMatrix66& );

        GlobalPoint getPropogatedTrack(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, double extZ);
        void fillTrackProperties(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, bool isMatched);

        virtual void endJob() {};


    private:
        edm::EDGetTokenT<reco::TrackCollection> track_token;
        edm::EDGetTokenT<reco::VertexCollection>          vtxToken_;
        edm::EDGetTokenT<TrackingParticleCollection>          tpToken_;
        edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
        edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
        edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;
        edm::EDGetTokenT<reco::SimToRecoCollection>          strToken_;

        edm::ParameterSet cfg_;

        TreeWriter tree;

        size gen_pdgid        ;
        size gen_pt           ;
        size gen_eta          ;
        size gen_phi          ;
        size gen_center_eta   ;
        size gen_center_phi   ;
        size gen_delta_eta    ;
        size gen_delta_phi    ;
        size gen_lay_eta      ;
        size gen_lay_phi      ;
        size gen_lay_lay      ;
        size rh_center_eta     ;
        size rh_center_phi     ;
        size rh_center_eta_e   ;
        size rh_center_phi_e   ;
        size rh_delta_eta      ;
        size rh_delta_phi      ;
        size rh_delta_eta_e    ;
        size rh_delta_phi_e    ;
        size rh_lay_eta        ;
        size rh_lay_phi        ;
        size rh_lay_eta_e      ;
        size rh_lay_phi_e      ;
        size rh_lay_lay         ;                                ;
        size track_isMatched    ;
        size track_pt           ;
        size track_eta          ;
        size track_phi          ;
        size track_center_eta   ;
        size track_center_phi   ;
        size track_center_eta_e ;
        size track_center_phi_e ;
        size track_delta_eta    ;
        size track_delta_phi    ;
        size track_delta_eta_e  ;
        size track_delta_phi_e  ;


        const double beginOfDet  = 527;
        const double centerOfDet = 539.5;
        const double endOfDet    = 552;

};


ME0TreeMaker::ME0TreeMaker(const edm::ParameterSet& iConfig) : cfg_(iConfig), tree (iConfig.getUntrackedParameter<std::string>("outFileName"),"Events","")
{
  track_token = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );
  vtxToken_   = consumes<reco::VertexCollection>( edm::InputTag("offlinePrimaryVertices") );
  tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );
  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
  dToken_     = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
//  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag("me0RecHits") );
  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag("newrh") );
  strToken_   = consumes<reco::SimToRecoCollection>( edm::InputTag("trackingParticleRecoTrackAsssociation") );

  gen_pdgid        = tree.add<int>(       "gen_pdgid"     ,"I",0);
  gen_pt           = tree.add<float>(     "gen_pt"        ,"F",0);
  gen_eta          = tree.add<float>(     "gen_eta"       ,"F",0);
  gen_phi          = tree.add<float>(     "gen_phi"       ,"F",0);
  gen_center_eta   = tree.add<float>(     "gen_center_eta","F",-99);
  gen_center_phi   = tree.add<float>(     "gen_center_phi","F",-99);
  gen_delta_eta    = tree.add<float>(     "gen_delta_eta" ,"F",-99);
  gen_delta_phi    = tree.add<float>(     "gen_delta_phi" ,"F",-99);
  gen_lay_eta      = tree.addMulti<float>("gen_lay_eta" , -99);
  gen_lay_phi      = tree.addMulti<float>("gen_lay_phi" , -99);
  gen_lay_lay      = tree.addMulti<size8>("gen_lay_lay" ,  9);
  rh_center_eta    = tree.add<float>(     "rh_center_eta"  ,"F",-99);
  rh_center_phi    = tree.add<float>(     "rh_center_phi"  ,"F",-99);
  rh_center_eta_e  = tree.add<float>(     "rh_center_eta_e","F",-99);
  rh_center_phi_e  = tree.add<float>(     "rh_center_phi_e","F",-99);
  rh_delta_eta     = tree.add<float>(     "rh_delta_eta"   ,"F",-99);
  rh_delta_phi     = tree.add<float>(     "rh_delta_phi"   ,"F",-99);
  rh_delta_eta_e   = tree.add<float>(     "rh_delta_eta_e" ,"F",-99);
  rh_delta_phi_e   = tree.add<float>(     "rh_delta_phi_e" ,"F",-99);
  rh_lay_eta       = tree.addMulti<float>("rh_lay_eta"     ,-99);
  rh_lay_phi       = tree.addMulti<float>("rh_lay_phi"     ,-99);
  rh_lay_eta_e     = tree.addMulti<float>("rh_lay_eta_e"   ,-99);
  rh_lay_phi_e     = tree.addMulti<float>("rh_lay_phi_e"   ,-99);
  rh_lay_lay       = tree.addMulti<size8>("rh_lay_lay"     , 9);
  track_isMatched     = tree.addMulti<size8>(     "track_isMatched"   ,9);
  track_pt            = tree.addMulti<float>(     "track_pt"          ,-99);
  track_eta           = tree.addMulti<float>(     "track_eta"         ,-99);
  track_phi           = tree.addMulti<float>(     "track_phi"         ,-99);
  track_center_eta    = tree.addMulti<float>(     "track_center_eta"  ,-99);
  track_center_phi    = tree.addMulti<float>(     "track_center_phi"  ,-99);
  track_center_eta_e  = tree.addMulti<float>(     "track_center_eta_e",-99);
  track_center_phi_e  = tree.addMulti<float>(     "track_center_phi_e",-99);
  track_delta_eta     = tree.addMulti<float>(     "track_delta_eta"   ,-99);
  track_delta_phi     = tree.addMulti<float>(     "track_delta_phi"   ,-99);
  track_delta_eta_e   = tree.addMulti<float>(     "track_delta_eta_e" ,-99);
  track_delta_phi_e   = tree.addMulti<float>(     "track_delta_phi_e" ,-99);


  tree.book();
}


ME0TreeMaker::~ME0TreeMaker() {tree.write();}





void ME0TreeMaker::getHits(edm::Handle<std::vector<PSimHit> >&  simHitH,edm::Handle<ME0DigiPreRecoCollection>& digis,edm::Handle<ME0RecHitCollection>& rechits,
    const SimTrack * simTrack, HitCollection& hits
){
  hits.clear();

  //find GEM Hits...should only work for the primary interaction
  for (const auto& simHit : (*simHitH)) {
    if(simHit.eventId() != simTrack->eventId()) continue;
    if(simHit.trackId() != simTrack->trackId()) continue;
    const ME0DetId detID = simHit.detUnitId();
    const ME0DigiPreRecoCollection::Range& range = digis->get(detID);

    unsigned int iD = 0;
    double minDR2 = -1;
    const ME0DigiPreReco * digi = 0;
    unsigned int digiInd = 0;
    for (ME0DigiPreRecoCollection::const_iterator idigi = range.first;
        idigi != range.second;idigi++) {
      double dr2 = (simHit.entryPoint().x() - idigi->x())*(simHit.entryPoint().x() - idigi->x()) + (simHit.entryPoint().y() - idigi->y())*(simHit.entryPoint().y() - idigi->y());
      if(minDR2 < 0 || dr2 < minDR2){
        digi = &(*idigi);
        digiInd = iD;
        minDR2 = dr2;
      }
      iD++;
    }
    if(minDR2 < 0 ) {std::cout << "Could not find matching digi!!!" <<std::endl;}

    auto recHitRange = rechits->get(detID);
    const ME0RecHit * rechit = 0;
    if(digi && (digiInd < recHitRange.second -  recHitRange.first)) rechit = &(*(recHitRange.first + digiInd));
    if(!rechit) {std::cout << "Could not find matching rechit!!!" <<std::endl;}
//    std::cout <<"("<< simHit.entryPoint().x() <<","<<simHit.entryPoint().y()<<","<<simHit.entryPoint().z() <<") ("<< digi->x() <<","<<digi->y() <<") ("<< rechit->localPosition().x() <<","<<rechit->localPosition().y()<<","<<rechit->localPosition().z() <<")"<<std::endl;
    hits.add(detID,&simHit,digi, rechit);
  }
  hits.sort();
}

void ME0TreeMaker::getSimInfo(const ME0Geometry* mgeom, const HitCollection& hits, float& gen_eta, float& gen_phi, float& gen_delta_eta, float& gen_delta_phi){
  gen_eta = -99;
  gen_phi = -99;
  gen_delta_eta = -99;
  gen_delta_phi = -99;
  if(hits.simHits.size() < 2) return;

  auto getCenter = [](double v1,double v2,double z1,double z2,double zc)->double{
    if(z2 == z1) return -99;
    double m = (v2 - v1)/(z2 - z1);
    double b =  (z2*v1 - z1*v2)/(z2 - z1);
    return m*zc+b;
  };
  auto getGlobal = [&](int ind) -> GlobalPoint {
    return mgeom->etaPartition(hits.detIDs[ind])->toGlobal(hits.simHits[ind]->entryPoint());
  };
  const int upInd = floor(hits.simHits.size()/2);

  //Do center first
  GlobalPoint centerUp = getGlobal(upInd);
  GlobalPoint centerDown = getGlobal(upInd-1);
  const double zc = centerUp.z() > 0 ? centerOfDet : -1*centerOfDet;
auto params = mgeom->etaPartition(hits.detIDs[upInd])->specs()->parameters();
//  std::cout << params[0] <<" "<< params[2] << std::endl;
//  std::cout <<mgeom->etaPartition(hits.detIDs[upInd-1])->toGlobal(LocalPoint(0,0,0)).phi() - mgeom->etaPartition(hits.detIDs[upInd])->toGlobal(LocalPoint(0,0,0)).phi()<<std::endl;

  double xc = getCenter(centerUp.x(),centerDown.x(),centerUp.z(),centerDown.z(),zc);
  double yc = getCenter(centerUp.y(),centerDown.y(),centerUp.z(),centerDown.z(),zc);
  GlobalPoint centerPt(xc,yc,zc);
  gen_eta = centerPt.eta();
  gen_phi = centerPt.phi();

  //now for deltas
  GlobalPoint up = getGlobal(hits.simHits.size()-1);
  GlobalPoint down = getGlobal(0);

  gen_delta_eta = up.eta() - down.eta();
  gen_delta_phi = deltaPhi(up.phi(), down.phi());
}


std::pair<ME0DetId,ME0Segment*> ME0TreeMaker::buildSegment(const ME0Geometry* mgeom, const HitCollection& hits) const
{
  std::vector<ME0DetId>         detIDs; detIDs.reserve(hits.recHits.size());
  std::vector<const ME0RecHit*> rechits; rechits.reserve(hits.recHits.size());
  for(unsigned int iH = 0; iH < hits.recHits.size(); ++iH) {
    detIDs.push_back(hits.detIDs[iH]);
    rechits.push_back(hits.recHits[iH]);
  }
  if(rechits.size() < 2) return std::pair<ME0DetId,ME0Segment*>(ME0DetId(),0);

  MuonSegFit::MuonRecHitContainer muonRecHits;

//  uint32_t refid = ensemble.second.begin()->first;
  const ME0EtaPartition * refPart = mgeom->etaPartition(detIDs[0]);
  // select hits from the ensemble and sort it
  for (unsigned int rh=0; rh < rechits.size();rh++){
    // for segFit - using local point in first partition frame
    const ME0EtaPartition * thePartition   =   mgeom->etaPartition(detIDs[rh]);
    GlobalPoint gp = thePartition->toGlobal(rechits[rh]->localPosition());
    const LocalPoint lp = refPart->toLocal(gp);
    ME0RecHit *newRH = rechits[rh]->clone();
    LocalError le(newRH->localPositionError().xx() <= 0 ? 0.01 : newRH->localPositionError().xx(),
        newRH->localPositionError().yy() <= 0 ? 0.01 : newRH->localPositionError().yy(),
            newRH->localPositionError().xy() <= 0 ? 0.01 : newRH->localPositionError().xy());

    newRH->setError(le);
    newRH->setPosition(lp);

//    std::cout<< "("<< lp.x() <<" "<<lp.y() <<","<<lp.z()<<") ";

    MuonSegFit::MuonRecHitPtr trkRecHit(newRH);
    muonRecHits.push_back(trkRecHit);
  }

  // The actual fit on all hits of the vector of the selected Tracking RecHits:
  MuonSegFit  * sfit_ = new MuonSegFit(muonRecHits);
  sfit_->fit();

  // obtain all information necessary to make the segment:
  LocalPoint protoIntercept      = sfit_->intercept();
  LocalVector protoDirection     = sfit_->localdir();
  AlgebraicSymMatrix protoErrors = sfit_->covarianceMatrix();
  double protoChi2               = sfit_->chi2();

//  std::cout<< " --> ("<< protoIntercept.x() <<" "<<protoIntercept.y() <<","<<protoIntercept.z()<<") ";


  // Calculate the central value and uncertainty of the segment time
  float averageTime=0.;
  for (auto rh=rechits.begin(); rh!=rechits.end(); ++rh){
    averageTime += (*rh)->tof();
  }
  if(rechits.size() != 0) averageTime=averageTime/(rechits.size());
  float timeUncrt=0.;
  for (auto rh=rechits.begin(); rh!=rechits.end(); ++rh){
    timeUncrt += pow((*rh)->tof()-averageTime,2);
  }
  if(rechits.size() > 1) timeUncrt=timeUncrt/(rechits.size()-1);
  timeUncrt = sqrt(timeUncrt);

  for (auto rh:muonRecHits) rh.reset();
  delete sfit_;

  return std::pair<ME0DetId,ME0Segment*>(detIDs[0], new ME0Segment(rechits, protoIntercept, protoDirection, protoErrors, protoChi2, averageTime, timeUncrt));
}

void ME0TreeMaker::getSegmentProperties(const ME0Segment * seg, const ME0EtaPartition* partition,
    float& cen_eta, float& cen_phi, float& delta_eta, float& delta_phi
)const{

  //get "local" z values
  auto layCen =  partition->toGlobal(LocalPoint(0,0,0));
//  double zSign = layCen.z() > 0 ? 1.0 : -1.0;

  auto segLP                   = seg->localPosition();
  auto segLD                   = seg->localDirection();

  auto extrap = [&] (double extZ) -> LocalPoint {
    double extX = segLP.x()+segLD.x()*extZ/segLD.z();
    double extY = segLP.y()+segLD.y()*extZ/segLD.z();
    return LocalPoint(extX,extY,extZ);

  };



  auto locLow  = extrap(beginOfDet - layCen.z());
  auto lochigh = extrap(endOfDet - layCen.z());
  auto loccen  = extrap(centerOfDet - layCen.z());

  auto globlow  = partition->toGlobal(locLow );
  auto globhigh = partition->toGlobal(lochigh);
  auto globcen  = partition->toGlobal(loccen );

  cen_eta   = globcen.eta();
  cen_phi   = globcen.phi();
  delta_eta = globhigh.eta()-globlow.eta();
  delta_phi = deltaPhi(globhigh.phi(), globlow.phi());
}
FreeTrajectoryState
ME0TreeMaker::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix55& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);

  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0TreeMaker::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix66& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);

  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}
void ME0TreeMaker::getFromFTS(const FreeTrajectoryState& fts,
            GlobalVector& p3, GlobalVector& r3,
            int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());

  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}


GlobalPoint ME0TreeMaker::getPropogatedTrack(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, double extZ){
  float zSign = track->pz() > 0 ? 1.0f : -1.0f;
  const float zValue = extZ * zSign;
  Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
  int chargeReco = track->charge();
  GlobalVector p3reco, r3reco;
  p3reco = GlobalVector(track->outerPx(), track->outerPy(), track->outerPz());
  r3reco = GlobalVector(track->outerX(), track->outerY(), track->outerZ());

  AlgebraicSymMatrix66 covReco;
  //This is to fill the cov matrix correctly
  AlgebraicSymMatrix55 covReco_curv;
  covReco_curv = track->outerStateCovariance();
  FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
  getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);
  //Now we propagate and get the propagated variables from the propagated state
  TrajectoryStateOnSurface lastrecostate;
  lastrecostate = ThisshProp->propagate(initrecostate,*plane);
  if (!lastrecostate.isValid()) return GlobalPoint();

  FreeTrajectoryState finalrecostate(*lastrecostate.freeTrajectoryState());
  AlgebraicSymMatrix66 covFinalReco;
  GlobalVector p3FinalReco_glob, r3FinalReco_globv;
  getFromFTS(finalrecostate, p3FinalReco_glob, r3FinalReco_globv, chargeReco, covFinalReco);

  return GlobalPoint(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z());
}



void ME0TreeMaker::fillTrackProperties(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, bool isMatched){

    tree.fillMulti( track_isMatched   , size8(isMatched));
    tree.fillMulti( track_pt          , float(track->pt()));
    tree.fillMulti( track_eta         , float(track->eta()));
    tree.fillMulti( track_phi         , float(track->phi()));

    GlobalPoint centerPoint = getPropogatedTrack(bField,ThisshProp,track,centerOfDet);
    GlobalPoint lowPoint = getPropogatedTrack(bField,ThisshProp,track,beginOfDet);
    GlobalPoint upPoint = getPropogatedTrack(bField,ThisshProp,track,endOfDet);

    tree.fillMulti( track_center_eta  , float(centerPoint.eta()));
    tree.fillMulti( track_center_phi  , float(centerPoint.phi()));
//    tree.fillMulti( track_center_eta_e, float(0));
//    tree.fillMulti( track_center_phi_e, float(0));
    tree.fillMulti( track_delta_eta   , float(upPoint.eta() - lowPoint.eta()));
    tree.fillMulti( track_delta_phi   , float(deltaPhi(upPoint.phi(), lowPoint.phi())));
//    tree.fillMulti( track_delta_eta_e , float(0));
//    tree.fillMulti( track_delta_phi_e , float(0));
}


void
ME0TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> me0g;
  iSetup.get<MuonGeometryRecord>().get(me0g);
  const ME0Geometry* mgeom = &*me0g;

//
//  std::cout << " Geometry node for GEMGeom is  " << &(*me0g) << std::endl;
//  std::cout << " detTypes       \t"              <<me0g->detTypes().size() << std::endl;
//  std::cout << " GeomDetUnit       \t"           <<me0g->detUnits().size() << std::endl;
//  std::cout << " GeomDet           \t"           <<me0g->dets().size() << std::endl;
//  std::cout << " GeomDetUnit DetIds\t"           <<me0g->detUnitIds().size() << std::endl;
//  std::cout << " eta partitions \t"              <<me0g->etaPartitions().size() << std::endl;
//  std::cout << " chambers       \t"              <<me0g->chambers().size() << std::endl;
////  std::cout << " super chambers  \t"             <<me0g->superChambers().size() << std::endl;
////  std::cout << " rings  \t\t"                    <<me0g->rings().size() << std::endl;
////  std::cout << " stations  \t\t"                 <<me0g->stations().size() << std::endl;
////  std::cout << " regions  \t\t"                  <<me0g->regions().size() << std::endl;
//
//  for(unsigned int iE = 0; iE < me0g->etaPartitions().size(); ++iE){
//    if( me0g->etaPartitions()[iE]->id().layer() != 1) continue;
//    std::cout << me0g->etaPartitions()[iE]->id().chamber() <<" "<< me0g->etaPartitions()[iE]->toGlobal(LocalPoint(0,0,0)).eta() <<" "<<me0g->etaPartitions()[iE]->toGlobal(LocalPoint(0,0,0)).phi()<<" "<<me0g->etaPartitions()[iE]->toGlobal(LocalPoint(0,0,0)).z()<<std::endl;
//  }


  TrackingParticle r;
  edm::Handle <reco::TrackCollection > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle <reco::VertexCollection > verts;
  iEvent.getByToken(vtxToken_,verts);

  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByToken(tpToken_,TPCollectionH);

  edm::Handle<std::vector<PSimHit> >  simHitH ;
  iEvent.getByToken(shToken_,simHitH);

  edm::Handle<ME0DigiPreRecoCollection> digisH;
  iEvent.getByToken(dToken_,digisH);

  edm::Handle<ME0RecHitCollection> rechitsH;
  iEvent.getByToken(rhToken_,rechitsH);

  edm::Handle<reco::SimToRecoCollection> simRecH;
  iEvent.getByToken(strToken_,simRecH);

  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  edm::ESHandle<Propagator> ThisshProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", ThisshProp);


  for (TrackingParticleCollection::size_type i=0; i<TPCollectionH->size(); i++) {
    TrackingParticleRef trpart(TPCollectionH, i);
    if(!trpart->genParticles().size() ||  !trpart->g4Tracks().size())continue;
    if(TMath::Abs(trpart->pdgId()) != 13) continue;

    const SimTrack * simTrack = 0;
    for(unsigned int iST = 0; iST <  trpart->g4Tracks().size(); ++iST){
      if(TMath::Abs(trpart->g4Tracks()[iST].type()) != 13) continue;
      if(simTrack) {std::cout << "FOUND AN EXTRA MUON SIM TRACK!!!" <<std::endl; continue;}
      simTrack =&trpart->g4Tracks()[iST];
      }
    if(!simTrack) continue;
    tree.reset();

    tree.fill(gen_pdgid      ,int(simTrack->type()));
    tree.fill(gen_pt         ,float(simTrack->momentum().pt()));
    tree.fill(gen_eta        ,float(simTrack->momentum().eta()));
    tree.fill(gen_phi        ,float(simTrack->momentum().phi()));

    HitCollection hits;
    getHits(simHitH,digisH,rechitsH,simTrack,hits);
    float g_eta; float g_phi; float g_delta_eta; float g_delta_phi;
    getSimInfo(mgeom,hits,g_eta,g_phi,g_delta_eta,g_delta_phi);

    tree.fill(gen_center_eta       ,g_eta);
    tree.fill(gen_center_phi       ,g_phi);
    tree.fill(gen_delta_eta        ,g_delta_eta);
    tree.fill(gen_delta_phi        ,g_delta_phi);

    for(unsigned int iH = 0; iH < hits.simHits.size(); ++iH){
      if(!hits.simHits[iH]) continue;
      float eta = mgeom->etaPartition(hits.detIDs[iH])->toGlobal(hits.simHits[iH]->entryPoint()).eta();
      float phi = mgeom->etaPartition(hits.detIDs[iH])->toGlobal(hits.simHits[iH]->entryPoint()).phi();
      tree.fillMulti(gen_lay_lay,size8(hits.detIDs[iH].layer()));
      tree.fillMulti(gen_lay_eta,eta);
      tree.fillMulti(gen_lay_phi,phi);
    }
    float s_eta      = -99;
    float s_phi      = -99;
    float s_delta_eta= -99;
    float s_delta_phi= -99;
    std::pair<ME0DetId,ME0Segment*>  seg = buildSegment(mgeom,hits);
    if(seg.second)
      getSegmentProperties(seg.second,mgeom->etaPartition(seg.first),s_eta, s_phi, s_delta_eta, s_delta_phi);
    tree.fill(rh_center_eta       ,s_eta);
    tree.fill(rh_center_phi       ,s_phi);
    tree.fill(rh_delta_eta        ,s_delta_eta);
    tree.fill(rh_delta_phi        ,s_delta_phi);

    for(unsigned int iH = 0; iH < hits.recHits.size(); ++iH){
      if(!hits.recHits[iH]) continue;
      float eta = mgeom->etaPartition(hits.detIDs[iH])->toGlobal(hits.recHits[iH]->localPosition()).eta();
      float phi = mgeom->etaPartition(hits.detIDs[iH])->toGlobal(hits.recHits[iH]->localPosition()).phi();
      tree.fillMulti(rh_lay_lay,size8(hits.detIDs[iH].layer()));
      tree.fillMulti(rh_lay_eta,eta);
      tree.fillMulti(rh_lay_phi,phi);
    }
    delete seg.second;

    edm::RefToBase<reco::Track> track;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > simRecAsso;
    if(simRecH->find(trpart) != simRecH->end()) {
          simRecAsso = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) (*simRecH)[trpart];
          if(simRecAsso.begin() != simRecAsso.end() ) track = simRecAsso.begin()->first;

    }
    if(track.isNonnull()){

      GlobalPoint centerPoint = getPropogatedTrack(bField,ThisshProp,track.get(),centerOfDet);
      GlobalPoint lowPoint = getPropogatedTrack(bField,ThisshProp,track.get(),beginOfDet);
      GlobalPoint upPoint = getPropogatedTrack(bField,ThisshProp,track.get(),endOfDet);
      tree.fillMulti( track_isMatched   , size8(true));
      tree.fillMulti( track_pt          , float(track->pt()));
      tree.fillMulti( track_eta         , float(track->eta()));
      tree.fillMulti( track_phi         , float(track->phi()));
      tree.fillMulti( track_center_eta  , float(centerPoint.eta()));
      tree.fillMulti( track_center_phi  , float(centerPoint.phi()));
      tree.fillMulti( track_delta_eta   , float(upPoint.eta() - lowPoint.eta()));
      tree.fillMulti( track_delta_phi   , float(deltaPhi(upPoint.phi(), lowPoint.phi())));
    }

    //Now fill all other tracks
    for(const auto& oTrack : *tracks){
      if(&oTrack == track.get()) continue;

      GlobalPoint centerPoint = getPropogatedTrack(bField,ThisshProp,&oTrack,centerOfDet);
      bool fill = false;
      if(TMath::Abs(centerPoint.eta() - g_eta) < 0.35 && TMath::Abs(deltaPhi(centerPoint.phi(), g_phi)) < 0.35    )fill =true;
      if(seg.second)
        if(TMath::Abs(centerPoint.eta() - s_eta) < 0.35 && TMath::Abs(deltaPhi(centerPoint.phi(), s_phi)) < 0.35    )fill =true;
      if(!fill) continue;

      GlobalPoint lowPoint = getPropogatedTrack(bField,ThisshProp,&oTrack,beginOfDet);
      GlobalPoint upPoint = getPropogatedTrack(bField,ThisshProp,&oTrack,endOfDet);
      tree.fillMulti( track_isMatched   , size8(false));
      tree.fillMulti( track_pt          , float(oTrack.pt()));
      tree.fillMulti( track_eta         , float(oTrack.eta()));
      tree.fillMulti( track_phi         , float(oTrack.phi()));
      tree.fillMulti( track_center_eta  , float(centerPoint.eta()));
      tree.fillMulti( track_center_phi  , float(centerPoint.phi()));
      tree.fillMulti( track_delta_eta   , float(upPoint.eta() - lowPoint.eta()));
      tree.fillMulti( track_delta_phi   , float(deltaPhi(upPoint.phi(), lowPoint.phi())));
    }



    tree.fillTree();
    }
}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0TreeMaker);
