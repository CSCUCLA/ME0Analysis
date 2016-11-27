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





class ME0TestFakeRate : public edm::EDAnalyzer {
  double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
  double deltaR2(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
  }



    public:
        explicit ME0TestFakeRate(const edm::ParameterSet&);
        ~ME0TestFakeRate();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);

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

        TString outFileName;
        HistGetter hists;



        const double beginOfDet  = 527;
        const double centerOfDet = 539.5;
        const double endOfDet    = 552;

};


ME0TestFakeRate::ME0TestFakeRate(const edm::ParameterSet& iConfig) : cfg_(iConfig),outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
  track_token = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );
  dToken_     = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );


}


ME0TestFakeRate::~ME0TestFakeRate() {  hists.write(outFileName);}



FreeTrajectoryState
ME0TestFakeRate::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix55& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);

  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0TestFakeRate::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix66& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);

  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}
void ME0TestFakeRate::getFromFTS(const FreeTrajectoryState& fts,
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


GlobalPoint ME0TestFakeRate::getPropogatedTrack(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, double extZ){
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



void
ME0TestFakeRate::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> me0g;
  iSetup.get<MuonGeometryRecord>().get(me0g);
  const ME0Geometry* mgeom = &*me0g;


  TrackingParticle r;
  edm::Handle <reco::TrackCollection > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle<ME0DigiPreRecoCollection> digisH;
  iEvent.getByToken(dToken_,digisH);

  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  edm::ESHandle<Propagator> ThisshProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", ThisshProp);

  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);

//  hists.getOrMake1D("nMuons",";# of muons",1,0,2)->Fill(1.0);
//  hists.getOrMake1D("nMuons",";# of muons",1,0,2)->Fill(1.0);

  int nMuons = 0;
  int nMatch = 0;
  int nMatch_3GeV = 0;
  int nMatch_tight = 0;
  int nMatch_tight_3GeV = 0;

  for(const auto& d: (*digisH) ){
   if(d.first.layer() != 1) continue;
   const auto* epart = mgeom->etaPartition(d.first);
   for (ME0DigiPreRecoCollection::const_iterator idigi = d.second.first;
       idigi != d.second.second;idigi++) {
     if(TMath::Abs(idigi->pdgid()) != 13) continue;
     nMuons++;

     bool match            = 0;
     bool match_3GeV       = 0;
     bool match_tight      = 0;
     bool match_tight_3GeV = 0;

     auto gl = epart->toGlobal(LocalPoint(idigi->x(),idigi->y(),0));
     for (std::vector<reco::Track>::const_iterator iTrack = tracks->begin();
         iTrack != tracks->end(); ++iTrack){
       GlobalPoint trackPoint = getPropogatedTrack(bField,ThisshProp,&(*iTrack),gl.z());
       bool match_this = true;
       bool matchTight_this = true;
       if(TMath::Abs(trackPoint.eta() - gl.eta() ) > 0.04) match_this = false;
       if(TMath::Abs(deltaPhi(trackPoint.phi(),gl.phi() ) )> 0.04) match_this = false;

       if(TMath::Abs(trackPoint.eta() - gl.eta() ) > 0.03) matchTight_this = false;
       if(TMath::Abs(deltaPhi(trackPoint.phi(),gl.phi() ) )> 0.03) matchTight_this = false;
       if(match_this) match = true;
       if(match_this && iTrack->pt() > 3) match_3GeV = true;
       if(matchTight_this) match_tight = true;
       if(matchTight_this && iTrack->pt() > 3) match_tight_3GeV = true;
     }

     if(match           ) nMatch           +=1;
     if(match_3GeV      ) nMatch_3GeV      +=1;
     if(match_tight     ) nMatch_tight     +=1;
     if(match_tight_3GeV) nMatch_tight_3GeV+=1;

   }
  }
  auto fill = [&](TString name, int number) {
    for(int iM = 0; iM < number; ++iM) hists.getOrMake1D(name,";# of muons",1,0,2)->Fill(1.0);
    hists.getOrMake1D(TString::Format("%s_perE",name.Data()),";# of muons",20,-.5,19.5)->Fill(number);
  };

  fill("nMuons",nMuons);
  fill("nMuons_match",nMatch);
  fill("nMuons_match_3GeV",nMatch_3GeV);
  fill("nMuons_match_tight",nMatch_tight);
  fill("nMuons_match_tight_3GeV",nMatch_tight_3GeV);
}


//define this as a plug-in
DEFINE_FWK_MODULE(ME0TestFakeRate);
