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

#include "../../AnalysisSupport/interface/HistGetter.h"
#include "../../AnalysisSupport/TreeWriter.h"
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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "../../AnalysisSupport/interface/ParticleUtilities.h"

#include <TVector2.h>
#include <TMath.h>





class ME0EstimateSampleFakeRate : public edm::EDAnalyzer {
  double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}
  double deltaR2(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
  }

    public:
        explicit ME0EstimateSampleFakeRate(const edm::ParameterSet&);
        ~ME0EstimateSampleFakeRate();


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

        void classification(const std::vector<std::pair<reco::GenParticleRef,reco::GenParticleRef> >& trueLeptonList,
            const std::vector<int>& matchedTP, edm::Handle <reco::TrackCollection >& tracks, const std::vector<int>& matchedTrack,const bool ambiguous  );
        float getTrackIso(edm::Handle <reco::TrackCollection >& tracks, const unsigned int thisTrackIDX, const double dr);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<reco::TrackCollection> track_token;
        edm::EDGetTokenT<reco::VertexCollection>          vtxToken_;
        edm::EDGetTokenT<TrackingParticleCollection>          tpToken_;
        edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
        edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
        edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;
        edm::EDGetTokenT<reco::SimToRecoCollection>          strToken_;
        edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken_;
        edm::EDGetTokenT<reco::GenParticleCollection>  genToken;

        edm::ParameterSet cfg_;

        TString outFileName;
        HistGetter hists;
        TreeWriter tree;

        size gen_only_type            ;
        size gen_only_pdgid           ;
        size gen_only_pt              ;
        size gen_only_eta             ;
        size gen_only_phi             ;


        size weight              ;
        size is_ambiguous        ;
        size track_pt            ;
        size track_eta           ;
        size track_phi           ;
        size track_iso           ;
        size track_proj_phi      ;
        size track_proj_eta      ;
        size track_proj_deltaPhi ;
        size gen_type            ;
        size gen_pdgid           ;
        size gen_pt              ;
        size gen_eta             ;
        size gen_phi             ;
        size segment_phi         ;
        size segment_eta         ;
        size segment_nHits       ;
        size segment_nEHits      ;
        size segment_deltaPhi    ;




        const double beginOfDet  = 527;
        const double centerOfDet = 539.5;
        const double endOfDet    = 552;





};


ME0EstimateSampleFakeRate::ME0EstimateSampleFakeRate(const edm::ParameterSet& iConfig) : cfg_(iConfig),outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
    tree (outFileName,"Events","")
{
  genToken     = consumes<reco::GenParticleCollection>( edm::InputTag("genParticles") );
  tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );
  strToken_   = consumes<reco::SimToRecoCollection>( edm::InputTag("trackingParticleRecoTrackAsssociation") );
  track_token = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );
  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
  genEvtInfoToken_  = consumes<GenEventInfoProduct>   (edm::InputTag("generator"));



  gen_only_type        = tree.addMulti<int>  (     "gen_only_type"           ,0);
  gen_only_pdgid       = tree.addMulti<int>  (     "gen_only_pdgid"          ,0);
  gen_only_pt          = tree.addMulti<float>(     "gen_only_pt"             ,0);
  gen_only_eta         = tree.addMulti<float>(     "gen_only_eta"            ,0);
  gen_only_phi         = tree.addMulti<float>(     "gen_only_phi"            ,0);

  weight    = tree.add<float> (       "weight"               ,"F",0);
  is_ambiguous           = tree.add<bool> (       "is_ambiguous"               ,"O",false);
  track_pt            = tree.addMulti<float>(     "track_pt"           ,0);
  track_eta           = tree.addMulti<float>(     "track_eta"          ,0);
  track_phi           = tree.addMulti<float>(     "track_phi"          ,0);
  track_iso           = tree.addMulti<float>(     "track_iso"          ,0);
  track_proj_phi      = tree.addMulti<float>(     "track_proj_phi"     ,0);
  track_proj_eta      = tree.addMulti<float>(     "track_proj_eta"     ,0);
  track_proj_deltaPhi = tree.addMulti<float>(     "track_proj_deltaPhi",0);
  gen_type            = tree.addMulti<int>  (     "gen_type"           ,0);
  gen_pdgid           = tree.addMulti<int>  (     "gen_pdgid"          ,0);
  gen_pt              = tree.addMulti<float>(     "gen_pt"             ,0);
  gen_eta             = tree.addMulti<float>(     "gen_eta"            ,0);
  gen_phi             = tree.addMulti<float>(     "gen_phi"            ,0);
  segment_eta         = tree.addMulti<float>(     "segment_eta"        ,0);
  segment_phi         = tree.addMulti<float>(     "segment_phi"        ,0);
  segment_nHits       = tree.addMulti<int>  (     "segment_nHits"      ,0);
  segment_nEHits      = tree.addMulti<int>  (     "segment_nEHits"     ,0);
  segment_deltaPhi    = tree.addMulti<float>(     "segment_deltaPhi"   ,0);


  tree.book();

}


ME0EstimateSampleFakeRate::~ME0EstimateSampleFakeRate() {
  tree.write();
//  hists.write("testHist.root");
}



FreeTrajectoryState
ME0EstimateSampleFakeRate::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix55& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);

  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0EstimateSampleFakeRate::getFTS(const GlobalVector& p3, const GlobalVector& r3,
         int charge, const AlgebraicSymMatrix66& cov,
         const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);

  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}
void ME0EstimateSampleFakeRate::getFromFTS(const FreeTrajectoryState& fts,
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


GlobalPoint ME0EstimateSampleFakeRate::getPropogatedTrack(const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp, const reco::Track * track, double extZ){
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


void ME0EstimateSampleFakeRate::classification(const std::vector<std::pair<reco::GenParticleRef,reco::GenParticleRef> >& trueLeptonList,
    const std::vector<int>& matchedTP, edm::Handle <reco::TrackCollection >& tracks, const std::vector<int>& matchedTrack,const bool ambiguous ) {

    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(0);

    int nZLeptons = 0;
    int nWithTP = 0;
    int nWithMatchedTrack = 0;
    int nWithMatchedTrack3GeV = 0;

    for(unsigned int iT = 0; iT < trueLeptonList.size(); ++iT ){
      if(trueLeptonList[iT].first.key() == trueLeptonList[iT].second.key() ) continue;
      nZLeptons++;
      if(ambiguous) continue;
      if(matchedTP[iT] < 0) continue;
      nWithTP++;
      if(matchedTrack[iT] <= 0) continue;
      nWithMatchedTrack++;
      if(tracks->at(matchedTrack[iT]).pt() < 3 )continue;
      nWithMatchedTrack3GeV++;
    }
    if(nZLeptons != 4) return;
    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(1);
    if(ambiguous) return;
    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(2);
    if(nWithTP != 4) return;
    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(3);
    if(nWithMatchedTrack != 4) return;
    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(4);
    if(nWithMatchedTrack3GeV != 4) return;
    hists.getOrMake1D("Events",";# of muons",6,-.5,5.5)->Fill(5);
}

float ME0EstimateSampleFakeRate::getTrackIso(edm::Handle <reco::TrackCollection >& tracks, const unsigned int thisTrackIDX, const double dr){
  double dr2 = dr*dr;
  double iso = 0;
  const auto& thisTrack = tracks->at(thisTrackIDX);
  for(unsigned int iT = 0 ; iT < tracks->size(); ++iT){
    if(iT == thisTrackIDX) continue;
    const auto& track = tracks->at(iT);
    if(deltaR2(thisTrack.eta(),thisTrack.phi(), track.eta(),track.phi()) > dr2 ) continue;
    iso += track.pt();
  }

  return iso/thisTrack.pt();

}



void
ME0EstimateSampleFakeRate::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::ESHandle<ME0Geometry> me0g;
  iSetup.get<MuonGeometryRecord>().get(me0g);
  const ME0Geometry* mgeom = &*me0g;

  edm::Handle<reco::GenParticleCollection> genH;
  iEvent.getByToken(genToken,genH);

  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByToken(tpToken_,TPCollectionH);

  edm::Handle<reco::SimToRecoCollection> simRecH;
  iEvent.getByToken(strToken_,simRecH);

  edm::Handle <reco::TrackCollection > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle<std::vector<PSimHit> >  simHitH ;
  iEvent.getByToken(shToken_,simHitH);

  edm::ESHandle<MagneticField> bField;
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  edm::ESHandle<Propagator> ThisshProp;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", ThisshProp);

  edm::Handle<GenEventInfoProduct>             genEvtInfo_;
  iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);


  std::vector<std::pair<reco::GenParticleRef,reco::GenParticleRef> > trueLeptonList;

  for(unsigned int iG = 0; iG < genH->size(); ++iG){
    const auto& gp = genH->at(iG);
    if(!ParticleInfo::isFinal(gp.status())) continue; //Final particle from pythia
    const int pdgid = TMath::Abs(gp.pdgId());
    if(pdgid != ParticleInfo::p_eminus && pdgid != ParticleInfo::p_muminus ) continue; //electron or muon
    //get first version of particle
//    reco::GenParticleRef gpr(genH, iG);
    auto ogp = ParticleUtilities::getOriginal(reco::GenParticleRef(genH, iG),genH);
    if(ogp->numberOfMothers() == 1 && TMath::Abs(ogp->mother(0)->pdgId()) == ParticleInfo::p_Z0 ){
      trueLeptonList.emplace_back( ogp,  reco::GenParticleRef( genH,ogp->motherRef(0).key()));
    } else {
      trueLeptonList.emplace_back( ogp, ogp);
    }
  }
//  std::cout << "----START--------" <<std::endl;
//  for(const auto& tl : trueLeptonList){
//    ParticleInfo::printGenParticleInfo(tl.first.get(),tl.first.key());
//    ParticleInfo::printGenParticleInfo(tl.second.get(),tl.second.key());
//    std::cout << "--" <<std::endl;
//  }

  bool ambiguous = false;
  std::vector<int> matchedTP(trueLeptonList.size(),-1);


  for (TrackingParticleCollection::size_type i=0; i<TPCollectionH->size(); i++) {
    TrackingParticleRef trpart(TPCollectionH, i);
    const auto& gps = trpart->genParticles();
    if(!gps.size()) continue;

    int iL = -1;
    int nL = 0;
    for(const auto& gp : gps  ){
      if(TMath::Abs(gp->pdgId()) != ParticleInfo::p_eminus && TMath::Abs(gp->pdgId()) != ParticleInfo::p_muminus  ) continue;
      nL++;
      auto ogp = ParticleUtilities::getOriginal(gp,genH);
      for(unsigned int iT = 0; iT < trueLeptonList.size(); ++iT ){
        if(ogp.key() != trueLeptonList[iT].first.key()) continue;
        if(matchedTP[iT] >= 0) ambiguous = true;
        iL = iT;
        break;
      }
    }
    if(nL > 1) ambiguous = true;
    if(nL == 1) matchedTP[iL] = i;
  }

  std::vector<int> matchedTrack(trueLeptonList.size(),-1);
  std::vector<int> trackToGP(tracks->size(),-1);
  for(unsigned int iT = 0; iT < trueLeptonList.size(); ++iT ){
    if(matchedTP[iT] < 0) continue;
    TrackingParticleRef trpart(TPCollectionH, matchedTP[iT]);

    edm::RefToBase<reco::Track> track;
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > simRecAsso;
    if(simRecH->find(trpart) != simRecH->end()) {
          simRecAsso = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) (*simRecH)[trpart];
          if(simRecAsso.begin() != simRecAsso.end() ) track = simRecAsso.begin()->first;

    }
    if(track.isNonnull()){
      trackToGP[track.key()] = iT;
      matchedTrack[iT] = track.key();
    }

  }

  classification(trueLeptonList,matchedTP,tracks, matchedTrack,ambiguous );





  tree.reset();

  //Gen only filling
  for(unsigned int iP = 0; iP < trueLeptonList.size(); ++iP){
    const auto& genInfo = trueLeptonList[iP];
    tree .fillMulti(gen_only_type         ,int( genInfo.first.key() == genInfo.second.key() ? 0 : genInfo.second.key() + 1));
    tree .fillMulti(gen_only_pdgid        ,int(genInfo.first->pdgId()));
    tree .fillMulti(gen_only_pt           ,float(genInfo.first->pt()));
    tree .fillMulti(gen_only_eta          ,float(genInfo.first->eta()));
    tree .fillMulti(gen_only_phi          ,float(genInfo.first->phi()));
  }


  tree.fill(weight,float(genEvtInfo_->weight()));
  tree.fill(is_ambiguous       ,ambiguous);

  for(unsigned int iT = 0; iT < tracks->size(); ++iT ){

    const auto& track = tracks->at(iT);

    if(trackToGP[iT] < 0 && track.pt() < 0.5) continue;

    tree.fillMulti(track_pt ,float(track.pt()));
    tree.fillMulti(track_eta,float(track.eta()));
    tree.fillMulti(track_phi,float(track.phi()));
    tree.fillMulti(track_iso,getTrackIso(tracks,iT,0.3));

    GlobalPoint centerPoint = getPropogatedTrack(bField,ThisshProp,&track,centerOfDet);
    GlobalPoint lowPoint = getPropogatedTrack(bField,ThisshProp,&track,beginOfDet);
    GlobalPoint upPoint = getPropogatedTrack(bField,ThisshProp,&track,endOfDet);

    tree.fillMulti( track_proj_eta        , centerPoint.z() != 0 ? float(centerPoint.eta())                       : float(-99));
    tree.fillMulti( track_proj_phi        , centerPoint.z() != 0 ? float(centerPoint.phi())                       : float(-99));
    tree.fillMulti( track_proj_deltaPhi   , upPoint.z() != 0  && lowPoint.z() != 0 ? float(deltaPhi(upPoint.phi(), lowPoint.phi())) : float(-99));

//    if(false){
    if(trackToGP[iT] >= 0 ){

      const auto& genInfo = trueLeptonList[trackToGP[iT]];
      const auto& tp = TPCollectionH->at(matchedTP[trackToGP[iT]]);

      tree .fillMulti(gen_type         ,int( genInfo.first.key() == genInfo.second.key() ? 0 : genInfo.second.key() + 1));
      tree .fillMulti(gen_pdgid        ,int(genInfo.first->pdgId()));
      tree .fillMulti(gen_pt           ,float(genInfo.first->pt()));
      tree .fillMulti(gen_eta          ,float(genInfo.first->eta()));
      tree .fillMulti(gen_phi          ,float(genInfo.first->phi()));

      float seg_phi = 0;
      float seg_eta = 0;
      int  seg_nHits = 0;
      int  seg_nEHits = 0;
      float seg_deltaPhi = 0;
//      std::cout << genInfo.first->pdgId()<<" "<< genInfo.first.key() <<" "<< genInfo.first->pt() <<" "<< genInfo.first->eta() <<" "<< genInfo.first->phi() << std::endl;
      const SimTrack * simTrack = 0;
      if( genInfo.first->pdgId() == ParticleInfo::p_muminus)
      for(unsigned int iST = 0; iST <  tp.g4Tracks().size(); ++iST){
//        std::cout << tp.g4Tracks()[iST].type()  <<" "<<tp.g4Tracks()[iST].genpartIndex() <<" "<< std::endl;
        if(tp.g4Tracks()[iST].type() != genInfo.first->pdgId()   ) continue;
        if(simTrack) {std::cout << "FOUND AN EXTRA MUON SIM TRACK!!!" <<std::endl; continue;}
        simTrack =&tp.g4Tracks()[iST];
        }

      std::vector<std::vector<const PSimHit*> > laysHit(6);



      if(simTrack){


        for (const auto& simHit : (*simHitH)) {
          if(simHit.eventId() != simTrack->eventId()) continue;
          if(simHit.trackId() != simTrack->trackId()) continue;
          const ME0DetId detID = simHit.detUnitId();
          laysHit[detID.layer()  - 1].push_back(&simHit);
        }
        std::vector<const PSimHit *> prunedHits;
        for(auto& l : laysHit){
          if(l.size() >0) {
            seg_nHits++;
            prunedHits.push_back(l.front());
          }
          if(l.size() > 1) seg_nEHits += l.size() -1;
        }

        if(seg_nHits >= 3){

          auto getGlobal = [&](const PSimHit * sh) -> GlobalPoint {
            return mgeom->etaPartition(sh->detUnitId())->toGlobal(sh->entryPoint());
          };
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
          seg_eta = centerPt.eta();
          seg_phi = centerPt.phi();
          seg_deltaPhi = deltaPhi(upPt.phi(), downPt.phi());
        }

      }
      tree .fillMulti(segment_phi      ,seg_phi);
      tree .fillMulti(segment_eta      ,seg_eta);
      tree .fillMulti(segment_nHits    ,seg_nHits);
      tree .fillMulti(segment_nEHits    ,seg_nEHits);
      tree .fillMulti(segment_deltaPhi ,seg_deltaPhi);


    } else {
      tree .fillMulti(gen_type         ,int(0));
      tree .fillMulti(gen_pdgid        ,int(0));
      tree .fillMulti(gen_pt           ,float(0));
      tree .fillMulti(gen_eta          ,float(0));
      tree .fillMulti(gen_phi          ,float(0));
      tree .fillMulti(segment_phi      ,float(0));
      tree .fillMulti(segment_eta      ,float(0));
      tree .fillMulti(segment_nHits    ,int(0));
      tree .fillMulti(segment_nEHits   ,int(0));
      tree .fillMulti(segment_deltaPhi ,float(0));
    }


  }
  tree.fillTree();






}


//define this as a plug-in
DEFINE_FWK_MODULE(ME0EstimateSampleFakeRate);
