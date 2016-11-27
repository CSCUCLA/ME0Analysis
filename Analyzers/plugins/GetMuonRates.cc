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


#include "../../ME0Analysis/plugins/HistGetter.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "../interface/ParticleUtilities.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class GetMuonRates : public edm::EDAnalyzer {
    public:
        explicit GetMuonRates(const edm::ParameterSet&);
        ~GetMuonRates();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<reco::GenParticleCollection> gen_token;
        edm::EDGetTokenT<pat::MuonCollection>            muonToken_;
        edm::EDGetTokenT<reco::VertexCollection>          vtxToken_;
        edm::EDGetTokenT<double>                         rhoToken_;
        edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken_;

        TString outFileName;
        HistGetter hists;
};


GetMuonRates::GetMuonRates(const edm::ParameterSet& iConfig)
   : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
  gen_token = consumes<reco::GenParticleCollection>( edm::InputTag("prunedGenParticles") );
  muonToken_ = consumes<pat::MuonCollection>( edm::InputTag("slimmedMuons") );
  vtxToken_ = consumes<reco::VertexCollection>( edm::InputTag("offlineSlimmedPrimaryVertices") );
  rhoToken_  =consumes<double>                        (edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
  genEvtInfoToken_  = consumes<GenEventInfoProduct>   (edm::InputTag("generator"));

}


GetMuonRates::~GetMuonRates() {
  hists.write(outFileName);
}


typedef std::map<unsigned int, std::vector<std::pair<const PSimHit*,ME0DetId> > > SimHCont;

void
GetMuonRates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  edm::Handle <reco::GenParticleCollection > genParts;
  iEvent.getByToken(gen_token,genParts);

  edm::Handle <pat::MuonCollection > muons;
  iEvent.getByToken(muonToken_,muons);

  edm::Handle <reco::VertexCollection > verts;
  iEvent.getByToken(vtxToken_,verts);

  edm::Handle<double>     rho_;
  iEvent.getByToken(rhoToken_,rho_);

  edm::Handle<GenEventInfoProduct>             genEvtInfo_;
  iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);

  if(verts->size() == 0 ) return;

  bool hasgoodvtx = false;
  if(!(*verts)[0].isFake() &&
      (*verts)[0].ndof() > 4. &&
      (*verts)[0].position().Rho() <= 2.0 &&
      fabs((*verts)[0].position().Z())<=24.0)
    hasgoodvtx = true;
  if(!hasgoodvtx) return;


  std::vector<std::pair<reco::GenParticleRef,reco::GenParticleRef> > trueLeptonList;
  for(unsigned int iG = 0; iG < genParts->size(); ++iG){
    const auto& gp = genParts->at(iG);
    if(!ParticleInfo::isFinal(gp.status())) continue; //Final particle from pythia
    const int pdgid = TMath::Abs(gp.pdgId());
    if(pdgid != ParticleInfo::p_eminus && pdgid != ParticleInfo::p_muminus ) continue; //electron or muon
    //get first version of particle
//    reco::GenParticleRef gpr(genH, iG);
    auto ogp = ParticleUtilities::getOriginal(reco::GenParticleRef(genParts, iG),genParts);
    if(ogp->numberOfMothers() == 1 && TMath::Abs(ogp->mother(0)->pdgId()) == ParticleInfo::p_Z0 ){
      trueLeptonList.emplace_back( ogp,  reco::GenParticleRef( genParts,ogp->motherRef(0).key()));
    } else {
      trueLeptonList.emplace_back( ogp, ogp);
    }
  }

  int nZMuons = 0;

  for(unsigned int iP = 0; iP < trueLeptonList.size(); ++iP){
    const auto& genInfo = trueLeptonList[iP];
    if(genInfo.first.key() == genInfo.second.key()) continue;
    if(TMath::Abs(genInfo.first->pdgId()) != ParticleInfo::p_muminus) continue;
    nZMuons++;
  }
  if(nZMuons != 2) return;

  int nLoose =0;
  int nTight = 0;
  int nIsoTight = 0;

  for(const auto& muon : (*muons)){
    if(muon.pt() < 5 || TMath::Abs(muon.eta()) > 2.4 ) continue;
    if(!((muon.isGlobalMuon() || (muon.isTrackerMuon() && muon.numberOfMatches() > 0 ) ) &&
        muon.muonBestTrackType() != 2 ) ) continue;

    double d0 = muon.innerTrack().isNonnull() ? -1.*muon.bestTrack()->dxy(verts->at(0).position()):0;
    double dz = muon.innerTrack().isNonnull() ? muon.bestTrack()->dz(verts->at(0).position()):0;
    if(TMath::Abs(d0) > 0.5) continue;
    if(TMath::Abs(dz) > 1) continue;
    nLoose++;
    hists.getOrMake1D("loose_pt",";p_{T} [GeV]",500,-.5,499.5)->Fill(muon.pt(),genEvtInfo_->weight());

    if (muon.pt()<=200.0){ if(! muon.isPFMuon()) continue;}
    else {
        if ( ! (muon.numberOfMatchedStations() > 1
                  && (muon.muonBestTrack()->ptError()/muon.muonBestTrack()->pt()) < 0.3
                  && std::abs(d0) < 0.2
                  && std::abs(dz) < 0.5
                  && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
                  && muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) || muon.isPFMuon() ) continue;
    }
    nTight++;
    hists.getOrMake1D("tight_pt",";p_{T} [GeV]",500,-.5,499.5)->Fill(muon.pt(),genEvtInfo_->weight());

    auto getIso = [&] (const pat::Muon& muon) -> double {
              double PUCorr = 0.5*muon.pfIsolationR03().sumPUPt;
              double isoCH = muon.pfIsolationR03().sumChargedHadronPt;
              double isoNH = muon.pfIsolationR03().sumNeutralHadronEt;
              double isoPhot = muon.pfIsolationR03().sumPhotonEt;
              double iso = (isoCH+std::max(isoPhot+isoNH-PUCorr,0.0))/muon.pt();
              return iso;
    };
    if(getIso(muon) > 0.35) continue;
    nIsoTight++;
    hists.getOrMake1D("isotight_pt",";p_{T} [GeV]",500,-.5,499.5)->Fill(muon.pt(),genEvtInfo_->weight());
  }




  hists.getOrMake1D("nEventsRaw",";# of events",1,0,2)->Fill(1.0);
  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0,genEvtInfo_->weight());
  hists.getOrMake1D("nLoose",";# of loose muons",6,-0.5,5.5)->Fill(nLoose,genEvtInfo_->weight());
  hists.getOrMake1D("nTight",";# of tight muons",6,-0.5,5.5)->Fill(nTight,genEvtInfo_->weight());
  hists.getOrMake1D("nTightISO",";# of isolated tight muons",6,-0.5,5.5)->Fill(nIsoTight,genEvtInfo_->weight());


}



//define this as a plug-in
DEFINE_FWK_MODULE(GetMuonRates);
