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
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"


#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "../interface/ME0Helper.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class GERHAnalyzer : public edm::EDAnalyzer {
    public:
        explicit GERHAnalyzer(const edm::ParameterSet&);
        ~GERHAnalyzer();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:

        edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
        edm::EDGetTokenT<std::vector<SimTrack>> track_token;
        edm::EDGetTokenT<TrackingParticleCollection>          tpToken_ ;
        edm::EDGetTokenT<GEMRecHitCollection>          rhToken_ ;


        TString outFileName;
        HistGetter hists;
};


GERHAnalyzer::GERHAnalyzer(const edm::ParameterSet& iConfig)
   : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
    shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
    track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
    tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );
    rhToken_    = consumes<GEMRecHitCollection>( edm::InputTag("gemRecHits") );
}

GERHAnalyzer::~GERHAnalyzer() {
  hists.write(outFileName);
}


void
GERHAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{




        edm::Handle <std::vector<SimTrack> > tracks;
      iEvent.getByToken(track_token,tracks);

      edm::Handle<std::vector<PSimHit> >  simHitH ;
      iEvent.getByToken(shToken_,simHitH);

      edm::Handle<TrackingParticleCollection>  TPCollectionH ;
      iEvent.getByToken(tpToken_,TPCollectionH);

      edm::Handle<GEMRecHitCollection> rechitsH;
      iEvent.getByToken(rhToken_,rechitsH);

      edm::ESHandle<GEMGeometry> gemg;
      iSetup.get<MuonGeometryRecord>().get(gemg);
      const GEMGeometry* ggeom = &*gemg;


      hists.getOrMake1D("incl_nEvents",";# of events",1,0,2)->Fill(1.0);

      auto simMuons = ME0Helper::fillSimMuonsByTP(TPCollectionH,*tracks,*simHitH);

      int nInAccept = 0;
      for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
          const auto& muon = simMuons[iM];
          const float absETA = TMath::Abs(muon.trPart->eta());
          if(absETA > 1.0 && absETA < 2.8) nInAccept++;

      }
//      if(nInAccept) return;

      hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);

      for(const auto& d: (*rechitsH) ){
          auto did = d.gemId();
          TString prefix = "";
          if(did.station() == 1)
                prefix = "ge11_";
            else
                prefix = "ge21_";

          const auto* epart = ggeom->etaPartition(did);
          auto glob = epart->toGlobal(d.localPosition());


          hists.getOrMake1D(TString::Format("%s_bx",prefix.Data()),";BX",41,-20.5,20.5)->Fill(d.BunchX());
          hists.getOrMake2D(TString::Format("%s_byLay_bx",prefix.Data()),";bx ;layer",41,-20.5,20.5,8,-0.5,7.5)->Fill(d.BunchX(),did.layer());
          hists.getOrMake2D(TString::Format("%s_byCH_bx",prefix.Data()),";bx ;CH",41,-20.5,20.5,50,-0.5,49.5)->Fill(d.BunchX(),did.chamber());
          hists.getOrMake1D(TString::Format("%s_hitsByRadius",prefix.Data()),";radius",1000,0,500)->Fill(glob.perp());
          hists.getOrMake1D(TString::Format("%s_hitsByZ",prefix.Data()),";radius",2000,-1000,1000)->Fill(glob.z());
          hists.getOrMake1D(TString::Format("%s_hitsByETA",prefix.Data()),";#eta",1000,-5,5)->Fill(glob.eta());
          hists.getOrMake1D(TString::Format("%s_hitsByPhi",prefix.Data()),";#phi",1000,-5,5)->Fill(glob.phi().phi());



      }
   }



//define this as a plug-in
DEFINE_FWK_MODULE(GERHAnalyzer);
