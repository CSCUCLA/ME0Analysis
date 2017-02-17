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


#include "AnalysisSupport/Utilities/interface/HistGetter.h"



using namespace std;

class ME0TrackDensity : public edm::EDAnalyzer {
    public:
        explicit ME0TrackDensity(const edm::ParameterSet&);
        ~ME0TrackDensity();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<reco::TrackCollection> track_token;
        edm::EDGetTokenT<reco::VertexCollection>          vtxToken_;
        edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;

        TString outFileName;
        HistGetter hists;
};


ME0TrackDensity::ME0TrackDensity(const edm::ParameterSet& iConfig) :
    outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
  track_token = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );
  vtxToken_   = consumes<reco::VertexCollection>( edm::InputTag("offlinePrimaryVertices") );
//  dToken_     = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
}


ME0TrackDensity::~ME0TrackDensity() {hists.write(outFileName);}


void
ME0TrackDensity::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle <reco::TrackCollection > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle <reco::VertexCollection > verts;
  iEvent.getByToken(vtxToken_,verts);

//  edm::Handle<ME0DigiPreRecoCollection> digisH;
//  iEvent.getByToken(dToken_,digisH);


//  edm::ESHandle<ME0Geometry> me0g;
//  iSetup.get<MuonGeometryRecord>().get(me0g);
//  const ME0Geometry* mgeom = &*me0g;

  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);

  int nV = verts->size();
  int nGV = 0;
  for (std::vector<reco::Vertex>::const_iterator iV = verts->begin();
      iV != verts->end(); ++iV){
    if(!iV->isFake() &&
        iV->ndof() > 4. &&
        iV->position().Rho() <= 2.0 &&
        fabs(iV->position().Z())<=24.0)
      nGV++;
  }


  hists.getOrMake1D("nVerts",";# of vert.",1000,-.5,999.5)->Fill(nV);
  hists.getOrMake1D("nGoodVerts",";# of good vert.",1000,-.5,999.5)->Fill(nGV);

  for (std::vector<reco::Track>::const_iterator iTrack = tracks->begin();
      iTrack != tracks->end(); ++iTrack){

    if(iTrack->pt() > 0.5){
      hists.getOrMake1D("tracks_gt0p5_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt0p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }
    if(iTrack->pt() > 1.5){
      hists.getOrMake1D("tracks_gt1p5_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt1p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }
    if(iTrack->pt() > 3.0){
      hists.getOrMake1D("tracks_gt3p0_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt3p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }

    if(iTrack->pt() > 5.0){
      hists.getOrMake1D("tracks_gt5p0_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt5p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }

    if(iTrack->pt() > 10.0){
      hists.getOrMake1D("tracks_gt10_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt10_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }
    if(iTrack->pt() > 20.0){
      hists.getOrMake1D("tracks_gt20_eta",";#eta",800,-4.0,4.0)->Fill(iTrack->eta());
      hists.getOrMake1D("tracks_gt20_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iTrack->eta()));
    }

  }

//  for(auto iD = digisH->begin(); iD != digisH->end(); ++iD){
//    auto detID = (*iD).first;
//    if(detID.layer() != 1) continue;
//    const ME0DigiPreRecoCollection::Range& range = (*iD).second;
//    for (ME0DigiPreRecoCollection::const_iterator idigi = range.first;
//        idigi != range.second;idigi++) {
//      float eta = mgeom->etaPartition(detID)->toGlobal(LocalPoint(idigi->x(),idigi->y(),0)).eta();
//      if(TMath::Abs(idigi->pdgid()) == 13){
//        hists.getOrMake1D("me0hits_muons_eta",";#eta",800,-4.0,4.0)->Fill(eta);
//        hists.getOrMake1D("me0hits_muons_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(eta));
//      }
//      else{
//        hists.getOrMake1D("me0hits_nonmuons_eta",";#eta",800,-4.0,4.0)->Fill(eta);
//        hists.getOrMake1D("me0hits_nonmuons_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(eta));
//      }
//
//    }
//
//  }


}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0TrackDensity);
