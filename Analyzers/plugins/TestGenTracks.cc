// user include files

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "../../AnalysisSupport/interface/HistGetter.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"


using namespace std;

class TreeContainer {
public:
  TreeContainer(TString fileName, TString treeName, TString treeTitle){
    file = new TFile(fileName, "RECREATE");
    tree = new TTree(treeName,treeTitle);
  }
  void write() {
    file->cd();
    tree->Write();
    file->Close();
    delete file;
  }

  template<class T>
  void    book(const char *name, T& var, const char *type) { tree->Branch(name, &var, TString(name).Append("/").Append(type).Data()); }


  void fill() {tree->Fill();}
  TFile * file;
  TTree * tree;

};

class TestGenTracks : public edm::EDAnalyzer {
    public:
        explicit TestGenTracks(const edm::ParameterSet&);
        ~TestGenTracks();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<reco::GenParticleCollection> gen_token;
        TString outFileName;
//        TreeContainer tree;
        HistGetter hists;


        int n_trks_pt_gt0p5      = 0;
        int n_trks_pt_gt3p0      = 0;
        int n_trks_pt_gt5p0      = 0;
        int n_trks_pt_gt10       = 0;
        int n_trks_pt_gt20       = 0;
        int n_muon_trks_pt_gt0p5 = 0;
        int n_muon_trks_pt_gt3p0 = 0;
        int n_muon_trks_pt_gt5p0 = 0;
        int n_muon_trks_pt_gt10  = 0;
        int n_muon_trks_pt_gt20  = 0;

};


TestGenTracks::TestGenTracks(const edm::ParameterSet& iConfig) :
    outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))//,tree(outFileName,"Events","")
{
  gen_token = consumes<reco::GenParticleCollection>( edm::InputTag("genParticles") );

//  tree.book("n_trks_pt_gt0p5",n_trks_pt_gt0p5      ,"I");
//  tree.book("n_trks_pt_gt3p0",n_trks_pt_gt3p0      ,"I");
//  tree.book("n_trks_pt_gt5p0",n_trks_pt_gt5p0      ,"I");
//  tree.book("n_trks_pt_gt10",n_trks_pt_gt10       ,"I");
//  tree.book("n_trks_pt_gt20",n_trks_pt_gt20       ,"I");
//  tree.book("n_muon_trks_pt_gt0p5",n_muon_trks_pt_gt0p5 ,"I");
//  tree.book("n_muon_trks_pt_gt3p0",n_muon_trks_pt_gt3p0 ,"I");
//  tree.book("n_muon_trks_pt_gt5p0",n_muon_trks_pt_gt5p0 ,"I");
//  tree.book("n_muon_trks_pt_gt10" ,n_muon_trks_pt_gt10  ,"I");
//  tree.book("n_muon_trks_pt_gt20" ,n_muon_trks_pt_gt20  ,"I");

}


TestGenTracks::~TestGenTracks() {
  hists.write(outFileName);
  //  tree.write();
}


void
TestGenTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle <reco::GenParticleCollection > genParts;
  iEvent.getByToken(gen_token,genParts);

  n_trks_pt_gt0p5      = 0;
  n_trks_pt_gt3p0      = 0;
  n_trks_pt_gt5p0      = 0;
  n_trks_pt_gt10       = 0;
  n_trks_pt_gt20       = 0;
  n_muon_trks_pt_gt0p5 = 0;
  n_muon_trks_pt_gt3p0 = 0;
  n_muon_trks_pt_gt5p0 = 0;
  n_muon_trks_pt_gt10  = 0;
  n_muon_trks_pt_gt20  = 0;

  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);

  for (std::vector<reco::GenParticle>::const_iterator iG = genParts->begin();
      iG != genParts->end(); ++iG){
    if(iG->status() != 1) continue;
    if(iG->charge() == 0) continue;

    if(iG->pt() > 0.5){
      hists.getOrMake1D("tracks_gt0p5_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt0p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }
    if(iG->pt() > 1.5){
      hists.getOrMake1D("tracks_gt1p5_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt1p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }
    if(iG->pt() > 3.0){
      hists.getOrMake1D("tracks_gt3p0_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt3p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }

    if(iG->pt() > 5.0){
      hists.getOrMake1D("tracks_gt5p0_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt5p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }

    if(iG->pt() > 10.0){
      hists.getOrMake1D("tracks_gt10_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt10_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }
    if(iG->pt() > 20.0){
      hists.getOrMake1D("tracks_gt20_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("tracks_gt20_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
    }

    if(TMath::Abs(iG->pdgId()) == 13){

      hists.getOrMake1D("muon_tracks_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
      hists.getOrMake1D("muon_tracks_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));

      if(iG->pt() > 0.5){
        hists.getOrMake1D("muon_tracks_gt0p5_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt0p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }
      if(iG->pt() > 1.5){
        hists.getOrMake1D("muon_tracks_gt1p5_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt1p5_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }
      if(iG->pt() > 3.0){
        hists.getOrMake1D("muon_tracks_gt3p0_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt3p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }

      if(iG->pt() > 5.0){
        hists.getOrMake1D("muon_tracks_gt5p0_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt5p0_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }

      if(iG->pt() > 10.0){
        hists.getOrMake1D("muon_tracks_gt10_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt10_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }
      if(iG->pt() > 20.0){
        hists.getOrMake1D("muon_tracks_gt20_eta",";#eta",800,-4.0,4.0)->Fill(iG->eta());
        hists.getOrMake1D("muon_tracks_gt20_abseta",";|#eta|",400,0,4.0)->Fill(TMath::Abs(iG->eta()));
      }

    }




////    if(TMath::Abs(iG->eta()) < 2 || TMath::Abs(iG->eta()) > 3 ) continue;
//    if(TMath::Abs(iG->eta()) > 1 ) continue;
//
//    if(iG->pt() > 0.5){
//           n_trks_pt_gt0p5++;
//           if(TMath::Abs(iG->pdgId()) == 13) n_muon_trks_pt_gt0p5++;
//    }
//    if(iG->pt() > 1.5){
//      n_trks_pt_gt0p5++;
//      if(TMath::Abs(iG->pdgId()) == 13)  n_muon_trks_pt_gt0p5++;
//    }
//    if(iG->pt() > 3.0){
//      n_trks_pt_gt3p0++;
//if(TMath::Abs(iG->pdgId()) == 13) n_muon_trks_pt_gt3p0++;
//    }
//
//    if(iG->pt() > 5.0){
//      n_trks_pt_gt5p0++;
//      if(TMath::Abs(iG->pdgId()) == 13) n_muon_trks_pt_gt5p0++;
//    }
//
//    if(iG->pt() > 10.0){
//      n_trks_pt_gt10++;
//      if(TMath::Abs(iG->pdgId()) == 13) n_muon_trks_pt_gt10++;
//    }
//    if(iG->pt() > 20.0){
//      n_trks_pt_gt20++;
//      if(TMath::Abs(iG->pdgId()) == 13) n_muon_trks_pt_gt20++;
//    }
  }
//  tree.fill();
}



//define this as a plug-in
DEFINE_FWK_MODULE(TestGenTracks);
