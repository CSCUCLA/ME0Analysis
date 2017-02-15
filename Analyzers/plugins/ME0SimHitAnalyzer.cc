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

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <TVector2.h>
#include <TMath.h>

#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"


using namespace std;


struct SimTrackProperties{
	double cenEta = -1;
	double lay1Eta = -1;
	double dPhi = -1;
	int nLaysHit = 0;
	bool inFirst= false;
};
SimTrackProperties getSimTrackProperties( const ME0Geometry* mgeom, const std::vector<std::pair<const PSimHit*,ME0DetId> >& simHits) {

	SimTrackProperties prop;

    std::vector<std::vector<const PSimHit*> > laysHit(6);
    for(const auto& h : simHits){
      laysHit[h.second.layer()  - 1].push_back(h.first);
    }
    std::vector<const PSimHit *> prunedHits;
    for(auto& l : laysHit){
      if(l.size() >0) {
    	  prop.nLaysHit++;
        prunedHits.push_back(l.front());
      }
    }

    auto getGlobal = [&](const PSimHit * sh) -> GlobalPoint {
    	return mgeom->etaPartition(sh->detUnitId())->toGlobal(sh->entryPoint());
    };

    if(laysHit[0].size()){
    	prop.inFirst = true;
    	prop.lay1Eta = getGlobal(laysHit[0].front()).eta();
    }


    if(prop.nLaysHit > 2){
        const double beginOfDet  = 527;
        const double centerOfDet = 539.5;
        const double endOfDet    = 552;
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
        prop.dPhi = TVector2::Phi_mpi_pi(upPt.phi() - downPt.phi());

    }
    return prop;
}


class ME0SimHitAnalyzer : public edm::EDAnalyzer {
    public:
        explicit ME0SimHitAnalyzer(const edm::ParameterSet&);
        ~ME0SimHitAnalyzer();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<std::vector<SimTrack>> track_token;
        edm::EDGetTokenT<std::vector<SimVertex>> vert_token;
        edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
        edm::EDGetTokenT<reco::GenParticleCollection> gen_token;


        TString outFileName;
        HistGetter hists;
};


ME0SimHitAnalyzer::ME0SimHitAnalyzer(const edm::ParameterSet& iConfig)
   : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
  vert_token = consumes<std::vector<SimVertex>>( edm::InputTag("g4SimHits","","SIM") );
  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
  gen_token = consumes<reco::GenParticleCollection>( edm::InputTag("genParticles") );


}


ME0SimHitAnalyzer::~ME0SimHitAnalyzer() {
  hists.write(outFileName);
}



typedef std::map<unsigned int, std::vector<std::pair<const PSimHit*,ME0DetId> > > SimHCont;



void
ME0SimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle <std::vector<SimTrack> > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle <std::vector<SimVertex> > verts;
  iEvent.getByToken(vert_token,verts);

  edm::Handle<std::vector<PSimHit> >  simHitH ;
  iEvent.getByToken(shToken_,simHitH);

  edm::Handle <reco::GenParticleCollection > genParts;
  iEvent.getByToken(gen_token,genParts);

  edm::ESHandle<ME0Geometry> me0g;
  iSetup.get<MuonGeometryRecord>().get(me0g);
  const ME0Geometry* mgeom = &*me0g;

  SimHCont hitcont;

  for (const auto & hit: *simHitH)
  {
    hitcont[hit.trackId()].emplace_back(&hit,ME0DetId(hit.detUnitId()));

  }



  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);


    for (SimHCont::iterator it=hitcont.begin(); it!=hitcont.end(); ++it){
      int vertT = -1;
      int idx = -1;
      for(unsigned int iT = 0; iT < tracks->size(); ++iT)
        if(tracks->at(iT).trackId() == it->first){
          idx = iT; break;
        }
      if(idx >=  0){
        const auto& track = tracks->at(idx);
         if(!track.noVertex()){
           const auto& vert = verts->at(track.vertIndex());
           vertT = vert.processType();
         }

      }

      int partID = TMath::Abs(it->second.front().first->particleType());
      hists.getOrMake2D("partTypes",";Vertex type;PDGID",303,-1.5,301.5,41,-.5,40.5)->Fill(vertT,partID < 38 ? partID : 40);

      //Filltype
      int type = -1;
      TString prefix;
      if(partID == 13){
        prefix = "muon";
        if(vertT == 0){
          type = 1;
        }
        else if(vertT >= 201)
          type = 2;
        else
          type = 3;
      } else if (partID == 11){
        prefix = "ele";
        if(vertT >= 1 && vertT <=23)
          type = 4;
        else
          type = 5;
      } else if (partID >= 38){
        prefix = "hadron";
        if(vertT == 0)
          type = 6;
        else if(vertT >= 201)
          type = 7;
        else if (vertT == 111 || vertT == 121)
          type = 8;
        else
          type = 9;
      } else {
        prefix = "other";
        type = 10;
      }


      hists.getOrMake1D("particleBreakdown",";Muons[123],Electron[45],Hadron[6-9],Other",10   ,.5,10.5)->Fill(type);


      //Numbers special plots;
      SimTrackProperties props = getSimTrackProperties(mgeom, it->second);
      if(props.inFirst)                                                                           hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(1);
      if(props.nLaysHit  && TMath::Abs(props.lay1Eta) > 2.4)                                      hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(2);
      if(props.nLaysHit >= 3)                                                                     hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(3);
      if(props.nLaysHit >= 3 && TMath::Abs(props.cenEta) > 2.4)                                   hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(4);
      if(props.nLaysHit >= 3 && TMath::Abs(props.dPhi) < 0.013)                                   hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(5);
      if(props.nLaysHit >= 3 && TMath::Abs(props.dPhi) < 0.013 && TMath::Abs(props.cenEta) > 2.4) hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI",10   ,.5,10.5)->Fill(6);
      if(props.nLaysHit >= 6 && TMath::Abs(props.dPhi) < 0.013)                                   hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI,6lDPHI",10   ,.5,10.5)->Fill(7);
      if(props.nLaysHit >= 6 && TMath::Abs(props.dPhi) < 0.013 && TMath::Abs(props.cenEta) > 2.4) hists.getOrMake1D(TString::Format("%s_specialYields",prefix.Data()),";1l,3l,3lDPHI,6lDPHI",10   ,.5,10.5)->Fill(8);



      double minETA = 0;
      double minPhi = 0;
      double maxPhi = 0;
      int minLay = -1;
      int maxLay = -1;

      std::vector<int> laysHit(6,0);
      for(const auto& h : it->second){
        laysHit[h.second.layer()  - 1] ++;
        auto gl = mgeom->etaPartition(h.second )->toGlobal(h.first->entryPoint());
//        std::vector<float> origPars = mgeom->etaPartition(h.second )->specs()->parameters();
//
//        cout << origPars[0] <<" "<< origPars[1] <<" "<< origPars[2] <<endl;


        if(minLay < 0 || h.second.layer() < minLay){
          minLay = h.second.layer();
          minPhi = gl.phi();
          minETA = gl.eta();
        } else if(maxLay < 0 || h.second.layer() > maxLay){
          maxLay = h.second.layer();
          maxPhi = gl.phi();
        }
      }
      int nLaysHit = 0;
      for(auto& l : laysHit){if(l >0) nLaysHit++;}
      hists.getOrMake1D(TString::Format("%s_nLaysHit",prefix.Data()),";# of hit layers",7,-0.5,6.5)->Fill(nLaysHit);

      for(unsigned int iL = 0; iL < laysHit.size(); ++iL){
        if(laysHit[iL]) hists.getOrMake1D(TString::Format("%s_hitLays",prefix.Data()),"; layer",7,-0.5,6.5)->Fill(iL);
        if(laysHit[iL]) hists.getOrMake2D(TString::Format("%s_xNhits_hitLays",prefix.Data()),"; layer ; nLays",7,-0.5,6.5,7,-0.5,6.5)->Fill(iL, nLaysHit);
      }

      if(idx >=  0){
        const auto& track = tracks->at(idx);
        hists.getOrMake1D(TString::Format("%s_trackPT",prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
        hists.getOrMake1D(TString::Format("%s_trackP",prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());

        if(nLaysHit == 1){
          hists.getOrMake1D(TString::Format("%s_1l_trackTheta",prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
          hists.getOrMake1D(TString::Format("%s_1l_trackPT",prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
          hists.getOrMake1D(TString::Format("%s_1l_trackP",prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());

          hists.getOrMake1D(TString::Format("%s_1l_pdgID",prefix.Data()),";track p",6000,-0.5,5999.5)->Fill(TMath::Abs(track.type()));
          if(!track.noVertex()){
          hists.getOrMake1D(TString::Format("%s_1l_posZ",prefix.Data()),";track p",600,0,600)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().z()));
          hists.getOrMake1D(TString::Format("%s_1l_posETA",prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));

          hists.getOrMake1D(TString::Format("%s_1l_posETA",prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));
          hists.getOrMake1D(TString::Format("%s_1l_shETA",prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(minETA));

          hists.getOrMake1D(TString::Format("%s_1l_posPHI",prefix.Data()),";track p",600,-6.5,6.5)->Fill(verts->at(track.vertIndex()).position().phi());
          }

        } else if (nLaysHit == 2){
          hists.getOrMake1D(TString::Format("%s_2l_trackPT",prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
          hists.getOrMake1D(TString::Format("%s_2l_trackP",prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());
          hists.getOrMake1D(TString::Format("%s_2l_trackTheta",prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
        } else {
          hists.getOrMake1D(TString::Format("%s_geq3l_trackPT",prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
          hists.getOrMake1D(TString::Format("%s_geq3l_trackP",prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());
          hists.getOrMake1D(TString::Format("%s_geq3l_trackTheta",prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
          hists.getOrMake1D(TString::Format("%s_geq3l_pdgID",prefix.Data()),";track p",6000,-0.5,5999.5)->Fill(TMath::Abs(track.type()));
          if(!track.noVertex()){
          hists.getOrMake1D(TString::Format("%s_geq3l_posZ",prefix.Data()),";track p",600,0,600)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().z()));
          hists.getOrMake1D(TString::Format("%s_geq3l_posETA",prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));
          hists.getOrMake1D(TString::Format("%s_geq3l_posPHI",prefix.Data()),";track p",600,-6.5,6.5)->Fill(verts->at(track.vertIndex()).position().phi());
          hists.getOrMake1D(TString::Format("%s_geq3l_shETA",prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(minETA));
          }
        }

      }

      if(nLaysHit >= 3){
        double dPhi = TMath::Abs(TVector2::Phi_mpi_pi(maxPhi - minPhi));
        hists.getOrMake1D(TString::Format("%s_nLays_geq3_dPhi",prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
        if(nLaysHit >= 4)hists.getOrMake1D(TString::Format("%s_nLays_geq4_dPhi",prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
        if(nLaysHit >= 5)hists.getOrMake1D(TString::Format("%s_nLays_geq5_dPhi",prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
        if(nLaysHit >= 6)hists.getOrMake1D(TString::Format("%s_nLays_geq6_dPhi",prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);

        if(type == 1) hists.getOrMake1D(TString::Format("%s_pi_nLays_geq3_dPhi",prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
      }




    }



//  cout << endl << "-------------------------New event-------------------------"<<endl;
//
//  for (SimHCont::iterator it=hitcont.begin(); it!=hitcont.end(); ++it){
//    int idx = -1;
//    for(unsigned int iT = 0; iT < tracks->size(); ++iT)
//      if(tracks->at(iT).trackId() == it->first){
//        idx = iT; break;
//      }
//    if(idx >=  0){
//      const auto& track = tracks->at(idx);
//       std::cout << "Track: " << track.type() << " ("<< track.momentum().pt() <<","<< track.momentum().eta()<<","<< track.momentum().phi() <<")" <<endl;
//       if(!track.noGenpart()){
//         const auto& gen =  genParts->at(track.genpartIndex());
//         std::cout << "GEN: " << gen.pdgId() <<" "<< gen.status() << "("<< gen.pt() <<","<< gen.eta()<<","<< gen.phi() <<")" <<endl;
//       }
//       if(!track.noVertex()){
//         const auto& vert = verts->at(track.vertIndex());
//         cout <<"Vert: "<< vert.vertexId() <<" "<< vert.processType() <<" ("<< vert.position().eta()<<","<< vert.position().phi() <<")" << endl;
//
//       }
//
//    } else {
//      cout <<"None: "<< it->first <<" ("<< tracks->size() <<")" <<endl;
//    }
//    for(const auto& h : it->second){
//      auto gl = mgeom->etaPartition(h.second )->toGlobal(h.first->entryPoint());
//      std::cout << h.second.chamber() <<" : " << h.second.layer() <<" "<< h.first->particleType() <<" ("<<gl.eta()<<","<<gl.phi()<<") " <<endl;
//    }
//  }

}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0SimHitAnalyzer);
