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
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"


#include "../../AnalysisSupport/interface/HistGetter.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <TVector2.h>
#include <TMath.h>

#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

using namespace std;

class GESimHitAnalyzer : public edm::EDAnalyzer {
    public:
	typedef std::map<unsigned int, std::vector<std::pair<const PSimHit*,GEMDetId> > > SimHCont;
        explicit GESimHitAnalyzer(const edm::ParameterSet&);
        ~GESimHitAnalyzer();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        void fillPlots(const GEMGeometry* ggeo, SimHCont& hits, edm::Handle <std::vector<SimTrack> >& tracks, edm::Handle <std::vector<SimVertex> >& verts, TString name);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<std::vector<SimTrack>> track_token;
        edm::EDGetTokenT<std::vector<SimVertex>> vert_token;
        edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
        edm::EDGetTokenT<reco::GenParticleCollection> gen_token;


        TString outFileName;
        HistGetter hists;
};


GESimHitAnalyzer::GESimHitAnalyzer(const edm::ParameterSet& iConfig)
   : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
  vert_token = consumes<std::vector<SimVertex>>( edm::InputTag("g4SimHits","","SIM") );
  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonGEMHits","SIM") );
  gen_token = consumes<reco::GenParticleCollection>( edm::InputTag("genParticles") );


}


GESimHitAnalyzer::~GESimHitAnalyzer() {
  hists.write(outFileName);
}




void GESimHitAnalyzer::fillPlots(const GEMGeometry* ggeo, SimHCont& hits, edm::Handle <std::vector<SimTrack> >& tracks, edm::Handle <std::vector<SimVertex> >& verts, TString name){
	for (SimHCont::iterator it=hits.begin(); it!=hits.end(); ++it){
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
		hists.getOrMake2D(name + "_partTypes",";Vertex type;PDGID",303,-1.5,301.5,41,-.5,40.5)->Fill(vertT,partID < 38 ? partID : 40);

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

		hists.getOrMake1D(name +"_particleBreakdown",";Muons[123],Electron[45],Hadron[6-9],Other",10   ,.5,10.5)->Fill(type);

		double minETA = 0;
		double minPhi = 0;
		double maxPhi = 0;
		int minLay = -1;
		int maxLay = -1;

		std::vector<int> laysHit(6,0);
		for(const auto& h : it->second){
			laysHit[h.second.layer()  - 1] ++;
				auto gl = ggeo->etaPartition(h.second  )->toGlobal(h.first->entryPoint());
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
		hists.getOrMake1D(TString::Format("%s_%s_nLaysHit",name.Data(),prefix.Data()),";# of hit layers",7,-0.5,6.5)->Fill(nLaysHit);

		for(unsigned int iL = 0; iL < laysHit.size(); ++iL){
			if(laysHit[iL]) hists.getOrMake1D(TString::Format("%s_%s_hitLays",name.Data(),prefix.Data()),"; layer",7,-0.5,6.5)->Fill(iL);
			if(laysHit[iL]) hists.getOrMake2D(TString::Format("%s_%s_xNhits_hitLays",name.Data(),prefix.Data()),"; layer ; nLays",7,-0.5,6.5,7,-0.5,6.5)->Fill(iL, nLaysHit);
		}

		if(idx >=  0){
			const auto& track = tracks->at(idx);
			hists.getOrMake1D(TString::Format("%s_%s_trackPT",name.Data(),prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
			hists.getOrMake1D(TString::Format("%s_%s_trackP",name.Data(),prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());

			if(nLaysHit == 1){
				hists.getOrMake1D(TString::Format("%s_%s_1l_trackTheta",name.Data(),prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
				hists.getOrMake1D(TString::Format("%s_%s_1l_trackPT",name.Data(),prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
				hists.getOrMake1D(TString::Format("%s_%s_1l_trackP",name.Data(),prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());

				hists.getOrMake1D(TString::Format("%s_%s_1l_pdgID",name.Data(),prefix.Data()),";track p",6000,-0.5,5999.5)->Fill(TMath::Abs(track.type()));
				if(!track.noVertex()){
					hists.getOrMake1D(TString::Format("%s_%s_1l_posZ",name.Data(),prefix.Data()),";track p",600,0,600)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().z()));
					hists.getOrMake1D(TString::Format("%s_%s_1l_posETA",name.Data(),prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));

					hists.getOrMake1D(TString::Format("%s_%s_1l_posETA",name.Data(),prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));
					hists.getOrMake1D(TString::Format("%s_%s_1l_shETA",name.Data(),prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(minETA));

					hists.getOrMake1D(TString::Format("%s_%s_1l_posPHI",name.Data(),prefix.Data()),";track p",600,-6.5,6.5)->Fill(verts->at(track.vertIndex()).position().phi());
				}

			} else if (nLaysHit == 2){
				hists.getOrMake1D(TString::Format("%s_%s_2l_trackPT",name.Data(),prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
				hists.getOrMake1D(TString::Format("%s_%s_2l_trackP",name.Data(),prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());
				hists.getOrMake1D(TString::Format("%s_%s_2l_trackTheta",name.Data(),prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
			} else {
				hists.getOrMake1D(TString::Format("%s_%s_geq3l_trackPT",name.Data(),prefix.Data()),";track p_{T}",100,0,10)->Fill(track.momentum().pt());
				hists.getOrMake1D(TString::Format("%s_%s_geq3l_trackP",name.Data(),prefix.Data()),";track p",100,0,10)->Fill(track.momentum().P());
				hists.getOrMake1D(TString::Format("%s_%s_geq3l_trackTheta",name.Data(),prefix.Data()),";track momentum eta",100,-6,6)->Fill(TMath::ASin(track.momentum().pt()/track.momentum().P()));
				hists.getOrMake1D(TString::Format("%s_%s_geq3l_pdgID",name.Data(),prefix.Data()),";track p",6000,-0.5,5999.5)->Fill(TMath::Abs(track.type()));
				if(!track.noVertex()){
					hists.getOrMake1D(TString::Format("%s_%s_geq3l_posZ",name.Data(),prefix.Data()),";track p",600,0,600)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().z()));
					hists.getOrMake1D(TString::Format("%s_%s_geq3l_posETA",name.Data(),prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(verts->at(track.vertIndex()).position().eta()));
					hists.getOrMake1D(TString::Format("%s_%s_geq3l_posPHI",name.Data(),prefix.Data()),";track p",600,-6.5,6.5)->Fill(verts->at(track.vertIndex()).position().phi());
					hists.getOrMake1D(TString::Format("%s_%s_geq3l_shETA",name.Data(),prefix.Data()),";track p",600,0,6)->Fill(TMath::Abs(minETA));
				}
			}

		}

		if(nLaysHit >= 3){
			double dPhi = TMath::Abs(TVector2::Phi_mpi_pi(maxPhi - minPhi));
			hists.getOrMake1D(TString::Format("%s_%s_nLays_geq3_dPhi",name.Data(),prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
			if(nLaysHit >= 4)hists.getOrMake1D(TString::Format("%s_%s_nLays_geq4_dPhi",name.Data(),prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
			if(nLaysHit >= 5)hists.getOrMake1D(TString::Format("%s_%s_nLays_geq5_dPhi",name.Data(),prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
			if(nLaysHit >= 6)hists.getOrMake1D(TString::Format("%s_%s_nLays_geq6_dPhi",name.Data(),prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);

			if(type == 1) hists.getOrMake1D(TString::Format("%s_%s_pi_nLays_geq3_dPhi",name.Data(),prefix.Data()),";#Delta#phi",20,0,0.05)->Fill(dPhi);
		}




	}

}

void
GESimHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle <std::vector<SimTrack> > tracks;
  iEvent.getByToken(track_token,tracks);

  edm::Handle <std::vector<SimVertex> > verts;
  iEvent.getByToken(vert_token,verts);

  edm::Handle<std::vector<PSimHit> >  simHitH ;
  iEvent.getByToken(shToken_,simHitH);

  edm::Handle <reco::GenParticleCollection > genParts;
  iEvent.getByToken(gen_token,genParts);

//  edm::ESHandle<ME0Geometry> me0g;
//  iSetup.get<MuonGeometryRecord>().get(me0g);
//  const ME0Geometry* mgeom = &*me0g;

  edm::ESHandle<GEMGeometry> gemg;
  iSetup.get<MuonGeometryRecord>().get(gemg);
  const GEMGeometry* ggeom = &*gemg;

  SimHCont ge11Hitcont;
  SimHCont ge21Hitcont;





  for (const auto & hit: *simHitH)
  {

	auto did = GEMDetId(hit.detUnitId());
//	auto gl = ggeom->etaPartition(did )->toGlobal(h.first->entryPoint());
	          std::vector<float> origPars = ggeom->etaPartition(did )->specs()->parameters();
	          auto p = ggeom->etaPartition(did )->position();


	cout << did.station() <<" "<<did.chamber()<<" "<< did.layer()<<" "<<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<< origPars[0] <<" "<< origPars[1] <<" "<< origPars[2] <<endl;
    if(did.station() == 1)
    	ge11Hitcont[hit.trackId()].emplace_back(&hit,GEMDetId(hit.detUnitId()));
    else
    	ge21Hitcont[hit.trackId()].emplace_back(&hit,GEMDetId(hit.detUnitId()));

  }



  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);
  fillPlots(ggeom,ge11Hitcont,tracks,verts,"ge11_");
  fillPlots(ggeom,ge21Hitcont,tracks,verts,"ge21_");
}



//define this as a plug-in
DEFINE_FWK_MODULE(GESimHitAnalyzer);
