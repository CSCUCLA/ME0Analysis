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
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"


#include "AnalysisSupport/Utilities/interface/HistGetter.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class GEDigiAnalyzer : public edm::EDAnalyzer {
    public:
        explicit GEDigiAnalyzer(const edm::ParameterSet&);
        ~GEDigiAnalyzer();


    private:
        virtual void beginJob() {};
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() {};


    private:
        edm::EDGetTokenT<GEMDigiCollection>          dToken_  ;



        TString outFileName;
        HistGetter hists;
};


GEDigiAnalyzer::GEDigiAnalyzer(const edm::ParameterSet& iConfig)
   : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName"))
{
	  dToken_     = consumes<GEMDigiCollection>( edm::InputTag("simMuonGEMDigis") );
}

GEDigiAnalyzer::~GEDigiAnalyzer() {
  hists.write(outFileName);
}


void
GEDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	  edm::Handle<GEMDigiCollection> digisH;
	  iEvent.getByToken(dToken_,digisH);

	  edm::ESHandle<GEMGeometry> gemg;
	  iSetup.get<MuonGeometryRecord>().get(gemg);
	  const GEMGeometry* ggeom = &*gemg;

  hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);


  for(const auto& d: (*digisH) ){
		auto did = d.first;
		TString prefix = "";
		if(did.station() == 1)
			prefix = "ge11_";
		else
			prefix = "ge21_";

   const auto* epart = ggeom->etaPartition(did);

   for (GEMDigiCollection::const_iterator idigi = d.second.first;
       idigi != d.second.second;idigi++) {


	      hists.getOrMake1D(TString::Format("%s_tof",prefix.Data()),";tof",41,-20.5,20.5)->Fill(idigi->bx());
	      hists.getOrMake2D(TString::Format("%s_byLay_tof",prefix.Data()),";tof ;layer",41,-20.5,20.5,8,-0.5,7.5)->Fill(idigi->bx(),d.first.layer());
	      hists.getOrMake1D(TString::Format("%s_hits",prefix.Data()),";layer",8,-0.5,7.5)->Fill(d.first.layer());


	      if(d.first.layer() == 1)
	    	  hists.getOrMake1D(TString::Format("%s_hitsByRadius",prefix.Data()),";radius",1000,0,500)->Fill(epart->toGlobal(epart->centreOfStrip(idigi->strip())).perp());







   }

   }





}



//define this as a plug-in
DEFINE_FWK_MODULE(GEDigiAnalyzer);
