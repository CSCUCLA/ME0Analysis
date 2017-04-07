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
#include "DataFormats/GEMRecHit/interface/ME0Segment.h"

#include "AnalysisSupport/Utilities/interface/HistGetter.h"
#include "../interface/ME0Helper.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class ME0BackgroundAnalyzer : public edm::EDAnalyzer {
public:
	explicit ME0BackgroundAnalyzer(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
	  runName(iConfig.getUntrackedParameter<std::string>("runName"))
	{
		  newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  segToken_    = consumes<ME0SegmentCollection>( edm::InputTag(iConfig.getParameter<std::string>("segmentCollection")) );

	}

	~ME0BackgroundAnalyzer() {
		hists.write(outFileName);
	}



private:
	virtual void beginJob() {};
	virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);



private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
    edm::EDGetTokenT<ME0SegmentCollection>          segToken_ ;



	TString outFileName;
	TString runName;
	HistGetter hists;



	void makeDigiPlots(const ME0Geometry* mgeom,  const ME0DigiPreRecoCollection& digis, TString prefix){


		std::vector<std::vector<std::vector<uint32_t>>> chMap(2, std::vector<std::vector<uint32_t>>(18,std::vector<uint32_t>(6,0) ) );


		for(const auto& ch_digis: digis ){
			const ME0DetId detID = ch_digis.first;
			const auto* epart = mgeom->etaPartition(detID);
			const unsigned int lay = detID.layer();
			chMap[detID.region() < 0 ? 0 : 1][detID.chamber() - 1][lay-1] =detID.rawId();

			for (ME0DigiPreRecoCollection::const_iterator idigi = ch_digis.second.first;
					idigi != ch_digis.second.second;idigi++) {
				if(idigi->tof() > 10 || idigi->tof() < -10) continue;
				bool isNeutron = !idigi->prompt();
				auto fillPlots = [&](   TString varName, TString xTitle, int nX, float minX, float maxX, float var) {
					if(isNeutron)hists.getOrMake1D(TString::Format("%s_neutrons_%s",prefix.Data(), varName.Data()),xTitle,nX,minX,maxX)->Fill(var);
					else hists.getOrMake1D(TString::Format("%s_pu_%s",prefix.Data(), varName.Data()),xTitle,nX,minX,maxX)->Fill(var);
					hists.getOrMake1D(TString::Format("%s_all_%s",prefix.Data(), varName.Data()),xTitle,nX,minX,maxX)->Fill(var);

					if(isNeutron)hists.getOrMake1D(TString::Format("%s_neutrons_lay%u_%s",prefix.Data(), lay, varName.Data()),xTitle,nX,minX,maxX)->Fill(var);
					else hists.getOrMake1D(TString::Format("%s_pu_lay%u_%s",prefix.Data(), lay,varName.Data()),xTitle,nX,minX,maxX)->Fill(var);
					hists.getOrMake1D(TString::Format("%s_all_lay%u_%s",prefix.Data(), lay, varName.Data()),xTitle,nX,minX,maxX)->Fill(var);

				};

				GlobalPoint glb = epart->toGlobal(LocalPoint(idigi->x(),idigi->y()));
				fillPlots("hitsByR",";radius [cm]",1000,0,200,glb.perp());
				float absE = std::fabs(glb.eta());
				fillPlots("hitsByEta",";#eta",1000,1.8,3.0,absE);
			}
		}

		for(unsigned int iE = 0; iE < 2; ++iE){
			for(unsigned int iC = 0; iC < 18; ++iC){
				for(unsigned int iL = 0; iL < 6; ++iL){
					int nNeutron = 0;
					int nPU = 0;

					int nNeutron1VFAT[] = {0,0,0};
					int nPU1VFAT[] = {0,0,0};

					uint32_t rawID = chMap[iE][iC][iL];
					if(rawID != 0){
						const auto& ch_digis = digis.get(ME0DetId(rawID));
						for (ME0DigiPreRecoCollection::const_iterator idigi = ch_digis.first;
								idigi != ch_digis.second;idigi++) {
							if(idigi->x() < -3.95 && idigi->y() < -39){
								if(std::fabs(idigi->tof()) < 27){
									int idx = 0;
									if (idigi->tof() < -12.5) idx = 0;
									else if(idigi->tof() < 12.5) idx = 1;
									else idx = 2;

									if(!idigi->prompt()) nNeutron1VFAT[idx]++;
												else nPU1VFAT[idx]++;
								}
							}
							if(std::fabs(idigi->tof()) > 10) continue;
							if(!idigi->prompt()) nNeutron++;
							else nPU++;
						}
					}

					hists.getOrMake1D(TString::Format("%s_lay%u_numberOfNeutronDigis",prefix.Data(), iL+1),";# of neutron digis",200,-0.5,199.5)->Fill(nNeutron);
					hists.getOrMake1D(TString::Format("%s_lay%u_numberOfPUDigis",prefix.Data(), iL+1),";# of pu digis",200,-0.5,199.5)->Fill(nPU);
					hists.getOrMake1D(TString::Format("%s_lay%u_numberOfDigis",prefix.Data(), iL+1),";# of digis",200,-0.5,199.5)->Fill(nNeutron+nPU);

					hists.getOrMake1D(TString::Format("%s_numberOfNeutronDigis",prefix.Data()),";# of neutron digis",200,-0.5,199.5)->Fill(nNeutron);
					hists.getOrMake1D(TString::Format("%s_numberOfPUDigis",prefix.Data()),";# of pu digis",200,-0.5,199.5)->Fill(nPU);
					hists.getOrMake1D(TString::Format("%s_numberOfDigis",prefix.Data()),";# of digis",200,-0.5,199.5)->Fill(nNeutron+nPU);


					int nVF10 =  ((nPU1VFAT[0] + nNeutron1VFAT[0]) > 0);
					int nVF11 =  ((nPU1VFAT[1] + nNeutron1VFAT[1]) > 0);
					int nVF12 =  ((nPU1VFAT[2] + nNeutron1VFAT[2]) > 0);
					int countVF1 = nVF10+nVF11+nVF12;

					hists.getOrMake1D(TString::Format("%s_numberOfVFAT1On",prefix.Data()),";VFAT ON",4,-0.5,3.5)->Fill(countVF1);
					hists.getOrMake1D(TString::Format("%s_numberOfNeutronDigisVFAT10",prefix.Data()),";# of neutron digis",200,-0.5,199.5)->Fill(nNeutron1VFAT[0]);
					hists.getOrMake1D(TString::Format("%s_numberOfPUDigisVFAT10",prefix.Data()),";# of pu digis",200,-0.5,199.5)->Fill(nPU1VFAT[0]);
					hists.getOrMake1D(TString::Format("%s_numberOfDigisVFAT10",prefix.Data()),";# of digis",200,-0.5,199.5)->Fill(nPU1VFAT[0] + nNeutron1VFAT[0] );
					hists.getOrMake1D(TString::Format("%s_numberOfNeutronDigisVFAT11",prefix.Data()),";# of neutron digis",200,-0.5,199.5)->Fill(nNeutron1VFAT[1]);
					hists.getOrMake1D(TString::Format("%s_numberOfPUDigisVFAT11",prefix.Data()),";# of pu digis",200,-0.5,199.5)->Fill(nPU1VFAT[1]);
					hists.getOrMake1D(TString::Format("%s_numberOfDigisVFAT11",prefix.Data()),";# of digis",200,-0.5,199.5)->Fill(nPU1VFAT[1] + nNeutron1VFAT[1] );
					hists.getOrMake1D(TString::Format("%s_numberOfNeutronDigisVFAT12",prefix.Data()),";# of neutron digis",200,-0.5,199.5)->Fill(nNeutron1VFAT[2]);
					hists.getOrMake1D(TString::Format("%s_numberOfPUDigisVFAT12",prefix.Data()),";# of pu digis",200,-0.5,199.5)->Fill(nPU1VFAT[2]);
					hists.getOrMake1D(TString::Format("%s_numberOfDigisVFAT12",prefix.Data()),";# of digis",200,-0.5,199.5)->Fill(nPU1VFAT[2] + nNeutron1VFAT[2] );

				}
			}

		}

	}


	void makeSegmentPlots(const ME0Geometry* mgeom,  const ME0SegmentCollection& segments, TString prefix){

		std::vector<std::vector<uint32_t>> chSegCount(2, std::vector<uint32_t>(18,0 ) );


		for(const auto& seg : segments ){
//			const auto* ch = mgeom->chamber(seg.me0DetId().chamberId());
//			GlobalPoint seggp = ch->toGlobal(seg.localPosition());
//			LocalPoint seglp = seg.localPosition();
			if(std::fabs(seg.time()) >= 11.0 ) continue;
			chSegCount[seg.me0DetId().chamberId().region() < 0 ? 0 : 1][seg.me0DetId().chamber() -1]++;
		}


		for(unsigned int iE = 0; iE < 2; ++iE){
			int sectN = 0;
			for(unsigned int iC = 0; iC < 18; ++iC){
				hists.getOrMake1D(TString::Format("%s_numberOfChamberSegments",prefix.Data()),";# of segments",40,-0.5,49.5)->Fill(chSegCount[iE][iC]);
				if(sectN == 2){
					hists.getOrMake1D(TString::Format("%s_numberOfSectorSegments",prefix.Data()),";# of segments",40,-0.5,49.5)->Fill(chSegCount[iE][iC]+chSegCount[iE][iC-1]+chSegCount[iE][iC-2]);
					sectN = 0;
				} else sectN++;
			}
	}

	}


};






void
ME0BackgroundAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


	edm::Handle<ME0DigiPreRecoCollection> digisNH;
	iEvent.getByToken(newDToken_,digisNH);

	  edm::Handle<ME0SegmentCollection >  segH ;
	  iEvent.getByToken(segToken_,segH);


	edm::ESHandle<ME0Geometry> me0g;
	iSetup.get<MuonGeometryRecord>().get(me0g);
	const ME0Geometry* mgeom = &*me0g;

	hists.getOrMake1D(TString::Format("%snEvents",runName.Data()),";# of events",1,0,2)->Fill(1.0);

	makeDigiPlots(mgeom,*digisNH,TString::Format("%s",runName.Data()));
	makeSegmentPlots(mgeom,*segH,TString::Format("%s",runName.Data()));
}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0BackgroundAnalyzer);
