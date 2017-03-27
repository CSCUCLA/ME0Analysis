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

#include <DataFormats/MuonReco/interface/ME0Muon.h>

#include <TVector2.h>
#include <TMath.h>

using namespace std;

class ME0TrackMatchingAnalyzer : public edm::EDAnalyzer {
public:
	explicit ME0TrackMatchingAnalyzer(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
	  runName(iConfig.getUntrackedParameter<std::string>("runName"))
	{
		dToken_        = consumes<ME0DigiPreRecoCollection>( edm::InputTag("simMuonME0Digis") );
		  newDToken_         = consumes<ME0DigiPreRecoCollection>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  oldToNewMapToken_  = consumes<ME0DigiPreRecoMap>( edm::InputTag(iConfig.getParameter<std::string>("newDigiCollection")) );
		  shToken_    = consumes<std::vector<PSimHit>>( edm::InputTag("g4SimHits","MuonME0Hits","SIM") );
		  track_token = consumes<std::vector<SimTrack>>( edm::InputTag("g4SimHits","","SIM") );
		  segToken_    = consumes<ME0SegmentCollection>( edm::InputTag(iConfig.getParameter<std::string>("segmentCollection")) );
		  rhToken_    = consumes<ME0RecHitCollection>( edm::InputTag(iConfig.getParameter<std::string>("recHitCollection")) );


		  strToken_   = consumes<reco::SimToRecoCollection>( edm::InputTag("trackingParticleRecoTrackAsssociation") );
		  tpToken_    = consumes<TrackingParticleCollection>( edm::InputTag("mix","MergedTrackTruth","HLT") );
		  retrack_token_ = consumes<reco::TrackCollection>( edm::InputTag("generalTracks") );

		  muonsToken = consumes<std::vector<reco::ME0Muon>>(iConfig.getParameter<edm::InputTag>("muonsTag"));


	}

	~ME0TrackMatchingAnalyzer() {
		hists.write(outFileName);
	}



private:
	virtual void beginJob() {};
	virtual void endJob() {};
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
template<typename Comp>
void plotTrackComp(const ME0Chamber * chamber, const Comp& comp, ME0Helper::PropogatedTrack& track, TString prefix){
	ME0Helper::LocalPropogatedTrack locTrack =  ME0Helper::getLocalPropogateTrack(track, chamber);

	LocalPoint compLocPos = comp.localPosition();
	LocalVector compLocDir = comp.localDirection();
	AlgebraicSymMatrix compLocCov = comp.parametersError();


	hists.getOrMake1D(TString::Format("%s_track_sigmax",prefix.Data()),   ";track #sigma_{x}",1000,-5,5)->   Fill(std::sqrt(locTrack.cov[3][3]));
	hists.getOrMake1D(TString::Format("%s_track_sigmay",prefix.Data()),   ";track #sigma_{y}",1000,-20,20)-> Fill(std::sqrt(locTrack.cov[4][4]));
	hists.getOrMake1D(TString::Format("%s_track_sigmaDXDZ",prefix.Data()),";track #sigma_{dxdz}",1000,-5,5)->Fill(std::sqrt(locTrack.cov[1][1]));
	hists.getOrMake1D(TString::Format("%s_track_sigmaDYDZ",prefix.Data()),";track #sigma_{dydz}",1000,-5,5)->Fill(std::sqrt(locTrack.cov[2][2]));

	hists.getOrMake1D(TString::Format("%s_comp_sigmax",prefix.Data()),   ";comp #sigma_{x}",1000,-5,5)->      Fill(std::sqrt(compLocCov[2][2]));
	hists.getOrMake1D(TString::Format("%s_comp_sigmay",prefix.Data()),   ";comp #sigma_{y}",1000,-20,20)->    Fill(std::sqrt(compLocCov[3][3]));
	hists.getOrMake1D(TString::Format("%s_comp_sigmaDXDZ",prefix.Data()),";comp #sigma_{dxdz}",1000,-5,5)->   Fill(std::sqrt(compLocCov[0][0]));
	hists.getOrMake1D(TString::Format("%s_comp_sigmaDYDZ",prefix.Data()),";comp #sigma_{dydz}",1000,-5,5)->   Fill(std::sqrt(compLocCov[1][1]));

	hists.getOrMake1D(TString::Format("%s_total_sigmax",prefix.Data()),   ";total #sigma_{x}",1000,-5,5)->   Fill(std::sqrt(locTrack.cov[3][3] + compLocCov[2][2]  )       );
	hists.getOrMake1D(TString::Format("%s_total_sigmay",prefix.Data()),   ";total #sigma_{y}",1000,-20,20)-> Fill(std::sqrt(locTrack.cov[4][4] + compLocCov[3][3]  )       );
	hists.getOrMake1D(TString::Format("%s_total_sigmaDXDZ",prefix.Data()),";total #sigma_{dxdz}",1000,-5,5)->Fill(std::sqrt(locTrack.cov[1][1] + compLocCov[0][0]) );
	hists.getOrMake1D(TString::Format("%s_total_sigmaDYDZ",prefix.Data()),";total #sigma_{dydz}",1000,-5,5)->Fill(std::sqrt(locTrack.cov[2][2] + compLocCov[1][1]));
	auto locTrackP = locTrack.ltp.vector();
	hists.getOrMake1D(TString::Format("%s_diff_x",prefix.Data()),   ";track - comp x",1000,-20,20)->   Fill(locTrackP[3]-compLocPos.x()                  );
	hists.getOrMake1D(TString::Format("%s_diff_y",prefix.Data()),   ";track - comp y",1000,-20,20)-> Fill(locTrackP[4]-compLocPos.y()                   );

	hists.getOrMake2D(TString::Format("%s_diff_xy",prefix.Data()),   ";track - comp x;track - comp y",1000,-20,20,1000,-20,20)->   Fill(locTrackP[3]-compLocPos.x(), locTrackP[4]-compLocPos.y());



	hists.getOrMake1D(TString::Format("%s_diff_DXDZ",prefix.Data()),";track - comp dxdz",1000,-5,5)->Fill(locTrackP[1]-compLocDir.x()/compLocDir.z());
	hists.getOrMake1D(TString::Format("%s_diff_DYDZ",prefix.Data()),";track - comp dydz",1000,-5,5)->Fill(locTrackP[2]-compLocDir.y()/compLocDir.z()) ;

//	double totChi2 = 0;
//	totChi2 += (locTrackP[3]-comp.x   )*(locTrackP[3]-comp.x   )/(locTrack.cov[3][3] + comp.cov_x_x         )       ;
//	totChi2 += (locTrackP[4]-comp.y   )*(locTrackP[4]-comp.y   )/(locTrack.cov[4][4] + comp.cov_y_y         )       ;
//	totChi2 += (locTrackP[1]-comp.dxdz)*(locTrackP[1]-comp.dxdz)/(locTrack.cov[1][1] + comp.cov_dxdz_dxdz) ;
//	totChi2 += (locTrackP[2]-comp.dydz)*(locTrackP[2]-comp.dydz)/(locTrack.cov[2][2] + comp.cov_dydz_dydz);
//	double chi2Prob = TMath::Prob(totChi2,4);
//	hists.getOrMake1D(TString::Format("%s_diff_chi2",prefix.Data()),";track - comp #Chi^{2}",1000,0,100)->Fill(totChi2);
//	hists.getOrMake1D(TString::Format("%s_diff_chi2Prob",prefix.Data()),";track - comp #Chi^{2} prob",1000,0,1)->Fill(chi2Prob);


	hists.getOrMake1D(TString::Format("%s_resid_x",prefix.Data()),   ";(track - comp)/sigma x",500,-10,10)->   Fill((locTrackP[3]-compLocPos.x()               )/std::sqrt(locTrack.cov[3][3] + compLocCov[2][2])       );
	hists.getOrMake1D(TString::Format("%s_resid_y",prefix.Data()),   ";(track - comp)/sigma y",500,-10,10)->   Fill((locTrackP[4]-compLocPos.y()               )/std::sqrt(locTrack.cov[4][4]   + compLocCov[3][3])         );
	hists.getOrMake1D(TString::Format("%s_resid_DXDZ",prefix.Data()),";(track - comp)/sigma dxdz",500,-10,10)->Fill((locTrackP[1]-compLocDir.x()/compLocDir.z())/std::sqrt(locTrack.cov[1][1] + compLocCov[0][0]) );
	hists.getOrMake1D(TString::Format("%s_resid_DYDZ",prefix.Data()),";(track - comp)/sigma dydz",500,-10,10)->Fill((locTrackP[2]-compLocDir.y()/compLocDir.z())/std::sqrt(locTrack.cov[2][2] + compLocCov[1][1]));

	hists.getOrMake1D(TString::Format("%s_trackresid_x",prefix.Data()),   ";(track - comp)/track sigma x",500,-10,10)->   Fill((locTrackP[3]-compLocPos.x()               )/std::sqrt(locTrack.cov[3][3] ));
	hists.getOrMake1D(TString::Format("%s_trackresid_y",prefix.Data()),   ";(track - comp)/track sigma y",500,-10,10)->   Fill((locTrackP[4]-compLocPos.y()               )/std::sqrt(locTrack.cov[4][4]   )  );
	hists.getOrMake1D(TString::Format("%s_trackresid_DXDZ",prefix.Data()),";(track - comp)/track sigma dxdz",500,-10,10)->Fill((locTrackP[1]-compLocDir.x()/compLocDir.z())/std::sqrt(locTrack.cov[1][1] ));
	hists.getOrMake1D(TString::Format("%s_trackresid_DYDZ",prefix.Data()),";(track - comp)/track sigma dydz",500,-10,10)->Fill((locTrackP[2]-compLocDir.y()/compLocDir.z())/std::sqrt(locTrack.cov[2][2] ));


	GlobalPoint segGlobPos(chamber->toGlobal(compLocPos));
    GlobalVector segGlobalVector =  chamber->toGlobal(compLocDir);
    float compDPhi = TVector2::Phi_mpi_pi(segGlobalVector.phi() - segGlobPos.phi());
    float trackDPhi = TVector2::Phi_mpi_pi(track.p3FinalReco_glob.phi() - track.r3FinalReco_globv.phi());


	hists.getOrMake1D(TString::Format("%s_diff_eta",prefix.Data()),   ";track - comp #eta",1000,-1,1)->   Fill(track.r3FinalReco_globv.eta() -segGlobPos.eta());
	hists.getOrMake1D(TString::Format("%s_diff_phi",prefix.Data()),   ";track - comp #phi",1000,-.5,.5)-> Fill(TVector2::Phi_mpi_pi(track.r3FinalReco_globv.phi() -segGlobPos.phi()));
	hists.getOrMake1D(TString::Format("%s_diff_DEta",prefix.Data()),   ";track - comp #Delta#eta",1000,-1,1)->   Fill(track.p3FinalReco_glob.eta() -segGlobalVector.eta());
	hists.getOrMake1D(TString::Format("%s_diff_DPhi",prefix.Data()),   ";track - comp #Delta#phi",1000,-2,2)-> Fill(TVector2::Phi_mpi_pi(trackDPhi - compDPhi));
	hists.getOrMake1D(TString::Format("%s_comp_DPhi",prefix.Data()),   ";comp #Delta#phi",1000,-2,2)-> Fill(TVector2::Phi_mpi_pi(compDPhi));
	hists.getOrMake1D(TString::Format("%s_track_DPhi",prefix.Data()),   ";track #Delta#phi",1000,-2,2)-> Fill(TVector2::Phi_mpi_pi(trackDPhi));

	auto extrap = [&] (const LocalPoint& point, const LocalVector& dir, double extZ) -> LocalPoint {
	    double extX = point.x()+extZ*dir.x()/dir.z();
	    double extY = point.y()+extZ*dir.y()/dir.z();
	    return LocalPoint(extX,extY,extZ);
	  };
	auto getDPhi2 = [&](const LocalPoint& point, const LocalVector& dir) -> float {
		float chamberZ = chamber->position().z();
		LocalPoint projHigh = extrap(point,dir, ME0Helper::endOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
		LocalPoint projLow = extrap(point,dir, ME0Helper::beginOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
	    auto globLow  = chamber->toGlobal(projLow );
		auto globHigh = chamber->toGlobal(projHigh);
		return TVector2::Phi_mpi_pi(globHigh.phi() - globLow.phi());
	};
	float compDPhi2 = getDPhi2(compLocPos,compLocDir);
	float trackDPhi2 = getDPhi2(LocalPoint(locTrackP[3],locTrackP[4],0),LocalVector(locTrackP[1],locTrackP[2],1));
	hists.getOrMake1D(TString::Format("%s_diff_DPhi2",prefix.Data()),   ";track - comp #Delta#phi 2",1000,-.5,.5)-> Fill(TVector2::Phi_mpi_pi(trackDPhi2 - compDPhi2));
	hists.getOrMake1D(TString::Format("%s_comp_DPhi2",prefix.Data()),   ";comp #Delta#phi 2",1000,-.5,.5)-> Fill(TVector2::Phi_mpi_pi(compDPhi2));
	hists.getOrMake1D(TString::Format("%s_track_DPhi2",prefix.Data()),   ";track #Delta#phi 2",1000,-.5,.5)-> Fill(TVector2::Phi_mpi_pi(trackDPhi2));
	hists.getOrMake2D(TString::Format("%s_diff_xy_dph2i",prefix.Data()),   ";track - comp max[x,y];track - comp #Delta#phi 2",500,0,-20,1000,-.5,.5)->   Fill( std::max(std::fabs(locTrackP[3]-compLocPos.x()),std::fabs(locTrackP[4]-compLocPos.y())),TVector2::Phi_mpi_pi(trackDPhi2 - compDPhi2));


//	float chamberZ = chamber->position().z();
//	LocalPoint projHigh = extrap(compLocPos,compLocDir, ME0Helper::endOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
//	LocalPoint projLow = extrap(compLocPos,compLocDir, ME0Helper::beginOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
//	LocalVector cdphi4( projHigh.x() - projLow.x(), projHigh.y() - projLow.y(), projHigh.z() - projLow.z()   );
//    auto globLow  = chamber->toGlobal(projLow );
//	auto globHigh = chamber->toGlobal(projHigh);
//	GlobalVector segGlobalVector2 =  chamber->toGlobal(cdphi4);

//	cout << endl << "C "<< compLocPos <<" "<< compLocDir<<" " << segGlobPos << " " << segGlobalVector<< "("<<compDPhi<<") :"<< locTrack.ltp.position() <<" "<< locTrack.ltp.direction()<<
////			" "<< track.r3FinalReco_globv <<" "<< track.p3FinalReco_glob <<"("<< trackDPhi <<") " <<" :: "<< chamber->toGlobal(locTrack.ltp.position()) <<" "<<chamber->toGlobal(locTrack.ltp.direction())  << compDPhi2 <<" "<<trackDPhi2<<std::endl;

//	cout << prefix <<" "<< compLocPos << " "<<compLocDir <<" "<< projHigh << " "<< projLow<< " "<<cdphi4<<" "<<segGlobalVector2<<" "<< segGlobalVector2.phi()<<" "<< globLow <<" "<< globHigh<<" "<< TVector2::Phi_mpi_pi(globHigh.phi() - globLow.phi()) <<endl;

//	cout << "TRK" <<" "<< LocalPoint(locTrackP[3],locTrackP[4],0) << " "<<LocalVector(locTrackP[1],locTrackP[2],1)  <<" "<<track.p3FinalReco_glob<<" "<< track.p3FinalReco_glob.phi()<<" "<< trackDPhi2 <<endl;

}

void getResPlots(const ME0Geometry* mgeom, const ME0SegmentCollection& segments, const ME0RecHitCollection& recHits, const std::vector<SimTrack>&  simTrackH, const std::vector<PSimHit>&  simHitH,
		const ME0DigiPreRecoCollection& oldDigis,const ME0DigiPreRecoCollection& newDigis,const ME0DigiPreRecoMap& digiMap,
		const edm::Handle<TrackingParticleCollection>& trParticles, const reco::SimToRecoCollection& simToReco , const edm::ESHandle<MagneticField>& bField,const edm::ESHandle<Propagator>& ThisshProp,
		const std::vector<reco::ME0Muon> & me0Muons, TString sname){

	auto simMuons = ME0Helper::fillSimMuons(simTrackH,simHitH);
	ME0Helper::DigiInfoMap digiInfo;
	ME0Helper::fillDigiInfoMap(newDigis,digiMap,oldDigis,simHitH,digiInfo);
	ME0Helper::associateSimMuons(simMuons,digiInfo);
	ME0Helper::SegmentCatMap segmentCats;
	ME0Helper::fillSegmentCategories(segments,recHits,digiInfo,simMuons,segmentCats);
	ME0Helper::associateSimMuonsToTracks(simMuons, trParticles, simToReco);
//	cout <<endl<< simMuons.size()<<" ";
	for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
		const auto& muon = simMuons[iM];
		const float pt = muon.track->momentum().pt();
		const float absETA = TMath::Abs(muon.track->momentum().eta());
//		cout <<"(" <<pt <<" "<< absETA <<") ";
		if(pt < 1) continue;
		if(absETA > 2.8 || absETA < 2.0 ) continue;
		hists.getOrMake1D(TString::Format("%sall_muon_pt",sname.Data()),";muon p_{T}",60,0,30)->Fill(pt);

		if(muon.recoTrack.isNull()) continue;

		ME0Helper::PropogatedTrack propTr = ME0Helper::propogateTrack(bField,ThisshProp,&*muon.recoTrack, 539.35);
		if(!propTr.isValid) continue;

		hists.getOrMake1D(TString::Format("%sgoodTrack_muon_pt",sname.Data()),";muon p_{T}",60,0,30)->Fill(pt);

		bool positiveEndcap = muon.recoTrack->pz() > 0;
		bool positiveCharge = muon.track->charge() > 0;

		auto simHitProp = ME0Helper::getSimTrackProperties(mgeom,muon.simHits);

		if(simHitProp.nPrimLaysHit >= 4){
			hists.getOrMake1D(TString::Format("%sgoodTrack_goodSimTrack_muon_pt",sname.Data()),";muon p_{T}",60,0,30)->Fill(pt);

			if(pt >= 1){
				std::string ptstr = "pteq1to3_";
				if(pt >= 3 && pt < 5 ) ptstr  = "pteq3to5_";
				else if(pt >= 5 && pt < 20 ) ptstr = "pteq5to20_";
				else  if(pt >= 20) ptstr = "ptgeq20_";
				plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_"+ptstr);
			}

			if(positiveCharge && positiveEndcap) plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_posQposE_");
			if(positiveCharge && !positiveEndcap) plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_posQnegE_");
			if(!positiveCharge && positiveEndcap) plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_negQposE_");
			if(!positiveCharge && !positiveEndcap) plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_negQnegE_");



				plotTrackComp(simHitProp.chamber,simHitProp,propTr,sname+"sim_");
		}


		int typeC = 3; //Perfect/MissOne/MissMultiple/Lost
		if(muon.segments.size()) {
			int nRD = (segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second)->nRecHits();
			int nSD = muon.segments[0].second;
			if(segmentCats[muon.segments[0].first.first][muon.segments[0].first.second].first == ME0Helper::MUON_COMP_PURE) typeC = 0;
			else if(nRD - nSD <= 1 ) typeC = 1;
			else if( nSD > (nRD - nSD)) typeC = 2;
		}
		if(typeC <= 2){
			const ME0Chamber * chamb = mgeom->chamber(muon.segments[0].first.first);
			const ME0Segment * segment = &*(segments.get(muon.segments[0].first.first).first + muon.segments[0].first.second);
			hists.getOrMake1D(TString::Format("%sgoodTrack_goodSegment_muon_pt",sname.Data()),";muon p_{T}",60,0,30)->Fill(pt);
//			cout << propTr.r3FinalReco_globv.z() <<" "<< chamb->position().z() <<" "<< segment->localPosition() <<" "<< segment->localDirection()<< " ";
//			ME0Helper::LocalPropogatedTrack locTrack =  ME0Helper::getLocalPropogateTrack(propTr, chamb);
//			   GlobalVector segGlobalVector =  chamb->toGlobal(segment->localDirection());

//			std::cout << locTrack.ltp.position() <<" "<< locTrack.ltp.direction()<<
//					" -> "<<  segGlobalVector<<" "<<propTr.p3FinalReco_glob <<" :: "<< segGlobalVector.phi() << " "<<propTr.p3FinalReco_glob.phi() <<std::endl;


			if(pt >= 1){
				std::string ptstr = "pteq1to3_";
				if(pt >= 3 && pt < 5 ) ptstr  = "pteq3to5_";
				else if(pt >= 5 && pt < 20 ) ptstr = "pteq5to20_";
				else  if(pt >= 20) ptstr = "ptgeq20_";
				plotTrackComp(chamb,*segment,propTr,sname+"seg_"+ptstr);
			}
			if(positiveCharge && positiveEndcap)   plotTrackComp(chamb,*segment,propTr,sname+"seg_posQposE_");
			if(positiveCharge && !positiveEndcap)  plotTrackComp(chamb,*segment,propTr,sname+"seg_posQnegE_");
			if(!positiveCharge && positiveEndcap)  plotTrackComp(chamb,*segment,propTr,sname+"seg_negQposE_");
			if(!positiveCharge && !positiveEndcap) plotTrackComp(chamb,*segment,propTr,sname+"seg_negQnegE_");
			plotTrackComp(chamb,*segment,propTr,sname+"seg_");
//			if(segment->localPosition().x() < -7 )
//				plotTrackComp(chamb,*segment,propTr,sname+"seg_ltm7_");
//			else if(segment->localPosition().x() > 7 )
//				plotTrackComp(chamb,*segment,propTr,sname+"seg_gt7_");
//			else
//				plotTrackComp(chamb,*segment,propTr,sname+"seg_m7to7_");
		}

	}



	//FAKES!
	int nFakes = 0;
	int nFakesGeq3 = 0;
	int nFakesGeq3PassDPhi = 0;
	for(const auto& me0Muon : me0Muons){
		const auto * chamber = mgeom->chamber(me0Muon.me0segment().me0DetId());
		auto ch_segs = segments.get(me0Muon.me0segment().me0DetId());
		auto& ch_cats = segmentCats[me0Muon.me0segment().me0DetId()];
		int idx = me0Muon.me0segid() - int(ch_segs.first - segments.begin());

		if(!(me0Muon.me0segment().localPosition() == (ch_segs.first + idx)->localPosition())){
			cout << endl<<"ERRORRR "<< me0Muon.me0segment() <<" <-> "<< *(ch_segs.first + idx)<<endl;
		}

		auto cat =ch_cats[idx].first;
		if(cat <= ME0Helper::MUON_MISS_DIRTY_NEUT) continue;
		ME0Helper::PropogatedTrack propTr = ME0Helper::propogateTrack(bField,ThisshProp,&*me0Muon.innerTrack(), 539.35);
		if(!propTr.isValid) continue;
		auto prop = ME0Helper::getSegmentProperties(chamber,&me0Muon.me0segment());
		const bool passDPhi = std::fabs(prop.dPhi) <0.013;
		const bool passTime = std::fabs(me0Muon.me0segment().time()) <11.0;

		float pt = me0Muon.innerTrack()->pt();

		hists.getOrMake1D(TString::Format("%sfakeSeg_pt",sname.Data()),";fake p_{T}",60,0,30)->Fill(pt);
		if(passDPhi) hists.getOrMake1D(TString::Format("%sfakeSeg_passDPhi_pt",sname.Data()),";fake p_{T}",60,0,30)->Fill(pt);
		if(passTime) hists.getOrMake1D(TString::Format("%sfakeSeg_passTime_pt",sname.Data()),";fake p_{T}",60,0,30)->Fill(pt);
		if(passDPhi && passTime) hists.getOrMake1D(TString::Format("%sfakeSeg_passDPhiTime_pt",sname.Data()),";fake p_{T}",60,0,30)->Fill(pt);

		if(!passTime) continue;

		nFakes++;

		if(pt >= 3) nFakesGeq3++;
		if(pt >= 3 && passDPhi) nFakesGeq3PassDPhi++;
		if(pt >= 1){
			std::string ptstr = "pteq1to3_";
			if(pt >= 3 && pt < 5 ) ptstr  = "pteq3to5_";
			else if(pt >= 5 && pt < 20 ) ptstr = "pteq5to20_";
			else  if(pt >= 20) ptstr = "ptgeq20_";
			plotTrackComp(chamber,me0Muon.me0segment(),propTr,sname+"fakeseg_"+ptstr);
			if(passDPhi) plotTrackComp(chamber,me0Muon.me0segment(),propTr,sname+"fakeseg_passDPhi_"+ptstr);
		}

		plotTrackComp(chamber,me0Muon.me0segment(),propTr,sname+"fakeseg_");
		plotTrackComp(chamber,me0Muon.me0segment(),propTr,sname+"fakeseg_passDPhi_");
	}
	hists.getOrMake1D(TString::Format("%sfake_nMuons",sname.Data()),";# of muons",100,-0.5,99.5)->Fill(nFakes);
	hists.getOrMake1D(TString::Format("%sfake_ptgeq3_nMuons",sname.Data()),";# of muons",100,-0.5,99.5)->Fill(nFakesGeq3);
	hists.getOrMake1D(TString::Format("%sfake_ptgeq3_passDPhi_nMuons",sname.Data()),";# of muons",100,-0.5,99.5)->Fill(nFakesGeq3PassDPhi);

}

private:
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          dToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoCollection>          newDToken_  ;
	edm::EDGetTokenT<ME0DigiPreRecoMap>                 oldToNewMapToken_  ;
    edm::EDGetTokenT<std::vector<PSimHit>>          shToken_ ;
    edm::EDGetTokenT<ME0SegmentCollection>          segToken_ ;
    edm::EDGetTokenT<std::vector<SimTrack>> track_token;
    edm::EDGetTokenT<ME0RecHitCollection>          rhToken_ ;

    edm::EDGetTokenT<reco::TrackCollection>           retrack_token_ ;
    edm::EDGetTokenT<reco::SimToRecoCollection>          strToken_ ;
    edm::EDGetTokenT<TrackingParticleCollection>          tpToken_ ;

    edm::EDGetTokenT<std::vector<reco::ME0Muon>>          muonsToken ;


	TString outFileName;
	TString runName;
	HistGetter hists;
};




void
ME0TrackMatchingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

//	edm::Handle<ME0DigiPreRecoCollection> digisH;
//	iEvent.getByToken(dToken_,digisH);
//
//	edm::Handle<ME0DigiPreRecoCollection> digisNH;
//	iEvent.getByToken(newDToken_,digisNH);
//
//	edm::Handle<ME0DigiPreRecoMap> digisM;
//	iEvent.getByToken(oldToNewMapToken_,digisM);
//
//		edm::Handle <std::vector<SimTrack> > tracks;
//	  iEvent.getByToken(track_token,tracks);
//
//	  edm::Handle<std::vector<PSimHit> >  simHitH ;
//	  iEvent.getByToken(shToken_,simHitH);

	  edm::Handle<ME0SegmentCollection >  segH ;
	  iEvent.getByToken(segToken_,segH);

//	  edm::Handle<ME0RecHitCollection> rechitsH;
//	  iEvent.getByToken(rhToken_,rechitsH);

	edm::ESHandle<ME0Geometry> me0g;
	iSetup.get<MuonGeometryRecord>().get(me0g);
	const ME0Geometry* mgeom = &*me0g;
	hists.getOrMake1D("nEvents",";# of events",1,0,2)->Fill(1.0);


//	  edm::ESHandle<MagneticField> bField;
//	  iSetup.get<IdealMagneticFieldRecord>().get(bField);
//
//	  edm::ESHandle<Propagator> ThisshProp;
//	  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAlong", ThisshProp);

	  edm::Handle <reco::TrackCollection > retracks;
	  iEvent.getByToken(retrack_token_,retracks);

//	  edm::Handle<reco::SimToRecoCollection> simRecH;
//	  iEvent.getByToken(strToken_,simRecH);
//
//	  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
//	  iEvent.getByToken(tpToken_,TPCollectionH);
//
//
//	  edm::Handle<std::vector<reco::ME0Muon>> muonCollectionH;
//	  iEvent.getByToken(muonsToken, muonCollectionH);

	  int nTracksPT1 = 0;
	  int nTracksPT2 = 0;

	  for(const auto& t: *retracks){
	      if (std::abs(t.eta()) < 1.8) continue;
	      if (std::abs(t.eta()) > 3.0) continue;
	      if(t.pt() < 1) continue;
	      nTracksPT1++;
	      if(t.pt() >= 2) nTracksPT2++;
	  }
	   int nSeg = 0;
	   int nSegDPhi = 0;

	for(const auto& s: *segH){
		if(std::fabs(s.time()) >= 11.0 ) continue;
		nSeg++;
    	auto extrap = [&] (const LocalPoint& point, const LocalVector& dir, double extZ) -> LocalPoint {
    	    double extX = point.x()+extZ*dir.x()/dir.z();
    	    double extY = point.y()+extZ*dir.y()/dir.z();
    	    return LocalPoint(extX,extY,extZ);
    	  };
	   	auto getDPhiDEta= [&](const LocalPoint& point, const LocalVector& dir, float& dPhi, float& dEta)  {
	   		const ME0Chamber * chamber = mgeom->chamber(s.me0DetId());
	   		float chamberZ = chamber->position().z();
	   		LocalPoint projHigh = extrap(point,dir, ME0Helper::endOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
	   		LocalPoint projLow = extrap(point,dir, ME0Helper::beginOfDet-   (chamberZ < 0 ? -1.0 : 1.0) * chamberZ);
	   	    auto globLow  = chamber->toGlobal(projLow );
	   		auto globHigh = chamber->toGlobal(projHigh);
	   		dPhi =  TVector2::Phi_mpi_pi(globHigh.phi() - globLow.phi());
	   		dEta =  globHigh.eta() - globLow.eta();
	   	};
    	float segDPhi, segDEta;
    	getDPhiDEta(s.localPosition(),s.localDirection(),segDPhi,segDEta);
		if(std::fabs(segDPhi) < 0.013) nSegDPhi++;
	}

	hists.getOrMake1D(TString::Format("%s_nTrPT1"  ,TString::Format("%s_",runName.Data()).Data()),";nTracks",1000,-0.5,999.5)->Fill(nTracksPT1);
	hists.getOrMake1D(TString::Format("%s_nTrPT2"  ,TString::Format("%s_",runName.Data()).Data()),";nTracks",1000,-0.5,999.5)->Fill(nTracksPT2);
	hists.getOrMake1D(TString::Format("%s_nSegDP2" ,TString::Format("%s_",runName.Data()).Data()),";nSegs",1000,-0.5,999.5)->Fill(nSeg);
	hists.getOrMake1D(TString::Format("%s_nSegDP13",TString::Format("%s_",runName.Data()).Data()),";nSegs",1000,-0.5,999.5)->Fill(nSegDPhi);


//	  getResPlots(mgeom,*segH,*rechitsH,*tracks,*simHitH,*digisH,*digisNH,*digisM,TPCollectionH, *simRecH,bField,ThisshProp, *muonCollectionH,TString::Format("%s_",runName.Data()));


}



//define this as a plug-in
DEFINE_FWK_MODULE(ME0TrackMatchingAnalyzer);
