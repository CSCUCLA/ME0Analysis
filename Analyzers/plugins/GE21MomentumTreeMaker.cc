// user include files

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
//#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
//#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "AnalysisSupport/Utilities/interface/Types.h"
#include "AnalysisSupport/TreeInterface/interface/TreeWriter.h"
#include "AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "AnalysisSupport/Utilities/interface/ParticleInfo.h"
#include "AnalysisSupport/CMSSWUtilities/interface/ParticleUtilities.h"

using namespace std;
class GE21MomentumTreeMaker : public edm::EDAnalyzer {
public:
    explicit GE21MomentumTreeMaker(const edm::ParameterSet& iConfig) : outFileName(iConfig.getUntrackedParameter<std::string>("outFileName")),
    tree (outFileName,"Events","")
    {

        genparticle_token         = consumes<std::vector<reco::GenParticle>>( edm::InputTag("genParticles") );
        recomuon_token            = consumes<std::vector<reco::Muon>>       ( edm::InputTag("muons") );
        picky_def_token           = consumes<reco::TrackToTrackMap>         ( edm::InputTag("tevMuons","picky") );
        picky_noGE21_token        = consumes<reco::TrackToTrackMap>         ( edm::InputTag("tevMuonsNoGE21","picky") );
        picky_noME21_token        = consumes<reco::TrackToTrackMap>         ( edm::InputTag("tevMuonsNoME21","picky") );
        picky_noGE11ME0_token     = consumes<reco::TrackToTrackMap>         ( edm::InputTag("tevMuonsNoGE11ME0","picky") );
        picky_noGE11ME0ME21_token = consumes<reco::TrackToTrackMap>         ( edm::InputTag("tevMuonsNoGE11ME0ME21","picky") );
        picky_oneGE21_token = consumes<reco::TrackToTrackMap>               ( edm::InputTag("tevMuonsOneGE21","picky") );






        nM_all                        = tree.add<ASTypes::size>   ("nM_all"                 ,"i",0);
        nM_good                       = tree.add<ASTypes::size>   ("nM_good"                ,"i",0);
        genMuon_pt                    = tree.add<float>  ("genMuon_pt"             ,"F",-1);
        genMuon_eta                   = tree.add<float>  ("genMuon_eta"            ,"F",-1);
        genMuon_q                     = tree.add<ASTypes::int8>   ("genMuon_q"              ,"B",0);
        genReco_dr                    = tree.add<float>  ("genReco_dr"             ,"F",-1);
        global_pt                     = tree.add<float>  ("global_pt"              ,"F",-1);
        global_eta                    = tree.add<float>  ("global_eta"             ,"F",-1);
        global_q                      = tree.add<ASTypes::int8>   ("global_q"               ,"B",0);

        global_v_d                      = tree.addMulti<ASTypes::int8>   ("global_v_d"               ,0);
        global_v_s                      = tree.addMulti<ASTypes::int8>   ("global_v_s"               ,0);
        global_v_r                      = tree.addMulti<ASTypes::int8>   ("global_v_r"               ,0);
        global_v_phie                   = tree.addMulti<float>           ("global_v_phie"            ,0);

        tracker_pt                    = tree.add<float>  ("tracker_pt"             ,"F",-1);
        tracker_eta                   = tree.add<float>  ("tracker_eta"            ,"F",-1);
        tracker_q                     = tree.add<ASTypes::int8>   ("tracker_q"              ,"B",0);
        standalone_pt                 = tree.add<float>  ("standalone_pt"          ,"F",-1);
        standalone_eta                = tree.add<float>  ("standalone_eta"         ,"F",-1);
        standalone_q                  = tree.add<ASTypes::int8>   ("standalone_q"           ,"B",0);
        picky_pt                      = tree.add<float>  ("picky_pt"               ,"F",-1);
        picky_eta                     = tree.add<float>  ("picky_eta"              ,"F",-1);
        picky_q                       = tree.add<ASTypes::int8>   ("picky_q"                ,"B",0);

        picky_v_d                      = tree.addMulti<ASTypes::int8>   ("picky_v_d"         ,0);
        picky_v_s                      = tree.addMulti<ASTypes::int8>   ("picky_v_s"         ,0);
        picky_v_r                      = tree.addMulti<ASTypes::int8>   ("picky_v_r"         ,0);
        picky_v_phie                   = tree.addMulti<float>  ("picky_v_phie"               ,0);

        picky_oneGE21_pt              = tree.add<float>  ("picky_oneGE21_pt"               ,"F",-1);
        picky_oneGE21_eta             = tree.add<float>  ("picky_oneGE21_eta"              ,"F",-1);
        picky_oneGE21_q               = tree.add<ASTypes::int8>   ("picky_oneGE21_q"                ,"B",0);

        picky_def_pt                  = tree.add<float>  ("picky_def_pt"           ,"F",-1);
        picky_def_eta                 = tree.add<float>  ("picky_def_eta"          ,"F",-1);
        picky_def_q                   = tree.add<ASTypes::int8>   ("picky_def_q"            ,"B",0);
        picky_noGE21_pt               = tree.add<float>  ("picky_noGE21_pt"        ,"F",-1);
        picky_noGE21_eta              = tree.add<float>  ("picky_noGE21_eta"       ,"F",-1);
        picky_noGE21_q                = tree.add<ASTypes::int8>   ("picky_noGE21_q"         ,"B",0);
        picky_noME21_pt               = tree.add<float>  ("picky_noME21_pt"        ,"F",-1);
        picky_noME21_eta              = tree.add<float>  ("picky_noME21_eta"       ,"F",-1);
        picky_noME21_q                = tree.add<ASTypes::int8>   ("picky_noME21_q"         ,"B",0);
        picky_noGE11ME0_pt            = tree.add<float>  ("picky_noGE11ME0_pt"     ,"F",-1);
        picky_noGE11ME0_eta           = tree.add<float>  ("picky_noGE11ME0_eta"    ,"F",-1);
        picky_noGE11ME0_q             = tree.add<ASTypes::int8>   ("picky_noGE11ME0_q"      ,"B",0);
        picky_noGE11ME0ME21_pt        = tree.add<float>  ("picky_noGE11ME0ME21_pt" ,"F",-1);
        picky_noGE11ME0ME21_eta       = tree.add<float>  ("picky_noGE11ME0ME21_eta","F",-1);
        picky_noGE11ME0ME21_q         = tree.add<ASTypes::int8>   ("picky_noGE11ME0ME21_q"  ,"B",0);
        tree.book();
    }

    ~GE21MomentumTreeMaker() {
        tree.write();

    }



private:
    virtual void beginJob() {};
    virtual void endJob() {};
    virtual void analyze(const edm::Event&, const edm::EventSetup&);


    bool isGoodMuon(const reco::Muon& recoMu){
        if(!recoMu.isGlobalMuon()) return false;
        if(recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() == 0) return false;
        if(recoMu.numberOfMatchedStations() <= 1) return false;
        if(recoMu.innerTrack()->hitPattern().numberOfValidPixelHits() == 0) return false;
        if(recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) return false;
        return true;
    }





    void fillTree(const reco::GenParticle& genMuon, const reco::Muon* recoMuon, const reco::TrackToTrackMap& picky_def, const reco::TrackToTrackMap& picky_noGE21, const reco::TrackToTrackMap& picky_noME21,
            const reco::TrackToTrackMap& picky_noGE11ME0, const reco::TrackToTrackMap& picky_noGE11ME0ME21, const reco::TrackToTrackMap& picky_oneGE21){
        tree.fill(genMuon_pt             ,float(genMuon.pt()));
        tree.fill(genMuon_eta            ,float(genMuon.eta()));
        tree.fill(genMuon_q              ,ASTypes::int8(genMuon.charge()));

        if(recoMuon ==0){
            return;
        }


        tree.fill(genReco_dr             ,float(PhysicsUtilities::deltaR(genMuon,*recoMuon->track())));
        tree.fill(global_pt              ,float(recoMuon->globalTrack()->pt()));
        tree.fill(global_eta             ,float(recoMuon->globalTrack()->eta()));
        tree.fill(global_q               ,ASTypes::int8(recoMuon->globalTrack()->charge()));
        tree.fill(tracker_pt             ,float(recoMuon->innerTrack()->pt()));
        tree.fill(tracker_eta            ,float(recoMuon->innerTrack()->eta()));
        tree.fill(tracker_q              ,ASTypes::int8(recoMuon->innerTrack()->charge()));
        tree.fill(standalone_pt          ,float(recoMuon->standAloneMuon()->pt()));
        tree.fill(standalone_eta         ,float(recoMuon->standAloneMuon()->eta()));
        tree.fill(standalone_q           ,ASTypes::int8(recoMuon->standAloneMuon()->charge()));
        tree.fill(picky_pt               ,float(recoMuon->pickyTrack().isNonnull() ? recoMuon->pickyTrack()->pt() : -1.0));
        tree.fill(picky_eta              ,float(recoMuon->pickyTrack().isNonnull() ? recoMuon->pickyTrack()->eta() : -1.0));
        tree.fill(picky_q                ,ASTypes::int8(recoMuon->pickyTrack().isNonnull() ? recoMuon->pickyTrack()->charge() : 0));


        auto getMatchedTrack =[&](const reco::TrackToTrackMap& trackMap) ->  reco::TrackRef  {
            reco::TrackToTrackMap::const_iterator tkRef =trackMap.find(recoMuon->combinedMuon());
            if (tkRef != trackMap.end()) return tkRef->val;
            return reco::TrackRef();
        };

        auto picky_def_trk           = getMatchedTrack(picky_def      );
        auto picky_noGE21_trk        = getMatchedTrack(picky_noGE21);
        auto picky_noME21_trk        = getMatchedTrack(picky_noME21);
        auto picky_noGE11ME0_trk     = getMatchedTrack(picky_noGE11ME0);
        auto picky_noGE11ME0ME21_trk = getMatchedTrack(picky_noGE11ME0ME21);
        auto picky_oneGE21_trk       = getMatchedTrack(picky_oneGE21);
//
//                std::cout << "Start >> "<<std::endl;
//                auto printTracks =[&](const reco::TrackRef& trk){
//                    if(!trk.isAvailable()) return;
//
//                    for(auto rh = trk->recHitsBegin(); rh != trk->recHitsEnd(); ++rh){
//                        if(!(**rh).isValid()) continue;
//                        DetId id = (**rh).geographicalId();
//                        if(id.det() != DetId::Muon) continue;
//
//                        int station = -999;
//                        int ring    = -999;
//                        if ( id.subdetId() == MuonSubdetId::DT ) {
//                            DTChamberId did(id.rawId());
//                            station = did.station();
//                            //UNUSED:   wheel = did.wheel();
//                        } else if ( id.subdetId() == MuonSubdetId::CSC ) {
//                            CSCDetId did(id.rawId());
//                            station = did.station();
//                            ring = did.ring();
//                        } else if ( id.subdetId() == MuonSubdetId::GEM ) {
//                            GEMDetId did(id.rawId());
//                            station = did.station();
//                            ring = did.ring();
//                        } else if ( id.subdetId() == MuonSubdetId::ME0 ) {
//                            ME0DetId did(id.rawId());
//                            station = did.station();
//                        }
//
//                        auto nrh = theMuonRecHitBuilder->build(&**rh);
//
//                        std::cout << "("<< id.subdetId()  <<","<<station<<","<<ring<<") --> ";
////                        if((*nrh).hasPositionAndError() ) std::cout << "["<< (*nrh).globalPositionError().phierr((*nrh).globalPosition()) <<".."<< (*nrh).globalPosition().eta()<<","<<(*nrh).globalPosition().phi()<<"..."<< (*nrh).globalPositionError().rerr((*nrh).globalPosition())<<"]";
//                        if((*nrh).hasPositionAndError() ) std::cout << "["<< (*nrh).globalPosition().phi()<<"]";
//                        std::cout <<" :: ";
//                    }
//                    std::cout << std::endl;
//
//
//
//                };

//                printTracks(recoMuon->globalTrack());
        //        printTracks(recoMuon->pickyTrack());

        auto addHits = [&](const reco::TrackRef& trk,const ASTypes::size v_d,const ASTypes::size v_s,const ASTypes::size v_r,const ASTypes::size v_phie){
            if(!trk.isAvailable()) return;
            for(auto rh = trk->recHitsBegin(); rh != trk->recHitsEnd(); ++rh){
                if(!(**rh).isValid()) continue;
                DetId id = (**rh).geographicalId();
                if(id.det() != DetId::Muon) continue;

                int station = -1;
                int ring    = -1;
                if ( id.subdetId() == MuonSubdetId::DT ) {
                    DTChamberId did(id.rawId());
                    station = did.station();
                    //UNUSED:   wheel = did.wheel();
                } else if ( id.subdetId() == MuonSubdetId::CSC ) {
                    CSCDetId did(id.rawId());
                    station = did.station();
                    ring = did.ring();
                } else if ( id.subdetId() == MuonSubdetId::GEM ) {
                    GEMDetId did(id.rawId());
                    station = did.station();
                    ring = did.ring();
                } else if ( id.subdetId() == MuonSubdetId::ME0 ) {
                    ME0DetId did(id.rawId());
                    station = did.station();
                }
                auto nrh = theMuonRecHitBuilder->build(&**rh);

                tree.fillMulti(v_d       ,ASTypes::convertTo<ASTypes::int8>(id.subdetId(),"subdet"));
                tree.fillMulti(v_s       ,ASTypes::convertTo<ASTypes::int8>(station,"station"));
                tree.fillMulti(v_r       ,ASTypes::convertTo<ASTypes::int8>(ring,"ring"));
                tree.fillMulti(v_phie    ,float((*nrh).globalPositionError().phierr((*nrh).globalPosition())));
            }
        };

        addHits(recoMuon->globalTrack(),global_v_d,global_v_s,global_v_r,global_v_phie);
        addHits(recoMuon->pickyTrack(),picky_v_d,picky_v_s,picky_v_r,picky_v_phie);

        tree.fill(picky_oneGE21_pt       ,float(picky_oneGE21_trk.isNonnull() ? picky_oneGE21_trk->pt() : -1.0));
        tree.fill(picky_oneGE21_eta      ,float(picky_oneGE21_trk.isNonnull() ? picky_oneGE21_trk->eta() : -1.0));
        tree.fill(picky_oneGE21_q        ,ASTypes::int8(picky_oneGE21_trk.isNonnull() ? picky_oneGE21_trk->charge() : 0));

        tree.fill(picky_def_pt           ,float(picky_def_trk.isNonnull() ? picky_def_trk->pt() : -1.0));
        tree.fill(picky_def_eta          ,float(picky_def_trk.isNonnull() ? picky_def_trk->eta() : -1.0));
        tree.fill(picky_def_q            ,ASTypes::int8(picky_def_trk.isNonnull() ? picky_def_trk->charge() : 0));

        tree.fill(picky_noGE21_pt        ,float(picky_noGE21_trk.isNonnull() ? picky_noGE21_trk->pt() : -1.0));
        tree.fill(picky_noGE21_eta       ,float(picky_noGE21_trk.isNonnull() ? picky_noGE21_trk->eta() : -1.0));
        tree.fill(picky_noGE21_q         ,ASTypes::int8(picky_noGE21_trk.isNonnull() ? picky_noGE21_trk->charge() : 0));

        tree.fill(picky_noME21_pt        ,float(picky_noME21_trk.isNonnull() ? picky_noME21_trk->pt() : -1.0));
        tree.fill(picky_noME21_eta       ,float(picky_noME21_trk.isNonnull() ? picky_noME21_trk->eta() : -1.0));
        tree.fill(picky_noME21_q         ,ASTypes::int8(picky_noME21_trk.isNonnull() ? picky_noME21_trk->charge() : 0));

        tree.fill(picky_noGE11ME0_pt     ,float(picky_noGE11ME0_trk.isNonnull() ? picky_noGE11ME0_trk->pt() : -1.0));
        tree.fill(picky_noGE11ME0_eta    ,float(picky_noGE11ME0_trk.isNonnull() ? picky_noGE11ME0_trk->eta() : -1.0));
        tree.fill(picky_noGE11ME0_q      ,ASTypes::int8(picky_noGE11ME0_trk.isNonnull() ? picky_noGE11ME0_trk->charge() : 0));

        tree.fill(picky_noGE11ME0ME21_pt ,float(picky_noGE11ME0ME21_trk.isNonnull() ? picky_noGE11ME0ME21_trk->pt() : -1.0));
        tree.fill(picky_noGE11ME0ME21_eta,float(picky_noGE11ME0ME21_trk.isNonnull() ? picky_noGE11ME0ME21_trk->eta() : -1.0));
        tree.fill(picky_noGE11ME0ME21_q  ,ASTypes::int8(picky_noGE11ME0ME21_trk.isNonnull() ? picky_noGE11ME0ME21_trk->charge() : 0));
    }

private:
    edm::EDGetTokenT<std::vector<reco::GenParticle> >      genparticle_token  ;
    edm::EDGetTokenT<std::vector<reco::Muon> >             recomuon_token  ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_def_token ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_noGE21_token ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_noME21_token ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_noGE11ME0_token ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_noGE11ME0ME21_token ;
    edm::EDGetTokenT<reco::TrackToTrackMap   >             picky_oneGE21_token ;

    TString outFileName;

    TreeWriter tree;
    ASTypes::size nM_all                        =-1;
    ASTypes::size nM_good                       =-1;
    ASTypes::size genMuon_pt                    =-1;
    ASTypes::size genMuon_eta                   =-1;
    ASTypes::size genMuon_q                     =-1;
    ASTypes::size genReco_dr                    =-1;
    ASTypes::size global_pt                     =-1;
    ASTypes::size global_eta                    =-1;
    ASTypes::size global_q                      =-1;
    ASTypes::size global_v_d                    = -1;
    ASTypes::size global_v_s                    = -1;
    ASTypes::size global_v_r                    = -1;
    ASTypes::size global_v_phie                 = -1;
    ASTypes::size tracker_pt                    =-1;
    ASTypes::size tracker_eta                   =-1;
    ASTypes::size tracker_q                     =-1;
    ASTypes::size standalone_pt                 =-1;
    ASTypes::size standalone_eta                =-1;
    ASTypes::size standalone_q                  =-1;
    ASTypes::size picky_pt                      =-1;
    ASTypes::size picky_eta                     =-1;
    ASTypes::size picky_q                       =-1;
    ASTypes::size picky_v_d                    = -1;
    ASTypes::size picky_v_s                    = -1;
    ASTypes::size picky_v_r                    = -1;
    ASTypes::size picky_v_phie                 = -1;
    ASTypes::size picky_oneGE21_pt              =-1;
    ASTypes::size picky_oneGE21_eta             =-1;
    ASTypes::size picky_oneGE21_q               =-1;
    ASTypes::size picky_def_pt                  =-1;
    ASTypes::size picky_def_eta                 =-1;
    ASTypes::size picky_def_q                   =-1;
    ASTypes::size picky_noGE21_pt               =-1;
    ASTypes::size picky_noGE21_eta              =-1;
    ASTypes::size picky_noGE21_q                =-1;
    ASTypes::size picky_noME21_pt               =-1;
    ASTypes::size picky_noME21_eta              =-1;
    ASTypes::size picky_noME21_q                =-1;
    ASTypes::size picky_noGE11ME0_pt            =-1;
    ASTypes::size picky_noGE11ME0_eta           =-1;
    ASTypes::size picky_noGE11ME0_q             =-1;
    ASTypes::size picky_noGE11ME0ME21_pt        =-1;
    ASTypes::size picky_noGE11ME0ME21_eta       =-1;
    ASTypes::size picky_noGE11ME0ME21_q         =-1;

    edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;
};




void
GE21MomentumTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    edm::Handle<std::vector<reco::GenParticle> >  genparticle_han;
    iEvent.getByToken(genparticle_token,genparticle_han);

    edm::Handle<std::vector<reco::Muon> >  recomuon_han;
    iEvent.getByToken(recomuon_token,recomuon_han);

    edm::Handle<reco::TrackToTrackMap >  picky_def_han;
    iEvent.getByToken(picky_def_token,picky_def_han);

    edm::Handle<reco::TrackToTrackMap >  picky_noGE21_han;
    iEvent.getByToken(picky_noGE21_token,picky_noGE21_han);

    edm::Handle<reco::TrackToTrackMap >  picky_noME21_han;
    iEvent.getByToken(picky_noME21_token,picky_noME21_han);

    edm::Handle<reco::TrackToTrackMap >  picky_noGE11ME0_han;
    iEvent.getByToken(picky_noGE11ME0_token,picky_noGE11ME0_han);

    edm::Handle<reco::TrackToTrackMap >  picky_noGE11ME0ME21_han;
    iEvent.getByToken(picky_noGE11ME0ME21_token,picky_noGE11ME0ME21_han);

    edm::Handle<reco::TrackToTrackMap >  picky_oneGE21_han;
    iEvent.getByToken(picky_oneGE21_token,picky_oneGE21_han);


    iSetup.get<TransientRecHitRecord>().get("MuonRecHitBuilder",theMuonRecHitBuilder);


    std::vector<const reco::Muon*> goodMuons;
    for(const auto& rm : *recomuon_han){
        if(!isGoodMuon(rm)) continue;
        goodMuons.push_back(&rm);
    }

    //    ParticleInfo::printGenInfo(*genparticle_han,-1);
    for(const auto& gp : *genparticle_han ){
        if(std::fabs(gp.pdgId()) != ParticleInfo::p_muminus) continue;
        if(!ParticleUtilities::isFirstInChain(&gp)) continue;
        //        if(!ParticleInfo::isDoc(gp.status())) continue;
        if(gp.numberOfMothers()) continue; //for pythia
        const auto* finalP = ParticleUtilities::getFinal(&gp);

        double nearDR2 = 99999;
        const reco::Muon* matched = 0;
        for(const auto * rm : goodMuons){
            const float dr2 = PhysicsUtilities::deltaR2(*finalP,*rm->innerTrack());
            if(dr2 >= nearDR2) continue;
            matched = rm;
            nearDR2 = dr2;
        }


        tree.reset();
        tree.fill(nM_all             ,ASTypes::size(recomuon_han->size()));
        tree.fill(nM_good            ,ASTypes::size(goodMuons.size()));
        fillTree(*finalP,matched==0 ? 0 : matched ,*picky_def_han,*picky_noGE21_han,*picky_noME21_han,*picky_noGE11ME0_han,*picky_noGE11ME0ME21_han,*picky_oneGE21_han);
        tree.fillTree();

    }



}



//define this as a plug-in
DEFINE_FWK_MODULE(GE21MomentumTreeMaker);
