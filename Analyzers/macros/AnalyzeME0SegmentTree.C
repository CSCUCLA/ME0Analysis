
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/PhysicsUtilities.h"
#include "HistoPlotting/include/Plotter.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TGraph.h"
#include "TFitter.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
using namespace std;
using namespace PhysicsUtilities;
using ASTypes::CylLorentzVectorF;
typedef std::vector<std::tuple<double,double,double>>  MeasureList;//pt DPhi eta

class Tau3MuAnalyzer{
public:
    TFile * file = 0;
    TTreeReader* reader = 0;
    TTreeReaderArray<double> * mu_pt  = 0;
    TTreeReaderArray<double> * mu_eta = 0;
    TTreeReaderArray<double> * mu_phi = 0;
    TTreeReaderArray<double> * mu_rec = 0;
    int nEvents = 0;
    TRandom3 * rand = 0;

    Tau3MuAnalyzer(TString filename ="tree_tau3mu_gen.root"){
        file = new TFile(filename,"Read");
        reader = new TTreeReader("demo/data", file);
        mu_pt = new TTreeReaderArray<double> (*reader, "mu_pt");
        mu_eta= new TTreeReaderArray<double> (*reader, "mu_eta");
        mu_phi= new TTreeReaderArray<double> (*reader, "mu_phi");
        mu_rec= new TTreeReaderArray<double> (*reader, "mu_rec");
        nEvents = reader->GetEntries(true);
        rand = new TRandom3(1234);
    }

    CylLorentzVectorF getRandomCentralDiMu() {
        while(true){
            int entry = rand->Uniform(0,nEvents);
            reader->SetEntry(entry);
            int nLT2p4 = 0;
            int n2p4to2p8 = 0;
            for(unsigned int iM = 0; iM < 3; ++iM){
                if (std::fabs((*mu_eta)[iM]) < 2.4 && (*mu_pt)[iM] > 2) nLT2p4++;
                else if (std::fabs((*mu_eta)[iM]) > 2.4 && std::fabs((*mu_eta)[iM]) < 2.8) n2p4to2p8++;
            }
            if(nLT2p4 != 2 || n2p4to2p8 != 1) continue;
            CylLorentzVectorF innerMuons;
            for(unsigned int iM = 0; iM < 3; ++iM){
                CylLorentzVectorF muon((*mu_pt)[iM],(*mu_eta)[iM],(*mu_phi)[iM],.105);
                if(std::fabs(muon.eta()) < 2.4)   innerMuons += muon;
            }
            return innerMuons;
        }
    }

    CylLorentzVectorF getRandomCentralSingleMu() {
        while(true){
            int entry = rand->Uniform(0,nEvents);
            reader->SetEntry(entry);
            int nLT2p4 = 0;
            int n2p4to2p8 = 0;
            for(unsigned int iM = 0; iM < 3; ++iM){
                if (std::fabs((*mu_eta)[iM]) < 2.4 && (*mu_pt)[iM] > 2) nLT2p4++;
                else if (std::fabs((*mu_eta)[iM]) > 2.4 && std::fabs((*mu_eta)[iM]) < 2.8) n2p4to2p8++;
            }
            if(nLT2p4 != 1 || n2p4to2p8 != 2) continue;
            CylLorentzVectorF innerMuons;
            for(unsigned int iM = 0; iM < 3; ++iM){
                CylLorentzVectorF muon((*mu_pt)[iM],(*mu_eta)[iM],(*mu_phi)[iM],.105);
                if(std::fabs(muon.eta()) < 2.4)   innerMuons += muon;
            }
            return innerMuons;
        }
    }

    ~Tau3MuAnalyzer(){
        if(file) file->Close();
        if(reader) delete reader;
        if(mu_pt ) delete mu_pt ;
        if(mu_eta) delete mu_eta;
        if(mu_phi) delete mu_phi;
        if(mu_rec) delete mu_rec;
        delete rand;
    }


};

class Analyzer : public BaseTupleAnalyzer{
public:
    TString glbPrefix = "";

    float minDPhi = 1;
    int   nStrips = 1;

    Analyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){
        setBranchAddress("simMuon_pt"               ,&simMuon_pt               );
        setBranchAddress("simMuon_eta"              ,&simMuon_eta              );
        setBranchAddress("simMuon_phi"              ,&simMuon_phi              );
        setBranchAddress("simMuon_q"                ,&simMuon_q                );
        setBranchAddress("simMuon_gen_eta"          ,&simMuon_gen_eta          );
        setBranchAddress("simMuon_gen_phi"          ,&simMuon_gen_phi          );
        setBranchAddress("simMuon_gen_dphi"         ,&simMuon_gen_dphi         );
        setBranchAddress("simMuon_gen_deta"         ,&simMuon_gen_deta         );
        setBranchAddress("simMuon_gen_x"            ,&simMuon_gen_x            );
        setBranchAddress("simMuon_gen_y"            ,&simMuon_gen_y            );
        setBranchAddress("simMuon_gen_dx"           ,&simMuon_gen_dx           );
        setBranchAddress("simMuon_gen_dy"           ,&simMuon_gen_dy           );
        setBranchAddress("simMuon_segmentIDX"       ,&simMuon_segmentIDX       );
        setBranchAddress("simMuon_segment_quality"  ,&simMuon_segment_quality  );
        setBranchAddress("simMuon_segment_nGoodHits",&simMuon_segment_nGoodHits);
        setBranchAddress("simMuon_segment_nBadHits" ,&simMuon_segment_nBadHits );
        setBranchAddress("segment_eta"              ,&segment_eta              );
        setBranchAddress("segment_phi"              ,&segment_phi              );
        setBranchAddress("segment_dphi"             ,&segment_dphi             );
        setBranchAddress("segment_deta"             ,&segment_deta             );
        setBranchAddress("segment_x"                ,&segment_x                );
        setBranchAddress("segment_y"                ,&segment_y                );
        setBranchAddress("segment_dx"               ,&segment_dx               );
        setBranchAddress("segment_dy"               ,&segment_dy               );
    }


    std::vector<float>* simMuon_pt               = new vector<float>;
    std::vector<float>* simMuon_eta              = new vector<float>;
    std::vector<float>* simMuon_phi              = new vector<float>;
    std::vector<int>  * simMuon_q                = new vector<int>  ;

    std::vector<float>* simMuon_gen_eta          = new vector<float>;
    std::vector<float>* simMuon_gen_phi          = new vector<float>;
    std::vector<float>* simMuon_gen_dphi         = new vector<float>;
    std::vector<float>* simMuon_gen_deta         = new vector<float>;
    std::vector<float>* simMuon_gen_x            = new vector<float>;
    std::vector<float>* simMuon_gen_y            = new vector<float>;
    std::vector<float>* simMuon_gen_dx           = new vector<float>;
    std::vector<float>* simMuon_gen_dy           = new vector<float>;

    std::vector<int>  * simMuon_segmentIDX       = new vector<int>  ;
    std::vector<int>  * simMuon_segment_quality  = new vector<int>  ;
    std::vector<int>  * simMuon_segment_nGoodHits= new vector<int>  ;
    std::vector<int>  * simMuon_segment_nBadHits = new vector<int>  ;
    std::vector<float>* segment_eta              = new vector<float>;
    std::vector<float>* segment_phi              = new vector<float>;
    std::vector<float>* segment_dphi             = new vector<float>;
    std::vector<float>* segment_deta             = new vector<float>;
    std::vector<float>* segment_x                = new vector<float>;
    std::vector<float>* segment_y                = new vector<float>;
    std::vector<float>* segment_dx               = new vector<float>;
    std::vector<float>* segment_dy               = new vector<float>;

    void write(TString fileName){ plotter.write(fileName);}
    HistGetter plotter;

    virtual void runAEvent() {
    }
};

class FindParams : public Analyzer {
public:

    MeasureList genEtaList;
    MeasureList genThetaList;
    MeasureList genSinThetaList;
    MeasureList recoEtaList;
    MeasureList recoThetaList;
    MeasureList recoSinThetaList;
    FindParams(std::string fileName, std::string treeName) : Analyzer(fileName,treeName) {}

    virtual void runAEvent() {
        for(unsigned int iS = 0; iS < simMuon_pt->size(); ++iS){
            CylLorentzVectorF muon(simMuon_pt->at(iS),simMuon_eta->at(iS),simMuon_phi->at(iS),0.105);
            if(muon.pt() < 2) continue;
            double genTheta = muon.theta() < TMath::PiOver2() ?  muon.theta() : TMath::Pi() - muon.theta() ;
            genEtaList     .emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),std::fabs(muon.eta()) );
            genThetaList   .emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),genTheta );
            genSinThetaList.emplace_back(muon.pt(),std::fabs(simMuon_gen_dphi->at(iS)),std::sin(genTheta) );
            if(simMuon_segment_quality->at(iS) == 3) continue;
            const int idx = simMuon_segmentIDX->at(iS);
            const double theta = PhysicsUtilities::etaToTheta(std::fabs(segment_eta->at(idx)));
            recoEtaList     .emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),std::fabs(segment_eta->at(idx)) );
            recoThetaList   .emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),theta );
            recoSinThetaList.emplace_back(muon.pt(),std::fabs(segment_dphi->at(idx)),std::sin(theta) );

            auto solPT = [&](float dPhi, float eta, float a, float b, bool gen = false) -> double { return (a + b*eta)/std::max(dPhi,(gen? float(0.0) :minDPhi)); };

            double genSolPT = solPT(std::fabs(simMuon_gen_dphi->at(iS)),genTheta,-0.0039,0.134,true );
            double recoSolPT = solPT(std::fabs(segment_dphi->at(idx)),theta,-0.0039,0.134   );

            double genSolPT2 = solPT(std::fabs(simMuon_gen_dphi->at(iS)),std::fabs(muon.eta()),0.07292,-0.02116,true    );
            double recoSolPT2 =  solPT(std::fabs(segment_dphi->at(idx)),std::fabs(segment_eta->at(idx)),0.07292,-0.02116  );

            plotter.getOrMake1D("gen_eta_PToPT",";reco pt / gen pt"       ,100,0,10)->Fill(  genSolPT2/muon.pt()  );
            plotter.getOrMake1D("gen_theta_PToPT",";reco pt / gen pt"     ,100,0,10)->Fill(  genSolPT/muon.pt()  );
            plotter.getOrMake1D("reco_eta_PToPT",";reco pt / gen pt"      ,100,0,10)->Fill(  recoSolPT2/muon.pt()  );
            plotter.getOrMake1D("reco_theta_PToPT",";reco pt / gen pt",100,0,10)->Fill(  recoSolPT/muon.pt()  );



            TString selString = "";
            if(muon.pt() < 5) selString = "pt1to5_";
            else if(muon.pt() < 15) selString = "pt5to15_";
            else if(muon.pt() < 18) selString = "pt15to18_";
            else selString = "ptgeq18_";
            //			if(std::fabs(muon.eta()) < 2.2) selString += "eta2to2p2";
            //			else if(std::fabs(muon.eta()) < 2.4) selString += "eta2p2to2p4";
            //			else if(std::fabs(muon.eta()) < 2.6) selString += "eta2p2to2p6";
            //			else selString += "eta2p6to2p8";

            if(std::fabs(muon.eta()) < 2.4) selString += "eta2to2p2";
            else selString += "eta2p6to2p8";


            plotter.getOrMake1D(TString::Format("gen_%s_PToPT",selString.Data()),";segment p_{T} / gen p_{T}",100,0,4)->Fill(  genSolPT/muon.pt()  );
            plotter.getOrMake1D(TString::Format("reco_%s_PToPT",selString.Data()),";segment p_{T} / gen p_{T}",100,0,4)->Fill(  recoSolPT/muon.pt()  );

            plotter.getOrMake2D(TString::Format("genPT_v_recoPT"),";gen p_{T};segment p_{T}",20,0,60,20,0,60)->Fill(  muon.pt(), recoSolPT  );

            plotter.getOrMake2D(TString::Format("reco_PToPT_v_nStrips"),";segment p_{T} / gen p_{T};dPhi [nStrips]",50,0,10,100,0,10)->Fill(  solPT(std::fabs(segment_dphi->at(idx)),theta,-0.0039,0.134,true)/muon.pt(),
                    std::fabs(segment_dphi->at(idx))/(20.0*TMath::Pi()/(180.0*384))  );



        }
    }

    void solve() {

        TFile * file = TFile::Open("profiles.root","RECREATE");
        auto makeSolution = [&](const MeasureList& list, TString name, int nBins, double min, double max){
            TProfile * n = new TProfile(TString::Format("%s_N",name.Data()),";|eta|",nBins,min,max);
            for(const auto& val : list ){
                if(std::get<1>(val) > .1) continue;
                if(std::get<1>(val) < 0.0008) continue;
                if(std::get<0>(val) > 10) continue;
                n->Fill(std::get<2>(val) ,  1/(std::get<0>(val)*std::get<1>(val)));
            }

            file->cd();
            n->Write();
        };
        makeSolution(genEtaList,"genEtaList",32,2,2.8);
        makeSolution(genThetaList,"genThetaList",32,.121,.269);
        makeSolution(genSinThetaList,"genSinThetaList",32,.121,.269);
        makeSolution(recoEtaList,"recoEtaList",32,2,2.8);
        makeSolution(recoThetaList,"recoThetaList",32,.121,.269);
        makeSolution(recoSinThetaList,"recoSinThetaList",32,.121,.269);
        file->Close();

    }

};

class TestParams : public Analyzer {
public:
    TRandom3 * rand = new TRandom3(12345);
    double getMeasuredPT(float dPhi, float theta, bool gen=false, float a = -0.0039, float b = 0.134 )const { return (a + b*theta)/std::max(dPhi,(gen? float(0.0) :minDPhi));}
    CylLorentzVectorF getMuon(unsigned int idx) const { return CylLorentzVectorF(simMuon_pt->at(idx),simMuon_eta->at(idx),simMuon_phi->at(idx),0.105);}
    CylLorentzVectorF getMeasuredGen(unsigned int idx) const {
        auto muon = getMuon(idx);
        double genTheta = muon.theta() < TMath::PiOver2() ?  muon.theta() : TMath::Pi() - muon.theta() ;
        return CylLorentzVectorF(getMeasuredPT(std::fabs(simMuon_gen_dphi->at(idx)), genTheta,true  ),simMuon_eta->at(idx),simMuon_phi->at(idx),0.105);
    }
    CylLorentzVectorF getMeasuredReco(unsigned int idx) const {
        const double theta = PhysicsUtilities::etaToTheta(std::fabs(segment_eta->at(idx)))   ;
        return CylLorentzVectorF(getMeasuredPT(std::fabs(segment_dphi->at(idx)), theta  ),segment_eta->at(idx),segment_phi->at(idx),0.105 );
    }
    int getRecoCharge(unsigned int idx) const {return segment_dphi->at(idx) > 0 ? -1 : 1;};
    int getGenCharge(unsigned int idx) const {return simMuon_gen_dphi->at(idx) > 0 ? -1 : 1;};
    CylLorentzVectorF getRandomTrack(const double& pt) { return CylLorentzVectorF(pt,rand->Uniform(-2.4,2.4),rand->Uniform(-TMath::Pi(),TMath::Pi()),.105);}

    TestParams(std::string fileName, std::string treeName) : Analyzer(fileName,treeName) {}

    void tau3MuStudies() {
        int muIDX1 = -1;
        int muIDX2 = -1;
        double bestMass = -1;
        for(unsigned int iM1 = 0; iM1 < simMuon_pt->size(); ++iM1){
            auto mu1 = getMuon(iM1);
            for(unsigned int iM2 = iM1+1; iM2 < simMuon_pt->size(); ++iM2){
                auto mu2 = getMuon(iM2);
                if(simMuon_q->at(iM1) == simMuon_q->at(iM2)) continue;
                double mass = (mu1+mu2).mass();
                if(bestMass < 0 || abs(mass - 91.1876)  <  abs(bestMass - 91.1876) ){
                    bestMass = mass;
                    muIDX1 = iM1; muIDX2 = iM2;
                }
            }
        }
        if(muIDX1 < 0 || muIDX2 < 0 ) return;
        const auto mu1 = getMuon(muIDX1); const float mu1AbsEta = std::fabs(mu1.eta());
        const auto mu2 = getMuon(muIDX2); const float mu2AbsEta = std::fabs(mu2.eta());

        if ((mu1AbsEta > 2.0) || (mu2AbsEta > 2.0) ) return;

        plotter.getOrMake1D(TString::Format("%s_tau3mu_nFakeEvents",glbPrefix.Data()),";nEvents",1,0,2)->Fill( 1 );

        //do 2 muons innner
        auto cenralDi = getRandomTrack(5);// tau3MuTree.getRandomCentralDiMu();
        bool foundTwoInnerPT1 = false;
        bool foundTwoInnerPT2 = false;
        CylLorentzVectorF minTwoInnerPT1;
        CylLorentzVectorF minTwoInnerPT2;
        plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_pt",glbPrefix.Data()),";di-muon inner p_{T} [GeV]",200,0,200)->Fill(cenralDi.pt() );
        plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_abseta",glbPrefix.Data()),";di-muon inner p_{T} [GeV]",50,0,5)->Fill(std::fabs(cenralDi.eta()) );
        for(unsigned int iS = 0; iS < segment_dx->size(); ++iS){
            auto fakeSeg = getMeasuredReco(iS);
            if(std::fabs(fakeSeg.eta()) < 2.4) continue;
            if(fakeSeg.pt() > 1){
                plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_pt1_bkg_invmass",glbPrefix.Data()),";2 inner muon + ME0Segment mass [GeV]",40,0,20)->Fill(  (cenralDi + fakeSeg).mass() );
                if(!foundTwoInnerPT1 || (cenralDi + fakeSeg).mass() < minTwoInnerPT1.mass())
                    minTwoInnerPT1 = (cenralDi + fakeSeg);
                foundTwoInnerPT1 =true;
            }

            if(fakeSeg.pt() > 2){
                plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_pt2_bkg_invmass",glbPrefix.Data()),";2 inner muon + ME0Segment mass [GeV]",40,0,20)->Fill(  (cenralDi + fakeSeg).mass() );
                if(!foundTwoInnerPT2 || (cenralDi + fakeSeg).mass() < minTwoInnerPT2.mass())
                    minTwoInnerPT2 = (cenralDi + fakeSeg);
                foundTwoInnerPT2 =true;
            }
        }
        if(foundTwoInnerPT1) plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_pt1_bkg_min_invmass",glbPrefix.Data()),";reconstructed #tau mass [GeV]",80,0,20)->Fill(  minTwoInnerPT1.mass() );
        if(foundTwoInnerPT2) plotter.getOrMake1D(TString::Format("%s_tau3mu_twoInner_pt2_bkg_min_invmass",glbPrefix.Data()),";reconstructed #tau mass [GeV]",80,0,20)->Fill(  minTwoInnerPT2.mass() );


        //Now do one inner
        auto singleDi = getRandomTrack(5);//tau3MuTree.getRandomCentralSingleMu();
        plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_pt",glbPrefix.Data()),";single-muon p_{T} [GeV]",40,0,20)->Fill(singleDi.pt() );
        plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_abseta",glbPrefix.Data()),";single-muon p_{T} [GeV]",50,0,5)->Fill(std::fabs(singleDi.eta()) );
        bool foundOneInnerPT1 = false;
        bool foundOneInnerPT2 = false;
        CylLorentzVectorF minOneInnerPT1;
        CylLorentzVectorF minOneInnerPT2;
        for(unsigned int iS1 = 0; iS1 < segment_dx->size(); ++iS1){
            auto fakeSeg1 = getMeasuredReco(iS1);
            if(std::fabs(fakeSeg1.eta()) < 2.4) continue;
            for(unsigned int iS2 = iS1+1; iS2 < segment_dx->size(); ++iS2){
                auto fakeSeg2 = getMeasuredReco(iS2);
                if(std::fabs(fakeSeg2.eta()) < 2.4) continue;
                if(fakeSeg1.pt() > 1 && fakeSeg2.pt() > 1){
                    plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_pt1_bkg_invmass",glbPrefix.Data()),";inner muon + 2 ME0Segments mass [GeV]",40,0,20)->Fill(  (singleDi + fakeSeg1 +fakeSeg2).mass() );
                    if(!foundOneInnerPT1 || (singleDi + fakeSeg1 + fakeSeg2).mass() < minOneInnerPT1.mass())
                        minOneInnerPT1 = (singleDi + fakeSeg1 + fakeSeg2);
                    foundOneInnerPT1 =true;
                }
                if(fakeSeg1.pt() > 2 && fakeSeg2.pt() > 2){
                    plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_pt2_bkg_invmass",glbPrefix.Data()),";inner muon + 2 ME0Segments mass [GeV]",40,0,20)->Fill(  (singleDi + fakeSeg1 +fakeSeg2).mass() );
                    if(!foundOneInnerPT2 || (singleDi + fakeSeg1 + fakeSeg2).mass() < minOneInnerPT2.mass())
                        minOneInnerPT2 = (singleDi + fakeSeg1 + fakeSeg2);
                    foundOneInnerPT2 =true;
                }
            }
        }
        if(foundOneInnerPT1) plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_pt1_bkg_min_invmass",glbPrefix.Data()),";reconstructed #tau mass [GeV]",40,0,20)->Fill(  minOneInnerPT1.mass() );
        if(foundOneInnerPT2) plotter.getOrMake1D(TString::Format("%s_tau3mu_oneInner_pt2_bkg_min_invmass",glbPrefix.Data()),";reconstructed #tau mass [GeV]",40,0,20)->Fill(  minOneInnerPT2.mass() );

        //		//Now do three inner
        //		for(unsigned int iS1 = 0; iS1 < segment_dx->size(); ++iS1){
        //			auto fakeSeg1 = getMeasuredReco(iS1);
        //			if(std::fabs(fakeSeg1.eta()) < 2.4) continue;
        //			for(unsigned int iS2 = iS1+1; iS2 < segment_dx->size(); ++iS2){
        //				auto fakeSeg2 = getMeasuredReco(iS2);
        //				if(std::fabs(fakeSeg2.eta()) < 2.4) continue;
        //				for(unsigned int iS3 = iS2+1; iS3 < segment_dx->size(); ++iS3){
        //					auto fakeSeg3 = getMeasuredReco(iS3);
        //					if(std::fabs(fakeSeg3.eta()) < 2.4) continue;
        //					plotter.getOrMake1D(TString::Format("%s_tau3mu_zeroInner_bkg_invmass",glbPrefix.Data()),";3 ME0Segments mass [GeV]",200,0,200)->Fill(  (fakeSeg3 + fakeSeg1 +fakeSeg2).mass() );
        //					if(fakeSeg1.pt() > 1 &&fakeSeg2.pt() > 1 && fakeSeg3.pt() > 1 )
        //						plotter.getOrMake1D(TString::Format("%s_tau3mu_zeroInner_pt1_bkg_invmass",glbPrefix.Data()),";3 ME0Segments mass [GeV]",200,0,200)->Fill(  (fakeSeg3 + fakeSeg1 +fakeSeg2).mass() );
        //					if(fakeSeg1.pt() > 2 &&fakeSeg2.pt() > 2 && fakeSeg3.pt() > 2 )
        //						plotter.getOrMake1D(TString::Format("%s_tau3mu_zeroInner_pt2_bkg_invmass",glbPrefix.Data()),";3 ME0Segments mass [GeV]",200,0,200)->Fill(  (fakeSeg3 + fakeSeg1 +fakeSeg2).mass() );
        //
        //				}
        //			}
        //		}

    }

    void zMuMuStudies() {
        int muIDX1 = -1;
        int muIDX2 = -1;
        double bestMass = -1;
        for(unsigned int iM1 = 0; iM1 < simMuon_pt->size(); ++iM1){
            auto mu1 = getMuon(iM1);
            for(unsigned int iM2 = iM1+1; iM2 < simMuon_pt->size(); ++iM2){
                auto mu2 = getMuon(iM2);
                if(simMuon_q->at(iM1) == simMuon_q->at(iM2)) continue;
                double mass = (mu1+mu2).mass();
                if(bestMass < 0 || abs(mass - 91.1876)  <  abs(bestMass - 91.1876) ){
                    bestMass = mass;
                    muIDX1 = iM1; muIDX2 = iM2;
                }
            }
        }
        if(muIDX1 < 0 ||muIDX2 < 0 ) return;
        const auto mu1 = getMuon(muIDX1); const float mu1AbsEta = std::fabs(mu1.eta());
        const auto mu2 = getMuon(muIDX2); const float mu2AbsEta = std::fabs(mu2.eta());
        //		plotter.getOrMake1D(TString::Format("%s_zmm_types",glbPrefix.Data()),";Incl/pt>2/2 <2.4/1 <2.4 1 2.4-28/2 2.4-28/Other",6,-0.5,5.5)->Fill(0 );

        //		if(mu1.pt() < 2 || mu2.pt() < 2) return;

        //		plotter.getOrMake1D(TString::Format("%s_zmm_types",glbPrefix.Data()),";Incl/pt>2/2 <2.4/1 <2.4 1 2.4-28/2 2.4-28/Other",6,-0.5,5.5)->Fill(1);


        int type = 5;
        int nLt2p4 = (mu1AbsEta < 2.4) + (mu2AbsEta < 2.4);
        int n2p4to2p8 = (mu1AbsEta > 2.4 && mu1AbsEta < 2.8 ) + (mu2AbsEta > 2.4 && mu2AbsEta < 2.8 );
        if(nLt2p4 == 2) type =2;
        if(nLt2p4 ==1 && n2p4to2p8 ==1  ) type =3;
        if(n2p4to2p8 ==2  ) type =4;

        //		plotter.getOrMake1D(TString::Format("%s_zmm_types",glbPrefix.Data()),";Incl/pt>2/2 <2.4/1 <2.4 1 2.4-28/2 2.4-28/Other",6,-0.5,5.5)->Fill(type);

        //pure fakes
        if(type == 2 && (mu1AbsEta < 2.0) && (mu2AbsEta < 2.0) ){
            plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_nFakeEvents",glbPrefix.Data()),";nEvents",1,0,2)->Fill( 1 );
            std::vector<double> pts;
            for(unsigned int iS = 0; iS < segment_dx->size(); ++iS){
                auto fakeSeg = getMeasuredReco(iS);
                if(TMath::Abs(fakeSeg.eta()) < 2.4) continue;
                pts.push_back(fakeSeg.pt());
                plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",200,0,200)->Fill(fakeSeg.pt() );
            }
            std::sort(pts.begin(),pts.end(),[](double a, double b) {return a > b;});
            if(pts.size() > 0)plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_max_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",80,0,40)->Fill(pts[0] );
            if(pts.size() > 1)plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_secmax_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",80,0,40)->Fill( pts[1]);



            //for dimuon
            bool twoMu = false;
            bool twoMuCharge = false;
            bool twoMuChargePT = false;
            bool twoMuChargePTMass = false;
            for(unsigned int iS1 = 0; iS1 < segment_dx->size(); ++iS1){
                auto fakeSeg1 = getMeasuredReco(iS1);
                if(std::fabs(fakeSeg1.eta()) < 2.4) continue;
                for(unsigned int iS2 = iS1+1; iS2 < segment_dx->size(); ++iS2){
                    auto fakeSeg2 = getMeasuredReco(iS2);
                    if(std::fabs(fakeSeg2.eta()) < 2.4) continue;
                    twoMu= true;

                    bool chargeCut = (std::fabs(segment_dphi->at(iS1)) < minDPhi) || (std::fabs(segment_dphi->at(iS2)) < minDPhi) ||
                            (getRecoCharge(iS1) != getRecoCharge(iS2) );
                    bool pTcut   = std::min(fakeSeg1.pt(),fakeSeg2.pt()) > 6;// (nStrips < 200 ? 6.0  : 10.0);
                    bool massCut = ( (fakeSeg1 + fakeSeg2).mass()  > (nStrips < 200 ? 15.0  : 30.0)  ) && ( (fakeSeg1 + fakeSeg2).mass()  < (nStrips < 200 ? 45.0  : 80.0)  );
                    if(chargeCut) twoMuCharge = true;
                    if(chargeCut&&pTcut) twoMuChargePT = true;
                    if(chargeCut&&pTcut&&massCut) twoMuChargePTMass = true;

                    plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_bkg_invmass",glbPrefix.Data()),";2 ME0Segment mass [GeV]",200,0,200)->Fill(  (fakeSeg1 + fakeSeg2).mass() );
                    if(chargeCut)  plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_passCharge_bkg_invmass",glbPrefix.Data()),";2 ME0Segment mass [GeV]",200,0,200)->Fill(  (fakeSeg1 + fakeSeg2).mass() );
                    if(chargeCut && pTcut)  plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_passChargeAndPT_bkg_invmass",glbPrefix.Data()),";2 ME0Segment mass [GeV]",200,0,200)->Fill(  (fakeSeg1 + fakeSeg2).mass() );

                }
            }
            if(twoMu) plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_bkgMult",glbPrefix.Data()),";2 Mu / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 0);
            if(twoMuCharge) plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_bkgMult",glbPrefix.Data()),";2 Mu / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 1);
            if(twoMuChargePT) plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_bkgMult",glbPrefix.Data()),";2 Mu / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 2);
            if(twoMuChargePTMass) plotter.getOrMake1D(TString::Format("%s_zmm_twoInner_diMu_bkgMult",glbPrefix.Data()),";2 Mu / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 3);






        }


        else if(type == 3){
            CylLorentzVectorF innerTrack;
            int genME0MuonIDX = -1;
            bool good = true;
            plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 0);
            if(mu1AbsEta < 2.4){
                innerTrack = mu1;
                if(simMuon_segment_quality->at(muIDX2) < 3 ) {genME0MuonIDX = muIDX2;}
                else good = false;
            }
            else {
                innerTrack = mu2;
                if(simMuon_segment_quality->at(muIDX1) < 3 ){genME0MuonIDX = muIDX1;}
                else good = false;
            }
            if(good){
                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 1);
                bool chargeCut = std::fabs(segment_dphi->at(simMuon_segmentIDX->at(genME0MuonIDX))) < minDPhi ||
                        getRecoCharge(simMuon_segmentIDX->at(genME0MuonIDX)) == simMuon_q->at(genME0MuonIDX);

                if(chargeCut) plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 2);

                auto me0Seg = getMeasuredReco(simMuon_segmentIDX->at(genME0MuonIDX));
                auto genSeg = getMeasuredGen(genME0MuonIDX);

                bool pTcut =me0Seg.pt() > 6;
                bool massCut =(innerTrack + me0Seg).mass()  > 40 && (innerTrack + me0Seg).mass()  < 100;

                if(chargeCut && pTcut)
                    plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 3);
                if(chargeCut && pTcut && massCut)
                    plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 4);

                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass",200,0,200)->Fill(  (innerTrack + me0Seg).mass() );
                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_simhit_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass",200,0,200)->Fill(  (innerTrack + genSeg).mass() );

                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_pt",glbPrefix.Data()),";ME0Segment  p_{T}",200,0,200)->Fill(me0Seg.pt() );
                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_simhit_pt",glbPrefix.Data()),";ME0Segment p_{T}",200,0,200)->Fill(genSeg.pt() );
                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_gen_pt",glbPrefix.Data()),";gen p_{T}",200,0,200)->Fill(simMuon_pt->at(genME0MuonIDX) );

                plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_nFakeEvents",glbPrefix.Data()),";nEvents",1,0,2)->Fill( 1 );
                for(unsigned int iS = 0; iS < segment_dx->size(); ++iS){
                    if(iS == simMuon_segmentIDX->at(genME0MuonIDX)) continue;
                    auto fakeSeg = getMeasuredReco(iS);
                    if(TMath::Abs(fakeSeg.eta()) < 2.4) continue;

                    bool chargeCut = std::fabs(segment_dphi->at(iS)) < minDPhi || getRecoCharge(iS) == simMuon_q->at(genME0MuonIDX);
                    bool pTcut =fakeSeg.pt() > 6;
                    bool massCut =(innerTrack + fakeSeg).mass()  > 40 && (innerTrack + fakeSeg).mass()  < 100;

                    plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_bkgMult",glbPrefix.Data()),";incl / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 0);
                    if(chargeCut) plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_bkgMult",glbPrefix.Data()),";incl / charge / pT cut / mass cut",4,-0.5,3.5)->Fill( 1);
                    if(chargeCut && pTcut ) plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_bkgMult",glbPrefix.Data()),";incl / charge / pT cut / mass cut",4,-0.5,3.5)->Fill(2);
                    if(chargeCut && pTcut && massCut ) plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_bkgMult",glbPrefix.Data()),";incl / charge / pT cut / mass cut",4,-0.5,3.5)->Fill(3);

                    plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_bkg_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass",200,0,200)->Fill(  (innerTrack + fakeSeg).mass() );
                    if(chargeCut)  plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_passCharge_bkg_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass",200,0,200)->Fill(  (innerTrack + fakeSeg).mass() );
                    if(chargeCut && pTcut)  plotter.getOrMake1D(TString::Format("%s_zmm_oneInnerOneME0_passChargeAndPT_bkg_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass",200,0,200)->Fill(  (innerTrack + fakeSeg).mass() );

                }

            }

        } else if(type == 4){
            plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 0);

            bool good = true;
            CylLorentzVectorF me0Seg1, me0Seg2;
            if(simMuon_segment_quality->at(muIDX1) < 3 ) me0Seg1 = getMeasuredReco(simMuon_segmentIDX->at(muIDX1));
            else good = false;
            if(simMuon_segment_quality->at(muIDX2) < 3 ) me0Seg2 = getMeasuredReco(simMuon_segmentIDX->at(muIDX2));
            else good = false;
            if(good){
                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 1);

                auto genSeg1 = getMeasuredGen(muIDX1);
                auto genSeg2 = getMeasuredGen(muIDX2);

                bool chargeCut = (std::fabs(segment_dphi->at(simMuon_segmentIDX->at(muIDX1))) < minDPhi) || (std::fabs(segment_dphi->at(simMuon_segmentIDX->at(muIDX2))) < minDPhi) ||
                        (getRecoCharge(simMuon_segmentIDX->at(muIDX1)) != getRecoCharge(simMuon_segmentIDX->at(muIDX2)) );
                bool pTcut   = std::min(me0Seg1.pt(),me0Seg2.pt()) > 6;// (nStrips < 200 ? 6.0  : 10.0);
                bool massCut = ( (me0Seg1 + me0Seg2).mass()  > (nStrips < 200 ? 15.0  : 30.0)  ) && ( (me0Seg1 + me0Seg2).mass()  < (nStrips < 200 ? 45.0  : 80.0)  );

                if(chargeCut) plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 2);
                if(chargeCut&&pTcut) plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 3);
                if(chargeCut&&pTcut&&massCut) plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_eff",glbPrefix.Data()),";total / seg reco / charge / pT cut / mass cut",5,-0.5,4.5)->Fill( 4);


                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_invmass",glbPrefix.Data()),";2 ME0Segments mass [GeV]",200,0,200)->Fill(  (me0Seg1 + me0Seg2).mass() );
                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_simhit_invmass",glbPrefix.Data()),";gen muon + ME0Segment mass [GeV]",200,0,200)->Fill(  (genSeg1 + genSeg2).mass() );

                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_min_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",200,0,200)->Fill(std::min(me0Seg1.pt(),me0Seg2.pt()) );
                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_max_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",200,0,200)->Fill(std::max(me0Seg1.pt(),me0Seg2.pt()) );
                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_simhit_min_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",200,0,200)->Fill(std::min(genSeg1.pt(),genSeg2.pt()) );
                plotter.getOrMake1D(TString::Format("%s_zmm_twoME0_simhit_max_pt",glbPrefix.Data()),";ME0Segment p_{T} [GeV]",200,0,200)->Fill(std::max(genSeg1.pt(),genSeg2.pt()) );
            }


        }

    }

    void fakeStudies() {


        int nGMuons = 0;
        for(const auto& qual : *simMuon_segment_quality ) if(qual < 3) nGMuons++;
        if(nGMuons != 2) return;
        plotter.getOrMake1D(TString::Format("%s_fs_nEvents",glbPrefix.Data()),";nEvents",1,0,2)->Fill(  1  );

        auto pt5SingleMu  = getRandomTrack(5);
        auto pt20SingleMu  = getRandomTrack(20);


        for(unsigned int idx = 0; idx < segment_eta->size(); ++idx){
            bool good = false; for(const auto& trueIDX : *simMuon_segmentIDX ) if(trueIDX == idx) good = true;
            if(good) continue;
            auto measuredMuon = getMeasuredReco(idx);
            plotter.getOrMake1D(TString::Format("%s_fs_incl_pt",glbPrefix.Data()),";bkg. segment p_{T}",60,0,30)->Fill(  measuredMuon.pt() );
            if(std::fabs(measuredMuon.eta()) < 2.4) continue;
            plotter.getOrMake1D(TString::Format("%s_fs_pt",glbPrefix.Data()),";bkg. segment p_{T}",60,0,30)->Fill(  measuredMuon.pt() );

            if(measuredMuon.pt() > 2) {
                plotter.getOrMake1D(TString::Format("%s_fs_singleMu5_me0Mugt2_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt5SingleMu).mass() );
                plotter.getOrMake1D(TString::Format("%s_fs_singleMu20_me0Mugt2_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt20SingleMu).mass() );

            }
            if(measuredMuon.pt() > 5) {
                plotter.getOrMake1D(TString::Format("%s_fs_singleMu5_me0Mugt5_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt5SingleMu).mass() );
                plotter.getOrMake1D(TString::Format("%s_fs_singleMu20_me0Mugt5_invmass",glbPrefix.Data()),";single muon + ME0Segment mass",200,0,200)->Fill(  (measuredMuon + pt20SingleMu).mass() );

            }




        }

    }

    void tdrPlot() {

        int nSM = 0;
        int nMOOAcc = 0;
        int nMIAcc = 0;
        for(unsigned int iM = 0; iM < simMuon_pt->size(); ++iM ){
            const float pt  = simMuon_pt->at(iM);
            const float eta = simMuon_eta->at(iM);
            const float absEta = std::fabs(eta);
            if(absEta < 1.5 || absEta > 3.5) nMOOAcc ++;
            else nMIAcc++;;
        }
        if(nMIAcc > 0) return;
        if(nMOOAcc != 2) return;
        TString sname = "p8s384__";
        plotter.getOrMake1D(TString::Format("%s_fake_nEvents",sname.Data()),";nEvents",1,0,2.0)->Fill(1);


        for(unsigned int idx = 0; idx < segment_eta->size(); ++idx){

            float dPhi = std::fabs(segment_dphi->at(idx));
            float eta  = std::fabs(segment_eta->at(idx));


            plotter.getOrMake1D(TString::Format("%s_fake_seg_eta",sname.Data()),";segment |#eta|",120,1.8,3.0)->Fill(eta);
            if(dPhi <0.013)  plotter.getOrMake1D(TString::Format("%s_fake_dPhi_lt0p013_seg_eta",sname.Data()),";segment |#eta|",120,1.8,3.0)->Fill(eta);
            if(dPhi <0.004)  plotter.getOrMake1D(TString::Format("%s_fake_dPhi_lt0p004_seg_eta",sname.Data()),";segment |#eta|",120,1.8,3.0)->Fill(eta);
            if(dPhi <0.002)  plotter.getOrMake1D(TString::Format("%s_fake_dPhi_lt0p002_seg_eta",sname.Data()),";segment |#eta|",120,1.8,3.0)->Fill(eta);


        }

    }

    virtual void runAEvent() {
        tdrPlot();
        fakeStudies();
        zMuMuStudies();
        tau3MuStudies();
    }


    Tau3MuAnalyzer tau3MuTree;
};

#endif

void AnalyzeME0SegmentTree(std::string fileName, std::string prefix, int nStrips, std::string outName){

    //	FindParams a (fileName,"Events");
    //	a.glbPrefix = prefix;
    //	a.minDPhi = 0.5*(20.0*TMath::Pi()/180.0)/float(nStrips);
    //	a.nStrips = nStrips;
    //	a.analyze();
    //	a.solve();
    //	a.write(outName);

    TestParams a (fileName,"Events");
    a.glbPrefix = prefix;
    a.minDPhi = 0.5*(20.0*TMath::Pi()/180.0)/float(nStrips);
    a.nStrips = nStrips;
    a.analyze();
    a.write(outName);
}

