
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../AnalsysisSupport/macros/BaseTupleAnalyzer.h"
#include "../AnalsysisSupport/interface/HistGetter.h"
#include "TVector2.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
using namespace std;
double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}

class Analyzer : public BaseTupleAnalyzer{
public:
  Analyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){

    gen_only_type            = new std::vector<int>  ;
    gen_only_pdgid           = new std::vector<int>  ;
    gen_only_pt              = new std::vector<float>;
    gen_only_eta             = new std::vector<float>;
    gen_only_phi             = new std::vector<float>;
    weight                   =   0;
    is_ambiguous             =   false;
    track_pt                 =   new std::vector<float>;
    track_eta                =   new std::vector<float>;
    track_phi                =   new std::vector<float>;
    track_iso                =   new std::vector<float>;
    track_proj_phi           =   new std::vector<float>;
    track_proj_eta           =   new std::vector<float>;
    track_proj_deltaPhi      =   new std::vector<float>;
    gen_type                 =   new std::vector<int>  ;
    gen_pdgid                =   new std::vector<int>  ;
    gen_pt                   =   new std::vector<float>;
    gen_eta                  =   new std::vector<float>;
    gen_phi                  =   new std::vector<float>;
    segment_phi              =   new std::vector<float>;
    segment_eta              =   new std::vector<float>;
    segment_nHits            =   new std::vector<int>  ;
    segment_nEHits           =   new std::vector<int>  ;
    segment_deltaPhi         =   new std::vector<float>;




    setBranchAddress(     "gen_only_type"    ,&gen_only_type  );
    setBranchAddress(     "gen_only_pdgid"   ,&gen_only_pdgid );
    setBranchAddress(     "gen_only_pt"      ,&gen_only_pt    );
    setBranchAddress(     "gen_only_eta"     ,&gen_only_eta   );
    setBranchAddress(     "gen_only_phi"     ,&gen_only_phi   );



    setBranchAddress("weight"             ,&weight);
    setBranchAddress("is_ambiguous"       ,&is_ambiguous        );
    setBranchAddress("track_pt"           ,&track_pt            );
    setBranchAddress("track_eta"          ,&track_eta           );
    setBranchAddress("track_phi"          ,&track_phi           );
    setBranchAddress("track_iso"          ,&track_iso           );
    setBranchAddress("track_proj_phi"     ,&track_proj_phi      );
    setBranchAddress("track_proj_eta"     ,&track_proj_eta      );
    setBranchAddress("track_proj_deltaPhi",&track_proj_deltaPhi );
    setBranchAddress("gen_type"           ,&gen_type            );
    setBranchAddress("gen_pdgid"          ,&gen_pdgid           );
    setBranchAddress("gen_pt"             ,&gen_pt              );
    setBranchAddress("gen_eta"            ,&gen_eta             );
    setBranchAddress("gen_phi"            ,&gen_phi             );
    setBranchAddress("segment_phi"        ,&segment_phi         );
    setBranchAddress("segment_eta"        ,&segment_eta         );
    setBranchAddress("segment_nHits"      ,&segment_nHits       );
    setBranchAddress("segment_nEHits"     ,&segment_nEHits      );
    setBranchAddress("segment_deltaPhi"   ,&segment_deltaPhi    );

  }



  std::vector<int>  *  gen_only_type   ;
  std::vector<int>  *  gen_only_pdgid  ;
  std::vector<float>*  gen_only_pt     ;
  std::vector<float>*  gen_only_eta    ;
  std::vector<float>*  gen_only_phi    ;

  float weight       ;
  bool is_ambiguous       ;
  std::vector<float>*  track_pt             ;
  std::vector<float>*  track_eta            ;
  std::vector<float>*  track_phi            ;
  std::vector<float>*  track_iso            ;
  std::vector<float>*  track_proj_phi       ;
  std::vector<float>*  track_proj_eta       ;
  std::vector<float>*  track_proj_deltaPhi  ;
  std::vector<int>  *  gen_type             ;
  std::vector<int>  *  gen_pdgid            ;
  std::vector<float>*  gen_pt               ;
  std::vector<float>*  gen_eta              ;
  std::vector<float>*  gen_phi              ;
  std::vector<float>*  segment_phi          ;
  std::vector<float>*  segment_eta          ;
  std::vector<int>  *  segment_nHits        ;
  std::vector<int>  *  segment_nEHits       ;
  std::vector<float>*  segment_deltaPhi     ;



  void basePlots(){
    auto fill = [&](int bin){
      plotter.getOrMake1D("EventCounts","",30,-0.5,29.5)->Fill(bin);
    };

    fill(0);
    vector<float> pts;
    vector<float> etas;

    for(unsigned int iT = 0; iT < gen_only_pt->size(); ++iT){
      if(gen_only_type->at(iT) == 0 ) continue;
      if(TMath::Abs(gen_only_pdgid->at(iT)) != 13 ) continue;
      pts.push_back(gen_only_pt->at(iT));
      etas.push_back(TMath::Abs(gen_only_eta->at(iT)));
    }

    if(pts.size() != 4) return;
    fill(1);

    int largestETA = -1;
    int lowestPT = -1;
    for(unsigned int iP = 0; iP < pts.size(); ++iP){
      if(largestETA < 0 || etas[iP] > etas[largestETA]) largestETA  = iP;
      if(lowestPT < 0 || pts[iP] < pts[lowestPT]) lowestPT  = iP;
    }

    plotter.getOrMake1D("incl_etadist",";|#eta|;a.u.",50,0,5.0)->Fill(etas[largestETA]);
    if(pts[lowestPT] < 3) plotter.getOrMake1D("pt_lt3_etadist",";|#eta|;a.u.",50,0,5.0)->Fill(etas[largestETA]);
    else if(pts[lowestPT] < 5) plotter.getOrMake1D("pt3to5_etadist",";|#eta|;a.u.",50,0,5.0)->Fill(etas[largestETA]);
    else if(pts[lowestPT] < 10) plotter.getOrMake1D("pt5to10_etadist",";|#eta|;a.u.",50,0,5.0)->Fill(etas[largestETA]);
    else if(pts[lowestPT] >= 10) plotter.getOrMake1D("ptgt10_etadist",";|#eta|;a.u.",50,0,5.0)->Fill(etas[largestETA]);

    float etaBins[] = {0,2.4,2.8,3.2};
    float ptBins[]  = {0,3,5,10,15};
    plotter.getOrMake2D("eta_pt_dist",";Maximum muon |#eta|;Minimum muon p_{T} [GeV]",3,etaBins,4,ptBins)->Fill(std::min(etas[largestETA],float(3.1)),std::min(pts[lowestPT],float(14)) );


    auto count = [&] (float minPT, float maxETA) ->int {
      int ct = 0;
      for(unsigned int iP = 0; iP < pts.size(); ++iP){
        if(etas[iP] > maxETA ) continue;
        if(pts[iP] < minPT) continue;
        ct++;
      }
      return ct;
    };

    bool passTwoCentral = count(17,2.4) > 0 && count(8,2.4) > 1;

    auto fillETA = [&] (double eta, int firstBin) -> int{
      if(count(-1,eta) >= 4 ) fill(firstBin);
      if(count(3,eta) >= 4 ) fill(firstBin+1);
      if(count(5,eta) >= 4 ) fill(firstBin+2);
      if(count(8,eta) >= 4 ) fill(firstBin+3);
      if(passTwoCentral && count(-1,eta) >= 4 ) fill(firstBin+4);
      if(passTwoCentral && count(3,eta) >= 4 ) fill(firstBin+5);
      if(passTwoCentral && count(5,eta) >= 4 ) fill(firstBin+6);
      if(passTwoCentral && count(8,eta) >= 4 ) fill(firstBin+7);
      return firstBin+8;
    };

    int firstBin = 2;
    firstBin = fillETA(2.4,firstBin);
    firstBin = fillETA(2.8,firstBin);
    firstBin = fillETA(4.0,firstBin);



    auto fill2 = [&](int bin){
      plotter.getOrMake1D("EventCounts2","",30,-0.5,29.5)->Fill(bin);
    };


    auto count2 = [&] (float minPT,float minETA, float maxETA) ->int {
      int ct = 0;
      for(unsigned int iP = 0; iP < pts.size(); ++iP){
        if(etas[iP] <= minETA ) continue;
        if(etas[iP] > maxETA ) continue;
        if(pts[iP] < minPT) continue;
        ct++;
      }
      return ct;
    };

    fill2(0);
    if(count2(-1,-1, 2.4) >= 4 ) fill2(1);
    if(count2(3,-1, 2.4) >= 4 )  fill2(2);
    if(count2(5,-1, 2.4) >= 4 )  fill2(3);
    if(count2(5,-1, 2.4) +count2(3,2.4, 2.8)  >= 4 )  fill2(4);
    if(count2(5,-1, 2.8)  >= 4 )  fill2(5);
    if(count2(5,-1, 2.4) +count2(10,2.4, 2.8)  >= 4 )  fill2(6);

    if(passTwoCentral){
      fill2(7);
      if(count2(-1,-1, 2.4) >= 4 ) fill2(8);
      if(count2(3,-1, 2.4) >= 4 ) fill2(9);
      if(count2(5,-1, 2.4) >= 4 )  fill2(10);
      if(count2(5,-1, 2.4) +count2(3,2.4, 2.8)  >= 4 )  fill2(11);
      if(count2(5,-1, 2.8)  >= 4 )  fill2(12);
      if(count2(5,-1, 2.4) +count2(10,2.4, 2.8)  >= 4 )  fill2(13);
    }




  }


  void checkTrackISO() {
    if(is_ambiguous) return;
    plotter.getOrMake1D("checkTrackIso_nEvents",";nEvents",1,0,2)->Fill(1);

    int nISOTrks = 0;
    for(unsigned int iT  = 0; iT < track_pt->size(); ++iT){
      if(track_pt->at(iT) < 5 || TMath::Abs(track_eta->at(iT)) > 2.4  ) continue;
      if(track_iso->at(iT) > 0.25) continue;
      nISOTrks++;
    }
    if(nISOTrks >= 3) plotter.getOrMake1D("checkTrackIso_3ISO_nEvents",";nEvents",1,0,2)->Fill(1);

    for(unsigned int iT  = 0; iT < track_pt->size(); ++iT){
      bool isMuon = false;
      bool isBKG = false;
      if(gen_pdgid->at(iT) == 0 ) isBKG = true;
      if(TMath::Abs(gen_pdgid->at(iT)) == 13 && gen_type->at(iT) != 0  ) isMuon = true;

      TString name = isMuon ? "zMuon_" : (isBKG ? "bkg_" : "other_");
      if(track_pt->at(iT) < 3) name += "ptlt3_";
      else if(track_pt->at(iT) < 5) name += "pteq3to5_";
      else if(track_pt->at(iT) < 10) name += "pteq5to10_";
      else if(track_pt->at(iT) < 20) name += "pteq10to20_";
      else name += "ptgeq20_";

      if(TMath::Abs(track_eta->at(iT)) < 1.4 ) name += "etalt1p4_";
      else if(TMath::Abs(track_eta->at(iT)) < 2.4 ) name += "etaeq1p4to2p4_";
      else if(TMath::Abs(track_eta->at(iT)) < 2.8 ) name += "etaeq2p4to2p8_";
      else name += "etageq2p8_";

      plotter.getOrMake1D(TString::Format("%siso",name.Data()),";relative track iso",100,0,5.0)->Fill(track_iso->at(iT));

      double pTs[] = {0,3,5,10,20,30};
      double etas[] = {0,1.4,2.4,2.8};
      TString name2D = isMuon ? "zMuon_" : (isBKG ? "bkg_" : "other_");
      if(track_iso->at(iT) < 0.25)
      plotter.getOrMake2D(TString::Format("%spass_2D_pass",name2D.Data()),";p_{T} [GeV]; |#eta|",5,pTs,3,etas)->Fill(track_pt->at(iT) > 30 ? 28 : track_pt->at(iT) , TMath::Abs(track_eta->at(iT))  );
      plotter.getOrMake2D(TString::Format("%sincl_2D_pass",name2D.Data()),";p_{T} [GeV]; |#eta|",5,pTs,3,etas)->Fill(track_pt->at(iT) > 30 ? 28 : track_pt->at(iT) , TMath::Abs(track_eta->at(iT))  );


      if(isBKG){

        auto fill = [&] (TString name) {
          plotter.getOrMake1D(TString::Format("bkg_%sincl_eta",name.Data()),";|#eta|",100,0,5.0)->Fill(TMath::Abs(track_eta->at(iT)));
          if(track_iso->at(iT) < 0.25)
          plotter.getOrMake1D(TString::Format("bkg_%siso_eta",name.Data()),";|#eta|",100,0,5.0)->Fill(TMath::Abs(track_eta->at(iT)));
        };

        if(track_pt->at(iT) >= 0.5) fill("ptgeq0p5_");
        if(track_pt->at(iT) >= 1) fill("ptgeq1_");
        if(track_pt->at(iT) >= 3) fill("ptgeq3_");
        if(track_pt->at(iT) >= 5) fill("ptgeq5_");
        if(track_pt->at(iT) >= 10) fill("ptgeq10_");

        if(nISOTrks >= 3 && track_pt->at(iT) >= 0.5) fill("ptgeq0p5_3ISO_");
        if(nISOTrks >= 3 && track_pt->at(iT) >= 1) fill("ptgeq1_3ISO_");
        if(nISOTrks >= 3 && track_pt->at(iT) >= 3) fill("ptgeq3_3ISO_");
        if(nISOTrks >= 3 && track_pt->at(iT) >= 5) fill("ptgeq5_3ISO_");
        if(nISOTrks >= 3 && track_pt->at(iT) >= 10) fill("ptgeq10_3ISO_");






      }

    }



  }







  virtual void runAEvent() {
    basePlots();
    checkTrackISO();



  }

  void write(TString fileName){ plotter.write(fileName);}
  HistGetter plotter;
};

#endif

void AnalyzeME0FakeRateTree(std::string fileName ="testout.root",std::string outFileName = "plots.root"){
  Analyzer a(fileName,"Events");
  a.analyze();
  a.write(outFileName);
}
