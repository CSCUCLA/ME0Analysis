
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "BaseTupleAnalyzer.h"
#include "HistGetter.h"
#include "TVector2.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
using namespace std;
double deltaPhi(const float phi1, const float phi2) {return TVector2::Phi_mpi_pi(phi1 - phi2);}

class Analyzer : public BaseTupleAnalyzer{
public:
  Analyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){




    gen_lay_eta         = new std::vector<float>;
    gen_lay_phi         = new std::vector<float>;
    gen_lay_lay         = new std::vector<size8>;
    rh_lay_eta          = new std::vector<float>;
    rh_lay_phi          = new std::vector<float>;
    rh_lay_eta_e        = new std::vector<float>;
    rh_lay_phi_e        = new std::vector<float>;
    rh_lay_lay          = new std::vector<size8>;
    track_isMatched     = new std::vector<size8>;
    track_pt            = new std::vector<float>;
    track_eta           = new std::vector<float>;
    track_phi           = new std::vector<float>;
    track_center_eta    = new std::vector<float>;
    track_center_phi    = new std::vector<float>;
    track_center_eta_e  = new std::vector<float>;
    track_center_phi_e  = new std::vector<float>;
    track_delta_eta     = new std::vector<float>;
    track_delta_phi     = new std::vector<float>;
    track_delta_eta_e   = new std::vector<float>;
    track_delta_phi_e   = new std::vector<float>;


    setBranchAddress("gen_pdgid"            ,&gen_pdgid         );
    setBranchAddress("gen_pt"               ,&gen_pt            );
    setBranchAddress("gen_eta"              ,&gen_eta           );
    setBranchAddress("gen_phi"              ,&gen_phi           );
    setBranchAddress("gen_center_eta"       ,&gen_center_eta    );
    setBranchAddress("gen_center_phi"       ,&gen_center_phi    );
    setBranchAddress("gen_delta_eta"        ,&gen_delta_eta     );
    setBranchAddress("gen_delta_phi"        ,&gen_delta_phi     );
    setBranchAddress("gen_lay_eta"          ,&gen_lay_eta       );
    setBranchAddress("gen_lay_phi"          ,&gen_lay_phi       );
    setBranchAddress("gen_lay_lay"          ,&gen_lay_lay       );
    setBranchAddress("rh_center_eta"        ,&rh_center_eta     );
    setBranchAddress("rh_center_phi"        ,&rh_center_phi     );
    setBranchAddress("rh_center_eta_e"      ,&rh_center_eta_e   );
    setBranchAddress("rh_center_phi_e"      ,&rh_center_phi_e   );
    setBranchAddress("rh_delta_eta"         ,&rh_delta_eta      );
    setBranchAddress("rh_delta_phi"         ,&rh_delta_phi      );
    setBranchAddress("rh_delta_eta_e"       ,&rh_delta_eta_e    );
    setBranchAddress("rh_delta_phi_e"       ,&rh_delta_phi_e    );
    setBranchAddress("rh_lay_eta"           ,&rh_lay_eta        );
    setBranchAddress("rh_lay_phi"           ,&rh_lay_phi        );
    setBranchAddress("rh_lay_eta_e"         ,&rh_lay_eta_e      );
    setBranchAddress("rh_lay_phi_e"         ,&rh_lay_phi_e      );
    setBranchAddress("rh_lay_lay"           ,&rh_lay_lay        );
    setBranchAddress("track_isMatched"      ,&track_isMatched   );
    setBranchAddress("track_pt"             ,&track_pt          );
    setBranchAddress("track_eta"            ,&track_eta         );
    setBranchAddress("track_phi"            ,&track_phi         );
    setBranchAddress("track_center_eta"     ,&track_center_eta  );
    setBranchAddress("track_center_phi"     ,&track_center_phi  );
    setBranchAddress("track_center_eta_e"   ,&track_center_eta_e);
    setBranchAddress("track_center_phi_e"   ,&track_center_phi_e);
    setBranchAddress("track_delta_eta"      ,&track_delta_eta   );
    setBranchAddress("track_delta_phi"      ,&track_delta_phi   );
    setBranchAddress("track_delta_eta_e"    ,&track_delta_eta_e );
    setBranchAddress("track_delta_phi_e"    ,&track_delta_phi_e );


  }


  int                 gen_pdgid           ;
  float               gen_pt              ;
  float               gen_eta             ;
  float               gen_phi             ;
  float               gen_center_eta      ;
  float               gen_center_phi      ;
  float               gen_delta_eta       ;
  float               gen_delta_phi       ;
  std::vector<float>* gen_lay_eta         ;
  std::vector<float>* gen_lay_phi         ;
  std::vector<size8>* gen_lay_lay         ;
  float               rh_center_eta       ;
  float               rh_center_phi       ;
  float               rh_center_eta_e     ;
  float               rh_center_phi_e     ;
  float               rh_delta_eta        ;
  float               rh_delta_phi        ;
  float               rh_delta_eta_e      ;
  float               rh_delta_phi_e      ;
  std::vector<float>* rh_lay_eta          ;
  std::vector<float>* rh_lay_phi          ;
  std::vector<float>* rh_lay_eta_e        ;
  std::vector<float>* rh_lay_phi_e        ;
  std::vector<size8>* rh_lay_lay          ;
  std::vector<size8>* track_isMatched     ;
  std::vector<float>* track_pt            ;
  std::vector<float>* track_eta           ;
  std::vector<float>* track_phi           ;
  std::vector<float>* track_center_eta    ;
  std::vector<float>* track_center_phi    ;
  std::vector<float>* track_center_eta_e  ;
  std::vector<float>* track_center_phi_e  ;
  std::vector<float>* track_delta_eta     ;
  std::vector<float>* track_delta_phi     ;
  std::vector<float>* track_delta_eta_e   ;
  std::vector<float>* track_delta_phi_e   ;


  void trackGenPlots(){
    if(track_pt->size() == 0 || !track_isMatched->at(0)) return;

    //plot gen info plots
    TString prefix = "";
    if(gen_pt >= 1 && gen_pt < 3){
      prefix = "gen_pt1to3";
    } else if(gen_pt >= 3 && gen_pt < 5){
      prefix = "gen_pt3to5";
    } else if(gen_pt >= 5 && gen_pt < 10){
      prefix = "gen_pt5to10";
    } else if(gen_pt >= 10 && gen_pt < 20){
      prefix = "gen_pt10to20";
    } else if (gen_pt >= 20){
      prefix = "gen_ptgeq20";
    } else
      return;

    plotter.getOrMake1D(TString::Format("%s_tg_dpt",prefix.Data()),";(track p_{T}-gen p_{T})/gen p_{T}",25,-.2,0.2)->Fill((track_pt->front()-gen_pt)/gen_pt);
    plotter.getOrMake1D(TString::Format("%s_tg_curv",prefix.Data()),";(track p_{T}-gen p_{T})/gen p_{T}^{2} [GeV^{-1}]",50,-0.05,0.05)->Fill((track_pt->front()-gen_pt)/(gen_pt*gen_pt));
    plotter.getOrMake1D(TString::Format("%s_tg_dphi",prefix.Data()),";track #phi - gen #phi (in ME0)",50,-.05,.05)->Fill(track_center_phi->front() - gen_center_phi);
    plotter.getOrMake1D(TString::Format("%s_tg_deta",prefix.Data()),";track #eta - gen #eta (in ME0)",50,-.05,.05)->Fill(track_center_eta->front() - gen_center_eta);
    plotter.getOrMake1D(TString::Format("%s_tg_ddphi",prefix.Data()),";track #Delta#phi - gen #Delta#phi (in ME0)",50,-.01,.01)->Fill(track_delta_phi->front() - gen_delta_phi);
    plotter.getOrMake1D(TString::Format("%s_tg_ddeta",prefix.Data()),";track #Delta#eta - gen #Delta#eta (in ME0)",50,-.01,.01)->Fill(track_delta_eta->front() - gen_delta_eta);


    plotter.getOrMake1D(TString::Format("%s_tg_abs_dphi",prefix.Data()) ,";|track #phi - gen #phi| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->front(),gen_center_phi)));
    plotter.getOrMake1D(TString::Format("%s_tg_abs_deta",prefix.Data()) ,";|track #eta - gen #eta| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(track_center_eta->front() - gen_center_eta));
    plotter.getOrMake1D(TString::Format("%s_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_phi->front() - gen_delta_phi  ));
    plotter.getOrMake1D(TString::Format("%s_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_eta->front() - gen_delta_eta  ));

    plotter.getOrMake2D(TString::Format("%s_tg_abs_dphi_deta",prefix.Data()) ,";|track #phi - gen #phi| (in ME0);|track #eta - gen #eta| (in ME0)" ,2000,-.01,.4,2000,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->front(),gen_center_phi)),TMath::Abs(track_center_eta->front() - gen_center_eta));


    if(TMath::Abs(deltaPhi(track_center_phi->front(),gen_center_phi)) < 0.04 && TMath::Abs(track_center_eta->front() - gen_center_eta) < 0.04 ){
      plotter.getOrMake1D(TString::Format("%s_tight_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_phi->front() - gen_delta_phi  ));
      plotter.getOrMake1D(TString::Format("%s_tight_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_eta->front() - gen_delta_eta  ));
      plotter.getOrMake1D(TString::Format("%s_tight_fine_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_phi->front() - gen_delta_phi  ));
      plotter.getOrMake1D(TString::Format("%s_tight_fine_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_eta->front() - gen_delta_eta  ));
    } else {
      plotter.getOrMake1D(TString::Format("%s_tight_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(10);
      plotter.getOrMake1D(TString::Format("%s_tight_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(10);
      plotter.getOrMake1D(TString::Format("%s_tight_fine_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(3.19);
      plotter.getOrMake1D(TString::Format("%s_tight_fine_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(3.19);
    }



    for(unsigned int iT = 0; iT < track_pt->size(); ++iT){
      if(track_isMatched->at(iT) ) continue;
      if(TMath::Abs(track_center_eta->at(iT) - gen_center_eta) >= 0.35 || TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)) >= 0.35    ) continue;
      plotter.getOrMake1D(TString::Format("%s_fake_tg_abs_dphi",prefix.Data()) ,";|track #phi - gen #phi| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)));
      plotter.getOrMake1D(TString::Format("%s_fake_tg_abs_deta",prefix.Data()) ,";|track #eta - gen #eta| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(track_center_eta->at(iT) - gen_center_eta));
      plotter.getOrMake1D(TString::Format("%s_fake_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_phi->at(iT)  - gen_delta_phi  ));
      plotter.getOrMake1D(TString::Format("%s_fake_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_eta->at(iT)  - gen_delta_eta  ));

      plotter.getOrMake2D(TString::Format("%s_fake_tg_abs_dphi_deta",prefix.Data()) ,";|track #phi - gen #phi| (in ME0);|track #eta - gen #eta| (in ME0)" ,2000,-.01,.4,2000,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)),TMath::Abs(track_center_eta->at(iT) - gen_center_eta));


      if(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)) < 0.04 && TMath::Abs(track_center_eta->at(iT) - gen_center_eta) < 0.04 ){
        plotter.getOrMake1D(TString::Format("%s_tight_fake_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_phi->at(iT) - gen_delta_phi  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fake_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_eta->at(iT) - gen_delta_eta  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_phi->at(iT) - gen_delta_phi  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_eta->at(iT) - gen_delta_eta  ));
      } else {
        plotter.getOrMake1D(TString::Format("%s_tight_fake_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(10);
        plotter.getOrMake1D(TString::Format("%s_tight_fake_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(10);

        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(3.19);
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(3.19);
      }

    }


    for(unsigned int iT = 0; iT < track_pt->size(); ++iT){
      if(track_isMatched->at(iT) ) continue;
      if(track_pt->at(iT) < 3 ) continue;
      if(TMath::Abs(track_center_eta->at(iT) - gen_center_eta) >= 0.35 || TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)) >= 0.35    ) continue;
      plotter.getOrMake1D(TString::Format("%s_fake_3GeV_tg_abs_dphi",prefix.Data()) ,";|track #phi - gen #phi| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)));
      plotter.getOrMake1D(TString::Format("%s_fake_3GeV_tg_abs_deta",prefix.Data()) ,";|track #eta - gen #eta| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(track_center_eta->at(iT) - gen_center_eta));
      plotter.getOrMake1D(TString::Format("%s_fake_3GeV_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_phi->at(iT)  - gen_delta_phi  ));
      plotter.getOrMake1D(TString::Format("%s_fake_3GeV_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(track_delta_eta->at(iT)  - gen_delta_eta  ));

      plotter.getOrMake2D(TString::Format("%s_fake_3GeV_tg_abs_dphi_deta",prefix.Data()) ,";|track #phi - gen #phi| (in ME0);|track #eta - gen #eta| (in ME0)" ,2000,-.01,.4,2000,-.01,.4)->Fill(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)),TMath::Abs(track_center_eta->at(iT) - gen_center_eta));


      if(TMath::Abs(deltaPhi(track_center_phi->at(iT),gen_center_phi)) < 0.04 && TMath::Abs(track_center_eta->at(iT) - gen_center_eta) < 0.04 ){
        plotter.getOrMake1D(TString::Format("%s_tight_fake_3GeV_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_phi->at(iT) - gen_delta_phi  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fake_3GeV_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(TMath::Abs(track_delta_eta->at(iT) - gen_delta_eta  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_3GeV_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_phi->at(iT) - gen_delta_phi  ));
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_3GeV_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(TMath::Abs(track_delta_eta->at(iT) - gen_delta_eta  ));
      } else {
        plotter.getOrMake1D(TString::Format("%s_tight_fake_3GeV_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",250,-.01,.1)->Fill(10);
        plotter.getOrMake1D(TString::Format("%s_tight_fake_3GeV_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",250,-.01,.1)->Fill(10);

        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_3GeV_tg_abs_ddphi",prefix.Data()),";|track #Delta#phi - gen #Delta#phi| (in ME0)",6420,-.01,3.2)->Fill(3.19);
        plotter.getOrMake1D(TString::Format("%s_tight_fine_fake_3GeV_tg_abs_ddeta",prefix.Data()),";|track #Delta#eta - gen #Delta#eta| (in ME0)",6420,-.01,3.2)->Fill(3.19);
      }

    }


  }

  void hitGenPlots(){
     if(gen_pt < 1) return;
     if(rh_lay_eta->size() < 6) return;
     if(rh_center_phi < -70) return;


     for(unsigned int iG = 0; iG < gen_lay_eta->size(); ++iG){
       for(unsigned int iH = 0; iH < rh_lay_eta->size(); ++iH){
         if(gen_lay_lay->at(iG) != rh_lay_lay->at(iG) ) continue;
         plotter.getOrMake1D("rg_rechit_deta",";rechit #eta - simhit #eta",50,-.3,.3)->Fill(rh_lay_eta->at(iH) - gen_lay_eta->at(iG));
         plotter.getOrMake1D("rg_rechit_dphi",";rechit #phi - simhit #phi",50,-.02,.02)->Fill(rh_lay_phi->at(iH) - gen_lay_phi->at(iG));
         break;
       }
     }

     //plot gen info plots
     TString prefix = "";
     if(gen_pt >= 1 && gen_pt < 3){
       prefix = "gen_pt1to3";
     } else if(gen_pt >= 3 && gen_pt < 5){
       prefix = "gen_pt3to5";
     } else if(gen_pt >= 5 && gen_pt < 10){
       prefix = "gen_pt5to10";
     } else{
       prefix = "gen_ptgeq10";
     }


     plotter.getOrMake1D(TString::Format("%s_rg_dphi",prefix.Data()),";segment #phi - gen #phi (in ME0)",50,-.3,.3)->Fill(deltaPhi(rh_center_phi,gen_center_phi));
     plotter.getOrMake1D(TString::Format("%s_rg_deta",prefix.Data()),";segment #eta - gen #eta (in ME0)",50,-.3,.3)->Fill(rh_center_eta - gen_center_eta);
     plotter.getOrMake1D(TString::Format("%s_rg_ddphi",prefix.Data()),";segment #Delta#phi - gen #Delta#phi (in ME0)",50,-.3,.3)->Fill(rh_delta_phi - gen_delta_phi);
     plotter.getOrMake1D(TString::Format("%s_rg_ddeta",prefix.Data()),";segment #Delta#eta - gen #Delta#eta (in ME0)",50,-.3,.3)->Fill(rh_delta_eta - gen_delta_eta);


     plotter.getOrMake1D(TString::Format("%s_rg_abs_dphi",prefix.Data()) ,";|segment #phi - gen #phi| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(deltaPhi(rh_center_phi,gen_center_phi)));
     plotter.getOrMake1D(TString::Format("%s_rg_abs_deta",prefix.Data()) ,";|segment #eta - gen #eta| (in ME0)"            ,70,-.01,.4)->Fill(TMath::Abs(rh_center_eta - gen_center_eta));
     plotter.getOrMake1D(TString::Format("%s_rg_abs_ddphi",prefix.Data()),";| #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(rh_delta_phi - gen_delta_phi  ));
     plotter.getOrMake1D(TString::Format("%s_rg_abs_ddeta",prefix.Data()),";| #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.4)->Fill(TMath::Abs(rh_delta_eta - gen_delta_eta  ));

     if(TMath::Abs(deltaPhi(rh_center_phi,gen_center_phi)) < 0.04 && TMath::Abs(rh_center_eta - gen_center_eta) < 0.04 ){
       plotter.getOrMake1D(TString::Format("%s_tight_rg_abs_ddphi",prefix.Data()),";|segment #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.1)->Fill(TMath::Abs(rh_delta_phi - gen_delta_phi  ));
       plotter.getOrMake1D(TString::Format("%s_tight_rg_abs_ddeta",prefix.Data()),";|segment #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.1)->Fill(TMath::Abs(rh_delta_eta - gen_delta_eta  ));
     } else {
       plotter.getOrMake1D(TString::Format("%s_tight_rg_abs_ddphi",prefix.Data()),";|segment #Delta#phi - gen #Delta#phi| (in ME0)",70,-.01,.1)->Fill(10);
       plotter.getOrMake1D(TString::Format("%s_tight_rg_abs_ddeta",prefix.Data()),";|segment #Delta#eta - gen #Delta#eta| (in ME0)",70,-.01,.1)->Fill(10);
     }

   }

  virtual void runAEvent() {
    if(gen_lay_eta->size() != 6) return;

    double norm = (gen_pdgid > 0 ? 1 : -1);

    //plot gen info plots
    TString prefix = "gen_";
    if(gen_pt >= 0.5 && gen_pt < 1){
      prefix += "pt0p5to1";
    }
    if(gen_pt >= 1 && gen_pt < 3){
      prefix += "pt1to3";
    }
    if(gen_pt >= 3 && gen_pt < 5){
      prefix += "pt3to5";
    }
    if(gen_pt >= 5 && gen_pt < 10){
      prefix += "pt5to10";
    }
    if(gen_pt >= 10 && gen_pt < 20){
      prefix += "pt10to20";
    }
    if(gen_pt >= 20){
      prefix += "ptgeq20";
    }
    plotter.getOrMake1D(TString::Format("%s_dphi",prefix.Data()),";|#phi_{lay 6}-#phi_{lay 1}| (in ME0)",50,0,.05)->Fill(TMath::Abs(gen_delta_phi));
    prefix += "_";
    if(TMath::Abs(gen_eta) >= 2 && TMath::Abs(gen_eta)  < 2.1) prefix += "eta2to2p1";
    if(TMath::Abs(gen_eta) >= 2.1 && TMath::Abs(gen_eta)  < 2.2) prefix += "eta2p1to2p2";
    if(TMath::Abs(gen_eta) >= 2.6 && TMath::Abs(gen_eta)  < 2.7) prefix += "eta2p6to2p7";
//    else if(TMath::Abs(gen_eta) >= 2.25 && TMath::Abs(gen_eta)  < 2.5) prefix += "eta2p25to2p5";
//    else if(TMath::Abs(gen_eta) >= 2.5 && TMath::Abs(gen_eta)  < 2.8) prefix += "eta2p25to2p8";
    plotter.getOrMake1D(TString::Format("%s_dphi",prefix.Data()),";|#phi_{lay 6}-#phi_{lay 1}| (in ME0)",50,0,.05)->Fill(TMath::Abs(gen_delta_phi));





    if(TMath::Abs(gen_delta_phi) < 0.001){
      plotter.getOrMake1D("gen_absdphilt0p001_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    } else if (TMath::Abs(gen_delta_phi) < 0.002)
      plotter.getOrMake1D("gen_absdphieq0p001to0p002_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    else if (TMath::Abs(gen_delta_phi) < 0.005)
      plotter.getOrMake1D("gen_absdphieq0p002to0p005_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    else if (TMath::Abs(gen_delta_phi) < 0.008)
      plotter.getOrMake1D("gen_absdphieq0p005to0p008_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    else if (TMath::Abs(gen_delta_phi) < 0.015)
      plotter.getOrMake1D("gen_absdphieq0p008to0p015_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    else if ( TMath::Abs(gen_delta_phi) >= 0.015)
      plotter.getOrMake1D("gen_absdphigeq0p015_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);
    plotter.getOrMake1D("gen_pt",";p_{T} [GeV]",59,0.5,30)->Fill(gen_pt);


    plotter.getOrMake2D("gen_pt_dphi",";p_{T} [GeV];#phi_{lay 6}-#phi_{lay 1} (in ME0)",59,0.5,30, 100,-0.01,0.1)->Fill(gen_pt, norm*gen_delta_phi);


    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> genP (gen_pt,gen_eta,gen_phi,.105658);
    plotter.getOrMake2D("gen_p_gen_pt",";p_{T} [GeV];|p| [GeV]",59,0.5,30, 100,0,200)->Fill(gen_pt, genP.P());



    trackGenPlots();
    hitGenPlots();




//    if(gen_eta < 0 && gen_pdgid < 0) plotter.getOrMake1D("gen_neg_negz",";|#phi|",100,-0.05,0.05)->Fill(norm*gen_delta_phi);
//    if(gen_eta < 0 && gen_pdgid > 0) plotter.getOrMake1D("gen_pos_negz",";|#phi|",100,-0.05,0.05)->Fill(norm*gen_delta_phi);
//    if(gen_eta > 0 && gen_pdgid < 0) plotter.getOrMake1D("gen_neg_posz",";|#phi|",100,-0.05,0.05)->Fill(norm*gen_delta_phi);
//    if(gen_eta > 0 && gen_pdgid > 0) plotter.getOrMake1D("gen_pos_posz",";|#phi|",100,-0.05,0.05)->Fill(norm*gen_delta_phi);




  }

  void write(TString fileName){ plotter.write(fileName);}
  HistGetter plotter;
};

#endif

void AnalyzeME0Tree(std::string fileName ="out.root",std::string treeName = "Events",std::string outFileName = "plots.root"){
  Analyzer a(fileName,treeName);
  a.analyze();
  a.write(outFileName);
}
