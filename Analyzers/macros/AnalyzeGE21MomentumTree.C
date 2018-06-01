

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/TreeInterface/interface/BaseTupleAnalyzer.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/HistGetter.h"
#include "/Users/nmccoll/Dropbox/Work/GitRepositories/AnalysisSupport/Utilities/interface/Types.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
using namespace std;

using ASTypes::int8;
using ASTypes::size;

class Analyzer : public BaseTupleAnalyzer{
public:
    TString glbPrefix = "";
    unsigned int nP = 0;
    unsigned int nS = 0;
    double alphaPhi2SC = 0;
    double maxPhi2SC=0;
    double alphaDPhi2SC = 0;
    double maxDPhi2SC=0;
    double phiC=0;
    double dPhiC=0;
    double etaC= 0;
    bool isPureBKG = false;

    Analyzer(std::string fileName, std::string treeName, std::string sampname) : BaseTupleAnalyzer(fileName,treeName),sampname(sampname){


        setBranchAddress("nM_all"                 ,&nM_all                 );
        setBranchAddress("nM_good"                ,&nM_good                );
        setBranchAddress("genMuon_pt"             ,&genMuon_pt             );
        setBranchAddress("genMuon_eta"            ,&genMuon_eta            );
        setBranchAddress("genMuon_q"              ,&genMuon_q              );
        setBranchAddress("genReco_dr"             ,&genReco_dr             );
        setBranchAddress("global_pt"              ,&global_pt              );
        setBranchAddress("global_eta"             ,&global_eta             );
        setBranchAddress("global_q"               ,&global_q               );
        setBranchAddress("tracker_pt"             ,&tracker_pt             );
        setBranchAddress("tracker_eta"            ,&tracker_eta            );
        setBranchAddress("tracker_q"              ,&tracker_q              );
        setBranchAddress("standalone_pt"          ,&standalone_pt          );
        setBranchAddress("standalone_eta"         ,&standalone_eta         );
        setBranchAddress("standalone_q"           ,&standalone_q           );
        setBranchAddress("picky_pt"               ,&picky_pt               );
        setBranchAddress("picky_eta"              ,&picky_eta              );
        setBranchAddress("picky_q"                ,&picky_q                );
        setBranchAddress("picky_def_pt"           ,&picky_def_pt           );
        setBranchAddress("picky_def_eta"          ,&picky_def_eta          );
        setBranchAddress("picky_def_q"            ,&picky_def_q            );
        setBranchAddress("picky_oneGE21_pt"       ,&picky_oneGE21_pt       );
        setBranchAddress("picky_oneGE21_eta"      ,&picky_oneGE21_eta      );
        setBranchAddress("picky_oneGE21_q"        ,&picky_oneGE21_q        );
        setBranchAddress("picky_noGE21_pt"        ,&picky_noGE21_pt        );
        setBranchAddress("picky_noGE21_eta"       ,&picky_noGE21_eta       );
        setBranchAddress("picky_noGE21_q"         ,&picky_noGE21_q         );
        setBranchAddress("picky_noME21_pt"        ,&picky_noME21_pt        );
        setBranchAddress("picky_noME21_eta"       ,&picky_noME21_eta       );
        setBranchAddress("picky_noME21_q"         ,&picky_noME21_q         );
        setBranchAddress("picky_noGE11ME0_pt"     ,&picky_noGE11ME0_pt     );
        setBranchAddress("picky_noGE11ME0_eta"    ,&picky_noGE11ME0_eta    );
        setBranchAddress("picky_noGE11ME0_q"      ,&picky_noGE11ME0_q      );
        setBranchAddress("picky_noGE11ME0ME21_pt" ,&picky_noGE11ME0ME21_pt );
        setBranchAddress("picky_noGE11ME0ME21_eta",&picky_noGE11ME0ME21_eta);
        setBranchAddress("picky_noGE11ME0ME21_q"  ,&picky_noGE11ME0ME21_q  );

        setBranchAddress("global_v_d"   ,&global_v_d     );
        setBranchAddress("global_v_s"   ,&global_v_s     );
        setBranchAddress("global_v_r"   ,&global_v_r     );
        setBranchAddress("global_v_phie",&global_v_phie  );

        setBranchAddress("picky_v_d"   ,&picky_v_d     );
        setBranchAddress("picky_v_s"   ,&picky_v_s     );
        setBranchAddress("picky_v_r"   ,&picky_v_r     );
        setBranchAddress("picky_v_phie",&picky_v_phie  );



    }

    std:: string sampname;


    size    nM_all                        =0;
    size    nM_good                       =0;
    float   genMuon_pt                    =0;
    float   genMuon_eta                   =0;
    int8    genMuon_q                     =0;
    float   genReco_dr                    =0;
    float   global_pt                     =0;
    float   global_eta                    =0;
    int8    global_q                      =0;
    float   tracker_pt                    =0;
    float   tracker_eta                   =0;
    int8    tracker_q                     =0;
    float   standalone_pt                 =0;
    float   standalone_eta                =0;
    int8    standalone_q                  =0;
    float   picky_pt                      =0;
    float   picky_eta                     =0;
    int8    picky_q                       =0;
    float   picky_oneGE21_pt              =0;
    float   picky_oneGE21_eta             =0;
    int8    picky_oneGE21_q               =0;
    float   picky_def_pt                  =0;
    float   picky_def_eta                 =0;
    int8    picky_def_q                   =0;
    float   picky_noGE21_pt               =0;
    float   picky_noGE21_eta              =0;
    int8    picky_noGE21_q                =0;
    float   picky_noME21_pt               =0;
    float   picky_noME21_eta              =0;
    int8    picky_noME21_q                =0;
    float   picky_noGE11ME0_pt            =0;
    float   picky_noGE11ME0_eta           =0;
    int8    picky_noGE11ME0_q             =0;
    float   picky_noGE11ME0ME21_pt        =0;
    float   picky_noGE11ME0ME21_eta       =0;
    int8    picky_noGE11ME0ME21_q         =0;

    std::vector<int8>          * global_v_d               = new std::vector<int8>;
    std::vector<int8>          * global_v_s               = new std::vector<int8>;
    std::vector<int8>          * global_v_r               = new std::vector<int8>;
    std::vector<float>         * global_v_phie            = new std::vector<float>;

    std::vector<int8>          * picky_v_d               = new std::vector<int8>;
    std::vector<int8>          * picky_v_s               = new std::vector<int8>;
    std::vector<int8>          * picky_v_r               = new std::vector<int8>;
    std::vector<float>         * picky_v_phie            = new std::vector<float>;


    virtual void runAEvent() {
        if(genReco_dr <0 || genReco_dr> 0.1) return;

        std::string prefix = sampname;
        plotter.getOrMake1DPre(prefix.c_str(),"nm_all",";n reco muons",10,-0.5,9.5)->Fill(nM_all);
        plotter.getOrMake1DPre(prefix.c_str(),"nm_good",";n reco muons",10,-0.5,9.5)->Fill(nM_good);

        auto makeResPlots = [&](std::string prefix, float q, float pt){
            if(pt <= 0) return;
            float genQoPT = genMuon_q/genMuon_pt;
            float recoQoPT = q/pt;
            plotter.getOrMake1DPre(prefix.c_str(),"q_o_pt",";q/#it{p}_{T} [GeV^{-1}]",1000,-0.003,0.003)->Fill(recoQoPT);
            plotter.getOrMake1DPre(prefix.c_str(),"rel_q_o_pt",";q/#it{p}_{T} (reco.)  - q/#it{p}_{T} (gen)  ",1000,-0.001,0.001)->Fill(recoQoPT-genQoPT);
        };

        auto makePlot = [&](std::string prefix, float q, float pt){
            const float eta = std::fabs(genMuon_eta);
            if(eta > 1.6 && eta < 2.4) makeResPlots(prefix+"_eta_1p6_2p4",q,pt);
            if(eta > 1.6 && eta < 1.8) makeResPlots(prefix+"_eta_1p6_1p8",q,pt);
            if(eta > 1.8 && eta < 2.0) makeResPlots(prefix+"_eta_1p8_2p0",q,pt);
            if(eta > 2.0 && eta < 2.2) makeResPlots(prefix+"_eta_2p0_2p2",q,pt);
            if(eta > 2.2 && eta < 2.4) makeResPlots(prefix+"_eta_2p2_2p4",q,pt);

        };

        makePlot(prefix + "_global",global_q,global_pt);
        makePlot(prefix + "_tracker",tracker_q,tracker_pt);
        makePlot(prefix + "_standalone",standalone_q,standalone_pt);
        makePlot(prefix + "_picky",picky_q,picky_pt);
        makePlot(prefix + "_noGE21_picky",picky_noGE21_q,picky_noGE21_pt);
        makePlot(prefix + "_noME21_picky",picky_noME21_q,picky_noME21_pt);
        makePlot(prefix + "_noGE11ME0_picky",picky_noGE11ME0_q,picky_noGE11ME0_pt);
        makePlot(prefix + "_noGE11ME0ME21_picky",picky_noGE11ME0ME21_q,picky_noGE11ME0ME21_pt);
        makePlot(prefix + "_oneGE21_picky",picky_oneGE21_q,picky_oneGE21_pt);


        auto fillInfo =[&](const std::string& prefix, const std::vector<int8>& v_d,const std::vector<int8>& v_s, const std::vector<int8>& v_r, const std::vector<float>& v_phie){
            bool filledME0 = false;
            bool filledME1 = false;
            bool filledME2 = false;
            bool filledME3 = false;
            bool filledME4 = false;
            bool filledGE1 = false;
            bool filledGE2 = false;

            double avgOtherError    =  0;
            int nInAverage = 0;
            double minOtherError    = -9;
            double minGE21Error     = -9;

            for(unsigned int irh = 0; irh < v_d.size(); ++irh){
                if(v_d[irh] == DT || v_d[irh] ==RPC) continue;
                bool isOther = false;
                if(v_d[irh] == GEM){
                    if(v_s[irh] == 1 ){
                        filledGE1 = true;
                        isOther = true;
                    } else if(v_s[irh] == 2 ){
                        filledGE2 = true;
                    }
                } else if(v_d[irh] == ME0){
                    filledME0 = true;
                    isOther = true;
                } else if(v_d[irh] == CSC){
                    isOther = true;
                    if(v_s[irh] == 1 ) filledME1 = true;
                    if(v_s[irh] == 2 ) filledME2 = true;
                    if(v_s[irh] == 3 ) filledME3 = true;
                    if(v_s[irh] == 4 ) filledME4 = true;
                }

                if(isOther){
                    minOtherError =  (minOtherError < 0 ||v_phie[irh] < minOtherError)  ?  v_phie[irh] : minOtherError;
                    nInAverage++;
                    avgOtherError +=  sqrt(v_phie[irh]);
                } else {
                    minGE21Error =  (minGE21Error < 0 ||v_phie[irh] < minGE21Error)  ?  v_phie[irh] : minGE21Error;
                }

            }
            minGE21Error = minGE21Error < 0 ? -1 : std::sqrt(minGE21Error);
            avgOtherError = nInAverage == 0 ? -1 :   avgOtherError/float(nInAverage);
            minOtherError = minOtherError < 0 ? -1 : std::sqrt(minOtherError);

            plotter.getOrMake1DPre(prefix.c_str(),"ge21_phierror",";#sigma_{#phi}",1000,0.00005,0.01)->Fill(minGE21Error < 0 ? -1.0 : minGE21Error);
            plotter.getOrMake1DPre(prefix.c_str(),"ge21_phierror_o_avg_phierror",";GE21 #sigma_{#phi}/ <#sigma_{#phi}>",10000,0.001,1000)->Fill((minGE21Error < 0 || avgOtherError <= 0 ) ? -1.0 : minGE21Error/avgOtherError);
            plotter.getOrMake1DPre(prefix.c_str(),"ge21_phierror_o_min_phierror",";GE21 #sigma_{#phi}/ <#sigma_{#phi}>",10000,0.001,1000)->Fill((minGE21Error < 0 || minOtherError <= 0 ) ? -1.0 : minGE21Error/minOtherError);

            int nOther = filledME0 + filledME1+filledME2+filledME3+filledME4+filledGE1;

            plotter.getOrMake1DPre(prefix.c_str(),"nOtherStations",";N. other stations",10,-0.5,9.5)->Fill(nOther);
        };

        auto makePlot2 = [&](const std::string& prefix, const std::vector<int8>& v_d,const std::vector<int8>& v_s, const std::vector<int8>& v_r, const std::vector<float>& v_phie){
            const float eta = std::fabs(genMuon_eta);
            if(eta > 1.6 && eta < 2.4) fillInfo(prefix+"_eta_1p6_2p4",v_d,v_s,v_r,v_phie);
            if(eta > 1.6 && eta < 1.8) fillInfo(prefix+"_eta_1p6_1p8",v_d,v_s,v_r,v_phie);
            if(eta > 1.8 && eta < 2.0) fillInfo(prefix+"_eta_1p8_2p0",v_d,v_s,v_r,v_phie);
            if(eta > 2.0 && eta < 2.2) fillInfo(prefix+"_eta_2p0_2p2",v_d,v_s,v_r,v_phie);
            if(eta > 2.2 && eta < 2.4) fillInfo(prefix+"_eta_2p2_2p4",v_d,v_s,v_r,v_phie);

        };

        makePlot2(prefix+"_globalrh",*global_v_d,*global_v_s,*global_v_r,*global_v_phie);
        makePlot2(prefix+"_pickyrh",*picky_v_d,*picky_v_s,*picky_v_r,*picky_v_phie);


    }

    const int DT= 1;
    const int CSC=2;
    const int RPC=3;
    const int GEM=4;
    const int ME0=5;

    void write(std::string fileName){ plotter.write(fileName);}
    HistGetter plotter;
};

#endif

void AnalyzeGE21MomentumTree(std::string fileName,std::string sampName){
    Analyzer a(fileName,"Events",sampName);
    a.analyze();
    a.write(std::string("granPlots_")+sampName+".root");
}
