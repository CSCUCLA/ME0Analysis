
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

using namespace std;
using namespace PhysicsUtilities;
using ASTypes::CylLorentzVectorF;
using ASTypes::size;
using ASTypes::size8;
const double CSV_LOOSE        = 0.460;
const double CSV_MEDIUM       = 0.800;
const double CSV_TIGHT        = 0.935;

template<typename Momentum1, typename Momentum2>
double transverseMass(const Momentum1& visible, const Momentum2& invisible)
{
  const double    cosDPhi   = TMath::Cos( PhysicsUtilities::deltaPhi(visible.phi(), invisible.phi()) );
  return TMath::Sqrt( 2 * visible.pt() * invisible.pt() * (1 - cosDPhi) );
}


CylLorentzVectorF getnz(const CylLorentzVectorF& met, const CylLorentzVectorF& vis) {
	const double mH = 125;
	const double a = mH*mH - vis.mass()*vis.mass() +2*vis.x()*met.x() +2*vis.y()*met.y();
	const double A = 4*(vis.E()*vis.E() - vis.z()*vis.z());
	const double B = -4*a* vis.z();
	const double C = 4*vis.E()*vis.E()*(met.x()*met.x() + met.y()*met.y()) - a*a;
	const double delta = B*B -4*A*C;

	double metZ = 0;
	if(delta < 0) {
		metZ= -B/(2*A);
	} else {
		const double pos = (-B + std::sqrt(delta))/(2*A);
		const double neg = (-B - std::sqrt(delta))/(2*A);
		if(std::fabs(pos) < std::fabs(neg)) metZ = pos;
		else metZ = neg;
	}
	ASTypes::CartLorentzVector neutrino(met.x(),met.y(),metZ,std::sqrt(met.x()*met.x()+met.y()*met.y()+metZ*metZ));
	return CylLorentzVectorF(neutrino);
};


class Jet : public CylLorentzVectorF{
public:
	Jet() : CylLorentzVectorF() {};
	Jet(const CylLorentzVectorF& mom) : CylLorentzVectorF(mom) {}
	float csv = -1;
};

class FatJet : public Jet{
public:
	FatJet() : Jet() {}
	FatJet(const CylLorentzVectorF& mom) : Jet(mom) {}
	float tau2otau1() const {return tau1 == 0 ? 99 : tau2/tau1;}
	float tau3otau1() const {return tau1 == 0 ? 99 : tau3/tau1;}
	float tau3otau2() const {return tau2 == 0 ? 99 : tau3/tau2;}
	float minSJCSV() const {return sj1.pt() > 0 && sj2.pt() > 0 ? std::min(sj1.csv,sj2.csv) : 0.0;}
	float maxSJCSV() const {return sj1.pt() > 0 && sj2.pt() > 0 ? std::max(sj1.csv,sj2.csv) : 0.0;}
	float sd_massFFJ=-1;
	float bbcsv    =-1;
	float tau1     =-1;
	float tau2     =-1;
	float tau3     =-1;
	CylLorentzVectorF sdMom;
	Jet sj1;
	Jet sj2;
};
class BaseCutAnalyzer : public BaseTupleAnalyzer{
public:
	TString glbPrefix = "";


	BaseCutAnalyzer(std::string fileName, std::string treeName) : BaseTupleAnalyzer(fileName,treeName){

		setBranchAddress("decayType"  ,&decayType  );
		setBranchAddress("hbb_pt"     ,&hbb_pt     );
		setBranchAddress("hbb_eta"    ,&hbb_eta    );
		setBranchAddress("hbb_phi"    ,&hbb_phi    );
		setBranchAddress("hbb_mass"   ,&hbb_mass   );
		setBranchAddress("b1_pt"      ,&b1_pt      );
		setBranchAddress("b1_eta"     ,&b1_eta     );
		setBranchAddress("b1_phi"     ,&b1_phi     );
		setBranchAddress("b1_mass"    ,&b1_mass    );
		setBranchAddress("b2_pt"      ,&b2_pt      );
		setBranchAddress("b2_eta"     ,&b2_eta     );
		setBranchAddress("b2_phi"     ,&b2_phi     );
		setBranchAddress("b2_mass"    ,&b2_mass    );
		setBranchAddress("hww_pt"     ,&hww_pt     );
		setBranchAddress("hww_eta"    ,&hww_eta    );
		setBranchAddress("hww_phi"    ,&hww_phi    );
		setBranchAddress("hww_mass"   ,&hww_mass   );
		setBranchAddress("w1d1_pt"    ,&w1d1_pt    );
		setBranchAddress("w1d1_eta"   ,&w1d1_eta   );
		setBranchAddress("w1d1_phi"   ,&w1d1_phi   );
		setBranchAddress("w1d1_mass"  ,&w1d1_mass  );
		setBranchAddress("w1d2_pt"    ,&w1d2_pt    );
		setBranchAddress("w1d2_eta"   ,&w1d2_eta   );
		setBranchAddress("w1d2_phi"   ,&w1d2_phi   );
		setBranchAddress("w1d2_mass"  ,&w1d2_mass  );
		setBranchAddress("w2d1_pt"    ,&w2d1_pt    );
		setBranchAddress("w2d1_eta"   ,&w2d1_eta   );
		setBranchAddress("w2d1_phi"   ,&w2d1_phi   );
		setBranchAddress("w2d1_mass"  ,&w2d1_mass  );
		setBranchAddress("w2d2_pt"    ,&w2d2_pt    );
		setBranchAddress("w2d2_eta"   ,&w2d2_eta   );
		setBranchAddress("w2d2_phi"   ,&w2d2_phi   );
		setBranchAddress("w2d2_mass"  ,&w2d2_mass  );

		setBranchAddress("t1_b_pt"    ,&t1_b_pt      );
		setBranchAddress("t1_b_eta"   ,&t1_b_eta     );
		setBranchAddress("t1_b_phi"   ,&t1_b_phi     );
		setBranchAddress("t1_w1_pt"   ,&t1_w1_pt     );
		setBranchAddress("t1_w1_eta"  ,&t1_w1_eta    );
		setBranchAddress("t1_w1_phi"  ,&t1_w1_phi    );
		setBranchAddress("t1_w1_pdgid",&t1_w1_pdgid  );
		setBranchAddress("t1_w2_pt"   ,&t1_w2_pt     );
		setBranchAddress("t1_w2_eta"  ,&t1_w2_eta    );
		setBranchAddress("t1_w2_phi"  ,&t1_w2_phi    );
		setBranchAddress("t1_w2_pdgid",&t1_w2_pdgid  );
		setBranchAddress("t2_b_pt"    ,&t2_b_pt      );
		setBranchAddress("t2_b_eta"   ,&t2_b_eta     );
		setBranchAddress("t2_b_phi"   ,&t2_b_phi     );
		setBranchAddress("t2_w1_pt"   ,&t2_w1_pt     );
		setBranchAddress("t2_w1_eta"  ,&t2_w1_eta    );
		setBranchAddress("t2_w1_phi"  ,&t2_w1_phi    );
		setBranchAddress("t2_w1_pdgid",&t2_w1_pdgid  );
		setBranchAddress("t2_w2_pt"   ,&t2_w2_pt     );
		setBranchAddress("t2_w2_eta"  ,&t2_w2_eta    );
		setBranchAddress("t2_w2_phi"  ,&t2_w2_phi    );
		setBranchAddress("t2_w2_pdgid",&t2_w2_pdgid  );

		setBranchAddress("weight"      ,&weight      );
		setBranchAddress("process"     ,&process     );
		setBranchAddress("ht"          ,&ht          );
		setBranchAddress("met"         ,&met         );
		setBranchAddress("met_phi"     ,&met_phi     );
		setBranchAddress("lepton_pt"   ,&lepton_pt   );
		setBranchAddress("lepton_eta"  ,&lepton_eta  );
		setBranchAddress("lepton_phi"  ,&lepton_phi  );
		setBranchAddress("lepton_muon" ,&lepton_muon );

		setBranchAddress("fj_pt"       ,&fj_pt       );
		setBranchAddress("fj_eta"      ,&fj_eta      );
		setBranchAddress("fj_phi"      ,&fj_phi      );
		setBranchAddress("fj_mass"     ,&fj_mass     );
		setBranchAddress("fj_sd_massFFJ",&fj_sd_massFFJ);
		setBranchAddress("fj_csv"      ,&fj_csv      );
		setBranchAddress("fj_bbcsv"    ,&fj_bbcsv    );
		setBranchAddress("fj_tau1",&fj_tau1     );
		setBranchAddress("fj_tau2",&fj_tau2     );
		setBranchAddress("fj_tau3",&fj_tau3     );
		setBranchAddress("fj_sd_pt"    ,&fj_sd_pt    );
		setBranchAddress("fj_sd_eta"   ,&fj_sd_eta   );
		setBranchAddress("fj_sd_phi"   ,&fj_sd_phi   );
		setBranchAddress("fj_sd_mass"  ,&fj_sd_mass  );
		setBranchAddress("fj_sj1_pt"   ,&fj_sj1_pt   );
		setBranchAddress("fj_sj1_eta"  ,&fj_sj1_eta   );
		setBranchAddress("fj_sj1_phi"  ,&fj_sj1_phi   );
		setBranchAddress("fj_sj1_mass" ,&fj_sj1_mass   );
		setBranchAddress("fj_sj1_csv"  ,&fj_sj1_csv  );
		setBranchAddress("fj_sj2_pt"   ,&fj_sj2_pt   );
		setBranchAddress("fj_sj2_eta"  ,&fj_sj2_eta   );
		setBranchAddress("fj_sj2_phi"  ,&fj_sj2_phi   );
		setBranchAddress("fj_sj2_mass" ,&fj_sj2_mass   );
		setBranchAddress("fj_sj2_csv"  ,&fj_sj2_csv  );

		setBranchAddress("fj1p5_pt"       ,&fj1p5_pt       );
		setBranchAddress("fj1p5_eta"      ,&fj1p5_eta      );
		setBranchAddress("fj1p5_phi"      ,&fj1p5_phi      );
		setBranchAddress("fj1p5_mass"     ,&fj1p5_mass     );
		setBranchAddress("fj1p5_sd_massFFJ"     ,&fj1p5_sd_massFFJ);
		setBranchAddress("fj1p5_csv"      ,&fj1p5_csv      );
		setBranchAddress("fj1p5_bbcsv"      ,&fj1p5_bbcsv    );
		setBranchAddress("fj1p5_tau1",&fj1p5_tau1     );
		setBranchAddress("fj1p5_tau2",&fj1p5_tau2     );
		setBranchAddress("fj1p5_tau3",&fj1p5_tau3     );
		setBranchAddress("fj1p5_sd_pt"    ,&fj1p5_sd_pt    );
		setBranchAddress("fj1p5_sd_eta"   ,&fj1p5_sd_eta   );
		setBranchAddress("fj1p5_sd_phi"   ,&fj1p5_sd_phi   );
		setBranchAddress("fj1p5_sd_mass"  ,&fj1p5_sd_mass  );
		setBranchAddress("fj1p5_sj1_pt"   ,&fj1p5_sj1_pt   );
		setBranchAddress("fj1p5_sj1_eta"  ,&fj1p5_sj1_eta   );
		setBranchAddress("fj1p5_sj1_phi"  ,&fj1p5_sj1_phi   );
		setBranchAddress("fj1p5_sj1_mass" ,&fj1p5_sj1_mass   );
		setBranchAddress("fj1p5_sj1_csv"  ,&fj1p5_sj1_csv  );
		setBranchAddress("fj1p5_sj2_pt"   ,&fj1p5_sj2_pt   );
		setBranchAddress("fj1p5_sj2_eta"  ,&fj1p5_sj2_eta   );
		setBranchAddress("fj1p5_sj2_phi"  ,&fj1p5_sj2_phi   );
		setBranchAddress("fj1p5_sj2_mass" ,&fj1p5_sj2_mass   );
		setBranchAddress("fj1p5_sj2_csv"  ,&fj1p5_sj2_csv  );

		setBranchAddress("jet_pt"      ,&jet_pt      );
		setBranchAddress("jet_eta"     ,&jet_eta     );
		setBranchAddress("jet_phi"     ,&jet_phi     );
		setBranchAddress("jet_mass"    ,&jet_mass    );
		setBranchAddress("jet_csv"     ,&jet_csv     );


	}


	size8        decayType      =0    ;
	float        hbb_pt         =0    ;
	float        hbb_eta        =0    ;
	float        hbb_phi        =0    ;
	float        hbb_mass       =0    ;
	float        b1_pt          =0   ;
	float        b1_eta         =0   ;
	float        b1_phi         =0   ;
	float        b1_mass        =0   ;
	float        b2_pt          =0   ;
	float        b2_eta         =0   ;
	float        b2_phi         =0   ;
	float        b2_mass        =0   ;
	float        hww_pt         =0    ;
	float        hww_eta        =0    ;
	float        hww_phi        =0    ;
	float        hww_mass       =0    ;
	float        w1d1_pt        =0     ;
	float        w1d1_eta       =0     ;
	float        w1d1_phi       =0     ;
	float        w1d1_mass      =0     ;
	float        w1d2_pt        =0     ;
	float        w1d2_eta       =0     ;
	float        w1d2_phi       =0     ;
	float        w1d2_mass      =0     ;
	float        w2d1_pt        =0     ;
	float        w2d1_eta       =0     ;
	float        w2d1_phi       =0     ;
	float        w2d1_mass      =0     ;
	float        w2d2_pt        =0     ;
	float        w2d2_eta       =0     ;
	float        w2d2_phi       =0     ;
	float        w2d2_mass      =0     ;

	float        t1_b_pt        =0;
	float        t1_b_eta       =0;
	float        t1_b_phi       =0;
	float        t1_w1_pt       =0;
	float        t1_w1_eta      =0;
	float        t1_w1_phi      =0;
	int          t1_w1_pdgid    =0;
	float        t1_w2_pt       =0;
	float        t1_w2_eta      =0;
	float        t1_w2_phi      =0;
	int          t1_w2_pdgid    =0;
	float        t2_b_pt        =0;
	float        t2_b_eta       =0;
	float        t2_b_phi       =0;
	float        t2_w1_pt       =0;
	float        t2_w1_eta      =0;
	float        t2_w1_phi      =0;
	int          t2_w1_pdgid    =0;
	float        t2_w2_pt       =0;
	float        t2_w2_eta      =0;
	float        t2_w2_phi      =0;
	int          t2_w2_pdgid    =0;

	float        weight         =0  ;
	size8        process        =0   ;
	float        ht             =0  ;
	float        met            =0  ;
	float        met_phi        =0         ;
	std::vector<float>  * lepton_pt    = new std::vector<float>   ;
	std::vector<float>  * lepton_eta   = new std::vector<float>    ;
	std::vector<float>  * lepton_phi   = new std::vector<float>    ;
	std::vector<size8>  * lepton_muon  = new std::vector<size8  >     ;

	std::vector<float>  * fj_pt        = new std::vector<float>    ;
	std::vector<float>  * fj_eta       = new std::vector<float>    ;
	std::vector<float>  * fj_phi       = new std::vector<float>    ;
	std::vector<float>  * fj_mass      = new std::vector<float>    ;
	std::vector<float>  * fj_sd_massFFJ= new std::vector<float>    ;
	std::vector<float>  * fj_csv       = new std::vector<float>    ;
	std::vector<float>  * fj_bbcsv     = new std::vector<float>    ;
	std::vector<float>  * fj_tau1      = new std::vector<float>    ;
	std::vector<float>  * fj_tau2      = new std::vector<float>    ;
	std::vector<float>  * fj_tau3      = new std::vector<float>    ;
	std::vector<float>  * fj_sd_pt     = new std::vector<float>    ;
	std::vector<float>  * fj_sd_eta    = new std::vector<float>    ;
	std::vector<float>  * fj_sd_phi    = new std::vector<float>    ;
	std::vector<float>  * fj_sd_mass   = new std::vector<float>    ;
	std::vector<float>  * fj_sj1_pt    = new std::vector<float>    ;
	std::vector<float>  * fj_sj1_eta   = new std::vector<float>    ;
	std::vector<float>  * fj_sj1_phi   = new std::vector<float>    ;
	std::vector<float>  * fj_sj1_mass  = new std::vector<float>    ;
	std::vector<float>  * fj_sj1_csv  = new std::vector<float>    ;
	std::vector<float>  * fj_sj2_pt    = new std::vector<float>    ;
	std::vector<float>  * fj_sj2_eta   = new std::vector<float>    ;
	std::vector<float>  * fj_sj2_phi   = new std::vector<float>    ;
	std::vector<float>  * fj_sj2_mass  = new std::vector<float>    ;
	std::vector<float>  * fj_sj2_csv   = new std::vector<float>    ;

	std::vector<float>  * fj1p5_pt        = new std::vector<float>    ;
	std::vector<float>  * fj1p5_eta       = new std::vector<float>    ;
	std::vector<float>  * fj1p5_phi       = new std::vector<float>    ;
	std::vector<float>  * fj1p5_mass      = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sd_massFFJ= new std::vector<float>    ;
	std::vector<float>  * fj1p5_csv       = new std::vector<float>    ;
	std::vector<float>  * fj1p5_bbcsv     = new std::vector<float>    ;
	std::vector<float>  * fj1p5_tau1      = new std::vector<float>    ;
	std::vector<float>  * fj1p5_tau2      = new std::vector<float>    ;
	std::vector<float>  * fj1p5_tau3      = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sd_pt     = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sd_eta    = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sd_phi    = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sd_mass   = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj1_pt    = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj1_eta   = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj1_phi   = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj1_mass  = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj1_csv  = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj2_pt    = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj2_eta   = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj2_phi   = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj2_mass  = new std::vector<float>    ;
	std::vector<float>  * fj1p5_sj2_csv   = new std::vector<float>    ;

	std::vector<float>  * jet_pt       = new std::vector<float>    ;
	std::vector<float>  * jet_eta      = new std::vector<float>    ;
	std::vector<float>  * jet_phi      = new std::vector<float>    ;
	std::vector<float>  * jet_mass     = new std::vector<float>    ;
	std::vector<float>  * jet_csv      = new std::vector<float>    ;

	const float lumi = 36.0;
	float nWeight = 0;
	float sigWeight = 0;

	void write(TString fileName){ plotter.write(fileName);}
	HistGetter plotter;

	CylLorentzVectorF getMET() const {return CylLorentzVectorF(met,0,met_phi,0);}
	CylLorentzVectorF getLep(const unsigned int idx = 0) const {return CylLorentzVectorF(lepton_pt->at(idx),lepton_eta->at(idx),lepton_phi->at(idx),lepton_muon->at(idx) ? 0.1 : 0);}
	CylLorentzVectorF getFJ(const unsigned int idx) const {return CylLorentzVectorF(fj_pt->at(idx),fj_eta->at(idx),fj_phi->at(idx),fj_mass->at(idx));}
	CylLorentzVectorF getSDFJ(const unsigned int idx) const {return CylLorentzVectorF(fj_sd_pt->at(idx),fj_sd_eta->at(idx),fj_sd_phi->at(idx),fj_sd_mass->at(idx));}
	CylLorentzVectorF getJet(const unsigned int idx) const {return CylLorentzVectorF(jet_pt->at(idx),jet_eta->at(idx),jet_phi->at(idx),jet_mass->at(idx));}

	struct genTop {
		CylLorentzVectorF top;
		CylLorentzVectorF b;
		CylLorentzVectorF W;
		CylLorentzVectorF w1;
		CylLorentzVectorF w2;
	};
	std::pair<genTop,genTop> getTTBar(){
		auto getMass = [](int pdgID) -> float {
			switch (std::abs(pdgID)){
			case 13 :
				return .1;
			case 15 :
				return 1.8;
			case 5:
				return 4.8;
			case 6:
				return 172.44;
			default:
				return 0;
			}
		};
		genTop ta;
		ta.b = CylLorentzVectorF(t1_b_pt,t1_b_eta,t1_b_phi,4.8);
		ta.w1 = CylLorentzVectorF(t1_w1_pt,t1_w1_eta,t1_w1_phi,getMass(t1_w1_pdgid));
		ta.w2 = CylLorentzVectorF(t1_w2_pt,t1_w2_eta,t1_w2_phi,getMass(t1_w2_pdgid));
		ta.W = ta.w1+ta.w2;
		ta.top = ta.b +ta.W;
		genTop tb;
		tb.b = CylLorentzVectorF(t2_b_pt,t2_b_eta,t2_b_phi,4.8);
		tb.w1 = CylLorentzVectorF(t2_w1_pt,t2_w1_eta,t2_w1_phi,getMass(t2_w1_pdgid));
		tb.w2 = CylLorentzVectorF(t2_w2_pt,t2_w2_eta,t2_w2_phi,getMass(t2_w2_pdgid));
		tb.W = tb.w1+tb.w2;
		tb.top = tb.b +tb.W;
		return std::make_pair(ta,tb);
	}

	void fillJets(std::vector<Jet>& jets){
		jets.clear();
		for(unsigned int idx = 0; idx < jet_pt->size(); ++idx ){
			jets.emplace_back(CylLorentzVectorF(jet_pt->at(idx),jet_eta->at(idx),jet_phi->at(idx),jet_mass->at(idx)));
			jets.back().csv = jet_csv->at(idx);
		}
	}

	void fillFJJets(std::vector<FatJet>& jets){
		jets.clear();
		for(unsigned int idx = 0; idx < fj_pt->size(); ++idx ){
			jets.emplace_back(CylLorentzVectorF(fj_pt->at(idx),fj_eta->at(idx),fj_phi->at(idx),fj_mass->at(idx)));
			auto & jet = jets.back();
			jet.sd_massFFJ =fj_sd_massFFJ->at(idx);
			jet.csv        =fj_csv       ->at(idx);
			jet.bbcsv      =fj_bbcsv     ->at(idx);
			jet.tau1       =fj_tau1      ->at(idx);
			jet.tau2       =fj_tau2      ->at(idx);
			jet.tau3       =fj_tau3      ->at(idx);
			jet.sj1      =CylLorentzVectorF(fj_sj1_pt->at(idx),fj_sj1_eta->at(idx),fj_sj1_phi->at(idx),fj_sj1_mass->at(idx)); jet.sj1.csv = fj_sj1_csv->at(idx);
			jet.sj2      =CylLorentzVectorF(fj_sj2_pt->at(idx),fj_sj2_eta->at(idx),fj_sj2_phi->at(idx),fj_sj2_mass->at(idx)); jet.sj2.csv = fj_sj2_csv->at(idx);
			jet.sdMom = jet.sj1 + jet.sj2;
		}
	}
	void fillFJ1p5Jets(std::vector<FatJet>& jets){
		jets.clear();
		for(unsigned int idx = 0; idx < fj1p5_pt->size(); ++idx ){
			jets.emplace_back(CylLorentzVectorF(fj1p5_pt->at(idx),fj1p5_eta->at(idx),fj1p5_phi->at(idx),fj1p5_mass->at(idx)));
			auto & jet = jets.back();
			jet.sd_massFFJ =fj1p5_sd_massFFJ->at(idx);
			jet.csv        =fj1p5_csv       ->at(idx);
			jet.bbcsv      =fj1p5_bbcsv     ->at(idx);
			jet.tau1       =fj1p5_tau1      ->at(idx);
			jet.tau2       =fj1p5_tau2      ->at(idx);
			jet.tau3       =fj1p5_tau3      ->at(idx);
			jet.sj1      =CylLorentzVectorF(fj1p5_sj1_pt->at(idx),fj1p5_sj1_eta->at(idx),fj1p5_sj1_phi->at(idx),fj1p5_sj1_mass->at(idx)); jet.sj1.csv = fj1p5_sj1_csv->at(idx);
			jet.sj2      =CylLorentzVectorF(fj1p5_sj2_pt->at(idx),fj1p5_sj2_eta->at(idx),fj1p5_sj2_phi->at(idx),fj1p5_sj2_mass->at(idx)); jet.sj2.csv = fj1p5_sj2_csv->at(idx);
			jet.sdMom = jet.sj1 + jet.sj2;
		}
	}

	bool isSignal() const {return process == 10;}
};
#endif

