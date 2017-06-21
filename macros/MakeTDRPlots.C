
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
#include "TString.h"
#include "TLine.h"
#include "HistoPlotting/include/Plotter.h"
#include "HistoPlotting/include/PlotTools.h"
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "HistoPlotting/include/StyleInfo.h"
using namespace std;

class Analyzer {
public:
	TString prefix = "p8s384";
	int dottedStyle = 1;
	int straightStyle = 1;
//	vector<int> colors = {800,1,634,857};
//	vector<int> syles  = {9,1,9,1};

	vector<int> colors = {800,634,1,857};
	vector<int> syles  = {9,9,1,1};

	Analyzer() {}

  void makeResolutionPlots(TString filename) {
	  TFile * fF = new TFile(filename,"READ");
	  TString pdfNames[] = {"etaRes.pdf","dPhiRes.pdf"};
	  TString vars[] = {"eta","dphi",""};
	  TString recos[] = {"track","segment","gen","",""};
	  TString recoNs[] = {"projected track","segment","gen","",""};


	  TString secVars[] = {"genmreco",""};

	  TString pts[] = {"pteq3to5","ptgeq20","","",""};
	  TString ptNs[] = {"p_{T} 3-5 GeV","p_{T} 20-30 GeV",""};



	        for(unsigned int iV = 0; vars[iV][0]; ++ iV){
	          for(unsigned int iS = 0; secVars[iS][0]; ++ iS){
	            Plotter * p = new Plotter();
	            int iC = 0;
	            for(unsigned int iR = 0; recos[iR][0]; ++ iR)
	            for(unsigned int iP = 0; pts[iP][0]; ++ iP){
	              TH1 *hf = 0;
	              fF->GetObject(TString::Format("%s_%s_real_muon_%s_%s_%s",prefix.Data(),pts[iP].Data(),recos[iR].Data(), secVars[iS].Data(),vars[iV].Data()),hf);
	              if(hf == 0) continue;
	              p->addHistLine(hf,TString::Format("%s, %s",recoNs[iR].Data(),ptNs[iP].Data()),colors[iC],syles[iC],4 );
	              iC++;
	            }
	            p->normalize();
	            p->rebin(4);
	            p->setYTitle("arbitrary units");
	            p->setCanvasSize(1024,726);
	            p->setMargins(0.08,.12,.12,.05);
	            p->setAxisTextSize(0);
	            p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
	            if(iV == 1) p->setMinMax(0.,0.46);
	            else p->setMinMax(0.,0.52);
	            p->setLegendPos(0.13,0.71,0.6,0.89);

	            p->draw(true,pdfNames[iV]);

	          }
	        }

	        fF->Close();
  }

  void makeLayersCrossedPlot(TString filename){

	      TFile * ff = new TFile(filename);
	      TH1 * hEv = 0;
	      ff->GetObject("nEvents",hEv);
	      float nEvents = hEv->GetBinContent(1);

	      TH1 * hM= 0;
	      TH1 * hE= 0;
	      TH1 * hH= 0;

	      ff->GetObject("muon_nLaysHit",hM);  hM->Scale(200./nEvents);
	      ff->GetObject("ele_nLaysHit",hE);   hE->Scale(200./nEvents);
	      ff->GetObject("hadron_nLaysHit",hH);hH->Scale(200./nEvents);
	      hE->GetXaxis()->SetRangeUser(0.5,6.5);

	      hE->SetXTitle("N. of ME0 layers crossed by particle");
	      hE->SetYTitle("<N. of bkg. particles> per BX / bin");


	      Plotter* p = new Plotter;
	      p->addStackHist(hE,"bkg. electrons");
	      p->addStackHist(hH,"bkg. hadrons");
	      p->addStackHist(hM,"bkg. muons",StyleInfo::defFillColors[3]);

          p->setCanvasSize(1024,726);
          p->setMargins(0.08,.12,.12,.05);
          p->setAxisTextSize(0);
          p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
          p->setLegendPos(0.6,0.72,1.0,0.89);
          p->addText("2.0 < |#eta| < 2.8",0.613,0.675,0.04);


	      TCanvas * c = p->draw(false,"bkgParticleLayerNumber.pdf");
	      p->xAxis()->SetRangeUser(0.5,6.5);
	      c->Print("bkgParticleLayerNumber.pdf");

	      ff->Close();


  }

  void makeDPhiPlot(TString fileBkg, TString fileSignal){
	    TFile * ff = new TFile(fileBkg);
	    TH1 * hfa = 0;
	    TH1 * hf3= 0;
	    TH1 * hf4= 0;
	    TH1 * hf4m= 0;
	    TH1 * hf5= 0;
	    TH1 * hf6= 0;
	    TH1 * hf6m= 0;
	    ff->GetObject("all_tr_signed_dPhi",hfa);
	    ff->GetObject("all_nLays_geq3_tr_signed_dPhi",hf3);
	    ff->GetObject("all_nLays_geq4_tr_signed_dPhi",hf4);
	    ff->GetObject("muon_nLays_geq4_tr_signed_dPhi",hf4m);
	    ff->GetObject("all_nLays_geq5_tr_signed_dPhi",hf5);
	    ff->GetObject("all_nLays_geq6_tr_signed_dPhi",hf6);
	    ff->GetObject("muon_nLays_geq6_tr_signed_dPhi",hf6m);

	    TH1 * hf6o = (TH1*)hf6->Clone();

	    TH1 * hf5o = (TH1*)hf5->Clone();
	    hf5o->Add(hf6,-1);

	    TH1 * hf4o = (TH1*)hf4->Clone();
	    hf4o->Add(hf5,-1);

	    TH1 * hf3o = (TH1*)hf3->Clone();
	    hf3o->Add(hf4,-1);

	    TH1 * hflt = (TH1*)hfa->Clone();
	    hflt->Add(hf3,-1);


	    TH1 * hf34o = (TH1*)hf3->Clone();
	    hf34o->Add(hf5,-1);

	    TH1 * hf56o = (TH1*)hf5->Clone();

	    TH1 * hflt4 = (TH1*)hfa->Clone();
	    hflt4->Add(hf4,-1);

	    TH1 * hf4nonM = (TH1*)hf4->Clone();
	    hf4nonM->Add(hf4m,-1);

	    TH1 * hf6nonM = (TH1*)hf6->Clone();
	    hf6nonM->Add(hf6m,-1);


	    TFile * fs = new TFile(fileSignal);
	    TH1 * hr0 =0;
	    TH1 * hr1 =0;
	    TH1 * hr3=0;
	    TH1 * hr5=0;
	    TH1 * hr20=0;
	    fs->GetObject("p8s384_ptleq2__tr_geq4_signed_dPhi",hr0);
	    fs->GetObject("p8s384_pteq2to3__tr_geq4_signed_dPhi",hr1);
	    fs->GetObject("p8s384_pteq3to5__tr_geq4_signed_dPhi",hr3);
	    fs->GetObject("p8s384_pteq10to15__tr_geq4_signed_dPhi",hr5);
	    fs->GetObject("p8s384_ptgeq20__tr_geq4_signed_dPhi",hr20);


	    Plotter * p = new Plotter;

	    p->addStackHist(hf4m,"bkg. muons",-1,1001,1,2);
	    p->addStackHist(hf4nonM,"bkg. electrons and hadrons",-1,1001,1,2);


	    p->addHistLine(hr1 ,"signal muons, p_{T} 2-3 GeV",colors[0],syles[0],4);
	    p->addHistLine(hr3 ,"signal muons, p_{T} 3-5 GeV",colors[1],syles[1],4);
	    p->addHistLine(hr5 ,"signal muons, p_{T} 10-15 GeV",colors[2],syles[2],4);
	    p->addHistLine(hr20,"signal muons, p_{T} 20-30 GeV",colors[3],syles[3],4);
	    p->rebin(1);
	    p->normalize();

	    p->setXTitle("sim -#it{q}#times#Delta#phi");
	    p->setYTitle("arbitrary units");

	    p->setMinMax(0,0.65);

        p->setCanvasSize(1024,726);
        p->setMargins(0.08,.12,.12,.05);
        p->setAxisTextSize(0);
        p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
        p->setLegendPos(0.4,0.55,.9,0.89);
        p->addText("N. of ME0 layers crossed by particle #geq 4",0.416,0.50,0.04);
	      TCanvas * c = p->draw(false,"generatedDPhi_signalAndBkg.pdf");
	      c->Print("generatedDPhi_signalAndBkg.pdf");

	    ff->Close();
	    fs->Close();
	}


  void makeNeutronPlot(){

	  vector<TString> xTitles{"7.5 (x1)","10 (x1.33)","22.5 (x3)","37.5 (x5)","75 (x10)"};
	  vector<TString> names{"4 ME0 layers per segment","5 ME0 layers per segment"};
	  vector<vector<double>> vals{{0.198,0.282,5.298,29.626,195.550},{0.008,0.023,0.214,1.748,39.374}};
	  Plotter * p = new Plotter();

	  TH1F * axisHist = new TH1F("axisHist",";neutron luminosity scale [10^{34} Hz/cm] (safety margin)",xTitles.size(),-0.5,float(xTitles.size()) -.5);
	  for (unsigned int i=1;i<=xTitles.size();i++) axisHist->GetXaxis()->SetBinLabel(i,xTitles[i-1]);

	  for(unsigned int iN = 0; iN < names.size(); ++iN){
	    TGraph * g = new TGraph();
	    for(unsigned int iP = 0; iP < vals[iN].size(); ++iP){
	      g->SetPoint(iP,iP,vals[iN][iP]);
	    }
	    g->SetLineColor  (StyleInfo::getLineColor(iN));
	    g->SetLineWidth  (4);
	    g->SetLineStyle  (1);
	    g->SetMarkerStyle(20);
	    g->SetMarkerColor(StyleInfo::getLineColor(iN));
	    g->SetMarkerSize (1);
	    Drawing::Drawable1D drawableF("P L",names[iN],Drawing::GRAPH,g,false);
	    drawableF.graphAxisHist = (TH1*)axisHist->Clone();
	    p->addDrawable(drawableF);
	  }
	  p->setYTitle("<N. of neutron bkg. segments> per BX");
	  p->setXTitle("neutron luminosity scale [10^{34}cm^{-2}s^{-1}] (safety margin)");

          p->setCanvasSize(1024,726);
          p->setMargins(0.08,.12,.12,.05);
          p->setAxisTextSize(0);
          p->setCMSLumi(0, "14 TeV", "Simulation preliminary",1.2 );
          p->setLegendPos(.13,0.78,0.65,0.89);
          p->addText("2.0 < |#eta| < 2.8",0.15,0.74,0.04);
          p->addText("|#Delta#phi| < 0.013",0.15,0.70,0.04);

	      TCanvas * c = p->draw(false,"neutronSegmentMultiplicity.pdf");
	      p->yAxis()->SetRangeUser(0.004,400);
	      c->SetLogy();



	      c->Print("neutronSegmentMultiplicity.pdf");
  }


  void makeSegmentBkgPlots(TString filename) {
	  TFile * f = new TFile(filename,"READ");

	  TString digi = "p8s384";

	  TString dPhis[] = {"_","_dPhi_lt0p013_","_dPhi_lt0p004_","_dPhi_lt0p002_","","",""};
	  TString dPhiNs[] = {"|#Delta#phi| < 0.02 (p_{T} > 2 GeV)","|#Delta#phi| < 0.013 (p_{T} > 3 GeV)","|#Delta#phi| < 0.004 (p_{T} > 10 GeV)","|#Delta#phi| < 0.002 (p_{T} > 20 GeV)","","",""};


	  TString vars[] ={ "seg_eta",""};
	  TH1 * hInt = 0;
	  f->GetObject(TString::Format("%s___fake_nEvents",digi.Data()),hInt);
	  float nEvents = hInt->GetBinContent(1);


	  for(unsigned int iV = 0; vars[iV][0]; ++iV){
		  Plotter * p = new Plotter;
		  for(unsigned int iD = 0; dPhis[iD][0]; ++iD){
			  TH1 * h;
			  f->GetObject(TString::Format("%s___fake%s%s",digi.Data(),dPhis[iD].Data(),vars[iV].Data()),h);
			  if(h == 0) continue;
			  h = (TH1*)h->Clone();
			  h->Rebin(10);
			  h->Scale(1./nEvents);
			  h->SetYTitle("<N. of bkg. segments> per BX / 0.1");
			  p->addHistLine(h,dPhiNs[iD]);
		  }

		    p->setCanvasSize(1024,726);
		    p->setMargins(0.08,.12,.12,.05);
		    p->setAxisTextSize(0);
		    p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
		    p->setLegendPos(.13,0.70,0.65,0.89);
		    p->addText("2.0 < |#eta| < 2.8",0.15,0.66,0.04);

		      TCanvas * c = p->draw(false,"segmentDensity_byETA.pdf");
		      c->Print("segmentDensity_byETA.pdf");

	  }
	  f->Close();
  }

  void makeMatchBkgPlots(TString filename) {
	  TString prefix = "p8s384";
	  TFile * fF = new TFile(filename,"READ");

	  TString extr = "_";
	  TString var = "eta";
	  TString xTitle = "pixel track |#eta|";

	  TString pts[] = {"ptgeq2","ptgeq3","ptgeq5","",""};
	  TString ptNs[] = {"pixel track p_{T} > 2 GeV","pixel track p_{T} > 3 GeV","pixel track p_{T} > 5 GeV","pixel track p_{T} > 20 GeV",""};

	  Plotter  * p = new Plotter;

	  TH1 * hnF = 0;
	  fF->GetObject(TString::Format("%s_nEvtsForFakesA",prefix.Data()),hnF);
	  if(hnF==0) cout << TString::Format("%s_nEvtsForFakesA",prefix.Data()) << endl;

	  auto addBKGGraph = [&](TString histName, TString title, int inColor) -> Drawing::Drawable1D {
		  TH1 *hf = 0;
		  fF->GetObject(histName,hf);
		  if(hf==0) cout << histName<< endl;
		  hf = (TH1*)hf->Clone();
		  // hf = hf->Rebin(nBins,"",bins);
		  PlotTools::toOverflow(hf);
		  PlotTools::toUnderflow(hf);
		  hf->Scale(1./hnF->GetBinContent(1));

		  TGraphAsymmErrors* gr = new TGraphAsymmErrors(hf);
		  gr->SetLineColor  (StyleInfo::getLineColor(inColor));
		  gr->SetLineWidth  (4);
		  gr->SetLineStyle  (1);
		  gr->SetMarkerStyle(inColor < 2 ? 20 : 21);
		  gr->SetMarkerColor(StyleInfo::getLineColor(inColor));
		  gr->SetMarkerSize (1);
		  Drawing::Drawable1D drawableF("P E 0",title,Drawing::GRAPH,gr,false);
		  drawableF.graphAxisHist = hf;
		  return drawableF;

	  };

	  for(unsigned int iP = 0; pts[iP][0]; ++ iP){
		  auto d1 = addBKGGraph(TString::Format("%s_fake_muon_%s_passME0Muon_%s",prefix.Data(),pts[iP].Data(),var.Data()),ptNs[iP],iP);
		  p->addDrawable(d1);
	  }
	  p->setCanvasSize(1024,726);
	  p->setMargins(0.08,.12,.12,.05);
	  p->setAxisTextSize(0);
	  p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
	  p->setLegendPos(.13,0.75,0.65,0.89);
	  p->addText("2.0 < |#eta| < 2.8",0.15,0.71,0.04);
	  p->setMinMax(0.0,.065);
	  p->setYTitle("<N. of bkg. matches> per BX / 0.1");
	  p->setXTitle(xTitle);
	  p->draw(true,"backgroundRate_eta.pdf");

	  fF->Close();
  }

  void makeTotalBkgPlots(TString trackFile, TString segFile, TString matchFile) {
	  Plotter  * p = new Plotter;


	  //matches
	  if(false){
		  TString digi = "p8s384";
		  TFile * ff = new TFile(matchFile,"READ");
//		  TString pts[] = {"ptgeq2","ptgeq3","ptgeq5","",""};
//		  TString ptNs[] = {"match p_{T} > 2 GeV","match p_{T} > 3 GeV","match p_{T} > 5 GeV",""};
//		  TString pts[] = {"ptgeq5","ptgeq3","ptgeq2","",""};
//		  TString ptNs[] = {"match p_{T} > 5 GeV","match p_{T} > 3 GeV","match p_{T} > 2 GeV",""};
		  TString pts[] = {"ptgeq2","",""};
		  TString ptNs[] = {"match p_{T} > 2 GeV",""};


		  TH1 * hn = 0;
		  ff->GetObject(TString::Format("%s_nEvtsForFakesA",digi.Data()),hn);
		  float nEvents = hn->GetBinContent(1);
		  for(unsigned int iP = 0; pts[iP][0]; ++iP){
			  TH1 * ho;
			  ff->GetObject(TString::Format("%s_fake_muon_%s_passME0Muon_eta",prefix.Data(),pts[iP].Data()),ho);
			  if(ho == 0) continue;
			  TH1 * h = new TH1F(TString::Format("rangeMatch_%s",pts[iP].Data()),"background object p_{T} [GeV]",10,1.9,2.9);
				  for(unsigned int iB = 1; iB <= ho->GetNbinsX(); ++iB){
					  double binLow = ho->GetBinLowEdge(iB);
					  if(binLow < 1.9) continue;
					  if(binLow >= 2.9) continue;
					  h->Fill(ho->GetBinCenter(iB),ho->GetBinContent(iB));
				  }
//			  h->Rebin(10);
			  h->Scale(1./nEvents);
			  p->addHistLine(h,ptNs[iP],StyleInfo::getLineColor(iP),1);
		  }

		  ff->Close();
	  }

	  //segments
	  {
		  TString digi = "p8s384";
		  TFile * ff = new TFile(segFile,"READ");
//		  TString pts[] = {"fake_seg","fake_pt_gt2_seg","fake_pt_gt3_seg","fake_pt_gt5_seg",""};
//		  TString ptNs[] = {"inclusive segments","segment p_{T} > 2 GeV","segment p_{T} > 3 GeV","segment p_{T} > 5 GeV",""};
		  TString pts[] = {"fake_pt_gt5_seg","fake_pt_gt3_seg","fake_pt_gt2_seg","fake_seg",""};
		  TString ptNs[] = {"segment p_{T} > 5 GeV","segment p_{T} > 3 GeV","segment p_{T} > 2 GeV","inclusive segments",""};
//		  TString pts[] = {"fake_pt_gt2_seg","fake_seg",""};
//		  TString ptNs[] = {"segment p_{T} > 2 GeV","inclusive segments",""};


		  TH1 * hn = 0;
		  ff->GetObject(TString::Format("%s___fake_nEvents",digi.Data()),hn);
		  float nEvents = hn->GetBinContent(1);
		  for(unsigned int iP = 0; pts[iP][0]; ++iP){
			  TH1 * ho;
			  ff->GetObject(TString::Format("%s___%s_eta",digi.Data(),pts[iP].Data()),ho);
			  if(ho == 0) continue;
			  TH1 * h = new TH1F(TString::Format("rangeSegs_%s",pts[iP].Data()),"background object p_{T} [GeV]",10,1.9,2.9);
				  for(unsigned int iB = 1; iB <= ho->GetNbinsX(); ++iB){
					  double binLow = ho->GetBinLowEdge(iB);
					  if(binLow < 1.9) continue;
					  if(binLow >= 2.9) continue;
					  h->Fill(ho->GetBinCenter(iB),ho->GetBinContent(iB));
				  }
			  h->Scale(1./nEvents);
			  p->addHistLine(h,ptNs[iP],StyleInfo::getLineColor(iP),9);
		  }

		  ff->Close();
	  }

	  //tracks
	  if(false){
		  TFile * ff = new TFile(trackFile,"READ");
//		  TString pts[] = {"gt0p5","gt2p0","gt3p0","gt5p0","",""};
//		  TString ptNs[] = {"track p_{T} > 0.5 GeV","track p_{T} > 2 GeV","track p_{T} > 3 GeV","track p_{T} > 5 GeV",""};
//		  TString pts[] = {"gt5p0","gt3p0","gt2p0","gt0p5","",""};
//		  TString ptNs[] = {"track p_{T} > 5 GeV","track p_{T} > 3 GeV","track p_{T} > 2 GeV","track p_{T} > 0.5 GeV",""};
		  TString pts[] = {"gt2p0","gt0p5","",""};
		  TString ptNs[] = {"track p_{T} > 2 GeV","track p_{T} > 0.5 GeV",""};


		  TH1 * hn = 0;
		  ff->GetObject(TString::Format("nEvents"),hn);
		  float nEvents = hn->GetBinContent(1);
		  for(unsigned int iP = 0; pts[iP][0]; ++iP){
			  TH1 * ho;
			  ff->GetObject(TString::Format("tracks_%s_abseta",pts[iP].Data()),ho);
			  if(ho == 0) continue;
			  TH1 * h = new TH1F(TString::Format("rangeTracks_%s",pts[iP].Data()),"background object p_{T} [GeV]",10,1.9,2.9);
			  for(unsigned int iB = 1; iB <= ho->GetNbinsX(); ++iB){
				  double binLow = ho->GetBinLowEdge(iB);
				  if(binLow < 1.9) continue;
				  if(binLow >= 2.9) continue;
				  h->Fill(ho->GetBinCenter(iB),ho->GetBinContent(iB));
			  }
			  h->Scale(1./nEvents);
			  p->addHistLine(h,ptNs[iP],StyleInfo::getLineColor(iP),7);
		  }
		  ff->Close();
	  }


	  p->setCanvasSize(1024,726);
	  p->setMargins(0.08,.12,.12,.05);
	  p->setAxisTextSize(0);
	  p->setCMSLumi(0, "14 TeV, 200 PU", "Simulation preliminary",1.2 );
	  p->setLegendPos(.14,0.75,0.92,0.89);
	  p->setLegendNColumns(3);
//	  p->addText("2.0 < |#eta| < 2.8",0.15,0.71,0.03);
	  p->setMinMax(0.001,10000);
	  p->setYTitle("<N. of bkg. objects> per BX / 0.1");
	  p->setXTitle("background object |#eta|");
	  TCanvas * c = p->draw(false,"backgroundRateThreeLevels_eta.pdf");
	  c->SetLogy();
	  c->Print("backgroundRateThreeLevels_eta.pdf");

  }




};

#endif

void MakeTDRPlots(){
	TH1::AddDirectory(false);
	Analyzer * a = new Analyzer();
//	a->makeResolutionPlots("/Users/nmccoll/Dropbox/Work/Projects/ME0/3_13_17_trackMuonMatching/TDRVersion/trackMatchingTree_p8s384_plots.root");
//	a->makeLayersCrossedPlot("/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/simHitAnalyzer.root");
//	a->makeDPhiPlot("/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/simHitTestForTDR.root","/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/digiTestForTDR_p8s384.root");
	a->makeNeutronPlot();
//	a->makeSegmentBkgPlots("/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/segmentAnalyzerForTDR_addPT_p8s384.root");
//	a->makeMatchBkgPlots("/Users/nmccoll/Dropbox/Work/Projects/ME0/3_13_17_trackMuonMatching/TDRVersion/trackMatchingTree_NU_p8s384_POGHP_plots.root");
//	a->makeTotalBkgPlots("/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/trackDensity_ForTDR_std.root",
//			"/Users/nmccoll/Dropbox/Work/Projects/ME0/2_15_17_updatedTruthPlots/segmentAnalyzerForTDR_addPT_p8s384.root",
//			"/Users/nmccoll/Dropbox/Work/Projects/ME0/3_13_17_trackMuonMatching/TDRVersion/trackMatchingTree_NU_p8s384_POGHP_plots.root");

}
