{

	Plotter * p = new Plotter();
	p->addHistLine(gen_eta_PToPT,"gen_eta_PToPT");
	p->addHistLine(gen_theta_PToPT,"gen_theta_PToPT");
	p->addHistLine(reco_eta_PToPT,"reco_eta_PToPT");
	p->addHistLine(reco_theta_PToPT,"reco_theta_PToPT");
	
	p->normalize();
	p->draw();
}


//pt resolution
{
	TFile * f = new TFile("plots_segmentTree_tenMu_p8s384.root","read");
	TString types[] = {"reco","gen",""};
	// TString etas[] = {"eta2to2p2","eta2p2to2p4","eta2p2to2p6","eta2p6to2p8",""};
		TString etas[] = {"eta2to2p2","eta2p6to2p8",""};
	TString pts[] = {"pt1to5","pt5to15","pt15to30",""};
	for(unsigned int iT = 0; types[iT][0]; ++iT){
	Plotter * p = new Plotter();
		for(unsigned int iP = 0; pts[iP][0]; ++iP){
		for(unsigned int iE = 0; etas[iE][0]; ++iE){
		TH1 * h = 0;
		f->GetObject(TString::Format("%s_%s_%s_PToPT",types[iT].Data(),pts[iP].Data(),etas[iE].Data()),h);
		if(h == 0) continue;
    int style = 1;
    if(iE) style = 9;
		p->addHistLine(h,TString::Format("%s, %s",pts[iP].Data(),etas[iE].Data()), StyleInfo::getLineColor(iP),style  );
    cout << TString::Format("%s_%s_%s_PToPT",types[iT].Data(),pts[iP].Data(),etas[iE].Data())<<" " << h->GetMean() <<" "<< h->GetStdDev() <<endl;
		}
		}
    // p->rebin(2);
    p->normalize();
    
		p->draw(false,types[iT]);
	}
}



//compare resolutions
{
	TFile * f1 = new TFile("plots_segmentTree_tenMu_p8s384.root","read");
  TFile * f2 = new TFile("plots_segmentTree_tenMu_p8s192.root","read");
  TString typeNames[] = {"sim-hit","384 strips","192 strips",""};
    std::vector<TFile*> files = {f2,f1,f2};
    TString types[] = {"gen","reco","reco",""};
  
  TString etas[] = {"eta2to2p2","eta2p6to2p8",""};
  TString pts[] = {"pt1to5","pt15to30",""};
  TString ptNs[] = {"p_{T} 2-5 GeV","p_{T} 15-30 GeV",""};

  
  		for(unsigned int iE = 0; etas[iE][0]; ++iE){
                   for(unsigned int iP = 0; pts[iP][0]; ++iP){
       Plotter * p = new Plotter();
       for(unsigned int iT = 0; typeNames[iT][0]; ++iT){      

         		TH1 * h = 0;
            files[iT]->cd();
         		files[iT]->GetObject(TString::Format("%s_%s_%s_PToPT",types[iT].Data(),pts[iP].Data(),etas[iE].Data()),h);
            if(h == 0) continue;
            int style = 1;
            if(iP) style = 9;
            // p->addHistLine(h,TString::Format("%s, %s",typeNames[iT].Data(),ptNs[iP].Data()), StyleInfo::getLineColor(iT),style  );  
            p->addHistLine(h,TString::Format("%s",typeNames[iT].Data()));  
            cout << TString::Format("%s_%s_%s_PToPT",types[iT].Data(),pts[iP].Data(),etas[iE].Data())<<" "<<files[iT]<<" " << h->GetMean() <<" "<< h->GetStdDev() <<endl;
                                
           }
           p->normalize();
           p->rebin(4);
           p->setYTitle("arbitrary units");
       		p->draw(false,etas[iE]+pts[iP]);
       }

      }
      
    }
    
    
//compare resolutions L1 /L2
{
	TFile * f1 = new TFile("plots_segmentTree_tenMu_p8s384.root","read");
  TFile * f2 = new TFile("plots_segmentTree_tenMu_p8s192.root","read");

  std::vector<TFile*> files = {f1,f2};
  TString filesNames[] = {"level-2","level-1",""};

  TString types[] = {"gen","reco","",""};
  TString genName = "sim-hits";

  TString pts[] = {"pt1to5","pt15to18",""};
  TString ptNs[] = {"p_{T} 2-5 GeV","p_{T} 15-18 GeV",""};

  TString etas[] = {"eta2to2p2","eta2p6to2p8",""};
  TString etaNs[] = {"2.0 < |#eta| < 2.2","2.6 < |#eta| < 2.8",""};

  for(unsigned int iF = 0; iF < files.size(); ++iF){
      for(unsigned int iP = 0; pts[iP][0]; ++iP){
        Plotter * p = new Plotter();
                  for(unsigned int iE = 0; etas[iE][0]; ++iE){
        for(unsigned int iT = 0; types[iT][0]; ++iT){

         		TH1 * h = 0;
            files[iF]->cd();
         		files[iF]->GetObject(TString::Format("%s_%s_%s_PToPT",types[iT].Data(),pts[iP].Data(),etas[iE].Data()),h);   
            if(h == 0) continue;
            int style = 1;
            if(iT == 0) style = 9;
            TString type = iT == 0 ? genName :filesNames[iF];
        		p->addHistLine(h,TString::Format("%s, %s",type.Data(),etaNs[iE].Data()));//, StyleInfo::getLineColor(iE),style  );
          }          
        }        
        p->normalize();
        p->rebin(4);
        p->setYTitle("arbitrary units");
    		p->draw(false,filesNames[iF]+pts[iP]);
      }              
  }

}
  
//integral pT 

//ZMM fake pt res
{
	TFile * f = new TFile("plots_segmentTree_ZMM.root","read");
  TString pres[] = {"p8s192","p8s384",""};
  TString prens[] = {"192 strips","384 strips",""};
  // TString vars[] = {"pt",""};
  // TString pre2 = "zmm_twoInner";
  // int rebin = -1;
  // TString ytitle = "<N. of bkg. segments> per BX / 1 GeV";
    // TString pre2 = "zmm_oneInnerOneME0";
    // TString vars[] = {"bkg_invmass","passChargeAndPT_bkg_invmass",""};
    // int rebin = 10;
    // TString ytitle = "<N. of bkg. segments> per BX / 10 GeV";
  
  // TString pre2 = "zmm_twoInner";
  // TString vars[] = {"diMu_bkg_invmass","diMu_passChargeAndPT_bkg_invmass",""};
  // int rebin = 10;
  // TString ytitle = "<N. of bkg. segments> per BX / 10 GeV";
  
  TString pre2 = "tau3mu";
  TString vars[] = {"twoInner_pt1_bkg_invmass","twoInner_pt2_bkg_invmass","twoInner_pt1_bkg_min_invmass","twoInner_pt2_bkg_min_invmass","oneInner_pt1_bkg_invmass","oneInner_pt2_bkg_invmass","oneInner_pt1_bkg_min_invmass","oneInner_pt2_bkg_min_invmass",""};
  int rebin = 1;
  TString ytitle = "<N. of objects> per BX / 0.5 GeV";
  
   
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * p = new Plotter();
    for(unsigned int iP = 0; pres[iP][0];++iP){
      TH1 * hn = 0;
      f->GetObject(TString::Format("%s_%s_nFakeEvents",pres[iP].Data(),pre2.Data()),hn);
      if(hn == 0) continue;
      double scale = 1./hn->GetBinContent(1);
      TH1 * h = 0;
      f->GetObject(TString::Format("%s_%s_%s",pres[iP].Data(),pre2.Data(),vars[iV].Data()),h);
      if(h == 0) continue;
      h->Scale(scale);
      p->addHistLine(h,prens[iP]);            
    }
    p->setYTitle(ytitle);
    // p->addText("2.4 < |#eta| < 2.8");
    // p->addText("PU 200");
    if(rebin > 0) p->rebin(rebin);
    p->draw(false,vars[iV]);

  }
  

}    

//Min for tau 
{
	TFile * f = new TFile("plots_segmentTree_ZMM.root","read");
  TString pres[] = {"p8s192","p8s384",""};
  TString prens[] = {"level-1","level-2",""};

  TString pre2 = "tau3mu";
  TString vars[] = {"twoInner_pt1_bkg_min_invmass","twoInner_pt2_bkg_min_invmass","oneInner_pt1_bkg_min_invmass","oneInner_pt2_bkg_min_invmass",""};
  int rebin = 1;
  TString ytitle = "fraction of events w/ M[#tau] < X" ;
  bool integrateDown = true;
  
  // bool integrateDown = true;
  // TString pre2 = "zmm_twoInner";
  // TString vars[] = {"pt",""};
  // int rebin = 1;
  // TString ytitle = "<N. of bkg. ME0Segments> per BX (p_{T} > X)" ;
  // bool integrateDown = false;
  
  // TString pre2 = "zmm_twoInner";
  // TString vars[] = {"max_pt",""};
  // int rebin = 1;
  // TString ytitle = "fraction of BX's w/ #geq1 bkg. ME0Segment (p_{T}>X)" ;
  // bool integrateDown = false;
  
  // TString pre2 = "zmm_twoInner";
  // TString vars[] = {"secmax_pt",""};
  // int rebin = 1;
  // TString ytitle = "fraction of BX's w/ #geq2 bkg. ME0Segment (p_{T}>X)" ;
  // bool integrateDown = false;
  
   
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * p = new Plotter();
    for(unsigned int iP = 0; pres[iP][0];++iP){
      TH1 * hn = 0;
      f->GetObject(TString::Format("%s_%s_nFakeEvents",pres[iP].Data(),pre2.Data()),hn);
      if(hn == 0) continue;
      double scale = 1./hn->GetBinContent(1);
      TH1 * h = 0;
      f->GetObject(TString::Format("%s_%s_%s",pres[iP].Data(),pre2.Data(),vars[iV].Data()),h);
      if(h == 0) continue;
      h->Scale(scale);
      TH1* hi = (TH1*)h->Clone();
      for(unsigned int iB = 1; iB <= h->GetNbinsX(); ++iB){
        if(integrateDown)
          hi->SetBinContent(iB,h->Integral(0,iB));
        else
          hi->SetBinContent(iB,h->Integral(iB,-1));
      }
      hi->SetBinContent(h->GetNbinsX()+1,0);
      p->addHistLine(hi,prens[iP]);            
    }
    p->setYTitle(ytitle);
    // p->addText("2.4 < |#eta| < 2.8");
    // p->addText("PU 200");
    if(rebin > 0) p->rebin(rebin);
    p->draw(false,vars[iV]);

  }
  

}    


//Print
{
  auto print =[](const TH1*h){
    for(unsigned int iB = 1; iB <= h->GetNbinsX(); ++iB)
      cout << h->GetBinContent(iB) <<"\t+/-\t"<< h->GetBinError(iB) <<"\t";
    cout << endl;    
  };
  p8s192_zmm_twoME0_eff->Scale(1./p8s192_zmm_twoME0_eff->GetBinContent(1.));
  print(p8s192_zmm_twoME0_eff);
  p8s384_zmm_twoME0_eff->Scale(1./p8s384_zmm_twoME0_eff->GetBinContent(1.));
  print(p8s384_zmm_twoME0_eff);
  
  p8s192_zmm_twoInner_diMu_bkgMult->Scale(1./p8s192_zmm_twoInner_nFakeEvents->GetBinContent(1.));
  print(p8s192_zmm_twoInner_diMu_bkgMult);
  p8s384_zmm_twoInner_diMu_bkgMult->Scale(1./p8s384_zmm_twoInner_nFakeEvents->GetBinContent(1.));
  print(p8s384_zmm_twoInner_diMu_bkgMult);


}


//single mu zmumu mass
{
	TFile * f = new TFile("plots_segmentTree_ZMM.root","read");
  TString pres[] = {"p8s192","p8s384","p8s384","",""};
  TString prens[] = {"192 strips","384 strips","sim-hit",""};
  TString prees[] = {"","","simhit",""};
  // TString vars[] = {"pt","invmass",""};
  // TString pre2 = "zmm_oneInnerOneME0";
  
  TString vars[] = {"min_pt","max_pt","invmass",""};
  TString pre2 = "zmm_twoME0";
   
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * p = new Plotter();
    for(unsigned int iP = 0; pres[iP][0];++iP){
      TH1 * h = 0;
      TString var = vars[iV];
      if(iP > 1) var = prees[iP]+"_" + var;
      f->GetObject(TString::Format("%s_%s_%s",pres[iP].Data(),pre2.Data(),var.Data()),h);
      if(h == 0) continue;
      p->addHistLine(h,prens[iP]);            
    }
    p->setYTitle("arbitrary units");
p->normalize();
p->rebin(4);
    p->draw(false,vars[iV]);

  }
  
}
    
{
	TFile * f = new TFile("plots_segmentTree_ZMM_p8s384.root","read");
	TString  prefix = "p8s384_fs";
	// TString etas[] = {"eta2to2p2","eta2p2to2p4","eta2p2to2p6","eta2p6to2p8",""};
		TString types[] = {"singleMu5_me0Mugt2","singleMu20_me0Mugt2","singleMu5_me0Mugt5","singleMu20_me0Mugt5", ""};
    		TString typeNs[] = {"p_{T} 5 GeV |#eta| < 2.4 muon + p_{T} > 2 GeV ME0Segment","p_{T} 20 GeV |#eta| < 2.4 muon + p_{T} > 2 GeV ME0Segment",
        "p_{T} 5 GeV |#eta| < 2.4 muon + p_{T} > 5 GeV ME0Segment","p_{T} 20 GeV |#eta| < 2.4 muon + p_{T} > 5 GeV ME0Segment" ""};
        
	TString vars[] = {"invmass",""};
  
  TH1 * hn = 0;
  f->GetObject(TString::Format("%s_nEvents",prefix.Data()),hn);
  double nEvents = hn->GetBinContent(1);
  
		for(unsigned int iV = 0; vars[iV][0]; ++iV){
	Plotter * p = new Plotter();
  	for(unsigned int iT = 0; types[iT][0]; ++iT){
		TH1 * h = 0;
		f->GetObject(TString::Format("%s_%s_%s",prefix.Data(),types[iT].Data(),vars[iV].Data()),h);
		if(h == 0) continue;
    h->Scale(1./nEvents);
		p->addHistLine(h,typeNs[iT]);
		}

    p->rebin(5);
		p->draw(false,vars[iV]);
	}
}