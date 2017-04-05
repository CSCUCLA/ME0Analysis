//Find cuts
{
  // TString filename = "trackMatchingTree_p6s512_plots.root";
  // TString prefix = "p6s512";
    // TString filename = "trackMatchingTree_p8s384_plots.root";
        // TString prefix = "p8s384";
        
        // TString filename = "trackMatchingTree_p8s256_plots.root";
            // TString prefix = "p8s256";
            
            TString filename = "trackMatchingTree_p8s384_plots.root";
                TString prefix = "p8s384";
  TFile * f = new TFile(filename,"READ");
  TString extr = "";
  // TString extr = "absDPhilt0p1_absPhiEtaLt0p1_";
  TString types [] = {"signal","fake",""};
  TString typeNs [] = {"signal","fake",""};
  TString pts [] = {"pt_eq2to3","pt_eq3to5","pt_eq5to10","pt_eq10to20","pt_geq20",""};
  TString ptNs [] = {"p_{T} 2-3 GeV","p_{T} 3-5 GeV","p_{T} 5-10 GeV","p_{T} 10-20 GeV","p_{T} >20 GeV",""};
  // TString vars [] = {"resid_x","pull_x","resid_y","pull_y","resid_dx","pull_dx","resid_dy","pull_dy","resid_eta","resid_phi","resid_deta","resid_dphi",""};
    // TString vars [] = {"resid_x","pull_x","resid_y","pull_y","resid_eta","resid_phi","resid_deta","resid_dphi",""};
        // TString vars [] = {"resid_eta","resid_phi","resid_deta","resid_dphi","resid_dphiS","resid_x","pull_x","resid_xs",""};
                TString vars [] = {"resid_eta","resid_phi","resid_dphi",""};
                
                                // TString vars [] = {"resid_phi","resid_phiS1","resid_phiS2","resid_dphi","resid_dphiS1","resid_dphiS2",""};
  
                  // TString vars [] = {"pull_x","pull_y","pull_dx","pull_dy",""};
        // TString vars [] = {"resid_seg_dphi",""};
  
  double values[] = {.999,.997,.995,.994,-1};
  
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * pDist = new Plotter;
    Plotter * pInt  = new Plotter;    
      std::vector<int> binS;    
    for(unsigned int iT = 0; types[iT][0]; ++iT){
          for(unsigned int iP = 0; pts[iP][0]; ++iP){
            TH1 * hd = 0;
            f->GetObject(TString::Format("%s_%s_%s_%s%s",prefix.Data(), types[iT].Data(), pts[iP].Data(), extr.Data(), vars[iV].Data()),hd);
            TH1 * hi = 0;
            f->GetObject(TString::Format("%s_%s_%s_%sabs_%s",prefix.Data(), types[iT].Data(), pts[iP].Data(), extr.Data(), vars[iV].Data()),hi);
            if(hd == 0 || hi == 0){
              std::cout << TString::Format("%s_%s_%s_%sabs_%s",prefix.Data(), types[iT].Data(), pts[iP].Data(), extr.Data(), vars[iV].Data())<<endl;
              std::cout << TString::Format("%s_%s_%s_%s%s",prefix.Data(), types[iT].Data(), pts[iP].Data(), extr.Data(), vars[iV].Data())<<endl;
              continue;
            }
            pDist->addHistLine(hd, TString::Format("%s, %s",typeNs[iT].Data(),ptNs[iP].Data()), StyleInfo::getLineColor(iP),iT + 1 );
            PlotTools::toUnderflow(hi);
            PlotTools::toOverflow(hi);
            PlotTools::normalize(hi);
            TH1 * hi2 = (TH1*)hi->Clone();
            for(unsigned int iB = 1; iB <= hi->GetNbinsX(); ++iB){
              hi2->SetBinContent(iB,hi->Integral(0,iB));
            }
            hi2->SetBinContent(0,0);
            hi2->SetBinContent(hi2->GetNbinsX()+1,0);
            
            
            if(iP == 3 && iT ==0 ){
              cout << vars[iV]<< "\t S \t";
               for(unsigned int iA = 0; values[iA] >= 0; ++iA){
                 for(unsigned int iB = 1; iB <= hi2->GetNbinsX(); ++iB){
                    if(hi2->GetBinContent(iB) < values[iA]  &&  hi2->GetBinContent(iB + 1)  >=  values[iA]  ){
                      binS.push_back(iB);
                      cout << hi2->GetBinLowEdge(iB) + hi2->GetBinWidth(iB) <<"\t";
                    }
                 }
               }
               cout << endl;
             }
             
             if(iP == 3 && iT ==1 ){
               cout << vars[iV]<< "\t B \t";
                for(unsigned int iA = 0; iA < binS.size(); ++iA){
                  cout << hi2->GetBinContent(binS[iA]) << "\t";
                }
                cout << endl;
              }
              
            pInt->addHistLine(hi2, TString::Format("%s, %s",typeNs[iT].Data(),ptNs[iP].Data()), StyleInfo::getLineColor(iP),iT + 1 );
          }
    }
    pDist->normalize();
    pDist->draw(false,vars[iV]);
    pInt->setMinMax(0,1.2);
    pInt->draw(false,vars[iV] +"_int");
  }
}

//Effplot
{
  TString filename = "trackMatchingTree_p8s384_plots.root";
  TString prefix = "p8s384";
  TFile * f = new TFile(filename,"READ");

  TString extr = "";
  
    // TString efftypes [] = {"incl","goodTrack","goodSegment","goodTrackAndSegment","goodME0Muon",""};
        // TString efftypes [] = {"incl","goodSegment"  ,"goodSegmentAndGoodPT","goodME0Muon",""};
        // TString efftypeNs[] = {"incl","segment reco.","+ gen-matched pixel track (20% p_{T} match)","+ exclusive track-segment match",""};
  
  // TString efftypes [] = {"incl","goodME0Muon","goodSegmentAndGoodPT","goodSegment"  ,""};
  // TString efftypeNs[] = {"incl","+ exclusive track-segment match","+ gen-matched pixel track (20% p_{T} match)","segment reco.",""};
  TString efftypes [] = {"goodTrackPT","goodME0Muon","goodSegmentAndGoodPT",""  ,""};
  TString efftypeNs[] = {"incl","ME0Muon","ME0Segment",""};
  
  double bins[] = {0,1,2,3,4,5,10,15,20,25,30};
  int nBins = 10;
    p = new Plotter;
        for(unsigned int iT = 0; efftypes[iT][0]; ++iT){
          TH1 * hd = 0;
          f->GetObject(TString::Format("%s_real_muon_%s_pt",prefix.Data(), efftypes[iT].Data()),hd);
          hd = hd->Rebin(nBins,"",bins);
          p->addHist(hd,efftypeNs[iT]);
        }
        

        // p->rebin(nBins,bins);
        p->drawRatio(0,"stack",true,false);
    
  
}

//Effplot++ : Add fakes for Anna
{
  TString filename = "trackMatchingTree_p8s384_plots.root";
  TString filenameF = "trackMatchingTree_NU_p8s384_plots.root";
  TString prefix = "p8s384";
  TFile * f = new TFile(filename,"READ");
  TFile * fF = new TFile(filenameF,"READ");

  TString extr = "";
  

  TString efftypes [] = {"passME0Muon","passSegment",""  ,""};
  TString efftypeNs[] = {"ME0Muon reconstruction efficiency","ME0Segment reconstruction efficiency",""};
  
  // double bins[] = {0,1,2,3,4,5,10,15,20,25,30};
  // int nBins = 10;
  
  // double bins[] = {2,3,4,5,7,10,15,20,25,30};
  // int nBins = 9;
  
    p = new Plotter;
    
    //eff
    TH1 * hd = 0;
    f->GetObject(TString::Format("%s_real_muon_passTrack_pt",prefix.Data()),hd);
    // hd = hd->Rebin(nBins,"",bins);
    PlotTools::toOverflow(hd);
    PlotTools::toUnderflow(hd);
    
        for(unsigned int iT = 0; efftypes[iT][0]; ++iT){
          TH1 * hn = 0;
          f->GetObject(TString::Format("%s_real_muon_%s_pt",prefix.Data(), efftypes[iT].Data()),hn);
          // hn = hn->Rebin(nBins,"",bins);
          PlotTools::toOverflow(hn);
          PlotTools::toUnderflow(hn);

          hn->SetLineColor(StyleInfo::getLineColor(iT));
          hn->SetLineWidth(3);
          hn->SetMarkerColor(StyleInfo::getLineColor(iT));
          hn->SetMarkerStyle(20);
          hn->SetMarkerSize(1);

          Drawing::Drawable1D drawable = Drawing::makeRatio(hn,hd,efftypeNs[iT],"E1X0 P ",true);
          p->addDrawable(drawable);
        }
        
        // FR
        TH1 * hnF = 0;
        fF->GetObject(TString::Format("%s_nEvtsForFakesA",prefix.Data()),hnF);
        if(hnF==0) cout << TString::Format("%s_nEvtsForFakesA",prefix.Data()) << endl;
        TH1 *hf = 0;
        fF->GetObject(TString::Format("%s_fake_muon_passME0Muon_pt",prefix.Data()),hf);
        if(hf==0) cout << TString::Format("%s_fake_muon_passME0Muon_pt",prefix.Data()) << endl;
        // hf = hf->Rebin(nBins,"",bins);
        PlotTools::toOverflow(hf);
        PlotTools::toUnderflow(hf);
        hf->Scale(1./hnF->GetBinContent(1));
        
        p->addStackHist(hf,"# of background muons per event");
        
        //         TGraphAsymmErrors* gr = new TGraphAsymmErrors(hf);
        // gr->SetLineColor  (StyleInfo::getLineColor(2));
        // gr->SetLineWidth  (3);
        // gr->SetLineStyle  (1);
        // gr->SetMarkerStyle(20);
        // gr->SetMarkerColor(StyleInfo::getLineColor(2));
        // gr->SetMarkerSize (1);
        //         Drawing::Drawable1D drawableF("P 0","Fakes",Drawing::GRAPH,gr,false);
        //         drawableF.graphAxisHist = hf;
        //         p->addDrawable(drawableF);
        p->setMinMax(0.0,1.0);
        p->draw(false,"tot");
        
    
  
}


{
  TString fN[] ={"trackMatchingTree_p6s512_plots.root","trackMatchingTree_p8s384_plots.root",""};
  TString fL[] ={"6x512","8x384",""};
  TString fpre[] ={"p6s512","p8s384",""};
  std::vector<TFile*> fls;
  for(unsigned int iF = 0; fN[iF][0]; ++iF){
    fls.push_back(new TFile(fN[iF],"READ"));
  }
  TString extr = "";     
  TString den = "noCuts";
  TString num = "eta_lt0p08_phi_lt0p05_dphi_lt0p0065";

  Plotter * pDphR = new Plotter;
  
  for(unsigned int iF = 0; fN[iF][0]; ++iF ){  
    TH1 * hd = 0;
    fls[iF]->GetObject(TString::Format("%s_real_muon_%s_pt",fpre[iF].Data(), den.Data()),hd);    
    TH1 * hn = 0;
    fls[iF]->GetObject(TString::Format("%s_real_muon_%s_pt",fpre[iF].Data(), num.Data()),hn);
    if(!hd || !hn) continue;
        PlotTools::toUnderflow(hd);PlotTools::toOverflow(hd);
        PlotTools::toUnderflow(hn);PlotTools::toOverflow(hn);
    hn->Divide(hd);
    pDphR->addHistLine(hn,fL[iF]);    
  }
  pDphR->draw(false,"fakeDPhi");
    
  
}


{
  TString filename = "trackMatchingTree_p6s512_plots.root";
  TFile * f = new TFile(filename,"READ");
  TString prefix = "p6s512";
  TString extr = "";
  
  double maxDPhi []  = {0,0.009,0.0065,0.0055,0.004,-1};
  double maxPhi  []  = {0,0.055,0.045,0.03,-1};
  double maxEta  []  = {0,0.08,0.07,0.06,-1};
  
  Plotter * pDphR = new Plotter;
  Plotter * pDphF = new Plotter;
  for(unsigned int iD = 0; maxDPhi[iD] >= 0; ++iD ){  
	  TString dPhiStr= TString::Format("%.4f",maxDPhi[iD]).ReplaceAll(".","p");
    TH1 * hm = 0;
    f->GetObject(TString::Format("%s_real_muon_dPhi_%s_phi_0p000_eta_0p00_pt",prefix.Data(), dPhiStr.Data()),hm);
    if(hm) pDphR->addHistLine(hm,dPhiStr);
    TH1 * hf = 0;
    f->GetObject(TString::Format("%s_fake_muon_dPhi_%s_phi_0p000_eta_0p00_pt",prefix.Data(), dPhiStr.Data()),hf);
    if(hf == 0) continue;
        hf = (TH1*)hf->Clone();
    hf->Scale(1./10000.);
    pDphF->addHistLine(hf,dPhiStr);    
    cout << dPhiStr.Data() <<" -> "<< hf->Integral(0,-1)<<endl;
    
  }
  pDphR->drawRatio(0,"stack",true,false,"realDPhi");
  pDphF->draw(false,"fakeDPhi");
  
  
  Plotter * pPhR = new Plotter;
  Plotter * pPhF = new Plotter;
  for(unsigned int iD = 0; maxPhi[iD] >= 0; ++iD ){  
	  TString dPhiStr= TString::Format("%.3f",maxPhi[iD]).ReplaceAll(".","p");
    TH1 * hm = 0;
    f->GetObject(TString::Format("%s_real_muon_dPhi_0p0000_phi_%s_eta_0p00_pt",prefix.Data(), dPhiStr.Data()),hm);
    if(hm) pPhR->addHistLine(hm,dPhiStr);
    TH1 * hf = 0;
    f->GetObject(TString::Format("%s_fake_muon_dPhi_0p0000_phi_%s_eta_0p00_pt",prefix.Data(), dPhiStr.Data()),hf);
    if(hf == 0) continue;
        hf = (TH1*)hf->Clone();
    hf->Scale(1./10000.);
    pPhF->addHistLine(hf,dPhiStr);    
    cout << dPhiStr.Data() <<" -> "<< hf->Integral(0,-1)<<endl;
    
  }
  pPhR->drawRatio(0,"stack",true,false,"realPhi");
  pPhF->draw(false,"fakePhi");
  
  Plotter * pEtR = new Plotter;
  Plotter * pEtF = new Plotter;
  for(unsigned int iD = 0; maxEta[iD] >= 0; ++iD ){  
	  TString dPhiStr= TString::Format("%.2f",maxEta[iD]).ReplaceAll(".","p");
    TH1 * hm = 0;
    f->GetObject(TString::Format("%s_real_muon_dPhi_0p0000_phi_0p000_eta_%s_pt",prefix.Data(), dPhiStr.Data()),hm);
    if(hm) pEtR->addHistLine(hm,dPhiStr);
    else cout << TString::Format("%s_real_muon_dPhi_0p0000_phi_0p000_eta_%s_pt",prefix.Data(), dPhiStr.Data()) <<endl;
    TH1 * hf = 0;
    f->GetObject(TString::Format("%s_fake_muon_dPhi_0p0000_phi_0p000_eta_%s_pt",prefix.Data(), dPhiStr.Data()),hf);
    if(hf == 0) continue;
    hf = (TH1*)hf->Clone();
    hf->Scale(1./10000.);
    pEtF->addHistLine(hf,dPhiStr);    
    cout << dPhiStr.Data() <<" -> "<< hf->Integral(0,-1)<<endl;
    
  }
  pEtR->drawRatio(0,"stack",true,false,"realEta");
  pEtF->draw(false,"fakeEta");
  

  
  
}


{
  TString filename = "trackMatchingTree_NU_p8s384_plots.root";
  TString prefix = "p8s384";
  // TString filename = "trackMatchingTree_NU_p6s512_plots.root";
  // TString prefix = "p6s512";
  
  // TString filename = "trackMatchingTree_NU_p6s384_plots.root";
  // TString prefix = "p6s384";
  // TString filename = "trackMatchingTree_NU_p8s256_plots.root";
  // TString prefix = "p8s256";
  
  TFile * f = new TFile(filename,"READ");

  TString extr = "";
  
  TH1 * hn = 0;
  f->GetObject(TString::Format("%s_nEvtsForFakes",prefix.Data()),hn);
  
  TString cuts []  = {"noCuts","dEta_lt0p08","eta_lt0p08_phi_lt0p05","eta_lt0p08_phi_lt0p05_dphi_lt0p0065","eta_lt0p08_phi_lt0p05_dphi_lt0p0065_tP","eta_lt0p08_phi_lt0p05_dphi_lt0p0065_tPtDP",""};
  TString cutNs []  = {"incl.","match #eta","match #eta, #phi","match #eta, #phi, #Delta#phi","match #eta, #phi, tight #phi","match #eta, #phi, tight #phi, #Delta#phi",""};
  double ptvals[] = {1,2,3,5,10,-1};
  Plotter * pDphR = new Plotter;
  Plotter * pDphF = new Plotter;
  for(unsigned int iD = 0; cuts[iD][0]; ++iD ){  
    TH1 * hm = 0;
    f->GetObject(TString::Format("%s_real_muon_%s_pt",prefix.Data(), cuts[iD].Data()),hm);
    if(hm) pDphR->addHistLine(hm,cutNs[iD]);
    TH1 * hf = 0;
    f->GetObject(TString::Format("%s_fake_muon_%s_pt",prefix.Data(), cuts[iD].Data()),hf);
    if(hf == 0) continue;
    hf = (TH1*)hf->Clone();
    hf->Scale(1./hn->GetBinContent(1));
    hf->SetYTitle("<# of BKG ME0Muons>/event");
    pDphF->addHistLine(hf,cutNs[iD]);    
    
    cout << cutNs[iD].Data() <<"\t";
    for(unsigned int iP=0; ptvals[iP] >= 0; ++iP)
      cout << hf->Integral(hf->FindFixBin(ptvals[iP]),-1) <<"\t";
    cout << endl;
    
  }
  // pDphR->drawRatio(0,"stack",true,false,"realDPhi");
  pDphF->draw(false,"fakeDPhi");
    
}




//compare all fakes
{
  TString filename = "trackMatchingTree_NU_plots.root";
  TFile * f = new TFile(filename,"READ");
  

  
  TString prefix[] = {"p8s384","p6s512","p6s384","p8s256",""};
  TString prefixNs[] = {"(8x384)","(6x512)","(6x384)","(8x256)",""};
  TString extr = "";
  TString vars [] = {"eta","eta_phi","eta_phi_dPhi","eta_phi_dPhi_S1","eta_phi_dPhi_S2","",""};

  // p8s384_fake_muon_eta_phi_dPhi_S2_Int
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * pDist = new Plotter;
    cout << vars[iV] <<endl;
    for(unsigned int iP = 0; prefix[iP][0]; ++iP){      
      TH1 * hm = 0;
      f->GetObject(TString::Format("%s_fake_muon_%s_Int",prefix[iP].Data(), vars[iV].Data()),hm);
      TH1 * hi = 0;
      f->GetObject(TString::Format("%s_nEvtsForFakes",prefix[iP].Data()),hi);
      cout <<prefixNs[iP]<<"\t";
      for(unsigned int iB = 1; iB <= hm->GetNbinsX(); ++iB)
        cout << TString::Format("%.4f+/-%.4f",hm->GetBinContent(iB)/hi->GetBinContent(1),hm->GetBinError(iB)/hi->GetBinContent(1))<<"\t";
      cout << endl;
      
      hm->Scale(1./hi->GetBinContent(1));
      pDist->addHist(hm, prefixNs[iP]);
    }
    pDist->draw(false,vars[iV]);
      }
      
    }
    
    //COmpare all effs
    
    {
      TString filename = "trackMatchingTree_plots.root";
  
      TString prefix[] = {"p8s384","p6s512","p6s384","p8s256",""};
      TString prefixNs[] = {"(8x384)","(6x512)","(6x384)","(8x256)",""};
      TFile * f = new TFile(filename,"READ");


      TString num = "goodME0Muon";
      TString den = "goodTrackPT";

  
      double bins[] = {0,1,2,3,4,5,10,15,20,25,30};
      int nBins = 10;
        p = new Plotter;
    for(unsigned int iP = 0; prefix[iP][0]; ++iP){   
              TH1 * hd = 0;
              f->GetObject(TString::Format("%s_real_muon_%s_pt",prefix[iP].Data(), den.Data()),hd);
              if(hd == 0) {
cout <<                TString::Format("%s_real_muon_%s_pt",prefix[iP].Data(), den.Data())<<endl;
                continue;
              }
              PlotTools::toUnderflow(hd);
              PlotTools::toOverflow(hd);
              hd = hd->Rebin(nBins,"",bins);
              TH1 * hn = 0;
              f->GetObject(TString::Format("%s_real_muon_%s_pt",prefix[iP].Data(), num.Data()),hn);
              if(hn == 0) {
cout <<                TString::Format("%s_real_muon_%s_pt",prefix[iP].Data(), num.Data())<<endl;
                continue;
              }
              PlotTools::toUnderflow(hn);
              PlotTools::toOverflow(hn);
              hn = hn->Rebin(nBins,"",bins);
              hn->Divide(hn,hd,1,1,"b");
              p->addHist(hn,prefixNs[iP]);
            }
        

            // p->rebin(nBins,bins);
            p->draw(false,"mon");
    //
  
    }
