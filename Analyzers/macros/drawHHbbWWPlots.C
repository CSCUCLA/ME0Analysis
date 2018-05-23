rr -b -q 'AnalyzeCutTree.C+("bkg_tree.root","NONE","bkg_plots.root")'                &
rr -b -q 'AnalyzeCutTree.C+("signal_1000_tree.root","1000","sig_1000_plots.root")'   &
rr -b -q 'AnalyzeCutTree.C+("signal_1600_tree.root","1600","sig_1600_plots.root")'   &
hadd -f compiled.root *_plots.root &
  
rr -b -q 'AnalyzeCutTree.C+("bkg_tree_lepFJSkim.root","NONE","bkg_plots_postHbb.root")'                &
rr -b -q 'AnalyzeCutTree.C+("signal_1000_tree.root","1000","sig_1000_plots_postHbb.root")'   &
rr -b -q 'AnalyzeCutTree.C+("signal_1600_tree.root","1600","sig_1600_plots_postHbb.root")'   &  
  hadd -f compiled_postHbb.root *_plots_postHbb.root &

{
  // TFile * f = new TFile("compiled_postHbb.root","read");
  TFile * f = new TFile("compiled.root","read");
  // TFile * f = new TFile("compiled.root","read");
  TString samples[] = {"t","wlnu","ttbar","signal_1000","signal_1600",""};
  TString sampleNs[] = {"single t","W#rightarrowl#nu","t#bar{t}","M(1000)","M(1600)"};
  int nBkg = 3;
  TString prefix = "";
  
 
  auto makePlots  = [&](TString * vars, TString name, TString prefix, int rebin, bool doNorm = true ){
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
    Plotter * p = new Plotter();
    for(unsigned int iS = 0; samples[iS][0]; ++iS){
      TH1 *h = 0;
      f->GetObject(TString::Format("%s_%s%s",samples[iS].Data(),prefix.Data(),vars[iV].Data()),h);
      if(h == 0) continue;
      if(iS < nBkg ) p->addStackHist(h,sampleNs[iS]);
      else{
        if(doNorm){
          float stackNorm = p->getStackIntegral();                  
          if(stackNorm) h->Scale(stackNorm/h->Integral(0,-1)); else PlotTools::normalize(h);
        } else if(p->getTotStack()){          
           h->Add(p->getTotStack());
         }
        
        p->addHistLine(h,sampleNs[iS]);
      } 
    }
    p->rebin(rebin);
    p->draw(false,TString::Format("%s_%s",name.Data(), vars[iV].Data()));
        // p->drawRatio(-1,"stack",false,false,TString::Format("%s_%s",name.Data(), vars[iV].Data()));
  }};
  
  auto makeRocs  = [&](TString * vars, TString name, TString prefix, bool doForward = true ){
    Plotter * p = new Plotter();
    Plotter * pSoB = new Plotter();
  for(unsigned int iV = 0; vars[iV][0]; ++iV){
      TH1 *hb = 0;
      Plotter * pI = new Plotter();      
    for(unsigned int iS = 0; samples[iS][0]; ++iS){
      TH1 *h = 0;
      f->GetObject(TString::Format("%s_%s%s",samples[iS].Data(),prefix.Data(),vars[iV].Data()),h);
      if(h == 0) continue;
      if(iS < nBkg ){
        if(hb == 0) hb = (TH1*)h->Clone();
        else hb->Add(h);
        TH1 *hI = PlotTools::getIntegral(h,doForward,false);                              
        pI->addStackHist(hI,sampleNs[iS].Data());
      }
      else{
        TGraph * roc = PlotTools::getRocCurve(h,hb,doForward,"signal","bkg");
        p->addGraphLine(roc,TString::Format("%s, %s",sampleNs[iS].Data(),vars[iV].Data()));        
        TH1 *hI = PlotTools::getIntegral(h,doForward,true);                              
        float stackNorm = hb->Integral(0,-1);
        hI->Scale(stackNorm);
        pI->addHistLine(hI,sampleNs[iS].Data());
      } 
    }
    pI->draw(false,TString::Format("%s_int_%s",name.Data(), vars[iV].Data()));    
  }
    p->draw(false,name.Data());
  };

  // TString genHbbVars[]= {"gen_Hbb_gen","gen_Hbb_reco","genHbbfj_inclrank","gen_dR_genHbb_hbb","gen_wjj_gen","gen_wjj_reco","genWjjfj_rank","gen_Hbb_recoComp",""};
  // makePlots(genHbbVars,"genHbbVars");
  
  // TString hbbVars[]= {"hbb_fj_mass","hbb_fj_sd_mass","hbb_fj_sdfj_mass","hbb_fj_tau2otau1","hbb_fj_bbcsv","hbb_fj_minsdcsv","hbb_fj_maxsdcsv","hbb_fj_maxMed_minsdcsv","hbb_fj_maxTight_minsdcsv",""};  
  // TString hbbVars[]= {"hbb_fj_oM_mass","hbb_fj_oM_sd_mass","hbb_fj_oM_sdfj_mass","hbb_fj_oM_rawmass","hbb_fj_oM_tau3otau2","hbb_fj_oM_tau3otau1","hbb_fj_oM_tau2otau1","hbb_fj_oM_maxTight_minsdcsv", "hbb_fj_oM_bbcsv","hbb_fj_oM_minsdcsv","hbb_fj_oM_maxsdcsv","hbb_fj_oM_maxMed_minsdcsv",""};
  // makePlots(hbbVars,"hbbVars","mass1400to1800_",5);
    // makePlots(hbbVars,"hbbVars","mass900to1100_",5);
      // makePlots(hbbVars,"hbbVars","",5);

  // TString hbbMass[] = {"hbb_fj_oM_tau2otau1","","",""};
    // makeRocs(hbbMass,"hbbMass","",false);
  // TString hbbMass[] = {"hbb_fj_oM_mass","hbb_fj_oM_sd_mass","hbb_fj_oM_sdfj_mass",""};
  // makeRocs(hbbMass,"hbbMass");
  
  // TString hbbCSV[] = {"hbb_fj_oM_csv","hbb_fj_oM_bbcsv","hbb_fj_oM_maxMed_minsdcsv","",""};
  // TString hbbCSV[] = {"hbb_fj_csv","hbb_fj_bbcsv","hbb_fj_maxMed_minsdcsv","hbb_fj_maxTight_minsdcsv",""};
  // makeRocs(hbbCSV,"hbbCSV");
  
  // TString wjjVars[]= {"wjj_fj_mass","wjj_fj_sd_mass","wjj_fj_sdfj_mass","wjj_fj_tau2otau1","wjj_fj_minsdcsv","wjj_fj_maxsdcsv",""};
  // TString wjjVars[]= {"wjj_fj_oM_mass","wjj_fj_oM_sd_mass","wjj_fj_oM_sdfj_mass","wjj_fj_oM_tau2otau1","wjj_fj_oM_minsdcsv","wjj_fj_oM_maxsdcsv","wjj_fj_oM_csv",""};
  // makePlots(wjjVars,"wjjVars","",5);
  // TString wjjROC[]= {"wjj_fj_oM_maxsdcsv","wjj_fj_oM_csv","wjj_fj_oM_tau2otau1",""};
  // makeRocs(wjjROC,"wjjROC","mass1400to1800_",false);
  
  // TString hhVars[]= {"hWW_mass","hh_mass","tightHbb_hWW_mass","tightHbb_hh_mass","looseHbb_hWW_mass","looseHbb_hh_mass","onShellWjj_hWW_mass","onShellWjj_hh_mass","offShellWjj_hWW_mass","offShellWjj_hh_mass",""};
  // TString hhVars[]= {"tightH_onShell_hWW_mass","tightH_onShell_hh_mass","looseH_onShell_hWW_mass","looseH_onShell_hh_mass","tightH_offShell_hWW_mass","tightH_offShell_hh_mass","looseH_offShell_hWW_mass","looseH_offShell_hh_mass","","",""};
    TString hhVars[]= {"hh_mass","tightHbb_hh_mass","looseHbb_hh_mass","DR_tightHbb_hh_mass","DR_looseHbb_hh_mass",""};
    // TString hhVars[]= {"hh_mass","tightHbb_hh_mass","","looseHbb_hh_mass","onShellWjj_hh_mass","offShellWjj_hh_mass",""};
      // TString hhVars[]= {"dR_wJJ_lepton","dR_hbb_lepton","wJJ_dRpto2","hbb_pt","wJJ_pt","hm_dR_wJJ_lepton","hm_wJJ_dPhipto2","hm_dR_hbb_lepton","hm_wJJ_pt","hm_hbb_pt",""};
  
      // TString hhVars[]= {"neutrino_eta","neutrino_lepton_dEta","neutrino_lepton_dR","wlnu_eta","wlnu_mass","wlnu_wjj_dEta","wlnu_wjj_dR","wlnu_wjj_rPTo2","hWW_eta","hWW_hbb_dEta","hWW_hbb_dR",""};
      // TString hhVars[]= {"neutrino_eta","neutrino_lepton_dEta","neutrino_lepton_dR","wlnu_eta","wlnu_mass","wlnu_wjj_dEta","wlnu_wjj_dR","wlnu_wjj_rPTo2","hWW_eta","hWW_hbb_dEta","hWW_hbb_dR",""};
  
  makePlots(hhVars,"hhVars","",10,false);
        // TString hhVars[]= {"dR_wJJ_lepton","wJJ_dPhipto2","hm_dR_wJJ_lepton","hm_wJJ_dPhipto2",""};
  // TString hhVars[]= {"wlnu_wjj_dR","wlnu_wjj_rPTo2","hm_wlnu_wjj_dR","hm_wlnu_wjj_rPTo2",""};
  // TString hhRocVars[]= {"tightH_onShell_hWW_mass","looseH_onShell_hWW_mass","tightH_offShell_hWW_mass","looseH_offShell_hWW_mass","",""};

  // makeRocs(hhRocVars,"hhRocVars",false);
        
  
  
  // TString lepVars[]= {"lepton_pt","nLeptons","goodLep_nLeptons","ht",""};
  // makePlots(lepVars,"lepVars");
  
  // TString jetVars[]= {"nFJ","maxFHDPhi","nhighDPHI","highDPHI_fj_pt","highDPHI_fj_mass","highDPHI_fj_sd_pt","highDPHI_fj_sd_mass","highDPHI_fj_tau2otau1","highDPHI_fj_csv","highDPHI_fj_minSJPT","highDPHI_fj_minSJCSV","highDPHI_fj_maxSJCSV","highDPHI_fj_oneM_sd_mass","highDPHI_fj_oneM_mass","highDPHI_fj_oneM_tau2otau1","highDPHI_fj_oneM_minSJCSV","highDPHI_fj_oneM_maxSJCSV","nHighDPhi_pass",""};
  // makePlots(jetVars,"jetVars");
  
  // TString genWjjVars[]= {"gen_WWjj_pt","gen_dPHI_lepWWjj","gen_dR_lepWWjj","gen_dPHI_lepWjj","gen_dR_lepWjj",
  // "gen_ptgeq200_dPHI_lepWWjj","gen_ptgeq200_dR_lepWWjj","gen_ptgeq200_dPHI_lepWjj","gen_ptgeq200_dR_lepWjj","",""};
  // makePlots(genWjjVars,"genWjjVars");
  
  // TString wjjVars[]= {"wJJ_cand","wjj_fj_pt"       ,"wjj_fj_mass"     ,"wjj_fj_sd_pt"    ,"wjj_fj_sd_mass"  ,"wjj_fj_tau2otau1","wjj_fj_csv"      ,"wjj_fj_maxSJCSV","wjj_fj_oneM_tau2otau1","wjj_fj_oneM_csv"      ,"wjj_fj_oneM_maxSJCSV"
      // ,"wjj_lepSub_fj_pt"       ,"wjj_lepSub_fj_mass"     ,"wjj_lepSub_fj_sd_pt"    ,"wjj_lepSub_fj_sd_mass"  ,"wjj_lepSub_fj_tau2otau1","wjj_lepSub_fj_csv"      ,"wjj_lepSub_fj_maxSJCSV" ,"wjj_lepSub_fj_oneM_tau2otau1","wjj_lepSub_fj_oneM_csv"      ,"wjj_lepSub_fj_oneM_maxSJCSV" ,""};
  // makePlots(wjjVars,"wjjVars");
  
    // TString metVars[]= {"solved_hWW_mass","solved_hh_mass","loose_solved_hWW_mass","tight_solved_hWW_mass","corr_hWW_mass","corr_hh_mass","tight_corr_hWW_mass","minLepMT","minHadMass","minLepMass","dPhi_met_lep","hh_mass_lt1200_hbb_fj_mass","hh_mass_geq1400_hbb_fj_mass","hh_mass_lt1200_hbb_fj_sd_mass","hh_mass_geq1400_hbb_fj_sd_mass","hh_mass_lt1200_hbb_fj_tau2otau1","hh_mass_geq1400_hbb_fj_tau2otau1",""};
  // TString metVars[]= {"met","dPhi_met_lep"       ,"mt_met_lep"     ,"mt_met_lepWjj" ,"hWW_mass","hWW_pt","hh_mass","hh_pt","hh_dPhi","hh_dEta","solved_hWW_mass","solved_hh_mass","minW_mass","maxW_mass","hWW_chi2",""};
  // TString metVars[] = {"hh_mass_lt1200_wjj_fj_sd_mass","hh_mass_lt1200_wjj_fj_tau2otau1","hh_mass_lt1200_wjj_fj_maxSJCSV","hh_mass_lt1200_wjj_fj_minSJCSV","hh_mass_lt1200_hbb_fj_sd_mass","hh_mass_lt1200_hbb_fj_tau2otau1","hh_mass_lt1200_hbb_fj_maxSJCSV","hh_mass_lt1200_hbb_fj_minSJCSV","hh_mass_lt1200_nFJ","hh_mass_lt1200_nAK4","hh_mass_lt1200_nExtraAK4","hh_mass_lt1200_max_extra_ak4", "hh_mass_geq1400_wjj_fj_sd_mass","hh_mass_geq1400_wjj_fj_tau2otau1","hh_mass_geq1400_wjj_fj_maxSJCSV","hh_mass_geq1400_wjj_fj_minSJCSV","hh_mass_geq1400_hbb_fj_sd_mass","hh_mass_geq1400_hbb_fj_tau2otau1","hh_mass_geq1400_hbb_fj_maxSJCSV","hh_mass_geq1400_hbb_fj_minSJCSV","hh_mass_geq1400_nFJ","hh_mass_geq1400_nAK4","hh_mass_geq1400_nExtraAK4","hh_mass_geq1400_max_extra_ak4",""};
  
  
    // TString metVars[]= {"solved_hh_mass","tight_solved_hh_mass","corr_hh_mass","tight_corr_hh_mass","corrPart_hh_mass","solved_hWW_mass","tight_corrPart_hh_mass",""};
    // makePlots(metVars,"metVars",false);
  //
    // TString ak4wjjVars[]= {"gen_min_Wj","gen_max_Wj","gen_dR_Wj1j2"       ,"gen_min_dR_lepWWj"  ,"gen_max_dR_lepWWj"  ,"gen_min_dR_lepWj"   ,"gen_max_dR_lepWj"   ,"gen_min_dPhi_lepWWj","gen_max_dPhi_lepWWj",""};
      // TString ak4wjjVars[]= {"nak4_nearW","nak4_pairs_nearW","ak4_1j_mass","ak4_1j_csv","ak4_2j_mass","ak4_2j_mincsv","ak4_2j_maxcsv","ak4_categories","wJJ_cand",""};
            // TString ak4wjjVars[]= {"ak42jpass_corr_hWW_mass","ak42jpass_corr_hh_mass","ak41jpass_corr_hWW_mass","ak41jpass_corr_hh_mass",""};
  // makePlots(ak4wjjVars,"ak4wjjVars",false);

}

{
  TH1 * hb = ttbar_met;
  TH1 * hs = signal_1600_met;
  Plotter * p = new Plotter();
  p->addHistLine(hs,"sig");
  p->addHistLine(hb,"bkg");
  p->normalize();
  p->draw(false,"dist");
  
  TGraph * roc = PlotTools::getRocCurve(hs,hb,true,"signal","t#bar{t}");
  TGraph * roc2 = PlotTools::getRocCurve(signal_1600_mt_met_lepWjj,ttbar_mt_met_lepWjj,true,"signal","t#bar{t}");
  Plotter *pr = new Plotter();
  pr->addGraphLine(roc,"met");
  pr->addGraphLine(roc2,"mt");
  pr->draw(false,"roc");
  
}


{
    TFile * f = new TFile("compiled_postHbb.root","read");
    TString samples[] = {"signal_1000","signal_1600",""};
    TString sampleNs[] = {"M(1000)","M(1600)"};

    // TString vars[] = {"solved_hh_mass","corrPart_hh_mass","corr_hh_mass",""};    
    // TString vars[] = {"tight_solved_hh_mass","tight_corrPart_hh_mass","tight_corr_hh_mass",""};
    // TString varNs[] = {"uncorrected","correct hbb","correct hbb and #slash{E}_{T}",""};
    TString vars[] = {"solved_hh_mass","corrPart_hh_mass","",""};    
    TString varNs[] = {"uncorrected","after h#rightarrowbb calib.","",""};

      Plotter * p = new Plotter();    
      for(unsigned int iS = 0; samples[iS][0]; ++iS){
    for(unsigned int iV = 0; vars[iV][0]; ++iV){
      TH1 *h = 0;
      f->GetObject(TString::Format("%s_%s",samples[iS].Data(),vars[iV].Data()),h);
      if(h == 0) continue;
      cout << TString::Format("%s_%s",samples[iS].Data(),vars[iV].Data()) <<" -> "<< h->GetMean() <<" "<< h->GetStdDev() <<endl;
      int lineStyle = 1; if(iS) lineStyle = 9;
      p->addHistLine(h,sampleNs[iS] +", "+ varNs[iV],-1,lineStyle );      
      
    }
  }
  p->normalize();
  p->rebin(10); 
  p->draw();

  
}

sd 100-175
  tau 2/1 < 0.7
    max subjet csv > 0.6
      