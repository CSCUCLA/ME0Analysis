{
        TFile * f = new TFile("test_plots_PU_p8s768Merge.root");
	TString samps[] = {"p8s768Merge","","p8s640Merge","p8s512Merge","p8s384Merge","p8s256Merge","p8s128Merge",""};
    TString sampNs[] =     {"768 strips, 8 part.","640 strips, 8 part.","512 strips, 8 part.","384 strips, 8 part.","256 strips, 8 part.","128 strips, 8 part."};
    // TString var = "ptgeq20_segmentcheck_simHit_newSeg_dphidiff";
    TString var = "pteq5to20_segmentcheck_simHit_newSeg_dphidiff";

        // TString var = "pteq3to5_segmentcheck_simHit_newSeg_dphidiff";

    Plotter * p = new Plotter;
    
    for(unsigned int iS = 0; samps[iS][0]; ++iS){
        TH1 * h;
        f->GetObject(TString::Format("%s_%s",samps[iS].Data(),var.Data()),h);
        h->Rebin(10);
        // h->SetXTitle("#Delta#phi[lay_{6}-lay_{1}](SimHit) - #Delta#phi[lay_{6}-lay_{1}](New seg.)");
        h->SetXTitle("#phi (SimHit) - #phi (New seg.)");
        p->addHist(h,sampNs[iS]);

    }
    p->draw();

}


{
        TString prefix = "p64s4Merge";
    TFile * f = new TFile("digiTester_plots_" + prefix+".root");
    // TString pts[] = {"ptleq3","pteq3to5","pteq5to20","ptgeq20",""};
    // TString ptNs[] =     {"p_{T}<3","3-p_{T}-5","5-p_{T}-20","p_{T}>20"};

	TString pts[] = {"pteq5to20","ptgeq20",""};
    TString ptNs[] =     {"5-p_{T}-20","p_{T}>20"};

    // TString hts[] = {"nHits3","nHits4to5","nHits6",""};
    // TString htNs[] =     {"3 hits","4-5 hits","6 hits",""};

	TString hts[] = {"nHits4to5","nHits6",""};
    TString htNs[] =     {"4-5 hits","6 hits",""};
    
    TString vars[] ={
        // "trueProperties_avgTOF","trueProperties_seedTOFDiff","trueProperties_TOFDiff",
                    // "trueProperties_seedETADiff","trueProperties_ETADiff",
                        // "trueProperties_seedYDiff","trueProperties_YDiff",
                                            // "trueProperties_seedRDiff","trueProperties_RDiff",
                    // "trueProperties_seedPHIDiff","trueProperties_PHIDiff",
            "trueProperties_seednETADiff","trueProperties_nETADiff","","trueProperties_nETACen",
                    ""};

    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        Plotter * p = new Plotter;
    
        for(unsigned int iP = 0; pts[iP][0]; ++iP){
                    for(unsigned int iH = 0; hts[iH][0]; ++iH){
            TH1 * h;
            f->GetObject(TString::Format("%s_%s_%s__%s",prefix.Data(), pts[iP].Data(),hts[iH].Data(),vars[iV].Data()),h);
            if(h == 0) continue;
            // h->Rebin(10);
            // h->SetXTitle("#Delta#phi[lay_{6}-lay_{1}](SimHit) - #Delta#phi[lay_{6}-lay_{1}](New seg.)");
            // h->SetXTitle("#phi (SimHit) - #phi (New seg.)");
            p->addHist(h,TString::Format("%s, %s",ptNs[iP].Data(),htNs[iH].Data()));

        }
    }
        p->draw();
    }
    
    
    
}



{
    // TFile * f = new TFile("eta_digiall.root");
    TFile * f = new TFile("phi_digiall.root");
    // TString pts[] = {"ptleq3","pteq3to5","pteq5to20","ptgeq20",""};
    // TString ptNs[] =     {"p_{T}<3","3-p_{T}-5","5-p_{T}-20","p_{T}>20"};
    // TString digis[] = {"p4s4Merge","p4s768Merge","p8s768Merge","p16s768Merge","p64s768Merge","p64s4Merge","p128s768Merge",""};
        TString digis[] = {"p8s128Merge","p8s256Merge","p8s512Merge","p4s512Merge","p8s1024Merge","","",""};
	TString pts[] = {"pteq3to5","pteq5to20",""};
    TString ptNs[] =     {"3-p_{T}-5","5-p_{T}-20"};
    
    // TString hts[] = {"nHits3","nHits4to5","nHits6",""};
    // TString htNs[] =     {"3 hits","4-5 hits","6 hits",""};

	TString hts[] = {"nHits4to5","nHits6",""};
    TString htNs[] =     {"4-5 hits","6 hits",""};
    
    TString vars[] ={ "trueProperties_seedPHINDiff"
            // "trueProperties_PHICenNDiff",
            // "trueProperties_ETASeedFitDiff",
                // "trueProperties_TOFDiff",
                    ""};

for(unsigned int iD = 0; digis[iD][0]; ++iD){
    TString prefix = digis[iD];
    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        Plotter * p = new Plotter;
    
        for(unsigned int iP = 0; pts[iP][0]; ++iP){
                    for(unsigned int iH = 0; hts[iH][0]; ++iH){
            TH1 * h;
            f->GetObject(TString::Format("%s_%s_%s__%s",prefix.Data(), pts[iP].Data(),hts[iH].Data(),vars[iV].Data()),h);
            if(h == 0) continue;
            // h->Rebin(10);
            // h->SetXTitle("#Delta#phi[lay_{6}-lay_{1}](SimHit) - #Delta#phi[lay_{6}-lay_{1}](New seg.)");
            // h->SetXTitle("#phi (SimHit) - #phi (New seg.)");
            h->SetTitle(digis[iD]);
            p->addHist(h,TString::Format("%s, %s",ptNs[iP].Data(),htNs[iH].Data()));
        
        }
    }
        p->draw(false,digis[iD]);
    }
}
    
    
    
}



{
    TFile * f = new TFile("testSegments.root");

        TString digis[] = {"p8s768","p12s256","p8s384","p6s512","p4s768","","",""};
        TString digiNs[] = {"8 part. 768 strips","12 part. 256 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};
        
        // TString digis[] = {"p8s768","p8s768M6","p8s768L3","p8s768M6L3","","","",""};
        // TString digiNs[] = {"6 layers, req. 4 layers per seg.","5 layers, req. 4 layers per seg.","6 layers, req. 3 layers per seg.","5 layers, req. 3 layers per seg.","",""};
    
    // TString digis[] = {"p8s768","p8s768L3","","","","",""};
    // TString digiNs[] = {"Req. 4 layers per seg.","Req. 3 layers per seg.","",""};
    
        // TString digis[] = {"p8s256","p8s384","p8s512","p8s640","p8s768","","",""};
        // TString digiNs[] = {"8 part. 256 strips","8 part. 384 strips","8 part. 512 strips","8 part. 640 strips","8 part. 768 strips","","",""};
    
    // TString digis[] = {"p2s768","p4s768","p6s768","p8s768","p10s768","","",""};
    // TString digiNs[] = {"2 part. 768 strips","4 part. 768 strips","6 part. 768 strips","8 part. 768 strips","10 part. 768 strips","","",""};

        // TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        // TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};
        


    TString vars[] ={ "fake_passDPhiTime_nSegs","","all_fake_seg_nNon0","all_fake_seg_tof" ,"all_fake_seg_nHits","all_fake_seg_dPhi","all_fake_seg_chi2ondof",
            // "trueProperties_PHICenNDiff",
            // "trueProperties_ETASeedFitDiff",
                // "trueProperties_TOFDiff",
                    ""};

    for(unsigned int iV = 0; vars[iV][0]; ++iV){
                Plotter * p = new Plotter;
                for(unsigned int iD = 0; digis[iD][0]; ++iD){
    TString prefix = digis[iD];
    TH1 * h;
    f->GetObject(TString::Format("%s__%s",prefix.Data(),vars[iV].Data()),h);
    if(h == 0) continue;
    if(iV == 0) h->SetXTitle("# of bkg. segments per event");
    if(iV == 0) h->SetYTitle("a.u.");
    cout << prefix.Data() << " "<< vars[iV].Data()<<" "<<h->GetMean() <<endl;
    // if(iV == 0) h->Rebin(2);
    h->SetTitle(digis[iD]);
    p->addHistLine(h,TString::Format("%s",digiNs[iD].Data()));
    }
        p->draw(false,vars[iV]);
    }
}


// Fake TDR plot by eta

{
    TFile * f = new TFile("segmentAnalyzerForTDR_p8s384.root");
    
    TString digi = "p8s384";
    
        // TString dPhis[] = {"_","_dPhi_lt0p013_","_dPhi_lt0p004_","_dPhi_lt0p002_","","","",""};
        // TString dPhiNs[] = {"|#Delta#phi| < 0.02 (p_{T} > 2 GeV)","|#Delta#phi| < 0.013 (p_{T} > 3 GeV)","|#Delta#phi| < 0.004 (p_{T} > 10 GeV)","|#Delta#phi| < 0.002 (p_{T} > 20 GeV)","","",""};
        
        TString dPhis[] = {"_","_dPhi_lt0p002scaled_","_dPhi_lt0p013_","_dPhi_lt0p004_","_dPhi_lt0p002_","","",""};
        TString dPhiNs[] = {"|#Delta#phi| < 0.02 (p_{T} > 2 GeV)","|#Delta#phi| < C(#eta) (p_{T} > 2 GeV)","|#Delta#phi| < 0.013 (p_{T} > 3 GeV)","|#Delta#phi| < 0.004 (p_{T} > 10 GeV)","|#Delta#phi| < 0.002 (p_{T} > 20 GeV)","","",""};
        


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
    h->Rebin(10);
    h->Scale(1./nEvents);
    // h->SetXTitle("# of bkg. segments per event");
    h->SetYTitle("<N. of bkg. segments> / event");
    cout << dPhis[iD].Data()<< " "<< vars[iV].Data()<<" "<<h->Integral(0,-1) <<endl;
    // if(iV == 0) h->Rebin(2);
    // h->SetTitle(digis[iD]);
    p->addHistLine(h,dPhiNs[iD]);
    }
        p->draw(false,vars[iV]);
    }
}


{
    TFile * f = new TFile("testTestSegments.root");

        // TString digis[] = {"p12s256","p8s768","p8s384","p6s512","p4s768","","",""};
        // TString digiNs[] = {"12 part. 256 strips","8 part. 768 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};

        TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};

        TString cats[] = {"all","perfect","oneFakeHit","twoFakeHits","lost","",""};
        // TString cats[] = {"all","MUON_COMP_PURE","MUON_COMP_DIRTY_TRACK","MUON_COMP_DIRTY_NEUT","MUON_MISS_PURE","MUON_MISS_DIRTY_TRACK","MUON_MISS_DIRTY_NEUT","MISSING_SEGMENT",""};


    
    TString vars[] ={ "pt",""};

    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
                Plotter * p = new Plotter();     int iH = 0;
                for(unsigned int iC = 0; cats[iC][0]; ++iC){
    TH1 * h = 0;
    f->GetObject(TString::Format("%s__%s_muon_%s",digis[iD].Data(),cats[iC].Data(),vars[iV].Data()),h);
    if(h == 0) continue;
    ++iH;
    h->Rebin(4);
    h->SetTitle(digis[iD]);
    p->addHist(h,TString::Format("%s",cats[iC].Data()));
    }
    // false,digis[iD])
        if(iH)p->drawRatio(0,"s",true,false,digis[iD]);
    // p->draw(false,digis[iD]);
    }
}

}


{
    TFile * f = new TFile("testSegments.root");

    TString digis[] = {"p8s768","p12s256","p8s384","p6s512","p4s768","","",""};
    TString digiNs[] = {"8 part. 768 strips","12 part. 256 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};

        // TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        // TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};
        TString den = "all";

        TString numcats[] = {"perfect","oneFakeHit","twoFakeHits","","",""};


    
    TString vars[] ={ "pt",""};

    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        Plotter * p = new Plotter();    
    TH1 * hd = 0;        
    f->GetObject(TString::Format("%s__%s_mCh_muon_%s",digis[0].Data(),den.Data(),vars[iV].Data()),hd);        
    hd->SetXTitle("muon p_{T} [GeV]");
    hd->SetYTitle("signal ME0Segment reconstruction efficiency");
    p->addHist(hd,"all");        
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            TH1 * hn = 0;
                int iH = 0;
                for(unsigned int iC = 0; numcats[iC][0]; ++iC){
    TH1 * h = 0;
    f->GetObject(TString::Format("%s__%s_mCh_passDPhiTime_muon_%s",digis[iD].Data(),numcats[iC].Data(),vars[iV].Data()),h);
    if(h == 0) continue;
    ++iH;
    if(hn) hn->Add(h); else hn = (TH1*)h->Clone();
    // h->SetTitle(digis[iD]);
    }
    if(hn) p->addHist(hn,digiNs[iD]);
    // false,digis[iD])
     
    // p->draw(false,digis[iD]);
    }
    p->rebin(2);
    p->drawRatio(0,"s",true,false);
}

}



{
    TFile * f = new TFile("testSegments.root");

        TString digis[] = {"p12s256","p8s768","p8s384","p6s512","p4s768","","",""};
        TString digiNs[] = {"12 part. 256 strips","8 part. 768 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};
        
        // TString cats[] = {"pteq1to3_perfect_muon","pteq3to5_perfect_muon","pteq5to20_perfect_muon","ptgeq20_perfect_muon","all_fake",""};

    
    TString vars[] ={ "dPhi","","tof","simMreco_phi","simMreco_eta","simMreco_dphi","simMreco_deta",""};
    TString fakeStr = "all_fake";
    TString ptstr[] = {"pteq1to3","pteq3to5","pteq5to20","ptgeq20","",""};
    TString ptstrNs[] =     {"#mu signal, p_{T} 1-3 GeV","#mu signal, p_{T} 3-5 GeV","#mu signal, p_{T} 5-20 GeV","#mu signal, p_{T} 20-30 GeV"};
    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
                Plotter * p = new Plotter();     
                
                TH1 * hf = 0;
                f->GetObject(TString::Format("%s__%s_seg_%s",digis[iD].Data(),fakeStr.Data(),vars[iV].Data()),hf);
                if(hf){
                    p->addStackHist(hf,"background ME0Segments");                    
                }
                
                
                
                for(unsigned int iP = 0; ptstr[iP][0]; ++iP){
    TH1 * hp = 0;
    f->GetObject(TString::Format("%s__%s_perfect_muon_seg_%s",digis[iD].Data(),ptstr[iP].Data(),vars[iV].Data()),hp);
    TH1 * hd = 0;
    f->GetObject(TString::Format("%s__%s_dirty_muon_seg_%s",digis[iD].Data(),ptstr[iP].Data(),vars[iV].Data()),hd);
    
    if(hp == 0) continue;
    if(hd) hp->Add(hd);

    hp->SetTitle(digis[iD]);
    p->addHistLine(hp,ptstrNs[iP].Data());
    }
    // false,digis[iD])
        // if(iH)p->drawRatio(0,"s",true,false,digis[iD]);
    p->rebin(2);
    p->normalize();
    p->setMinMax(0.001,1);
    p->draw(false,digis[iD] +"_"+ vars[iV]);
    }
}

}

//Table

{
    // TFile * f = new TFile("testSegments_justNeut_112p5_fineP.root");
        TFile * f = new TFile("testSegments_TDR_justNeut_112p5_fineP.root");

        // TString digis[] = {"p8s768","p8s768M6","p8s768L3","p8s768M6L3","","","",""};
        // TString digiNs[] = {"6 layers, req. 4 layers per seg.","5 layers, req. 4 layers per seg.","6 layers, req. 3 layers per seg.","5 layers, req. 3 layers per seg.","",""};
    
        // TString digis[] = {"p8s256","p8s384","p8s512","p8s640","p8s768","","",""};
        // TString digiNs[] = {"8 part. 256 strips","8 part. 384 strips","8 part. 512 strips","8 part. 640 strips","8 part. 768 strips","","",""};
    
    // TString digis[] = {"p2s768","p4s768","p6s768","p8s768","p10s768","","",""};
    // TString digiNs[] = {"2 part. 768 strips","4 part. 768 strips","6 part. 768 strips","8 part. 768 strips","10 part. 768 strips","","",""};
    

        // TString digis[] = {"p8s768","p12s256","p8s384","p6s512","p4s768","","",""};
        // TString digiNs[] = {"8 part. 768 strips","12 part. 256 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};
        
       //   TString digis[] = {"p6s51210NL3",
       //  "p6s5127p5NL3"    ,
       //  "p6s51210NL4"     ,
       //  "p6s5127p5NL4"    ,
       //  "p6s51210NL5"     ,
       //  "p6s5127p5NL5"    ,
       //  "p6s51210NL6"     ,
       //  "p6s5127p5NL6"    ,
       //  "p6s51210NM6L3"   ,
       //  "p6s5127p5NM6L3"  ,
       //  "p6s51210NM6L4"   ,
       //  "p6s5127p5NM6L4"  ,
       //  "p6s51210NM6L5"   ,
       //  "p6s5127p5NM6L5"  ,
       //  "p6s51210NM56L3"  ,
       //  "p6s5127p5NM56L3" ,
       //  "p6s51210NM56L4"  ,
       //  "p6s5127p5NM56L4" ,
       //  "p6s51210NM456L3" ,
       //  "p6s5127p5NM456L3",""};
       //  TString digiNs[] = {"10N L3",
       // "7p5N L3"    ,
       // "10N L4"     ,
       // "7p5N L4"    ,
       // "10N L5"     ,
       // "7p5N L5"    ,
       // "10N L6"     ,
       // "7p5N L6"    ,
       // "10N M6L3"   ,
       // "7p5N M6L3"  ,
       // "10N M6L4"   ,
       // "7p5N M6L4"  ,
       // "10N M6L5"   ,
       // "7p5N M6L5"  ,
       // "10N M56L3"  ,
       // "7p5N M56L3" ,
       // "10N M56L4"  ,
       // "7p5N M56L4" ,
       // "10N M456L3" ,
       // "7p5N M456L3",""};
    
  //   TString digis[] = {
  //  "p6s5127p5NL4"    ,
  //  "p6s51210NL4"     ,
  //  "p6s51237p5NL4"     ,
  //  "p6s51275NL4"     ,
  //  ""};
  //  TString digiNs[] = {
  //    "p6s5127p5NL4"    ,
  //    "p6s51210NL4"     ,
  //    "p6s51237p5NL4"     ,
  //    "p6s51275NL4"     ,
  // ""};
        
        TString digis[] = {
          "p8s3847p5NL4"  ,
          "p8s38410NL4"   ,
          "p8s38422p5NL4" ,
          "p8s38437p5NL4" ,     
          "p8s38475NL4"   ,
          "p8s3847p5NL5"  ,
          "p8s38410NL5"   ,
          "p8s38422p5NL5" ,
          "p8s38437p5NL5" ,
          "p8s38475NL5"   ,
       ""};
       TString digiNs[] = {     
         "p8s3847p5NL4"  ,
         "p8s38410NL4"   ,
         "p8s38422p5NL4" ,
         "p8s38437p5NL4" ,     
         "p8s38475NL4"   ,
         "p8s3847p5NL5"  ,
         "p8s38410NL5"   ,
         "p8s38422p5NL5" ,
         "p8s38437p5NL5" ,
         "p8s38475NL5"   ,
      ""};






        

        // TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        // TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};


        cout <<" \t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            cout << digiNs[iD] <<"\t";
        }
        cout <<endl;
        
        //Do fake rate
        cout << "Fake Rate (mean,sigma)\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            TH1 * h = 0;
            f->GetObject(TString::Format("%s__fake_passDPhiTime_ext_nSegs",digis[iD].Data()),h);
            float mean = 0;
            float std = 0;
            if(h){
                mean = h->GetMean();
                std  = h->GetStdDev();
            }
                 cout << TString::Format("%.3f,%.3f\t",mean,std);
        }
        cout <<endl;
        //fake proportions
        TString faketypes[] = { "FAKE_TRACK_PURE", "FAKE_NEUT_PURE", "FAKE_MUON_MIX", "FAKE_OTHER",""};
        TString fakeden = "all";
            for(unsigned int iN = 0; faketypes[iN][0]; ++iN){
                cout << faketypes[iN] <<" (%)\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            
            TH1 * h = 0;
            TH1 * hd = 0;
            f->GetObject(TString::Format("%s__%s_fake_seg_nHits",digis[iD].Data(),fakeden.Data()),hd);
            f->GetObject(TString::Format("%s__%s_fake_seg_nHits",digis[iD].Data(),faketypes[iN].Data()),h);
            float num = 0;
            float tot = 0;
            float avg = 0;
            if(h && hd){
                tot = hd->Integral();
                num  = h->Integral();
                avg = num/tot * 100.0;
            }
                 cout << TString::Format("%.0f%%\t",avg);
        }
        cout << endl;
    }
    cout <<endl;
        

        //efficiency
        TString pts[] = {"pteq1to3","pteq3to5","pteq5to20","ptgeq20",""};
        TString ptNs[] = {"1-pT-3","3-pT-5","5-pT-20","20-pT-30",""};
        TString effnums[] = {"perfect","lost",""};
        TString effquals[] = {"seg_passDPhiTime_nHits","seg_passDPhiTime_nHits",""};
        TString effnumNs[] = {"perfect (%)","lost (%)",""};
        for(unsigned int iP = 0; pts[iP][0]; ++iP){
            for(unsigned int iN = 0; effnums[iN][0]; ++iN){
            cout << ptNs[iP] <<" "<<effnumNs[iN]<< "\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){


            TH1 * hd = 0;
            f->GetObject(TString::Format("%s__%s_all_muon_seg_nHits",digis[iD].Data(),pts[iP].Data()),hd);
            if(!hd) {
              cout << TString::Format("%.1f%%\t",0.0); continue;
            }
            float tot = hd->Integral();
            
            float num = 0;
            if(iN){
                TH1 * hp = 0;
                f->GetObject(TString::Format("%s__%s_perfect_muon_%s",digis[iD].Data(),pts[iP].Data(),effquals[iN].Data()),hp);
            
                TH1 * hdirty = 0;
                f->GetObject(TString::Format("%s__%s_dirty_muon_%s",digis[iD].Data(),pts[iP].Data(),effquals[iN].Data()),hdirty);  
                num = hd->Integral();
                if(hp) num -= hp->Integral();
                if(hdirty) num -= hdirty->Integral();
                
            } else {

                TH1 * h = 0;
                f->GetObject(TString::Format("%s__%s_%s_muon_%s",digis[iD].Data(),pts[iP].Data(),effnums[iN].Data(),effquals[iN].Data()),h);
                if(h)num = h->Integral();
            }
            float eff = num/tot * 100.0;
            
                 cout << TString::Format("%.1f%%\t",eff);
                 // cout << TString::Format("%.1f%%\t",TMath::Sqrt(num)/tot * 100);
                 
                 // cout << TString::Format("(%f,%f)\t",num,tot);
        }
        cout << endl;
    }

}
cout << endl;


        //lost muon causes
        TString effLostCausesNums[] = {"MUON_COMP_DIRTY_TRACK","MUON_COMP_DIRTY_NEUT",
        "MUON_MISS_PURE","MUON_MISS_DIRTY_TRACK","MUON_MISS_DIRTY_NEUT","MISSING_SEGMENT",""};

        for(unsigned int iP = 0; pts[iP][0]; ++iP){
            for(unsigned int iN = 0; effLostCausesNums[iN][0]; ++iN){
            cout << ptNs[iP] <<" "<<effLostCausesNums[iN]<< "\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){

            TH1 * hd = 0;
            for(unsigned int iN2 = 0; effLostCausesNums[iN2][0]; ++iN2){
                TH1 * h = 0;
                f->GetObject(TString::Format("%s__%s_%s_muon_seg_nHits",digis[iD].Data(),pts[iP].Data(),effLostCausesNums[iN2].Data()),h);
                if(!h) continue;
                if(hd) hd->Add(h);
                else hd = (TH1*)h->Clone();
            }
            TH1 * h = 0;
            f->GetObject(TString::Format("%s__%s_%s_muon_seg_nHits",digis[iD].Data(),pts[iP].Data(),effLostCausesNums[iN].Data()),h);
            float eff = 0;
            float num = 0;
            float tot = 0;
            if(h && hd){
                tot = hd->Integral();
                num = h->Integral();
                eff = num/tot * 100.0;
            }
                 cout << TString::Format("%.1f%%\t",eff);
                 // cout << TString::Format("(%f,%f)\t",num,tot);
        }
        cout << endl;
    }

}
cout << endl;


        //distributions
        TString distTypes[] = {"perfect","dirty",""};
        TString distNs[] = {"perfect","dirty",""};
        TString vars[] = {"phi","eta","dphi","deta",""};
        TString varNs[] = {"SmS(phi)","SmS(eta)","SmS(dphi)","SmS(deta)"};

        for(unsigned int iV = 0; vars[iV][0]; ++iV){
        for(unsigned int iP = 0; pts[iP][0]; ++iP){
            for(unsigned int iN = 0; distNs[iN][0]; ++iN){
            cout << ptNs[iP] <<" "<<distNs[iN]<< " "<< varNs[iV]<<" (sigma)\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            TH1 * h = 0;
            f->GetObject(TString::Format("%s__%s_%s_muon_seg_simMreco_%s",digis[iD].Data(),pts[iP].Data(),distTypes[iN].Data(),vars[iV].Data()),h);
            float mean = 0;
            float std = 0;
            if(h){
                mean = h->GetMean();
                std  = h->GetStdDev();
            }
                 cout << TString::Format("%.1e\t",std);
        }
        cout << endl;
    }

}
cout << endl;

}
}



////FAKE MULTIPLICITY TABLE
{
    TFile * f = new TFile("testSegments_justNeut_112p5_0p1P.root");
    TString STLAY[] = {"p6s512","p8s768","p12s512"};
    TString GEOM[] = {"L3","L4" ,"L5","L6","M6L3"  ,"M6L4"   ,"M6L5"   ,"M56L3"  ,"M56L4"  ,"M456L3" ,"L4T0" ,"L5T0","M6L3T0","M6L4T0","L4" ,"L5", "L4","L5" ,""};
    TString GEOMN[] = {"6 lay., 3 per seg.","6 lay., 4 per seg." ,"6 lay., 5 per seg.","6 lay., 6 per seg." ,"5 lay., 3 per seg."  ,"5 lay., 4 per seg."   ,"5 lay., 5 per seg."   ,"4 lay., 3 per seg."  ,"4 lay., 4 per seg."  ,"3 lay., 3 per seg.","6 lay., 4 per seg. 1BX" ,"6 lay., 5 per seg. 1BX","5 lay., 3 per seg. 1BX"  ,"5 lay., 4 per seg. 1BX", "48S 6 lay., 4 per seg. ","48S 6 lay., 5 per seg. ","48E 6 lay., 4 per seg. "  ,"48E 6 lay., 5 per seg. "  ,""};
    TString LUM[] = {"7p5N","10N","22p5N","37p5N","75N",""};
    TString LUMN[] = {"7.5x10^34","10x10^34","22.5x10^34","37.5x10^34","75x10^34",""};
    TString VARS[] = {"Ineff","BKG"};
    int nVars = 2;
    
    cout <<"\t";
    int midV = nVars/2;
    for(unsigned int iL = 0; LUM[iL][0]; ++iL){
      for(unsigned int iV = 0; iV < nVars; ++iV){
        if(midV == iV) cout << LUMN[iL];
        cout <<"\t";
      }
    }
    cout <<endl<<"\t";
    for(unsigned int iL = 0; LUM[iL][0]; ++iL){
      for(unsigned int iV = 0; iV < nVars; ++iV){
        cout << VARS[iV] <<"\t";
      }
    }
    cout <<endl;
    
    for(unsigned int iG = 0; GEOM[iG][0]; ++iG){
     cout << GEOMN[iG]<<"\t";
         for(unsigned int iL = 0; LUM[iL][0]; ++iL){
           TString thisStrip = STLAY[0];
           if(iG == 14 || iG == 15) thisStrip =STLAY[1];
           if(iG == 16 || iG == 17) thisStrip =STLAY[2];
                TString prefix = thisStrip + LUM[iL] + GEOM[iG];
                
           //EFFICIENCY
           int effBin = 7; //3 GeV;
           TH1 * hd = 0;
           f->GetObject(TString::Format("%s__all_mCh_muon_pt",prefix.Data()),hd);
           TH1 * h1 = 0;
           f->GetObject(TString::Format("%s__perfect_mCh_passDPhiTime_muon_pt",prefix.Data()),h1);
           TH1 * h2 = 0;
           f->GetObject(TString::Format("%s__oneFakeHit_mCh_passDPhiTime_muon_pt",prefix.Data()),h2);
           TH1 * h3 = 0;
           f->GetObject(TString::Format("%s__twoFakeHits_mCh_passDPhiTime_muon_pt",prefix.Data()),h3);
           float den = 0;
           float num = 0;
           if(hd) den = hd->Integral(effBin,-1);
           if(h1) num += h1->Integral(effBin,-1);
           if(h2) num += h2->Integral(effBin,-1);
           if(h3) num += h3->Integral(effBin,-1);
           num = den - num;
           if(den){
             float eff = num/den*100;
              cout << TString::Format("%.1f%%\t",eff);                                                         
           } else {
             cout << TString::Format("-\t");                                                         
           }
           
           
           //Do fake rate

               TH1 * hf = 0;
               f->GetObject(TString::Format("%s__fake_passDPhiTime_ext_nSegs",prefix.Data()),hf);
               float mean = 0;
               float std = 0;
               
               if(hf){
                   mean = hf->GetMean();
                   std  = hf->GetStdDev();
                   cout << TString::Format("%.1f\t",mean);
                   
               } else {
                 cout << TString::Format("-\t");
               }
           
         }
         
         
         
         
    cout << endl;     
    }

       
       
}



//MAKE NEUTRON PLOT

{
  vector<TString> xTitles{"7.5 (x1)","10 (x1.33)","22.5 (x3)","37.5 (x5)","75 (x10)"};
  vector<TString> names{"4 out of 6 layers per segment","5 out of 6 layers per segment"};
  vector<vector<double>> vals{{0.198,0.282,5.298,29.626,195.550},{0.008,0.023,0.214,1.748,39.374}};
  Plotter * p = new Plotter();
  
  TH1F * axisHist = new TH1F("axisHist",";neutron luminosity scale [10^{34} Hz/cm] (safety margin)",xTitles.size(),-0.5,float(xTitles.size()) -.5);
  for (i=1;i<=xTitles.size();i++) axisHist->GetXaxis()->SetBinLabel(i,xTitles[i-1]);
  
  for(unsigned int iN = 0; iN < names.size(); ++iN){
    TGraph * g = new TGraph();
    for(unsigned int iP = 0; iP < vals[iN].size(); ++iP){
      g->SetPoint(iP,iP,vals[iN][iP]);
    }
    g->SetLineColor  (StyleInfo::getLineColor(iN));
    g->SetLineWidth  (3);
    g->SetLineStyle  (1);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(StyleInfo::getLineColor(iN));
    g->SetMarkerSize (1);
    Drawing::Drawable1D drawableF("P L",names[iN],Drawing::GRAPH,g,false);
    drawableF.graphAxisHist = (TH1*)axisHist->Clone();
    p->addDrawable(drawableF);     
  }
  p->setYTitle("<N. of neutron bkg. segments> per event");
  p->setXTitle("neutron luminosity scale [10^{34}cm^{-2}s^{-1}] (safety margin)");
  p->draw(false,"tot");


}


//Make bkg composition table
  
