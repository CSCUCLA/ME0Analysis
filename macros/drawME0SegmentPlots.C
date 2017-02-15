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

        // TString digis[] = {"p12s256","p8s768","p8s384","p6s512","p4s768","","",""};
        // TString digiNs[] = {"12 part. 256 strips","8 part. 768 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};

        TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};
        


    TString vars[] ={ "fake_nSegs","all_fake_seg_nHits","all_fake_seg_dPhi","all_fake_seg_chi2ondof",
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
    cout << prefix.Data() << " "<< vars[iV].Data()<<" "<<h->GetMean() <<endl;
    if(iV == 0) h->Rebin(2);
    h->SetTitle(digis[iD]);
    p->addHist(h,TString::Format("%s",digiNs[iD].Data()));
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

        TString digis[] = {"p12s256","p8s768","p8s384","p6s512","p4s768","","",""};
        TString digiNs[] = {"12 part. 256 strips","8 part. 768 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};
        
        TString cats[] = {"pteq1to3_perfect_muon","pteq3to5_perfect_muon","pteq5to20_perfect_muon","ptgeq20_perfect_muon","all_fake",""};

    
    TString vars[] ={ "dPhi","simMreco_phi","simMreco_eta","simMreco_dphi","simMreco_deta",""};
    
    

    for(unsigned int iV = 0; vars[iV][0]; ++iV){
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
                Plotter * p = new Plotter();     int iH = 0;
                for(unsigned int iC = 0; cats[iC][0]; ++iC){
    TH1 * h = 0;
    f->GetObject(TString::Format("%s__%s_seg_%s",digis[iD].Data(),cats[iC].Data(),vars[iV].Data()),h);
    if(h == 0) continue;
    ++iH;
    // h->Rebin(4);
    h->SetTitle(digis[iD]);
    p->addHist(h,TString::Format("%s",cats[iC].Data()));
    }
    // false,digis[iD])
        // if(iH)p->drawRatio(0,"s",true,false,digis[iD]);
    p->normalize();
    p->setMinMax(0.001,1);
    p->draw(false,digis[iD] +"_"+ vars[iV]);
    }
}

}

//Table

{
    TFile * f = new TFile("testSegments.root");

        // TString digis[] = {"p8s768","p12s256","p8s384","p6s512","p4s768","","",""};
        // TString digiNs[] = {"8 part. 768 strips","12 part. 256 strips","8 part. 384 strips","6 part. 512 strips","4 part. 768 strips","","",""};

        TString digis[] = {"p6s512","p6s512M1","p6s512M2","p6s512M3","p6s512M4","p6s512M5","p6s512M6","","",""};
        TString digiNs[] = {"6 part. 512 strips","6 part. 512 strips, no lay 1","6 part. 512 strips, no lay 2","6 part. 512 strips, no lay 3","6 part. 512 strips, no lay 4","6 part. 512 strips, no lay 5","6 part. 512 strips, no lay 6","","",""};


        cout <<" \t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            cout << digiNs[iD] <<"\t";
        }
        cout <<endl;
        
        //Do fake rate
        cout << "Fake Rate (mean,sigma)\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            TH1 * h = 0;
            f->GetObject(TString::Format("%s__fake_nSegs",digis[iD].Data()),h);
            float mean = 0;
            float std = 0;
            if(h){
                mean = h->GetMean();
                std  = h->GetStdDev();
            }
                 cout << TString::Format("%.0f,%.0f\t",mean,std);
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
        TString effnumNs[] = {"perfect (%)","lost (%)",""};
        for(unsigned int iP = 0; pts[iP][0]; ++iP){
            for(unsigned int iN = 0; effnums[iN][0]; ++iN){
            cout << ptNs[iP] <<" "<<effnumNs[iN]<< "\t";
        for(unsigned int iD = 0; digis[iD][0]; ++iD){
            TH1 * h1 = 0;
            f->GetObject(TString::Format("%s__%s_dirty_muon_seg_nHits",digis[iD].Data(),pts[iP].Data()),h1);
            TH1 * h2 = 0;
            f->GetObject(TString::Format("%s__%s_perfect_muon_seg_nHits",digis[iD].Data(),pts[iP].Data()),h2);
            TH1 * h3 = 0;
            f->GetObject(TString::Format("%s__%s_lost_muon_seg_nHits",digis[iD].Data(),pts[iP].Data()),h3);

            TH1 * h = 0;
            f->GetObject(TString::Format("%s__%s_%s_muon_seg_nHits",digis[iD].Data(),pts[iP].Data(),effnums[iN].Data()),h);
            float eff = 0;
            float num = 0;
            float tot = 0;
            if(h && h1){
                tot = h1->Integral()+h2->Integral()+h3->Integral();
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
                 cout << TString::Format("%.2e\t",std);
        }
        cout << endl;
    }

}
cout << endl;

}
}
