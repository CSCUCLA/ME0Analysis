{
    TFile * f = new TFile("plots_testGEMGeo_clean.root");
    TString vars[] = {"deltaY","deltaYoSig",""};

    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iE =  0 ; iE < 8; ++iE){
         TH1 * h;
         f->GetObject(TString::Format("%s_eta_%u",vars[iV].Data(),iE),h);
         p->addHistLine(h,TString::Format("#eta part. %u",iE),-1,1,2);         
         }  
         p->normalize();
         p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));               
    }    
}


{
    TFile * f = new TFile("puStudy.root");
    TString etas[] = {"tracks_gt0p5_abseta","tracks_gt3p0_abseta","tracks_gt5p0_abseta","tracks_gt10_abseta","","tracks_gt20_abseta",""};
    TString etaNames[]= {"p_{T} > 0.5 GeV","p_{T} > 3.0 GeV","p_{T} > 5.0 GeV","p_{T} > 10 GeV","p_{T} > 20 GeV"};
    TH1 * hNorm;
    f->GetObject("nEvents",hNorm);
    double nE = hNorm->GetBinContent(1);
    Plotter * p = new Plotter;
    



     for(unsigned int iE =  0 ; etas[iE][0]; ++iE){
         TH1 * h;
         f->GetObject(etas[iE],h);

         h->Rebin(20);
         for(unsigned int iB = 0; iB <= h->GetNbinsX(); ++iB){
             double count = h->GetBinContent(iB);
             double error = h->GetBinError(iB);
             double width = h->GetBinWidth(iB);
             // double norm = width*2*TMath::Pi()*nE*2;
             double norm = nE;
             // double norm = width*nE*2;
             h->SetBinContent(iB,count/norm);
             h->SetBinError(iB,error/norm);
             h->SetYTitle("d^{2}N/d#phid#eta");
             // h->SetYTitle("dN/d#eta");
             h->SetXTitle("#eta");
         }
         p->addHist(h,etaNames[iE]);
         }

                  p->draw();
         // p->normalize();
         // p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));
}

//PU plot without abs
{
    TFile * f = new TFile("trackDensity.root");
    TString etas[] = {"tracks_gt0p5_eta","tracks_gt3p0_eta","tracks_gt5p0_eta","","tracks_gt10_eta","","tracks_gt20_eta",""};
    TString etaNames[]= {"p_{T} > 0.5 GeV","p_{T} > 3.0 GeV","p_{T} > 5.0 GeV","p_{T} > 10 GeV","p_{T} > 20 GeV"};
    TH1 * hNorm;
    f->GetObject("nEvents",hNorm);
    double nE = hNorm->GetBinContent(1);
    Plotter * p = new Plotter;
    



     for(unsigned int iE =  0 ; etas[iE][0]; ++iE){
         TH1 * h;
         f->GetObject(etas[iE],h);

         h->Rebin(20);
         for(unsigned int iB = 0; iB <= h->GetNbinsX(); ++iB){
             double count = h->GetBinContent(iB);
             double error = h->GetBinError(iB);
             double width = h->GetBinWidth(iB);
             // double norm = width*2*TMath::Pi()*nE*2;
             double norm = nE;
             // double norm = width*nE*2;
             h->SetBinContent(iB,count/norm);
             h->SetBinError(iB,error/norm);
             h->SetYTitle("# of inner tracks from PU per BX / 0.2");
             // h->SetYTitle("dN/d#eta");
             h->SetXTitle("#eta");
         }
         p->addHist(h,etaNames[iE]);
         }

                  p->draw();
         // p->normalize();
         // p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));
}


//Gen DPHI

{
    TFile * f = new TFile("plots.root");

    TString var = "dphi";    
    // TString var = "eta2p6to2p7_dphi";
    // TString var = "eta2p1to2p2_dphi";
    

    TString plots[] = {"gen_pt0p5to1","gen_pt1to3","gen_pt3to5","gen_pt5to10","gen_pt10to20","gen_ptgeq20","","",""};
    TString plotNames[]= {"0.5 < p_{T} <  1 GeV","1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} <  20 GeV","p_{T} > 20 GeV"};
    
    // TString plots[] = {"gen_absdphilt0p001_pt","gen_absdphieq0p001to0p002_pt","gen_absdphieq0p001to0p002_pt","gen_absdphieq0p002to0p005_pt","gen_absdphieq0p005to0p015_pt","gen_absdphigeq0p015_pt","","",""};
    // TString plotNames[]= {"0.5 < p_{T} <  1 GeV","1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} <  20 GeV","p_{T} > 20 GeV"};
    
    
    Plotter * p = new Plotter;
    for(unsigned int iE =  0 ; plots[iE][0]; ++iE){
        TH1 * h;
        f->GetObject(TString::Format("%s_%s",plots[iE].Data(),var.Data()),h);
        h->SetYTitle("a.u.");
        // p->addHist(h,plotNames[iE]);
                p->addHistLine(h,plotNames[iE]);
    }
     p->normalize();    
    p->draw();
    
}

{
    TFile * f = new TFile("plots.root");

    // TString plots[] = {"gen_pt0p5to1_dphi","gen_pt1to3_dphi","gen_pt3to5_dphi","gen_pt5to10_dphi","gen_pt10to20_dphi","gen_ptgeq20_dphi","","",""};    
    // TString plotNames[]= {"0.5 < p_{T} <  1 GeV","1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} <  20 GeV","p_{T} > 20 GeV"};
    
    TString plots[] = {"gen_pt","gen_absdphilt0p001_pt","gen_absdphieq0p001to0p002_pt","gen_absdphieq0p002to0p005_pt","gen_absdphieq0p005to0p008_pt","gen_absdphieq0p008to0p015_pt","gen_absdphigeq0p015_pt","","",""};    
    TString plotNames[]= {"all","|#Delta#phi| < 0.001","0.001 < |#Delta#phi| < 0.002","0.002 < |#Delta#phi| < 0.005","0.005 < |#Delta#phi| < 0.008","0.008 < |#Delta#phi| < 0.015","0.015 < |#Delta#phi|",""};
    
    Plotter * p = new Plotter;
    for(unsigned int iE =  0 ; plots[iE][0]; ++iE){
        TH1 * h;
        f->GetObject(plots[iE],h);
        h->SetYTitle("a.u.");
        // p->addHist(h,plotNames[iE]);
                p->addHistLine(h,plotNames[iE]);
    }
     // p->normalize();
    p->drawRatio(0,"stack",true);
    
}


{
    TFile * f = new TFile("plots.root");
    TString plots[] = {"gen_pt1to3","gen_pt3to5","gen_pt5to10","gen_pt10to20","gen_ptgeq20","","",""};
    TString plotNames[]= {"1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} < 20 GeV","p_{T} > 20 GeV"};
    TString vs[] = {"tg_dpt","tg_dphi","tg_deta","tg_ddphi","tg_ddeta","tg_curv",""};
    
    // TString vs[] = {"rg_dphi","rg_deta","rg_ddphi",""};

    for(unsigned int iV =  0 ; vs[iV][0]; ++iV){    
    Plotter * p = new Plotter;
    for(unsigned int iE =  0 ; plots[iE][0]; ++iE){
        TH1 * h;
        f->GetObject(TString::Format("%s_%s",plots[iE].Data(),vs[iV].Data()) ,h);
        h->SetYTitle("a.u.");
        p->addHist(h,plotNames[iE]);

    }
     p->normalize();    
    p->draw();
}
}

//Effs
{
    TFile * f = new TFile("plots.root");
    // TString plots[] = {"gen_pt1to3","gen_pt3to5","gen_pt5to10","gen_ptgeq10","","",""};    
    // TString plotNames[]= {"1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","p_{T} > 10 GeV"};
    
    // TString plots[] = {"gen_pt3to5","","gen_ptgeq10","","",""};
    TString plotNames[]= {"3 < p_{T} <  5 GeV","","p_{T} > 10 GeV"};
    TString plots[] = {"gen_pt5to10_tight","","gen_ptgeq10_tight","","",""};
    // TString plotNames[]= {"3 < p_{T} <  5 GeV","p_{T} > 10 GeV"};

    // TString vs[] = {"tg_abs_dphi","tg_abs_deta","tg_abs_ddphi",""};
    TString vs[] = {"tg_abs_ddphi",""};
    
    // TString vs[] = {"tg_abs_dphi","tg_abs_deta","tg_abs_ddphi",""};


    for(unsigned int iV =  0 ; vs[iV][0]; ++iV){    
        for(unsigned int iE =  0 ; plots[iE][0]; ++iE){            
        TH1 * hf;
        TH1 * hr;
        f->GetObject(TString::Format("%s_%s",plots[iE].Data(),vs[iV].Data()) ,hr);
        f->GetObject(TString::Format("%s_fake_%s",plots[iE].Data(),vs[iV].Data()) ,hf);
        // PlotTools::toUnderflow(hf);PlotTools::toOverflow(hf);
        // PlotTools::toUnderflow(hr);PlotTools::toOverflow(hr);
        double totNum = hr->Integral(0,-1);
        TH1 * hfInt = (TH1*)hf->Clone();
        TH1 * hrInt = (TH1*)hr->Clone();        
        // for(unsigned int iB = 0; iB <= hf->GetNbinsX(); ++iB ){
        //     hfInt->SetBinContent(iB, hf->Integral(0,iB) );
        //     hrInt->SetBinContent(iB, hr->Integral(0,iB) );
        // }
        // hfInt->SetBinContent(hf->GetNbinsX(), hf->Integral(0,hf->GetNbinsX()+1) );
        // hrInt->SetBinContent(hf->GetNbinsX(), hr->Integral(0,hf->GetNbinsX()+1) );
        // hfInt->SetBinContent(hf->GetNbinsX()+1, 0);
        // hrInt->SetBinContent(hf->GetNbinsX()+1, 0);
        // hfInt->SetBinContent(0, 0);
        // hrInt->SetBinContent(0, 0);
        
        hfInt->Scale(1/totNum);
        hrInt->Scale(1/totNum);
        hfInt->SetTitle(plotNames[iE]);
        hrInt->SetTitle(plotNames[iE]);
        Plotter * p = new Plotter;
        p->addHist(hfInt,"(PU tracks)/(# of muon tracks)");
        p->addHist(hrInt,"Normalized muon tracks");
        p->setMinMax(0.,1.);
        p->draw();
    }
    
}
}



{

    TFile * f = new TFile("plots.root");
    TString plots[] = {"gen_pt1to3","gen_pt3to5","gen_pt5to10","gen_pt10to20","gen_ptgeq20","","",""};
    TString plotNames[]= {"1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} < 20 GeV","p_{T} > 20 GeV"};
    TString v = "tg_abs_dphi_deta";
    double cuts[] = {-1,0.05,0.04,.03,.02,-1};
    
    
    
    for(unsigned int iE =  0 ; plots[iE][0]; ++iE){            
        double nMuons = 0;
    TH2 * hf;
    TH2 * hr;
    f->GetObject(TString::Format("%s_%s",plots[iE].Data(),v.Data()) ,hr);
    f->GetObject(TString::Format("%s_fake_3GeV_%s",plots[iE].Data(),v.Data()) ,hf);
    cout << plotNames[iE] <<"\t";
    for(unsigned int iC = 0; !iC || cuts[iC] >= 0; ++iC){
        double nF = 0;
        double nR = 0;
        double eF = 0;
        double eRu = 0;
        double eRd = 0;
        if(iC){
        int binNumber = hf->GetXaxis()->FindFixBin(cuts[iC]);        
        nF = hf->Integral(0,binNumber,0,binNumber);
        nR = hr->Integral(0,binNumber,0,binNumber);            
        eRu = TEfficiency::ClopperPearson((unsigned int)nMuons, (unsigned int)nR, 0.683, true) - nR/nMuons;
        eRd = nR/nMuons - TEfficiency::ClopperPearson((unsigned int)nMuons, (unsigned int)nR, 0.683, false);        
        } else {
            nF = hf->Integral(0,-1,0,-1);
            nR = hr->Integral(0,-1,0,-1);  
            nMuons = nR;
        }
        eF = TMath::Sqrt(nF);
        // cout << nR << "\t";        
        // cout<< TString::Format("%.1f +%.1f -%.1f\t",100.*nR/nMuons,100.*eRu,100.*eRd);
        cout << TString::Format("%.2f +/- %.2f\t",nF/nMuons,eF/nMuons);
        // cout << nF/nMuons <<"\t";

        
        
    }
    cout << endl;
}


}


{

    TFile * f = new TFile("plots.root");
    // TString plots[] = {"gen_pt1to3","gen_pt3to5","gen_pt5to10","gen_pt10to20","gen_ptgeq20","","",""};
    TString plots[] = {"gen_pt1to3_tight_fine","gen_pt3to5_tight_fine","gen_pt5to10_tight_fine","gen_pt10to20_tight_fine","gen_ptgeq20_tight_fine","","",""};
    TString plotNames[]= {"1 < p_{T} <  3 GeV","3 < p_{T} <  5 GeV","5 < p_{T} <  10 GeV","10 < p_{T} < 20 GeV","p_{T} > 20 GeV"};
    TString v = "tg_abs_ddphi";
    double cuts[] = {-1,3.18,.007,.006,.005,-1};
    
    for(unsigned int iE =  0 ; plots[iE][0]; ++iE){            
        double nMuons = 0;
    TH1 * hf;
    TH1 * hr;
    f->GetObject(TString::Format("%s_%s",plots[iE].Data(),v.Data()) ,hr);
    f->GetObject(TString::Format("%s_fake_3GeV_%s",plots[iE].Data(),v.Data()) ,hf);
    cout << plotNames[iE] <<"\t";
    for(unsigned int iC = 0; !iC || cuts[iC] >= 0; ++iC){
        double nF = 0;
        double nR = 0;
        double eF = 0;
        double eRu = 0;
        double eRd = 0;
        if(iC){
        int binNumber = hf->GetXaxis()->FindFixBin(cuts[iC]);        
        nF = hf->Integral(0,binNumber);
        nR = hr->Integral(0,binNumber);            
        eRu = TEfficiency::ClopperPearson((unsigned int)nMuons, (unsigned int)nR, 0.683, true) - nR/nMuons;
        eRd = nR/nMuons - TEfficiency::ClopperPearson((unsigned int)nMuons, (unsigned int)nR, 0.683, false);   
        } else {
            nF = hf->Integral(0,-1);
            nR = hr->Integral(0,-1);  
            nMuons = nR;
        }
        eF = TMath::Sqrt(nF);
        
        // cout << nR <<"\t";
        // cout << nR/nMuons <<"\t";
        // cout << nF/nMuons <<"\t";
        // cout<< TString::Format("%.1f +%.1f -%.1f\t",100.*nR/nMuons,100.*eRu,100.*eRd);
        cout << TString::Format("%.2f +/- %.2f\t",nF/nMuons,eF/nMuons);

        
        
    }
    cout << endl;
}


}



{
    TFile * f = new TFile("plots.root");

    TH1 * h;
    // f->GetObject("rg_rechit_deta",h);
    f->GetObject("rg_rechit_dphi",h);
    h->SetYTitle("a.u.");
    // h->SetXTitle("rechit position - simhit position");
    p = new Plotter;
    // p->addHist(h,"#eta");
    // TH1 * h2;

    // h->SetYTitle("a.u.");
    // h->SetXTitle("rechit position - simhit position");
    p->addHist(h,"#phi");
    p->normalize();
    p->draw();
  
    
}





/////////////// Fake study;
{
double nE = nEvents->GetBinContent(1);    
TH1 * h = particleBreakdown;
h->Scale(1./nE);
h->Draw();
TString names[] = {"mu from p.i.","mu from decay","other mu","e from photon", "other e","hadron from p.i.","hadron from decay","hadron from scattering","other hadron","other"};
cout << endl;
for(unsigned int iB =0 ; iB < h->GetNbinsX();++iB){
    cout << names[iB] <<"\t";
}
cout << endl;
// for(unsigned int iB =1 ; iB <= h->GetNbinsX();++iB){
//     cout << TString::Format("%.3f +/- %.3f\t",h->GetBinContent(iB),h->GetBinError(iB));
// }
for(unsigned int iB =1 ; iB <= h->GetNbinsX();++iB){
    cout << TString::Format("%.1f +/- %.1f\t",200*h->GetBinContent(iB),200*h->GetBinError(iB));
    // cout << TString::Format("%.1f +/- %.1f\t",h->GetBinContent(iB),h->GetBinError(iB));
}
cout << endl;
}


{
    TFile * f = new TFile("simHitAnalyzer_MinBias.root");
	
	TH1 *nH;     f->GetObject("nEvents",nH);

double nE = nH->GetBinContent(1);    
TString names[] = {"muons","hadrons","electrons"};
TString sampT[] = {"muon","hadron","ele",""};
cout << endl;

vector<TH1*> hs;
for(unsigned int iS = 0; sampT[iS][0]; ++iS){
	TH1 * h;
    f->GetObject(TString::Format("%s_hitLays",sampT[iS].Data()),h);
	h->Scale(1./nE);
	hs.push_back(h);
    cout << names[iS] <<"\t";
}
cout << endl;

for(unsigned int iB = 1; iB <= 6; ++iB){
	for(unsigned int iH = 0; iH < hs.size(); ++iH){
	    // cout << TString::Format("%.1f +/- %.1f\t",200*hs[iH]->GetBinContent(iB),200*hs[iH]->GetBinError(iB));
	    cout << TString::Format("%.2f +/- %.2f\t",200*hs[iH]->GetBinContent(iB)/36,200*hs[iH]->GetBinError(iB)/36);		
		
	}
	cout << endl;
}
cout << endl;
}

{
    TFile * f = new TFile("digiout_NUNU_0.root");
    TString hs[] = {"ele_hitsByRadius","neutron_hitsByRadius","photon_hitsByRadius",""};
    TString hs2[] = {"ele_hits","neutron_hits","photon_hits",""};
    TString hsn[] = {"electrons","neutrons","photons",""};

    Plotter * p = new Plotter;
    
    for(unsigned int iS = 0; hs[iS][0]; ++iS){
        TH1 * h2;
        f->GetObject(hs2[iS],h2);
        double scale = 100.0*9*6;
        cout <<h2 ->Integral()/scale <<" +/- "<< sqrt(h2->Integral())/scale << endl;

        TH1 * h;
        f->GetObject(hs[iS],h);
        h->Rebin(10);
        
        

        h->GetYaxis()->SetTitle("Rate [Hz/cm^{2}]");        
        h->GetXaxis()->SetTitle("radius [cm]");
    for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB){
        double low = h->GetBinLowEdge(iB);
        double up = h->GetBinLowEdge(iB) + h->GetBinWidth(iB);
        double area = 2*TMath::Pi() *(up*up -low*low);
        double timeW = 100.0*25.0e-9*1.0*9.0;
        double newC = h->GetBinContent(iB)/(area*timeW);
        double newE = h->GetBinError(iB)/(area*timeW);
        h->SetBinContent(iB,newC);
        h->SetBinError(iB,newE);
    }
    p->addHist(h,hsn[iS]);
}
    p->draw();
}


{
    TFile * f = new TFile("gedigiout_MB.root"); //NUNU for neutrons
    TString hs[] = {"ge11__hitsByRadius","ge21__hitsByRadius","",""};
    TString hs2[] = {"ge11__hits","ge21__hits","",""};
    TString hsn[] = {"ge11","ge21","",""};  
    
    TH1 * hn;
    f->GetObject("nEvents",hn);
    double nE = hn->GetBinContent(1);    
    nE /= (9.0*200); //Remove for neutrons

    Plotter * p = new Plotter;
    
    for(unsigned int iS = 0; hs[iS][0]; ++iS){
        TH1 * h2;
        f->GetObject(hs2[iS],h2);
        double scale =  1;
        if(iS == 0){
            scale = nE *9.0*2.0*72.0; 
        } else {
            scale = nE *9.0*2.0*36.0;            
        }
        cout <<h2 ->Integral()/scale <<" +/- "<< sqrt(h2->Integral())/scale << endl;
        TH1 * h;
        f->GetObject(hs[iS],h);
        h->Rebin(100);
        
        

        h->GetYaxis()->SetTitle("Rate [Hz/cm^{2}]");        
        h->GetXaxis()->SetTitle("radius [cm]");
    for(unsigned int iB =1; iB <= h->GetNbinsX(); ++iB){
        double low = h->GetBinLowEdge(iB);
        double up = h->GetBinLowEdge(iB) + h->GetBinWidth(iB);
        double area = 2*TMath::Pi() *(up*up -low*low);
        double timeW = nE*25.0e-9*1.0*9.0;
        double newC = h->GetBinContent(iB)/(area*timeW);
        double newE = h->GetBinError(iB)/(area*timeW);
        h->SetBinContent(iB,newC);
        h->SetBinError(iB,newE);
    }
    p->addHist(h,hsn[iS]);
}
    p->draw();
}



{
    TFile * f = new TFile("simHitAnalyzer_MinBias.root");
    TString vars[] = {"trackPT","trackP","nLaysHit","hitLays","nLays_geq3_dPhi","","nLays_geq4_dPhi","nLays_geq5_dPhi","nLays_geq6_dPhi",""};
    TString samps[] = {"hadron","ele","muon",""};
        TString sampNs[] = {"Hadrons","Electrons","Muons",""};
    
        TH1 * hn;
        f->GetObject("nEvents",hn);
        double nE = hn->GetBinContent(1);    
    
    
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iS =  0 ; samps[iS][0]; ++iS){
         TH1 * h;
         f->GetObject(TString::Format("%s_%s",samps[iS].Data(),vars[iV].Data()),h);
         h->GetYaxis()->SetTitle("# per MinBias event");
         h->Scale(1./nE);
         p->addHistLine(h,sampNs[iS]);         
         }  
         p->draw();               
    }    
}

1.4((.139*.0634*27500 + 703)/703)

{
double nE = nEvents->GetBinContent(1);    

hadron_trackPT->Scale(1./nE);
ele_trackPT->Scale(1./nE);
muon_trackPT->Scale(1./nE);

hadron_trackPT->GetYaxis()->SetTitle("<#> per MinBias event");

Plotter * p = new Plotter;
p->addHistLine(hadron_trackPT,"Hadrons");
p->addHistLine(ele_trackPT,"Electrons");
p->addHistLine(muon_trackPT,"Muons");

p->draw();
}

{
double nE = nEvents->GetBinContent(1);    

hadron_trackP->Scale(1./nE);
ele_trackP->Scale(1./nE);
muon_trackP->Scale(1./nE);

hadron_trackP->GetYaxis()->SetTitle("<#> per MinBias event");

Plotter * p = new Plotter;
p->addHistLine(hadron_trackP,"Hadrons");
p->addHistLine(ele_trackP,"Electrons");
p->addHistLine(muon_trackP,"Muons");

p->draw();
}

{
double nE = nEvents->GetBinContent(1);    

hadron_trackP->Scale(1./nE);
ele_trackP->Scale(1./nE);
muon_trackP->Scale(1./nE);

hadron_trackP->GetYaxis()->SetTitle("<#> per MinBias event");

Plotter * p = new Plotter;
p->addHistLine(hadron_trackP,"Hadrons");
p->addHistLine(ele_trackP,"Electrons");
p->addHistLine(muon_trackP,"Muons");

p->draw();
}



/// HZZ

{
    TFile * f = new TFile("plots.root");
    TString vars[] = {"etadist","",""};
    TString samps[] = {"incl","pt_lt3","pt3to5","pt5to10","ptgt10",""};
        TString sampNs[] = {"Inclusive","p_{T} < 3 GeV","3 < p_{T} < 5 GeV","5 < p_{T} < 10 GeV","p_{T} > 10 GeV",""};
              
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iS =  0 ; samps[iS][0]; ++iS){
         TH1 * h;
         f->GetObject(TString::Format("%s_%s",samps[iS].Data(),vars[iV].Data()),h);
         // h->GetYaxis()->SetTitle("# per MinBias event");
         // h->Scale(1./nE);
         p->addHistLine(h,sampNs[iS]);
         }  
         p->draw();               
    }    
}

{
        TFile * f = new TFile("plots.root");
        TH2 * h;
        f->GetObject("eta_pt_dist",h);
        double integral = h->Integral(0,-1,0,-1);
        cout << integral << endl;
        
        for(unsigned int iX = 1; iX <= h->GetNbinsX(); ++iX)
            for(unsigned int iY = 1; iY <= h->GetNbinsY(); ++iY){
                double p = h->GetBinContent(iX,iY)/integral;
                double e = TMath::Sqrt(p*(1-p)/integral);
                h->SetBinContent(iX,iY,p*100);
                h->SetBinError(iX,iY,e*100);
            }
            gStyle->SetPaintTextFormat(".1f%%");
            h->Draw("COLZTEXTE");
            
        
    
}

{
    TFile * f = new TFile("plots.root");
    TString vars[] = {"iso","",""};
    TString pts[] = {"pteq3to5","pteq5to10","pteq10to20","ptgeq20",""};
    TString ptsNs[] = {"3 < p_{T} < 5 GeV","5 < p_{T} < 10 GeV","10 < p_{T} < 20 GeV","p_{T} > 20 GeV",""};
    
    TString samps[] = {"zMuon","bkg",""};
    TString sampNs[] = {"muon","bkg.",""};
    
    TString etas[] = {"etalt1p4","etaeq1p4to2p4","etaeq2p4to2p8",""};
    TString etaNs[] = {"|#eta| < 1.4","1.4 < |#eta| < 2.4","2.4 < |#eta| < 2.8",""};
              
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){

     for(unsigned int iS =  0 ; samps[iS][0]; ++iS){
              for(unsigned int iE =  0 ; etas[iE][0]; ++iE){
                  Plotter * p = new Plotter;
              for(unsigned int iP =  0 ; pts[iP][0]; ++iP){
                  // KEY: TH1F    zMuon_ptlt3_etaeq2p4to2p8_iso;1
                  
         TH1 * h;
         cout << TString::Format("%s_%s_%s_%s",samps[iS].Data(),pts[iP].Data(),etas[iE].Data(),vars[iV].Data()) << endl;
         f->GetObject(TString::Format("%s_%s_%s_%s",samps[iS].Data(),pts[iP].Data(),etas[iE].Data(),vars[iV].Data()),h);
         
         TH1 *hi = (TH1*)h->Clone();
         for(unsigned int iB = 1; iB <= hi->GetNbinsX()+1; ++iB){
             double val = h->Integral(0,iB);
             if(iB == hi->GetNbinsX())val = h->Integral(0,-1);
             if(iB == hi->GetNbinsX()+1)val = 0;
             hi->SetBinContent(iB, val/h->Integral(0,-1));
         }
         // h = hi;
         
         p->addHistLine(h,TString::Format("%s, %s, %s",sampNs[iS].Data(),ptsNs[iP].Data(),etaNs[iE].Data()));         
         }  
         p->normalize();
         p->draw();               
    }    
}
}

}

{
    // Bool_t Divide(const TH1* h1, const TH1* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option = "")     // *MENU*
    // KEY: TH2F    zMuon_pass__2D_pass;1
    // KEY: TH2F    bkg_pass__2D_pass;1
    // KEY: TH2F    other_pass__2D_pass;1
    // KEY: TH2F    bkg_incl__2D_pass;1
    // KEY: TH2F    zMuon_incl__2D_pass;1
    // KEY: TH2F    other_incl__2D_pass;1
    gStyle->SetPaintTextFormat(".3f");
    zMuon_pass_2D_pass->Divide(zMuon_pass_2D_pass,zMuon_incl_2D_pass,1,1,"b");
    zMuon_pass_2D_pass->Draw("COLZTEXTE");
    new TCanvas();
    bkg_pass_2D_pass->Divide(bkg_pass_2D_pass,bkg_incl_2D_pass,1,1,"b");
    bkg_pass_2D_pass->Draw("COLZTEXTE");
    
}


{
    TFile * f = new TFile("zmumu_sampleFakeRate_plots.root");
    // TString etas[] = {"bkg_ptgeq3_incl_eta","bkg_ptgeq3_iso_eta",""};
    // TString etaNames[]= {"p_{T} > 3 GeV","isolated p_{T} > 3 GeV",""};
    TString etas[] = {"bkg_ptgeq3_incl_eta","bkg_ptgeq3_iso_eta","","bkg_ptgeq3_3ISO_incl_eta","bkg_ptgeq3_3ISO_iso_eta",""};
    TString etaNames[]= {"p_{T} > 3 GeV","isolated p_{T} > 3 GeV","p_{T} > 3 GeV (w/ISO)","isolated p_{T} > 3 GeV (w/ISO)",""};
    
    TH1 * hNorm;
    f->GetObject("checkTrackIso_nEvents",hNorm);
    double nE = hNorm->GetBinContent(1);
    
    TH1 * hNorm2;
    f->GetObject("checkTrackIso_3ISO_nEvents",hNorm2);
    double nE2 = hNorm2->GetBinContent(1);
    
    Plotter * p = new Plotter;
    



     for(unsigned int iE =  0 ; etas[iE][0]; ++iE){
         TH1 * h;
         f->GetObject(etas[iE],h);

         h->Rebin(5);
         for(unsigned int iB = 0; iB <= h->GetNbinsX(); ++iB){
             double count = h->GetBinContent(iB);
             double error = h->GetBinError(iB);
             double width = h->GetBinWidth(iB);
             // double norm = width*2*TMath::Pi()*nE*2;
             // double norm = nE;
             double norm = width*nE*2;
             if(iE > 1 ) norm = width*nE2*2;
             h->SetBinContent(iB,count/norm);
             h->SetBinError(iB,error/norm);
             // h->SetYTitle("d^{2}N/d#phid#eta");
             h->SetYTitle("dN/d#eta");
             h->SetXTitle("#eta");
         }
         p->addHist(h,etaNames[iE]);
         }

                  p->draw();
         // p->normalize();
         // p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));
}

{
    Plotter * p = new Plotter;
    nLoose->SetXTitle("# of reconstructed muons");
    p->addHist(nLoose,"Loose ID");
    nTight->SetXTitle("# of reconstructed muons");
    p->addHist(nTight,"Tight ID");
    nTightISO->SetXTitle("# of reconstructed muons");
    p->addHist(nTightISO,"Tight ID + ISO");
    p->draw();
}

{
    vector<TH1*> hists; hists.push_back(nLoose);hists.push_back(nTight);hists.push_back(nTightISO);
    for(iH = 0; iH < hists.size(); ++iH){
        double one = hists[iH]->GetBinContent(4)/hists[iH]->GetBinContent(3);
        double two = (hists[iH]->GetBinContent(5)+hists[iH]->GetBinContent(6))/hists[iH]->GetBinContent(4);
        cout << one <<"  " << two <<" "<< hists[iH]->GetBinError(5)/hists[iH]->GetBinContent(4) <<endl;
    }
}


// factor of 10 plot
{
    TFile * ff = new TFile("simHitAnalyzer.root");
    TH1 * hfa = 0;
    TH1 * hf3= 0; 
    TH1 * hf4= 0;
    TH1 * hf4m= 0;
    TH1 * hf5= 0;
    TH1 * hf6= 0;
    TH1 * hf6m= 0;
    ff->GetObject("all_tr_dPhi",hfa);
    ff->GetObject("all_nLays_geq3_tr_dPhi",hf3);
    ff->GetObject("all_nLays_geq4_tr_dPhi",hf4);
    ff->GetObject("muon_nLays_geq4_tr_dPhi",hf4m);
    ff->GetObject("all_nLays_geq5_tr_dPhi",hf5);
    ff->GetObject("all_nLays_geq6_tr_dPhi",hf6);
    ff->GetObject("muon_nLays_geq6_tr_dPhi",hf6m);

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
    
    
    TFile * fs = new TFile("digiAnalyzer_p8s768.root");
    TH1 * hr1 =0;
    TH1 * hr3=0;
    TH1 * hr5=0;
    TH1 * hr20=0;
    fs->GetObject("p8s768_ptleq3__tr_dPhi",hr1);
    fs->GetObject("p8s768_pteq3to5__tr_dPhi",hr3);
    fs->GetObject("p8s768_pteq5to20__tr_dPhi",hr5);
    fs->GetObject("p8s768_ptgeq20__tr_dPhi",hr20);
    
    hflt->SetXTitle("SimHit #Delta#phi");
    hflt->SetYTitle("a.u.");
    
    hf4nonM->SetXTitle("SimHit #Delta#phi");
    hf4nonM->SetYTitle("a.u.");
    
    p = new Plotter;
    // p->addStackHist(hflt ,"bkg, < 3 layers hit");
    // p->addStackHist(hf34o,"bkg, 3-4 layers hit");
    // p->addStackHist(hf56o,"bkg, 5-6 layers hit");
        // p->addStackHist(hfa,"bkg, 5-6 layers hit");
    
    p->addStackHist(hf4nonM,"non-#mu PU bkg, #geq 4 layers hit");
        p->addStackHist(hf4m,"#mu PU bkg, #geq 4 layers hit");
    
    
        // p->addHistLine(hf4,"bkg, >= 4 layers hit");
        // p->addHistLine(hflt4,"bkg, < 4 layers hit");
    
    // p->addHistLine(hr1 ,"muon, p_{T} 1-3");
    p->addHistLine(hr3 ,"#mu signal, p_{T} 3-5 GeV");
    // p->addHistLine(hr5 ,"muon, p_{T} 5-20");
    // p->addHistLine(hr20,"muon, p_{T} 20-30");
    p->rebin(10);
    p->normalize();

    p->draw();
        
}