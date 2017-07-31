//Hits per layer
{
  TFile * f = new TFile("bkgTest.root");
  
  TString prefix = "p8s384";
  TString settings[] = {"pu140n5","pu200n7p14","pu250n8p93","pu140n15","pu200n21p43","pu250n26p79",""};
  TString settingNs[] = {"140 PU, x1 neut.","200 PU, x1 neut.","250 PU, x1 neut.","140 PU, x3 neut.","200 PU, x3 neut.","250 PU, x3 neut.",""};
  // TString settings[] = {"pu200n5l","pu200n10","pu200n7p14","pu200n21p43",""};
  vector<double> overflows;
  double ofV = 22;//value at which you start to lose data
  TString var = "numberOfDigis";
  Plotter * p = new Plotter();
  cout <<"Average number of hits per layer"<<endl;
  for(unsigned int iS = 0; settings[iS][0]; ++iS){
    int color =  StyleInfo::getLineColor(iS);
    int linestyle = 1;
    if(iS > 2){
      linestyle = 2;
      color = StyleInfo::getLineColor(iS -3);
    }
    
    if(iS % 3 == 0){
      cout << endl;
    } 
    TH1 * hN = 0;
    f->GetObject(TString::Format("%s%s_nEvents",prefix.Data(),settings[iS].Data()),hN);
    if(hN == 0) continue;
    TH1 * hV = 0;
    f->GetObject(TString::Format("%s%s__%s",prefix.Data(),settings[iS].Data(),var.Data()),hV);
    if(hV == 0) continue; 
    PlotTools::normalize(hV);
    cout << TString::Format("%.2f",hV->GetMean()) <<"\t";
    overflows.push_back(hV->Integral(hV->FindFixBin(ofV),-1));
    p->addHistLine(hV, settingNs[iS],color,linestyle);
    
    
  }
    cout<<endl <<"Overflow >>"<<endl;
  for(unsigned int iS = 0; iS < overflows.size();++iS){
    if(iS % 3 == 0) cout << endl;
    cout << TString::Format("%.2e", overflows[iS]) <<"\t"; 
  }
  cout << endl;
  p->setYTitle("fraction of BXs");
  p->setXTitle("Number of hits per ME0 layer");
  p->draw(false,"something");
}

//VFats ON
{
  TFile * f = new TFile("bkgTest.root");
  
  TString prefix = "p8s384";
  TString settings[] = {"pu140n5","pu200n7p14","pu250n8p93","pu140n15","pu200n21p43","pu250n26p79",""};
  TString settingNs[] = {"140 PU, x1 neut.","200 PU, x1 neut.","250 PU, x1 neut.","140 PU, x3 neut.","200 PU, x3 neut.","250 PU, x3 neut.",""};
  // TString settings[] = {"pu200n5l","pu200n10","pu200n7p14","pu200n21p43",""};
  vector<double> overflows;
  double ofV = 1;//value at which you start to lose data
  TString var = "numberOfVFAT1On";
  Plotter * p = new Plotter();
  cout <<"Average number of vfats on"<<endl;
  for(unsigned int iS = 0; settings[iS][0]; ++iS){
    if(iS % 3 == 0){
      cout << endl;

    } 
    TH1 * hN = 0;
    f->GetObject(TString::Format("%s%s_nEvents",prefix.Data(),settings[iS].Data()),hN);
    if(hN == 0) continue;
    TH1 * hV = 0;
    f->GetObject(TString::Format("%s%s__%s",prefix.Data(),settings[iS].Data(),var.Data()),hV);
    if(hV == 0) continue; 
    PlotTools::normalize(hV);
    cout << TString::Format("%.2f",hV->GetMean()) <<"\t";
    overflows.push_back(hV->Integral(hV->FindFixBin(ofV),-1));
    p->addHistLine(hV, settingNs[iS]);
    
    
  }
    cout <<endl<<">= 1 >>"<<endl;
  for(unsigned int iS = 0; iS < overflows.size();++iS){
    if(iS % 3 == 0) cout << endl;
    cout << TString::Format("%.0f%%", overflows[iS]*100) <<"\t"; 
  }
  cout << endl;
  p->setYTitle("fraction of L1As");
  p->setXTitle("Number of BXs with a hit");
  p->draw(false,"something");
}

//Segments
{
  TFile * f = new TFile("bkgTest.root");
  double bins[] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5};
  int nBins = 11;
  TString prefix = "p8s384";
  TString settings[] = {"pu140n5","pu200n7p14","pu250n8p93","pu140n15","pu200n21p43","pu250n26p79",""};
  TString settingNs[] = {"140 PU, x1 neut.","200 PU, x1 neut.","250 PU, x1 neut.","140 PU, x3 neut.","200 PU, x3 neut.","250 PU, x3 neut.",""};
  // TString settings[] = {"pu200n5l","pu200n10","pu200n7p14","pu200n21p43",""};
  vector<double> overflows;
  double ofV = 9;//value at which you start to lose data
  TString var = "numberOfSectorSegments";
  Plotter * p = new Plotter();
  cout <<"Average number of vfats on"<<endl;
  for(unsigned int iS = 0; settings[iS][0]; ++iS){
    int color =  StyleInfo::getLineColor(iS);
    int linestyle = 1;
    if(iS > 2){
      linestyle = 2;
      color = StyleInfo::getLineColor(iS -3);
    }
    
    if(iS % 3 == 0) cout << endl;
    TH1 * hN = 0;
    f->GetObject(TString::Format("%s%s_nEvents",prefix.Data(),settings[iS].Data()),hN);
    if(hN == 0) continue;
    TH1 * hV = 0;
    f->GetObject(TString::Format("%s%s__%s",prefix.Data(),settings[iS].Data(),var.Data()),hV);
    if(hV == 0) continue; 
    PlotTools::normalize(hV);
    cout << TString::Format("%.2f",hV->GetMean()) <<"\t";
    overflows.push_back(hV->Integral(hV->FindFixBin(ofV),-1));
    // hV = PlotTools::rebin(hV,nBins,bins);
    p->addHistLine(hV, settingNs[iS],color,linestyle);
    
    
  }
    cout <<endl<<">= 1 >>"<<endl;
  for(unsigned int iS = 0; iS < overflows.size();++iS){
    if(iS % 3 == 0) cout << endl;
    cout << TString::Format("%.2e", overflows[iS]) <<"\t"; 
  }
  cout << endl;
  p->setYTitle("fraction of BXs");
  p->setXTitle("Number of segments per CTP7");
  p->draw(false,"something");
}

