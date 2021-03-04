void plotEfficiency_multi_Ele(TString triggerFile){

//TFile *f1 = TFile::Open("DY2018_PT_eff_tightid_nominal.root");
//TFile *f2 = TFile::Open("DY2021_PT_eff_tightid_nominal.root");


//TFile *f1 = TFile::Open("DY2018_ETA_eff_tightid_nominal.root");
//TFile *f2 = TFile::Open("DY2021_ETA_eff_tightid_nominal.root");

//TFile *f1 = TFile::Open("DY2018_PU_eff_tightid_nominal.root");
//TFile *f2 = TFile::Open("DY2021_PU_eff_tightid_nominal.root");

TFile *f1 = TFile::Open("Zprime2018_efficiency_tightid_nominal.root");
TFile *f2 = TFile::Open("Zprime2021_efficiency_tightid_nominal.root");

     ifstream inFile;
     inFile.open(triggerFile.Data());
     TString trigger_name;

  while(!inFile.eof())
  {

   TH1F *h_total_2018A,*h_total_2018B,*h_total_2018C,*h_total_2018D;
   TH1F *h_pass_2018A,*h_pass_2018B,*h_pass_2018C,*h_pass_2018D;
  
   TGraphAsymmErrors *gra_2018A=0;
   TGraphAsymmErrors *gra_2018B=0;

   inFile>>trigger_name;
   if(inFile.eof()) break;
   TString total = trigger_name.Data();
   total += "_total";
   if(total.Contains("_Iso") ) total.ReplaceAll("_Iso","");
   TString pass = trigger_name.Data();
   pass += "_pass"; 

  h_total_2018A = (TH1F*)f1->Get(total.Data());
  h_total_2018B = (TH1F*)f2->Get(total.Data());


  h_pass_2018A = (TH1F*)f1->Get(pass.Data());
  h_pass_2018B = (TH1F*)f2->Get(pass.Data());

for(int i = 1; i < h_pass_2018A->GetNbinsX(); i ++) {
cout << "total =  " << h_total_2018A->GetBinContent(i) << "  " << "pass = " <<  h_pass_2018A->GetBinContent(i) << endl;
}

  TMultiGraph *mg = new TMultiGraph();

  gra_2018A = new TGraphAsymmErrors (h_pass_2018A,h_total_2018A);
  gra_2018B = new TGraphAsymmErrors (h_pass_2018B,h_total_2018B);

  mg->Add(gra_2018A);
  mg->Add(gra_2018B);
//  mg->Add(gra_2018C);
//  mg->Add(gra_2018D);


  gra_2018A->SetLineColor(1);
  gra_2018B->SetLineColor(2);

  gra_2018A->SetMarkerStyle(20);
  gra_2018B->SetMarkerStyle(20);

  gra_2018A->SetMarkerColor(1);
  gra_2018B->SetMarkerColor(2);


  mg->SetTitle(h_total_2018A->GetTitle());
  TCanvas *c = new TCanvas("c","",200,10,600,500);
//  c->SetLogx();
  TPad *pad1 = new TPad("pad1","pad1",0,0.2,1,1);
   pad1->SetBottomMargin(0.01);
   pad1->SetLogy(0);
   pad1->Draw();
   pad1->cd();
  
  mg->Draw("APE");
  mg->SetMinimum(0.0);
  mg->SetMaximum(1.05);
  mg->GetXaxis()->CenterTitle(1);
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle(1);
  mg->GetYaxis()->SetTitleOffset(1.4);

   TLegend *legend = new TLegend(0.5,0.3,0.7,0.5);
   legend->AddEntry(gra_2018A,"2018","p");
   legend->AddEntry(gra_2018B,"2021","p");
   legend->Draw();

  c->cd();
TH1F *h1 = (TH1F *)h_pass_2018A->Clone();
TH1F *h2 = (TH1F *)h_total_2018A->Clone();
   
h1->Divide(h2);
   
TH1F *h3 = (TH1F *)h_pass_2018B->Clone();
TH1F *h4 = (TH1F *)h_total_2018B->Clone();
  
h3->Divide(h4);
  
h3->Divide(h1);
  
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.2);
  pad2->SetGridy(1);
  pad2->SetPad(0,0.0,1.0,0.2);
  pad2->SetTopMargin(0.);
  pad2->SetBottomMargin(0.5);
  pad2->Draw();
  pad2->cd();
  float yscale = (1.0-0.2)/(0.18-0);

  h3->SetMarkerStyle(21);
  h3->SetMarkerSize(0.8);
  h3->SetStats(0);
  h3->SetMinimum(0.7);
  h3->SetMaximum(1.1);
  h3->SetTitle("");
  h3->GetXaxis()->SetTitle("");
  h3->GetYaxis()->SetTitle("2021/2018");
  h3->GetXaxis()->SetTitleOffset(1.3);
  h3->GetXaxis()->SetLabelSize(0.033*yscale);
  h3->GetXaxis()->SetTitleSize(0.036*yscale);
  h3->GetXaxis()->SetTickLength(0.03*yscale);
  h3->GetYaxis()->SetTitleOffset(0.3);
  h3->GetYaxis()->SetNdivisions(5);
  h3->GetYaxis()->SetLabelSize(0.020*yscale);
  h3->GetYaxis()->SetTitleSize(0.036*yscale);
  h3->Draw("");


  TString pngFileName = trigger_name.Data() ;
  pngFileName +=  "_efficiency_highPT.png";
  c->SaveAs(pngFileName.Data());
  c->Clear();

        delete c;

}//inFile while look
}// end of program
