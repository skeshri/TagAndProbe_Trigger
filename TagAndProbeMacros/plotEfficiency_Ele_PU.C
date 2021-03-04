void plotEfficiency_Ele_PU(TString triggerFile){

TFile *f1 = TFile::Open("DY2018_PU_eff_tightid_nominal.root");
TFile *f2 = TFile::Open("DY2021_PU_eff_tightid_nominal.root");
TFile *f3 = TFile::Open("DY2021_PU_eff_tightid_nominal.root");

     ifstream inFile;
     inFile.open(triggerFile.Data());
     TString trigger_name;

  while(!inFile.eof())
  {

   TH1F *h_total_2018A,*h_total_2018B,*h_total_2018C;
   TH1F *h_pass_2018A,*h_pass_2018B,*h_pass_2018C;
  
   TGraphAsymmErrors *gra_2018A=0;
   TGraphAsymmErrors *gra_2018B=0;
   TGraphAsymmErrors *gra_2018C=0;

   inFile>>trigger_name;
   if(inFile.eof()) break;
   TString total = trigger_name.Data();
   total += "_total";
   cout << total << endl;
   if(total.Contains("_Iso") ) total.ReplaceAll("_Iso","");
   TString pass = trigger_name.Data();
   pass += "_pass"; 

  h_total_2018A = (TH1F*)f1->Get(total.Data());
  h_total_2018B = (TH1F*)f2->Get(total.Data());
  h_total_2018C = (TH1F*)f3->Get(total.Data());


  h_pass_2018A = (TH1F*)f1->Get(pass.Data());
  h_pass_2018B = (TH1F*)f2->Get(pass.Data());
  h_pass_2018C = (TH1F*)f3->Get(pass.Data());

cout << "h_total_2018A = " << h_total_2018A->Integral() << "h_pass_2018A = " << h_pass_2018A->Integral() << endl;
cout << "h_total_2018B = " << h_total_2018B->Integral() << "h_pass_2018B = " << h_pass_2018B->Integral() << endl;
cout << "h_total_2018C = " << h_total_2018C->Integral() << "h_pass_2018C = " << h_pass_2018C->Integral() << endl;

  TMultiGraph *mg = new TMultiGraph();

  gra_2018A = new TGraphAsymmErrors (h_pass_2018A,h_total_2018A);
  gra_2018B = new TGraphAsymmErrors (h_pass_2018B,h_total_2018B);
  gra_2018C = new TGraphAsymmErrors (h_pass_2018C,h_total_2018C);

  mg->Add(gra_2018A);
  mg->Add(gra_2018B);
//  mg->Add(gra_2018C);


  gra_2018A->SetLineColor(1);
  gra_2018B->SetLineColor(2);
  gra_2018C->SetLineColor(3);

  gra_2018A->SetMarkerStyle(20);
  gra_2018B->SetMarkerStyle(20);
  gra_2018C->SetMarkerStyle(20);

  gra_2018A->SetMarkerColor(1);
  gra_2018B->SetMarkerColor(2);
  gra_2018C->SetMarkerColor(3);


  mg->SetTitle(h_total_2018A->GetTitle());
  TCanvas *c = new TCanvas("c","",200,10,600,500);
//  c->SetLogx();
  mg->Draw("APE");
  mg->SetMinimum(0.2);
  mg->SetMaximum(1.02);
  mg->GetXaxis()->CenterTitle(1);
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle(1);
  mg->GetYaxis()->SetTitleOffset(1.4);

   TLegend *legend = new TLegend(0.5,0.3,0.7,0.5);
   legend->AddEntry(gra_2018A,"2018","p");
   legend->AddEntry(gra_2018B,"2021","p");
//   legend->AddEntry(gra_2018C,"2023","p");
   legend->Draw();

  TString pngFileName = trigger_name.Data() ;
  pngFileName +=  "_efficiency.png";
  c->SaveAs(pngFileName.Data());
  c->Clear();
//  gSystem->Exec("mkdir -p results_Ele");
//  TString dirName = "results_Ele/";
//  dirName+=trigger_name.Data();
 // gSystem->Exec("mkdir -p  "+ dirName);
//  gSystem->Exec("mv " + trigger_name+"*.png  "+dirName);

        delete c;

}//inFile while look
}// end of program
