#define Scouting_Efficiency_cxx
#include "Scouting_Efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Scouting_Efficiency::Loop(TString output)
{
//   In a ROOT session, you can do:
//      root> .L Scouting_Efficiency.C
//      root> Scouting_Efficiency t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile *file = new TFile(output.Data(),"RECREATE");
   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   std::cout<< "total events : "<<nentries<<endl;
   double eta_bins[11] = {-2.5,-2.1,-1.6,-1.4,-0.8,0,0.8,1.4,1.6,2.1,2.5};
   double pt_bins_Ele30[18] = {0,10,20,25,27,29,30,31,33,33,35,40,45,50,60,100,200,1000};
   double pt_bins_Ele16[20] = {0,10,12,14,15,16,17,19,21,23,25,30,35,40,45,50,60,100,200,1000};
   double pt_bins_Ele12[18] = {0,10,11,12,13,14,15,20,25,30,35,40,45,50,60,100,200,1000};

// HLT Ele30
   TH1F *h_Ele30_pt_total = new TH1F("Ele30_pt_total","Ele30_pt",17,pt_bins_Ele30);
   TH1F *h_Ele30_eta_total = new TH1F("Ele30_eta_total","Ele30_eta",10,eta_bins);
   TH2F *h_Ele30_pt_eta_total = new TH2F("Ele30_pt_eta_total","Ele30_pt_eta",10,eta_bins,17,pt_bins_Ele30);
   TH1F *h_Ele30_pt_pass = new TH1F("Ele30_pt_pass","Ele30_pt",17,pt_bins_Ele30);
   TH1F *h_Ele30_eta_pass = new TH1F("Ele30_eta_pass","Ele30_eta",10,eta_bins);
   TH2F *h_Ele30_pt_eta_pass = new TH2F("Ele30_pt_eta_pass","Ele30_pt_eta",10,eta_bins,17,pt_bins_Ele30);

   h_Ele30_pt_total->Sumw2();
   h_Ele30_eta_total->Sumw2();
   h_Ele30_pt_eta_total->Sumw2();
   h_Ele30_pt_pass->Sumw2();
   h_Ele30_eta_pass->Sumw2();
   h_Ele30_pt_eta_pass->Sumw2();

// HLT Ele16_Ele12 Ele16 leg
   TH1F *h_Ele16_pt_total = new TH1F("Ele16_pt_total","Ele16_pt",19,pt_bins_Ele16);
   TH1F *h_Ele16_eta_total = new TH1F("Ele16_eta_total","Ele16_eta",10,eta_bins);
   TH2F *h_Ele16_pt_eta_total = new TH2F("Ele16_pt_eta_total","Ele16_pt_eta",10,eta_bins,19,pt_bins_Ele16);
   TH1F *h_Ele16_pt_pass = new TH1F("Ele16_pt_pass","Ele16_pt",19,pt_bins_Ele16);
   TH1F *h_Ele16_eta_pass = new TH1F("Ele16_eta_pass","Ele16_eta",10,eta_bins);
   TH2F *h_Ele16_pt_eta_pass = new TH2F("Ele16_pt_eta_pass","Ele16_pt_eta",10,eta_bins,19,pt_bins_Ele16);

   h_Ele16_pt_total->Sumw2();
   h_Ele16_eta_total->Sumw2();
   h_Ele16_pt_eta_total->Sumw2();
   h_Ele16_pt_pass->Sumw2();
   h_Ele16_eta_pass->Sumw2();
   h_Ele16_pt_eta_pass->Sumw2();

// HLT Ele16_Ele12 Ele12 leg
   TH1F *h_Ele12_pt_total = new TH1F("Ele12_pt_total","Ele12_pt",17,pt_bins_Ele12);
   TH1F *h_Ele12_eta_total = new TH1F("Ele12_eta_total","Ele12_eta",10,eta_bins);
   TH2F *h_Ele12_pt_eta_total = new TH2F("Ele12_pt_eta_total","Ele12_pt_eta",10,eta_bins,17,pt_bins_Ele12);
   TH1F *h_Ele12_pt_pass = new TH1F("Ele12_pt_pass","Ele12_pt",16,pt_bins_Ele12);
   TH1F *h_Ele12_eta_pass = new TH1F("Ele12_eta_pass","Ele12_eta",10,eta_bins);
   TH2F *h_Ele12_pt_eta_pass = new TH2F("Ele12_pt_eta_pass","Ele12_pt_eta",10,eta_bins,16,pt_bins_Ele12);

   h_Ele12_pt_total->Sumw2();
   h_Ele12_eta_total->Sumw2();
   h_Ele12_pt_eta_total->Sumw2();
   h_Ele12_pt_pass->Sumw2();
   h_Ele12_eta_pass->Sumw2();
   h_Ele12_pt_eta_pass->Sumw2();


// HLT Ele16_Ele12
   TH2F *h_Ele16_Ele12_pt_total = new TH2F("Ele_16_Ele12_pt_total","Ele16_Ele12_pt",19,pt_bins_Ele16,17,pt_bins_Ele12);
   TH2F *h_Ele16_Ele12_pt_pass = new TH2F("Ele_16_Ele12_pt_pass","Ele16_Ele12_pt",19,pt_bins_Ele16,17,pt_bins_Ele12);
   TH2F *h_Ele16_Ele12_eta_total = new TH2F("Ele_16_Ele12_eta_total","Ele16_Ele12_eta",10,eta_bins,10,eta_bins);
   TH2F *h_Ele16_Ele12_eta_pass = new TH2F("Ele_16_Ele12_eta_pass","Ele16_Ele12_eta",10,eta_bins,10,eta_bins);

   h_Ele16_Ele12_pt_total->Sumw2();
   h_Ele16_Ele12_pt_pass->Sumw2();
   h_Ele16_Ele12_eta_total->Sumw2();
   h_Ele16_Ele12_eta_pass->Sumw2();

   Long64_t nbytes = 0, nb = 0;
   int pass_nEle=0;
   int pass_RefTrig=0;
   int pass_TargetTrig=0;
   int pass_Id = 0;
   int pass_Id_RefTrig = 0;
   int pass_denominator=0;
   int pass_filter=0;
   int pass_numerator=0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%50000 == 0)cout<<"number of events processed : "<<jentry<<endl;
      if(nEle!=2) continue;
      pass_nEle++;

      bool denominator = false;
      bool numerator = false;
      bool passEle1Id = passEleIdLoose->at(0);
      bool passEle2Id = passEleIdLoose->at(1);

      bool passRefTrig = false;
      bool passTargetTrig = false;
      int i_path=0;
      for(auto p: *triggerPath)
      {
          TString path(p);
          if (path.Contains("DST_HLTMuon_Run3_PFScoutingPixelTracking") && triggerDecision->at(i_path)) 
          //if (path.Contains("DST_Run3_JetHT_PFScoutingPixelTracking") && triggerDecision->at(i_path)) 
          {
              passRefTrig = true;
              pass_RefTrig++;
          //    if (passEle1Id && passEle2Id) pass_Id_RefTrig++;
          }
          if (path.Contains("DST_Run3_EG16_EG12_PFScoutingPixelTracking") && triggerDecision->at(i_path)) 
          {
              passTargetTrig = true;
              pass_TargetTrig++;
          }

          i_path++;
      }
      if (passEle1Id && passEle2Id) pass_Id++;

      denominator = passEle1Id && passEle2Id && passRefTrig;
      //denominator = passRefTrig;
      
      bool matchEle1Fil = passFilterDSTEG16EG12->at(0);
      bool matchEle2Fil = passFilterDSTEG16EG12->at(1);

      numerator = denominator && matchEle1Fil && matchEle2Fil && passTargetTrig;
      if (matchEle1Fil && matchEle2Fil) pass_filter++;
      //numerator = denominator && passTargetTrig;
      
      if (denominator)
      {
          pass_denominator++;
          h_Ele16_pt_total->Fill(ele_pt->at(0));
          h_Ele16_eta_total->Fill(ele_eta->at(0));
          h_Ele16_pt_eta_total->Fill(ele_eta->at(0),ele_pt->at(0));
          h_Ele12_pt_total->Fill(ele_pt->at(1));
          h_Ele12_eta_total->Fill(ele_eta->at(1));
          h_Ele12_pt_eta_total->Fill(ele_eta->at(1),ele_pt->at(1));
          h_Ele16_Ele12_pt_total->Fill(ele_pt->at(0),ele_pt->at(1));
          h_Ele16_Ele12_eta_total->Fill(ele_eta->at(0),ele_eta->at(1));
      }
      
      if (numerator)
      {
          pass_numerator++;
          h_Ele16_pt_pass->Fill(ele_pt->at(0));
          h_Ele16_eta_pass->Fill(ele_eta->at(0));
          h_Ele16_pt_eta_pass->Fill(ele_eta->at(0),ele_pt->at(0));
          h_Ele12_pt_pass->Fill(ele_pt->at(1));
          h_Ele12_eta_pass->Fill(ele_eta->at(1));
          h_Ele12_pt_eta_pass->Fill(ele_eta->at(1),ele_pt->at(1));
          h_Ele16_Ele12_pt_pass->Fill(ele_pt->at(0),ele_pt->at(1));
          h_Ele16_Ele12_eta_pass->Fill(ele_eta->at(0),ele_eta->at(1));
      }


   }
   file->Write();
   std::cout<<"pass_nEle: "<<pass_nEle<<std::endl;
   std::cout<<"pass_RefTrig: "<<pass_RefTrig<<std::endl;
   std::cout<<"pass_Id: "<<pass_Id<<std::endl;
   std::cout<<"pass_Id_RefTrig: "<<pass_Id_RefTrig<<std::endl;
   std::cout<<"pass_denominator: "<<pass_denominator<<std::endl;
   std::cout<<"pass_TargetTrig: "<<pass_TargetTrig<<std::endl;
   std::cout<<"pass_filter: "<<pass_filter<<std::endl;
   std::cout<<"pass_numerator: "<<pass_numerator<<std::endl;
}
