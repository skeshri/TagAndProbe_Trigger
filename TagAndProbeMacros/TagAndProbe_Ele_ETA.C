#define TagAndProbe_cxx
#include "TagAndProbe.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TagAndProbe::Loop(TString output)
{
   if (fChain == 0) return;
  bool RunSystematic=false;
   vector<TString> systematicVar =  {"nominal"};
   
   if (RunSystematic)
   {
	systematicVar.push_back("TagPt_up");
	systematicVar.push_back("TagPt_down");
	systematicVar.push_back("Zmass_up");
	systematicVar.push_back("Zmass_down");
   }

for(int i=0; i<systematicVar.size();i++)
{

   double ptTag = 37; 
   double zMassL = 60; 
   double zMassR = 120;

   if(systematicVar.at(i).Contains("TagPt_up")) ptTag = 41;
   if(systematicVar.at(i).Contains("TagPt_down")) ptTag = 33;
   if(systematicVar.at(i).Contains("Zmass_up")) {zMassL = 50; zMassR = 130;}
   if(systematicVar.at(i).Contains("Zmass_down")) {zMassL = 70; zMassR = 110;}
 
   cout<<systematicVar.at(i).Data()<<endl;
   cout<<" ptTag : "<<ptTag<<" , zMassL : "<<zMassL<<" , zMassR : "<<zMassR<<endl;
   output +="_ETA_eff_";
   output +="tightid_";
   output += systematicVar.at(i);
   output +=".root";

   std::cout<< "total events : "<<fChain->GetEntries()<<endl;
   TFile *file = new TFile(output.Data(),"RECREATE");
   Long64_t nentries = fChain->GetEntriesFast();
//   std::cout<< "total events : "<<nentries<<endl;
   double eta_bins[19] = {-2.5,-2.4,-2.3,-2.2,-2.1,-1.566,-1.4442,-0.8,-0.4,0,0.4,0.8,1.4442,1.566,2.1,2.2,2.3,2.4,2.5};
   double eta_bins_115[11] = {-2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5};


// HLT Ele25

   TH1F *h_Ele25_eta_total = new TH1F("Ele25_eta_total","Ele25_eta",18,eta_bins);
   TH1F *h_Ele25_eta_pass = new TH1F("Ele25_eta_pass","Ele25_eta",18,eta_bins);

   h_Ele25_eta_total->Sumw2();
   h_Ele25_eta_pass->Sumw2();

// HLT Ele27

   TH1F *h_Ele27_eta_total = new TH1F("Ele27_eta_total","Ele27_eta",18,eta_bins);
   TH1F *h_Ele27_eta_pass = new TH1F("Ele27_eta_pass","Ele27_eta",18,eta_bins);

   h_Ele27_eta_total->Sumw2();
   h_Ele27_eta_pass->Sumw2();

// HLT Ele32

   TH1F *h_Ele32_eta_total = new TH1F("Ele32_eta_total","Ele32_eta",18,eta_bins);
   TH1F *h_Ele32_eta_pass = new TH1F("Ele32_eta_pass","Ele32_eta",18,eta_bins);

   h_Ele32_eta_total->Sumw2();
   h_Ele32_eta_pass->Sumw2();

// HLT Ele50

   TH1F *h_Ele50_eta_total = new TH1F("Ele50_eta_total","Ele50_eta",18,eta_bins);
   TH1F *h_Ele50_eta_pass = new TH1F("Ele50_eta_pass","Ele50_eta",18,eta_bins);

   h_Ele50_eta_total->Sumw2();
   h_Ele50_eta_pass->Sumw2();


// HLT Ele115

   TH1F *h_Ele115_eta_total = new TH1F("Ele115_eta_total","Ele115_eta",10,eta_bins_115);
   TH1F *h_Ele115_eta_pass = new TH1F("Ele115_eta_pass","Ele115_eta",10,eta_bins_115);

   h_Ele115_eta_total->Sumw2();
   h_Ele115_eta_pass->Sumw2();

// HLT Ele23_Ele12 Ele23 leg
   TH1F *h_Ele23_Ele12_leg1_eta_total = new TH1F("Ele23_Ele12_leg1_eta_total","Ele23_Ele12_leg1_eta",18,eta_bins);
   TH1F *h_Ele23_Ele12_leg1_eta_pass = new TH1F("Ele23_Ele12_leg1_eta_pass","Ele23_Ele12_leg1_eta",18,eta_bins);
      
   h_Ele23_Ele12_leg1_eta_total->Sumw2();
   h_Ele23_Ele12_leg1_eta_pass->Sumw2();

// HLT Ele23_Ele12 Ele12 leg
   TH1F *h_Ele23_Ele12_leg2_eta_total = new TH1F("Ele23_Ele12_leg2_eta_total","Ele23_Ele12_leg2_eta",18,eta_bins);
   TH1F *h_Ele23_Ele12_leg2_eta_pass = new TH1F("Ele23_Ele12_leg2_eta_pass","Ele23_Ele12_leg2_eta",18,eta_bins);

   h_Ele23_Ele12_leg2_eta_total->Sumw2();
   h_Ele23_Ele12_leg2_eta_pass->Sumw2();


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%50000 == 0)cout<<"number of events processed : "<<jentry<<endl;
      if(nEle!=2) continue;
      int first  = rand()%2;
      int second = (first+1)%2;      
      if(ele_charge->at(first) * ele_charge->at(second)>0)continue;
      bool tag_EleId = passEleIdTight->at(first);
      bool tag_EleKin = ele_pt->at(first)>ptTag && fabs(ele_etaSC->at(first))<2.5;
      bool tag_TriggerMatch = passFilterEle32->at(first);

      if(!(tag_EleId && tag_EleKin && tag_TriggerMatch))continue;

      bool probe_EleId = HWW_Electron_NewDef(second, ele_etaSC->at(second));
      bool probe_EleKin = fabs(ele_etaSC->at(second))<2.5;

      if(!(probe_EleId && probe_EleKin)) continue;

      TLorentzVector tag_eleLV, probe_eleLV, Z_candLV;
      tag_eleLV.SetPtEtaPhiM(ele_pt->at(first), ele_etaSC->at(first), ele_phi->at(first),0.);
      probe_eleLV.SetPtEtaPhiM(ele_pt->at(second), ele_etaSC->at(second), ele_phi->at(second),0.);
      Z_candLV = tag_eleLV + probe_eleLV;

      if (Z_candLV.M()<zMassL || Z_candLV.M() > zMassR) continue;

      if(!(genElectron_fromZ->at(0) == 1 && genElectron_fromZ->at(1) == 1)) continue;


      if(ele_pt->at(second)>33)h_Ele25_eta_total->Fill(ele_etaSC->at(second));
      if(ele_pt->at(second)>35)h_Ele27_eta_total->Fill(ele_etaSC->at(second));
      if(ele_pt->at(second)>40)h_Ele32_eta_total->Fill(ele_etaSC->at(second));
      if(ele_pt->at(second)>58)h_Ele50_eta_total->Fill(ele_etaSC->at(second)); 
      if(ele_pt->at(second)>123)h_Ele115_eta_total->Fill(ele_etaSC->at(second));
      if(ele_pt->at(second)>25)h_Ele23_Ele12_leg1_eta_total->Fill(ele_etaSC->at(second)); 
      if(ele_pt->at(second)>15)h_Ele23_Ele12_leg2_eta_total->Fill(ele_etaSC->at(second)); 


      if (passFilterEle25->at(second)){
      if(ele_pt->at(second)>33)h_Ele25_eta_pass->Fill(ele_etaSC->at(second));
      }

      if (passFilterEle27->at(second)){
      if(ele_pt->at(second)>35)h_Ele27_eta_pass->Fill(ele_etaSC->at(second));
      }

      if (passFilterEle32->at(second)){
      if(ele_pt->at(second)>40)h_Ele32_eta_pass->Fill(ele_etaSC->at(second));
      }

      if (passFilterEle50->at(second)){
      if(ele_pt->at(second)>58)h_Ele50_eta_pass->Fill(ele_etaSC->at(second)); 
      }
 
      if (passFilterEle115->at(second)){
      if(ele_pt->at(second)>123)h_Ele115_eta_pass->Fill(ele_etaSC->at(second));
      }

      if (passFilterEle23_12_leg1->at(second)){
      if(ele_pt->at(second)>25)h_Ele23_Ele12_leg1_eta_pass->Fill(ele_etaSC->at(second)); 
      } 
      if (passFilterEle23_12_leg2->at(second)){
      if(ele_pt->at(second)>15)h_Ele23_Ele12_leg2_eta_pass->Fill(ele_etaSC->at(second)); 
      } 

      // if (Cut(ientry) < 0) continue;
    }
file->Write(); 
  }
 
}

bool TagAndProbe::HWW_Electron_NewDef(int i, double eta)
{
bool mvaid = false;
bool tightid = false;
bool mediumid= false;
double sieie = ele_full5x5_sigmaIetaIeta->at(i);
double eInvMinusPInv = ele_ooEmooP->at(i);
bool convVeto = ele_passConversionVeto->at(i);
double relIso = ele_relCombIsoWithEA->at(i);
mvaid = passMVAIsoWP90->at(i);
tightid = passEleIdTight->at(i);
mediumid = passEleIdMedium->at(i);

if(fabs(eta) < 1.479){    // Barrel
//if (!convVeto) return false;
//if (relIso >= 0.06) return false;
//if (!mvaid) return false;
if (!tightid) return false;
//if (!mediumid) return false;
if (fabs(ele_d0->at(i))>=0.05) return false;
if (fabs(ele_dz->at(i))>=0.1) return false;
}

else{
//if (sieie >= 0.03) return false;
//if (fabs(eInvMinusPInv) >= 0.014) return false;
//if (!convVeto) return false;
//if (relIso >= 0.06) return false;
//if (!mvaid) return false;
if (!tightid) return false;
//if (!mediumid) return false;
if (fabs(ele_d0->at(i))>=0.1) return false;
if (fabs(ele_dz->at(i))>=0.2) return false;
}
return true;
}
