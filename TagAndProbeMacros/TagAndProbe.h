//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  2 12:51:27 2021 by ROOT version 6.14/09
// from TTree EventTree/Event data
// found on file: /eos/cms/store/group/phys_egamma/arun/TriggerEff_Run3/DYEE_Run3_v1/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/DY2021/210301_180846/0000/TnP_ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef TagAndProbe_h
#define TagAndProbe_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "iostream"
#include "fstream"
#include "cstdlib"

using namespace std;
class TagAndProbe {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           pvNTracks;
   Int_t           good_vertices;
   Int_t           nPV;
   Int_t           nPU;
   Int_t           nPUTrue;
   Float_t         rho;
   Float_t         genWeight;
   UInt_t          genParticles_n;
   vector<double>  *genElectron_pt;
   vector<double>  *genElectron_eta;
   vector<double>  *genElectron_phi;
   vector<double>  *genElectron_energy;
   vector<bool>    *genElectron_fromZ;
   Int_t           nEle;
   vector<double>  *ele_pt;
   vector<double>  *ele_eta;
   vector<double>  *ele_etaSC;
   vector<double>  *ele_phi;
   vector<double>  *ele_tricharge;
   vector<double>  *ele_phiSC;
   vector<double>  *ele_energy;
   vector<double>  *ele_energySC;
   vector<double>  *ele_charge;
   vector<double>  *ele_dEtaIn;
   vector<double>  *ele_dEtaSeed;
   vector<double>  *ele_dPhiIn;
   vector<double>  *ele_hOverE;
   vector<double>  *ele_full5x5_sigmaIetaIeta;
   vector<double>  *ele_isoChargedHadrons;
   vector<double>  *ele_isoNeutralHadrons;
   vector<double>  *ele_isoPhotons;
   vector<double>  *ele_relCombIsoWithEA;
   vector<double>  *ele_isoChargedFromPU;
   vector<double>  *ele_ooEmooP;
   vector<double>  *ele_dr03TkSumPt;
   vector<double>  *ele_dr03EcalRecHitSumEt;
   vector<double>  *ele_dr03HcalDepth1TowerSumEt;
   vector<double>  *ele_d0;
   vector<double>  *ele_dz;
   vector<double>  *ele_SIP;
   vector<int>     *ele_expectedMissingInnerHits;
   vector<int>     *ele_passConversionVeto;
   vector<bool>    *passEleIdLoose;
   vector<bool>    *passEleIdMedium;
   vector<bool>    *passEleIdTight;
   vector<bool>    *passMVAnoIsoWP90;
   vector<bool>    *passMVAnoIsoWP80;
   vector<bool>    *passMVAIsoWP90;
   vector<bool>    *passMVAIsoWP80;
   vector<bool>    *hasMatchedToZ;
   vector<bool>    *passL1EG10;
   vector<bool>    *passL1EG17;
   vector<bool>    *passL1EG23;
   vector<bool>    *passL1EG23Iso;
   vector<bool>    *passL1EG20Iso;
   vector<string>  *triggerPath;
   vector<bool>    *triggerDecision;
   vector<bool>    *passFilterEle32;
   vector<bool>    *passFilterEle23_12_leg1;
   vector<bool>    *passFilterEle23_12_leg2;
   vector<bool>    *passFilterEle115;
   vector<bool>    *passFilterEle50;
   vector<bool>    *passFilterEle25;
   vector<bool>    *passFilterEle27;
   vector<bool>    *passFilterMu12_Ele23_legEle;
   vector<bool>    *passFilterMu23_Ele12_legEle;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_pvNTracks;   //!
   TBranch        *b_good_vertices;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_nPUTrue;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genParticles_n;   //!
   TBranch        *b_genElectron_pt;   //!
   TBranch        *b_genElectron_eta;   //!
   TBranch        *b_genElectron_phi;   //!
   TBranch        *b_genElectron_energy;   //!
   TBranch        *b_genElectron_fromZ;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_ele_pt;   //!
   TBranch        *b_ele_eta;   //!
   TBranch        *b_ele_etaSC;   //!
   TBranch        *b_ele_phi;   //!
   TBranch        *b_ele_tricharge;   //!
   TBranch        *b_ele_phiSC;   //!
   TBranch        *b_ele_energy;   //!
   TBranch        *b_ele_energySC;   //!
   TBranch        *b_ele_charge;   //!
   TBranch        *b_ele_dEtaIn;   //!
   TBranch        *b_ele_dEtaSeed;   //!
   TBranch        *b_ele_dPhiIn;   //!
   TBranch        *b_ele_hOverE;   //!
   TBranch        *b_ele_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_ele_isoChargedHadrons;   //!
   TBranch        *b_ele_isoNeutralHadrons;   //!
   TBranch        *b_ele_isoPhotons;   //!
   TBranch        *b_ele_relCombIsoWithEA;   //!
   TBranch        *b_ele_isoChargedFromPU;   //!
   TBranch        *b_ele_ooEmooP;   //!
   TBranch        *b_ele_dr03TkSumPt;   //!
   TBranch        *b_ele_dr03EcalRecHitSumEt;   //!
   TBranch        *b_ele_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_ele_d0;   //!
   TBranch        *b_ele_dz;   //!
   TBranch        *b_ele_SIP;   //!
   TBranch        *b_ele_expectedMissingInnerHits;   //!
   TBranch        *b_ele_passConversionVeto;   //!
   TBranch        *b_passEleIdLoose;   //!
   TBranch        *b_passEleIdMedium;   //!
   TBranch        *b_passEleIdTight;   //!
   TBranch        *b_passMVAnoIsoWP90;   //!
   TBranch        *b_passMVAnoIsoWP80;   //!
   TBranch        *b_passMVAIsoWP90;   //!
   TBranch        *b_passMVAIsoWP80;   //!
   TBranch        *b_hasMatchedToZ;   //!
   TBranch        *b_passL1EG10;   //!
   TBranch        *b_passL1EG17;   //!
   TBranch        *b_passL1EG23;   //!
   TBranch        *b_passL1EG23Iso;   //!
   TBranch        *b_passL1EG20Iso;   //!
   TBranch        *b_triggerPath;   //!
   TBranch        *b_triggerDecision;   //!
   TBranch        *b_passFilterEle32;   //!
   TBranch        *b_passFilterEle23_12_leg1;   //!
   TBranch        *b_passFilterEle23_12_leg2;   //!
   TBranch        *b_passFilterEle115;   //!
   TBranch        *b_passFilterEle50;   //!
   TBranch        *b_passFilterEle25;   //!
   TBranch        *b_passFilterEle27;   //!
   TBranch        *b_passFilterMu12_Ele23_legEle;   //!
   TBranch        *b_passFilterMu23_Ele12_legEle;   //!

   TagAndProbe(TTree *tree=0);
   virtual ~TagAndProbe();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual void     Loop(TString output);      
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool    HWW_Electron_NewDef(int i, double eta);

};

#endif

#ifdef TagAndProbe_cxx
TagAndProbe::TagAndProbe(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_egamma/arun/TriggerEff_Run3/DYEE_Run3_v1/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/DY2021/210301_180846/0000/TnP_ntuple_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/group/phys_egamma/arun/TriggerEff_Run3/DYEE_Run3_v1/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/DY2021/210301_180846/0000/TnP_ntuple_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/cms/store/group/phys_egamma/arun/TriggerEff_Run3/DYEE_Run3_v1/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/DY2021/210301_180846/0000/TnP_ntuple_1.root:/ntupler");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

TagAndProbe::~TagAndProbe()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TagAndProbe::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TagAndProbe::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TagAndProbe::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ele_pt = 0;
   genElectron_pt = 0;
   genElectron_eta = 0;
   genElectron_phi = 0;
   genElectron_energy = 0;
   genElectron_fromZ = 0;
   ele_eta = 0;
   ele_etaSC = 0;
   ele_phi = 0;
   ele_tricharge = 0;
   ele_phiSC = 0;
   ele_energy = 0;
   ele_energySC = 0;
   ele_charge = 0;
   ele_dEtaIn = 0;
   ele_dEtaSeed = 0;
   ele_dPhiIn = 0;
   ele_hOverE = 0;
   ele_full5x5_sigmaIetaIeta = 0;
   ele_isoChargedHadrons = 0;
   ele_isoNeutralHadrons = 0;
   ele_isoPhotons = 0;
   ele_relCombIsoWithEA = 0;
   ele_isoChargedFromPU = 0;
   ele_ooEmooP = 0;
   ele_dr03TkSumPt = 0;
   ele_dr03EcalRecHitSumEt = 0;
   ele_dr03HcalDepth1TowerSumEt = 0;
   ele_d0 = 0;
   ele_dz = 0;
   ele_SIP = 0;
   ele_expectedMissingInnerHits = 0;
   ele_passConversionVeto = 0;
   passEleIdLoose = 0;
   passEleIdMedium = 0;
   passEleIdTight = 0;
   passMVAnoIsoWP90 = 0;
   passMVAnoIsoWP80 = 0;
   passMVAIsoWP90 = 0;
   passMVAIsoWP80 = 0;
   hasMatchedToZ = 0;
   passL1EG10 = 0;
   passL1EG17 = 0;
   passL1EG23 = 0;
   passL1EG23Iso = 0;
   passL1EG20Iso = 0;
   triggerPath = 0;
   triggerDecision = 0;
   passFilterEle32 = 0;
   passFilterEle23_12_leg1 = 0;
   passFilterEle23_12_leg2 = 0;
   passFilterEle115 = 0;
   passFilterEle50 = 0;
   passFilterEle25 = 0;
   passFilterEle27 = 0;
   passFilterMu12_Ele23_legEle = 0;
   passFilterMu23_Ele12_legEle = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("pvNTracks", &pvNTracks, &b_pvNTracks);
   fChain->SetBranchAddress("good_vertices", &good_vertices, &b_good_vertices);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("nPUTrue", &nPUTrue, &b_nPUTrue);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genParticles_n", &genParticles_n, &b_genParticles_n);
   fChain->SetBranchAddress("genElectron_pt", &genElectron_pt, &b_genElectron_pt);
   fChain->SetBranchAddress("genElectron_eta", &genElectron_eta, &b_genElectron_eta);
   fChain->SetBranchAddress("genElectron_phi", &genElectron_phi, &b_genElectron_phi);
   fChain->SetBranchAddress("genElectron_energy", &genElectron_energy, &b_genElectron_energy);
   fChain->SetBranchAddress("genElectron_fromZ", &genElectron_fromZ, &b_genElectron_fromZ);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
   fChain->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
   fChain->SetBranchAddress("ele_etaSC", &ele_etaSC, &b_ele_etaSC);
   fChain->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
   fChain->SetBranchAddress("ele_tricharge", &ele_tricharge, &b_ele_tricharge);
   fChain->SetBranchAddress("ele_phiSC", &ele_phiSC, &b_ele_phiSC);
   fChain->SetBranchAddress("ele_energy", &ele_energy, &b_ele_energy);
   fChain->SetBranchAddress("ele_energySC", &ele_energySC, &b_ele_energySC);
   fChain->SetBranchAddress("ele_charge", &ele_charge, &b_ele_charge);
   fChain->SetBranchAddress("ele_dEtaIn", &ele_dEtaIn, &b_ele_dEtaIn);
   fChain->SetBranchAddress("ele_dEtaSeed", &ele_dEtaSeed, &b_ele_dEtaSeed);
   fChain->SetBranchAddress("ele_dPhiIn", &ele_dPhiIn, &b_ele_dPhiIn);
   fChain->SetBranchAddress("ele_hOverE", &ele_hOverE, &b_ele_hOverE);
   fChain->SetBranchAddress("ele_full5x5_sigmaIetaIeta", &ele_full5x5_sigmaIetaIeta, &b_ele_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("ele_isoChargedHadrons", &ele_isoChargedHadrons, &b_ele_isoChargedHadrons);
   fChain->SetBranchAddress("ele_isoNeutralHadrons", &ele_isoNeutralHadrons, &b_ele_isoNeutralHadrons);
   fChain->SetBranchAddress("ele_isoPhotons", &ele_isoPhotons, &b_ele_isoPhotons);
   fChain->SetBranchAddress("ele_relCombIsoWithEA", &ele_relCombIsoWithEA, &b_ele_relCombIsoWithEA);
   fChain->SetBranchAddress("ele_isoChargedFromPU", &ele_isoChargedFromPU, &b_ele_isoChargedFromPU);
   fChain->SetBranchAddress("ele_ooEmooP", &ele_ooEmooP, &b_ele_ooEmooP);
   fChain->SetBranchAddress("ele_dr03TkSumPt", &ele_dr03TkSumPt, &b_ele_dr03TkSumPt);
   fChain->SetBranchAddress("ele_dr03EcalRecHitSumEt", &ele_dr03EcalRecHitSumEt, &b_ele_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("ele_dr03HcalDepth1TowerSumEt", &ele_dr03HcalDepth1TowerSumEt, &b_ele_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("ele_d0", &ele_d0, &b_ele_d0);
   fChain->SetBranchAddress("ele_dz", &ele_dz, &b_ele_dz);
   fChain->SetBranchAddress("ele_SIP", &ele_SIP, &b_ele_SIP);
   fChain->SetBranchAddress("ele_expectedMissingInnerHits", &ele_expectedMissingInnerHits, &b_ele_expectedMissingInnerHits);
   fChain->SetBranchAddress("ele_passConversionVeto", &ele_passConversionVeto, &b_ele_passConversionVeto);
   fChain->SetBranchAddress("passEleIdLoose", &passEleIdLoose, &b_passEleIdLoose);
   fChain->SetBranchAddress("passEleIdMedium", &passEleIdMedium, &b_passEleIdMedium);
   fChain->SetBranchAddress("passEleIdTight", &passEleIdTight, &b_passEleIdTight);
   fChain->SetBranchAddress("passMVAnoIsoWP90", &passMVAnoIsoWP90, &b_passMVAnoIsoWP90);
   fChain->SetBranchAddress("passMVAnoIsoWP80", &passMVAnoIsoWP80, &b_passMVAnoIsoWP80);
   fChain->SetBranchAddress("passMVAIsoWP90", &passMVAIsoWP90, &b_passMVAIsoWP90);
   fChain->SetBranchAddress("passMVAIsoWP80", &passMVAIsoWP80, &b_passMVAIsoWP80);
   fChain->SetBranchAddress("hasMatchedToZ", &hasMatchedToZ, &b_hasMatchedToZ);
   fChain->SetBranchAddress("passL1EG10", &passL1EG10, &b_passL1EG10);
   fChain->SetBranchAddress("passL1EG17", &passL1EG17, &b_passL1EG17);
   fChain->SetBranchAddress("passL1EG23", &passL1EG23, &b_passL1EG23);
   fChain->SetBranchAddress("passL1EG23Iso", &passL1EG23Iso, &b_passL1EG23Iso);
   fChain->SetBranchAddress("passL1EG20Iso", &passL1EG20Iso, &b_passL1EG20Iso);
   fChain->SetBranchAddress("triggerPath", &triggerPath, &b_triggerPath);
   fChain->SetBranchAddress("triggerDecision", &triggerDecision, &b_triggerDecision);
   fChain->SetBranchAddress("passFilterEle32", &passFilterEle32, &b_passFilterEle32);
   fChain->SetBranchAddress("passFilterEle23_12_leg1", &passFilterEle23_12_leg1, &b_passFilterEle23_12_leg1);
   fChain->SetBranchAddress("passFilterEle23_12_leg2", &passFilterEle23_12_leg2, &b_passFilterEle23_12_leg2);
   fChain->SetBranchAddress("passFilterEle115", &passFilterEle115, &b_passFilterEle115);
   fChain->SetBranchAddress("passFilterEle50", &passFilterEle50, &b_passFilterEle50);
   fChain->SetBranchAddress("passFilterEle25", &passFilterEle25, &b_passFilterEle25);
   fChain->SetBranchAddress("passFilterEle27", &passFilterEle27, &b_passFilterEle27);
   fChain->SetBranchAddress("passFilterMu12_Ele23_legEle", &passFilterMu12_Ele23_legEle, &b_passFilterMu12_Ele23_legEle);
   fChain->SetBranchAddress("passFilterMu23_Ele12_legEle", &passFilterMu23_Ele12_legEle, &b_passFilterMu23_Ele12_legEle);
   Notify();
}

Bool_t TagAndProbe::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TagAndProbe::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TagAndProbe::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TagAndProbe_cxx
