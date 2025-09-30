//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 30 16:35:14 2025 by ROOT version 6.32.06
// from TTree T/T
// found on file: 72867/run3auau_new_newcalib_v008-00072867-000000.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

// Header file for the classes stored in the TTree if any.

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event_id;
   Float_t         zvrtx;
   Int_t           centrality;
   Float_t         rho_M;
   Float_t         rho_M_sigma;
   Float_t         rho_A;
   Float_t         rho_A_sigma;
   Float_t         mbd_total_charge_south;
   Float_t         mbd_total_charge_north;
   Float_t         mbd_mean_time_south;
   Float_t         mbd_mean_time_north;
   Int_t           mbd_total_npmt_south;
   Int_t           mbd_total_npmt_north;
   Int_t           gl1_scaled_trigger_vec[64];
   Int_t           gl1_live_trigger_vec[64];
   Float_t         emcal_energy[24][64];
   Float_t         emcal_time[24][64];
   Int_t           emcal_isgood[24][64];
   Float_t         hcalin_energy[24][64];
   Float_t         hcalin_time[24][64];
   Int_t           hcalin_isgood[24][64];
   Float_t         hcalout_energy[24][64];
   Float_t         hcalout_time[24][64];
   Float_t         sepd_energy[744];
   Int_t           sepd_isgood[744];
   Int_t           sepd_arm[744];
   Float_t         sepd_radius[744];
   Float_t         sepd_phi[744];

   // List of branches
   TBranch        *b_event_id;   //!
   TBranch        *b_zvrtx;   //!
   TBranch        *b_centrality;   //!
   TBranch        *b_rho_M;   //!
   TBranch        *b_rho_M_sigma;   //!
   TBranch        *b_rho_A;   //!
   TBranch        *b_rho_A_sigma;   //!
   TBranch        *b_mbd_total_charge_south;   //!
   TBranch        *b_mbd_total_charge_north;   //!
   TBranch        *b_mbd_mean_time_south;   //!
   TBranch        *b_mbd_mean_time_north;   //!
   TBranch        *b_mbd_total_npmt_south;   //!
   TBranch        *b_mbd_total_npmt_north;   //!
   TBranch        *b_gl1_scaled_trigger_vec;   //!
   TBranch        *b_gl1_live_trigger_vec;   //!
   TBranch        *b_emcal_energy;   //!
   TBranch        *b_emcal_time;   //!
   TBranch        *b_emcal_isgood;   //!
   TBranch        *b_hcalin_energy;   //!
   TBranch        *b_hcalin_time;   //!
   TBranch        *b_hcalin_isgood;   //!
   TBranch        *b_hcalout_energy;   //!
   TBranch        *b_hcalout_time;   //!
   TBranch        *b_sepd_energy;   //!
   TBranch        *b_sepd_isgood;   //!
   TBranch        *b_sepd_arm;   //!
   TBranch        *b_sepd_radius;   //!
   TBranch        *b_sepd_phi;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("72867/run3auau_new_newcalib_v008-00072867-000000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("72867/run3auau_new_newcalib_v008-00072867-000000.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
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

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("zvrtx", &zvrtx, &b_zvrtx);
   fChain->SetBranchAddress("centrality", &centrality, &b_centrality);
   fChain->SetBranchAddress("rho_M", &rho_M, &b_rho_M);
   fChain->SetBranchAddress("rho_M_sigma", &rho_M_sigma, &b_rho_M_sigma);
   fChain->SetBranchAddress("rho_A", &rho_A, &b_rho_A);
   fChain->SetBranchAddress("rho_A_sigma", &rho_A_sigma, &b_rho_A_sigma);
   fChain->SetBranchAddress("mbd_total_charge_south", &mbd_total_charge_south, &b_mbd_total_charge_south);
   fChain->SetBranchAddress("mbd_total_charge_north", &mbd_total_charge_north, &b_mbd_total_charge_north);
   fChain->SetBranchAddress("mbd_mean_time_south", &mbd_mean_time_south, &b_mbd_mean_time_south);
   fChain->SetBranchAddress("mbd_mean_time_north", &mbd_mean_time_north, &b_mbd_mean_time_north);
   fChain->SetBranchAddress("mbd_total_npmt_south", &mbd_total_npmt_south, &b_mbd_total_npmt_south);
   fChain->SetBranchAddress("mbd_total_npmt_north", &mbd_total_npmt_north, &b_mbd_total_npmt_north);
   fChain->SetBranchAddress("gl1_scaled_trigger_vec", gl1_scaled_trigger_vec, &b_gl1_scaled_trigger_vec);
   fChain->SetBranchAddress("gl1_live_trigger_vec", gl1_live_trigger_vec, &b_gl1_live_trigger_vec);
   fChain->SetBranchAddress("emcal_energy", emcal_energy, &b_emcal_energy);
   fChain->SetBranchAddress("emcal_time", emcal_time, &b_emcal_time);
   fChain->SetBranchAddress("emcal_isgood", emcal_isgood, &b_emcal_isgood);
   fChain->SetBranchAddress("hcalin_energy", hcalin_energy, &b_hcalin_energy);
   fChain->SetBranchAddress("hcalin_time", hcalin_time, &b_hcalin_time);
   fChain->SetBranchAddress("hcalin_isgood", hcalin_isgood, &b_hcalin_isgood);
   fChain->SetBranchAddress("hcalout_energy", hcalout_energy, &b_hcalout_energy);
   fChain->SetBranchAddress("hcalout_time", hcalout_time, &b_hcalout_time);
   fChain->SetBranchAddress("sepd_energy", sepd_energy, &b_sepd_energy);
   fChain->SetBranchAddress("sepd_isgood", sepd_isgood, &b_sepd_isgood);
   fChain->SetBranchAddress("sepd_arm", sepd_arm, &b_sepd_arm);
   fChain->SetBranchAddress("sepd_radius", sepd_radius, &b_sepd_radius);
   fChain->SetBranchAddress("sepd_phi", sepd_phi, &b_sepd_phi);
   Notify();
}

bool test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
