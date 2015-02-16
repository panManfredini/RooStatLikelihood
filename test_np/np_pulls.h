//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 10 18:11:10 2015 by ROOT version 6.02/04
// from TTree toyTree/NP post fitted
// found on file: ../mytoys.root
//////////////////////////////////////////////////////////

#ifndef np_pulls_h
#define np_pulls_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class np_pulls {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        T_Acc;
   Double_t        T_Acc_errhi;
   Double_t        T_Acc_errlo;
   Double_t        T_Bkg_norm_factor;
   Double_t        T_Bkg_norm_factor_errhi;
   Double_t        T_Bkg_norm_factor_errlo;
   Double_t        T_eff_1;
   Double_t        T_eff_10;
   Double_t        T_eff_10_errhi;
   Double_t        T_eff_10_errlo;
   Double_t        T_eff_11;
   Double_t        T_eff_11_errhi;
   Double_t        T_eff_11_errlo;
   Double_t        T_eff_12;
   Double_t        T_eff_12_errhi;
   Double_t        T_eff_12_errlo;
   Double_t        T_eff_13;
   Double_t        T_eff_13_errhi;
   Double_t        T_eff_13_errlo;
   Double_t        T_eff_14;
   Double_t        T_eff_14_errhi;
   Double_t        T_eff_14_errlo;
   Double_t        T_eff_15;
   Double_t        T_eff_15_errhi;
   Double_t        T_eff_15_errlo;
   Double_t        T_eff_16;
   Double_t        T_eff_16_errhi;
   Double_t        T_eff_16_errlo;
   Double_t        T_eff_17;
   Double_t        T_eff_17_errhi;
   Double_t        T_eff_17_errlo;
   Double_t        T_eff_18;
   Double_t        T_eff_18_errhi;
   Double_t        T_eff_18_errlo;
   Double_t        T_eff_19;
   Double_t        T_eff_19_errhi;
   Double_t        T_eff_19_errlo;
   Double_t        T_eff_1_errhi;
   Double_t        T_eff_1_errlo;
   Double_t        T_eff_2;
   Double_t        T_eff_20;
   Double_t        T_eff_20_errhi;
   Double_t        T_eff_20_errlo;
   Double_t        T_eff_2_errhi;
   Double_t        T_eff_2_errlo;
   Double_t        T_eff_3;
   Double_t        T_eff_3_errhi;
   Double_t        T_eff_3_errlo;
   Double_t        T_eff_4;
   Double_t        T_eff_4_errhi;
   Double_t        T_eff_4_errlo;
   Double_t        T_eff_5;
   Double_t        T_eff_5_errhi;
   Double_t        T_eff_5_errlo;
   Double_t        T_eff_6;
   Double_t        T_eff_6_errhi;
   Double_t        T_eff_6_errlo;
   Double_t        T_eff_7;
   Double_t        T_eff_7_errhi;
   Double_t        T_eff_7_errlo;
   Double_t        T_eff_8;
   Double_t        T_eff_8_errhi;
   Double_t        T_eff_8_errlo;
   Double_t        T_eff_9;
   Double_t        T_eff_9_errhi;
   Double_t        T_eff_9_errlo;
   Double_t        T_n_Co_1;
   Double_t        T_n_Co_10;
   Double_t        T_n_Co_10_errhi;
   Double_t        T_n_Co_10_errlo;
   Double_t        T_n_Co_11;
   Double_t        T_n_Co_11_errhi;
   Double_t        T_n_Co_11_errlo;
   Double_t        T_n_Co_12;
   Double_t        T_n_Co_12_errhi;
   Double_t        T_n_Co_12_errlo;
   Double_t        T_n_Co_13;
   Double_t        T_n_Co_13_errhi;
   Double_t        T_n_Co_13_errlo;
   Double_t        T_n_Co_14;
   Double_t        T_n_Co_14_errhi;
   Double_t        T_n_Co_14_errlo;
   Double_t        T_n_Co_15;
   Double_t        T_n_Co_15_errhi;
   Double_t        T_n_Co_15_errlo;
   Double_t        T_n_Co_16;
   Double_t        T_n_Co_16_errhi;
   Double_t        T_n_Co_16_errlo;
   Double_t        T_n_Co_17;
   Double_t        T_n_Co_17_errhi;
   Double_t        T_n_Co_17_errlo;
   Double_t        T_n_Co_18;
   Double_t        T_n_Co_18_errhi;
   Double_t        T_n_Co_18_errlo;
   Double_t        T_n_Co_19;
   Double_t        T_n_Co_19_errhi;
   Double_t        T_n_Co_19_errlo;
   Double_t        T_n_Co_1_errhi;
   Double_t        T_n_Co_1_errlo;
   Double_t        T_n_Co_2;
   Double_t        T_n_Co_20;
   Double_t        T_n_Co_20_errhi;
   Double_t        T_n_Co_20_errlo;
   Double_t        T_n_Co_2_errhi;
   Double_t        T_n_Co_2_errlo;
   Double_t        T_n_Co_3;
   Double_t        T_n_Co_3_errhi;
   Double_t        T_n_Co_3_errlo;
   Double_t        T_n_Co_4;
   Double_t        T_n_Co_4_errhi;
   Double_t        T_n_Co_4_errlo;
   Double_t        T_n_Co_5;
   Double_t        T_n_Co_5_errhi;
   Double_t        T_n_Co_5_errlo;
   Double_t        T_n_Co_6;
   Double_t        T_n_Co_6_errhi;
   Double_t        T_n_Co_6_errlo;
   Double_t        T_n_Co_7;
   Double_t        T_n_Co_7_errhi;
   Double_t        T_n_Co_7_errlo;
   Double_t        T_n_Co_8;
   Double_t        T_n_Co_8_errhi;
   Double_t        T_n_Co_8_errlo;
   Double_t        T_n_Co_9;
   Double_t        T_n_Co_9_errhi;
   Double_t        T_n_Co_9_errlo;
   Double_t        mu;
   Double_t        mu_errhi;
   Double_t        mu_errlo;

   // List of branches
   TBranch        *b_T_Acc;   //!
   TBranch        *b_T_Acc_errhi;   //!
   TBranch        *b_T_Acc_errlo;   //!
   TBranch        *b_T_Bkg_norm_factor;   //!
   TBranch        *b_T_Bkg_norm_factor_errhi;   //!
   TBranch        *b_T_Bkg_norm_factor_errlo;   //!
   TBranch        *b_T_eff_1;   //!
   TBranch        *b_T_eff_10;   //!
   TBranch        *b_T_eff_10_errhi;   //!
   TBranch        *b_T_eff_10_errlo;   //!
   TBranch        *b_T_eff_11;   //!
   TBranch        *b_T_eff_11_errhi;   //!
   TBranch        *b_T_eff_11_errlo;   //!
   TBranch        *b_T_eff_12;   //!
   TBranch        *b_T_eff_12_errhi;   //!
   TBranch        *b_T_eff_12_errlo;   //!
   TBranch        *b_T_eff_13;   //!
   TBranch        *b_T_eff_13_errhi;   //!
   TBranch        *b_T_eff_13_errlo;   //!
   TBranch        *b_T_eff_14;   //!
   TBranch        *b_T_eff_14_errhi;   //!
   TBranch        *b_T_eff_14_errlo;   //!
   TBranch        *b_T_eff_15;   //!
   TBranch        *b_T_eff_15_errhi;   //!
   TBranch        *b_T_eff_15_errlo;   //!
   TBranch        *b_T_eff_16;   //!
   TBranch        *b_T_eff_16_errhi;   //!
   TBranch        *b_T_eff_16_errlo;   //!
   TBranch        *b_T_eff_17;   //!
   TBranch        *b_T_eff_17_errhi;   //!
   TBranch        *b_T_eff_17_errlo;   //!
   TBranch        *b_T_eff_18;   //!
   TBranch        *b_T_eff_18_errhi;   //!
   TBranch        *b_T_eff_18_errlo;   //!
   TBranch        *b_T_eff_19;   //!
   TBranch        *b_T_eff_19_errhi;   //!
   TBranch        *b_T_eff_19_errlo;   //!
   TBranch        *b_T_eff_1_errhi;   //!
   TBranch        *b_T_eff_1_errlo;   //!
   TBranch        *b_T_eff_2;   //!
   TBranch        *b_T_eff_20;   //!
   TBranch        *b_T_eff_20_errhi;   //!
   TBranch        *b_T_eff_20_errlo;   //!
   TBranch        *b_T_eff_2_errhi;   //!
   TBranch        *b_T_eff_2_errlo;   //!
   TBranch        *b_T_eff_3;   //!
   TBranch        *b_T_eff_3_errhi;   //!
   TBranch        *b_T_eff_3_errlo;   //!
   TBranch        *b_T_eff_4;   //!
   TBranch        *b_T_eff_4_errhi;   //!
   TBranch        *b_T_eff_4_errlo;   //!
   TBranch        *b_T_eff_5;   //!
   TBranch        *b_T_eff_5_errhi;   //!
   TBranch        *b_T_eff_5_errlo;   //!
   TBranch        *b_T_eff_6;   //!
   TBranch        *b_T_eff_6_errhi;   //!
   TBranch        *b_T_eff_6_errlo;   //!
   TBranch        *b_T_eff_7;   //!
   TBranch        *b_T_eff_7_errhi;   //!
   TBranch        *b_T_eff_7_errlo;   //!
   TBranch        *b_T_eff_8;   //!
   TBranch        *b_T_eff_8_errhi;   //!
   TBranch        *b_T_eff_8_errlo;   //!
   TBranch        *b_T_eff_9;   //!
   TBranch        *b_T_eff_9_errhi;   //!
   TBranch        *b_T_eff_9_errlo;   //!
   TBranch        *b_T_n_Co_1;   //!
   TBranch        *b_T_n_Co_10;   //!
   TBranch        *b_T_n_Co_10_errhi;   //!
   TBranch        *b_T_n_Co_10_errlo;   //!
   TBranch        *b_T_n_Co_11;   //!
   TBranch        *b_T_n_Co_11_errhi;   //!
   TBranch        *b_T_n_Co_11_errlo;   //!
   TBranch        *b_T_n_Co_12;   //!
   TBranch        *b_T_n_Co_12_errhi;   //!
   TBranch        *b_T_n_Co_12_errlo;   //!
   TBranch        *b_T_n_Co_13;   //!
   TBranch        *b_T_n_Co_13_errhi;   //!
   TBranch        *b_T_n_Co_13_errlo;   //!
   TBranch        *b_T_n_Co_14;   //!
   TBranch        *b_T_n_Co_14_errhi;   //!
   TBranch        *b_T_n_Co_14_errlo;   //!
   TBranch        *b_T_n_Co_15;   //!
   TBranch        *b_T_n_Co_15_errhi;   //!
   TBranch        *b_T_n_Co_15_errlo;   //!
   TBranch        *b_T_n_Co_16;   //!
   TBranch        *b_T_n_Co_16_errhi;   //!
   TBranch        *b_T_n_Co_16_errlo;   //!
   TBranch        *b_T_n_Co_17;   //!
   TBranch        *b_T_n_Co_17_errhi;   //!
   TBranch        *b_T_n_Co_17_errlo;   //!
   TBranch        *b_T_n_Co_18;   //!
   TBranch        *b_T_n_Co_18_errhi;   //!
   TBranch        *b_T_n_Co_18_errlo;   //!
   TBranch        *b_T_n_Co_19;   //!
   TBranch        *b_T_n_Co_19_errhi;   //!
   TBranch        *b_T_n_Co_19_errlo;   //!
   TBranch        *b_T_n_Co_1_errhi;   //!
   TBranch        *b_T_n_Co_1_errlo;   //!
   TBranch        *b_T_n_Co_2;   //!
   TBranch        *b_T_n_Co_20;   //!
   TBranch        *b_T_n_Co_20_errhi;   //!
   TBranch        *b_T_n_Co_20_errlo;   //!
   TBranch        *b_T_n_Co_2_errhi;   //!
   TBranch        *b_T_n_Co_2_errlo;   //!
   TBranch        *b_T_n_Co_3;   //!
   TBranch        *b_T_n_Co_3_errhi;   //!
   TBranch        *b_T_n_Co_3_errlo;   //!
   TBranch        *b_T_n_Co_4;   //!
   TBranch        *b_T_n_Co_4_errhi;   //!
   TBranch        *b_T_n_Co_4_errlo;   //!
   TBranch        *b_T_n_Co_5;   //!
   TBranch        *b_T_n_Co_5_errhi;   //!
   TBranch        *b_T_n_Co_5_errlo;   //!
   TBranch        *b_T_n_Co_6;   //!
   TBranch        *b_T_n_Co_6_errhi;   //!
   TBranch        *b_T_n_Co_6_errlo;   //!
   TBranch        *b_T_n_Co_7;   //!
   TBranch        *b_T_n_Co_7_errhi;   //!
   TBranch        *b_T_n_Co_7_errlo;   //!
   TBranch        *b_T_n_Co_8;   //!
   TBranch        *b_T_n_Co_8_errhi;   //!
   TBranch        *b_T_n_Co_8_errlo;   //!
   TBranch        *b_T_n_Co_9;   //!
   TBranch        *b_T_n_Co_9_errhi;   //!
   TBranch        *b_T_n_Co_9_errlo;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_mu_errhi;   //!
   TBranch        *b_mu_errlo;   //!

   TString nameFile;
   np_pulls(TTree *tree=0, TString name="unconditional_mu1.root");
   virtual ~np_pulls();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef np_pulls_cxx
np_pulls::np_pulls(TTree *tree, TString name) : fChain(0), nameFile(name) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../mytoys.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../mytoys.root");
      }
      f->GetObject("toyTree",tree);

   }
   Init(tree);
}

np_pulls::~np_pulls()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t np_pulls::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t np_pulls::LoadTree(Long64_t entry)
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

void np_pulls::Init(TTree *tree)
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

   fChain->SetBranchAddress("T_Acc", &T_Acc, &b_T_Acc);
   fChain->SetBranchAddress("T_Acc_errhi", &T_Acc_errhi, &b_T_Acc_errhi);
   fChain->SetBranchAddress("T_Acc_errlo", &T_Acc_errlo, &b_T_Acc_errlo);
   fChain->SetBranchAddress("T_Bkg_norm_factor", &T_Bkg_norm_factor, &b_T_Bkg_norm_factor);
   fChain->SetBranchAddress("T_Bkg_norm_factor_errhi", &T_Bkg_norm_factor_errhi, &b_T_Bkg_norm_factor_errhi);
   fChain->SetBranchAddress("T_Bkg_norm_factor_errlo", &T_Bkg_norm_factor_errlo, &b_T_Bkg_norm_factor_errlo);
   fChain->SetBranchAddress("T_eff_1", &T_eff_1, &b_T_eff_1);
   fChain->SetBranchAddress("T_eff_10", &T_eff_10, &b_T_eff_10);
   fChain->SetBranchAddress("T_eff_10_errhi", &T_eff_10_errhi, &b_T_eff_10_errhi);
   fChain->SetBranchAddress("T_eff_10_errlo", &T_eff_10_errlo, &b_T_eff_10_errlo);
   fChain->SetBranchAddress("T_eff_11", &T_eff_11, &b_T_eff_11);
   fChain->SetBranchAddress("T_eff_11_errhi", &T_eff_11_errhi, &b_T_eff_11_errhi);
   fChain->SetBranchAddress("T_eff_11_errlo", &T_eff_11_errlo, &b_T_eff_11_errlo);
   fChain->SetBranchAddress("T_eff_12", &T_eff_12, &b_T_eff_12);
   fChain->SetBranchAddress("T_eff_12_errhi", &T_eff_12_errhi, &b_T_eff_12_errhi);
   fChain->SetBranchAddress("T_eff_12_errlo", &T_eff_12_errlo, &b_T_eff_12_errlo);
   fChain->SetBranchAddress("T_eff_13", &T_eff_13, &b_T_eff_13);
   fChain->SetBranchAddress("T_eff_13_errhi", &T_eff_13_errhi, &b_T_eff_13_errhi);
   fChain->SetBranchAddress("T_eff_13_errlo", &T_eff_13_errlo, &b_T_eff_13_errlo);
   fChain->SetBranchAddress("T_eff_14", &T_eff_14, &b_T_eff_14);
   fChain->SetBranchAddress("T_eff_14_errhi", &T_eff_14_errhi, &b_T_eff_14_errhi);
   fChain->SetBranchAddress("T_eff_14_errlo", &T_eff_14_errlo, &b_T_eff_14_errlo);
   fChain->SetBranchAddress("T_eff_15", &T_eff_15, &b_T_eff_15);
   fChain->SetBranchAddress("T_eff_15_errhi", &T_eff_15_errhi, &b_T_eff_15_errhi);
   fChain->SetBranchAddress("T_eff_15_errlo", &T_eff_15_errlo, &b_T_eff_15_errlo);
   fChain->SetBranchAddress("T_eff_16", &T_eff_16, &b_T_eff_16);
   fChain->SetBranchAddress("T_eff_16_errhi", &T_eff_16_errhi, &b_T_eff_16_errhi);
   fChain->SetBranchAddress("T_eff_16_errlo", &T_eff_16_errlo, &b_T_eff_16_errlo);
   fChain->SetBranchAddress("T_eff_17", &T_eff_17, &b_T_eff_17);
   fChain->SetBranchAddress("T_eff_17_errhi", &T_eff_17_errhi, &b_T_eff_17_errhi);
   fChain->SetBranchAddress("T_eff_17_errlo", &T_eff_17_errlo, &b_T_eff_17_errlo);
   fChain->SetBranchAddress("T_eff_18", &T_eff_18, &b_T_eff_18);
   fChain->SetBranchAddress("T_eff_18_errhi", &T_eff_18_errhi, &b_T_eff_18_errhi);
   fChain->SetBranchAddress("T_eff_18_errlo", &T_eff_18_errlo, &b_T_eff_18_errlo);
   fChain->SetBranchAddress("T_eff_19", &T_eff_19, &b_T_eff_19);
   fChain->SetBranchAddress("T_eff_19_errhi", &T_eff_19_errhi, &b_T_eff_19_errhi);
   fChain->SetBranchAddress("T_eff_19_errlo", &T_eff_19_errlo, &b_T_eff_19_errlo);
   fChain->SetBranchAddress("T_eff_1_errhi", &T_eff_1_errhi, &b_T_eff_1_errhi);
   fChain->SetBranchAddress("T_eff_1_errlo", &T_eff_1_errlo, &b_T_eff_1_errlo);
   fChain->SetBranchAddress("T_eff_2", &T_eff_2, &b_T_eff_2);
   fChain->SetBranchAddress("T_eff_20", &T_eff_20, &b_T_eff_20);
   fChain->SetBranchAddress("T_eff_20_errhi", &T_eff_20_errhi, &b_T_eff_20_errhi);
   fChain->SetBranchAddress("T_eff_20_errlo", &T_eff_20_errlo, &b_T_eff_20_errlo);
   fChain->SetBranchAddress("T_eff_2_errhi", &T_eff_2_errhi, &b_T_eff_2_errhi);
   fChain->SetBranchAddress("T_eff_2_errlo", &T_eff_2_errlo, &b_T_eff_2_errlo);
   fChain->SetBranchAddress("T_eff_3", &T_eff_3, &b_T_eff_3);
   fChain->SetBranchAddress("T_eff_3_errhi", &T_eff_3_errhi, &b_T_eff_3_errhi);
   fChain->SetBranchAddress("T_eff_3_errlo", &T_eff_3_errlo, &b_T_eff_3_errlo);
   fChain->SetBranchAddress("T_eff_4", &T_eff_4, &b_T_eff_4);
   fChain->SetBranchAddress("T_eff_4_errhi", &T_eff_4_errhi, &b_T_eff_4_errhi);
   fChain->SetBranchAddress("T_eff_4_errlo", &T_eff_4_errlo, &b_T_eff_4_errlo);
   fChain->SetBranchAddress("T_eff_5", &T_eff_5, &b_T_eff_5);
   fChain->SetBranchAddress("T_eff_5_errhi", &T_eff_5_errhi, &b_T_eff_5_errhi);
   fChain->SetBranchAddress("T_eff_5_errlo", &T_eff_5_errlo, &b_T_eff_5_errlo);
   fChain->SetBranchAddress("T_eff_6", &T_eff_6, &b_T_eff_6);
   fChain->SetBranchAddress("T_eff_6_errhi", &T_eff_6_errhi, &b_T_eff_6_errhi);
   fChain->SetBranchAddress("T_eff_6_errlo", &T_eff_6_errlo, &b_T_eff_6_errlo);
   fChain->SetBranchAddress("T_eff_7", &T_eff_7, &b_T_eff_7);
   fChain->SetBranchAddress("T_eff_7_errhi", &T_eff_7_errhi, &b_T_eff_7_errhi);
   fChain->SetBranchAddress("T_eff_7_errlo", &T_eff_7_errlo, &b_T_eff_7_errlo);
   fChain->SetBranchAddress("T_eff_8", &T_eff_8, &b_T_eff_8);
   fChain->SetBranchAddress("T_eff_8_errhi", &T_eff_8_errhi, &b_T_eff_8_errhi);
   fChain->SetBranchAddress("T_eff_8_errlo", &T_eff_8_errlo, &b_T_eff_8_errlo);
   fChain->SetBranchAddress("T_eff_9", &T_eff_9, &b_T_eff_9);
   fChain->SetBranchAddress("T_eff_9_errhi", &T_eff_9_errhi, &b_T_eff_9_errhi);
   fChain->SetBranchAddress("T_eff_9_errlo", &T_eff_9_errlo, &b_T_eff_9_errlo);
   fChain->SetBranchAddress("T_n_Co_1", &T_n_Co_1, &b_T_n_Co_1);
   fChain->SetBranchAddress("T_n_Co_10", &T_n_Co_10, &b_T_n_Co_10);
   fChain->SetBranchAddress("T_n_Co_10_errhi", &T_n_Co_10_errhi, &b_T_n_Co_10_errhi);
   fChain->SetBranchAddress("T_n_Co_10_errlo", &T_n_Co_10_errlo, &b_T_n_Co_10_errlo);
   fChain->SetBranchAddress("T_n_Co_11", &T_n_Co_11, &b_T_n_Co_11);
   fChain->SetBranchAddress("T_n_Co_11_errhi", &T_n_Co_11_errhi, &b_T_n_Co_11_errhi);
   fChain->SetBranchAddress("T_n_Co_11_errlo", &T_n_Co_11_errlo, &b_T_n_Co_11_errlo);
   fChain->SetBranchAddress("T_n_Co_12", &T_n_Co_12, &b_T_n_Co_12);
   fChain->SetBranchAddress("T_n_Co_12_errhi", &T_n_Co_12_errhi, &b_T_n_Co_12_errhi);
   fChain->SetBranchAddress("T_n_Co_12_errlo", &T_n_Co_12_errlo, &b_T_n_Co_12_errlo);
   fChain->SetBranchAddress("T_n_Co_13", &T_n_Co_13, &b_T_n_Co_13);
   fChain->SetBranchAddress("T_n_Co_13_errhi", &T_n_Co_13_errhi, &b_T_n_Co_13_errhi);
   fChain->SetBranchAddress("T_n_Co_13_errlo", &T_n_Co_13_errlo, &b_T_n_Co_13_errlo);
   fChain->SetBranchAddress("T_n_Co_14", &T_n_Co_14, &b_T_n_Co_14);
   fChain->SetBranchAddress("T_n_Co_14_errhi", &T_n_Co_14_errhi, &b_T_n_Co_14_errhi);
   fChain->SetBranchAddress("T_n_Co_14_errlo", &T_n_Co_14_errlo, &b_T_n_Co_14_errlo);
   fChain->SetBranchAddress("T_n_Co_15", &T_n_Co_15, &b_T_n_Co_15);
   fChain->SetBranchAddress("T_n_Co_15_errhi", &T_n_Co_15_errhi, &b_T_n_Co_15_errhi);
   fChain->SetBranchAddress("T_n_Co_15_errlo", &T_n_Co_15_errlo, &b_T_n_Co_15_errlo);
   fChain->SetBranchAddress("T_n_Co_16", &T_n_Co_16, &b_T_n_Co_16);
   fChain->SetBranchAddress("T_n_Co_16_errhi", &T_n_Co_16_errhi, &b_T_n_Co_16_errhi);
   fChain->SetBranchAddress("T_n_Co_16_errlo", &T_n_Co_16_errlo, &b_T_n_Co_16_errlo);
   fChain->SetBranchAddress("T_n_Co_17", &T_n_Co_17, &b_T_n_Co_17);
   fChain->SetBranchAddress("T_n_Co_17_errhi", &T_n_Co_17_errhi, &b_T_n_Co_17_errhi);
   fChain->SetBranchAddress("T_n_Co_17_errlo", &T_n_Co_17_errlo, &b_T_n_Co_17_errlo);
   fChain->SetBranchAddress("T_n_Co_18", &T_n_Co_18, &b_T_n_Co_18);
   fChain->SetBranchAddress("T_n_Co_18_errhi", &T_n_Co_18_errhi, &b_T_n_Co_18_errhi);
   fChain->SetBranchAddress("T_n_Co_18_errlo", &T_n_Co_18_errlo, &b_T_n_Co_18_errlo);
   fChain->SetBranchAddress("T_n_Co_19", &T_n_Co_19, &b_T_n_Co_19);
   fChain->SetBranchAddress("T_n_Co_19_errhi", &T_n_Co_19_errhi, &b_T_n_Co_19_errhi);
   fChain->SetBranchAddress("T_n_Co_19_errlo", &T_n_Co_19_errlo, &b_T_n_Co_19_errlo);
   fChain->SetBranchAddress("T_n_Co_1_errhi", &T_n_Co_1_errhi, &b_T_n_Co_1_errhi);
   fChain->SetBranchAddress("T_n_Co_1_errlo", &T_n_Co_1_errlo, &b_T_n_Co_1_errlo);
   fChain->SetBranchAddress("T_n_Co_2", &T_n_Co_2, &b_T_n_Co_2);
   fChain->SetBranchAddress("T_n_Co_20", &T_n_Co_20, &b_T_n_Co_20);
   fChain->SetBranchAddress("T_n_Co_20_errhi", &T_n_Co_20_errhi, &b_T_n_Co_20_errhi);
   fChain->SetBranchAddress("T_n_Co_20_errlo", &T_n_Co_20_errlo, &b_T_n_Co_20_errlo);
   fChain->SetBranchAddress("T_n_Co_2_errhi", &T_n_Co_2_errhi, &b_T_n_Co_2_errhi);
   fChain->SetBranchAddress("T_n_Co_2_errlo", &T_n_Co_2_errlo, &b_T_n_Co_2_errlo);
   fChain->SetBranchAddress("T_n_Co_3", &T_n_Co_3, &b_T_n_Co_3);
   fChain->SetBranchAddress("T_n_Co_3_errhi", &T_n_Co_3_errhi, &b_T_n_Co_3_errhi);
   fChain->SetBranchAddress("T_n_Co_3_errlo", &T_n_Co_3_errlo, &b_T_n_Co_3_errlo);
   fChain->SetBranchAddress("T_n_Co_4", &T_n_Co_4, &b_T_n_Co_4);
   fChain->SetBranchAddress("T_n_Co_4_errhi", &T_n_Co_4_errhi, &b_T_n_Co_4_errhi);
   fChain->SetBranchAddress("T_n_Co_4_errlo", &T_n_Co_4_errlo, &b_T_n_Co_4_errlo);
   fChain->SetBranchAddress("T_n_Co_5", &T_n_Co_5, &b_T_n_Co_5);
   fChain->SetBranchAddress("T_n_Co_5_errhi", &T_n_Co_5_errhi, &b_T_n_Co_5_errhi);
   fChain->SetBranchAddress("T_n_Co_5_errlo", &T_n_Co_5_errlo, &b_T_n_Co_5_errlo);
   fChain->SetBranchAddress("T_n_Co_6", &T_n_Co_6, &b_T_n_Co_6);
   fChain->SetBranchAddress("T_n_Co_6_errhi", &T_n_Co_6_errhi, &b_T_n_Co_6_errhi);
   fChain->SetBranchAddress("T_n_Co_6_errlo", &T_n_Co_6_errlo, &b_T_n_Co_6_errlo);
   fChain->SetBranchAddress("T_n_Co_7", &T_n_Co_7, &b_T_n_Co_7);
   fChain->SetBranchAddress("T_n_Co_7_errhi", &T_n_Co_7_errhi, &b_T_n_Co_7_errhi);
   fChain->SetBranchAddress("T_n_Co_7_errlo", &T_n_Co_7_errlo, &b_T_n_Co_7_errlo);
   fChain->SetBranchAddress("T_n_Co_8", &T_n_Co_8, &b_T_n_Co_8);
   fChain->SetBranchAddress("T_n_Co_8_errhi", &T_n_Co_8_errhi, &b_T_n_Co_8_errhi);
   fChain->SetBranchAddress("T_n_Co_8_errlo", &T_n_Co_8_errlo, &b_T_n_Co_8_errlo);
   fChain->SetBranchAddress("T_n_Co_9", &T_n_Co_9, &b_T_n_Co_9);
   fChain->SetBranchAddress("T_n_Co_9_errhi", &T_n_Co_9_errhi, &b_T_n_Co_9_errhi);
   fChain->SetBranchAddress("T_n_Co_9_errlo", &T_n_Co_9_errlo, &b_T_n_Co_9_errlo);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("mu_errhi", &mu_errhi, &b_mu_errhi);
   fChain->SetBranchAddress("mu_errlo", &mu_errlo, &b_mu_errlo);
   Notify();
}

Bool_t np_pulls::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void np_pulls::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t np_pulls::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef np_pulls_cxx
