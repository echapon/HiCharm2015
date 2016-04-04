//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr  3 18:12:23 2016 by ROOT version 5.34/18
// from TTree plotParamEvolution/fit results
// found on file: MCSignalParameters.root
//////////////////////////////////////////////////////////

#ifndef plotParamEvolution_h
#define plotParamEvolution_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class plotParamEvolution {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         ptmin;
   Float_t         ptmax;
   Float_t         ymin;
   Float_t         ymax;
   Float_t         centmin;
   Float_t         centmax;
   Char_t          jpsiName[23];
   Char_t          psipName[5];
   Char_t          bkgName;
   Char_t          collSystem[5];
   Float_t         alpha_Jpsi_val;
   Float_t         alpha_Jpsi_err;
   Float_t         f_Jpsi_val;
   Float_t         f_Jpsi_err;
   Float_t         n_Jpsi_val;
   Float_t         n_Jpsi_err;
   Float_t         sigma1_Jpsi_val;
   Float_t         sigma1_Jpsi_err;
   Float_t         rSigma21_Jpsi_val;
   Float_t         rSigma21_Jpsi_err;

   // List of branches
   TBranch        *b_ptmin;   //!
   TBranch        *b_ptmax;   //!
   TBranch        *b_ymin;   //!
   TBranch        *b_ymax;   //!
   TBranch        *b_centmin;   //!
   TBranch        *b_centmax;   //!
   TBranch        *b_jpsiName;   //!
   TBranch        *b_psipName;   //!
   TBranch        *b_bkgName;   //!
   TBranch        *b_collSystem;   //!
   TBranch        *b_alpha_Jpsi_val;   //!
   TBranch        *b_alpha_Jpsi_err;   //!
   TBranch        *b_f_Jpsi_val;   //!
   TBranch        *b_f_Jpsi_err;   //!
   TBranch        *b_n_Jpsi_val;   //!
   TBranch        *b_n_Jpsi_err;   //!
   TBranch        *b_sigma1_Jpsi_val;   //!
   TBranch        *b_sigma1_Jpsi_err;   //!
   TBranch        *b_rSigma21_Jpsi_val;   //!
   TBranch        *b_rSigma21_Jpsi_err;   //!

   plotParamEvolution(TTree *tree=0);
   virtual ~plotParamEvolution();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef plotParamEvolution_cxx
plotParamEvolution::plotParamEvolution(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MCSignalParameters.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MCSignalParameters.root");
      }
      f->GetObject("fitresults",tree);

   }
   Init(tree);
}

plotParamEvolution::~plotParamEvolution()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t plotParamEvolution::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plotParamEvolution::LoadTree(Long64_t entry)
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

void plotParamEvolution::Init(TTree *tree)
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

   fChain->SetBranchAddress("ptmin", &ptmin, &b_ptmin);
   fChain->SetBranchAddress("ptmax", &ptmax, &b_ptmax);
   fChain->SetBranchAddress("ymin", &ymin, &b_ymin);
   fChain->SetBranchAddress("ymax", &ymax, &b_ymax);
   fChain->SetBranchAddress("centmin", &centmin, &b_centmin);
   fChain->SetBranchAddress("centmax", &centmax, &b_centmax);
   fChain->SetBranchAddress("jpsiName", jpsiName, &b_jpsiName);
   fChain->SetBranchAddress("psipName", psipName, &b_psipName);
   fChain->SetBranchAddress("bkgName", &bkgName, &b_bkgName);
   fChain->SetBranchAddress("collSystem", collSystem, &b_collSystem);
   fChain->SetBranchAddress("alpha_Jpsi_val", &alpha_Jpsi_val, &b_alpha_Jpsi_val);
   fChain->SetBranchAddress("alpha_Jpsi_err", &alpha_Jpsi_err, &b_alpha_Jpsi_err);
   fChain->SetBranchAddress("f_Jpsi_val", &f_Jpsi_val, &b_f_Jpsi_val);
   fChain->SetBranchAddress("f_Jpsi_err", &f_Jpsi_err, &b_f_Jpsi_err);
   fChain->SetBranchAddress("n_Jpsi_val", &n_Jpsi_val, &b_n_Jpsi_val);
   fChain->SetBranchAddress("n_Jpsi_err", &n_Jpsi_err, &b_n_Jpsi_err);
   fChain->SetBranchAddress("sigma1_Jpsi_val", &sigma1_Jpsi_val, &b_sigma1_Jpsi_val);
   fChain->SetBranchAddress("sigma1_Jpsi_err", &sigma1_Jpsi_err, &b_sigma1_Jpsi_err);
   fChain->SetBranchAddress("rSigma21_Jpsi_val", &rSigma21_Jpsi_val, &b_rSigma21_Jpsi_val);
   fChain->SetBranchAddress("rSigma21_Jpsi_err", &rSigma21_Jpsi_err, &b_rSigma21_Jpsi_err);
   Notify();
}

Bool_t plotParamEvolution::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotParamEvolution::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotParamEvolution::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plotParamEvolution_cxx
