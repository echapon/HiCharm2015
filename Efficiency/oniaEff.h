//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 22 10:54:36 2016 by ROOT version 6.06/00
// from TTree myTree/My TTree of dimuons
// found on file: root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root
//////////////////////////////////////////////////////////

#ifndef oniaEff_h
#define oniaEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "../Fitter/Macros/Utilities/EVENTUTILS.h"

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TLorentzVector.h"

class oniaEff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventNb;
   UInt_t          runNb;
   UInt_t          LS;
   Float_t         zVtx;
   Int_t           Centrality;
   Int_t           nTrig;
   Int_t           trigPrescale[43];   //[nTrig]
   ULong64_t       HLTriggers;
   Int_t           Npix;
   Int_t           NpixelTracks;
   Int_t           Ntracks;
   Float_t         SumET_HF;
   Float_t         SumET_HFplus;
   Float_t         SumET_HFminus;
   Float_t         SumET_HFplusEta4;
   Float_t         SumET_HFminusEta4;
   Float_t         SumET_ET;
   Float_t         SumET_EE;
   Float_t         SumET_EB;
   Float_t         SumET_EEplus;
   Float_t         SumET_EEminus;
   Float_t         SumET_ZDC;
   Float_t         SumET_ZDCplus;
   Float_t         SumET_ZDCminus;
   Int_t           nEP;
   Float_t         rpAng[29];   //[nEP]
   Float_t         rpSin[29];   //[nEP]
   Float_t         rpCos[29];   //[nEP]
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[36];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   TClonesArray    *Reco_QQ_mupl_4mom;
   TClonesArray    *Reco_QQ_mumi_4mom;
   ULong64_t       Reco_QQ_trig[36];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mupl_trig[36];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mumi_trig[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauTrue[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauTrue3D[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[36];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_isGoodMuon[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_isGoodMuon[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_highPurity[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_highPurity[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TrkMuArb[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TrkMuArb[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TMOneStaTight[36];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TMOneStaTight[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixValHits[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixValHits[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nMuValHits[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nMuValHits[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkHits[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkHits[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_global[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_global[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixWMea[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixWMea[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkWMea[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkWMea[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_StationsMatched[36];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_StationsMatched[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxy[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxy[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxyErr[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxyErr[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dz[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dz[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dzErr[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dzErr[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_global[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_global[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_inner[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_global[36];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_global[36];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[9];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[9];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[9];   //[Reco_mu_size]
   Bool_t          Reco_mu_isGoodMuon[9];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[9];   //[Reco_mu_size]
   Bool_t          Reco_mu_TrkMuArb[9];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[9];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[9];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[9];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[9];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[9];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[9];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[9];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[9];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[9];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[9];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[9];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[9];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[9];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_inner[9];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_global[9];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[9];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_global[9];   //[Reco_mu_size]
   Int_t           Gen_QQ_size;
   Int_t           Gen_QQ_type[3];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_4mom;
   Int_t           Gen_QQ_momId[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_mupl_4mom;
   TClonesArray    *Gen_QQ_mumi_4mom;
   Int_t           Gen_mu_size;
   Int_t           Gen_mu_type[14];   //[Gen_mu_size]
   Int_t           Gen_mu_charge[14];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_LS;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_Centrality;   //!
   TBranch        *b_nTrig;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_Npix;   //!
   TBranch        *b_NpixelTracks;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_SumET_HF;   //!
   TBranch        *b_SumET_HFplus;   //!
   TBranch        *b_SumET_HFminus;   //!
   TBranch        *b_SumET_HFplusEta4;   //!
   TBranch        *b_SumET_HFminusEta4;   //!
   TBranch        *b_SumET_ET;   //!
   TBranch        *b_SumET_EE;   //!
   TBranch        *b_SumET_EB;   //!
   TBranch        *b_SumET_EEplus;   //!
   TBranch        *b_SumET_EEminus;   //!
   TBranch        *b_SumET_ZDC;   //!
   TBranch        *b_SumET_ZDCplus;   //!
   TBranch        *b_SumET_ZDCminus;   //!
   TBranch        *b_nEP;   //!
   TBranch        *b_rpAng;   //!
   TBranch        *b_rpSin;   //!
   TBranch        *b_rpCos;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_4mom;   //!
   TBranch        *b_Reco_QQ_mumi_4mom;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_mupl_trig;   //!
   TBranch        *b_Reco_QQ_mumi_trig;   //!
   TBranch        *b_Reco_QQ_isCowboy;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_ctauErr;   //!
   TBranch        *b_Reco_QQ_ctau3D;   //!
   TBranch        *b_Reco_QQ_ctauErr3D;   //!
   TBranch        *b_Reco_QQ_ctauTrue;   //!
   TBranch        *b_Reco_QQ_ctauTrue3D;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_vtx;   //!
   TBranch        *b_Reco_QQ_Ntrk;   //!
   TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mupl_highPurity;   //!
   TBranch        *b_Reco_QQ_mumi_highPurity;   //!
   TBranch        *b_Reco_QQ_mupl_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mumi_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mupl_dxy;   //!
   TBranch        *b_Reco_QQ_mumi_dxy;   //!
   TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
   TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
   TBranch        *b_Reco_QQ_mupl_dz;   //!
   TBranch        *b_Reco_QQ_mumi_dz;   //!
   TBranch        *b_Reco_QQ_mupl_dzErr;   //!
   TBranch        *b_Reco_QQ_mumi_dzErr;   //!
   TBranch        *b_Reco_QQ_mupl_pt_inner;   //!
   TBranch        *b_Reco_QQ_mumi_pt_inner;   //!
   TBranch        *b_Reco_QQ_mupl_pt_global;   //!
   TBranch        *b_Reco_QQ_mumi_pt_global;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_charge;   //!
   TBranch        *b_Reco_mu_4mom;   //!
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_isGoodMuon;   //!
   TBranch        *b_Reco_mu_highPurity;   //!
   TBranch        *b_Reco_mu_TrkMuArb;   //!
   TBranch        *b_Reco_mu_TMOneStaTight;   //!
   TBranch        *b_Reco_mu_nPixValHits;   //!
   TBranch        *b_Reco_mu_nMuValHits;   //!
   TBranch        *b_Reco_mu_nTrkHits;   //!
   TBranch        *b_Reco_mu_normChi2_inner;   //!
   TBranch        *b_Reco_mu_normChi2_global;   //!
   TBranch        *b_Reco_mu_nPixWMea;   //!
   TBranch        *b_Reco_mu_nTrkWMea;   //!
   TBranch        *b_Reco_mu_StationsMatched;   //!
   TBranch        *b_Reco_mu_dxy;   //!
   TBranch        *b_Reco_mu_dxyErr;   //!
   TBranch        *b_Reco_mu_dz;   //!
   TBranch        *b_Reco_mu_dzErr;   //!
   TBranch        *b_Reco_mu_pt_inner;   //!
   TBranch        *b_Reco_mu_pt_global;   //!
   TBranch        *b_Reco_mu_ptErr_inner;   //!
   TBranch        *b_Reco_mu_ptErr_global;   //!
   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_type;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_momId;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_4mom;   //!
   TBranch        *b_Gen_QQ_mumi_4mom;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_type;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!

   oniaEff(TTree *tree=0);
   virtual ~oniaEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char* fname, bool ispbpb=false, bool isPsip=false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Bool_t isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit);
   Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon);
   Bool_t areMuonsInAcceptance2015 (Int_t iRecoQQ);
   Bool_t passQualityCuts2015 (Int_t iRecoQQ);
};

#endif

#ifdef oniaEff_cxx
oniaEff::oniaEff(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root:/hionia");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

oniaEff::~oniaEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t oniaEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t oniaEff::LoadTree(Long64_t entry)
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

void oniaEff::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_QQ_4mom = 0;
   Reco_QQ_mupl_4mom = 0;
   Reco_QQ_mumi_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;
   Gen_QQ_4mom = 0;
   Gen_QQ_mupl_4mom = 0;
   Gen_QQ_mumi_4mom = 0;
   Gen_mu_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("LS", &LS, &b_LS);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
   fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
   fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
   fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
   fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
   fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
   fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
   fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
   fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
   fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
   fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
   fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
   fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
   fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
   fChain->SetBranchAddress("nEP", &nEP, &b_nEP);
   fChain->SetBranchAddress("rpAng", rpAng, &b_rpAng);
   fChain->SetBranchAddress("rpSin", rpSin, &b_rpSin);
   fChain->SetBranchAddress("rpCos", rpCos, &b_rpCos);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
   fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
   fChain->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
   fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   fChain->SetBranchAddress("Reco_QQ_ctauTrue", Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
   fChain->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
   fChain->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
   fChain->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
   fChain->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
   fChain->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
   fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
   fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
   fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
   fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
   fChain->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
   fChain->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
   fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
   fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
   fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
   fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
   fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
   fChain->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
   fChain->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
   fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
   fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
   fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
   fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
   fChain->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
   fChain->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
   fChain->SetBranchAddress("Reco_QQ_mupl_pt_inner", Reco_QQ_mupl_pt_inner, &b_Reco_QQ_mupl_pt_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_pt_inner", Reco_QQ_mumi_pt_inner, &b_Reco_QQ_mumi_pt_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_pt_global", Reco_QQ_mupl_pt_global, &b_Reco_QQ_mupl_pt_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_pt_global", Reco_QQ_mumi_pt_global, &b_Reco_QQ_mumi_pt_global);
   fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_inner", Reco_QQ_mupl_ptErr_inner, &b_Reco_QQ_mupl_ptErr_inner);
   fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_inner", Reco_QQ_mumi_ptErr_inner, &b_Reco_QQ_mumi_ptErr_inner);
   fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
   fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
   fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
   fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
   fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
   fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
   fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
   fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   Notify();
}

Bool_t oniaEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void oniaEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t oniaEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
Bool_t oniaEff::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit) 
{
   Bool_t cond = true;
   cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) ); 
   cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
   return cond;
};
Bool_t oniaEff::isGlobalMuonInAccept2015 (TLorentzVector* Muon) 
{
   return (fabs(Muon->Eta()) < 2.4 &&
         ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
          (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
          (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
};

Bool_t oniaEff::areMuonsInAcceptance2015 (Int_t iRecoQQ)
{
   TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
   TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
   return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
};  

Bool_t oniaEff::passQualityCuts2015 (Int_t iRecoQQ) 
{
   Bool_t cond = true;

   // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
   cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
   cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
   cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
   cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
   cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);

   // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
   cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
   cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
   cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
   cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
   cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);

   cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);

   return cond;
}; 
#endif // #ifdef oniaEff_cxx
