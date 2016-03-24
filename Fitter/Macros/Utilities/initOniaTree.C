//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 16 22:46:58 2016 by ROOT version 5.34/32
// from TTree myTree/My TTree of dimuons
// found on file: OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262620_263511_noCUT.root
//////////////////////////////////////////////////////////

#ifndef initOniaTree_C
#define initOniaTree_C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

Int_t           fCurrent; //!current Tree number in a TChain

// Declaration of leaf types
UInt_t          eventNb;
UInt_t          runNb;
UInt_t          LS;
Float_t         zVtx;
Int_t           Centrality;
Int_t           nTrig;
Int_t           trigPrescale[99];   //[nTrig]
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
Float_t         rpAng[99];   //[nEP]
Float_t         rpSin[99];   //[nEP]
Float_t         rpCos[99];   //[nEP]
Int_t           Reco_QQ_size;
Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
TClonesArray    *Reco_QQ_4mom;
TClonesArray    *Reco_QQ_mupl_4mom;
TClonesArray    *Reco_QQ_mumi_4mom;
ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
ULong64_t       Reco_QQ_mupl_trig[99];   //[Reco_QQ_size]
ULong64_t       Reco_QQ_mumi_trig[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauTrue[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauTrue3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
TClonesArray    *Reco_QQ_vtx;
Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_isGoodMuon[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_isGoodMuon[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_highPurity[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_highPurity[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_TrkMuArb[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_TrkMuArb[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_TMOneStaTight[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_TMOneStaTight[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nPixValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nPixValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nMuValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nMuValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nTrkHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkHits[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_normChi2_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_normChi2_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_normChi2_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_normChi2_global[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nPixWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nPixWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nTrkWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_StationsMatched[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_StationsMatched[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dxy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dxyErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxyErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dz[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dz[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dzErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dzErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_pt_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_pt_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_pt_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_pt_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_ptErr_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_ptErr_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_ptErr_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_ptErr_global[99];   //[Reco_QQ_size]
Int_t           Reco_mu_size;
Int_t           Reco_mu_type[99];   //[Reco_mu_size]
Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
TClonesArray    *Reco_mu_4mom;
ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
Bool_t          Reco_mu_isGoodMuon[99];   //[Reco_mu_size]
Bool_t          Reco_mu_highPurity[99];   //[Reco_mu_size]
Bool_t          Reco_mu_TrkMuArb[99];   //[Reco_mu_size]
Bool_t          Reco_mu_TMOneStaTight[99];   //[Reco_mu_size]
Int_t           Reco_mu_nPixValHits[99];   //[Reco_mu_size]
Int_t           Reco_mu_nMuValHits[99];   //[Reco_mu_size]
Int_t           Reco_mu_nTrkHits[99];   //[Reco_mu_size]
Float_t         Reco_mu_normChi2_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_normChi2_global[99];   //[Reco_mu_size]
Int_t           Reco_mu_nPixWMea[99];   //[Reco_mu_size]
Int_t           Reco_mu_nTrkWMea[99];   //[Reco_mu_size]
Int_t           Reco_mu_StationsMatched[99];   //[Reco_mu_size]
Float_t         Reco_mu_dxy[99];   //[Reco_mu_size]
Float_t         Reco_mu_dxyErr[99];   //[Reco_mu_size]
Float_t         Reco_mu_dz[99];   //[Reco_mu_size]
Float_t         Reco_mu_dzErr[99];   //[Reco_mu_size]
Float_t         Reco_mu_pt_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_pt_global[99];   //[Reco_mu_size]
Float_t         Reco_mu_ptErr_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_ptErr_global[99];   //[Reco_mu_size]
Int_t           Gen_QQ_size;
Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_4mom;
Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_mupl_4mom;
TClonesArray    *Gen_QQ_mumi_4mom;
Int_t           Gen_mu_size;
Int_t           Gen_mu_type[99];   //[Gen_mu_size]
Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
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

string TreeName("hionia/myTree");

void initOniaTree(TChain *tree)
{
   std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   TChain   *fChain;   //!pointer to the analyzed TTree or TChain

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
   if (fChain->GetBranch("eventNb")) fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   if (fChain->GetBranch("runNb")) fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   if (fChain->GetBranch("LS")) fChain->SetBranchAddress("LS", &LS, &b_LS);
   if (fChain->GetBranch("zVtx")) fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   if (fChain->GetBranch("Centrality")) fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   if (fChain->GetBranch("nTrig")) fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   if (fChain->GetBranch("trigPrescale")) fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   if (fChain->GetBranch("HLTriggers")) fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   if (fChain->GetBranch("Npix")) fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
   if (fChain->GetBranch("NpixelTracks")) fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
   if (fChain->GetBranch("Ntracks")) fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   if (fChain->GetBranch("SumET_HF")) fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
   if (fChain->GetBranch("SumET_HFplus")) fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
   if (fChain->GetBranch("SumET_HFminus")) fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
   if (fChain->GetBranch("SumET_HFplusEta4")) fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
   if (fChain->GetBranch("SumET_HFminusEta4")) fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
   if (fChain->GetBranch("SumET_ET")) fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
   if (fChain->GetBranch("SumET_EE")) fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
   if (fChain->GetBranch("SumET_EB")) fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
   if (fChain->GetBranch("SumET_EEplus")) fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
   if (fChain->GetBranch("SumET_EEminus")) fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
   if (fChain->GetBranch("SumET_ZDC")) fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
   if (fChain->GetBranch("SumET_ZDCplus")) fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
   if (fChain->GetBranch("SumET_ZDCminus")) fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
   if (fChain->GetBranch("nEP")) fChain->SetBranchAddress("nEP", &nEP, &b_nEP);
   if (fChain->GetBranch("rpAng")) fChain->SetBranchAddress("rpAng", rpAng, &b_rpAng);
   if (fChain->GetBranch("rpSin")) fChain->SetBranchAddress("rpSin", rpSin, &b_rpSin);
   if (fChain->GetBranch("rpCos")) fChain->SetBranchAddress("rpCos", rpCos, &b_rpCos);
   if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   if (fChain->GetBranch("Reco_QQ_type")) fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   if (fChain->GetBranch("Reco_QQ_mupl_4mom")) fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
   if (fChain->GetBranch("Reco_QQ_mumi_4mom")) fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
   if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   if (fChain->GetBranch("Reco_QQ_mupl_trig")) fChain->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
   if (fChain->GetBranch("Reco_QQ_mumi_trig")) fChain->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
   if (fChain->GetBranch("Reco_QQ_isCowboy")) fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
   if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   if (fChain->GetBranch("Reco_QQ_ctauTrue")) fChain->SetBranchAddress("Reco_QQ_ctauTrue", Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
   if (fChain->GetBranch("Reco_QQ_ctauTrue3D")) fChain->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);
   if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
   if (fChain->GetBranch("Reco_QQ_mupl_isGoodMuon")) fChain->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
   if (fChain->GetBranch("Reco_QQ_mumi_isGoodMuon")) fChain->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
   if (fChain->GetBranch("Reco_QQ_mupl_highPurity")) fChain->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
   if (fChain->GetBranch("Reco_QQ_mumi_highPurity")) fChain->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
   if (fChain->GetBranch("Reco_QQ_mupl_TrkMuArb")) fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
   if (fChain->GetBranch("Reco_QQ_mumi_TrkMuArb")) fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
   if (fChain->GetBranch("Reco_QQ_mupl_TMOneStaTight")) fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
   if (fChain->GetBranch("Reco_QQ_mumi_TMOneStaTight")) fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
   if (fChain->GetBranch("Reco_QQ_mupl_nPixValHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
   if (fChain->GetBranch("Reco_QQ_mumi_nPixValHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
   if (fChain->GetBranch("Reco_QQ_mupl_nMuValHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
   if (fChain->GetBranch("Reco_QQ_mumi_nMuValHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
   if (fChain->GetBranch("Reco_QQ_mupl_nTrkHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
   if (fChain->GetBranch("Reco_QQ_mumi_nTrkHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
   if (fChain->GetBranch("Reco_QQ_mupl_normChi2_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
   if (fChain->GetBranch("Reco_QQ_mumi_normChi2_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
   if (fChain->GetBranch("Reco_QQ_mupl_normChi2_global")) fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
   if (fChain->GetBranch("Reco_QQ_mumi_normChi2_global")) fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
   if (fChain->GetBranch("Reco_QQ_mupl_nPixWMea")) fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
   if (fChain->GetBranch("Reco_QQ_mumi_nPixWMea")) fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
   if (fChain->GetBranch("Reco_QQ_mupl_nTrkWMea")) fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
   if (fChain->GetBranch("Reco_QQ_mumi_nTrkWMea")) fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
   if (fChain->GetBranch("Reco_QQ_mupl_StationsMatched")) fChain->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
   if (fChain->GetBranch("Reco_QQ_mumi_StationsMatched")) fChain->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
   if (fChain->GetBranch("Reco_QQ_mupl_dxy")) fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
   if (fChain->GetBranch("Reco_QQ_mumi_dxy")) fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
   if (fChain->GetBranch("Reco_QQ_mupl_dxyErr")) fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
   if (fChain->GetBranch("Reco_QQ_mumi_dxyErr")) fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
   if (fChain->GetBranch("Reco_QQ_mupl_dz")) fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
   if (fChain->GetBranch("Reco_QQ_mumi_dz")) fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
   if (fChain->GetBranch("Reco_QQ_mupl_dzErr")) fChain->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
   if (fChain->GetBranch("Reco_QQ_mumi_dzErr")) fChain->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
   if (fChain->GetBranch("Reco_QQ_mupl_pt_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_pt_inner", Reco_QQ_mupl_pt_inner, &b_Reco_QQ_mupl_pt_inner);
   if (fChain->GetBranch("Reco_QQ_mumi_pt_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_pt_inner", Reco_QQ_mumi_pt_inner, &b_Reco_QQ_mumi_pt_inner);
   if (fChain->GetBranch("Reco_QQ_mupl_pt_global")) fChain->SetBranchAddress("Reco_QQ_mupl_pt_global", Reco_QQ_mupl_pt_global, &b_Reco_QQ_mupl_pt_global);
   if (fChain->GetBranch("Reco_QQ_mumi_pt_global")) fChain->SetBranchAddress("Reco_QQ_mumi_pt_global", Reco_QQ_mumi_pt_global, &b_Reco_QQ_mumi_pt_global);
   if (fChain->GetBranch("Reco_QQ_mupl_ptErr_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_inner", Reco_QQ_mupl_ptErr_inner, &b_Reco_QQ_mupl_ptErr_inner);
   if (fChain->GetBranch("Reco_QQ_mumi_ptErr_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_inner", Reco_QQ_mumi_ptErr_inner, &b_Reco_QQ_mumi_ptErr_inner);
   if (fChain->GetBranch("Reco_QQ_mupl_ptErr_global")) fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
   if (fChain->GetBranch("Reco_QQ_mumi_ptErr_global")) fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
   if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   if (fChain->GetBranch("Reco_mu_isGoodMuon")) fChain->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
   if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
   if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
   if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
   if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
   if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
   if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   if (fChain->GetBranch("Gen_QQ_type")) fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   if (fChain->GetBranch("Gen_QQ_momId")) fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
   if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   if (fChain->GetBranch("Gen_mu_type")) fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
   if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
}
#endif // #ifndef initOniaTree_C
