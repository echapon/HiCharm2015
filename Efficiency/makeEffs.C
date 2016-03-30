#include "oniaEff.C"
#include "TChain.h"

#include <iostream>

using namespace std;

void makeEffs() {
   // PP
   TChain *tch_jpsi_pp = new TChain("hionia/myTree");
   tch_jpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
   TChain *tch_psi2s_pp = new TChain("hionia/myTree");
   tch_psi2s_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_Psi2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
   TChain *tch_npjpsi_pp = new TChain("hionia/myTree");
   tch_npjpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");

   // PbPb
   TChain *tch_jpsi_pbpb = new TChain("hionia/myTree");
   // tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   // tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   TChain *tch_psi2s_pbpb = new TChain("hionia/myTree");
   // tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_psi2s_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_Psi2SMM_ptPsi2_15_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   TChain *tch_npjpsi_pbpb = new TChain("hionia/myTree");
   // tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   // tch_npjpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_BJpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");

   // make the efficiency histos
   cout << "Efficiencies for pp prompt Jpsi" << endl;
   oniaEff obj_jpsi_pp(tch_jpsi_pp);
   obj_jpsi_pp.Loop("files/histos_jpsi_pp.root",false);
   cout << "Efficiencies for pp prompt Psi(2S)" << endl;
   oniaEff obj_psi2s_pp(tch_psi2s_pp);
   obj_psi2s_pp.Loop("files/histos_psi2s_pp.root",false);
   cout << "Efficiencies for pp non-prompt Jpsi" << endl;
   oniaEff obj_npjpsi_pp(tch_npjpsi_pp);
   obj_npjpsi_pp.Loop("files/histos_npjpsi_pp.root",false);
   cout << "Efficiencies for pbpb prompt Jpsi" << endl;
   oniaEff obj_jpsi_pbpb(tch_jpsi_pbpb);
   obj_jpsi_pbpb.Loop("files/histos_jpsi_pbpb.root",true);
   cout << "Efficiencies for pbpb prompt Psi(2S)" << endl;
   oniaEff obj_psi2s_pbpb(tch_psi2s_pbpb);
   obj_psi2s_pbpb.Loop("files/histos_psi2s_pbpb.root",true);
   cout << "Efficiencies for pbpb non-prompt Jpsi" << endl;
   oniaEff obj_npjpsi_pbpb(tch_npjpsi_pbpb);
   obj_npjpsi_pbpb.Loop("files/histos_npjpsi_pbpb.root",true);
}
