# muonIDcutsOptimization
Macros for muon ID cuts optimization studies

//==================== muIDCutsOptim.C =====================\\

This macro is used to analyse the OniaTrees, and produce the basic mu-ID distributions for data and MC. The file CentralityMap_PbPb2015.txt contains the Ncoll values corresponding to a given centrality bin (binLowEdge binUpEdge Ncoll), which is used for the weighs of histos.

1) Ensure to acomodate the sizes of variables in muIDCutsOptim.h, to that of the tree you are reading.

2) In muIDCutsOptim.C select if the sample is "DATA" or "MC", and the name of the particle (just for the name of the final .root). Set fIncludeHighPurity to kFALSE if you want to remove the high purity cut from the default soft muon cuts.

3) Set the kinematical cuts and invariant mass regions for the signal and sidebands. The dimuon p_T is restricted so we analyse the same region in DATA and MC.

4) In a root session:

- Example of usage in DATA:

  > TFile* f = TFile::Open("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262548_263757_noCUT.root");

  > TDirectory* dir = f->GetDirectory("hionia")

  > TTree* tree = (TTree*)dir->Get("myTree");

  > .L muIDCutsOptim.C+

  > muIDCutsOptim toto(tree)

  > toto.Loop()



- Example of usage in MC:

  > TChain * chain = new TChain("myTree","");

  > chain->Add("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_ptJpsi36_noCUT.root/myTree");

  > chain->Add("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_ptJpsi69_noCUT.root/myTree");

  > chain->Add("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_ptJpsi912_noCUT.root/myTree");

  > .L muIDCutsOptim.C+

  > muIDCutsOptim toto(chain)

  > toto.Loop()

//===============================================================\\


//==================== muIDplots.C =====================\\

This macro is used to produce muon ID cuts plots (variables distributions, efficiencies...) from the Data and MC .root files produced by muIDCutsOptim.C 

1) If you changed the invariant mass regions for the signal and sidebands in muIDCutsOptim.C (step 2 above), you will need to set accordingly the varialbe bkgNorm. bkgNorm = (wide of signal region) / (wide of sidebands).

2) In a root session:

  > .L muIDplots.C+

  > muIDplots("histos_MC_JPsi_HPincl.root","histos_DATA_JPsi_HPincl.root")

Note that the signal (from MC) goes as first argument and background (from Data) as second. The naming of the histos comparing data and MC is made following this order. WARNING: Take care about this if you change the order of the input files!!

//===============================================================\\

