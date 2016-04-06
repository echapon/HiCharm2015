# HiCharm2015
Repository for code regarding charmoniun in 5.02TeV data.

DO NOT PUSH ANYTHING HERE ANYMORE. THE CODE HAS BEEN MOVED TO https://github.com/CMS-HIN-dilepton/DimuonCADIs

## Fitter
* See documentation in https://twiki.cern.ch/twiki/bin/view/CMS/HiCharm2015Fitter
* fitter.C: master file
 * The fitter can be run with the following command:  root -l -b -q fitter.C+'("name of work dir")' 
* clean_wd.sh: source this script to clean the working directory (RooDatasets, plots)
* simplePrintResults.C: minimalistic macro to print the fit result from a list of files
* Input/createInput.py: automatically create the input (configuration) files
* plotResults.C: plotting macro containing two main functions. It is largely automated.
 * plotPt("name of work dir"): make plot vs pt
 * plotCent("name of work dir"): make plot cs centrality
* results2tree.C: parse an output directory and put the results into a TTree
 * results2tree("name of work dir","var1,var2,var3")
 * NB: the names of the variables ("var1", etc. in the example above) should be the actual name of the variable in the workspace, without "\_PP" or "\_PbPb". e.g. "sigma1\_Jpsi"
* plotParamEvolution.C: macro to plot the evolution of the signal parameters as a function of pt and centrality. it takes the TTree from the results2tree.C macro as input.
* plotVars.C: alternative macro to plot the evolution of variables.
 * plotPt(const char\* workDirName, const char\* varname, const char* collTag="", bool plotErr=true, bool isMC=false)
   where workDirName is the tag that was given to fitter.C, varname is the name of the variable, 
   collTag can be PP or PbPb (leave empty for both), plotErr tells if errors should be plotted, and isMC is false for
   data, true for MC.
 * plotRap, plotCent work the same
 * plotFiles(const char* workDirNames, const char* varname, const char* xaxis, float rapmin, float rapmax, 
   float ptmin, float ptmax, int centmin, int centmax, const char* collTag="PP", bool plotErr=true, bool isMC=false)
   compare different working directories (workDirNames should be a comma-separated list "dir1,dir2,..,dirN")
   for the bins defined by rapmin, rapmax, ptmin, ptmax, centmin, centmax.


## Efficiency
* makeEffs.C: make the histograms (numerators and denominators) from the onia trees
* plotEffs.C: make the efficiency plots from the histograms we just created
