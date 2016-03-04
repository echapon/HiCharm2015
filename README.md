# HiCharm2015
Repository for code regarding charmoniun in 5.02TeV data.

## Fitter
* See documentation in https://twiki.cern.ch/twiki/bin/view/CMS/HiCharm2015Fitter
* fitter.C: master file
* The fitter can be run with the following command:  root -l -b -q fitter.C+'("name of work dir")' 
* clean_wd.sh: source this script to clean the working directory (RooDatasets, plots)
* simplePrintResults.C: minimalistic macro to print the fit result from a list of files
* Input/createInput.py: automatically create the input (configuration) files
