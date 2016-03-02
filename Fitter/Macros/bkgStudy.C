#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooWorkspace.h"
#include "RooAbsPdf.h"

#include "TFile.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <iostream>

using namespace std;
using namespace RooFit;

int cFactor = 2; // 1 for regular LLR test, 2 for Akaike Information Criterion (AIC)

void bkgStudy(bool isPbPb = false) {
   vector<string> filenames;
   // Background functions with the smallest number of parameters should go first
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_Exponential.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_FirstOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_SecondOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_ThirdOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_FourthOrderChebychev.root");
   
   const char* fcnname = isPbPb ? "pdfMASS_Bkg_PbPb" : "pdfMASS_Bkg_PP";
   const char* fcnname_full = isPbPb ? "pdfMASS_Tot_PbPb" : "pdfMASS_Tot_PP";
   const char* dataname = isPbPb ? "dOS_DATA_PbPb" : "dOS_DATA_PP";

   vector<string> fcntags;
   fcntags.push_back("Exponential");
   fcntags.push_back("FirstOrderChebychev");
   fcntags.push_back("SecondOrderChebychev");
   fcntags.push_back("ThirdOrderChebychev");
   fcntags.push_back("FourthOrderChebychev");

   vector<double> nlls;
   vector<int> npars;


   // get the bkg functions
   for (unsigned int i=0; i<filenames.size(); i++) {
      TFile *tfile = new TFile(filenames[i].c_str());
      RooWorkspace *ws = (RooWorkspace*) tfile->Get("myws");
      RooAbsReal *nll =ws->pdf(fcnname_full)->createNLL(*ws->data(dataname));
      // ws->pdf(fcnname)->fitTo(*ws->data("dOS_DATA_PP"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow"),NormRange("MassWindow"));
      nlls.push_back(nll->getVal());
      int npar = ws->pdf(fcnname)->getParameters(*ws->data(dataname))->getSize();
      npars.push_back(npar);
   }

   // print results
   for (unsigned int ibkg=0; ibkg<filenames.size(); ibkg++) {
      cout.setf(ios::fixed, ios::floatfield);
      // cout << fcntags[ibkg] << " (" << npars[ibkg] << ") " << setprecision(2) << ((ibkg==0) ? nlls[ibkg] : 2.*(nlls[ibkg]-nlls[0]));
      cout << fcntags[ibkg] << " (" << npars[ibkg] << ") " << setprecision(2) << nlls[ibkg];
      cout << setprecision(4);
      cout.unsetf(ios::floatfield);
      for (unsigned int ibkg2=0; ibkg2<ibkg; ibkg2++) {
         double deltanll = -2.*(nlls[ibkg]-nlls[ibkg2]);
         int deltanpar = abs(npars[ibkg]-npars[ibkg2]);
         double prob = deltanll > 0 ? 100.*TMath::Prob(deltanll,cFactor*deltanpar) : 100;
         cout << "\t" <<  prob;
      } // for (ibkg2)
      cout << endl;
   } // for (ibkg)
}

