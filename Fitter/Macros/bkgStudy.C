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

void bkgStudy() {
   vector<string> filenames;
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_FirstOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_SecondOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_ThirdOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_FourthOrderChebychev.root");
   filenames.push_back("Workspace/DATA/ws_DATA_Psi2SJpsi_PPPrompt_pt65300_rap1624_cent0200_262620_263757_Exponential.root");
   
   const char* fcnname = "pdfMASS_Bkg_PP";
   const char* fcnname_full = "pdfMASS_Tot_PP";

   vector<string> fcntags;
   fcntags.push_back("FirstOrderChebychev");
   fcntags.push_back("SecondOrderChebychev");
   fcntags.push_back("ThirdOrderChebychev");
   fcntags.push_back("FourthOrderChebychev");
   fcntags.push_back("Exponential");

   vector<double> nlls;
   vector<int> npars;


   // get the bkg functions
   for (unsigned int i=0; i<filenames.size(); i++) {
      TFile *tfile = new TFile(filenames[i].c_str());
      RooWorkspace *ws = (RooWorkspace*) tfile->Get("myws");
      RooAbsReal *nll =ws->pdf(fcnname_full)->createNLL(*ws->data("dOS_DATA_PP"));
      // ws->pdf(fcnname)->fitTo(*ws->data("dOS_DATA_PP"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow"),NormRange("MassWindow"));
      nlls.push_back(nll->getVal());
      int npar = ws->pdf(fcnname)->getParameters(*ws->data("dOS_DATA_PP"))->getSize();
      npars.push_back(npar);
   }

   // print results
   for (unsigned int ibkg=0; ibkg<filenames.size(); ibkg++) {
      cout << fcntags[ibkg] << " (" << npars[ibkg] << ") " << ((ibkg==0) ? nlls[ibkg] : 2.*(nlls[ibkg]-nlls[0]));
      for (unsigned int ibkg2=0; ibkg2<ibkg; ibkg2++) {
         double deltanll = 2.*fabs(nlls[ibkg]-nlls[ibkg2]);
         int deltanpar = abs(npars[ibkg]-npars[ibkg2]);
         double prob = 100.*TMath::Prob(deltanll,deltanpar);
         cout << "\t" << 100.*TMath::Prob(deltanll,deltanpar);
      } // for (ibkg2)
      cout << endl;
   } // for (ibkg)
}
