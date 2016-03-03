// a macro that simply extracts the number of jpsis from a list of files containing fit results

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include <iostream>

using namespace std;

void simplePrintResults() {
   vector<string> filenames;
   filenames.push_back("Output/Test2/result/FIT_DATA_Psi2SJpsi_PPPrompt_Bkg_SecondOrderChebychev_pt65300_rap016_cent0200_262620_263757.root");

   const char* parname = "N_Jpsi_PP";

   vector<string>::iterator it = filenames.begin();
   for (it; it<filenames.end(); it++) {
      TFile *f = new TFile(it->c_str());
      if (!f) {
         cout << "Error, " << *it << " not found" << endl;
         continue;
      }
      RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
      if (!ws) {
         cout << "Error, workspace not found in " << *it << endl;
         continue;
      }

      RooRealVar *var = ws->var(parname);
      if (!ws) {
         cout << "Error, variable " << parname << " not found in " << *it << endl;
         continue;
      }
      cout << *it << " " << var->getVal() << " +- " << var->getError() << endl;
   }
}
