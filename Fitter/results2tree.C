#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <vector>
#include <cstring>

#include "Macros/Utilities/resultUtils.h"

using namespace std;

struct poi {
   Char_t name[64];
   float val;
   float err;
};

void results2tree(const char* workDirName, const char* outputFileName, const char* thePoiNames) {
   // workDirName: usual tag where to look for files in Output
   // outFileName: will create a file with this name
   // thePoiNames: comma-separated list of parameters to store ("par1,par2,par3")
   TFile *f = new TFile(outputFileName,"RECREATE");
   TTree *tr = new TTree("fitresults","fit results");


   // bin edges
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   // model names
   Char_t jpsiName[128], psipName[128], bkgName[128];
   // collision system
   Char_t collSystem[8];
   // parameters to store: make it a vector
   vector<poi> thePois;
   TString thePoiNamesStr(thePoiNames);
   TString t; Int_t from = 0;
   while (thePoiNamesStr.Tokenize(t, from , ",")) {
      poi p; strcpy(p.name, t.Data());
      cout << p.name << endl;
      thePois.push_back(p);
   }

   // create tree branches
   tr->Branch("ptmin",&ptmin,"ptmin/F");
   tr->Branch("ptmax",&ptmax,"ptmax/F");
   tr->Branch("ymin",&ymin,"ymin/F");
   tr->Branch("ymax",&ymax,"ymax/F");
   tr->Branch("centmin",&centmin,"centmin/F");
   tr->Branch("centmax",&centmax,"centmax/F");
   tr->Branch("jpsiName",jpsiName,"jpsiName/C");
   tr->Branch("psipName",psipName,"psipName/C");
   tr->Branch("bkgName",bkgName,"bkgName/C");
   tr->Branch("collSystem",collSystem,"collSystem/C");

   for (vector<poi>::iterator it=thePois.begin(); it!=thePois.end(); it++) {
      tr->Branch(Form("%s_val",it->name),&(it->val),Form("%s_val/F",it->name));
      tr->Branch(Form("%s_err",it->name),&(it->err),Form("%s_err/F",it->name));
   }

   // list of files
   vector<TString> theFiles = fileList(workDirName);

   int cnt=0;
   for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      cout << "Parsing file " << cnt << " / " << theFiles.size() << ": " << *it << endl;

      // parse the file name to get info
      anabin thebin = binFromFile(*it);
      ptmin = thebin.ptbin().low();
      ptmax = thebin.ptbin().high();
      ymin = thebin.rapbin().low();
      ymax = thebin.rapbin().high();
      centmin = thebin.centbin().low();
      centmax = thebin.centbin().high();
      strcpy(collSystem, (it->Index("PbPb")) ? "PbPb" : "PP");

      // get the model names
      from = 0;
      bool catchjpsi=false, catchpsip=false, catchbkg=false;
      while (it->Tokenize(t, from, "_")) {
         if (catchjpsi) {strcpy(jpsiName, t.Data()); catchjpsi=false;}
         if (catchpsip) {strcpy(psipName, t.Data()); catchpsip=false;}
         if (catchbkg) {strcpy(bkgName, t.Data()); catchbkg=false;}
         if (t=="Jpsi") catchjpsi=true;
         if (t=="Psi2S") catchpsip=true;
         if (t=="Bkg") catchbkg=true;
      }

      // get the POIs
      for (vector<poi>::iterator itpoi=thePois.begin(); itpoi!=thePois.end(); itpoi++) {
         RooRealVar *thevar = poiFromFile(*it, Form("_%s",collSystem), itpoi->name);
         itpoi->val = thevar ? thevar->getVal() : 0;
         itpoi->err = thevar ? thevar->getError() : 0;
      }

      // fill the tree
      tr->Fill();
      cnt++;
   } // loop on the files

   f->Write();
   f->Close();
}
