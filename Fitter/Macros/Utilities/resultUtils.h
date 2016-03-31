#ifndef resultUtils_h
#define resultUtils_h

#include "bin.h"

#include "RooRealVar.h"
#include "TString.h"
#include "TFile.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

RooRealVar* poiFromFile(const char* filename, const char* token, const char* thepoiname);
vector<TString> fileList(const char* input, const char* token="");
RooRealVar* ratioVar(RooRealVar *num, RooRealVar *den, bool usedenerror=true);
anabin binFromFile(const char* filename);
bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck);

RooRealVar* poiFromFile(const char* filename, const char* token, const char* thepoiname) {
   TFile *f = new TFile(filename);
   if (!f) {
      cout << "Error, file " << filename << " does not exist." << endl;
      return NULL;
   }
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   if (!ws) {
      cout << "Error, file " << filename << " is bad." << endl;
      return NULL;
   }
   TString poiname_and_token = TString(thepoiname) + TString(token);
   RooRealVar *ans = (RooRealVar*) ws->var(poiname_and_token);
   if (!ans) return NULL;
   RooRealVar* ansc = new RooRealVar(*ans,poiname_and_token + Form("_from_%s",filename));
   f->Close(); delete f;
   return ansc;
}

vector<TString> fileList(const char* input, const char* token) {
   vector<TString> ans;

   TString basedir(Form("Output/%s/result/DATA/",input));
   TSystemDirectory dir(input,basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".root") && (TString(token) == "" || fname.Index(token) != kNPOS)) {
            ans.push_back(basedir+fname);
         }
      }
   }

   return ans;
}

RooRealVar* ratioVar(RooRealVar *num, RooRealVar *den, bool usedenerror) {
   double n = num->getVal();
   double d = den->getVal();
   double dn = num->getError();
   double dd = den->getError();

   double r = d!=0 ? n/d : 0;
   double dr = n!=0 && d!=0 ? r * sqrt(pow(dn/n,2) + pow(dd/d,2)) : 0;
   if (!usedenerror && n!=0) dr = (dn/n)*r;
   RooRealVar *ans = new RooRealVar(Form("%s_over_%s",num->GetName(),den->GetName()), Form("%s / %s",num->GetTitle(),den->GetTitle()), r);
   ans->setError(dr);

   return ans;
}

anabin binFromFile(const char* filename) {
   TFile *f = new TFile(filename);
   if (!f) {
      cout << "Error, file " << filename << " does not exist." << endl;
      return anabin(0,0,0,0,0,0);
   }
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   if (!ws) {
      cout << "Error, file " << filename << " is bad." << endl;
      return anabin(0,0,0,0,0,0);
   }
   RooRealVar *pt = (RooRealVar*) ws->var("pt");
   RooRealVar *rap = (RooRealVar*) ws->var("rap");
   RooRealVar *cent = (RooRealVar*) ws->var("cent");
   if (!pt || !rap || !cent) {
      cout << "Error, file " << filename << " is bad." << endl;
      return anabin(0,0,0,0,0,0);
   }
   anabin ans(rap->getMin(),rap->getMax(),pt->getMin(),pt->getMax(),cent->getMin(),cent->getMax());
   f->Close(); delete f;
   return ans;
}

bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck) {
   bool ok=false;

   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      if (xaxis=="pt" && it->rapbin()==tocheck.rapbin() && it->centbin()==tocheck.centbin()
            && ! (it->ptbin()==tocheck.ptbin())) {
         ok=true;
         tocheck.setptbin(it->ptbin());
         break;
      } else if (xaxis=="cent" && it->rapbin()==tocheck.rapbin() && it->ptbin()==tocheck.ptbin()
            && ! (it->centbin()==tocheck.centbin())) {
         ok=true;
         tocheck.setcentbin(it->centbin());
         break;
      }
   }

   return ok;
}


#endif // ifndef resultUtils_h
