#ifndef syst_h
#define syst_h

#include "../Macros/Utilities/bin.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include <map>
#include <string>
#include <fstream>
#include <iostream>

struct syst {
   string name;
   double value;
};

using namespace std;

vector<TString> fileList_syst(const char* token) {
   vector<TString> ans;

   TString basedir("Systematics/csv/");
   TSystemDirectory dir(token,basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".csv") && (TString(token) == "" || fname.Index(token) != kNPOS)) {
            ans.push_back(basedir+fname);
         }
      }
   }

   return ans;
};


map<anabin, syst> readSyst(const char* systfile) {
   map<anabin, syst> ans;

   ifstream file(systfile);
   if (!(file.good())) return ans;

   string systname; getline(file,systname);

   string line;
   double rapmin, rapmax, ptmin, ptmax, centmin, centmax, value;

   while (file.good()) {
      getline(file,line);
      if (line.size()==0) break;
      TString tline(line.c_str());
      TString t; Int_t from = 0, cnt=0;
      while (tline.Tokenize(t, from , ",")) {
         t.Strip(TString::kBoth,' ');
         value = atof(t.Data());
         if (cnt==0) rapmin = atof(t.Data());
         else if (cnt==1) rapmax = value;
         else if (cnt==2) ptmin = value;
         else if (cnt==3) ptmax = value;
         else if (cnt==4) centmin = value;
         else if (cnt==5) centmax = value;
         else if (cnt>6) {
            cout << "Warning, too many fields, I'll take the last one." << endl;
            continue;
         }
         cnt++;
      }
      anabin thebin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);
      syst thesyst; thesyst.value = value; thesyst.name = systname;
      ans[thebin] = thesyst;
   }

   return ans;
};

map<anabin, syst> combineSyst(vector< map<anabin, syst> > theSysts, string name="total") {
   map<anabin, syst> ans;

   vector< map<anabin, syst> >::const_iterator it;
   for (it=theSysts.begin(); it!=theSysts.end(); it++) {
      map<anabin, syst>::const_iterator it2;
      for (it2=it->begin(); it2!=it->end(); it2++) {
         anabin thebin = it2->first;
         syst thesyst = it2->second;
         thesyst.name = name;

         // if we already have a syst for this bin, sum quadractically the existing syst and the new syst
         if (ans.find(thebin) != ans.end()) thesyst.value = sqrt(pow(thesyst.value,2) + pow(ans[thebin].value,2));
         ans[thebin] = thesyst;
      }
   }

   return ans;
};

map<anabin, syst> readSyst_all(const char* token) {
   // token should be PP or PbPb

   vector< map<anabin, syst> > systmap_all;
   vector<TString> filelist = fileList_syst(token);

   for (vector<TString>::const_iterator it=filelist.begin(); it!=filelist.end(); it++) {
      cout << "Reading file " << *it << endl;
      map<anabin,syst> systmap = readSyst(it->Data());
      systmap_all.push_back(systmap);
   }

   map<anabin,syst> ans = combineSyst(systmap_all,token);
   return ans;
}


#endif // ifndef syst_h
