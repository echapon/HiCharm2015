#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/resultUtils.h"

#include "TString.h"
#include "RooRealVar.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

const char* poiname = "RFrac2Svs1S";
const bool  doratio = true; // true -> look for separate PP and PbPb files, false -> input files are with simultaneous pp-PbPb fits


//////////////////
// DECLARATIONS //
//////////////////

// methods for computing the systematic
double rms(vector<double> v);
double maxdiff(vector<double> v);
// other
RooRealVar* poiFromFile(const char* filename, const char* token="");

///////////////////
// MAIN FUNCTION //
///////////////////

void results2syst(const char* workDirNames, const char* systFileName, const char* systHeader, int method) {
// workDirNames: of the form "dir1,dir2,dir3,..."
// systFileName: "someDir/syst_blabla.csv"
// systHeader: this will be the header of the systematics file. A few words describing what this systematic is.
// method: 0 -> RMS, 1 -> max difference to the first work dir (= nominal)

   map<anabin, vector<double> > mapvals;
   TString workDirNamesStr(workDirNames);
   TString workDirName; Int_t from = 0;
   while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
      // list of files
      vector<TString> theFiles, theFiles2;
      if (!doratio) {
         theFiles = fileList(workDirName.Data());
      } else {
         theFiles = fileList(workDirName.Data(),"PbPb");
         theFiles2 = fileList(workDirName.Data(),"PP");
      }

      cout << "Now parsing files in " << workDirName << endl;

      RooRealVar* theVar=NULL;
      vector<TString>::const_iterator it,it2;
      for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
         anabin thebin = binFromFile(it->Data());
         if (!doratio) {
            theVar = poiFromFile(it->Data(),"_PbPbvsPP");
         } else {
            RooRealVar *num = poiFromFile(it->Data(),"_PbPb");

            // force the centrality to 0-200 for pp
            anabin thebinpp = thebin;
            binI centbin(0,200);
            thebinpp.setcentbin(centbin);

            // look for the pp file corresponding to the pbpb one
            bool found=false;
            for (it2=theFiles2.begin(); it2!=theFiles2.end(); it2++) {
               if (binFromFile(it2->Data()) == thebinpp) {
                  found=true;
                  break;
               }
            }
            if (!found) {
               cout << "Error! I did not find the PP file for " << *it << endl;
               return;
            }
            RooRealVar *den = poiFromFile(it2->Data(),"_PP");
            theVar = ratioVar(num,den);
         }
         if (theVar) mapvals[thebin].push_back(theVar->getVal());
      }
   }

   cout << "Done parsing files. Now computing the systematic." << endl;

   ofstream file(systFileName);
   file << systHeader << endl;
   map<anabin, vector<double> >::const_iterator it;
   for (it=mapvals.begin(); it!=mapvals.end(); it++) {
      anabin thebin = it->first;
      vector<double> v = it->second;
      double value;
      if (method==0) value=rms(v);
      else if (method==1) value=maxdiff(v);
      else value=0;
      file << thebin.rapbin().low() << ", " << thebin.rapbin().high() << ", " 
         << thebin.ptbin().low() << ", " << thebin.ptbin().high() << ", " 
         << thebin.centbin().low() << ", " << thebin.centbin().high() << ", " 
         << value << endl;
   }
   file.close();
   cout << "Closed " << systFileName << endl;
}

double rms(vector<double> v) {
   if (v.size()==0 || v[0]==0) return 0;
    double s=0,s2=0;
    for (unsigned int i=0; i<v.size(); i++) {
       s+=v[i];
       s2+=v[i]*v[i];
    }
    return sqrt(s2-(s*s))/v[0];
 }

double maxdiff(vector<double> v) {
   if (v.size()==0 || v[0]==0) return 0;
   double maxdiff=0;
    for (unsigned int i=1; i<v.size(); i++) {
       maxdiff=max(maxdiff,fabs(v[i]-v[0]));
    }
    return maxdiff/v[0];
}

RooRealVar* poiFromFile(const char* filename, const char* token) {
   return poiFromFile(filename,token,poiname);
}
