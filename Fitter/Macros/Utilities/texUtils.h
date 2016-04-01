#ifndef texutils_h
#define texutils_h

#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include <vector>
#include <string>
#include <fstream>

using namespace std;

void printHist(vector<TH1F*> hist, const char* filename);
void printHist(TH1F* hist, const char* filename);
void printGraph(vector<TGraphAsymmErrors*> tg, const char* filename);
void printGraph(TGraphAsymmErrors* tg, const char* filename);
void printGraph(vector<TGraphAsymmErrors*> tg, vector<TGraphAsymmErrors*> tg_syst, const char* filename);
void printGraph(TGraphAsymmErrors* tg, TGraphAsymmErrors* tg_syst, const char* filename);
void inittex(const char* filename, const char* xname, vector<string> yname);
void inittex(const char* filename, const char* xname, string yname);
void closetex(const char* filename);
void addline(const char* filename, string line, int n=1);

void printHist(vector<TH1F*> hist, const char* filename) {
   ofstream file(filename, ofstream::app);
   file.precision(3);
   if (hist.size()==0) return;
   int nbins = hist[0]->GetNbinsX();
   for (int ibin=1; ibin<nbins+1; ibin++) {
      file << "$[" << hist[0]->GetXaxis()->GetBinLowEdge(ibin) << "-" << hist[0]->GetXaxis()->GetBinUpEdge(ibin) << "]$";
      file.setf(ios::fixed);
      for (vector<TH1F*>::const_iterator ith=hist.begin(); ith!=hist.end(); ith++) {
         file << " & $" << (*ith)->GetBinContent(ibin) << " \\pm " << (*ith)->GetBinError(ibin) << "$";
      }
      file.unsetf(ios::fixed);
      file << "\\\\" << endl;
   }
   file << "\\hline" << endl;
   file.close();
}

void printHist(TH1F* hist, const char* filename) {
   vector<TH1F*> v; v.push_back(hist);
   printHist(v, filename);
}

void printGraph(vector<TGraphAsymmErrors*> tg, const char* filename) {
   ofstream file(filename, ofstream::app);
   file.precision(3);
   if (tg.size()==0) return;
   int nbins = tg[0]->GetN();
   for (int ibin=0; ibin<nbins; ibin++) {
      file << "$[" << tg[0]->GetX()[ibin]-tg[0]->GetErrorXlow(ibin) << "-" << tg[0]->GetX()[ibin]+tg[0]->GetErrorXhigh(ibin) << "]$";
      file.setf(ios::fixed);
      for (vector<TGraphAsymmErrors*>::const_iterator itg=tg.begin(); itg!=tg.end(); itg++) {
         double errlow = (*itg)->GetErrorYlow(ibin);
         double errhigh = (*itg)->GetErrorYhigh(ibin);
         if (fabs(errlow-errhigh)>1e-3) {
            file << " & $" << (*itg)->GetY()[ibin] << "_{-" << errlow << "}^{+" << errhigh << "} $";
         } else {
            file << " & $" << (*itg)->GetY()[ibin] << " \\pm " << max(errlow,errhigh) << "$";
         }
      }
      file.unsetf(ios::fixed);
      file << "\\\\" << endl;
   }
   file << "\\hline" << endl;
   file.close();
}

void printGraph(TGraphAsymmErrors* tg, const char* filename) {
   vector<TGraphAsymmErrors*> v; v.push_back(tg);
   printGraph(v, filename);
}

void printGraph(vector<TGraphAsymmErrors*> tg, vector<TGraphAsymmErrors*> tg_syst, const char* filename) {
   ofstream file(filename, ofstream::app);
   file.precision(3);
   if (tg.size()==0) return;
   if (tg.size() != tg_syst.size()) return;
   int nbins = tg[0]->GetN();
   for (int ibin=0; ibin<nbins; ibin++) {
      file << "$[" << tg[0]->GetX()[ibin]-tg[0]->GetErrorXlow(ibin) << "-" << tg[0]->GetX()[ibin]+tg[0]->GetErrorXhigh(ibin) << "]$";
      file.setf(ios::fixed);
      vector<TGraphAsymmErrors*>::const_iterator itg=tg.begin(), itg_syst=tg_syst.begin();
      for (; itg!=tg.end(); itg++) {
         double errlow = (*itg)->GetErrorYlow(ibin);
         double errhigh = (*itg)->GetErrorYhigh(ibin);
         double errlow_syst = (*itg_syst)->GetErrorYlow(ibin);
         double errhigh_syst = (*itg_syst)->GetErrorYhigh(ibin);
         if (fabs(errlow-errhigh)>1e-3) {
            file << " & $" << (*itg)->GetY()[ibin] << "_{-" << errlow << "}^{+" << errhigh << "}\\textrm{ (stat.)} ";
         } else {
            file << " & $" << (*itg)->GetY()[ibin] << " \\pm" << max(errlow,errhigh) << "\\textrm{ (stat.)} ";
         }
         if (fabs(errlow_syst-errhigh_syst)>1e-3) {
            file << "_{-" << errlow_syst << "}^{+" << errhigh_syst << "}\\textrm{ (syst.)} $";
         } else {
            file << " \\pm " << max(errlow_syst,errhigh_syst) << "\\textrm{ (syst.)} $";
         }
         itg_syst++;
      }
      file.unsetf(ios::fixed);
      file << "\\\\" << endl;
   }
   file << "\\hline" << endl;
   file.close();
}

void printGraph(TGraphAsymmErrors* tg, TGraphAsymmErrors* tg_syst, const char* filename) {
   vector<TGraphAsymmErrors*> v, vs;
   v.push_back(tg);
   vs.push_back(tg_syst);
   printGraph(v,vs,filename);
}

void inittex(const char* filename, const char* xname, vector<string> yname) {
   ofstream file(filename);
   file << "\\begin{tabular}{c"; 
   for (unsigned int i=0; i<yname.size(); i++) file << "c";
   file << "}" << endl;
   file << "\\hline" << endl;
   file << xname;
   for (unsigned int i=0; i<yname.size(); i++) file << " & " << yname[i];
   file<< "\\\\" << endl;
   file << "\\hline" << endl;
   file.close();
}

void inittex(const char* filename, const char* xname, string yname) {
   vector<string> v; v.push_back(yname);
   inittex(filename, xname, v);
}

void closetex(const char* filename) {
   ofstream file(filename, ofstream::app);
   file << "\\end{tabular}" << endl;
   file.close();
}

void addline(const char* filename, string line, int n) {
   ofstream file(filename, ofstream::app);
   file << "& \\multicolumn{" << n << "}{c}{" << line << "}\\\\" << endl;
   file << "\\hline" << endl;
   file.close();
}

#endif // ifndef texutils_h
