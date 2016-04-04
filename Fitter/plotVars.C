#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/bin.h"
#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Systematics/syst.h"
#include "results2tree.C"

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;


/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

// will plot the dependence of varname in workDirName, as a function of pt, centrality or rapidity. 
// collTag should be PP or PbPb. If collTag="", then both PP and PbPb are plotted.
void plotPt(const char* workDirName, const char* varname, const char* collTag="", bool plotErr=true, bool isMC=false);
void plotCent(const char* workDirName, const char* varname, const char* collTag="", bool plotErr=true, bool isMC=false);
void plotRap(const char* workDirName, const char* varname, const char* collTag="", bool plotErr=true, bool isMC=false);

// will plot the dependence of varname as a function to xaxis (=pt, cent or rap) for the file in workDirNames (of the form "dir1,dir2,dir3,...")
void plotFiles(const char* workDirNames, const char* varname, const char* xaxis, float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax, 
      const char* collTag="PP", bool plotErr=true, bool isMC=false);

////////////////////////
// OTHER DECLARATIONS //
////////////////////////

TGraphErrors* plotVar(TTree *tr, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true);
vector<TGraphErrors*> plotVar(TTree *tr, const char* varname, vector<anabin> theBin, string xaxis, string collTag, bool plotErr=true);
TGraphErrors* plotVar(const char* filename, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true);
void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags, const char* workDirName, const char* basename="");



////////////////////
// IMPLEMENTATION //
////////////////////

void plotPt(const char* workDirName, const char* varname, const char* collTag, bool plotErr, bool isMC) {
   string xaxis = "pt";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   TFile *f = new TFile(treeFileName(workDirName,isMC));
   if (!f) {
      results2tree(workDirName,isMC);
      f = new TFile(treeFileName(workDirName,isMC));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg;
   if (string(collTag) != "") tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr);
   else {
      // plot for pp and pbpb and concatenate the results to tg
      vector<TGraphErrors*> tgpp =  plotVar(tr, varname, theCats, xaxis, "PP", plotErr);
      vector<TGraphErrors*> tgpbpb =  plotVar(tr, varname, theCats, xaxis, "PbPb", plotErr);
      tg.insert(tg.end(), tgpp.begin(), tgpp.end());
      tg.insert(tg.end(), tgpbpb.begin(), tgpbpb.end());
   }
   vector<string> tags;
   vector<string> collTags;
   if (string(collTag) != "") collTags.push_back(collTag);
   else {
      collTags.push_back("PP");
      collTags.push_back("PbPb");
   }
   vector<anabin>::const_iterator it;
   vector<string>::const_iterator itc;
   for (itc=collTags.begin(); itc!=collTags.end(); itc++) {
      for (it=theCats.begin(); it!=theCats.end(); it++) {
         ostringstream oss;
         oss.precision(1);
         oss.setf(ios::fixed);
         oss << *itc << ": "
            << it->rapbin().low() << " < |y| < " << it->rapbin().high() << ", " 
            << it->ptbin().low() << " < p_{T} < " << it->ptbin().high() << " GeV/c";
         tags.push_back(oss.str());
      }
   }
   plotGraphs(tg, tags, workDirName);
}

void plotCent(const char* workDirName, const char* varname, const char* collTag, bool plotErr, bool isMC) {
   string xaxis = "cent";
   vector<anabin> theCats;

   // centrality dependence
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   TFile *f = new TFile(treeFileName(workDirName,isMC));
   if (!f) {
      results2tree(workDirName,isMC);
      f = new TFile(treeFileName(workDirName,isMC));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr);
   vector<string> tags;
   vector<anabin>::const_iterator it;
   for (it=theCats.begin(); it!=theCats.end(); it++) {
      ostringstream oss;
      oss.precision(1);
      oss.setf(ios::fixed);
      oss << it->rapbin().low() << " < |y| < " << it->rapbin().high() << ", " 
         << it->centbin().low()/2. << "%-" << it->centbin().high()/2. << "%";
      tags.push_back(oss.str());
   }
   plotGraphs(tg, tags, workDirName);
}

void plotRap(const char* workDirName, const char* varname, const char* collTag, bool plotErr, bool isMC) {
   string xaxis = "rap";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,-200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,-200));

   TFile *f = new TFile(treeFileName(workDirName,isMC));
   if (!f) {
      results2tree(workDirName,isMC);
      f = new TFile(treeFileName(workDirName,isMC));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr);
   vector<string> tags;
   vector<anabin>::const_iterator it;
   for (it=theCats.begin(); it!=theCats.end(); it++) {
      ostringstream oss;
      oss.precision(1);
      oss.setf(ios::fixed);
      oss << it->rapbin().low() << " < |y| < " << it->rapbin().high() << ", " 
         << it->ptbin().low() << " < p_{T} < " << it->ptbin().high() << " GeV/c";
      tags.push_back(oss.str());
   }
   plotGraphs(tg, tags, workDirName);
}

void plotFiles(const char* workDirNames, const char* varname, const char* xaxis, float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax, 
      const char* collTag, bool plotErr, bool isMC) {

   vector<TGraphErrors*> tg;
   vector<string> tags;
   anabin theBin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);

   TString workDirNamesStr(workDirNames);
   TString workDirName; Int_t from = 0;
   while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
      TGraphErrors *tgg = plotVar(treeFileName(workDirName,isMC), varname, theBin, xaxis, collTag, plotErr);
      if (!tgg) {
         results2tree(workDirName,isMC);
         tgg = plotVar(treeFileName(workDirName,isMC), varname, theBin, xaxis, collTag, plotErr);
      }
      if (tgg) {
         tg.push_back(tgg);
         tags.push_back(workDirName.Data());
      }
   }

   ostringstream oss; oss.precision(0); oss.setf(ios::fixed);
   oss << collTag << "_pt" << ptmin*10. << ptmax*10. << "_rap" << rapmin*10. << rapmax*10. << "_cent" << centmin << centmax;

   if (tg.size()>0) plotGraphs(tg, tags, tags[0].c_str(), oss.str().c_str());
}

TGraphErrors* plotVar(TTree *tr, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr) {
   if (!tr) return NULL;
   vector<double> x, ex, y, ey;
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   float val, val_err=0;
   float valmax=0, valmin=0;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);
   if (string(varname)=="nll" || string(varname)=="chi2" || string(varname)=="normchi2") {
      tr->SetBranchAddress(varname,&val);
   } else {
      tr->SetBranchAddress(Form("%s_val",varname),&val);
   }
   if (plotErr) tr->SetBranchAddress(Form("%s_err",varname),&val_err);
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);
      anabin trbin(ymin, ymax, ptmin, ptmax, centmin, centmax);
      if (!binok(theBin, xaxis, trbin)) continue;
      if (string(collSystem) != collTag) continue;

      if (xaxis=="pt") {
         x.push_back((ptmin+ptmax)/2.);
         ex.push_back((ptmax-ptmin)/2.);
      } else if (xaxis=="cent") {
         x.push_back((centmin+centmax)/4.);
         ex.push_back((centmax-centmin)/4.);
      } else { // if (xaxis=="rap") 
         x.push_back((ymin+ymax)/2.);
         ex.push_back((ymax-ymin)/2.);
      }
      y.push_back(val);
      ey.push_back(val_err);

      // min and max
      valmax = max(valmax, (float) 1.4*(val+val_err));
      valmin = min(valmin, (float) (val>0 ? 0 : 1.1*(val-val_err)));
   }

   int n = x.size();

   TString name = Form("%s_%.1f_%.1f_%.1f_%.1f_%i_%i",collTag.c_str(),
         theBin.rapbin().low(),theBin.rapbin().high(),
         theBin.ptbin().low(),theBin.ptbin().high(),
         theBin.centbin().low(),theBin.centbin().high());
   TString hname = "haxes_" + name;
   TString gname = "gr_" + name;
   TGraphErrors *ans = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
   ans->SetName(gname);
   TH1F *haxes=NULL;
   if (xaxis=="pt") {
      haxes = new TH1F(hname,Form(";p_{T} (GeV/c);%s",varname),1,0,30);
   } else if (xaxis=="cent") {
      haxes = new TH1F(hname,Form(";Centrality bin;%s",varname),1,0,100);
   } else { // if (xaxis=="rap")
      haxes = new TH1F(hname,Form(";|y|;%s",varname),1,0,2.4);
   }
   haxes->GetYaxis()->SetLimits(valmin, valmax);
   haxes->GetYaxis()->SetRangeUser(valmin, valmax);
   haxes->GetYaxis()->SetTitleOffset(2);
   ans->SetHistogram(haxes);

   ans->Sort();

   return ans;
}

vector<TGraphErrors*> plotVar(TTree *tr, const char* varname, vector<anabin> theBin, string xaxis, string collTag, bool plotErr) {
   vector<TGraphErrors*> ans;

   vector<anabin>::const_iterator it;
   for (it=theBin.begin(); it!=theBin.end(); it++) {
      ans.push_back(plotVar(tr, varname, *it, xaxis, collTag, plotErr));
   }

   return ans;
}

TGraphErrors* plotVar(const char* filename, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr) {
   TFile *f = new TFile(filename);
   if (!f) return NULL;
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return NULL;
   return plotVar(tr, varname, theBin, xaxis, collTag, plotErr);
}

void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags, const char* workDirName, const char* basename) {
   if (graphs.size() != tags.size()) {
      cout << "Different number of graphs and legends" << endl;
      return;
   }
   if (graphs.size() == 0) return;

   setTDRStyle();
   TCanvas *c1 = new TCanvas("c1","c1",600,600);

   TLegend *tleg = new TLegend(0.18,0.73,0.52,0.89);
   tleg->SetBorderSize(0);
   tleg->SetTextSize(0.03);

   TH1 *haxes = NULL;
   double ymin=0, ymax=0;

   for (unsigned int i=0; i<graphs.size(); i++) {
      graphs[i]->SetLineColor(1+i);
      graphs[i]->SetMarkerColor(1+i);
      graphs[i]->SetMarkerStyle(20+i);
      if (i==0) {
         graphs[i]->Draw("AP");
         haxes = graphs[i]->GetHistogram();
         ymin = haxes->GetYaxis()->GetXmin();
         ymax = haxes->GetYaxis()->GetXmax();
      }
      else {
         graphs[i]->Draw("P");
         ymin = min(ymin, graphs[i]->GetYaxis()->GetXmin());
         ymax = max(ymax, graphs[i]->GetYaxis()->GetXmax());
      }
      tleg->AddEntry(graphs[i],tags[i].c_str(),"LP");
   }

   if (haxes) haxes->GetYaxis()->SetRangeUser(ymin, ymax);
   // c1->Update();

   tleg->Draw();

   string yaxis = haxes->GetYaxis()->GetTitle();
   string xaxis = "rap";
   TString txaxis(haxes->GetXaxis()->GetTitle());
   if (txaxis.Index("Centrality") != kNPOS) xaxis = "cent";
   if (txaxis.Index("p_{T}") != kNPOS) xaxis = "pt";

   gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", workDirName), kTRUE); 
   c1->SaveAs(Form("Output/%s/plot/RESULT/root/plot_%s_%s_vs_%s.root",workDirName, basename, yaxis.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", workDirName), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/png/plot_%s_%s_vs_%s.png",workDirName, basename, yaxis.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", workDirName), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/plot_%s_%s_vs_%s.pdf",workDirName, basename, yaxis.c_str(), xaxis.c_str()));

   // print the result tables
   string xname = "$|y|$";
   if (xaxis=="cent") xname = "Centrality";
   if (xaxis=="pt") xname = "\\pt";
   string yname = latexSafe(yaxis); // make the name safe for LaTeX
   gSystem->mkdir(Form("Output/%s/tex/", workDirName), kTRUE); 
   char texname[2048]; 
   sprintf(texname, "Output/%s/tex/result_%s_%s_vs_%s.tex",workDirName, basename, yaxis.c_str(), xaxis.c_str());
   vector<string> tags_fixed;
   for (vector<string>::const_iterator it=tags.begin(); it!=tags.end(); it++) {
      string tagfixed=latexSafe(*it);
      tags_fixed.push_back(tagfixed);
   }
   bool samesize=true;
   for (unsigned int i=0; i<graphs.size(); i++) {
      if (graphs[i]->GetN() != graphs[0]->GetN()) samesize=false;
   }
   if (samesize) {
      inittex(texname, xname.c_str(), tags_fixed);
      addline(texname,yname,graphs.size());
      printGraph(graphs, texname);
   } else {
      inittex(texname, xname.c_str(), yname);
      for (unsigned int i=0; i<graphs.size(); i++) {
         addline(texname,tags_fixed[i]);
         printGraph(graphs[i], texname);
      }
   }
   closetex(texname);
   cout << "Closed " << texname << endl;
   cout << "It is advised that you check the contents of " << texname << " as it may not compile nicely as is." << endl;
}
