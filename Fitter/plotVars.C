#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/bin.h"
#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Systematics/syst.h"

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLegend.h"
#include "TCanvas.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

//////////////////
// DECLARATIONS //
//////////////////

TGraphErrors* plotVar(TTree *tr, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true);
vector<TGraphErrors*> plotVar(TTree *tr, const char* varname, vector<anabin> theBin, string xaxis, string collTag, bool plotErr=true);
TGraphErrors* plotVar(const char* filename, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true);
void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags);


/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

void plotPt(const char* filename, const char* varname, const char* collTag="", bool plotErr=true);
void plotCent(const char* filename, const char* varname, const char* collTag="", bool plotErr=true);
void plotRap(const char* filename, const char* varname, const char* collTag="", bool plotErr=true);

void plotPt(const char* filename, const char* varname, const char* collTag, bool plotErr) {
   string xaxis = "pt";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   TFile *f = new TFile(filename);
   if (!f) return;
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
   plotGraphs(tg, tags);
}

void plotCent(const char* filename, const char* varname, const char* collTag, bool plotErr) {
   string xaxis = "cent";
   vector<anabin> theCats;

   // centrality dependence
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   TFile *f = new TFile(filename);
   if (!f) return;
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
   plotGraphs(tg, tags);
}

void plotRap(const char* filename, const char* varname, const char* collTag, bool plotErr) {
   string xaxis = "rap";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,-200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,-200));

   TFile *f = new TFile(filename);
   if (!f) return;
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
   plotGraphs(tg, tags);
}


////////////////////
// IMPLEMENTATION //
////////////////////

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
   tr->SetBranchAddress(Form("%s_val",varname),&val);
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
      valmax = max(valmax, (float) 1.1*(val+val_err));
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

void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags) {
   if (graphs.size() != tags.size()) {
      cout << "Different number of graphs and legends" << endl;
      return;
   }

   setTDRStyle();
   TCanvas *c1 = new TCanvas("c1","c1",600,600);

   TLegend *tleg = new TLegend(0.19,0.73,0.53,0.89);
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
         cout << ymin << " " << ymax << endl;
      }
      else {
         graphs[i]->Draw("P");
         ymin = min(ymin, graphs[i]->GetYaxis()->GetXmin());
         ymax = max(ymax, graphs[i]->GetYaxis()->GetXmax());
         cout << ymin << " " << ymax << endl;
      }
      tleg->AddEntry(graphs[i],tags[i].c_str(),"LP");
   }

   if (haxes) haxes->GetYaxis()->SetRangeUser(ymin, ymax);
   // c1->Update();

   tleg->Draw();

   c1->SaveAs("plot.pdf");
   c1->SaveAs("plot.png");
   c1->SaveAs("plot.root");
}
