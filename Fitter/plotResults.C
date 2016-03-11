#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"

#include "Macros/Utilities/initClasses.h"

#include <vector>
#include <map>
#include <string>
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace std;

const char* ylabel = "(#Psi(2S)/J/#Psi)_{PbPb} / (#Psi(2S)/J/#Psi)_{pp}";
const char* poiname = "RFrac2Svs1S_PbPbvsPP";

// function to get the param of interest from a workspace
RooRealVar* poiFromFile(const char* filename);
// function to get the analysis bin from a file
anabin binFromFile(const char* filename);
// check
bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck);
// plot
void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, string outputDir);
void plot(vector<anabin> thecats, string xaxis, string outputDir);

void plotPt(string outputDir) {
   string xaxis = "pt";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   plot(theCats,xaxis,outputDir);
};

void plotCent(string outputDir) {
   string xaxis = "cent";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   plot(theCats,xaxis,outputDir);
};


void plot(vector<anabin> thecats, string xaxis, string outputDir) {
   // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
   // xaxis is the variable to be plotted. "pt", "rap" or "cent"

   // list of files
   vector<const char*> theFiles;
   theFiles.push_back("");
   theFiles.push_back("");

   map<anabin, RooRealVar*> theVars;

   for (vector<const char*>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      anabin thebin = binFromFile(*it);
      theVars[thebin] = poiFromFile(*it);
   }

   map<anabin, vector<anabin> > theBins;
   map<anabin, vector<RooRealVar*> > theVarsBinned;
   map<anabin, TGraphAsymmErrors* > theGraphs;

   // initialize the maps
   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      theBins[*it] = vector<anabin>();
      theVarsBinned[*it] = vector<RooRealVar*>();
   }

   for (map<anabin, RooRealVar*>::const_iterator it=theVars.begin(); it!=theVars.end(); it++) {
      anabin thebin = it->first;
      if (!binok(thecats,xaxis,thebin)) continue;
      theBins[thebin].push_back(it->first);
      theVarsBinned[thebin].push_back(it->second);
   }

   // make TGraphAsymmErrors
   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      int n = theBins[*it].size();
      if(n==0) {
         cout << "Error, nothing found for category" << endl;
         theGraphs[*it] = NULL;
         continue;
      }
      theGraphs[*it] = new TGraphAsymmErrors(n);
      for (int i=0; i<n; i++) {
         double x, exl, exh, y, eyl, eyh;
         double low = theBins[*it][i].ptbin.low;
         double high = theBins[*it][i].ptbin.high;
         if (xaxis=="pt") {
            x = (low+high)/2.;
            exh = (high-low)/2.;
            exl = (high-low)/2.;
         }
         if (xaxis=="cent") {
            x = HI::findNcollAverage(low,high);
            exl = 0.;
            exh = 0.;
         }
         y = theVarsBinned[*it][i]->getVal();
         eyl = theVarsBinned[*it][i]->getErrorLo();
         eyh = theVarsBinned[*it][i]->getErrorHi();
         theGraphs[*it]->SetPoint(i,x,y);
         theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
      }
   }

   // plot
   plotGraph(theGraphs, xaxis, outputDir);
}

anabin binFromFile(const char* filename) {
   TFile *f = new TFile(filename);
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   RooRealVar *pt = (RooRealVar*) ws->var("pt");
   RooRealVar *rap = (RooRealVar*) ws->var("rap");
   RooRealVar *cent = (RooRealVar*) ws->var("cent");
   anabin ans(pt->getMin(),pt->getMax(),rap->getMin(),rap->getMax(),cent->getMin(),cent->getMax());
   f->Close(); delete f;
   return ans;
}

RooRealVar* poiFromFile(const char* filename) {
   TFile *f = new TFile(filename);
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   RooRealVar *ans = (RooRealVar*) ws->var(poiname);
   RooRealVar* ansc = new RooRealVar(*ans,Form("%s_from_%s",poiname,filename));
   f->Close(); delete f;
   return ansc;
}

bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck) {
   bool ok=false;

   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      if (xaxis=="pt" && it->rapbin==tocheck.rapbin && it->centbin==tocheck.centbin) {
         ok=true;
         tocheck.ptbin = it->ptbin;
         break;
      } else if (xaxis=="cent" && it->rapbin==tocheck.rapbin && it->ptbin==tocheck.ptbin) {
         ok=true;
         tocheck.centbin = it->centbin;
         break;
      }
   }

   return ok;
}

void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, string outputDir) {
   TCanvas *c1 = new TCanvas("c1","c1",600,600);

   // the axes
   TH1F *haxes; TLine line;
   if (xaxis=="pt") {
      haxes = new TH1F("haxes","haxes",1,0,30);
      line = TLine(0,1,30,1);
   }
   if (xaxis=="cent") {
      haxes = new TH1F("haxes","haxes",1,0,420);
      line = TLine(0,1,420,1);
   }
   haxes->GetYaxis()->SetRangeUser(0,1.5);
   haxes->GetYaxis()->SetTitle(ylabel);
   const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "N_{part}";
   haxes->GetXaxis()->SetTitle(xlabel);
   haxes->Draw();
   line.Draw();

   TLegend *tleg = new TLegend(0.4,0.4,0.6,0.6);

   int cnt=0;
   for (map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin(); it!=theGraphs.end(); it++) {
      TGraphAsymmErrors* tg = it->second;
      if (!tg) continue;

      if (cnt==0) {
         tg->SetMarkerStyle(kFullSquare);
         tg->SetMarkerColor(kRed);
         tg->SetLineColor(kRed);
      } else if (cnt==1) {
         tg->SetMarkerStyle(kFullCircle);
         tg->SetMarkerColor(kBlue);
         tg->SetLineColor(kBlue);
      } else if (cnt==2) {
         tg->SetMarkerStyle(kFullTriangleUp);
         tg->SetMarkerColor(kGreen);
         tg->SetLineColor(kGreen);
      }
      tg->Draw("AP");      

      TString raplabel = Form("%.1f < |y| < %.1f ; ",it->first.rapbin.low,it->first.rapbin.high);
      TString otherlabel = "BWAA";
      if (xaxis == "pt") otherlabel = Form("%i%-%i%",(int) (it->first.centbin.low/2.), (int) (it->first.centbin.high/2.));
      tleg->AddEntry(tg, raplabel + otherlabel, "p");

      cnt++;
   }

   tleg->Draw();

   int iPos = 33;
   CMS_lumi( (TPad*) gPad, 99, iPos, "Preliminary pp " );

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   gSystem->mkdir(Form("%splot/RESULT/root/", outputDir.c_str()), kTRUE); 
   c1->SaveAs(Form("%splot/RESULT/root/result_%s.root",outputDir.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("%splot/RESULT/png/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("%splot/RESULT/png/result_%s.png",outputDir.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("%splot/RESULT/pdf/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("%splot/RESULT/pdf/result_%s.pdf",outputDir.c_str(), xaxis.c_str()));

   delete c1;
}
