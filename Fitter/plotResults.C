#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"

#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Systematics/syst.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

using namespace std;

////////////////
// PARAMETERS //
////////////////

const char* ylabel = "(#Psi(2S)/J/#Psi)_{PbPb} / (#Psi(2S)/J/#Psi)_{pp}";
const char* poiname = "RFrac2Svs1S";
const bool  doratio = true; // true -> look for separate PP and PbPb files, false -> input files are with simultaneous pp-PbPb fits


///////////////
// CONSTANTS //
///////////////

const double xfrac = 0.8; // for the 2-panel plot for the centrality dependence.


//////////////////
// DECLARATIONS //
//////////////////

RooRealVar* poiFromFile(const char* filename, const char* token="");
// plot
void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst);
void plot(vector<anabin> thecats, string xaxis, string workDirName);
void centrality2npart(TGraphAsymmErrors* tg, bool issyst=false, bool isMB=false, double xshift=0.);



/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

void plotPt(string workDirName) {
   string xaxis = "pt";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   plot(theCats,xaxis,workDirName);
};

void plotCent(string workDirName) {
   string xaxis = "cent";
   vector<anabin> theCats;

   // centrality dependence
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   // for minimum bias
   theCats.push_back(anabin(0,1.6,6.5,30,0,-200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,-200));

   plot(theCats,xaxis,workDirName);
};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plot(vector<anabin> thecats, string xaxis, string outputDir) {
   // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
   // xaxis is the variable to be plotted. "pt", "rap" or "cent"

   // list of files
   vector<TString> theFiles, theFiles2;
   if (!doratio) {
      theFiles = fileList(outputDir.c_str());
   } else {
      theFiles = fileList(outputDir.c_str(),"PbPb");
      theFiles2 = fileList(outputDir.c_str(),"PP");

      // if (theFiles.size() != theFiles2.size()) {
      //    cout << "I found a different number of pp and PbPb files!! This is no good, I give up." << endl;
      //    return;
      // }
   }

   map<anabin, RooRealVar*> theVars;
   map<anabin,syst> stat_PP; // statistical uncertainty on PP fits (useful for the case of the centrality dependence)

   vector<TString>::const_iterator it,it2;
   for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      anabin thebin = binFromFile(it->Data());
      if (!doratio) {
         theVars[thebin] = poiFromFile(it->Data(),"_PbPbvsPP");
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
         theVars[thebin] = ratioVar(num,den,!(xaxis=="cent")); // if plotting the centrality dependence, do not put the pp stat
         syst thestat_PP;
         thestat_PP.name = "stat_PP";
         thestat_PP.value = den->getVal() != 0 ? den->getError()/den->getVal() : 0;
         stat_PP[thebinpp] = thestat_PP;
      }
   }

   map<anabin, vector<anabin> > theBins;
   map<anabin, vector<RooRealVar*> > theVarsBinned;
   map<anabin, TGraphAsymmErrors* > theGraphs;
   map<anabin, TGraphAsymmErrors* > theGraphs_syst;

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

   // systematics
   map<anabin, syst> syst_PP = readSyst_all("PP");
   if (doratio && xaxis=="cent") { // put the PP stat error into the PP syst, that will go into a box at 1
      vector< map<anabin, syst> > all_PP;
      all_PP.push_back(syst_PP);
      all_PP.push_back(stat_PP);
      syst_PP = combineSyst(all_PP,"statsyst_PP");
   }
   map<anabin, syst> syst_PbPb = readSyst_all("PbPb");

   // make TGraphAsymmErrors
   int cnt=0;
   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      int n = theBins[*it].size();
      if(n==0) {
         cout << "Error, nothing found for category" << endl;
         theGraphs[*it] = NULL;
         continue;
      }
      bool isMB = (it->centbin().low()<=0 && it->centbin().high()<=0);
      if (isMB && n!=1) cout << "Warning, I have " << n << " bins for a MinBias category, I expected 1." << endl;

      theGraphs[*it] = new TGraphAsymmErrors(n);
      theGraphs[*it]->SetName(Form("bin_%i",cnt));
      theGraphs_syst[*it] = new TGraphAsymmErrors(n);
      theGraphs_syst[*it]->SetName(Form("bin_%i_syst",cnt));

      for (int i=0; i<n; i++) {
         double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
         double exsyst=0, eysyst=0;
         double low=0, high=0; 
         anabin thebin = theBins[*it][i];
         if (xaxis=="pt") {
            low= thebin.ptbin().low();
            high = thebin.ptbin().high();
            x = (low+high)/2.;
            exh = (high-low)/2.;
            exl = (high-low)/2.;
            exsyst = 0.5;
            eysyst = sqrt(pow(syst_PP[thebin].value,2) + pow(syst_PbPb[thebin].value,2)); // quadratic sum of PP and PbPb systs
         }
         if (xaxis=="cent") {
            low= thebin.centbin().low();
            high = thebin.centbin().high();
            // x = isMB ? 150 + (150./1.6)*it->rapbin().low() : HI::findNpartAverage(low,high);
            // exl = 0.;
            // exh = 0.;
            x = (low+high)/2./2.;
            exh = (high-low)/2./2.;
            exl = (high-low)/2./2.;
            exsyst = !isMB ? 5 : 5./(1.-xfrac);
            eysyst = syst_PbPb[thebin].value; // only PbPb syst: the PP one will go to a dedicated box
            // also add
         }
         y = theVarsBinned[*it][i]->getVal();
         eyl = fabs(theVarsBinned[*it][i]->getErrorLo());
         eyh = theVarsBinned[*it][i]->getErrorHi();
         eysyst = y*eysyst;
         theGraphs[*it]->SetPoint(i,x,y);
         theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
         theGraphs_syst[*it]->SetPoint(i,x,y);
         theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst,eysyst);
         cout << x << " " << y << " " << eyl << " " << eyh << " " << eysyst << endl;
      }
      cnt++;
   }

   // plot
   plotGraph(theGraphs, theGraphs_syst, xaxis, outputDir, syst_PP);
}

RooRealVar* poiFromFile(const char* filename, const char* token) {
   return poiFromFile(filename,token,poiname);
}


void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst) {
   setTDRStyle();


   TCanvas *c1 = NULL;
   if (xaxis=="cent") c1 = new TCanvas("c1","c1",600/xfrac,600);
   else c1 = new TCanvas("c1","c1",600,600);

   // in the case of the centrality dependence, we need the minimum bias panel on the right
   TPad *padl=NULL, *padr=NULL;
   if (xaxis=="cent") {
         padl = new TPad("padl_left","",0,0,xfrac,1);
         padl->SetLeftMargin(gStyle->GetPadLeftMargin()/xfrac);
         padl->SetRightMargin(0);
         // padl->SetBottomMargin(gStyle->GetPadBottomMargin());
         // padl->SetTopMargin(gStyle->GetPadTopMargin());
         padl->SetFrameBorderMode(0);
         padl->SetBorderMode(0);
         padl->SetBorderSize(0);
         padl->Draw();

         padr = new TPad("padr_right","",xfrac,0,1,1);
         padr->SetRightMargin(gStyle->GetPadRightMargin()/(1.-xfrac));
         padr->SetLeftMargin(0);
         // padr->SetBottomMargin(gStyle->GetPadBottomMargin());
         // padr->SetTopMargin(gStyle->GetPadTopMargin());
         padr->SetFrameBorderMode(0);
         padr->SetBorderMode(0);
         padr->SetBorderSize(0);
         padr->Draw();

         padl->cd();
   }

   // the axes
   TH1F *haxes=NULL, *haxesr=NULL; TLine line, liner;
   if (xaxis=="pt") {
      haxes = new TH1F("haxes","haxes",1,0,30);
      line = TLine(0,1,30,1);
   }
   if (xaxis=="cent") {
      haxes = new TH1F("haxesl","haxesl",1,0,420);
      haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X")/xfrac);
      // haxes->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")*xfrac);
      line = TLine(0,1,420,1);
      haxesr = new TH1F("haxesr","haxesr",1,0,420);
      haxesr->GetXaxis()->SetTickLength(0);
      haxesr->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/(1.-xfrac));
      haxesr->GetYaxis()->SetRangeUser(0,1.5);
      haxesr->GetXaxis()->SetTitleSize(0);
      haxesr->GetXaxis()->SetLabelSize(0);
      liner = TLine(0,1,420,1);
   }
   haxes->GetYaxis()->SetRangeUser(0,1.5);
   haxes->GetYaxis()->SetTitle(ylabel);
   const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "N_{part}";
   haxes->GetXaxis()->SetTitle(xlabel);
   haxes->Draw();
   line.Draw();
   if (xaxis=="cent") {
      padr->cd();
      haxesr->Draw();
      liner.Draw();
      padl->cd();
   }

   double xshift=0.;
   if (xaxis=="cent") xshift=0.05;
   TLegend *tleg = new TLegend(0.16+xshift,0.73,0.50,0.89);
   tleg->SetBorderSize(0);
   tleg->SetTextSize(0.03);

   // prepare for the printing of the result tables
   const char* xname = (xaxis=="cent") ? "Centrality" : "\\pt";
   gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE); 
   char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s.tex",outputDir.c_str(),xaxis.c_str());
   string yname("\\doubleratio");
   inittex(texname, xname, yname);

   int cnt=0;
   map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
   map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
   for (; it!=theGraphs.end(); it++) {
      anabin thebin = it->first;
      bool isMB = (thebin.centbin().low()<=0 && thebin.centbin().high()<=0);
      TGraphAsymmErrors* tg = it->second;
      TGraphAsymmErrors* tg_syst = it_syst->second;
      if (!tg || !tg_syst) continue;

      if (thebin.rapbin() == binD(0.,1.6)) {
         tg->SetMarkerStyle(kFullSquare);
         tg->SetMarkerColor(kRed);
         tg->SetLineColor(kRed);
         tg_syst->SetFillColorAlpha(kRed, 0.5);
      } else if (thebin.rapbin() == binD(1.6,2.4)) {
         tg->SetMarkerStyle(kFullCircle);
         tg->SetMarkerColor(kBlue);
         tg->SetLineColor(kBlue);
         tg_syst->SetFillColorAlpha(kBlue, 0.5);
      } else {
         tg->SetMarkerStyle(kFullTriangleUp);
         tg->SetMarkerColor(kGreen);
         tg->SetLineColor(kGreen);
         tg_syst->SetFillColorAlpha(kGreen, 0.5);
      }

      if (xaxis=="cent") {
         if (thebin.centbin().low()<=0 && thebin.centbin().high()<=0) padr->cd();
         else padl->cd();
      }
      tg_syst->Draw("2");      
      tg->Draw("P");      

      TString raplabel = Form("%.1f < |y| < %.1f, ",it->first.rapbin().low(),it->first.rapbin().high());
      TString otherlabel = "BWAA";
      if (xaxis == "pt") otherlabel.Form("%i%s-%i%s",(int) (it->first.centbin().low()/2.), "%", (int) (it->first.centbin().high()/2.), "%");
      if (xaxis == "cent") otherlabel.Form("%.1f < p_{T} < %.1f GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
      if (!isMB) {
         tleg->AddEntry(tg, (raplabel + otherlabel), "p");
      }

      // print tex
      ostringstream oss;
      oss.precision(1); oss.setf(ios::fixed);
      oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
      if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
      if (xaxis == "cent") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";

      addline(texname,oss.str());
      printGraph(tg, tg_syst, texname);

      // for the centrality dependence: we want Npart plotted, not the centrality
      if (xaxis == "cent") {
         centrality2npart(tg, isMB, false);
         centrality2npart(tg_syst, isMB, true, (150./1.6)*it->first.rapbin().low());
      }

      // in the case where the centrality dependence is plotted: treat the PP uncertainties as global systematics
      if (xaxis == "cent" && !isMB) {
         double x, dx, y, dy;
         dx = 10;
         x = 2*dx*(thebin.rapbin().low()/1.6) + dx;
         y = 1;
         anabin thebin(it->first.rapbin().low(),
               it->first.rapbin().high(),
               it->first.ptbin().low(),
               it->first.ptbin().high(),
               0,200);
         dy = gsyst[thebin].value;
         cout << "global syst: " << dy << endl;
         TBox *tbox = new TBox(x-dx,y-dy,x+dx,y+dy);
         if (thebin.rapbin() == binD((double) 0,1.6)) tbox->SetFillColorAlpha(kRed, 0.5);
         else if (thebin.rapbin() == binD(1.6,2.4)) tbox->SetFillColorAlpha(kBlue, 0.5);
         else tbox->SetFillColorAlpha(kGreen, 0.5);
         tbox->Draw();
      }

      cnt++;
      it_syst++;
   }

   if (xaxis=="cent") padl->cd();
   tleg->Draw();

   int iPos = 33;
   CMS_lumi( (TPad*) gPad, 105, iPos, "PbPb-PP" );

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE); 
   c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s.root",outputDir.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s.png",outputDir.c_str(), xaxis.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s.pdf",outputDir.c_str(), xaxis.c_str()));

   delete c1;

   // close tex
   closetex(texname);
   cout << "Closed " << texname << endl;
}

void centrality2npart(TGraphAsymmErrors* tg, bool issyst, bool isMB, double xshift) {
   int n = tg->GetN();
   for (int i=0; i<n; i++) {
      double x, y, exl, exh, eyl, eyh;
      x = tg->GetX()[i];
      y = tg->GetY()[i];
      exl = tg->GetErrorXlow(i);
      exh = tg->GetErrorXhigh(i);
      eyl = tg->GetErrorYlow(i);
      eyh = tg->GetErrorYhigh(i);
      x = isMB ? 150 + xshift : HI::findNpartAverage(x-exl,x+exh);
      if (!issyst) {
         exl = 0.;
         exh = 0.;
      } else {
         exl = !isMB ? 5 : 5./(1.-xfrac);
         exh = exl;
      }
      tg->SetPoint(i,x,y);
      tg->SetPointError(i,exl,exh,eyl,eyh);
   }
}
