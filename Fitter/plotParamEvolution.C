#define plotParamEvolution_cxx
#include "plotParamEvolution.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>

void plotParamEvolution::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L plotParamEvolution.C
//      Root > plotParamEvolution t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

  //___________ Create Histos
  int parN = 5;
  const char* parName[5] = {"alpha", "n", "sigma1", "rSigma21", "f"};
  int varN = 4;
  const char* varname[4] = {"IntPt","IntCent","Pt", "Cent"};
  int rapN = 2;
  const char* rapname[2] = {"MidY", "FwdY"};
  int sysN = 2;
  const char* sysname[2] = {"PP", "PbPb"};
  TObjArray* histosArr = new TObjArray();
  histosArr->SetOwner(kTRUE);
  for ( int i = 0 ; i < parN ; i++ )
  {
    for ( int j = 0 ; j < varN ; j++ )
    {
      for ( int k = 0 ; k < sysN ; k++ )
      {
        for ( int l = 0 ; l < rapN ; l++ )
        {
          int nbins = 0;
          double xmin = 0.;
          double xmax = 0.;
          TString xName("");
          TString yName("");
          if ( !strcmp(varname[j],"IntPt") || !strcmp(varname[j],"Pt") )
          {
            nbins = 60;
            xmin = 0.;
            xmax = 30.;
            xName += "p_{T}";
          }
          else
          {
            nbins = 200;
            xmin = 0.;
            xmax = 200.;
            xName += "cent";
          }
          TString hName = Form("h%s_%s_%s_%s",parName[i],rapname[l],varname[j],sysname[k]);
          TH1* h = new TH1D(hName.Data(),hName.Data(), nbins, xmin, xmax);
          h->GetXaxis()->SetTitle(xName.Data());
          h->GetYaxis()->SetTitle(parName[i]);
          h->SetStats(0);
          h->GetXaxis()->SetLabelSize(0.06);
          h->GetXaxis()->SetTitleSize(0.06);
          h->GetXaxis()->SetTitleOffset(0.76);
          h->GetYaxis()->SetLabelSize(0.06);
          h->GetYaxis()->SetTitleSize(0.06);
          h->SetMarkerStyle(20);
          h->SetMarkerColor(1);
          h->SetLineColor(1);
          histosArr->Add(h);
        }
      }
    }
  }
  //___________
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  //___________ Fill histos
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    if ( strcmp(jpsiName,"DoubleCrystalBall") ) continue;
    
    TH1* hf(0x0);
    TH1* hs(0x0);
    TH1* hr(0x0);
    TH1* ha(0x0);
    TH1* hn(0x0);
    
    if ( (ymin ==  0.) && (ymax > 1.59 && ymax < 1.61)) // Mid-Y
    {
      if ( (centmin == 0.) && (centmax > 198 && centmax < 201) ) // Cent integrated
      {
        if ( (ptmin == 6.5) && (ptmax > 29.7 && ptmax < 30.1) ) // Cent-Pt Integrated
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_MidY_IntPt_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_MidY_IntPt_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_MidY_IntPt_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_MidY_IntPt_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_MidY_IntPt_%s",collSystem)));
        }
        else // Pt differential
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_MidY_Pt_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_MidY_Pt_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_MidY_Pt_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_MidY_Pt_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_MidY_Pt_%s",collSystem)));
        }
        
        if ( !hf || !hs || !hr || !ha || !hn )
        {
          std::cout << "Error: Some histogram is missing" << std::endl;
          std::cout << collSystem << " ; " << ptmin << " ; " << ptmax << " ; " << centmin << " ; " << centmax << " ; " << ymin << " ; " << ymax << " ; " << jpsiName << std::endl;
          continue;
        }
        
        Int_t bin = hf->FindBin(ptmin);
        hf->SetBinContent(bin,f_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hf->SetBinError(bin,f_Jpsi_err);
        hs->SetBinContent(bin,sigma1_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hs->SetBinError(bin,sigma1_Jpsi_err);
        hr->SetBinContent(bin,rSigma21_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hr->SetBinError(bin,rSigma21_Jpsi_err);
        ha->SetBinContent(bin,alpha_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) ha->SetBinError(bin,alpha_Jpsi_err);
        hn->SetBinContent(bin,n_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hn->SetBinError(bin,n_Jpsi_err);
      }
      if ( (ptmin == 6.5) && (ptmax > 29.7 && ptmax < 30.1) ) // Pt integrated
      {
        if ( (centmin == 0.) && (centmax > 198 && centmax < 201.) ) // Cent differential
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_MidY_IntCent_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_MidY_IntCent_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_MidY_IntCent_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_MidY_IntCent_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_MidY_IntCent_%s",collSystem)));
        }
        else
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_MidY_Cent_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_MidY_Cent_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_MidY_Cent_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_MidY_Cent_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_MidY_Cent_%s",collSystem)));
        }
        
        if ( !hf || !hs || !hr || !ha || !hn )
        {
          std::cout << "Error: Some histogram is missing" << std::endl;
          std::cout << collSystem << " ; " << ptmin << " ; " << ptmax << " ; " << centmin << " ; " << centmax << " ; " << ymin << " ; " << ymax << " ; " << jpsiName << std::endl;
          continue;
        }
        
        Int_t bin = hf->FindBin(centmin);
        hf->SetBinContent(bin,f_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hf->SetBinError(bin,f_Jpsi_err);
        hs->SetBinContent(bin,sigma1_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hs->SetBinError(bin,sigma1_Jpsi_err);
        hr->SetBinContent(bin,rSigma21_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hr->SetBinError(bin,rSigma21_Jpsi_err);
        ha->SetBinContent(bin,alpha_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) ha->SetBinError(bin,alpha_Jpsi_err);
        hn->SetBinContent(bin,n_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hn->SetBinError(bin,n_Jpsi_err);
      }
    }
    else if ( (ymin > 1.59 && ymin < 1.61) && (ymax > 2.39 && ymax < 2.405) ) // Fwd-Y
    {
      if ( (centmin == 0.) && (centmax > 198 && centmax < 201) ) // Cent integrated
      {
        if ( (ptmin == 3.0) && (ptmax > 29.7 && ptmax < 30.1) ) // Cent-Pt Integrated
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_FwdY_IntPt_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_FwdY_IntPt_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_FwdY_IntPt_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_FwdY_IntPt_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_FwdY_IntPt_%s",collSystem)));
        }
        else // Pt differential
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_FwdY_Pt_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_FwdY_Pt_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_FwdY_Pt_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_FwdY_Pt_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_FwdY_Pt_%s",collSystem)));
        }
        
        if ( !hf || !hs || !hr || !ha || !hn )
        {
          std::cout << "Error: Some histogram is missing" << std::endl;
          std::cout << collSystem << " ; " << ptmin << " ; " << ptmax << " ; " << centmin << " ; " << centmax << " ; " << ymin << " ; " << ymax << " ; " << jpsiName << std::endl;
          continue;
        }
        
        Int_t bin = hf->FindBin(ptmin);
        hf->SetBinContent(bin,f_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hf->SetBinError(bin,f_Jpsi_err);
        hs->SetBinContent(bin,sigma1_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hs->SetBinError(bin,sigma1_Jpsi_err);
        hr->SetBinContent(bin,rSigma21_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hr->SetBinError(bin,rSigma21_Jpsi_err);
        ha->SetBinContent(bin,alpha_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) ha->SetBinError(bin,alpha_Jpsi_err);
        hn->SetBinContent(bin,n_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hn->SetBinError(bin,n_Jpsi_err);
      }
      
      if ( (ptmin == 3.0) && (ptmax > 29.7 && ptmax < 30.1) ) // Pt integrated
      {
        if ( (centmin == 0.) && (centmax > 198 && centmax < 201) ) // Cent differential
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_FwdY_IntCent_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_FwdY_IntCent_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_FwdY_IntCent_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_FwdY_IntCent_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_FwdY_IntCent_%s",collSystem)));
        }
        else
        {
          hf = static_cast<TH1*>(histosArr->FindObject(Form("hf_FwdY_Cent_%s",collSystem)));
          hs = static_cast<TH1*>(histosArr->FindObject(Form("hsigma1_FwdY_Cent_%s",collSystem)));
          hr = static_cast<TH1*>(histosArr->FindObject(Form("hrSigma21_FwdY_Cent_%s",collSystem)));
          ha = static_cast<TH1*>(histosArr->FindObject(Form("halpha_FwdY_Cent_%s",collSystem)));
          hn = static_cast<TH1*>(histosArr->FindObject(Form("hn_FwdY_Cent_%s",collSystem)));
        }
        
        if ( !hf || !hs || !hr || !ha || !hn )
        {
          std::cout << "Error: Some histogram is missing" << std::endl;
          std::cout << collSystem << " ; " << ptmin << " ; " << ptmax << " ; " << centmin << " ; " << centmax << " ; " << ymin << " ; " << ymax << " ; " << jpsiName << std::endl;
          continue;
        }
        
        Int_t bin = hf->FindBin(centmin);
        hf->SetBinContent(bin,f_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hf->SetBinError(bin,f_Jpsi_err);
        hs->SetBinContent(bin,sigma1_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hs->SetBinError(bin,sigma1_Jpsi_err);
        hr->SetBinContent(bin,rSigma21_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hr->SetBinError(bin,rSigma21_Jpsi_err);
        ha->SetBinContent(bin,alpha_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) ha->SetBinError(bin,alpha_Jpsi_err);
        hn->SetBinContent(bin,n_Jpsi_val);
        if( !strcmp(collSystem,"PP") ) hn->SetBinError(bin,n_Jpsi_err);
      }
    }
  }
  //___________
  
  
  //___________ Plot results
  TObjArray* aResults = new TObjArray();
  aResults->SetOwner(kTRUE);
  
  int varN2 = 2;
  const char* varname2[2] = {"Pt", "Cent"};

  for ( int i = 0 ; i < parN ; i++ )
  {
    for ( int j = 0 ; j < varN2 ; j++ )
    {
      for ( int k = 0 ; k < sysN ; k++ )
      {
        for ( int l = 0 ; l < rapN ; l++ )
        {
          TString hName = Form("h%s_%s_%s_%s",parName[i],rapname[l],varname2[j],sysname[k]);
          TString hIntName = Form("h%s_%s_Int%s_%s",parName[i],rapname[l],varname2[j],sysname[k]);
          
          TCanvas* c1 = new TCanvas(hName.Data(),"canvas",81,98,700,504);
          c1->Range(-5.547577,-0.004024593,31.93896,0.03073326);
          c1->SetFillColor(0);
          c1->SetBorderMode(0);
          c1->SetBorderSize(2);
          c1->SetLeftMargin(0.1479885);
          c1->SetRightMargin(0.05172414);
          c1->SetTopMargin(0.08421053);
          c1->SetBottomMargin(0.1157895);
          c1->SetFrameBorderMode(0);
          c1->SetFrameBorderMode(0);
        
          TH1* haInt = static_cast<TH1*>(histosArr->FindObject(hIntName.Data())->Clone());
          haInt->SetMarkerStyle(21);
          haInt->SetMarkerColor(2);
          haInt->SetLineColor(2);
          TH1* ha = static_cast<TH1*>(histosArr->FindObject(hName.Data())->Clone());
          ha->GetYaxis()->SetTitleOffset(1.15);
          ha->GetXaxis()->SetTitleOffset(0.85);
          
          //  TLegend* l1 = new TLegend();
          //  l1->AddEntry(haInt,"Integrated","lp");
          //  l1->AddEntry(ha,"p_{T} diferential","lp");
          
          TLatex *  text = new TLatex(0.4913793,0.2494759,Form("%s %s",sysname[k], strcmp(rapname[l],"MidY") ? "1.6 < |y| < 2.4" : "0 < |y| < 1.6"));
          text->SetNDC();
          text->SetTextFont(42);
          text->SetTextSize(0.06708595);
          text->SetLineWidth(2);
          
          ha->Draw("p");
          haInt->Draw("psame");
          //  l1->Draw("same");
          text->Draw("same");
          c1->Update();
          aResults->Add(c1);
          
          c1->SaveAs(Form("%s.pdf",hName.Data()));
          
        }
      }
    }
  }
  
  
  TFile* f = new TFile("fitParamEvol.root","RECREATE");
  histosArr->Write("paramEvolHistos", TObject::kOverwrite | TObject::kSingleKey);
  aResults->Write("plots", TObject::kOverwrite | TObject::kSingleKey);
  f->Close();
  
}
