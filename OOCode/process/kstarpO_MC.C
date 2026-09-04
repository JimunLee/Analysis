#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TF1Convolution.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TLine.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace TMath;

// Draw function

void Shreck();
TLatex* DrawTLatex(double x, double y, const TString& text, double textSize, int textFont, int textColor) {
  TLatex* txtInfo = new TLatex(x, y, text);
  txtInfo->SetNDC(kTRUE);
  txtInfo->SetTextSize(textSize);
  txtInfo->SetTextFont(textFont);
  txtInfo->SetTextColor(textColor);
  txtInfo->Draw("SAME");
  return txtInfo;
}


void HistoAxisTitles(TH1* hist, const TString& xTitle, const TString& yTitle,
		     double xTitleSize=0.05, int xTitleFont=62, double xTitleOffset=1.0,
		     double yTitleSize=0.05, int yTitleFont=62, double yTitleOffset=1.0)
{
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetXaxis()->SetTitleFont(xTitleFont);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);

  hist->GetYaxis()->SetTitle(yTitle);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetYaxis()->SetTitleFont(yTitleFont);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
}

// Main Drawfunction

void DrawBW()
{
  cout.precision(6);
  gStyle->SetOptStat(0);
  
  // pTBins setting 
  const Int_t nPtBins = 3;
  Double_t xBins[nPtBins+1] = {0.0, 2.5, 5.0, 10.0};

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  
  TString localpath = "MC/";

  // Data file inputs 
  TString mainfile = "../kstar-in-po/Result_rootfile/pOMC_AnalysisResults.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString USSName("kstar-in-o-o/hMC_USS;1");
  TString MixedUSSName("kstar-in-o-o/hMC_USS_Mix;1");
  TString Eventdata("kstar-in-o-o/nEvents_MC;1");
  
  // eventdata counts
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(1); //number 1 is clear number

  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  THnSparseT<TArrayF> *hUSS_MC = (THnSparseT<TArrayF>*)data->Get(USSName);
  THnSparseT<TArrayF> *hUSS_Mix_MC = (THnSparseT<TArrayF>*)data->Get(MixedUSSName);

  TH1D* hUSS [nPtBins];
  TH1D* hUSS_Mixed [nPtBins];
  TH1D* hUSS_Mixed_SUB [nPtBins];
  TH1D* hUSSMC_SUB [nPtBins];

  TH1* hUSSFit [nPtBins];
  TH1* hUSSBg [nPtBins];
  TH1* hUSSkstarSig [nPtBins];
  
  Double_t ptlow = 0;
  Double_t pthigh = 0;
  Double_t ReBinNum = 4;

  TH1D *hProj_USS;
  TH1D *hProj_USS_Mixed;
  for (Int_t i=0; i<nPtBins; i++)
    {
      Double_t nBins = 300.0;
      hUSSMC_SUB[i] = new TH1D(Form("hUSSMC_SUB_%i", i), Form("hUSSMC_SUB_%i", i), nBins, imlow, imtop);
      hUSS_Mixed_SUB[i] = new TH1D(Form("hUSS_Mixed_SUB_%i", i), Form("hUSS_Mixed_SUB_%i", i), nBins, imlow, imtop);

      hUSSBg[i] = new TH1D(Form("hUSSBg_%i", i), Form("hUSSBg_%i", i), nBins, imlow, imtop);
      hUSSkstarSig[i] = new TH1D(Form("hUSSkstarSig_%i", i), Form("hUSSkstarSig_%i", i), nBins, imlow, imtop);
      hUSSFit[i] = new TH1D(Form("hUSSFit_%i", i), Form("hUSSFit_%i", i), nBins, imlow, imtop);
      
      ptlow = xBins[i];
      pthigh = xBins[i+1];
      (hUSS_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hUSS_Mix_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);

      hUSS[i] = (TH1D*)hUSS_MC->Projection(2);
      hUSS_Mixed[i] = (TH1D*)hUSS_Mix_MC->Projection(2);
      hUSS[i]->SetName(Form("hUSSMC_%i", i));
      hUSS_Mixed[i]->SetName(Form("hUSSMC_Mixed_%i", i));
    }//for1

  for (Int_t j=0; j<nPtBins; j++)
    {
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      Double_t Nlow = 1.16;
      Double_t Nhigh = 1.19;

      Double_t Scale = hUSS_Mixed[j]->Integral(hUSS_Mixed[j]->FindBin(Nlow), hUSS_Mixed[j]->FindBin(Nhigh))
	/ hUSS[j]->Integral(hUSS[j]->FindBin(Nlow), hUSS[j]->FindBin(Nhigh));
      
      double delta =1.0;
      hUSS_Mixed[j]->Scale(1.0/Scale);

      HistoAxisTitles(hUSS[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS[j]->SetTitle(Form(""));

      hUSS[j]->SetMarkerStyle(4); //circle
      hUSS[j]->SetMarkerSize(1.2);
      hUSS[j]->SetMarkerColor(4);
      hUSS[j]->SetLineColor(4);

      hUSS_Mixed[j]->SetMarkerStyle(26); //triangle
      hUSS_Mixed[j]->SetMarkerSize(1.2);
      hUSS_Mixed[j]->SetMarkerColor(2);
      hUSS_Mixed[j]->SetLineColor(2);

      hUSS[j]->Rebin(5);
      hUSS_Mixed[j]->Rebin(5);
      hUSS[j]->Scale(0.2);
      hUSS_Mixed[j]->Scale(0.2);

      // 1D projection plots
      TCanvas* plots = new TCanvas(Form("plots_%d", j), "", 940, 800);
      plots->SetLeftMargin(0.16);
      plots->SetRightMargin(0.06);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);

      HistoAxisTitles(hUSS[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS[j]->SetTitle(Form(""));

      hUSS[j]->Draw();
      hUSS_Mixed[j]->Draw("SAME");

      TLatex* txtInfoA = DrawTLatex(0.641, 0.310, "O-O yield", 0.04, 62, 1);
      TLatex* txtInfoB = DrawTLatex(0.172, 0.935, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoC = DrawTLatex(0.641, 0.265, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
      TLatex* txtInfoD = DrawTLatex(0.641, 0.221, "Min. Bias", 0.04, 42, 1);
      TLatex* txtInfoF = DrawTLatex(0.641, 0.176, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);

      TLegend* legend = new TLegend(0.621,0.381,0.855,0.498); //LEFT x, y, Right x, y
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.04);
      legend->AddEntry(hUSS[j], "K^{+}#pi^{-}", "lpf");
      legend->AddEntry(hUSS_Mixed[j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
      legend->Draw("SAME");
      
      plots->SaveAs(Form("../kstar-in-po/Plot/"+localpath+"/MC_projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");
    
      // ==========================
      // Drawing Substraction
      
      hUSSMC_SUB[j] = (TH1D*)hUSS[j]->Clone(Form("hUSS_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hUSS_Mixed_SUB[j] = (TH1D*)hUSS_Mixed[j]->Clone(Form("hUSS_Mixed_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      TCanvas* substraction = new TCanvas(Form("MC_substraction_%d", j), "", 550, 440);
      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.08);
      substraction->SetBottomMargin(0.14);

      hUSSMC_SUB[j]->Add(hUSS_Mixed_SUB[j], -1);
      /* hUSSMC_SUB[j]->Rebin(5); */
      /* hUSSMC_SUB[j]->Scale(1.0/5.0); */
      
      /* hUSSMC_SUB[j]->Draw(); */
      /* hUSSMCTrue[j]->Rebin(5); */
      /* hUSSMCTrue[j]->Scale(1.0/5.0); */
      HistoAxisTitles(hUSSMC_SUB[j], "#it{M}_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]", "Counts");
      hUSSMC_SUB[j]->SetTitle(Form(""));

      TLatex* txtInfoa = DrawTLatex(0.55, 0.88, "O-O yield", 0.04, 62, 1);
      TLatex* txtInfob = DrawTLatex(0.55, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoc = DrawTLatex(0.55, 0.83, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
      TLatex* txtInfoe = DrawTLatex(0.55, 0.741, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);

      substraction->SaveAs(Form("../kstar-in-po/Plot/"+localpath+"/MCRec_reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh)); 
    }//for2
  
  TF1* fkstar [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  TF1* fBg_POST [nPtBins];
  TF1* fkstar_POST [nPtBins];

  double kstarMass = 0.89581;
  double kstarwidth = 0.0474;
  
  double kstarwidthMin = 0.0434;
  double kstarwidthMax = 0.046;
  
  double BgMinRange = 0.60;
  double kstarMinRange = 0.76;
  double kstarMaxRange = 1.1;
  double BgMaxRange = 1.20;

  for (int k=0; k<nPtBins; k++)
    {
      hUSSkstarSig[k]->Add(hUSSMC_SUB[k]);
      hUSSBg[k]->Add(hUSSMC_SUB[k]);
      hUSSFit[k]->Add(hUSSMC_SUB[k]);

      ptlow = xBins[k];
      pthigh = xBins[k+1];
      // ==============
      // Function
      // ==============
      /* fBg[k] = new TF1(Form("kstar_Bg_%i",k), "pol3(0)", BgMinRange, BgMaxRange); */
      /* fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange); */
      /* fkstar[k]->SetParameter(1, kstarMass); */
      /* fkstar[k]->SetParameter(2, kstarwidth); */
      
      fFit[k] = new TF1(Form("Fullfit_%i", k), "pol3(0) + [4]*BreitWignerRelativistic(x, [5], [6])", BgMinRange, BgMaxRange);
      fFit[k]->SetParameter(5, kstarMass);
      fFit[k]->SetParameter(6, kstarwidth);

      if(k==0){
	/* for (int bin = hUSSBg[k]->FindBin(0.82); bin <= hUSSBg[k]->FindBin(0.96); bin++) */
	/*   hUSSBg[k]->SetBinError(bin, 0); */

        /* hUSSBg[k]->Fit(fBg[k], "SR0", "", BgMinRange, BgMaxRange); */
        /* hUSSkstarSig[k]->Fit(fkstar[k], "SR0", "", BgMinRange, BgMaxRange); */
	    
        fFit[k]->SetParameter(0, fFit[k]->GetParameter(0));
        fFit[k]->SetParameter(1, fFit[k]->GetParameter(1));
        fFit[k]->SetParameter(2, fFit[k]->GetParameter(2));
        fFit[k]->SetParameter(3, fFit[k]->GetParameter(3));
        fFit[k]->SetParameter(4, fFit[k]->GetParameter(4));
        fFit[k]->FixParameter(5, fFit[k]->GetParameter(5));
        fFit[k]->FixParameter(6, fFit[k]->GetParameter(6));
      } else{
        fFit[k]->SetParameter(0, fBg_POST[k-1]->GetParameter(0));
        fFit[k]->SetParameter(1, fBg_POST[k-1]->GetParameter(1));
        fFit[k]->SetParameter(2, fBg_POST[k-1]->GetParameter(2));
        fFit[k]->SetParameter(3, fBg_POST[k-1]->GetParameter(3));
	fFit[k]->SetParameter(4, fkstar_POST[k-1]->GetParameter(0));
        fFit[k]->FixParameter(5, fkstar_POST[k-1]->GetParameter(1));
        fFit[k]->FixParameter(6, fkstar_POST[k-1]->GetParameter(2));
      }
      hUSSFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange);
	
      /* if(k==0){ */
      /* 	for (int bin = hUSSBg[k]->FindBin(0.82); bin <= hUSSBg[k]->FindBin(0.96); bin++) */
      /* 	  hUSSBg[k]->SetBinError(bin, 0); */

      /*   hUSSBg[k]->Fit(fBg[k], "SR0", "", BgMinRange, BgMaxRange); */
      /*   /\* hUSSkstarSig[k]->Fit(fkstar[k], "SR0", "", BgMinRange, BgMaxRange); *\/ */
	    
      /*   fFit[k]->SetParameter(0, fBg[k]->GetParameter(0)); */
      /*   fFit[k]->SetParameter(1, fBg[k]->GetParameter(1)); */
      /*   fFit[k]->SetParameter(2, fBg[k]->GetParameter(2)); */
      /*   fFit[k]->SetParameter(3, fBg[k]->GetParameter(3)); */
      /*   fFit[k]->SetParameter(4, fkstar[k]->GetParameter(0)); */
      /*   fFit[k]->FixParameter(5, fkstar[k]->GetParameter(1)); */
      /*   fFit[k]->FixParameter(6, fkstar[k]->GetParameter(2)); */
      /* } else{ */
      /*   fFit[k]->SetParameter(0, fBg_POST[k-1]->GetParameter(0)); */
      /*   fFit[k]->SetParameter(1, fBg_POST[k-1]->GetParameter(1)); */
      /*   fFit[k]->SetParameter(2, fBg_POST[k-1]->GetParameter(2)); */
      /*   fFit[k]->SetParameter(3, fBg_POST[k-1]->GetParameter(3)); */
      /* 	fFit[k]->SetParameter(4, fkstar[k-1]->GetParameter(0)); */
      /*   fFit[k]->FixParameter(5, fkstar[k]->GetParameter(1)); */
      /*   fFit[k]->FixParameter(6, fkstar[k]->GetParameter(2)); */
      /* } */
      /* hUSSFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange); */
      //hUSSFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange);
      
      // ================================
      // for fitting each one
      // ================================
      fBg_POST[k] = new TF1(Form("kstar_Bg_POST_%i",k), "pol3(0)", BgMinRange, BgMaxRange);
      fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(0));
      fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(1));
      fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(2));
      fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(3));

      fkstar_POST[k] = new TF1(Form("kstar_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
      fkstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(4));
      fkstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(5));
      fkstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(6));

      TFitResultPtr ptr = hUSSFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange);
      double Chi2 = ptr->Chi2();
      int NDF = ptr->Ndf();
      double chi2Ndf = Chi2/NDF;
      string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);

      // Draw
      TCanvas* lines = new TCanvas(Form("line_%d",k), "", 940, 800);
      lines->SetLeftMargin(0.14);
      lines->SetRightMargin(0.1);
      lines->SetTopMargin(0.08);
      lines->SetBottomMargin(0.14);
      HistoAxisTitles(hUSSFit[k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
		      "d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);
	
      hUSSFit[k]->SetTitle(Form(""));
      hUSSFit[k]->SetLineColor(1);
      hUSSFit[k]->GetYaxis()->SetRangeUser(0.0001, hUSSFit[k]->GetMaximum()*1.1);
      hUSSFit[k]->Draw();

      // total
      fFit[k]->SetLineColor(2);
      fFit[k]->SetLineWidth(4);
      fFit[k]->Draw("SAME");

      // background
      fBg_POST[k]->SetLineColor(4);
      fBg_POST[k]->SetLineWidth(4);
      fBg_POST[k]->Draw("SAME");

      // kstar
      fkstar_POST[k]->SetLineColor(6);
      fkstar_POST[k]->SetLineWidth(4);
      fkstar_POST[k]->Draw("SAME");
      
      TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, "p-O yield", 0.04, 62, 1);
      TLatex* txtInfoCC = DrawTLatex(0.56, 0.833, "p-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
      TLatex* txtInfoDD = DrawTLatex(0.56, 0.789, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
      TLatex* txtInfoEE = DrawTLatex(0.56, 0.74, chi2NdfText, 0.04, 42, 1);

      TLegend *legendd = new TLegend(0.55,0.54,0.66,0.64);
      legendd->SetFillColor(4000);
      legendd->SetBorderSize(0);
      legendd->SetTextSize(0.03);
      legendd->AddEntry(hUSSFit[k],"MC", "lep");
      legendd->AddEntry(fBg_POST[k], "Bkg. MCFit", "l");
      legendd->AddEntry(fkstar_POST[k], "BreitWignerRelativistic. MCFit", "l");
      legendd->AddEntry(fFit[k], "Total. Fit", "l");
      legendd->Draw("SAME");
	
      lines->SaveAs(Form("../kstar-in-po/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
      // lines->SaveAs(Form("../kstar-in-po/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.pdf", ptlow, pthigh));
    } //for
  
}//DrawBW


void Shreck() {
  cout << R"(⢀⡴⠑⡄⠀⠀⠀⠀⠀⠀⠀⣀⣀⣤⣤⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
    ⠸⡇⠀⠿⡀⠀⠀⠀⣀⡴⢿⣿⣿⣿⣿⣿⣿⣿⣷⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠑⢄⣠⠾⠁⣀⣄⡈⠙⣿⣿⣿⣿⣿⣿⣿⣿⣆⠀⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⢀⡀⠁⠀⠀⠈⠙⠛⠂⠈⣿⣿⣿⣿⣿⠿⡿⢿⣆⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⢀⡾⣁⣀⠀⠴⠂⠙⣗⡀⠀⢻⣿⣿⠭⢤⣴⣦⣤⣹⠀⠀⠀⢀⢴⣶⣆ 
    ⠀⠀⢀⣾⣿⣿⣿⣷⣮⣽⣾⣿⣥⣴⣿⣿⡿⢂⠔⢚⡿⢿⣿⣦⣴⣾⠁⠸⣼⡿ 
    ⠀⢀⡞⠁⠙⠻⠿⠟⠉⠀⠛⢹⣿⣿⣿⣿⣿⣌⢤⣼⣿⣾⣿⡟⠉⠀⠀⠀⠀⠀ 
    ⠀⣾⣷⣶⠇⠀⠀⣤⣄⣀⡀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀ 
    ⠀⠉⠈⠉⠀⠀⢦⡈⢻⣿⣿⣿⣶⣶⣶⣶⣤⣽⡹⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⠀⠉⠲⣽⡻⢿⣿⣿⣿⣿⣿⣿⣷⣜⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣷⣶⣮⣭⣽⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⣀⣀⣈⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⠀⠹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠟⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠻⠿⠿⠿⠿⠛⠉
    )" << endl;

    }
