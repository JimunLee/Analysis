#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

void DrawBW(TString Kind, TString USS, TString LSS)
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  // pTBins setting 
  const Int_t nPtBins = 9;
  Double_t xBins[nPtBins+1];

  // TLatex setting
  TLatex* SpecialLatex_plot = nullptr;
  TLatex* SpecialLatex_fit = nullptr;
  TLatex* SpecialLatex_plot2 = nullptr;
  TLatex* SpecialLatex_fit2 = nullptr;

  if (Kind == "Inclusive"){
    Double_t inclusive_xBins[nPtBins+1] = {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(inclusive_xBins), end(inclusive_xBins), xBins);
    
    SpecialLatex_plot = DrawTLatex(0.657, 0.775, "Min. Bias", 0.04, 42, 1);
    SpecialLatex_fit = DrawTLatex(0.557, 0.79, "Min. Bias", 0.04, 42, 1);
  }
  else if (Kind == "Outside"){
    Double_t outside_xBins[nPtBins+1] = {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(outside_xBins), end(outside_xBins), xBins);
    
    SpecialLatex_plot = DrawTLatex(0.547, 0.285, "Min. Bias", 0.04, 42, 1);
    SpecialLatex_fit = DrawTLatex(0.557, 0.79, "Min. Bias", 0.04, 42, 1);
  }
  else if (Kind == "Inside"){
    Double_t inside_xBins[nPtBins+1] = {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(inside_xBins), end(inside_xBins), xBins);

    SpecialLatex_plot = DrawTLatex(0.657, 0.775, "#Delta #it{R}_{ jet}^{ K*} < 0.4", 0.04, 42, 1);
    SpecialLatex_plot2 = DrawTLatex(0.657, 0.665, "#it{p}_{T}^{ Jet ch uncorr} > 8 GeV/c", 0.04, 42, 1);
    SpecialLatex_fit = DrawTLatex(0.558, 0.785, "#Delta #it{R}_{ jet}^{ K*} < 0.4", 0.04, 42, 1);
    SpecialLatex_fit2 = DrawTLatex(0.558, 0.630, "#it{p}_{T}^{ Jet ch uncorr} > 8 GeV/c", 0.04, 42, 1);
  }
  else{
    cout << "Please specify either 'Inside', 'Outside' or 'Inclusive'." << endl;
    return;
  }

  // invariant mass range
  Double_t imlow = 0.75;
  Double_t imtop = 1.25;
  
  TString widelocal = "kstarInjets/";
  TString localpath = widelocal + Kind;

  // Data file inputs
  TString mainfile = "../kstar-in-jets/Result_rootfile/tpctof_two.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString SignalName("phi-in-jets/"+USS+";1");
  TString BackgroundName("phi-in-jets/"+LSS+";1");
  TString Eventdata("phi-in-jets/nEvents;1");

  // eventdata counts
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(2); //number 2 is clear number

  // Drawing JetPtHistogram
  TH1D* FJpt = (TH1D*)data->Get("phi-in-jets/FJptHistogram;1");
  TCanvas* FJc = new TCanvas("FJc", "fastjetpt", 940, 880);
  FJc->SetLeftMargin(0.15);
  FJc->SetRightMargin(0.05);
  FJc->SetTopMargin(0.10);
  FJc->SetBottomMargin(0.15);

  FJpt->SetMarkerStyle(20);
  FJpt->SetMarkerSize(1.4);
  FJpt->Sumw2();
  FJpt->Draw();

  HistoAxisTitles(FJpt, "#it{p}_{T} [GeV/#it{c}]", "#it{N }_{jet}^{Raw}");
  FJpt->SetTitle("");

  TLatex* textpTa = DrawTLatex(0.526, 0.915, "ALICE #bf{WIP}", 0.05, 62, 1);
  TLatex* textpTb = DrawTLatex(0.526, 0.854, "pp #sqrt{s} = 13.6 TeV", 0.05, 42, 1);
  TLatex* textpTc = DrawTLatex(0.526, 0.790, "#it{p}_{T}^{ jet ch uncorr} > 8 GeV/#it{c}", 0.05, 42, 1);
  TLatex* textpTd = DrawTLatex(0.526, 0.725, "#it{R} = 0.4", 0.05, 42, 1);

  FJc->SetLogy();
  FJc->SaveAs("../kstar-in-jets/Plot/"+widelocal+"FastJetpt.png", "RECREATE");
  
  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  TH3D* hLSS = (TH3D*)data->Get(BackgroundName);
  TH3D* hUSS = (TH3D*)data->Get(SignalName);
  
  TH1D* hLSS_Y [nPtBins];
  TH1D* hUSS_Y [nPtBins];
  
  TH1D* hLSS_Y_SUB [nPtBins];
  TH1D* hUSS_Y_SUB [nPtBins];

  TH1* hUSS_Y_Bg [nPtBins];
  TH1* hUSS_Y_ExtraBg [nPtBins];
  TH1* hUSS_Y_kstarSig [nPtBins];

  TH1* hUSS_Y_BgTmp [nPtBins];
  TH1* hUSS_Y_ExtraBgTmp [nPtBins];
  TH1* hUSS_Y_kstarSigTmp [nPtBins];
  
  TH1* hUSS_Y_BgTmp2 [nPtBins];
  TH1* hUSS_Y_ExtraBgTmp2 [nPtBins];
  TH1* hUSS_Y_kstarSigTmp2 [nPtBins];
  
  TH1* hUSS_Y_Fit [nPtBins];
  TH1* hUSS_Y_FitTmp [nPtBins];
  TH1* hUSS_Y_FitTmp2 [nPtBins];

  if(Kind == "Inside"){
    for (Int_t i=0; i<nPtBins; i++)
      {
	Double_t nBins = 100.0;
	
	hLSS_Y[i] = new TH1D(Form("hLSS_Y_%i", i), Form("hLSS_Y_%i", i), nBins, imlow, imtop);
	hUSS_Y[i] = new TH1D(Form("hUSS_Y_%i", i), Form("hUSS_Y_%i", i), nBins, imlow, imtop);

	hUSS_Y_Bg[i] = new TH1D(Form("hUSS_Y_Bg_%i", i), Form("hUSS_Y_Bg_%i", i), nBins, imlow, imtop);
	hUSS_Y_ExtraBg[i] = new TH1D(Form("hUSS_Y_ExtraBg_%i", i), Form("hUSS_Y_ExtraBg_%i", i), nBins, imlow, imtop);
	hUSS_Y_kstarSig[i] = new TH1D(Form("hUSS_Y_kstarSig_%i", i), Form("hUSS_Y_kstarSig_%i", i), nBins, imlow, imtop);      
	hUSS_Y_Fit[i] = new TH1D(Form("hUSS_Y_Fit_%i", i), Form("hUSS_Y_Fit_%i", i), nBins, imlow, imtop);

	hUSS_Y_BgTmp[i] = new TH1D(Form("hUSS_Y_BgTmp_%i", i), Form("hUSS_Y_BgTmp_%i", i), nBins, imlow, imtop);
	hUSS_Y_ExtraBgTmp[i] = new TH1D(Form("hUSS_Y_ExtraBgTmp_%i", i), Form("hUSS_Y_ExtraBgTmp_%i", i), nBins, imlow, imtop);
	hUSS_Y_kstarSigTmp[i] = new TH1D(Form("hUSS_Y_kstarSigTmp_%i", i), Form("hUSS_Y_kstarSigTmp_%i", i), nBins, imlow, imtop);      
	hUSS_Y_FitTmp[i] = new TH1D(Form("hUSS_Y_FitTmp_%i", i), Form("hUSS_Y_FitTmp_%i", i), nBins, imlow, imtop);

	hUSS_Y_BgTmp2[i] = new TH1D(Form("hUSS_Y_BgTmp2_%i", i), Form("hUSS_Y_BgTmp2_%i", i), nBins, imlow, imtop);
	hUSS_Y_ExtraBgTmp2[i] = new TH1D(Form("hUSS_Y_ExtraBgTmp2_%i", i), Form("hUSS_Y_ExtraBgTmp2_%i", i), nBins, imlow, imtop);
	hUSS_Y_kstarSigTmp2[i] = new TH1D(Form("hUSS_Y_kstarSigTmp2_%i", i), Form("hUSS_Y_kstarSigTmp2_%i", i), nBins, imlow, imtop);      
	hUSS_Y_FitTmp2[i] = new TH1D(Form("hUSS_Y_FitTmp2_%i", i), Form("hUSS_Y_FitTmp2_%i", i), nBins, imlow, imtop);
      }
  }
  else{
    for (Int_t i=0; i<nPtBins; i++)
      {
	Double_t nBins = 500.0;
	
	hLSS_Y[i] = new TH1D(Form("hLSS_Y_%i", i), Form("hLSS_Y_%i", i), nBins, imlow, imtop);
	hUSS_Y[i] = new TH1D(Form("hUSS_Y_%i", i), Form("hUSS_Y_%i", i), nBins, imlow, imtop);

	hUSS_Y_Bg[i] = new TH1D(Form("hUSS_Y_Bg_%i", i), Form("hUSS_Y_Bg_%i", i), nBins, imlow, imtop);
	hUSS_Y_ExtraBg[i] = new TH1D(Form("hUSS_Y_ExtraBg_%i", i), Form("hUSS_Y_ExtraBg_%i", i), nBins, imlow, imtop);
	hUSS_Y_kstarSig[i] = new TH1D(Form("hUSS_Y_kstarSig_%i", i), Form("hUSS_Y_kstarSig_%i", i), nBins, imlow, imtop);      
	hUSS_Y_Fit[i] = new TH1D(Form("hUSS_Y_Fit_%i", i), Form("hUSS_Y_Fit_%i", i), nBins, imlow, imtop);
      }
  }

  Double_t ptlow = 0;
  Double_t pthigh = 0;
  Double_t ReBinNum = 4;
  
  for (Int_t j=0; j<nPtBins; j++)
    {
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      (hLSS->GetYaxis())->SetRangeUser(ptlow, pthigh);
      (hUSS->GetYaxis())->SetRangeUser(ptlow, pthigh);

      hLSS_Y[j]->Add(hLSS->Project3D("z o"), 1);
      hUSS_Y[j]->Add(hUSS->Project3D("z o"), 1);
      // =======================================================================================
      // Normalization 
      Double_t Nlow = 1.19;
      Double_t Nhigh = 1.24;

      /* hUSS_Y[j]->Scale(1.0/nEvents); */
      /* hLSS_Y[j]->Scale(1.0/nEvents); */

      Double_t Scale = hUSS_Y[j]->Integral(hUSS_Y[j]->FindBin(Nlow), hUSS_Y[j]->FindBin(Nhigh))
       	/ hLSS_Y[j]->Integral(hLSS_Y[j]->FindBin(Nlow), hLSS_Y[j]->FindBin(Nhigh));
      hUSS_Y[j]->Scale(Scale);
      hLSS_Y[j]->Scale(Scale);
      
      // ========================================================================================
      // 1D projection plots
      TCanvas* plots = new TCanvas(Form("plots_%d", j), "", 940, 800);
      plots->SetLeftMargin(0.14);
      plots->SetRightMargin(0.08);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);

      HistoAxisTitles(hUSS_Y[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "1/#it{N}_{ ev} d#it{N}^{ K#pi}/d#it{M}_{ inv}");
      hUSS_Y[j]->SetTitle(Form(""));

      hUSS_Y[j]->SetMarkerStyle(26); //triangle
      hUSS_Y[j]->SetMarkerSize(1.2);
      hUSS_Y[j]->SetMarkerColor(2);
      hUSS_Y[j]->SetLineColor(2);
      hUSS_Y[j]->Draw();

      hLSS_Y[j]->SetMarkerStyle(4); //circle
      hLSS_Y[j]->SetMarkerSize(1.2);
      hLSS_Y[j]->SetMarkerColor(4);
      hLSS_Y[j]->SetLineColor(4);
      hLSS_Y[j]->Draw("SAME");
      
      TLatex* txtInfoA = DrawTLatex(0.657, 0.876, Kind+" yield", 0.04, 62, 1);
      TLatex* txtInfoB = DrawTLatex(0.657, 0.929, "ALICE #bf{WIP}", 0.05, 62, 1);
      TLatex* txtInfoC = DrawTLatex(0.657, 0.825, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
      if (SpecialLatex_plot != nullptr) { SpecialLatex_plot->Draw("SAME"); }
      if (SpecialLatex_plot2 != nullptr) { SpecialLatex_plot2->Draw("SAME"); }
      TLatex* txtInfoE = DrawTLatex(0.657, 0.723, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
      TLegend* legend = new TLegend(0.42,0.17,0.72,0.32); //LEFT x, y, Right x, y
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.05);
      legend->AddEntry(hUSS_Y[j], "K^{+}#pi^{-}", "lpf");
      legend->AddEntry(hLSS_Y[j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}+#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
      legend->Draw("SAME");

      plots->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");
      // ========================================================================================
      // Drawing Substraction
      hUSS_Y_SUB[j] = (TH1D*)hUSS_Y[j]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_SUB[j] = (TH1D*)hLSS_Y[j]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      TCanvas* substraction = new TCanvas(Form("substraction_%d", j), "", 550, 440);
      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.08);
      substraction->SetBottomMargin(0.14);

      hUSS_Y_SUB[j]->Add(hLSS_Y_SUB[j], -1);
      hUSS_Y_SUB[j]->Draw();
      
      HistoAxisTitles(hUSS_Y_SUB[j], "#it{M}_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]", "Counts");
      hUSS_Y_SUB[j]->SetTitle(Form(""));

      TLatex* txtInfoa = DrawTLatex(0.55, 0.88, Kind+" yield", 0.04, 62, 1);
      TLatex* txtInfob = DrawTLatex(0.55, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoc = DrawTLatex(0.55, 0.83, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1); 
      if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
      TLatex* txtInfoe = DrawTLatex(0.55, 0.741, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);

      substraction->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh));
    }
  
  // fitting of K*(892)0
  //=======
  // main
  //=======
  TF1* fExtraBg [nPtBins];
  TF1* fkstar [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];

  TF1* fExtraBgTmp [nPtBins];
  TF1* fkstarTmp [nPtBins];
  TF1* fBgTmp [nPtBins];
  TF1* fFit_Tmp [nPtBins];

  TF1* fExtraBgTmp2 [nPtBins];
  TF1* fkstarTmp2 [nPtBins];
  TF1* fBgTmp2 [nPtBins];
  TF1* fFit_Tmp2 [nPtBins];
  
  //=======
  // POST
  //=======
  TF1* fBg_POST [nPtBins];
  TF1* fkstar_POST [nPtBins];
  
  TF1* fBgTmp_POST [nPtBins];
  TF1* fkstarTmp_POST [nPtBins];
  
  TF1* fBgTmp_POST2 [nPtBins];
  TF1* fkstarTmp_POST2 [nPtBins];
  
  //=======
  // Range
  //=======
  double kstarMass = 0.89581;  
  double ExtraMinRange = 0.75;
  double ExtraMaxRange = 0.84;
  double kstarMinRange = 0.75;
  double kstarMaxRange = 1.05;
  double BgMinRange = 0.83;
  double BgMaxRange = 1.25;

  TH1D* hYield = new TH1D("binYield","binYield;pT;N", nPtBins,xBins); 

  if (Kind == "Inclusive"){
    // Part 1. At low pT
    for (int k=0; k<5; k++)
      {
	hUSS_Y_kstarSig[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Bg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_ExtraBg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Fit[k]->Add(hUSS_Y_SUB[k]);

	if (k==0){
	  for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.77); bin <= hUSS_Y_ExtraBg[k]->FindBin(1.02); bin++)
	    hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	}
	
	else{
	  for (int bin = hUSS_Y_Bg[k]->FindBin(0.82); bin <= hUSS_Y_Bg[k]->FindBin(1.06); bin++)
	    hUSS_Y_Bg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.77); bin <= hUSS_Y_ExtraBg[k]->FindBin(0.98); bin++)
	    hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.77); bin <= hUSS_Y_kstarSig[k]->FindBin(0.825); bin++)
	    hUSS_Y_kstarSig[k]->SetBinError(bin, 0);
	}

	ptlow = xBins[k];
	pthigh = xBins[k+1];
	
	// ==============
	// 1. Extra fit
	fBg[k] = new TF1(Form("kstar_ExtraBg_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)", ExtraMinRange, BgMaxRange);
	hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
      
	// 2. Background fit
	fExtraBg[k] = new TF1(Form("kstar_bg_%i",k), "pol4", ExtraMinRange, BgMaxRange);
	hUSS_Y_ExtraBg[k]->Fit(fExtraBg[k], "SR0", "", ExtraMinRange, BgMaxRange);

	// 3. kstar
	fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
	fkstar[k]->SetParameter(1, kstarMass);
	fkstar[k]->SetParameter(2, 0.0474);
	hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	// ===============
	// fitting entire
	fFit[k] = new TF1(Form("Fullfit_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);
	
	if(k==0 || k==1){
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(4, fExtraBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(5, fExtraBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(6, fExtraBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(7, fExtraBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(8, fExtraBg[k]->GetParameter(4));		      
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));
	}
	else{
	  fFit[k]->SetParameter(0, fBg_POST[k-1]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg_POST[k-1]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg_POST[k-1]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg_POST[k-1]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBg_POST[k-1]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBg_POST[k-1]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBg_POST[k-1]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBg_POST[k-1]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBg_POST[k-1]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstar_POST[k-1]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar_POST[k-1]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar_POST[k-1]->GetParameter(2));
	}
	
	hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	// =================================================================================
	// for fitting each onea
	fBg_POST[k] = new TF1(Form("kstar_ExtraBg_POST_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
	fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(0));
	fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(1));
	fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(2));
	fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(3));
	fBg_POST[k]->SetParameter(4, fFit[k]->GetParameter(4));
	fBg_POST[k]->SetParameter(5, fFit[k]->GetParameter(5));
	fBg_POST[k]->SetParameter(6, fFit[k]->GetParameter(6));
	fBg_POST[k]->SetParameter(7, fFit[k]->GetParameter(7));
	fBg_POST[k]->SetParameter(8, fFit[k]->GetParameter(8));

	fkstar_POST[k] = new TF1(Form("kstar_kstar_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", ExtraMinRange, BgMaxRange);
	fkstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(9));
	fkstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(10));
	fkstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(11));

	TFitResultPtr ptr = hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, 1.25);
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

	HistoAxisTitles(hUSS_Y_Fit[k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
			"d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);
	
	hUSS_Y_Fit[k]->SetTitle(Form(""));

	hUSS_Y_Fit[k]->SetLineColor(1);
	hUSS_Y_Fit[k]->Draw();

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
      
	TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{WIP}", 0.05, 62, 1);
	TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, Kind+" yield", 0.04, 62, 1);
	TLatex* txtInfoCC = DrawTLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
	if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
	if (SpecialLatex_fit2 != nullptr) { SpecialLatex_fit2->Draw("SAME"); }
	TLatex* txtInfoDD = DrawTLatex(0.56, 0.74, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	TLatex* txtInfoEE = DrawTLatex(0.56, 0.68, chi2NdfText, 0.04, 42, 1);
	
	lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));

	//=============================================
	// Yield Main
	//=============================================
	double HistY = 0;
	double HistYerr = 0;
	double lowM = hUSS_Y_Fit[k]->FindBin(0.85);
	double highM = hUSS_Y_Fit[k]->FindBin(0.95);
	HistY=hUSS_Y_Fit[k]->IntegralAndError(lowM, highM, HistYerr);
	double binW = hUSS_Y_Fit[k]->GetBinWidth(50); //Any number
	double BgY = fBg_POST[k]->Integral(0.85,0.95)/binW;
	
	TMatrixDSym covTot(fFit[k]->GetNpar());
	TMatrixDSym covGG(fkstar_POST[k]->GetNpar());
	TMatrixDSym covGG_err(fBg_POST[k]->GetNpar());
	covTot = ptr->GetCovarianceMatrix();
	covGG_err = covTot.GetSub(0,8,0,8);// these have to be fBG parameters,
	Double_t Bg_err = (1/binW)*(fBg_POST[k]->IntegralError(0.85,0.95,fBg_POST[k]->GetParameters(),covGG_err.GetMatrixArray()));

	double HistCorrY = HistY-BgY;
	double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	hYield->SetBinContent(k+1, HistCorrY);
	hYield->SetBinContent(k+1, hYield->GetBinContent(k+1)/hYield->GetBinWidth(k+1));
	hYield->SetBinError(k+1, HistCorrYerr);
	Shreck();
	cout<<HistCorrY<<endl;
	cout<<Bg_err<<endl;
	
	TCanvas* numbin = new TCanvas();
	numbin->SetRightMargin(0.05);
	numbin->SetLeftMargin(0.17);
	numbin->SetBottomMargin(0.12);
	
	hYield->SetLineWidth(2);
	hYield->SetLineColor(1);
	hYield->SetMarkerStyle(20); //circle
	hYield->SetMarkerSize(1);
	hYield->SetMarkerColor(1);

	double yMin = hYield->GetMinimum()*(-0.7);
	double yMax = hYield->GetMaximum()*1.1;
	
	TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
	hdummy->SetMinimum(yMin);
	hdummy->SetMaximum(yMax);

	HistoAxisTitles(hdummy, "#it{p}_{T} [GeV/#it{c}]", "#frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}",
			0.05, 62, 1.0, 0.05, 62, 1.4);
	hdummy->SetLineWidth(2);
	hdummy->SetLineColor(1);
	
	hdummy->Draw();
	hYield->Draw("SAME");
	numbin->SaveAs("../kstar-in-jets/Plot/"+localpath+"/h_Yield.png","RECREATE");
      }

     for (int k=5; k<nPtBins; k++)
      {
	hUSS_Y_kstarSig[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Bg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_ExtraBg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Fit[k]->Add(hUSS_Y_SUB[k]);
	
	for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.83); bin <= hUSS_Y_ExtraBg[k]->FindBin(0.99); bin++)
	  hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.77); bin <= hUSS_Y_kstarSig[k]->FindBin(0.81); bin++)
	  hUSS_Y_kstarSig[k]->SetBinError(bin, 0);

	ptlow = xBins[k];
	pthigh = xBins[k+1];
	
	// ==============
	// 1. Extra fit
	fBg[k] = new TF1(Form("kstar_ExtraBg_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)", ExtraMinRange, BgMaxRange);
	hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
      
	// 2. Background fit
	fExtraBg[k] = new TF1(Form("kstar_bg_%i",k), "pol4", ExtraMinRange, BgMaxRange);
	hUSS_Y_ExtraBg[k]->Fit(fExtraBg[k], "SR0", "", ExtraMinRange, BgMaxRange);

	// 3. kstar
	fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
	fkstar[k]->SetParameter(1, kstarMass);
	fkstar[k]->SetParameter(2, 0.0474);
	hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	// ===============
	// fitting entire
	fFit[k] = new TF1(Form("Fullfit_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);
	
	if(k==5){
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	
	  fFit[k]->SetParameter(4, fExtraBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(5, fExtraBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(6, fExtraBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(7, fExtraBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(8, fExtraBg[k]->GetParameter(4));
			      
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));
	}
	else{
	  fFit[k]->SetParameter(0, fBg_POST[k-1]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg_POST[k-1]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg_POST[k-1]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg_POST[k-1]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBg_POST[k-1]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBg_POST[k-1]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBg_POST[k-1]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBg_POST[k-1]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBg_POST[k-1]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstar_POST[k-1]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar_POST[k-1]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar_POST[k-1]->GetParameter(2));
	}
	
	hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	// =================================================================================
	// for fitting each onea
	fBg_POST[k] = new TF1(Form("kstar_ExtraBg_POST_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
	fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(0));
	fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(1));
	fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(2));
	fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(3));
	fBg_POST[k]->SetParameter(4, fFit[k]->GetParameter(4));
	fBg_POST[k]->SetParameter(5, fFit[k]->GetParameter(5));
	fBg_POST[k]->SetParameter(6, fFit[k]->GetParameter(6));
	fBg_POST[k]->SetParameter(7, fFit[k]->GetParameter(7));
	fBg_POST[k]->SetParameter(8, fFit[k]->GetParameter(8));

	fkstar_POST[k] = new TF1(Form("kstar_kstar_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", ExtraMinRange, BgMaxRange);
	fkstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(9));
	fkstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(10));
	fkstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(11));

	TFitResultPtr ptr = hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, 1.25);
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

	HistoAxisTitles(hUSS_Y_Fit[k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
			"d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);
	
	hUSS_Y_Fit[k]->SetTitle(Form(""));

	hUSS_Y_Fit[k]->SetLineColor(1);
	hUSS_Y_Fit[k]->Draw();

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
      
	TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{WIP}", 0.05, 62, 1);
	TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, Kind+" yield", 0.04, 62, 1);
	TLatex* txtInfoCC = DrawTLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
	if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
	if (SpecialLatex_fit2 != nullptr) { SpecialLatex_fit2->Draw("SAME"); }
	TLatex* txtInfoDD = DrawTLatex(0.56, 0.74, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	TLatex* txtInfoEE = DrawTLatex(0.56, 0.68, chi2NdfText, 0.04, 42, 1);
	
	lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
	//=============================================
	// Yield Main
	//=============================================
	double HistY = 0;
	double HistYerr = 0;
	double lowM = hUSS_Y_Fit[k]->FindBin(0.85);
	double highM = hUSS_Y_Fit[k]->FindBin(0.95);
	HistY=hUSS_Y_Fit[k]->IntegralAndError(lowM, highM, HistYerr);
	double binW = hUSS_Y_Fit[k]->GetBinWidth(50); //Any number
	double BgY = fBg_POST[k]->Integral(0.85,0.95)/binW;
	
	TMatrixDSym covTot(fFit[k]->GetNpar());
	TMatrixDSym covGG(fkstar_POST[k]->GetNpar());
	TMatrixDSym covGG_err(fBg_POST[k]->GetNpar());
	covTot = ptr->GetCovarianceMatrix();
	covGG_err = covTot.GetSub(0,8,0,8);// these have to be fBG parameters,
	Double_t Bg_err = (1/binW)*(fBg_POST[k]->IntegralError(0.85,0.95,fBg_POST[k]->GetParameters(),covGG_err.GetMatrixArray()));

	double HistCorrY = HistY-BgY;
	double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	hYield->SetBinContent(k+1, HistCorrY);
	hYield->SetBinContent(k+1, hYield->GetBinContent(k+1)/hYield->GetBinWidth(k+1));
	hYield->SetBinError(k+1, HistCorrYerr);
	Shreck();
	cout<<HistCorrY<<endl;
	cout<<Bg_err<<endl;
	
	TCanvas* numbin = new TCanvas();
	numbin->SetRightMargin(0.05);
	numbin->SetLeftMargin(0.17);
	numbin->SetBottomMargin(0.12);
	
	hYield->SetLineWidth(2);
	hYield->SetLineColor(1);
	hYield->SetMarkerStyle(20); //circle
	hYield->SetMarkerSize(1);
	hYield->SetMarkerColor(1);

	double yMin = hYield->GetMinimum()*(-0.7);
	double yMax = hYield->GetMaximum()*1.1;
	
	TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
	hdummy->SetMinimum(yMin);
	hdummy->SetMaximum(yMax);

	HistoAxisTitles(hdummy, "#it{p}_{T} [GeV/#it{c}]", "#frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}",
			0.05, 62, 1.0, 0.05, 62, 1.4);
	hdummy->SetLineWidth(2);
	hdummy->SetLineColor(1);
	
	hdummy->Draw();
	hYield->Draw("SAME");
	numbin->SaveAs("../kstar-in-jets/Plot/"+localpath+"/h_Yield.png","RECREATE");
      }
  }
  // Jet
  else if (Kind == "Inside"){
    for (int k=6; k>-1; k--)
      {
	hUSS_Y_kstarSig[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Bg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_ExtraBg[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_Fit[k]->Add(hUSS_Y_SUB[k]);

	hUSS_Y_kstarSigTmp[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_BgTmp[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_ExtraBgTmp[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_FitTmp[k]->Add(hUSS_Y_SUB[k]);

	hUSS_Y_kstarSigTmp2[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_BgTmp2[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_ExtraBgTmp2[k]->Add(hUSS_Y_SUB[k]);
	hUSS_Y_FitTmp2[k]->Add(hUSS_Y_SUB[k]);

	ptlow = xBins[k];
	pthigh = xBins[k+1];

	//======
	//Origin
	//======
	fBg[k] = new TF1(Form("kstar_ExtraBg_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)", ExtraMinRange, BgMaxRange);
	fExtraBg[k]= new TF1(Form("kstar_ExtraBg_%i",k), "pol4", ExtraMinRange, BgMaxRange);
	fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
	fkstar[k]->SetParameter(1, kstarMass);
	fkstar[k]->SetParameter(2, 0.0474);

	//======
	//Tmp1
	//======
	fBgTmp[k] = new TF1(Form("kstar_ExtraBgTmp_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)", ExtraMinRange, BgMaxRange);
	fExtraBgTmp[k]= new TF1(Form("kstar_ExtraBgTmp_%i",k), "pol4", ExtraMinRange, BgMaxRange);
	fkstarTmp[k] = new TF1(Form("kstarTmp_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
	fkstarTmp[k]->SetParameter(1, kstarMass);
	fkstarTmp[k]->SetParameter(2, 0.0474);
	//======
	//Tmp2
	//======
	fBgTmp2[k] = new TF1(Form("kstar_ExtraBgTmp2_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)", ExtraMinRange, BgMaxRange);
	fExtraBgTmp2[k]= new TF1(Form("kstar_ExtraBgTmp2_%i",k), "pol4", ExtraMinRange, BgMaxRange);
	fkstarTmp2[k] = new TF1(Form("kstarTmp2_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
	fkstarTmp2[k]->SetParameter(1, kstarMass);
	fkstarTmp2[k]->SetParameter(2, 0.0474);
	

	fFit[k] = new TF1(Form("Fullfit_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);
	fFit_Tmp[k] = new TF1(Form("FullfitTmp_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);
	fFit_Tmp2[k] = new TF1(Form("FullfitTmp2_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);


	
	if (k==8){
	  for (int bin = hUSS_Y_Bg[k]->FindBin(0.90); bin <= hUSS_Y_Bg[k]->FindBin(0.96); bin++)
	    hUSS_Y_Bg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.84); bin <= hUSS_Y_ExtraBg[k]->FindBin(0.94); bin++)
	    hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.765); bin <= hUSS_Y_kstarSig[k]->FindBin(0.835); bin++)
	    hUSS_Y_kstarSig[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_ExtraBg[k]->Fit(fExtraBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(4, fExtraBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(5, fExtraBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(6, fExtraBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(7, fExtraBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(8, fExtraBg[k]->GetParameter(4));
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));
	}

	else if (k==5){
	  for (int bin = hUSS_Y_Bg[k]->FindBin(0.76); bin <= hUSS_Y_Bg[k]->FindBin(0.99); bin++)
	    hUSS_Y_Bg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.84); bin <= hUSS_Y_ExtraBg[k]->FindBin(0.94); bin++)
	    hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.77); bin <= hUSS_Y_kstarSig[k]->FindBin(0.83); bin++)
	    hUSS_Y_kstarSig[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_ExtraBg[k]->Fit(fExtraBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(4, fExtraBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(5, fExtraBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(6, fExtraBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(7, fExtraBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(8, fExtraBg[k]->GetParameter(4));
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));
	}

	else if (k==4){
	  for (int bin = hUSS_Y_Bg[k]->FindBin(0.76); bin <= hUSS_Y_Bg[k]->FindBin(0.99); bin++)
	    hUSS_Y_Bg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBg[k]->FindBin(0.84); bin <= hUSS_Y_ExtraBg[k]->FindBin(0.94); bin++)
	    hUSS_Y_ExtraBg[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.77); bin <= hUSS_Y_kstarSig[k]->FindBin(0.84); bin++)
	    hUSS_Y_kstarSig[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_ExtraBg[k]->Fit(fExtraBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(4, fExtraBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(5, fExtraBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(6, fExtraBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(7, fExtraBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(8, fExtraBg[k]->GetParameter(4));
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));

	  //===================
	  //Temp1 : 1.5~1.8
	  //===================
	  for (int bin = hUSS_Y_BgTmp[k]->FindBin(0.76); bin <= hUSS_Y_BgTmp[k]->FindBin(0.99); bin++)
	    hUSS_Y_BgTmp[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBgTmp[k]->FindBin(0.83); bin <= hUSS_Y_ExtraBgTmp[k]->FindBin(0.94); bin++)
	    hUSS_Y_ExtraBgTmp[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSigTmp[k]->FindBin(0.77); bin <= hUSS_Y_kstarSigTmp[k]->FindBin(0.84); bin++)
	    hUSS_Y_kstarSigTmp[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_BgTmp[k]->Fit(fBgTmp[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_ExtraBgTmp[k]->Fit(fExtraBgTmp[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSigTmp[k]->Fit(fkstarTmp[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	  fFit_Tmp[k]->SetParameter(0, fBgTmp[k]->GetParameter(0));
	  fFit_Tmp[k]->SetParameter(1, fBgTmp[k]->GetParameter(1));
	  fFit_Tmp[k]->SetParameter(2, fBgTmp[k]->GetParameter(2));
	  fFit_Tmp[k]->SetParameter(3, fBgTmp[k]->GetParameter(3));
	  fFit_Tmp[k]->SetParameter(4, fExtraBgTmp[k]->GetParameter(0));
	  fFit_Tmp[k]->SetParameter(5, fExtraBgTmp[k]->GetParameter(1));
	  fFit_Tmp[k]->SetParameter(6, fExtraBgTmp[k]->GetParameter(2));
	  fFit_Tmp[k]->SetParameter(7, fExtraBgTmp[k]->GetParameter(3));
	  fFit_Tmp[k]->SetParameter(8, fExtraBgTmp[k]->GetParameter(4));
	  fFit_Tmp[k]->SetParameter(9, fkstarTmp[k]->GetParameter(0));
	  fFit_Tmp[k]->SetParameter(10, fkstarTmp[k]->GetParameter(1));
	  fFit_Tmp[k]->SetParameter(11, fkstarTmp[k]->GetParameter(2));

	  //===================
	  //Temp2 :1.8~2.2
	  //===================
	  for (int bin = hUSS_Y_BgTmp2[k]->FindBin(0.76); bin <= hUSS_Y_BgTmp2[k]->FindBin(0.99); bin++)
	    hUSS_Y_BgTmp2[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_ExtraBgTmp2[k]->FindBin(0.82); bin <= hUSS_Y_ExtraBgTmp2[k]->FindBin(0.94); bin++)
	    hUSS_Y_ExtraBgTmp2[k]->SetBinError(bin, 0);
	  for (int bin = hUSS_Y_kstarSigTmp2[k]->FindBin(0.76); bin <= hUSS_Y_kstarSigTmp2[k]->FindBin(0.79); bin++)
	    hUSS_Y_kstarSigTmp2[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_BgTmp2[k]->Fit(fBgTmp2[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_ExtraBgTmp2[k]->Fit(fExtraBgTmp2[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSigTmp2[k]->Fit(fkstarTmp2[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	  fFit_Tmp2[k]->SetParameter(0, fBgTmp2[k]->GetParameter(0));
	  fFit_Tmp2[k]->SetParameter(1, fBgTmp2[k]->GetParameter(1));
	  fFit_Tmp2[k]->SetParameter(2, fBgTmp2[k]->GetParameter(2));
	  fFit_Tmp2[k]->SetParameter(3, fBgTmp2[k]->GetParameter(3));
	  fFit_Tmp2[k]->SetParameter(4, fExtraBgTmp2[k]->GetParameter(0));
	  fFit_Tmp2[k]->SetParameter(5, fExtraBgTmp2[k]->GetParameter(1));
	  fFit_Tmp2[k]->SetParameter(6, fExtraBgTmp2[k]->GetParameter(2));
	  fFit_Tmp2[k]->SetParameter(7, fExtraBgTmp2[k]->GetParameter(3));
	  fFit_Tmp2[k]->SetParameter(8, fExtraBgTmp2[k]->GetParameter(4));
	  fFit_Tmp2[k]->SetParameter(9, fkstarTmp2[k]->GetParameter(0));
	  fFit_Tmp2[k]->SetParameter(10, fkstarTmp2[k]->GetParameter(1));
	  fFit_Tmp2[k]->SetParameter(11, fkstarTmp2[k]->GetParameter(2));
	}
	
	else if (k==3){
	  //===================
	  //Temp1 :1.5~1.8
	  //===================
	  fFit_Tmp[k]->SetParameter(0, fBgTmp_POST[k+1]->GetParameter(0));
	  fFit_Tmp[k]->SetParameter(1, fBgTmp_POST[k+1]->GetParameter(1));
	  fFit_Tmp[k]->SetParameter(2, fBgTmp_POST[k+1]->GetParameter(2));
	  fFit_Tmp[k]->SetParameter(3, fBgTmp_POST[k+1]->GetParameter(3));
	  fFit_Tmp[k]->SetParameter(4, fBgTmp_POST[k+1]->GetParameter(4));
	  fFit_Tmp[k]->SetParameter(5, fBgTmp_POST[k+1]->GetParameter(5));
	  fFit_Tmp[k]->SetParameter(6, fBgTmp_POST[k+1]->GetParameter(6));
	  fFit_Tmp[k]->SetParameter(7, fBgTmp_POST[k+1]->GetParameter(7));
	  fFit_Tmp[k]->SetParameter(8, fBgTmp_POST[k+1]->GetParameter(8));
	  fFit_Tmp[k]->SetParameter(9, fkstarTmp_POST[k+1]->GetParameter(0));
	  fFit_Tmp[k]->SetParameter(10, fkstarTmp_POST[k+1]->GetParameter(1));
	  fFit_Tmp[k]->SetParameter(11, fkstarTmp_POST[k+1]->GetParameter(2));

	  //===================
	  //Temp2 : 1.8~2.2
	  //===================
	  fFit[k]->SetParameter(0, fBgTmp_POST2[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBgTmp_POST2[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBgTmp_POST2[k+1]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBgTmp_POST2[k+1]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBgTmp_POST2[k+1]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBgTmp_POST2[k+1]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBgTmp_POST2[k+1]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBgTmp_POST2[k+1]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBgTmp_POST2[k+1]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstarTmp_POST2[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstarTmp_POST2[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstarTmp_POST2[k+1]->GetParameter(2));
	}

	else if (k==2){
	  //===================
	  //Temp1 :1.5~1.8
	  //===================
	  fFit[k]->SetParameter(0, fBgTmp_POST[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBgTmp_POST[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBgTmp_POST[k+1]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBgTmp_POST[k+1]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBgTmp_POST[k+1]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBgTmp_POST[k+1]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBgTmp_POST[k+1]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBgTmp_POST[k+1]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBgTmp_POST[k+1]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstarTmp_POST[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstarTmp_POST[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstarTmp_POST[k+1]->GetParameter(2));
	}

	
	else if (k==0 || k==1){
	  fBg[k] = new TF1(Form("kstar_ExtraBg_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)+pol4(4)", ExtraMinRange, BgMaxRange);
	  
	  for (int bin = hUSS_Y_kstarSig[k]->FindBin(0.75); bin <= hUSS_Y_kstarSig[k]->FindBin(0.85); bin++)
	    hUSS_Y_kstarSig[k]->SetBinError(bin, 0);
	  
	  hUSS_Y_Bg[k]->Fit(fBg[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  hUSS_Y_kstarSig[k]->Fit(fkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);
	  
	  fFit[k]->SetParameter(0, fBg[k]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBg[k]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBg[k]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBg[k]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBg[k]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBg[k]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar[k]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar[k]->GetParameter(2));
	}
	
	else{
	  fFit[k]->SetParameter(0, fBg_POST[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(1, fBg_POST[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(2, fBg_POST[k+1]->GetParameter(2));
	  fFit[k]->SetParameter(3, fBg_POST[k+1]->GetParameter(3));
	  fFit[k]->SetParameter(4, fBg_POST[k+1]->GetParameter(4));
	  fFit[k]->SetParameter(5, fBg_POST[k+1]->GetParameter(5));
	  fFit[k]->SetParameter(6, fBg_POST[k+1]->GetParameter(6));
	  fFit[k]->SetParameter(7, fBg_POST[k+1]->GetParameter(7));
	  fFit[k]->SetParameter(8, fBg_POST[k+1]->GetParameter(8));
	  fFit[k]->SetParameter(9, fkstar_POST[k+1]->GetParameter(0));
	  fFit[k]->SetParameter(10, fkstar_POST[k+1]->GetParameter(1));
	  fFit[k]->SetParameter(11, fkstar_POST[k+1]->GetParameter(2));
	}
	
	hUSS_Y_FitTmp2[k]->Fit(fFit_Tmp2[k], "SR0", "", ExtraMinRange, BgMaxRange);
	hUSS_Y_FitTmp[k]->Fit(fFit_Tmp[k], "SR0", "", ExtraMinRange, BgMaxRange);
	hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);
	
	// First
	fBg_POST[k] = new TF1(Form("kstar_ExtraBg_POST_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
	fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(0));
	fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(1));
	fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(2));
	fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(3));
	fBg_POST[k]->SetParameter(4, fFit[k]->GetParameter(4));
	fBg_POST[k]->SetParameter(5, fFit[k]->GetParameter(5));
	fBg_POST[k]->SetParameter(6, fFit[k]->GetParameter(6));
	fBg_POST[k]->SetParameter(7, fFit[k]->GetParameter(7));
	fBg_POST[k]->SetParameter(8, fFit[k]->GetParameter(8));

	fkstar_POST[k] = new TF1(Form("kstar_kstar_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", ExtraMinRange, BgMaxRange);
	fkstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(9));
	fkstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(10));
	fkstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(11));

	//=============================================================
	//Tmp
	//=============================================================
	fBgTmp_POST[k] = new TF1(Form("kstar_ExtraBgTMP_POST_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
	fBgTmp_POST[k]->SetParameter(0, fFit_Tmp[k]->GetParameter(0));
	fBgTmp_POST[k]->SetParameter(1, fFit_Tmp[k]->GetParameter(1));
	fBgTmp_POST[k]->SetParameter(2, fFit_Tmp[k]->GetParameter(2));
	fBgTmp_POST[k]->SetParameter(3, fFit_Tmp[k]->GetParameter(3));
	fBgTmp_POST[k]->SetParameter(4, fFit_Tmp[k]->GetParameter(4));
	fBgTmp_POST[k]->SetParameter(5, fFit_Tmp[k]->GetParameter(5));
	fBgTmp_POST[k]->SetParameter(6, fFit_Tmp[k]->GetParameter(6));
	fBgTmp_POST[k]->SetParameter(7, fFit_Tmp[k]->GetParameter(7));
	fBgTmp_POST[k]->SetParameter(8, fFit_Tmp[k]->GetParameter(8));

	fkstarTmp_POST[k] = new TF1(Form("kstar_kstarTMP_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", ExtraMinRange, BgMaxRange);
	fkstarTmp_POST[k]->SetParameter(0, fFit_Tmp[k]->GetParameter(9));
	fkstarTmp_POST[k]->SetParameter(1, fFit_Tmp[k]->GetParameter(10));
	fkstarTmp_POST[k]->SetParameter(2, fFit_Tmp[k]->GetParameter(11));

	//=============================================================
	//Tmp2
	//=============================================================
	fBgTmp_POST2[k] = new TF1(Form("kstar_ExtraBgTMP_POST2_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
	fBgTmp_POST2[k]->SetParameter(0, fFit_Tmp2[k]->GetParameter(0));
	fBgTmp_POST2[k]->SetParameter(1, fFit_Tmp2[k]->GetParameter(1));
	fBgTmp_POST2[k]->SetParameter(2, fFit_Tmp2[k]->GetParameter(2));
	fBgTmp_POST2[k]->SetParameter(3, fFit_Tmp2[k]->GetParameter(3));
	fBgTmp_POST2[k]->SetParameter(4, fFit_Tmp2[k]->GetParameter(4));
	fBgTmp_POST2[k]->SetParameter(5, fFit_Tmp2[k]->GetParameter(5));
	fBgTmp_POST2[k]->SetParameter(6, fFit_Tmp2[k]->GetParameter(6));
	fBgTmp_POST2[k]->SetParameter(7, fFit_Tmp2[k]->GetParameter(7));
	fBgTmp_POST2[k]->SetParameter(8, fFit_Tmp2[k]->GetParameter(8));

	fkstarTmp_POST2[k] = new TF1(Form("kstar_kstarTMP_POST2_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", ExtraMinRange, BgMaxRange);
	fkstarTmp_POST2[k]->SetParameter(0, fFit_Tmp2[k]->GetParameter(9));
	fkstarTmp_POST2[k]->SetParameter(1, fFit_Tmp2[k]->GetParameter(10));
	fkstarTmp_POST2[k]->SetParameter(2, fFit_Tmp2[k]->GetParameter(11));


	

	TFitResultPtr ptr = hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, 1.25);
	double Chi2 = ptr->Chi2();
	int NDF = ptr->Ndf();
	double chi2Ndf = Chi2/NDF;
	string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);
      
	// Draw
	TCanvas* lines = new TCanvas(Form("line_%d",k), "", 940, 800);
	lines->SetLeftMargin(0.12);
	lines->SetRightMargin(0.1);
	lines->SetTopMargin(0.08);
	lines->SetBottomMargin(0.14);

	HistoAxisTitles(hUSS_Y_Fit[k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
			"d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.5);
	hUSS_Y_Fit[k]->SetTitle(Form(""));
	hUSS_Y_Fit[k]->SetLineColor(1);
	hUSS_Y_Fit[k]->Draw();

	// total
	fFit[k]->SetLineColor(2);
	fFit[k]->SetLineWidth(4);
	fFit[k]->Draw("SAME");

	//expo
	fBg_POST[k]->SetLineColor(4);
	fBg_POST[k]->SetLineWidth(4);
	fBg_POST[k]->Draw("SAME");

	//kstar
	fkstar_POST[k]->SetLineColor(6);
	fkstar_POST[k]->SetLineWidth(4);
	fkstar_POST[k]->Draw("SAME");
      
	TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{WIP}", 0.05, 62, 1);
	TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, Kind+" yield", 0.04, 62, 1);
	TLatex* txtInfoCC = DrawTLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
	if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
	if (SpecialLatex_fit2 != nullptr) { SpecialLatex_fit2->Draw("SAME"); }
	TLatex* txtInfoDD = DrawTLatex(0.56, 0.74, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	TLatex* txtInfoEE = DrawTLatex(0.56, 0.68, chi2NdfText, 0.04, 42, 1);

	lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));

	double HistY = 0;
	double HistYerr = 0;
	double lowM = hUSS_Y_Fit[k]->FindBin(0.85);
	double highM = hUSS_Y_Fit[k]->FindBin(0.95);
	HistY=hUSS_Y_Fit[k]->IntegralAndError(lowM, highM, HistYerr);
	double binW = hUSS_Y_Fit[k]->GetBinWidth(50); //Any number
	double BgY = fBg_POST[k]->Integral(0.85,0.95)/binW;
	
	TMatrixDSym covTot(fFit[k]->GetNpar());
	TMatrixDSym covGG(fkstar_POST[k]->GetNpar());
	TMatrixDSym covGG_err(fBg_POST[k]->GetNpar());
	covTot = ptr->GetCovarianceMatrix();
	covGG_err = covTot.GetSub(0,8,0,8);// these have to be fBG parameters,
	Double_t Bg_err = (1/binW)*(fBg_POST[k]->IntegralError(0.85,0.95,fBg_POST[k]->GetParameters(),covGG_err.GetMatrixArray()));

	double HistCorrY = HistY-BgY;
	double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	hYield->SetBinContent(k+1, HistCorrY);
	hYield->SetBinContent(k+1, hYield->GetBinContent(k+1)/hYield->GetBinWidth(k+1));
	hYield->SetBinError(k+1, HistCorrYerr);
	Shreck();
	cout<<HistCorrY<<endl;
	cout<<Bg_err<<endl;

	TCanvas* numbin = new TCanvas();
	numbin->SetRightMargin(0.05);
	numbin->SetLeftMargin(0.17);
	numbin->SetBottomMargin(0.12);
	
	hYield->SetLineWidth(2);
	hYield->SetLineColor(2);
	hYield->SetMarkerStyle(23);
	hYield->SetMarkerSize(1);
	hYield->SetMarkerColor(2);

	double yMin = hYield->GetMinimum()*(-0.7);
	double yMax = hYield->GetMaximum()*1.1;
	
	TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
	hdummy->SetMinimum(yMin);
	hdummy->SetMaximum(yMax);

	HistoAxisTitles(hdummy, "#it{p}_{T} [GeV/#it{c}]", "#frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}",
			0.05, 62, 1.0, 0.05, 62, 1.4);
	hdummy->SetLineWidth(2);
	hdummy->SetLineColor(1);
	
	hdummy->Draw();
	hYield->Draw("SAME");
	numbin->SaveAs("../kstar-in-jets/Plot/"+localpath+"/h_Yield.png","RECREATE");
      }
    /* } */
 
    
  }//jet
}//DrawBW


void Shreck() {
  const char* shrekAscii = R"(
    ⢀⡴⠑⡄⠀⠀⠀⠀⠀⠀⠀⣀⣀⣤⣤⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
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
    )";

    std::cout << shrekAscii << std::endl;
}
