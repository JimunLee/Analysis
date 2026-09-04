#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/PdfFuncMathCore.h"
#include <TLine.h>
#include <iostream>
#include <fstream>

using namespace std;

// Draw function

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

// crystsalball_function

double crystalball_function(double *x, double *par) {
  double alpha = par[0];
  double n = par[1];
  double mean = par[2];
  double sigma = par[3];

  if (sigma < 0.) return 0.;
  double z = (x[0] - mean) / sigma;
  if (alpha < 0) z = -z;
  double abs_alpha = abs(alpha);

  if (z > -abs_alpha)
    return std::exp(-0.5 * z * z);  // Gaussian

  double nDivAlpha = n / abs_alpha;
  double AA = exp(-0.5 * abs_alpha * abs_alpha);
  double B = nDivAlpha - abs_alpha;
  double arg = nDivAlpha / (B - z);

  return AA * pow(arg, n);
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
  Double_t nBins = 500.0;

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
  
  TH1* hUSS_Y_Bkg [nPtBins];
  TH1* hUSS_Y_Peak [nPtBins];
  TH1* hUSS_Y_Fit [nPtBins];

  for (Int_t i=0; i<nPtBins; i++)
    {
      hLSS_Y[i] = new TH1D(Form("hLSS_Y_%i", i), Form("hLSS_Y_%i", i), nBins, imlow, imtop);
      hUSS_Y[i] = new TH1D(Form("hUSS_Y_%i", i), Form("hUSS_Y_%i", i), nBins, imlow, imtop);

      hUSS_Y_Bkg[i] = new TH1D(Form("hUSS_Y_Bkg_%i", i), Form("hUSS_Y_Bkg_%i", i), nBins, imlow, imtop); 
      hUSS_Y_Peak[i] = new TH1D(Form("hUSS_Y_Peak_%i", i), Form("hUSS_Y_Peak_%i", i), nBins, imlow, imtop);
      hUSS_Y_Fit[i] = new TH1D(Form("hUSS_Y_Fit_%i", i), Form("hUSS_Y_Fit_%i", i), nBins, imlow, imtop);
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
      Double_t Nlow = 1.2;
      Double_t Nhigh = 1.25;
      // Scale
      Double_t Scale = hUSS_Y[j]->Integral(hUSS_Y[j]->FindBin(Nlow), hUSS_Y[j]->FindBin(Nhigh))
	/ hLSS_Y[j]->Integral(hLSS_Y[j]->FindBin(Nlow), hLSS_Y[j]->FindBin(Nhigh));
      hLSS_Y[j]->Scale(Scale);

      hUSS_Y[j]->Rebin(ReBinNum);
      hUSS_Y[j]->Scale(1./ReBinNum);
      hLSS_Y[j]->Rebin(ReBinNum);
      hLSS_Y[j]->Scale(1./ReBinNum);

      hUSS_Y[j]->Scale(1.0/nEvents);
      hLSS_Y[j]->Scale(1.0/nEvents);
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
  TF1* fBW_kstar [nPtBins];
  TF1* fCB [nPtBins];
  TF1Convolution* fconv [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit_kstar [nPtBins];
  TF1* fFit [nPtBins];
  
  double Min=0.0;
  double Max=10.0;

  double kstarMassMin = 0.89581-0.00019;
  double kstarMassMax = 0.89581+0.00019;
  
  double kstarMinRange = 0.76;
  double kstarMaxRange = 1.04;
  
  double Bglow = 0.75;
  double Bghigh = 1.25;

  for (int k=1; k<nPtBins+1; k++)
    {
      hUSS_Y_Fit[nPtBins-k]->Add(hUSS_Y_SUB[nPtBins-k]);
      
      // add
      for(int bin = hUSS_Y_Fit[nPtBins-k]->FindBin(0.758); bin<=hUSS_Y_Fit[nPtBins-k]->FindBin(0.827); bin++)
	hUSS_Y_Fit[nPtBins-k]->SetBinError(bin,0);
      
      // fill bin
      int startBin = hUSS_Y_Fit[nPtBins-k]->FindBin(0.87);
      int endBin = hUSS_Y_Fit[nPtBins-k]->FindBin(0.79);
      double refValue = hUSS_Y_Fit[nPtBins-k]->GetBinContent(startBin);
      
      for(int bin = startBin -1; bin >= endBin; --bin){
	double currentValue = hUSS_Y_Fit[nPtBins-k]->GetBinContent(bin);

	if (currentValue >= refValue) {
	  hUSS_Y_Fit[nPtBins-k]->SetBinError(bin, 0);
	}
      }
      ptlow = xBins[nPtBins-k];
      pthigh = xBins[nPtBins-k+1];
      //===========================================================================================================================
      // ==============
      // fitting kstar
      fBW_kstar[nPtBins-k] = new TF1(Form("kstar_fit_%i", nPtBins-k), "[0]*TMath::BreitWigner(x, [1], [2])", kstarMinRange, kstarMaxRange);
      fBW_kstar[nPtBins-k]->SetParLimits(0, Min, Max);
      fBW_kstar[nPtBins-k]->SetParLimits(1, kstarMassMin, kstarMassMax);
      fBW_kstar[nPtBins-k]->FixParameter(2, 0.0474);

      // plus function 
      // fCB_kstar[nPtBins-k] = new TF1(Form("fFit_kstar_%i", nPtBins-k), crystalball_function, kstarMinRange, kstarMaxRange, 4);
      // fCB_kstar[nPtBins-k]->SetParLimits(2, kstarMassMin, kstarMassMax);
      //  fCB_kstar[nPtBins-k]->FixParameter(3, 0.0474); 


      fCB[nPtBins-k] = new TF1(Form("fCB_fiti_%i", nPtBins-k), "gaus", kstarMinRange, kstarMaxRange);
      fCB[nPtBins-k]->GetParameters();
      // convolution
      TF1Convolution *fconv = new TF1Convolution(fBW_kstar[nPtBins-k], fCB[nPtBins-k], Bglow, Bghigh, true);
      fFit_kstar[nPtBins-k] = new TF1(Form("Total_kstar_%i", nPtBins-k), fCB[nPtBins-k], Bglow, Bghigh, fconv->GetNpar());

      // histogram fitting
      hUSS_Y_Fit[nPtBins-k]->Fit(fFit_kstar[nPtBins-k], "0", "", kstarMinRange, kstarMaxRange);
      
      
      /* fFit_kstar[nPtBins-k]->SetParameter(0, fBW_kstar[nPtBins-k]->GetParameter(0)); */
      /* fFit_kstar[nPtBins-k]->SetParameter(1, fBW_kstar[nPtBins-k]->GetParameter(1)); */
      /* fFit_kstar[nPtBins-k]->SetParameter(2, fBW_kstar[nPtBins-k]->GetParameter(2)); */
      /* fFit_kstar[nPtBins-k]->SetParLimits(3, Min, Max); */
      
      // =====================
      // fitting background
      /* fBg[nPtBins-k] = new TF1(Form("BgFit_%i", nPtBins-k), "pol3", Bglow, Bghigh); */

      /* hUSS_Y_Fit[nPtBins-k]->Fit(fBg[nPtBins-k], "0", "", Bglow, Bghigh); */

      // ===============
      // fitting entire
      /* fFit[nPtBins-k] = new TF1(Form("Fullfit_%i", nPtBins-k), "[0]*TMath::BreitWigner(x, [1], [2]) + [3]*fConv + pol3(4)", Bglow, Bghigh); */
      
      /* fFit[nPtBins-k]->SetParameter(0, fFit_kstar[nPtBins-k]->GetParameter(0)); */
      /* fFit[nPtBins-k]->SetParameter(1, fFit_kstar[nPtBins-k]->GetParameter(1)); */
      /* fFit[nPtBins-k]->SetParameter(2, fFit_kstar[nPtBins-k]->GetParameter(2)); */
      /* fFit[nPtBins-k]->SetParameter(3, fFit_kstar[nPtBins-k]->GetParameter(3)); */
       
      /* fFit[nPtBins-k]->SetParameter(4, fBg[nPtBins-k]->GetParameter(0)); */
      /* fFit[nPtBins-k]->SetParameter(5, fBg[nPtBins-k]->GetParameter(1)); */
      /* fFit[nPtBins-k]->SetParameter(6, fBg[nPtBins-k]->GetParameter(2)); */
      /* fFit[nPtBins-k]->SetParameter(7, fBg[nPtBins-k]->GetParameter(3)); */

      /* hUSS_Y_Fit[nPtBins-k]->GetXaxis()->SetRangeUser(0.75, 1.15); */
      // hUSS_Y_Fit[nPtBins-k]->Fit(fFit[nPtBins-k], "SR0", "", Bglow, Bghigh); */

      // =================
      // =================================================================================
      // for fitting each one
      fFit_kstar[nPtBins-k]->SetParameter(0, fFit[nPtBins-k]->GetParameter(0));
      fFit_kstar[nPtBins-k]->SetParameter(1, fFit[nPtBins-k]->GetParameter(1));
      fFit_kstar[nPtBins-k]->SetParameter(2, fFit[nPtBins-k]->GetParameter(2));
      fFit_kstar[nPtBins-k]->SetParameter(3, fFit[nPtBins-k]->GetParameter(3));
      
      /* fBg[nPtBins-k]->SetParameter(0, fFit[nPtBins-k]->GetParameter(4)); */
      /* fBg[nPtBins-k]->SetParameter(1, fFit[nPtBins-k]->GetParameter(5)); */
      /* fBg[nPtBins-k]->SetParameter(2, fFit[nPtBins-k]->GetParameter(6)); */
      /* fBg[nPtBins-k]->SetParameter(3, fFit[nPtBins-k]->GetParameter(7)); */

      //===========================
      //===========================================================================================================================
      
      TFitResultPtr ptr = hUSS_Y_Fit[nPtBins-k]->Fit(fFit[nPtBins-k], "SR0", "", Bglow, Bghigh);
      
      double Chi2 = ptr->Chi2();
      int NDF = ptr->Ndf();
      double chi2Ndf = Chi2/NDF;
      string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);
      
      // Draw
      TCanvas* lines = new TCanvas(Form("line_%d",nPtBins-k), "", 940, 800);
      lines->SetLeftMargin(0.12);
      lines->SetRightMargin(0.1);
      lines->SetTopMargin(0.08);
      lines->SetBottomMargin(0.14);

      HistoAxisTitles(hUSS_Y_Fit[nPtBins-k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
		      "1/#it{N}_{ ev} d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.0);
      hUSS_Y_Fit[nPtBins-k]->SetTitle(Form(""));
      
      hUSS_Y_Fit[nPtBins-k]->SetLineColor(1);
      hUSS_Y_Fit[nPtBins-k]->Draw();

      /* fFit[nPtBins-k]->SetLineColor(2); */
      /* fFit[nPtBins-k]->Draw("SAME"); */
      
      fFit_kstar[nPtBins-k]->SetLineColor(6);
      fFit_kstar[nPtBins-k]->Draw("SAME");

      /* fBg[nPtBins-k]->SetLineColor(4); */
      /* fBg[nPtBins-k]->SetLineStyle(2); */
      /* fBg[nPtBins-k]->Draw("SAME"); */

      TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{WIP}", 0.05, 62, 1);
      TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, Kind+" yield", 0.04, 62, 1);
      TLatex* txtInfoCC = DrawTLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
      if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
      if (SpecialLatex_fit2 != nullptr) { SpecialLatex_fit2->Draw("SAME"); }
      TLatex* txtInfoDD = DrawTLatex(0.56, 0.74, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
      TLatex* txtInfoEE = DrawTLatex(0.56, 0.68, chi2NdfText, 0.04, 42, 1);
      
      lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
    }
  
}
