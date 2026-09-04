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
  const Int_t nPtBins = 11;
  Double_t xBins[nPtBins+1] = {0.0, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0};

  // Centrality Bins setting
  const Int_t nCentBins = 3;
  Double_t CentBins[nCentBins+1] = {0.0, 0.05, 0.1, 0.2};

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  
  TString localpath = "MC/";

  // Data file inputs 
  TString mainfile = "../kstar-in-oo/Result_rootfile/OOMC_v1.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString USSKPiName("kstar-in-o-o/hMC_USS_KPi;1");
  TString USSPiKName("kstar-in-o-o/hMC_USS_PiK;1");
  TString LSSKPiName("kstar-in-o-o/hMC_LSS_KPi;1");
  TString LSSPiKName("kstar-in-o-o/hMC_LSS_PiK;1");
 
  TString Eventdata("kstar-in-o-o/nEvents_MC;1");
  TString TrueName("kstar-in-o-o/hMC_USS_True;1");
  
  // eventdata counts
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(1); //number 1 is clear number

  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  THnSparseT<TArrayF> *hUSS_KPi_MC = (THnSparseT<TArrayF>*)data->Get(USSKPiName);
  THnSparseT<TArrayF> *hUSS_PiK_MC = (THnSparseT<TArrayF>*)data->Get(USSPiKName);
  THnSparseT<TArrayF> *hLSS_KPi_MC = (THnSparseT<TArrayF>*)data->Get(LSSKPiName);
  THnSparseT<TArrayF> *hLSS_PiK_MC = (THnSparseT<TArrayF>*)data->Get(LSSPiKName);
  
  THnSparseT<TArrayF> *hUSS_True = (THnSparseT<TArrayF>*)data->Get(TrueName);

  TH1D* hUSSKPi [nPtBins][nCentBins];
  TH1D* hUSSPiK [nPtBins][nCentBins];
  TH1D* hUSS [nPtBins][nCentBins];

  TH1D* hLSSKPi [nPtBins][nCentBins];
  TH1D* hLSSPiK [nPtBins][nCentBins];
  TH1D* hLSS [nPtBins][nCentBins];
  
  TH1D* hLSSMC_SUB [nPtBins][nCentBins];
  TH1D* hUSSMC_SUB [nPtBins][nCentBins];
  TH1D* hUSSMC_True [nPtBins][nCentBins];

  TH1D* hCent [nPtBins][nCentBins];

  TH1* hUSSFit [nPtBins][nCentBins];
  TH1* hUSSBg [nPtBins][nCentBins];
  TH1* hUSSkstarSig [nPtBins][nCentBins];
  
  double ptlow = 0;
  double pthigh = 0;
  double centlow = 0;
  double centhigh = 0;
  
  Double_t ReBinNum = 5;
  Double_t nBins = 300.0;
 
  for (Int_t i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];

      for (Int_t j=0; j<nCentBins; j++)
	{
	  centlow = CentBins[j];
	  centhigh = CentBins[j+1];
	  
	  hUSSMC_SUB[i][j] = new TH1D(Form("hUSSMC_SUB_%i_%i", i,j), Form("hUSSMC_SUB_%i_%i", i,j), nBins, imlow, imtop);
	  hLSSMC_SUB[i][j] = new TH1D(Form("hLSSMC_SUB_%i_%i", i,j), Form("hLSSMC_SUB_%i_%i", i,j), nBins, imlow, imtop);
	  hUSSMC_True[i][j] = new TH1D(Form("hUSSMC_True_%i_%i", i,j), Form("hUSSMC_True_%i_%i", i,j), nBins, imlow, imtop);
	  hUSSFit[i][j] = new TH1D(Form("hUSSFit_%i_%i", i,j), Form("hUSSFit_%i_%i", i,j), nBins/5, imlow, imtop); 

	  
	  (hUSS_KPi_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hUSS_PiK_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hUSS_KPi_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  (hUSS_PiK_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  hUSSKPi[i][j] = (TH1D*)hUSS_KPi_MC->Projection(2);
	  hUSSPiK[i][j] = (TH1D*)hUSS_PiK_MC->Projection(2);
      
	  hUSS[i][j] = (TH1D*)hUSSKPi[i][j]->Clone(Form("hUSS_KPi%i_%i", i,j));
	  hUSS[i][j]->Add(hUSSPiK[i][j]);
	  hUSS[i][j]->Scale(0.5);
  
	  (hLSS_KPi_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hLSS_PiK_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hLSS_KPi_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  (hLSS_PiK_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  
	  hLSSKPi[i][j] = (TH1D*)hLSS_KPi_MC->Projection(2);
	  hLSSPiK[i][j] = (TH1D*)hLSS_PiK_MC->Projection(2);

	  hLSS[i][j] = (TH1D*)hLSSKPi[i][j]->Clone(Form("hLSS_KPi_%i_%i", i,j));
	  hLSS[i][j]->Add(hLSSPiK[i][j]);
	  hLSS[i][j]->Scale(0.5);

	  // True
	  (hUSS_True->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hUSS_True->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  hUSSMC_True[i][j] = (TH1D*)hUSS_True->Projection(2);
      
	  hUSS[i][j]->SetName(Form("hUSSMC_%i_%i", i,j));
	  hLSS[i][j]->SetName(Form("hLSSMC_%i_%i", i,j));
	  hUSSMC_True[i][j]->SetName(Form("hUSSMC_True_%i_%i", i,j));
	  
	  hCent[i][j] = (TH1D*)hUSS_KPi_MC->Projection(0);
	  hCent[i][j]->SetName(Form("hUSSMC_Cent_%i_%i", i,j));

	  delete hUSSKPi[i][j];
	  delete hUSSPiK[i][j];
	  delete hLSSKPi[i][j];
	  delete hLSSPiK[i][j];
      
	}//for1
    }//for2
  for (Int_t i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];
      
      for (Int_t j=0; j<nCentBins; j++)
	{
	  centlow = CentBins[j];
	  centhigh = CentBins[j+1];

	  Double_t Nlow = 1.12;
	  Double_t Nhigh = 1.199;
	  Double_t Scale = hLSS[i][j]->Integral(hLSS[i][j]->FindBin(Nlow), hLSS[i][j]->FindBin(Nhigh))
	    / hUSS[i][j]->Integral(hUSS[i][j]->FindBin(Nlow), hUSS[i][j]->FindBin(Nhigh));
      
	  double delta;
	  delta = 0.986;
	  hLSS[i][j]->Scale(1.0/Scale*delta);
	  
	  HistoAxisTitles(hUSS[i][j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
	  hUSS[i][j]->SetTitle(Form(""));

	  hUSS[i][j]->SetMarkerStyle(4); //circle
	  hUSS[i][j]->SetMarkerSize(1.2);
	  hUSS[i][j]->SetMarkerColor(4);
	  hUSS[i][j]->SetLineColor(4);

	  hLSS[i][j]->SetMarkerStyle(26); //triangle
	  hLSS[i][j]->SetMarkerSize(1.2);
	  hLSS[i][j]->SetMarkerColor(2);
	  hLSS[i][j]->SetLineColor(2);

	  hUSS[i][j]->Rebin(5);
	  hLSS[i][j]->Rebin(5);
	  hUSSMC_True[i][j]->Rebin(5);
      
	  hUSS[i][j]->Scale(0.2);
	  hLSS[i][j]->Scale(0.2);
	  hUSSMC_True[i][j]->Scale(0.2);

	  // 1D projection plots
	  TCanvas* plots = new TCanvas(Form("plots_%d_%d", i,j), "", 940, 800);
	  plots->SetLeftMargin(0.16);
	  plots->SetRightMargin(0.06);
	  plots->SetTopMargin(0.08);
	  plots->SetBottomMargin(0.14);

	  hUSS[i][j]->Draw();
	  hLSS[i][j]->Draw("SAME");
	  
	  TH1D* hShade = (TH1D*)hLSS[i][j]->Clone(Form("hShade_%d_%d", i, j));
	  int binLow = hShade->FindBin(Nlow);
	  int binHigh = hShade->FindBin(Nhigh);
      
	  for (int bin=1; bin <= hShade->GetNbinsX(); ++bin){
	    if (bin < binLow || bin > binHigh) {
	      hShade->SetBinContent(bin, 0);
	    }
	  }
	  hShade->SetFillStyle(3003);
	  hShade->SetFillColor(kBlue);
	  hShade->Draw("SAME HIST F");

	  TLatex* txtInfoA = DrawTLatex(0.526, 0.392, "O-O yield", 0.04, 62, 1);
	  TLatex* txtInfoB = DrawTLatex(0.526, 0.445, "ALICE #bf{Performance}", 0.05, 62, 1);
	  TLatex* txtInfoC = DrawTLatex(0.526, 0.341, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
	  TLatex* txtInfoD = DrawTLatex(0.526, 0.286, "Min. Bias", 0.04, 42, 1);
	  TLatex* txtInfoF = DrawTLatex(0.526, 0.241, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  TLatex* txtInfoG = DrawTLatex(0.526, 0.187, Form("Centrality %d-%d %%", (int)(centlow*100), (int)(centhigh*100)), 0.04, 42, 1);
 
	  TLegend* legend = new TLegend(0.265,0.167,0.498,0.285); //LEFT x, y, Right x, y
	  legend->SetFillColor(4000);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.04);
	  legend->AddEntry(hUSS[i][j], "K^{+}#pi^{-}", "lpf");
	  legend->AddEntry(hLSS[i][j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
	  legend->Draw("SAME");

	  plots->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"proj_%1.1f_%1.1f_cent_%1.1f_%1.1f.png", ptlow, pthigh, centlow, centhigh),"RECREATE");
  
	  // ==========================
	  // Drawing Substraction
	  hUSSMC_SUB[i][j] = (TH1D*)hUSS[i][j]->Clone(Form("hUSS_Minv_%1.1f_%1.1f_cent_%1.1f_%1.1f", ptlow, pthigh, centlow, centhigh));
	  hLSSMC_SUB[i][j] = (TH1D*)hLSS[i][j]->Clone(Form("hLSS_Minv_%1.1f_%1.1f_cent_%1.1f_%1.1f", ptlow, pthigh, centlow, centhigh));
	  hUSSMC_SUB[i][j]->Add(hLSSMC_SUB[i][j], -1);
	}//for2
    }//for1

  TF1* fkstar [nPtBins][nCentBins];
  TF1* fBg [nPtBins][nCentBins];
  TF1* fBg2 [nPtBins][nCentBins];
  TF1* fFit [nPtBins][nCentBins];
  TF1* fBg_POST [nPtBins][nCentBins];
  TF1* fkstar_POST [nPtBins][nCentBins];
  TF1* gausbg [nPtBins][nCentBins];
  TF1* f0bg [nPtBins][nCentBins];
  TF1* pol2bg [nPtBins][nCentBins];
  TF1* bw [nPtBins][nCentBins];
  TF1* totbg [nPtBins][nCentBins];

  double kstarMass = 0.89581;
  double kstarwidth = 0.0464;

  double BgMinRange = 0.69;
  double BgMaxRange = 1.18;

  TH1D* hYield[3];
  hYield[0] = new TH1D("binYield_0_5","binYield;pT;N", nPtBins,xBins);
  hYield[1] = new TH1D("binYield_5_10","binYield;pT;N", nPtBins,xBins);
  hYield[2] = new TH1D("binYield_10_20","binYield;pT;N", nPtBins,xBins);

  int red0 = TColor::GetColor("#FF0000");     // pink and red
  int violet1 = TColor::GetColor("#BF00FF");    // hot pink
  int blue3 = TColor::GetColor("#0096FF");  // orange
  
  for (int i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];
    
      for (int j=0; j<nCentBins; j++)
	{
	  centlow = CentBins[j];
	  centhigh = CentBins[j+1];
	  
	  hUSSFit[i][j]->Add(hUSSMC_SUB[i][j]);

	  // ==============
	  // Function
	  // ==============
	  fkstar[i][j] = new TF1(Form("kstar_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", 0.84, 0.95);
	  fBg2[i][j] = new TF1(Form("kstar_GausBg_%i_cent_%i",i,j), "gaus(0)+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6)", BgMinRange, BgMaxRange);
	  fFit[i][j] = new TF1(Form("kstar_fullFit_%i_cent_%i",i,j), "gaus(0)+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6) + [9]*BreitWignerRelativistic(x, [10], [11])", BgMinRange, BgMaxRange);


	  TCanvas* lines = new TCanvas(Form("line_%d_cent_%d",i,j), "", 940, 800);
     
	  fkstar[i][j]->FixParameter(1, kstarMass); 
	  fkstar[i][j]->FixParameter(2, kstarwidth);

	  //=================================================================
	  fFit[i][j]->SetParLimits(0, 0, 40000);
	  fFit[i][j]->SetParLimits(1, 0.725, 0.82);
	  fFit[i][j]->SetParLimits(2, 0.0363, 0.126);

	  fFit[i][j]->SetParLimits(4, 0.940, 0.98);
	  fFit[i][j]->SetParLimits(5, 0.09, 0.187);
	  
	  fFit[i][j]->SetParameter(6, fBg2[i][j]->GetParameter(6));
	  fFit[i][j]->SetParameter(7, fBg2[i][j]->GetParameter(7));
	  fFit[i][j]->SetParameter(8, fBg2[i][j]->GetParameter(8));
	  fFit[i][j]->SetParameter(9, fkstar[i][j]->GetParameter(0)); //DO NOT FIX      
	  fFit[i][j]->FixParameter(10, fkstar[i][j]->GetParameter(1));
	  fFit[i][j]->FixParameter(11, fkstar[i][j]->GetParameter(2));

	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->SetTitle(Form(""));
	  hUSSFit[i][j]->Draw();
	  hUSSMC_True[i][j]->Draw("SAME");
      
	  gausbg[i][j] = new TF1(Form("gausbg_%i_cent_%i",i,j), "gaus(0)", BgMinRange, BgMaxRange);
	  f0bg[i][j] = new TF1(Form("f0bg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
	  pol2bg[i][j] = new TF1(Form("pol2bg_%i_cent_%i",i,j), "pol2(0)", BgMinRange, BgMaxRange);
	  bw[i][j] = new TF1(Form("bw_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
	  totbg[i][j] = new TF1(Form("totbg_%i_cent_%i",i,j), "gaus(0) + [3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6)", BgMinRange, BgMaxRange);

	  gausbg[i][j]->FixParameter(0, fFit[i][j]->GetParameter(0));
	  gausbg[i][j]->FixParameter(1, fFit[i][j]->GetParameter(1));
	  gausbg[i][j]->FixParameter(2, fFit[i][j]->GetParameter(2));

	  f0bg[i][j]->FixParameter(0, fFit[i][j]->GetParameter(3));
	  f0bg[i][j]->FixParameter(1, fFit[i][j]->GetParameter(4));
	  f0bg[i][j]->FixParameter(2, fFit[i][j]->GetParameter(5));

	  pol2bg[i][j]->FixParameter(0, fFit[i][j]->GetParameter(6));
	  pol2bg[i][j]->FixParameter(1, fFit[i][j]->GetParameter(7));
	  pol2bg[i][j]->FixParameter(2, fFit[i][j]->GetParameter(8));
   
	  bw[i][j]->FixParameter(0, fFit[i][j]->GetParameter(9));
	  bw[i][j]->FixParameter(1, fFit[i][j]->GetParameter(10));
	  bw[i][j]->FixParameter(2, fFit[i][j]->GetParameter(11));

	  totbg[i][j]->FixParameter(0, fFit[i][j]->GetParameter(0));
	  totbg[i][j]->FixParameter(1, fFit[i][j]->GetParameter(1));
	  totbg[i][j]->FixParameter(2, fFit[i][j]->GetParameter(2));
	  totbg[i][j]->FixParameter(3, fFit[i][j]->GetParameter(3));
	  totbg[i][j]->FixParameter(4, fFit[i][j]->GetParameter(4));
	  totbg[i][j]->FixParameter(5, fFit[i][j]->GetParameter(5));
	  totbg[i][j]->FixParameter(6, fFit[i][j]->GetParameter(6));
	  totbg[i][j]->FixParameter(7, fFit[i][j]->GetParameter(7));
	  totbg[i][j]->FixParameter(8, fFit[i][j]->GetParameter(8));

	  TFitResultPtr ptr = hUSSFit[i][j]->Fit(fFit[i][j], "SR0", "", BgMinRange, BgMaxRange);
	  double Chi2 = ptr->Chi2();
	  int NDF = ptr->Ndf();
	  double chi2Ndf = Chi2/NDF;
	  string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);

	  // Draw
	  lines->SetLeftMargin(0.14);
	  lines->SetRightMargin(0.1);
	  lines->SetTopMargin(0.08);
	  lines->SetBottomMargin(0.14);
	  HistoAxisTitles(hUSSFit[i][j], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
			  "d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);

	  gausbg[i][j]->SetLineColor(6); //magenta
	  gausbg[i][j]->Draw("SAME");
      
	  f0bg[i][j]->SetLineColor(4); //blue
	  f0bg[i][j]->Draw("SAME");

	  pol2bg[i][j]->SetLineColor(3); //green
	  pol2bg[i][j]->Draw("SAME");

	  bw[i][j]->SetLineColor(1);
	  bw[i][j]->Draw("SAME");

	  totbg[i][j]->SetLineColor(7);
	  totbg[i][j]->Draw("SAME");

	  fFit[i][j]->SetLineColor(2);
	  fFit[i][j]->Draw("SAME");

	  TLatex* txtInfoAA = DrawTLatex(0.167, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
	  TLatex* txtInfoBB = DrawTLatex(0.173, 0.873, "O-O yield", 0.04, 62, 1);
	  TLatex* txtInfoCC = DrawTLatex(0.173, 0.833, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
	  TLatex* txtInfoDD = DrawTLatex(0.173, 0.789, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  TLatex* txtInfoEE = DrawTLatex(0.173, 0.74, Form("Centrality %d-%d %%", (int)(centlow*100), (int)(centhigh*100)), 0.04, 42, 1);
	  TLatex* txtInfoFF = DrawTLatex(0.173, 0.70, chi2NdfText, 0.04, 42, 1);

	  TLegend *legendd = new TLegend(0.77,0.74,0.89,0.90);
	  legendd->SetFillColor(4000);
	  legendd->SetBorderSize(0);
	  legendd->SetTextSize(0.03);

	  legendd->AddEntry(totbg[i][j], "Tot Bkg", "l");
	  legendd->AddEntry(bw[i][j], "kstar", "l");
	  legendd->AddEntry(gausbg[i][j], "gaus", "l");
	  legendd->AddEntry(f0bg[i][j], "f0bw", "l");
	  legendd->AddEntry(pol2bg[i][j], "pol2", "l");
	  legendd->Draw("SAME");
      
	  lines->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"two_fit_%1.1f_%1.1f_cent_%1.1f_%1.1f.png", ptlow, pthigh, centlow, centhigh), "RECREATE");

	  //=============================================
	  // Yield Main
	  //=============================================
	  double HistY = 0;
	  double HistYerr = 0;
	  double lowM = hUSSFit[i][j]->FindBin(0.85);
	  double highM = hUSSFit[i][j]->FindBin(0.95);
	  HistY=hUSSFit[i][j]->IntegralAndError(lowM, highM, HistYerr);
	  double binW = hUSSFit[i][j]->GetBinWidth(50); //Any number
	  double BgY = totbg[i][j]->Integral(0.85,0.95)/binW;
	
	  TMatrixDSym covTot(fFit[i][j]->GetNpar());
	  TMatrixDSym covGG(bw[i][j]->GetNpar());
	  TMatrixDSym covGG_err(totbg[i][j]->GetNpar());
	  covTot = ptr->GetCovarianceMatrix();
	  covGG_err = covTot.GetSub(0,8,0,8);// these have to be totbg parameters
	  
	  Double_t Bg_err = (1/binW)*(totbg[i][j]->IntegralError(0.85, 0.95, totbg[i][j]->GetParameters(),covGG_err.GetMatrixArray()));

	  double HistCorrY = HistY-BgY;
	  double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	  hYield[j]->SetBinContent(i+1, HistCorrY);
	  hYield[j]->SetBinContent(i+1, hYield[j]->GetBinContent(i+1)/(hYield[j]->GetBinWidth(i+1)));
	  hYield[j]->SetBinError(i+1, HistCorrYerr);

	  double CLow = hCent[i][j]->FindBin(CentBins[j]);
	  double CHigh = hCent[i][j]->FindBin(CentBins[j+1]);
	  double nEvents = 0.0;
	  double totEvents = hCent[i][j]->Integral(1, hCent[i][j]->FindBin(0.2));
	  if (j==0)
	    {
	      nEvents = hCent[i][j]->Integral(1, hCent[i][j]->FindBin(0.05));
	    }
	  else if (j!=0){
	    nEvents = hCent[i][j]->Integral(CLow, CHigh);
	  }
	  // hYield[j]->Scale(1.0 / nEvents);
	  
	  Shreck();
	  cout<<"HistCorrY: "<<HistCorrY<<endl;
	  cout<<"Bg_err: "<<Bg_err<<endl;
	  cout<<"nEvents: "<<nEvents<<endl;
	
	  //Extra work
	  //0-5%
	  if(j==0){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(blue3);
	    hYield[j]->SetMarkerStyle(43); //circle
	    hYield[j]->SetMarkerSize(3);
	    hYield[j]->SetMarkerColor(blue3);
	  }
	  //5-10%
	  if(j==1){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(violet1);
	    hYield[j]->SetMarkerStyle(93); //hospital
	    hYield[j]->SetMarkerSize(3);
	    hYield[j]->SetMarkerColor(violet1);
	  }

	  //10-20%
	  if(j==2){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(red0);
	    hYield[j]->SetMarkerStyle(20); //circle
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(red0);
	  }
	  cout << "//==========================================================//" << endl;
	  cout << "Total events in all centrality bins: " << totEvents << endl;
	} //for j
    }//for i
  
  TCanvas* numbin = new TCanvas("","",1200,800);
  numbin->SetRightMargin(0.05);
  numbin->SetLeftMargin(0.17);
  numbin->SetBottomMargin(0.12);

  double yMin = hYield[0]->GetMinimum();
  double yMax = hYield[0]->GetMaximum();
	
  TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
  hdummy->SetMinimum(yMin*0.5);
  hdummy->SetMaximum(yMax*1.5);


  HistoAxisTitles(hdummy, "#it{p}_{T} [G.qeV/#it{c}]", "#frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}",
		  0.05, 62, 1.0, 0.05, 62, 1.4);
  hdummy->SetLineWidth(2);
  hdummy->SetLineColor(1);
  hdummy->Draw();

  for(int j=0; j<3; j++)
    hYield[j]->Draw("SAME");

  TLatex* txtInfoAAA = DrawTLatex(0.66, 0.845, "ALICE #bf{Performance}", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.66, 0.789, "O-O #sqrt{s_{NN}} = 5.38 TeV", 0.05, 62, 1);
  TLegend* legenddd = new TLegend(0.64, 0.54, 0.97, 0.74);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);

  legenddd->AddEntry(hYield[2], "Centrality 10-20%", "pl");
  legenddd->AddEntry(hYield[1], "Centrality 5-10%", "pl");
  legenddd->AddEntry(hYield[0], "Centrality 0-5%", "pl"); //pl : point line

  legenddd->Draw("SAME");

  numbin->SetLogy();
  numbin->SaveAs("../kstar-in-oo/Plot/"+localpath+"h_Yield.png", "RECREATE");
      
  
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
