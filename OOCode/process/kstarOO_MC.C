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

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  
  TString localpath = "MC/";

  // Data file inputs 
  TString mainfile = "../kstar-in-oo/Result_rootfile/GenResults.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString USSName("kstar-in-o-o/hMC_USS;1");
  TString LSSName("kstar-in-o-o/hMC_LSS;1");
  TString MixedUSSName("kstar-in-o-o/hMC_USS_Mix;1");
  TString Eventdata("kstar-in-o-o/nEvents_MC;1");
  TString TrueName("kstar-in-o-o/hMC_USS_True;1");
  
  // eventdata counts
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(1); //number 1 is clear number

  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  THnSparseT<TArrayF> *hUSS_MC = (THnSparseT<TArrayF>*)data->Get(USSName);
  THnSparseT<TArrayF> *hLSS_MC = (THnSparseT<TArrayF>*)data->Get(LSSName);
  THnSparseT<TArrayF> *hUSS_Mix_MC = (THnSparseT<TArrayF>*)data->Get(MixedUSSName);
  THnSparseT<TArrayF> *hUSS_True = (THnSparseT<TArrayF>*)data->Get(TrueName);
  
  TH1D* hUSS [nPtBins];
  TH1D* hLSS [nPtBins];
  TH1D* hUSS_Mixed [nPtBins];
  TH1D* hUSS_Mixed_SUB [nPtBins];
  TH1D* hUSSMC_SUB [nPtBins];
  TH1D* hLSSMC_SUB [nPtBins];
  TH1D* hTrue [nPtBins];

  TH1* hFit [nPtBins];
  TH1* hBg [nPtBins];
  TH1* hSig [nPtBins];
  
  Double_t ptlow = 0;
  Double_t pthigh = 0;
  Double_t ReBinNum = 4;
  Double_t nBins = 300.0;
 
  for (Int_t i=0; i<nPtBins; i++)
    {
      hUSSMC_SUB[i] = new TH1D(Form("hUSSMC_SUB_%i", i), Form("hUSSMC_SUB_%i", i), nBins, imlow, imtop);
      hLSSMC_SUB[i] = new TH1D(Form("hLSSMC_SUB_%i", i), Form("hLSSMC_SUB_%i", i), nBins, imlow, imtop);
      hUSS_Mixed_SUB[i] = new TH1D(Form("hUSS_Mixed_SUB_%i", i), Form("hUSS_Mixed_SUB_%i", i), nBins, imlow, imtop);
      hTrue[i] = new TH1D(Form("hTrue_%i", i), Form("hTrue_%i", i), nBins, imlow, imtop);
      
      hBg[i] = new TH1D(Form("hBg_%i", i), Form("hBg_%i", i), nBins, imlow, imtop);
      hSig[i] = new TH1D(Form("hSig_%i", i), Form("hSig_%i", i), nBins, imlow, imtop);
      hFit[i] = new TH1D(Form("hFit_%i", i), Form("hFit_%i", i), nBins, imlow, imtop);
      
      ptlow = xBins[i];
      pthigh = xBins[i+1];
      
      (hUSS_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hLSS_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hUSS_Mix_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hUSS_True->GetAxis(1))->SetRangeUser(ptlow, pthigh);

      hUSS[i] = (TH1D*)hUSS_MC->Projection(2);
      hLSS[i] = (TH1D*)hLSS_MC->Projection(2);
      hUSS_Mixed[i] = (TH1D*)hUSS_Mix_MC->Projection(2);
      hTrue[i] = (TH1D*)hUSS_True->Projection(2);
      
      hUSS[i]->SetName(Form("hUSSMC_%i", i));
      hLSS[i]->SetName(Form("hLSSMC_%i", i));
      hUSS_Mixed[i]->SetName(Form("hUSSMC_Mixed_%i", i));
      hTrue[i]->SetName(Form("hTrue_%i", i));
    }//for1

  for (Int_t j=0; j<nPtBins; j++)
    {
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      Double_t Nlow = 1.1;
      Double_t Nhigh = 1.19;
      
      Double_t Scale = hLSS[j]->Integral(hLSS[j]->FindBin(Nlow), hLSS[j]->FindBin(Nhigh))
	/ hUSS[j]->Integral(hUSS[j]->FindBin(Nlow), hUSS[j]->FindBin(Nhigh));
      
      double delta = 1.006;
      hLSS[j]->Scale(1.0/(Scale*delta));
      
      HistoAxisTitles(hUSS[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS[j]->SetTitle(Form(""));

      hUSS[j]->SetMarkerStyle(4); //circle
      hUSS[j]->SetMarkerSize(1.2);
      hUSS[j]->SetMarkerColor(4);
      hUSS[j]->SetLineColor(4);

      hLSS[j]->SetMarkerStyle(26); //triangle
      hLSS[j]->SetMarkerSize(1.2);
      hLSS[j]->SetMarkerColor(2);
      hLSS[j]->SetLineColor(2);

      hUSS[j]->Rebin(5);
      hLSS[j]->Rebin(5);
      
      hUSS[j]->Scale(0.2); // 1/5 = 0.2
      hLSS[j]->Scale(0.2);

      // 1D projection plots
      TCanvas* plots = new TCanvas(Form("plots_%d", j), "", 940, 800);
      plots->SetLeftMargin(0.16);
      plots->SetRightMargin(0.06);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);

      HistoAxisTitles(hUSS[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS[j]->SetTitle(Form(""));

      hUSS[j]->Draw();
      hLSS[j]->Draw("SAME");

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
      legend->AddEntry(hLSS[j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
      legend->Draw("SAME");
      
      plots->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"/MC_projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");

      // ==========================
      // Drawing Substraction
      hUSSMC_SUB[j] = (TH1D*)hUSS[j]->Clone(Form("hUSS_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSSMC_SUB[j] = (TH1D*)hLSS[j]->Clone(Form("hLSS_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      TCanvas* substraction = new TCanvas(Form("MC_substraction_%d", j), "", 550, 440);
      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.08);
      substraction->SetBottomMargin(0.14);

      hUSSMC_SUB[j]->Add(hLSSMC_SUB[j], -1);
      HistoAxisTitles(hUSSMC_SUB[j], "#it{M}_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]", "Counts");
      hUSSMC_SUB[j]->SetTitle(Form(""));

      TLatex* txtInfoa = DrawTLatex(0.55, 0.88, "O-O yield", 0.04, 62, 1);
      TLatex* txtInfob = DrawTLatex(0.55, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoc = DrawTLatex(0.55, 0.83, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1);
      TLatex* txtInfoe = DrawTLatex(0.55, 0.741, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);

      substraction->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"/MCRec_reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh)); 
    }//for2

  TF1* fkstar [nPtBins];
  TF1* fgaus [nPtBins];
  TF1* ff0 [nPtBins];
  TF1* fpol3 [nPtBins];
  TF1* fpol4 [nPtBins];
  
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];

  TF1* gausbg [nPtBins];
  TF1* pol3bg [nPtBins];
  TF1* f0bg [nPtBins];
  TF1* bw [nPtBins];
  
  double kstarMass = 0.89581;
  double kstarwidth = 0.0464;
  
  double BgMinRange = 0.63;
  double BgMaxRange = 1.20;
  
  for (int k=0; k<nPtBins; k++)
    {
      hSig[k]->Add(hUSSMC_SUB[k]);
      hBg[k]->Add(hUSSMC_SUB[k]);
      hFit[k]->Add(hUSSMC_SUB[k]);

      ptlow = xBins[k];
      pthigh = xBins[k+1];
      // ==============
      // Function
      // ==============
      fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
      fgaus[k] =  new TF1(Form("gaus_%i",k), "gaus(0)", BgMinRange, BgMaxRange);
      fpol3[k] =  new TF1(Form("pol3_%i",k), "pol3(0)", BgMinRange, BgMaxRange);
      fBg[k] = new TF1(Form("Bg_%i",k), "gaus(0)+[3]*BreitWignerRelativistic(x, [4], [5]) + pol3(6)", BgMinRange, BgMaxRange);
      fFit[k] = new TF1(Form("fullFit_%i",k), "gaus(0)+[3]*BreitWignerRelativistic(x, [4], [5]) + pol3(6) + [10]*BreitWignerRelativistic(x, [11], [12])", BgMinRange, BgMaxRange);

      TCanvas* c = new TCanvas(Form("c_%d",k),"", 940, 800);
      //==============================================================================================================
      //playground of fitting parameter

      TH1D* CloneHist = new TH1D("test", "test", nBins/5, imlow, imtop);
      for (int b=0; b<nBins; b++){
	if(b>22 && b<42)
	  continue;
	CloneHist->SetBinContent(b, hBg[k]->GetBinContent(b));
	CloneHist->SetBinError(b, hBg[k]->GetBinError(b));
      }
      CloneHist->Fit(fpol3[k], "SR0", "", 0.60, 1.20);

      fgaus[k]->SetParLimits(1, 0.695, 0.780);
      TH1D* CloneHist3 = new TH1D("test3", "test3", nBins/5, imlow, imtop);
      for (int b=0; b<nBins; b++){
	if(b>26)
	  continue;
	CloneHist3->SetBinContent(b, hBg[k]->GetBinContent(b));
	CloneHist3->SetBinError(b, hBg[k]->GetBinError(b));
      }
      CloneHist3->Fit(fgaus[k], "SR0", "", 0.60, 0.82);
      /* CloneHist3->Draw(); */
      /* fgaus[k]->Draw("SAME"); */
		     
      /* continue; */
      //==============================================================================================================
      //=====
      // 1.
      // Up to this, I find and set the parameter tuned for background
      // Because of the insane residual background shape.
      //====
      
      //======
      // 2.
      // Here, I have the tuning parameter for background, and
      // Now we can merge the part of signal.
      //======

      fBg[k]->SetParameter(0, fgaus[k]->GetParameter(0));
      fBg[k]->SetParameter(1, fgaus[k]->GetParameter(1));
      fBg[k]->SetParameter(2, fgaus[k]->GetParameter(2));
      
      /* fBg[k]->SetParameter(4, 0.968); */
      /* fBg[k]->SetParLimits(5, 0.054, 0.124); */

      fBg[k]->SetParameter(6, fpol3[k]->GetParameter(0));
      fBg[k]->SetParameter(7, fpol3[k]->GetParameter(1));
      fBg[k]->SetParameter(8, fpol3[k]->GetParameter(2));
      fBg[k]->SetParameter(9, fpol3[k]->GetParameter(3));
      
      fkstar[k]->FixParameter(1, 0.892);
      fkstar[k]->FixParameter(2, 0.0464);
      
      TH1D*CloneHist2 = new TH1D("test2", "test2", nBins/5, imlow, imtop);
      for (int b=0; b<nBins; b++){
	if(b<26 || b>36)
	  continue;
	CloneHist2->SetBinContent(b, hBg[k]->GetBinContent(b));
	CloneHist2->SetBinError(b, hBg[k]->GetBinError(b));
      }
      
      CloneHist2->Fit(fkstar[k], "SR0", "", 0.82, 1.01);
      /* CloneHist2->Draw(); */
      /* fkstar[k]->Draw("SAME"); */
      
      fFit[k]->SetParLimits(0, 0, fBg[k]->GetParameter(0));
      fFit[k]->SetParameter(1, fBg[k]->GetParameter(1));
      fFit[k]->SetParameter(2, fBg[k]->GetParameter(2));

      fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));
      fFit[k]->SetParameter(4, fBg[k]->GetParameter(4));
      fFit[k]->SetParameter(5, fBg[k]->GetParameter(5));
      
      fFit[k]->SetParameter(6, fBg[k]->GetParameter(6));
      fFit[k]->SetParameter(7, fBg[k]->GetParameter(7));
      fFit[k]->SetParameter(8, fBg[k]->GetParameter(8));
      fFit[k]->SetParameter(9, fBg[k]->GetParameter(9));
      
      fFit[k]->SetParameter(10, fkstar[k]->GetParameter(0));
      fFit[k]->SetParameter(11, fkstar[k]->GetParameter(1));
      fFit[k]->SetParameter(12, fkstar[k]->GetParameter(2));

      
      hFit[k]->Fit(fFit[k], "SR", "", BgMinRange, BgMaxRange);
      hFit[k]->Fit(fFit[k], "SR", "", BgMinRange, BgMaxRange);
      hFit[k]->Draw();
      hTrue[k]->Draw("SAME");

      gausbg[k] = new TF1(Form("gausbg_%i",k), "gaus(0)", BgMinRange, BgMaxRange);
      f0bg[k] = new TF1(Form("f0bg_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
      pol3bg[k] = new TF1(Form("pol3bg_%i",k), "pol3(0)", BgMinRange, BgMaxRange);
      bw[k] = new TF1(Form("bw_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);

      gausbg[k]->FixParameter(0, fFit[k]->GetParameter(0));
      gausbg[k]->FixParameter(1, fFit[k]->GetParameter(1));
      gausbg[k]->FixParameter(2, fFit[k]->GetParameter(2));

      f0bg[k]->FixParameter(0, fFit[k]->GetParameter(3));
      f0bg[k]->FixParameter(1, fFit[k]->GetParameter(4));
      f0bg[k]->FixParameter(2, fFit[k]->GetParameter(5));

      pol3bg[k]->FixParameter(0, fFit[k]->GetParameter(6));
      pol3bg[k]->FixParameter(1, fFit[k]->GetParameter(7));
      pol3bg[k]->FixParameter(2, fFit[k]->GetParameter(8));
      pol3bg[k]->FixParameter(3, fFit[k]->GetParameter(9));
      
      bw[k]->FixParameter(0, fFit[k]->GetParameter(10));
      bw[k]->FixParameter(1, fFit[k]->GetParameter(11));
      bw[k]->FixParameter(2, fFit[k]->GetParameter(12));

      gausbg[k]->SetLineColor(6); //magenta
      gausbg[k]->Draw("SAME");
      
      f0bg[k]->SetLineColor(4); //blue
      f0bg[k]->Draw("SAME");

      pol3bg[k]->SetLineColor(3); //green
      pol3bg[k]->Draw("SAME");

      bw[k]->SetLineColor(1);
      bw[k]->Draw("SAME");

      fFit[k]->SetLineColor(2);
      fFit[k]->Draw("SAME");

      continue;

   
      /* fBg[k]->FixParameter(8, kstarMass); */
      /* fBg[k]->FixParameter(9, kstarwidth); */
      /* fBg[k]->SetParameter(0,250000); */
      /* fBg[k]->SetParLimits(0,1.2340700e+09/2,1.2340700e+09); */

      /* fBg[k]->SetParameter(0,536538); */
      /* fBg[k]->FixParameter(1,0.632662); */
      /* fBg[k]->FixParameter(2,0.1588); */
      /* hBg[k]->Fit(fBg[k], "SR", "", BgMinRange, BgMaxRange); */
      /* continue; */
      /* hBg[k]->Fit(fBg[k], "SR0", "", BgMinRange, BgMaxRange); */
	
      /* fFit[k] = new TF1(Form("Fullfit_%i", k), "pol3(0) + [4]*BreitWignerRelativistic(x, [5], [6])", BgMinRange, BgMaxRange); */
	
      /* if(k==0){ */
      /* 	hBg[k]->Fit(fBg[k], "SR0", "", BgMinRange, BgMaxRange); */
      /* 	hSig[k]->Fit(fkstar[k], "SR0", "", BgMinRange, BgMaxRange); */
	
      /* 	fFit[k]->SetParameter(0, fBg[k]->GetParameter(0)); */
      /* 	fFit[k]->SetParameter(1, fBg[k]->GetParameter(1)); */
      /* 	fFit[k]->SetParameter(2, fBg[k]->GetParameter(2)); */
      /* 	fFit[k]->SetParameter(3, fBg[k]->GetParameter(3));	 */
      /* 	fFit[k]->SetParameter(4, fkstar[k]->GetParameter(0)); */
      /* 	fFit[k]->FixParameter(5, fkstar[k]->GetParameter(1)); */
      /* 	fFit[k]->FixParameter(6, fkstar[k]->GetParameter(2)); */
      /* } */
      /* else{ */
      /* 	fFit[k]->SetParameter(0, fBg_POST[k-1]->GetParameter(0)); */
      /* 	fFit[k]->SetParameter(1, fBg_POST[k-1]->GetParameter(1)); */
      /* 	fFit[k]->SetParameter(2, fBg_POST[k-1]->GetParameter(2)); */
      /* 	fFit[k]->SetParameter(3, fBg_POST[k-1]->GetParameter(3)); */
      /* 	fFit[k]->SetParameter(4, fkstar_POST[k]->GetParameter(0)); */
      /* 	fFit[k]->FixParameter(5, fkstar_POST[k]->GetParameter(1)); */
      /* 	fFit[k]->FixParameter(6, fkstar_POST[k]->GetParameter(2)); */
      /* } */
      /* hFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange); */
      /* //hFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange); */
      
      /* // ================================ */
      /* // for fitting each one */
      /* // ================================ */
      /* fBg_POST[k] = new TF1(Form("kstar_Bg_POST_%i",k), "pol3(0)", BgMinRange, BgMaxRange); */
      /* fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(0)); */
      /* fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(1)); */
      /* fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(2)); */
      /* fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(3)); */
      /* /\* fBg_POST[k]->SetParameter(4, fFit[k]->GetParameter(4)); *\/ */
      /* /\* fBg_POST[k]->SetParameter(5, fFit[k]->GetParameter(5)); *\/ */
      /* /\* fBg_POST[k]->SetParameter(6, fFit[k]->GetParameter(6)); *\/ */

      /* fkstar_POST[k] = new TF1(Form("kstar_POST_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange); */
      /* fkstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(4)); */
      /* fkstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(5)); */
      /* fkstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(6)); */

      /* TFitResultPtr ptr = hFit[k]->Fit(fFit[k], "SR0", "", BgMinRange, BgMaxRange); */
      /* double Chi2 = ptr->Chi2(); */
      /* int NDF = ptr->Ndf(); */
      /* double chi2Ndf = Chi2/NDF; */
      /* string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf); */

      /* // Draw */
      /* TCanvas* lines = new TCanvas(Form("line_%d",k), "", 940, 800); */
      /* lines->SetLeftMargin(0.14); */
      /* lines->SetRightMargin(0.1); */
      /* lines->SetTopMargin(0.08); */
      /* lines->SetBottomMargin(0.14); */
      /* HistoAxisTitles(hFit[k], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]", */
      /* 		      "d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7); */
      
      /* hFit[k]->SetTitle(Form("")); */
      /* hFit[k]->SetLineColor(1); */
      /* hFit[k]->Draw(); */
     
      /* // total */
      /* fFit[k]->SetLineColor(2); */
      /* fFit[k]->SetLineWidth(4); */
      /* fFit[k]->Draw("SAME"); */
      
      /* // background */
      /* fBg_POST[k]->SetLineColor(4); */
      /* fBg_POST[k]->SetLineWidth(4); */
      /* fBg_POST[k]->Draw("SAME"); */

      /* // kstar */
      /* fkstar_POST[k]->SetLineColor(6); */
      /* fkstar_POST[k]->SetLineWidth(4); */
      /* fkstar_POST[k]->Draw("SAME"); */
      
      /* TLatex* txtInfoAA = DrawTLatex(0.555, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1); */
      /* TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, "O-O yield", 0.04, 62, 1); */
      /* TLatex* txtInfoCC = DrawTLatex(0.56, 0.833, "o-o #sqrt{s_{NN}} = 5.38 TeV", 0.04, 42, 1); */
      /* TLatex* txtInfoDD = DrawTLatex(0.56, 0.789, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1); */
      /* TLatex* txtInfoEE = DrawTLatex(0.56, 0.74, chi2NdfText, 0.04, 42, 1); */

      /* TLegend *legendd = new TLegend(0.55,0.54,0.66,0.64); */
      /* legendd->SetFillColor(4000); */
      /* legendd->SetBorderSize(0); */
      /* legendd->SetTextSize(0.03); */
      /* legendd->AddEntry(hFit[k],"MC", "lep"); */
      /* legendd->AddEntry(fBg_POST[k], "Bkg. MCFit", "l"); */
      /* legendd->AddEntry(fkstar_POST[k], "BreitWignerRelativistic. MCFit", "l"); */
      /* legendd->AddEntry(fFit[k], "Total. Fit", "l"); */
      /* legendd->Draw("SAME"); */
	
      /* lines->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh)); */
      /* // lines->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.pdf", ptlow, pthigh)); */
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
