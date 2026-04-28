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
} // Main Drawfunction


void DrawBW()
{
  cout.precision(6);
  gStyle->SetOptStat(0);
  
  // pTBins setting 
  const Int_t nPtBins = 9;
  Double_t xBins[nPtBins+1] = {0.0, 1.2, 1.8, 2.4, 3.2, 4.0, 5.0, 6.0, 7.0, 8.0};

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;

  TString localpath = "Data/Inside/JetHighPt/";

  // Data file inputs 
  TString mainfile = "../kstar-in-jets/Results_rootfile/MB_Jets_pass1_ao_skimmed.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString USSKPiName("kstar-in-o-o_id47421/hUSS_KPi_INSIDE;1");
  TString USSPiKName("kstar-in-o-o_id47421/hUSS_PiK_INSIDE;1");
  TString LSSKPiName("kstar-in-o-o_id47421/hLSS_KPi_INSIDE;1");
  TString LSSPiKName("kstar-in-o-o_id47421/hLSS_PiK_INSIDE;1");
  TString CentName("kstar-in-o-o_id47421/nEvents;1");
  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  THnSparseT<TArrayF> *hUSS_KPi = (THnSparseT<TArrayF>*)data->Get(USSKPiName);
  THnSparseT<TArrayF> *hUSS_PiK = (THnSparseT<TArrayF>*)data->Get(USSPiKName);
  THnSparseT<TArrayF> *hLSS_KPi = (THnSparseT<TArrayF>*)data->Get(LSSKPiName);
  THnSparseT<TArrayF> *hLSS_PiK = (THnSparseT<TArrayF>*)data->Get(LSSPiKName);

  TH1D* hUSSKPi [nPtBins];
  TH1D* hUSSPiK [nPtBins];
  TH1D* hUSS [nPtBins];

  TH1D* hLSSKPi [nPtBins];
  TH1D* hLSSPiK [nPtBins];
  TH1D* hLSS [nPtBins];
  
  TH1D* hLSS_SUB [nPtBins];
  TH1D* hUSS_SUB [nPtBins];

  TH1D* hCent;
  hCent = (TH1D*)data->Get(CentName);

  TH1* hUSSFit [nPtBins];
  TH1* hUSSBg [nPtBins];
  TH1* hUSSkstarSig [nPtBins];
  
  double ptlow = 0;
  double pthigh = 0;

  Double_t ReBinNum = 5;
  Double_t nBins = 300.0;
 
  for (Int_t i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];

      hUSS_SUB[i] = new TH1D(Form("hUSS_SUB_%i", i), Form("hUSS_SUB_%i", i), nBins, imlow, imtop);
      hLSS_SUB[i] = new TH1D(Form("hLSS_SUB_%i", i), Form("hLSS_SUB_%i", i), nBins, imlow, imtop);
      hUSSFit[i] = new TH1D(Form("hUSSFit_%i", i), Form("hUSSFit_%i", i), nBins/5, imlow, imtop); 
	  
      (hUSS_KPi->GetAxis(2))->SetRangeUser(ptlow, pthigh);
      (hUSS_PiK->GetAxis(2))->SetRangeUser(ptlow, pthigh);

      hUSSKPi[i] = (TH1D*)hUSS_KPi->Projection(3); // Mass
      hUSSPiK[i] = (TH1D*)hUSS_PiK->Projection(3);
      
      hUSS[i] = (TH1D*)hUSSKPi[i]->Clone(Form("hUSS_KPi%i", i));
      hUSS[i]->Add(hUSSPiK[i]); //KPi+PiK
      hUSS[i]->Scale(0.5);

      (hLSS_KPi->GetAxis(2))->SetRangeUser(ptlow, pthigh);
      (hLSS_PiK->GetAxis(2))->SetRangeUser(ptlow, pthigh);

      hLSSKPi[i] = (TH1D*)hLSS_KPi->Projection(3);
      hLSSPiK[i] = (TH1D*)hLSS_PiK->Projection(3);

      hLSS[i] = (TH1D*)hLSSKPi[i]->Clone(Form("hLSS_KPi_%i", i));
      hLSS[i]->Add(hLSSPiK[i]);
      hLSS[i]->Scale(0.5);

      hUSS[i]->SetName(Form("hUSS_%i", i));
      hLSS[i]->SetName(Form("hLSS_%i", i));

      delete hUSSKPi[i];
      delete hUSSPiK[i];
      delete hLSSKPi[i];
      delete hLSSPiK[i];
    }//for

  for (Int_t i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];

      Double_t Nlow = 1.10;
      Double_t Nhigh = 1.1999;
      Double_t Scale = hLSS[i]->Integral(hLSS[i]->FindBin(Nlow), hLSS[i]->FindBin(Nhigh))
	/ hUSS[i]->Integral(hUSS[i]->FindBin(Nlow), hUSS[i]->FindBin(Nhigh));

      double delta;
      delta = 1.0;
      hLSS[i]->Scale(1.0/Scale*delta);

      HistoAxisTitles(hUSS[i], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS[i]->SetTitle(Form(""));

      hUSS[i]->SetMarkerStyle(4); //circle
      hUSS[i]->SetMarkerSize(1.2);
      hUSS[i]->SetMarkerColor(4);
      hUSS[i]->SetLineColor(4);

      hLSS[i]->SetMarkerStyle(26); //triangle
      hLSS[i]->SetMarkerSize(1.2);
      hLSS[i]->SetMarkerColor(2);
      hLSS[i]->SetLineColor(2);

      hUSS[i]->Rebin(5);
      hLSS[i]->Rebin(5);
      
      hUSS[i]->Scale(0.2);
      hLSS[i]->Scale(0.2);

      // 1D projection plots
      TCanvas* plots = new TCanvas(Form("plots_%d", i), "", 940, 800);
      plots->SetLeftMargin(0.16);
      plots->SetRightMargin(0.06);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);

      hUSS[i]->Draw();
      hLSS[i]->Draw("SAME");
      TH1D* hShade = (TH1D*)hLSS[i]->Clone(Form("hShade_%d", i));
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

      TLatex* txtInfoA = DrawTLatex(0.521, 0.318, "Jets", 0.04, 62, 1);
      TLatex* txtInfoB = DrawTLatex(0.339, 0.365, "ALICE Work In Progress", 0.05, 62, 1);
      TLatex* txtInfoC = DrawTLatex(0.522, 0.265, "pp, #sqrt{s_{NN}} = 13.6 TeV", 0.04, 42, 1);
      TLatex* txtInfoD = DrawTLatex(0.522, 0.216, "Jet Ch High pT", 0.04, 42, 1);
      TLatex* txtInfoF = DrawTLatex(0.522, 0.169, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
 
      TLegend* legend = new TLegend(0.265,0.167,0.498,0.285); //LEFT x, y, Right x, y
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.04);
      legend->AddEntry(hUSS[i], "K^{+}#pi^{-}", "lpf");
      legend->AddEntry(hLSS[i], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
      legend->Draw("SAME");

      plots->SaveAs(Form("../kstar-in-jets/Plot/kstarInjets/"+localpath+"proj_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");
  
      // ==========================
      // Drawing Substraction
      hUSS_SUB[i] = (TH1D*)hUSS[i]->Clone(Form("hUSS_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_SUB[i] = (TH1D*)hLSS[i]->Clone(Form("hLSS_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hUSS_SUB[i]->Add(hLSS_SUB[i], -1);
    }//for

  TF1* fkstar [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  TF1* pol3bg [nPtBins];
  TF1* bw [nPtBins];

  double kstarMass;
  double kstarwidth;

  double BgMinRange = 0.70;
  double BgMaxRange = 1.16;
  
  TH1D* hYield = new TH1D("binYield","binYield;pT;N", nPtBins,xBins);
  int red = TColor::GetColor("#ff0400");
  
  double nEvents = 0.0;
  
  for (int i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];

      hUSSFit[i]->Add(hUSS_SUB[i]);
      // ==============
      // Function
      // ==============
      fkstar[i] = new TF1(Form("kstar_%i",i), "[0]*BreitWignerRelativistic(x, [1], [2])", 0.84, 0.95);
      fBg[i] = new TF1(Form("Bg_%i",i), "pol3(0)", BgMinRange, BgMaxRange);
      fFit[i] = new TF1(Form("kstar_fullFit_%i",i), "[0]*BreitWignerRelativistic(x, [1], [2]) + pol3(3)", BgMinRange, BgMaxRange);


      TCanvas* lines = new TCanvas(Form("line_%d",i), "", 940, 800);

      fkstar[i]->SetParameter(1, 0.890);
      fkstar[i]->SetParameter(2, 0.0524);
      //=================================================================
      
      fFit[i]->SetParameter(0, fkstar[i]->GetParameter(0)); //DO NOT FIX
      fFit[i]->SetParameter(1, fkstar[i]->GetParameter(1));
      fFit[i]->FixParameter(2, fkstar[i]->GetParameter(2));

      fFit[i]->SetParameter(3, fBg[i]->GetParameter(0));
      fFit[i]->SetParameter(4, fBg[i]->GetParameter(1));
      fFit[i]->SetParameter(5, fBg[i]->GetParameter(2));
      fFit[i]->SetParameter(6, fBg[i]->GetParameter(3));
  
      hUSSFit[i]->Fit(fFit[i], "SR", "", BgMinRange, BgMaxRange);
      hUSSFit[i]->Fit(fFit[i], "SR", "", BgMinRange, BgMaxRange);
      hUSSFit[i]->SetTitle(Form(""));
      hUSSFit[i]->Draw();

      pol3bg[i] = new TF1(Form("pol3bg_%i",i), "pol3(0)", BgMinRange, BgMaxRange);
      bw[i] = new TF1(Form("bw_%i",i), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);

      bw[i]->FixParameter(0, fFit[i]->GetParameter(0));
      bw[i]->FixParameter(1, fFit[i]->GetParameter(1));
      bw[i]->FixParameter(2, fFit[i]->GetParameter(2));
      
      pol3bg[i]->FixParameter(0, fFit[i]->GetParameter(3));
      pol3bg[i]->FixParameter(1, fFit[i]->GetParameter(4));
      pol3bg[i]->FixParameter(2, fFit[i]->GetParameter(5));
      pol3bg[i]->FixParameter(3, fFit[i]->GetParameter(6));

      TFitResultPtr ptr = hUSSFit[i]->Fit(fFit[i], "SR0", "", BgMinRange, BgMaxRange);
      double Chi2 = ptr->Chi2();
      int NDF = ptr->Ndf();
      double chi2Ndf = Chi2/NDF;
      string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);

      // Draw
      lines->SetLeftMargin(0.14);
      lines->SetRightMargin(0.1);
      lines->SetTopMargin(0.08);
      lines->SetBottomMargin(0.14);
      HistoAxisTitles(hUSSFit[i], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
		      "d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);

      pol3bg[i]->SetLineColor(3); //green
      pol3bg[i]->Draw("SAME");

      bw[i]->SetLineColor(1);
      bw[i]->Draw("SAME");

      fFit[i]->SetLineColor(2);
      fFit[i]->Draw("SAME");

      TLatex* txtInfoAA = DrawTLatex(0.167, 0.93, "ALICE Work In Progress", 0.05, 62, 1);
      TLatex* txtInfoBB = DrawTLatex(0.173, 0.873, "JetChHighPt", 0.04, 62, 1);
      TLatex* txtInfoCC = DrawTLatex(0.173, 0.820, "pp #sqrt{s_{NN}} = 13.6 TeV", 0.04, 42, 1);
      TLatex* txtInfoDD = DrawTLatex(0.173, 0.767, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
      TLatex* txtInfoFF = DrawTLatex(0.173, 0.714, chi2NdfText, 0.04, 42, 1);

      TLegend *legendd = new TLegend(0.717,0.824,0.837,0.903);
      legendd->SetFillColor(4000);
      legendd->SetBorderSize(0);
      legendd->SetTextSize(0.03);

      legendd->AddEntry(bw[i], "K^{*0}", "l");
      legendd->AddEntry(pol3bg[i], "Residual bkg.", "l");
      legendd->Draw("SAME");
      
      lines->SaveAs(Form("../kstar-in-jets/Plot/kstarInjets/"+localpath+"two_fit_%1.1f_%1.1f.png", ptlow, pthigh), "RECREATE");
      
      //=============================================
      // Yield Main
      //=============================================
      // Recon Yield
      double HistY = 0;
      double HistYerr = 0;
      double lowM = hUSSFit[i]->FindBin(BgMinRange);
      double highM = hUSSFit[i]->FindBin(BgMaxRange);
      HistY=hUSSFit[i]->IntegralAndError(lowM, highM, HistYerr);
      double binW = hUSSFit[i]->GetBinWidth(50); //Any bin number
      double BgY = (1/binW)*pol3bg[i]->Integral(BgMinRange, BgMaxRange);
	
      TMatrixDSym covTot(fFit[i]->GetNpar());
      TMatrixDSym covGG(bw[i]->GetNpar());
      TMatrixDSym covGG_err(pol3bg[i]->GetNpar());
      covTot = ptr->GetCovarianceMatrix();
      covGG_err = covTot.GetSub(0, 3, 0, 3);// these have to be pol3bg parameters

      Double_t Bg_err = (1/binW)*pol3bg[i]->IntegralError(BgMinRange, BgMaxRange, pol3bg[i]->GetParameters(), covGG_err.GetMatrixArray());
	 
      double HistCorrY = HistY-BgY;
      double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

      hYield->SetBinContent(i+1, HistCorrY);
      hYield->SetBinContent(i+1, hYield->GetBinContent(i+1)/hYield->GetBinWidth(i+1));
      hYield->SetBinError(i+1, HistCorrYerr);

      nEvents = hCent->GetBinContent(3); //binx=3: Inclusive(nEvents), binx=4: Has Jets(nTrigers)
    }//for i
  TFile* fout = TFile::Open("../kstar-in-jets/Results_draw/JetUncorrpTspectra_trg1_JetHighPt.root", "RECREATE");
  fout->cd();
  hCent->Write("nEvents_Jets");
  hYield->Write("BC_Jets_spectra");
  
  hYield->Scale(1/nEvents);
  hYield->Write("Jets_spectra");
  fout->Close();
  
  
  TCanvas* numbin = new TCanvas();
  numbin->SetRightMargin(0.05);
  numbin->SetLeftMargin(0.17);
  numbin->SetBottomMargin(0.12);

  HistoAxisTitles(hYield, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}", 0.05, 62, 1.0, 0.05, 62, 1.2);
  hYield->SetLineWidth(2);
  hYield->SetLineColor(red);
  hYield->SetMarkerStyle(43); //circle
  hYield->SetMarkerSize(3);
  hYield->SetMarkerColor(red);
  hYield->SetTitle("");
  hYield->Draw("SAME");

  TLatex* txtInfoAAA = DrawTLatex(0.58, 0.845, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.59, 0.789, "JetChHighpT, #sqrt{s} = 13.6 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.59, 0.725, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.59, 0.665, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.17, 0.78, 0.51, 0.907);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hYield, "Reconstructed yield", "pl"); //pl : point line
  legenddd->Draw("SAME");

  numbin->SetLogy();
  numbin->SaveAs("../kstar-in-jets/Plot/kstarInjets/"+localpath+"h_Yield_Jets.png", "RECREATE");
}//DrawBW
