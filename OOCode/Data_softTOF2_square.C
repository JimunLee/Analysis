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
  Double_t xBins[nPtBins+1] = {0.0, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 8.0};
  //Double_t xBins[nPtBins+1] = {0.0, 8.0};
  // Centrality Bins setting
  const Int_t nCentBins = 5;
  Double_t CentBins[nCentBins+1] = {0.0, 10.0, 30.0, 50.0, 70.0, 100.0};

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  
  TString localpath = "DATA/";

  // Data file inputs 
  TString mainfile = "../kstar-in-oo/Result_rootfile/softTOF2_square.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString USSKPiName("kstar-in-o-o/hUSS_KPi;1");
  TString USSPiKName("kstar-in-o-o/hUSS_PiK;1");
  TString LSSKPiName("kstar-in-o-o/hLSS_KPi;1");
  TString LSSPiKName("kstar-in-o-o/hLSS_PiK;1");
 
  TString CentName("kstar-in-o-o/hcentFT0C_AC;1");
  
  // eventdata counts
  //TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  // Double_t nEvents = hVtx->GetBinContent(1); //number 1 is clear number

  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  THnSparseT<TArrayF> *hUSS_KPi_MC = (THnSparseT<TArrayF>*)data->Get(USSKPiName);
  THnSparseT<TArrayF> *hUSS_PiK_MC = (THnSparseT<TArrayF>*)data->Get(USSPiKName);
  THnSparseT<TArrayF> *hLSS_KPi_MC = (THnSparseT<TArrayF>*)data->Get(LSSKPiName);
  THnSparseT<TArrayF> *hLSS_PiK_MC = (THnSparseT<TArrayF>*)data->Get(LSSPiKName);
  //THnSparseT<TArrayF> *hCentFT0C = (THnSparseT<TArrayF>*)data->Get(CentName);

  TH1D* hUSSKPi [nPtBins][nCentBins];
  TH1D* hUSSPiK [nPtBins][nCentBins];
  TH1D* hUSS [nPtBins][nCentBins];

  TH1D* hLSSKPi [nPtBins][nCentBins];
  TH1D* hLSSPiK [nPtBins][nCentBins];
  TH1D* hLSS [nPtBins][nCentBins];
  
  TH1D* hLSSMC_SUB [nPtBins][nCentBins];
  TH1D* hUSSMC_SUB [nPtBins][nCentBins];

  TH1D* hCent [nCentBins];

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
	  
	  hUSSMC_SUB[i][j] = new TH1D(Form("hUSS_SUB_%i_%i", i,j), Form("hUSS_SUB_%i_%i", i,j), nBins, imlow, imtop);
	  hLSSMC_SUB[i][j] = new TH1D(Form("hLSS_SUB_%i_%i", i,j), Form("hLSS_SUB_%i_%i", i,j), nBins, imlow, imtop);
	  hUSSFit[i][j] = new TH1D(Form("hUSSFit_%i_%i", i,j), Form("hUSSFit_%i_%i", i,j), nBins/5, imlow, imtop); 

	  
	  (hUSS_KPi_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hUSS_PiK_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hUSS_KPi_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);
	  (hUSS_PiK_MC->GetAxis(0))->SetRangeUser(centlow, centhigh);

	  hUSSKPi[i][j] = (TH1D*)hUSS_KPi_MC->Projection(2); // Mass
	  hUSSPiK[i][j] = (TH1D*)hUSS_PiK_MC->Projection(2);
      
	  hUSS[i][j] = (TH1D*)hUSSKPi[i][j]->Clone(Form("hUSS_KPi%i_%i", i,j));
	  hUSS[i][j]->Add(hUSSPiK[i][j]); //KPi+PiK
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

	  hUSS[i][j]->SetName(Form("hUSS_%i_%i", i,j));
	  hLSS[i][j]->SetName(Form("hLSS_%i_%i", i,j));
	  
	  hCent[j] = (TH1D*)data->Get(CentName);
	  hCent[j]->SetName(Form("hUSS_Cent_%i_%i", i,j));
	  hCent[j]->GetXaxis()->SetRangeUser(centlow, centhigh);

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
	  
	  Double_t Nlow = 1.10;
	  Double_t Nhigh = 1.1999;
	  Double_t Scale = hLSS[i][j]->Integral(hLSS[i][j]->FindBin(Nlow), hLSS[i][j]->FindBin(Nhigh))
	    / hUSS[i][j]->Integral(hUSS[i][j]->FindBin(Nlow), hUSS[i][j]->FindBin(Nhigh));
      
	  double delta;
	  if (i<4){
	    if(j==0)
	      delta = 0.999;
	    else if(j==3){
	      if(i==0)
		delta= 0.999;
	      else
		delta = 0.996;
	    }
	    else if(i==3&&j==4){
	      delta = 0.99;
	    }
	    else if((i<4)&&j==1)
	      delta = 0.999;
	    else if(j!=4){
	      delta = 1.0;
	    }
	    else{
	      delta = 0.992;
	    }
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if (i==4){
	    if (j==3||j==0)
	      delta = 0.994;
	    else if (j==2)
	      delta = 0.996;
	    else
	      delta = 0.992;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if ((i==5||i==6) && j==2){
	    delta = 0.996;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==6&&j==1){
	    delta = 0.998;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==5&&j==1){
	    delta = 0.994;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==5&&j==0){
	    delta = 0.992;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==6&&j==0){
	    delta = 0.997;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else{
	    delta = 0.986;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
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
      
	  hUSS[i][j]->Scale(0.2);
	  hLSS[i][j]->Scale(0.2);
	  
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

	  TLatex* txtInfoA = DrawTLatex(0.521, 0.358, "O-O yield", 0.04, 62, 1);
	  TLatex* txtInfoB = DrawTLatex(0.339, 0.409, "ALICE Work In Progress", 0.05, 62, 1);
	  TLatex* txtInfoC = DrawTLatex(0.522, 0.313, "o-o #sqrt{s_{NN}} = 5.36 TeV", 0.04, 42, 1);
	  TLatex* txtInfoD = DrawTLatex(0.522, 0.264, "Min. Bias", 0.04, 42, 1);
	  TLatex* txtInfoF = DrawTLatex(0.522, 0.169, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  TLatex* txtInfoG = DrawTLatex(0.526, 0.216, Form("Centrality %d-%d %%", (int)centlow, (int)(centhigh)), 0.04, 42, 1);
 
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

  double kstarMass;
  double kstarwidth;

  double BgMinRange = 0.70;
  double BgMaxRange = 1.16;

  TH1D* hYield[5];
  hYield[0] = new TH1D("binYield_0_10","binYield;pT;N", nPtBins,xBins);
  hYield[1] = new TH1D("binYield_10_30","binYield;pT;N", nPtBins,xBins);
  hYield[2] = new TH1D("binYield_30_50","binYield;pT;N", nPtBins,xBins);
  hYield[3] = new TH1D("binYield_50_70","binYield;pT;N", nPtBins,xBins);
  hYield[4] = new TH1D("binYield_70_100","binYield;pT;N", nPtBins,xBins);

  int red = TColor::GetColor("#FF0000");
  int orange = TColor::GetColor("#ffa200");
  int green = TColor::GetColor("#09b347");
  int violet = TColor::GetColor("#BF00FF");
  int blue = TColor::GetColor("#176ceb");

  double nEvents = 0.0;
  double nEventsArr [nCentBins];
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
	  fBg2[i][j] = new TF1(Form("kstar_GausBg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6)", BgMinRange, BgMaxRange);
	  fFit[i][j] = new TF1(Form("kstar_fullFit_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6) + [9]*BreitWignerRelativistic(x, [10], [11])", BgMinRange, BgMaxRange);


	  TCanvas* lines = new TCanvas(Form("line_%d_cent_%d",i,j), "", 940, 800);

	  fkstar[i][j]->SetParameter(1, 0.890);
	  fkstar[i][j]->SetParameter(2, 0.0524);
	  //=================================================================
	  if((i==4||i==5)&&j<4){
	    for (int bin = hUSSFit[i][j]->FindBin(1.05); bin <= hUSSFit[i][j]->FindBin(1.20); bin++)
	      hUSSFit[i][j]->SetBinError(bin, 0);
	  }
	  if(i==3&&j==0){
	    for (int bin = hUSSFit[i][j]->FindBin(1.02); bin <= hUSSFit[i][j]->FindBin(1.1); bin++)
	      hUSSFit[i][j]->SetBinError(bin, 0);
	  }
	  if(i==4&&j==0){
	     for (int bin = hUSSFit[i][j]->FindBin(1.03); bin <= hUSSFit[i][j]->FindBin(1.20); bin++)
	    hUSSFit[i][j]->SetBinError(bin, 0);
	  }
	  //=================================================================
	  if((i==8&&j>2)||i==9||i==10){
	    fFit[i][j]->SetParLimits(0, 0, 20);
	    fFit[i][j]->SetParLimits(1, 0.785, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if(j==4&&(i<6)){
	    fFit[i][j]->SetParLimits(0, 0, 300);
	    fFit[i][j]->SetParLimits(1, 0.75, 0.81);
	    fFit[i][j]->SetParLimits(2, 0.070, 0.120);
	  }
	  else if(i==6&&j==2){
	    fFit[i][j]->SetParLimits(0, 0, 800);
	    fFit[i][j]->SetParLimits(1, 0.725, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if (i==6&&j==1){
	    fFit[i][j]->SetParLimits(0, 0, 3500);
	    fFit[i][j]->SetParLimits(1, 0.775, 0.81);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if (i==6&&j==0){
	    fFit[i][j]->SetParLimits(0, 0, 2200);
	    fFit[i][j]->SetParLimits(1, 0.74, 0.82);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.150);
	  }
	  else if (i==5&&j==1){
	    fFit[i][j]->SetParLimits(0, 0, 3600);
	    fFit[i][j]->SetParLimits(1, 0.746, 0.802);
	    fFit[i][j]->SetParLimits(2, 0.040, 0.090);
	  }
	  else if (i==5&&j==0){
	    fFit[i][j]->SetParLimits(0, 0, 2800);
	    fFit[i][j]->SetParLimits(1, 0.74, 0.802);
	    fFit[i][j]->SetParLimits(2, 0.040, 0.090);
	  }
	  else if((i==6&&j!=2)||i==7){
	    fFit[i][j]->SetParLimits(0, 0, 700);
	    fFit[i][j]->SetParLimits(1, 0.796, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.120);
	  }
	  else if((i==3||i==4||i==5)&&j==3){
	    fFit[i][j]->SetParLimits(0, 0, 500);
	    fFit[i][j]->SetParLimits(1, 0.75, 0.82);
	    fFit[i][j]->SetParLimits(2, 0.070, 0.120);
	  }
	  else if((i==1||i==2)&&j==3){
	    fFit[i][j]->SetParLimits(0, 0, 300);
	    fFit[i][j]->SetParLimits(1, 0.75, 0.81);
	    fFit[i][j]->SetParLimits(2, 0.070, 0.120);
	  }
	  else if (i==5&&j==2){
	    fFit[i][j]->SetParLimits(0, 0, 2000);
	    fFit[i][j]->SetParLimits(1, 0.765, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.060, 0.150);
	  }
	  else if (i==4&&j==2){
	    fFit[i][j]->SetParLimits(0, 0, 1300);
	    fFit[i][j]->SetParLimits(1, 0.746, 0.804);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.150);
	  }
	  else if(i==4&&j==0){
	    fFit[i][j]->SetParLimits(0, 0, 3000);
	    fFit[i][j]->SetParLimits(1, 0.735, 0.802);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.090);
	  }
	  else if(i==4&&j==1){
	    fFit[i][j]->SetParLimits(0, 0, 5000);
	    fFit[i][j]->SetParLimits(1, 0.735, 0.802);
	    fFit[i][j]->SetParLimits(2, 0.040, 0.090);
	  }
	  else if (i<4&&j==2){
	    fFit[i][j]->SetParLimits(0, 0, 1000);
	    fFit[i][j]->SetParLimits(1, 0.745, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if ((i==2||i==3)&&j==1){
	    fFit[i][j]->SetParLimits(0, 0, 3800);
	    fFit[i][j]->SetParLimits(1, 0.736, 0.810);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if ((i==0||i==1)&&(j==1||j==0)){
	    fFit[i][j]->SetParLimits(0, 0, 2800);
	    fFit[i][j]->SetParLimits(1, 0.736, 0.810);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }
	  else if(i==2&&j==0){
	    fFit[i][j]->SetParLimits(0, 0, 5000);
	    fFit[i][j]->SetParLimits(1, 0.71, 0.80);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.150);
	  }
	  else if (i==3&&j==0){
	    fFit[i][j]->SetParLimits(0, 0, 5000);
	    fFit[i][j]->SetParLimits(1, 0.735, 0.802);
	    fFit[i][j]->SetParLimits(2, 0.040, 0.090);
	  }
	  else{
	    fFit[i][j]->SetParLimits(0, 0, 20000);
	    fFit[i][j]->SetParLimits(1, 0.725, 0.814);
	    fFit[i][j]->SetParLimits(2, 0.030, 0.150);
	  }

	  //===============================================
	  if((i==8&&j>2)||i==9||i==10){
	    fFit[i][j]->SetParLimits(3, 0, 10);
	    fFit[i][j]->SetParLimits(4, 0.95, 0.98);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if(j==4&&(i<6)){
	    fFit[i][j]->SetParLimits(3, 0, 500);
	    fFit[i][j]->SetParLimits(4, 0.965, 0.980);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if(i==6&&j==2){
	    fFit[i][j]->SetParLimits(3, 0, 1000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.180);
	  }
	  else if (i==6&&j==1){
	    fFit[i][j]->SetParLimits(3, 0, 2800);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.155);
	  }
	  else if (i==6&&j==0){
	    fFit[i][j]->SetParLimits(3, 0, 5000);
	    fFit[i][j]->SetParLimits(4, 0.96, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.180);
	  }
	  else if(i==6||i==7){
	    fFit[i][j]->SetParLimits(3, 0, 300);
	    fFit[i][j]->SetParLimits(4, 0.96, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if((i==3||i==5)&&j==3){
	     fFit[i][j]->SetParLimits(3, 0, 600);
	     fFit[i][j]->SetParLimits(4, 0.95, 0.980);
	     fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	 
	  else if((i<3)&&j==3){
	    fFit[i][j]->SetParLimits(3, 0, 600);
	    fFit[i][j]->SetParLimits(4, 0.95, 0.980);
	    fFit[i][j]->SetParLimits(5, 0.070, 0.120);
	  }
	  else if(i==5&&j==2){
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.180);
	  }
	  else if(i==5&&j==1){
	    fFit[i][j]->SetParLimits(3, 0, 4500);
	    fFit[i][j]->SetParLimits(4, 0.93, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.080, 0.180);
	  }
	  else if(i==5&&j==0){
	    fFit[i][j]->SetParLimits(3, 0, 4500);
	    fFit[i][j]->SetParLimits(4, 0.93, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.080, 0.180);
	  }
	  else if(i==4&&j==2){
	    fFit[i][j]->SetParLimits(3, 0, 2400);
	    fFit[i][j]->SetParLimits(4, 0.95, 0.995);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if(i<4&&j==2){
	    fFit[i][j]->SetParLimits(3, 0, 1000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.0);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if(i==4&&(j==0||j==1)){
	    fFit[i][j]->SetParLimits(3, 0, 4500);
	    fFit[i][j]->SetParLimits(4, 0.92, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.180);
	  }
	  else if(i==4&&j==3){
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.95, 0.980);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  else if ((i==2||i==3)&&j==1){
	    fFit[i][j]->SetParLimits(3, 0, 3000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.04);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.120);
	  }
	
	  else if ((i==0||i==1)&&(j==1||j==0)){
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.04);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.120);
	  }
	  else if(i==2&&j==0){
	    fFit[i][j]->SetParLimits(3, 0, 5000);
	    fFit[i][j]->SetParLimits(4, 0.94, 1.03);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.100);
	  }
	  else if(i==3&&j==0){
	    fFit[i][j]->SetParLimits(3, 0, 5000);
	    fFit[i][j]->SetParLimits(4, 0.93, 1.00);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.050);
	  }
	  else{
	    fFit[i][j]->SetParLimits(3, 0, 30000);
	    fFit[i][j]->SetParLimits(4, 0.95, 1.02);
	    fFit[i][j]->SetParLimits(5, 0.020, 0.120);
	  }
	  
	  fFit[i][j]->SetParameter(6, fBg2[i][j]->GetParameter(6));
	  fFit[i][j]->SetParameter(7, fBg2[i][j]->GetParameter(7));
	  fFit[i][j]->SetParameter(8, fBg2[i][j]->GetParameter(8));

	  fFit[i][j]->SetParameter(9, fkstar[i][j]->GetParameter(0)); //DO NOT FIX
	  fFit[i][j]->SetParameter(10, fkstar[i][j]->GetParameter(1));
	  fFit[i][j]->SetParameter(11, fkstar[i][j]->GetParameter(2));
       
	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->SetTitle(Form(""));
	  hUSSFit[i][j]->Draw();
      
	  gausbg[i][j] = new TF1(Form("gausbg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
	  f0bg[i][j] = new TF1(Form("f0bg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
	  pol2bg[i][j] = new TF1(Form("pol2bg_%i_cent_%i",i,j), "pol2(0)", BgMinRange, BgMaxRange);
	  bw[i][j] = new TF1(Form("bw_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", BgMinRange, BgMaxRange);
	  totbg[i][j] = new TF1(Form("totbg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2]) + [3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6)", BgMinRange, BgMaxRange);

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

	  TLatex* txtInfoAA = DrawTLatex(0.167, 0.93, "ALICE Work In Progress", 0.05, 62, 1);
	  TLatex* txtInfoBB = DrawTLatex(0.173, 0.873, "O-O yield", 0.04, 62, 1);
	  TLatex* txtInfoCC = DrawTLatex(0.173, 0.820, "o-o #sqrt{s_{NN}} = 5.36 TeV", 0.04, 42, 1);
	  TLatex* txtInfoDD = DrawTLatex(0.173, 0.767, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  TLatex* txtInfoEE = DrawTLatex(0.173, 0.665, Form("Centrality %d-%d %%", (int)(centlow), (int)(centhigh)), 0.04, 42, 1);
	  TLatex* txtInfoFF = DrawTLatex(0.173, 0.714, chi2NdfText, 0.04, 42, 1);

	  TLegend *legendd = new TLegend(0.72,0.74,0.84,0.90);
	  legendd->SetFillColor(4000);
	  legendd->SetBorderSize(0);
	  legendd->SetTextSize(0.03);

	  legendd->AddEntry(totbg[i][j], "Tot Bkg", "l");
	  legendd->AddEntry(bw[i][j], "K^{*0}", "l");
	  legendd->AddEntry(gausbg[i][j], "#rho", "l");
	  legendd->AddEntry(f0bg[i][j], "f0", "l");
	  legendd->AddEntry(pol2bg[i][j], "Residual bkg.", "l");
	  legendd->Draw("SAME");
      
	  lines->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"two_fit_%1.1f_%1.1f_cent_%1.1f_%1.1f.png", ptlow, pthigh, centlow, centhigh), "RECREATE");
	  
	  //=============================================
	  // Yield Main
	  //=============================================
	  double HistY = 0;
	  double HistYerr = 0;
	  double lowM = hUSSFit[i][j]->FindBin(BgMinRange);
	  double highM = hUSSFit[i][j]->FindBin(BgMaxRange);
	  HistY=hUSSFit[i][j]->IntegralAndError(lowM, highM, HistYerr);
	  double binW = hUSSFit[i][j]->GetBinWidth(50); //Any number
	  double BgY = totbg[i][j]->Integral(BgMinRange, BgMaxRange)/binW;
	
	  TMatrixDSym covTot(fFit[i][j]->GetNpar());
	  TMatrixDSym covGG(bw[i][j]->GetNpar());
	  TMatrixDSym covGG_err(totbg[i][j]->GetNpar());
	  covTot = ptr->GetCovarianceMatrix();
	  covGG_err = covTot.GetSub(0,8,0,8);// these have to be totbg parameters
	  
	  Double_t Bg_err = (1/binW)*(totbg[i][j]->IntegralError(BgMinRange, BgMaxRange, totbg[i][j]->GetParameters(),covGG_err.GetMatrixArray()));

	  double HistCorrY = HistY-BgY;
	  double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	  hYield[j]->SetBinContent(i+1, HistCorrY);
	  hYield[j]->SetBinContent(i+1, hYield[j]->GetBinContent(i+1)/(hYield[j]->GetBinWidth(i+1)));
	  hYield[j]->SetBinError(i+1, HistCorrYerr);

	  double CLow = hCent[j]->FindBin(centlow);
	  double CHigh = hCent[j]->FindBin(centhigh);

	  nEvents = hCent[j]->Integral(CLow, CHigh);
	  
	  cout<<"HistCorrY: "<<HistCorrY<<endl;
	  cout<<"Bg_err: "<<Bg_err<<endl;
	  cout << "nEvents1: " << nEvents << endl;
	  Shreck();
	  nEventsArr[j]=nEvents;
	
	  //Extra work
	  //0-10%
	  if(j==0){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(red);
	    hYield[j]->SetMarkerStyle(43);
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(red);
	  }
	  //10-30%
	  if(j==1){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(orange);
	    hYield[j]->SetMarkerStyle(93); //hospital
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(orange);
	  }
	  //30-50%
	  if(j==2){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(green);
	    hYield[j]->SetMarkerStyle(94); //star
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(green);
	  }
	  //50-70%
	  if(j==3){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(blue);
	    hYield[j]->SetMarkerStyle(101); //dia
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(blue);
	  }
	  //70-100%
	  if(j==4){
	    hYield[j]->SetLineWidth(2);
	    hYield[j]->SetLineColor(violet);
	    hYield[j]->SetMarkerStyle(92); //realdia
	    hYield[j]->SetMarkerSize(2);
	    hYield[j]->SetMarkerColor(violet);
	  }
	  
	} //for j  
    }//for i

  TCanvas* numbin = new TCanvas("","",1200,800);
  numbin->SetRightMargin(0.05);
  numbin->SetLeftMargin(0.17);
  numbin->SetBottomMargin(0.12);

  double yMin = 1e-6;
  double yMax = 1e-1;
	
  TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
  hdummy->SetMinimum(yMin*0.5);
  hdummy->SetMaximum(yMax*1.5);

  HistoAxisTitles(hdummy, "#it{p}_{T} [GeV/#it{c}]", "1/N_{evt}*d^{2}N_{K^{*0}}/d#it{p}_{T}d#eta",
		  0.05, 62, 1.0, 0.05, 62, 1.4);
  hdummy->SetLineWidth(2);
  hdummy->SetLineColor(1);
  hdummy->Draw();
  

  for (int j=0; j<nCentBins; j++){
    hYield[j]->Scale(1/nEventsArr[j]);
    cout<<"nEvents: "<<nEventsArr[j]<<endl;
    hYield[j]->Draw("SAME");
  }
  
  
  TLatex* txtInfoAAA = DrawTLatex(0.58, 0.845, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.59, 0.789, "o-o, #sqrt{s_{NN}} = 5.36 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.59, 0.725, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.59, 0.665, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.64, 0.54, 0.97, 0.74);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->Draw("SAME");

  TLegend* cLegend = new TLegend(0.189,0.147,0.423,0.341); //LEFT x, y, Right x, y
  cLegend->SetFillColor(4000);
  cLegend->SetBorderSize(0);
  cLegend->SetTextSize(0.04);
  cLegend->AddEntry(hYield[0], "0-10% Cent", "lpf");
  cLegend->AddEntry(hYield[1], "10-30% Cent", "lpf");
  cLegend->AddEntry(hYield[2], "30-50% Cent", "lpf");
  cLegend->AddEntry(hYield[3], "50-70% Cent", "lpf");
  cLegend->AddEntry(hYield[4], "70-100% Cent", "lpf");
  cLegend->Draw("SAME");

  numbin->SetLogy();
  numbin->SaveAs("../kstar-in-oo/Plot/"+localpath+"h_Yield.png", "RECREATE");

  /* TFile *fout = TFile::Open("../kstar-in-oo/Result_rootfile/uncorrected_spectra_kstar.root", "RECREATE"); */
  /* numbin->Write(); */
  /* fout->Close(); */

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
