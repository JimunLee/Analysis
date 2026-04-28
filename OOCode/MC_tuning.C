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
#include <TGraphErrors.h>
#include <cmath>

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

  // Centrality Bins setting
  const Int_t nCentBins = 1;
  Double_t CentBins[nCentBins+1] = {0.0, 100.0};

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  
  TString localpath = "MC/";

  // Data file inputs 
  TString mainfile = "../kstar-in-oo/Result_rootfile/MC_offLumi.root";
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
	  hUSSKPi[i][j] = (TH1D*)hUSS_KPi_MC->Projection(2);
	  hUSSPiK[i][j] = (TH1D*)hUSS_PiK_MC->Projection(2);

	  //Merge KPi+PiK
	  hUSS[i][j] = (TH1D*)hUSSKPi[i][j]->Clone(Form("hUSS_KPi%i_%i", i,j));
	  hUSS[i][j]->Add(hUSSPiK[i][j]);
	  hUSS[i][j]->Scale(0.5); //Counted double yield so set to control the scale
  
	  (hLSS_KPi_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  (hLSS_PiK_MC->GetAxis(1))->SetRangeUser(ptlow, pthigh);
	  hLSSKPi[i][j] = (TH1D*)hLSS_KPi_MC->Projection(2);
	  hLSSPiK[i][j] = (TH1D*)hLSS_PiK_MC->Projection(2);

	  hLSS[i][j] = (TH1D*)hLSSKPi[i][j]->Clone(Form("hLSS_KPi_%i_%i", i,j));
	  hLSS[i][j]->Add(hLSSPiK[i][j]);
	  hLSS[i][j]->Scale(0.5);

	  // Generated K*0
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
  //=================
  // M distribution
  //=================
  for (Int_t i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];
      
      for (Int_t j=0; j<nCentBins; j++)
	{
	  centlow = CentBins[j];
	  centhigh = CentBins[j+1];

	  Double_t Nlow = 1.1;
	  Double_t Nhigh = 1.199;
	  Double_t Scale = hLSS[i][j]->Integral(hLSS[i][j]->FindBin(Nlow), hLSS[i][j]->FindBin(Nhigh))
	    / hUSS[i][j]->Integral(hUSS[i][j]->FindBin(Nlow), hUSS[i][j]->FindBin(Nhigh));
      
	  double delta;
	  if(i>4){
	    delta=0.986;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==3||i==4){
	    delta = 0.992;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else if(i==1||i==2){
	    delta = 0.995;
	    hLSS[i][j]->Scale(1.0/Scale*delta);
	  }
	  else{
	    delta = 0.99; //i==0
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
	  hUSSMC_True[i][j]->Rebin(5);
      
	  hUSS[i][j]->Scale(0.2); // 1/5
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

	  TLatex* txtInfoA = DrawTLatex(0.521, 0.314, "OO yield", 0.04, 62, 1);
	  TLatex* txtInfoB = DrawTLatex(0.367, 0.357, "ALICE Work In Progress", 0.05, 62, 1);
	  TLatex* txtInfoC = DrawTLatex(0.522, 0.264, "OO #sqrt{s_{NN}} = 5.36 TeV", 0.04, 42, 1);
	  TLatex* txtInfoD = DrawTLatex(0.522, 0.214, "Min. Bias", 0.04, 42, 1);
	  TLatex* txtInfoF = DrawTLatex(0.522, 0.169, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
 
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

  TF1* fGen [nPtBins][nCentBins];
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

  double kstarMass = 0.89281;
  double kstarwidth = 0.0464;

  double BgMinRange = 0.70;
  double BgMaxRange = 1.16;

  TH1D* hYield[1];
  hYield[0] = new TH1D("binYield","binYield;pT;N", nPtBins,xBins);

  TH1D* hTrueYield[1];
  hTrueYield[0] = new TH1D("GenbinYield","GenbinYield;pT;N", nPtBins,xBins);


  int red0 = TColor::GetColor("#FF0000");     // pink and red
  int violet1 = TColor::GetColor("#BF00FF");    // hot pink
  int blue3 = TColor::GetColor("#0096FF");  // orange

  
  for (int i=0; i<nPtBins; i++)
    {
      ptlow = xBins[i];
      pthigh = xBins[i+1];

      Shreck();
      std::cout<<"*************************************************"<<std::endl;
      std::cout<<"*************************************************"<<std::endl;
      std::cout<<"THIS IS START OF FITTING FOR pTBins "<<ptlow<<" - "<<pthigh<<" with corresponding i value "<<i<<std::endl;
      std::cout<<"*************************************************"<<std::endl;
      std::cout<<"*************************************************"<<std::endl;
      
      for (int j=0; j<nCentBins; j++)
	{
	  centlow = CentBins[j];
	  centhigh = CentBins[j+1];
	  
	  hUSSFit[i][j]->Add(hUSSMC_SUB[i][j]);
	  //===================================================
	  // Gen_kstar
	  //===================================================
	  TCanvas* gen_kstar = new TCanvas(Form("gen_pT_%d_%d", i, i+1), "", 940, 800);
	  gen_kstar->SetLeftMargin(0.16);
	  gen_kstar->SetRightMargin(0.06);
	  gen_kstar->SetTopMargin(0.08);
	  gen_kstar->SetBottomMargin(0.14);
	  
	  fGen[i][j] = new TF1(Form("Genkstar_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", 0.84, 0.96);
	  fGen[i][j]->SetParameter(1, kstarMass); 
	  fGen[i][j]->SetParameter(2, kstarwidth);

	  hUSSMC_True[i][j]->Fit(fGen[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSMC_True[i][j]->SetTitle(Form(""));
	  HistoAxisTitles(hUSSMC_True[i][j], "#it{M}_{ inv}^{ K^{+}#pi^{-} }[GeV/#it{c}^{2}]",
			   "d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.7);

	  TLatex* txtInfogenA = DrawTLatex(0.183, 0.938, "ALICE Work In Progress", 0.05, 62, 1);
	  TLatex* txtInfogenB = DrawTLatex(0.191, 0.874, "OO yield", 0.04, 62, 1);
	  TLatex* txtInfogenC = DrawTLatex(0.191, 0.820, "oo #sqrt{s_{NN}} = 5.36 TeV", 0.04, 42, 1);
	  TLatex* txtInfogenD = DrawTLatex(0.191, 0.765, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  
	  gen_kstar->SaveAs(Form("../kstar-in-oo/Plot/"+localpath+"hGen_%1.1f_%1.1f.png", ptlow, pthigh), "RECREATE");
	  
	  // ==============
	  // Function
	  // ==============
	  fkstar[i][j] = new TF1(Form("kstar_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])", 0.84, 0.95);
	  fBg2[i][j] = new TF1(Form("kstar_GausBg_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6)", BgMinRange, BgMaxRange);
	  fFit[i][j] = new TF1(Form("kstar_fullFit_%i_cent_%i",i,j), "[0]*BreitWignerRelativistic(x, [1], [2])+[3]*BreitWignerRelativistic(x, [4], [5]) + pol2(6) + [9]*BreitWignerRelativistic(x, [10], [11])", BgMinRange, BgMaxRange);


	  TCanvas* lines = new TCanvas(Form("line_%d_cent_%d",i,j), "", 940, 800);
	  if(i==5){
	    for (int bin = hUSSFit[i][j]->FindBin(1.09); bin <= hUSSFit[i][j]->FindBin(1.20); bin++)
	      hUSSFit[i][j]->SetBinError(bin, 0);
	  }
	  if(i==1||i==2){
	    for (int bin = hUSSFit[i][j]->FindBin(0.60); bin <= hUSSFit[i][j]->FindBin(0.76); bin++)
	      hUSSFit[i][j]->SetBinError(bin, 0);
	  }
     
	  //=================================================================
	  if(i>6){
	    fFit[i][j]->SetParLimits(0, 0, 1000);
	    fFit[i][j]->SetParLimits(1, 0.765, 0.81);
	    fFit[i][j]->SetParLimits(2, 0.0363, 0.126);
	  }
	  else if(i==6){
	    fFit[i][j]->SetParLimits(0, 0, 600);
	    fFit[i][j]->SetParLimits(1, 0.775, 0.82);
	    fFit[i][j]->SetParLimits(2, 0.0363, 0.180);
	  }
	  else if(i==5){
	    fFit[i][j]->SetParLimits(0, 0, 1200);
	    fFit[i][j]->SetParLimits(1, 0.735, 0.78);
	    fFit[i][j]->SetParLimits(2, 0.074, 0.120);
	  }
	  else if(i==4){
	    fFit[i][j]->SetParLimits(0, 0, 1000);
	    fFit[i][j]->SetParLimits(1, 0.725, 0.78);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.120);
	  }
	  else if(i==3){
	    fFit[i][j]->SetParLimits(0, 0, 800);
	    fFit[i][j]->SetParLimits(1, 0.725, 0.78);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.120);
	  }
	  else if(i==1||i==2){
	    fFit[i][j]->SetParLimits(0, 0, 10000*hUSSFit[i][j]->GetXaxis()->GetBinWidth(5));
	    fFit[i][j]->SetParLimits(1, 0.765, 0.80);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.120);
	  }
	  else{
	    fFit[i][j]->SetParLimits(0, 0, 2000);
	    fFit[i][j]->SetParLimits(1, 0.725, 0.82);
	    fFit[i][j]->SetParLimits(2, 0.020, 0.120);
	  }

	  //=================================================================
	  if(i==9||i==10){
	    fFit[i][j]->SetParLimits(3, 0, 100);
	    fFit[i][j]->SetParLimits(4, 0.950, 0.98);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }
	  else if(i==8){
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.950, 0.98);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }
	  else if(i==7){
	    fFit[i][j]->SetParLimits(3, 0, 1000);
	    fFit[i][j]->SetParLimits(4, 0.940, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.11, 0.187);
	  }
	  else if(i==6){
	    fFit[i][j]->SetParLimits(3, 0, 1000);
	    fFit[i][j]->SetParLimits(4, 0.950, 1.04);
	    fFit[i][j]->SetParLimits(5, 0.11, 0.187);
	  }
	  else if(i==5){
	    fFit[i][j]->SetParLimits(3, 0, 3500);
	    fFit[i][j]->SetParLimits(4, 0.938, 0.98);
	    fFit[i][j]->SetParLimits(5, 0.050, 0.187);
	  }
	  else if(i==4){
	    fFit[i][j]->SetParLimits(3, 0, 1500);
	    fFit[i][j]->SetParLimits(4, 0.938, 0.948);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }
	  else if(i==3){
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.938, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }
	  else if(i==1||i==2){
	    fFit[i][j]->SetParLimits(3, 0, 800);
	    fFit[i][j]->SetParLimits(4, 0.955, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }
	  else{
	    fFit[i][j]->SetParLimits(3, 0, 2000);
	    fFit[i][j]->SetParLimits(4, 0.938, 0.99);
	    fFit[i][j]->SetParLimits(5, 0.030, 0.187);
	  }

	  fFit[i][j]->SetParameter(6, fBg2[i][j]->GetParameter(6));
	  fFit[i][j]->SetParameter(7, fBg2[i][j]->GetParameter(7));
	  fFit[i][j]->SetParameter(8, fBg2[i][j]->GetParameter(8));
	  
	  fFit[i][j]->SetParameter(9, fGen[i][j]->GetParameter(0)); //DO NOT FIX      
	  fFit[i][j]->FixParameter(10, fGen[i][j]->GetParameter(1));
	  fFit[i][j]->FixParameter(11, fGen[i][j]->GetParameter(2));

	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->Fit(fFit[i][j], "SR", "", BgMinRange, BgMaxRange);
	  hUSSFit[i][j]->SetTitle(Form(""));
	  hUSSFit[i][j]->Draw();
	  hUSSMC_True[i][j]->Draw("SAME");
	  
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
	  TLatex* txtInfoBB = DrawTLatex(0.173, 0.873, "OO yield", 0.04, 62, 1);
	  TLatex* txtInfoCC = DrawTLatex(0.173, 0.820, "oo #sqrt{s_{NN}} = 5.36 TeV", 0.04, 42, 1);
	  TLatex* txtInfoDD = DrawTLatex(0.173, 0.767, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	  TLatex* txtInfoFF = DrawTLatex(0.173, 0.713, chi2NdfText, 0.04, 42, 1);

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
	  // Recon Yield
	  double HistY = 0;
	  double HistYerr = 0;
	  double lowM = hUSSFit[i][j]->FindBin(0.80);
	  double highM = hUSSFit[i][j]->FindBin(1.10);
	  HistY=hUSSFit[i][j]->IntegralAndError(lowM, highM, HistYerr);
	  double binW = hUSSFit[i][j]->GetBinWidth(50); //Any bin number
	  double BgY = (1/binW)*totbg[i][j]->Integral(0.80,1.10);
	
	  TMatrixDSym covTot(fFit[i][j]->GetNpar());
	  TMatrixDSym covGG(bw[i][j]->GetNpar());
	  TMatrixDSym covGG_err(totbg[i][j]->GetNpar());
	  covTot = ptr->GetCovarianceMatrix();
	  covGG_err = covTot.GetSub(0, 8, 0, 8);// these have to be totbg parameters


	  cout <<"//============================================================\\" <<endl;
	  Double_t Bg_err = (1/binW)*totbg[i][j]->IntegralError(0.80, 1.10, totbg[i][j]->GetParameters(), covGG_err.GetMatrixArray());
	  //  Double_t Bg_err = (1/binW)*totbg[i][j]->IntegralError(0.80, 1.10, totbg[i][j]->GetParameters(), covGG_err.GetMatrixArray());
	 
	  double HistCorrY = HistY-BgY;
	  double HistCorrYerr = TMath::Sqrt(HistYerr*HistYerr + Bg_err*Bg_err);

	  hYield[j]->SetBinContent(i+1, HistCorrY);
	  hYield[j]->SetBinContent(i+1, hYield[j]->GetBinContent(i+1)/(hYield[j]->GetBinWidth(i+1)));
	  hYield[j]->SetBinError(i+1, HistCorrYerr);
	  /* hYield[j]->SetBinContent(i, HistCorrY); */
	  /* hYield[j]->SetBinContent(i, hYield[j]->GetBinContent(i)/(hYield[j]->GetBinWidth(i))); */
	  /* hYield[j]->SetBinError(i, HistCorrYerr); */
	   
	  hYield[j]->SetLineWidth(2);
	  hYield[j]->SetLineColor(blue3);
	  hYield[j]->SetMarkerStyle(43); //circle
	  hYield[j]->SetMarkerSize(3);
	  hYield[j]->SetMarkerColor(blue3);

	  hTrueYield[j]->SetBinContent(i+1, hUSSMC_True[i][j]->Integral());
	  
	  hTrueYield[j]->SetBinContent(i+1, hTrueYield[j]->GetBinContent(i+1)/(hTrueYield[j]->GetBinWidth(i+1)));
	  hTrueYield[j]->SetBinError(i+1, TMath::Sqrt(hTrueYield[i]->GetBinContent(i+1)));
	  
	  hTrueYield[j]->SetLineWidth(2);
	  hTrueYield[j]->SetLineColor(red0);
	  hTrueYield[j]->SetMarkerStyle(43); //circle
	  hTrueYield[j]->SetMarkerSize(3);
	   hTrueYield[j]->SetMarkerColor(red0);

	  for (int b = 1; b <= hYield[0]->GetNbinsX(); ++b) {
	    cout << "bin " << b
		 << " Yield= " << hYield[0]->GetBinContent(b)
		 << " True= " << hTrueYield[0]->GetBinContent(b)
		 << " Ratio= "<< hYield[0]->GetBinContent(b)/hTrueYield[0]->GetBinContent(b)
		 << endl;
	  }
	  cout<<"HistCorrY: "<<HistCorrY<<endl;
	  cout<<"Bg_err: "<<Bg_err<<endl;
	
	  // Gen Yield
	  double GenHistY = 0;
	  double GenHistYerr = 0;
	  double GenlowM = hUSSMC_True[i][j]->FindBin(0.70);
	  double GenhighM = hUSSMC_True[i][j]->FindBin(1.10);
	  GenHistY=hUSSMC_True[i][j]->Integral(GenlowM, GenhighM);
	  double GenbinW = hUSSMC_True[i][j]->GetBinWidth(50); //Any bin number
	  double GenBgY = (1/GenbinW)*fGen[i][j]->Integral(0.70,1.10);

	  
	  cout << "//==========================================================//" << endl;
	} //for j
    }//for i
  
  
  TCanvas* numbin = new TCanvas("","",1200,1800);
  TPad* p1 = new TPad("","",0,0.3,1,1);
  p1->SetLeftMargin(0.1);
  p1->Draw();
  p1->cd();
  p1->SetLogy();
  p1->SetLeftMargin(0.15);
  p1->SetBottomMargin(0.001);
  
  double yMin = hYield[0]->GetMinimum();
  double yMax = hYield[0]->GetMaximum();
	
  TH1D* hdummy = new TH1D("hdummy", "", 1, 0,9);
  HistoAxisTitles(hdummy, "#it{p}_{T} [GeV/#it{c}]", "#frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}",
		  0.05, 62, 1.0, 0.05, 62, 1.2);
  hdummy->SetMinimum(yMin*0.4);
  hdummy->SetMaximum(yMax*1.5);
  hdummy->SetLineWidth(2);
  hdummy->SetLineColor(1);
  hdummy->Draw();
  
  hYield[0]->Draw("SAME");


  TLatex* txtInfoAAA = DrawTLatex(0.58, 0.845, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.59, 0.789, "o-o, #sqrt{s_{NN}} = 5.36 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.59, 0.725, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.59, 0.665, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.57, 0.51, 0.904, 0.629);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hYield[0], "Reconstructed yield", "pl"); //pl : point line
  legenddd->Draw("SAME");

  //=======
  // pad2
  //=======
  numbin->cd();
  TPad* p2 = new TPad("", "", 0,0,1,0.3);
  p2->SetLeftMargin(0.15);
  p2->SetBottomMargin(0.4);
  p2->SetTopMargin(0.0001);
  p2->Draw();
  p2->cd();
  
  TH1D* hdummy2 = new TH1D("", "", 1, 0,9);
  HistoAxisTitles(hdummy2, "#it{p}_{T} [GeV/#it{c}]", "",
		  0.1, 62, 1.4);
  hdummy2->GetXaxis()->SetLabelSize(0.10);
  hdummy2->GetYaxis()->SetLabelSize(0.10);
  
  hdummy2->SetMinimum(0.5);
  hdummy2->SetMaximum(1.6);
  hdummy2->GetYaxis()->SetNdivisions(303);
  hdummy2->Draw();

  TLine* l1 = new TLine(0, 1.0, 9, 1.0); //xmin, ymin, xmax,ymax
  l1->SetLineStyle(2);
  l1->SetLineWidth(1);
  l1->SetLineColor(kBlack);
  l1->Draw("SAME");
  
  TH1D* hYieldRatio = (TH1D*)hYield[0]->Clone("hYieldRatio");
  hYieldRatio->Divide((TH1D*)hTrueYield[0]);
  hYieldRatio->Draw("SAME");

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
