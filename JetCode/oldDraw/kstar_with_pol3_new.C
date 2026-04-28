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
  TH1* hUSS_Y_Extra [nPtBins];

  for (Int_t i=0; i<nPtBins; i++)
    {
      hLSS_Y[i] = new TH1D(Form("hLSS_Y_%i", i), Form("hLSS_Y_%i", i), nBins, imlow, imtop);
      hUSS_Y[i] = new TH1D(Form("hUSS_Y_%i", i), Form("hUSS_Y_%i", i), nBins, imlow, imtop);

      hUSS_Y_Bkg[i] = new TH1D(Form("hUSS_Y_Bkg_%i", i), Form("hUSS_Y_Bkg_%i", i), nBins, imlow, imtop); 
      hUSS_Y_Peak[i] = new TH1D(Form("hUSS_Y_Peak_%i", i), Form("hUSS_Y_Peak_%i", i), nBins, imlow, imtop);
      hUSS_Y_Fit[i] = new TH1D(Form("hUSS_Y_Fit_%i", i), Form("hUSS_Y_Fit_%i", i), nBins, imlow, imtop);
      hUSS_Y_Extra[i] = new TH1D(Form("hUSS_Y_Extra_%i", i), Form("hUSS_Y_Extra_%i", i), nBins, imlow, imtop);
    }

  Double_t ptlow = 0;
  Double_t pthigh = 0;
  Double_t ReBinNum = 4;
  
  for (Int_t j=0; j<nPtBins; j++)
    {
      if(j!=1)
	continue;
 
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      (hLSS->GetYaxis())->SetRangeUser(ptlow, pthigh);
      (hUSS->GetYaxis())->SetRangeUser(ptlow, pthigh);

      hLSS_Y[j]->Add(hLSS->Project3D("z o"), 1);
      hUSS_Y[j]->Add(hUSS->Project3D("z o"), 1);
      // =======================================================================================
      // Normalization 
      Double_t Nlow = 1.20;
      Double_t Nhigh = 1.24;
      
      // if(j==4){
      // 	Nlow=0.8;
      // 	Nlow=0.83;
      // }

      hUSS_Y[j]->Scale(1.0/nEvents);
      hLSS_Y[j]->Scale(1.0/nEvents);

      Double_t Scale = hUSS_Y[j]->Integral(hUSS_Y[j]->FindBin(Nlow), hUSS_Y[j]->FindBin(Nhigh))
       	/ hLSS_Y[j]->Integral(hLSS_Y[j]->FindBin(Nlow), hLSS_Y[j]->FindBin(Nhigh));
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
  TF1* fBW_kstar [nPtBins];
  TF1* fExtra [nPtBins];
  TF1* fExtraExtra [nPtBins];

  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  TF1* fSuperFit [nPtBins];

  TF1* fBg_POST [nPtBins];
  TF1* fBW_kstar_POST [nPtBins];
  TF1* fExtra_POST [nPtBins];

  
  double Min=0.0;
  double Max=10.0;

  /* double kstarMassMin = 0.89581-0.00019; */
  /* double kstarMassMax = 0.89581+0.00019; */

  /* double omegaMassMin = 0.78265-0.00012; */
  /* double omegaMassMax = 0.78265+0.00012; */

  /* double rhoMassMin = 0.77526-0.00025; */
  /* double rhoMassMax = 0.77526+0.00025; */

  double kstarMass = 0.89581;
  double omegaMass = 0.78265;
  double rhoMass = 0.77526;
    
  double kstarMinRange = 0.75;
  double kstarMaxRange = 1.15;

  double omegaMinRange = 0.75;
  double omegaMaxRange = 1.15;

  double rhoMinRange = 0.75;
  double rhoMaxRange = 1.15;
  
  double Bglow = 0.75;
  double Bghigh = 1.25;
    

  for (int k=1; k<nPtBins+1; k++)
    {
      if(k!=1)
	continue;
      hUSS_Y_Extra[k]->Add(hUSS_Y_SUB[k]);
      hUSS_Y_Fit[k]->Add(hUSS_Y_SUB[k]);
      hUSS_Y_Bkg[k]->Add(hUSS_Y_SUB[k]);
      
      for (int bin = hUSS_Y_Bkg[k]->FindBin(0.76); bin <= hUSS_Y_Bkg[k]->FindBin(0.819); bin++){
	hUSS_Y_Bkg[k]->SetBinError(bin, 0);
     }
      for (int bin = hUSS_Y_Bkg[k]->FindBin(0.83); bin <= hUSS_Y_Bkg[k]->FindBin(0.962); bin++){
	hUSS_Y_Bkg[k]->SetBinError(bin, 0);
	hUSS_Y_Peak[k]->SetBinError(bin, 0);
	hUSS_Y_Extra[k]->SetBinError(bin, 0);
      }
	
      ptlow = xBins[k];
      pthigh = xBins[k+1];

      // ==============
      // fitting kstar
      //kstar
      fBW_kstar[k] = new TF1(Form("kstar_bw_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
      fBW_kstar[k]->SetParameter(1, kstarMass);
      fBW_kstar[k]->FixParameter(2, 0.0474);
      hUSS_Y_Fit[k]->Fit(fBW_kstar[k], "0", "", kstarMinRange, kstarMaxRange);

      
      fExtra[k] = new TF1(Form("kstar_Extra_%i",k), "pol2 + gaus(3)", 0.75, 0.85);
      hUSS_Y_Extra[k]->Fit(fExtra[k], "0", "", 0.75, 0.82);
      hUSS_Y_Extra[k]->Fit(fExtra[k], "0", "", 0.75, 0.82);
      hUSS_Y_Extra[k]->Fit(fExtra[k], "", "", 0.75, 0.82);

  

      // =====================
      // fitting background
      fBg[k] = new TF1(Form("BgFit_%i", k), "pol3(0)", Bglow, Bghigh);
      hUSS_Y_Bkg[k]->Fit(fBg[k], "0", "", Bglow, Bghigh);
      
      
      // ===============
      // fitting entire
      // if(k==1) {
      // fFit[k] = new TF1(Form("Fullfit_%i", k), "[0]*BreitWignerRelativistic(x, [1], [2]) + pol3(3)", Bglow, Bghigh);
       // fFit[k] = new TF1(Form("Fullfit_%i", k), Form("[0]*BreitWignerRelativistic(x, [1], [2]) + pol3(3)+kstar_Extra_%i(7)",k), Bglow, Bghigh);
      fFit[k] = new TF1(Form("Fullfit_%i", k), Form("kstar_bw_%i + BgFit_%i + kstar_Extra_%i",k,k,k), Bglow, Bghigh);
      
      fFit[k]->SetParameter(0, fBW_kstar[k]->GetParameter(0));
      fFit[k]->SetParameter(1, fBW_kstar[k]->GetParameter(1));
      fFit[k]->FixParameter(2, fBW_kstar[k]->GetParameter(2));

      fFit[k]->SetParameter(3, fBg[k]->GetParameter(0));
      fFit[k]->SetParameter(4, fBg[k]->GetParameter(1));
      fFit[k]->SetParameter(5, fBg[k]->GetParameter(2));
      fFit[k]->SetParameter(6, fBg[k]->GetParameter(3));

      fFit[k]->SetParameter(7, fExtra[k]->GetParameter(0));
      fFit[k]->SetParameter(8, fExtra[k]->GetParameter(1));
      fFit[k]->SetParameter(9, fExtra[k]->GetParameter(2));
      fFit[k]->SetParameter(10, fExtra[k]->GetParameter(3));
      fFit[k]->SetParameter(11, fExtra[k]->GetParameter(4));
      fFit[k]->SetParameter(12, fExtra[k]->GetParameter(5));
	  

	// if(k==4) {
	 // fFit[k] = new TF1(Form("Fullfit_%i", k), "[0]*BreitWignerRelativistic(x, [1], [2]) + pol3(3)+gaus(7)", Bglow, Bghigh);
      
	 // fFit[k]->SetParameter(0, fBW_kstar[k]->GetParameter(0));
	 // fFit[k]->FixParameter(1, fBW_kstar[k]->GetParameter(1));
	 // fFit[k]->FixParameter(2, fBW_kstar[k]->GetParameter(2));

	 // fFit[k]->FixParameter(3, fBg[k]->GetParameter(0));
	 // fFit[k]->FixParameter(4, fBg[k]->GetParameter(1));
	 // fFit[k]->FixParameter(5, fBg[k]->GetParameter(2));
	 // fFit[k]->FixParameter(6, fBg[k]->GetParameter(3));

	 // fFit[k]->SetParLimits(7, 0, 1);
	 // fFit[k]->SetParameter(8, fExtra[k]->GetParameter(1));
	 // fFit[k]->SetParameter(9, fExtra[k]->GetParameter(2));

      // }      

      // hUSS_Y_Fit[k]->GetXaxis()->SetRangeUser(0.75, 1.15);
      hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", Bglow, Bghigh);
      hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", Bglow, Bghigh);
      hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", Bglow, Bghigh);

      // =================================================================================
      // for fitting each one
      fBW_kstar_POST[k] = new TF1(Form("kstar_POST_bw_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", Bglow, Bghigh);
      fBg_POST[k] = new TF1(Form("Bg_POSTFit_%i", k), "pol3(0)", Bglow, Bghigh);
      fExtra_POST[k] = new TF1(Form("kstar_Extra_POST_%i",k), "gaus", 0.75, 0.85);

      fBW_kstar_POST[k]->SetParameter(0, fFit[k]->GetParameter(0));
      fBW_kstar_POST[k]->SetParameter(1, fFit[k]->GetParameter(1));
      fBW_kstar_POST[k]->SetParameter(2, fFit[k]->GetParameter(2));

      fBg_POST[k]->SetParameter(0, fFit[k]->GetParameter(3));
      fBg_POST[k]->SetParameter(1, fFit[k]->GetParameter(4));
      fBg_POST[k]->SetParameter(2, fFit[k]->GetParameter(5));
      fBg_POST[k]->SetParameter(3, fFit[k]->GetParameter(6));

       // fExtra_POST[k]->SetParameter(0, fFit[k]->GetParameter(7));
       // fExtra_POST[k]->SetParameter(1, fFit[k]->GetParameter(8));
       // fExtra_POST[k]->SetParameter(2, fFit[k]->GetParameter(9));
      
      TFitResultPtr ptr = hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", Bglow, Bghigh);
      
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
		      "1/#it{N}_{ ev} d#it{N}^{ K^{+}#pi^{-}}_{USS-LSS}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.04, 62, 1.0);
      hUSS_Y_Fit[k]->SetTitle(Form(""));
      
      hUSS_Y_Fit[k]->SetLineColor(1);
      hUSS_Y_Fit[k]->Draw();

      fFit[k]->SetLineColor(2);
      fFit[k]->Draw("SAME");
      
      // fBW_kstar_POST[k]->SetLineColor(6);
      // fBW_kstar_POST[k]->Draw("SAME");
      fBg_POST[k]->SetLineColor(4);
      fBg_POST[k]->SetLineStyle(2);
      fBg_POST[k]->Draw("SAME");
      // fExtra_POST[k]->Draw("SAME");

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
