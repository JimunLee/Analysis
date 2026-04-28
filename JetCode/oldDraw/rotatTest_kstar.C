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

void DrawBW(TString Kind, TString USS, TString LSS)
{
  cout.precision(6);
  gStyle->SetOptStat(0);
  
  // pTBins setting 
  const Int_t nPtBins = 8;
  Double_t xBins[nPtBins+1];

  // TLatex setting
  TLatex* SpecialLatex_plot = nullptr;
  TLatex* SpecialLatex_fit = nullptr;
  TLatex* SpecialLatex_plot2 = nullptr;
  TLatex* SpecialLatex_fit2 = nullptr;

  if (Kind == "Inclusive"){
    Double_t inclusive_xBins[nPtBins+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.0, 8.0};
    copy(begin(inclusive_xBins), end(inclusive_xBins), xBins);

    SpecialLatex_fit = DrawTLatex(0.557, 0.79, "Min. Bias", 0.04, 42, 1);
  }
  else if(Kind == "Inside"){
    Double_t inside_xBins[nPtBins+1] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.0, 8.0};
    copy(begin(inside_xBins), end(inside_xBins), xBins);

    SpecialLatex_fit = DrawTLatex(0.558, 0.785, "#Delta #it{R}_{ jet}^{ K*} < 0.4", 0.04, 42, 1);
    SpecialLatex_fit2 = DrawTLatex(0.558, 0.630, "#it{p}_{T}^{ Jet ch uncorr} > 8 GeV/c", 0.04, 42, 1);
  }
  else{
    cout << "Please specify either 'Inside', 'Outside' or 'Inclusive'." << endl;
    return;
  }

  // invariant mass range
  float imlow = 0.60;
  float imtop = 1.20;
  /* Double_t imlow = 0.55; */
  /* Double_t imtop = 1.25; */
  
  TString widelocal = "kstarInjets/MC/";
  TString localpath = widelocal + Kind;

  // Data file inputs
  TString mainfile = "../kstar-in-jets/Result_rootfile/MB_fullrange.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString TrueName("phi-in-jets/hMCRecTrue_"+USS+";1");
  TString USSName("phi-in-jets/hMCRec_"+USS+";1");
  TString LSSName("phi-in-jets/hMCRec_"+LSS+";1");
  TString R_LSSName("phi-in-jets/hMCRec_"+USS+"_Rotational;1");
  TString Eventdata("phi-in-jets/nEvents_MCRec;1");
  
  // eventdata counts
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(2); //number 2 is clear number

  //============================================================================
  // 1D projection plots : raw yields plot of invariant mass following each pT
  
  THnSparseT<TArrayF> *hLSS_MCRecRot = (THnSparseT<TArrayF>*)data->Get(R_LSSName);
  THnSparseT<TArrayF> *hUSS_MCRecTrue = (THnSparseT<TArrayF>*)data->Get(TrueName);
  THnSparseT<TArrayF> *hUSS_MCRec = (THnSparseT<TArrayF>*)data->Get(USSName);
  THnSparseT<TArrayF> *hLSS_MCRec = (THnSparseT<TArrayF>*)data->Get(LSSName);


  TH1D* hLSS_Y_MCRot [nPtBins];
  TH1D* hUSS_Y_MCTrue [nPtBins];
  TH1D* hUSS_Y_MC [nPtBins];
  TH1D* hLSS_Y_MC [nPtBins];
  TH1D* hUSS_Y_MC_SUB [nPtBins];
  TH1D* hLSS_Y_MC_SUB [nPtBins];
  TH1D* hSub_USSLSS [nPtBins];

  TH1* hUSS_Y_Fit [nPtBins];
  TH1* hUSS_Y_Bg [nPtBins];
  TH1* hUSS_Y_kstarSig [nPtBins];
  
  Double_t ptlow = 0;
  Double_t pthigh = 0;
  Double_t ReBinNum = 4;

  TH1D *hProj_Rot;
  TH1D *hProj_True;
  TH1D *hProj_USS;
  TH1D *hProj_LSS;
  
    for (Int_t i=0; i<nPtBins; i++)
    {
      Double_t nBins = 700.0;
      hLSS_Y_MCRot[i] = new TH1D(Form("hLSS_Y_MCRot%i", i), Form("hLSS_Y_MCRot%i", i), nBins, imlow, imtop);
      hUSS_Y_MCTrue[i] = new TH1D(Form("hUSS_Y_MCTrue%i", i), Form("hUSS_Y_MCTrue%i", i), nBins, imlow, imtop);
      hUSS_Y_MC[i] = new TH1D(Form("hUSS_Y_MC%i", i), Form("hUSS_Y_MC%i", i), nBins, imlow, imtop);
      hLSS_Y_MC[i] = new TH1D(Form("hLSS_Y_MC%i", i), Form("hLSS_Y_MC%i", i), nBins, imlow, imtop);
      hSub_USSLSS[i] = new TH1D(Form("hSub_USSLSS%i", i), Form("hSub_USSLSS%i", i), nBins, imlow, imtop);

      hUSS_Y_Bg[i] = new TH1D(Form("hUSS_Y_Bg_%i", i), Form("hUSS_Y_Bg_%i", i), nBins, imlow, imtop);
      hUSS_Y_kstarSig[i] = new TH1D(Form("hUSS_Y_kstarSig_%i", i), Form("hUSS_Y_kstarSig_%i", i), nBins, imlow, imtop);
      hUSS_Y_Fit[i] = new TH1D(Form("hUSS_Y_Fit_%i", i), Form("hUSS_Y_Fit_%i", i), nBins, imlow, imtop);
      
      ptlow = xBins[i];
      pthigh = xBins[i+1];
      
      (hLSS_MCRecRot->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hUSS_MCRecTrue->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hUSS_MCRec->GetAxis(1))->SetRangeUser(ptlow, pthigh);
      (hLSS_MCRec->GetAxis(1))->SetRangeUser(ptlow, pthigh);

      hProj_Rot = hLSS_MCRecRot->Projection(2);
      hProj_True = hUSS_MCRecTrue->Projection(2);
      hProj_USS = hUSS_MCRec->Projection(2);
      hProj_LSS = hLSS_MCRec->Projection(2);

      hLSS_Y_MCRot[i]->Add(hProj_Rot);
      hUSS_Y_MCTrue[i]->Add(hProj_True);
      hUSS_Y_MC[i]->Add(hProj_USS);
      hLSS_Y_MC[i]->Add(hProj_LSS);

      delete hProj_Rot;
      delete hProj_True;
      delete hProj_USS;
      delete hProj_LSS;
    }//for1


 
  for (Int_t j=0; j<nPtBins; j++)
    {
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      Double_t Nlow = 1.19;
      Double_t Nhigh = 1.24;

      Double_t Scale = hLSS_Y_MC[j]->Integral(hLSS_Y_MC[j]->FindBin(Nlow), hLSS_Y_MC[j]->FindBin(Nhigh))
	/ hUSS_Y_MC[j]->Integral(hUSS_Y_MC[j]->FindBin(Nlow), hUSS_Y_MC[j]->FindBin(Nhigh));
      
      double delta =1.0;
      hLSS_Y_MC[j]->Scale(1.0);
      hUSS_Y_MC[j]->Scale(Scale*delta);
    
      // 1D projection plots
      TCanvas* plots = new TCanvas(Form("plots_%d", j), "", 940, 800);
      plots->SetLeftMargin(0.16);
      plots->SetRightMargin(0.06);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);

      HistoAxisTitles(hUSS_Y_MC[j], "#it{M}_{ inv}^{ K^{+}#pi^{-}} [GeV/#it{c}^{2}]", "d#it{N}^{ K#pi}/d#it{M}_{ inv}", 0.06, 62, 1.0, 0.06, 62, 1.2);
      hUSS_Y_MC[j]->SetTitle(Form(""));

      hUSS_Y_MC[j]->SetMarkerStyle(4); //circle
      hUSS_Y_MC[j]->SetMarkerSize(1.2);
      hUSS_Y_MC[j]->SetMarkerColor(4);
      hUSS_Y_MC[j]->SetLineColor(4); 

      hLSS_Y_MC[j]->SetMarkerStyle(26); //triangle
      hLSS_Y_MC[j]->SetMarkerSize(1.2);
      hLSS_Y_MC[j]->SetMarkerColor(2);
      hLSS_Y_MC[j]->SetLineColor(2);

      /* double maxY = max(hUSS_Y_MC[j]->GetMaximum(), hLSS_Y_MCRot[j]->GetMaximum()); */
      /* double minY = min(hUSS_Y_MC[j]->GetMinimum(), hLSS_Y_MCRot[j]->GetMinimum()); */
      
      /* hUSS_Y_MC[j]->SetMaximum(maxY * 1.1);  */
      /* hUSS_Y_MC[j]->SetMinimum(minY * 0.9);   */

      hUSS_Y_MC[j]->Draw();
      hLSS_Y_MC[j]->Draw("SAME");

      if (Kind == "Inclusive"){
	TLatex* txtInfoA = DrawTLatex(0.641, 0.876, Kind+" yield", 0.04, 62, 1);
	TLatex* txtInfoB = DrawTLatex(0.641, 0.929, "ALICE #bf{Performance}", 0.05, 62, 1);
	TLatex* txtInfoC = DrawTLatex(0.641, 0.825, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
	TLatex* txtInfoD = DrawTLatex(0.641, 0.775, "Min. Bias", 0.04, 42, 1);
	TLatex* txtInfoF = DrawTLatex(0.641, 0.723, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
	if (j==0 || j==1){
	  TLegend* legend = new TLegend(0.18,0.17,0.48,0.32); //LEFT x, y, Right x, y
	  legend->SetFillColor(4000);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.05);
	  legend->AddEntry(hLSS_Y_MCRot[j], "K^{+}#pi^{-}", "lpf");
	  legend->AddEntry(hLSS_Y_MC[j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
	  legend->Draw("SAME");
	}
	else {
	  TLegend* legend = new TLegend(0.53,0.16,0.83,0.31); //LEFT x, y, Right x, y
	  legend->SetFillColor(4000);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.05);
	  legend->AddEntry(hLSS_Y_MCRot[j], "K^{+}#pi^{-}", "lpf");
	  legend->AddEntry(hLSS_Y_MC[j], "2#sqrt{#it{M}_{inv}^{ K^{+}#pi^{+}}*#it{M}_{inv}^{ K^{-}#pi^{-}}}", "lpf");
	  legend->Draw("SAME");
	}
      }
      plots->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/MCRec_projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");

      // ==========================
      // Drawing Substraction
      hUSS_Y_MC_SUB[j] = (TH1D*)hUSS_Y_MC[j]->Clone(Form("hUSS_Y_MC_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_MC_SUB[j] = (TH1D*)hLSS_Y_MC[j]->Clone(Form("hLSS_Y_MC_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      TCanvas* substraction = new TCanvas(Form("MC_substraction_%d", j), "", 550, 440);
      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.08);
      substraction->SetBottomMargin(0.14);

      hUSS_Y_MC_SUB[j]->Add(hLSS_Y_MC_SUB[j], -1);
      hUSS_Y_MC_SUB[j]->Rebin(5);
      hUSS_Y_MC_SUB[j]->Scale(1.0/5.0);
      
      hUSS_Y_MC_SUB[j]->Draw();
      hUSS_Y_MCTrue[j]->Rebin(5);
      hUSS_Y_MCTrue[j]->Scale(1.0/5.0);
      hUSS_Y_MCTrue[j]->Draw("SAME");      
      HistoAxisTitles(hUSS_Y_MC_SUB[j], "#it{M}_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]", "Counts");
      hUSS_Y_MC_SUB[j]->SetTitle(Form(""));

      TLatex* txtInfoa = DrawTLatex(0.55, 0.88, Kind+" yield", 0.04, 62, 1);
      TLatex* txtInfob = DrawTLatex(0.55, 0.93, "ALICE #bf{Performance}", 0.05, 62, 1);
      TLatex* txtInfoc = DrawTLatex(0.55, 0.83, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1); 
      if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
      TLatex* txtInfoe = DrawTLatex(0.55, 0.741, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);

      substraction->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/MCRec_reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh));
      
    }//for2
  TF1* fkstar [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  TF1* fBg_POST [nPtBins];
  TF1* fkstar_POST [nPtBins];
  TF1* fMCkstar [nPtBins];

  double kstarMass = 0.89581;
  double kstarwidth = 0.0474;
  
  double kstarMassMin = 0.87;
  double kstarMassMax = 0.90;
  
  double kstarwidthMin = 0.0434;
  double kstarwidthMax = 0.046;
  
  double ExtraMinRange = 0.55;
  double kstarMinRange = 0.80;
  double kstarMaxRange = 1.02;
  double BgMaxRange = 1.25;

  for (int k=0; k<nPtBins; k++)
    {
      hUSS_Y_kstarSig[k]->Add(hUSS_Y_MC_SUB[k]);
      hUSS_Y_Bg[k]->Add(hUSS_Y_MC_SUB[k]);
      hUSS_Y_Fit[k]->Add(hUSS_Y_MC_SUB[k]);

      ptlow = xBins[k];
      pthigh = xBins[k+1];
      // ==============
      // Function
      // ==============
      fBg[k] = new TF1(Form("kstar_DoubleBg_%i",k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4)", ExtraMinRange, BgMaxRange);
      fkstar[k] = new TF1(Form("kstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
      fFit[k] = new TF1(Form("Fullfit_%i", k), "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + pol4(4) + [9]*BreitWignerRelativistic(x, [10], [11])", ExtraMinRange, BgMaxRange);
	
      fMCkstar[k] = new TF1(Form("MCkstar_%i",k), "[0]*BreitWignerRelativistic(x, [1], [2])", kstarMinRange, kstarMaxRange);
      fMCkstar[k]->SetParLimits(1, kstarMassMin, kstarMassMax);
      fMCkstar[k]->SetParLimits(2, kstarwidthMin, kstarwidthMax);
      hUSS_Y_MCTrue[k]->Fit(fMCkstar[k], "SR0", "", ExtraMinRange, BgMaxRange);

      fkstar[k]->FixParameter(1, fMCkstar[k]->GetParameter(1));
      fkstar[k]->FixParameter(2, fMCkstar[k]->GetParameter(2));	  	
	
      if(k==0){
	/* for (int bin = hUSS_Y_Bg[k]->FindBin(0.87); bin <= hUSS_Y_Bg[k]->FindBin(0.93); bin++) */
	/*   hUSS_Y_Bg[k]->SetBinError(bin, 0); */
      
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
	fFit[k]->FixParameter(10, fkstar[k]->GetParameter(1));
	fFit[k]->FixParameter(11, fkstar[k]->GetParameter(2));
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
	fFit[k]->SetParameter(9, fkstar[k]->GetParameter(0));
	fFit[k]->FixParameter(10, fkstar[k]->GetParameter(1));
	fFit[k]->FixParameter(11, fkstar[k]->GetParameter(2));
      }
      hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);
      hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);

      // ================================
      // for fitting each one
      // ================================
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

      TFitResultPtr ptr = hUSS_Y_Fit[k]->Fit(fFit[k], "SR0", "", ExtraMinRange, BgMaxRange);
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
      hUSS_Y_Fit[k]->GetYaxis()->SetRangeUser(0.0001, hUSS_Y_Fit[k]->GetMaximum()*1.1);
      hUSS_Y_Fit[k]->Draw();

      
      // True
      hUSS_Y_Fit[k]->SetLineColor(1);
      hUSS_Y_MCTrue[k]->Draw("SAME");

      // MC
      fMCkstar[k]->SetLineColor(1);
      fMCkstar[k]->SetLineWidth(4);
      fMCkstar[k]->Draw("SAME");

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
      TLatex* txtInfoBB = DrawTLatex(0.56, 0.88, Kind+" yield", 0.04, 62, 1);
      TLatex* txtInfoCC = DrawTLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
      if (SpecialLatex_fit != nullptr) { SpecialLatex_fit->Draw("SAME"); }
      if (SpecialLatex_fit2 != nullptr) { SpecialLatex_fit2->Draw("SAME"); }
      TLatex* txtInfoDD = DrawTLatex(0.56, 0.74, Form("%.1f < #it{p}_{T}^{ K*}< %.1f GeV/#it{c}", ptlow, pthigh), 0.04, 42, 1);
      TLatex* txtInfoEE = DrawTLatex(0.56, 0.68, chi2NdfText, 0.04, 42, 1);

      TLegend *legendd = new TLegend(0.55,0.54,0.66,0.64);
      legendd->SetFillColor(4000);
      legendd->SetBorderSize(0);
      legendd->SetTextSize(0.03);
      legendd->AddEntry(hUSS_Y_Fit[k],"MC", "lep");
      legendd->AddEntry(fBg_POST[k], "Bkg. MCFit", "l");
      legendd->AddEntry(fkstar_POST[k], "BreitWignerRelativistic. MCFit", "l");
      legendd->AddEntry(fFit[k], "Total. Fit", "l");
      legendd->Draw("SAME");
	
      lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
      // lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.pdf", ptlow, pthigh));
      cout << "===============================" << endl;
      cout << "nth :"<< k << endl; 
    } //for

}//DrawBW
