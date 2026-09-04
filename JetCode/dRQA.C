#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TPaveStats.h>
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

TLatex* DrawTLatex(double x, double y, const TString& text, double textSize, int textFont, int textColor) {
  TLatex* txtInfo = new TLatex(x, y, text);
  txtInfo->SetNDC(kTRUE);
  txtInfo->SetTextSize(textSize);
  txtInfo->SetTextFont(textFont);
  txtInfo->SetTextColor(textColor);
  txtInfo->Draw("SAME");
  return txtInfo;
}

// Draw function
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

void DrawQA()
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString mainfile_PYTHIA = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_PYTHIA.root";
  TString mainfile_HERWIG = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_HERWIG.root";
  TString mainfile_EPOS = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_EPOS.root";
  
  TFile* data_PYTHIA = TFile::Open(mainfile_PYTHIA);
  if(!data_PYTHIA)   return;
  TFile* data_HERWIG = TFile::Open(mainfile_HERWIG);
  if(!data_HERWIG)   return;
  TFile* data_EPOS = TFile::Open(mainfile_EPOS);
  if(!data_EPOS)   return;

  TString dR_jet_daughters("kstar-in-o-o/dR_taggedjet_all;1");
  TString dR_jet_daughters_6_8("kstar-in-o-o/dR_taggedjet_all_6_8;1");

  TString recoveredJet("kstar-in-o-o/recoveredJetpT_6_8to8;1");

  TString origin_6_8("kstar-in-o-o/Kstar_pT_INJet_6_8;1");
  TString origin_8_10("kstar-in-o-o/Kstar_pT_INJet_8_10;1");
  TString origin_10_12("kstar-in-o-o/Kstar_pT_INJet_10_12;1");
  TString recovered_8_12("kstar-in-o-o/Kstar_pT_6_8to8_IN_recoveredJet;1");
  //------------------------------------------------
  // dR Study

  // jet pT > 8 GeV
  
  TH2D* h_dR_kpi_all = (TH2D*)data_PYTHIA->Get(dR_jet_daughters);
  TCanvas* dRall = new TCanvas("dRall", "", 880, 700);
  dRall->SetLeftMargin(0.13);
  dRall->SetBottomMargin(0.12);
  dRall->SetRightMargin(0.14);
  dRall->SetTopMargin(0.06);

  HistoAxisTitles(h_dR_kpi_all, "dR_{ jet to kpi}", "#it{p}_{T}^{ K^{*0}} (GeV/#it{c})",
		  0.045, 42, 1.10,
		  0.045, 42, 1.10);
  
  h_dR_kpi_all->SetTitle("");
  h_dR_kpi_all->Draw("COLZ");

  TLatex* txt_A = DrawTLatex(0.54, 0.885, "PYTHIA", 0.04, 62, 1);
  TLatex* txt_B = DrawTLatex(0.54, 0.830, "pp, #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
  TLatex* txt_C = DrawTLatex(0.54, 0.770, "K*(892)^{0} #rightarrow K^{+}#pi^{-}", 0.04, 42, 1);
  TLatex* txt_D = DrawTLatex(0.54, 0.718, "|#eta_{track}^{ch}| < 0.8,", 0.04, 42, 1);
  TLatex* txt_E = DrawTLatex(0.71, 0.718, "|#eta_{jet}^{ch} | < 0.5", 0.04, 42, 1);
  TLatex* txt_F = DrawTLatex(0.54, 0.658, "#it{p}_{T, jet}^{ ch} > 8.0 GeV/#it{c}", 0.04, 42, 1);
  TLatex* txt_G = DrawTLatex(0.54, 0.602, "Fastjet anti-k_{T} R=0.4", 0.04, 42, 1);

  // 6 < jet pT < 8 GeV/c
  
  TH2D* h_dR_kpi_6_8 = (TH2D*)data_PYTHIA->Get(dR_jet_daughters_6_8);
  TCanvas* dR_6_8 = new TCanvas("dR_6_8", "", 880, 700);
  dR_6_8->SetLeftMargin(0.13);
  dR_6_8->SetBottomMargin(0.12);
  dR_6_8->SetRightMargin(0.14);
  dR_6_8->SetTopMargin(0.06);

  HistoAxisTitles(h_dR_kpi_6_8, "dR_{ jet to kpi}", "#it{p}_{T}^{ K^{*0}} (GeV/#it{c})",
		  0.045, 42, 1.10,
		  0.045, 42, 1.10);
  
  h_dR_kpi_6_8->SetTitle("");
  h_dR_kpi_6_8->Draw("COLZ");

  TLatex* txt_AA = DrawTLatex(0.54, 0.885, "PYTHIA", 0.04, 62, 1);
  TLatex* txt_BB = DrawTLatex(0.54, 0.830, "pp, #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
  TLatex* txt_CC = DrawTLatex(0.54, 0.770, "K*(892)^{0} #rightarrow K^{+}#pi^{-}", 0.04, 42, 1);
  TLatex* txt_DD = DrawTLatex(0.54, 0.718, "|#eta_{track}^{ch}| < 0.8,", 0.04, 42, 1);
  TLatex* txt_EE = DrawTLatex(0.71, 0.718, "|#eta_{jet}^{ch} | < 0.5", 0.04, 42, 1);
  TLatex* txt_FF = DrawTLatex(0.54, 0.658, "6.0 < #it{p}_{T, jet}^{ ch} < 8.0 GeV/#it{c}", 0.04, 42, 1);
  TLatex* txt_GG = DrawTLatex(0.54, 0.602, "Fastjet anti-k_{T} R=0.4", 0.04, 42, 1);

  // recovered jet

  TH1D* h_recoveredJet = (TH1D*)data_PYTHIA->Get(recoveredJet);
  HistoAxisTitles(h_recoveredJet, "#it{p}_{T, Recovered jet} (GeV/#it{c})", "Counts",
		  0.045, 42, 1.10,
		  0.045, 42, 1.10);
  h_recoveredJet->GetXaxis()->SetRangeUser(0., 20.);
  h_recoveredJet->SetTitle("");
  h_recoveredJet->Draw();
  
  // N_K*0 in jet pT > 8 GeV/c

  TH1D* h_kstar_8_10_PY = (TH1D*)data_PYTHIA->Get(origin_8_10);
  TH1D* h_kstar_10_12_PY = (TH1D*)data_PYTHIA->Get(origin_10_12);
  TH1D* h_kstar_8_12_PY = (TH1D*)h_kstar_8_10_PY->Clone("8_10_origin_PY");
  h_kstar_8_12_PY->Add(h_kstar_10_12_PY);

  // gStyle->SetOptStat(1110);
  TCanvas* c_8_12 = new TCanvas("c", "", 880, 700);
  
  h_kstar_8_12_PY->SetStats(1);
  
  h_kstar_8_12_PY->SetTitle("");
  HistoAxisTitles(h_kstar_8_12_PY, "K^{*0} #it{p}_{T} (GeV/#it{c})","Counts",
		  0.045, 42, 1.0,
		  0.045, 42, 1.05);
  
  h_kstar_8_12_PY->GetXaxis()->SetRangeUser(0., 20.);
  h_kstar_8_12_PY->Draw();
  TLatex* txt_AAA = DrawTLatex(0.53, 0.829, "N_{ K^{*0}} in 8 < #it{p}_{T, jet} < 12 (GeV/#it{c})", 0.04, 42, 1);
  /* TPaveStats* st = (TPaveStats*)h_kstar_8_12_PY->FindObject("stats"); */
  /* if (st) { */
  /*   st->SetX1NDC(0.55); */
  /*   st->SetY1NDC(0.65); */
  /*   st->SetX2NDC(0.87); */
  /*   st->SetY2NDC(0.88); */

  /*   st->SetTextSize(0.045); */
  /*   st->SetBorderSize(1); */
  /* } */
  /* c_8_12->Modified(); */
  /* c_8_12->Update(); */

  // N_K*0 in recovered jet pT
  
  TH1D* h_recovered_kstar_8_12_PY = (TH1D*)data_PYTHIA->Get(recovered_8_12);

  TCanvas* c_recovered = new TCanvas("c", "", 880, 700);
  
  h_recovered_kstar_8_12_PY->SetStats(1);
  
  h_recovered_kstar_8_12_PY->SetTitle("");
  HistoAxisTitles(h_recovered_kstar_8_12_PY, "K^{*0} #it{p}_{T} (GeV/#it{c})","Counts",
		  0.045, 42, 1.0,
		  0.045, 42, 1.05);
  
  h_recovered_kstar_8_12_PY->GetXaxis()->SetRangeUser(0., 20.);
  h_recovered_kstar_8_12_PY->Draw();
  TLatex* txt_AAAA = DrawTLatex(0.53, 0.829, "N_{ K^{*0}} in shifted jet", 0.04, 42, 1);
  // TPaveStats* st2 = (TPaveStats*)h_recovered_kstar_8_12_PY->FindObject("stats");
  /* if (st) { */
  /*   st2->SetX1NDC(0.55); */
  /*   st2->SetY1NDC(0.65); */
  /*   st2->SetX2NDC(0.87); */
  /*   st2->SetY2NDC(0.88); */

  /*   st2->SetTextSize(0.045); */
  /*   st2->SetBorderSize(1); */
  /* } */
  /* c_recovered->Modified(); */
  /* c_recovered->Update(); */
  
  //------------------------------------------------
  // PYTHIA Corr
  // for correction
  TH1D* h_kstar_8_12_PY_corr = (TH1D*)h_kstar_8_12_PY->Clone("h_kstar_8_12_PY_corr");
  TH1D* h_kstar_8_12_PY_copy = (TH1D*)h_kstar_8_12_PY->Clone("h_kstar_8_12_PY_copy");
  
  // denominator
  TH1D* denominator_PY = (TH1D*)h_kstar_8_12_PY->Clone("8_12_Truth_deno_PY");
  denominator_PY->Add(h_recovered_kstar_8_12_PY);

  
  double A_PY = h_recovered_kstar_8_12_PY->Integral();
  double B_PY = h_kstar_8_12_PY->Integral();
    
  cout << "-----------------------------------------" << endl;
  cout << "           PYTHIA           " << endl;
  cout << "A(nK*0 in [8,12]->[6,8]): " << A_PY << endl;
  cout << "B(nK*0 in [8,12]): " << B_PY << endl;
  cout << "A+B: " << A_PY+B_PY << endl;
  cout << "X(%): " << A_PY/(A_PY+B_PY) << endl;
  cout << "*****************************************" << endl;

  double X_PY = A_PY/(A_PY+B_PY);
  cout << "Scaled B: (1-X)(A+B) : " << (1-X_PY)*(A_PY+B_PY) << endl; 
  
  // numerator

  h_kstar_8_12_PY->Scale(1.0/(1.0-X_PY));
  TH1D* numerator_PY = (TH1D*)h_kstar_8_12_PY->Clone("8_12_Truth_numerator_PY");

  cout << "Scaled A+B: " << numerator_PY->Integral() << endl;
  
  //------------------------------------------------
  // HERWIG Corr

  // denominator
  
  TH1D* h_kstar_8_10_HW = (TH1D*)data_HERWIG->Get(origin_8_10);
  TH1D* h_kstar_10_12_HW = (TH1D*)data_HERWIG->Get(origin_10_12);
  TH1D* h_kstar_8_12_HW = (TH1D*)h_kstar_8_10_HW->Clone("8_10_origin_HW");
  h_kstar_8_12_HW->Add(h_kstar_10_12_HW);

  TH1D* h_recovered_kstar_8_12_HW = (TH1D*)data_HERWIG->Get(recovered_8_12);

  // for correction
  TH1D* h_kstar_8_12_HW_corr = (TH1D*)h_kstar_8_12_HW->Clone("h_kstar_8_12_HW_corr");
  TH1D* h_kstar_8_12_HW_copy = (TH1D*)h_kstar_8_12_HW->Clone("h_kstar_8_12_HW_copy");
  

  TH1D* denominator_HW = (TH1D*)h_kstar_8_12_HW->Clone("8_12_Truth_deno_HW");
  denominator_HW->Add(h_recovered_kstar_8_12_HW);

  double A_HW = h_recovered_kstar_8_12_HW->Integral();
  double B_HW = h_kstar_8_12_HW->Integral();
    
  cout << "-----------------------------------------" << endl;
  cout << "           HERWIG           " << endl;
  cout << "A(nK*0 in [8,12]->[6,8]): " << A_HW << endl;
  cout << "B(nK*0 in [8,12]): " << B_HW << endl;
  cout << "A+B: " << A_HW+B_HW << endl;
  cout << "X(%): " << A_HW/(A_HW+B_HW) << endl;
  cout << "*****************************************" << endl;

  double X_HW = A_HW/(A_HW+B_HW);
  cout << "Scaled B: (1-X)(A+B) : " << (1-X_HW)*(A_HW+B_HW) << endl; 
  
  // numerator

  h_kstar_8_12_HW->Scale(1.0/(1.0-X_HW));
  TH1D* numerator_HW = (TH1D*)h_kstar_8_12_HW->Clone("8_12_Truth_nume_HW");

  cout << "Scaled A+B: " << numerator_HW->Integral() << endl;
  
  //------------------------------------------------
  // EPOS Corr

  // denominator
  
  TH1D* h_kstar_8_10_EPOS = (TH1D*)data_EPOS->Get(origin_8_10);
  TH1D* h_kstar_10_12_EPOS = (TH1D*)data_EPOS->Get(origin_10_12);
  TH1D* h_kstar_8_12_EPOS = (TH1D*)h_kstar_8_10_EPOS->Clone("8_10_origin_EPOS");
  h_kstar_8_12_EPOS->Add(h_kstar_10_12_EPOS);
  
  TH1D* h_recovered_kstar_8_12_EPOS = (TH1D*)data_EPOS->Get(recovered_8_12);

  // for correction
  TH1D* h_kstar_8_12_EPOS_corr = (TH1D*)h_kstar_8_12_EPOS->Clone("h_kstar_8_12_EPOS_corr");
  TH1D* h_kstar_8_12_EPOS_copy = (TH1D*)h_kstar_8_12_EPOS->Clone("h_kstar_8_12_EPOS_copy");

  TH1D* denominator_EPOS = (TH1D*)h_kstar_8_12_EPOS->Clone("8_12_Truth_deno_EPOS");
  denominator_EPOS->Add(h_recovered_kstar_8_12_EPOS);

  double A_EPOS = h_recovered_kstar_8_12_EPOS->Integral();
  double B_EPOS = h_kstar_8_12_EPOS->Integral();
    
  cout << "-----------------------------------------" << endl;
  cout << "           EPOS           " << endl;
  cout << "A(nK*0 in [8,12]->[6,8]): " << A_EPOS << endl;
  cout << "B(nK*0 in [8,12]): " << B_EPOS << endl;
  cout << "A+B: " << A_EPOS+B_EPOS << endl;
  cout << "X(%): " << A_EPOS/(A_EPOS+B_EPOS) << endl;
  cout << "*****************************************" << endl;

  double X_EPOS = A_EPOS/(A_EPOS+B_EPOS);
  cout << "Scaled B: (1-X)(A+B) : " << (1-X_EPOS)*(A_EPOS+B_EPOS) << endl; 
  
  // numerator

  h_kstar_8_12_EPOS->Scale(1.0/(1.0-X_EPOS));
  TH1D* numerator_EPOS = (TH1D*)h_kstar_8_12_EPOS->Clone("8_12_Truth_nume_EPOS");

  cout << "Scaled A+B: " << numerator_EPOS->Integral() << endl;
  cout << "-----------------------------------------" << endl;

  //================================================
  // Common binning
  //================================================
  TH1D* clone_denominator_PY = (TH1D*)denominator_PY->Clone("corr_test_PY");
  TH1D* clone_denominator_HW = (TH1D*)denominator_HW->Clone("corr_test_HW");
  TH1D* clone_denominator_EPOS = (TH1D*)denominator_EPOS->Clone("corr_test_EPOS");
  
  numerator_PY->Rebin(4);
  denominator_PY->Rebin(4);

  // Make ratios
  numerator_PY->Divide(denominator_PY);
  numerator_HW->Divide(denominator_HW);
  numerator_EPOS->Divide(denominator_EPOS);

  TCanvas* p = new TCanvas("p", "PYTHIA recovery correction", 800, 700);
  p->SetLeftMargin(0.13);
  p->SetBottomMargin(0.12);
  p->SetRightMargin(0.04);
  p->SetTopMargin(0.06);
  numerator_PY->Draw();

  TLegend* leg1 = new TLegend(0.17, 0.81, 0.43, 0.90);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.04);

  leg1->AddEntry(numerator_PY,   "PYTHIA", "lep");

  leg1->Draw("SAME");

  TCanvas* models = new TCanvas("models", "Correction comparison", 800, 700);

  models->SetLeftMargin(0.13);
  models->SetBottomMargin(0.12);
  models->SetRightMargin(0.04);
  models->SetTopMargin(0.06);

  // PYTHIA
  numerator_PY->SetLineColor(kBlack);
  numerator_PY->SetMarkerColor(kBlack);
  numerator_PY->SetMarkerStyle(20);
  numerator_PY->SetMarkerSize(1.0);
  numerator_PY->SetLineWidth(2);

  // HERWIG
  numerator_HW->SetLineColor(kRed+1);
  numerator_HW->SetMarkerColor(kRed+1);
  numerator_HW->SetMarkerStyle(21);
  numerator_HW->SetMarkerSize(1.0);
  numerator_HW->SetLineWidth(2);

  // EPOS
  numerator_EPOS->SetLineColor(kBlue+1);
  numerator_EPOS->SetMarkerColor(kBlue+1);
  numerator_EPOS->SetMarkerStyle(33);
  numerator_EPOS->SetMarkerSize(1.3);
  numerator_EPOS->SetLineWidth(2);

  HistoAxisTitles(numerator_PY,"#it{p}_{T}^{ K^{*0}} (GeV/#it{c})","Recovery ratio",
		  0.045, 42, 1.10,
		  0.045, 42, 1.25);
  
  numerator_PY->SetTitle("");
  numerator_PY->GetXaxis()->SetRangeUser(0.0, 10.0);
  
  numerator_PY->Draw();
  numerator_HW->Draw("SAME");
  numerator_EPOS->Draw("SAME");

  TLegend* leg = new TLegend(0.17, 0.75, 0.43, 0.91);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(numerator_PY,   "PYTHIA", "lep");
  leg->AddEntry(numerator_HW,   "HERWIG", "lep");
  leg->AddEntry(numerator_EPOS, "EPOS",   "lep");

  leg->Draw("SAME");
  models->SaveAs("/Users/jimun/Analysis/kstar-in-jets/Plot/kstarInjets/MC/Correction/model_dependency.png", "RECREATE");

  //=============================
  //
  //   Fit through all models
  //
  //============================
  TH1D* numerator_PY_common = (TH1D*)numerator_HW->Clone("numerator_PY_common");

  numerator_PY_common->Reset("ICES");

  for (int i = 1; i <= numerator_PY_common->GetNbinsX(); i++) {
    double pt = numerator_PY_common->GetBinCenter(i);
    int binPY = numerator_PY->FindBin(pt);
    
    numerator_PY_common->SetBinContent(i, numerator_PY->GetBinContent(binPY));
    numerator_PY_common->SetBinError(i,numerator_PY->GetBinError(binPY));
  }

 
  TH1D* h_Ratio_Avg = (TH1D*)numerator_PY_common->Clone("h_Ratio_Avg");
  h_Ratio_Avg->Add(numerator_HW);
  h_Ratio_Avg->Add(numerator_EPOS);

  h_Ratio_Avg->Scale(1.0/3.0);
  h_Ratio_Avg->SetLineColor(kMagenta+2);
 
  HistoAxisTitles(h_Ratio_Avg, "#it{p}_{T}^{ K^{*0}} (GeV/#it{c})", "Average recovery ratio",
		  0.045, 42, 1.10, 0.045, 42, 1.25);
  
  TCanvas* c_avg = new TCanvas("c_avg", "", 800, 700);
  c_avg->SetLeftMargin(0.13);
  c_avg->SetBottomMargin(0.12);
  c_avg->SetRightMargin(0.04);
  c_avg->SetTopMargin(0.06);
  
  h_Ratio_Avg->SetTitle("");
  h_Ratio_Avg->GetXaxis()->SetRangeUser(0., 6.);
  h_Ratio_Avg->Draw();


  TF1* fgaus = new TF1("fgaus", "[0]-gaus(1)", 0, 3.0);
  fgaus->SetParameter(0, 1.0);
  fgaus->SetParameter(1, 0.12);
  fgaus->SetParameter(2, 1.50);
  fgaus->SetParameter(3, 0.80);

  TF1* flogi = new TF1("flogi", "[0]+[1]/(1.0 + exp(-(x-[2])/[3]))", 3.0, 6.0);
  flogi->SetParameter(0, 0.95);
  flogi->SetParameter(1, 0.30);
  flogi->SetParameter(2, 3.8);
  flogi->SetParameter(3, 0.6);
  
  TF1* fCorr = new TF1("fCorr", "[0]-gaus(1) + [4]+[5]/(1.0 + exp(-(x-[6])/[7]))", 0.0, 6.0);
  fCorr->SetParameter(0, fgaus->GetParameter(0));
  fCorr->SetParameter(1, fgaus->GetParameter(1));
  fCorr->SetParameter(2, fgaus->GetParameter(2));
  fCorr->SetParameter(3, fgaus->GetParameter(3));
  fCorr->SetParameter(4, flogi->GetParameter(0));
  fCorr->SetParameter(5, flogi->GetParameter(1));
  fCorr->SetParameter(6, flogi->GetParameter(2));
  fCorr->SetParameter(7, flogi->GetParameter(3));
  
  h_Ratio_Avg->Fit(fCorr, "R", "", 0.0, 6.0); //USE the fit value i gave the info.
  h_Ratio_Avg->Draw();
  fCorr->Draw("SAME");

  TLegend* leg_fit = new TLegend(0.16, 0.82, 0.41, 0.93);
  leg_fit->SetBorderSize(0);
  leg_fit->SetFillStyle(0);
  leg_fit->SetTextFont(42);
  leg_fit->SetTextSize(0.04);

  leg_fit->AddEntry(h_Ratio_Avg, "Average value from models", "lep");
  leg_fit->Draw("SAME");

  //============================
  //
  //      Correction
  //
  //============================
  double X_avg = (X_PY+X_HW+X_EPOS)/3.0;
  cout << "X_avg: " << X_avg << endl;

  //==========
  // PYTHIA
  //==========
  double corr_global_PY = 1.0/(1.0-X_PY);
  for (int i=1; i<=h_kstar_8_12_PY_corr->GetNbinsX(); i++){
    
    double pt = h_kstar_8_12_PY_corr->GetBinCenter(i);
    double yield = h_kstar_8_12_PY_copy->GetBinContent(i);
    double error = h_kstar_8_12_PY_copy->GetBinError(i);

    double corr_tot = corr_global_PY;

    if(pt >=0.0 && pt <= 6.0){
      double closure = fCorr->Eval(pt);
      if (closure <= 0.0)
	continue;

      double corr_shape = 1.0/closure;
      corr_tot *=corr_shape;
    }
    
    h_kstar_8_12_PY_corr->SetBinContent(i, yield*corr_tot);
    h_kstar_8_12_PY_corr->SetBinError(i, error*corr_tot);
  }
  
  TCanvas* c_corr = new TCanvas("c_corr", "", 800, 700);
  TPad* p1_corr = new TPad("PYTHIA_p1", "", 0, 0.3, 1, 1);
  p1_corr->SetTopMargin(0.05);
  p1_corr->SetLeftMargin(0.15);
  p1_corr->Draw();
  p1_corr->cd();
  p1_corr->SetBottomMargin(0.01);

  h_kstar_8_12_PY_copy->Rebin(4);
  h_kstar_8_12_PY_corr->Rebin(4);
  clone_denominator_PY->Rebin(4);

  h_kstar_8_12_PY_copy->Scale(1.0/4.0);
  h_kstar_8_12_PY_corr->Scale(1.0/4.0);
  clone_denominator_PY->Scale(1.0/4.0);
  
  h_kstar_8_12_PY_copy->SetTitle("");
  clone_denominator_PY->SetTitle("");
  
  h_kstar_8_12_PY_copy->GetXaxis()->SetRangeUser(0., 6.);
  h_kstar_8_12_PY_corr->GetXaxis()->SetRangeUser(0., 6.);
  clone_denominator_PY->GetXaxis()->SetRangeUser(0., 6.);

  clone_denominator_PY->GetYaxis()->SetRangeUser(0., 680);

  clone_denominator_PY->SetLineColor(kBlue);
  clone_denominator_PY->SetMarkerColor(kBlue);
  clone_denominator_PY->SetMarkerStyle(24);
  
  h_kstar_8_12_PY_copy->SetLineColor(kBlack);
  h_kstar_8_12_PY_copy->SetMarkerColor(kBlack);
  h_kstar_8_12_PY_copy->SetMarkerStyle(20);

  h_kstar_8_12_PY_corr->SetLineColor(kRed+1);
  h_kstar_8_12_PY_corr->SetMarkerColor(kRed+1);
  h_kstar_8_12_PY_corr->SetMarkerStyle(21);

  clone_denominator_PY->Draw();
  h_kstar_8_12_PY_copy->Draw("SAME");
  h_kstar_8_12_PY_corr->Draw("SAME");
  
  TLatex* txt = DrawTLatex(0.61, 0.87, "PYTHIA", 0.06, 62, 1);
  TLegend* leg_corr = new TLegend(0.60, 0.67, 0.86, 0.84);
  leg_corr->SetBorderSize(0);
  leg_corr->SetFillStyle(0);
  leg_corr->SetTextFont(42);
  leg_corr->SetTextSize(0.04);
  
  leg_corr->AddEntry(h_kstar_8_12_PY_copy, "Measured K^{*0}", "lep");
  leg_corr->AddEntry(clone_denominator_PY, "Inclusive K^{*0}", "lep");
  leg_corr->AddEntry(h_kstar_8_12_PY_corr, "Migration-corrected K^{*0}", "lep");
 
  leg_corr->Draw("SAME");

  c_corr->cd();
  TPad* p2_corr = new TPad("PYTHIA_p2", "", 0, 0, 1, 0.3);
  p2_corr->SetLeftMargin(0.15);
  p2_corr->SetBottomMargin(0.3);
  p2_corr->SetTopMargin(0.0001);
  p2_corr->Draw();
  p2_corr->cd();

  TH1D* h_ratio_corr = (TH1D*)h_kstar_8_12_PY_corr->Clone("h_ratio_corr");

  HistoAxisTitles(h_ratio_corr, "K^{*0} #it{p}_{T} (GeV/#it{c})", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  h_ratio_corr->GetXaxis()->SetRangeUser(0., 6.);
  h_ratio_corr->GetXaxis()->SetRangeUser(0., 6.);
  
  h_ratio_corr->SetTitle("");
  h_ratio_corr->GetXaxis()->SetLabelSize(0.09);
  h_ratio_corr->GetYaxis()->SetLabelSize(0.08);

  h_ratio_corr->SetLineColor(4);

  h_ratio_corr->Divide(clone_denominator_PY);
  h_ratio_corr->Draw("SAME");

  TLegend* leg_below = new TLegend(0.52, 0.35, 0.78, 0.45);
  leg_below->SetBorderSize(0);
  leg_below->SetFillStyle(0);
  leg_below->SetTextFont(42);
  leg_below->SetTextSize(0.1);

  leg_below->AddEntry(h_ratio_corr, "Migration-corrected / Inclusive", "lep");
  leg_below->Draw("SAME");

  
  //====================
  // HERWIG correction
  //====================
 
  double corr_global_HW = 1.0/(1.0-X_HW);
  for (int i=1; i<=h_kstar_8_12_HW_corr->GetNbinsX(); i++){
    
    double pt_HW = h_kstar_8_12_HW_corr->GetBinCenter(i);
    double yield_HW = h_kstar_8_12_HW_copy->GetBinContent(i);
    double error_HW = h_kstar_8_12_HW_copy->GetBinError(i);

    double corr_tot_HW = corr_global_HW;

    if(pt_HW >=0.0 && pt_HW <= 6.0){
      double closure_HW = fCorr->Eval(pt_HW);
      if (closure_HW <= 0.0)
	continue;

      double corr_shape_HW = 1.0/closure_HW;
      corr_tot_HW *=corr_shape_HW;
    }
    
    h_kstar_8_12_HW_corr->SetBinContent(i, yield_HW*corr_tot_HW);
    h_kstar_8_12_HW_corr->SetBinError(i, error_HW*corr_tot_HW);
  }
  
  TCanvas* c_corr_HW = new TCanvas("c_corr_HW", "", 800, 700);
  TPad* p1_corr_HW = new TPad("HERWIG_p1", "", 0, 0.3, 1, 1);
  p1_corr_HW->SetTopMargin(0.05);
  p1_corr_HW->SetLeftMargin(0.15);
  p1_corr_HW->Draw();
  p1_corr_HW->cd();
  p1_corr_HW->SetBottomMargin(0.01);

  h_kstar_8_12_HW_copy->Scale(1.0);
  h_kstar_8_12_HW_corr->Scale(1.0);
  clone_denominator_HW->Scale(1.0);
  
  h_kstar_8_12_HW_copy->SetTitle("");
  clone_denominator_HW->SetTitle("");
  
  h_kstar_8_12_HW_copy->GetXaxis()->SetRangeUser(0., 6.);
  h_kstar_8_12_HW_corr->GetXaxis()->SetRangeUser(0., 6.);
  clone_denominator_HW->GetXaxis()->SetRangeUser(0., 6.);

  clone_denominator_HW->GetYaxis()->SetRangeUser(0., 7500);
  HistoAxisTitles(clone_denominator_HW, "", "Counts",
		  0.045, 42, 1.0, 0.045, 42, 1.2);

  clone_denominator_HW->SetLineColor(kBlue);
  clone_denominator_HW->SetMarkerColor(kBlue);
  clone_denominator_HW->SetMarkerStyle(24);
  
  h_kstar_8_12_HW_copy->SetLineColor(kBlack);
  h_kstar_8_12_HW_copy->SetMarkerColor(kBlack);
  h_kstar_8_12_HW_copy->SetMarkerStyle(20);

  h_kstar_8_12_HW_corr->SetLineColor(kRed+1);
  h_kstar_8_12_HW_corr->SetMarkerColor(kRed+1);
  h_kstar_8_12_HW_corr->SetMarkerStyle(21);

  clone_denominator_HW->Draw();
  h_kstar_8_12_HW_copy->Draw("SAME");
  h_kstar_8_12_HW_corr->Draw("SAME");
  
  TLatex* txt_HW = DrawTLatex(0.61, 0.87, "HERWIG", 0.06, 62, 1);
  TLegend* leg_corr_HW = new TLegend(0.60, 0.67, 0.86, 0.84);
  leg_corr_HW->SetBorderSize(0);
  leg_corr_HW->SetFillStyle(0);
  leg_corr_HW->SetTextFont(42);
  leg_corr_HW->SetTextSize(0.04);
  
  leg_corr_HW->AddEntry(h_kstar_8_12_HW_copy, "Measured K^{*0}", "lep");
  leg_corr_HW->AddEntry(clone_denominator_HW, "Inclusive K^{*0}", "lep");
  leg_corr_HW->AddEntry(h_kstar_8_12_HW_corr, "Migration-corrected K^{*0}", "lep");
 
  leg_corr_HW->Draw("SAME");

  c_corr_HW->cd();
  TPad* p2_corr_HW = new TPad("HERWIG_p2", "", 0, 0, 1, 0.3);
  p2_corr_HW->SetLeftMargin(0.15);
  p2_corr_HW->SetBottomMargin(0.3);
  p2_corr_HW->SetTopMargin(0.0001);
  p2_corr_HW->Draw();
  p2_corr_HW->cd();

  TH1D* h_ratio_corr_HW = (TH1D*)h_kstar_8_12_HW_corr->Clone("h_ratio_corr_HW");

  HistoAxisTitles(h_ratio_corr_HW, "K^{*0} #it{p}_{T} (GeV/#it{c})", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  h_ratio_corr_HW->GetXaxis()->SetRangeUser(0., 6.);
  h_ratio_corr_HW->GetXaxis()->SetRangeUser(0., 6.);
  
  h_ratio_corr_HW->SetTitle("");
  h_ratio_corr_HW->GetXaxis()->SetLabelSize(0.09);
  h_ratio_corr_HW->GetYaxis()->SetLabelSize(0.08);

  h_ratio_corr_HW->SetLineColor(4);

  h_ratio_corr_HW->Divide(clone_denominator_HW);
  h_ratio_corr_HW->Draw("SAME");

  TLegend* leg_below_HW = new TLegend(0.52, 0.35, 0.78, 0.45);
  leg_below_HW->SetBorderSize(0);
  leg_below_HW->SetFillStyle(0);
  leg_below_HW->SetTextFont(42);
  leg_below_HW->SetTextSize(0.1);

  leg_below_HW->AddEntry(h_ratio_corr_HW, "Migration-corrected / Inclusive", "lep");
  leg_below_HW->Draw("SAME");


  //====================
  // EPOS correction
  //====================
 
  double corr_global_EPOS = 1.0/(1.0-X_EPOS);
  for (int i=1; i<=h_kstar_8_12_EPOS_corr->GetNbinsX(); i++){
    
    double pt_EPOS = h_kstar_8_12_EPOS_corr->GetBinCenter(i);
    double yield_EPOS = h_kstar_8_12_EPOS_copy->GetBinContent(i);
    double error_EPOS = h_kstar_8_12_EPOS_copy->GetBinError(i);

    double corr_tot_EPOS = corr_global_EPOS;

    if(pt_EPOS >=0.0 && pt_EPOS <= 6.0){
      double closure_EPOS = fCorr->Eval(pt_EPOS);
      if (closure_EPOS <= 0.0)
	continue;

      double corr_shape_EPOS = 1.0/closure_EPOS;
      corr_tot_EPOS *=corr_shape_EPOS;
    }
    
    h_kstar_8_12_EPOS_corr->SetBinContent(i, yield_EPOS*corr_tot_EPOS);
    h_kstar_8_12_EPOS_corr->SetBinError(i, error_EPOS*corr_tot_EPOS);
  }
  
  TCanvas* c_corr_EPOS = new TCanvas("c_corr_EPOS", "", 800, 700);
  TPad* p1_corr_EPOS = new TPad("EPOS_p1", "", 0, 0.3, 1, 1);
  p1_corr_EPOS->SetTopMargin(0.05);
  p1_corr_EPOS->SetLeftMargin(0.15);
  p1_corr_EPOS->Draw();
  p1_corr_EPOS->cd();
  p1_corr_EPOS->SetBottomMargin(0.01);

  h_kstar_8_12_EPOS_copy->Scale(1.0);
  h_kstar_8_12_EPOS_corr->Scale(1.0);
  clone_denominator_EPOS->Scale(1.0);
  
  h_kstar_8_12_EPOS_copy->SetTitle("");
  clone_denominator_EPOS->SetTitle("");
  
  h_kstar_8_12_EPOS_copy->GetXaxis()->SetRangeUser(0., 6.);
  h_kstar_8_12_EPOS_corr->GetXaxis()->SetRangeUser(0., 6.);
  clone_denominator_EPOS->GetXaxis()->SetRangeUser(0., 6.);

  clone_denominator_EPOS->GetYaxis()->SetRangeUser(0., 900);
  HistoAxisTitles(clone_denominator_EPOS, "", "Counts",
		  0.045, 42, 1.0, 0.045, 42, 1.2);

  clone_denominator_EPOS->SetLineColor(kBlue);
  clone_denominator_EPOS->SetMarkerColor(kBlue);
  clone_denominator_EPOS->SetMarkerStyle(24);
  
  h_kstar_8_12_EPOS_copy->SetLineColor(kBlack);
  h_kstar_8_12_EPOS_copy->SetMarkerColor(kBlack);
  h_kstar_8_12_EPOS_copy->SetMarkerStyle(20);

  h_kstar_8_12_EPOS_corr->SetLineColor(kRed+1);
  h_kstar_8_12_EPOS_corr->SetMarkerColor(kRed+1);
  h_kstar_8_12_EPOS_corr->SetMarkerStyle(21);

  clone_denominator_EPOS->Draw();
  h_kstar_8_12_EPOS_copy->Draw("SAME");
  h_kstar_8_12_EPOS_corr->Draw("SAME");
  
  TLatex* txt_EPOS = DrawTLatex(0.61, 0.87, "EPOS", 0.06, 62, 1);
  TLegend* leg_corr_EPOS = new TLegend(0.60, 0.67, 0.86, 0.84);
  leg_corr_EPOS->SetBorderSize(0);
  leg_corr_EPOS->SetFillStyle(0);
  leg_corr_EPOS->SetTextFont(42);
  leg_corr_EPOS->SetTextSize(0.04);
  
  leg_corr_EPOS->AddEntry(h_kstar_8_12_EPOS_copy, "Measured K^{*0}", "lep");
  leg_corr_EPOS->AddEntry(clone_denominator_EPOS, "Inclusive K^{*0}", "lep");
  leg_corr_EPOS->AddEntry(h_kstar_8_12_EPOS_corr, "Migration-corrected K^{*0}", "lep");
 
  leg_corr_EPOS->Draw("SAME");

  c_corr_EPOS->cd();
  TPad* p2_corr_EPOS = new TPad("EPOS_p2", "", 0, 0, 1, 0.3);
  p2_corr_EPOS->SetLeftMargin(0.15);
  p2_corr_EPOS->SetBottomMargin(0.3);
  p2_corr_EPOS->SetTopMargin(0.0001);
  p2_corr_EPOS->Draw();
  p2_corr_EPOS->cd();

  TH1D* h_ratio_corr_EPOS = (TH1D*)h_kstar_8_12_EPOS_corr->Clone("h_ratio_corr_EPOS");

  HistoAxisTitles(h_ratio_corr_EPOS, "K^{*0} #it{p}_{T} (GeV/#it{c})", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  h_ratio_corr_EPOS->GetXaxis()->SetRangeUser(0., 6.);
  h_ratio_corr_EPOS->GetXaxis()->SetRangeUser(0., 6.);
  
  h_ratio_corr_EPOS->SetTitle("");
  h_ratio_corr_EPOS->GetXaxis()->SetLabelSize(0.09);
  h_ratio_corr_EPOS->GetYaxis()->SetLabelSize(0.08);

  h_ratio_corr_EPOS->SetLineColor(4);

  h_ratio_corr_EPOS->Divide(clone_denominator_EPOS);
  h_ratio_corr_EPOS->Draw("SAME");

  TLegend* leg_below_EPOS = new TLegend(0.52, 0.35, 0.78, 0.45);
  leg_below_EPOS->SetBorderSize(0);
  leg_below_EPOS->SetFillStyle(0);
  leg_below_EPOS->SetTextFont(42);
  leg_below_EPOS->SetTextSize(0.1);

  leg_below_EPOS->AddEntry(h_ratio_corr_EPOS, "Migration-corrected / Inclusive", "lep");
  leg_below_EPOS->Draw("SAME");
  
}
