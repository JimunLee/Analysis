#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TString.h>
#include <TLegend.h>
#include <TColor.h>
#include <TStyle.h>

using namespace std;
using namespace TMath;


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


void triggers()
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString MB = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_MB.root";
  TString trigger1 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg1_JetHighPt.root";
  TString trigger2 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg2_JetLowPt.root";
  TString trigger3 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg3_TrackHighPt.root";
  TString trigger4 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg4_TrackLowPt.root";

  TFile* minbias = TFile::Open(MB);
  if(!minbias) return;
  TFile* trg1 = TFile::Open(trigger1);
  if(!trg1) return;
  TFile* trg2 = TFile::Open(trigger2);
  if(!trg2) return;
  TFile* trg3 = TFile::Open(trigger3);
  if(!trg3) return;
  TFile* trg4 = TFile::Open(trigger4);
  if(!trg4) return;
  
  TString trgName("Inclusive_spectra;1");

  TH1D* hMB = (TH1D*)minbias->Get(trgName);
  TH1D* htrg1 = (TH1D*)trg1->Get(trgName);
  TH1D* htrg2 = (TH1D*)trg2->Get(trgName);
  TH1D* htrg3 = (TH1D*)trg3->Get(trgName);
  TH1D* htrg4 = (TH1D*)trg4->Get(trgName);

  int blue3 = TColor::GetColor("#0096FF");
  int red = TColor::GetColor("#ff0400");
  int orange = TColor::GetColor("#ff9900");
  int green = TColor::GetColor("#23b028");
  int purple = TColor::GetColor("#7d20c9");  
  int black = TColor::GetColor("#2c2d33");  
  
  
  TCanvas* c = new TCanvas();
  /* TPad* p = new TPad("", "", 0, 0.3, 1, 1); */
  /* p->SetLeftMargin(0.15); */
  /* p->SetBottomMargin(0.001); */
  /* p->SetLogy(); */
  /* p->Draw(); */
  /* p->cd(); */
  
  HistoAxisTitles(hMB, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}", 0.05, 62, 1.0, 0.05, 62, 1.2);
  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  c->SetBottomMargin(0.12);

  hMB->SetLineWidth(3);
  hMB->SetLineColor(blue3);
  hMB->SetMarkerStyle(43);
  hMB->SetMarkerSize(4);
  hMB->SetMarkerColor(blue3);

  htrg1->SetLineWidth(2);
  htrg1->SetLineColor(red);
  htrg1->SetMarkerStyle(89);
  htrg1->SetMarkerSize(2);
  htrg1->SetMarkerColor(red);

  htrg2->SetLineWidth(2);
  htrg2->SetLineColor(orange);
  htrg2->SetMarkerStyle(89);
  htrg2->SetMarkerSize(2);
  htrg2->SetMarkerColor(orange);

  htrg3->SetLineWidth(2);
  htrg3->SetLineColor(green);
  htrg3->SetMarkerStyle(89);
  htrg3->SetMarkerSize(2);
  htrg3->SetMarkerColor(green);

  htrg4->SetLineWidth(2);
  htrg4->SetLineColor(purple);
  htrg4->SetMarkerStyle(89);
  htrg4->SetMarkerSize(2);
  htrg4->SetMarkerColor(purple);

  hMB->Draw("SAME");
  htrg1->Draw("SAME");
  htrg2->Draw("SAME");
  htrg3->Draw("SAME");
  htrg4->Draw("SAME");

  TLegend* legenddd = new TLegend(0.66, 0.61, 0.99, 0.86);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hMB, "Min Bias", "pl"); //pl : point line
  legenddd->AddEntry(htrg1, "Jet High Pt", "pl"); //pl : point line
  legenddd->AddEntry(htrg2, "Jet Low Pt", "pl"); //pl : point line
  legenddd->AddEntry(htrg3, "Track High Pt", "pl"); //pl : point line
  legenddd->AddEntry(htrg4, "Track Low Pt", "pl"); //pl : point line
  legenddd->Draw("SAME");

  c->SetLogy();
  c->cd();

  /* TPad* p2 = new TPad("", "", 0, 0, 1, 0.3); */
  /* p2->SetLeftMargin(0.15); */
  /* p2->SetBottomMargin(0.4); */
  /* p2->SetTopMargin(0.0001); */
  /* p2->Draw(); */
  /* p2->cd(); */

  /* TH1D* hdummy2 = new TH1D("", "", 1, 0,9); */
  /* HistoAxisTitles(hdummy2, "#it{p}_{T} [GeV/#it{c}]", "", */
  /* 		  0.1, 62, 1.4); */
  /* hdummy2->GetXaxis()->SetLabelSize(0.10); */
  /* hdummy2->GetYaxis()->SetLabelSize(0.10); */
  
  /* hdummy2->SetMinimum(0.0); */
  /* hdummy2->SetMaximum(20); */
  /* hdummy2->GetYaxis()->SetNdivisions(303); */
  /* hdummy2->Draw(); */

  /* TH1D* hRatio = (TH1D*)hskimmed->Clone("hRatio"); */
  /* hRatio->Divide((TH1D*)hMB); */
  /* hRatio->Draw("SAME"); */

  c->SaveAs("../kstar-in-jets/Plot/kstarInjets/Data/Inclusive/Uncorr_triggers_MB_skimmed.png", "RECREATE");
}

