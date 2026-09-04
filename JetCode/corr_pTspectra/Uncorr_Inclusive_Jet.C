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
#include <TLatex.h>

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

TLatex* DrawTLatex(double x, double y, const TString& text, double textSize, int textFont, int textColor) {
  TLatex* txtInfo = new TLatex(x, y, text);
  txtInfo->SetNDC(kTRUE);
  txtInfo->SetTextSize(textSize);
  txtInfo->SetTextFont(textFont);
  txtInfo->SetTextColor(textColor);
  txtInfo->Draw("SAME");
  return txtInfo;
}

void Uncorr()
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString Inclusive_PATH = "/Users/jimun/Analysis/kstar-in-jets/Results_draw/Uncorr_Inclusive.root";
  TString Jet_PATH = "/Users/jimun/Analysis/kstar-in-jets/Results_draw/Uncorr_Jet.root";
  
  TFile* inclusive = TFile::Open(Inclusive_PATH);
  if(!inclusive) return;
  TFile* jet = TFile::Open(Jet_PATH);
  if(!jet) return;
  
  TString InclusiveName("Inclusive_spectra;1");
  TString JetName("Jets_spectra;1");

  TH1D* hInclusive = (TH1D*)inclusive->Get(InclusiveName);
  TH1D* hJets = (TH1D*)jet->Get(JetName);

  int blue3 = TColor::GetColor("#0096FF");
  int red = TColor::GetColor("#ff0400");
  
  TCanvas* c = new TCanvas("", "", 940, 700);
  HistoAxisTitles(hInclusive, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta}", 0.05, 62, 1.0, 0.05, 62, 1.2);
  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  c->SetBottomMargin(0.12);

  hInclusive->SetLineWidth(3);
  hInclusive->SetLineColor(blue3);
  hInclusive->SetMarkerStyle(43);
  hInclusive->SetMarkerSize(4);
  hInclusive->SetMarkerColor(blue3);

  hJets->SetLineWidth(2);
  hJets->SetLineColor(red);
  hJets->SetMarkerStyle(89);
  hJets->SetMarkerSize(2);
  hJets->SetMarkerColor(red);
 
  hInclusive->Draw();
  hJets->Draw("SAME");

  TLatex* txt_A = DrawTLatex(0.61, 0.845, "ALICE Work In Progress", 0.04, 62, 1);
  TLatex* txt_B = DrawTLatex(0.61, 0.789, "pp, #sqrt{s} = 13.6 TeV", 0.04, 42, 1);
  TLatex* txt_C = DrawTLatex(0.61, 0.730, "K*(892)^{0} #rightarrow K^{+}#pi^{-}", 0.04, 42, 1);
  TLatex* txt_D = DrawTLatex(0.61, 0.675, "|#eta_{track}^{ch}| < 0.8,", 0.04, 42, 1);
  TLatex* txt_E = DrawTLatex(0.78, 0.675, "|#eta_{jet}^{ch} | < 0.5", 0.04, 42, 1);
  TLatex* txt_F = DrawTLatex(0.61, 0.610, "#it{p}_{T, jet}^{ ch} > 8.0 GeV/#it{c}", 0.04, 42, 1);
  TLatex* txt_G = DrawTLatex(0.61, 0.545, "Fastjet anti-k_{T} R=0.4", 0.04, 42, 1);
 
  TLegend* legenddd = new TLegend(0.18, 0.43, 0.52, 0.52);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hInclusive, "Inclusive", "pl"); //pl : point line
  legenddd->AddEntry(hJets, "In Jet", "pl"); //pl : point line
  legenddd->Draw("SAME");

  hInclusive->SetMaximum(7e-2);
  hInclusive->SetMinimum(2e-6);
  
  c->SetLogy();
  c->cd();

  c->SaveAs("/Users/jimun/Analysis/kstar-in-jets/Plot/kstarInjets/Data/Uncorr_Inclusive_Jet.png", "RECREATE");
}

