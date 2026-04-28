#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
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

void corrected()
{  
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString localpath = "Data/Inclusive/";
  int blue3 = TColor::GetColor("#0096FF");
  int red = TColor::GetColor("#ff0400");
  int orange = TColor::GetColor("#ff9900");
  int green = TColor::GetColor("#23b028");
  int purple = TColor::GetColor("#7d20c9");  
  int black = TColor::GetColor("#2c2d33");  
  //==================
  // efficiency data
  //==================
  TString mainfile1 = "../kstar-in-jets/Results_draw/eff_Inclusive_MB.root";
  TString mainfile1_jet = "../kstar-in-jets/Results_draw/eff_Jet_MB.root";
  
  TFile* data1 = TFile::Open(mainfile1);
  TFile* data1_jet = TFile::Open(mainfile1_jet);
  
  if(!data1)   return;
  if(!data1_jet)   return;
  //==================
  // uncorrected data
  //==================
  TString MB = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_MB.root";
  TString Jet = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_MB.root";
 
  TFile* data2 = TFile::Open(MB); 
  TFile* data2_jet = TFile::Open(Jet); 

  if(!data2)   return;
  if(!data2_jet)   return;

  //==================
  // Gen data
  //==================
  TString mainfile3 = "../kstar-in-jets/Results_rootfile/Eff_MB_hyperloop.root";
  TFile* data3 = TFile::Open(mainfile3);
  if(!data3)   return;
  TString uncorrMBName("Inclusive_spectra;1");
  TString uncorrJetName("Jets_spectra;1");
  TString effName("eff_MB;1");
 
  TH1D* heff = (TH1D*)data1->Get(effName);
  TH1D* heff_jet = (TH1D*)data1_jet->Get(effName);
  
  TH1D* huncorr = (TH1D*)data2->Get(uncorrMBName);
  TH1D* huncorr_jet = (TH1D*)data2_jet->Get(uncorrJetName);
  
  TH1D* hcorr = (TH1D*)huncorr->Clone("hcorr");
  TH1D* hcorr_jet = (TH1D*)huncorr_jet->Clone("hcorr_jet");
 
  HistoAxisTitles(hcorr, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta} #frac{1}{eff x Acc} #frac{1}{BR}", 0.04, 62, 1.0, 0.04, 62, 1.7);
  TCanvas* c = new TCanvas("", "", 1020, 800);
  c->SetLeftMargin(0.16);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.14);
  c->SetLogy();

  TFile* fout = TFile::Open("/Users/jimun/Analysis/kstar-in-jets/Results_draw/corr_Inclusive_Jet_MB.root", "RECREATE");
  fout->cd();
  hcorr->Write("Uncorr_Inclusive_MB");
  hcorr_jet->Write("Uncorr_Jets_MB");

  hcorr->Divide(heff);
  hcorr_jet->Divide(heff_jet);

  hcorr->Write("Eff_corr_Inclusive_MB");
  hcorr_jet->Write("Eff_corr_Jets_MB");

  //==================
  // signal Loss
  //==================
  /* const Int_t nPtBins = 11; */
  /* Double_t xBins[nPtBins+1] = {0.0, 0.8, 1.0, 1.2, 1.4, 1.8, 2.3, 2.8, 3.4, 4.0, 5.0, 8.0}; */

  /* TString kstarGenName("kstar-in-o-o/hEffGen_pT;1"); */
  /* TString kstarRawName("kstar-in-o-o/hEffGen_pT_Raw;1"); */

  /* TH1D* hGen = (TH1D*)data3->Get(kstarGenName); */
  /* TH1D* hRaw = (TH1D*)data3->Get(kstarRawName); */
  
  /* TH1D* hGenR = new TH1D("hGenRebin", "hGenRebin", nPtBins, xBins); */
  /* TH1D* hRawR = new TH1D("hRawRebin", "hRawRebin", nPtBins, xBins); */
  
  /* hGenR = (TH1D*)hGen->Rebin(nPtBins, "hGenRebin", xBins); */
  /* hRawR = (TH1D*)hRaw->Rebin(nPtBins, "hRawRebin", xBins); */

  /* double nGenR = hGenR->Integral(); */
  /* double nRawR = hRawR->Integral(); */

  
  /* double signalLoss = nGenR/nRawR; */
  
  /* cout << "nGenR: " << nGenR << endl; */
  /* cout << "nRawR: " << nRawR << endl; */
  /* cout << "Signal Loss Ratio: " << signalLoss <<endl; */

  /* hcorr->Scale(1.0/signalLoss); */
  
  /* hcorr->Write("signal_corrected_pT"); */

  //==================
  //|| Event Loss
  // ==================
  /* TString RawEventName("kstar-in-o-o/nEvents_MC;1"); */
  /* TString GenEventName("kstar-in-o-o/nEvents_MC_True;1"); */

  /* TH1D* hEvRaw = (TH1D*)data3->Get(RawEventName); */
  /* TH1D* hEvGen = (TH1D*)data3->Get(GenEventName); */
  
  /* double numRaw = hEvRaw->GetBinContent(3); */
  /* double numGen = hEvGen->GetBinContent(1); */
  
  /* double eventLoss = numRaw/numGen; */
  
  /* cout<<"number of EvRaw: " << numRaw <<endl; */
  /* cout<<"number of EvGen: " << numGen <<endl; */
  /* cout<<"Event Loss Ratio: " << eventLoss <<endl; */

  /* hcorr->Scale(eventLoss); */
  /* hcorr->Write("event_corrected_pT"); */
  
  fout->Close();

  hcorr->SetLineWidth(2);
  hcorr->SetLineColor(blue3);
  hcorr->SetMarkerStyle(43); //circle
  hcorr->SetMarkerSize(3);
  hcorr->SetMarkerColor(blue3);
  hcorr->Draw("SAME");
  
  hcorr_jet->SetLineWidth(2);
  hcorr_jet->SetLineColor(red);
  hcorr_jet->SetMarkerStyle(89);
  hcorr_jet->SetMarkerSize(2);
  hcorr_jet->SetMarkerColor(red);
  hcorr_jet->Draw("SAME");

  TLatex* txtInfoAAA = DrawTLatex(0.555, 0.868, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.555, 0.810, "Jets, #sqrt{s} = 13.6 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.555, 0.748, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.555, 0.695, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.543, 0.589, 0.877, 0.670);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hcorr, "Inclusive", "pl"); //pl : point line
  legenddd->AddEntry(hcorr_jet, "Jet", "pl"); //pl : point line
  legenddd->Draw("SAME");
  
  c->SaveAs("../kstar-in-jets/Plot/kstarInjets/"+localpath+"correctedpT_Inclusive_Jets.png", "RECREATE");
}


