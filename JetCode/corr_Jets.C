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

  TString localpath = "Data/Inside/";
  int blue3 = TColor::GetColor("#0096FF");
  int red = TColor::GetColor("#ff0400");
  int orange = TColor::GetColor("#ff9900");
  int green = TColor::GetColor("#23b028");
  int purple = TColor::GetColor("#7d20c9");  
  int black = TColor::GetColor("#2c2d33");  
  //==================
  // efficiency data
  //==================
  TString mainfile1 = "../kstar-in-jets/Results_draw/eff_Jet_MB.root";
  TFile* data1 = TFile::Open(mainfile1);
  if(!data1)   return;

  //==================
  // uncorrected data
  //==================
  TString MB = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_MB.root";
  TString trigger1 = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_trg1_JetHighPt.root";
  TString trigger2 = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_trg2_JetLowPt.root";
  TString trigger3 = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_trg3_TrackHighPt.root";
  TString trigger4 = "../kstar-in-jets/Results_draw/JetUncorrpTspectra_trg4_TrackLowPt.root";
  
  TFile* data2 = TFile::Open(MB); 
  TFile* data2_trg1 = TFile::Open(trigger1);
  TFile* data2_trg2 = TFile::Open(trigger2);
  TFile* data2_trg3 = TFile::Open(trigger3);
  TFile* data2_trg4 = TFile::Open(trigger4);
  
  if(!data2)   return;
  if(!data2_trg1)   return;
  if(!data2_trg2)   return;
  if(!data2_trg3)   return;
  if(!data2_trg4)   return;

  TString uncorrName("Jets_spectra;1");
  TString effName("eff_MB;1");
 
  TH1D* heff = (TH1D*)data1->Get(effName);
  TH1D* huncorr = (TH1D*)data2->Get(uncorrName);
  TH1D* huncorr_trg1 = (TH1D*)data2_trg1->Get(uncorrName);
  TH1D* huncorr_trg2 = (TH1D*)data2_trg2->Get(uncorrName);
  TH1D* huncorr_trg3 = (TH1D*)data2_trg3->Get(uncorrName);
  TH1D* huncorr_trg4 = (TH1D*)data2_trg4->Get(uncorrName);
  
  TH1D* hcorr = (TH1D*)huncorr->Clone("hcorr");
  TH1D* hcorr_trg1 = (TH1D*)huncorr_trg1->Clone("hcorr_trg1");
  TH1D* hcorr_trg2 = (TH1D*)huncorr_trg2->Clone("hcorr_trg2");
  TH1D* hcorr_trg3 = (TH1D*)huncorr_trg3->Clone("hcorr_trg3");
  TH1D* hcorr_trg4 = (TH1D*)huncorr_trg4->Clone("hcorr_trg4");
 
  HistoAxisTitles(hcorr, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta} #frac{1}{eff x Acc} #frac{1}{BR}", 0.04, 62, 1.0, 0.04, 62, 1.7);
  hcorr->SetTitle("");
  
  
  hcorr->Divide(heff);
  hcorr_trg1->Divide(heff);
  hcorr_trg2->Divide(heff);
  hcorr_trg3->Divide(heff);
  hcorr_trg4->Divide(heff);

  TFile* fout = TFile::Open("/Users/jimun/Analysis/kstar-in-jets/Results_draw/corr_Jets.root", "RECREATE");
  fout->cd();
  hcorr->Write("Eff_corrected_pT_MB");
  hcorr_trg1->Write("Eff_corrected_pT_JetChHighpT");
  hcorr_trg2->Write("Eff_corrected_pT_JetChLowpT");
  hcorr_trg3->Write("Eff_corrected_pT_TrackHighpT");
  hcorr_trg4->Write("Eff_corrected_pT_TrackLowpT");
  

  //==================
  //
  //    signal Loss
  //
  //==================
  
  //==================
  // Gen data
  //==================
  TString mainfile3 = "../kstar-in-jets/Results_rootfile/Eff_MB.root";
  TFile* data3 = TFile::Open(mainfile3);
  if(!data3)   return;
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
  
  double yMin = 1e-6;
  double yMax = 1e-1;
  hcorr->SetMinimum(yMin);
  hcorr->SetMaximum(yMax);

  TCanvas* c = new TCanvas("", "", 1020, 800);
  c->SetLeftMargin(0.16);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.14);


  hcorr->SetLineWidth(2);
  hcorr->SetLineColor(blue3);
  hcorr->SetMarkerStyle(43); //circle
  hcorr->SetMarkerSize(3);
  hcorr->SetMarkerColor(blue3);
  hcorr->Draw();
  
  hcorr_trg1->SetLineWidth(2);
  hcorr_trg1->SetLineColor(red);
  hcorr_trg1->SetMarkerStyle(89);
  hcorr_trg1->SetMarkerSize(2);
  hcorr_trg1->SetMarkerColor(red);
  hcorr_trg1->Draw("SAME");
  
  hcorr_trg2->SetLineWidth(2);
  hcorr_trg2->SetLineColor(orange);
  hcorr_trg2->SetMarkerStyle(89);
  hcorr_trg2->SetMarkerSize(2);
  hcorr_trg2->SetMarkerColor(orange);
  hcorr_trg2->Draw("SAME");
  
  hcorr_trg3->SetLineWidth(2);
  hcorr_trg3->SetLineColor(green);
  hcorr_trg3->SetMarkerStyle(89);
  hcorr_trg3->SetMarkerSize(2);
  hcorr_trg3->SetMarkerColor(green);
  hcorr_trg3->Draw("SAME");

  hcorr_trg4->SetLineWidth(2);
  hcorr_trg4->SetLineColor(purple);
  hcorr_trg4->SetMarkerStyle(89);
  hcorr_trg4->SetMarkerSize(2);
  hcorr_trg4->SetMarkerColor(purple);
  hcorr_trg4->Draw("SAME");

  c->SetLogy();

  TLatex* txtInfoAAA = DrawTLatex(0.555, 0.868, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.555, 0.810, "Jets, #sqrt{s} = 13.6 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.555, 0.748, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.555, 0.695, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.166, 0.694, 0.5, 0.904);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hcorr, "MB", "pl"); //pl : point line
  legenddd->AddEntry(hcorr_trg1, "JetHighpT skim", "pl"); //pl : point line
  legenddd->AddEntry(hcorr_trg2, "JetLowpT skim", "pl"); //pl : point line
  legenddd->AddEntry(hcorr_trg3, "TrackHighpT skim", "pl"); //pl : point line
  legenddd->AddEntry(hcorr_trg4, "TrackLowpT skim", "pl"); //pl : point line
  
  legenddd->Draw("SAME");
  
  c->SaveAs("../kstar-in-jets/Plot/kstarInjets/"+localpath+"corr_Jets.png", "RECREATE");
}


