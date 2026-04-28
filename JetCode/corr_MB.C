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
  TFile* data1 = TFile::Open(mainfile1);
  if(!data1)   return;

  //==================
  // uncorrected data
  //==================
  TString MB = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_MB.root";
  TString trigger1 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg1_JetHighPt.root";
  TString trigger2 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg2_JetLowPt.root";
  TString trigger3 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg3_TrackHighPt.root";
  TString trigger4 = "../kstar-in-jets/Results_draw/UncorrectedpTspectra_trg4_TrackLowPt.root";
  
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
  //==================
  // Gen data
  //==================
  TString mainfile3 = "../kstar-in-jets/Results_rootfile/signalLoss_test.root";
  TFile* data3 = TFile::Open(mainfile3);
  if(!data3)   return;

  TString uncorrName("Inclusive_spectra;1");
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

  HistoAxisTitles(hcorr, "#it{p}_{T} [GeV/#it{c}]", "#frac{1}{N_{evt}} #frac{d^{2}N_{K^{*0}}}{d#it{p}_{T}d#eta} #frac{1}{eff x Acc x BR} #frac{#epsilon_{EL}}{#epsilon_{SL}}", 0.04, 62, 1.0, 0.04, 62, 1.7);
  TCanvas* c = new TCanvas("", "", 1020, 800);
  c->SetLeftMargin(0.16);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.14);
  
  hcorr->Divide(heff);
  hcorr_trg1->Divide(heff);
  hcorr_trg2->Divide(heff);
  hcorr_trg3->Divide(heff);
  hcorr_trg4->Divide(heff);

  //hcorr->Draw();
  /* hcorr_trg1->Draw("SAME"); */
  /* hcorr_trg2->Draw("SAME"); */
  /* hcorr_trg3->Draw("SAME"); */
  /* hcorr_trg4->Draw("SAME"); */
  
  c->SetLogy();

  TFile* fout = TFile::Open("/Users/jimun/Analysis/kstar-in-jets/Results_draw/corr_Inclusive.root", "RECREATE");
  fout->cd();
  hcorr->Write("Eff_corrected_pT_MB");
  hcorr_trg1->Write("Eff_corrected_pT_JetChHighpT");
  hcorr_trg2->Write("Eff_corrected_pT_JetChLowpT");
  hcorr_trg3->Write("Eff_corrected_pT_TrackHighpT");
  hcorr_trg4->Write("Eff_corrected_pT_TrackLowpT");

  //==================
  //| Event Loss
  // ==================
  TString GenEventName("kstar-in-o-o/nEvents_Gen;1");

  TH1D* hEvGen = (TH1D*)data3->Get(GenEventName);
  
  double RecColl = hEvGen->GetBinContent(4);
  double GenColl = hEvGen->GetBinContent(1);
  
  double eventLoss = RecColl/GenColl;
  
  cout<<"number of Recon.Coll: " << RecColl <<endl;
  cout<<"number of Gen.Coll: " << GenColl <<endl;
  cout<<"Event Loss: " << eventLoss <<endl;

  hcorr->Scale(eventLoss);
  hcorr->Write("event_corrected_pT");

  //==================
  // signal Loss
  //==================
  TString mainfile4 = "../kstar-in-jets/Results_draw/SignalLoss_Inclusive_MB.root";
  TFile* data4 = TFile::Open(mainfile4);
  if(!data4)   return;

  TString SLName("signalLoss_MB;1");
  TH1D* hsigLoss = (TH1D*)data4->Get(SLName);

  hcorr->Divide(hsigLoss);
  
  fout->Close();

  hcorr->SetLineWidth(2);
  hcorr->SetLineColor(blue3);
  hcorr->SetMarkerStyle(43); //circle
  hcorr->SetMarkerSize(3);
  hcorr->SetMarkerColor(blue3);
  hcorr->Draw("SAME");
  
  /* hcorr_trg1->SetLineWidth(2); */
  /* hcorr_trg1->SetLineColor(red); */
  /* hcorr_trg1->SetMarkerStyle(89); */
  /* hcorr_trg1->SetMarkerSize(2); */
  /* hcorr_trg1->SetMarkerColor(red); */
  /* hcorr_trg1->Draw("SAME"); */
  
  /* hcorr_trg2->SetLineWidth(2); */
  /* hcorr_trg2->SetLineColor(orange); */
  /* hcorr_trg2->SetMarkerStyle(89); */
  /* hcorr_trg2->SetMarkerSize(2); */
  /* hcorr_trg2->SetMarkerColor(orange); */
  /* hcorr_trg2->Draw("SAME"); */
  
  /* hcorr_trg3->SetLineWidth(2); */
  /* hcorr_trg3->SetLineColor(green); */
  /* hcorr_trg3->SetMarkerStyle(89); */
  /* hcorr_trg3->SetMarkerSize(2); */
  /* hcorr_trg3->SetMarkerColor(green); */
  /* hcorr_trg3->Draw("SAME"); */

  /* hcorr_trg4->SetLineWidth(2); */
  /* hcorr_trg4->SetLineColor(purple); */
  /* hcorr_trg4->SetMarkerStyle(89); */
  /* hcorr_trg4->SetMarkerSize(2); */
  /* hcorr_trg4->SetMarkerColor(purple); */
  /* hcorr_trg4->Draw("SAME"); */

  TLatex* txtInfoAAA = DrawTLatex(0.555, 0.868, "ALICE Work In Progress", 0.05, 62, 1);
  TLatex* txtInfoBBB = DrawTLatex(0.555, 0.810, "Inclusive, #sqrt{s} = 13.6 TeV", 0.05, 42, 1);
  TLatex* txtInfoCCC = DrawTLatex(0.555, 0.748, "K^{*0}(892) #rightarrow K^{+}#pi^{-}", 0.05, 42, 1);
  TLatex* txtInfoDDD = DrawTLatex(0.555, 0.695, "|#eta| < 0.8", 0.05, 42, 1);
  TLegend* legenddd = new TLegend(0.166, 0.152, 0.5, 0.363);
  legenddd->SetBorderSize(0);
  legenddd->SetFillStyle(0);
  legenddd->SetTextSize(0.045);
  legenddd->SetTextFont(42);
  legenddd->AddEntry(hcorr, "Inclusive", "pl"); //pl : point line
  /* legenddd->AddEntry(hcorr_trg1, "JetHighpT skim", "pl"); //pl : point line */
  /* legenddd->AddEntry(hcorr_trg2, "JetLowpT skim", "pl"); //pl : point line */
  /* legenddd->AddEntry(hcorr_trg3, "TrackHighpT skim", "pl"); //pl : point line */
  /* legenddd->AddEntry(hcorr_trg4, "TrackLowpT skim", "pl"); //pl : point line */
  
  legenddd->Draw("SAME");
  
  c->SaveAs("../kstar-in-jets/Plot/kstarInjets/"+localpath+"correctedpT.png", "RECREATE");
}


