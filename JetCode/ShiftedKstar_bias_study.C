#include <TH1.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>

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

void ShiftedKstar()
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString main_PY = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_PYTHIA.root";
  TFile* data_PY = TFile::Open(main_PY);
  if(!data_PY)
    return;

  TString main_EPOS = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_EPOS.root";
  TFile* data_EPOS = TFile::Open(main_EPOS);
  if(!data_EPOS)
    return;

  TString main_HW = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_HERWIG.root";
  TFile* data_HW = TFile::Open(main_HW);
  if(!data_HW)
    return;


  TString recovered_8_12("kstar-in-o-o/Kstar_pT_6_8to8_IN_recoveredJet;1");
  TString origin_6_8("kstar-in-o-o/Kstar_pT_INJet_6_8;1");
  TString origin_8_10("kstar-in-o-o/Kstar_pT_INJet_8_10;1");
  TString origin_10_12("kstar-in-o-o/Kstar_pT_INJet_10_12;1");

  //=================================
  //|
  //|        PYTHIA
  //|
  TH1D* h_recovered_8_12_PY = (TH1D*)data_PY->Get(recovered_8_12);
  TH1D* h_origin_6_8_PY = (TH1D*)data_PY->Get(origin_6_8);
  TH1D* h_origin_8_10_PY = (TH1D*)data_PY->Get(origin_8_10);
  TH1D* h_origin_10_12_PY = (TH1D*)data_PY->Get(origin_10_12);

  TH1D* h_origin_8_12_PY = (TH1D*)h_origin_8_10_PY->Clone("h_clone_8_10");
  h_origin_8_12_PY->Add(h_origin_10_12_PY);
  
  HistoAxisTitles(h_origin_6_8_PY, "K^{*0} #it{p}_{T} [GeV/#it{c}]", "Normalized distribution", 0.05, 62, 1.0, 0.05, 62, 1.0);
  h_origin_6_8_PY->SetTitle(Form(""));
  h_origin_6_8_PY->SetLineColor(1);
  h_recovered_8_12_PY->SetLineColor(2);
  h_origin_8_12_PY->SetLineColor(4);

  TCanvas* c1 = new TCanvas("", "", 880, 600);
  TPad* p1 = new TPad("", "", 0, 0.3, 1, 1);
  p1->SetTopMargin(0.05);
  p1->SetLeftMargin(0.15);
  p1->Draw();
  p1->cd();
  p1->SetBottomMargin(0.01); 

  h_origin_6_8_PY->Scale(1.0/h_origin_6_8_PY->Integral());
  h_recovered_8_12_PY->Scale(1.0/h_recovered_8_12_PY->Integral());
  h_origin_8_12_PY->Scale(1.0/h_origin_8_12_PY->Integral());

  h_origin_6_8_PY->Rebin(4);
  h_recovered_8_12_PY->Rebin(4);
  h_origin_8_12_PY->Rebin(4);

  h_origin_6_8_PY->GetXaxis()->SetRangeUser(0., 14.);
  h_origin_6_8_PY->GetYaxis()->SetRangeUser(0., 0.12);
  
  h_origin_6_8_PY->Draw();
  h_recovered_8_12_PY->Draw("SAME");
  h_origin_8_12_PY->Draw("SAME");


  
  TLegend* L1 = new TLegend(0.46, 0.58, 0.74, 0.81);
  L1->SetBorderSize(0);
  L1->SetTextSize(0.05);
  L1->AddEntry(h_origin_6_8_PY, "K^{*0} in 6-8 GeV/c jets", "l");
  L1->AddEntry(h_recovered_8_12_PY, "K^{*0} in migrated jets (8-12 #rightarrow 6-8 GeV/c)", "l");
  L1->AddEntry(h_origin_8_12_PY, "K^{*0} in 8-12 GeV/c jets", "l");
  L1->Draw("SAME");

  TLatex* txtInfoA = DrawTLatex(0.471, 0.858, "PYTHIA", 0.06, 62, 1);
  
  
  c1->cd();
  TPad* p2 = new TPad("", "", 0, 0, 1, 0.3);
  p2->SetLeftMargin(0.15);
  p2->SetBottomMargin(0.3);
  p2->SetTopMargin(0.0001);
  p2->Draw();
  p2->cd();


  TH1D* h_Ratio_recOver8_12_PY = (TH1D*)h_recovered_8_12_PY->Clone("hRatio_rec_Over_8_12");
  TH1D* h_Ratio_recOver6_8_PY = (TH1D*)h_recovered_8_12_PY->Clone("hRatio_rec_Over_6_8");
  HistoAxisTitles(h_Ratio_recOver8_12_PY, "K^{*0} #it{p}_{T} [GeV/#it{c}]", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  
  h_Ratio_recOver8_12_PY->GetXaxis()->SetRangeUser(0., 14.);
  h_Ratio_recOver6_8_PY->GetXaxis()->SetRangeUser(0., 14.);
  

  h_Ratio_recOver8_12_PY->SetTitle("");
  h_Ratio_recOver8_12_PY->GetXaxis()->SetLabelSize(0.09);
  h_Ratio_recOver8_12_PY->GetYaxis()->SetLabelSize(0.08);

  h_Ratio_recOver8_12_PY->SetMinimum(-0.1);
  h_Ratio_recOver8_12_PY->SetMaximum(2.5);
  h_Ratio_recOver8_12_PY->GetYaxis()->SetNdivisions(408);
  
  
  h_Ratio_recOver8_12_PY->SetLineColor(4);
  h_Ratio_recOver6_8_PY->SetLineColor(1);
  
  h_Ratio_recOver8_12_PY->Divide(h_origin_8_12_PY);
  h_Ratio_recOver8_12_PY->Draw();

  h_Ratio_recOver6_8_PY->Divide(h_origin_6_8_PY);
  h_Ratio_recOver6_8_PY->Draw("SAME");

  TLegend* L2 = new TLegend(0.65, 0.68, 0.79, 0.92);
  L2->SetBorderSize(0);
  L2->SetTextSize(0.1);
  L2->AddEntry(h_Ratio_recOver8_12_PY, "Migrated / 8-12 GeV/c jets", "l");
  L2->AddEntry(h_Ratio_recOver6_8_PY, "Migrated / 6-8 GeV/c jets", "l");
  L2->Draw();

  //=================================
  //|
  //|        HERWIG
  //|
  TH1D* h_recovered_8_12_HW = (TH1D*)data_HW->Get(recovered_8_12);
  TH1D* h_origin_6_8_HW = (TH1D*)data_HW->Get(origin_6_8);
  TH1D* h_origin_8_10_HW = (TH1D*)data_HW->Get(origin_8_10);
  TH1D* h_origin_10_12_HW = (TH1D*)data_HW->Get(origin_10_12);

  TH1D* h_origin_8_12_HW = (TH1D*)h_origin_8_10_HW->Clone("h_clone_8_10");
  h_origin_8_12_HW->Add(h_origin_10_12_HW);
  
  
  HistoAxisTitles(h_origin_6_8_HW, "K*0 pT [GeV/#it{c}]", "Counts", 0.05, 62, 1.0, 0.05, 62, 1.0);
  h_origin_6_8_HW->SetTitle(Form(""));
  h_origin_6_8_HW->GetXaxis()->SetRangeUser(0., 14.);
  h_origin_8_12_HW->GetXaxis()->SetRangeUser(0., 14.);
  h_recovered_8_12_HW->GetXaxis()->SetRangeUser(0., 14.);

  h_origin_6_8_HW->SetLineColor(1);
  h_recovered_8_12_HW->SetLineColor(2);
  h_origin_8_12_HW->SetLineColor(4);

  TCanvas* c1_HW = new TCanvas("", "", 880, 600);
  TPad* p1_HW = new TPad("", "", 0, 0.3, 1, 1);
  p1_HW->SetTopMargin(0.05);
  p1_HW->SetLeftMargin(0.15);
  p1_HW->Draw();
  p1_HW->cd();
  p1_HW->SetBottomMargin(0.01);

  h_origin_6_8_HW->Scale(1.0/h_origin_6_8_HW->Integral());
  h_recovered_8_12_HW->Scale(1.0/h_recovered_8_12_HW->Integral());
  h_origin_8_12_HW->Scale(1.0/h_origin_8_12_HW->Integral());

  h_origin_6_8_HW->GetYaxis()->SetRangeUser(0., 0.12);

  h_origin_6_8_HW->Draw();
  h_recovered_8_12_HW->Draw("SAME");
  h_origin_8_12_HW->Draw("SAME");

  TLegend* L1_HW = new TLegend(0.51, 0.58, 0.79, 0.81);
  L1_HW->SetBorderSize(0);
  L1_HW->SetTextSize(0.05);
  L1_HW->AddEntry(h_origin_6_8_HW, "K^{*0} in 6-8 GeV/c jets", "l");
  L1_HW->AddEntry(h_recovered_8_12_HW, "K^{*0} in migrated jets (8-12 #rightarrow 6-8 GeV/c)", "l");
  L1_HW->AddEntry(h_origin_8_12_HW, "K^{*0} in 8-12 GeV/c jets", "l");
  L1_HW->Draw();

  TLatex* txtInfoB = DrawTLatex(0.525, 0.858, "HERWIG", 0.06, 62, 1);

  c1_HW->cd();
  TPad* p2_HW = new TPad("", "", 0, 0, 1, 0.3);
  p2_HW->SetLeftMargin(0.15);
  p2_HW->SetBottomMargin(0.3);
  p2_HW->SetTopMargin(0.0001);
  p2_HW->Draw();
  p2_HW->cd();


  TH1D* h_Ratio_recOver8_12_HW = (TH1D*)h_recovered_8_12_HW->Clone("hRatio_rec_Over_8_12");
  TH1D* h_Ratio_recOver6_8_HW = (TH1D*)h_recovered_8_12_HW->Clone("hRatio_rec_Over_6_8");
  HistoAxisTitles(h_Ratio_recOver8_12_HW, "K*0 pT [GeV/#it{c}]", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  
  h_Ratio_recOver8_12_HW->GetXaxis()->SetRangeUser(0., 14.);
  h_Ratio_recOver6_8_HW->GetXaxis()->SetRangeUser(0., 14.);
  
  h_Ratio_recOver8_12_HW->SetTitle("");
  h_Ratio_recOver8_12_HW->GetXaxis()->SetLabelSize(0.09);
  h_Ratio_recOver8_12_HW->GetYaxis()->SetLabelSize(0.08);

  h_Ratio_recOver8_12_HW->SetMinimum(-0.1);
  h_Ratio_recOver8_12_HW->SetMaximum(2.5);
  h_Ratio_recOver8_12_HW->GetYaxis()->SetNdivisions(408);
  
  
  h_Ratio_recOver8_12_HW->SetLineColor(4);
  h_Ratio_recOver6_8_HW->SetLineColor(1);
  
  h_Ratio_recOver8_12_HW->Divide(h_origin_8_12_HW);
  h_Ratio_recOver8_12_HW->Draw();

  h_Ratio_recOver6_8_HW->Divide(h_origin_6_8_HW);
  h_Ratio_recOver6_8_HW->Draw("SAME");

  TLegend* L2_HW = new TLegend(0.65, 0.68, 0.79, 0.92);
  L2_HW->SetBorderSize(0);
  L2_HW->SetTextSize(0.1);
  L2_HW->AddEntry(h_Ratio_recOver8_12_HW, "Migrated / 8-12 GeV/c jets", "l");
  L2_HW->AddEntry(h_Ratio_recOver6_8_HW, "Migrated / 6-8 GeV/c jets", "l");
  L2_HW->Draw();


    
  //=================================
  //|
  //|        EPOS
  //|
  TH1D* h_recovered_8_12_EPOS = (TH1D*)data_EPOS->Get(recovered_8_12);
  TH1D* h_origin_6_8_EPOS = (TH1D*)data_EPOS->Get(origin_6_8);
  TH1D* h_origin_8_10_EPOS = (TH1D*)data_EPOS->Get(origin_8_10);
  TH1D* h_origin_10_12_EPOS = (TH1D*)data_EPOS->Get(origin_10_12);

  TH1D* h_origin_8_12_EPOS = (TH1D*)h_origin_8_10_EPOS->Clone("h_clone_8_10");
  h_origin_8_12_EPOS->Add(h_origin_10_12_EPOS);
  
  
  HistoAxisTitles(h_origin_6_8_EPOS, "K*0 pT [GeV/#it{c}]", "Counts", 0.05, 62, 1.0, 0.05, 62, 1.0);
  h_origin_6_8_EPOS->SetTitle(Form(""));
  h_origin_6_8_EPOS->GetXaxis()->SetRangeUser(0., 14.);
  h_origin_8_12_EPOS->GetXaxis()->SetRangeUser(0., 14.);
  h_recovered_8_12_EPOS->GetXaxis()->SetRangeUser(0., 14.);

  h_origin_6_8_EPOS->SetLineColor(1);
  h_recovered_8_12_EPOS->SetLineColor(2);
  h_origin_8_12_EPOS->SetLineColor(4);

  TCanvas* c1_EPOS = new TCanvas("", "", 880, 600);
  TPad* p1_EPOS = new TPad("", "", 0, 0.3, 1, 1);
  p1_EPOS->SetTopMargin(0.05);
  p1_EPOS->SetLeftMargin(0.15);
  p1_EPOS->Draw();
  p1_EPOS->cd();
  p1_EPOS->SetBottomMargin(0.01);

  h_origin_6_8_EPOS->Scale(1.0/h_origin_6_8_EPOS->Integral());
  h_recovered_8_12_EPOS->Scale(1.0/h_recovered_8_12_EPOS->Integral());
  h_origin_8_12_EPOS->Scale(1.0/h_origin_8_12_EPOS->Integral());

  h_origin_6_8_EPOS->GetYaxis()->SetRangeUser(0., 0.12);

  h_origin_6_8_EPOS->Draw("SAME");
  h_recovered_8_12_EPOS->Draw("SAME");
  h_origin_8_12_EPOS->Draw("SAME");

  TLegend* L1_EPOS = new TLegend(0.51, 0.58, 0.79, 0.81);
  L1_EPOS->SetBorderSize(0);
  L1_EPOS->SetTextSize(0.05);

 L1_HW->AddEntry(h_origin_6_8_HW, "K^{*0} in 6-8 GeV/c jets", "l");
  L1_HW->AddEntry(h_recovered_8_12_HW, "K^{*0} in migrated jets (8-12 #rightarrow 6-8 GeV/c)", "l");
  L1_HW->AddEntry(h_origin_8_12_HW, "K^{*0} in 8-12 GeV/c jets", "l");

  
  L1_EPOS->AddEntry(h_origin_6_8_EPOS,  "K^{*0} in 6-8 GeV/c jets", "l");
  L1_EPOS->AddEntry(h_recovered_8_12_EPOS, "K^{*0} in migrated jets (8-12 #rightarrow 6-8 GeV/c)", "l");
  L1_EPOS->AddEntry(h_origin_8_12_EPOS, "K^{*0} in 8-12 GeV/c jets", "l");
  L1_EPOS->Draw("SAME");

  TLatex* txtInfoC = DrawTLatex(0.525, 0.858, "EPOS", 0.06, 62, 1);

  c1_EPOS->cd();
  TPad* p2_EPOS = new TPad("", "", 0, 0, 1, 0.3);
  p2_EPOS->SetLeftMargin(0.15);
  p2_EPOS->SetBottomMargin(0.3);
  p2_EPOS->SetTopMargin(0.0001);
  p2_EPOS->Draw();
  p2_EPOS->cd();


  TH1D* h_Ratio_recOver8_12_EPOS = (TH1D*)h_recovered_8_12_EPOS->Clone("hRatio_rec_Over_8_12");
  TH1D* h_Ratio_recOver6_8_EPOS = (TH1D*)h_recovered_8_12_EPOS->Clone("hRatio_rec_Over_6_8");
  HistoAxisTitles(h_Ratio_recOver8_12_EPOS, "K*0 pT [GeV/#it{c}]", "", 0.12, 62, 0.8, 0.1, 62, 1.0);
  
  h_Ratio_recOver8_12_EPOS->GetXaxis()->SetRangeUser(0., 14.);
  h_Ratio_recOver6_8_EPOS->GetXaxis()->SetRangeUser(0., 14.);
  
  h_Ratio_recOver8_12_EPOS->SetTitle("");
  h_Ratio_recOver8_12_EPOS->GetXaxis()->SetLabelSize(0.09);
  h_Ratio_recOver8_12_EPOS->GetYaxis()->SetLabelSize(0.08);

  h_Ratio_recOver8_12_EPOS->SetMinimum(-0.1);
  h_Ratio_recOver8_12_EPOS->SetMaximum(2.5);
  h_Ratio_recOver8_12_EPOS->GetYaxis()->SetNdivisions(408);
  
  
  h_Ratio_recOver8_12_EPOS->SetLineColor(4);
  h_Ratio_recOver6_8_EPOS->SetLineColor(1);
  
  h_Ratio_recOver8_12_EPOS->Divide(h_origin_8_12_EPOS);
  h_Ratio_recOver8_12_EPOS->Draw();

  h_Ratio_recOver6_8_EPOS->Divide(h_origin_6_8_EPOS);
  h_Ratio_recOver6_8_EPOS->Draw("SAME");

  TLegend* L2_EPOS = new TLegend(0.65, 0.68, 0.79, 0.92);
  L2_EPOS->SetBorderSize(0);
  L2_EPOS->SetTextSize(0.1);
  L2_EPOS->AddEntry(h_Ratio_recOver8_12_EPOS, "Migrated / 8-12 GeV/c jets", "l");
  L2_EPOS->AddEntry(h_Ratio_recOver6_8_EPOS, "Migrated / 6-8 GeV/c jets", "l");
  L2_EPOS->Draw();

  
}
