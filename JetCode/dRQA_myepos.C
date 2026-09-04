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
  //TString mainfile_EPOS = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/dR_EPOS.root";
  TString mainfile_EPOS = "/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/myepos.root";
  
  TFile* data_PYTHIA = TFile::Open(mainfile_PYTHIA);
  if(!data_PYTHIA)   return;
  TFile* data_HERWIG = TFile::Open(mainfile_HERWIG);
  if(!data_HERWIG)   return;
  TFile* data_EPOS = TFile::Open(mainfile_EPOS);
  if(!data_EPOS)   return;
  

  TString origin_8_10("kstar-in-o-o/Kstar_pT_INJet_8_10;1");
  TString origin_10_12("kstar-in-o-o/Kstar_pT_INJet_10_12;1");
  TString recovered_8_12("kstar-in-o-o/Kstar_pT_6_8to8_IN_recoveredJet;1");

  //------------------------------------------------
  // PYTHIA
  
  // denominator
  
  TH1D* h_kstar_8_10_PY = (TH1D*)data_PYTHIA->Get(origin_8_10);
  TH1D* h_kstar_10_12_PY = (TH1D*)data_PYTHIA->Get(origin_10_12);
  TH1D* h_kstar_8_12_PY = (TH1D*)h_kstar_8_10_PY->Clone("8_10_origin");
  h_kstar_8_12_PY->Add(h_kstar_10_12_PY);

  TH1D* h_recovered_kstar_8_12_PY = (TH1D*)data_PYTHIA->Get(recovered_8_12);

  TH1D* denominator_PY = (TH1D*)h_kstar_8_12_PY->Clone("8_12_Truth");
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
  TH1D* numerator_PY = (TH1D*)h_kstar_8_12_PY->Clone("8_12_Truth");

  cout << "Scaled A+B: " << numerator_PY->Integral() << endl;
  
  //------------------------------------------------
  // HERWIG

  // denominator
  
  TH1D* h_kstar_8_10_HW = (TH1D*)data_HERWIG->Get(origin_8_10);
  TH1D* h_kstar_10_12_HW = (TH1D*)data_HERWIG->Get(origin_10_12);
  TH1D* h_kstar_8_12_HW = (TH1D*)h_kstar_8_10_HW->Clone("8_10_origin");
  h_kstar_8_12_HW->Add(h_kstar_10_12_HW);

  TH1D* h_recovered_kstar_8_12_HW = (TH1D*)data_HERWIG->Get(recovered_8_12);

  TH1D* denominator_HW = (TH1D*)h_kstar_8_12_HW->Clone("8_12_Truth");
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
  TH1D* numerator_HW = (TH1D*)h_kstar_8_12_HW->Clone("8_12_Truth");

  cout << "Scaled A+B: " << numerator_HW->Integral() << endl;
  
  //------------------------------------------------
  // EPOS

  // denominator
  
  TH1D* h_kstar_8_10_EPOS = (TH1D*)data_EPOS->Get(origin_8_10);
  TH1D* h_kstar_10_12_EPOS = (TH1D*)data_EPOS->Get(origin_10_12);
  TH1D* h_kstar_8_12_EPOS = (TH1D*)h_kstar_8_10_EPOS->Clone("8_10_origin");
  h_kstar_8_12_EPOS->Add(h_kstar_10_12_EPOS);

  TH1D* h_recovered_kstar_8_12_EPOS = (TH1D*)data_EPOS->Get(recovered_8_12);

  TH1D* denominator_EPOS = (TH1D*)h_kstar_8_12_EPOS->Clone("8_12_Truth");
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
  TH1D* numerator_EPOS = (TH1D*)h_kstar_8_12_EPOS->Clone("8_12_Truth");

  cout << "Scaled A+B: " << numerator_EPOS->Integral() << endl;
  cout << "-----------------------------------------" << endl;

  double X_average = (X_PY + X_HW + X_EPOS)/3.0;
  cout << "X_average: " << X_average << endl;

  //================================================
  // Common binning
  //================================================
  numerator_PY->Rebin(4);
  denominator_PY->Rebin(4);

  //myepos
  numerator_EPOS->Rebin(4);
  denominator_EPOS->Rebin(4);

  // Make ratios
  numerator_PY->Divide(denominator_PY);
  numerator_HW->Divide(denominator_HW);
  numerator_EPOS->Divide(denominator_EPOS);

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
  numerator_HW->SetLineColor(kRed + 1);
  numerator_HW->SetMarkerColor(kRed + 1);
  numerator_HW->SetMarkerStyle(21);
  numerator_HW->SetMarkerSize(1.0);
  numerator_HW->SetLineWidth(2);

  // EPOS
  numerator_EPOS->SetLineColor(kBlue + 1);
  numerator_EPOS->SetMarkerColor(kBlue + 1);
  numerator_EPOS->SetMarkerStyle(33);
  numerator_EPOS->SetMarkerSize(1.3);
  numerator_EPOS->SetLineWidth(2);

  HistoAxisTitles(numerator_PY,"#it{p}_{T}^{K^{*0}} (GeV/#it{c})","Correction ratio",
		  0.045, 42, 1.10,
		  0.045, 42, 1.25);
  
  numerator_PY->SetTitle("");
  numerator_PY->GetXaxis()->SetRangeUser(0.0, 10.0);
  
  numerator_PY->Draw();
  numerator_HW->Draw("SAME");
  numerator_EPOS->Draw("SAME");

  TLegend* leg = new TLegend(0.26, 0.75, 0.52, 0.91);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(numerator_PY,   "PYTHIA", "lep");
  leg->AddEntry(numerator_HW,   "HERWIG", "lep");
  leg->AddEntry(numerator_EPOS, "EPOS",   "lep");

  leg->Draw("SAME");
  models->SaveAs("/Users/jimun/Analysis/kstar-in-jets/Plot/kstarInjets/MC/Correction/model_dependency.png", "RECREATE");
}
