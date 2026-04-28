#include <TH1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>

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

void signalLoss()
{
  cout.precision(6);
  gStyle->SetOptStat(0);

  TString mainfile3 = "../kstar-in-jets/Results_rootfile/signalLoss_test.root";
  TFile* data3 = TFile::Open(mainfile3);
  if(!data3)   return;

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
  
  //==================
  // signal Loss
  //==================
  const Int_t nPtBins = 11;
  Double_t xBins[nPtBins+1] = {0.0, 0.8, 1.0, 1.2, 1.4, 1.8, 2.3, 2.8, 3.4, 4.0, 5.0, 8.0};
  int deepblue = TColor::GetColor("#167ec7");
  int yellorange = TColor::GetColor("#f05902");

  TString RecName("kstar-in-o-o/hEffGen_pT;1");
  TString GenName("kstar-in-o-o/hGen_pT_Raw;1");

  TH1D* hRecR = new TH1D("hRecRebin", "hRecRebin", nPtBins, xBins);
  TH1D* hGenR = new TH1D("hGenRebin", "hGenRebin", nPtBins, xBins);

  TH1D* hRec = (TH1D*)data3->Get(RecName);
  TH1D* hGen = (TH1D*)data3->Get(GenName);

  hRecR = (TH1D*)hRec->Rebin(nPtBins, "hRecRebin", xBins);
  hGenR = (TH1D*)hGen->Rebin(nPtBins, "hGenRebin", xBins);

  hRecR->Sumw2();
  hGenR->Sumw2();
  hRecR->Divide(hGenR);
 
  TFile* fout = TFile::Open("../kstar-in-jets/Results_draw/SignalLoss_Inclusive_MB.root", "RECREATE");
  hRecR->Write("signalLoss_MB");

  fout->Write();
  fout->Close();

  HistoAxisTitles(hRecR, "#it{p}_{T} [GeV/#it{c}]", "Signal Loss", 0.04, 62, 1.0, 0.04, 62, 1.1);
  hRecR->SetLineWidth(2);
  hRecR->SetLineColor(deepblue);
  hRecR->SetMarkerStyle(43); //circle
  hRecR->SetMarkerSize(3);
  hRecR->SetMarkerColor(deepblue);

  hRecR->SetMinimum(0.69);
  hRecR->SetMaximum(0.78);
  
  hRecR->SetTitle("");
  hRecR->Draw();

  TLine* EL = new TLine(0, eventLoss, 8.0, eventLoss);
  EL->SetLineStyle(2);
  EL->SetLineWidth(3);
  EL->SetLineColor(yellorange);
  EL->Draw("SAME");

  TLegend* legend = new TLegend(0.67, 0.789, 0.938, 0.92);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(EL, "Event loss", "l");
  legend->Draw("SAME");
  
}
