#include <TH1.h>
#include <TStyle.h>
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

void eff()
{
  cout.precision(6);
  gStyle->SetOptStat(0);
  
  TString mainfile = "../kstar-in-jets/Results_rootfile/Eff_MB_server.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  const Int_t nPtBins = 9;
  Double_t xBins[nPtBins+1] = {0.0, 1.2, 1.8, 2.4, 3.2, 4.0, 5.0, 6.0, 7.0, 8.0};
  
  int deepblue = TColor::GetColor("#167ec7");
  
  TString effRecName("kstar-in-o-o/hEffRec_pT;1");
  TString effGenName("kstar-in-o-o/hEffGen_pT;1");
  TString RecEventName("kstar-in-o-o/nEvents_MC;1");
  TString GenEventName("kstar-in-o-o/nEvents_MC_True;1");
  
  /* TString RecEventName("kstar-in-o-o/nEvents;1"); */
  /* TString GenEventName("kstar-in-o-o/nEvents_Gen;1"); */

  TH1D* hRecR = new TH1D("hRecRebin", "hRecRebin", nPtBins, xBins);
  TH1D* hGenR = new TH1D("hGenRebin", "hGenRebin", nPtBins, xBins);
  
  TH1D* hRec = (TH1D*)data->Get(effRecName);
  TH1D* hGen = (TH1D*)data->Get(effGenName);

  TH1D* nEvRec = (TH1D*)data->Get(RecEventName);
  TH1D* nEvGen = (TH1D*)data->Get(GenEventName);

  hRec->Scale(1.0/(nEvRec->GetBinContent(3))); // the number of events of MB
  hGen->Scale(1.0/(nEvGen->GetBinContent(1))); // the number of events of Jets

  hRecR = (TH1D*)hRec->Rebin(nPtBins,"hRecRebin",xBins);
  hGenR = (TH1D*)hGen->Rebin(nPtBins,"hGenRebin",xBins);

  hRec->Divide(hGen);
  hRecR->Divide(hGenR);
 
  hRec->Scale(1.0/0.66); //BR
  hRecR->Scale(1.0/0.66); //BR

  TFile* fout = TFile::Open("../kstar-in-jets/Results_draw/eff_Jet_MB.root", "RECREATE");
  hRecR->Write("eff_MB");

  fout->Write();
  fout->Close();

  hRec->GetXaxis()->SetRangeUser(0, 8);
  HistoAxisTitles(hRecR, "#it{p}_{T} [GeV/#it{c}]", "Efficiency x Acceptance", 0.04, 62, 1.0, 0.04, 62, 1.1);
  hRecR->SetLineWidth(2);
  hRecR->SetLineColor(deepblue);
  hRecR->SetMarkerStyle(43); //circle
  hRecR->SetMarkerSize(3);
  hRecR->SetMarkerColor(deepblue);
  
  hRecR->SetTitle("");
  hRecR->Draw();
  
}

