#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace TMath;

void signalLoss()
{
  TString mainfile = "../kstar-in-jets/Results_rootfile/MB_EffwTpcTof.root";
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;

  TString mainfile2 = "../kstar-in-jets/Results_draw/MBuncorrectedpTspectra.root";
  TFile* data2 = TFile::Open(mainfile2);
  if(!data2)   return;

  const Int_t nPtBins = 11;
  Double_t xBins[nPtBins+1] = {0.0, 0.8, 1.0, 1.2, 1.4, 1.8, 2.3, 2.8, 3.4, 4.0, 5.0, 8.0};

  TString kstarGenName("kstar-in-o-o/hEffGen_pT;1");
  TString kstarRawName("BC_MBspectra;1");

  // The number of generated kstar 
  TH1D* hGen = (TH1D*)data->Get(kstarGenName);
  TH1D* hGenR = new TH1D("hGenRebin", "hGenRebin", nPtBins, xBins);
  hGenR = (TH1D*)hGen->Rebin(nPtBins, "hGenRebin", xBins);
  double nGenR = hGenR->Integral();
  cout << "nGenR: " << nGenR << endl;

  // The number of raw kstar
  TH1D* hRaw = (TH1D*)data2->Get(kstarRawName);
  double nRaw = hRaw->Integral();
  cout << "nRaw: " << nRaw << endl;

  // Signal Loss Ratio
  double Loss_ratio = nGenR/nRaw;
  cout << "Loss Ratio: " << Loss_ratio <<endl;

  
  
}
