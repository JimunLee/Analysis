#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TLine.h>
#include <iostream>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
/*
HOW TO USE:

1. .L Drawer_Newer.C++

2. Drawer_Newer("Phi_pp16l_FineBinning_01_06_18.root", "MCPHI.root","Phi_Yield_MB_27bin_04_06_18.root", 0, 2)

That should do it! Fitting starts at Line 448 (Ch.3)
*/

//Optimal bins
      const Int_t nPtBins = 12;
      Double_t xBins[nPtBins+1] = 
        {0.9, 
         1.2, 1.4, 1.6,  
         1.8, 2.0, 2.2, 
         2.6, 3.0, 3.4,
         4.0, 5.0, 8.0};

      const Int_t nOldPtBins = 14;
      Double_t xOldBins[nOldPtBins+1] = 
        {0.5, 0.7, 0.9, 
         1.2, 1.4, 1.6,  
         1.8, 2.0, 2.2, 
         2.6, 3.0, 3.4,
         4.0, 5.0, 8.0};


Double_t imlow = 0.95;
Double_t imtop = 1.20;
void ShapeMyHist(TH1D **hist);
Double_t rap_correction(Double_t* x);
Double_t mybw(Double_t* x, Double_t* par);
Double_t Background(Double_t* x, Double_t* par);
Double_t allFit(Double_t* x, Double_t* par);
void Drawer_Voig(TString mainfile, TString outputfile,  const Double_t Fitlow, const Double_t Fithigh)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(01111);
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);                
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.08); 

               
  //TABLE OF CONTENTS
  // Ch1: Loading & Booking
  // Ch2: Reweighing
  // Ch3: Fitting
  // Ch4: Yield Extraction
  mainfile = "DataInputFiles/"+mainfile;

/////////////////////////////////////////////////////////////////////////////
//Ch1: Loading & Booking
////////////////////////////////////////////////////////////////////////////
//////////////////////////Loading Files//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  TFile* data = TFile::Open(mainfile);   
  if(!data)
    return;
/////////////////////////////////////////////////////////////////////////////
//////////////////////////Create Binning/////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//////////////////Selecting Data + Fetching Histograms///////////////////////
///////////////////////////////////////////////////////////////////////////// 
  char* charSelect1 = new char [30];
  char* charSelect2 = new char [30];
  TString SignalName("phi-in-jets/hUSS_INSIDE;1");
  TString BackgroundName("phi-in-jets/hLSS_INSIDE;1");

  
  TH3D* hLSS = (TH3D*)data->Get(BackgroundName);
  TH3D* hUSS = (TH3D*)data->Get(SignalName);

    TH1D* hCalib1 = new TH1D("hCalib1", "hCalib1", nPtBins, xBins);
  TH1D* hCalib2 = new TH1D("hCalib2", "hCalib2", nPtBins, xBins);
  TH1D* hCalib3 = new TH1D("hCalib3", "hCalib3", nPtBins, xBins);
  TH1D* hCalib5 = new TH1D("hCalib5", "hCalib5", nPtBins, xBins);
  TH1D* hCalib6 = new TH1D("hCalib6", "hCalib6", nPtBins, xBins);
  TH1D* h_binYield = new TH1D("_binYield", "binPhi Yield", nPtBins, xBins);

      
/////////////////////////////////////////////////////////////////////////////
//////////////////////Booking New Histograms, constants and characters///////
/////////////////////////////////////////////////////////////////////////////
  TH1D* hLSS_Y [nPtBins];
  TH1D* hUSS_Y [nPtBins];
  TH1D* hLSS_Y_Copy [nPtBins];
  TH1D* hUSS_Y_Copy [nPtBins];
  TH1D* hUSS_Y_Bkg [nPtBins];
  TH1D* hFitVData [nPtBins];


  TF1* fBreit [nPtBins];
  TF1* fBreitCopy[nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];

  for (Int_t i=0;i<nPtBins;i++) 
    {
      hLSS_Y[i]=new TH1D(Form("hLSS_Y_%i",i), Form("hLSS_Y_%i",i), 400, 0.95, 1.35);
      hUSS_Y[i]=new TH1D(Form("hUSS_Y_%i",i), Form("hUSS_Y_%i",i), 400, 0.95, 1.35);
      hUSS_Y_Bkg[i]=new TH1D(Form("hUSS_Y_Bkg_%i",i), Form("hUSS_Y_Bkg_%i",i), 400, 0.95, 1.35);


      hLSS_Y_Copy[i]=new TH1D(Form("hLSS_Y_Copy_%i",i), Form("hLSS_Y_Copy_%i",i), 400, 0.95, 1.35);
      hUSS_Y_Copy[i]=new TH1D(Form("hUSS_Y_Copy_%i",i), Form("hUSS_Y_Copy_%i",i), 400, 0.95, 1.35);

      hFitVData[i] = new TH1D(Form("hFitVData_%i",i), "(Fit-Data)/Peak Value", 400, 0.95,1.35);
      hFitVData[i]->Sumw2();


    }



   
//Ch2: Reweighing
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  //gROOT->SetBatch(kTRUE);
  Double_t ptlow = 0;
  Double_t pthigh = 0;
  //Toggle to 1 for single run
  Int_t CanvasCounter = 0;
  Int_t PadCounter = 0;


  
  for(Int_t j=0;j<nPtBins;j++) //THIS SELECTS THE OVERARCHING PT-BINNING
    {
      ptlow = xBins[j];
      pthigh = xBins[j+1];

      cout<<"ptlow"<<endl;

      (hUSS->GetYaxis())->SetRangeUser(ptlow, pthigh);
      (hLSS->GetYaxis())->SetRangeUser(ptlow, pthigh);

       hLSS_Y[j]->Add(hLSS->Project3D("z o"),1);
       hUSS_Y[j]->Add(hUSS->Project3D("z o"),1);

       Double_t Nlow = 0.97;
       Double_t Nhigh =1.0;

       Double_t Scale1 = hUSS_Y[j]->Integral(hUSS_Y[j]->FindBin(Nlow), hUSS_Y[j]->FindBin(Nhigh)) / hLSS_Y[j]->Integral(hLSS_Y[j]->FindBin(Nlow), hLSS_Y[j]->FindBin(Nhigh));
       hLSS_Y[j]->Scale(Scale1);     
       // hUSS_Y[j]->Rebin(2);
       // hUSS_Y[j]->Scale(0.5);
       // hLSS_Y[j]->Rebin(2);
       // hLSS_Y[j]->Scale(0.5);
       // hUSS_Y_Copy[j]->Rebin(2);
       // hUSS_Y_Copy[j]->Scale(0.5);
       // hLSS_Y_Copy[j]->Rebin(2);
       // hLSS_Y_Copy[j]->Scale(0.5);


       hUSS_Y_Copy[j]->Add(hUSS_Y[j]);
       hLSS_Y_Copy[j]->Add(hLSS_Y[j]);       

       hUSS_Y[j]->Add(hLSS_Y[j], -1);

       TCanvas* c1 = new TCanvas("c1","c1", 1920, 1080);
       hUSS_Y_Copy[j]->Draw();
       hLSS_Y_Copy[j]->SetLineColor(2);
       hLSS_Y_Copy[j]->Draw("SAME");

       c1->SaveAs(Form("BackgroundPlots/Minv_Plot_%f_to_%f.png",ptlow,pthigh),"RECREATE");
        
    }
  Double_t Bglow =0.95;
  Double_t Bghigh = 1.35;
  Double_t Siglow = Fitlow;
  Double_t Sighigh = Fithigh;


  for(int f=0; f<nPtBins;f++){


    ptlow = xBins[f];
    pthigh = xBins[f+1];

    TCanvas* c2 = new TCanvas(Form("c_%d", f),Form("c_%d", f),1920,1080);

    Int_t binz = 0;
    while(binz<25)
      {
	if(hUSS_Y[f]->GetBinCenter(binz)<0.995)
	  hUSS_Y[f]->SetBinError(binz,0);
	else
	  break;
	binz++;
      }

    // hUSS_Y_Bkg[f]->Rebin(2);    
    // hUSS_Y_Bkg[f]->Scale(0.5);

    hUSS_Y_Bkg[f]->Add(hUSS_Y[f]);
    Int_t binx = hUSS_Y_Bkg[f]->FindBin(1.007);
    while(binx<hUSS_Y_Bkg[f]->FindBin(1.033))
      {
	hUSS_Y_Bkg[f]->SetBinError(binx,0);
	binx++;
      }


    
  
    
    fBreit[f] = new TF1(Form("BwFit_%i",f),"[0]*TMath::Voigt(x-[2], [1], [3])",Siglow, Sighigh);
    fBreit[f]->SetLineColor(2);
    fBreit[f]->SetParameter(0,1);
    fBreit[f]->SetParameter(1,1);
    fBreit[f]->SetParameter(2,1.019455);
    fBreit[f]->SetParameter(3,0.004266);
    fBreitCopy[f] = new TF1(Form("BwFit_Copy_%i",f),"[0]*TMath::Voigt(x-[2], [1], [3])",0, 1000000);
    fBreitCopy[f]->SetLineColor(2);
    fFit[f] = new TF1(Form("FullFit_%i",f), "[0]*TMath::Voigt(x-[2], [1], [3]) + pol2(4)", Siglow, Sighigh);

    
    fBg[f] = new TF1(Form("BkgFit_%i",f), "pol2", Bglow, Bghigh);
    fBg[f]->SetLineColor(1);
    fBg[f]->SetLineStyle(2);
    fBg[f]->SetLineWidth(4);
    fBg[f]->SetRange(Siglow, Sighigh);
    // Fit background in order to estimate its initial values for the complete fit    
    hUSS_Y_Bkg[f]->Fit(fBg[f], "0", "", Bglow, Bghigh);
    hUSS_Y_Bkg[f]->Fit(fBg[f], "0", "", Bglow, Bghigh);
    hUSS_Y_Bkg[f]->Fit(fBg[f], "0", "", Bglow, Bghigh);

    // Extract parameters from previous BG and MC Signal fits, and plug them into combined fit.
    fFit[f]->SetLineColor(2);
    fFit[f]->SetParLimits(1,0.0001,0.003);
    fFit[f]->SetParameter(3,0.004266);
    fFit[f]->SetParameter(4, fBg[f]->GetParameter(0));
    fFit[f]->SetParameter(5, fBg[f]->GetParameter(1));
    fFit[f]->SetParameter(6, fBg[f]->GetParameter(2));
    fFit[f]->SetParameter(0, fBreit[f]->GetParameter(0));
    fFit[f]->SetParameter(1, fBreit[f]->GetParameter(1));
    fFit[f]->SetParameter(2, fBreit[f]->GetParameter(2));




    Bool_t fCalib=kTRUE;
    if(fCalib){
      TFile* CalibFile = TFile::Open("CalibFile.root");
      TH1D* hCalibPar1 = (TH1D*)CalibFile->Get("hCalib1");
      TH1D* hCalibPar2 = (TH1D*)CalibFile->Get("hCalib2");
      TH1D* hCalibPar3 = (TH1D*)CalibFile->Get("hCalib3");
      TH1D* hCalibPar5 = (TH1D*)CalibFile->Get("hCalib5");
      TH1D* hCalibPar6 = (TH1D*)CalibFile->Get("hCalib6");

      fFit[f]->FixParameter(1,hCalibPar1->GetBinContent(f+1));
      fFit[f]->FixParameter(2,hCalibPar2->GetBinContent(f+1));
      fFit[f]->FixParameter(3,hCalibPar3->GetBinContent(f+1));
      // fFit[f]->FixParameter(5,hCalibPar5->GetBinContent(f+1));
      // fFit[f]->FixParameter(6,hCalibPar6->GetBinContent(f+1));

      // fFit[f]->SetParameter(1,hCalibPar1->GetBinContent(f+1));
      // fFit[f]->SetParameter(2,hCalibPar2->GetBinContent(f+1));
      // fFit[f]->SetParameter(3,hCalibPar3->GetBinContent(f+1));
    }


    hUSS_Y[f]->Fit(fFit[f], "SR0", "", Siglow, Sighigh);
    hUSS_Y[f]->Fit(fFit[f], "SR0", "", Siglow, Sighigh);


    TFitResultPtr ptr;
    ptr = hUSS_Y[f]->Fit(fFit[f], "SR", "", Siglow, Sighigh);

    if(!fCalib){
      hCalib1->SetBinContent(f+1, fFit[f]->GetParameter(1));
      hCalib1->SetBinError(f+1, fFit[f]->GetParError(1));
      hCalib2->SetBinContent(f+1, fFit[f]->GetParameter(2));
      hCalib2->SetBinError(f+1, fFit[f]->GetParError(2));
      hCalib3->SetBinContent(f+1, fFit[f]->GetParameter(3));
      hCalib3->SetBinError(f+1, fFit[f]->GetParError(3));
      hCalib5->SetBinContent(f+1, fFit[f]->GetParameter(5));
      hCalib5->SetBinError(f+1, fFit[f]->GetParError(5));
      hCalib6->SetBinContent(f+1, fFit[f]->GetParameter(6));
      hCalib6->SetBinError(f+1, fFit[f]->GetParError(6));
 
    }    


    fBg[f]->FixParameter(0, fFit[f]->GetParameter(4));
    fBg[f]->FixParameter(1, fFit[f]->GetParameter(5));
    fBg[f]->FixParameter(2, fFit[f]->GetParameter(6));
    fBg[f]->SetRange(Siglow, Sighigh);

    fBreit[f]->SetParameters(fFit[f]->GetParameters());
    fBreitCopy[f]->SetParameters(fFit[f]->GetParameters());

    ShapeMyHist(&hUSS_Y[f]);
    hUSS_Y[f]->SetTitle("");
    hUSS_Y[f]->GetXaxis()->SetTitle("M_{inv} [Gev/c^{2}]");
    hUSS_Y[f]->GetYaxis()->SetTitle("Entries");
    hUSS_Y[f]->Draw();
    fFit[f]->Draw("SAME");
    fBg[f]->Draw("SAME");





    hUSS_Y[f]->GetXaxis()->SetRangeUser(0.99,1.07);
    Double_t I_min = 1.01;
    Double_t I_max = 1.03;
     
    TLine* Upp = new TLine(I_max, hUSS_Y[f]->GetMinimum()-50, I_max,hUSS_Y[f]->GetMaximum()+50);
    TLine* Ner = new TLine(I_min, hUSS_Y[f]->GetMinimum()-50, I_min,hUSS_Y[f]->GetMaximum()+50);
    TLine* Noll = new TLine(0.95, 0, 1.20, 0);
    Upp->Draw("SAME");
    Ner->Draw("SAME");
    Noll->Draw("SAME");

    c2->SaveAs(Form("FitPlots/Fit_Plot_%f_to_%f.png",ptlow,pthigh),"RECREATE");

    Double_t Peak = 0;
    for(Int_t bin = 1; bin<=400; bin++)
      {
	Double_t Xaxis = hUSS_Y[f]->GetXaxis()->GetBinCenter(bin);
	if(Xaxis>=Bglow)
	  {
	    Double_t Data = hUSS_Y[f]->GetBinContent(bin);
	    Double_t Fit = fFit[f]->Eval(Xaxis); //extracting Xaxis's value
	    Double_t Ratio = Fit-Data;
	    Double_t erData = hUSS_Y[f]->GetBinError(bin);
	    if(Data>Peak)
	      Peak = Data;
	    hFitVData[f]->SetBinContent(bin, Ratio);
	  }
	if(Xaxis>Bghigh)
	  break;
      }


    Double_t binw = hUSS_Y[f]->GetXaxis()->GetBinWidth(1);
    Double_t Sint = fBreit[f]->Integral(I_min, I_max);
    Double_t yield = Sint/binw;
    Double_t Sint_top = fFit[f]->GetMaximum();
    Double_t Tails = 0;

    Double_t binyield = 0;
    Double_t binpiece1 = 0;
    Double_t binpiece2= 0;
    cout<<hUSS_Y[f]->GetBinCenter(hUSS_Y[f]->FindBin(I_max))<<endl;
    for(Int_t binz = hUSS_Y[f]->FindBin(I_min); binz < hUSS_Y[f]->FindBin(I_max); binz++)
      {
	binpiece2 +=  fBg[f]->Eval(hUSS_Y[f]->GetBinCenter(binz));
      }
    if(binpiece2<0)
      binpiece2=0;

    Double_t YieldError = 0;
    binyield = hUSS_Y[f]->IntegralAndError(hUSS_Y[f]->FindBin(I_min), hUSS_Y[f]->FindBin(I_max)-1, YieldError) - binpiece2;

    Double_t FPars[4] = {fFit[f]->GetParameter(0),fFit[f]->GetParameter(1),fFit[f]->GetParameter(2),fFit[f]->GetParameter(3)};
    Double_t BPar[3] = {fFit[f]->GetParameter(4),fFit[f]->GetParameter(5), fFit[f]->GetParameter(6)};
    TMatrixDSym covTot,covGG,covGG_err;
    covTot = ptr->GetCovarianceMatrix();
    covGG = covTot.GetSub(0,3,0,3);
    covGG_err = covTot.GetSub(4,6,4,6);
    Double_t error_bin1 = (1/binw)*(fBreit[f]->IntegralError(Siglow,I_min,fBreit[f]->GetParameters(),covGG.GetMatrixArray()));
    Double_t error_bin2 = (1/binw)*(fBreit[f]->IntegralError(I_max,Sighigh,fBreit[f]->GetParameters(),covGG.GetMatrixArray()));
    Double_t Bg_err = (1/binw)*(fBg[f]->IntegralError(I_min,I_max,fBg[f]->GetParameters(),covGG_err.GetMatrixArray()));

    // if(binyield<0)
    //   return;
    Tails = TMath::Abs(fBreitCopy[f]->Integral(0.987354,I_min)/binw);
    Tails+= TMath::Abs(fBreitCopy[f]->Integral(I_max,TMath::Infinity())/binw);

    Double_t Error = TMath::Sqrt(YieldError*YieldError +  error_bin1* error_bin1 + error_bin2* error_bin2 +Bg_err*Bg_err);
    binyield+=Tails;
    yield+=Tails;
    cout<<"Peak Value: "<<Sint_top<<endl;
    cout<<"Yield: "<<yield<<endl;
    cout<<"Bins: "<<xBins[f]<< " - " <<xBins[f+1]<<endl;
    cout<<"f: "<<f<<endl;
    cout<<"Ndf: "<<fFit[f]->GetNDF()<<endl;
    cout<<"Bin Yield: "<<binyield<<endl;


    h_binYield->SetBinContent(f+1, binyield/h_binYield->GetXaxis()->GetBinWidth(f+1));  
    h_binYield->SetBinError(f+1,Error/h_binYield->GetXaxis()->GetBinWidth(f+1));


  }


  
  //    TFile* file = new TFile("CalibFile.root", "RECREATE");
   // hCalib1->Write();
   // hCalib2->Write();
   // hCalib3->Write();
   // hCalib5->Write();
   // hCalib6->Write();

  TFile* file = new TFile("ResultsFile.root", "RECREATE");
 ShapeMyHist(&h_binYield);
 h_binYield->GetYaxis()->SetTitle("Counts");
 h_binYield->SetTitle("#frac{d^{2}N_{#phi}}{dp_{T}d#eta}");
 h_binYield->GetXaxis()->SetTitle("p_{T} [GeV/c]");



  
  TH1D* hVtx = (TH1D*)data->Get("phi-in-jets/nEvents"); 
  Double_t nEvents = hVtx->GetBinContent(1);
  TH1D* hCorr_binYield = new TH1D("hCorr_binYield", "hCorr_binYield",nPtBins,xBins);
  // hCorr_binYield->Scale(1.0/1.6);
  hCorr_binYield->Add(h_binYield);
  hCorr_binYield->Scale(1.0/nEvents);

  for(int i=1; i<15; i++){
    hCorr_binYield->SetBinContent(i, hCorr_binYield->GetBinContent(i)/ hCorr_binYield->GetBinWidth(i));
    hCorr_binYield->SetBinError(i, hCorr_binYield->GetBinError(i)/ hCorr_binYield->GetBinWidth(i));
  }

  hCorr_binYield->Write();

  // Double_t nVtxOk = hVtx->GetBinContent(3)/(hVtx->GetBinContent(2)+hVtx->GetBinContent(3));
  // Double_t nNoVtx = nVtxOk*hVtx->GetBinContent(1);
  // cout<<"nNoVtx : "<<nNoVtx<<endl;
  // cout<<"nEvents : "<<nEvents<<endl;
  //   nEvents += nNoVtx;
 // hCorr_binYield->GetYaxis()->SetTitle("Counts");
 // hCorr_binYield->SetTitle("#frac{d^{2}N_{#phi}}{dp_{T}dy} #frac{1}{e_{corr}} #frac{1}{N_{vtx<10} + (N_{noVtx} #scale[0.5]{#bullet} e_{vtxcut})}");
 // hCorr_binYield->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////  
  
}//End of void loop
void ShapeMyHist(TH1D **hist)
{
  TH1D *htemp = *hist;
  htemp->SetLineWidth(4);
  htemp->SetLabelSize(0.035,"xyz");
  htemp->SetLabelOffset(0.0,"y");
  htemp->SetLabelOffset(0.01,"x");
  htemp->SetLabelColor(kBlack,"xyz");
  htemp->SetTitleSize(0.039,"xyz");
  htemp->SetTitleOffset(1.4,"y");
  htemp->SetTitleOffset(1.2,"x");
  htemp->SetLabelFont(62,"xy");
  htemp->SetTitleFont(62,"xy");
}
Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 2/(TMath::Pi()); // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*(arg2/(arg3 + arg4));
}


Double_t Background(Double_t* x, Double_t*par)
{

  // return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]) ;
  return (par[0]+par[1]*x[0]) ;
}

Double_t allFit(Double_t* x, Double_t* par)
{
  return mybw(x,par) + Background(x, &par[3]);
}
Double_t rap_correction(Double_t* x)
{
  Double_t pt = x[0];  
  Double_t eta  = 0.8;
  Double_t mass = 1.019;
  const Double_t mt = TMath::Sqrt(pt*pt + mass*mass);
  const Double_t rap = TMath::ASinH(pt/mt*TMath::SinH(eta));
  //  return rap/eta;
  return rap/eta;
}


