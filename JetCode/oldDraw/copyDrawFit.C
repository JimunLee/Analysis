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

//=========================================
// Guide of function
// => DrawVoigt("mainfile_name.root", "Inclusive", 1)
//=========================================
Double_t mybw(Double_t* x, Double_t* par);
void Shreck();
void DrawVoigt(TString Kind, TString DataNum)
{
  std::cout.precision(6); //keep the rest of the code unchanged
  gStyle->SetOptStat(0);

  const Int_t nPtBins = 9;
  Double_t xBins[nPtBins+1];

  if (Kind == "Outside"){
    //Outside range
    Double_t outside_xBins[nPtBins+1] =
      {0.8, 1.2, 1.5, 1.8, 2.4, 3.0, 3.6, 4.5, 5.5, 8.0};
    copy(begin(outside_xBins), end(outside_xBins), xBins);
  } else if (Kind == "Inclusive"){
    //Inclusive range
    Double_t inclusive_xBins[nPtBins+1] =
      {0.8, 1.2, 1.5, 1.8, 2.4, 3.0, 3.6, 4.5, 5.5, 8.0};
    copy(begin(inclusive_xBins), end(inclusive_xBins), xBins);
  } else if (Kind == "Inside"){
    //Inside range
    Double_t inside_xBins[nPtBins+1] = 
      {0.8, 1.2, 1.5, 1.8, 2.4, 3.0, 3.6, 4.5, 5.5, 8.0};
    copy(begin(inside_xBins), end(inside_xBins), xBins);
  }else{
    cout << "Please specify either 'Inside', 'Outside', or 'Inclusive'." << endl;
    return;
  }
  
  Double_t imlow = 0.75;
  Double_t imtop = 1.15;
  Double_t nBins = 400.0;
  
  TString localpath = "Tot_combined/"+Kind; 

  //=========================================
  // Data file inputs
  //=========================================

  TString mainfile = "../phi-in-jets/Result_rootfile/pidtest_combined.root";
  
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;
  //=========================================
  //| Information : Data names              | 
  //| combined data ;1  single data ;2      |
  //| hUSS_INSIDE;1 hUSS_OUTSIDE;1  hUSS;1  |
  //=========================================
  TString SignalName("phi-in-jets/hUSS;"+DataNum);
  TString BackgroundName("phi-in-jets/hLSS;"+DataNum);
  if(DataNum != 1 && DataNum != 2) return;
  
  TH3D* hLSS = (TH3D*)data->Get(BackgroundName);
  TH3D* hUSS = (TH3D*)data->Get(SignalName);
  
  //=========================================
  // 1-D projection plots
  //=========================================
  TH1D* hLSS_Y [nPtBins];
  TH1D* hUSS_Y [nPtBins];
  
  TH1* hUSS_Y_Bkg [nPtBins];
  TH1* hUSS_Y_Peak [nPtBins];
  TH1D* h_binYield = new TH1D("_binYield", "binPhi Yield", nPtBins, xBins);
  
  for (Int_t jpt=0; jpt<nPtBins; jpt++)
    {
      hLSS_Y[jpt]=new TH1D(Form("hLSS_Y_%i", jpt), Form("hLSS_Y_%i", jpt), nBins, imlow, imtop);
      hUSS_Y[jpt]=new TH1D(Form("hUSS_Y_%i", jpt), Form("hUSS_Y_%i", jpt), nBins, imlow, imtop);
      
      hUSS_Y_Bkg[jpt]=new TH1D(Form("hUSS_Y_Bkg_%i", jpt), Form("hUSS_Y_Bkg_%i", jpt), nBins, imlow, imtop);
      hUSS_Y_Peak[jpt]=new TH1D(Form("hUSS_Y_Peak_%i", jpt), Form("hUSS_Y_Peak_%i", jpt), nBins, imlow, imtop);
    }

  //=========================================
  // Loop over pT bins
  //=========================================
  TH1D* hLSS_Y_COPY [nPtBins];
  TH1D* hUSS_Y_COPY [nPtBins];
  
  Double_t ReBinNum = 4;
  Double_t ptlow = 0;
  Double_t pthigh = 0;
  
  for(Int_t ipt=0; ipt<nPtBins; ipt++)
    {
      ptlow = xBins[ipt];
      pthigh = xBins[ipt+1];

      (hLSS->GetYaxis())->SetRangeUser(ptlow, pthigh);
      (hUSS->GetYaxis())->SetRangeUser(ptlow, pthigh);

      hLSS_Y[ipt]->Add(hLSS->Project3D("z o"),1);
      hUSS_Y[ipt]->Add(hUSS->Project3D("z o"),1);
      
      Double_t Nlow = 1.05;
      Double_t Nhigh = 1.1;
      Double_t Scale = hUSS_Y[ipt]->Integral(hUSS_Y[ipt]->FindBin(Nlow), hUSS_Y[ipt]->FindBin(Nhigh))
	/ hLSS_Y[ipt]->Integral(hLSS_Y[ipt]->FindBin(Nlow), hLSS_Y[ipt]->FindBin(Nhigh));
      hLSS_Y[ipt]->Scale(Scale);

      hUSS_Y[ipt]->Rebin(ReBinNum);
      hUSS_Y[ipt]->Scale(1./ReBinNum);
      hLSS_Y[ipt]->Rebin(ReBinNum);
      hLSS_Y[ipt]->Scale(1./ReBinNum);
      
      hLSS_Y[ipt]->SetLineColor(2); // 2 : red

      //=========================================
      // Draw 1-D projection plots (hUSS, hLSS)
      //=========================================
      
      TCanvas* plots = new TCanvas(Form("plots_%d",ipt), Form("plots_%d",ipt), 550, 440);

      plots->SetLeftMargin(0.12);
      plots->SetRightMargin(0.1);
      plots->SetTopMargin(0.1);
      plots->SetBottomMargin(0.12);

      hUSS_Y_COPY[ipt] = (TH1D*)hUSS_Y[ipt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[ipt] = (TH1D*)hLSS_Y[ipt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hUSS_Y[ipt]->Add(hLSS_Y[ipt], -1);
      
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitle("M_{inv}[GeV/c^{2}]");
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitle("Counts");
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->SetTitle(Form(""));
      hUSS_Y_COPY[ipt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[ipt]->Draw();
      hUSS_Y_COPY[ipt]->SetLineColor(kRed);
      hLSS_Y_COPY[ipt]->Draw("SAME");
      hLSS_Y_COPY[ipt]->SetLineColor(kBlue);

      TLatex* txtInfo = new TLatex(0.64, 0.31, Form("%.1f < pT < %.1f GeV/c", ptlow, pthigh));
      txtInfo->SetNDC(kTRUE); // NDC 좌표 사용
      txtInfo->SetTextSize(0.04);
      txtInfo->SetTextFont(62);
      txtInfo->SetTextColor(1);
      txtInfo->Draw("SAME");
      TLegend* legend = new TLegend(0.63,0.17,0.83,0.29);
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.05);
      legend->AddEntry(hUSS_Y_COPY[ipt], "Signal", "lpf");
      legend->AddEntry(hLSS_Y_COPY[ipt], "Background", "lpf");
      legend->Draw("SAME");
  
      //   plots->SaveAs(Form("../phi-in-jets/Plot/"+localpath+"/projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");
      
      //=========================================
      // Draw 1-D projection hUSS-hLSS 
      //=========================================
      TCanvas* substraction = new TCanvas(Form("substraction_%d",ipt), "", 550, 440);

      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.1);
      substraction->SetBottomMargin(0.12);
      
      hUSS_Y_COPY[ipt] = (TH1D*)hUSS_Y[ipt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[ipt] = (TH1D*)hLSS_Y[ipt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitle("M_{inv}[GeV/c^{2}]");
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitle("Counts");
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->SetTitle(Form(""));
      hUSS_Y_COPY[ipt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[ipt]->Draw();
      hUSS_Y[ipt]->Draw("SAME");

      TLegend* legend_sub = new TLegend(0.14,0.16,0.46,0.20);
      legend_sub->SetFillColor(4000);
      legend_sub->SetBorderSize(0);
      legend_sub->SetTextSize(0.04);
      legend_sub->AddEntry(hUSS_Y_COPY[ipt], "hUSS - hLSS", "lpf");
      legend_sub->Draw("SAME");      
      TLatex* txtInfo_sub = new TLatex(0.1624, 0.2385, Form("%.1f < pT < %.1f GeV/c", ptlow, pthigh));
      txtInfo_sub->SetNDC(kTRUE); // NDC 좌표 사용
      txtInfo_sub->SetTextSize(0.04);
      txtInfo_sub->SetTextFont(62);
      txtInfo_sub->SetTextColor(1);
      txtInfo_sub->Draw("SAME");
      // substraction->SaveAs(Form("../phi-in-jets/Plot/"+localpath+"/reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh));
    }
  
  Double_t Siglow = 0.85;
  Double_t Sighigh = 0.95;
  Double_t Bglow = 0.75;
  Double_t Bghigh = 1.15;

  TF1* fBreit [nPtBins];
  TF1* fBreitCopy [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  
  TH1D* hYield_Bin = new TH1D("binYield","binYield", nPtBins,xBins);
  TH1D* hYield_Peak = new TH1D("fitYield","fitYield;pT;N", nPtBins,xBins);  

  TF1* fBreit_only [nPtBins];
  TF1* fBg_only [nPtBins];
  
  for(int fpt=0; fpt<nPtBins; fpt++)
    {
      hUSS_Y_Bkg[fpt]->Add(hUSS_Y[fpt]);
      hUSS_Y_Peak[fpt]->Add(hUSS_Y[fpt]);
      
      for(int bin = hUSS_Y_Bkg[fpt]->FindBin(Siglow); bin<=hUSS_Y_Bkg[fpt]->FindBin(Sighigh); bin++)
	hUSS_Y_Bkg[fpt]->SetBinError(bin,0);
      
      ptlow = xBins[fpt];
      pthigh = xBins[fpt+1];
      
      //=========================================
      // Fitting Peak Function 
      //=========================================
      //| Information : K*                      | 
      //| gamma(1,3) : 0.0508                   |
      //| mass(2) : 0.89166                     |
      //=========================================
      fBreit[fpt] = new TF1(Form("BwFit_%i",fpt),"[0]*TMath::Voigt(x-[2], [1], [3])", Siglow, Sighigh);
      fBreit[fpt]->SetLineColor(2); //Red:2

      fBreit[fpt]->SetParameter(0,1);
      fBreit[fpt]->SetParameter(1,0.0508);
      fBreit[fpt]->SetParameter(2,0.89166);
      fBreit[fpt]->SetParameter(3,0.0508);
      fBreitCopy[fpt] = new TF1(Form("BwFit_Copy_%i",fpt),"[0]*TMath::Voigt(x-[2], [1], [3])", Siglow, Sighigh);
      fBreitCopy[fpt]->SetLineColor(2);
      
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      Shreck();
      cout<<fBreit[fpt]->GetParameter(0)<<endl;
      cout<<fBreit[fpt]->GetParameter(1)<<endl;
      cout<<fBreit[fpt]->GetParameter(2)<<endl;
      Shreck();
      // add part : Draw only fBreit_only line
      fBreit_only[fpt] = new TF1(Form("BwFit_%i_only",fpt),"[0]*TMath::Voigt(x-[2], [1], [3])", Siglow, Sighigh);
      fBreit_only[fpt]->SetLineColor(5);
      fBreit_only[fpt]->FixParameter(0, fBreit[fpt]->GetParameter(4));
      fBreit_only[fpt]->FixParameter(1, fBreit[fpt]->GetParameter(5));
      fBreit_only[fpt]->FixParameter(2, fBreit[fpt]->GetParameter(6));

      TCanvas* test = new TCanvas();
      fBreit_only[fpt]->Draw();

      continue;
      
      //=========================================
      // Fitting Background Function 
      //=========================================
      if(fpt==4 || fpt==5 || fpt==6){
	fBg[fpt] = new TF1(Form("BkgFit_%i",fpt), "pol3", Bglow, Bghigh); //you have to be careful! This is 4 par
	fBg_only[fpt] = new TF1(Form("BgFitonly_%i",fpt), "pol3", Bglow, Bghigh);
      }
      else{
	fBg[fpt] = new TF1(Form("BkgFit_%i",fpt), "pol2", Bglow, Bghigh); //this has 3 par
	fBg_only[fpt] = new TF1(Form("BgFitonly_%i",fpt), "pol2", Bglow, Bghigh);
      }
      //fBg[fpt] = new TF1(Form("BkgFit_%i",fpt), "pol3", Bglow, Bghigh);
      fBg[fpt]->SetLineColor(1);//Black:1
      fBg[fpt]->SetLineStyle(2);
      fBg[fpt]->SetLineWidth(4);
      fBg[fpt]->SetRange(Bglow, Bghigh);
      
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);

       // add part : Draw only fBg_only
      fBg_only[fpt]->SetParameter(0, fBg[fpt]->GetParameter(4));
      fBg_only[fpt]->SetParameter(1, fBg[fpt]->GetParameter(5));
      fBg_only[fpt]->SetParameter(2, fBg[fpt]->GetParameter(6));
      if(fpt==4 || fpt==5 || fpt ==6){
	fBg_only[fpt]->SetParameter(3, fBg[fpt]->GetParameter(7));
      }
      //=========================================
      // Combine Peak + Background
      //=========================================
      fFit[fpt] = new TF1(Form("FullFit_%i",fpt), "[0]*TMath::Voigt(x-[2], [1], [3]) + pol3(4)", Bglow, Bghigh);
      fFit[fpt]->SetLineColor(2);
      
      //=========================================
      // Costrain Peak + Background  
      //=========================================
      fFit[fpt]->SetParLimits(1, 0.0001, 0.0003);
      
      // BACKGROUND
      fFit[fpt]->SetParameter(4, fBg[fpt]->GetParameter(0));
      fFit[fpt]->SetParameter(5, fBg[fpt]->GetParameter(1));
      fFit[fpt]->SetParameter(6, fBg[fpt]->GetParameter(2));
      if(fpt==4 || fpt==5|| fpt==6){
	fFit[fpt]->SetParameter(7, fBg[fpt]->GetParameter(3));
      }
      // PEAK
      fFit[fpt]->SetParameter(0, fBreit[fpt]->GetParameter(0));
      fFit[fpt]->SetParameter(1, fBreit[fpt]->GetParameter(1));
      fFit[fpt]->SetParameter(2, fBreit[fpt]->GetParameter(2));
      fFit[fpt]->SetParameter(3, fBreit[fpt]->GetParameter(3));
      
      
      // OK, WE NOW HAVE GOOD STARTING VALUES
      hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Siglow, Sighigh);
      hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Siglow, Sighigh);

      // Now we take the newely fitted background values in the full fit, and store them in the background-only fit.
      fBg[fpt]->SetParameter(0, fFit[fpt]->GetParameter(4));
      fBg[fpt]->SetParameter(1, fFit[fpt]->GetParameter(5));
      fBg[fpt]->SetParameter(2, fFit[fpt]->GetParameter(6));
      if(fpt==4 || fpt==5|| fpt ==6){
      fBg[fpt]->SetParameter(3, fFit[fpt]->GetParameter(7));
      }
      
      //==========================
      // Draw
      //==========================   

      TCanvas * lines = new TCanvas(Form("line_%d",fpt),"" , 800, 660);
      lines->SetLogy();
      lines->SetLeftMargin(0.12);
      lines->SetRightMargin(0.1);
      lines->SetTopMargin(0.1);
      lines->SetBottomMargin(0.12);

      hUSS_Y_COPY[fpt] = (TH1D*)hUSS_Y[fpt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[fpt] = (TH1D*)hLSS_Y[fpt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitle("M_{inv}[GeV/c^{2}]");
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitle("Entries");
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleSize(0.05);
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleOffset(1.0);
      hUSS_Y_COPY[fpt]->SetTitle(Form(""));
      hUSS_Y_COPY[fpt]->SetTitleSize(0.05, "t");
      hUSS_Y_COPY[fpt]->SetTitle(Form(""));
      hUSS_Y_COPY[fpt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[fpt]->SetLineColor(1);
      hUSS_Y_COPY[fpt]->Draw();
      
      fFit[fpt]->SetLineColor(2);
      fFit[fpt]->Draw("SAME");
      
      fBg[fpt]->SetLineColor(4);
      fBg[fpt]->Draw("SAME");

      fBreit_only[fpt]->SetLineColor(5);
      fBreit_only[fpt]->Draw("SAME");
  
      fBg_only[fpt]->SetLineColor(6);
      fBg_only[fpt]->Draw("SAME");

      TLatex* txtInfoB = new TLatex(0.18, 0.34, Form("%.1f < pT < %.1f GeV/c", ptlow, pthigh));
      txtInfoB->SetNDC(kTRUE); // use NDC
      txtInfoB->SetTextSize(0.04);
      txtInfoB->SetTextFont(62);
      txtInfoB->SetTextColor(1);
      txtInfoB->Draw("SAME");
      TLegend *legend = new TLegend(0.18,0.18,0.40,0.32);
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.05);
      legend->AddEntry(fBg[fpt], "Background Fit", "lpf");
      legend->AddEntry(fFit[fpt], "Signal Fit", "lpf");
      legend->Draw("SAME");
      
      // lines->SaveAs(Form("../phi-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
      
      //==========================
      // Reflecting an error
      //========================== 
      Int_t Siglow_bin = hUSS_Y[fpt]->FindBin(Siglow);
      Int_t Sighigh_bin = hUSS_Y[fpt]->FindBin(Sighigh);
      
      double binW = hUSS_Y[fpt]->GetBinWidth(1); //any number
      double bin_int = hUSS_Y[fpt]->Integral(Siglow_bin, Sighigh_bin);
      double peak_int = fFit[fpt]->Integral(Siglow,Sighigh)/binW;
      double bks_int = fBg[fpt]->Integral(Siglow,Sighigh)/binW;

      double dif_one = peak_int - bks_int;
      double dif_two = bin_int - bks_int;

      cout<<"BIN INT "<< bin_int<<endl;
      cout<<"PEAK INT : "<<peak_int<<endl;
      cout<<"PEAK YIELD : "<< dif_one<<endl;
      cout<<"BIN YIELD : "<< dif_two<<endl;

      hYield_Peak->SetBinContent(fpt+1, dif_one);
      hYield_Bin->SetBinContent(fpt+1, dif_two);
      
      TFitResultPtr ptr;
      ptr = hUSS_Y[fpt]->Fit(fFit[fpt], "SR", "", Siglow, Sighigh);
      fBreitCopy[fpt]->SetParameters(fFit[fpt]->GetParameters());
      
      TMatrixDSym covTot(fFit[fpt]->GetNpar() );
      TMatrixDSym covGG(fBreit[fpt]->GetNpar() );
      TMatrixDSym covGG_err(fBg[fpt]->GetNpar() );
      
      covTot = ptr->GetCovarianceMatrix();
      covGG = covTot.GetSub(0,3,0,3);
      covGG_err = covTot.GetSub(4,7,4,7);

      Double_t error_fit = (1/binW)*(fBreit[fpt]->IntegralError(Siglow,Sighigh,fBreit[fpt]->GetParameters(),covGG.GetMatrixArray()));
      Double_t Bg_err = (1/binW)*(fBg[fpt]->IntegralError(Siglow,Sighigh,fBg[fpt]->GetParameters(),covGG_err.GetMatrixArray()));
   
      cout<<"error_fit : "<<error_fit<<endl;
      cout<<"Bg_err : "<<Bg_err<<endl;
      
      Double_t Tails = 0;
      Tails = TMath::Abs(fBreitCopy[fpt]->Integral(0.749735,Siglow)/binW);
      Tails += TMath::Abs(fBreitCopy[fpt]->Integral(Sighigh,TMath::Infinity())/binW);

      Int_t binz = 0;
      while(binz<25)
	{

      cout<<"Bin Yield : "<<binyield<<endl;
      
      h_binYield->SetBinContent(fpt+1, binyield/h_binYield->GetXaxis()->GetBinWidth(fpt+1));
      h_binYield->SetBinError(fpt+1,Error/h_binYield->GetXaxis()->GetBinWidth(fpt+1));
    }

  //================================
  // Draw # of bin of (peak - bkg)
  //================================
  
  TCanvas* numbin = new TCanvas();
  hYield_Peak->Draw();
  hYield_Bin->SetLineColor(3);
  hYield_Peak->GetXaxis()->SetTitle("pT[GeV/c]");
  hYield_Bin->Draw("SAME");
  numbin->SaveAs("../phi-in-jets/Plot/"+localpath+"/hYield.png","RECREATE");

  TString Eventdata("phi-in-jets/nEvents;"+DataNum);
  TH1D* hVtx = (TH1D*)data->Get(Eventdata);
  Double_t nEvents = hVtx->GetBinContent(1);
  TH1D* hCorr_binYield = new TH1D("hCorr_binYield", "hCorr_binYield", nPtBins,xBins);
  
  hCorr_binYield->Add(h_binYield);
  hCorr_binYield->Scale(1.0/nEvents);
    
  // hCorr_binYield->GetYaxis()->SetTitle("#frac{d^{2}N_{K*}}{dpT dy}");
  hCorr_binYield->GetXaxis()->SetTitle("pT[GeV/c]");

  for(int i=1; i<10; i++)
    {
      hCorr_binYield->SetBinContent(i, hCorr_binYield->GetBinContent(i)/ hCorr_binYield->GetBinWidth(i));
      hCorr_binYield->SetBinError(i, hCorr_binYield->GetBinError(i)/ hCorr_binYield->GetBinWidth(i));
    }
  
  TCanvas* logbin = new TCanvas();
  hCorr_binYield->Draw();
  logbin->SetLogy();
  //logbin->SaveAs("../phi-in-jets/Plot/"+localpath+"/logbinhYield.png");
  
}


Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 2/(TMath::Pi()); // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*(arg2/(arg3 + arg4));
}

 void Shreck()
{
  std::cout << "⢀⡴⠑⡄⠀⠀⠀⠀⠀⠀⠀⣀⣀⣤⣤⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠸⡇⠀⠿⡀⠀⠀⠀⣀⡴⢿⣿⣿⣿⣿⣿⣿⣷⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠑⢄⣠⠾⠁⣀⣄⡈⠙⣿⣿⣿⣿⣿⣿⣿⣆⠀⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⢀⡀⠁⠀⠀⠈⠙⠛⠂⠈⣿⣿⣿⣿⣿⠿⡿⢿⣆⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⢀⡾⣁⣀⠀⠴⠂⠙⣗⡀⠀⢻⣿⣿⠭⢤⣴⣦⣤⣹⠀⠀⠀⢀⢴⣶⣆\n";
  std::cout << "⠀⠀⢀⣾⣿⣿⣿⣷⣮⣽⣾⣿⣥⣴⣿⣿⡿⢂⠔⢚⡿⢿⣿⣦⣴⣾⠁⠸⣼⡿\n";
  std::cout << "⠀⢀⡞⠁⠙⠻⠿⠟⠉⠀⠛⢹⣿⣿⣿⣿⣌⢤⣼⣿⣾⣿⡟⠉⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⣾⣷⣶⠇⠀⠀⣤⣄⣀⡀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠉⠈⠉⠀⠀⢦⡈⢻⣿⣿⣿⣶⣶⣶⣶⣤⣽⡹⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠉⠲⣽⡻⢿⣿⣿⣿⣿⣿⣷⣜⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣷⣶⣮⣭⣽⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⣀⣀⣈⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠹⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠟⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀\n";
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⠻⠿⠿⠿⠛⠉";

}

