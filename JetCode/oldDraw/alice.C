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
// => DrawVoigt("Inclusive", 1, "hUSS;", "hLSS;")
//=========================================

void Shreck();
void DrawVoigt(TString Kind, TString DataNum, TString USS, TString LSS)
{
  std::cout.precision(6); //keep the rest of the code unchanged
  gStyle->SetOptStat(0);

  const Int_t nPtBins = 9;
  Double_t xBins[nPtBins+1];
  TLatex* SpecialLatex_plot = nullptr;
  TLatex* SpecialLatex_fit = nullptr;

  if (Kind == "Outside"){
    Double_t outside_xBins[nPtBins+1] =
      {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(outside_xBins), end(outside_xBins), xBins);
  }
  else if (Kind == "Inclusive"){
    Double_t inclusive_xBins[nPtBins+1] =
      {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(inclusive_xBins), end(inclusive_xBins), xBins);

    SpecialLatex_plot = new TLatex(0.547, 0.243, "Min. Bias");
    SpecialLatex_plot->SetNDC(kTRUE);
    SpecialLatex_plot->SetTextSize(0.04);
    SpecialLatex_plot->SetTextFont(42);
    SpecialLatex_plot->SetTextColor(1);
    SpecialLatex_plot->Draw("SAME");
    
    SpecialLatex_fit = new TLatex(0.557, 0.79, "Min. Bias");
    SpecialLatex_fit->SetNDC(kTRUE);
    SpecialLatex_fit->SetTextSize(0.04);
    SpecialLatex_fit->SetTextFont(42);
    SpecialLatex_fit->SetTextColor(1);
    SpecialLatex_fit->Draw("SAME");
    
  }
  else if (Kind == "Inside"){
    Double_t inside_xBins[nPtBins+1] = 
      {0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.2, 4.0, 5.0, 8.0};
    copy(begin(inside_xBins), end(inside_xBins), xBins);

    SpecialLatex_plot = new TLatex(0.547, 0.243, "#Delta R_{jet}^{K*}< 0.4,  p_{T}^{jet} : 0-0.4%");
    SpecialLatex_plot->SetNDC(kTRUE);
    SpecialLatex_plot->SetTextSize(0.04);
    SpecialLatex_plot->SetTextFont(42);
    SpecialLatex_plot->SetTextColor(1);
    SpecialLatex_plot->Draw("SAME");
    
    SpecialLatex_fit = new TLatex(0.558, 0.785, "#Delta R_{jet}^{K*}< 0.4,  p_{T}^{jet} : 0-0.4%");
    SpecialLatex_fit->SetNDC(kTRUE);
    SpecialLatex_fit->SetTextSize(0.04);
    SpecialLatex_fit->SetTextFont(42);
    SpecialLatex_fit->SetTextColor(1);
    SpecialLatex_fit->Draw("SAME");
  }
  else{
    cout << "Please specify either 'Inside', 'Outside', or 'Inclusive'." << endl;
    return;
  }

  Double_t imlow = 0.65;
  Double_t imtop = 1.15;
  Double_t nBins = 400.0;

  TString widelocal = "tem/";
  TString localpath = widelocal + Kind; 

  //=========================================
  // Data file inputs
  //=========================================

  TString mainfile = "../kstar-in-jets/Result_rootfile/alice1.root";
  
  TFile* data = TFile::Open(mainfile);
  if(!data)   return;
  //=========================================
  //| Information : Data names              | 
  //| combined data ;1  single data ;2      |
  //| hUSS_INSIDE;1 hUSS_OUTSIDE;1  hUSS;1  |
  //=========================================
  TString SignalName("kstar-in-jets/"+USS+DataNum);
  TString BackgroundName("kstar-in-jets/"+LSS+DataNum);
  if(DataNum != 1 && DataNum != 2) return;
  
  TH3D* hLSS = (TH3D*)data->Get(BackgroundName);
  TH3D* hUSS = (TH3D*)data->Get(SignalName);

  //=========================================
  // FJptHistogram
  //=========================================

  TH1D* FJpt = (TH1D*)data->Get("kstar-in-jets/FJptHistogram");
  TCanvas* FJc = new TCanvas("FJc","fastjetpt", 940, 880);
  FJc->SetLeftMargin(0.15);
  FJc->SetRightMargin(0.05);
  FJc->SetTopMargin(0.10);
  FJc->SetBottomMargin(0.15);
  
  FJpt->Sumw2();
  FJpt->Draw();
  
  FJpt->GetXaxis()->SetTitle("#it{p}_{T}[GeV/#it{c}]");
  FJpt->GetXaxis()->SetTitleSize(0.06);
  FJpt->GetXaxis()->SetTitleFont(62);
  FJpt->GetXaxis()->SetTitleOffset(1.0);

  FJpt->GetYaxis()->SetTitle("#it{N }_{jet}^{Raw}");
  FJpt->GetYaxis()->SetTitleSize(0.06);
  FJpt->GetYaxis()->SetTitleFont(62);
  FJpt->GetYaxis()->SetTitleOffset(1.0);
  
  FJpt->SetTitle("");
  
  TLatex* textpta = new TLatex(0.526, 0.915, "ALICE #bf{Performance}");
  textpta->SetNDC(kTRUE);
  textpta->SetTextSize(0.05);
  textpta->SetTextFont(62);
  textpta->SetTextColor(1);
  textpta->Draw("SAME");
  TLatex* textptb = new TLatex(0.526, 0.854, "pp @ 13.6 TeV");
  textptb->SetNDC(kTRUE);
  textptb->SetTextSize(0.05);
  textptb->SetTextFont(42);
  textptb->SetTextColor(1);
  textptb->Draw("SAME");
  TLatex* textptc = new TLatex(0.526, 0.790, "#it{p}_{T}^{jet} > 8 GeV/#it{c}");
  textptc->SetNDC(kTRUE);
  textptc->SetTextSize(0.05);
  textptc->SetTextFont(42);
  textptc->SetTextColor(1);
  textptc->Draw("SAME");
  TLatex* textptd = new TLatex(0.526, 0.726, "1.05695 x 10^{8} Events");
  textptd->SetNDC(kTRUE);
  textptd->SetTextSize(0.05);
  textptd->SetTextFont(42);
  textptd->SetTextColor(1);
  textptd->Draw("SAME");
  FJc->SetLogy();
  FJc->SaveAs("../kstar-in-jets/Plot/"+widelocal+"FastJetpt.png","RECREATE");
  
  //=========================================
  // 1-D projection plots
  //=========================================
  TH1D* hLSS_Y [nPtBins];
  TH1D* hUSS_Y [nPtBins];
  
  TH1* hUSS_Y_Bkg [nPtBins];
  TH1* hUSS_Y_Peak [nPtBins];
  
  TH1D* hFitVData [nPtBins];
  
  for (Int_t jpt=0; jpt<nPtBins; jpt++)
    {
      hLSS_Y[jpt]=new TH1D(Form("hLSS_Y_%i", jpt), Form("hLSS_Y_%i", jpt), nBins, imlow, imtop);
      hUSS_Y[jpt]=new TH1D(Form("hUSS_Y_%i", jpt), Form("hUSS_Y_%i", jpt), nBins, imlow, imtop);
      
      hUSS_Y_Bkg[jpt]=new TH1D(Form("hUSS_Y_Bkg_%i", jpt), Form("hUSS_Y_Bkg_%i", jpt), nBins, imlow, imtop);
      hUSS_Y_Peak[jpt]=new TH1D(Form("hUSS_Y_Peak_%i", jpt), Form("hUSS_Y_Peak_%i", jpt), nBins, imlow, imtop);

      hFitVData[jpt] = new TH1D(Form("hFitVData_%i",jpt), "(Fit-Data)/Peak Value", nBins, imlow, imtop);
      hFitVData[jpt]->Sumw2();
    }

  //=========================================
  // Loop over pT bins
  //=========================================
  TH1D* hLSS_Y_COPY [nPtBins];
  TH1D* hUSS_Y_COPY [nPtBins];
  
  Double_t ReBinNum = 4;
  Double_t ptlow = 0;
  Double_t pthigh = 0;

  /* TCanvas* plots = new TCanvas("plots", "Combined Plots", 1280, 1080); */
  /* plots->Divide(3,3); */
  /* TCanvas* substraction = new TCanvas("substraction", "Combined substractions", 1280, 1080); */
  /* substraction->Divide(3,3); */
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
      //plots->cd(ipt +1);
      TCanvas* plots = new TCanvas(Form("plots_%d",ipt), "", 940, 800);
      
      plots->SetLeftMargin(0.14);
      plots->SetRightMargin(0.08);
      plots->SetTopMargin(0.08);
      plots->SetBottomMargin(0.14);
      
      hUSS_Y_COPY[ipt] = (TH1D*)hUSS_Y[ipt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[ipt] = (TH1D*)hLSS_Y[ipt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hUSS_Y[ipt]->Add(hLSS_Y[ipt], -1);
      
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitle("M_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]");
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitle("Counts");
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->SetTitle(Form(""));
      hUSS_Y_COPY[ipt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[ipt]->Draw();
      hUSS_Y_COPY[ipt]->SetLineColor(kRed);
      hUSS_Y_COPY[ipt]->SetMarkerStyle(26);
      hUSS_Y_COPY[ipt]->SetMarkerSize(1.2);
      hUSS_Y_COPY[ipt]->SetMarkerColor(kRed); 
      hLSS_Y_COPY[ipt]->Draw("SAME");
      hLSS_Y_COPY[ipt]->SetLineColor(kBlue);
      hLSS_Y_COPY[ipt]->SetMarkerStyle(4);
      hLSS_Y_COPY[ipt]->SetMarkerSize(1.2);
      hLSS_Y_COPY[ipt]->SetMarkerColor(kBlue);

      TLatex* txtInfoA = new TLatex(0.55, 0.35, Kind+" yield");
      txtInfoA->SetNDC(kTRUE); // use NDC
      txtInfoA->SetTextSize(0.04);
      txtInfoA->SetTextFont(62);
      txtInfoA->SetTextColor(1);
      txtInfoA->Draw("SAME");
      TLatex* txtInfoB = new TLatex(0.55, 0.40, "ALICE #bf{Performance}");
      txtInfoB->SetNDC(kTRUE); // use NDC
      txtInfoB->SetTextSize(0.05);
      txtInfoB->SetTextFont(62);
      txtInfoB->SetTextColor(1);
      txtInfoB->Draw("SAME");
      TLatex* txtInfoC = new TLatex(0.55, 0.30, "pp #sqrt{s} = 13.6 TeV");
      txtInfoC->SetNDC(kTRUE); // use NDC
      txtInfoC->SetTextSize(0.04);
      txtInfoC->SetTextFont(42);
      txtInfoC->SetTextColor(1);
      txtInfoC->Draw("SAME");
      if (SpecialLatex_plot != nullptr) {
	SpecialLatex_plot->Draw("SAME");
      }
      TLatex* txtInfoE = new TLatex(0.543, 0.185, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh));
      txtInfoE->SetNDC(kTRUE); // use NDC
      txtInfoE->SetTextSize(0.04);
      txtInfoE->SetTextFont(42);
      txtInfoE->SetTextColor(1);
      txtInfoE->Draw("SAME");
      TLegend* legend = new TLegend(0.20,0.17,0.45,0.29);
      legend->SetFillColor(4000);
      legend->SetBorderSize(0);
      legend->SetTextSize(0.05);
      legend->AddEntry(hUSS_Y_COPY[ipt], "Signal", "lpf");
      legend->AddEntry(hLSS_Y_COPY[ipt], "Background", "lpf");
      legend->Draw("SAME");

      //plots->SaveAs("../kstar-in-jets/Plot/" + localpath + "/combined_plots.png");
      plots->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/projection_%1.1f_%1.1f.png", ptlow, pthigh),"RECREATE");
      
      //=========================================
      // Draw 1-D projection hUSS-hLSS 
      //=========================================
      
      TCanvas* substraction = new TCanvas(Form("substraction_%d",ipt), "", 550, 440);
      
      substraction->SetLeftMargin(0.12);
      substraction->SetRightMargin(0.1);
      substraction->SetTopMargin(0.08);
      substraction->SetBottomMargin(0.14);
      //substraction->cd(ipt +1);
      
      hUSS_Y_COPY[ipt] = (TH1D*)hUSS_Y[ipt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[ipt] = (TH1D*)hLSS_Y[ipt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitle("M_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]");
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitle("Counts");
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[ipt]->GetYaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[ipt]->SetTitle(Form(""));
      hUSS_Y_COPY[ipt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[ipt]->Draw();
      TLatex* txtInfoa = new TLatex(0.55, 0.88, Kind+" yield");
      txtInfoa->SetNDC(kTRUE); // use NDC
      txtInfoa->SetTextSize(0.04);
      txtInfoa->SetTextFont(62);
      txtInfoa->SetTextColor(1);
      txtInfoa->Draw("SAME");
      TLatex* txtInfob = new TLatex(0.55, 0.93, "ALICE #bf{Performance}");
      txtInfob->SetNDC(kTRUE); // use NDC
      txtInfob->SetTextSize(0.05);
      txtInfob->SetTextFont(62);
      txtInfob->SetTextColor(1);
      txtInfob->Draw("SAME");
      TLatex* txtInfoc = new TLatex(0.55, 0.83, "pp #sqrt{s} = 13.6 TeV");
      txtInfoc->SetNDC(kTRUE); // use NDC
      txtInfoc->SetTextSize(0.04);
      txtInfoc->SetTextFont(42);
      txtInfoc->SetTextColor(1);
      txtInfoc->Draw("SAME");
      if (SpecialLatex_fit != nullptr) {
	SpecialLatex_fit->Draw("SAME");
      }
      TLatex* txtInfoe = new TLatex(0.55, 0.71, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh));
      txtInfoe->SetNDC(kTRUE); // use NDC
      txtInfoe->SetTextSize(0.04);
      txtInfoe->SetTextFont(42);
      txtInfoe->SetTextColor(1);
      txtInfoe->Draw("SAME");    

      substraction->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/reduced_projection_%1.1f_%1.1f.png", ptlow, pthigh));
      // substraction->SaveAs("../kstar-in-jets/Plot/" + localpath + "/combined_substractions.png");
    }
  
  Double_t Siglow = 0.85;
  Double_t Sighigh = 0.95;
  Double_t Bglow = 0.65;
  Double_t Bghigh = 1.15;

  TF1* fBreit [nPtBins];
  TF1* fBreitCopy [nPtBins];
  TF1* fBg [nPtBins];
  TF1* fFit [nPtBins];
  
  TH1D* hYield_Bin = new TH1D("binYield","binYield", nPtBins,xBins);
  TH1D* h_binYield = new TH1D("fitYield","fitYield;pT;N", nPtBins,xBins);  

  TF1* fBreit_only [nPtBins];
  TF1* fBg_only [nPtBins];

  /* TCanvas* lines = new TCanvas("lines", "Combined lines", 1280, 1080); */
  /* lines->Divide(3,3); */
  for(int fpt=0; fpt<nPtBins; fpt++) //nPtBins+1 MUST DELETE +! When i wanna draw one histogram
    {
      hUSS_Y_Bkg[fpt]->Add(hUSS_Y[fpt]);
      hUSS_Y_Peak[fpt]->Add(hUSS_Y[fpt]);
      
      for(int bin = hUSS_Y_Bkg[fpt]->FindBin(Siglow); bin<=hUSS_Y_Bkg[fpt]->FindBin(Sighigh); bin++)
	hUSS_Y_Bkg[fpt]->SetBinError(bin,0);
      
      ptlow = xBins[fpt];
      pthigh = xBins[fpt+1];
      
      //=========================================
      // Peak Function 
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
      
      //add
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
      hUSS_Y_Peak[fpt]->Fit(fBreit[fpt], "0", "", Siglow, Sighigh);
          
      fBreitCopy[fpt] = new TF1(Form("BwFit_Copy_%i",fpt),"[0]*TMath::Voigt(x-[2], [1], [3])", Siglow, Sighigh);
      fBreitCopy[fpt]->SetLineColor(2);

      //=========================================
      // Background Function 
      //=========================================
      fBg[fpt] = new TF1(Form("BkgFit_%i",fpt), "pol3", Bglow, Bghigh); //you have to be careful! This is 4 par

      fBg[fpt]->SetLineColor(1); //Black:1
      fBg[fpt]->SetLineStyle(2);
      fBg[fpt]->SetLineWidth(4);
      fBg[fpt]->SetRange(Bglow, Bghigh);
      
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);
      hUSS_Y_Bkg[fpt]->Fit(fBg[fpt], "0", "", Bglow, Bghigh);    

      //=========================================
      // Combine Peak + Background
      //=========================================
      fFit[fpt] = new TF1(Form("FullFit_%i",fpt), "[0]*TMath::Voigt(x-[2], [1], [3]) + pol3(4)", Bglow, Bghigh);
      fFit[fpt]->SetLineColor(2);
      // SIGNAL
      fFit[fpt]->SetParameter(0, fBreit[fpt]->GetParameter(0));
      fFit[fpt]->SetParameter(1, fBreit[fpt]->GetParameter(1));
      fFit[fpt]->SetParameter(2, fBreit[fpt]->GetParameter(2));
      fFit[fpt]->SetParameter(3, fBreit[fpt]->GetParameter(3));
      
      // BACKGROUND
      fFit[fpt]->SetParameter(4, fBg[fpt]->GetParameter(0));
      fFit[fpt]->SetParameter(5, fBg[fpt]->GetParameter(1));
      fFit[fpt]->SetParameter(6, fBg[fpt]->GetParameter(2));
      fFit[fpt]->SetParameter(7, fBg[fpt]->GetParameter(3));
      
      // OK, WE NOW HAVE GOOD STARTING VALUES
      hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Bglow, Bghigh);
      hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Bglow, Bghigh);
      hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Bglow, Bghigh);
      
      fBg[fpt]->SetRange(Bglow, Bghigh);

      fBreit[fpt]->SetParameters(fFit[fpt]->GetParameters());
      fBreitCopy[fpt]->SetParameters(fFit[fpt]->GetParameters());

      //==========================
      // Draw
      //==========================   
      TCanvas * lines = new TCanvas(Form("line_%d",fpt),"" , 940, 800);

      lines->SetLeftMargin(0.12);
      lines->SetRightMargin(0.1);
      lines->SetTopMargin(0.08);
      lines->SetBottomMargin(0.14);
      // lines->cd(fpt +1);

      hUSS_Y_COPY[fpt] = (TH1D*)hUSS_Y[fpt]->Clone(Form("hUSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      hLSS_Y_COPY[fpt] = (TH1D*)hLSS_Y[fpt]->Clone(Form("hLSS_Y_Minv_%1.1f_%1.1f", ptlow, pthigh));
      
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitle("M_{inv}^{K^{+}#pi^{-}}[GeV/#it{c}^{2}]");
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[fpt]->GetXaxis()->SetTitleOffset(1.0);
      
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitle("Entries");
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleSize(0.06);
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleFont(62);  
      hUSS_Y_COPY[fpt]->GetYaxis()->SetTitleOffset(1.0);
      hUSS_Y_COPY[fpt]->SetTitle(Form(""));
      hUSS_Y_COPY[fpt]->SetTitleSize(0.05, "t");

      hUSS_Y_COPY[fpt]->SetLineColor(1);
      hUSS_Y_COPY[fpt]->Draw();
      
      fFit[fpt]->SetLineColor(2);
      fFit[fpt]->Draw("SAME");

      fBg[fpt]->SetParameter(0, fFit[fpt]->GetParameter(4));
      fBg[fpt]->SetParameter(1, fFit[fpt]->GetParameter(5));
      fBg[fpt]->SetParameter(2, fFit[fpt]->GetParameter(6));
      fBg[fpt]->SetParameter(3, fFit[fpt]->GetParameter(7));

      fBg[fpt]->SetLineColor(4);
      fBg[fpt]->Draw("SAME");

      TFitResultPtr ptr = hUSS_Y[fpt]->Fit(fFit[fpt], "SR0", "", Bglow, Bghigh);
      double Chi2 = ptr->Chi2();
      int NDF = ptr->Ndf();
      double chi2Ndf = Chi2 / NDF;
      cout << "chi2/ndf :"<<chi2Ndf<<endl;
      string chi2NdfText = Form("#chi^{2} / NDF = %.2f", chi2Ndf);
      
      TLatex* txtInfoAA = new TLatex(0.55, 0.88, Kind+" yield");
      txtInfoAA->SetNDC(kTRUE); // use NDC
      txtInfoAA->SetTextSize(0.04);
      txtInfoAA->SetTextFont(62);
      txtInfoAA->SetTextColor(1);
      txtInfoAA->Draw("SAME");
      TLatex* txtInfoBB = new TLatex(0.55, 0.93, "ALICE #bf{Performance}");
      txtInfoBB->SetNDC(kTRUE); // use NDC
      txtInfoBB->SetTextSize(0.05);
      txtInfoBB->SetTextFont(62);
      txtInfoBB->SetTextColor(1);
      txtInfoBB->Draw("SAME");
      TLatex* txtInfoCC = new TLatex(0.56, 0.84, "pp #sqrt{s} = 13.6 TeV");
      txtInfoCC->SetNDC(kTRUE); // use NDC
      txtInfoCC->SetTextSize(0.04);
      txtInfoCC->SetTextFont(42);
      txtInfoCC->SetTextColor(1);
      txtInfoCC->Draw("SAME");
      if (SpecialLatex_fit != nullptr) {
	SpecialLatex_fit->Draw("SAME");
      }
      TLatex* txtInfoEE = new TLatex(0.56, 0.74, Form("%.1f < p_{T}^{K*}< %.1f GeV/#it{c}", ptlow, pthigh));
      txtInfoEE->SetNDC(kTRUE); // use NDC
      txtInfoEE->SetTextSize(0.04);
      txtInfoEE->SetTextFont(42);
      txtInfoEE->SetTextColor(1);
      txtInfoEE->Draw("SAME");
      TLatex* txtInfoFF = new TLatex(0.56, 0.69, chi2NdfText.c_str());
      txtInfoFF->SetNDC(kTRUE); // use NDC
      txtInfoFF->SetTextSize(0.04);
      txtInfoFF->SetTextFont(42);
      txtInfoFF->SetTextColor(1);
      txtInfoFF->Draw("SAME");
      TLegend *legendd = new TLegend(0.15,0.79,0.26,0.89);
      legendd->SetFillColor(4000);
      legendd->SetBorderSize(0);
      legendd->SetTextSize(0.03);
      legendd->AddEntry(hUSS_Y_COPY[fpt],"Data. M_{inv}^{K^{+}#pi^{-}}", "lpf");
      legendd->AddEntry(fBg[fpt], "Bkg. Fit", "lpf");
      legendd->AddEntry(fFit[fpt], "Voigt. Fit", "lpf");
      legendd->Draw("SAME");
      
      lines->SaveAs(Form("../kstar-in-jets/Plot/"+localpath+"/two_fit_%1.1f_%1.1f.png", ptlow, pthigh));
      //lines->SaveAs("../kstar-in-jets/Plot/" + localpath + "/combined_lines.png");

      //==========================
      // Reflecting an error
      //==========================

      /*     //================================ */
      /*     // Draw # of bin of (peak - bkg) */
      /*     //================================ */
  
      /*     TCanvas* numbin = new TCanvas(); */
      /*     h_binYield->Draw(); */
      /*     h_binYield->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]"); */
      /*     h_binYield->GetYaxis()->SetTitle("Counts"); */
      /*     h_binYield->SetTitle("#frac{d^{2}N_{#K^*}}{dp_{T}d#eta}"); */
      /*     h_binYield->Draw(); */
      /*     numbin->SaveAs("../phi-in-jets/Plot/"+localpath+"/h_binYield.png","RECREATE"); */

      /*     Double_t Peak = 0; */
      /*     for(Int_t bin = 1; bin<=nBins; bin++) */
      /* 	{ */
      /* 	  Double_t Xaxis = hUSS_Y[fpt]->GetXaxis()->GetBinCenter(bin); */
      /* 	  if(Xaxis>=Bglow) */
      /* 	    { */
      /* 	      Double_t Data = hUSS_Y[fpt]->GetBinContent(bin); */
      /* 	      Double_t Fit = fFit[fpt]->Eval(Xaxis); */
      /* 	      Double_t Ratio = Fit-Data; */
      /* 	      Double_t erData = hUSS_Y[fpt]->GetBinError(bin); */
      /* 	      if(Data>Peak) */
      /* 		Peak = Data; */
      /* 	      hFitVData[fpt]->SetBinContent(bin, Ratio); */
      /* 	    } */
      /* 	  if(Xaxis>Bghigh) */
      /* 	    break; */
      /* 	} */

      /*     Double_t binW = hUSS_Y[fpt]->GetBinWidth(1); //any number */
      /*     Double_t Sint = fBreit[fpt]->Integral(Siglow, Sighigh); */
      /*     Double_t yield = Sint/binW; */
      /*     Double_t Tails = 0; */
      /*     Double_t binyield = 0; */
      /*     Double_t binpiece = 0; */
  
      /*     for(Int_t binz = hUSS_Y[fpt]->FindBin(Siglow); binz < hUSS_Y[fpt]->FindBin(Sighigh); binz++) */
      /* 	{ */
      /* 	  binpiece += fBg[fpt]->Eval(hUSS_Y[fpt]->GetBinCenter(binz)); */
      /* 	} */
      /*     if(binpiece<0) */
      /* 	binpiece=0; */
      
      /*     Double_t YieldError = 0; */
      /*     binyield = hUSS_Y[fpt]->IntegralAndError(hUSS_Y[fpt]->FindBin(Siglow), hUSS_Y[fpt]->FindBin(Sighigh)-1, YieldError) - binpiece; */

      /*     TFitResultPtr ptr; */
      /*     ptr = hUSS_Y[fpt]->Fit(fFit[fpt], "SR", "", Siglow, Sighigh); */

      /*     TMatrixDSym covTot(fFit[fpt]->GetNpar()); */
      /*     TMatrixDSym covGG(fBreit[fpt]->GetNpar()); */
      /*     TMatrixDSym covGG_err(fBg[fpt]->GetNpar()); */
      /*     covTot = ptr->GetCovarianceMatrix(); */
      /*     covGG = covTot.GetSub(0,3,0,3); */
      /*     covGG_err = covTot.GetSub(4,7,4,7); */
 
      /*     Double_t error_fit = (1/binW)*(fBreit[fpt]->IntegralError(Siglow,Sighigh,fBreit[fpt]->GetParameters(),covGG.GetMatrixArray())); */
      /*     Double_t Bg_err = (1/binW)*(fBg[fpt]->IntegralError(Siglow,Sighigh,fBg[fpt]->GetParameters(),covGG_err.GetMatrixArray())); */
      
      /*     Tails = TMath::Abs(fBreitCopy[fpt]->Integral(0.749735,Siglow)/binW); */
      /*     Tails += TMath::Abs(fBreitCopy[fpt]->Integral(Sighigh,TMath::Infinity())/binW); */
      
      /*     //Double_t Error = TMath::Sqrt(YieldError*YieldError +  error_fit*error_fit + Bg_err*Bg_err); */
      /*     Double_t Error = TMath::Sqrt(YieldError*YieldError); */
      
      /*     binyield += Tails; */
      /*     yield += Tails; */
      
      /*     h_binYield->SetBinContent(fpt+1, binyield/h_binYield->GetXaxis()->GetBinWidth(fpt+1)); */
      /*     h_binYield->SetBinError(fpt+1,Error/h_binYield->GetXaxis()->GetBinWidth(fpt+1)); */
      
      /*   } */
  
      /* TString Eventdata("phi-in-jets/nEvents;"+DataNum); */
      /* TH1D* hVtx = (TH1D*)data->Get(Eventdata); */
      /* Double_t nEvents = hVtx->GetBinContent(1); */
      /* TH1D* hCorr_binYield = new TH1D("hCorr_binYield", "hCorr_binYield", nPtBins,xBins); */
  
      /* hCorr_binYield->Add(h_binYield); */
      /* hCorr_binYield->Scale(1.0/nEvents); */
    
      /* // hCorr_binYield->GetYaxis()->SetTitle("#frac{d^{2}N_{K*}}{dp_{T} dy}"); */
      /* hCorr_binYield->GetXaxis()->SetTitle("p_{T}[GeV/c]"); */

      /* for(int i=1; i<10; i++) */
      /*   { */
      /*     hCorr_binYield->SetBinContent(i, hCorr_binYield->GetBinContent(i)/ hCorr_binYield->GetBinWidth(i)); */
      /*     hCorr_binYield->SetBinError(i, hCorr_binYield->GetBinError(i)/ hCorr_binYield->GetBinWidth(i)); */
      /*   } */
  
      /* TCanvas* logbin = new TCanvas(); */
      /* hCorr_binYield->Draw(); */
      /* logbin->SetLogy(); */
      /* logbin->SaveAs("../phi-in-jets/Plot/"+localpath+"/logbinhYield.png"); */

      /* TFile* outFile=new TFile("outfile_outside.root", "RECREATE"); */
      /* outFile->cd(); */
      /* hCorr_binYield->Write(); */
      /* outFile->Write(); */
      /* outFile->Close(); */
  
    }
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

