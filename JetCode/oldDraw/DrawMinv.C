void DrawMinv()
{
  std::cout.precision(6);
  gStyle->SetOptStat(0);

  const Int_t nPtBins = 12;
  Double_t xBins [nPtBins+1] = { 0.8, 1.2, 1.9, 2.1, 2.3, 2.8, 3.0, 3.2, 3.4, 3.9,
				 4.2, 4.4, 8.0};
  Double_t imlow = 0.75;
  Double_t imtop = 1.15;
  Double_t nBins = 400.0;

  TString firstpath = "../phi-in-jets/Result_rootfile/pidtest2.root";
  TString secondpath = "../phi-in-jets/Result_rootfile/output_totcombined.root";
  TString localpath = "PID2test";
  
  //==========================
  // Data file inputs
  //==========================
  TFile* dataA = TFile::Open(firstpath);
  TFile* dataB = TFile::Open(secondpath);
  if(!dataA || !dataB)  return;
  
  TString SignalName("phi-in-jets/hUSS;2");
  TString BackgroundName("phi-in-jets/hLSS;2");
  TH3D* hLSSA = (TH3D*)dataA->Get(BackgroundName);
  TH3D* hUSSA = (TH3D*)dataA->Get(SignalName);

  TString SignalNameB("phi-in-jets/hUSS;1");
  TString BackgroundNameB("phi-in-jets/hLSS;1");
  TH3D* hLSSB = (TH3D*)dataB->Get(BackgroundNameB);
  TH3D* hUSSB = (TH3D*)dataB->Get(SignalNameB);

  //==========================
  // 1-D projection plots
  //==========================
  TH1D* hLSS_YA[nPtBins];
  TH1D* hUSS_YA[nPtBins];
  TH1D* hLSS_YB[nPtBins];
  TH1D* hUSS_YB[nPtBins];

  for(Int_t i=0; i<nPtBins; i++)
  {
    hLSS_YA[i] = new TH1D(Form("hLSS_Y_%i", i), Form("hLSS_Y_%i", i), 400, 0.75, 1.15);
    hUSS_YA[i] = new TH1D(Form("hUSS_Y_%i", i), Form("hUSS_Y_%i", i), 400, 0.75, 1.15);
    hLSS_YA[i]->Sumw2();
    hUSS_YA[i]->Sumw2();
    
    hLSS_YB[i] = new TH1D(Form("hLSS_YB_%i", i), Form("hLSS_YB_%i", i), 400, 0.75, 1.15);
    hUSS_YB[i] = new TH1D(Form("hUSS_YB_%i", i), Form("hUSS_YB_%i", i), 400, 0.75, 1.15);
    hLSS_YB[i]->Sumw2();
    hUSS_YB[i]->Sumw2();
  }

  Double_t ptlow = 0;
  Double_t pthigh = 0;

  TH1F* hDummy = new TH1F("hDummy", "hDummy", 250, 0.7, 3);
  hDummy->SetTitle("");
  hDummy->GetXaxis()->SetTitle("#it{M}_{inv} (GeV/#it{c}^{2})");
  hDummy->GetYaxis()->SetTitle("1/N dN/d#it{M}_{inv}");
  hDummy->GetXaxis()->SetTitleSize(0.07);
  hDummy->GetXaxis()->SetTitleOffset(0.8);
  hDummy->GetYaxis()->SetTitleSize(0.07);
  hDummy->GetYaxis()->SetTitleOffset(0.85);
  hDummy->GetYaxis()->SetRangeUser(0.5, 4.5);
  hDummy->GetXaxis()->SetRangeUser(0.75, 1.15);

  TCanvas* cP;
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.01);

  Int_t Nrebin = 2;

  //==========================
  // Loop over pT bins
  //==========================
  for(Int_t ipt=0; ipt<1; ipt++)//nPtBins; a++)
  {
    ptlow =xBins[ipt];
    pthigh = xBins[ipt+1];

    //==========================
    // PID cut 1
    //==========================
    hUSSA->GetYaxis()->SetRangeUser(ptlow, pthigh);
    hLSSA->GetYaxis()->SetRangeUser(ptlow, pthigh);
       
    hLSS_YA[ipt]->Add(hLSSA->Project3D("z o"),1);
    hUSS_YA[ipt]->Add(hUSSA->Project3D("z o"),1);
    
    hUSS_YA[ipt]->Scale(1.0/hUSS_YA[ipt]->Integral(1,400));
    hUSS_YA[ipt]->Rebin(Nrebin);
    hUSS_YA[ipt]->Scale(1.0/hUSS_YA[ipt]->GetXaxis()->GetBinWidth(1));
    
    hLSS_YA[ipt]->Scale(1.0/hLSS_YA[ipt]->Integral(1,400));
    hLSS_YA[ipt]->Rebin(Nrebin);
    hLSS_YA[ipt]->Scale(1.0/hLSS_YA[ipt]->GetXaxis()->GetBinWidth(1));
    
    //hUSS_YA[ipt]->Add(hLSS_YA[ipt], -1);

    //==========================
    // PID cut 2
    //==========================
    hUSSB->GetYaxis()->SetRangeUser(ptlow, pthigh);
    hLSSB->GetYaxis()->SetRangeUser(ptlow, pthigh);
    
    hLSS_YB[ipt]->Add(hLSSB->Project3D("z o"),1);
    hUSS_YB[ipt]->Add(hUSSB->Project3D("z o"),1);
    
    hUSS_YB[ipt]->Scale(1.0/hUSS_YB[ipt]->Integral(1,400));
    hUSS_YB[ipt]->Rebin(Nrebin);
    hUSS_YB[ipt]->Scale(1.0/hUSS_YB[ipt]->GetXaxis()->GetBinWidth(1));
    
    hLSS_YB[ipt]->Scale(1.0/hLSS_YB[ipt]->Integral(1,400));
    hLSS_YB[ipt]->Rebin(Nrebin);
    hLSS_YB[ipt]->Scale(1.0/hLSS_YB[ipt]->GetXaxis()->GetBinWidth(1));
    
    //hUSS_YB[ipt]->Add(hLSS_YB[ipt], -1);

    //==========================
    // Draw
    //==========================
    cP = new TCanvas(Form("cP_%d", ipt), "", 550, 400);
    cP->cd();
    hDummy->Draw();
    
    hUSS_YA[ipt]->SetLineColor(kGreen+2);
    hUSS_YA[ipt]->Draw("same");

    hUSS_YB[ipt]->SetLineColor(kRed);
    hUSS_YB[ipt]->Draw("same");

    TLegend *legend = new TLegend(0.18,0.72,0.40,0.86);
    legend->SetFillColor(4000);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(hUSS_YB[ipt], "PID cut 1", "lpf");
    legend->AddEntry(hUSS_YA[ipt], "PID cut 2", "lpf");
    legend->Draw("same");
    
    TLatex* txta = new TLatex(0.18, 0.92, Form("%1.1f < p_{T} < %1.1f GeV/c", ptlow, pthigh));
    txta->SetNDC(kTRUE); // <- use NDC coordinate
    txta->SetTextSize(0.06);
    txta->SetTextFont(62);
    txta->SetTextColor(1);
    txta->Draw("same");
     
    cP->SaveAs(Form("../phi-in-jets/Plot/testpid/"+localpath+"/difPID_%f_%f.png", xBins[ipt], xBins[ipt+1]));

    // 첫 번째 히스토그램을 두 번째 히스토그램으로 나눕니다.
    //hUSS_YA[ipt]->Divide(hUSS_YB[ipt]);
    
    // 나눈 결과인 hUSS_YA[ipt]를 초록색으로 설정하여 그립니다.
    //hUSS_YA[ipt]->SetLineColor(kGreen);
    //hUSS_YA[ipt]->Draw("same");
  }
}


/* void Drawtotal(){ */
  
/*   TString fullpath = "../phi-in-jets/Result_rootfile/pidtest2.root"; */
/*   TString localpath = "PID2test"; */
/*   //////////////////////////////////////////////////////////////////////// */
/*   //////////////////////////////////////////////////////////////////////// */
/*   TFile* data = TFile::Open(fullpath); */
/*   if(!data) */
/*     return; */
/*   //////////////////////////////////////////////////////////////////////// */
/*   //////////////////////////////////////////////////////////////////////// */
/*   //hUSS_INSIDE;1 */
/*   //hUSS_OUTSIDE;1 */
/*   //hUSS;1 */
  
/*   TString SignalName("phi-in-jets/hUSS;1"); */
/*   TString BackgroundName("phi-in-jets/hLSS;1"); */
  
/*   TH3D* hLSS = (TH3D*)data->Get(BackgroundName); */
/*   TH3D* hUSS = (TH3D*)data->Get(SignalName); */

/*   // 초기값 0으로 start  */
/*   Double_t ptlow = 0; */
/*   Double_t pthigh = 0; */
  
  
/*   (hUSS->GetYaxis())->SetRangeUser(ptlow, pthigh); */
/*   (hLSS->GetYaxis())->SetRangeUser(ptlow, pthigh); */
  
  
  
/*   TString SignalName2("phi-in-jets/hUSS;2"); */
/*   TString BackgroundName2("phi-in-jets/hLSS;2"); */
  
/*   TH3D* hLSS2 = (TH3D*)data->Get(BackgroundName2); */
/*   TH3D* hUSS2 = (TH3D*)data->Get(SignalName2); */
 
/*   TH1D* hRatio = new TH1D("hRatio", "hRatio", 400, 0.75, 1.15); */
  
/*   //projection */
/*   hLSS_Y->Add(hLSS->Project3D("z o"),1); */
/*   hUSS_Y->Add(hUSS->Project3D("z o"),1); */
  
/*   //scale */
/*   hUSS_Y->Scale(1.0/hUSS_Y->Integral(1,400)); */
/*   hUSS_Y->Scale(1.0/hUSS_Y->GetXaxis()->GetBinWidth(1)); */

/*   hLSS_Y->Scale(1.0/hLSS_Y->Integral(1,400)); */
/*   hLSS_Y->Scale(1.0/hLSS_Y->GetXaxis()->GetBinWidth(1)); */

/*   hUSS->SetLineColor(kRed); */
/*   hUSS_Y->Add(hLSS_Y, -1); */
  
/*   //projection2 */
/*   hLSS_YB->Add(hLSS2->Project3D("z o"),1); */
/*   hUSS_YB->Add(hUSS2->Project3D("z o"),1); */
  
/*   //scale2 */
/*   hUSS_YB->Scale(1.0/hUSS_YB->Integral(1,400)); */
/*   hUSS_YB->Scale(1.0/hUSS_YB->GetXaxis()->GetBinWidth(1)); */

/*   hLSS_YB->Scale(1.0/hLSS_YB->Integral(1,400)); */
/*   hLSS_YB->Scale(1.0/hLSS_YB->GetXaxis()->GetBinWidth(1)); */

/*   hUSS2->SetLineColor(kBlue); */
/*   hUSS_YB->Add(hLSS_YB, -1); */

/*   //ratio */
/*   hRatio->Divide(hUSS_Y, hUSS_YB); */

  
/*   TCanvas *cp_inclusive = new TCanvas("name", "title", 440, 500); */

/*   // 첫 번째 히스토그램인 hUSS_YA[a]를 파란색으로 설정하여 그립니다. */
/*   hUSS_Y->Draw(); */
      
/*   // 두 번째 히스토그램인 hUSS_YB[a]를 빨간색으로 설정하여 겹쳐서 그립니다. */
/*   hUSS_YB->Draw("same"); */
      
      
/*   // 나눈 결과인 hUSS_Y[a]를 초록색으로 설정하여 그립니다. */
/*   hRatio->SetLineColor(kGreen); */
/*   hRatio->Draw("same"); */

/*   cp_inclusive->SaveAs("../phi-in-jets/Plot/Tot_combined2/"+localpath+"/pid2_Inclusive.png"); */

/* } */
