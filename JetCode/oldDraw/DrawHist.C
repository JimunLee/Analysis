#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TObject.h>
#include <TMath.h>
#include <TTimeStamp.h>
#include <TParticle.h>
#include <TLeaf.h>
#include <TLorentzVector.h>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>

using namespace std;

void DrawHist()
{  
  gStyle->SetOptStat(0); // hide the statistics box by 0 

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //// Step 1 ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  TString folpath = "/Users/jimun/Analysis/phi-in-jets/Plot/testpid2/" ;
  
  TFile* fin = TFile::Open("../phi-in-jets/Result_rootfile/pidtest2.root");

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // INCLUSIVE

  TH1D *hUSS = (TH1D*)fin->Get("phi-in-jets/hUSS_1D_2_3;2");
  TH1D *hLSS = (TH1D*)fin->Get("phi-in-jets/hLSS_1D_2_3;2");
  TCanvas *cp_inclusive = new TCanvas("name", "title", 440, 500);
  hUSS->SetLineColor(2); // Red:2
  hLSS->SetLineColor(4); // Blue:4
  hUSS->GetYaxis()->SetRangeUser(1500, 5000);
  hUSS->Draw();
  hLSS->Draw("same");
  cp_inclusive->SaveAs(folpath+"Inclusive_Graph.png");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // OUTSIDE
  
  /* TH1D *hUSS = (TH1D*)fin->Get("phi-in-jets/hUSS_OUTSIDE_1D"); */
  /* TH1D *hLSS = (TH1D*)fin->Get("phi-in-jets/hLSS_OUTSIDE_1D"); */

  /* // means of parenthesis : (name of canvas, title of canvas, width, height) */
  /* TCanvas *cp_outside = new TCanvas("name", "title", 440, 500); */
  
  /* // canvas3 (normalization) */
  /* /\* double norm = hUSS->Integral(350, 400)/hLSS->Integral(350,400); *\/ */
  /* /\* hLSS->Scale(norm); // h1->Scale(x)=h1->GetBincontetn(i) * x *\/ */
  
  /* // Color */
  /* hUSS->SetLineColor(2); // Red:2 */
  /* hLSS->SetLineColor(4); // Blue:4 */

  /* // canvas1 (Draw each one) */
  /* hUSS->GetYaxis()->SetRangeUser(0, 12000); */
  /* hUSS->Draw(); */
  /* hLSS->Draw("same"); //It can be used to draw over a canvas that has already been drawn. */

  /* //// file save */
  /* cp_outside->SaveAs(folpath+"outside_Graph.png"); // 1.Draw the graph each one */

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // INSIDE
  
  /* TH1D *hUSS = (TH1D*)fin->Get("phi-in-jets/hUSS_INSIDE_1D"); */
  /* TH1D *hLSS = (TH1D*)fin->Get("phi-in-jets/hLSS_INSIDE_1D"); */

  /* TCanvas *cp_inside = new TCanvas("name", "title", 440, 500); */
  /* hUSS->SetLineColor(2); // Red:2 */
  /* hLSS->SetLineColor(4); // Blue:4 */
  
  /* hUSS->GetYaxis()->SetRangeUser(2500, 5000); */
  /* hUSS->Draw(); */
  /* hLSS->Draw("same"); //It can be used to draw over a canvas that has already been drawn. */

  /* //// file save */
  /* cp_inside->SaveAs("/Users/jimun/Analysis/phi-in-jets/Plot/Tot_combined/inside_Graph.png"); */

  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////s////////////////////////////////
  //// Step 2 //////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  TH1D *hDummy = new TH1D("hDummy", "hDummy", 10, 0.75, 1.15);
  TLatex TL;
  hDummy->SetTitle("");
  //in
  // hDummy->GetYaxis()->SetRangeUser(0, 2000);
  //out
  hDummy->GetYaxis()->SetRangeUser(0, 1500);
  
  hDummy->GetYaxis()->SetTitle("N_{K^{+}\\pi^{-}}^{Pairs}");
  hDummy->GetXaxis()->SetTitle("#it{M}_{inv}");
  
  TH1D *hSubtract = (TH1D*)hUSS->Clone("USSCLone");
  hSubtract->Add(hLSS, -1);

  TCanvas *hSub = new TCanvas("hSub", "hSub", 1200, 1200);
  hSub->SetLeftMargin(0.15);
  hSub->SetBottomMargin(0.10);
  hDummy->Draw();
  hSubtract->Draw("same");

  //TLatex TL;
  TL.SetTextSize(0.04);
  TL.SetTextAlign(22);
  
  //INCLUSIVE
  TL.DrawLatexNDC(0.506678, 0.935961, "pp @ 13.6 TeV, Inclusive Jet Cone, R < 0.4");
  TL.DrawLatex(0.50667, 0.86, "#splitline{#it{p}_{T}^{jet}>10 GeV/c}{#bf{1.74014*e^{8}} Events}");
  hSub->SaveAs(folpath+"inclusive_USS-LSS.png");
  
  //OUT
  /* TL.DrawLatexNDC(0.506678, 0.935961, "pp @ 13.6 TeV, Outside Jet Cone, R < 0.4"); */
  /* TL.DrawLatex(0.50667, 0.86, "#splitline{#it{p}_{T}^{jet}>10 GeV/c}{#bf{1.74014*e^{8}} Events}"); */
  /* hSub->SaveAs("/Users/jimun/Analysis/phi-in-jets/Plot/Tot_combined/outside_USS-LSS.png"); */
  
  //IN
  /* TL.DrawLatexNDC(0.480801, 0.935319, "pp @ 13.6 TeV, Inside Jet Cone, R < 0.4"); */
  /* TL.DrawLatex(0.850749, 1832.98, "#splitline{#it{p}_{T}^{jet}>10 GeV/c}{#bf{1.74014*e^{8}} Events}"); */
  /* hSub->SaveAs("/Users/jimun/Analysis/phi-in-jets/Plot/Tot_combined/inside_USS-LSS.png"); */
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //// Step 3 //////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  /* TH1D *hJetPt = (TH1D*)fin->Get("phi-in-jets/FJptHistogram"); */
  /* hJetPt->Sumw2(); */
  /* TCanvas *FJpt = new TCanvas("FJpt", "FJpt", 1000, 1000); */
  /* FJpt->SetLogy(); */
  /* hJetPt->SetTitle(""); */
  /* hJetPt->GetYaxis()->SetTitle("N_{Jets}^{Raw}"); */
  /* hJetPt->GetXaxis()->SetTitle("#it{p}_{T} [GeV/c]"); */

  /* hJetPt->Draw(); */
 
  /* TL.SetTextSize(0.0483936); */
  /* TL.DrawLatex(110.02, 27274, "pp @ 13.6 TeV"); */
  /* TL.DrawLatex(109.77, 5508.67, "#splitline{#it{p}_{T}^{Jet}>10 GeV/c}{#bf{1.74014*e^{8 Event}}"); */

  /* //// file save */
  /* FJpt->SaveAs("/Users/jimun/Analysis/phi-in-jets/Plot/Tot_combined/FJpt_kaonstar.png"); */
  }
