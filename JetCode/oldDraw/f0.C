Double_t threshold = 0.139570*2.;
double PhaseSpaceFactor(double x, double T, double pT){
 return ( x / ( sqrt( pow(x,2) + pow(pT,2) ) ) ) * exp( -sqrt( pow(x,2) + pow(pT,2) )/T );
}

double Width(double x, double M0, double Gam0, int Spin){
 double w = pow( (x*x - threshold*threshold)/(M0*M0 - threshold*threshold),0.5+(double)Spin )*Gam0*M0/x;
 return w;
}

double breitWigner(double x, double Amp, double M0, double Gam0, int Spin){
 double br = Amp*x*M0*Width(x, M0, Gam0, Spin);
 br /= ( pow(M0*M0-x*x,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2) );
 return br;
}
double background(double x, double ind, double b1, double b2){
 double bg = pow(x-threshold,ind)*exp( b1*(x-threshold) + b2*pow(x-threshold,2) );
 return bg;
}
double breitWigner_f0(double *x, double* par){
 return breitWigner(x[0], par[1], par[0], par[2], 0);
}
double sfit(double *x, double* par){

 double bgA = par[0];
 double bgind = par[1];
 double bgb1 = par[2];
 double bgb2 = par[3];

 double f0m = par[4];
 double f0A = par[5];
 double f0g = par[6];
 int f0s = 0;

 double f2m = par[7];
 double f2A = par[8];
 double f2g = par[9];
 int f2s = 2;

 double rhom = par[10];
 double rhoA = par[11];
 double rhog = par[12];
 int rhos = 1;

 double T = par[13];
 double pT = par[14];

 double breitWigner_f0  = breitWigner(x[0],f0A,f0m,f0g,f0s);
 double breitWigner_f2  = breitWigner(x[0],f2A,f2m,f2g,f2s);
 double breitWigner_rho  = breitWigner(x[0],rhoA,rhom,rhog,rhos);

 double bgfunc = bgA*background(x[0],bgind,bgb1,bgb2);

 return ( breitWigner_f0 + breitWigner_f2 + breitWigner_rho )*PhaseSpaceFactor(x[0],T,pT) + bgfunc;
// return breitWigner_f0 + breitWigner_f2 + breitWigner_rho + bgfunc;
}


void FitLP(int LP=0, int Rebinning=1){
 TFile* fin = new TFile("InvMassOut.root","read");

 const int nmult = 4;
 const int npt = 4;
 const int nlpt = 2;
 const int nrt = 4;

 double m_min[nmult] = {
     0,  0,  20, 50 };
 double m_max[nmult] = {
     90, 20, 50, 80 };

 double p_min[npt] = {
    0.4, 1.0, 3.0, 0.0 };
 double p_max[npt] = {
    1.0, 3.0, 8.0, 8.0 };

 double lpt_min[nlpt] = {
    0, 5 };

 const char trnsName[nrt][1000] = {
	"Toward",
	"Transverse",
	"Opposite",
	"In-jet"};


 TH1D* hMassSig[nrt][nmult][npt];
 TF1* fitBin[nrt][nmult][npt];
 TF1* fitBinComp[nrt][nmult][npt][4];

 TFitResultPtr fitResults[nrt][nmult][npt];

 TGraphErrors* gf0Width[nrt][nmult];
 TGraphErrors* gf0Mass[nrt][nmult];
 TGraphErrors* gChi2OverNDF[nrt][nmult];
 for(int r=0;r<nrt;r++){
	for(int i=0;i<nmult;i++){
		gf0Width[r][i] = new TGraphErrors();
		gf0Mass[r][i] = new TGraphErrors();
		gChi2OverNDF[r][i] = new TGraphErrors();
	}
 }


 double Intgr_RawY_cntl[nrt][nmult][npt];
 double Intgr_RawY_stat[nrt][nmult][npt];
 TH1D* hPtStat[nrt][nmult];
 double PtBinnings[npt] = {
	0.4, 1.0, 3.0, 8.0 };

 int MagIndex[4] = {0,5,8,11};
 int lcolor[4] = {46, 30, 41, 28};

 double BGAmpMin = 1e2*2.0*1e1;
 double BGAmpMax = 3e7*2.0*1e1;

 double BGSlopeMin = -2.5;
 double BGSlopeMax = 3.0;

 double BGInd1Min = -5.0;
 double BGInd1Max = 5.0;

 double BGInd2Min = -0.5;
 double BGInd2Max = 5.0;

 double f0MassMinRange = 0.95;
 double f0MassMaxRange = 1.0;

 double f0AmpMin = 1.0*2.0;
 double f0AmpMax = 1e6*2.0;

 double f0WMin = 0.01;
 double f0WMax = 0.10;

 double f2MassMinRange = 1.2755-0.0012*0.002;
 double f2MassMaxRange = 1.2755+0.0012*0.002;

 double f2AmpMin = 1.0*2.0;
 double f2AmpMax = 1e6*2.0;

 double f2WMin = 0.1867-0.0024*0.02;
 double f2WMax = 0.1867+0.0029*0.02;

 double rhoAmpMin = 1.0*2.0;
 double rhoAmpMax = 1e6*2.0;

 double rhoMassMin = 0.7752;
 double rhoMassMax = 0.7754;

 double rhoWMin = 0.1491-0.00009;
 double rhoWMax = 0.1491+0.00009;

 TCanvas* c = new TCanvas("c","c",800,600);
 gPad->SetLeftMargin(0.14);
 gPad->SetBottomMargin(0.14);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.05);
 gPad->SetTicks();
 gStyle->SetOptStat(0);

 TLegend* leg = new TLegend(0.45,0.4,0.95,0.95);
 leg->SetTextFont(43);
 leg->SetTextSize(26);
 leg->SetLineWidth(0.0);
 leg->SetFillStyle(0);

 for(int i=0;i<nmult;i++){
	for(int j=0;j<npt;j++){
		for(int r=0;r<3;r++){
			hMassSig[r][i][j] = (TH1D*)fin->Get(Form("hProjInvMassSub_%d_%d_%d_%d",LP,r,i,j));
		}
		hMassSig[3][i][j] = (TH1D*)fin->Get(Form("hProjInvMassSub_InJet_%d_%d_%d",LP,i,j));
	}
 }

 double nsigIntgr = 3.;

 for(int r=0;r<nrt;r++){
	for(int i=0;i<nmult;i++){
		for(int j=0;j<npt;j++){
			fitBin[r][i][j] = new TF1("ffit",sfit,0.7,1.8,15);
			for(int k=0;k<4;k++){
				fitBinComp[r][i][j][k] = new TF1("ffit",sfit,0.7,1.8,15);
			}

			fitBin[r][i][j]->FixParameter(13, 0.160 );
			fitBin[r][i][j]->FixParameter(14, (p_min[j]+p_max[j])/2. );

			fitBin[r][i][j]->SetParameter(1, 2.0 );
			fitBin[r][i][j]->SetParameter(2, -4.0 );
			fitBin[r][i][j]->SetParameter(3, 0.0 );

			fitBin[r][i][j]->SetParLimits(0,hMassSig[r][i][j]->GetMaximum()*0.1,BGAmpMax);
			fitBin[r][i][j]->SetParLimits(1,BGSlopeMin,BGSlopeMax);
			fitBin[r][i][j]->SetParLimits(2,BGInd1Min,BGInd1Max);
			fitBin[r][i][j]->SetParLimits(3,BGInd2Min,BGInd2Max);

			fitBin[r][i][j]->SetParLimits(4,f0MassMinRange,f0MassMaxRange);
			fitBin[r][i][j]->SetParLimits(6,0.054,0.056);
			if( LP==0 ) fitBin[r][i][j]->SetParLimits(6,0.01, 0.1);

			fitBin[r][i][j]->SetParLimits(7,f2MassMinRange,f2MassMaxRange);
			fitBin[r][i][j]->SetParLimits(9,f2WMin,f2WMax);

			fitBin[r][i][j]->SetParLimits(10,rhoMassMin,rhoMassMax);
			fitBin[r][i][j]->SetParLimits(12,rhoWMin,rhoWMax);

			fitBin[r][i][j]->FixParameter(7, (f2MassMinRange+f2MassMaxRange)/2. );
			fitBin[r][i][j]->FixParameter(9, (f2WMin+f2WMax)/2. );
			fitBin[r][i][j]->FixParameter(10, (rhoMassMin+rhoMassMax)/2. );
			fitBin[r][i][j]->FixParameter(12, (rhoWMin+rhoWMax)/2. );

			if( j==0 ){
				f0AmpMin = 1.0*2.0*1e1;
				f2AmpMin = 1.0*2.0*1e1;
				rhoAmpMin = 1.0*2.0*1e1;

				rhoAmpMax = hMassSig[r][i][j]->GetMaximum()*1e4;
				f0AmpMax = hMassSig[r][i][j]->GetMaximum()*1e3;
				f2AmpMax = hMassSig[r][i][j]->GetMaximum()*1e2;

			} else if( j>0 ){
			    f0AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			    f0AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			    f2AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			    f2AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			    rhoAmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			    rhoAmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			}

			fitBin[r][i][j]->SetParLimits(5,f0AmpMin, f0AmpMax );
			fitBin[r][i][j]->SetParLimits(8,f2AmpMin, f2AmpMax);
			fitBin[r][i][j]->SetParLimits(11,rhoAmpMin, rhoAmpMax);

			if( r ){
				fitBin[r][i][j]->SetParLimits(1, fitBin[0][i][j]->GetParameter(0)*0.1, fitBin[0][i][j]->GetParameter(0)*10. );
				fitBin[r][i][j]->SetParLimits(1, fitBin[0][i][j]->GetParameter(1)-0.2, fitBin[0][i][j]->GetParameter(1)+0.2 );
				fitBin[r][i][j]->SetParLimits(2, fitBin[0][i][j]->GetParameter(2)-0.2, fitBin[0][i][j]->GetParameter(2)+0.2 );
				fitBin[r][i][j]->SetParLimits(3, fitBin[0][i][j]->GetParameter(3)-0.2, fitBin[0][i][j]->GetParameter(3)+0.2 );
			}

			if( Rebinning>1 ) hMassSig[r][i][j]->Rebin( Rebinning );
			if( j>8 ) hMassSig[r][i][j]->Rebin( 2 );

			fitResults[r][i][j] = (TFitResultPtr)hMassSig[r][i][j]->Fit(fitBin[r][i][j],"sbq0","",0.8,1.55);

			for(int k=1;k<4;k++){
				if( k==7 || k==9 || k==10 || k==12 ) continue;
				fitBin[r][i][j]->FixParameter( k, fitBin[r][i][j]->GetParameter(k) );

			}
			fitResults[r][i][j] = (TFitResultPtr)hMassSig[r][i][j]->Fit(fitBin[r][i][j],"sb","",0.8,1.55);

			if (j<npt-1){
				gf0Width[r][i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[r][i][j]->GetParameter(6) );
				gf0Width[r][i]->SetPointError( j, 0, fitBin[r][i][j]->GetParError(6) );

				gf0Mass[r][i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[r][i][j]->GetParameter(4) );
				gf0Mass[r][i]->SetPointError( j, 0, fitBin[r][i][j]->GetParError(4) );

				gChi2OverNDF[r][i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[r][i][j]->GetChisquare()/fitBin[r][i][j]->GetNDF() );
			}

	        c->cd();

			for(int k=0;k<4;k++){
				for(int p=0;p<15;p++){
					fitBinComp[r][i][j][k]->SetParameter( p, fitBin[r][i][j]->GetParameter(p) );
				}
				for(int p=0;p<4;p++){
					if( k==p ) continue;
					fitBinComp[r][i][j][k]->SetParameter( MagIndex[p], 0 );
				}
			}

	        hMassSig[r][i][j]->GetXaxis()->SetRangeUser(0.8,1.85);
	        hMassSig[r][i][j]->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
	        hMassSig[r][i][j]->GetYaxis()->SetTitle(Form("Counts / %.1lf MeV",hMassSig[r][i][j]->GetBinWidth(1)*1e3));
	        hMassSig[r][i][j]->SetTitle("");

			hMassSig[r][i][j]->GetXaxis()->SetTitleFont(43);
			hMassSig[r][i][j]->GetXaxis()->SetLabelFont(43);
			hMassSig[r][i][j]->GetYaxis()->SetTitleFont(43);
			hMassSig[r][i][j]->GetYaxis()->SetLabelFont(43);
 
			hMassSig[r][i][j]->GetXaxis()->SetTitleSize(32);
			hMassSig[r][i][j]->GetYaxis()->SetTitleSize(32);
			hMassSig[r][i][j]->GetXaxis()->SetLabelSize(28);
			hMassSig[r][i][j]->GetYaxis()->SetLabelSize(28);

			hMassSig[r][i][j]->SetMarkerStyle(24);
			hMassSig[r][i][j]->SetMarkerColor(1);
			hMassSig[r][i][j]->SetLineColor(1);
			hMassSig[r][i][j]->GetFunction("ffit")->SetLineColor(kBlue);
			hMassSig[r][i][j]->SetMinimum(0);
			hMassSig[r][i][j]->SetMaximum( hMassSig[r][i][j]->GetBinContent( hMassSig[r][i][j]->GetXaxis()->FindBin(0.8) )*1.1 );
			hMassSig[r][i][j]->Draw();
			for(int k=0;k<4;k++){
			    fitBinComp[r][i][j][k]->SetLineWidth(2.0);
			    fitBinComp[r][i][j][k]->SetLineColor( lcolor[k] );
			    fitBinComp[r][i][j][k]->Draw("same");
			}

			leg->Clear();

			leg->AddEntry( (TObject*)0, "ALICE work in progress", "");
			leg->AddEntry( (TObject*)0, "LHC22m apass4", "");
			leg->AddEntry( (TObject*)0, Form(" FT0M %.0lf#font[122]{-}%.0lf%%",m_min[i],m_max[i]), "");
			leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.5", "");
			leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
			leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",p_min[j],p_max[j]), "");
			leg->AddEntry( fitBinComp[r][i][j][0], "Residual bkg.", "l");
			leg->AddEntry( fitBinComp[r][i][j][1], "f_{0}(980)", "l");
			leg->AddEntry( fitBinComp[r][i][j][3], "#rho(770)^{0}", "l");
			leg->AddEntry( fitBinComp[r][i][j][2], "f_{2}(1270)", "l");

			leg->Draw();

			c->SaveAs(Form("figs/Fit_%dLP/fitout_%d_%d_%d.pdf",LP,r,i,j));

			fitBin[r][i][j]->SetParameter(0, 0);
			fitBin[r][i][j]->SetParameter(8, 0);
			fitBin[r][i][j]->SetParameter(11, 0);

			Intgr_RawY_cntl[r][i][j] = fitBin[r][i][j]->Integral(
				fitBin[r][i][j]->GetParameter(4) - fitBin[r][i][j]->GetParameter(6)*nsigIntgr,
				fitBin[r][i][j]->GetParameter(4) + fitBin[r][i][j]->GetParameter(6)*nsigIntgr );
			Intgr_RawY_cntl[r][i][j] /= hMassSig[r][i][j]->GetBinWidth(1);

			Intgr_RawY_stat[r][i][j] = fitBin[r][i][j]->IntegralError(
				fitBin[r][i][j]->GetParameter(4) - fitBin[r][i][j]->GetParameter(6)*nsigIntgr,
                fitBin[r][i][j]->GetParameter(4) + fitBin[r][i][j]->GetParameter(6)*nsigIntgr,
				fitBin[r][i][j]->GetParameters(), fitResults[r][i][j]->GetCovarianceMatrix().GetMatrixArray(), 1e-6 );
			Intgr_RawY_stat[r][i][j] /= hMassSig[r][i][j]->GetBinWidth(1);
		}
	}
 }


 for(int r=0;r<3;r++){
    for(int i=0;i<nmult;i++){
		hPtStat[r][i] = new TH1D(Form("hPtStat_%d_%d",r,i),"",npt-1,PtBinnings);
        for(int j=0;j<npt-1;j++){
			hPtStat[r][i]->SetBinContent( j+1, Intgr_RawY_cntl[r][i][j] );
			hPtStat[r][i]->SetBinError( j+1, Intgr_RawY_stat[r][i][j] );
		}
	}
 }

 TFile* fout = new TFile(Form("F0IntgrOut_%dLP.root",LP),"recreate");
 for(int r=0;r<3;r++){
    for(int i=0;i<nmult;i++){
		hPtStat[r][i]->Write();
	}
 }

}

void Fit(){
 FitLP(0, 1);
 FitLP(1, 2);
}
/*

 TH1D* hMassSig[nmult][npt];
 TF1* fitBin[nmult][npt];
 TF1* fitBinComp[nmult][npt][4];
 int MagIndex[4] = {0,5,8,11};
 int lcolor[4] = {46, 30, 41, 28};
 double BGAmpMin = 1e2*2.0*1e2;
 double BGAmpMax = 3e7*2.0*1e2;

 double BGSlopeMin = -2.5;
 double BGSlopeMax = 3.0;

 double BGInd1Min = -5.0;
 double BGInd1Max = 5.0;

 double BGInd2Min = -0.5;
 double BGInd2Max = 5.0;

 double f0MassMinRange = 0.95;
 double f0MassMaxRange = 1.0;

 double f0AmpMin = 1.0*2.0;
 double f0AmpMax = 1e6*2.0;

 double f0WMin = 0.01;
 double f0WMax = 0.10;

 double f2MassMinRange = 1.2755-0.0012*0.002;
 double f2MassMaxRange = 1.2755+0.0012*0.002;

 double f2AmpMin = 1.0*2.0;
 double f2AmpMax = 1e6*2.0;

 double f2WMin = 0.1867-0.0024*0.02;
 double f2WMax = 0.1867+0.0029*0.02;

 double rhoAmpMin = 1.0*2.0;
 double rhoAmpMax = 1e6*2.0;

 double rhoMassMin = 0.7752;
 double rhoMassMax = 0.7754;

 double rhoWMin = 0.1491-0.00009;
 double rhoWMax = 0.1491+0.00009;

 TCanvas* c = new TCanvas("c","c",800,600);
 gPad->SetLeftMargin(0.14);
 gPad->SetBottomMargin(0.14);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();
 gStyle->SetOptStat(0);

 TLegend* leg = new TLegend(0.45,0.4,0.95,0.95);
 leg->SetTextFont(43);
 leg->SetTextSize(26);
 leg->SetLineWidth(0.0);
 leg->SetFillStyle(0);

 for(int i=0;i<nmult;i++){
	for(int j=0;j<npt;j++){
		hMassSig[i][j] = (TH1D*)fin->Get(Form("hInvMass_sub_%d_%d",i,j));
		hMassSig[i][j]->Rebin(2);

		fitBin[i][j] = new TF1("f1",sfit,0.7,1.8,15);
		for(int k=0;k<4;k++){
			fitBinComp[i][j][k] = new TF1("f1",sfit,0.7,1.8,15);
		}

		fitBin[i][j]->FixParameter(13, 0.160 );
		fitBin[i][j]->FixParameter(14, (p_min[j]+p_max[j])/2. );

//		fitBin[i][j]->FixParameter( 6, 0.055 );
		fitBin[i][j]->SetParameter( 4, 0.98 );

        fitBin[i][j]->SetParameter(1, 2.0 );
        fitBin[i][j]->SetParameter(2, -4.0 );
        fitBin[i][j]->SetParameter(3, 0.0 );

        fitBin[i][j]->SetParLimits(0,hMassSig[i][j]->GetMaximum()*0.1,BGAmpMax);
        fitBin[i][j]->SetParLimits(1,BGSlopeMin,BGSlopeMax);
        fitBin[i][j]->SetParLimits(2,BGInd1Min,BGInd1Max);
        fitBin[i][j]->SetParLimits(3,BGInd2Min,BGInd2Max);

        fitBin[i][j]->SetParLimits(4,f0MassMinRange,f0MassMaxRange);
        fitBin[i][j]->SetParLimits(6,0.04,0.06);

        fitBin[i][j]->SetParLimits(7,f2MassMinRange,f2MassMaxRange);
        fitBin[i][j]->SetParLimits(9,f2WMin,f2WMax);

        fitBin[i][j]->SetParLimits(10,rhoMassMin,rhoMassMax);
        fitBin[i][j]->SetParLimits(12,rhoWMin,rhoWMax);

		if( j==0 ){
			f0AmpMin = 1.0*2.0*1e2;
			f0AmpMax = 1e6*2.0*1e2;
			f2AmpMin = 1.0*2.0*1e2;
			f2AmpMax = 1e6*2.0*1e2;
			rhoAmpMin = 1.0*2.0*1e2;
			rhoAmpMax = 1e6*2.0*1e2;
		} else if( j>0 ){
            f0AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
            f0AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
            f2AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
            f2AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
            rhoAmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
            rhoAmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
        }

		fitBin[i][j]->SetParLimits(5,f0AmpMin, f0AmpMax );
		fitBin[i][j]->SetParLimits(8,f2AmpMin, f2AmpMax);
		fitBin[i][j]->SetParLimits(11,rhoAmpMin, rhoAmpMax);

		hMassSig[i][j]->Fit(fitBin[i][j], "sb","",0.8,1.76);

		for(int k=0;k<4;k++){
			for(int p=0;p<15;p++){
				fitBinComp[i][j][k]->SetParameter( p, fitBin[i][j]->GetParameter(p) );
			}
			for(int p=0;p<4;p++){
				if( k==p ) continue;
				fitBinComp[i][j][k]->SetParameter( MagIndex[p], 0 );
			}
		}

		c->cd();

		hMassSig[i][j]->GetXaxis()->SetRangeUser(0.8,1.85);
		hMassSig[i][j]->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
		hMassSig[i][j]->GetYaxis()->SetTitle(Form("Counts / %.1lf MeV",hMassSig[i][j]->GetBinWidth(1)*1e3));
		hMassSig[i][j]->SetTitle("");

		hMassSig[i][j]->GetXaxis()->SetTitleFont(43);
		hMassSig[i][j]->GetXaxis()->SetLabelFont(43);
		hMassSig[i][j]->GetYaxis()->SetTitleFont(43);
		hMassSig[i][j]->GetYaxis()->SetLabelFont(43);
 
		hMassSig[i][j]->GetXaxis()->SetTitleSize(32);
		hMassSig[i][j]->GetYaxis()->SetTitleSize(32);
		hMassSig[i][j]->GetXaxis()->SetLabelSize(28);
		hMassSig[i][j]->GetYaxis()->SetLabelSize(28);

		hMassSig[i][j]->SetMarkerStyle(24);
		hMassSig[i][j]->SetMarkerColor(1);
		hMassSig[i][j]->SetLineColor(1);
		hMassSig[i][j]->GetFunction("f1")->SetLineColor(kBlue);
		hMassSig[i][j]->SetMinimum(0);
		hMassSig[i][j]->Draw();
		for(int k=0;k<4;k++){
			fitBinComp[i][j][k]->SetLineWidth(2.0);
		    fitBinComp[i][j][k]->SetLineColor( lcolor[k] );
		    fitBinComp[i][j][k]->Draw("same");
		}

		leg->Clear();

		leg->AddEntry( (TObject*)0, "ALICE work in progress", "");
		leg->AddEntry( (TObject*)0, "LHC22o apass4", "");
		leg->AddEntry( (TObject*)0, Form(" FT0M %.0lf#font[122]{-}%.0lf%%",m_min[i],m_max[i]), "");
		leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.5", "");
		leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",p_min[j],p_max[j]), "");
		leg->AddEntry( fitBinComp[i][j][0], "Residual bkg.", "l");
		leg->AddEntry( fitBinComp[i][j][1], "f_{0}(980)", "l");
		leg->AddEntry( fitBinComp[i][j][3], "#rho(770)^{0}", "l");
		leg->AddEntry( fitBinComp[i][j][2], "f_{2}(1270)", "l");


		leg->Draw();

		c->SaveAs(Form("figs/Fit/fitout_%d_%d.pdf",i,j));

	}
 }


 hMassSig[0][1]->GetXaxis()->SetRangeUser(0.5,1.85);
 hMassSig[0][1]->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
 hMassSig[0][1]->GetYaxis()->SetTitle(Form("Counts / %.1lf MeV",hMassSig[0][1]->GetBinWidth(1)*1e3));
 hMassSig[0][1]->SetTitle("");

 hMassSig[0][1]->GetXaxis()->SetTitleFont(43);
 hMassSig[0][1]->GetXaxis()->SetLabelFont(43);
 hMassSig[0][1]->GetYaxis()->SetTitleFont(43);
 hMassSig[0][1]->GetYaxis()->SetLabelFont(43);
 
 hMassSig[0][1]->GetXaxis()->SetTitleSize(32);
 hMassSig[0][1]->GetYaxis()->SetTitleSize(32);
 hMassSig[0][1]->GetXaxis()->SetLabelSize(28);
 hMassSig[0][1]->GetYaxis()->SetLabelSize(28);

 hMassSig[0][1]->SetMarkerStyle(24);
 hMassSig[0][1]->SetMarkerColor(1);
 hMassSig[0][1]->SetLineColor(1);
 hMassSig[0][1]->GetFunction("f1")->SetLineColor(kBlue);
 hMassSig[0][1]->Draw();
 for(int k=0;k<4;k++){
	fitBinComp[0][1][k]->SetLineWidth(2.0);
	fitBinComp[0][1][k]->SetLineColor( lcolor[k] );
	fitBinComp[0][1][k]->Draw("same");
 }

 TLegend* leg = new TLegend(0.4,0.4,0.95,0.95);
 leg->SetTextFont(43);
 leg->SetTextSize(26);
 leg->SetLineWidth(0.0);
 leg->SetFillStyle(0);


 leg->AddEntry( (TObject*)0, "ALICE work in progress", "");
 leg->AddEntry( (TObject*)0, "LHC22e pass3, FT0M 0#font[122]{-}80%", "");
 leg->AddEntry( (TObject*)0, "pp 900 GeV , |#it{y}|<0.5", "");
 leg->AddEntry( (TObject*)0, "1.5 < #it{p}_{T} < 6.0 GeV/#it{c}", "");
 leg->AddEntry( fitBinComp[0][1][0], "Residual bkg.", "l");
 leg->AddEntry( fitBinComp[0][1][1], "f_{0}(980)", "l");
 leg->AddEntry( fitBinComp[0][1][3], "#rho(770)^{0}", "l");
 leg->AddEntry( fitBinComp[0][1][2], "f_{2}(1270)", "l");


 leg->Draw();

 c->SaveAs("fitres.pdf");
*/


