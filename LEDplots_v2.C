//Second attempt at plotting ave int charge vs HV setting for LED pulse data - sseeds 12.10.20

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const int kNrows = 24; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int kNLED = 6; //Number of led channels = 5 leds + off = 6 channels
const int runT = 13; //Total number of runs
double chiSqr = 0.0;

TF1 *t1;
TF1 *t2;
TF1 *t3;

TFile *LEDpulse_HVcal = new TFile("LEDpulse_HVCal.root","RECREATE");
//TGraph *G1 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt for int val
//TGraph *G2 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt for int val

//TF1 *pFit = new TF1("pFit","[0]*pow((x-[1]),[2])+[3]"); //Produces remarkably bad fits
//TF1 *pFit = new TF1("pFit","[0]*pow((x-[1])/1200,[2])",0,2000); //Normalize x
//TF1 *pFit = new TF1("pFit","pow(x,[0])",0,2000);
//TF1 *pFit;

double intMean[kNrows][kNcols][kNLED][runT];
double maxMean[kNrows][kNcols][kNLED][runT];
double intSD[kNrows][kNcols][kNLED][runT];
double maxSD[kNrows][kNcols][kNLED][runT];
double HVset[kNrows][kNcols][runT];
double NPE[kNrows][kNcols][kNLED][runT]; //Just look at LED 1 for now using max ADC

//TH1F *intAlpha = new TH1F("intAlpha","intAlpha",50,0,0);
//TH1F *maxAlpha = new TH1F("intAlpha","intAlpha",50,0,0);

double powFit(double *X,double *p) 
{
  //Double_t fitval = par[0]*pow(X[0],par[1]);
  double fitval = p[0]*pow((X[0]-p[1])/1200,p[2]) + p[3]; //X normalized to lowest HV value
  return fitval;
}

void LEDplots_v2(int run1 = 1286, int run2 = 1288, int run3 = 1289, int run4 = 1290, int run5 = 1291, int run6 = 1292, int run7 = 1293, int run8 = 1294, int run9 = 1297, int run10 = 1300, int run11 = 1301, int run12 = 1302, int run13 = 1303){ 

  cout << "Beginning loop over runs " << run1 << ", " << run2 << ", " << run3 << ", " << run4 << ", " << run5 << ", " << run6 << ", " << run7 << ", " << run8 << ", " << run9 << ", " << run10 << ", " << run11 << ", " << run12 << ", and " << run13 << ".." << endl;  
  int runs[runT]={run1,run2,run3,run4,run5,run6,run7,run8,run9,run10,run11,run12,run13};

  cout << "HV settings and LED intensities hardcoded from HCAL wiki for now.." << endl;

  int runHV_jlab[runT] = {1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050};
  int runHV_cmu[runT] = {1350,1400,1425,1450,1475,1500,1525,1550,1575,1600,1625,1650,1675};
  
  double ledInten[kNLED]= {0,1,2,4,8,16};  //Hardcoded again as led intensity relative to led 1.
  
  //Loop over all runs and fill master array
  cout << "Filling master array from analysis histograms.." << endl;
  for(int i=0; i<runT; i++){

    TFile *f = TFile::Open(Form("outFiles/ledSpectFile_run%d.root",runs[i]));  //Read in file
    
    for(int r=12; r<kNrows; r++){
      for(int c=0; c<kNcols; c++){

	TH1F *h3 = (TH1F*)f->Get(Form("Pedestal Spect R%d C%d",r,c));
	
	for(int l=1; l<kNLED; l++){ //Ignoring led = 0, as no led is on for this setting
	  	  
	  TH1F *h1 = (TH1F*)f->Get(Form("Int ADC Spect R%d C%d LED%d",r,c,l));  //Save all integrated spectra to master histo
	  TH1F *h2 = (TH1F*)f->Get(Form("Max ADC Spect R%d C%d LED%d",r,c,l));  //Save all integrated spectra to master histo
	  
	  if(h1){
	    h1->Fit("gaus","Q");
	    t1=h1->GetFunction("gaus");
	    intMean[r][c][l][i] = t1->GetParameter(1);

	    intSD[r][c][l][i] = t1->GetParameter(2);
	  }

	  if(h2){
	    h2->Fit("gaus","Q");
	    t2=h2->GetFunction("gaus");
	    maxMean[r][c][l][i] = t2->GetParameter(1);

	    maxSD[r][c][l][i] = t2->GetParameter(2);
	    if(h3){
	      h3->Fit("gaus","Q");
	      t3=h3->GetFunction("gaus");
	      NPE[r][c][l][i]=pow(t2->GetParameter(1)/pow(pow(t2->GetParameter(2),2.0)-pow(t3->GetParameter(2),2.0),0.5),2.0);
	    }
	  }	  
	}	
      }
    }

    f->Close();

    cout << "Ending run " << runs[i] << ".." << endl;
  }

  cout << "Creating file to hold fit parameters.." << endl;
  ofstream outFile;
  outFile.open("PMTCalibrationData.txt");
  outFile << "#Alpha parameter obtained from integrated pulse data." << endl;
  outFile << "Data taken from runs " << run1 << ", " << run2 << ", " << run3 << ", " << run4 << ", " << run5 << ", " << run6 << ", and " << run7 << "." << endl;  
  outFile << "#Fit to [0]*pow(x-[1],[2])+[3]. Data represents average over all LED settings." << endl;
  outFile << "#Row  Column <intAlpha> intsigma <maxAlpha> maxsigma" << endl;
  
  //Now use intMean to draw TGraphs
  cout << "Using master array to draw graphs and obtaining fit parameters.." << endl;
  for(int r=12; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      TH1F *intAlpha = new TH1F(Form("R%dC%dIntAlpha",r,c),Form("R%dC%dIntAlpha",r,c),50,0,0);
      for(int l=1; l<kNLED; l++){ //LED zero is where all LEDs are off
	double HV[runT];
	double Sum[runT];
	double SD[runT];
	for(int i=0; i<runT; i++){
	  if(c>3&&c<8){
	    HV[i]=runHV_jlab[i];
	    HVset[r][c][i]=runHV_jlab[i];
	  }else{
	    HV[i]=runHV_cmu[i];
	    HVset[r][c][i]=runHV_cmu[i];  
	  }
	  Sum[i]=intMean[r][c][l][i];
	  SD[i]=intSD[r][c][l][i]/sqrt(5000);  //Error is (std dev)/sqrt(N). Hardcoded for now; can use h1->GetEntries() to obtain N
	}
	TGraph *G1 = new TGraphErrors(runT,HV,Sum,0,SD);
	TF1 *pFit = new TF1("pFit",powFit,1100,1500,4);

	//pFit->SetParLimits(0,0,0.01);
	if(c>3&&c<8){
	  pFit->SetParameters(0,10,8,0);
	  pFit->SetParLimits(1,0,1500);
	  pFit->SetParLimits(2,1,100);
	}else{
	  //pFit->SetParameters(0,1200,12,0);
	  pFit->SetParameters(0,800,10,10);
	  //pFit->SetParLimits(1,0,1500);
	  //pFit->SetParLimits(2,1,100);
	}
	
	//G1->Fit("pol2","QR"); //options M to improve minimum and fit
	G1->Fit("pFit","Q");
	//G1->SetMinimum(0.0);
	G1->SetTitle(Form("SUM R%d-C%d LED%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	//G1->SetTitle(Form("SUM R%d-C%d LED%d: coeff=%g shift=%f alpha=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2)));
	//G1->SetTitle(Form("%d-%d LED%d: alpha=%f",r,c,l,pFit->GetParameter(0)));
	G1->GetYaxis()->SetTitle("Average Integrated Signal (sRAU)");  //1 RAU approx 0.5 mV
	G1->GetYaxis()->CenterTitle();
	G1->GetXaxis()->SetTitle("PMT HV Setting (-V)");
	G1->GetXaxis()->CenterTitle();
	//G1->GetXaxis()->SetLimits(0.0,1500.0);
	G1->SetMinimum(0.0);
	G1->SetLineColor(kWhite);
	G1->SetMarkerStyle(8);
	G1->SetMarkerSize(1);
	//chiSqr = pFit->GetChisquare();
	LEDpulse_HVcal->cd();
	G1->Write(Form("SUM R%d-C%d LED%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	//G1->Write(Form("SUM R%d-C%d LED%d: coeff=%g shift=%f alpha=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2)));
	//G1->Write(Form("%d-%d LED%d: alpha=%f",r,c,l,pFit->GetParameter(0)));
	intAlpha->Fill(pFit->GetParameter(2));
	G1->Draw("AP");
      }
      TH1F *maxAlpha = new TH1F(Form("R%dC%dMaxAlpha",r,c),Form("R%dC%dMaxAlpha",r,c),50,0,0);
      for(int l=1; l<kNLED; l++){
	double HV[runT];
	double Max[runT];
	double SD[runT];
	for(int i=0; i<runT; i++){
	  if(c>3&&c<8){
	    HV[i]=runHV_jlab[i];
	  }else{
	    HV[i]=runHV_cmu[i];
	  }	  
	  Max[i]=maxMean[r][c][l][i];
	  SD[i]=maxSD[r][c][l][i]/sqrt(5000);  //Error is (std dev)/sqrt(N)
	}
	TGraph *G1 = new TGraphErrors(runT,HV,Max,0,SD);	
	TF1 *pFit = new TF1("pFit",powFit,1100,1500,4);

	//pFit->SetParLimits(0,0,0.01);
	if(c>3&&c<8){
	  pFit->SetParameters(0,10,8,0);
	  pFit->SetParLimits(1,0,1500);
	  pFit->SetParLimits(2,1,100);
	}else{
	  pFit->SetParameters(0,800,10,10);
	  //pFit->SetParLimits(1,0,1500);
	  //pFit->SetParLimits(2,1,100);
	}
	
	//G1->Fit("pol2","QR"); //options M to improve minimum and fit
	G1->Fit("pFit","Q");
	//G1->SetMinimum(0.0);
	G1->SetTitle(Form("MAX R%d-C%d LED%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	//G1->SetTitle(Form("MAX R%d-C%d LED%d: coeff=%g shift=%f alpha=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2)));
	//G1->SetTitle(Form("%d-%d LED%d: alpha=%f",r,c,l,pFit->GetParameter(0)));
	G1->GetYaxis()->SetTitle("Average fADC Amplitude (RAU)");  //1 RAU approx 0.5 mV
	G1->GetYaxis()->CenterTitle();
	G1->GetXaxis()->SetTitle("PMT HV Setting (-V)");
	G1->GetXaxis()->CenterTitle();
	//G1->GetXaxis()->SetLimits(0.0,1500.0);
	G1->SetMinimum(0.0);
	G1->SetLineColor(kWhite);
	G1->SetMarkerStyle(8);
	G1->SetMarkerSize(1);
	//chiSqr = pFit->GetChisquare();
	LEDpulse_HVcal->cd();
	G1->Write(Form("MAX R%d-C%d LED%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	//G1->Write(Form("MAX R%d-C%d LED%d: coeff=%g shift=%f alpha=%f",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2)));
	//G1->Write(Form("%d-%d LED%d: alpha=%f",r,c,l,pFit->GetParameter(0)));
	maxAlpha->Fill(pFit->GetParameter(2));
	G1->Draw("AP");
      }
      /*
      for(int hv=0; hv<runT; hv++){
	double LI[kNLED];
	double SumI[kNLED];
	for(int i=0; i<kNLED; i++){
	  LI[i]=ledInten[i];
	  SumI[i]=intMean[r][c][i][hv];
	}
	TGraph *G2 = new TGraph(kNLED,LI,SumI);
	//G2->Fit("pFit","QR");
	//G2->SetMinimum(0.0);
	G2->SetTitle(Form("R%d-C%d HV%d",r,c,hv));
	G2->GetYaxis()->SetTitle("Average Integrated Signal ");  //Need actual units here
	G2->GetXaxis()->SetTitle("LED Intensity Setting Relative to LED 1");
	G2->SetLineColor(kWhite);
	G2->SetMarkerStyle(8);
	G2->SetMarkerSize(1);
	LEDpulse_HVcal->cd();
	G2->Write(Form("R%d-C%d HV%d",r,c,hv));
	G2->Draw("AP");
      }
      */
      //outFile << r << "  " << c << "  " << intAlpha->GetMean() << "  " << intAlpha->GetRMS() << "  " << maxAlpha->GetMean() << "  " << maxAlpha->GetRMS() << "  " << pow(maxAlpha->GetMean()/maxAlpha->GetRMS(),2) << endl;
      outFile << r << "  " << c << "  " << intAlpha->GetMean() << "  " << intAlpha->GetRMS() << "  " << maxAlpha->GetMean() << "  " << maxAlpha->GetRMS() << endl;
      
    }
  }

  // Draw NPE v HV graphs to look for plateau region
  for(int r=12; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=1; l<kNLED; l++){
	TGraph *G3 = new TGraph(runT, HVset[r][c], NPE[r][c][l]);
	G3->SetTitle(Form("R%d-C%d LED %d NPE v HV",r,c,l));
	G3->GetYaxis()->SetTitle("NPE");
	G3->GetXaxis()->SetTitle("HV Setting");
	G3->SetMinimum(0.0);
	G3->SetLineColor(kWhite);
	G3->SetMarkerStyle(8);
	G3->SetMarkerSize(1);
	LEDpulse_HVcal->cd();
	G3->Write(Form("R%d-C%d LED %d NPE v HV",r,c,l));
	G3->Draw("AP");
      }
    }
  }
  //%g, %e
  
  cout << "Completed loop over runs. Histograms written to LEDpulse_HVcal.root. Fit parameters written to file PMTCalibrationData.txt." << endl;
}
