//Third attempt at plotting ave int charge vs HV setting - sseeds 12.10.20

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const int kNrows = 12; //Number of pmt rows
const int kNcols = 12; //Number of pmt columns
const int kNLED = 6; //Number of led channels = 5 leds + off = 6 channels
//const int HVt = 5;
//const int PMTt = 144;
const int runT = 7; //Seven total runs read in by function

TFile *LEDpulse_HVcal = new TFile("LEDpulse_HVCal.root","RECREATE");
TGraph *G1 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt for int val
//TGraph *G2 = new TGraph(); //For max val

//TCanvas *C1 = new TCanvas("C1","PMT intPulse",900,700); //First 12 x 6 grid
//TCanvas *C2 = new TCanvas("C2","PMT intPulse",900,700); //Temporary for max val graphs

TF1 *pFit = new TF1("pFit","[0]*pow(x,[1])+[2]",0,2000); //Produces remarkably bad fits
//TF1 *testFit;
//TF1 *pFitted;
TH1F *HI[kNrows][kNcols][kNLED][runT];
TH1F *HP[kNrows][kNcols][kNLED][runT];
double intMean[kNrows][kNcols][kNLED][runT];


void LEDplots(int run1 = 1198, int run2 = 1200, int run3 = 1201, int run4 = 1202, int run5 = 1203, int run6 = 1204, int run7 = 1205){ 

  int runs[runT]={run1,run2,run3,run4,run5,run6,run7};

  int runHV[runT]={1200,1250,1300,1350,1400,1450,1425};  //Will hardcode these corresponding to defaults for now. Will be better to read in conversion file from text. Later.

  //Loop over all runs and fill master array
  for(int i=0; i<runT; i++){
    for(int r=0; r<kNrows; r++){
      for(int c=0; c<kNcols; c++){
	for(int l=0; l<kNLED; l++){
	  
	  TFile *f = TFile::Open(Form("ledIntFile_run%d.root",runs[i]));  //Read in file
	  f->ls();
	  
	  TH1F *h1 = (TH1F*)f->Get(Form("Int ADC Spect R%d C%d LED%d",r,c,l));  //Save all integrated spectra to master histo

	  intMean[r][c][l][i] = h1->GetMean(); //Get mean from the integrated specra and save to master array

	  //HI[r][c][l][i]=h1;

	  //TH1F *h2 = (TH1F*)f->Get("Pedestal Spect R%d C%d LED%d",r,c,l); //Save all pedestals to master histo
	  //HP[r][c][l][i]=h1;
     
	}
      }
    }
  }

  //Now use intMean to draw TGraphs
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=0; l<kNLED; l++){
	double HV[runT];
	double Sum[runT];
	for(int i=0; i<runT; i++){
	  HV[i]=runHV[i];
	  Sum[i]=intMean[r][c][l][i];
	}
	TGraph *G1 = new TGraph(runT,HV,Sum);
	G1->Fit("pFit","QR");
	G1->SetTitle(Form("%d-%d LED%d: coeff%d alpha%d offset%d",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(0)));
	G1->GetYaxis()->SetTitle("Average Integrated Signal");  //Need actual units here
	G1->GetXaxis()->SetTitle("PMT HV Setting");
	LEDpulse_HVcal->cd();
	G1->Write(Form("%d-%d LED%d: coeff%d alpha%d offset%d",r,c,l,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(0)));
	
      }
    }
  }
  
  cout << "Completed loop over runs." << endl;


  
}
