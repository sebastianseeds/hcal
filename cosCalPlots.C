//First attempt at plotting ave int charge vs HV setting for LED pulse data - sseeds 12.10.20

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
const int runT = 6; //Six total runs read in by function
double chiSqr = 0.0;

TFile *F1 = new TFile("cosmicPresPlots.root","RECREATE");

double pmtHV[kNrows][kNcols][runT];
double avgADC[kNrows][kNcols][runT];
double avgTDC[kNrows][kNcols][runT];

TH1F *ChannelADC[runT];
TH1F *ChannelTDC[runT];

TH1F *ADCvChannel[runT];
TH1F *TDCvChannel[runT];

TH1F *ChannelPed[runT];

TGraph *ADCvHV[kNrows][kNcols];
TGraph *TDCvHV[kNrows][kNcols];

double powFit(double *X,double *p) 
{
  //Double_t fitval = par[0]*pow(X[0],par[1]);
  double fitval = p[0]*pow((X[0]-p[1])/1200,(1/p[2])) + p[3]; //X normalized to lowest HV value
  return fitval;
}

void cosCalPlots(int run1 = 1235, int run2 = 1239, int run3 = 1263, int run4 = 1265, int run5 = 1269, int run6 = 1271){ 
  
  cout << "Beginning loop over runs " << run1 << ", " << run2 << ", " << run3 << ", " << run4 << ", " << run5 << ", and " << run6 << ".." << endl;  
  int runs[runT]={run1,run2,run3,run4,run5,run6};
  
  cout << "Reading in HV settings per pmt for each run.." << endl;

  for(int i=1; i<runT; i++){
    
    ifstream file(Form("HV_run%d.txt",runs[i]));

    cout << "Getting HV settings for each pmt for run " << runs[i] << "." << endl;
    
    int n1=0;
    double d1;
    
    int rval, cval;
    string line;
    
    while (getline(file,line)){
      if(line.at(0) == '#'){
	continue;
      }
      
      stringstream ss(line);
      ss >> d1;
      
      rval = floor(n1/kNcols);
      cval = n1 % kNcols;
      
      pmtHV[rval][cval][i-1] = -d1; 
      
      cout << "run = " << runs[i] << ", row = " << rval << ", col = " << cval << ", array val = " << pmtHV[rval][cval][i] << endl;
      
      n1++;
    }
  }

  //Loop over runs and populate ADC and TDC histograms
  cout<< "Filling histograms with average ADC and TDC values by run.." << endl;
  for(int i=0; i<runT; i++){
    TFile *f = TFile::Open(Form("cosmicHistograms_run%d_TEST.root",runs[i]));  //Read in file

    TH1F *h0 = (TH1F*)f->Get("pedvChannel");
    
    TH1F *h1 = (TH1F*)f->Get("ADCvChannel_l"); //Get ADC Histogram
    //ADCvChannel[i] = (TH1F*)f->Get("ADCvChannel_l");
    TH1F *h2 = (TH1F*)f->Get("TDCvChannel_l"); //Get TDC Histogram
    //TDCvChannel[i] = (TH1F*)f->Get("TDCvChannel_l");

    ChannelPed[i] = new TH1F(Form("Avg Pedestal, All Channels, Run %d",runs[i]),Form("Avg Pedestal, All Channels, Run %d",runs[i]),50,0,300);
    
    ChannelADC[i] = new TH1F(Form("Avg ADC, All Channels, Run %d",runs[i]),Form("Avg ADC, All Channels, Run %d",runs[i]),50,0,200);
    ChannelTDC[i] = new TH1F(Form("Avg TDC, All Channels, Run %d",runs[i]),Form("Avg TDC, All Channels, Run %d",runs[i]),50,0,200);

    ADCvChannel[i] = new TH1F(Form("ADC vs Channel Run %d",runs[i]),Form("ADC vs Channel Run %d",runs[i]),kNrows*kNcols,0,288);
    TDCvChannel[i] = new TH1F(Form("TDC vs Channel Run %d",runs[i]),Form("TDC vs Channel Run %d",runs[i]),kNrows*kNcols,0,288);

    for(int j=0;j<h1->GetEntries();j++){

      int row = floor(j/kNcols);
      int col = (j%kNcols);

      //cout << "r " << row << " c " << col << " run " << runs[i] << " idx " << j << " content " << h1->GetBinContent(j+1) << endl;

      TH1F *h3 = (TH1F*)f->Get(Form("Max ADC Spect R%d C%d",row,col));
      TH1F *h4 = (TH1F*)f->Get(Form("Max ADC Spect TDC R%d C%d",row,col));
    
      if(row==11&&col==1){
	//cout << "DOING THE THINGS" << endl;

	
	F1->cd();
	h3->Write(Form("A ChannelADC R%d C%d Run%d",row,col,runs[i]));
	h3->Draw("AP");
	F1->cd();
	h4->Write(Form("A ChannelTDC R%d C%d Run%d",row,col,runs[i]));
	h4->Draw("AP");    
      }

      ChannelPed[i]->Fill(h0->GetBinContent(j+1));

      if(h1->GetBinContent(j+1)<250&&h1->GetBinContent(j+1)>20){
	
	ChannelADC[i]->Fill(h1->GetBinContent(j+1));
	ChannelTDC[i]->Fill(h2->GetBinContent(j+1));

	ADCvChannel[i]->SetBinContent(j+1,h1->GetBinContent(j+1));
	TDCvChannel[i]->SetBinContent(j+1,h2->GetBinContent(j+1));

	//cout << ADCvChannel[i]->GetBinContent(j+1) << endl;
	if(i!=0){
	  avgADC[row][col][i]=h1->GetBinContent(j+1);
	  avgTDC[row][col][i]=h2->GetBinContent(j+1);
	}
	
      }
      
      if(h1->GetBinContent(j+1)<20){

	//cout << "Content less than 20" << endl;
	/*
	F1->cd();
	h3->Write(Form("ChannelADC R%d C%d Run%d idx%d",row,col,runs[i],j));
	h3->Draw("AP");
	F1->cd();
	h4->Write(Form("ChannelTDC R%d C%d Run%d idx%d",row,col,runs[i],j));
	h4->Draw("AP");
	*/
	ChannelADC[i]->Fill(h3->GetMean());
	//cout << runs[i] << " ADC Mean " << j << " " << h3->GetMean() << endl;
	ChannelTDC[i]->Fill(h2->GetMean());
	//cout << runs[i] << " TDC Mean " << j << " " << h4->GetMean() << endl;
	  
	ADCvChannel[i]->SetBinContent(j+1,h3->GetMean());
	TDCvChannel[i]->SetBinContent(j+1,h4->GetMean());
	if(i!=0){
	  avgADC[row][col][i-1]=h3->GetMean();
	  avgTDC[row][col][i-1]=h4->GetMean();
	}
      }

      if(h1->GetBinContent(j+1)>250){
	/*
	F1->cd();
	h3->Write(Form("ChannelADC R%d C%d Run%d idx%d",row,col,runs[i],j));
	h3->Draw("AP");
	F1->cd();
	h4->Write(Form("ChannelTDC R%d C%d Run%d idx%d",row,col,runs[i],j));
	h4->Draw("AP");
	*/
	ChannelADC[i]->Fill(h3->GetMaximumBin()*10+5);
	//cout << "ADC Max Bin " << h3->GetMaximumBin()*10+5 << endl;
	ChannelTDC[i]->Fill(h2->GetMaximumBin()*10+5);
	//cout << "TDC Max Bin " << h4->GetMaximumBin()*10+5 << endl;
	
	ADCvChannel[i]->SetBinContent(j+1,h3->GetMaximumBin()*10+5);
	TDCvChannel[i]->SetBinContent(j+1,h4->GetMaximumBin()*10+5);
	if(i!=0){
	  avgADC[row][col][i]=h3->GetMaximumBin()*10+5;
	  avgTDC[row][col][i]=h4->GetMaximumBin()*10+5;
	}
      }
	
      //cout << "r " << row << " c " << col << " run " << i << "." << endl;
    }
  }
  
  //Draw gain curve tgraphs by pmt
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      

      TF1 *pFit = new TF1("pFit",powFit,1100,1500,4); //Define fit function
	
      //Set limits for fit based on pmt type (jlab or cmu)
      if(c>3&&c<8){
	pFit->SetParameters(0,10,8,0);
	pFit->SetParLimits(1,0,1500);
	pFit->SetParLimits(2,1,100);
      }else{
	pFit->SetParameters(0,800,10,10);
      }
      //Draw ADC graph
      ADCvHV[r][c] = new TGraph(runT,pmtHV[r][c],avgADC[r][c]); 
      ADCvHV[r][c]->Fit("pFit","Q");
      //if(pFit->GetParameter(2)<20&&pFit->GetParameter(2)>5){

      //for(int i=1; i<runT; i++){
      //  cout << pmtHV[r][c][i] << "  " << avgADC[r][c][i] << endl;
      //}
	ADCvHV[r][c]->SetTitle(Form("Avg ADC vs HV R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	ADCvHV[r][c]->GetYaxis()->SetTitle("Avg Max ADC (RAU)");  
	ADCvHV[r][c]->GetYaxis()->CenterTitle();
	ADCvHV[r][c]->GetXaxis()->SetTitle("PMT HV Setting (-V)");
	ADCvHV[r][c]->GetXaxis()->CenterTitle();
	ADCvHV[r][c]->SetMinimum(0.0);
	ADCvHV[r][c]->SetLineColor(kWhite);
	ADCvHV[r][c]->SetMarkerStyle(8);
	ADCvHV[r][c]->SetMarkerSize(1);
	F1->cd();
	ADCvHV[r][c]->Write(Form("Avg ADC R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
	ADCvHV[r][c]->Draw("AP");
	//}
      /*
      //Draw TDC graph
      TDCvHV[r][c] = new TGraph(runT,pmtHV[r][c],avgTDC[r][c]); 
      TDCvHV[r][c]->Fit("pFit","Q");
      TDCvHV[r][c]->SetTitle(Form("Avg TDC vs HV R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      TDCvHV[r][c]->GetYaxis()->SetTitle("Avg Max ADC (RAU)"); 
      TDCvHV[r][c]->GetYaxis()->CenterTitle();
      TDCvHV[r][c]->GetXaxis()->SetTitle("PMT HV Setting (-V)");
      TDCvHV[r][c]->GetXaxis()->CenterTitle();
      TDCvHV[r][c]->SetMinimum(0.0);
      TDCvHV[r][c]->SetLineColor(kWhite);
      TDCvHV[r][c]->SetMarkerStyle(8);
      TDCvHV[r][c]->SetMarkerSize(1);
      F1->cd();
      TDCvHV[r][c]->Write(Form("Avg TDC R%d-C%d: coeff=%g shift=%f alpha=%f offset=%f",r,c,pFit->GetParameter(0),pFit->GetParameter(1),pFit->GetParameter(2),pFit->GetParameter(3)));
      TDCvHV[r][c]->Draw("AP");
      */
    }
  }
  
  //Draw ADC and TDC Histograms
  for(int i=0; i<runT; i++){
    ChannelADC[i]->SetTitle(Form("Avg ADC, All Channels, Run %d",runs[i]));
    ChannelADC[i]->GetYaxis()->SetTitle("Counts");  
    ChannelADC[i]->GetYaxis()->CenterTitle();
    ChannelADC[i]->GetXaxis()->SetTitle("Avg Max ADC (RAU)");
    ChannelADC[i]->GetXaxis()->CenterTitle();
    ChannelADC[i]->SetMinimum(0.0);
    ChannelADC[i]->SetMaximum(150.0);

    F1->cd();
    ChannelADC[i]->Write(Form("ChannelADC Run%d",runs[i]));
    ChannelADC[i]->Draw("AP");

    ChannelPed[i]->SetTitle(Form("Avg Pedestal, All Channels, Run %d",runs[i]));
    ChannelPed[i]->GetYaxis()->SetTitle("Counts");  
    ChannelPed[i]->GetYaxis()->CenterTitle();
    ChannelPed[i]->GetXaxis()->SetTitle("Avg Pedestal (RAU)");
    ChannelPed[i]->GetXaxis()->CenterTitle();
    ChannelPed[i]->SetMinimum(0.0);
    ChannelPed[i]->SetMaximum(150.0);

    F1->cd();
    ChannelPed[i]->Write(Form("ChannelPed Run%d",runs[i]));
    ChannelPed[i]->Draw("AP");

    ChannelTDC[i]->SetTitle(Form("Avg TDC, All Channels, Run %d",runs[i]));
    ChannelTDC[i]->GetYaxis()->SetTitle("Counts");  
    ChannelTDC[i]->GetYaxis()->CenterTitle();
    ChannelTDC[i]->GetXaxis()->SetTitle("Avg Max TDC (RAU)");
    ChannelTDC[i]->GetXaxis()->CenterTitle();
    ChannelTDC[i]->SetMinimum(0.0);
    ChannelTDC[i]->SetMaximum(150.0);
    
    F1->cd();
    ChannelTDC[i]->Write(Form("ChannelTDC Run%d",runs[i]));
    ChannelTDC[i]->Draw("AP");
    
    ADCvChannel[i]->SetTitle(Form("Avg ADC v Channels, Run %d",runs[i]));
    ADCvChannel[i]->SetMaximum(350.0);
    F1->cd();
    ADCvChannel[i]->Write(Form("ADCvChannel Run%d",runs[i]));
    ADCvChannel[i]->Draw("AP");
    
    TDCvChannel[i]->SetTitle(Form("Avg TDC v Channels, Run %d",runs[i]));
    TDCvChannel[i]->SetMaximum(350.0);
    F1->cd();
    TDCvChannel[i]->Write(Form("TDCvChannel Run%d",runs[i]));
    TDCvChannel[i]->Draw("AP");

  }

  
  cout << "Completed loop over runs. Histograms written to cosmicPresPlots.root. Fit parameters written to file cosPMTCalibrationData.txt." << endl;
}



//Scratch code for salvage
/*
//Now use intMean to draw TGraphs
  cout << "Using master array to draw graphs and obtaining fit parameters.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      TH1F *intAlpha = new TH1F(Form("R%dC%dIntAlpha",r,c),Form("R%dC%dIntAlpha",r,c),50,0,0);
      for(int l=1; l<kNLED; l++){ //LED zero is off
	double HV[runT];
	double Sum[runT];
	double SD[runT];
	for(int i=0; i<runT; i++){
	  HV[i]=runHV[i];
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
	G1->GetYaxis()->SetTitle("Average Integrated Signal (sRAU)");  //Need actual units here
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
	  HV[i]=runHV[i];
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
	G1->GetYaxis()->SetTitle("Average fADC Amplitude (RAU)");  //Need actual units here
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
      //outFile << r << "  " << c << "  " << intAlpha->GetMean() << "  " << intAlpha->GetRMS() << "  " << maxAlpha->GetMean() << "  " << maxAlpha->GetRMS() << "  " << pow(maxAlpha->GetMean()/maxAlpha->GetRMS(),2) << endl;
      outFile << r << "  " << c << "  " << intAlpha->GetMean() << "  " << intAlpha->GetRMS() << "  " << maxAlpha->GetMean() << "  " << maxAlpha->GetRMS() << endl;
      
    }
  }
  */
  //%g, %e
