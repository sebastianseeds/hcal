#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 12;
const Int_t kNcols = 12;

const Int_t DISP_MIN_SAMPLE = 85*0;
const Int_t DISP_MAX_SAMPLE = 135*0+200*0+150*0+20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

const Int_t minADC = 0;
const Int_t maxADC = 4000;

Int_t gCurrentEntry = -1;

TChain *T = 0;
TCanvas *C1 = new TCanvas();

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, int, int);
void getPedestal(Int_t);
void integratePulse(TH1F*, TF1*, int, int);

//Augmenting with parameters necessary to obtain integrated landau fit values from good pulses
Double_t nhit = 0;
TH1F *histos[kNrows][kNcols];
vector<double> intValues;
vector<int> intValues_row;
vector<int> intValues_col;
int signalTotal = 1;
TFile *HistosFile = new TFile("HistosFile.root","RECREATE");  //File for checking histogram fits and integrity
int checkHistos = 0;  //Iterator to save only some histograms for a check
TF1 *fFit;  //Declaration for landau fit check
double pedestals[kNrows][kNcols];
Bool_t gSaturated[kNrows][kNcols];
int badFit = 0; //Counts up per landau fit failing chi^2 test
int DRAWNHISTO = 0; //Counts up per histo drawn


TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
      TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void displayEvent(Int_t entry = -1){
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);

  Int_t r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
    }
  }

  Float_t peak[kNrows][kNcols];
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }

  
  for(Int_t m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      std::cerr << "Why is row negative? Or col?" << std::endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool displayed = false;
    for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      displayed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    if(!displayed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(Int_t s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(gSaturated[r][c])
        histos[r][c]->SetLineColor(kRed+1);
      else
        histos[r][c]->SetLineColor(kBlue+1);
      if(tdc[r][c]!=0){
        histos[r][c]->SetLineColor(kGreen+1);
      }
      if(tdc[r][c]!=0 && !gSaturated[r][c]){ //Must check to make sure that the cols are cols and rows are rows
	
	goodHistoTest(histos[r][c],r,c);
	checkHistos++;
	
      }
    }
  }
}

//Open pedestal text file and read out all pedestal values by HV run
void getPedestal(Int_t run = 989){  //is 0-12 in file along the row or column?---------------------------------------*
  Int_t lineNumber = 0;
  ifstream pedFile(Form("Pedestals_%i.txt",run));

  cout << "Getting pedestal values from run " << run << "." << endl;
  
  int n1;
  double d1, d2;

  int rval, cval;
  string line;
  
  while (getline(pedFile,line))
    {
      if(line.at(0) == '#'){
	continue;
      }

      stringstream ss(line);
      ss >> n1 >> d1 >> d2;
      
      rval = floor(n1/kNrows);
      cval = n1 % kNcols;
      
      pedestals[rval][cval] = d1; //Are the numbers in pedFile counting rows first then columns?---------------*
      
      std::cout << "row = " << rval << ", col = " << cval << ", array val = " << pedestals[rval][cval] << std::endl;
    }

  cout << "Pedestal values extracted." << endl;
  
}

//Saturation and TDC test already applied in display(). Test now for chi^2 goodness of fit to landau.
void goodHistoTest(TH1F *testHisto, int row, int col){

  //Try to improve Landau fits
  //int binmax = testHisto->GetMaximumBin();
  //testHisto->SetBinContent(20,3*(testHisto->GetBinContent(binmax)));

  //Get pedestal value for the run and pmt
  
  double pedVal = pedestals[row][col];

  //Subtract pedestal value from each bin
  for(int b=0 ; b<testHisto->GetNbinsX() ; ++b) {  //Is this how I should be subtracting the pedestal?----------*
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
    if(testHisto->GetBinContent(b) < 0) testHisto->SetBinContent(b,0.0);
  }  
  
  testHisto->Fit("landau", "Q", 0, 19);  //Fit histogram with landau, limits window minus overflow on 20
  fFit=testHisto->GetFunction("landau");

  int dof = 18;  //20 bins -1 per histo, removing last bin to erase overflow -> 19 - 1 = 18
  double rChiSqr = (fFit->GetChisquare())/dof;

  if (rChiSqr > 20.0){  //Is this a good limit?---------------------------------------------------------------*
    if(badFit % 10000 == 0){
      HistosFile->cd();
      testHisto->Write(Form("badFitHisto%d_rChiSqr_gt5",badFit));
      testHisto->Draw("AP");
      DRAWNHISTO++;
    }
    badFit++;
  } else if (rChiSqr < 1){
    if(badFit % 10000 == 0){
      HistosFile->cd();
      testHisto->Write(Form("badFitHisto%d_rChiSqr_lt1",badFit));
      testHisto->Draw("AP");
      DRAWNHISTO++;
      cout << "Chisqr <1" << endl;
    }
    badFit++;
  } else {

    //if(row == 0 && col == 7) cout << "Event passes cuts" << endl;
    
    integratePulse(testHisto, fFit, row, col); //Call integratePulse() for histograms that pass all tests
  }
  /*
  if (checkHistos % 10000 == 0){
  HistosFile->cd();
  testHisto->Write(Form("checkHisto%d_r%d_c%d",iter,row,col));
  testHisto->Draw("AP");
  DRAWNHISTO++;
  }
  */
}


//integrate good signals, return pulse integrated value by PMT and output 2xnevents array
void integratePulse(TH1F *intHisto, TF1 *fFit, int row, int col){
  
  //Should I be integrating the remaining bins or the landau fit?---------------------------------------------*
  //intHisto->Fit("landau", "Q");  //Fit histogram with landau, will need to verify the limits - "Q" not working
  //fFit=intHisto->GetFunction("landau");

  intValues.push_back(fFit->Integral(0,DISP_FADC_SAMPLES));  //Not sure if this integral is working correctly - from 0 to 20 bins - Perhaps integrate to inf
  intValues_row.push_back(row);
  intValues_col.push_back(col);

  //if(row == 0 && col == 0) cout << "ROW AND COL ZERO" << endl; //Not sure why r=0,c=0 is always empty

  if (signalTotal % 10000 == 0){
    HistosFile->cd();
    intHisto->Write(Form("goodFitHisto%d_r%d_c%d",signalTotal,row,col));
    intHisto->Draw("AP");
    DRAWNHISTO++;
  }
  signalTotal++;
}

Int_t cosCal(Int_t run = 989, Int_t event = -1)
{

  getPedestal(run);
  
  //cout << pedestals[6][6] << "!!!!" << endl;

  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.row",hcalt::row);
    T->SetBranchAddress("sbs.hcal.col",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.row",1);
    T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
    for(Int_t r = 0; r < kNrows; r++) {
      for(Int_t c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
        gSaturated[r][c] = false;
      }
    }
  }

  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;
  
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    displayEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals integrated = " << signalTotal << "." << endl;
  cout << "Total good signals integrated (intValues.size()) = " << intValues.size() << "..." << endl;
  
  
  //Now take in three vectors and get average integrated pulse value per pmt (unique r c value)
  double pmtIntTotal[kNrows*kNcols]={0};
  double pmtEvTotal[kNrows*kNcols]={0};

  for(int ir=0; ir<kNrows; ir++){
    
    for(int ic=0; ic<kNcols; ic++){
      
      for(int iev=0; iev<intValues.size(); iev++){
	
	int evRow = intValues_row[iev];
	int evCol = intValues_col[iev];

	if(evRow == ir && evCol == ic){

	  //if(evRow == 0 && evCol == 7) cout << "Event is processed" << endl;
	  
	  pmtIntTotal[ir*kNrows+(ic)] += intValues[iev]; //build value for int tot, one val per pmt
	  pmtEvTotal[ir*kNrows+(ic)] += 1.0;
	  
	}
      }
    }
  }

  cout << "pmtIntTotal[0] = " <<  pmtIntTotal[0]  << "." << endl;

  //cout << "pmtIntTotal.size = " << pmtIntTotal.size << "." << endl;

  cout << "pmtEvTotal[0] = " << pmtEvTotal[0] << "." << endl;

  //cout << "pmtEvTotal.size = " << pmtEvTotal.size << "." << endl;
  
  cout << "Integrated totals and total events per PMT calculated." << endl;

  cout << "Checked histograms = " << checkHistos << "." << endl;
  
  cout << "Total bad fits = " << badFit << "." << endl;

  cout << "Total histograms drawn = " << DRAWNHISTO << "." << endl;

  ofstream outFile;

  outFile.open(Form("intPulseRun_%i.txt",run));

  outFile << "#Integrated pulse data for run " << run << "." << endl;
  outFile << "#Total number of integrated events is " << intValues.size() << "." << endl;
  outFile << "#Row  Column  pmtIntValue  EvTotal" << endl;
  
  for(int e=0; e<(kNrows*kNcols); e++){
    if(pmtEvTotal[e] == 0){
      pmtIntTotal[e]=-1.0;
    } else {
    pmtIntTotal[e] = pmtIntTotal[e]/pmtEvTotal[e];
    }
    cout << "PMT row = " << floor(e/kNrows) << ", col = " << e % kNcols << ": Int Pulse Avg = " << pmtIntTotal[e] << " with " << pmtEvTotal[e] << " events averaged over." << endl;
    outFile << floor(e/kNrows) << "  " << e % kNcols << "  " << pmtIntTotal[e] << "  " << pmtEvTotal[e] << endl;
  }
  
  return 0;
}

