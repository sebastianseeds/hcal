//SSeeds 2020 - Cosmic PMT HV calibration code for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

const int kNrows = 24;
const int kNcols = 12;

const int DISP_MIN_SAMPLE = 0.0;
const int DISP_MAX_SAMPLE = 30.0;
const int DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

int gCurrentEntry = -1;
TChain *T = 0;

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(int);

TH1F *ADCvChannel;
TH1F *TDCvChannel;
TH1F *intADCvChannel;
TH1F *intTDCvChannel;
TH1F *NEVvChannel;
TH1F *pedvChannel;
//TH1F *intPedvChannel;
TH1F *histos[kNrows][kNcols];
TH1F *pedSpec[kNrows][kNcols];
TH1F *intPedSpec[kNrows][kNcols];
TH1F *PMTIntSpec[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *PMTIntSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits, cut on tdc coincidence
TH1F *PMTMaxSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits, cut on tdc coincidence
TF1 *f1;
TF1 *f2;

int signalTotal = 1;
TFile *HistosFile = new TFile("HistosFile.root","RECREATE");  //File for checking histogram fits and integrity
TFile *HistosFilePedestal = new TFile("HistosFilePedestal.root","RECREATE");  //File for checking histogram fits and integrity
double pedestals[kNrows][kNcols];
//double intPedestals[kNrows][kNcols];
double pedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols];
//bool gPulse[kNrows+2][kNcols+2]; //Needs to be larger for false-valued buffers (1 per side)
//bool gPulseTDC[kNrows+2][kNcols+2]; //Needs to be larger for false-valued buffers (1 per side)
bool gPulse[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
bool gPulseTDC[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
int DRAWNHISTO = 0; //Counts up per histo drawn
int rejectNote = 0;
int passNote = 0;
int rejectNoteTDC = 0;
int passNoteTDC = 0;
int histNote = 0;

double targetRAU = 61.425; //Taken from Scott B.
int runN[9] = {989,988,987,984,983,981,980,978,1235}; //All run numbers
//double jlab_lHV[9] = {1600.0,1650.0,1750.0,1850.0,2000.0,1900.0,1800.0,1700.0,1600.0}; //For all runs, left half jlab pmt HV
//double cmu_lHV[9] = {1500.0,1550.0,1650.0,1750.0,1900.0,1800.0,1700.0,1600.0,1500.0}; //For all runs, left half cmu pmt HV

double pmtHV[kNrows][kNcols];

//double jlab_rHV = 1600.0; //Not used since lHV is the same and referenced instead
//double cmu_rHV = 1500.0; //Not used since lHV is the same and referenced instead
//double jlab_lHVset, cmu_lHVset;
double cmuExp = 10.5; //Exp parameter per gain equation on data sheet. Starting values from LED runs
double jlabExp = 8.0; //Exp parameter per gain equation on data sheet. Starting values from LED runs
double targetHV[kNrows][kNcols];
double targetHVTDC[kNrows][kNcols];

double powFit(double *X,double *p) 
{
  double fitval = p[0]*pow(X[0],p[1]);
  //double fitval = p[0]*pow((X[0]-p[1])/1200,p[2]) + p[3]; //X normalized to lowest HV value
  return fitval;
}

//Create skewed normal distribution to fit the timing resolution.
double skewFit(double *X,double *p) 
{
  double fitval = p[0]*exp(-0.5*pow((X[0]-p[1])/p[2],2.0))*0.5*(1+erf(p[3]*((X[0]-p[1])/p[2])));
  return fitval;
}


TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void processEvent(int entry = -1){
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  T->GetEntry(gCurrentEntry);

  //histNote=0;
  
  int r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
      //gPulse[r+1][c+1] = false;
      //gPulseTDC[r+1][c+1] = false;
      gPulse[r+2][c+2] = false;
      gPulseTDC[r+2][c+2] = false;
    }
  }
  
  float peak[kNrows][kNcols];
  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      peak[r][c] = 0.0;
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  for(int m = 0; m < hcalt::ndata; m++) {
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
    bool processed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    if(!processed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,-404);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(TString::Format("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(!gSaturated[r][c]){
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
      }
    }
  }

  //Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. Only building functionality for straight vertical line test, leaving diagonal test machinery for potential implementation.
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gPulse[r+2][c+2]==true){
	if((gPulse[r][c+2]==true&&gPulse[r+1][c+2]==true)||(gPulse[r+1][c+2]==true&&gPulse[r+3][c+2]==true)||(gPulse[r+3][c+2]==true&&gPulse[r+4][c+2]==true)){
	  if(passNote<10) cout << "Pulse r" << r << " c" << c << " passes verticality cut." << endl;
	  passNote++;
	  continue;
	}else{
	  gPulse[r+2][c+2]=false;
	  if(rejectNote<10) cout << "Pulse r" << r << " c" << c << " REJECTED on verticality cut." << endl;
	  rejectNote++;
	}
      }
    }
  }
  
  /*
  //Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. This checks for two adjacent pulses in separate rows including diagonal tracks.
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gPulse[r+1][c+1]==true){
	if(gPulse[r+2][c]==false&&gPulse[r+2][c+1]==false&&gPulse[r+2][c+2]==false){
	  if(gPulse[r][c]==false&&gPulse[r][c+1]==false&&gPulse[r][c+2]==false){
	    gPulse[r+1][c+1]=false;
	    
	    //if(rejectNote<10) cout << "Pulse r" << r << " c" << c << " REJECTED on verticality cut." << endl;
	    rejectNote++;
	  }
	}
      }
      if(gPulse[r+1][c+1]==true&&passNote<10){
	//cout << "Pulse r" << r << " c" << c << " passes verticality cut." << endl;
	passNote++;
      }
      if(gPulseTDC[r+1][c+1]==true){
	if(gPulseTDC[r+2][c]==false&&gPulseTDC[r+2][c+1]==false&&gPulseTDC[r+2][c+2]==false){
	  if(gPulseTDC[r][c]==false&&gPulseTDC[r][c+1]==false&&gPulseTDC[r][c+2]==false){
	    gPulseTDC[r+1][c+1]=false;
	    //if(rejectNoteTDC<10) cout << "TDC cut pulse r" << r << " c" << c << " REJECTED on verticality cut." << endl;
	    rejectNoteTDC++;
	  }
	}
      }
      if(gPulse[r+1][c+1]==true&&passNoteTDC<10){
	//cout << "TDC cut pulse r" << r << " c" << c << " passes verticality cut." << endl;
	passNoteTDC++;
      }
    }
  }
  */  

  //Now, if pulse passes verticality test, pedestal subtract and fill spectra histograms
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      //if(gPulse[r+1][c+1]==true){

      if(gPulse[r+2][c+2]==true){
	
	PMTIntSpec[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	PMTMaxSpec[r][c]->Fill(histos[r][c]->GetMaximum());

	//Only fill TDC if both ADC and TDC pulses present
	if(gPulseTDC[r+2][c+2]==true){
	  PMTIntSpecTDC[r][c]->Fill(histos[r][c]->Integral(1,DISP_MAX_SAMPLE)); //Integral from total bin content
	  PMTMaxSpecTDC[r][c]->Fill(histos[r][c]->GetMaximum());
	}	
      }	
    }
  }
}

//Acquire pedestal from first 4 bins of each histogram for each module
void getPedestal(int entry=-1){ 
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }
  
  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }
  
  T->GetEntry(gCurrentEntry);
  
  int r,c,idx,n,sub;
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
    }
  }
  
  double adc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      adc[r][c] = 0.0;
    }
  }
  
  
  for(int m = 0; m < hcalt::ndata; m++) {
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
    bool processed = false;
    for(int s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
    }
  
    if(!processed) {
      cerr << "Skipping empty module: " << m << endl;
      for(int s = 0;  s < DISP_FADC_SAMPLES; s++) {
        histos[r][c]->SetBinContent(s+1,0);
      }
    }
  }
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      
      double sum = 0.0;


      /*
      for(int b=1 ; b<5 ; b++) {
	sum+=histos[r][c]->GetBinContent(b);
      }
      
      pedSpec[r][c]->Fill(sum/4); //Fill histogram by pmt with the average of the first 4 bins
      //intPedSpec[r][c]->Fill(sum);
      */

      //cout << "Max value is " << histos[r][c]->GetMaximum() << "." << endl;
      
      if(histos[r][c]->GetMaximum()<320.0){  //eliminating adc pulses from calc of pedestal

	for(int b=0 ; b<DISP_MAX_SAMPLE ; b++) {
	  pedSpec[r][c]->Fill(histos[r][c]->GetBinContent(b+1));
	}
	
	/*
	for(int b=0 ; b<DISP_MAX_SAMPLE ; b++) {
	  sum+=histos[r][c]->GetBinContent(b+1);
	}
	
	pedSpec[r][c]->Fill(sum/DISP_MAX_SAMPLE); //Fill histogram by pmt with the average of the first 4 bins
	//intPedSpec[r][c]->Fill(sum);
	
	if(sum/DISP_MAX_SAMPLE>320) cout << "Sum OVER range on " << r << "  " << c << "  " << "." << endl;

	if(sum/DISP_MAX_SAMPLE<120) cout << "Sum UNDER range on " << r << "  " << c << "  " << "." << endl;

	//cout << "Huh?" << endl;

	*/
      }
      
    }
  }
}

//Pedestal subtract then verify that max value is greater than twice sigma of pedestals and lies w/in bounds of histo
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col){
  
  //Get pedestal value for the run and pmt 
  double pedVal = pedestals[row][col];
  
  //Subtract pedestal value from each bin
  for(int b=1 ; b<=testHisto->GetNbinsX() ; b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
  }
 
  //if(testHisto->GetMaximum()>12.0 && testHisto->GetMaximum()>(3.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<28){  //Cut out noise by ensuring all pulses are greater than twice the sigma of the ped and that the max lies in the center of the histogram
  if(testHisto->GetMaximum()>(6.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<28){  //Cut out noise by ensuring all pulses are greater than twice the sigma of the ped and that the max lies in the center of the histogram
    
    if(tdcVal!=0) gPulseTDC[row+2][col+2]=true;
    
    gPulse[row+2][col+2]=true;
    
    if (signalTotal % 5000 == 0){
      cout << "Writing a reference histogram to file.." << endl;
      cout << "Pedestal for same histogram = " << pedestals[row][col] << "." << endl;
      cout << "Pedestal sigma for same histogram = " << pedSigma[row][col] << "." << endl;
      cout << "Maximum value for same histogram = " << testHisto->GetMaximum() << "." << endl;
      cout << "TDC value for same histogram = " << tdcVal << "." << endl << endl;
      
      HistosFile->cd();
      testHisto->SetTitle(Form("Sample R%d C%d",row,col));
      testHisto->GetYaxis()->SetTitle("RAU");
      testHisto->GetYaxis()->CenterTitle();
      testHisto->GetXaxis()->SetTitle("ADC Sample");
      testHisto->GetXaxis()->CenterTitle();
      testHisto->Write(Form("sampleHisto%d_r%d_c%d",signalTotal,row,col));
      testHisto->Draw("AP");
      DRAWNHISTO++;
    }
    
    signalTotal++;
    
  }
}

int cosCal_v4(int run = 989, int event = -1)
{

  //Declare outfile
  TFile *cosAggFile = new TFile(Form("cosmicHistograms_run%d.root",run),"RECREATE");
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      PMTIntSpec[r][c] = new TH1F(Form("Int ADC Spect R%d C%d",r,c),Form("Int ADC Spect R%d C%d",r,c),150,0,4000);
      PMTIntSpec[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),150,0,1000);
      PMTMaxSpec[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      PMTIntSpecTDC[r][c] = new TH1F(Form("Int ADC Spect R%d C%d, TDC Cut",r,c),Form("Int ADC Spect R%d C%d",r,c),150,0,4000);
      PMTIntSpecTDC[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      PMTMaxSpecTDC[r][c] = new TH1F(Form("Max ADC Spect R%d C%d, TDC Cut",r,c),Form("Max ADC Spect R%d C%d",r,c),150,0,1000);
      PMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();
      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,120,320);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
      //intPedSpec[r][c] = new TH1F(Form("Sum Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,550,1200);
      //intPedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      //intPedSpec[r][c]->GetXaxis()->CenterTitle();
      
    }
  }

  //Build extra analysis histograms
  ADCvChannel = new TH1F("ADCvChannel", "ADCvChannel", 288, 0, 287);
  TDCvChannel = new TH1F("TDCvChannel", "TDCvChannel", 288, 0, 287);
  intADCvChannel = new TH1F("intADCvChannel", "intADCvChannel", 288, 0, 287);
  intTDCvChannel = new TH1F("intTDCvChannel", "intTDCvChannel", 288, 0, 287);
  NEVvChannel = new TH1F("NEVvChannel", "NEVvChannel", 288, 0, 287);
  pedvChannel = new TH1F("pedvChannel", "pedvChannel", 288, 0, 287);
  //intPedvChannel = new TH1F("intPedvChannel", "intPedvChannel", 288, 0, 287);
  
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d_0.root",run));
    
    T->Add(TString::Format("fadc_f1tdc_%d_1.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_2.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_3.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_4.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_5.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_6.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_7.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_8.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_9.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_10.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_11.root",run));
    T->Add(TString::Format("fadc_f1tdc_%d_12.root",run));
    
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
    for(int r = 0; r < kNrows; r++) {
      for(int c = 0; c < kNcols; c++) {
        histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
        gSaturated[r][c] = false;
      }
    }
  }
  /*
  //Set appropriate HVs for run
  for(int i=0; i<9; i++){
    if(run==runN[i]){
      jlab_lHVset = jlab_lHV[i];
      cmu_lHVset = cmu_lHV[i];
    }
  }
  */

  //Set appropriate HVs for run. Settings from HCAL wiki. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream file(Form("HV_run%d.txt",run));

  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while (getline(file,line))
    {
      if(line.at(0) == '#'){
	continue;
      }

      stringstream ss(line);
      ss >> d1;
      
      rval = floor(n1/kNcols);
      cval = n1 % kNcols;
      
      pmtHV[rval][cval] = -d1; //Are the numbers in pedFile counting rows first then columns?---------------*
      
      cout << "row = " << rval << ", col = " << cval << ", array val = " << pmtHV[rval][cval] << endl;

      n1++;
    }
  
  //Set default values for pulse check histograms. Arrays contain false-valued buffer on all sides
  for(int r = 0; r < kNrows+4; r++) {
    for(int c = 0; c < kNcols+4; c++) {
      gPulse[r][c] = false;
      gPulseTDC[r][c] = false;
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;

  cout << "Filling pedestal histograms.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current pedestal event: " << gCurrentEntry << endl;
    }
  }
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      HistosFilePedestal->cd();
      pedSpec[r][c]->SetTitle(Form("Sample R%d C%d",r,c));
      pedSpec[r][c]->Write(Form("sampleHisto_r%d_c%d",r,c));
      pedSpec[r][c]->Draw("AP");
    }
  }
  
  cout << "Fitting pedestal histos and extracting mean/sigma.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      pedSpec[r][c]->Fit("gaus","Q");
      //intPedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      //f2=intPedSpec[r][c]->GetFunction("gaus");
      pedestals[r][c]=f1->GetParameter(1);
      pedSigma[r][c]=f1->GetParameter(2);
      //intPedestals[r][c]=f2->GetParameter(1);
    }
  }

  gCurrentEntry = event; //Resetting entries for pulse analysis
  
  cout << "Processing pedestal subtracted pulses.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ //Not sure about GetEntries() here
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }
  
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;

  cout << "Creating file to hold fit parameters.." << endl;
  ofstream outFile;
  outFile.open(Form("cosCalOut_run%d.txt",run));
  outFile << "#Target HV settings for run " << run << "." << endl;
  outFile << "#Row Col targetHV targetHV_TDCCut" << endl;

  cout << "Writing spectrum histograms and calibration constants to file.." << endl;

  TF1 *skewFitMax[kNrows][kNcols] = {};
  TF1 *skewFitMaxTDC[kNrows][kNcols] = {};
  TF1 *skewFitInt[kNrows][kNcols] = {};
  TF1 *skewFitIntTDC[kNrows][kNcols] = {};
  
  double fitMaxX;
  double TDCfitMaxX;
  double fitIntX;
  double TDCfitIntX;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      skewFitMax[r][c] = new TF1();
      skewFitMaxTDC[r][c] = new TF1();
      skewFitInt[r][c] = new TF1();
      skewFitIntTDC[r][c] = new TF1();
      
      //Create skew fit function per r and c
      skewFitMax[r][c] = new TF1(Form("SkewFitMax r%d c%d",r,c),skewFit, 0, 1000, 4);
      skewFitMax[r][c]->SetLineColor(4);
      skewFitMax[r][c]->SetNpx(1000);
      skewFitMaxTDC[r][c] = new TF1(Form("SkewFitMaxTDC r%d c%d",r,c),skewFit, 0, 1000, 4);
      skewFitMaxTDC[r][c]->SetLineColor(4);
      skewFitMaxTDC[r][c]->SetNpx(1000);
      skewFitInt[r][c] = new TF1(Form("SkewFitInt r%d c%d",r,c),skewFit, 0, 4000, 4);
      skewFitInt[r][c]->SetLineColor(4);
      skewFitInt[r][c]->SetNpx(1000);
      skewFitIntTDC[r][c] = new TF1(Form("SkewFitIntTDC r%d c%d",r,c),skewFit, 0, 4000, 4);
      skewFitIntTDC[r][c]->SetLineColor(4);
      skewFitIntTDC[r][c]->SetNpx(1000);

      //Set parameters for function from spectrum histograms
      skewFitMax[r][c]->SetRange(PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())-0.8*PMTMaxSpec[r][c]->GetStdDev(),PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin())+1.3*PMTMaxSpec[r][c]->GetStdDev());
      skewFitMax[r][c]->SetParameter(0,PMTMaxSpec[r][c]->GetBinContent(PMTMaxSpec[r][c]->GetMaximumBin()));;
      skewFitMax[r][c]->SetParameter(1,PMTMaxSpec[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpec[r][c]->GetMaximumBin()));
      skewFitMax[r][c]->SetParameter(2,PMTMaxSpec[r][c]->GetStdDev());
      skewFitMax[r][c]->SetParameter(3,2);
      skewFitMaxTDC[r][c]->SetRange(PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())-0.8*PMTMaxSpecTDC[r][c]->GetStdDev(),PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin())+1.3*PMTMaxSpecTDC[r][c]->GetStdDev());
      skewFitMaxTDC[r][c]->SetParameter(0,PMTMaxSpecTDC[r][c]->GetBinContent(PMTMaxSpecTDC[r][c]->GetMaximumBin()));;
      skewFitMaxTDC[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      skewFitMaxTDC[r][c]->SetParameter(2,PMTMaxSpecTDC[r][c]->GetStdDev());
      skewFitMaxTDC[r][c]->SetParameter(3,2);
      skewFitInt[r][c]->SetRange(PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())-0.5*PMTIntSpec[r][c]->GetStdDev(),PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin())+1.3*PMTIntSpec[r][c]->GetStdDev());
      skewFitInt[r][c]->SetParameter(0,PMTIntSpec[r][c]->GetBinContent(PMTIntSpec[r][c]->GetMaximumBin()));;
      skewFitInt[r][c]->SetParameter(1,PMTIntSpec[r][c]->GetXaxis()->GetBinCenter(PMTIntSpec[r][c]->GetMaximumBin()));
      skewFitInt[r][c]->SetParameter(2,PMTIntSpec[r][c]->GetStdDev());
      skewFitInt[r][c]->SetParameter(3,2);
      skewFitIntTDC[r][c]->SetRange(PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())-0.5*PMTIntSpecTDC[r][c]->GetStdDev(),PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin())+1.3*PMTIntSpecTDC[r][c]->GetStdDev());
      skewFitIntTDC[r][c]->SetParameter(0,PMTIntSpecTDC[r][c]->GetBinContent(PMTIntSpecTDC[r][c]->GetMaximumBin()));;
      skewFitIntTDC[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      skewFitIntTDC[r][c]->SetParameter(2,PMTIntSpecTDC[r][c]->GetStdDev());
      skewFitIntTDC[r][c]->SetParameter(3,2);
      
      int iHV=0; 
      if(c>3&&c<8) iHV=1;

      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      
      cosAggFile->cd();
      pedSpec[r][c]->SetTitle(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->GetXaxis()->SetTitle("RAU");
      pedSpec[r][c]->GetXaxis()->CenterTitle();
      pedSpec[r][c]->Write(Form("Pedestal Spect R%d C%d",r,c));
      pedSpec[r][c]->Draw("AP");
      /*
      intPedSpec[r][c]->Fit("gaus","Q");
      f1=intPedSpec[r][c]->GetFunction("gaus");
      
      cosAggFile->cd();
      intPedSpec[r][c]->SetTitle(Form("Sum Pedestal 4 Bins Spect R%d C%d",r,c));
      intPedSpec[r][c]->GetXaxis()->SetTitle("sRAU");
      intPedSpec[r][c]->GetXaxis()->CenterTitle();
      intPedSpec[r][c]->Write(Form("Sum Pedestal 4 Bins Spect R%d C%d",r,c));
      intPedSpec[r][c]->Draw("AP");
      */

      //Spectra without TDC cut
      if( PMTIntSpec[r][c]->GetEntries() > 0){
	//cout << "PMTIntSpec entries for row " << r << ", col " << c << " = " << PMTIntSpec[r][c]->GetEntries() << "." << endl;
	PMTIntSpec[r][c]->Fit(skewFitInt[r][c],"RQ");
	fitIntX = skewFitInt[r][c]->GetMaximumX();
	cosAggFile->cd();
	PMTIntSpec[r][c]->SetTitle(Form("Int ADC Spect R%d C%d MaxX%f",r,c,fitIntX));
	PMTIntSpec[r][c]->Write(Form("Int ADC Spect R%d C%d",r,c));
	PMTIntSpec[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpec[r][c]->GetEntries() > 0){
	//cout << "PMTMaxSpec entries for row " << r << ", col " << c << " = " << PMTMaxSpec[r][c]->GetEntries() << "." << endl;
	PMTMaxSpec[r][c]->Fit(skewFitMax[r][c],"RQ");
	fitMaxX = skewFitMax[r][c]->GetMaximumX();
	cosAggFile->cd();
	PMTMaxSpec[r][c]->SetTitle(Form("Max ADC Spect R%d C%d MaxX%f",r,c,fitMaxX));
	PMTMaxSpec[r][c]->Write(Form("Max ADC Spect R%d C%d",r,c));
	PMTMaxSpec[r][c]->Draw("AP");
      }
      
      //Calculate target HV from gain equation on data sheet
      
      if(iHV==1){
	targetHV[r][c] = pmtHV[r][c]/pow(fitMaxX/targetRAU,1.0/jlabExp);
      }else{
	targetHV[r][c] = pmtHV[r][c]/pow(fitMaxX/targetRAU,1.0/cmuExp);
      }
      
      
      //Debugging
      /*
	if(iHV==1){
	targetHV[r][c] = pmtHV[r][c];
	}else{
	targetHV[r][c] = pmtHV[r][c];
	}
      */
      cout << "fitMaxX = " << fitMaxX << endl;
      
      //Spectra with TDC cut
      if( PMTIntSpecTDC[r][c]->GetEntries() > 0){
	//cout << "PMTIntSpecTDC entries for row " << r << ", col " << c << " = " << PMTIntSpecTDC[r][c]->GetEntries() << "." << endl;
	PMTIntSpecTDC[r][c]->Fit(skewFitIntTDC[r][c],"RQ");
	TDCfitIntX = skewFitIntTDC[r][c]->GetMaximumX();
	cosAggFile->cd();
	PMTIntSpecTDC[r][c]->SetTitle(Form("Int ADC Spect TDCcut R%d C%d MaxX%f",r,c,TDCfitIntX));
	PMTIntSpecTDC[r][c]->Write(Form("Int ADC Spect TDC R%d C%d",r,c));
	PMTIntSpecTDC[r][c]->Draw("AP");
      }
      
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
	//cout << "PMTMaxSpecTDC entries for row " << r << ", col " << c << " = " << PMTMaxSpecTDC[r][c]->GetEntries() << "." << endl;
	PMTMaxSpecTDC[r][c]->Fit(skewFitMaxTDC[r][c],"RQ");
	TDCfitMaxX = skewFitMaxTDC[r][c]->GetMaximumX();
	cosAggFile->cd();
	PMTMaxSpecTDC[r][c]->SetTitle(Form("Max ADC Spect TDCcut R%d C%d MaxX%f",r,c,TDCfitMaxX));
	PMTMaxSpecTDC[r][c]->Write(Form("Max ADC Spect TDC R%d C%d",r,c));
	PMTMaxSpecTDC[r][c]->Draw("AP");
      }

      //Calculate target HV with TDC cut from gain equation on data sheet

      if(iHV==1){
	targetHVTDC[r][c] = pmtHV[r][c]/pow(TDCfitMaxX/targetRAU,1.0/jlabExp);
      }else{
	targetHVTDC[r][c] = pmtHV[r][c]/pow(TDCfitMaxX/targetRAU,1.0/cmuExp);
      }

      //Fill Max ADC vs Channel Histogram
      ADCvChannel->SetBinContent(kNcols*r+c+1,fitMaxX);
      TDCvChannel->SetBinContent(kNcols*r+c+1,TDCfitMaxX);
      intADCvChannel->SetBinContent(kNcols*r+c+1,fitIntX);
      intTDCvChannel->SetBinContent(kNcols*r+c+1,TDCfitIntX);

      //Fill Number of events per module histogram
      //int NEV = PMTIntSpec[r][c]->GetEntries();
      NEVvChannel->SetBinContent(kNcols*r+c+1,PMTIntSpec[r][c]->GetEntries());

      //Fill average pedestal value per module histograms
      pedvChannel->SetBinContent(kNcols*r+c+1,pedestals[r][c]);
      //intPedvChannel->SetBinContent(kNcols*r+c,intPedestals[r][c]);

      cout << "For row " << r << " and col " << c << ", target HV is " << targetHV[r][c] << " and target HV (TDC cut) is " << targetHVTDC[r][c] << "." << endl;
      outFile << r << "  " << c << "  " << targetHV[r][c] << "  " << targetHVTDC[r][c] << endl;
      
    }
  }

  //Draw analysis histograms
  cosAggFile->cd();
  ADCvChannel->SetTitle("Average Max ADC Val vs Channel");
  ADCvChannel->Write("ADCvChannel");
  ADCvChannel->Draw("AP");
  
  cosAggFile->cd();
  TDCvChannel->SetTitle("Average Max ADC Val with TDC Cut vs Channel");
  TDCvChannel->Write("TDCvChannel");
  TDCvChannel->Draw("AP");
  
  cosAggFile->cd();
  intADCvChannel->SetTitle("Average Int ADC Val vs Channel");
  intADCvChannel->Write("intADCvChannel");
  intADCvChannel->Draw("AP");
  
  cosAggFile->cd();
  intTDCvChannel->SetTitle("Average Int ADC Val with TDC Cut vs Channel");
  intTDCvChannel->Write("intTDCvChannel");
  intTDCvChannel->Draw("AP");
  
  cosAggFile->cd();
  NEVvChannel->SetTitle("Number of Event Pulses vs Channel");
  NEVvChannel->Write("NEVvChannel");
  NEVvChannel->Draw("AP");

  cosAggFile->cd();
  pedvChannel->SetTitle("Avg Pedestal vs Channel");
  pedvChannel->Write("pedvChannel");
  pedvChannel->Draw("AP");
  /*
  cosAggFile->cd();
  intPedvChannel->SetTitle("Avg Sum Pedestal (4 first bins) vs Channel");
  intPedvChannel->Write("intPedvChannel");
  intPedvChannel->Draw("AP");
  */
  outFile << endl << endl;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      
      outFile << "Pedestal for row " << r << " and col " << c << " = " << pedestals[r][c] << endl;

    }
  }
  
  cout << "ADC spectrum histograms written to file cosSpectFile_run" << run << ".root" << endl;
  cout << "Sample histograms drawn to file HistosFile.root" << endl;
  cout << "Target HV settings written to cosCalOut_run" << run << ".txt" << endl;
  
  return 0;
}
  
