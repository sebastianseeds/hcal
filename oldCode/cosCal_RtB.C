//SSeeds April 2021 - Cosmic PMT HV calibration code version 3 for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include "hcal.h"

//Detector parameters and flags
const int kNrows = 24;
const int kNcols = 12;

const int minSample = 0.0;
const int maxSample = 50.0;
const int totSample = (maxSample-minSample); //Should be the total fADC window size with each samp = 4ns

//Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int);
void getPedestal(int);

//Declare necessary histograms
TH1F *histos[kNrows][kNcols];
TH1F *pedSpec[kNrows][kNcols];
TH1F *retToBase[kNrows][kNcols]; //=new TH1F("","",totSample,0,totSample); //Limits set by histos[][]
TH1F *PMTIntSpec[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *PMTIntSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,5000); //Empirical limits, cut on tdc coincidence
TH1F *PMTMaxSpecTDC[kNrows][kNcols]; //=new TH1F("","",20,0,1000); //Empirical limits, cut on tdc coincidence

//Declare fit function
TF1 *f1;

//Create marker to track number of signals that pass cuts
int signalTotal = 1;

//Create files to keep temporary histograms for integrity checks
TFile *HistosFile = new TFile("outFiles/CosHistosFile.root","RECREATE");  //File for checking histogram fits and integrity
TFile *HistosFilePedestal = new TFile("outFiles/CosHistosFilePedestal.root","RECREATE");  //File for checking histogram fits and integrity

//Declare necessary arrays
double pedestals[kNrows][kNcols];
//double retToBase[kNrows][kNcols];
double pedSigma[kNrows][kNcols];
bool gSaturated[kNrows][kNcols];
bool gPulse[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
bool gPulseTDC[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)
bool gVert[kNrows+4][kNcols+4]; //Needs to be larger for false-valued buffers (2 per side)

//double targetRAU = 48.15; //4095*3.0/2.0*0.39*14/700 or (digMaxADC)*(dynRanAmp/dynRanADC)*(long cable attenuation)*(cosEdep/maxGMnEdep). Attenuation average over all channels current April 2021 pending effective gain measurements in Hall.

double targetRAU = 78.355; //Previous value of 61.425*(2.55/2) or targetRAU_prev*(attLongCable_newMeas/attLongCable_oldMeas)

double pmtHV[kNrows][kNcols];
double alphas[kNrows][kNcols];
double targetHV[kNrows][kNcols];

//Landau fit. May improve with convolution with dying exponential to deal with low amplitude EM effects.
double landauFit(double *x, double *par) {
  double amp = par[0];
  double mpv = par[1];
  double sigma = par[2];
  double offset = par[3];
  double ADC = x[0];
  return amp * TMath::Landau(ADC,mpv,sigma) + offset;
}

//Create generic histogram function
TH1F* MakeHisto(int row, int col, int bins){
  TH1F *h = new TH1F(Form("h%02d%02d",row,col),Form("%d-%d",row+1,col+1),bins,minSample,maxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//
void processEvent(int entry = -1){
  //Check event increment and increment
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  //Get the event from the TTree
  T->GetEntry(gCurrentEntry);
  
  int r,c,idx,n,sub;
  // Clear old signal histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
      gSaturated[r][c] = false;
      gPulse[r+2][c+2] = false;
      gPulseTDC[r+2][c+2] = false;
      gVert[r+2][c+2] = false;
    }
  }

  //Reset signal peak, adc, and tdc arrays
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

  //Process event with m data
  for(int m = 0; m < hcalt::ndata; m++) {
    //Define row and column
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      cerr << "Why is row negative? Or col?" << endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols) continue;

    //Define index, number of samples, fill adc and tdc arrays, and switch processed marker for error reporting
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;

    //Fill signal histogram from samps and mark saturated array if applicable
    for(int s = minSample; s < maxSample && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]);
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c] = true;
      }
    }
    //Report error if module is empty
    if(!processed) {
      std::cerr << "Skipping empty module: " << m << std::endl;
      for(int s = minSample;  s < totSample; s++) {
        histos[r][c]->SetBinContent(s+1-minSample,-404);
      }
    }
  }

  //Pass unsaturated signal histos to goodHistoTest to see if a signal pulse exists there. gPulse array value  marked true if test passed.
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->SetTitle(Form("%d-%d (ADC=%g,TDC=%g)",r+1,c+1,adc[r][c],tdc[r][c]));
      if(!gSaturated[r][c]){
	goodHistoTest(histos[r][c],tdc[r][c],r,c);
      }
    }
  }

  //Vertical line test - remember to shift all values up by one due to intialization of gPulse/gPulseTDC and 'false' buffer. Only building functionality for straight vertical line test.
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gPulse[r+2][c+2]==true){
	if((gPulse[r][c+2]==true&&gPulse[r+1][c+2]==true)||(gPulse[r+1][c+2]==true&&gPulse[r+3][c+2]==true)||(gPulse[r+3][c+2]==true&&gPulse[r+4][c+2]==true)){ //Checks if two pulses exist above, below, or one above and below every given pulse to ensure a track exists.
	  if(gPulse[r+2][c+1]==false&&gPulse[r+2][c+3]==false) //Check on each side of track for pulses to exclude events which clip the corner of the module.
	    gVert[r+2][c+2]=true; //Flip the bool to determine a good cosmic pulse.
	}else{ //Else excluded. Diagonal tracks also excluded.
	  continue;
	}
      }
    }
  }

  //Now, if pulse passes verticality test, pedestal subtract and fill spectra histograms
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gVert[r+2][c+2]==true){ 
	
	//Get return to baseline parameter (in bins (ADC samples, 4ns per sample))
	double baseline = 0.0;
	//int return = 0.0;
	for(int s = 0; s < 4; s++){ //Get baseline from first four bins
	  
	  baseline += histos[r][c]->GetBinContent(s+1-minSample);
	  
	}
	baseline /= 4.0; //Get average value of the baseline per event

	for(int s = histos[r][c]->GetMaximumBin(); s < totSample; s++){
	  double binContent = histos[r][c]->GetBinContent(s);
	  if(binContent<baseline+3*pedSigma[r][c]){
	    //cout << s-histos[r][c]->GetMaximumBin() << endl;
	    retToBase[r][c]->Fill(s-histos[r][c]->GetMaximumBin());
	  }else if(s==totSample-1){
	    retToBase[r][c]->Fill(50);
	  }else{
	    continue;
	  }
	}
	
	
	PMTIntSpec[r][c]->Fill(histos[r][c]->Integral(1,maxSample)); //Integral from total bin content
	PMTMaxSpec[r][c]->Fill(histos[r][c]->GetMaximum());
	
	//Only fill TDC if both ADC and TDC pulses present
	if(gPulseTDC[r+2][c+2]==true){
	  PMTIntSpecTDC[r][c]->Fill(histos[r][c]->Integral(1,maxSample)); //Integral from total bin content
	  PMTMaxSpecTDC[r][c]->Fill(histos[r][c]->GetMaximum());
	}	
      }	
    }
  }
}

//Acquire pedestal for each module on events where no cosmic passed through said module (TDC==0). See comments for processEvent().
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
  
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c]->Reset("ICES M");
    }
  }

  double adc[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      adc[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }
  
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      cerr << "Why is row negative? Or col?" << endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols) continue;
    
    idx = hcalt::samps_idx[m];
    n = hcalt::nsamps[m];
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s = minSample; s < maxSample && s < n; s++) {
      processed = true;
      histos[r][c]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]);
    }
  
    if(!processed) {
      cerr << "Skipping empty module: " << m << endl;
      for(int s = minSample;  s < totSample; s++) {
        histos[r][c]->SetBinContent(s+1-minSample,0);
      }
    }
  }
  
  //With a histogram populated for each module on this event, fill pedestal histogram per module over all events with contents
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(tdc[r][c]==0){ //eliminating adc pulses from coincident tdc measurement in pedestal calculation	
	for(int b=0 ; b<maxSample ; b++) {	  
	  pedSpec[r][c]->Fill(histos[r][c]->GetBinContent(b+1));
	}
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

  //Primary signal cut. All defined signals are greater than twice the sigma of the pedestal and the maximum value is not in the first or last bin
  if(testHisto->GetMaximum()>(6.0*pedSigma[row][col]) && testHisto->GetMaximumBin()>2 && testHisto->GetMaximumBin()<28){

    //Switch if both ADC signal and TDC signal exists
    if(tdcVal!=0) gPulseTDC[row+2][col+2]=true;

    //Switch if ADC signal exists
    gPulse[row+2][col+2]=true;

    //Write out sample event histograms for independent verification
    if (signalTotal % 50000 == 0){
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
    }
    
    signalTotal++;
    
  }
}

//Main
int cosCal_RtB(int run = 1725, int event = -1){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);
  
  //Declare outfile
  TFile *cosHistFile = new TFile(Form("outFiles/cosHistRtB_run%d.root",run),"RECREATE");
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building ADC and TDC spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      PMTIntSpec[r][c] = new TH1F(Form("Int ADC Spect R%d C%d",r,c),Form("Int ADC Spect R%d C%d",r,c),100,0,4000);
      PMTIntSpec[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      
      PMTMaxSpec[r][c] = new TH1F(Form("Max ADC Spect R%d C%d",r,c),Form("Max ADC Spect R%d C%d",r,c),100,0,1000);
      PMTMaxSpec[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      
      PMTIntSpecTDC[r][c] = new TH1F(Form("Int ADC Spect R%d C%d, TDC Cut",r,c),Form("Int ADC Spect R%d C%d",r,c),100,0,4000);
      PMTIntSpecTDC[r][c]->GetXaxis()->SetTitle("sRAU");
      PMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      PMTMaxSpecTDC[r][c] = new TH1F(Form("Max ADC Spect R%d C%d, TDC Cut",r,c),Form("Max ADC Spect R%d C%d",r,c),100,0,1000);
      PMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle("RAU");
      PMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();

      pedSpec[r][c] = new TH1F(Form("Pedestal Spect R%d C%d",r,c),Form("Pedestal Spect R%d C%d",r,c),200,120,320);
      pedSpec[r][c]->GetXaxis()->SetTitle("<RAU>");
      pedSpec[r][c]->GetXaxis()->CenterTitle();

      retToBase[r][c] = new TH1F(Form("Return to Baseline R%d C%d",r,c),Form("Return to Baseline R%d C%d",r,c),totSample,minSample,maxSample);
      retToBase[r][c]->GetXaxis()->SetTitle("Samples");
      retToBase[r][c]->GetXaxis()->CenterTitle();
    }
  }

  //Read in data produced by analyzer in root format
  cout << "Reading raw data from analyzer.." << endl;
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("rootFiles/cosmic/fadc_f1tdc_%d_1.root",run)); //DEBUGGING
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
        histos[r][c] = MakeHisto(r,c,totSample);
        gSaturated[r][c] = false;
      }
    }
  }

  //Set appropriate HV and alphas for run. HV settings from HCAL wiki. Alphas from LED analysis. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream file(Form("setFiles/HV_run%d.txt",run));

  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while(getline(file,line)){
    if(line.at(0) == '#'){
      continue;
    }
    
    stringstream ss(line);
    ss >> d1;
    
    rval = floor(n1/kNcols);
    cval = n1 % kNcols;
    
    pmtHV[rval][cval] = -d1; 
        
    n1++;
  }

  ifstream file2("setFiles/alphas.txt");
  
  n1=0;
  string line2;

  while(getline(file2,line2)){
    if(line2.at(0)=='#'){
      continue;
    }

    stringstream ss(line2);
    ss >> d1;

    rval = floor(n1/kNcols);
    cval = n1 % kNcols;

    alphas[rval][cval] = d1;
    
    n1++;
  }
  
  //Set default values for pulse check bool arrays. Arrays contain false-valued buffer for possible future diagonal track implementation.
  cout << "Resetting particle track test flags.." << endl;
  for(int r = 0; r < kNrows+4; r++) {
    for(int c = 0; c < kNcols+4; c++) {
      gPulse[r][c] = false;
      gPulseTDC[r][c] = false;
      gVert[r][c] = false;
    }
  }
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;

  //Fill pedestal histograms from data where f1TDC did not fire
  cout << "Filling pedestal histograms.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
    getPedestal(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current pedestal event: " << gCurrentEntry << endl;
    }
  }

  /*
  //Write pedestal histograms to file
  cout << "Writing pedestal histograms to file.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){  
      HistosFilePedestal->cd();
      pedSpec[r][c]->SetTitle(Form("Pedestal R%d C%d",r,c));
      pedSpec[r][c]->Write(Form("pedHisto_r%d_c%d",r,c));
      pedSpec[r][c]->Draw("AP");
    }
  }
  */

  //Fitting pedestal histograms and getting mean to subtract event pulses by module
  cout << "Processing pedestals and saving to file.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      pedSpec[r][c]->Fit("gaus","Q");
      f1=pedSpec[r][c]->GetFunction("gaus");
      pedestals[r][c]=f1->GetParameter(1);
      pedSigma[r][c]=f1->GetParameter(2);
      HistosFilePedestal->cd();
      pedSpec[r][c]->SetTitle(Form("Pedestal R%d C%d",r,c));
      pedSpec[r][c]->Write(Form("pedHisto_r%d_c%d",r,c));
      pedSpec[r][c]->Draw("AP");
    }
  }

  gCurrentEntry = event; //Resetting entries for pulse analysis

  //Pedestal subtract each pulse by module and populate fADC and f1TDC spectra histograms
  cout << "Pedestal subtracting and processing signals.." << endl;
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
    processEvent(gCurrentEntry);
    gCurrentEntry++;
    
    //Keep count of events processed for monitoring
    if (gCurrentEntry%10000 == 0){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }

  //Fit ADC spectra with landau to extract mpv for HV calibration
  cout << "Writing spectrum histograms and calibration constants to file.." << endl;
  ofstream outFile;
  outFile.open(Form("outFiles/HVTargetsV3_run%d.txt",run)); //Create text file to hold target HVs
  time_t now = time(0); 
  char *dt = ctime(&now);
  outFile << "#Target HV settings for run " << run << " from cosCal_v3 on " << dt << "#" << endl;
  outFile << "#Row Col targetHV" << endl;

  //Make array of fit functions to pass initial conditions to and to fit ADC spectra per channel
  TF1 *landauFitMax[kNrows][kNcols] = {};
  TF1 *landauFitInt[kNrows][kNcols] = {};
  TF1 *landauFitMaxTrue[kNrows][kNcols] = {};
  TF1 *landauFitIntTrue[kNrows][kNcols] = {};
  

  double fitMaxX;
  double fitIntX;
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      //Fit twice - first across entire ADC range, second restricted with parameters from first fit

      //First Fit

      //Create landau fit function per r and c
      landauFitMax[r][c] = new TF1(Form("landauFitMax r%d c%d",r,c),landauFit, 0, 1000, 4);  //Can use user defined landauFit here, but worst results
      landauFitMax[r][c]->SetLineColor(4);
      landauFitMax[r][c]->SetNpx(1000);

      landauFitInt[r][c] = new TF1(Form("landauFitInt r%d c%d",r,c),landauFit, 0, 4000, 4);
      landauFitInt[r][c]->SetLineColor(4);
      landauFitInt[r][c]->SetNpx(1000);

      //MAX
      //Set parameters for landau fit with maximized range
      landauFitMax[r][c]->SetRange(0,1000);
      landauFitMax[r][c]->SetParameter(0,PMTMaxSpecTDC[r][c]->GetBinContent(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      //Use TDC to get good guess at mpv for landau fit
      landauFitMax[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));

      //landauFitMax[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetMaximumBin());
      landauFitMax[r][c]->SetParameter(2,PMTMaxSpecTDC[r][c]->GetStdDev());
      landauFitMax[r][c]->SetParameter(3,0);

      //INT
      landauFitInt[r][c]->SetRange(0,4000);
      landauFitInt[r][c]->SetParameter(0,PMTIntSpecTDC[r][c]->GetBinContent(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      //Use TDC to get good guess at mpv for landau fit
      landauFitInt[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      //landauFitInt[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetMaximumBin());
      landauFitInt[r][c]->SetParameter(2,PMTIntSpecTDC[r][c]->GetStdDev());
      landauFitInt[r][c]->SetParameter(3,0);

      //FIT
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
	PMTMaxSpecTDC[r][c]->Fit(landauFitMax[r][c],"+RQ");
	if(landauFitMax[r][c]->GetParameter(1)>0&&landauFitMax[r][c]->GetParameter(1)<2000){
	  fitMaxX = landauFitMax[r][c]->GetParameter(1);
	  landauFitMax[r][c]->SetLineColor(kBlue);
	}else{
	  fitMaxX = PMTMaxSpecTDC[r][c]->GetMean();
	  landauFitMax[r][c]->SetLineColor(kRed);
	}
	cosHistFile->cd();
	PMTMaxSpecTDC[r][c]->SetTitle(Form("Max TDC Spect R%d C%d MaxX%f",r,c,fitMaxX));
	PMTMaxSpecTDC[r][c]->Write(Form("Max TDC Spect R%d C%d",r,c));
	PMTMaxSpecTDC[r][c]->Draw("AP");
      }


      //Second fit
      //Create landau fit function per r and c
      landauFitMaxTrue[r][c] = new TF1(Form("landauFitMax r%d c%d",r,c),landauFit, 0, 1000, 4);
      landauFitMaxTrue[r][c]->SetLineColor(4);
      landauFitMaxTrue[r][c]->SetNpx(1000);

      landauFitIntTrue[r][c] = new TF1(Form("landauFitInt r%d c%d",r,c),landauFit, 0, 4000, 4);
      landauFitIntTrue[r][c]->SetLineColor(4);
      landauFitIntTrue[r][c]->SetNpx(1000);

      //MAX
      //Set parameters for landau fit from previous fit
      landauFitMaxTrue[r][c]->SetRange(landauFitMax[r][c]->GetParameter(1)-landauFitMax[r][c]->GetParameter(2),landauFitMax[r][c]->GetParameter(1)+landauFitMax[r][c]->GetParameter(2));
      landauFitMaxTrue[r][c]->SetParameter(0,PMTMaxSpecTDC[r][c]->GetBinContent(PMTMaxSpecTDC[r][c]->GetMaximumBin()));
      //Use TDC to get good guess at mpv for landau fit
      landauFitMaxTrue[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTMaxSpecTDC[r][c]->GetMaximumBin()));

      //landauFitMaxTrue[r][c]->SetParameter(1,PMTMaxSpecTDC[r][c]->GetMaximumBin());
      landauFitMaxTrue[r][c]->SetParameter(2,PMTMaxSpecTDC[r][c]->GetStdDev());
      landauFitMaxTrue[r][c]->SetParameter(3,0);

      //INT
      landauFitIntTrue[r][c]->SetRange(landauFitInt[r][c]->GetParameter(1)-landauFitInt[r][c]->GetParameter(2),landauFitInt[r][c]->GetParameter(1)+landauFitInt[r][c]->GetParameter(2));
      landauFitIntTrue[r][c]->SetParameter(0,PMTIntSpecTDC[r][c]->GetBinContent(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      //Use TDC to get good guess at mpv for landau fit
      landauFitIntTrue[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetXaxis()->GetBinCenter(PMTIntSpecTDC[r][c]->GetMaximumBin()));
      //landauFitIntTrue[r][c]->SetParameter(1,PMTIntSpecTDC[r][c]->GetMaximumBin());
      landauFitIntTrue[r][c]->SetParameter(2,PMTIntSpecTDC[r][c]->GetStdDev());
      landauFitIntTrue[r][c]->SetParameter(3,0);

      //FIT
      if( PMTMaxSpecTDC[r][c]->GetEntries() > 0){
	PMTMaxSpecTDC[r][c]->Fit(landauFitMaxTrue[r][c],"+RQ");
	if(landauFitMaxTrue[r][c]->GetParameter(1)>0&&landauFitMaxTrue[r][c]->GetParameter(1)<2000){
	  fitMaxX = landauFitMaxTrue[r][c]->GetParameter(1);
	  landauFitMaxTrue[r][c]->SetLineColor(kBlue);
	}else{
	  fitMaxX = PMTMaxSpecTDC[r][c]->GetMean();
	  landauFitMaxTrue[r][c]->SetLineColor(kRed);
	}
	cosHistFile->cd();
	PMTMaxSpecTDC[r][c]->SetTitle(Form("Max TDC Spect R%d C%d MaxX%f",r,c,fitMaxX));
	PMTMaxSpecTDC[r][c]->Write(Form("Max TDC Spect R%d C%d",r,c));
	PMTMaxSpecTDC[r][c]->Draw("AP");
      }

      
      //Calculate target HV
      targetHV[r][c] = pmtHV[r][c]/pow(fitMaxX/targetRAU,1.0/alphas[r][c]);

      cout << "For row " << r << " and col " << c << ", target HV for next cosmic run is " << targetHV[r][c] <<  "." << endl;
      outFile << r << "  " << c << "  " << targetHV[r][c] << "  " << endl;

      //Draw return to baseline histograms
      
      cosHistFile->cd();
      retToBase[r][c]->SetTitle(Form("Return to Baseline R%d C%d (Samples)",r,c));
      retToBase[r][c]->Write(Form("Return to Baseline R%d C%d (Samples)",r,c));
      retToBase[r][c]->Draw("AP");
      
    }
  }
  
  
  
  //Post analysis reporting
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;
  cout << "Sample histograms drawn to file HistosFile.root." << endl;
  cout << "Target HV settings written to HVTargetsV3_run" << run << ".txt." << endl;

  st->Stop();

  //Sent time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}
  
