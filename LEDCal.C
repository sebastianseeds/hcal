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

const Int_t minSample = 0; //Initialize histogram minimum
const Int_t maxSample = 30; //Initialize histogram maximum - cannot exceed hcalt::nsamps per pmt per event
const Int_t totSample = (maxSample-minSample);

const Int_t minADC = 0;
const Int_t maxADC = 4000;

//LED variables
const int kNLED = 6; //For 5 LEDs and an off position at zero
//const int kLED[kNLED] = {16,32};
int gFoundLED;
int gEntries;

Int_t gCurrentEntry = -1;

TChain *T = 0;
TCanvas *C1 = new TCanvas();
TCanvas *C2 = new TCanvas();
TCanvas *C3 = new TCanvas();
TCanvas *C4 = new TCanvas();
TCanvas *C5 = new TCanvas();

//Augment fn's for cosmic Calibration
void goodHistoTest(TH1F*, double, int, int, int);
//void getPedestal(Int_t);
void getPedestal(TH1F*, int, int, int);
void pulseSpec(TH1F*, TF1*, double, int, int, int);

//Augmenting with parameters necessary to obtain integrated landau fit values from good pulses
Double_t nhit = 0;
TH1F *histos[kNrows][kNcols][kNLED];
TH1F *PMTIntSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,5000); //Empirical limits
TH1F *PMTMaxSpec[kNrows][kNcols][kNLED]; //=new TH1F("","",20,0,1000); //Empirical limits
TH1F *pedSpec[kNrows][kNcols][kNLED]; 
TH1F *passHisto = new TH1F("passHisto","passHisto", 19, 0, 19);
vector<double> intValues;
vector<double> maxValues;
vector<int> intValues_led;
vector<int> intValues_row;
vector<int> intValues_col;
int signalTotal = 1;
TFile *LEDHistosFile = new TFile("LEDHistosFile.root","RECREATE");  //File for checking histogram fits and integrity
int checkHistos = 0; //Iterator to save only some histograms for a check
TF1 *fFit; //Declaration for landau fit check
TF1 *fFitPed; //Declaration for gaussian fit to pedestal spectra
double pedestals[kNrows][kNcols][kNLED];
double pedSigma[kNrows][kNcols][kNLED];
int pedProcessed[kNrows][kNcols][kNLED];
Bool_t gSaturated[kNrows][kNcols][kNLED];
int badFit = 0; //Counts up per landau fit failing chi^2 test
int chisqr = 0; //Counts up per chi sqr <1 fit
int DRAWNHISTO = 0; //Counts up per histo drawn
double con; //var for landau fit. Needed for maximum value of each histo
double mu; //var for landau fit. Needed for maximum value of each histo
double sig; //var for landau fit. Needed for maximum value of each histo
double xmaxy; //var for landau fit. Needed for maximum value of each histo
double maxy; //var for landau fit. Needed for maximum value of each histo
int ledNumber; //to convert from ledbit to led on
double pedSig; //std dev of pedestals histogram
double sampSum; //var for sum of samples from pedestal bins
double sampSumN; //Normalized sum of samples from ped bins
bool P {true}; //Pass into processEvent to fill pedestal spectra
bool NP {false}; //Pass into processEvent to fill int/max spectra

int testEvents = 1500; //Limit tested events to small amount for debugging starting from event 1000
int histoInc = 0; //For sampling random histos[][][] in processEvent

int INC = 0;

TH1F* MakeHisto(int row, int col, int led, int bins){
  TH1F *h = new TH1F(TString::Format("h%02d%02d LED%d",row,col,led),TString::Format("%d-%d LED-%d",row+1,col+1,led),bins,minSample,maxSample);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Rows and Columns start at zero and go to kNrows-1 and kNcols-1
void processEvent(int entry = -1, bool pedbool = 1){
  if(entry == -1) {
    gCurrentEntry++;
    //cout << "HERE1" << endl;
  } else {
    gCurrentEntry = entry;
    //cout << "HERE2" << endl;
  }

  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
    //cout << "HERE3" << endl;
  }

  T->GetEntry(gCurrentEntry);

  if(gCurrentEntry%1000==0){
  cout << "Looping over entry " << gCurrentEntry << "." << endl;
  }
  
  //General vars by event
  int r,c,l,idx,n,sub;
  int ledbit = hcalt::ledbit;  //Attempt 1 at including ledbit info in output histograms. One ledbit per event, no loops over ledbit.
  
  ledNumber = (int) (ceil(log2(ledbit+1)));  //Verified map

  //cout << "ledNumber is " << ledNumber << "." << endl;
  
  // Clear old histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      for(l = 0; l < kNLED; l++){
	histos[r][c][l]->Reset("ICES M");
	gSaturated[r][c][l] = false;
      }
    }
  }

  //Initializing peak, adc, and tdc arrays. No ledbit dimension as redefined by event.
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

  //Get rows and columns from Tree for the event and fill histos and gsaturated
  for(int m = 0; m < hcalt::ndata; m++) {
    r = hcalt::row[m]-1;
    c = hcalt::col[m]-1;
    if(r < 0 || c < 0) {
      cout << "Why is row negative? Or col?" << endl;
      continue;
    }
    
    if(r >= kNrows || c >= kNcols)
      continue;
    idx = hcalt::samps_idx[m];
    //cout << "Samps_idx is " << hcalt::samps_idx[m] << "." << endl;
    n = hcalt::nsamps[m];  
    adc[r][c] = hcalt::a[m];
    tdc[r][c] = hcalt::tdc[m];
    bool processed = false;
    for(int s = minSample; s < maxSample && s < n; s++) {
      processed = true;
      histos[r][c][ledNumber]->SetBinContent(s+1-minSample,hcalt::samps[idx+s]);
      //cout << "R " << r << " C " << c << " LED " << ledNumber << " Samps " << hcalt::samps[idx+s] << "." << endl;
      //cout << "ledNumber is " << ledNumber << "." << endl;
      
      if(peak[r][c]<hcalt::samps[idx+s])
        peak[r][c]=hcalt::samps[idx+s];
      if(peak[r][c]>4095) {
        gSaturated[r][c][ledNumber] = true;
	//cout << "SATURATION STATION!" << endl;
      }
    }
    if(!processed) {
      cout << "Skipping empty module: " << m << endl;
      for(int s = 0;  s < totSample; s++) {
        histos[r][c][ledNumber]->SetBinContent(s+1,-404);
      }
    }
  }

  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      histos[r][c][ledNumber]->SetTitle(TString::Format("%d-%d LED-%d (ADC=%g,TDC=%g)",r+1,c+1,ledNumber,adc[r][c],tdc[r][c]));     
      if(gSaturated[r][c][ledNumber])
        histos[r][c][ledNumber]->SetLineColor(kRed+1);
      else
        histos[r][c][ledNumber]->SetLineColor(kBlue+1);
      if(tdc[r][c]!=0){
        histos[r][c][ledNumber]->SetLineColor(kGreen+1);
      }
      
      //if(!gSaturated[r][c][ledNumber]&&!pedbool){
      if(!pedbool){
	goodHistoTest(histos[r][c][ledNumber],tdc[r][c],r,c,ledNumber);
	checkHistos++;
	//}else if(!gSaturated[r][c][ledNumber]&&pedbool){
      }else if(pedbool){

	//C2->cd();
	//histos[r][c][ledNumber]->Draw();
	
	getPedestal(histos[r][c][ledNumber],r,c,ledNumber);
      }else{
	//cout << "A CELEBRATION!" << endl;
	continue;
      }
    }
  }
}

void getPedestal(TH1F *pedHisto, int row, int col, int ledNumber){

  sampSum = 0.0;
  sampSumN = 0.0;

  //C3->cd();
  //pedHisto->Draw();
  
  for (int i=minSample; i<(minSample+5); i++){  //Obtain pedestal from first five bins only where no pulse exists. Will expect pulse from every pmt on every event.
    sampSum+=pedHisto->GetBinContent(i+1); //Data starts at bin 1, not zero (underflow)
  }
  sampSumN = sampSum/(minSample+5); //Normalize by number of bins used per event
  
  pedSpec[row][col][ledNumber]->Fill(sampSumN); //Fill histogram with pedestal value from current event

  //cout << sampSumN << endl;

  //C1->cd();
  //pedSpec[row][col][ledNumber]->Draw();
  
}

//Title a bit of a misnomer here. Simply subtract pedestal and fill histograms. Will need landau to compute maxvalue per histogram.
void goodHistoTest(TH1F *testHisto, double tdcVal, int row, int col, int ledNumber){

  //Get pedestal value and std dev for the event 
  double pedVal = pedestals[row][col][ledNumber];
  double pedSig = pedSigma[row][col][ledNumber];
  
  //Subtract pedestal value from each bin
  for(int b=1 ; b<=testHisto->GetNbinsX() ; b++) {
    testHisto->SetBinContent(b,testHisto->GetBinContent(b)-pedVal);
    testHisto->SetBinError(b,pedSig);
    //if(testHisto->GetBinContent(b) < 0) testHisto->SetBinContent(b,0.0); //Remove negative values after pedestal subtraction
    //if(testHisto->GetBinContent(b) < -100) cout << "LOW BIN CONTENT " << testHisto->GetBinContent(b) << "." << endl; //Report if ped sub is failing
  }
  
  //testHisto->Fit("landau", "Q", "", 1, 30);  //Fit histogram with landau, limits window minus overflow on 20
  testHisto->Fit("landau", "QR");  //Fit histogram with landau, limits window minus overflow on 20
  fFit=testHisto->GetFunction("landau"); //Set error of bin content
  
  int dof = 27;  //30 bins -3 for fit per histo (30 bins - 3 fit parameters (con, mu, sigma))
  double rChiSqr = (fFit->GetChisquare())/dof;
  
  pulseSpec(testHisto, fFit, tdcVal, row, col, ledNumber); //Call integratePulse() for histograms that pass all tests
  
}

//integrate good signals, get max from good signals, return spectrums and averages
void pulseSpec(TH1F *intHisto, TF1 *fFit, double tdcVal, int row, int col, int ledNumber){
  
  intValues.push_back(intHisto->Integral(1,30));  //Just getting total bin content instead of integrating landau
  
  //Using TMath and obtaining parameters from fit function - similar to method in analyze_cosmics.C l:377. May wish to save these params to root tree output.
  con = fFit->GetParameter(0);
  mu = fFit->GetParameter(1);
  sig = fFit->GetParameter(2);
  xmaxy = fFit->GetMaximumX(mu*0.5,mu*1.5);
  maxy  = TMath::Landau(xmaxy,mu,sig)*con; //Calculate maximum value in histogram from maximum of landau fit

  //Fill vectors
  maxValues.push_back(maxy);
  intValues_row.push_back(row);
  intValues_col.push_back(col);
  intValues_led.push_back(ledNumber);

  //Fill spectra histograms
  //PMTIntSpec[row][col]->Fill(fFit->Integral(0,1000)); //Integral from landau fit
  PMTIntSpec[row][col][ledNumber]->Fill(intHisto->Integral(1,30)); //Integral from total bin content
  PMTMaxSpec[row][col][ledNumber]->Fill(maxy);
  
  if (signalTotal % 10000 == 0){
    cout << "Writing a reference good fit histogram to file.." << endl;
    cout << "Maximum y value for this histogram = " << maxy << "." << endl;
    cout << "Pedestal for same histogram = " << pedestals[row][col][ledNumber] << "." << endl;
    cout << "TDC value for same histogram = " << tdcVal << "." << endl;
    cout << "LED number for same histogram = " << ledNumber << "." << endl << endl;
    LEDHistosFile->cd();
    intHisto->SetTitle(Form("Good Fit R%d C%d LED%d",row,col,ledNumber));
    intHisto->Write(Form("sampleHisto%d_r%d_c%d_led%d",signalTotal,row,col,ledNumber));
    intHisto->Draw("AP");
    DRAWNHISTO++;
  }
  signalTotal++;
}

Int_t LEDCal(int run = 1205, int event = 1000)  //Start after LEDs warm up ~1000 events
{
  //Build spectra histograms. Empirical limits.
  cout << "Building spectrum histograms.." << endl;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=0; l<kNLED; l++){
	PMTIntSpec[r][c][l] = new TH1F(Form("Int ADC Spect R%d C%d LED%d",r,c,l),Form("Int ADC Spect R%d C%d LED%d",r,c,l),30,0,0); 
	PMTMaxSpec[r][c][l] = new TH1F(Form("Max ADC Spect R%d C%d LED%d",r,c,l),Form("Max ADC Spect R%d C%d LED%d",r,c,l),30,0,0);
	pedSpec[r][c][l] = new TH1F(Form("Pedestal Spect R%d C%d LED%d",r,c,l),Form("Pedestal Spect R%d C%d LED%d",r,c,l),30,0,0);
      }
    }
  }

  //cout << "1" << endl;
  
  if(!T) { 
    T = new TChain("T");
    T->Add(TString::Format("fadc_f1tdc_%d.root",run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
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
	for(int l = 0; l < kNLED; l++) {
	  histos[r][c][l] = MakeHisto(r,c,l,totSample);
	  gSaturated[r][c][l] = false;
	}
      }
    }
  }

  //cout << "2" << endl;
  
  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << endl;

  //Obtain pedestal histograms from all entries for each led at each ledbit using first three bins
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
  //for (int i = gCurrentEntry; i < testEvents; i++){ 
    processEvent(gCurrentEntry, P);
    gCurrentEntry++;

    //cout << "INC " << INC << endl;
    //INC++;
    
    //Keep count of events processed for pedestal
    if (gCurrentEntry%10000 == 0){
      cout << "Pedestal -> processing event: " << gCurrentEntry << endl;
    }
  }

  //C4->cd();
  //pedSpec[0][0][1]->Draw();

  //cout << "3" << endl;
  
  //Get mean and std dev from each pedestal spectra histo and put in array
  for(int r = 0; r < kNrows; r++) {
    for(int c = 0; c < kNcols; c++) {
      for(int l = 0; l < kNLED; l++){
	if(pedSpec[r][c][l]->GetEntries()==0){
	  continue;
	}else{
	  pedestals[r][c][l]=pedSpec[r][c][l]->GetMean(); //Get the mean
	  pedSigma[r][c][l]=pedSpec[r][c][l]->GetRMS(); //Get std dev
	  
	  cout << "Pedestals: row col LED mean sigma -> " << r << " " << c << " " << l << " " << pedestals[r][c][l] << " " << pedSigma[r][c][l] << endl;
	}
      }
    }
  }

  gCurrentEntry = event; //Reset event count for Int/Max calculation
  
  //Obtain average max/int value over histograms for all entries by led at each ledbit
  for (int i = gCurrentEntry; i < T->GetEntries(); i++){ 
  //for (int i = gCurrentEntry; i < testEvents; i++){ 
    processEvent(gCurrentEntry, NP);
    gCurrentEntry++;
    
    //Keep count of events processed for max/int avg values
    if (gCurrentEntry%10000 == 0){
      cout << "Max/Int Value -> processing event: " << gCurrentEntry << endl;
    }
  }
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << signalTotal << "." << endl;
  cout << "Obtaining averages for max value and integrated pulses.." << endl;
  
  //Now take in three vectors and get average integrated and max value pulse per pmt (unique r c value)
  double pmtIntTotal[kNrows*kNcols]={0};
  double pmtEvTotal[kNrows*kNcols]={0};
  double pmtMaxTotal[kNrows*kNcols]={0};

  for(int ir=0; ir<kNrows; ir++){
    
    for(int ic=0; ic<kNcols; ic++){
      
      for(int iev=0; iev<intValues.size(); iev++){
	
	int evRow = intValues_row[iev];
	int evCol = intValues_col[iev];

	if(evRow == ir && evCol == ic){
	  
	  pmtIntTotal[ir*kNrows+ic] += intValues[iev]; //build value for int tot, one val per pmt
	  pmtMaxTotal[ir*kNrows+ic] += maxValues[iev];
	  pmtEvTotal[ir*kNrows+ic] += 1.0;
	  
	}
      }
    }
  }

  cout << "Writing spectrum histograms to file.." << endl;

  TFile *ledAggFile = new TFile(Form("ledIntFile_run%d.root",run),"RECREATE");
  
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      for(int l=0; l<kNLED; l++){
	LEDHistosFile->cd();
	PMTIntSpec[r][c][l]->Write(Form("Int ADC Spect R%d C%d LED%d",r,c,l));
	PMTIntSpec[r][c][l]->Draw("AP");

	ledAggFile->cd();
	PMTIntSpec[r][c][l]->Write(Form("Int ADC Spect R%d C%d LED%d",r,c,l));
	PMTIntSpec[r][c][l]->Draw("AP");
	
	LEDHistosFile->cd();
	PMTMaxSpec[r][c][l]->Write(Form("Max ADC Spect R%d C%d LED%d",r,c,l));
	PMTMaxSpec[r][c][l]->Draw("AP");

	LEDHistosFile->cd();
	pedSpec[r][c][l]->Write(Form("Pedestal Spect R%d C%d LED%d",r,c,l));
	pedSpec[r][c][l]->Draw("AP");

	ledAggFile->cd();
	pedSpec[r][c][l]->Write(Form("Pedestal Spect R%d C%d LED%d",r,c,l));
	pedSpec[r][c][l]->Draw("AP");
      }
    }
  }  

  
  ofstream outFile;

  outFile.open(Form("AvePulseRun_%i.txt",run));

  outFile << "#Averaged int pulse and pulse max data for run " << run << "." << endl;
  outFile << "#Total number of pulse events is " << intValues.size() << "." << endl;
  outFile << "#Row  Column  pmtIntValue  pmtMaxValue  EvTotal" << endl;
  
  for(int e=0; e<(kNrows*kNcols); e++){
    if(pmtEvTotal[e] == 0){
      pmtIntTotal[e]=-1.0; //Make empty bins obvious
      pmtMaxTotal[e]=-1.0; //Make empty bins obvious
    } else {
    pmtIntTotal[e] = pmtIntTotal[e]/pmtEvTotal[e];
    pmtMaxTotal[e] = pmtMaxTotal[e]/pmtEvTotal[e];
    }
    cout << "PMT row = " << floor(e/kNrows) << ", col = " << e % kNcols << ": Int Pulse Avg = " << pmtIntTotal[e] << " and Max Pulse Avg = " << pmtMaxTotal[e] << " with " << pmtEvTotal[e] << " events averaged over." << endl;
    outFile << floor(e/kNrows) << "  " << e % kNcols << "  " << pmtIntTotal[e] << "  " << pmtMaxTotal[e] << "  " << pmtEvTotal[e] << endl;
  }

  cout << "Integrated averages, maximum averages, and total events per PMT written to file AvePulseRun_" << run << ".txt." << endl;

  cout << "ADC spectrum histograms written to file LEDHistosFile.root" << endl;

  cout << "Checked histograms = " << checkHistos << " written to file LEDHistosFile.root." << endl;
  
  cout << "Total bad fits = " << badFit << " rejected." << endl;

  cout << "Total fitted histograms drawn to file LEDHistosFile.root = " << DRAWNHISTO << "." << endl;
  
  return 0;
}

