//SSeeds August 2021 - Cosmic PMT HV calibration code for use with SBS HCAL experiments. 24 rows x 12 cols HCAL modules. Relevant information can be found here: hallaweb.jlab.org/wiki/index.php/HCAL_Cosmics.

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include <TStopwatch.h>
#include <iomanip>
#include <ctime>
#include "hcal.h"

// Detector parameters and flags
const int kNrows = 24;
const int kNcols = 12;

const double kTargetRAU = 78.355; // Previous value of 61.425*(2.55/2) or kTargetRAU_prev*(attLongCable_newMeas/attLongCable_oldMeas)

// Counter to keep track of T tree entry for processing
int gCurrentEntry = -1;
TChain *T = 0;

// Augment fn's for cosmic Calibration and recording
string getDate( );

// Declare necessary histograms
TH1F *gPMTIntSpec[kNrows][kNcols]; 
TH1F *gPMTMaxSpec[kNrows][kNcols]; 
TH1F *gPMTIntSpecTDC[kNrows][kNcols]; 
TH1F *gPMTMaxSpecTDC[kNrows][kNcols]; 
TH1F *gPedSpec[kNrows][kNcols];
TH1F *gADCvChannel;
TH1F *gAmpvChannel;
TH1F *gNEVvChannel;
TH1F *gPedvChannel;

// Declare fit function
TF1 *g1;

// Create marker to track number of signals that pass cuts
int gSignalTotal = 0;

// Declare necessary arrays
double gPMTHV[kNrows][kNcols];
double gAlphas[kNrows][kNcols];

// Make machine-based date function
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);

  string yyyy = std::to_string(1900 + ltm.tm_year);
  string mm = std::to_string(1 + ltm.tm_mon);
  string dd = std::to_string(ltm.tm_mday);
  string date = mm + '/' + dd + '/' + yyyy;

  return date;
} 

void processEvent(int entry = -1, double cut = 5, int diagPlots = 0 ){
  // Check event increment and increment
  if( entry == -1 ) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }

  if( gCurrentEntry < 0 ) {
    gCurrentEntry = 0;
  }

  // Get the event from the TTree
  T->GetEntry( gCurrentEntry );
  
  int r,c;

  // Declare signal amplitude, adc, and tdc arrays for this event
  double adc[kNrows][kNcols] = {0.0};
  double tdc[kNrows][kNcols] = {0.0};
  double adc_p[kNrows][kNcols] = {0.0};
  double amp[kNrows][kNcols] = {0.0};
  double amp_p[kNrows][kNcols] = {0.0};
  double ped[kNrows][kNcols] = {0.0};
  bool saturated[kNrows][kNcols] = { {false} };

  // Process event with m data
  for( int m = 0; m < hcalt::ndata; m++ ) {
    // Define row and column
    r = hcalt::row[m];
    c = hcalt::col[m];
    if( r < 0 || c < 0 ) {
      cerr << "Error: row or col negative." << endl;
      continue;
    }
    
    if( r >= kNrows || c >= kNcols ) continue;

    // Fill adc and tdc arrays
    adc[r][c] = hcalt::a[m];
    adc_p[r][c] = hcalt::a_p[m];
    amp[r][c] = hcalt::a_amp[m];
    amp_p[r][c] = hcalt::a_amp_p[m];
    tdc[r][c] = hcalt::tdc[m];
    if( diagPlots==1 ) ped[r][c] = hcalt::ped[m];
    // Mark saturated array when amplitude meets max RAU
    if( amp[r][c] > 4094 ) {
      saturated[r][c] = true;
    }
  }

  // Vertical line test - accept only if three hits in vertical track. Reject if hit adjacent in row. Assuming pedestal sigma reasonable for all channels all events (sigma << 5 mV or 10 RAU (observed roughly 1 mV or 2 RAU))
 
  for( r = 0; r < kNrows; r++ ) {
    for( c = 0; c < kNcols; c++ ) {

      if( diagPlots==1 ) gPedSpec[r][c]->Fill( ped[r][c] );

      if( amp_p[r][c]>cut && saturated[r][c]==false ) { // Only examine channels with max > 5 sigma pedestal for hits and no saturation
	if( r==0 && c==0 ) { // Logic for top left channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 && c==11 ) { // Logic for top right channel
	  if( amp_p[r+1][c] > cut && amp_p[r+2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==0 ) { // Logic for bottom left channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 && c==11 ) { // Logic for bottom right channel
	  if( amp_p[r-1][c] > cut && amp_p[r-2][c] > cut ){
	    if( amp_p[r][c-1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==0 || r==1) { // Logic for other row 1 or 2 channels
	  if( (amp_p[r+1][c] > cut && amp_p[r+2][c] > cut) ||
	      (r==1 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else if( r==23 || r==22 ) { // Logic for other row 22 or 23 channels
	  if( (amp_p[r-1][c] > cut && amp_p[r-2][c] > cut) ||
	      (r==22 && amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
	else { // Logic for all other channels
	  if( (amp_p[r-2][c] > cut && amp_p[r-1][c] > cut) ||
	      (amp_p[r-1][c] > cut && amp_p[r+1][c] > cut) ||
	      (amp_p[r+1][c] > cut && amp_p[r+2][c] > cut) ) {
	    if( amp_p[r][c-1] < 6 && amp_p[r][c+1] < 6 ) {
	      if( tdc[r][c] != 0 ) {
		gPMTIntSpecTDC[r][c]->Fill( adc_p[r][c] );
		gPMTMaxSpecTDC[r][c]->Fill( amp_p[r][c] );
	      }
	      gPMTIntSpec[r][c]->Fill( adc_p[r][c] );
	      gPMTMaxSpec[r][c]->Fill( amp_p[r][c] );
	      gSignalTotal++;
	    }
	  }
	}
      }
    }
  }
}

// Main
int hcal_gain_match(int run = -1, int event = -1){
    
  // Initialize function with user commands
  bool diagPlots = 0;
  bool vertCut = 0;
  string date = getDate();
  double adcCut = 5;

  cout << "Enter run number for analysis." << endl;
  cin >> run;
  cout << "Cut on vertical tracks? Enter 0 for NO and 1 for YES." << endl;
  cin >> vertCut;
  cout << "Print diagnostic plots? Enter 0 for NO and 1 for YES." << endl;
  cin >> diagPlots;
  cout << "Enter ADC pulse threshold over pedestal (in mV)." << endl;
  cin >> adcCut;
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  //Declare outfile
  TFile *cosHistFile = new TFile( Form( "outFiles/gainMatchResults_run%d.root", run ), "RECREATE" );
  ofstream fitData;
  fitData.open( Form( "outFiles/HCal_%d_FitData.txt", run ) );
  fitData << "Run Number: " << run << " Desired Peak Position: " << kTargetRAU << endl;
  fitData << "Module " << " " << " Target HV " << " " << " Stat " << " " << " ErrStat " << " " << " Peak Pos " << " " << " ErrPPos " << " " << " Peak Width " << " " << " ErrPWid " << " " << " NGoodEvents " << " " <<  " " << " Flag " << std::endl;

  //Set spectrum histogram limits
  int PMTSpecBins = 100;
  //INT
  double PMTIntSpec_min = 0.0, PMTIntSpec_max = 250.0;
  //MAX
  double PMTMaxSpec_min = 0.0, PMTMaxSpec_max = 500.0;
  
  //Build spectrum histograms. Empirical limits.
  cout << "Building ADC and TDC spectrum histograms.." << endl;
  for( int r=0; r<kNrows; r++ ){
    for( int c=0; c<kNcols; c++ ){  
      gPMTIntSpec[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d", r, c ), Form( "Int ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpec[r][c]->GetXaxis()->SetTitle( "sRAU" );
      gPMTIntSpec[r][c]->GetXaxis()->CenterTitle();
      
      gPMTMaxSpec[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d", r, c ), Form( "Max ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpec[r][c]->GetXaxis()->SetTitle( "RAU" );
      gPMTMaxSpec[r][c]->GetXaxis()->CenterTitle();
      
      gPMTIntSpecTDC[r][c] = new TH1F( Form( "Int ADC Spect R%d C%d, TDC Cut", r, c ), Form( "Int ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpecTDC[r][c]->GetXaxis()->SetTitle( "sRAU" );
      gPMTIntSpecTDC[r][c]->GetXaxis()->CenterTitle();
      
      gPMTMaxSpecTDC[r][c] = new TH1F( Form( "Max ADC Spect R%d C%d, TDC Cut", r, c ), Form( "Max ADC Spect R%d C%d", r, c ), PMTSpecBins, PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpecTDC[r][c]->GetXaxis()->SetTitle( "RAU" );
      gPMTMaxSpecTDC[r][c]->GetXaxis()->CenterTitle();

      gPedSpec[r][c] = new TH1F( Form( "Pedestal Spect R%d C%d, TDC Cut", r, c ), Form( "Pedestal Spect R%d C%d", r, c ), PMTSpecBins, PMTMaxSpec_min, PMTMaxSpec_max );
      gPedSpec[r][c]->GetXaxis()->SetTitle( "RAU" );
      gPedSpec[r][c]->GetXaxis()->CenterTitle();
    }
  }

  // Build diagnostic histograms
  if( diagPlots==1 ){
    gADCvChannel = new TH1F( "ADCvChannel", "ADC vs Channel", kNcols*kNrows, 0, kNcols*kNrows-1 );
    gAmpvChannel = new TH1F( "AmpvChannel", "Amplitude vs Channel", kNcols*kNrows, 0, kNcols*kNrows-1 );
    gNEVvChannel = new TH1F( "NEVvChannel", "Number of events vs Channel", kNcols*kNrows, 0, kNcols*kNrows-1 );
    gPedvChannel = new TH1F( "PedvChannel", "Pedestal vs Channel", kNcols*kNrows, 0, kNcols*kNrows-1 );
  }
  
  // Read in data produced by analyzer in root format
  cout << "Reading replayed root file.." << endl;
  if( !T ) { 
    T = new TChain( "T" );
    T->Add( Form( "/volatile/halla/sbs/seeds/rootFiles/hcal_%d_*.root", run ) );
    T->SetBranchStatus( "*", 0 );
    T->SetBranchStatus( "sbs.hcal.*", 1 );
    T->SetBranchAddress( "sbs.hcal.a", hcalt::a );
    T->SetBranchAddress( "sbs.hcal.a_amp", hcalt::a_amp );
    T->SetBranchAddress( "sbs.hcal.a_p", hcalt::a_p );
    T->SetBranchAddress( "sbs.hcal.a_amp_p", hcalt::a_amp_p );
    T->SetBranchAddress( "sbs.hcal.tdc", hcalt::tdc );
    T->SetBranchAddress( "sbs.hcal.ped", hcalt::ped );
    T->SetBranchAddress( "sbs.hcal.adcrow", hcalt::row );
    T->SetBranchAddress( "sbs.hcal.adccol", hcalt::col );
    T->SetBranchStatus( "Ndata.sbs.hcal.adcrow", 1 );
    T->SetBranchAddress( "Ndata.sbs.hcal.adcrow", &hcalt::ndata );
    cout << "Opened up tree with nentries=" << T->GetEntries() << endl;
  }

  // Set appropriate HV and alphas for run. HV settings from HCAL wiki. Alphas from LED analysis. Must have accompanying text file, one double per line by module number. Assuming negative voltage inputs.

  ifstream HVFile( Form( "setFiles/HV_run%d.txt", run ) );

  if( !HVFile ){
    cerr << "No HV settings from run present -> setFiles/HV_run" << run << ".txt expected." << endl;
    return 0;
  }
  
  cout << "Getting HV settings for each pmt for run " << run << "." << endl;

  int n1=0;
  double d1;
  
  int rval, cval;
  string line;
    
  while( getline( HVFile, line ) ){
    if( line.at(0) == '#' ) {
      continue;
    }
    
    stringstream ss( line );
    ss >> d1;
    
    rval = floor( n1/kNcols );
    cval = n1 % kNcols;
    
    gPMTHV[rval][cval] = -d1; 
        
    n1++;
  }

  ifstream alphaFile( "setFiles/alphas.txt" );

  if( !alphaFile ){
    cerr << "No PMT alphas file present -> setFiles/alphas.txt expected." << endl;
    return 0;
  }
  
  n1=0;
  string line2;

  while( getline( alphaFile, line2 ) ){
    if( line2.at( 0 )=='#' ) {
      continue;
    }

    stringstream ss( line2 );
    ss >> d1;

    rval = floor( n1/kNcols );
    cval = n1 % kNcols;

    gAlphas[rval][cval] = d1;
    
    n1++;
  }

  gCurrentEntry = event;

  cout << "Total events to process " << T->GetEntries() << ".." << endl;

  cout << "Processing hits by event.." << endl;

  for ( int i = gCurrentEntry; i < T->GetEntries(); i++ ){ 
    processEvent( gCurrentEntry, adcCut, diagPlots );
    gCurrentEntry++;
    
    // Keep count of events processed for monitoring
    if ( gCurrentEntry%10000 == 0 ){
      cout << "Current event: " << gCurrentEntry << endl;
    }
  }

  // Fit ADC spectra with landau to extract mpv for HV calibration
  cout << "Writing spectrum histograms and calibration constants to file.." << endl;
  ofstream outFile;
  outFile.open( Form( "outFiles/HVTargets_run%d.txt", run ) ); // Create text file to hold target HVs
  time_t now = time(0); 
  char *dt = ctime(&now);
  outFile << "#Target HV settings for run " << run << " from hcal_gain_match.C on " << dt << "#" << endl;
  outFile << "#Row Col targetHV" << endl;

  // Implement two-fit scheme
  TF1 *landauInt[kNrows][kNcols] = {};
  TF1 *landauMax[kNrows][kNcols] = {};

  // Fit parameter array for built-in landau
  double parsInt[3];
  double parErrInt[3];
  double parsMax[3];
  double parErrMax[3];

  // Array to store target high voltage per module as it's calculated
  double targetHV[kNrows][kNcols];

  // Fit all spectra histograms
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){

      // Create fit functions for each module
      // INT
      landauInt[r][c] = new TF1(Form("landauInt_r%d_c%d",r,c),"landau");
      landauInt[r][c]->SetLineColor(2);
      landauInt[r][c]->SetNpx(1000);
      // MAX
      landauMax[r][c] = new TF1(Form("landauMax_r%d_c%d",r,c),"landau");;
      landauMax[r][c]->SetLineColor(4);
      landauMax[r][c]->SetNpx(1000);

      // First fit
      // INT
      int maxBinInt = gPMTIntSpec[r][c]->GetMaximumBin();
      double maxBinCenterInt = gPMTIntSpec[r][c]->GetXaxis()->GetBinCenter( maxBinInt );
      double maxCountInt = gPMTIntSpec[r][c]->GetMaximum();
      double binWidthInt = gPMTIntSpec[r][c]->GetBinWidth( maxBinInt );
      double stdDevInt = gPMTIntSpec[r][c]->GetStdDev();
      // MAX
      int maxBinMax = gPMTMaxSpec[r][c]->GetMaximumBin();
      double maxBinCenterMax = gPMTMaxSpec[r][c]->GetXaxis()->GetBinCenter( maxBinMax );
      double maxCountMax = gPMTMaxSpec[r][c]->GetMaximum();
      double binWidthMax = gPMTMaxSpec[r][c]->GetBinWidth( maxBinMax );
      double stdDevMax = gPMTMaxSpec[r][c]->GetStdDev();
      // INT
      landauInt[r][c]->SetParameters( maxCountInt, maxBinCenterInt, stdDevInt );
      landauInt[r][c]->SetRange( PMTIntSpec_min, PMTIntSpec_max );
      gPMTIntSpec[r][c]->Fit(landauInt[r][c],"No+RQ");
      landauInt[r][c]->GetParameters( parsInt );
      // MAX
      landauMax[r][c]->SetParameters( maxCountMax, maxBinCenterMax, stdDevMax );
      landauMax[r][c]->SetRange( PMTMaxSpec_min, PMTMaxSpec_max );
      gPMTMaxSpec[r][c]->Fit(landauMax[r][c],"No+RQ");
      landauMax[r][c]->GetParameters( parsMax );

  
      // Second fit with tailored range
      // INT
      int lowBinCentInt = ( maxBinInt*binWidthInt ) - ( 4.5*parsInt[2] );
      int upBinCentInt = ( maxBinInt*binWidthInt ) + ( 4.5*parsInt[2] );
      landauInt[r][c]->SetParameters( parsInt[0], parsInt[1], parsInt[2] );
      landauInt[r][c]->SetRange( lowBinCentInt, upBinCentInt );
      gPMTIntSpec[r][c]->Fit( landauInt[r][c], "+RQ" );
      landauInt[r][c]->GetParameters( parsInt );
      for ( int i=0; i<3; i++ ) parErrInt[i] = landauInt[r][c]->GetParError( i );
      // MAX
      int lowBinCentMax = ( maxBinMax*binWidthMax ) - ( 4.5*parsMax[2] );
      int upBinCentMax = ( maxBinMax*binWidthMax ) + ( 4.5*parsMax[2] );
      landauMax[r][c]->SetParameters( parsMax[0], parsMax[1], parsMax[2] );
      landauMax[r][c]->SetRange( lowBinCentMax, upBinCentMax );
      gPMTMaxSpec[r][c]->Fit( landauMax[r][c], "+RQ" );
      landauMax[r][c]->GetParameters( parsMax );
      for ( int i=0; i<3; i++ ) parErrMax[i] = landauMax[r][c]->GetParError( i ); 

      // Count up the good events
      // INT
      double goodEvInt = gPMTIntSpec[r][c]->Integral( gPMTIntSpec[r][c]->FindFixBin( lowBinCentInt ), gPMTIntSpec[r][c]->FindFixBin( upBinCentInt ), "" );
      // MAX
      double goodEvMax = gPMTMaxSpec[r][c]->Integral( gPMTMaxSpec[r][c]->FindFixBin( lowBinCentMax ), gPMTMaxSpec[r][c]->FindFixBin( upBinCentMax ), "" );

      // Calculate target HV
      targetHV[r][c] = gPMTHV[r][c]/pow(parsMax[1]/kTargetRAU,1.0/gAlphas[r][c]);
      outFile << r << " " << c << " " << targetHV[r][c] << endl; 

      // Checking on goodness of fit for max spectra
      string Flag = "Good"; // Flag to keep track of fit status
      if( gPMTMaxSpec[r][c]->GetEntries() == 0 ){
	cout << "**Module " << r << " " << c << " histo empty." << endl;
	Flag = "No_Data";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] > 60.0 ){
	cout << "**Module " << r << " " << c << " fit too wide." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Wide";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] < 1.0 ){
	cout << "**Module " << r << " " << c << " fit too narrow." << endl;
	cout << "parsMax[2] = " << parsMax[2] << endl;
	Flag = "Narrow";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parErrMax[1] > 20.0 || parErrMax[2] > 20.0 ){
	cout << "**Module " << r << " " << c << " error bar too high." << endl;
	cout << "parErrMax[1] = " << parErrMax[1] << " parErrMax[2] = " << parErrMax[2] << endl;
	Flag = "Big_error";
	for( int i=0; i<4; i++ ) { parsMax[i] = 0.0; parErrMax[i] = 0.0; }
	targetHV[r][c] = 1.0;
	goodEvMax = 0.0;
      }else if( parsMax[2] < 5.0 ){
	cout << "**Warning: Module " << r << " " << c << " appears narrow." << endl;
      }

      if( Flag != "Good" )
	cout << "Problem fit to module " << r << " " << c << " -> " << Flag << endl;

      // Write fitted histograms to file

      cout << "Writing fitted histograms to file.." << endl;
      
      cosHistFile->cd();
      // INT
      gPMTIntSpec[r][c]->SetTitle( Form( "Int Spect R%d C%d FitMean %f", r, c, parsInt[1] ) );
      gPMTIntSpec[r][c]->Write( Form( "Int Spect R%d C%d", r, c ) );
      gPMTIntSpec[r][c]->Draw( "AP" );
      // MAX
      gPMTMaxSpec[r][c]->SetTitle( Form( "Max Spect R%d C%d FitMean %f", r, c, parsMax[1] ) );
      gPMTMaxSpec[r][c]->Write( Form( "Max Spect R%d C%d", r, c ) );
      gPMTMaxSpec[r][c]->Draw( "AP" );

      if( diagPlots==1 ){
	cout << "Writing diagnostic plots to file.." << endl;
	gADCvChannel->SetBinContent( kNcols*r+c, parsInt[1] );
	gAmpvChannel->SetBinContent( kNcols*r+c, parsMax[1] );
	gPedvChannel->SetBinContent( kNcols*r+c, gPedSpec[r][c]->GetMean() );
	gNEVvChannel->SetBinContent( kNcols*r+c, gPMTIntSpec[r][c]->GetEntries() );
      }
    }
  }
    
  if( diagPlots==1 ){
    cosHistFile->cd();
    gADCvChannel->SetTitle( "Average ADC vs Channel" );
    gADCvChannel->Write( "ADCvChannel" );
    gADCvChannel->Draw( "AP" );

    gAmpvChannel->SetTitle( "Average Amplitude vs Channel" );
    gAmpvChannel->Write( "AmpvChannel" );
    gAmpvChannel->Draw( "AP" );

    gPedvChannel->SetTitle( "Average Pedestal vs Channel" );
    gPedvChannel->Write( "PedvChannel" );
    gPedvChannel->Draw( "AP" );

    gNEVvChannel->SetTitle( "Numver of Hits vs Channel" );
    gNEVvChannel->Write( "NEVvChannel" );
    gNEVvChannel->Draw( "AP" );
  }

  // Post analysis reporting
  cout << "Finished loop over run " << run << "." << endl;
  cout << "Total good signals = " << gSignalTotal << "." << endl;
  cout << "Target HV settings written to HVTargets_run" << run << ".txt." << endl;

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

  return 0;
}
  
