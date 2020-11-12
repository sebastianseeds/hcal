//Second attempt at plotting ave int charge vs HV setting - sseeds 11.12.20

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TSystem.h>
#include "hcal.h"

const Int_t kNrows = 12; //Number of pmt rows
const Int_t kNcols = 12; //Number of pmt columns

TFile *CosmicCalPlots = new TFile("CosmicCalPlots.root","RECREATE");
TGraph *G1 = new TGraph(); //Single graph is filled then overwritten after writing to canvas for each pmt

TCanvas *C1 = new TCanvas("C1","PMT intPulse",900,700); //First 12 x 6 grid
//TCanvas *C2 = new TCanvas("C2","PMT 72-143 intPulse",900,700); //Second 12 x 6 grid

void HCALplots(int run1 = 978, int run2 = 980, int run3 = 987, int run4 = 988, int run5 = 989){ //Can pass 5 HV runs at a time.
  
  int runs[5] = {run1, run2, run3, run4, run5}; //Configure runs for loop
  
  G1->GetXaxis()->SetTitle("HV Setting (units)"); //Same for all
  G1->GetYaxis()->SetTitle("Average Integrated Charge (units)"); //Same for all
  
  //C1->Divide(kNcols,ceil(kNrows/2)); //Divide up the canvas into two
  //C2->Divide(kNcols,floor(kNrows/2));

  //C1->Divide(kNcols,kNrows); //Just one huge canvas until this works
  
  cout << ceil(kNrows/2) << "  " << floor(kNrows/2) << endl;

  int row, col, evTot;
  double intVal;
  string line;

  vector<double> intVals = {0};
  vector<int> rowVals = {0};
  vector<int> colVals = {0};
  vector<int> HVVals = {0};
  
  for (int R=0; R<5; R++){

    int lineNumber = 0;
    ifstream inFile1(Form("intPulseRun_%d.txt",runs[R])); //Read in file for first run

    cout << "Reading in " << Form("intPulseRun_%d.txt",runs[R]) << "..." << endl;
    
    //intPulseRun must be of the format <row  col  intVal  evTot> after lines which start with #

    
    while (getline(inFile1,line))
      {
	
	if(line.at(0) == '#'){
	  continue;
	}

	lineNumber++;
	stringstream ss(line);
	ss >> row >> col >> intVal >> evTot;

	intVals.push_back(intVal);
       	rowVals.push_back(row);
       	colVals.push_back(col);
       	HVVals.push_back(runs[R]);
       	
	/*
	if(row<ceil(kNrows/2)) {
	  C1->cd(kNrows*row+col+1);
	  cout << "Draw at position " << kNrows*row+col+1 << "." << endl;
	} else {
	  C2->cd(kNrows*(row-ceil(kNrows/2))+col+1);
	  cout << "Draw at position " << kNrows*(row-ceil(kNrows/2))+col+1 << "." << endl;	  
	}
	
	G1->SetTitle(Form("IntAve Run%d Row%d Col%d",runs[R],row,col));
	cout << Form("IntAve Run%d Row%d Col%d",runs[R],row,col) << endl;
	G1->SetMinimum(0.0);
	G1->SetMaximum(5000.0);
	G1->SetPoint(R,runs[R],intVal);
	//cout << R << "  " << runs[R] << "  " << intVal << endl;
	G1->SetMarkerStyle(119);
	G1->SetMarkerSize(3);
	G1->SetMarkerColor(2);
	G1->Draw("AP same");

	CosmicCalPlots->cd();
	G1->Write(Form("Run%d,Row%d,Col%d",runs[R],row,col));
	G1->Draw("AP");
	*/
	//cout << "row = " << row << ", col = " << col << ", intVal = " << intVal << endl;
      }            
  }

  //TGraph *G[5] = new TGraph();

  int HVct = 0;
  
  for (int cL=0; cL<(kNcols*kNrows); cL++){
    C1->Clear();
    C1->cd();
    HVct = 0;
    
    for (int r=0; r<kNrows; r++){
      for (int c=0; c<kNcols; c++){
	for (int i=1; i<intVals.size()+1; i++){
	  //if (rowVals[i]*kNcols==cL && colVals[i]*kNrows==cL){
	  if (rowVals[i]*kNcols+colVals[i]==cL){
	    G1->SetPoint(HVct, HVVals[i], intVals[i]);
	    G1->SetTitle(Form("IntAve Row%d Col%d",r,c));

	    //if(HVVals[i] == 0) cout << "HV Borked at " << cL << endl;
	    //if(intVals[i] == 0) cout << "Int Borked at " << cL << endl;
	  
	    
	    //cout << "HV " << HVVals[i] << " point Set at row " << rowVals[i] << ", col " << colVals[i] << " with value " << intVals[i] << "." << endl;
	    
	    G1->SetMinimum(0.0);
	    //G1->SetMaximum(5000.0);
	    G1->SetMarkerStyle(8);
	    G1->SetMarkerSize(2);
	    //G1->Draw("AP same");

	    
	    
	    HVct++;
	    
	  }
	      
	}
      }
    }

    double rw = floor(cL/kNcols);
    double cl = cL%kNrows;
    
    CosmicCalPlots->cd();
    G1->Write(Form("Row %d,Col %d",rw,cl));
    G1->Draw("AP");
   
    cout << cL << " TGraph drawn." << endl;
    
  }

  cout << "Extraction complete." << endl;

}

