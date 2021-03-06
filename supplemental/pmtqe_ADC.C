#include "hcal.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

TChain *T = 0;

TCanvas *tmpCanv = 0;

const Int_t kSkip = 10;
const Int_t kLEDCount = 1000;
//const Int_t kLEDCount = 500.0;
Bool_t kPedSubtractQuad = false;

const Int_t kNLED = 63;
//const Int_t kLED[kNLED] = { 8, 16, 32 };
const Int_t kLED[kNLED] = { 16, 32 };
Int_t gFoundLED;
Int_t gEntries;

const Double_t kMinNPE[kNLED] = {  0.0,  0.0,  /*0.0*/ };
const Double_t kMaxNPE[kNLED] = { 10000.0, 10000.0, /*10000.0*/ };
const Int_t kMinLEDCount[kNLED] = { 150, 600 } ;
Double_t kMinADC[7] = {  0.0, 0.0, 0.0, 0.0, 0.0, 12000., 20000.  /*0.0*/ };
Double_t kMaxADC[7] = {  0.0, 0.0, 0.0, 0.0, 0.0, 30000., 30000. /*10000.0*/ };

const Int_t kPedMinSamp = 0;
const Int_t kPedMaxSamp = 20;
const Bool_t kLogY = true;

const Int_t kNcols = 6;
const Int_t kNrows = 6;
const Int_t kColOffset = 0;
const Int_t kRowOffset = 0;
TH1F *hNPE[kNrows][kNcols][kNLED] = { 0 };
TGraphErrors *gNPE[kNrows][kNcols][kNLED] = { 0 };
TGraphErrors *gADC[kNrows][kNcols][kNLED] = { 0 };
//std::vector<Double_t> vNPE[kNrows][kNcols][kNLED][4];
Int_t gTotalEntries[kNrows][kNcols][kNLED];

struct Value_t {
  Double_t avg; ///< Value (avg if taken from a number of values)
  Double_t sig; ///< Standard deviation of this average
  Double_t err; ///< Error (if any) associated to this value
  Int_t n; ///< Entries that went into computing this average
  Value_t() : avg(0.), sig(-1),err(-1),n(0) {}
  Value_t(Double_t p_avg, Double_t p_sig = -1., Double_t p_err = -1.,
      Int_t p_n = 1) {
    avg = p_avg;
    sig = p_sig;
    err = p_err;
    n   = p_n;
  }

  Value_t& operator=(Value_t const&v) {
    avg = v.avg;
    sig = v.sig;
    err = v.err;
    n   = v.n;
    return *this;
  }
};

struct GraphYData_t {
  std::vector<Double_t> xv;
  std::vector<Double_t> xe;
  std::vector<Double_t> yv;
  std::vector<Double_t> ye;
  size_t size() { return yv.size(); };
  void push_back(Double_t y, Double_t e = 0.0) {
    xv.push_back(xv.size());
    yv.push_back(y);
    ye.push_back(e);
    xe.push_back(0.0);
    //std::cerr << "Pushing back: " << xv.size()-1 << " " << 0 << " " << y << " " << e;
  }
  void push_back(Value_t y) {
    push_back(y.avg,y.err);
  }
};

struct ResLED_t {
  Int_t row;
  Int_t col;
  Int_t led;
  Value_t npe;
  ResLED_t(Int_t r = -1, Int_t c = -1, Int_t l = -1) : row(r), col(c), led(l) {};
  ResLED_t(Int_t r, Int_t c, Int_t l,Value_t v) :
    row(r), col(c), led(l), npe(v) {};
};

std::vector<Value_t> vNPE[kNrows][kNcols][kNLED];
std::vector<Value_t> vADC[kNrows][kNcols][kNLED];
std::vector<ResLED_t> results[kNrows][kNcols];
std::vector<ResLED_t> resultsADC[kNrows][kNcols];
GraphYData_t toGraphNPE[kNrows][kNcols][kNLED];
GraphYData_t toGraphADC[kNrows][kNcols][kNLED];


Value_t ComputeStats(std::vector<Value_t> &vec, Bool_t weighted = kFALSE)
{
  Value_t res;
  res.avg = res.err = res.sig = 0.0;
  res.n = vec.size();
  if(res.n<=0) {
    res.sig = res.err = -1;
    return res;
  } else if ( res.n == 1) {
    res.avg = vec[0].avg;
    res.sig = 0;
    res.err = TMath::Abs(res.avg);
    return res;
  }

  for(Int_t i = 0; i < res.n; i++) {
    res.avg += vec[i].avg;
  }
  res.avg /= Double_t(res.n);
  Double_t sig2 = 0.0;
  Double_t sigi = 0.0;
  for(Int_t i = 0; i < res.n; i++) {
    sigi  = vec[i].avg-res.avg;
    sig2 += sigi*sigi;
    //sig2 += TMath::Power(vec[i].avg-res.avg,2.0);
  }
  res.sig = TMath::Sqrt(sig2/Double_t(res.n-1));
  res.err = res.sig/TMath::Sqrt(res.n);

  return res;
}

Value_t ComputeStatsADC(std::vector<Value_t> &vec, Bool_t weighted = kFALSE)
{
  Value_t res;
  res.avg = res.err = res.sig = 0.0;
  res.n = vec.size();
  if(res.n<=0) {
    res.sig = res.err = -1;
    return res;
  } else if ( res.n == 1) {
    res.avg = vec[0].avg;
    res.sig = 0;
    res.err = TMath::Abs(res.avg);
    return res;
  }

  for(Int_t i = 0; i < res.n; i++) {
    res.avg += vec[i].avg;
  }
  res.avg /= Double_t(res.n);
  for(Int_t j = 0; j < res.n; j++) {
    res.sig += vec[j].sig;
    //printf("%d - %f; sum = %f;\n", j, vec[j].sig, res.sig);
  }
  res.sig /= Double_t(res.n);
  //printf("divided by: %d = %f\n", res.n, res.sig);
  return res;
}


Value_t ComputeNPE(Value_t adc)
{
  Value_t npe(0.0);
  if(adc.sig > 0) {
    npe.avg = TMath::Power(adc.avg/adc.sig,2.0);
    npe.n   = adc.n;
    // npe.err = 1./TMath::Sqrt(npe.n);
    /* npe.err = TMath::Sqrt(
                   (TMath::Power(2*adc.avg*adc.sig,2.0))
                 + (TMath::Power(adc.sig/TMath::Sqrt(adc.n),2.0)) ) */
    npe.err = 2*TMath::Sqrt( (npe.avg+TMath::Power(npe.avg, 2.0)) / npe.n );
    npe.sig = adc.sig;
  }
  return npe;
}

Int_t getLEDIdx(Int_t led)
{
  for(Int_t l; l < kNLED; l++) {
    if(kLED[l]==led)
     return l;
  }
  return -1;
}
//std::set<Int_t> minadc_p[kNrows][kNcols][kNLED];
//std::set<Int_t> maxadc_p[kNrows][kNcols][kNLED];
//std::set<Int_t> minadc_r[kNrows][kNcols][kNLED];
//std::set<Int_t> maxadc_r[kNrows][kNcols][kNLED];
Double_t minadc_p[kNrows][kNcols][kNLED];
Double_t maxadc_p[kNrows][kNcols][kNLED];
Double_t minadc_r[kNrows][kNcols][kNLED];
Double_t maxadc_r[kNrows][kNcols][kNLED];
Bool_t gSaturated[kNrows][kNcols][kNLED];
Bool_t gFoundBits[kNrows][kNcols][kNLED];

TH1F *hADC[kNrows][kNcols][kNLED];

void ProcessNextLED(Int_t &entry)
{
  T->GetEntry(entry);
  Int_t led = hcalt::ledbit;
  if(led<0||led>=64) {
    entry++;
    return;
  }
  //Int_t iled = getLEDIdx(led);
  Int_t iled = led;
  Bool_t isGood = true;
  std::vector<Value_t> a_r[kNrows][kNcols];
  std::vector<Value_t> a_p[kNrows][kNcols];
  std::vector<Value_t> ped[kNrows][kNcols];


  // Common use variables
  Int_t r,c;
  Int_t nsamps;
  Int_t idx,is;
  std::vector<Value_t> temp_ped;
  temp_ped.resize(kPedMaxSamp-kPedMinSamp);
  Value_t lped,la_r,la_p;
  Int_t lidx = 0;
  Int_t ecount = 0;
  while(entry < gEntries) {
    T->GetEntry(entry++);
    ecount++;
    if(hcalt::ledbit != led) {
      //std::cerr << "Entry: " << entry-1 << " "
      //  << " (" << led << " != " << hcalt::ledbit << ")" << std::endl;
      entry--; // Return it to the previous entry
      break;
    } else if ( ecount < kSkip ) {
      continue;
    } else if ( ecount >= kLEDCount ) {
      break;
    }
    for(Int_t m = 0; m < hcalt::ndata; m++) {
      r = hcalt::row[m]-kRowOffset-1;
      c = hcalt::col[m]-kColOffset-1;
      nsamps = hcalt::nsamps[m];
      idx = hcalt::samps_idx[m];
      if(r>=0 && c>=0 && nsamps>kPedMaxSamp && r<kNrows && c < kNcols) {
        la_r = Value_t(hcalt::a[m]);
	//	printf(">>>%d %d %d %f\n",r,c,led,la_r.avg);
        if(!gFoundBits[r][c][iled]) {
          gFoundBits[r][c][iled] = true;
          a_r[r][c].reserve(kLEDCount);
          a_p[r][c].reserve(kLEDCount);
          ped[r][c].reserve(kLEDCount);
        }
        is = 0;
        for(Int_t s = kPedMinSamp; s < kPedMaxSamp; s++,is++) {
          temp_ped[is] = Value_t(hcalt::samps[idx+s]);
          if(hcalt::samps[idx+s]>4095)
            gSaturated[r][c][iled] = true;
        }
        // Check for gSaturated ADC
        for(Int_t s = kPedMaxSamp; s < nsamps; s++) {
          if(hcalt::samps[idx+s]>4095)
            gSaturated[r][c][iled] = true;
        }
        if(!gSaturated[r][c][iled]) {
          lped = ComputeStats(temp_ped,r==0&&c==0&&led==32);
          // Since fA in SBSCalorimeter is now the *average* ADC value of all samples, don't subtract out
          // more than 1 sample from this value
	  lped.avg *= Double_t(nsamps);
          if(lped.sig > 0) {
            la_p  = Value_t(la_r.avg-lped.avg);
            a_r[r][c].push_back(la_r);
            a_p[r][c].push_back(la_p);
            ped[r][c].push_back(lped);
            hADC[r][c][iled]->Fill(la_p.avg,1.0);
          }
        }
      } // end if good r,c, nsamps
    } // end for loop over modules
  } // end while loop over entries

  // Finally, compute the NPE for all the found Modules for this LED setting
  Value_t avg_p;
  Value_t npe;
  Value_t tst_adc;
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(gFoundBits[r][c][iled]) {
        tst_adc = ComputeStats(a_p[r][c]);
        //printf("[%2d,%2d]\tADCavg = %f\tADCsig = %f\n", r+1,c+1,tst_adc.avg, tst_adc.sig);
        // tst_adc.err = tst_adc.sig;
        npe = ComputeNPE(ComputeStats(a_p[r][c]));
        vADC[r][c][iled].push_back(tst_adc);
        vNPE[r][c][iled].push_back(npe);
        toGraphNPE[r][c][iled].push_back(npe);
        toGraphADC[r][c][iled].push_back(tst_adc);
	//printf("r,c=%d %d patt=%2d npe=%7.1f adc=%8.1f adc_sig=%5.1f\n",r,c,iled,npe.avg,tst_adc.avg,tst_adc.sig);
        if(r==0&&c==0&&iled==32) {
          Value_t tmp2 = ComputeStats(a_p[r][c]);
          /*for(size_t k = 0; k < a_p[r][c].size(); k++) {
            std::cerr << "k: " << k
              << ", a_r: " << a_r[r][c][k].avg << " +/- " << a_r[r][c][k].err
              << ", ped: " << ped[r][c][k].avg << " +/- " << ped[r][c][k].err
              << ", a_p: " << a_p[r][c][k].avg << " +/- " << a_p[r][c][k].err
              << std::endl;
          }*/
          /*
          std::cerr << "avg: " << tmp2.avg
            << ", sig: " << tmp2.sig << ", err: " << tmp2.err << std::endl;
          std::cerr << std::endl << std::endl << std::endl
            << "DATA DATA DATA vvvvv" << std::endl << std::endl << std::endl;
          for(size_t j = 0; j < a_p[r][c].size(); j++) {
            std::cerr << "j: " << j << "\t" << a_p[r][c][j].avg << "\t" << a_r[r][c][j].avg << "\t" << ped[r][c][j].avg << std::endl;
          }
          std::cerr << std::endl << std::endl << std::endl
            << "DATA DATA DATA ^^^^^" << std::endl << std::endl << std::endl;
            */
          a_r[r][c].clear();
          ped[r][c].clear();
          a_p[r][c].clear();
        }
      }
    }
  }
}

/*
void processLED2(Int_t iled, Int_t led, Int_t &entry)
{
  Int_t maxcount = (kLEDCount-kSkip);
  Bool_t isGood=true;
  Double_t ped[kNrows][kNcols][kNLED];
  Double_t pedEntries[kNrows][kNcols][kNLED];
  //std::vector<Double_t> vals_r[kNrows][kNcols];
  std::vector<Double_t> vals_p[kNrows][kNcols];
  Double_t lavg[kNrows][kNcols];
  Int_t r,c;
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      lavg[r][c] = 0.0;
    }
  }
  Double_t adc_r,adc_p;
  Int_t ecount = 0;
  while(isGood && entry < gEntries) {
    if(ecount<kSkip) {
      ecount++;
      entry++;
      continue;
    }
    T->GetEntry(entry++);
    if(hcalt::ledbit != led) {
      isGood = false;
      continue;
    }
    if(!gFoundBits[iled]) {
      gFoundBits[iled]=1;
      gFoundLED++;
    }

    for(Int_t m = 0; m < hcalt::ndata; m++) {
      r = hcalt::row[m]-1;
      c = hcalt::col[m]-1;
      adc_r = hcalt::a[m];
      if(!hNPE[r][c][iled]) {
        hNPE[r][c][iled] = new TH1F(TString::Format("hLEDNPE%02d%02d%02d",r+1,c+1,led),
         TString::Format("NPE [%02d,%02d] LED %02d",r+1,c+1,led),
	  100, kMinNPE[iled], kMaxNPE[iled]);
        gSaturated[r][c][iled] = 0;
        //vals_r[r][c].reserve(kLEDCount);
        vals_p[r][c].reserve(kLEDCount);
        lavg[r][c] = 0.0;
        gTotalEntries[r][c][iled] = 0;
      }
      Int_t nsamps = hcalt::nsamps[m];
      Int_t idx = hcalt::samps_idx[m];
      if(nsamps>kPedMaxSamp) {
        ped[r][c][iled] = 0;
        pedEntries[r][c][iled] = 0;
        for(Int_t s = kPedMinSamp; s < kPedMaxSamp; s++) {
          ped[r][c][iled] += hcalt::samps[idx+s];
          pedEntries[r][c][iled]++;
          if(hcalt::samps[idx+s]>4095)
            gSaturated[r][c][iled] = true;
        }
        // Check for gSaturated ADC
        for(Int_t s = kPedMaxSamp; s < nsamps; s++) {
          if(hcalt::samps[idx+s]>4095)
            gSaturated[r][c][iled] = true;
        }
        if(pedEntries[r][c][iled]>0) {
          Double_t tmp2 = ped[r][c][iled];
          ped[r][c][iled] /= Double_t(pedEntries[r][c][iled]);
          adc_p = adc_r - ped[r][c][iled]*Double_t(nsamps);
          vals_p[r][c].push_back(adc_p);
          Double_t tmp3 = lavg[r][c];
          lavg[r][c] += adc_p;
          if(false&&(r!=4&&r!=0)&&(entry==4409||tmp3!=tmp3 || lavg[r][c] != lavg[r][c] || adc_p < 1||lavg[r][c]<1)) {
              std::cerr << "Strange: " << r+1 << "-" << c+1
                << " " << tmp2
                << " " << Double_t(pedEntries[r][c][iled])
                << " " << ped[r][c][iled]
                << " " << adc_r
                << " " << adc_p
                << " " << tmp3
                << " " << lavg[r][c]
                << " " << entry-1
                << " " << led
                << " " << hcalt::ledcount
                << std::endl;
                exit(-1);
          }
          if(adc_r<10) {
            std::cerr << "Entry: " << entry << ", adc_r: " << adc_r << std::endl;
          }
        }
      }
    }
  }
  Double_t sig2;
  Double_t npe;
  Double_t sig;
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      if(vals_p[r][c].size()>kMinLEDCount[iled]) {
        Double_t tmp = lavg[r][c];
        lavg[r][c] /= Double_t(vals_p[r][c].size());
        gTotalEntries[r][c][iled] += Int_t(vals_p[r][c].size());
        sig2 = 0;
        for(size_t e = 0; e < vals_p[r][c].size(); e++) {
          sig2 += TMath::Power(vals_p[r][c][e]-lavg[r][c],2.);
        }
        sig = TMath::Sqrt(sig2/(Double_t(vals_p[r][c].size())-1.));
        npe = TMath::Power(lavg[r][c]/sig,2.);
        if(npe != npe) {
          std::cerr << "NAN: " << r+1 << "-" << c+1
            << " " << tmp
            << " " << Double_t(vals_p[r][c].size())
            << " " << lavg[r][c]
            << " " << sig2
            << " " << sig
            << " " << ped[r][c][iled]
            << std::endl;
        }
        hNPE[r][c][iled]->Fill(npe);
        vNPE[r][c][iled][0].push_back(npe);
        vNPE[r][c][iled][1].push_back(lavg[r][c]);
        vNPE[r][c][iled][2].push_back(sig);
        vNPE[r][c][iled][3].push_back(vNPE[r][c][iled][0].size());
      }
    }
  }
}
*/

void pmtqe_ADC(Int_t run=301, Bool_t saveCanvas = kFALSE)
{
  Bool_t readLimsFromFile=true;
  gStyle->SetLabelSize(0.05,"XY");
  if(!T) {
    T = new TChain("T");
    T->Add(TString::Format("%s/fadcAERO_%d.root",getenv("HCAL_ROOTFILES"),run));
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("sbs.hcal.*",1);
    T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
    T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
    T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);
    T->SetBranchAddress("sbs.hcal.a",hcalt::a);
    T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);
    T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx);
    T->SetBranchAddress("sbs.hcal.row",hcalt::row);
    T->SetBranchAddress("sbs.hcal.col",hcalt::col);
    T->SetBranchStatus("Ndata.sbs.hcal.row",1);
    T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);
  }
  gStyle->SetLabelSize(0.1,"XY");
  gStyle->SetNdivisions(5);
  // Initialize things...
  for(Int_t r = 0; r < kNrows; r++) {
    for(Int_t c = 0; c < kNcols; c++) {
      for(Int_t iled = 0; iled < kNLED; iled++) {
        gSaturated[r][c][iled] = gFoundBits[r][c][iled] = false;
        int ledbit = (int) (log2(iled+1)+1);
        hADC[r][c][iled] = new TH1F(TString::Format("hADC%d%d%d",r+1,c+1,iled),
          TString::Format("ADC [%d,%d]",r+1,c+1),100,kMinADC[ledbit],kMaxADC[ledbit]);
      }
    }
  }

  gEntries = T->GetEntries();
  Int_t led,r,c,ledcount,iled;
  Int_t foundLED=0;
  Double_t adc_r = 0;
  Double_t adc_p = 0;
  //for(Int_t entry  = 0; entry < entries; entry++) {
  Int_t entry = 6000;
  while(entry<gEntries) {
    /* if(entry < 100) { 
      std::cout << "Processing entry: " << entry << std::endl;
      } */
    if(entry%1000 == 0) {
      printf("%d\n", entry);
    }
    //T->GetEntry(entry);
    ProcessNextLED(entry);
  }
  
  std::vector<TCanvas*> canv;
  Double_t ravg,rsig,rent;
  ResLED_t res;
  ResLED_t resADC;
  ofstream outlim;
  ofstream outsummary(TString::Format("amplitudes/%d.dat",run),std::ios::out);;
  if(!readLimsFromFile) {
    outlim.open(TString::Format("amplitudes/limits_%d.dat",run),std::ios::out);
  }
  for(iled = 0; iled < kNLED; iled++) {
    //led=kLED[iled];
    led = iled;
    ravg = rsig = rent = 0;
    TCanvas *canvas =  0;
    TCanvas *canvasADC =  0;
    for(r = 0; r < kNrows; r++) {
      for(c = 0; c < kNcols; c++) {
        if(gFoundBits[r][c][iled]) {
          if(!canvasADC && iled==16) {
            std::cerr << "Creating canvas " << TString::Format("canvas%02d",led)
              << std::endl;
            
            canvas = new TCanvas(TString::Format("canvas%02d",led),
                TString::Format("canvas%02d",led),400*kNcols,400*kNrows);
            canvas->Divide(kNrows,kNcols);
            
            canvasADC = new TCanvas(TString::Format("canvasADC%02d",led),
                TString::Format("canvasADC%02d",led),400*kNcols,400*kNrows);
            canvasADC->Divide(kNrows,kNcols);

            //canv.push_back(canvas);
          }
          if(canvas){
          canvas->cd(c+r*kNcols+1);
          //hNPE[r][c][iled]->Draw();
          //gPad->SetLogy(kLogY);
          gNPE[r][c][iled] = new TGraphErrors(toGraphNPE[r][c][iled].size(),
            toGraphNPE[r][c][iled].xv.data(),toGraphNPE[r][c][iled].yv.data(),
            toGraphNPE[r][c][iled].xe.data(),toGraphNPE[r][c][iled].ye.data());
          //TString tmp = TString::Format("[%2d,%2d,%2d]",r+1,c+1,led);
          gNPE[r][c][iled]->SetMarkerStyle(20);
          gNPE[r][c][iled]->SetTitle(TString::Format("%02d-%02d LED %02d",r+1+kRowOffset,c+1+kColOffset,led));
          gNPE[r][c][iled]->Draw("AP");
          gPad->SetLogy(false);
          }

          // Also graph the ADC (where the error is the standard deviation so we can see them both visually
          if(canvasADC) {
            canvasADC->cd(c+r*kNcols+1);
            hADC[r][c][iled]->Draw();
            Double_t hmean = hADC[r][c][iled]->GetMean();
            printf("Now fitting %d, %d...\n", r+1, c+1);
            hADC[r][c][iled]->Fit("gaus","","",hmean*0.75,hmean*1.25);
          }
          /*
          gADC[r][c][iled] = new TGraphErrors(toGraphADC[r][c][iled].size(),
            toGraphADC[r][c][iled].xv.data(),toGraphADC[r][c][iled].yv.data(),
            toGraphADC[r][c][iled].xe.data(),toGraphADC[r][c][iled].ye.data());
          TString tmp = TString::Format("[%2d,%2d,%2d]",r+1+kRowOffset,c+1+kColOffset,led);
          gADC[r][c][iled]->SetMarkerStyle(20);
          gADC[r][c][iled]->SetTitle(TString::Format("%02d-%02d LED %02d",r+1+kRowOffset,c+1+kColOffset,led));
          gADC[r][c][iled]->Draw("AP");
*/
          if(led>=0) {
            resADC = ResLED_t(r+1+kRowOffset,c+1+kColOffset,led,ComputeStatsADC(vADC[r][c][iled]));
            //res = ResLED_t(r+1+kRowOffset,c+1+kColOffset,led,ComputeStats(vNPE[r][c][iled]));
            //res.npeE = hNPE[r][c][iled]->GetMeanError();;
            //res.npe  = hNPE[r][c][iled]->GetMean();
            //res.npeE = hNPE[r][c][iled]->GetRMS();
            // New technique, just compute average from graph myself
            //ComputeStats(vNPE[r][c][iled][0],res.npe,res.npeE,res.npeS,res.npeN);
            //std::cerr << "r: " << res.row << ", c: " << res.col << ", npe: " << res.npe << std::endl;
            //if(res.npe.avg != res.npe.avg) {
            //  std::cerr << "Something happened to " << res.row << "-" << res.col <<  std::endl;
            //}
            //res.nentries = hNPE[r][c][iled]->GetEntries();
            //results[r][c].push_back(res);
            resultsADC[r][c].push_back(resADC);
            // Output a summary of results
            if(readLimsFromFile) {
            outsummary << resADC.row << " " << resADC.col << " " << resADC.led
              << " " << resADC.npe.avg << " " << resADC.npe.err
              << " " << resADC.npe.sig << " " << resADC.npe.n
              << std::endl;
            }
          }
        }
      }
    }
    if(saveCanvas&&canvas) {
      canvas->SaveAs(TString::Format("plots/%d_LED%02d_npe.png",run,led));
      canvasADC->SaveAs(TString::Format("plots/%d_LED%02d_adc.png",run,led));
    }
  }
  outsummary.close();

  // print result
  std::cout << "#Results for Run: " << run << std::endl;
  printf("\t");
  for(size_t l = 0; l < resultsADC[0][0].size(); l++) {
    printf("%7d ", resultsADC[0][0][l].led);
  }
  printf("\n");
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) { 
      //if(r==4&&c==1)
      //  continue; // Skip row 5-2 which is bad anyways!
      //std::cout << TString::Format("LED %2d",results[r][c][l].led);
      printf("[%02d,%02d] ",r+1+kRowOffset,c+1+kColOffset);
/*      for(size_t l = 0; l < results[r][c].size(); l++) {
        printf("%7.1f ", results[r][c][l].npe.avg);
      }
      printf("\n\t");
      for(size_t l = 0; l < results[r][c].size(); l++) {
        printf("%7.1f ", results[r][c][l].npe.err);
      }
      printf("\n\t");*/
      for(size_t l = 0; l < resultsADC[r][c].size(); l++) {
        printf("%7.1f ", resultsADC[r][c][l].npe.avg);
      }
      printf("\n\t");
      for(size_t l = 0; l < resultsADC[r][c].size(); l++) {
        printf("%7.1f ", resultsADC[r][c][l].npe.sig);
      }
      printf("\n\n");
    }
  }
  /*
   printf("\n\n=== ADC Values ===\n");
   printf("\t");
  for(size_t l = 0; l < resultsADC[0][0].size(); l++) {
    printf("%7d ", resultsADC[0][0][l].led);
  }
  printf("\n");
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) { 
      //if(r==4&&c==1)
      //  continue; // Skip row 5-2 which is bad anyways!
      //std::cout << TString::Format("LED %2d",results[r][c][l].led);
      printf("[%02d,%02d] ",r+1+kRowOffset,c+1+kColOffset);
      for(size_t l = 0; l < resultsADC[r][c].size(); l++) {
        printf("%7.1f ", resultsADC[r][c][l].npe.avg);
      }
      printf("\n\t");
      for(size_t l = 0; l < resultsADC[r][c].size(); l++) {
        printf("%7.1f ", resultsADC[r][c][l].npe.err);
      }
      printf("\n");
    }
  }
*/

      /* std::cout << TString::Format("[%02d,%02d] = ",r+1,c+1);
      for(size_t l = 0; l < results[r][c].size(); l++) {
        std::cout << TString::Format("  (%2d) %10.3f %10.3f %10d",results[r][c][l].led,
            results[r][c][l].npe.avg,results[r][c][l].npe.err,
            results[r][c][l].npe.n);
      }
      std::cout << std::endl; */
    //std::cout << TString::Format("LED %2d %10.3f %10.3f %10.3f",led,npe,avg,sig) << std::endl;
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      for(iled = 0; iled < kNLED; iled++) {
        if(gSaturated[r][c][iled]) {
           led = iled;
           std::cout << "*** Saturated *** ADC [" << r+1+kRowOffset << ", " << c+1+kColOffset << "] LED " << led << std::endl; 
        }
      }
    }
  }
}
