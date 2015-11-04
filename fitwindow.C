#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>  

#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TImage.h"
#include <TSystem.h>
#include <TMarker.h>
#include "TMath.h"
TCanvas *fcanvas;
vector<TMarker*> fmarkers;
vector<double> fpoints;
Double_t fbinw;
Int_t fnpeaks;
TF1 *ffit;
TH1F* fh;
bool fcommonwidth = false;
bool fnobg = false;
Double_t multgausbg(Double_t *x, Double_t *par);
void ClickFit();
void Clear();
void Fit();
void Zoom();
void UnZoom();
void SetCommonWidth(bool set =true){
  fcommonwidth = set;
}
void UnSetCommonWidth(){
  fcommonwidth = false;
}
void SetNoBG(bool set =true){
  fnobg = set;
}
void UnSetNoBG(){
  fnobg = false;
}
void window(){
  fmarkers.clear();
  fpoints.clear();
  fcanvas = new TCanvas("window","window",600,400);
  fcanvas->AddExec("fit","ClickFit()");
  //TFile *f = new TFile("gamma_new.root");
  // fh = (TH1F*)f->Get("hgamma");
  // fh->Draw();
  fh = NULL;
  fnpeaks =1;
}
void ClickFit(){
  TObject *select = gPad->GetSelected();
  if(!select)
    return;

  TIter next (fcanvas->GetListOfPrimitives());
  TObject *obj;
  while((obj=next())){
    //cout << "reading: " << obj->GetName() << endl;
    if(obj->InheritsFrom("TH1")){
      //cout << "histo: " << obj->GetName() << endl;
      fh = (TH1F*)obj;
    }
  }

  gPad->GetCanvas()->FeedbackMode(kTRUE);
  //Only act on middle-clicks.
  if(gPad->GetEvent() == 12){
    //cout << "middle button " << endl;
    //Find the location of the click.
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();
    //cout << xp <<"\t" << yp << endl;
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);
    cout << "point: " << fpoints.size() << " selected at "<< xp << endl;
    fpoints.push_back(xp);
    fmarkers.push_back(new TMarker(xp,yp,2));
    fmarkers.back()->SetMarkerSize(2);
    fmarkers.back()->SetMarkerColor(2);
    for(UShort_t i=0;i<fmarkers.size()-1;i++){
      fmarkers[i]->SetMarkerColor(1);
      fmarkers[i]->Draw("same");
    }
    fmarkers.back()->Draw("same");
    gPad->Modified();
    gPad->Update();
    gSystem->ProcessEvents();

  }
  if(gPad->GetEvent() ==kKeyPress){
    switch(gPad->GetEventX()){
    case 'f': //fit
      Fit();
      break;
    case 'c': //clear
      Clear();
      break;
    case 'z': //zoom
      Zoom();
      break;
    case 'u': //unzoom
      UnZoom();
      break;
    default:
      break;
    };
  }
}
void Addpoint(double xp){
  double yp = 0;
  yp = gPad->AbsPixeltoY((int)yp);
  cout << "point: " << fpoints.size() << " selected at "<< xp << endl;
  fpoints.push_back(xp);
  fmarkers.push_back(new TMarker(xp,yp,2));
  fmarkers.back()->SetMarkerSize(2);
  fmarkers.back()->SetMarkerColor(2);
  for(UShort_t i=0;i<fmarkers.size()-1;i++){
    fmarkers[i]->SetMarkerColor(1);
    fmarkers[i]->Draw("same");
    }
  fmarkers.back()->Draw("same");
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();

}
void Zoom(){
  if(fpoints.size() >2){
    cout << "too many points for zooming! Clearing!" << endl;
    Clear();
    return;    
  }
  sort(fpoints.begin(), fpoints.end());
  fh->GetXaxis()->SetRangeUser(fpoints[0],fpoints[1]);
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
  Clear();
}
void UnZoom(){
  fh->GetXaxis()->UnZoom();
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}
void Fit(){
  sort(fpoints.begin(), fpoints.end());
  fbinw = fh->GetBinWidth(1);
  //cout << fbinw << endl;
  double slope;
  double bg;
  double sigma;
  double content;
  double mean;
  int npar =5;
  if(fpoints.size()==4)
    npar = 8;
  if(fpoints.size()==5)
    npar = 11;
  ffit = new TF1("fit",multgausbg,fpoints[0],fpoints.back(),npar);
  ffit->SetLineColor(3);
  ffit->SetLineWidth(1);


  if(fpoints.size()==3){
    //cout << fpoints[0] << "\t"<<fpoints[1] << "\t"<<fpoints[2] <<endl;

    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[2]));
    slope/=(fpoints[0]-fpoints[2]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[2]));
    bg-=slope*(fpoints[0]+fpoints[2]);
    bg/=2;

    sigma = fbinw*5;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[2]));
    mean = fpoints[1];
    fnpeaks = 1;
    ffit->SetParameters(bg,slope,content,mean,sigma);
    ffit->SetParName(0,"BgConstant");
    ffit->SetParName(1,"BgSlope   ");
    ffit->SetParName(2,"Content   ");
    ffit->SetParName(3,"Mean      ");
    ffit->SetParName(4,"Sigma     ");
  }
  else if(fpoints.size()==2){
    //cout << fpoints[0] << "\t"<<fpoints[1] <<endl;
    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[1]));
    slope/=(fpoints[0]-fpoints[1]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[1]));
    bg-=slope*(fpoints[0]+fpoints[1]);
    bg/=2;

    sigma = fbinw*5;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[1]));
    mean = fpoints[0]+fpoints[1];
    mean/=2;
    fnpeaks = 1;
    ffit->SetParameters(bg,slope,content,mean,sigma);
    ffit->SetParName(0,"BgConstant");
    ffit->SetParName(1,"BgSlope   ");
    ffit->SetParName(2,"Content   ");
    ffit->SetParName(3,"Mean      ");
    ffit->SetParName(4,"Sigma     ");

  }
  else if(fpoints.size()==4){
    //cout << fpoints[0] << "\t"<<fpoints[1] <<endl;
    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[3]));
    slope/=(fpoints[0]-fpoints[3]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[3]));
    bg-=slope*(fpoints[0]+fpoints[3]);
    bg/=2;

    sigma = fbinw;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[3]));
    content/= 2;
    fnpeaks = 2;
    ffit->SetParameters(bg,slope,content,fpoints[1],sigma,content,fpoints[2],sigma);
    ffit->SetParName(0,"BgConstant");
    ffit->SetParName(1,"BgSlope   ");
    ffit->SetParName(2,"Content0  ");
    ffit->SetParName(3,"Mean0     ");
    ffit->SetParName(4,"Sigma0    ");
    ffit->SetParName(5,"Content1  ");
    ffit->SetParName(6,"Mean1     ");
    ffit->SetParName(7,"Sigma1    ");

  }
  else if(fpoints.size()==5){
    //cout << fpoints[0] << "\t"<<fpoints[1] <<endl;
    slope = fh->GetBinContent(fh->FindBin(fpoints[0])) - fh->GetBinContent(fh->FindBin(fpoints[4]));
    slope/=(fpoints[0]-fpoints[4]);
    bg = fh->GetBinContent(fh->FindBin(fpoints[0])) + fh->GetBinContent(fh->FindBin(fpoints[4]));
    bg-=slope*(fpoints[0]+fpoints[4]);
    bg/=2;

    sigma = fbinw;
    content = fh->Integral(fh->FindBin(fpoints[0]),fh->FindBin(fpoints[4]));
    content/= 2;
    fnpeaks = 3;
    ffit->SetParameters(bg,slope,content,fpoints[1],sigma,content,fpoints[2],sigma,content,fpoints[3],sigma);
    ffit->SetParName(0,"BgConstant");
    ffit->SetParName(1,"BgSlope   ");
    ffit->SetParName(2,"Content0  ");
    ffit->SetParName(3,"Mean0     ");
    ffit->SetParName(4,"Sigma0    ");
    ffit->SetParName(5,"Content1  ");
    ffit->SetParName(6,"Mean1     ");
    ffit->SetParName(7,"Sigma1    ");
    ffit->SetParName(8,"Content2  ");
    ffit->SetParName(9,"Mean2     ");
    ffit->SetParName(10,"Sigma2    ");

  }
  else if(fpoints.size()==1){
    cout << "too few points for fitting! " << endl;
    return;
  }
  else{
    cout << "too many points for fitting! Clearing!" << endl;
    Clear();
    return;
  }
  if(fnobg){
    ffit->SetParameter(0,0);
    ffit->SetParameter(1,0);
    ffit->FixParameter(0,0);
    ffit->FixParameter(1,0);
  }
  fh->Fit(ffit,"R");
  //ffit->Draw("same");
  cout << "      Chi Square: " << ffit->GetChisquare() << endl;
  cout << "      FWHM:       " << 2*ffit->GetParameter(4)*sqrt(2*log(2)) << "\t" <<2*ffit->GetParameter(4)*sqrt(2*log(2))/ffit->GetParameter(3) << endl;

  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
  Clear();

}
void Clear(){
  for(UShort_t i=0;i<fmarkers.size();i++){
    fmarkers[i]->Delete();
  }
  fmarkers.clear();
  fpoints.clear();
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();      

}
Double_t multgausbg(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss0 content
  par[3]   gauss0 mean
  par[4]   gauss0 width
  par[5]   gauss1 content
  par[6]   gauss1 mean
  par[7]   gauss1 width
  ......
  */
  Double_t result = par[0] + par[1]*x[0];

  for (Int_t p=0;p<fnpeaks;p++) {
    Double_t norm  = par[3*p+2];
    Double_t mean  = par[3*p+3];
    Double_t sigma = par[3*p+4];
    if(fcommonwidth==true)
      sigma  = par[4];
    arg = (x[0]-mean)/(sqrt2*sigma);
    result += fbinw/(sqrt2pi*sigma) * norm * exp(-arg*arg);
  }

  return result;
}
