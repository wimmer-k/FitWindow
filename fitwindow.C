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
#include "TH2F.h"
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
TF1 *fbg;
TF1* fpeak[5];
TH1F* fh;
TH2F* fh2;
bool fcommonwidth = false;
bool fnobg = false;
bool flog = false;
bool frebinpressed = false;
int flinew = 1;
Double_t multgausbg(Double_t *x, Double_t *par);
void ClickFit();
void LiveFit();
void Clear();
void Fit();
void Zoom();
void UnZoom();
void CrossHair();
void LogY();
void Project(int axis);
void Rebin(int reb);
void SetLineWidth(int lw){
  flinew = lw;
}
void SetCommonWidth(bool set =true){
  fcommonwidth = set;
  cout << " common width ";
  if(fcommonwidth)
    cout << "ON " << endl;
  else
    cout << "OFF " << endl;
}
void UnSetCommonWidth(){
  fcommonwidth = false;
  cout << " common width OFF" << endl;
}
void SetNoBG(bool set =true){
  fnobg = set;
  cout << " no background ";
  if(fnobg)
    cout << "ON " << endl;
  else
    cout << "OFF " << endl;
}
void SetLogLikelyhood(bool set =true){
  flog= set;
  cout << " log likelyhood ";
  if(flog)
    cout << "ON " << endl;
  else
    cout << "OFF " << endl;
}
void UnSetNoBG(){
  fnobg = false;
  cout << " no background OFF" << endl;
}
void window(){
  fmarkers.clear();
  fpoints.clear();
  fcanvas = new TCanvas("window","window",600,400);
  fcanvas->AddExec("fit","ClickFit()");
  //fcanvas->AddExec("fit","LiveFit()");
  fh = NULL;
  fh2 = NULL;
  fnpeaks =1;
}
void LiveFit(){
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
  //cout << "clicked " << gPad->GetEvent() << endl;
  //Only act on rightclick-clicks.
  if(gPad->GetEvent() == 1){
    Clear();
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);
    if(yp>fh->GetMinimum() && yp<fh->GetMaximum()){
      cout << "point: " << fpoints.size() << " selected at "<< xp << "\t" << yp<< endl;
      fpoints.push_back(xp);
      gPad->SetEditable(0);
    }
  }
  if(gPad->GetEvent() == 21){
    double xp = gPad->GetEventX();
    double yp = gPad->GetEventY();
    xp = gPad->AbsPixeltoX((int)xp);
    yp = gPad->AbsPixeltoY((int)yp);
    if(yp>fh->GetMinimum() && yp<fh->GetMaximum()){
      cout << "point: " << fpoints.size() << " selected at "<< xp << "\t" << yp<< endl;
      fpoints.push_back(xp);
      Fit();
      fpoints.pop_back();
      gPad->Modified();
      gPad->Update();
      gSystem->ProcessEvents();      
    }
  }
  if(gPad->GetEvent() == 11){
    gPad->SetEditable(1);
    Clear();
  }
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
      fh2 = NULL;
    }
    if(obj->InheritsFrom("TH2")){
      //cout << "histo: " << obj->GetName() << endl;
      fh2 = (TH2F*)obj;
      fh = NULL;
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
    cout << "point: " << fpoints.size() << " selected at "<< xp << "\t" << yp<< endl;
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
  //cout << "event " << gPad->GetEvent() << endl;
  if(gPad->GetEvent() ==kKeyPress){
    //cout << "pressed " << gPad->GetEventX() << endl;
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
    case 'h': //crosshair
      CrossHair();
      break;
    case 'l': //log Y
      LogY();
      break;
    case 'x': //project x
      Project(0);
      break;
    case 'y': //project x
      Project(1);
      break;
    case 'r': //rebin
      cout << "rebinning " << flush;
      frebinpressed = true;    
      break;
    default:
      break;
    };
    if(frebinpressed && gPad->GetEventX()>48 && gPad->GetEventX()<57){
      cout << "factor " << gPad->GetEventX()-48 << endl;
      Rebin(gPad->GetEventX()-48);
      frebinpressed = false;
    }
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
  cout << "Zoom" <<endl;
  if(fpoints.size() >2){
    cout << "too many points for zooming! Clearing!" << endl;
    Clear();
    return;    
  }
  sort(fpoints.begin(), fpoints.end());
  cout << fpoints[0] <<"\t"<<fpoints[1]<<endl;
  fh->GetXaxis()->SetRangeUser(fpoints[0],fpoints[1]);
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
  Clear();
}
void UnZoom(){
  cout << "UnZoom" <<endl;
  if(fh){
    fh->GetXaxis()->UnZoom();
    fh->GetYaxis()->UnZoom();
  }
  if(fh2){
    fh2->GetXaxis()->UnZoom();
    fh2->GetYaxis()->UnZoom();
    fh2->GetZaxis()->UnZoom();
  }
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}
void CrossHair(){
  cout << "CrossHair" <<endl;
  if(gPad->HasCrosshair())
    gPad->SetCrosshair(0);
  else
    gPad->SetCrosshair(1);
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}
void LogY(){
  TIter next (fcanvas->GetListOfPrimitives());
  TObject *obj;
  if(fh){
    cout << "LogY" <<endl;
    if(gPad->GetLogy())
      gPad->SetLogy(0);
    else
      gPad->SetLogy(1);
  }
  else if(fh2){
    cout << "LogZ" <<endl;
    if(gPad->GetLogy())
      gPad->SetLogy(0);
    if(gPad->GetLogz())
      gPad->SetLogz(0);
    else
      gPad->SetLogz(1);
  }
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}
void Project(int axis){
  cout << "Project on ";
  if(axis == 0)
    cout << "X";
  if(axis == 1)
    cout << "Y";
  cout << " axis" <<endl;
  if(!fh2){
    cout << " not a 2D histogram" << endl;
  }
  if(axis == 0)
    fh2->ProjectionX()->Draw();
  if(axis == 1)
    fh2->ProjectionY()->Draw();
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}

void Rebin(int reb){
  //fh = (TH1F*)fh->Rebin(reb,"reb");
  cout << gPad->GetUxmin() <<"\t"<< gPad->GetUxmax() << endl;
  double min =  gPad->GetUxmin();
  double max =  gPad->GetUxmax();
  fh->Rebin(reb);
  fh->GetXaxis()->SetRangeUser(min,max);
  fh->Draw();
  gPad->Modified();
  gPad->Update();
  gSystem->ProcessEvents();
}
void Fit(){
  if(fpoints.size()<2){
    cout << "too few points for fitting! " << endl;
    return;
  }

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
  ffit->SetLineWidth(flinew);


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
  if(flog)
    fh->Fit(ffit,"RL");
  else
    fh->Fit(ffit,"R");
  //ffit->Draw("same");
  cout << "      Chi Square: " << ffit->GetChisquare() << endl;
  cout << "      FWHM:       " << 2*ffit->GetParameter(4)*sqrt(2*log(2)) << "\t" <<2*ffit->GetParameter(4)*sqrt(2*log(2))/ffit->GetParameter(3) << endl;

  fbg = new TF1("bg",multgausbg,fpoints[0],fpoints.back(),npar);
  fbg->SetLineColor(2);
  fbg->SetLineWidth(flinew);
  for(int i=0;i<2;i++){
    fbg->SetParameter(i,ffit->GetParameter(i));
  }
  for(int i=2;i<npar;i++){
    fbg->SetParameter(i,0);
  }
  fbg->Draw("same");
  
  for(int j=0;j<fnpeaks;j++){
    cout << "peak " << j << endl;
    fpeak[j] = new TF1(Form("peak_%d",j),multgausbg,fpoints[0],fpoints.back(),npar);
    fpeak[j]->SetLineColor(4);
    fpeak[j]->SetLineWidth(flinew);
    for(int i=0;i<2;i++)
      fpeak[j]->SetParameter(i,0);
    for(int i=2;i<npar;i++)
      fpeak[j]->SetParameter(i,ffit->GetParameter(i));
    for(int i=0;i<fnpeaks;i++){
      if(i!=j)
	fpeak[j]->SetParameter(2+i*3,0);
      if(fcommonwidth)
	fpeak[j]->SetParameter(2+i*3+2,ffit->GetParameter(4));
    }
    // for(int i=0;i<npar;i++){
    //   cout << "setting par " << i << " to " << fpeak[j]->GetParameter(i) << endl;
    // }    
    fpeak[j]->Draw("same");
  }

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
    if(sigma==0)
      continue;
    arg = (x[0]-mean)/(sqrt2*sigma);
    result += fbinw/(sqrt2pi*sigma) * norm * exp(-arg*arg);
  }
  

  return result;
}
