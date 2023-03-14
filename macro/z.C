#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
void z()
{
  TFile *f[2];
  f[0]=new TFile("totFAZIAPRE_4840.root");
  f[1]=new TFile("totFAZIAPRE_4825.root");
    TH2F *h[2];
 TH1D *pz[2];
 for(int j=0;j<2;j++)
   {
     h[j]=(TH2F*)f[j]->Get("hzv");
     pz[j]=(TH1D*)h[j]->ProjectionY(Form("pz%d",j),-1,-1);
    }
 TCanvas *c=new TCanvas("c","c");
 c->Draw();
 pz[0]->SetLineColor(2);
 pz[0]->Draw();
 pz[0]->GetXaxis()->SetTitle("Z");
 pz[1]->Draw("same");
 TLegend *leg=new TLegend(0.55,0.6,0.9,0.9);
 leg->SetTextSize(0.04);
 leg->AddEntry(pz[0],"48Ca+12C@40AMeV","l");
 leg->AddEntry(pz[1],"48Ca+12C@25AMeV","l");
 leg->Draw();
 gPad->SetGridx(kFALSE);
 gPad->SetGridy(kFALSE);

}
