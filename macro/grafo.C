#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void grafo()
{
  gStyle->SetGridWidth(0);


  TFile *f[3];
  char tit[3][100]={"4025","4825","4840"};
  TH1F *haa_erel[3];
  TH1F *haa_emix[3];
  TH1F *hda_erel[3];
  TH1F *hda_emix[3];
  TH1F *haa[3];
  TH1F *hda[3];


  int col[3]={1,2,3};
  int mark[3]={20,20,24};
  int marke[2]={20,24};
  for(int j=0;j<3;j++)
    {
      f[j]=new TFile(Form("feb4/totFAZIAPRE_%s.root",tit[j]));

      haa_erel[j]=(TH1F*)f[j]->Get("haa_erel");
      haa_emix[j]=(TH1F*)f[j]->Get("haa_emix");
      hda_erel[j]=(TH1F*)f[j]->Get("hda_erel");
      hda_emix[j]=(TH1F*)f[j]->Get("hda_emix");
      haa[j]=(TH1F*)haa_erel[j]->Clone(Form("haa_%d",j));
      hda[j]=(TH1F*)hda_erel[j]->Clone(Form("hda_%d",j));
      haa[j]->Divide(haa_emix[j]);
      hda[j]->Divide(hda_emix[j]);


    }

  TCanvas *caa=new TCanvas("caa","caa");
  caa->Draw();

  TLegend *legaa=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      haa[j]->SetLineColor(col[j]);
      if(j==0)haa[j]->Draw();
      if(j>0) haa[j]->Draw("same");
      legaa->AddEntry(haa[j],tit[j],"l");
    }
  legaa->Draw();
  caa->Print("aa.pdf");
  TCanvas *cda=new TCanvas("cda","cda");
  cda->Draw();

  TLegend *legda=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      hda[j]->SetLineColor(col[j]);
      if(j==0)hda[j]->Draw();
      if(j>0) hda[j]->Draw("same");
      legda->AddEntry(hda[j],tit[j],"l");
    }
  legda->Draw();
  cda->Print("da.pdf");


  TCanvas *cda4=new TCanvas("cda4","cda4");
  cda4->Divide(2,2);
  cda4->Draw();

  TLegend *legda4=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      cda4->cd(j+1);

      hda[j]->Draw();
      
      legda4->AddEntry(hda[j],tit[j],"l");
    }
  cda4->cd(4);
  legda4->Draw();
  cda4->Print("da4.pdf");

  TCanvas *caa4=new TCanvas("caa4","caa4");
  caa4->Divide(2,2);
  caa4->Draw();

  TLegend *legaa4=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      caa4->cd(j+1);

      haa[j]->Draw();
      
      legaa4->AddEntry(haa[j],tit[j],"l");
    }
  caa4->cd(4);
  legaa4->Draw();
  caa4->Print("aa4.pdf");


}
