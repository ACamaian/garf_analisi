#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void cfr()
{
  gStyle->SetGridWidth(0);


  TFile *f[3];
  char tit[3][100]={"4025","4825","4840"};
  TH1F *hmolt[3];
  TH2F*hthetarelvrel[3];
  TH2F*hthetarelvrelge5[3];
  TH2F *hthetarelvrelsumztge5[3];
  TGraphErrors *gnzbig[3];
  TGraphErrors *gsigmabig[3];

  TGraphErrors *gzbig[3][27];
  TGraphErrors *gnzbiggest[3];
  TGraphErrors *gnzbiggestm2[3];
  TGraphErrors *gnzsecond[3];

  int col[3]={1,2,3};
  int mark[3]={20,20,24};
  int marke[2]={20,24};
  for(int j=0;j<3;j++)
    {
      f[j]=new TFile(Form("feb4/totFAZIAPRE_%s.root",tit[j]));
      hmolt[j]=(TH1F*)f[j]->Get("hm");

hthetarelvrel[j]=(TH2F*)f[j]->Get("hthetarelvrel");
hthetarelvrelge5[j]=(TH2F*)f[j]->Get("hthetarelvrelge5");
 hthetarelvrelsumztge5[j] =(TH2F*)f[j]->Get("hthetarelvrelsumztge5");
      gnzbig[j]=(TGraphErrors*)f[j]->Get("gnzbig");
      gnzbiggest[j]=(TGraphErrors*)f[j]->Get("gnzbiggesttutti");
      gnzbiggestm2[j]=(TGraphErrors*)f[j]->Get("gnzbiggest");
      gsigmabig[j]=(TGraphErrors*)f[j]->Get("gsigmanzbig");
      gnzsecond[j]=(TGraphErrors*)f[j]->Get("gnzsecond");

 gnzbig[j]->SetMarkerColor(col[j]);
 gnzbig[j]->SetLineColor(col[j]);
 gnzsecond[j]->SetMarkerColor(col[j]);
 gnzsecond[j]->SetLineColor(col[j]);

 gnzbig[j]->SetMarkerStyle(mark[j]);
 gnzbiggest[j]->SetMarkerColor(col[j]);
 gnzbiggest[j]->SetLineColor(col[j]);

 gnzbiggest[j]->SetMarkerStyle(mark[j]);
 gnzbiggestm2[j]->SetMarkerColor(col[j]);
 gnzbiggestm2[j]->SetLineColor(col[j]);

 gnzbiggestm2[j]->SetMarkerStyle(20);
 gnzsecond[j]->SetMarkerStyle(24);

 gsigmabig[j]->SetMarkerColor(col[j]);
 gsigmabig[j]->SetLineColor(col[j]);
 gsigmabig[j]->SetMarkerStyle(mark[j]);
      for(int iz=10;iz<27;iz++)
	{
	  gzbig[j][iz]=0;
	  gzbig[j][iz]=(TGraphErrors*)f[j]->Get(Form("gz_%d",iz));
	  if(gzbig[j][iz]!=0)
	    {
	      if(gzbig[j][iz]->GetN()>0)
		{
		  gzbig[j][iz]->SetMarkerColor(col[j]);
		  gzbig[j][iz]->SetLineColor(col[j]);
		  gzbig[j][iz]->SetMarkerStyle(mark[j]);

		}
	    }
	}
    }

  TCanvas *cnz=new TCanvas("cnz","cnz");
  cnz->Draw();
  TH2F *hnz=new TH2F("hnz","hnz",20,9,27,100,1.,1.4);
  hnz->Draw();
  hnz->GetXaxis()->SetTitle("Z");
  hnz->GetYaxis()->SetTitle("<N>/Z");
  hnz->SetStats(kFALSE);
  TLegend *legnz=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      gnzbig[j]->Draw("pl");
      legnz->AddEntry(gnzbig[j],tit[j],"p");
    }
  legnz->Draw();
  cnz->Print("nzbig.pdf");

  TCanvas *cnzbiggest=new TCanvas("cnzbiggest","cnzbiggest");
  cnz->Draw();

  TH2F *hnzbiggest=new TH2F("hnzbiggest","hnzbiggest",20,9,27,100,1.,1.4);
  hnzbiggest->Draw();
  hnzbiggest->GetXaxis()->SetTitle("Z");
  hnzbiggest->GetYaxis()->SetTitle("<N>/Z");
  hnzbiggest->SetStats(kFALSE);
  TLegend *legnzbiggest=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      gnzbiggest[j]->Draw("pl");
      legnzbiggest->AddEntry(gnzbiggest[j],tit[j],"p");
    }
  legnzbiggest->Draw();
  cnzbiggest->Print("nzbiggest.pdf");


  TCanvas *cnzbiggestm2=new TCanvas("cnzbiggestm2","cnzbiggestm2");
  cnzbiggestm2->Draw();

  TH2F *hnzbiggestm2=new TH2F("hnzbiggestm2","hnzbiggestm2",20,9,27,100,1.,1.4);
  hnzbiggestm2->Draw();
  hnzbiggestm2->GetXaxis()->SetTitle("Z");
  hnzbiggestm2->GetYaxis()->SetTitle("<N>/Z");
  hnzbiggestm2->SetStats(kFALSE);
  TLegend *legnzbiggestm2=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      gnzbiggestm2[j]->Draw("pl");
      legnzbiggestm2->AddEntry(gnzbiggestm2[j],tit[j],"p");
    }
  legnzbiggestm2->Draw();
  cnzbiggestm2->Print("nzbiggestm2.pdf");

 TCanvas *cnz12=new TCanvas("cnz12","cnz12");
  cnz12->Draw();

  TH2F *hnz12=new TH2F("hnz12","hnz12",20,0,27,100,0.1,1.4);
  hnz12->Draw();
  hnz12->GetXaxis()->SetTitle("Z");
  hnz12->GetYaxis()->SetTitle("<N>/Z");
  hnz12->SetStats(kFALSE);
  TLegend *legnz12=new TLegend(0.6,0.1,0.9,0.3);
  for(int j=0;j<3;j++)
    {
      gnzbiggestm2[j]->Draw("pl");
      gnzsecond[j]->Draw("pl");
      legnz12->AddEntry(gnzbiggestm2[j],Form("%s_biggest",tit[j]),"p");
      legnz12->AddEntry(gnzsecond[j],Form("%s_second",tit[j]),"p");
    }
  legnz12->Draw();
  cnz12->Print("nz12.pdf");



  TCanvas *csigma=new TCanvas("csigma","csigma");
  csigma->Draw();
  TH2F *hsigma=new TH2F("hsigma","hsigma",20,9,27,100,0.,3.2);
  hsigma->Draw();
  hsigma->GetXaxis()->SetTitle("Z");
  hsigma->GetYaxis()->SetTitle("sigma");
  hsigma->SetStats(kFALSE);
  TLegend *legsigma=new TLegend(0.2,0.7,0.5,0.9);
  for(int j=0;j<3;j++)
    {
      gsigmabig[j]->Draw("pl");
      legsigma->AddEntry(gsigmabig[j],tit[j],"p");
    }
  legsigma->Draw();
  csigma->Print("sigmabig.pdf");

  TCanvas *cn1 =new TCanvas("cn1","cn1",0,0,2000,1000);
  cn1->Divide(4,2);
  cn1->Draw();

  TH2F *hn1[8];
  for(int j=0;j<8;j++)
    {cn1->cd(j+1);

      if(j+10<14)hn1[j]=new TH2F(Form("hn1_%d",j),Form("Z=%d",10+j),20,8,20,100,0,1);
if(j+10>=14)hn1[j]=new TH2F(Form("hn1_%d",j),Form("Z=%d",10+j),20,12,25,100,0,1);
      hn1[j]->GetXaxis()->SetTitle("N");
      hn1[j]->SetStats(kFALSE);
      hn1[j]->Draw();
      
      for(int k=0;k<3;k++)
	{
	  gzbig[k][j+10]->Draw("pl");
	}
gPad->SetLogy(kTRUE);
    }
  cn1->cd(1);
  TLegend *legn1=new TLegend(0.5,0.5,0.9,0.9);
  for(int j=0;j<3;j++)
    {
      legn1->AddEntry(gzbig[j][10],tit[j],"p");
    }
  legn1->Draw();
  cn1->Print("cn1.pdf");
  TCanvas *cn2 =new TCanvas("cn2","cn2",0,0,2000,1000);
  cn2->Divide(4,2);
  cn2->Draw();

  TH2F *hn2[8];
  for(int j=0;j<5;j++)
    {cn2->cd(j+1);

      if(j+18<20)hn2[j]=new TH2F(Form("hn2_%d",j),Form("Z=%d",18+j),20,16,30,100,0,1);
if(j+18>=20)hn2[j]=new TH2F(Form("hn2_%d",j),Form("Z=%d",18+j),20,18,40,100,0,1);
      hn2[j]->GetXaxis()->SetTitle("N");
      hn2[j]->SetStats(kFALSE);
      hn2[j]->Draw();
      
      for(int k=0;k<3;k++)
	{
	  gzbig[k][j+18]->Draw("pl");
	}
gPad->SetLogy(kTRUE);
    }
  cn2->cd(1);
  TLegend *legn2=new TLegend(0.5,0.5,0.9,0.9);
  for(int j=0;j<3;j++)
    {
      legn2->AddEntry(gzbig[j][18],tit[j],"p");
    }
  legn2->Draw();
  cn2->Print("cn2.pdf");



  TCanvas *cmolt=new TCanvas("cmolt","cmolt");
  cmolt->Draw();
   
  TLegend *legmolt=new TLegend(0.7,0.7,0.9,0.9);
  for(int j=2;j>=0;j--)
    {
      hmolt[j]->SetLineColor(col[j]);
hmolt[j]->SetStats(kFALSE);
 hmolt[j]->GetXaxis()->SetRangeUser(0,15);
 if(j==2)hmolt[j]->Draw();
 else hmolt[j]->Draw("same");
 gPad->SetLogy(kTRUE);
      legmolt->AddEntry(hmolt[j],tit[j],"l");
    }
  legmolt->Draw();
  cmolt->Print("molt.pdf");

//   TFile *fout[3];
//   for(int j=0;j<3;j++)
//     {
//       fout[j]=new TFile(Form("out%s.root",tit[j]),"RECREATE");
//       fout[j]->cd();
//       fout[j]->Add(gnzbig[j]);
//       fout[j]->Add(gnzbiggest[j]);
//       fout[j]->Add(gnzbiggestm2[j]);
//       fout[j]->Add(gnzsecond[j]);
//       fout[j]->Write();
      
//       fout[j]->Close();
//     }

  TCanvas *c4 =new TCanvas("c4","c4",0,0,500,1000);
  c4->Divide(1,3);
  c4->Draw();
  for(int j=0;j<3;j++)
    {
      c4->cd(j+1);
      hthetarelvrel[j]->SetTitle(Form("%s",tit[j]));
      hthetarelvrel[j]->Draw("zcol");
      gPad->SetLogz(kTRUE);
    }
  c4->Print("thetarelvrel.pdf");
  TCanvas *c5 =new TCanvas("c5","c5",0,0,500,1000);
  c5->Divide(1,3);
  c5->Draw();
  for(int j=0;j<3;j++)
    {
      c5->cd(j+1);
      hthetarelvrelge5[j]->SetTitle(Form("%s",tit[j]));
      hthetarelvrelge5[j]->Draw("zcol");
      gPad->SetLogz(kTRUE);
    }
  c5->Print("thetarelvrelge5.pdf");
  TCanvas *c6 =new TCanvas("c6","c6",0,0,500,1000);
  c6->Divide(1,3);
  c6->Draw();
  for(int j=0;j<3;j++)
    {
      c6->cd(j+1);
      hthetarelvrelsumztge5[j]->SetTitle(Form("%s",tit[j]));
      hthetarelvrelsumztge5[j]->Draw("zcol");
      gPad->SetLogz(kTRUE);
    }
  c6->Print("thetarelvrelsumztge5.pdf");
}
