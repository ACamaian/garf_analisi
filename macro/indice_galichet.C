#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TPaletteAxis.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPolyLine.h>
#include <TMath.h>
void indice_galichet()
{
  TFile *f[3];

  char tit[3][100]={"4025","4825","4840"}; 
  for(int j=0;j<3;j++)
    {
            f[j]=new TFile(Form("feb9/totFAZIAPRE_%s.root",tit[j]));
    }

  TH2F *hsololcp[3];
  TH2F *h[3];
  TH1D *hsololcp1d[3][2];
  TH1D *h1d[3][2];

  TH1F *hneventi[3];
  float norm[3];
  TGraphErrors *g[3];
  for(int j=0;j<3;j++)
    {
      hsololcp[j]=(TH2F*)f[j]->Get("hnzaltrivparmenovbigsololcp");
      h[j]=(TH2F*)f[j]->Get("hnzaltrivparmenovbig");
      hsololcp1d[j][0]=(TH1D*) hsololcp[j]->ProjectionY(Form("hsololcp1d0_%d",j),-1,-1);
      hsololcp1d[j][1]=(TH1D*) hsololcp[j]->ProjectionY(Form("hsololcp1d1_%d",j),hsololcp[j]->GetXaxis()->FindBin(0.),hsololcp[j]->GetXaxis()->GetNbins());
      h1d[j][0]=(TH1D*) h[j]->ProjectionY(Form("h1d0_%d",j),-1,-1);
      h1d[j][1]=(TH1D*) h[j]->ProjectionY(Form("h1d1_%d",j),h[j]->GetXaxis()->FindBin(0.),h[j]->GetXaxis()->GetNbins());
      hneventi[j]=(TH1F*)f[j]->Get("hneventi");
      norm[j]=hneventi[j]->GetBinContent(hneventi[j]->GetXaxis()->FindBin(45.));
  cout<<"norm="<<norm[j]<<endl;
  g[j]=(TGraphErrors*)f[j]->Get("gnzbiggestconlcp");
    }



  float np=0;
  float nn=0;
  float itot=0;
  float dnn=0;
  float dnp=0;
  for(int j=0;j<3;j++)
    {
      for(int k=0;k<2;k++)
	{
      itot=0;
      dnn=0;
      dnp=0;
      np=0;
      nn=0;
      for(int i=0;i<h1d[j][k]->GetXaxis()->GetNbins();i++)
	{
	  if(h1d[j][k]->GetXaxis()->GetBinCenter(i+1)<500 &&h1d[j][k]->GetBinContent(i+1)>0)
	    {
	      float iz=floor(h1d[j][k]->GetXaxis()->GetBinCenter(i+1)/100);
	      float in=floor(h1d[j][k]->GetXaxis()->GetBinCenter(i+1)-iz*100);
	      if(in!=0)
		{
	      np=np+h1d[j][k]->GetBinContent(i+1)*iz;
	      nn=nn+h1d[j][k]->GetBinContent(i+1)*in;
	      //printf("%d %f %f %f\n",i+1,h1[j]->GetXaxis()->GetBinCenter(i+1),iz,in);
	      //cout<<i+1<<" "<<h1[j]->GetXaxis()->GetBinCenter(i+1)<<" "<<iz<<" "<<in<<" "<<h1[j]->GetBinContent(i+1)<<endl;
	      itot=itot+h1d[j][k]->GetBinContent(i+1);
	      dnn=dnn+pow(in*sqrt(h1d[j][k]->GetBinContent(i+1)),2);
	      dnp=dnp+pow(iz*sqrt(h1d[j][k]->GetBinContent(i+1)),2);

		}
	    }
	}
      float  Er=sqrt(pow(sqrt(dnn)/nn,2)+pow(sqrt(dnp)/np,2));
      Er=Er*nn/np;
      cout<<itot<<endl;
      cout<<"tutti="<<j<<" "<<k<<" "<<nn/np<<" "<<Er<<endl;
	}
      for(int k=0;k<2;k++)
	{
      itot=0;
      dnn=0;
      dnp=0;
      np=0;
      nn=0;
      for(int i=0;i<hsololcp1d[j][k]->GetXaxis()->GetNbins();i++)
	{
	  if(hsololcp1d[j][k]->GetXaxis()->GetBinCenter(i+1)<500 &&hsololcp1d[j][k]->GetBinContent(i+1)>0)
	    {
	      float iz=floor(hsololcp1d[j][k]->GetXaxis()->GetBinCenter(i+1)/100);
	      float in=floor(hsololcp1d[j][k]->GetXaxis()->GetBinCenter(i+1)-iz*100);
	      

	      if(in!=0)
		{
	      np=np+hsololcp1d[j][k]->GetBinContent(i+1)*iz;
	      nn=nn+hsololcp1d[j][k]->GetBinContent(i+1)*in;
	      //printf("%d %f %f %f\n",i+1,hsololcp1[j]->GetXaxis()->GetBinCenter(i+1),iz,in);
	      //cout<<i+1<<" "<<hsololcp1[j]->GetXaxis()->GetBinCenter(i+1)<<" "<<iz<<" "<<in<<" "<<hsololcp1[j]->GetBinContent(i+1)<<endl;
	      itot=itot+hsololcp1d[j][k]->GetBinContent(i+1);
	      dnn=dnn+pow(in*sqrt(hsololcp1d[j][k]->GetBinContent(i+1)),2);
	      dnp=dnp+pow(iz*sqrt(hsololcp1d[j][k]->GetBinContent(i+1)),2);

		}
	    }
	}
      float Er=sqrt(pow(sqrt(dnn)/nn,2)+pow(sqrt(dnp)/np,2));
      Er=Er*nn/np;
      cout<<itot<<endl;
      cout<<"solo LCP "<<j<<" "<<k<<" "<<nn/np<<" "<<Er<<endl;


	}


	}
  float mp[3][2];
  float mpsololcp[3][2];

  for(int j=0;j<3;j++)
    {
      for(int k=0;k<2;k++)
	{
	  mp[j][k]=h1d[j][k]->GetBinContent(101)/norm[j];
	
	  mpsololcp[j][k]=hsololcp1d[j][k]->GetBinContent(101)/norm[j];
	  cout<<"mp="<<j<<" "<<k<<" "<<mp[j][k]<<endl;
	  cout<<"mpsololcp="<<j<<" "<<k<<" "<<mpsololcp[j][k]<<endl;
	}
    }
  int col[3]={1,2,2};
  int mark[3]={20,20,24};
  for(int j=0;j<3;j++)
    {
	  g[j]->SetMarkerStyle(mark[j]);
	  g[j]->SetMarkerColor(col[j]);
	  g[j]->SetLineColor(col[j]);

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
      g[j]->Draw("pl");
      legnz->AddEntry(g[j],tit[j],"p");
    }
  legnz->Draw();
  cnz->Print("nzbigconlcp.pdf");



}
