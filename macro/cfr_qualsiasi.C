#include <stdio.h>
#include <iostream>
#include <string.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TStyle.h>
int nmax;
TFile *f[10];
void apri(int nm=3,char *f1="jan27_23/totINDRAFAZIA_64Ni64Ni_32.root",char *f2="jan27_23/totFAZIASYM_484835.root",char *f3="jan27_23/totFAZIAPRE_4840.root",char *f4="a4.root",char *f5="a5.root",char *f6="a6.root", char *f7="a7.root", char *f8="a8.root", char *f9="a9.root", char*f10="a10.root")
{
  nmax=nm;
 char title[10][1000];
  sprintf(title[0],"%s",f1);
  sprintf(title[1],"%s",f2);
  sprintf(title[2],"%s",f3);
  sprintf(title[3],"%s",f4);
  sprintf(title[4],"%s",f5);
  sprintf(title[5],"%s",f6);
  sprintf(title[6],"%s",f7);
  sprintf(title[7],"%s",f8);
  sprintf(title[8],"%s",f9);
  sprintf(title[9],"%s",f10);

 

  for(int j=0;j<nmax;j++)
    {
      f[j]=new TFile(Form("%s",title[j]));
      if(f[j]->IsZombie())
	{
	  cout<<"Il file "<<f[j]->GetName()<<" non esiste"<<endl;
	  return;
	}
    }

}
void spettro(char *nome)
{


  TH1F *h1d[nmax];
  TH2F *h2d[nmax];
  TGraphErrors *g[nmax];
  for(int j=0;j<nmax;j++)
    {
      h1d[j]=0;
      h2d[j]=0;
      g[j]=0;

    }

  TObject *obj[nmax];
  for(int j=0;j<nmax;j++)
    {
      obj[j]=0;
      obj[j]=(TObject*)f[j]->Get(Form("%s",nome));
      if(obj[j]!=0)
	{
	  if(obj[j]->InheritsFrom("TH2"))
	    {
	      cout<<"2d"<<endl;
	      h2d[j]=(TH2F*)obj[j];
	      h2d[j]->SetTitle(Form("%s_%d",nome,j));
	      h2d[j]->SetName(Form("%s_%d",nome,j));


	    }
	  if(obj[j]->InheritsFrom("TH1")&&!obj[j]->InheritsFrom("TH2"))
	    {
	      cout<<"1d"<<endl;
	      h1d[j]=(TH1F*)obj[j];
	      h1d[j]->SetTitle(Form("%s_%d",nome,j));
	      h1d[j]->SetName(Form("%s_%d",nome,j));
	      h1d[j]->SetLineColor(j+1);
	    }
	  if(obj[j]->InheritsFrom("TGraph"))
	    {
	      cout<<"graph"<<endl;
	      g[j]=(TGraphErrors*)obj[j];
	      g[j]->SetTitle(Form("%s_%d",nome,j));
	      g[j]->SetName(Form("%s_%d",nome,j));
	      g[j]->SetLineColor(j+1);
	      g[j]->SetMarkerColor(j+1);

	    }

	}
      else
	{
	  cout<<nome<< " in "<<f[j]->GetName()<<" non esiste"<<endl;
	}
    }
  TCanvas *csomma;
  TCanvas *csommanorm;
  TCanvas *c[nmax];
  for(int j=0;j<nmax;j++)
    {
      if((h1d[j]!=0&&h1d[j]->Integral()>0) || (h2d[j]!=0&&h2d[j]->Integral()>0)||(g[j]!=0&&g[j]->GetN()>0))
	{
      c[j]=new TCanvas(Form("c%d",j),Form("%s-%s",nome,f[j]->GetName()));
      c[j]->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      if(h1d[j]!=0)
	{
	  h1d[j]->Draw();
	}
      if(h2d[j]!=0)
	{
	  gPad->SetLogz(kTRUE);
	  h2d[j]->Draw("zcol");
	}
      if(g[j]!=0)
	{
	  g[j]->Draw("pal");
	}

	}
      if((h1d[j]!=0&&h1d[j]->Integral()==0) || (h2d[j]!=0&&h2d[j]->Integral()==0))
	{
cout<<nome<< " in "<<f[j]->GetName()<<" e' vuoto"<<endl;
	}
    }
  if(nmax>1)
    {
      if(h1d[0]!=0)
	{
	  csomma=new TCanvas("csomma","tutti");
	       gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
	  float inte=h1d[0]->GetMaximum();
	  int jinte=0;
   for(int j=1;j<nmax;j++)
    {
      if(h1d[j]!=0)
	{
	  if(inte<h1d[j]->GetMaximum())
	    {
	      inte=h1d[j]->GetMaximum();
	      jinte=j;
	    }
	}
    }
   h1d[jinte]->Draw();
   for(int j=0;j<nmax;j++)
    {
      if(j!=inte)
	{
	  if(h1d[j]!=0)
	    {
	      h1d[j]->Draw("same");
	    }
	}
    }

 csommanorm=new TCanvas("csommanorm","tuttinorm");
 csommanorm->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
 TH1F *h1dnorm[nmax];
 for(int j=0;j<nmax;j++)
   {
     if(h1d[j]!=0)
       {
	 h1dnorm[j]=(TH1F*)h1d[j]->Clone();
	 h1dnorm[j]->SetName(Form("%s_norm",h1d[j]->GetName()));
	 h1dnorm[j]->Scale(1./h1d[j]->Integral());
       }
   }
 int nn=0;
   for(int j=0;j<nmax;j++)
    {

	  if(h1dnorm[j]!=0)
	    {
	      nn++;
	      if(nn==1)
		{
	      h1dnorm[j]->Draw();
		}
	      else
		{
		  h1dnorm[j]->Draw("same");
		}
	    }

    }

	}

      if(g[0]!=0)
	{
	  csomma=new TCanvas("csomma","tutti");
	  float min=100000;
	  float max=-100000;
	  float xmin=100000;
	  float xmax=-100000;
	  double *x;
	  double *y;
	  for(int k=0;k<nmax;k++)
	    {
	      x=g[k]->GetX();
	      y=g[k]->GetY();
	      for(int j=0;j<g[k]->GetN();j++)
		{
		  if(x[j]<xmin)
		    {
		      xmin=x[j];
		    }
		  if(x[j]>xmax)
		    {
		      xmax=x[j];
		    }
		  if(y[j]>max)
		    {
		      max=y[j];
		    }
		  if(y[j]<min)
		    {
		      min=y[j];
		    }

		}
	    }

	  TH2F* hsum=new TH2F("hsum","",100,xmin-0.5,xmax+0.5,100,min-0.02,max+0.02);
	  hsum->Draw();
	  hsum->SetStats(kFALSE); 
   for(int j=0;j<nmax;j++)
    {
      if(g[j]!=0)
	{
	  g[j]->Draw("pl");
	}
	}
	}
    }
}

