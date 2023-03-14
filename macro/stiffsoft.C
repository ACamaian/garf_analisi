#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLegend.h>

void tz(TH1D *h,TGraphErrors *g[27],TGraphErrors *gt);
void stiffsoft()
{
  char tit[3][200]={"totFAZIAPRE_4840.root","totFAZIAPRE_4840_amd_geo_stiff.root","totFAZIAPRE_4840_amd_geo_soft.root"};

  TH2F *h[3];
  TFile *f[3];
  TGraphErrors *gg[3][27];
  TGraphErrors *ggr[3][27];
   TGraphErrors *ggs[3][27];
   TGraphErrors *gt[3];
   TH1F *hnev[3];
 int rif[27];
  rif[1]=1;
  rif[2]=4;
  rif[3]=7;
  rif[4]=9;
  rif[5]=10;
  rif[6]=12;
  rif[7]=15;
  rif[8]=16;
  TH1D *px[3];
  for(int j=0;j<3;j++)
    {
      f[j]=new TFile(tit[j]);
      h[j]=(TH2F*)f[j]->Get("hivalmolt");
      hnev[j]=(TH1F*)f[j]->Get("hneventi");

      px[j]=(TH1D*)h[j]->ProjectionX(Form("px%d",j),0,-1);
	  gt[j]=new TGraphErrors();
	  gt[j]->SetMarkerStyle(20);
	  gt[j]->SetMarkerColor(j+1);
      for(int iz=0;iz<27;iz++)
	{
	  gg[j][iz]=new TGraphErrors();
	  gg[j][iz]->SetMarkerStyle(20);
	  gg[j][iz]->SetMarkerColor(j+1);
	  ggr[j][iz]=new TGraphErrors();
	  ggr[j][iz]->SetMarkerStyle(20);
	  ggr[j][iz]->SetMarkerColor(j+1);
	  ggs[j][iz]=new TGraphErrors();
	  ggs[j][iz]->SetMarkerStyle(20);
	  ggs[j][iz]->SetMarkerColor(j+1);
	}
      tz(px[j],gg[j],gt[j]);
	
	  for(int iz=1;iz<9;iz++)
	    {
	      double *x,*y,*ey;
	      if(gg[j][iz]->GetN()>0)
		{
		  x=gg[j][iz]->GetX();
		  y=gg[j][iz]->GetY();
		  ey=gg[j][iz]->GetEY();
		  int k0=-1;
		  for(int k=0;k<gg[j][iz]->GetN();k++)
		    {
		      if(x[k]==rif[iz])
			{
			  k0=k;
			  break;
			}
		    }
		  if(k0>=0)
		    {
		  for(int k=0;k<gg[j][iz]->GetN();k++)
		    {
	float	      rap=y[k]/y[k0];
		float      erap=rap*(ey[k]/y[j]+ey[k0]/y[k0]);
		
		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
		      ggr[j][iz]->SetPointError(ggr[j][iz]->GetN()-1,0.0001,erap);


		    }

		    }


		}

	    }


    }

  for(int iz=1;iz<9;iz++)
	    {
	      double *x[3],*y[3],*ey[3];
	      for(int j=0;j<3;j++)
		{
	      if(gg[j][iz]->GetN()>0)
		{
		  x[j]=gg[j][iz]->GetX();
		  y[j]=gg[j][iz]->GetY();
		  ey[j]=gg[j][iz]->GetEY();
		}
		}
		  for(int k=0;k<gg[0][iz]->GetN();k++)
		    {
		      float a0=x[0][k];

		      for(int j=1;j<3;j++)
			{
			  int kk0=-1;
			  if(gg[j][iz]->GetN()>0)
			    {
			      for(int kk=0;kk<gg[j][iz]->GetN();kk++)
				{
				  if(x[j][kk]==a0)
				    {
				      kk0=kk;
				      break;
				    }
				}


			    }
			  if(kk0>=0)
			    {
			      float rap=y[0][k]*hnev[j]->Integral(2,2)/(y[j][kk0]*hnev[0]->Integral(2,2));
			      float      erap=rap*(ey[0][k]/y[0][k]+ey[j][kk0]/y[j][kk0]+1/sqrt(hnev[j]->Integral(2,2))+1/sqrt(hnev[0]->Integral(2,2)));
		
		      ggs[j][iz]->SetPoint(ggs[j][iz]->GetN(),x[0][k],rap);
		      ggs[j][iz]->SetPointError(ggs[j][iz]->GetN()-1,0.0001,erap);


			    }
			}

		    }
	    }

    
		 
		 
		 
  TCanvas *c[27];
  TH2F *hf[27];
  float amin[27]={0,0.5,2.5,5.5,6,9,10,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float amax[27]={0,3.5,6.5,10,14,14,17,20,25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float ymin[27]={0,1e-1,1e-3,1e-2,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  TLegend *leg[27];
  for(int iz=1;iz<9;iz++)
    {
      c[iz]=new TCanvas(Form("c%d",iz),Form("c%d",iz));
      c[iz]->Draw();
      hf[iz]=new TH2F(Form("hf%d",iz),Form("Z=%d",iz),100,amin[iz],amax[iz],100,ymin[iz],5);
      hf[iz]->GetXaxis()->SetTitle("A");
      hf[iz]->GetYaxis()->SetTitle("Y(A)/Y(Aref)");
      hf[iz]->SetStats(kFALSE);
      hf[iz]->Draw();
      ggr[0][iz]->Draw("p");
      ggr[1][iz]->Draw("p");
      ggr[2][iz]->Draw("p");
      gPad->SetLogy(1);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      leg[iz]=new TLegend(0.3,0.9,0.8,0.99);
      leg[iz]->AddEntry(ggr[0][iz],"exp","p");
      leg[iz]->AddEntry(ggr[1][iz],"STIFF","p");
      leg[iz]->AddEntry(ggr[2][iz],"SOFT","p");
      leg[iz]->Draw();
    }

  TCanvas *cmod[27];
  TH2F *hmod[27];

  float ymin2[27]={0,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float ymax2[27]={0,1e2,1e2,1e2,1e2,1e2,1e2,1e2,1e2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
TLegend *leg2[27];
  for(int iz=1;iz<9;iz++)
    {
      cmod[iz]=new TCanvas(Form("cmod%d",iz),Form("cmod%d",iz));
      cmod[iz]->Draw();
      hmod[iz]=new TH2F(Form("hmod%d",iz),Form("Z=%d",iz),100,amin[iz],amax[iz],100,ymin2[iz],ymax2[iz]);
      hmod[iz]->GetXaxis()->SetTitle("A");
hmod[iz]->GetYaxis()->SetTitle("Y(A)EXP/Y(A)SIM");
      hmod[iz]->SetStats(kFALSE);
      hmod[iz]->Draw();
      
      ggs[1][iz]->Draw("p");
      ggs[2][iz]->Draw("p");
      gPad->SetLogy(1);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      leg2[iz]=new TLegend(0.3,0.9,0.8,0.99);
      
      leg2[iz]->AddEntry(ggs[1][iz],"exp/STIFF","p");
      leg2[iz]->AddEntry(ggs[2][iz],"exp/SOFT","p");
      leg2[iz]->Draw();
    }

  
}
void tz(TH1D *h,TGraphErrors *g[27],TGraphErrors *gt)
{
  int totz[100];
  int totzn[100][100];


  
      for(int iz=0;iz<100;iz++)
	{
	      totz[iz]=0;

	  for(int ia=0;ia<100;ia++)
	    {
	      totzn[iz][ia]=0;
	 
	    }
	}
    

  for(int k=0;k<h->GetXaxis()->GetNbins();k++)
       {
	 int icount=h->GetBinContent(k+1);
	
	 if(icount>0)
	   {

	     int iz=(k+1)/100;
	     int in=k+1-iz*100-1;
	   
	            	 totz[iz]+=icount;
			 totzn[iz][in]+=icount; 
	   }
       }
  for(int iz=1;iz<26;iz++)
    {
      if(totz[iz]>0)
	{
	  gt->SetPoint(gt->GetN(),(float)iz,totz[iz]);
	  gt->SetPointError(gt->GetN()-1,0.001,sqrt((float)totz[iz]));

	  for(int in=0;in<100;in++)
	    {
	      if(totzn[iz][in]>0)
		{

		  g[iz]->SetPoint(g[iz]->GetN(),(float)(iz+in),(float)totzn[iz][in]);
		  g[iz]->SetPointError(g[iz]->GetN()-1,0.0001,sqrt((float)totzn[iz][in]));
  }
  }
	}
    }



}
