#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2);
void nsuz0(TH1F *h,TGraphErrors *gnz);
void cfr()
{
  gStyle->SetGridWidth(0);

  int col[3]={1,2,2};
  int mark[3]={20,20,24};
  int marke[2]={20,24};
  TFile *f[3];
  char tit[3][100]={"4025","4825","4840"};
  float ebeam[3]={25,25,40};
  float ap[3]={40,48,48};
  float at=12.;

  float vplab[3];
  float vpcm[3];
  for(int j=0;j<3;j++)
    {
      vplab[j]=300*sqrt(2*ebeam[j]/931.5);
      vpcm[j]=vplab[j]-ap[j]*vplab[j]/(ap[j]+at);
      cout<<vplab[j]<<" "<<vpcm[j]<<endl;
    }

  TH2F*hnzbigvl[3];
  TH2F*hnzbigvcm[3];
 
  int ny;

  for(int j=0;j<3;j++)
    {
            f[j]=new TFile(Form("feb5/totFAZIAPRE_%s.root",tit[j]));
      //f[j]=new TFile(Form("/home/piantell/prove/test_2020/%s.root",tit[j]));
      

hnzbigvl[j]=(TH2F*)f[j]->Get("hnzbigvlnorm");
hnzbigvcm[j]=(TH2F*)f[j]->Get("hnzbigvcmnorm");
cout<<j<<" "<<hnzbigvl[j]<<" "<<hnzbigvcm[j]<<endl;
    }

  ny=hnzbigvl[0]->GetYaxis()->GetNbins();
  //  int div=10;
  int div=2;
  ny=ny/div;
  TH1D *hl[3][ny];
  TH1D *hcm[3][ny];
  float vl[3][ny],vcm[3][ny];
  for(int j=0;j<3;j++)
    {
      for(int k=0;k<ny;k++)
	{
	  hl[j][k]=(TH1D*)hnzbigvl[j]->ProjectionX(Form("hl_%d_%d",j,k),k*div+1,(k+1)*div);
	  vl[j][k]=(hnzbigvl[j]->GetYaxis()->GetBinCenter(k*div+1)+hnzbigvl[j]->GetYaxis()->GetBinCenter((k+1)*div))/2;
	  hcm[j][k]=(TH1D*)hnzbigvcm[j]->ProjectionX(Form("hcm_%d_%d",j,k),k*div+1,(k+1)*div);
	  vcm[j][k]=(hnzbigvcm[j]->GetYaxis()->GetBinCenter(k*div+1)+hnzbigvcm[j]->GetYaxis()->GetBinCenter((k+1)*div))/2;
	  //	  cout<<hnzbigvcm[j]->GetYaxis()->GetBinCenter(k*div+1)<<" "<<hnzbigvcm[j]->GetYaxis()->GetBinCenter((k+1)*div)<<" "<<k*div+1<<" "<<(k+1)*div<<endl;
	}
    }


 

 //  cout<<ny<<" "<<endl;
  TGraphErrors *gl[3][ny];
  TGraphErrors *gcm[3][ny];
  int igoodl[3][ny];
  for(int j=0;j<3;j++)
    {
  for(int k=0;k<ny;k++)
    {
      gl[j][k]=0;
      igoodl[j][k]=0;

      if(hl[j][k]->Integral()>0)
	{

	  gl[j][k]=new TGraphErrors();
	  nsuz0((TH1F*)hl[j][k],gl[j][k]);
	  igoodl[j][k]=1;
	  cout<<"l "<<k<<" "<<vl[j][k]<<endl;
	  gl[j][k]->SetMarkerStyle(mark[j]);
	  gl[j][k]->SetMarkerColor(col[j]);
	  gl[j][k]->SetLineColor(col[j]);

	}
 gcm[j][k]=0;
      if(hcm[j][k]->Integral()>0)
	{
     

  gcm[j][k]=new TGraphErrors();
	  nsuz0((TH1F*)hcm[j][k],gcm[j][k]);
	  gcm[j][k]->SetMarkerStyle(mark[j]);
	  gcm[j][k]->SetMarkerColor(col[j]);
    gcm[j][k]->SetLineColor(col[j]);
	  //  cout<<"cm "<<k<<" "<<vcm[j][k]<<endl;
	}


    }
    }

  

  int imin[3],imax[3];
      for(int j=0;j<3;j++)
	{
	  for(int i=0;i<ny;i++)
	    {
	      if(igoodl[j][i]>0)
		{
		  imin[j]=i;
		  break;
		}
	    }
	  for(int i=ny-1;i>=0;i--)
	    {
	      if(igoodl[j][i]>0)
		{
		  imax[j]=i;
		  break;
		}
	    }

	}
      int Imin,Imax,Jmin,Jmax;
      Imin=1000;
      Imax=0;
      for(int j=0;j<3;j++)
	{
	  for(int i=imin[j];i<=imax[j];i++)
	    {
	     
	      if(Imin>imin[j]&&igoodl[j][i]>0)
		{
		  
		  Imin=imin[j];
		  Jmin=j;
		}
	     
	      if(Imax<imax[j]&&igoodl[j][i]>0)
		{
		  
		  Jmax=j;
		  Imax=imax[j];
		}

	    }
	}



       cout<<Imin<<" "<<Jmin<<" "<<Imax<<" "<<Jmax<<" "<<vl[Jmin][Imin]<<" "<<vl[Jmax][Imax]<<endl;
   TCanvas *cl=new TCanvas("cl","cl",1000,1000);
   cl->Divide(3,4);
   cl->Draw();
   TH2F *hkl[20];
   TLegend *leg=new TLegend(0.5,0.5,0.9,0.9);
  

   int n0=0;
  for(int k=Imin;k<=Imax;k++)
    {
      if(n0>12)continue;
      cl->cd(k-Imin+1);
hkl[k-Imin]=new TH2F(Form("hkl_%d",k),Form("vlab/vp=%f",vl[0][k]),28,9,28,100,0.8,1.3);
hkl[k-Imin]->Draw();
      hkl[k-Imin]->GetXaxis()->SetTitle("Z");
      hkl[k-Imin]->GetYaxis()->SetTitle("N/Z");
      for(int j=0;j<3;j++)
	{
	  if(igoodl[j][k]>0)
	    {
	      if(gl[j][k]!=0)
	    {
	  gl[j][k]->Draw("pl");
	    }
	    }



	}
      n0++;
    }
  for(int j=0;j<3;j++)
    {
      for(int k=0;k<ny;k++)
	{
	  if(gl[j][k]!=0)
	    {
	      leg->AddEntry(gl[j][k],tit[j],"p");
	      break;
	    }
	}
    }
  cl->cd(12);
  leg->Draw();
  cl->Print("gnzvl.pdf");

}

void nsuz0(TH1F *h,TGraphErrors *gnz)
{
  int totz[100];
  int totzn[100][100];
float nzval[100];
 float 	 enzval[100];
 float	 err[100];
  
      for(int iz=0;iz<100;iz++)
	{
	      totz[iz]=0;
nzval[iz]=0;
 	 enzval[iz]=0;
 	 err[iz]=0;
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


 for(int iz=0;iz<100;iz++)
    {
      if(totz[iz]>10)//ridurre per aumentare la statistica
	{

 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in]>0)
 	       {
		 		
 		 nzval[iz]=nzval[iz]+totzn[iz][in]*(float)in;
		
 	       }
 	   }
 	 if(nzval[iz]>0)
 	   {
 	     float nsuz=nzval[iz]/totz[iz];
 	     for(int in=0;in<100;in++)
 	       {
 		 if(totzn[iz][in]>0)
 		   {

 			 enzval[iz]=enzval[iz]+totzn[iz][in]*(float)in*(float)in;


 		   }
 	       }

 	     enzval[iz]=enzval[iz]/totz[iz];
 	     float nzvaldue=nsuz*nsuz;
 	     enzval[iz]=sqrt((enzval[iz]-nzvaldue)/totz[iz]);


 	      err[iz]=enzval[iz]/((float)iz);
	      
	     

 	     nsuz=nsuz/(float)iz;
	     if(gnz!=0)
	       {
 	 gnz->SetPoint(gnz->GetN(),(float)iz,nsuz);
 	 gnz->SetPointError(gnz->GetN()-1,0.001,err[iz]);
	       }
 	   }
		 		 
	}//totz>10
    }//iz=1,80
 
 cout<<gnz->GetName()<<endl;
 gnz->GetXaxis()->SetTitle("Z");
 gnz->GetYaxis()->SetTitle("N/Z");

}



void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100])
{
  int totz[100];
  int totzn[100][100];
float nzval[100];
 float 	 enzval[100];
 float	 err[100];
  
      for(int iz=0;iz<100;iz++)
	{
	      totz[iz]=0;
nzval[iz]=0;
 	 enzval[iz]=0;
 	 err[iz]=0;
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


 for(int iz=0;iz<100;iz++)
    {
      if(totz[iz]>10)//ridurre per aumentare la statistica
	{

 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in]>0)
 	       {
		 if(gz[iz]!=0)
		   {
 		 gz[iz]->SetPoint(gz[iz]->GetN(),(float)in,(float)totzn[iz][in]/(float)totz[iz]);
 		 gz[iz]->SetPointError(gz[iz]->GetN()-1,0.0001,((float)totzn[iz][in]/(float)totz[iz])*(sqrt((float)totzn[iz][in])/(float)totzn[iz][in]+sqrt((float)totz[iz])/(float)totz[iz]));
	

		   }		
 		 nzval[iz]=nzval[iz]+totzn[iz][in]*(float)in;
		
 	       }
 	   }
 	 if(nzval[iz]>0)
 	   {
 	     float nsuz=nzval[iz]/totz[iz];
 	     for(int in=0;in<100;in++)
 	       {
 		 if(totzn[iz][in]>0)
 		   {

 			 enzval[iz]=enzval[iz]+totzn[iz][in]*(float)in*(float)in;


 		   }
 	       }

 	     enzval[iz]=enzval[iz]/totz[iz];
 	     float nzvaldue=nsuz*nsuz;
 	     enzval[iz]=sqrt((enzval[iz]-nzvaldue)/totz[iz]);


 	      err[iz]=enzval[iz]/((float)iz);
	      
	     

 	     nsuz=nsuz/(float)iz;
	     if(gnz!=0)
	       {
 	 gnz->SetPoint(gnz->GetN(),(float)iz,nsuz);
 	 gnz->SetPointError(gnz->GetN()-1,0.001,err[iz]);
	       }
 	   }
		 if(gz[iz]!=0)
		   {
 gz[iz]->GetXaxis()->SetTitle("N");		
     gz[iz]->SetTitle(Form("Z=%d",iz));
		
		   }		 
	}//totz>10
    }//iz=1,80
 
 cout<<gnz->GetName()<<endl;
 gnz->GetXaxis()->SetTitle("Z");
 gnz->GetYaxis()->SetTitle("N/Z");

}
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2)
{
  int totz[100];
  int totzn[100][100];
float nzval[100];
 
     for(int iz=0;iz<100;iz++)
	{
	      totz[iz]=0;
nzval[iz]=0;
 	 
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
  double*z;
  double *nmed;
  double *emed;
  z=gnz->GetX();
  nmed=gnz->GetY();
  emed=gnz->GetEY();
  double nnmed[100],eemed[100];
  float err,err2,somman,sumerr,ey,nzval0;
  for(int j=0;j<gnz->GetN();j++)
    {
   
      //      gnz->GetPoint(j,z,nmed);
      //int iz=(int)z;
      int iz=(int)z[j];
      nnmed[j]=nmed[j]*z[j];
      eemed[j]=emed[j]*z[j];
      err=0;
      err2=0;
      somman=0;
      sumerr=0;
      ey=0;
      nzval0=0;
 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in]>0)
	       {
		 nzval[iz]=nzval[iz]+totzn[iz][in]*pow((float)in-nnmed[j],2);
		 err=err+totzn[iz][in];
		 // err=err+sqrt(totzn[iz][in])*pow((float)in-nnmed[j],2)+2*((float)in-nnmed[j])*eemed[j]*totzn[iz][in];
		 
		 //err2=err2+sqrt(totzn[iz][in]);
		 
	       }
	   }
	 
	 // err2=err2/totz[iz];
	 //err=err/nzval[iz];
	 nzval[iz]=nzval[iz]/totz[iz];
	 nzval0=nzval[iz];
	 //ey=nzval0*(err+err2);
	 nzval[iz]=sqrt(nzval0);
	 //ey=0.5*ey/sqrt(nzval0);
	 ey=nzval[iz]/sqrt(2*(err-1));
	 g2->SetPoint(g2->GetN(),z[j],nzval[iz]);
	 g2->SetPointError(g2->GetN()-1,0.0001,ey);
    }

}
