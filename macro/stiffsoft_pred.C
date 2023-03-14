#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <TLatex.h>

void tz(TH1D *h,TGraphErrors *g[27],TGraphErrors *gt);
void stiffsoft_tutti(TString s="feb24_23/totFAZIAPRE_4840_amd_geo_stiff.root,feb24_23/totFAZIAPRE_4840_amd_geo_soft.root,feb24_23/totFAZIAPRE_4840.root,feb24_23/totFAZIASYM_484835_amd_geo_stiff.root,feb24_23/totFAZIASYM_484835_amd_geo_soft.root,feb24_23/totFAZIASYM_484835.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_stiff.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_soft.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32.root")
{
 
  TFile *f[100];
  char name[10000];
  sprintf(name,"%s",s.Data());
  
  int jl[100];
  jl[0]=0;
  int nf=1;
  for(int j=0;j<strlen(name);j++)
    {
      if(name[j]==',')
	{
	  jl[nf]=j+1;
	  nf++;
	}
    }
  cout<<nf<<endl;
  TString sf[100];
  TString stit[100];
  for(int j=0;j<nf;j++)
    {
      if(j<nf-1)
	{
      sf[j]=s(jl[j],jl[j+1]-jl[j]-1);
	}
      else
	{
	  sf[j]=s(jl[j],s.Length()-jl[j]);
	}
      cout<<j<<" "<<sf[j].Data()<<endl;

      //      f[j]=new TFile(Form("%s",sf[j].Data()));
      f[j]=new TFile(sf[j].Data());
      if(sf[j].Contains("/"))
	{
	  int kk=sf[j].Last('/');
	  stit[j]=sf[j](kk+1,sf[j].Length()-kk);
	}
      else
	{
	  stit[j]=sf[j];
	}
      cout<<stit[j].Data()<<endl;
    }

  int * col=(int*)malloc(nf*sizeof (int));
  int * sym=(int*)malloc(nf*sizeof (sym));
 
  char ** tit=(char**)malloc(nf*sizeof(char*));
  int * tipo=(int*)malloc(nf*sizeof (int));
 int * etipo=(int*)malloc(nf*sizeof (int));
 
 //  int * serv=(int*)malloc(nf*sizeof (int));
  for(int j=0;j<nf;j++)
    {
      tit[j]=(char*)malloc(100*sizeof(char));
      if(stit[j].Contains("64")&&stit[j].Contains("32")) //6464_32
	{
	  
	  col[j]=1;
	  if(stit[j].Contains("stiff"))
	    {
	      sym[j]=26;
	  
	  sprintf(tit[j],"6464_32ST");
	    }
	  else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
	  

	  sprintf(tit[j],"6464_32SO");
	    }
	  else
	    {
	  sym[j]=20;
	  
	  sprintf(tit[j],"6464_32E");
	  etipo[j]=0;
	    }
	  tipo[j]=0;

	}



      if(stit[j].Contains("64")&&stit[j].Contains("52")) //6464_52
	{

	  col[j]=4;
	  if(stit[j].Contains("stiff"))
	    {
	      sym[j]=26;
	     
	  sprintf(tit[j],"6464_52ST");
	    }
	  else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
	       
	  sprintf(tit[j],"6464_52SO");
	    }
	  else
	    {
	  sym[j]=20;
	
	  sprintf(tit[j],"6464_52E");
	    }
	  tipo[j]=1;
	}


      if(stit[j].Contains("FAZIASYM")) //4848_35
	{
	  
	  col[j]=2;
	  if(stit[j].Contains("stiff"))
	    {
	      sym[j]=26;
	     
	     sprintf(tit[j],"4848_35ST");
	     
	    }
	 else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
	     
	     sprintf(tit[j],"4848_35SO");
	    }
	 else
	   {
	     sym[j]=21;
	     
	     sprintf(tit[j],"4848_35E");
	   }
	  tipo[j]=2;
	}
      if(stit[j].Contains("FAZIAPRE")&&stit[j].Contains("4840")) //4812_40
	{
	
	  col[j]=kGreen+2;
	
	  if(stit[j].Contains("stiff"))
	    {
	      sym[j]=26;
 
     sprintf(tit[j],"4812_40ST");
	    }
	else  if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;

     sprintf(tit[j],"4812_40SO");
	    }
	else
	  {
  sym[j]=34;

     sprintf(tit[j],"4812_40E");
	  }

	  tipo[j]=3;
    }



     if(stit[j].Contains("FAZIAPRE")&&stit[j].Contains("4825")) //4812_25
	{
	
	  col[j]=6;
	
	  if(stit[j].Contains("stiff"))
	    {
	      sym[j]=26;

     sprintf(tit[j],"4812_40ST");
	    }
	else  if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
 
     sprintf(tit[j],"4812_40SO");
	    }
	else
	  {
  sym[j]=29;

     sprintf(tit[j],"4812_40E");
	  }

	  tipo[j]=4;
    }
    }





  TH2F ***h=(TH2F***)malloc(nf*sizeof(TH2F**));
 
  TH1F **hnev=(TH1F**)malloc(nf*sizeof(TH1F*));
  TH2F **hpz=(TH2F**)malloc(nf*sizeof(TH2F*));

  int rif[27];
   rif[1]=1;
   rif[2]=4;
   rif[3]=7;
   rif[4]=9;
   rif[5]=10;
   rif[6]=12;
   rif[7]=15;
   rif[8]=16;



   TH1D ***px=(TH1D***)malloc(nf*sizeof(TH1D**));
 
   TGraphErrors ***gt=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
   TGraphErrors ****gg=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
   TGraphErrors ****ggr=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
   TGraphErrors ****ggs=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
  TGraphErrors ****ggy=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
  TGraphErrors ****ggt=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
  TGraphErrors ****ggt2=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));

  TGraphErrors ****ggp=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));

 
   for(int j=0;j<nf;j++)
     {
       hnev[j]=(TH1F*)f[j]->Get("hneventi");
       h[j]=(TH2F**)malloc(10*sizeof(TH2F*));
       hpz[j]=(TH2F*)f[j]->Get("hpz");

       px[j]=(TH1D**)malloc(10*sizeof(TH1D*));
       gt[j]=(TGraphErrors**)malloc(10*sizeof(TGraphErrors*));
	  gg[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
	  ggr[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
	  ggs[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
	  ggy[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
	  ggt[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
	  ggt2[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));

	  ggp[j]=(TGraphErrors***)malloc(50*sizeof(TGraphErrors**));

       for(int k=0;k<10;k++)
	 {
	   h[j][k]=(TH2F*)f[j]->Get(Form("h%d",100+k));
	   px[j][k]=(TH1D*)h[j][k]->ProjectionX(Form("px%d_%d",j,k),1,h[j][k]->GetYaxis()->FindBin(0.));
	   gt[j][k]=new TGraphErrors(); 
	   gt[j][k]->SetMarkerStyle(sym[j]);
	   gt[j][k]->SetMarkerColor(col[j]);
	  gg[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
	  ggr[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
	  ggs[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
	  ggy[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
	  ggt[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
	  ggt2[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));

       for(int iz=0;iz<27;iz++)
 	{
 	  gg[j][k][iz]=new TGraphErrors();
 	  gg[j][k][iz]->SetMarkerStyle(sym[j]);
 	  gg[j][k][iz]->SetMarkerColor(col[j]);
	  gg[j][k][iz]->SetLineColor(col[j]);
	  gg[j][k][iz]->SetMarkerSize(1.2);

 	  ggr[j][k][iz]=new TGraphErrors();
 	  ggr[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggr[j][k][iz]->SetMarkerColor(col[j]);
 	  ggr[j][k][iz]->SetLineColor(col[j]);
	  ggr[j][k][iz]->SetMarkerSize(1.2);

 	  ggs[j][k][iz]=new TGraphErrors();
 	  ggs[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggs[j][k][iz]->SetMarkerColor(col[j]);
	  ggs[j][k][iz]->SetLineColor(col[j]);
	  ggs[j][k][iz]->SetMarkerSize(1.2);

 	  ggy[j][k][iz]=new TGraphErrors();
 	  ggy[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggy[j][k][iz]->SetMarkerColor(col[j]);
	  ggy[j][k][iz]->SetLineColor(col[j]);
	  ggy[j][k][iz]->SetMarkerSize(1.2);

 	  ggt[j][k][iz]=new TGraphErrors();
 	  ggt[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggt[j][k][iz]->SetMarkerColor(col[j]);
	  ggt[j][k][iz]->SetLineColor(col[j]);
	  ggt[j][k][iz]->SetMarkerSize(1.2);

 	  ggt2[j][k][iz]=new TGraphErrors();
 	  ggt2[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggt2[j][k][iz]->SetMarkerColor(col[j]);
	  ggt2[j][k][iz]->SetLineColor(col[j]);
	  ggt2[j][k][iz]->SetMarkerSize(1.2);

 	}

       tz(px[j][k],gg[j][k],gt[j][k]);
 
	 }
    for(int k=0;k<50;k++)
	 {
	   ggp[j][k]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
      for(int iz=0;iz<27;iz++)
 	{
 	  ggp[j][k][iz]=new TGraphErrors();
 	  ggp[j][k][iz]->SetMarkerStyle(sym[j]);
 	  ggp[j][k][iz]->SetMarkerColor(col[j]);
	  ggp[j][k][iz]->SetLineColor(col[j]);
	  ggp[j][k][iz]->SetMarkerSize(1.2);
	}
	 }
     }




 int *** ainf=(int***)malloc(nf*sizeof (int**));
  int *** asup=(int***)malloc(nf*sizeof (int**));
 
  for(int i=0;i<nf;i++)
    {
      ainf[i]=(int**)malloc(10*sizeof(int*));
      asup[i]=(int**)malloc(10*sizeof(int*));
      for(int k=0;k<10;k++)
	{
     ainf[i][k]=(int*)malloc(27*sizeof(int));
      asup[i][k]=(int*)malloc(27*sizeof(int));
 
	}
    }
    
  for(int j=0;j<nf;j++)
    {
      for(int k=0;k<10;k++)
	{
      for(int iz=0;iz<27;iz++)
	{
	  ainf[j][k][iz]=-1;
	  asup[j][k][iz]=-1;
	}

	}
    }
  
  for(int j=0;j<nf;j++)
    {
      for(int m=0;m<10;m++)
	{
     
      if(gg[j][m][5]->GetN()>0)
	{
	  
	  int i1=1;
	  while(i1==1)
	    {
	      i1=0;
	      double *x=gg[j][m][5]->GetX();
	  for(int k=0;k<gg[j][m][5]->GetN();k++)
	    {
	      if(x[k]<10)
		{
		 
		gg[j][m][5]->RemovePoint(k);
		i1=1;
		break;
		}
	    }
	    }
	}
    }	  
    }



    for(int j=0;j<nf;j++)
      {
        if(!stit[j].Contains("stiff")&&!stit[j].Contains("soft"))
 	 {
	 
      for(int m=0;m<10;m++)
	{

  for(int iz=0;iz<27;iz++)
	     {
	       if(gg[j][m][iz]->GetN()>0)
		 {
	
	   double *x=gg[j][m][iz]->GetX();
 	   ainf[j][m][iz]=1000;
 	   asup[j][m][iz]=-1000;
	   for(int i=0;i<gg[j][m][iz]->GetN();i++)
	     {
	       if(x[i]<ainf[j][m][iz])ainf[j][m][iz]=x[i];
	       if(x[i]>asup[j][m][iz])asup[j][m][iz]=x[i];
	     }
	   
 	 }
	 }
      }
      }
      }
    for(int j=0;j<nf;j++)
      {
	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
	  {
	    for(int m=0;m<10;m++)
	      {
	    for(int k=0;k<nf;k++)
	      {
		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
		  {
		for(int iz=0;iz<27;iz++)
	     {
	       ainf[j][m][iz]= ainf[k][m][iz];
	       asup[j][m][iz]= asup[k][m][iz];
	       	       
	       
	     }    

		    break;
		  }
	      }
	  }
      }
      }


  




  for(int j=0;j<nf;j++)
     {
  for(int m=0;m<10;m++)
	      {
 	  for(int iz=1;iz<9;iz++)
 	    {
 	      double *x,*y,*ey;
 	      if(gg[j][m][iz]->GetN()>0)
 		{
 		  x=gg[j][m][iz]->GetX();
 		  y=gg[j][m][iz]->GetY();
 		  ey=gg[j][m][iz]->GetEY();
 		  int k0=-1;
 		  for(int k=0;k<gg[j][m][iz]->GetN();k++)
 		    {
 		      if(x[k]==rif[iz])
 			{
 			  k0=k;
 			  break;
 			}
 		    }
 		  if(k0>=0) 		    
{
 		  for(int k=0;k<gg[j][m][iz]->GetN();k++)
 		    {

		      //		      if(iz==1 && x[k]>3)continue;
		      //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
		      if(x[k]<ainf[j][m][iz])continue;
		      if(x[k]>asup[j][m][iz])continue;

		      
		      float rap=y[k]/y[k0];
		      float erap=rap*(ey[k]/y[k]+ey[k0]/y[k0]);
		      
		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
 		      ggr[j][m][iz]->SetPoint(ggr[j][m][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
 		      ggr[j][m][iz]->SetPointError(ggr[j][m][iz]->GetN()-1,0.0001,erap);
		      

 		    }
 		    }


 		}

 	    }


     


     }
     }

 


   for(int j=0;j<nf;j++)
     {
  for(int m=0;m<10;m++)
	      {
       if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
	 {
       for(int iz=1;iz<9;iz++)
	 {
	   double *x,*y,*ey; //sim
 	      if(ggr[j][m][iz]->GetN()>0)
 		{
 		  x=ggr[j][m][iz]->GetX();
 		  y=ggr[j][m][iz]->GetY();
 		  ey=ggr[j][m][iz]->GetEY();
		  

       for(int iq=0;iq<nf;iq++)
	 {
	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
	     {
	   if(tipo[j]==tipo[iq])
	     {

	 
	       double *xq,*yq,*eyq; //exp
 	      if(ggr[iq][m][iz]->GetN()>0)
 		{
 		  xq=ggr[iq][m][iz]->GetX();
 		  yq=ggr[iq][m][iz]->GetY();
 		  eyq=ggr[iq][m][iz]->GetEY();
 		  
 		  for(int k=0;k<ggr[j][m][iz]->GetN();k++)
 		    {
		      int k0=-1;
		      for(int l=0;l<ggr[iq][m][iz]->GetN();l++)
			{

 		      if(x[k]==xq[l])
 			{
 			  k0=l;
 			  break;
 			}
			}
 		  if(k0>=0) 		    
		    {
		      if(x[k]<ainf[j][m][iz])continue;
		      if(x[k]>asup[j][m][iz])continue;
  
		      float rap=yq[k0]/y[k];
		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		
		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
 		      ggs[j][m][iz]->SetPoint(ggs[j][m][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
 		      ggs[j][m][iz]->SetPointError(ggs[j][m][iz]->GetN()-1,0.0001,erap);
		      
		      
 		    }
 		    }


 	

 	    }


	     }
	     }
	 }
		}
	 }
	 }
     } 



     }




   for(int j=0;j<nf;j++)
     {
       if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
	 {
  for(int m=0;m<10;m++)
	      {
       for(int iz=1;iz<9;iz++)
	 {
	   double *x,*y,*ey; //sim
 	      if(gg[j][m][iz]->GetN()>0)
 		{
 		  x=gg[j][m][iz]->GetX();
 		  y=gg[j][m][iz]->GetY();
 		  ey=gg[j][m][iz]->GetEY();
		  

       for(int iq=0;iq<nf;iq++)
	 {
	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
	     {
	   if(tipo[j]==tipo[iq])
	     {

	 
	       double *xq,*yq,*eyq; //exp
 	      if(gg[iq][m][iz]->GetN()>0)
 		{
 		  xq=gg[iq][m][iz]->GetX();
 		  yq=gg[iq][m][iz]->GetY();
 		  eyq=gg[iq][m][iz]->GetEY();
 		  
 		  for(int k=0;k<gg[j][m][iz]->GetN();k++)
 		    {
		      int k0=-1;
		      for(int l=0;l<gg[iq][m][iz]->GetN();l++)
			{
 		      if(x[k]==xq[l])
 			{
 			  k0=l;
 			  break;
 			}
			}
 		  if(k0>=0) 		    
		    {
		      if(x[k]<ainf[j][m][iz])continue;
		      if(x[k]>asup[j][m][iz])continue;
 		  
		      //		      float rap=yq[k0]*hnev[j]->Integral(102,102)/(y[k]*hnev[iq]->Integral(102,102));
		      float rap=yq[k0]*hpz[j]->Integral(hpz[j]->GetYaxis()->FindBin(0.1*k),hpz[j]->GetYaxis()->FindBin(0.1*(k+1))-1)/(y[k]*hpz[iq]->Integral(hpz[iq]->GetYaxis()->FindBin(0.1*k),hpz[iq]->GetYaxis()->FindBin(0.1*(k+1))-1));
		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]+1/sqrt(hpz[j]->Integral(hpz[j]->GetYaxis()->FindBin(0.1*k),hpz[j]->GetYaxis()->FindBin(0.1*(k+1))-1))+1/sqrt(hpz[iq]->Integral(hpz[iq]->GetYaxis()->FindBin(0.1*k),hpz[iq]->GetYaxis()->FindBin(0.1*(k+1))-1)));
		

 		      ggy[j][m][iz]->SetPoint(ggy[j][m][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
 		      ggy[j][m][iz]->SetPointError(ggy[j][m][iz]->GetN()-1,0.0001,erap);
		      

		    }
 		    }


 	

 	    }


	     }
	     }
	 }
		}
	 }
	 }
     } 

     }

   for(int j=0;j<nf;j++)
     {
       for(int m=0;m<10;m++)
	 {
	   for(int iz=1;iz<9;iz++)
	     {
	       double *x,*y,*ey;
	       if(gg[j][m][iz]->GetN()>0)
		 {
		   x=gg[j][m][iz]->GetX();
		   y=gg[j][m][iz]->GetY();
		   ey=gg[j][m][iz]->GetEY();
		   float tot=0;
		   for(int k=0;k<gg[j][m][iz]->GetN();k++)
		     {

		       //		      if(iz==1 && x[k]>3)continue;
		       //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
		       if(x[k]<ainf[j][m][iz])continue;
		       if(x[k]>asup[j][m][iz])continue;
		       tot=tot+y[k];
		     }
		   for(int k=0;k<gg[j][m][iz]->GetN();k++)
		     {
		       if(x[k]<ainf[j][m][iz])continue;
		       if(x[k]>asup[j][m][iz])continue;

		      
		       float rap=y[k]/tot;
		       float erap=rap*(ey[k]/y[k]+1/sqrt(tot));
		
		       // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
		       ggt[j][m][iz]->SetPoint(ggt[j][m][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
		       ggt[j][m][iz]->SetPointError(ggt[j][m][iz]->GetN()-1,0.0001,erap);
		      

		     }
 		    


		 }

	     }


     


	 }
     }


   int iamin[27],iamax[27];
   for(int iz=0;iz<27;iz++)
     {
       iamax[iz]=-1;
       iamin[iz]=1000;

     }


   for(int j=0;j<nf;j++)
     {
       for(int m=0;m<10;m++)
	 {
	  
	   for(int iz=1;iz<9;iz++)
	     {
	      
	       double *x,*y,*ey;
		   x=ggt[j][m][iz]->GetX();
		   y=ggt[j][m][iz]->GetY();
		   ey=ggt[j][m][iz]->GetEY();
	       
	       if(ggt[j][m][iz]->GetN()>0)
		 {
		   for(int k=0;k<ggt[j][m][iz]->GetN();k++)
		     {
		       int ia=x[k];
		       if(iamin[iz]>ia)iamin[iz]=ia;
		       if(iamax[iz]<ia)iamax[iz]=ia;
		       cout<<j<<" "<<ia<<" "<<iz<<" "<<m<<" "<<y[k]<<endl;
		       ggp[j][ia][iz]->SetPoint(ggp[j][ia][iz]->GetN(),(float)(m+0.5),y[k]);
		       ggp[j][ia][iz]->SetPointError(ggp[j][ia][iz]->GetN()-1,0.0001,ey[k]);

		     }
		 }

	     }


	 }
     }
  



   for(int j=0;j<nf;j++)
     {
       if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
	 {
	   for(int m=0;m<10;m++)
	     {
	       for(int iz=1;iz<9;iz++)
		 {
		   double *x,*y,*ey; //sim
		   if(ggt[j][m][iz]->GetN()>0)
		     {
		       x=ggt[j][m][iz]->GetX();
		       y=ggt[j][m][iz]->GetY();
		       ey=ggt[j][m][iz]->GetEY();
		  

		       for(int iq=0;iq<nf;iq++)
			 {
			   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
			     {
			       if(tipo[j]==tipo[iq])
				 {

	 
				   double *xq,*yq,*eyq; //exp
				   if(ggt[iq][m][iz]->GetN()>0)
				     {
				       xq=ggt[iq][m][iz]->GetX();
				       yq=ggt[iq][m][iz]->GetY();
				       eyq=ggt[iq][m][iz]->GetEY();
 		  
				       for(int k=0;k<ggt[j][m][iz]->GetN();k++)
					 {
					   int k0=-1;
					   for(int l=0;l<ggt[iq][m][iz]->GetN();l++)
					     {
					       if(x[k]==xq[l])
						 {
						   k0=l;
						   break;
						 }
					     }
					   if(k0>=0) 		    
					     {
 		  
					       float rap=yq[k0]/y[k];
					       float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		

					       ggt2[j][m][iz]->SetPoint(ggt2[j][m][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
					       ggt2[j][m][iz]->SetPointError(ggt2[j][m][iz]->GetN()-1,0.0001,erap);
		      

					     }
					 }


 	

				     }


				 }
			     }
			 }
		     }
		 }
	     }
	 } 


     }


   TCanvas **cpz=(TCanvas**)malloc(nf*sizeof(TCanvas*));
   for(int j=0;j<nf;j++)
     {
       cpz[j]=new TCanvas(Form("cpz_%d",j),Form("%s",tit[j]));
       cpz[j]->Draw();
       hpz[j]->SetTitle(Form("%s",tit[j]));
       hpz[j]->Draw("zcol");
       gPad->SetLogz(kTRUE);
     }




   TLine *l=new TLine(0,1,30,1);
   l->SetLineStyle(2);
   l->SetLineWidth(2);

   l->SetLineColor(kBlue-7);
   
   TCanvas *c[27];
   TH2F *hf[27][10];
   float amin[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float amax[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float ymin[27]={0,1e-1,1e-5,1e-3,1e-5,1e-5,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *leg[27];
   
   for(int iz=1;iz<9;iz++)
     {

       c[iz]=new TCanvas(Form("c%d",iz),Form("c%d",iz),0,0,1000,1000);
       c[iz]->Divide(3,3);
       c[iz]->Draw();
	   leg[iz]=new TLegend(0.62,0.49,0.97,0.98);
       for(int m=1;m<10;m++)
	 {
	   c[iz]->cd(m);
	   hf[iz][m]=new TH2F(Form("hf%d_%d",iz,m),Form("Z=%d pred=%d -%d",iz,m,m+1),100,amin[iz],amax[iz],100,ymin[iz],10);
       hf[iz][m]->GetXaxis()->SetTitle("A");
       hf[iz][m]->GetYaxis()->SetTitle(Form("Y(A)/Y(%d)",rif[iz]));
       hf[iz][m]->SetStats(kFALSE);
       hf[iz][m]->Draw();
	   

       for(int k=0;k<nf;k++)
	 {
       ggr[k][m][iz]->Draw("pl");
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);


       
       if(m==1)
	 {

	   leg[iz]->AddEntry(ggr[k][m][iz],tit[k],"p");
	
	 }

	 }

       if(m==1)leg[iz]->Draw();

	 }
     }


   TCanvas *cb[27];
   TH2F *hfb[27][10];
   float aminb[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float amaxb[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   float yminb[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *legb[27];
   for(int iz=1;iz<9;iz++)
     {

       cb[iz]=new TCanvas(Form("cb%d",iz),Form("cb%d",iz),0,0,1000,1000);
       cb[iz]->Divide(3,3);
       cb[iz]->Draw();
       legb[iz]=new TLegend(0.62,0.49,0.97,0.98);
       for(int m=1;m<10;m++)
	 {
	   cb[iz]->cd(m);      
	   hfb[iz][m]=new TH2F(Form("hfb%d_%d",iz,m),Form("Z=%d pred=%d -%d",iz,m,m+1),100,aminb[iz],amaxb[iz],100,yminb[iz],100);
       hfb[iz][m]->GetXaxis()->SetTitle("A");
       hfb[iz][m]->GetYaxis()->SetTitle("(Y/Yrif)EXP/(Y/Yrif)Sim");
       hfb[iz][m]->SetStats(kFALSE);
       hfb[iz][m]->Draw();
       l->Draw();

       for(int k=0;k<nf;k++)
	 {
	   if(ggs[k][m][iz]->GetN()>0)
	     {
       ggs[k][m][iz]->Draw("p");

	     }
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);
       if(m==1)
	 {     
       legb[iz]->AddEntry(ggs[k][m][iz],tit[k],"p");
	     }
	     }
	 
           if(m==1) legb[iz]->Draw();
	 }
     }

  TCanvas *cy[27];
   TH2F *hfy[27][10];
   float aminy[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float amaxy[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   float yminy[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *legy[27];
   for(int iz=1;iz<9;iz++)
     {

       cy[iz]=new TCanvas(Form("cy%d",iz),Form("cy%d",iz),0,0,1000,1000);
       cy[iz]->Divide(3,3);
       cy[iz]->Draw();
       legy[iz]=new TLegend(0.62,0.49,0.97,0.98);
      for(int m=1;m<10;m++)
	 {
	   cy[iz]->cd(m);   
	   hfy[iz][m]=new TH2F(Form("hfy%d_%d",iz,m),Form("Z=%d pred=%d -%d",iz,m,m+1),100,aminy[iz],amaxy[iz],100,yminy[iz],100);
       hfy[iz][m]->GetXaxis()->SetTitle("A");
       hfy[iz][m]->GetYaxis()->SetTitle("YEXP/YSIM (norm neventi)");
       hfy[iz][m]->SetStats(kFALSE);
       hfy[iz][m]->Draw();
       l->Draw();

       for(int k=0;k<nf;k++)
	 {
	   if(ggy[k][m][iz]->GetN()>0)
	     {
       ggy[k][m][iz]->Draw("p");
	     }
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);
	     
	   if(m==1)
	     {
       legy[iz]->AddEntry(ggy[k][m][iz],tit[k],"p");
	     }
	 }

	 
      if(m==1)legy[iz]->Draw();
	 }
     }




  TCanvas *ct[27];
   TH2F *hft[27][10];
   float amint[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float amaxt[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float ymint[27]={0,1e-1,1e-5,1e-3,1e-5,1e-6,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *legt[27];
   for(int iz=1;iz<9;iz++)
     {

       ct[iz]=new TCanvas(Form("ct%d",iz),Form("ct%d",iz),0,0,1000,1000);
        ct[iz]->Divide(3,3);
       ct[iz]->Draw();
       legt[iz]=new TLegend(0.62,0.49,0.97,0.98);
      for(int m=1;m<10;m++)
	 {   
	   ct[iz]->cd(m);        
	   hft[iz][m]=new TH2F(Form("hft%d_%d",iz,m),Form("Z=%d pred=%d -%d",iz,m,m+1),100,amint[iz],amaxt[iz],100,ymint[iz],10);
       hft[iz][m]->GetXaxis()->SetTitle("A");
       hft[iz][m]->GetYaxis()->SetTitle("Y(A)(norm 1)");
       hft[iz][m]->SetStats(kFALSE);
       hft[iz][m]->Draw();
       
       for(int k=0;k<nf;k++)
	 {
	   if(!stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
	     {
       ggt[k][m][iz]->Draw("pl");
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);
       
       if(m==1)
	 {
       legt[iz]->AddEntry(ggt[k][m][iz],tit[k],"p");
	     }
	     }
	 }
       if(m==1)legt[iz]->Draw();
	 }
     }

  TCanvas *ct2[27];
   TH2F *hft2[27][10];
   float amint2[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float amaxt2[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

   float ymint2[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *legt2[27];
   for(int iz=1;iz<9;iz++)
     {

       ct2[iz]=new TCanvas(Form("ct2_%d",iz),Form("ct2_%d",iz),0,0,1000,1000);
       ct2[iz]->Divide(3,3);
       ct2[iz]->Draw();
       legt2[iz]=new TLegend(0.62,0.49,0.97,0.98);
       for(int m=1;m<10;m++)
	 {
	   ct2[iz]->cd(m);
	   hft2[iz][m]=new TH2F(Form("hft2_%d_%d",iz,m),Form("Z=%d pred=%d -%d",iz,m,m+1),100,amint2[iz],amaxt2[iz],100,ymint2[iz],100);
       hft2[iz][m]->GetXaxis()->SetTitle("A");
       hft2[iz][m]->GetYaxis()->SetTitle("YEXP/YSIM (norm1)");
       hft2[iz][m]->SetStats(kFALSE);
       hft2[iz][m]->Draw();
       l->Draw();

       for(int k=0;k<nf;k++)
	 {
	   if(ggt2[k][m][iz]->GetN()>0)
	     {
       ggt2[k][m][iz]->Draw("p");
	     }
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);
	     
	   if(m==1)
	     {
       legt2[iz]->AddEntry(ggy[k][m][iz],tit[k],"p");
	     }
	 }
       if(m==1)legt2[iz]->Draw();
	 }
     }






  TCanvas *cp[27];
   TH2F *hfp[27][50];
   //   float aminp[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   //float amaxp[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   float yminp[27]={0,1e-1,1e-5,1e-3,1e-5,1e-6,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   TLegend *legp[27];
   for(int iz=1;iz<9;iz++)
     {

       cp[iz]=new TCanvas(Form("cp%d",iz),Form("cp%d",iz),0,0,1000,1000);
        cp[iz]->Divide(3,3);
       cp[iz]->Draw();
       legp[iz]=new TLegend(0.62,0.49,0.97,0.98);
       int nia=0;
       for(int ia=iamin[iz];ia<=iamax[iz];ia++)
	 {   

	   nia++;
	   if(nia<10)
	     {
	   cp[iz]->cd(nia);        
	   hfp[iz][nia-1]=new TH2F(Form("hfp%d_%d",iz,ia),Form("Z=%d A=%d",iz,ia),110,0,11,100,yminp[iz],10);
       hfp[iz][nia-1]->GetXaxis()->SetTitle("pred");
       hfp[iz][nia-1]->GetYaxis()->SetTitle("Y(norm 1)");
       hfp[iz][nia-1]->SetStats(kFALSE);
       hfp[iz][nia-1]->Draw();
       
       for(int k=0;k<nf;k++)
	 {
	   if(!stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
	     {
	      
	       if(ggp[k][ia][iz]->GetN()>0)ggp[k][ia][iz]->Draw("pl");
       gPad->SetLogy(1);
       gPad->SetGridx(kFALSE);
       gPad->SetGridy(kFALSE);
       
       if(nia==1)
	 {
       legp[iz]->AddEntry(ggp[k][ia][iz],tit[k],"p");
	     }
	     }
	 }
       if(nia==1)legp[iz]->Draw();
	     }
	   else
	     {
	       cout<<"ci sono altri isotopi "<<" Z="<<iz<<" A="<<ia<<endl;
	     }
	 }
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
