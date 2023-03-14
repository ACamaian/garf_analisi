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
float sum(TH1D *p);
void stiffsoft_tutti(TString s="mar03_23/totFAZIAPRE_4840_amd_geo_stiff.root,mar03_23/totFAZIAPRE_4840_amd_geo_soft.root,mar03_23/totFAZIAPRE_4840.root,mar03_23/totFAZIASYM_484835_amd_geo_stiff.root,mar03_23/totFAZIASYM_484835_amd_geo_soft.root,mar03_23/totFAZIASYM_484835.root,mar03_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_stiff.root,mar03_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_soft.root,mar03_23/totINDRAFAZIA_64Ni64Ni_32.root")
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





  TH2F **h=(TH2F**)malloc(nf*sizeof(TH2F*));
  TH2F **hsub=(TH2F**)malloc(nf*sizeof(TH2F*));
  TH1F **hnev=(TH1F**)malloc(nf*sizeof(TH1F*));

  int rif[27];
   rif[1]=1;
   rif[2]=4;
   rif[3]=7;
   rif[4]=9;
   rif[5]=10;
   rif[6]=12;
   rif[7]=15;
   rif[8]=16;



   TH1D **px=(TH1D**)malloc(nf*sizeof(TH1D*));
   TH1D **pn=(TH1D**)malloc(nf*sizeof(TH1D*));
   TH1D **pnz=(TH1D**)malloc(nf*sizeof(TH1D*));
   TH1D **pz=(TH1D**)malloc(nf*sizeof(TH1D*));
   TH1D **pxsub=(TH1D**)malloc(nf*sizeof(TH1D*));
   TGraphErrors **gt=(TGraphErrors**)malloc(nf*sizeof(TGraphErrors*));
   TGraphErrors ***gg=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
   TGraphErrors ***ggr=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
   TGraphErrors ***ggs=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggy=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggt=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggt2=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));


   TGraphErrors **gtsub=(TGraphErrors**)malloc(nf*sizeof(TGraphErrors*));
   TGraphErrors ***ggsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
   TGraphErrors ***ggrsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
   TGraphErrors ***ggssub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggysub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggtsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
  TGraphErrors ***ggt2sub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));

float *** tot=(float ***)malloc(nf*sizeof (float**));
float **** tia=(float ****)malloc(nf*sizeof (float***));

 for(int j=0;j<nf;j++)
   {
     tot[j]=(float**)malloc(60*sizeof(float*));
     tia[j]=(float***)malloc(60*sizeof(float**));
     for(int iz=0;iz<60;iz++)
       {
	 tot[j][iz]=(float*)malloc(10*sizeof(float*));
	 tia[j][iz]=(float**)malloc(10*sizeof(float**));
     for(int ias=0;ias<10;ias++)
       {
	 tot[j][iz][ias]=0;
	 tia[j][iz][ias]=(float*)malloc(100*sizeof(float*));
     for(int ia=0;ia<100;ia++)
       {
	 tia[j][iz][ias][ia]=0;
       }
       }
       }
   }


   for(int j=0;j<nf;j++)
     {
     
       h[j]=(TH2F*)f[j]->Get("hpz");
       hnev[j]=(TH1F*)f[j]->Get("hneventi");
       
       cout<<j<<endl;
       hsub[j]=(TH2F*)f[j]->Get("hivalvpar");

       px[j]=(TH1D*)h[j]->ProjectionX(Form("px%d",j),1,h[j]->GetYaxis()->FindBin(0.));
       pxsub[j]=(TH1D*)hsub[j]->ProjectionX(Form("pxsub%d",j),1,hsub[j]->GetYaxis()->GetNbins());

       pz[j]=(TH1D*)px[j]->Clone(Form("pz_%d",j));

       pz[j]->Reset();
	    for(int ix=0;ix<pxsub[j]->GetXaxis()->GetNbins();ix++)
	      {
		int ival1=pxsub[j]->GetXaxis()->GetBinCenter(ix+1);
		int iz1=ival1/100;
		int in1=ival1-iz1*100;
		int ia1=iz1+in1;
		int ias=0;
		if(pxsub[j]->GetBinContent(ix+1)>0)
		  {
		//cout<<asy<<" "<<ias<<endl;
		tot[j][iz1][ias]=tot[j][iz1][ias]+pxsub[j]->GetBinContent(ix+1);
		tia[j][iz1][ias][ia1]=tia[j][iz1][ias][ia1]+pxsub[j]->GetBinContent(ix+1);
		pz[j]->Fill((float)iz1,pxsub[j]->GetBinContent(ix+1));
		      }
	      }
px[j]->Sumw2();
pz[j]->Sumw2();
	    pz[j]->Scale(1/pz[j]->Integral());
	    px[j]->Scale(1/px[j]->Integral());
	    cout<<j<<endl;

  //      cout<<j<<endl;
  	  gt[j]=new TGraphErrors();
  	  gt[j]->SetMarkerStyle(sym[j]);
  	  gt[j]->SetMarkerColor(col[j]);
// 	  gg[j]=(TGraphErrors**)malloc(2*sizeof(TGraphErrors*));
	  //	   	  ggr[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggs[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggy[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggt[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggt2[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));



//  	  gtsub[j]=new TGraphErrors();
//  	  gtsub[j]->SetMarkerStyle(sym[j]);
//  	  gtsub[j]->SetMarkerColor(col[j]);
// 	  ggsub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggrsub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggssub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggysub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggtsub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));
// 	  ggt2sub[j]=(TGraphErrors**)malloc(27*sizeof(TGraphErrors*));

// 	  for(int ib=0;ib<2;ib++)
// 	    {
// 	      gg[j][iz]=new TGraphErrors();
// 	      gg[j][iz]->SetMarkerStyle(sym[j]);
//   	  gg[j][iz]->SetMarkerColor(col[j]);
//  	  gg[j][iz]->SetLineColor(col[j]);
//  	  gg[j][iz]->SetMarkerSize(1.2);

// 	    }

//        for(int iz=0;iz<27;iz++)
//  	{
//  	  
//  	  


//  	  ggr[j][iz]=new TGraphErrors();
//  	  ggr[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggr[j][iz]->SetMarkerColor(col[j]);
//  	  ggr[j][iz]->SetLineColor(col[j]);
// 	  ggr[j][iz]->SetMarkerSize(1.2);

//  	  ggs[j][iz]=new TGraphErrors();
//  	  ggs[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggs[j][iz]->SetMarkerColor(col[j]);
// 	  ggs[j][iz]->SetLineColor(col[j]);
// 	  ggs[j][iz]->SetMarkerSize(1.2);

//  	  ggy[j][iz]=new TGraphErrors();
//  	  ggy[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggy[j][iz]->SetMarkerColor(col[j]);
// 	  ggy[j][iz]->SetLineColor(col[j]);
// 	  ggy[j][iz]->SetMarkerSize(1.2);

//  	  ggt[j][iz]=new TGraphErrors();
//  	  ggt[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggt[j][iz]->SetMarkerColor(col[j]);
// 	  ggt[j][iz]->SetLineColor(col[j]);
// 	  ggt[j][iz]->SetMarkerSize(1.2);

//  	  ggt2[j][iz]=new TGraphErrors();
//  	  ggt2[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggt2[j][iz]->SetMarkerColor(col[j]);
// 	  ggt2[j][iz]->SetLineColor(col[j]);
// 	  ggt2[j][iz]->SetMarkerSize(1.2);


//  	  ggsub[j][iz]=new TGraphErrors();
//  	  ggsub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggsub[j][iz]->SetMarkerColor(col[j]);
// 	  ggsub[j][iz]->SetLineColor(col[j]);
// 	  ggsub[j][iz]->SetMarkerSize(1.2);

//  	  ggrsub[j][iz]=new TGraphErrors();
//  	  ggrsub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggrsub[j][iz]->SetMarkerColor(col[j]);
//  	  ggrsub[j][iz]->SetLineColor(col[j]);
// 	  ggrsub[j][iz]->SetMarkerSize(1.2);

//  	  ggssub[j][iz]=new TGraphErrors();
//  	  ggssub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggssub[j][iz]->SetMarkerColor(col[j]);
// 	  ggssub[j][iz]->SetLineColor(col[j]);
// 	  ggssub[j][iz]->SetMarkerSize(1.2);

//  	  ggysub[j][iz]=new TGraphErrors();
//  	  ggysub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggysub[j][iz]->SetMarkerColor(col[j]);
// 	  ggysub[j][iz]->SetLineColor(col[j]);
// 	  ggysub[j][iz]->SetMarkerSize(1.2);

//  	  ggtsub[j][iz]=new TGraphErrors();
//  	  ggtsub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggtsub[j][iz]->SetMarkerColor(col[j]);
// 	  ggtsub[j][iz]->SetLineColor(col[j]);
// 	  ggtsub[j][iz]->SetMarkerSize(1.2);

//  	  ggt2sub[j][iz]=new TGraphErrors();
//  	  ggt2sub[j][iz]->SetMarkerStyle(sym[j]);
//  	  ggt2sub[j][iz]->SetMarkerColor(col[j]);
// 	  ggt2sub[j][iz]->SetLineColor(col[j]);
// 	  ggt2sub[j][iz]->SetMarkerSize(1.2);


//  	}

//        tz(px[j],gg[j],gt[j]);
//        tz(pxsub[j],ggsub[j],gtsub[j]);
     }

 

     for(int j=0;j<nf;j++)
      {
	pn[j]=0;
	pnz[j]=0;
	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
	  {
	    for(int k=0;k<nf;k++)
	      {
		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
		  {
		    pn[j]=(TH1D*)px[k]->Clone(Form("pn_%d",j));
		    pn[j]->Divide(px[j]);

		    pnz[j]=(TH1D*)pz[k]->Clone(Form("pnz_%d",j));
		    pnz[j]->Divide(pz[j]);

		    break;
		  }
	      }
	  }
      }
	    


//  int ** ainf=(int**)malloc(nf*sizeof (int*));
//   int ** asup=(int**)malloc(nf*sizeof (int*));
//  int ** ainfsub=(int**)malloc(nf*sizeof (int*));
//   int ** asupsub=(int**)malloc(nf*sizeof (int*));
//   for(int i=0;i<nf;i++)
//     {
//       ainf[i]=(int*)malloc(27*sizeof(int));
//       asup[i]=(int*)malloc(27*sizeof(int));
//       ainfsub[i]=(int*)malloc(27*sizeof(int));
//       asupsub[i]=(int*)malloc(27*sizeof(int));

//     }
//   for(int j=0;j<nf;j++)
//     {
//       for(int iz=0;iz<27;iz++)
// 	{
// 	  ainf[j][iz]=-1;
// 	  asup[j][iz]=-1;
// 	  ainfsub[j][iz]=-1;
// 	  asupsub[j][iz]=-1;

// 	}
//     }
  
//   for(int j=0;j<nf;j++)
//     {
     
//       if(gg[j][5]->GetN()>0)
// 	{
	  
// 	  int i1=1;
// 	  while(i1==1)
// 	    {
// 	      i1=0;
// 	      double *x=gg[j][5]->GetX();
// 	  for(int k=0;k<gg[j][5]->GetN();k++)
// 	    {
// 	      if(x[k]<10)
// 		{
		 
// 		gg[j][5]->RemovePoint(k);
// 		i1=1;
// 		break;
// 		}
// 	    }
// 	    }
// 	}
//     }	  


//  for(int j=0;j<nf;j++)
//     {
     
//       if(ggsub[j][5]->GetN()>0)
// 	{
	  
// 	  int i1=1;
// 	  while(i1==1)
// 	    {
// 	      i1=0;
// 	      double *x=ggsub[j][5]->GetX();
// 	  for(int k=0;k<ggsub[j][5]->GetN();k++)
// 	    {
// 	      if(x[k]<10)
// 		{
		 
// 		ggsub[j][5]->RemovePoint(k);
// 		i1=1;
// 		break;
// 		}
// 	    }
// 	    }
// 	}
//     }	 


//     for(int j=0;j<nf;j++)
//       {
//         if(!stit[j].Contains("stiff")&&!stit[j].Contains("soft"))
//  	 {
// 	   for(int iz=0;iz<27;iz++)
// 	     {
// 	       if(gg[j][iz]->GetN()>0)
// 		 {
	
// 	   double *x=gg[j][iz]->GetX();
//  	   ainf[j][iz]=1000;
//  	   asup[j][iz]=-1000;
// 	   for(int i=0;i<gg[j][iz]->GetN();i++)
// 	     {
// 	       if(x[i]<ainf[j][iz])ainf[j][iz]=x[i];
// 	       if(x[i]>asup[j][iz])asup[j][iz]=x[i];
// 	     }
	   
//  	 }
// 	 }
//       }
//       }

//     for(int j=0;j<nf;j++)
//       {
// 	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
// 	  {
// 	    for(int k=0;k<nf;k++)
// 	      {
// 		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
// 		  {
// 		for(int iz=0;iz<27;iz++)
// 	     {
// 	       ainf[j][iz]= ainf[k][iz];
// 	       asup[j][iz]= asup[k][iz];
	       
	       
// 	     }    

// 		    break;
// 		  }
// 	      }
// 	  }
//       }
	    


//    for(int j=0;j<nf;j++)
//       {
//         if(!stit[j].Contains("stiff")&&!stit[j].Contains("soft"))
//  	 {
// 	   for(int iz=0;iz<27;iz++)
// 	     {
// 	       if(ggsub[j][iz]->GetN()>0)
// 		 {
	
// 	   double *x=ggsub[j][iz]->GetX();
//  	   ainfsub[j][iz]=1000;
//  	   asupsub[j][iz]=-1000;
// 	   for(int i=0;i<ggsub[j][iz]->GetN();i++)
// 	     {
// 	       if(x[i]<ainfsub[j][iz])ainfsub[j][iz]=x[i];
// 	       if(x[i]>asupsub[j][iz])asupsub[j][iz]=x[i];
// 	     }
	   
//  	 }
// 	 }
//       }
//       }

//     for(int j=0;j<nf;j++)
//       {
// 	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
// 	  {
// 	    for(int k=0;k<nf;k++)
// 	      {
// 		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
// 		  {
// 		for(int iz=0;iz<27;iz++)
// 	     {
// 	       ainf[j][iz]= ainf[k][iz];
// 	       asup[j][iz]= asup[k][iz];
	       
	       
// 	     }    

// 		    break;
// 		  }
// 	      }
// 	  }
//       }

//     for(int j=0;j<nf;j++)
//       {
// 	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
// 	  {
// 	    for(int k=0;k<nf;k++)
// 	      {
// 		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
// 		  {
// 		for(int iz=0;iz<27;iz++)
// 	     {
// 	       ainfsub[j][iz]= ainfsub[k][iz];
// 	       asupsub[j][iz]= asupsub[k][iz];
	       
	       
// 	     }    

// 		    break;
// 		  }
// 	      }
// 	  }
//       }


//   for(int j=0;j<nf;j++)
//      {
//  	  for(int iz=1;iz<9;iz++)
//  	    {
//  	      double *x,*y,*ey;
//  	      if(gg[j][iz]->GetN()>0)
//  		{
//  		  x=gg[j][iz]->GetX();
//  		  y=gg[j][iz]->GetY();
//  		  ey=gg[j][iz]->GetEY();
//  		  int k0=-1;
//  		  for(int k=0;k<gg[j][iz]->GetN();k++)
//  		    {
//  		      if(x[k]==rif[iz])
//  			{
//  			  k0=k;
//  			  break;
//  			}
//  		    }
//  		  if(k0>=0) 		    
// {
//  		  for(int k=0;k<gg[j][iz]->GetN();k++)
//  		    {

// 		      //		      if(iz==1 && x[k]>3)continue;
// 		      //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
// 		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;

		      
// 		      float rap=y[k]/y[k0];
// 		      float erap=rap*(ey[k]/y[k]+ey[k0]/y[k0]);
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggr[j][iz]->SetPointError(ggr[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


//  		}

//  	    }


     


//      }


//   for(int j=0;j<nf;j++)
//      {
//  	  for(int iz=1;iz<9;iz++)
//  	    {
//  	      double *x,*y,*ey;
//  	      if(ggsub[j][iz]->GetN()>0)
//  		{
//  		  x=ggsub[j][iz]->GetX();
//  		  y=ggsub[j][iz]->GetY();
//  		  ey=ggsub[j][iz]->GetEY();
//  		  int k0=-1;
//  		  for(int k=0;k<ggsub[j][iz]->GetN();k++)
//  		    {
//  		      if(x[k]==rif[iz])
//  			{
//  			  k0=k;
//  			  break;
//  			}
//  		    }
//  		  if(k0>=0) 		    
// {
//  		  for(int k=0;k<ggsub[j][iz]->GetN();k++)
//  		    {

// 		      //		      if(iz==1 && x[k]>3)continue;
// 		      //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
// 		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;

		      
// 		      float rap=y[k]/y[k0];
// 		      float erap=rap*(ey[k]/y[k]+ey[k0]/y[k0]);
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggrsub[j][iz]->SetPoint(ggrsub[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggrsub[j][iz]->SetPointError(ggrsub[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


//  		}

//  	    }


     


//      }




//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(ggr[j][iz]->GetN()>0)
//  		{
//  		  x=ggr[j][iz]->GetX();
//  		  y=ggr[j][iz]->GetY();
//  		  ey=ggr[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(ggr[iq][iz]->GetN()>0)
//  		{
//  		  xq=ggr[iq][iz]->GetX();
//  		  yq=ggr[iq][iz]->GetY();
//  		  eyq=ggr[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<ggr[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<ggr[iq][iz]->GetN();l++)
// 			{

//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
// 		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;
  
// 		      float rap=yq[k0]/y[k];
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggs[j][iz]->SetPoint(ggs[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggs[j][iz]->SetPointError(ggs[j][iz]->GetN()-1,0.0001,erap);
		      
		      
//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 



//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(ggrsub[j][iz]->GetN()>0)
//  		{
//  		  x=ggrsub[j][iz]->GetX();
//  		  y=ggrsub[j][iz]->GetY();
//  		  ey=ggrsub[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(ggrsub[iq][iz]->GetN()>0)
//  		{
//  		  xq=ggrsub[iq][iz]->GetX();
//  		  yq=ggrsub[iq][iz]->GetY();
//  		  eyq=ggrsub[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<ggrsub[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<ggrsub[iq][iz]->GetN();l++)
// 			{

//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
// 		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;

// 		      float rap=yq[k0]/y[k];
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggssub[j][iz]->SetPoint(ggssub[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggssub[j][iz]->SetPointError(ggssub[j][iz]->GetN()-1,0.0001,erap);
		      
		      
//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 





//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(gg[j][iz]->GetN()>0)
//  		{
//  		  x=gg[j][iz]->GetX();
//  		  y=gg[j][iz]->GetY();
//  		  ey=gg[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(gg[iq][iz]->GetN()>0)
//  		{
//  		  xq=gg[iq][iz]->GetX();
//  		  yq=gg[iq][iz]->GetY();
//  		  eyq=gg[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<gg[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<gg[iq][iz]->GetN();l++)
// 			{
//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
//  		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;

// 		      float rap=yq[k0]*hnev[j]->Integral(102,102)/(y[k]*hnev[iq]->Integral(102,102));
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]+1/sqrt(hnev[j]->Integral(102,102))+1/sqrt(hnev[iq]->Integral(102,102)));
		

//  		      ggy[j][iz]->SetPoint(ggy[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggy[j][iz]->SetPointError(ggy[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 



//   for(int j=0;j<nf;j++)
//      {
//  	  for(int iz=1;iz<9;iz++)
//  	    {
//  	      double *x,*y,*ey;
//  	      if(gg[j][iz]->GetN()>0)
//  		{
//  		  x=gg[j][iz]->GetX();
//  		  y=gg[j][iz]->GetY();
//  		  ey=gg[j][iz]->GetEY();
// 		  float tot=0;
//  		  for(int k=0;k<gg[j][iz]->GetN();k++)
//  		    {

// 		      //		      if(iz==1 && x[k]>3)continue;
// 		      //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
// 		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;
// 		      tot=tot+y[k];
// 		    }
// 		  for(int k=0;k<gg[j][iz]->GetN();k++)
// 		    {
// 		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;
 
// 		      float rap=y[k]/tot;
// 		      float erap=rap*(ey[k]/y[k]+1/sqrt(tot));
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggt[j][iz]->SetPoint(ggt[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggt[j][iz]->SetPointError(ggt[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
 		    


//  		}

//  	    }


     


//      }



//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(ggsub[j][iz]->GetN()>0)
//  		{
//  		  x=ggsub[j][iz]->GetX();
//  		  y=ggsub[j][iz]->GetY();
//  		  ey=ggsub[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(ggsub[iq][iz]->GetN()>0)
//  		{
//  		  xq=ggsub[iq][iz]->GetX();
//  		  yq=ggsub[iq][iz]->GetY();
//  		  eyq=ggsub[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<ggsub[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<ggsub[iq][iz]->GetN();l++)
// 			{
//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
//  		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;

// 		      float rap=yq[k0]*hnev[j]->Integral(102,102)/(y[k]*hnev[iq]->Integral(102,102));
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]+1/sqrt(hnev[j]->Integral(102,102))+1/sqrt(hnev[iq]->Integral(102,102)));
		

//  		      ggysub[j][iz]->SetPoint(ggysub[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggysub[j][iz]->SetPointError(ggysub[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 



//   for(int j=0;j<nf;j++)
//      {
//  	  for(int iz=1;iz<9;iz++)
//  	    {
//  	      double *x,*y,*ey;
//  	      if(ggsub[j][iz]->GetN()>0)
//  		{
//  		  x=ggsub[j][iz]->GetX();
//  		  y=ggsub[j][iz]->GetY();
//  		  ey=ggsub[j][iz]->GetEY();
// 		  float tot=0;
//  		  for(int k=0;k<ggsub[j][iz]->GetN();k++)
//  		    {

// 		      //		      if(iz==1 && x[k]>3)continue;
// 		      //if(iz==2 && (x[k]>6|| x[k]<3))continue;

		      
// 		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;
// 		      tot=tot+y[k];
// 		    }
// 		  for(int k=0;k<ggsub[j][iz]->GetN();k++)
// 		    {
// 		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;

// 		      float rap=y[k]/tot;
// 		      float erap=rap*(ey[k]/y[k]+1/sqrt(tot));
		
// 		      // 		      ggr[j][iz]->SetPoint(ggr[j][iz]->GetN(),x[k],rap);
//  		      ggtsub[j][iz]->SetPoint(ggtsub[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggtsub[j][iz]->SetPointError(ggtsub[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
 		    


//  		}

//  	    }


     


//      }




//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(ggt[j][iz]->GetN()>0)
//  		{
//  		  x=ggt[j][iz]->GetX();
//  		  y=ggt[j][iz]->GetY();
//  		  ey=ggt[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(ggt[iq][iz]->GetN()>0)
//  		{
//  		  xq=ggt[iq][iz]->GetX();
//  		  yq=ggt[iq][iz]->GetY();
//  		  eyq=ggt[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<ggt[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<ggt[iq][iz]->GetN();l++)
// 			{
//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
//  		      if(x[k]<ainf[j][iz])continue;
// 		      if(x[k]>asup[j][iz])continue;

// 		      float rap=yq[k0]/y[k];
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		

//  		      ggt2[j][iz]->SetPoint(ggt2[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggt2[j][iz]->SetPointError(ggt2[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 





//    for(int j=0;j<nf;j++)
//      {
//        if(stit[j].Contains("stiff") || stit[j].Contains("soft"))
// 	 {
//        for(int iz=1;iz<9;iz++)
// 	 {
// 	   double *x,*y,*ey; //sim
//  	      if(ggtsub[j][iz]->GetN()>0)
//  		{
//  		  x=ggtsub[j][iz]->GetX();
//  		  y=ggtsub[j][iz]->GetY();
//  		  ey=ggtsub[j][iz]->GetEY();
		  

//        for(int iq=0;iq<nf;iq++)
// 	 {
// 	   if(iq!=j &&!stit[iq].Contains("stiff") &&!stit[iq].Contains("soft"))
// 	     {
// 	   if(tipo[j]==tipo[iq])
// 	     {

	 
// 	       double *xq,*yq,*eyq; //exp
//  	      if(ggtsub[iq][iz]->GetN()>0)
//  		{
//  		  xq=ggtsub[iq][iz]->GetX();
//  		  yq=ggtsub[iq][iz]->GetY();
//  		  eyq=ggtsub[iq][iz]->GetEY();
 		  
//  		  for(int k=0;k<ggtsub[j][iz]->GetN();k++)
//  		    {
// 		      int k0=-1;
// 		      for(int l=0;l<ggt[iq][iz]->GetN();l++)
// 			{
//  		      if(x[k]==xq[l])
//  			{
//  			  k0=l;
//  			  break;
//  			}
// 			}
//  		  if(k0>=0) 		    
// 		    {
//  		      if(x[k]<ainfsub[j][iz])continue;
// 		      if(x[k]>asupsub[j][iz])continue;

// 		      float rap=yq[k0]/y[k];
// 		      float erap=rap*(ey[k]/y[k]+eyq[k0]/yq[k0]);
		

//  		      ggt2sub[j][iz]->SetPoint(ggt2sub[j][iz]->GetN(),x[k]+tipo[j]*0.1,rap);
//  		      ggt2sub[j][iz]->SetPointError(ggt2sub[j][iz]->GetN()-1,0.0001,erap);
		      

//  		    }
//  		    }


 	

//  	    }


// 	     }
// 	     }
// 	 }
// 		}
// 	 }
// 	 }
//      } 












//    TLine *l=new TLine(0,1,30,1);
//    l->SetLineStyle(2);
//    l->SetLineWidth(2);

//    l->SetLineColor(kBlue-7);
		 
//    TCanvas *c[27];
//    TH2F *hf[27];
//    float amin[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amax[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float ymin[27]={0,1e-1,1e-5,1e-3,1e-5,1e-5,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *leg[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        c[iz]=new TCanvas(Form("c%d",iz),Form("c%d",iz));
//        c[iz]->Draw();
//        hf[iz]=new TH2F(Form("hf%d",iz),Form("Z=%d",iz),100,amin[iz],amax[iz],100,ymin[iz],10);
//        hf[iz]->GetXaxis()->SetTitle("A");
//        hf[iz]->GetYaxis()->SetTitle(Form("Y(A)/Y(%d)",rif[iz]));
//        hf[iz]->SetStats(kFALSE);
//        hf[iz]->Draw();
// leg[iz]=new TLegend(0.79,0.65,0.97,0.98);
//        for(int k=0;k<nf;k++)
// 	 {
//        ggr[k][iz]->Draw("pl");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        leg[iz]->AddEntry(ggr[k][iz],tit[k],"p");

// 	 }
//        leg[iz]->Draw();
//      }




//   TCanvas *cb[27];
//    TH2F *hfb[27];
//    float aminb[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxb[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//    float yminb[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legb[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        cb[iz]=new TCanvas(Form("cb%d",iz),Form("cb%d",iz));
//        cb[iz]->Draw();
//        hfb[iz]=new TH2F(Form("hfb%d",iz),Form("Z=%d",iz),100,aminb[iz],amaxb[iz],100,yminb[iz],100);
//        hfb[iz]->GetXaxis()->SetTitle("A");
//        hfb[iz]->GetYaxis()->SetTitle("(Y/Yrif)EXP/(Y/Yrif)Sim");
//        hfb[iz]->SetStats(kFALSE);
//        hfb[iz]->Draw();
//        l->Draw();
// legb[iz]=new TLegend(0.12,0.14,0.3,0.47);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(ggs[k][iz]->GetN()>0)
// 	     {
//        ggs[k][iz]->Draw("p");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legb[iz]->AddEntry(ggs[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legb[iz]->Draw();
//      }


//   TCanvas *cy[27];
//    TH2F *hfy[27];
//    float aminy[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxy[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//    float yminy[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legy[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        cy[iz]=new TCanvas(Form("cy%d",iz),Form("cy%d",iz));
//        cy[iz]->Draw();
//        hfy[iz]=new TH2F(Form("hfy%d",iz),Form("Z=%d",iz),100,aminy[iz],amaxy[iz],100,yminy[iz],100);
//        hfy[iz]->GetXaxis()->SetTitle("A");
//        hfy[iz]->GetYaxis()->SetTitle("YEXP/YSIM (norm neventi)");
//        hfy[iz]->SetStats(kFALSE);
//        hfy[iz]->Draw();
//        l->Draw();
// legy[iz]=new TLegend(0.12,0.14,0.3,0.47);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(ggy[k][iz]->GetN()>0)
// 	     {
//        ggy[k][iz]->Draw("p");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legy[iz]->AddEntry(ggy[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legy[iz]->Draw();
//      }





//   TCanvas *ct[27];
//    TH2F *hft[27];
//    float amint[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxt[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float ymint[27]={0,1e-1,1e-5,1e-3,1e-5,1e-6,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legt[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        ct[iz]=new TCanvas(Form("ct%d",iz),Form("ct%d",iz));
//        ct[iz]->Draw();
//        hft[iz]=new TH2F(Form("hft%d",iz),Form("Z=%d",iz),100,amint[iz],amaxt[iz],100,ymint[iz],10);
//        hft[iz]->GetXaxis()->SetTitle("A");
//        hft[iz]->GetYaxis()->SetTitle("Y(A)(norm 1)");
//        hft[iz]->SetStats(kFALSE);
//        hft[iz]->Draw();
// legt[iz]=new TLegend(0.79,0.65,0.97,0.98);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(!stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
// 	     {
//        ggt[k][iz]->Draw("pl");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legt[iz]->AddEntry(ggt[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legt[iz]->Draw();
//      }

//   TCanvas *ct2[27];
//    TH2F *hft2[27];
//    float amint2[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxt2[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//    float ymint2[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legt2[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        ct2[iz]=new TCanvas(Form("ct2_%d",iz),Form("ct2_%d",iz));
//        ct2[iz]->Draw();
//        hft2[iz]=new TH2F(Form("hft2_%d",iz),Form("Z=%d",iz),100,amint2[iz],amaxt2[iz],100,ymint2[iz],100);
//        hft2[iz]->GetXaxis()->SetTitle("A");
//        hft2[iz]->GetYaxis()->SetTitle("YEXP/YSIM (norm1)");
//        hft2[iz]->SetStats(kFALSE);
//        hft2[iz]->Draw();
//        l->Draw();
// legt2[iz]=new TLegend(0.12,0.14,0.3,0.47);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(ggt2[k][iz]->GetN()>0)
// 	     {
//        ggt2[k][iz]->Draw("p");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legt2[iz]->AddEntry(ggy[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legt2[iz]->Draw();
//      }





//   TCanvas *csub[27];
//    TH2F *hfsub[27];
//    float aminsub[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxsub[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float yminsub[27]={0,1e-1,1e-5,1e-3,1e-5,1e-5,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legsub[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        csub[iz]=new TCanvas(Form("csub%d",iz),Form("csub%d",iz));
//        csub[iz]->Draw();
//        hfsub[iz]=new TH2F(Form("hfsub%d",iz),Form("Z=%d",iz),100,aminsub[iz],amaxsub[iz],100,yminsub[iz],10);
//        hfsub[iz]->GetXaxis()->SetTitle("A");
//        hfsub[iz]->GetYaxis()->SetTitle(Form("Y(A)/Y(%d) SUB",rif[iz]));
//        hfsub[iz]->SetStats(kFALSE);
//        hfsub[iz]->Draw();
// legsub[iz]=new TLegend(0.79,0.65,0.97,0.98);
//        for(int k=0;k<nf;k++)
// 	 {
//        ggrsub[k][iz]->Draw("pl");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legsub[iz]->AddEntry(ggrsub[k][iz],tit[k],"p");

// 	 }
//        legsub[iz]->Draw();
//      }


//  TCanvas *ctsub[27];
//    TH2F *hftsub[27];
//    float amintsub[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxtsub[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float ymintsub[27]={0,1e-1,1e-5,1e-3,1e-5,1e-6,1e-4,1e-5,1e-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legtsub[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        ctsub[iz]=new TCanvas(Form("ctsub%d",iz),Form("ctsub%d",iz));
//        ctsub[iz]->Draw();
//        hftsub[iz]=new TH2F(Form("hftsub%d",iz),Form("Z=%d",iz),100,amintsub[iz],amaxtsub[iz],100,ymintsub[iz],10);
//        hftsub[iz]->GetXaxis()->SetTitle("A");
//        hftsub[iz]->GetYaxis()->SetTitle("Y(A)(norm 1)");
//        hftsub[iz]->SetStats(kFALSE);
//        hftsub[iz]->Draw();
// legtsub[iz]=new TLegend(0.79,0.65,0.97,0.98);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(!stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
// 	     {
//        ggtsub[k][iz]->Draw("pl");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legtsub[iz]->AddEntry(ggtsub[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legtsub[iz]->Draw();
//      }

//   TCanvas *ct2sub[27];
//    TH2F *hft2sub[27];
//    float amint2sub[27]={0,0.5,2.5,5.5,6,6,8,11,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    float amaxt2sub[27]={0,3.5,8.5,10,14,15,19,20,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

//    float ymint2sub[27]={0,1e-1,1e-3,1e-3,1e-5,1e-1,1e-4,1e-5,1e-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//    TLegend *legt2sub[27];
//    for(int iz=1;iz<9;iz++)
//      {

//        ct2sub[iz]=new TCanvas(Form("ct2sub_%d",iz),Form("ct2sub_%d",iz));
//        ct2sub[iz]->Draw();
//        hft2sub[iz]=new TH2F(Form("hft2sub_%d",iz),Form("Z=%d",iz),100,amint2sub[iz],amaxt2sub[iz],100,ymint2sub[iz],100);
//        hft2sub[iz]->GetXaxis()->SetTitle("A");
//        hft2sub[iz]->GetYaxis()->SetTitle("YEXP/YSIM (norm1)");
//        hft2sub[iz]->SetStats(kFALSE);
//        hft2sub[iz]->Draw();
//        l->Draw();
// legt2sub[iz]=new TLegend(0.12,0.14,0.3,0.47);
//        for(int k=0;k<nf;k++)
// 	 {
// 	   if(ggt2sub[k][iz]->GetN()>0)
// 	     {
//        ggt2sub[k][iz]->Draw("p");
//        gPad->SetLogy(1);
//        gPad->SetGridx(kFALSE);
//        gPad->SetGridy(kFALSE);
       
//        legt2sub[iz]->AddEntry(ggysub[k][iz],tit[k],"p");
// 	     }
// 	 }
//        legt2sub[iz]->Draw();
//      }

      float m2=0;
      int jm2;
      char stile[2][10]={"pE","hist"};
  for(int j=0;j<nf;j++)
    {
      if(px[j]->GetMaximum()>m2)
	{
	  m2=px[j]->GetMaximum();
	  jm2=j;
  }

  px[j]->GetXaxis()->SetRangeUser(8,29);
      px[j]->SetMarkerStyle(sym[j]);
 px[j]->SetMarkerColor(col[j]);
px[j]->SetLineColor(col[j]);
    if(stit[j].Contains("soft"))px[j]->SetLineStyle(2);

      px[j]->GetXaxis()->SetNdivisions(505);      
    }
  for(int j=0;j<nf;j++)
    {


px[j]->GetYaxis()->SetNdivisions(505);
pz[j]->GetYaxis()->SetNdivisions(505);

      px[j]->GetXaxis()->SetLabelSize(0.08);
      px[j]->GetYaxis()->SetLabelSize(0.08);
      px[j]->GetYaxis()->SetLabelFont(42);
      px[j]->GetXaxis()->SetTitleSize(0.08);
      px[j]->GetYaxis()->SetTitleSize(0.08);


      pz[j]->GetXaxis()->SetLabelSize(0.08);
      pz[j]->GetYaxis()->SetLabelSize(0.08);
      pz[j]->GetYaxis()->SetLabelFont(42);
      pz[j]->GetXaxis()->SetTitleSize(0.08);
      pz[j]->GetYaxis()->SetTitleSize(0.08);



     px[j]->SetTitle("");
      px[j]->SetStats(kFALSE);
 pz[j]->SetStats(kFALSE);

      px[j]->GetYaxis()->SetTitle("P(Z)");
 pz[j]->GetYaxis()->SetTitle("P(Z)");
      

      px[j]->GetYaxis()->SetTitleFont(42);

     px[j]->GetYaxis()->SetTitleOffset(1.35);
     px[j]->GetXaxis()->SetTitle("Z");
     pz[j]->GetXaxis()->SetTitle("Z");
pz[j]->GetYaxis()->SetTitleOffset(2);
    

     px[j]->GetYaxis()->SetLabelOffset(0.008);
     pz[j]->GetYaxis()->SetLabelOffset(0.008);

    }

  TCanvas *c2=new TCanvas("c2","c2",0,0,1200,600);
  c2->Draw();
  c2->Divide(2,1);
  TLatex *l2[2];
//   char let2[2][4]={"c)","d)"};
//   TLegend *leg1=new TLegend(0.4,0.6,0.8,0.8);
//   leg1->AddEntry(h1x[0],"Exp","p");
//   leg1->AddEntry(h1x[1],"Sim","l");
//   leg1->SetTextSize(0.08);
//   leg1->SetBorderSize(0);
//   leg1->SetFillStyle(0);

//   TLegend *leg1a=new TLegend(0.25,0.6,0.45,0.8);
//   leg1a->AddEntry(h1y[0],"Exp","p");
//   leg1a->AddEntry(h1y[1],"Sim","l");
//   leg1a->SetTextSize(0.08);
//   leg1a->SetBorderSize(0);
//   leg1a->SetFillStyle(0);


      c2->cd(1);
    int is=1;
 if(!stit[jm2].Contains("stiff") && !stit[jm2].Contains("soft"))
   {
     is=0;
   }
  px[jm2]->Draw(Form("%s",stile[is]));
  for(int j=0;j<nf;j++)
    {


   if(j==jm2)continue;

      is=1;
 if(!stit[j].Contains("stiff") && !stit[j].Contains("soft"))
   {
     is=0;
   }
      px[j]->Draw(Form("%ssame",stile[is]));  
    


      //      px[0]->SetMarkerStyle(20);
      // px[0]->SetMarkerColor(1);

      //px[j]->SetLineColor(1);




     

   
	  
    }

     
      gPad->SetLogy(kTRUE);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      gPad->SetBottomMargin(0.17);
      gPad->SetLeftMargin(0.21);
      gPad->SetRightMargin(0.01);
 gPad->Modified();
      gPad->Update();
     
 //      l2[j]=new TLatex(0.8,0.8,let2[j]);
//       l2[j]->SetNDC();
//       l2[j]->SetTextSize(0.08);
//       l2[j]->Draw();
    
//   c2->cd(2);
//   leg1->Draw();
// c2->cd(1);
//  leg1a->Draw();


      c2->cd(2);
      float m1=0;
      int jm1;
  
  for(int j=0;j<nf;j++)
    {
      if(pz[j]->GetMaximum()>m1)
	{
	  m1=pz[j]->GetMaximum();
	  jm1=j;
  }

     
      pz[j]->GetXaxis()->SetRangeUser(0,19);
      pz[j]->SetMarkerStyle(sym[j]);
 pz[j]->SetMarkerColor(col[j]);
pz[j]->SetLineColor(col[j]);


    pz[j]->GetXaxis()->SetNdivisions(505);     
    pz[j]->GetYaxis()->SetNdivisions(503);     
    if(stit[j].Contains("soft"))pz[j]->SetLineStyle(2);
    }
  cout<<m1<<" "<<jm1<<endl;
    is=1;
 if(!stit[jm1].Contains("stiff") && !stit[jm1].Contains("soft"))
   {
     is=0;
   }
  pz[jm1]->Draw(Form("%s",stile[is]));

  for(int j=0;j<nf;j++)
    {
      if(j==jm1)continue;

      is=1;
 if(!stit[j].Contains("stiff") && !stit[j].Contains("soft"))
   {
     is=0;
   }
      pz[j]->Draw(Form("%ssame",stile[is]));


    }
     
      gPad->SetLogy(kTRUE);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      gPad->SetBottomMargin(0.17);
      gPad->SetLeftMargin(0.21);
      gPad->SetRightMargin(0.01);
 gPad->Modified();
      gPad->Update();
     
 //      l2[j]=new TLatex(0.8,0.8,let2[j]);
//       l2[j]->SetNDC();
//       l2[j]->SetTextSize(0.08);
//       l2[j]->Draw();
    
//   c2->cd(2);
//   leg1->Draw();
// c2->cd(1);
//  leg1a->Draw();


      cout<<"qui"<<endl;


      float m3=0;
      int jm3;

  for(int j=0;j<nf;j++)
    {
      if(pn[j]!=0)
	{
      if(pn[j]->GetMaximum()>m3)
	{
	  m3=pn[j]->GetMaximum();
	  jm3=j;
	}
	

  pn[j]->GetXaxis()->SetRangeUser(8,29);
      pn[j]->SetMarkerStyle(sym[j]);
 pn[j]->SetMarkerColor(col[j]);
pn[j]->SetLineColor(col[j]);
  

      pn[j]->GetXaxis()->SetNdivisions(505);      
	}



    }
  
  for(int j=0;j<nf;j++)
    {
 if(pn[j]!=0)
	{     

pn[j]->GetYaxis()->SetNdivisions(505);
   pn[j]->GetXaxis()->SetLabelSize(0.08);
      pn[j]->GetYaxis()->SetLabelSize(0.08);
      pn[j]->GetYaxis()->SetLabelFont(42);
      pn[j]->GetXaxis()->SetTitleSize(0.08);
      pn[j]->GetYaxis()->SetTitleSize(0.08);
     pn[j]->SetTitle("");
      pn[j]->SetStats(kFALSE);
      pn[j]->GetYaxis()->SetTitle("P(Z)");
      pn[j]->GetYaxis()->SetTitleFont(42);

     pn[j]->GetYaxis()->SetTitleOffset(1.35);
     pn[j]->GetXaxis()->SetTitle("Z");
     pn[j]->GetYaxis()->SetLabelOffset(0.008);
	}
 if(pnz[j]!=0)
   {
pnz[j]->GetYaxis()->SetNdivisions(505);

   


      pnz[j]->GetXaxis()->SetLabelSize(0.08);
      pnz[j]->GetYaxis()->SetLabelSize(0.08);
      pnz[j]->GetYaxis()->SetLabelFont(42);
      pnz[j]->GetXaxis()->SetTitleSize(0.08);
      pnz[j]->GetYaxis()->SetTitleSize(0.08);




 pnz[j]->SetStats(kFALSE);


 pnz[j]->GetYaxis()->SetTitle("P(Z)");
      


     pnz[j]->GetXaxis()->SetTitle("Z");
pnz[j]->GetYaxis()->SetTitleOffset(2);
    


     pnz[j]->GetYaxis()->SetLabelOffset(0.008);
   }
    }
cout<<"qqqx"<<endl;

  TCanvas *cn2=new TCanvas("cn2","cn2",0,0,1200,600);
  cn2->Draw();
  cn2->Divide(2,1);
  TLatex *ln2[2];
//   char let2[2][4]={"c)","d)"};
//   TLegend *leg1=new TLegend(0.4,0.6,0.8,0.8);
//   leg1->AddEntry(h1x[0],"Exp","p");
//   leg1->AddEntry(h1x[1],"Sim","l");
//   leg1->SetTextSize(0.08);
//   leg1->SetBorderSize(0);
//   leg1->SetFillStyle(0);

//   TLegend *leg1a=new TLegend(0.25,0.6,0.45,0.8);
//   leg1a->AddEntry(h1y[0],"Exp","p");
//   leg1a->AddEntry(h1y[1],"Sim","l");
//   leg1a->SetTextSize(0.08);
//   leg1a->SetBorderSize(0);
//   leg1a->SetFillStyle(0);


      cn2->cd(1);

  pn[jm3]->Draw("pE");
  for(int j=0;j<nf;j++)
    {


   if(j==jm3)continue;


 if(pn[j]!=0)
   {
      pn[j]->Draw("pEsame");  
   }


      //      pn[0]->SetMarkerStyle(20);
      // pn[0]->SetMarkerColor(1);

      //pn[j]->SetLineColor(1);




     

   
	  
    }

     
      gPad->SetLogy(kTRUE);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      gPad->SetBottomMargin(0.17);
      gPad->SetLeftMargin(0.21);
      gPad->SetRightMargin(0.01);
 gPad->Modified();
      gPad->Update();
     
 //      l2[j]=new TLatex(0.8,0.8,let2[j]);
//       l2[j]->SetNDC();
//       l2[j]->SetTextSize(0.08);
//       l2[j]->Draw();
    
//   c2->cd(2);
//   leg1->Draw();
// c2->cd(1);
//  leg1a->Draw();


      cn2->cd(2);
      float m4=0;
      int jm4;
  
  for(int j=0;j<nf;j++)
    {
       if(pnz[j]!=0)
	{
	  if(pnz[j]->GetMaximum()>m1)
	{
	  m4=pnz[j]->GetMaximum();
	  jm4=j;
	}

     
      pnz[j]->GetXaxis()->SetRangeUser(0,19);
      pnz[j]->SetMarkerStyle(sym[j]);
 pnz[j]->SetMarkerColor(col[j]);
pnz[j]->SetLineColor(col[j]);


    pnz[j]->GetXaxis()->SetNdivisions(505);     
    pnz[j]->GetYaxis()->SetNdivisions(503);     
   
	}
    }
  
       pnz[jm4]->Draw("pE");

  for(int j=0;j<nf;j++)
    {
      if(j==jm4)continue;

 if(pnz[j]!=0)
   {
      pnz[j]->Draw("pEsame");
   }

    }
     
      gPad->SetLogy(kTRUE);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      gPad->SetBottomMargin(0.17);
      gPad->SetLeftMargin(0.21);
      gPad->SetRightMargin(0.01);
 gPad->Modified();
      gPad->Update();
     
 //      l2[j]=new TLatex(0.8,0.8,let2[j]);
//       l2[j]->SetNDC();
//       l2[j]->SetTextSize(0.08);
//       l2[j]->Draw();
    
//   c2->cd(2);
//   leg1->Draw();
// c2->cd(1);
//  leg1a->Draw();


      for(int j=0;j<nf;j++)
	{
	  if(pn[j]!=0)
	    {
	      
	      cout<<sum(pn[j])<<" "<<stit[j]<<endl;
	    }
	}



  
}
float sum(TH1D *p)
{
  float out=0;
  int n=0;
  for(int ix=0;ix<p->GetXaxis()->GetNbins();ix++)
    {
      if(p->GetBinContent(ix+1)>0)
	{
      out=out+pow(p->GetBinContent(ix+1)-1,2);
      n++;
	}

    }
      out=sqrt(out);
      out=out/(float)n;
  return out;
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
