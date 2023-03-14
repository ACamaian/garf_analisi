#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TString.h>
#include <TLine.h>
#include <TMath.h>

void tz(TH1D *h,TGraphErrors *g[27],TGraphErrors *gt);
void tz(TH1F *h,TGraphErrors *g[27],TGraphErrors *gt);
// void stiffsoft_tutti(TString s="feb24_23/totFAZIAPRE_4840_amd_geo_stiff.root,feb24_23/totFAZIAPRE_4840_amd_geo_soft.root,feb24_23/totFAZIAPRE_4840.root,feb24_23/totFAZIASYM_484835_amd_geo_stiff.root,feb24_23/totFAZIASYM_484835_amd_geo_soft.root,feb24_23/totFAZIASYM_484835.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_stiff.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_soft.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32.root,feb24_23/totINDRAFAZIA_64Ni64Ni_52_amd_geo_stiff.root,feb24_23/totINDRAFAZIA_64Ni64Ni_52_amd_geo_soft.root,feb24_23/totINDRAFAZIA_64Ni64Ni_52.root")
// void stiffsoft_tutti(TString s="feb24_23/totFAZIAPRE_4840.root,feb24_23/totFAZIASYM_484835.root,feb24_23/totINDRAFAZIA_64Ni64Ni_32.root")
 void stiffsoft_tutti(TString s="feb24_23/totFAZIASYM_484835.root")

//void stiffsoft_tutti(TString s="feb22_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_stiff.root,feb22_23/totINDRAFAZIA_64Ni64Ni_32_amd_geo_soft.root,feb22_23/totINDRAFAZIA_64Ni64Ni_32.root,feb22_23/totINDRAFAZIA_64Ni64Ni_52_amd_geo_stiff.root,feb22_23/totINDRAFAZIA_64Ni64Ni_52_amd_geo_soft.root,feb22_23/totINDRAFAZIA_64Ni64Ni_52.root")

//void stiffsoft_tutti(TString s="feb24_23/totFAZIASYM_484835.root")
//void stiffsoft_tutti(TString s="feb22_23/totFAZIASYM_484835.root,feb22_23/totFAZIAPRE_4840.root")
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
 int * sym2b=(int*)malloc(nf*sizeof (sym));
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
 sym2b[j]=22;
	  sprintf(tit[j],"6464_32ST");
	    }
	  else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
    sym2b[j]=23;
	  sprintf(tit[j],"6464_32SO");
	    }
	  else
	    {
	  sym[j]=20;
sym2b[j]=24;
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
  sym2b[j]=22;
	  sprintf(tit[j],"6464_52ST");
	    }
	  else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
sym2b[j]=23;
	  sprintf(tit[j],"6464_52SO");
	    }
	  else
	    {
	  sym[j]=20;
sym2b[j]=24;
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
 sym2b[j]=22;
	     sprintf(tit[j],"4848_35ST");
	    }
	 else if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
  sym2b[j]=23;
	     sprintf(tit[j],"4848_35SO");
	    }
	 else
	   {
	     sym[j]=21;
 sym2b[j]=25;
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
sym2b[j]=22;
     sprintf(tit[j],"4812_40ST");
	    }
	else  if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
 sym2b[j]=23;
     sprintf(tit[j],"4812_40SO");
	    }
	else
	  {
  sym[j]=34;
 sym2b[j]=28;
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
 sym2b[j]=22;
     sprintf(tit[j],"4812_40ST");
	    }
	else  if(stit[j].Contains("soft"))
	    {
	      sym[j]=32;
sym2b[j]=23;
     sprintf(tit[j],"4812_40SO");
	    }
	else
	  {
  sym[j]=29;
 sym2b[j]=30;
     sprintf(tit[j],"4812_40E");
	  }

	  tipo[j]=4;
    }
    }






  TH2F **h1=(TH2F**)malloc(nf*sizeof(TH2F*));

  TH1F **hnev=(TH1F**)malloc(nf*sizeof(TH1F*));


  //TH2F ****h=(TH2F****)malloc(nf*sizeof(TH2F***));

  int rif[57];
   rif[1]=1;
   rif[2]=4;
   rif[3]=7;
   rif[4]=9;
   rif[5]=10;
   rif[6]=12;
   rif[7]=15;
   rif[8]=16;



//     TH1D **px=(TH1D**)malloc(nf*sizeof(TH1D*));
//     TH1D **py=(TH1D**)malloc(nf*sizeof(TH1D*));

//    TH1D **pxsub=(TH1D**)malloc(nf*sizeof(TH1D*));
//    TGraphErrors **gt=(TGraphErrors**)malloc(nf*sizeof(TGraphErrors*));
//    TGraphErrors ***ggqp=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//    TGraphErrors ***ggr=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//    TGraphErrors ***ggs=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggy=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggt=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggt2=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));


//    TGraphErrors **gtsub=(TGraphErrors**)malloc(nf*sizeof(TGraphErrors*));
    TGraphErrors ****ggrec=(TGraphErrors****)malloc(nf*sizeof(TGraphErrors***));
//    TGraphErrors ***ggrsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//    TGraphErrors ***ggssub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggysub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggtsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
//   TGraphErrors ***ggt2sub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));
    //  TGraphErrors ***gtsub=(TGraphErrors***)malloc(nf*sizeof(TGraphErrors**));

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
       cout<<j<<endl;
            h1[j]=(TH2F*)f[j]->Get("hival1ival2");
	    ggrec[j]=(TGraphErrors***)malloc(10*sizeof(TGraphErrors**));
        for(int ias=0;ias<10;ias++)
 	 {
 	   ggrec[j][ias]=(TGraphErrors**)malloc(57*sizeof(TGraphErrors*));
       for(int iz=0;iz<57;iz++)
 	{
 	  ggrec[j][ias][iz]=new TGraphErrors();
 	  ggrec[j][ias][iz]->SetMarkerStyle(sym[j]);
	  // 	  ggrec[j][ias][iz]->SetMarkerColor(col[j]);
	  //ggrec[j][ias][iz]->SetLineColor(col[j]);
 	  ggrec[j][ias][iz]->SetMarkerColor(ias+1);
	  ggrec[j][ias][iz]->SetLineColor(ias+1);
	  ggrec[j][ias][iz]->SetMarkerSize(1.2);
	  ggrec[j][ias][iz]->SetName(Form("%d_A%d",j,ias));

	}
	 }


	    for(int ix=0;ix<h1[j]->GetXaxis()->GetNbins();ix++)
	      {
		int ival1=h1[j]->GetXaxis()->GetBinCenter(ix+1);
		int iz1=ival1/100;
		int in1=ival1-iz1*100;
		int ia1=iz1+in1;
		for(int iy=0;iy<h1[j]->GetYaxis()->GetNbins();iy++)
		  {
		    if(h1[j]->GetBinContent(ix+1,iy+1)>0)
		      {
		int ival2=h1[j]->GetYaxis()->GetBinCenter(iy+1);
		int iz2=ival2/100;
		int in2=ival2-iz2*100;
		int ia2=iz2+in2;

		//cout<<ix+1<<" "<<ival1<<" "<<iy+1<<" "<<ival2<<" "<<iz1<<" "<<ia1<<" "<<iz2<<" "<<ia2<<endl;
		float asy=(float)(iz1-iz2)/(iz1+iz2);
		int zt=iz1+iz2;
		int at=ia1+ia2;
		int ias=TMath::Nint(asy*10);
		//cout<<asy<<" "<<ias<<endl;
		tot[j][zt][ias]=tot[j][zt][ias]+h1[j]->GetBinContent(ix+1,iy+1);
		tia[j][zt][ias][at]=tia[j][zt][ias][at]+h1[j]->GetBinContent(ix+1,iy+1);
		//		cout<<j<<" "<<zt<<" "<<ias<<" "<<at<<" "<<tia[j][zt][ias][at]<<endl;
		//if(j==0 &&zt==25&&ias==4&&at==54)cout<<tia[0][25][4][54]<<endl;
		if(j==0 && zt==17 && at==36)cout<<" "<<zt<<" "<<at<<" "<<iz1<<" "<<iz2<<" "<<asy<<" "<<ias<<endl;
		      }
		  }

	      }
	    //	    cout<<"F"<<tia[0][25][4][54]<<endl;
	    

	    for(int iz=0;iz<57;iz++)
	      {
		for(int ias=0;ias<10;ias++)
		  {
		    if(tot[j][iz][ias]>0)
		      {
		    for(int ia=0;ia<100;ia++)
		      {
			if(tia[j][iz][ias][ia]>0)
			  {
			float ertia=(tia[j][iz][ias][ia]/tot[j][iz][ias])*(1/sqrt(tia[j][iz][ias][ia])+1/sqrt(tot[j][iz][ias]));
			//cout<<tia[j][iz][ias][ia]<<" "<<tot[j][iz][ias]<<endl;
			tia[j][iz][ias][ia]=tia[j][iz][ias][ia]/tot[j][iz][ias];	

			ggrec[j][ias][iz]->SetPoint(ggrec[j][ias][iz]->GetN(),(float)ia,tia[j][iz][ias][ia]);
			ggrec[j][ias][iz]->SetPointError(ggrec[j][ias][iz]->GetN()-1,0.001,ertia);
			  }
		      }
		      }
		  }
	      }
	    

           

     }



   //cout<<"FB"<<tia[0][25][4][54]<<" "<<ggrec[0][4][25]->GetN()<<endl;

   //return;


 int ** ainfsub=(int**)malloc(nf*sizeof (int*));
  int ** asupsub=(int**)malloc(nf*sizeof (int*));
  for(int i=0;i<nf;i++)
    {

      ainfsub[i]=(int*)malloc(57*sizeof(int));
      asupsub[i]=(int*)malloc(57*sizeof(int));

    }
  for(int j=0;j<nf;j++)
    {
      for(int iz=0;iz<57;iz++)
	{
	
	  ainfsub[j][iz]=-1;
	  asupsub[j][iz]=-1;

	}
    }

   for(int j=0;j<nf;j++)
      {
        if(!stit[j].Contains("stiff")&&!stit[j].Contains("soft"))
 	 {
	   for(int iz=0;iz<57;iz++)
	     {
 	   ainfsub[j][iz]=1000;
 	   asupsub[j][iz]=-1000;
	       for(int ip=0;ip<10;ip++)
		 {

	       if(ggrec[j][ip][iz]->GetN()>0)
		 {
	
	   double *x=ggrec[j][ip][iz]->GetX();
	   double *y=ggrec[j][ip][iz]->GetY();

	   for(int i=0;i<ggrec[j][ip][iz]->GetN();i++)
	     {
	       //cout<<"a="<<j<<" "<<iz<<" "<<x[i]<<" "<<y[i]<<endl;
	       if(x[i]<ainfsub[j][iz])ainfsub[j][iz]=x[i];
	       if(x[i]>asupsub[j][iz])asupsub[j][iz]=x[i];
	       //cout<<asupsub[j][iz]<<endl;
	     }
	   
 	 }
		 }
	 }
      }
      }
   //   return;
  

    for(int j=0;j<nf;j++)
      {
	if(stit[j].Contains("stiff")||stit[j].Contains("soft"))
	  {
	    for(int k=0;k<nf;k++)
	      {
		if(tipo[j]==tipo[k]&&k!=j&& !stit[k].Contains("stiff")&&!stit[k].Contains("soft"))
		  {
		for(int iz=0;iz<57;iz++)
	     {
	       ainfsub[j][iz]= ainfsub[k][iz];
	       asupsub[j][iz]= asupsub[k][iz];
	       
	       
	     }    

		    break;
		  }
	      }
	  }
      }

    cout<<"qui"<<endl;
   TLine *l=new TLine(0,1,30,1);
   l->SetLineStyle(2);
   l->SetLineWidth(2);

   l->SetLineColor(kBlue-7);
		 
   TCanvas *c[57];
   TH2F *hf[57];


   TLegend *leg[57];
   for(int iz=1;iz<57;iz++)
     {
       hf[iz]=0;
       int io=0;
       for(int k=0;k<nf;k++)
	 {
	   for(int ip=0;ip<10;ip++)
	     {
	   if(ggrec[k][ip][iz]->GetN()>0)
	     {
	       io=1;
	     }
	     }
	 }
       if(io==1)
	 {
	   c[iz]=new TCanvas(Form("c%d",iz),Form("c%d",iz));
	   c[iz]->Draw();
	   float am=1000;
	   float aM=-1000;
	   
	   for(int j=0;j<nf;j++)
	     {
	       if(ainfsub[j][iz]<am && ainfsub[j][iz]>0)am=ainfsub[j][iz];
	       if(asupsub[j][iz]>aM)aM=asupsub[j][iz];
	     }

	   //cout<<iz<<" "<<am<<" "<<aM<<endl;
	   hf[iz]=new TH2F(Form("hf%d",iz),Form("Z=%d",iz),100,am-1,aM+1,100,1e-5,10);
	   hf[iz]->GetXaxis()->SetTitle("A");
	   hf[iz]->GetYaxis()->SetTitle(Form("Y(A)Norm"));
	   hf[iz]->SetStats(kFALSE);
	   hf[iz]->Draw();
	   leg[iz]=new TLegend(0.79,0.65,0.97,0.98);
	   for(int k=0;k<nf;k++)
	     {
	         if(!stit[k].Contains("stiff") &&!stit[k].Contains("soft"))
		 {
		   for(int ip=0;ip<10;ip++)
		     {
		    
		       if( ggrec[k][ip][iz]->GetN()>0)
			 {
	       ggrec[k][ip][iz]->Draw("pl");
			 }
	       gPad->SetLogy(1);
	       gPad->SetGridx(kFALSE);
	       gPad->SetGridy(kFALSE);
		     }
	      
	       leg[iz]->AddEntry(ggrec[k][0][iz],Form("%s", tit[k]),"p");
		 }
	     }
	   leg[iz]->Draw();
	 }

     }
    TCanvas *cv=new TCanvas("cv","cv",0,0,2000,2000);
   cv->Divide(4,4);
   for(int j=1;j<16;j++)
     {
       cv->cd(j);
       if(hf[j+11]!=0)
	 {
       hf[j+11]->Draw();
	   for(int k=0;k<nf;k++)
	     {
	       if(!stit[k].Contains("stiff") &&!stit[k].Contains("soft"))
		 {
		   for(int ip=0;ip<10;ip++)
		     {

if( ggrec[k][ip][j+11]->GetN()>0)
			 {
	       ggrec[k][ip][j+11]->Draw("pl");
			 }
		     }
	       gPad->SetLogy(1);
	       gPad->SetGridx(kFALSE);
	       gPad->SetGridy(kFALSE);
		 }
	     }
leg[j+11]->Draw();       
     }
     } 
}
void tz(TH1D *h,TGraphErrors *g[57],TGraphErrors *gt)
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
  for(int iz=1;iz<57;iz++)
    {
      if(totz[iz]>0)
	{
	  gt->SetPoint(gt->GetN(),(float)iz,totz[iz]);
	  gt->SetPointError(gt->GetN()-1,0.001,sqrt((float)totz[iz]));

	  for(int in=0;in<100;in++)
	    {
	      if(totzn[iz][in]>0)
		{
		  float rr=(float)totzn[iz][in]/(float)totz[iz];
		  g[iz]->SetPoint(g[iz]->GetN(),(float)(iz+in),(float)totzn[iz][in]/(float)totz[iz]);
		  g[iz]->SetPointError(g[iz]->GetN()-1,0.0001,rr*(1/sqrt((float)totzn[iz][in])+1/sqrt((float)totz[iz])));
  }
  }
	}
    }



}
void tz(TH1F *h,TGraphErrors *g[57],TGraphErrors *gt)
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
  for(int iz=1;iz<57;iz++)
    {
      if(totz[iz]>0)
	{
	  gt->SetPoint(gt->GetN(),(float)iz,totz[iz]);
	  gt->SetPointError(gt->GetN()-1,0.001,sqrt((float)totz[iz]));
	  //cout<<iz<<" "<<totz[iz]<<endl;
	  for(int in=0;in<100;in++)
	    {
	      if(totzn[iz][in]>0)
		{

		  float rr=(float)totzn[iz][in]/(float)totz[iz];
		  //cout<<iz<<" "<<totz[iz]<<" "<<in<<" "<<totzn[iz][in]<<" "<<iz+in<<" "<<rr<<endl;
		  g[iz]->SetPoint(g[iz]->GetN(),(float)(iz+in),(float)totzn[iz][in]/(float)totz[iz]);
		  g[iz]->SetPointError(g[iz]->GetN()-1,0.0001,rr*(1/sqrt((float)totzn[iz][in])+1/sqrt((float)totz[iz])));
  }
  }
	}
    }



}
