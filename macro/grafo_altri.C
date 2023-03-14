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
void nsuz0(TH1D *h,TGraphErrors *gnz);
void grafo_altri()
{

  TFile *fout=new TFile("out19mar.root","RECREATE");

  gStyle->SetGridWidth(0);

  int col[4]={1,2,3,4};
  int mark[4]={20,20,24,25};
  int marke[2]={20,24};
  TFile *f[4];
  char tit[4][100]={"4025","4825","4840","4835"};
  char nome[4][100]={"mar19/totFAZIAPRE","mar19/totFAZIAPRE","mar19/totFAZIAPRE","mar19/totFAZIASYM"};
  float ebeam[4]={25,25,40};
  float ap[4]={40,48,48};
  float at=12.;

  float vplab[4];
  float vpcm[4];
  for(int j=0;j<4;j++)
    {
      vplab[j]=300*sqrt(2*ebeam[j]/931.5);
      vpcm[j]=vplab[j]-ap[j]*vplab[j]/(ap[j]+at);
      cout<<vplab[j]<<" "<<vpcm[j]<<endl;
    }
  TH2F *h100[4][30];
  TH1D *h100pos[4][30];
  TH2F *haltri[4];
  TH1D *haltripos[4];

  TGraphErrors *g[4][30];
  TGraphErrors *gz[4][5];
  TGraphErrors *gbig[4];
  TGraphErrors *gzevap[4];
  TGraphErrors *gzprim[4];
TGraphErrors *grap[4];
  TH1F *hb[4];
  float tot[4][30];
  float totzevap[4][30];
  float zevap[4][30];
  float zprim[4][30];
  for(int j=0;j<4;j++)
    {
	  gbig[j]=new TGraphErrors();
	  gbig[j]->SetMarkerStyle(20);
	  gbig[j]->SetMarkerColor(col[j]);
	  gbig[j]->SetLineColor(col[j]);
 gbig[j]->SetName(Form("%s",tit[j]));

	  gzevap[j]=new TGraphErrors();
	  gzevap[j]->SetMarkerStyle(20);
	  gzevap[j]->SetMarkerColor(col[j]);
	  gzevap[j]->SetLineColor(col[j]);

	  gzprim[j]=new TGraphErrors();
	  gzprim[j]->SetMarkerStyle(20);
	  gzprim[j]->SetMarkerColor(col[j]);
	  gzprim[j]->SetLineColor(col[j]);


      for(int iz=0;iz<5;iz++)
	{
	  gz[j][iz]=new TGraphErrors();
	  gz[j][iz]->SetMarkerStyle(20);
	  gz[j][iz]->SetMarkerColor(col[j]);
	  gz[j][iz]->SetLineColor(col[j]);
	}

  for(int iz=0;iz<30;iz++)
    {
      g[j][iz]=0;
      h100pos[j][iz]=0;
      tot[j][iz]=0;
      totzevap[j][iz]=0;
      zevap[j][iz]=0;
      zprim[j][iz]=0;
    }
    }
  int ny;

  for(int j=0;j<4;j++)
    {
      //            f[j]=new TFile(Form("feb5/totFAZIAPRE_%s.root",tit[j]));
      f[j]=new TFile(Form("%s_%s.root",nome[j],tit[j]));
      //f[j]=new TFile(Form("/home/piantell/prove/test_2020/%s.root",tit[j]));
      haltri[j]=(TH2F*)f[j]->Get("hnzaltrivparmenovbigsololcp");
      haltripos[j]=(TH1D*)haltri[j]->ProjectionY(Form("haltripos%d",j),haltri[j]->GetXaxis()->FindBin(0.),haltri[j]->GetXaxis()->GetNbins());


      hb[j]=(TH1F*)f[j]->Get("hnzbiggestconlcp");
      nsuz0(hb[j],gbig[j]);
     
      for(int iz=10;iz<26;iz++)
	{
	  tot[j][iz]=hb[j]->Integral(hb[j]->GetXaxis()->FindBin(iz*100),hb[j]->GetXaxis()->FindBin(iz+1*100)-1);
	}


       for(int iz=10;iz<26;iz++)
	 {
	   h100[j][iz]=(TH2F*)f[j]->Get(Form("h%d",100+iz));
	   
	   h100pos[j][iz]=(TH1D*)h100[j][iz]->ProjectionY(Form("z%d_%d",iz,j),h100[j][iz]->GetXaxis()->FindBin(0.),h100[j][iz]->GetXaxis()->GetNbins());
	   if(h100pos[j][iz]!=0)
	     {
	       g[j][iz]=new TGraphErrors();
	       g[j][iz]->SetMarkerStyle(iz+10);
	       g[j][iz]->SetName(Form("g%d_%d",iz,j));
	   nsuz0(h100pos[j][iz],g[j][iz]);

	   if(g[j][iz]->GetN()>0)
	     {
	       double *x=g[j][iz]->GetX();
	       double *y=g[j][iz]->GetY();
	       double *ey=g[j][iz]->GetEY();
	       for(int k=0;k<g[j][iz]->GetN();k++)
		 {
		   gz[j][(int)x[k]]->SetPoint(gz[j][(int)x[k]]->GetN(),(float)iz,y[k]);
		   gz[j][(int)x[k]]->SetPointError(gz[j][(int)x[k]]->GetN()-1,0.001,ey[k]);
		 }

	     }

	   if(tot[j][iz]>10)
	     {
	   for(int k=0;k<h100pos[j][iz]->GetXaxis()->GetNbins();k++)
	     {
	       float icount=h100pos[j][iz]->GetBinContent(k+1);
	       
	       int izz=h100pos[j][iz]->GetXaxis()->GetBinLowEdge(k+1)/100;
	  
	       if(icount>0)
		 {
		   //   cout<<k+1<<" "<<h100pos[j][iz]->GetXaxis()->GetBinLowEdge(k+1)<<" "<<izz<<endl;
	       totzevap[j][iz]=totzevap[j][iz]+(float)izz*icount;
		 }
	     }
	   zevap[j][iz]=2*totzevap[j][iz]/tot[j][iz];
	   float e=2*zevap[j][iz]*sqrt(1/totzevap[j][iz]+1/tot[j][iz]);

	   zprim[j][iz]=zevap[j][iz]+(float)iz;
	  
	   gzevap[j]->SetPoint(gzevap[j]->GetN(),(float)iz,zevap[j][iz]);
	   gzprim[j]->SetPoint(gzprim[j]->GetN(),(float)iz,zprim[j][iz]);


 gzevap[j]->SetPointError(gzevap[j]->GetN()-1,0.001,e);
gzprim[j]->SetPointError(gzprim[j]->GetN()-1,0.001,e);
	 
	     }

	     }
	 }

    }




  for(int j=0;j<4;j++)
    { 
grap[j]=new TGraphErrors();
      
       grap[j]->SetMarkerStyle(mark[j]);
       grap[j]->SetMarkerColor(col[j]);
       grap[j]->SetLineColor(col[j]);

        grap[j]->SetMarkerSize(3);
       grap[j]->SetPoint(grap[j]->GetN(),1.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(101))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(101))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(101)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(101))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100))));
       grap[j]->SetPoint(grap[j]->GetN(),2.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(102)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(100)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(100))));

       grap[j]->SetPoint(grap[j]->GetN(),3.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(201)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(201)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(102)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(102))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(201)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(201))));


       grap[j]->SetPoint(grap[j]->GetN(),4.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(204))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(202)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(204))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(202)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(204)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(204))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(202)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(202))));
       grap[j]->SetPoint(grap[j]->GetN(),5.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(305))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(303)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(305))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(303)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(305)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(305))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(303)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(303))));
       grap[j]->SetPoint(grap[j]->GetN(),6.5,haltripos[j]->GetBinContent(haltripos[j]->FindBin(304))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(403)));
       grap[j]->SetPointError(grap[j]->GetN()-1,0.0001,(haltripos[j]->GetBinContent(haltripos[j]->FindBin(304))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(403)))*(sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(304)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(304))+sqrt(haltripos[j]->GetBinContent(haltripos[j]->FindBin(403)))/haltripos[j]->GetBinContent(haltripos[j]->FindBin(403))));


    }






  TCanvas *c[4];
  TH2F *h[4];
  TLegend *leg[4];
  for(int j=0;j<4;j++)
    {
      c[j]=new TCanvas(Form("c%d",j),Form("c%d",j));
      c[j]->Draw();
      h[j]=new TH2F(Form("h%d",j),Form("%s",tit[j]),10,0,10,100,0.,1.4);
      h[j]->Draw();
      h[j]->GetXaxis()->SetTitle("Zaltri");
      h[j]->GetYaxis()->SetTitle("N/Z altri");

      leg[j]=new TLegend(0.5,0.3,0.9,0.9);
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);

      for(int iz=10;iz<26;iz++)
	{
	  if(g[j][iz]->GetN()>0)
	    {
	      g[j][iz]->SetMarkerSize(2);
	  g[j][iz]->Draw("p");
	  
	  leg[j]->AddEntry(g[j][iz],Form("z%d",iz),"p");
	    }
	}
      leg[j]->Draw();
      c[j]->Print(Form("nsuzaltrivszp%d.pdf",j));

    }
  TCanvas *cz[4];
 TH2F *hz[4];
  for(int k=0;k<4;k++)
    {
      cz[k]=new TCanvas(Form("cz%d",k+1),Form("cz%d",k+1));
      cz[k]->Draw();
      if(k==0)hz[k]=new TH2F(Form("hz%d",k+1),Form("Z=%d",k+1),30,9,24,100,0.,1);
      if(k==1)hz[k]=new TH2F(Form("hz%d",k+1),Form("Z=%d",k+1),30,9,24,100,0.94,1.02);
      if(k==2)hz[k]=new TH2F(Form("hz%d",k+1),Form("Z=%d",k+1),30,9,24,100,1.,1.4);
      if(k==3)hz[k]=new TH2F(Form("hz%d",k+1),Form("Z=%d",k+1),30,9,24,100,0.8,1.4);
      hz[k]->GetXaxis()->SetTitle("Zbig");
      hz[k]->GetYaxis()->SetTitle(Form("N/Z altri per Z=%d",k+1));

      hz[k]->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);

      for(int j=0;j<4;j++)
	{
	  gz[j][k+1]->Draw("pl");
	}
      cz[k]->Print(Form("zaltri%d.pdf",k+1));
    }


  TCanvas *cbig=new TCanvas("cbig","cbig");
  TH2F *hbig=new TH2F("hbig","hbig",30,9,25,100,0.8,1.4);
  hbig->GetXaxis()->SetTitle("zbig");
      cbig->Draw();
		      hbig->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);

      for(int j=0;j<4;j++)
	{
	  gbig[j]->Draw("pl");
	  fout->Add(gbig[j]);
	}
      cbig->Print("zbig.pdf");
    

  TCanvas *czevap=new TCanvas("czevap","czevap");
  TH2F *hzevap=new TH2F("hzevap","hzevap",30,9,25,100,0.001,2);
  hzevap->GetXaxis()->SetTitle("zbig");
      czevap->Draw();
		      hzevap->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);

      for(int j=0;j<4;j++)
	{
	  gzevap[j]->Draw("pl");
	}
      gPad->SetLogy(kTRUE);
      czevap->Print("zevap.pdf");
    

  TCanvas *czprim=new TCanvas("czprim","czprim");
  TH2F *hzprim=new TH2F("hzprim","hzprim",30,9,25,100,0,30);
  hzprim->GetXaxis()->SetTitle("zbig");
      czprim->Draw();
		      hzprim->Draw();
      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);

      for(int j=0;j<4;j++)
	{
	  gzprim[j]->Draw("pl");
	}
      czprim->Print("zprim.pdf");
  

  TCanvas *c2r=new TCanvas("c2r","c2r",0,0,1200,600);
  c2r->Draw();


    TH2F *h2r=new TH2F("h2r","",20,0.5,7,100,1e-4,30);

     h2r->GetXaxis()->SetLabelSize(0.15);
 h2r->GetXaxis()->SetLabelOffset(0.03);
      h2r->GetYaxis()->SetLabelSize(0.1);
      h2r->GetXaxis()->SetTitleSize(0.08);
      h2r->GetYaxis()->SetTitleSize(0.12);
      h2r->SetTitle("");
       h2r->SetStats(kFALSE);
      h2r->GetYaxis()->SetTitle("Yield ratio");
      h2r->GetXaxis()->SetTitle("");
      h2r->GetXaxis()->SetNdivisions(505);
  h2r->GetYaxis()->SetNdivisions(505);
    h2r->GetXaxis()->SetTitleOffset(0.9);
    h2r->GetYaxis()->SetTitleOffset(0.6);

//     h2r->GetYaxis()->SetRangeUser(0,60);

  h2r->Draw();
 
    gPad->SetLogy(kTRUE);
  

  h2r->GetXaxis()->SetBinLabel(h2r->FindBin(1.5)-1,"d/p");
  h2r->GetXaxis()->SetBinLabel(h2r->FindBin(2.5)-1,"t/p");
  h2r->GetXaxis()->SetBinLabel(h2r->FindBin(3.5)-1,"t/^{3}He");
  h2r->GetXaxis()->SetBinLabel(h2r->FindBin(4.5)-1,"^{6}He/#alpha");
  h2r->GetXaxis()->SetBinLabel(h2r->FindBin(5.5)-1,"^{8}Li/^{6}Li");
 h2r->GetXaxis()->SetBinLabel(h2r->FindBin(6.5)-1,"^{7}Li/^{7}Be");


  for(int j=0;j<4;j++)
    {
      
      grap[j]->Draw("pl");
 

      gPad->SetGridx(kFALSE);
      gPad->SetGridy(kFALSE);
      gPad->SetBottomMargin(0.2);
            gPad->SetLeftMargin(0.14);
	    gPad->SetTopMargin(0.11);
	    //gPad->SetRightMargin(0.15);
	     gPad->Modified();
      gPad->Update();

    }
  
 
  c2r->Print("rap.pdf");

  fout->Write();
  fout->Close();

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
 
 //cout<<gnz->GetName()<<endl;
 gnz->GetXaxis()->SetTitle("Z");
 gnz->GetYaxis()->SetTitle("N/Z");

}

void nsuz0(TH1D *h,TGraphErrors *gnz)
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
 
 //cout<<gnz->GetName()<<endl;
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
 
 //cout<<gnz->GetName()<<endl;
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
