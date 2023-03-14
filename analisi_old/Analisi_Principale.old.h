#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>
#define RCO(j) (((evento.coderiv[j]>=1000000)&&(evento.coderiv[j]<10000000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
#define GARF(j) (((evento.coderiv[j]>=1000)&&(evento.coderiv[j]<10000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))

void fillh(string hname,double X,double peso) {
	//if(fabs(peso)<0.001) return;
	TH1D *h1temp=(TH1D *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
if(h1temp==0){cout<<hname<<endl;}
	h1temp->Fill(X,peso);
	return;
}

void fillh(string hname,double X,double Y,double peso) {
	//if(fabs(peso)<0.001) return;
	TH2F *h2temp=(TH2F *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
if(h2temp==0){cout<<hname<<endl;}
	h2temp->Fill(X,Y,peso);
	return;
}
void Classe_evento::AnalisiPrincipale()
{

//Scelta trigger (solo exp)
 int bitpat[8];
if(Classe_analisi::Getanalisi()->tipo_analisi>20)
  {
 for(int j=0;j<8;j++)
   {
     bitpat[j]=0;
     int ival=pow(2,j);
     bitpat[j]=expevent.trig&ival;
     if(bitpat[j]>0){bitpat[j]=1;}
   }
 // se il trigger e' 132, bitpat e' 0 0 1 0 0 0 0 1 (2^2+2^7)
  }

  //<<<<<<< Analisi_Principale.h
  //printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);
  //cout<<"molt="<<evento.moltepl<<endl;

 vector <int> isec(evento.moltepl);
 vector <int> istrip(evento.moltepl);
 vector <int> icsi(evento.moltepl);


 for(int j=0;j<evento.moltepl;j++)
    {
      isec[j]=-1;
      istrip[j]=-1;
      icsi[j]=-1;
      if(evento.coderiv[j]>=1000000)
	{
	 isec[j]=(evento.coderiv[j]-1000000)/100;
	 istrip[j]=((evento.coderiv[j]-1000000)-isec[j]*100)/10;
	 icsi[j]=(evento.coderiv[j]-1000000)-isec[j]*100-istrip[j]*10;
		}
      if(evento.coderiv[j]>=1000 && evento.coderiv[j]<1000000)
	{
	  istrip[j]=-1;
	  isec[j]=(evento.coderiv[j]-1000)/10;
	  icsi[j]=(evento.coderiv[j]-1000)-10*isec[j];
	}
    }
int nfusbis=-1;
int lista_fusbis[30];

TH1F*hz=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hz");
TH2F*hros=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hros");
TH2F*hros2=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hros2");
TH1F*hcode=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hcode");
  for(int j=0;j<evento.moltepl;j++)
    {
      //  cout<<"z="<<evento.z[j]<<endl;
      hz->Fill(evento.z[j]);
      
      hros->Fill(evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296));
  hros2->Fill(evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296));
	
      hcode->Fill((float)evento.coderiv[j]-10000000);



      if(Classe_geo::Getgeo()->ncuts>0)
	{
	  int nc0=-1;
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	      if(strcmp("taglio2",Form("%s",Classe_analisi::Getanalisi()->gggcuts[k]->GetName()))==0)
		{
		  nc0=k;
		  break;
		}
	    }
	  if(nc0>=0)
	    {
	  if(Classe_analisi::Getanalisi()->gggcuts[nc0]->IsInside(evento.vpartcm[j],evento.z[j]))
	    {
	      
	      lista_fusbis[nfusbis]=j;
	      nfusbis++;
	    }
	    }
	}
	
    }
  hz=0;
  hros=0;
  hcode=0;
  hros2=0;
  int nheavy=0;
  int indice[100];
  int indice_ori[100];
  int origine[100];
  // cout<<evento.moltepl<<endl;
  if(Classe_analisi::Getanalisi()->tipo_analisi==1 || Classe_analisi::Getanalisi()->tipo_analisi==11)
    {
  for(int j=0;j<evento.moltepl;j++)
    {
      if(evento.z[j]>=10)
	{
	  indice[nheavy]=j;
	  origine[nheavy]=mcevent.origine[evento.indice_originale[j]];
	  indice_ori[nheavy]=evento.indice_originale[j];
	  	  
	  nheavy++;
	  
	}
    }
  if(evento.moltepl>2)
    {
      if(nheavy>2)
	{
  for(int j=0;j<nheavy;j++)
    {
      //      cout<<j<<" "<<nheavy<<" "<<evento.z[indice[j]]<<" "<<origine[j]<<endl;
    }
	}

    }
  //<<<<<<< Analisi_Principale.h
  if(nheavy>2)
    {
  for(int j=0;j<nheavy-1;j++)
    {
      for(int k=j+1;k<nheavy;k++)
	{
	  if(origine[j]==origine[k])
	    {
	      // 	      cout<<Classe_analisi::Getanalisi()->nentry<<" "<<origine[j]<<" "<<origine[k]<<" "<<evento.z[indice[j]]<<" "<<evento.z[indice[k]]<<endl;
		      TH2F *hthecm=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hthecm");
		      hthecm->Fill(evento.thetacm[indice[j]],evento.thetacm[indice[k]],1);
		      hthecm=0;
TH2F *hthe=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hthe");
		      hthe->Fill(evento.thetalab[indice[j]],evento.thetalab[indice[k]],1);
		      hthe=0;
	    }
	}
      
      //cout<<indice_ori[j]<<" "<<evento.z[indice[j]]<<" "<<mcevent.z[indice_ori[j]]<<" "<<origine[j]<<endl;
    }
    }//nheavy>2
  if(nheavy==3)
    {
      float vmin=1000;
      int jmin=-1;
      for(int j=0;j<nheavy;j++)
	{
	  if(evento.vpartlab[indice[j]]<vmin)
	    {
	      vmin=evento.vpartlab[indice[j]];
	      jmin=j;
	    }
	}

 for(int j=0;j<nheavy-1;j++)
    {
      for(int k=j+1;k<nheavy;k++)
	{
	  if(j!=jmin&& k!=jmin && origine[j]!=origine[k])
	    {

TH2F *hthe2=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hthe2");
		      hthe2->Fill(evento.thetalab[indice[j]],evento.thetalab[indice[k]],1);
		      hthe2=0;
	    }
	}

    }
    }//NHEAVY=3

    }//hipse



  if(Classe_analisi::Getanalisi()->tipo_analisi==0 || Classe_analisi::Getanalisi()->tipo_analisi==10)
    {
      if(mcevent.fissione==1)
	{
for(int j=0;j<evento.moltepl;j++)
    {
      if(evento.z[j]>=10)
	{
	  indice[nheavy]=j;
	  origine[nheavy]=mcevent.origine[evento.indice_originale[j]];
	 
	  	  
	  nheavy++;
	  
	}
    }
 if(nheavy>=2)
   {
     for(int j=0;j<nheavy-1;j++)
       {
	 for(int k=j+1;k<nheavy;k++)
	   {
	    float vrel=sqrt(pow((evento.vpartcm_x[indice[j]]-evento.vpartcm_x[indice[k]]),2)+pow((evento.vpartcm_y[indice[j]]-evento.vpartcm_y[indice[k]]),2)+pow((evento.vpartcm_z[indice[j]]-evento.vpartcm_z[indice[k]]),2));
TH2F *hvrel=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hvrel");
	     hvrel->Fill(vrel,1);
	     hvrel=0;
TH2F *hthecm=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hthecm");
		      hthecm->Fill(evento.thetacm[indice[j]],evento.thetacm[indice[k]],1);
		      hthecm=0;
TH2F *hthe=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hthe");
		      hthe->Fill(evento.thetalab[indice[j]],evento.thetalab[indice[k]],1);
		      hthe=0;

	   }
       }

   }



 if(nheavy>2)
   {
     for(int j=0;j<nheavy;j++)
       {
	 cout<<nheavy<<" "<<j<<" "<<evento.z[indice[j]]<<endl;


       }
   }
	}
    }//gemini

  for(int j=0;j<evento.moltepl;j++)
    {
		if(RCO(j)) 
		  {
fillh("hzcode_RCO",evento.z[j],evento.rcocode[j],1.);
		  }
    }
}
#endif
void Classe_analisi::RoutineFinale()
{
}
