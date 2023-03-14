#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>
#define RCO(j) (((evento.coderiv[j]>=1000000)&&(evento.coderiv[j]<10000000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
#define GARF(j) (((evento.coderiv[j]>=1000)&&(evento.coderiv[j]<10000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
//evento.coderiv.push_back(1000000+(isec+1)*100+(istrip+1)*10+icsi+1); RCO (1000111-1000897)
//evento.coderiv.push_back((isec+1)*10+icsi+1+1000); garfield (1011-1248)
//tipo_analisi=210 =>odie

void fillh(string hname,double X,double peso)
{
Classe_formule::fillh(hname,X,peso);
 return;

}
void fillh(string hname,double X,double Y,double peso)
{
  Classe_formule::fillh(hname,X,Y,peso);
  return;
}
TH1F *geth1(string hname)
{
  return Classe_formule::geth1(hname);
}
TH2F *geth2(string hname)
{
  return Classe_formule::geth2(hname);
}
TGraphErrors *getg(string hname)
{
  return Classe_formule::getg(hname);
}

void Classe_evento::AnalisiPrincipale()
{
  //  printf("Evento %lld\n",Classe_analisi::Getanalisi()->nentry);
  //if(Classe_analisi::Getanalisi()->tipo_analisi==0 || Classe_analisi::Getanalisi()->tipo_analisi==100)
  // {  if(mcevent.isresidue!=1)
  //    {
  //	return;
  //    }
  // }
 fillh("hneventi",30.,1);  
//Scelta trigger (solo exp)
 int bitpat[8];
if(Classe_analisi::Getanalisi()->tipo_analisi>200)
  {
 for(int j=0;j<8;j++)
   {
     bitpat[j]=0;
     int ival=pow(2,j);
     bitpat[j]=expevent.trig&ival;
    
     if(bitpat[j]>0)
       {
	 bitpat[j]=1;
	 fillh("hbit",(float)j,1);
       }
 
 // se il trigger e' 132, bitpat e' 0 0 1 0 0 0 0 1 (2^2+2^7)
   }
 if(bitpat[2]==0){return;}//se OR SIRCO strip1-5==0 E OR SIRCO && OR CSI Garf==0 si butta l'evento
  }
//mixatore(0,2,2,4,-1,-1,-1);
// mixatore(1,1,2,4,1,1,2);
/* for(UInt_t i=0; i<mixing.Erel.at(0).size(); i++){ */
/*   fillh("herel",mixing.Erel.at(0).at(i),1); */
/* 	 } */

/*  for(UInt_t i=0; i<mixing.Prel.at(0).size(); i++){ */
/*    fillh("hprel",Classe_formule::modulo(mixing.Prel.at(0).at(i)),1); */
/* 	 } */

 fillh("hneventi",31.,1);

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

if(Classe_analisi::Getanalisi()->tipo_analisi<200 && Classe_analisi::Getanalisi()->tipo_analisi>=100)//Mc geo
  {
    fillh("hmoltprimari",1.,(float)mcevent.moltprimari,1); 
  }

 //*****************************************
 if(evento.moltepl<=1)
   {
 fillh("hneventi",28.,1);
     return;

   }
 //*****************************************



 float vpcm[500][3];


 fillh("hneventi",1.,1.);
 fillh("hmtot",(float)evento.moltepl,1);
  char contfus[100];

int nc0fus=-1;
 if(Classe_analisi::Getanalisi()->tipo_analisi>=200)//exp
  {
 sprintf(contfus,"fus");
  



      

      if(Classe_geo::Getgeo()->ncuts>0)
	{
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	      if(strcmp(contfus,Form("%s",Classe_analisi::Getanalisi()->gggcuts[k]->GetName()))==0)
		{
		  nc0fus=k;
		  break;
		}
	    }	  
	}	  
  }


      int nfus=0;
 

      int lista_fus[10];

  for(int j=0;j<evento.moltepl;j++)
    {
//cout<<Classe_geo::Getgeo()->D[TMath::Nint(evento.z[j])][TMath::Nint(evento.a[j])-TMath::Nint(evento.z[j])]<<endl;// Difetto di massa
      Classe_formule::da_xyz_a1(evento.vpartcm_x[j],evento.vpartcm_y[j],evento.vpartcm_z[j],vpcm[j]);

      fillh("hzthe",evento.thetalab[j],evento.z[j],1);
     
 
      fillh("hros",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      fillh("hros2",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      

         fillh("hzv",evento.vpartcm[j],evento.z[j],1);
 
      //******** FUSIONE
	  if(nc0fus>=0)
	    {
      	      	  if(Classe_analisi::Getanalisi()->gggcuts[nc0fus]->IsInside(evento.vpartcm[j],evento.z[j])==1)
		    {
		     
		      fillh("htest",evento.vpartcm[j],evento.z[j],1);
		        lista_fus[nfus]=j; 
		        nfus++; 
		    }
	    }
	
    
	
      else
	{
	  if(evento.z[j]>16&&evento.vpartcm[j]<10)
	{
	  lista_fus[nfus]=j;
	  nfus++;
	}
	}

	 if(Classe_analisi::Getanalisi()->tipo_analisi<200 && Classe_analisi::Getanalisi()->tipo_analisi>=100)//Mc geo
	   {
	     fillh("he1e2",evento.epartcm[j],mcevent.epartcm[evento.indice_originale[j]],1);

	   }

    }//giro su moltepl

 
}
#endif

void Classe_analisi::RoutineFinale()
{
}
