#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>
void Classe_evento::AnalisiPrincipale()
{
  // CORREZIONE SUL TEMPO DI VOLO DEI PHOS
  //  if(Classe_analisi::Getanalisi()->tipo_analisi>=2)//solo per analisi exp
  // {
  //   for(int j=0;j<evento.moltepl;j++)
  //	{
  //	  if(evento.coderiv[j]<1000) 
  //	    { 
  //	      int ipp=evento.coderiv[j]/10; 
  //	      int ihh=evento.coderiv[j]%10; 
  //	      Correggi_Tempo_Phos(ipp,ihh,j);

  //    }
  //}//fine giro sulle particelle
  // }
  // FINE CORREZIONE SUL TEMPO DI VOLO DEI PHOS

  int flag=0;
  if(Classe_analisi::Getanalisi()->tipo_analisi>=2)//experimental analysis
    {
if(expevent.trig&16)
  {
    flag=1;
  }
    }
  else//MonteCarlo
    {
      flag=1;
    }
  // flag=1;//all triggers
  if(flag==0){return;}
  TH1F*hgamma_mult=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgamma_mult");
  TH1F*hgamma_mult_er=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgamma_mult_er");
  TH1F *hgamma_er=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgamma_er");
 
 TH1F *htbaf=(TH1F*)gROOT->GetListOfSpecials()->FindObject("htbaf");
 TH1F *hgammabaf=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgammabaf");
 TH2F *halphagamma=(TH2F*) gROOT->GetListOfSpecials()->FindObject("halphagamma");
 TH2F *halphaenergy=(TH2F*) gROOT->GetListOfSpecials()->FindObject("halphaenergy");
  TH2F *hprotongamma=(TH2F*) gROOT->GetListOfSpecials()->FindObject("hprotongamma");
 TH2F *hprotonenergy=(TH2F*) gROOT->GetListOfSpecials()->FindObject("hprotonenergy");
  TH1F *hgamma_syfiss=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgamma_syfiss");
  TH1F *hgamma_asyfiss=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hgamma_asyfiss");
		 double EBeam=Classe_analisi::Getanalisi()->reazione.ebeam;
		 double tMin,tMax;
		 if(EBeam<7.5)tMin=70;
		 if(EBeam>8. && EBeam<11)tMin=60;
		 if(EBeam>11)tMin=50;
		 if(EBeam<7.5)tMax=120;
		 if(EBeam>8. && EBeam<11)tMax=95;
		 if(EBeam>11)tMax=70;
  if(flag==1)//good trigger
    {
	int residue=0;
	vector <int> jresidue;
	int heavy=0;
	vector <int> jheavy;
	int lcp=0;
	vector <int> jlcp;
	int imftot=0;
	vector <int> jimftot;
	int nalpha=0;
	vector <int> jalpha;
	int nproton=0;
	vector <int> jproton;
	int ngamma=0;
	vector <int> jgamma;
	vector <int> flagdoppler;
	int nores=0;
	vector <int> jnores;

	for(int j=0;j<evento.moltepl;j++)
	  {
	    if(evento.coderiv[j]>10000)//BaF2
	      {
		if(TMath::Nint(evento.z[j])==0&&TMath::Nint(evento.a[j])==0)
		  {
		    ngamma++;
		    jgamma.push_back(j);
		    if(Classe_analisi::Getanalisi()->tipo_analisi<2)//no doppler correction for Mcarlo
		      {
		    flagdoppler.push_back(1);
		      }
		    else
		      {
flagdoppler.push_back(0);
		      }
		  }
	      }//BaF2

	    //	    if(evento.coderiv[j]<1000)//phos
	    if(evento.coderiv[j]<50)//phos
	      {
		if(TMath::Nint(evento.z[j])==100)
		  {
		    if(evento.tvolo[j]>=tMin && evento.tvolo[j]<=tMax)
		      {
		    residue++;
		    jresidue.push_back(j);
		      }
		  }
		if(evento.z[j]>=50)
		  {
		    heavy++;
		    jheavy.push_back(j);
		  }	       
	      }//phos
	    if(TMath::Nint(evento.z[j])==2)
	      {
		nalpha++;
		jalpha.push_back(j);
	      }
	      if( (TMath::Nint(evento.z[j]) == 1 &&(TMath::Nint(evento.a[j]) ==
	      1))){
		      	nproton++;
			jproton.push_back(j);
	      }
	    if(evento.z[j]>2 && evento.z[j]<50)
	      {
		imftot++;
		jimftot.push_back(j);
	      }

	  }//loop on particle multiplicity
	if(ngamma>0)
	  {
	    for(int j=0;j<ngamma;j++)
	      {
		htbaf->Fill(evento.tvolo[jgamma[j]]);
		if(evento.coderiv[jgamma[j]]==10001)
		  {
		hgammabaf->Fill(evento.epartlab[jgamma[j]]);
		  }
		TH1F *h=(TH1F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",10+
									   evento.coderiv[jgamma[j]]-10000));
		h->Fill(evento.tvolo[jgamma[j]]);
		h=0;
			}
	  }
	if(residue==1 && heavy==1 && imftot==0)
	  {
	    float vv=evento.vpartlab[jresidue[0]];
	hgamma_mult_er->Fill(ngamma);	
	for(int j=0;j<ngamma;j++)
	  {
	    float ed=Classe_formule::doppler(evento.epartlab[jgamma[j]],Classe_analisi::Getanalisi()->reazione.vcm,evento.thetalab[jgamma[j]]);
	    // float ed=Classe_formule::doppler(evento.epartlab[jgamma[j]],evento.vpartlab[jresidue[0]],evento.thetalab[jgamma[j]]);
	    if(evento.coderiv[jgamma[j]]<10008)
	      {
	    evento.epartlab[jgamma[j]]=ed;	
	    flagdoppler[j]=1;
	    hgamma_er->Fill(evento.epartlab[jgamma[j]]);
	    halphagamma->Fill(evento.epartlab[jgamma[j]],nalpha);
	      }
	  }
	
	for(int j=0;j<nalpha;j++)
	  {
	    if(evento.coderiv[jalpha[j]]>1000 && evento.coderiv[jalpha[j]]<9999)
	      {
		int ncsi=(evento.coderiv[jalpha[j]]-1000)%10;
	    halphaenergy->Fill(evento.epartcm[jalpha[j]],ngamma);
		TH1F *h3=(TH1F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",100+ncsi));
		h3->Fill(evento.epartcm[jalpha[j]],ngamma);
		h3=0;
	      }
	  }
	for(int j=0;j<nproton;j++)
	  {
	    if(evento.coderiv[jproton[j]]>1000 && evento.coderiv[jproton[j]]<9999)
	      {
		int ncsi=(evento.coderiv[jproton[j]]-1000)%10;
	    hprotonenergy->Fill(evento.epartcm[jproton[j]],ngamma);
		TH1F *h2=(TH1F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",200+ncsi));
		h2->Fill(evento.epartcm[jproton[j]],ngamma);
		h2=0;
	      }
	  }


	  }//residue condition

	//Asymmetric fission condition
  if((heavy==1 && imftot==1&&(evento.z[jimftot[0]]>2&&evento.z[jimftot[0]]<10))) //selezione fissione asimmetrica PER ANALISI SPECIALE
    {
      float v1[3],v2[3];
      v1[0]=evento.vpartcm_x[jheavy[0]];
      v1[1]=evento.vpartcm_y[jheavy[0]];
      v1[2]=evento.vpartcm_z[jheavy[0]];
      v2[0]=evento.vpartcm_x[jimftot[0]];
      v2[1]=evento.vpartcm_y[jimftot[0]];
      v2[2]=evento.vpartcm_z[jimftot[0]];


      float scp;
      Classe_formule::scaprod(v1,v2,&scp);
      float thetarel=57.296*TMath::ACos(scp/(Classe_formule::modulo(v1)*Classe_formule::modulo(v2)));
     
     
      float phi1=90.-evento.phicm[jheavy[0]];
      if (phi1>180) phi1=phi1-360.;
      float phi2=90.-evento.phicm[jimftot[0]];
      if (phi2>180) phi2=phi2-360.;
      int ip1=evento.coderiv[jheavy[0]]/10;
      int ih1=evento.coderiv[jheavy[0]]-ip1*10;
      int ip2=evento.coderiv[jimftot[0]]/10;
      int ih2=evento.coderiv[jimftot[0]]-ip2*10;     
      float vrel=TMath::Sqrt(pow((v1[0]-v2[0]),2)+pow((v1[1]-v2[1]),2)+pow((v1[2]-v2[2]),2));
     if(thetarel>160 && TMath::Abs(phi1-phi2)>140 && TMath::Abs(phi1-phi2)<220 && vrel>15
	//	 && (Classe_analisi::Getanalisi()->tipo_analisi>=1 ||(Classe_analisi::Getanalisi()->tipo_analisi<1 && mcevent.fissione==1))
	 )
	{	  
	for(int j=0;j<ngamma;j++)
	  {
	    if(evento.coderiv[jgamma[j]]<10008)
	      {
	      hgamma_asyfiss->Fill(evento.epartlab[jgamma[j]]);
	   
	      }
	  }
	}
    }
     //end of asymmetric fission
     //Symmetric fission
     if(heavy==2)
       {
      float v1[3],v2[3];
      v1[0]=evento.vpartcm_x[jheavy[0]];
      v1[1]=evento.vpartcm_y[jheavy[0]];
      v1[2]=evento.vpartcm_z[jheavy[0]];
      v2[0]=evento.vpartcm_x[jheavy[1]];
      v2[1]=evento.vpartcm_y[jheavy[1]];
      v2[2]=evento.vpartcm_z[jheavy[1]];
      float scp;
      Classe_formule::scaprod(v1,v2,&scp);
      float thetarel=57.296*TMath::ACos(scp/(Classe_formule::modulo(v1)*Classe_formule::modulo(v2)));
     

     
      float phi1=90.-evento.phicm[jheavy[0]];
      if (phi1>180) phi1=phi1-360.;
      float phi2=90.-evento.phicm[jheavy[1]];
      if (phi2>180) phi2=phi2-360.;
      int ip1=evento.coderiv[jheavy[0]]/10;
      int ih1=evento.coderiv[jheavy[0]]-ip1*10;
      int ip2=evento.coderiv[jheavy[1]]/10;
      int ih2=evento.coderiv[jheavy[1]]-ip2*10;     
      float vrel=TMath::Sqrt(pow((v1[0]-v2[0]),2)+pow((v1[1]-v2[1]),2)+pow((v1[2]-v2[2]),2));
      if(thetarel>160 && TMath::Abs(phi1-phi2)>140 && TMath::Abs(phi1-phi2)<220 && vrel>15
	 //&& (Classe_analisi::Getanalisi()->tipo_analisi>=1 ||(Classe_analisi::Getanalisi()->tipo_analisi<1 && mcevent.fissione==1)))
	 )    
	{
	for(int j=0;j<ngamma;j++)
	  {
	    if(evento.coderiv[jgamma[j]]<10008)
	      {
	      hgamma_syfiss->Fill(evento.epartlab[jgamma[j]]);
	   
	      }
	  }
	}


	
       }//end of symmetric fission
	  

	    hgamma_mult->Fill(ngamma);

    }//if(flag==1)
hgamma_mult=0;
hgamma_mult_er=0;
hgamma_er=0;
 htbaf=0;
 hgammabaf=0;
}
#endif
void Classe_analisi::RoutineFinale()
{
}
