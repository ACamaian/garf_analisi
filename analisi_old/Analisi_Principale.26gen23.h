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

float Alpha(float *vp1,float *vp2, float a1,float a2) //1 e' il piu' grosso
{
  float vcmc[3];
  for(int j=0;j<3;j++) //cm della coppia 1 2
    {
      vcmc[j]=(a1*vp1[j]+a2*vp2[j])/(a1+a2);
    }

  float vcc[3];
  for(int j=0;j<3;j++)//velocita' del piu' grosso rispetto al cm della coppia
    {
      vcc[j]=vp1[j]-vcmc[j];
    }
  float sca;
  Classe_formule::scaprod(vcc,vcmc,&sca); // prodotto scalare fra velocita' del piu' grosso rispetto al cm della coppia e c.m. della coppia
  float alpha=57.296*TMath::ACos(sca/(Classe_formule::modulo(vcc)*Classe_formule::modulo(vcmc)));
  if(Classe_formule::modulo(vcc)==0)
    {
      cout<<"vcc=0"<<" "<<vp1[0]<<" "<<vp2[0]<<" "<<vp1[1]<<" "<<vp2[1]<<" "<<vp1[2]<<" "<<vp2[2]<<endl;
      cout<<Classe_analisi::Getanalisi()->nentry<<endl;
      return -1000;
 }
  return alpha;
				 
}


float viola(int Z1,int Z2,int A1, int A2)
{
  float a1=A1;
  float z1=Z1;
  float a2=A2;
  float z2=Z2;

  float z;
  float a;
  z=z1+z2;
  a=a1+a2;

  //  float Ek=(0.1189+0.0011)*z*z/pow(a,0.333)+7.3+1.5;
  float Ek=(0.1189)*z*z/pow(a,0.333)+7.3;
float mi=a1*a2/(a1+a2);
 float viola=300*sqrt(2*Ek/(mi*931.5));
 return viola;
}

float pred(float a,float v[3])
{
float pbeam=(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
  float p=a*Classe_formule::amu*Classe_formule::modulo(v)/Classe_formule::cluce;
  p=p/pbeam;
  return p;
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

 int iela=0;  
//Scelta trigger (solo exp)
 int bitpat[8];
/* if(Classe_analisi::Getanalisi()->tipo_analisi>200) */
/*   { */
/*  for(int j=0;j<8;j++) */
/*    { */
/*      bitpat[j]=0; */
/*      int ival0=pow(2,j); */
/*      bitpat[j]=expevent.trig&ival0; */
    
/*      if(bitpat[j]>0) */
/*        { */
/* 	 bitpat[j]=1; */
/* 	 fillh("hbit",(float)j,1); */
/*        } */
 
/*  // se il trigger e' 132, bitpat e' 0 0 1 0 0 0 0 1 (2^2+2^7) */
/*    } */
/*  //if(bitpat[2]==0){return;}//se OR SIRCO strip1-5==0 E OR SIRCO && OR CSI Garf==0 si butta l'evento */
/*   } */
//mixatore(0,2,2,4,-1,-1,-1);
// mixatore(1,1,2,4,1,1,2);
/* for(UInt_t i=0; i<mixing.Erel.at(0).size(); i++){ */
/*   fillh("herel",mixing.Erel.at(0).at(i),1); */
/* 	 } */

/*  for(UInt_t i=0; i<mixing.Prel.at(0).size(); i++){ */
/*    fillh("hprel",Classe_formule::modulo(mixing.Prel.at(0).at(i)),1); */
/* 	 } */





 int ival[100];
float nsuzval[100];

if(Classe_analisi::Getanalisi()->tipo_analisi<200 && Classe_analisi::Getanalisi()->tipo_analisi>=100)//Mc geo
  {
    fillh("hmoltprimari",1.,(float)mcevent.moltprimari,1); 
  }

 //*****************************************
 if(evento.moltepl<=1)
   {
 fillh("hneventi",28.,1);
 //    return;

   }
 //*****************************************
     float asse[3]={0.,0.,1.};

float ztot=0;
float ptotz=0;
 float pbeam=0;
 for(int j=0;j<evento.moltepl;j++)
    {
 ztot=ztot+evento.z[j];
 ptotz=ptotz+Classe_formule::amu*evento.a[j]*evento.vpartlab_z[j]/Classe_formule::cluce;	  
    }
 ptotz=ptotz/(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);

pbeam=(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
 fillh("hztotptotz",ptotz,ztot,1);

  if(ztot>(Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt))return;

  if(ptotz>1.1)return;


 fillh("hneventi",31.,1);
 float vpcm[500][3];


 fillh("hneventi",1.,1.);
 fillh("hmtot",(float)evento.moltepl,1);
 
 if(Classe_analisi::Getanalisi()->tipo_analisi<200)
   {
     fillh("hb",mcevent.par_urto,1);
   }

float vplab[100][3]; 

 
 float zmax=0;
 int malpha=0;
 int md=0;
 int mt=0;

 int mC=0;
 int mH=0;
 int m7Li=0;
 int m6Li=0;

int m12C=0;
int m16O=0;
int m20Ne=0;
int m24Mg=0;
int m28Si=0;
int m32S=0;
int m36Ar=0;

 int m7Be=0;
 int m9Be=0;
 int m10Be=0;
 int m40Ca=0;
 int m48Ca=0;

 int m6He=0;
 int m3He=0;
 int mp=0;
 int m27Al=0;

 int m10B=0;
 int m11B=0;
 int m12B=0;
 int m13C=0;
 int m14N=0;

 int m11C=0;
 
 int jimf[500];

 int j6He[500];

 int j7Li[500];
 int j6Li[500];

 int j12C[500];
 int j11C[500];
 int j13C[500];
 int j16O[500];
 int j20Ne[500];
 int j24Mg[500];
 int j28Si[500];
 int j32S[500];
 int j36Ar[500];

 int jp[500];
 int j27Al[500];

 int j7Be[100];
 int j9Be[100];
 int j10Be[100];

 int j40Ca[100];
 int j48Ca[100];

 int j3He[100];

 int j10B[100];
 int j11B[100];
 int j12B[100];


 int j14N[100];

 int jmax=-1;
 int zlcp=0;
 int nbig=0;
 int nimf=0;
 float vpar,vperp;
 int jalpha[100];
 int jd[100];
 int jt[100];
 int index[500];
 int vec[500];

 for(int j=0;j<500;j++)
   {
     index[j]=-1;
   }

  for(int j=0;j<evento.moltepl;j++)
    {
vec[j]=evento.z[j];
     ival[j]=evento.z[j]*100+evento.a[j]-evento.z[j];
     nsuzval[j]=(float)(evento.a[j]-evento.z[j])/(float)evento.z[j];

      vpcm[j][0]=evento.vpartcm_x[j];
      vpcm[j][1]=evento.vpartcm_y[j];
      vpcm[j][2]=evento.vpartcm_z[j];

      vplab[j][0]=evento.vpartlab_x[j];
      vplab[j][1]=evento.vpartlab_y[j];
      vplab[j][2]=evento.vpartlab_z[j];

//cout<<Classe_geo::Getgeo()->D[TMath::Nint(evento.z[j])][TMath::Nint(evento.a[j])-TMath::Nint(evento.z[j])]<<endl;// Difetto di massa
      Classe_formule::da_xyz_a1(evento.vpartcm_x[j],evento.vpartcm_y[j],evento.vpartcm_z[j],vpcm[j]);

      fillh("hzthe",evento.thetalab[j],evento.z[j],1);
     
 
      fillh("hros",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      fillh("hros2",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      

         fillh("hzv",evento.vpartcm[j],evento.z[j],1);
	

 fillh("hzvlab",evento.vpartlab[j],evento.z[j],1);
 fillh("hzvz",vpcm[j][2],evento.z[j],1);
 //rcocode è idtype cioè PSA, DeltaE -E etc


 if(evento.z[j]>zmax)
   {
zmax=evento.z[j];
 jmax=j;
   }
 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4&&evento.rcoqf[j]==1)
   {
     jalpha[malpha]=j;
malpha++;
   } 
 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==2&&evento.rcoqf[j]==1)
   {
     jd[md]=j;
md++;
   } 
 if(TMath::Nint(evento.z[j])==3 &&TMath::Nint(evento.a[j])==7&&evento.rcoqf[j]==1)
   {
     j7Li[m7Li]=j;
m7Li++;
   } 

 if(TMath::Nint(evento.z[j])==3 &&TMath::Nint(evento.a[j])==6&&evento.rcoqf[j]==1)
   {
     j6Li[m6Li]=j;
m6Li++;
   } 
 if(TMath::Nint(evento.z[j])==6 &&TMath::Nint(evento.a[j])==12&&evento.rcoqf[j]==1)
   {
     j12C[m12C]=j;
m12C++;
   } 
 if(TMath::Nint(evento.z[j])==6 &&TMath::Nint(evento.a[j])==11&&evento.rcoqf[j]==1)
   {
     j11C[m11C]=j;
m11C++;
   }


 
 if(TMath::Nint(evento.z[j])==8 &&TMath::Nint(evento.a[j])==16&&evento.rcoqf[j]==1)
   {
     j16O[m16O]=j;
m16O++;
   } 

 if(TMath::Nint(evento.z[j])==10 &&TMath::Nint(evento.a[j])==20&&evento.rcoqf[j]==1)
   {
     j20Ne[m20Ne]=j;
m20Ne++;
   } 
 if(TMath::Nint(evento.z[j])==12 &&TMath::Nint(evento.a[j])==24&&evento.rcoqf[j]==1)
   {
     j24Mg[m24Mg]=j;
m24Mg++;
   } 
 if(TMath::Nint(evento.z[j])==14 &&TMath::Nint(evento.a[j])==28&&evento.rcoqf[j]==1)
   {
     j28Si[m28Si]=j;
m28Si++;
   } 

 if(TMath::Nint(evento.z[j])==16 &&TMath::Nint(evento.a[j])==32&&evento.rcoqf[j]==1)
   {
     j32S[m32S]=j;
m32S++;
   } 
 if(TMath::Nint(evento.z[j])==18 &&TMath::Nint(evento.a[j])==36&&evento.rcoqf[j]==1)
   {
     j36Ar[m36Ar]=j;
m36Ar++;
   } 

 if(TMath::Nint(evento.z[j])==13 &&TMath::Nint(evento.a[j])==27&&evento.rcoqf[j]==1)
   {
     j27Al[m27Al]=j;
m27Al++;
   } 
 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==1&&evento.rcoqf[j]==1)
   {
     jp[mp]=j;
mp++;
   } 
 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==3&&evento.rcoqf[j]==1)
   {
     jt[mt]=j;
mt++;
   } 

 if(TMath::Nint(evento.z[j])==4 &&TMath::Nint(evento.a[j])==7&&evento.rcoqf[j]==1)
   {
     j7Be[m7Be]=j;
m7Be++;
   }

 if(TMath::Nint(evento.z[j])==4 &&TMath::Nint(evento.a[j])==9&&evento.rcoqf[j]==1)
   {
     j9Be[m9Be]=j;
m9Be++;
   }
 if(TMath::Nint(evento.z[j])==4 &&TMath::Nint(evento.a[j])==10&&evento.rcoqf[j]==1)
   {
     j10Be[m10Be]=j;
m10Be++;
   }

 if(TMath::Nint(evento.z[j])==20 &&TMath::Nint(evento.a[j])==40&&evento.rcoqf[j]==1)
   {
     j40Ca[m40Ca]=j;
m40Ca++;
   }

 if(TMath::Nint(evento.z[j])==20 &&TMath::Nint(evento.a[j])==48&&evento.rcoqf[j]==1)
   {
     j48Ca[m48Ca]=j;
m48Ca++;
   }

 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==6&&evento.rcoqf[j]==1)
   {
     j6He[m6He]=j;
m6He++;
   }

 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==3&&evento.rcoqf[j]==1)
   {
     j3He[m3He]=j;
m3He++;
   }

 if(TMath::Nint(evento.z[j])==5 &&TMath::Nint(evento.a[j])==10&&evento.rcoqf[j]==1)
   {
     j10B[m10B]=j;
m10B++;
   }

 if(TMath::Nint(evento.z[j])==5 &&TMath::Nint(evento.a[j])==11&&evento.rcoqf[j]==1)
   {
     j11B[m11B]=j;
m11B++;
   }

 if(TMath::Nint(evento.z[j])==5 &&TMath::Nint(evento.a[j])==12&&evento.rcoqf[j]==1)
   {
     j12B[m12B]=j;
m12B++;
   }

  if(TMath::Nint(evento.z[j])==6 &&TMath::Nint(evento.a[j])==13&&evento.rcoqf[j]==1)
   {
     j13C[m13C]=j;
m13C++;
   }
  if(TMath::Nint(evento.z[j])==7 &&TMath::Nint(evento.a[j])==14&&evento.rcoqf[j]==1)
   {
     j14N[m14N]=j;
m14N++;
   }

  

if(TMath::Nint(evento.z[j])==6)mC++;
 if(TMath::Nint(evento.z[j])==1)mH++;
 if(TMath::Nint(evento.z[j])==1 || TMath::Nint(evento.z[j])==2)zlcp=zlcp+evento.z[j];
 if(evento.z[j]>=10)nbig++;

 if(evento.z[j]>=3 && evento.z[j]<10)
    {
      jimf[nimf]=j;
nimf++;
    }
 if(evento.rcoqf[j]==1)
   {
 fillh("hivalmolt",(float)ival[j],(float)evento.moltepl,1);
   }

    }//giro su moltepl


     


 TMath::Sort(evento.moltepl,vec,index);


 float ae1,athe2,ae2;
 if(evento.moltepl==1)
   {
     if(TMath::Nint(evento.z[0])==Classe_analisi::Getanalisi()->reazione.zp&&TMath::Nint(evento.a[0])==Classe_analisi::Getanalisi()->reazione.ap &&evento.rcoqf[0]==1)
       {
float Ap=Classe_analisi::Getanalisi()->reazione.ap;
float At=Classe_analisi::Getanalisi()->reazione.at;
 float aein=Classe_analisi::Getanalisi()->reazione.ebeam*Ap/Classe_formule::amu;
 // cout<<Classe_analisi::Getanalisi()->reazione.ebeam<<" "<<Ap<<" "<<Classe_formule::amu<<endl;

 Classe_formule::rel_kin(&Ap,&At,&aein,&evento.thetalab[0],&ae1,&athe2,&ae2);
 if(TMath::Abs(ae1-evento.epartlab[0])<0.01*evento.epartlab[0])iela=1;
       }
   }
 if(evento.moltepl==2)
   {
     if((TMath::Nint(evento.z[index[0]])==Classe_analisi::Getanalisi()->reazione.zp&&TMath::Nint(evento.a[index[0]])==Classe_analisi::Getanalisi()->reazione.ap &&TMath::Nint(evento.z[index[1]])==Classe_analisi::Getanalisi()->reazione.zt&&TMath::Nint(evento.a[index[1]])==Classe_analisi::Getanalisi()->reazione.at)&&evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1)
       {
float Ap=Classe_analisi::Getanalisi()->reazione.ap;
float At=Classe_analisi::Getanalisi()->reazione.at;
 float aein=Classe_analisi::Getanalisi()->reazione.ebeam*Ap/Classe_formule::amu;
 // cout<<Classe_analisi::Getanalisi()->reazione.ebeam<<" "<<Ap<<" "<<Classe_formule::amu<<endl;

 Classe_formule::rel_kin(&Ap,&At,&aein,&evento.thetalab[index[0]],&ae1,&athe2,&ae2);
 if(TMath::Abs(ae1-evento.epartlab[index[0]])<0.01*evento.epartlab[index[0]])iela=1;
       }
   }
 if(iela==1)fillh("hneventi",40.,1);
 float massvec[500],vcmvec[500][3];
 int nbuone=0;
 for(int j=0;j<evento.moltepl;j++)
   {
   
     massvec[nbuone]=evento.a[j];
     for(int k=0;k<3;k++)
       {
	 vcmvec[nbuone][k]=vpcm[j][k];
       }
     nbuone++;
       
   }

float tflow=-1;
 if(nbuone>=2)
   {
Classe_formule::theflow(&tflow,nbuone,massvec,vcmvec);
   }
 if(tflow>0)
   {
     fillh("htflowztot",tflow,ztot,1);
   }
 
 if(evento.moltepl>1)
   {
     if(evento.z[index[0]]>=18 && evento.z[index[1]]>=18)
       {
	 fillh("hvv",vpcm[index[0]][2],vpcm[index[1]][2],1);
       }
   }


 

  float assez[3]={0.,0.,1.};



  fillh("himf",(float)nimf,1);
  fillh("himfzb",(float)nimf,(float)evento.z[index[0]],1);

  fillh("hzindex0",(float)evento.z[index[0]],1);
  float vparrc,vperprc,v1cm,v2cm;
  float vrec[3];
  float pr;
  int iqp=0;
  int aqp=(Classe_analisi::Getanalisi()->reazione.zp)/2;
  //cout<<aqp<<endl;
  if(evento.z[index[0]]>=aqp && vpcm[index[0]][2]>=0)
    {
      if(evento.moltepl==1)
	{
	  iqp=1;
	}
      else
	{
	  if(evento.z[index[1]]<10||(evento.z[index[1]]>=10 && vpcm[index[1]][2]<0))
	    {
      iqp=1;
	    }

	}
    }

      if(evento.z[index[0]]>=aqp)
	{
	  pr=pred(evento.a[index[0]],vplab[index[0]]);
	  fillh("hpz",evento.z[index[0]],pr,1);

	  if(iqp==1)
	    {
	      fillh("hpz_iqp",evento.z[index[0]],pr,1);
	    }

	  if(evento.rcoqf[index[0]]==1)
	    {
	      fillh("hpz_aid",evento.z[index[0]],pr,1);
	      if(iqp==1)
		{
		  fillh("hpz_aid_iqp",evento.z[index[0]],pr,1);

		}
	    }
	}


      if(iqp==1)
	{
	
	  for(int j=0;j<evento.moltepl;j++)
	    {
	      if(j!=index[0]&&evento.rcoqf[j]==1)
		{
		  fillh("hvzivalaltri",(float)ival[j],vpcm[j][2],1);
		  fillh("hvivalaltri",(float)ival[j],Classe_formule::modulo(vpcm[j]),1);
		  fillh("hivalaltrizb",(float)ival[j],evento.z[index[0]],1);

		  if(evento.rcoqf[index[0]]==1)
		    {
		      fillh("hivalaltriivalbig",(float)ival[j],(float)ival[index[0]],1);
		      fillh("hvzivalaltri2",(float)ival[j],vpcm[j][2],1);
		      fillh("hvivalaltri2",(float)ival[j],Classe_formule::modulo(vpcm[j]),1);
		    }


		}
	    }

	}




  for(int j=0;j<evento.moltepl;j++)
    {
if(evento.rcoqf[j]==1)
   {

      fillh("hivalvz",(float)ival[j],vpcm[j][2],1);
      
      if(evento.z[index[0]]>=aqp)
	{
	  int isto=100+(int)(pr*10);
	 

	  if(j!=index[0])
	    {
	     
	  fillh("hivalvz_zbig",(float)ival[j],vpcm[j][2],1);
	  fillh("hivalzbig",(float)ival[j],evento.z[index[0]],1);

	      Classe_formule::vparvperprecoil(vpcm[j],vpcm[index[0]],evento.a[j],evento.a[index[0]],&vparrc,&vperprc,&v1cm,&v2cm,vrec);
	      fillh("hivalvpar",(float)ival[j],vparrc,1);
	      if(vpcm[index[0]][2]>0)
	    {
	      fillh("hivalvparvpos",(float)ival[j],vparrc,1);
	      if(v1cm>0)
		{
		  int npp=0;
		  int npp2=0;
		  int mpp=0;
		  int mpp2=0;
		 
		  fillh("hivalvparvposvcmpos",(float)ival[j],vparrc,1);

		  fillh("hivalmolt2",(float)ival[j],evento.moltepl,1);

		  if(vparrc<0)
		    {
		      fillh("hivalmolt3",(float)ival[j],evento.moltepl,1);
		    }

		  for(int k=0;k<evento.moltepl;k++)
		    {
		      if(k!=j && k!=index[0])
			{
			  
			  npp++;
			  fillh("hivalzaltri",(float)ival[j],evento.z[k],1);
			  if(vparrc<0)
			    {
			      mpp++;
			      fillh("hivalzaltri3",(float)ival[j],evento.z[k],1);
			    }
			  if(vpcm[k][2]>=0)
			    {
			  fillh("hivalzaltrivcmpos",(float)ival[j],evento.z[k],1);
			  //if(evento.z[k]==20)cout<<evento.z[index[0]]<<" "<<vpcm[k][2]<<" "<<vpcm[index[0]][2]<<endl;
			  npp2++;
			  if(vparrc<0)
			    {
			      mpp2++;
			      fillh("hivalzaltrivcmpos3",(float)ival[j],evento.z[k],1);
			    }
			    }

			  if(evento.rcoqf[k]==1)
			    {
			  fillh("hivalivalaltri",(float)ival[j],(float)ival[k],1);
			  if(vparrc<0)
			    {
			       fillh("hivalivalaltri3",(float)ival[j],(float)ival[k],1);
			    }
			  if(vpcm[k][2]>=0)
			    {
			  fillh("hivalivalaltrivcmpos",(float)ival[j],(float)ival[k],1);
			  if(vparrc<0)
			    {
			      fillh("hivalivalaltrivcmpos3",(float)ival[j],(float)ival[k],1);
			    }
			    }
			    }
			}
		    }
		  
		  fillh("hivalnpp",(float)ival[j],(float)npp,1);
		  fillh("hivalmpp",(float)ival[j],(float)mpp,1);
		  if(vparrc<0)
		    {
		  fillh("hivalnpp2",(float)ival[j],(float)npp2,1);
		  
		  fillh("hivalmpp2",(float)ival[j],(float)mpp2,1);
 
		    }
		  if(vparrc>0)
		    {
		      if(Classe_formule::modulo(vrec)-vparrc>0)
			{
		      fillh("hsub",(float)ival[j],-vparrc,-1);
			}
		    }
		  else
		    {
		      fillh("hsub",(float)ival[j],vparrc,1);
		    }


		  if(iqp==1)
		    {
		  fillh("hivalvparvposvcmpos_iqp",(float)ival[j],vparrc,1);
		   fillh(Form("h%d",isto),(float)ival[j],vparrc,1);

		   if(vparrc>0)
		     {
		       if(Classe_formule::modulo(vrec)-vparrc>0)
			{
		       fillh("hsub2",(float)ival[j],-vparrc,-1);
		       fillh(Form("h%d",isto+100),(float)ival[j],-vparrc,-1);
			}
		       
		     }
		   else
		     {
		   fillh(Form("h%d",isto+100),(float)ival[j],vparrc,1);
		   
		   fillh("hsub2",(float)ival[j],vparrc,1);
		     }


		  if(evento.rcoqf[index[0]]==1)
		    {
		  fillh("hivalvparvposvcmpos_iqp_abig",(float)ival[j],vparrc,1);
		  
		  

		  fillh(Form("h%d_aid",isto),(float)ival[j],vparrc,1);
		  

		  if(vparrc>0)
		     {
		       if(Classe_formule::modulo(vrec)-vparrc>0)
			{
		       fillh("hsub3",(float)ival[j],-vparrc,-1);
		       fillh(Form("h%d_aid",isto+100),(float)ival[j],-vparrc,-1);
			}
		     }
		  else
		    {
		      fillh("hsub3",(float)ival[j],vparrc,1);
		      fillh(Form("h%d_aid",isto+100),(float)ival[j],vparrc,1);
		    }

		    }//rcoqf_big=1

		    }//iqp==1

		  
		}//v>vcm
	      fillh("hivalvparcm",(float)ival[j],v1cm,1);
	      if(vparrc>0)
		{
	      fillh("hivalvparcmvparrcpos",(float)ival[j],v1cm,1);
		}

	    }
	    }

	  if(vpcm[index[0]][2]>0)
	    {

	      fillh("hivalvz_zbigvpos",(float)ival[j],vpcm[j][2],1);

	    }


	}

    }

    }
}






#endif
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2);
void nsuz0(TH1F *h,TGraphErrors *gnz);
void Classe_analisi::RoutineFinale()
{
 

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




