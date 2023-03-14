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

  float Ek=(0.1189+0.0011)*z*z/pow(a,0.333)+7.3+1.5;
float mi=a1*a2/(a1+a2);
 float viola=300*sqrt(2*Ek/(mi*931.5));
 return viola;
}


void Classe_evento::AnalisiPrincipale()
{
   float Delta[30][60];
   Delta[1][2]=13.136;
   Delta[1][3]=14.9498;
   Delta[1][1]=7.289;
   Delta[2][4]=2.425;
   Delta[2][6]=17.594;

   Delta[2][3]=14.931;


   Delta[6][12]=0;
   Delta[3][5]=11.68;
   Delta[3][9]=24.9539;
   Delta[4][8]=4.942;

   Delta[5][9]=12.416;
   Delta[5][10]=12.051;
   Delta[5][11]=8.668;
   Delta[3][7]=14.908;
   Delta[3][6]=14.086;


   Delta[8][16]=-4.737;
   Delta[10][20]=-7.042;
   Delta[12][24]=-13.933;
   Delta[13][27]=-17.197;

   Delta[14][28]=-21.493;
   Delta[16][32]=-26.016;
   Delta[18][36]=-30.230;
   Delta[20][40]=-34.846;
   Delta[20][48]=-44.215;
   Delta[4][7]=15.769;
   Delta[4][9]=11.348;
   Delta[4][10]=12.607;

   Delta[6][11]=10.650;
   Delta[6][13]=3.125;
   Delta[6][14]=3.020;

   Delta[22][44]=-37.548;
   Delta[21][41]=-28.642;

   Delta[22][52]=-49.464;
   Delta[21][49]=-46.552;

   Delta[3][4]=25.3;


   Delta[8][15]=2.85;
   Delta[7][14]=2.86;
   Delta[5][12]=13.37;


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


   for(int j=0;j<evento.moltepl;j++)
   {
     if(TMath::Nint(evento.z[j])==17&&TMath::Nint(evento.a[j])==35&& evento.rcoqf[j]==1)
       {
	 fillh("he_17_35",faziakali.esi1[evento.indice_originale[j]],(float)evento.coderiv[j]-10000000,1);
       }
   }


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
 for(int j=0;j<evento.moltepl;j++)
    {
 ztot=ztot+evento.z[j];
 ptotz=ptotz+Classe_formule::amu*evento.a[j]*evento.vpartlab_z[j]/Classe_formule::cluce;	  
    }
 ptotz=ptotz/(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);

 fillh("hztotptotz",ptotz,ztot,1);

  if(ztot>(Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt))return;

  if(ptotz>1.1)return;


 fillh("hneventi",31.,1);
 float vpcm[500][3];


 fillh("hneventi",1.,1.);
 fillh("hmtot",(float)evento.moltepl,1);
 


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


 int j6He[500];

 int j7Li[500];
 int j6Li[500];

 int j12C[500];
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

 int j13C[100];
 int j14N[100];

 int jmax=-1;
 int zlcp=0;
 int nbig=0;
 float vpar,vperp;
 int jalpha[100];
 int jd[100];
 int jt[100];
 int index[500];
 int vec[500];
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
	 if(evento.rcocode[j]==11)fillh("hzv11",evento.vpartcm[j],evento.z[j],1);
	 if(evento.rcocode[j]==12)fillh("hzv12",evento.vpartcm[j],evento.z[j],1);
	 if(evento.rcocode[j]==23)fillh("hzv23",evento.vpartcm[j],evento.z[j],1);
	 if(evento.rcocode[j]==33)fillh("hzv33",evento.vpartcm[j],evento.z[j],1);

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

 fillh("hivalmolt",(float)ival[j],(float)evento.moltepl,1);


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
 

  float assez[3]={0.,0.,1.};
  fillh("hztotzmax",ztot,zmax,1);

  fillh("hzmaxmalpha",zmax,(float)malpha,1);
  fillh("hzmaxmC",zmax,(float)mC,1);
  fillh("hzmaxmH",zmax,(float)mH,1);
 fillh("hzmaxzlcp",zmax,(float)zlcp,1);

 fillh("hmalpha",(float)malpha,1);
 fillh("hnbig",(float)nbig,1);

 if(jmax>=0)
   {
    if(vpcm[jmax][2]>0&&zmax>=10 && nbig==1&&iela==0)
   {
     if(evento.rcoqf[jmax]==1)
       {
fillh("hnzbig",(float)ival[jmax],1);
 fillh("hnzbigriv",(float)ival[jmax],(float)evento.coderiv[jmax]-10000000,1);
 fillh("hnzbigvl",(float)ival[jmax],evento.vpartlab[jmax],1);
 fillh("hnzbigvcm",(float)ival[jmax],evento.vpartcm[jmax],1);

 fillh("hnzbigvlnorm",(float)ival[jmax],evento.vpartlab[jmax]/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzbigvcmnorm",(float)ival[jmax],evento.vpartcm[jmax]/Classe_analisi::Getanalisi()->reazione.vp_cm,1);


 fillh("hnzbige",(float)ival[jmax],evento.epartlab[jmax],1);
 // cout<<evento.epartlab[jmax]<<endl;
       }
     //  if(TMath::Nint(evento.z[jmax])==20)cout<<evento.a[jmax]<<" "<<evento.rcoqf[jmax]<<" "<<evento.coderiv[jmax]<<endl;

     fillh("hzvpcmbig",evento.vpartcm[jmax],(float)zmax,1);
     for(int k=0;k<malpha;k++)
       {
     Classe_formule::vparvperp(vpcm[jalpha[k]],assez,&vpar,&vperp);
     fillh("hocchialphafascio",vpar,vperp,1);
     Classe_formule::vparvperp(vpcm[jalpha[k]],vpcm[jmax],&vpar,&vperp);
     fillh("hocchialpha",vpar,vperp,1);
       }
   /*   if(malpha>1) */
/*        { */
/* 	 fillh("herel",erel,1); */
/*        } */
     
   }
   }


  if(evento.moltepl>=2)
    {
      if(vpcm[index[0]][2]>0 && evento.z[index[0]]>=10 && iela==0)
	{
if(evento.rcoqf[index[0]]==1)
   {
fillh("hnz30",(float)ival[index[0]],1);
 if(evento.z[index[1]]<5)
   {
 for(int j=1;j<evento.moltepl;j++)
   {
     if(evento.rcoqf[index[j]]==1)
       {
	 fillh("hnz31",(float)ival[index[j]],1);
       }
   }

fillh("hnz32",(float)ival[index[0]],1);

   }
   }
	}
    }

 float thetarel;
 float vrel[3];

  if(evento.moltepl>1)
    {

      if(evento.z[index[1]]>=5)
	{
  thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);
		      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);

		      fillh("hthetarelvrel_nosel",Classe_formule::modulo(vrel),thetarel,1);

	  float etaz=(float)(evento.z[index[0]]-evento.z[index[1]])/(evento.z[index[0]]+evento.z[index[1]]);

	  fillh("hetaz",etaz,1);

		    

							     
	  fillh("hz1z2",evento.z[index[0]],evento.z[index[1]],1);
	  fillh("hv1v2lab",evento.vpartlab[index[0]],evento.vpartlab[index[1]],1);
	  fillh("hv1zv2zcm",evento.vpcm[index[0]][2],evento.vpcm[index[1]][2],1);
	  fillh("hz1v1lab",evento.vpartlab[index[0]],evento.z[index[0]],1);
	  fillh("hz2v2lab",evento.vpartlab[index[1]],evento.z[index[1]],1);
		fillh("hz1v1zcm",evento.vpcm[index[0]][2],evento.z[index[0]],1);
		      fillh("hz2v2zcm",evento.vpcm[index[1]][2],evento.z[index[1]],1);
		      float ztt,att,ivt;
		      float etaa;
		      int ieta=etaz*10;
int ieta_a;
 float mi,tke;
 float vvp[3];
float vparc,vperpc;
 float vpart,vperpt;
 float vpars,vperps;
if(evento.rcoqf[index[1]]==1&&evento.rcoqf[index[0]]==1)
  {
	 ztt=evento.z[index[0]]+evento.z[index[1]];
	 att=evento.a[index[0]]+evento.a[index[1]];
	 ivt=100*ztt+att-ztt;
 etaa=(float)(evento.a[index[0]]-evento.a[index[1]])/(evento.a[index[0]]+evento.a[index[1]]);
if(etaa<0)etaa=-etaa;
ieta_a=etaa*10; 

	 mi=(float)evento.a[index[0]]*(float)evento.a[index[1]]/(evento.a[index[0]]+evento.a[index[1]]);
	 tke=0.5*mi*931.5*pow(Classe_formule::modulo(vrel),2)/pow(300,2);
	 for(int k=0;k<3;k++)
	   {
	     vvp[k]=((float)evento.a[index[0]]*vpcm[index[0]][k]+(float)evento.a[index[1]]*vpcm[index[1]][k])/att;
	   }
 }

if(evento.rcoqf[index[1]]==1&&evento.rcoqf[index[0]]==1)
  {
	 fillh("hindex1index2",(float)ival[index[0]],(float)ival[index[1]],1);

	 
	 fillh("hindextot",ivt,1);
	 fillh("hvrelindextot",Classe_formule::modulo(vrel),ivt,1);
	 fillh("hvnormindextot",Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);

	
	 fillh(Form("h%d",80+ieta),Classe_formule::modulo(vrel),ivt,1);
	 fillh(Form("h%d",180+ieta),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);

	 fillh("htkeindextot",tke,ivt,1);
	   fillh(Form("h%d",60+ieta),tke,ivt,1);
	
	 
	 

	 fillh(Form("h%d",30+ieta_a),tke,ivt,1);

	 fillh(Form("h%d",(int)ivt+4000),(float)ival[index[0]],Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vrel_zt%d_%d",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vnorm_zt%d_%d",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),1);

	




	 fillh("hzvricos",Classe_formule::modulo(vvp),ztt,1);
	 fillh("hzvzricos",vvp[2],ztt,1);

	 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
	   {
	     fillh("hzvricoslt90",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90",vvp[2],ztt,1);
	   }
	 else
	   {
	     fillh("hzvricosgt90",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90",vvp[2],ztt,1);
	   }
	 if((int)ztt==20 && ieta>=2)
	   {
 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
   {
	     fillh("hzvricoslt90z20",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90z20",vvp[2],ztt,1);

	     fillh("hzvzbiglt90z20",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmalllt90z20",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }
 else
   {

     fillh("hzvricosgt90z20",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90z20",vvp[2],ztt,1);

	     fillh("hzvzbiggt90z20",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmallgt90z20",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }
 

	     for(int m=2;m<evento.moltepl;m++)
	       {

		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[0]],&vparc,&vperpc);
		    Classe_formule::vparvperp(vpcm[index[m]],vvp,&vpart,&vperpt);
		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[1]],&vpars,&vperps);

		 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
		   {
		 fillh("hocchialt90",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);
		 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchialt90p",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

 fillh("hocchilt90p_zbig",vparc,vperpc,1);
 fillh("hocchilt90p_zt",vpart,vperpt,1);
 fillh("hocchilt90p_zs",vpars,vperps,1);

		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchialt90a",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

		    fillh("hocchilt90a_zbig",vparc,vperpc,1);
fillh("hocchilt90a_zt",vpart,vperpt,1);
fillh("hocchilt90a_zs",vpars,vperps,1);
		   }


		   }
		 else
		   {
fillh("hocchiagt90",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);

 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchiagt90p",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 
fillh("hocchigt90p_zbig",vparc,vperpc,1);
fillh("hocchigt90p_zt",vpart,vperpt,1);
fillh("hocchigt90p_zs",vpars,vperps,1);
		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchiagt90a",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

fillh("hocchigt90a_zbig",vparc,vperpc,1);
fillh("hocchigt90a_zt",vpart,vperpt,1);
fillh("hocchigt90a_zs",vpars,vperps,1);
		   }		   }
	       }

	   }

  }




		         if(thetarel>=120)		      
			{


	  fillh("hetazgt120",etaz,1);

		    

							     
	  fillh("hz1z2gt120",evento.z[index[0]],evento.z[index[1]],1);
	  fillh("hv1v2labgt120",evento.vpartlab[index[0]],evento.vpartlab[index[1]],1);
	  fillh("hv1zv2zcmgt120",evento.vpcm[index[0]][2],evento.vpcm[index[1]][2],1);
	  fillh("hz1v1labgt120",evento.vpartlab[index[0]],evento.z[index[0]],1);
	  fillh("hz2v2labgt120",evento.vpartlab[index[1]],evento.z[index[1]],1);
		fillh("hz1v1zcmgt120",evento.vpcm[index[0]][2],evento.z[index[0]],1);
		      fillh("hz2v2zcmgt120",evento.vpcm[index[1]][2],evento.z[index[1]],1);


if(evento.rcoqf[index[1]]==1&&evento.rcoqf[index[0]]==1)
  {
	 fillh("hindex1index2gt120",(float)ival[index[0]],(float)ival[index[1]],1);

	 
	 fillh("hindextotgt120",ivt,1);
	 fillh("hvrelindextotgt120",Classe_formule::modulo(vrel),ivt,1);
	 fillh("hvnormindextotgt120",Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);

	 fillh(Form("h%dgt120",80+ieta),Classe_formule::modulo(vrel),ivt,1);
	 fillh(Form("h%dgt120",180+ieta),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);



	 fillh("htkeindextotgt120",tke,ivt,1);
	   fillh(Form("h%dgt120",60+ieta),tke,ivt,1);


	 fillh(Form("h%dgt120",30+ieta_a),tke,ivt,1);

	 fillh(Form("h%dgt120",(int)ivt+4000),(float)ival[index[0]],Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vrel_zt%d_%dgt120",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vnorm_zt%d_%dgt120",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),1);




	 fillh("hzvricosgt120",Classe_formule::modulo(vvp),ztt,1);
	 fillh("hzvzricosgt120",vvp[2],ztt,1);
	 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
	   {
	     fillh("hzvricoslt90gt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90gt120",vvp[2],ztt,1);
	   }
	 else
	   {
	     fillh("hzvricosgt90gt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90gt120",vvp[2],ztt,1);
	   }
	 if((int)ztt==20 && ieta>=2)
	   {
 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
   {
	     fillh("hzvricoslt90z20gt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90z20gt120",vvp[2],ztt,1);

	     fillh("hzvzbiglt90z20gt120",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmalllt90z20gt120",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }
 else
   {

     fillh("hzvricosgt90z20gt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90z20gt120",vvp[2],ztt,1);

	     fillh("hzvzbiggt90z20gt120",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmallgt90z20gt120",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }


	     for(int m=2;m<evento.moltepl;m++)
	       {

		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[0]],&vparc,&vperpc);
		    Classe_formule::vparvperp(vpcm[index[m]],vvp,&vpart,&vperpt);
		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[1]],&vpars,&vperps);

		 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
		   {
		 fillh("hocchialt90gt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);
		 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchialt90pgt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

 fillh("hocchilt90p_zbiggt120",vparc,vperpc,1);
 fillh("hocchilt90p_ztgt120",vpart,vperpt,1);
 fillh("hocchilt90p_zsgt120",vpars,vperps,1);

		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchialt90agt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

		    fillh("hocchilt90a_zbiggt120",vparc,vperpc,1);
fillh("hocchilt90a_ztgt120",vpart,vperpt,1);
fillh("hocchilt90a_zsgt120",vpars,vperps,1);
		   }


		   }
		 else
		   {
fillh("hocchiagt90gt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);

 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchiagt90pgt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 
fillh("hocchigt90p_zbiggt120",vparc,vperpc,1);
fillh("hocchigt90p_ztgt120",vpart,vperpt,1);
fillh("hocchigt90p_zsgt120",vpars,vperps,1);
		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchiagt90agt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

fillh("hocchigt90a_zbiggt120",vparc,vperpc,1);
fillh("hocchigt90a_ztgt120",vpart,vperpt,1);
fillh("hocchigt90a_zsgt120",vpars,vperps,1);
		   }		   }
	       }

	   }

  }




			}//thetarel>=120




		         if(thetarel<120)		      
			{


	  fillh("hetazlt120",etaz,1);

		    

							     
	  fillh("hz1z2lt120",evento.z[index[0]],evento.z[index[1]],1);
	  fillh("hv1v2lablt120",evento.vpartlab[index[0]],evento.vpartlab[index[1]],1);
	  fillh("hv1zv2zcmlt120",evento.vpcm[index[0]][2],evento.vpcm[index[1]][2],1);
	  fillh("hz1v1lablt120",evento.vpartlab[index[0]],evento.z[index[0]],1);
	  fillh("hz2v2lablt120",evento.vpartlab[index[1]],evento.z[index[1]],1);
		fillh("hz1v1zcmlt120",evento.vpcm[index[0]][2],evento.z[index[0]],1);
		      fillh("hz2v2zcmlt120",evento.vpcm[index[1]][2],evento.z[index[1]],1);


if(evento.rcoqf[index[1]]==1&&evento.rcoqf[index[0]]==1)
  {
	 fillh("hindex1index2lt120",(float)ival[index[0]],(float)ival[index[1]],1);

	 
	 fillh("hindextotlt120",ivt,1);
	 fillh("hvrelindextotlt120",Classe_formule::modulo(vrel),ivt,1);
	 fillh("hvnormindextotlt120",Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);

	 fillh(Form("h%dlt120",80+ieta),Classe_formule::modulo(vrel),ivt,1);
	 fillh(Form("h%dlt120",180+ieta),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),ivt,1);



	 fillh("htkeindextotlt120",tke,ivt,1);
	   fillh(Form("h%dlt120",60+ieta),tke,ivt,1);


	 fillh(Form("h%dlt120",30+ieta_a),tke,ivt,1);

	 fillh(Form("h%dlt120",(int)ivt+4000),(float)ival[index[0]],Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vrel_zt%d_%dlt120",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel),1);

	 fillh(Form("halpha_vnorm_zt%d_%dlt120",(int)ztt,ieta),Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]),Classe_formule::modulo(vrel)/viola((int)evento.z[index[0]],(int)evento.z[index[1]],(int)evento.a[index[0]],(int)evento.a[index[1]]),1);



	 fillh("hzvricoslt120",Classe_formule::modulo(vvp),ztt,1);
	 fillh("hzvzricoslt120",vvp[2],ztt,1);
	 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
	   {
	     fillh("hzvricoslt90lt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90lt120",vvp[2],ztt,1);
	   }
	 else
	   {
	     fillh("hzvricosgt90lt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90lt120",vvp[2],ztt,1);
	   }
	 if((int)ztt==20 && ieta>=2)
	   {
 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
   {
	     fillh("hzvricoslt90z20lt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricoslt90z20lt120",vvp[2],ztt,1);

	     fillh("hzvzbiglt90z20lt120",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmalllt90z20lt120",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }
 else
   {

     fillh("hzvricosgt90z20lt120",Classe_formule::modulo(vvp),ztt,1);
	     fillh("hzvzricosgt90z20lt120",vvp[2],ztt,1);

	     fillh("hzvzbiggt90z20lt120",vpcm[index[0]][2],(float)evento.z[index[0]],1); 
	     fillh("hzvzsmallgt90z20lt120",vpcm[index[1]][2],(float)evento.z[index[1]],1); 

   }
 float vparc,vperpc;
 float vpart,vperpt;
 float vpars,vperps;

	     for(int m=2;m<evento.moltepl;m++)
	       {

		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[0]],&vparc,&vperpc);
		    Classe_formule::vparvperp(vpcm[index[m]],vvp,&vpart,&vperpt);
		    Classe_formule::vparvperp(vpcm[index[m]],vpcm[index[1]],&vpars,&vperps);

		 if(Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]])<90)
		   {
		 fillh("hocchialt90lt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);
		 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchialt90plt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

 fillh("hocchilt90p_zbiglt120",vparc,vperpc,1);
 fillh("hocchilt90p_ztlt120",vpart,vperpt,1);
 fillh("hocchilt90p_zslt120",vpars,vperps,1);

		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchialt90alt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

		    fillh("hocchilt90a_zbiglt120",vparc,vperpc,1);
fillh("hocchilt90a_ztlt120",vpart,vperpt,1);
fillh("hocchilt90a_zslt120",vpars,vperps,1);
		   }


		   }
		 else
		   {
fillh("hocchiagt90lt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1);

 if(TMath::Nint(evento.z[index[m]])==1)
		   {
		    fillh("hocchiagt90plt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 
fillh("hocchigt90p_zbiglt120",vparc,vperpc,1);
fillh("hocchigt90p_ztlt120",vpart,vperpt,1);
fillh("hocchigt90p_zslt120",vpars,vperps,1);
		   }
		 if(TMath::Nint(evento.z[index[m]])==2)
		   {
		    fillh("hocchiagt90alt120",vpcm[index[m]][2],sqrt(pow(vpcm[index[m]][0],2)+pow(vpcm[index[m]][1],2)),1); 

fillh("hocchigt90a_zbiglt120",vparc,vperpc,1);
fillh("hocchigt90a_ztlt120",vpart,vperpt,1);
fillh("hocchigt90a_zslt120",vpars,vperps,1);
		   }		   }
	       }

	   }

  }




			}//thetarel<120




	}//z[1]>=5
      if(evento.z[index[1]]<5)
	{
	  if(evento.rcoqf[index[0]]==1)
	    {
	  fillh("hindexb",(float)ival[index[0]],1);
	  fillh("hzvbb",Classe_formule::modulo(vpcm[index[0]]),(float)evento.z[index[0]],1);

	  if(evento.z[index[0]]>5)
	    {
fillh("hzvbb5",Classe_formule::modulo(vpcm[index[0]]),(float)evento.z[index[0]],1);


	    }

	    }


	}


     if(vpcm[index[0]][2]>0 && evento.z[index[0]]>=10 && iela==0)
	{
	  fillh("hneventi",48.,1);
 if(evento.rcoqf[index[0]]==1)
   {
fillh("hnzbiggesttuttinotf",(float)ival[index[0]],1);
fillh("hzbiggestnotf",evento.z[index[0]],1);

 int mpos=0;
 for(int j=1;j<evento.moltepl;j++)
   {
     if(evento.rcoqf[index[j]]==1)
       {
	 fillh("hnzaltrinotf",(float)ival[index[j]],1);
	 fillh("hivalbigivalaltrinotf",(float)ival[index[0]],(float)ival[index[j]],1);

	 Classe_formule::vparvperp(vpcm[index[j]],vpcm[index[0]],&vpar,&vperp);
	 fillh("hnzaltrivparnotf",vpar,(float)ival[index[j]],1);
	 fillh("hnzaltrivparmenovbignotf",vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);

	 if(evento.z[index[1]]<5)
	   {
	 fillh("hnzaltrisololcpnotf",(float)ival[index[j]],1);

	 fillh("hnzaltrivparsololcpnotf",vpar,(float)ival[index[j]],1);
	 fillh("hnzaltrivparmenovbigsololcpnotf",vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);

	 
	 	 fillh(Form("h%d_notf",100+(int)evento.z[index[0]]),vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);
	 
		 if(vpar-Classe_formule::modulo(vpcm[index[0]])>=0)mpos++;
	   }

       }

   }


 if(evento.z[index[1]]<5)
   {
     if(mpos>0)fillh("hnzbiggestconlcpposnotf",(float)ival[index[0]],1);
     fillh("hnzbiggestconlcpnotf",(float)ival[index[0]],1);




   }

   }// if(evento.rcoqf[index[0]]==1)
	}//if(vpcm[index[0]][2]>0 && evento.z[index[0]]>=10 && iela==0
    }//if(evento.moltepl>1)



  if(tflow<0)return;


   if(ztot>=14&&tflow<30)
    {
      fillh("hm",(float)evento.moltepl,1);
      if(vpcm[index[0]][2]>0 && evento.z[index[0]]>=10 && iela==0)
	{
	  fillh("hneventi",45.,1);
	  if(evento.moltepl>1)
	    {
 if(evento.rcoqf[index[0]]==1)
   {
fillh("hnzbiggesttutti",(float)ival[index[0]],1);
 fillh("hnzbiggesttuttivl",(float)ival[index[0]],evento.vpartlab[index[0]],1);
 fillh("hnzbiggesttuttivcm",(float)ival[index[0]],evento.vpartcm[index[0]],1);
 fillh("hnzbiggesttuttivlnorm",(float)ival[index[0]],evento.vpartlab[index[0]]/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzbiggesttuttivcmnorm",(float)ival[index[0]],evento.vpartcm[index[0]]/Classe_analisi::Getanalisi()->reazione.vp_cm,1);
 fillh("hzbiggest",evento.z[index[0]],1);
 int mpos=0;
 for(int j=1;j<evento.moltepl;j++)
   {
     if(evento.rcoqf[index[j]]==1)
       {
	 fillh("hnzaltri",(float)ival[index[j]],1);
	 fillh("hivalbigivalaltri",(float)ival[index[0]],(float)ival[index[j]],1);

	 Classe_formule::vparvperp(vpcm[index[j]],vpcm[index[0]],&vpar,&vperp);
	 fillh("hnzaltrivpar",vpar,(float)ival[index[j]],1);
	 fillh("hnzaltrivparmenovbig",vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);

	 if(evento.z[index[1]]<5)
	   {
	 fillh("hnzaltrisololcp",(float)ival[index[j]],1);

	 fillh("hnzaltrivparsololcp",vpar,(float)ival[index[j]],1);
	 fillh("hnzaltrivparmenovbigsololcp",vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);

	 
	 	 fillh(Form("h%d",100+(int)evento.z[index[0]]),vpar-Classe_formule::modulo(vpcm[index[0]]),(float)ival[index[j]],1);
	 
		 if(vpar-Classe_formule::modulo(vpcm[index[0]])>=0)mpos++;
	   }

       }

   }


 if(evento.z[index[1]]<5)
   {
     if(mpos>0)fillh("hnzbiggestconlcppos",(float)ival[index[0]],1);
     fillh("hnzbiggestconlcp",(float)ival[index[0]],1);
 fillh("hnzbiggestconlcpvlnorm",(float)ival[index[0]],evento.vpartlab[index[0]]/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzbiggestconlcpvcmnorm",(float)ival[index[0]],evento.vpartcm[index[0]]/Classe_analisi::Getanalisi()->reazione.vp_cm,1);



   }


   }


      
	    }
	  if(evento.moltepl==2)
	    {
	      float vcmpaircm[3];
	      float vcmpairlab[3];
	      for(int k=0;k<3;k++)
		{
		  vcmpaircm[k]=((float)evento.a[index[0]]*vpcm[index[0]][k]+(float)evento.a[index[1]]*vpcm[index[1]][k])/((float)evento.a[index[0]]+(float)evento.a[index[1]]);
		  vcmpairlab[k]=((float)evento.a[index[0]]*vplab[index[0]][k]+(float)evento.a[index[1]]*vplab[index[1]][k])/((float)evento.a[index[0]]+(float)evento.a[index[1]]);


		}

	      if(evento.rcoqf[index[0]]==1)
		{
fillh("hnzbiggest",(float)ival[index[0]],1);
 fillh("hnzbiggestvl",(float)ival[index[0]],evento.vpartlab[index[0]],1);
 fillh("hnzbiggestvcm",(float)ival[index[0]],evento.vpartcm[index[0]],1);
 fillh("hnzbiggestvlpair",(float)ival[index[0]],Classe_formule::modulo(vcmpairlab),1);
 fillh("hnzbiggestvcmpair",(float)ival[index[0]],Classe_formule::modulo(vcmpaircm),1);
 fillh("hnzbiggestvlnorm",(float)ival[index[0]],evento.vpartlab[index[0]]/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzbiggestvcmnorm",(float)ival[index[0]],evento.vpartcm[index[0]]/Classe_analisi::Getanalisi()->reazione.vp_cm,1);
 fillh("hnzbiggestvlpairnorm",(float)ival[index[0]],Classe_formule::modulo(vcmpairlab)/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzbiggestvcmpairnorm",(float)ival[index[0]],Classe_formule::modulo(vcmpaircm)/Classe_analisi::Getanalisi()->reazione.vp_cm,1);


		}
	      if(evento.rcoqf[index[1]]==1)
		{
fillh("hnzsecond",(float)ival[index[1]],1);
 fillh("hnzsecondvl",(float)ival[index[1]],evento.vpartlab[index[0]],1);
 fillh("hnzsecondvcm",(float)ival[index[1]],evento.vpartcm[index[0]],1);
 fillh("hnzsecondvlpair",(float)ival[index[0]],Classe_formule::modulo(vcmpairlab),1);
 fillh("hnzsecondvcmpair",(float)ival[index[0]],Classe_formule::modulo(vcmpaircm),1);
 fillh("hnzsecondvlnorm",(float)ival[index[1]],evento.vpartlab[index[0]]/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzsecondvcmnorm",(float)ival[index[1]],evento.vpartcm[index[0]]/Classe_analisi::Getanalisi()->reazione.vp_cm,1);
 fillh("hnzsecondvlpairnorm",(float)ival[index[0]],Classe_formule::modulo(vcmpairlab)/Classe_analisi::Getanalisi()->reazione.vplab,1);
 fillh("hnzsecondvcmpairnorm",(float)ival[index[0]],Classe_formule::modulo(vcmpaircm)/Classe_analisi::Getanalisi()->reazione.vp_cm,1);

		}
	      float v12cm[3];
	      for(int k=0;k<3;k++)
		{
		  v12cm[k]=((float)evento.a[index[0]]*vpcm[index[0]][k]+(float)evento.a[index[1]]*vpcm[index[1]][k])/((float)evento.a[index[0]]+(float)evento.a[index[1]]);
		}
	      float ivec[3],jvec[3],kvec[3];
	      Classe_formule::sdroutofplane(v12cm,ivec,jvec,kvec);
	      
	      float vvrel0[3];
	      float vvrel1[3];
	      Classe_formule::vrel(vpcm[index[0]],v12cm,vvrel0);
	      Classe_formule::vrel(vpcm[index[1]],v12cm,vvrel1);
	      float theout1;
	      Classe_formule::outofplane(vvrel1,kvec,&theout1);
	      if(TMath::Nint(evento.z[index[1]])==2 &&TMath::Nint(evento.a[index[1]])==4)
		{

		  fillh("houtalpha",theout1,1);
		  fillh("houtalphazbig",(float)evento.z[index[0]],theout1,1);
		}
	      if(evento.z[index[1]]>=5)
		{
		 fillh("houtzgt5",theout1,1); 
		 fillh("houtzgt5zbig",(float)evento.z[index[0]],theout1,1); 
		} 

		      thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);
		      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);
		      fillh("hthetarelvrel",Classe_formule::modulo(vrel),thetarel,1);
		      if(evento.z[index[1]]>=5)fillh("hthetarelvrelge5",Classe_formule::modulo(vrel),thetarel,1);

		      if(thetarel>140)
			{fillh("hm2therelgt140",(float)evento.z[index[0]],(float)evento.z[index[1]],1);
			}
		      if(thetarel<60)
			{fillh("hm2therellt60",(float)evento.z[index[0]],(float)evento.z[index[1]],1);
			}
		      if(TMath::Nint(evento.z[index[0]])+TMath::Nint(evento.z[index[1]])==Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt)
			{
			  fillh("hthetarelvrelsumzt",Classe_formule::modulo(vrel),thetarel,1);
			  if(evento.z[index[1]]>=5) fillh("hthetarelvrelsumztge5",Classe_formule::modulo(vrel),thetarel,1);  
	    }
	    }//evento.moltepl==2

      for(int j=2;j<6;j++)
	{

	  if(evento.moltepl==j)
	    {
	      for(int k=1;k<j;k++)
		{
		
		  if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[k]]==1)fillh(Form("h%dm%d",k,j),(float)ival[index[0]],(float)ival[index[k]],1);
		  fillh(Form("hz%dm%d",k,j),(float)evento.z[index[0]],(float)evento.z[index[k]],1);
		
		}
	    }
	}
	}// 1grosso


    }//periferiche		  
		 
   for(int j=0;j<evento.moltepl;j++)
     {
       if(evento.rcoqf[j]==1)
	 {
	   fillh("hival",(float)ival[j],1);
	 }
     }

   if(malpha>=3)
     {
	 float mu=4.*4./(4.+4.);
       float erel12=10000;
       int j1=-1;
       int j2=-1;
       int j3=-1;
      float vrel12[3]; 
 float vrel13[3]; 
float vrel23[3]; 

	  for(int k=0;k<malpha-1;k++)
	    {
	      for(int i=k+1;i<malpha;i++)
		{
		  for(int pp=0;pp<3;pp++)
		    {
		  vrel12[pp]=vpcm[jalpha[k]][pp]-vpcm[jalpha[i]][pp];
		    }
		  
		   float erelloc=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel12),2)/pow(Classe_formule::cluce,2);

		   if(erelloc<erel12)
		     {
		       erel12=erelloc;
		       j1=k;
		       j2=i;
		     }
		}
	    }//j1 j2 coppia di erel minima
	  //ora si fa j1 con tutte le altre che non siano j2

	  int j13=-1;
	  float erel13=1000;
	  for(int k=0;k<malpha;k++)
	    {
	      if(k!=j1&&k!=j2)
		{
		   for(int pp=0;pp<3;pp++)
		    {
		  vrel13[pp]=vpcm[jalpha[j1]][pp]-vpcm[jalpha[k]][pp];
		    }
		   float erelloc=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel13),2)/pow(Classe_formule::cluce,2);   
		   if(erelloc<erel13)
		     {
		       erel13=erelloc;
		       j13=k;

		     }

		}
	    }

	  //ora si fa j2 con tutte le altre che non siano j1

	  int j23=-1;
	  float erel23=1000;
	  for(int k=0;k<malpha;k++)
	    {
	      if(k!=j1&&k!=j2)
		{
		   for(int pp=0;pp<3;pp++)
		    {
		  vrel23[pp]=vpcm[jalpha[j2]][pp]-vpcm[jalpha[k]][pp];
		    }
		   float erelloc=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel23),2)/pow(Classe_formule::cluce,2);   
		   if(erelloc<erel23)
		     {
		       erel23=erelloc;
		       j23=k;

		     }

		}
	    }

	  //ora si confrontano erel23 e erel13; se erel23<erel13 allora come 3 si prende quello della coppia 23, altrimenti quello della coppia 13
	  j3=j13;
	  if(erel23<erel13)j3=j23;

		   for(int pp=0;pp<3;pp++)
		    {
		  vrel23[pp]=vpcm[jalpha[j2]][pp]-vpcm[jalpha[j3]][pp];
		    }	  
		   erel23=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel23),2)/pow(Classe_formule::cluce,2); 

		   for(int pp=0;pp<3;pp++)
		    {
		  vrel13[pp]=vpcm[jalpha[j1]][pp]-vpcm[jalpha[j3]][pp];
		    }	  
		   erel13=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel13),2)/pow(Classe_formule::cluce,2); 


			 float xD=sqrt(3.)*(erel23-erel12)/2;
			 float yD=(2*erel13-erel23-erel12)/2;
			 
			 fillh("hDalitztutte",xD,yD,1);
     }//malpha>=3


      if(malpha>1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=1;
	     
	    }
	  int indice=0;
	  mixatore(indice,1,2,2,4,-1,-1,-1);
	  // cout<<"malpha="<<malpha<<" "<<mixing.Prel.at(indice).size()<<endl;
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];

	  //mixing.Prel.at(indice).at(i).at(0) componente 0 di prel
	  //mixing.Prel.at(indice).at(i).at(1) componente 1 di prel
	  //mixing.Prel.at(indice).at(i).at(2) componente 2 di prel
	  //mixing.Prel.at(indice).at(i).at(3) indice in evento della prima particella della coppia su cui è calcolato prel
	  //mixing.Prel.at(indice).at(i).at(4) indice in evento della seconda particella della coppia

	  //mixing.Vcm.at(indice).at(i).at(0) componente 0 della velocità nel lab del cm della coppia
	  //mixing.Vcm.at(indice).at(i).at(1) componente 1 della velocità nel lab del cm della coppia
	  //mixing.Vcm.at(indice).at(i).at(2) componente 2 della velocità nel lab del cm della coppia
	  

	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*4/8.;
		float mtot=4.+4.;
		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		float therel=mixing.Thetarel.at(indice).at(i);
		fillh("haa_vcm_erel",vvcm,erel,1);
		fillh("haa_thetarel_erel",therel,erel,1);
		fillh("haa_erel",erel,1);
		float Q=Delta[4][8]-2*Delta[2][4];
		fillh("h8Be_estar",erel-Q,1);
			
		
		if(erel<erelmin)
		  {
		    erelmin=erel;
		    imin=i;
		  }
		
	
	}

	    vcmc[0]=mixing.Vcm.at(indice).at(imin).at(0);
	    vcmc[1]=mixing.Vcm.at(indice).at(imin).at(1);
	    vcmc[2]=mixing.Vcm.at(indice).at(imin).at(2)-Classe_analisi::Getanalisi()->reazione.vcm;

	if(erelmin<0.3&&mp>=1)
	  {
	    
	    float erelminp=1000;
	    int imp=-1;
	    float vectrel9B[3];
	    float mu=1.*8./(1.+8.);
	    for(int k=0;k<mp;k++)
	      {
		for(int i=0;i<3;i++)
		  {
		    vectrel9B[i]=vpcm[jp[k]][i]-vcmc[i];
		  }

		 float vrel9B=Classe_formule::modulo(vectrel9B);
		     float erel9B=0.5*mu*Classe_formule::amu*pow(vrel9B,2)/pow(Classe_formule::cluce,2);
		     
		     if(erelminp>erel9B)
		       {
			 erelminp=erel9B;
			 imp=jp[k];
		       }


	      }
	    if(imp>=0)
	      {
		     fillh("h8Bep_erel",erelminp,1);
		     float Q=Delta[5][9]-Delta[4][8]-Delta[1][1];
		     fillh("h9B_estar",erelminp-Q,1);

	      }

	  }//if(erelmin<0.3&&mp>=1)


	if(erelmin<0.3 && malpha>=3)
	  {

		 float erela8Bemin=1000;
		 int ja8Bemin=-1;
		     float vectrelC[3];
		     float VC12[3];
		     float mu=4.*8./(4.+8.);
		 for(int k=0;k<malpha;k++)
		   {
		     if(jalpha[k]!=mixing.Prel.at(indice).at(imin).at(3) && jalpha[k]!=mixing.Prel.at(indice).at(imin).at(4))
		       {


		     for(int i=0;i<3;i++)
		       {
			 vectrelC[i]=vpcm[jalpha[k]][i]-vcmc[i];
			 
		       }
		     float vrelC=Classe_formule::modulo(vectrelC);
		     float erelC=0.5*mu*Classe_formule::amu*pow(vrelC,2)/pow(Classe_formule::cluce,2);
		     if(erelC<erela8Bemin)
		       {
			 erela8Bemin=erelC;
			 ja8Bemin=jalpha[k];
		       }

		       }

		   }
		 if(ja8Bemin>=0)
		   {
		     fillh("hBea_erelmge3tutti",erela8Bemin,1);
		     float Q=Delta[6][12]-Delta[4][8]-Delta[2][4];
		     fillh("h12C_estarmge3tutti",erela8Bemin-Q,1);
			 float VC12[3];
		     for(int k=0;k<3;k++)
		       {
			 VC12[k]=(4.*vpcm[ja8Bemin][k]+8.*vcmc[k])/(4.+8.);

		       } 
		     float vparC,vperpC;
		     Classe_formule::vparvperp(VC12,asse,&vparC,&vperpC);
		     fillh("h799",vparC,vperpC,1);
		     if(evento.z[index[0]]>=10)
		       {
			 fillh("h700",vparC,vperpC,1);
			 fillh("he700",erela8Bemin-Q,1);
		       }
		     if(evento.z[index[0]]>=15)
		       {
			 fillh("h701",vparC,vperpC,1);
			 fillh("he701",erela8Bemin-Q,1);
		       }
		     if(evento.z[index[0]]>=10&&evento.z[index[0]]<15)
		       {
			 fillh("h702",vparC,vperpC,1);
			 fillh("he702",erela8Bemin-Q,1);
		       }
		    
		     if(evento.z[index[0]]<10)
		       {
			 fillh("h703",vparC,vperpC,1);
			 fillh("he703",erela8Bemin-Q,1);
		       }
		     if(erela8Bemin-Q<8)//Hoyle state
		       {
			 float mu=4.*4./(4.+4.);
			 float erel12=erelmin;//1-2
			 float vrel13[3];
			 float vrel23[3];
			 int j1=mixing.Prel.at(indice).at(imin).at(3);
 			 int j2=mixing.Prel.at(indice).at(imin).at(4); 
			 for(int k=0;k<3;k++)
			   {
			     vrel13[k]=vpcm[j1][k]-vpcm[ja8Bemin][k];//1-3 
			     vrel23[k]=vpcm[j2][k]-vpcm[ja8Bemin][k];//2-3
			   }
			 float erel13=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel13),2)/pow(Classe_formule::cluce,2);
			 float erel23=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel23),2)/pow(Classe_formule::cluce,2);
			 float xD=sqrt(3.)*(erel23-erel12)/2;
			 float yD=(2*erel13-erel23-erel12)/2;
			 
			 fillh("hDalitz",xD,yD,1);

		       }//if(erela8Bemin-Q<8)
		   }


	  }//erelmin<0.3 malpha>=3



	if(erelmin<10 && malpha>=3)
	  {

		 float erela8Bemin=1000;
		 int ja8Bemin=-1;
		     float vectrelC[3];
		     float VC12[3];
		     float mu=4.*8./(4.+8.);
		 for(int k=0;k<malpha;k++)
		   {
		     if(jalpha[k]!=mixing.Prel.at(indice).at(imin).at(3) && jalpha[k]!=mixing.Prel.at(indice).at(imin).at(4))
		       {


		     for(int i=0;i<3;i++)
		       {
			 vectrelC[i]=vpcm[jalpha[k]][i]-vcmc[i];
			 
		       }
		     float vrelC=Classe_formule::modulo(vectrelC);
		     float erelC=0.5*mu*Classe_formule::amu*pow(vrelC,2)/pow(Classe_formule::cluce,2);
		     if(erelC<erela8Bemin)
		       {
			 erela8Bemin=erelC;
			 ja8Bemin=jalpha[k];
		       }

		       }

		   }
		 if(ja8Bemin>=0)
		   {
		     fillh("hBea_erel10mge3tutti",erela8Bemin,1);
		     float Q=Delta[6][12]-Delta[4][8]-Delta[2][4];
		     fillh("h12C_estar10mge3tutti",erela8Bemin-Q,1);

		     
			 float mu=4.*4./(4.+4.);
			 float erel12=erelmin;//1-2
			 float vrel13[3];
			 float vrel23[3];
			 int j1=mixing.Prel.at(indice).at(imin).at(3);
 			 int j2=mixing.Prel.at(indice).at(imin).at(4); 
			 for(int k=0;k<3;k++)
			   {
			     vrel13[k]=vpcm[j1][k]-vpcm[ja8Bemin][k];//1-3 
			     vrel23[k]=vpcm[j2][k]-vpcm[ja8Bemin][k];//2-3
			   }
			 float erel13=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel13),2)/pow(Classe_formule::cluce,2);
			 float erel23=0.5*mu*Classe_formule::amu*pow(Classe_formule::modulo(vrel23),2)/pow(Classe_formule::cluce,2);
			 float xD=sqrt(3.)*(erel23-erel12)/2;
			 float yD=(2*erel13-erel23-erel12)/2;
			 
			 fillh("hDalitze10",xD,yD,1);
			 float VC12[3];
		     for(int k=0;k<3;k++)
		       {
			 VC12[k]=(4.*vpcm[ja8Bemin][k]+8.*vcmc[k])/(4.+8.);

		       } 
		     float vparC,vperpC;
		     Classe_formule::vparvperp(VC12,asse,&vparC,&vperpC);
		     fillh("h699",vparC,vperpC,1);
		     if(evento.z[index[0]]>=10)
		       {
			 fillh("h600",vparC,vperpC,1);
			 fillh("he600",erela8Bemin-Q,1);
		       }
		     if(evento.z[index[0]]>=15)
		       {
			 fillh("h601",vparC,vperpC,1);
			 fillh("he601",erela8Bemin-Q,1);
		       }
		     if(evento.z[index[0]]>=10&&evento.z[index[0]]<15)
		       {
			 fillh("h602",vparC,vperpC,1);
			 fillh("he602",erela8Bemin-Q,1);
		       }
		     
		     if(evento.z[index[0]]<10)
		       {
			 fillh("h603",vparC,vperpC,1);
			 fillh("he603",erela8Bemin-Q,1);
		       }

		   }


	  }//erelmin<10 malpha>=3



	if(evento.z[index[0]]>=10)
	  {
	    float vpar,vperp;
	    



	     Classe_formule::vparvperp(vcmc,vpcm[index[0]],&vpar,&vperp);
	     
	     fillh(Form("h%d",200+(int)evento.z[index[0]]),vpar,erelmin,1);
	     fillh("h200",vpar,erelmin,1);
	     fillh(Form("h%d",(int)ztot*100+(int)evento.z[index[0]]),vpar,erelmin,1);
	     if(erelmin<0.3)
	       {
	     fillh(Form("h%d",300+(int)evento.z[index[0]]),vpar,vperp,1);
	     fillh("h300",vpar,vperp,1);


	     if(malpha>=3)
	       {
		 float erela8Bemin=1000;
		 int ja8Bemin=-1;
		     float vectrelC[3];
		     float VC12[3];
		     float mu=4.*8./(4.+8.);
		 for(int k=0;k<malpha;k++)
		   {
		     if(jalpha[k]!=mixing.Prel.at(indice).at(imin).at(3) && jalpha[k]!=mixing.Prel.at(indice).at(imin).at(4))
		       {


		     for(int i=0;i<3;i++)
		       {
			 vectrelC[i]=vpcm[jalpha[k]][i]-vcmc[i];
			 
		       }
		     float vrelC=Classe_formule::modulo(vectrelC);
		     float erelC=0.5*mu*Classe_formule::amu*pow(vrelC,2)/pow(Classe_formule::cluce,2);
		     if(erelC<erela8Bemin)
		       {
			 erela8Bemin=erelC;
			 ja8Bemin=jalpha[k];
		       }

		       }

		   }
		 if(ja8Bemin>=0)
		   {
		     fillh("hBea_erelmge3",erela8Bemin,1);
		     float Q=Delta[6][12]-Delta[4][8]-Delta[2][4];
		     fillh("h12C_estarmge3",erela8Bemin-Q,1);
		   }


	       }//malpha>=3


	     if(malpha==3)
	       {
		 int ja=-1;
		 for(int k=0;k<malpha;k++)
		   {
		     if(jalpha[k]!=mixing.Prel.at(indice).at(imin).at(3) && jalpha[k]!=mixing.Prel.at(indice).at(imin).at(4))
		       {
			 ja=jalpha[k];
			 break;
		       }

		   }

		 if(ja>=0)
		   {
		     float vpara,vperpa;
		     Classe_formule::vparvperp(vpcm[ja],vpcm[index[0]],&vpara,&vperpa);
		     fillh("h401",vpara,vperpa,1);
		     float mu=4.*8./(4.+8.);
		     float vectrelC[3];
		     float VC12[3];
		     for(int k=0;k<3;k++)
		       {
			 vectrelC[k]=vpcm[ja][k]-vcmc[k];
			 VC12[k]=(4.*vpcm[ja][k]+8.*vcmc[k])/(4.+8.);
		       }
		     float vrelC=Classe_formule::modulo(vectrelC);
		     float erelC=0.5*mu*Classe_formule::amu*pow(vrelC,2)/pow(Classe_formule::cluce,2);
fillh("hBea_erel",erelC,1);
 float Q=Delta[6][12]-Delta[4][8]-Delta[2][4];
 fillh("h12C_estar",erelC-Q,1);
 
 float vpar2,vperp2;
 Classe_formule::vparvperp(VC12,vpcm[index[0]],&vpar2,&vperp2);

fillh("h400",vpar2,vperp2,1);

fillh("h405",vpar2,erelC-Q,1);

// cout<<vpar2<<" "<<erelC-Q<<endl;


float thetarelq=Classe_formule::thetarel(VC12,vpcm[index[0]]);
 float vrelq[3];
Classe_formule::vrel(VC12,vpcm[index[0]],vrelq);
 fillh("h402",Classe_formule::modulo(vrelq),thetarelq,1);

 if(vpar2>-8)
   {
     fillh("h403",Classe_formule::modulo(vrelq),thetarelq,1);
   }
 else
   {
     fillh("h404",Classe_formule::modulo(vrelq),thetarelq,1);
   }
 if(erelC-Q<15)
   {
       fillh("h406",Classe_formule::modulo(vrelq),thetarelq,1);
       if(Classe_analisi::Getanalisi()->reazione.ebeam>25)
	 {
	   if(vpar2>8)
	     {
fillh("h407",Classe_formule::modulo(vrelq),thetarelq,1);
	       
	     }
	   else
	     {
fillh("h408",Classe_formule::modulo(vrelq),thetarelq,1);
	     }
	 }
       else
	 {
	   if(vpar2>0)
	     {
fillh("h407",Classe_formule::modulo(vrelq),thetarelq,1);
	     }
	   else
	     {
fillh("h408",Classe_formule::modulo(vrelq),thetarelq,1);
	     }

	 }

   }

		   }





	       }//malpha==3
	       }//erelmin <0.3
	  } //eventoz[0]>10

 
	       



/* 	if(malpha>2) */
/* 	  { */
/* 	    for(int k=0;k<malpha-1;k++) */
/* 	      { */
/* 		for(int k1=k+1;k1<malpha;k1++) */
/* 		  { */
/* 		   float  erel=0.5*931.5*2*(pow((vpcm[jalpha[k]][0]- vpcm[jalpha[k1]][0]),2)+  pow((vpcm[jalpha[k]][1]- vpcm[jalpha[k1]][1]),2)+pow((vpcm[jalpha[k]][2]- vpcm[jalpha[k1]][2]),2))/(300.*300.);  */
/* 		    cout<<k<<"k "<<k1<<" "<<erel<<endl; */
/* 		  } */
		
/* 	      } */
/* 	  } */

 	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){ 
 		float mu=4*4/8.; 
		float mtot=4.+4.;
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
	
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		float therel=mixing.Thetamix.at(indice).at(i);
		fillh("haa_vcm_emix",vvcm,erel,1);
		fillh("haa_thetarel_emix",therel,erel,1);
		fillh("haa_emix",erel,1);
		float Q=Delta[4][8]-2*Delta[2][4];
		fillh("h8Be_estarmix",erel-Q,1);

 		} 

	} //malpha>1
      if(malpha>=1 && md>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=2;
	    }
	  for(int k=0;k<md;k++)
	    {
	      evento.isformix[jd[k]]=2;
	    }

	  int indice=1;
 mixatore(indice,2,1,2,4,1,1,2);
	float Q=Delta[3][6]-Delta[1][2]-Delta[2][4];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){
		float mu=2*4/6.;
		
		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
	
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		float therel=mixing.Thetarel.at(indice).at(i);
		fillh("hda_vcm_erel",vvcm,erel,1);
		fillh("hda_thetarel_erel",therel,erel,1);
		fillh("hda_erel",erel,1);
	
		fillh("h6Li_estar",erel-Q,1);

	}

	



	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*2/6.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		float therel=mixing.Thetamix.at(indice).at(i);
		fillh("hda_vcm_emix",vvcm,erel,1);
		fillh("hda_thetarel_emix",therel,erel,1);
		fillh("hda_emix",erel,1);
fillh("h6Li_estarmix",erel-Q,1);
		}



	}

      if(malpha>=1 && m6Li>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=3;
	     
	    }
	  for(int k=0;k<m6Li;k++)
	    {
	      evento.isformix[j6Li[k]]=3;
	     
	    }

	  int indice=2;
mixatore(indice,3,1,2,4,1,3,6);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[5][10]-Delta[2][4]-Delta[3][6];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*6/10.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("ha6Li_erel",erel,1);
		
		fillh("h10B_estar",erel-Q,1);
		
	}

	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*6/10.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("ha6Li_emix",erel,1);
		fillh("h10B_estarmix",erel-Q,1);

		}



	}//6Li+alpha

      if(malpha>=1 && m7Li>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=4;
	     
	    }
	  for(int k=0;k<m7Li;k++)
	    {
	      evento.isformix[j7Li[k]]=4;
	     
	    }

	  int indice=3;
	  mixatore(indice,4,1,2,4,1,3,7);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
	float Q=Delta[5][11]-Delta[2][4]-Delta[3][7];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*7/11.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("ha7Li_erel",erel,1);
	
		fillh("h11B_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*7/11.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("ha7Li_emix",erel,1);
		fillh("h11B_estarmix",erel-Q,1);

		}
	}//7Li+alpha



      if(malpha>=1 && m16O>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=5;
	     
	    }
	  for(int k=0;k<m16O;k++)
	    {
	      evento.isformix[j16O[k]]=5;
	     
	    }

	  int indice=4;
	  mixatore(indice,5,1,2,4,1,8,16);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[10][20]-Delta[2][4]-Delta[8][16];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*16/20.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h16Oa_erel",erel,1);
		
		fillh("h20Ne_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*16/20.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h16Oa_emix",erel,1);
		fillh("h20Ne_estarmix",erel-Q,1);

		}
	}//16O+alpha

      if(malpha>=1 && m12C>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=6;
	     
	    }
	  for(int k=0;k<m12C;k++)
	    {
	      evento.isformix[j12C[k]]=6;
	     
	    }

	  int indice=5;
	  mixatore(indice,6,1,2,4,1,6,12);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
	float Q=Delta[8][16]-Delta[2][4]-Delta[6][12];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*12/16.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h12Ca_erel",erel,1);
	
		fillh("h16O_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*12/16.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h12Ca_emix",erel,1);
		fillh("h16O_estarmix",erel-Q,1);

		}
	}//12C+alpha


      if(malpha>=1 && m20Ne>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=7;
	     
	    }
	  for(int k=0;k<m20Ne;k++)
	    {
	      evento.isformix[j20Ne[k]]=7;
	     
	    }

	  int indice=6;
	  mixatore(indice,7,1,2,4,1,10,20);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
	float Q=Delta[12][24]-Delta[2][4]-Delta[10][20];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*20/24.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h20Nea_erel",erel,1);
	
		fillh("h24Mg_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*20/24.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h20Nea_emix",erel,1);
		fillh("h24Mg_estarmix",erel-Q,1);

		}
	}//20Ne+alpha


      if(malpha>=1 && m24Mg>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=8;
	     
	    }
	  for(int k=0;k<m24Mg;k++)
	    {
	      evento.isformix[j24Mg[k]]=8;
	     
	    }

	  int indice=7;
	  mixatore(indice,8,1,2,4,1,12,24);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[14][28]-Delta[2][4]-Delta[12][24];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*24/28.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h24Mga_erel",erel,1);
		
		fillh("h28Si_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*24/28.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h24Mga_emix",erel,1);
		fillh("h28Si_estarmix",erel-Q,1);

		}
	}//24Mg+alpha


      if(malpha>=1 && m28Si>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=9;
	     
	    }
	  for(int k=0;k<m28Si;k++)
	    {
	      evento.isformix[j28Si[k]]=9;
	     
	    }

	  int indice=8;
	  mixatore(indice,9,1,2,4,1,14,28);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
	float Q=Delta[16][32]-Delta[2][4]-Delta[14][28];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*28/32.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h28Sia_erel",erel,1);
	
		fillh("h32S_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*28/32.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h28Sia_emix",erel,1);
		fillh("h32S_estarmix",erel-Q,1);

		}
	}//28Si+alpha


      if(malpha>=1 && m32S>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=10;
	     
	    }
	  for(int k=0;k<m32S;k++)
	    {
	      evento.isformix[j32S[k]]=10;
	     
	    }

	  int indice=9;
	  mixatore(indice,10,1,2,4,1,16,32);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[18][36]-Delta[2][4]-Delta[16][32];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*32/36.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h32Sa_erel",erel,1);
		
		fillh("h36Ar_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*32/36.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h32Sa_emix",erel,1);
		fillh("h36Ar_estarmix",erel-Q,1);

		}
	}//32S+alpha



      if(malpha>=1 && m36Ar>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=11;
	     
	    }
	  for(int k=0;k<m36Ar;k++)
	    {
	      evento.isformix[j36Ar[k]]=11;
	     
	    }

	  int indice=10;
	  mixatore(indice,11,1,2,4,1,18,36);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[20][40]-Delta[2][4]-Delta[18][36];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*36/40.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h36Ara_erel",erel,1);
		
		fillh("h40Ca_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*36/40.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h36Ara_emix",erel,1);
		fillh("h40Ca_estarmix",erel-Q,1);

		}
	}//36Ar+alpha

      if(mp>=1 && m27Al>=1)
	{
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=12;
	     
	    }
	  for(int k=0;k<m27Al;k++)
	    {
	      evento.isformix[j27Al[k]]=12;
	     
	    }

	  int indice=11;
	  mixatore(indice,12,1,1,1,1,13,27);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
	float Q=Delta[14][28]-Delta[1][1]-Delta[13][27];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=1*27/28.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h27Alp_erel",erel,1);
	
		fillh("h28Sip_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=1*27/28.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h27Alp_emix",erel,1);
		fillh("h28Sip_estarmix",erel-Q,1);

		}
	}//27Al+p







      if(malpha>=1 && m40Ca>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=13;
	     
	    }
	  for(int k=0;k<m40Ca;k++)
	    {
	      evento.isformix[j40Ca[k]]=13;
	     
	    }

	  int indice=12;
	  mixatore(indice,13,1,2,4,1,20,40);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[22][44]-Delta[2][4]-Delta[20][40];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*40/44.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h40Caa_erel",erel,1);
		
		fillh("h44Ti_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*40/44.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h40Caa_emix",erel,1);
		fillh("h44Ti_estarmix",erel-Q,1);

		}
	}//40Ca+alpha


      if(mp>=1 && m40Ca>=1)
	{
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=14;
	     
	    }
	  for(int k=0;k<m40Ca;k++)
	    {
	      evento.isformix[j40Ca[k]]=14;
	     
	    }

	  int indice=13;
	  mixatore(indice,14,1,1,1,1,20,40);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[21][41]-Delta[1][1]-Delta[20][40];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=1*40/41.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h40Cap_erel",erel,1);
		
		fillh("h41Sc_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=1*40/41.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h40Cap_emix",erel,1);
		fillh("h41Sc_estarmix",erel-Q,1);

		}
	}//40Ca+p


      if(malpha>=1 && m48Ca>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=15;
	     
	    }
	  for(int k=0;k<m48Ca;k++)
	    {
	      evento.isformix[j48Ca[k]]=15;
	     
	    }

	  int indice=14;
	  mixatore(indice,15,1,2,4,1,20,48);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[22][52]-Delta[2][4]-Delta[20][48];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*48/52.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h48Caa_erel",erel,1);
		
		fillh("h52Ti_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*48/52.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h48Caa_emix",erel,1);
		fillh("h52Ti_estarmix",erel-Q,1);

		}
	}//48Ca+alpha


      if(mp>=1 && m48Ca>=1)
	{
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=16;
	     
	    }
	  for(int k=0;k<m48Ca;k++)
	    {
	      evento.isformix[j48Ca[k]]=16;
	     
	    }

	  int indice=15;
	  mixatore(indice,16,1,1,1,1,20,48);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[21][49]-Delta[1][1]-Delta[20][48];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=1*48/49.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h48Cap_erel",erel,1);
		
		fillh("h49Sc_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=1*48/49.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h48Cap_emix",erel,1);
		fillh("h49Sc_estarmix",erel-Q,1);

		}
	}//48Ca+p

      if(malpha>=1 && m7Be>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=17;
	     
	    }
	  for(int k=0;k<m7Be;k++)
	    {
	      evento.isformix[j7Be[k]]=17;
	     
	    }

	  int indice=16;
	  mixatore(indice,17,1,2,4,1,4,7);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][11]-Delta[2][4]-Delta[4][7];
 float estarC11=0;
 int n11C=0;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*7/11.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h7Bea_erel",erel,1);
		
		fillh("h11C_estar",erel-Q,1);
		estarC11=erel-Q;
	if(estarC11<15)
	  {
	    n11C=1;

	    if(evento.z[index[0]]>=10)
	      {
		float v11C[3];
		int i3=mixing.Prel.at(indice).at(i).at(3);
		int i4=mixing.Prel.at(indice).at(i).at(3);


		for(int kk=0;kk<3;kk++)
		  {
		    v11C[kk]=((float)evento.a[i3]*vpcm[i3][kk]+(float)evento.a[i4]*vpcm[i4][kk])/((float)evento.a[i3]+(float)evento.a[i4]);
		  }

		Classe_formule::vparvperp(v11C,vpcm[index[0]],&vpar,&vperp);
		fillh("hocchi11C",vpar,vperp,1);
		float vrel11C[3];
		Classe_formule::vrel(v11C,vpcm[index[0]],vrel11C);
		float thetarel11C=Classe_formule::thetarel(v11C,vpcm[index[0]]);
		fillh("hthetarelvrel11C",Classe_formule::modulo(vrel11C),thetarel11C,1);
		fillh("hvz11C",v11C[2],1);
	      }

	  }
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*7/11.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h7Bea_emix",erel,1);
		fillh("h11C_estarmix",erel-Q,1);

		}
	if(n11C==1)
	  {
	    for(int kk=0;kk<evento.moltepl;kk++)
	      {
		
		if(TMath::Nint(evento.z[kk])!=2 && TMath::Nint(evento.a[kk])!=4 && TMath::Nint(evento.z[kk])!=4 && TMath::Nint(evento.a[kk])!=7)
		  {
		fillh("hzvz_11C",vpcm[kk][2],(float)evento.z[kk],1);
		  }

	      }
	  }





	}//7Be+a




      if(malpha>=1 && m9Be>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=18;
	     
	    }
	  for(int k=0;k<m9Be;k++)
	    {
	      evento.isformix[j9Be[k]]=18;
	     
	    }

	  int indice=17;
	  mixatore(indice,18,1,2,4,1,4,9);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][13]-Delta[2][4]-Delta[4][9];
 float estar13C=0;
 int n13C=0;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*9/13.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h9Bea_erel",erel,1);
		
		fillh("h13C_estar",erel-Q,1);
		estar13C=erel-Q;

	if(estar13C<15)
	  {
	    n13C=1;

	    if(evento.z[index[0]]>=10)
	      {
		float v13C[3];
		int i3=mixing.Prel.at(indice).at(i).at(3);
		int i4=mixing.Prel.at(indice).at(i).at(3);


		for(int kk=0;kk<3;kk++)
		  {
		    v13C[kk]=((float)evento.a[i3]*vpcm[i3][kk]+(float)evento.a[i4]*vpcm[i4][kk])/((float)evento.a[i3]+(float)evento.a[i4]);
		  }

		Classe_formule::vparvperp(v13C,vpcm[index[0]],&vpar,&vperp);
		fillh("hocchi13C",vpar,vperp,1);
		float vrel13C[3];
		Classe_formule::vrel(v13C,vpcm[index[0]],vrel13C);
		float thetarel13C=Classe_formule::thetarel(v13C,vpcm[index[0]]);
		fillh("hthetarelvrel13C",Classe_formule::modulo(vrel13C),thetarel13C,1);
		fillh("hvz13C",v13C[2],1);
	      }

	  }


		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*9/13.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h9Bea_emix",erel,1);
		fillh("h13C_estarmix",erel-Q,1);

		}
	if(n13C==1)
	  {
	    for(int kk=0;kk<evento.moltepl;kk++)
	      {
		
		if(TMath::Nint(evento.z[kk])!=2 && TMath::Nint(evento.a[kk])!=4 && TMath::Nint(evento.z[kk])!=4 && TMath::Nint(evento.a[kk])!=9)
		  {
		fillh("hzvz_13C",vpcm[kk][2],(float)evento.z[kk],1);
		  }

	      }
	  }
	}//9Be+a

      if(malpha>=1 && m10Be>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=19;
	     
	    }
	  for(int k=0;k<m10Be;k++)
	    {
	      evento.isformix[j10Be[k]]=19;
	     
	    }

	  int indice=18;
	  mixatore(indice,19,1,2,4,1,4,10);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][14]-Delta[2][4]-Delta[4][10];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*10/14.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h10Bea_erel",erel,1);
		
		fillh("h14C_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*10/14.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h10Bea_emix",erel,1);
		fillh("h14C_estarmix",erel-Q,1);

		}
	}//10Be+a



      if(malpha>=1 && m6He>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=20;
	     
	    }
	  for(int k=0;k<m6He;k++)
	    {
	      evento.isformix[j6He[k]]=20;
	     
	    }

	  int indice=19;
	  mixatore(indice,20,1,2,4,1,2,6);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[4][10]-Delta[2][4]-Delta[2][6];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*6/10.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h6Hea_erel",erel,1);
		
		fillh("h10Be_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*6/10.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h6Hea_emix",erel,1);
		fillh("h10Be_estarmix",erel-Q,1);

		}
	}//6He+a

      if(malpha>=1 && mp>=1)
	{
	  for(int k=0;k<malpha;k++)
	    {
	      evento.isformix[jalpha[k]]=21;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=21;
	     
	    }

	  int indice=20;
	  mixatore(indice,21,1,2,4,1,1,1);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[3][5]-Delta[2][4]-Delta[1][1];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=4*1/5.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("hpa_erel",erel,1);
		
		fillh("h5Li_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=4*1/5.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("hpa_emix",erel,1);
		fillh("h5Li_estarmix",erel-Q,1);

		}
	}//alpha+p










      if(m3He>=1 && md>=1)
	{
	  for(int k=0;k<m3He;k++)
	    {
	      evento.isformix[j3He[k]]=22;
	     
	    }
	  for(int k=0;k<md;k++)
	    {
	      evento.isformix[jd[k]]=22;
	     
	    }

	  int indice=21;
	  mixatore(indice,22,1,2,3,1,1,2);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[3][5]-Delta[2][3]-Delta[1][2];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=3*2/5.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("hd3He_erel",erel,1);
		
		fillh("h5Li_d3He_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=3*2/5.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("hd3He_emix",erel,1);
		fillh("h5Li_d3He_estarmix",erel-Q,1);

		}
	}//3He+d



      if(m9Be>=1 && mp>=1)
	{
	  for(int k=0;k<m9Be;k++)
	    {
	      evento.isformix[j9Be[k]]=23;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=23;
	     
	    }

	  int indice=22;
	  mixatore(indice,23,1,1,1,1,4,9);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[5][10]-Delta[1][1]-Delta[4][9];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=1*9/10.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("hp9Be_erel",erel,1);
		
		fillh("h10B_p9Be_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=1*9/10.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("hp9Be_emix",erel,1);
		fillh("h10B_p9Be_estarmix",erel-Q,1);

		}
	}//9Be+p



      if(m3He>=1 && mp>=1)
	{
	  for(int k=0;k<m3He;k++)
	    {
	      evento.isformix[j3He[k]]=24;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=24;
	     
	    }

	  int indice=23;
	  mixatore(indice,24,1,1,1,1,2,3);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[3][4]-Delta[2][3]-Delta[1][1];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=3*1/4.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("hp3He_erel",erel,1);
		
		fillh("h4Li_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=3*1/4.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("hp3He_emix",erel,1);
		fillh("h4Li_estarmix",erel-Q,1);

		}
	}//3He+p



      if(m7Li>=1 && md>=1)
	{
	  for(int k=0;k<m7Li;k++)
	    {
	      evento.isformix[j7Li[k]]=25;
	     
	    }
	  for(int k=0;k<md;k++)
	    {
	      evento.isformix[jd[k]]=25;
	     
	    }

	  int indice=24;
	  mixatore(indice,25,1,1,2,1,3,7);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[4][9]-Delta[3][7]-Delta[1][2];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=7*2/9.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("hd7Li_erel",erel,1);
		
		fillh("h9Be_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=7*2/9.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("hd7Li_emix",erel,1);
		fillh("h9Be_estarmix",erel-Q,1);

		}
	}//7Li+d


      if(m6He>=1 && mt>=1)
	{
	  for(int k=0;k<m6He;k++)
	    {
	      evento.isformix[j6He[k]]=26;
	     
	    }
	  for(int k=0;k<mt;k++)
	    {
	      evento.isformix[jt[k]]=26;
	     
	    }

	  int indice=25;
	  mixatore(indice,26,1,1,3,1,2,6);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[3][9]-Delta[2][6]-Delta[1][3];
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	
		float mu=6*3/9.;

		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h6Het_erel",erel,1);
		
		fillh("h9Li_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		float mu=6*3/9.;
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h6Het_emix",erel,1);
		fillh("h9Li_estarmix",erel-Q,1);

		}
	}//6He+t







      if(m14N>=1 && mp>=1)
	{
	  for(int k=0;k<m14N;k++)
	    {
	      evento.isformix[j14N[k]]=27;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=27;
	     
	    }

	  int indice=26;
	  mixatore(indice,27,1,1,1,1,7,14);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[8][15]-Delta[1][1]-Delta[7][14];
		float mu=14*1/15.;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	


		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h14Np_erel",erel,1);
		
		fillh("h15O_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h14Np_emix",erel,1);
		fillh("h15O_estarmix",erel-Q,1);

		}
	}//14N+p



     if(m13C>=1 && mp>=1)
	{
	  for(int k=0;k<m13C;k++)
	    {
	      evento.isformix[j13C[k]]=28;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=28;
	     
	    }

	  int indice=27;
	  mixatore(indice,28,1,1,1,1,6,13);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[7][14]-Delta[1][1]-Delta[6][13];
		float mu=13*1/14.;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	


		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h13Cp_erel",erel,1);
		
		fillh("h14N_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h13Cp_emix",erel,1);
		fillh("h14N_estarmix",erel-Q,1);

		}
	}//13C+p


     if(m12B>=1 && mp>=1)
	{
	  for(int k=0;k<m12B;k++)
	    {
	      evento.isformix[j12B[k]]=29;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=29;
	     
	    }

	  int indice=28;
	  mixatore(indice,29,1,1,1,1,5,12);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][13]-Delta[1][1]-Delta[5][12];
		float mu=12*1/13.;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	


		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h12Bp_erel",erel,1);
		
		fillh("h13Cdap_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h12Bp_emix",erel,1);
		fillh("h13Cdap_estarmix",erel-Q,1);

		}
	}//12B+p


     if(m11B>=1 && mp>=1)
	{
	  for(int k=0;k<m11B;k++)
	    {
	      evento.isformix[j11B[k]]=30;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=30;
	     
	    }

	  int indice=29;
	  mixatore(indice,30,1,1,1,1,5,11);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][12]-Delta[1][1]-Delta[5][11];
		float mu=11*1/12.;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	


		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h11Bp_erel",erel,1);
		
		fillh("h12Cdap_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h11Bp_emix",erel,1);
		fillh("h12Cdap_estarmix",erel-Q,1);

		}
	}//11B+p

     if(m10B>=1 && mp>=1)
	{
	  for(int k=0;k<m10B;k++)
	    {
	      evento.isformix[j10B[k]]=31;
	     
	    }
	  for(int k=0;k<mp;k++)
	    {
	      evento.isformix[jp[k]]=31;
	     
	    }

	  int indice=30;
	  mixatore(indice,31,1,1,1,1,5,10);
	  float erelmin=1000;
	  int imin=-1;
	  float vcmc[3];
float Q=Delta[6][11]-Delta[1][1]-Delta[5][10];
		float mu=10*1/11.;
	for(UInt_t i=0; i<mixing.Prel.at(indice).size(); i++){

	


		float vvrel=Classe_formule::modulo(mixing.Prel.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		float vvcm=Classe_formule::modulo(mixing.Vcm.at(indice).at(i));
		 
		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		fillh("h10Bp_erel",erel,1);
		
		fillh("h11Cdap_estar",erel-Q,1);
		
	}
	for(UInt_t i=0; i<mixing.Pmix.at(indice).size(); i++){
		
		
		float vvrel=Classe_formule::modulo(mixing.Pmix.at(indice).at(i))/mu/Classe_formule::amu*Classe_formule::cluce;
		
float vvcm=Classe_formule::modulo(mixing.Vcm_mix.at(indice).at(i));

		float erel=0.5*mu*Classe_formule::amu*pow(vvrel/Classe_formule::cluce,2); //- Qval (con il segno) per estar
		
		
		fillh("h10Bp_emix",erel,1);
		fillh("h11Cdap_estarmix",erel-Q,1);

		}
	}//11B+p

}






#endif
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2);
void nsuz0(TH1F *h,TGraphErrors *gnz);
void Classe_analisi::RoutineFinale()
{
  TGraphErrors *gz[100];
  for(int j=0;j<100;j++)
    {
      gz[j]=0;
    }
  for(int j=1;j<27;j++)
    {
      gz[j]=getg(Form("gz_%d",j));
    }

  nsuz(geth1("hnzbig"),getg("gnzbig"),gz);
  ordinisup(geth1("hnzbig"),getg("gnzbig"),getg("gsigmanzbig"));
  nsuz0(geth1("hnzbiggest"),getg("gnzbiggest"));
  ordinisup(geth1("hnzbiggest"),getg("gnzbiggest"),getg("gsigmanzbiggest"));
  nsuz0(geth1("hnzsecond"),getg("gnzsecond"));
  ordinisup(geth1("hnzsecond"),getg("gnzsecond"),getg("gsigmanzsecond"));
  nsuz0(geth1("hnzbiggesttutti"),getg("gnzbiggesttutti"));
  nsuz0(geth1("hnzbiggesttuttinotf"),getg("gnzbiggesttuttinotf"));
 nsuz0(geth1("hnzbiggestconlcp"),getg("gnzbiggestconlcp"));
 nsuz0(geth1("hnzbiggestconlcpnotf"),getg("gnzbiggestconlcpnotf"));
 nsuz0(geth1("hnzaltri"),getg("gnzaltri"));
 nsuz0(geth1("hnzaltrinotf"),getg("gnzaltrinotf"));
 nsuz0(geth1("hnzbiggestconlcppos"),getg("gnzbiggestconlcppos"));
 nsuz0(geth1("hnzbiggestconlcpposnotf"),getg("gnzbiggestconlcpposnotf"));


 nsuz0(geth1("hnz30"),getg("gnz30"));
       nsuz0(geth1("hnz31"),getg("gnz31"));
       nsuz0(geth1("hnz32"),getg("gnz32"));

       nsuz0(geth1("hindexb"),getg("gindexb"));
       nsuz0(geth1("hindextot"),getg("gindextot"));

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




