#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>

#define RCO(j) (((evento.coderiv[j]>=1000000)&&(evento.coderiv[j]<10000000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
#define GARF(j) (((evento.coderiv[j]>=1000)&&(evento.coderiv[j]<10000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))

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


void Classe_evento::AnalisiPrincipale()
{
//	if (evento.run<21000) {return;}
//printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);
//cout<<"molt="<<evento.moltepl<<endl;

  fillh("hneventi",1.,1);
  fillh("hmolt",evento.moltepl,1); 
  if(Classe_analisi::Getanalisi()->tipo_analisi<20)
    {
      fillh("hbfiltrato",mcevent.par_urto,1);
    }


  vector <int> ibloc(evento.moltepl);
  vector <int> iqua(evento.moltepl);
  vector <int> itel(evento.moltepl);

 for(int j=0;j<evento.moltepl;j++)
    {
      ibloc[j]=-1;
     iqua[j]=-1;
     itel[j]=-1;
     if(evento.coderiv[j]>=10000000)
       {
	 ibloc[j]=evento.coderiv[j]-10000000;
	 
	 ibloc[j]=ibloc[j]/100;
	 int ival0=evento.coderiv[j]-10000000-ibloc[j]*100;
	 iqua[j]=ival0/4+1;
	 itel[j]=ival0-(iqua[j]-1)*4+1;
	 fillh("hcode2",(float)ibloc[j]*100+iqua[j]*10+itel[j],1);

	 //cout<<ibloc[j]<<" "<<iqua[j]<<" "<<itel[j]<<" "<<evento.coderiv[j]-10000000<<endl;
	 
       }
     
    }
float vpcm[100][3];
float vplab[100][3]; 
int ival[100];
float nsuzval[100];
float sumz=0;
float ptotz=0;
int nfrag=0;
int listafrag[100];

 int index[500];
 int vec[500];

 int ithe;

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

     if(evento.z[j]>=3)
       {
	 listafrag[nfrag]=j;
	 nfrag++;
       }
     fillh("hzv",evento.vpartlab[j],evento.z[j],1);
     ithe=evento.thetalab[j]/2;
     if(ithe<=10&&ithe>=0)
       { 
	 fillh(Form("hzvthe%d",ithe),evento.vpartlab[j],evento.z[j],1);}
     if(evento.moltepl>1)
       {
	 fillh("hzvm",evento.vpartlab[j],evento.z[j],1);
	 fillh("hzvcm",evento.vpartcm[j],evento.z[j],1);
       }

 sumz=sumz+evento.z[j];
 ptotz=ptotz+Classe_formule::amu*evento.a[j]*evento.vpartlab_z[j]/Classe_formule::cluce;
   }//giro sulla moltepl
 
 TMath::Sort(evento.moltepl,vec,index);

 if(evento.moltepl>=2)
   {
 if(TMath::Nint(evento.z[index[0]])==TMath::Nint(evento.z[index[1]])) //se il piu' grosso e il secondo piu' grosso sono uguali, si prende come piu' grosso quello a theta piu' piccolo
   {
     fillh("hneventi",31.,1);
     if(evento.thetacm[index[0]]>evento.thetacm[index[1]])
       {
	 int serv=index[0];
	 index[0]=index[1];
	 index[1]=serv;
       }
   }
   }
 if(evento.moltepl>=2)
   {
 fillh("hztotthebiggest",evento.thetacm[index[0]],sumz,1);
   }
 ptotz=ptotz/(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
 fillh("hztotptotz",ptotz,(float)sumz,1);


fillh("hmfrag",(float)nfrag,1);


 

 //calcolo thetaflow
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
     fillh("htflowztot",tflow,sumz,1);
   }

 if(evento.moltepl>1)
   {
fillh("hmfragmgt2",(float)nfrag,1);
   }

 if(nfrag==0) //se non ci sono frammenti con Z>=3 si esce
   {
     fillh("hneventi",10.,1);
     return;
   }
fillh("hneventi",11.,1);

 if(evento.moltepl<2)
   {
     fillh("hneventi",2.,1);
     return;
   }
if(tflow>0)
   {
     fillh("htflowztotgood",tflow,sumz,1);
     fillh("htflowtbiggest",tflow,evento.thetacm[index[0]],1);
//Gius 1/9/17
     fillh("htflowtbiggestmolt",evento.moltepl, tflow-evento.thetacm[index[0]],1);
     fillh("htflowtheorico",tflow,mcevent.par_urto,1);
     if(mcevent.par_urto<3)(fillh("htflowtheorico1",tflow,mcevent.par_urto,1));
     if(mcevent.par_urto>3&&mcevent.par_urto<8)(fillh("htflowtheorico2",tflow,mcevent.par_urto,1));
     if(mcevent.par_urto>8&&mcevent.par_urto<11)(fillh("htflowtheorico3",tflow,mcevent.par_urto,1));
   }


fillh("hneventi",3.,1);

fillh("hthecmvcmtot",evento.vpartcm[index[0]],evento.thetacm[index[0]],1);

//Giu 5/6/17
fillh("hthevlabvlabtot",evento.vpartlab[index[0]],evento.thetalab[index[0]],1);
if(evento.rcoqf[index[0]]==1 && evento.a[0]>14 && evento.a[0]<40) fillh("hthecmtot",evento.thetacm[index[0]],1);
//finish

 for(int j=0;j<evento.moltepl;j++)
   {
     if(evento.rcoqf[j]==1)
       {
	 fillh("hnztutti",(float)ival[j],1);
       }
   }
     float thetarel=-1;
     float vrel[3];
     float thetarel2=-1;
     float vrel2[3];
     int idic=-1;
     int ifis=-1;
     int idic2=-1;
     int ifis2=-1;
     int idic218=-1;
     int ifis218=-1;

     int idic1=-1;
     int ifis1=-1;
     int idic18=-1;
     int ifis118=-1;
     int idic118=-1;
     int ifis18=-1;
     int iqp=-1;
     int iqp18=-1;
     int icent=-1;
     int icent18=-1;

     float vpar=-1000;
     float vperp=-1000;
     float vcmcoppiacm[3];
     float vcmcoppialab[3];
     float zcoppia=0;

if(tflow>=50 && tflow<=90 && sumz>=12)
   {
     icent=1;
fillh("hneventi",21.,1);
fillh("hmoltgrossone",evento.moltepl,1);
fillh("hfraggrossone",nfrag,1);

 fillh("hthecmvcmicent",evento.vpartcm[index[0]],evento.thetacm[index[0]],1);

  if(Classe_analisi::Getanalisi()->tipo_analisi<20)
    {
      fillh("hbicent",mcevent.par_urto,1);
    }

     fillh("hzvgrossone",evento.vpartlab[index[0]],evento.z[index[0]],1);
 thetarel=-1;
	  if(nfrag>=2)
	    {
	     
	      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);
	      thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);  
	      fillh("hvreltherelselgrossone",Classe_formule::modulo(vrel),thetarel,1);
	      if(evento.z[index[1]]>=5)
		   {
		     fillh("hvreltherelselgrossonegt5",Classe_formule::modulo(vrel),thetarel,1);
		   }
	    }
     if(evento.rcoqf[index[0]]==1)
	   {
	     fillh("hzngrossone",(float)ival[index[0]],1);
	     fillh("hmassagrossone",evento.a[index[0]],1);
	   }
     for(int j=0;j<evento.moltepl;j++)
       {
	 fillh("hzvgrossonet",evento.vpartlab[j],evento.z[j],1);
	 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4)
	   {
	 Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	 fillh("hocchiagrossone",vpar,vperp,1);
	   }
       }

   }//tflow>=50 && tflow<=90 && sumz>=12

 if(tflow>=50 && tflow<=90 && sumz>18)
   {
     icent18=1;
fillh("hneventi",23.,1);
     fillh("hzvgrossone18",evento.vpartlab[index[0]],evento.z[index[0]],1);
 thetarel=-1;
	  if(nfrag>=2)
	    {
	     
	      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);
	      thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);  
	      fillh("hvreltherelselgrossone18",Classe_formule::modulo(vrel),thetarel,1);
	      if(evento.z[index[1]]>=5)
		   {
		     fillh("hvreltherelselgrossonegt518",Classe_formule::modulo(vrel),thetarel,1);
		   }
	    }
     if(evento.rcoqf[index[0]]==1)
	   {
	     fillh("hzngrossone18",(float)ival[index[0]],1);
	   }
     for(int j=0;j<evento.moltepl;j++)
       {
	 fillh("hzvgrossonet18",evento.vpartlab[j],evento.z[j],1);
       }

   }//tflow>=50 && tflow<=90 && sumz>18



if(tflow>=8 && tflow<=30 && sumz>=12)
   {
     fillh("hzvsel1big",evento.vpartlab[index[0]],evento.z[index[0]],1);
	 if(evento.rcoqf[index[0]]==1 &&evento.vpartcm_z[index[0]]>=0)
	   {
	     fillh("hznbig",(float)ival[index[0]],1);
	   }
       for(int j=0;j<evento.moltepl;j++)
       {
	 fillh("hzvsel1",evento.vpartlab[j],evento.z[j],1);
       }
     

     if(((nfrag==1 &&evento.z[index[0]]>=12) ||(nfrag>=2 && evento.z[index[1]]<5 &&evento.z[index[0]]>=12))&&evento.vpartcm_z[index[0]]>0)
       {
	 iqp=1;
       }
     if(((nfrag==1 &&evento.z[index[0]]>18) ||(nfrag>=2 && evento.z[index[1]]<5&&evento.z[index[0]]>18))&&evento.vpartcm_z[index[0]]>0)
//     if(((nfrag==1 &&evento.z[index[0]]>23) ||(nfrag>=2 && evento.z[index[1]]<5 &&evento.z[index[0]]>23))&&evento.vpartcm_z[index[0]]>0)
       {
	 iqp18=1;
       }
     thetarel=-1;
     thetarel2=-1;

	  if(nfrag>=2)
	    {
	     
	      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);
	      thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);  
	      fillh("hvreltherelsel1",Classe_formule::modulo(vrel),thetarel,1);
	      if(nfrag>=3 &&TMath::Nint(evento.z[index[2]])==TMath::Nint(evento.z[index[1]])&&TMath::Nint(evento.z[index[1]])>=5)
		{
	      Classe_formule::vrel(vpcm[index[2]],vpcm[index[0]],vrel2);
	      thetarel2=Classe_formule::thetarel(vpcm[index[2]],vpcm[index[0]]);
	      fillh("hneventi",30.,1);
		}//if(nfrag>=3 &&TMath::Nint(evento.z[index[2]])==TMath::Nint(evento.z[index[1]])&&TMath::Nint(evento.z[index[1]])>=5)
    if(evento.vpartcm_z[index[0]]>=0)
      {
	fillh("hzvsel1bigvzpos",evento.vpartlab[index[0]],evento.z[index[0]],1);
	 for(int j=0;j<evento.moltepl;j++)
	   {
	     if(j!=index[0])
	       {
	     fillh("hzvsel1altri",evento.vpartlab[j],evento.z[j],1);
	       }
	   }
	if(thetarel>=0&& evento.z[index[1]]>=5)
	  {
	    fillh("hvreltherelsel1vposgt5",Classe_formule::modulo(vrel),thetarel,1); 
	    fillh("htherelvrel01",Classe_formule::modulo(vrel),thetarel,1);
	    fillh("hz12therelsel1",evento.z[index[0]]+evento.z[index[1]],thetarel,1);

	    if(thetarel>=160 &&evento.z[index[0]]>=12)
	      {
		idic1=1;
	      }
	    if(thetarel>=160 &&evento.z[index[0]]>=18)
	      {
		idic118=1;
	      }
	    if(thetarel>=40 && thetarel<=100)
	      {
		for(int k=0;k<3;k++)
		   {

		   vcmcoppiacm[k]=(evento.a[index[0]]*vpcm[index[0]][k]+evento.a[index[1]]*vpcm[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
		   vcmcoppialab[k]=(evento.a[index[0]]*vplab[index[0]][k]+evento.a[index[1]]*vplab[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
			    }
		zcoppia=evento.z[index[0]]+evento.z[index[1]];
if((evento.z[index[0]]+evento.z[index[1]])>=12 &&vcmcoppiacm[2]>=0)
  {
    ifis1=1;
  }
 if((evento.z[index[0]]+evento.z[index[1]])>18&&vcmcoppiacm[2]>=0)
  {
    ifis118=1;
  }
	      }
	  }//if(thetarel>=0&& evento.z[index[1]]>=5)
	if(thetarel2>=0&& evento.z[index[2]]>=5)
	  {
	    fillh("htherelvrel02",Classe_formule::modulo(vrel2),thetarel2,1);
	    if(thetarel2>=160 &&evento.z[index[0]]>=12)
	      {
		idic2=1;
	      }
	    if(thetarel2>=160 &&evento.z[index[0]]>=18)
	      {
		idic218=1;
	      }

	    if(thetarel2>=40 && thetarel2<=100)
	      {
		for(int k=0;k<3;k++)
		   {

		   vcmcoppiacm[k]=(evento.a[index[0]]*vpcm[index[0]][k]+evento.a[index[2]]*vpcm[index[2]][k])/(evento.a[index[0]]+evento.a[index[2]]);
		   vcmcoppialab[k]=(evento.a[index[0]]*vplab[index[0]][k]+evento.a[index[2]]*vplab[index[2]][k])/(evento.a[index[0]]+evento.a[index[2]]);
			    }
		zcoppia=evento.z[index[0]]+evento.z[index[2]];
if((evento.z[index[0]]+evento.z[index[2]])>=12 &&vcmcoppiacm[2]>=0)
  {
    ifis2=1;
  }
if((evento.z[index[0]]+evento.z[index[2]])>=18 &&vcmcoppiacm[2]>=0)
  {
    ifis218=1;
  }

	      }


	  }//if(thetarel2>=0&& evento.z[index[2]]>=5)



      }//if(evento.vpartcm_z[index[0]]>=0)


	    }//if(nfrag>=2)
   }//if(tflow>=8 && tflow<=30 && sumz>=12)

 if(ifis1==1 && ifis2<0 && idic2<0)
   {
     ifis=1;
   } 
 if(ifis118==1 && ifis218<0 && idic218<0)
   {
     ifis18=1;
   } 
 if(ifis1<0 && idic1<0 && ifis2==1)
   {
     ifis=1;
     int serv=index[2];
     index[2]=index[1];
     index[1]=serv;
   } 
 if(ifis118<0 && idic118<0 && ifis218==1)
   {
     ifis18=1;
     int serv=index[2];
     index[2]=index[1];
     index[1]=serv;
   } 
 if(idic1==1 && ifis2<0 && idic2<0)
   {
     idic=1;
   } 
 if(idic118==1 && ifis218<0 && idic218<0)
   {
     idic18=1;
   } 
 if(ifis1<0 && idic1<0 && idic2==1)
   {
     idic=1;
     int serv=index[2];
     index[2]=index[1];
     index[1]=serv;
   } 
 if(ifis118<0 && idic118<0 && idic218==1)
   {
     idic18=1;
     int serv=index[2];
     index[2]=index[1];
     index[1]=serv;
   } 




  for(int j=0;j<evento.moltepl;j++) 
   { 
      if(evento.rcocode[j]==11&&evento.rcoqf[j]==1) 
        { 
 fillh(Form("hab%dq%dt%d",ibloc[j],iqua[j],itel[j]),faziakali.esi1[evento.indice_originale[j]],evento.z[j],1); 
        } 
      if(evento.rcocode[j]==11) 
        { 
 fillh(Form("hzb%dq%dt%d",ibloc[j],iqua[j],itel[j]),faziakali.esi1[evento.indice_originale[j]],evento.z[j],1); 

        } 
    } 



 if(iqp==1)
   {
 fillh("hneventi",20.,1);
	 fillh("hzvqpsolo",evento.vpartlab[index[0]],evento.z[index[0]],1);
	 if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzniqp",(float)ival[index[0]],1);
			  }
	 for(int j=0;j<evento.moltepl;j++)
	   {
	     if(j!=index[0])
	       {
	     fillh("hzvqpsoloaltri",evento.vpartlab[j],evento.z[j],1);
	       }
	   }
       }//iqp==1
 if(iqp18==1)
       {
fillh("hneventi",22.,1);

if(evento.moltepl>1){
		if(TMath::Nint(evento.z[index[1]])>=3)
		{
		float vvrel1[3];
		Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vvrel1);
		
fillh("hvvreliqpztutti",Classe_formule::modulo(vvrel1),1);
		if(TMath::Nint(evento.z[index[1]])==3)fillh("hvvreliqpz3",Classe_formule::modulo(vvrel1),1);
if(TMath::Nint(evento.z[index[1]])==4)fillh("hvvreliqpz4",Classe_formule::modulo(vvrel1),1);
}
}

	 fillh("hzvqpsolo18",evento.vpartlab[index[0]],evento.z[index[0]],1);
	 if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzniqp18",(float)ival[index[0]],1);
			  }
	 for(int j=0;j<evento.moltepl;j++)
	   {
	     if(j!=index[0])
	       {
	     fillh("hzvqpsoloaltri18",evento.vpartlab[j],evento.z[j],1);
	       }
	   }
       }//iqp18


 if(idic==1)
   {
     fillh("hneventi",18.,1);
     fillh("hmfragdic",(float)nfrag,1);
 fillh("hthecmvcmidic",evento.vpartcm[index[0]],evento.thetacm[index[0]],1);

    float mi=evento.a[index[0]]*evento.a[index[1]]/(evento.a[index[0]]+evento.a[index[1]]);
    float tke=0.5*mi*Classe_formule::amu*pow(Classe_formule::modulo(vrel),2)/pow(Classe_formule::cluce,2);
fillh("hzvseldic",evento.vpartlab[index[0]],evento.z[index[0]],1);
			if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzndic",(float)ival[index[0]],1);
			  }  
			for(int j=0;j<evento.moltepl;j++)
			  {
			    if(j!=index[0])
			      {
			    fillh("hzvseldicaltri",evento.vpartlab[j],evento.z[j],1);
			      }
			    fillh("hzvseldictutti",evento.vpartlab[j],evento.z[j],1);
			  }
			fillh("hzvseldicqpqt",evento.vpartlab[index[0]],evento.z[index[0]],1);
			fillh("hzvseldicqpqt",evento.vpartlab[index[1]],evento.z[index[1]],1);
			
//Giu 5/6/17
			fillh("hzvseldicqpqtneck",evento.vpartlab[index[0]],evento.z[index[0]],1);
			fillh("hzvseldicqpqtneck",evento.vpartlab[index[1]],evento.z[index[1]],1);
			if(evento.moltepl>2)fillh("hzvseldicqpqtneck",evento.vpartlab[index[2]],evento.z[index[2]],1);		
			if(evento.moltepl>3)fillh("hzvseldicqpqtneck",evento.vpartlab[index[3]],evento.z[index[3]],1);
//finish

			fillh("hdiff",evento.z[index[0]],tke,1);
			fillh("hwilc",evento.thetacm[index[0]],tke,1);
			fillh("hz1z2dic",evento.z[index[0]],evento.z[index[1]],1);
    for(int j=0;j<evento.moltepl;j++)
       {
	 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchiadic",vpar,vperp,1);
	   }
	 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==1)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchipdic",vpar,vperp,1);
	   }
	 if(j!=index[0]&&j!=index[1])
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hvparzidic",evento.z[j],vpar,1);
	   }
       }

     if(nfrag>=3)
       {
	 Classe_formule::vparvperp(vpcm[index[2]],vpcm[index[0]],&vpar,&vperp);
	 if(evento.z[index[2]]<7)
	   {
	 fillh(Form("hocchiz%ddic",(int)evento.z[index[2]]),vpar,vperp,1);
	   }
	 if(evento.rcoqf[index[2]]==1)
	   {
	     fillh("hnzimfdic",(float)ival[index[2]],1);
	  
	     fillh("hznvimfdic",(float)ival[index[2]],vpar,1);
	   }
       }

   }//idic==1

 if(idic18==1)
   {	
     fillh("hneventi",17.,1);
     fillh("hmfragdic18",(float)nfrag,1);
   float mi=evento.a[index[0]]*evento.a[index[1]]/(evento.a[index[0]]+evento.a[index[1]]);
    float tke=0.5*mi*Classe_formule::amu*pow(Classe_formule::modulo(vrel),2)/pow(Classe_formule::cluce,2);
		fillh("hzvseldic18",evento.vpartlab[index[0]],evento.z[index[0]],1);
		if(evento.moltepl>2)
		{
		if(TMath::Nint(evento.z[index[2]])>=3)
		{
		float vvrel1[3],vvrel2[3];
		Classe_formule::vrel(vpcm[index[2]],vpcm[index[0]],vvrel1);
		Classe_formule::vrel(vpcm[index[2]],vpcm[index[1]],vvrel2);
fillh("hvvreldicztutti",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
		if(TMath::Nint(evento.z[index[2]])==3)fillh("hvvrelzdic3",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==4)fillh("hvvrelzdic4",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==5)fillh("hvvrelzdic5",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==6)fillh("hvvrelzdic6",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==7)fillh("hvvrelzdic7",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==8)fillh("hvvrelzdic8",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==9)fillh("hvvrelzdic9",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==10)fillh("hvvrelzdic10",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
		
		}
		}
			if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzndic18",(float)ival[index[0]],1);
			  }  
			if(evento.rcoqf[index[1]]==1)
			  {
			fillh("hzndic18qt",(float)ival[index[1]],1);
			  }  
			for(int j=0;j<evento.moltepl;j++)
			  {
			    if(j!=index[0])
			      {
			    fillh("hzvseldicaltri18",evento.vpartlab[j],evento.z[j],1);
			      }
			    fillh("hzvseldictutti18",evento.vpartlab[j],evento.z[j],1);
			  }
			fillh("hzvseldicqpqt18",evento.vpartlab[index[0]],evento.z[index[0]],1);
			fillh("hzvseldicqpqt18",evento.vpartlab[index[1]],evento.z[index[1]],1);

			fillh("hdiff18",evento.z[index[0]],tke,1);
			fillh("hwilc18",evento.thetacm[index[0]],tke,1);
			fillh("hz1z2dic18",evento.z[index[0]],evento.z[index[1]],1);
     for(int j=0;j<evento.moltepl;j++)
       {
	 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchiadic18",vpar,vperp,1);
	   }
	 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==1)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchipdic18",vpar,vperp,1);
	   }

       }

     if(nfrag>=3)
       {
	 Classe_formule::vparvperp(vpcm[index[2]],vpcm[index[0]],&vpar,&vperp);
	 if(evento.z[index[2]]<7)
	   {
	 fillh(Form("hocchiz%ddic18",(int)evento.z[index[2]]),vpar,vperp,1);
	   }
	 if(evento.rcoqf[index[2]]==1)
	   {
	     fillh("hnzimfdic18",(float)ival[index[2]],1);
	  
	     fillh("hznvimfdic18",(float)ival[index[2]],vpar,1);
	   }
       }

   }//idic18==1
 float zsomma;
 float eta,alpha;
 if(ifis==1 || ifis18==1)
   {

     zsomma=evento.z[index[0]]+evento.z[index[1]];
    for(int k=0;k<3;k++)
     {
			    vcmcoppialab[k]=(evento.a[index[0]]*vplab[index[0]][k]+evento.a[index[1]]*vplab[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
			    vcmcoppiacm[k]=(evento.a[index[0]]*vpcm[index[0]][k]+evento.a[index[1]]*vpcm[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
			    
   }

     eta=(evento.z[index[0]]-evento.z[index[1]])/(evento.z[index[0]]+evento.z[index[1]]);
     alpha=Alpha(vpcm[index[0]],vpcm[index[1]],evento.a[index[0]],evento.a[index[1]]);
   }

 if(ifis==1)
   {

  if(Classe_analisi::Getanalisi()->tipo_analisi<20)
    {
      fillh("hbifis",mcevent.par_urto,1);
    }


 fillh("hneventi",19.,1);
 fillh("hthecmvcmifis",evento.vpartcm[index[0]],evento.thetacm[index[0]],1);

     if(nfrag>2)
       {
     fillh("hzv3",Classe_formule::modulo(vcmcoppialab),zsomma,1);
     fillh("hzv3",evento.vpartlab[index[2]],evento.z[index[2]],1);
       }



fillh("hzvfiscoppiavzpos",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);

     fillh("hzalphaifiss_small",alpha,evento.z[index[1]],1);
     fillh("hzalphaifiss_big",alpha,evento.z[index[0]],1);
     fillh("hzalphaifiss_eta",alpha,eta,1);

    if(evento.rcoqf[index[0]]==1)
      {
	fillh("hnzf1",(float)ival[index[0]],1);
	fillh("hnzf1alpha",(float)ival[index[0]],alpha,1);
	if(alpha<90)
	  {
	    fillh("hnzf1lt90",(float)ival[index[0]],1);
	    if(zsomma>=25 && zsomma<=35)
	      {
	    if(evento.z[index[1]]<5)
	      {
		fillh("hnzf1lt90b1",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==5)
	      {
		fillh("hnzf1lt90f5",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==6)
	      {
		fillh("hnzf1lt90f6",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==7)
	      {
		fillh("hnzf1lt90f7",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==8)
	      {
		fillh("hnzf1lt90f8",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==9)
	      {
		fillh("hnzf1lt90f9",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==10)
	      {
		fillh("hnzf1lt90f10",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==11)
	      {
		fillh("hnzf1lt90f11",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==12)
	      {
		fillh("hnzf1lt90f12",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==13)
	      {
		fillh("hnzf1lt90f13",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==14)
	      {
		fillh("hnzf1lt90f14",(float)ival[index[0]],1);	
	      }

	      }

	  }
      }
    if(evento.rcoqf[index[1]]==1)
      {
	fillh("hnzf2",(float)ival[index[1]],1);
	fillh("hnzf2alpha",(float)ival[index[1]],alpha,1);
	if(alpha<90)
	  {
	    fillh("hnzf2lt90",(float)ival[index[1]],1);
	  }
      }

	if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1){
			    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);

	fillh("hznfis",(float)ivalc,1);
			    fillh("hdeltaealfis",evento.z[index[0]]+evento.z[index[1]],evento.a[index[0]]+evento.a[index[1]]-Classe_formule::EAL(evento.z[index[0]]+evento.z[index[1]]),1);
			    fillh("hcorrdeltaeal",evento.a[index[0]]-Classe_formule::EAL(evento.z[index[0]]),evento.a[index[1]]-Classe_formule::EAL(evento.z[index[1]]),1);

			    float deltaeal=evento.a[index[0]]-Classe_formule::EAL(evento.z[index[0]])+evento.a[index[1]]-Classe_formule::EAL(evento.z[index[1]]);
			    float deltaeal2=evento.a[index[0]]-Classe_formule::EAL(evento.z[index[0]])-evento.a[index[1]]+Classe_formule::EAL(evento.z[index[1]]);
			    
			    fillh("hdeltaealeta",eta,deltaeal,1);
			    fillh("hdeltaealztot",evento.z[index[0]]+evento.z[index[1]],deltaeal,1);
			    fillh("hdeltaeal2eta",eta,deltaeal2,1);
			    fillh("hdeltaeal2ztot",evento.z[index[0]]+evento.z[index[1]],deltaeal2,1);

			  }
			
			fillh("hz1z2fis",evento.z[index[0]],evento.z[index[1]],1);
			    for(int j=0;j<evento.moltepl;j++)
			      {
				if(j!=index[0]&&j!=index[1])
				  {
				    if(TMath::Nint(evento.z[j])==2&&TMath::Nint(evento.a[j])==4)
			      {
			    Classe_formule::vparvperp(vpcm[j],vcmcoppiacm,&vpar,&vperp);
			    fillh("hocchiafis",vpar,vperp,1);
			      }
				  }
				if(j!=index[0])
				  {
				    fillh("hzalphaifiss_tutti",Alpha(vpcm[index[0]],vpcm[j],evento.a[index[0]],evento.a[j]),evento.z[j],1);
				  }

			      }






   

   }//ifis==1

 if(ifis18==1)
   {
fillh("hneventi",24.,1);

//Giu 19/10/17
for(int jj=0;jj<evento.moltepl;jj++){
	if(evento.rcoqf[index[jj]]==1){
		fillh("hnfistuttitot",(float)ival[index[jj]],1);
	}
	}
//
		


//Giu
fillh("hmolt3frag",(float)nfrag,1);
	if(evento.moltepl>2&&evento.rcoqf[index[2]]==1){
	     	Classe_formule::vparvperp(vpcm[index[2]],vpcm[index[0]],&vpar,&vperp);

		fillh("hnfissuz3frag",(float)ival[index[2]],1);
		fillh("hnfissuz3fragv",(float)ival[index[2]],vpar,1);
	}
//
	if(evento.moltepl>2){
		if(TMath::Nint(evento.z[index[2]])>=3)
		{
		float vvrel1[3];
		float vvrel2[3];
		Classe_formule::vrel(vpcm[index[2]],vpcm[index[0]],vvrel1);
		Classe_formule::vrel(vpcm[index[2]],vpcm[index[1]],vvrel2);
fillh("hvvrelifisztutti",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
		if(TMath::Nint(evento.z[index[2]])==3)fillh("hvvrelifisz3",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==4)fillh("hvvrelifisz4",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==5)fillh("hvvrelifisz5",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==6)fillh("hvvrelifisz6",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);
if(TMath::Nint(evento.z[index[2]])==7)fillh("hvvrelifisz7",Classe_formule::modulo(vvrel1),Classe_formule::modulo(vvrel2),1);

}
}

     fillh("hzalphaifiss18_small",alpha,evento.z[index[1]],1);
     fillh("hzalphaifiss18_big",alpha,evento.z[index[0]],1);
     //prova
//     if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1)
     fillh("hzalphaifiss18_eta",TMath::Cos(alpha/57.296),eta,1);

//Giu 3/10/17
	fillh("hz1vsz1z2",evento.z[index[0]]+evento.z[index[1]],evento.z[index[0]],1);
	fillh("hz2vsz1z2",evento.z[index[0]]+evento.z[index[1]],evento.z[index[1]],1);
	fillh("hz1z2vsvlab18",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);
	fillh("hz1z2vsvcm18",Classe_formule::modulo(vcmcoppiacm),evento.z[index[0]]+evento.z[index[1]],1);
	fillh("hetavsvcm",Classe_formule::modulo(vcmcoppiacm),eta,1);
	fillh("hetavsvlab",Classe_formule::modulo(vcmcoppialab),eta,1);
	fillh("hzvf1",evento.vpartlab[index[0]],evento.z[index[0]],1);
	fillh("hzvf2",evento.vpartlab[index[1]],evento.z[index[1]],1);

//

    if(evento.rcoqf[index[0]]==1)
      {
	fillh("hnz18f1",(float)ival[index[0]],1);
	fillh("hnz18f1alpha",(float)ival[index[0]],alpha,1);
//Giu 5/10/17
	if(alpha<60)fillh("hnz18f1alphaleq60",(float)ival[index[0]],1);
	if(alpha<100)fillh("hnz18f1alphageq100",(float)ival[index[0]],1);
	if(eta<0.25)fillh("hnz18f1etaleq03",(float)ival[index[0]],1);
	if(eta>0.45)fillh("hnz18f1etageq04",(float)ival[index[0]],1);	
//


//8/10/17
	if(TMath::Nint(evento.z[index[1]])==5){
		if(eta>0.5){
			fillh("hnz18fbig5eta06",(float)ival[index[0]],alpha,1);
			if(evento.rcoqf[index[1]]==1)
				fillh("hnz18fsmall5eta6",(float)ival[index[1]],alpha,1);
		}
		if(eta<0.3){
			fillh("hnz18fbig5eta02",(float)ival[index[0]],alpha,1);
			if(evento.rcoqf[index[1]]==1)
				fillh("hnz18fsmall5eta2",(float)ival[index[1]],alpha,1);
		}
	}
	if(TMath::Nint(evento.z[index[1]])==7){
		if(eta>0.5){
			fillh("hnz18fbig6eta06",(float)ival[index[0]],alpha,1);
			if(evento.rcoqf[index[1]]==1)
				fillh("hnz18fsmall6eta6",(float)ival[index[1]],alpha,1);
		}
		if(eta<0.3){
			fillh("hnz18fbig6eta02",(float)ival[index[0]],alpha,1);
		}
			if(evento.rcoqf[index[1]]==1)
				fillh("hnz18fsmall6eta2",(float)ival[index[1]],alpha,1);
	}
//
//20/10/17
	if((TMath::Nint(evento.z[index[1]])==5&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==17&&evento.rcoqf[index[0]]==1) || (TMath::Nint(evento.z[index[1]])==5&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==18&&evento.rcoqf[index[0]]==1)){
				    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
		fillh("hnz18asym_sm1",(float)ival[index[1]],alpha,1);
		fillh("hnz18asym_bg1",(float)ival[index[0]],alpha,1);
		fillh("hnz18asym_cp",(float)ivalc,1);
		//cout<<"5-17,18   "<<Classe_analisi::Getanalisi()->nentry<<endl;
}		
if((TMath::Nint(evento.z[index[1]])==5&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==20&&evento.rcoqf[index[0]]==1) || (TMath::Nint(evento.z[index[1]])==5&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==21&&evento.rcoqf[index[0]]==1)){
				    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
		fillh("hnz18asym_sm2",(float)ival[index[1]],alpha,1);
		fillh("hnz18asym_bg2",(float)ival[index[0]],alpha,1);
		fillh("hnz18asym_cp",(float)ivalc,1);
		//	cout<<"5-20,21   "<<Classe_analisi::Getanalisi()->nentry<<endl;
}		

if((TMath::Nint(evento.z[index[1]])==10&&evento.rcoqf[index[1]]==1&&((TMath::Nint(evento.z[index[0]]))>=15&&(TMath::Nint(evento.z[index[0]]))<=18)&&evento.rcoqf[index[0]]==1) || (TMath::Nint(evento.z[index[1]])==11&&evento.rcoqf[index[1]]==1&&((TMath::Nint(evento.z[index[0]]))>=15&&(TMath::Nint(evento.z[index[0]]))<=18)&&evento.rcoqf[index[0]]==1)){
				    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
		fillh("hnz18sym_cp",(float)ivalc,1);
		fillh("hnz18sym_sm",(float)ival[index[1]],alpha,1);
		fillh("hnz18sym_bg",(float)ival[index[0]],alpha,1);
		//cout<<"10,11-15,18   "<<Classe_analisi::Getanalisi()->nentry<<endl;
}
		
if((TMath::Nint(evento.z[index[1]])==9&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==16&&evento.rcoqf[index[0]]==1) || (TMath::Nint(evento.z[index[1]])==9&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==17&&evento.rcoqf[index[0]]==1)){
				    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
		fillh("hnz18alan_sm1",(float)ival[index[1]],alpha,1);
		fillh("hnz18alan_bg1",(float)ival[index[0]],alpha,1);
		fillh("hnz18alan_cp",(float)ivalc,1);
		//cout<<"9-16,17   "<<Classe_analisi::Getanalisi()->nentry<<endl;
}		
if((TMath::Nint(evento.z[index[1]])==9&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==20&&evento.rcoqf[index[0]]==1) || (TMath::Nint(evento.z[index[1]])==9&&evento.rcoqf[index[1]]==1&&(TMath::Nint(evento.z[index[0]]))==21&&evento.rcoqf[index[0]]==1)){
				    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
		fillh("hnz18alan_sm2",(float)ival[index[1]],alpha,1);
		fillh("hnz18alan_bg2",(float)ival[index[0]],alpha,1);
		fillh("hnz18alan_cp",(float)ivalc,1);
		//cout<<"9-20,21   "<<Classe_analisi::Getanalisi()->nentry<<endl;
}		
//

//Giu 26/10/17
if(evento.rcoqf[index[1]]==1&&evento.rcoqf[index[0]]==1){
int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);

	if(TMath::Nint(evento.z[index[1]])==5)fillh("hnz18cp_z5",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==6)fillh("hnz18cp_z6",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==7)fillh("hnz18cp_z7",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==8)fillh("hnz18cp_z8",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==9)fillh("hnz18cp_z9",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==10)fillh("hnz18cp_z10",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==11)fillh("hnz18cp_z11",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==12)fillh("hnz18cp_z12",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==13)fillh("hnz18cp_z13",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==14)fillh("hnz18cp_z14",(float)ivalc,1);
	if(TMath::Nint(evento.z[index[1]])==15)fillh("hnz18cp_z15",(float)ivalc,1);
}		
//


	if(alpha<90)
	  {
	    fillh("hnz18f1lt90",(float)ival[index[0]],1);
	    if(zsomma>=25 && zsomma<=35)
	      {
	    if(evento.z[index[1]]<5)
	      {
		fillh("hnz18f1lt90b1",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==5)
	      {
		fillh("hnz18f1lt90f5",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==6)
	      {
		fillh("hnz18f1lt90f6",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==7)
	      {
		fillh("hnz18f1lt90f7",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==8)
	      {
		fillh("hnz18f1lt90f8",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==9)
	      {
		fillh("hnz18f1lt90f9",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[index[1]])==10)
	      {
		fillh("hnz18f1lt90f10",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==11)
	      {
		fillh("hnz18f1lt90f11",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==12)
	      {
		fillh("hnz18f1lt90f12",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==13)
	      {
		fillh("hnz18f1lt90f13",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[index[1]])==14)
	      {
		fillh("hnz18f1lt90f14",(float)ival[index[0]],1);	
	      }

	      }

	  }
      }

    if(evento.rcoqf[index[1]]==1)
      {
	fillh("hnz18f2",(float)ival[index[1]],1);
	fillh("hnz18f2alpha",(float)ival[index[1]],alpha,1);

//Giu 5/10/17
	if(alpha<60)fillh("hnz18f2alphaleq60",(float)ival[index[1]],1);
	if(alpha<100)fillh("hnz18f2alphageq100",(float)ival[index[1]],1);
	if(eta<0.3)fillh("hnz18f2etaleq03",(float)ival[index[1]],1);
	if(eta>0.4)fillh("hnz18f2etageq04",(float)ival[index[1]],1);	
//



	if(alpha<90)
	  {
	    fillh("hnz18f2lt90",(float)ival[index[1]],1);
	  }
      }


			if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1)
			  {
			    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
			    fillh("hznfis18",(float)ivalc,1);
			  }
			
			fillh("hz1z2fis18",evento.z[index[0]],evento.z[index[1]],1);
			

			    fillh("hzvfiscoppiavzpos18",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);


	if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1){
			    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
//Giu 13/06/17
if(eta<0.3)
fillh("hznfiscoppiasmall",(float)ivalc,1);
else if(eta>0.3&&eta<0.6)
fillh("hznfiscoppiamid",(float)ivalc,1);
else if(eta>0.6)
fillh("hznfiscoppiabig",(float)ivalc,1);
//finish

//Giu 5/10/17
	if(alpha<60)fillh("hnz18coppiaalphaleq60",(float)ivalc,1);
	if(alpha<100)fillh("hnz18coppiaalphageq100",(float)ivalc,1);
	if(eta<0.25)fillh("hnz18coppiaetaleq03",(float)ivalc,1);
	if(eta>0.45)fillh("hnz18coppiaetageq04",(float)ivalc,1);	
//

}


   }//ifis18==1




 

 //cout<<iqp<<" "<<idic<<" "<<ifis<<endl;
 fillh("hcond",(float)iqp,(float)idic,1);
 fillh("hcond2",(float)iqp,(float)ifis,1);





 
 


 float nsuzlcp=0;
 float nt=0;
 float zt=0;
 float ntneg=0;
 float ztneg=0;
 if(iqp==1)
   {
  if(Classe_analisi::Getanalisi()->tipo_analisi<20)
    {
      fillh("hbiqp",mcevent.par_urto,1);
    }
   }
 if(idic==1)
   {
  if(Classe_analisi::Getanalisi()->tipo_analisi<20)
    {
      fillh("hbidic",mcevent.par_urto,1);
    }
   }

 if(iqp==1 || idic==1)
   {
     if(evento.rcoqf[index[0]]==1)
       {
	 fillh("hdeltaealqp",evento.z[index[0]],evento.a[index[0]]-Classe_formule::EAL(evento.z[index[0]]),1);
       }


     for(int j=0;j<evento.moltepl;j++)
       {

	 Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	 if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0])
	   {
	    
	     fillh("hvparzidiciqp",evento.z[j],vpar,1);


	   }

	     if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1)
	       //	 if(j!=index[0] && evento.rcoqf[j]==1)
	   {

//		if(TMath::Nint(evento.z[j])==3&&TMath::Nint(evento.a[j])==9) printf("Entro nell'evento %lld %f %d %d\n",Classe_analisi::Getanalisi()->nentry,mcevent.par_urto, mcevent.moltepl,TMath::Nint(evento.a[j]));
	     fillh("hnzconqp",(float)ival[j],1);
	     fillh("hnzconqpv",(float)ival[j],vpar,1);
//Giu
//		printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);
	if(evento.z[j]==3 && vpar>10 && vpar<40){ 
		fillh("htest",evento.rcocode[j], ibloc[j]*100+iqua[j]*10+itel[j],1);
		fillh("htest2",evento.moltepl, ibloc[j]*100+iqua[j]*10+itel[j],1);
		}
//	     
	     
	     float alphapic=Alpha(vpcm[index[0]],vpcm[j],evento.a[index[0]],evento.a[j]);
	     fillh("hnzconqpalpha",(float)ival[j],alphapic,1);
	     if(evento.rcoqf[index[0]]==1)
	       {
		 if(alphapic<90)
		   {
		     fillh("hnzqplt90",(float)ival[index[0]],1);
		     if(TMath::Nint(evento.z[j])==3) 
	      {
		fillh("hnzqplt90f3",(float)ival[index[0]],1);	
	      }
		     if(TMath::Nint(evento.z[j])==4)
	      {
		fillh("hnzqplt90f4",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==5)
	      {
		fillh("hnzqplt90f5",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==6)
	      {
		fillh("hnzqplt90f6",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==7)
	      {
		fillh("hnzqplt90f7",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==8)
	      {
		fillh("hnzqplt90f8",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==9)
	      {
		fillh("hnzqplt90f9",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==10)
	      {
		fillh("hnzqplt90f10",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[11])==11)
	      {
		fillh("hnzqplt90f11",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==12)
	      {
		fillh("hnzqplt90f12",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==13)
	      {
		fillh("hnzqplt90f13",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==14)
	      {
		fillh("hnzqplt90f14",(float)ival[index[0]],1);	
	      }


		   }

	       } 

	   }
	    

	 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchiadictot",vpar,vperp,1);
	   }
	 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==1)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchipdictot",vpar,vperp,1);
	   }
if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar>vpcm[index[0]][2])
  //	 if(j!=index[0] && vpar>vpcm[index[0]][2])
	   {
	     fillh("hnzvmagqp",(float)ival[j],1);

	     if(evento.z[j]<3)
	       {
		 nt=nt+evento.a[j]-evento.z[j];
		 zt=zt+evento.z[j];
	       }
	   }

if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<vpcm[index[0]][2])
 
	   {
	     fillh("hnzvminqp",(float)ival[j],1);
	     if(evento.z[j]<3)
	       {
		 ntneg=ntneg+evento.a[j]-evento.z[j];
		 ztneg=ztneg+evento.z[j];
	       }
	   }

if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<vpcm[index[0]][2]&&vpar>=0)
  {
fillh("hnzvminqpvpos",(float)ival[j],1);
  }
 if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<0)
  {
fillh("hnzvneg",(float)ival[j],1);
  }


       }


     if(zt>0)
       {
     fillh("hnsuzlcp",nt/zt,1);
     //   cout<<nt/zt<<endl;
       }
     if(ztneg>0)
       {
     fillh("hnsuzlcpneg",nt/zt,1);
     //   cout<<nt/zt<<endl;
       }
   }//iqp || idic



 if(iqp18==1 || idic18==1)
   {

//Giu
fillh("hmoltdic3frag",(float)nfrag,1);
	if(idic18==1 && evento.moltepl>2){
	     	Classe_formule::vparvperp(vpcm[index[2]],vpcm[index[0]],&vpar,&vperp);
		fillh("hnzdic18frag3",(float)ival[index[2]],1);
		fillh("hnzdic18frag3v",(float)ival[index[2]],vpar,1);
	}
//



     for(int j=0;j<evento.moltepl;j++)
       {

	 Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     
	     if(((idic18==1 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1)
	       //	 if(j!=index[0] && evento.rcoqf[j]==1)
	   {
	     fillh("hvparzidiciqp18",evento.z[j],vpar,1);

	     fillh("hnzconqp18",(float)ival[j],1);
	     fillh("hnzconqp18v",(float)ival[j],vpar,1);
	     float alphapic=Alpha(vpcm[index[0]],vpcm[j],evento.a[index[0]],evento.a[j]);
	     fillh("hnzconqp18alpha",(float)ival[j],alphapic,1);
	     if(evento.rcoqf[index[0]]==1)
	       {
		 if(alphapic<90)
		   {
		     fillh("hnzqp18lt90",(float)ival[index[0]],1);
		     if(TMath::Nint(evento.z[j])==3) 
	      {
		fillh("hnzqp18lt90f3",(float)ival[index[0]],1);	
	      }
		     if(TMath::Nint(evento.z[j])==4)
	      {
		fillh("hnzqp18lt90f4",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==5)
	      {
		fillh("hnzqp18lt90f5",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==6)
	      {
		fillh("hnzqp18lt90f6",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==7)
	      {
		fillh("hnzqp18lt90f7",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==8)
	      {
		fillh("hnzqp18lt90f8",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==9)
	      {
		fillh("hnzqp18lt90f9",(float)ival[index[0]],1);	
	      }

	    if(TMath::Nint(evento.z[j])==10)
	      {
		fillh("hnzqp18lt90f10",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[11])==11)
	      {
		fillh("hnzqp18lt90f11",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==12)
	      {
		fillh("hnzqp18lt90f12",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==13)
	      {
		fillh("hnzqp18lt90f13",(float)ival[index[0]],1);	
	      }
	    if(TMath::Nint(evento.z[j])==14)
	      {
		fillh("hnzqp18lt90f14",(float)ival[index[0]],1);	
	      }


		   }
	       }
	   }
	 if(TMath::Nint(evento.z[j])==2 &&TMath::Nint(evento.a[j])==4)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchiadictot18",vpar,vperp,1);
	   }
	 if(TMath::Nint(evento.z[j])==1 &&TMath::Nint(evento.a[j])==1)
	   {
	     Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     fillh("hocchipdictot18",vpar,vperp,1);
	   }
if(((idic18==1 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar>vpcm[index[0]][2])
  //	 if(j!=index[0] && vpar>vpcm[index[0]][2])
	   {
	     fillh("hnzvmagqp18",(float)ival[j],1);

	     
	   }

if(((idic18==1 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<vpcm[index[0]][2])
 
	   {
	     fillh("hnzvminqp18",(float)ival[j],1);
	     
	   }

if(((idic==18 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<vpcm[index[0]][2]&&vpar>=0)
  {
fillh("hnzvminqpvpos18",(float)ival[j],1);
  }
 if(((idic18==1 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1&& vpar<0)
  {
fillh("hnzvneg18",(float)ival[j],1);
  }


       }


 
   }//iqp18 || idic18

 if(iqp==1)
   {
     fillh("hthecmvcmiqp",evento.vpartcm[index[0]],evento.thetacm[index[0]],1);

     for(int j=0;j<evento.moltepl;j++)
       {
	 if(j!=index[0])
	   {
	     fillh("hzalphaiqp",Alpha(vpcm[index[0]],vpcm[j],evento.a[index[0]],evento.a[j]),evento.z[j],1);
	   }
       }
   }
 if(idic==1)
   {
      for(int j=0;j<evento.moltepl;j++)
       {
	 if(j!=index[0] && j!=index[1])
	   {
	     fillh("hzalphaidic",Alpha(vpcm[index[0]],vpcm[j],evento.a[index[0]],evento.a[j]),evento.z[j],1);
	   }
       }
      if(nfrag==3)
	{
	  fillh("hzalphaidicneck",Alpha(vpcm[index[0]],vpcm[index[2]],evento.a[index[0]],evento.a[index[2]]),evento.z[index[2]],1);
	}    
   }



}
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void nsuzv(TH2F *h,TGraphErrors *gnz[100]);
void nsuzJed(TH2F *h,TGraphErrors *gnz[100]);
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2);
void Classe_analisi::RoutineFinale()
{
  cout<<"entro nella routine finale"<<endl;
TGraphErrors *gm[100];
TGraphErrors *gznvimfdic[100];
TGraphErrors *gznvimf18dic[100];
TGraphErrors *gnzconqpv[100];
TGraphErrors *gnzconqp18v[100];
TGraphErrors *gnzconqpalpha[100];
TGraphErrors *gnzconqp18alpha[100];
TGraphErrors *gnzf1alpha[100];
TGraphErrors *gnzf2alpha[100];
TGraphErrors *gnz18f1alpha[100];
TGraphErrors *gnz18f2alpha[100];

//Giu
TGraphErrors *gnfissuz3fragv[100];
TGraphErrors *gnzdic18frag3v[100];
TGraphErrors *gnz18fsmall5eta6[100];
TGraphErrors *gnz18fsmall6eta6[100];
TGraphErrors *gnz18fsmall5eta2[100];
TGraphErrors *gnz18fsmall6eta2[100];
TGraphErrors *gnz18fbig5eta06[100];
TGraphErrors *gnz18fbig5eta02[100];
TGraphErrors *gnz18fbig6eta06[100];
TGraphErrors *gnz18fbig6eta02[100];
TGraphErrors *gnz18sym_bg[100];
TGraphErrors *gnz18sym_sm[100];
TGraphErrors *gnz18asym_bg1[100];
TGraphErrors *gnz18asym_sm1[100];
TGraphErrors *gnz18asym_bg2[100];
TGraphErrors *gnz18asym_sm2[100];
TGraphErrors *gnz18alan_bg1[100];
TGraphErrors *gnz18alan_sm1[100];
TGraphErrors *gnz18alan_bg2[100];
TGraphErrors *gnz18alan_sm2[100];

//

  for(int j=0;j<100;j++) 
    {
      gm[j]=0;
      gznvimfdic[j]=0;
      gnzconqpv[j]=0;
      gnzconqp18v[j]=0;
      gznvimf18dic[j]=0;
     gnzconqpalpha[j]=0;
     gnzconqp18alpha[j]=0;
     gnzf2alpha[j]=0;
     gnzf1alpha[j]=0;
     gnz18f2alpha[j]=0;
     gnz18f1alpha[j]=0;
	gnz18sym_bg[j]=0;
	gnz18sym_sm[j]=0;
	gnz18asym_bg1[j]=0;
	gnz18asym_sm1[j]=0;
	gnz18asym_bg2[j]=0;
	gnz18asym_sm2[j]=0;
	gnz18alan_bg1[j]=0;
	gnz18alan_sm1[j]=0;
	gnz18alan_bg2[j]=0;
	gnz18alan_sm2[j]=0;

//Giu
gnfissuz3fragv[j]=0;gnfissuz3fragv[j]=getg(Form("gnfissuz3fragv_l%d",j));
gnfissuz3fragv[j]=0;gnfissuz3fragv[j]=getg(Form("gnfissuz3fragv_l%d",j));
gnfissuz3fragv[j]=0;gnfissuz3fragv[j]=getg(Form("gnfissuz3fragv_l%d",j));
gnzdic18frag3v[j]=0;gnzdic18frag3v[j]=getg(Form("gnzdic18frag3v_l%d",j));
gnz18fsmall5eta6[j]=0;gnz18fsmall5eta6[j]=getg(Form("gnz18fsmall5eta6_l%d",j));
gnz18fsmall6eta6[j]=0;gnz18fsmall6eta6[j]=getg(Form("gnz18fsmall6eta6_l%d",j));
gnz18fsmall5eta2[j]=0;gnz18fsmall5eta2[j]=getg(Form("gnz18fsmall5eta2_l%d",j));
gnz18fsmall6eta2[j]=0;gnz18fsmall6eta2[j]=getg(Form("gnz18fsmall6eta2_l%d",j));
gnz18fbig5eta06[j]=0;gnz18fbig5eta06[j]=getg(Form("gnz18fbig5eta06_l%d",j));
gnz18fbig5eta02[j]=0;gnz18fbig5eta02[j]=getg(Form("gnz18fbig5eta02_l%d",j));
gnz18fbig6eta06[j]=0;gnz18fbig6eta06[j]=getg(Form("gnz18fbig6eta06_l%d",j));
gnz18fbig6eta02[j]=0;gnz18fbig6eta02[j]=getg(Form("gnz18fbig6eta02_l%d",j));
//
    }

  for(int j=1;j<26;j++) 
    { 
      gznvimfdic[j]=getg(Form("gznvimfdic%d",j));
      gznvimf18dic[j]=getg(Form("gznvimf18dic%d",j));
gnzconqpv[j]=getg(Form("gnzconqpv%d",j));
gnzconqpalpha[j]=getg(Form("gnzconqpalpha%d",j));
gnzconqp18alpha[j]=getg(Form("gnzconqp18alpha%d",j));
gnzconqp18v[j]=getg(Form("gnzconqp18v%d",j));
 gnzf1alpha[j]=getg(Form("gnzf1alpha%d",j));
 gnzf2alpha[j]=getg(Form("gnzf2alpha%d",j));
 gnz18f1alpha[j]=getg(Form("gnz18f1alpha%d",j));
 gnz18f2alpha[j]=getg(Form("gnz18f2alpha%d",j));
	gnz18sym_bg[j]=getg(Form("gnz18sym_bg_l%d",j));
	gnz18sym_sm[j]=getg(Form("gnz18sym_sm_l%d",j));
	gnz18asym_bg1[j]=getg(Form("gnz18asym_bg1_l%d",j));
	gnz18asym_sm1[j]=getg(Form("gnz18asym_sm1_l%d",j));
	gnz18asym_bg2[j]=getg(Form("gnz18asym_bg2_l%d",j));
	gnz18asym_sm2[j]=getg(Form("gnz18asym_sm2_l%d",j));
	gnz18alan_bg1[j]=getg(Form("gnz18alan_bg1_l%d",j));
	gnz18alan_sm1[j]=getg(Form("gnz18alan_sm1_l%d",j));
	gnz18alan_bg2[j]=getg(Form("gnz18alan_bg2_l%d",j));
	gnz18alan_sm2[j]=getg(Form("gnz18alan_sm2_l%d",j));

    }

nsuz(geth1("hznbig"),getg("gnzbig"),gm);
nsuz(geth1("hzngrossone"),getg("gnzgrossone"),gm);
nsuz(geth1("hzngrossone18"),getg("gnzgrossone18"),gm);
nsuz(geth1("hnzimfdic"),getg("gnzimfdic"),gm);
nsuz(geth1("hzndic"),getg("gnzdic"),gm);
nsuz(geth1("hzndic18"),getg("gnzdic18"),gm);
nsuz(geth1("hzndic18qt"),getg("gnzdic18qt"),gm);
nsuz(geth1("hzndic0"),getg("gnzdic0"),gm);
nsuz(geth1("hznfis"),getg("gnzfis"),gm);
nsuz(geth1("hznfis18"),getg("gnzfis18"),gm);
nsuz(geth1("hnz18alan_cp"),getg("gnz18alan_cp"),gm);
nsuz(geth1("hnz18asym_cp"),getg("gnz18asym_cp"),gm);
nsuz(geth1("hnz18sym_cp"),getg("gnz18sym_cp"),gm);
nsuz(geth1("hnz18cp_z5"),getg("gnz18cp_z5"),gm);
nsuz(geth1("hnz18cp_z6"),getg("gnz18cp_z6"),gm);
nsuz(geth1("hnz18cp_z7"),getg("gnz18cp_z7"),gm);
nsuz(geth1("hnz18cp_z8"),getg("gnz18cp_z8"),gm);
nsuz(geth1("hnz18cp_z9"),getg("gnz18cp_z9"),gm);
nsuz(geth1("hnz18cp_z10"),getg("gnz18cp_z10"),gm);
nsuz(geth1("hnz18cp_z11"),getg("gnz18cp_z11"),gm);
nsuz(geth1("hnz18cp_z12"),getg("gnz18cp_z12"),gm);
nsuz(geth1("hnz18cp_z13"),getg("gnz18cp_z13"),gm);
nsuz(geth1("hnz18cp_z14"),getg("gnz18cp_z14"),gm);
nsuz(geth1("hnz18cp_z15"),getg("gnz18cp_z15"),gm);
//Giu 13/06/17
nsuz(geth1("hznfiscoppiasmall"),getg("gnzfiscoppiasmall"),gm);
nsuz(geth1("hznfiscoppiamid"),getg("gnzfiscoppiamid"),gm);
nsuz(geth1("hznfiscoppiabig"),getg("gnzfiscoppiabig"),gm);
nsuz(geth1("hnfissuz3frag"),getg("gnfissuz3frag"),gm);
nsuz(geth1("hnfistuttitot"),getg("gnfistuttitot"),gm);
nsuz(geth1("hnzdic18frag3"),getg("gnzdic18frag3"),gm);
nsuzv(geth2("hnfissuz3fragv"),gnfissuz3fragv);
nsuzv(geth2("hnzdic18frag3v"),gnzdic18frag3v);
//
//Giu 13/06/17
nsuzJed(geth2("hnz18asym_sm1"),gnz18asym_sm1);
nsuzJed(geth2("hnz18asym_bg1"),gnz18asym_bg1);
nsuzJed(geth2("hnz18asym_sm2"),gnz18asym_sm2);
nsuzJed(geth2("hnz18asym_bg2"),gnz18asym_bg2);
nsuzJed(geth2("hnz18sym_sm"),gnz18sym_sm);
nsuzJed(geth2("hnz18sym_bg"),gnz18sym_bg);
nsuzJed(geth2("hnz18alan_sm1"),gnz18alan_sm1);
nsuzJed(geth2("hnz18alan_bg1"),gnz18alan_bg1);
nsuzJed(geth2("hnz18alan_sm2"),gnz18alan_sm2);
nsuzJed(geth2("hnz18alan_bg2"),gnz18alan_bg2);
nsuzv(geth2("hnz18fsmall5eta6"),gnz18fsmall5eta6);
nsuzv(geth2("hnz18fsmall6eta6"),gnz18fsmall6eta6);
nsuzv(geth2("hnz18fsmall5eta2"),gnz18fsmall5eta2);
nsuzv(geth2("hnz18fsmall6eta2"),gnz18fsmall6eta2);
nsuzv(geth2("hnz18fbig5eta06"), gnz18fbig5eta06);
nsuzv(geth2("hnz18fbig5eta02"), gnz18fbig5eta02);
nsuzv(geth2("hnz18fbig6eta06"), gnz18fbig6eta06);
nsuzv(geth2("hnz18fbig6eta02"), gnz18fbig6eta02);
nsuz(geth1("hnz18f1alphaleq60"),getg("gnz18f1alphaleq60"),gm);
nsuz(geth1("hnz18f1alphageq100"),getg("gnz18f1alphageq100"),gm);
nsuz(geth1("hnz18f1etaleq03"),getg("gnz18f1etaleq03"),gm);
nsuz(geth1("hnz18f1etageq04"),getg("gnz18f1etageq04"),gm);
nsuz(geth1("hnz18f2alphaleq60"),getg("gnz18f2alphaleq60"),gm);
nsuz(geth1("hnz18f2alphageq100"),getg("gnz18f2alphageq100"),gm);
nsuz(geth1("hnz18f2etaleq03"),getg("gnz18f2etaleq03"),gm);
nsuz(geth1("hnz18f2etageq04"),getg("gnz18f2etageq04"),gm);
nsuz(geth1("hnz18coppiaalphaleq60"),getg("gnz18coppiaalphaleq60"),gm);
nsuz(geth1("hnz18coppiaalphageq100"),getg("gnz18coppiaalphageq100"),gm);
nsuz(geth1("hnz18coppiaetaleq03"),getg("gnz18coppiaetaleq03"),gm);
nsuz(geth1("hnz18coppiaetageq04"),getg("gnz18coppiaetageq04"),gm);
//
 nsuz(geth1("hzniqp"),getg("gnziqp"),gm);
 nsuz(geth1("hzniqp18"),getg("gnziqp18"),gm);
 nsuz(geth1("hnzconqp"),getg("gnzconqp"),gm);
  nsuz(geth1("hnzconqp18"),getg("gnzconqp18"),gm);
 nsuz(geth1("hnzvmagqp"),getg("gnzvmagqp"),gm);
 nsuz(geth1("hnzvminqp"),getg("gnzvminqp"),gm);
 nsuz(geth1("hnzvminqpvpos"),getg("gnzvminqpvpos"),gm);
 nsuz(geth1("hnzvneg"),getg("gnzvneg"),gm);

 nsuz(geth1("hnzvmagqp18"),getg("gnzvmagqp18"),gm);
 nsuz(geth1("hnzvminqp18"),getg("gnzvminqp18"),gm);
 nsuz(geth1("hnzvminqpvpos18"),getg("gnzvminqpvpos18"),gm);
 nsuz(geth1("hnzvneg18"),getg("gnzvneg18"),gm);

 nsuz(geth1("hnztutti"),getg("gnztutti"),gm);
 nsuz(geth1("hnzf1"),getg("gnzf1"),gm);
  nsuz(geth1("hnzf2"),getg("gnzf2"),gm);
 nsuz(geth1("hnzf1lt90"),getg("gnzf1lt90"),gm);
  nsuz(geth1("hnzf2lt90"),getg("gnzf2lt90"),gm);
  nsuz(geth1("hnzf1lt90b1"),getg("gnzf1lt90b1"),gm);
  nsuz(geth1("hnzf1lt90f5"),getg("gnzf1lt90f5"),gm);
  nsuz(geth1("hnzf1lt90f6"),getg("gnzf1lt90f6"),gm);
  nsuz(geth1("hnzf1lt90f7"),getg("gnzf1lt90f7"),gm);
  nsuz(geth1("hnzf1lt90f8"),getg("gnzf1lt90f8"),gm);
  nsuz(geth1("hnzf1lt90f9"),getg("gnzf1lt90f9"),gm);
  nsuz(geth1("hnzf1lt90f10"),getg("gnzf1lt90f10"),gm);
  nsuz(geth1("hnzf1lt90f11"),getg("gnzf1lt90f11"),gm);
  nsuz(geth1("hnzf1lt90f12"),getg("gnzf1lt90f12"),gm);
  nsuz(geth1("hnzf1lt90f13"),getg("gnzf1lt90f13"),gm);
  nsuz(geth1("hnzf1lt90f14"),getg("gnzf1lt90f14"),gm);


 nsuz(geth1("hnz18f1"),getg("gnz18f1"),gm);
  nsuz(geth1("hnz18f2"),getg("gnz18f2"),gm);
 nsuz(geth1("hnz18f1lt90"),getg("gnz18f1lt90"),gm);
  nsuz(geth1("hnz18f2lt90"),getg("gnz18f2lt90"),gm);
  nsuz(geth1("hnz18f1lt90b1"),getg("gnz18f1lt90b1"),gm);
  nsuz(geth1("hnz18f1lt90f5"),getg("gnz18f1lt90f5"),gm);
  nsuz(geth1("hnz18f1lt90f6"),getg("gnz18f1lt90f6"),gm);
  nsuz(geth1("hnz18f1lt90f7"),getg("gnz18f1lt90f7"),gm);
  nsuz(geth1("hnz18f1lt90f8"),getg("gnz18f1lt90f8"),gm);
  nsuz(geth1("hnz18f1lt90f9"),getg("gnz18f1lt90f9"),gm);
  nsuz(geth1("hnz18f1lt90f10"),getg("gnz18f1lt90f10"),gm);
  nsuz(geth1("hnz18f1lt90f11"),getg("gnz18f1lt90f11"),gm);
  nsuz(geth1("hnz18f1lt90f12"),getg("gnz18f1lt90f12"),gm);
  nsuz(geth1("hnz18f1lt90f13"),getg("gnz18f1lt90f13"),gm);
  nsuz(geth1("hnz18f1lt90f14"),getg("gnz18f1lt90f14"),gm);

  nsuz(geth1("hnzqplt90"),getg("gnzqplt90"),gm);
  nsuz(geth1("hnzqplt90f3"),getg("gnzqplt90f3"),gm);
  nsuz(geth1("hnzqplt90f4"),getg("gnzqplt90f4"),gm);

  nsuz(geth1("hnzqplt90f5"),getg("gnzqplt90f5"),gm);
  nsuz(geth1("hnzqplt90f6"),getg("gnzqplt90f6"),gm);
  nsuz(geth1("hnzqplt90f7"),getg("gnzqplt90f7"),gm);
  nsuz(geth1("hnzqplt90f8"),getg("gnzqplt90f8"),gm);
  nsuz(geth1("hnzqplt90f9"),getg("gnzqplt90f9"),gm);
  nsuz(geth1("hnzqplt90f10"),getg("gnzqplt90f10"),gm);
  nsuz(geth1("hnzqplt90f11"),getg("gnzqplt90f11"),gm);
  nsuz(geth1("hnzqplt90f12"),getg("gnzqplt90f12"),gm);
  nsuz(geth1("hnzqplt90f13"),getg("gnzqplt90f13"),gm);
  nsuz(geth1("hnzqplt90f14"),getg("gnzqplt90f14"),gm);

  nsuz(geth1("hnzqp18lt90"),getg("gnzqp18lt90"),gm);
  nsuz(geth1("hnzqp18lt90f3"),getg("gnzqp18lt90f3"),gm);
  nsuz(geth1("hnzqp18lt90f4"),getg("gnzqp18lt90f4"),gm);

  nsuz(geth1("hnzqp18lt90f5"),getg("gnzqp18lt90f5"),gm);
  nsuz(geth1("hnzqp18lt90f6"),getg("gnzqp18lt90f6"),gm);
  nsuz(geth1("hnzqp18lt90f7"),getg("gnzqp18lt90f7"),gm);
  nsuz(geth1("hnzqp18lt90f8"),getg("gnzqp18lt90f8"),gm);
  nsuz(geth1("hnzqp18lt90f9"),getg("gnzqp18lt90f9"),gm);
  nsuz(geth1("hnzqp18lt90f10"),getg("gnzqp18lt90f10"),gm);
  nsuz(geth1("hnzqp18lt90f11"),getg("gnzqp18lt90f11"),gm);
  nsuz(geth1("hnzqp18lt90f12"),getg("gnzqp18lt90f12"),gm);
  nsuz(geth1("hnzqp18lt90f13"),getg("gnzqp18lt90f13"),gm);
  nsuz(geth1("hnzqp18lt90f14"),getg("gnzqp18lt90f14"),gm);



nsuzv(geth2("hznvimfdic"),gznvimfdic);

nsuzv(geth2("hznvimfdic18"),gznvimf18dic);
 
nsuzv(geth2("hnzconqpv"),gnzconqpv);
nsuzv(geth2("hnzconqp18v"),gnzconqp18v);

nsuzv(geth2("hnzconqpalpha"),gnzconqpalpha);
nsuzv(geth2("hnzf1alpha"),gnzf1alpha);
nsuzv(geth2("hnzf2alpha"),gnzf2alpha);

nsuzv(geth2("hnz18f1alpha"),gnz18f1alpha);
nsuzv(geth2("hnz18f2alpha"),gnz18f2alpha);

 for(int j=0;j<100;j++)
   {
     if(gnzconqpalpha[j]!=0)
       {
	 gnzconqpalpha[j]->GetXaxis()->SetTitle("alpha"); 
       }
     if(gnzf1alpha[j]!=0)
       {
	 gnzf1alpha[j]->GetXaxis()->SetTitle("alpha"); 
       }
     if(gnzf2alpha[j]!=0)
       {
	 gnzf2alpha[j]->GetXaxis()->SetTitle("alpha"); 
       }
    if(gnz18f1alpha[j]!=0)
       {
	 gnz18f1alpha[j]->GetXaxis()->SetTitle("alpha"); 
       }
     if(gnz18f2alpha[j]!=0)
       {
	 gnz18f2alpha[j]->GetXaxis()->SetTitle("alpha"); 
       }

   }


 TGraphErrors *grap[7];
 for(int j=0;j<7;j++) 
    { 
     
      grap[j]=getg(Form("grap%d",j));
    }
  for(int j=0;j<geth2("hnzconqp18v")->GetYaxis()->GetNbins();j++)
   {
   
     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(102.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1)>0)
       {
	 grap[0]->SetPoint(grap[0]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(102.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1));
grap[0]->SetPointError(grap[0]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(102.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(102.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1))));

       }
     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(101.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1)>0)
       {
	 grap[1]->SetPoint(grap[1]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(101.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1));
grap[1]->SetPointError(grap[1]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(101.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(101.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(100.),j+1))));
       }
     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(202.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(204.),j+1)>0)
       {
     grap[2]->SetPoint(grap[2]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(204.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(202.),j+1));
     grap[2]->SetPointError(grap[2]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(204.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(202.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(204.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(202.),j+1))));
       }

     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(304.),j+1)>0)
       {
     grap[3]->SetPoint(grap[3]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(304.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1));
     grap[3]->SetPointError(grap[3]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(304.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(304.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1))));
       }

     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(405.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1)>0)
       {
     grap[4]->SetPoint(grap[4]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(405.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1));
     grap[4]->SetPointError(grap[4]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(405.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(405.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(403.),j+1))));
       }

     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(303.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(206.),j+1)>0)
       {
     grap[5]->SetPoint(grap[5]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(306.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(206.),j+1));
 grap[5]->SetPointError(grap[5]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(306.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(206.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(306.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(206.),j+1))));
       }

     if(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(103.),j+1)>0 && geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(203.),j+1)>0)
       {
     grap[6]->SetPoint(grap[6]->GetN(),geth2("hnzconqp18v")->GetYaxis()->GetBinCenter(j+1),geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(103.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(203.),j+1));
 grap[6]->SetPointError(grap[6]->GetN()-1,0.001,(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(103.),j+1)/geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(203.),j+1))*(1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(103.),j+1))+1/sqrt(geth2("hnzconqp18v")->GetBinContent(geth2("hnzconqp18v")->GetXaxis()->FindBin(203.),j+1))));
       }

   }
 grap[0]->GetXaxis()->SetTitle("vpar");
 grap[0]->GetYaxis()->SetTitle("t/p");
 grap[1]->GetXaxis()->SetTitle("vpar");
 grap[1]->GetYaxis()->SetTitle("d/p");
 grap[2]->GetXaxis()->SetTitle("vpar");
 grap[2]->GetYaxis()->SetTitle("6He/alpha");
 grap[3]->GetXaxis()->SetTitle("vpar");
 grap[3]->GetYaxis()->SetTitle("Li7/Be7");
 grap[4]->GetXaxis()->SetTitle("vpar");
 grap[4]->GetYaxis()->SetTitle("Be9/Be7");
 grap[5]->GetXaxis()->SetTitle("vpar");
 grap[5]->GetYaxis()->SetTitle("Li6/He6");
 grap[6]->GetXaxis()->SetTitle("vpar");
 grap[6]->GetYaxis()->SetTitle("t/He3");

 float rap=geth1("hneventi")->GetBinContent(geth1("hneventi")->FindBin(24.))/(geth1("hneventi")->GetBinContent(geth1("hneventi")->FindBin(24.))+geth1("hneventi")->GetBinContent(geth1("hneventi")->FindBin(17.))+geth1("hneventi")->GetBinContent(geth1("hneventi")->FindBin(22.)));
  
 ordinisup(geth1("hnfistuttitot"),getg("gnfistuttitot"),getg("gnfistuttitotsigma"));
 ordinisup(geth1("hzniqp18"),getg("gnziqp18"),getg("gnziqpsigma18"));
 ordinisup(geth1("hzngrossone"),getg("gnzgrossone"),getg("gnzgrossonesigma"));
 ordinisup(geth1("hzngrossone18"),getg("gnzgrossone18"),getg("gnzgrossonesigma18"));
 ordinisup(geth1("hnzf1"),getg("gnzf1"),getg("gnzf1sigma"));
 ordinisup(geth1("hnzf2"),getg("gnzf2"),getg("gnzf2sigma"));
 ordinisup(geth1("hzndic"),getg("gnzdic"),getg("gnzdicsigma"));
 ordinisup(geth1("hzndic18"),getg("gnzdic18"),getg("gnzdicsigma18"));
 ordinisup(geth1("hzndic18qt"),getg("gnzdic18qt"),getg("gnzdicsigma18qt"));
 ordinisup(geth1("hnztutti"),getg("gnztutti"),getg("gnztuttisigma"));
 ordinisup(geth1("hznfis"),getg("gnzfis"),getg("gnzfissigma"));
 ordinisup(geth1("hznfis18"),getg("gnzfis18"),getg("gnzfissigma18"));
 ordinisup(geth1("hnzconqp"),getg("gnzconqp"),getg("gnzconqpsigma"));
 ordinisup(geth1("hnzconqp18"),getg("gnzconqp18"),getg("gnzconqpsigma18"));
 cout<<"rapporto fissioni/non fissioni="<<rap<<endl;

  return;
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


 for(int iz=0;iz<=100;iz++)
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




void nsuzv(TH2F *h,TGraphErrors *gnz[100])
{
  int nbx=h->GetYaxis()->GetNbins();
  int totz[100][nbx];
  int totzn[100][100][nbx];
float nzval[100][nbx];
 float 	 enzval[100][nbx];
 float	 err[100][nbx];
  
      for(int iz=0;iz<100;iz++)
	{
	  for(int iv=0;iv<nbx;iv++)
	    {
	      totz[iz][iv]=0;
nzval[iz][iv]=0;
 	 enzval[iz][iv]=0;
 	 err[iz][iv]=0;
	  for(int ia=0;ia<100;ia++)
	    {
	      totzn[iz][ia][iv]=0;
	 
	    }
	    }
	}
    

  for(int k=0;k<h->GetXaxis()->GetNbins();k++)
       {
	 for(int iv=0;iv<h->GetYaxis()->GetNbins();iv++)
	   {
	     int icount=h->GetBinContent(k+1,iv+1);
	
	 if(icount>0)
	   {

	     int iz=(k+1)/100;
	     int in=k+1-iz*100-1;
	   
	            	 totz[iz][iv]+=icount;
			 totzn[iz][in][iv]+=icount; 
	   }
       }
       }


 for(int iz=0;iz<=25;iz++)
    {
      for(int iv=0;iv<nbx;iv++)
	{
	  float velo=h->GetYaxis()->GetBinCenter(iv+1);
      if(totz[iz][iv]>5)//cambiato da 10 a 5 il minimo
	{

 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in][iv]>0)
 	       {

 		 nzval[iz][iv]=nzval[iz][iv]+totzn[iz][in][iv]*(float)in;
		
 	       }
 	   }
 	 if(nzval[iz][iv]>0)
 	   {
 	     float nsuz=nzval[iz][iv]/totz[iz][iv];
 	     for(int in=0;in<100;in++)
 	       {
 		 if(totzn[iz][in][iv]>0)
 		   {

 			 enzval[iz][iv]=enzval[iz][iv]+totzn[iz][in][iv]*(float)in*(float)in;


 		   }
 	       }

 	     enzval[iz][iv]=enzval[iz][iv]/totz[iz][iv];
 	     float nzvaldue=nsuz*nsuz;
 	     enzval[iz][iv]=sqrt((enzval[iz][iv]-nzvaldue)/totz[iz][iv]);


 	      err[iz][iv]=enzval[iz][iv]/((float)iz);
	      
	     

 	     nsuz=nsuz/(float)iz;
	     if(gnz[iz]!=0)
	       {
		 gnz[iz]->SetPoint(gnz[iz]->GetN(),velo,nsuz);
 	 gnz[iz]->SetPointError(gnz[iz]->GetN()-1,0.001,err[iz][iv]);
	       }
 	   }

	}//totz>10
	}//iv
      if(gnz[iz]!=0)
	{
gnz[iz]->GetXaxis()->SetTitle("vpar");
 gnz[iz]->GetYaxis()->SetTitle("N/Z");
 gnz[iz]->SetTitle(Form("Z=%d",iz));
	}
    }//iz=1,80

 

}


void nsuzJed(TH2F *h,TGraphErrors *gnz[100])
{
  int nbx=h->GetYaxis()->GetNbins();
  int totz[100][nbx];
  int totzn[100][100][nbx];
float nzval[100][nbx];
 float 	 enzval[100][nbx];
 float	 err[100][nbx];
  
      for(int iz=0;iz<100;iz++)
	{
	  for(int iv=0;iv<nbx;iv++)
	    {
	      totz[iz][iv]=0;
nzval[iz][iv]=0;
 	 enzval[iz][iv]=0;
 	 err[iz][iv]=0;
	  for(int ia=0;ia<100;ia++)
	    {
	      totzn[iz][ia][iv]=0;
	 
	    }
	    }
	}
    

  for(int k=0;k<h->GetXaxis()->GetNbins();k++)
       {
	 for(int iv=0;iv<h->GetYaxis()->GetNbins();iv++)
	   {
	     int icount=h->GetBinContent(k+1,iv+1);
	
	 if(icount>0)
	   {

	     int iz=(k+1)/100;
	     int in=k+1-iz*100-1;
	   
	            	 totz[iz][iv]+=icount;
			 totzn[iz][in][iv]+=icount; 
	   }
       }
       }


 for(int iz=0;iz<=25;iz++)
    {
      for(int iv=0;iv<nbx;iv++)
	{
	  float velo=h->GetYaxis()->GetBinCenter(iv+1);
      if(totz[iz][iv]>5)//cambiato da 10 a 5 il minimo
	{

 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in][iv]>0)
 	       {

 		 nzval[iz][iv]=nzval[iz][iv]+totzn[iz][in][iv]*(float)in;
		
 	       }
 	   }
 	 if(nzval[iz][iv]>0)
 	   {
 	     float nsuz=nzval[iz][iv]/totz[iz][iv];
 	     for(int in=0;in<100;in++)
 	       {
 		 if(totzn[iz][in][iv]>0)
 		   {

 			 enzval[iz][iv]=enzval[iz][iv]+totzn[iz][in][iv]*(float)in*(float)in;


 		   }
 	       }

 	     enzval[iz][iv]=enzval[iz][iv]/totz[iz][iv];
 	     float nzvaldue=nsuz*nsuz;
 	     enzval[iz][iv]=sqrt((enzval[iz][iv]-nzvaldue)/totz[iz][iv]);


 	      err[iz][iv]=enzval[iz][iv]/((float)iz);
	      
	     

 	     nsuz=nsuz/(float)iz;
	     if(gnz[iz]!=0)
	       {
		 gnz[iz]->SetPoint(gnz[iz]->GetN(),velo,(nsuz-1)/(nsuz+1));
 	 	 gnz[iz]->SetPointError(gnz[iz]->GetN()-1,0.001,2*(err[iz][iv]/pow((nsuz+1),2)));
	       }
 	   }

	}//totz>10
	}//iv
      if(gnz[iz]!=0)
	{
 gnz[iz]->GetXaxis()->SetTitle("#alpha");
 gnz[iz]->GetYaxis()->SetTitle("N/Z");
 gnz[iz]->SetTitle(Form("Z=%d",iz));
	}
    }//iz=1,80

 

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


#endif
