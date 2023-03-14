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


float Alpha(float *v1,float *v2)
{
  float vcc[3];
  for(int j=0;j<3;j++)
    {
      vcc[j]=v1[j]-v2[j];
    }
  float sca;
  Classe_formule::scaprod(vcc,v2,&sca);
  float alpha=57.296*TMath::ACos(sca/(Classe_formule::modulo(vcc)*Classe_formule::modulo(v2)));
  if(Classe_formule::modulo(vcc)==0)
    {
      cout<<"vcc=0"<<" "<<v1[0]<<" "<<v2[0]<<" "<<v1[1]<<" "<<v2[1]<<" "<<v1[2]<<" "<<v2[2]<<endl;
      cout<<Classe_analisi::Getanalisi()->nentry<<endl;
      return -1000;
 }
  return alpha;
				 
}


void Classe_evento::AnalisiPrincipale()
{
  //printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);
  // cout<<"molt="<<evento.moltepl<<endl;
  fillh("hneventi",1.,1);
  fillh("hmolt",evento.moltepl,1); 
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
   }


fillh("hneventi",3.,1);
     float thetarel=-1;
     float vrel[3];
     int idic=-1;
     int ifis=-1;
     int idic18=-1;
     int ifis18=-1;
     int iqp=-1;
     int iqp18=-1;
     int icent=-1;
     int icent18=-1;

     float vpar=-1000;
     float vperp=-1000;
 if(tflow>=8 && tflow<=30 && sumz>=12)
   {
     fillh("hzvsel1big",evento.vpartlab[index[0]],evento.z[index[0]],1);

     if(((nfrag==1 &&evento.z[index[0]]>=12) ||(nfrag>=2 && evento.z[index[1]]<5 &&evento.z[index[0]]>=12))&&evento.vpartcm_z[index[0]]>0)
       {
	 iqp=1;
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
       }//iqp
     if(((nfrag==1 &&evento.z[index[0]]>18) ||(nfrag>=2 && evento.z[index[1]]<5 &&evento.z[index[0]]>18))&&evento.vpartcm_z[index[0]]>0)
       {
	 iqp18=1;
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

     thetarel=-1;
	  if(nfrag>=2)
	    {
	     
	      Classe_formule::vrel(vpcm[index[1]],vpcm[index[0]],vrel);
	      thetarel=Classe_formule::thetarel(vpcm[index[1]],vpcm[index[0]]);  
	      fillh("hvreltherelsel1",Classe_formule::modulo(vrel),thetarel,1);
	      if(evento.z[index[1]]>=5)
		   {
		    fillh("hvreltherelsel1gt5",Classe_formule::modulo(vrel),thetarel,1); 
		   }
	    }//nfrag>=2

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
	 if(thetarel>=0) //vero con nfrag>=2
	   {
	     fillh("hvreltherelsel1vpos",Classe_formule::modulo(vrel),thetarel,1);

	      if(evento.z[index[1]]>=5)
		   {
		    fillh("hvreltherelsel1vposgt5",Classe_formule::modulo(vrel),thetarel,1); 
		    fillh("hz12therelsel1",evento.z[index[0]]+evento.z[index[1]],thetarel,1);
		    if(thetarel>=160)
		      {
			float mi=evento.a[index[0]]*evento.a[index[1]]/(evento.a[index[0]]+evento.a[index[1]]);
			float tke=0.5*mi*Classe_formule::amu*pow(Classe_formule::modulo(vrel),2)/pow(Classe_formule::cluce,2);
			if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzndic0",(float)ival[index[0]],1);
			  }
			if(evento.z[index[0]]>=12)
			  {
			idic=1;
			
			fillh("hzvseldic",evento.vpartlab[index[0]],evento.z[index[0]],1);
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

			if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzndic",(float)ival[index[0]],1);
			  }
		
			fillh("hdiff",evento.z[index[0]],tke,1);
			fillh("hwilc",evento.thetacm[index[0]],tke,1);
			fillh("hz1z2dic",evento.z[index[0]],evento.z[index[1]],1);
			  }//idic1

			if(evento.z[index[0]]>18)
			  {
			idic18=1;
			
			fillh("hzvseldic18",evento.vpartlab[index[0]],evento.z[index[0]],1);
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

			if(evento.rcoqf[index[0]]==1)
			  {
			fillh("hzndic18",(float)ival[index[0]],1);
			  }
		
			fillh("hdiff18",evento.z[index[0]],tke,1);
			fillh("hwilc18",evento.thetacm[index[0]],tke,1);
			fillh("hz1z2dic18",evento.z[index[0]],evento.z[index[1]],1);
			  }//idic18



		      }
		    if(thetarel>=40 && thetarel<=100)
		      {
			float vcmcoppialab[3];
			float vcmcoppiacm[3];
			for(int k=0;k<3;k++)
			  {
			    vcmcoppialab[k]=(evento.a[index[0]]*vplab[index[0]][k]+evento.a[index[1]]*vplab[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
			    vcmcoppiacm[k]=(evento.a[index[0]]*vpcm[index[0]][k]+evento.a[index[1]]*vpcm[index[1]][k])/(evento.a[index[0]]+evento.a[index[1]]);
			    
			  }
//************** ifis 1
			if((evento.z[index[0]]+evento.z[index[1]])>=12)
			  {
			ifis=1;
			if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1)
			  {
			    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
			    fillh("hznfis",(float)ivalc,1);
			  }
			
			fillh("hz1z2fis",evento.z[index[0]],evento.z[index[1]],1);
	
			fillh("hzvfiscoppia",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);
			if(vcmcoppiacm[2]>=0)
			  {
			    fillh("hzvfiscoppiavzpos",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);
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
			      }

			  }
			  }//ifis
 //************** ifis 18
			if((evento.z[index[0]]+evento.z[index[1]])>18)
			  {
			ifis18=1;
			if(evento.rcoqf[index[0]]==1&&evento.rcoqf[index[1]]==1)
			  {
			    int ivalc=(evento.z[index[0]]+evento.z[index[1]])*100+(evento.a[index[0]]+evento.a[index[1]])-(evento.z[index[0]]+evento.z[index[1]]);
			    fillh("hznfis18",(float)ivalc,1);
			  }
			
			fillh("hz1z2fis18",evento.z[index[0]],evento.z[index[1]],1);
			
			fillh("hzvfiscoppia18",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);
			if(vcmcoppiacm[2]>=0)
			  {
			    fillh("hzvfiscoppiavzpos18",Classe_formule::modulo(vcmcoppialab),evento.z[index[0]]+evento.z[index[1]],1);
			  }
			  }//ifis18

		      }


		   }

	   }//thetarel>=0
	 if(evento.rcoqf[index[0]]==1)
	   {
	     fillh("hznbig",(float)ival[index[0]],1);
	   }
       }
     for(int j=0;j<evento.moltepl;j++)
       {
	 fillh("hzvsel1",evento.vpartlab[j],evento.z[j],1);
       }




   }//tflow>=8 && tflow<=30 && sumz>=12

 if(tflow>=50 && tflow<=90 && sumz>=12)
   {
     icent=1;
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


 //cout<<iqp<<" "<<idic<<" "<<ifis<<endl;
 fillh("hcond",(float)iqp,(float)idic,1);
 fillh("hcond2",(float)iqp,(float)ifis,1);

 if(idic==1)
   {
     fillh("hneventi",18.,1);
     fillh("hmfragdic",(float)nfrag,1);
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


   }//idic


 if(idic18==1)
   {
     fillh("hneventi",17.,1);
     fillh("hmfragdic18",(float)nfrag,1);
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


   }//idic18





 if(ifis==1)
   {
     fillh("hneventi",19.,1);
   }


 if(ifis18==1)
   {
     fillh("hneventi",24.,1);
   }
 if(iqp==1)
   {
     fillh("hneventi",20.,1);
   }
 if(iqp18==1)
   {
     fillh("hneventi",22.,1);
   }

 if(icent==1)
   {
     fillh("hneventi",21.,1);
   }
 if(icent18==1)
   {
     fillh("hneventi",23.,1);
   }
 float nsuzlcp=0;
 float nt=0;
 float zt=0;
 float ntneg=0;
 float ztneg=0;
 if(iqp==1 || idic==1)
   {
     for(int j=0;j<evento.moltepl;j++)
       {

	 Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     if(((idic==1 && j!=index[1])||(iqp==1)) && j!=index[0]&&evento.rcoqf[j]==1)
	       //	 if(j!=index[0] && evento.rcoqf[j]==1)
	   {

	     fillh("hnzconqp",(float)ival[j],1);
	     fillh("hnzconqpv",(float)ival[j],vpar,1);
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
     for(int j=0;j<evento.moltepl;j++)
       {

	 Classe_formule::vparvperp(vpcm[j],vpcm[index[0]],&vpar,&vperp);
	     if(((idic18==1 && j!=index[1])||(iqp18==1)) && j!=index[0]&&evento.rcoqf[j]==1)
	       //	 if(j!=index[0] && evento.rcoqf[j]==1)
	   {

	     fillh("hnzconqp18",(float)ival[j],1);
	     fillh("hnzconqp18v",(float)ival[j],vpar,1);
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


}
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void nsuzv(TH2F *h,TGraphErrors *gnz[100]);
void Classe_analisi::RoutineFinale()
{
  cout<<"entro nella routine finale"<<endl;
TGraphErrors *gm[100];
TGraphErrors *gznvimfdic[100];
TGraphErrors *gznvimf18dic[100];
TGraphErrors *gnzconqpv[100];
TGraphErrors *gnzconqp18v[100];
  for(int j=0;j<100;j++) 
    {
      gm[j]=0;
      gznvimfdic[j]=0;
      gnzconqpv[j]=0;
      gnzconqp18v[j]=0;
    }

  for(int j=1;j<26;j++) 
    { 
      gznvimfdic[j]=getg(Form("gznvimfdic%d",j));
      gznvimf18dic[j]=getg(Form("gznvimf18dic%d",j));
gnzconqpv[j]=getg(Form("gnzconqpv%d",j));
gnzconqp18v[j]=getg(Form("gnzconqp18v%d",j));
    }
nsuz(geth1("hznbig"),getg("gnzbig"),gm);
nsuz(geth1("hzngrossone"),getg("gnzgrossone"),gm);
nsuz(geth1("hzngrossone18"),getg("gnzgrossone18"),gm);
nsuz(geth1("hnzimfdic"),getg("gnzimfdic"),gm);
nsuz(geth1("hzndic"),getg("gnzdic"),gm);
nsuz(geth1("hzndic18"),getg("gnzdic18"),gm);
nsuz(geth1("hzndic0"),getg("gnzdic0"),gm);
nsuz(geth1("hznfis"),getg("gnzfis"),gm);
nsuz(geth1("hznfis18"),getg("gnzfis18"),gm);
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


nsuzv(geth2("hznvimfdic"),gznvimfdic);
nsuzv(geth2("hznvimfdic18"),gznvimf18dic);
nsuzv(geth2("hnzconqpv"),gnzconqpv);
nsuzv(geth2("hnzconqp18v"),gnzconqp18v);

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


 for(int iz=3;iz<=100;iz++)
    {
      if(totz[iz]>10)
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


 for(int iz=3;iz<=25;iz++)
    {
      for(int iv=0;iv<nbx;iv++)
	{
	  float velo=h->GetYaxis()->GetBinCenter(iv+1);
      if(totz[iz][iv]>10)
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






#endif
