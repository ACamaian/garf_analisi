#ifndef FORMULE
#define FORMULE
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <TMath.h>
#include "Classe_analisi.h"
//static struct
//{
// float pvettore[200][3];
//}pvec_;
//#ifndef pvec
//#define pvec
//extern "C"
//{struct
//{
// float pvettore[200][3];
//}perflowpvec_;
//}
//#endif
#ifndef _ODIE_
extern struct perflow
{
 float pvettore[500][3];
}perflowpvec_;

// struct vettore
//{float pvettore[500][3];
//};
//extern struct vettore *pvec_;
 extern "C"{
   void thetaflow_(float *tflow,int *npart);
 }

extern "C"{
void relkin_(float *apr,float *atar,float *aein, float *athe, float *ae1, float *athe2,float *ae2);

}
 #endif

class Classe_formule
{
 public:
 static constexpr float  cluce=300.;
  static constexpr float amu=931.5;
  static constexpr float me=0.511;

  Classe_formule()
    {

    }

  static float  vcm_classica(float vp,int ap, int at)
    {
      float vcm=ap*vp/(at+ap);
      return vcm;
    }
  static float  vplab_classica(float ebeam)
    {
      float m=1.;
      float vplab=e2v(ebeam,m);

      return vplab;
    }
  static float e2v(float e,float m)
    {

      float v=cluce*TMath::Sqrt(2*e/(amu*m));
      return v;
    }
  static float v2e(float v,float m)
    {
      float e=0.5*m*amu*v*v/(cluce*cluce);

      return e;
    }

  static void polcar(float r, float the, float phi,float *out)
      {
	//printf("r=%f the=%f phi=%f\n",r,the,phi);
	//	if(the==0.)
	//  {
	//   out[2]=0.;
	//  }
	//	else
	//  {
	out[2]=r*cos(the/57.296);
	//  }
	out[0]=r*sin(the/57.296)*cos(phi/57.296);
	out[1]=r*sin(the/57.296)*sin(phi/57.296);


    }
  static float modulo(float *input)
    {
      float out=TMath::Sqrt(input[0]*input[0]+input[1]*input[1]+input[2]*input[2]);
      return out;
    }
static float modulo(vector<float> input)    {
      float out=TMath::Sqrt(input[0]*input[0]+input[1]*input[1]+input[2]*input[2]);
      return out;
    }

  static void carpol(float *r,float *the, float *phi,float *input)
    {
      *r=modulo(input);
    
      float comp=input[2]/(*r);
      if(comp>1){comp=1.;}
      if(comp<-1){comp=-1.;}
      (*the)=57.296*TMath::ACos(comp);

 /*  if(input[0]==0. && input[1]==0.) */
/* 	{ */
/* 	  (*phi)=0.; */
/* 	} */
/*       else */
/* 	{ */
	  (*phi)=57.296*TMath::ATan2(input[1],input[0]);
/* 	} */

  
    }
  static void cm2lab_classica(float *velocm,float *velolab, float vcm)
    {
      for(int k=0;k<3;k++)
	{
	  velolab[k]=velocm[k];
	  if(k==2){velolab[k]=velolab[k]+vcm;}
	}
    }
  static void lab2cm_classica(float *velocm,float *velolab, float vcm)
    {
      for(int k=0;k<3;k++)
	{
	  velocm[k]=velolab[k];
	  if(k==2){velocm[k]=velocm[k]-vcm;}
	}
    }

  static void vecprod(float *a,float *b,float *vout)
    {
      vout[0]=a[1]*b[2]-a[2]*b[1];
      vout[1]=-a[0]*b[2]+a[2]*b[0];
      vout[2]=a[0]*b[1]-a[1]*b[0];

    }
  static void scaprod(float *v1, float *v2, float *vout)
    {
      *vout=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

    }
  static float EAL(float z)
    {
      float a=2.072*z + 2.32E-03 * z*z;
      return a;
    }
	
	static int QualeA(int iz) {
		if(iz<0) return -1;
		if(iz>200) return -1;
		switch(iz) {
			case 0: return 1;
			case 1: return 1;
			case 2: return 4;
			case 3: return 7;
			case 4: return 9;
			case 5: return 11;
			case 6: return 13;
			case 7: return 15;
			case 8: return 16;
			default: return (int)(0.5+2.072*(double)iz + 2.32E-03 * (double)(iz*iz)); //da 9 in su => EAL arrotondata
		}
	}
	
	static double QualeA(double z) {
		return (double)QualeA((int)(z+0.5));
	}
	
	static float QualeA(float z) {
		return (float)QualeA((int)(z+0.5));
	}
	
	static void vparvperp(float *v,float *asse,float *vpar,float *vperp)
    {
      float cosal;
      float scp;
      scaprod(v,asse,&scp);
      cosal=scp/(modulo(v)*modulo(asse));
      if(cosal>1){cosal=1;}
      *vpar=modulo(v)*cosal;      
      float sinal=TMath::Sqrt(1-cosal*cosal);
      *vperp=modulo(v)*sinal;
      return;
      
    }

	//si corregge per il rinculo correlato
	//vpar e vperp sono rispetto all'emettitore (origine e direzione, corretti per il rinculo), vparcm e vperpcm sono con origine nel cm ma con asse orientata lungo l'emettitore (con direzione corretta per il rinculo)
	static void vparvperprecoil(float *vp,float *vsource,float apart,float asource,float *vpar,float *vperp,float *vparcm,float *vperpcm,float *vrec)
    {
      float asse[3];
      float v[3];
      for(int k=0;k<3;k++)
	{
	  asse[k]=(vp[k]*apart+vsource[k]*asource)/(apart+asource);
	  v[k]=vp[k]-asse[k];
	}
      *vrec=modulo(asse);
      float cosal;
      float scp;
      scaprod(v,asse,&scp);
      cosal=scp/(modulo(v)*modulo(asse));
      if(cosal>1){cosal=1;}
      *vpar=modulo(v)*cosal;      
      float sinal=TMath::Sqrt(1-cosal*cosal);
      *vperp=modulo(v)*sinal;

      //cout<<sinal<<" sinal 1 "<<scp<<" "<<modulo(v)*sinal<<endl;

      scaprod(vp,asse,&scp);
      cosal=scp/(modulo(vp)*modulo(asse));
      if(cosal>1){cosal=1;}
      *vparcm=modulo(vp)*cosal;      
      sinal=TMath::Sqrt(1-cosal*cosal);
      *vperpcm=modulo(vp)*sinal;
      //cout<<*vperpcm<<" "<<*vperp<<" "<<modulo(vp)<<" "<<modulo(v)<<endl;

      //cout<<sinal<<" sinal 2 "<<scp<<" "<<modulo(vp)*sinal<<endl;

      return;
      
    }

static  float doppler(float ein,float velo,float theta)
  {float beta=velo/cluce;
    float eout=ein*(1-beta*TMath::Cos(theta/57.296))/(TMath::Sqrt(1-beta*beta));
    return eout;
  }
 static void p2v_cart(float *p,float a,float *v)
   {
     for(int j=0;j<3;j++)
       {
	 v[j]=cluce*p[j]/(a*amu);	 
       }
     return;
   }
 static float p2e(float *p,float a)
   {
     float cp=modulo(p);
     float e=cp*cp/(2*a*931.5);
  return e;
   }
 static float calctkel(float vqp)
   {
     float atot=Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at;
     float vrel=vqp*atot/((float)Classe_analisi::Getanalisi()->reazione.at);
        float tkel=Classe_analisi::Getanalisi()->reazione.Ecm-0.5*Classe_analisi::Getanalisi()->reazione.mi*amu*vrel*vrel/(cluce*cluce);
     return tkel;
   }
 static float v_ela(float theta)
 {
   float vscat;
   
    float vscat1=(pow(Classe_analisi::Getanalisi()->reazione.ap,2)*cos(theta/57.296)*Classe_analisi::Getanalisi()->reazione.vplab+sqrt(pow(Classe_analisi::Getanalisi()->reazione.ap,4)*pow(Classe_analisi::Getanalisi()->reazione.vplab,2)*pow(cos(theta/57.296),2)-(pow(Classe_analisi::Getanalisi()->reazione.ap,4)-pow(Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.at,2))*pow(Classe_analisi::Getanalisi()->reazione.vplab,2)))/(pow(Classe_analisi::Getanalisi()->reazione.ap,2)+Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.at); 
 float vscat2=(pow(Classe_analisi::Getanalisi()->reazione.ap,2)*cos(theta/57.296)*Classe_analisi::Getanalisi()->reazione.vplab-sqrt(pow(Classe_analisi::Getanalisi()->reazione.ap,4)*pow(Classe_analisi::Getanalisi()->reazione.vplab,2)*pow(cos(theta/57.296),2)-(pow(Classe_analisi::Getanalisi()->reazione.ap,4)-pow(Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.at,2))*pow(Classe_analisi::Getanalisi()->reazione.vplab,2)))/(pow(Classe_analisi::Getanalisi()->reazione.ap,2)+Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.at); 
  vscat=vscat1; 
  if(vscat<0){vscat=vscat2;} 
   return vscat;
 }


 static float thetarel(float *v1,float *v2)
 {
   float thetarel;
   float vvrel[3];
   for(int j=0;j<3;j++)
     {
       vvrel[j]=v1[j]-v2[j];
     }
   thetarel=57.296*TMath::ACos((pow(modulo(v1),2)+pow(modulo(v2),2)-pow(modulo(vvrel),2))/(2*modulo(v1)*modulo(v2)));

   return thetarel;
 }

  static float thetarel(vector<float> v1, vector<float> v2)
 {
   float thetarel;
   float vv1[3],vv2[3],vvrel[3];
   for(int j=0;j<3;j++)
     {
     	vv1[j]=v1.at(j);
	vv2[j]=v2.at(j);
       vvrel[j]=vv1[j]-vv2[j];
     }
   thetarel=57.296*TMath::ACos((pow(modulo(vv1),2)+pow(modulo(vv2),2)-pow(modulo(vvrel),2))/(2*modulo(vv1)*modulo(vv2)));

   return thetarel;
 }

 static void da_xyz_a1(float vx,float vy,float vz,float *vout)
 {
   vout[0]=vx;
   vout[1]=vy;
   vout[2]=vz;
 }
 static void vrel(float *v1,float *v2, float *vvrel)
 {
   for(int j=0;j<3;j++)
     {
       vvrel[j]=v1[j]-v2[j];
     }
 }

static void fillh(string hname,double X,double peso) {
	//if(fabs(peso)<0.001) return;
	TH1D *h1temp=(TH1D *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
	if(h1temp==0){cout<<"Manca "<<hname<<endl;}
	h1temp->Fill(X,peso);
	return;
}

static void fillh(string hname,double X,double Y,double peso) {
	//if(fabs(peso)<0.001) return;
	TH2F *h2temp=(TH2F *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
	if(h2temp==0){cout<<"Manca "<<hname<<endl;}
	h2temp->Fill(X,Y,peso);
	return;
}

 static TH1F *geth1(string hname)
   {
     TH1F *h1temp=(TH1F *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
     if(h1temp==0){cout<<"Manca "<<hname<<endl;}
     return h1temp;
   }
 
 static TH2F *geth2(string hname)
   {
     TH2F *h1temp=(TH2F *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
     if(h1temp==0){cout<<"Manca "<<hname<<endl;}
     return h1temp;
   }

 static TGraphErrors *getg(string hname)
   {
         TGraphErrors *h1temp=(TGraphErrors *)gROOT->GetListOfSpecials()->FindObject(hname.c_str());
     if(h1temp==0){cout<<"Manca "<<hname<<endl;}
     return h1temp;
   }

 static void sdroutofplane(float *v1,float *ivec,float *jvec, float *kvec)
 {
   float assez[3]={0.,0.,1.};
   float vout[3];
   vecprod(v1,assez,vout);
 
   for(int j=0;j<3;j++)
     {
       kvec[j]=vout[j]/modulo(vout);
     }
 
   vecprod(kvec,assez,jvec);
 
   vecprod(jvec,kvec,ivec);


 }
 static void outofplane(float *v1,float *kvec,float *the)
   {
     float vout;
     scaprod(v1,kvec,&vout);
     *the=57.296*TMath::ACos(vout/(modulo(v1)*modulo(kvec)));
   }

 static  float Alpha(float *vp1,float *vp2, float a1,float a2) //1 e' il piu' grosso
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


static float viola(int Z1,int Z2,int A1, int A2)
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



#ifndef _ODIE_
static void theta_flow(float *tflow,int *npart)
 {
   thetaflow_(tflow,npart);
   return;
 }

static  void theflow(float *tflow,int npart,float massvec[500],float vcmvec[500][3])
 {

   for(int j=0;j<npart;j++)
     {
       float modvcmvec=sqrt(pow(vcmvec[j][0],2)+pow(vcmvec[j][1],2)+pow(vcmvec[j][2],2));

       float gamma=TMath::Sqrt(1/(1-pow(modvcmvec,2)/(cluce*cluce)));
       //cout<<gamma<<endl;
       for(int k=0;k<3;k++)
	 {

 perflowpvec_.pvettore[j][k]=massvec[j]*amu*gamma*vcmvec[j][k]/cluce;
 //cout<<npart<<" "<<j<<" "<<k<<" "<<vcmvec[j][k]<<" "<<perflowpvec_.pvettore[j][k]<<endl;
	   
	 }
     }

   theta_flow(tflow,&npart);
 }
 

 static void rel_kin(float *apr, float *atar,float *aein, float *athe, float *ae1, float *athe2, float * ae2)
   {
   
     relkin_(apr,atar,aein,athe,ae1,athe2,ae2);
   } 
 #endif
};
#endif
