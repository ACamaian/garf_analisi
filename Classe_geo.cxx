#include "Classe_geo.h"
#include "Classe_formule.h"
#include "Classe_analisi.h"
#include <TRandom3.h>
#include <TNamed.h>
#include <TF1.h>
#include <sys/file.h>
#include <stdio.h>

//energia persa minima in almeno un rivelatore per scrivere la riga in tabella
#define ETHR 0.01
//energia iniziale minima
#define EMIN 0.1
//passo in energia per nucleone (attenzione: influenza il numero di energie in tabella e quindi la sua dimensione)
#define ESTPN 0.05
//passo in energia MAX (attenzione: influenza il numero di energie in tabella e quindi la sua dimensione)
#define ESTP 0.2
//dimensione massima del vettore di energie (massimo numero di energie), non superare mai 10000.
#define NMAX 2000

#define BARBUI 2
#define SCHWALM 4
#define VEDALOSS 5

#define TABGARF VEDALOSS
//#define TABGARF SCHWALM
//#define TABGARF BARBUI
#define TABRCO BARBUI

#define GGASRES 0.08
#define GCSIRES 0.05
#define RGASRES 0.05
#define RSIRES 0.004
#define RCSIRES 0.04


#define GGASMINRES 0.2
#define GCSIMINRES 0.5
#define RGASMINRES 0.2
#define RSIMINRES 0.03
#define RCSIMINRES 0.5


extern "C" {
 void vedaloss_();
}
extern "C" {
 void ecorr_veda_(float *eingresso,float *zpr,float *apr,float *atar, float *eout,float *elost,int *mate,float *thick,int *idir, int *icod, float *pressione);
}

extern "C" {
  void de_vedaloss_(float *z,float *a,float *atar,float *de,float *spess,float *e,int *mate,float *pressione);
}



double Classe_geo::E2Luce_csi(double z,double a,double E)//FORMULA MICHELA
{
	//   28/6/2016 disabilitata la linearizzazione!!!!
        //11/11/2016 riabilitata la linearizzazione
	double ei,coeff,esp;
	double p[8]={12.2424988,288.645213,0.420237154,-0.805045025,1.58972328,0.51370864,0.165660342,0.0181196633};
	
	ei=10.5*a;
	coeff=(p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]));
	esp=p[4]-p[5]*exp(-p[6]*z);
	return (E<=ei)?(coeff*pow(E,esp)):(coeff*((E/ei-1.)*esp+1.)*pow(ei,esp));
	//return (p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]))*pow(E,p[4]-p[5]*exp(-p[6]*z));
}

double Classe_geo::Luce2E_csi(double z,double a,double lo) //formula Michela INVERTITA
{
  // 28/6/2016 disabilitata la linearizzazione!!!!
	//11/11/2016 riabilitata la linearizzazione
	double ei,loi,coeff,esp;
	double p[8]={12.2424988,288.645213,0.420237154,-0.805045025,1.58972328,0.51370864,0.165660342,0.0181196633};
	
	ei=10.5*a;
	coeff=(p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]));
	esp=p[4]-p[5]*exp(-p[6]*z);
	loi=(coeff*pow(ei,esp));
	return (lo<=loi)?(pow(lo/coeff,1./esp)):(((lo/(coeff*pow(ei,esp))-1.)/esp+1)*ei);
	//return pow(lo/((p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]))),1./(p[4]-p[5]*exp(-p[6]*z)));
}

void Classe_geo::vedaloss()
{
  vedaloss_();
  return;
}

void Classe_geo::ecorr_veda(float *eingresso,float *zpr,float *apr,float *atar, float *eout,float *elost,int *mate,float *thick,int *idir, int *icod, float *pressione)
{
  ecorr_veda_(eingresso,zpr,apr,atar,eout,elost,mate,thick,idir,icod,pressione);
  return;
}

int Classe_geo::Soglie_Garfield(float zpart,float apart,float epart,float thetapart,int codice,float *estrip1,float *estrip2,float *ecsi, float *luce,float *zout,float *aout,float *eout,int *iginocchio,int code_micro) {
	*iginocchio=0;
	*estrip1=0;
	*estrip2=0;
	*ecsi=0;
	*luce=0;
	*eout=-1;
	*zout=-1;
	*aout=-1;
	int isec=codice/10;
	int icsi=codice%10;
	int iz=(int)(zpart+0.5);
	int ia=(int)(apart+0.5);
	double epass,eprima,alpha,am,phi,calp,spess;
	double s_tar=target.spessore*10/target.mat.r; //µm
	double s_m1=garf.mylar1;      //µm
	double s_d1=garf.dead1;   //µm
	double s_s1=garf.strip1;   //µm
	double s_d2=garf.dead2;    //µm
	double s_s2=garf.strip2;   //µm
	double s_br=garf.braccio;   //µm
	double s_m2=garf.mylar2;     //µm
	double d_d3;
	
	//roba rotta (gia' tolto in InsideGarf)
	//	if(garf.codice[icsi-1][isec-1]==1) {return 1;}
	//nella camera indietro si prendono solo z=1,2
	if(icsi<5 && zpart>2.5) {return 1;}
	//modifica fatta 2-12-11 per conteggio eventi buttati dalla bl rejection
	if(gRandom->Uniform(0.,1.)>garf.bl_rej[icsi-1][isec-1]) {return 1;}
	//i neutroni non si rivelano
	if(iz<=0) {return 1;}
	
	MATE mylar={"mylar",3,{{6,-1,8},{1,-1,10},{8,-1,4}},1395.,-1.,-1.};
	MATE cf4  ={"CF4",2,{{6,-1,1},{9,-1,4}},-1,garf.pressione,300.};
	if(icsi<5) cf4.P=garf.pressione_back;
	
	//MEZZO TARGET
	epass=eloss((double)epart,iz,ia,&(target.mat),fabs(s_tar/(2.*(double)cos(thetapart/57.296))),SCHWALM);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	if(thetapart<90) {alpha=90.-garf.alpha0-thetapart;}
	else {alpha=thetapart-90.-garf.alpha0;}
	calp=1./cos(alpha/57.296);
	
	//MYLAR 1
	epass=eloss(epass,iz,ia,&mylar,s_m1*calp,TABGARF);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	//DEAD LAYER 1
	epass=eloss(epass,iz,ia,&cf4,s_d1*calp,TABGARF);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	//STRIP DOWN
	eprima=epass;
	epass=eloss(epass,iz,ia,&cf4,s_s1*calp,TABGARF);
	*estrip1=(float)(eprima-epass);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	//DEAD LAYER 2
	epass=eloss(epass,iz,ia,&cf4,s_d2*calp,TABGARF);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	if(icsi>4) {am=90.-garf.alpha0-garf.thecsi[icsi-1];}
	else {am=garf.thecsi[icsi-1]-90.-garf.alpha0;}
	phi=alpha-am;
	d_d3=1000.*garf.dist[icsi-1]/cos(phi/57.296)-(s_br+s_m1+s_d1+s_s1+s_d2+s_s2)*calp;
	
	//STRIP UP
	spess=s_s2*calp;
	if(d_d3<0) {
		spess+=d_d3;
		*iginocchio=1;
	}
	eprima=epass;
	epass=eloss(epass,iz,ia,&cf4,spess,TABGARF);
	*estrip2=(float)(eprima-epass);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	
	//DEAD LAYER 3, SE C'È
	if(d_d3>0) {
		epass=eloss(epass,iz,ia,&cf4,d_d3,TABGARF);
		if((epass<=0)||(!isfinite(epass))) {return 1;}
	}
	
	//MYLAR 2
	epass=eloss(epass,iz,ia,&mylar,s_m2/cos(phi/57.296),TABGARF);
	if((epass<=0)||(!isfinite(epass))) {return 1;}
	*ecsi=(float)epass;
	
	//	if(*ecsi>=0&&iz==20&&ia==48 &&icsi==8){cout<<epart<<" "<<300*sqrt(2*epart/(931.5*apart))<<endl;}

	//spalmatura energia CsI
	*ecsi=gRandom->Gaus(*ecsi,*ecsi*GCSIRES+GCSIMINRES); //risoluz cesio 3%
	if(*ecsi<0.001) {*ecsi=0.001;}
	//spalmatura energia strip down
	*estrip1=gRandom->Gaus(*estrip1,*estrip1*GGASRES+GGASMINRES);//risoluz micro 6%
	if(*estrip1<0.001) {*estrip1=0.001;}
	//spalmatura energia strip up
	*estrip2=gRandom->Gaus(*estrip2,*estrip2*GGASRES+GGASMINRES);//risoluz micro 6%
	if(*estrip2<0.001) {*estrip2=0.001;}
	//niente strip indietro
	if(icsi<5) {*estrip1=0; *estrip2=0;}
	
	//TAGLI IN LUCE CSI DOVUTI A SOGLIE
	double aeff=((iz==1)?(1.):((double)apart));
	*luce=(float)E2Luce_csi((double)zpart,aeff,(double)(*ecsi));
	if(iz==1) {
		if(ia==3 && *luce<garf.cut_csi_t[icsi-1][isec-1]) {apart=2.; ia=2;}
		if(ia==2 && *luce<garf.cut_csi_d[icsi-1][isec-1]) {apart=1.; ia=1;}
		if(*luce<garf.cut_csi_p[icsi-1][isec-1]) return 1;
	}
	if(iz==2) {
		if(ia==3 && *luce<garf.cut_csi_3he[icsi-1][isec-1]) {apart=4.; ia=4;}
		if(*luce<garf.cut_csi_a[icsi-1][isec-1]) return 1;
	}
	
	//RISOLUZIONE ISOTOPICA SOLO DA FAST-SLOW
	*zout=zpart;
	switch(iz) {
		case 1: *aout=apart; break;
		case 2: *aout=(apart<3.5)?3.:4.; break;
		//case 2: *aout=4.; break;
		  //		case 3: *aout=7.; break;
		  //case 4: *aout=7.; break;
		default:
		  *aout=Classe_formule::QualeA(zpart);
		  //	if(iz<=12) *aout=2.*zpart;
		  //	else *aout=Classe_formule::EAL(zpart);
	}
	
	//if(iz==8&&epart>25&&epart<26) printf("%f, %f, %f\n",*estrip1,*estrip2,*ecsi);
	
	//UN UNICO CalcoloEGarf PER GARFIELD FORWARD E BACK!
	*eout=CalcoloEGarf(iz,TMath::Nint(*aout),*estrip1,*estrip2,*ecsi,icsi-1);
	
	//if(fabs(epart-*eout)/epart>0.2&&epart/ *aout>1.) {
	//	printf("Z=%2d A=%2d->%2d riv=(%2d,%d), E=%f->%f\n",iz,ia,TMath::Nint(*aout),isec,icsi,epart,*eout); getchar();
	//}
	
	return 0;
}

int Classe_geo::Soglie_RCO(float zpart,float apart,float epart,float thetapart,float phipart,int *codice,float *egas,float *esi,float *ecsi, float *luce,float *aout, float *zout,float *eout, float *qf)
{
	*egas=0;
	*esi=0;
	*ecsi=0;
	*luce=0;
	*eout=-1;
	*zout=-1;
	*aout=-1;
	int isec0=*codice/100;
	int istrip0=(*codice%100)/10;
	if(istrip0<=0) {istrip0=9;}
	int icsi0=*codice%10;
	if(icsi0<=0) {icsi0=7;}
	int iz=(int)(zpart+0.5);
	int ia=(int)(apart+0.5);
	int bit=0,stop=0,working=7;
	double epass,eprima,calp,range;
	double s_tar=target.spessore*10/target.mat.r; //µm
	double s_m1=rco.spesmylar;       //µm
	double s_gas=rco.spesgas;        //µm
	double s_m2=rco.spesmylar_mezzo; //µm
	double s_m3=rco.spesmylar;       //µm
	double s_si;
	if(istrip0<9) s_si=rco.strip_spess[isec0-1][istrip0-1]; //µm;
	else s_si=rco.strip_spess[isec0-1][0]; //µm;
	double s_mcs=rco.spesmylarcsi;   //µm
	
	//i neutroni non si rivelano
	if(iz<=0) {return 0;}
	
	MATE mylar={"mylar",3,{{6,-1,8},{1,-1,10},{8,-1,4}},1395.,-1.,-1.};
	MATE cf4  ={"CF4",2,{{6,-1,1},{9,-1,4}},-1,rco.pressione,300.};
	MATE si   ={"Si",1,{{14,-1,1}},0,-1.,-1.};
	
	calp=1./cos(thetapart/57.296);
	
	//MEZZO TARGET
	epass=eloss((double)epart,iz,ia,&(target.mat),s_tar*calp/2.,SCHWALM);
	if((epass<=0)||(!isfinite(epass))) {return 0;}
	
	//MYLAR 1
	epass=eloss(epass,iz,ia,&mylar,s_m1*calp,TABRCO);
	if((epass<=0)||(!isfinite(epass))) {return 0;}
	
	//GAS 1
	eprima=epass;
	epass=eloss(epass,iz,ia,&cf4,s_gas*calp,TABRCO);
	*egas=eprima-epass;
	bit+=1;
	stop=1;
	
	//MYLAR 2
	if((epass>0)&&isfinite(epass)) {
		epass=eloss(epass,iz,ia,&mylar,s_m2*calp,TABRCO);
	}
	
	//GAS 2
	if((epass>0)&&isfinite(epass)) {
		eprima=epass;
		epass=eloss(epass,iz,ia,&cf4,s_gas*calp,TABRCO);
		*egas+=eprima-epass;
	}
	
	//MYLAR 3
	if((epass>0)&&isfinite(epass)) {
		epass=eloss(epass,iz,ia,&mylar,s_m3*calp,TABRCO);
	}
	
	//SILICIO
	range=0;
	if((epass>0)&&isfinite(epass)) {
		eprima=epass;
		epass=eloss(epass,iz,ia,&si,s_si*calp,VEDALOSS);
		*esi=eprima-epass;
		bit+=2;
		stop=2;
		//Se non buca il Silicio calcolo il range per PSA
		if((epass<=0)||(!isfinite(epass))) {
			range=e2range(eprima,iz,ia,&si,SCHWALM);
		}
		else range=s_si*calp;
	}
	
	//MYLAR CSI
	if((epass>0)&&isfinite(epass)) {
		epass=eloss(epass,iz,ia,&mylar,s_mcs*calp,TABRCO);
		if(isfinite(epass)) {
			*ecsi=epass;
			*luce=E2Luce_csi((double)zpart,(double)apart,(double)(*ecsi));
			bit+=4;
			stop=3;
		}
		else *ecsi=0;
	}
	
	//SOGLIE "SPERIMENTALI" IN ENERGIA, PER ORA MESSE A SPANNE GUARDANDO LE CORRELAZIONI
	if((*luce<60)&&(bit&4)) {
		bit-=4;
		stop=2;
	}
	//	if((*esi<1.2)&&(bit&2)) {
	if((*esi<0.6)&&(bit&2)) { //soglia in canali: 5 chan
		bit-=2;
		if(stop==2) stop=1;
	}
	if((*egas<0.2)&&(bit&1)) {
		bit-=1;
		if(stop==1) stop=0;
	}
	
	*luce=0;
	
	//La bitpattern "bit" indica i rivelatori che hanno parlato
	//La bitpattern "working" indica i rivelatori funzionanti, anche se non hanno parlato!
	
	//costruisco working
	//tolgo il bit del gas
	if(rco.codice_gas[isec0-1]>0) working&=6;
	//tolgo il bit del Si
	if(istrip0<9) {
		if(rco.codice_strip[isec0-1][istrip0-1]>0) working&=5;
	}
	else working&=5;
	//tolgo il bit del CsI
	if(icsi0<7) {
		if(rco.codice_csi[isec0-1][icsi0-1]>0) working&=3;
	}
	else working&=3;
	
	//costruisco bit
	//se la particella arriverebbe in CsI ma il CsI è guasto la butto!
	if((bit&4)&&((working&4)==0)) return 5;
	//se la particella arriverebbe in Si ma il Si è guasto la butto!
	if((bit&2)&&((working&2)==0)&&(stop==2)) return 6;
	bit&=working;
	//se nessuno ha parlato esco!
	if(bit==0) return 0;
	//se ho solo Si posso fare solo PSA, quindi devo avere range>30µm!
	if((bit==2)&&(range<30)) {return 0;}
	//se ho solo CsI posso fare solo fast-slow, quindi devo avere z<3 e superare le soglie (per controllo soglie devo prima spalmare E)!
	if((bit==4)&&(iz>2)) {return 0;}
	//Se la particella arriva in CsI, non la vedo in Si, ed è con z>2: non la posso identificare! La butto.
	if((bit==5)&&(iz>2)) {return 7;}
	
	//costruisco codice
	*codice=isec0*100;
	if(bit&2) *codice+=istrip0*10;
	else *codice+=90;
	if(bit&4) *codice+=icsi0;
	else *codice+=7;
	
	//SPALMATURE IN ENERGIA
	if(bit&1) {
		*egas=gRandom->Gaus(*egas,*egas*RGASRES+RGASMINRES); //risoluzione gas 6%
		if(*egas<0.001) *egas=0.001;
	}
	else {*egas=0;}
	if(bit&2) {
		*esi=gRandom->Gaus(*esi,*esi*RSIRES+RSIMINRES); //risoluzione Si 0.3%
		if(*esi<0.001) *esi=0.001;
	}
	else {*esi=0;}
	if(bit&4) {
		*ecsi=gRandom->Gaus(*ecsi,*ecsi*RCSIRES+RCSIMINRES); //risoluzione CsI 3%
		if(*ecsi<0.001) *ecsi=0.001;
		*luce=E2Luce_csi((double)zpart,(double)apart,(double)(*ecsi));
	}
	else {*ecsi=0;}
	
	//SOGLIE CSI
	//	if((bit&4)&&((bit&2)==0)) {
	  //	if(((iz==1)&&(*luce<rco.soglie_csi_p[isec0-1][icsi0-1]))||((iz==2)&&(*luce<rco.soglie_csi_a[isec0-1][icsi0-1]))) {
	//		return 8;
	//	}
	//	}
	
	//ASSEGNAZIONE Z,A
	if(bit==1) {
		*zout=100;
		*aout=100;
		//ACHTUNG! Tenere conto in analisi che quando c'è solo gas eout ha un valore fisso.
		*eout=10000;
		return (stop==1)?1:4;
	}
	*zout=zpart;
	if(bit&4)//se ha parlato il cesio (cioe' c'e' il bit 4, non necessariamente da solo)
	  {
	    
	    *aout=apart;
	    *qf=1200;
	    if((bit&2) && rco.intervalli==1) //se ha parlato anche il silicio 
	      {
		if(iz<rco.zinf[isec0-1][istrip0-1][icsi0-1] || iz >rco.zsup[isec0-1][istrip0-1][icsi0-1])//se siamo nella zona in cui le masse non sono cliccate
		  { *aout=Classe_formule::QualeA(zpart);
		    *qf=200;
		  }
		else
		  {
		    if(ia <rco.ainf[isec0-1][istrip0-1][icsi0-1][iz-1] || ia>rco.asup[isec0-1][istrip0-1][icsi0-1][iz-1])//se siamo fuori dagli intervallini		
		  {
		    *qf=-1;
		    *zout=-1;
		    *aout=-1;
		    return 0;
		  }
		  }
	      }
	  }
	else {*aout=Classe_formule::QualeA(zpart);}
	
	//13/4/2015 Se EAL dà A maggiore della massa totale iniziale prendo quest'ultima!!
	if(*aout>Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at) *aout=Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at;
	
	//if(iz==8&&((bit&4)==0)) printf("e= %f %f %f\n",*egas,*esi,*ecsi);
	
	*eout=CalcoloERCO(iz,(int)(*aout+0.5),*egas,*esi,*ecsi,bit,*codice);
	
	//if(fabs(epart-*eout)/epart>0.2&&fabs(epart-*eout)>1.&&(iz==1)) {
	//if(iz==8&&((bit&4)==0)) {
		//printf("Z=%2d A=%2d->%2d riv=(%d,%d,%d), bit=%d, codice=%d, E=%f->%f\n",iz,ia,TMath::Nint(*aout),isec0,istrip0,icsi0,bit,*codice,epart,*eout); getchar();
	//}
	
	return stop;
}

void Classe_geo::prova()
{
  printf("pippo\n");
}

Classe_geo *Classe_geo::geo=0;

Classe_geo * Classe_geo::Getgeo()
   {
     if(geo!=0)
       {return geo;}
     geo=new Classe_geo();
     return geo; 
   }

void Classe_geo::Leggigeo(char *filegeo)
{
	double rho_el[86]={-1,-1,535,1848,2460,2260,-1,-1,-1,-1, //H->Ne
		968,1738,2700,2330,1823,1960,-1,-1, //Na->Ar
		856,1550,2985,4507,6110,7140,7470,7874,8900,8908,8920,7140,5904,5323,5727,4819,3120,-1, //K->Kr
		1532,2630,4472,6511,8570,10280,11500,12370,12450,12023,10490,8650,7310,7310,6697,6240,4940,-1,//Rb->Xe
		1879,3510,6146,6689,6640,7010,7264,7353,5244,7901,8219,8551,8795,9066,9321,6570,9841,13310,16650,19250,21020,22610,22650,21090,19300,13534,11850,11340,9780,9196,-1,-1//Cs->Rn
	};
	char nom_el[86][3]={"H","He","Li","Be","B","C","N","O","F","Ne", //H->Ne
		"Na","Mg","Al","Si","P","S","Cl","Ar", //Na->Ar
		"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", //K->Kr
		"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe", //Rb->Xe
		"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"//Cs->Rn
	};
  vedaloss();

  garf.pressione=0;
  garf.pressione_back=0;
  rco.pressione=0;

  FILE *apri=fopen(filegeo,"r");
  if(apri==0)
    {
      printf("File della geometria %s non esistente\n",filegeo);
      return;
    }
  
   if(flock(fileno(apri),LOCK_EX)) perror("flock");
  target.spessore=Classe_analisi::Getanalisi()->reazione.spess_target;
  target.zt=Classe_analisi::Getanalisi()->reazione.zt;
  target.at=Classe_analisi::Getanalisi()->reazione.at;

	sprintf(target.mat.nome,"%s",nom_el[target.zt-1]);
	target.mat.ne=1;
	target.mat.e[0][0]=target.zt;
	target.mat.e[0][1]=target.at;
	target.mat.e[0][2]=1;
	target.mat.r=rho_el[target.zt-1];
	target.mat.P=-1;
	target.mat.T=-1;
  
  if(target.zt==14){target.materiale=1; target.densita=0.233;}//Si
  if(target.zt==28){target.materiale=4; target.densita=0.8902;}//Ni
  if(target.zt==6){target.materiale=6; target.densita=0.190;}//C
  if(target.zt==47){target.materiale=7; target.densita=1.05;}//Ag
  if(target.zt==50){target.materiale=8; target.densita=0.713;}//Sn
  if(target.zt==79){target.materiale=10; target.densita=1.93;}//Au
  if(target.zt==92){target.materiale=11; target.densita=1.895;}//U
  if(target.zt==41){target.materiale=13; target.densita=0.857;}//Nb
  if(target.zt==73){target.materiale=14; target.densita=1.6654;}//Ta
  if(target.zt==23){target.materiale=15; target.densita=0.611;}//V
  if(target.zt==13){target.materiale=18; target.densita=0.26989;}//Al
  if(target.zt==82){target.materiale=19; target.densita=1.135;}//Pb
  if(target.zt==32){target.materiale=22; target.densita=0.5323;}//Ge
  if(target.zt==20){target.materiale=23; target.densita=0.155;}//Ca
  if(target.zt==29){target.materiale=24; target.densita=0.896;}//Cu
  if(target.zt==22){target.materiale=25; target.densita=0.454;}//Ti
  if(target.zt==83){target.materiale=26; target.densita=0.9747;}//Bi
  if(target.zt==12){target.materiale=27; target.densita=0.1738;}//Mg
  if(target.zt==3){target.materiale=28; target.densita=0.0534;}//Li
  if(target.zt==30){target.materiale=29; target.densita=0.7133;}//Zn

 if(target.zt==5){target.materiale=6; target.densita=0.190;}//il B si tratta come C

  if(target.materiale==-1)
    {
      printf("Materiale del target non presente nel data base; si prende Au\n");
      target.materiale=10;
    }

 char riga[1000];
  int zl,al;
  double ddelta;
  char iso[10];
 
  FILE *aprimasse=fopen("awm95.txt","r");
  if(aprimasse==NULL)
    {
      cout<<"File dei difetti di massa non trovato"<<endl;
    }
  else
    {
        while(fscanf(aprimasse,"\n%[^\n]",riga)!=EOF)
    {
      if(riga[0]!='#')
	{
	  sscanf(riga,"%s %d %lg",iso,&zl,&ddelta);
	  sscanf(iso,"%d",&al);
	  //printf("%s %d %f %d\n",iso,z,ddelta,a);
	  int in=al-zl;
	  // D[zl][in]=ddelta*0.001;
	  D[zl][in]=ddelta*0.001;
	  // cout<<"Z="<<zl<<" N="<<in<<" A="<<zl+in<<" "<<D[zl][in]<<endl;
	}
    }
  fclose(aprimasse);
//  double me=0.511;
//   double amu=931.49;
  double mp=D[1][0]+Classe_formule::amu-Classe_formule::me;
  double mn=D[0][1]+Classe_formule::amu;
  for(int iz=0;iz<111;iz++)
    {
      for(int in=0;in<163;in++)
	{
	  if(D[iz][in]>-100000)
	    {
	      B[iz][in]=(double)iz*mp+(double)in*mn+(double)iz*Classe_formule::me-(double)(in+iz)*Classe_formule::amu-D[iz][in];

	      
	    }
	}
    }
  B[0][1]=0;
  B[1][0]=0;
    }

 
  int iflag=0;
  int iphos=0;
  int igarf_f=0;
  int igarf_b=0;
  int ihector=0;
  //<<<<<<< Classe_geo.cxx
  int irco_strip=0;
  int irco_gas=0;
  int irco_csi=0;
  int ifazietto=0;
  int ispefazietto=0;
  int icl=0;
  int icls=0;
  int iclg=0;
	      
  char riv[200],tipo[200],denominaz[200];
  float d,t,f,spp;
  int ip,ih; 
  int isettore; 
	      float thm,thM;
	      int istr;
		  
		  
  while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
    {
      iflag=0;
      if(riga[0]!='#')
	{
	  sscanf(riga,"%s %s",riv,tipo);

	  if(strcmp(riv,"RIV")==0)
	    {iflag=1;
	      iphos=0;
	      igarf_f=0;
	      igarf_b=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;

	      if(strcmp(tipo,"PHOS")==0)
		{
		  iphos=1;
		  ThereIsPhos=1;
		}
	      if(strcmp(tipo,"GARFIELD_FORWARD")==0)
		{
		  igarf_f=1;
		  //<<<<<<< Classe_geo.cxx
		  ThereIsGarf=1;
		}
	      if(strcmp(tipo,"GARFIELD_BACKWARD")==0)
		{
		  igarf_b=1;
		  ThereIsGarf=1;
		}
	      if(strcmp(tipo,"RCO_STRIP")==0)
		{
		  irco_strip=1;
		  ThereIsRCO=1;
		}
	      if(strcmp(tipo,"RCO_GAS")==0)
		{
		  irco_gas=1;
		  ThereIsRCO=1;
		}
	      if(strcmp(tipo,"RCO_CSI")==0)
		{
		  irco_csi=1;
		  ThereIsRCO=1;

		}
	      //<<<<<<< Classe_geo.cxx


	      if(strcmp(tipo,"HECTOR")==0)
		{
		  ihector=1;
		  //<<<<<<< Classe_geo.cxx
		  ThereIsHector=1;
		}
	      if(strcmp(tipo,"FAZIETTINO")==0)
		{
		  ifazietto=1;
		  ThereIsBlocchiFazia=1;
	
		}
	    }
	    
	    
	  //	  printf("%s\n",tipo);
	  if(iphos==1 && iflag==0)
	    {
	      sscanf(riga,"%d %d %f %f %f %f",&ip,&ih,&d,&t,&f,&spp);
	      phos.dist[ip-1][ih-1]=d;
	      phos.theta[ip-1][ih-1]=t;
	      phos.phi[ip-1][ih-1]=f;
	      phos.spess[ip-1][ih-1]=spp;
	      	Classe_formule::polcar(phos.dist[ip-1][ih-1],phos.theta[ip-1][ih-1],phos.phi[ip-1][ih-1],phos.veccentro[ip-1][ih-1]);
float vkm=Classe_formule::modulo(phos.veccentro[ip-1][ih-1]);
              for(int nn=0;nn<3;nn++)
		{
		  phos.vk[ip-1][ih-1][nn]=phos.veccentro[ip-1][ih-1][nn]/vkm;
		}
	      Classe_formule::vecprod(phos.vk[ip-1][ih-1],assex,phos.vj[ip-1][ih-1]);
	      Classe_formule::vecprod(phos.vj[ip-1][ih-1],phos.vk[ip-1][ih-1],phos.vi[ip-1][ih-1]);

	      float vim=Classe_formule::modulo(phos.vi[ip-1][ih-1]);
	      float vjm=Classe_formule::modulo(phos.vj[ip-1][ih-1]);
	      for(int nn=0;nn<3;nn++)
		{
		  phos.vi[ip-1][ih-1][nn]=phos.vi[ip-1][ih-1][nn]/vim;
		  phos.vj[ip-1][ih-1][nn]=phos.vj[ip-1][ih-1][nn]/vjm;
		}

	    }
	  if(igarf_f==1 && iflag==0)
	    {
	      sscanf(riga,"%d %f %f %f %f",&ip,&d,&t,&thm,&thM);
	      garf.thecsi[ip-1]=t;
	      garf.dist[ip-1]=d;
	      garf.themincsi[ip-1]=thm;
	      garf.themaxcsi[ip-1]=thM;
	      if(garf.gcsi[ip-1]==0)
		{
		  garf.gcsi[ip-1]=new TGraph();
		  //		  		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2,-garf.h[ip-1]/2);
		  // garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l1[ip-1]/2,-garf.h[ip-1]/2);
		  //garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l2[ip-1]/2,garf.h[ip-1]/2);
		  //garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l2[ip-1]/2,garf.h[ip-1]/2);
		  // garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2,-garf.h[ip-1]/2);

		  		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2+garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l1[ip-1]/2-garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l2[ip-1]/2-garf.spess_fasciatura,garf.h[ip-1]/2-garf.hup[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l2[ip-1]/2+garf.spess_fasciatura,garf.h[ip-1]/2-garf.hup[ip-1]);
		   garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2+garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);

	
		}
	      for(int isec=0;isec<24;isec++)
		{
		  //		  if(isec<6){garf.phi[ip-1][isec]=82.5-isec*garf.dphi;}
		  if(isec<6){garf.phi[ip-1][isec]=-82.5+isec*garf.dphi;}
		  if(isec>=6 && isec<=17){garf.phi[ip-1][isec]=7.5+(isec-6)*garf.dphi;}

 if(isec>=18){garf.phi[ip-1][isec]=-172.5+(isec-18)*garf.dphi;}
	      	Classe_formule::polcar(garf.dist[ip-1],garf.thecsi[ip-1],garf.phi[ip-1][isec],garf.veccentro[isec][ip-1]);
float vkm=Classe_formule::modulo(garf.veccentro[isec][ip-1]);

              for(int nn=0;nn<3;nn++)
		{
		  garf.vk[isec][ip-1][nn]=garf.veccentro[isec][ip-1][nn]/vkm;
		}
	      //	      Classe_formule::vecprod(garf.vk[isec][ip-1],vfascio,garf.vi[isec][ip-1]);
	      Classe_formule::vecprod(vfascio,garf.vk[isec][ip-1],garf.vi[isec][ip-1]);

	      //Classe_formule::vecprod(garf.vk[isec][ip-1],garf.vi[isec][ip-1],garf.vj[isec][ip-1]);
	      Classe_formule::vecprod(garf.vi[isec][ip-1],garf.vk[isec][ip-1],garf.vj[isec][ip-1]);

	      float vim=Classe_formule::modulo(garf.vi[isec][ip-1]);
	      float vjm=Classe_formule::modulo(garf.vj[isec][ip-1]);
	      for(int nn=0;nn<3;nn++)
		{
		  garf.vi[isec][ip-1][nn]=garf.vi[isec][ip-1][nn]/vim;
		  garf.vj[isec][ip-1][nn]=garf.vj[isec][ip-1][nn]/vjm;
		}
		}
	    }
	  if(igarf_b==1 && iflag==0)
	    {
	      sscanf(riga,"%d %f %f %f %f",&ip,&d,&t,&thm,&thM);
	      garf.thecsi[ip-1]=t;
	      garf.dist[ip-1]=d;
	      garf.themincsi[ip-1]=thm;
	      garf.themaxcsi[ip-1]=thM;
	      if(garf.gcsi[ip-1]==0)
		{
		  garf.gcsi[ip-1]=new TGraph();
// 		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2,-garf.h[ip-1]/2);
// 		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l1[ip-1]/2,-garf.h[ip-1]/2);
// 		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l2[ip-1]/2,garf.h[ip-1]/2);
// 		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l2[ip-1]/2,garf.h[ip-1]/2);
// 		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2,-garf.h[ip-1]/2);
		  		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2+garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l1[ip-1]/2-garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),garf.l2[ip-1]/2-garf.spess_fasciatura,garf.h[ip-1]/2-garf.hup[ip-1]);
		  garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l2[ip-1]/2+garf.spess_fasciatura,garf.h[ip-1]/2-garf.hup[ip-1]);
		   garf.gcsi[ip-1]->SetPoint(garf.gcsi[ip-1]->GetN(),-garf.l1[ip-1]/2+garf.spess_fasciatura,-garf.h[ip-1]/2+garf.hdown[ip-1]);

		}
	      for(int isec=2;isec<23;isec++)
		{
		  //		  if(isec<6){garf.phi[ip-1][isec]=82.5-isec*garf.dphi;}
		  if(isec<6){garf.phi[ip-1][isec]=-82.5+isec*garf.dphi;}
		  if(isec>=6 && isec<=17){garf.phi[ip-1][isec]=7.5+(isec-6)*garf.dphi;}

 if(isec>=18){garf.phi[ip-1][isec]=-172.5+(isec-18)*garf.dphi;}
	      	Classe_formule::polcar(garf.dist[ip-1],garf.thecsi[ip-1],garf.phi[ip-1][isec],garf.veccentro[isec][ip-1]);
float vkm=Classe_formule::modulo(garf.veccentro[isec][ip-1]);

              for(int nn=0;nn<3;nn++)
		{
		  garf.vk[isec][ip-1][nn]=garf.veccentro[isec][ip-1][nn]/vkm;
		}
	      //	      Classe_formule::vecprod(garf.vk[isec][ip-1],vfascio,garf.vi[isec][ip-1]);
	      Classe_formule::vecprod(vfascio,garf.vk[isec][ip-1],garf.vi[isec][ip-1]);

	      //	      Classe_formule::vecprod(garf.vk[isec][ip-1],garf.vi[isec][ip-1],garf.vj[isec][ip-1]);
	      Classe_formule::vecprod(garf.vi[isec][ip-1],garf.vk[isec][ip-1],garf.vj[isec][ip-1]);


	      float vim=Classe_formule::modulo(garf.vi[isec][ip-1]);
	      float vjm=Classe_formule::modulo(garf.vj[isec][ip-1]);
	      for(int nn=0;nn<3;nn++)
		{
		  garf.vi[isec][ip-1][nn]=garf.vi[isec][ip-1][nn]/vim;
		  garf.vj[isec][ip-1][nn]=garf.vj[isec][ip-1][nn]/vjm;
		}
		}
	    }
	  if(irco_strip==1 && iflag==0)
	    {
	      float spestr;
	      sscanf(riga,"%d %d %f %f %f %f",&isettore,&istr,&d,&thm,&thM,&spestr);
	      if(icls==0)
		{
	      rco.strip_veccentro[0]=0.;
	      rco.strip_veccentro[1]=0.;
	      rco.strip_veccentro[2]=d;
	      icls=1;
		}
		  rco.strip_dist[isettore-1][istr-1]=d;
		  rco.strip_themin[isettore-1][istr-1]=thm;
		  rco.strip_themax[isettore-1][istr-1]=thM;
		  rco.strip_the[isettore-1][istr-1]=(thm+thM)/2;
		  rco.strip_phimin[isettore-1][istr-1]=(isettore-1)*rco.dphi-rco.dphi/2;
		  rco.strip_phimax[isettore-1][istr-1]=rco.strip_phimin[isettore-1][istr-1]+rco.dphi;
		  rco.strip_phi[isettore-1][istr-1]=rco.strip_phimin[isettore-1][istr-1]+rco.dphi/2;
		  while(rco.strip_phimin[isettore-1][istr-1]>=180)
		    {rco.strip_phimin[isettore-1][istr-1]=rco.strip_phimin[isettore-1][istr-1]-360;}
		  while(rco.strip_phimax[isettore-1][istr-1]>180)
		    {rco.strip_phimax[isettore-1][istr-1]=rco.strip_phimax[isettore-1][istr-1]-360;}
		  while(rco.strip_phi[isettore-1][istr-1]>180)
		    {rco.strip_phi[isettore-1][istr-1]=rco.strip_phi[isettore-1][istr-1]-360;}

		  rco.strip_spess[isettore-1][istr-1]=spestr;
	    
	    }
	  if(irco_gas==1 && iflag==0)
	    {
	      sscanf(riga,"%d %f %f %f",&isettore,&d,&thm,&thM);
	      if(iclg==0)
		{
	      rco.gas_veccentro[0]=0.;
	      rco.gas_veccentro[1]=0.;
	      rco.gas_veccentro[2]=d;
	      iclg=1;
		}
	      rco.gas_dist[isettore-1]=d;
	      rco.gas_themin[isettore-1]=thm;
	      rco.gas_themax[isettore-1]=thM;
	      rco.gas_the[isettore-1]=(thm+thM)/2;
	      rco.gas_phimin[isettore-1]=(isettore-1)*rco.dphi-rco.dphi/2;
	      rco.gas_phimax[isettore-1]=rco.gas_phimin[isettore-1]+rco.dphi;
	      rco.gas_phi[isettore-1]=rco.gas_phimin[isettore-1]+rco.dphi/2;
	      	  while(rco.gas_phimin[isettore-1]>=180)
		    {
		  rco.gas_phimin[isettore-1]=rco.gas_phimin[isettore-1]-360;
		    }
	      while(rco.gas_phimax[isettore-1]>180)
		{
		  rco.gas_phimax[isettore-1]=rco.gas_phimax[isettore-1]-360;
		}
	      while(rco.gas_phi[isettore-1]>180)
		{
		  rco.gas_phi[isettore-1]=rco.gas_phi[isettore-1]-360;
		}
	      

	      rco.gas_vk[isettore-1][0]=0.;
	      rco.gas_vk[isettore-1][1]=0.;
	      rco.gas_vk[isettore-1][2]=1.;
	      rco.gas_vi[isettore-1][0]=TMath::Cos(rco.gas_phi[isettore-1]/57.296);
	      rco.gas_vi[isettore-1][1]=-TMath::Sin(rco.gas_phi[isettore-1]/57.296);
	      rco.gas_vi[isettore-1][2]=0;
	      rco.gas_vj[isettore-1][0]=TMath::Sin(rco.gas_phi[isettore-1]/57.296);
	      rco.gas_vj[isettore-1][1]=TMath::Cos(rco.gas_phi[isettore-1]/57.296);
	      rco.gas_vj[isettore-1][2]=0;

	     

	    }
	  if(irco_csi==1 && iflag==0)
	    {
	      float th;

	      sscanf(riga,"%d %f %f %f %f",&ip,&d,&th,&thm,&thM);
	      if(icl==0)
		{
	      rco.csi_veccentro[0]=0.;
	      rco.csi_veccentro[1]=0.;
	      rco.csi_veccentro[2]=d;
	      icl=1;
		}
	      for(int k=0;k<8;k++)
		{
		  rco.csi_dist[k][ip-1]=d;
		  rco.csi_the[k][ip-1]=th;
		  rco.csi_themin[k][ip-1]=thm;
		  rco.csi_themax[k][ip-1]=thM;
		  if(ip%2==1)//dispari
		    {
		      rco.csi_phimin[k][ip-1]=k*rco.dphi-rco.dphi/2;
		      rco.csi_phimax[k][ip-1]=rco.csi_phimin[k][ip-1]+rco.dphi/2;
		      rco.csi_phi[k][ip-1]=rco.csi_phimin[k][ip-1]+rco.dphi/4;
		    }
		  else
		    {
		    rco.csi_phimin[k][ip-1]=k*rco.dphi;
		    rco.csi_phimax[k][ip-1]=rco.csi_phimin[k][ip-1]+rco.dphi/2;
		    rco.csi_phi[k][ip-1]=rco.csi_phimin[k][ip-1]+rco.dphi/4;
		    }
		  		  
		  while(rco.csi_phimin[k][ip-1]>=180)
		    {rco.csi_phimin[k][ip-1]=rco.csi_phimin[k][ip-1]-360;}
		  while(rco.csi_phimax[k][ip-1]>180)
		    {rco.csi_phimax[k][ip-1]=rco.csi_phimax[k][ip-1]-360;}
		  while(rco.csi_phi[k][ip-1]>180)
		    {rco.csi_phi[k][ip-1]=rco.csi_phi[k][ip-1]-360;}

		  //cout<<"set="<<k<<" csi="<<ip-1<<" "<<rco.csi_phimin[k][ip-1]<<" "<<rco.csi_phimax[k][ip-1]<<endl;


		}
	    }
	  if(ihector==1 && iflag==0)
	    {
	      sscanf(riga,"%d %f %f %f",&ip,&d,&t,&f);
	      hector.theta[ip-1]=t;
	      hector.phi[ip-1]=f;
	      hector.dist[ip-1]=d;
	      Classe_formule::polcar(hector.dist[ip-1],hector.theta[ip-1],hector.phi[ip-1],hector.veccentro[ip-1]);	      
float vkm=Classe_formule::modulo(hector.veccentro[ip-1]);
              for(int nn=0;nn<3;nn++)
		{
		  hector.vk[ip-1][nn]=hector.veccentro[ip-1][nn]/vkm;
		}
	      Classe_formule::vecprod(hector.vk[ip-1],assex,hector.vj[ip-1]);
	      Classe_formule::vecprod(hector.vj[ip-1],hector.vk[ip-1],hector.vi[ip-1]);
	      float vim=Classe_formule::modulo(hector.vi[ip-1]);
	      float vjm=Classe_formule::modulo(hector.vj[ip-1]);
	      for(int nn=0;nn<3;nn++)
		{
		  hector.vi[ip-1][nn]=hector.vi[ip-1][nn]/vim;
		  hector.vj[ip-1][nn]=hector.vj[ip-1][nn]/vjm;
		}

	     	    }
	  if(ifazietto==1 && iflag==0)
	    {
	      sscanf(riga,"%d %f %f %f",&ip,&d,&t,&f);//si legge il numero del blocco, la distanza del centro, theta e phi
	      fazietto.dist[ip]=d;
	      fazietto.theta[ip]=t;
	      fazietto.phi[ip]=f;

	      Classe_formule::polcar(d,t,f,fazietto.vcentroblocco[ip]);
	      float vkm=Classe_formule::modulo(fazietto.vcentroblocco[ip]);
	      
	      for(int nn=0;nn<3;nn++)
		{
		  fazietto.vkblocco[ip][nn]=fazietto.vcentroblocco[ip][nn]/vkm;
		}
	      //	      Classe_formule::vecprod(fazietto.vkblocco[ip],assex,fazietto.vjblocco[ip]);
	      // Classe_formule::vecprod(fazietto.vjblocco[ip],fazietto.vkblocco[ip],fazietto.viblocco[ip]);
	      Classe_formule::vecprod(assey,fazietto.vkblocco[ip],fazietto.vjblocco[ip]);
	      Classe_formule::vecprod(fazietto.vkblocco[ip],fazietto.vjblocco[ip],fazietto.viblocco[ip]);
	      //coordinate sul piano del blocco dei centri dei 16 rivelatori del blocco:

	         //Quadrifoglio 4
	      fazietto.vx[ip][15]=-fazietto.dx-fazietto.dx/2;
	      fazietto.vy[ip][15]=-fazietto.dy-fazietto.dy/2;
	      fazietto.vx[ip][12]=-fazietto.dx-fazietto.dx/2;
	      fazietto.vy[ip][12]=-fazietto.dy/2;
	      fazietto.vx[ip][13]=-fazietto.dx/2;
	      fazietto.vy[ip][13]=-fazietto.dy/2;
	      fazietto.vx[ip][14]=-fazietto.dx/2;
	      fazietto.vy[ip][14]=-fazietto.dy-fazietto.dy/2;
	      //Quadrifoglio 0
	      fazietto.vx[ip][3]=-fazietto.dx-fazietto.dx/2;
	      fazietto.vy[ip][3]=fazietto.dy/2;
	      fazietto.vx[ip][0]=-fazietto.dx-fazietto.dx/2;
	      fazietto.vy[ip][0]=fazietto.dy+fazietto.dy/2;
	      fazietto.vx[ip][1]=-fazietto.dx/2;
	      fazietto.vy[ip][1]=fazietto.dy+fazietto.dy/2;
	      fazietto.vx[ip][2]=-fazietto.dx/2;
	      fazietto.vy[ip][2]=fazietto.dy/2;
	      //Quadrifoglio 2
	      fazietto.vx[ip][7]=fazietto.dx/2;
	      fazietto.vy[ip][7]=fazietto.dy/2;
	      fazietto.vx[ip][4]=fazietto.dx/2;
	      fazietto.vy[ip][4]=fazietto.dy+fazietto.dy/2;
	      fazietto.vx[ip][5]=fazietto.dx+fazietto.dx/2;
	      fazietto.vy[ip][5]=fazietto.dy+fazietto.dy/2;
	      fazietto.vx[ip][6]=fazietto.dx+fazietto.dx/2;
	      fazietto.vy[ip][6]=fazietto.dy/2;
	      //Quadrifoglio 3
	      fazietto.vx[ip][11]=fazietto.dx/2;
	      fazietto.vy[ip][11]=-fazietto.dy-fazietto.dy/2;
	      fazietto.vx[ip][8]=fazietto.dx/2;
	      fazietto.vy[ip][8]=-fazietto.dy/2;
	      fazietto.vx[ip][9]=fazietto.dx+fazietto.dx/2;
	      fazietto.vy[ip][9]=-fazietto.dy/2;
	      fazietto.vx[ip][10]=fazietto.dx+fazietto.dx/2;
	      fazietto.vy[ip][10]=-fazietto.dy-fazietto.dy/2;

	     
	    }

	  if(strcmp(riv,"CODICE_PHOS")==0)
	    {
#ifndef _ODIE_
	      TNamed *name1=new TNamed(riga,"s=1000");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name1->Write();
#endif
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=2;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int ip,ih,ic1,ic2,ic3;

	      char riga2[1000];
	      int codex;
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d %d %d",&ip,&ih,&ic1,&ic2,&ic3);
			  
			  if(ic1>0){ic1=1;}
			  if(ic2>0){ic2=1;}
			  if(ic3>0){ic3=1;}
			  codex=ic1+ic2*2+ic3*4;

			      phos.codice[ip-1][ih-1]=codex;
			     
			    
			}
		    }
		    }
		  fclose(ap);
		}

		
	  if(strcmp(riv,"CODICE_GARF")==0)
	    {
#ifndef _ODIE_
	      TNamed *name4=new TNamed(riga,"s=1003");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name4->Write();
#endif
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=5;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int icsi,isec;

	      char riga2[1000];
	      int codex;
	      float code_bl_rej;

	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d %f",&icsi,&isec,&codex,&code_bl_rej);
			  garf.codice[icsi-1][isec-1]=codex;
			  garf.bl_rej[icsi-1][isec-1]=code_bl_rej;
			     
			    
			}
		    }
		    

	      fclose(ap);
		}

		}
	  if(strcmp(riv,"MICRO_ROTTE")==0)//Micro rotte
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=5;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("micro rotte=%s\n",tipo);

	      int icsi,isec;

	      char riga2[1000];
	      int cl1,cl2,cl3,cl4;
	     
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d %d %d %d",&isec,&icsi,&cl1,&cl2,&cl3,&cl4);
			    garf.micro_rotte[isec-1][icsi-1][0]=cl1;//up left	1 buona, 0 rotta		  			     
			  garf.micro_rotte[isec-1][icsi-1][1]=cl2;//up right			  			     
			  garf.micro_rotte[isec-1][icsi-1][2]=cl3;//down left			  			     
			  garf.micro_rotte[isec-1][icsi-1][3]=cl4;//down right			  			     
			   // printf("TROVATE MICRO_ROTTE: s%2d c%d: %d %d %d %d\n",isec,icsi,cl1,cl2,cl3,cl4);
			}
		    }
		    

	      fclose(ap);
          //getchar();
		}

	    }
if(strcmp(riv,"FAZIA_TEL_ROTTI")==0)//tele fazia rotti
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=15;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("tele fazia rotti=%s\n",tipo);

	      int ibl,iqua,itel;

	      char riga2[1000];
	      	      int cl1,cl2,cl3;
	     
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d %d %d %d",&ibl,&iqua,&itel,&cl1,&cl2,&cl3);
			  fazietto.rotti[ibl][iqua-1][itel-1]=cl3+cl2*2+cl1*4;
			  // cout<<ibl<<" "<<iqua<<" "<<itel<<" "<<fazietto.rotti[ibl][iqua-1][itel-1]<<" "<<(fazietto.rotti[ibl][iqua-1][itel-1]&1)<<" "<<(fazietto.rotti[ibl][iqua-1][itel-1]&2)<<" "<<(fazietto.rotti[ibl][iqua-1][itel-1]&4)<<endl;
			}
		    }
		    

	      fclose(ap);
          //getchar();
		}

	    }

if(strcmp(riv,"FAZIA_GRID_Si1Si2")==0)//griglie Si1Si2
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=215;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("Griglie Si1Si2 per fazia=%s\n",tipo);
	      cout<<"Si leggono i limiti di Z per identificazione in massa per Si1 Si2"<<endl;
	      int ibl,iqua,itel;

	      char riga2[1000];
	      //	      int cl1,cl2,cl3,cl4;
	      TString slocal;
	      TString sloc;
	      int num,zval;
	      if(ap!=0)
		{
		  int Ik2=-2;
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      
		      if(riga2[0]!='#')
			{
			  slocal=riga2;
		
			  if(slocal.Contains("<PARAMETER> IDTelescopes"))
			    {
			       Ik2=-2;
			      // cout<<slocal.Data()<<endl;
			  int ik=slocal.Index("SI2_");
			  sscanf(&riga2[ik+4],"%d",&num);
			  
			  if(num>100)
			    {
			      ibl=num/100;
			    }
			  else
			    {
			      ibl=0;
			    }
			  iqua=(num-ibl*100)/10;
			  itel=num-ibl*100-iqua*10;
			  //  cout<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<endl;

			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      zval=-1;
			      if(sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				 
				  if(sloc.Contains("-"))
				    {
				      
				      Ik2=-1;
				      break;
				    }
				  else
				    {
				  Ik2=sloc.Index("PIDRANGE=")+9;
				    
				  sscanf(&riga2[Ik2],"%d",&zval);
				  break;
				    }
				}
			      if(Ik2==-2 && sloc.Contains("ID:Z="))
				{
				  break;
				}
			    }//for(;;)
			  if(Ik2==-2)continue;
			  if(zval>0)
			    {
			  fazietto.alimSi1Si2[ibl][iqua-1][itel-1]=zval;
			  //			  cout<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<" "<<aval<<endl;
			    }
			  if(Ik2>=-1)
			    {
			  int iz,alom,aloM;
			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      
			      if(sloc.Contains("<PARAMETER> PIDRANGE")&&!sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				 
				  int ik2=sloc.Index("PIDRANGE")+8;
				  sscanf(&riga2[ik2],"%d",&iz);
				  int ik3=sloc.Index("=")+1;
				  sscanf(&riga2[ik3],"%d",&alom);
				  if(sloc.Contains("|"))
				    {
				    
				  int ik4=sloc.Last('|')+1;
				  sscanf(&riga2[ik4],"%d",&aloM);
				    
				  fazietto.infmassesi1si2[ibl][iqua-1][itel-1][iz]=alom;
				  fazietto.supmassesi1si2[ibl][iqua-1][itel-1][iz]=aloM;
				    }

				  //				  cout<<ibl<<" "<<iqua<<" "<<itel<<" "<<zval<<" "<<iz<<" "<<alom<<" "<<aloM<<endl;
				  if(iz==zval)
				    {break;}
				}
if((sloc.Contains("+KVIDZALine")||sloc.Contains("+KVIDCutContour"))&&zval==-1)
				{
				  int ii=-1;
				  for(int izz=iz;izz>=0;izz--)
				    {
				      if(fazietto.infmassesi1si2[ibl][iqua-1][itel-1][izz]>=0)
					{
					  ii=izz;
					  break;
					}
				    }
				  //				  zval=iz;
				  zval=ii;
				  fazietto.alimSi1Si2[ibl][iqua-1][itel-1]=zval;
				  break;
				}
			    }//for(;;)

			    }

			    }
			  
			}
		    }
		    

	      fclose(ap);
          //getchar();
		}
	      //cout<<"silici-silici"<<endl;
	      for(int ibl=0;ibl<12;ibl++)
		{
		   if(fazietto.dist[ibl]>0)
		     {
		       for(int iqua=0;iqua<4;iqua++)
			 {
			   for(int itel=0;itel<4;itel++)
			     {
			       //cout<<"ibl="<<ibl<<" "<<iqua+1<<" "<<itel+1<<" "<<fazietto.alimSi1Si2[ibl][iqua][itel]<<endl;
			       for(int iz=1;iz<=fazietto.alimSi1Si2[ibl][iqua][itel];iz++)
				 {
				   //  cout<<iz<<" "<<fazietto.infmassesi1si2[ibl][iqua][itel][iz]<<" "<<fazietto.supmassesi1si2[ibl][iqua][itel][iz]<<endl;
				 }
			     }
			 }
		     }
		}
	    }


if(strcmp(riv,"FAZIA_GRID_Si2CsI")==0)//griglie Si2CsI
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=215;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("Griglie Si2CsI per fazia=%s\n",tipo);
	      cout<<"Si leggono i limiti di Z per identificazione in massa per Si2 CsI"<<endl;
	      int ibl,iqua,itel;

	      char riga2[1000];
	      //	      int cl1,cl2,cl3,cl4;
	      TString slocal;
	      TString sloc;
	      int num,zval;
	      if(ap!=0)
		{
		  int Ik2=-2;
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  slocal=riga2;
			  if(slocal.Contains("<PARAMETER> IDTelescopes"))
			    {
			      Ik2=-2;
			      // cout<<slocal.Data()<<endl;
			  int ik=slocal.Index("CSI_");
			  sscanf(&riga2[ik+4],"%d",&num);
			  
			  if(num>100)
			    {
			      ibl=num/100;
			    }
			  else
			    {
			      ibl=0;
			    }
			  iqua=(num-ibl*100)/10;
			  itel=num-ibl*100-iqua*10;
			  //cout<<"N="<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<endl;
			  
			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      zval=-1;
			      if(sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				  //cout<<"qui "<<Ik2<<endl;
				  if(sloc.Contains("-"))
				    {
				      //Ik2=sloc.Index("-")+1;
				      Ik2=-1;
				      break;
				    }
				  else
				    {
				  Ik2=sloc.Index("PIDRANGE=")+9;
				    
				  sscanf(&riga2[Ik2],"%d",&zval);
				  break;
				    }
				}
			      //cout<<"E="<<Ik2<<endl;
			      if(Ik2==-2 && sloc.Contains("ID:Z="))
				{
				  break;
				}
			    }//for(;;)
			  //cout<<"Ik2="<<Ik2<<endl;
			  if(Ik2==-2)continue;
			  if(zval>0)
			    {
			  fazietto.alimSi2CsI[ibl][iqua-1][itel-1]=zval;
			  //cout<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<" "<<zval<<endl;
			    }
			  if(Ik2>=-1)
			    {
			  int iz,alom,aloM;
			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      //cout<<riga2<<endl;
			      
			      if(sloc.Contains("<PARAMETER> PIDRANGE")&&!sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				 
				  int ik2=sloc.Index("PIDRANGE")+8;
				  sscanf(&riga2[ik2],"%d",&iz);
				  int ik3=sloc.Index("=")+1;
				  sscanf(&riga2[ik3],"%d",&alom);
				  if(sloc.Contains("|"))
				    {

				  int ik4=sloc.Last('|')+1;
				  sscanf(&riga2[ik4],"%d",&aloM);

				  fazietto.infmassesi2csi[ibl][iqua-1][itel-1][iz]=alom;
				  fazietto.supmassesi2csi[ibl][iqua-1][itel-1][iz]=aloM;
				    }
				  //cout<<ibl<<" "<<iqua<<" "<<itel<<" "<<zval<<" "<<iz<<" "<<alom<<" "<<aloM<<endl;
				  if(iz==zval)
				    {break;}
				}
			      if((sloc.Contains("+KVIDZALine")||sloc.Contains("+KVIDCutContour"))&&zval==-1)
				{
				  int ii=-1;
				  for(int izz=iz;izz>=0;izz--)
				    {
				      if(fazietto.infmassesi2csi[ibl][iqua-1][itel-1][izz]>=0)
					{
					  ii=izz;
					  break;
					}
				    }
				  //				  zval=iz;
				  zval=ii;

				 
				  fazietto.alimSi2CsI[ibl][iqua-1][itel-1]=zval;
				  break;
				}

			    }//for(;;)
			  //cout<<"FF="<<fazietto.alimSi2CsI[ibl][iqua-1][itel-1]<<endl;
			    }
			    }//if(slocal.Contains("<PARAMETER> IDTelescopes"))
			 
			  
			}
		    }
		    

	      fclose(ap);
          //getchar();
		}
	      //cout<<"silici-cesi"<<endl;
	      for(int ibl=0;ibl<12;ibl++)
		{
		   if(fazietto.dist[ibl]>0)
		     {
		       for(int iqua=0;iqua<4;iqua++)
			 {
			   for(int itel=0;itel<4;itel++)
			     {
			       //cout<<"ibl="<<ibl<<" "<<iqua+1<<" "<<itel+1<<" "<<fazietto.alimSi2CsI[ibl][iqua][itel]<<endl;
			       for(int iz=1;iz<=fazietto.alimSi2CsI[ibl][iqua][itel];iz++)
				 {
				   //cout<<iz<<" "<<fazietto.infmassesi2csi[ibl][iqua][itel][iz]<<" "<<fazietto.supmassesi2csi[ibl][iqua][itel][iz]<<endl;
				 }
			     }
			 }
		     }
		}

	    }


if(strcmp(riv,"FAZIA_GRID_Si1PSA")==0)//griglie PSA corrente Si1
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=215;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("Griglie Si1 PSA corrente per fazia=%s\n",tipo);
	      cout<<"Si leggono i limiti di Z per identificazione in massa per PSA in corrente in Si1"<<endl;
	      int ibl,iqua,itel;

	      char riga2[1000];
	      //int cl1,cl2,cl3,cl4;
	      TString slocal;
	      TString sloc;
	      int num,z1,z2;
	      if(ap!=0)
		{

		  int Ik2=-2;
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  slocal=riga2;
			  if(slocal.Contains("<PARAMETER> IDTelescopes"))
			    {
			      Ik2=-2;
			      // cout<<slocal.Data()<<endl;
			  int ik=slocal.Index("SI1_");
			  sscanf(&riga2[ik+4],"%d",&num);
			  
			  if(num>100)
			    {
			      ibl=num/100;
			    }
			  else
			    {
			      ibl=0;
			    }
			  iqua=(num-ibl*100)/10;
			  itel=num-ibl*100-iqua*10;
			  //cout<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<endl;
			    z2=-1;
			    z1=-1;
			    //cout<<"Z2="<<z2<<endl;
			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      z1=-1;
			      z2=-1;
			      //cout<<riga2<<endl;
			      if(sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				  int ik2=sloc.Index("PIDRANGE=")+9;
				  sscanf(&riga2[ik2],"%d",&z1);
				  int ik3=sloc.Index("-")+1;
				  sscanf(&riga2[ik3],"%d",&z2);
				  break;
				}
			      if(z2==-1&&(sloc.Contains("OK:")||sloc.Contains("ID:Z"))) break;
			    }
			  //			  			    cout<<"bZ2="<<z2<<endl;
			  //cout<<z2<<" "<<z1<<endl;
			  if(z2==-1)continue;
			  if(z2<200&&z2>=0)
			    {
			  fazietto.aliminfSi1PSA[ibl][iqua-1][itel-1]=z1;
			  fazietto.alimsupSi1PSA[ibl][iqua-1][itel-1]=z2;
			    }
			  //cout<<num<<" "<<ibl<<" "<<iqua<<" "<<itel<<" "<<z1<<" "<<z2<<endl;
			  int iz,alom,aloM;
			
			  for(;;)
			    {
			      int oo= fscanf(ap,"\n%[^\n]",riga2);
			      if(oo==EOF)
				{
				  break;
				}
			      sloc=riga2;
			      //cout<<riga2<<endl;
			      if(sloc.Contains("<PARAMETER> PIDRANGE")&&!sloc.Contains("<PARAMETER> PIDRANGE="))
				{
				  //cout<<"dentro"<<endl;
				  int ik2=sloc.Index("PIDRANGE")+8;
				  sscanf(&riga2[ik2],"%d",&iz);
				  int ik3=sloc.Index("=")+1;
				  sscanf(&riga2[ik3],"%d",&alom);
				  int ik4=sloc.Last('|')+1;
				  sscanf(&riga2[ik4],"%d",&aloM);
				  fazietto.infmassepsa[ibl][iqua-1][itel-1][iz]=alom;
				  fazietto.supmassepsa[ibl][iqua-1][itel-1][iz]=aloM;
				  //cout<<alom<<" "<<aloM<<" "<<iz<<endl;

				  //	  				  cout<<ibl<<" "<<iqua<<" "<<itel<<" "<<z2<<" "<<iz<<" "<<alom<<" "<<aloM<<endl;
				  if(iz==z2)
				    {break;}
				}


			      if((sloc.Contains("+KVIDZALine")||sloc.Contains("+KVIDCutContour")||sloc.Contains("+KVIDZoneContour"))&&z2==200)
				{

				  //cout<<"EEE"<<endl;
				  int ii=-1;
				  for(int izz=iz;izz>=0;izz--)
				    {
				      if(fazietto.infmassepsa[ibl][iqua-1][itel-1][izz]>=0)
					{
					  ii=izz;
					  break;
					}
				    }
				  //				  zval=iz;
				  z2=ii;
				  //cout<<z2<<endl;
				  ii=-1;
				  for(int izz=1;izz<=iz;izz++)
				    {
				      if(fazietto.infmassepsa[ibl][iqua-1][itel-1][izz]>=0)
					{
					  ii=izz;
					  break;
					}
				    }
				  //				  zval=iz;
				  z1=ii;


			 
				 
				  fazietto.alimsupSi1PSA[ibl][iqua-1][itel-1]=z2;
				  fazietto.aliminfSi1PSA[ibl][iqua-1][itel-1]=z1;
				  break;
				}

			    }//for(;;)
			    }//if(slocal.Contains("<PARAMETER> IDTelescopes"))


			  
			}
		    }
		    

	      fclose(ap);
          //getchar();
		}
	      //cout<<"silicio psa"<<endl;
	      for(int ibl=0;ibl<12;ibl++)
		{
		   if(fazietto.dist[ibl]>0)
		     {
		       for(int iqua=0;iqua<4;iqua++)
			 {
			   for(int itel=0;itel<4;itel++)
			     {
			       //cout<<"ibl="<<ibl<<" "<<iqua+1<<" "<<itel+1<<" "<<fazietto.aliminfSi1PSA[ibl][iqua][itel]<<" "<<fazietto.alimsupSi1PSA[ibl][iqua][itel]<<endl;
			       for(int iz=fazietto.aliminfSi1PSA[ibl][iqua][itel];iz<=fazietto.alimsupSi1PSA[ibl][iqua][itel];iz++)
				 {
				   //cout<<iz<<" "<<fazietto.infmassepsa[ibl][iqua][itel][iz]<<" "<<fazietto.supmassepsa[ibl][iqua][itel][iz]<<endl;
				 }
			     }
			 }
		     }
		}
	    }


if(strcmp(riv,"FAZIA_LIMITI_PSA")==0)//limiti in energia lasciata in Si1 per l'identificazione in Z e A per la PSA in Si1 (Imax-E)
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=3215;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("Limiti di energia in Si1 per identificazione in Z e A da PSA Imax-E in Si1=%s\n",tipo);
	      cout<<"Si leggono i limiti di energia in Si1 per identificazione in Z e A da PSA Imax-E in Si1"<<endl;
	      int ibl,iqua,itel;

	      char riga2[1000];

	      float ezinf,ezsup,eainf,easup;
	      int zval;
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d %d %f %f %f %f",&ibl,&iqua,&itel,&zval,&ezinf,&ezsup,&eainf,&easup);
			  fazietto.zesi1inf[ibl][iqua-1][itel-1][zval]=ezinf;
			  fazietto.zesi1sup[ibl][iqua-1][itel-1][zval]=ezsup;
			  fazietto.aesi1inf[ibl][iqua-1][itel-1][zval]=eainf;
			  fazietto.aesi1sup[ibl][iqua-1][itel-1][zval]=easup;
			  
			}


			  
		    }
		  fclose(ap);
		}
		    

	    
          //getchar();
	    
	    }


if(strcmp(riv,"FAZIA_SPESSORI")==0)//spessori dei rivelatori di Fazietto
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=1;
	      iflag=155;
	      
		}


	  if(ispefazietto==1 && iflag==0)
	    {
	      int ib,iqua,itel,spe1,spe2;
	      sscanf(riga,"%d %d %d %d %d",&ib,&iqua,&itel,&spe1,&spe2);
	      fazietto.spes_si1[ib][iqua-1][itel-1]=spe1;
	      fazietto.spes_si2[ib][iqua-1][itel-1]=spe2;


	    }



		//!Aggiunto per inserire nel MC i tagli in LO sui Cesi
		if(strcmp(riv,"SOGLIE_GARF")==0)
		{

			sscanf(riga,"%s %s",riv,tipo);
			FILE *ap=fopen(tipo,"r");
			printf("soglie=%s\n",tipo);
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=51;			
			int icsi,isec;
			
			char riga2[1000];
			float soglia[5];
			int n;
			if(ap!=0)
			{
				while(fgets(riga2,1000,ap)!=NULL)
				{
					if(riga2[0]!='#')
					{
						n=sscanf(riga2,"%d %d %f %f %f %f %f",&isec,&icsi,soglia,soglia+1,soglia+2,soglia+3,soglia+4);
						if(n>6) {
							garf.cut_csi_p[icsi-1][isec-1]=soglia[0];
							garf.cut_csi_d[icsi-1][isec-1]=soglia[1];
							garf.cut_csi_t[icsi-1][isec-1]=soglia[2];
							garf.cut_csi_3he[icsi-1][isec-1]=soglia[3];
							garf.cut_csi_a[icsi-1][isec-1]=soglia[4];
						}
						else {
							if(n>3) {
								garf.cut_csi_p[icsi-1][isec-1]=soglia[0];
								garf.cut_csi_a[icsi-1][isec-1]=soglia[1];
							}
						}
					}
				}
			}
			fclose(ap);
		}
		if(strcmp(riv,"SOGLIE_RCO_CSI")==0)
		{
			sscanf(riga,"%s %s",riv,tipo);
			FILE *ap=fopen(tipo,"r");
			printf("soglie=%s\n",tipo);
			
			int icsi,isec;
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=52;			
			char riga2[1000];
			float soglia_p,soglia_a;
			if(ap!=0)
			{
				while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
				{
					if(riga2[0]!='#')
					{
						sscanf(riga2,"%d %d %f %f",&isec,&icsi,&soglia_p,&soglia_a);
						rco.soglie_csi_p[isec-1][icsi-1]=soglia_p;
						rco.soglie_csi_a[isec-1][icsi-1]=soglia_a;
					}
				}
			}
			fclose(ap);
		}
	  if(strcmp(riv,"CODICE_RCO_GAS")==0)
	    {
#ifndef _ODIE_
	     TNamed *name4=new TNamed(riga,"s=1300");
	     Classe_analisi::Getanalisi()->Getfout()->cd();
	     name4->Write();
#endif
			iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=15;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int isec;

	      char riga2[1000];
	      int codex;
	      
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d ",&isec,&codex);
			  rco.codice_gas[isec-1]=codex;
			 
			     
			    
			}
		    }
		    }
		  fclose(ap);
		}
	  if(strcmp(riv,"CODICE_RCO_STRIP")==0)
	    {
#ifndef _ODIE_
	      TNamed *name4=new TNamed(riga,"s=1301");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name4->Write();
#endif
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=16;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int isec,istrip;

	      char riga2[1000];
	      int codex;
	      
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d",&isec,&istrip,&codex);
			  rco.codice_strip[isec-1][istrip-1]=codex;
			 
			     
			    
			}
		    }
		    }
		  fclose(ap);
		}

	  if(strcmp(riv,"INTERVALLI_MASSE")==0)
	    {
#ifndef _ODIE_
	      TNamed *name44=new TNamed(riga,"s=13301");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name44->Write();
#endif
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=166;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int isec,istrip,icsi;

	      char riga2[1000];
	      //	      int codex;
	      
	      if(ap!=0)
		{
		  rco.intervalli=1;
		  int zmi,zma;
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  //			  cout<<riga2<<endl;
			  sscanf(riga2,"%d %d %d %d %d",&isec,&istrip,&icsi,&zmi,&zma);
			  rco.zinf[isec-1][istrip-1][icsi-1]=zmi;
			  rco.zsup[isec-1][istrip-1][icsi-1]=zma;
			  // cout<<isec<<" "<<istrip<<" "<<icsi<<" "<<zmi<<" "<<zma<<endl;
			  int ij=0;
			  int kj=0;
			  int kv=0;
			  int amin,amax;
			  for(int k=0;k<strlen(riga2);k++)
			    {
			      if(riga2[k]==' ' && riga2[k+1]!=' ')
				{
				  ij++;
				}
			      if(ij==4)
				{
				  kj=k+1;
				  break;
				}
			    }
			  for(int k=kj;k<strlen(riga2);k++)
			    {
			      if(riga2[k]==' ')
				{
				  kv=k;
				  break;
				}
			    }
			  for(int i=zmi;i<=zma;i++)
			    {
			      ij=0;
			      sscanf(&riga2[kv],"%d %d",&amin,&amax);
			      //  cout<<i<<" "<<amin<<" "<<amax<<endl;
			      rco.ainf[isec-1][istrip-1][icsi-1][i-1]=amin;
			      rco.asup[isec-1][istrip-1][icsi-1][i-1]=amax;
			      
			  
			      for(int k=kv;k<strlen(riga2);k++)
				{
				  if(riga2[k]==' ' && riga2[k+1]!=' ')
				{
				  ij++;
				}
				  if(ij==3)
				    {
				      kv=k;
				      break;
				    }
				}
			    }
			  //			  rco.codice_strip[isec-1][istrip-1]=codex;
			 
			     
			    
			}
		    }
		    }
		  fclose(ap);
		}
	  
	  if(strcmp(riv,"CODICE_RCO_CSI")==0)
	    {
#ifndef _ODIE_
	     TNamed *name4=new TNamed(riga,"s=1302");
	     Classe_analisi::Getanalisi()->Getfout()->cd();
	     name4->Write();
#endif
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=17;
	      sscanf(riga,"%s %s",riv,tipo);
	      FILE *ap=fopen(tipo,"r");
	      printf("codice=%s\n",tipo);

	      int isec,icsi;

	      char riga2[1000];
	      int codex;
	      
	      if(ap!=0)
		{
		  while(fscanf(ap,"\n%[^\n]",riga2)!=EOF)
		    {
		      if(riga2[0]!='#')
			{
			  sscanf(riga2,"%d %d %d",&isec,&icsi,&codex);
			  rco.codice_csi[isec-1][icsi-1]=codex;
			 
			     
			    
			}
		    }
		    }
		  fclose(ap);
		}

		
	  if(strcmp(riv,"PRESSIONE_GARFIELD")==0)
	    {iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=3;
	      sscanf(riga,"%s %f",riv,&garf.pressione);
#ifndef _ODIE_
	      TNamed *name2=new TNamed(riga,"s=1001");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name2->Write();
#endif
	    }
	  //<<<<<<< Classe_geo.cxx
	  if(strcmp(riv,"PRESSIONE_GARFIELD_BACK")==0)
	    {iphos=0;
	      igarf_f=0;
	      ihector=0;
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=6;
	      sscanf(riga,"%s %f",riv,&garf.pressione_back);
#ifndef _ODIE_
		  TNamed *name4=new TNamed(riga,"s=1004");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name4->Write();
#endif
	    }
	  if(strcmp(riv,"CUTG")==0)
	    {

	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=27;
	      sscanf(riga,"%s %s %s",riv,tipo,denominaz);
	      TFile *filecut=0;
	      filecut=new TFile(Form("%s",tipo));
	      if(!filecut->IsZombie())
		{
		  Classe_analisi::Getanalisi()->gggcuts[ncuts]=(TGraph*)filecut->Get(Form("%s",denominaz));
		  Classe_analisi::Getanalisi()->gggcuts[ncuts]->SetName(Form("%s",denominaz));
		  Classe_analisi::Getanalisi()->gggcuts[ncuts]->SetTitle(Form("%s",denominaz));

		  ncuts++;
		  cout<<"carico il cut "<<denominaz<<" da "<<tipo<<endl;
		}

	    
		}

	  if(strcmp(riv,"TABELLA_PERDITE_GARFIELD")==0)
	    {
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      //<<<<<<< Classe_geo.cxx
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=4;
	      char file_perdite[200];
	      sscanf(riga,"%s %s",riv,file_perdite);
#ifndef _ODIE_
		  if(Classe_analisi::Getanalisi()->tipo_analisi>=100 && Classe_analisi::Getanalisi()->tipo_analisi<200)
		  {
#endif
	      TabellaPerdite(file_perdite);
#ifndef _ODIE_
	      TNamed *name3=new TNamed(riga,"s=1002");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name3->Write();
		  }
		  else
		  {
			  printf("Non carico la tabella %s perche' analisi exp o mcarlo 4p\n",file_perdite);
		  }
#endif
	    }
	  //<<<<<<< Classe_geo.cxx
	  if(strcmp(riv,"TABELLA_PERDITE_RCO")==0)
	    {
	      iphos=0;
	      igarf_f=0;
	      ihector=0;
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=7;
	      char file_perdite[200];
	      sscanf(riga,"%s %s",riv,file_perdite);
#ifndef _ODIE_
		  if(Classe_analisi::Getanalisi()->tipo_analisi>=100&&Classe_analisi::Getanalisi()->tipo_analisi<200)
		  {
#endif
			  TabellaPerditeRCO(file_perdite);
#ifndef _ODIE_
			  TNamed *name73=new TNamed(riga,"s=1073");
			  Classe_analisi::Getanalisi()->Getfout()->cd();
			  name73->Write();
		  }
		  else
		  {
			  printf("Non carico la tabella %s perche' analisi exp o mcarlo 4p\n",file_perdite);
		  }
#endif
	    }
	  if(strcmp(riv,"PRESSIONE_RCO")==0)
	    {iphos=0;
	      igarf_f=0;
	      ihector=0;
	      igarf_b=0;
	      irco_strip=0;
	      irco_gas=0;
	      irco_csi=0;
	      ifazietto=0;
	      ispefazietto=0;
	      iflag=3;
	      sscanf(riga,"%s %f",riv,&rco.pressione);
#ifndef _ODIE_
		  TNamed *name20=new TNamed(riga,"s=1020");
	      Classe_analisi::Getanalisi()->Getfout()->cd();
	      name20->Write();
#endif
	    }


	}
    }
    
    if(flock(fileno(apri),LOCK_UN)) perror("flock");
  fclose(apri);


  
}

int Classe_geo::InsidePhos(float thetapart,float phipart,float vpart,float *tvolo,int *code)
{

  int inside=0;

	  *code=-1;
  for(int ip=0;ip<6;ip++)
    {
      for(int ih=0;ih<9;ih++)
	{

	  if(phos.dist[ip][ih]>0) 
	    {
	      if(thetapart-phos.theta[ip][ih]<45)
	    {

	      float t=(pow(phos.veccentro[ip][ih][0],2)+pow(phos.veccentro[ip][ih][1],2)+pow(phos.veccentro[ip][ih][2],2))/(phos.veccentro[ip][ih][0]*sin(thetapart/57.296)*cos(phipart/57.296)+phos.veccentro[ip][ih][1]*sin(thetapart/57.296)*sin(phipart/57.296)+phos.veccentro[ip][ih][2]*cos(thetapart/57.296));
	      float punto[3];
	      punto[0]=t*sin(thetapart/57.296)*cos(phipart/57.296);
	      punto[1]=t*sin(thetapart/57.296)*sin(phipart/57.296);
	      punto[2]=t*cos(thetapart/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-phos.veccentro[ip][ih][nn];
		}
	      float pi,pj;
	      Classe_formule::scaprod(punto_c,phos.vi[ip][ih],&pi);
	      Classe_formule::scaprod(punto_c,phos.vj[ip][ih],&pj);
	      // printf("%f %f %f %f %f\n",pi,pj,phos.vi[ip][ih][0],phos.vi[ip][ih][1],phos.vi[ip][ih][2]);
	      if(TMath::Abs(pi)-phos.dx/2<=0 && TMath::Abs(pj)-phos.dy/2<=0)
		{
		  
		  inside=1;
		  float basevolo=TMath::Sqrt(pow(punto[0],2)+pow(punto[1],2)+pow(punto[2],2));
		  *tvolo=basevolo/vpart;
		  //		  if(phos.codice[ip][ih]>5|| phos.codice[ip][ih]<0) //phos cattivi
		  if(phos.codice[ip][ih]<=0) //phos cattivi
		    {
		      inside=2;
		    }
		  *code=(ip+1)*10+(ih+1);
		  return inside;
		}

	    }	


	    }
	}
    }


  return inside;
}

int Classe_geo::InsideFazietto(float thetapart,float phipart,float vpart,float *distanza,int *code)
{
  *code=-1;
  *distanza=-1;
  int inside=0;
 
  for(int j=0;j<12;j++)
    {
      //      cout<<"blocco===="<<j<<" "<<fazietto.dist[j]<<endl;
      if(fazietto.dist[j]>0)
	{
	if(thetapart-fazietto.theta[j]<90)//c'e' il blocco
	{
	      float t=(pow(fazietto.vcentroblocco[j][0],2)+pow(fazietto.vcentroblocco[j][1],2)+pow(fazietto.vcentroblocco[j][2],2))/(fazietto.vcentroblocco[j][0]*sin(thetapart/57.296)*cos(phipart/57.296)+fazietto.vcentroblocco[j][1]*sin(thetapart/57.296)*sin(phipart/57.296)+fazietto.vcentroblocco[j][2]*cos(thetapart/57.296));
	      float punto[3];
	      punto[0]=t*sin(thetapart/57.296)*cos(phipart/57.296);
	      punto[1]=t*sin(thetapart/57.296)*sin(phipart/57.296);
	      punto[2]=t*cos(thetapart/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-fazietto.vcentroblocco[j][nn];
		}
	      float pi,pj;
	      Classe_formule::scaprod(punto_c,fazietto.viblocco[j],&pi);
	      Classe_formule::scaprod(punto_c,fazietto.vjblocco[j],&pj);
	      // cout<<"geo blocco="<<j<<" "<<pi<<" "<<pj<<" "<<thetapart<<" "<<phipart<<endl;
	      	      for(int k=0;k<16;k++)
	      

		{
		  
  if(TMath::Abs(fazietto.vx[j][k]-pi)<fazietto.dx_active/2 && TMath::Abs(fazietto.vy[j][k]-pj)<fazietto.dy_active/2)
		 //  float xmin=fazietto.vx[j][k]-fazietto.dx_active/2;
// 		  float xmax=fazietto.vx[j][k]+fazietto.dx_active/2;
// 		  float ymin=fazietto.vy[j][k]-fazietto.dy_active/2;
// 		  float ymax=fazietto.vy[j][k]+fazietto.dy_active/2;
// 		  if(pi>=xmin && pi<=xmax && pj>=ymin && pj<=ymax)

		{
		  int iqua=k/4;
		  int itel=k-iqua*4;
		  // cout<<"inside "<<k<<" bl="<<j<<"q="<<iqua+1<<" t="<<itel+1<<" "<<fazietto.rotti[j][iqua][itel]<<" evento="<<Classe_analisi::Getanalisi()->nentry<<endl;
		  // cout<<j<<" "<<k<<" "<<iqua<<" "<<itel<<" "<<fazietto.rotti[j][iqua][itel]<<endl;
		  if(fazietto.rotti[j][iqua][itel]==0)
		    {
		      inside=0;
		    
		      return inside;
		    }

		  inside=1;
		 *code=j*100+k;
		 
		 // cout<<"inside"<<j<<" "<<k<<" "<<iqua<<" "<<itel<<endl;
		 *distanza=sqrt(pi*pi+pj*pj+fazietto.dist[j]*fazietto.dist[j]);
		 return inside;
		}

		}
	}
	}
    }


  return inside;

}

//da provare se ci sara' un montecarlo con i gamma
int Classe_geo::InsideHector(float thetapart,float phipart, float *tvolo, int *code)
{
  int inside=0;
  *code=-1;
  for(int ip=0;ip<8;ip++)
    {
      if(hector.dist[ip]>0)
	{
	  if(thetapart-hector.theta[ip]<45)
	    {
	      float t=(pow(hector.veccentro[ip][0],2)+pow(hector.veccentro[ip][1],2)+pow(hector.veccentro[ip][2],2))/(hector.veccentro[ip][0]*sin(thetapart/57.296)*cos(phipart/57.296)+hector.veccentro[ip][1]*sin(thetapart/57.296)*sin(phipart/57.296)+hector.veccentro[ip][2]*cos(thetapart/57.296));
	      float punto[3];
	      punto[0]=t*sin(thetapart/57.296)*cos(phipart/57.296);
	      punto[1]=t*sin(thetapart/57.296)*sin(phipart/57.296);
	      punto[2]=t*cos(thetapart/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-hector.veccentro[ip][nn];
		}

	      float pi,pj;
	      Classe_formule::scaprod(punto_c,hector.vi[ip],&pi);
	      Classe_formule::scaprod(punto_c,hector.vj[ip],&pj);
	      if(TMath::Sqrt(pow(pi,2)+pow(pj,2))<=hector.raggio)
		{

		  inside=1;
		  float basevolo=TMath::Sqrt(pow(punto[0],2)+pow(punto[1],2)+pow(punto[2],2));
		  *tvolo=basevolo/Classe_formule::cluce;

		  *code=ip+1;
		  return inside;

		}



	    }
	}


    }
  return inside;

}

int Classe_geo::InsideGarf(float thetapart,float phipart,float zpart,int *code,int *code_micro)
{
  //microstrip
  //code_micro=1 left, 2 right
  //si mette 0 per Z=1,2 per i quali la micro non si usa (neppure l'info se e' dx o sx)
  int inside=0;
  *code_micro=-1;
  *code=-1;
  int ip0=-1;
  int isec0=-1;
  if(thetapart<garf.themincsi[7] && garf.thecsi[7]>=0)
    {
      return inside;
    }
  if(thetapart>garf.themaxcsi[0] && garf.thecsi[0]>=0)
    {
      return inside;
    }
  if(thetapart>garf.themaxcsi[4] && garf.thecsi[4]>=0 && garf.thecsi[0]<0)
    {
      return inside;
    }
  if(thetapart<garf.themincsi[3] && garf.thecsi[3]>=0 && garf.thecsi[7]<0)
    {
      return inside;
    }
    
//si cercano le coordinate della particella sul piano del rivelatore
  int ipmin=0;
  int ipmax=8;
  if(thetapart<90){ipmin=4;}
  if(thetapart>=90){ipmax=4;}
  //  for(int ip=0;ip<8;ip++)
  for(int ip=ipmin;ip<ipmax;ip++)
    {
            if(garf.thecsi[ip]>=0 && (thetapart-garf.thecsi[ip]<45))
      // if(garf.thecsi[ip]>=0 && thetapart>=garf.themincsi[ip]&&thetapart<garf.themaxcsi[ip])
	{
	  int isecmin=0;
	  int isecmax=24;
	  if(ip<4){isecmin=2;isecmax=23;}
	  for(int isec=isecmin;isec<isecmax;isec++)
	    {

	      float punto[3];

	      	      float t=(pow(garf.veccentro[isec][ip][0],2)+pow(garf.veccentro[isec][ip][1],2)+pow(garf.veccentro[isec][ip][2],2))/(garf.veccentro[isec][ip][0]*sin(thetapart/57.296)*cos(phipart/57.296)+garf.veccentro[isec][ip][1]*sin(thetapart/57.296)*sin(phipart/57.296)+garf.veccentro[isec][ip][2]*cos(thetapart/57.296));



	      	      punto[0]=t*sin(thetapart/57.296)*cos(phipart/57.296);
	       punto[1]=t*sin(thetapart/57.296)*sin(phipart/57.296);
	      punto[2]=t*cos(thetapart/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-garf.veccentro[isec][ip][nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y ortogonale alle due basi del trapezio
	      Classe_formule::scaprod(punto_c,garf.vi[isec][ip],&pi);
	      Classe_formule::scaprod(punto_c,garf.vj[isec][ip],&pj);	   
	      //double *xpcsi;
	      //double *ypcsi;
	      //xpcsi=garf.gcsi[ip]->GetX();
	      //ypcsi=garf.gcsi[ip]->GetY();


	      if(garf.gcsi[ip]->IsInside(pi,pj)==1)
		{
		  ip0=ip;
		  isec0=isec;

		  break;
		}
	    }
	}
    }

    
  if(ip0>=0 && isec0>=0)
  {
	//cesi rotti
	 if(garf.codice[ip0][isec0]==1)
	   {
	     inside=0;
	     return inside;
	   }

      if(garf.micro_rotte[isec0][ip0][0]==0&&garf.micro_rotte[isec0][ip0][1]==0&&garf.micro_rotte[isec0][ip0][2]==0&&garf.micro_rotte[isec0][ip0][3]==0 && zpart>2.5)//tutte le micro rotte e siamo sopra le alfa
	{
	  inside=0;
	  return inside;
	}

	

	if(ip0<4 &&zpart>2.5) //la micro non e' utilizzabile nella camera back e dal solo cesio non ci si fa niente
	  {
	    inside=0;
	    return inside;
	  }
	  
      inside=1;
      *code=(isec0+1)*10+ip0+1;


       if(zpart<=2.5)//per Z=1,2 l'info delle micro non si usa
       	{
       	  *code_micro=0; 
       	}
       else
       	{

      if(phipart<garf.phi[ip0][isec0])
	{
	  *code_micro=1;//sx
	}
      else
	{
	  *code_micro=2;//dx
	 
	  
	}
	 } //zpart <=2.5

      if(*code_micro==1 && garf.micro_rotte[isec0][ip0][0]==0)//up left rotta
	{
	  inside=0;
	  return inside;
	}
      if(*code_micro==2 && garf.micro_rotte[isec0][ip0][1]==0) //up right rotta
	{
	  inside=0;
	  return inside;
	}

      return inside;
     }
  


//   for(int ip=0;ip<8;ip++)
//     {
//       if(garf.thecsi[ip]>=0)
// 	{
// 	  if(thetapart>=garf.themincsi[ip]&&thetapart<=garf.themaxcsi[ip])
// 	    {
// 	      int isecmin=0;
// 	      int isecmax=24;
// 	      if(ip<4){isecmin=2;isecmax=23;}
//       for(int isec=isecmin;isec<isecmax;isec++)
// 	{
// 	  if(phipart>=garf.phi[ip][isec]-garf.dphi_corr/2 && phipart<garf.phi[ip][isec]+garf.dphi_corr/2)
// 	    {
//       *code=(isec+1)*10+ip+1;    
//       if(phipart<garf.phi[ip][isec])
// 	{
// 	  *code_micro=1;
// 	  if(ip<4){*code_micro=291;if(zpart>2.5){inside=0;return inside;}}//la micro non e' utilizzabile mai nella camera back
// 	}
//       else
// 	{
// 	  *code_micro=2;
// 	  if(ip<4){*code_micro=292;if(zpart>2.5){inside=0;return inside;}}//la micro non e' utilizzabile mai nella camera back
// 	} 
//       if(ip>=4)
// 	{
// 	  if(isec==21&&phipart<garf.phi[ip][isec] &&zpart>=2.5){*code_micro=301;inside=0; return inside;}
//       if(isec==21&&phipart>=garf.phi[ip][isec]&&zpart>=2.5){*code_micro=302;}
// 	}
//       inside=1;

//       return inside;
// 	    }

// 	}
// 	} 
//     }
//     }
  return inside;

}
//<<<<<<< Classe_geo.cxx

void Classe_geo::InizializzaRCO()
{
 if(rco.gas_a0<0) // si prendono i parametri del settore 0, tanto sono tutti uguali
    {
	      rco.gas_a8=rco.gas_dist[0]*TMath::Tan(rco.strip_themin[0][7]/57.296);
	      rco.gas_a0=rco.gas_dist[0]*TMath::Tan(rco.strip_themax[0][0]/57.296);
	      rco.gas_acc=rco.gas_a0/TMath::Cos(rco.dphi/(2*57.296));
	      rco.gas_gcc=rco.gas_a8/TMath::Cos(rco.dphi/(2*57.296));
	      rco.gas_xa=rco.gas_acc*TMath::Cos(3.14/2-rco.dphi/(2*57.296));
	      rco.gas_ya=rco.gas_acc*TMath::Sin(3.14/2-rco.dphi/(2*57.296));
	      rco.gas_yg=rco.gas_gcc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.strip_a8=rco.strip_dist[0][7]*TMath::Tan(rco.strip_themin[0][7]/57.296);
	      rco.strip_gcc=rco.strip_a8/TMath::Cos(rco.dphi/(2*57.296));
	      rco.strip_yg=rco.strip_gcc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_a0=rco.csi_dist[0][0]*TMath::Tan(rco.strip_themax[0][0]/57.296);
	      rco.csi_acc=rco.csi_a0/TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_xa=rco.csi_acc*TMath::Sin(rco.dphi/(2*57.296));
	      rco.csi_ya=rco.csi_acc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_yb=rco.csi_ya;
	      rco.csi_xb=-rco.csi_xa;
	      rco.csi_xA=rco.spe_morto_csi;
	      rco.csi_xB=-(float)(rco.spe_morto_csi);
	      rco.csi_yA=rco.csi_ya;
	      rco.csi_yB=rco.csi_yb;
	      rco.csi_a2=rco.csi_dist[0][0]*TMath::Tan(rco.strip_themax[0][2]/57.296);
	      rco.csi_ccc=rco.csi_a2/TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_xc=rco.csi_ccc*TMath::Sin(rco.dphi/(2*57.296));
	      rco.csi_yc=rco.csi_ccc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_yd=rco.csi_yc;
	      rco.csi_xd=-rco.csi_xc;
	      rco.csi_yD=rco.csi_yd;
	      rco.csi_yC=rco.csi_yc;	      
	      rco.csi_xC=rco.spe_morto_csi;
	      rco.csi_xD=-(float)(rco.spe_morto_csi);
	      rco.csi_a4=rco.csi_dist[0][0]*TMath::Tan(rco.strip_themax[0][4]/57.296);
	      rco.csi_ecc=rco.csi_a4/TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_xe=rco.csi_ecc*TMath::Sin(rco.dphi/(2*57.296));
	      rco.csi_ye=rco.csi_ecc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_yf=rco.csi_ye;
	      rco.csi_xf=-rco.csi_xe;
	      rco.csi_yF=rco.csi_yf;
	      rco.csi_yE=rco.csi_ye;	      
	      rco.csi_xE=rco.spe_morto_csi;
	      rco.csi_xF=-(float)(rco.spe_morto_csi);	      
	      rco.csi_a8=rco.csi_dist[0][0]*TMath::Tan(rco.strip_themin[0][7]/57.296);
	      rco.csi_gcc=rco.csi_a8/TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_xg=rco.csi_gcc*TMath::Sin(rco.dphi/(2*57.296));
	      rco.csi_yg=rco.csi_gcc*TMath::Cos(rco.dphi/(2*57.296));
	      rco.csi_yh=rco.csi_yg;
	      rco.csi_xh=-rco.csi_xg;
	      rco.csi_yH=rco.csi_yh;
	      rco.csi_yG=rco.csi_yg;	      
	      rco.csi_xG=rco.spe_morto_csi;
	      rco.csi_xH=-(float)(rco.spe_morto_csi);	      

	      if(rco.csi[0]==0)
		{
		  rco.csi[0]=new TGraph();
		  rco.csi[0]->SetPoint(0,rco.csi_xb,rco.csi_yb);
		  rco.csi[0]->SetPoint(1,rco.csi_xB,rco.csi_yB);
		  rco.csi[0]->SetPoint(2,rco.csi_xD,rco.csi_yD+rco.spe_morto_csi);
		  rco.csi[0]->SetPoint(3,rco.csi_xd,rco.csi_yd+rco.spe_morto_csi);
		}
	      if(rco.csi[1]==0)
		{
		  rco.csi[1]=new TGraph();
		  rco.csi[1]->SetPoint(0,rco.csi_xa,rco.csi_ya);
		  rco.csi[1]->SetPoint(1,rco.csi_xA,rco.csi_yA);
		  rco.csi[1]->SetPoint(2,rco.csi_xC,rco.csi_yC+rco.spe_morto_csi);
		  rco.csi[1]->SetPoint(3,rco.csi_xc,rco.csi_yc+rco.spe_morto_csi);
		}
	      if(rco.csi[2]==0)
		{
		  rco.csi[2]=new TGraph();
		  rco.csi[2]->SetPoint(0,rco.csi_xd,rco.csi_yd-rco.spe_morto_csi);
		  rco.csi[2]->SetPoint(1,rco.csi_xD,rco.csi_yD-rco.spe_morto_csi);
		  rco.csi[2]->SetPoint(2,rco.csi_xF,rco.csi_yF+rco.spe_morto_csi);
		  rco.csi[2]->SetPoint(3,rco.csi_xf,rco.csi_yf+rco.spe_morto_csi);
		}
	      if(rco.csi[3]==0)
		{
		  rco.csi[3]=new TGraph();
		  rco.csi[3]->SetPoint(0,rco.csi_xc,rco.csi_yc-rco.spe_morto_csi);
		  rco.csi[3]->SetPoint(1,rco.csi_xC,rco.csi_yC-rco.spe_morto_csi);
		  rco.csi[3]->SetPoint(2,rco.csi_xE,rco.csi_yE+rco.spe_morto_csi);
		  rco.csi[3]->SetPoint(3,rco.csi_xe,rco.csi_ye+rco.spe_morto_csi);
		}
	      if(rco.csi[4]==0)
		{
		  rco.csi[4]=new TGraph();
		  rco.csi[4]->SetPoint(0,rco.csi_xf,rco.csi_yf-rco.spe_morto_csi);
		  rco.csi[4]->SetPoint(1,rco.csi_xF,rco.csi_yF-rco.spe_morto_csi);
		  rco.csi[4]->SetPoint(2,rco.csi_xH,rco.csi_yH+rco.spe_morto_csi);
		  rco.csi[4]->SetPoint(3,rco.csi_xh,rco.csi_yh+rco.spe_morto_csi);
		}
	      if(rco.csi[5]==0)
		{
		  rco.csi[5]=new TGraph();
		  rco.csi[5]->SetPoint(0,rco.csi_xe,rco.csi_ye-rco.spe_morto_csi);
		  rco.csi[5]->SetPoint(1,rco.csi_xE,rco.csi_yE-rco.spe_morto_csi);
		  rco.csi[5]->SetPoint(2,rco.csi_xG,rco.csi_yG+rco.spe_morto_csi);
		  rco.csi[5]->SetPoint(3,rco.csi_xg,rco.csi_yg+rco.spe_morto_csi);
		}

	      for(int j=0;j<6;j++)
		{
		  double a,b;
		  rco.csi[j]->GetPoint(0,a,b);
		  rco.csi[j]->SetPoint(rco.csi[j]->GetN(),a,b);
		}
	      int totstriprotte[8];
	      rco.stripcopertamax=-1;
	      for(int j=7;j>=0;j--)
		{
		  totstriprotte[j]=0;
		  for(int k=0;k<8;k++)
		    {
		      if(rco.codice_strip[j][k]>0)
			{
			  totstriprotte[j]=totstriprotte[j]+1;
			}
		    }
		  if(totstriprotte[j]==8)
		    {
		      rco.stripcopertamax=j;
		    }
		}
	      


    }
 return;
}

int Classe_geo::InsideRCO(float thetapart,float phipart,float zpart,int *code)
{
  InizializzaRCO();

  //<<<<<<< Classe_geo.cxx
  int inside=0;
  *code=-1;
  int isec0=-1;
  int istrip0=-1;
  int icsi0=-1;

  if(rco.stripcopertamax>=0)//fino alla strip coperta non si vede niente (c'e' uno schermo di alluminio)
    {
      if(thetapart<rco.strip_themax[0][rco.stripcopertamax])
	{
	  inside=0;
      return inside;
	}
    }

  for(int isec=0;isec<8;isec++)
    {
      if(rco.gas_the[isec]>0)
	{
	  if(thetapart>=rco.gas_themin[isec] && thetapart<rco.gas_themax[isec])
	    {
	  if(rco.gas_phimin[isec]*rco.gas_phimax[isec]>=0)
	    {
	      //if(phipart>=rco.gas_phimin[isec]+rco.dphimorto && phipart<rco.gas_phimax[isec]-rco.dphimorto)
if(phipart>=rco.gas_phimin[isec] && phipart<rco.gas_phimax[isec])
	    {



      isec0=isec;
      break;


	}
	    }//phimin*phimax>=0
	  else
	    {
	      if(rco.gas_phimin[isec]>0 && rco.gas_phimax[isec]<0)
		{

		  if((phipart>=rco.gas_phimin[isec] && phipart<rco.gas_phimin[isec]+rco.dphi)||(phipart<rco.gas_phimax[isec] && phipart>=rco.gas_phimax[isec]-rco.dphi))
		    {
	  

      isec0=isec;
      break;

		    }


		}//phimin>0, phimax<0
	      if(rco.gas_phimin[isec]<0 && rco.gas_phimax[isec]>0)
		{
		
if(phipart>=rco.gas_phimin[isec] && phipart<rco.gas_phimax[isec])
	    {
	      


      isec0=isec;
      break;


	  
	    }
		}//phimin<0 phimax>0

	    }//phimin*phimax<0
	
	    }//thetapart=themin-themax
	}//rco_the>0
    }//isec 0-7

  if(isec0<0) // SE NON SIAMO DENTRO IL GAS, SI BUTTA VIA TUTTO
    {
      inside=0;
      return inside;
    }
  
  if(isec0>=0)//si cercano le coordinate della particella sul piano del rivelatore
    {
      //ci va pi/2-phi perche' il phi parte dall'asse y (quindi l'angolo rispetto a x e' pi/2-phi)
	      float punto[3];
	      float t=(pow(rco.gas_veccentro[0],2)+pow(rco.gas_veccentro[1],2)+pow(rco.gas_veccentro[2],2))/(rco.gas_veccentro[0]*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296)+rco.gas_veccentro[1]*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296)+rco.gas_veccentro[2]*cos(thetapart/57.296));

	      punto[0]=t*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296);
	      punto[1]=t*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296);
	      punto[2]=t*cos(thetapart/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-rco.gas_veccentro[nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_c,rco.gas_vi[isec0],&pi);
	      Classe_formule::scaprod(punto_c,rco.gas_vj[isec0],&pj);
	      float yrif=rco.gas_ya/rco.gas_xa*pi+rco.dx_spe*sqrt(1+pow(rco.gas_ya/rco.gas_xa,2));
	      float xb=-rco.gas_xa;
	      float yrifb=rco.gas_ya/xb*pi+rco.sx_spe*sqrt(1+pow(rco.gas_ya/xb,2));

	      if(pi>=0 &&pj<yrif)
		{
		  inside=0;
		  return inside;
		}
	      if(pi<0&&pj<yrifb)
		{
		  inside=0;
		  return inside;
		}

	      if(pj<rco.gas_yg)//la strip 8 in basso e' dritta (e non tonda)
		{
		  inside=0;
		  return inside;
		}
	      *code=(isec0+1)*100;
	      inside=1;

	      for(int istrip=0;istrip<7;istrip++)
		{
		  if(thetapart>=rco.strip_themin[isec0][istrip]&& thetapart<=rco.strip_themax[isec0][istrip])
		    {
		      istrip0=istrip;

		      break;
		    }
		}
	      if(istrip0<0)
		{

		  if(thetapart<=rco.strip_themax[isec0][7]) 
		{
	      float ts=(pow(rco.strip_veccentro[0],2)+pow(rco.strip_veccentro[1],2)+pow(rco.strip_veccentro[2],2))/(rco.strip_veccentro[0]*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296)+rco.strip_veccentro[1]*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296)+rco.strip_veccentro[2]*cos(thetapart/57.296));

	      float puntos[3];
	      puntos[0]=ts*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296);
	      puntos[1]=ts*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296);
	      puntos[2]=ts*cos(thetapart/57.296);
	      float punto_cs[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_cs[nn]=puntos[nn]-rco.strip_veccentro[nn];
		}
	      float pis,pjs;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_cs,rco.gas_vi[isec0],&pis);
	      Classe_formule::scaprod(punto_cs,rco.gas_vj[isec0],&pjs);
	     
	     
	      if(pjs>=rco.strip_yg) // va verificato se ci sono casi in cui questa condizione sarebbe verificata ma pj>ybb non lo e' (in questo caso si perde potenzialmente roba che potrebbe essere identificata da Si-CsI
		    {
		      istrip0=7;
		         
		    }
		}
		}
	      if(istrip0>=0)
		{
		  *code=*code+(istrip0+1)*10;

		}

	      float tcsi=(pow(rco.csi_veccentro[0],2)+pow(rco.csi_veccentro[1],2)+pow(rco.csi_veccentro[2],2))/(rco.csi_veccentro[0]*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296)+rco.csi_veccentro[1]*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296)+rco.csi_veccentro[2]*cos(thetapart/57.296));
	      float puntocsi[3];
	      puntocsi[0]=tcsi*sin(thetapart/57.296)*cos(TMath::Pi()/2-phipart/57.296);
	      puntocsi[1]=tcsi*sin(thetapart/57.296)*sin(TMath::Pi()/2-phipart/57.296);
	      puntocsi[2]=tcsi*cos(thetapart/57.296);
	      float punto_ccsi[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_ccsi[nn]=puntocsi[nn]-rco.csi_veccentro[nn];
		}
	      float picsi,pjcsi;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vi[isec0],&picsi);
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vj[isec0],&pjcsi);


    
  for(int icsi=0;icsi<6;icsi++)
    {
      if(rco.csi[icsi]->IsInside(picsi,pjcsi)==1)
	{
	  icsi0=icsi;
	  break;
	}
    }
  if(icsi0>=0)
    {
      *code=*code+(icsi0+1);
    }  


    }
  return inside;

}
int Classe_geo::Soglie_Phos(float zpart, float apart, float thetapart, float epart,float tvolo,float codephos,float *zout,float *aout,float *tout)
{
  //ritorna 3 se la roba arriva in cesio; ritorna 2 se si ferma nel secondo plastico; ritorna 1 se si ferma nel primo plastico; ritorna 0 se si ferma nel target o nel mylar di ingresso dei phoswich

  int sottosoglia=3;
  if(TMath::Nint(zpart)==0)//i neutroni si rivelano ma solo come conteggi
    {
      *zout=0;
      *aout=1;
      *tout=-1;
      return sottosoglia;
    }

  float dummy_pressione=0;
  int idir=1;
 int icod;
  float atloc=(float)target.at;
  float spesst=target.spessore/(2.*cos(thetapart/57.296)*1000*target.densita);
  float el1,e1;
  ecorr_veda(&epart,&zpart,&apart,&atloc,&e1,&el1,&target.materiale,&spesst,&idir,&icod,&dummy_pressione);
  if(e1<=0)
    {
      sottosoglia=0;
      return sottosoglia;
    }
  float emylar1;//energia residua dopo mylar1
  float elmylar1;//energia persa in mylar1
  int matemylar=16;
  atloc=0;
  float spesmylar=1.5;//micron
  spesmylar=spesmylar/TMath::Cos(thetapart/57.296);
  ecorr_veda(&e1,&zpart,&apart,&atloc,&emylar1,&elmylar1,&matemylar,&spesmylar,&idir,&icod,&dummy_pressione);
  if(emylar1<=0)
    {
      sottosoglia=0;
      return sottosoglia;
    }

  int mate=3;
  int ip=codephos/10-1;
  int ih=codephos-(ip+1)*10-1;
  float aplastico=6.5095;
  float spess=phos.spess[ip][ih];
  float spess2=5000;//micron spessore plastico spesso
  spess=spess/TMath::Cos(thetapart/57.296);
  spess2=spess2/TMath::Cos(thetapart/57.296);
  float epostplastico1,epostplastico2,el2;
 
 ecorr_veda(&emylar1,&zpart,&apart,&aplastico,&epostplastico1,&el1,&mate,&spess,&idir,&icod,&dummy_pressione);
 int stop1=0;
 int stop2=0;

 float vpostmezzotarget=Classe_formule::cluce*TMath::Sqrt(2*e1/(Classe_formule::amu*apart));
 float tpostmezzotarget=phos.dist[ip][ih]/vpostmezzotarget;
 
if(epostplastico1<=0)
   {
     stop1=1;
     sottosoglia=1;
   }
 else
   {
 ecorr_veda(&epostplastico1,&zpart,&apart,&aplastico,&epostplastico2,&el2,&mate,&spess2,&idir,&icod,&dummy_pressione);
 if(epostplastico2<=0)
   {
     stop2=1;
     sottosoglia=2;
   }
   }
 if(stop1==1)
   {
     //     *tout=tvolo+gRandom->Gaus(0.,1.);
*tout=tpostmezzotarget+gRandom->Gaus(0.,1.);
	 // si ha il tempo di volo
     if(TMath::Nint(zpart)==1)
       {
	 *zout=zpart;
	 *aout=1;
       }
     if(TMath::Nint(zpart)==2)
       {
	 *zout=zpart;
	 *aout=4;
       }
     if(TMath::Nint(zpart)>12)//modificato per pulire gli eventi ER da fissione come suggerito da MC 9-11-11
       {
	 // frammenti heavy
	 //  {
	     *zout=50;
	     *aout=50;
	     //  }
		 double EBeam=Classe_analisi::Getanalisi()->reazione.ebeam;
		 double tMin,tMax;
		 if(EBeam<7.5)tMin=70;
		 if(EBeam>8. && EBeam<11)tMin=60;
		 if(EBeam>11)tMin=50;
		 if(EBeam<7.5)tMax=120;
		 if(EBeam>8. && EBeam<11)tMax=95;
		 if(EBeam>11)tMax=70;

	     //	     if(tvolo>=50 && tvolo<=70&&TMath::Nint(zpart)>=20)//regione attesa per residuo di fusione (ma sicuramente inquinato da ff asimmetriche) ho aggiunto la condizione Z>25 solo sul residuo ripristinando z>12 per gli "z50" stoppati in primo plastico perche con z>25 non c'era mai la fiss simmetrica :) visto che ztot=42. 10-11-11
		 if((tvolo>=tMin && tvolo <=tMax && TMath::Nint(zpart)>=20)&&(codephos<50))
	       {
		 //		 if(codephos>50)// printf("codephos=%f zpart=%f\n",codephos,zpart);
	     *zout=100;
	     *aout=100;
	   }
       }
     if(TMath::Nint(zpart)<=12 && TMath::Nint(zpart)>2)
       {
	 *zout=50; 
	 *aout=50; //messo 50 invece che -1 per uniformarsi al taglio triangolare 50 andante messo per i FF "simmetrici" 21-11-11
       }

   }

 if(stop2==1)
   {
     //     *tout=tvolo+gRandom->Gaus(0.,1.);
*tout=tpostmezzotarget+gRandom->Gaus(0.,1.);
     //si ha il tempo di volo
    
    

*zout=zpart;

     if(TMath::Nint(zpart)==1)
       {
	 *aout=apart;

       }     
    

	 if(TMath::Nint(zpart)==2)
	   {
	     *aout=4;
	   }
	   if(TMath::Nint(zpart)==3)
	   {
	     *aout=7;
	   }
	 if(TMath::Nint(zpart)==4)
	   {
	     *aout=7;
	   }
	 if(TMath::Nint(zpart)==5)
	   {
	     *aout=10;
	   }
	 if(TMath::Nint(zpart)>=6 &&TMath::Nint(zpart)<=12)
	   {
	 *aout=2*zpart;
	   }
	 if(TMath::Nint(zpart)>12)
	   {
	     *aout=Classe_formule::EAL(zpart);
	     
	   }
     }


 if(stop1==0&&stop2==0)// id in cesio
   {
     //     *tout=-1;
     //non c'e' tempo di volo
     //     *tout=tvolo+gRandom->Gaus(0.,1.);
     *tout=tpostmezzotarget+gRandom->Gaus(0.,1.);


     *zout=zpart;
     if(TMath::Nint(zpart)==1)
       {
	 *aout=apart;
       }     
     else
       {
	 if(TMath::Nint(zpart)==2)
	   {
	     *aout=4;
	   }
	   if(TMath::Nint(zpart)==3)
	   {
	     *aout=7;
	   }
	 if(TMath::Nint(zpart)==4)
	   {
	     *aout=7;
	   }
	 if(TMath::Nint(zpart)==5)
	   {
	     *aout=10;
	   }
	 if(TMath::Nint(zpart)>=6 &&TMath::Nint(zpart)<=12)
	   {
	 *aout=2*zpart;
	   }
	 else
	   {
	     *aout=Classe_formule::EAL(zpart);
	     
	   }
       }
       }     

 //phos per i quali non si corregge i tvolo per p e/o alfa:

 // if(ip+1==2 && ih+1==4 && (zpart==1 ||zpart==2)){*tout=-1;}
//  if(ip+1==3 && ih+1==6 && (zpart==1 ||zpart==2)){*tout=-1;}
//  if(ip+1==5 && (zpart==1 ||zpart==2)){*tout=-1;}
//  if(ip+1==6 && (zpart==1 ||zpart==2)){*tout=-1;}
//  if(ip+1==1 && ih+1==2 && zpart==1){*tout=-1;}
//  if(ip+1==1 && ih+1==6 && zpart==1){*tout=-1;}
//  if(ip+1==1 && ih+1==8 && zpart==1){*tout=-1;}
//  if(ip+1==1 && ih+1==9 && zpart==1){*tout=-1;}
//  if(ip+1==2 && ih+1==3 && zpart==1){*tout=-1;}
//  if(ip+1==2 && ih+1==8 && zpart==1){*tout=-1;}
//  if(ip+1==3 && ih+1==2 && zpart==1){*tout=-1;}
//  if(ip+1==3 && ih+1==7 && zpart==1){*tout=-1;}
//  if(ip+1==3 && ih+1==8 && zpart==1){*tout=-1;}
//  if(ip+1==4 && ih+1==3 && zpart==1){*tout=-1;}
//  if(ip+1==4 && ih+1==5 && zpart==1){*tout=-1;}
//  if(ip+1==4 && ih+1==7 && zpart==1){*tout=-1;}
//  if(ip+1==4 && ih+1==9 && zpart==1){*tout=-1;}

 return sottosoglia;
}



int Classe_geo::Soglie_Fazietto(float zpart, float apart, float thetapart, float epart, int code,float *zout, float *aout, float *eout, int *aid)
{
  *zout=-1;
  *aout=-1;
  *eout=-1;
  *aid=0;
  int zll;
  int zaa;
  int rivelato=-1;
  //0 roba non identificata (per es primi 30 micron del silicio 1)
  //11 roba stoppata in Si1; 12 roba stoppata in Si2; 23 roba stoppata in CsI;
  if(TMath::Nint(zpart)==0 && TMath::Nint(apart)==1)
    {
      *zout=0;
      *aout=0;
      *eout=-1;
      rivelato=33;
      return rivelato;
    }
  int iblocco=code/100;//parte da 0
  int  iqua=(code-iblocco*100)/4; // parte da 0
  int itel=(code-iblocco*100)-iqua*4; // parte da 0
 //cout<<code<<" " <<iblocco<<" "<<iqua<<" "<<itel<<endl;
 
  float dummy_pressione=0;
  int idir=1;
 int icod;
int idirt=2;
       float elostt;
  float atloct=(float)target.at;
float spesst=target.spessore/(2.*cos(thetapart/57.296)*1000*target.densita);
 float el1,e1;
  ecorr_veda(&epart,&zpart,&apart,&atloct,&e1,&el1,&target.materiale,&spesst,&idir,&icod,&dummy_pressione);
  if(e1<=0)//roba che muore nel target
    {
      rivelato=-1;
      return rivelato;
    }
int  matesi=1;
float atloc=0;
float spes_si1=fazietto.spes_si1[iblocco][iqua][itel]/TMath::Cos(thetapart/57.296);
 float epostsi1,esi1;
ecorr_veda(&e1,&zpart,&apart,&atloc,&epostsi1,&esi1,&matesi,&spes_si1,&idir,&icod,&dummy_pressione);
 if(epostsi1<=0)
   {

     //     float range=Range(zpart,apart,e1,matesi);
     // range=range*TMath::Cos(thetapart/57.296);
     //if(range>50)//micron
     if((fazietto.rotti[iblocco][iqua][itel]&4)!=0)
       {

     if(fazietto.zesi1inf[iblocco][iqua][itel][(int)zpart]>0 && e1>fazietto.zesi1inf[iblocco][iqua][itel][(int)zpart])
       {

     *zout=zpart;
     zaa=1000;
     zll=-1000;
     float eainf=1000;
     float easup=0;
       if(fazietto.aliminfSi1PSA[iblocco][iqua][itel]>0)
	 {
	   zaa=fazietto.aliminfSi1PSA[iblocco][iqua][itel];
	 }
       if(fazietto.alimsupSi1PSA[iblocco][iqua][itel]>0)
	 {
	   zll=fazietto.alimsupSi1PSA[iblocco][iqua][itel];
	 }
       if(fazietto.aesi1inf[iblocco][iqua][itel][(int)zpart]>0)
	 {
	   eainf=fazietto.aesi1inf[iblocco][iqua][itel][(int)zpart];
	 }
       if(fazietto.aesi1sup[iblocco][iqua][itel][(int)zpart]>0)
	 {
	   easup=fazietto.aesi1sup[iblocco][iqua][itel][(int)zpart];
	 }
       if(zpart>=zaa && zpart<=zll)
	 {
	   if(apart<fazietto.infmassepsa[iblocco][iqua][itel][(int)zpart]||apart > fazietto.supmassepsa[iblocco][iqua][itel][(int)zpart])
	     {
	       //     cout<<"qui"<<" "<<zpart<<" "<<apart<<" "<<iblocco<<" "<<iqua<<" "<<itel<<endl;
     *zout=-1;
     *aout=-1;
     *eout=-1;
	       return 0;
	     }
	 }
     if(zpart<=zll&&zpart>=zaa && e1>=eainf && e1<=easup)
       {
	 *aout=apart;
	 *aid=1;
       }
	 else
	   {
     *aout=Classe_formule::EAL(zpart);
     *aid=0;
	   }
     // *eout=esi1+gRandom->Gaus(0,0.005*esi1);
     
  
     esi1=esi1+gRandom->Gaus(0,0.006*esi1);
     ecorr_veda(&esi1,&zpart,&apart,&atloct,eout,&elostt,&target.materiale,&spesst,&idirt,&icod,&dummy_pressione);//L'energia in uscita tiene conto delle perdite in mezzo target (idirt=2, ecorr_veda usata al rovescio), quindi e' esi1 + qualcosa. eout e' il risultato
  
     rivelato=11;//PSA in Si1

     return rivelato;
       }
     else
       {
	 *zout=-1;
	 *aout=-1;
	 *aid=0;
	 *eout=esi1+gRandom->Gaus(0,0.006*esi1);
	 return 0;	 
       }
       }
     else //Se si1 e' rotto
       {
	 *zout=-1;
	 *aout=-1;
	 *aid=0;
	 *eout=-1;
	 return 0;
       }
   }//epostsi1<=0
 //spessore morto di 1 micron fra Si1 e Si2
 float spemorto=1./TMath::Cos(thetapart/57.296);
 float epostmorto,emorto;
ecorr_veda(&epostsi1,&zpart,&apart,&atloc,&epostmorto,&emorto,&matesi,&spemorto,&idir,&icod,&dummy_pressione);
 if(epostmorto<=0) //roba che muore nello strato morto fra si1 e si2
   {
     *zout=-1;
     *aout=-1;
     *eout=-1;
     return 0;
   }
float spes_si2=fazietto.spes_si2[iblocco][iqua][itel]/TMath::Cos(thetapart/57.296);
 float epostsi2,esi2;
ecorr_veda(&epostmorto,&zpart,&apart,&atloc,&epostsi2,&esi2,&matesi,&spes_si2,&idir,&icod,&dummy_pressione);
//ecorr_veda(&epostsi1,&zpart,&apart,&atloc,&epostsi2,&esi2,&matesi,&spes_si2,&idir,&icod,&dummy_pressione);
 float esi2_0=esi2;
 if(epostsi2<=0)
   {
     if(esi2<0.2)//soglia in si2
       {
     *zout=-1;
     *aout=-1;
     *eout=-1;
     return 0;
       }
     if((fazietto.rotti[iblocco][iqua][itel]&4)!=0 && (fazietto.rotti[iblocco][iqua][itel]&2)!=0)
       {
     if(fazietto.alimSi1Si2[iblocco][iqua][itel]<0)
       {
     *zout=-1;
     *aout=-1;
     *eout=-1;
     return 0;
       }
     *zout=zpart;
     zll=fazietto.alimSi1Si2[iblocco][iqua][itel];

     if(zpart<=zll)
       {
	 if(apart>fazietto.supmassesi1si2[iblocco][iqua][itel][(int)zpart]||apart<fazietto.infmassesi1si2[iblocco][iqua][itel][(int)zpart])
	   {
	     *zout=-1;
	     *aout=-1;
	     *eout=-1;
	     //cout<<"qui"<<" "<<zpart<<" "<<apart<<" "<<iblocco<<" "<<iqua<<" "<<itel<<endl;
	     return 0;
	   }
	 *aout=apart;
	 *aid=1;
       }
     else
       {
	 *aout=Classe_formule::EAL(zpart);
	 *aid=0;
       }
     //*eout=esi1+esi2+gRandom->Gaus(0,0.006*esi1)+gRandom->Gaus(0,0.007*esi2);
      esi1=esi1+gRandom->Gaus(0,0.006*esi1);
      esi2=esi2+gRandom->Gaus(0,0.007*esi2);
      float esi12=esi1+esi2;
     ecorr_veda(&esi12,&zpart,&apart,&atloct,eout,&elostt,&target.materiale,&spesst,&idirt,&icod,&dummy_pressione);//L'energia in uscita tiene conto delle perdite in mezzo target (idirt=2, ecorr_veda usata al rovescio), quindi e' esi1+esi2 + qualcosa. eout e' il risultato

 rivelato=12;//da Si1-Si2
     return rivelato;
   }
     else
       {
     *zout=-1;
     *aout=-1;
     *eout=-1;
     return 0;
       }
   }
 //si arriva in cesio
 rivelato=23;//da Si2-CsI

if((fazietto.rotti[iblocco][iqua][itel]&2)==0 && (fazietto.rotti[iblocco][iqua][itel]&1)==0)
  {
	*zout=-1;
	*aout=-1;
	*eout=-1;
	return 0;
  }

    if((fazietto.rotti[iblocco][iqua][itel]&1)!=0 && (fazietto.rotti[iblocco][iqua][itel]&2)!=0)
      {
     if(fazietto.alimSi2CsI[iblocco][iqua][itel]<0)
       {
	 if(zpart>=5)
	   {
	*zout=-1;
	*aout=-1;
	*eout=-1;
	return 0;
	   }

       }
     if(apart<fazietto.infmassesi2csi[iblocco][iqua][itel][(int)zpart]||apart>fazietto.supmassesi2csi[iblocco][iqua][itel][(int)zpart])
       {
	*zout=-1;
	*aout=-1;
	*eout=-1;
	return 0;	 
       }

     *zout=zpart;
     zll=fazietto.alimSi2CsI[iblocco][iqua][itel];

     if(zpart<=zll)
       {
	 *aout=apart;
	 *aid=1;
       }
     else
       {
	 *aout=Classe_formule::EAL(zpart);
	 *aid=0;
       }
      }

    if((fazietto.rotti[iblocco][iqua][itel]&2)==0 && (fazietto.rotti[iblocco][iqua][itel]&1)!=0)
      {
	if(zpart<5)
	  {
	    if(TMath::Nint(zpart)==1 && (TMath::Nint(apart)<1 ||TMath::Nint(apart)>3))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==2 && (TMath::Nint(apart)<3 ||TMath::Nint(apart)>6))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==3 && (TMath::Nint(apart)<6 ||TMath::Nint(apart)>9))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==4 && (TMath::Nint(apart)<7 ||TMath::Nint(apart)>11))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }

	    *zout=zpart;
	    *aout=apart;
	    *aid=1;
	    rivelato=33;
	  }
	else
	  {
	*zout=-1;
	*aout=-1;
	*eout=-1;
	return 0;
	  }
      }
    //    if(*zout<5 && *zout>0 && esi2_0<1)
    if(*zout<5 && *zout>0 && (fazietto.alimSi2CsI[iblocco][iqua][itel]<0 || (esi2_0<1&&fazietto.alimSi2CsI[iblocco][iqua][itel]>0)))
      {
	    if(TMath::Nint(zpart)==1 && (TMath::Nint(apart)<1 ||TMath::Nint(apart)>3))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==2 && (TMath::Nint(apart)<3 ||TMath::Nint(apart)>6))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==3 && (TMath::Nint(apart)<6 ||TMath::Nint(apart)>9))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	    if(TMath::Nint(zpart)==4 && (TMath::Nint(apart)<7 ||TMath::Nint(apart)>11))
	      {
		*zout=-1;
		*aout=-1;
		*eout=-1;
		return 0;
	      }
	*aout=apart;
	*aid=1;
	rivelato=33;
      }
    

    //  if(iblocco==0 &&iqua==1 && itel==0)cout<<rivelato<<" "<<zpart<<" "<<esi2_0<<endl;
     float ecsi=epostsi2;
     float ett=ecsi+gRandom->Gaus(0,0.03*ecsi)+esi1+gRandom->Gaus(0,0.006*esi1)+esi2+gRandom->Gaus(0,0.007*esi2);
     //*eout=ett;
     ecorr_veda(&ett,&zpart,&apart,&atloct,eout,&elostt,&target.materiale,&spesst,&idirt,&icod,&dummy_pressione);//L'energia in uscita tiene conto delle perdite in mezzo target (idirt=2, ecorr_veda usata al rovescio), quindi e' esi1 + qualcosa. eout e' il risultato

     float elost=esi2+esi1+gRandom->Gaus(0,0.006*esi1)+gRandom->Gaus(0,0.007*esi2);
     float spessi=(fazietto.spes_si2[iblocco][iqua][itel]+fazietto.spes_si1[iblocco][iqua][itel])/TMath::Cos(thetapart/57.296);
//    *eout= DE2E(zpart,*aout,elost,matesi,spessi); //E totale da Delta E in Si1 e Si2
   

  return rivelato;
}
void Classe_geo::Spalma_Fazietto(float codice, float *theout, float *phiout, float *dout)
{
  int iblocco=codice/100;
  int iriv=codice-iblocco*100;

  //cout<<"spalma"<<iblocco<<" "<<iriv<<endl;
  float pi=fazietto.vx[iblocco][iriv]+(1-2*gRandom->Rndm(0))*fazietto.dx_active/2; // componente rispetto al versore vi nel piano del rivelatore
  float pj=fazietto.vy[iblocco][iriv]+(1-2*gRandom->Rndm(0))*fazietto.dy_active/2;// idem per vj
  float pnew[3];
  pnew[0]=pi*fazietto.viblocco[iblocco][0]+pj*fazietto.vjblocco[iblocco][0]+fazietto.vcentroblocco[iblocco][0];
  pnew[1]=pi*fazietto.viblocco[iblocco][1]+pj*fazietto.vjblocco[iblocco][1]+fazietto.vcentroblocco[iblocco][1];
  pnew[2]=pi*fazietto.viblocco[iblocco][2]+pj*fazietto.vjblocco[iblocco][2]+fazietto.vcentroblocco[iblocco][2];
  Classe_formule::carpol(dout,theout,phiout,pnew);

  return;
}


// da provare se ci sara' un mcarlo con gamma
int Classe_geo::Soglie_Hector(float zpart,float apart,float epart,float tvolo,float codehector,float *zout,float *aout,float *tout,float *eout)
{
  //rivelatore per gamma
  int sottosoglia=-1;
  *zout=-1;
  *aout=-1;
  *tout=-1;
  *eout=-1;
  if(TMath::Nint(zpart)>0)//roba carica
    {
      *zout=-1;
      *aout=-1;
      *tout=-1;
      *eout=-1;
      return sottosoglia;
    }
  if(TMath::Nint(apart)>0)//neutroni
    {
      *zout=-1;
      *aout=-1;
      *tout=-1;
      *eout=-1;
      return sottosoglia;
    }
  //  if(TMath::Nint(zpart)==0 && TMath::Nint(apart)==0 && epart>0.08)
  if(TMath::Nint(zpart)==0 && TMath::Nint(apart)==0 && epart>4.5)
    {
      *zout=0;
      *aout=0;
      sottosoglia=1;
 //     int ip=codehector-1;    a cosa serve????
      *tout=tvolo;//da sporcare con risoluz temporale
      *eout=epart;//da sporcare con risoluz energetica
    }

  return sottosoglia;

}


void Classe_geo::Spalma_Garfield(float codice, float *theout, float *phiout,int code_micro)
{
  int isec=codice/10;
  int icsi=codice-isec*10;
  int dentro=-1;

  if(garf.dist[icsi-1]<=0)
    {
      cout<<"cesio "<<icsi<<" non contenuto nella geometria attualmente caricata"<<endl;
      *theout=-1000;
      *phiout=-1000;
      exit(1);

    }

  while(dentro<0)
    {
      float pi,pj;
      float pk=0;
      pj=-garf.h[icsi-1]/2+garf.hdown[icsi-1]+(garf.h[icsi-1]/2-garf.hup[icsi-1]+garf.h[icsi-1]/2-garf.hdown[icsi-1])*gRandom->Rndm(0);
      float pimin,pimax;
      if(garf.l1[icsi-1]>garf.l2[icsi-1]) 
	{
	  pimin=-garf.l1[icsi-1]/2+garf.spess_fasciatura;
	  pimax=garf.l1[icsi-1]/2-garf.spess_fasciatura;
	}
      else
	{
	  pimin=-garf.l2[icsi-1]/2+garf.spess_fasciatura;
	  pimax=garf.l2[icsi-1]/2-garf.spess_fasciatura;
	}
      if(code_micro==1)//solo sx
	{
	   pimax=0;
	  
	}
      if(code_micro==2)//solo dx
	{
	  pimin=0;
	  
	}
      pi=pimin+(pimax-pimin)*gRandom->Rndm(0);
	      if(garf.gcsi[icsi-1]->IsInside(pi,pj)==1)
		{


		  dentro=1;

		  float v[3];
	      v[0]=pi*garf.vi[isec-1][icsi-1][0]+pj*garf.vj[isec-1][icsi-1][0]+pk*garf.vk[isec-1][icsi-1][0];
	      v[1]=pi*garf.vi[isec-1][icsi-1][1]+pj*garf.vj[isec-1][icsi-1][1]+pk*garf.vk[isec-1][icsi-1][1];
	      v[2]=pi*garf.vi[isec-1][icsi-1][2]+pj*garf.vj[isec-1][icsi-1][2]+pk*garf.vk[isec-1][icsi-1][2];
	      for(int j=0;j<3;j++)
		{
		  v[j]=v[j]+garf.veccentro[isec-1][icsi-1][j];
		}
	      float r;
	      Classe_formule::carpol(&r,theout,phiout,v);
	

	      return;
		}
    }  
 




//   *theout=garf.themincsi[icsi-1]+(garf.themaxcsi[icsi-1]-garf.themincsi[icsi-1])*gRandom->Rndm(0);
 

//   *phiout=garf.phi[icsi-1][isec-1]-garf.dphi_corr/2+gRandom->Rndm(0)*garf.dphi_corr;

//   if(code_micro==301)
//     {
// *phiout=garf.phi[icsi-1][isec-1]-garf.dphi_corr/2*gRandom->Rndm(0);
//     }
//   if(code_micro==302)
//     {
// *phiout=garf.phi[icsi-1][isec-1]+garf.dphi_corr/2*gRandom->Rndm(0);
//     }

}
void Classe_geo::Spalma_Phos(float codice, float *theout, float *phiout, float *dout)
{
 int ip=codice/10-1;
  int ih=codice-(ip+1)*10-1;
  float pi=(1-2*gRandom->Rndm(0))*phos.dx/2; // componente rispetto al versore vi nel piano del rivelatore
  float pj=(1-2*gRandom->Rndm(0))*phos.dx/2;// idem per vj
  float pnew[3];
  pnew[0]=pi*phos.vi[ip][ih][0]+pj*phos.vj[ip][ih][0]+phos.veccentro[ip][ih][0];
  pnew[1]=pi*phos.vi[ip][ih][1]+pj*phos.vj[ip][ih][1]+phos.veccentro[ip][ih][1];
  pnew[2]=pi*phos.vi[ip][ih][2]+pj*phos.vj[ip][ih][2]+phos.veccentro[ip][ih][2];
  Classe_formule::carpol(dout,theout,phiout,pnew);
  return;
}

void Classe_geo::Spalma_Hector(int codice,float *theout, float *phiout,float *dout)
{
  int ip=codice-1;
  float ri=gRandom->Rndm(0)*hector.raggio;
  float ther=gRandom->Rndm(0)*360;
  float pi=ri*TMath::Cos(ther/57.296);
  float pj=ri*TMath::Sin(ther/57.296);
  float pnew[3];
  pnew[0]=pi*hector.vi[ip][0]+pj*hector.vj[ip][0]+hector.veccentro[ip][0];
  pnew[1]=pi*hector.vi[ip][1]+pj*hector.vj[ip][1]+hector.veccentro[ip][1];
  pnew[2]=pi*hector.vi[ip][2]+pj*hector.vj[ip][2]+hector.veccentro[ip][2];
  Classe_formule::carpol(dout,theout,phiout,pnew);

}


int Classe_geo::calcolaGARF(int z,int a,int icsi,double *ein,double *e1,double *e2,double *e3,double *loo) {
	int np=0,ie;
	double ei,emin,emax,estp,spess,epass,eprima,lo;
	double alpha,calp;
	
	double s_tar=target.spessore*10/target.mat.r; //µm
	double s_m1=garf.mylar1;      //µm
	double s_d1=garf.dead1;   //µm
	double s_s1=garf.strip1;   //µm
	double s_d2=garf.dead2;    //µm
	double s_s2=garf.strip2;   //µm
	double s_br=garf.braccio;   //µm
	double s_m2=garf.mylar2;     //µm
	double d_d3;
	
	MATE mylar={"mylar",3,{{6,-1,8},{1,-1,10},{8,-1,4}},1395.,-1.,-1.};
	MATE cf4  ={"CF4",2,{{6,-1,1},{9,-1,4}},-1,garf.pressione,300.};
	if(icsi<5) cf4.P=garf.pressione_back;
	//printf("Z=%d A=%d CsI%d",z,a,icsi+5);
	
	//emin=EMIN*(double)a; estp=ESTP*(double)a;
	double eb=Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.ebeam;
	double Ma=Classe_analisi::Getanalisi()->reazione.ap*931.5;
	double Mb=Classe_analisi::Getanalisi()->reazione.at*931.5;
	double Mtot=Ma+Mb;
	double Mc=931.5*(double)a;
	double Md=Mtot-Mc;
	emin=EMIN; estp=ESTPN*(double)a;
	if(estp>ESTP) estp=ESTP;
	if(emin<estp) emin=estp;
	emax=eb*(Mb*eb+Ma*Mc+Mb*Md+sqrt((eb+2.*Ma)*(Mb*Mb*eb+2.*Mb*Mc*Md)))/(Mtot*Mtot+2.*Mb*eb);
	emax*=1.3; //si sta un po' larghi per sicurezza (nel caso in cui si selezioni l'isotopo sbagliato)
	if((emax-emin)/estp>(double)NMAX) {
		estp=ceil(10.*(emax-emin)/(double)NMAX)/10.;
		//printf("z=%d a=%d icsi=%d => emin=%lf emax=%lf estp=%lf\n",z,a,icsi,emin,emax,estp);
	}
	for(ie=(int)(emin/estp);ie<(int)(emax/estp);ie++) {
		ei=estp*(double)ie;
		ein[np]=ei;
		
// 		if(z==12) {
// 		getchar();
// 		printf("\nStrato     Ein   Elost\n");
// 		printf("----------------------\n");
// 		}
		
		//MEZZO TARGET
		epass=eloss(ei,z,a,&(target.mat),fabs(s_tar/(2.*cos(garf.thecsi[icsi-1]/57.296))),SCHWALM);
// 		if(z==12) printf("½ target  %5.1lf %5.1lf\n",ei,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
        if(icsi>4) {alpha=90.-garf.alpha0-garf.thecsi[icsi-1];}
        else {alpha=garf.thecsi[icsi-1]-90.-garf.alpha0;}
		calp=1./cos(alpha/57.296);
		
		//MYLAR 1
// 		eprima=epass;
		epass=eloss(epass,z,a,&mylar,s_m1*calp,TABGARF);
// 		if(z==12) printf("Mylar 1   %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		//DEAD LAYER 1
// 		eprima=epass;
		epass=eloss(epass,z,a,&cf4,s_d1*calp,TABGARF);
// 		if(z==12) printf("Dead l 1  %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		//STRIP DOWN
		eprima=epass;
		epass=eloss(epass,z,a,&cf4,s_s1*calp,TABGARF);
// 		if(z==12) printf("Strip dn  %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		e1[np]=eprima-epass;
		
		//DEAD LAYER 2
// 		eprima=epass;
		epass=eloss(epass,z,a,&cf4,s_d2*calp,TABGARF);
// 		if(z==12) printf("Dead l 2  %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		//phi=alpha-(82.-ang[icsi]); //È sempre 0 perché prendo solo i centri dei cesi!
		d_d3=1000.*garf.dist[icsi-1]-(s_br+s_m1+s_d1+s_s1+s_d2+s_s2)*calp;
		
		//STRIP UP
		spess=s_s2*calp;
		if(d_d3<0) {
			spess+=d_d3;
		}
		eprima=epass;
		epass=eloss(epass,z,a,&cf4,spess,TABGARF);
// 		if(z==12) printf("Strip up  %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		e2[np]=eprima-epass;
		
		//DEAD LAYER 3, SE C'È
// 		eprima=epass;
		if(d_d3>0) {
			epass=eloss(epass,z,a,&cf4,d_d3,TABGARF);
// 			if(z==12) printf("Dead l 3  %5.1lf %5.1lf\n",eprima,epass);
			if((epass<=0)||(!isfinite(epass))) continue;
		}
		
		//MYLAR 2
// 		eprima=epass;
		epass=eloss(epass,z,a,&mylar,s_m2,TABGARF);
// 		if(z==12) printf("Mylar 2   %5.1lf %5.1lf\n",eprima,epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		if(epass<ETHR) continue;
		lo=E2Luce_csi((double)z,(double)a,epass);
		//if(lo>LMAX) break;
		e3[np]=epass;
		loo[np]=lo;
		np++;
		//if(np>=NMAX) {
		//	printf("[WARNING] Raggiunto il limite di punti in GARF per Z=%d A=%d CsI %d ad E=%lf\n",z,a,icsi,ei);
		//	break;
		//}
	}
	return np;
}


int Classe_geo::calcolaRCO(int z,int a,int isec,int istrip,double *ein,double *e1,double *e2,double *e3,double *loo,int *k) {
	int np=0,ie;
	double ei,emin,emax,estp,epass,eprima;
	double calp;
	
	double s_tar=target.spessore*10/target.mat.r; //µm
	double s_m1=rco.spesmylar;       //µm
	double s_gas=rco.spesgas;        //µm
	double s_m2=rco.spesmylar_mezzo; //µm
	double s_m3=rco.spesmylar;       //µm
	double s_si=rco.strip_spess[isec][istrip]; //µm
	double s_mcs=rco.spesmylarcsi;   //µm
	
	
	MATE mylar={"mylar",3,{{6,-1,8},{1,-1,10},{8,-1,4}},1395.,-1.,-1.};
	MATE cf4  ={"CF4",2,{{6,-1,1},{9,-1,4}},-1,rco.pressione,300.};
	MATE si   ={"Si",1,{{14,-1,1}},0,-1.,-1.};
	
	*k=0;
	//emin=EMIN*(double)a; estp=ESTP*(double)a;
	double eb=Classe_analisi::Getanalisi()->reazione.ap*Classe_analisi::Getanalisi()->reazione.ebeam;
	double Ma=Classe_analisi::Getanalisi()->reazione.ap*931.5;
	double Mb=Classe_analisi::Getanalisi()->reazione.at*931.5;
	double Mtot=Ma+Mb;
	double Mc=931.5*(double)a;
	double Md=Mtot-Mc;
	emin=EMIN; estp=ESTPN*(double)a;
	if(estp>ESTP) estp=ESTP;
	if(emin<estp) emin=estp;
	emax=eb*(Mb*eb+Ma*Mc+Mb*Md+sqrt((eb+2.*Ma)*(Mb*Mb*eb+2.*Mb*Mc*Md)))/(Mtot*Mtot+2.*Mb*eb);
	emax*=1.3; //si sta un po' larghi per sicurezza (nel caso in cui si selezioni l'isotopo sbagliato)
	if((emax-emin)/estp>(double)NMAX) {
		estp=ceil(10.*(emax-emin)/(double)NMAX)/10.;
		//printf("z=%d a=%d isec=%d istrip=%d => emin=%lf emax=%lf estp=%lf\n",z,a,isec,istrip,emin,emax,estp);
	}
	for(ie=(int)(emin/estp);ie<(int)(emax/estp);ie++) {
		ei=estp*(double)ie;
		ein[np]=ei;
		
		//if(z==7) {
		//getchar();
		//printf("\nStrato     Ein   Elost\n");
		//printf("----------------------\n");
		//}
		
		calp=1./cos(rco.strip_the[isec][istrip]/57.296);
		e1[np]=0;
		e2[np]=0;
		e3[np]=0;
		loo[np]=0;
		
		//MEZZO TARGET
		epass=eloss(ei,z,a,&(target.mat),s_tar*calp/2.,SCHWALM);
		//if(z==7) printf("½ target  %5.1lf %5.1lf\n",ei,ei-epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		//MYLAR 1
		//eprima=epass;
		epass=eloss(epass,z,a,&mylar,s_m1*calp,TABRCO);
		//if(z==7) printf("Mylar 1   %5.1lf %5.1lf\n",eprima,eprima-epass);
		if((epass<=0)||(!isfinite(epass))) continue;
		
		//GAS 1
		eprima=epass;
		epass=eloss(epass,z,a,&cf4,s_gas*calp,TABRCO);
		//if(z==7) printf("Mezzo gas %5.1lf %5.1lf\n",eprima,eprima-epass);
		e1[np]=eprima-epass;
		
		//MYLAR 2
		if((epass>0)&&isfinite(epass)) {
			//eprima=epass;
			epass=eloss(epass,z,a,&mylar,s_m2*calp,TABRCO);
			//if(z==7) printf("Mylar 2   %5.1lf %5.1lf\n",eprima,eprima-epass);
		}
		
		//GAS 2
		if((epass>0)&&isfinite(epass)) {
			eprima=epass;
			epass=eloss(epass,z,a,&cf4,s_gas*calp,TABRCO);
			//if(z==7) printf("Mezzo gas %5.1lf %5.1lf\n",eprima,eprima-epass);
			e1[np]+=eprima-epass;
		}
		
		if(((epass<=0)||(!isfinite(epass)))&&(e1[np]<ETHR)) continue;
		
		//MYLAR 3
		if((epass>0)&&isfinite(epass)) {
			//eprima=epass;
			epass=eloss(epass,z,a,&mylar,s_m3*calp,TABRCO);
			if(epass<ETHR) epass=0.;
			//if(z==7) printf("Mylar 3   %5.1lf %5.1lf\n",eprima,eprima-epass);
		}
		
		//SILICIO
		if((epass>0)&&isfinite(epass)) {
			eprima=epass;
			epass=eloss(epass,z,a,&si,s_si*calp,VEDALOSS);
			//if(z==7) printf("Silicio   %5.1lf %5.1lf\n",eprima,eprima-epass);
			e2[np]=eprima-epass;
			if(*k==0) *k=np+1;
		}
		
		//MYLAR CSI
		if((epass>0)&&isfinite(epass)) {
			//eprima=epass;
			epass=eloss(epass,z,a,&mylar,s_mcs*calp,TABRCO);
			if(epass<ETHR) epass=0.;
			//printf("Mylar 3   %5.1lf %5.1lf\n",eprima,epass);
			if(isfinite(epass)) e3[np]=epass;
			else e3[np]=0;
		}
		
		if((e1[np]<ETHR)&&(e2[np]<ETHR)&&(e3[np]<ETHR)) continue;
		
		if((epass>0)&&isfinite(epass)) {
			loo[np]=E2Luce_csi((double)z,(double)a,epass);
			if(*k<10000) *k+=(np+1)*10000;
		}
		//if(loo[np]>LMAX) break;
		np++;
		//if(np>=NMAX) {
		//	printf("[WARNING] Raggiunto il limite di punti in RCO per Z=%d A=%d Sec %d, strip %d ad E=%lf\n",z,a,isec,istrip,ei);
		//	break;
		//}
	}
	return np;
}


void scrivi_ric(FILE *f,nodo *radice,int l,int z,int a,int riv,int app) {
	int i;
	energie *en;
	
	fprintf(f,"%d %d\n",radice->ne,radice->min);
	if(radice->min<0) {
		fprintf(f,"# TABELLA PER Z=%d A=%d ",z,a);
		if(app) fprintf(f,"Sett %d, strip %d\n",1+riv/8,1+riv%8);
		else fprintf(f,"CsI %d\n",riv);
		en=(energie *)(radice->son);
		for(i=0;i<radice->ne;i++) {
			//fprintf(f,"%7.2lf %8.3lf %8.3lf %8.3lf %8.2lf\n",en[i].ei,en[i].e1,en[i].e2,en[i].e3,en[i].lo);
			//fprintf(f,"%7.2lf %8.4lf %8.4lf %8.4lf %8.2lf\n",en[i].ei,en[i].e1,en[i].e2,en[i].e3,en[i].lo);
			if(app) fprintf(f,"%7.2lf %7.4lf %8.3lf %8.3lf %6.0lf\n",en[i].ei,en[i].e1,en[i].e2,en[i].e3,en[i].lo);
			else fprintf(f,"%7.2lf %7.4lf %7.4lf %8.3lf %6.0lf\n",en[i].ei,en[i].e1,en[i].e2,en[i].e3,en[i].lo);
		}
	}
	else {
		for(i=0;i<radice->ne;i++) {
			switch(l) {
				case 0: z=i+radice->min; break;
				case 1: a=i+radice->min; break;
				case 2: riv=i+radice->min;
			}
			scrivi_ric(f,(nodo *)(radice->son)+i,l+1,z,a,riv,app);
		}
	}
	return;
}

void leggi_ric(FILE *f,nodo *radice) {
	int i,k;
	energie *en;
	char rigo[1000];
	
	for(;;) {
		if(fgets(rigo,1000,f)==NULL) return;
		for(i=0;i<1000;i++) {
			if(rigo[i]!=' ') break;
		}
		if(rigo[i]!='#') break;
	}
	if(sscanf(rigo," %d %d ",&(radice->ne),&(radice->min))<2) return;
	//printf("Leggo %d elementi, minimo %d",radice->ne,radice->min); getchar();
	if(radice->min<0) {
		en=(energie *)calloc(radice->ne,sizeof(energie));
		for(k=0;k<radice->ne;k++) {
			for(;;) {
				if(fgets(rigo,1000,f)==NULL) return;
				for(i=0;i<1000;i++) {
					if(rigo[i]!=' ') break;
				}
				if(rigo[i]!='#') break;
			}
			if(sscanf(rigo," %lg %lg %lg %lg %lg ",&(en[k].ei),&(en[k].e1),&(en[k].e2),&(en[k].e3),&(en[k].lo))<5) return;
		}
		radice->son=en;
	}
	else {
		radice->son=(void *)calloc(radice->ne,sizeof(nodo));
		for(i=0;i<radice->ne;i++) {
			leggi_ric(f,(nodo *)(radice->son)+i);
		}
	}
	return;
}

void Classe_geo::TabellaPerdite(char *fileperdite) {
	int i,k,z,a,np;
	int zmin,zmax,amin,amax,rmin,rmax;
	nodo *radice,*nodz,*noda,*nodr;
	double ei[NMAX],e1[NMAX],e2[NMAX],e3[NMAX],lo[NMAX];
	
	zmin=1;
	zmax=Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt;
	//	if(zmax>20) zmax=20;
	if(zmax>30)zmax=30;
	//printf("zmax=%d",zmax); getchar();
	rmax=8;
	FILE *f=fopen(fileperdite,"r");
	if(f==0) {
		printf("Il file %s non esiste, lo creo\n",fileperdite);
		radice=(nodo *)malloc(sizeof(nodo));
		radice->ne=zmax-zmin+1;
		radice->min=zmin;
		nodz=(nodo *)calloc(radice->ne,sizeof(nodo));
		radice->son=(void *)nodz;
		for(z=zmin;z<=zmax;z++) {
			if(z<=2) rmin=1; //NELLA CAMERA BACK SOLO FINO A Z=2
			else rmin=5; //IMF SOLO IN AVANTI
			switch(z) {
				case 1: amin=1; amax=3; break;
				case 2: amin=3; amax=4; break;
				default: amin=Classe_formule::QualeA(z); amax=amin;
			}
			if(amax>Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at) {
				amax=Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at;
				if(amin>amax) amin=amax;
			}
			nodz[z-zmin].ne=amax-amin+1;
			nodz[z-zmin].min=amin;
			noda=(nodo *)calloc(amax-amin+1,sizeof(nodo));
			nodz[z-zmin].son=(void *)noda;
			for(a=amin;a<=amax;a++) {
				//printf("        Z=%d A=%d\n",z,a);
				noda[a-amin].ne=rmax-rmin+1;
				noda[a-amin].min=rmin;
				nodr=(nodo *)calloc(rmax-rmin+1,sizeof(nodo));
				noda[a-amin].son=(void *)nodr;
				for(i=rmin;i<=rmax;i++) {
					np=calcolaGARF(z,a,i,ei,e1,e2,e3,lo);
					nodr[i-rmin].ne=np;
					nodr[i-rmin].min=-1;
					nodr[i-rmin].son=(void *)calloc(np,sizeof(energie));
					for(k=0;k<np;k++) {
						((energie *)nodr[i-rmin].son)[k].ei=ei[k];
						((energie *)nodr[i-rmin].son)[k].e2=e2[k];
						((energie *)nodr[i-rmin].son)[k].e3=e3[k];
						((energie *)nodr[i-rmin].son)[k].e1=e1[k];
						((energie *)nodr[i-rmin].son)[k].lo=lo[k];
					}
				}
				printf("Creata tabella per Z=%d A=%d. emin=%lf emax=%lf estp=%lf. N=%d\n",z,a,ei[0],ei[np-1],ei[1]-ei[0],np);
			}
		}
		f=fopen(fileperdite,"w");
		fprintf(f,"# TABELLA GARFIELD; Zmin=%d Zmax=%d\n",zmin,zmax);
		fprintf(f,"# Pressione forward  = %f mbar\n",garf.pressione);
		fprintf(f,"# Pressione backward = %f mbar\n",garf.pressione_back);
		fprintf(f,"# Target Z=%d A=%d, Materiale=%s, Spessore = %f µg/cm²\n",target.zt,target.at,target.mat.nome,target.spessore);
		fprintf(f,"############################################\n");
		fprintf(f,"#  Ein | µs down |  µs up  | Ecsi  |  L.O. #\n");
		fprintf(f,"############################################\n");
		printf("\nScrivo su disco la tabella di GARFIELD\n");
		scrivi_ric(f,radice,0,0,0,0,0);
		fclose(f);
	}
	else {
	  printf("Leggo da disco la tabella di GARFIELD...%s\n",fileperdite);
		radice=(nodo *)malloc(sizeof(nodo));
		leggi_ric(f,radice);
	}
	garf.TabGarf=radice;
	fHaveTable = 1;
	return;
}

void Classe_geo::TabellaPerditeRCO(char *fileperdite) {
	int i,k,z,a,np;
	int zmin,zmax,amin,amax,rmin,rmax;
	nodo *radice,*nodz,*noda,*nodr;
	double ei[NMAX],e1[NMAX],e2[NMAX],e3[NMAX],lo[NMAX];
	
	zmin=1;
	zmax=Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt;
	//	if(zmax>30) zmax=30;
	
	rmin=0;
	rmax=63;
	FILE *f=fopen(fileperdite,"r");
	if(f==0) {
		printf("Il file %s non esiste, lo creo\n",fileperdite);
		radice=(nodo *)malloc(sizeof(nodo));
		radice->ne=zmax-zmin+1;
		radice->min=zmin;
		nodz=(nodo *)calloc(radice->ne,sizeof(nodo));
		radice->son=(void *)nodz;
		for(z=zmin;z<=zmax;z++) {
			switch(z) {
				case 1: amin=1; amax=3; break;
				case 2: amin=3; amax=6; break;
				case 3: amin=6; amax=8; break;
				case 4: amin=7; amax=10; break;
				case 5: amin=10; amax=13; break;
				case 6: amin=11; amax=15; break;
				case 7: amin=13; amax=17; break;
				case 8: amin=15; amax=20; break;
				case 9: amin=17; amax=22; break;
				case 10: amin=19; amax=23; break;
				case 11: amin=21; amax=25; break;
				case 12: amin=23; amax=27; break;
				case 13: amin=26; amax=29; break;
				case 14: amin=28; amax=30; break;
				default:
					amin=Classe_formule::QualeA(z);
					amax=amin;
			}
			if(amax>Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at) {
				amax=Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at;
				if(amin>amax) amin=amax;
			}
			nodz[z-zmin].ne=amax-amin+1;
			nodz[z-zmin].min=amin;
			noda=(nodo *)calloc(amax-amin+1,sizeof(nodo));
			nodz[z-zmin].son=(void *)noda;
			for(a=amin;a<=amax;a++) {
				//printf("        Z=%d A=%d\n",z,a);
				noda[a-amin].ne=rmax-rmin+1;
				noda[a-amin].min=rmin;
				nodr=(nodo *)calloc(rmax-rmin+1,sizeof(nodo));
				noda[a-amin].son=(void *)nodr;
				for(i=rmin;i<=rmax;i++) {
					np=calcolaRCO(z,a,i/8,i%8,ei,e1,e2,e3,lo,&k);
					nodr[i-rmin].ne=np;
					nodr[i-rmin].min=-k;
					nodr[i-rmin].son=(void *)calloc(np,sizeof(energie));
					for(k=0;k<np;k++) {
						((energie *)nodr[i-rmin].son)[k].ei=ei[k];
						((energie *)nodr[i-rmin].son)[k].e2=e2[k];
						((energie *)nodr[i-rmin].son)[k].e3=e3[k];
						((energie *)nodr[i-rmin].son)[k].e1=e1[k];
						((energie *)nodr[i-rmin].son)[k].lo=lo[k];
					}
				}
				printf("Creata tabella per Z=%d A=%d. emin=%lf emax=%lf estp=%lf. N=%d\n",z,a,ei[0],ei[np-1],ei[1]-ei[0],np);
			}
		}
		f=fopen(fileperdite,"w");
		fprintf(f,"# TABELLA Ring Counter; Zmin=%d Zmax=%d\n",zmin,zmax);
		fprintf(f,"# Pressione IC = %f mbar\n",rco.pressione);
		fprintf(f,"# Target Z=%d A=%d, Materiale=%s, Spessore = %f µg/cm²\n",target.zt,target.at,target.mat.nome,target.spessore);
		fprintf(f,"############################################\n");
		fprintf(f,"#  Ein | µs down |  µs up  | Ecsi  |  L.O. #\n");
		fprintf(f,"############################################\n");
		printf("\nScrivo su disco la tabella del Ring Counter\n");
		scrivi_ric(f,radice,0,0,0,0,1);
		fclose(f);
	}
	else {
	  printf("Leggo da disco la tabella del Ring Counter...%s\n",fileperdite);
		radice=(nodo *)malloc(sizeof(nodo));
		leggi_ric(f,radice);
	}
	rco.TabRCO=radice;
	fHaveTableRCO = 1;
	return;
}

float Classe_geo::CalcoloEGarf(int z,int a,float estrip1,float estrip2,float ecsi,int i0) {
	float emin,emintab,estp;
	nodo *radice;
	energie *en;
	int am,ie;
	
	//if(z>2) {printf("z=%d a=%d estrip1=%f estrip2=%f ecsi=%f theta=%f, i0=%d",z,a,estrip1,estrip2,ecsi,theta,i0); getchar();}
	radice=garf.TabGarf;
	if(z<2) {estrip1=0.; estrip2=0.;}
	if(estrip1<0.) estrip1=0.;
	if(estrip2<0.) estrip2=0.;
	emin=ecsi+estrip1+estrip2;
	
	if((z<radice->min)||(z>=radice->ne+radice->min)) return -1.; //z out of range
	radice=(nodo *)radice->son+(z-radice->min);
	
	if(a<=0) a=Classe_formule::QualeA(z);
	if(a<=0) {am=radice->min+radice->ne/2;}
	else {
		if(a<radice->min) {am=radice->min;}
		else {
			if(a>=radice->ne+radice->min) {am=radice->min+radice->ne-1;}
			else am=a;
		}
	}
	radice=(nodo *)radice->son+(am-radice->min);
	
	if((i0+1<radice->min)||(i0+1>=radice->ne+radice->min)) return -1.; //i0 out of range, in pratica se si cerca di vedere uno Z>2 in back.
	radice=(nodo *)(radice->son)+(i0+1-radice->min);
	en=(energie *)(radice->son);
	emintab=en[0].ei;
	if(emin>en[radice->ne-1].ei) {
		printf("[TABELLA CORTA!] z=%d a=%d GARF csi %d estrip+ecsi=%lf\n",z,a,i0+1,emin);
		//return emin; //tabella troppo corta, esco con la somma delle energie gas+csi
		return -1; //  28/6/2016 se succede qualcosa di strano meglio uscire con -1
	}
	estp=en[1].ei-en[0].ei; //28/4/2015 definizione universale di estp!
	ie=(emin-en[0].ei)/estp-1;
	if(ie<0) ie=0;
	
	//15/4/2015 uso ecsi+estrip come Estart
	if(ie>radice->ne-3) ie=radice->ne-3;
	for(;(en[ie+2].e1+en[ie+2].e2+en[ie+2].e3<emin)&&(ie<radice->ne-3);ie++);
	
	float dminecsi;
	float dminestrip1;
	float dminestrip2;
	float riscsi=GCSIRES*ecsi+GCSIMINRES;
	float risstrip1=GGASRES*estrip1+GGASMINRES;
	float risstrip2=GGASRES*estrip2+GGASMINRES;
	float d1=-1.,d2=-1.,d3=-1.;//,E1,E2,E3;
	TRandom3 rando; //29/6/16 via l'interpolazione: si spalma!
	rando.SetSeed(0);
	
	dminecsi=pow(TMath::Abs(ecsi-en[ie].e3)/riscsi,2);
	if(estrip1>0.001) dminestrip1=pow(TMath::Abs(estrip1-en[ie].e1)/risstrip1,2); else dminestrip1=0.;
	if(estrip2>0.001) dminestrip2=pow(TMath::Abs(estrip2-en[ie].e2)/risstrip2,2); else dminestrip2=0.;
	
	d3=TMath::Sqrt(dminecsi+dminestrip1+dminestrip2);
	for(ie++;ie<radice->ne;ie++) {
		d1=d2; d2=d3;
		dminecsi=pow(TMath::Abs(ecsi-en[ie].e3)/riscsi,2);
		if(estrip1>0.001) dminestrip1=pow(TMath::Abs(estrip1-en[ie].e1)/risstrip1,2); else dminestrip1=0.;
		if(estrip2>0.001) dminestrip2=pow(TMath::Abs(estrip2-en[ie].e2)/risstrip2,2); else dminestrip2=0.;
		d3=TMath::Sqrt(dminecsi+dminestrip1+dminestrip2);
		if(d3>d2) break;
	}
	if(d3<=d2) ie--;
	
	if((d1>=0)&&(d2>=0)) {
		
		//   28/06/2016  RIMOSSA PARABOLA
		return rando.Uniform(en[ie-1].ei-estp/2.,en[ie-1].ei+estp/2.);
		
// 		E1=en[ie-2].ei;
// 		E2=en[ie-1].ei;
// 		E3=en[ie].ei;
// 		return (E1*E1*(d3-d2)+E2*E2*(d1-d3)+E3*E3*(d2-d1))/(2.*(d1*(E2-E3)+d2*(E3-E1)+d3*(E1-E2)));
	}
	if(d2<0) { //Questa condizione non puo' e non deve verificarsi mai!!
		printf("DEBUG [5] z=%d a=%d GARF csi %d estrip+ecsi=%lf\n",z,a,i0+1,emin);
		//return emin;
		return -1; //  28/6/2016 se succede qualcosa di strano meglio uscire con -1
	}
	//Sono partito troppo avanti, torno indietro!
	for(;ie>=2;ie--) {
		if(d1>=0) {d3=d2; d2=d1;}
		dminecsi=pow(TMath::Abs(ecsi-en[ie-2].e3)/riscsi,2);
		if(estrip1>0.001) dminestrip1=pow(TMath::Abs(estrip1-en[ie-2].e1)/risstrip1,2); else dminestrip1=0.;
		if(estrip2>0.001) dminestrip2=pow(TMath::Abs(estrip2-en[ie-2].e2)/risstrip2,2); else dminestrip2=0.;
		d1=TMath::Sqrt(dminecsi+dminestrip1+dminestrip2);
		if(d1>d2) break;
	}
	if((d1<=d2)&&(d1>=0)) ie++;
	if(d1>=0) {
		
		//   28/06/2016  RIMOSSA PARABOLA
		return rando.Uniform(en[ie-1].ei-estp/2.,en[ie-1].ei+estp/2.);
		
// 		E1=en[ie-2].ei;
// 		E2=en[ie-1].ei;
// 		E3=en[ie].ei;
// 		return (E1*E1*(d3-d2)+E2*E2*(d1-d3)+E3*E3*(d2-d1))/(2.*(d1*(E2-E3)+d2*(E3-E1)+d3*(E1-E2)));
	}
	//RESTANO I CASI IN CUI LA PRIMA RIGA DELLA TABELLA e' Gia' TROPPO GRANDE
	//printf("DEBUG [1] z=%d a=%d GARF csi %d estrip+ecsi=%lf\n",z,a,i0+1,emin);
	return rando.Uniform(emintab-estp/2.,emintab+estp/2.);
}


float Classe_geo::CalcoloERCO(int z,int a,float egas,float esi,float ecsi,int bit,int codice) {
	float emin,emintab,estp;
	nodo *radice,*radamax;
	energie *en,*enmax;
	int am,ie,iemin,iestart,iesi,iecsi,istop;
	int isec=codice/100;
	int istrip=(codice%100)/10;
	int icsi=codice%10;
	
	//if(z==19&&a==40) printf("CalcoloERCO:\n");
	//if(z==19&&a==40) printf(" inizio: z=%3d a=%3d egas=%f esi=%f ecsi=%f bit=%d codice=%d\n",z,a,egas,esi,ecsi,bit,codice);
	if(egas<=0) {bit&=6;}
	if(esi<=0) {bit&=5;}
	if(ecsi<=0) {bit&=3;}
	if(bit==0) return -1.;
	if((bit&1)==0) egas=0;
	if((bit&2)==0) esi=0;
	if((bit&4)==0) ecsi=0;
	emin=egas+esi+ecsi;
	
	radice=rco.TabRCO;
	if((z<radice->min)||(z>=radice->ne+radice->min)) return -1.; //z out of range
	radice=(nodo *)radice->son+(z-radice->min);
	
	if(a<=0) a=Classe_formule::QualeA(z);
	if(a<=0) {am=radice->min+radice->ne/2;}
	else {
		if(a<radice->min) {am=radice->min;}
		else {
			if(a>=radice->ne+radice->min) {am=radice->min+radice->ne-1;}
			else am=a;
		}
	}
	radamax=(nodo *)radice->son+(radice->ne-1);
	radice=(nodo *)radice->son+(am-radice->min);
	//VECCHIA DEFINIZIONE DI ESTP!!
	//estp=am*ESTPN;
	//if(estp>ESTP) estp=ESTP;
	
	if((icsi<=0)||(icsi>6)) {
		if((istrip<=0)||(istrip>8)) istop=1;
		else istop=2;
	}
	else istop=3;
	if((isec<=0)||(isec>8)) isec=1; //settore non specificato, prendo l'1
	if((istrip<=0)||(istrip>=8)) {
		//istrip out of range, in pratica la strip non funziona.
		//prendo una strip ragionevole!
		if((icsi<=0)||(icsi>6)) istrip=4;
		else {
			if((icsi==1)||(icsi==2)) istrip=1;
			if((icsi==3)||(icsi==4)) istrip=4;
			if((icsi==5)||(icsi==6)) istrip=7;
		}
	}
	int iriv=8*isec+istrip-9; //sarebbe 8*(isec-1)+(istrip-1)
	if((iriv<radice->min)||(iriv>=radice->ne+radice->min)) return -1.; //iriv out of range, non dovrebbe mai accadere!
	radice=(nodo *)(radice->son)+(iriv-radice->min);
	if((iriv<radamax->min)||(iriv>=radamax->ne+radamax->min)) radamax=radice;
	else radamax=(nodo *)(radamax->son)+(iriv-radamax->min);
	en=(energie *)(radice->son);
	enmax=(energie *)(radamax->son);
	emintab=en[0].ei;
	if(emin>en[radice->ne-1].ei) {
		printf("[TABELLA CORTA!] z=%d a=%d RCO codice %d egas+esi+ecsi=%lf\n",z,a,codice,emin);
		//return emin; //tabella troppo corta, esco con la somma delle energie gas+si+csi
		return -1; //  28/6/2016 se succede qualcosa di strano meglio uscire con -1
	}
	//CERCO UN'ENERGIA DI PARTENZA RAGIONEVOLE
	estp=en[1].ei-en[0].ei; //28/4/2015 definizione universale di estp!
	ie=(emin-emintab)/estp;
	if(ie<0) ie=0;
	if(ie>radice->ne-3) ie=radice->ne-3;
	iemin=0; iestart=0;
	iesi=-radice->min%10000-1;
	iecsi=-radice->min/10000-1;
	if(iesi>0) {
		if((en[iesi-1].e2>0)||(en[iesi].e2<=0)) printf("DEBUG [2]: z=%d a=%d codice=%d iesi=%d\n",z,a,codice,iesi);
	}
	if(iecsi>0) {
		if((en[iecsi-1].e3>0)||(en[iecsi].e3<=0)) printf("DEBUG [3]: z=%d a=%d codice=%d iecsi=%d\n",z,a,codice,iecsi);
		if((esi>en[iecsi-1].e2)&&(esi>en[iecsi].e2)) {
			//Esi misurato e' maggiore della massima energia possibile!
			//sperimentalmente si vede che questa condizione si verifica sempre e solo quando la camera
			//soffre e raccoglie meno carica; quindi viene identificato uno Z sbagliato e chiaramente per quello
			//Z l'energia della cuspide e' diversa. L'unico modo per pulire i dati e' uscire in questo caso!
			//printf("DEBUG [4]: z=%d a=%d codice=%d esi=%f\n",z,a,codice,esi);
			//5/10/2015 NON SI SPALMA, SI ESCE!
			//TRandom3 rando;
			//rando.SetSeed(0);
			//return rando.Uniform(en[iecsi-1].ei,en[iecsi].ei);
			iecsi=-radamax->min/10000-1;
			if(iecsi<=0) return -1;
			if((esi>enmax[iecsi-1].e2)&&(esi>enmax[iecsi].e2)) return -1;
			//printf("DEBUG [4]: z=%d a=%d codice=%d esi=%f RECUPERATO\n",z,a,codice,esi); getchar();
			radice=radamax; en=enmax;
			emintab=en[0].ei;
			if(emin>en[radice->ne-1].ei) {
				printf("[TABELLA CORTA!] z=%d a=%d RCO codice %d egas+esi+ecsi=%lf\n",z,a,codice,emin);
				//return emin; //tabella troppo corta, esco con la somma delle energie gas+si+csi
				return -1; //  28/6/2016 se succede qualcosa di strano meglio uscire con -1
			}
			estp=en[1].ei-en[0].ei;
			ie=(emin-emintab)/estp;
			if(ie<0) ie=0;
			if(ie>radice->ne-3) ie=radice->ne-3;
			iemin=0; iestart=0;
			iesi=-radice->min%10000-1;
		}
	}
	if((istop==1)&&(ie>=iesi)&&(iesi>0)) ie=iesi-1;
	if(istop==2) {
		if((ie>=iecsi)&&(iecsi>0)) ie=iecsi-1;
		for(;ie>iesi;ie--) {
			if((en[ie+1].e2<esi)||esi==0.) break;
		}
		if(iesi>0) {
			iemin=iesi-1;
			for(iestart=iemin;iestart<radice->ne-2;iestart++) {
				if(en[iestart].e1>en[iestart+1].e1) break;
			}
			//CUSPIDE DELLA PERDITA IN CAMERA
		}
	}
	if((istop==3)&&(iecsi>0)) {
		for(;ie>iecsi;ie--) {
			if((en[ie+1].e2>esi)&&(en[ie+1].e3<ecsi||ecsi==0)) break;
		}
		if(iecsi>0) {
			iemin=iecsi-1;
			for(iestart=iemin;iestart<radice->ne-2;iestart++) {
				if(en[iestart].e2>en[iestart+1].e2) break;
			}
			//CUSPIDE DELLA PERDITA IN Si
		}
	}
	if(iestart>ie) ie=iestart;
	
	//if(z==19&&a==40) printf("                     egas=%f esi=%f ecsi=%f => emin=%f\n",egas,esi,ecsi,emin);
	//if(z==19&&a==40) printf("  RIGA %5d         egas=%f esi=%f ecsi=%f => etot=%f\n",ie,en[ie].e1,en[ie].e2,en[ie].e3,en[ie].ei);
	
	float dminecsi;
	float dminesi;
	float dminegas;
	float riscsi=RCSIRES*ecsi+RCSIMINRES;
	float risgas=RGASRES*egas+RGASMINRES;
	float rissi=RSIRES*esi+RSIMINRES;
	float d1=-1.,d2=-1.,d3=-1.;//,E1,E2,E3;
	TRandom3 rando; //29/6/16 via l'interpolazione: si spalma!
	rando.SetSeed(0);
	
	if(egas>0.) dminegas=pow(TMath::Abs(egas-en[ie].e1)/risgas,2); else dminegas=0.;
	if(esi>0.) dminesi=pow(TMath::Abs(esi-en[ie].e2)/rissi,2); else dminesi=0.;
	if(ecsi>0.) dminecsi=pow(TMath::Abs(ecsi-en[ie].e3)/riscsi,2); else dminecsi=0.;
	
	d3=TMath::Sqrt(dminecsi+dminesi+dminegas);
	for(ie++;ie<radice->ne;ie++) {
		d1=d2; d2=d3;
		if(egas>0.) dminegas=pow(TMath::Abs(egas-en[ie].e1)/risgas,2); else dminegas=0.;
		if(esi>0.) dminesi=pow(TMath::Abs(esi-en[ie].e2)/rissi,2); else dminesi=0.;
		if(ecsi>0.) dminecsi=pow(TMath::Abs(ecsi-en[ie].e3)/riscsi,2); else dminecsi=0.;
		d3=TMath::Sqrt(dminecsi+dminesi+dminegas);
		if(d3>d2) break;
	}
	if(d3<=d2) ie--;
	
	//if(z==19&&a==40) printf("  RIGA %5d (fw)    egas=%f esi=%f ecsi=%f => etot=%f\n",ie-1,en[ie-1].e1,en[ie-1].e2,en[ie-1].e3,en[ie-1].ei);
	
	if((d1>=0)&&(d2>=0)) {
		
		//   28/06/2016  RIMOSSA PARABOLA
		return rando.Uniform(en[ie-1].ei-estp/2.,en[ie-1].ei+estp/2.);
		
// 		E1=en[ie-2].ei;
// 		E2=en[ie-1].ei;
// 		E3=en[ie].ei;
// 		//if(z==8&&((bit&4)==0)) printf("(1) % 5.1f %8.3f  (2) % 5.1f %8.3f  (3) % 5.1f %8.3f\n",E1,d1,E2,d2,E3,d3);
// 		return (E1*E1*(d3-d2)+E2*E2*(d1-d3)+E3*E3*(d2-d1))/(2.*(d1*(E2-E3)+d2*(E3-E1)+d3*(E1-E2)));
	}
	if(d2<0) { //Questa condizione non puo' e non deve verificarsi mai!!
		printf("DEBUG [5] z=%d a=%d RCO codice %d egas+esi+ecsi=%lf\n",z,a,codice,emin);
		//return emin;
		return -1; //  28/6/2016 se succede qualcosa di strano meglio uscire con -1
	}
	//Sono partito troppo avanti, torno indietro!
	for(;ie-2>=iemin;ie--) {
		if(d1>=0) {d3=d2; d2=d1;}
		if(egas>0.) dminegas=pow(TMath::Abs(egas-en[ie-2].e1)/risgas,2); else dminegas=0.;
		if(esi>0.) dminesi=pow(TMath::Abs(esi-en[ie-2].e2)/rissi,2); else dminesi=0.;
		if(ecsi>0.) dminecsi=pow(TMath::Abs(ecsi-en[ie-2].e3)/riscsi,2); else dminecsi=0.;
		d1=TMath::Sqrt(dminecsi+dminesi+dminegas);
		if(d1>d2) break;
	}
	if((d1<=d2)&&(d1>=0)) ie++;
	
	//if(z==19&&a==40) printf("  RIGA %5d (bw)    egas=%f esi=%f ecsi=%f => etot=%f\n",ie-1,en[ie-1].e1,en[ie-1].e2,en[ie-1].e3,en[ie-1].ei);
	
	if(d1>=0) {
		
		//   28/06/2016  RIMOSSA PARABOLA
		return rando.Uniform(en[ie-1].ei-estp/2.,en[ie-1].ei+estp/2.);
		
// 		E1=en[ie-2].ei;
// 		E2=en[ie-1].ei;
// 		E3=en[ie].ei;
// 		return (E1*E1*(d3-d2)+E2*E2*(d1-d3)+E3*E3*(d2-d1))/(2.*(d1*(E2-E3)+d2*(E3-E1)+d3*(E1-E2)));
	}
	
	//RESTANO I CASI IN CUI LA PRIMA RIGA UTILE DELLA TABELLA e' GIa' TROPPO GRANDE
	//printf("DEBUG [1] z=%d a=%d RCO codice %d egas+esi+ecsi=%lf\n",z,a,codice,emin);
	if(iemin==0) return rando.Uniform(emintab-estp/2.,emintab+estp/2.);
	else return rando.Uniform(en[iemin].ei,en[iemin+1].ei);
}





void Classe_geo::SpalmaRCO(int isec,int istrip,int icsi,float *theout,float *phiout)
{

  InizializzaRCO();
 
  if(istrip<8 && icsi<6)//strip definita, cesio definito=>phi da cesio, theta da strip
    {
      int dentro=-1;
      while(dentro<0)
	{
	  if(icsi%2==0)//icsi pari (0,2,4)
	    {
	  *phiout=rco.gas_phimin[isec]+rco.dphi/2*gRandom->Rndm(0);
	    }
	  else//icsi dispari (1,3,5)
	    {
	  *phiout=rco.gas_phimax[isec]-rco.dphi/2*gRandom->Rndm(0);

	    }

*theout=rco.strip_themin[isec][istrip]+(rco.strip_themax[isec][istrip]-rco.strip_themin[isec][istrip])*gRandom->Rndm(0);
 dentro=0;
 if(istrip==7)
   {
    float ts=(pow(rco.strip_veccentro[0],2)+pow(rco.strip_veccentro[1],2)+pow(rco.strip_veccentro[2],2))/(rco.strip_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.strip_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.strip_veccentro[2]*cos(*theout/57.296));
    float pis,pjs;
	      float puntos[3];
	      puntos[0]=ts*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      puntos[1]=ts*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      puntos[2]=ts*cos(*theout/57.296);
	      float punto_cs[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_cs[nn]=puntos[nn]-rco.strip_veccentro[nn];
		}
	     
	      Classe_formule::scaprod(punto_cs,rco.gas_vi[isec],&pis);
	      Classe_formule::scaprod(punto_cs,rco.gas_vj[isec],&pjs);
	      if(pjs<rco.strip_yg)
		{
		  dentro=-1;
		}
   }//strip 7

 if(dentro>=0)
   {
float t=(pow(rco.gas_veccentro[0],2)+pow(rco.gas_veccentro[1],2)+pow(rco.gas_veccentro[2],2))/(rco.gas_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[2]*cos(*theout/57.296));
	      float punto[3];
	      punto[0]=t*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      punto[1]=t*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      punto[2]=t*cos(*theout/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-rco.gas_veccentro[nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_c,rco.gas_vi[isec],&pi);
	      Classe_formule::scaprod(punto_c,rco.gas_vj[isec],&pj);
 

	      float yrif=rco.gas_ya/rco.gas_xa*pi+rco.dx_spe*sqrt(1+pow(rco.gas_ya/rco.gas_xa,2));
	      float xb=-rco.gas_xa;
	      float yrifb=rco.gas_ya/xb*pi+rco.sx_spe*sqrt(1+pow(rco.gas_ya/xb,2));

	      if(pj<yrif)
		{
		  dentro=-1;
		}
	      if(pj<yrifb)
		{
		  dentro=-1;
		}

	      if(pj<rco.gas_yg)//la strip 8 in basso e' dritta (e non tonda)
		{
		  dentro=-1;
		}
	      if(dentro>=0)
		{
	      float tcsi=(pow(rco.csi_veccentro[0],2)+pow(rco.csi_veccentro[1],2)+pow(rco.csi_veccentro[2],2))/(rco.csi_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.csi_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.csi_veccentro[2]*cos(*theout/57.296));
	      float puntocsi[3];
	      puntocsi[0]=tcsi*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      puntocsi[1]=tcsi*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      puntocsi[2]=tcsi*cos(*theout/57.296);
	      float punto_ccsi[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_ccsi[nn]=puntocsi[nn]-rco.csi_veccentro[nn];
		}
	      float picsi,pjcsi;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vi[isec],&picsi);
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vj[isec],&pjcsi);
	    
	    
	      if(rco.csi[icsi]->IsInside(picsi,pjcsi)==0)
		{
		  dentro=-1;
		}
   }
   }
	}//dentro<0
    }
  if(istrip==8 && icsi<6)//strip non definita, cesio definito => theta e phi da cesio
    {
      int dentro=-1;
      while(dentro<0)
	{
	  if(icsi%2==0)
	    {
	  *phiout=rco.gas_phimin[isec]+rco.dphi/2*gRandom->Rndm(0);
	    }
	  else
	    {
	  *phiout=rco.gas_phimax[isec]-rco.dphi/2*gRandom->Rndm(0);
	    }

  if(rco.stripcopertamax>=0)//fino alla strip coperta non si vede niente (c'e' uno schermo di alluminio)
    {
*theout=rco.strip_themax[isec][rco.stripcopertamax]+(rco.gas_themax[isec]-rco.strip_themax[isec][rco.stripcopertamax])*gRandom->Rndm(0);
    }
  else
    {
*theout=rco.gas_themin[isec]+(rco.gas_themax[isec]-rco.gas_themin[isec])*gRandom->Rndm(0);
    }
 dentro=0;
	      float tcsi=(pow(rco.csi_veccentro[0],2)+pow(rco.csi_veccentro[1],2)+pow(rco.csi_veccentro[2],2))/(rco.csi_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.csi_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.csi_veccentro[2]*cos(*theout/57.296));
	      float puntocsi[3];
	      puntocsi[0]=tcsi*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      puntocsi[1]=tcsi*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      puntocsi[2]=tcsi*cos(*theout/57.296);
	      float punto_ccsi[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_ccsi[nn]=puntocsi[nn]-rco.csi_veccentro[nn];
		}
	      float picsi,pjcsi;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vi[isec],&picsi);
	      Classe_formule::scaprod(punto_ccsi,rco.gas_vj[isec],&pjcsi);
	      if(rco.csi[icsi]->IsInside(picsi,pjcsi)==0)
		{
		  dentro=-1;
		}
	      if(dentro>=0)
		{
float t=(pow(rco.gas_veccentro[0],2)+pow(rco.gas_veccentro[1],2)+pow(rco.gas_veccentro[2],2))/(rco.gas_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[2]*cos(*theout/57.296));
	      float punto[3];
	      punto[0]=t*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      punto[1]=t*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      punto[2]=t*cos(*theout/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-rco.gas_veccentro[nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_c,rco.gas_vi[isec],&pi);
	      Classe_formule::scaprod(punto_c,rco.gas_vj[isec],&pj);
 

	      float yrif=rco.gas_ya/rco.gas_xa*pi+rco.dx_spe*sqrt(1+pow(rco.gas_ya/rco.gas_xa,2));
	      float xb=-rco.gas_xa;
	      float yrifb=rco.gas_ya/xb*pi+rco.sx_spe*sqrt(1+pow(rco.gas_ya/xb,2));
	      
	      if(pj<yrif)
		{
		  dentro=-1;
		}
	      if(pj<yrifb)
		{
		  dentro=-1;
		}

	      if(pj<rco.gas_yg)//la strip 8 in basso e' dritta (e non tonda)
		{
		  dentro=-1;
		}
		}
	}//dentro<0
    }
  if(istrip<8 && icsi==6)//strip definita, cesio no => phi da gas e theta dalle strip
    {
      int dentro=-1;
      while(dentro<0)
	{
      if(rco.gas_phimin[isec]*rco.gas_phimax[isec]>=0)//min e max tutti e 2 positivi o tutti e 2 negativi
	{
	  *phiout=rco.gas_phimin[isec]+rco.dphi*gRandom->Rndm(0);

	}
      else//min positivo max negativo o viceversa
	{
	  if(gRandom->Rndm(0)<0.5)
	    {
	      *phiout=rco.gas_phimin[isec]+rco.dphi/2*gRandom->Rndm(0);

	    }
	  else
	    {
	      *phiout=rco.gas_phimax[isec]-rco.dphi/2*gRandom->Rndm(0);

	    }
	}

*theout=rco.strip_themin[isec][istrip]+(rco.strip_themax[isec][istrip]-rco.strip_themin[isec][istrip])*gRandom->Rndm(0);

float t=(pow(rco.gas_veccentro[0],2)+pow(rco.gas_veccentro[1],2)+pow(rco.gas_veccentro[2],2))/(rco.gas_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[2]*cos(*theout/57.296));
	      float punto[3];
	      punto[0]=t*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      punto[1]=t*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      punto[2]=t*cos(*theout/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-rco.gas_veccentro[nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_c,rco.gas_vi[isec],&pi);
	      Classe_formule::scaprod(punto_c,rco.gas_vj[isec],&pj);
 

	      float yrif=rco.gas_ya/rco.gas_xa*pi+rco.dx_spe*sqrt(1+pow(rco.gas_ya/rco.gas_xa,2));
	      float xb=-rco.gas_xa;
	      float yrifb=rco.gas_ya/xb*pi+rco.sx_spe*sqrt(1+pow(rco.gas_ya/xb,2));
	      dentro=0;
	      if(pj<yrif)
		{
		  dentro=-1;
		}
	      if(pj<yrifb)
		{
		  dentro=-1;
		}

	      if(pj<rco.gas_yg)//la strip 8 in basso e' dritta (e non tonda)
		{
		  dentro=-1;
		}
	      if(istrip==7 && dentro==0)
		{

     float ts=(pow(rco.strip_veccentro[0],2)+pow(rco.strip_veccentro[1],2)+pow(rco.strip_veccentro[2],2))/(rco.strip_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.strip_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.strip_veccentro[2]*cos(*theout/57.296));

	      float puntos[3];
	      puntos[0]=ts*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      puntos[1]=ts*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      puntos[2]=ts*cos(*theout/57.296);
	      float punto_cs[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_cs[nn]=puntos[nn]-rco.strip_veccentro[nn];
		}
	      float pis,pjs;
	      Classe_formule::scaprod(punto_cs,rco.gas_vi[isec],&pis);
	      Classe_formule::scaprod(punto_cs,rco.gas_vj[isec],&pjs);
	      if(pjs<rco.strip_yg)
		{
		  dentro=-1;
		}
		}//strip 7 (da 0)

	}//while dentro<0


    }
  if(istrip==8 && icsi==6)//ne' strip ne' cesio definiti, si prende theta e phi dal gas
    {
      int dentro=-1;
      while(dentro<0)
	{

      if(rco.gas_phimin[isec]*rco.gas_phimax[isec]>=0)//min e max tutti e 2 positivi o tutti e 2 negativi
	{
	  *phiout=rco.gas_phimin[isec]+rco.dphi*gRandom->Rndm(0);

	}
      else//min positivo max negativo o viceversa
	{
	  if(gRandom->Rndm(0)<0.5)
	    {
	      *phiout=rco.gas_phimin[isec]+rco.dphi/2*gRandom->Rndm(0);

	    }
	  else
	    {
	      *phiout=rco.gas_phimax[isec]-rco.dphi/2*gRandom->Rndm(0);

	    }
	}
  if(rco.stripcopertamax>=0)//fino alla strip coperta non si vede niente (c'e' uno schermo di alluminio)
    {
*theout=rco.strip_themax[isec][rco.stripcopertamax]+(rco.gas_themax[isec]-rco.strip_themax[isec][rco.stripcopertamax])*gRandom->Rndm(0);
    }
  else
    {
*theout=rco.gas_themin[isec]+(rco.gas_themax[isec]-rco.gas_themin[isec])*gRandom->Rndm(0);
    }
float t=(pow(rco.gas_veccentro[0],2)+pow(rco.gas_veccentro[1],2)+pow(rco.gas_veccentro[2],2))/(rco.gas_veccentro[0]*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[1]*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296)+rco.gas_veccentro[2]*cos(*theout/57.296));
	      float punto[3];
	      punto[0]=t*sin(*theout/57.296)*cos(TMath::Pi()/2-*phiout/57.296);
	      punto[1]=t*sin(*theout/57.296)*sin(TMath::Pi()/2-*phiout/57.296);
	      punto[2]=t*cos(*theout/57.296);
	      float punto_c[3];
	      for(int nn=0;nn<3;nn++)
		{
	      punto_c[nn]=punto[nn]-rco.gas_veccentro[nn];
		}
	      float pi,pj;//coordinate sul piano del riv in un SDR che ha l'asse y come il centro dello spicchio
	      Classe_formule::scaprod(punto_c,rco.gas_vi[isec],&pi);
	      Classe_formule::scaprod(punto_c,rco.gas_vj[isec],&pj);
 

	      float yrif=rco.gas_ya/rco.gas_xa*pi+rco.dx_spe*sqrt(1+pow(rco.gas_ya/rco.gas_xa,2));
	      float xb=-rco.gas_xa;
	      float yrifb=rco.gas_ya/xb*pi+rco.sx_spe*sqrt(1+pow(rco.gas_ya/xb,2));
	      dentro=0;
	      if(pj<yrif)
		{
		  dentro=-1;
		}
	      if(pj<yrifb)
		{
		  dentro=-1;
		}

	      if(pj<rco.gas_yg)//la strip 8 in basso e' dritta (e non tonda)
		{
		  dentro=-1;
		}


	}
    }


  return;
}

float Classe_geo::Range(float zpart,float apart,float epart,int mate)
{
  float spess=1;//micron
  int idir=1;
  int icod;
  float dummy_pressione=0;
  float atloc=0;
  float elost,epost;
  while(1)
    {
  ecorr_veda(&epart,&zpart,&apart,&atloc,&epost,&elost,&mate,&spess,&idir,&icod,&dummy_pressione);    
  if(epost<=0)
    {
      return spess-1;
    }
      spess++;
      if(spess>1000000)
	{
	  cout<<"Piu' di 1000000 iterazioni in range "<<zpart<<" "<<apart<<" "<<epart<<" "<<mate<<" "<<spess<<endl;
	}
    }
  return 0;
}
 float Classe_geo::DE2E(float zpart,float apart,float elost,int mate,float spess)
 {
 //int idir=1;
  //int icod;
  float dummy_pressione=0;
  float atloc=0;
  float e=-1;
  #ifndef _ODIE_
      de_vedaloss_(&zpart,&apart,&atloc,&elost,&spess,&e,&mate,&dummy_pressione);
  #endif
      return e;
 }

float Classe_geo::Epunch(float zpart,float apart,float spess,int mate)
{
  float epart=0.1;
  int idir=1;
  int icod;
  float dummy_pressione=0;
  float atloc=0;
  float elost,epost;
  while(1)
    {
  ecorr_veda(&epart,&zpart,&apart,&atloc,&epost,&elost,&mate,&spess,&idir,&icod,&dummy_pressione);    
  if(epost>0)
    {
      return epart-0.1;
    }
      epart=epart+0.1;
      if(epart>10000)
	{
	  cout<<"Piu' di 100000 iterazioni in Epunch "<<zpart<<" "<<apart<<" "<<epart<<" "<<mate<<" "<<spess<<endl;
	}
    }
  return 0;
}
