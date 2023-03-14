#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "perdenc.h"

#define EMIN 0.01
#define EMAX 1000
#define PREC 0.0001
#define ITER 1000

typedef struct Tabededx {
	double e,dd;
	struct Tabededx *next;
} tabededx;

typedef struct tabDB {
	char nome[32];
	int z;
	tabededx *tsrim;
	tabededx *tnist;
	tabededx *tbarb;
	tabededx *tschw;
	struct tabDB *next;
} tabdb;

extern void ecorr_veda_(float *eingresso,float *zpr,float *apr,float *atar, float *eout,float *elost,int *mate,float *thick,int *idir, int *icod, float *pressione);
extern void de_vedaloss_(float *zpr,float *apr,float *atar,float *de,float *thick,float *e,int *mate,float *pressione);

tabdb *db=NULL;

const double M[86]={1.008,4.003,6.941,9.012,10.81,12.01,14.01,16.00,19.00,20.18, //H->Ne
	22.99,24.31,26.98,28.09,30.97,32.07,35.45,39.95, //Na->Ar
	39.10,40.08,44.96,47.87,50.94,52.00,54.94,55.85,58.93,58.69,63.55,65.41,69.72,72.64,74.92,78.96,79.90,83.80, //K->Kr
	85.47,87.62,88.91,91.22,92.91,95.94,98.00,101.1,102.9,106.4,107.9,112.4,114.8,118.7,121.8,127.6,126.9,131.3, //Rb->Xe
	132.9,137.3,138.9,140.1,140.9,144.2,144.9,150.4,152.0,157.3,158.9,162.5,164.9,167.3,168.9,173.1,175.0,178.5,181.0,183.8,186.2,190.2,192.2,195.1,197.0,200.6,204.4,207.2,209.0,209.0,210.0,222.0//Cs->Rn
};

double schwalm(double e,double Z1,double Z2,double A2,double rho) {
	//double dedxnu,dedxel;
	//double z,eps;
	double v,vq,y,xi,h,r,AL,g;
	double H[101]={8.080,5.706,0.609,-6.214,-10.646,-9.119,-3.171,2.205,3.656,1.097,-3.536,-7.364,-7.897,-4.762,0.120,4.370,6.883,7.880,8.016,7.764,
		7.327,6.766,6.100,5.347,4.534,3.701,2.898,2.182,1.606,1.219,1.054,1.125,1.426,1.930,2.595,3.368,4.193,5.016,5.793,6.490,
		7.087,7.575,7.957,8.234,8.384,8.239,7.350,5.386,3.265,2.891,4.604,6.653,7.722,7.937,7.806,7.575,7.303,6.998,6.661,6.294,
		5.901,5.487,5.056,4.617,4.177,3.746,3.334,2.950,2.604,2.305,2.063,1.883,1.771,1.730,1.762,1.866,2.038,2.273,2.565,2.906,
		3.286,3.696,4.125,4.564,5.004,5.436,5.853,6.249,6.618,6.959,7.269,7.546,7.791,8.004,8.188,8.345,8.477,8.586,8.676,8.749,8.807};
		
		if(e<EMIN) e=EMIN;
		
		//electronic dedx
		vq=e/469.;
		v=sqrt(vq);
		xi=vq/Z2;
		h=(Z2>100)?9.:H[(int)Z2];
		//h=0;
		r=v/pow(Z1,0.509);
		AL=log(xi*(4.444*Z2*Z2+922.2*A2/rho));
		y=3.3e-4*log(1.+54721.*xi*(1.76e-2*(A2/rho)*(1.-exp(-pow(0.054*Z2,5)))+1.-exp(-0.3*Z2)))+(2.*sqrt(xi)/(Z2*(1.+1.e4*sqrt(xi)))-1.32e-5*h*AL*exp(-0.32*AL*AL))/(1.+1.e12*xi*xi*xi);
		g=(Z1<=2)?(1.-exp(-116.79*r-3350.4*r*r)):(1.-pow(1.035-0.4*exp(-0.16*Z1),1.-exp(-137.*v))*exp(-120.4*v/pow(Z1,0.65)));
		
		return g*g*Z1*Z1*Z2*y/(A2*vq);
}

double dedx(int ne,int zp,int *zt,int *at,double *w,double rho,double e) {
	double dEdx=0.,wtot=0.;
	int j;
	
	for(j=0;j<ne;j++) {
		dEdx+=schwalm(e,(double)zp,(double)(zt[j]),(double)(at[j]),rho/1000.)*w[j];
		wtot+=w[j];
	}
	return dEdx/wtot;
}

double nudedx(int ne,int zp,int ap,int *zt,int *at,double *w,double e) {
	double dEdx=0.,wtot=0.,z,eps,dedxnu;
	int j;
	
	for(j=0;j<ne;j++) {
		z=(ap+at[j])*sqrt(pow(zp,0.6666667)+pow(zt[j],0.6666667));
		eps=3.25e4*e*ap*at[j]/(z*zp*zt[j]);
		dedxnu=8.6785*(sqrt(eps)*log(eps+2.1718)/(1.+6.8*eps+3.4*pow(eps,1.5)))*zp*zt[j]*ap/(z*at[j]);
		dEdx+=dedxnu*w[j];
		wtot+=w[j];
	}
	return dEdx/wtot;
}

double leggimate(MATE *mat,int *zt,int *at,double *w) {
	int i;
	double Mtot=0,r=-1;
	
	double rho[86]={-1,-1,535,1848,2460,2260,-1,-1,-1,-1, //H->Ne
		968,1738,2700,2330,1823,1960,-1,-1, //Na->Ar
		856,1550,2985,4507,6110,7140,7470,7874,8900,8908,8920,7140,5904,5323,5727,4819,3120,-1, //K->Kr
		1532,2630,4472,6511,8570,10280,11500,12370,12450,12023,10490,8650,7310,7310,6697,6240,4940,-1,//Rb->Xe
		1879,3510,6146,6689,6640,7010,7264,7353,5244,7901,8219,8551,8795,9066,9321,6570,9841,13310,16650,19250,21020,22610,22650,21090,19300,13534,11850,11340,9780,9196,-1,-1//Cs->Rn
	};
	
	for(i=0;i<mat->ne;i++) {
		w[i]=M[mat->e[i][0]-1]*(double)(mat->e[i][2]);
		Mtot+=w[i];
		zt[i]=mat->e[i][0];
		if(mat->e[i][1]<=0) at[i]=(int)(M[mat->e[i][0]-1]+0.5);
		else at[i]=mat->e[i][1];
		//printf("Z=%3d A=%3d w=%lf; ",zt[i],at[i],w[i]);
	}
	
	if(mat->r<0) r=Mtot*mat->P/(83.144621*mat->T); //in mg/cm??
	else {
		if(mat->ne==1) r=rho[mat->e[0][0]-1];
		else r=mat->r;
	}
	//printf("r=%lf\n",r);
	
	return r;
}

void integrale(tabededx *primo) {
	long double buff,intg;
	double e1,dd1,e2,dd2;
	tabededx *scorri=NULL;
	
	//CALCOLO L'INTEGRALE (TRAPEZI):
	buff=0; intg=0;
	e1=0;
	dd1=primo->dd;
	primo->dd=0;
	for(scorri=primo->next;scorri!=NULL;scorri=scorri->next) {
		e2=scorri->e;
		dd2=scorri->dd;
		buff+=(long double)((dd1+dd2)*(e2-e1)/2.);
		if(intg*0.01<buff) {
			intg+=buff;
			buff=0;
		}
		scorri->dd=(double)(intg+buff);
		e1=e2;
		dd1=dd2;
	}
	return;
}

void creatab(int Zp,MATE *mat,int *zt,int *at,double *w,double rho,tabdb *tab) {
	tabededx *primo=NULL,*scorri=NULL,*tmp;
	double emin,emax,dd0,e1,e2,ec,dd1,dd2,ddc;
	int pos;
	int ne=mat->ne,srim=0;
	FILE *f;
	double P0,P1,p2,p3;
	char fn[1000];
	
	//CREO SEMPRE TUTTE LE TABELLE!!
	//TABELLA DA NIST
	primo=NULL;
	sprintf(fn,"data/nist_%d_%d_%s.txt",(int)(M[Zp-1]+0.5),Zp,mat->nome);
	f=fopen(fn,"r");
	if(f!=NULL) {
		for(;fscanf(f," %lg %lg %lg ",&e1,&dd1,&dd2)==3;) {
			e1/=(double)((int)(M[Zp-1]+0.5));
			if(e1<EMIN) continue;
			dd1+=dd2;
			if(primo==NULL) {
				primo=(tabededx *)malloc(sizeof(tabededx));
				primo->e=0.;
				primo->dd=1./dd1;
				primo->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=primo->next;
			}
			else {
				scorri->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=scorri->next;
			}
			scorri->e=e1;
			scorri->dd=1./dd1;
			scorri->next=NULL;
		}
		fclose(f);
	}
	if(primo!=NULL) integrale(primo);
	tab->tnist=primo;
	
	//TABELLA DA SRIM
	primo=NULL;
	sprintf(fn,"data/srim_%d_%s.txt",Zp,mat->nome);
	f=fopen(fn,"r");
	if(f!=NULL) {
		for(;fscanf(f," %lg %lg %lg ",&e1,&dd1,&dd2)==3;) {
			//e1/=(double)((int)(M[Zp-1]+0.5));
			if(e1<EMIN) continue;
			dd1+=dd2;
			if(primo==NULL) {
				primo=(tabededx *)malloc(sizeof(tabededx));
				primo->e=0.;
				primo->dd=1./dd1;
				primo->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=primo->next;
			}
			else {
				scorri->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=scorri->next;
			}
			scorri->e=e1;
			scorri->dd=1./dd1;
			scorri->next=NULL;
		}
		fclose(f);
	}
	if(primo!=NULL) integrale(primo);
	tab->tsrim=primo;
	
	//TABELLA BARBUI (solo per gas, mylar ed elementi)
	primo=NULL;
	if((strcmp(mat->nome,"mylar")==0)||(mat->r<0)||(ne==1)) {
		sprintf(fn,"data/nist_1_1_%s.txt",mat->nome);
		f=fopen(fn,"r");
		if(f==NULL) {
			sprintf(fn,"data/srim_1_%s.txt",mat->nome);
			f=fopen(fn,"r");
			srim=1;
		}
	}
	else f=NULL;
	if(f!=NULL) {
		if(mat->r<0) { //GAS
			P0=1.468-0.08301*log((double)Zp);
			P1=4.615-0.71544*log((double)Zp);
			p2=0.4039;
			p3=0.2965;
		}
		else {
			P0=1.45-0.07*log((double)Zp);
			P1=6.;
			if(strcmp(mat->nome,"mylar")==0) { //MYLAR
				p2=0.5089;
				p3=0.62;
			}
			else { //ELEMENTO SOLIDO
				p2=0.4641-0.0011*zt[0];
				p3=0.6;
			}
		}
		//e2=1-P0*exp(-P1*pow(GTH,p2)/pow((double)Zp,p3));
		//dd2=P0*P1*p2*pow(GTH,p2-1)*exp(-P1*pow(GTH,p2)/pow((double)Zp,p3))/pow((double)Zp,p3);
		for(;fscanf(f," %lg %lg %lg ",&e1,&dd1,&dd2)==3;) {
			if(e1<EMIN) continue;
			if(srim) {
				if(strcmp(mat->nome,"CF4")==0) {
					if(e1<2) dd1*=1.06-0.0413*e1;
					else dd1*=0.9744;
				}
				if(strcmp(mat->nome,"mylar")==0) {
					if(e1<0.2) dd1*=1.0868;
					else  {
						if(e1<1.5268) dd1*=1.0398+0.3286*e1-0.5029*e1*e1+0.1726*e1*e1*e1;
						else dd1*=0.9835;
					}
				}
				if(strcmp(mat->nome,"C4H10")==0) {
					if(e1<0.1203) dd1*=0.7851+1.7868*e1;
				}
			}
			
			ec=1-P0*exp(-P1*pow(e1,p2)/pow((double)Zp,p3));
			
			if(ec>0) dd1*=pow(ec*(double)Zp,2.);
			else dd1=0.;
			
			if(e1<1.) {
				dd1=e1*dd1+(1.-e1)*dedx(ne,Zp,zt,at,w,rho,e1);
				//if(dd1/M[Zp-1]<0.05) dd1=M[Zp-1]*0.05;
			}
			
			//nuclear dedx
			dd1+=nudedx(ne,Zp,(int)(M[Zp-1]+0.5),zt,at,w,e1);
			
			if(primo==NULL) {
				primo=(tabededx *)malloc(sizeof(tabededx));
				primo->e=0.;
				primo->dd=1./dd1;
				primo->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=primo->next;
			}
			else {
				scorri->next=(tabededx *)malloc(sizeof(tabededx));
				scorri=scorri->next;
			}
			scorri->e=e1;
			scorri->dd=1./dd1;
			scorri->next=NULL;
		}
		fclose(f);
	}
	if(primo!=NULL) integrale(primo);
	tab->tbarb=primo;
	
	//TABELLA SCHWALM
	primo=NULL;
	emin=EMIN;
	emax=EMAX;
	dd0=1./dedx(ne,Zp,zt,at,w,rho,EMIN);
	primo=(tabededx *)malloc(sizeof(tabededx));
	primo->e=0;
	primo->dd=dd0;
	primo->next=(tabededx *)malloc(sizeof(tabededx));
	tmp=primo->next;
	tmp->e=emin;
	tmp->dd=dd0;
	tmp->next=(tabededx *)malloc(sizeof(tabededx));
	tmp=tmp->next;
	ec=(emax+emin)/2.;
	tmp->e=ec;
	tmp->dd=1./dedx(ne,Zp,zt,at,w,rho,ec);
	tmp->next=(tabededx *)malloc(sizeof(tabededx));
	tmp=tmp->next;
	tmp->e=emax;
	tmp->dd=1./dedx(ne,Zp,zt,at,w,rho,EMAX);
	tmp->next=NULL;
	for(scorri=primo->next;scorri->next->next!=NULL;scorri=scorri->next) {
		ec=scorri->next->e; ddc=scorri->next->dd;
		e1=scorri->e; dd1=scorri->dd;
		e2=scorri->next->next->e; dd2=scorri->next->next->dd;
		for(;fabs(ddc-dd1-(ec-e1)*(dd2-dd1)/(e2-e1))>PREC;) {
			if(ec-e1<e2-ec) pos=1;
			else pos=0;
			if(pos==0 && ec<emin) pos=1;
			tmp=(tabededx *)malloc(sizeof(tabededx));
			tmp->e=(pos?((ec+e2)/2.):((ec+e1)/2.));
			tmp->dd=1./dedx(ne,Zp,zt,at,w,rho,tmp->e);
			if(pos) {
				tmp->next=scorri->next->next;
				scorri->next->next=tmp;
			}
			else {
				tmp->next=scorri->next;
				scorri->next=tmp;
			}
			ec=scorri->next->e; ddc=scorri->next->dd;
			e2=scorri->next->next->e; dd2=scorri->next->next->dd;
		}
	}
	if(primo!=NULL) integrale(primo);
	tab->tschw=primo;
	
	return;
}

//form: 0=auto
//      1=nist
//      2=barbui
//      3=srim
//      4=schwalm

double elosscore(int opz,double ein,int Zp,MATE *mat,double thickness,int form) {
	int zt[5],at[5];
	double w[5],rho,rhosch,t=0,tres,e1,e2,e3,r1,r2,r3,a,b,c,d,range,ee1,ee2;
	tabededx *primo=NULL,*scorri;
	tabdb *t1;
	
	if(ein<0||Zp<0||mat==NULL||thickness<0) {
		printf("[elosscore] ERROR: un parametro inserito ?? negativo!\n");
		exit(1);
	}
	
	rho=leggimate(mat,zt,at,w);
	rhosch=rho;
	if(opz!=4) {
		t=rho*thickness/10000.; //in mg/cm??
		if(mat->ne>1) rhosch=100000.;
	}
	if(db==NULL) {
		//printf("[elosscore] Creazione tabella perdite per Z=%d in %s\n",Zp,mat->nome);
		db=(tabdb *)malloc(sizeof(tabdb));
		strcpy(db->nome,mat->nome);
		db->z=Zp;
		creatab(Zp,mat,zt,at,w,rhosch,db);
		db->next=NULL;
		t1=db;
	}
	else {
		for(t1=db;t1->next!=NULL;t1=t1->next) {
			if((strcmp(mat->nome,t1->nome)==0)&&(Zp==t1->z)) break;
		}
		if((strcmp(mat->nome,t1->nome)!=0)||(Zp!=t1->z)) {
			//printf("[elosscore] Creazione tabella perdite per Z=%d in %s\n",Zp,mat->nome);
			t1->next=(tabdb *)malloc(sizeof(tabdb));
			strcpy(t1->next->nome,mat->nome);
			t1->next->z=Zp;
			creatab(Zp,mat,zt,at,w,rhosch,t1->next);
			t1->next->next=NULL;
			t1=t1->next;
		}
	}
	switch(form) {
		case 0: case 1:
			primo=t1->tnist;
			if(primo!=NULL) break;
			if(form==1) printf("[elosscore] WARNING: impossibile creare la tabella NIST per Z=%d in %s\n",Zp,mat->nome);
		case 2:
			if((mat->r<0)||(strcmp(mat->nome,"mylar")==0)||(form==2)) primo=t1->tbarb;
			else primo=t1->tsrim;
			if(primo!=NULL) break;
			if(form==2) printf("[elosscore] WARNING: impossibile creare la tabella Barbui per Z=%d in %s\n",Zp,mat->nome);
		case 3:
			primo=t1->tsrim;
			if(primo!=NULL) break;
			if(form==3) printf("[elosscore] WARNING: impossibile creare la tabella SRIM per Z=%d in %s\n",Zp,mat->nome);
		case 4:
			primo=t1->tschw;
	}
	if(primo==NULL) {
		printf("[elosscore] ERROR: non ho tabelle e Schwalm non funziona per Z=%d in %s\n",Zp,mat->nome);
		return -1;
	}
	if(opz==3) {
		ein*=rho/10000.;
		for(scorri=primo;scorri->next->next!=NULL;scorri=scorri->next) {
			if(fabs(ein-scorri->next->dd)<fabs(ein-scorri->next->next->dd)) break;
			if(scorri->next->next->next==NULL) break;
		}
		if(ein>scorri->next->next->dd) {
			printf("[elosscore] ERROR: too high input range!\n");
			return -1;
		}
		e1=scorri->dd;  e2=scorri->next->dd;  e3=scorri->next->next->dd;
		r1=scorri->e; r2=scorri->next->e; r3=scorri->next->next->e;
	}
	else {
		for(scorri=primo;scorri->next->next!=NULL;scorri=scorri->next) {
			if(fabs(ein-scorri->next->e)<fabs(ein-scorri->next->next->e)) break;
			if(scorri->next->next->next==NULL) break;
		}
		if(ein>scorri->next->next->e) {
			printf("[elosscore] ERROR: too high input energy!\n");
			return -1;
		}
		e1=scorri->e;  e2=scorri->next->e;  e3=scorri->next->next->e;
		r1=scorri->dd; r2=scorri->next->dd; r3=scorri->next->next->dd;
	}
	d=e1*e1*(e2-e3)+e2*e2*(e3-e1)+e3*e3*(e1-e2);
	a=r1*(e2-e3)+r2*(e3-e1)+r3*(e1-e2);
	b=e1*e1*(r2-r3)+e2*e2*(r3-r1)+e3*e3*(r1-r2);
	c=e1*e1*(e2*r3-e3*r2)+e2*e2*(e3*r1-e1*r3)+e3*e3*(e1*r2-e2*r1);
	
	range=(a*ein*ein+b*ein+c)/d;
	
	if(opz==2) return range*10000./rho;
	if(opz==3) return range;
	
	if(opz==4) {
		for(scorri=primo;scorri->next->next!=NULL;scorri=scorri->next) {
			if(fabs(thickness-scorri->next->e)<fabs(thickness-scorri->next->next->e)) break;
			if(scorri->next->next->next==NULL) break;
		}
		if(e1!=scorri->e) {
			e1=scorri->e;  e2=scorri->next->e;  e3=scorri->next->next->e;
			r1=scorri->dd; r2=scorri->next->dd; r3=scorri->next->next->dd;
			d=e1*e1*(e2-e3)+e2*e2*(e3-e1)+e3*e3*(e1-e2);
			a=r1*(e2-e3)+r2*(e3-e1)+r3*(e1-e2);
			b=e1*e1*(r2-r3)+e2*e2*(r3-r1)+e3*e3*(r1-r2);
			c=e1*e1*(e2*r3-e3*r2)+e2*e2*(e3*r1-e1*r3)+e3*e3*(e1*r2-e2*r1);
		}
		tres=(a*thickness*thickness+b*thickness+c)/d;
		return (range-tres)*10000./rho;;
	}
	
	if(opz==0) {
		if(t>=range) return 0;
		tres=range-t;
		scorri=primo;
	}
	else { //SOLO opz==1
		tres=range+t;
	}
	
	for(;scorri->next->next!=NULL;scorri=scorri->next) {
		if(fabs(tres-scorri->next->dd)<fabs(tres-scorri->next->next->dd)) break;
		if(scorri->next->next->next==NULL) break;
	}
	if(e1!=scorri->e) {
		e1=scorri->e;  e2=scorri->next->e;  e3=scorri->next->next->e;
		r1=scorri->dd; r2=scorri->next->dd; r3=scorri->next->next->dd;
		d=e1*e1*(e2-e3)+e2*e2*(e3-e1)+e3*e3*(e1-e2);
		a=r1*(e2-e3)+r2*(e3-e1)+r3*(e1-e2);
		b=e1*e1*(r2-r3)+e2*e2*(r3-r1)+e3*e3*(r1-r2);
		c=e1*e1*(e2*r3-e3*r2)+e2*e2*(e3*r1-e1*r3)+e3*e3*(e1*r2-e2*r1);
	}
	
	
	if(fabs(a/d)<1.e-10) return (tres-c)/b;
	
	ee1=b*(-1.+sqrt(1-4*a*(c-tres*d)/(b*b)))/(2.*a);
	ee2=b*(-1.-sqrt(1-4*a*(c-tres*d)/(b*b)))/(2.*a);
	
	return ((fabs(ee1-e2)<fabs(ee2-e2))?ee1:ee2);
}

void mateveda(MATE *mat,int *materiale,float *atloc,float *pressione) {
	char nmat[29][10]={
		"Si","C10H8O4","CH2","Ni","C3F8","C","Ag","Sn","CsI","Au",
		"U","air","Nb","Ta","V","CF4","C4H10","Al","Pb","PbS",
		"KCl","Ge","Ca","Cu","Ti","Bi","Mg","Li","Zn"
	};
	
	//aggiunta del 7/6/16 per associare all'elemento giusto gli ioni in cui A ?? specificato (es: 40Ca -> Ca)
	int istart=0;
	for(;mat->nome[istart]>='0'&&mat->nome[istart]<='9';istart++);
	
	for(*materiale=0;*materiale<29;(*materiale)++) {
		if(strcmp(mat->nome+istart,nmat[*materiale])==0) break;
	}
	if(*materiale==29) *materiale=10;
	else (*materiale)++;
	if(strcmp(mat->nome,"mylar")==0) *materiale=2;
	if(strcmp(mat->nome,"NE102")==0) *materiale=3;
	if(strcmp(mat->nome,"isobutane")==0) *materiale=17;
	
	*pressione=0;
	if(mat->e[0][1]>=0) *atloc=(float)(mat->e[0][1]);
	else {
		if((mat->e[0][0]>=1)&&(mat->e[0][0]<=86)) *atloc=(float)(M[mat->e[0][0]-1]);
		else *atloc=0;
	}
	switch(*materiale) {
		case 5: case 12: case 16: case 17:
			*pressione=(float)(mat->P); //niente break!! Anche in questi casi atloc=0!
		case 2: case 3: case 9: case 20: case 21:
			*atloc=0;
	}
	return;
}

double eloss(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in ??m*/,int form) {
	if(form<5) return (double)Ap*elosscore(0,ein/(double)Ap,Zp,mat,thickness/(double)Ap,form);
	else {
		float e1,el1;
		float epart=(float)ein;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess=(float)thickness;
		int idir=1,icod,materiale;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		ecorr_veda_(&epart,&zpart,&apart,&atloc,&e1,&el1,&materiale,&spess,&idir,&icod,&pressione);
		double eout=(double)e1;
		return ((eout>ein)?ein:eout);
	}
}

double eorig(double epass /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in ??m*/,int form) {
	if(form<5) return (double)Ap*elosscore(1,epass/(double)Ap,Zp,mat,thickness/(double)Ap,form);
	else {
		float e1,el1;
		float epart=(float)epass;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess=(float)thickness;
		int idir=2,icod,materiale;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		ecorr_veda_(&epart,&zpart,&apart,&atloc,&e1,&el1,&materiale,&spess,&idir,&icod,&pressione);
		double ein=(double)e1;
		return ((ein<epass)?epass:ein);
	}
}

double e2range(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,int form) {
	if(form<5) return (double)Ap*elosscore(2,ein/(double)Ap,Zp,mat,0,form);
	else {
		float e1,el1,e2,el2;
		float epart=(float)ein;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess1=20000.,spess2,step=10000.;
		int idir=1,icod,materiale,n;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		for(n=0;n<ITER;n++) {
			spess2=spess1+step;
			ecorr_veda_(&epart,&zpart,&apart,&atloc,&e1,&el1,&materiale,&spess1,&idir,&icod,&pressione);
			ecorr_veda_(&epart,&zpart,&apart,&atloc,&e2,&el2,&materiale,&spess2,&idir,&icod,&pressione);
			if(e1>0) {
				if(e2>0) spess1+=step;
				else {
					if(step/spess1<0.001) break;
					step/=2.;
				}
			}
			else {
				while(spess1<1.5*step) step/=2;
				spess1-=step;
			}
		}
		if(n>=ITER) printf("[elosscore] WARNING: too many iterations while searching range\n");
		return (((double)spess1)+((double)spess2))/2.;
	}
}

double range2e(double range /*in ??m*/,int Zp,int Ap,MATE *mat,int form) {
	if(form<5) return (double)Ap*elosscore(3,range/(double)Ap,Zp,mat,0,form);
	else {
		float e1,el1,e2,el2;
		float epart1=200.,epart2,step=100.;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess=(float)range;
		int idir=1,icod,materiale,n;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		for(n=0;n<ITER;n++) {
			epart2=epart1+step;
			ecorr_veda_(&epart1,&zpart,&apart,&atloc,&e1,&el1,&materiale,&spess,&idir,&icod,&pressione);
			ecorr_veda_(&epart2,&zpart,&apart,&atloc,&e2,&el2,&materiale,&spess,&idir,&icod,&pressione);
			if(e1>0) {
				while(epart1<1.5*step) step/=2;
				epart1-=step;
			}
			else {
				if(e2>0) {
					if(step/epart1<0.001) break;
					step/=2.;
				}
				else epart1+=step;
			}
		}
		if(n>=ITER) printf("[elosscore] WARNING: too many iterations while searching range\n");
		return (((double)epart1)+((double)epart2))/2.;
	}
}

double e2spes(double ein /*in MeV*/,int Zp,int Ap,MATE *mat,double epass /*in MeV*/,int form) {
	if(form<5) return (double)Ap*elosscore(4,ein/(double)Ap,Zp,mat,epass/(double)Ap,form);
	else {
		float e1,el1,e2,el2;
		float epart=(float)ein;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess1=20000.,spess2,step=10000.;
		int idir=1,icod,materiale,n;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		for(n=0;n<ITER;n++) {
			spess2=spess1+step;
			ecorr_veda_(&epart,&zpart,&apart,&atloc,&e1,&el1,&materiale,&spess1,&idir,&icod,&pressione);
			ecorr_veda_(&epart,&zpart,&apart,&atloc,&e2,&el2,&materiale,&spess2,&idir,&icod,&pressione);
			if((double)e1>epass) {
				if((double)e2>epass) spess1+=step;
				else {
					if(step/spess1<0.001) break;
					step/=2.;
				}
			}
			else {
				while(spess1<1.5*step) step/=2;
				spess1-=step;
			}
		}
		if(n>=ITER) printf("[elosscore] WARNING: too many iterations while searching range\n");
		return (((double)spess1)+((double)spess2))/2.;
	}
}

double spes2e(double el /*in MeV*/,int Zp,int Ap,MATE *mat,double thickness /*in ??m*/,int form) {
	double ein=100.;
	if(form<5) {
		double epass;
		double step=100.;
		int n;
		
		epass=eloss(ein,Zp,Ap,mat,thickness,form);
		//printf("  0) ein = %lf  el_st = %lf  el_in = %lf  step = %lf\n",ein,ein-epass,el,step);
		for(n=0;(fabs(ein-epass-el)>0.001)&&(n<ITER);n++) {
			if((ein-epass-el>0)||(epass==0)) {
				if(step<0) step*=-0.5;
			}
			else {
				if(step>0) step*=-0.5;
				while(ein+step<0.01) step*=0.5;
			}
			ein+=step;
			epass=eloss(ein,Zp,Ap,mat,thickness,form);
			//printf("%3d) ein = %lf  el_st = %lf  el_in = %lf  step = %lf\n",n+1,ein,ein-epass,el,step);
		}
		if(n>=ITER) printf("[elosscore] WARNING: too many iterations while searching energy\n");
		if(fabs(ein-epass-el)>0.001) {printf("[elosscore] WARNING: energia persa inserita troppo elevata\n"); return -ein;}
	}
	else {
		float epart;
		float de=(float)el;
		float zpart=(float)Zp;
		float apart=(float)Ap;
		float spess=(float)thickness;
		int materiale;
		float atloc,pressione;
		mateveda(mat,&materiale,&atloc,&pressione);
		de_vedaloss_(&zpart,&apart,&atloc,&de,&spess,&epart,&materiale,&pressione);
		if(epart>0) ein=(double)epart;
		else ein=-range2e(thickness,Zp,Ap,mat,form);
	}
	return ein;
}
