#ifndef GEO
#define GEO
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <TGraph.h>
#include "perdenc.h"

typedef struct {
	double ei,e1,e2,e3,lo;
} energie;

typedef struct Nodo {
	int ne;
	int min;
	void *son;
} nodo;

using namespace std;

// fatto con la classe container
class Classe_geo{
 public:
  //<<<<<<< Classe_geo.h

  int ThereIsPhos;
  int ThereIsHector;
  int ThereIsGarf;
  int ThereIsRCO;
  int ThereIsBlocchiFazia;

  double D[111][163],B[111][163];

  struct PHOS
  {
    float dist[6][9];
    float theta[6][9];
    float phi[6][9];
    float spess[6][9];
  float veccentro[6][9][3];     

    float vi[6][9][3];//versore sul piano del riv.
    float vj[6][9][3];//versore sul piano del riv.
    float vk[6][9][3];//versore ortogonale al riv e passante per il suo centro (NB APPROSSIMAZIONE!)
     static constexpr float dx=64;
     static constexpr float  dy=64;
    
    int codice[6][9];    
 
  }phos;
  float vfascio[3];
  float assex[3];
  float assey[3];

struct GARF
{
  float thecsi[8];
  float themincsi[8];
  float themaxcsi[8];
  float dist[8];
  float phi[8][24];

  float veccentro[24][8][3];     

    float vi[24][8][3];//versore sul piano del riv.
    float vj[24][8][3];//versore sul piano del riv.
    float vk[24][8][3];//versore ortogonale al riv e passante per il suo centro


  static constexpr float dphi=15; 
  //  static constexpr float dphi_corr=13.5;
  static constexpr float dphi_corr=14.5;

  static constexpr float spess_fasciatura=0.;//0.25;//mm

  float pressione;
  float pressione_back;
  static constexpr float alpha0=8;//inclinazione del piano di garfield su cui sono impiantate le strip rispetto a 90 gradi
  static constexpr float mylar1=6;//micron mylar di ingresso
  static constexpr float dead1=11500;//micron di gas morto 1
  static constexpr float strip1=36000;//micron lunghezza strip 1 
  static constexpr float strip2=35000;//micron lunghezza strip 2 

 static constexpr float dead2=4000;//micron di gas morto 2
  
  static constexpr float braccio=60400;//micron lunghezza del percorso dal target all'ingresso in garfield
 
static constexpr float mylar2=3.5;//micron mylar di uscita + mylar cesio
//   vector <float> Edown;//estrip1
//   vector <float> Eup;//estrip2
//   vector <float> Ecsi;//ecsi
//   vector <int> indiciz;
//   vector <int> indicie;
//   vector <int> indicitheta;

// new table organization
  
  nodo *TabGarf;
  

  int codice[8][24];
  float bl_rej[8][24];
  float cut_csi_p[8][24];
  float cut_csi_d[8][24];
  float cut_csi_t[8][24];
  float cut_csi_3he[8][24];
  float cut_csi_a[8][24];
  TGraph *gcsi[8];

  float h[8],l1[8],l2[8];
  float hdown[8],hup[8];
  int micro_rotte[24][8][4];

}garf;

 struct RCO{
   float strip_dist[8][8];//settore strip
   float strip_themin[8][8];
   float strip_themax[8][8];
   float strip_the[8][8];
   float strip_phi[8][8];
   float strip_phimin[8][8];
   float strip_phimax[8][8];
   float strip_spess[8][8];
   float strip_veccentro[3];
   static constexpr float dphi=45; //apertura in phi di un settore
   float gas_dist[8];//settore
   float gas_themin[8];
   float gas_themax[8];
   float gas_the[8];
   float gas_phi[8];
   float gas_phimin[8];
   float gas_phimax[8];  
   float gas_veccentro[3];
   float gas_vk[8][3];
   float gas_vj[8][3];
   float gas_vi[8][3];
   float csi_dist[8][6];
   float csi_phi[8][6];
   float csi_phimin[8][6];
   float csi_phimax[8][6];
   float csi_themin[8][6];
   float csi_themax[8][6];
   float csi_the[8][6];
   float csi_veccentro[3];
   TGraph *csi[6];

   float gas_a0;
   float gas_a8;
   float gas_acc;
   float gas_gcc;
   float gas_xa;
   float gas_ya;
   float gas_yg;
   float strip_a8;
   float strip_gcc;
   float strip_yg;
   float csi_a0;
   float csi_acc;
   float csi_xa;
   float csi_ya;
   float csi_xb;
   float csi_yb;
   float csi_xA;
   float csi_yA;
   float csi_xB;
   float csi_yB;
   float csi_a2;
   float csi_ccc;
   float csi_xc;
   float csi_yc;
   float csi_xd;
   float csi_yd;
   float csi_xC;
   float csi_yC;
   float csi_xD;
   float csi_yD;
   float csi_a4;
   float csi_ecc;
   float csi_xe;
   float csi_ye;
   float csi_xf;
   float csi_yf;
   float csi_xE;
   float csi_yE;
   float csi_xF;
   float csi_yF;
   float csi_a8;
   float csi_gcc;
   float csi_xg;
   float csi_yg;
   float csi_xh;
   float csi_yh;
   float csi_xH;
   float csi_yH;
   float csi_xG;
   float csi_yG;

   float pressione;
   static constexpr float spesgas=30000;//micron spessore mezza camera
   static constexpr float spesmylar=3.;//micron spessore finestra ingresso e uscita della camera
   static constexpr float spesmylar_mezzo=1.5;//micron spessore del catodo centrale della camera
   static constexpr float spesmylarcsi=1.5;//micron spessore finestra ingresso cesi
   //      static constexpr float dphimorto=0.26;//gradi morti da entrambe le parti di un settore (2.5mm)
   // static constexpr float dphimortocsi=0.13;//gradi morti fra i cesi dietro lo stesso settore (1.4mm)
   //spazio in theta fra i cesi 1.4mm

      static constexpr float dx_spe=0.8;//mm spessore morto a dx fra un settore e l'altro
    static constexpr float sx_spe=1.7;//mm spessore morto a sx fra un settore e l'altro
   static constexpr float spe_morto_csi=0.7;//mm spessore morto fra i cesi dietro lo stesso settore 


   int codice_gas[8];//0 buono, 1 rotto
   int codice_strip[8][9];
   int codice_csi[8][7];
   float soglie_csi_p[8][6];
   float soglie_csi_a[8][6];
   
   
   int stripcopertamax;
   
   nodo *TabRCO;

   int zinf[8][8][6];
   int zsup[8][8][6];
   int ainf[8][8][6][100];
   int asup[8][8][6][100];
   int intervalli;
   
 }rco;

 struct HECTOR{
   float theta[8];
   float phi[8];
   float dist[8];
   static constexpr float raggio=14.6/2;
   float veccentro[8][3];
   float vi[8][3];
   float vj[8][3];
   float vk[8][3];

 }hector;

 struct FAZIETTO{
   float dist[12];//blocchi
   float theta[12];
   float phi[12];
   //Attenzione: si fa una geo approssimata; si suppone che i riv siano su un piano ortogonale alla congiungente centro del blocco-target
   float vkblocco[12][3];
   float vjblocco[12][3];
   float viblocco[12][3];
   float vcentroblocco[12][3];
   float vx[12][16];//posizioni dei centri dei singoli riv rispetto a un SDR centrato nel centro del blocco
   float vy[12][16];
   //schema della numerazione dei riv nel blocco
   //0 1   |  4 5     Q1T1 Q1T2 | Q2T1 Q2T2
   //3 2   |  7 6     Q1T4 Q1T3 | Q2T4 Q2T3
   //-----------      ----------------------
   //12 13 | 8 9      Q4T1 Q4T2 | Q3T1 Q3T2
   //15 14 |11 10     Q4T4 Q4T3 | Q3T4 Q3T3

 int rotti[12][4][4];
 float spes_si1[12][4][4];//spessore silicio 1 del telescopio, micron
 float spes_si2[12][4][4];//spessore silicio 2 del telescopio, micron

   
           static constexpr float dx_active=20;
    static constexpr float  dy_active=20;
    

    
     static constexpr float dx=22;
     static constexpr float  dy=22;

   static constexpr float spes_csi=100000;//spessore csi del telescopio, micron

   int alimSi1Si2[12][4][4]; //limiti di Z per l'identificazione in massa
   int alimSi2CsI[12][4][4];
   int aliminfSi1PSA[12][4][4];
   int alimsupSi1PSA[12][4][4];
   float zesi1inf[12][4][4][100];
   float zesi1sup[12][4][4][100];
   float aesi1inf[12][4][4][100];
   float aesi1sup[12][4][4][100];

   int infmassesi2csi[12][4][4][100];
   int supmassesi2csi[12][4][4][100];
   int infmassesi1si2[12][4][4][100];
   int supmassesi1si2[12][4][4][100];
   int infmassepsa[12][4][4][100];
   int supmassepsa[12][4][4][100];


 }fazietto;

struct TARGET
{
  float spessore;
  int zt,at;
  int materiale; //da rimuovere una volta abolita vedaloss.
  float densita;
  MATE mat;
} target;

 int ncuts;

  static Classe_geo *geo;


  Classe_geo()
    {
      fHaveTable = 0;
      fHaveTableRCO = 0;

      for(int ip=0;ip<6;ip++)
	{
	  for(int ih=0;ih<9;ih++)
	    {
	      phos.dist[ip][ih]=-1;
	      phos.theta[ip][ih]=-1;
	      phos.phi[ip][ih]=-1;
	      phos.codice[ip][ih]=-1;
	    }
	}
      garf.h[0]=45;
      garf.h[1]=43;
      garf.h[2]=41;
      garf.h[3]=41;
      garf.h[4]=41;
      garf.h[5]=41;
      garf.h[6]=43;
      garf.h[7]=45;

      // Attenzione: i cesi 5-8 sono messi con la base piu' stretta verso il theta=0, mentre i cesi 1-4 hanno la base piu' larga verso il theta=0

      garf.l2[0]=26.3;
      garf.l2[1]=32.2;
      garf.l2[2]=35.3;
      garf.l2[3]=35.4;

      garf.l1[4]=35.4;
      garf.l1[5]=35.3;
      garf.l1[6]=32.2;
      garf.l1[7]=26.3;

      garf.l1[0]=36.0;
      garf.l1[1]=40.0;
      garf.l1[2]=40.7;
      garf.l1[3]=38.3;

      garf.l2[4]=38.3;
      garf.l2[5]=40.7;
      garf.l2[6]=40.0;
      garf.l2[7]=36.0;

      garf.hdown[0]=2.66;
      garf.hdown[1]=2.82;
      garf.hdown[2]=0.905;
      garf.hdown[3]=1.42;
      garf.hdown[4]=1.42;
      garf.hdown[5]=-0.6;
      garf.hdown[6]=1.15;
      garf.hdown[7]=2.3;

      garf.hup[0]=2.66;
      garf.hup[1]=2.85;
      garf.hup[2]=0.92;
      garf.hup[3]=1.42;
      garf.hup[4]=1.42;
      garf.hup[5]=2.43;
      garf.hup[6]=4.53;
      garf.hup[7]=3.05;

      for(int ip=0;ip<8;ip++)
	{
	  garf.thecsi[ip]=-1;
	  garf.themincsi[ip]=-1;
	  garf.themaxcsi[ip]=-1;
	  garf.dist[ip]=-1;
	  garf.gcsi[ip]=0;
	  for(int isec=0;isec<24;isec++)
	    {
	      garf.phi[ip][isec]=-1;
	      garf.codice[ip][isec]=0;
		  garf.bl_rej[ip][isec]=1;
		  garf.cut_csi_3he[ip][isec]=0;
		  garf.cut_csi_a[ip][isec]=0;
		  garf.cut_csi_p[ip][isec]=0;
		  garf.cut_csi_d[ip][isec]=0;
		  garf.cut_csi_t[ip][isec]=0;
		  for(int icsi=0;icsi<8;icsi++)
		    {
		      for(int ik=0;ik<4;ik++)
			{
		      garf.micro_rotte[isec][icsi][ik]=0;
			}
		    }	
	    }
	  hector.dist[ip]=-1;
	  hector.theta[ip]=-1;
	  hector.phi[ip]=-1;
	  //<<<<<<< Classe_geo.h
	  for(int istrip=0;istrip<8;istrip++)
	    {
	      rco.strip_dist[ip][istrip]=-1;
	      rco.strip_themin[ip][istrip]=-1;
	      rco.strip_themax[ip][istrip]=-1;
	      rco.strip_the[ip][istrip]=-1;
	      rco.strip_phimin[ip][istrip]=-1;
	      rco.strip_phimax[ip][istrip]=-1;
	      rco.strip_phi[ip][istrip]=-1;
	    }
	  rco.gas_dist[ip]=-1;
	  rco.gas_phimin[ip]=-1;
	  rco.gas_phimax[ip]=-1;
	  rco.gas_phi[ip]=-1;
	  rco.gas_themin[ip]=-1;
	  rco.gas_themax[ip]=-1;
	  rco.gas_the[ip]=-1;
	  for(int ic=0;ic<6;ic++)
	    {
	      rco.csi_dist[ip][ic]=-1;
	      rco.csi_the[ip][ic]=-1;
	      rco.csi_themin[ip][ic]=-1;
	      rco.csi_themax[ip][ic]=-1;
	      rco.csi_phi[ip][ic]=-1;
	      rco.csi_phimin[ip][ic]=-1;
	      rco.csi_phimax[ip][ic]=-1;

	    }
	}
      for(int j=0;j<3;j++)
	{
	  rco.gas_veccentro[j]=0;
	  rco.csi_veccentro[j]=0;
	  rco.strip_veccentro[j]=0;

	}
      for(int j=0;j<6;j++)
	{
	  rco.csi[j]=0;
	}
   rco.gas_a0=-1;
   rco.gas_a8=-1;
   rco.gas_acc=-1;
   rco.gas_gcc=-1;
   rco.gas_xa=-1;
   rco.gas_ya=-1;
   rco.gas_yg=-1;
   rco.strip_a8=-1;
   rco.strip_gcc=-1;
   rco.strip_yg=-1;

   rco.csi_a0=-1;
   rco.csi_acc=-1;
   rco.csi_xa=-1;
   rco.csi_ya=-1;
   rco.csi_xb=-1;
   rco.csi_yb=-1;
   rco.csi_xA=-1;
   rco.csi_yA=-1;
   rco.csi_xB=-1;
   rco.csi_yB=-1;
   rco.csi_a2=-1;
   rco.csi_ccc=-1;
   rco.csi_xc=-1;
   rco.csi_yc=-1;
   rco.csi_xd=-1;
   rco.csi_yd=-1;
   rco.csi_xC=-1;
   rco.csi_yC=-1;
   rco.csi_xD=-1;
   rco.csi_yD=-1;
   rco.csi_a4=-1;
   rco.csi_ecc=-1;
   rco.csi_xe=-1;
   rco.csi_ye=-1;
   rco.csi_xf=-1;
   rco.csi_yf=-1;
   rco.csi_xE=-1;
   rco.csi_yE=-1;
   rco.csi_xF=-1;
   rco.csi_yF=-1;
   rco.csi_a8=-1;
   rco.csi_gcc=-1;
   rco.csi_xg=-1;
   rco.csi_yg=-1;
   rco.csi_xh=-1;
   rco.csi_yh=-1;
   rco.csi_xH=-1;
   rco.csi_yH=-1;
   rco.csi_xG=-1;
   rco.csi_yG=-1;

   rco.stripcopertamax=-1;

      for(int j=0;j<12;j++)
	{
	  fazietto.dist[j]=-1;
	  fazietto.theta[j]=-1;
	  fazietto.phi[j]=-1;
	    for(int iqua=0;iqua<4;iqua++)
		{
		  for(int itele=0;itele<4;itele++)
		    {
		      fazietto.rotti[j][iqua][itele]=1+2+4;//2^0+2^1+2^2
		      fazietto.spes_si1[j][iqua][itele]=300;
		      fazietto.spes_si2[j][iqua][itele]=500;
		      fazietto.alimSi1Si2[j][iqua][itele]=-1;
		      fazietto.alimSi2CsI[j][iqua][itele]=-1;
		      fazietto.aliminfSi1PSA[j][iqua][itele]=-1;
		      fazietto.alimsupSi1PSA[j][iqua][itele]=-1;
		      for(int iz=0;iz<100;iz++)
			{
		      fazietto.zesi1inf[j][iqua][itele][iz]=-1;
		      fazietto.zesi1sup[j][iqua][itele][iz]=-1;
		      fazietto.aesi1inf[j][iqua][itele][iz]=-1;
		      fazietto.aesi1sup[j][iqua][itele][iz]=-1;
		      fazietto.infmassesi2csi[j][iqua][itele][iz]=-1;
		      fazietto.supmassesi2csi[j][iqua][itele][iz]=-1;
		      fazietto.infmassesi1si2[j][iqua][itele][iz]=-1;
		      fazietto.supmassesi1si2[j][iqua][itele][iz]=-1;
		      fazietto.infmassepsa[j][iqua][itele][iz]=-1;
		      fazietto.supmassepsa[j][iqua][itele][iz]=-1;

			}
		    }
		}

	    
	}
      for(int j=0;j<8;j++)
	{
	  rco.codice_gas[j]=0;
	  for(int k=0;k<6;k++)
	    {

	      rco.soglie_csi_p[j][k]=0;
	      rco.soglie_csi_a[j][k]=0;
	    }
	  for(int k=0;k<9;k++)
	    {
	      rco.codice_strip[j][k]=0;
	    }
	  for(int k=0;k<7;k++)
	    {
	      rco.codice_csi[j][k]=0;
	    }
	}
      
      rco.intervalli=0;
      for(int isec=0;isec<8;isec++)
	{
	  for(int istrip=0;istrip<8;istrip++)
	    {
	      for(int icsi=0;icsi<6;icsi++)
		{
		  rco.zinf[isec][istrip][icsi]=1;
		  rco.zsup[isec][istrip][icsi]=0;
		  for(int iz=0;iz<100;iz++)
		    {
		  rco.ainf[isec][istrip][icsi][iz]=0;
		  rco.asup[isec][istrip][icsi][iz]=0;
		    }
		}
	    }
	}
      
      vfascio[0]=0.;
      vfascio[1]=0.;
      vfascio[2]=1.;
      assex[0]=1.;
      assex[1]=0.;
      assex[2]=0.;
      assey[0]=0.;
      assey[1]=1.;
      assey[2]=0.;
      target.materiale=-1;
      //<<<<<<< Classe_geo.h
      ThereIsPhos=0;
      ThereIsGarf=0;
      ThereIsHector=0;
      ThereIsRCO=0;
      ThereIsBlocchiFazia=0;

      ncuts=0;

        for(int j=0;j<111;j++)
    {
      for(int k=0;k<163;k++)
	{
	  D[j][k]=-100000;
	  B[j][k]=-100000;
	}
    }

    }

static void prova();

 int HaveTable(){return fHaveTable;};
 int HaveTableRCO(){return fHaveTableRCO;};

 void Leggigeo(char *filegeo);
 static Classe_geo *Getgeo();
 int InsidePhos(float thetapart,float phipart,float vpart, float *tvolo,int *code);
int InsideGarf(float thetapart,float phipart,float zpart,int *code,int *code_micro);
int InsideRCO(float thetapart,float phipart,float zpart,int *code);
 void vedaloss();
 void ecorr_veda(float *eingresso,float *zpr,float *apr,float *atar, float *eout,float *elost,int *mate,float *thick,int *idir, int *icod, float *pressione); 
int Soglie_Garfield(float zpart,float apart,float epart,float thetapart,int codice,float *estrip1,float *estrip2,float *ecsi, float *luce,float *aout, float *zout,float *eout,int *iginocchio,int code_micro);
//<<<<<<< Classe_geo.h
 int Soglie_RCO(float zpart,float apart,float epart,float thetapart,float phipart,int *codice,float *egas,float *esi,float *ecsi, float *luce,float *aout, float *zout,float *eout, float *qf);
 int Soglie_Phos(float zpart, float apart, float thetapart, float epart,float tvolo,float codephos,float *zout,float *aout, float *tout);

void Spalma_Garfield(float codice, float *theout, float *phiout,int code_micro);
 void Spalma_Phos(float codice, float *theout, float *phiout, float *dout);
 void Spalma_Fazietto(float codice, float *theout, float *phiout, float *dout);

 void TabellaPerdite(char *fileperdite);
 void TabellaPerditeBack(char *fileperdite);

float CalcoloEGarf(int z,int a,float estrip1,float estrip2, float ecsi,int i0);
 void Spalma_Hector(int codice,float *theout,float *phiout,float *dout);
 int InsideHector(float thetapart,float phipart,float *tvolo, int *code);
 int InsideFazietto(float thetapart,float phipart,float vpart,float *distanza,int *code);
 int Soglie_Hector(float zpart, float apart, float epart,float tvolo,float codephos,float *zout,float *aout, float *tout,float *eout);
//<<<<<<< Classe_geo.h
 void SpalmaRCO(int isec,int istrip,int icsi,float *theout,float *phiout);
 float CalcoloERCO(int z,int a,float egas,float esi, float ecsi,int bit, int codice);
 int Soglie_Fazietto(float zpart, float apart, float thetapart, float epart, int codice,float *zout, float *aout, float *eout,int *aid);
 void TabellaPerditeRCO(char *fileperdite);
 float Range(float zpart,float apart,float epart,int mate);
 float DE2E(float zpart,float apart,float elost,int mate,float spess);
float Epunch(float zpart,float apart,float spess,int mate);
 void InizializzaRCO();

 double Luce2E_csi(double z, double a, double lo);
 double E2Luce_csi(double z, double a, double E);

 protected:
  int fHaveTable;
  int fHaveTableRCO;
 
private:
	int calcolaGARF(int z,int a,int icsi,double *ein,double *e1,double *e2,double *e3,double *loo);
	int calcolaRCO(int, int, int, int, double*, double*, double*, double*, double*, int*);

};

#endif
