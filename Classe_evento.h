#ifndef EVENTO
#define EVENTO

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>

#include <TROOT.h>
#include <stdlib.h>
#include <vector>
using namespace std;


class Classe_evento
{
 public:
	int tempo;
  struct eventomc {
    //variabili originali lette dal montecarlo
    UInt_t moltepl;
    float z[500];
    float a[500];
    float phicm[500]; // il phi del montecarlo va da -180 a 180 come deve
    float thetacm[500];
    float epartcm[500];
    UInt_t moltgamma;
    float egamma[500];
        UInt_t fissione;
     float spin;
    int moltcharged;
    int moltneutr;
    UInt_t CN;
    UInt_t isresidue;
    float zprimari[500];
    float aprimari[500];
    float vcm_primari[500];
    float thecm_primari[500];
    float phi_primari[500];
    float estar_primari[500];
    float spin_primari[500];
    float tkel;
    // per hipse
    //    int mch_mneutr;
    UInt_t mch_mneutr;
    float par_urto;
    float estar;
    float Px[500];
    float Py[500];
    float Pz[500];
    float Pxprim[500];
    float Pyprim[500];
    float Pzprim[500];
    float jx[500];
    float jy[500];
    float jz[500];

    UInt_t origine[500];
    //fine per hipse
    //per baiocco
    UInt_t discreto[500];
    float exc[500];
    //fine per baiocco
    //per AMD sdecay
    UInt_t nevorig;
    UInt_t zamd[500];
    UInt_t aamd[500];
    float zsecb[500];
    float asecb[500];
    float vxsecamd[500];
    float vysecamd[500];
    float vzsecamd[500];
    UInt_t moltprimari;
    UInt_t zamdprim[500];
    UInt_t aamdprim[500];
        float vxprimamd[500];
    float vyprimamd[500];
    float vzprimamd[500];
    //fine per AMD sdecay
    //per twingo gemini++
    UInt_t isresidue_prim[500];
       float pxprim[500];
    float pyprim[500];
    float pzprim[500];
    //fine per twingo gemini++

    //per langevin
    float eservizio[500];
    //fine per langevin
    
    // variabili calcolate (da valori originali)
    float vpartcm[500];
    float vpartcm_cart[500][3];
    float vpartlab_cart[500][3];
    float vpartlab[500];
    float thetalab[500];
    float philab[500];
    float epartlab[500];
    // geometria
    int index_phos[500];//box*10+phos (partono da 1)
    int index_garf[500];//isec*10+cesio (isec=1-24,cesio=5-8)
    int index_baf[500];//nbaf (1-8)
    //<<<<<<< Classe_evento.h
    int index_rco[500];
    int index_fazietto[500];
    int soglia[500];

    float estrip[500][4];//4 strip per ogni settore UL UR DL DR; il settore si ricava da index_garf; Left o Right da philab (0-7.5 Left; 7.5-15 Right) e ce lo dice il codice ritornato da Classe_geo::InsideGarf    
    float ecsi[500];
    float tvolo[500];
    //RCO
    float rco_gas[500];
    float rco_si[500];
    float rco_csi[500];
    float qf[500];
    // variabili ricalcolate dopo effetti sperimentali (uguali alle variabili originali se non c'e' geometria)
    float zexp[500];
    float aexp[500];
    float epartlabexp[500];
    float thetalabexp[500];
    float philabexp[500];
    float vpartlabexp[500];
    float vpartlab_cartexp[500][3];
    float vpartcm_cartexp[500][3];
    float vpartcmexp[500];
    float epartcmexp[500];
    float thetacmexp[500];
    float phicmexp[500];
    float tvoloexp[500];
    int codphos[500];
    int partbuona[500];

    int array[500]; //per output kvsim
    int ntele[500];
    int idcode[500];
    int ecode[500];
    int Ameasured[500];
  


  }mcevent;

struct EXP
{
 Double_t phos_tof[6][9];
 Double_t  phos_traw[6][9];
  Double_t phos_z[6][9];
  Double_t phos_a[6][9];
  Double_t phos_qf[6][9];//100 se identif. da contorno; diverso se id da ga gb con linee Z
  Double_t phos_cod[6][9];//1=id da gatof; 2=id da gagb; 3=id da gbgc

  Double_t phos_ga[6][9];//ga utile per avere un'idea dell'energia
  Double_t garf_z[24][8];
  Double_t garf_a[24][8];
  Double_t garf_qf[24][8];
  Double_t garf_theta[24][8];
  Double_t garf_phi[24][8];
  Double_t garf_epart[24][8];
  Double_t garf_tvolo[24][8];
  
  Double_t garf_raw_egas[24][8];
 Double_t garf_raw_egashg[24][8];
 
  Double_t garf_raw_fast[24][8];
  Double_t garf_raw_lo[24][8];
  Double_t garf_raw_slowpsa[24][8];
  Double_t garf_cal_egas[24][8];  
  Double_t garf_cal_ecsi[24][8];  


  Double_t rco_z[8][9][7];//8 settori, 9 strip, 7 cesi; la nona strip e l'ottavo cesio vengono messi quando non si possono definire (es. roba stoppata nel silicio-> niente csi)
  Double_t rco_a[8][9][7];
  Double_t rco_epart[8][9][7];
  Double_t rco_code[8][9][7];
  Double_t rco_qf[8][9][7];
  Double_t rco_tvolo[8][9][7];
  Double_t rco_tcode[8][9][7]; //tempo dato da camera (1), Si (2), CsI (3)

  Double_t rco_raw_egas[8][9][7];
  Double_t rco_raw_esi[8][9][7];
  Double_t rco_raw_trise[8][9][7];
  Double_t rco_raw_fast[8][9][7];
  Double_t rco_raw_slow[8][9][7];
  Double_t rco_raw_slowpsa[8][9][7];
  Double_t rco_cal_egas[8][9][7];
  Double_t rco_cal_esi[8][9][7];
  Double_t rco_cal_ecsi[8][9][7];

  int trig;
  float tplast;
Long64_t run;
} expevent;

struct EXP_COMP
{
  UInt_t value_N;//n. di valori nell'evento: 5*n. di phos buoni
  Float_t value_val[1000];//per ora 5 valori per ogni phos dell'evento
  UInt_t value_worker_id[1000];//n. phos 1000*ip+100*ih
  UInt_t value_worker_class_code[1000];//15328 per il worker phosid
  char run_num[200];
  UInt_t event_num; 
  Long64_t run;

}expcomp;

 struct EXP_BO
{
  Int_t molt;
  Float_t codiceriv[1000];
  Float_t vx[1000];
  Float_t vy[1000];
  Float_t vz[1000];
  Int_t z[1000];
  Float_t a[1000];
  Int_t qf[1000];
  Float_t e[1000];
  Float_t theta[1000];
  Float_t phi[1000];
  Int_t trigger[8];
}expbo;

struct FAZIAKALI
{
  int run;
  int mtot;
  int idtel[1000];//nome rivelatore 123 blocco 1 qua 2 tel 3 (per blocco 0, solo 2 numeri)
  int blocco[1000];
  int qua[1000];
  int tel[1000];
  int z[1000];
  int a[1000];
  int aid[1000]; //identificazione in massa; 1 massa misurata, 0 massa non misurata (da betastability)
  int idcode[1000];// qualita' dell'identificazione; 0 tutto ok, 4 male (1,2,3 cosi' e cosi'); aid=1 idcode=0 masse e zeta giuste; aid=0 idcode=0 zeta giuste
  int ecode[1000];// se 0, energia OK (calibrata per Si1,Si2,CsI); se 1 una delle energie non e' calibrata indipendentemente ma e' ricostruita; -1 non c'e' l'energia , 2-3 ci sono incoerenze nella calibrazione
  int idtype[1000];// come e' stata identificata la particella: 11 PSA in Si1; 12 DE-E Si1-Si2; 22 PSA Si2 (solo se 12 e' impossibile, es. manca Si1); 23 De-E Si2-CsI; 33 fast slow in CsI (solo se 23 non e' possibile, es particella non vista in Si o Si2 rotto)
  float esi1[1000]; //energia in Si1
  float esi2[1000];//energia in Si2
  float ecsi[1000];//energia in CsI (slow)
  float etot[1000];// energia totale della particella (dopo perdite in target, penso)
  float chsi1[1000];//energia in Si1 in canali
  float chsi2[1000];//energia in Si2 in canali
  float chcsi[1000];//energia in csi in canali
  float theta[1000];
  float phi[1000];
}faziakali;

struct INDRAFAZIA
{
  int run;
  int mtot;
  int idtel[1000];//nome rivelatore 123 blocco 1 qua 2 tel 3 (per blocco 0, solo 2 numeri)//fazia
  int blocco[1000];//fazia
  int qua[1000];//fazia
  int tel[1000];//fazia
  int ring[1000];//indra
  int module[1000];//indra
  int z[1000];
  int a[1000];
  int aid[1000]; //identificazione in massa; 1 massa misurata, 0 massa non misurata (da betastability)
  int idcode[1000];// qualita' dell'identificazione; 0 tutto ok, 4 male (1,2,3 cosi' e cosi'); aid=1 idcode=0 masse e zeta giuste; aid=0 idcode=0 zeta giuste
  int ecode[1000];// se 0, energia OK (calibrata per Si1,Si2,CsI); se 1 una delle energie non e' calibrata indipendentemente ma e' ricostruita; -1 non c'e' l'energia , 2-3 ci sono incoerenze nella calibrazione
  int idtype[1000];// come e' stata identificata la particella: 11 PSA in Si1; 12 DE-E Si1-Si2; 22 PSA Si2 (solo se 12 e' impossibile, es. manca Si1); 23 De-E Si2-CsI; 33 fast slow in CsI (solo se 23 non e' possibile, es particella non vista in Si o Si2 rotto)
  float esi1[1000]; //energia in Si1
  float esi2[1000];//energia in Si2
  float ecsi[1000];//energia in CsI (slow)
  float etot[1000];// energia totale della particella (dopo perdite in target, penso)
  float chsi1[1000];//energia in Si1 in canali
  float chsi2[1000];//energia in Si2 in canali
  float chcsi[1000];//energia in csi in canali
  float chcsifast[1000];//csifast in canali
  float theta[1000];
  float phi[1000];
  int array[1000]; // =1 se colpisco fazia, =0 se colpisco indra 
  float desi1[1000];
  float desi2[1000];
  int gt_dt[1000];
}indrafazia;


struct EXP_HECTOR
{
  float ebaf[8];
  float tbaf[8];

}exphector;

struct MixPar{
  vector<vector<float> > Erel;
  vector<vector<float> > Emix;
  vector<vector<vector<float> > > Prel;
 vector<vector<vector<float> > > Pcm; //mixatore 29-01-21
 vector<vector<vector<float> > > Pmix;
  vector<vector<vector<float> > > Pcm_mix;//mixatore 29-01-21
  vector<vector<float> > Thetarel;//mixatore 29-01-21
  vector<vector<float> > Thetamix;//mixatore 29-01-21
  
 vector<vector<vector<float> > > Vcm; //mixatore 29-01-21
 vector<vector<vector<float> > > Vcm_mix;//mixatore 29-01-21

  void reset(int i) {
    
    Erel.at(i).clear();
    Emix.at(i).clear();
    Prel.at(i).clear();
    Pmix.at(i).clear();
    Pcm.at(i).clear();//mixatore 29-01-21
    Pcm_mix.at(i).clear();//mixatore 29-01-21
    Thetarel.at(i).clear();//mixatore 29-01-21
    Thetamix.at(i).clear();//mixatore 29-01-21

    Vcm.at(i).clear();//mixatore 29-01-21
    Vcm_mix.at(i).clear();//mixatore 29-01-21

  }
  void build()
  {
    Erel.resize(40);

    Emix.resize(40);

    Prel.resize(40);

    Pmix.resize(40);
    Pcm.resize(40);//mixatore 29-01-21
    Pcm_mix.resize(40);//mixatore 29-01-21
    Thetarel.resize(40);//mixatore 29-01-21
    Thetamix.resize(40);//mixatore 29-01-21

    Vcm.resize(40);//mixatore 29-01-21
    Vcm_mix.resize(40);//mixatore 29-01-21
  }
} mixing;


struct EVENTO
{
  int moltepl;
  Long64_t run;
  vector <int> isformix;//per mixatore 29-01-21
  vector <float> z;
  vector <float> a;
  vector <float> epartcm;
  vector <float> thetacm;
  vector <float> phicm;
  vector <float> vpartcm;
  vector <float> vpartcm_x ;
  vector <float> vpartcm_y;
  vector <float> vpartcm_z;
  vector <vector <float> > vpcm;
  vector <float> epartlab;
  vector <float> thetalab;
  vector <float> philab;
  vector <float> vpartlab;
  vector <float> vpartlab_x ;
  vector <float> vpartlab_y;
  vector <float> vpartlab_z;
  vector <vector <float> > vplab;
  vector <int> indice_originale;//solo per montecarlo 
  vector <float> tvolo;//solo per sperimentale
  int bunch;
  vector <float> tbunch;//tempo per trovare il bunch riportato
  vector <int> codphos;//solo per phoswich mc ed exp 1 gatof, 2 gagb, 3 gbgc
  vector <int> phosqf;//solo per sperimentale phoswich 100 per id da gate, altro per linee Z
  vector <int> coderiv; // per i phos box*10+phos (parte da 1); per garfield isec*10+cesio+1000; BaF2 is 10000+nbaf
  vector <float> phosga;// ga (solo per i phos exp)
  vector <int> garfqf;//solo per exp garfield
  int trig;//trigger exp
  vector <float> tplast;// tempo del plastichino (solo exp)
  vector <int> rcoqf;//solo per exp RCo
  vector <int> rcocode;

}evento;


 

  Classe_evento()
    {
      evento.moltepl=0; 
      mixing.build();
                  
    }
  void Leggievento();
  void AnalizzaEvento();
 
  void AzzeraEvento();
  void NoeffettiExp();
  void Copia_in_Evento(int tipo_analisi);
  void AnalisiPrincipale();
  void Riporto();
  void RiportoGarfRCo(); //30/9/2015
void Escludi_Multiple_Hits();
void Correggi_Tempo_Phos(int ip,int ih,int j);


// versione vecchia di mixatore (fino a 29-01-21)
/* vector<float> erel2part(vector<vector<float> > mixingprel, int a1, int a2); */
/* vector<float> emix(vector<vector<float> > mixingprel, int a1, int a2); */

/* vector<vector <float> > prel2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2,int a1,int a2); */

/* vector<vector<float> > pmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart,int a1,int a2); */
/* void mixatore(int index, int n1, int z1, int a1, int n2, int z2, int a2); */

vector<vector <float> > prel2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2);//mixatore 29-01-21
vector<vector <float> > pcm2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2);//mixatore 29-01-21
vector<vector <float> > Vcm2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2);//mixatore 29-01-21
vector<vector<float> > pmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart,vector<vector<float> > ppart2);//mixatore 29-01-21
vector<vector<float> > pcmmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart,vector<vector<float> > ppart2);//mixatore 29-01-21
vector<vector<float> > Vcmmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart,vector<vector<float> > ppart2);//mixatore 29-01-21
vector<float>  thetarelmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart,vector<vector<float> > ppart2);//mixatore 29-01-21
vector<float>  thetarel2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2);//mixatore 29-01-21
void mixatore(int index, int codice, int n1, int z1, int a1, int n2, int z2, int a2);//mixatore 29-01-21

};
#endif
