#ifndef ANALISI
#define ANALISI
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <TChain.h>
#include <stdlib.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include "Classe_evento.h"

using namespace std;

class Classe_analisi
{
 public:
  int tipo_analisi;//0 MC Gemini 4p, 1 MC Hipse 4p, MC Baiocco 4p, 10 MC Gemini geo, 11 MC Hipse geo, 12 MC Baiocco geo
  //20 Exp forma leggibile, 21 exp compatto
  //Convenzione: 0-9 simulazioni 4p; 10-19 le stesse simulazioni geo; 20 exp human readable; 21 exp compatto
  TChain *ch;
  int numero_files;
  int *lista_eventi;
    Long64_t eventi_letti;
  // int eventi_letti;
  // int nentry;
  // long long int nentry;
         Long64_t nentry;
	 Long64_t ghost_part;
struct sistema
{
  int zp,ap,zt,at;
  float ebeam;//energia del fascio per nucleone
  float vcm;
  float vplab;
  float spess_target;
  float Ecm;//Ecm
  float mi;//massa ridotta
  float vp_cm;
  float bgrazing;
  float thegrazingcm;
  float thegrazinglab;
  float betacm;
  float gammacm;

  // Veri parametri della reazione (serve per MC analizzato come reazione diversa da quella che e' davvero)
  int zpv,apv,ztv,atv;
  float ebeamv;//energia del fascio per nucleone
  float vcmv;
  float vplabv;
  float spess_targetv;
  float Ecmv;//Ecm
  float miv;//massa ridotta
  float vp_cmv;
  float bgrazingv;
  float thegrazingcmv;
  float thegrazinglabv;
  float betacmv;
  float gammacmv;

}reazione;

struct trplast
{
  int Trigplast;
}triggerplastichino;

 TH1D *h1d[10000];
 TH2F *h2d[10000];
 TGraph *gggcuts[100];
TGraphErrors *g[10000];
 TTree *ntupla;
 int n1d,n2d;
int ng;
 int ntupl;


struct Tree
{
  int moltepl;
  float z[200];
  float a[200];
 float vxcm[200];
  float vycm[200];
  float vzcm[200];
  float vlabmod[200];
  float thelab[200];
  float philab[200];

  //float vxlab[200];
  //float vylab[200];
  // float vzlab[200];
  //float vcmmod[200];
  // float thecm[200];
  //float phicm[200];
  

}tree;

float D[111][163]; //Difetti di massa Wapstra z,n
float Delta[111][274]; //Difetti di massa Wapstra z,a

 struct sistema *reac;

 static Classe_analisi *ana;

  Classe_evento *Evento;
  TFile *fout; 
 
//int buffer;


 vector<vector<vector<vector <float> > > > tableP1;

  vector<vector<vector<vector <float> > > > tableP2;

vector <int> iokmix;
vector <int> jokmix;
vector <int> ioktimesmix;



  Classe_analisi()
    {
      ch=0;
      tipo_analisi=-1;
      numero_files=1;
      lista_eventi=0;
      Evento=0;
      reac=&reazione;
      fout=0;
      n1d=0;
      n2d=0;
      ng=0;
      eventi_letti=0;
      ghost_part=0;
       for(int j=0;j<10000;j++)
      	{
      	  h1d[j]=0;
      	  h2d[j]=0;

      	}
       for(int j=0;j<10000;j++)
	 {
	   g[j]=0;
	 }

       for(int j=0;j<100;j++)
	 {
	   gggcuts[j]=0;
	 }
      //      printf("qui %s\n",ch->GetName());
       nentry=0;
       ntupl=0;       
       ntupla=0;

       tableP1.resize(40);
       for(int i=0;i<tableP1.size();i++)
	 {
	   tableP1.at(i).resize(100);
	   //  cout<<tableP1.at(i).size()<<endl;
	 }
       //       tableP1.at(tableP1.size()-1).resize(100);
       tableP2.resize(40);
       for(int i=0;i<tableP2.size();i++)
	 {
	   tableP2.at(i).resize(100);
	 }
       //       tableP2.at(tableP2.size()-1).resize(100);

       iokmix.resize(40);
       jokmix.resize(40);
       ioktimesmix.resize(40);

       //       buffer=100;


 for(int j=0;j<111;j++)
    {
      for(int k=0;k<163;k++)
	{
	  D[j][k]=-100000;//Difetti di massa di Wapstra z,n
	  

	}
     for(int k=0;k<274;k++)
	{
	  Delta[j][k]=-100000;//Difetti di massa di Wapstra z,a
	  

	}
    }



    }
  void Set_Reazione(int Zp,int Ap,int Zt,int At,float Ebeam,float spess);
  void Set_ReazioneVera(int Zp,int Ap,int Zt,int At,float Ebeam,float spess);
  void  Loop_eventi();
  void Set_Eventi(int *eventi_totali,int *eventi_richiesti,int *evento_iniziale,int nfile);
  void ApriFileUscita(char *file);
  void DefinisciHisto(char *file);
  void ScriviFileUscita();
  void CreaNtupla();
void Leggi_Wapstra();

 static Classe_analisi *Getanalisi();
 void SetChain(TChain *chain)
 {
   ch=chain;
 }
 void SetTipoAnalisi(int itip)
 {
   tipo_analisi=itip;
 }
 TFile *Getfout(){return fout;}
 void Set_TrigPlastichino(int trigplast);
void RoutineFinale();
};


#endif
