#include <iostream>
#include <stdio.h>
#include <string.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TObject.h>
#include <TCollection.h>
#include "Classe_analisi.h"
#include "Classe_formule.h"
#include <signal.h>
#include "Classe_geo.h"
#include <TNamed.h>


void mysignal(int a)
{

printf("segnale%d control c Scrivo il file di uscita e chiudo\n",a);

 Classe_analisi::Getanalisi()->RoutineFinale();
  Classe_analisi::Getanalisi()->ScriviFileUscita();

  exit(1);


  //  return;

}

int main(int argv,char *argc[])
{

  signal(2,mysignal);

  FILE *lista_comandi=0;
  if(argv>=2)
    {
      lista_comandi=fopen(argc[1],"r");
      
    }
  if(lista_comandi==0)
    {
      printf("Manca il file dei comandi\n");
      return 1;
    }

  Classe_analisi::Getanalisi()->Leggi_Wapstra();
  TNamed *str[1000];
  int nstr=0;
 
  printf("Prendo i comandi dal file %s\n",argc[1]);
 
  str[nstr]=new TNamed(Form("file comandi=%s",argc[1]),Form("s=%d",nstr));;
  nstr++;
 char tipo[1000];
 char filename[1000],file_histo[1000];
 char file_uscita[1000];
 char file_geo[1000];
 char ntu[100];
 
 char geo[10];
 int iexp=-1;
 int nev=0;
 TFile *f[200];

 int eventi_richiesti[200];
 int evento_iniziale[200];
 int eventi_totali[200];

 TChain *chain=0;
 int zp,ap,zt,at;
 int zpv,apv,ztv,atv;
 int trigplast;
 float ebeam,spess;
 float ebeamv,spessv;
 zpv=-1;
 apv=-1;
 atv=-1;
 ztv=-1;
 ebeamv=-1;
 spessv=-1;
 for(int j=0;j<200;j++)
   {
     eventi_richiesti[j]=0;
     evento_iniziale[j]=0;
     eventi_totali[j]=0;
     f[j]=0;
   }
 int nfile=0;
 char riga[1000];
 int iniz=0;
 while(fscanf(lista_comandi,"\n%[^\n]",riga)!=EOF)
   {

     if(riga[0]!='#')
       {

     sscanf(riga,"%s",tipo);
    
     str[nstr]=new TNamed(riga,Form("s=%d",nstr));
     // printf("%s %d\n",riga,nstr);
     nstr++;

     //     if(strcmp(tipo,"MC")==0 || strcmp(tipo,"EXPCOMP")==0||strcmp(tipo,"EXP-READ")==0)
     if(strcmp(tipo,"MC")==0 || strcmp(tipo,"EXPCOMP")==0||strcmp(tipo,"EXPBO")==0||strcmp(tipo,"EXP-READ")==0 || strcmp(tipo,"MC_HIPSE")==0|| strcmp(tipo,"MC_BAIOCCO")==0||strcmp(tipo,"MC_AMD")==0||strcmp(tipo,"MC_G++_TWINGO")==0||strcmp(tipo,"MC_G++_AMD")==0||strcmp(tipo,"MC_G++_BLOB")==0||strcmp(tipo,"MC_f90_TWINGO")==0||strcmp(tipo,"MC_f90_AMD")==0||strcmp(tipo,"LANGE")==0||strcmp(tipo,"EXP-FAZIA")==0||strcmp(tipo,"MC_SIMON_AMD")==0||strcmp(tipo,"MC_HFL_AMD")==0||strcmp(tipo,"MC_BAIOCCO_PAR")==0||strcmp(tipo,"COMD")==0||strcmp(tipo,"MC_HFL_HIPSE")==0||strcmp(tipo,"EXP-INDRAFAZIA")==0||strcmp(tipo,"MC_G++_AMD_KALISIM")==0)
       {

	 if(strcmp(tipo,"MC")==0 || strcmp(tipo,"MC_HIPSE")==0||strcmp(tipo,"MC_BAIOCCO")==0||strcmp(tipo,"MC_AMD")==0||strcmp(tipo,"MC_G++_TWINGO")==0||strcmp(tipo,"MC_G++_AMD")==0||strcmp(tipo,"MC_G++_BLOB")==0||strcmp(tipo,"MC_f90_TWINGO")==0||strcmp(tipo,"MC_f90_AMD")==0||strcmp(tipo,"LANGE")==0||strcmp(tipo,"MC_SIMON_AMD")==0||strcmp(tipo,"MC_HFL_AMD")==0||strcmp(tipo,"MC_BAIOCCO_PAR")==0||strcmp(tipo,"COMD")==0||strcmp(tipo,"MC_HFL_HIPSE")==0||strcmp(tipo,"MC_G++_AMD_KALISIM")==0)
       {
	 printf("Analisi Simulazione\n");

	 sscanf(riga,"%s %s %s %d",tipo,filename,geo,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %s %d %d",tipo,filename,geo,&nev,&iniz);
	   }
	 if(strcmp(tipo,"MC")==0)
	   {
	     printf("Gemini\n");
	     iexp=0;
	   }
	 if(strcmp(tipo,"MC_HIPSE")==0)
	   {
	     printf("HIPSE\n");
	     iexp=1;
	   }
	 if(strcmp(tipo,"MC_BAIOCCO")==0)
	   {
	     printf("BAIOCCO\n");
	     iexp=2;
	   }
	 if(strcmp(tipo,"MC_AMD")==0)
	   {
	     printf("AMD con sdecay\n");
	     iexp=3;
	   }
	 if(strcmp(tipo,"MC_G++_TWINGO")==0)
	   {
	     printf("Twingo con Gemini++\n");
	     iexp=4;

	   }
	 if(strcmp(tipo,"MC_G++_AMD")==0)
	   {
	     printf("AMD con Gemini++\n");
	     iexp=5;

	   }
	 if(strcmp(tipo,"MC_f90_TWINGO")==0)
	   {
	     printf("TWINGO con Geminif90\n");
	     iexp=6;
	   }
	 if(strcmp(tipo,"MC_f90_AMD")==0)
	   {
	     printf("AMD con geminif90\n");
	     iexp=7;
	   }
	 if(strcmp(tipo,"LANGE")==0)
	   {
	     printf("Langevin4D\n");
	     iexp=8;
	   }
	 if(strcmp(tipo,"MC_G++_BLOB")==0)
	   {
	     printf("BLOB con Gemini++\n");
	     iexp=9;

	   }
	 if(strcmp(tipo,"MC_SIMON_AMD")==0)
	   {
	     printf("AMD con Simon\n");
	     iexp=10;

	   }
	 if(strcmp(tipo,"MC_HFL_AMD")==0)
	   {
	     printf("AMD con HFL\n");
	     iexp=11;

	   }
	 if(strcmp(tipo,"MC_BAIOCCO_PAR")==0)
	   {
	     printf("HFL parallelo\n");
	     iexp=12;

	   }
	 if(strcmp(tipo,"COMD")==0)
	   {
	     printf("COMD\n");
	     iexp=13;

	   }
	 if(strcmp(tipo,"MC_HFL_HIPSE")==0)
	   {
	     printf("HIPSE con HFL\n");
	     iexp=14;

	   }

	 if(strcmp(tipo,"MC_G++_AMD_KALISIM")==0)
	   {
	     printf("AMD+GEMINI++ filtrato con KALIVEDASIM per INDRA FAZIA\n");
	     iexp=15;
	   }	 
	 if(strcmp(geo,"geo")==0)
	   {
	     iexp=iexp+100;
	     printf("Analisi in geometria\n");
	   }
	 else
	   {
	     printf("Analisi in 4pi\n");
	   }
	 //<<<<<<< faianalisi.C

	 str[nstr]=new TNamed(Form("%s %s",tipo,geo),Form("%s %d",filename,nev));
 nstr++;

       }
     if(strcmp(tipo,"EXPCOMP")==0)
       {
	 printf("Analisi Dati Sperimentali Forma Compatta\n");
	 iexp=210;
	 sscanf(riga,"%s %s %d",tipo,filename,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %d %d",tipo,filename,&nev,&iniz);
	   }
	 str[nstr]=new TNamed(tipo,Form("%s %d",filename,nev));
 nstr++;
       }
          if(strcmp(tipo,"EXPBO")==0)
       {
	 printf("Analisi Dati Sperimentali Ntupla di Bologna\n");
	 iexp=220;
	 sscanf(riga,"%s %s %d",tipo,filename,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %d %d",tipo,filename,&nev,&iniz);
	   }
	 str[nstr]=new TNamed(tipo,Form("%s %d",filename,nev));
 nstr++;
       }

          if(strcmp(tipo,"EXP-FAZIA")==0)
       {
	 printf("Analisi Dati Sperimentali Fazietto (uscita di Kaliveda)\n");
	 iexp=230;
	 sscanf(riga,"%s %s %d",tipo,filename,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %d %d",tipo,filename,&nev,&iniz);
	   }
	 str[nstr]=new TNamed(tipo,Form("%s %d",filename,nev));
 nstr++;
       }



          if(strcmp(tipo,"EXP-INDRAFAZIA")==0)
       {
	 printf("Analisi Dati Sperimentali INDRA FAZIA (uscita di Kaliveda)\n");
	 iexp=240;
	 sscanf(riga,"%s %s %d",tipo,filename,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %d %d",tipo,filename,&nev,&iniz);
	   }
	 str[nstr]=new TNamed(tipo,Form("%s %d",filename,nev));
 nstr++;
       }


     if(strcmp(tipo,"EXP-READ")==0)
       {
	 printf("Analisi Dati Sperimentali Forma Leggibile\n");
	 iexp=200;
	 sscanf(riga,"%s %s %d",tipo,filename,&nev);
	 if(nev!=-1)
	   {
	     sscanf(riga,"%s %s %d",tipo,filename,&nev);
	   }
	 //<<<<<<< faianalisi.C
	 str[nstr]=new TNamed(tipo,Form("%s %d %d",filename,nev,iniz));
 nstr++;

       }



     f[nfile]=new TFile(filename);

     if(!f[nfile]->IsZombie())
       {
     printf("Apro il file %s\n",f[nfile]->GetName());
     eventi_richiesti[nfile]=nev;
     evento_iniziale[nfile]=iniz;
     nfile++;
       }
       }

     if(strcmp(tipo,"GEO")==0)
       {
	 sscanf(riga,"%s %s",tipo,file_geo);
	 printf("File geometria=%s\n",file_geo);
       }

     if(strcmp(tipo,"HISTO")==0)
       {
	 sscanf(riga,"%s %s",tipo,file_histo);
       }
     if(strcmp(tipo,"REAC")==0)
       {
	 sscanf(riga,"%s %d %d %d %d %f %f %d",tipo,&zp,&ap,&zt,&at,&ebeam,&spess,&trigplast);
       }
     if(strcmp(tipo,"REAC_VERA")==0)
       {
	 sscanf(riga,"%s %d %d %d %d %f %f",tipo,&zpv,&apv,&ztv,&atv,&ebeamv,&spessv);
       }

     if(strcmp(tipo,"USCITA")==0)
       {
	 sscanf(riga,"%s %s",tipo,file_uscita);
	 Classe_analisi::Getanalisi()->ApriFileUscita(file_uscita);
       }

     if(strcmp(tipo,"NTUPLA")==0)
       {
	 sscanf(riga,"%s %s",tipo,ntu);
	 if(strcmp(ntu,"yes")==0)
	   {
	     cout<<"inserisco nel file di uscita anche ntupla"<<endl;
	     Classe_analisi::Getanalisi()->CreaNtupla();
	   }
       }

       }
   }
 fclose(lista_comandi);
//Non funziona; dovrebbe permettere di trovare automaticamente il nome del TTree
// if(nfile>0)
//    {
//      TIter next(f[0]->GetListOfKeys());
//      TObject *obj;
//      while(obj=(TObject*)next())
//        {
// 	 printf("NOME=%s %d\n",obj->GetName(),obj->InheritsFrom("TTree"));
// 	 if(obj->InheritsFrom("TTree"))
// 	   {
	    
// 	     printf("%s\n",Form("%s",obj->GetName()));
// 	   }
//        }
//    }

 if(iexp==0||iexp==100)
   {
     chain=new TChain("geminitree","geminitree");
   }
if(iexp==4||iexp==104)
  {
    chain=new TChain("twingo_geminiC","twingo_geminiC");
  }
 if(iexp==1||iexp==101)
   {
     chain=new TChain("hipse","hipse");
   }
 if(iexp==2||iexp==102)
   {
     chain=new TChain("baiocco","baiocco");
   }
 if(iexp==3||iexp==103)
   {
     chain=new TChain("amd_sec","amd_sec");
   }
 if(iexp==5||iexp==105 || iexp==15 || iexp==115)
   {
     chain=new TChain("amd_geminiC","amd_geminiC");

     if(iexp==115)
       {
	 iexp=15;
	 printf("ATTENZIONE: dati gia' filtrati da KALIVEDASIM, sono solo in geometria (ma si devono analizzare come se fossero 4pi)\n");
       }

   }
 if(iexp==6||iexp==106)
   {
     chain=new TChain("twingo_geminif90","twingo_geminif90");
   }
 if(iexp==7||iexp==107)
   {
     chain=new TChain("amd_geminif90","amd_geminif90");
   }
if(iexp==8||iexp==108)
   {
     chain=new TChain("langevin","langevin");
   }
if(iexp==9||iexp==109)
   {
     chain=new TChain("blob_geminiC","blob_geminiC");
   }

if(iexp==10||iexp==110)
   {
     chain=new TChain("amd_simon","amd_simon");
   }
if(iexp==11||iexp==111)
   {
     chain=new TChain("amd_sec_hfl","amd_sec_hfl");
   }
if(iexp==12||iexp==112)
   {
     chain=new TChain("baiocco","baiocco");
   }

if(iexp==13||iexp==113)
   {
     chain=new TChain("comd","comd");
   }
if(iexp==14||iexp==114)
   {
     chain=new TChain("hipse_sec_hfl","hipse_sec_hfl");
   }

  if(iexp==200||iexp==210)
   {
     chain=new TChain("odietree","odietree");
   }

 if(iexp==220)
   {
     chain=new TChain("pi","pi");
   }
 if(iexp==230||iexp==240)
   {
     chain=new TChain("events","events");
   }
     for(int j=0;j<nfile;j++)
       {
	 chain->Add(f[j]->GetName());
	 printf("Aggiungo il file %s\n",f[j]->GetName());       
	 if(j==0){eventi_totali[j]=chain->GetEntries();}
	 else{eventi_totali[j]=chain->GetEntries()-eventi_totali[j-1];}
	 if(eventi_richiesti[j]==-1)
	   {
	     eventi_richiesti[j]=eventi_totali[j];
	   }
       }

   
 printf("%lld\n",chain->GetEntries());
 str[nstr]=new TNamed(Form("Entries totali=%lld",chain->GetEntries()),Form("s=%d",nstr));
 nstr++;

 if(Classe_analisi::Getanalisi()->Getfout()==0)
   {
     printf("non c'e' file di uscita definito; si apre uscita_default.root\n");
     sprintf(file_uscita,"uscita_default.root");
     Classe_analisi::Getanalisi()->ApriFileUscita(file_uscita);
   }

 Classe_analisi::Getanalisi()->SetChain(chain);
 Classe_analisi::Getanalisi()->SetTipoAnalisi(iexp);

Classe_analisi::Getanalisi()->Set_Reazione(zp,ap,zt,at,ebeam,spess);
 if(zpv>0)
   {
Classe_analisi::Getanalisi()->Set_ReazioneVera(zpv,apv,ztv,atv,ebeamv,spessv);
   }
 else
   {
Classe_analisi::Getanalisi()->Set_ReazioneVera(zp,ap,zt,at,ebeam,spess);
   }
 Classe_analisi::Getanalisi()->Set_TrigPlastichino(trigplast);
Classe_geo::Getgeo()->Leggigeo(file_geo); 





 Classe_analisi::Getanalisi()->Set_Eventi(eventi_totali,eventi_richiesti,evento_iniziale,nfile);

Classe_analisi::Getanalisi()->DefinisciHisto(file_histo);


 for(int j=0;j<nstr;j++)
   {
     Classe_analisi::Getanalisi()->Getfout()->cd();
    
     str[j]->Write();
   }
if(Classe_geo::Getgeo()->ncuts>0)
	{
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	  Classe_analisi::Getanalisi()->gggcuts[k]->Write();
	    }
	}
Classe_analisi::Getanalisi()->Loop_eventi();


 
Classe_analisi::Getanalisi()->RoutineFinale();
Classe_analisi::Getanalisi()->ScriviFileUscita();
  return 0;
}

