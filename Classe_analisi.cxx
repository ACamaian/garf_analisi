#include "Classe_analisi.h"
#include "Classe_formule.h"
#include <time.h>
#include <TROOT.h>

Classe_analisi *Classe_analisi::ana=0;
Classe_analisi *Classe_analisi::Getanalisi()
{
  if(ana!=0)
    {
      return ana;
    }
  ana=new Classe_analisi();
  return ana;
}

void Classe_analisi::Loop_eventi()
{
  //	time_t end;

  //	int h,m,s;

       const char up[4] = {0x1B,'[','A',0};
        double ev_al_sec;
        time_t start,now;
        double diff=-1;


  nentry=0;
  printf("n. tot di eventi=%lld\n",ch->GetEntries());
ghost_part=0;
  for(int j=0;j<numero_files;j++)
    {
     for(Long64_t k=lista_eventi[j];k<lista_eventi[j+numero_files];k++)
	      //for(int k=lista_eventi[j];k<lista_eventi[j+numero_files];k++)
	{

	  if(k<ch->GetEntries())
	    {
	      nentry=k;

              //                   if(eventi_letti%1000000==1){printf("eventi letti %lld\n",eventi_letti);}
                   if(eventi_letti%100000==1)
{
 //printf("eventi letti %lld\n",eventi_letti);
if(diff<0) {
                                  start=time(NULL); diff=0;
                                  printf("Inizio analisi: %s",ctime(&start));
                        }
                          else {
                                  now=time(NULL);
                                  ev_al_sec=(double)k/((double)now-(double)start);
                                  printf("Now event %7lld of %lld (%5.1f%%). %7.2f ev/sec. ETA: %3d'%2.2d\" TT: %d'\n%s",
eventi_letti-1, ch->GetEntries(), 100.*(eventi_letti)/(float)ch->GetEntries(),
                 ev_al_sec ,
((int)((ch->GetEntries()-eventi_letti-1)/ev_al_sec))/60,((int)((ch->GetEntries()-eventi_letti-1)/ev_al_sec))%60,(now-start)/60,
up);
                        }



 }
     
      Evento=new Classe_evento();

      Evento->AzzeraEvento();


      	   Evento->Leggievento();


	   Evento->AnalizzaEvento();
	   delete Evento;
	  Evento=0;
	  eventi_letti++;
	  //    printf("eventi letti %lld nentry %lld\n",eventi_letti,nentry);	      
	 
	    }
	}
    }



}
void Classe_analisi::Set_Eventi(int *eventi_totali,int *eventi_richiesti,int *evento_iniziale,int nfile)
{

  numero_files=nfile;
  lista_eventi=(int *)malloc(2*numero_files*sizeof(int));
  for(int j=0;j<nfile;j++)
    {
      int sumj=0;
      for(int k=0;k<j;k++)
	{
	  sumj=sumj+eventi_totali[k];
	}
      lista_eventi[j]=sumj+evento_iniziale[j];
      lista_eventi[j+nfile]=lista_eventi[j]+eventi_richiesti[j];
    }
  //  printf("nfile=%d\n",nfile);
  //for(int j=0;j<nfile;j++)
  // {      
      //  printf("j%d %d %d\n",j,lista_eventi[j],lista_eventi[j+nfile]);
  //  }

}

void Classe_analisi::Set_Reazione(int Zp,int Ap,int Zt,int At,float Ebeam,float spess)
{

  reazione.zp=Zp;
  reazione.ap=Ap;
  reazione.zt=Zt;
  reazione.at=At;
  reazione.ebeam=Ebeam;
  reazione.spess_target=spess;
  printf("ebeam tot=%f %f %d\n",reazione.ebeam*reazione.ap,reazione.ebeam,reazione.ap);
 reazione.vplab=Classe_formule::vplab_classica(reazione.ebeam);
 reazione.vcm=Classe_formule::vcm_classica(reazione.vplab,reazione.ap,reazione.at); 
 reazione.vp_cm=reazione.vplab-reazione.vcm;
 printf("vcm=%f mm/ns vp_lab=%f mm/ns vp_cm=%f mm/ns\n",reazione.vcm,reazione.vplab,reazione.vp_cm);
 reazione.mi=(float)Ap*(float)At/(float)(Ap+At);
 reazione.Ecm=reazione.mi*reazione.ebeam;
 cout<<"Massa ridotta="<<reazione.mi<<" Ebeam="<<reazione.ebeam<<" Ecm="<<reazione.Ecm<<endl;
 float r0=1.7445-0.04388*log(reazione.zp*reazione.zt);
  float R12=r0*(pow(reazione.ap,0.3333)+pow(reazione.at,0.3333));
  float Vcoul=1.44*reazione.zp*reazione.zt/R12;
   reazione.bgrazing=R12*sqrt(1-Vcoul/reazione.Ecm);
   reazione.thegrazingcm=57.296*2*atan(1.44*reazione.zt*reazione.zp/(2*reazione.Ecm*reazione.bgrazing));
 
  reazione.thegrazinglab=57.296*atan(reazione.vp_cm*sin(reazione.thegrazingcm/57.296)/(reazione.vcm+reazione.vp_cm*cos(reazione.thegrazingcm/57.296)));
   printf("bgrazing (fm)=%f thetagrazing_cm=%f thetagrazing_lab=%f\n",reazione.bgrazing,reazione.thegrazingcm,reazione.thegrazinglab);
   reazione.betacm=reazione.vcm/Classe_formule::cluce;
   reazione.gammacm=1/sqrt(1-pow(reazione.betacm,2));
}


//Serve quando si analizza il MC come se fosse una reazione diversa da quella che e'
void Classe_analisi::Set_ReazioneVera(int Zp,int Ap,int Zt,int At,float Ebeam,float spess)
{

  reazione.zpv=Zp;
  reazione.apv=Ap;
  reazione.ztv=Zt;
  reazione.atv=At;
  reazione.ebeamv=Ebeam;
  reazione.spess_targetv=spess;
  
 reazione.vplabv=Classe_formule::vplab_classica(reazione.ebeamv);
 reazione.vcmv=Classe_formule::vcm_classica(reazione.vplabv,reazione.apv,reazione.atv); 
 reazione.vp_cmv=reazione.vplabv-reazione.vcmv;

 reazione.miv=(float)Ap*(float)At/(float)(Ap+At);
 reazione.Ecmv=reazione.miv*reazione.ebeamv;

 float r0=1.7445-0.04388*log(reazione.zpv*reazione.ztv);
  float R12=r0*(pow(reazione.apv,0.3333)+pow(reazione.atv,0.3333));
  float Vcoul=1.44*reazione.zpv*reazione.ztv/R12;
   reazione.bgrazingv=R12*sqrt(1-Vcoul/reazione.Ecmv);
   reazione.thegrazingcmv=57.296*2*atan(1.44*reazione.ztv*reazione.zpv/(2*reazione.Ecmv*reazione.bgrazingv));
 
  reazione.thegrazinglabv=57.296*atan(reazione.vp_cmv*sin(reazione.thegrazingcmv/57.296)/(reazione.vcmv+reazione.vp_cmv*cos(reazione.thegrazingcmv/57.296)));

   reazione.betacmv=reazione.vcmv/Classe_formule::cluce;
   reazione.gammacmv=1/sqrt(1-pow(reazione.betacmv,2));
   if(tipo_analisi<200 && (TMath::Abs(reazione.vcmv-reazione.vcm)>0.01||reazione.zpv!=reazione.zp || reazione.ztv!=reazione.zt ||reazione.apv!=reazione.ap||reazione.atv!=reazione.at||TMath::Abs(reazione.ebeam-reazione.ebeamv)>0.1))
     {
       cout<<"Per il MC quando si analizza un sistema con una cinematica diversa da quella vera (per fare il fondo dovuto a un contaminante):"<<endl;
 cout<<"Vera reazione: Massa ridotta="<<reazione.miv<<" Ebeam="<<reazione.ebeamv<<" Ecm="<<reazione.Ecmv<<endl;
printf("Vera reazione: ebeam tot=%f %f %d\n",reazione.ebeamv*reazione.apv,reazione.ebeamv,reazione.apv);
 printf("Vera reazione: vcm=%f mm/ns vp_lab=%f mm/ns vp_cm=%f mm/ns\n",reazione.vcmv,reazione.vplabv,reazione.vp_cmv);
   printf("Vera reazione: bgrazing (fm)=%f thetagrazing_cm=%f thetagrazing_lab=%f\n",reazione.bgrazingv,reazione.thegrazingcmv,reazione.thegrazinglabv);
     }
}


void Classe_analisi::Set_TrigPlastichino(int trigplast)
{
  triggerplastichino.Trigplast=trigplast;
}

void Classe_analisi::ApriFileUscita(char *file)
{
  fout=new TFile(file,"RECREATE");
  printf("File di uscita %s\n",fout->GetName());
}

void Classe_analisi::DefinisciHisto(char *file_uscita)
{
  FILE *apri=fopen(file_uscita,"r");
  if(apri==0)
    {
      printf("Manca lista istogrammi\n");
      return;
    }
  char num[100];
  int chx,chy;
  float x1,x2,y1,y2;
  char riga[1000];

  while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
    {
      if(riga[0]!='#')
	{
	  sscanf(riga,"%s",num);
	  if(strlen(num)==strlen(riga))
	    {
	      
	      g[ng]=new TGraphErrors();
	      g[ng]->SetName(num);
	      g[ng]->SetMarkerStyle(20);
	      gROOT->GetListOfSpecials()->Add(g[ng]);
	      ng++;
	      
	    }
	  else
	    {
	      
	  sscanf(riga,"%s %d %f %f %d %f %f",num,&chx,&x1,&x2,&chy,&y1,&y2);
	  if(chy>1)
	    {
	      h2d[n2d]=new TH2F(Form("h%s",num),Form("h%s",num),chx,x1,x2,chy,y1,y2);
	      h2d[n2d]->SetContour(128);
	      gROOT->GetListOfSpecials()->Add(h2d[n2d]);
	      n2d++;
	    }
	  else
	    {
	      h1d[n1d]=new TH1D(Form("h%s",num),Form("h%s",num),chx,x1,x2);
	      gROOT->GetListOfSpecials()->Add(h1d[n1d]);
	      n1d++;
	    }
	    
	    }
	}
    }
  fclose(apri);

}
void Classe_analisi::ScriviFileUscita()
{
     printf("Eventi analizzati=%lld\n",eventi_letti);
     // printf("Eventi analizzati=%d\n",eventi_letti);

  fout->cd();

  for(int j=0;j<n1d;j++)
    {
      h1d[j]->Write();
    }
  for(int j=0;j<n2d;j++)
    {
      h2d[j]->Write();
    }
  if(ntupl==1)
    {
  ntupla->Write();
    }
  for(int j=0;j<ng;j++)
    {
      g[j]->Write();
    }

  fout->Close();

  return;
}


void Classe_analisi::CreaNtupla()
{
  ntupl=1;
  ntupla=new TTree();
  ntupla->SetName("ntupla");
  ntupla->SetTitle("ntupla");
      
	  ntupla->Branch("evento",&nentry,"evento/L");
	    ntupla->Branch("moltepl",&tree.moltepl,"moltepl/I");
	    ntupla->Branch("z",tree.z,"z[moltepl]/F");
	    ntupla->Branch("a",tree.a,"a[moltepl]/F");
	    ntupla->Branch("vxcm",tree.vxcm,"vxcm[moltepl]/F");
	    ntupla->Branch("vycm",tree.vycm,"vycm[moltepl]/F");
	    ntupla->Branch("vzcm",tree.vzcm,"vzcm[moltepl]/F");

	    //	    ntupla->Branch("vxlab",tree.vxlab,"vxlab[moltepl]/F");
	    //ntupla->Branch("vylab",tree.vylab,"vylab[moltepl]/F");
	    // ntupla->Branch("vzlab",tree.vzlab,"vzlab[moltepl]/F");
	    
	    ntupla->Branch("vlabmod",tree.vlabmod,"vlabmod[moltepl]/F");
	    ntupla->Branch("thelab",tree.thelab,"thelab[moltepl]/F");
	    ntupla->Branch("philab",tree.philab,"philab[moltepl]/F");
	    //ntupla->Branch("vcmmod",tree.vcmmod,"vcmmod[moltepl]/F");
	    // ntupla->Branch("thecm",tree.thecm,"thecm[moltepl]/F");
	    // ntupla->Branch("phicm",tree.phicm,"phicm[moltepl]/F");



}
void Classe_analisi::Leggi_Wapstra()
{
  FILE *wapstra=fopen("awm95.txt","r");
  if(wapstra==0)
    {
      cout<<"Manca awm95.txt"<<endl;
      exit(0);
    }
  int zl,al;
  double ddelta;
  char iso[10];

  char riga[1000];
 
  cout<<"Sto leggendo awm95.txt"<<endl;
  while(fscanf(wapstra,"\n%[^\n]",riga)!=EOF)
    {
      if(riga[0]!='#')
	{
	  sscanf(riga,"%s %d %lg",iso,&zl,&ddelta);
	  sscanf(iso,"%d",&al);
	  //printf("%s %d %f %d\n",iso,z,ddelta,a);
	  int in=al-zl;
	  D[zl][in]=ddelta*0.001;
	  Delta[zl][in+zl]=D[zl][in];
	  // cout<<zl<<" "<<in+zl<<" "<<Delta[zl][in+zl]<<endl;
	  
	}
    }

  fclose(wapstra);
  double me=0.511;
  double amu=931.49;
    double mp=D[1][0]+amu-me;
  double mn=D[0][1]+amu;


  return;

}
