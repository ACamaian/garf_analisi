#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <TH2.h>
#include <TControlBar.h>
#include <TFile.h>
#include <TList.h>
#include <TLatex.h>
#include <TMath.h>
#include <TObject.h>
#include <TCollection.h>
#include <unistd.h>
//#include "Mygraf.h"
#ifndef ROOT_TNamed
#include <TGraph.h>
#include <TPad.h>
#endif
#include "Mybar.h"
#include <TF1.h>
#include <TCutG.h>
#include "Mydialogo.h"
#include <TStyle.h>
#include <math.h>
#include <TROOT.h>
#include <TTree.h>
#include <TStyle.h>
#include <TColor.h>
#include <TChain.h>
#include <stdlib.h>
#include <TSpline.h>
#include <TChain.h>
#include <Getline.h>
#include <TKey.h>

TH2F *h;
TGraph *linea[500];
TSpline3 *spline[500];
TF1 *funz[500];
TCanvas *c;
int nlinee;
int indice;
int indicec;
//void color();
int ica;
char nomefil[500];
char nomefilc[500];
double vecx[500];
double vecy[500];
char run1[13];
char run2[13];
int sec,strip,csi;

int zval[500],aval[500];
float mm0;
char fapri[500];
char ftit[500];
char Nfil[500];
double iper(double *x,double *par);
double retta(double *x,double *par);

int ordina(int value);
void mettil(Int_t event, Int_t x, Int_t y, TObject *selected);
float pid(float fast,float slow,int *code);
float pidzlines(float x,float y);
double spezzata(double *x, double *par);
void copia(int indice);
float scaprod(float *px,float *py);
float calcola_val(int jlinea,float x);
void mettic(Int_t event, Int_t x, Int_t y, TObject *selected);
TGraph *gcont[100];
int ncont;
//int icod;
char icod[100];
char histoname[1000];
//ih=0 Si1-Si2 ih=1 Si-CsI
void metti_z(char *fileroot="histopsa_allgood.root",char *nomehisto="ID_SI1_11")

//void metti_z(int isec=1,int istrip=1,int icsi=1,char *nfil="RCO_SivsCsi_17ott_ca40",char *rini="111207161336",char *rfin="111211073832",char *file="lineez.dat",char *filec="contorno")
//void metti_z(int isec=1,int istrip=1,int icsi=1,char *nfil="RCO_test_SivsCsi_11mar14_ca48_histo",char *rini="111203002353",char *rfin="111207161335",char *file="lineez.dat",char *filec="contorno")
{
 //  sprintf(Nfil,"%s",nfil);
//   sec=isec;
//   strip=istrip;
//   csi=icsi; 
   sprintf(run1,"000000000000");
   sprintf(run2,"999999999999");
   sprintf(histoname,"%s",nomehisto);

//   sprintf(fapri,"rco_set%d_strip%d_csi%d_sivscsi_devse",isec,istrip,icsi);  
//   sprintf(ftit,"event->ringco->settore%d->strip%d->csi%d->sivscsi|dEvsE",isec,istrip,icsi);

  char filec[100];
  char file[100];
  char nfil[500];

  sprintf(fapri,"%s",nomehisto);
  TString ss=fileroot;
  int kk=ss.Index(".root"); 
  TString ss2=ss(0,kk);
  sprintf(Nfil,"%s",ss2.Data());
sprintf(nfil,"%s",ss2.Data());
  
  sprintf(file,"lineez.dat");
  sprintf(filec,"contorno");

 TIter next(gROOT->GetListOfSpecials());
 TObject *obj;
 while(obj=(TObject*)next())
   {if(obj->InheritsFrom("Mybar"))
     {Mybar *bar=(Mybar *)obj;
   bar->~Mybar();
     }
   }
for(int j=0;j<500;j++)
   {

     if(linea[j]!=0)
       {
	 delete linea[j];
       }
     linea[j]=0;
     if(funz[j]!=0)
       {
	 delete funz[j];
       }
     funz[j]=0;

     if(spline[j]!=0)
       {
	 delete spline[j];
       }
     spline[j]=0;

   }
 nlinee=0;
for(int j=0;j<100;j++)
   {

     if(gcont[j]!=0)
       {
	 delete gcont[j];
       }
     gcont[j]=0;
   }
 ncont=0;
TCanvas *cc=(TCanvas*) gROOT->FindObject("c");
 if(cc!=0)
   {
     delete cc;
     cc=0;
   }
 cc=(TCanvas *)gROOT->FindObject("cloc");
 if(cc!=0)
   {
 delete cc;
     cc=0;
   }

 cc=(TCanvas *)gROOT->FindObject("cloc2");
 if(cc!=0)
   {
 delete cc;
     cc=0;
   }

c=new TCanvas("c","c",200,0,500,500);
 c->Draw();
 // color();
 sprintf(nomefil,"%s_%s_%s",nfil,fapri,file);
 sprintf(nomefilc,"%s_%s_%s",nfil,fapri,filec);

 TFile *f=new TFile(Form("%s",fileroot));


		    //  TFile *f=new TFile(Form("%s.root",nfil));
 


     h=(TH2F*)f->Get(fapri);



 if(h==0)
   {
     cout<<"l'istogramma non esiste"<<endl;
     return;
   }
 // h->Rebin2D(2,2);
     sprintf(ftit,"%s",h->GetTitle());	
	  h->SetContour(128);
	  h->SetMaximum(50);
	  h->SetMinimum(1);

	  h->SetStats(kFALSE);
	  gPad->SetLogz(kTRUE);
	  h->Draw("zcol");
	  mm0=h->GetMaximum();
	  gPad->Modified();
	  gPad->Update();

 Mybar *bar = new Mybar("bar","vertical","zeta",10,0);
 bar->AddButton("mettilinea","lineaz()","mette una linea","button");
 bar->AddButton("duplicalinea","duplica()","duplica una linea","button");

 bar->AddButton("scrivilinee","scrivilinee()","scrive tutte le linee","button");
 bar->AddButton("caricalinee","caricalinee()","carica le linee da file","button");
 bar->AddButton("identifica","identifica()","fa le identificazioni","button");
 bar->AddButton("unzoom","unzoom()","Unzoom","button");
bar->AddButton("salva_formato_clicca_zlines","formato()","salva nel formato di clicca_z_lines","button");
 bar->AddButton("intervallo_run_validita","intervallo()","button");
bar->AddButton("aggiungi_punto_iniziale","punto_iniziale()","aggiunge a tutte le linee un punto fino x=0","button");
bar->AddButton("prolunga_lungo_iperbole","iperbole()","prolunga lungo un'iperbole tutte le linee","button");
bar->AddButton("Recupera_formato_clicca_zlines","recupera()","recupera dal formato di clicca_z_lines","button");
 bar->AddButton("stiracchia in x e y","stirax()","stira in x e y e shifta","button");
bar->AddButton("rimuovi_punto_iniziale","rimuovi_punto_iniziale()","toglie a tutte le linee il primo punto","button");
 bar->AddButton("metticontorno","contorno()","mette un contorno","button");
 bar->AddButton("leggicontorno","leggicontorno()","carica un contorno","button");
 bar->AddButton("scrivicontorni","scrivicontorni()","scrive tutti i contorni","button");
 bar->AddButton("carica_tutti_i_contorni","leggitutticontorni()","carica tutti i contorni","button");
 bar->AddButton("togli_contorno_da_schermo","togli_contorno()","toglie un contorno dal video");
 bar->AddButton("identifica_clicca_zlines","idzlines()","identifica in modo simile a clicca_z_lines","button");
 bar->AddButton("recupera_da_fileroot","recuperaroot()","recupera da file root linee e contorni","button");
 bar->AddButton("recupera_da_kaliveda","recuperakali()","recupera da kaliveda","button");
 bar->AddButton("scrivi_come_kaliveda","scrivikali()","scrivi kaliveda","button");
bar->AddButton("zoombassissimo","zoom01()","zoom","button");

bar->AddButton("zoombasso","zoom1()","zoom","button");
bar->AddButton("zoomstrettosux","zoom2()","zoom","button");
bar->AddButton("zoomalto","zoom3()","zoom","button");
bar->AddButton("zoombassoapertox","zoom4()","zoom","button");
 bar->AddButton("zoomquasialtissimo","zoom6()","button");
bar->AddButton("zoomaltissimo","zoom5()","zoom","button");


bar->AddButton("ripristinacol","ripristina()","ripristina","button");

//bar->AddButton("daicodice","codice()","Codice","button");

 bar->Show();

}

void lineaz()
{
  static char nomec[500];
  Mydialogo *input1=new Mydialogo("PI della linea?",Form("PI %d",(nlinee+1)*10),nomec);
  
  cout<<nomec<<endl;
  TString sc=nomec;
  int pival,pic;
sscanf(&sc.Data()[sc.Index("PI")+2],"%d",&pic);
  indice=nlinee;

   cout<<"Mettere la linea "<<nomec<<endl;
   for(int j=0;j<nlinee;j++)
     {
       if(linea[j]!=0)
	 {
	   TString sl=linea[j]->GetName();
	   sscanf(&sl.Data()[sl.Index("PI")+2],"%d",&pival);
	   if(pival==pic)
	 {
	   indice=j;
	   break;
	 }
	 }
     }

   cout<<"indice="<<indice<<" nlinee= "<<nlinee<<endl;
   cout<<aval[indice]<<endl;
  ica=0;
  int flag=0;
  if(linea[indice]==0)
    {
      cout<<"entro qui"<<endl;
      linea[indice]=new TGraph();
      linea[indice]->SetMarkerStyle(20);
      linea[indice]->SetMarkerColor(2);
      linea[indice]->SetMarkerSize(0.5);
      linea[indice]->SetName(Form("%s",nomec));
      //   spline[indice]=new TSpline3(Form("%s_spline",nomec),linea[indice]);
      // cout<<"spline qui="<<spline[indice]<<endl;
      flag=1;
      //      char aa[20];
      // sscanf(nomec,"%s %d",aa,&zval[indice]);
      // aval[indice]=2*zval[indice];
      // cout<<nomec<<" Z="<<zval[indice]<<" A="<<aval[indice]<<endl;
      zval[indice]=pic/10;
      aval[indice]=-1;
    }

      for(int j=0;j<linea[indice]->GetN();j++)
	{
	  vecx[j]=linea[indice]->GetX()[j];
	  vecy[j]=linea[indice]->GetY()[j];

	}
      cout<<"clic sx per mettere punti; c per uscire dalla linea; d per cancellare linea; u per togliere l'ultimo punto messo; x per togliere il primo punto della linea"<<endl;
      cout<<"n aggiunge un punto in fondo lungo una retta; f aggiunge un punto in cima lungo una retta"<<endl;
      cout<<"p prolunga una linea lungo una retta; m prolunga lungo un'iperbole; a prolunga all'inizio della linea lungo un'iperbole"<<endl;
      cout<<"g prolunga con fit iperbolico proiettando ortogonalmente alla tangente alla riga; k prolunga con fit iperbolico; q prolunga con fit rettilineo; r prolunga con fit rettilineo proiettando ortogonalmente alla tangente alla riga"<<endl;

  cout<<"+ per allungare in x; - per restringere in x; s per allungare in y; t per accorciare in y"<<endl;
  cout<<"> per shiftare a dx (in x); < per shiftare a sx (in x); h per shiftare in su; b per shiftare in giu"<<endl;

  c->cd();
c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
	      "mettil(Int_t,Int_t,Int_t,TObject*)");//connette il canvas a un evento grafico (non chiarissimo come funziona)
 while(ica!=1 && ica!=2 &&ica!=3)
   {usleep(100);
   gClient->HandleInput();//fondamentale, se no non funziona

   }

c->Disconnect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)");

// cout<<spline[nlinee]->GetName()<<endl;
// spline[nlinee]->Draw("same");
// gPad->Modified();
// gPad->Update();
 if(ica!=3)
 {
  cout<<"ho messo la linea "<<nomec<<endl;
 }
  if(flag==1)
    {
  nlinee++;
    }

  cout<<"nlinee="<<nlinee<<endl;
}

void scrivilinee()
{
  if(nlinee==0)
    {
      cout<<"Non ci sono linee da salvare"<<endl;
      return;
    }
  c->cd();
  static char nomec[500];
  Mydialogo *input2=new Mydialogo("Nome del file delle linee?",Form("%s",nomefil),nomec);
  if(system(Form("dir %s",nomec))==0)
    {
  static char nomec2[500];
       Mydialogo *input2=new Mydialogo("Il file esiste gia'. Sovrascrivo?","yes",nomec2);
       if(strcmp(nomec2,"yes")!=0)
	 {
return;
	 }
       system(Form("cp %s %s.old",nomec,nomec));//backup!
    }
   int ip,ih;
  char v1[2];
  char v2[2];
  int listapi[500];
  char pi[10];
  int pival[500];
  for(int j=0;j<nlinee;j++)
    {
      listapi[j]=j;
      char nomel[100];
      if(linea[j]!=0)
	{
      sprintf(nomel,"%s",linea[j]->GetName());
      sscanf(nomel,"%s %d",pi,&pival[j]);
	}
      else
	{
	  pival[j]=1000;
	}
    } 
  for(int j=0;j<nlinee-1;j++)
    {
      for(int k=j+1;k<nlinee;k++)
	{
	  if(pival[k]<pival[j])
	    {
	      int lj=listapi[j];
	      int pj=pival[j];
	      listapi[j]=listapi[k];
	      listapi[k]=lj;
	      pival[j]=pival[k];
	      pival[k]=pj;
	    }
	}
    }
 
char stringa[400];

  sprintf(stringa,"%s",nomefil);

  FILE *apri=fopen(nomec,"w");
  fprintf(apri,"# %s\n",stringa);
  for(int j=0;j<nlinee;j++)
    {
      int jj=listapi[j];
      if(linea[jj]!=0)
	{
      int flag=ordina(jj);
      if(flag==1)
	{
	  gPad->Modified();
	  gPad->Update();
	}
      if(linea[jj]->GetN()>0)
	{
      cout<<"salvo la linea "<<linea[jj]->GetName()<<endl;
      fprintf(apri,"%s\n",linea[jj]->GetName());
      double *xxj;
      double *yyj;
  xxj=linea[jj]->GetX();
  yyj=linea[jj]->GetY();


  for(int k=0;k<linea[jj]->GetN();k++)
    {

printf("%d %f %f\n",k,xxj[k],yyj[k]);      
fprintf(apri,"%d %f %f\n",k,xxj[k],yyj[k]);
      
    }
	}
	}
    }
  fclose(apri);

  TFile *froot=new TFile(Form("%s.root",nomefil),"RECREATE");
  froot->cd();
  TGraph *gr[nlinee];
  for(int j=0;j<nlinee;j++)
    {
      gr[j]=(TGraph*)linea[j]->Clone();
      TString sl=linea[j]->GetName();
      int ikk=sl.Index("PI");
      int pival;
      sscanf(&sl.Data()[ikk+2],"%d",&pival);

      gr[j]->SetName(Form("linea%d",pival));
      gr[j]->Write();
    } 
  froot->Write();
  froot->Close();
  for(int j=0;j<nlinee;j++)
    {
      delete gr[j];
      gr[j]=0;
    }
  cout<<"Ho scritto"<<endl;
}


void mettil(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  // cout<<event<<" "<<gPad->GetEventX()<<endl;
  c->cd();

  if(event==1)
    {
float px=gPad->AbsPixeltoX(x);
     float py=gPad->AbsPixeltoY(y);
     py=gPad->PadtoY(py);
     float uymin=gPad->GetUymin();
     float uymax=gPad->GetUymax();

     if(px>=gPad->GetUxmin() && px<=gPad->GetUxmax() && py>=gPad->PadtoY(uymin) && py<=gPad->PadtoY(uymax))
       {
	 // cout<<px<<" "<<py<<endl;
	 int n=linea[indice]->GetN();
	 // cout<<n<<endl;
	 linea[indice]->Set(n+1);
	 linea[indice]->SetPoint(n,px,py);


	 vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[n];
	 vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[n];

	 if(linea[indice]->GetN()==1)
	   {
	     linea[indice]->Draw("pl");

	   }


	 //	 gPad->Modified();
	 // gPad->Update();
       }
       }
  if(event==24 && gPad->GetEventX()==99)//c per uscire
    {
      copia(indice);
      gPad->Modified();
      gPad->Update();
      //      spline[nlinee]->Draw("same");
      
      int go=0;
go=ordina(indice);
    if(go==1)
      {

	     gPad->Modified();
      	     gPad->Update();

      }
  cout<<"fine linea "<<linea[indice]->GetName()<<endl;
      ica=2;
    }

  if(event==24 && gPad->GetEventX()==100) // d per cancellare linea
    {
      cout<<"Cancello la linea "<<linea[indice]->GetName()<<endl;
      delete linea[indice];
      linea[indice]=0;
      gPad->Modified();
      gPad->Update();

      ica=3;
    }

  if(event==24 && gPad->GetEventX()==117)// u per togliere l'ultimo punto
    {
      if(linea[indice]->GetN()>1)
	{
	  linea[indice]->RemovePoint(linea[indice]->GetN()-1);
	  gPad->Modified();
	  gPad->Update();
	  cout<<"tolgo l'ultimo punto messo"<<endl;
	  

	}
    }
  if(event==24 && gPad->GetEventX()==120)// x per togliere il primo
    {
      if(linea[indice]->GetN()>1)
	{
	  linea[indice]->RemovePoint(0);
	  gPad->Modified();
	  gPad->Update();
	  cout<<"tolgo il primo punto della linea"<<endl;
      for(int j=0;j<linea[indice]->GetN();j++)
	{
	  vecx[j]=linea[indice]->GetX()[j];
	  vecy[j]=linea[indice]->GetY()[j];

	}

	}
    }

  //a=97;b=98;c=99;f=102;g=103;n=110;
  if(event==24 &&gPad->GetEventX()==110)//n aggiungo un punto in fondo lungo la derivata
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=2)
	{
	  double xa,ya,xb,yb;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xa,ya);
	    linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	    double x0=xb;
	    double y0=yb;
	    double delta=(ya-yb)/(xa-xb);
	    double offset=(xa*yb-ya*xb)/(xa-xb);
	    if(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
	      {
			double xn=x0+(xb-xa);
		double yn=delta*xn+offset;
		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
		  {
		int np=linea[indice]->GetN();
		linea[indice]->Set(np+1);
		linea[indice]->SetPoint(np,xn,yn);
		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	
		  }
		x0=xn;
		y0=yn;
	      }
	    
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);	
	}
      else
	{
	  cout<<"occorrono almeno 2 punti"<<endl;
	}	
    }


 if(event==24 &&gPad->GetEventX()==102)//f aggiungo un punto in cima lungo la derivata
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=2)
	{
	  double xa,ya,xb,yb;
	    linea[indice]->GetPoint(0,xa,ya);
	    linea[indice]->GetPoint(1,xb,yb);
	    double x0=xa;
	    double y0=ya;
	    double delta=(ya-yb)/(xa-xb);
	    double offset=(xa*yb-ya*xb)/(xa-xb);
	    
	    if(x0>h->GetXaxis()->GetXmin() && y0>h->GetYaxis()->GetXmin())
	      {
			double xn=x0-(xb-xa);
		double yn=delta*xn+offset;
		//	cout<<xn<<" "<<yn<<endl;
		if(xn>h->GetXaxis()->GetXmin() && yn>h->GetYaxis()->GetXmin())
		  {
		int np=linea[indice]->GetN();
		linea[indice]->Set(np+1);
		linea[indice]->SetPoint(np,xn,yn);
		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	
		  }
		x0=xn;
		y0=yn;
	      
	    
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
	    ordina(indice);
	      }
	    if(x0-(xb-xa)<h->GetXaxis()->GetXmin())
	      {
		
		double xn=h->GetXaxis()->GetXmin();
		double yn=delta*xn+offset;
		int np=linea[indice]->GetN();
		linea[indice]->Set(np+1);
		linea[indice]->SetPoint(np,xn,yn);
		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
	    ordina(indice);
	    ica=2;
	      }
	
	}
      else
	{
	  cout<<"occorrono almeno 2 punti"<<endl;
	}	
    }


  if(event==24 && gPad->GetEventX()==112)//p per prolungare una linea lungo una retta
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=2)
	{
	  double xa,ya,xb,yb;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xa,ya);
	    linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	    double x0=xb;
	    double y0=yb;
	    double delta=(ya-yb)/(xa-xb);
	    double offset=(xa*yb-ya*xb)/(xa-xb);
	    while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
	      {
		double xn=x0+(xb-xa);
		double yn=delta*xn+offset;
		
		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
		  {
		int np=linea[indice]->GetN();
		linea[indice]->Set(np+1);
		linea[indice]->SetPoint(np,xn,yn);
		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	
		  }
		x0=xn;
		y0=yn;
	      }
	    
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
	}
      else
	{
	  cout<<"occorrono almeno 2 punti"<<endl;
	}

      ica=2;
    }

  if(event==24 && gPad->GetEventX()==109)//m per prolungare una linea lungo una iperbole
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      int incro=0;
      if(linea[indice]->GetN()>=5)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	  linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
 	  TF1 *fiper=new TF1("fiper",iper,xa,xb,2);
	  linea[indice]->Fit(fiper,"R");
	  char stringhetta[20];
	  sprintf(stringhetta,"%s",linea[indice]->GetName());
	  char pi[10];
	  int pival;
	  int pivalb;
	  sscanf(stringhetta,"%s %d",pi,&pival);
	  for(int jj=0;jj<nlinee;jj++)
	    {
	     if(jj!=indice)
	       {
		 if(linea[jj]!=0)
		   {

		     sprintf(stringhetta,"%s",linea[jj]->GetName());
		     sscanf(stringhetta,"%s %d",pi,&pivalb);
		     // cout<<"Check con linea PI"<<pivalb<<endl;
		     double xpp,ypp;
		     linea[jj]->GetPoint(linea[jj]->GetN()-1,xpp,ypp);
		     //cout<<TMath::Abs(xpp-h->GetXaxis()->GetXmax())/h->GetXaxis()->GetXmax()<<"ypp= "<<ypp<<" "<<fiper->Eval(h->GetXaxis()->GetXmax())<<endl;
		     if(((TMath::Abs(xpp-h->GetXaxis()->GetXmax())/h->GetXaxis()->GetXmax()<0.1)||(xpp>h->GetXaxis()->GetXmax())) &&pivalb<pival)
		       {
		     if(fiper->Eval(xpp)<ypp)	  
		       {
			 cout<<"Si incrocia la linea di PI inferiore="<<pivalb<<endl;
			 incro=1;
			 
			 delete fiper;
			 fiper=0;
			 gPad->Modified();
			 gPad->Update();
			 break;
		       }
		       }
		     if(((TMath::Abs(xpp-h->GetXaxis()->GetXmax())/h->GetXaxis()->GetXmax()<0.1||xpp>h->GetXaxis()->GetXmax()) )&&pivalb>pival)
		       {
		     if(fiper->Eval(xpp)>ypp)	  
		       {
			 cout<<"Si incrocia la linea di PI superiore="<<pivalb<<endl;
			 incro=1;
			 delete fiper;
			 fiper=0;
			 gPad->Modified();
			 gPad->Update();
			 break;
		       }
		       }

		   }
	       }
	    }
	  if(incro==0)
	    {

 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);

 	    while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
 		double xn=x0+(xb-xc);
		double yn=fiper->Eval(xn);
 		
		//cout<<xn<<" "<<yn<<endl;
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	
 		  }
 		x0=xn;
 		y0=yn;
 	      }
	    double xn=h->GetXaxis()->GetXmax();
	    double yn=fiper->Eval(xn);
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	    
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
	    delete fiper;
	    fiper=0;
	    }
	}
      else
	{
	  cout<<"occorrono almeno 5 punti"<<endl;
	}
      if(incro==0)
	{
      ica=2;
	}
    }

if(event==24 && gPad->GetEventX()==97)//a per prolungare una linea all'inizio lungo una iperbole
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=5)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(4,xb,yb);
	  linea[indice]->GetPoint(0,xa,ya);
 	  TF1 *fiper=new TF1("fiper",iper,xa,xb,2);
	  linea[indice]->Fit(fiper,"R");
	  

 	    double x0=xa;
 	    double y0=ya;
	    double xc,yc;
	    linea[indice]->GetPoint(1,xc,yc);

 	    while(x0>h->GetXaxis()->GetXmin() && y0>h->GetYaxis()->GetXmin())
 	      {
 		double xn=x0-(xc-xa);
		double yn=fiper->Eval(xn);
 		
		//cout<<xn<<" "<<yn<<endl;
 		if(xn>h->GetXaxis()->GetXmin() && yn>h->GetYaxis()->GetXmin())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
	
 		  }
 		x0=xn;
 		y0=yn;
 	      }
	    double xn=h->GetXaxis()->GetXmin();
	    double yn=fiper->Eval(xn);
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];

	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
	    ordina(indice);
	    delete fiper;
	    fiper=0;
	}
      else
	{
	  cout<<"occorrono almeno 5 punti"<<endl;
	}

      ica=2;
    }
 if(event==24 && gPad->GetEventX()==103)//g per prolungare una linea con un fit iperbolico proiettando ortogonalmente alla tangente alla riga
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=5)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	  linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
 	  
	  TCanvas *cnp=new TCanvas("cnp","cnp",800,0,250,250);	  
			
			cnp->Draw();

 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
	    float yn0=-1;
	    float delta0=xb-xc;
	     while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
		TF1 *fiper=new TF1("fiper",iper,xa,xb,2);
		linea[indice]->Fit(fiper,"R");
 		double xn=x0+(xb-xc);
		double yn=fiper->Eval((xn+x0)/2);
		double yderi=fiper->Derivative((xn+x0)/2);
		float vers1[2]={TMath::Cos(TMath::ATan(yderi)),TMath::Sin(TMath::ATan(yderi))};
		float vers2[2]={TMath::Cos(TMath::ATan(-1/yderi)),TMath::Sin(TMath::ATan(-1/yderi))};
		double bwx=h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin((xn+x0)/2));
		double bwy=h->GetYaxis()->GetBinWidth(h->GetYaxis()->FindBin(yn));
		double chx=(xn-x0)/bwx;
		double chy=(30.)/bwy;
		TH2F *hll=new TH2F("hll","hll",chx,(x0-xn)/2,(xn-x0)/2,chy,-15,+15);
		int ixmin=h->GetXaxis()->FindBin(x0);
		int ixmax=h->GetXaxis()->FindBin(xn);
		int iymin=h->GetYaxis()->FindBin(yn-15);
		int iymax=h->GetYaxis()->FindBin(yn+15);
		for(int ix=ixmin;ix<=ixmax;ix++)
		  {
		    float valx=h->GetXaxis()->GetBinCenter(ix)-(xn+x0)/2;
		    for(int iy=iymin;iy<=iymax;iy++)
		      {
		      float valy=h->GetYaxis()->GetBinCenter(iy)-yn;
			float hc=h->GetBinContent(ix,iy);
			float pcoord[2]={valx,valy};
			float newx=scaprod(pcoord,vers1);
			float newy=scaprod(pcoord,vers2);
			
			hll->Fill(newx,newy,hc);
		      }
		  }

		//cout<<xn<<" "<<yn<<" "<<x0<<" "<<xa<<" "<<xb<<" "<<xc<<endl;
		//cout<<yn0<<endl;


		//cout<<"integrale="<<proj->Integral()<<" chi_fit="<<fiper->GetChisquare()/fiper->GetNDF()<<endl;
		TH1D *proj=(TH1D*)hll->ProjectionY("proj",1,chx);

		if(proj->Integral()>20 && fiper->GetChisquare()/fiper->GetNDF()<100)
		  {


		    TF1 *fgaus=new TF1("fgaus","gaus",-15,+15);

		fgaus->SetParameter(1,0.);

		fgaus->SetParameter(0,proj->GetBinContent(proj->GetXaxis()->FindBin(0.)));
		int bmax=proj->GetMaximumBin();
		int jmezzo=-1;
		//		int blim=proj->GetXaxis()->FindBin(yn-15);
		int blim=proj->GetXaxis()->FindBin(-15.);
		for(int nm=bmax;nm>=blim;nm--)
		  {
		   
if(proj->GetBinContent(nm)<=0.5*proj->GetMaximum()
 		     &&proj->GetBinContent(nm-1)<proj->GetBinContent(nm)&&proj->GetBinContent(nm-2)<proj->GetBinContent(nm)&& proj->GetBinContent(nm-3)<proj->GetBinContent(nm))
 		   {
 		   jmezzo=nm;
 		    break;
 		   }
		  }


		if(jmezzo!=-1)
		  {
		    float delta=2*(proj->GetXaxis()->GetBinCenter(bmax)-proj->GetXaxis()->GetBinCenter(jmezzo));
		    	    fgaus->SetParameter(2,delta/2.35);
			    fgaus->SetParLimits(2,0.8*delta/2.35,1.2*delta/2.35);//FAZIA
			    if(1.2*delta/2.35>2)
			      {
			fgaus->SetParLimits(2,0.8,2.);//FAZIA	
			      }
			    cout<<"LLLLL "<<0.8*delta/2.35<<" "<<1.2*delta/2.35<<endl;
		  }
		//cout<<fgaus<<endl;
		cnp->cd();
		proj->Draw();
		//cout<<cnp<<endl;
		cnp->Modified();
		cnp->Update();
		//proj->Fit(fgaus,"R");
		proj->Fit(fgaus,"RB");
		cnp->Modified();
		cnp->Update();


		yn=fgaus->GetParameter(1)*vers2[1]+fiper->Eval((xn+x0)/2);
		xn=fgaus->GetParameter(1)*vers2[0]+(xn+x0)/2;

		//cout<<xn<<" "<<yn<<" "<<"chi="<<fgaus->GetChisquare()/fgaus->GetNDF()<<endl;

		//cout<<xn<<" "<<yn<<endl;
 
		//		if(fgaus->GetChisquare()/fgaus->GetNDF()>100 || (yn0>0 &&TMath::Abs(yn-yn0)/yn0>0.1)||(TMath::Abs(yn-fiper->Eval(xn))/fiper->Eval(xn)>0.01))
		//	if(fgaus->GetChisquare()/fgaus->GetNDF()>200 ||(TMath::Abs(yn-fiper->Eval(xn))/fiper->Eval(xn)>0.02))
		//			if(fgaus->GetChisquare()/fgaus->GetNDF()>500 ||(TMath::Abs(yn-fiper->Eval(xn))/fiper->Eval(xn)>0.02))
	if(fgaus->GetChisquare()/fgaus->GetNDF()>500 ||(TMath::Abs(yn-fiper->Eval(xn))/fiper->Eval(xn)>0.05))
		  {
		    cout<<"chisquare_gaus="<<fgaus->GetChisquare()/fgaus->GetNDF()<<" (yn-yn0)/yn0="<<TMath::Abs(yn-yn0)/yn0<<" yn-fiper->Eval(xn))/fiper->Eval(xn)="<<TMath::Abs(yn-fiper->Eval(xn))/fiper->Eval(xn)<<endl;
		 delete fiper;
		fiper=0;
		delete fgaus;
		fgaus=0;
		delete proj;
		proj=0;
		delete cnp;
		cnp=0;
		delete hll;
		hll=0;

		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;
		    
		    break;
		  }
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
		c->cd();
		gPad->Modified();
		gPad->Update();
		copia(indice);	
 		  }
		//	x0=xn+(xb-xc)/2;
 		x0=xn+delta0/2;
		//x0=xn+delta0;
		cout<<x0<<" "<<xb<<" "<<xc<<endl;
 		y0=yn;
		yn0=yn;
		linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
		linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
		linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
		delete fgaus;
		fgaus=0;

		  }
		
		else
		  {
		 delete fiper;
		fiper=0;
		delete proj;
		proj=0;
		if(cnp!=0)
		  {
		    delete cnp;
		    cnp=0;
		  }

		delete hll;
		hll=0;
	
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;

		    break;
		  }
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
		delete fiper;
		fiper=0;
		delete proj;
		proj=0;

		delete hll;
		hll=0;
	
	// 	if(cnp!=0)
// 		  {
// 		    delete cnp;
// 		    cnp=0;
// 		  }
	      }

	    


	}
      else
	{
	  cout<<"occorrono almeno 5 punti"<<endl;
	}

      ica=2;
    }
 if(event==24 && gPad->GetEventX()==114)//r per prolungare una linea con un fit rettilineo proiettando ortogonalmente alla tangente alla riga
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=4)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	  linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
 	  
	  TCanvas *cnp=new TCanvas("cnp","cnp",800,0,250,250);	  
			
			cnp->Draw();

 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
	    float yn0=-1;
	    float delta0=xb-xc;
	     while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
		TF1 *fretta=new TF1("fretta",retta,xa,xb,2);
		linea[indice]->Fit(fretta,"R");
 		double xn=x0+(xb-xc);
		double yn=fretta->Eval((xn+x0)/2);
		double yderi=fretta->Derivative((xn+x0)/2);
		float vers1[2]={TMath::Cos(TMath::ATan(yderi)),TMath::Sin(TMath::ATan(yderi))};
		float vers2[2]={TMath::Cos(TMath::ATan(-1/yderi)),TMath::Sin(TMath::ATan(-1/yderi))};
		double bwx=h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin((xn+x0)/2));
		double bwy=h->GetYaxis()->GetBinWidth(h->GetYaxis()->FindBin(yn));
		double chx=(xn-x0)/bwx;
		double chy=(30.)/bwy;
		TH2F *hll=new TH2F("hll","hll",chx,(x0-xn)/2,(xn-x0)/2,chy,-15,+15);
		int ixmin=h->GetXaxis()->FindBin(x0);
		int ixmax=h->GetXaxis()->FindBin(xn);
		int iymin=h->GetYaxis()->FindBin(yn-15);
		int iymax=h->GetYaxis()->FindBin(yn+15);
		for(int ix=ixmin;ix<=ixmax;ix++)
		  {
		    float valx=h->GetXaxis()->GetBinCenter(ix)-(xn+x0)/2;
		    for(int iy=iymin;iy<=iymax;iy++)
		      {
		      float valy=h->GetYaxis()->GetBinCenter(iy)-yn;
			float hc=h->GetBinContent(ix,iy);
			float pcoord[2]={valx,valy};
			float newx=scaprod(pcoord,vers1);
			float newy=scaprod(pcoord,vers2);
			
			hll->Fill(newx,newy,hc);
		      }
		  }

		//cout<<xn<<" "<<yn<<" "<<x0<<" "<<xa<<" "<<xb<<" "<<xc<<endl;
		//cout<<yn0<<endl;


		//cout<<"integrale="<<proj->Integral()<<" chi_fit="<<fretta->GetChisquare()/fretta->GetNDF()<<endl;
		TH1D *proj=(TH1D*)hll->ProjectionY("proj",1,chx);

		if(proj->Integral()>20 && fretta->GetChisquare()/fretta->GetNDF()<100)
		  {


		    TF1 *fgaus=new TF1("fgaus","gaus",-15,+15);

		fgaus->SetParameter(1,0.);

		fgaus->SetParameter(0,proj->GetBinContent(proj->GetXaxis()->FindBin(0.)));
		int bmax=proj->GetMaximumBin();
		int jmezzo=-1;
		//		int blim=proj->GetXaxis()->FindBin(yn-15);
		int blim=proj->GetXaxis()->FindBin(-15.);
		for(int nm=bmax;nm>=blim;nm--)
		  {
		   
if(proj->GetBinContent(nm)<=0.5*proj->GetMaximum()
 		     &&proj->GetBinContent(nm-1)<proj->GetBinContent(nm)&&proj->GetBinContent(nm-2)<proj->GetBinContent(nm)&& proj->GetBinContent(nm-3)<proj->GetBinContent(nm))
 		   {
 		   jmezzo=nm;
 		    break;
 		   }
		  }


		if(jmezzo!=-1)
		  {
		    float delta=2*(proj->GetXaxis()->GetBinCenter(bmax)-proj->GetXaxis()->GetBinCenter(jmezzo));
		    fgaus->SetParameter(2,delta/2.35);

		  }
		//cout<<fgaus<<endl;
		cnp->cd();
		proj->Draw();
		//cout<<cnp<<endl;
		cnp->Modified();
		cnp->Update();
		proj->Fit(fgaus,"R");
		cnp->Modified();
		cnp->Update();


		yn=fgaus->GetParameter(1)*vers2[1]+fretta->Eval((xn+x0)/2);
		xn=fgaus->GetParameter(1)*vers2[0]+(xn+x0)/2;

		//cout<<xn<<" "<<yn<<" "<<"chi="<<fgaus->GetChisquare()/fgaus->GetNDF()<<endl;

		//cout<<xn<<" "<<yn<<endl;
 
		//		if(fgaus->GetChisquare()/fgaus->GetNDF()>100 || (yn0>0 &&TMath::Abs(yn-yn0)/yn0>0.1)||(TMath::Abs(yn-fretta->Eval(xn))/fretta->Eval(xn)>0.01))
		//		if(fgaus->GetChisquare()/fgaus->GetNDF()>200 ||(TMath::Abs(yn-fretta->Eval(xn))/fretta->Eval(xn)>0.02))
		if(fgaus->GetChisquare()/fgaus->GetNDF()>500 ||(TMath::Abs(yn-fretta->Eval(xn))/fretta->Eval(xn)>0.02))
		  {
		    cout<<"chisquare_gaus="<<fgaus->GetChisquare()/fgaus->GetNDF()<<" (yn-yn0)/yn0="<<TMath::Abs(yn-yn0)/yn0<<" yn-fretta->Eval(xn))/fretta->Eval(xn)="<<TMath::Abs(yn-fretta->Eval(xn))/fretta->Eval(xn)<<endl;
		 delete fretta;
		fretta=0;
		delete fgaus;
		fgaus=0;
		delete proj;
		proj=0;
		delete cnp;
		cnp=0;
		delete hll;
		hll=0;

		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;
		    
		    break;
		  }
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
		c->cd();
		gPad->Modified();
		gPad->Update();
		copia(indice);	
 		  }
 		//x0=xn+(xb-xc)/2;
 		x0=xn+delta0/2;
		//x0=xn+delta0;
 		y0=yn;
		yn0=yn;
		linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
		linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
		linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
		delete fgaus;
		fgaus=0;

		  }
		
		else
		  {
		 delete fretta;
		fretta=0;
		delete proj;
		proj=0;
		if(cnp!=0)
		  {
		    delete cnp;
		    cnp=0;
		  }

		delete hll;
		hll=0;
	
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;

		    break;
		  }
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
		delete fretta;
		fretta=0;
		delete proj;
		proj=0;

		delete hll;
		hll=0;
	
	// 	if(cnp!=0)
// 		  {
// 		    delete cnp;
// 		    cnp=0;
// 		  }
	      }

	    


	}
      else
	{
	  cout<<"occorrono almeno 4 punti"<<endl;
	}

      ica=2;
    }

 if(event==24 && gPad->GetEventX()==107)//k per prolungare una linea con un fit
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=5)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	  linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
 	  
	  TCanvas *cnp=new TCanvas("cnp","cnp",800,0,250,250);	  
			
			cnp->Draw();
 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
	    float yn0=-1;
	     while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
		TF1 *fiper=new TF1("fiper",iper,xa,xb,2);
		linea[indice]->Fit(fiper,"R");
	

 		double xn=x0+(xb-xc);
		double yn=fiper->Eval((xn+x0)/2);
		//cout<<xn<<" "<<yn<<" "<<x0<<" "<<xa<<" "<<xb<<" "<<xc<<endl;
		//cout<<yn0<<endl;
		TH1D *proj=h->ProjectionY("proj",h->GetXaxis()->FindBin(x0),h->GetXaxis()->FindBin(xn));
		proj->GetXaxis()->SetRangeUser(yn-15,yn+15);
		cout<<"integrale="<<proj->Integral()<<" chi_fit="<<fiper->GetChisquare()/fiper->GetNDF()<<endl;
	
		if(proj->Integral()>20 && fiper->GetChisquare()/fiper->GetNDF()<100)
		  {

		    TF1 *fgaus=new TF1("fgaus","gaus",yn-10,yn+10);
		fgaus->SetParameter(1,yn);
		fgaus->SetParameter(0,proj->GetBinContent(proj->GetXaxis()->FindBin(yn)));
		int bmax=proj->GetMaximumBin();
		int jmezzo=-1;
		int blim=proj->GetXaxis()->FindBin(yn-15);
		for(int nm=bmax;nm>=blim;nm--)
		  {
		   
if(proj->GetBinContent(nm)<=0.5*proj->GetMaximum()
 		     &&proj->GetBinContent(nm-1)<proj->GetBinContent(nm)&&proj->GetBinContent(nm-2)<proj->GetBinContent(nm)&& proj->GetBinContent(nm-3)<proj->GetBinContent(nm))
 		   {
 		   jmezzo=nm;
 		    break;
 		   }
		  }


		if(jmezzo!=-1)
		  {
		    float delta=2*(proj->GetXaxis()->GetBinCenter(bmax)-proj->GetXaxis()->GetBinCenter(jmezzo));
		    fgaus->SetParameter(2,delta/2.35);
		  }
		//cout<<fgaus<<endl;
		cnp->cd();
		proj->Draw();
		//cout<<cnp<<endl;
		cnp->Modified();
		cnp->Update();
		proj->Fit(fgaus,"R");
		cnp->Modified();
		cnp->Update();
		xn=(xn+x0)/2;
		yn=fgaus->GetParameter(1);
		//cout<<xn<<" "<<yn<<" "<<"chi="<<fgaus->GetChisquare()/fgaus->GetNDF()<<endl;

		//cout<<xn<<" "<<yn<<endl;
		//const char *inp=Getline("qqq\n");
		
		if(fgaus->GetChisquare()/fgaus->GetNDF()>100 || (yn0>0 &&TMath::Abs(yn-yn0)/yn0>0.1))
		  {
		    cout<<"chisquare_gaus="<<fgaus->GetChisquare()/fgaus->GetNDF()<<" (yn-yn0)/yn0="<<TMath::Abs(yn-yn0)/yn0<<endl;
		 delete fiper;
		fiper=0;
		delete fgaus;
		fgaus=0;
		delete proj;
		proj=0;
		delete cnp;
		cnp=0;
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;
		    
		    break;
		  }
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
		c->cd();
		gPad->Modified();
		gPad->Update();
		copia(indice);	
 		  }
 		x0=xn+(xb-xc)/2;
 		y0=yn;
		yn0=yn;
		linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
		linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
		linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
		delete fgaus;
		fgaus=0;

		  }
		
		else
		  {
		 delete fiper;
		fiper=0;
		delete proj;
		proj=0;
		if(cnp!=0)
		  {
		    delete cnp;
		    cnp=0;
		  }
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;

		    break;
		  }
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
		delete fiper;
		fiper=0;
		delete proj;
		proj=0;
	// 	if(cnp!=0)
	//	  {
	//	    delete cnp;
	//	    cnp=0;
	//	  }
	      }

	    


	}
      else
	{
	  cout<<"occorrono almeno 5 punti"<<endl;
	}

      ica=2;
    }
if(event==24 && gPad->GetEventX()==113)//q per prolungare una linea con un fit lungo una retta
    {
      cout<<"prolungo la linea "<<linea[indice]->GetN()<<endl;
      if(linea[indice]->GetN()>=4)
	{

	  double xa,ya,xb,yb;
	  linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
	  linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
 	  
	  TCanvas *cnp=new TCanvas("cnp","cnp",800,0,250,250);	  
			
			cnp->Draw();
 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
	    float yn0=-1;
	     while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
		TF1 *fretta=new TF1("fretta",retta,xa,xb,2);
		linea[indice]->Fit(fretta,"R");
	

 		double xn=x0+(xb-xc);
		double yn=fretta->Eval((xn+x0)/2);
		//cout<<xn<<" "<<yn<<" "<<x0<<" "<<xa<<" "<<xb<<" "<<xc<<endl;
		//cout<<yn0<<endl;
		TH1D *proj=h->ProjectionY("proj",h->GetXaxis()->FindBin(x0),h->GetXaxis()->FindBin(xn));
		proj->GetXaxis()->SetRangeUser(yn-15,yn+15);
		cout<<"integrale="<<proj->Integral()<<" chi_fit="<<fretta->GetChisquare()/fretta->GetNDF()<<endl;
	
		if(proj->Integral()>20 && fretta->GetChisquare()/fretta->GetNDF()<100)
		  {

		    TF1 *fgaus=new TF1("fgaus","gaus",yn-10,yn+10);
		fgaus->SetParameter(1,yn);
		fgaus->SetParameter(0,proj->GetBinContent(proj->GetXaxis()->FindBin(yn)));
		int bmax=proj->GetMaximumBin();
		int jmezzo=-1;
		int blim=proj->GetXaxis()->FindBin(yn-15);
		for(int nm=bmax;nm>=blim;nm--)
		  {
		   
if(proj->GetBinContent(nm)<=0.5*proj->GetMaximum()
 		     &&proj->GetBinContent(nm-1)<proj->GetBinContent(nm)&&proj->GetBinContent(nm-2)<proj->GetBinContent(nm)&& proj->GetBinContent(nm-3)<proj->GetBinContent(nm))
 		   {
 		   jmezzo=nm;
 		    break;
 		   }
		  }


		if(jmezzo!=-1)
		  {
		    float delta=2*(proj->GetXaxis()->GetBinCenter(bmax)-proj->GetXaxis()->GetBinCenter(jmezzo));
		    fgaus->SetParameter(2,delta/2.35);
		  }
		//cout<<fgaus<<endl;
		cnp->cd();
		proj->Draw();
		//cout<<cnp<<endl;
		cnp->Modified();
		cnp->Update();
		proj->Fit(fgaus,"R");
		cnp->Modified();
		cnp->Update();
		xn=(xn+x0)/2;
		yn=fgaus->GetParameter(1);
		//cout<<xn<<" "<<yn<<" "<<"chi="<<fgaus->GetChisquare()/fgaus->GetNDF()<<endl;

		//cout<<xn<<" "<<yn<<endl;
		//const char *inp=Getline("qqq\n");
		
		if(fgaus->GetChisquare()/fgaus->GetNDF()>100 || (yn0>0 &&TMath::Abs(yn-yn0)/yn0>0.1))
		  {
		    cout<<"chisquare_gaus="<<fgaus->GetChisquare()/fgaus->GetNDF()<<" (yn-yn0)/yn0="<<TMath::Abs(yn-yn0)/yn0<<endl;
		 delete fretta;
		fretta=0;
		delete fgaus;
		fgaus=0;
		delete proj;
		proj=0;
		delete cnp;
		cnp=0;
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;
		    
		    break;
		  }
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[indice]->GetN();
 		linea[indice]->Set(np+1);
 		linea[indice]->SetPoint(np,xn,yn);
 		vecx[linea[indice]->GetN()-1]=linea[indice]->GetX()[np];
 		vecy[linea[indice]->GetN()-1]=linea[indice]->GetY()[np];
		c->cd();
		gPad->Modified();
		gPad->Update();
		copia(indice);	
 		  }
 		x0=xn+(xb-xc)/2;
 		y0=yn;
		yn0=yn;
		linea[indice]->GetPoint(linea[indice]->GetN()-1,xb,yb);
		linea[indice]->GetPoint(linea[indice]->GetN()-5,xa,ya);
		linea[indice]->GetPoint(linea[indice]->GetN()-2,xc,yc);
		delete fgaus;
		fgaus=0;

		  }
		
		else
		  {
		 delete fretta;
		fretta=0;
		delete proj;
		proj=0;
		if(cnp!=0)
		  {
		    delete cnp;
		    cnp=0;
		  }
		linea[indice]->RemovePoint(linea[indice]->GetN()-1);
		      for(int j=0;j<linea[indice]->GetN();j++)
			{
			  vecx[j]=linea[indice]->GetX()[j];
			  vecy[j]=linea[indice]->GetY()[j];

			}
		    ica=2;

		    break;
		  }
	    gPad->Modified();
	    gPad->Update();
	    copia(indice);
		delete fretta;
		fretta=0;
		delete proj;
		proj=0;
	// 	if(cnp!=0)
// 		  {
// 		    delete cnp;
// 		    cnp=0;
// 		  }
	      }

	    


	}
      else
	{
	  cout<<"occorrono almeno 4 punti"<<endl;
	}

      ica=2;
    }

 if(event==24 && gPad->GetEventX()==43)//+ per stirare (a crescere) la linea su x
   {
     if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  xx[k]=(xx[k]-xx[0])*1.001;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		 	 vecx[k]=xx[k];
			 vecy[k]=yy[k]; 
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==45)//- per restringere (a diminuire) la linea su x
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  xx[k]=(xx[k]-xx[0])/1.001;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==115)//s per stirare (a crescere) la linea su y
   {
     if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  yy[k]=yy[k]*1.001;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		 	 vecx[k]=xx[k];
			 vecy[k]=yy[k]; 
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==116)//t per restringere (a diminuire) la linea su y
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  yy[k]=yy[k]/1.001;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }

 if(event==24 && gPad->GetEventX()==104)//h per shiftare in su
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  yy[k]=yy[k]+5.;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==98)//b per shiftare in giu'
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  yy[k]=yy[k]-5.;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==62)//> per shiftare a dx in x
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  xx[k]=xx[k]+5.;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==60)//< per shiftare a sx in x
   {
  if(linea[indice]->GetN()>0)
       {
	 double *xx=linea[indice]->GetX();
	      double *yy=linea[indice]->GetY();

	      for(int k=0;k<linea[indice]->GetN();k++)
		{
		  xx[k]=xx[k]-5.;
		  linea[indice]->SetPoint(k,xx[k],yy[k]);
		  vecx[k]=xx[k];
		  vecy[k]=yy[k]; 
		  
		}
	      gPad->Modified();
	      gPad->Update();
       }
   }
 if(event==24 && gPad->GetEventX()==119)//w dare z e a
  {
    cout<<linea[indice]->GetName()<<endl;
  static char nomec[500];
  Mydialogo *input1=new Mydialogo("Z,A?",Form("%d,%d",zval[indice],aval[indice]),nomec);

  
  sscanf(nomec,"%d,%d",&zval[indice],&aval[indice]);
    if(aval[indice]==-1)
      {
	aval[indice]=2*zval[indice];
      }
     cout<<"Z="<<zval[indice]<<" A="<<aval[indice]<<endl;

     int pl;
     char aal[100];
     sscanf(linea[indice]->GetName(),"%s %d",aal,&pl);
     linea[indice]->SetName(Form("PI %d Z=%d A=%d",pl,zval[indice],aval[indice]));
    c->cd();
  }

}




void copia(int indice)
{
  for(int j=0;j<linea[indice]->GetN();j++)
    {
     
      linea[indice]->SetPoint(j,vecx[j],vecy[j]);
    }
  return;
}




int ordina(int value)
{

  float vx[500],vy[500],servx,servy;
  for(int j=0;j<linea[value]->GetN();j++)
    {
      double xx,yy;
      linea[value]->GetPoint(j,xx,yy);
      vx[j]=xx;
      vy[j]=yy;
    }
  int iflag=0;
  for(int j=0;j<linea[value]->GetN()-1;j++)
    {
      for(int k=j;k<linea[value]->GetN();k++)
	{
	  if(vx[j]>vx[k])
	    {
	      iflag=1;
	      float servx=vx[j];
	      float servy=vy[j];
	      vx[j]=vx[k];
	      vy[j]=vy[k];
	      vx[k]=servx;
	      vy[k]=servy;
	    }
	}

    }
  int out=0;
  if(iflag==1)
    {
      for(int j=0;j<linea[value]->GetN();j++)
	{
	  linea[value]->SetPoint(j,vx[j],vy[j]);
	  // cout<<j<<" "<<vx[j]<<" "<<vy[j]<<endl;
	}

      out=1;
      cout<<"riordino la linea "<<linea[value]->GetName()<<endl;
    }

  return out;
}
void caricalinee()
{
 static char nomec[500];
  if(nlinee>0)
    {
  Mydialogo *input1=new Mydialogo("Attenzione, cancello tutte le linee presenti; continuare?","yes",nomec);
  if(strcmp(nomec,"yes")!=0)
    {
      return;
    }     

  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
      delete linea[j];
      linea[j]=0;
	}
      if(funz[j]!=0)
	{
	  delete funz[j];
	  funz[j]=0;
	}
    }
  nlinee=0;
    }

  Mydialogo *input2=new Mydialogo("Nome del file delle linee?",Form("%s",nomefil),nomec);
  FILE *apri=fopen(nomec,"r");
  if(apri==0)
    {
      cout<<nomefil<<" non esiste"<<endl;
      return;
    }
  char riga[500];
  int valori[500],npunti[500];
  char pi[10];
  int pival;
  float xx[500][500];
  float yy[500][500];
  for(int j=0;j<500;j++)
    {
      valori[j]=0;
      npunti[j]=0;
      for(int k=0;k<500;k++)
	{
	  xx[j][k]=0;
	  yy[j][k]=0;
	}
    }
  int iflag=0;
  int ip;
  fscanf(apri,"\n%[^\n]",riga);
  TString sr;
  while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
    {
      sr=riga;
      sscanf(riga,"%s",pi);
      //      cout<<riga<<" "<<pi<<endl;
      if(strcmp(pi,"PI")==0)
	{
	  sscanf(riga,"%s %d",pi,&pival);
	  //  cout<<"pi="<<pi<<" "<<pival<<endl;
	  valori[nlinee]=pival;
	  zval[nlinee]=-1;
	  aval[nlinee]=-1;
	  if(sr.Contains("Z="))
	    {

	      sscanf(&sr.Data()[sr.Index("Z=")+2],"%d",&zval[nlinee]);
	      sscanf(&sr.Data()[sr.Index("A=")+2],"%d",&aval[nlinee]);
	    }
	  
	  nlinee++;
	}
      else
	{
	  sscanf(riga,"%d %f %f",&ip,&xx[nlinee-1][npunti[nlinee-1]],&yy[nlinee-1][npunti[nlinee-1]]);

	  //	  if(ipp==112){yy[nlinee-1][npunti[nlinee-1]]=yy[nlinee-1][npunti[nlinee-1]]/1.005;}

	  npunti[nlinee-1]++;
	}

    }
  fclose(apri);
  cout<<"nlinee="<<nlinee<<endl;
  int colonna[500];
  for(int j=0;j<500;j++)
    {
      colonna[j]=j;
    }
  for(int j=0;j<nlinee-1;j++)
    {
      for(int k=j+1;k<nlinee;k++)
	{
	  if(valori[j]>valori[k])
	    {

	      int li=colonna[k];
	      colonna[k]=colonna[j];
	      colonna[j]=li;
	      	      float serv=valori[k];
	       valori[k]=valori[j];
	      valori[j]=serv;
	    }
	}
    }
  c->cd();
  for(int j=0;j<nlinee;j++)
    {
      
      linea[j]=new TGraph(npunti[colonna[j]],xx[colonna[j]],yy[colonna[j]]);
      if(zval[j]<0)
	{
      linea[j]->SetName(Form("PI %d",valori[j]));
	}
      else
	{
	  linea[j]->SetName(Form("PI %d Z=%d A=%d",valori[j],zval[j],aval[j]));
	}
      cout<<linea[j]->GetName()<<endl;
      linea[j]->SetMarkerStyle(20);
      linea[j]->SetMarkerColor(2);
      linea[j]->SetMarkerSize(0.5);
      linea[j]->Draw("pl");
    }
  gPad->Modified();
  gPad->Update();
}

void identifica()
{
  cout<<"perche' funzioni bisogna che le linee siano in ordine di PI crescente! se non lo sono, passare dal ciclo scrivi - leggi -scrivi"<<endl;
static char nomec[500];
  Mydialogo *input1=new Mydialogo("Sei passato dal ciclo scrivi - leggi -scrivi?","yes",nomec);
  if(strcmp(nomec,"yes")!=0)
    {
      return;
    }

  int conta=0;
  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
	  conta++;
	}
    }
  if(conta<2)
    {
  Mydialogo *input1=new Mydialogo("Non ci sono almeno 2 linee in canna","OK?",nomec);
  return;
    }
  c->cd();
  double xx,yy,a;
  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
	  funz[j]=new TF1(Form("f %s",linea[j]->GetName()),spezzata,0,4000,2*linea[j]->GetN()+1);
	  cout<<funz[j]->GetName()<<endl;


	  a=linea[j]->GetN();
	  funz[j]->SetParameter(0,a);
	  for(int k=0;k<linea[j]->GetN();k++)
	    {
	      linea[j]->GetPoint(k,xx,yy);
	      funz[j]->SetParameter(k+1,xx);
	      funz[j]->SetParameter(k+1+linea[j]->GetN(),yy);
	    }	  
	  funz[j]->SetLineColor(5);
	  //funz[j]->Draw("same");
	}
    }
  gPad->Modified();
  gPad->Update();




TCanvas *cloc=(TCanvas*) gROOT->FindObject("cloc");
 if(cloc!=0)
   {
     delete cloc;
     cloc=0;
   }

  cloc=new TCanvas("cloc","cloc",500,0,500,500);
  cloc->Draw();
  cloc->Divide(2,2);
  cloc->cd(1);
  gPad->SetLogy(kTRUE);
TCanvas *cloc2=(TCanvas*) gROOT->FindObject("cloc2");
 if(cloc2!=0)
   {
     delete cloc2;
     cloc2=0;
   }

  cloc2=new TCanvas("cloc2","cloc2",1000,0,500,500);
  cloc2->Draw();
  cloc2->Divide(2,2);
  TH1D *hpi=(TH1D*)gROOT->FindObject("hpi");
  if(hpi!=0)
    {
      delete hpi;
      hpi=0;
    }
  TH2F *hnoid=(TH2F*)gROOT->FindObject("hnoid");
  if(hnoid!=0)
    {
      delete hnoid;
      hnoid=0;
    }
  TH2F *hid=(TH2F*)gROOT->FindObject("hid");
  if(hid!=0)
    {
      delete hid;
      hid=0;
    }
  TH2F *hidcode0=(TH2F*)gROOT->FindObject("hidcode0");
  if(hidcode0!=0)
    {
      delete hidcode0;
      hidcode0=0;
    }
  TH2F *hpie=(TH2F*)gROOT->FindObject("hpie");
  if(hpie!=0)
    {
      delete hpie;
      hpie=0;
    }
  TH2F *hli=(TH2F*)gROOT->FindObject("hli");
  if(hli!=0)
    {
      delete hli;
      hli=0;
    }
  TH2F *hBe=(TH2F*)gROOT->FindObject("hBe");
  if(hBe!=0)
    {
      delete hBe;
      hBe=0;
    }
  TH2F *halfa=(TH2F*)gROOT->FindObject("halfa");
  if(halfa!=0)
    {
      delete halfa;
      halfa=0;
    }


  TCutG *gcut=(TCutG*)gROOT->FindObject(Form("%s_%s_contorno999.dat",Nfil,fapri));
  if(gcut!=0)
    {
      delete gcut;
      gcut=0;
    }


  FILE *ap=fopen(Form("%s_%s_contorno999.dat",Nfil,fapri),"r");
		 if(ap==0)
		   {
		     cout<<"non esiste il contorno per gli stoppati"<<endl;
		   }
		 else
		   {
		     gcut=new TCutG();
		     gcut->SetName(Form("%s_%s_contorno999.dat",Nfil,fapri));
		     gcut->SetLineColor(2);
		     char rr[400];
		     fscanf(ap,"\n%[^\n]",rr);
		     int bnp;
		     float bnx,bny;
		     while(fscanf(ap,"%d %f %f",&bnp,&bnx,&bny)!=EOF)
		       {
			 cout<<bnp<<" "<<bnx<<" "<<bny<<endl;
			 gcut->Set(bnp+1);
			 gcut->SetPoint(bnp,bnx,bny);
		       }
		     fclose(ap);
		   }

  hpi=new TH1D("hpi","hpi",2500,0,250);
  hnoid=new TH2F("hnoid","hnoid",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  hli=new TH2F("hli","hli",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  halfa=new TH2F("halfa","halfa",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  hBe=new TH2F("hBe","hBe",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());

  hid=new TH2F("hid","hid",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
   hidcode0=new TH2F("hidcode0","hidcode0",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  hpie=new TH2F("hpie","hpie",hpi->GetXaxis()->GetNbins(),hpi->GetXaxis()->GetXmin(),hpi->GetXaxis()->GetXmax(),hid->GetYaxis()->GetNbins(),hid->GetYaxis()->GetXmin(),hid->GetYaxis()->GetXmax());

  int code;
  float pi;

  for(int j=0;j<h->GetXaxis()->GetNbins();j++)
    {
	  float gb=h->GetXaxis()->GetBinCenter(j+1);
      for(int k=0;k<h->GetYaxis()->GetNbins();k++)
	{
	  if(h->GetBinContent(j+1,k+1)>0)
	    {	  
	  float ga=h->GetYaxis()->GetBinCenter(k+1);
	  if((gcut!=0&& gcut->IsInside(gb,ga)==0) || gcut==0)
	    {

	  pi=pid(ga,gb,&code);
	 
	  if(pi<0)
	    {
	      hnoid->Fill(gb,ga,h->GetBinContent(j+1,k+1));
	    }
	  else
	    {
	      //	      if(slowaperto>100)
	      //	{
	      hpi->Fill(pi,h->GetBinContent(j+1,k+1));
	      hpie->Fill(pi,ga,h->GetBinContent(j+1,k+1));
	      if(pi>25 && pi<=35){hli->SetBinContent(j+1,k+1,h->GetBinContent(j+1,k+1));}
	      if(pi>18 && pi<=25){halfa->SetBinContent(j+1,k+1,h->GetBinContent(j+1,k+1));}
	      if(pi>35 && pi<=45){hBe->SetBinContent(j+1,k+1,h->GetBinContent(j+1,k+1));}

	      //	}
	      if(code==1)
		{
		  hid->Fill(gb,ga,h->GetBinContent(j+1,k+1));
		}
	      //	      if(code==0)
	      if(code==2)
		{
		  hidcode0->Fill(gb,ga,h->GetBinContent(j+1,k+1));
		}

	    }
	}
	    }
    }
    }
  cloc->cd();
  cloc->cd(1);
  hpi->Draw();
  gPad->Modified();
  gPad->Update();
  cloc->cd(2);
  hid->Draw("zcol");
//   for(int n=0;n<nlinee;n++)
//     {
//       if(linea[n]!=0)
// 	{
// 	  linea[n]->Draw();
// 	}
//     }
  gPad->Modified();
  gPad->Update();
  cloc->cd(3);
  hidcode0->Draw("zcol");
  gPad->Modified();
  gPad->Update();

  cloc->cd(4);
  hnoid->Draw("zcol");
  gPad->Modified();
  gPad->Update();
  cloc2->cd();
  cloc2->cd(1);
  hpie->Draw("zcol");
  gPad->Modified();
  gPad->Update();
  cloc2->cd(2);
  halfa->Draw("zcol");
  for(int j=0;j<5;j++)
    {
      if(linea[j]!=0)
	{
	  linea[j]->Draw();
	}
    }
  gPad->Modified();
  gPad->Update();
  cloc2->cd(3);
  hli->Draw("zcol");
  for(int j=0;j<5;j++)
    {
      if(linea[j]!=0)
	{
	  linea[j]->Draw();
	}
    }
  gPad->Modified();
  gPad->Update();
  cloc2->cd(4);
  hBe->Draw("zcol");
  for(int j=0;j<5;j++)
    {
      if(linea[j]!=0)
	{
	  linea[j]->Draw();
	}
    }
  gPad->Modified();
  gPad->Update();


  c->cd();


}
float pid(float fast,float slow, int *code)
{
  *code=-1;
  float pi=-1;
  if(fast<funz[0]->Eval(slow)||fast>funz[nlinee-1]->Eval(slow))
    {
      //      pi=-1;
      //  return pi;
    }
  double xm1,ym1,xM1,yM1,xf1,yf1,xF1,yF1;
  linea[0]->GetPoint(0,xm1,ym1);
  linea[nlinee-1]->GetPoint(0,xM1,yM1);
  linea[0]->GetPoint(linea[0]->GetN()-1,xf1,yf1);
  linea[nlinee-1]->GetPoint(linea[nlinee-1]->GetN()-1,xF1,yF1);

  if(slow<xm1 || slow<xM1)
    {
      // pi=-1;
      // return pi;
    }

  if(slow>xf1 || slow>xF1)
    {
      // pi=-1;
      // return pi;
    }



  for(int j=0;j<nlinee-1;j++)
    {
      double ya=funz[j]->Eval(slow);
      double yb=funz[j+1]->Eval(slow);
      double yc=-1000;
      double yd=-1000;
      double p0,pn,yp0,ypn;
      linea[j]->GetPoint(0,p0,yp0);
      linea[j]->GetPoint(linea[j]->GetN()-1,pn,ypn);
      double P0,Pn,yP0,yPn;
      linea[j+1]->GetPoint(0,P0,yP0);
      linea[j+1]->GetPoint(linea[j+1]->GetN()-1,Pn,yPn);
      if(fast>=ya && fast<yb && slow>=p0 && slow<=pn && slow>=P0 && slow<=Pn)
	//  if(fast>=ya && fast<yb && slow>=funz[j]->GetParameter(1) && slow<=funz[j]->GetParameter(funz[j]->GetParameter(0)))
	{
	  char aa[2];
	  char bb[3];
	  int pi1,pi2;
	  if(j>0 && j<nlinee-2)
	    {
	      yc=funz[j-1]->Eval(slow);
	      double PPn,yPPn,PP0,yPP0;
	      linea[j-1]->GetPoint(0,PP0,yPP0);
	      linea[j-1]->GetPoint(linea[j-1]->GetN()-1,PPn,yPPn);

	      float distya=fast-ya;
	      float distyb=yb-fast;
	      if(fabs(yb-ya)>fabs(ya-yc) && fabs(distya)<fabs(distyb) && slow>=p0 && slow<=pn && slow>=PP0 && slow<=PPn)
		{

	
	  sscanf(funz[j]->GetName(),"%1s %2s %d",aa,bb,&pi1);
	  sscanf(funz[j-1]->GetName(),"%1s %2s %d",aa,bb,&pi2);
	  pi=pi2+(fast-yc)/(ya-yc)*(pi1-pi2);
	  *code=2;
	  return pi;
		}
	      yd=funz[j+2]->Eval(slow);
	      float distyd=fast-yd;
	      linea[j+2]->GetPoint(0,PP0,yPP0);
	      linea[j+2]->GetPoint(linea[j-1]->GetN()-1,PPn,yPPn);

	      if(fabs(yb-ya)>fabs(ya-yd) && fabs(distya)>fabs(distyb) && slow>=p0 && slow<=pn && slow>=PP0 && slow<=PPn)
		{
	  sscanf(funz[j+1]->GetName(),"%1s %2s %d",aa,bb,&pi1);
	  sscanf(funz[j+2]->GetName(),"%1s %2s %d",aa,bb,&pi2);
	  pi=pi1+(fast-yb)/(yd-yb)*(pi2-pi1);
	  *code=2;
	  return pi;
		}


	    }



	  sscanf(funz[j]->GetName(),"%1s %2s %d",aa,bb,&pi1);
	  sscanf(funz[j+1]->GetName(),"%1s %2s %d",aa,bb,&pi2);

	  pi=pi1+(fast-ya)/(yb-ya)*(pi2-pi1);
	  //	  cout<<pi1<<" "<<pi2<<" "<<fast<<" "<<slow<<" "<<ya<<" "<<yb<<endl;
	  // cout<<pi<<endl;
	  *code=1;
	  return pi;
	}
    }
  int jmin=-1;
  for(int j=0;j<nlinee-1;j++)
    {
      double ya=funz[j]->Eval(slow);
      double x0,y0,xf,yf;
      linea[j]->GetPoint(0,x0,y0);
      linea[j]->GetPoint(linea[j]->GetN()-1,xf,yf);

      if(fast>=ya &&slow>=x0 && slow <=xf)
	{
	  jmin=j;
	  //break;
	}
    }
  int jmax=-1;
  if(jmin>=0)
    {
  for(int j=jmin+1;j<nlinee;j++)
    {
      double ya=funz[j]->Eval(slow);
      double x0,y0,xf,yf;
      linea[j]->GetPoint(0,x0,y0);
      linea[j]->GetPoint(linea[j]->GetN()-1,xf,yf);

      if(fast<ya &&slow>=x0 && slow<=xf)
	{
	  jmax=j;
	    break;
	}

    }
  if(jmax>=0)
    {
      double ya=funz[jmin]->Eval(slow);
      double yb=funz[jmax]->Eval(slow);

	  char aa[2];
	  char bb[3];
	  int pi1,pi2;
	  sscanf(funz[jmin]->GetName(),"%1s %2s %d",aa,bb,&pi1);
	  sscanf(funz[jmax]->GetName(),"%1s %2s %d",aa,bb,&pi2);

	  pi=pi1+(fast-ya)/(yb-ya)*(pi2-pi1);
	  //	  cout<<pi1<<" "<<pi2<<" "<<fast<<" "<<slow<<" "<<ya<<" "<<yb<<endl;
	  // cout<<pi<<endl;
	  *code=0;
	  return pi;
	}

    }


    



  return pi;
}

double spezzata(double *x,double *par)
{
  int npoint=(int)par[0];
  double out=0;
  double xx[npoint],yy[npoint];
  for(int j=0;j<npoint;j++)
    {
      xx[j]=par[j+1];
      yy[j]=par[j+1+npoint];
    }
  for(int j=0;j<npoint-1;j++)
    {
      
      if(x[0]>=par[j+1] && x[0]<par[j+2])
	{
	  double xa=par[j+1];
	  double xb=par[j+2];
	  double ya=par[j+1+npoint];
	  double yb=par[j+2+npoint];

	  out=x[0]*(ya-yb)/(xa-xb)+(xa*yb-ya*xb)/(xa-xb);

  break;
	}
    }
  if(x[0]<par[1])
    {
	  double xa=par[1];
	  double xb=par[2];
	  double ya=par[1+npoint];
	  double yb=par[2+npoint];
	  out=x[0]*(ya-yb)/(xa-xb)+(xa*yb-ya*xb)/(xa-xb);
    }
  if(x[0]>=par[npoint])
    {
	  double xa=par[npoint-1];
	  double xb=par[npoint];
	  double ya=par[2*npoint-1];
	  double yb=par[2*npoint];
	  out=x[0]*(ya-yb)/(xa-xb)+(xa*yb-ya*xb)/(xa-xb);
    }

  //  cout<<npoint<<endl;
  return out;
}

void unzoom()
{
h->GetXaxis()->UnZoom();
h->GetYaxis()->UnZoom();
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom1()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax/4);
  h->GetYaxis()->SetRangeUser(0,ymax/16);
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom01()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax/8);
  h->GetYaxis()->SetRangeUser(0,ymax/20);
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom2()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax/7);
  h->GetYaxis()->SetRangeUser(0,ymax/10);
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom3()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax);
  h->GetYaxis()->SetRangeUser(ymax/10,ymax/2);
 c->cd();
 gPad->Modified();
 gPad->Update();
}

void zoom4()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax/2);
  h->GetYaxis()->SetRangeUser(ymax/16,ymax/4);
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom5()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax);
  h->GetYaxis()->SetRangeUser(ymax/2,ymax);
 c->cd();
 gPad->Modified();
 gPad->Update();
}
void zoom6()
{
  float xmax=h->GetXaxis()->GetXmax();
  float ymax=h->GetYaxis()->GetXmax();

  h->GetXaxis()->SetRangeUser(0,xmax);
  h->GetYaxis()->SetRangeUser(ymax/3.2,ymax/1.333);
 c->cd();
 gPad->Modified();
 gPad->Update();
}

void ripristina()
{

  h->SetMinimum(1);
 
 c->cd();
 gPad->Modified();
 gPad->Update();
}

// void codice()
// {

//   char riga[1000];
//   char stringa[6][9][1000];
//   for(int j=0;j<6;j++)
//     {
//       for(int i=0;i<9;i++)
// 	{
// 	  sprintf(stringa[j][i],"-");
// 	}
//     }

//    int ip,ih;
//   char v1[2];
//   char v2[2];
//   sscanf(h->GetName(),"%1s%d%1s%d",v1,&ip,v2,&ih);
//   int ipl,ihl;
// static char nomec[500];
//  char value[500];
//  FILE *apri=fopen(Form("codici_qualita_gagb_%s-%s.txt",run1,run2),"r");
//   if(apri==0)
//     {
//   Mydialogo *input2=new Mydialogo("Codice di qualita'?","0 buono",nomec);
//   apri=fopen(Form("codici_qualita_gagb_%s-%s.txt",run1,run2),"w");
//       fprintf(apri,"%s hector->phosbox%d->phos%d->gagb|zquality %s %s C=%s\n",h->GetName(),ip,ih,run1,run2,nomec);
//       fclose(apri);
//     }
//   else
//     {
//       int iflag=0;
//       while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
// 	{
// 	  sscanf(riga,"%1s%d%1s%d",v1,&ipl,v2,&ihl);
// 	  if(ipl==ip &&ihl==ih)
// 	    {
// 	      for(int i=0;i<strlen(riga)-1;i++)
// 		{
// 		  if(riga[i]=='C'&&riga[i+1]=='=')
// 		    {
// 		      sscanf(&riga[i+2],"%[^\n]",value);
// 		      break;
// 		    }
// 		}

	 
	 
// Mydialogo *input2=new Mydialogo("Codice di qualita' (gia' presente)?",value,nomec);
//  sprintf(stringa[ip-1][ih-1],"%s hector->phosbox%d->phos%d->gagb|zquality %s %s C=%s",h->GetName(),ip,ih,run1,run2,nomec);
      
      
//  iflag=1;
// 	    }
// 	  else
// 	    {
// 	  sprintf(stringa[ipl-1][ihl-1],"%s",riga);
// 	    }
// 	}
//       fclose(apri);
//       if(iflag==0)
// 	{
//   Mydialogo *input2=new Mydialogo("Codice di qualita'?","0 buono",nomec);
//   sprintf(stringa[ip-1][ih-1],"%s hector->phosbox%d->phos%d->gagb|zquality %s %s C=%s",h->GetName(),ip,ih,run1,run2,nomec);
// 	}
//       apri=fopen(Form("codici_qualita_gagb_%s-%s.txt",run1,run2),"w");
//       for(int j=0;j<6;j++)
// 	{
// 	  for(int i=0;i<9;i++)
// 	    {
// 	      if(strcmp(stringa[j][i],"-")!=0)
// 		{
// 	      fprintf(apri,"%s\n",stringa[j][i]);
// 		}
// 	    }
// 	}
//       fclose(apri);
//     }
// }

void formato()
{
// static char nomec[500];
//   Mydialogo *input2=new Mydialogo("Nome del file delle linee?",Form("%s",nomefil),nomec);
//   FILE *apri=fopen(nomec,"r");
//   if(apri==0)
//     {
//       cout<<nomefil<<" non esiste"<<endl;
//       return;
//     }
char nomec[500];
 sprintf(nomec,"%s_%s_lineez.dat",Nfil,fapri);
 FILE *apri=fopen(nomec,"r");
 int flaglinee=0;

  if(apri==0)
    {
      cout<<nomec<<" non esiste"<<endl;
      flaglinee=1;
      //      return;
    }
  int out=system(Form("ls %s_%s_contorno*.dat>oloc",Nfil,fapri));
  int ncont=0; 
  char listacont[500][500];
  if(out==0)
    {

      FILE *aloc=fopen("oloc","r");
      while(fscanf(aloc,"%s",listacont[ncont])!=EOF)
	{
	  cout<<listacont[ncont]<<endl;
	  ncont++;
	}
      fclose(aloc);
      system("rm -f oloc");
      cout<<"Ci sono anche "<<ncont<<" contorni"<<endl;
    }


  FILE *aclicca=fopen(Form("DBaseToTxt_clicca_zlines_%s_%s.txt",Nfil,fapri),"w");
  


  char riga[500];
  int valori[500],npunti[500];
  char pi[10];
  int pival;
  float xx[500][500];
  float yy[500][500];
  int numlinee=0;
  for(int j=0;j<500;j++)
    {
      valori[j]=0;
      npunti[j]=0;
      for(int k=0;k<500;k++)
	{
	  xx[j][k]=0;
	  yy[j][k]=0;
	}
    }
  int iflag=0;
  int ip;
  if(flaglinee==0)
    {
  fscanf(apri,"\n%[^\n]",riga);

  while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
    {
      
      sscanf(riga,"%s",pi);
      //      cout<<riga<<" "<<pi<<endl;
      if(strcmp(pi,"PI")==0)
	{
	  sscanf(riga,"%s %d",pi,&pival);
	  //  cout<<"pi="<<pi<<" "<<pival<<endl;
	  valori[numlinee]=pival;

	  
	  numlinee++;
	}
      else
	{
	  sscanf(riga,"%d %f %f",&ip,&xx[numlinee-1][npunti[numlinee-1]],&yy[numlinee-1][npunti[numlinee-1]]);

	  npunti[numlinee-1]++;
	}

    }
  fclose(apri);
    }//flaglinee==0
  cout<<"numlinee="<<numlinee<<endl;

      if(ncont>0)
	{
	  for(int kk=0;kk<ncont;kk++)
	    {
	      FILE *apricont=fopen(listacont[kk],"r");
	      cout<<"apro "<<listacont[kk]<<endl;
  while(fscanf(apricont,"\n%[^\n]",riga)!=EOF)
    {
      TString ns=riga;

      if(ns.Contains("PI"))
	{
	  int kl=ns.Index("=");
	  sscanf(&riga[kl+1],"%d",&pival);
	  //  cout<<"pi="<<pi<<" "<<pival<<endl;
	  valori[numlinee]=pival;

	  
	  numlinee++;
	}
      else
	{
	  sscanf(riga,"%d %f %f",&ip,&xx[numlinee-1][npunti[numlinee-1]],&yy[numlinee-1][npunti[numlinee-1]]);

	  npunti[numlinee-1]++;
	}

    }
  fclose(apricont);
	    }//kk
	}//ncont>0

  int ntp=0;
  //  numlinee=2;
  cout<<numlinee<<endl;
  for(int j=0;j<numlinee;j++)
    {
      //     ntp=ntp+npunti[j];
     ntp=ntp+npunti[j]*2;
    }
  ntp=ntp+2*numlinee+2;
  float enumlinee=numlinee-ncont;
  fprintf(aclicca,"%s %s %s %d %e ",ftit,run1,run2,ntp,enumlinee);
  for(int j=0;j<numlinee-ncont;j++)
    {
            float eval=valori[j];
      //      float eval=valori[j]/10.;
      float enpunti=npunti[j];
      fprintf(aclicca,"%e %e ",eval,enpunti);
      for(int k=0;k<npunti[j];k++)
	{
	  fprintf(aclicca,"%e %e ",xx[j][k],yy[j][k]);
	}
    }//j=0,numlinee-ncont

  fprintf(aclicca,"%e ",(float)ncont);
  if(ncont>0)
    {
      for(int j=numlinee-ncont;j<numlinee;j++)
	{
      float eval=valori[j];
      float enpunti=npunti[j];
      fprintf(aclicca,"%e %e ",eval,enpunti);
      for(int k=0;k<npunti[j];k++)
	{
	  fprintf(aclicca,"%e %e ",xx[j][k],yy[j][k]);
	}

	    }//j
    }//ncont>0
  fprintf(aclicca,"\n");
  fclose(aclicca);
  system(Form("ls DBaseToTxt_clicca_zlines_%s_*.txt>o",Nfil));
  FILE *oo=fopen("o","r");
  char afil[500];
  TString unione;
  unione="cat ";
  while(fscanf(oo,"%s",afil)!=EOF)
    {
      unione.Append(Form("%s ",afil));
    }
  fclose(oo);
  cout<<unione.Data()<<endl;
  unione.Append(Form(" > DBaseToTxt_clicca_zlines_%s.txt",Nfil));
  system(unione.Data());
  cout<<"creo il file "<<Form("DBaseToTxt_clicca_zlines_%s.txt",Nfil)<<endl;
}

double iper(double *x,double *par)
{
  double out=par[0]/(x[0]+par[1]);
  return out;
}
double retta(double *x,double *par)
{
  double out=par[0]+x[0]*par[1];
  return out;
}
void punto_iniziale()
{
    for(int j=0;j<nlinee;j++)
  
    {
      if(linea[j]!=0)
	{
	  if(linea[j]->GetN()>=2)
	    {
	  double xa,ya,xb,yb;
	    linea[j]->GetPoint(0,xa,ya);
	    linea[j]->GetPoint(1,xb,yb);

	    if(xa>h->GetXaxis()->GetXmin())
	      {

	    double delta=(ya-yb)/(xa-xb);
	    double offset=(xa*yb-ya*xb)/(xa-xb);
	    
	   
		
		double xn=h->GetXaxis()->GetXmin();
		double yn=delta*xn+offset;
		int np=linea[j]->GetN();
		linea[j]->Set(np+1);
		linea[j]->SetPoint(np,xn,yn);
	
	    gPad->Modified();
	    gPad->Update();
		    ordina(j);
	
	      }
	      }
	
	}

	    }
    cout<<"ho finito"<<endl;
}



void recupera()
{
  
  //static char nomec[500];
  // Mydialogo *input2=new Mydialogo("Nome del file delle linee?",Form("%s",nomefil),nomec);
  

  FILE *aclicca=fopen(Form("DBaseToTxt_clicca_zlines_%s.txt",Nfil),"r");
  if(aclicca==0)
    {
      cout<<"Il Dbase.txt"<<" "<<Form("DBaseToTxt_clicca_zlines_%s.txt",Nfil)<<" non esiste"<<endl;
      return;
    }
 char Riga[1000000];
  char Nome[500];
  TString Tit;
  
  int isec,istrip,icsi;
  FILE *Out;
  while(fscanf(aclicca,"\n%[^\n]",Riga)!=EOF)
    {
      sscanf(Riga,"%s",Nome);
      Tit=Nome;
      if(Tit.Contains("->")||Tit.Contains("|"))
	   {
	     
 Tit.ReplaceAll(" ","");
  Tit.ReplaceAll("event->","");
  Tit.ReplaceAll("->value","");
  Tit.ReplaceAll("camere","ca");
  Tit.ReplaceAll("spicchio_","");
  Tit.ReplaceAll("ringco","rco");
  Tit.ReplaceAll("settore","set");
  Tit.ToLower();
  Tit.ReplaceAll("->","_");  
 
	
	
      if(Tit.Contains("|"))
	{
	  TString sfapri=fapri;
	  int jh=sfapri.Last('_');
	  TString pezzofin=sfapri(jh+1,sfapri.Length());
	  cout<<pezzofin.Data()<<endl;
      int jk=Tit.Index("|");
     
      Tit.Replace(jk,1,"_");
       int jlast=Tit.Last('_');
       Tit.Replace(jlast+1,pezzofin.Length(),pezzofin.Data());
	}
	
	 }
      Tit.Append(".txt");
      cout<<Tit.Data()<<endl;

	 //      sscanf(&Nome[Tit.Index("strip")+strlen("strip")],"%d",&istrip);
	 // sscanf(&Nome[Tit.Index("settore")+strlen("settore")],"%d",&isec);
	 // sscanf(&Nome[Tit.Index("csi")+strlen("csi")],"%d",&icsi);
      Out=fopen(Form("DBaseToTxt_clicca_zlines_%s_%s",Nfil,Tit.Data()),"w");
	 //      Out=fopen(Form("DBaseToTxt_clicca_zlines_%s_rco_set%d_strip%d_csi%d_sivscsi_devse.txt",Nfil,isec,istrip,icsi),"w");
      fprintf(Out,"%s\n",Riga);
      fclose(Out);

      //      cout<<isec<<endl;
    }
fclose(aclicca);

 aclicca=fopen(Form("DBaseToTxt_clicca_zlines_%s_%s.txt",Nfil,fapri),"r");
  char riga[1000000];
  int spazi[1000000];
  char nomef[100];
  int ra,rb,ilung;  
  float aval,npp,nl;
  int valori;
  int ic=0;
  int n0,n1;

   while(fscanf(aclicca,"\n%[^\n]",riga)!=EOF)
    {
      for(int nn=0;nn<strlen(riga)-1;nn++)
	{
	  if(riga[nn]==' '&&riga[nn+1]!=' ')
	    {
	      spazi[ic]=nn;
	      ic++;
	    }
	}

      sscanf(&riga[spazi[3]+1],"%f",&nl);
      nlinee=nl;
      // if(nlinee==0)
      //	{
      //	  cout<<"Non ci sono linee di Z"<<endl;
      //	  return;
      //	}
      //cout<<nlinee<<endl;
      if(nlinee>0)
	{
      int ic0=4;
      for(int k=0;k<nlinee;k++)
	{
	  sscanf(&riga[spazi[ic0]+1],"%f %f",&aval,&npp);
	  //  aval=aval*10.;
	  valori=TMath::Nint(aval);

	  linea[k]=new TGraph();
	  linea[k]->SetName(Form("PI %d",valori));
	  linea[k]->SetMarkerStyle(20);
	  linea[k]->SetMarkerColor(2);
	  linea[k]->SetMarkerSize(0.5);
	 
	  int inpp=TMath::Nint(npp);
	  //cout<<"valori="<<valori<<" "<<aval<<" "<<inpp<<endl;
	  float xa,ya;
	  for(int i=0;i<inpp;i++)
	    {
	      sscanf(&riga[spazi[ic0+2+2*i]+1],"%f %f",&xa,&ya);
	      linea[k]->SetPoint(linea[k]->GetN(),xa,ya);
	      //      if(k==0){cout<<xa<<" "<<ya<<endl;}
	    }
	  ic0=ic0+2+2*inpp;
	  
	}
sscanf(&riga[spazi[ic0]+1],"%f",&nl);

      ncont=nl;
      cout<<"ncont="<<ncont<<" "<<spazi[ic0]+1<<" "<<ic0<<endl;
      ic0=ic0+1;
      for(int k=0;k<ncont;k++)
	{
	  sscanf(&riga[spazi[ic0]+1],"%f %f",&aval,&npp);
 	  //aval=aval*10;
 	  valori=TMath::Nint(aval);
	  cout<<aval<<" "<<npp<<endl;
 	  gcont[k]=new TGraph();
 	  gcont[k]->SetName(Form("%s%d.dat",nomefilc,valori));
 	  gcont[k]->SetMarkerStyle(20);
 	  gcont[k]->SetMarkerColor(3);
 	  gcont[k]->SetLineColor(3);
	  //	  sscanf(riga,"%f",&npp);
	  int inpp=TMath::Nint(npp);
	  //cout<<"valori="<<valori<<" "<<aval<<" "<<inpp<<endl;
 	  float xa,ya;
 	  for(int i=0;i<inpp;i++)
 	    {
 	      sscanf(&riga[spazi[ic0+2+2*i]+1],"%f %f",&xa,&ya);
 	      gcont[k]->SetPoint(gcont[k]->GetN(),xa,ya);
 	      //      if(k==0){cout<<xa<<" "<<ya<<endl;}
 	    }
	  ic0=ic0+2+2*inpp;
	  
	}
	}//nlinee >0
      else
	{
	  cout<<"Non ci sono linee"<<endl;
     sscanf(&riga[spazi[4]+1],"%f",&nl);
     ncont=nl;
     int ic0=5;
      for(int k=0;k<ncont;k++)
	{
	  sscanf(&riga[spazi[ic0]+1],"%f %f",&aval,&npp);
 	  //aval=aval*10;
 	  valori=TMath::Nint(aval);

 	  gcont[k]=new TGraph();
 	  gcont[k]->SetName(Form("%s%d.dat",nomefilc,valori));
 	  gcont[k]->SetMarkerStyle(20);
 	  gcont[k]->SetMarkerColor(3);
 	  gcont[k]->SetLineColor(3);

	  //sscanf(riga,"%f",&npp);
	  int inpp=TMath::Nint(npp);
	  //cout<<"valori="<<valori<<" "<<aval<<" "<<inpp<<endl;
 	  float xa,ya;
 	  for(int i=0;i<inpp;i++)
 	    {
 	      sscanf(&riga[spazi[ic0+2+2*i]+1],"%f %f",&xa,&ya);
 	      gcont[k]->SetPoint(gcont[k]->GetN(),xa,ya);
 	      //      if(k==0){cout<<xa<<" "<<ya<<endl;}
 	    }
	  ic0=ic0+2+2*inpp;
	  
	}
	}//nlinee=0
    }
   fclose(aclicca);
   for(int j=0;j<nlinee;j++)
     {
       linea[j]->Draw("pl");
     }
   gPad->Modified();
   gPad->Update();
  for(int j=0;j<ncont;j++)
     {
       gcont[j]->Draw("pl");
     }
   gPad->Modified();
   gPad->Update();
}

void stirax()
{
  if(nlinee>0)
    {
  c->cd();
 cout<<"+ per allungare in x; - per restringere in x; s per allungare in y; t per accorciare in y"<<endl;
  cout<<"> per shiftare a dx (in x); < per shiftare a sx (in x); h per shiftare in su; b per shiftare in giu"<<endl; 
  ica=0;
c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
	      "stiratorex(Int_t,Int_t,Int_t,TObject*)");//connette il canvas a un evento grafico (non chiarissimo come funziona)
 while(ica!=1 && ica!=2 &&ica!=3)
   {usleep(100);
   gClient->HandleInput();//fondamentale, se no non funziona

   }

c->Disconnect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)");
}
}
void stiratorex(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  //   cout<<event<<" "<<gPad->GetEventX()<<endl;
  c->cd();
  if(event==24 && gPad->GetEventX()==43)//+
    {
      cout<<"allungo di 0.01 in x "<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		  xx[k]=(xx[k]-xx[0])*1.005;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
  if(event==24 && gPad->GetEventX()==45)//-
    {
      cout<<"accorcio di 0.01 in x"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		  xx[k]=(xx[k]-xx[0])/1.005;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
  if(event==24 && gPad->GetEventX()==62)//>
    {
      cout<<"shifto di 5 a dx"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		 
		  xx[k]=xx[k]+5.;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
 if(event==24 && gPad->GetEventX()==60)//<
    {
      cout<<"shifto di 5 a sx"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		  
		  xx[k]=xx[k]-5.;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
 if(event==24 && gPad->GetEventX()==104)//h
    {
      cout<<"shifto di 5 in su"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		 
		  yy[k]=yy[k]+5.;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
 if(event==24 && gPad->GetEventX()==98)//b
    {
      cout<<"shifto di 5 in giu"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		 
		  yy[k]=yy[k]-5.;
		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }

 if(event==24 && gPad->GetEventX()==115)//s
    {
      cout<<"allungo di 0.01 su y"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		  yy[k]=yy[k]*1.005;

		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }
  if(event==24 && gPad->GetEventX()==116)//t
    {
      cout<<"accorcio di 0.01 in y"<<endl;
      for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      double *xx=linea[j]->GetX();
	      double *yy=linea[j]->GetY();

	      for(int k=0;k<linea[j]->GetN();k++)
		{
		  yy[k]=yy[k]/1.005;

		  linea[j]->SetPoint(k,xx[k],yy[k]);
		  
		}
	    }
	}
      gPad->Modified();
      gPad->Update();
    }

 if(event==24 && gPad->GetEventX()==99)//c
    {
      cout<<"Esco da stiratore"<<endl;
      ica=2;
    }

}

void iperbole()
{
    for(int j=0;j<nlinee;j++)
  
    {
      if(linea[j]!=0)
	{
	  if(linea[j]->GetN()>=5)
	    {
	  double xa,ya,xb,yb;
	  linea[j]->GetPoint(linea[j]->GetN()-1,xb,yb);
	  if(xb<h->GetXaxis()->GetXmax())
	    {
	  linea[j]->GetPoint(linea[j]->GetN()-5,xa,ya);
 	  TF1 *fiper=new TF1("fiper",iper,xa,xb,2);
	  linea[j]->Fit(fiper,"R");
 	    double x0=xb;
 	    double y0=yb;
	    double xc,yc;
	    linea[j]->GetPoint(linea[j]->GetN()-2,xc,yc);
	    cout<<x0<<" "<<h->GetXaxis()->GetXmax()<<endl;

 	    while(x0<h->GetXaxis()->GetXmax() && y0<h->GetYaxis()->GetXmax())
 	      {
		cout<<"entro"<<endl;
 		double xn=x0+(xb-xc);
		double yn=fiper->Eval(xn);
 		
		//cout<<xn<<" "<<yn<<endl;
 		if(xn<h->GetXaxis()->GetXmax() && yn<h->GetYaxis()->GetXmax())
 		  {
 		int np=linea[j]->GetN();
 		linea[j]->Set(np+1);
 		linea[j]->SetPoint(np,xn,yn);
 		vecx[linea[j]->GetN()-1]=linea[j]->GetX()[np];
 		vecy[linea[j]->GetN()-1]=linea[j]->GetY()[np];
	
 		  }
 		x0=xn;
 		y0=yn;
 	      }
	    linea[j]->GetPoint(linea[j]->GetN()-1,xc,yc);
	    if(xc<h->GetXaxis()->GetXmax())
	      {
	    double xn=h->GetXaxis()->GetXmax();
	    double yn=fiper->Eval(xn);
 		int np=linea[j]->GetN();
 		linea[j]->Set(np+1);
 		linea[j]->SetPoint(np,xn,yn);
	      }
		gPad->Modified();
		gPad->Update();
		delete fiper;
		fiper=0;
	    }
	    }
	}
    }
    cout<<"ho finito"<<endl;
}




float scaprod(float *v1,float *v2)
{
  float out=v1[0]*v2[0]+v1[1]*v2[1];
  return out;
}
void rimuovi_punto_iniziale()
{
     for(int j=0;j<nlinee;j++)
	{
	  if(linea[j]->GetN()>0)
	    {
	      linea[j]->RemovePoint(0);
	      gPad->Modified();
	      gPad->Update();
	    }
	  
	}
}
void idzlines()
{
  cout<<"perche' funzioni bisogna che le linee siano in ordine di PI crescente! se non lo sono, passare dal ciclo scrivi - leggi -scrivi"<<endl;
static char nomec[500];
  Mydialogo *input1=new Mydialogo("Sei passato dal ciclo scrivi - leggi -scrivi?","yes",nomec);
  if(strcmp(nomec,"yes")!=0)
    {
      return;
    }

  int conta=0;
  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
	  conta++;
	}
    }
  if(conta<2)
    {
  Mydialogo *input1=new Mydialogo("Non ci sono almeno 2 linee in canna","OK?",nomec);
  return;
    }
  c->cd();
  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
	  spline[j]=new TSpline3(Form("spline_%s",linea[j]->GetName()),linea[j]);
	  spline[j]->SetTitle(Form("spline_%s",linea[j]->GetName()));
	  spline[j]->SetName(Form("spline_%s",linea[j]->GetName()));
	  spline[j]->SetLineColor(2);
	  spline[j]->SetMarkerColor(2);
	  spline[j]->SetMarkerStyle(20);
	 
	  // spline[j]->Draw("same");
	}
    }
	  gPad->Modified();
	  gPad->Update();

TCanvas *cloc=(TCanvas*) gROOT->FindObject("cloc");
 if(cloc!=0)
   {
     delete cloc;
     cloc=0;
   }

  cloc=new TCanvas("cloc","cloc",500,0,800,800);
  cloc->Draw();
  cloc->Divide(1,2);
  cloc->cd(1);
 gPad->SetLogy(kTRUE);
 cloc->cd(2);
 gPad->SetLogz(kTRUE);
  TH1D *hpi=(TH1D*)gROOT->FindObject("hpi");
  if(hpi!=0)
    {
      delete hpi;
      hpi=0;
    }
  TH2F *hid=(TH2F*)gROOT->FindObject("hid");
  if(hid!=0)
    {
      delete hid;
      hid=0;
    }
  TCutG *gcut=(TCutG*)gROOT->FindObject(Form("%s_%s_contorno999.dat",Nfil,fapri));
  if(gcut!=0)
    {
      delete gcut;
      gcut=0;
    }


  FILE *ap=fopen(Form("%s_%s_contorno999.dat",Nfil,fapri),"r");
		 if(ap==0)
		   {
		     cout<<"non esiste il contorno per gli stoppati"<<endl;
		   }
		 else
		   {
		     gcut=new TCutG();
		     gcut->SetName(Form("%s_%s_contorno999.dat",Nfil,fapri));
		     gcut->SetLineColor(2);
		     char rr[400];
		     fscanf(ap,"\n%[^\n]",rr);
		     int bnp;
		     float bnx,bny;
		     while(fscanf(ap,"%d %f %f",&bnp,&bnx,&bny)!=EOF)
		       {
			 cout<<bnp<<" "<<bnx<<" "<<bny<<endl;
			 gcut->Set(bnp+1);
			 gcut->SetPoint(bnp,bnx,bny);
		       }
		     fclose(ap);
		   }
		 float pi=-1;
  hpi=new TH1D("hpi","hpi",2800,0,280);
 hid=new TH2F("hid","hid",h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),h->GetYaxis()->GetNbins(),h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
 for(int j=0;j<h->GetXaxis()->GetNbins();j++)
    {
	  float gb=h->GetXaxis()->GetBinCenter(j+1);
      for(int k=0;k<h->GetYaxis()->GetNbins();k++)
	{
	  if(h->GetBinContent(j+1,k+1)>0)
	    {
	  float ga=h->GetYaxis()->GetBinCenter(k+1);
	  if((gcut!=0&& gcut->IsInside(gb,ga)==0) || gcut==0)
	    {

	  pi=pidzlines(gb,ga);
	  if(pi>0)
	    {
	      hpi->Fill(pi,h->GetBinContent(j+1,k+1));
	      hid->Fill(gb,ga,h->GetBinContent(j+1,k+1));
	    }
	    }
	    }
	}
    }
 cloc->cd(1);

  hpi->Draw();
  gPad->Modified();
  gPad->Update();
 cloc->cd(2);

  hid->Draw("zcol");
  gPad->Modified();
  gPad->Update();
  for(int j=0;j<nlinee;j++)
    {
       if(spline[j]!=0)
 	{
 	  spline[j]->Draw("same");
 	}
       //      if(linea[j]!=0)
       //	{
       //	  linea[j]->Draw("pl");
       //	}
    }
  gPad->Modified();
  gPad->Update();

  c->cd();
}
float pidzlines(float x,float y)
{
  float pi=-1;
  float y1=calcola_val(0,x);
  float y2=calcola_val(1,x);
  float pi1,pi2;
  const char * nome1;
  const char *nome2;
      nome1=linea[0]->GetName();
      nome2=linea[1]->GetName();
      sscanf(&nome1[3],"%f",&pi1);
      sscanf(&nome2[3],"%f",&pi2);

   if(y1>y && y2>y)//sotto la prima linea
    {

      pi=(y-y1)/(y2-y1)*(pi2-pi1)+pi1;
      
     
      return pi;
    }
  //tra 2 linee
  int tra_linee=0;
  float d1,d2;
  int j0=-1;
  for(int j=0;j<nlinee-1;j++)
    {
      if(linea[j]!=0 && linea[j+1]!=0)
	{
	  y1=calcola_val(j,x);
	  y2=calcola_val(j+1,x);
	  d1=y1-y;
	  d2=y2-y;
	  if(d1<0 && d2>0)
	    {
	      j0=j;
	      tra_linee=1;
	      break;
	    }
	}
    }
    if(j0>-1)
  
    {
            nome1=linea[j0]->GetName();
      nome2=linea[j0+1]->GetName();
      sscanf(&nome1[3],"%f",&pi1);
      sscanf(&nome2[3],"%f",&pi2);
     pi=(y-y1)/(y2-y1)*(pi2-pi1)+pi1;
  
      return pi;
    }
  //sopra l'ultima linea
  y1=calcola_val(nlinee-2,x);
  y2=calcola_val(nlinee-1,x);
            nome1=linea[nlinee-2]->GetName();
      nome2=linea[nlinee-1]->GetName();
      sscanf(&nome1[3],"%f",&pi1);
      sscanf(&nome2[3],"%f",&pi2);
      if((y-y1)>0 && (y-y2)>0)
	{
     pi=(y-y1)/(y2-y1)*(pi2-pi1)+pi1;
     return pi;
	}
  return pi;
}
 float calcola_val(int jlinea,float x)
 {
   float val=-1;
   double xlinea0,ylinea0;
   linea[jlinea]->GetPoint(0,xlinea0,ylinea0);
   if(x<xlinea0)
     {
       val=ylinea0+(x-xlinea0)*spline[jlinea]->Derivative(xlinea0);
       return val;
     }
   linea[jlinea]->GetPoint(linea[jlinea]->GetN()-1,xlinea0,ylinea0);
   if(x>xlinea0)
     {
 val=ylinea0+(x-xlinea0)*spline[jlinea]->Derivative(xlinea0);
       return val;
     }
val=spline[jlinea]->Eval(x);
   return val;
 }
void contorno()
{
  char  ico[100];
  sprintf(ico,"999");
  static char nomec[100];

  Mydialogo *input1=new Mydialogo("codice contorno (1/999)?",Form("%s",ico),nomec);
  sscanf(nomec,"%s",icod);

  //  Mydialogo *input1=new Mydialogo("nome del contorno?",Form("%s_contorno%d.dat",h->GetName(),icod),nomec);
  sprintf(nomec,"%s%s.dat",nomefilc,icod);
  cout<<nomec<<endl;
  if(system(Form("dir %s",nomec))==0)
    {
      static char answ[10];
Mydialogo *input2=new Mydialogo("Il file esiste gia', continuare?","yes",answ);
 if(strcmp(answ,"no")==0)
   {
     return;
   }
    }
  int ij0=-1;
  for(int j=0;j<ncont;j++)
    {
      if(strcmp(nomec,gcont[j]->GetName())==0)
	{
	  ij0=j;
	  break;
	}
    }
  if(ij0>-1)
    {
      cout<<"contorno gia' caricato; procedere con le modifiche"<<endl;
      indicec=ij0;
    }
  else
    {

  cout<<"Mettere il contorno "<<nomec<<endl;
  indicec=ncont;
    }
  cout<<"ncont="<<ncont<<endl;
  ica=0;
  if(gcont[indicec]==0)
    {
      gcont[indicec]=new TGraph();
      gcont[indicec]->SetMarkerStyle(20);
      gcont[indicec]->SetMarkerColor(3);
      gcont[indicec]->SetLineColor(3);
      gcont[indicec]->SetName(Form("%s",nomec));

    }

  c->cd();
c->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
	      "mettic(Int_t,Int_t,Int_t,TObject*)");//connette il canvas a un evento grafico (non chiarissimo come funziona)
 while(ica!=1 && ica!=2 &&ica!=3)
   {usleep(100);
   gClient->HandleInput();//fondamentale, se no non funziona
   }

c->Disconnect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)");

 if(ica!=3)
   {
  cout<<"ho messo il contorno "<<nomec<<endl;
  FILE *apri=fopen(nomec,"w");
  double *xxj;
  double *yyj;

  xxj=gcont[indicec]->GetX();
  yyj=gcont[indicec]->GetY();


  fprintf(apri,"# PI=%s\n",icod);

  for(int j=0;j<gcont[indicec]->GetN();j++)
    {

printf("%d %f %f\n",j,xxj[j],yyj[j]);      
fprintf(apri,"%d %f %f\n",j,xxj[j],yyj[j]);
      
    }
  fclose(apri);
  TFile *froot=new TFile(Form("%s.root",nomec),"RECREATE");
  froot->cd();
  TGraph *gr;
  gr=(TGraph*)gcont[indicec]->Clone();
  gr->SetName(Form("contorno%s",icod));
  gr->Write();
  froot->Write();
  froot->Close();
  delete gr;
  gr=0;

  cout<<"ho salvato il contorno "<<gcont[indicec]->GetName()<<endl;
  if(ij0<0){  ncont++;}
   }


}
void scrivicontorni()
{
  for(int j=0;j<ncont;j++)
    {
      cout<<"salvo il contorno "<<gcont[j]->GetName()<<endl;
      FILE *apri=fopen(gcont[j]->GetName(),"w");
      double *xxj;
      double *yyj;
  xxj=gcont[j]->GetX();
  yyj=gcont[j]->GetY();

  char stringa[400];
  sprintf(stringa,"%s",gcont[j]->GetName());
  for(int nn=0;nn<strlen(stringa)-8;nn++)
    {
     
      if(stringa[nn]=='c'&&stringa[nn+1]=='o'&&stringa[nn+2]=='n' && stringa[nn+3]=='t' && stringa[nn+4]=='o'&&stringa[nn+5]=='r'&&stringa[nn+6]=='n'&& stringa[nn+7]=='o')
	{
	  	  sscanf(&stringa[nn+8],"%s",icod);
	  break;
	}
    }

  
  cout<<"icod="<<icod<<endl; 

  fprintf(apri,"# PI=%s\n",icod);

  for(int k=0;k<gcont[j]->GetN();k++)
    {

printf("%d %f %f\n",k,xxj[k],yyj[k]);      
fprintf(apri,"%d %f %f\n",k,xxj[k],yyj[k]);
      
    }
  fclose(apri);
    }


  TGraph *gr;
  for(int j=0;j<ncont;j++)
    {
      TFile *froot=new TFile(Form("%s.root",gcont[j]->GetName()),"RECREATE");
  froot->cd(); 
      gr=(TGraph*)gcont[j]->Clone();
      TString sl=gcont[j]->GetName();
      int jk=sl.Index("contorno");
      int icd;
      sscanf(&gcont[j]->GetName()[jk+8],"%d",&icd);
      
      gr->SetName(Form("contorno%d",icd));
      gr->Write();
  froot->Write();
  froot->Close();
      delete gr;
      gr=0;
    }


}
void leggicontorno()
{
  static char nomec[100];
int  ico=999;
  Mydialogo *input1=new Mydialogo("codice contorno (1/999)?",Form("%d",ico),nomec);
 sscanf(nomec,"%s",icod);
  //  Mydialogo *input1=new Mydialogo("nome del contorno?",Form("%s_contorno%d.dat",h->GetName(),ncont),nomec);
sprintf(nomec,"%s%s.dat",nomefilc,icod);

  cout<<nomec<<endl;

  FILE *apri=fopen(nomec,"r");
  if(apri==0)
    {
      cout<<"il file "<<nomec<<" non esiste"<<endl;
      return;
    }
  int iflag0=-1;
int ij0=-1;
  for(int j=0;j<ncont;j++)
    {
      if(strcmp(nomec,gcont[j]->GetName())==0)
	{
	  cout<<"questo contorno e' gia' caricato"<<endl;
	  ij0=j;
	  iflag0=1;
	  break;

	}
    }
  if(iflag0==1)
    {
 static char nomec33[100];
  Mydialogo *input33=new Mydialogo("Attenzione, contorno gia' caricato; lo cancello e lo ricarico?","yes",nomec33);
  if(strcmp(nomec33,"yes")!=0)
    {
      return;
    }
  else
    {
      delete gcont[ij0];
      gcont[ij0]=0;
      gPad->Modified();
      gPad->Update();
 
      //      gPad->GetListOfPrimitives()->FindObject(nomec)->Delete();
 
      //gPad->Modified();
      //gPad->Update();
 
      for(int k=ij0;k<ncont-1;k++)
	{
	  gcont[k]=gcont[k+1];
	}
      gcont[ncont-1]=0;
 
      ncont--;
 
    }
    }

 
  gcont[ncont]=new TGraph();
      gcont[ncont]->SetMarkerStyle(20);
 	  gcont[ncont]->SetMarkerColor(3);
 	  gcont[ncont]->SetLineColor(3);
      gcont[ncont]->SetName(Form("%s",nomec));
  int np;
  float x,y;
  char riga[400];
  fscanf(apri,"\n%[^\n]",riga);

  while(fscanf(apri,"%d %f %f",&np,&x,&y)!=EOF)
    {
      gcont[ncont]->Set(np+1);
      gcont[ncont]->SetPoint(np,x,y);
      
    }
  gcont[ncont]->Draw("pl");
  gPad->Modified();
  gPad->Update();
  ncont++;
}

void leggitutticontorni()
{

  static char nomec[100];

 Mydialogo *input1=new Mydialogo("Nome files contorni?",Form("%s*.dat",nomefilc),nomec);


  cout<<nomec<<endl;
  system(Form("ls %s >oo",nomec));
  FILE *ooo=fopen("oo","r");
  char nomecont[200];
  FILE *apri=0;
  while(fscanf(ooo,"%s",&nomecont)!=EOF)
    {
      int iflag0=-1;
      int ij0=-1;
      for(int j=0;j<ncont;j++)
	{
	  cout<<nomecont<<" "<<gcont[j]->GetName()<<endl;
       if(strcmp(nomecont,gcont[j]->GetName())==0)
	{
	  cout<<"questo contorno e' gia' caricato"<<endl;
	  ij0=j;
	  iflag0=1;
	  break;
	}	  
	}
  if(iflag0==1)
    {
 static char nomec33[100];

 Mydialogo *input33=new Mydialogo(Form("contorno %s gia' caricato. Cancello e ricarico?",nomecont),"yes",nomec33);
  if(strcmp(nomec33,"yes")!=0)
    {
      return;
    }
  else
    {
      delete gcont[ij0];
      gcont[ij0]=0;
      gPad->Modified();
      gPad->Update();
     
      //  gPad->GetListOfPrimitives()->FindObject(nomecont)->Delete();
      //gPad->Modified();
      //gPad->Update();
      for(int k=ij0;k<ncont-1;k++)
	{
	  gcont[k]=gcont[k+1];
	}
      gcont[ncont-1]=0;
      ncont--;
    }
    }

  gcont[ncont]=new TGraph();
      gcont[ncont]->SetMarkerStyle(20);
 	  gcont[ncont]->SetMarkerColor(6);
 	  gcont[ncont]->SetLineColor(6);
	   gcont[ncont]->SetMarkerSize(0.6);
      gcont[ncont]->SetName(Form("%s",nomecont));
  int np;
  float x,y;
  char riga[400];
  apri=fopen(nomecont,"r");
  fscanf(apri,"\n%[^\n]",riga);

  while(fscanf(apri,"%d %f %f",&np,&x,&y)!=EOF)
    {
      gcont[ncont]->Set(np+1);
      gcont[ncont]->SetPoint(np,x,y);
      
    }
  gcont[ncont]->Draw("pl");
  gPad->Modified();
  gPad->Update();
  ncont++;
  fclose(apri);
      
    }//while fscanf
  fclose(ooo);


}

void togli_contorno()
{
 static char nomec[100];
int  ico=999;
  Mydialogo *input1=new Mydialogo("codice contorno (1/999)?",Form("%d",ico),nomec);
 sscanf(nomec,"%s",icod);
sprintf(nomec,"%s%s.dat",nomefilc,icod);

  cout<<nomec<<endl;
  int ij0=-1;
  for(int j=0;j<ncont;j++)
    {
      if(strcmp(nomec,gcont[j]->GetName())==0)
	{
	  ij0=j;
	  break;
	}
    }
  if(ij0==-1)
    {
      cout<<"questo contorno non e' caricato"<<endl;
      return;
    }
      delete gcont[ij0];
      gcont[ij0]=0;
      gPad->Modified();
      gPad->Update();
      for(int k=ij0;k<ncont-1;k++)
	{
	  gcont[k]=gcont[k+1];
	}
      gcont[ncont-1]=0;
      ncont--;
    
return;
}


void mettic(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  //  cout<<event<<endl;
  c->cd();
  if(event==1)
    {
float px=gPad->AbsPixeltoX(x);
     float py=gPad->AbsPixeltoY(y);
     py=gPad->PadtoY(py);
     float uymin=gPad->GetUymin();
     float uymax=gPad->GetUymax();
     if(px>=gPad->GetUxmin() && px<=gPad->GetUxmax() && py>=gPad->PadtoY(uymin) && py<=gPad->PadtoY(uymax))
       {
	 int n=gcont[indicec]->GetN();
	 gcont[indicec]->Set(n+1);
	 gcont[indicec]->SetPoint(n,px,py);
	 if(gcont[indicec]->GetN()==1)
	   {
	     gcont[indicec]->Draw("pl");

	   }

	 //	 gPad->Modified();
	 //gPad->Update();
       }
       }
  if(event==24 && gPad->GetEventX()==99)//c per uscire
    {
  int np=gcont[indicec]->GetN();

  gcont[indicec]->Set(np+1);

  double xl,yl;
  gcont[indicec]->GetPoint(0,xl,yl);
  gcont[indicec]->SetPoint(np,xl,yl);
  cout<<"fine contorno "<<gcont[indicec]->GetName()<<endl;
      ica=2;
    }

  if(event==24 && gPad->GetEventX()==117)// u per togliere l'ultimo punto
    {
      if(gcont[indicec]->GetN()>1)
	{
	  gcont[indicec]->RemovePoint(gcont[indicec]->GetN()-1);
	  gPad->Modified();
	  gPad->Update();
	  cout<<"tolgo l'ultimo punto messo"<<endl;
	  

	}
    }

  if(event==24 && gPad->GetEventX()==120)// x per togliere il primo punto
    {
      if(gcont[indicec]->GetN()>1)
	{
	  gcont[indicec]->RemovePoint(0);
	  gPad->Modified();
	  gPad->Update();
	  cout<<"tolgo il primo punto messo"<<endl;
	  

	}
    }
 if(event==24 && gPad->GetEventX()==100) // d per cancellare contorno
    {
      cout<<"Cancello il contorno "<<gcont[indicec]->GetName()<<endl;
      delete gcont[indicec];
      gcont[indicec]=0;
      gPad->Modified();
      gPad->Update();
      for(int k=indicec;k<ncont-1;k++)
	{
	  gcont[k]=gcont[k+1];
	}
      gcont[ncont-1]=0;
      ncont--;
      ica=3;
    }
 if(event==24 && gPad->GetEventX()==105)//i inverto l'ordine dei punti 
   {
     if(gcont[indicec]->GetN()>1)
       {
	 cout<<"Inverto l'ordine dei punti"<<endl;
	 double xc[500];
	 double yc[500];
	 for(int k=0;k<gcont[indicec]->GetN();k++)
	   {double a,b;
	     gcont[indicec]->GetPoint(k,a,b);
	     xc[k]=a;
	     yc[k]=b;
	     
	   }

	 for(int k=0;k<gcont[indicec]->GetN();k++)
	   {
	     gcont[indicec]->SetPoint(k,xc[gcont[indicec]->GetN()-k-1],yc[gcont[indicec]->GetN()-k-1]);
	   }

	 gPad->Modified();
	 gPad->Update();

       }
   }
}

void intervallo()
{ 
static char nomec[500];
 Mydialogo *input2=new Mydialogo("Intervallo di validita'? run1,run2","000000000000,999999999999",nomec);
 cout<<nomec<<endl;
 
   sscanf(nomec,"%12s,%12s",run1,run2);
 cout<<"run1="<<run1<<endl;
 cout<<"run2="<<run2<<endl;

}

void recuperaroot()
{
 static char nomec[100];

 Mydialogo *input1=new Mydialogo("Nome files contorni?",Form("%s*.dat.root",nomefilc),nomec);

  cout<<nomec<<endl;
  system(Form("ls %s >oo",nomec));
  FILE *ooo=fopen("oo","r");
  char nomecont[200];

  while(fscanf(ooo,"%s",&nomecont)!=EOF)
    {
      TString nct=nomecont;
      int jk=nct.Index(".root");
      TString nct2=nct(0,jk);
      sprintf(nomecont,"%s",nct2.Data());
    
      int iflag0=-1;
      int ij0=-1;
      for(int j=0;j<ncont;j++)
	{
	  cout<<nomecont<<" "<<gcont[j]->GetName()<<endl;
       if(strcmp(nomecont,gcont[j]->GetName())==0)
	{
	  cout<<"questo contorno e' gia' caricato"<<endl;
	  ij0=j;
	  iflag0=1;
	  break;
	}	  
	}
  if(iflag0==1)
    {
 static char nomec33[100];

 Mydialogo *input33=new Mydialogo(Form("contorno %s gia' caricato. Cancello e ricarico?",nomecont),"yes",nomec33);
  if(strcmp(nomec33,"yes")!=0)
    {
      return;
    }
  else
    {
      delete gcont[ij0];
      gcont[ij0]=0;
      gPad->Modified();
      gPad->Update();
     
      //  gPad->GetListOfPrimitives()->FindObject(nomecont)->Delete();
      //gPad->Modified();
      //gPad->Update();
      for(int k=ij0;k<ncont-1;k++)
	{
	  gcont[k]=gcont[k+1];
	}
      gcont[ncont-1]=0;
      ncont--;
    }
    }

  TFile *ffapri=new TFile(Form("%s.root",nomecont));
  TIter next(ffapri->GetListOfKeys());
  TKey *key;

  while(key=(TKey*)next())
    {
      
      cout<<key->GetName()<<" "<<key->GetClassName()<<endl;
      
      if(strcmp(key->GetClassName(),"TGraph")==0 ||strcmp(key->GetClassName(),"TCutG")==0)
	{
	  cout<<"qui"<<endl;
	  gcont[ncont]=(TGraph*)key->ReadObj();
	  cout<<"dopo"<<" "<<gcont[ncont]<<endl;


      gcont[ncont]->SetMarkerStyle(20);
 	  gcont[ncont]->SetMarkerColor(6);
 	  gcont[ncont]->SetLineColor(6);
	   gcont[ncont]->SetMarkerSize(0.6);
	   // gcont[ncont]->SetName(Form("%s",nomecont));
gcont[ncont]->Draw("pl");
  gPad->Modified();
  gPad->Update();
  ncont++;
	} 
    }
      
    }//while fscanf
  fclose(ooo);
 static char nomecx[500];
  if(nlinee>0)
    {
  Mydialogo *input11=new Mydialogo("Attenzione, cancello tutte le linee presenti; continuare?","yes",nomecx);
  if(strcmp(nomecx,"yes")!=0)
    {
      return;
    }     

  for(int j=0;j<nlinee;j++)
    {
      if(linea[j]!=0)
	{
      delete linea[j];
      linea[j]=0;
	}
      if(funz[j]!=0)
	{
	  delete funz[j];
	  funz[j]=0;
	}
    }
  nlinee=0;
    }

  Mydialogo *input21=new Mydialogo("Nome del file delle linee?",Form("%s.root",nomefil),nomec);
  FILE *apri=fopen(nomec,"r");
  if(apri==0)
    {
      cout<<Form("%s.root",nomefil)<<" non esiste"<<endl;
      return;
    }
  fclose(apri);

 
  int valori[500];
  char pi[10];
  int pival;

  for(int j=0;j<500;j++)
    {
      valori[j]=0;

    }
  TFile *ffroot=new TFile(Form("%s.root",nomefil));
  TIter next(gFile->GetListOfKeys());
  TObject *obj;
  nlinee=0;
  TGraph *grlinee[500];
  char lll[20];
  while(obj=(TObject*)next())
    {cout<<obj->IsA()->GetName()<<endl;
      if(obj->InheritsFrom("TKey"))
	{
	  cout<<obj->GetName()<<endl;
	  grlinee[nlinee]=(TGraph*)ffroot->Get(obj->GetName());
	  int k0=-1;
	  for(int k=0;k<strlen(grlinee[nlinee]->GetName());k++)
	    {
	      if(grlinee[nlinee]->GetName()[k]>='0'&&grlinee[nlinee]->GetName()[k]<='9')
		{
		  k0=k;
		  break;
		}
	    }
	  sscanf(&grlinee[nlinee]->GetName()[k0],"%d",&pival);

	  //sscanf(grlinee[nlinee]->GetName(),"%s%d",lll,&pival);
	  cout<<pival<<endl;
	  valori[nlinee]=pival;
	  nlinee++;
	}
    }
  int colonna[500];
  for(int j=0;j<500;j++)
    {
      colonna[j]=j;
    }
  for(int j=0;j<nlinee-1;j++)
    {
      for(int k=j+1;k<nlinee;k++)
	{
	  if(valori[j]>valori[k])
	    {

	      int li=colonna[k];
	      colonna[k]=colonna[j];
	      colonna[j]=li;
	      	      float serv=valori[k];
	       valori[k]=valori[j];
	      valori[j]=serv;
	    }
	}
    }  
  c->cd();
  for(int j=0;j<nlinee;j++)
    {
      linea[j]=(TGraph*)grlinee[colonna[j]]->Clone();
      linea[j]->SetName(Form("PI %d",valori[j]));
     cout<<linea[j]->GetName()<<endl;
      linea[j]->SetMarkerStyle(20);
      linea[j]->SetMarkerColor(2);
      linea[j]->SetMarkerSize(0.5);
      linea[j]->Draw("pl");
    }
  gPad->Modified();
  gPad->Update();

  for(int j=0;j<nlinee;j++)
    {
      delete grlinee[j];
     
    }
  for(int j=0;j<500;j++)
    {
      grlinee[j]=0;
    }
  ffroot->Close();

}


 void recuperakali()
 {
    static char nomec[100];


  static char nomecx[500];
   if(nlinee>0 || ncont>0)
     {
   Mydialogo *input11=new Mydialogo("Attenzione, cancello tutte le linee  e i contorni presenti; continuare?","yes",nomecx);
   if(strcmp(nomecx,"yes")!=0)
     {
       return;
     }     

   for(int j=0;j<nlinee;j++)
     {
       if(linea[j]!=0)
 	{
       delete linea[j];
       linea[j]=0;
 	}
       if(funz[j]!=0)
 	{
 	  delete funz[j];
 	  funz[j]=0;
 	}
     }
   nlinee=0;

   for(int j=0;j<ncont;j++)
     {
       if(gcont[j]!=0)
	 {
	   delete gcont[j];
	   gcont[j]=0;
      gPad->Modified();
      gPad->Update();	   
	 }
     }
   ncont=0;
     }
   char riga[1000];
  Mydialogo *input21=new Mydialogo("Nome del file delle linee?",Form("%s.dat",histoname),nomec);
   FILE *apri=fopen(nomec,"r");
   if(apri==0)
     {
       cout<<Form("%s",nomec)<<" non esiste"<<endl;
       return;
     }
   else
     {
       int inew=0;
       int npunti=0;
       while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
	 {
	  
	   TString sriga=riga;
	   cout<<riga<<endl;
	   if(sriga.Contains("+KVIDZALine"))
	     {
	       cout<<"parte una nuova riga"<<endl;
	       linea[nlinee]=new TGraph();
	       linea[nlinee]->SetMarkerStyle(20);
	       linea[nlinee]->SetMarkerColor(2);
	       linea[nlinee]->SetMarkerSize(0.5);
	       nlinee++;
	       
	       npunti=0;
	       fscanf(apri,"\n%[^\n]",riga);

		  fscanf(apri,"%d %d",&zval[nlinee-1],&aval[nlinee-1]);
		  linea[nlinee-1]->SetName(Form("PI %d Z=%d A=%d",zval[nlinee-1]*100+aval[nlinee-1]-zval[nlinee-1],zval[nlinee-1],aval[nlinee-1]));
		  fscanf(apri,"%d",&npunti);

		  float a,b;
		  for(int k=0;k<npunti;k++)
		    {
		  fscanf(apri,"%f %f",&a,&b);
		  linea[nlinee-1]->SetPoint(linea[nlinee-1]->GetN(),a,b);
		}
		  linea[nlinee-1]->Draw("pl");
		  gPad->Modified();
		  gPad->Update();
		}
	   if(sriga.Contains("+KVIDCutContour")||sriga.Contains("+KVIDCutLine"))
	     {
	       cout<<"parte un nuovo contorno"<<endl;
	       gcont[ncont]=new TGraph();
	       gcont[ncont]->SetMarkerStyle(20);
	       gcont[ncont]->SetMarkerColor(6);
	       gcont[ncont]->SetLineColor(6);
	       gcont[ncont]->SetMarkerSize(0.6);
	       ncont++;
	       fscanf(apri,"\n%[^\n]",riga);
	       char nome[100];
	       sscanf(&riga[3],"%s",nome);
	       gcont[ncont-1]->SetName(Form("%s%s.dat",nomefilc,nome));
	       fscanf(apri,"\n%[^\n]",riga);
	       int npt=0;
	       fscanf(apri,"%d",&npt);
	       float a,b;
	       for(int j=0;j<npt;j++)
		 {
		   fscanf(apri,"%f %f",&a,&b);
		   gcont[ncont-1]->SetPoint(gcont[ncont-1]->GetN(),a,b);
		 }
	       
	       
	   gcont[ncont-1]->Draw("pl");
	   gPad->Modified();
	   gPad->Update();
	     }
	 }
   fclose(apri);
     }
 

 }

void scrivikali()
{
static char nomec[100];
    if(nlinee==0 && ncont==0)
    {
      cout<<"Non ci sono linee o contorni da salvare"<<endl;
      return;
    }
  Mydialogo *input2=new Mydialogo("Nome del file delle linee e contorni?",Form("%s.dat",histoname),nomec);
  if(system(Form("dir %s",nomec))==0)
    {
  static char nomec2[500];
       Mydialogo *input2=new Mydialogo("Il file esiste gia'. Sovrascrivo?","yes",nomec2);
       if(strcmp(nomec2,"yes")!=0)
	 {
return;
	 }
       system(Form("cp %s %s.old",nomec,nomec));//backup!
    }

  FILE *apri=fopen(nomec,"w");
  fprintf(apri,"# ASCII file generated by KVIDZAGrid::WriteToAsciiFile\n");
  TString sh=histoname;
  char nomet[1000];
  TString si="SI";
  TString tele="T";
    TString qua="Q";
    TString blo="B";
    int nsil,ntele,nqua,nblo,ntelesc;
 int k0=-1;
 for(int k=strlen(histoname)-1;k>=0;k--)
   {
     if(histoname[k]=='_')
       {
	 k0=k;
	 break;
       }

   }

 if(k0>=0)
   {
     sscanf(&histoname[k0+1],"%d",&ntelesc);
   }
 cout<<"ntele="<<ntelesc<<endl;
 nblo=ntelesc/100;
 nqua=(ntelesc-nblo*100)/10;
 ntele=ntelesc-nblo*100-nqua*10;



 //  sscanf(&sh.Data()[sh.Index(si.Data())+si.Length()],"%d",&nsil);
 // sscanf(&sh.Data()[sh.Index(tele.Data())+tele.Length()],"%d",&ntele);
 // sscanf(&sh.Data()[sh.Index(qua.Data())+qua.Length()],"%d",&nqua);
 // sscanf(&sh.Data()[sh.Index(blo.Data())+blo.Length()],"%d",&nblo);
  sprintf(nomet,"%s",histoname);
  //sprintf(nomet,"ID_%s%d-%s%d-%s%d-%s%03d",si.Data(),nsil,tele.Data(),ntele,qua.Data(),nqua,blo.Data(),nblo);
 fprintf(apri,"# ID Graph Name : %s\n",nomet);
 fprintf(apri,"# This file can be read using KVIDZAGrid::ReadFromAsciiFile\n");
 fprintf(apri,"#\n");
 fprintf(apri,"++KVIDZAGrid\n");
 fprintf(apri,"<VARX>\n");
 fprintf(apri,"<VARY>\n");
 fprintf(apri,"<PARAMETER> IDTelescopes=%s\n",nomet);
 fprintf(apri,"<PARAMETER> Runlist=\n");
 fprintf(apri,"OnlyZId 0\n");

 for(int j=0;j<nlinee;j++)
   {
     cout<<"scrivo la linea "<<j<<" "<<linea[j]->GetName()<<endl;
     fprintf(apri,"+KVIDZALine\n");
     fprintf(apri,"ID:Z=%d A=%d\n",zval[j],aval[j]);
     fprintf(apri,"%d\t%d\n",zval[j],aval[j]);
     fprintf(apri,"%d\n",linea[j]->GetN());
     double *xp;
     double *yp;
     xp=linea[j]->GetX();
     yp=linea[j]->GetY();
     for(int k=0;k<linea[j]->GetN();k++)
       {
	 fprintf(apri,"%f   %f\n",xp[k],yp[k]);
       }
   }
 for(int j=0;j<ncont;j++)
   {
     cout<<"scrivo il contorno "<<j<<" "<<gcont[j]->GetName()<<endl;
     fprintf(apri,"+KVIDCutContour\n");
     TString sc=gcont[j]->GetName();
     TString scr="contorno";
     int ikk=sc.Index(scr.Data())+scr.Length();
     int ik2=sc.Index(".dat");
     TString sc2=sc(ikk,ik2-ikk);
     //cout<<sc2.Data()<<" "<<scr.Data()<<" "<<ikk<<" "<<" "<<ik2<<" "<<sc.Data()<<endl;
     fprintf(apri,"OK:%s\n",sc2.Data());
     fprintf(apri,"0\n");
     fprintf(apri,"%d\n",gcont[j]->GetN());

     double *xp;
     double *yp;
     xp=gcont[j]->GetX();
     yp=gcont[j]->GetY();
     for(int k=0;k<gcont[j]->GetN();k++)
       {
	 fprintf(apri,"%f   %f\n",xp[k],yp[k]);
       }
   }

 fprintf(apri,"!\n");
  fclose(apri);
}

void duplica()
{
  if(nlinee>=1)
    {
	   TString sl=linea[nlinee-1]->GetName();
	   int pival;
	   sscanf(&sl.Data()[sl.Index("PI")+2],"%d",&pival);
  static char nomec[500];
  Mydialogo *input2=new Mydialogo("PI linea da duplicare?",Form("PI %d",pival),nomec);


  TString sc=nomec;
  int pival2,pic;
sscanf(&sc.Data()[sc.Index("PI")+2],"%d",&pic);
  indice=-1;

   cout<<"Duplicare la linea "<<nomec<<endl;
   for(int j=0;j<nlinee;j++)
     {
       if(linea[j]!=0)
	 {
	   TString sl=linea[j]->GetName();
	   sscanf(&sl.Data()[sl.Index("PI")+2],"%d",&pival2);
	   if(pival2==pic)
	 {
	   indice=j;
	   break;
	 }
	 }
     }
   cout<<"Indice linea da duplicare="<<indice<<endl;
   if(indice>=0)
     {
       Mydialogo *input3=new Mydialogo("PI linea duplicata?",Form("PI %d Z=%d A=%d",(nlinee+1)*10,nlinee+1,2*(nlinee+1)),nomec);
  TString sc2=nomec;

sscanf(&sc2.Data()[sc2.Index("PI")+2],"%d",&pic);

      linea[nlinee]=new TGraph();
      linea[nlinee]->SetMarkerStyle(20);
      linea[nlinee]->SetMarkerColor(2);
      linea[nlinee]->SetMarkerSize(0.5);
      linea[nlinee]->SetName(Form("%s",nomec));
    
    zval[nlinee]=nlinee+1;
      aval[nlinee]=2*zval[nlinee];

      double *xp=linea[indice]->GetX();
      double *yp=linea[indice]->GetY();
      cout<<indice<<endl;
      double xx,yy;
      if(indice>0)
	{
	  funz[indice-1]=new TF1(Form("f %s",linea[indice-1]->GetName()),spezzata,0,4000,2*linea[indice-1]->GetN()+1);
	  cout<<funz[indice-1]->GetName()<<endl;


	 int a=linea[indice-1]->GetN();
	  funz[indice-1]->SetParameter(0,a);
	  for(int k=0;k<linea[indice-1]->GetN();k++)
	    {
	      linea[indice-1]->GetPoint(k,xx,yy);
	      funz[indice-1]->SetParameter(k+1,xx);
	      funz[indice-1]->SetParameter(k+1+linea[indice-1]->GetN(),yy);
	    }	  


	  //  spline[indice-1]=new TSpline3(Form("spline_%s",linea[indice-1]->GetName()),linea[indice-1]);
	  //spline[indice-1]->SetTitle(Form("spline_%s",linea[indice-1]->GetName()));
	  //spline[indice-1]->SetName(Form("spline_%s",linea[indice-1]->GetName()));

	  for(int j=0;j<linea[indice]->GetN();j++)
	    {
	      //	      linea[nlinee]->SetPoint(linea[nlinee]->GetN(),xp[j],yp[j]+yp[j]-spline[indice-1]->Eval(xp[j]));
  linea[nlinee]->SetPoint(linea[nlinee]->GetN(),xp[j],yp[j]+yp[j]-funz[indice-1]->Eval(xp[j]));
	    }

	}
      else
	{
  Mydialogo *input4=new Mydialogo("shift in y?",Form("%d",100),nomec);
  int shift;
  sscanf(nomec,"%d",&shift);
  cout<<"shift="<<shift<<endl;
	  for(int j=0;j<linea[indice]->GetN();j++)
	    {
	      linea[nlinee]->SetPoint(linea[nlinee]->GetN(),xp[j],yp[j]+shift);
	    }
	}

      linea[nlinee]->Draw("pl");
  
 nlinee++;
cout<<"ho messo la linea "<<nomec<<endl;
cout<<"nlinee="<<nlinee<<endl;
 gPad->Modified();
 gPad->Update();
     }
   else
     {
       cout<<"Linea da duplicare non esistente"<<endl;
       return;
     }

    }
  else
    {
      cout<<"non ci sono linee da duplicare"<<endl;
    }
  return;
}
