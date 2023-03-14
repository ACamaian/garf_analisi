#include <stdio.h>
#include <iostream>
#include <TCanvas.h>
#include <string.h>
#include <TH2.h>
#include <TH1.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMath.h>
float modulo(float *input)
    {
      float out=TMath::Sqrt(input[0]*input[0]+input[1]*input[1]+input[2]*input[2]);
      return out;
    }

void carpol(float *r,float *the, float *phi,float *input)
    {
      *r=modulo(input);
    
      float comp=input[2]/(*r);
      if(comp>1){comp=1.;}
      if(comp<-1){comp=-1.;}
      (*the)=57.296*TMath::ACos(comp);

  if(input[0]==0. && input[1]==0.)
	{
	  (*phi)=0.;
	}
      else
	{
	  (*phi)=57.296*TMath::ATan2(input[1],input[0]);
	}
  return;
  
    }
void polcar(float r, float the, float phi,float *out)
      {
	//printf("r=%f the=%f phi=%f\n",r,the,phi);
	if(the==0.)
	  {
	    out[2]=0.;
	  }
	else
	  {
	out[2]=r*cos(the/57.296);
	  }
	out[0]=r*sin(the/57.296)*cos(phi/57.296);
	out[1]=r*sin(the/57.296)*sin(phi/57.296);


    }

void baiocco(char *file="FI_24Mg_2.6_Jtp12.dat")
{
  TString s=file;
  
  TTree *baiocco=new TTree("baiocco","baiocco");
  FILE *apri=fopen(file,"r");
  


		UInt_t mtot;
		float spin;
		float a[100],z[100],Px[100],Py[100],Pz[100],exc[100];
		UInt_t discreto[100];
		baiocco->Branch("mtot",&mtot,"mtot/i");
		baiocco->Branch("spin",&spin,"spin/F");
		baiocco->Branch("z",z,"z[mtot]/F");	
	baiocco->Branch("a",a,"a[mtot]/F");
	baiocco->Branch("exc",exc,"exc[mtot]/F");
	baiocco->Branch("Px",Px,"Px[mtot]/F");
	baiocco->Branch("Py",Py,"Py[mtot]/F");
	baiocco->Branch("Pz",Pz,"Pz[mtot]/F");
      
		baiocco->Branch("discreto",discreto,"discreto[mtot]/i");
		char riga[400];
	  float PP[3];
  int nev=0;
  int ia,iz;
  float pthe[100],pphi[100],pmod[100];
    while(1)
  // while(nev<10)
    {

        if(fscanf(apri,"\n%[^\n]",riga)!=EOF)
      {
      //      cout<<riga<<endl;
	sscanf(riga,"%d %f",&mtot,&spin);


      for(int j=0;j<mtot;j++)//loop sull'evento
	{
	  fscanf(apri,"%d %d %f %f %f %f %d",&ia,&iz,&exc[j],&Px[j],&Py[j],&Pz[j],&discreto[j]);
	  a[j]=(float)ia;
          z[j]=(float)iz;

	  PP[0]=Px[j];
	  PP[1]=Py[j];
	  PP[2]=Pz[j];

	  carpol(&pmod[j],&pthe[j],&pphi[j],PP);
	
	}
      float phirand= 180-360*gRandom->Rndm();
      for(int j=0;j<mtot;j++)
	{
	  pphi[j]=pphi[j]+phirand;
	  if(pphi[j]>180)
	    {
	      pphi[j]=pphi[j]-360;
	    }
	  if(pphi[j]<-180)
	    {
	      pphi[j]=360+pphi[j];
	    }
	  polcar(pmod[j],pthe[j],pphi[j],PP);
	  Px[j]=PP[0];
	  Py[j]=PP[1];
	  Pz[j]=PP[2];
	}

      baiocco->Fill();
      nev++;
      //  cout<<nev<<endl;
      }
	else
      {
        break;
      }
    }
  fclose(apri);
   cout<<nev<<endl;

   TFile *fou=new TFile(Form("%s.root",file),"RECREATE");
   fou->cd();
   baiocco->Write();
   fou->Close();
}
