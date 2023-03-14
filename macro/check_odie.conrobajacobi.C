#include <stdio.h>
#include <string.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TObject.h>
#include <TCollection.h>
#include <signal.h>
#include <TNamed.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>

using namespace std;

TFile *fout;
TTree *tree;
void mysignal(int a)
{

printf("segnale%d control c Scrivo il file di uscita e chiudo\n",a);
 tree->Write();
  fout->Write();
  fout->Close();



exit(1);

 
//    return;

}
int main(int argv,char *argc[])
{

  signal(2,mysignal);
  printf("%d\n",argv);
  if(argv<2)
    {
      cout<<"Manca il file da leggere"<<endl;
      return 0;
    }

  TString slist[50];
      for(int k=0;k<argv;k++)
	{
	
	  slist[k]=argc[k];
	  cout<<slist[k].Data()<<" "<<slist[k].IsDigit()<<endl;
	  
	}
 int neventi=10;
 int flaglast=0;
      if(slist[argv-1].IsDigit())
	{
	  neventi=atoi(slist[argv-1].Data());
	  flaglast=1;
	}
      TString smeno;
     
      if(slist[argv-1][0]=='-')
	{
	
	  smeno=slist[argv-1](1,slist[argv-1].Length());
	  if(smeno.IsDigit())
	    {
	      flaglast=1;
	      neventi=-1;
	    }
	}
      for(int k=1;k<argv-flaglast;k++)
	{
	  //	  cout<<"list file="<<argc[k]<<endl;

TFile *apri=new TFile(argc[k]);
		  if(apri->IsZombie()==1)

 	    {
 	      cout<<"Il file "<<argc[k]<<" non esiste"<<endl;
 	      return 0;
 	    }
		  apri->Close();
	}
 //  TFile *f=new TFile(argc[1]);
  TChain *ch=new TChain("odietree","odietree");
  TFile *f[20];

  for(int j=1;j<argv-flaglast;j++)
    {
      f[j-1]=new TFile(argc[j]);
      ch->Add(f[j-1]->GetName());
      cout<<"Aggiungo i files "<<f[j-1]->GetName()<<endl;
    }

  //ch->Add(f->GetName());

  cout<<"neventi="<<ch->GetEntries()<<endl;
  char str[2000];
  TString ss2;
  for(int j=1;j<argv-flaglast;j++)
    {
      TString ss=f[j-1]->GetName();
      TString ss1=ss(0,ss.Index(".root"));
		     if(j==1)
		       {
			 ss2=ss1;
		       }
		     else
		       {
			 ss2.Append("_");
			 ss2.Append(ss1.Data());
		       }
    }
  ss2.Append("_check.root");

  fout=new TFile(Form("%s",ss2.Data()),"RECREATE");

  cout<<"Sto creando il file di test "<<fout->GetName()<<endl;
  
  TNamed *name[20];
  
  for(int j=1;j<argv-flaglast;j++)
    {
      name[j-1]=new TNamed(Form("%s",f[j-1]->GetName()),Form("listafiles_%d",j-1));
      name[j-1]->Write();
    }  




  
  int molt,trigger;
  float z[500],a[500];
  int coderiv[500];
  float tempo[500],e[500];
  float pi[500],up[500],down[500],Ecsi[500];
  int icsi[500],isec[500],lato[500];
  //char run[13];
  Long64_t run;
  tree=new TTree("odie","odie");
   tree->Branch("molt",&molt,"molt/i");
   tree->Branch("isec",isec,"isec[molt]/i");
   tree->Branch("icsi",icsi,"icsi[molt]/i");
   tree->Branch("lato",lato,"lato[molt]/i");
   
   //  tree->Branch("z",z,"z[molt]/F");
//   tree->Branch("a",a,"a[molt]/F");
//   tree->Branch("tempo",tempo,"tempo[molt]/F");
//   tree->Branch("e",e,"e[molt]/F");
//   tree->Branch("coderiv",coderiv,"coderiv[molt]/i");
//   tree->Branch("trigger",&trigger,"trigger/i");
   tree->Branch("pi",pi,"pi[molt]/F");
   tree->Branch("up",up,"up[molt]/F");
   tree->Branch("down",down,"down[molt]/F");
   tree->Branch("Ecsi",Ecsi,"Ecsi[molt]/F");
   //   tree->Branch("run",run,"run/C");
   //   tree->Branch("run",&run,"run/Long64_t");
   tree->Branch("run",&run,"run/Long64_t");
   if(neventi<0)
    {
      neventi=ch->GetEntries();
    }
  cout<<"Eventi da leggere="<<neventi<<endl; 


  UInt_t value_N;//n. di valori nell'evento
  Float_t value_val[1000];
  int value_worker_id[1000];
  int value_worker_class_code[1000];
  char run_num[200];
	      ch->SetBranchAddress("value_N",&value_N);
	      ch->SetBranchAddress("value_val",value_val);
	      ch->SetBranchAddress("value_worker_id",value_worker_id);
	      ch->SetBranchAddress("value_worker_class_code",value_worker_class_code);
	      ch->SetBranchAddress("run_num",run_num);



  for(int iev=0;iev<neventi;iev++)
    {
      ch->GetEntry(iev);
      molt=0;
      //      sprintf(run,"%s",run_num);
      sscanf(run_num,"%lld",&run);
      for(int j=0;j<value_N;j++)
	{
	  //  cout<<value_N<<endl;
	  if(value_worker_class_code[j]==14688)
	    {
	      isec[molt]=value_worker_id[j]/10;
	      icsi[molt]=value_worker_id[j]-isec[molt]*10;
	      up[molt]=value_val[j];
	      j++;
	      down[molt]=value_val[j];
	      j++;
	      Ecsi[molt]=value_val[j];
	      j++;
	      lato[molt]=value_val[j];
	      j++;
	      z[molt]=value_val[j];
	      j++;
	      pi[molt]=value_val[j];
	      molt++;
	    }//14688

// 	  if(value_worker_class_code[j]==11584)
// 	    {
// 	      trigger=value_val[j];


// 	    }
// 	  if(value_worker_class_code[j]==18400)//plastichino: 9999
// 	    {
// 	      coderiv[molt]=9999;
// 	      tempo[molt]=value_val[j];
// 	      e[molt]=-1;
// 	      molt++;
// 	    }

// 	  if(value_worker_class_code[j]==15328)//phoswich: ip*10+ih (ip e ih da 1) (11-69)
// 	    {
// 	      int ip=value_worker_id[j]/1000;
// 	      int ih=(value_worker_id[j]-ip*1000)/100;
// 	      z[molt]=value_val[j];
// 	      j++;
// 	      a[molt]=value_val[j];
// 	      j++;
// 	      // expevent.phos_qf[ip-1][ih-1]=value_val[j];
// 	      j++;
// 	      // expevent.phos_cod[ip-1][ih-1]=value_val[j];
// 	      j++;
// 	      tempo[molt]=value_val[j];

// 	      //*********** DA TOGLIERE SE NON c'E' ga nell'Ntupla
// 			      j++;
// 	      //expevent.phos_ga[ip-1][ih-1]=value_val[j];
// 	      //*********** FINE
// 		  e[molt]=-1;
// 		  coderiv[molt]=ip*10+ih;
// 		  molt++;  
// 		  }//classe 15328
// 	  if(value_worker_class_code[j]==24880)//garfield: isec*10+icsi+1000 (ip e ih da 1) (1011-1248)
// 	    {
// 	      int isec=value_worker_id[j]/10;
// 	      int icsi=value_worker_id[j]-isec*10;
	      
// 	      z[molt]=value_val[j];
// 	      j++;
// 	      a[molt]=value_val[j];
// 	      j++;
// 	      //expevent.garf_qf[isec-1][icsi-1]=value_val[j];
// 	      j++;
// 	      //expevent.garf_theta[isec-1][icsi-1]=value_val[j];
// 	      j++;
// 	      //expevent.garf_phi[isec-1][icsi-1]=value_val[j];
// 	      j++;
// 	      e[molt]=value_val[j];
// 	      tempo[molt]=-1;
// 	      coderiv[molt]=isec*10+icsi+1000;
// 	      molt++;
// 	    }//classe 24880
// 	  if(value_worker_class_code[j]==3728)//BaF2 ibaf+10000 (ibaf da 1) (10001-10008)
//    {
// 	      int ibaf=value_worker_id[j]/10;

	      
// 	      e[molt]=value_val[j];
// 	      j++;
// 	      tempo[molt]=value_val[j];
// 	      coderiv[molt]=ibaf+10000;

	 
	 
// 	      molt++;
//    }//classe 3728
// 	  //<<<<<<< Classe_evento.cxx
// 	  if(value_worker_class_code[j]==24416)//RCO isec*100+istrip*10+icsi+1000000 (isec,istrip,icsi da 1) (1000111-1000897)
// 	    {
// 	      int settore=value_worker_id[j]/1000-1;
// 	      z[molt]=value_val[j];
// 	      j++;
// 	      a[molt]=value_val[j];
// 	      j++;
// 	      //	      Double_t theloc=value_val[j];
// 	      j++;
// 	       Double_t philoc=value_val[j];
// 	      j++;	  
// 	      // Double_t qfloc=value_val[j];
// 	      j++;
// 	       Double_t codeloc=value_val[j];
// 	      j++;

// 	      e[molt]=value_val[j];
// 	      int istrip0,icsi0;
// 	      istrip0=(value_worker_id[j]-1000*(settore+1))/100-1;

// 	      if(TMath::Nint(codeloc)==1 || TMath::Nint(codeloc)==2)//identificazione da camera-Si o da PSA in Si
// 		{

// 		  icsi0=6;//cesio non definito
// 		}
// 	      if(TMath::Nint(codeloc)==3)//identificazione da Si-Csi
// 		{

// 		  icsi0=philoc-(settore+1)*10-1;
// 		}
// 	      tempo[molt]=-1;
// 	      coderiv[molt]=(settore+1)*100+(istrip0+1)*10+icsi0+1+1000000;

// 		  molt++;

// 	    }//classe 24420
//   if(value_worker_class_code[j]==17344)//gamma e neutroni RCO isec*100+istrip*10+icsi+1000000 (isec,istrip,icsi da 1) (1000191-1000896)
//  	    {
//  	      int settore=value_worker_id[j]/1000-1;
// 	      z[molt]=value_val[j];

//  	      j++;

// 	      a[molt]=value_val[j];
//  	      j++;
// // 	      Double_t theloc=value_val[j];
//  	      j++;
// // 	      Double_t philoc=value_val[j];
//  	      j++;	  
// 	      //  	      Double_t qfloc=value_val[j];
//  	      j++;
//  	      Double_t codeloc=value_val[j];
//  	  j++;
// 	  e[molt]=value_val[j];
// 	  tempo[molt]=-1;

//  	      int istrip0,icsi0;
//  	      icsi0=(value_worker_id[j]-1000*(settore+1))/10-1;

//  	      if(TMath::Nint(codeloc)==4)//identificazione da Csi Fast -Slow
//  		{
//  		  istrip0=8;//strip non definita

//  		}

// 	 coderiv[molt]=(settore+1)*100+(istrip0+1)*10+icsi0+1+1000000;	

// 	  molt++;

// 	    }//17344
    

    }
  tree->Fill();
    }
  tree->Write();
  fout->Write();
  fout->Close();
  return 0;
}
