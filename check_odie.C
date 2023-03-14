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
  float z[500],a[500],tplast;
  int coderiv[500];
  float e[500];
  int isec[500],istrip[500],icsi[500];
  int rco_code[500];
  int value[500];
  Long64_t run;
  tree=new TTree("odie","odie");
   tree->Branch("molt",&molt,"molt/i");
   tree->Branch("z",z,"z[molt]/F");
   tree->Branch("a",a,"a[molt]/F");
  tree->Branch("e",e,"e[molt]/F");
  //   tree->Branch("coderiv",coderiv,"coderiv[molt]/i");
  tree->Branch("isec",isec,"isec[molt]/I");
  tree->Branch("istrip",istrip,"istrip[molt]/I");
  tree->Branch("icsi",icsi,"icsi[molt]/I");
  tree->Branch("rco_code",rco_code,"rco_code[molt]/I");
   tree->Branch("trigger",&trigger,"trigger/I");
   tree->Branch("tplast",&tplast,"tplast/F");
   tree->Branch("run",&run,"run/Long64_t");
   tree->Branch("value_worker_class_code",value,"value[molt]/I");

   //per fare micro-csi di garfield
   // tree->Branch("lato",lato,"lato[molt]/i");
   // tree->Branch("pi",pi,"pi[molt]/F");
   // tree->Branch("up",up,"up[molt]/F");
   // tree->Branch("down",down,"down[molt]/F");
   // tree->Branch("Ecsi",Ecsi,"Ecsi[molt]/F")

   if(neventi<0)
    {
      neventi=ch->GetEntries();
    }
  cout<<"Eventi da leggere="<<neventi<<endl; 


  UInt_t value_N;//n. di valori nell'evento
  Float_t value_val[1000];
  UInt_t value_worker_id[1000];
  UInt_t value_worker_class_code[1000];
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
      tplast=-1;
      //      sprintf(run,"%s",run_num);
      sscanf(run_num,"%lld",&run);
      for(int j=0;j<value_N;j++)
	{
	  //  cout<<value_N<<endl;
	  if(value_worker_class_code[j]==11584)
	    {
	      trigger=value_val[j];
	      //  printf("%d trig %f\n",trig,expcomp.value_val[j]);
	    }
	  if(value_worker_class_code[j]==18400)
	    {
	      tplast=value_val[j];
	    }

	  if(value_worker_class_code[j]==24880)//GARFIELD
	    {
	      value[molt]=24880;
	      isec[molt]=value_worker_id[j]/10;
	      int isec0=isec[molt];
	      icsi[molt]=value_worker_id[j]-isec0*10;
	      istrip[molt]=-1;
	      z[molt]=value_val[j];
	      j++;
	      a[molt]=value_val[j];
	      j++;
	      // garf_qf[isec-1][icsi-1]=value_val[j];
	      j++;
	      //garf_theta[isec-1][icsi-1]=value_val[j];
	      j++;
	      //garf_phi[isec-1][icsi-1]=value_val[j];
	      j++;
	      e[molt]=value_val[j];

	      if(TMath::Nint(z[molt])==50 && TMath::Nint(a[molt])==50)
		{
		  z[molt]=0;
		  a[molt]=0;
		}

	      if(z[molt]>0 && a[molt]<0)
		{
		  if(TMath::Nint(z[molt])==1)
		    {
		      a[molt]=1; 
		    }
		  else
		    {
		      a[molt]=2* z[molt];
		    }
		}

	      rco_code[molt]=-1;
	      molt++;
	    }//classe 24880

	  //<<<<<<< Classe_evento.cxx
	  if(value_worker_class_code[j]==24416)//RCO 
	    {
	      value[molt]=24416;
	      //+1000000
	      int settore=value_worker_id[j]/1000;
	      isec[molt]=settore;
	      z[molt]=value_val[j];
	      j++;
	      a[molt]=value_val[j];
	      j++;
	      //Double_t theloc=value_val[j];
	      j++;
	      double philoc=value_val[j];
	      j++;	  
 	      //qfloc=value_val[j];
	      j++;
	        Double_t codeloc=value_val[j];
	      j++;
	      
	      e[molt]=value_val[j];
	      //cout<<settore<<" "<<zloc<<" "<<aloc<<" "<<theloc<<" "<<philoc<<" "<<qfloc<<" "<<codeloc<<endl;
	      int istrip0,icsi0=0;
	      istrip0=(value_worker_id[j]-1000*(settore))/100;
	      istrip[molt]=istrip0;
	      if(TMath::Nint(codeloc)==1 || TMath::Nint(codeloc)==2)//identificazione da camera-Si o da PSA in Si
		{

		  icsi0=7;//cesio non definito
		}
	      if(TMath::Nint(codeloc)==3)//identificazione da Si-Csi
		{

		  icsi0=philoc-(settore)*10;
		}
	      icsi[molt]=icsi0;
	      rco_code[molt]=codeloc;
	      molt++;
		


	    }//classe 24420
  if(value_worker_class_code[j]==17344)//(solo fast-slow)
 	    {
//+1000000
	      value[molt]=17344;
 	      int settore=value_worker_id[j]/1000;
	      isec[molt]=settore;
	      z[molt]=value_val[j];
 	      j++;
 	      a[molt]=value_val[j];
 	      j++;
 	      //Double_t theloc=value_val[j];
 	      j++;
 	      //Double_t philoc=value_val[j];
 	      j++;	  
  	      //Double_t qfloc=value_val[j];
 	      j++;
 	      Double_t codeloc=value_val[j];
 	  j++;
 	  e[molt]=value_val[j];//occhio e' luce!
// 	      //cout<<settore<<" "<<zloc<<" "<<aloc<<" "<<theloc<<" "<<philoc<<" "<<qfloc<<" "<<codeloc<<endl;
 	      int istrip0,icsi0;
 	      icsi0=(value_worker_id[j]-1000*(settore))/10;
	      icsi[molt]=icsi0;
 	      if(TMath::Nint(codeloc)==4)//identificazione da Csi Fast -Slow
 		{
 		  istrip0=9;//strip non definita

 		}
	      if(TMath::Nint(z[molt])==10 && TMath::Nint(a[molt])==10)
		{
 		  z[molt]=0.;
 		  a[molt]=0.;
		}
	     		
 		  rco_code[molt]=codeloc;
 		 

		  istrip[molt]=istrip0;

// 		  //il theta e il phi si assegnano con la spalmatura
		
		  molt++;

 	    }//classe 17344	  

  //per micro-csi di Garfield
  // 	  if(value_worker_class_code[j]==14688)
  //	    {
  //	      isec[molt]=value_worker_id[j]/10;
  //	      icsi[molt]=value_worker_id[j]-isec[molt]*10;
  //	      up[molt]=value_val[j];
  //	      j++;
  //	      down[molt]=value_val[j];
  //	      j++;
  //	      Ecsi[molt]=value_val[j];
  //	      j++;
  //	      lato[molt]=value_val[j];
  //	      j++;
  //	      z[molt]=value_val[j];
  //	      j++;
  //	      pi[molt]=value_val[j];
  //	      molt++;
  //	    }//14688   

    }
  tree->Fill();
    }
  tree->Write();
  fout->Write();
  fout->Close();
  return 0;
}
