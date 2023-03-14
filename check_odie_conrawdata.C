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
double E2Luce_csi(double z,double a,double E);//FORMULA MICHELA
double Luce2E_csi(double z,double a,double lo); //formula Michela INVERTITA
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

  int moltrawg,moltrawr,moltrawr2;
  int isecrawg[500],istriprawr[500],icsirawg[500],isecrawr[500],icsirawr[500],isecrawr2[500],icsirawr2[500];
  float ic[500],icraw[500],sie[500],sieraw[500],sitrise[500],csislowr[500],csislowpsar[500],csiluceg[500],ecsipsar2[500],ecsi[500];
  float csislowg[500],csifastg[500],csislowpsag[500],micro[500],csifastr2[500],csislowpsar2[500],microchan[500],microchanhg[500];
  float Sie[500],Ic[500],Sitrise[500],Csislow[500],Csifast[500],Csislowpsa[500],Micro[500],Sieraw[500],Icraw[500],Microchan[500],Microchanhg[500],Csiluce[500],CsiE[500];
  //char rivelatore[500][20];
  // struct nome{
  //  char rivelatore[500][20];
  // }tipo;
  char rivelatore[10000];
  for(int j=0;j<500;j++)
    {

      isecrawg[j]=-1;
 isecrawr[j]=-1;
 isecrawr2[j]=-1;
      istriprawr[j]=-1;
      icsirawg[j]=-1;
      icsirawr2[j]=-1;
      icsirawr[j]=-1;
      ic[j]=-1;
      icraw[j]=-1;

      sie[j]=-1;
      sieraw[j]=-1;
      sitrise[j]=-1;
      csiluceg[j]=-1;
      csifastg[j]=-1;
      csifastr2[j]=-1;
      ecsipsar2[j]=-1;
      csislowpsag[j]=-1;
      csislowpsar2[j]=-1;
      csislowr[j]=-1;
      ecsi[j]=-1;
     
      micro[j]=-1;
      microchan[j]=-1;
      microchanhg[j]=-1;
      
      Sie[j]=-1;
      Sieraw[j]=-1;
      Micro[j]=-1;
      Microchan[j]=-1;
      Microchanhg[j]=-1;
      Ic[j]=-1;
      Icraw[j]=-1;
      Sitrise[j]=-1;
      Csislow[j]=-1;
      Csifast[j]=-1;
      Csislowpsa[j]=-1;
      CsiE[j]=-1;
      Csiluce[j]=-1;
      //      tipo.rivelatore[j][0]='\0';

    }
  rivelatore[0]='\0';

  int indexrco[500];

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
   //tree->Branch("rivelatore",rivelatore,"rivelatore[molt]/C");
   //tree->Branch("rivelatore",tipo.rivelatore,"rivelatore[molt]/C");
   tree->Branch("rivelatore",rivelatore,"rivelatore/C");
   tree->Branch("Sie",Sie,"Sie[molt]/F");
 tree->Branch("Sieraw",Sieraw,"Sieraw[molt]/F");
   tree->Branch("Ic",Ic,"Ic[molt]/F");
   tree->Branch("Icraw",Icraw,"Icraw[molt]/F");
   tree->Branch("Sitrise",Sitrise,"Sitrise[molt]/F");
   tree->Branch("Csislow",Csislow,"Csislow[molt]/F");
   tree->Branch("Csifast",Csifast,"Csifast[molt]/F");
   tree->Branch("Csislowpsa",Csislowpsa,"Csislowpsa[molt]/F");
   tree->Branch("Micro",Micro,"Micro[molt]/F");
  tree->Branch("Microchan",Microchan,"Microchan[molt]/F");
 tree->Branch("Microchanhg",Microchanhg,"Microchanhg[molt]/F");
tree->Branch("CsiE",CsiE,"CsiE[molt]/F");
 tree->Branch("Csiluce",Csiluce,"Csiluce[molt]/F");
   

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
  for(int j=0;j<500;j++)
    {

      isecrawg[j]=-1;
 isecrawr[j]=-1;
 isecrawr2[j]=-1;
      istriprawr[j]=-1;
      icsirawg[j]=-1;
      icsirawr2[j]=-1;
      icsirawr[j]=-1;
      ic[j]=-1;
      icraw[j]=-1;

      sie[j]=-1;
      sieraw[j]=-1;
      sitrise[j]=-1;
      csiluceg[j]=-1;
      csifastg[j]=-1;
      csifastr2[j]=-1;
      ecsipsar2[j]=-1;
      csislowpsag[j]=-1;
      csislowpsar2[j]=-1;
      csislowr[j]=-1;
     
     
      micro[j]=-1;
      microchan[j]=-1;
microchanhg[j]=-1;

      Sie[j]=-1;
      Sieraw[j]=-1;
      Micro[j]=-1;
      Microchan[j]=-1;
      Microchanhg[j]=-1;
      Ic[j]=-1;
      Icraw[j]=-1;
      Sitrise[j]=-1;
      Csislow[j]=-1;
      Csifast[j]=-1;
      Csislowpsa[j]=-1;
      CsiE[j]=-1;
      Csiluce[j]=-1;
    }
rivelatore[0]='\0';
      ch->GetEntry(iev);
      //cout<<"iev="<<iev<<endl;
      molt=0;
      moltrawg=0;
      moltrawr=0;
      moltrawr2=0;
      tplast=-1;
      //      sprintf(run,"%s",run_num);
      sscanf(run_num,"%lld",&run);
      for(int j=0;j<value_N;j++)
	{

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
	      //istrip[molt]=-1;
	      isec[molt]=value_worker_id[j]/10;
	      int isec0=isec[molt];
	      icsi[molt]=value_worker_id[j]-isec0*10;
	      
	      z[molt]=value_val[j];
	      j++;
	      a[molt]=value_val[j];
	      j++;
	      float qf=value_val[j];
	      j++;
	      //garf_theta[isec-1][icsi-1]=value_val[j];
	      j++;
	      //garf_phi[isec-1][icsi-1]=value_val[j];
	      istrip[molt]=value_val[j]*10;//per garfield il numero di strip e' 10 per la sinistra, 20 per la destra; -10 se e' fast slow
	      j++;


	      e[molt]=value_val[j];
	      char val[20];
	

	      if(qf>=190)
		{
		  sprintf(val,"GARF Micro-CsI;");

		}
	      else
		{
		 sprintf(val,"GARF Fast Slow CsI;");
		}

	      if(strlen(rivelatore)>0)
		{
	      int kx=-1;
		  for(int k=strlen(rivelatore);k>0;k--)
		    {
		      if(rivelatore[k]==';')
			{
			  kx=k;
			  break;
			}
		    }
		  sprintf(&rivelatore[kx+1],"%s",val);
		}
	      else
		{
		  sprintf(rivelatore,"%s",val);
		}
	      // if(TMath::Nint(z[molt])==50 && TMath::Nint(a[molt])==50)
	      // 	{
	      // 	  z[molt]=0;
	      // 	  a[molt]=0;
	      // 	}

	      // if(z[molt]>0 && a[molt]<0)
	      // 	{
	      // 	  if(TMath::Nint(z[molt])==1)
	      // 	    {
	      // 	      a[molt]=1; 
	      // 	    }
	      // 	  else
	      // 	    {
	      // 	      a[molt]=2* z[molt];
	      // 	    }
	      // 	}

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
	      char val[20];

	      if(rco_code[molt]==1)
		{
		  	      sprintf(val,"RCO IC-Si;");
		}
	      if(rco_code[molt]==2)
		{
		  sprintf(val,"RCO PSA Si;");
		}
	      if(rco_code[molt]==3)
		{
		  sprintf(val,"RCO Si-CsI;");
		}
	      if(strlen(rivelatore)>0)
		{
	      int kx=-1;
		  for(int k=strlen(rivelatore);k>0;k--)
		    {
		      if(rivelatore[k]==';')
			{
			  kx=k;
			  break;
			}
		    }
		  sprintf(&rivelatore[kx+1],"%s",val);
		}
	      else
		{
		  sprintf(rivelatore,"%s",val);
		}
	      // if(TMath::Nint(a[molt]==-1))
	      // 	{
	      // 	  if(TMath::Nint(z[molt]==1))
	      // 	    {
	      // 	      a[molt]=1;
	      // 	    }
	      // 	  else
	      // 	    {
	      // 	      a[molt]=2*z[molt];
	      // 	    }
	      // 	}
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
 	  e[molt]=value_val[j];
// 	      //cout<<settore<<" "<<zloc<<" "<<aloc<<" "<<theloc<<" "<<philoc<<" "<<qfloc<<" "<<codeloc<<endl;
 	      int istrip0,icsi0;
 	      icsi0=(value_worker_id[j]-1000*(settore))/10;
	      icsi[molt]=icsi0;
 	      if(TMath::Nint(codeloc)==4)//identificazione da Csi Fast -Slow
 		{
 		  istrip0=9;//strip non definita

 		}
	      // if(TMath::Nint(z[molt])==10 && TMath::Nint(a[molt])==10)//eventuali gamma
	      // 	{
 	      // 	  z[molt]=0.;
 	      // 	  a[molt]=0.;
	      // 	}

 		  rco_code[molt]=codeloc;
		  char val[20];
 		  sprintf(val,"RCO Fast-Slow CsI;");
	      if(strlen(rivelatore)>0)
		{
	      int kx=-1;
		  for(int k=strlen(rivelatore);k>0;k--)
		    {
		      if(rivelatore[k]==';')
			{
			  kx=k;
			  break;
			}
		    }
		  sprintf(&rivelatore[kx+1],"%s",val);
		}
	      else
		{
		  sprintf(rivelatore,"%s",val);
		}

		  istrip[molt]=istrip0;

// 		  //il theta e il phi si assegnano con la spalmatura
		
		  molt++;

 	    }//classe 17344	  
    
                    if(value_worker_class_code[j]==10052) { //tvolo GARF
			int isec=value_worker_id[j]/10;
			int icsi=value_worker_id[j]%10;
			isecrawg[moltrawg]=isec;
			icsirawg[moltrawg]=icsi;
			//	expevent.garf_tvolo[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			//RAW DATA
			microchan[moltrawg]=value_val[j];
			j++;

			micro[moltrawg]=value_val[j];
			j++;
			csifastg[moltrawg]=value_val[j];
			j++;
			csiluceg[moltrawg]=value_val[j];
			j++;
			csislowpsag[moltrawg]=value_val[j];
			moltrawg++;
			


		}//classe 10052


	if(value_worker_class_code[j]==2000) { //raw GARF
			int isec=value_worker_id[j]/10;
			int icsi=value_worker_id[j]%10;	
			isecrawg[moltrawg]=isec;
			icsirawg[moltrawg]=icsi;	
			//RAW DATA
			microchan[moltrawg]=value_val[j];
			j++;
			micro[moltrawg]=value_val[j];
			j++;
			csifastg[moltrawg]=value_val[j];
			j++;
			csiluceg[moltrawg]=value_val[j];
			j++;
			csislowpsag[moltrawg]=value_val[j];
		
			
			moltrawg++;

		}//classe 2000 (sostituisce 10052)
 	if(value_worker_class_code[j]==2001) { //raw GARF
			int isec=value_worker_id[j]/10;
			int icsi=value_worker_id[j]%10;	
			isecrawg[moltrawg]=isec;
			icsirawg[moltrawg]=icsi;	
			//RAW DATA
			microchan[moltrawg]=value_val[j];
			j++;
			microchanhg[moltrawg]=value_val[j];
			j++;			
			micro[moltrawg]=value_val[j];
			j++;
			csifastg[moltrawg]=value_val[j];
			j++;
			csiluceg[moltrawg]=value_val[j];
			j++;
			csislowpsag[moltrawg]=value_val[j];
		
			
			moltrawg++;

		}//classe 2001 (sostituisce 10052 e 2000)
 	    





 	    if(value_worker_class_code[j]==10053) { //tvolo RCO strip
			int isec=value_worker_id[j]/1000;
			int isi=(value_worker_id[j]%1000)/100;
			isecrawr[moltrawr]=isec;
			istriprawr[moltrawr]=isi;
			

			//			double tvolo=value_val[j];
			j++;
			//			double tcode=value_val[j];
			j++;
			
			int icsi=(int)(value_val[j]+0.5);
		
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
				icsirawr[moltrawr]=icsi;

	//			expevent.rco_tvolo[isec-1][isi-1][icsi-1]=tvolo;
			//	expevent.rco_tcode[isec-1][isi-1][icsi-1]=tcode;//1 tempo da IC; 2 tempo da Si; 3 tempo da CsI
			//RAW DATA
			icraw[moltrawr]=value_val[j];
			j++;
			ic[moltrawr]=value_val[j];
			j++;

			sieraw[moltrawr]=value_val[j];
			j++;
			sie[moltrawr]=value_val[j];
			j++;
			sitrise[moltrawr]=value_val[j];
			j++;
			csislowr[moltrawr]=value_val[j];
			moltrawr++;
			//	tempo=1;
		}//classe 10053



	    if(value_worker_class_code[j]==1000) { //raw RCO strip
			int isec=value_worker_id[j]/1000;
			int isi=(value_worker_id[j]%1000)/100;
			isecrawr[moltrawr]=isec;
			istriprawr[moltrawr]=isi;
			int icsi=(int)(value_val[j]+0.5);
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			icsirawr[moltrawr]=icsi;
			//RAW DATA
			icraw[moltrawr]=value_val[j];
			
			j++;
			ic[moltrawr]=value_val[j];
			j++;

			sieraw[moltrawr]=value_val[j];
			j++;
			sie[moltrawr]=value_val[j];
			j++;
			sitrise[moltrawr]=value_val[j];
			j++;
			csislowr[moltrawr]=value_val[j];
			moltrawr++;
			
		}//classe 1000 (sostituisce 10053)

	    if(value_worker_class_code[j]==1001) { //raw RCO strip
			int isec=value_worker_id[j]/1000;
			int isi=(value_worker_id[j]%1000)/100;
			isecrawr[moltrawr]=isec;
			istriprawr[moltrawr]=isi;
			int icsi=(int)(value_val[j]+0.5);
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			icsirawr[moltrawr]=icsi;
			//RAW DATA
			icraw[moltrawr]=value_val[j];
			
			j++;
			ic[moltrawr]=value_val[j];
			j++;

			sieraw[moltrawr]=value_val[j];
			j++;
			sie[moltrawr]=value_val[j];
			j++;
			sitrise[moltrawr]=value_val[j];
			j++;
			csislowr[moltrawr]=value_val[j];
			j++;
			ecsi[moltrawr]=value_val[j];
			moltrawr++;
			
		}//classe 1001 (sostituisce 10053 e 1000)


		if(value_worker_class_code[j]==10054) { //tvolo RCO exit
			int isec=value_worker_id[j]/1000;
			int icsi=(value_worker_id[j]%100)/10;
			isecrawr2[moltrawr2]=isec;
			icsirawr2[moltrawr2]=icsi;
			
			//			expevent.rco_tvolo[isec-1][8][icsi-1]=value_val[j];
			j++;
			//			expevent.rco_tcode[isec-1][8][icsi-1]=3.;
			//RAW DATA
			csifastr2[moltrawr2]=value_val[j];
			j++;
			ecsipsar2[moltrawr2]=value_val[j];
			j++;
			csislowpsar2[moltrawr2]=value_val[j];
			moltrawr2++;
			//			tempo=1;
		}//classe 10054
  		if(value_worker_class_code[j]==1500) { //RCO exit raw
			int isec=value_worker_id[j]/1000;
			int icsi=(value_worker_id[j]%100)/10;
			isecrawr2[moltrawr2]=isec;
			icsirawr2[moltrawr2]=icsi;
			
			//RAW DATA
			csifastr2[moltrawr2]=value_val[j];
			j++;
			ecsipsar2[moltrawr2]=value_val[j];
			j++;
			csislowpsar2[moltrawr2]=value_val[j];
			moltrawr2++;	
		}//classe 1500 (sostituisce 10054)
 	    	

	}
      //      cout<<moltrawr<<" "<<moltrawg<<" "<<moltrawr2<<endl;

      if(moltrawr>0||moltrawr2>0 ||moltrawg>0)
	{
	  if(molt!=(moltrawr+moltrawr2+moltrawg))
	    {
	      cout<<"Molt raw diversa da molt cal"<<" "<<iev<<endl;
	    }
	  

	  for(int k=0;k<molt;k++)
	    {

	      if(value[k]==24880)//garfield
		{
		  
		  for(int i=0;i<moltrawg;i++)
		    {
		  
		      if(isecrawg[i]==isec[k]&&icsirawg[i]==icsi[k])
			{
			  Micro[k]=micro[i];
			  Microchan[k]=microchan[i];
			   Microchanhg[k]=microchanhg[i];
			  
			  Csifast[k]=csifastg[i];
			  Csislowpsa[k]=csislowpsag[i];
			  Csiluce[k]=csiluceg[i];
			  CsiE[k]=Luce2E_csi((double)z[k],(double)a[k],csiluceg[i]);

			  break;
			}

		    }

		}

	      if(value[k]==24416)//RCO
		{

		  for(int i=0;i<moltrawr;i++)
		    {

		      if(isecrawr[i]==isec[k]&&icsirawr[i]==icsi[k]&&istriprawr[i]==istrip[k])
			{

			  Icraw[k]=icraw[i];
			  
			  Ic[k]=ic[i];
			  
			  Sieraw[k]=sieraw[i];
			  Sie[k]=sie[i];
			  Sitrise[k]=sitrise[i];
			  Csislow[k]=csislowr[i];
			  CsiE[k]=ecsi[i];

			  break;
			}
		    }		  
		}

	      if(value[k]==17344)//fast slow RCO
		{
		  for(int i=0;i<moltrawr2;i++)
		    {
		      if(isecrawr2[i]==isec[k]&&icsirawr2[i]==icsi[k]&&istrip[k]==9)
			{
		
			  Csifast[k]=csifastr2[i];
			  Csislowpsa[k]=csislowpsar2[i];
			  CsiE[k]=ecsipsar2[i];
			  Csiluce[k]=E2Luce_csi((double)z[k],(double)a[k],CsiE[k]);
			  break;
			}
		    }		  
		}




	    }

	}





  tree->Fill();
}

  tree->Write();
  fout->Write();
  fout->Close();
  return 0;
}
double E2Luce_csi(double z,double a,double E)//FORMULA MICHELA
{
	double ei,coeff,esp;
	double p[8]={12.2424988,288.645213,0.420237154,-0.805045025,1.58972328,0.51370864,0.165660342,0.0181196633};
	
	ei=10.5*a;
	coeff=(p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]));
	esp=p[4]-p[5]*exp(-p[6]*z);
	return (E<=ei)?(coeff*pow(E,esp)):(coeff*((E/ei-1.)*esp+1.)*pow(ei,esp));
}

double Luce2E_csi(double z,double a,double lo) //formula Michela INVERTITA
{
	double ei,loi,coeff,esp;
	double p[8]={12.2424988,288.645213,0.420237154,-0.805045025,1.58972328,0.51370864,0.165660342,0.0181196633};
	
	ei=10.5*a;
	coeff=(p[0]+p[1]*exp(-p[2]*z))*(1.+p[3]*pow(a*z*z,p[7]));
	esp=p[4]-p[5]*exp(-p[6]*z);
	loi=(coeff*pow(ei,esp));
	return (lo<=loi)?(pow(lo/coeff,1./esp)):(((lo/(coeff*pow(ei,esp))-1.)/esp+1)*ei);
}
