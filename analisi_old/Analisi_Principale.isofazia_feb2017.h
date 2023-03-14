#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>
#define RCO(j) (((evento.coderiv[j]>=1000000)&&(evento.coderiv[j]<10000000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
#define GARF(j) (((evento.coderiv[j]>=1000)&&(evento.coderiv[j]<10000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))


void fillh(string hname,double X,double peso)
{
Classe_formule::fillh(hname,X,peso);
 return;

}
void fillh(string hname,double X,double Y,double peso)
{
  Classe_formule::fillh(hname,X,Y,peso);
  return;
}
TH1F *geth1(string hname)
{
  return Classe_formule::geth1(hname);
}
TH2F *geth2(string hname)
{
  return Classe_formule::geth2(hname);
}
TGraphErrors *getg(string hname)
{
  return Classe_formule::getg(hname);
}

void Classe_evento::AnalisiPrincipale()
{
  

  //<<<<<<< Analisi_Principale.h
  //printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);
  // cout<<"molt="<<evento.moltepl<<endl;
  char contQP[100];

  if(Classe_analisi::Getanalisi()->tipo_analisi==23)
    {
      sprintf(contQP,"contornoQPexp");
  if(faziakali.run>21000)
    {
      return;
    }
    }
  else
    {
      sprintf(contQP,"contornoQPsim");
    }
    
  int cqp=-1;
      if(Classe_geo::Getgeo()->ncuts>0)
	{
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	      if(strcmp(contQP,Form("%s",Classe_analisi::Getanalisi()->gggcuts[k]->GetName()))==0)
		{
		  cqp=k;
		  break;
		}	      
	    }

	}



  fillh("hneventi",1.,1);

  vector <int> ibloc(evento.moltepl);
  vector <int> iqua(evento.moltepl);
  vector <int> itel(evento.moltepl);

 for(int j=0;j<evento.moltepl;j++)
    {
      ibloc[j]=-1;
     iqua[j]=-1;
     itel[j]=-1;
     if(evento.coderiv[j]>=10000000)
       {
	 ibloc[j]=evento.coderiv[j]-10000000;
	 
	 ibloc[j]=ibloc[j]/100;
	 int ival=evento.coderiv[j]-10000000-ibloc[j]*100;
	 iqua[j]=ival/4+1;
	 itel[j]=ival-(iqua[j]-1)*4+1;
	 fillh("hcode2",(float)ibloc[j]*100+iqua[j]*10+itel[j],1);
	 // 	 cout<<"bl="<<ibloc[j]<<" q="<<iqua[j]<<" t="<<itel[j]<<" "<<evento.coderiv[j]-10000000<<" "<<faziakali.idtel[evento.indice_originale[j]]<<"evento="<<Classe_analisi::Getanalisi()->nentry<<endl;
	 //cout<<ibloc[j]<<" "<<iqua[j]<<" "<<itel[j]<<" "<<evento.coderiv[j]-10000000<<endl;
	 
       }
     
    }


 for(int j=0;j<evento.moltepl;j++)
   {
     if(evento.rcocode[j]==11&&evento.rcoqf[j]==1)
       {
	 // fillh(Form("hab%dq%dt%d",ibloc[j],iqua[j],itel[j]),faziakali.esi1[evento.indice_originale[j]],evento.z[j],1);
       }
     if(evento.rcocode[j]==11)
       {
	 // fillh(Form("hzb%dq%dt%d",ibloc[j],iqua[j],itel[j]),faziakali.esi1[evento.indice_originale[j]],evento.z[j],1);
	
       }
   }
 int iqp=-1;
 int nqp=0;
 if(cqp>=0)
   {
     for(int j=0;j<evento.moltepl;j++)
       {
	 if(Classe_analisi::Getanalisi()->gggcuts[cqp]->IsInside(evento.vpartlab[j],evento.z[j])==1)
	   {
	     iqp=j;
	     nqp++;
	   }
       }
   }

  for(int j=0;j<evento.moltepl;j++)
    {

      fillh("hros",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      fillh("hzv",evento.vpartlab[j],evento.z[j],1);
      fillh("hzvcm",evento.vpartcm[j],evento.z[j],1);
      
      fillh("hcode",(float)evento.coderiv[j]-10000000,1);
if(evento.rcoqf[j]==1)
    {
      int ival=evento.z[j]*100+evento.a[j]-evento.z[j];
      fillh("hzn",(float)ival,1);
    } 
      if(evento.moltepl>1)
	{
      fillh("hzvm",evento.vpartlab[j],evento.z[j],1);    
      if(j!=iqp && nqp==1)
	{
if(evento.rcoqf[j]==1)
    {

 int ivalqp=evento.z[j]*100+evento.a[j]-evento.z[j];
      fillh("hznqp",(float)ivalqp,1);
    }

	}

	}

   }

  fillh("hmolt",evento.moltepl,1); 




  
}
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void Classe_analisi::RoutineFinale()
{
  cout<<"entro nella routine finale"<<endl;
TGraphErrors *gtz[26];
 for(int j=0;j<26;j++)
   {
     gtz[j]=0;
   }
 for(int j=1;j<26;j++)
   {
     gtz[j]=getg(Form("gtz%d",j));
   }

 TGraphErrors *gqpz[26];
 for(int j=0;j<26;j++)
   {
     gqpz[j]=new TGraphErrors();
   }
 nsuz(geth1("hzn"),getg("gtnz"),gtz); 
 nsuz(geth1("hznqp"),getg("gqpnz"),gqpz);
 


  return;
}

void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100])
{
  int totz[100];
  int totzn[100][100];
float nzval[100];
 float 	 enzval[100];
 float	 err[100];
  
      for(int iz=0;iz<100;iz++)
	{
	      totz[iz]=0;
nzval[iz]=0;
 	 enzval[iz]=0;
 	 err[iz]=0;
	  for(int ia=0;ia<100;ia++)
	    {
	      totzn[iz][ia]=0;
	 
	    }
	}
    

  for(int k=0;k<h->GetXaxis()->GetNbins();k++)
       {
	 int icount=h->GetBinContent(k+1);
	
	 if(icount>0)
	   {

	     int iz=(k+1)/100;
	     int in=k+1-iz*100-1;
	   
	            	 totz[iz]+=icount;
			 totzn[iz][in]+=icount; 
	   }
       }


 for(int iz=3;iz<=25;iz++)
    {
      if(totz[iz]>10)
	{

 	 for(int in=0;in<100;in++)
 	   {
 	     if(totzn[iz][in]>0)
 	       {
		 if(gz[iz]!=0)
		   {
 		 gz[iz]->SetPoint(gz[iz]->GetN(),(float)in,(float)totzn[iz][in]/(float)totz[iz]);
 		 gz[iz]->SetPointError(gz[iz]->GetN()-1,0.0001,((float)totzn[iz][in]/(float)totz[iz])*(sqrt((float)totzn[iz][in])/(float)totzn[iz][in]+sqrt((float)totz[iz])/(float)totz[iz]));
	

		   }
 		 nzval[iz]=nzval[iz]+totzn[iz][in]*(float)in;
		
 	       }
 	   }
 	 if(nzval[iz]>0)
 	   {
 	     float nsuz=nzval[iz]/totz[iz];
 	     for(int in=0;in<100;in++)
 	       {
 		 if(totzn[iz][in]>0)
 		   {

 			 enzval[iz]=enzval[iz]+totzn[iz][in]*(float)in*(float)in;


 		   }
 	       }

 	     enzval[iz]=enzval[iz]/totz[iz];
 	     float nzvaldue=nsuz*nsuz;
 	     enzval[iz]=sqrt((enzval[iz]-nzvaldue)/totz[iz]);


 	      err[iz]=enzval[iz]/((float)iz);
	      
	     

 	     nsuz=nsuz/(float)iz;
	     if(gnz!=0)
	       {
 	 gnz->SetPoint(gnz->GetN(),(float)iz,nsuz);
 	 gnz->SetPointError(gnz->GetN()-1,0.001,err[iz]);
	       }
 	   }
		 if(gz[iz]!=0)
		   {
 gz[iz]->GetXaxis()->SetTitle("N");		
     gz[iz]->SetTitle(Form("Z=%d",iz));
		
		   }
	}//totz>10
    }//iz=1,80
 cout<<gnz->GetName()<<endl;
 gnz->GetXaxis()->SetTitle("Z");
 gnz->GetYaxis()->SetTitle("N/Z");

}
#endif
