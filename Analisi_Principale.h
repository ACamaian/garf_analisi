#ifndef ANALISI_PRINCIPALE
#define ANALISI_PRINCIPALE
#include "Classe_evento.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include <TRandom.h>
#define RCO(j) (((evento.coderiv[j]>=1000000)&&(evento.coderiv[j]<10000000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
#define GARF(j) (((evento.coderiv[j]>=1000)&&(evento.coderiv[j]<10000))||(Classe_analisi::Getanalisi()->tipo_analisi==0))
//evento.coderiv.push_back(1000000+(isec+1)*100+(istrip+1)*10+icsi+1); RCO (1000111-1000897)
//evento.coderiv.push_back((isec+1)*10+icsi+1+1000); garfield (1011-1248)
//tipo_analisi=210 =>odie

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

  fillh("hneventi",30.,1);

  float asse[3]={0.,0.,1.};

  float ztot=0;
  float ptotz=0;
  float pbeam=0;
  for(int j=0;j<evento.moltepl;j++)
    {
      ztot=ztot+evento.z[j];
      ptotz=ptotz+Classe_formule::amu*evento.a[j]*evento.vpartlab_z[j]/Classe_formule::cluce;	  
    }
  ptotz=ptotz/(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
  pbeam=(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
  fillh("hztotptotz",ptotz,ztot,1);

  if(ztot>(Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt))return;

  if(ptotz>1.1)return;

  fillh("hneventi",31.,1);
 
  fillh("hneventi",1.,1.);
  fillh("hmtot",(float)evento.moltepl,1);
 
  if(Classe_analisi::Getanalisi()->tipo_analisi<200)
    {
      fillh("hb",mcevent.par_urto,1);
    }

  
  int ival[100];
  float nsuzval[100];
  float vplab[100][3];  
  float zmax=0;
  int zlcp=0;
  int nbig=0;
  float vpar,vperp;
  int index[500];
  int vec[500];
  float vpcm[500][3];
  for(int j=0;j<evento.moltepl;j++)
    {
      vec[j]=evento.z[j];
      ival[j]=evento.z[j]*100+evento.a[j]-evento.z[j];
      nsuzval[j]=(float)(evento.a[j]-evento.z[j])/(float)evento.z[j];

      vpcm[j][0]=evento.vpartcm_x[j];
      vpcm[j][1]=evento.vpartcm_y[j];
      vpcm[j][2]=evento.vpartcm_z[j];

      vplab[j][0]=evento.vpartlab_x[j];
      vplab[j][1]=evento.vpartlab_y[j];
      vplab[j][2]=evento.vpartlab_z[j];

      //cout<<Classe_geo::Getgeo()->D[TMath::Nint(evento.z[j])][TMath::Nint(evento.a[j])-TMath::Nint(evento.z[j])]<<endl;// Difetto di massa
      Classe_formule::da_xyz_a1(evento.vpartcm_x[j],evento.vpartcm_y[j],evento.vpartcm_z[j],vpcm[j]);

      fillh("hztheta",evento.thetalab[j],evento.z[j],1);
      
      fillh("hros",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      fillh("hros2",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      
      fillh("hzvcm",evento.vpartcm[j],evento.z[j],1);
      
      //rcocode è idtype cioè PSA, DeltaE -E etc
      if(evento.rcocode[j]==11)fillh("hzvcm11",evento.vpartcm[j],evento.z[j],1);
      if(evento.rcocode[j]==12)fillh("hzvcm12",evento.vpartcm[j],evento.z[j],1);
      if(evento.rcocode[j]==23)fillh("hzvcm23",evento.vpartcm[j],evento.z[j],1);
      if(evento.rcocode[j]==33)fillh("hzvcm33",evento.vpartcm[j],evento.z[j],1);

      fillh("hzvlab",evento.vpartlab[j],evento.z[j],1);
      fillh("hzvzcm",vpcm[j][2],evento.z[j],1);

      if(evento.z[j]>zmax)  zmax=evento.z[j];
      
      fillh("hivalmolt",(float)ival[j],(float)evento.moltepl,1);

    }//giro su moltepl

  TMath::Sort(evento.moltepl,vec,index);

 

}






#endif
void nsuz(TH1F *h,TGraphErrors *gnz,TGraphErrors *gz[100]);
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2);
void nsuz0(TH1F *h,TGraphErrors *gnz);
void Classe_analisi::RoutineFinale()
{
 

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


  for(int iz=0;iz<100;iz++)
    {
      if(totz[iz]>10)//ridurre per aumentare la statistica
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
void ordinisup(TH1F *h,TGraphErrors *gnz,TGraphErrors *g2)
{
  int totz[100];
  int totzn[100][100];
  float nzval[100];
 
  for(int iz=0;iz<100;iz++)
    {
      totz[iz]=0;
      nzval[iz]=0;
 	 
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
  double*z;
  double *nmed;
  double *emed;
  z=gnz->GetX();
  nmed=gnz->GetY();
  emed=gnz->GetEY();
  double nnmed[100],eemed[100];
  float err,err2,somman,sumerr,ey,nzval0;
  for(int j=0;j<gnz->GetN();j++)
    {
   
      //      gnz->GetPoint(j,z,nmed);
      //int iz=(int)z;
      int iz=(int)z[j];
      nnmed[j]=nmed[j]*z[j];
      eemed[j]=emed[j]*z[j];
      err=0;
      err2=0;
      somman=0;
      sumerr=0;
      ey=0;
      nzval0=0;
      for(int in=0;in<100;in++)
	{
	  if(totzn[iz][in]>0)
	    {
	      nzval[iz]=nzval[iz]+totzn[iz][in]*pow((float)in-nnmed[j],2);
	      err=err+totzn[iz][in];
	      // err=err+sqrt(totzn[iz][in])*pow((float)in-nnmed[j],2)+2*((float)in-nnmed[j])*eemed[j]*totzn[iz][in];
		 
	      //err2=err2+sqrt(totzn[iz][in]);
		 
	    }
	}
	 
      // err2=err2/totz[iz];
      //err=err/nzval[iz];
      nzval[iz]=nzval[iz]/totz[iz];
      nzval0=nzval[iz];
      //ey=nzval0*(err+err2);
      nzval[iz]=sqrt(nzval0);
      //ey=0.5*ey/sqrt(nzval0);
      ey=nzval[iz]/sqrt(2*(err-1));
      g2->SetPoint(g2->GetN(),z[j],nzval[iz]);
      g2->SetPointError(g2->GetN()-1,0.0001,ey);
    }

}




void nsuz0(TH1F *h,TGraphErrors *gnz)
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


  for(int iz=0;iz<100;iz++)
    {
      if(totz[iz]>10)//ridurre per aumentare la statistica
	{

	  for(int in=0;in<100;in++)
	    {
	      if(totzn[iz][in]>0)
		{
		 		
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
		 		 
	}//totz>10
    }//iz=1,80
 
  cout<<gnz->GetName()<<endl;
  gnz->GetXaxis()->SetTitle("Z");
  gnz->GetYaxis()->SetTitle("N/Z");

}




