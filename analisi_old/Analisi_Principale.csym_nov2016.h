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
//tipo_analisi=21 =>odie


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
  //if(Classe_analisi::Getanalisi()->tipo_analisi==0 || Classe_analisi::Getanalisi()->tipo_analisi==10)
  // {  if(mcevent.isresidue!=1)
  //    {
  //	return;
  //    }
  // }
 fillh("hneventi",30.,1);  
//Scelta trigger (solo exp)
 int bitpat[8];
if(Classe_analisi::Getanalisi()->tipo_analisi>20)
  {
 for(int j=0;j<8;j++)
   {
     bitpat[j]=0;
     int ival=pow(2,j);
     bitpat[j]=expevent.trig&ival;
    
     if(bitpat[j]>0)
       {
	 bitpat[j]=1;
	 fillh("hbit",(float)j,1);
       }
 
 // se il trigger e' 132, bitpat e' 0 0 1 0 0 0 0 1 (2^2+2^7)
   }
 if(bitpat[2]==0){return;}//se OR SIRCO strip1-5==0 E OR SIRCO && OR CSI Garf==0 si butta l'evento
  }


 fillh("hneventi",31.,1);

 vector <int> isec(evento.moltepl);
 vector <int> istrip(evento.moltepl);
 vector <int> icsi(evento.moltepl);

 for(int j=0;j<evento.moltepl;j++)
    {
      isec[j]=-1;
      istrip[j]=-1;
      icsi[j]=-1;
      if(evento.coderiv[j]>=1000000)
	{
	 isec[j]=(evento.coderiv[j]-1000000)/100;
	 istrip[j]=((evento.coderiv[j]-1000000)-isec[j]*100)/10;
	 icsi[j]=(evento.coderiv[j]-1000000)-isec[j]*100-istrip[j]*10;
		}
      if(evento.coderiv[j]>=1000 && evento.coderiv[j]<1000000)
	{
	  istrip[j]=-1;
	  isec[j]=(evento.coderiv[j]-1000)/10;
	  icsi[j]=(evento.coderiv[j]-1000)-10*isec[j];
	  fillh("hrip",15.,1);
	}
    }

if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)//Mc geo
  {
    fillh("hmoltprimari",1.,(float)mcevent.moltprimari,1); 
  }

 //*****************************************
 if(evento.moltepl<=1)
   {
 fillh("hneventi",28.,1);
     return;

   }
 //*****************************************

if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)//Mc geo
  {
    fillh("hmoltprimari",2.,(float)mcevent.moltprimari,1); 
  }



 if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)//Mc geo
   {
     int irco=0;
     int igarf=0;
     for(int j=0;j<evento.moltepl;j++)
       {
	 if(RCO(j))
	   {  
	    
	     irco++;
	   }
	 if(GARF(j))
	   {
	    
	     igarf++;
	   }
       }
     if(irco==0 || igarf==0)
       {
	 return;
       }
   }
 fillh("hneventi",27.,1);

if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)//Mc geo
  {
    fillh("hmoltprimari",3.,(float)mcevent.moltprimari,1); 
  }

 int tagbutta[evento.moltepl];
 for(int j=0;j<evento.moltepl;j++)
   {
     tagbutta[j]=0;
   }

 float ztot;
 float vpcm[500][3],vref[3],vvrel[3];
 float vreffus[3];
 float vpar,vperp;	  

 fillh("hneventi",1.,1.);
 fillh("hmtot",(float)evento.moltepl,1);
  char contfus[100];

if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)
  {
    sprintf(contfus,"contorno100"); //twingo gemini++ geo
    if(Classe_analisi::Getanalisi()->tipo_analisi==10)//gemini geo
   {
     sprintf(contfus,"contorno18");
   }
  }
 if(Classe_analisi::Getanalisi()->tipo_analisi>=20)//exp
  {
 sprintf(contfus,"fus");
  }
if(Classe_analisi::Getanalisi()->tipo_analisi<10)
  {
    sprintf(contfus,"contorno500"); //twingo gemini++ 4pi
 if(Classe_analisi::Getanalisi()->tipo_analisi==0)
   {
     sprintf(contfus,"contorno8");//gemini 4pi
   }
  }


      int nc0fus=-1;
      int nc0break=-1;
      char contbreak[100];
      sprintf(contbreak,"contorno200");
      if(Classe_geo::Getgeo()->ncuts>0)
	{
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	      if(strcmp(contfus,Form("%s",Classe_analisi::Getanalisi()->gggcuts[k]->GetName()))==0)
		{
		  nc0fus=k;
		  break;
		}
	    }	  
	  for(int k=0;k<Classe_geo::Getgeo()->ncuts;k++)
	    {
	      if(strcmp(contbreak,Form("%s",Classe_analisi::Getanalisi()->gggcuts[k]->GetName()))==0)
		{
		  nc0break=k;
		  break;
		}
	    }	  
	}


      int nfus=0;
      int nfussolo=0;
      int nfrag=0;
      int nheavy=0;
      int nbig=0;
      int mlcp=0;
      ztot=0;
      int lista_frag[50],lista_heavy[10],lista_big[10],lista_fus[10],lista_qp[10];
      float zmax=0;
      int nqp=0;
      int ival;
      int jmax=-1;
      float ptotz=0;
      float ptotperp=0;
      int nfis=0;
      int Iqp=-1;
      int nalpha=0;
      int lista_alpha[100];
      int nbreak=0;
  for(int j=0;j<evento.moltepl;j++)
    {
      Classe_formule::da_xyz_a1(evento.vpartcm_x[j],evento.vpartcm_y[j],evento.vpartcm_z[j],vpcm[j]);

      fillh("hzthe",evento.thetalab[j],evento.z[j],1);
      //cout<<Classe_geo::Getgeo()->D[TMath::Nint(evento.z[j])][TMath::Nint(evento.a[j])-TMath::Nint(evento.z[j])]<<endl;
    
      if(Classe_analisi::Getanalisi()->tipo_analisi>20)
	{
	  // if(evento.rcocode[j]==4)
	  // {
	  //   tagbutta[j]=1;
	  //  }
	  if(evento.rcocode[j]==3 && evento.z[j]>15)
	    {
	      fillh("hneventi",181.,1);
	      tagbutta[j]=1;
	      //cout<<"roba di Z "<<evento.z[j]<<" che arriva in cesio???? Si butta"<<endl; 
	    }
	  if(evento.rcocode[j]==1 && evento.z[j]==1) //aggiunto  20 mag 2016
	    {
	      	      fillh("hneventi",183.,1);
	      tagbutta[j]=1;
	    }
	  if(evento.z[j]>(Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt))
	    {
	      fillh("hneventi",180.,1);
	      tagbutta[j]=1;
	      cout<<"Z>Zp+Zt!! "<<evento.z[j]<<endl;
	    }
	  if(GARF(j)&& evento.z[j]>16)
	    {
	      tagbutta[j]=1;
	    }

	}
      if(tagbutta[j]==0)
	{

	  ztot=ztot+evento.z[j];
	  ptotz=ptotz+Classe_formule::amu*evento.a[j]*evento.vpartlab_z[j]/Classe_formule::cluce;
	  ptotperp=ptotperp+sqrt(pow((Classe_formule::amu*evento.a[j]*evento.vpartlab_x[j]/Classe_formule::cluce),2)+pow((Classe_formule::amu*evento.a[j]*evento.vpartlab_y[j]/Classe_formule::cluce),2));
	}
    }
 ptotz=ptotz/(Classe_analisi::Getanalisi()->reazione.ap*Classe_formule::amu*Classe_analisi::Getanalisi()->reazione.vplab/Classe_formule::cluce);
 fillh("hztotptotz",ptotz,ztot,1);
 fillh("hptotperpz",ptotz,ptotperp,1); 

fillh("hneventi",50.,1);

 if((ptotz>1.1)||(ztot>(Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt)))
   {

     fillh("hneventi",51.,1);
     return;
   }
 fillh("hneventi",52.,1);
 if(ptotz<0.2&&ztot<10)//10 febbraio 2016
   {
     fillh("hneventi",53.,1);
     return;
   }
 if(ztot<10)//12 febbraio 2016
   {
      fillh("hneventi",56.,1);
     return;
   }
 fillh("hneventi",54.,1);

 if(Classe_analisi::Getanalisi()->tipo_analisi<20)
   {
     if(Classe_analisi::Getanalisi()->tipo_analisi==0 ||Classe_analisi::Getanalisi()->tipo_analisi==10)
       {
	 fillh("hspin",mcevent.spin,1);
       }
     else
       {
	 fillh("hb",mcevent.par_urto,1);
       }
   }


 int mnobutta=0;
 for(int j=0;j<evento.moltepl;j++)
   {
     if(tagbutta[j]==0)
       {
	 mnobutta++;
      fillh("hzthebis",evento.thetalab[j],evento.z[j],1);

	  fillh("hz",evento.z[j],1);
      ival=evento.z[j]*100+(evento.a[j]-evento.z[j]);
      fillh("hzn",(float)ival,1);
      if(Classe_analisi::Getanalisi()->tipo_analisi>=10 &&evento.coderiv[j]>1000000)//exp o geo
	{
	  if(evento.rcocode[j]==3)
	    {
	      fillh("hqf",(float)evento.rcoqf[j],1);
	    }
	  if(evento.rcocode[j]==3&&evento.rcoqf[j]>1000){fillh("hzn2",(float)ival,1);}
	}

      if(Classe_analisi::Getanalisi()->tipo_analisi<10 && evento.thetalab[j]<17 && evento.thetalab[j]>5.3)//4pi
	{
	 
	  fillh("hzn2",(float)ival,1);
	}
      fillh("hros",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      fillh("hros2",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
      
      if(evento.z[j]>=2.5)
	{
	  
	  fillh("hros3",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
	}
      if(evento.z[j]<2.5)
	{
	  fillh("hros4",evento.thetalab[j]*sin(evento.philab[j]/57.296),evento.thetalab[j]*cos(evento.philab[j]/57.296),1);
	}
      if(evento.coderiv[j]>10000000)
	{
  fillh("hcoderco",(float)evento.coderiv[j]-10000000,1);
	}
      else
	{
  fillh("hcodegarf",(float)evento.coderiv[j]-1000,1);
	}
         fillh("hzv",evento.vpartcm[j],evento.z[j],1);
      fillh("hzvlab",evento.vpartlab[j],evento.z[j],1);
      if(istrip[j]==-1)
	{fillh("hzvgarf",evento.vpartlab[j],evento.z[j],1);}
      if(istrip[j]>=0)
	{fillh("hzvrco",evento.vpartlab[j],evento.z[j],1);}
      if(evento.z[j]>zmax)
	{
	  zmax=evento.z[j];
	  jmax=j;
	}
      if(evento.z[j]<2.5)
	{
	  mlcp++;
	}

      //      if(evento.z[j]>5)
      if(evento.z[j]>3)
	{
	  lista_frag[nfrag]=j;

	  nfrag++;

	}
      if(evento.z[j]>7)
	{
	  lista_qp[nqp]=j;

	  nqp++;

	}

      if(evento.z[j]>=10)
	{
	  lista_heavy[nheavy]=j;
	  nheavy++;
	}
      if(evento.z[j]>18)
	{
	  lista_big[nbig]=j;
	  nbig++;
	}
     
	  if(nc0fus>=0)
	    {
      	      	  if(Classe_analisi::Getanalisi()->gggcuts[nc0fus]->IsInside(evento.vpartcm[j],evento.z[j])==1)
		    {
		      // cout<<"entro2 "<<Classe_analisi::Getanalisi()->nentry<<endl;
		      fillh("hpro",evento.vpartcm[j],evento.z[j],1);
		        lista_fus[nfus]=j; 
		        nfus++; 
		    }
	    }
	
    
	
      else
	{
	  if(evento.z[j]>16&&evento.vpartcm[j]<10)
	{
	  lista_fus[nfus]=j;
	  nfus++;
	}
	}
    

	if(TMath::Nint(evento.z[j])==2&&TMath::Nint(evento.a[j])==4)
	  {
	    lista_alpha[nalpha]=j;
	    nalpha++;
	  }

	}//ibutta==0
    }//giro su moltepl

 float massvec[500],vcmvec[500][3];
 int nbuone=0;
 for(int j=0;j<evento.moltepl;j++)
   {
     if(tagbutta[j]==0)
       {
     massvec[nbuone]=evento.a[j];
     for(int k=0;k<3;k++)
       {
	 vcmvec[nbuone][k]=vpcm[j][k];
       }
     nbuone++;
       }
   }
 fillh("hneventi",48.,1);
 if(nbuone<2)
   {
     fillh("hneventi",49.,1);
     return;
   }
 fillh("hneventi",47.,1);
 float tflow;
 // cout<<evento.moltepl<<endl;
 // Classe_formule::theflow(&tflow,evento.moltepl,massvec,vpcm);
 Classe_formule::theflow(&tflow,nbuone,massvec,vcmvec);
 // cout<<tflow<<endl;
  fillh("htflow",tflow,1); 
  fillh("htflowztot",tflow,ztot,1); 


fillh("hcostflowztot",cos(tflow/57.296),ztot,1);

  //if(Classe_analisi::Getanalisi()->tipo_analisi<20 && Classe_analisi::Getanalisi()->tipo_analisi>=10)//Mc geo
if(Classe_analisi::Getanalisi()->tipo_analisi<20 )//Mc
  {
    if(mcevent.moltprimari==1)
      {
	 fillh("htflowztotM1",tflow,ztot,1);
	 /* if(ztot<20) */
	 /*   { */
	 /*     cout<<tflow<<" "<<ztot<<" "<<evento.moltepl<<endl; */
	 /*     for(int j=0;j<evento.moltepl;j++)  */
	 /*       { */
	 /*   	 cout<<j<<" "<<evento.z[j]<<" "<<evento.vpartcm[j]<<" "<<evento.thetacm[j]<<" "<<Classe_analisi::Getanalisi()->nentry<<endl; */
	 /*       } */
	 /*   } */

fillh("hcostflowztotM1",cos(tflow/57.296),ztot,1);
      }
    if(mcevent.moltprimari==2)
      {
	 fillh("htflowztotM2",tflow,ztot,1);
	 
	 fillh("hcostflowztotM2",cos(tflow/57.296),ztot,1);
	 //if(ztot>22 && tflow<40)
	 //  if(ztot>22 && tflow>60)
	 //if(ztot<22 && tflow>60)
	 // {  
	      //	      cout<<tflow<<" "<<ztot<<" "<<evento.moltepl<<endl; 
	 //   for(int j=0;j<evento.moltepl;j++) 
	 //     { 
		  // cout<<j<<" "<<evento.z[j]<<" "<<evento.vpartcm[j]<<" "<<evento.thetacm[j]<<" "<<Classe_analisi::Getanalisi()->nentry<<endl; 
	 //     } 
	 // } 
      }

  }
  int nf0=0;
  int nf1=0;
  int jfus;
  float zfus;
    for(int j=0;j<evento.moltepl;j++)
      {
	if(tagbutta[j]==0)
	  {
        fillh("htflowz",tflow,evento.z[j],1);
	// if(ztot>=25 && tflow>60)
if(ztot>22 && tflow>60)
   {
	fillh("hzv1",evento.vpartcm[j],evento.z[j],1);
	nf0++;
   }
// if(ztot<25 && tflow<40)
 if(ztot<=22 && tflow<40)
   {
     fillh("hzv2",evento.vpartcm[j],evento.z[j],1);
     nf1++;
   }
	  }
      }
    if(nf0>0 &&nfrag==1) //FUSSOLO
      {

	jfus=-1;
	zfus=0;
	for(int j=0;j<evento.moltepl;j++)
	  {
	    if(tagbutta[j]==0)
	      {
	    fillh("hzv3",evento.vpartcm[j],evento.z[j],1);
	    if(evento.z[j]>zfus)
	      {
		zfus=evento.z[j];
		jfus=j;
	      }
	      }
	  }
	if(jfus>=0)
	  {
	fillh("hzv4",evento.vpartcm[jfus],evento.z[jfus],1);//eventualmente aggiungere condizione vcm fus <15
	if(evento.vpartcm[jfus]<15)
	  {
	fillh("hzfus",evento.z[jfus],1);
fillh("hzv6",evento.vpartcm[jfus],evento.z[jfus],1);
	nfussolo=1;
fillh("hneventi",13.,1);
fillh("htheflowfus",tflow,1);
if(Classe_analisi::Getanalisi()->tipo_analisi<20&&Classe_analisi::Getanalisi()->tipo_analisi>=10)
   {
     fillh("hmprimarifus",(float)mcevent.moltprimari,1); 
     fillh("hmprimarifuszfus",(float)mcevent.moltprimari,evento.z[jfus],1);
     if(mcevent.moltprimari==2)
       {
	 fillh("hdicfuszthe",evento.z[jfus],evento.thetacm[jfus],1);
       }
/*       if(mcevent.moltprimari==2)  */
/*  {  */
/* cout<<"evento="<<Classe_analisi::Getanalisi()->nentry<<" jfus="<<jfus<<" "<<zfus<<" "<<ztot<<" "<<tflow<<endl; */
/*  for(int j=0;j<evento.moltepl;j++) */
/*    { */
/*      cout<<j<<" "<<evento.z[j]<<" "<<evento.vpartcm[j]<<endl; */
/*    } */
/* } */
   }
 for(int k=0;k<3;k++)
   {
     vreffus[k]=vpcm[jfus][k];
   }
 float vpfus[3];
for(int j=0;j<evento.moltepl;j++)
	    {
	      if(j!=jfus)
		{
		  if(tagbutta[j]==0)
		{
		   Classe_formule::vrel(vpcm[j],vreffus,vpfus);
		   float vpfusmod=Classe_formule::modulo(vpfus);
		  float epfus=Classe_formule::v2e(vpfusmod,evento.a[j]);
		  if(TMath::Nint(evento.z[j])==1 && TMath::Nint(evento.a[j])==1)
		    {
		      fillh("hefusp",epfus,1);
		    }
		  if(TMath::Nint(evento.z[j])==2 && TMath::Nint(evento.a[j])==4)
		    {
		      fillh("hefusa",epfus,1);
		    }		  
		}
		}
	    }
	  }
	  }
	  
      }//FINE FUSSOLO
    if(nfrag==2)//Due frammenti
      {
 int igrosso,ipiccolo;
 if(evento.z[lista_frag[0]]>evento.z[lista_frag[1]])
   {
 igrosso=lista_frag[0];
 ipiccolo=lista_frag[1];
 

   }
 else
   {
     igrosso=lista_frag[1];
     ipiccolo=lista_frag[0];

   }

fillh("hznfrag2",evento.z[igrosso],evento.z[ipiccolo],1); 

float vgrosso[3],vpiccolo[3];

 Classe_formule::da_xyz_a1(evento.vpartcm_x[igrosso],evento.vpartcm_y[igrosso],evento.vpartcm_z[igrosso],vgrosso);
 Classe_formule::da_xyz_a1(evento.vpartcm_x[ipiccolo],evento.vpartcm_y[ipiccolo],evento.vpartcm_z[ipiccolo],vpiccolo);
 Classe_formule::vrel(vgrosso,vpiccolo,vvrel);

 float vcmf[3];
 for(int k=0;k<3;k++)
   {
     vcmf[k]=(evento.a[igrosso]*vgrosso[k]+evento.a[ipiccolo]*vpiccolo[k])/(evento.a[igrosso]+evento.a[ipiccolo]);
   } 
     
      float vvrelmod=Classe_formule::modulo(vvrel);
      float thetarel=Classe_formule::thetarel(vgrosso,vpiccolo);

      float zt12=evento.z[igrosso]+evento.z[ipiccolo];
      fillh("hvrelztotnfrag2",Classe_formule::modulo(vvrel),zt12,1);
      fillh("hvrelthetarelnfrag2",Classe_formule::modulo(vvrel),thetarel,1);
     fillh("hztotthetarelnfrag2",zt12,thetarel,1);
      if(nc0break>=0)
	 {
	   if(Classe_analisi::Getanalisi()->gggcuts[nc0break]->IsInside(zt12,thetarel))
	     {
	       nbreak++;
	     }
	 }
      if(nbreak==1)//NORMALIZZAZIONE PER SOTTRAZIONE CARBONIO
	{
	  fillh("htestbreak",zt12,thetarel,1);
	}
      if(nf0>0)//NFRAG=2 E SELEZIONE CENTRALI
      {
	for(int j=0;j<evento.moltepl;j++)
	  {
	    fillh("hzv5",evento.vpartcm[j],evento.z[j],1);
	  }

 fillh("hneventi",4.,1.);
      fillh("hvrelztotnfis0",Classe_formule::modulo(vvrel),zt12,1);
      fillh("hvrelthetarelnfis0",Classe_formule::modulo(vvrel),thetarel,1);
     fillh("hztotthetarelnfis0",zt12,thetarel,1);

     if(thetarel>=120&& vvrelmod>20 && vvrelmod<35) //SELEZIONE FISSIONE
       {

	 nfis++;
fillh("hneventi",5.,1.);
fillh("htheflowfis",tflow,1);
      fillh("hvrelztotnfis",Classe_formule::modulo(vvrel),zt12,1);
      fillh("hvrelthetarelnfis",Classe_formule::modulo(vvrel),thetarel,1);
     fillh("hztotthetarelnfis",zt12,thetarel,1);
     fillh("hzbigzsmall",evento.z[igrosso],evento.z[ipiccolo],1);

if(Classe_analisi::Getanalisi()->tipo_analisi<20&&Classe_analisi::Getanalisi()->tipo_analisi>=10)
   {
     fillh("hmprimarifiss",(float)mcevent.moltprimari,1); 
   }

       }//FINE SELEZIONE FISSIONE
 

      }//FINE NFRAG=2 E SELEZIONE CENTRALI
      }//FINE NFRAG=2
    if(nf1>0 && nqp==1)//SELEZIONE PERIFERICHE E ZQP>7
      {
fillh("hzv0",evento.vpartcm[lista_qp[0]],evento.z[lista_qp[0]],1);
	for(int j=0;j<evento.moltepl;j++)
	  {	
	    if(tagbutta[j]==0)
	      {
fillh("hzv7",evento.vpartcm[j],evento.z[j],1);
	      }
	  }
	if(evento.z[lista_qp[0]]<(Classe_analisi::Getanalisi()->reazione.zp+1) && evento.thetacm[lista_qp[0]]<60)//altri tagli per qp
   {
fillh("hzv9",evento.vpartcm[lista_qp[0]],evento.z[lista_qp[0]],1);
float asyst= Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at;
float aqt=Classe_analisi::Getanalisi()->reazione.ap+Classe_analisi::Getanalisi()->reazione.at-evento.a[lista_qp[0]];
float zsyst=Classe_analisi::Getanalisi()->reazione.zp+Classe_analisi::Getanalisi()->reazione.zt;
float vqt=evento.a[lista_qp[0]]*evento.vpartcm[lista_qp[0]]/aqt;

 float mi=evento.a[lista_qp[0]]*aqt/(aqt+evento.a[lista_qp[0]]);
 float vrel=evento.vpartcm[lista_qp[0]]+vqt;
 float eviola=1.022*pow(zsyst,2)/pow(asyst,0.3333);
 float vviola=Classe_formule::cluce*sqrt(2*eviola/(Classe_formule::amu*(asyst)));
 float vmaxattesa=1.2*vviola*aqt/asyst;
 // if(evento.vpartcm[lista_qp[0]]>=vmaxattesa)
 if(evento.vpartcm[lista_qp[0]]>=0.8*vmaxattesa)//SELEZIONE QP
   {
fillh("hzv8",evento.vpartcm[lista_qp[0]],evento.z[lista_qp[0]],1);
 Iqp=lista_qp[0];
fillh("hneventi",7.,1);
fillh("htheflowdic",tflow,1);
 if(Classe_analisi::Getanalisi()->tipo_analisi<20&&Classe_analisi::Getanalisi()->tipo_analisi>=10)
   {
     fillh("hmprimariDIC",(float)mcevent.moltprimari,1); 
   }
fillh("htheflowselect1",tflow,1);
int ztsel1=0;
int ivalqp=evento.z[Iqp]*100+(evento.a[Iqp]-evento.z[Iqp]);
    if(evento.rcocode[Iqp]==3&&evento.rcoqf[Iqp]>1000)
	{
	  fillh("hznqpselect1",(float)ivalqp,1);
	}
		  vref[0]=evento.vpartcm_x[Iqp];
		  vref[1]=evento.vpartcm_y[Iqp];
		  vref[2]=evento.vpartcm_z[Iqp];
	  for(int j=0;j<evento.moltepl;j++)
	    {
	      if(j!=Iqp)
		{
	      if(tagbutta[j]==0)
		{

	  Classe_formule::vparvperp(vpcm[j],vref,&vpar,&vperp);
	  float vpqp[3];
	  Classe_formule::vrel(vpcm[j],vref,vpqp);
	  float vpqpm,theqp,phiqp;
	  Classe_formule::carpol(&vpqpm,&theqp,&phiqp,vpqp);
	  float epqp=Classe_formule::v2e(vpqpm,evento.a[j]);
	  

	  if(vpar>=0)
	    {ztsel1=ztsel1+evento.z[j];}
	  
		  fillh("hzselect1",evento.z[j],vpar,1);
		  if(TMath::Nint(evento.z[j])<3)
		    {
		      fillh("hlcpselect1",evento.z[j]*10+evento.a[j],vpar,1);
		      if(vpar>Classe_formule::modulo(vref))
			{
			  fillh("hlcpselect1forw",evento.z[j]*10+evento.a[j],vpar,1);
			}

		    }
		  ival=evento.z[j]*100+(evento.a[j]-evento.z[j]);
		  if(Classe_analisi::Getanalisi()->tipo_analisi>=10 &&evento.coderiv[j]>1000000)//exp o geo
		    {
		      if(evento.rcocode[j]==3&&evento.rcoqf[j]>1000){fillh("hznselect1",(float)ival,vpar,1);}
		    }
	  if(TMath::Nint(evento.z[j])==2&&TMath::Nint(evento.a[j])==4)
	    {
	      fillh("hocchiaselect1",vpar,vperp,1);
	      fillh("haeqptheqp",epqp,theqp,1);
	      fillh("hocchiafascioselect1",evento.vpartcm[j]*cos(evento.thetacm[j]/57.296),evento.vpartcm[j]*sin(evento.thetacm[j]/57.296),1);
	    } 
	  if(TMath::Nint(evento.z[j])==2&&TMath::Nint(evento.a[j])==3)
	    {
	      fillh("hocchiHe3select1",vpar,vperp,1);
	      fillh("hHe3eqptheqp",epqp,theqp,1);

	    } 
	  if(TMath::Nint(evento.z[j])==2&&TMath::Nint(evento.a[j])==6)
	    {
	      fillh("hocchiHe6select1",vpar,vperp,1);
	      fillh("hHe6eqptheqp",epqp,theqp,1);

	    } 

	  if(TMath::Nint(evento.z[j])==1&&TMath::Nint(evento.a[j])==1)
	    {
	      fillh("hocchipselect1",vpar,vperp,1);
	      fillh("hpeqptheqp",epqp,theqp,1);
	      fillh("hocchipfascioselect1",evento.vpartcm[j]*cos(evento.thetacm[j]/57.296),evento.vpartcm[j]*sin(evento.thetacm[j]/57.296),1);

	    } 
	  if(TMath::Nint(evento.z[j])==1&&TMath::Nint(evento.a[j])==2)
	    {
	      fillh("hocchidselect1",vpar,vperp,1);
	      fillh("hdeqptheqp",epqp,theqp,1);

	    } 
	  if(TMath::Nint(evento.z[j])==1&&TMath::Nint(evento.a[j])==3)
	    {
	      fillh("hocchitselect1",vpar,vperp,1);
	      fillh("hteqptheqp",epqp,theqp,1);

	    } 

	  if(TMath::Nint(evento.z[j])>2)
	    {
	      fillh("hocchiIMFselect1",vpar,vperp,1);
	      float vrelimf[3];
	      Classe_formule::vrel(vpcm[Iqp],vpcm[j],vrelimf);

	      fillh("hvrelimfselect1",Classe_formule::modulo(vrelimf),1);
	      fillh("hzimfvrelimfselect1",Classe_formule::modulo(vrelimf),evento.z[j],1);
      if(TMath::Nint(evento.z[j])==3)
	{
	  fillh("hocchiz3select1",vpar,vperp,1);
	}

      if(TMath::Nint(evento.z[j])==4)
	{
	  fillh("hocchiz3select1",vpar,vperp,1);
	}
      if(TMath::Nint(evento.z[j])==5)
	{
	  fillh("hocchiz5select1",vpar,vperp,1);
	}
      if(TMath::Nint(evento.z[j])==6)
	{
	  fillh("hocchiz6select1",vpar,vperp,1);
	}
      if(TMath::Nint(evento.z[j])==7)
	{
	  fillh("hocchiz7select1",vpar,vperp,1);
	}
      if(TMath::Nint(evento.z[j])==8)
	{
	  fillh("hocchiz8select1",vpar,vperp,1);
	}

	    }
	  
		}//tagbutta==0
		}//j!selected
	    }//moltepl
  fillh("hztsel1",(float)ztsel1,1);
   }//FINE SELEZIONE QP

   }//FINE ALTRI TAGLI PER QP

      }//FINE SELEZIONE PERIFERICHE E ZQP>7
    int na=0;
    int lista_a[100];
 int naf1=0;
    int lista_af1[100];
 int naf0=0;
    int lista_af0[100];

	for(int j=0;j<evento.moltepl;j++)
	  {
	    if(TMath::Nint(evento.z[j]==2) && TMath::Nint(evento.a[j]==4))
	      {
		lista_a[na]=j;
		na++;
		if(nf1==1)
		  {
		    lista_af1[na]=j;
		    naf1++;
		  }
		if(nf0==1)
		  {
		    lista_af0[na]=j;
		    naf0++;
		  }
	      }
	  }
	//	if(naf1>0 || naf0>0)cout<<na<<" "<<naf1<<" "<<naf0<<endl;
	if(na>4)
	  {
	    fillh("htheflowztot_alpha",tflow,ztot,1);
	  }
     

}
#endif
void Classe_analisi::RoutineFinale()
{
}
