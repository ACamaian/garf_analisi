#include "Classe_evento.h"
#include "Classe_analisi.h"
#include "Classe_geo.h"
#include "Classe_formule.h"
#include "Analisi_Principale.h"
#include <TRandom.h>


void Classe_evento::Leggievento()
{

  // printf("Entro nell'evento %lld\n",Classe_analisi::Getanalisi()->nentry);


 //Convenzione tipo_analisi: 0-9 simulazioni 4p; 100-109 le stesse simulazioni geo; 200 exp human readable; 210 exp compatto  
  if(Classe_analisi::Getanalisi()->tipo_analisi<200) //MC geo o 4pi
    {
      if(Classe_analisi::Getanalisi()->tipo_analisi==0 || Classe_analisi::Getanalisi()->tipo_analisi==100)//Gemini 4p (0) o geo (100)
	{
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Nproducts",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodZ",mcevent.z);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodA",mcevent.a);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodphiCM",mcevent.phicm);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodthetaCM",mcevent.thetacm);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodECM",mcevent.epartcm);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("isfission",&mcevent.fissione);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("fJ",&mcevent.spin);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("isresidue",&mcevent.isresidue);
 if(Classe_analisi::Getanalisi()->ch->GetBranch("Ngamma")!=0)
   {
  Classe_analisi::Getanalisi()->ch->SetBranchAddress("Ngamma",&mcevent.moltgamma);
  Classe_analisi::Getanalisi()->ch->SetBranchAddress("gammaECM",mcevent.egamma);
   }
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("CN",&mcevent.CN);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("Estar",&mcevent.estar);

 //Solo per GeminiDIC
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
 if(Classe_analisi::Getanalisi()->ch->GetBranch("tkel")!=0)
   {
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("zprimari",mcevent.zprimari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("aprimari",mcevent.aprimari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("vcm_primari",mcevent.vcm_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("thecm_primari",mcevent.thecm_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("phi_primari",mcevent.phi_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar_primari",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("spin_primari",mcevent.spin_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("tkel",&mcevent.tkel);
   }
 //fine per GeminiDIC
	}
      if(Classe_analisi::Getanalisi()->tipo_analisi==1||Classe_analisi::Getanalisi()->tipo_analisi==101)//Hipse 4p (1) o geo (101)
	{
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mch",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.z);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.a);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Px",mcevent.Px);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Py",mcevent.Py);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pz",mcevent.Pz);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mtot",&mcevent.mch_mneutr);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("exc",&mcevent.estar);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("origine",mcevent.origine);
      //Origine:
      //  = 0 -> fusion of the QP and QT 
      //  = 1 -> QP
      // = 2 -> QT
      //  > 2 -> other
      // After-decay if a fragments is emitted from another
      // fragment of origin "i" it also have the
      // origin "i"
	}
      if(Classe_analisi::Getanalisi()->tipo_analisi==2||Classe_analisi::Getanalisi()->tipo_analisi==102)//Baiocco 4p (2) o geo (102)
	{
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mtot",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.z);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.a);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Px",mcevent.Px);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Py",mcevent.Py);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pz",mcevent.Pz);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("exc",mcevent.exc);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("discreto",mcevent.discreto);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("spin",&mcevent.spin);
	}


      if(Classe_analisi::Getanalisi()->tipo_analisi==11 || Classe_analisi::Getanalisi()->tipo_analisi==111||Classe_analisi::Getanalisi()->tipo_analisi==12 || Classe_analisi::Getanalisi()->tipo_analisi==112)
	{

Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zsecb);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.asecb);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jx",mcevent.jx);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jy",mcevent.jy);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jz",mcevent.jz);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamdprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vyprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzprimamd);


 if(Classe_analisi::Getanalisi()->tipo_analisi==12 || Classe_analisi::Getanalisi()->tipo_analisi==112)
   {
Classe_analisi::Getanalisi()->ch->SetBranchAddress("discreto",mcevent.discreto);
   }

	} //AMD+HFL (11-111) o baiocco parallelo 12-112
      if(Classe_analisi::Getanalisi()->tipo_analisi==14 || Classe_analisi::Getanalisi()->tipo_analisi==114)
	{

Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zsecb);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.asecb);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);


Classe_analisi::Getanalisi()->ch->SetBranchAddress("mch",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamdprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Px",mcevent.Pxprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Py",mcevent.Pyprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pz",mcevent.Pzprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("spinx",mcevent.jx);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("spiny",mcevent.jy);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("spinz",mcevent.jz);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("molttotsec",&mcevent.mch_mneutr);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("discreto",mcevent.discreto);


	} //HIPSE+HFL (14-114) 


      if(Classe_analisi::Getanalisi()->tipo_analisi==3 || Classe_analisi::Getanalisi()->tipo_analisi==103||Classe_analisi::Getanalisi()->tipo_analisi==5 || Classe_analisi::Getanalisi()->tipo_analisi==105||Classe_analisi::Getanalisi()->tipo_analisi==7 || Classe_analisi::Getanalisi()->tipo_analisi==107||Classe_analisi::Getanalisi()->tipo_analisi==9 || Classe_analisi::Getanalisi()->tipo_analisi==109||Classe_analisi::Getanalisi()->tipo_analisi==10 || Classe_analisi::Getanalisi()->tipo_analisi==110)//AMD sdecay 4p (3) o geo (103) oppure gemini++ decay (5=105) oppure geminif90 decay (7-107) o blob gemini++ (9-109) o AMD+SIMON (10-110)
	{
	
Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zamd);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.aamd);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("j",mcevent.spin_primari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamdprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vyprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzprimamd);

 if(Classe_analisi::Getanalisi()->tipo_analisi==5||Classe_analisi::Getanalisi()->tipo_analisi==105)
{
  if(Classe_analisi::Getanalisi()->ch->GetBranch("isresidue")!=0)
    {
Classe_analisi::Getanalisi()->ch->SetBranchAddress("isresidue",mcevent.isresidue_prim);
    }
}	



	}


      if(Classe_analisi::Getanalisi()->tipo_analisi==15 || Classe_analisi::Getanalisi()->tipo_analisi==115)//AMD+GEMINI++ filtrato con KALIVEDASIM (da usare solo come fosse 4pi ma in realtà è geometria)
	{
	
Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zamd);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.aamd);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("j",mcevent.spin_primari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamdprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vyprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzprimamd);


  if(Classe_analisi::Getanalisi()->ch->GetBranch("isresidue")!=0)
    {
Classe_analisi::Getanalisi()->ch->SetBranchAddress("isresidue",mcevent.isresidue_prim);
    }

 if(Classe_analisi::Getanalisi()->ch->GetBranch("array")!=0) Classe_analisi::Getanalisi()->ch->SetBranchAddress("array",mcevent.array);
  if(Classe_analisi::Getanalisi()->ch->GetBranch("ntele")!=0) Classe_analisi::Getanalisi()->ch->SetBranchAddress("ntele",mcevent.ntele);
  if(Classe_analisi::Getanalisi()->ch->GetBranch("idcode")!=0) Classe_analisi::Getanalisi()->ch->SetBranchAddress("idcode",mcevent.idcode);
  if(Classe_analisi::Getanalisi()->ch->GetBranch("ecode")!=0) Classe_analisi::Getanalisi()->ch->SetBranchAddress("ecode",mcevent.ecode);
  if(Classe_analisi::Getanalisi()->ch->GetBranch("Ameasured")!=0) Classe_analisi::Getanalisi()->ch->SetBranchAddress("Ameasured",mcevent.Ameasured);



	}//tipo_analisi=15-115 AMD+GEMINI++ filtrato con KALIVEDASIM (da usare solo come fosse 4pi ma in realtà è geometria)


      if(Classe_analisi::Getanalisi()->tipo_analisi==4 || Classe_analisi::Getanalisi()->tipo_analisi==104)//twingo+ gemini++
	{

Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zamd);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.aamd);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("j",mcevent.spin_primari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamdprim);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vyprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzprimamd);

// Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
// Classe_analisi::Getanalisi()->ch->SetBranchAddress("eventotwingo",&mcevent.nevorig);

//      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Nproducts",&mcevent.moltepl);
//       Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodZ",mcevent.z);
//            Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodA",mcevent.a);
//        Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodphiCM",mcevent.phicm);
//        Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodthetaCM",mcevent.thetacm);
//        Classe_analisi::Getanalisi()->ch->SetBranchAddress("prodECM",mcevent.epartcm);
//  Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);
// Classe_analisi::Getanalisi()->ch->SetBranchAddress("Zfragprimari",mcevent.zprimari);
// Classe_analisi::Getanalisi()->ch->SetBranchAddress("Afragprimari",mcevent.aprimari);
// Classe_analisi::Getanalisi()->ch->SetBranchAddress("Nfragprimari",&mcevent.moltprimari);
// Classe_analisi::Getanalisi()->ch->SetBranchAddress("Estarfragprimari",mcevent.estar_primari);
//  Classe_analisi::Getanalisi()->ch->SetBranchAddress("isresidue",mcevent.isresidue_prim);
//  Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pxcmfragprimari",mcevent.pxprim);
//  Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pycmfragprimari",mcevent.pyprim);
//  Classe_analisi::Getanalisi()->ch->SetBranchAddress("Pzcmfragprimari",mcevent.pzprim);

	}


      if(Classe_analisi::Getanalisi()->tipo_analisi==13 || Classe_analisi::Getanalisi()->tipo_analisi==113)//COMD da solo fino a tempi lunghi
	{
	
Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltepl);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamd);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamd);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzsecamd);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("estar",mcevent.estar_primari);

 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jx",mcevent.jx);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jy",mcevent.jy);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("jz",mcevent.jz);

	}
      if(Classe_analisi::Getanalisi()->tipo_analisi==6 || Classe_analisi::Getanalisi()->tipo_analisi==106)//twingo+geminif90
	{
Classe_analisi::Getanalisi()->ch->SetBranchAddress("nev",&mcevent.nevorig);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltsec",&mcevent.moltepl);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("b",&mcevent.par_urto);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("zsec",mcevent.zamd);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("asec",mcevent.aamd);
	   
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vxsec",mcevent.vxsecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vysec",mcevent.vysecamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vzsec",mcevent.vzsecamd);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("sorgente",mcevent.origine);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("moltprim",&mcevent.moltprimari);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("zprim",mcevent.zamdprim);
           Classe_analisi::Getanalisi()->ch->SetBranchAddress("aprim",mcevent.aamdprim);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vyprimamd);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzprimamd);

Classe_analisi::Getanalisi()->ch->SetBranchAddress("estarprim",mcevent.estar_primari);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress("j",mcevent.spin_primari);




 

	}  

     if(Classe_analisi::Getanalisi()->tipo_analisi==8 || Classe_analisi::Getanalisi()->tipo_analisi==108)//Langevin4D
       {
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("z",mcevent.zamd);
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("a",mcevent.aamd);
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx",mcevent.vxsecamd);
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy",mcevent.vysecamd);
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz",mcevent.vzsecamd);
	 Classe_analisi::Getanalisi()->ch->SetBranchAddress("molt",&mcevent.moltepl);
	  Classe_analisi::Getanalisi()->ch->SetBranchAddress("rand",&mcevent.nevorig);
	  Classe_analisi::Getanalisi()->ch->SetBranchAddress("epart",mcevent.eservizio);
	}

    }

  if(Classe_analisi::Getanalisi()->tipo_analisi==200)//Exp LETTURA DI ODIE FORMA LEGGIBILE
    {
      if(Classe_geo::Getgeo()->ThereIsPhos==1)
	{
      for(int ip=0;ip<6;ip++)
	{
	  for(int ih=0;ih<9;ih++)
	    {

    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_traw",ip+1,ih+1),&expevent.phos_traw[ip][ih]);
Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_z",ip+1,ih+1),&expevent.phos_z[ip][ih]);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_a",ip+1,ih+1),&expevent.phos_a[ip][ih]);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_qf",ip+1,ih+1),&expevent.phos_qf[ip][ih]);
 Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_cod",ip+1,ih+1),&expevent.phos_cod[ip][ih]);
 //Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("hector_phosbox%d_phos%d_phosid_ga",ip+1,ih+1),&expevent.phos_ga[ip][ih]); //NOn ci sono entuple non compatte con ga dentro

	    }
	}
      //<<<<<<< Classe_evento.cxx
	}//ThereIsPhos
      if(Classe_geo::Getgeo()->ThereIsGarf==1)
	{
    for(int isec=0;isec<24;isec++)
	{
	  for(int icsi=0;icsi<8;icsi++)
	    {

    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_z",isec+1,icsi+1),&expevent.garf_z[isec][icsi]);
    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_a",isec+1,icsi+1),&expevent.garf_a[isec][icsi]);
    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_qf",isec+1,icsi+1),&expevent.garf_qf[isec][icsi]);
    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_theta",isec+1,icsi+1),&expevent.garf_theta[isec][icsi]);
    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_phi",isec+1,icsi+1),&expevent.garf_phi[isec][icsi]);
    Classe_analisi::Getanalisi()->ch->SetBranchAddress(Form("ca_old%d_csi%d_part_mul1_epart",isec+1,icsi+1),&expevent.garf_epart[isec][icsi]);
	    }
	}
	}//ThereIsGarf
    }//tipo_analisi==200
  if(Classe_analisi::Getanalisi()->tipo_analisi==210)//Exp Lettura di ODIE Forma Compatta

  {
	      Classe_analisi::Getanalisi()->ch->SetBranchAddress("value_N",&expcomp.value_N);
	      Classe_analisi::Getanalisi()->ch->SetBranchAddress("value_val",expcomp.value_val);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("value_worker_id",expcomp.value_worker_id);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("value_worker_class_code",expcomp.value_worker_class_code);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("run_num",expcomp.run_num);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("event_num",&expcomp.event_num);

  }//tipo_analisi==210

  if(Classe_analisi::Getanalisi()->tipo_analisi==220)//Exp Lettura di Ntupla Bologna
    {
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mult_1",&expbo.molt);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("numtel_1",expbo.codiceriv);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("vx_1",expbo.vx);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("vy_1",expbo.vy);
       Classe_analisi::Getanalisi()->ch->SetBranchAddress("vz_1",expbo.vz);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("z_1",expbo.z);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("a_1",expbo.a);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("cod_1",expbo.qf);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("ei_1",expbo.e);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("th_1",expbo.theta);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("ph_1",expbo.phi);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("trg_1",expbo.trigger);


       
    }//tipo_analisi==220

  if(Classe_analisi::Getanalisi()->tipo_analisi==230)//Exp Lettura albero Fazietto da Kaliveda
    {
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mtot",&faziakali.mtot);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("run",&faziakali.run);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("idtel",&faziakali.idtel);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Blk",&faziakali.blocco);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Qua",&faziakali.qua);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Tel",&faziakali.tel);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Z",&faziakali.z);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("A",&faziakali.a);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Aid",&faziakali.aid);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("IDCode",&faziakali.idcode);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("ECode",&faziakali.ecode);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("IDType",&faziakali.idtype);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Esi1",&faziakali.esi1);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Esi2",&faziakali.esi2);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Ecsi",&faziakali.ecsi);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("Etot",&faziakali.etot);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHsi1",&faziakali.chsi1);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHsi2",&faziakali.chsi2);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHcsi",&faziakali.chcsi);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("theta",&faziakali.theta);
Classe_analisi::Getanalisi()->ch->SetBranchAddress("phi",&faziakali.phi);




    }//tipo_analisi==230

  if(Classe_analisi::Getanalisi()->tipo_analisi==240)//Exp lettura tree di INDRAFAZIA (output kaliveda)
    {
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("mtot",&indrafazia.mtot);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("run",&indrafazia.run);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("idtel",&indrafazia.idtel);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Blk",&indrafazia.blocco);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Qua",&indrafazia.qua);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Tel",&indrafazia.tel);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Ring",&indrafazia.ring);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Module",&indrafazia.module);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Z",&indrafazia.z);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("A",&indrafazia.a);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Aid",&indrafazia.aid);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("IDCode",&indrafazia.idcode);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("ECode",&indrafazia.ecode);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("IDType",&indrafazia.idtype);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Esi1",&indrafazia.esi1);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Esi2",&indrafazia.esi2);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Ecsi",&indrafazia.ecsi);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("Etot",&indrafazia.etot);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHsi1",&indrafazia.chsi1);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHsi2",&indrafazia.chsi2);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHcsi",&indrafazia.chcsi);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("CHcsifast",&indrafazia.chcsifast);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("theta",&indrafazia.theta);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("phi",&indrafazia.phi);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("array",&indrafazia.array);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("DESi1",&indrafazia.desi1);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("DESi2",&indrafazia.desi2);
      Classe_analisi::Getanalisi()->ch->SetBranchAddress("GT_DT",&indrafazia.gt_dt);



    }//tipo_analisi==240

  //*********************
  Classe_analisi::Getanalisi()->ch->GetEntry(Classe_analisi::Getanalisi()->nentry);
  //*********************
 
      }
void Classe_evento::AnalizzaEvento()
{
  //cout<<"in"<<endl;
  // printf("cluce=%f\n",Classe_formule::cluce);
  if(Classe_analisi::Getanalisi()->tipo_analisi<200)
    {
      if(Classe_analisi::Getanalisi()->tipo_analisi==1 || Classe_analisi::Getanalisi()->tipo_analisi==101)//HIPSE+SIMON
	{
	  mcevent.fissione=-1;
	  
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      float ppp[3];
	      ppp[0]=mcevent.Px[j];
	      ppp[1]=mcevent.Py[j];
	      ppp[2]=mcevent.Pz[j];
	      mcevent.epartcm[j]=Classe_formule::p2e(ppp,mcevent.a[j]);
	      float vv[3];
	      Classe_formule::p2v_cart(ppp,mcevent.a[j],vv);
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],vv);
	     
	      
	    }
	  mcevent.moltneutr=mcevent.mch_mneutr-mcevent.moltepl;
	  mcevent.moltgamma=0;
	}//HIPSE
      if(Classe_analisi::Getanalisi()->tipo_analisi==14 || Classe_analisi::Getanalisi()->tipo_analisi==114)//HIPSE+HFL
	{
	  mcevent.moltneutr=mcevent.mch_mneutr-mcevent.moltepl;//attenzione: qui dentro ci finiscono solo i neutroni primari (i secondari sono nel file delle particelle evaporate)
mcevent.moltgamma=0;


	  mcevent.fissione=-1;
	  mcevent.CN=-1;
	  mcevent.spin=-1;

	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      mcevent.z[j]=mcevent.zsecb[j];
	      mcevent.a[j]=mcevent.asecb[j];
	      mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j];
	      mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j];
	      mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j];
	      
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],mcevent.vpartcm_cart[j]);
	      mcevent.epartcm[j]=Classe_formule::v2e(mcevent.vpartcm[j],mcevent.a[j]);

	    if(TMath::Nint(mcevent.z[j])==0 &&TMath::Nint(mcevent.a[j])==1)
		{
		  mcevent.moltneutr++;
		}
	    }
	    }//HIPSE+HFL


      if(Classe_analisi::Getanalisi()->tipo_analisi==2 || Classe_analisi::Getanalisi()->tipo_analisi==102)//BAIOCCO
	{
	  mcevent.fissione=-1;
	  
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      float ppp[3];
	      ppp[0]=mcevent.Px[j];
	      ppp[1]=mcevent.Py[j];
	      ppp[2]=mcevent.Pz[j];
	      mcevent.epartcm[j]=Classe_formule::p2e(ppp,mcevent.a[j]);
	      float vv[3];
	      Classe_formule::p2v_cart(ppp,mcevent.a[j],vv);
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],vv);
	    if(TMath::Nint(mcevent.z[j])==0 &&TMath::Nint(mcevent.a[j])==1)
		{
		  mcevent.moltneutr++;
		}  
	      
	    }

	  mcevent.moltgamma=0;
	}//BAIOCCO


      if(Classe_analisi::Getanalisi()->tipo_analisi==12 || Classe_analisi::Getanalisi()->tipo_analisi==112)//BAIOCCO PAR
	{

	  mcevent.fissione=-1;
	   mcevent.spin=sqrt(pow(mcevent.jx[0],2)+pow(mcevent.jy[0],2)+pow(mcevent.jz[0],2));
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      mcevent.z[j]=mcevent.zsecb[j];
	      mcevent.a[j]=mcevent.asecb[j];
	      float ppp[3];
	      ppp[0]=mcevent.vxsecamd[j];
	      ppp[1]=mcevent.vysecamd[j];
	      ppp[2]=mcevent.vzsecamd[j];
	      
	      mcevent.epartcm[j]=Classe_formule::v2e(Classe_formule::modulo(ppp),mcevent.asecb[j]);

	     

	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],ppp);
	      if(TMath::Nint(mcevent.zsecb[j])==0 &&TMath::Nint(mcevent.asecb[j])==1)
		{
		  mcevent.moltneutr++;
		}  
	      
	    }

	  mcevent.moltgamma=0;
	}//BAIOCCO PAR



      // if(Classe_analisi::Getanalisi()->tipo_analisi==0 || Classe_analisi::Getanalisi()->tipo_analisi==100||Classe_analisi::Getanalisi()->tipo_analisi==4||Classe_analisi::Getanalisi()->tipo_analisi==104)//Gemini 
      // 	{
      // 		for(unsigned j=0;j<mcevent.moltepl;j++)
      // 	    {
      // 	      if(TMath::Nint(mcevent.z[j])==0 &&TMath::Nint(mcevent.a[j])==1)
      // 		{
      // 		  mcevent.moltneutr++;
      // 		} 
      // 	    }
      // 	}//Gemini

      if(Classe_analisi::Getanalisi()->tipo_analisi==13 || Classe_analisi::Getanalisi()->tipo_analisi==113)//COMD solo
	{
	  mcevent.fissione=-1;
	  mcevent.CN=-1;
	  mcevent.spin=-1;

	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      mcevent.spin_primari[j]=sqrt(pow(mcevent.jx[j],2)+pow(mcevent.jy[j],2)+pow(mcevent.jz[j],2));

	      mcevent.z[j]=(float)mcevent.zamd[j];
	      mcevent.a[j]=(float)mcevent.aamd[j];
	      mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j];
	      mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j];
	      mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j];
	     
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],mcevent.vpartcm_cart[j]);
	      mcevent.epartcm[j]=Classe_formule::v2e(mcevent.vpartcm[j],mcevent.a[j]);
	      if(mcevent.zamd[j]==0 &&mcevent.aamd[j]==1){mcevent.moltneutr++;}
	    }
	  
	  mcevent.moltgamma=0;

	}//COMD solo

     if(Classe_analisi::Getanalisi()->tipo_analisi==3 || Classe_analisi::Getanalisi()->tipo_analisi==103||Classe_analisi::Getanalisi()->tipo_analisi==4 || Classe_analisi::Getanalisi()->tipo_analisi==104||Classe_analisi::Getanalisi()->tipo_analisi==5 || Classe_analisi::Getanalisi()->tipo_analisi==105||Classe_analisi::Getanalisi()->tipo_analisi==6 || Classe_analisi::Getanalisi()->tipo_analisi==106||Classe_analisi::Getanalisi()->tipo_analisi==7 || Classe_analisi::Getanalisi()->tipo_analisi==107||Classe_analisi::Getanalisi()->tipo_analisi==9 || Classe_analisi::Getanalisi()->tipo_analisi==109||Classe_analisi::Getanalisi()->tipo_analisi==10 || Classe_analisi::Getanalisi()->tipo_analisi==110 ||Classe_analisi::Getanalisi()->tipo_analisi==15||Classe_analisi::Getanalisi()->tipo_analisi==115)//AMD (3-103, 5-105, 7-107) oppure twingo+geminif90 oppure twingo+gemini++ oppure blob oppure AMD+SIMON oppure AMD+GEMINI++ filtrato da KALIVEDA SIM
	{
	  mcevent.fissione=-1;
	  mcevent.CN=-1;
	  mcevent.spin=-1;

	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      mcevent.z[j]=(float)mcevent.zamd[j];
	      mcevent.a[j]=(float)mcevent.aamd[j];
	      mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j];
	      mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j];
	      mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j];
	      //mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j]+gRandom->Gaus(0,2.);//prova di Ono per aumentare lo spred Delta p=sqrt(10)*80MeV/c
	      //mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j]+gRandom->Gaus(0,2.);
	      //mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j]+gRandom->Gaus(0,2.);
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],mcevent.vpartcm_cart[j]);
	      mcevent.epartcm[j]=Classe_formule::v2e(mcevent.vpartcm[j],mcevent.a[j]);
	      if(mcevent.zamd[j]==0 &&mcevent.aamd[j]==1){mcevent.moltneutr++;}
	    }
	  
	  mcevent.moltgamma=0;
	}//AMD sdecay
     if(Classe_analisi::Getanalisi()->tipo_analisi==11 || Classe_analisi::Getanalisi()->tipo_analisi==111)
       {
	  mcevent.fissione=-1;
	  mcevent.CN=-1;
	  mcevent.spin=-1;

	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      mcevent.z[j]=mcevent.zsecb[j];
	      mcevent.a[j]=mcevent.asecb[j];
	      mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j];
	      mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j];
	      mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j];
	      //mcevent.vpartcm_cart[j][0]=mcevent.vxsecamd[j]+gRandom->Gaus(0,2.);//prova di Ono per aumentare lo spred Delta p=sqrt(10)*80MeV/c
	      //mcevent.vpartcm_cart[j][1]=mcevent.vysecamd[j]+gRandom->Gaus(0,2.);
	      //mcevent.vpartcm_cart[j][2]=mcevent.vzsecamd[j]+gRandom->Gaus(0,2.);
	      Classe_formule::carpol(&mcevent.vpartcm[j],&mcevent.thetacm[j],&mcevent.phicm[j],mcevent.vpartcm_cart[j]);
	      mcevent.epartcm[j]=Classe_formule::v2e(mcevent.vpartcm[j],mcevent.a[j]);
	      if(mcevent.zamd[j]==0 &&mcevent.aamd[j]==1){mcevent.moltneutr++;}
	    }
	  
	  mcevent.moltgamma=0;
       }


      // printf("Moltepl=%d\n",mcevent.moltepl);
     if(Classe_analisi::Getanalisi()->tipo_analisi!=8 && Classe_analisi::Getanalisi()->tipo_analisi!=108)
       {
     for(unsigned j=0;j<mcevent.moltepl;j++)
    {
      if(TMath::Nint(mcevent.epartcm[j])!=-1)
	{

    mcevent.vpartcm[j]=  Classe_formule::e2v(mcevent.epartcm[j],mcevent.a[j]);

    Classe_formule::polcar(mcevent.vpartcm[j],mcevent.thetacm[j],mcevent.phicm[j],mcevent.vpartcm_cart[j]);

    //    Classe_formule::cm2lab_classica(mcevent.vpartcm_cart[j],mcevent.vpartlab_cart[j],Classe_analisi::Getanalisi()->reazione.vcm);  
    Classe_formule::cm2lab_classica(mcevent.vpartcm_cart[j],mcevent.vpartlab_cart[j],Classe_analisi::Getanalisi()->reazione.vcmv);  
    Classe_formule::carpol(&mcevent.vpartlab[j],&mcevent.thetalab[j],&mcevent.philab[j],mcevent.vpartlab_cart[j]);
       mcevent.epartlab[j]=Classe_formule::v2e(mcevent.vpartlab[j],mcevent.a[j]);

	}//epartcm>0
    }//giro sulla molteplicita' (va fatto sia per mcarlo 4p sia per mcarlo geo

      // printf("Moltepl=%d\n",mcevent.moltepl);
  mcevent.moltcharged=mcevent.moltepl;
  //in hipse in coda alla roba carica ci sono i neutroni (sono contati in mcevent.moltcharged)

  //*************GAMMA
  //********for gammas: APPROXIMATED! elab=ecm=egemini without doppler shift and without geometry filter (all gammas are taken and assigned to Baf 1)
  for(unsigned j=0;j<mcevent.moltgamma;j++)
    {
     mcevent.vpartcm[j+mcevent.moltcharged]= Classe_formule::cluce;
      mcevent.vpartlab[j+mcevent.moltcharged]= Classe_formule::cluce;
      mcevent.epartcm[j+mcevent.moltcharged]=mcevent.egamma[j];
      mcevent.thetacm[j+mcevent.moltcharged]=180*gRandom->Rndm();
     mcevent.phicm[j+mcevent.moltcharged]=180-360*gRandom->Rndm();
     mcevent.philab[j+mcevent.moltcharged]= mcevent.phicm[j+mcevent.moltcharged];
     mcevent.vpartcm_cart[j+mcevent.moltcharged][0]=Classe_formule::cluce*TMath::Sin(mcevent.thetacm[j+mcevent.moltcharged]/57.296)*TMath::Cos(mcevent.phicm[j+mcevent.moltcharged]/57.296);
     mcevent.vpartcm_cart[j+mcevent.moltcharged][1]=Classe_formule::cluce*TMath::Sin(mcevent.thetacm[j+mcevent.moltcharged]/57.296)*TMath::Sin(mcevent.phicm[j+mcevent.moltcharged]/57.296);
     mcevent.vpartcm_cart[j+mcevent.moltcharged][2]=Classe_formule::cluce*TMath::Cos(mcevent.thetacm[j+mcevent.moltcharged]/57.296);
     
    	float pzcm=mcevent.egamma[j]*TMath::Cos(mcevent.thetacm[j+mcevent.moltcharged]/57.296)/Classe_formule::cluce;
     //	mcevent.epartlab[j+mcevent.moltcharged]=(Classe_analisi::Getanalisi()->reazione.gammacm*mcevent.egamma[j]/Classe_formule::cluce+Classe_analisi::Getanalisi()->reazione.gammacm*Classe_analisi::Getanalisi()->reazione.betacm*pzcm)*Classe_formule::cluce;
	mcevent.epartlab[j+mcevent.moltcharged]=(Classe_analisi::Getanalisi()->reazione.gammacmv*mcevent.egamma[j]/Classe_formule::cluce+Classe_analisi::Getanalisi()->reazione.gammacmv*Classe_analisi::Getanalisi()->reazione.betacmv*pzcm)*Classe_formule::cluce;
	
      float plab= mcevent.epartlab[j+mcevent.moltcharged]/Classe_formule::cluce;
      
     //      float pzlab=Classe_analisi::Getanalisi()->reazione.gammacm*pzcm+Classe_analisi::Getanalisi()->reazione.gammacm*Classe_analisi::Getanalisi()->reazione.betacm*mcevent.egamma[j]/Classe_formule::cluce;
      float pzlab=Classe_analisi::Getanalisi()->reazione.gammacmv*pzcm+Classe_analisi::Getanalisi()->reazione.gammacmv*Classe_analisi::Getanalisi()->reazione.betacmv*mcevent.egamma[j]/Classe_formule::cluce;
      
      mcevent.thetalab[j+mcevent.moltcharged]=57.296*TMath::ACos(pzlab/plab);
      
   
       mcevent.vpartlab_cart[j+mcevent.moltcharged][0]=Classe_formule::cluce*TMath::Sin(mcevent.thetalab[j+mcevent.moltcharged]/57.296)*TMath::Cos(mcevent.philab[j+mcevent.moltcharged]/57.296);
     mcevent.vpartlab_cart[j+mcevent.moltcharged][1]=Classe_formule::cluce*TMath::Sin(mcevent.thetalab[j+mcevent.moltcharged]/57.296)*TMath::Sin(mcevent.philab[j+mcevent.moltcharged]/57.296);
     mcevent.vpartlab_cart[j+mcevent.moltcharged][2]=Classe_formule::cluce*TMath::Cos(mcevent.thetalab[j+mcevent.moltcharged]/57.296);
     
      mcevent.z[j+mcevent.moltcharged]=0;
      mcevent.a[j+mcevent.moltcharged]=0;
      mcevent.phicm[j+mcevent.moltcharged]=-28.65;
      mcevent.thetacm[j+mcevent.moltcharged]=110.;
      mcevent.moltepl++;
 
    }//giro sulla molteplicita' di gamma (va fatto sia per mcarlo 4p sia per mcarlo geo
  // FINE GAMMA
  //     printf("Spin=%f\n",mcevent.spin);
       }//tutti i casi tranne langevin
 if(Classe_analisi::Getanalisi()->tipo_analisi==8 || Classe_analisi::Getanalisi()->tipo_analisi==108)
       {
	 mcevent.moltgamma=0;
	 mcevent.moltcharged=0;
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      if(mcevent.aamd[j]>0)
		{
	      mcevent.z[mcevent.moltcharged]=(float)mcevent.zamd[j];
	      mcevent.a[mcevent.moltcharged]=(float)mcevent.aamd[j];
	      mcevent.vpartcm_cart[mcevent.moltcharged][0]=mcevent.vxsecamd[j];
	      mcevent.vpartcm_cart[mcevent.moltcharged][1]=mcevent.vysecamd[j];
	      mcevent.vpartcm_cart[mcevent.moltcharged][2]=mcevent.vzsecamd[j];
	      mcevent.epartcm[mcevent.moltcharged]=mcevent.eservizio[j];
	      if(mcevent.aamd[j]>0)
		{
	      Classe_formule::carpol(&mcevent.vpartcm[mcevent.moltcharged],&mcevent.thetacm[mcevent.moltcharged],&mcevent.phicm[mcevent.moltcharged],mcevent.vpartcm_cart[mcevent.moltcharged]);
	      
		}
	      if(mcevent.zamd[j]==0 &&mcevent.aamd[j]==1){mcevent.moltneutr++;}
		}
	      mcevent.moltcharged++;
	    }
     for(unsigned j=0;j<mcevent.moltcharged;j++)
    {
      if(TMath::Nint(mcevent.epartcm[j])!=-1)
	{


	  //    Classe_formule::cm2lab_classica(mcevent.vpartcm_cart[j],mcevent.vpartlab_cart[j],Classe_analisi::Getanalisi()->reazione.vcm);  
    Classe_formule::cm2lab_classica(mcevent.vpartcm_cart[j],mcevent.vpartlab_cart[j],Classe_analisi::Getanalisi()->reazione.vcmv);  
    Classe_formule::carpol(&mcevent.vpartlab[j],&mcevent.thetalab[j],&mcevent.philab[j],mcevent.vpartlab_cart[j]);
       mcevent.epartlab[j]=Classe_formule::v2e(mcevent.vpartlab[j],mcevent.a[j]);

	}//epartcm>0
    }//giro sulla molteplicita' dei carichi (va fatto sia per mcarlo 4p sia per mcarlo geo


	  
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	    {
	      if(mcevent.aamd[j]==0)
		{
		  mcevent.z[mcevent.moltcharged+mcevent.moltgamma]=0;
	      mcevent.a[mcevent.moltcharged+mcevent.moltgamma]=0;
	      mcevent.epartcm[mcevent.moltcharged+mcevent.moltgamma]=mcevent.eservizio[j];
      mcevent.vpartcm_cart[mcevent.moltcharged+mcevent.moltgamma][0]=mcevent.vxsecamd[j];
      mcevent.vpartcm_cart[mcevent.moltcharged+mcevent.moltgamma][1]=mcevent.vysecamd[j];
      mcevent.vpartcm_cart[mcevent.moltcharged+mcevent.moltgamma][2]=mcevent.vzsecamd[j];
      
	mcevent.thetacm[mcevent.moltcharged+mcevent.moltgamma]=57.296*TMath::ACos(mcevent.vzsecamd[j]/Classe_formule::cluce);
      mcevent.phicm[mcevent.moltcharged+mcevent.moltgamma]=57.296*TMath::ATan2(mcevent.vysecamd[j],mcevent.vxsecamd[j]);
	float pzcm=mcevent.eservizio[j]*TMath::Cos(mcevent.thetacm[mcevent.moltcharged+mcevent.moltgamma]/57.296)/Classe_formule::cluce;

	//	mcevent.epartlab[mcevent.moltcharged+mcevent.moltgamma]=(Classe_analisi::Getanalisi()->reazione.gammacm*mcevent.eservizio[j]/Classe_formule::cluce+Classe_analisi::Getanalisi()->reazione.gammacm*Classe_analisi::Getanalisi()->reazione.betacm*pzcm)*Classe_formule::cluce;
	mcevent.epartlab[mcevent.moltcharged+mcevent.moltgamma]=(Classe_analisi::Getanalisi()->reazione.gammacmv*mcevent.eservizio[j]/Classe_formule::cluce+Classe_analisi::Getanalisi()->reazione.gammacmv*Classe_analisi::Getanalisi()->reazione.betacmv*pzcm)*Classe_formule::cluce;
	
      float plab= mcevent.epartlab[mcevent.moltcharged+mcevent.moltgamma]/Classe_formule::cluce;
      
      //      float pzlab=Classe_analisi::Getanalisi()->reazione.gammacm*pzcm+Classe_analisi::Getanalisi()->reazione.gammacm*Classe_analisi::Getanalisi()->reazione.betacm*mcevent.eservizio[j]/Classe_formule::cluce;
      float pzlab=Classe_analisi::Getanalisi()->reazione.gammacmv*pzcm+Classe_analisi::Getanalisi()->reazione.gammacmv*Classe_analisi::Getanalisi()->reazione.betacmv*mcevent.eservizio[j]/Classe_formule::cluce;
      
      mcevent.thetalab[mcevent.moltcharged+mcevent.moltgamma]=57.296*TMath::ACos(pzlab/plab);
mcevent.philab[mcevent.moltcharged+mcevent.moltgamma]=57.296*TMath::ATan2(mcevent.vysecamd[j],mcevent.vxsecamd[j]);
 
 mcevent.vpartlab_cart[mcevent.moltcharged+mcevent.moltgamma][0]=Classe_formule::cluce*TMath::Sin(mcevent.thetalab[mcevent.moltcharged+mcevent.moltgamma]/57.296)*TMath::Cos(mcevent.phicm[mcevent.moltcharged+mcevent.moltgamma]/57.296);
      mcevent.vpartlab_cart[mcevent.moltcharged+mcevent.moltgamma][1]=Classe_formule::cluce*TMath::Sin(mcevent.thetalab[mcevent.moltcharged+mcevent.moltgamma]/57.296)*TMath::Sin(mcevent.phicm[mcevent.moltcharged+mcevent.moltgamma]/57.296);
      mcevent.vpartlab_cart[mcevent.moltcharged+mcevent.moltgamma][2]=Classe_formule::cluce*TMath::Cos(mcevent.thetalab[mcevent.moltcharged+mcevent.moltgamma]/57.296);
          mcevent.vpartlab[mcevent.moltcharged+mcevent.moltgamma]= Classe_formule::cluce;

      


		  mcevent.moltgamma++;
		}
	    }
	 
       }//solo langevin

      if(Classe_analisi::Getanalisi()->tipo_analisi<200)
    {
      //   Classe_formule::fillh("hb",mcevent.par_urto,1);
    }
  if(Classe_analisi::Getanalisi()->tipo_analisi<100){NoeffettiExp();}//mcarlo 4p

  if(Classe_analisi::Getanalisi()->tipo_analisi>=100 &&Classe_analisi::Getanalisi()->tipo_analisi<200 ) //mcarlo geo
      {
 //**********************
	  //  int nogarf=0;//fa garfield 
  //      nogarf=1;//non fa garfield
	//	int nobaf=1;//non fa i BaF2
	  //	int	nobaf=0; //fa i BaF2
 //*********************************
		  for(unsigned j=0;j<mcevent.moltepl;j++)
      {
  
  int iphos=0;//geometria phos; 0 fuori phos; 1 in phos buono; 2 in phos rotto
  int igarf=0;//geometria garf; 0 fuori garf; 1 in garf
  int ihector=0;//geometria hector; 0 fuori hector; 1 in hector
  //<<<<<<< Classe_evento.cxx
  int irco=0;//geometria RCO; 0 fuori RCO; 1 in RCO
  int ifazietto=0;//geometria Blocchi Fazia; 0 fuori blocchi; 1 nei blocchi


  int sottosoglia_garf=1;//1 sotto soglia per i cesi di garfield e quindi non identificato; 0 sopra soglia per i cesi di garfield e quindi identificato
  int soglia_rco=0;//1 stoppato in camera; 2 stoppato in silicio; 3 stoppato in cesio
  int soglia_phos=0;//1 stoppato nel primo plastico; 2 stoppato nel secondo plastico ; 3 stoppato in cesio; 0 non rivelato (morto nel target o nel mylar di ingresso)

  int sotto_hector=0;//0 non identificato in hector; 1 identificato in hector
  int rivelato_fazietto=-1;

  if((int)j<mcevent.moltcharged) //in Phos e Garfield cerca solo i carichi
    {
      if(Classe_geo::Getgeo()->ThereIsPhos==1)
	{
      iphos=Classe_geo::Getgeo()->InsidePhos(mcevent.thetalab[j],mcevent.philab[j],mcevent.vpartlab[j],&mcevent.tvolo[j],&mcevent.index_phos[j]);
      //           if(iphos!=0)
      if(iphos==1) //esclude quelli rotti
	{
	  soglia_phos=Classe_geo::Getgeo()->Soglie_Phos(mcevent.z[j],mcevent.a[j],mcevent.thetalab[j],mcevent.epartlab[j],mcevent.tvolo[j],mcevent.index_phos[j],&mcevent.zexp[j],&mcevent.aexp[j],&mcevent.tvoloexp[j]);
	  mcevent.soglia[j]=soglia_phos;
	  if(soglia_phos>0)
	    {
	  //********************************
	  mcevent.partbuona[j]=1;
	  //***********************
	  mcevent.codphos[j]=soglia_phos;

	  float dout;	  
	  Classe_geo::Getgeo()->Spalma_Phos(mcevent.index_phos[j],&mcevent.thetalabexp[j],&mcevent.philabexp[j],&dout);
	  if(mcevent.tvoloexp[j]>0) // se non c'e' il tempo di volo le variabili sperimentali legate alla velocita' rimangono a -1
	    {
	  mcevent.vpartlabexp[j]=dout/mcevent.tvoloexp[j];
	  Classe_formule::polcar(mcevent.vpartlabexp[j],mcevent.thetalabexp[j],mcevent.philabexp[j],mcevent.vpartlab_cartexp[j]);
	  if(mcevent.aexp[j]<50)
	    {
	  mcevent.epartlabexp[j]=Classe_formule::v2e(mcevent.vpartlabexp[j],mcevent.aexp[j]);
	    }
	  else
	    {
	      mcevent.epartlabexp[j]=-1;
	    }
	    Classe_formule::lab2cm_classica(mcevent.vpartcm_cartexp[j],mcevent.vpartlab_cartexp[j],Classe_analisi::Getanalisi()->reazione.vcm);
	  Classe_formule::carpol(&mcevent.vpartcmexp[j],&mcevent.thetacmexp[j],&mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);
	  if(mcevent.aexp[j]<50)
	    {
	  mcevent.epartcmexp[j]= Classe_formule::v2e(mcevent.vpartcmexp[j],mcevent.aexp[j]);
	    }
	  else
	    {
	      mcevent.epartcmexp[j]=-1;
	    }
	    }//if tvolo exp>0
	  
	  //<<<<<<< Classe_evento.cxx
	    }//soglia_phos>0

	}//if iphos==1
	}//ThereIsPhos==1
    



      //	 if(iphos==0 && nogarf!=1)
      if(iphos==0 && Classe_geo::Getgeo()->ThereIsGarf==1)
	   {
	     int code_micro=-1;
	     //geometria Garfield; 0 fuori Garfield; 1 in Garfield
	     igarf=Classe_geo::Getgeo()->InsideGarf(mcevent.thetalab[j],mcevent.philab[j],mcevent.z[j],&mcevent.index_garf[j],&code_micro);

	       if(igarf==1)
		 {

		   float estrip1,estrip2,ecsi,luce;
		   int iginocchio;
		   sottosoglia_garf=Classe_geo::Getgeo()->Soglie_Garfield(mcevent.z[j],mcevent.a[j],mcevent.epartlab[j],mcevent.thetalab[j],mcevent.index_garf[j],&estrip1,&estrip2,&ecsi,&luce,&mcevent.zexp[j],&mcevent.aexp[j],&mcevent.epartlabexp[j],&iginocchio,code_micro);
		   mcevent.soglia[j]=sottosoglia_garf;
			   if(sottosoglia_garf==0)
			     {
			       if(code_micro==0)//la micro non si usa per Z=1,2
				 {
				   mcevent.estrip[j][0]=-1;
				   mcevent.estrip[j][1]=-1;
				   mcevent.estrip[j][2]=-1;
				   mcevent.estrip[j][3]=-1;

				   
				 }	   
			   if(code_micro==1)
			     {
			       mcevent.estrip[j][0]=estrip2; //up left
			       mcevent.estrip[j][2]=estrip1; //down left
			     }
			   if(code_micro==2)
			     {
			       mcevent.estrip[j][1]=estrip2; //up right
			       mcevent.estrip[j][3]=estrip1; //down right
			     }
// 			   if(code_micro==101)
// 			     {
// 			       mcevent.estrip[j][0]=-1; //up left
// 			       mcevent.estrip[j][2]=estrip1; //down left
// 			     }
// 			   if(code_micro==102)
// 			     {
// 			       mcevent.estrip[j][1]=-1; //up right
// 			       mcevent.estrip[j][3]=estrip1; //down right
// 			     }
// 			   if(code_micro==201)
// 			     {
// 			       mcevent.estrip[j][0]=estrip2; //up left
// 			       mcevent.estrip[j][2]=-1; //down left
// 			     }
// 			   if(code_micro==202)
// 			     {
// 			       mcevent.estrip[j][1]=estrip2; //up right
// 			       mcevent.estrip[j][3]=-1; //down right
// 			     }

// 			   if(code_micro==301)
// 			     {
// 			       mcevent.estrip[j][0]=-1; //up left
// 			       mcevent.estrip[j][2]=-1; //down left
// 			     }
// 			   if(code_micro==302)
// 			     {
// 			       mcevent.estrip[j][1]=-1; //up right
// 			       mcevent.estrip[j][3]=-1; //down right
// 			     }
// 			   if(code_micro==291)//per camera back
// 			     {
// 			       mcevent.estrip[j][0]=-1;
// 			       mcevent.estrip[j][2]=-1;
// 			     }
// 			   if(code_micro==292)//per camera back
// 			     {
// 			       mcevent.estrip[j][1]=-1;
// 			       mcevent.estrip[j][3]=-1;
// 			     }

			   mcevent.ecsi[j]=ecsi;
			   Classe_geo::Getgeo()->Spalma_Garfield(mcevent.index_garf[j],&mcevent.thetalabexp[j],&mcevent.philabexp[j],code_micro);
	

			   int isec=mcevent.index_garf[j]/10;
			   int icsi=mcevent.index_garf[j]-isec*10;
 			   //int index=100+icsi;


			   if(mcevent.epartlabexp[j]>0)
			     {
			       //****************************
			       mcevent.partbuona[j]=1;
			       //**********************************
			   mcevent.vpartlabexp[j]=Classe_formule::e2v(mcevent.epartlabexp[j],mcevent.aexp[j]);
			   Classe_formule::polcar(mcevent.vpartlabexp[j],mcevent.thetalabexp[j],mcevent.philabexp[j],mcevent.vpartlab_cartexp[j]);
Classe_formule::lab2cm_classica(mcevent.vpartcm_cartexp[j],mcevent.vpartlab_cartexp[j],Classe_analisi::Getanalisi()->reazione.vcm);
	  Classe_formule::carpol(&mcevent.vpartcmexp[j],&mcevent.thetacmexp[j],&mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);
 mcevent.epartcmexp[j]= Classe_formule::v2e(mcevent.vpartcmexp[j],mcevent.aexp[j]);
 mcevent.tvoloexp[j]=Classe_geo::Getgeo()->garf.dist[icsi-1]/mcevent.vpartlabexp[j];
 


			     }//epartlabexp>0
			     }//sottosoglia_garf==0
		 }//igarf==1
	       //<<<<<<< Classe_evento.cxx
	   }//iphos==0 && ThereIsGarf
      if(iphos==0 && igarf==0 && Classe_geo::Getgeo()->ThereIsRCO==1)
	{
	  //geometria RCO; 0 fuori RCO; 1 in RCO
	  irco=Classe_geo::Getgeo()->InsideRCO(mcevent.thetalab[j],mcevent.philab[j],mcevent.z[j],&mcevent.index_rco[j]);
	  if(irco==1)//dentro il RCO
	    {

	      float luce;
	      soglia_rco=Classe_geo::Getgeo()->Soglie_RCO(mcevent.z[j],mcevent.a[j],mcevent.epartlab[j],mcevent.thetalab[j],mcevent.philab[j],&mcevent.index_rco[j],&mcevent.rco_gas[j],&mcevent.rco_si[j],&mcevent.rco_csi[j],&luce,&mcevent.aexp[j],&mcevent.zexp[j],&mcevent.epartlabexp[j],&mcevent.qf[j]);
	      mcevent.soglia[j]=soglia_rco;
	     
	      // nei casi 4 e 5 mancano sia il cesio sia la strip perche' si cade fuori dallo spessore attivo; l'energia in camera non e' l'energia totale.
	      //      	      if(mcevent.soglia_rco[j]>=1 && mcevent.soglia_rco[j]<=3)//1 solo camera ->buona solo E (pero' non corretta per perdite in strati morti perche' lo Z non si sa); 2-3 evento identificabile in Z
			if(soglia_rco==2 || soglia_rco==3)//evento identificabile in Z perche' sono definiti cesio e/ strip
		{
  int isec=mcevent.index_rco[j]/100-1;
  int istrip=(mcevent.index_rco[j]-(isec+1)*100)/10-1;
  int icsi=mcevent.index_rco[j]-(isec+1)*100-(istrip+1)*10-1;


   Classe_geo::Getgeo()->SpalmaRCO(isec,istrip,icsi,&mcevent.thetalabexp[j],&mcevent.philabexp[j]);
      //   mcevent.thetalabexp[j]=mcevent.thetalab[j];//senza spalmatura
      // mcevent.philabexp[j]=mcevent.philab[j];//senza spalmatura
    
		  if(mcevent.epartlabexp[j]>0)
		    {

			       //****************************
			       mcevent.partbuona[j]=1;
			       //**********************************
				   if(soglia_rco==2 || soglia_rco==3)//definiti cesio e/o strip-> si puo' fare identificazione in Z
				   {
			   mcevent.vpartlabexp[j]=Classe_formule::e2v(mcevent.epartlabexp[j],mcevent.aexp[j]);
			   Classe_formule::polcar(mcevent.vpartlabexp[j],mcevent.thetalabexp[j],mcevent.philabexp[j],mcevent.vpartlab_cartexp[j]);
Classe_formule::lab2cm_classica(mcevent.vpartcm_cartexp[j],mcevent.vpartlab_cartexp[j],Classe_analisi::Getanalisi()->reazione.vcm);
	  Classe_formule::carpol(&mcevent.vpartcmexp[j],&mcevent.thetacmexp[j],&mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);
 mcevent.epartcmexp[j]= Classe_formule::v2e(mcevent.vpartcmexp[j],mcevent.aexp[j]);
 mcevent.tvoloexp[j]=Classe_geo::Getgeo()->rco.gas_dist[isec]/mcevent.vpartlabexp[j];
				   }//soglia_rco=2,3
		    }//epartlab>0

		}//soglia_rco==1-3


	    }//irco==1
	}//ThereIs Rco ==1
      
      if(iphos==0 && igarf==0 && irco==0 && Classe_geo::Getgeo()->ThereIsBlocchiFazia==1)
	{
	  float distanza_riv;
	 
	  ifazietto=Classe_geo::Getgeo()->InsideFazietto(mcevent.thetalab[j],mcevent.philab[j],mcevent.vpartlab[j],&distanza_riv,&mcevent.index_fazietto[j]);
	
	    if(ifazietto==1)
	      {
		
		int aid=0;
		rivelato_fazietto=Classe_geo::Getgeo()->Soglie_Fazietto(mcevent.z[j],mcevent.a[j],mcevent.thetalab[j],mcevent.epartlab[j],mcevent.index_fazietto[j],&mcevent.zexp[j],&mcevent.aexp[j],&mcevent.epartlabexp[j],&aid);
		mcevent.qf[j]=aid;
mcevent.soglia[j]=rivelato_fazietto;

	      if(rivelato_fazietto>=0)
		{
	      float dout;
	      Classe_geo::Getgeo()->Spalma_Fazietto(mcevent.index_fazietto[j],&mcevent.thetalabexp[j],&mcevent.philabexp[j],&dout);

	      	      if(mcevent.epartlabexp[j]>0 && mcevent.zexp[j]>0)
	      //   if(mcevent.epartlabexp[j]>0 && mcevent.zexp[j]>0&&rivelato_fazietto>=2)
		    {

			       //****************************
			       mcevent.partbuona[j]=1;
			       //**********************************	      
			   mcevent.vpartlabexp[j]=Classe_formule::e2v(mcevent.epartlabexp[j],mcevent.aexp[j]);
	  Classe_formule::polcar(mcevent.vpartlabexp[j],mcevent.thetalabexp[j],mcevent.philabexp[j],mcevent.vpartlab_cartexp[j]);
Classe_formule::lab2cm_classica(mcevent.vpartcm_cartexp[j],mcevent.vpartlab_cartexp[j],Classe_analisi::Getanalisi()->reazione.vcm);
	  Classe_formule::carpol(&mcevent.vpartcmexp[j],&mcevent.thetacmexp[j],&mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);
 mcevent.epartcmexp[j]= Classe_formule::v2e(mcevent.vpartcmexp[j],mcevent.aexp[j]);
mcevent.tvoloexp[j]=distanza_riv/mcevent.vpartlabexp[j];
 
				   }//epartlabexp[j]>0;
		}//rivelato_fazietto>=0
	      }//fazietto==1
	}//ThereIsBlocchiFazia


    }//fine j<moltcharged
  else // per i soli gamma (j>=moltcharged)
    {
      //	 if(nobaf!=1)
      if(Classe_geo::Getgeo()->ThereIsHector==1)
	   {
	     // ihector=Classe_geo::Getgeo()->InsideHector(mcevent.thetalab[j],mcevent.philab[j],&mcevent.tvolo[j],&mcevent.index_baf[j]);
	     ihector=1; // per questioni di statistica si prendono tutti i gamma
	     mcevent.index_baf[j]=1; // si mettono tutti nel baf 1;
	     mcevent.tvolo[j]=1;//1 ns fisso
	     if(ihector==1)
	       {
		 sotto_hector=Classe_geo::Getgeo()->Soglie_Hector(mcevent.z[j],mcevent.a[j],mcevent.epartlab[j],mcevent.tvolo[j],mcevent.index_baf[j],&mcevent.zexp[j],&mcevent.aexp[j],&mcevent.tvoloexp[j],&mcevent.epartlabexp[j]);
		 if(sotto_hector==1)
		   {
		     //***************************
		     mcevent.partbuona[j]=1;
		     //***************************
			 float dout;
			 Classe_geo::Getgeo()->Spalma_Hector(mcevent.index_baf[j],&mcevent.thetalabexp[j],&mcevent.philabexp[j],&dout);
			 mcevent.vpartlabexp[j]=Classe_formule::cluce;
			 Classe_formule::polcar(mcevent.vpartlabexp[j],mcevent.thetalabexp[j],mcevent.philabexp[j],mcevent.vpartlab_cartexp[j]);


  float plab=mcevent.epartlabexp[j]/Classe_formule::cluce;

 float pzlab=plab*TMath::Cos(mcevent.thetalabexp[j]/57.296);
 mcevent.epartcmexp[j]=mcevent.epartlabexp[j]*Classe_analisi::Getanalisi()->reazione.gammacm/Classe_formule::cluce-pzlab*Classe_analisi::Getanalisi()->reazione.betacm*Classe_analisi::Getanalisi()->reazione.gammacm;
 float pcm=mcevent.epartcmexp[j]/Classe_formule::cluce;
 float pzcm=-Classe_analisi::Getanalisi()->reazione.betacm*Classe_analisi::Getanalisi()->reazione.gammacm*mcevent.epartlabexp[j]/Classe_formule::cluce+Classe_analisi::Getanalisi()->reazione.gammacm*pzlab;
 mcevent.thetacmexp[j]=57.296*TMath::ACos(pzcm/pcm);
 mcevent.phicmexp[j]=mcevent.philabexp[j];
mcevent.vpartcmexp[j]=Classe_formule::cluce;
Classe_formule::polcar(mcevent.vpartcmexp[j],mcevent.thetacmexp[j],mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);

   




		   }//sotto_hector==1
 
	       }//ihector==1

	   }//ThereIsHector
	      }//else per j>=moltcharged

      }//giro sulla molteplicita' di particelle

      }//tipo_analisi==100-109

    }//tipo_analisi<200

  //**********************ODIE COMPATTO***********************
  if(Classe_analisi::Getanalisi()->tipo_analisi==210)//ODIE COMPATTO
    {

      sscanf(expcomp.run_num,"%lld",&expcomp.run);
      expevent.run=expcomp.run;
		for(unsigned j=0;j<expcomp.value_N;j++)
	{
	  //   printf("j=%d %d\n",j,expcomp.value_worker_class_code[j]);

	  if(expcomp.value_worker_class_code[j]==11584)
	    {
	      expevent.trig=expcomp.value_val[j];
	      //  printf("%d trig %f\n",expevent.trig,expcomp.value_val[j]);
	    }
	  if(expcomp.value_worker_class_code[j]==18400)
	    {
	      expevent.tplast=expcomp.value_val[j];
	    }

	  if(expcomp.value_worker_class_code[j]==15328)
	    {
	      int ip=expcomp.value_worker_id[j]/1000;
	      int ih=(expcomp.value_worker_id[j]-ip*1000)/100;
	      expevent.phos_z[ip-1][ih-1]=expcomp.value_val[j];
	      j++;
	      expevent.phos_a[ip-1][ih-1]=expcomp.value_val[j];
	      j++;
	      expevent.phos_qf[ip-1][ih-1]=expcomp.value_val[j];
	      j++;
	      expevent.phos_cod[ip-1][ih-1]=expcomp.value_val[j];
	      j++;
	      expevent.phos_traw[ip-1][ih-1]=expcomp.value_val[j];
	      //*********** DA TOGLIERE SE NON c'E' ga nell'Ntupla
			      j++;
		 expevent.phos_ga[ip-1][ih-1]=expcomp.value_val[j];
	      //*********** FINE
		  //	   printf("%d zphos=%d %d %f %f %f %f %f %f\n",j,ip,ih, expevent.phos_z[ip-1][ih-1],expevent.phos_a[ip-1][ih-1],expevent.phos_qf[ip-1][ih-1],expevent.phos_cod[ip-1][ih-1],expevent.phos_traw[ip-1][ih-1],expevent.phos_ga[ip-1][ih-1]);
		  
		  }//classe 15328
	  if(expcomp.value_worker_class_code[j]==24880)
	    {
	      int isec=expcomp.value_worker_id[j]/10;
	      int icsi=expcomp.value_worker_id[j]-isec*10;
	      expevent.garf_z[isec-1][icsi-1]=expcomp.value_val[j];
	 
	      j++;
	      expevent.garf_a[isec-1][icsi-1]=expcomp.value_val[j];
	      j++;
	      expevent.garf_qf[isec-1][icsi-1]=expcomp.value_val[j];
	      j++;
	      expevent.garf_theta[isec-1][icsi-1]=expcomp.value_val[j];
	      j++;
	      expevent.garf_phi[isec-1][icsi-1]=expcomp.value_val[j];
	      j++;
	      expevent.garf_epart[isec-1][icsi-1]=expcomp.value_val[j];

	      if(TMath::Nint(expevent.garf_z[isec-1][icsi-1])==50 && TMath::Nint(expevent.garf_a[isec-1][icsi-1])==50)
		{
		  expevent.garf_z[isec-1][icsi-1]=0;
		  expevent.garf_a[isec-1][icsi-1]=0;
		}

	      if(expevent.garf_z[isec-1][icsi-1]>0 && expevent.garf_a[isec-1][icsi-1]<0)
		{
		  if(TMath::Nint(expevent.garf_z[isec-1][icsi-1])==1)
		    {
		      expevent.garf_a[isec-1][icsi-1]=1; 
		    }
		  else
		    {
		      expevent.garf_a[isec-1][icsi-1]=2* expevent.garf_z[isec-1][icsi-1];
		    }
		}
	      //  printf("%d zgarf %d %d %f %f %f %f %f %f\n",j,isec,icsi,expevent.garf_z[isec-1][icsi-1],expevent.garf_a[isec-1][icsi-1],expevent.garf_qf[isec-1][icsi-1],expevent.garf_theta[isec-1][icsi-1],expevent.garf_phi[isec-1][icsi-1],expevent.garf_epart[isec-1][icsi-1]);

	    }//classe 24880
	
	
	
	
	
	
	
	  if(expcomp.value_worker_class_code[j]==3728)//BaF2
   {
	      int ibaf=expcomp.value_worker_id[j]/10;

	      exphector.ebaf[ibaf-1]=expcomp.value_val[j];
	      j++;
	      exphector.tbaf[ibaf-1]=expcomp.value_val[j];
	 
	      // printf("baf=%d %f %f\n",ibaf,exphector.ebaf[ibaf-1],exphector.tbaf[ibaf-1]);
   }//classe 3728
	  //<<<<<<< Classe_evento.cxx
	  if(expcomp.value_worker_class_code[j]==24416)//RCO strip
	    {
	      int settore=expcomp.value_worker_id[j]/1000-1;
	      Double_t zloc=expcomp.value_val[j];
	      j++;
	      Double_t aloc=expcomp.value_val[j];
	      j++;
	      //Double_t theloc=expcomp.value_val[j];
	      j++;
	      Double_t philoc=expcomp.value_val[j];
	      j++;	  
 	      Double_t qfloc=expcomp.value_val[j];
	      j++;
	      Double_t codeloc=expcomp.value_val[j];
	      j++;
	      Double_t eloc=expcomp.value_val[j];

	      int istrip0,icsi0=0;
	      istrip0=(expcomp.value_worker_id[j]-1000*(settore+1))/100-1;

	      if(TMath::Nint(codeloc)==1 || TMath::Nint(codeloc)==2)//identificazione da camera-Si o da PSA in Si
		{

		  icsi0=6;//cesio non definito
		}
	      if(TMath::Nint(codeloc)==3)//identificazione da Si-Csi
		{

		  icsi0=philoc-(settore+1)*10-1;
		}
	      //	      TH2F *hhc=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hhc");
	     
	      // if(expcomp.run==111203223004){hhc->Fill((float)(settore+1)*100+((float)istrip0+1)*10+(float)icsi0+1,zloc);}
	     
	      if(istrip0==4 &&(icsi0==2 || icsi0==3))// sono eventi in combinazioni di strip-cesio che la geo NON riconosce
		{
                Classe_analisi::Getanalisi()->ghost_part++;
                 if (Classe_analisi::Getanalisi()->ghost_part%100==1) printf(" --> %Ld ghost particles \n", Classe_analisi::Getanalisi()->ghost_part);

		}
	      else
	      {
		  expevent.rco_z[settore][istrip0][icsi0]=zloc;
		  expevent.rco_a[settore][istrip0][icsi0]=aloc;
		  expevent.rco_code[settore][istrip0][icsi0]=codeloc;//1 IC-Si, 2 PSA, 3 Si-CsI (solo csi in 17344)
		  expevent.rco_qf[settore][istrip0][icsi0]=qfloc;
		  //energia da ottobre 2013
		  expevent.rco_epart[settore][istrip0][icsi0]=eloc;
	      }
		  //il theta e il phi si assegnano con la spalmatura
		


	    }//classe 24420
  if(expcomp.value_worker_class_code[j]==17344)
 	    {
 	      int settore=expcomp.value_worker_id[j]/1000-1;

 	      Double_t zloc=expcomp.value_val[j];
 	      j++;
 	      Double_t aloc=expcomp.value_val[j];
 	      j++;
 	      //Double_t theloc=expcomp.value_val[j];
 	      j++;
 	      //Double_t philoc=expcomp.value_val[j];
 	      j++;	  
  	      Double_t qfloc=expcomp.value_val[j];
 	      j++;
 	      Double_t codeloc=expcomp.value_val[j];
 	  j++;
 	  Double_t eloc=expcomp.value_val[j];

 	      int istrip0,icsi0;
 	      icsi0=(expcomp.value_worker_id[j]-1000*(settore+1))/10-1;

 	      if(TMath::Nint(codeloc)==4)//identificazione da Csi Fast -Slow
 		{
 		  istrip0=8;//strip non definita

 		}
	    
	      //	      TH2F *hhc=(TH2F*)gROOT->GetListOfSpecials()->FindObject("hhc");
	      //  if(expcomp.run==111203223004){ hhc->Fill((float)(settore+1)*100+((float)istrip0+1)*10+(float)icsi0+1,zloc);}

	      if(TMath::Nint(zloc)==10 && TMath::Nint(aloc)==10)//eventuali gamma
		{
 		  expevent.rco_z[settore][istrip0][icsi0]=0.;
 		  expevent.rco_a[settore][istrip0][icsi0]=0.;
		}
	      else
		{
		  expevent.rco_z[settore][istrip0][icsi0]=zloc;
 		  expevent.rco_a[settore][istrip0][icsi0]=aloc;
		}
 		  expevent.rco_code[settore][istrip0][icsi0]=codeloc;
 		  expevent.rco_qf[settore][istrip0][icsi0]=qfloc;

 		  expevent.rco_epart[settore][istrip0][icsi0]=eloc; 

// 		  //il theta e il phi si assegnano con la spalmatura
		


 	    }//classe 17344
 	    
	if(expcomp.value_worker_class_code[j]==10052) { //tvolo GARF
			int isec=expcomp.value_worker_id[j]/10;
			int icsi=expcomp.value_worker_id[j]%10;
			expevent.garf_tvolo[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			//RAW DATA
			expevent.garf_raw_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_cal_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_fast[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_lo[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_slowpsa[isec-1][icsi-1]=expcomp.value_val[j];
			tempo=1;
			


		}//classe 10052
 	if(expcomp.value_worker_class_code[j]==2052) { //tvolo GARF
			int isec=expcomp.value_worker_id[j]/10;
			int icsi=expcomp.value_worker_id[j]%10;
			expevent.garf_tvolo[isec-1][icsi-1]=expcomp.value_val[j];
			tempo=1;
	}//classe 2052 (sostituisce 10052)
	if(expcomp.value_worker_class_code[j]==2000) { //raw GARF
			int isec=expcomp.value_worker_id[j]/10;
			int icsi=expcomp.value_worker_id[j]%10;		
			//RAW DATA
			expevent.garf_raw_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_cal_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_fast[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_lo[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_slowpsa[isec-1][icsi-1]=expcomp.value_val[j];
		
			


		}//classe 2000 (sostituisce 10052)
 	if(expcomp.value_worker_class_code[j]==2001) { //raw GARF
			int isec=expcomp.value_worker_id[j]/10;
			int icsi=expcomp.value_worker_id[j]%10;		
			//RAW DATA
			expevent.garf_raw_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_egashg[isec-1][icsi-1]=expcomp.value_val[j];
			j++;			
			expevent.garf_cal_egas[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_fast[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_lo[isec-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.garf_raw_slowpsa[isec-1][icsi-1]=expcomp.value_val[j];
		
			


		}//classe 2001 (sostituisce 10052 e 2000)
 	    	    
 	    if(expcomp.value_worker_class_code[j]==10053) { //tvolo RCO strip
			int isec=expcomp.value_worker_id[j]/1000;
			int isi=(expcomp.value_worker_id[j]%1000)/100;
			
			double tvolo=expcomp.value_val[j];
			j++;
			double tcode=expcomp.value_val[j];
			j++;
			
			int icsi=(int)(expcomp.value_val[j]+0.5);
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			expevent.rco_tvolo[isec-1][isi-1][icsi-1]=tvolo;
			expevent.rco_tcode[isec-1][isi-1][icsi-1]=tcode;//1 tempo da IC; 2 tempo da Si; 3 tempo da CsI
			//RAW DATA
			expevent.rco_raw_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;

			expevent.rco_raw_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_trise[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_slow[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			
			tempo=1;
		}//classe 10053

	    if(expcomp.value_worker_class_code[j]==1052) { //tvolo RCO strip
			int isec=expcomp.value_worker_id[j]/1000;
			int isi=(expcomp.value_worker_id[j]%1000)/100;
			
			double tvolo=expcomp.value_val[j];
			j++;
			double tcode=expcomp.value_val[j];
			j++;
			
			int icsi=(int)(expcomp.value_val[j]+0.5);
			
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			expevent.rco_tvolo[isec-1][isi-1][icsi-1]=tvolo;
			expevent.rco_tcode[isec-1][isi-1][icsi-1]=tcode;//1 tempo da IC; 2 tempo da Si; 3 tempo da CsI
			tempo=1;

	    }//classe 1052 (sostituisce 10053
	    if(expcomp.value_worker_class_code[j]==1000) { //raw RCO strip
			int isec=expcomp.value_worker_id[j]/1000;
			int isi=(expcomp.value_worker_id[j]%1000)/100;
			int icsi=(int)(expcomp.value_val[j]+0.5);
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			//RAW DATA
			expevent.rco_raw_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;

			expevent.rco_raw_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_trise[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_slow[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			
		
		}//classe 1000 (sostituisce 10053)

 if(expcomp.value_worker_class_code[j]==1001) { //raw RCO strip
			int isec=expcomp.value_worker_id[j]/1000;
			int isi=(expcomp.value_worker_id[j]%1000)/100;
			int icsi=(int)(expcomp.value_val[j]+0.5);
			j++;
			if(icsi<=0) icsi=7;
			if((icsi<1)||(icsi>7)) {
				printf("[DEBUG] ODIE ritorna un codice cesio errato: %d!\n",icsi);
				icsi=7;
			}
			//RAW DATA
			expevent.rco_raw_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_egas[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;

			expevent.rco_raw_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_esi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_trise[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_slow[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_ecsi[isec-1][isi-1][icsi-1]=expcomp.value_val[j];
		
		}//classe 1001 (sostituisce 10053 e 1000)

		
		if(expcomp.value_worker_class_code[j]==10054) { //tvolo RCO exit
			int isec=expcomp.value_worker_id[j]/1000;
			int icsi=(expcomp.value_worker_id[j]%100)/10;
			expevent.rco_tvolo[isec-1][8][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_tcode[isec-1][8][icsi-1]=3.;
			//RAW DATA
			expevent.rco_raw_fast[isec-1][8][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_ecsi[isec-1][8][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_slowpsa[isec-1][8][icsi-1]=expcomp.value_val[j];
			tempo=1;
		}//classe 10054
 	    	
 		if(expcomp.value_worker_class_code[j]==1552) { //tvolo RCO exit
			int isec=expcomp.value_worker_id[j]/1000;
			int icsi=(expcomp.value_worker_id[j]%100)/10;
			expevent.rco_tvolo[isec-1][8][icsi-1]=expcomp.value_val[j];
			expevent.rco_tcode[isec-1][8][icsi-1]=3.;
			tempo=1;
		}//classe 1552 (sostituisce 10054)
		if(expcomp.value_worker_class_code[j]==1500) { //RCO exit raw
			int isec=expcomp.value_worker_id[j]/1000;
			int icsi=(expcomp.value_worker_id[j]%100)/10;

			
			//RAW DATA
			expevent.rco_raw_fast[isec-1][8][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_cal_ecsi[isec-1][8][icsi-1]=expcomp.value_val[j];
			j++;
			expevent.rco_raw_slowpsa[isec-1][8][icsi-1]=expcomp.value_val[j];
			
		}//classe 1500 (sostituisce 10054)
 	    	
 	    
	}//giro su expcompvalue
    }//tipo_analisi==210


    //**********************NTUPLA di BOLOGNA***********************
  if(Classe_analisi::Getanalisi()->tipo_analisi==220)//NTUPLA di BOLOGNA
    {
      expevent.trig=0;
      for(int j=0;j<8;j++)
	{
	  expevent.trig=expevent.trig+pow(2,j)*expbo.trigger[j];
	}
	
      
      for(int j=0;j<expbo.molt;j++)
	{
	 
	  if(expbo.codiceriv[j]<3000)//garfield
	    {
	      int isec=expbo.codiceriv[j]/10;
	      int icsi=expbo.codiceriv[j]-isec*10;
	      expevent.garf_z[isec-1][icsi-1]=expbo.z[j];
	      expevent.garf_a[isec-1][icsi-1]=TMath::Nint(expbo.a[j]);
	      expevent.garf_theta[isec-1][icsi-1]=expbo.theta[j];
	      expevent.garf_phi[isec-1][icsi-1]=expbo.phi[j];
	      expevent.garf_epart[isec-1][icsi-1]=expbo.e[j];
	      expevent.garf_qf[isec-1][icsi-1]=expbo.qf[j];
	      if(expevent.garf_a[isec-1][icsi-1]<=0)
		{
		  expevent.garf_a[isec-1][icsi-1]=Classe_formule::QualeA(expevent.garf_z[isec-1][icsi-1]);
		  expevent.garf_qf[isec-1][icsi-1]=expevent.garf_qf[isec-1][icsi-1]+1000000;

		}
	
	    }//garfield
	  if(expbo.codiceriv[j]>=3000)//RCO
	    {
	      int isec=(expbo.codiceriv[j]-3000)/100;
	      int istrip=((expbo.codiceriv[j]-3000)-isec*100)/10;
	      int icsi=(expbo.codiceriv[j]-3000)-isec*100-istrip*10;


	      if(icsi==0)
		{
		  icsi=7;
		}
	      if(istrip==0)
		{
		  istrip=9;
		}
	      	      if(istrip==5 &&(icsi==3 || icsi==4))// sono eventi in combinazioni di strip-cesio che la geo NON riconosce
		{
                Classe_analisi::Getanalisi()->ghost_part++;
                 if (Classe_analisi::Getanalisi()->ghost_part%100==1) printf(" --> %Ld ghost particles \n", Classe_analisi::Getanalisi()->ghost_part);

		}
		      else
			{
	      expevent.rco_z[isec-1][istrip-1][icsi-1]=expbo.z[j];
	      expevent.rco_a[isec-1][istrip-1][icsi-1]=TMath::Nint(expbo.a[j]);
	      expevent.rco_epart[isec-1][istrip-1][icsi-1]=expbo.e[j];
	      expevent.rco_qf[isec-1][istrip-1][icsi-1]=expbo.qf[j];
	      	   //    TH1F *hzdopo=(TH1F*)gROOT->GetListOfSpecials()->FindObject("hzdopo");
	      // hzdopo->Fill((float)expbo.z[j],1);
	      if(expevent.rco_a[isec-1][istrip-1][icsi-1]<=0)
		{
		  expevent.rco_a[isec-1][istrip-1][icsi-1]=Classe_formule::QualeA(expevent.rco_z[isec-1][istrip-1][icsi-1]);
		  expevent.rco_qf[isec-1][istrip-1][icsi-1]=expevent.rco_qf[isec-1][istrip-1][icsi-1]+1000000;
		}

	 
			}
	    }//RCO
	}//giro su molteplicita'

	

      
    }//tipo analisi==220



    //**********************NTUPLA DI KALIVEDA PER ISOFAZIA***********************
  if(Classe_analisi::Getanalisi()->tipo_analisi==230)//NTUPLA di KALIVEDA per ISOFAZIA
    {
      expevent.trig=-1;

    }//tipo analisi==230

  if(Classe_analisi::Getanalisi()->tipo_analisi==240)//NTUPLA di KALIVEDA per INDRA-FAZIA
    {
      expevent.trig=-1;

    }//tipo analisi==240




if(Classe_analisi::Getanalisi()->tipo_analisi>=100 &&Classe_analisi::Getanalisi()->tipo_analisi<200 ) //mcarlo geo Si buttano i multiple hits
  {
    Escludi_Multiple_Hits();

  }

    Copia_in_Evento(Classe_analisi::Getanalisi()->tipo_analisi);




  if(evento.moltepl>0)
    {
      // cout<<"chiamo "<<evento.moltepl<<endl;
      AnalisiPrincipale();
    }


}

void Classe_evento::AzzeraEvento()
{
  evento.moltepl=0;//19-01-21

if(Classe_analisi::Getanalisi()->tipo_analisi<200)
  {
    mcevent.moltepl=0;
    mcevent.moltcharged=0;
    mcevent.moltgamma=0;
mcevent.fissione=-1;
 mcevent.isresidue=-1;
mcevent.spin=-1;
 mcevent.par_urto=-1;
 mcevent.estar=-1;
 mcevent.moltneutr=0;
 mcevent.mch_mneutr=0;
 mcevent.CN=-1;
 mcevent.tkel=-1;
 for(int j=0;j<2;j++)
   {
     mcevent.zprimari[j]=-1;
     mcevent.aprimari[j]=-1;
     mcevent.vcm_primari[j]=-1;
     mcevent.thecm_primari[j]=-1;
     mcevent.phi_primari[j]=-1;
     mcevent.estar_primari[j]=-1;
     mcevent.spin_primari[j]=-1;

   }

    for(int j=0;j<500;j++)
      {
	mcevent.z[j]=-1;
	mcevent.a[j]=-1;
	mcevent.phicm[j]=-1;
	mcevent.thetacm[j]=-1;
	mcevent.epartcm[j]=-1;
	mcevent.Px[j]=-1;
	mcevent.Py[j]=-1;
	mcevent.Pz[j]=-1;
	mcevent.origine[j]=-1;
	mcevent.discreto[j]=-1;
	mcevent.exc[j]=-1;
	mcevent.index_phos[j]=-1;
	mcevent.index_garf[j]=-1;
	mcevent.index_baf[j]=-1;
	//<<<<<<< Classe_evento.cxx
	mcevent.index_rco[j]=-1;
	mcevent.index_fazietto[j]=-1;
	mcevent.soglia[j]=-1;

	mcevent.ecsi[j]=-1;
	mcevent.rco_gas[j]=-1;
	mcevent.rco_csi[j]=-1;
	mcevent.rco_si[j]=-1;
	mcevent.qf[j]=-1;
	mcevent.zexp[j]=-1;
	mcevent.aexp[j]=-1;
	mcevent.epartlabexp[j]=-1;
	mcevent.thetalabexp[j]=-1;
	mcevent.philabexp[j]=-1;
	mcevent.vpartlabexp[j]=-1;
	mcevent.tvolo[j]=-1;
	mcevent.tvoloexp[j]=-1;
	mcevent.codphos[j]=-1;
	mcevent.partbuona[j]=-1;


        mcevent.array[j]=-1; //per output KALIVEDASIM
        mcevent.ntele[j]=-1;
        mcevent.idcode[j]=-1;
        mcevent.ecode[j]=-2;
        mcevent.Ameasured[j]=-1;


	for(int k=0;k<3;k++)
	  {
	    mcevent.vpartlab_cartexp[j][k]=-1;
	    mcevent.vpartcm_cartexp[j][k]=-1;
	    
	  }
	mcevent.vpartcmexp[j]=-1;
	mcevent.epartcmexp[j]=-1;
	mcevent.thetacmexp[j]=-1;
	mcevent.phicmexp[j]=-1;

	for(int k=0;k<4;k++)
	  {
	    mcevent.estrip[j][k]=-1;
	  }

      }
  }

if(Classe_analisi::Getanalisi()->tipo_analisi>=200)
  {

faziakali.mtot=0;
indrafazia.mtot=0;

    expcomp.value_N=0;//dovrebbe bastare azzerare solo questo
	tempo=0; //30/9/2015 aggiunto tvolo GARF & RCO
    for(int ip=0;ip<6;ip++)
      {
	for(int ih=0;ih<9;ih++)
	  {
	    expevent.phos_tof[ip][ih]=-1;
	    expevent.phos_traw[ip][ih]=-1;
	    expevent.phos_z[ip][ih]=-1;
	    expevent.phos_a[ip][ih]=-1;
	    expevent.phos_qf[ip][ih]=-1;
	    expevent.phos_cod[ip][ih]=-1;
	    expevent.phos_ga[ip][ih]=-1;
	    
	  }
      }

    for(int isec=0;isec<24;isec++)
      {
	for(int icsi=0;icsi<8;icsi++)
	  {
	    expevent.garf_z[isec][icsi]=-1;
	    expevent.garf_a[isec][icsi]=-1;
	    expevent.garf_qf[isec][icsi]=-1;
	    expevent.garf_theta[isec][icsi]=-1;
	    expevent.garf_phi[isec][icsi]=-1;
	    expevent.garf_epart[isec][icsi]=-1;
		expevent.garf_tvolo[isec][icsi]=-1; //30/9/2015 aggiunto tvolo
		expevent.garf_raw_egas[isec][icsi]=-1; 
		expevent.garf_raw_egashg[isec][icsi]=-1; 
		
		expevent.garf_raw_fast[isec][icsi]=-1; 
		expevent.garf_raw_lo[isec][icsi]=-1; 
		expevent.garf_raw_slowpsa[isec][icsi]=-1; 
		expevent.garf_cal_egas[isec][icsi]=-1; 
		expevent.garf_cal_ecsi[isec][icsi]=-1;
	  }
      }

    for(int ip=0;ip<8;ip++)
      {
	exphector.ebaf[ip]=-1;
	exphector.tbaf[ip]=-1;
      }
    expevent.trig=-1;
    expevent.tplast=-1;
    expevent.run=-1;
    for(int ip=0;ip<8;ip++)
      {
	for(int istrip=0;istrip<9;istrip++)
	  {
	    for(int icsi=0;icsi<7;icsi++)
	      {
		expevent.rco_z[ip][istrip][icsi]=-1;
		expevent.rco_a[ip][istrip][icsi]=-1;
		expevent.rco_epart[ip][istrip][icsi]=-1;
		expevent.rco_code[ip][istrip][icsi]=-1;
		expevent.rco_qf[ip][istrip][icsi]=-1;
		expevent.rco_tvolo[ip][istrip][icsi]=-1; //30/9/2015 aggiunto tvolo
		expevent.rco_tcode[ip][istrip][icsi]=-1; //30/9/2015 aggiunto tvolo
		expevent.rco_raw_egas[ip][istrip][icsi]=-1;
	        expevent.rco_cal_egas[ip][istrip][icsi]=-1;
		expevent.rco_raw_esi[ip][istrip][icsi]=-1;
		expevent.rco_cal_esi[ip][istrip][icsi]=-1;
		expevent.rco_raw_trise[ip][istrip][icsi]=-1;
		expevent.rco_raw_fast[ip][istrip][icsi]=-1;
		expevent.rco_raw_slowpsa[ip][istrip][icsi]=-1;
		expevent.rco_cal_ecsi[ip][istrip][icsi]=-1;
	
		expevent.rco_raw_slow[ip][istrip][icsi]=-1;
	
		
	      }
	  }
      }
  }



}


void Classe_evento::NoeffettiExp()
{
  //copia nelle variabili post geometria le variabili originali; da usare in 4p
  if(Classe_analisi::Getanalisi()->tipo_analisi<100)
    {
		for(unsigned j=0;j<mcevent.moltepl;j++)
	{
	  mcevent.zexp[j]=mcevent.z[j];
	  mcevent.aexp[j]=mcevent.a[j];
	  mcevent.epartlabexp[j]=mcevent.epartlab[j];
	  mcevent.thetalabexp[j]=mcevent.thetalab[j];
	  mcevent.philabexp[j]=mcevent.philab[j];
	  mcevent.vpartlabexp[j]=mcevent.vpartlab[j];
	  if(TMath::Abs(Classe_analisi::Getanalisi()->reazione.vcmv-Classe_analisi::Getanalisi()->reazione.vcm)>0.5)//se si sta analizzando il MC con una vcm diversa da quella con cui e' stato prodotto
	    {
	      Classe_formule::lab2cm_classica(mcevent.vpartcm_cartexp[j],mcevent.vpartlab_cart[j],Classe_analisi::Getanalisi()->reazione.vcm);
	      
	      Classe_formule::carpol(&mcevent.vpartcmexp[j],&mcevent.thetacmexp[j],&mcevent.phicmexp[j],mcevent.vpartcm_cartexp[j]);
	      mcevent.epartcmexp[j]= Classe_formule::v2e(mcevent.vpartcmexp[j],mcevent.aexp[j]);
	    }
	  else
	    {

	  mcevent.vpartcmexp[j]=mcevent.vpartcm[j];
	  mcevent.epartcmexp[j]=mcevent.epartcm[j];
	  mcevent.thetacmexp[j]=mcevent.thetacm[j];
	  mcevent.phicmexp[j]=mcevent.phicm[j];
	  for(int k=0;k<3;k++)
	    {
	      mcevent.vpartcm_cartexp[j][k]=mcevent.vpartcm_cart[j][k];
	    }
	    }
	  mcevent.tvoloexp[j]=mcevent.tvolo[j];
	  mcevent.qf[j]=1;//per fazietto
	  //********************
	      mcevent.partbuona[j]=1; //tutte buone
	  //******************
	  for(int k=0;k<3;k++)
	    {
	      	  mcevent.vpartlab_cartexp[j][k]=mcevent.vpartlab_cart[j][k];
		

	    }
	}//fine giro particelle

    }
  else
    {
      return;
    }


}

void Classe_evento::Escludi_Multiple_Hits()
{
	for(unsigned j=0;j<mcevent.moltepl-1;j++)
    {
      if(mcevent.zexp[j]>0)
	{
		for(unsigned k=j+1;k<mcevent.moltepl;k++)
	{
	  if(mcevent.zexp[k]>0)
	    {
	  if(mcevent.index_phos[j]>0 && mcevent.index_phos[k]>0)
	    {
	      if(mcevent.index_phos[j]==mcevent.index_phos[k])
		{
		  mcevent.partbuona[j]=-1;
		  mcevent.partbuona[k]=-1;
		}
	    }
	  if(mcevent.index_garf[j]>0 && mcevent.index_garf[k]>0)
	    {
	      if(mcevent.index_garf[j]==mcevent.index_garf[k])
		{
		  mcevent.partbuona[j]=-1;
		  mcevent.partbuona[k]=-1;
		}
	    }

	  if(mcevent.index_baf[j]>0 && mcevent.index_baf[k]>0)
	    {
	      if(mcevent.index_baf[j]==mcevent.index_baf[k])
		{
		  mcevent.partbuona[j]=-1;
		  mcevent.partbuona[k]=-1;
		}
	    }
	  if(mcevent.index_rco[j]>0 && mcevent.index_rco[k]>0)
	    {
	      if(mcevent.index_rco[j]==mcevent.index_rco[k])
		{
		  mcevent.partbuona[j]=-1;
		  mcevent.partbuona[k]=-1;
		}
	    }
	  if(mcevent.index_fazietto[j]>=0 && mcevent.index_fazietto[k]>=0)
	    {
	      if(mcevent.index_fazietto[j]==mcevent.index_fazietto[k])
		{
		  mcevent.partbuona[j]=-1;
		  mcevent.partbuona[k]=-1;
		 
		}
	    }


	    }
	}//k
	}

    }//j

  return;
}

void Classe_evento::Copia_in_Evento(int tipo_analisi)
{
	      vector <float> vjl;
 
  if(tipo_analisi<200)
    {
      evento.trig=-1;
	  evento.bunch=-1;
	  evento.run=-1;
	  for(unsigned j=0;j<mcevent.moltepl;j++)
	{
	  if(mcevent.partbuona[j]>0)
	    {
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
	      evento.z.push_back(mcevent.zexp[j]);
	      evento.a.push_back(mcevent.aexp[j]);
	      evento.epartcm.push_back(mcevent.epartcmexp[j]);
	      evento.thetacm.push_back(mcevent.thetacmexp[j]);
	      evento.phicm.push_back(mcevent.phicmexp[j]);
	      evento.vpartcm.push_back(mcevent.vpartcmexp[j]);
	      evento.vpartcm_x.push_back(mcevent.vpartcm_cartexp[j][0]);
	      evento.vpartcm_y.push_back(mcevent.vpartcm_cartexp[j][1]);
	      evento.vpartcm_z.push_back(mcevent.vpartcm_cartexp[j][2]);
	      
	      vjl.push_back(mcevent.vpartcm_cartexp[j][0]);
	      vjl.push_back(mcevent.vpartcm_cartexp[j][1]);
	      vjl.push_back(mcevent.vpartcm_cartexp[j][2]);
		evento.vpcm.push_back(vjl);
	      	vjl.clear();
	      
	      evento.epartlab.push_back(mcevent.epartlabexp[j]);
	      evento.thetalab.push_back(mcevent.thetalabexp[j]);
	      evento.philab.push_back(mcevent.philabexp[j]);
	      evento.vpartlab.push_back(mcevent.vpartlabexp[j]);
	      evento.vpartlab_x.push_back(mcevent.vpartlab_cartexp[j][0]);
	      evento.vpartlab_y.push_back(mcevent.vpartlab_cartexp[j][1]);
	      evento.vpartlab_z.push_back(mcevent.vpartlab_cartexp[j][2]);
	      vjl.push_back(mcevent.vpartlab_cartexp[j][0]);
	      vjl.push_back(mcevent.vpartlab_cartexp[j][1]);
	      vjl.push_back(mcevent.vpartlab_cartexp[j][2]);
		evento.vplab.push_back(vjl);
	      	vjl.clear();
	      evento.indice_originale.push_back(j);
	      evento.tvolo.push_back(mcevent.tvoloexp[j]);
		  evento.tbunch.push_back(-1);

	      if(mcevent.Ameasured[j]>-1)//per output KALIVEDASIM	      
	       {
		 if(mcevent.array[j]>-1) evento.codphos.push_back(mcevent.array[j]); //1=fazia, 0=indra	       
		 if(mcevent.ntele[j]>-1) evento.coderiv.push_back(mcevent.ntele[j]); //modulo colpito
		 if(mcevent.idcode[j]>-1) evento.rcocode.push_back(mcevent.idcode[j]); //IDCODE
		 if(mcevent.ecode[j]==-1){mcevent.ecode[j]=0;}

		 if(mcevent.ecode[j]>-2) evento.phosqf.push_back(mcevent.ecode[j]);//equality -> non calibrazione della parte indietro
		 if(mcevent.Ameasured[j]>-1) evento.rcoqf.push_back(mcevent.Ameasured[j]);//1 se massa identificata, 0 se ricalcolata da EAL
	       }

	      else if(mcevent.index_phos[j]>0)
		{
	      evento.coderiv.push_back(mcevent.index_phos[j]);
	      evento.rcocode.push_back(-1);
		}
	      else if(mcevent.index_garf[j]>0)
		{
		  evento.coderiv.push_back(mcevent.index_garf[j]+1000);
		  evento.rcocode.push_back(-1);
		}
	      else if(mcevent.index_baf[j]>0)
		{
		  evento.coderiv.push_back(mcevent.index_baf[j]+10000);
		  evento.rcocode.push_back(-1);
		}
	      //<<<<<<< Classe_evento.cxx
	      else if(mcevent.index_rco[j]>0)
		{
		  evento.coderiv.push_back(mcevent.index_rco[j]+1000000);
		  int ivva=-1;
		  if(mcevent.soglia[j]==2){ivva=1;}
		  if(mcevent.soglia[j]==3){ivva=3;}
		  evento.rcocode.push_back(ivva);

		}
	      else if(mcevent.index_fazietto[j]>=0)
		{
		  evento.coderiv.push_back(mcevent.index_fazietto[j]+10000000);
		  //evento.rcocode.push_back(-1);
		  evento.rcocode.push_back(mcevent.soglia[j]);

		}

	      else
		{
		  //		  printf("Attenzione non risulta ne' index_phos ne' index_garf\n");
		  evento.coderiv.push_back(-1); //dovrebbe essere il caso del montecarlo 4p
		  evento.rcocode.push_back(-1);
		}

	      if(mcevent.Ameasured[j]<0)//caso NON KALIVEDASIM
		{
	      evento.codphos.push_back(mcevent.codphos[j]);
	      //cose che non vengono settate per il montecarlo 


	      evento.phosqf.push_back(-1);
//	      evento.rcoqf.push_back(-1);
	      evento.rcoqf.push_back(mcevent.qf[j]);
		}

	      evento.garfqf.push_back(-1);
	      
	      
	     
	      //evento.phosga.push_back(-1);
	      evento.phosga.push_back(mcevent.epartlab[j]);//energia teorica del mc; serve per fare figure ga/tof

	      evento.tplast.push_back(-1);
	      evento.moltepl++;
	    }
	}
    }
  if(tipo_analisi>=200&& tipo_analisi<220)//EXP leggibile EXP compatto; uscite di ODIE
    {
		
      evento.trig=expevent.trig;
      evento.run=expevent.run;
      if(Classe_geo::Getgeo()->ThereIsPhos==1)
	{
      Riporto();

      for(int ip=0;ip<6;ip++)
	{
	  for(int ih=0;ih<9;ih++)
	    {
	      //	      if((expevent.phos_z[ip][ih]>0.1||expevent.phos_a[ip][ih]>0.1)&&(expevent.phos_tof[ip][ih]>0.))
	      if((expevent.phos_z[ip][ih]>0.1||expevent.phos_a[ip][ih]>0.1)&&(expevent.phos_tof[ip][ih]>0.) && (Classe_geo::Getgeo()->phos.codice[ip][ih]>0)) //si buttano i phos cattivi
		{
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
		  evento.z.push_back(expevent.phos_z[ip][ih]);
		  evento.a.push_back(expevent.phos_a[ip][ih]);
		  float dout;	  
		  int idphos=(ip+1)*10+ih+1;
		  float the,phi;
		  Classe_geo::Getgeo()->Spalma_Phos(idphos,&the,&phi,&dout);
		  evento.thetalab.push_back(the);
		  evento.philab.push_back(phi);

		  evento.tvolo.push_back(expevent.phos_tof[ip][ih]);
		  evento.tbunch.push_back(-1);
		  //		  		  printf("tempo volo=%d %d %f %f\n",ip,ih,expevent.phos_tof[ip][ih],expevent.phos_traw[ip][ih]);
		  evento.codphos.push_back(expevent.phos_cod[ip][ih]);
		  evento.phosqf.push_back(expevent.phos_qf[ip][ih]);
		  evento.vpartlab.push_back(dout/expevent.phos_tof[ip][ih]);

		  float vl[3];

		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],the,phi,vl);
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();

		  evento.coderiv.push_back((ip+1)*10+ih+1);

		  if(evento.a[evento.moltepl]<50)
		    {
		      
		      evento.epartlab.push_back(Classe_formule::v2e(evento.vpartlab[evento.moltepl],evento.a[evento.moltepl]));
		    }
		  else
		    {
		      evento.epartlab.push_back(-1);
		    }

		  float vc[3];
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);
     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
		  if(evento.a[evento.moltepl]<50)
		    {
		      evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));
		    }
		  else
		    {
		      evento.epartcm.push_back(-1);
		    }

		  evento.phosga.push_back(expevent.phos_ga[ip][ih]);

		  // cose che non vengono settate per i phos exp:
		  evento.indice_originale.push_back(-1);
		  evento.garfqf.push_back(-1);
		  evento.rcoqf.push_back(-1);
		  evento.rcocode.push_back(-1);

		  evento.tplast.push_back(-1);
		 


		  evento.moltepl++;
		}
	    }
	}//fine phos
	}//ThereIsPhos
	
	evento.bunch=-1;
	if(tempo) RiportoGarfRCo();
	
      if(Classe_geo::Getgeo()->ThereIsGarf==1)
	{
      for(int isec=0;isec<24;isec++)
	{
	  for(int icsi=0;icsi<8;icsi++)
	    {
	      	      
	      if((expevent.garf_z[isec][icsi]>0. && expevent.garf_a[isec][icsi]>0.)&&expevent.garf_epart[isec][icsi]>0&&Classe_geo::Getgeo()->garf.codice[icsi][isec]==0)
		{
		  int code_micro;
		  if(expevent.garf_phi[isec][icsi]<=0||TMath::Nint(expevent.garf_phi[isec][icsi])==3)//3 vuol dire condivisa; 0 o -1 micro non definita
		    {
		      code_micro=0;
		    }
		  else//1 sx 2 dx
		    {
		      code_micro=TMath::Nint(expevent.garf_phi[isec][icsi]);
		    }
		  

		   int non_accetta=0;
		   if(expevent.garf_z[isec][icsi]>2.5 && code_micro>0 && Classe_geo::Getgeo()->garf.micro_rotte[isec][icsi][code_micro-1]==0)
		     {
		       non_accetta=1;
		     }
		   if(code_micro==0 && expevent.garf_z[isec][icsi]>2.5)
		     {
		       non_accetta=1;
		     }
		   
		  if(non_accetta==0)

		    {
		 evento.isformix.push_back(-1); //inserito per mixatore 29-01-21 
		  evento.z.push_back(expevent.garf_z[isec][icsi]);
		  evento.a.push_back(expevent.garf_a[isec][icsi]);
		  evento.epartlab.push_back(expevent.garf_epart[isec][icsi]);
		  
		  if(expevent.garf_raw_lo[isec][icsi]>0)
		    {
		      expevent.garf_cal_ecsi[isec][icsi]=Classe_geo::Getgeo()->Luce2E_csi(expevent.garf_z[isec][icsi],expevent.garf_a[isec][icsi],expevent.garf_raw_lo[isec][icsi]);
		    
		    }


		   float thetout,phiout;
		   int codice=(isec+1)*10+icsi+1;
		 Classe_geo::Getgeo()->Spalma_Garfield(codice,&thetout,&phiout,code_micro);
		 
		  if(expevent.garf_theta[isec][icsi]!=-1)
		    {
		      if((expevent.garf_theta[isec][icsi]<Classe_geo::Getgeo()->garf.themincsi[icsi]||expevent.garf_theta[isec][icsi]>Classe_geo::Getgeo()->garf.themaxcsi[icsi])||(expevent.garf_z[isec][icsi]<2.5))//i p e le alpha sono sempre spalmati (non si usa il tdrift)
			{
		
		      evento.thetalab.push_back(thetout);//prende quello spalmato
		      
			}
		      else
			{
		      evento.thetalab.push_back(expevent.garf_theta[isec][icsi]);//prende quello che arriva da odie
			}      
		    }
		  else
		    {
		      // evento.thetalab.push_back(Classe_geo::Getgeo()->garf.thecsi[icsi]);//angolo centrale
		    
		      evento.thetalab.push_back(thetout);//prende quello spalmato
		    }
		 
		  evento.garfqf.push_back(expevent.garf_qf[isec][icsi]);
		  evento.rcoqf.push_back(-1);
	      evento.rcocode.push_back(-1);

		  evento.coderiv.push_back((isec+1)*10+icsi+1+1000);
		  evento.philab.push_back(phiout);
		  

		  float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
	
		 
		  evento.vpartlab.push_back(vv);
		evento.tvolo.push_back(Classe_geo::Getgeo()->garf.dist[icsi]/vv);
		evento.tbunch.push_back(expevent.garf_tvolo[isec][icsi]-200.*(float)(evento.bunch+1));
		  float vl[3];
		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();
		  float vc[3];
	
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);

     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
 evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));
 
		  //cose che non vengono settate per garfield exp
 evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);
		  evento.tplast.push_back(-1);



evento.moltepl++;
		    }//non_accetta=0
		}//evento buono
	      if(TMath::Nint(expevent.garf_z[isec][icsi])==0 &&TMath::Nint(expevent.garf_a[isec][icsi])==0 && expevent.garf_epart[isec][icsi]>0&&Classe_geo::Getgeo()->garf.codice[icsi][isec]==0)// evento gamma in garfield
		{
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
	      evento.z.push_back(0.);
	      evento.a.push_back(0.);
	      int codice=(isec+1)*10+icsi+1;
	      float thetout,phiout;
	      int code_micro=0;
	      Classe_geo::Getgeo()->Spalma_Garfield(codice,&thetout,&phiout,code_micro);
	      evento.thetalab.push_back(thetout);
	      evento.philab.push_back(phiout);
	      evento.tvolo.push_back(-1);
		  evento.tbunch.push_back(expevent.garf_tvolo[isec][icsi]-200.*(float)(evento.bunch+1));
	      evento.epartcm.push_back(-1);
	      evento.thetacm.push_back(-1);
	      evento.phicm.push_back(-1);
	      evento.vpartcm.push_back(-1);
	      evento.vpartcm_x.push_back(-1);
	      evento.vpartcm_y.push_back(-1);
	      evento.vpartcm_z.push_back(-1);

	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
	      evento.epartlab.push_back(expevent.garf_epart[isec][icsi]);
	      evento.vpartlab.push_back(-1);
	      evento.vpartlab_x.push_back(-1);
	      evento.vpartlab_y.push_back(-1);
	      evento.vpartlab_z.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vplab.push_back(vjl);
		vjl.clear();
	      evento.indice_originale.push_back(-1);
	      evento.codphos.push_back(-1);
	      evento.phosqf.push_back(-1);
	      evento.phosga.push_back(-1);
	      
	      evento.rcoqf.push_back(-1);
	      evento.rcocode.push_back(-1);
	      evento.garfqf.push_back(expevent.garf_qf[isec][icsi]);
	      evento.tplast.push_back(-1);

	      evento.coderiv.push_back((isec+1)*10+icsi+1+1000);





	      evento.moltepl++;


		}//evento gamma in garfield

	    }
	}//fine garf
	}//ThereIsGarf==1
      if(Classe_geo::Getgeo()->ThereIsRCO==1)
	{
	  for(int isec=0;isec<8;isec++)
	    {
	      for(int istrip=0;istrip<9;istrip++)
		{
		  for(int icsi=0;icsi<7;icsi++)
		    {
		      int non_accetta=0;
		      if(expevent.rco_z[isec][istrip][icsi]>0.)//Se esiste uno Z identificato
			{
			  if(expevent.rco_code[isec][istrip][icsi]==1&&(Classe_geo::Getgeo()->rco.codice_gas[isec]==1||Classe_geo::Getgeo()->rco.codice_strip[isec][istrip]==1))
			    {
			      non_accetta=1;
			    }
			  if(expevent.rco_code[isec][istrip][icsi]==2&&Classe_geo::Getgeo()->rco.codice_strip[isec][istrip]==1)
			{
			  non_accetta=1;
			}
			  if(expevent.rco_code[isec][istrip][icsi]==4 &&Classe_geo::Getgeo()->rco.codice_csi[isec][icsi]==1)
			    {
			      non_accetta=1;
			    }

			  if(expevent.rco_code[isec][istrip][icsi]==3&&(Classe_geo::Getgeo()->rco.codice_csi[isec][icsi]==1||Classe_geo::Getgeo()->rco.codice_strip[isec][istrip]==1))
			    {
			      non_accetta=1;
			    }
		

			  if(non_accetta==0)
			    {
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
			  evento.z.push_back(expevent.rco_z[isec][istrip][icsi]);
			  if(expevent.rco_a[isec][istrip][icsi]<=0)
			    {

			      if(TMath::Nint(expevent.rco_z[isec][istrip][icsi])==1)
				{
				  evento.a.push_back(1.); 
				}
			      else
				{
			      evento.a.push_back(2*expevent.rco_z[isec][istrip][icsi]);
				}
				}
			  else
			    {
			  evento.a.push_back(expevent.rco_a[isec][istrip][icsi]);
			    }
			  evento.epartlab.push_back(expevent.rco_epart[isec][istrip][icsi]);
			  float theout,phiout;
			  Classe_geo::Getgeo()->SpalmaRCO(isec,istrip,icsi,&theout,&phiout);
			  evento.thetalab.push_back(theout);
			  evento.philab.push_back(phiout);
			  
			  evento.rcoqf.push_back(expevent.rco_qf[isec][istrip][icsi]);
			  evento.coderiv.push_back(1000000+(isec+1)*100+(istrip+1)*10+icsi+1);
			  evento.rcocode.push_back(expevent.rco_code[isec][istrip][icsi]);
			  
		  float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
	
		 
		  evento.vpartlab.push_back(vv);
			evento.tvolo.push_back(Classe_geo::Getgeo()->rco.gas_dist[isec]/vv);
			evento.tbunch.push_back(expevent.rco_tvolo[isec][istrip][icsi]-200.*(float)(evento.bunch+1));
		  float vl[3];
		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();
		  float vc[3];
	
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);

     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
 evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));
 
		  //cose che non vengono settate per RCO exp
 evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);
evento.tplast.push_back(-1);
evento.garfqf.push_back(-1);
evento.moltepl++;
			    }//non accetta==0			  
			}//z identificato

 if(TMath::Nint(expevent.rco_z[isec][istrip][icsi])==0 &&TMath::Nint(expevent.rco_a[isec][istrip][icsi])==0 && expevent.rco_epart[isec][istrip][icsi]>0)// evento gamma in garfield
		{
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
	      evento.z.push_back(0.);
	      evento.a.push_back(0.);
	      
	      float theout,phiout;
	      Classe_geo::Getgeo()->SpalmaRCO(isec,istrip,icsi,&theout,&phiout);
	      
	      evento.thetalab.push_back(theout);
	      evento.philab.push_back(phiout);
	      evento.tvolo.push_back(-1);
		  evento.tbunch.push_back(expevent.rco_tvolo[isec][istrip][icsi]-200.*(float)(evento.bunch+1));
	      evento.epartcm.push_back(-1);
	      evento.thetacm.push_back(-1);
	      evento.phicm.push_back(-1);
	      evento.vpartcm.push_back(-1);
	      evento.vpartcm_x.push_back(-1);
	      evento.vpartcm_y.push_back(-1);
	      evento.vpartcm_z.push_back(-1);

	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
	      evento.epartlab.push_back(expevent.rco_epart[isec][istrip][icsi]);
	      evento.vpartlab.push_back(-1);
	      evento.vpartlab_x.push_back(-1);
	      evento.vpartlab_y.push_back(-1);
	      evento.vpartlab_z.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vplab.push_back(vjl);
		vjl.clear();
	      evento.indice_originale.push_back(-1);
	      evento.codphos.push_back(-1);
	      evento.phosqf.push_back(-1);
	      evento.phosga.push_back(-1);
	      
	      evento.rcoqf.push_back(expevent.rco_qf[isec][istrip][icsi]);
	      evento.rcocode.push_back(expevent.rco_code[isec][istrip][icsi]);
	      evento.garfqf.push_back(-1);
	      evento.tplast.push_back(-1);

	      evento.coderiv.push_back(1000000+(isec+1)*100+(istrip+1)*10+icsi+1);





	      evento.moltepl++;


		}//evento gamma in RCO


		    }//cesio
		}//strip
	    }//settore
	}//thereisRCO==1

      if(Classe_geo::Getgeo()->ThereIsHector==1)
	{
       float ddout,tthe,pphi;

      for(int ip=0;ip<8;ip++)
	{
       float tbaf1=0;
       float tbaf2=2;
       //if(ip==1){tbaf1=32.5;tbaf2=34.5;}
       //  if(ip==2){tbaf1=-32.;tbaf2=-30.;}
       //  if(ip==3){tbaf1=22.;tbaf2=24.;}
       //  if(ip==5){tbaf1=38.5;tbaf2=40.5;}
	  //	  if(exphector.ebaf[ip]>0 && exphector.tbaf[ip]>0)
	  	  if(exphector.ebaf[ip]>0 && exphector.tbaf[ip]<=tbaf2 &&exphector.tbaf[ip]>=tbaf1)
		    //if(exphector.ebaf[ip]>0)
	    {
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
	      evento.z.push_back(0.);
	      evento.a.push_back(0.);
	      Classe_geo::Getgeo()->Spalma_Hector(ip+1,&tthe,&pphi,&ddout);
	      evento.thetalab.push_back(tthe);
	      evento.philab.push_back(pphi);
	      evento.tvolo.push_back(exphector.tbaf[ip]);
		  evento.tbunch.push_back(-1);
	      evento.epartcm.push_back(-1);
	      evento.thetacm.push_back(-1);
	      evento.phicm.push_back(-1);
	      evento.vpartcm.push_back(-1);
	      evento.vpartcm_x.push_back(-1);
	      evento.vpartcm_y.push_back(-1);
	      evento.vpartcm_z.push_back(-1);

	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
	      evento.epartlab.push_back(exphector.ebaf[ip]/1000.);
	      evento.vpartlab.push_back(-1);
	      evento.vpartlab_x.push_back(-1);
	      evento.vpartlab_y.push_back(-1);
	      evento.vpartlab_z.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vplab.push_back(vjl);
		vjl.clear();
	      evento.indice_originale.push_back(-1);
	      evento.codphos.push_back(-1);
	      evento.phosqf.push_back(-1);
	      evento.phosga.push_back(-1);
	      evento.garfqf.push_back(-1);
	      evento.rcoqf.push_back(-1);
	      evento.rcocode.push_back(-1);

	      evento.tplast.push_back(-1);
	      int codeb=10000+ip+1;
	      evento.coderiv.push_back(codeb);




	      evento.moltepl++;
	      // printf("in strut %d %f %f\n",ip+1,evento.epartlab[evento.moltepl-1],evento.tvolo[evento.moltepl-1]);
	    }
	}//fine baf
	}//thereishector




      if(expevent.trig==Classe_analisi::Getanalisi()->triggerplastichino.Trigplast)//evento di plastichino
	{
	  evento.run=expevent.run;	
evento.tplast.push_back(expevent.tplast);
 evento.tvolo.push_back(-1);
 evento.tbunch.push_back(-1);
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
 evento.z.push_back(-1);
 evento.a.push_back(-1);
 evento.epartlab.push_back(-1);
 evento.thetalab.push_back(-1);
 evento.garfqf.push_back(-1);
 evento.rcoqf.push_back(-1);
 evento.rcocode.push_back(-1);

 evento.coderiv.push_back(9999);

 evento.philab.push_back(-1);
 evento.vpartlab.push_back(-1);
 evento.vpartlab_x.push_back(-1);
 evento.vpartlab_y.push_back(-1);
 evento.vpartlab_z.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vplab.push_back(vjl);
		vjl.clear();
 evento.vpartcm.push_back(-1);
 evento.thetacm.push_back(-1);
 evento.phicm.push_back(-1);
 evento.vpartcm_x.push_back(-1);
 evento.vpartcm_y.push_back(-1);
 evento.vpartcm_z.push_back(-1);

	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      evento.vpcm.push_back(vjl);
	vjl.clear();
	      
 evento.epartcm.push_back(-1);
evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);


evento.moltepl++;
	}

    }//if(tipo_analisi>=200&&tipo_analisi<220)

  //ntupla di BOLOGNA***************
 if(tipo_analisi==220)//EXP Ntupla di Bologna
    {
      
      evento.trig=expevent.trig;
      evento.run=-1;
	
      if(Classe_geo::Getgeo()->ThereIsGarf==1)
	{
      for(int isec=0;isec<24;isec++)
	{
	  for(int icsi=0;icsi<8;icsi++)
	    {
	      if(expevent.garf_z[isec][icsi]>0.&&expevent.garf_a[isec][icsi]>0.&&expevent.garf_epart[isec][icsi]>0)
		{
		  int non_accetta=0;
		  


		  if(non_accetta==0)

		    {

	evento.isformix.push_back(-1); //inserito per mixatore 29-01-21	      	
		  evento.z.push_back(expevent.garf_z[isec][icsi]);
		  evento.a.push_back(expevent.garf_a[isec][icsi]);
		  evento.epartlab.push_back(expevent.garf_epart[isec][icsi]);
		  
		  int code_micro=TMath::Nint(1+gRandom->Rndm());
		   float thetout,phiout;
		   int codice=(isec+1)*10+icsi+1;
		 Classe_geo::Getgeo()->Spalma_Garfield(codice,&thetout,&phiout,code_micro);
		

		
		      evento.thetalab.push_back(thetout);//prende quello spalmato
		      evento.philab.push_back(phiout);//prende quello spalmato
		 
		  evento.garfqf.push_back(expevent.garf_qf[isec][icsi]);
		  evento.rcoqf.push_back(-1);
	      evento.rcocode.push_back(-1);

		  evento.coderiv.push_back((isec+1)*10+icsi+1+1000);
		  
		  

		  float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
		  
		 
		  evento.vpartlab.push_back(vv);
		evento.tvolo.push_back(Classe_geo::Getgeo()->garf.dist[icsi]/vv);
		evento.tbunch.push_back(-1);
		  float vl[3];
		  
		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
		  
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();
		  float vc[3];
		  
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);
     
     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
 evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));
 
		  //cose che non vengono settate per garfield exp
 evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);
		  evento.tplast.push_back(-1);

 

evento.moltepl++;
		    }//non_accetta=0
		}//evento buono

	    }
	}//fine garf
	}//ThereIsGarf==1
      if(Classe_geo::Getgeo()->ThereIsRCO==1)
	{
	  int non_accetta;
	  for(int isec=0;isec<8;isec++)
	    {
	      for(int istrip=0;istrip<9;istrip++)
		{
		  for(int icsi=0;icsi<7;icsi++)
		    {
		      non_accetta=1;

			  		      if(expevent.rco_z[isec][istrip][icsi]>0. &&expevent.rco_a[isec][istrip][icsi]>0 && expevent.rco_epart[isec][istrip][icsi]>0)//Se esiste uno Z identificato
			    	{
			  non_accetta=0;
				}
			  if(non_accetta==0)
			    {
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
			  evento.z.push_back(expevent.rco_z[isec][istrip][icsi]);

			  evento.a.push_back(expevent.rco_a[isec][istrip][icsi]);
			    
			  evento.epartlab.push_back(expevent.rco_epart[isec][istrip][icsi]);
			  float theout,phiout;
			  Classe_geo::Getgeo()->SpalmaRCO(isec,istrip,icsi,&theout,&phiout);
			  evento.thetalab.push_back(theout);
			  evento.philab.push_back(phiout);
			  
			  evento.rcoqf.push_back(expevent.rco_qf[isec][istrip][icsi]);
			  evento.coderiv.push_back(1000000+(isec+1)*100+(istrip+1)*10+icsi+1);
			  evento.rcocode.push_back(-1);
			  
		  float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
	
		 
		  evento.vpartlab.push_back(vv);
			evento.tvolo.push_back(Classe_geo::Getgeo()->rco.gas_dist[isec]/vv);
			evento.tbunch.push_back(-1);
		  float vl[3];
		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();
		  float vc[3];
	
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);

     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	vjl.clear();
	      
 evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));
 
		  //cose che non vengono settate per RCO exp
 evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);
evento.tplast.push_back(-1);
evento.garfqf.push_back(-1);
evento.moltepl++;
			    }//non accetta==0			  
		
		    }//cesio
		}//strip
	    }//settore
	}//thereisRCO==1

     


      if(expevent.trig==Classe_analisi::Getanalisi()->triggerplastichino.Trigplast)//evento di plastichino
	{
evento.run=expevent.run;	
evento.tplast.push_back(expevent.tplast);
 evento.tvolo.push_back(-1);
 evento.tbunch.push_back(-1);
evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
 evento.z.push_back(-1);
 evento.a.push_back(-1);
 evento.epartlab.push_back(-1);
 evento.thetalab.push_back(-1);
 evento.garfqf.push_back(-1);
 evento.rcoqf.push_back(-1);
 evento.rcocode.push_back(-1);

 evento.coderiv.push_back(9999);

 evento.philab.push_back(-1);
 evento.vpartlab.push_back(-1);
 evento.vpartlab_x.push_back(-1);
 evento.vpartlab_y.push_back(-1);
 evento.vpartlab_z.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vplab.push_back(vjl);
		vjl.clear();
 evento.vpartcm.push_back(-1);
 evento.thetacm.push_back(-1);
 evento.phicm.push_back(-1);
 evento.vpartcm_x.push_back(-1);
 evento.vpartcm_y.push_back(-1);
 evento.vpartcm_z.push_back(-1);

	      vjl.push_back(-1);
	      vjl.push_back(-1);
	      vjl.push_back(-1);
		evento.vpcm.push_back(vjl);
	vjl.clear();

 evento.epartcm.push_back(-1);
evento.indice_originale.push_back(-1);
 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);


evento.moltepl++;
	}

    }//if(tipo_analisi==220)


  //ntupla di KALIVEDA per ISOFAZIA***************
 if(tipo_analisi==230)//EXP ntupla KALIVEDA per ISOFAZIA
    {
      
      evento.trig=expevent.trig;
      evento.run=faziakali.run;
	  float atloc=0;
			  float epost,elost;
			  int mate=1;
			  int idir=1;
			  int icod;
			  float dummy_pressione=0;
      if(Classe_geo::Getgeo()->ThereIsBlocchiFazia==1 && faziakali.mtot>0)
	{
	  int scarta=0;
	  for(int j=0;j<faziakali.mtot;j++)
	    {


	      scarta=1;
	      	      if((faziakali.z[j]<=2 && faziakali.idcode[j]<3 && faziakali.ecode[j]<3)||(faziakali.z[j]>2 && faziakali.idcode[j]<3 && faziakali.ecode[j]!=3)) //da faziasym, 16feb2021, condizione meno restrittiva
			//	if(faziakali.idcode[j]==0 && faziakali.ecode[j]>=0&&faziakali.ecode[j]<3)//identificazione buona e energia buona, condizione fino a 16feb2021
		{
		  scarta=0;
		  // cout<<faziakali.idtype[j]<<endl;
		  //	  if(Classe_geo::Getgeo()->fazietto.alimSi1Si2[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]>=0&&Classe_geo::Getgeo()->fazietto.alimSi2CsI[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]>=0&&Classe_geo::Getgeo()->fazietto.aliminfSi1PSA[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]>=0&&Classe_geo::Getgeo()->fazietto.alimsupSi1PSA[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]>=0)//ci sono tutte le griglie
		  //	    {
		      if(Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]>0)//rivelatori funzionanti
			{
			  float zz=faziakali.z[j];
			  float aa=faziakali.a[j];
			  //		       if(faziakali.esi2[j]>2)
			  if((Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&4)&&(faziakali.ecsi[j]>10))
 			{
		

			  float spe=Classe_geo::Getgeo()->fazietto.spes_si1[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1];
			  Classe_geo::Getgeo()->ecorr_veda(&faziakali.etot[j],&zz,&aa,&atloc,&epost,&elost,&mate,&spe,&idir,&icod,&dummy_pressione);
			  

			  if(faziakali.esi1[j]-elost>50)
			    {
			      // cout<<"el "<<faziakali.esi1[j]<<" "<<elost<<" "<<faziakali.etot[j]<<endl;
			      scarta=1;
			    }
 			}
			  if((faziakali.ecsi[j]>10)&&(Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&4)&&(Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&2))
 			{
		

			  float spe=Classe_geo::Getgeo()->fazietto.spes_si1[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]+Classe_geo::Getgeo()->fazietto.spes_si2[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1];
			  Classe_geo::Getgeo()->ecorr_veda(&faziakali.etot[j],&zz,&aa,&atloc,&epost,&elost,&mate,&spe,&idir,&icod,&dummy_pressione);
			  // cout<<"spe="<<Classe_geo::Getgeo()->fazietto.spes_si1[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]<<" "<<Classe_geo::Getgeo()->fazietto.spes_si2[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]<<" "<<spe<<endl;
			  // cout<<elost<<" "<<faziakali.esi1[j]+faziakali.esi2[j]<<" "<<faziakali.esi1[j]+faziakali.esi2[j]-elost<<" "<<spe<<endl;
			  if((faziakali.esi1[j]+faziakali.esi2[j])-elost>50)
			    {
			      //cout<<"el2"<<" "<<faziakali.esi1[j]<<" "<<elost<<" "<<faziakali.etot[j]<<endl;
			      scarta=1;
			    }
 			}

			  if((Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&4)==0 &&  (faziakali.idtype[j]==11 || faziakali.idtype[j]==12))
			    {
			      scarta=1;
			    }
			  if((Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&2)==0 &&  (faziakali.idtype[j]==12 || faziakali.idtype[j]==23))
			    {
			      scarta=1;
			    }
			  if((Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]&1)==0 &&  (faziakali.idtype[j]==33 || faziakali.idtype[j]==23))
			    {
			      scarta=1;
			    }
			  // cout<<Classe_geo::Getgeo()->fazietto.rotti[faziakali.blocco[j]][faziakali.qua[j]-1][faziakali.tel[j]-1]<<" "<<faziakali.idtype[j]<<" "<<scarta<<endl;
			  if(faziakali.z[j]==4 &&faziakali.a[j]==8)// si butta il Be8
			    {
			      scarta=1;
			    }
			
		       if(scarta==0)
			 {

			   //   if(Classe_analisi::Getanalisi()->reazione.zp==20 && Classe_analisi::Getanalisi()->reazione.ap==48 && TMath::Nint(Classe_analisi::Getanalisi()->reazione.ebeam)==40 && faziakali.z[j]>=17 && faziakali.z[j]<=19)
			   //	{
			   //faziakali.a[j]=faziakali.a[j]+1;
			   //	}



evento.isformix.push_back(-1); //inserito per mixatore 29-01-21
		  evento.z.push_back((float)faziakali.z[j]);
		  //		  if(faziakali.aid[j]==1)
		  //  {
		  evento.a.push_back((float)faziakali.a[j]);
		  // }
		  // else
		  // {
		  // evento.a.push_back(Classe_formule::QualeA((float)faziakali.z[j]));
		      // }
		  evento.indice_originale.push_back(j);//posizione della particella in faziakali
		  evento.rcoqf.push_back(faziakali.aid[j]);//1 se la massa e' misurata, 0 se e' la beta stability
		  evento.epartlab.push_back(faziakali.etot[j]);
		  

		  float dout=Classe_geo::Getgeo()->fazietto.dist[faziakali.blocco[j]];
	     float theout,phiout; 
	     int code=faziakali.blocco[j]*100+(faziakali.qua[j]-1)*4+faziakali.tel[j]-1;
	     evento.coderiv.push_back(code+10000000);
	     // Classe_geo::Getgeo()->Spalma_Fazietto(code,&theout,&phiout,&dout);//per spalmare
	      theout=faziakali.theta[j];
	      phiout=faziakali.phi[j];
	      evento.thetalab.push_back(theout);
	      evento.philab.push_back(phiout);
	      evento.rcocode.push_back(faziakali.idtype[j]);// come e' stata identificata la particella: 11 PSA in Si1; 12 DE-E Si1-Si2; 22 PSA Si2 (solo se 12 e' impossibile, es. manca Si1); 23 De-E Si2-CsI; 33 fast slow in CsI (solo se 23 non e' possibile, es particella non vista in Si o Si2 rotto)
float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
evento.vpartlab.push_back(vv);

 evento.tvolo.push_back(dout/vv);
float vl[3];
		  Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
		  evento.vpartlab_x.push_back(vl[0]);
		  evento.vpartlab_y.push_back(vl[1]);
		  evento.vpartlab_z.push_back(vl[2]);
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
		evento.vplab.push_back(vjl);
		vjl.clear();
		  float vc[3];
	
     Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);

     float VC,THEC,PHIC;
     Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
     evento.vpartcm.push_back(VC);
     evento.thetacm.push_back(THEC);
     evento.phicm.push_back(PHIC);
     evento.vpartcm_x.push_back(vc[0]);
     evento.vpartcm_y.push_back(vc[1]);
     evento.vpartcm_z.push_back(vc[2]);

	      vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
		evento.vpcm.push_back(vjl);
	
		vjl.clear();
	
	
 evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));

		  //cose che non vengono settate per FAZIETTO exp

 evento.phosga.push_back(-1);
 evento.codphos.push_back(-1);
 evento.phosqf.push_back(-1);
evento.tplast.push_back(-1);
evento.garfqf.push_back(-1);

 evento.tbunch.push_back(-1);


		  evento.moltepl++;
			 }//scarta==0
			}//i rivelatori sono funzionanti (da tel rotti)
		      // }//ci sono tutte le griglie
		}
	    }

	}//thereisblocchifazia
    }// fine tipo_analisi==230 ntupla di KALIVEDA per ISOFAZIA


 //-----------------------------------------ntupla per INDRAFAZIA-----------------
 if(tipo_analisi==240)//EXP ntupla KALIVEDA per INDRAFAZIA
    {
      
      evento.trig=expevent.trig;
      evento.run=indrafazia.run;
      if(indrafazia.mtot>0)
       {
        int scarta = 0;  // lo faccio diventare 1 per scartare roba (particelle) che non voglio
        for(int j=0; j<indrafazia.mtot;j++)
         {
          scarta=0;
          if(indrafazia.idcode[j]==0) //controllo il successo dell'identificazione
          {
           //if(indrafazia.ecode[j]>0 && indrafazia.ecode[j]<3) //controllo il successo della calibrazione
           if(true) //per ora prendo anche le cose non cali
           {
            float zz=indrafazia.z[j];
	    float aa=indrafazia.a[j];
            if(indrafazia.z[j]>indrafazia.a[j]) scarta=1;
	    //if(indrafazia.z[j]==4 && indrafazia.a[j]==8) scarta=1;// ad esempio, se si vuole buttare il Be8... 
	    //if(indrafazia.array[j]==1 && (indrafazia.gt_dt[j]<64 ||indrafazia.gt_dt[j]>72)) scarta=1;

            //if(indrafazia.array[j]==1 && indrafazia.z[j]<6 && indrafazia.gt_dt[j]>82) scarta=1;
            //if(indrafazia.array[j]==1 && (indrafazia.z[j]>=6 && indrafazia.z[j]<9) && (indrafazia.gt_dt[j]<40 ||indrafazia.gt_dt[j]>82)) scarta=1;
            //if(indrafazia.array[j]==1 && (indrafazia.z[j]>=9 && indrafazia.z[j]<20) && (indrafazia.gt_dt[j]<66 ||indrafazia.gt_dt[j]>83)) scarta=1;
            //if(indrafazia.array[j]==1 && (indrafazia.z[j]>=20 && indrafazia.z[j]<=24) && (indrafazia.gt_dt[j]<71 ||indrafazia.gt_dt[j]>83)) scarta=1;
	    //if(indrafazia.array[j]==1 && indrafazia.z[j]>24 && (indrafazia.gt_dt[j]<75 ||indrafazia.gt_dt[j]>81)) scarta=1;


	    if(scarta==0)
	     {
	      evento.indice_originale.push_back(j);//posizione della particella in indrafazia
              evento.z.push_back((float)indrafazia.z[j]);
	      evento.a.push_back((float)indrafazia.a[j]);
	      evento.rcoqf.push_back(indrafazia.aid[j]);//1 se la massa e' misurata, 0 se e' la beta stability
	      evento.epartlab.push_back(indrafazia.etot[j]);
	      float theout,phiout; 
	      theout=indrafazia.theta[j];
	      phiout=indrafazia.phi[j];
	      evento.thetalab.push_back(theout);
	      evento.philab.push_back(phiout);
	      float dout=100.; //----------------------------------------------------------ATTENZIONE: DISTANZA DEL RIVELATORE. SERVE SOLO PER IL TEMPO DI VOLO, CHE ORA Ãˆ ERRATO!!!!!!! 
	      float vv=Classe_formule::e2v(evento.epartlab[evento.moltepl],evento.a[evento.moltepl]);
              evento.vpartlab.push_back(vv);
              if(vv!=0)evento.tvolo.push_back(dout/vv);
              float vl[3];
	      Classe_formule::polcar(evento.vpartlab[evento.moltepl],evento.thetalab[evento.moltepl],evento.philab[evento.moltepl],vl);
	      evento.vpartlab_x.push_back(vl[0]);
	      evento.vpartlab_y.push_back(vl[1]);
	      evento.vpartlab_z.push_back(vl[2]);  
	      vjl.push_back(vl[0]);
	      vjl.push_back(vl[1]);
	      vjl.push_back(vl[2]);
	      evento.vplab.push_back(vjl);
	      vjl.clear();        
              float vc[3];
              Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);
              float VC,THEC,PHIC;
              Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
              evento.vpartcm.push_back(VC);
              evento.thetacm.push_back(THEC);
              evento.phicm.push_back(PHIC);
              evento.vpartcm_x.push_back(vc[0]);
              evento.vpartcm_y.push_back(vc[1]);
              evento.vpartcm_z.push_back(vc[2]);
              vjl.push_back(vc[0]);
	      vjl.push_back(vc[1]);
	      vjl.push_back(vc[2]);
	      evento.vpcm.push_back(vjl);
	      vjl.clear();
              evento.epartcm.push_back(Classe_formule::v2e(evento.vpartcm[evento.moltepl],evento.a[evento.moltepl]));              
              
	      evento.rcocode.push_back(indrafazia.idtype[j]);
	      
	      evento.codphos.push_back(indrafazia.array[j]);
	      int code=0;
	      if(evento.codphos[evento.moltepl]==1) code = indrafazia.idtel[j]; //fazia
	      if(evento.codphos[evento.moltepl]==0) code = indrafazia.ring[j]*100 + indrafazia.module[j];//indra
	      evento.coderiv.push_back(code); 
              evento.tbunch.push_back((float)indrafazia.gt_dt[j]);
              evento.phosqf.push_back(indrafazia.ecode[j]);
		 //cose che non vengono settate per INDRAFAZIA exp
                 evento.phosga.push_back(-1);
                 evento.tplast.push_back(-1);
                 evento.garfqf.push_back(-1);
		 //   evento.tbunch.push_back(-1);
              evento.moltepl++;
	     }//scarta==0
		    
           } //controllo ecode
          } //controllo idcode
         } //fine giro su molt
       } //chiude if molt>0
       	
    }// fine tipo_analisi==240 ntupla di KALIVEDA per INDRAFAZIA




  
  if(evento.moltepl>0)
    {
  if(Classe_analisi::Getanalisi()->ntupl==1)
    {
      Classe_analisi::Getanalisi()->tree.moltepl=evento.moltepl;
        for(int jk=0;jk<evento.moltepl;jk++)
      	{
      	  Classe_analisi::Getanalisi()->tree.z[jk]=evento.z[jk];
      	  Classe_analisi::Getanalisi()->tree.a[jk]=evento.a[jk];
      	  Classe_analisi::Getanalisi()->tree.vxcm[jk]=evento.vpartcm_x[jk];
      	  Classe_analisi::Getanalisi()->tree.vycm[jk]=evento.vpartcm_y[jk];
      	  Classe_analisi::Getanalisi()->tree.vzcm[jk]=evento.vpartcm_z[jk];

	  //      	  Classe_analisi::Getanalisi()->tree.vxlab[jk]=evento.vpartlab_x[jk];
      	  //Classe_analisi::Getanalisi()->tree.vylab[jk]=evento.vpartlab_y[jk];
      	  //Classe_analisi::Getanalisi()->tree.vzlab[jk]=evento.vpartlab_z[jk];

      	  Classe_analisi::Getanalisi()->tree.vlabmod[jk]=evento.vpartlab[jk];
      	  Classe_analisi::Getanalisi()->tree.thelab[jk]=evento.thetalab[jk];
      	  Classe_analisi::Getanalisi()->tree.philab[jk]=evento.philab[jk];

	  //      	  Classe_analisi::Getanalisi()->tree.vcmmod[jk]=evento.vpartcm[jk];
	  //              Classe_analisi::Getanalisi()->tree.thecm[jk]=evento.thetacm[jk];
      	  //                  Classe_analisi::Getanalisi()->tree.phicm[jk]=evento.phicm[jk];



      	}

 


  Classe_analisi::Getanalisi()->ntupla->Fill();
    }
    }
}//Fine copia in evento

//ga0=offset _1 box 1 [phos]
//ga1=slope _1 box 1 [phos]
static float ga0_1[9]={999,-209.51,-746.39,-399.33,-494.82,-1264.96,-798.07,-750.72,-1055.39};
static float ga1_1[9]={999,0.79,1.80,0.92,1.32,2.88,1.86,1.83,2.53};

static float ga0_2[9]={999,-319.03,-1059.92,-834.14,-306.84,-343.84,-727.84,-405.42,-677.52};
static float ga1_2[9]={999,0.83,2.42,1.83,0.78,0.82,1.77,1.22,1.51};

static float ga0_3[9]={999,-198.79,-669.77,-503.13,-1246.79,999,-421.31,-369.57,-443.01};
static float ga1_3[9]={999,0.71,1.63,1.31,2.79,999,0.97,1.10,1.00};

static float ga0_4[9]={999,-567.44,-408.06,-224.73,-869.92,-1386.17,-1488.01,-673.21,-1314.44};
static float ga1_4[9]={999,1.47,0.96,0.57,2.03,3.12,3.35,1.44,3.01};

static float ga0_5[9]={-192.35,-276.81,-98.81,999,999,999,999,999,999};
static float ga1_5[9]={0.75,0.92,0.54};

static float ga0_6[9]={-1365.38,-172.21,999,999,999,-2132.36,999,999,999};
static float ga1_6[9]={3.05,0.72,999,999,999,4.58,999,999,999};

//protoni box phos [7 parametri]; alpha box phos [7 parametri]

static float p1_2[8]={-4.0040349121e+03,6.6696594238e+01,-4.5359614491e-01,1.6307237092e-03,-3.2617094803e-06,3.4471252519e-09,-1.5059508545e-12,0.0000000000e+00};
static float a1_2[8]={-3.2560925293e+02,3.1972057819e+00,-1.0083393194e-02,1.4533080503e-05,-6.6153207356e-09,-5.9601074076e-12,7.7676487272e-15,-2.4189081353e-18};
static float p1_3[8]={-1.2916474609e+03,1.7531852722e+01,-9.4318613410e-02,2.6687918580e-04,-4.1534855200e-07,3.3752112216e-10,-1.1209232192e-13,0.0000000000e+00};
static float a1_3[8]={-2.6929119873e+02,2.7354657650e+00,-1.0626775213e-02,2.3766764571e-05,-3.1534884926e-08,2.4340286731e-11,-1.0053735595e-14,1.7145215541e-18};
static float p1_4[8]={-6.3446733398e+03,2.4734960938e+02,-4.0295138359e+00,3.5816133022e-02,-1.8746456772e-04,5.7860637526e-07,-9.7658492404e-10,6.9629057294e-13};
static float a1_4[8]={-1.8662678528e+02,3.2524726391e+00,-1.6781792045e-02,4.3474250560e-05,-5.5960356349e-08,2.7304561390e-11,7.4551003459e-15,-8.6413606447e-18};
static float p1_5[8]={-9.7997187500e+03,1.7856428528e+02,-1.3321653605e+00,5.2235778421e-03,-1.1336226635e-05,1.2919074699e-08,-6.0468539897e-12,0.0000000000e+00};
static float a1_5[8]={3.6187835693e+01,-1.0334856510e+00,8.9321872219e-03,-2.9133034332e-05,4.8683016729e-08,-4.4309812763e-11,2.0955417238e-14,-4.0416363194e-18};
static float p1_6[8]={-2.8530798340e+03,4.4555339813e+01,-2.8196987510e-01,9.3929102877e-04,-1.7291735048e-06,1.6713020967e-09,-6.6451656460e-13,0.0000000000e+00};
static float a1_6[8]={-4.6038104248e+02,5.3241400719e+00,-2.5350714102e-02,7.0610920375e-05,-1.1987451387e-07,1.2127621130e-10,-6.6839634142e-14,1.5395287037e-17};
static float p1_7[8]={-1.3371432495e+02,-5.8384245634e-01,2.1829849109e-02,-1.2725990382e-04,3.4643332469e-07,-4.9923121193e-10,3.6909890292e-13,-1.1047957272e-16};
static float a1_7[8]={-3.7500033569e+02,3.6669754982e+00,-1.3995140791e-02,2.9793804060e-05,-3.7152481980e-08,2.6900266736e-11,-1.0461348178e-14,1.6892679368e-18};
static float p1_8[8]={7.3878271484e+02,-2.5100494385e+01,2.9384341836e-01,-1.6938686604e-03,5.4665979405e-06,-1.0078491819e-08,9.9405171594e-12,-4.0760013163e-15};
static float a1_8[8]={-1.4839857483e+02,1.2449846268e+00,-2.8714104556e-03,3.2484269923e-06,-1.7525767504e-09,3.6049470158e-13,0.0000000000e+00,0.0000000000e+00};
static float p1_9[8]={-2.2164562988e+03,2.8144714355e+01,-1.4365707338e-01,3.8563759881e-04,-5.7177186363e-07,4.4458073334e-10,-1.4184989717e-13,0.0000000000e+00};
static float a1_9[8]={-9.5876571655e+01,3.0784016848e-01,1.6878369497e-03,-7.0700266406e-06,1.0970405917e-08,-8.5399456604e-12,3.3421904160e-15,-5.2473546326e-19};
static float p2_2[8]={-3.9262687988e+03,1.1797425079e+02,-1.4687540531e+00,9.9133299664e-03,-3.9050755731e-05,8.9898385625e-08,-1.1222315438e-10,5.8725431202e-14};
static float a2_2[8]={-1.1239528656e+02,1.4112386703e+00,-1.2180192862e-03,-2.1090592782e-05,9.1014342729e-08,-1.6050615004e-10,1.3411449414e-13,-4.3791269192e-17};
static float p2_3[8]={-1.1680218750e+04,1.9386967468e+02,-1.3564451933e+00,5.2039045841e-03,-1.1813462152e-05,1.5875890824e-08,-1.1700802653e-11,3.6500585167e-15};
static float a2_3[8]={-8.2155950928e+02,8.2457876205e+00,-3.4489572048e-02,8.1780017354e-05,-1.1615338735e-07,9.7588000181e-11,-4.4616718496e-14,8.5424993266e-18};
static float p2_4[8]={-6.5386254883e+02,1.3450217247e+01,-1.0646681488e-01,4.5471161138e-04,-1.1208626347e-06,1.6008744330e-09,-1.2324017674e-12,3.9606450638e-16};
static float a2_4[8]={-6.8391441345e+01,9.1534274817e-01,-2.9261705931e-03,6.0324209699e-06,-8.4645614962e-09,7.5601330751e-12,-3.7407393230e-15,7.6771542550e-19};
static float p2_5[8]={-1.1645745117e+04,4.0500595093e+02,-5.8929586411e+00,4.6579383314e-02,-2.1551860846e-04,5.8340344822e-07,-8.5517865012e-10,5.2329603324e-13};
static float a2_5[8]={-3.9623370361e+02,6.5204939842e+00,-3.7425398827e-02,1.1112225911e-04,-1.7898817362e-07,1.4935679082e-10,-5.0736639960e-14,0.0000000000e+00};
static float p2_6[8]={-2.9668454590e+03,1.0126020050e+02,-1.4277074337e+00,1.0949739255e-02,-4.9219343055e-05,1.2987905507e-07,-1.8665397006e-10,1.1289673216e-13};
static float a2_6[8]={-1.7620214844e+02,2.8545460701e+00,-1.3641007245e-02,3.3513566450e-05,-4.4249119924e-08,3.0019053215e-11,-8.2288700254e-15,0.0000000000e+00};
static float p2_7[8]={-1.6630506592e+03,2.3203887939e+01,-1.2148306519e-01,2.8242517146e-04,-1.4698615303e-07,-5.7649407470e-10,1.0812458784e-12,-5.7100121513e-16};
static float a2_7[8]={-1.2400426483e+02,1.0057499409e+00,-2.1263370290e-03,2.2335136691e-06,-1.1420904222e-09,2.3062083384e-13,0.0000000000e+00,0.0000000000e+00};
static float p2_8[8]={-5.0777519531e+03,7.5092536926e+01,-4.6310895681e-01,1.5572846169e-03,-3.0817545849e-06,3.5926124298e-09,-2.2825352719e-12,6.0744819582e-16};
static float a2_8[8]={-2.7033358765e+02,1.4245300293e+00,6.9537211675e-04,-1.3951862456e-05,3.1963796943e-08,-3.3267753136e-11,1.6888483410e-14,-3.3945147282e-18};
static float p2_9[8]={-4.4696865234e+03,1.3881864929e+02,-1.7745459080e+00,1.2124995701e-02,-4.7189889301e-05,1.0289155483e-07,-1.1259453786e-10,4.3747134143e-14};
static float a2_9[8]={-1.2293434143e+02,1.3154457808e+00,6.1033701058e-04,-3.0247651011e-05,1.1224834395e-07,-1.8536618074e-10,1.4784896397e-13,-4.6387418558e-17};
static float p3_2[8]={-5.0800698242e+03,9.7213661194e+01,-7.7637398243e-01,3.3817607909e-03,-8.6758609541e-06,1.3126649989e-08,-1.0855804970e-11,3.7884662414e-15};
static float a3_2[8]={-1.8541755676e+02,1.3019353151e+00,-8.0412556417e-04,-7.7943859651e-06,2.2021254154e-08,-2.5387967770e-11,1.3951230647e-14,-3.0149598022e-18};
static float p3_3[8]={-7.6833120117e+03,1.3338256836e+02,-9.7074037790e-01,3.8621623535e-03,-9.0632802312e-06,1.2551217843e-08,-9.5018055923e-12,3.0343559849e-15};
static float a3_3[8]={-2.3726199341e+02,1.9116970301e+00,-5.0892247818e-03,7.5138964348e-06,-6.8085479477e-09,3.9463198334e-12,-1.3853705602e-15,2.2455001630e-19};
static float p3_4[8]={-2.9713525391e+03,5.2374820709e+01,-3.7506663799e-01,1.4152438380e-03,-2.9534719488e-06,3.2349649626e-09,-1.4556718461e-12,0.0000000000e+00};
static float a3_4[8]={-1.6189535522e+02,1.4831246138e+00,-3.8762583863e-03,5.0343205658e-06,-3.1893385710e-09,7.9303252333e-13,0.0000000000e+00,0.0000000000e+00};
static float p3_5[8]={-9.4561022949e+02,1.2084351540e+01,-5.9936791658e-02,1.5517199063e-04,-2.1809326256e-07,1.5766572770e-10,-4.5829144516e-14,0.0000000000e+00};
static float a3_5[8]={-4.2049597168e+02,4.1047601700e+00,-1.5770556405e-02,3.3847893064e-05,-4.2628908403e-08,3.1211297435e-11,-1.2281873550e-14,2.0074378929e-18};
static float p3_7[8]={-1.5744055176e+02,1.5818041563e+00,6.0510113835e-02,-1.1516949162e-03,8.6276686488e-06,-3.2764145175e-08,6.2664956557e-11,-4.8029556542e-14};
static float a3_7[8]={-7.3276832581e+01,1.5884759426e+00,-6.8838326260e-03,1.2870877072e-05,-3.9873282454e-09,-2.0031139941e-11,2.7554295176e-14,-1.1029402949e-17};
static float p3_8[8]={-1.8005867188e+04,3.8181851196e+02,-3.4218053818e+00,1.6805520281e-02,-4.8790261644e-05,8.3759211122e-08,-7.8783049529e-11,3.1345324963e-14};
static float a3_8[8]={-7.9168098450e+01,5.0873268396e-02,6.0779163614e-03,-3.0868432077e-05,6.9698693039e-08,-8.2497529641e-11,4.9947008741e-14,-1.2229465828e-17};
static float p3_9[8]={-6.5716296387e+02,2.1180049896e+01,-2.3424954712e-01,9.9341501482e-04,1.2928419437e-06,-2.6196193659e-08,8.4604032113e-11,-9.0666975880e-14};
static float a3_9[8]={-6.9682151794e+01,1.4015486240e+00,-4.4616563246e-03,-6.8938680897e-07,3.4626292944e-08,-7.8850995178e-11,7.3385640284e-14,-2.5411014887e-17};
static float p4_2[8]={-2.6565891113e+03,3.9216537476e+01,-2.3441153765e-01,7.3934294051e-04,-1.2932002846e-06,1.1904247588e-09,-4.5107084842e-13,0.0000000000e+00};
static float a4_2[8]={-1.4119215393e+02,1.1188602448e+00,-2.4509031791e-03,2.7181624773e-06,-1.5063978953e-09,3.3557241729e-13,0.0000000000e+00,0.0000000000e+00};
static float p4_3[8]={-3.6449926758e+03,1.2605455017e+02,-1.7714353800e+00,1.3079551049e-02,-5.3257950640e-05,1.1342685013e-07,-9.8845896745e-11,0.0000000000e+00};
static float a4_3[8]={-2.8516937256e+02,4.9398727417e+00,-2.8168287128e-02,8.2841223048e-05,-1.3156876832e-07,1.0787973292e-10,-3.5907095439e-14,0.0000000000e+00};
static float p4_4[8]={-3.3134838867e+03,1.2235340118e+02,-1.8276492357e+00,1.4325841330e-02,-6.1941049353e-05,1.4023646600e-07,-1.3011018651e-10,0.0000000000e+00};
static float a4_4[8]={-1.5988375854e+02,2.8553576469e+00,-1.4578028582e-02,3.8079746446e-05,-5.2932055183e-08,3.7476240861e-11,-1.0631796618e-14,0.0000000000e+00};
static float p4_5[8]={-5.5084228516e+02,7.6232442856e+00,-3.9970427752e-02,1.1097604147e-04,-1.6870615127e-07,1.3343631755e-10,-4.2970472693e-14,0.0000000000e+00};
static float a4_5[8]={-1.0177314758e+02,9.0369194746e-01,-2.0166968461e-03,2.3825666631e-06,-1.4540856297e-09,3.6359790504e-13,0.0000000000e+00,0.0000000000e+00};
static float p4_6[8]={-1.9292406250e+04,3.3384024048e+02,-2.4376173019e+00,9.7489971668e-03,-2.3049818992e-05,3.2232257752e-08,-2.4700743187e-11,8.0082687647e-15};
static float a4_6[8]={-1.8309881592e+02,1.2771540880e+00,-2.7786684223e-03,4.3285026550e-06,-6.6474807880e-09,7.8342896406e-12,-5.0069879341e-15,1.2565152678e-18};
static float p4_7[8]={-1.8412921143e+03,2.2004608154e+01,-1.0603181273e-01,2.7165596839e-04,-3.9235649751e-07,3.1123820188e-10,-1.1776468472e-13,1.2943761103e-17};
static float a4_7[8]={-1.4274264526e+02,8.8965970278e-01,-1.4694449492e-03,1.1801188293e-06,-4.5158815576e-10,6.5976920921e-14,0.0000000000e+00,0.0000000000e+00};
static float p4_8[8]={-2.9562646484e+02,7.7519502640e+00,-3.4535158426e-02,-5.0222617574e-04,6.9584812081e-06,-3.4738739885e-08,8.0252082757e-11,-7.1766987427e-14};
static float a4_8[8]={-6.8183746338e+01,1.4072570801e+00,-5.7451124303e-03,1.1562948202e-05,-1.1106941145e-08,4.1080329242e-12,0.0000000000e+00,0.0000000000e+00};
static float p4_9[8]={-1.4796660156e+03,1.9983949661e+01,-1.0701104999e-01,3.0022350256e-04,-4.6184547386e-07,3.7025818522e-10,-1.2122122967e-13,0.0000000000e+00};
static float a4_9[8]={-1.3549725342e+02,1.0708279610e+00,-2.2006253712e-03,2.1938940336e-06,-1.0366790759e-09,1.8645297166e-13,0.0000000000e+00,0.0000000000e+00};
static float p5_1[8]={-2.4554224609e+04,4.9399911499e+02,-4.1004753113e+00,1.7988244072e-02,-4.3941872718e-05,5.6681699334e-08,-3.0175782012e-11,0.0000000000e+00};
static float a5_1[8]={-8.0506964111e+02,9.4262208939e+00,-4.1530162096e-02,9.4891474873e-05,-1.1827143709e-07,7.6596971821e-11,-2.0214700478e-14,0.0000000000e+00};
static float p5_2[8]={-8.8243681641e+03,1.5747789001e+02,-1.1531716585e+00,4.4458983466e-03,-9.4915130830e-06,1.0639223724e-08,-4.8954846253e-12,0.0000000000e+00};
static float a5_2[8]={-3.0819207764e+02,2.8870427608e+00,-9.0159010142e-03,1.4307859601e-05,-1.1930595178e-08,4.9326059730e-12,-7.7800526777e-16,0.0000000000e+00};
static float p5_3[8]={-4.6519929688e+04,1.1207482910e+03,-1.1447717667e+01,6.4303115010e-02,-2.1451656357e-04,4.2524962396e-07,-4.6412523824e-10,2.1527888521e-13};
static float a5_3[8]={-5.6711853027e+02,6.6946702003e+00,-2.7789141983e-02,5.4892116168e-05,-4.5136008708e-08,-6.4443702111e-12,3.3032177525e-14,-1.4635349591e-17};
static float p6_1[8]={-5.0804394531e+02,4.7811341286e+00,-1.5493327752e-02,2.2063095457e-05,-8.0226492116e-09,-9.7796371695e-12,7.2493399494e-15,0.0000000000e+00};
static float a6_1[8]={-2.7245166016e+02,1.7862100601e+00,-4.0790708736e-03,4.7454032028e-06,-2.7148334691e-09,6.0153330884e-13,0.0000000000e+00,0.0000000000e+00};
static float p6_2[8]={-4.5153920898e+03,7.6649002075e+01,-5.2973771095e-01,1.9251039485e-03,-3.8719040276e-06,4.0917322863e-09,-1.7771822425e-12,0.0000000000e+00};
static float a6_2[8]={-3.9333441162e+02,4.0230326653e+00,-1.4329838566e-02,2.6473770049e-05,-2.6608010018e-08,1.3870537531e-11,-2.9444059563e-15,0.0000000000e+00};
static float p6_6[8]={-2.8213708496e+02,2.5735435486e+00,-7.8197475523e-03,1.1915668438e-05,-8.8851557223e-09,2.5953096009e-12,0.0000000000e+00,0.0000000000e+00};
static float a6_6[8]={2.9427927246e+03,-3.2528903961e+01,1.5111996233e-01,-3.8087458233e-04,5.6586748087e-07,-4.9696613491e-10,2.3930217617e-13,-4.8797535004e-17};

static float sogliap_1[9]={100,215,238,100,215,230,239,224,282};
static float sogliaa_1[9]={100,235,311,140,252,287,340,276,342};
static float sigmap_1[9]={2.00,2.59,3.10,2.32,2.64,3.17,3.18,6.67,2.77};
static float sigmaa_1[9]={2.00,2.92,2.11,1.99,2.16,3.08,2.18,4.23,3.17};
static float sogliap_2[9]={100,123,258,100,120,110,236,254,120};
static float sogliaa_2[9]={100,124,378,100,149,160,300,300,160};
static float sigmap_2[9]={2.00,3.71,2.96,2.00,2.39,2.49,2.22,2.43,2.48};
static float sigmaa_2[9]={2.00,3.73,2.52,2.00,2.06,1.85,1.73,1.96,2.46};
static float sogliap_3[9]={100,211,229,205,250,100,100,228,91};
static float sogliaa_3[9]={100,265,319,270,338,100,115,259,114};
static float sigmap_3[9]={2.00,3.02,3.00,2.74,2.95,2.00,1.98,3.41,2.32};
static float sigmaa_3[9]={2.00,2.19,2.56,3.04,2.54,2.00,2.15,3.28,2.59};
static float sogliap_4[9]={100,243,110,100,230,270,300,88,274};
static float sogliaa_4[9]={100,300,147,130,285,346,395,130,300};
static float sigmap_4[9]={2.00,1.61,2.68,1.88,2.67,2.68,4.81,2.35,5.79};
static float sigmaa_4[9]={2.00,3.93,2.95,1.79,2.26,2.30,1.92,2.21,3.46};
static float sogliap_5[9]={205,223,190,100,100,100,100,100,100};
static float sogliaa_5[9]={220,257,210,100,100,100,100,100,100};
static float sigmap_5[9]={4.40,5.26,4.65,2.00,2.00,2.00,2.00,2.00,2.00};
static float sigmaa_5[9]={2.72,4.04,2.99,2.00,2.00,2.00,2.00,2.00,2.00};
static float sogliap_6[9]={328,220,100,100,100,317,100,100,100};
static float sogliaa_6[9]={455,246,100,100,100,428,100,100,100};
static float sigmap_6[9]={6.55,4.44,2.00,2.00,2.00,6.85,2.00,2.00,2.00};
static float sigmaa_6[9]={10.76,3.44,2.00,2.00,2.00,7.49,2.00,2.00,2.00};


static float xminp_1[9]={0,200,216,100,200,209,224,200,250};
static float xmina_1[9]={0,220,274,130,240,260,300,260,330};
static float xmaxp_1[9]={1000,490,701,280,480,536,800,530,705};
static float xmaxa_1[9]={1000,1100,1400,790,1170,1000,1400,1080,1490};
static float xminp_2[9]={0,112,242,123,112,100,213,240,115};
static float xmina_2[9]={0,144,333,182,148,144,250,300,152};
static float xmaxp_2[9]={1000,390,597,673,289,350,586,590,305};
static float xmaxa_2[9]={1000,800,1180,1276,718,941,1180,1176,823};
static float xminp_3[9]={0,190,220,186,220,0,82,207,78};
static float xmina_3[9]={0,250,290,244,304,0,102,230,102};
static float xmaxp_3[9]={1000,600,580,500,727,1000,330,460,280};
static float xmaxa_3[9]={1000,1160,1400,995,1400,1000,800,990,800};
static float xminp_4[9]={0,210,100,95,200,242,270,75,225};
static float xmina_4[9]={0,265,140,130,244,300,340,100,260};
static float xmaxp_4[9]={1000,600,270,250,600,600,800,300,727};
static float xmaxa_4[9]={1000,1200,752,800,1210,1200,1500,790,1400};
static float xminp_5[9]={195,205,180,0,0,0,0,0,0};
static float xmina_5[9]={225,245,200,0,0,0,0,0,0};
static float xmaxp_5[9]={400,500,400,1000,1000,1000,1000,1000,1000};
static float xmaxa_5[9]={920,1160,860,1000,1000,1000,1000,1000,1000};
static float xminp_6[9]={260,200,0,0,0,267,0,0,0};
static float xmina_6[9]={400,230,0,0,0,350,0,0,0};
static float xmaxp_6[9]={800,530,1000,1000,1000,1000,1000,1000,1000};
static float xmaxa_6[9]={1100,1100,1000,1000,1000,950,1000,1000,1000};




void Classe_evento::Correggi_Tempo_Phos(int ip,int ih,int j)
{
  if(evento.z[j]>2){return;}
  int iflag=0;
  
  if(ip==1 && ih>=2 && ih<=9){iflag=1;}
  if(ip==2 && ih>=2 && ih<=3){iflag=1;}
  if(ip==2 && ih>=5 && ih<=9){iflag=1;}
  if(ip==3 && ih>=2 && ih<=5){iflag=1;}
  if(ip==3 &&ih>6){iflag=1;}
  if(ip==4 && ih>=2 && ih<=9){iflag=1;}
  if(ip==5 &&ih<4){iflag=1;}
  if(ip==6 &&ih<3){iflag=1;}
  if(ip==6 &&ih==6){iflag=1;}
  //il 6.3 sarebbe buono ma non si puo' correggere perche' mancano i dati di Fiasco

  // Quando non si puo' correggere il tempo, si mettono la v e la e =-1

  if(ip==2 && ih==4 && (evento.z[j]==1||evento.z[j]==2))
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    } 
  if(ip==3 && ih==6 && (evento.z[j]==1||evento.z[j]==2))
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    } 
  if(ip==5 && ih>3 && (evento.z[j]==1||evento.z[j]==2))
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    } 
  if(ip==6 && ih>2 && ih<6 && (evento.z[j]==1||evento.z[j]==2))
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    } 
  if(ip==6 && ih>6 && (evento.z[j]==1||evento.z[j]==2))
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    } 


  if(iflag==0){return;}

  // Correzioni non riuscite
  if(ip==1 && ih==2 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==1 && ih==6 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==1 && ih==8 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==1 && ih==9 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==2 && ih==3 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==2 && ih==8 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==3 && ih==2 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==3 && ih==7 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==3 && ih==8 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==4 && ih==3 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==4 && ih==5 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==4 && ih==7 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==4 && ih==9 &&evento.z[j]==1)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }
  if(ip==5||ip==6)
    {
      evento.vpartlab[j]=-1;
      evento.epartlab[j]=-1;
      evento.vpartcm[j]=-1;
      evento.epartcm[j]=-1;
      return;
    }





  float newv=-1;
  float newga=-1;
  float newt=-1;
  float newe=-1;


	      //int indice=ip*10+ih;
	      //    TH2F *hnewgavp=(TH2F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",3000+indice));
	      // TH2F *hnewgava=(TH2F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",3100+indice));
	      // TH2F *hnewganewv=(TH2F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",3200+indice));
	      // TH1F *hnewep=(TH1F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",3300+indice));
	      // TH1F *hnewea=(TH1F*)gROOT->GetListOfSpecials()->FindObject(Form("h%d",3400+indice));

	      if(ip==1){newga=ga0_1[ih-1]+ga1_1[ih-1]*evento.phosga[j];}
	      if(ip==2){newga=ga0_2[ih-1]+ga1_2[ih-1]*evento.phosga[j];}
	      if(ip==3){newga=ga0_3[ih-1]+ga1_3[ih-1]*evento.phosga[j];}
	      if(ip==4){newga=ga0_4[ih-1]+ga1_4[ih-1]*evento.phosga[j];}
	      if(ip==5){newga=ga0_5[ih-1]+ga1_5[ih-1]*evento.phosga[j];}
	      if(ip==6){newga=ga0_6[ih-1]+ga1_6[ih-1]*evento.phosga[j];}
	      float sogliap=0;
	      float sogliaa=0;
	      if(ip==1){sogliap=sogliap_1[ih-1];sogliaa=sogliaa_1[ih-1];}
	      if(ip==2){sogliap=sogliap_2[ih-1];sogliaa=sogliaa_2[ih-1];}
	      if(ip==3){sogliap=sogliap_3[ih-1];sogliaa=sogliaa_3[ih-1];}
	      if(ip==4){sogliap=sogliap_4[ih-1];sogliaa=sogliaa_4[ih-1];}
	      if(ip==5){sogliap=sogliap_5[ih-1];sogliaa=sogliaa_5[ih-1];}
	      if(ip==6){sogliap=sogliap_6[ih-1];sogliaa=sogliaa_6[ih-1];}
	      float sigmap=0,sigmaa=0;
	      if(ip==1){sigmap=sigmap_1[ih-1];sigmaa=sigmaa_1[ih-1];}
	      if(ip==2){sigmap=sigmap_2[ih-1];sigmaa=sigmaa_2[ih-1];}
	      if(ip==3){sigmap=sigmap_3[ih-1];sigmaa=sigmaa_3[ih-1];}
	      if(ip==4){sigmap=sigmap_4[ih-1];sigmaa=sigmaa_4[ih-1];}
	      if(ip==5){sigmap=sigmap_5[ih-1];sigmaa=sigmaa_5[ih-1];}
	      if(ip==6){sigmap=sigmap_6[ih-1];sigmaa=sigmaa_6[ih-1];}
	      
	      float xminp=0,xmaxp=0,xmina=0,xmaxa=0;
	      if(ip==1){xminp=xminp_1[ih-1];xmaxp=xmaxp_1[ih-1];}
	      if(ip==1){xmina=xmina_1[ih-1];xmaxa=xmaxa_1[ih-1];}
	      if(ip==2){xminp=xminp_2[ih-1];xmaxp=xmaxp_2[ih-1];}
	      if(ip==2){xmina=xmina_2[ih-1];xmaxa=xmaxa_2[ih-1];}
	      if(ip==3){xminp=xminp_3[ih-1];xmaxp=xmaxp_3[ih-1];}
	      if(ip==3){xmina=xmina_3[ih-1];xmaxa=xmaxa_3[ih-1];}
	      if(ip==4){xminp=xminp_4[ih-1];xmaxp=xmaxp_4[ih-1];}
	      if(ip==4){xmina=xmina_4[ih-1];xmaxa=xmaxa_4[ih-1];}
	      if(ip==5){xminp=xminp_5[ih-1];xmaxp=xmaxp_5[ih-1];}
	      if(ip==5){xmina=xmina_5[ih-1];xmaxa=xmaxa_5[ih-1];}
	      if(ip==6){xminp=xminp_6[ih-1];xmaxp=xmaxp_6[ih-1];}
	      if(ip==6){xmina=xmina_6[ih-1];xmaxa=xmaxa_6[ih-1];}




	      float cp[8],ca[8];
	      if(ip==1 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_2[i];
		      ca[i]=a1_2[i];
		    }
		}
	      if(ip==1 && ih==3)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_3[i];
		      ca[i]=a1_3[i];
		    }
		}
	      if(ip==1 && ih==4)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_4[i];
		      ca[i]=a1_4[i];
		    }
		}
	      if(ip==1 && ih==5)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_5[i];
		      ca[i]=a1_5[i];
		    }
		}
	      if(ip==1 && ih==6)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_6[i];
		      ca[i]=a1_6[i];
		    }
		}
	      if(ip==1 && ih==7)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_7[i];
		      ca[i]=a1_7[i];
		    }
		}
	      if(ip==1 && ih==8)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_8[i];
		      ca[i]=a1_8[i];
		    }
		}
	      if(ip==1 && ih==9)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p1_9[i];
		      ca[i]=a1_9[i];
		    }
		}


	      if(ip==2 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_2[i];
		      ca[i]=a2_2[i];
		    }
		}
	      if(ip==2 && ih==3)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_3[i];
		      ca[i]=a2_3[i];
		    }
		}
	      if(ip==2 && ih==4)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_4[i];
		      ca[i]=a2_4[i];
		    }
		}
	      if(ip==2 && ih==5)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_5[i];
		      ca[i]=a2_5[i];
		    }
		}
	      if(ip==2 && ih==6)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_6[i];
		      ca[i]=a2_6[i];
		    }
		}
	      if(ip==2 && ih==7)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_7[i];
		      ca[i]=a2_7[i];
		    }
		}
	      if(ip==2 && ih==8)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_8[i];
		      ca[i]=a2_8[i];
		    }
		}
	      if(ip==2 && ih==9)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p2_9[i];
		      ca[i]=a2_9[i];
		    }
		}

	      if(ip==3 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_2[i];
		      ca[i]=a3_2[i];
		    }
		}
	      if(ip==3 && ih==3)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_3[i];
		      ca[i]=a3_3[i];
		    }
		}
	      if(ip==3 && ih==4)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_4[i];
		      ca[i]=a3_4[i];
		    }
		}
	      if(ip==3 && ih==5)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_5[i];
		      ca[i]=a3_5[i];
		    }
		}

	      if(ip==3 && ih==7)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_7[i];
		      ca[i]=a3_7[i];
		    }
		}
	      if(ip==3 && ih==8)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_8[i];
		      ca[i]=a3_8[i];
		    }
		}
	      if(ip==3 && ih==9)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p3_9[i];
		      ca[i]=a3_9[i];
		    }
		}
	      if(ip==4 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_2[i];
		      ca[i]=a4_2[i];
		    }
		}
	      if(ip==4 && ih==3)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_3[i];
		      ca[i]=a4_3[i];
		    }
		}
	      if(ip==4 && ih==4)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_4[i];
		      ca[i]=a4_4[i];
		    }
		}
	      if(ip==4 && ih==5)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_5[i];
		      ca[i]=a4_5[i];
		    }
		}
	      if(ip==4 && ih==6)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_6[i];
		      ca[i]=a4_6[i];
		    }
		}
	      if(ip==4 && ih==7)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_7[i];
		      ca[i]=a4_7[i];
		    }
		}
	      if(ip==4 && ih==8)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_8[i];
		      ca[i]=a4_8[i];
		    }
		}
	      if(ip==4 && ih==9)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p4_9[i];
		      ca[i]=a4_9[i];
		    }
		}
	      if(ip==5 && ih==1)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p5_1[i];
		      ca[i]=a5_1[i];
		    }
		}
	      if(ip==5 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p5_2[i];
		      ca[i]=a5_2[i];
		    }
		}
	      if(ip==5 && ih==3)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p5_3[i];
		      ca[i]=a5_3[i];
		    }
		}
	      if(ip==6 && ih==1)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p6_1[i];
		      ca[i]=a6_1[i];
		    }
		}
	      if(ip==6 && ih==2)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p6_2[i];
		      ca[i]=a6_2[i];
		    }
		}
	      if(ip==6 && ih==6)
		{
		  for(int i=0;i<8;i++)
		    {
		      cp[i]=p6_6[i];
		      ca[i]=a6_6[i];
		    }
		}


	      int ispecial=0;
	      if(ip==4 && ih==4){ispecial=1;}
	      if(ip==2 && ih==9){ispecial=1;}
	      if(ispecial==0)
		{
	      
if(evento.z[j]==1 && evento.a[j]==1)
  {
    //hnewgavp->Fill(evento.vpartlab[j],newga,1);    
 if(newga<sogliap)
   {newv=evento.vpartlab[j];}
 else
   {
     newv=cp[0]+cp[1]*newga+cp[2]*pow(newga,2)+cp[3]*pow(newga,3)+cp[4]*pow(newga,4)+cp[5]*pow(newga,5)+cp[6]*pow(newga,6)+cp[7]*pow(newga,7);

     if(newga>xmaxp)
       {
	 float x1=xmaxp-50;
	 float x2=xmaxp;
	 float y1=cp[0]+cp[1]*x1+cp[2]*pow(x1,2)+cp[3]*pow(x1,3)+cp[4]*pow(x1,4)+cp[5]*pow(x1,5)+cp[6]*pow(x1,6)+cp[7]*pow(x1,7);
float y2=cp[0]+cp[1]*x2+cp[2]*pow(x2,2)+cp[3]*pow(x2,3)+cp[4]*pow(x2,4)+cp[5]*pow(x2,5)+cp[6]*pow(x2,6)+cp[7]*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
       }
     if(newga<xminp)
       {
	 float x1=xminp+10;
	 float x2=xminp;
	 float y1=cp[0]+cp[1]*x1+cp[2]*pow(x1,2)+cp[3]*pow(x1,3)+cp[4]*pow(x1,4)+cp[5]*pow(x1,5)+cp[6]*pow(x1,6)+cp[7]*pow(x1,7);
float y2=cp[0]+cp[1]*x2+cp[2]*pow(x2,2)+cp[3]*pow(x2,3)+cp[4]*pow(x2,4)+cp[5]*pow(x2,5)+cp[6]*pow(x2,6)+cp[7]*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);	 
       }

     newv=gRandom->Gaus(newv,sigmap);
   }

  }
 if(evento.z[j]==2)
   {
     //hnewgava->Fill(evento.vpartlab[j],newga,1); 
 if(newga<sogliaa)
   {
     newv=evento.vpartlab[j];
   } 
 else
   {
     newv=ca[0]+ca[1]*newga+ca[2]*pow(newga,2)+ca[3]*pow(newga,3)+ca[4]*pow(newga,4)+ca[5]*pow(newga,5)+ca[6]*pow(newga,6)+ca[7]*pow(newga,7);

     if(newga>xmaxa)
       {
	 float x1=xmaxa-50;
	 float x2=xmaxa;
	 float y1=ca[0]+ca[1]*x1+ca[2]*pow(x1,2)+ca[3]*pow(x1,3)+ca[4]*pow(x1,4)+ca[5]*pow(x1,5)+ca[6]*pow(x1,6)+ca[7]*pow(x1,7);
	 float y2=ca[0]+ca[1]*x2+ca[2]*pow(x2,2)+ca[3]*pow(x2,3)+ca[4]*pow(x2,4)+ca[5]*pow(x2,5)+ca[6]*pow(x2,6)+ca[7]*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
       }
     if(newga<xmina)
       {
	 float x1=xmina+10;
	 float x2=xmina;
	 float y1=ca[0]+ca[1]*x1+ca[2]*pow(x1,2)+ca[3]*pow(x1,3)+ca[4]*pow(x1,4)+ca[5]*pow(x1,5)+ca[6]*pow(x1,6)+ca[7]*pow(x1,7);
	 float y2=ca[0]+ca[1]*x2+ca[2]*pow(x2,2)+ca[3]*pow(x2,3)+ca[4]*pow(x2,4)+ca[5]*pow(x2,5)+ca[6]*pow(x2,6)+ca[7]*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
       }


     newv=gRandom->Gaus(newv,sigmaa);
   }
   }
		}//ispecial=0



   if(ip==4 && ih==4)
     {
       newga=57.37+0.6236*(evento.phosga[j]-500);
     
       if(evento.z[j]==1 && evento.a[j]==1)
 	{
	  //hnewgavp->Fill(evento.vpartlab[j],newga,1);    

 		  if(newga<100)
 		    {newv=evento.vpartlab[j];}
 		  else
 		    {
 		      newv=-2310.59+86.825*newga-1.33823*pow(newga,2)+0.0111604*pow(newga,3)-5.40545e-05*pow(newga,4)+1.52066e-07*pow(newga,5)-2.30209e-10*pow(newga,6)+1.44677e-13*pow(newga,7);
 if(newga>xmaxp)
   {
 	 float x1=xmaxp-50;
	 float x2=xmaxp;
	 float y1=-2310.59+86.825*x1-1.33823*pow(x1,2)+0.0111604*pow(x1,3)-5.40545e-05*pow(x1,4)+1.52066e-07*pow(x1,5)-2.30209e-10*pow(x1,6)+1.44677e-13*pow(x1,7);
float y2=-2310.59+86.825*x2-1.33823*pow(x2,2)+0.0111604*pow(x2,3)-5.40545e-05*pow(x2,4)+1.52066e-07*pow(x2,5)-2.30209e-10*pow(x2,6)+1.44677e-13*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
   }
 newv=gRandom->Gaus(newv,2.);
 		    }

 	}//z=1
       if(evento.z[j]==2)
 	{
	  //hnewgava->Fill(evento.vpartlab[j],newga,1); 

 		  if(newga<150)
 		    {newv=evento.vpartlab[j];}
 		  else
 		    {
 		      newv=-127.266+2.2369*newga-0.0100762*newga*newga+2.19455e-05*pow(newga,3)-2.26366e-08*pow(newga,4)+8.93669e-12*pow(newga,5);
		      if(newga>xmaxa)
   {
 	 float x1=xmaxa-50;
	 float x2=xmaxa;
	 float y1=-127.266+2.2369*x1-0.0100762*x1*x1+2.19455e-05*pow(x1,3)-2.26366e-08*pow(x1,4)+8.93669e-12*pow(x1,5);
	 float y2=-127.266+2.2369*x2-0.0100762*x2*x2+2.19455e-05*pow(x2,3)-2.26366e-08*pow(x2,4)+8.93669e-12*pow(x2,5);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
   } 

 		      		      newv=gRandom->Gaus(newv,1.3);

 		    }	  

 	}//z==2
     }//4.4

   if(ip==2 && ih==9)
     {

 newga=80.09+1.576*(evento.phosga[j]-500);
	      
 if(evento.z[j]==1 &&evento.a[j]==1) 
   {

     //hnewgavp->Fill(evento.vpartlab[j],newga,1);
 		  if(newga<125)
 		    {newv=evento.vpartlab[j];}
 		  else
 		    {
 newv=-5646.18+184.16*newga-2.51077*pow(newga,2)+0.0186537*pow(newga,3)-8.13425e-05*pow(newga,4)+2.08309e-07*pow(newga,5)-2.90438e-10*pow(newga,6)+1.7031e-13*pow(newga,7);

 if(newga>xmaxp)
   {
 	 float x1=xmaxp-50;
	 float x2=xmaxp; 
   float y1=-5646.18+184.16*x1-2.51077*pow(x1,2)+0.0186537*pow(x1,3)-8.13425e-05*pow(x1,4)+2.08309e-07*pow(x1,5)-2.90438e-10*pow(x1,6)+1.7031e-13*pow(x1,7);
   float y2=-5646.18+184.16*x2-2.51077*pow(x2,2)+0.0186537*pow(x2,3)-8.13425e-05*pow(x2,4)+2.08309e-07*pow(x2,5)-2.90438e-10*pow(x2,6)+1.7031e-13*pow(x2,7);
newv=newga*(y1-y2)/(x1-x2)+(x1*y2-x2*y1)/(x1-x2);
   }


 newv=gRandom->Gaus(newv,2.4);
 		    }

   }//z=1
  if(evento.z[j]==2)
    {
      //hnewgava->Fill(evento.vpartlab[j],newga,1); 

 		  if(newga<175)
 		    {newv=evento.vpartlab[j];}
 		  else
 		    {
 		      newv=-137.18+1.614*newga-0.00192151*newga*newga-1.89748e-05*pow(newga,3)+ 8.37152e-08*pow(newga,4)-1.44155e-10*pow(newga,5)+1.16262e-13*pow(newga,6)-3.64328e-17*pow(newga,7);
		      if(newga>xmaxa)
			{
 	 //float x1=xmaxa-50;
	 //float x2=xmaxa;
	 //float y1=-137.18+1.614*x1-0.00192151*x1*x1-1.89748e-05*pow(x1,3)+ 8.37152e-08*pow(x1,4)-1.44155e-10*pow(x1,5)+1.16262e-13*pow(x1,6)-3.64328e-17*pow(x1,7);
	 //float y2=-137.18+1.614*x2-0.00192151*x2*x2-1.89748e-05*pow(x2,3)+ 8.37152e-08*pow(x2,4)-1.44155e-10*pow(x2,5)+1.16262e-13*pow(x2,6)-3.64328e-17*pow(x2,7);

			}

 		      		      newv=gRandom->Gaus(newv,1.96);

 		    }


    }//z=2

     }//2.9

      if(newv>0)
	{
newt=Classe_geo::Getgeo()->phos.dist[ip-1][ih-1]/newv;
 newe=-1;
 if(evento.a[j]<50)
   {
newe=Classe_formule::v2e(newv,evento.a[j]);
   }

 evento.tvolo[j]=newt;
 evento.vpartlab[j]=newv;
 evento.epartlab[j]=newe;
 float vl[3];
 Classe_formule::polcar(evento.vpartlab[j],evento.thetalab[j],evento.philab[j],vl);
 evento.vpartlab_x[j]=vl[0];
 evento.vpartlab_y[j]=vl[1];
 evento.vpartlab_z[j]=vl[2];
 float vc[3];
 Classe_formule::lab2cm_classica(vc,vl,Classe_analisi::Getanalisi()->reazione.vcm);
 float VC,THEC,PHIC;
 Classe_formule::carpol(&VC,&THEC,&PHIC,vc);
 evento.vpartcm[j]=VC;
 evento.thetacm[j]=THEC;
 evento.phicm[j]=PHIC;
 evento.vpartcm_x[j]=vc[0];
 evento.vpartcm_y[j]=vc[1];
 evento.vpartcm_z[j]=vc[2];
 if(evento.a[j]<50)
   {
     evento.epartcm[j]=Classe_formule::v2e(evento.vpartcm[j],evento.a[j]);
   }
 else
   {
     evento.epartcm[j]=-1;
   }

	}//newv>0

      if(newv>0)
	{
	  //	  hnewganewv->Fill(newv,newga,1);
	  // hnewganewv->Fill(evento.vpartlab[j],newga,1);
if(evento.z[j]==1 &&evento.a[j]==1)
  {
    //    hnewep->Fill(newe,1);
    // hnewep->Fill(evento.epartlab[j],1);
  }	  
if(evento.z[j]==2&&evento.a[j]==4)
  {
    //    hnewea->Fill(newe,1);
    // hnewea->Fill(evento.epartlab[j],1);
  }


	}//newv>0

      //      hnewgavp=0;
      // hnewgava=0;
      // hnewganewv=0;
      //hnewep=0;
      // hnewea=0;

  return;
}

void Classe_evento::Riporto()
{
  int bmin=10;
  int bunch[6][9];
  float t0;
for(int ip=0;ip<6;ip++)
	{
	  for(int ih=0;ih<9;ih++)
	    {
	      bunch[ip][ih]=0;
	      if((expevent.phos_z[ip][ih]>0.1||expevent.phos_a[ip][ih]>0.1) && (expevent.phos_traw[ip][ih]>0)&& (Classe_geo::Getgeo()->phos.codice[ip][ih]>0)) //<=0 phos cattivi
		{
		  t0=expevent.phos_traw[ip][ih]-10;

		  if(t0<0)
		    {
		      bunch[ip][ih]=-1;
		    }
		  else
		    {
		      bunch[ip][ih]=t0/200;
		    }
		  if(bunch[ip][ih]<bmin && bunch[ip][ih]>=0)
		    {
		      bmin=bunch[ip][ih];
		    }

		}
	    }
	}

 for(int ip=0;ip<6;ip++)
	{
	  for(int ih=0;ih<9;ih++)
	    {
	      expevent.phos_tof[ip][ih]=-1;
	      if((expevent.phos_z[ip][ih]>0.1||expevent.phos_a[ip][ih]>0.1) && (expevent.phos_traw[ip][ih]>0)&& (Classe_geo::Getgeo()->phos.codice[ip][ih]>0))//phos cattivi <=0
		{
		  if(bunch[ip][ih]==bmin)
		    {
		      expevent.phos_tof[ip][ih]= expevent.phos_traw[ip][ih]-bmin*200;
		    }
		  else
		    {
		      expevent.phos_tof[ip][ih]=-1; 
		    }
		}
	    }
	}
 if(bmin==10){bmin=0;}
for(int ip=0;ip<8;ip++)
	{

	  if(exphector.ebaf[ip]>0)
	    {
	      float TBAF=exphector.tbaf[ip];
	       exphector.tbaf[ip]=TBAF-200*bmin;
	    }
	}
	 


}

void Classe_evento::RiportoGarfRCo() {
	int isgarf=0,isrco=0;
	int s1=0,s2=0,c1=0,c2=0,b;
	double tmin[2]={400,600};
	double tmax[2]={600,800};
	double tadd[3]={100,0,100};
	
	if(Classe_geo::Getgeo()->ThereIsRCO) isrco=1;
	if(Classe_geo::Getgeo()->ThereIsGarf) isgarf=1;
	
	
	//Osservo in quale bunch stanno le particelle
	for(int isec=0;(isec<8)&&isrco;isec++) {
		for(int isi=0;isi<9;isi++) {
			for(int icsi=0;icsi<7;icsi++) {
				if(expevent.rco_z[isec][isi][icsi]<0) continue;
				int tc=(int)(expevent.rco_tcode[isec][isi][icsi]+0.5);
				if((tc<1)||(tc>3)) continue;
				if(tc==2) {
					if((expevent.rco_tvolo[isec][isi][icsi]>=400)&&(expevent.rco_tvolo[isec][isi][icsi]<=600)) s1++;
					if((expevent.rco_tvolo[isec][isi][icsi]>=600)&&(expevent.rco_tvolo[isec][isi][icsi]<=800)) s2++;
				}
				else {
					if((expevent.rco_tvolo[isec][isi][icsi]>=350)&&(expevent.rco_tvolo[isec][isi][icsi]<=650)) c1++;
					if((expevent.rco_tvolo[isec][isi][icsi]>=550)&&(expevent.rco_tvolo[isec][isi][icsi]<=850)) c2++;
				}
			}
		}
	}
	for(int isec=0;(isec<24)&&isgarf;isec++) {
		for(int icsi=0;icsi<8;icsi++) {
			if(expevent.garf_z[isec][icsi]<0) continue;
			if((expevent.garf_tvolo[isec][icsi]>=350)&&(expevent.garf_tvolo[isec][icsi]<=650)) c1++;
			if((expevent.garf_tvolo[isec][icsi]>=550)&&(expevent.garf_tvolo[isec][icsi]<=850)) c2++;
		}
	}
	if(s1+s2+c1+c2==0) return; //se non ho info valide sul tempo esco!
	// SELEZIONE DEL BUNCH
	b=((s1==s2)?((c2>c1)?1:0):((s2>s1)?1:0));
	evento.bunch=b;
	// SCARTO LE PARTICELLE FUORI TEMPO!
	for(int isec=0;(isec<8)&&isrco;isec++) {
		for(int isi=0;isi<9;isi++) {
			for(int icsi=0;icsi<7;icsi++) {
				if(expevent.rco_z[isec][isi][icsi]<0) continue;
				int tc=(int)(expevent.rco_tcode[isec][isi][icsi]+0.5);
				if((tc<1)||(tc>3)) {
					expevent.rco_z[isec][isi][icsi]=-1;
					continue;
				}
				if((expevent.rco_tvolo[isec][isi][icsi]<tmin[b]-tadd[tc-1])||(expevent.rco_tvolo[isec][isi][icsi]>tmax[b]+tadd[tc-1])) {
					expevent.rco_z[isec][isi][icsi]=-1;
				}
			}
		}
	}
	for(int isec=0;(isec<24)&&isgarf;isec++) {
		for(int icsi=0;icsi<8;icsi++) {
			if(expevent.garf_z[isec][icsi]<0) continue;
			if((expevent.garf_tvolo[isec][icsi]<tmin[b]-100.)||(expevent.garf_tvolo[isec][icsi]>tmax[b]+100.)) {
				expevent.garf_z[isec][icsi]=-1;
			}
		}
	}
	return;
}

//Nuova versione del mixatore dal 29-01-21
vector<vector <float> > Classe_evento::prel2part(vector<vector<float> > ppart1, vector<vector<float> > ppart2){

	//if(ppart1.size()!=ppart2.size()) {cout << "prel2part: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<vector<float> > prel_vect;
  if(ppart1!=ppart2){
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=0; j<ppart2.size(); j++){
	vector<float> prel_temp;
	float ared=(float)ppart1.at(i).at(5)*(float)ppart2.at(j).at(5)/((float)ppart1.at(i).at(5)+(float)ppart2.at(j).at(5));
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	for(int k=0;k<3;k++){
	  prel_temp.push_back(ared*(ppart1.at(i).at(k)/(float)a1-ppart2.at(j).at(k)/(float)a2));
	}
	prel_temp.push_back(ppart1.at(i).at(6));
	prel_temp.push_back(ppart2.at(j).at(6));

	prel_vect.push_back(prel_temp);
	prel_temp.clear();
      }
    }
  }
  else {
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=i+1; j<ppart2.size(); j++){
	vector<float> prel_temp;
	float ared=(float)ppart1.at(i).at(5)*(float)ppart2.at(j).at(5)/((float)ppart1.at(i).at(5)+(float)ppart2.at(j).at(5));
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	for(int k=0;k<3;k++){
	  prel_temp.push_back(ared*(ppart1.at(i).at(k)/(float)a1-ppart2.at(j).at(k)/(float)a2));
	}
	prel_temp.push_back(ppart1.at(i).at(6));
	prel_temp.push_back(ppart2.at(j).at(6));

	prel_vect.push_back(prel_temp);
	prel_temp.clear();
      }
    }
  }
  return prel_vect;
}


vector<vector<float> > Classe_evento::pmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart, vector<vector<float> > ppart2){
//if(ppart.size()!=ppart2.size()) {cout << "pmix: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}


  vector<vector<float> > pmix_vect;
  for(UInt_t i=0; i<ppart.size(); i++){
    for(UInt_t j=0; j<tableP.at(index).size(); j++){
      if(!tableP.at(index).at(j).empty()){
	for(UInt_t q=0; q<tableP.at(index).at(j).size(); q++){
	  if(!tableP.at(index).at(j).at(q).empty()){
	    vector<float> p_temp;
	    float ared=(float)ppart.at(i).at(5)*(float)tableP.at(index).at(j).at(q).at(5)/((float)ppart.at(i).at(5)+(float)tableP.at(index).at(j).at(q).at(5));
	    float a1=(float)ppart.at(i).at(5);
	    float a2=(float)tableP.at(index).at(j).at(q).at(5);
	    
	    if(ppart.at(i).at(3)!=tableP.at(index).at(j).at(q).at(3)){ //no part in stesso riv
	    	//if(ppart2.at(q).at(4)== tableP.at(index).at(j).at(q).at(4) && ppart2.at(q).at(5)== tableP.at(index).at(j).at(q).at(5) ){ //mixo ppart1 solo con ioni di stesso Z,A di ppart2 in Table || ppart2 puo avere size diverso di TableP 
	    		for(int k=0; k<3; k++){
		  	p_temp.push_back(ared*(ppart.at(i).at(k)/(float)a1-tableP.at(index).at(j).at(q).at(k)/(float)a2));
	    		}
	    		pmix_vect.push_back(p_temp);
	    		p_temp.clear();
	    		//} 
			//else continue;
		}
	    else continue;
	      }
	}
      }
    }
  }
  return pmix_vect;
}

vector<vector<float> > Classe_evento::pcmmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart, vector<vector<float> > ppart2){
//if(ppart.size()!=ppart2.size()) {cout << "pcmmix: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<vector<float> > pmix_vect;
  for(UInt_t i=0; i<ppart.size(); i++){
    for(UInt_t j=0; j<tableP.at(index).size(); j++){
      if(!tableP.at(index).at(j).empty()){
	for(UInt_t q=0; q<tableP.at(index).at(j).size(); q++){
	  if(!tableP.at(index).at(j).at(q).empty()){
	    vector<float> p_temp;
	    float ared=(float)ppart.at(i).at(5)*(float)tableP.at(index).at(j).at(q).at(5)/((float)ppart.at(i).at(5)+(float)tableP.at(index).at(j).at(q).at(5));
	    float a1=(float)ppart.at(i).at(5);
	    float a2=(float)tableP.at(index).at(j).at(q).at(5);
	    if(ppart.at(i).at(3)!=tableP.at(index).at(j).at(q).at(3)){ //no part in stesso riv
	    	//if(ppart2.at(q).at(4)== tableP.at(index).at(j).at(q).at(4) && ppart2.at(q).at(5)== tableP.at(index).at(j).at(q).at(5) ){ //mixo ppart1 solo con ioni di stesso Z,A di ppart2 in Table
	    		for(int k=0; k<3; k++){
		  	p_temp.push_back((a1*ppart.at(i).at(k)+a2*tableP.at(index).at(j).at(q).at(k))/(a1+a2));
	    		}
	    		pmix_vect.push_back(p_temp);
	    		p_temp.clear();
	    		//} 
			//else continue;
		}
	    else continue;
	      }
	}
      }
    }
  }
  return pmix_vect;
}



vector<vector<float> > Classe_evento::Vcmmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart, vector<vector<float> > ppart2){
//if(ppart.size()!=ppart2.size()) {cout << "pcmmix: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<vector<float> > vmix_vect;
  for(UInt_t i=0; i<ppart.size(); i++){
    for(UInt_t j=0; j<tableP.at(index).size(); j++){
      if(!tableP.at(index).at(j).empty()){
	for(UInt_t q=0; q<tableP.at(index).at(j).size(); q++){
	  if(!tableP.at(index).at(j).at(q).empty()){
	    vector<float> v_temp;
	    float ared=(float)ppart.at(i).at(5)*(float)tableP.at(index).at(j).at(q).at(5)/((float)ppart.at(i).at(5)+(float)tableP.at(index).at(j).at(q).at(5));
	    float a1=(float)ppart.at(i).at(5);
	    float a2=(float)tableP.at(index).at(j).at(q).at(5);
	    if(ppart.at(i).at(3)!=tableP.at(index).at(j).at(q).at(3)){ //no part in stesso riv
	    	//if(ppart2.at(q).at(4)== tableP.at(index).at(j).at(q).at(4) && ppart2.at(q).at(5)== tableP.at(index).at(j).at(q).at(5) ){ //mixo ppart1 solo con ioni di stesso Z,A di ppart2 in Table
	    		for(int k=0; k<3; k++){

			  //		  	p_temp.push_back((a1*ppart.at(i).at(k)+a2*tableP.at(index).at(j).at(q).at(k))/(a1+a2));
			  v_temp.push_back((a1*Classe_formule::cluce*ppart.at(i).at(k)/(a1*931.5)+a2*Classe_formule::cluce*tableP.at(index).at(j).at(q).at(k)/(a2*931.5))/(a1+a2));

	    		}
	    		vmix_vect.push_back(v_temp);
	    		v_temp.clear();
	    		//} 
			//else continue;
		}
	    else continue;
	      }
	}
      }
    }
  }
  return vmix_vect;
}



vector<float>  Classe_evento::thetarelmix(int index, vector<vector<vector<vector<float> > > > tableP, vector<vector<float> > ppart, vector<vector<float> > ppart2){
//if(ppart.size()!=ppart2.size()) {cout << "thetarelmix: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<float>  theta_vect;
  for(UInt_t i=0; i<ppart.size(); i++){
    for(UInt_t j=0; j<tableP.at(index).size(); j++){
      if(!tableP.at(index).at(j).empty()){
	for(UInt_t q=0; q<tableP.at(index).at(j).size(); q++){
	  if(!tableP.at(index).at(j).at(q).empty()){
	    float thetarel_temp=0;
	    float ared=(float)ppart.at(i).at(5)*(float)tableP.at(index).at(j).at(q).at(5)/((float)ppart.at(i).at(5)+(float)tableP.at(index).at(j).at(q).at(5));
	    float a1=(float)ppart.at(i).at(5);
	    float a2=(float)tableP.at(index).at(j).at(q).at(5);
	    if(ppart.at(i).at(3)!=tableP.at(index).at(j).at(q).at(3)){ //no part in stesso riv
	    	//if(ppart2.at(q).at(4)== tableP.at(index).at(j).at(q).at(4) && ppart2.at(q).at(5)== tableP.at(index).at(j).at(q).at(5) ){ //mixo ppart1 solo con ioni di stesso Z,A di ppart2 in Table
		thetarel_temp=Classe_formule::thetarel(ppart.at(i),tableP.at(index).at(j).at(q));
	    	theta_vect.push_back(thetarel_temp);
	    		//} 
			//else continue;
		}
	    else continue;
	      }
	}
      }
    }
  }
  return theta_vect;
}


vector<vector<float> > Classe_evento::pcm2part(vector<vector<float> >ppart1, vector<vector<float> >ppart2){
		
	//if(ppart1.size()!=ppart2.size()) {cout << "pcm2part: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<vector<float> > pcm_vect;
  if(ppart1!=ppart2){
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=0; j<ppart2.size(); j++){
	vector<float> pcm_temp;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	for(int k=0;k<3;k++){
	  pcm_temp.push_back((a1*ppart1.at(i).at(k)+a2*ppart2.at(j).at(k))/atot);
	}
	pcm_vect.push_back(pcm_temp);
	pcm_temp.clear();
      }
    }
  }
  else {
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=i+1; j<ppart2.size(); j++){
	vector<float> pcm_temp;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	for(int k=0;k<3;k++){
	  pcm_temp.push_back((a1*ppart1.at(i).at(k)+a2*ppart2.at(j).at(k))/atot);
	}	
	pcm_vect.push_back(pcm_temp);
	pcm_temp.clear();
      }
    }
  }
  return pcm_vect;
	
}

vector<vector<float> > Classe_evento::Vcm2part(vector<vector<float> >ppart1, vector<vector<float> >ppart2){
		
	//if(ppart1.size()!=ppart2.size()) {cout << "pcm2part: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<vector<float> > vcm_vect;
  if(ppart1!=ppart2){
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=0; j<ppart2.size(); j++){
	vector<float> vcm_temp;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	for(int k=0;k<3;k++){
	  vcm_temp.push_back((a1*Classe_formule::cluce*ppart1.at(i).at(k)/(a1*931.5)+a2*Classe_formule::cluce*ppart2.at(j).at(k)/(a2*931.5))/atot);
	}
	vcm_vect.push_back(vcm_temp);
	vcm_temp.clear();
      }
    }
  }
  else {
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=i+1; j<ppart2.size(); j++){
	vector<float> vcm_temp;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	for(int k=0;k<3;k++){

	  vcm_temp.push_back((a1*Classe_formule::cluce*ppart1.at(i).at(k)/(a1*931.5)+a2*Classe_formule::cluce*ppart2.at(j).at(k)/(a2*931.5))/atot);
	}	
	vcm_vect.push_back(vcm_temp);
	vcm_temp.clear();
      }
    }
  }
  return vcm_vect;
	
}

vector<float>  Classe_evento::thetarel2part(vector<vector<float> >ppart1, vector<vector<float> >ppart2){
		
	//if(ppart1.size()!=ppart2.size()) {cout << "thetarel2part: ppart1 e ppart2 non hanno le stesse dimensioni" << endl;}

  vector<float>  theta_vect;
  if(ppart1!=ppart2){
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=0; j<ppart2.size(); j++){
	float theta_temp=0;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	
	theta_temp=Classe_formule::thetarel(ppart1.at(i),ppart2.at(j));
	theta_vect.push_back(theta_temp);
	}
    }
  }
  else {
    for(UInt_t i=0; i<ppart1.size(); i++){
      for(UInt_t j=i+1; j<ppart2.size(); j++){
	float theta_temp;
	float a1=(float)ppart1.at(i).at(5);
	float a2=(float)ppart2.at(j).at(5);
	float atot=a1+a2;
	
	theta_temp=Classe_formule::thetarel(ppart1.at(i),ppart2.at(j));
	theta_vect.push_back(theta_temp);
      }
    }
  }
  return theta_vect;
	
}

//index è il tipo di correlazione: per es alpha-alpha è 0, d-alpha è 1 etc. Il numero lo definisce l'utente
//codice: è il valore di evento.isfoxmix[j] definito dall'utente. Per esempio, se si definisce che le particelle di garfield hanno isformix=1 e quelle del RCO hanno 0, si può fare il mixatore solo su quelle con  1 o con 0
//se si fa il mixing fra 2 particelle uguali si deve mettere n1=2, z1, a1, -1, -1, -1
//se sono 2 diverse si deve fare n1=1 z1 a1, n2=1, z2, a2 

void Classe_evento::mixatore(int index, int codice, int n1, int z1, int a1, int n2, int z2, int a2)
{
  
  mixing.reset(index);
  
  
  //  cout<<"sono in mixatore"<<endl;
  int buffer=100;
  if(n1>1&&n2==-1){
   
         
      if(Classe_analisi::Getanalisi()->iokmix.at(index)>0&&Classe_analisi::Getanalisi()->iokmix.at(index)%buffer==0) Classe_analisi::Getanalisi()->ioktimesmix.at(index)++;
      Classe_analisi::Getanalisi()->jokmix.at(index)=Classe_analisi::Getanalisi()->iokmix.at(index)-Classe_analisi::Getanalisi()->ioktimesmix.at(index)*buffer;
      
      
      vector< vector<float> > Ppart;
     
     
      for(UInt_t m=0; m<evento.moltepl; m++){
	//	cout<<TMath::Nint(evento.a[m])<<" "<<a1<<" "<<TMath::Nint(evento.z[m])<<" "<<z1<<" "<<evento.isformix[m]<<" "<<codice<<endl;

	if(TMath::Nint(evento.a[m])==a1&&TMath::Nint(evento.z[m])==z1&&evento.isformix[m]==codice) {
	
	  vector<float> ppart;
	  vector<float> vpart;
	  
	  ppart.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_x[m]/Classe_formule::cluce);
	  ppart.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_y[m]/Classe_formule::cluce);
	  ppart.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_z[m]/Classe_formule::cluce);
	 
	  //metto come in quarta posizione il coderiv
	  ppart.push_back(evento.coderiv[m]);
	  ppart.push_back(evento.z[m]);
	  ppart.push_back(evento.a[m]);
	  ppart.push_back(m);
 
	  Ppart.push_back(ppart); 
	  ppart.clear();
	 
	  
	}
      }
       if(Ppart.size()<=0) return;
      //calcolo Prel
      
      mixing.Prel.at(index) = prel2part(Ppart,Ppart);
      mixing.Pcm.at(index) = pcm2part(Ppart,Ppart);
      mixing.Vcm.at(index) = Vcm2part(Ppart,Ppart);

      mixing.Thetarel.at(index)=thetarel2part(Ppart,Ppart);
     
      //Mi calcolo il mix con gli elementi messi nella tabella (primo giro no)
      if(Classe_analisi::Getanalisi()->iokmix.at(index)>0) {
	mixing.Pmix.at(index)=pmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart,Ppart);
	mixing.Pcm_mix.at(index)=pcmmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart,Ppart);
	mixing.Vcm_mix.at(index)=Vcmmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart,Ppart);

	mixing.Thetamix.at(index)=thetarelmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart,Ppart);
      }
      
      //aggiungo alla tabella con buffer circolare il vettore precedente
      //	 Classe_analisi::Getanalisi()->tableP1.at(index).resize(Classe_analisi::Getanalisi()->tableP1.at(index).size()+1);
      // cout<<Classe_analisi::Getanalisi()->tableP1.at(index).size()<<" "<<Classe_analisi::Getanalisi()->tableP1.at(index).at(Classe_analisi::Getanalisi()->jokmix.at(index)).size()<<" "<<index<<" "<<Classe_analisi::Getanalisi()->jokmix.at(index)<<endl;
      Classe_analisi::Getanalisi()->tableP1.at(index).at(Classe_analisi::Getanalisi()->jokmix.at(index))=Ppart;
      //svuoto tutti i vector (esclusa la tabella) e incremento per ibe per il successivo evento scelto
      
      
      Ppart.clear();
      Classe_analisi::Getanalisi()->iokmix.at(index)++;
    
    
  }//due part uguali
  
  if(n1>0&&n2>0){
  
    
      if(Classe_analisi::Getanalisi()->iokmix.at(index)>0&&Classe_analisi::Getanalisi()->iokmix.at(index)%buffer==0) Classe_analisi::Getanalisi()->ioktimesmix.at(index)++;
      Classe_analisi::Getanalisi()->jokmix.at(index)=Classe_analisi::Getanalisi()->iokmix.at(index)-Classe_analisi::Getanalisi()->ioktimesmix.at(index)*buffer;
      
      vector<vector<float> > Ppart1;
      vector<vector<float> > Ppart2; 
      for(UInt_t m=0; m<evento.moltepl; m++){
	if(TMath::Nint(evento.a[m])==a1&&TMath::Nint(evento.z[m])==z1&&evento.isformix[m]==codice) {
	  vector<float> ppart1;
	  ppart1.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_x[m]/Classe_formule::cluce);
	  ppart1.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_y[m]/Classe_formule::cluce);
	  ppart1.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_z[m]/Classe_formule::cluce);
	  ppart1.push_back(evento.coderiv[m]);
	  ppart1.push_back(evento.z[m]);
	  ppart1.push_back(evento.a[m]);
	  ppart1.push_back(m);
	  Ppart1.push_back(ppart1);
	  ppart1.clear();
	 
	}
	
	
	if(TMath::Nint(evento.a[m])==a2&&TMath::Nint(evento.z[m])==z2&&evento.isformix[m]==codice) {
	  vector<float> ppart2;
	  ppart2.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_x[m]/Classe_formule::cluce);
	  ppart2.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_y[m]/Classe_formule::cluce);
	  ppart2.push_back(evento.a[m]*Classe_formule::amu*evento.vpartlab_z[m]/Classe_formule::cluce);
	  ppart2.push_back(evento.coderiv[m]);
	  ppart2.push_back(evento.z[m]);
	  ppart2.push_back(evento.a[m]);
	  ppart2.push_back(m);
	  Ppart2.push_back(ppart2);
	  ppart2.clear();
	  
	}
      }
      
      if(Ppart1.size()<=0 || Ppart2.size()<=0) return;
      
      /*cout << "MIX-Ppart1:\n";
	for(int i=0; i<Ppart1.size(); i++){
	for(int k=0; k<Ppart1.at(i).size(); k++){
	cout << Ppart1.at(i).at(k) << " ";
	}
	}
	
	cout << "\nMIX-Ppart2:\n";
	for(int i=0; i<Ppart2.size(); i++){
	for(int k=0; k<.size(); k++){
	cout << Ppart2.at(i).at(k) << " ";
	}
	}*/
      
      //calcolo Prel
      mixing.Prel.at(index) = prel2part(Ppart1,Ppart2);	 
      mixing.Pcm.at(index) = pcm2part(Ppart1,Ppart2);
      mixing.Vcm.at(index) = Vcm2part(Ppart1,Ppart2);

      mixing.Thetarel.at(index) = thetarel2part(Ppart1,Ppart2);
      
      
     /* cout << "\nMIX-Prel:\n";
	for(int i=0; i<mixing.Prel.at(index).size(); i++){
	for(int k=0; k<mixing.Prel.at(index).at(i).size(); k++){
	cout << mixing.Prel.at(index).at(i).at(k) << " ";
	}
	}
	cout << " \n";
	cout << a1*a2/((float)a1+(float)a2) << " MIX-Erel:\n";
	for(int i=0; i<mixing.Erel.at(index).size(); i++) cout << mixing.Erel.at(index).at(i) << " ";
	cout << " \n\n";*/
      
      //Mi calcolo il mix con gli elementi messi nella tabella (primo giro no)
      if(Classe_analisi::Getanalisi()->iokmix.at(index)>0) {
	
	
	/*cout << "Table 1:\n";
	  for(int i=0;i<Classe_analisi::Getanalisi()->tableP1.at(index).size(); i++){
	  //cout << i << endl;
	  if(!Classe_analisi::Getanalisi()->tableP1.at(index).at(i).empty()){
	  for(int j=0; j<Classe_analisi::Getanalisi()->tableP1.at(index).at(i).size(); j++){
	  //cout << j << endl;
	  if(!Classe_analisi::Getanalisi()->tableP1.at(index).at(i).at(j).empty()){
	  for(int q=0; q<Classe_analisi::Getanalisi()->tableP1.at(index).at(i).at(j).size(); q++){
	  //cout << q << endl;
	  cout << Classe_analisi::Getanalisi()->tableP1.at(index).at(i).at(j).at(q) << " " ;		
	  
	  }
	  }p
	  }
	  cout << " \n";	 
	  }
	  }
	  
	  cout << " \n";*/
	
	mixing.Pmix.at(index)=pmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart1,Ppart2);
	vector<vector<float> > pstep;
	mixing.Pcm_mix.at(index)=pcmmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart1,Ppart2);
	mixing.Vcm_mix.at(index)=Vcmmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart1,Ppart2);
	vector<vector<float> > pcmstep;
	mixing.Thetamix.at(index)=thetarelmix(index,Classe_analisi::Getanalisi()->tableP1,Ppart1,Ppart2);
	vector<float> thetastep; 

	
	/*   cout << "Table 2:\n";
	     for(int i=0;i<Classe_analisi::Getanalisi()->tableP2.at(index).size(); i++){
	     //cout << i << endl;
	     if(!Classe_analisi::Getanalisi()->tableP2.at(index).at(i).empty()){
	     for(int j=0; j<Classe_analisi::Getanalisi()->tableP2.at(index).at(i).size(); j++){
	     //cout << j << endl;
	     if(!Classe_analisi::Getanalisi()->tableP2.at(index).at(i).at(j).empty()){
	     for(int q=0; q<Classe_analisi::Getanalisi()->tableP2.at(index).at(i).at(j).size(); q++){
	     //cout << q << endl;
	     cout << Classe_analisi::Getanalisi()->tableP2.at(index).at(i).at(j).at(q) << " " ;		
	     
	     }
	     }
	     }
	     cout << " \n";	 
	     }
	     }
	     
	     cout << " \n";*/
	
	
	
	
	//pstep=pmix(index,Classe_analisi::Getanalisi()->tableP2,Ppart2,Ppart1);
	//mixing.Pmix.at(index).insert(mixing.Pmix.at(index).end(), pstep.begin(), pstep.end() );
	//pcmstep=pcmmix(index,Classe_analisi::Getanalisi()->tableP2,Ppart2,Ppart1);
	//mixing.Pcm_mix.at(index).insert(mixing.Pcm_mix.at(index).end(), pstep.begin(), pstep.end() );
	//thetastep=thetarelmix(index,Classe_analisi::Getanalisi()->tableP2,Ppart2,Ppart1);
	//mixing.Thetamix.at(index).insert(mixing.Thetamix.at(index).end(),thetastep.begin(), thetastep.end());
	
	/* cout << "Mix:\n";
	   for(int i=0; i<mixing.Pmix.at(index).size(); i++){
	   for(int k=0; k<mixing.Pmix.at(index).at(i).size(); k++){
	   cout << mixing.Pmix.at(index).at(i).at(k) << " ";
	   }
	   cout << " \n";
	   }
	   cout << " \n\n";*/
	
	
      }
      
      //	 Classe_analisi::Getanalisi()->tableP1.at(index).resize(Classe_analisi::Getanalisi()->tableP1.at(index).size()+1);
      Classe_analisi::Getanalisi()->tableP1.at(index).at(Classe_analisi::Getanalisi()->jokmix.at(index))=Ppart2;
      //	 Classe_analisi::Getanalisi()->tableP2.at(index).resize(Classe_analisi::Getanalisi()->tableP2.at(index).size()+1);
      Classe_analisi::Getanalisi()->tableP2.at(index).at(Classe_analisi::Getanalisi()->jokmix.at(index))=Ppart1;	 
      Ppart1.clear();
      Ppart2.clear();
     
      
      Classe_analisi::Getanalisi()->iokmix.at(index)++;
  }//fine alpha-deuton like
  
  return;
}
