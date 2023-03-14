#include <stdio.h>
#include <iostream>
#include <TCanvas.h>
#include <TH2.h>
void pac()
{
  FILE *apri=fopen("limitizid_FAZIAPRE4825.txt","r");
  char riga[1000];
  int ib,iq,it,iz;
  float e1,e2,ea1,ea2,ea10;
TH1F *hz=new TH1F("hz","hz",20,5,25);
TH1F *hemin=new TH1F("hemin","hemin",200,0,1000);

 int z0=-1;
 int st=0;
  while(fscanf(apri,"\n%[^\n]",riga)!=EOF)
    {
      if(riga[0]!='#')
	{
	  sscanf(riga,"%d %d %d %d %f %f %f %f",&ib,&iq,&it,&iz,&e1,&e2,&ea1,&ea2);
	  if(ea1>-1 && iz>5){z0=iz;st=1;ea10=ea1;}
	  if(ea1<0&&st==1&&iz>5)
	    {
	  hz->Fill((float)z0);

	  if(z0==20)
	    {
	      hemin->Fill(ea10);
	      cout<<ea10<<" "<<ib<<" "<<iq<<" "<<it<<" "<<iz<<endl;
	    }

	  st=0;
	    }
	  if(ea1<0 &&iz>5){st=0;}
	}
    }
  fclose(apri);
TCanvas *c=new TCanvas("c","c");
 c->Draw();
hz->Draw();
TCanvas *cemin=new TCanvas("cemin","cemin");
 cemin->Draw();
hemin->Draw();
}
