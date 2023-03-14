          subroutine relkin(apr,atar,aein,athe,ae1,athe2,ae2)
          implicit double precision (b-h,o-z)
c          write(6,*)aein
               radgr=57.29578d0
               euma=931.494d0
               bapr=apr
               batar=atar
               baein=aein
               bthe=athe
               CTH=dcos(bthe/radgr)
c               TGTH=tan(bthe/radgr)
               ecin=baein
               etopr=ecin+bapr
               pcpr=dsqrt(etopr*etopr-bapr*bapr)
               etot=ecin+bapr+batar
               bfac=etot*etot-pcpr*pcpr*cth*cth
               cfac=2.d0*(batar*ecin+bapr*bapr+bapr*batar)
               crad1=cfac*cfac-4.d0*bapr*bapr*bfac
               crad=dsqrt(crad1)
               eprtot=((etot*cfac+pcpr*cth*crad)/(2.d0*bfac))
               E1B=eprtot-bapr
               pcprn=dsqrt(eprtot*eprtot-bapr*bapr)
               E1=E1B*euma
               ae1=e1
               etar=(etot-e1b-bapr)
               pctgn=dsqrt(etar*etar-batar*batar)
               ectar=(etar-batar)*euma
               the2=dasin((pcprn/pctgn)*dsin(bthe/radgr))
               athe2=radgr*the2
               ae2=ectar
c         write(6,*)'    theta,     en cin(Mev), en tot(amu)'
c         write(6,*)'pr',athe,e1,eprtot
c         write(6,*)'tg',athe2,ectar,etar
         return
         end
