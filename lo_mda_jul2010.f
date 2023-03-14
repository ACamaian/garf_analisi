

**********la chiamata alla function nel main program e':

c         elux_csi=lo_mda_tg(iz,ia,eres)
c modificata in subroutine lo_mda_tg(iz,ia,eres,aluce)
* dove iz,ia sono ovviamente z,a dello ione di cui si stanno calcolando delta_e
* ed e_residua (eres in cesio in MeV)

************* segue function usata per i tabelloni di
************* commissioning,isomix,dipolino,jacobi
************* luglio 2010: richiesta di Faby e Sandro di proseguire 
************* con una retta le uscite di luce vs. energia
************* ricordo che l'introduzione della retta e' arbitraria 
************* non si basa su nessuna info exp
************* se e quando ci saranno dati ulteriori si cerchera' di fare meglio
************* PER ORA CON RETTA ARBITRARIA LE USCITE DI LUCE VS.ENERGIA 
************* PER Z=1,2 SI INCROCIANO (Sandro luglio 2010).
c-------------------------------------------
c segue mda dec09-jan10+modifica per tangente 5/7/2010 Legnaro (Faby-Sandro)
c-------------------------------------------

************* UPGRADING 2011: nuovi parametri per la funzione per evitare 
*************                 il crossing tra p,d e t



      subroutine lo_mda_tg(iz,ia,eres,aluce)
      implicit none
      real eamin,ei,rloi,tga
      real aluce
      integer iz,ia
      real*4 par(8)
**************************************    ok 5/7/2010 per tangente
c---------------- retta tangente da e=42 MeV, lo= circa 1250
c segue 19/5/2010 param. usati per commissioning,isomix,dipolino,jacobi
c      data par/19.0147656,243.959916,0.422541642,-0.858159425,1.5065705,
c     &0.288082746,0.275077286,0.0121368638/

c *********** Nuovi parametri 2011 **************************************
      data par/12.2424988,288.645213,0.420237154,-0.805045025,1.58972328
     &,0.51370864,0.165660342,0.0181196633/



*************************************************
      real*4 xz,amas,eres,zc
      data eamin/10.5/
*
      xz  = iz 
      amas= ia
c------------------------------------------
      ei=eamin*amas
c-------------------------------------------
c eamin=ene/A da cui far partire retta tangente; 
c lux da cui far partire retta tangente,  viene calcolata da ei in poi
*
      zc= (amas*xz**2)**par(8)
      if(eres.le.ei)then
      aluce   = (par(1)+par(2)*exp(-par(3)*xz))*(1.+par(4)*zc)*
     &    eres**(par(5)-par(6)*exp(-par(7)*xz))
      else
        rloi= (par(1)+par(2)*exp(-par(3)*xz))*
     &      (1.+par(4)*zc)*
     &      ei**(par(5)-par(6)*exp(-par(7)*xz))
        tga= (par(1)+par(2)*exp(-par(3)*xz))*(1.+par(4)*zc)*
     &          (par(5)-par(6)*exp(-par(7)*xz))*
     &         ei**(par(5)-par(6)*exp(-par(7)*xz)-1)
        aluce = eres*tga+(rloi-tga*ei)
c come dire:        deriv = (ene-ei)*tga+rloi  cioe' ei=e_intercetto
      endif
      return
      end
c-------------------------------------------
