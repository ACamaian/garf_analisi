      subroutine thetaflow(tflow,npart)
c      common/pvec/pvettore(0:2,0:499)
      common/perflowpvec/pvettore(3,500)
      common /eigen/W,X

      real*8 W(3),X(3)
      real p(500,3),pmod(500)
c      dimension vlpcm(3),vlablcp(3)
      real*8 t(3,3),vec(3,3),eigenval(3)
      real fascio(3),vett3(3),xx(3),amatrix(3,3)
      tflow=-1

      fascio(1)=0
      fascio(2)=0
      fascio(3)=1
c      write(6,*)'entro',npart

      do j=1,npart

c         write(6,*)j,pvettore(1,j),pvettore(2,j),pvettore(3,j)
           

        pmod(j)=sqrt(pvettore(1,j)**2+pvettore(2,j)**2+pvettore(3,j)**2)

      enddo
      icall=0
      do i=1,3
         do j=1,3
            t(i,j)=0
            do k=1,npart
c               t(i,j)=t(i,j)+p(k,i)*p(k,j)/pmod(k)               
               t(i,j)=t(i,j)+pvettore(i,k)*pvettore(j,k)/pmod(k)               


            enddo
            amatrix(i,j)=t(i,j)
               if(t(i,j).ne.0)icall=1
            vec(i,j)=0
         enddo
         eigenval(i)=0
      enddo
      
      in=3
      istep=1000
      if(icall.eq.1)then
      call FINDEIGEN(amatrix,in,xx,eigen,istep)
c      write(6,*)'Piu grosso=',eigen,xx
      costflow=abs(scaprod(fascio,xx)/
     c                      (sqrt(xx(1)**2+xx(2)**2+xx(3)**2)))
      tflow=57.296*acos(costflow)
c      write(6,*)tflow,costflow
      endif

c*************Metodo alternativo (+generale) per calcolare tutti gli autovalori e autovettori; spesso non converge. Nella maggior parte dei casi da' risultato uguale a al metodo precedente che calcola solo l'autovettore relativo all'autovalore piu' grosso e dato che serve solo quell'autovettore li' per calcolare il thetaflow, e' inutile complicare i calcoli

c      it=0 !aggiunto perche' non funzionava
c      do k=1,3
c         W(k)=0
c         X(k)=0
c      enddo
c      if(icall.eq.1)call calcola_autovec(t,vec,eigenval,it)
c      if(it.ge.0)then
c      if(eigenval(1).ne.0.or.eigenval(2).ne.0.or.eigenval(3).ne.0)then
c         write(6,*)'E=',eigenval
c         write(6,*)'vec=',vec
c      emax=eigenval(1)
c         do k=1,3
c            vett3(k)=vec(1,k)
c         enddo
c      if(eigenval(2).gt.emax)then
c         emax=eigenval(2)
c         do k=1,3
c            vett3(k)=vec(2,k)
c         enddo
c      endif
c      if(eigenval(3).gt.emax)then
c         emax=eigenval(3)
c         do k=1,3
c            vett3(k)=vec(3,k)
c         enddo
c      endif
c         tflow=57.296*acos(scaprod(vett3,fascio))
c         write(6,*)tflow,cos(tflow/57.296)
c      endif
c      endif

      return
      end



      SUBROUTINE FINDEIGEN(Matrix, n, x, EigenVal, steps)
      real Matrix(3,3)
      real x(3)
      integer steps
c        IMPLICIT NONE
c        INTEGER, INTENT(IN) :: n, steps  !n = order of matrix, steps = number of iterations
c        REAL, INTENT(IN), DIMENSION(n,n) :: Matrix(n,n)  !Input Matrix
c        REAL, INTENT(INOUT), DIMENSION(n) :: x !Eigenvector
c        REAL, INTENT(INOUT) :: EigenVal !Eigenvalue
c        INTEGER :: i, j
        
c        x  = 1 !Initialize eigen vector to any value.
        do i=1,3
           x(i)=1
        enddo
c        n=3
        DO i = 1, steps
                CALL MULMATRIX(Matrix, x, n)       !Multiply input matrix by eigenvector
                CALL FINDLARGEST(x, n, EigenVal)   !Find eigenvalue
                if(eigenval.eq.0.)return
c                IF(EigenVal == 0)return                        
                DO j = 1, n                        !Find eigenvector
                        x(j) = x(j)/EigenVal
                ENDDO        
        ENDDO
        return
        end
cEND SUBROUTINE FINDEIGEN
 
      SUBROUTINE MULMATRIX(a, b, n)
c        IMPLICIT NONE
c        INTEGER, INTENT(IN) :: n !matrix size
c        REAL, INTENT(IN), DIMENSION(n,n) :: a  !Matrix of order > 1
c        REAL, INTENT(INOUT), DIMENSION(n) :: b !1x1 matrix
      real a(3,3)
      real b(3),temp(3)
c        INTEGER i, j
c        REAL, DIMENSION(n) :: temp !temporary matrix
 
c        temp = 0
        do j=1,3
           temp(j)=0
        enddo
        !These two loops to the multiplication
        DO i = 1, n
                DO j = 1, n
                        temp(i) = temp(i) + a(i,j)*b(j)
                END DO
        END DO
        do i=1,3
           b(i)=temp(i)
        enddo
c        b = temp
        return
        end
c END SUBROUTINE MULMATRIX
 
      SUBROUTINE FINDLARGEST(x, n, l)
      real x(3),l
c        IMPLICIT NONE
c        INTEGER, INTENT(IN) :: n
c        REAL, INTENT(IN), DIMENSION(n) :: x
c        REAL, INTENT(INOUT) :: l !Largest value
c        
c        INTEGER :: i
c        !Algorithm is easy
c        !Let the largest number be the first one.
c        !If you find a number larger than it, store this number and then continue
c        l = ABS(x(1))
c        DO i = 2, n
c                IF (ABS(x(i)) > l) l = ABS(x(i))
c        END DO
c
      l=abs(x(1))
      do i=2,n
         if(abs(x(i)).gt.l)l=abs(x(i))
      enddo
      return
      end

                
cEND SUBROUTINE FINDLARGEST














!********************************************************
!* Eigenvalues and eigenvectors of a real square matrix *
!* by Rutishauser's method and inverse iteration method *
!* ---------------------------------------------------- *
!* Reference:                                           *
!*                                                      *
!*   "ALGEBRE Algorithmes et programmes en Pascal       *
!*    de Jean-Louis Jardrin - Dunod BO-PRE 1988."       *
!*    [BIBLI 10].                                       *
!*                                                      *
!*                    F90 Release By J-P Moreau, Paris. *
!* ---------------------------------------------------- *
!* SAMPLE RUN:                                          *
!*                                                      *
!* Input file: elpro.dat                                *
!*                                                      *
!* 5                                                    *
!*   1   2   3  -7    12                                *
!*   2   4   7   3    -1                                *
!*   3   7  10   8     4                                *
!*  -7   3   8  -0.75 -9                                *
!*  12  -1   4  -9    10                                *
!*                                                      *
!* Output to screen:                                    *
!*                                                      *
!*         *** EIGENVALUES AND EIGENVECTORS ***         *
!*              OF A REAL SQUARE MATRIX                 *
!*              BY RUTISHAUSER'S METHOD                 *
!*            AND INVERSE ITERATION METHOD              *
!*                                                      *
!*  Matrix A:                                           *
!*                                                      *
!*   1.0000   2.0000   3.0000  -7.0000  12.000000       *
!*   2.0000   4.0000   7.0000   3.0000  -1.000000       *
!*   3.0000   7.0000  10.0000   8.0000   4.000000       *
!*  -7.0000   3.0000   8.0000  -0.7500  -9.000000       *
!*  12.0000  -1.0000   4.0000  -9.0000  10.000000       *
!*                                                      *
!*                                                      *
!*  Eigenvalue 1:   23.75595                            *
!*                                                      *
!*   Eigenvector:                                       *
!*                                                      *
!*    0.70582  -0.01268   0.13149  -0.52750   1.00000   *
!*                                                      *
!*  Eigenvalue 2:  -10.48655                            *
!*                                                      *
!*   Eigenvector:                                       *
!*                                                      *
!*    0.46717  -0.00795  -0.50783   1.00000   0.26443   *
!*                                                      *
!*  Eigenvalue 3:   0.46335                             *
!*                                                      *
!*   Eigenvector:                                       *
!*                                                      *
!*    0.20881   1.00000  -0.47803  -0.27500  -0.21691   *
!*                                                      *
!*  Eigenvalue 4:   -7.77458                            *
!*                                                      *
!*   Eigenvector:                                       *
!*                                                      *
!*    1.00000  -0.32610   0.20962  -0.14769  -0.81542   *
!*                                                      *
!*  Eigenvalue 5:   0.46335                             *
!*                                                      *
!*   Eigenvector:                                       *
!*                                                      *
!*    0.20881   1.00000  -0.47803  -0.27500  -0.21691   *
!*                                                      * 
!********************************************************
c     PROGRAM Elpro
      subroutine calcola_autovec(A,VX,R,it)
cinteger, parameter :: NMAX = 20
      common /eigen/W,X
      real*8 d1,d2,eps
c    real*8 A(NMAX,NMAX), VX(NMAX,NMAX), R(NMAX)
      real*8 A(3,3),VX(3,3),R(3)
      real*8 W(3),X(3)
c    print *,' '
c    print *,'    *** EIGENVALUES AND EIGENVECTORS ***'
c    print *,'          OF A REAL SQUARE MATRIX'
c    print *,'          BY RUTISHAUSER''S METHOD'
c    print *,'        AND INVERSE ITERATION METHOD'
c	print *,' '
      nmax=3
      n=3
!   read data from input text file
c      write(6,*)'A=',A
      call Data(n, A, eps, d1, d2, m)

      call EPMRI(eps,d1,d2,m,n,A,it,R,VX)

      if(it.eq.-1)then
c    if (it==-1) then
         write(6,*)'  No convergence!'
c      print *,'  No convergence!'
c	  print *,' '
      else if(it.eq.-2)then
         write(6,*)'CASO NAN'
        else
         
	  do j=1, n
c	    write(*,10) j, R(j)
c	    print *,' Eigenvector:'
c	    write(*,20)  (VX(i,j),i=1, n)
c        print *,' '
c		print *,' '
      end do
      end if
c    print *,' '  
c      stop
      return
 10   format(' Eigenvalue ',I1,': ',F20.5)
 20   format(5(F20.5)) 

      END


      Subroutine Data(n, A, eps, d1, d2, m)
c  integer, parameter :: NMAX = 20
c  real*8 A(NMAX,NMAX)
      common /eigen/W,X
      real*8 W(3),X(3)
      real*8 A(3,3)
      real*8 eps,d1,d2
      nmax=3
      n=nmax
c	open(unit=1,file='elpro.dat')
c	open(unit=1,file='input.dat')
c        do i=1,n
c           read(1,*)a(i,1),a(i,2),a(i,3)
c        enddo
c    read(1,*)  n
c	print *,' '
c    print *,' Matrix A:'
c	print *,' '
c    do i=1, n
c	  read(1,*) (A(i,j),j=1,n)
c	  write(*,10) (A(i,j),j=1,n)
c    end do  
c    print *,' ';
c       close(1)
      eps=1.d-10
      d1=1.d-8
      d2=1.d-8
      m=600
c      m=10000
	return
c10  format(5F10.4)
        return
        End

  !***************************************************************
  !* Subroutine DECCRM determines the lower triangular matrix and*
  !* the upper triangukar matrix of Crout's decomposition of a   *
  !* given square real matrix, A.                                *
  !* ----------------------------------------------------------- *
  !* INPUTS:                                                     *
  !*         eps: required precision (double)                    *
  !*          n : size of matrix A (integer)                     *
  !*          A : input matrix (n x n)                           *
  !* OUTPUTS:                                                    *
  !*          it: flag, =0 if method does not apply              *
  !*                    =1 if method is ok.                      *
  !*           U: lower triangular matrix.                       *
  !*           V: upper triangular matrix.                       *
  !***************************************************************
      Subroutine DECCRM(eps, n, A, it, U, V)
      common /eigen/W,X
      real*8 W(3),X(3)
c    integer, parameter :: NMAX = 20
      real*8 a(3,3),u(3,3),v(3,3)
c    real*8 A(NMAX,NMAX), U(NMAX,NMAX), V(NMAX,NMAX)
      real*8 eps, s
      nmax=3

c    if (dabs(A(1,1)) < eps)  then
      if (dabs(A(1,1)).lt.eps)  then
	  it=0
       else
      do i=1, n
	    U(i,1)=A(i,1)
      end do 
      V(1,1)=1.d0
      do j=2, n
	    V(1,j)=A(1,j)/U(1,1)
      end do
c      it=1; k=2
      it=1
      k=2
c      do while (it.ne.0.and.k<=n)
      do while (it.ne.0.and.k.le.n)
	    do i=1, n
c	      if (i < k) then  
	      if (i.lt.k) then  
		    U(i,k)=0.d0
		  else
	        s=0.d0
	        do j=1, k-1
			  s = s + U(i,j)*V(j,k)
            end do
	        U(i,k) = A(i,k) - s
		  end if
        end do
c	    if (dabs(U(k,k)) < eps) then 
	    if (dabs(U(k,k)).lt.eps) then 
          it=0
		else
	      do j=1, n
c	        if (j < k) then 
	        if (j.lt.k) then 
			  V(k,j)=0.d0
c	        else if (j==k) then
	        else if (j.eq.k) then
			  V(k,j)=1.d0
            else
	          s=0.d0
	          do i=1, k-1
			    s = s + U(k,i)*V(i,j)
              end do
	          V(k,j) = A(k,j)/U(k,k);
			end if
          end do
	      k = k + 1
		end if
      end do !while
      end if
	return
        End

!  void MatPrint(char *title, MAT A, int n) {
!  int i,j;
!   printf("\n %s\n", title);
!for (i=1; i<=n; i++) {
!     for (j=1; j<=n; j++)
!       printf(" %f", A[i][j]);
!	  printf("\n");
!   }
! }

! void VecPrint(char *title, VEC X, int n) {
!   int i;
!   printf("\n %s\n", title);
!   for (i=1; i<=n; i++) printf(" %f", X[i]);
!	printf("\n");
! }

  !*********************************************************
  !* Calculate the eigenvalues of a real square matrix by  *
  !* Rutishauser's Method.                                 *
  !* ----------------------------------------------------- *
  !* INPUTS:                                               *
  !*        eps: absolute precision (double)               *
  !*        dta: relative precision (double)               *
  !*         m : maximum number of iterations (integer)    *
  !*         n : size of matrix A (integer)                *
  !*         A : input real square matrix (n x n)          *
  !* OUTPUTS:                                              *
  !*         it: flag, =-1 if convergence is not obtained  *
  !*                   =1 if convergence is ok.            *
  !*         R : contains in output the n eigenvalues of A *
  !*                                                       *         
  !*********************************************************
      Subroutine VAMR(eps, dta, m, n, A, it, R)
      common /eigen/W,X
      real*8 W(3),X(3)
c    integer, parameter :: NMAX = 20
c    real*8 A(NMAX,NMAX), R(NMAX)
      real*8 eps, dta
      real*8 phi,s,t0
c    real*8 U(NMAX,NMAX), V(NMAX,NMAX)
      real*8 a(3,3),r(3),u(3,3),v(3,3)
      nmax=3
      t0=0.d0
      l=1
c      do while (l<=m.and.it.ne.1)
      do while (l.le.m.and.it.ne.1)
      do i=1, n
	    R(i)=A(i,i)
c            write(6,*)'i,R(i)',R(i)
      end do

      call DECCRM(eps, n, A, it, U, V)

c      if (it==0) then
      if (it.eq.0) then
	    do i=1, n
		  A(i,i)=A(i,i) + 1.d0
        end do
	    t0 = t0 + 1.d0
      else
	    do i=1, n
		  do j=1, n
	        s=0.d0
	        do k=1, n
			  s = s + V(i,k) * U(k,j)
            end do
	        A(i,j) = s
		  end do
        end do
	    phi=0.d0
        do i=1, n
	      s= dabs(A(i,i)-R(i))
c	      if (s > phi)  phi=s
	      if (s.gt.phi)  phi=s
		end do
c	    if (phi < dta) then
	    if (phi.lt.dta) then
	      do i=1, n
		    R(i) = A(i,i) - t0
          end do
		else
	      l = l + 1
	      it=-1
		end if
      end if
      end do                    !while
c      write(6,*)'RRR',R
      return
      End

  !************************************************************
  !* Procedure IIM calculates a real eigenvalue and the asso- *
  !* ciated eigenvector of a real square matrix the inverse   *
  !* iteration method.                                        *
  !* -------------------------------------------------------- *
  !* INPUTS:                                                  *
  !*         eps : absolute precision (double)                *
  !*         dta : relative precision (double)                *
  !*          m  : maximum number of iterations (integer)     *
  !*          n  : size of matrix A                           *
  !*          A  : input real square matrix (n x n)           *
  !* OUTPUTS:                                                 *
  !*          it : flag, =-1 if convergence is not obtained   *
  !*                     =1 if convergence is ok.             *
  !*        Gamma: starting value for the eigenvalue as input *
  !*               approximation of the eigenvalue with preci-*
  !*               sion dta in output.                        *
  !*          X1 : contains in output the associated eigen-   *
  !*               vector.                                    *
  !*                                                          *
  !************************************************************
      Subroutine IIM(eps, dta, m, n, A, it, gamma, X1)
      common /eigen/W,X
      real*8 X(3)
c      integer, parameter :: NMAX = 20
c      real*8 A(NMAX,NMAX), X1(NMAX)
      real*8 a(3,3),x1(3)
      real*8 eps, dta, gamma
      real*8  p0,phi,s,t0
c      real*8 W(NMAX), X0(NMAX)
      real*8 w(3),x0(3)
c      integer LP(NMAX)
      integer lp(3)
      nmax=3
      if(gamma.ne.gamma)then
         return
      endif
c      write(6,*)'inizio',W,A
c      write(6,*)'gamma,gamma+1',gamma,gamma+1
c      if(gamma.ne.gamma)write(6,*)'GAMMA is NAN',gamma
      do i=1, n
	  A(i,i) = A(i,i) - gamma
       end do
c       write(6,*)'GAMMA',gamma
       do k=1, n-1
          p0=A(k,k)
          l0=k
c          write(6,*)'P0=',p0
c          write(6,*)'Akk',k,A
c          p0=A(k,k); l0=k
          do i=k+1, n
c		if (dabs(A(i,k)) > dabs(p0)) then
		if (dabs(A(i,k)).gt.dabs(p0)) then
c	      p0=A(i,k); l0=i
	      p0=A(i,k)
              l0=i
c              write(6,*)'P0dentro=',p0
		end if
      end do 
      LP(k)=l0
c      write(6,*)'EPS',eps,p0
c      if (dabs(p0) < eps) then
      if (dabs(p0).lt.eps) then
c	    p0=eps; A(l0,k)=eps
	    p0=eps
            A(l0,k)=eps
      end if

      if (l0.ne.k) then
	    do j=k, n
c	      t0=A(k,j); A(k,j)=A(l0,j); A(l0,j)=t0
	      t0=A(k,j)
              A(k,j)=A(l0,j)
              A(l0,j)=t0
		end do
      end if
	  do i=k+1, n
             if(p0.ne.0.d0)then
	    A(i,k)=A(i,k)/p0
            else
               A(i,k)=0.
c               write(6,*)'AIK'
            endif
c            write(6,*)'AAAAAA','p0=',p0,i,k,A
	    do j=k+1, n
	      A(i,j)=A(i,j)-A(i,k)*A(k,j)
        end do
      end do
	end do !k loop

c    if (dabs(A(n,n)) < eps)  A(n,n)=eps
        if (dabs(A(n,n)).lt.eps)  A(n,n)=eps
        do i=1, n
	  X0(i)=1.d0/sqrt(1.d0*i)
       end do 

c       write(6,*)'AAdopo',A

c    it=-1; l=1
       it=-1
      l=1
c    do while (it==-1.and.l<=m)
      do while (it.eq.-1.and.l.le.m)
      do i=1, n
	    W(i)=X0(i)
      end do
c      write(6,*)'WWa=',W
      do k=1, n-1
        l0=LP(k)
        if (l0.ne.k) then
c          t0=W(k); W(k)=W(l0); W(l0)=t0
          t0=W(k)
          W(k)=W(l0)
          W(l0)=t0
        end if
c        write(6,*)'l0,Lp',k,l0,LP(k),t0
        do i=k+1, n
		  W(i)=W(i)-A(i,k)*W(k)
c                  write(6,*)'WiAi',i,k,W(i),W(k),A(i,k)
        end do
      end do
c      write(6,*)'WWW=',W
      if(A(n,n).ne.0.d0)then
      X1(n)=W(n)/A(n,n)
      else
         X1(n)=0
c         write(6,*)'Ann',A(n,n)
      endif

c      write(6,*)'X1=',X1,W
      do i=n-1, 1, -1
        s=0.d0
        do j=i+1, n
		  s = s + A(i,j)*X1(j)
        end do
        if(A(i,i).ne.0.d0)then
	    X1(i)=(W(i)-s)/A(i,i)
        else
c           write(6,*)'Aii=0',A(i,i)
           X1(i)=0
        endif
      end do
      p0=0.d0
      do i=1, n
c        if (dabs(X1(i)) > dabs(p0))  p0=X1(i)
        if (dabs(X1(i)).gt.dabs(p0))  p0=X1(i)
      end do
c      if(p0.eq.0.d0)write(6,*)'X1(i)',X1
      do i=1, n
         if(p0.ne.0.d0)then
        X1(i)=X1(i)/p0
        else
c           write(6,*)'p0=0',p0
           X1(i)=0
        endif
      end do
      phi=0.d0
      do i=1, n
        s=dabs(X1(i)-X0(i))
c        if (s > phi)  phi=s
        if (s.gt.phi)  phi=s
      end do

c      if (phi < dta) then
      if (phi.lt.dta) then
         if(p0.ne.0.d0)then
        gamma = gamma + 1.d0/p0
        else
c           write(6,*)'p0=0bis',p0
           gamma=gamma
        endif
        it=1
      else
        do i=1, n
		  X0(i)=X1(i)
        end do 
        l = l + 1
      end if
	end do !while
	return
        End                     !IIM  	  						    


  !******************************************************
  !* INPUTS:                                            *
  !* EPS : precision (Double)                           *
  !* D1  : precision d1 (Double)                        *
  !* D2  : precision d2 (Double)                        *
  !* M   : maximum number of iterations (integer)       *
  !* N   : order of matrix A (integer)                  *
  !* A   : input matrix to study (of MAT type)          *
  !* -------------------------------------------------- *
  !* OUTPUTS:                                           *
  !* IT  : -1 if no convergence is obtained (integer)   *
  !* R   : table of eigenvalues (of VEC type)           *
  !* VX  : table of eigenvectors (of MAT type)          *
  !******************************************************
      Subroutine EPMRI(eps, d1, d2, m, n, A, it, R, VX)
      common /eigen/W,X
c    integer, parameter :: NMAX = 20
c    real*8 A(NMAX,NMAX), R(NMAX), VX(NMAX,NMAX)
c    real*8 X(NMAX), A1(NMAX,NMAX)
      real*8 W(3)
      real*8 a(3,3),r(3),vx(3,3)
      real*8 x(3),a1(3,3)
      real*8 eps,d2,d1
      nmax=20
      do i=1, n
      do j=1, n
        A1(i,j) = A(i,j)
      end do
      end do

      call VAMR(eps,d2,m,n,A1,it,R)
c      write(6,*)R
      if(R(1).ne.R(1).or.R(2).ne.R(2).or.R(3).ne.R(3))then
c         write(6,*)'CASO NAN'
         it=-2
         return
      endif
! restore A1 after VAMR
      do i=1, n
      do j=1, n
        A1(i,j) = A(i,j)
      end do
	end do   

        j=1

c        X(1)=0 !aggiunta
c        X(2)=0 !aggiunta
c        X(3)=0!aggiunta

c    do while (it==1.and.j<=n)
        do while (it.eq.1.and.j.le.n)

           call IIM(eps,d1,m,n,A1,it,R(j),X)
          
! restore A1 after IIM
      do i=1, n
        do k=1, n
          A1(i,k) = A(i,k)
        end do
        VX(i,j) = X(i)
      end do
      j = j + 1     
      end do                    !while
	
	return
        End

! End of file elpro.f90


c      subroutine correlazioni
c      common /varcm/velolablcp(100,3),velocmlcp(100,3),thetacmlcp(100),
c     c                   phicmlcp(100),vcmlcp(100)
c            common /qp/jmaxsicsi,zmaxsicsi
c      dimension vlistap(5,200,3),vlistaa(5,200,3)
c      dimension inda(50),indp(50)
c      np=0
c      na=0
c      do j=1,npart
c         if(zpart(j).eq.1.and.apart(j).eq.1)then
c            np=np+1
c            indp(np)=j
c         endif
c         if(zpart(j).eq.2.and.apart(j).eq.4)then
c            na=na+1
c            inda(na)=j
c         endif
c      enddo
c      return
c      end

      FUNCTION SCAPROD(U,V)                                             
      DIMENSION U(3),V(3)                                               
      SCAPROD=0.                                                        
      DO I=1,3                                                          
         SCAPROD=SCAPROD+ U(I)*V(I)                                     
      END DO                                                            
      RETURN                                                            
      END        
