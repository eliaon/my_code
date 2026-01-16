      Program Jpsi_tot_gbw

      IMPLICIT DOUBLE PRECISION(A-H,L-Z)
      common/dados/alfem,pi
      common/q2/q2
      common/xbjorken/xbj
      external ampl,lnampl!,gammln
      open (unit=22,file='jpsi_gbw_gammap_bg.dat',
     # status='unknown')

      pi=3.1416d0
      alfem=1.d0/137.d0
      q2 = 0.d0  !não adianta fazer q2 diferente de 0, pois o programa 
		 !só considera contribuição da polarização transversal



	mv = 3.097d0 !J/Psi



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	WA = mv
	WB = 200.d0
	WC = 9000.d0

	W = WA

	j1 = 80
	aj1 = j1

	j2 = 40
	aj2 = j2

	step1 = (WB - WA) / aj1

	step2 = (WC - WB) / aj2

	j = j1 + j2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	


      do 30 i=0,j

	xbj = (q2 + mv**2.d0)/(w*w)

	amp = ampl(xbj)
	amp2 = amp**2.d0

ccccccccccccccccccccccccccccccccccccccccccccccc
c	valores de Bv
ccccccccccccccccccccccccccccccccccccccccccccccc

	b0 = 4.99d0	! GeV^-2
	W0 = 90.d0	! GeV
	al = 0.25d0	! GeV^-2
	bv = b0 + 4.d0 * al * dlog(W / W0)
cccccccccccccccccccccccccccccccccccccccccccccc
	
	y = dlog(1.d0/xbj)

	dsigtdt = amp2/(16.d0*pi)

      sigtot = (dsigtdt/bv)*0.389d0! mb


      write(*,*)rg2,w,sigtot*1.D6,i
      write(22,100)w,sigtot*1.D6! nb



	if(W .lt. 200.d0) then

	W = W + step1 !1.20d0

	else

	W = W + step2 !94.d0

	end if


 30      continue
 100     format(2x,6(E10.4, 2x))
      end
c
c===================================================================
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC    
CCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       double precision function ampl(xb)
       IMPLICIT DOUBLE PRECISION (A-H,L-Z)
       common/dados/alfem,pi
       common/q2/q2
       common/xbjorken/xbj
       external rint
	lir=0.d0
	lsr=50.d0!100.d0
	ampl = sgs0(lir,lsr,0.001d0,rint)
       return
       end



      DOUBLE PRECISION function rint(rd)
      IMPLICIT DOUBLE PRECISION(A-H,L-Z)
      common/dados/alfem,pi
      common/q2/q2
      common/xbjorken/xbj
      common/rd2/rd2
      external zint
	liz=0.d0
	lsz=1.d0
	rd2=rd*rd
	rint=2.d0*pi*rd*sgs0(liz,lsz,0.001d0,zint)
      return
      end


      DOUBLE PRECISION function zint(zz)
      IMPLICIT DOUBLE PRECISION(A-H,L-Z)
      common/dados/alfem,pi
      common/q2/q2
      common/xbjorken/xbj
      common/rd2/rd2
      external wfgt,sigdipgbw

	rdt = dsqrt(rd2)
	zzt = zz
	q2t = q2
        prodwf = wfgt(rdt,zzt,q2t)/(4.d0*pi) !divisão por (4*pi)

        sdip = sigdipgbw(rd2) ! sigiim(rd2)!!!

        zint=prodwf*sdip
      return
      end



      DOUBLE PRECISION function SIGDIPGBW(r2)
      IMPLICIT DOUBLE PRECISION(A-H,L-Z)
      common/dados/alfem,pi
      common/q2/q2
      common/xbjorken/xbj
       sig0=23.0d0/0.3894d0 !23.03d0/0.389d0 ! GeV^-2
       lambexpgw=0.29d0   !0.288d0
       x0= 3.0d-4   !3.04d-4
       qs2gw=((x0/xbj)**lambexpgw)
       xx=r2*qs2gw
       nq = 1.d0-dexp(-(1.d0/4.d0)*xx)  !(1.d0/4.d0)*xx!
       sigdipgbw=sig0*nq*((1.d0-xbj)**5.d0)
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     produto funcoes de onda foton - meson
C
C     Contribuicao transversa
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     Calcula a função de overlap méson-fóton
c     Modelo Boosted - Gaussian do artigo
c     PRD 74, 074016 (2006) - Kowalski
c
c     Esta função de overlap tem que ser dividida
c     por (4*pi) antes de ser inserida nas fórmulas
c     das seções de choque, pois Kowalski usou 
c     uma normalização diferente para elas, 
c     multiplicando a função original por um fator 
c     (4*pi).
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

       FUNCTION WFGT(rdt,zzt,q2t)
       IMPLICIT DOUBLE PRECISION (A-H,L-Z)
       common/dados/alfem,pi


c	Meson jpsi 
	Mf2 = (1.4d0)**2.d0      !massa do quark (GeV)
	ehat = 2.d0 / 3.d0
	NT = 0.578d0
	R2 = 2.3d0      !GeV-2	


	Nc = 3.d0
	ANORM = ehat*(dsqrt(4.d0*pi*alfem))*Nc / (pi*zzt*(1.d0-zzt))

	EPS2 = zzt*(1.D0-zzt)*q2t + Mf2
	EPS = DSQRT(EPS2)
	mK0 = dbesk0(EPS*rdt)
	mK1 = dbesk1(EPS*rdt)

C####################################################################
	PHITpart1 = NT*zzt*(1.d0-zzt)
	
	ARGpart1 =  -(Mf2*R2)/(8.d0*zzt*(1.D0-zzt))
	
	ARGpart2 = -(2.d0*zzt*(1.D0-zzt)*rdt*rdt)/(R2)
	
	ARGpart3 = (Mf2*R2)/(2.d0)
	
	ARGexp = ARGpart1 + ARGpart2 + ARGpart3
	
	PHIT = PHITpart1*dexp(ARGexp)

	ZZZ = zzt**2.d0 + (1.d0-zzt)**2.d0

	DPHIT = -((4.d0*rdt*zzt*(1.D0-zzt))/(R2))*PHIT
C####################################################################

        WF = ANORM*((Mf2*mK0*PHIT)-(ZZZ*EPS*mK1*DPHIT))
	WFGT = wf


	RETURN
        END






CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Rotina de Integracao
c
c     sgs0(li,ls,precisao,integrando)
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION SGS0 (A,B,EPS,F )
      DOUBLE PRECISION A,B,S,U,V,SF,C,SL,SG
      DOUBLE PRECISION SP,SA,F,SGS8,EPS,abb
      EXTERNAL F
      S = 0.d0
      U = A
  1   V = B
      IF ( U .LT. B ) THEN
      SF = SGS8 ( F,U,V )
  2   C  = (U+V)/2
      SL = SGS8 ( F,U,C )
      SG = SGS8 ( F,C,V )
      SP = SL+SG
      abb=abs(sf)
      if(abb.ne.0.0) go to 5
      abb=1.
  5   SA = ABS(SP-SF)/(abb*EPS)
      IF (  SA.GE.1.0 ) THEN
      V  = C
      SF = SL
      GOTO 2
      END IF
      U = V
      S = S+SP
      GOTO 1
      END IF
      SGS0=S
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SGS8 ( F,A,B )
      DOUBLE PRECISION A,B,H,S,C,X,F
      EXTERNAL F
      C = (A+B)/2
      H = (B-A)/2
      X = .96028985E0*H
      S = .10122853E0*(F(C+X)+F(C-X))
      X = .79666647E0*H
      S = S + .22238103E0*(F(C+X)+F(C-X))
      X = .52553240E0*H
      S = S + .31370664E0*(F(C+X)+F(C-X))
      X = .18343464E0*H
      S = S + .36268378E0*(F(C+X)+F(C-X))
      SGS8 = S * H
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cc
cc
cc
cc
cc
c	****************************************************************
c
c	*********           MODIFIED BESSEL FUNCTIONS         **********
c
c	****************************************************************
c
c
c
c	%%%%%%%%%	DOUBLE PRECISION K0 BESSEL FUNCTION	%%%%%%%%
c
c
      FUNCTION dbesk0(x)
      DOUBLE PRECISION bessk0,x,dbesk0
CU    USES bessi0
      DOUBLE PRECISION bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.d0) then
        y=x*x/4.d0
        bessk0=(-dlog(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
	dbesk0=bessk0
      else
        y=(2.d0/x)
        bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
	dbesk0=bessk0
      endif
      return
      END

      FUNCTION bessi0(x)
      double precision bessi0,x
      double precision ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
      endif
      return
      END
c
c
c
c	%%%%%%%%%	DOUBLE PRECISION K1 BESSEL FUNCTION	%%%%%%%%
c
c
      FUNCTION dbesk1(x)
      DOUBLE PRECISION bessk1,x,dbesk1
CU    USES bessi1
      DOUBLE PRECISION bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.d0) then
        y=x*x/4.d0
        bessk1=(log(x/2.d0)*bessi1(x))+(1.d0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     *(p5+y*(p6+y*p7))))))
	dbesk1=bessk1
      else
        y=2.d0/x
        bessk1=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
	dbesk1=bessk1
      endif
      return
      END

      FUNCTION bessi1(x)
      double precision bessi1,x
      double precision ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     *0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     *-0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1,
     *0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi1=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
        if(x.lt.0.)bessi1=-bessi1
      endif
      return
      END





