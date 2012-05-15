!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

function rslf(p,t)

!     This function calculates the liquid saturation vapor mixing ratio as
!     a function of pressure and Kelvin temperature

implicit none
real esl,rslf,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)

x=max(-80.,t-273.16)

esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rslf=.622*esl/(p-esl)

return
end

!     ******************************************************************

function rsif(p,t)

!     This function calculates the ice saturation vapor mixing ratio as a
!     function of pressure and Kelvin temperature

implicit none
real esi,rsif,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01)
parameter (c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05)
parameter (c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12)

x=max(-80.,t-273.16)
esi=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rsif=.622*esi/(p-esi)

return
end

!     ******************************************************************

function eslf(t)

!     This function calculates the liquid saturation vapor pressure as a
!     function of Celcius temperature

implicit none
real eslf,x,t,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)

x=max(-80.,t)
eslf=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************

function esif(t)

!     This function calculates the ice saturation vapor pressure as a
!     function of Celsius temperature

implicit none
real esif,x,t,c0,c1,c2,c3,c4,c5,c6,c7,c8
parameter (c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01)
parameter (c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05)
parameter (c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12)

x=max(-80.,t)
esif=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end

!     ******************************************************************

function eslpf(t)

!     This function calculates the partial derivative of liquid saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature

implicit none
real eslpf,x,t,d0,d1,d2,d3,d4,d5,d6,d7,d8
parameter (d0= .4443216e+02,d1= .2861503e+01,d2= .7943347e-01)
parameter (d3= .1209650e-02,d4= .1036937e-04,d5= .4058663e-07)
parameter (d6=-.5805342e-10,d7=-.1159088e-11,d8=-.3189651e-14)

x=max(-80.,t)
eslpf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end

!     ******************************************************************

function esipf(t)

!     This function calculates the partial derivative of ice saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature

implicit none
real esipf,x,t,d0,d1,d2,d3,d4,d5,d6,d7,d8
parameter (d0= .5036342e+02,d1= .3775758e+01,d2= .1269736e+00)
parameter (d3= .2503052e-02,d4= .3163761e-04,d5= .2623881e-06)
parameter (d6= .1392546e-08,d7= .4315126e-11,d8= .5961476e-14)

x=max(-80.,t)
esipf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end

!     ******************************************************************

!     This function calculates the partial derivative of liquid saturation vapor
!     mixing ratio with respect to temperature as a function of pressure and
!     Kelvin temperature

function rslfp(p,t)

implicit none
real rslfp,eslpf,x,t,d0,d1,d2,d3,d4,d5,d6,d7,d8,p
parameter (d0= .4443216e+02,d1= .2861503e+01,d2= .7943347e-01)
parameter (d3= .1209650e-02,d4= .1036937e-04,d5= .4058663e-07)
parameter (d6=-.5805342e-10,d7=-.1159088e-11,d8=-.3189651e-14)

x=max(-80.,t-273.16)
eslpf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
rslfp=.622*eslpf/(p-eslpf)

return
end

!     ******************************************************************

!     This function calculates the partial derivative of ice saturation vapor
!     mixing ratio with respect to temperature as a function of pressure and
!     Kelvin temperature

function rsifp(p,t)

implicit none
real rsifp,esipf,x,t,d0,d1,d2,d3,d4,d5,d6,d7,d8,p
parameter (d0= .5036342e+02,d1= .3775758e+01,d2= .1269736e+00)
parameter (d3= .2503052e-02,d4= .3163761e-04,d5= .2623881e-06)
parameter (d6= .1392546e-08,d7= .4315126e-11,d8= .5961476e-14)

x=max(-80.,t-273.16)
esipf=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
rsifp=.622*esipf/(p-esipf)

return
end

!     ******************************************************************

subroutine mrsl(n1,p,t,rsl)

implicit none
integer n,n1
real rsl(n1),rslf,p(n1),t(n1)

do n=1,n1
   rsl(n)=rslf(p(n),t(n))
enddo

return
end

!     ******************************************************************

subroutine mrsi(n1,p,t,rsi)

implicit none
integer n,n1
real rsi(n1),rsif,p(n1),t(n1)

do n=1,n1
   rsi(n)=rsif(p(n),t(n))
enddo

return
end

!     ******************************************************************

subroutine thvtoth(nn,theta,rv,rtp)

implicit none
integer nn,k
real theta(nn),rv(nn),rtp(nn)

do k=1,nn
  theta(k)=theta(k)*(1.+rtp(k))/(1.+1.61*rv(k))
enddo

return
end

!     ******************************************************************

function td(p,rs)

implicit none
real rr,rs,es,esln,p,td

rr=rs+1e-8
es=p*rr/(.622+rr)
esln=log(es)
td=(35.86*esln-4947.2325)/(esln-23.6837)

return
end

!     ******************************************************************

FUNCTION RS(P,T)

ES=610.78*EXP(17.269*(T-273.16)/(T-35.86))
RS=.622*ES/(P-ES)

RETURN
END

!     ******************************************************************

SUBROUTINE THETAE(P,T,RV,THE)
DATA CP/1004./,G/9.8/,R/287./,ALVL/2.35E6/

CPG=CP/G
PIT=P
TUPO=T
TTD=TD(P,RV)
DZ=CPG*(T-TTD)
IF(DZ.LE.0.) GOTO 20
DO ITTER=1,50
   TUPN=T-G/CP*DZ
   TMN=(TUPN+T)*.5*(1.+.61*RV)
   PIT=P*EXP(-G*DZ/(R*TMN))
   IF(ABS(TUPN-TUPO).LT.0.001) GOTO 20
   TTD=TD(PIT,RV)
   TUPO=TUPN
   DZ=DZ+CPG*(TUPN-TTD)
ENDDO
STOP 10
20 CONTINUE
THE=TUPO*(1E5/PIT)**.286*EXP(ALVL*RV/(CP*TUPO))

RETURN
END

!     ******************************************************************

FUNCTION  TW( RVP,THET,P)

!     ABS IS ABSOLUTE VALUE
!     ALL ARGUMENTS AND TW (KELVIN)

PRESS=P*1.E-2
RVAP=RVP*1.E3
PITER =  PRESS
DO ID=1,10
   TEMPER=THET*(PITER*1.E-3)**.286
   X  =  .02*( TMR(RVAP,PITER) - TEMPER)
   IF( ABS(X).LT. 0.01  ) GOTO 5
   PITER = PITER* ( 2.**(X)  )
ENDDO
5    TEMPER=THET*(PITER*1.E-3)**.286

AOS  =   OS(TEMPER,PITER)
TW   =  TSA( AOS,PRESS)

RETURN
END

!     ******************************************************************

FUNCTION   OS(T,P)

!     OS AND T (KELVIN) , P (MILLIBARS )

OS=T*((1000./P)**.286)/(EXP(-2.6518986*W(T,P)/T))

RETURN
END

!     ******************************************************************

FUNCTION TSA(OS,P)

!     TSA AND OS(KELVIN),P(MILLIBARS)
!     SIGN(A,B) RREPLACES THE ALGEBRETIC SIGN OF A WITH THE SIGN OF B

A  =  OS
TQ =  253.16
D  =  120
DO ID= 1,12
   D = D/2.
!     IF THE TEMPERATURE DIFFERENCE,X, IS SMALL,EXIT THIS LOOP
   X=A*EXP(-2.6518986*W(TQ,P)/TQ)-TQ*((1000./P)**.286)
   IF(ABS(X).LT.0.01) GOTO 2
   TQ = TQ + SIGN(D,X)
ENDDO
2    TSA=TQ

RETURN
END

!     ******************************************************************

FUNCTION  ESAT(T)

!     ESAT(MILLIBARS),T(KELVIN)

DATA ABZ/273.16/

TC=T-ABZ
ESAT=6.1078*EXP((17.2693882*TC)/(TC+237.3))

RETURN
END

!     ******************************************************************

FUNCTION  W(T,P)

!     W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )

IF(T.GE.999.) GOTO 10
X  =  ESAT(T)
W  =  622.*X/(P-X)
RETURN
10 W=0.0

RETURN
END

!     ******************************************************************

FUNCTION  TMR(W,P)

!     TMR(KELVIN),W(GRAMS WATER VAPOR/KILOGRAM DRY AIR),P(MILLIBAR)
!     LOG10  15   LOG TO THE BASE  TEN.

X =  LOG10(   W*P/(622.+ W)  )
TMR=10.**(.0498646455*X+2.4082965)-7.07475+38.9114*((10.**(  &
  .0915*X ) - 1.2035 )**2 )

RETURN
END

!     ******************************************************************

SUBROUTINE THE2T(THE,P,TH,T,R)

DATA CP/1004./,ALVL/2.350E6/

PI=(P*1E-5)**.286
TO=THE/EXP(ALVL*.012/(CP*295.))*PI
DO ITTER=1,50
   R=RS(P,TO)
   TH=THE/EXP(ALVL*R/(CP*TO))
   TN=TH*PI
   IF(ABS(TO-TN).LT.0.005) GOTO 12
   TO=TO+(TN-TO)*.3
ENDDO
WRITE(6,1) THE,P,TO,TN,TH,R
1 FORMAT(' STOP IN ROUTINE THE2T '/' THE,P,TO,TN,TH,R',6E15.6)
STOP 10
12 CONTINUE
T=TN

RETURN
END

!     ******************************************************************

subroutine qtk(q,tempk,fracliq)
parameter (r4186=1./4186.,r2093=1./2093.,r334000=1./334000.)

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!       tempk    temperature [K]
!       fracliq  liquid fraction [dimensionless]
!     Local Constants:
!       4186     specific heat of liquid [J/(kg K)]
!       2093     specific heat of ice [J/(kg K)]
!       334000   latent heat of fusion [J/kg]
!       273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempk = q * r2093 + 273.15
elseif (q .ge. 334000.) then
   fracliq = 1.
   tempk = q * r4186 + 193.36
else
   fracliq = q * r334000
   tempk = 273.15
endif

return
end

!     ******************************************************************

subroutine qtc(q,tempc,fracliq)
parameter (r4186=1./4186.,r2093=1./2093.,r334000=1./334000.)

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!        tempc    temperature [C]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempc = q * r2093
elseif (q .ge. 334000.) then
   fracliq = 1.
   tempc = q * r4186 - 80.
else
   fracliq = q * r334000
   tempc = 0.
endif

return
end

!     ******************************************************************

subroutine qwtk(qw,w,dryhcap,tempk,fracliq)
parameter (r4186=1./4186.,r2093=1./2093.,r334000=1./334000.)

!     Inputs:
!        qw       internal energy [J/m^2] or [J/m^3]
!        w        mass [kg/m^2] or [kg/m^3]
!        dryhcap  heat capacity of nonwater part [J/(m^2 K)] or [J/(m^3 K)]
!     Outputs:
!        tempk    temperature [K]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

qwliq0 = w * 334000.
if (qw .le. 0.) then
   fracliq = 0.
   tempk = qw / (2093. * w + dryhcap) + 273.15
elseif (qw .ge. qwliq0) then
   fracliq = 1.
   tempk = (qw - qwliq0) / (4186. * w + dryhcap) + 273.15
else
   fracliq = qw / qwliq0
   tempk = 273.15
endif
return
end

