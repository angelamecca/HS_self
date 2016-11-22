      double precision function fint(r,f,nr,hr,x,ichoic)
      implicit real*8 (a-h,o-z)
C                                                                       OMF00030
C     ICHOIC=1   INTERPOLATES IN X THE FUNC. F GIVEN IN NR POINTS R     OMF00040
C     ICHOIC=2 CALCULATES THE DERIVATIVE OF F IN X                      OMF00050
C     ICHOIC=3 CALCULATES THE ABSCISSA CORRESPONDING TO THE VALUE X     OMF00060
C               THE FUNCTION F (TWO-POINT FORMULA)                      OMF00070
C                                                                       OMF00080
      DIMENSION R(1:*),F(1:*)                                           OMF00090
C                                                                       OMF00100
      HRS=HR*HR                                                         OMF00110
      HRC=HRS*HR                                                        OMF00120
      I0=1                                                              OMF00130
      IF(ICHOIC.EQ.3) GOTO 400                                          OMF00140
C                                                                       OMF00150
      NR1=NR-1                                                          OMF00160
      IF(NR1)2,2,6                                                      OMF00170
    2 IF(ICHOIC-2)3,4,400                                               OMF00180
    3 FINT=F(1)                                                         OMF00190
      RETURN                                                            OMF00200
    4 FINT=0.                                                           OMF00210
      RETURN                                                            OMF00220
    6 IF(NR1.GT.3) GOTO 100                                             OMF00230
      IF(NR1-2) 10,20,30                                                OMF00240
   10 IF(ICHOIC-2) 11,12,400                                            OMF00250
   11 FINT=((X-R(1))*F(2)-(X-R(2))*F(1))/HR                             OMF00260
      RETURN                                                            OMF00270
   12 FINT=(F(2)-F(1))/HR                                               OMF00280
      RETURN                                                            OMF00290
C                                                                       OMF00300
   20 A1=X-R(1)                                                         OMF00310
      A2=X-R(2)                                                         OMF00320
      A3=X-R(3)                                                         OMF00330
      F2=2.*F(2)                                                        OMF00340
      IF(ICHOIC-2) 21,22,400                                            OMF00350
   21 FINT=0.5*(A2*A3*F(1)-A1*A3*F2+A1*A2*F(3))/HRS                     OMF00360
      RETURN                                                            OMF00370
   22 FINT=0.5*(F(1)*(A2+A3)-F2*(A1+A3)+F(3)*(A1+A2))/HRS               OMF00380
      RETURN                                                            OMF00390
C                                                                       OMF00400
   30 A1=X-R(I0)                                                        OMF00410
      A2=X-R(I0+1)                                                      OMF00420
      A3=X-R(I0+2)                                                      OMF00430
      A4=X-R(I0+3)                                                      OMF00440
      F1=F(I0)/3.                                                       OMF00450
      F4=F(I0+3)/3.                                                     OMF00460
      IF(ICHOIC-2) 31,32,32                                             OMF00470
   31 FINT=0.5*(A3*A4*(-A2*F1+A1*F(I0+1))+A1*A2*(-A4*F(I0+2)            OMF00480
     1    +A3*F4))/HRC                                                  OMF00490
      RETURN                                                            OMF00500
   32 FINT=0.5*(A3*A4*(-F1+F(I0+1))+(A3+A4)*(-A2*F1+A1*F(I0+1))         OMF00510
     1    +A1*A2*(-F(I0+2)+F4)+(A1+A2)*(-A4*F(I0+2)                     OMF00520
     2    +A3*F4))/HRC                                                  OMF00530
      RETURN                                                            OMF00540
C                                                                       OMF00550
  100 N=ABS(X-R(1))/HR                                                  OMF00560
      N2=N+2                                                            OMF00570
      I0=N                                                              OMF00580
      IF(N.LE.0) I0=1                                                   OMF00590
      IF(N2 - NR) 30,110,120                                            OMF00600
  110 I0=N-1                                                            OMF00610
      GOTO 30                                                           OMF00620
  120 I0=NR-3                                                           OMF00630
      GOTO 30                                                           OMF00640
C                                                                       OMF00650
C     INVERSE INTERPOLATION - X IS NOW A VALUE OF F                     OMF00660
C                                                                       OMF00670
  400 SOLD=1.                                                           OMF00680
      IF(X.LE.F(1)) SOLD=-1.                                            OMF00690
      DO 420 I=2,NR                                                     OMF00700
      SNEW=1.                                                           OMF00710
      IF(X.LE.F(I)) SNEW=-1.                                            OMF00720
      ISS=SNEW*SOLD                                                     OMF00730
      IF(ISS) 410,410,415                                               OMF00740
  410 FINT=R(I-1)+HR*(X-F(I-1))/(F(I)-F(I-1))                           OMF00750
      GOTO 510                                                          OMF00760
  415 SOLD=SNEW                                                         OMF00770
  420 CONTINUE                                                          OMF00780
  500 WRITE(6,1000 ) X,F(1),F(NR)                                       OMF00790
 1000 FORMAT (5X,' POINTS OUT OF RANGE '/5X,'X,F(1),F(N)=',3E15.6)      OMF00800
      STOP                                                              OMF00810
  510 RETURN                                                            OMF00820
      END                                                               OMF00830
