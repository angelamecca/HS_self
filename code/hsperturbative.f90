!
!=======================================================================
!
! effective mass mstar
!
!=======================================================================
!

double precision function  mstar_pert(ki,kF,a,deg)
  implicit none
  integer nk
  double precision ki, mstar
  double precision kF, a, deg
  double precision c, x, Fx, one, two, three, four, pi, mstarpert
  !
  pi = acos(-1.d0)
  c = kF*a
  x = ki/kF
  if(x.eq.1.d0) then
     mstarpert = 8.d0*(deg-1.d0)/(15.d0*pi**2)* (7.d0*log(2.d0) -1.d0)*c**2
     mstar_pert = 1.d0 + mstarpert
     return
  end if
  
  if (x.le.sqrt(2.d0)) then
     one = 4.d0/x
     two = 8.d0*x**3*log(x**2/abs(x**2-1.d0))
     three= -10.d0*(1.d0 + 1.d0/x**2)* log(abs((x+1.d0)/(x-1.d0)))
     four = (1.d0 + x*(2.d0-x**2)**0.5)/(1.d0 - x*(2.d0-x**2)**0.5)
     four = 2.d0*(2.d0+1/x**2)*(2.d0-x**2)**1.5*log(abs(four))
  else 
     one = (-2.d0*x**5 + 4.d0*x**4 + 8.d0*x**3 -20.d0*x**2 -8.d0*x + 20.d0)
     one = one/x/(1.d0-x**2)
     two= x**2/(x**2-1.d0)
     two= 8.d0*x**3*log(two)
     three= (x+1.d0)/(x-1.d0)
     three= -10.d0*(1.d0 + 1/x**2)*log(three)
     four= (x**2-2.d0)**(-0.5)
     four = atan(four)
     four= -4.d0*(2.d0+1/x**2)*(x**2-2.d0)**1.5*four
  end if
  Fx= one + two + three + four
  mstar = 1.d0 - 2.d0*(deg-1.d0)/15.d0/pi**2 * c**2 * fx/x
  mstar_pert = mstar
  return
end function mstar_pert
!
!=======================================================================
!
!momentum distribution nk
!
!=======================================================================
!
subroutine nk_pert(ki,nk,kF,a,deg,rho1,rho2)
  implicit none
  integer nk
  double precision ki(nk), rho1(nk),rho2(nk)
  double precision deg,a,kF
  double precision pi, c, x
  integer k
  double precision l1, l2, l3, all
  double precision t1,t2,t3
  double precision one, two, three, tot
  pi = 4.*atan(1.d0)
  c = kF*a
  
  do k = 1, nk
     x = ki(k)/kF         
     if(x.lt.1.d0) then
        one = (7.d0*log(2.d0)-8.d0)*x**3 + (10.d0 -3.0*log(2.0))*x
        two = 2.d0* log( (1.d0+x)/(1.d0-x))- 2.d0*(2.0d0-x**2)**1.5* &
             log(((2.d0-x**2)**0.5 + x)/( (2.d0-x**2)**0.5 - x))
        tot= 1.d0 - (deg-1.d0)/(3.0*pi**2*x)*c**2*(one+two)
     else if(x.eq.1.0) then
        tot = 1- 2.d0/(3.d0*pi**2)*(deg-1.d0)*c**2*(3.d0*log(2.0)+1.d0)
     else
        tot =0.d0
     end if
     rho1(k) = tot
  end do
  
  do k = 1, nk
     x = ki(k)/kF         
     if(x.lt.1.d0) then
        tot = 0.d0
     else if(x.eq.1) then
        tot =  2.d0/(3.d0*pi**2)*(deg-1.d0)*c**2* (3.d0*log(2.d0)-1.d0)
     else if(c.lt.3.0) then
        one = (7.d0*x**3-3.d0*x-6.d0)*log((x-1.d0)/(x+1.d0))
        two = (7.d0*x**3-3.d0*x+2.d0)*log(2.d0)-8.d0*x**3+22.d0*x**2 + 6.d0*x-24.d0           
        if(x.lt.sqrt(2.d0)) then
           l1 = (2.0d0+x+(2.d0-x**2)**0.5)/(2.d0+x-(2.d0-x**2)**0.5)
           l1 = log(l1)
           l2 = (1.d0+(2.d0-x**2)**0.5)/(1.d0-(2-x**2)**0.5)
           l2 = log(l2)
           l3 = (x+ (2.d0-x**2)**0.5)/(x-(2.d0-x**2)**0.5)
           l3 = -2.d0*log(l3)
           all= l1+l2+l3
           three= 2.d0*(2.d0-x**2)**1.5* all
        else 
           t1 = (x+2.d0)/((x**2-2.d0)**0.5)
           t1 = atan(t1)
           t2 = (x**2-2.d0)**(-0.5)
           t2 = atan(t2)
           t3 = x*(x**2-2.d0)**(-0.5)
           t3 = -2.d0*atan(t3)
           all = t1+t2+t3
           three = -4.d0*(x**2-2.d0)**1.5* all
        end if
        tot = (deg-1.d0)/(6.d0*pi**2*x)*c**2*(one+two+three)
     else
        t1 = x*(x**2-2.)**(-0.5)
        t1 = 2.d0*atan(t1)
        t2 = (x-2.d0)*(x**2-2.d0)**(-0.5)
        t2 = -atan(t2)
        t3 = (x+2.d0)*(x**2-2.d0)**(-0.5)
        t3 = -atan(t3)
        all = t1+t2+t3
        one = 2.d0*log((x+1.d0)/(x-1.d0))
        two = -2.d0*x
        three = (x**2-2.d0)**1.5 *all
        tot = 2.d0*(deg-1.d0)/(3.d0*pi**2*x)*c**2*(one+two+three)
     end if
     rho2(k) = tot
  end do
  
end subroutine nk_pert


