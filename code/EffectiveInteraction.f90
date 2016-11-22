!
!=======================================================================
!
MODULE Parameters
  implicit none
  double precision kF, EkF, deg, a, d, rmax,  m, chempot
  double precision htm, hbarc, pi
end MODULE Parameters
!
!=======================================================================
!



MODULE EffectiveInteraction
  double precision, dimension(:), allocatable:: r,fc,dfc,ddfc,veff
  double precision hr
  integer na, nd, nr, nt
  integer np
  double precision dp, pcut
  double precision, dimension(:), allocatable:: p, Mp, dMp, ddMp
  
contains


!
!=======================================================================
!

  double precision function T_2bc(na,nd,ntot,h,r,deg,m,htm,rho, sl,d_fr)
!.....total energy
!
    implicit none
    double precision rho, deg
    integer na, nd, ntot
    double precision  sl(ntot), d_fr(ntot)
    double precision  ha, hd, hr, r(ntot), h
    double precision  k_F, a, d, htm, m
    double precision  pi, veff(ntot), E, E1, vint
    integer i
    
    pi = 4.d0*atan(1.d0)
    E1 = 0.d0
    
    do i = 1,ntot
       veff(i) = d_fr(i)**2
    end do
    
    do i = 1,na+nd
       vint = veff(i)*( 1.d0 - sl(i)*sl(i)/deg )*r(i)*r(i)
       E1 = E1 + vint
    end do
    
    E1 = 2.d0*pi* rho*E1*h*htm/m
    T_2bc= E1
    write(6,*) 'E2bc:::: T2bc', E1
    return
  end function T_2bc
  !
  !=======================================================================
  !
  double precision function V_2bc(na,nd,ntot,h,r,rho,deg,sl,fr2,vr)
    !.....total energy
    !     
    implicit none
    double precision rho, deg
    integer na, nd, nr, ntot
    double precision sl(ntot), fr2(ntot), vr(ntot)
    double precision ha, hd, hr, r(ntot), h
    double precision k_F, a, d, htm, m
    double precision pi, veff(ntot), E2, vint
    integer i
    
    pi = 4.*atan(1.d0)
    E2 = 0.d0
    
    do i = 1, ntot
       E2 = E2+vr(i)*fr2(i)*( 1.d0 - sl(i)*sl(i)/deg )*r(i)*r(i)
    end do
    
    E2 = 2.d0*pi* rho*h*E2
    write(6,*) 'E2bc:::: V2bc', E2
    V_2bc= E2
    return
  end function V_2bc
  !
!=======================================================================
  !
  double precision function  spe_2bc(na,nd,ntot,h,r,deg,m,htm,k,rho, sl,dr_fr)
    !.....single partticle energy
    implicit none
    integer na, nd, ntot
    double precision ha, hd, hr, k_F, a, d, rho, deg, htm, m
    double precision r(ntot), h
    double precision k, sl(ntot), dr_fr(ntot)
    double precision pi, veff(ntot), vint, E 
    integer i      
    
    pi = 4.*atan(1.d0)
    E=0.d0
    
    do i = 1,ntot
       veff(i) = dr_fr(i)**2
    end do
    
    do i=1,na+nd
       vint = veff(i)* ( 1.d0 - sl(i)/deg * sin(k*r(i))/(k*r(i)) ) *r(i)*r(i)
       E = E + vint
    end do
    E = (4.d0*pi* rho*E*h + k*k/2.d0)*htm/m
    
    spe_2bc= E
    return
    
  end function spe_2bc
!
!=======================================================================
!
end MODULE EffectiveInteraction
