!=============================================================================================================
!
! fermionic hard&soft spheres, 
! self-energy( 2nd order) on shell 
! single particle spectrum
! effective mass
! momentum distribution
!
! input: name of the working directory
! (the path has to be ../Data/directory)
! ( create also the directory ../Data/directory/Results)
!
! reads input paramters from unit 5 (file ../Data/directory/param.in)
! kF    : Fermi momentum [fm^-1]
! deg   : degeneracy of momentum eigenstates
! a     : hard-core diameter [fm]
! d     : healing distance [fm]
! rmax  : maximum distance [fm]
! opt   : correlation type(1=fcsin, 2=fceulero)
! htm   :  
! m     : mass
! type  : 1 = HARD spheres, 2 = SOFT spheres
! alpha : V0 (soft spheres) 
!
! reads correlations from ../Data/directory/fcel.in or ../Data/directory/fcsin.in 
!
!
! output in ../Data/directory/Results
!=============================================================================================================
!
program main
  use parameters
  use EffectiveInteraction
  use boundary
  use selfenergy
  
  implicit none
  
  interface

     subroutine im_SE_ek(ImSE,the_err,Ene,kk,dove)      
       double precision ImSE, Ene, kk, the_err
       character*64 dove
     end subroutine im_SE_ek

     subroutine real_SE_ek(ReSE,the_err,Ene,kk,dove)      
       double precision ReSE, Ene, kk, the_err
       character*64 dove
     end subroutine real_SE_ek

     subroutine der5p_err(x,Nx,dx,f,f_er,df,df_er,n0)
       integer Nx, n0
       double precision x(Nx), f(Nx), f_er(Nx), df, df_er
       double precision dx
     end subroutine der5p_err

     subroutine deriv5p(x,Nx,dx,f,df,n0)
       implicit none
       integer Nx, n0
       double precision x(Nx), f(Nx), df
       double precision dx
     end subroutine deriv5p

      subroutine sigma_hf(r, nr, dr, ell, vr, ki, rho, deg, ReSE)
      implicit none
      integer nr, i
      double precision r(nr), vr(nr), ell(nr)
      double precision dr, ki, rho, deg, ReSE
    end subroutine sigma_hf
    
    double precision function  mstar_pert(ki,kF,a,deg)
      implicit none
      double precision ki, kF,a, deg
    end function mstar_pert


    double precision function fint(r,f,nr,hr,x,ichoic)
      implicit real*8 (a-h,o-z)
      dimension R(1:*),F(1:*)  
    end function fint
 end interface



  character*64 arg, arg2,home, dove, filename 
  integer opt, type
  double precision dens
  character temp_str
  integer i, j, k , n, nkp, nkh, i0, k0, j0
  double precision y, E2bc, T2bc, V2bc, TF, Nyq
  double precision constE1, constE2, constk, kcoef, dE1, dE2
  double precision, parameter:: ahs = 1.d0
  double precision, dimension(:), allocatable:: ell, dell, ddell, fcq, vr
  double precision, dimension(:), allocatable:: Ene
  double precision, dimension(:), allocatable:: Reh, Rep, Reh_er, Rep_er
  double precision, dimension(:), allocatable:: Imh, Imp, Imh_er, Imp_er
  double precision, dimension(:), allocatable:: Reh_de, Reh_de_er, Reh_dk, Reh_dk_er
  double precision, dimension(:), allocatable:: Rep_de, Rep_de_er, Rep_dk, Rep_dk_er
  double precision, dimension(:), allocatable:: Rehf, Rehf_dk
  double precision, dimension(:), allocatable:: spe, Retot, Eoff
  double precision, dimension(:), allocatable:: kin, pot
  double precision e0
  integer, parameter :: jmax = 5
  integer, parameter :: stencil = 5
  integer ste_tmp
  double precision  val, dkk, dee, mstarpert
  double precision, dimension(stencil):: Etmp, ktmp
  double precision, dimension(stencil)::  tmp_eh, dtmp_eh, er_eh, der_eh
  double precision, dimension(stencil)::  tmp_ep, dtmp_ep, er_ep, der_ep
  double precision, dimension(stencil)::  tmp_kh, dtmp_kh, er_kh, der_kh
  double precision, dimension(stencil)::  tmp_kp, dtmp_kp, er_kp, der_kp
  double precision, dimension(stencil)::  tmp_hf, dtmp_hf
  character*64 filenum, name, name1, name2
  double precision, dimension(:), allocatable:: Emh,Emh_er,Emp,Emp_er,Em, Em_er
  double precision, dimension(:), allocatable:: kmh, kmh_er,kmp,kmp_er, kmhf, km, km_er
  double precision, dimension(:), allocatable:: mmh, mmh_er,  mmp,mmp_er, mm, mm_er
  double precision  alpha, tmp_dp, spe_hf, zeta, zeta_er



  
  mah = .false.
  if(mah) write(6,*)'------>> MAHAUX matrix elements' 
  
  if(iargc().lt.2) then
     write(*,*) 'name of the directory  with input parameters required'
     write(*,*) ' dk value required'
     stop
  end if
  call getarg(1,arg)
  call getarg(2,arg2)
  
  read(arg2,*) val

  home = '../Data/'
  dove = trim((home))//trim((arg))//'/'
  write(6,*) 'working directory: ', dove
  
  open(5,file=trim(dove)//'param.in') 
  
  ! read system parameters
  !
     
  read(5,6100)kF
  read(5,6100)deg
  read(5,6300)a,d,rmax
  read(5,*)opt
  read(5,6100)htm
  read(5,6100)m
  read(5,*)type,alpha
  
  a = a*ahs
  d = d*ahs
  rmax = rmax* ahs

  hbarc  = 197.3d0
  pi = 4.d0*atan(1.d0)
  dens= deg*kF**3/(6.d0*pi*pi)
  EkF =  kF**2/(2.d0*m)
  write(6,*) 'EkF [fm^-1] = ', EkF, htm*EkF, 'mu [fm^-1]=' , chempot/hbarc
  
  if(opt.eq.1) then
     filename=trim(dove)//"fcsin.in"
     write(6,*) "Correlation ====> sin, deg = ", deg
  else
     filename=trim(dove)//"fcel.in"
     write(6,*) "Correlation ====> Eulero-Lagrange, deg = ", deg
  end if

  if(type.eq.1) write(6,*) 'HARD-SPHERES system'
  if(type.eq.2) write(6,*) 'SOFT-SPHERES system'
  
  !.....read correlations
  open(39,file=filename)
  read(39,*)temp_str,nt,na,nd     
  allocate(r(nt),fc(nt),dfc(nt),ddfc(nt))
  read(39,6400) (r(i),fc(i),dfc(i),ddfc(i), i=1,nt)
  ! define all the extra  parameters (hr, nr) for correlations
  hr = r(2)-r(1)
  nr = nt - nd
  
  allocate(fcq(nt),ell(nt),dell(nt),ddell(nt),veff(nt), vr(nt),kin(nt), pot(nt))
  
  do i = 1, nt
     if(r(i).eq.0.d0) then
        ell(i) = 1.d0
        dell(i) = 0.0d0
        ddell(i) = -kF*kF * 0.2d0
     else
        y        =  kF*r(i)
        ell(i)   =  3./(y**3)*( sin(y)- y*cos(y) )
        dell(i)  =  -(3.*kF/y)*( ell(i)- sin(y)/y)
        ddell(i) =  kF**2 *( 12./(y**2) * (ell(i) - sin(y)/y )- ell(i))
     end if
     fcq(i)   =  fc(i)**2
  end do

  do i = 1, na
     vr(i)  = alpha* (dble(type)-1.d0)
  end do
  do i = na+1, nt
     vr(i)  = 0.d0
  end do
  
  

  open(13,file=trim(adjustl(dove))//'veff_r.dat')
  read(13,*) temp_str
  read(13,6400)( r(i), veff(i), kin(i), pot(i), i = 1, nt)
  close(13)
  
  !     
  !
  TF= 3.d0/5.d0*(kF*kF)/2.*htm/m
  T2bc =  T_2bc(na,nd,nt,hr,r,deg,m,htm,dens,ell,dfc)
  V2bc =  V_2bc(na,nd,nt,hr,r,dens,deg,ell,fcq,vr)
  E2bc =  T2bc + V2bc
  write(6,*) ' '
  write(6,*) 'Energy 2body E2bc = ', TF+E2bc, '   TF =', TF, 'T2bc', E2bc, 'V2bc = ', V2bc, '(T+V) 2bc', E2bc
  write(6,*) ' '
  

  open(2,file=trim(adjustl(dove))//'kgrid.dat')
  read(2,8200) temp_str, nk, nkF
  allocate(kk(nk))
  read(2,8300) temp_str, kmin, kmax, dk
  read(2,6100) ( kk(k), k = 1, nk)
  close(2)
  write(6,*) ' nk,nkF,  dk, kmin, kmax', nk, nkF,  dk, kmin, kmax
  

  open(11,file=trim(adjustl(dove))//'Egridp.dat')
  read(1,8200)  temp_str, nEh, neFh
  read(11,8200) temp_str, nEp, neFp
  allocate(EEh(nEh),dEh(nEh),EEp(nEp),dEp(nEp))
  read(1,8220)  temp_str, Emin, Emax
  read(11,8220) temp_str, Emin, Emax
  read(1,6150) ( EEh(n), dEh(n), n = 1, nEh)      
  read(11,6150) (EEp(n), dEp(n), n = 1, nEp)      
  close(1)
  close(11)
  
  !
  Nyq =  1.d0/(2.d0*hr)
  pcut = 15.d0*kF

  open(47,file=trim(adjustl(dove))//'veff_q.dat')
  read(47,*) temp_str,np
  allocate(p(np),Mp(np),dMp(np),ddMp(np))
  read(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)
  dp   =  p(2)-p(1)


  write(filenum, '(6f5.3)')  kF

  !===========================================
  ! read e(k) = e0(k) + Re(k, e(k))
  !===========================================
  allocate(Ene(nk),spe(nk))

  
  name=trim(adjustl(dove))//"Results/spe_new"//trim(filenum)//".dat"

  open(22,file=name)
  read(22,*) temp_str
  do k = 1, nk
     read(22,*) tmp_dp, e0,spe_hf, spe(k)
     write(6,*) 'spe', k, spe(k)
  end do
  close(22)

     
  
  !===========================================
  ! read  ON SHELL  Sigma2p2h_Im & Sigma2p2h_Re
  !===========================================
  
  allocate(Reh(nk),Rep(nk),Reh_er(nk),Rep_er(nk), Rehf(nk))
  allocate(Imh(nk),Imp(nk),Imh_er(nk),Imp_er(nk))
    
!  name=trim(adjustl(dove))//"Results/self_new"//trim(filenum)//".dat"
!  open(22,file=name)
!  read(22,*) temp_str
!  do k = 1, nk
!     read(22,7000) tmp_dp, Imh(k), Imh_er(k), Reh(k), Reh_er(k), Imp(k), Imp_er(k), Rep(k), Rep_er(k), Rehf(k)
!  end do
!  close(22)
  


!============================================  
! EFFECTIVE MASS 
! DERIVATIVES d/d(omega) &  d/dk at e = spe(k)
!============================================  

  allocate(Reh_de(nk), Reh_de_er(nk), Reh_dk(nk), Reh_dk_er(nk))
  allocate(Rep_de(nk), Rep_de_er(nk), Rep_dk(nk), Rep_dk_er(nk))
  allocate(Rehf_dk(nk))

  !  step dkk & dee
  !dkk = 0.125d0*kF
  dkk = val * kF
  dee = 0.250d0*EkF


  name=trim(adjustl(dove))//"der"//trim(filenum)//".dat"
  open(32,file=name)
  write(32,*) '# dkk = ', dkk, 'dee = ', dee
  write(32,*) '# kk, Reh_de, Reh_de_er, Rep_de, Rep_de_er, Reh_dk, Reh_dk_er, Rep_dk, Rep_dk_er, Rehf_dk'
  flush(32)

  write(6,*) 'effective mass : derivatives in e(k)'
  do k = 1, nk
     Ene(k) = spe(k)
     i0 = stencil/2 + 1
     
     if (kk(k).lt.dkk) then
        k0 = 1
     else if(kk(k).lt.2.d0*dkk) then
        k0 = 2
     else
        k0 = stencil/2 + 1
     end if
     
     do i = 1, stencil
        Etmp(i) = Ene(k) + dble(i-i0)*dee
        ktmp(i) = kk(k)  + dble(i-k0)*dkk
        
        write(6,*) 'derivative E ', i, k, kk(k), Ene(k), kk(k), Etmp(i), ktmp(i)
        part = .true.
        flush(32)
        call real_SE_ek(tmp_ep(i),er_ep(i),Etmp(i),kk(k),dove)  
        flush(32)
        call real_SE_ek(tmp_kp(i),er_kp(i),Ene(k),ktmp(i),dove) 
        part = .false.
        flush(32)
        call real_SE_ek(tmp_eh(i),er_eh(i),Etmp(i),kk(k),dove)
        flush(32)
        call real_SE_ek(tmp_kh(i),er_kh(i),Ene(k),ktmp(i),dove)
       ! flush(32)
        call sigma_hf(r,nt,hr,ell,veff,ktmp(i),dens,deg,tmp_hf(i))
        write(27,*) i, k, Etmp(i), ktmp(i), tmp_eh(i),er_eh(i), tmp_ep(i),er_ep(i),tmp_kh(i),er_kh(i), tmp_kp(i),er_kp(i), tmp_hf(i)
        flush(27)
     end do
     
     call der5p_err(Etmp,stencil,dee,tmp_eh,er_eh,dtmp_eh(i0),der_eh(i0),i0)
     call der5p_err(Etmp,stencil,dee,tmp_ep,er_ep,dtmp_ep(i0),der_ep(i0),i0)
     call der5p_err(ktmp,stencil,dkk,tmp_kh,er_kh,dtmp_kh(k0),der_kh(k0),k0)
     call der5p_err(ktmp,stencil,dkk,tmp_kp,er_kp,dtmp_kp(k0),der_kp(k0),k0)
     call deriv5p(ktmp,stencil,dkk,tmp_hf,dtmp_hf(k0),k0)
     
     
     
     Rehf(k)      =   tmp_hf(k0)
     Rehf_dk(k)   =  dtmp_hf(k0)
     
     Reh(k)    =  tmp_eh(i0)
     Reh_er(k) =   er_eh(i0)
     Rep(k)    =  tmp_ep(i0)
     Rep_er(k) =   er_ep(i0) 
     
     Reh_de(k)    =  dtmp_eh(i0)
     Reh_de_er(k) =   der_eh(i0)
     Rep_de(k)    =  dtmp_ep(i0)
     Rep_de_er(k) =   der_ep(i0)
     
     Reh_dk(k)    =  dtmp_kh(k0)
     Reh_dk_er(k) =   der_kh(k0)
     Rep_dk(k)    =  dtmp_kp(k0)
     Rep_dk_er(k) =   der_kp(k0)
     
     write(32,7000) kk(k), Reh_de(k),Reh_de_er(k), Rep_de(k), Rep_de_er(k), Reh_dk(k),  Reh_dk_er(k), Rep_dk(k), Rep_dk_er(k), Rehf_dk(k)
     flush(32)
  end do
  close(32)

!===========================================
!   effective mass
!===========================================

  allocate(Emh(nk),Emh_er(nk),Emp(nk),Emp_er(nk),Em(nk), Em_er(nk))
  allocate(kmh(nk), kmh_er(nk),kmp(nk),kmp_er(nk), kmhf(nk), km(nk), km_er(nk))
  allocate(mmh(nk), mmh_er(nk),  mmp(nk),mmp_er(nk), mm(nk), mm_er(nk))
  
  name  = 'Results/Emass'//trim(filenum)//'.dat'
  name1 = 'Results/kmass'//trim(filenum)//'.dat'
  name2 = 'Results/Mass'//trim(filenum)//'.dat'
!
  open(1,file=trim(adjustl(dove))//name)
  open(2,file=trim(adjustl(dove))//name1)
  open(3,file=trim(adjustl(dove))//name2)
  write(1,*) '# kk, kk/kF, Emh, Emh_er, Emp, Emp_er, Em, Em_er'
  write(2,*) '# kk, kk/kF, kmh, kmh_er, kmp, kmp_er, km, km_er, kmhf'
  write(3,*) '# kk, kk/kF, mmh, mmh_er, mmp, mmp_er, mm, mm_er,  mstarpert'
     
  do k = 1, nk
     Emh(k)    =  1.d0 - Reh_de(k)
     Emh_er(k) = Reh_de_er(k)
     Emp(k)    =  1.d0 - Rep_de(k)
     Emp_er(k) = Rep_de_er(k)
     Em(k)     =  Emh(k) + Emp(k) - 1.d0
     Em_er(k)  = sqrt(Emh_er(k)**2 + Emp_er(k)**2 )

     kmh(k)    =  1.0d0/(1.d0 + Reh_dk(k)* m/kk(k) )
     kmh_er(k) = kmh(k)**2 * Reh_dk_er(k)* m/kk(k)
     kmp(k)    =  1.0d0/(1.d0 + Rep_dk(k)* m/kk(k) )
     kmp_er(k) = kmh(k)**2 * Rep_dk_er(k)* m/kk(k)
     kmhf(k)   =  1.0d0/(1.d0 + Rehf_dk(k) *m/kk(k) )     
     km(k)     = (Reh_dk(k) + Rep_dk(k) + Rehf_dk(k) ) * m /kk(k)
     km(k)     =  1.0d0/( 1.d0 + km(k) )
     km_er(k)  = km(k)**2* m/kk(k)* sqrt(Reh_dk_er(k)**2 + Rep_dk_er(k)**2)

     mmh(k)    =  Emh(k)*kmh(k)
     mmh_er(k) = sqrt(Emh_er(k)**2*kmh(k)**2+Emh(k)**2*kmh_er(k)**2 )
     mmp(k)    =  Emp(k)*kmp(k)
     mmp_er(k) = sqrt(Emp_er(k)**2*kmp(k)**2+Emp(k)**2*kmp_er(k)**2 )
     mm(k)     =  Em(k)*km(k)
     mm_er(k)  = sqrt(Em_er(k)**2*km(k)**2+Em(k)**2*km_er(k)**2 )
     
     mstarpert = mstar_pert(kk(k),kF,ahs,deg)

     write(1,6800) kk(k), kk(k)/kF, Emh(k), Emh_er(k), Emp(k), Emp_er(k), Em(k), Em_er(k)
     write(2,6900) kk(k), kk(k)/kF, kmh(k), kmh_er(k), kmp(k), kmp_er(k), km(k), km_er(k), kmhf(k)
     write(3,6900) kk(k), kk(k)/kF, mmh(k), mmh_er(k), mmp(k), mmp_er(k), mm(k), mm_er(k),  mstarpert   
  end do

  close(1)
  close(2)

  
  mstarpert = 8.d0*(deg-1.d0)/(15.d0*pi**2)* (7.d0*log(2.d0) -1.d0)*(kF*ahs)**2
  mstarpert = 1.d0 + mstarpert

  write(3,*) '# Perturbative value m*(kF)/m = ', mstarpert
  write(3,*) '# Effective interaction  m*(kF)/m = ', fint(kk,mm,nk,dk,kF,1)

  close(3)

  write(6,*) 'Perturbative value m*(kF)/m = ', mstarpert
  write(6,*) 'Effective interaction  m*(kF)/m = ', fint(kk,mm,nk,dk,kF,1)



!===========================================
!   zeta
!===========================================
  name = trim(adjustl(dove))//"Results/zeta_spe"//trim(filenum)//".dat"
  open(2,file=name)
  write(2,*) '# k z(k) z_er(k)  (nk(kF-) - n(kF+))  k/kF'
  do k = 1, nk
     zeta =  Reh_de(k)  + Rep_de(k)
     zeta = 1.d0 / (1.d0 - zeta)
     zeta_er = zeta**2 * sqrt(Reh_de_er(k)**2 + Rep_de_er(k)**2)
     write(2,6400) kk(k), zeta, zeta_er, kk(k)/kF
  end do


! deallocate
  deallocate(r,fc,dfc,ddfc)
  deallocate(fcq,ell,dell,ddell,veff)
  deallocate(p,Mp,dMp,ddMp)
  deallocate(kk,EEh,dEh,EEp,dEp)
  deallocate(Ene)
  deallocate(Reh,Rep,Reh_er,Rep_er)
  deallocate(Imh,Imp,Imh_er,Imp_er)
  deallocate(Reh_de, Reh_de_er, Reh_dk, Reh_dk_er)
  deallocate(Rep_de, Rep_de_er, Rep_dk, Rep_dk_er)
  deallocate(Rehf, Rehf_dk)
  deallocate(Emh,Emh_er,Emp,Emp_er,Em, Em_er)
  deallocate(kmh, kmh_er,kmp,kmp_er, kmhf, km, km_er)
  deallocate(mmh, mmh_er,mmp,mmp_er, mm, mm_er)
  deallocate(spe)
  deallocate(kin, pot)

  
2222 format(e25.15)
6100 format(e12.4)
6150 format(2e12.4)
6200 format(2e20.10e3)
6300 format(3e12.4)
6400 format(4e20.10)    
6500 format(5e17.8)
6700 format(7e17.8)
7100 format(11e17.8)
8200 format(a,2i5.1)
8220 format(a,2e12.4)
8300 format(a,3e12.4)      
6800 format(8e17.8)
6900 format(9e17.8)
7000 format(10e17.8)
7600 format(6e17.8)
  
  
  
  
  stop
end program main
