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

    ! subroutine im_SE_ek(ImSE,the_err,Ene,kk,dove)      
    !   double precision ImSE, Ene, kk, the_err
    !   character*64 dove
    ! end subroutine im_SE_ek

     !subroutine real_SE_ek(ReSE,the_err,Ene,kk,dove)      
      ! double precision ReSE, Ene, kk, the_err
       !character*64 dove
     !end subroutine real_SE_ek

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

     ! subroutine sigma_hf(r, nr, dr, ell, vr, ki, rho, deg, ReSE)
     ! implicit none
     ! integer nr, i
     ! double precision r(nr), vr(nr), ell(nr)
     ! double precision dr, ki, rho, deg, ReSE
    !end subroutine sigma_hf
    
    double precision function  mstar_pert(ki,kF,a,deg)
      implicit none
      double precision ki, kF,a, deg
    end function mstar_pert

    subroutine nk_pert(ki,nk,kF,a,deg,rho1,rho2)
      implicit none
      integer nk
      double precision ki(nk), rho1(nk),rho2(nk)
      double precision deg,a,kF
    end subroutine

    double precision function fint(r,f,nr,hr,x,ichoic)
      implicit real*8 (a-h,o-z)
      dimension R(1:*),F(1:*)  
    end function fint
 end interface



  character*128 arg, home, dove, filename, varname, workdir
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
  double precision e0, spein, ein, eout, ris
  integer, parameter :: jmax = 5
  integer, parameter :: stencil = 5
  integer ste_tmp
  double precision  dkk, dee, mstarpert
  double precision, dimension(stencil):: Etmp, ktmp
  double precision, dimension(stencil)::  tmp_eh, dtmp_eh, er_eh, der_eh
  double precision, dimension(stencil)::  tmp_ep, dtmp_ep, er_ep, der_ep
  double precision, dimension(stencil)::  tmp_kh, dtmp_kh, er_kh, der_kh
  double precision, dimension(stencil)::  tmp_kp, dtmp_kp, er_kp, der_kp
  double precision, dimension(stencil)::  tmp_hf, dtmp_hf
  character*128 filenum, name, name1, name2
  double precision, dimension(:), allocatable:: Emh,Emh_er,Emp,Emp_er,Em, Em_er
  double precision, dimension(:), allocatable:: kmh, kmh_er,kmp,kmp_er, kmhf, km, km_er
  double precision, dimension(:), allocatable:: mmh, mmh_er,  mmp,mmp_er, mm, mm_er
  double precision, dimension(:), allocatable:: nmag, nmin, rho1, rho2
  double precision zeta, zeta_er, zeta1, alpha
  double precision Ehf, Uk, Ek
  double precision s0, s, s1, ss, ss1, sss, sss1, delta , delta1 , norm
  double precision t0, t, t1, tt, tt1, ttt, ttt1, deltat, deltat1, normt
  double precision thf, tp, th, normself, ekin, tderp, tderh
  integer  maxit
  double precision small, realh, realh_er, realp, realp_er, dE
  double precision asoft
  integer status_val
  
  mah = .false.
  if(mah) write(6,*)'------>> MAHAUX matrix elements' 
  
  if(iargc().lt.1) then
     write(6,*) 'name of the directory  with input parameters required'
     stop
  end if
  call getarg(1,arg)
  
 CALL GET_ENVIRONMENT_VARIABLE("WORK", varname, STATUS= status_val)

  if (status_val.eq.2) then
  write(*,*) 'WARNING: the processor does not support environment variables '
  else if (status_val.eq.-1) then
     write(*,*) 'WARNING: varname is too short: varname is truncated'
  else if(status_val .eq.1) then
     write(*,*) 'WARNING: $WORK environment variable does not exist'
  end if

  if(status_val.eq.0) then
     workdir = trim(varname)//'/HS_selfenergy/'
  else
     workdir = '../'
  end if

  home = trim(workdir)//'Data/'

  write(*,*) 'FINAL RESULTS IN : ',  home


!  home = '../Data/'
  dove = trim((home))//trim((arg))//'/'
  write(6,*) 'working directory: ', dove
  maxit = 200
  small = 1.e-9
  
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
  open(39,file=filename)

  if(type.eq.1) write(6,*) 'HARD-SPHERES system'
  if(type.eq.2) write(6,*) 'SOFT-SPHERES system'
  
  !.....read correlations
  read(39,*)temp_str,nt,na,nd     
  allocate(r(nt),fc(nt),dfc(nt),ddfc(nt))
  do i = 1, nt
     read(39,6400) r(i),fc(i),dfc(i),ddfc(i)
  end do
  ! define all the extra  parameters (hr, nr) for correlations
  hr = r(2)-r(1)
  nr = nt - nd
  
  allocate(fcq(nt),ell(nt),dell(nt),ddell(nt),veff(nt), vr(nt))
  
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

    if(type.eq.2 ) then
       asoft = 0.d0
       do i = 1, nt
          asoft = asoft + vr(i)*r(i)*r(i)
       end do
       
       asoft = asoft * hr * m**2/htm
       write(6,*) 'a  soft =',  asoft
    end if
    
  !.....compute the effective interaction [ (hbar**2/m)  ]
  
    allocate(kin(i), pot(i))
    do i = 1,nt
       !     veff(i) = htm/m * dfc(i)**2 + fcq(i)*vr(i)
       kin(i) = htm/m * dfc(i)**2
       pot(i) = fcq(i)*vr(i)
       veff(i) = kin(i) + pot(i)
  end do
  !     
  open(13,file=trim(adjustl(dove))//'veff_r.dat')
  write(13,*) ' #', nt,kF
  do i = 1, nt
     write(13,6400) r(i), veff(i), kin(i), pot(i)
end do
  close(13)
  !
  TF= 3.d0/5.d0*(kF*kF)/2.*htm/m
  T2bc =  T_2bc(na,nd,nt,hr,r,deg,m,htm,dens,ell,dfc)
  V2bc =  V_2bc(na,nd,nt,hr,r,dens,deg,ell,fcq,vr)
  E2bc =  T2bc + V2bc
  write(6,*) ' '
  write(6,*) 'Energy 2body E2bc = ', TF+E2bc, '   TF =', TF, 'T2bc', E2bc, 'V2bc = ', V2bc, '(T+V) 2bc', E2bc
  write(6,*) ' '
  
  
  constE1 = 0.1d0
  constE2 = 0.5d0
!  constk  = 0.075d0
  constk  = 0.15d0
  dk  = constk*kF
  dE1 = constE1*EkF
  dE2 = constE2*EkF
  kcoef = 5.0d0
  kmax =  1.5d0*kF
  kmin =  0.50d0*kF
  Emax =  kcoef*kcoef*EkF
  Emin = -kcoef*kcoef*EkF
  
  !nk = int((kmax - kmin)/dk + 0.5d0) + 1 
  nk = 1
  !nEh = int(EkF/dE1 + 0.5d0)
  !nEp = int(EkF/dE2 + 0.5d0)
  nEh = 10
  nEp = 10
  write(6,*) 'nEh, nEp, dEh, dEp, Emin, Emax', nEh, nEp, dE1, dE2, Emax, Emin
  
  ! allocate grid E & k 
  allocate(kk(nk),EEh(nEh),dEh(nEh),EEp(nEp),dEp(nEp))
  write(6,*) ' nk, dk, kmin, kmax', nk, dk, kmin, kmax
  

  if(nk.eq.1) then
     kk(nk) = 1.0d0*kF
     nkF = nk
     write(6,*) '------>  nk = 1' , kk(nk), kF
     !     
  else if (kF.gt.kmin) then
     nkh = int((kF-kmin)/dk) + 1
     nkp = nk - nkh
     nkF = nkh
     if(kF.gt.kmax) then
        do k = 1, nk
           kk(k) = kmin + (dble(k))*dk
           write(6,*)  'kk',k,kk(k), dk
        end do
     else
        do k = 1, nkh
           kk(k) = kF - dble(nkh-k)*dk
           write(6,*) 'kk', k,kk(k)
        end do
        do k = nkh + 1 , nk
           kk(k) = kF + dble(k-nkh)*dk
           write(6,*) k,kk(k)
        end do
     end if
  else
     nkh = 0
     nkp = nk
     do k = 1, nkp
        kk(k) = kmin + dble(k)*dk
        write(6,*) k,kk(k)
     end do
  end if
  !
  do n = 1, nEh
     dEh(n) = dE1
  end do
  do n = 1, nEp
     dEp(n) = dE2
  end do
  
  if (Emin.gt.EkF) then
     write(6,*) '-----> Emin > EkF'
     do n = 1, nEh
        EEh(n)  = Emin + (dble(n)) * dEh(n)
     end do
     do n = 1, nEp
        EEp(n)  = Emin + (dble(n)) * dEp(n)
     end do
     
  else     
     neFh = int((EkF  - Emin)/dE1 + 0.5d0)
     neFp = int((EkF  - Emin)/dE2 + 0.5d0)
     do n =  1, nEh
        EEh(n) = EkF + (dble(n-neFh)) *dEh(n)
     end do
     do n =  1, nEp
        EEp(n) = EkF + (dble(n-neFp)) *dEp(n)
     end do
  end if


  open(11,file=trim(adjustl(dove))//'Egridp.dat')
  write(1,8200)  '#', nEh, neFh
  write(11,8200) '#', nEp, neFp
  write(1,8220)  '#', Emin, Emax
  write(11,8220) '#', Emin, Emax
  write(1,6150) ( EEh(n), dEh(n), n = 1, nEh)      
  write(11,6150) (EEp(n), dEp(n), n = 1, nEp)      
  close(1)
  close(11)
  
  open(2,file=trim(adjustl(dove))//'kgrid.dat')
  write(2,8200) '#', nk, nkF
  write(2,8300) '#', kmin, kmax, dk
  write(2,6100) ( kk(k), k = 1, nk)
  close(2)
  !
  ! calcola la trasf do fourier del potenziale
  Nyq =  1.d0/(2.d0*hr)
  pcut = 15.d0*kF
  np = 5000
  allocate(p(np),Mp(np),dMp(np),ddMp(np))
  dp   =  Pcut/dble(np)
  do  i = 1, np
     p(i) = (dble(i)-0.5d0)*dp
  end do
  !
  call Fourier(r,veff,nt,hr,p,Mp,np,dp)
  call der(p,1,np,dp,Mp,dMp,ddMp)
  !
  open(47,file=trim(adjustl(dove))//'veff_q.dat')
  write(47,*)'#',np,kF
  write(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)
  
  write(filenum, '(6f5.3)')  kF



  !===========================================
  ! solve e(k) = e0(k) + Re(k, e(k))
  !===========================================
  allocate(Ene(nk),spe(k))
  allocate(Reh(nk),Rep(nk),Reh_er(nk),Rep_er(nk), Rehf(nk))
  allocate(Imh(nk),Imp(nk),Imh_er(nk),Imp_er(nk))
  allocate(Retot(jmax), Eoff(jmax))

  
  name=trim(adjustl(dove))//"Results/spe_new"//trim(filenum)//".dat"
  name1=trim(adjustl(dove))//"Results/converg"//trim(filenum)//".dat"

  open(22,file=name)
  write(22,*) '#  k, e0, spe_old, spe_new, diff, diff %'
  open(23, file=name1)



  j0 = int((jmax+1)/2.d0)
  write(6,*) 'single particle energy e(k) = e0(k) + Re(k, e(k))'
  do k = 1, nk
     Ene(k) = (kk(k) )**2/(2.d0*M)
     call sigma_hf(r,nt,hr,ell,veff,kk(k),dens,deg,Rehf(k))
     dE = 2.d0*Rehf(k)/dble(jmax)
     write(6,*) 'k, ReHF', k, Rehf(k)
     do j = 1, jmax
        Eoff(j) = Ene(k) + Rehf(k) + dble(j-j0)*dE
        write(6,*) k, j, Eoff(j)
        part = .false.
        flush(23)
        flush(22)
        call real_SE_ek(realh,realh_er,Eoff(j),kk(k),dove)      
        part = .true.
        flush(23)
        flush(22)
        call real_SE_ek(realp,realp_er,Eoff(j),kk(k),dove)
        Retot(j) = realp + realh + Rehf(k)
     end do
     e0 = Ene(k)
     spein = e0 + Retot(j0)    
     ein = e0
     eout = e0
     !write(6,*) 'start', ein, eout
     do i = 1, maxit
        ein  =  e0 +  fint(Eoff,Retot,jmax,dE,ein,1)
        eout =  e0 +  fint(Eoff,Retot,jmax,dE,ein,1)
        !write(6,*) 'in,out', ein, eout
        if(abs(ein-eout).lt.small) then
           ris = eout
           write(23,*), 'convergence!', k, i , abs(ein-eout), eout
           exit
        else
           ein = eout
        end if
     end do
     spe(k) = ris
     write(22,7600) kk(k), e0, spein, ris,(ris-spein),(ris-spein)/spein
  end do
  close(22)
  close(23)

     



!===========================================
! momentum distribution
! DERIVATIVES d/d(omega) at  e = e_0(k)
!============================================  

  allocate(Reh_de(nk), Reh_de_er(nk), Reh_dk(nk), Reh_dk_er(nk))
  allocate(Rep_de(nk), Rep_de_er(nk), Rep_dk(nk), Rep_dk_er(nk))
  allocate(Rehf_dk(nk))

!  step dkk & dee
  dkk = 0.125d0*kF
  dee = 0.250d0*EkF

  name=trim(adjustl(dove))//"der_e0"//trim(filenum)//".dat"
  open(32,file=name)
  write(32,*) '# dkk = ', dkk, 'dee = ', dee
  write(32,*) '# kk, Reh_de, Reh_de_er, Rep_de, Rep_de_er'
  flush(32)

  write(6,*) 'momentum distribution: d Re/dE in e0(k)'
  do k = 1, nk
     Ene(k) = kk(k)**2/(2.d0*m)
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
        write(6,*) 'derivative in e0 ', i, k, kk(k), Ene(k),  Etmp(i)        
        part = .true.
        flush(32)
        call real_SE_ek(tmp_ep(i),er_ep(i),Etmp(i),kk(k),dove)  
        part = .false.
        flush(32)
        call real_SE_ek(tmp_eh(i),er_eh(i),Etmp(i),kk(k),dove)
        flush(32)
        write(27,*) i, k, Etmp(i), tmp_eh(i),er_eh(i), tmp_ep(i),er_ep(i)
        flush(27)
     end do
     
     call der5p_err(Etmp,stencil,dee,tmp_eh,er_eh,dtmp_eh(i0),der_eh(i0),i0)
     call der5p_err(Etmp,stencil,dee,tmp_ep,er_ep,dtmp_ep(i0),der_ep(i0),i0)

     
     Reh(k)    =  tmp_eh(i0)
     Reh_er(k) =   er_eh(i0)
     Rep(k)    =  tmp_ep(i0)
     Rep_er(k) =   er_ep(i0) 
     
     Reh_de(k)    =  dtmp_eh(i0)
     Reh_de_er(k) =   der_eh(i0)
     Rep_de(k)    =  dtmp_ep(i0)
     Rep_de_er(k) =   der_ep(i0)
          
     write(32,7000) kk(k), Reh_de(k),Reh_de_er(k), Rep_de(k), Rep_de_er(k)
     flush(32)
  end do
  close(32)


!===========================================
!    momentum distribution
!===========================================




  allocate(nmag(nk), nmin(nk), rho1(nk), rho2(nk))
  
  
  !.....pertubative n(k)
  call nk_pert(kk,nk,kF,ahs,deg,rho1,rho2)
  
  !....effective interaction n(k)
  do k = 1, nk
     if(kk(k).le.kF) then
        nmin(k) = 1.d0 + Rep_de(k)
     else 
        nmin(k) = 0.d0
     end if
     
     if (kk(k).ge.kF) then
        nmag(k) =   - Reh_de(k)
     else
        nmag(k) = 0.d0
     end if
  end do
  
  if(nk.gt.(2*nkF-1)) then
     s    = 0.d0
     s1   = 0.d0
     ss   = 0.d0
     ss1  = 0.d0
     sss  = 0.d0
     sss1 = 0.d0
     
     t    = 0.d0
     t1   = 0.d0
     tt   = 0.d0
     tt1  = 0.d0
     ttt  = 0.d0
     ttt1 = 0.d0
     thf  = 0.d0
     t0   = 0.d0
     th   = 0.d0
     ekin = 0.d0


     do k = 1, nkF-1
        s = s + (kk(k)**2*nmin(k)+kk(k+1)**2*nmin(k+1))/2.d0    * (kk(k+1) - kk(k))
        s1= s1+ (kk(k)**2*rho1(k) +kk(k+1)**2* rho1(k+1))/2.d0  * (kk(k+1) - kk(k)) 
        t = t + (kk(k)**4*nmin(k)+kk(k+1)**4*nmin(k+1))/2.d0    * (kk(k+1) - kk(k))
        t1= t1 +(kk(k)**4*rho1(k) + kk(k+1)**4* rho1(k+1))/2.d0 * (kk(k+1) - kk(k))
        thf = thf +( kk(k)**2 * Rehf(k) +  kk(k+1)**2 * Rehf(k+1)) / 2.d0 * (kk(k+1) - kk(k)) 
        tp =  tp + ( kk(k)**2 * Rep(k) + kk(k+1)**2 * Rep(k+1)) / 2.d0 *  (kk(k+1) - kk(k))
        ekin = ekin + (kk(k)**4 + kk(k+1)**4)/2.d0 *(kk(k+1) - kk(k))
     end do
     
     do k = nkF, 2*nkF-1 
        ss = ss + (kk(k)**2*nmag(k)+kk(k+1)**2*nmag(k+1))/2.d0 * (kk(k+1) - kk(k))
        ss1= ss1+ (kk(k)**2*rho2(k) +kk(k+1)**2* rho2(k+1))/2.d0 * (kk(k+1) - kk(k))
        tt = tt + (kk(k)**4*nmag(k)+kk(k+1)**4*nmag(k+1))/2.d0 * (kk(k+1) - kk(k))
        tt1= tt1 +(kk(k)**4*rho2(k) + kk(k+1)**4* rho2(k+1))/2.d0 * (kk(k+1) - kk(k))
        th =  th + ( kk(k)**2 * Reh(k) + kk(k+1)**2 * Reh(k+1)) / 2.d0 *  (kk(k+1) - kk(k))

     end do
     
     do k = 2*nkF, nk-1
        sss = sss + (kk(k)**2*nmag(k) + kk(k+1)**2*nmag(k+1))/2.d0 * (kk(k+1) - kk(k))
        sss1= sss1+ (kk(k)**2*rho2(k) + kk(k+1)**2* rho2(k+1))/2.d0 * (kk(k+1) - kk(k))
        ttt = ttt + (kk(k)**4*nmag(k) + kk(k+1)**4*nmag(k+1))/2.d0 * (kk(k+1) - kk(k))
        ttt1= ttt1 +(kk(k)**4*rho2(k) + kk(k+1)**4* rho2(k+1))/2.d0 * (kk(k+1) - kk(k))
        th =  th + ( kk(k)**2 * Reh(k) + kk(k+1)**2 * Reh(k+1)) / 2.d0 *  (kk(k+1) - kk(k))
     end do
     
     s0 = 4.d0*pi
     t0 = 4.d0*pi/(2.d0*pi)**3
     s  = s*s0      
     s1 = s1*s0     
     ss = ss*s0     
     ss1 = ss1*s0   
     sss = sss*s0   
     sss1 = sss1*s0 
     
     t    = t    * t0 
     t1   = t1   * t0 
     tt   = tt   * t0  
     tt1  = tt1  * t0 
     ttt  = ttt  * t0 
     ttt1 = ttt1 * t0
     thf = thf * t0
     tp = tp * t0
     th = th * t0
     ekin = ekin * t0
     
     tderp = t - ekin
     tderh = -(tt + ttt)

     norm = 4.d0*pi*kF**3/3.d0
     normt= deg/(2.d0*m* dens)
     normself= deg/(2.d0* dens)
     
     delta  = (kk(1))**2*nmin(1)*dk*4.d0*pi
     delta1 = (kk(1))**2*rho1(1)*dk*4.d0*pi
     deltat = (kk(1))**4*nmin(1)*dk*4.d0*pi
     deltat1= (kk(1))**4*rho1(1)*dk*4.d0*pi
     
     
     name = trim(adjustl(dove))//"Results/nk_norm"//trim(filenum)//".dat"

     open(16,file=name)
     write(16,*) ' '
     write(16,*)'------- n(k) Normalization-----' 
     write(16,*) 'n(k) eff  =',(s+ss+sss)/norm,    'kmax/kF =', kk(nk)/kF
     write(16,*) 'n(k) pert =',(s1+ss1+sss1)/norm, 'kmax/kF =',kk(nk)/kF
     write(16,*)' '
     write(16,*) 'k<kF         : n(k) eff =' ,s/norm, 'n(k) pert =',s1/norm
     write(16,*) 'kF < k < 2kF : n(k) eff =',ss/norm, 'n(k) pert =', ss1/norm
     write(16,*) '2kF < k < 5kF: n(k) eff =',sss/norm,'n(k) pert =',sss1/norm
     write(16,*)' '
     write(16,*) ' k < k(1) : delta eff = '  , delta/norm,'delta pert = ',delta1/norm
     
     write(16,*) ' '
     write(16,*)'------- n(k) Kinetic Energy-----' 
     write(16,*) 'T_2bc       = ', TF+E2bc, '   TF =', TF, 'E2bc=', E2bc
     write(16,*) 'T eff       =',    (t+tt+ttt)*normt, 'kmax/kF =', kk(nk)/kF
     write(16,*) 'T pert      =', (t1+tt1+ttt1)*normt, 'kmax/kF =', kk(nk)/kF
     write(16,*) 'T eff self  =', (thf + tp-th)*normself, 'kmax/kF =', kk(nk)/kF
     write(16,*) 'T eff total =', (t+tt+ttt)*normt + (thf+tp-th)*normself, 'kmax/kF =', kk(nk)/kF
     write(16,*)' '
     write(16,*) 'k < kF       : T eff=',t*normt,  'T pert =', t1*normt,   'diff % =', (t-t1)/t1
     write(16,*) 'kF< k < 2kF  : T eff=',tt*normt, 'T pert =', tt1*normt,  'diff % =',  (tt-tt1)/tt1
     write(16,*) '2kf < k < 5kF: T eff=',ttt*normt,'T pert =', ttt1*normt, 'diff % =', (ttt-ttt1)/ttt1
     write(16,*) 'k < kF       : E total eff=',t*normt + thf*normself
     write(16,*) ' Self HF = ', thf*normself
     write(16,*) ' Self p  = ', tp*normself, 'Self  h =',  th*normself, ' Self 2order', (tp-th)*normself
     write(16,*) ' '
     write(16,*) ' E fg     = ', ekin*normt
     write(16,*) ' E HF      = ', thf*normself
     write(16,*) ' E rep   = ', tp*normself,  ' E reh  = ', th*normself, ' E rep+reh  = ', (tp-th)*normself
     write(16,*) ' E derp  = ', tderp*normt,  ' E derh  = ', tderh*normt, ' E derp+derh  = ', (tderp-tderh)*normt
     write(16,*) ' E total 1ord    = ', ekin*normt + thf*normself
     write(16,*) ' E total 2nd ord = ', ekin*normt + thf*normself + (tp-th)*normself + (tderp-tderh)*normt, &
          'delta E = ', (tp-th)*normself + (tderp-tderh)*normt


     write(16,*) ' '
     write(16,*) 'k < k(1) :  delta eff = ', deltat*normt, 'delta  pert = ',deltat1*normt
     write(16,*) ' '
     write(16,*) ' '  
     write(16,*) ' (Teff-E2bc)   % ',   (t+tt+ttt)*normt/ (TF+E2bc ) -  1.d0
     write(16,*) ' (Tpert-E2bc)  % ',  (t1+tt1+ttt1)*normt/ (TF+E2bc ) -  1.d0
     write(16,*) ' (T+self-E2bc) % ',   ((t+tt+ttt)*normt + (thf+tp-th)*normself) /(TF+E2bc ) -  1.d0


     close(16)

     
  
  end if
  !write results n(k)

  name = trim(adjustl(dove))//"Results/nk"//trim(filenum)//".dat"
  open(2,file=name)
  write(2,8220) '# kF, kmax :', kF, kk(nk)
  write(2,8220) '# effective    nk_norm, kinetic energy :',(s+ss+sss)/norm,(t+tt+ttt)*normt
  write(2,8220) '# perturbative nk_norm, kinetic energy :',(s1+ss1+sss1)/norm,(t1+tt1+ttt1)*normt
  write(2,*) '# nk_eff nk_pert k kF k/kF'
  do k = 1, nk
     if(kk(k).le.kF) then
        write(2,6500) kk(k), nmin(k), Rep_de_er(k), rho1(k), kk(k)/kF
     end if
     if(kk(k).ge.kF) then
        write(2,6500) kk(k), nmag(k), Reh_de_er(k), rho2(k),  kk(k)/KF
     end if
  end do
  close(2)

     
     


!============================================  
! EFFECTIVE MASS 
! DERIVATIVES d/d(omega) &  d/dk at e = spe(k)
!============================================  


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
        flush(32)
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
     
     write(32,7000) kk(k), Reh_de(k),Reh_de_er(k), Rep_de(k), Rep_de_er(k), Reh_dk(k), &
          Reh_dk_er(k), Rep_dk(k), Rep_dk_er(k), Rehf_dk
     flush(32)
  end do
  close(32)
  
  !===========================================
  !  ON SHELL  Sigma2p2h_Im & Sigma2p2h_Re
  !===========================================
  
    
  name=trim(adjustl(dove))//"Results/self_new"//trim(filenum)//".dat"
  open(22,file=name)
  write(22,*) '# kk, Imh, Imh_er, Reh,Reh_er, Imp, Imp_er, Rep, Rep_er, Rehf '
  flush(22)
     

  do k = 1, nk
     part = .false.
     call real_SE_ek(Reh(k),Reh_er(k),spe(k),kk(k),dove)      
     call im_SE_ek(Imh(k),Imh_er(k),spe(k),kk(k),dove)      
     part = .true.
     call real_SE_ek(Rep(k),Rep_er(k),spe(k),kk(k),dove)      
     call im_SE_ek(Imp(k),Imp_er(k),spe(k),kk(k),dove)      
     write(22,7000) kk(k), Imh(k), Imh_er(k), Reh(k), Reh_er(k), Imp(k), Imp_er(k), Rep(k), Rep_er(k), Rehf(k)
     flush(22)
  end do

  close(22)

!===========================================
!   Optical potential
!===========================================

  name=trim(adjustl(dove))//"Results/Optical"//trim(filenum)//".dat"
  open(21,file=name)
  write(21,*) '# k E E_hf Uk Ek Reh Rep Rehf'
  do k = 1, nk
     e0 = kk(k)**2/(2.d0*m)
     EhF  = e0 + Rehf(k)
     Ek   = e0 + Rehf(k) + Rep(k) + Reh(k)
     Uk   = Rehf(k) + Rep(k) + Reh(k)
     write(21,6800) kk(k), e0, EhF, Uk, Ek , Reh(k), Rep(k), Rehf     
  end do
  close(21)


  


  
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
     Emh(k) =  1.d0 - Reh_de(k)
     Emh_er(k) = Reh_de_er(k)
     Emp(k) =  1.d0 - Rep_de(k)
     Emp_er(k) = Rep_de_er(k)
     Em(k)  =  Emh(k) + Emp(k) - 1.d0
     Em_er(k) = sqrt(Emh_er(k)**2 + Emp_er(k)**2 )

     kmh(k) =  1.0d0/(1.d0 + Reh_dk(k)* m/kk(k) )
     kmh_er(k) = kmh(k)**2 * Reh_dk_er(k)* m/kk(k)
     kmp(k) =  1.0d0/(1.d0 + Rep_dk(k)* m/kk(k) )
     kmp_er(k) = kmh(k)**2 * Rep_dk_er(k)* m/kk(k)
     kmhf(k)=  1.0d0/(1.d0 + Rehf_dk(k) *m/kk(k) )     
     km(k) = (Reh_dk(k) + Rep_dk(k) + Rehf_dk(k) ) * m /kk(k)
     km(k) =  1.0d0/(1.d0 + km(k))
     km_er(k) = km(k)**2* m/kk(k)* sqrt(Reh_dk_er(k)**2 + Rep_dk_er(k)**2)

     mmh(k) =  Emh(k)*kmh(k)
     mmh_er(k) =sqrt(Emh_er(k)**2*kmh(k)**2+Emh(k)**2*kmh_er(k)**2 )
     mmp(k) =  Emp(k)*kmp(k)
     mmp_er(k) =sqrt(Emp_er(k)**2*kmp(k)**2+Emp(k)**2*kmp_er(k)**2 )
     mm(k)  =  Em(k)*km(k)
     mm_er(k) =sqrt(Em_er(k)**2*km(k)**2+Em(k)**2*km_er(k)**2 )
     
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
  name = trim(adjustl(dove))//"Results/zeta"//trim(filenum)//".dat"
  open(2,file=name)
  write(2,*) '# k z(k) z_er(k)  (nk(kF-) - n(kF+))  k/kF'
  do k = 1, nk
     zeta =  Reh_de(k)  + Rep_de(k)
     zeta = 1.d0 / (1.d0 - zeta)
     zeta1 = nmin(k) - nmag(k)
     zeta_er = zeta**2 * sqrt(Reh_de_er(k)**2 + Rep_de_er(k)**2)
     write(2,*) kk(k), zeta, zeta_er, zeta1, kk(k)/kF
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
  deallocate(nmag,nmin,rho1,rho2)
  deallocate(Emh,Emh_er,Emp,Emp_er,Em, Em_er)
  deallocate(kmh, kmh_er,kmp,kmp_er, kmhf, km, km_er)
  deallocate(mmh, mmh_er,mmp,mmp_er, mm, mm_er)
  deallocate(spe,Retot, Eoff)
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
