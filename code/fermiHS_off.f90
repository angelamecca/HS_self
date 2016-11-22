!=============================================================================================================
!
! fermionic hard&soft spheres, 
! self-energy( 2nd order) off shell
! spectral functions 
!
! input: name of the working directory
! (the path has to be ../Data/directory)
! ( create also the directory ../Data/directory/Spectral/ for output)
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
! output in ../Data/directory/ and ../Data/directory/Spectral/
!=============================================================================================================
!

program main
  use parameters
  use EffectiveInteraction
  use boundary
  use selfenergy
  
  implicit none
  interface
     !subroutine imaginary_SE(ImSE,ImSE_er,EE,nE,dE,kk,nk,dk,dove)   
     !  integer nE, nk
     !  double precision dk
     !  double precision ImSE(nE,nk), ImSE_er(nE,nk),  EE(nE), dE(nE), kk(nk)
     !  character*64 dove
     !end subroutine imaginary_SE
     
     !subroutine real_SE(ReSE,ReSE_er,EE,nE,dE,kk,nk,dk,dove)   
     !  integer nE, nk
     ! double precision dk
      ! double precision ReSE(nE,nk), ReSE_er(nE,nk),  EE(nE), dE(nE), kk(nk)
     !  character*64 dove
     !end subroutine real_SE
     
     !subroutine real_SE_ek(ReSE,the_err,Ene,kk,dove)      
     !!  double precision ReSE, Ene, kk, the_err
     !  character*64 dove
     !end subroutine real_SE_ek
     
     !subroutine sigma_hf(r, nr, dr, ell, vr, ki, rho, deg, ReSE)
     !  implicit none
      ! integer nr, i
      ! double precision r(nr), vr(nr), ell(nr)
       !double precision dr, ki, rho, deg, ReSE
     !end subroutine sigma_hf
     
     
  end interface

  
  
  character*128 arg, home, dove, filename, filenum, filenum1, name, name1, varname, workdir 
  integer opt, type
  double precision dens, alpha
  character temp_str
  integer i, k , n, nkp, nkh
  double precision y, E2bc, T2bc, TF, Nyq
  double precision constE1, constE2, constk, kcoef, dE1, dE2
  double precision, parameter:: ahs = 1.d0
  double precision, dimension(:), allocatable:: ell, dell, ddell, fcq, vr
  double precision, dimension(:), allocatable:: kin, pot
  double precision, dimension(:,:), allocatable:: PP, PP_er, HH, HH_er
  double precision, dimension(:,:), allocatable:: PP1, PP1_er, HH1, HH1_er
  double precision Ene, Reh, Reh_er, Rep, Rep_er, Rehf, retot, retot_er
  double precision denom
  double precision sumh, sumh1, sump, sump1,  mdis, mdis1
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
  
  
  open(5,file=trim(dove)//'param.in') 
  
  ! read system parameters
  
  read(5,6100)kF
  read(5,6100)deg
  read(5,6300)a,d,rmax
  read(5,*)opt
  read(5,6100)htm
  read(5,6100)m
  read(5,*)type, alpha

  a = a*ahs
  d = d*ahs
  rmax = rmax* ahs

  
  if(type.eq.1) write(6,*) 'HARD-SPHERES system'
  if(type.eq.2) write(6,*) 'SOFT-SPHERES system'


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


  !.....read correlations
  read(39,*)temp_str,nt,na,nd     
  allocate(r(nt),fc(nt),dfc(nt),ddfc(nt))
  do i = 1, nt
     read(39,6400) r(i),fc(i),dfc(i),ddfc(i)
  end do

  ! define all the extra  parameters (hr, nr) for correlations
  hr = r(2)-r(1)
  nr = nt - nd
  
  allocate(vr(nt),fcq(nt),ell(nt),dell(nt),ddell(nt),veff(nt))
  


  vr(1:nt) = alpha

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
  
  !.....compute the effective interaction [ (hbar**2/m)  ]
  allocate(kin(i), pot(i))
  do i = 1,nt
     kin(i) = htm/m * dfc(i)**2
     pot(i) = fcq(i)*vr(i)
     veff(i) = kin(i) + pot(i)
  end do
  
  open(13,file=trim(adjustl(dove))//'veff_r.dat')
  write(13,*) ' #', nt,kF
  write(13,6200) (r(i), veff(i),i=1, nt)
  close(13)
  
  TF= 3.d0/5.d0*(kF*kF)/2.*htm/m
  E2bc =  T_2bc(na,nd,nt,hr,r,deg,m,htm,dens,ell,dfc)
  write(6,*) ' '
  write(6,*) 'Energy 2body:', TF+E2bc, '   TF=', TF
  write(6,*) ' '
  
  
  constE1 = 0.1d0
  constE2 = 0.5d0
  constk  = 0.5d0
  dk  = constk*kF
  dE1 = constE1*EkF
  dE2 = constE2*EkF
  kcoef = 5.0d0
  kmax =  1.50d0*kF
  kmin = 0.50d0*kF
  Emax =  kcoef*kcoef*EkF
  Emin = -kcoef*kcoef*EkF



  !nk = int((kmax - kmin)/dk + 0.5d0)
  nk = 1
  nEh = int((Emax  - Emin)/dE1 + 0.5d0)
  nEp = int((Emax  - Emin)/dE2 + 0.5d0)
  write(6,*) 'nEh, nEp, dEh, dEp, Emin, Emax', nEh, nEp, dE1, dE2, Emax, Emin
  

  ! allocate grid E & k 
  allocate(kk(nk),EEh(nEh),dEh(nEh),EEp(nEp),dEp(nEp))
  
  if(nk.eq.1) then
     kk(nk) = kF !0.1d0
     nkF = nk
     write(6,*) '------>  nk = 1' , kk(nk), kF
     !     
  else if (kF.gt.kmin) then
     nkh = int((kF-kmin)/dk +0.5d0)
     nkp = nk - nkh
     nkF = nkh
     !do k = 1, nk
     !kk(k) = kmin + (dble(k))*dk
     !write(6,*) k,kk(k)
     !end do
     
     do k = 1, nkh
        kk(k) = kF - dble(nkh-k)*dk
     end do
     do k = nkh + 1 , nk
        kk(k) = kF + dble(k-nkh)*dk
     end do
  else
     nkh = 0
     nkp = nk
     do k = 1, nkp
        kk(k) = kmin + dble(k)*dk
        write(6,*) k,kk(k)
     end do
  end if
  
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
  
  write(6,*) ' nk, dk, kmin, kmax', nk, dk, kmin, kmax
  
  ! write Egrid & kgrid
  open(1,file=trim(adjustl(dove))//'Egridh.dat')
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
  
  ! calcola la trasf do fourier del potenziale
  Nyq =  1.d0/(2.d0*hr)
  pcut = 5.d0
  np = 5000
  allocate(p(np),Mp(np),dMp(np),ddMp(np))
  dp   =  Pcut/dble(np)
  do  i = 1, np
     p(i) = (dble(i)-0.5d0)*dp
  end do
  
  call Fourier(r,veff,nt,hr,p,Mp,np)
  call der(p,1,np,dp,Mp,dMp,ddMp)
  
  open(47,file=trim(adjustl(dove))//'veff_q.dat')
  write(47,*)'#',np,kF
  write(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)
  
  write(filenum, '(6f5.3)')  kF
  

  
  
  
  
!===========================
! SELF ENERGY off shell
!===========================
  
  allocate(ImSEh(nEh,nk),ImSEh_er(nEh,nk),ImSEp(nEp,nk),ImSEp_er(nEp,nk))
  allocate(ReSEh(nEh,nk),ReSEh_er(nEh,nk),ReSEp(nEp,nk),ReSEp_er(nEp,nk))
  
  
  if(mah)  then
     write(6,*) 'aggiungere mahaux'
     
  else
     part = .false.
     !call imaginary_SE(ImSEh,ImSEh_er,EEh,nEh,dEh,kk,nk,dk,dove)
     call real_SE(ReSEh,ReSEh_er,EEh,nEh,dEh,kk,nk,dk,dove)      
     part = .true.
     !call imaginary_SE(ImSEp,ImSEp_er,EEp,nEp,dEp,kk,nk,dk,dove)
     call real_SE(ReSEp,ReSEp_er,EEp,nEp,dEp,kk,nk,dk,dove)      
  end if

!========================
!   SPECTRAL FUNCTIONS
!========================

  allocate(HH(nEh,nk), HH_er(nEh, nk), PP(nEp,nk), PP_er(nEp,nk))
  allocate(HH1(nEh,nk), HH1_er(nEh, nk), PP1(nEp,nk), PP1_er(nEp,nk))

  do k = 1,  nk
     Ene = htm* kk(k)**2/(2.d0*M)
     call sigma_hf(r,nt,hr,ell,veff,kk(k),dens,deg,Rehf)
     part = .true.
     call real_SE_ek(Rep,Rep_er,Ene,kk(k),dove)     
     part = .false.
     call real_SE_ek(Reh,Reh_er,Ene,kk(k),dove)     
     retot    = Rehf + Reh + Rep
     retot_er = sqrt(Reh_er**2 + Rep_er**2)


     write(filenum1, '(6f5.3)')  kk(k)
     name   = trim(adjustl(dove))//'/Spectral/H'//trim(filenum1)//'.dat'
     name1  = trim(adjustl(dove))//'/Spectral/P'//trim(filenum1)//'.dat'
     open(31,file=name)         
     open(32,file=name1)         

     do n = 1, nEh
        denom    = (EEh(n) - Ene - retot)**2 + ImSEh(n,k)**2
        HH(n,k) = ImSEh(n,k)/denom/pi
        !H_er
        denom = (EEh(n) - Ene)**2 + ImSEh(n,k)**2
        HH1(n,k) = ImSEh(n,k)/denom/pi
        !H1_er
        write(31,6700) HH(n,k), HH_er(n,k), HH1(n,k), HH1_er(n,k), kk(k), kk(k)/kF, EEh(n)
     end do

     do n = 1, nEp
        denom = (EEp(n) - Ene - retot)**2 + ImSEp(n,k)**2
        PP(n,k) = ImSEp(n,k)/denom/pi
        !P_er
        denom = (EEp(n) - Ene)**2 + ImSEp(n,k)**2
        PP1(n,k) = ImSEp(n,k)/denom/pi
        !P1_er
        write(32,6700) PP(n,k), PP_er(n,k), PP1(n,k), PP1_er(n,k), kk(k), kk(k)/kF, EEp(n)
     end do
     flush(31)
     flush(32)
     close(31)
     close(32)
  end do






!====================================================
! momentum distribution using spectral functions
!===================================================
  
  
 
  name  = trim(adjustl(dove))//'Spectral/nk_spectral'//trim(filenum)//'.dat'
  open(32,file=name) 
  write(32,*) 'k, nk, norm, nk1, norm1, kk/kF'
  
  dE1 = dEh(1)
  dE2 = dEp(1)     
  do k = 1, nk
     sumh = 0.d0
     sump =0.d0
     sumh1 = 0.d0
     sump1 =0.d0
     do n= 1,nEh
        sumh  = sumh  + HH(n,k)
        sumh1 = sumh1 + HH1(n,k)
     end do
     do n= 1,nEp
        sump  = sump  + PP(n,k)
        sump1 = sump1 + PP1(n,k)
     end do
     
     sumh = sumh * dE1
     sumh1 = sumh1 * dE1

     mdis  = sumh
     mdis1 = sumh1
     sump = sump*dE2
     sump1 = sump1*dE2
     
     write(32,6600) kk(k), mdis,  sumh+sump, mdis1, sumh1+sump1,  kk(k)/kF
     
  end do
  close(32)
  
  
  
  
  
  


  !deallocate
  deallocate(r,fc,dfc,ddfc)
  deallocate(fcq,ell,dell,ddell,veff)
  deallocate(p,Mp,dMp,ddMp)
  deallocate(kk,EEh,dEh,EEp,dEp)
  
  deallocate(ImSEh,ImSEh_er,ImSEp,ImSEp_er)
  deallocate(ReSEh,ReSEh_er,ReSEp,ReSEp_er)
  deallocate(kin, pot)
  deallocate( PP, PP_er, HH, HH_er)
  deallocate( PP1, PP1_er, HH1, HH1_er)
  
  
2222 format(e25.15)
6100 format(e12.4)
6150 format(2e12.4)
6200 format(2e20.10e3)
6300 format(3e12.4)
6400 format(4e20.10)    
6500 format(5e17.8)
6600 format(6e17.8)
6700 format(7e17.8)    
8200 format(a,2i5.1)
8220 format(a,2e12.4)
8300 format(a,3e12.4)      

      stop
  end program main
  
