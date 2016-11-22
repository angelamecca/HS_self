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
  use Parameters
  use EffectiveInteraction
  use Boundary
  use selfenergy
  
  implicit none
  
  interface
     double precision function fint(r,f,nr,hr,x,ichoic)
       implicit real*8 (a-h,o-z)
       dimension R(1:*),F(1:*)  
     end function fint
     
     subroutine interpolint(xin,yin,nin,x,y,n,np)
       integer*4 :: nin,n,np,npp,i,j
       real*8 :: xin(nin),yin(nin),x(n),y(n),dy
     end subroutine interpolint

  end interface

  
  
  character*64 arg, home, dove, filename 
  integer opt, type
  double precision dens, alpha
  character temp_str
  integer i, j, k , n, nkp, nkh, i0, k0, j0
  double precision y, E2bc, T2bc, V2bc, TF, Nyq
  double precision constE1, constE2, constk, kcoef, dE1, dE2
  double precision, parameter:: ahs = 1.d0
  double precision, dimension(:), allocatable:: ell, dell, ddell, fcq, vr,  vrnew
  double precision, dimension(:), allocatable:: kin, pot
  integer ste_tmp
  character*64 filenum, name, name1, name2
  integer na_in, nd_in, nr_in, ntot_in
  double precision, dimension(:), allocatable :: inr, inveff_CW, inveff_JF, inveff_BP
  double precision, dimension(:), allocatable :: veff_CW, veff_JF, veff_BP
  double precision :: error
  

  interface
     subroutine FilonFourier(rr,fr,qq,fq,opt)
       implicit none
       logical opt
       real(kind(1.d0)), dimension(:) :: rr, fr, qq, fq 
     end subroutine FilonFourier
  end interface

  mah = .false.
  if(mah) write(6,*)'------>> MAHAUX matrix elements' 
  
  if(iargc().lt.1) then
     write(6,*) 'name of the directory  with input parameters required'
     stop
  end if
  call getarg(1,arg)
  
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
  
  allocate(fcq(nt),ell(nt),dell(nt),ddell(nt),veff(nt),vr(nt))
  
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

  
  !.....compute the effective interaction [ (hbar**2/m)  ]

    open(21,file=trim(adjustl(dove))//'prova_veff.dat')
    read(21,*)
    read(21,*) temp_str, na_in, nd_in, nr_in, ntot_in
    write(*,*) temp_str, na_in, nd_in, nr_in, ntot_in

    allocate(inr(ntot_in), inveff_CW(ntot_in), inveff_JF(ntot_in), inveff_BP(ntot_in))
    read(21,'(1p,4e20.8e3)') (inr(i), inveff_CW(i), inveff_JF(i), inveff_BP(i),  i=1,ntot_in)

    
    allocate(veff_CW(nt), veff_JF(nt),veff_BP(nt))
    
    allocate(kin(nt), pot(nt))
    !do i = 1,nt
    !veff(i) = htm/m * dfc(i)**2 + fcq(i)*vr(i)
    !kin(i) = htm/m * dfc(i)**2
    !pot(i) = fcq(i)*vr(i)
    !veff(i) = kin(i) + pot(i)
    !veff_CW(i) = fint(inr,inveff_CW,ntot_in,r(i),2)
    !veff_JF(i) = fint(inr,inveff_JF,ntot_in,r(i),2)
    !veff_BP(i) = fint(inr,inveff_BP,ntot_in,r(i),2)
    !end do
    call interpolint(inr,inveff_CW,ntot_in,r,veff_CW,nt,4)
    call interpolint(inr,inveff_JF,ntot_in,r,veff_JF,nt,4)
    call interpolint(inr,inveff_BP,ntot_in,r,veff_BP,nt,4)
    
  !
  pot(:) = 0.d0

  ! scrivi veff derivata da CW
  open(13,file=trim(adjustl(dove))//'CW/veff_r.dat')
  write(13,*) ' #', nt,kF
  kin(:) = veff_CW(:)
  do i = 1, nt
     write(13,6400) r(i), veff_CW(i), kin(i), pot(i)
  end do
  close(13)
  !
  !scrivi veff derivata da JF
  open(13,file=trim(adjustl(dove))//'JF/veff_r.dat')
  write(13,*) ' #', nt,kF
  kin(:) = veff_JF(:)
  do i = 1, nt
     write(13,6400) r(i), veff_JF(i), kin(i), pot(i)
  end do
  close(13)

  ! scrivi veff derivata da BP
  open(13,file=trim(adjustl(dove))//'BP/veff_r.dat')
  write(13,*) ' #', nt,kF
  kin(:) = veff_BP(:)
  do i = 1, nt
     write(13,6400) r(i), veff_BP(i), kin(i), pot(i)
  end do
  close(13)


  TF= 3.d0/5.d0*(kF*kF)/2.*htm/m
  T2bc =  T_2bc(na,nd,nt,hr,r,deg,m,htm,dens,ell,dfc)
  V2bc =  V_2bc(na,nd,nt,hr,r,dens,deg,ell,fcq,vr)
  E2bc =  T2bc + V2bc
  write(6,*) ' '
  write(6,*) 'Energy 2body E2bc = ', TF+E2bc, '   TF =', TF, 'T2bc', E2bc, 'V2bc = ', V2bc, '(T+V) 2bc', E2bc
  write(6,*) ' '
  
  
  constE1 = 0.1d0
  constE2 = 0.5d0
  !constk  = 0.075d0
  !constk  = 0.15d0
  constk = 0.10d0
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
  nEh = 1
  nEp = 1
  write(6,*) 'nEh, nEp, dEh, dEp, Emin, Emax', nEh, nEp, dE1, dE2, Emax, Emin
  
  ! allocate grid E & k 
  allocate(kk(nk),EEh(nEh),dEh(nEh),EEp(nEp),dEp(nEp))
  write(6,*) ' nk, dk, kmin, kmax', nk, dk, kmin, kmax
  

  if(nk.eq.1) then
     kk(nk) = kF
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
  ! calcola la trasf di fourier del potenziale
  Nyq =  1.d0/(2.d0*hr)
  pcut = 15.d0*kF
  np = 5000
  allocate(p(np),Mp(np),dMp(np),ddMp(np))
  dp   =  min(Pcut/dble(np),Nyq)
  do  i = 1, np
     p(i) = (dble(i)-0.5d0)*dp
  end do
  !

  write(6,*) 'size', size(r), size(veff), size(p), size(Mp)

  ! transformata di Fouriere di veff_CW
  call FilonFourier(r, veff_CW, p, Mp, .true.)
  call der(p,1,np,dp,Mp,dMp,ddMp)
  open(47,file=trim(adjustl(dove))//'CW/veff_q.dat')
  write(47,*)'#',np,kF
  write(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)

  ! transformata di Fouriere di veff_JF
  call FilonFourier(r, veff_JF, p, Mp, .true.)
  call der(p,1,np,dp,Mp,dMp,ddMp)
  open(47,file=trim(adjustl(dove))//'JF/veff_q.dat')
  write(47,*)'#',np,kF
  write(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)

  ! transformata di Fouriere di veff_BP
  call FilonFourier(r, veff_BP, p, Mp, .true.)
  call der(p,1,np,dp,Mp,dMp,ddMp)
  open(47,file=trim(adjustl(dove))//'BP/veff_q.dat')
  write(47,*)'#',np,kF
  write(47,6400) (p(i),Mp(i),dMp(i),ddMp(i), i=1,np)
  close(47)

  
  allocate(vrnew(nt))

  !   call Fourier(p,Mp,np,dp,r,vrnew,nt)
  call FilonFourier(p, Mp, r, vrnew,.false.)
  
  open(49,file=trim(adjustl(dove))//'antitrasformata.dat')
  write(49,*)'#',np,kF
  do i = 1, nt
     write(49,*) r(i),vrnew(i)
  end do
  close(49)

  deallocate(vrnew)
  


! deallocate
  deallocate(r,fc,dfc,ddfc)
  deallocate(fcq,ell,dell,ddell,veff)
  deallocate(p,Mp,dMp,ddMp)
  deallocate(kk,EEh,dEh,EEp,dEp)


  
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
