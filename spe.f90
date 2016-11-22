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
       character*128 dove
     end subroutine im_SE_ek

     subroutine real_SE_ek(ReSE,the_err,Ene,kk,dove)      
       double precision ReSE, Ene, kk, the_err
       character*128 dove
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
  double precision, parameter:: ahs = 1.d0
  double precision, dimension(:), allocatable:: ell, dell, ddell, fcq, vr
  double precision, dimension(:), allocatable:: Ene
  double precision, dimension(:), allocatable:: Reh, Rep, Reh_er, Rep_er
  double precision, dimension(:), allocatable:: Imh, Imp, Imh_er, Imp_er
  double precision, dimension(:), allocatable:: Rehf
  double precision, dimension(:), allocatable:: spe, Retot, Eoff
  double precision, dimension(:), allocatable:: kin, pot
  double precision e0, spein, ein, eout, ris
  integer, parameter :: jmax = 5
  integer, parameter :: stencil = 5
  integer ste_tmp
  character*128 filenum, name, name1, name2

  double precision  alpha
  double precision Ehf, Uk, Ek
  integer  maxit
  double precision small, realh, realh_er, realp, realp_er, dE
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
  ! solve e(k) = e0(k) + Re(k, e(k))
  !===========================================
  allocate(Ene(nk),spe(k))
  allocate(Reh(nk),Rep(nk),Reh_er(nk),Rep_er(nk), Rehf(nk))
  allocate(Imh(nk),Imp(nk),Imh_er(nk),Imp_er(nk))
  allocate(Retot(jmax), Eoff(jmax))

  
  name=trim(adjustl(dove))//"Results/spe_new"//trim(filenum)//".dat"
  name1=trim(adjustl(dove))//"Results/converg"//trim(filenum)//".dat"

  open(22,file=name)
  open(23,file=name1)

  write(22,*) '#  k, e0, spe_old, spe_new, diff, diff %'

  realh=0.d0
  realh_er = 0.d0
  realp = 0.d0
  realp_er = 0.d0

  
  j0 = int(dble(jmax+1)/2.d0)
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
        call flush(23)
        call flush(22)
        call real_SE_ek(realh,realh_er,Eoff(j),kk(k),dove)      
        part = .true.
        call flush(23)
        call flush(22)
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
  !  ON SHELL  Sigma2p2h_Im & Sigma2p2h_Re
  !===========================================
  
    
  name=trim(adjustl(dove))//"Results/self_new"//trim(filenum)//".dat"
  open(22,file=name)
  write(22,*) '# kk, Imh, Imh_er, Reh,Reh_er, Imp, Imp_er, Rep, Rep_er, Rehf '
  flush(22)
     

  do k = 1, nk
     part = .false.
     call real_SE_ek(Reh(k),Reh_er(k),spe(k),kk(k),dove)   
     flush(22)
     call im_SE_ek(Imh(k),Imh_er(k),spe(k),kk(k),dove)      
     flush(22)
     part = .true.
     call real_SE_ek(Rep(k),Rep_er(k),spe(k),kk(k),dove) 
     flush(22)     
     call im_SE_ek(Imp(k),Imp_er(k),spe(k),kk(k),dove) 
     flush(22)     
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


  


  





  


! deallocate
  deallocate(r,fc,dfc,ddfc)
  deallocate(fcq,ell,dell,ddell,veff)
  deallocate(p,Mp,dMp,ddMp)
  deallocate(kk,EEh,dEh,EEp,dEp)
  deallocate(Ene)
  deallocate(Reh,Rep,Reh_er,Rep_er)
  deallocate(Imh,Imp,Imh_er,Imp_er)
  deallocate(Rehf)
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
