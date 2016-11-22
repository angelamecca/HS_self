!=======================================================================
!
!    Hard & Soft Spheres
!    Program to  determine correlation functions by solving
!    Eulero-Lagrange equations
!
!     Input file: "../Data/hs_el.in" with:
!     k_F : fermi momentum
!     deg : degeneracy     
!     alpha : V0 (for soft spheres)
!     a,d,rmax : core, healing distance, maximum distance
!     n1 : # steps per unit lenght
!     ht 
!     m
!     opt : 1 hard spheres, 2 soft spheres, .....
!     Output files:
!     "../Data/fcel.in" with correlations and its first and second derivatives
!  
!     "../Data/param.in" with parameters for my fhnc code:
!      k_F, deg, a,d,rmax,corr_type = 2 (1=euler, 2= sin),
!      htm, m, opt(1=hard, 2=soft), alpha (0 for hard, V0 for soft)
!
!     "../Data/vpot.in" with potential
!
!=======================================================================

program correlation_EL
  use euler

  implicit none
  integer :: n1, na, nd, nc, nt, nr
  integer :: opt, ir, i, corr_type
  double precision :: dens ,r_0, k_F, small, alpha, deg, TF
  double precision :: x, y, z
  double precision :: phi00, phip00, phidp00
  double precision :: a, d, rmax
  double precision :: ahs, htma, ht
  double precision, dimension(:), allocatable:: v00, vnn, fc00
  double precision, dimension(:), allocatable:: df, d2f
  double precision, dimension(:), allocatable:: slf, slpf, sldpf 
  double precision :: slat, slat1, zEper, Eper, c, teff, teff1, veff, veff1
  double precision :: e2bc, t2bc, v2bc, z2bc
  character*128 varname, workdir, home, dove
  integer status_val
!.....set constants
!
  corr_type = 2 !(eulero)
  pi = 4.d0*atan(1.d0)    
  small = 1.e-15

  !output



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
  
  open(49,file=trim(home)//'fcel.in')
  open(48,file=trim(home)//'param.in')
  open(47,file=trim(home)//'vpot.dat')
  
  !.....read input from unit 5
!     
  open(5,file=trim(home)//'param_el.in')
  !
  read(5,*)k_F
  read(5,*)deg      
  read(5,*)alpha
  read(5,*)a,d,rmax
  read(5,*)n1
  read(5,*)ht
  read(5,*)m
  read(5,*)opt,ahs !hard = 1, soft = 2

  ! all the distances in units of ahs
  ! hard spheres: set a = 1; soft spheres set a = range of the interaction
  
  a = a*ahs
  d = d*ahs
  rmax = rmax* ahs

  htm = ht*ht/m
  htma = ht*ht/m/(ahs**2.d0)

  
!     n1: # steps per unit lenght------------> na = n1*a
     
  h = 1.d0/dble(n1)
  na = int(a*dble(n1)+0.5d0)
  nd = int((d-a)*dble(n1)+0.5d0)
  nr = int((rmax-d)*dble(n1)+0.5d0)
  nc = na + nd
  nt = na + nd + nr
  write(6,*) na, nc, nt
  dens = deg * k_F**3/ (6.d0*pi*pi)
  !r_0 = (3./(4.*pi*dens))**(1./3.)


  ! allocate
  allocate(r(nt), vnn(nt))
  allocate(a00(nt), ap00(nt), adp00(nt), F00(nt), v00(nt), fc00(nt) )
  allocate (df(nt), d2f(nt))
  allocate(slf(nt), slpf(nt), sldpf(nt))    
  
  write(6,*)
  write(6,*) ' PARAMETERS : ahs = ', ahs
  write(6,*) ' core =', a, ' healind distance = ', d, 'maximum range =' , rmax, 'c = ', k_F*ahs
  write(6,*) ' kF = ', k_F, '  dens = ', dens, 'deg = ', deg
  write(6,*)
  
  !set grid & potential vnn(r)
  
  do ir= 1, nt     
     r(ir) = dble(ir)*h
  end do

  if (opt.eq.1) then
     write(6,*) ' '
     write(6,*) 'Solving Euler-Lagrange equation for HARD SPHERES'
     write(6,*) ' '
     do ir = 1, nt
        vnn(ir) = 0.d0
     end do
  else if(opt.eq.2) then
     !.... define step-function potential for soft spheres V0 = alpha * theta(a-r)
     write(6,*) ' '
     write(6,*) ' Solving Euler-Lagrange equation for SOFT SPHERES'
     write(6,*) ' '     
     do ir= 1, na     
        vnn(ir) = alpha
     end do
     do ir= na+1, nt     
        vnn(ir) = 0.d0
     end do
  else
     write(6,*) ' error! set opt value  = 1 for  hard spheres,  opt = 2  for soft spheres'
     return
  end if
  
  
  !slater and derivatives  
  do ir= 1, nt
     x = k_F*r(ir)
     if(x.eq.0.) then
        slf(ir)   = 1.d0
        slpf(ir)  = 0.d0
        sldpf(ir) = -k_F*k_F *0.2d0
     else
        y=sin(x)
        z=cos(x)
        slf(ir)   = 3.d0*(y-x*z)/x**3.
        slpf(ir)  = 3.d0*(y/x-slf(ir))/r(ir)
        sldpf(ir) = 3.d0*(-slpf(ir) + (slf(ir)+z-2.d0*y/x )/r(ir))/r(ir)
     end if
     !a(r), see FFPR nuovo cimento
     
     phi00     =  sqrt( 1.d0 - slf(ir)**2./deg )
     phip00    =  -slf(ir)*slpf(ir)/phi00 /deg
     phidp00   = -( ((slf(ir)*slpf(ir)/phi00)**2.)/deg + slpf(ir)**2 + slf(ir)*sldpf(ir) ) /phi00/deg
     a00(ir)   =  r(ir)*phi00
     ap00(ir)  =  phi00 + r(ir)*phip00
     adp00(ir) =  2.d0*phip00 + r(ir)*phidp00
          
!....set potential     
     v00(ir) =   vnn(ir)

!.....initialize FT0 , F1T1, F2T1 and F3T1
     F00(ir) = adp00(ir)/a00(ir) + v00(ir)*m/htm   

     write(47,*), r(ir), v00(ir)

  end do

     
     
!.....Solve Eulero-Lagrange equation in 0 < r < d   
!.....Nb: sng_chann solves EL equation for r in [r(nstart+1), r(nfin)]
  
  if (opt.eq.1) then
     !f(r) = 0  if r <= a:  r(nstart+1) =r(na+1), nstart must be = na     
     do ir = 1,na
        fc00(ir) = 0.d0
        df(ir) = 0.d0
        d2f(ir) = 0.d0
     end do
     call sngl_chann(fc00, na, nc) 
     call der(r,na+1,nc,h,fc00,df,d2f)  
     
  else if(opt.eq.2) then
     ! for SOFT SPHERES r(nstart+1) = r(1), nstart must be 0 
     call sngl_chann(fc00, 0, nc) 
     call der(r,1,nc,h,fc00,df,d2f)
  end if
  
  !f(r) = 1  if r > d      
  
  do ir = nc+1,nt
     fc00(ir) = 1.d0
     df(ir)   = 0.d0
     d2f(ir)  = 0.d0
  end do
  write(6,*) ' '
  write(6,*)'check boundary conditions for correlation function'
  write(6,*) ' '

  if (opt.eq.1) then
     write(6,*) "fc(a), fc(a + h),fc(d), fc'(d)", fc00(na),fc00(na+1), fc00(nc), df(nc)
  else
     write(6,*) "fc(h), fc(d), fc'(d)", fc00(1), fc00(nc), df(nc)
  end if
  
  !     
  !write input parameters  for g_fhnc(r) in param.in       
  !     
  write(48,'(1p,e12.4)')k_F
  write(48,'(1p,e12.4)')deg
  write(48,'(1p,3e12.4)')a/ahs,d/ahs,rmax/ahs
  write(48,*)corr_type
  write(48,'(1p,e12.4)')ht
  write(48,'(1p,e12.4)')m
  write(48,*)opt,ahs,alpha
  
  !write f(r) in fcel.in
  write(49,*)'#',nt,na,nc
  write(49,'(1p,4e20.10)')(r(i),fc00(i),df(i),d2f(i),i=1,nt)
  
  !.... compute energy: TF
  TF= 3.d0/5.d0*(k_F*k_F)/2.d0*htm
  
  ! perturbative g.s. energy
  if(opt.eq.1) then
     c = k_F * a 
  else
     c = k_F * ahs
  end if

  
  if(deg.eq.2.) then
     zEper=  5.d0/3.d0 *( 2.d0*c/(3.d0*pi) + 4.d0/(35.d0*pi**2)*(11.d0- 2.d0*log(2.))*c**2+ 0.230*c**3 )
  else
     !.....deg=4      
     zEper= 5./3. *( 2.*c/pi + 12./(35.*pi**2)*(11.- 2.*log(2.))*c**2+ 0.780*c**3 + &
          32./(9.*pi**3)*(4.*pi-3.*3**0.5)*c**4*log(c) )
  end if

  Eper = TF*(1.d0 + zEper)
  write(6,*) ''
  write(6,*) 'Low density espansion in c=' ,c, ' correction z=', zEper, ' E/TF = ', (1.d0+ zEper)
  write(6,*) ''

  

! energy 2 body cluster  
  t2bc = 0.d0
  v2bc = 0.d0

  do i = 1,nc-1
     slat  = (1.d0 - slf(i)*slf(i)   /deg)* r(i) *r(i)
     slat1 = (1.d0 - slf(i+1)*slf(i+1)/deg)*r(i+1)*r(i+1) 

     veff = v00(i)* fc00(i)**2 * slat
     teff = df(i)**2  * slat 

     veff1 = v00(i+1)* fc00(i+1)**2 *     slat1
     teff1 = df(i+1)**2   *slat1

     v2bc = v2bc + (veff + veff1)/2.d0 
     t2bc = t2bc + (teff + teff1)/2.d0

  end do

  v2bc= 2.d0*pi*dens*v2bc*h*htm
  t2bc= 2.d0*pi*dens*t2bc*h*htm
  e2bc = t2bc + v2bc

  z2bc = (TF+e2bc)/TF - 1.d0



  write(6,*) ' '
  write(6,*) ' a, d, kF, rho, c = (kF*a)', a/ahs,  d/ahs, k_F, dens, c
  write(6,*) ' '
  write(6,*) ' =========  2 body ================'
  write(6,*) ' '
  write(6,*) ' t 2bc =  ',t2bc, 'v 2bc =  ', v2bc, ' e_tot 2bc =   ', e2bc, 'z 2bc = ', z2bc
  write(6,*) ' '
  write(6,*) ' ========= low-density  =========='
  write(6,*) ' '
  write(6,*) ' E pert = ' , Eper, 'z pert =', zEper
  write(6,*) ' '
  write(6,*) ' ========= energy   =========='
  write(6,*) ' '
  write(6,*) 'TF = ', TF, '  E tot 2bc  =    ' , TF + e2bc, '  E pert = ', Eper
  write(6,*) ''
  write(6,*) ' ========= energy  in unit ht^2/ma^2  =========='
  write(6,*) ' '
  write(6,*) 'c = ' , c, 'TF = ', TF, '  T tot 2bc  =    ' , (TF + e2bc)*htma/htm, 'Etot 2bc = ', (TF + e2bc + v2bc)*htma/htm, '  E pert = ', Eper*htma/htm
  write(6,*) ''
  write(6,*) ''
  write(6,*) ' ========= z (c)  =========='
  write(6,*) ''
  write(6,*) ' c= ', c, 'z_2bc = ', z2bc, 'z_pert', zEper
  write(6,*) ''
  write(6,*) ' ======== a, d, kF,  c=kF*a, deg, TF, e2bc, T2bc, Epert, z2bc, zEper ======== '
  write(6,*) ' '
  write(6,'(1p,12e20.8e2)')  a/ahs, d/ahs, k_F, c, deg, TF, e2bc*htma/htm, T2bc*htma/htm, Eper*htma/htm, V2bc*htma/htm, z2bc, zEper
  write(6,*) ' '




  

  stop
  
  return     

end program correlation_EL
!***********************************************************************

subroutine delsq(f,r,hr,nr,delsq_f)
  
  implicit real*8(a-h,o-z)
  
  dimension r(100000),f(100000),delsq_f(100000),df(100000), d2f(100000)
  
  call der (r,1,nr,hr,f,df,d2f)
  do  i = 1,nr
     delsq_f(i) = d2f(i) + 2.*df(i)/r(i)
  end do
  return
end subroutine delsq

!***********************************************************************

subroutine der(x,n_0,n,h,f,df,d2f)
      
  !     returns first (df) and second (d2f) derivative of the 
  !     function f. The input function is given in the n points
  !     x(n_0),....,x(n) with steplength h.
  
  implicit real*8(a-h,o-z)
  
  dimension x(100000),f(100000),df(100000),d2f(100000)
  
  nm3 = n - 3
  n0p3 = n_0 + 3
  hr = 1./60./h
  h2r = hr/3./h
  hqr = 60.*h2r
  hrm = 10.*hr
  
  f0 = f(n_0)
  f1 = f(n_0 + 1)
  f2 = f(n_0 + 2)
  f3 = f(n_0 + 3)
  f4 = f(n_0 + 4)
  f5 = f(n_0 + 5)
  f6 = f(n_0 + 6)
  df(n_0) = hr*(-147.*f0+360.*f1-450.*f2+400.*f3-225.*f4 &
       +72.*f5-10.*f6)
  d2f(n_0) = h2r*(812.*f0-3132.*f1+5265.*f2-5080.*f3+2970.*f4 &
       -972.*f5+137.*f6)
  df(n_0 + 1) = hr*(-10.*f0-77.*f1+150.*f2-100.*f3+50.*f4 &
       -15.*f5+2.*f6)
  d2f(n_0 +1) = h2r*(137.*f0-147.*f1-255.*f2+470.*f3-285.*f4 &
       +93.*f5-13.*f6)
  df(n_0 + 2) = hr*(2.*f0-24.*f1-35.*f2+80.*f3-30.*f4 &
       +8.*f5-f6)
  d2f(n_0 +2) = h2r*(-13.*f0+228.*f1-420.*f2+200.*f3+15.*f4  &
       -12.*f5+2.*f6)
  
  do  i = n0p3,nm3
     f1 = f(i + 1)
     f2 = f(i + 2)
     f3 = f(i + 3)
     f4 = f(i - 1)
     f5 = f(i - 2)
     f6 = f(i - 3)
     f7 = f(i)
     
     df(i) = hr*(45.*(f1-f4)+9.*(f5-f2)+f3-f6)
     d2f(i) = h2r*(-490.*f7+270.*(f1+f4)-27.*(f2+f5) &
          +2.*(f3+f6))
  end do
  f1 = f(n - 6)
  f2 = f(n - 5)
  f3 = f(n - 4)
  f4 = f(n - 3)
  f5 = f(n - 2)
  f6 = f(n - 1)
  f7 = f(n)
  
  df(n-2) = hrm*(-f7+6.*f6-3.*f5-2.*f4)
  d2f(n-2) = hqr*(3.*f6-6.*f5+3.*f4)
  df(n-1) = hrm*(2.*f7+3.*f6-6.*f5+f4)
  d2f(n-1) = hqr*(3.*f7-6.*f6+3.*f5)
  df(n) = hrm*(11.*f7-18.*f6+9.*f5-2.*f4)
  d2f(n) = hqr*(6.*f7-15.*f6+12.*f5-3.*f4)
  
  return
end subroutine der
!***********************************************************************
