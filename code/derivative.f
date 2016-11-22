c     
c=======================================================================
c
      subroutine derivative1(x, fx, dx_fx,dx,Nx)
      integer Nx
      real*8 x(Nx), fx(Nx), dx_fx(Nx)
      real*8 dx
      integer n
c     
c     dx= p(2)-p(1)
c     n0= int(1.+ (k_F-p(1))/dp)
c     write(6,*) dp, n0, p(n0), f_K
      do n = 3, Nx-2 
         dx_fx(n)= (-fx(n+2)+8.*fx(n+1)-8.*fx(n-1)+fx(n-2))/(12.*dx)
      end do
      dx_fx(1)= (fx(2)-fx(1))/dx  
      dx_fx(2) = (fx(3)-fx(1))/(2.*dx)
      dx_fx(Nx-1) = (fx(Nx)-fx(Nx-2))/(2.*dx)
      dx_fx(Nx)= (fx(Nx-1)-fx(Nx))/(-dx)
      end subroutine derivative1
c     
c=======================================================================
c
      subroutine delsq(f,r,hr,nr,delsq_f)
c     
      implicit real*8(a-h,o-z)
c     
      dimension r(nr),f(nr),delsq_f(nr),df(nr),d2f(nr)
c     
      call der(r,1,nr,hr,f,df,d2f)
c     
      do 1000 i = 1,nr
 1000    delsq_f(i) = d2f(i) + 2.*df(i)/r(i)
c
      return
      end
c     
c=======================================================================
c
      subroutine der(x,n_0,n,h,f,df,d2f)
c
c     returns first (df) and second (d2f) derivative of the 
c     function f. The input function is given in the n points
c     x(n_0),....,x(n) with steplength h.
c
      implicit real*8(a-h,o-z)
c
      dimension x(n),f(n),df(n),d2f(n)
c
      nm3 = n - 3
      n0p3 = n_0 + 3
      hr = 1./60./h
      h2r = hr/3./h
      hqr = 60.*h2r
      hrm = 10.*hr
c
      f0 = f(n_0)
      f1 = f(n_0 + 1)
      f2 = f(n_0 + 2)
      f3 = f(n_0 + 3)
      f4 = f(n_0 + 4)
      f5 = f(n_0 + 5)
      f6 = f(n_0 + 6)
c
      df(n_0) = hr*(-147.*f0+360.*f1-450.*f2+400.*f3-225.*f4
     $         +72.*f5-10.*f6)
      d2f(n_0) = h2r*(812.*f0-3132.*f1+5265.*f2-5080.*f3+2970.*f4
     $          -972.*f5+137.*f6)
      df(n_0 + 1) = hr*(-10.*f0-77.*f1+150.*f2-100.*f3+50.*f4
     $         -15.*f5+2.*f6)
      d2f(n_0 +1) = h2r*(137.*f0-147.*f1-255.*f2+470.*f3-285.*f4
     $          +93.*f5-13.*f6)
      df(n_0 + 2) = hr*(2.*f0-24.*f1-35.*f2+80.*f3-30.*f4
     $         +8.*f5-f6)
      d2f(n_0 +2) = h2r*(-13.*f0+228.*f1-420.*f2+200.*f3+15.*f4
     $          -12.*f5+2.*f6)
c
      do 1000 i = n0p3,nm3
c
      f1 = f(i + 1)
      f2 = f(i + 2)
      f3 = f(i + 3)
      f4 = f(i - 1)
      f5 = f(i - 2)
      f6 = f(i - 3)
      f7 = f(i)
c
      df(i) = hr*(45.*(f1-f4)+9.*(f5-f2)+f3-f6)
      d2f(i) = h2r*(-490.*f7+270.*(f1+f4)-27.*(f2+f5)
     $        +2.*(f3+f6))
 1000 continue
c
      f1 = f(n - 6)
      f2 = f(n - 5)
      f3 = f(n - 4)
      f4 = f(n - 3)
      f5 = f(n - 2)
      f6 = f(n - 1)
      f7 = f(n)
c
      df(n-2) = hrm*(-f7+6.*f6-3.*f5-2.*f4)
      d2f(n-2) = hqr*(3.*f6-6.*f5+3.*f4)
      df(n-1) = hrm*(2.*f7+3.*f6-6.*f5+f4)
      d2f(n-1) = hqr*(3.*f7-6.*f6+3.*f5)
      df(n) = hrm*(11.*f7-18.*f6+9.*f5-2.*f4)
      d2f(n) = hqr*(6.*f7-15.*f6+12.*f5-3.*f4)
c
      return
      end
c***********************************************************************
      subroutine der_3p(n_0,n,h,x,f,df)
      implicit none
      integer n_0, n
      real*8 h, f(n),df(n),x(n)
      integer i
c     
      do i=n_0+1,n-1
         df(i) = ( f(i+1)-f(i-1) )/(x(i+1)-x(i-1))
      end do
      df(n_0)  = ( 4.d0*f(n_0+1)-3.d0*f(n_0)-f(n_0+2) )/(2.d0*h)
      df(n)    = (-4.d0*f(n-1)  +3.d0*f(n)  +f(n-2) )/(2.d0*h)
c     
      return
      end subroutine
c***********************************************************************
      subroutine der3p(x,n_0,n,h,f,df,ddf)
      implicit none
      integer n_0, n
      real*8 h, f(n),df(n),ddf(n),x(n)
      integer i
c     
      do i=n_0+1,n-1
         df(i)  = ( f(i+1)-f(i-1) )/(2.d0*h)
         ddf(i) = ( f(i+1) + f(i-1) -2.d0*f(i))/h**2
      end do
      df(n_0)  = ( 4.d0*f(n_0+1)-3.d0*f(n_0)-f(n_0+2) )/(2.d0*h)
      df(n)    = (-4.d0*f(n-1)  +3.d0*f(n)  +f(n-2) )/(2.d0*h)
c
      ddf(n_0) = (f(n_0+2) - 2.d0*f(n_0+1) + f(n_0))/h**2
      ddf(n)   = (f(n)    -  2.d0*f(n-1)   + f(n-2) )/h**2
c     
      return
      end subroutine
c***********************************************************************
      subroutine der3p_asym(x,n0,n1,nt,h1,h2,f,df,ddf)
      implicit none
      integer n0, n1, nt
      real*8 h1,h2, f(nt),df(nt),ddf(nt),x(nt)
      integer i
c     
      df(n0)  = ( 4.d0*f(n0+1)-3.d0*f(n0)-f(n0+2) )/(2.d0*h1)
      ddf(n0) = (f(n0+2) - 2.d0*f(n0+1) + f(n0))/h1**2
c
      do i=n0+1,n1-1
         df(i)  = ( f(i+1)-f(i-1) )/(2.d0*h1)
         ddf(i) = ( f(i+1) + f(i-1) -2.d0*f(i))/h1**2
      end do
c
      df(n1) = ( f(n1+1)-f(n1-1) )/(h1+h2)
      ddf(n1) = ( f(n1+1)*h1 + f(n1-1)*h2 - f(n1)*(h1+h1) )
     &     *2.d0 / ( h1 * h2 * (h1+h2))
c
      do i = n1+1,nt-1
         df(i)  = ( f(i+1)-f(i-1) )/(2.d0*h2)
         ddf(i) = ( f(i+1) + f(i-1) -2.d0*f(i))/h2**2
      end do     
c
      df(nt)    = (-4.d0*f(nt-1)  +3.d0*f(nt)  +f(nt-2) )/(2.d0*h2)
      ddf(nt)   = (f(nt)    -  2.d0*f(nt-1)   + f(nt-2) )/h2**2
c     
      return
      end subroutine
c***********************************************************************
