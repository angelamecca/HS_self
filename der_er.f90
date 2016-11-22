!
!=======================================================================
!
!     5 points der
!
!=======================================================================
!
subroutine der5p_err(x,Nx,dx,f,f_er,df,df_er,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), f_er(Nx), df, df_er
  double precision dx, hr
  integer i
  double precision f0, f1, f2, f3, f4, fe0, fe1, fe2, fe3, fe4
  !     
  hr = 1.d0/12.d0/dx

  if (n0.eq.3) then 
     f0 = f(n0)
     f1 = f(n0 + 2)
     f2 = f(n0 + 1)
     f3 = f(n0 - 1)
     f4 = f(n0 - 2)
     fe0 = f_er(n0)
     fe1 = f_er(n0 + 2)
     fe2 = f_er(n0 + 1)
     fe3 = f_er(n0 - 1)
     fe4 = f_er(n0 - 2)
     !
     df    = hr*( -f1 + 8.d0*(f2-f3) + f4 )
     df_er = hr*sqrt(fe1**2 + fe4**2 + 64.d0*(fe2**2 + fe3**2))
     
  else if(n0.eq.1) then
     f0  = f(n0)
     f1  = f(n0 + 1)
     f2  = f(n0 + 2)
     f3  = f(n0 + 3)
     f4  = f(n0 + 4)
     fe0 = f_er(n0)
     fe1 = f_er(n0 + 1)
     fe2 = f_er(n0 + 2)
     fe3 = f_er(n0 + 3)
     fe4 = f_er(n0 + 4)
     
     df    = hr *(-25.d0*f0 + 48.d0*f1 - 36.d0*f2 + 16.d0*f3 - 3.d0*f4)
     df_er = hr*sqrt( 625.d0*fe0**2 + 2304.d0*fe1**2 + 1296.d0*fe2**2 + 256.d0*fe3**2 + 9.d0*fe4**2)

  else if(n0.eq.2) then
     f0  = f(n0 - 1)
     f1  = f(n0)
     f2  = f(n0 + 1)
     f3  = f(n0 + 2)
     f4  = f(n0 + 3)
     fe0 = f_er(n0 - 1)
     fe1 = f_er(n0)
     fe2 = f_er(n0 + 1)
     fe3 = f_er(n0 + 2)
     fe4 = f_er(n0 + 3)
     df    = hr *(-3.d0*f0 - 10.d0*f1 + 18.d0*f2 - 6.d0*f3 + f4)
     df_er = hr*sqrt( 9.d0*fe0**2 + 100.d0*fe1**2 + 324.d0*fe2**2 + 36.d0*fe3**2 + fe4**2)
  else if (n0.eq.4) then
     f0  = f(n0 - 3)
     f1  = f(n0 - 2)
     f2  = f(n0 - 1)
     f3  = f(n0)
     f4  = f(n0 + 1)
     fe0 = f_er(n0 - 3)
     fe1 = f_er(n0 - 2)
     fe2 = f_er(n0 - 1)
     fe3 = f_er(n0)
     fe4 = f_er(n0 + 1)
     
     df    = hr *(-f0 + 6.d0*f1 - 18.d0*f2 + 10.d0*f3 + 3.d0*f4)
     df_er = hr*sqrt( fe0**2 + 36.d0*fe1**2 + 324.d0*fe2**2 + 100.d0*fe3**2 + 9.d0* fe4**2)
  else
     f0  = f(n0)
     f1  = f(n0 - 1)
     f2  = f(n0 - 2)
     f3  = f(n0 - 3)
     f4  = f(n0 - 4)
     fe0 = f_er(n0)
     fe1 = f_er(n0 - 1)
     fe2 = f_er(n0 - 2)
     fe3 = f_er(n0 - 3)
     fe4 = f_er(n0 - 4)
     
     df    = hr *(25.d0*f0 - 48.d0*f1 + 36.d0*f2 - 16.d0*f3 + 3.d0*f4)
     df_er = hr*sqrt( 625.d0*fe0**2 + 2304.d0*fe1**2 + 1296.d0*fe2**2 + 256.d0*fe3**2 + 9.d0*fe4**2)
     
  end if
  return  
end subroutine der5p_err
!
!=======================================================================
!
subroutine deriv5p(x,Nx,dx,f,df,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), df
  double precision dx, hr
  integer i
  double precision f0, f1, f2,f3,f4
  
  hr = 1.d0/12.d0/dx

  hr = 1.d0/12.d0/dx

  if (n0.eq.3) then 
     f0 = f(n0)
     f1 = f(n0 + 2)
     f2 = f(n0 + 1)
     f3 = f(n0 - 1)
     f4 = f(n0 - 2)
     df    = hr*( -f1 + 8.d0*(f2-f3) + f4 )
  else if(n0.eq.1) then
     f0  = f(n0)
     f1  = f(n0 + 1)
     f2  = f(n0 + 2)
     f3  = f(n0 + 3)
     f4  = f(n0 + 4)
     df    = hr *(-25.d0*f0 + 48.d0*f1 - 36.d0*f2 + 16.d0*f3 - 3.d0*f4)
  else if(n0.eq.2) then
     f0  = f(n0 - 1)
     f1  = f(n0)
     f2  = f(n0 + 1)
     f3  = f(n0 + 2)
     f4  = f(n0 + 3)
     df    = hr *(-3.d0*f0 - 10.d0*f1 + 18.d0*f2 - 6.d0*f3 + f4)
  else if (n0.eq.4) then
     f0  = f(n0 - 3)
     f1  = f(n0 - 2)
     f2  = f(n0 - 1)
     f3  = f(n0)
     f4  = f(n0 + 1)
     df    = hr *(-f0 + 6.d0*f1 - 18.d0*f2 + 10.d0*f3 + 3.d0*f4)
  else
     f0  = f(n0)
     f1  = f(n0 - 1)
     f2  = f(n0 - 2)
     f3  = f(n0 - 3)
     f4  = f(n0 - 4)
     df    = hr *(25.d0*f0 - 48.d0*f1 + 36.d0*f2 - 16.d0*f3 + 3.d0*f4)     
  end if
  
  return
end subroutine deriv5p
!
!=======================================================================
!
!     7 points der
!
!=======================================================================
!

subroutine der7p_err(x,Nx,dx,f,f_er,df,df_er,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), f_er(Nx), df, df_er
  double precision dx, hr
  integer i
  double precision f1, f2, f3, f4, f5, f6, f7
  double precision fe1, fe2, fe3, fe4, fe5, fe6, fe7
  !     
  hr = 1.d0/60.d0/dx
  f1 = f(n0 + 1)
  f2 = f(n0 + 2)
  f3 = f(n0 + 3)
  f4 = f(n0 - 1)
  f5 = f(n0 - 2)
  f6 = f(n0 - 3)
  f7 = f(n0)
  !
  fe1 = f_er(n0 + 1)
  fe2 = f_er(n0 + 2)
  fe3 = f_er(n0 + 3)
  fe4 = f_er(n0 - 1)
  fe5 = f_er(n0 - 2)
  fe6 = f_er(n0 - 3)
  fe7 = f_er(n0)
  !
  df = hr*(45.d0*(f1-f4)+9.d0*(f5-f2) + f3-f6 )
  df_er = hr*sqrt(45.d0**2 *( fe1**2 + fe4**2 ) + 9.d0**2*( fe5**2 + fe2**2) &
       +  fe3**2 + fe6**2)
  return
end subroutine der7p_err
!
!=======================================================================
!

subroutine deriv7p(x,Nx,dx,f,df,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), df
  double precision dx, hr
  integer i
  double precision f1, f2, f3, f4, f5, f6, f7
  !     
  hr = 1.d0/60.d0/dx
  f1 = f(n0 + 1)
  f2 = f(n0 + 2)
  f3 = f(n0 + 3)
  f4 = f(n0 - 1)
  f5 = f(n0 - 2)
  f6 = f(n0 - 3)
  f7 = f(n0)
  !
  !
  df = hr*(45.d0*(f1-f4)+9.d0*(f5-f2) + f3-f6 )
  
  return
end subroutine deriv7p
!
!=======================================================================

!     3 points der
!
!=======================================================================
!

subroutine der3p_err(x,Nx,dx,f,f_er,df,df_er,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), f_er(Nx), df, df_er
  double precision dx, hr
  integer i
  double precision f1, f2, f3
  double precision fe1, fe2, fe3
  !     
  hr = 1.d0/2.d0/dx
  if(n0.eq.1) then
     f1  = f(n0)
     f2  = f(n0 + 1)
     f3  = f(n0 + 2)
     fe1 = f_er(n0)
     fe2 = f_er(n0 + 1)
     fe3 = f_er(n0 + 2)
     df = hr*( -3.d0* f1 + 4.d0*f2  - f3)
     df_er = hr*sqrt(9.d0 * fe1**2 + 16.d0 * fe2**2 +  fe3**2)
  else
     f1 = f(n0 + 1)
     f2 = f(n0 - 1)
     f3 = f(n0)
     
     fe1 = f_er(n0 + 1)
     fe2 = f_er(n0 - 1)
     fe3 = f_er(n0)
     
     df = hr*( f1 - f2)
     df_er = hr*sqrt(fe1**2 + fe2**2)
  end if
  return
end subroutine der3p_err
!
!=======================================================================
!

subroutine deriv3p(x,Nx,dx,f,df,n0)
  implicit none
  integer Nx, n0
  double precision x(Nx), f(Nx), df
  double precision dx, hr
  integer i
  double precision f1, f2, f3
  !     
  hr = 1.d0/2.d0/dx
  if(n0.eq.1) then
     f1  = f(n0)
     f2  = f(n0 + 1)
     f3  = f(n0 + 2)
     df = hr*( -3.d0* f1 + 4.d0*f2  - f3)
  else
     f1 = f(n0 + 1)
     f2 = f(n0 - 1)
     f3 = f(n0)
     df = hr*( f1-f2 )
  end if
  return
end subroutine deriv3p
!
!=======================================================================
!
