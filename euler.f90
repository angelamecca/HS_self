module euler
  implicit none
  double precision h
  double precision, dimension(:), allocatable:: r, a00, ap00, adp00, F00
  double precision pi, htm, m
  
contains
 
!===============================================================================
!     Integration of Euler-Lagrange equations for spin-singlet channels
!     using Numerov's algorithm in r(n0+1) < r < r(nc)
!
!     Returns the correlation functions f_00 and f_10   
!
!     Roma - January 2006-----> 2013 modified for hard spheres 
!===============================================================================
!
  subroutine sngl_chann(fc00, n0, nc)
    
    implicit none
    
    integer jr,ilambda, nc, n0
    integer nmaxr 
    double precision lamold,lamnew, lam00, hq12, small
    double precision delta00, deltaold, dg00, diff, fac 
    

    double precision fc00(nc),X00(nc),g00(nc),zz(nc)
    
    !.....set constants
    !
    hq12 = h*h/12.d0
    small = 1.d-12    
    !
    lam00 = 0.d0
    lamold = lam00

    do  ilambda = 1,40
       
       do  jr = n0+1,nc
          X00(jr) = F00(jr) - lam00/htm
       end do
       g00(n0+1) = h
       zz(n0+1) = (1.d0 - hq12*X00(n0+1))*g00(n0+1)
       g00(n0+2) = 12.d0*(g00(n0+1)-zz(n0+1)) + 2.d0*g00(n0+1)
       zz(n0+2) = (1.d0-hq12*X00(n0+2))*g00(n0+2)
       
       do  jr = n0+3,nc
          zz(jr) = 10.d0*(g00(jr-1) - zz(jr-1)) + 2.d0*g00(jr-1) - zz(jr-2)
          g00(jr) = zz(jr)/(1.-hq12*X00(jr))
       end do
       
       fac = a00(nc)/g00(nc)

       do  jr = n0+1,nc
          g00(jr) = fac*g00(jr)
       end do
       dg00 = (3.*g00(nc) - 4.*g00(nc-1) + g00(nc-2))/2.d0/h
       delta00 = a00(nc)*dg00 - ap00(nc)*g00(nc)
       
       
       if(ilambda.eq.1) then
          deltaold = delta00
          lamold = lam00
          lam00 =  0.0001
       else
          diff = abs(delta00-deltaold)
          if(diff.le.small) then 
             write(8,*)' convergence reached in channel 00 '
             write(8,*)' # of iterations ',ilambda,' lambda ',lam00
             write(8,*) 'g e dg, a e da :', g00(nc),dg00,a00(nc),ap00(nc)
            go to 5000
          endif
          lamnew = (delta00*lamold - deltaold*lam00)/(delta00 - deltaold)
          deltaold = delta00
          lamold = lam00
          lam00 = lamnew
       endif
       
    end do
    
 5000 continue

    do  jr = n0+1,nc
       fc00(jr) = g00(jr)/a00(jr)
    end do

    return
  end subroutine sngl_chann
  
!===============================================================================
  




end module euler
