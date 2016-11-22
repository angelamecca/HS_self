!=======================================================================
!
!     Fourier Tranform
!     
!=======================================================================
!     
      subroutine Fourier(r,fr,nr,dr,q,fq,nq)
      implicit none
      integer nr, nq, n, i
      double precision r(nr),fr(nr),q(nq), fq(nq)
      double precision dr, dq, Iq, pi
!
      
      pi = 4.*atan(1.d0)
      
      do  n = 1, nq
         Iq = 0.d0
         if(q(n).ne.0.d0) then
            do i = 1, nr
               Iq = Iq + fr(i)*r(i)*( sin( q(n)*r(i)) / q(n) ) * dr
            end do
         else
            do i = 1, nr
               Iq = Iq + fr(i)*r(i)**2 * dr
            end do
         end if
         fq(n)=  4.d0*pi*Iq
      end do      
      
      end subroutine Fourier
!
!=======================================================================
!
