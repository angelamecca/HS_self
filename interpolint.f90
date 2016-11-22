      subroutine interpolint(xin,yin,nin,x,y,n,np)
      integer*4 :: nin,n,np,npp,i,j
      real*8 :: xin(nin),yin(nin),x(n),y(n),dy
      
      do i=1,n
         call locate(xin,nin,x(i),j)
         if(j.eq.0) j=1
         npp=min(np,iabs(nin-j))
         if (npp.eq.nin) y(i)=yin(nin)
         call polint(xin(j:j+npp),yin(j:j+npp),npp+1,x(i),y(i),dy)
      enddo
      end subroutine interpolint

!===================================================================================
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a value y,
! and an error estimate dy. If P (x) is the polynomial of degree N − 1 such that 
!P(xai) = yai,i = 1,...,n, then the returned value y = P(x).
!===================================================================================
      subroutine polint(xa,ya,n,x,y,dy)
      integer*4,parameter :: NMAX=10
      integer*4 :: n,i,m,ns
      real*8 :: dy,x,y,xa(n),ya(n)
      real*8 :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.)stop 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo
      return
      end subroutine polint

!===================================================================================
!Given an array xx(1:n), and given a value x, returns a value j such that x is between xx(j)
! and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0 or j=n is 
!returned to indicate that x is out of range. 
!===================================================================================
      subroutine locate(xx,n,x,j) 
      integer*4 :: j,n,jl,jm,ju
      real*8 :: x,xx(n)

      jl=0 
      ju=n+1
      do while (ju-jl.gt.1)  
         jm=(ju+jl)/2
         if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
            jl=jm 
         else
            ju=jm 
         endif
      enddo
      if (x.eq.xx(1))then
         j=1
      else if(x.eq.xx(n))then
         j=n-1 
      else
         j=jl 
      endif
      return
      end subroutine locate
