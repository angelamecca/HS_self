!=======================================================================
!
!     Fourier Transform
!     
!=======================================================================
!     
! opt = true  Fourier transform
! opt = false Inverse Fourier transform
!=======================================================================
!

subroutine FilonFourier(rr,fr,qq,fq,opt)
  implicit none
  logical opt
!  real(kind(1.d0)), dimension(:) :: rr(1:nr), fr(1:nr), qq(1:nq), fq(1:nq) 
  real(kind(1.d0)), dimension(:) :: rr, fr, qq, fq 
  integer(kind=4) nr, nq
  real(kind(1.d0)), dimension(:), allocatable:: gr
  real(kind(1.d0)):: h, theta, even, odd, q, r_1, r_n, small
  real(kind(1.d0)):: alpha, beta, gamma, snt, cst, sn2t
  real(kind(1.d0)):: one, two, three, pi
  integer ::  i, j


  small = 1.d-5
  pi = acos(-1.d0)
  nr = size(fr)
  nq = size(fq)
  h = abs(rr(2) -rr(1))

  
  allocate(gr(nr))


  gr(:) = rr(:) *  fr(:)

  r_1 = rr(1)
  r_n = rr(nr)

  do i = 1,nq
  
     q = qq(i)
     theta = h * q
     snt = sin(theta)
     sn2t = sin(2.d0*theta)
     cst = cos(theta)
     
     if(theta.lt.small) then
        alpha = 2.d0 * theta**3.d0/ 45.d0 - 2.d0 * theta**5.d0/315.d0 + 2.d0 * theta**7.d0/4725.d0
        beta  = 2.d0/3.d0  + 2.d0 * theta**2.d0/15.d0  - 4.d0 * theta**4.d0/105.d0 + 2.d0 * theta**6.d0/567
        gamma = 4.d0/3.d0  - 2.d0 * theta**2/ 15.d0 + theta**4.d0 / 210.d0 - theta**6.d0/11340.d0  
     else
        
        alpha = 1.d0 / theta + sn2t/ 2.d0/theta**2 - 2.d0 * snt*snt / theta**3
        beta  = 2.d0 *  ( (1.d0 + cst * cst)/theta**2  - sn2t/theta**3 )
        gamma = 4.d0 * ( snt / theta**3  - cst /theta**2)
     end if
     
     even = 0.d0
     odd  = 0.d0
     do j = 1, nr/2
        even = even + gr(2*j)   * sin(q*rr(2*j))
        odd  = odd  + gr(2*j-1) * sin(q*rr(2*j-1))
     end do
     even = even - 1.d0/2.d0 * ( gr(nr)* sin(q*r_n) + gr(1) * sin(q*r_1) )

     one   = alpha * (gr(1) * cos(q *r_1) - gr(nr) * cos(q * r_n) ) 
     two   = beta * even 
     three = gamma * odd
     fq(i) =   (one + two + three ) / q
  end do

  if(opt) then
     fq = 4.d0 * pi * h * fq
  else 
     fq = 4.d0 * pi * h * fq/(2.d0*pi)**3
  end if
  
  deallocate(gr)


  return
  
end subroutine FilonFourier
