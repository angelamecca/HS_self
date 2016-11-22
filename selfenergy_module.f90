
!=======================================================================
!
MODULE Boundary
  implicit none
  double precision lowerIM(10),upperIM(10)
  double precision lowerRE(10),upperRE(10)
  logical part
  logical mah
END MODULE Boundary
!
!=======================================================================
!


MODULE selfenergy
  
  implicit none
  
  double precision, dimension(:,:), allocatable:: ImSEh, ImSEh_er, ImSEp, ImSEp_er
  double precision, dimension(:,:), allocatable:: ReSEh, ReSEh_er, ReSEp, ReSEp_er
  double precision, dimension(:), allocatable:: ReSE1
  integer nEh, nEp, nk, nkF, neFh, neFp
  double precision Emax, Emin, kmax, kmin, dk
  double precision, dimension(:), allocatable::  EEh, EEp, dEh, dEp, kk
  double precision kval, Eval
  integer kprime, nprime
  
  
contains
  
  subroutine myvegas(ndim, ncomp, integrand, userdata, nvec,epsrel, epsabs, verbose, seed,&
       mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile,&
       neval, fail, integral, error, prob)
    implicit none
    integer ndim, ncomp, nvec, last, seed, mineval, maxeval
    double precision epsrel, epsabs, userdata
    integer verbose, nregions, fail
    integer nstart, nincrease, nbatch, gridno
    character*(*) statefile
    double precision integral(ncomp), error(ncomp), prob(ncomp)
    integer neval
    integer spin
    parameter (spin = -1)
    
    
    interface
       integer function integrand(ndim, xx, ncomp, ff)
         implicit none
         integer ndim, ncomp
         double precision xx(ndim), ff(ncomp)
       end function integrand
    end interface
    
    ! cuba 4.1
    call vegas(ndim, ncomp, integrand, userdata, nvec,epsrel, epsabs, verbose, seed,&
         mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, spin, &
         neval, fail, integral, error, prob)
    
    ! cuba 3.3
    !call myvegas(ndim, ncomp, integrand, userdata, nvec,epsrel, epsabs, verbose, seed,&
    !     mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, &
    !     neval, fail, integral, error, prob)

    
    return
  end subroutine myvegas
  
  
end module selfenergy


!=======================================================================

subroutine real_SE(ReSE,RESE_er,EE,nE,dE,qq,nq,dq,dove)      
  !use selfenergy, ONLY:kval, Eval, Emax
  use Boundary
  use EffectiveInteraction
  use Parameters
  use selfenergy
    
  implicit none
  !.....SET VEGAS param
  integer ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision epsrel, epsabs, userdata
  parameter (ndim = 5)
  parameter (ncomp = 1)
  parameter (userdata = 0.d0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-9)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval =  1000000000)
  
  integer nstart, nincrease, nbatch, gridno
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno = 0)
  parameter (statefile = "")
  
  
  interface
     integer function integrandREAL(ndim, xx, ncomp, ff)
       implicit none
       integer ndim, ncomp
       double precision xx(ndim), ff(ncomp)
     end function integrandREAL
  end interface
  
  double precision integral1(ncomp), error1(ncomp), prob1(ncomp)
  double precision integral2(ncomp), error2(ncomp), prob2(ncomp)
  double precision integral3(ncomp), error3(ncomp), prob3(ncomp)
  
  integer verbose, nregions, fail1, fail2, fail3
  integer neval, neval1, neval2, neval3
  character*16 env
  
  integer c
  !....END SET VEGAS
  
  double precision  pmin, pmax, const
  integer n, k, nE, nq, i
  double precision  dq, kcut, kpole
  double precision ReSE(nE,nq), ReSE_er(nE,nq),  EE(nE), dE(nE), qq(nq)
  integer  failed(nE)
  character*128 dove, filenum,nameV,nameD, filenum2
  logical im
  
  
  im      = .false.
  verbose = 0
  const   =  2.d0*pi  * M / (2.d0*pi)**6
  kcut    = 12.5d0*kF
  
  do k = 1, nq
     kval = qq(k)
     write(filenum,'(6f5.3)')  qq(k)   
     if(mah) then
        if (part) then
           write(6,*) '---------> perturbative RE Particle', qq(k), k
           nameV ='pert_repV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(nameV)//'.dat')
           open(18,file=trim(adjustl(dove))//trim(nameV)//'_1.dat')
           open(19,file=trim(adjustl(dove))//trim(nameV)//'_2.dat')
           open(20,file=trim(adjustl(dove))//trim(nameV)//'_3.dat')
        else
           write(6,*) '---------> perturbative RE Hole    ', qq(k), k
           nameV ='pert_rehV'//trim(filenum)//'.dat'
           open(17,file=trim(adjustl(dove))//nameV)            
        end if
     else
        if (part) then
           write(6,*) '---------> RE Particle', qq(k), k
           nameV ='repV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(nameV)//'.dat')
           open(18,file=trim(adjustl(dove))//trim(nameV)//'_1.dat')
           open(19,file=trim(adjustl(dove))//trim(nameV)//'_2.dat')
           open(20,file=trim(adjustl(dove))//trim(nameV)//'_3.dat')
        else
           write(6,*) '---------> RE Hole    ', qq(k), k
           nameV ='rehV'//trim(filenum)//'.dat'
           open(17,file=trim(adjustl(dove))//nameV)            
        end if
     end if
     
     do n = 1, nE
        Eval = EE(n)
        kpole = max(kF,kF+sqrt(max (0.d0,(kF)**2 + 2.d0*M*EE(n)) ))
        
        !----------------> PARTICLES  
        if(part) then
           
           lowerIM(1) =  kF
           upperIM(1) =  kpole
           lowerIM(2) =  kF
           upperIM(2) =  kpole
           lowerIM(3) = -1.d0
           upperIM(3) =  1.d0
           lowerIM(4) = -1.d0
           upperIM(4) =  1.d0
           lowerIM(5) =  0.d0
           upperIM(5) =  2.d0*pi
           
           !                  print *, "---------- Vegas RE pole pole-----------"
           !                  write(6,*) 'pole pole', upperIM(1),  upperIM(2)
           call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
                mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile,&
                neval1, fail1, integral1, error1, prob1)
           
           !                  print *, EE(n), n , EE(n)/EkF
           !                  print *, "neval    =", neval1
           !                  print *, "fail     =", fail1
           !                  print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
           !     &                 (integral1(c), error1(c), prob1(c), c = 1, ncomp)                 
           !           print *, " "
           
           write(18,155) const*integral1(1),const*error1(1),fail1,qq(k)/kF,EE(n)/EkF,neval1,kpole
           flush(18)
           
           
           lowerIM(1) =  kF
           upperIM(1) =  kpole
           lowerIM(2) =  kpole
           upperIM(2) =  kcut
           lowerIM(3) = -1.d0
           upperIM(3) =  1.d0
           lowerIM(4) = -1.d0
           upperIM(4) =  1.d0
           lowerIM(5) =  0.d0
           upperIM(5) =  2.d0*pi
           
           !                  print *, "---------- Vegas RE pole cut-----------"
           !                  write(6,*) 'pole cut',  upperIM(1),  upperIM(2)
           
           call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
                mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile, &
                neval2, fail2, integral2, error2, prob2)
           
           !                  print *, EE(n), n , EE(n)/EkF
           !                  print *, "neval    =", neval2
           !                  print *, "fail     =", fail2
           !                  print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
           !     &                 (integral2(c), error2(c), prob2(c), c = 1, ncomp)
           
           write(19,155) const*integral2(1),const*error2(1),fail2,qq(k)/kF,EE(n)/EkF,neval2,kpole
           flush(19)
           
           !           print *, " "
           
           lowerIM(1) =  kpole
           upperIM(1) =  kcut
           lowerIM(2) =  kpole
           upperIM(2) =  kcut
           lowerIM(3) = -1.d0
           upperIM(3) =  1.d0
           lowerIM(4) = -1.d0
           upperIM(4) =  1.d0
           lowerIM(5) =  0.d0
           upperIM(5) =  2.d0*pi
           
           !                  print *, "---------- Vegas RE cut cut -----------"
           !                  write(6,*) 'cut cut',  upperIM(1),  upperIM(2)
           
           call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
                mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile, &
                neval3, fail3, integral3, error3, prob3)
           
           !                  print *, EE(n), n , EE(n)/EkF
           !                  print *, "neval    =", neval3
           !                  print *, "fail     =", fail3
           !                  print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
           !     &                 (integral3(c), error3(c), prob3(c), c = 1, ncomp)
           
           !           print *, " "
           write(20,155) const*integral3(1),const*error3(1),fail3,qq(k)/kF,EE(n)/EkF,neval3,kpole
           flush(20)
           !                                                     
           ReSE(n, k) = const*(integral1(1)+2.d0*integral2(1)+integral3(1))
           ReSE_er(n, k) = const*sqrt(error1(1)**2+error2(1)**2+error3(1)**2)
           Failed(n) = fail1 + fail2 + fail3
           neval     = neval1+ neval2 + neval3
           
           write(6,156) ReSE(n, k),ReSE_er(n,k),Failed(n),qq(k)/kF,EE(n)/EkF,neval
           
           write(17,156) ReSE(n, k),ReSE_er(n,k),Failed(n),qq(k)/kF,EE(n)/EkF,neval
           flush(17)     
           
           !----------------> HOLES  
        else
           
           lowerIM(1) =  0.d0
           upperIM(1) =  kF
           lowerIM(2) =  0.d0
           upperIM(2) =  kF
           lowerIM(3) = -1.d0
           upperIM(3) =  1.d0
           lowerIM(4) = -1.d0
           upperIM(4) =  1.d0
           lowerIM(5) =  0.d0
           upperIM(5) =  2.d0*pi
           
           !                  print *, "---------- Vegas RE hole-----------"
           call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
                mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile,&
                neval, fail1, integral1, error1, prob1)
           
           !                  print *, EE(n), n , EE(n)/EkF
           !                  print *, "neval    =", neval
           !                  print *, "fail     =", fail1
           !                  print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
           !     &                 (integral1(c), error1(c), prob1(c), c = 1, ncomp)
           
           !           print *, " "
           ReSE(n, k) = const*integral1(1)
           ReSE_er(n,k)= const*error1(1)
           Failed(n)= fail1
           neval = neval1
           write(17,156) ReSE(n, k),ReSE_er(n,k),Failed(n),qq(k)/kF,EE(n)/EkF,neval
           flush(17)
        end if
        
     end do
     close(17)
     close(18)
     close(19) 
     close(20)
     !     
     call printMC(dove,ReSe,EE,nE,qq,nq,k,kF,EkF,part,im,ReSE_er,failed)
  end do
  !the_err
  return
  
155 format(2e18.8e3,i12.1,2e18.8e3,i12.1,e18.8e3) 
156 format(2e18.8e3,i12.1,2e18.8e3,i12.1)     
  
end subroutine real_SE

!===============================================================================



!=======================================================================

subroutine real_SE_ek(ReSE,the_err,Ene,qq,dove)      
  !use selfenergy, ONLY:kval, Eval, Emax
  use Boundary
  use EffectiveInteraction
  use Parameters
  use selfenergy
  
  implicit none
  
  !.....SET VEGAS param
  integer ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision epsrel, epsabs, userdata
  parameter (ndim = 5)
  parameter (ncomp = 1)
  parameter (userdata = 0.d0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-9)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval =  1000000000)
  
  integer nstart, nincrease, nbatch, gridno
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno = 0)
  parameter (statefile = "")
  
  interface
     integer function integrandREAL(ndim, xx, ncomp, ff)
       implicit none
       integer ndim, ncomp
       double precision xx(ndim), ff(ncomp)
     end function integrandREAL
  end interface
  
  double precision integral1(ncomp), error1(ncomp), prob1(ncomp)
  double precision integral2(ncomp), error2(ncomp), prob2(ncomp)
  double precision integral3(ncomp), error3(ncomp), prob3(ncomp)
  
  integer verbose, nregions, fail1, fail2, fail3
  integer neval, neval1, neval2, neval3
  character*16 env
  
  integer c
  !....END SET VEGAS
  
  double precision  pmin, pmax, const
  integer n, k, nE,  i
  double precision   kcut, kpole
  double precision ReSE, Ene, qq, the_err
  integer  failed
  character*128 dove, filenum,nameV,nameD, filenum2
  logical im
  
  im      = .false.
  verbose = 0
  const   = 2.d0*pi  * M / (2.d0*pi)**6
  kcut    = 12.5d0*kF
  
  kval = qq
  Eval = Ene
  
  
  !----------------> PARTICLES  
  if(part) then
     kpole = max(kF, kF+ sqrt(max (0.d0, (kF)**2 + 2.d0*M*Eval)))
     !write(6,*) 'kcut kpole', kcut, kpole
     
     
     lowerIM(1) =  kF
     upperIM(1) =  kpole
     lowerIM(2) =  kF
     upperIM(2) =  kpole
     lowerIM(3) = -1.d0
     upperIM(3) =  1.d0
     lowerIM(4) = -1.d0
     upperIM(4) =  1.d0
     lowerIM(5) =  0.d0
     upperIM(5) =  2.d0*pi
     
     !            print *, "---------- Vegas RE pole pole-----------"
     !            write(6,*) 'pole pole', upperIM(1),  upperIM(2)
     call myvegas(ndim, ncomp, integrandREAL, userdata, nvec, &
          epsrel, epsabs, verbose, seed,mineval, maxeval, nstart, nincrease, nbatch, &
          gridno, statefile, neval1, fail1, integral1, error1, prob1)
     
     !            print *, Ene, Ene/EkF
     !            print *, "neval    =", neval1
     !            print *, "fail     =", fail1
     
     !            print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
     !     &           (integral1(c), error1(c), prob1(c), c = 1, ncomp)
     !            write(6,155) const*integral1(1),const*error1(1),fail1
     !     &           ,qq/kF,Ene/EkF,neval1,kpole
     !     print *, " "
     
     lowerIM(1) =  kF
     upperIM(1) =  kpole
     lowerIM(2) =  kpole
     upperIM(2) =  kcut
     lowerIM(3) = -1.d0
     upperIM(3) =  1.d0
     lowerIM(4) = -1.d0
     upperIM(4) =  1.d0
     lowerIM(5) =  0.d0
     upperIM(5) =  2.d0*pi
     !     
     !            print *, "---------- Vegas RE pole cut-----------"
     !            write(6,*) 'pole cut',  upperIM(1),  upperIM(2)
     
     call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
          mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile,&
          neval2, fail2, integral2, error2, prob2)
     !     
     !            print *, Ene, Ene/EkF
     !            print *, "neval    =", neval2
     !            print *, "fail     =", fail2
     
       !            print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
     !     &           (integral2(c), error2(c), prob2(c), c = 1, ncomp)
     !            write(6,155) const*integral2(1),const*error2(1),fail2
     !     &           ,qq/kF,Ene/EkF,neval2,kpole
     !            
     !     print *, " "
     !     
     lowerIM(1) =  kpole
     upperIM(1) =  kcut
     lowerIM(2) =  kpole
     upperIM(2) =  kcut
     lowerIM(3) = -1.d0
     upperIM(3) =  1.d0
     lowerIM(4) = -1.d0
     upperIM(4) =  1.d0
     lowerIM(5) =  0.d0
     upperIM(5) =  2.d0*pi
     !     
     !            print *, "---------- Vegas RE cut cut -----------"
     !            write(6,*) 'cut cut',  upperIM(1),  upperIM(2)
     
     call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
          mineval, maxeval, nstart, nincrease, nbatch,gridno, statefile,&
          neval3, fail3, integral3, error3, prob3)
     !     
     !            print *, Ene,Ene/EkF
     !            print *, "neval    =", neval3
     !            print *, "fail     =", fail3
     
     !           print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
     !     &           (integral3(c), error3(c), prob3(c), c = 1, ncomp)            
     !            write(6,155) const*integral3(1),const*error3(1),fail3
     !     &           ,qq/kF,Ene/EkF,neval3,kpole
     !     print *, " "
     !     
     
     ReSE = const*(integral1(1)+2.d0*integral2(1)+integral3(1))
     the_err = const*sqrt(error1(1)**2 + 4.d0*error2(1)**2 + error3(1)**2)
     Failed = fail1 + fail2 + fail3
     neval     = neval1+ neval2 + neval3
     write(77,*) ReSE,the_err,Failed, Ene,Ene/EkF,neval
     flush(77)
     !     
     
     !
     !----------------> HOLES  
  else
     !     
     lowerIM(1) =  0.d0
     upperIM(1) =  kF
     lowerIM(2) =  0.d0
     upperIM(2) =  kF
     lowerIM(3) = -1.d0
     upperIM(3) =  1.d0
     lowerIM(4) = -1.d0
     upperIM(4) =  1.d0
     lowerIM(5) =  0.d0
     upperIM(5) =  2.d0*pi
     !     
     !            print *, "---------- Vegas RE hole-----------"
     call myvegas(ndim, ncomp, integrandREAL, userdata, nvec,epsrel, epsabs, verbose, seed,&
          mineval, maxeval, nstart, nincrease, nbatch, gridno, statefile, &
          neval, fail1, integral1, error1, prob1)
     !     
     !            print *, Ene, Ene/EkF
     !            print *, "neval    =", neval
     !            print *, "fail     =", fail1
     !            print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
     !     &           (integral1(c), error1(c), prob1(c), c = 1, ncomp)
     !            write(6,*) const*integral2(1),const*error2(1),fail2
     !     &           ,qq/kF,Ene/EkF,neval2,kpole
     
     !     print *, " "
     ReSE = const*integral1(1)
     the_err= const*error1(1)
     Failed= fail1
     neval = neval1
     write(78,*) ReSE,the_err,Failed, Ene,qq,neval
     flush(78)
     
     
    end if
    !     
    close(17)
    !     close(18)
    !     close(19) 
    !     close(20)
    !     
    
    return
    
155 format(2e18.8e3,i12.1,2e18.8e3,i12.1,e18.8e3) 
156 format(2e18.12e3,i12.1,2e18.12e3,i12.1)     
    
  end subroutine real_SE_ek
  
!===============================================================================
  
  
  
  
!===============================================================================
  
  integer function integrandREAL(ndim, xx, ncomp, ff)
    use Boundary
    use Parameters
    
    implicit none
    integer i,ndim, ncomp
    double precision xx(ndim), ff(ncomp)
    double precision x(ndim),range(ndim)
    double precision jacobian, ffint_realp, ffint_realh
    double precision ffint_rp_pert, ffint_rh_pert
    !     
    !
    jacobian = 1.
    do i = 1,ndim
       range(i) = upperIM(i) - lowerIM(i)
       jacobian = jacobian*range(i) 
       x(i) = lowerIM(i) + range(i)*xx(i)
    end do
    !  
    if(mah) then
       if(part) then
          ff(1) = jacobian*ffint_rp_pert(x(1),x(2),x(3),x(4),x(5))   
       else
          ff(1) = jacobian*ffint_rh_pert(x(1),x(2),x(3),x(4),x(5))   
       end if
    else
       if(part) then
          ff(1) = jacobian*ffint_realp(x(1),x(2),x(3),x(4),x(5))   
       else
          ff(1) = jacobian*ffint_realh(x(1),x(2),x(3),x(4),x(5))   
       end if
    end if
    !     
    !     
    integrandREAL = 0
    !     
    return
  end function integrandREAL
  
!===============================================================================
!  
!     Real part Sigma HOLES
!
!===============================================================================

  double precision function ffint_realh(q1,q2,cq1,cq2,phi)
    use selfenergy
    use EffectiveInteraction
    use Parameters

    implicit none
    double precision q1, q2, cq1, cq2, phi
    double precision small, fint
    double precision ff1, ff3, M1, M2, Muv, u, v
    double precision sq1, sq2, k1, k2
    double precision den
    
    ffint_realh = 0.d0
    small = 1.e-3
    k1 = kval
         
    sq1 = sqrt(1.d0 - cq1**2)
    sq2 = sqrt(1.d0 - cq2**2)   
    
    k2 =   sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*(q1*cq1 + q2*cq2 ) &
         + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
    
    if (k2.lt.kF) return
    den = 2.d0*m*Eval -  q1**2 - q2**2 + k2**2 
    ff1 = den /  (den**2 + small**2)
    !write(6,*), 'realh', Eval, q1, q2, k2, den**2, small**2
    !write(6,*), 'realh', den**2, small**2
    
    ff3 = q1**2 * q2**2
    
    u = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
    v = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
    M1 = fint(p,Mp,np,dp,u,1)
    M2 = fint(p,Mp,np,dp,v,1)      
    
    Muv = deg*(M1**2 + M2**2) - 2.d0*M1*M2      
         
    ffint_realh = ff1*ff3*Muv
    return
  end  function ffint_realh
  
!===============================================================================
!     
!     Real part Sigma PARTICLES
!     
!===============================================================================
  
  double precision function ffint_realp(q1,q2,cq1,cq2,phi)
    use selfenergy
    use EffectiveInteraction
    use Parameters

    implicit none
    double precision q1, q2, cq1, cq2, phi
    double precision small, fint
    double precision ff1, ff3, M1, M2, Muv, u, v
    double precision sq1, sq2, k1, k2
    double precision den
    
    ffint_realp = 0.d0
    small = 1.e-3
    k1 = kval
    !      write(6,*) ' kk nella funzione', kval, k1
    !     
    sq1 = sqrt(1.d0 - cq1**2)
    sq2 = sqrt(1.d0 - cq2**2)   
    !     
    k2 =   sqrt(q1**2 + q2**2 + k1**2  - 2.d0*k1*(q1*cq1 + q2*cq2) &
         + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
    !     
    if (k2.gt.kF) return
    den =  -2.d0*m*Eval +  q1**2 + q2**2 - k2**2 
    ff1 = den /  (den**2 + small**2)
    !write(6,*), 'realp', Eval, q1, q2, k2, den**2, small**2
    !write(6,*), 'realp',  den**2, small**2
    
    
    ff3 = q1**2 * q2**2
    
    u = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
    v = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
    M1 = fint(p,Mp,np,dp,u,1)
    M2 = fint(p,Mp,np,dp,v,1)      
    
    Muv = deg*(M1**2 + M2**2) - 2.d0*M1*M2      
    
    ffint_realp = - ff1*ff3*Muv
    return
  end  function ffint_realp
  
!===============================================================================
  
  
  
  
  
  


!===============================================================================
!
!     PERTURBATIVE Real part Sigma HOLES
!
!===============================================================================
!
  double precision function ffint_rh_pert(q1,q2,cq1,cq2,phi)
    use selfenergy
    use EffectiveInteraction
    use Parameters

    implicit none
    double precision q1, q2, cq1, cq2, phi
    double precision small, fint
    double precision ff1, ff3, M1, M2, Muv, u, v
    double precision sq1, sq2, k1, k2
    double precision den
    !
    ffint_rh_pert = 0.d0
    small = 1.e-3
    k1 = kval
    !     
    sq1 = sqrt(1.d0 - cq1**2)
    sq2 = sqrt(1.d0 - cq2**2)   
    !     
    k2 =   sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*(q1*cq1 + q2*cq2 ) &
         + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
    !     
    if (k2.lt.kF) return
    den = 2.d0*m*Eval -  q1**2 - q2**2 + k2**2 
    ff1 = den /  (den**2 + small**2)
    !     
    ff3 = q1**2 * q2**2
    !     
    Muv = 2.d0 *(deg-1.d0)*(4.d0*pi*a*htm/m)**2
    !     
    ffint_rh_pert = ff1*ff3*Muv
    return
  end  function ffint_rh_pert
  !     
!===============================================================================
!     
!     PERTURBATIVE Real part Sigma PARTICLES
!     
!===============================================================================
!     
  double precision function ffint_rp_pert(q1,q2,cq1,cq2,phi)
    use selfenergy
    use EffectiveInteraction
    use Parameters

    implicit none
    double precision q1, q2, cq1, cq2, phi
    double precision small, fint
    double precision ff1, ff3, M1, M2, Muv, u, v
    double precision sq1, sq2, k1, k2
    double precision den
    
    ffint_rp_pert = 0.d0
    small = 1.e-3
    k1 = kval
    
    sq1 = sqrt(1.d0 - cq1**2)
    sq2 = sqrt(1.d0 - cq2**2)   
    
    k2 =   sqrt(q1**2 + q2**2 + k1**2  - 2.d0*k1*(q1*cq1 + q2*cq2) &
         + 2.d0*q1*q2*(sq1*sq2*sin(phi) + cq1*cq2) )
    
    if (k2.gt.kF) return
    den =  -2.d0*m*Eval +  q1**2 + q2**2 - k2**2 
    ff1 = den /  (den**2 + small**2)
    
      
    ff3 = q1**2 * q2**2
    
    Muv = 2.d0 *(deg-1.d0)*(4.d0*pi*a*htm/m)**2
    
    ffint_rp_pert = - ff1*ff3*Muv
    return
  end  function ffint_rp_pert
     
!===============================================================================











!
!=======================================================================
!
  subroutine imaginary_SE(ImSE,ImSE_er,EE,nE,dE,qq,nq,dq,dove)      
    !use selfenergy, ONLY:kval, Eval
    use Boundary
    use EffectiveInteraction
    use Parameters
    use selfenergy
    implicit none
    !
    !.....SET VEGAS param
    integer ndim, ncomp, nvec, last, seed, mineval, maxeval
    double precision epsrel, epsabs, userdata
    parameter (ndim = 4)
    parameter (ncomp = 1)
    parameter (userdata = 0.d0)
    parameter (nvec = 1)
    parameter (epsrel = 1D-3)
    parameter (epsabs = 1D-9)
    parameter (last = 4)
    parameter (seed = 0)
    parameter (mineval = 0)
    parameter (maxeval = 1000000000)
    
  integer nstart, nincrease, nbatch, gridno, gridno1, gridno2
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno1 = 10)
  parameter (gridno2 = 0)
  parameter (statefile = "")
  
    interface
       integer function integrandIM(ndim, xx, ncomp, ff)
         implicit none
         integer ndim, ncomp
         double precision xx(ndim), ff(ncomp)
         
       end function integrandIM
    end interface
  
  double precision integral(ncomp), error(ncomp), prob(ncomp)
  integer verbose, nregions, neval, fail
  character*16 env
  integer c
!....END SET VEGAS
  
  double precision plim_csieta, pp_csieta, pmin, pmax
  integer n, k, nE, nq
  double precision  dq
  double precision ImSE(nE,nq),ImSE_er(nE,nq), EE(nE), dE(nE), qq(nq), the_err(nE)
  double precision const
  integer  failed(nE)
  character*128 dove, name, filenum
  logical im

  im = .true.
  verbose = 0
  const =  2.d0*pi * pi * M / (2.d0*pi)**6

  open(7,file='fail.dat')            
  do k = 1, nq
     kval = qq(k)

     write(filenum,'(6f5.3)')  qq(k)   
     if(mah) then
        if (part) then
           write(6,*) '---------> perturbative Im Particle', qq(k), k
           name ='pert_impV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(name)//'.dat')
        else
           write(6,*) '---------> perturbative Im Hole    ', qq(k), k
           name ='pert_imhV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(name)//'.dat')            
        end if
     else
        if (part) then
           write(6,*) '---------> Im Particle', qq(k), k
           name ='impV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(name)//'.dat')
        else
           write(6,*) '---------> Im Hole    ', qq(k), k
           name ='imhV'//trim(filenum)
           open(17,file=trim(adjustl(dove))//trim(name)//'.dat')            
        end if
     end if
     
     do n = 1, nE
        Eval = EE(n)
        
        pp_csieta   = sqrt(2.d0*M*EE(n))    
                  
        if(part) then
           !               write(6,*) '---------> Particle', qq(k), EE(n), k, n
           pmin = kF
           pmax = max(kF,pp_csieta)
           !pmax = pp_csieta
           gridno = gridno1
        else
           !               write(6,*) '---------> Hole', qq(k), EE(n), k, n
           pmin = 0.d0
           pmax = kF
           gridno = gridno2
        end if
        !     
        !     
        lowerIM(1) =  pmin
        upperIM(1) =  pmax
        lowerIM(2) =  pmin
        upperIM(2) =  pmax
        lowerIM(3) = -1.d0
        upperIM(3) =  1.d0
        lowerIM(4) = -1.d0
        upperIM(4) =  1.d0
!          
        call myvegas(ndim, ncomp, integrandIM, userdata, nvec,  &
             epsrel, epsabs, verbose, seed, mineval, maxeval, &
             nstart, nincrease, nbatch,gridno, statefile,     &
             neval, fail, integral, error, prob)
        !     
        !     print *, EE(n), n 
        !     print *, "neval    =", neval
        !     print *, "fail     =", fail
        !     print '(F25.12," +- ",F25.12,"   p = ",F8.3)',
        !     &           (integral(c), error(c), prob(c), c = 1, ncomp)
        
        !     print *, " "
        IMSE(n, k)    = const*integral(1)
        IMSE_er(n, k) = const*error(1)
        Failed(n)     = fail

        write(17,156) ImSE(n, k),ImSE_er(n,k),Failed(n),qq(k)/kF,EE(n)/EkF,neval
        flush(17)     
     end do
     close(17)
     
     call printMC(dove,ImSe,EE,nE,qq,nq,k,kF,EkF,part,im,ImSE_er,failed)
     
  end do
  !  

156 format(2e18.8e3,i12.1,2e18.8e3,i12.1)     
  
  return
end subroutine imaginary_SE
!     
!===============================================================================
!     

!
!===============================================================================
!
subroutine im_SE_ek(ImSE,the_err, Ene,qq,dove)  
  use Boundary
  use EffectiveInteraction
  use Parameters
  use selfenergy
  
  !use selfenergy, ONLY:kval, Eval
  implicit none
  !.....SET VEGAS param
  
  integer ndim, ncomp, nvec, last, seed, mineval, maxeval
  double precision epsrel, epsabs, userdata
  parameter (ndim = 4)
  parameter (ncomp = 1)
  parameter (userdata = 0.d0)
  parameter (nvec = 1)
  parameter (epsrel = 1D-3)
  parameter (epsabs = 1D-9)
  parameter (last = 4)
  parameter (seed = 0)
  parameter (mineval = 0)
  parameter (maxeval =  1000000000)
  
  integer nstart, nincrease, nbatch, gridno
  character*(*) statefile
  parameter (nstart = 50000)
  parameter (nincrease = 500)
  parameter (nbatch = 1000)
  parameter (gridno = 0)
  parameter (statefile = "")
  
    interface
       integer function integrandIM(ndim, xx, ncomp, ff)
         implicit none
         integer ndim, ncomp
         double precision xx(ndim), ff(ncomp)
         
       end function integrandIM
    end interface
  
  double precision integral(ncomp), error(ncomp), prob(ncomp)
  integer verbose, nregions, neval, fail
  character*16 env
  integer c
  !....END SET VEGAS
  
  double precision plim_csieta, pp_csieta, pmin, pmax
  double precision ImSE,the_err, Ene, qq
  double precision const
  character*128 dove, name


  verbose = 0
  const =  2.d0*pi * pi * M / (2.d0*pi)**6
  kval = qq
  Eval = Ene
  pp_csieta   = sqrt(2.d0*M*Ene)    
!     
!     
  if(part) then
     pmin = kF
     pmax = max(kF,pp_csieta)
     write(6,*) '---------> Im Particle E k', qq, Ene, pmin, pmax

  else
     pmin = 0.d0
     pmax = kF
     write(6,*) '---------> Im  Hole E k ', qq, Ene, pmin, pmax

  end if
  
  lowerIM(1) =  pmin
  upperIM(1) =  pmax
  lowerIM(2) =  pmin
  upperIM(2) =  pmax
  lowerIM(3) = -1.d0
  upperIM(3) =  1.d0
  lowerIM(4) = -1.d0
  upperIM(4) =  1.d0
  
  call myvegas(ndim, ncomp, integrandIM, userdata, nvec,&
       epsrel, epsabs, verbose, seed, mineval, maxeval, &
       nstart, nincrease, nbatch,gridno, statefile,neval,&
       fail, integral, error, prob)
  
  !print *, "neval    =", neval
  !print *, "fail     =", fail
  !print '(F25.12," +- ",F25.12,"   p = ",F8.3)', (const*integral(c), const*error(c), prob(c), c = 1, ncomp)
  
  !print *, " "
  IMSE = const* integral(1)
  the_err= const*error(1)
  
  !     
  !     
  return
end subroutine im_SE_ek
!     
!===============================================================================
!
!     
!===============================================================================
!
integer function integrandIM(ndim, xx, ncomp, ff)
  use Boundary
  use Parameters

  implicit none
  integer i,ndim, ncomp
  double precision xx(*), ff(*)
  double precision x(ndim),range(10)
  double precision jacobian, ffint_imh,ffint_imp 
  double precision ffint_imh_pert,ffint_imp_pert
  
  jacobian = 1.
  do i = 1,ndim
     range(i) = upperIM(i) - lowerIM(i)
     jacobian = jacobian*range(i) 
     x(i) = lowerIM(i) + range(i)*xx(i)
  end do
  !    
  if(mah) then
     if(part) then
        ff(1) = jacobian*ffint_imp_pert(x(1),x(2),x(3),x(4))        
     else
        ff(1) = jacobian*ffint_imh_pert(x(1),x(2),x(3),x(4))        
     end if

  else
     if(part) then
        ff(1) = jacobian*ffint_imp(x(1),x(2),x(3),x(4))        
     else
        ff(1) = jacobian*ffint_imh(x(1),x(2),x(3),x(4))        
     end if
  end if

  integrandIM = 0
  
  return
end function integrandIM

!===============================================================================
!
!     Imaginary part Sigma HOLES
!
!===============================================================================


double precision function ffint_imh(q1,q2,cq1,cq2)
  use selfenergy
  use EffectiveInteraction
  use Parameters

  implicit none
  double precision q1, q2, cq1, cq2
  double precision k1, k2, sq1, sq2, phi, fint
  double precision num, den, sphi0, cphi0
  double precision ff1, ff3,  M1, M2, u, v, Muv
  
  ffint_imh = 0.d0
  k1 = kval   
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)
  num = -2.d0*M*Eval- k1**2 +  2.d0*k1*(q1*cq1 + q2*cq2)
  den =  2.d0*q1*sq1*q2*sq2
  sphi0 = num/den - (cq1*cq2)/(sq1 *sq2 )
  if (abs(sphi0).gt.1.d0) return
  k2 = sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ) )     
  if(k2.lt.kF) return
  
  cphi0 = sqrt(1.d0-sphi0**2)
  ff1 =  2.d0/cphi0/den
  ff3 = q1**2 * q2**2
  u = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  v = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
  M1 = fint(p,Mp,np,dp,u,1)
  M2 = fint(p,Mp,np,dp,v,1)           
  Muv = deg*(M1**2 + M2**2) - 2.d0*M1*M2 
  
  !      const =  2.d0*pi * pi * M / (2.d0*pi)**6
  
  ffint_imh =  ff1*ff3*Muv
  return
end  function ffint_imh
!     
!===============================================================================
!     





!===============================================================================
!
!     Imaginary part Sigma PARTICLES
!
!===============================================================================


double precision function ffint_imp(q1,q2,cq1,cq2)
  use selfenergy
  use EffectiveInteraction
  use Parameters

  implicit none
  double precision q1, q2, cq1, cq2
  double precision k1, k2, sq1, sq2, phi,small, fint
  double precision num, den, sphi0, cphi0
  double precision ff1, ff3,  M1, M2, u, v, Muv
  
  ffint_imp = 0.d0
  small = 1.e-8  
  k1 = kval       
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)
  
  num = -2.d0*M*Eval- k1**2 +  2.d0*k1*(q1*cq1 + q2*cq2)
  den =  2.d0*q1*sq1*q2*sq2
  sphi0 = num/den - cq1/sq1* cq2/sq2
  if (abs(sphi0).gt.1.) return
  k2 = sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ))      
  if(k2.gt.kF) return
  cphi0 = sqrt(1.d0-sphi0**2)
  ff1 =  2.d0/(cphi0*den)
  ff3 = q1**2 * q2**2
  u = sqrt( q1**2 + k1**2 - 2.d0*k1*q1*cq1)
  v = sqrt( q2**2 + k1**2 - 2.d0*k1*q2*cq2)
  M1 = fint(p,Mp,np,dp,u,1)
  M2 = fint(p,Mp,np,dp,v,1)           
  Muv = deg*(M1**2 + M2**2) - 2.d0*M1*M2 
  
  !      const =  2.d0*pi * pi * M / (2.d0*pi)**6
  
  ffint_imp =  ff1*ff3*Muv
  return
end  function ffint_imp
!     
!===============================================================================
!     


!===============================================================================
!
!     Perturbative Imaginary part Sigma HOLES
!
!===============================================================================


double precision function ffint_imh_pert(q1,q2,cq1,cq2)
  use selfenergy
  use EffectiveInteraction
  use Parameters

  implicit none
  double precision q1, q2, cq1, cq2
  double precision k1, k2, sq1, sq2, phi, fint
  double precision num, den, sphi0, cphi0
  double precision ff1, ff3,  M1, M2, u, v, Muv
  
  ffint_imh_pert = 0.d0
  k1 = kval   
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)
  num = -2.d0*M*Eval- k1**2 +  2.d0*k1*(q1*cq1 + q2*cq2)
  den =  2.d0*q1*sq1*q2*sq2
  sphi0 = num/den - (cq1*cq2)/(sq1 *sq2 )
  if (abs(sphi0).gt.1.d0) return
  k2 = sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ) )     
  if(k2.lt.kF) return
  
  cphi0 = sqrt(1.d0-sphi0**2)
  ff1 =  2.d0/cphi0/den
  ff3 = q1**2 * q2**2
  Muv = 2.d0 *(deg-1.d0)*(4.d0*pi*a*htm/m)**2

  
  !      const =  2.d0*pi * pi * M / (2.d0*pi)**6
  
  ffint_imh_pert =  ff1*ff3*Muv
  return
end  function ffint_imh_pert
!     
!===============================================================================
!     





!===============================================================================
!
!     Perturbative Imaginary part Sigma PARTICLES
!
!===============================================================================


double precision function ffint_imp_pert(q1,q2,cq1,cq2)
  use selfenergy
  use EffectiveInteraction
  use Parameters

  implicit none
  double precision q1, q2, cq1, cq2
  double precision k1, k2, sq1, sq2, phi,small, fint
  double precision num, den, sphi0, cphi0
  double precision ff1, ff3,  M1, M2, u, v, Muv
  
  ffint_imp_pert = 0.d0
  small = 1.e-8  
  k1 = kval       
  
  sq1 = sqrt(1.d0 - cq1**2)
  sq2 = sqrt(1.d0 - cq2**2)
  
  num = -2.d0*M*Eval- k1**2 +  2.d0*k1*(q1*cq1 + q2*cq2)
  den =  2.d0*q1*sq1*q2*sq2
  sphi0 = num/den - cq1/sq1* cq2/sq2
  if (abs(sphi0).gt.1.) return
  k2 = sqrt(q1**2 + q2**2 + k1**2 - 2.d0*k1*( q1 * cq1 + q2 * cq2 ) &
       + 2.d0*q1*q2*( sq1*sq2*sphi0 + cq1*cq2 ))      
  if(k2.gt.kF) return
  cphi0 = sqrt(1.d0-sphi0**2)
  ff1 =  2.d0/(cphi0*den)
  ff3 = q1**2 * q2**2
  Muv = 2.d0 *(deg-1.d0)*(4.d0*pi*a*htm/m)**2

  
  !      const =  2.d0*pi * pi * M / (2.d0*pi)**6
  
  ffint_imp_pert =  ff1*ff3*Muv
  return
end  function ffint_imp_pert
!     
!===============================================================================
!     



!=======================================================================
!
!    Sigma REAL PART 1 order 
!
!=======================================================================
!
subroutine sigma_hf(r, nr, dr, ell, vr, ki, rho, deg, ReSE)
  use Boundary, ONLY:mah
  use Parameters, ONLY: pi, a, htm, m
  
  implicit none
  integer nr, i
  double precision r(nr), vr(nr), ell(nr)
  double precision dr, ki, rho, deg, ReSE
  double precision I0, I1, y
  
  if (mah) then
     ReSE = rho * (1.d0 -1.d0/deg) * (4.d0*pi*a*htm/m)  
  else
     I0=0.d0
     I1=0.d0
     if (ki.eq.0.) then
        do i=1, nr
           I0 = I0 + r(i)*r(i)*vr(i)
           I1 = I1 + r(i) *r(i) * ell(i)*vr(i)
        end do
     else
        do i=1, nr
           I0 = I0 + r(i)*r(i)*vr(i)
           I1 = I1 + r(i) * ell(i)*vr(i)* sin(ki*r(i))/ki
        end do
     end if
     
     ReSE = 4.d0*pi*dr*rho* (I0 - I1/deg )
  end if
  return
end subroutine sigma_hf
!     
!=======================================================================
!




  
  
