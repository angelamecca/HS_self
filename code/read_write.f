!
!=======================================================================
!
      subroutine print(dove,ImSeh,ImSep,ReSEh,ReSEp,VD,ReSE1,EEh,nEh,
     &     EEp,nEp,kk,nk,kF,EkF,a,d,deg)
      implicit none
      integer nEh, nEp, nk, n, k
      double precision ImSEh(nEh,nk),ImSEp(nEp,nk)
      double precision ReSEh(nEh,nk),ReSEp(nEp,nk),VD(nEp,nk),ReSE1(nk)
      double precision EEh(nEh), EEp(nEp), kk(nk)
      double precision kF,EkF,a,d,deg,pi,delta
      character*64 filenum, nameh, namep, dove
      character*2 filenum1
      pi = acos(-1.d0)
!
!
      do k = 1, nk
         write(filenum,'(6f5.3)')  kk(k)
         nameh = 'seh'//trim(filenum)//'.dat'
         namep = 'sep'//trim(filenum)//'.dat'
         open(11,file=trim(adjustl(dove))//nameh)
         open(12,file=trim(adjustl(dove))//namep)
!
         write(11,154) '#',a,d,kF,deg
         write(12,154) '#',a,d,kF,deg

         do n = 1, nEh
            write(11,158) ImSeh(n,k),ReSeh(n,k),ReSE1(k),EEh(n),
     &           EEh(n)/EkF,kk(k),kk(k)/kF,ImSeh(n,k)/(EEh(n)-EkF)
         end do
!     
         do n = 1, nEp
            delta = ImSEp(nEp,K)*log((EEp(n)-EEp(nEp))/
     &           (EEp(n)-2.d0*EEp(nEp))) /pi
            write(12,160) ImSep(n,k),ReSep(n,k),ReSE1(k),VD(n,k),EEp(n),
     &           EEp(n)/EkF,kk(k),kk(k)/kF,ImSep(n,k)/(EEp(n)-EkF),delta
         end do
         close(11) 
         close(12)
      end do
!
      return
 154  format(a,4e18.6e3)     
 157  format(7e18.6e3)     
 158  format(8e18.6e3)     
 159  format(9e18.6e3)     
 160  format(10e18.6e3)     

!
      end subroutine print
!
!=======================================================================
!
      subroutine read(dove,ImSeh,ImSep,ReSEh,ReSEp,VD,ReSE1,EEh,nEh,
     &     EEp,nEp,kk,nk)
      implicit none
      integer nEh, nEp, nk, n, k
      double precision ImSEh(nEh,nk),ImSEp(nEp,nk)
      double precision ReSEh(nEh,nk),ReSEp(nEp,nk),VD(nEp,nk),ReSE1(nk)
      double precision EEh(nEh), EEp(nEp), kk(nk)
      double precision tmp1, tmp2, tmp3, tmp4, tmp5
      character*64 filenum, nameh, namep, dove
      character*2 filenum1
      character str
!     
      write(6,*) 'leggo', nk, nEp, nEh

      do k = 1, nk
         write(filenum,'(6f5.3)')  kk(k)
         nameh = 'seh'//trim(filenum)//'.dat'
         namep = 'sep'//trim(filenum)//'.dat'
         open(11,file=trim(adjustl(dove))//nameh)
         open(12,file=trim(adjustl(dove))//namep)
!     
         read(11,154)str,tmp1,tmp2,tmp3,tmp4
         read(12,154)str,tmp1,tmp2,tmp3,tmp4

         do n = 1, nEh
            read(11,158) ImSeh(n,k),ReSeh(n,k),ReSE1(k)
         end do
         do n = 1, nEp
            read(12,159) ImSep(n,k),ReSep(n,k),ReSE1(k),VD(n,k)
         end do
         close(11) 
         close(12)
      end do
!
      return
 154  format(a,4e18.6e3)     
 157  format(7e18.6e3)     
 158  format(8e18.6e3)     
 159  format(9e18.6e3)     

      end subroutine read
!
!=======================================================================
!
      subroutine printIM(dove,ImSe,EE,nE,kk,nk,k,kF,EkF,part)
      implicit none
      integer nE, n, nk, k
      double precision ImSE(nE,nk)
      double precision EE(nE),kk(nk)
      double precision kF,EkF,a,d,deg
      character*64 filenum, name, dove
      logical part
!
!     
      write(filenum,'(6f5.3)')  kk(k)
      if(part) then
         name = 'imp'//trim(filenum)//'.dat'
      else
         name = 'imh'//trim(filenum)//'.dat'
      end if
      open(11,file=trim(adjustl(dove))//name)
!          
      do n = 1, nE
         write(11,155) ImSE(n,k),EE(n),kk(k),kk(k)/kF,EE(n)/EkF
      end do
!     
      close(11) 
!      
      return
 154  format(a,4e18.6e3)     
 155  format(5e18.6e3)     
!
      end subroutine printIM
!
!=======================================================================
!
      subroutine readIM(dove,ImSe,EE,nE,kk,nk,k,kF,EkF,part)
      implicit none
      integer nE, n, nk, k
      double precision ImSE(nE,nk)
      double precision EE(nE),kk(nk)
      double precision kF,EkF
      character*64 filenum, name, dove
      character str
      logical part
!
!     
      write(filenum,'(6f5.3)')  kk(k)
      if(part) then
         name = 'imp'//trim(filenum)//'.dat'
      else
         name = 'imh'//trim(filenum)//'.dat'
      end if
      open(11,file=trim(adjustl(dove))//name)
!     
!     
      do n = 1, nE
         read(11,152) ImSE(n,k),EE(n)
      end do
!     
      close(11) 
!      
      return
 154  format(a,4e18.6e3)     
 152  format(5e18.6e3)     

!
      end subroutine readIM
!
!=======================================================================
!
      subroutine printMC(dove,Ris,EE,nE,kk,nk,k,kF,EkF,part,im,Ris_er,
     &     failed)
      implicit none
      integer nE, n, nk, k
      double precision Ris(nE,nk), Ris_er(nE,nk)
      double precision EE(nE),kk(nk),the_err(nE)
      integer failed(nE)
      double precision kF,EkF,a,d,deg
      character*64 filenum, name, dove
      logical part
      logical im
!
!     
      write(filenum,'(6f5.3)')  kk(k)
      if(part) then
         if(im) then
            name = 'imp'//trim(filenum)//'.dat'
         else
            name = 'rep'//trim(filenum)//'.dat'
         end if
      else
         if(im) then
            name = 'imh'//trim(filenum)//'.dat'
         else
            name = 'reh'//trim(filenum)//'.dat'
         end if
      end if
      open(11,file=trim(adjustl(dove))//name)
!          
      do n = 1, nE
         write(11,155) Ris(n,k),Ris_er(n,k),failed(n),kk(k)/kF,EE(n)/EkF
      end do
!     
      close(11) 
!      
      return
 154  format(a,4e18.6e3)     
 155  format(2e18.6e3,i5.1,2e18.6e3)     
!
      end subroutine printMC
!
!=======================================================================
!

      subroutine readMC(dove,Ris,EE,nE,kk,nk,k,part,im)
      implicit none
      integer nE, n, nk, k
      double precision Ris(nE,nk)
      double precision EE(nE),kk(nk),the_err(nE)
      integer failed(nE)
      double precision kF,EkF,a,d,deg
      character*64 filenum, name, dove
      logical part
      logical im
!
!     
      write(filenum,'(6f5.3)')  kk(k)
      if(im) then
         if(part) then
            name = 'imp'//trim(filenum)//'.dat'
         else
            name = 'imh'//trim(filenum)//'.dat'
         end if
      else
         if(part) then
            name = 'rep'//trim(filenum)//'.dat'
         else
            name = 'reh'//trim(filenum)//'.dat'
         end if
      end if
      open(11,file=trim(adjustl(dove))//name)
!          
      do n = 1, nE
         read(11,155) Ris(n,k)
      end do
!     
      close(11) 
!      
      return
 154  format(a,4e18.6e3)     
 155  format(1e18.6e3)
!
      end subroutine readMC
!
!=======================================================================
!
