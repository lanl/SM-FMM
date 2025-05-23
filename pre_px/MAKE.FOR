      parameter(npts=64)
      dimension phg(10000),thg(10000),omg(10000)
      dimension indaux(npts,npts,npts)
c
      open(10,file='make.in',status='old')
      open(12,file='table.dat',status='old')
      open(11,file='microstr.dat',status='unknown')
c
      read(10,*) ngrain
      do i=1,ngrain
      read(10,*) phg(i),thg(i),omg(i)
      enddo
c
      read(12,*) npts1,npts2,npts3,ng
c
      if(ng.ne.ngrain) then
      write(*,*) 'DIFF # OF ORIENTATIONS IN MAKE.IN & TABLE'
      stop
      endif
c
      do k=1,npts3
cver      read(12,'(a)') dum
      read(12,*)
      do i=1,npts1
      read(12,*) (indaux(i,j,k),j=1,npts2)
      enddo
      enddo

      ig=0
      do k=1,npts3
      do j=1,npts2
      do i=1,npts1
      ig=ig+1
c
      iph=1
cw      if(indaux(i,j,k).lt.0) iph=2
c
      write(11,'(3(F9.2,1x),1x,5(i4,1x))') 
     #      phg(abs(indaux(i,j,k))),thg(abs(indaux(i,j,k))),
     #      omg(abs(indaux(i,j,k))),i,j,k,abs(indaux(i,j,k)),iph
c
      enddo
      enddo
      enddo

      end
