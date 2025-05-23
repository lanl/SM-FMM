      parameter (nptsmax=128)
      parameter (ngrmax=5000)
c
      dimension kgrain(0:nptsmax+1,0:nptsmax+1,0:nptsmax+1)
      dimension x(27*ngrmax),y(27*ngrmax),z(27*ngrmax),kgr(27*ngrmax)
      dimension aspect(3)
cz
      open(15,file='voro.in',status='old')
cz
      open(16,file='table.dat',status='unknown')
c
cz      write(*,*) ' HOW MANY GRAINS -> '
cz      read(*,*) ngr
      read(15,*) ngr
cz
cc      write(5,*)'NGR  = ',ngr
cz      write(*,*) ' HOW MANY POINTS X,Y,Z -> '
cz      read(*,*) npts1,npts2,npts3
      read(15,*) npts1,npts2,npts3
      read(15,*) aspect(1),aspect(2),aspect(3)
cc      write(5,*)'NPTS = ',npts
c
cz      write(6,*) npts1,npts2,npts3,ngr
      write(16,*) npts1,npts2,npts3,ngr
c
cz      write(*,*) ' RANDOM SEED FOR GRAIN GENERATION  -> '
cz      read(*,*) idum
      read(15,*) idum
cc      write(5,*)'SEED FOR GRAIN GENERATION  = ',idum
c
c     discard first random number
c
      dum=ran(idum)
c
      write(*,*) 'GENERATING SET OF POINTS'
c
      do i=1,ngr
      kgr(i)=i
      x(i)=ran2(idum)*npts1
      y(i)=ran2(idum)*npts2
      z(i)=ran2(idum)*npts3
      enddo
c
      jgr=ngr
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)
      z(jgr)=z(i)-float(npts3)
      enddo
c
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)+float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)+float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)+float(npts3)
      enddo
c
      do i=1,ngr
      jgr=jgr+1
      kgr(jgr)=i
      x(jgr)=x(i)-float(npts1)
      y(jgr)=y(i)-float(npts2)
      z(jgr)=z(i)-float(npts3)
      enddo
c
cc      do j=1,jgr
cc      write(1,*) x(j),y(j),kgr(j)
cc      enddo
c
      write(*,*) 'EVALUATING DISTANCES AND GENERATING GRAIN TABLE'
      write(*,*)
c
      do i=0,npts1+1
c
      write(*,'(1H+,a,i7)') 'X COORDINATE = ',i
c
      do j=0,npts2+1
      do k=0,npts3+1
c
      xg=float(i)-0.5
      yg=float(j)-0.5
      zg=float(k)-0.5
c
cc      write(2,*) xg,yg
c
      dist=1000000.
      do kk=1,jgr
      dist1=sqrt(((xg-x(kk))/aspect(1))**2+((yg-y(kk))/aspect(2))**2+
     # ((zg-z(kk))/aspect(3))**2)
      if(dist1.lt.dist) then
      dist=dist1
      kgrain(i,j,k)=kgr(kk)
      endif
      enddo
c
      enddo
      enddo
      enddo
cv
cv      do kv=1,1000
c
cv      i=int(ran(idum)*npts1)+1
cv      j=int(ran(idum)*npts2)+1
cv      k=int(ran(idum)*npts3)+1
cv      do i1=0,1
cv      do j1=0,1
cv      do k1=0,1
cv      ii=i+i1
cv      jj=j+j1
cv      kk=k+k1
cv      if(kgrain(ii,jj,kk).gt.0) then
cv       kgrain(ii,jj,kk)=-kgrain(ii,jj,kk)
cv      endif
cv      enddo
cv      enddo
cv      enddo
c
cv      enddo
cv
      do k=1,npts3
cz      write(6,*)'K= ',k
      write(16,*)'K= ',k
      do i=1,npts1
cz      write(6,'(514i5)')(kgrain(i,j,k),j=1,npts2)
      write(16,'(10000i5)')(kgrain(i,j,k),j=1,npts2)
      enddo
      enddo
c
      end
c
      function ran2(idum)
C---> Taken from 'Numerical Recipes'
      integer ia,im,iq,ir,idum,mask,k
      real*4 am,ran2
      parameter (ia=16807,im=2147483647,am=1./im,
     &          iq=127773,ir=2836,mask=123459876)
      idum=ieor(idum,mask)
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran2=am*idum
      idum=ieor(idum,mask)
      return
      end
