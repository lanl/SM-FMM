module various_functions
  use global
  implicit none

contains

  FUNCTION ran0(idum)
    INTEGER :: idum,IA,IM,IQ,IR,MASK
    DOUBLE PRECISION ::  ran0,AM
    PARAMETER (IA=16807,IM=2147483647,AM=1./IM, IQ=127773,IR=2836,MASK=123459876)
    ! Minimal" random number generator of Park and Miller. Returns a uniform random deviate
    ! between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
    ! to initialize the sequence; idum must not be altered between calls for successive deviates
    ! in a sequence.
    INTEGER k
    idum=ieor(idum,MASK)        ! XORing with MASK allows use of zero and other simple
    k=idum/IQ                   !   bit patterns for idum.
    idum=IA*(idum-k*IQ)-IR*k    ! Compute idum=mod(IA*idum,IM) without overflows by
    if (idum.lt.0) idum=idum+IM ! Schrage's method.
    ran0=AM*idum                ! Convert idum to a floating result.
    idum=ieor(idum,MASK)        ! Unmask before return.
    return
  END
 
  subroutine form_G
    use tensor_functions
    use fourier_functions
    implicit none
 
    double precision :: xk(3),g1(3,3,3,3),D_dft(3,3)
    double precision :: a(3,3)
    double precision :: g1tmp(3,3,3,3),aux3333(3,3,3,3)
  
    integer :: i,j,k,l,m,n,o,p,q,ii,itt,jj,kxx,kyy,kzz
    double precision :: dum,r,rold,tol,xisum,xjsum,xksum,xmix
  
    do kzz = 1,npts3
    do kyy = 1,npts2
    do kxx = 1,npts1
  
      xk = spatial_freq(kxx,kyy,kzz,npts1,npts2,npts3)
  
!      if(((mod(npts1,2)==0 .and. kxx==npts1/2+1).or. &
!          (mod(npts2,2)==0 .and. kyy==npts2/2+1).or. &
!         ((mod(npts3,2)==0 .and. kzz==npts3/2+1).and.npts3>1)).and.igamma.eq.0) then
!  
!        g1 = -s0
!  
!      else
  
        ! Build discrete and continuous D operator
        if(igamma.eq.0)then
  
          do ii = 1,3
          do jj = 1,3
            D_dft(ii,jj) = xk(ii)*xk(jj)
          enddo
          enddo
  
        elseif(igamma.eq.1)then
  
          xisum = 2.*pi*float(kxx-1)/float(npts1)
          xjsum = 2.*pi*float(kyy-1)/float(npts2)
          xksum = 2.*pi*float(kzz-1)/float(npts3)
  
          D_dft(1,1) = 2.0*(cos(xisum)-1.0)
          D_dft(2,2) = 2.0*(cos(xjsum)-1.0)
          D_dft(3,3) = 2.0*(cos(xksum)-1.0)
          D_dft(2,1) = -sin(xisum)*sin(xjsum) ! 0.5*(cos(xisum+xjsum)-cos(xisum-xjsum))
          D_dft(3,1) = -sin(xisum)*sin(xksum) ! 0.5*(cos(xisum+xksum)-cos(xisum-xksum))
          D_dft(3,2) = -sin(xjsum)*sin(xksum) ! 0.5*(cos(xjsum+xksum)-cos(xjsum-xksum))
          D_dft(1,2) = D_dft(2,1)
          D_dft(1,3) = D_dft(3,1)
          D_dft(2,3) = D_dft(3,2)
  
        endif
  
        ! Compute acoustic tensor inverse
        if (kxx==1 .and. kyy==1 .and. kzz==1) then
          a = 0.
        else
          do i = 1,3
          do k = 1,3
            a(i,k)=0.
            do j = 1,3
            do l = 1,3
              a(i,k) = a(i,k) + c0(i,j,k,l)*D_dft(j,l)!*xk(j)*xk(l)
            enddo
            enddo
          enddo
          enddo
  
          call lu_inverse(a, 3)
        endif
  
        do p = 1,3
        do q = 1,3
        do i = 1,3
        do j = 1,3
          g1(p,q,i,j) = -a(p,i)*D_dft(q,j) !*xk(q)*xk(j)
        enddo
        enddo
        enddo
        enddo
  
!      endif
  
!      if(igamma.eq.1)then
!  
!  !      c0n=sqrt(sum(c03333**2))
!        g1tmp = g1
!        tol = 1.0e-50
!        r = 2.0*tol
!        rold = 2.0e2
!        itt = 0
!        xmix = 1.0
!        do while(r.gt.tol.and.itt.lt.4)
!          itt=itt+1
!  
!          do i = 1,3
!          do j = 1,3
!          do k = 1,3
!          do l = 1,3
!            dum = 0.0
!            do m = 1,3
!            do n = 1,3
!            do o = 1,3
!            do p = 1,3
!              dum = dum - g1(i,j,m,n)*c0(m,n,o,p)*g1tmp(o,p,k,l)
!            enddo
!            enddo
!            enddo
!            enddo
!            aux3333(i,j,k,l) = dum
!          enddo
!          enddo
!          enddo
!          enddo
!  
!          g1tmp = xmix*aux3333 + (1.0-xmix)*g1tmp
!  
!        enddo
!        g1 = g1tmp
!  
!      endif
  
      Goper(:,:,:,:,kxx,kyy,kzz) = g1
  
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
        fourgrid3333(kzz,kyy,kxx,l,k,j,i) = cmplx(Goper(i,j,k,l,kxx,kyy,kzz))
      enddo
      enddo
      enddo
      enddo
  
    enddo
    enddo
    enddo
  
    call ifft_tensor3333(iplan_advanced3333, fourgrid3333)
 
    do kzz = 1,npts3
    do kyy = 1,npts2
    do kxx = 1,npts1
    do i = 1,3
    do j = 1,3
    do k = 1,3
    do l = 1,3
      Goperreal(i,j,k,l,kxx,kyy,kzz) = real(fourgrid3333(kzz,kyy,kxx,l,k,j,i))
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  
  end

  subroutine form_G_box
    use tensor_functions
    use fourier_functions
    implicit none
 
    double precision :: xk(3),g1(3,3,3,3),D_dft(3,3)
    double precision :: a(3,3)
    double precision :: g1tmp(3,3,3,3),aux3333(3,3,3,3)
  
    integer :: i,j,k,l,m,n,o,p,q,ii,itt,jj,kxx,kyy,kzz
    double precision :: dum,r,rold,tol,xisum,xjsum,xksum,xmix
  
    do kzz = 1,nbox3
    do kyy = 1,nbox2
    do kxx = 1,nbox1
  
      xk = spatial_freq(kxx,kyy,kzz,nbox1,nbox2,nbox3)
  
!      if(((mod(npts1,2)==0 .and. kxx==npts1/2+1).or. &
!          (mod(npts2,2)==0 .and. kyy==npts2/2+1).or. &
!         ((mod(npts3,2)==0 .and. kzz==npts3/2+1).and.npts3>1)).and.igamma.eq.0) then
!  
!        g1 = -s0
!  
!      else
  
        ! Build discrete and continuous D operator
        if(igamma.eq.0)then
  
          do ii = 1,3
          do jj = 1,3
            D_dft(ii,jj) = xk(ii)*xk(jj)
          enddo
          enddo
  
        elseif(igamma.eq.1)then
  
          xisum = 2.*pi*float(kxx-1)/float(npts1)
          xjsum = 2.*pi*float(kyy-1)/float(npts2)
          xksum = 2.*pi*float(kzz-1)/float(npts3)
  
          D_dft(1,1) = 2.0*(cos(xisum)-1.0)
          D_dft(2,2) = 2.0*(cos(xjsum)-1.0)
          D_dft(3,3) = 2.0*(cos(xksum)-1.0)
          D_dft(2,1) = -sin(xisum)*sin(xjsum) ! 0.5*(cos(xisum+xjsum)-cos(xisum-xjsum))
          D_dft(3,1) = -sin(xisum)*sin(xksum) ! 0.5*(cos(xisum+xksum)-cos(xisum-xksum))
          D_dft(3,2) = -sin(xjsum)*sin(xksum) ! 0.5*(cos(xjsum+xksum)-cos(xjsum-xksum))
          D_dft(1,2) = D_dft(2,1)
          D_dft(1,3) = D_dft(3,1)
          D_dft(2,3) = D_dft(3,2)
  
        endif
  
        ! Compute acoustic tensor inverse
        if (kxx==1 .and. kyy==1 .and. kzz==1) then
          a = 0.
        else
          do i = 1,3
          do k = 1,3
            a(i,k)=0.
            do j = 1,3
            do l = 1,3
              a(i,k) = a(i,k) + c0(i,j,k,l)*D_dft(j,l)!*xk(j)*xk(l)
            enddo
            enddo
          enddo
          enddo
  
          call lu_inverse(a, 3)
        endif
  
        do p = 1,3
        do q = 1,3
        do i = 1,3
        do j = 1,3
          g1(p,q,i,j) = -a(p,i)*D_dft(q,j) !*xk(q)*xk(j)
        enddo
        enddo
        enddo
        enddo
  
!      endif
  
!      if(igamma.eq.1)then
!  
!  !      c0n=sqrt(sum(c03333**2))
!        g1tmp = g1
!        tol = 1.0e-50
!        r = 2.0*tol
!        rold = 2.0e2
!        itt = 0
!        xmix = 1.0
!        do while(r.gt.tol.and.itt.lt.4)
!          itt=itt+1
!  
!          do i = 1,3
!          do j = 1,3
!          do k = 1,3
!          do l = 1,3
!            dum = 0.0
!            do m = 1,3
!            do n = 1,3
!            do o = 1,3
!            do p = 1,3
!              dum = dum - g1(i,j,m,n)*c0(m,n,o,p)*g1tmp(o,p,k,l)
!            enddo
!            enddo
!            enddo
!            enddo
!            aux3333(i,j,k,l) = dum
!          enddo
!          enddo
!          enddo
!          enddo
!  
!          g1tmp = xmix*aux3333 + (1.0-xmix)*g1tmp
!  
!        enddo
!        g1 = g1tmp
!  
!      endif
  
      Goper_box(:,:,:,:,kxx,kyy,kzz) = g1
  
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
        fourgrid3333box(kzz,kyy,kxx,l,k,j,i) = cmplx(Goper_box(i,j,k,l,kxx,kyy,kzz))
      enddo
      enddo
      enddo
      enddo
  
    enddo
    enddo
    enddo
  
    call ifft_tensor3333box(iplan_advanced3333box, fourgrid3333box)
 
    do kzz = 1,nbox3
    do kyy = 1,nbox2
    do kxx = 1,nbox1
    do i = 1,3
    do j = 1,3
    do k = 1,3
    do l = 1,3
      Goperreal_box(i,j,k,l,kxx,kyy,kzz) = real(fourgrid3333box(kzz,kyy,kxx,l,k,j,i))
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  
  end

  subroutine integrate_G_intpt(iel1,iel2,elements,G)
    use types
    use tensor_functions

    integer, intent(in) :: iel1, iel2
    type (element_type), allocatable, intent (in) :: elements(:)
    type (element_type) :: element1
    type (element_type) :: element2

    double precision, intent(out) :: G(3,3,3,3)

    integer :: ip1, ip2
    double precision :: dx(3), v

    element1 = elements(iel1)
    element2 = elements(iel2)
    G = 0.0
    do ip1 = 1,size(element1%w)
      do ip2 = 1,size(element2%w)
        ! write(*,*)'ip1,ip2',ip1,ip2
        dx = element1%xintpt(:,ip1) - element2%xintpt(:,ip2)
        ! write(*,*)'dx',dx
        ! write(*,*)enforce_periodic(nint(dx(1)) + 1, 1, npts1)
        ! write(*,*)enforce_periodic(nint(dx(2)) + 1, 1, npts2)
        ! write(*,*)enforce_periodic(nint(dx(3)) + 1, 1, npts3)
        v = element1%w(ip1)*element1%V*element2%w(ip2)*element2%V
        ! write(*,*)'v',v
        G = G + get_G(dx)*V
        ! write(*,*)'get_G(dx)',get_G(dx)
      enddo
    enddo
    G = G/element1%V

  end

  function get_G(dx) result(G)
    use tensor_functions

    double precision, intent(in) :: dx(3)
    double precision :: G(3,3,3,3)
    integer :: i, j, k

    i = enforce_periodic(nint(dx(1)) + 1, 1, npts1)
    j = enforce_periodic(nint(dx(2)) + 1, 1, npts2)
    k = enforce_periodic(nint(dx(3)) + 1, 1, npts3)
    G = Goperreal(:,:,:,:,i,j,k) 

  end

  subroutine form_Gmat_full_intpt(elements)
    use global
    use types
    use tensor_functions

    type (element_type), allocatable, intent (in) :: elements(:)
    double precision :: G(3,3,3,3)
    integer :: iel1, iel2

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(iel1, iel2, G) &
    !$OMP SHARED(elements, plot_el_scalar, nel, Gmat)
    !$OMP DO COLLAPSE(2)
    do iel1 = 1, nel
      do iel2 = 1, nel
        ! write(*,*) iel1, iel2
        call integrate_G_intpt(iel1,iel2,elements,G)
        Gmat(:,:,:,:,iel1,iel2) = G
        plot_el_scalar(iel2) = norm2(G)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  end

  subroutine form_Gmat_full(elements)
    use global
    use types
    use tensor_functions

    type (element_type), allocatable, intent (in) :: elements(:)
    double precision :: G(3,3,3,3), dx(3)
    double precision :: r = 1.0e5
    integer :: iel1, iel2

    open(111,file='full_Gamma',form='unformatted',access='sequential',status='unknown')

    if (iload == 0) then

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(iel1, iel2, G, dx) &
      !$OMP SHARED(elements, plot_el_scalar, nel, Gmat, npts1,npts2,npts3,r)
      !$OMP DO COLLAPSE(2)
      do iel1 = 1, nel
        do iel2 = 1, nel
          ! write(*,*) iel1, iel2
          dx = distance_periodic(elements(iel1)%x,elements(iel2)%x,npts1,npts2,npts3)
          if (norm2(dx) <= r) then
            call integrate_G(iel1,iel2,elements,G)
            write(*,*) norm2(G)
            Gmat(:,:,:,:,iel1,iel2) = G
          else
            Gmat(:,:,:,:,iel1,iel2) = 0.0
          endif
          ! plot_el_scalar(iel2) = norm2(G)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
  
      write(111)Gmat

    else

      read(111) Gmat

    endif

    close(111)

  end


  subroutine integrate_G(iel1,iel2,elements,G)
    use types
    use tensor_functions

    integer, intent(in) :: iel1, iel2
    type (element_type), allocatable, intent (in) :: elements(:)
    type (element_type) :: element1
    type (element_type) :: element2

    double precision, intent(out) :: G(3,3,3,3)

    integer :: ip1, ip2
    double precision :: dx(3), v

    element1 = elements(iel1)
    element2 = elements(iel2)
    G = 0.0
    do ip1 = 1,element1%ngrid_ind
      do ip2 = 1,element2%ngrid_ind
        ! write(*,*)'ip1,ip2',ip1,ip2
        dx = float(element1%grid_ind(:,ip1)) - float(element2%grid_ind(:,ip2))
        ! write(*,*)'element1%grid_ind(:,ip1)',element1%grid_ind(:,ip1)
        ! write(*,*)'element2%grid_ind(:,ip2)',element2%grid_ind(:,ip2)
        ! write(*,*)enforce_periodic(nint(dx(1)) + 1, 1, npts1)
        ! write(*,*)enforce_periodic(nint(dx(2)) + 1, 1, npts2)
        ! write(*,*)enforce_periodic(nint(dx(3)) + 1, 1, npts3)
        v = 1.0
        ! write(*,*)'v',v
        G = G + get_G(dx)*V
        ! write(*,*)'get_G(dx)',get_G(dx)
      enddo
    enddo
    if (element1%ngrid_ind > 0) G = G/float(element1%ngrid_ind)

  end

  subroutine grains_to_elements(grains, elements)
    use global
    use types
    use tensor_functions

    type (grain_type), allocatable, intent (in) :: grains(:)
    type (element_type), allocatable, intent (inout) :: elements(:)

    integer :: iel, id

    do iel = 1, nel

      id = elements(iel)%grain_id
      elements(iel)%C = grains(id)%C
      elements(iel)%Q = grains(id)%Q

    enddo
 
  end

  subroutine initial_guess(elements)
    use global
    use types
    use tensor_functions

    type (element_type), allocatable, intent (inout) :: elements(:)

    integer :: iel

    polaravg = 0.0
    do iel = 1, nel

      elements(iel)%strain = Eapp
      elements(iel)%stress = TijklTkl(elements(iel)%C, elements(iel)%strain)
      elements(iel)%polar = elements(iel)%stress - TijklTkl(c0(:,:,:,:), elements(iel)%strain)
      polaravg = polaravg + elements(iel)%polar*elements(iel)%v

    enddo
    polaravg = polaravg/vtot
 
  end

  function volume_tetra(xnode) result(v)
    use tensor_functions
    double precision, intent(in) :: xnode(3,4)
    double precision :: v, volmat(3,3)
    double precision, parameter ::  one_sixth = 0.166666666666667

    volmat(:,1) = xnode(:,1) - xnode(:,4)
    volmat(:,2) = xnode(:,2) - xnode(:,4)
    volmat(:,3) = xnode(:,3) - xnode(:,4)

    v = (determinant33(volmat))*one_sixth

  end

  function point_in_tetra(x,xnode) result(point_in)
    double precision, intent(in) :: x(3), xnode(3,4)
    double precision :: xnode_tmp(3,4), xi(4), v
    logical :: point_in
    integer :: i

    ! barycentric
    v = volume_tetra(xnode)
    do i = 1, 4
      xnode_tmp = xnode
      xnode_tmp(:,i) = x
      xi(i) = volume_tetra(xnode_tmp)/v
    enddo
    if (xi(1) > 0 .and. xi(2) > 0 .and. xi(3) > 0 .and. xi(4) > 0) then
      point_in = .true.
    else
      point_in = .false.
    endif

  end

  subroutine points_to_elements(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: xmin(3), xmax(3), x(3), xref(3), dx(3), dxavg(3), xrefold(3)
  double precision, allocatable :: x_grid(:,:)
  integer, allocatable :: grid_ind(:,:)
  integer :: iel, i, gridmin(3), gridmax(3), ip1, ip2, ip3, ic, iteravg

  do ip3 = 1,npts3
    do ip2 = 1,npts2
      do ip1 = 1,npts1

        iel = el_id(ip1,ip2,ip3)
        elements(iel)%ngrid_ind = elements(iel)%ngrid_ind + 1

      enddo
    enddo
  enddo

  vtot = 0.0
  do iel = 1, nel

    elements(iel)%v = float(elements(iel)%ngrid_ind)
    allocate(elements(iel)%grid_ind(3,elements(iel)%ngrid_ind))
    allocate(elements(iel)%x_grid(3,elements(iel)%ngrid_ind))
    elements(iel)%ngrid_ind = 0 ! reset

    vtot = vtot + elements(iel)%v

  enddo

  do ip3 = 1,npts3
    do ip2 = 1,npts2
      do ip1 = 1,npts1

        iel = el_id(ip1,ip2,ip3)
        elements(iel)%ngrid_ind = elements(iel)%ngrid_ind + 1
        elements(iel)%grid_ind(:,elements(iel)%ngrid_ind) = [ip1,ip2,ip3]
        elements(iel)%x_grid(:,elements(iel)%ngrid_ind) = float([ip1,ip2,ip3])

      enddo
    enddo
  enddo
  
  do iel = 1, nel

    ! write(*,*)  elements(iel)%ngrid_ind
    xref = elements(iel)%x_grid(:,1)
    ! xrefold = xref
    err = 1.0
    do while(err > 1.0e-4)
      dxavg = 0.0
      do i = 1,elements(iel)%ngrid_ind
        dx = distance_periodic(elements(iel)%x_grid(:,i),xref,npts1,npts2,npts3)
        elements(iel)%x_grid(:,i) = xref + dx
        dxavg = dxavg + dx
      enddo
      dxavg = dxavg/float(elements(iel)%ngrid_ind)
  
      xrefold = xref
      xref = xref + dxavg*1.0
      if (xref(1)<0.5)xref(1) = xref(1) + float(npts1)
      if (xref(1)>float(npts1)+0.5)xref(1) = xref(1) - float(npts1)
      if (xref(2)<0.5)xref(2) = xref(2) + float(npts2)
      if (xref(2)>float(npts2)+0.5)xref(2) = xref(2) - float(npts2)
      if (xref(3)<0.5)xref(3) = xref(3) + float(npts3)
      if (xref(3)>float(npts3)+0.5)xref(3) = xref(3) - float(npts3)
      err = norm2(xrefold - xref)
    enddo
    elements(iel)%x = xref

    xrefold = xref
    xref = 0.0
    do i = 1,elements(iel)%ngrid_ind
      xref = xref + elements(iel)%x_grid(:,i)
    enddo
    xref = xref/float(elements(iel)%ngrid_ind)
    write(*,*) norm2(xrefold - xref)
  enddo

  write(*,*) 'unassigned'
  do ip1 = 1,npts1
    do ip2 = 1,npts2
      do ip3 = 1,npts3

        if (el_id(ip1,ip2,ip3) == 0) then 
          write(*,*) ip1,ip2,ip3
        endif

      enddo
    enddo
  enddo

end

function distance_periodic(x1,x2,npts1,npts2,npts3) result(dxmin)
  implicit none
  
  double precision, intent(in) :: x1(3), x2(3)
  integer, intent(in) :: npts1,npts2,npts3
  
  double precision :: dx(3), x2per(3), dxmin(3), dxnmin, dxn
  integer :: iper1,iper2,iper3
  
  dxnmin = 1.0e10
  do iper1 = -1,1
    x2per(1) = x2(1) + float(iper1*npts1)
    do iper2 = -1,1
      x2per(2) = x2(2) + float(iper2*npts2)
      do iper3 = -1,1
        x2per(3) = x2(3) + float(iper3*npts3)
  
        dx = x1 - x2per
        dxn = norm2(dx)
  
        if (dxn < dxnmin) then
          dxnmin = dxn
          dxmin = dx
        endif
  
      enddo
    enddo
  enddo

end

subroutine points_to_boxes(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: edge_lngth(3), xstart(3), xend(3)
  integer :: ib, ib1, ib2, ib3

  ! grid points in each box (for interaction calculation analogous to element approach)
  edge_lngth(1) = (float(npts1)/float(nbox1)) ! length of each box edge (domain 0.5 to npts + 0.5)
  edge_lngth(2) = (float(npts2)/float(nbox2)) 
  edge_lngth(3) = (float(npts3)/float(nbox3)) 
  ib = 0
  do ib1 = 1,nbox1
   do ib2 = 1,nbox2
    do ib3 = 1,nbox3
     ib = ib + 1

     xstart(1) = 0.5 + edge_lngth(1)*float(ib1 - 1)
     xstart(2) = 0.5 + edge_lngth(2)*float(ib2 - 1)
     xstart(3) = 0.5 + edge_lngth(3)*float(ib3 - 1)
     xend = xstart + edge_lngth
     boxes(ib)%ipstart = ceiling(xstart)
     boxes(ib)%ipend = floor(xend)
     boxes(ib)%ngrid_ind = (boxes(ib)%ipend(1) - boxes(ib)%ipstart(1))*(boxes(ib)%ipend(2) - boxes(ib)%ipstart(2))* &
       (boxes(ib)%ipend(3) - boxes(ib)%ipstart(3))
     boxes(ib)%edge_lngth = edge_lngth
     boxes(ib)%xstart = xstart
     boxes(ib)%xend = xend
     boxes(ib)%x = (xstart + xend)*0.5
     boxes(ib)%ib_grid = [ib1,ib2,ib3]

   enddo
  enddo
 enddo

end


subroutine elements_to_boxes(elements,boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: xc(3)
  integer, allocatable :: element_ids(:)
  integer :: ib, itot, ic, iel


  allocate(element_ids(nel))

  ! elements in each box
  itot = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ic,iel,xc,element_ids) &
  !$OMP SHARED(nbox, boxes, elements, nel) &
  !$OMP REDUCTION(+:itot)
  !$OMP DO 
  do ib = 1, nbox

     ic = 0
     do iel = 1,nel
       xc = elements(iel)%x
       if (xc(1).ge.boxes(ib)%xstart(1).and.xc(1).lt.boxes(ib)%xend(1).and. &
           xc(2).ge.boxes(ib)%xstart(2).and.xc(2).lt.boxes(ib)%xend(2).and. &
           xc(3).ge.boxes(ib)%xstart(3).and.xc(3).lt.boxes(ib)%xend(3).and.elements(iel)%v > 0.0) then
         ic = ic + 1
         itot = itot + 1
         element_ids(ic) = iel
         elements(iel)%box_id = ib
       endif
     enddo
     boxes(ib)%nel = ic
     allocate(boxes(ib)%element_ids(ic))
     boxes(ib)%element_ids(1:ic) = element_ids(1:ic)

  enddo 
  !$OMP END DO
  !$OMP END PARALLEL

  if (itot.ne.nel) then
    write(*,*)'itot,nel',itot,nel
    ! stop
  endif

  deallocate(element_ids)

end

subroutine boxes_near_far(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: x(3)
  integer :: ib1, ib2, dib(3), icnear, icfar
  integer, allocatable :: boxnear(:), boxfar(:)

  allocate(boxnear(nbox))
  allocate(boxfar(nbox))

  do ib1 = 1,nbox

    icnear = 0
    icfar = 0
    do ib2 = 1,nbox

      dib = abs(boxes(ib1)%ib_grid - boxes(ib2)%ib_grid)
      if (dib(1).eq.nbox1 - 1) dib(1) = 1
      if (dib(2).eq.nbox2 - 1) dib(2) = 1
      if (dib(3).eq.nbox3 - 1) dib(3) = 1

      if (maxval(dib).le.1) then ! near field
        icnear = icnear + 1
        boxnear(icnear) = ib2
      else ! far field
        icfar = icfar + 1
        boxfar(icfar) = ib2
      endif

    enddo

    boxes(ib1)%nnear = icnear
    boxes(ib1)%nfar = icfar
    allocate(boxes(ib1)%near_box_id(icnear))
    allocate(boxes(ib1)%far_box_id(icfar))
    boxes(ib1)%near_box_id = boxnear(1:icnear)
    boxes(ib1)%far_box_id = boxfar(1:icfar)

  enddo 

 deallocate(boxnear)
 deallocate(boxfar)

end

subroutine element_moments(elements,boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dx(3), x(3), dxdx(3,3)
  integer :: iel, ib, ip, i, j
  
  ! element momments within boxes
  do iel = 1,nel ! loop over elements
    ib = elements(iel)%box_id
 
    elements(iel)%dx = 0.0
    elements(iel)%dxdx = 0.0
    do ip = 1,elements(iel)%ngrid_ind ! loop over points in each element
      x(1) = elements(iel)%x_grid(1,ip)
      x(2) = elements(iel)%x_grid(2,ip)
      x(3) = elements(iel)%x_grid(3,ip)
      dx = - boxes(ib)%x + x
      elements(iel)%dx = elements(iel)%dx + dx
      do i = 1,3
       do j = 1,3
         elements(iel)%dxdx(i,j) = elements(iel)%dxdx(i,j) + dx(i)*dx(j)
       enddo
      enddo
    enddo
    if (elements(iel)%v > 0) elements(iel)%dx = elements(iel)%dx/elements(iel)%v
    if (elements(iel)%v > 0) elements(iel)%dxdx = elements(iel)%dxdx/elements(iel)%v

  enddo

!  iel = 1401
!  ib = elements(iel)%box_id
!  write(*,*) iel 
!  write(*,*) elements(iel)%dx(1)**2
!  write(*,*) elements(iel)%dxdx(1,1)
!  x = 0.0
!  do ip = 1,elements(iel)%ngrid_ind
!    x = x + elements(iel)%x_grid(:,ip) - boxes(ib)%x
!  enddo
!  x = x/float(elements(iel)%ngrid_ind)
!  write(*,*) x(1)**2
!  dxdx = 0.0
!  do ip = 1,elements(iel)%ngrid_ind
!    ! write(*,*) elements(iel)%x_grid(:,ip)
!    do i = 1,3
!      do j = 1,3
!        ! dxdx(i,j) = dxdx(i,j) + (elements(iel)%x_grid(i,ip) - x(i))*(elements(iel)%x_grid(j,ip) - x(j))
!        dxdx(i,j) = dxdx(i,j) + (elements(iel)%x_grid(i,ip) - boxes(ib)%x(i))*(elements(iel)%x_grid(j,ip) - boxes(ib)%x(j))
!      enddo
!    enddo
!  enddo
!  dxdx = dxdx/float(elements(iel)%ngrid_ind)
!  write(*,*) dxdx(1,1)

end

subroutine element_close_inter_ids(elements,boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (in) :: boxes(:)

  integer, allocatable :: close_inter_el_ids(:)
  integer :: ib, iel_in_box1, iel1, ic, ib_in_near, ibnear, iel2, iel_in_box2

  allocate(close_inter_el_ids(nel))
 
  ! close interaction element sets
  nclose_inter_el_ids_mx = 0
  do ib = 1,nbox ! loop over boxes
    do iel_in_box1 = 1,boxes(ib)%nel ! loop over elements in each box
      iel1 = boxes(ib)%element_ids(iel_in_box1)

      ic = 0
      do ib_in_near = 1,boxes(ib)%nnear ! loop over near field boxes fr this box
        ibnear = boxes(ib)%near_box_id(ib_in_near) 
        do iel_in_box2 = 1,boxes(ibnear)%nel ! loop over elements in near field boxes

          iel2 = boxes(ibnear)%element_ids(iel_in_box2)
          if (iel1 .ne. iel2) then
            ic = ic + 1
            close_inter_el_ids(ic) = iel2
          endif
        enddo
      enddo

      if (ic > nclose_inter_el_ids_mx) nclose_inter_el_ids_mx = ic
      elements(iel1)%nel_close = ic
      if (allocated(elements(iel1)%close_inter_el_ids)) then
        write(*,*) 'ib, iel1:', ib, iel1
      endif
      allocate(elements(iel1)%close_inter_el_ids(ic))
      elements(iel1)%close_inter_el_ids(1:ic) = close_inter_el_ids(1:ic)
    enddo
  enddo

  deallocate(close_inter_el_ids)

  write(*,*) 'nclose_inter_el_ids_mx', nclose_inter_el_ids_mx

end

subroutine element_close_inter_Gamma(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: G(3,3,3,3)
  integer :: iel1, iel2, iel_in_close

  do iel1 = 1, nel
    allocate(elements(iel1)%Gamma(3,3,3,3,elements(iel1)%nel_close))
  enddo

!  ! interactions
!  !$OMP PARALLEL DEFAULT(NONE) &
!  !$OMP PRIVATE(iel1, iel2, G, iel_in_close) &
!  !$OMP SHARED(elements, nel)
!  !$OMP DO
!  do iel1 = 1, nel
!    do iel_in_close = 1,elements(iel1)%nel_close
!      iel2 = elements(iel1)%close_inter_el_ids(iel_in_close)
!      ! write(*,*) iel1,iel2
!
!      call integrate_G(iel1,iel2,elements,G)
!      elements(iel1)%Gamma(:,:,:,:,iel_in_close) = G
!
!    enddo
!  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL

  ! interactions
  open(111,file='elements_Gamma',form='unformatted',access='sequential',status='unknown')

  if (iload == 0) then

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(iel1, iel2, G, iel_in_close) &
    !$OMP SHARED(elements, nel, nclose_inter_el_ids_mx)
    !$OMP DO COLLAPSE(2)
    do iel1 = 1, nel
      do iel_in_close = 1,nclose_inter_el_ids_mx
        ! write(*,*) iel1
        if (iel_in_close <= size(elements(iel1)%close_inter_el_ids))then
          iel2 = elements(iel1)%close_inter_el_ids(iel_in_close)
          ! write(*,*) iel1,iel2
    
          call integrate_G(iel1,iel2,elements,G)
          elements(iel1)%Gamma(:,:,:,:,iel_in_close) = G
        endif
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    ! self interaction
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(iel1, G) &
    !$OMP SHARED(elements, nel)
    !$OMP DO
    do iel1 = 1, nel
      call integrate_G(iel1,iel1,elements,G)
      elements(iel1)%Gamma_self = G
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    do iel1 = 1, nel
      write(111) elements(iel1)%Gamma
      write(111) elements(iel1)%Gamma_self
    enddo

  else

    do iel1 = 1, nel
      read(111) elements(iel1)%Gamma
      read(111) elements(iel1)%Gamma_self
    enddo

  endif

  close(111)


end

subroutine box_far_inter_Gamma(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dx(3), dir1_inc(3), dir2_inc(3), dir3_inc(3), inc
  integer :: ib, i, ibfar, ind(3)

  do ib = 1,nbox
    allocate(boxes(ib)%Gamma(3,3,3,3,boxes(ib)%nfar))
    if (order > 0) allocate(boxes(ib)%dGamma_dX(3,3,3,3,3,boxes(ib)%nfar))
    if (order > 1) allocate(boxes(ib)%d2Gamma_dX2(3,3,3,3,3,3,boxes(ib)%nfar))
  enddo

  ! far field interactions
  ! inc = float(nint(boxes(1)%edge_lngth(1)/2.0))/1.0! 1.0
  inc = 1.0
  write(*,*) 'inc for derivatives of Gamma is ', inc
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,i, ibfar, dx, ind) &
  !$OMP SHARED(nbox, boxes, inc, dir1_inc, dir2_inc, dir3_inc,order,nbox1,nbox2,nbox3, Goperreal_box)
  !$OMP DO 
  do ib = 1,nbox
    do i = 1,boxes(ib)%nfar

      ibfar = boxes(ib)%far_box_id(i)
      ! write(*,*) ib,ibfar
      dx = - boxes(ibfar)%x + boxes(ib)%x
      dx = float(nint(dx))
        
      boxes(ib)%Gamma(:,:,:,:,i) = get_G(dx)

      ind = - boxes(ibfar)%ib_grid + boxes(ib)%ib_grid + 1
      ind(1) = enforce_periodic(ind(1),1,nbox1)
      ind(2) = enforce_periodic(ind(2),1,nbox2)
      ind(3) = enforce_periodic(ind(3),1,nbox3)

      ! boxes(ib)%Gamma(:,:,:,:,i) = Goperreal_box(:,:,:,:,ind(1),ind(2),ind(3))/ &
      !  float((boxes(ib)%ipend(1)-boxes(ib)%ipstart(1))*(boxes(ib)%ipend(2)-boxes(ib)%ipstart(2))* &
      !  (boxes(ib)%ipend(3)-boxes(ib)%ipstart(3)))
  
      if (order > 0) then
        boxes(ib)%dGamma_dX(:,:,:,:,1,i) = (get_G(dx + dir1_inc) - get_G(dx - dir1_inc))/(2.0*inc)
        boxes(ib)%dGamma_dX(:,:,:,:,2,i) = (get_G(dx + dir2_inc) - get_G(dx - dir2_inc))/(2.0*inc)
        boxes(ib)%dGamma_dX(:,:,:,:,3,i) = (get_G(dx + dir3_inc) - get_G(dx - dir3_inc))/(2.0*inc)
      endif
 
      if (order > 1) then
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,1,i) = &
         (get_G(dx + dir1_inc) - 2.0*get_G(dx) + get_G(dx - dir1_inc))/(inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,2,2,i) = &
         (get_G(dx + dir2_inc) - 2.0*get_G(dx) + get_G(dx - dir2_inc))/(inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,3,i) = &
         (get_G(dx + dir3_inc) - 2.0*get_G(dx) + get_G(dx - dir3_inc))/(inc**2)
  
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,2,i) = &
         (get_G(dx + dir1_inc + dir2_inc) - get_G(dx + dir1_inc - dir2_inc) - &
          get_G(dx - dir1_inc + dir2_inc) + get_G(dx - dir1_inc - dir2_inc))/(4.0*inc**2)
  
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,3,i) = &
         (get_G(dx + dir1_inc + dir3_inc) - get_G(dx + dir1_inc - dir3_inc) - &
          get_G(dx - dir1_inc + dir3_inc) + get_G(dx - dir1_inc - dir3_inc))/(4.0*inc**2)
        
        boxes(ib)%d2Gamma_dX2(:,:,:,:,2,3,i) = &
         (get_G(dx + dir2_inc + dir3_inc) - get_G(dx + dir2_inc - dir3_inc) - &
          get_G(dx - dir2_inc + dir3_inc) + get_G(dx - dir2_inc - dir3_inc))/(4.0*inc**2)
  
        boxes(ib)%d2Gamma_dX2(:,:,:,:,2,1,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,1,2,i)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,1,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,1,3,i)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,2,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,2,3,i)
      endif

      ! boxes(ib)%dGamma_dX = 0.0
      ! boxes(ib)%d2Gamma_dX2 = 0.0
  
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine box_far_inter_Gamma_FFT(boxes)
  use global
  use types
  use tensor_functions
  use fourier_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dx(3), dir1_inc(3), dir2_inc(3), dir3_inc(3), inc
  integer :: ib, i, ibfar, ind(3), j, k, l, m, n, ib1, ib2, ib3, ibfartmp

  inc = 1.0
  write(*,*) 'inc for derivatives of Gamma is ', inc
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc
  ib = 1
  Gamma_box = 0.0
  dGamma_dX_box = 0.0
  d2Gamma_dX2_box = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i, ibfar, dx, ind) &
  !$OMP SHARED(nbox, boxes, inc, dir1_inc, dir2_inc, dir3_inc,order,nbox1,nbox2,nbox3) &
  !$OMP SHARED(Gamma_box, dGamma_dX_box, d2Gamma_dX2_box, ib)
  !$OMP DO 
  do i = 1,boxes(ib)%nfar

    ibfar = boxes(ib)%far_box_id(i)
    dx = boxes(ibfar)%x - boxes(ib)%x ! from ib to ibfar
    dx = float(nint(dx))
    ind = boxes(ibfar)%ib_grid 

    Gamma_box(:,:,:,:,ind(1),ind(2),ind(3)) = get_G(dx)

    ! ! verify
    ! do j = 1,boxes(ibfar)%nfar
    !   ibfartmp = boxes(ibfar)%far_box_id(j)
    !   if (ibfartmp == ib) exit
    ! enddo
    ! write(*,*)i,norm2(boxes(ibfar)%Gamma(:,:,:,:,j)-Gamma_box(:,:,:,:,ind(1),ind(2),ind(3)))

    if (order > 0) then
      dGamma_dX_box(:,:,:,:,1,ind(1),ind(2),ind(3)) = (get_G(dx + dir1_inc) - get_G(dx - dir1_inc))/(2.0*inc)
      dGamma_dX_box(:,:,:,:,2,ind(1),ind(2),ind(3)) = (get_G(dx + dir2_inc) - get_G(dx - dir2_inc))/(2.0*inc)
      dGamma_dX_box(:,:,:,:,3,ind(1),ind(2),ind(3)) = (get_G(dx + dir3_inc) - get_G(dx - dir3_inc))/(2.0*inc)
    endif

    if (order > 1) then
      d2Gamma_dX2_box(:,:,:,:,1,1,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir1_inc) - 2.0*get_G(dx) + get_G(dx - dir1_inc))/(inc**2)
      d2Gamma_dX2_box(:,:,:,:,2,2,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir2_inc) - 2.0*get_G(dx) + get_G(dx - dir2_inc))/(inc**2)
      d2Gamma_dX2_box(:,:,:,:,3,3,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir3_inc) - 2.0*get_G(dx) + get_G(dx - dir3_inc))/(inc**2)

      d2Gamma_dX2_box(:,:,:,:,1,2,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir1_inc + dir2_inc) - get_G(dx + dir1_inc - dir2_inc) - &
        get_G(dx - dir1_inc + dir2_inc) + get_G(dx - dir1_inc - dir2_inc))/(4.0*inc**2)

      d2Gamma_dX2_box(:,:,:,:,1,3,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir1_inc + dir3_inc) - get_G(dx + dir1_inc - dir3_inc) - &
        get_G(dx - dir1_inc + dir3_inc) + get_G(dx - dir1_inc - dir3_inc))/(4.0*inc**2)
      
      d2Gamma_dX2_box(:,:,:,:,2,3,ind(1),ind(2),ind(3)) = &
        (get_G(dx + dir2_inc + dir3_inc) - get_G(dx + dir2_inc - dir3_inc) - &
        get_G(dx - dir2_inc + dir3_inc) + get_G(dx - dir2_inc - dir3_inc))/(4.0*inc**2)

      d2Gamma_dX2_box(:,:,:,:,2,1,ind(1),ind(2),ind(3)) = d2Gamma_dX2_box(:,:,:,:,1,2,ind(1),ind(2),ind(3))
      d2Gamma_dX2_box(:,:,:,:,3,1,ind(1),ind(2),ind(3)) = d2Gamma_dX2_box(:,:,:,:,1,3,ind(1),ind(2),ind(3))
      d2Gamma_dX2_box(:,:,:,:,3,2,ind(1),ind(2),ind(3)) = d2Gamma_dX2_box(:,:,:,:,2,3,ind(1),ind(2),ind(3))
    endif

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  do ib3 = 1,nbox3
    do ib2 = 1,nbox2
      do ib1 = 1,nbox1
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                fourgrid3333box(ib3,ib2,ib1,l,k,j,i) = cmplx(Gamma_box(i,j,k,l,ib1,ib2,ib3))
                if (order > 0) then
                  do m = 1,3
                    fourgrid33333box(ib3,ib2,ib1,m,l,k,j,i) = cmplx(dGamma_dX_box(i,j,k,l,m,ib1,ib2,ib3))
                    if (order > 1) then
                      do n = 1,3
                        fourgrid333333box(ib3,ib2,ib1,n,m,l,k,j,i) = cmplx(d2Gamma_dX2_box(i,j,k,l,m,n,ib1,ib2,ib3))
                      enddo
                    endif
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
        ! write(*,*)fourgrid33333box(ib3,ib2,ib1,1,1,1,1,1)
      enddo
    enddo
  enddo

  call fft_tensor3333box(plan_advanced3333box, fourgrid3333box)
  if (order > 0) call fft_tensor33333box(plan_advanced33333box, fourgrid33333box)
  if (order > 1) call fft_tensor333333box(plan_advanced333333box, fourgrid333333box)

  do ib3 = 1,nbox3
    do ib2 = 1,nbox2
      do ib1 = 1,nbox1
        Gamma_box_hat(ib3,ib2,ib1,:,:,:,:) = fourgrid3333box(ib3,ib2,ib1,:,:,:,:)
        if (order > 0) dGamma_dX_box_hat(ib3,ib2,ib1,:,:,:,:,:) = fourgrid33333box(ib3,ib2,ib1,:,:,:,:,:)
        if (order > 1) d2Gamma_dX2_box_hat(ib3,ib2,ib1,:,:,:,:,:,:) = fourgrid333333box(ib3,ib2,ib1,:,:,:,:,:,:)
      enddo
    enddo
  enddo

end

subroutine outgoing_expansion(elements,boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: ib, iel_in_box, i, j, iel

  ! outgoing expansions
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,iel_in_box,iel,i,j) &
  !$OMP SHARED(boxes,elements,nbox,order)
  !$OMP DO 
  do ib = 1,nbox

    boxes(ib)%outgoing0 = 0.0
    boxes(ib)%outgoing1 = 0.0
    boxes(ib)%outgoing2 = 0.0
    do iel_in_box = 1,boxes(ib)%nel
      iel = boxes(ib)%element_ids(iel_in_box)
      boxes(ib)%outgoing0 = boxes(ib)%outgoing0 + elements(iel)%polar*elements(iel)%v
      if (order > 0) then
        do i = 1,3
          boxes(ib)%outgoing1(:,:,i) = boxes(ib)%outgoing1(:,:,i) + elements(iel)%polar*elements(iel)%dx(i)* &
          elements(iel)%v
          if (order > 1) then
            do j = 1,3
              boxes(ib)%outgoing2(:,:,i,j) = boxes(ib)%outgoing2(:,:,i,j) + elements(iel)%polar* &
              elements(iel)%dx(i)*elements(iel)%dx(j)*elements(iel)%v
            enddo
          endif
        enddo
      endif
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine incoming_expansion(elements,boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: ib, ibfar, i

  ! incoming expansions
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ibfar,i) &
  !$OMP SHARED(nbox,boxes,order)
  !$OMP DO 
  do ib = 1,nbox

    boxes(ib)%incoming0 = 0.0
    boxes(ib)%incoming1 = 0.0
    boxes(ib)%incoming2 = 0.0
    do i = 1,boxes(ib)%nfar
      ibfar = boxes(ib)%far_box_id(i)
      if (order == 0) then

        boxes(ib)%incoming0 = boxes(ib)%incoming0 + TijklTkl(boxes(ib)%Gamma(:,:,:,:,i), boxes(ibfar)%outgoing0)

      elseif (order == 1) then

        boxes(ib)%incoming0 = boxes(ib)%incoming0 + TijklTkl(boxes(ib)%Gamma(:,:,:,:,i), boxes(ibfar)%outgoing0) - &
          TijklmTklm(boxes(ib)%dGamma_dX(:,:,:,:,:,i), boxes(ibfar)%outgoing1)
  
        boxes(ib)%incoming1 = boxes(ib)%incoming1 + TijklmTkl(boxes(ib)%dGamma_dX(:,:,:,:,:,i), boxes(ibfar)%outgoing0)

      elseif (order == 2) then

        boxes(ib)%incoming0 = boxes(ib)%incoming0 + TijklTkl(boxes(ib)%Gamma(:,:,:,:,i), boxes(ibfar)%outgoing0) - &
          TijklmTklm(boxes(ib)%dGamma_dX(:,:,:,:,:,i), boxes(ibfar)%outgoing1) + &
          0.5*TijklmnTklmn(boxes(ib)%d2Gamma_dX2(:,:,:,:,:,:,i), boxes(ibfar)%outgoing2)
  
        boxes(ib)%incoming1 = boxes(ib)%incoming1 + TijklmTkl(boxes(ib)%dGamma_dX(:,:,:,:,:,i), boxes(ibfar)%outgoing0) - &
          0.5*TijklmnTkln(boxes(ib)%d2Gamma_dX2(:,:,:,:,:,:,i), boxes(ibfar)%outgoing1) - &
          0.5*TijklmnTklm(boxes(ib)%d2Gamma_dX2(:,:,:,:,:,:,i), boxes(ibfar)%outgoing1)
  
        boxes(ib)%incoming2 = boxes(ib)%incoming2 + &
          0.5*TijklmnTkl(boxes(ib)%d2Gamma_dX2(:,:,:,:,:,:,i), boxes(ibfar)%outgoing0)

      endif

    enddo
    
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine incoming_expansion_FFT(boxes)
  use global
  use types
  use tensor_functions
  use fourier_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: ib, i, ib1, ib2, ib3, j, k, l, m
  double complex :: aux33c(3,3), dumc, aux333c(3,3,3), outgoing0(3,3), outgoing1(3,3,3)
  double complex :: outgoing2(3,3,3,3), aux3333c(3,3,3,3)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,i,j,k,l,ib1,ib2,ib3) &
  !$OMP SHARED(nbox,boxes,order,fourgrid33box,fourgrid333box,fourgrid3333box)
  !$OMP DO 
  do ib = 1,nbox
    ib1 = boxes(ib)%ib_grid(1)
    ib2 = boxes(ib)%ib_grid(2)
    ib3 = boxes(ib)%ib_grid(3)
    do i = 1,3
      do j = 1,3
        fourgrid33box(ib3,ib2,ib1,j,i) = cmplx(boxes(ib)%outgoing0(i,j))
        if (order > 0) then
          do k = 1,3
            fourgrid333box(ib3,ib2,ib1,k,j,i) = cmplx(boxes(ib)%outgoing1(i,j,k))
            if (order > 1) then
              do l = 1,3
                fourgrid3333box(ib3,ib2,ib1,l,k,j,i) = cmplx(boxes(ib)%outgoing2(i,j,k,l))
              enddo
            endif
          enddo
        endif
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call fft_tensor33box(plan_advanced33box, fourgrid33box)
  if (order > 0) call fft_tensor333box(plan_advanced333box, fourgrid333box)
  if (order > 1) call fft_tensor3333box(plan_advanced3333box, fourgrid3333box)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,outgoing0,outgoing1,outgoing2,aux33c,aux333c,aux3333c) &
  !$OMP SHARED(nbox1,nbox2,nbox3,boxes,order,fourgrid33box,fourgrid333box,fourgrid3333box) &
  !$OMP SHARED(Gamma_box_hat,dGamma_dX_box_hat,d2Gamma_dX2_box_hat)
  !$OMP DO COLLAPSE(3)
  do ib1 = 1,nbox1
    do ib2 = 1,nbox2
      do ib3 = 1,nbox3

        ! do i = 1,3
        !   do j = 1,3
        !     dumc = (0.0,0.0)
        !     do k = 1,3
        !       do l = 1,3
        !         dumc = dumc + fourgrid3333box(ib3,ib2,ib1,l,k,j,i)*fourgrid33box(ib3,ib2,ib1,l,k)
        !       enddo
        !     enddo
        !     aux33c(j,i) = dumc
        !   enddo
        ! enddo
        ! do i = 1,3
        !   do j = 1,3
        !     fourgrid33box(ib3,ib2,ib1,j,i) = aux33c(j,i)
        !   enddo
        ! enddo

        if (order == 0) then

          outgoing0 = fourgrid33box(ib3,ib2,ib1,:,:)

          aux33c = TlkjiTlk_cmplx(Gamma_box_hat(ib3,ib2,ib1,:,:,:,:),outgoing0)
          fourgrid33box(ib3,ib2,ib1,:,:) = aux33c

        elseif (order == 1) then

          outgoing0 = fourgrid33box(ib3,ib2,ib1,:,:)
          outgoing1 = fourgrid333box(ib3,ib2,ib1,:,:,:)

          aux33c = TlkjiTlk_cmplx(Gamma_box_hat(ib3,ib2,ib1,:,:,:,:),outgoing0) - &
           TmlkjiTmlk_cmplx(dGamma_dX_box_hat(ib3,ib2,ib1,:,:,:,:,:),outgoing1)
          fourgrid33box(ib3,ib2,ib1,:,:) = aux33c

          aux333c = TmlkjiTlk_cmplx(dGamma_dX_box_hat(ib3,ib2,ib1,:,:,:,:,:),outgoing0)
          fourgrid333box(ib3,ib2,ib1,:,:,:) = aux333c

        elseif (order == 2) then

          outgoing0 = fourgrid33box(ib3,ib2,ib1,:,:)
          outgoing1 = fourgrid333box(ib3,ib2,ib1,:,:,:)
          outgoing2 = fourgrid3333box(ib3,ib2,ib1,:,:,:,:)
  
          aux33c = TlkjiTlk_cmplx(Gamma_box_hat(ib3,ib2,ib1,:,:,:,:),outgoing0) - &
           TmlkjiTmlk_cmplx(dGamma_dX_box_hat(ib3,ib2,ib1,:,:,:,:,:),outgoing1) + &
           0.5*TnmlkjiTnmlk_cmplx(d2Gamma_dX2_box_hat(ib3,ib2,ib1,:,:,:,:,:,:),outgoing2)
          fourgrid33box(ib3,ib2,ib1,:,:) = aux33c
  
          aux333c = TmlkjiTlk_cmplx(dGamma_dX_box_hat(ib3,ib2,ib1,:,:,:,:,:),outgoing0) - &
           0.5*TnmlkjiTnlk_cmplx(d2Gamma_dX2_box_hat(ib3,ib2,ib1,:,:,:,:,:,:), outgoing1) - &
           0.5*TnmlkjiTmlk_cmplx(d2Gamma_dX2_box_hat(ib3,ib2,ib1,:,:,:,:,:,:), outgoing1)
          fourgrid333box(ib3,ib2,ib1,:,:,:) = aux333c

          aux3333c = 0.5*TnmlkjiTlk_cmplx(d2Gamma_dX2_box_hat(ib3,ib2,ib1,:,:,:,:,:,:), outgoing0)
          fourgrid3333box(ib3,ib2,ib1,:,:,:,:) = aux3333c

        endif

      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call ifft_tensor33box(iplan_advanced33box, fourgrid33box)
  if (order > 0) call ifft_tensor333box(iplan_advanced333box, fourgrid333box)
  if (order > 1) call ifft_tensor3333box(iplan_advanced3333box, fourgrid3333box)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,i,j,k,l,ib1,ib2,ib3) &
  !$OMP SHARED(nbox,boxes,order,fourgrid33box,fourgrid333box,fourgrid3333box)
  !$OMP DO 
  do ib = 1,nbox
    ib1 = boxes(ib)%ib_grid(1)
    ib2 = boxes(ib)%ib_grid(2)
    ib3 = boxes(ib)%ib_grid(3)
    ! write(*,*) ib
    ! write(*,*) boxes(ib)%incoming0(:,:)
    do i = 1,3
      do j = 1,3
        boxes(ib)%incoming0(i,j) = real(fourgrid33box(ib3,ib2,ib1,j,i))
        if (order > 0) then
          do k = 1,3
            boxes(ib)%incoming1(i,j,k) = real(fourgrid333box(ib3,ib2,ib1,k,j,i))
            if (order > 1) then
              do l = 1,3
                boxes(ib)%incoming2(i,j,k,l) = real(fourgrid3333box(ib3,ib2,ib1,l,k,j,i))
              enddo
            endif
          enddo
        endif
      enddo
    enddo
    ! write(*,*) boxes(ib)%incoming0(:,:)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine convolution(elements, boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (in) :: boxes(:)

  double precision :: strainacc(3,3)
  integer :: iel1, i, iel2, ib

  strainavg = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel1, iel2, i, ib, strainacc) &
  !$OMP SHARED(elements, boxes, nel, Eapp, ifull, Gmat, iJacobi, order) &
  !$OMP REDUCTION(+:strainavg)
  !$OMP DO 
  do iel1 = 1,nel
    if (elements(iel1)%ngrid_ind > 0) then
      ib = elements(iel1)%box_id
   
      if (ifull == 1) then
        strainacc = 0.0
        do iel2 = 1,nel
          if (iel2 .ne. iel1) then
            strainacc = strainacc + TijklTkl(Gmat(:,:,:,:,iel1,iel2),elements(iel2)%polar(:,:))
          endif
        enddo
        elements(iel1)%strain = Eapp + strainacc
        if (iJacobi == 1) then
          elements(iel1)%strain = TijklTkl(elements(iel1)%I_GdC_inv,elements(iel1)%strain)
        else
          elements(iel1)%strain = elements(iel1)%strain + TijklTkl(elements(iel1)%Gamma_self,elements(iel1)%polar)
        endif
      elseif (ifull == 0) then
        strainacc = 0.0 
        do i = 1,elements(iel1)%nel_close
          iel2 = elements(iel1)%close_inter_el_ids(i)
          strainacc = strainacc + TijklTkl(elements(iel1)%Gamma(:,:,:,:,i),elements(iel2)%polar(:,:))
        enddo
   
        if (order == 0) then
         elements(iel1)%strain = Eapp + strainacc + boxes(ib)%incoming0
        elseif (order == 1) then
         elements(iel1)%strain = Eapp + strainacc + boxes(ib)%incoming0 + &
          TijkTk(boxes(ib)%incoming1, elements(iel1)%dx)
        elseif (order == 2) then
         elements(iel1)%strain = Eapp + strainacc + boxes(ib)%incoming0 + &
          TijkTk(boxes(ib)%incoming1, elements(iel1)%dx) + &
          TijklTkl(boxes(ib)%incoming2, elements(iel1)%dxdx)
        endif
   
        if (iJacobi == 1) then
          elements(iel1)%strain = TijklTkl(elements(iel1)%I_GdC_inv,elements(iel1)%strain)
        else
          elements(iel1)%strain = elements(iel1)%strain + TijklTkl(elements(iel1)%Gamma_self,elements(iel1)%polar)
        endif
      endif
      strainavg = strainavg + elements(iel1)%strain(:,:)*elements(iel1)%v
    endif
 
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  strainavg = strainavg/vtot

  ! do iel1 = 1,nel
  !  elements(iel1)%strain = elements(iel1)%strain - strainavg + Eapp
  ! enddo

end

subroutine box_far_inter_Gamma_old(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dx(3)
  integer :: ib, i, ibfar, inc, ip1, ip2, ip3, ip1p, ip2p, ip3p, ip1m, ip2m, ip3m

  do ib = 1,nbox
    allocate(boxes(ib)%Gamma(3,3,3,3,boxes(ib)%nfar))
    if (order > 0) allocate(boxes(ib)%dGamma_dX(3,3,3,3,3,boxes(ib)%nfar))
    if (order > 1) allocate(boxes(ib)%d2Gamma_dX2(3,3,3,3,3,3,boxes(ib)%nfar))
  enddo

  ! far field interactions
  inc = 1
  do ib = 1,nbox
    do i = 1,boxes(ib)%nfar
      ibfar = boxes(ib)%far_box_id(i)
      dx = - boxes(ibfar)%x + boxes(ib)%x
  
      ip1 = enforce_periodic(nint(dx(1)) + 1, 1, npts1)
      ip2 = enforce_periodic(nint(dx(2)) + 1, 1, npts2)
      ip3 = enforce_periodic(nint(dx(3)) + 1, 1, npts3)
  
      ip1p = enforce_periodic(ip1 + inc, 1, npts1)
      ip1m = enforce_periodic(ip1 - inc, 1, npts1)
      ip2p = enforce_periodic(ip2 + inc, 1, npts2)
      ip2m = enforce_periodic(ip2 - inc, 1, npts2)
      ip3p = enforce_periodic(ip3 + inc, 1, npts3)
      ip3m = enforce_periodic(ip3 - inc, 1, npts3)
  
      boxes(ib)%Gamma(:,:,:,:,i) = Goperreal(:,:,:,:,ip1,ip2,ip3)

      if (order > 0) then
        boxes(ib)%dGamma_dX(:,:,:,:,:,i) = 0.0
        boxes(ib)%dGamma_dX(:,:,:,:,1,i) = (Goperreal(:,:,:,:,ip1p,ip2,ip3) - Goperreal(:,:,:,:,ip1m,ip2,ip3))/float(2*inc)
        boxes(ib)%dGamma_dX(:,:,:,:,2,i) = (Goperreal(:,:,:,:,ip1,ip2p,ip3) - Goperreal(:,:,:,:,ip1,ip2m,ip3))/float(2*inc)
        boxes(ib)%dGamma_dX(:,:,:,:,3,i) = (Goperreal(:,:,:,:,ip1,ip2,ip3p) - Goperreal(:,:,:,:,ip1,ip2,ip3m))/float(2*inc)
      endif

      if (order > 1) then
        boxes(ib)%d2Gamma_dX2(:,:,:,:,:,:,i) = 0.0
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,1,i) = (Goperreal(:,:,:,:,ip1p,ip2,ip3) - 2.0*Goperreal(:,:,:,:,ip1,ip2,ip3) + &
         Goperreal(:,:,:,:,ip1m,ip2,ip3))/float(inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,2,2,i) = (Goperreal(:,:,:,:,ip1,ip2p,ip3) - 2.0*Goperreal(:,:,:,:,ip1,ip2,ip3) + &
         Goperreal(:,:,:,:,ip1,ip2m,ip3))/float(inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,3,i) = (Goperreal(:,:,:,:,ip1,ip2,ip3p) - 2.0*Goperreal(:,:,:,:,ip1,ip2,ip3) + &
         Goperreal(:,:,:,:,ip1,ip2,ip3m))/float(inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,2,i) = (Goperreal(:,:,:,:,ip1p,ip2p,ip3) - Goperreal(:,:,:,:,ip1p,ip2m,ip3) - &
         Goperreal(:,:,:,:,ip1m,ip2p,ip3) + Goperreal(:,:,:,:,ip1m,ip2m,ip3))/float(4*inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,1,3,i) = (Goperreal(:,:,:,:,ip1p,ip2,ip3p) - Goperreal(:,:,:,:,ip1p,ip2,ip3m) - &
         Goperreal(:,:,:,:,ip1m,ip2,ip3p) + Goperreal(:,:,:,:,ip1m,ip2,ip3m))/float(4*inc**2)
         boxes(ib)%d2Gamma_dX2(:,:,:,:,2,3,i) = (Goperreal(:,:,:,:,ip1,ip2p,ip3p) - Goperreal(:,:,:,:,ip1,ip2p,ip3m) - &
         Goperreal(:,:,:,:,ip1,ip2m,ip3p) + Goperreal(:,:,:,:,ip1,ip2m,ip3m))/float(4*inc**2)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,2,1,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,1,2,i)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,1,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,1,3,i)
        boxes(ib)%d2Gamma_dX2(:,:,:,:,3,2,i) = boxes(ib)%d2Gamma_dX2(:,:,:,:,2,3,i)
      endif
  
    enddo
  enddo

end

subroutine integrate_G_Taylor(iel1,iel2,elements,boxes,G)
  use global
  use types
  use tensor_functions

  integer, intent (in) :: iel1,iel2
  type (element_type), allocatable, intent (in) :: elements(:)
  type (box_type), allocatable, intent (in) :: boxes(:)
  double precision, intent (out) :: G(3,3,3,3)

  integer :: ib1, ib2, i, j, k, ib

  G = 0.0
  ib1 = elements(iel1)%box_id
  ib2 = elements(iel2)%box_id
  do i = 1,boxes(ib1)%nfar
    ib = boxes(ib1)%far_box_id(i)
    if (ib == ib2) then
      G = boxes(ib1)%Gamma(:,:,:,:,i)*elements(iel2)%v
      do j = 1,3
        G = G + boxes(ib1)%dGamma_dX(:,:,:,:,j,i)*(elements(iel1)%dx(j) - elements(iel2)%dx(j))*elements(iel2)%v
      enddo

      do j = 1,3
        do k = 1,3
          G = G + boxes(ib1)%d2Gamma_dX2(:,:,:,:,j,k,i)* &
           (elements(iel1)%dxdx(j,k) + elements(iel2)%dxdx(j,k) - &
            elements(iel1)%dx(j)*elements(iel2)%dx(k) - elements(iel2)%dx(j)*elements(iel1)%dx(k))*elements(iel2)%v
        enddo
      enddo
      exit
    endif
  enddo

end

function tensor4_to_tensor2(A4) result(A2)
  implicit none
  double precision, intent (in) :: A4(3,3,3,3)

  double precision :: A2(9,9)
  integer :: i, j, k, l, m, n

  m = 0
  do i = 1,3
    do j = 1,3
      m = m + 1
      n = 0
      do k = 1,3
        do l = 1,3
          n = n + 1
          A2(m,n) = A4(i,j,k,l)
        enddo
      enddo
    enddo
  enddo

end

function tensor2_to_tensor4(A2) result(A4)
  implicit none
  double precision, intent (in) :: A2(9,9)

  double precision :: A4(3,3,3,3)
  integer :: i, j, k, l, m, n

  m = 0
  do i = 1,3
    do j = 1,3
      m = m + 1
      n = 0
      do k = 1,3
        do l = 1,3
          n = n + 1
          A4(i,j,k,l) = A2(m,n)
        enddo
      enddo
    enddo
  enddo

end

subroutine calc_I_GdC_inv(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: id4(3,3,3,3), I_GdC(3,3,3,3), I_GdC_inv99(9,9)
  integer :: i, iel, j, k, l

  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          id4(i,j,k,l) = id3(i,k)*id3(j,l)
        enddo
      enddo
    enddo
  enddo

  do iel = 1,nel

    I_GdC = id4 - TijklTklmn(elements(iel)%Gamma_self, elements(iel)%c - c0)
    I_GdC_inv99 = tensor4_to_tensor2(I_GdC)
    call lu_inverse(I_GdC_inv99,9)
    elements(iel)%I_GdC_inv = tensor2_to_tensor4(I_GdC_inv99)

  enddo

end

subroutine form_G_approx(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: G(3,3,3,3),x_ref(3),xtilde_ref(3), y(3), ytilde(3), dx(3)
  integer :: i, ib, ibfar, j, k, l, ib_ref, ip1, ip2, ip3

  Goperreal_approx = Goperreal

  ib_ref = 1
  x_ref = 1.0
  xtilde_ref = x_ref - boxes(ib_ref)%x
  do ib = 1,nbox
    do i = 1,boxes(ib)%nfar

      ibfar = boxes(ib)%far_box_id(i)

      if (ibfar == ib_ref) then ! reference is in far field

        do ip1 = boxes(ib)%ipstart(1),boxes(ib)%ipend(1)
          do ip2 = boxes(ib)%ipstart(2),boxes(ib)%ipend(2)
            do ip3 = boxes(ib)%ipstart(3),boxes(ib)%ipend(3)

              y = float([ip1,ip2,ip3])
              ytilde = y - boxes(ib)%x

              dx = ytilde - xtilde_ref ! for Taylor expansion

              G = boxes(ib)%Gamma(:,:,:,:,i)
              Goperreal_approx(:,:,:,:,ip1,ip2,ip3) = G

              if (order > 0) then
                G = 0.0
                do j = 1,3
                  G = G + boxes(ib)%dGamma_dX(:,:,:,:,j,i)*dx(j)
                enddo
                Goperreal_approx(:,:,:,:,ip1,ip2,ip3) = Goperreal_approx(:,:,:,:,ip1,ip2,ip3) + G
              endif

              if (order > 1) then
                G = 0.0
                do j = 1,3
                  do k = 1,3
                    G = G + boxes(ib)%d2Gamma_dX2(:,:,:,:,j,k,i)*dx(j)*dx(k)*0.5
                  enddo
                enddo
                Goperreal_approx(:,:,:,:,ip1,ip2,ip3) = Goperreal_approx(:,:,:,:,ip1,ip2,ip3) + G
              endif

            enddo
          enddo
        enddo

      
      endif

    enddo
  enddo

end

end module various_functions
