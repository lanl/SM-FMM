module IO_functions
  implicit none

contains


  subroutine load_input(input_file, elements, grains, boxes)
    use global
    use fourier_functions
    use types
    use various_functions

    character(len = 500), intent (in) :: input_file
    type (element_type), allocatable, intent (inout) :: elements(:)
    type (grain_type), allocatable, intent (inout) :: grains(:)
    type (box_type), allocatable, intent (inout) :: boxes(:)

    character(len = 500) :: mesh_file, tex_file, elas_file
    integer :: iFFT

    open (unit = 11, file = input_file, status = 'old')

    read (11,'(A)') mesh_file
    read (11,'(A)') tex_file
    read (11,*) elas_file
    read (11,*) npts1, npts2, npts3
    read (11,*) nbox1, nbox2, nbox3
    read (11,*) 
    read (11,*) Eapp(1,1), Eapp(1,2), Eapp(1,3)
    read (11,*) Eapp(2,1), Eapp(2,2), Eapp(2,3)
    read (11,*) Eapp(3,1), Eapp(3,2), Eapp(3,3)
    read (11,*) ifull
    read (11,*) iload
    read (11,*) iJacobi
    read (11,*) xc0
    read (11,*) step
    read (11,*) order
    read (11,*) iFFT

    close (11)

    call allocate_globals

    ! call load_vtk(mesh_file, elements)
    call load_micro(mesh_file, elements)
    if (ifull == 1) allocate(Gmat(3,3,3,3,nel,nel))
    if (ifull == 1) Gmat = 0
    call data_crystal_elast(elas_file)
    call data_grain(tex_file,grains)

    !$ call setup_openmp()
    call initialize_fftw(npts1,npts2,npts3,nbox1,nbox2,nbox3,nthreads)

    write(*,*) 'form_G start'
    call form_G
    call form_G_box
    write(*,*) 'form_G end'

    ! test box for same nbox as npts
    ! write(*,*) sqrt(sum((Goperreal-Goperreal_box)**2)) 

    ! call enforce_domain(elements)
    ! call init_integration_pt(elements)

    nbox = nbox1*nbox2*nbox3
    allocate(boxes(nbox))

    if (iFFT == 0) then
      FFT = .false.
    elseif (iFFT == 1) then
      FFT = .true.
    endif

  end

  subroutine load_vtk(mesh_file, elements)
    use global
    use types

    character(len = 500), intent (in) :: mesh_file
    type (element_type), allocatable, intent (inout) :: elements(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel
    integer :: nh_lines = 4
    double precision, allocatable :: x_node(:,:)
    logical, allocatable :: tetra(:)

    open (unit = 33, file = mesh_file, status = 'old')
     
    do i = 1,nh_lines
      read(33,'(A)') line
    enddo
    read(33,*) str1, nnode, str2
    if (trim(str1) .ne. 'POINTS') stop 'trim(str1) .ne. POINTS'

    allocate(x_node(3,nnode))
    do i = 1,nnode

      read(33,*) x_node(1,i), x_node(2,i), x_node(3,i)
      x_node(1,i) = x_node(1,i)*float(npts1) + 0.5 ! position is in domain 0.5 - npts + 0.5 (like in standard FFT)
      x_node(2,i) = x_node(2,i)*float(npts2) + 0.5
      x_node(3,i) = x_node(3,i)*float(npts3) + 0.5

    enddo

    read(33,'(A)')line

    read(33,*) str1, neltot, i
    if (trim(str1) .ne. 'CELLS') stop 'trim(str1) .ne. CELLS'

    allocate(ind(20,neltot))
    allocate(tetra(neltot))
    ind = 0
    do i = 1,neltot

      read(33,*) icell, (ind(j,i), j = 1,icell)

    enddo

    read(33,'(A)')line

    read(33,*) str1, neltot
    if (trim(str1) .ne. 'CELL_TYPES') stop 'trim(str1) .ne. CELL_TYPES'

    nel = 0
    do i = 1,neltot

      read(33,*) icell
      if (icell == 10) then ! tetra
        nel = nel + 1
        tetra(i) = .true.
      else
        tetra(i) = .false.
      endif

    enddo

    allocate(elements(nel))
    allocate(plot_el_scalar(nel))
    i = 0
    do iel = 1,neltot

      if (tetra(iel)) then
        i = i + 1
        do j = 1,4
          elements(i)%xnode(:,j) = x_node(:,ind(j,iel) + 1)
          elements(i)%ind_node(j) = ind(j,iel) + 1
        enddo
      endif
      
    enddo

    read(33,'(A)')line
    read(33,'(A)') line
    read(33,'(A)') line
    read(33,'(A)') line

    i = 0
    do iel = 1,neltot

      if (tetra(iel)) then
        i = i + 1
        read(33,*) elements(i)%grain_id
      else
        read(33,*)
      endif
      
    enddo

    deallocate(x_node)
    deallocate(ind)
    deallocate(tetra)

    close (33)

  end


  subroutine load_micro(mesh_file, elements)
    use global
    use types

    character(len = 500), intent (in) :: mesh_file
    type (element_type), allocatable, intent (inout) :: elements(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel, ip1, ip2, ip3, idum
    integer :: nh_lines = 4
    double precision, allocatable :: x_node(:,:)
    double precision :: dum
    logical, allocatable :: tetra(:)

    open (unit = 33, file = mesh_file, status = 'old')

    nel = -1
    do ip3 = 1,npts3
      do ip2 = 1,npts2
        do ip1 = 1,npts1
          read(33,*) dum,dum,dum,idum,idum,idum,el_id(ip1,ip2,ip3),idum
          if (el_id(ip1,ip2,ip3) > nel) nel = el_id(ip1,ip2,ip3)
        enddo
      enddo
    enddo

    allocate(elements(nel))
    do iel = 1,nel
      elements(iel)%grain_id = iel
    enddo
    write(*,*) nel

    close (33)

  end

  subroutine write_vtk(vtk_file, elements)
    use global
    use types

    character(len = 500), intent (in) :: vtk_file
    type (element_type), allocatable, intent (inout) :: elements(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel, in, number_fields, ic
    integer :: nh_lines = 4
    double precision, allocatable :: x_node(:,:)
    double precision :: dxdx(3,3)
    logical, allocatable :: tetra(:)

    open (unit = 33, file = vtk_file, status = 'unknown')
     
    write(33,'(A)')'# vtk DataFile Version 3.0'
    write(33,'(A)')'microstructure data'
    write(33,'(A)')'ASCII'
    write(33,'(A)')'DATASET UNSTRUCTURED_GRID'

    ! allocate(x_node(3,nnode))
    allocate(x_node(3,nel*4))

    ic = 0
    do iel = 1, nel

      do i = 1,4
        ic = ic + 1
        j = elements(iel)%ind_node(i)

        ! x_node(:,j) = elements(iel)%xnode(:,i)
        x_node(:,ic) = elements(iel)%xnode(:,i)
        elements(iel)%ind_node(i) = ic
      enddo

    enddo

    ! write(33,'(A,I10,A)')'POINTS',nnode,' double'
    write(33,'(A,I10,A)')'POINTS',ic,' double'
    do in = 1, ic ! nnode
      write(33,'(3E18.8E3)') x_node(:,in)
    enddo

    deallocate(x_node)

    write(33,'(A,I10,I10)')'CELLS',nel,nel*5
    do iel = 1, nel
      write(33,'(5I10)')4, elements(iel)%ind_node - 1
    enddo

    write(33,'(A,I10)')'CELL_TYPES',nel
    do iel = 1, nel
      write(33,'(I5)') 10 
    enddo

    number_fields = 6
    write(33,'(A,I10)')'CELL_DATA',nel
    write(33,'(A,I10)')'FIELD FieldData',number_fields

    write(33,'(A,I3,I10,A)')'grain_id',1,nel,' int'
    do iel = 1, nel
      write(33,'(I10)') elements(iel)%grain_id
    enddo
    write(33,'(A,I3,I10,A)')'plot_el_scalar',1,nel,' double'
    do iel = 1, nel
      write(33,'(100E18.8E3)') plot_el_scalar(iel)
    enddo
    write(33,'(A,I3,I10,A)')'el_id',1,nel,' int'
    do iel = 1, nel
      write(33,'(I10)') iel
    enddo

    write(33,'(A,I3,I10,A)')'stress',6,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%stress(1,1),elements(iel)%stress(2,2), &
      elements(iel)%stress(3,3),elements(iel)%stress(1,2),elements(iel)%stress(2,3),&
      elements(iel)%stress(1,3)
    enddo

    write(33,'(A,I3,I10,A)')'dx',3,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%dx(1),elements(iel)%dx(2), &
      elements(iel)%dx(3)
    enddo

    write(33,'(A,I3,I10,A)')'dxdx-dx*dx',9,nel,' double'
    do iel = 1, nel
      ! write(33,'(6E18.8E3)')elements(iel)%dxdx(1,1),elements(iel)%dxdx(2,1), &
      ! elements(iel)%dxdx(3,1),elements(iel)%dxdx(1,2),elements(iel)%dxdx(2,2),&
      ! elements(iel)%dxdx(3,2),elements(iel)%dxdx(1,3),elements(iel)%dxdx(2,3),&
      ! elements(iel)%dxdx(3,3)
      do i = 1,3
        do j = 1,3
          dxdx(i,j) = elements(iel)%dxdx(i,j) - elements(iel)%dx(i)*elements(iel)%dx(j)
        enddo
      enddo
      write(33,'(6E18.8E3)')dxdx(1,1),dxdx(2,1), &
      dxdx(3,1),dxdx(1,2),dxdx(2,2),&
      dxdx(3,2),dxdx(1,3),dxdx(2,3),&
      dxdx(3,3)

    enddo

    close (33)

  end

  subroutine write_grid_vtk(vtk_file, elements)
    use global
    use types

    character(len = 500), intent (in) :: vtk_file
    type (element_type), allocatable, intent (in) :: elements(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel, in, number_fields, ic, ip1, ip2, ip3
    integer :: nh_lines = 4
    double precision, allocatable :: x_node(:,:)
    logical, allocatable :: tetra(:)

    open (unit = 33, file = vtk_file, status = 'unknown')
     
    write(33,'(A)')'# vtk DataFile Version 3.0'
    write(33,'(A)')'microstructure data'
    write(33,'(A)')'ASCII'
    write(33,'(A)')'DATASET STRUCTURED_POINTS'
    write(33,'(A,3I10)')'DIMENSIONS', npts1+1, npts2+1, npts3+1
    write(33,'(A,3E18.8E3)')'ORIGIN', 0.5, 0.5, 0.5
    write(33,'(A,3E18.8E3)')'SPACING', 1.0, 1.0, 1.0

    number_fields = 3
    write(33,'(A,I10)')'CELL_DATA',npts1*npts2*npts3
    write(33,'(A,I10)')'FIELD FieldData',number_fields

    write(33,'(A,I3,I10,A)')'el_id',1,npts1*npts2*npts3,' int'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          write(33,'(I10)') el_id(ip1,ip2,ip3)
        enddo
      enddo
    enddo

    write(33,'(A,I3,I10,A)')'stress',6,npts1*npts2*npts3,' double'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          iel = el_id(ip1,ip2,ip3)
          write(33,'(6E18.8E3)')elements(iel)%stress(1,1),elements(iel)%stress(2,2), &
          elements(iel)%stress(3,3),elements(iel)%stress(1,2),elements(iel)%stress(2,3),&
          elements(iel)%stress(1,3)
        enddo
      enddo
    enddo

    write(33,'(A,I3,I10,A)')'strain',6,npts1*npts2*npts3,' double'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          iel = el_id(ip1,ip2,ip3)
          write(33,'(6E18.8E3)')elements(iel)%strain(1,1),elements(iel)%strain(2,2), &
          elements(iel)%strain(3,3),elements(iel)%strain(1,2),elements(iel)%strain(2,3),&
          elements(iel)%strain(1,3)
        enddo
      enddo
    enddo

!    write(33,'(A,I3,I10,A)')'G_approx',2,npts1*npts2*npts3,' double'
!    do ip3 = 1, npts3
!      do ip2 = 1, npts2
!        do ip1 = 1, npts1
!          write(33,'(813E18.8E3)') Goperreal_approx(1,1,1,1,ip1,ip2,ip3),Goperreal_approx(1,2,1,2,ip1,ip2,ip3)
!        enddo
!      enddo
!    enddo
!
!    write(33,'(A,I3,I10,A)')'G',2,npts1*npts2*npts3,' double'
!    do ip3 = 1, npts3
!      do ip2 = 1, npts2
!        do ip1 = 1, npts1
!          write(33,'(813E18.8E3)') Goperreal(1,1,1,1,ip1,ip2,ip3),Goperreal(1,2,1,2,ip1,ip2,ip3)
!        enddo
!      enddo
!    enddo
!
!    write(33,'(A,I3,I10,A)')'point',3,npts1*npts2*npts3,' int'
!    do ip3 = 1, npts3
!      do ip2 = 1, npts2
!        do ip1 = 1, npts1
!          write(33,'(3I10)') ip1,ip2,ip3
!        enddo
!      enddo
!    enddo

    close (33)

  end

  subroutine data_crystal_elast(elas_file)
    use tensor_functions
    use global
    implicit none
  
    character(len = 500), intent (in) :: elas_file

    double precision :: dde(3,3), xid4(3,3,3,3)
    double precision :: cc66v(6,6), ccaux(3,3,3,3)
    double precision :: aux6(6), aux33(3,3)
  
    integer :: i, iso, j, k, l
    double precision :: tla, tmu, tnu, young
  
    open (unit = 33, file = elas_file, status = 'old')

    ! unitary tensors
    do i = 1,3
    do j = 1,3
      dde(i,j) = 0.d0
      if (i.eq.j) dde(i,j) = 1.d0
    enddo
    enddo
  
    do i = 1,3
    do j = 1,3
    do k = 1,3
    do l = 1,3
      xid4(i,j,k,l) = (dde(i,k)*dde(j,l) + dde(i,l)*dde(j,k)) / 2.d0
    enddo
    enddo
    enddo
    enddo
  
    read(33,*) iso
  
    if(iso.eq.0) then
  
      do i=1,6
        read(33,*)(cc66v(i,j), j=1,6)
      enddo
  
      call voigt_vpsc(aux6,aux33,cc66v,ccaux,3)
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
        C(i,j,k,l) = ccaux(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
  
    else
  
      read(33,*) young, tnu
      tmu = young / (2.d0*(1.+tnu))
      tla = 2.d0*tmu*tnu / (1.d0 - 2.d0*tnu)
  
      do i = 1,3
      do j = 1,3
      do k = 1,3
      do l = 1,3
        C(i,j,k,l) = tla*dde(i,j)*dde(k,l) + 2.d0*tmu*xid4(i,j,k,l)
      enddo
      enddo
      enddo
      enddo
  
    endif

    close (33)
  
    return
  end
  
  subroutine data_grain(tex_file,grains)
    use global
   use tensor_functions
   use types
   implicit none
  
   character(len = 500), intent (in) :: tex_file
   type (grain_type), allocatable, intent (inout) :: grains(:)

   double precision :: aa(3,3)
   double precision :: caux3333(3,3,3,3),caux66(6,6)
   double precision :: aux6(6),aux33(3,3)
  
   double precision :: dum, om, ph, th
   integer :: i1, i2, ii, j, j1, j2, jgr, jj, jph
   integer :: k, k1, k2, kk, kkk, l1, l2, nph1
  
   open (unit = 33, file = tex_file, status = 'old')
   
   nph1 = 0
   read(33,*) ng
   allocate(grains(ng))

   c066 = 0.0
   do kkk = 1,ng
  
     read(33,*) ph,th,om
  
     grains(kkk)%id = kkk
     grains(kkk)%phase = 1
     grains(kkk)%eul = [ph,th,om]
  
     ! calculates the transformation matrix aa which transforms from
     ! sample to crystal. stores ag, which transforms from crystal to sample.
     call euler(2,ph*pi/180.0,th*pi/180.0,om*pi/180.0,aa)
  
     do j=1,3
     do k=1,3
      grains(kkk)%Q(j,k) = aa(k,j)
     enddo
     enddo
  
     do i1=1,3
     do j1=1,3
     do k1=1,3
     do l1=1,3
       dum=0.
       do i2=1,3
       do j2=1,3
       do k2=1,3
       do l2=1,3
         dum=dum+aa(i2,i1)*aa(j2,j1)*aa(k2,k1)*aa(l2,l1)*C(i2,j2,k2,l2)
       enddo
       enddo
       enddo
       enddo
       grains(kkk)%C(i1,j1,k1,l1) = dum
     enddo
     enddo
     enddo
     enddo
  
     call chg_basis(aux6,aux33,grains(kkk)%C66,grains(kkk)%C,4,6)

     c066 = c066 + grains(kkk)%C66
  
   enddo
   c066 = c066/float(ng)*xc0
   s066 = c066
   call lu_inverse(s066,6)
   call chg_basis(aux6,aux33,c066,c0,3,6)
   call chg_basis(aux6,aux33,s066,s0,3,6)

   write(*,*) 'c066 = '
   write(*,*) c066

   close (33)

   return
  end

  subroutine enforce_domain(elements)
    use global
    use types
    use tensor_functions

    type (element_type), allocatable, intent (inout) :: elements(:)
    
    integer :: iel, idir
    double precision, parameter ::  one_sixth = 0.166666666666667
    double precision :: dir_range(3)

    dir_range(1) = float(npts1)
    dir_range(2) = float(npts2)
    dir_range(3) = float(npts3)
    
    do iel = 1,nel

      elements(iel)%x = &
        elements(iel)%xnode(:,1)*0.25 + elements(iel)%xnode(:,2)*0.25 + &
        elements(iel)%xnode(:,3)*0.25 + elements(iel)%xnode(:,4)*0.25

      ! do idir = 1,3
      !   if (elements(iel)%x(idir) <= 0.0) then
      !     elements(iel)%x(idir) = elements(iel)%x(idir) + 1.0
      !     elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) + 1.0
      !   elseif (elements(iel)%x(idir) > 1.0) then
      !     elements(iel)%x(idir) = elements(iel)%x(idir) - 1.0
      !     elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) - 1.0
      !   endif
      ! enddo

      do idir = 1,3
        if (elements(iel)%x(idir) <= 0.5) then
          elements(iel)%x(idir) = elements(iel)%x(idir) + dir_range(idir)
          elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) + dir_range(idir)
          ! if (elements(iel)%x(idir) <= 1.0) then
          !   elements(iel)%x(idir) = elements(iel)%x(idir) + dir_range(idir)
          !   elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) + dir_range(idir)
          ! endif
        elseif (elements(iel)%x(idir) > dir_range(idir) + 0.5) then
          elements(iel)%x(idir) = elements(iel)%x(idir) - dir_range(idir)
          elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) - dir_range(idir)
          ! if (elements(iel)%x(idir) > dir_range(idir)) then
          !   elements(iel)%x(idir) = elements(iel)%x(idir) - dir_range(idir)
          !   elements(iel)%xnode(idir,:) = elements(iel)%xnode(idir,:) - dir_range(idir)
          ! endif
        endif
      enddo

    enddo

  end

  subroutine init_integration_pt(elements)
    use global
    use types
    use tensor_functions
    use various_functions

    type (element_type), allocatable, intent (inout) :: elements(:)
    
    integer :: iel
    double precision :: volmat(3,3)
    double precision, parameter ::  one_sixth = 0.166666666666667
    double precision, parameter ::  alpha = 0.58541020
    double precision, parameter ::  beta = 0.13819660

    ! vtot = 0.0
    do iel = 1,nel

!       ! cubic
!       elements(iel)%xintpt(:,1) = &
!         elements(iel)%xnode(:,1)*0.25 + elements(iel)%xnode(:,2)*0.25 + &
!         elements(iel)%xnode(:,3)*0.25 + elements(iel)%xnode(:,4)*0.25
! 
!       elements(iel)%xintpt(:,2) = &
!         elements(iel)%xnode(:,1)*0.5 + elements(iel)%xnode(:,2)*one_sixth + &
!         elements(iel)%xnode(:,3)*one_sixth + elements(iel)%xnode(:,4)*one_sixth
! 
!       elements(iel)%xintpt(:,3) = &
!         elements(iel)%xnode(:,1)*one_sixth + elements(iel)%xnode(:,2)*0.5 + &
!         elements(iel)%xnode(:,3)*one_sixth + elements(iel)%xnode(:,4)*one_sixth
! 
!       elements(iel)%xintpt(:,4) = &
!         elements(iel)%xnode(:,1)*one_sixth + elements(iel)%xnode(:,2)*one_sixth + &
!         elements(iel)%xnode(:,3)*0.5 + elements(iel)%xnode(:,4)*one_sixth
! 
!       elements(iel)%xintpt(:,5) = &
!         elements(iel)%xnode(:,1)*one_sixth + elements(iel)%xnode(:,2)*one_sixth + &
!         elements(iel)%xnode(:,3)*one_sixth + elements(iel)%xnode(:,4)*0.5
! 
!       elements(iel)%w = [-0.8, 0.45, 0.45, 0.45, 0.45]

      ! quadratic
      elements(iel)%xintpt(:,1) = &
        elements(iel)%xnode(:,1)*alpha + elements(iel)%xnode(:,2)*beta + &
        elements(iel)%xnode(:,3)*beta + elements(iel)%xnode(:,4)*beta

      elements(iel)%xintpt(:,2) = &
        elements(iel)%xnode(:,1)*beta + elements(iel)%xnode(:,2)*alpha + &
        elements(iel)%xnode(:,3)*beta + elements(iel)%xnode(:,4)*beta

      elements(iel)%xintpt(:,3) = &
      elements(iel)%xnode(:,1)*beta + elements(iel)%xnode(:,2)*beta + &
      elements(iel)%xnode(:,3)*alpha + elements(iel)%xnode(:,4)*beta

      elements(iel)%xintpt(:,4) = &
      elements(iel)%xnode(:,1)*beta + elements(iel)%xnode(:,2)*beta + &
      elements(iel)%xnode(:,3)*beta + elements(iel)%xnode(:,4)*alpha

      elements(iel)%w = [0.25, 0.25, 0.25, 0.25]


      ! volmat(:,1) = elements(iel)%xnode(:,1) - elements(iel)%xnode(:,4)
      ! volmat(:,2) = elements(iel)%xnode(:,2) - elements(iel)%xnode(:,4)
      ! volmat(:,3) = elements(iel)%xnode(:,3) - elements(iel)%xnode(:,4)

      ! elements(iel)%v = abs(determinant33(volmat))*one_sixth
      ! elements(iel)%v = abs(volume_tetra(elements(iel)%xnode))
      ! write(*,*) iel! , elements(iel)%v

      ! verify
      ! elements(iel)%v = abs(dot_product(elements(iel)%xnode(:,1) - elements(iel)%xnode(:,4), &
      !   cross_prod(elements(iel)%xnode(:,2) - elements(iel)%xnode(:,4),elements(iel)%xnode(:,3) - elements(iel)%xnode(:,4))))*one_sixth

      ! vtot = vtot + elements(iel)%v

    enddo
    ! write(*,*) 'vtot = ', vtot

  end

  subroutine write_txfft(elements, grains)
    use global
    use types
    use tensor_functions
  
    type (element_type), allocatable, intent (in) :: elements(:)
    type (grain_type), allocatable, intent (in) :: grains(:)

    integer :: iel, ip1, ip2, ip3, gr_id

    open (unit = 33, file = 'txfft', status = 'unknown')
    do ip3 = 1,npts3
      do ip2 = 1,npts2
        do ip1 = 1,npts1
          
          iel = el_id(ip1,ip2,ip3)
          gr_id = elements(iel)%grain_id
          write(33,*) grains(gr_id)%eul, ip1, ip2, ip3, gr_id, grains(gr_id)%phase

        enddo
      enddo
    enddo
    close(33)

  end
end module IO_functions
