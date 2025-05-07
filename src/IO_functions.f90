module IO_functions
  implicit none

contains


  subroutine load_input(input_file, elements, grains, boxes, ffts)
    use global
    use fourier_functions_r2c
    use types
    use various_functions

    character(len = 500), intent (in) :: input_file
    type (element_type), allocatable, intent (inout) :: elements(:)
    type (grain_type), allocatable, intent (inout) :: grains(:)
    type (box_type), allocatable, intent (inout) :: boxes(:)
    type (fftw_type), intent (inout) :: ffts

    character(len = 500) :: mesh_file, tex_file, elas_file, eigenstrain_file, slip
    integer :: iFFT, iperiodic, ieigenstrain_flag, iplast

    open (unit = 11, file = input_file, status = 'old')

    read (11,*) 
    read (11,'(A)') mesh_file
    read (11,'(A)') tex_file
    read (11,'(A)') elas_file
    read (11,*) npts1, npts2, npts3
    read (11,*) nbox1, nbox2, nbox3
    read (11,*) id_gas
    read (11,*) 
    read (11,*) Edotapp(1,1), Edotapp(1,2), Edotapp(1,3)
    read (11,*) Edotapp(2,1), Edotapp(2,2), Edotapp(2,3)
    read (11,*) Edotapp(3,1), Edotapp(3,2), Edotapp(3,3)
    read (11,*) time_inc, nincr
    read (11,*) 
    read (11,*) ieigenstrain_flag
    read (11,'(A)') eigenstrain_file
    read (11,*) 
    read (11,*) ifull
    read (11,*) iload
    read (11,*) iJacobi
    read (11,*) xc0
    read (11,*) step
    read (11,*) order
    read (11,*) iFFT
    read (11,*) iperiodic
    read (11,*) dx_min
    read (11,*) tolerance, itmx, itmin
    read (11,*) tolerance_out, itout_mx, k_penalty, k_penalty_incr
    read (11,*) dx_th
    read (11,*) iplast, slip, tauc, rate_exp, nit_non_loc

    close (11)

    if (iFFT == 0) then
      FFT = .false.
    elseif (iFFT == 1) then
      FFT = .true.
    endif

    if (iperiodic == 0) then
      periodic = .false.
    elseif (iperiodic == 1) then
      periodic = .true.
    endif

    if (iplast == 0) then
      plasticity = .false.
    elseif (iplast == 1) then
      plasticity = .true.
    endif

    call load_vtk(mesh_file, elements)
    call data_crystal_elast(elas_file)
    call read_and_init_grain(tex_file,grains)
    if (ieigenstrain_flag == 1) call data_grain_eigen(eigenstrain_file,grains)

    call renumber_el_grain_ids(elements)

    call allocate_globals
    !$ call setup_openmp()
    call allocate_and_plan(ffts,nbox1,nbox2,nbox3,nthreads)

    ! test box for same nbox as npts
    ! write(*,*) sqrt(sum((Goperreal-Goperreal_box)**2)) 

    if (periodic) call enforce_domain(elements)
    call init_integration_pt(elements)

    nbox = nbox1*nbox2*nbox3
    allocate(boxes(nbox))

    if (plasticity .and. slip(1:3) == 'fcc') then
      call init_slip_fcc
    endif

    call grains_to_elements(grains,elements)

    if (plasticity) then
      call allocate_slip_arr(elements)
      call calc_Schmid(elements)
      call assign_tauc(elements)
    endif

    ! stop

  end

  subroutine load_vtk(mesh_file, elements)
    use global
    use types
    use various_functions

    character(len = 500), intent (in) :: mesh_file
    type (element_type), allocatable, intent (inout) :: elements(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel
    integer :: nh_lines = 4
    double precision, allocatable :: x_node(:,:)
    double precision :: dum
    logical, allocatable :: tetra(:)
    logical :: gmsh = .true.

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

    if (gmsh) read(33,'(A)')line

    read(33,*) str1, neltot, i
    if (trim(str1) .ne. 'CELLS') stop 'trim(str1) .ne. CELLS'

    allocate(ind(20,neltot))
    allocate(tetra(neltot))
    ind = 0
    do i = 1,neltot

      read(33,*) icell, (ind(j,i), j = 1,icell)

    enddo

    if (gmsh) read(33,'(A)')line

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
        elements(i)%v = abs(volume_tetra(elements(i)%xnode))
        elements(i)%x = &
        elements(i)%xnode(:,1)*0.25 + elements(i)%xnode(:,2)*0.25 + &
        elements(i)%xnode(:,3)*0.25 + elements(i)%xnode(:,4)*0.25
      endif
      
    enddo

    if (gmsh) read(33,'(A)')line
    read(33,'(A)') line
    read(33,'(A)') line
    read(33,'(A)') line

    i = 0
    do iel = 1,neltot

      if (tetra(iel)) then
        i = i + 1
        if (gmsh) then
        read(33,*) elements(i)%grain_id
        else
          read(33,*) dum
          elements(i)%grain_id = nint(dum)
        endif
      else
        read(33,*)
      endif
      
    enddo

    deallocate(x_node)
    deallocate(ind)
    deallocate(tetra)

    close (33)

  end

  subroutine write_vtk(vtk_file, elements)
    use tensor_functions
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

    number_fields = 11
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

    write(33,'(A,I3,I10,A)')'volume',1,nel,' double'
    do iel = 1, nel
      write(33,'(1E18.8E3)')elements(iel)%v
    enddo

    write(33,'(A,I3,I10,A)')'stress',6,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%stresst(1),elements(iel)%stresst(2), &
      elements(iel)%stresst(3),elements(iel)%stresst(6)/sqrt2,elements(iel)%stresst(4)/sqrt2,&
      elements(iel)%stresst(5)/sqrt2
    enddo

    write(33,'(A,I3,I10,A)')'stressvm',1,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')stress_vm(elements(iel)%stresst)
    enddo

    write(33,'(A,I3,I10,A)')'strain',6,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%straint(1),elements(iel)%straint(2), &
      elements(iel)%straint(3),elements(iel)%straint(6)/sqrt2,elements(iel)%straint(4)/sqrt2,&
      elements(iel)%straint(5)/sqrt2
    enddo

    write(33,'(A,I3,I10,A)')'strain_pl_vm',1,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%strain_pl_vm
    enddo

    write(33,'(A,I3,I10,A)')'edotp',6,nel,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%edotp(1),elements(iel)%edotp(2), &
      elements(iel)%edotp(3),elements(iel)%edotp(6)/sqrt2,elements(iel)%edotp(4)/sqrt2,&
      elements(iel)%edotp(5)/sqrt2
    enddo

    write(33,'(A,I3,I10,A)')'distorted',1,nel,' int'
    do iel = 1, nel
      if (allocated(elements(iel)%iel_interp)) then
        write(33,'(I10)') 1
      else
        write(33,'(I10)') 0
      endif
    enddo

    write(33,'(A,I3,I10,A)')'dpolar',1,nel,' double'
    do iel = 1, nel
      write(33,'(E18.8E3)') norm2(elements(iel)%polar - elements(iel)%polarold)
    enddo

!    write(33,'(A,I3,I10,A)')'dx',3,nel,' double'
!    do iel = 1, nel
!      write(33,'(6E18.8E3)')elements(iel)%dx(1),elements(iel)%dx(2), &
!      elements(iel)%dx(3)
!    enddo
!
!    write(33,'(A,I3,I10,A)')'dxdx-dx*dx',9,nel,' double'
!    do iel = 1, nel
!      ! write(33,'(6E18.8E3)')elements(iel)%dxdx(1,1),elements(iel)%dxdx(2,1), &
!      ! elements(iel)%dxdx(3,1),elements(iel)%dxdx(1,2),elements(iel)%dxdx(2,2),&
!      ! elements(iel)%dxdx(3,2),elements(iel)%dxdx(1,3),elements(iel)%dxdx(2,3),&
!      ! elements(iel)%dxdx(3,3)
!      do i = 1,3
!        do j = 1,3
!          dxdx(i,j) = elements(iel)%dxdx(i,j) - elements(iel)%dx(i)*elements(iel)%dx(j)
!        enddo
!      enddo
!      write(33,'(6E18.8E3)')dxdx(1,1),dxdx(2,1), &
!      dxdx(3,1),dxdx(1,2),dxdx(2,2),&
!      dxdx(3,2),dxdx(1,3),dxdx(2,3),&
!      dxdx(3,3)
!
!    enddo

    close (33)

  end

  subroutine write_grid_vtk(vtk_file, boxes)
    use global
    use types
    use various_functions

    character(len = 500), intent (in) :: vtk_file
    type (box_type), allocatable, intent (in) :: boxes(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel, in, number_fields, ic, ip1, ip2, ip3, ib, ib1, ib2, ib3
    integer :: nh_lines = 4
    double precision :: G(3,3), dx(3)
    double precision, allocatable :: x_node(:,:), G_grid_box(:,:,:,:,:)
    logical, allocatable :: tetra(:)

    open (unit = 33, file = vtk_file, status = 'unknown')
     
    write(33,'(A)')'# vtk DataFile Version 3.0'
    write(33,'(A)')'microstructure data'
    write(33,'(A)')'ASCII'
    write(33,'(A)')'DATASET STRUCTURED_POINTS'
    write(33,'(A,3I10)')'DIMENSIONS', npts1, npts2, npts3
    write(33,'(A,3E18.8E3)')'ORIGIN', 1.0, 1.0, 1.0
    write(33,'(A,3E18.8E3)')'SPACING', 1.0, 1.0, 1.0

    number_fields = 5
    write(33,'(A,I10)')'POINT_DATA',npts1*npts2*npts3
    write(33,'(A,I10)')'FIELD FieldData',number_fields

    write(33,'(A,I3,I10,A)')'el_id',1,npts1*npts2*npts3,' int'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          write(33,'(I10)') el_id(ip1,ip2,ip3)
        enddo
      enddo
    enddo

    write(33,'(A,I3,I10,A)')'point',3,npts1*npts2*npts3,' int'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          write(33,'(3I10)') ip1,ip2,ip3
        enddo
      enddo
    enddo


    write(33,'(A,I3,I10,A)')'Gper',9,npts1*npts2*npts3,' double'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          write(33,'(9E18.8E3)') &
          Gper(1,1,ip1,ip2,ip3),Gper(2,1,ip1,ip2,ip3),Gper(3,1,ip1,ip2,ip3), &
          Gper(1,2,ip1,ip2,ip3),Gper(2,2,ip1,ip2,ip3),Gper(3,2,ip1,ip2,ip3), &
          Gper(1,3,ip1,ip2,ip3),Gper(2,3,ip1,ip2,ip3),Gper(3,3,ip1,ip2,ip3)
        enddo
      enddo
    enddo

    write(33,'(A,I3,I10,A)')'G',9,npts1*npts2*npts3,' double'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          dx(1) = float(ip1 - 1)
          dx(2) = float(ip2 - 1)
          dx(3) = float(ip3 - 1)
          dx = per_dx(dx)
          if (norm2(dx) > 0.0) then
            G = Green_function(dx, lambda, mu)
          else
            G = 0.0
          endif
          write(33,'(9E18.8E3)') &
          G(1,1),G(2,1),G(3,1), &
          G(1,2),G(2,2),G(3,2), &
          G(1,3),G(2,3),G(3,3)
        enddo
      enddo
    enddo

    allocate(G_grid_box(3,3,npts1,npts2,npts3))
    do ib = 1, nbox
      ib1 = boxes(ib)%ib_grid(1)
      ib2 = boxes(ib)%ib_grid(2)
      ib3 = boxes(ib)%ib_grid(3)
      do ip1 = boxes(ib)%ipstart(1),boxes(ib)%ipend(1)
        do ip2 = boxes(ib)%ipstart(2),boxes(ib)%ipend(2)
          do ip3 = boxes(ib)%ipstart(3),boxes(ib)%ipend(3)
            G_grid_box(:,:,ip1,ip2,ip3) = G_box(:,:,ib1,ib2,ib3)
          enddo
        enddo
      enddo
    enddo

    write(33,'(A,I3,I10,A)')'Gper_box',9,npts1*npts2*npts3,' double'
    do ip3 = 1, npts3
      do ip2 = 1, npts2
        do ip1 = 1, npts1
          write(33,'(9E18.8E3)') &
          G_grid_box(1,1,ip1,ip2,ip3),G_grid_box(2,1,ip1,ip2,ip3),G_grid_box(3,1,ip1,ip2,ip3), &
          G_grid_box(1,2,ip1,ip2,ip3),G_grid_box(2,2,ip1,ip2,ip3),G_grid_box(3,2,ip1,ip2,ip3), &
          G_grid_box(1,3,ip1,ip2,ip3),G_grid_box(2,3,ip1,ip2,ip3),G_grid_box(3,3,ip1,ip2,ip3)
        enddo
      enddo
    enddo

    deallocate(G_grid_box)

    close (33)

  end


  subroutine write_vtk_all_cells(vtk_file, elements, faces, edges)
    use global
    use types

    character(len = 500), intent (in) :: vtk_file
    type (element_type), allocatable, intent (inout) :: elements(:)
    type (face_type), allocatable, intent (in) :: faces(:)
    type (edge_type), allocatable, intent (in) :: edges(:)

    character(len = 500) :: line, str1, str2
    integer, allocatable :: ind(:,:)
    integer :: i, neltot, icell, j, iel, in, number_fields, ic, ifc, iedg
    integer :: ind_node_face(3), ind_node_edge(2)
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

    write(33,'(A,I10,I10)')'CELLS',nel + nfc + nedg,nel*5 + nfc*4 + 3*nedg
    do iel = 1, nel
      write(33,'(5I10)')4, elements(iel)%ind_node - 1
    enddo

    do ifc = 1, nfc
      iel = faces(ifc)%el(1)
      i = faces(ifc)%el_face(1)
      ind_node_face(1) = elements(iel)%ind_node(el_face_nodes(1,i))
      ind_node_face(2) = elements(iel)%ind_node(el_face_nodes(2,i))
      ind_node_face(3) = elements(iel)%ind_node(el_face_nodes(3,i))
      write(33,'(4I10)')3, ind_node_face - 1
    enddo

    do iedg = 1, nedg
      ifc = edges(iedg)%face(1)
      iel = faces(ifc)%el(1)
      i = faces(ifc)%el_face(1)
      ind_node_face(1) = elements(iel)%ind_node(el_face_nodes(1,i))
      ind_node_face(2) = elements(iel)%ind_node(el_face_nodes(2,i))
      ind_node_face(3) = elements(iel)%ind_node(el_face_nodes(3,i))
      i = edges(iedg)%face_edge(1)
      ind_node_edge(1) = ind_node_face(face_edge_nodes(1,i))
      ind_node_edge(2) = ind_node_face(face_edge_nodes(2,i))
      write(33,'(3I10)')2, ind_node_edge - 1
    enddo

    write(33,'(A,I10)')'CELL_TYPES',nel + nfc + nedg
    do iel = 1, nel
      write(33,'(I5)') 10 
    enddo
    do ifc = 1, nfc
      write(33,'(I5)') 5 
    enddo
    do iedg = 1, nedg
      write(33,'(I5)') 3 
    enddo

    number_fields = 12
    write(33,'(A,I10)')'CELL_DATA',nel + nfc + nedg
    write(33,'(A,I10)')'FIELD FieldData',number_fields

    write(33,'(A,I3,I10,A)')'cell_type',1,nel + nfc + nedg,' int'
    do iel = 1, nel
      write(33,'(I10)') 10 
    enddo
    do ifc = 1, nfc
      write(33,'(I10)') 5 
    enddo
    do iedg = 1, nedg
      write(33,'(I10)') 3 
    enddo

    write(33,'(A,I3,I10,A)')'ind',1,nel + nfc + nedg,' int'
    do iel = 1, nel
      write(33,'(I10)') iel 
    enddo
    do ifc = 1, nfc
      write(33,'(I10)') ifc
    enddo
    do iedg = 1, nedg
      write(33,'(I10)') iedg
    enddo

    write(33,'(A,I3,I10,A)')'box_id',1,nel + nfc + nedg,' int'
    do iel = 1, nel
      write(33,'(I10)') elements(iel)%box_id
    enddo
    do ifc = 1, nfc
      write(33,'(I10)') 0
    enddo
    do iedg = 1, nedg
      write(33,'(I10)') edges(iedg)%box_id(1)
    enddo

    write(33,'(A,I3,I10,A)')'polar',6,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%polar(1),elements(iel)%polar(2), &
      elements(iel)%polar(3),elements(iel)%polar(6)/sqrt2,elements(iel)%polar(4)/sqrt2,&
      elements(iel)%polar(5)/sqrt2
    enddo
    do ifc = 1, nfc
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo

    write(33,'(A,I3,I10,A)')'disgrad',1,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(E18.8E3)') norm2(elements(iel)%disgrad) 
    enddo
    do ifc = 1, nfc
      write(33,'(E18.8E3)') 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(E18.8E3)') 0.0
    enddo

    write(33,'(A,I3,I10,A)')'stress',6,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%stress(1),elements(iel)%stress(2), &
      elements(iel)%stress(3),elements(iel)%stress(6)/sqrt2,elements(iel)%stress(4)/sqrt2,&
      elements(iel)%stress(5)/sqrt2
    enddo
    do ifc = 1, nfc
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo

    write(33,'(A,I3,I10,A)')'strain',6,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(6E18.8E3)')elements(iel)%strain(1),elements(iel)%strain(2), &
      elements(iel)%strain(3),elements(iel)%strain(6)/sqrt2,elements(iel)%strain(4)/sqrt2,&
      elements(iel)%strain(5)/sqrt2
    enddo
    do ifc = 1, nfc
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(6E18.8E3)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    enddo

    write(33,'(A,I3,I10,A)')'face_F',3,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do ifc = 1, nfc
      write(33,'(3E18.8E3)') faces(ifc)%F
    enddo
    do iedg = 1, nedg
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo

    write(33,'(A,I3,I10,A)')'edge_F',3,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do ifc = 1, nfc
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(3E18.8E3)') edges(iedg)%f
    enddo

    write(33,'(A,I3,I10,A)')'face_u',3,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do ifc = 1, nfc
      write(33,'(3E18.8E3)') faces(ifc)%u
    enddo
    do iedg = 1, nedg
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo

    write(33,'(A,I3,I10,A)')'edge_u',3,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do ifc = 1, nfc
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do iedg = 1, nedg
      write(33,'(3E18.8E3)') edges(iedg)%u
    enddo

    write(33,'(A,I3,I10,A)')'u_tot',3,nel + nfc + nedg,' double'
    do iel = 1, nel
      write(33,'(3E18.8E3)') 0.0, 0.0, 0.0
    enddo
    do ifc = 1, nfc
      write(33,'(3E18.8E3)') faces(ifc)%u + matmul(Eapp, faces(ifc)%x)
    enddo
    do iedg = 1, nedg
      write(33,'(3E18.8E3)') edges(iedg)%u + matmul(Eapp, edges(iedg)%xc)
    enddo

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
  
  subroutine read_and_init_grain(tex_file,grains)
   use global
   use tensor_functions
   use types
   use various_functions
   implicit none
  
   character(len = 500), intent (in) :: tex_file
   type (grain_type), allocatable, intent (inout) :: grains(:)

   double precision :: aa(3,3)
   double precision :: caux3333(3,3,3,3),caux66(6,6)
   double precision :: aux6(6),aux33(3,3), c066b(6,6)
  
   double precision :: dum, om, ph, th
   integer :: i1, i2, ii, j, j1, j2, jgr, jj, jph
   integer :: k, k1, k2, kk, kkk, l1, l2, nph1, i
  
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
     if (kkk == id_gas) grains(kkk)%C = grains(kkk)%C*0.01
     ! write(*,*) 'scaled stiffness of grain', id_gas
     ! if (kkk == 2) grains(kkk)%C = grains(1)%C
     ! write(*,*) 'force same'
     ! read(*,*)
  
     ! call chg_basis(aux6,aux33,grains(kkk)%C66,grains(kkk)%C,4,6)
     ! grains(kkk)%S66 = grains(kkk)%C66
     ! call lu_inverse(grains(kkk)%S66, 6)
     ! call chg_basis(aux6,aux33,grains(kkk)%S66,grains(kkk)%S,3,6)

     grains(kkk)%C66 = mandel_t4_to_t2(grains(kkk)%C)
     grains(kkk)%S66 = grains(kkk)%C66
     call lu_inverse(grains(kkk)%S66, 6)
     grains(kkk)%S = mandel_t2_to_t4(grains(kkk)%S66)

     c066 = c066 + grains(kkk)%C66
  
   enddo
   c066 = c066/float(ng)*xc0
   ! read(*,*)

   ! c066 = force_isotropic(c066)
   c0 = mandel_t2_to_t4(c066)
   call chg_basis(aux6,aux33,c066b,c0,4,6)
   c066b = force_isotropic(c066b)
   call chg_basis(aux6,aux33,c066b,c0,3,6)
   c066 = mandel_t4_to_t2(c0)
   write(*,*) 'forcing isotropic stiffness in b-basis, to be changed'
   ! read(*,*)

   s066 = c066
   call lu_inverse(s066,6)
   
   ! call chg_basis(aux6,aux33,c066,c0,3,6)
   ! call chg_basis(aux6,aux33,s066,s0,3,6)
   ! call get_lambda_mu(c0, lambda, mu)

   c0 = mandel_t2_to_t4(c066)
   s0 = mandel_t2_to_t4(s066)
   call get_lambda_mu(c0, lambda, mu)

   write(*,*) 'c066 = '
   write(*,'(6E16.8)') ((c066(i,j),j=1,6),i=1,6)
   write(*,*) 'lambda = ', lambda
   write(*,*) 'mu = ', mu

   close (33)

   return
  end


  subroutine data_grain_eigen(eigenstrain_file, grains)
    use global
    use tensor_functions
    use types
    use various_functions
    implicit none
   
    character(len = 500), intent (in) :: eigenstrain_file
    type (grain_type), allocatable, intent (inout) :: grains(:)

    double precision :: dum, om, ph, th
    integer :: kkk
   
    open (unit = 33, file = eigenstrain_file, status = 'old')
    read(33,*)
    do kkk = 1,ng
   
      read(33,*) grains(kkk)%eigenstrain(1,1), grains(kkk)%eigenstrain(2,2), grains(kkk)%eigenstrain(3,3), &
        grains(kkk)%eigenstrain(2,3), grains(kkk)%eigenstrain(1,3), grains(kkk)%eigenstrain(1,2)

      grains(kkk)%eigenstrain(2,1) = grains(kkk)%eigenstrain(1,2)
      grains(kkk)%eigenstrain(3,1) = grains(kkk)%eigenstrain(1,3)
      grains(kkk)%eigenstrain(3,2) = grains(kkk)%eigenstrain(2,3)

    enddo
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

  subroutine renumber_el_grain_ids(elements)
    use global
    use types
    use tensor_functions
  
    type (element_type), allocatable, intent (inout) :: elements(:)

    integer :: iel 

    do iel = 1,nel
      elements(iel)%grain_id = mod(elements(iel)%grain_id,ng)
      if(elements(iel)%grain_id == 0) elements(iel)%grain_id = ng
    enddo

  end

  subroutine calc_write_averages(elements)
    use global
    use types
  
    type (element_type), allocatable, intent (inout) :: elements(:)

    double precision :: stress_avg(6), strain_avg(6)
    integer :: iel

    stress_avg = 0.0
    strain_avg = 0.0
    vtot = 0.0
    do iel = 1, nel
      stress_avg = stress_avg + elements(iel)%stress*elements(iel)%v
      strain_avg = strain_avg + elements(iel)%strain*elements(iel)%v
      vtot = vtot + elements(iel)%v
    enddo
    stress_avg = stress_avg/vtot
    strain_avg = strain_avg/vtot

    write(77,'(I5,100E16.8)') incr, stress_avg, strain_avg

  end


end module IO_functions
