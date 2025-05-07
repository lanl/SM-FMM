program main
  !$ use omp_lib
  use types
  use IO_functions
  use tensor_functions
  use global
  use fourier_functions_r2c
  use various_functions
  implicit none

  type (edge_type), allocatable :: edges(:)
  type (face_type), allocatable :: faces(:)
  type (element_type), allocatable :: elements(:)
  type (box_type), allocatable :: boxes(:)
  type (grain_type), allocatable :: grains(:)
  type (intg_pt_type) :: intg_pt(nintg_pt)
  type (intg_pt_type) :: intg_pt_quad(nintg_pt_quad)
  type (fftw_type) :: ffts

  integer :: iel, i, iel1, iel2, ib, ib1, ib2, iel_print, ifc, iedg, iedg_print, j, ifc_print, ibfar, ic, itertot
  character(len = 500) :: input_file, vtk_file
  double precision :: G(3,3,3,3), Gtmp(3,3,3,3), strainacc(3,3), x(3), t1, t2, tolerance_tmp
  double precision :: aux33(3,3), dum, dum1, avg_disp_b, k_penalty_init, Eapp_tot(3,3)
  double precision :: uface_norm_avg, uedge_norm_avg, dtconv, t1conv, t2conv, f_0, errold
  double precision :: t_outgoing, t_incoming, t_near, err_bc_old

  logical :: skip_next = .false.
  character(len=32) :: arg
  character(len=5) :: str, str_incr

  ! call test_nD
  ! call test_nD_many
  ! call allocate_and_plan(ffts,3,3,3,nthreads)
  ! ! call test_fft_arr(ffts)
  ! call test_conv_arr(ffts)
  ! stop

  ! call initialize_fftw(10,10,10,3,3,3,.false.,1)
  ! call test_conv_arr_cmplx
  ! stop

  ! ! test mandel
  ! call test_Mandel

  ! ! test loop
  ! x = 0.0
  ! write(*,*) 128**3
  ! !$ t1 = omp_get_wtime()
  ! do i = 1,128**3
  !   if (mod(i,10000) == 0) write(*,*)float(i)/float(128**3)
  !   do j = 1,128**3
  !     dum = norm2(x)
  !   enddo
  ! enddo
  ! !$ t2 = omp_get_wtime()
  ! !$ write(*,*) t2 - t1
  ! stop

  ! Read optional command line arguments
  do i = 1, command_argument_count()
    if (skip_next) then
      skip_next = .false.
      cycle
    endif
    call get_command_argument(i, arg)
  
    select case (arg)

      case ('--nthreads')
        skip_next = .true.
        call get_command_argument(i+1, arg)
        read(arg, *) nthreads

      case default
        write(*,'(3A)') 'Unknown command-line option: "', trim(arg), '" use -h or --help for argument list'
        stop

    end select
  end do

  ! initialize mapping arrays
  call init_map_outgoing_exp_arr_to_1d

  ! initialize integration points
  call initialize_integration_points(intg_pt)
  call initialize_integration_points_quad(intg_pt_quad)

  ! ! test
  ! call test_intgpt(intg_pt)

  ! load file
  input_file = 'fmm.in'
  write(*,*)'load_input'
  !$ t1 = omp_get_wtime()
  call load_input(input_file, elements, grains, boxes, ffts)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  ! ! test
  ! call test_plastic_rate(elements)

  ! regularization coeficients
  f_0 = 1.0/dx_min*2.0
  coef_reg = calc_coef(dx_min, f_0)
  ! write(*,*) coef_reg
  ! write(*,*) f_0
  ! x = 0.0
  ! x(1) = dx_min + 1.0e-4
  ! write(*,*) x
  ! write(*,*) Green_function(x, lambda, mu)
  ! x(1) = dx_min - 1.0e-4
  ! write(*,*) x
  ! write(*,*) Green_function(x, lambda, mu)
  ! x(1) = dx_min/10.0
  ! write(*,*) x
  ! write(*,*) Green_function(x, lambda, mu)
  ! x(1) = 2.0e-2
  ! write(*,*) x
  ! write(*,*) Green_function(x, lambda, mu)
  ! x(1) = 0.0
  ! write(*,*) x
  ! write(*,*) Green_function(x, lambda, mu)
  ! stop

  ! test
  Sapp = 0.0
  Sapp(3,3) = 1000.0

  ! call test_Green_function_integration(intg_pt)
  call test_Green_function_integration_quad(intg_pt, intg_pt_quad)
  call test_Green_function_integration_nearly_sing(intg_pt_quad)

  write(*,*)'number of elements:', nel
  write(*,*)'suggested number of boxes (total and per dim):', nint(float(nel)/10.0), nint((float(nel)/10.0)**0.33333333)

  ! write(*,*)'write_vtk'
  ! !$ t1 = omp_get_wtime()
  ! vtk_file = 'init.vtk'
  ! call write_vtk(vtk_file, elements)
  ! !$ t2 = omp_get_wtime()
  ! !$ write(*,*) t2 - t1

  write(*,*)'points_to_elements'
  !$ t1 = omp_get_wtime()
  call points_to_elements(elements)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'element_aspect_ratio'
  !$ t1 = omp_get_wtime()
  call element_aspect_ratio(elements, faces)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  if (periodic) then
    write(*,*)'create_faces_per'
    !$ t1 = omp_get_wtime()
    call create_faces_per(elements, faces)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  else
    write(*,*)'create_faces'
    !$ t1 = omp_get_wtime()
    call create_faces(elements, faces)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  endif

  if (periodic) then
    write(*,*)'create_edges_per'
    !$ t1 = omp_get_wtime()
    call create_edges_per(elements, faces, edges)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  else
    write(*,*)'create_edges'
    !$ t1 = omp_get_wtime()
    call create_edges(elements, faces, edges)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  endif

  write(*,*)'nel, nfc, nedg:', nel, nfc, nedg

  write(*,*)'integrate_Green_function_edges'
  !$ t1 = omp_get_wtime()
  call integrate_Green_function_edges(elements, faces, edges, intg_pt, intg_pt_quad)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'write_txfft'
  !$ t1 = omp_get_wtime()
  call write_txfft(elements, grains)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'points_to_boxes'
  !$ t1 = omp_get_wtime()
  call points_to_boxes(boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'elements_to_boxes'
  !$ t1 = omp_get_wtime()
  call elements_to_boxes(elements,boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'boxes_near_far'
  !$ t1 = omp_get_wtime()
  call boxes_near_far(boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'edges_faces_to_boxes'
  !$ t1 = omp_get_wtime()
  call edges_faces_to_boxes(edges,faces,elements,boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  ib = 1
  write(*,*) boxes(ib)%nfar
  ! ! write(*,*) 'near_edge_unique_ids', boxes(ib)%near_edge_unique_ids
  ! write(*,*) size(boxes(ib)%near_edge_unique_ids), nedg
  ! do i = 1, size(boxes(ib)%near_edge_unique_ids)
  !   do j = 1,i-1
  !     if (boxes(ib)%near_edge_unique_ids(j) == boxes(ib)%near_edge_unique_ids(i)) stop 'issue'
  !   enddo
  ! enddo
  ! stop

  write(*,*)'element_moments'
  !$ t1 = omp_get_wtime()
  call element_moments(elements,boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'element_moments'
  !$ t1 = omp_get_wtime()
  call edge_moments(edges,boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  ! write(*,*)'write_vtk_all_cells'
  ! !$ t1 = omp_get_wtime()
  ! vtk_file = 'init_all.vtk'
  ! call write_vtk_all_cells(vtk_file, elements, faces, edges)
  ! !$ t2 = omp_get_wtime()
  ! !$ write(*,*) t2 - t1

  if (solution_method == 2 .and. (.not.plasticity)) then
    !$ t1 = omp_get_wtime()
    call calc_I_C0S_inv(elements)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  endif

  if (ifull == 0) then

    ! write(*,*)'element_close_inter_ids'
    ! !$ t1 = omp_get_wtime()
    ! call element_close_inter_ids(elements,boxes)
    ! !$ t2 = omp_get_wtime()
    ! !$ write(*,*) t2 - t1

    if (FFT) then

      if (periodic) then

        ! write(*,*)'box_far_inter_G_FFT_per'
        ! !$ t1 = omp_get_wtime()
        ! call box_far_inter_G_FFT_per(boxes)
        ! !$ t2 = omp_get_wtime()
        ! !$ write(*,*) t2 - t1
        stop

      else 

        ! if (order == 4) stop 'order 4 not implemented in FFT'

        write(*,*)'box_far_inter_G_FFT'
        !$ t1 = omp_get_wtime()
        call box_far_inter_G_FFT(boxes, ffts)
        !$ t2 = omp_get_wtime()
        !$ write(*,*) t2 - t1

!        ! test
!        write(*,*)'box_far_inter_G'
!        !$ t1 = omp_get_wtime()
!        call box_far_inter_G(boxes)
!        !$ t2 = omp_get_wtime()
!        !$ write(*,*) t2 - t1
!
!        ib = 1
!        do i = 1,boxes(ib)%nfar
!          ibfar = boxes(ib)%far_box_id(i)
!          write(*,*) boxes(ibfar)%ib_grid
!          write(*,*) norm2(boxes(ib)%G(:,:,i) - &
!           G_box(:,:,boxes(ibfar)%ib_grid(1),boxes(ibfar)%ib_grid(2),boxes(ibfar)%ib_grid(3)))
!          write(*,*) norm2(boxes(ib)%dG_dX(:,:,:,i) + &
!           dG_dX_box(:,:,:,boxes(ibfar)%ib_grid(1),boxes(ibfar)%ib_grid(2),boxes(ibfar)%ib_grid(3)))
!          write(*,*) norm2(boxes(ib)%d2G_dX2(:,:,:,:,i) - &
!           d2G_dX2_box(:,:,:,:,boxes(ibfar)%ib_grid(1),boxes(ibfar)%ib_grid(2),boxes(ibfar)%ib_grid(3)))
!          write(*,*) norm2(boxes(ib)%d3G_dX3(:,:,:,:,:,i) + &
!           d3G_dX3_box(:,:,:,:,:,boxes(ibfar)%ib_grid(1),boxes(ibfar)%ib_grid(2),boxes(ibfar)%ib_grid(3)))
!          write(*,*) boxes(ib)%d3G_dX3(1,1,1,1,1,i), &
!           d3G_dX3_box(1,1,1,1,1,boxes(ibfar)%ib_grid(1),boxes(ibfar)%ib_grid(2),boxes(ibfar)%ib_grid(3))
!        enddo
!        stop

      endif

    else

     ! if (periodic) stop 'only FFT for periodic'

     write(*,*)'box_far_inter_G'
     !$ t1 = omp_get_wtime()
     call box_far_inter_G(boxes)
     !$ t2 = omp_get_wtime()
     !$ write(*,*) t2 - t1

    endif

  endif

  if (periodic .and. ifull == 0 .and. precalc_G) then
    write(*,*)'calc_G'
    !$ t1 = omp_get_wtime()
    call calc_G(edges, boxes)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1
  endif

  write(*,*)'integrate_Green_function_edges_close'
  !$ t1 = omp_get_wtime()
  ! call integrate_Green_function_edges_close(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  call integrate_Green_function_edges_close_test(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  ! call integrate_Green_function_edges_dist_el(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  ! call test_geometry(faces, elements)
  if (iJacobi == 1) then
    write(*,*)'self_interaction_Gamma'
    !$ t1 = omp_get_wtime()
    call self_interaction_Gamma(elements, faces, edges, intg_pt)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    if ((.not.plasticity) .and. solution_method == 1) then
      write(*,*)'calc_I_GdC_inv'
      !$ t1 = omp_get_wtime()
      call calc_I_GdC_inv(elements)
      !$ t2 = omp_get_wtime()
      !$ write(*,*) t2 - t1
    endif
  endif

  ! open output file
  open (unit = 77, file = 'averages.out', status = 'unknown')
  write(77,'(A)')'incr, stress_avg, strain_avg - Mandel notation'

  ! total load
  call total_load(boxes, faces, edges)

  ! vtk_file = 'grid.vtk'
  ! call write_grid_vtk(vtk_file, boxes)
  ! ! stop

  ! ! force eigenstrain
  ! do iel = 1,nel
  !   if (elements(iel)%grain_id == 2) elements(iel)%eigen = id3*0.1
  ! enddo

  ! G coeficients
  G_coef1 = (lambda + 3.0*mu)/(8.0*pi*mu*(lambda + 2.0*mu))
  G_coef2 = (lambda + mu)/(8.0*pi*mu*(lambda + 2.0*mu))

  ! eigenstrain in elements
  ! open(111, file = 'eigenstrain', form='unformatted',access='sequential', status = 'unknown')
  do iel = 1,nel
    elements(iel)%eigen = mandel_t2_to_t1(grains(elements(iel)%grain_id)%eigenstrain)
    ! read(111) elements(iel)%eigen 
  enddo
  ! close(111)

  ! iterate
  iel_print = 1503
  iedg_print = 704
  ifc_print = 6665

  k_penalty_init = k_penalty
  ! Eapp_tot = Eapp
  Eapp = 0.0
  do incr = 1,nincr

    Eapp = Edotapp*float(incr)*time_inc

    write(*,*) 'Increment:', incr
    write(*,*) 'Applied strain:'
    write(*,'(3E16.8)') (Eapp(i,:),i=1,3)

    if (incr == 1) then
      call initial_guess(elements)
      call calc_F_faces(elements, faces)
      call calc_F_faces_box(elements, faces, boxes)
      call calc_f_edges(faces, edges)
      call calc_f_edges_box(faces, edges, boxes)
      call calc_near_f_box(boxes, edges)
    endif

    ! write(*,*) 'elements(iel_print)%Gamma_self = ', elements(iel_print)%Gamma_self
    ! write(*,*) 'elements(iel_print)%stress = ', elements(iel_print)%stress
    ! write(*,*) 'elements(iel_print)%strain = ', elements(iel_print)%strain
    ! write(*,*) 'elements(iel_print)%eigen = ', elements(iel_print)%eigen
    ! write(*,*) 'elements(iel_print)%grain_id = ', elements(iel_print)%grain_id
    ! write(*,*) 'elements(iel_print)%C = ', elements(iel_print)%C
    ! write(*,*) 'elements(iel_print)%polar = ', elements(iel_print)%polar
    ! do i = 1,4
    !   ifc = elements(iel_print)%face(i)
    !   write(*,*) 'i = ', i
    !   write(*,*) 'ifc = ', ifc
    !   write(*,*) 'faces(ifc)%u = ', faces(ifc)%u
    !   write(*,*) 'faces(ifc)%area = ', faces(ifc)%area
    !   write(*,*) 'faces(ifc)%normal = ', faces(ifc)%normal
    ! enddo
    ! write(*,*) 'edges(iedg_print)%u = ', edges(iedg_print)%u
    ! ! read(*,*)

    !$ dtconv = 0.0
    !$ t_outgoing = 0.0
    !$ t_incoming = 0.0
    !$ t_near = 0.0

    err_bc = 1.0
    err_bc_old = 1.0
    iter_out = 0 
    tolerance_tmp = tolerance
    k_penalty = k_penalty_init
    itertot = 0
    do while ((err_bc.gt.tolerance_out.and.iter_out.lt.itout_mx)) ! .or. (iter_out.lt.5))
  
      iter_out = iter_out + 1
  
      if (iter_out.gt.1) then
        write(*,*) 'err_bc/err_bc_old', err_bc/err_bc_old
      !   if (err_bc/err_bc_old >= 0.9 .and. err_bc/err_bc_old < 0.99) k_penalty = k_penalty*k_penalty_incr
      endif
      k_penalty = k_penalty*k_penalty_incr
      write(*,*) 'k_penalty', k_penalty
    
      iter = 0
      err = 1.0
      errold = 1.0
  
      do while ((err.gt.tolerance_tmp.and.iter.lt.itmx).or.iter<itmin)
       iter = iter + 1
       itertot = itertot + 1
    
       ! strain as convolution
       do iel1 = 1,nel
         elements(iel1)%strainold = elements(iel1)%strain
       enddo
    
       !$ t1conv = omp_get_wtime()
       if (ifull == 0) then
         write(*,*) 'outgoing_expansion'
         !$ t1 = omp_get_wtime()
         call outgoing_expansion_f(edges,faces,elements,boxes)
         !$ t2 = omp_get_wtime()
         !$ t_outgoing = t_outgoing + t2 - t1
         !$ write(*,*) t2 - t1
    
         write(*,*) 'incoming_expansion'
         !$ t1 = omp_get_wtime()
         if (FFT) then
           ! stop 'zeroth order term not implemented'
           if (periodic) then
             ! call incoming_expansion_f_FFT_per(boxes)
            stop
           else
             ! call incoming_expansion_f_FFT(boxes)
             call incoming_expansion_f_FFT_1d(boxes, ffts)
           endif
         else
           call incoming_expansion_f(boxes)
         endif
         !$ t2 = omp_get_wtime()
         !$ t_incoming = t_incoming + t2 - t1
         !$ write(*,*) t2 - t1
       endif
    
       write(*,*) 'convolution_edges'
       !$ t1 = omp_get_wtime()
       if (ifull == 1) then
         ! call convolution_edges(faces, edges)
         call convolution_edges_box(elements, faces, edges, boxes)
         call calc_u_faces(edges, faces)
         call calc_disgrad_elements(faces, elements, boxes)
       else
         ! call convolution_all(boxes, elements, faces, edges)
         call convolution_all_box(boxes, elements, faces, edges)
       endif
       !$ t2 = omp_get_wtime()
       !$ t_near = t_near + t2 - t1
       !$ t2conv = omp_get_wtime()
       !$ write(*,*) t2 - t1
       !$ dtconv = dtconv + (t2conv - t1conv)
    
       ! write(*,*) 'faces(ifc_print)%F = ', faces(ifc_print)%F
       ! write(*,*) 'faces(ifc_print)%u = ', faces(ifc_print)%u
       ! write(*,*) 'faces(ifc_print)%normal = ', faces(ifc_print)%normal
       ! write(*,*) 'faces(ifc_print)%el = ', faces(ifc_print)%el
       ! write(*,*) 'elements(faces(ifc_print)%el(1))%polar = ', elements(faces(ifc_print)%el(1))%polar
       ! write(*,*) 'elements(faces(ifc_print)%el(2))%polar = ', elements(faces(ifc_print)%el(2))%polar
       ! do i = 1,3
       !   iedg = faces(ifc_print)%edge(i)
       !   write(*,*) 'i = ', i
       !   write(*,*) 'iedg = ', iedg
       !   write(*,*) 'edges(iedg)%u = ', edges(iedg)%u
       !   write(*,*) 'edges(iedg)%f = ', edges(iedg)%f
       ! enddo
    
       ! write(str,'(I5.5)')iter
       ! vtk_file = 'output_all_'//str//'.vtk'
       ! call write_vtk_all_cells(vtk_file, elements, faces, edges)
       ! ! stop
    
       if (plasticity .and. nit_non_loc > 0) then
          call non_local_hyd_correction(elements, faces, grains, boxes)
       endif
    
       errold = err
       err = 0.0
       dum = 0.0
       do iel1 = 1,nel
         if (elements(iel1)%grain_id .ne. id_gas) then
           err = err + sum((elements(iel1)%strainold - elements(iel1)%strain)**2)
         endif
         dum = dum + norm2(elements(iel1)%strain)
         vtot = vtot + elements(iel1)%v
       enddo
       dum = dum/float(nel)
       err = sqrt(err/float(nel))/dum ! norm2(Eapp)
    
    !   if (err > errold .and. iter > 2 .and. step > 0.01) then
    !    write(*,*) 'cut step'
    !    step = step*0.2
    !    if (step < 0.01) step = 0.01
    !   endif
    !
       do iel1 = 1,nel
         elements(iel1)%strain = elements(iel1)%strainold*(1.0-step) + elements(iel1)%strain*step
       enddo
  
       write(*,*) 'calculate stress and polarization'
       !$ t1 = omp_get_wtime()
       if (.not.plasticity) then
         call calc_stress_polar_elastic(elements)
       else
        if (solution_method == 1) then
          stop 'not implemented basic for plastic'
        else
          call calc_stress_polar_plastic_AL(elements)
        endif
       endif
       !$ t2 = omp_get_wtime()
       !$ t2conv = omp_get_wtime()
       !$ write(*,*) t2 - t1

       write(*,*) 'iter, err = ', iter, err, err_polar
  
       write(*,*) 'calculate forces'
       !$ t1 = omp_get_wtime()
       call calc_F_faces(elements, faces)
       call calc_F_faces_box(elements, faces, boxes)
       call calc_f_edges(faces, edges)
       call calc_f_edges_box(faces, edges, boxes)
       call calc_near_f_box(boxes, edges)
       !$ t2 = omp_get_wtime()
       !$ t2conv = omp_get_wtime()
       !$ write(*,*) t2 - t1
       
       ! write(*,*) 'elements(iel_print)%stress = ', elements(iel_print)%stress
       ! write(*,*) 'elements(iel_print)%strain = ', elements(iel_print)%strain
       ! write(*,*) 'elements(iel_print)%grain_id = ', elements(iel_print)%grain_id
       ! write(*,*) 'elements(iel_print)%polar = ', elements(iel_print)%polar
       ! do i = 1,4
       !   ifc = elements(iel_print)%face(i)
       !   write(*,*) 'i = ', i
       !   write(*,*) 'ifc = ', ifc
       !   write(*,*) 'faces(ifc)%u = ', faces(ifc)%u
       !   write(*,*) 'faces(ifc)%area = ', faces(ifc)%area
       !   write(*,*) 'faces(ifc)%normal = ', faces(ifc)%normal
       ! enddo
       ! write(*,*) 'edges(iedg_print)%u = ', edges(iedg_print)%u
       ! ! read(*,*)
    
       ! dum = 0.0
       ! ic = 0
       ! do ifc = 1,nfc
       !   if (faces(ifc)%boundary) then
       !    iel1 = faces(ifc)%el(1)
       !    dum = dum + norm2(matmul(elements(iel1)%stress - stressavg, faces(ifc)%normal))
       !    ic = ic + 1
       !   endif
       ! enddo
       ! write(*,*) 'traction err',dum/float(ic), dum/norm2(Sapp)/float(ic)
    
      enddo
    
      err_bc_old = err_bc
      err_bc = 0.0
      avg_disp_b = 0.0
      ic = 0
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(ifc,iel1) &
      !$OMP SHARED(nfc,faces,k_penalty,Eapp,elements,id_gas) &
      !$OMP REDUCTION(+:err_bc,avg_disp_b,ic)
      !$OMP DO 
      do ifc = 1,nfc
        iel1 = faces(ifc)%el(1)
        if (faces(ifc)%boundary .and. elements(iel1)%grain_id .ne. id_gas) then
         ! faces(ifc)%Fbc = faces(ifc)%Fbc + matmul(elements(iel1)%stress - stressavg, faces(ifc)%normal)*faces(ifc)%area*k_penalty
         faces(ifc)%Fbc = faces(ifc)%Fbc - faces(ifc)%u*k_penalty
         ! faces(ifc)%Fbc = faces(ifc)%F
         err_bc = err_bc + norm2(faces(ifc)%u)
         avg_disp_b = avg_disp_b + norm2(matmul(Eapp,faces(ifc)%x))
         ic = ic + 1
        endif
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      err_bc = err_bc/float(ic)
      avg_disp_b = avg_disp_b/float(ic)
      err_bc = err_bc/avg_disp_b
      tolerance_tmp = err_bc*0.5
      if (tolerance_tmp < tolerance) tolerance_tmp = tolerance
      write(*,*) 'iterout, err disp surf', iter_out, err_bc
  
    enddo

    call calc_write_averages(elements)

    call update_state(elements)

    write(str_incr,'(I5.5)')incr
    vtk_file = 'output_incr_'//str_incr//'.vtk'
    call write_vtk(vtk_file, elements)

  enddo

  !$ write(*,*) 'avg_conv_time', dtconv/float(itertot)
  !$ write(*,*) 'avg_outgoing', t_outgoing/float(itertot)
  !$ write(*,*) 'avg_incoming', t_incoming/float(itertot)
  !$ write(*,*) 'avg_near', t_near/float(itertot)

!   stressavg = 0.0
!   strainavg = 0.0
!   do iel1 = 1,nel
!     stressavg = stressavg + elements(iel1)%stress(:,:)*elements(iel1)%v
!     strainavg = strainavg + elements(iel1)%strain(:,:)*elements(iel1)%v
!   enddo
!   stressavg = stressavg/vtot
!   strainavg = strainavg/vtot
! 
!   write(*,*) 'stressavg = ', stressavg
!   write(*,*) 'strainavg = ', strainavg
! 
!   uface_norm_avg = 0.0
!   do ifc = 1,nfc
!     uface_norm_avg = uface_norm_avg + norm2(faces(ifc)%u)
!   enddo
!   uface_norm_avg = uface_norm_avg/float(nfc)
!   write(*,*) 'uface_norm_avg = ', uface_norm_avg
! 
!   uedge_norm_avg = 0.0
!   do iedg = 1,nedg
!     uedge_norm_avg = uedge_norm_avg + norm2(edges(iedg)%u)
!   enddo
!   uedge_norm_avg = uedge_norm_avg/float(nedg)
!   write(*,*) 'uedge_norm_avg = ', uedge_norm_avg

  ! vtk_file = 'output.vtk'
  ! call write_vtk(vtk_file, elements)
  ! vtk_file = 'output_all.vtk'
  ! call write_vtk_all_cells(vtk_file, elements, faces, edges)

!   ! analytical Eshelby
!   call analytical_sol_thermal_exp(elements,grains)
! 
!   vtk_file = 'output_analytical.vtk'
!   call write_vtk(vtk_file, elements)

  write(*,*) 'elements array size:', sizeof(elements)
  write(*,*) 'faces array size:', sizeof(faces)
  write(*,*) 'edges array size:', sizeof(edges)
  write(*,*) 'boxes array size:', sizeof(boxes)

  ! open(111, file = 'eigenstrain', form='unformatted',access='sequential', status = 'unknown')
  ! do iel = 1,nel
  !   write(111) elements(iel)%strain_pl
  ! enddo
  ! close(111)

  close(77)

end
