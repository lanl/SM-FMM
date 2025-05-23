program main
  !$ use omp_lib
  use types
  use IO_functions
  use tensor_functions
  use global
  use fourier_functions
  use various_functions
  implicit none

  type (element_type), allocatable :: elements(:)
  type (box_type), allocatable :: boxes(:)
  type (grain_type), allocatable :: grains(:)

  integer :: iel, i, iel1, iel2, ib, ib1, ib2
  character(len = 500) :: input_file, vtk_file
  double precision :: G(3,3,3,3), Gtmp(3,3,3,3), strainacc(3,3), x(3), t1, t2

  logical :: skip_next = .false.
  character(len=32) :: arg

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

  ! load file
  input_file = 'fmm.in'
  write(*,*)'load_input'
  !$ t1 = omp_get_wtime()
  call load_input(input_file, elements, grains, boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'number of elements:', nel
  write(*,*)'suggested number of boxes (total and per dim):', nint(float(nel)/10.0), nint((float(nel)/10.0)**0.33333333)

  ! write(*,*)'write_vtk'
  ! !$ t1 = omp_get_wtime()
  ! vtk_file = 'init.vtk'
  ! call write_vtk(vtk_file, elements)
  ! !$ t2 = omp_get_wtime()
  ! !$ write(*,*) t2 - t1

  write(*,*)'grains_to_elements'
  !$ t1 = omp_get_wtime()
  call grains_to_elements(grains,elements)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  write(*,*)'points_to_elements'
  !$ t1 = omp_get_wtime()
  call points_to_elements(elements)
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

  write(*,*)'element_moments'
  !$ t1 = omp_get_wtime()
  call element_moments(elements,boxes)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) t2 - t1

  ! vtk_file = 'test.vtk'
  ! call write_vtk(vtk_file, elements)
  ! return

  if (ifull == 0) then

    write(*,*)'element_close_inter_ids'
    !$ t1 = omp_get_wtime()
    call element_close_inter_ids(elements,boxes)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    write(*,*)'element_close_inter_Gamma'
    !$ t1 = omp_get_wtime()
    call element_close_inter_Gamma(elements)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    if (iJacobi == 1) then
      write(*,*)'calc_I_GdC_inv'
      !$ t1 = omp_get_wtime()
      call calc_I_GdC_inv(elements)
      !$ t2 = omp_get_wtime()
      !$ write(*,*) t2 - t1
    endif
  
    write(*,*)'box_far_inter_Gamma'
    !$ t1 = omp_get_wtime()
    call box_far_inter_Gamma(boxes)
    ! call box_far_inter_Gamma_old(boxes)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    write(*,*)'box_far_inter_Gamma_FFT'
    !$ t1 = omp_get_wtime()
    call box_far_inter_Gamma_FFT(boxes)
    ! call box_far_inter_Gamma_old(boxes)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    write(*,*)'form_G_approx'
    call form_G_approx(boxes)

    vtk_file = 'grid.vtk'
    call write_grid_vtk(vtk_file,elements)

  elseif (ifull == 1) then

    write(*,*) 'form_Gmat_full'
    !$ t1 = omp_get_wtime()
    ! call form_Gmat_full_intpt(elements)
    call form_Gmat_full(elements)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

    do iel1 = 1,nel
      elements(iel1)%Gamma_self = Gmat(:,:,:,:,iel1,iel1)
    enddo
    write(*,*)'calc_I_GdC_inv'
    !$ t1 = omp_get_wtime()
    call calc_I_GdC_inv(elements)
    !$ t2 = omp_get_wtime()
    !$ write(*,*) t2 - t1

  endif

!   iel1 = 1257
!   iel2 = 1132
! 
!   ib1 = 15
!   ib2 = boxes(ib1)%far_box_id(1)
!   do i = 1,boxes(ib1)%nel
!     iel1 = boxes(ib1)%element_ids(i)
!     if (norm2(elements(iel1)%dx) <= 3.0) exit
!   enddo
!   do i = 1,boxes(ib2)%nel
!     iel2 = boxes(ib2)%element_ids(i)
!     if (norm2(elements(iel2)%dx) <= 3.0) exit
!   enddo
!   write(*,*)iel1,iel2,norm2(elements(iel1)%dx),norm2(elements(iel2)%dx)
!   write(*,*)norm2(elements(iel1)%x - elements(iel2)%x)
! 
!   call integrate_G(iel1,iel2,elements,G)
!   write(*,*)'G(1,1,1,1)',G(1,1,1,1)
! 
!   call integrate_G_Taylor(iel1,iel2,elements,boxes,Gtmp)
!   write(*,*)'Gtmp(1,1,1,1)',Gtmp(1,1,1,1)
!   write(*,*)'norm2(Gtmp-G)/norm2(G)',norm2(Gtmp-G)/norm2(Goperreal(:,:,:,:,1,1,1))
! 
!   return

  iel = 100
  write(*,*) 'iel = ', iel
  do i = 1,4
   write(*,*) 'elements(iel)%xnode = ', elements(iel)%xnode(:,i)
  enddo
  write(*,*) 'elements(iel)%x = ', elements(iel)%x
  write(*,*) '<elements(iel)%x> = ', elements(iel)%dx + boxes(elements(iel)%box_id)%x
  write(*,*) 'elements(iel)%grain_id = ', elements(iel)%grain_id

  do i = 1,5
   write(*,*) 'elements(iel)%xintpt = ', elements(iel)%xintpt(:,i)
  enddo
  write(*,*) 'elements(iel)%w = ', elements(iel)%w(:)
  write(*,*) 'elements(iel)%v = ', elements(iel)%v
  write(*,*) 'elements(iel)%dx = ', elements(iel)%dx
  write(*,*) 'elements(iel)%dxdx = ', elements(iel)%dxdx

  write(*,*) 'elements(iel)%nel_close = ', elements(iel)%nel_close
  do i = 1,elements(iel)%nel_close
    write(*,*) 'elements(iel)%close_inter_el_ids(i) = ', elements(iel)%close_inter_el_ids(i)
    write(*,*) 'elements(iel)%Gamma(1,1,1,1,i) = ', elements(iel)%Gamma(1,1,1,1,i)
  enddo

  write(*,*) 'elements(iel)%Gamma_self(1,1,1,1) = ', elements(iel)%Gamma_self(1,1,1,1)
  write(*,*) 'elements(iel)%I_GdC_inv(1,1,1,1) = ', elements(iel)%I_GdC_inv(1,1,1,1)


  if (ifull == 0) then
    ib = 100
    write(*,*) 'ib = ', ib
    write(*,*) 'boxes(ib)%x = ', boxes(ib)%x
    write(*,*) 'boxes(ib)%xstart = ', boxes(ib)%xstart
    write(*,*) 'boxes(ib)%xend = ', boxes(ib)%xend
    write(*,*) 'boxes(ib)%nel = ', boxes(ib)%nel
    do i = 1,boxes(ib)%nel
     write(*,*) '  iel = ', boxes(ib)%element_ids(i)
     write(*,*) '  elements(iel)%x = ', elements(boxes(ib)%element_ids(i))%x
    enddo
    write(*,*) 'boxes(ib)%ib_grid = ', boxes(ib)%ib_grid
    write(*,*) 'boxes(ib)%nnear = ', boxes(ib)%nnear
    do i = 1,boxes(ib)%nnear
      write(*,*) '  boxes(ib)%near_box_id(i) = ', boxes(ib)%near_box_id(i)
    enddo
    write(*,*) 'boxes(ib)%nnear = ', boxes(ib)%nfar
    do i = 1,boxes(ib)%nfar
      write(*,*) '  boxes(ib)%far_box_id(i) = ', boxes(ib)%far_box_id(i)
      write(*,*) '  boxes(ib)%Gamma(1,1,1,1,i) = ', boxes(ib)%Gamma(1,1,1,1,i)
      if (order > 0) write(*,*) '  boxes(ib)%dGamma_dX(1,1,1,1,1,i) = ', boxes(ib)%dGamma_dX(1,1,1,1,1,i)
      if (order > 1) write(*,*) '  boxes(ib)%d2Gamma_dX2(1,1,1,1,1,1,i) = ', boxes(ib)%d2Gamma_dX2(1,1,1,1,1,1,i)
    enddo
  endif

  call initial_guess(elements)

  ! iterate
  iter = 0
  err = 1.0
  do while (err.gt.5.0e-4.and.iter.lt.500)
   iter = iter + 1

   ! strain as convolution
   do iel1 = 1,nel
     elements(iel1)%strainold = elements(iel1)%strain
   enddo

   if (ifull == 0) write(*,*) 'outgoing_expansion'
   !$ if (ifull == 0) t1 = omp_get_wtime()
   if (ifull == 0) call outgoing_expansion(elements, boxes)
   !$ if (ifull == 0) t2 = omp_get_wtime()
   !$ if (ifull == 0) write(*,*) t2 - t1

   if (ifull == 0) write(*,*) 'incoming_expansion'
   !$ if (ifull == 0) t1 = omp_get_wtime()
   if (FFT) then
    if (ifull == 0) call incoming_expansion_FFT(boxes)
   else
    if (ifull == 0) call incoming_expansion(elements, boxes)
   endif
   !$ if (ifull == 0) t2 = omp_get_wtime()
   !$ if (ifull == 0) write(*,*) t2 - t1

   ! write(*,*)boxes(1)%incoming0
   ! write(*,*)boxes(1)%incoming1
   ! write(*,*)boxes(1)%incoming2
   ! stop

   write(*,*) 'convolution'
   !$ t1 = omp_get_wtime()
   call convolution(elements, boxes)
   !$ t2 = omp_get_wtime()
   !$ write(*,*) t2 - t1

   err = 0.0
   do iel1 = 1,nel
     err = err + sum(elements(iel1)%strainold - elements(iel1)%strain)**2
   enddo
   err = sqrt(err)/norm2(Eapp)
   write(*,*) 'iter, err = ', iter, err

   do iel1 = 1,nel
     elements(iel1)%strain = elements(iel1)%strainold*(1.0-step) + elements(iel1)%strain*step
   enddo

   ! stress and polarization
   polaravg = 0.0
   do iel1 = 1,nel
     elements(iel1)%stress(:,:) = TijklTkl(elements(iel1)%C(:,:,:,:), elements(iel1)%strain(:,:))
     elements(iel1)%polar(:,:) = elements(iel1)%stress(:,:) - TijklTkl(c0(:,:,:,:), elements(iel1)%strain(:,:))
     polaravg = polaravg + elements(iel1)%polar(:,:)*elements(iel1)%v
   enddo
   polaravg = polaravg/vtot

  enddo

  vtk_file = 'grid.vtk'
  call write_grid_vtk(vtk_file,elements)
  ! vtk_file = 'output.vtk'
  ! call write_vtk(vtk_file, elements)

end
