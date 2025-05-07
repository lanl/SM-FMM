module fourier_functions_r2c
  !$ use omp_lib
  use FFTW3
  implicit none

type fftw_type

  integer :: nbox1, nbox2, nbox3, nbox3ext, nbox3cmplx

  real(C_DOUBLE), pointer :: fourr(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc(:,:,:)
  type(C_PTR) :: plan, data, iplan
  double precision :: normfact

  real(C_DOUBLE), pointer :: fourr3(:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc3(:,:,:,:)
  type(C_PTR) :: plan3, data3, iplan3
  ! double precision :: normfact3

  real(C_DOUBLE), pointer :: fourr33(:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc33(:,:,:,:,:)
  type(C_PTR) :: plan33, data33, iplan33

  real(C_DOUBLE), pointer :: fourr333(:,:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc333(:,:,:,:,:,:)
  type(C_PTR) :: plan333, data333, iplan333

  real(C_DOUBLE), pointer :: fourr3333(:,:,:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc3333(:,:,:,:,:,:,:)
  type(C_PTR) :: plan3333, data3333, iplan3333

  real(C_DOUBLE), pointer :: fourr33333(:,:,:,:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc33333(:,:,:,:,:,:,:,:)
  type(C_PTR) :: plan33333, data33333, iplan33333

  real(C_DOUBLE), pointer :: fourr333333(:,:,:,:,:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc333333(:,:,:,:,:,:,:,:,:)
  type(C_PTR) :: plan333333, data333333, iplan333333

  real(C_DOUBLE), pointer :: fourr36(:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc36(:,:,:,:,:)
  type(C_PTR) :: plan36, data36, iplan36

  real(C_DOUBLE), pointer :: fourr310(:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc310(:,:,:,:,:)
  type(C_PTR) :: plan310, data310, iplan310

  real(C_DOUBLE), pointer :: fourr315(:,:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: fourc315(:,:,:,:,:)
  type(C_PTR) :: plan315, data315, iplan315

end type fftw_type

contains

subroutine allocate_and_plan(ffts,nbox1in,nbox2in,nbox3in,nthreads)

  type (fftw_type) :: ffts

  integer, intent(in) :: nbox1in,nbox2in,nbox3in,nthreads
  integer :: nbox1,nbox2,nbox3,nbox_shp(3)
  integer :: dimtrue(3), dimr(3), dimc(3), sizer, sizec
  !$ integer :: iret
  logical :: openmp = .false.

  
  nbox1 = 2*nbox1in - 1
  nbox2 = 2*nbox2in - 1
  nbox3 = 2*nbox3in - 1

  ffts%nbox1 = nbox1
  ffts%nbox2 = nbox2
  ffts%nbox3 = nbox3
  ffts%nbox3ext = 2*(nbox3/2+1)
  ffts%nbox3cmplx = nbox3/2+1

  !$ openmp = .true.
  !$ iret = fftw_init_threads()
  !$ if (iret==0) stop '-> Error initializing FFTW threads'
  !$ call fftw_plan_with_nthreads(nthreads)

  dimtrue = [nbox1,nbox2,nbox3]
  dimr = [nbox1,nbox2,2*(nbox3/2+1)]
  dimc = [nbox1,nbox2,(nbox3/2+1)]
  sizer = product(dimr)
  sizec = product(dimc)

  ffts%data = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1, C_SIZE_T))
  call c_f_pointer(ffts%data, ffts%fourr, [2*(nbox3/2+1),nbox2,nbox1])
  call c_f_pointer(ffts%data, ffts%fourc, [nbox3/2+1,nbox2,nbox1])
  ffts%plan = fftw_plan_dft_r2c_3d(nbox1,nbox2,nbox3, ffts%fourr, ffts%fourc, FFTW_ESTIMATE)
  ffts%iplan = fftw_plan_dft_c2r_3d(nbox1,nbox2,nbox3, ffts%fourc,ffts%fourr, FFTW_ESTIMATE)
  ffts%normfact = 1.0/float(nbox1*nbox2*nbox3)

  ffts%data3 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data3, ffts%fourr3, [2*(nbox3/2+1),nbox2,nbox1,3])
  call c_f_pointer(ffts%data3, ffts%fourc3, [nbox3/2+1,nbox2,nbox1,3])
  ! ffts%plan3 = fftw_plan_dft_r2c(4,[3,nbox1,nbox2,nbox3], ffts%fourr3, ffts%fourc3, FFTW_ESTIMATE)
  ! ffts%iplan3 = fftw_plan_dft_c2r(4,[3,nbox1,nbox2,nbox3], ffts%fourc3,ffts%fourr3, FFTW_ESTIMATE)
  ffts%plan3 = fftw_plan_many_dft_r2c(3, dimtrue, 3,  &
                  ffts%fourr3, dimr, 1, sizer,  &
                  ffts%fourc3, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan3 = fftw_plan_many_dft_c2r(3, dimtrue, 3, &
                  ffts%fourc3, dimc, 1, sizec,  &
                  ffts%fourr3, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data33 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data33, ffts%fourr33, [2*(nbox3/2+1),nbox2,nbox1,3,3])
  call c_f_pointer(ffts%data33, ffts%fourc33, [nbox3/2+1,nbox2,nbox1,3,3])
  ffts%plan33 = fftw_plan_many_dft_r2c(3, dimtrue, 3*3,  &
                  ffts%fourr33, dimr, 1, sizer,  &
                  ffts%fourc33, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan33 = fftw_plan_many_dft_c2r(3, dimtrue, 3*3, &
                  ffts%fourc33, dimc, 1, sizec,  &
                  ffts%fourr33, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data333 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 3 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data333, ffts%fourr333, [2*(nbox3/2+1),nbox2,nbox1,3,3,3])
  call c_f_pointer(ffts%data333, ffts%fourc333, [nbox3/2+1,nbox2,nbox1,3,3,3])
  ffts%plan333 = fftw_plan_many_dft_r2c(3, dimtrue, 3*3*3,  &
                  ffts%fourr333, dimr, 1, sizer,  &
                  ffts%fourc333, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan333 = fftw_plan_many_dft_c2r(3, dimtrue, 3*3*3, &
                  ffts%fourc333, dimc, 1, sizec,  &
                  ffts%fourr333, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data3333 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 3 * 3 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data3333, ffts%fourr3333, [2*(nbox3/2+1),nbox2,nbox1,3,3,3,3])
  call c_f_pointer(ffts%data3333, ffts%fourc3333, [nbox3/2+1,nbox2,nbox1,3,3,3,3])
  ffts%plan3333 = fftw_plan_many_dft_r2c(3, dimtrue, 3*3*3*3,  &
                  ffts%fourr3333, dimr, 1, sizer,  &
                  ffts%fourc3333, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan3333 = fftw_plan_many_dft_c2r(3, dimtrue, 3*3*3*3, &
                  ffts%fourc3333, dimc, 1, sizec,  &
                  ffts%fourr3333, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data33333 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 3 * 3 * 3 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data33333, ffts%fourr33333, [2*(nbox3/2+1),nbox2,nbox1,3,3,3,3,3])
  call c_f_pointer(ffts%data33333, ffts%fourc33333, [nbox3/2+1,nbox2,nbox1,3,3,3,3,3])
  ffts%plan33333 = fftw_plan_many_dft_r2c(3, dimtrue, 3*3*3*3*3,  &
                  ffts%fourr33333, dimr, 1, sizer,  &
                  ffts%fourc33333, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan33333 = fftw_plan_many_dft_c2r(3, dimtrue, 3*3*3*3*3, &
                  ffts%fourc33333, dimc, 1, sizec,  &
                  ffts%fourr33333, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data333333 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 3 * 3 * 3 * 3 * 3, C_SIZE_T))
  call c_f_pointer(ffts%data333333, ffts%fourr333333, [2*(nbox3/2+1),nbox2,nbox1,3,3,3,3,3,3])
  call c_f_pointer(ffts%data333333, ffts%fourc333333, [nbox3/2+1,nbox2,nbox1,3,3,3,3,3,3])
  ffts%plan333333 = fftw_plan_many_dft_r2c(3, dimtrue, 3*3*3*3*3*3,  &
                  ffts%fourr333333, dimr, 1, sizer,  &
                  ffts%fourc333333, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan333333 = fftw_plan_many_dft_c2r(3, dimtrue, 3*3*3*3*3*3, &
                  ffts%fourc333333, dimc, 1, sizec,  &
                  ffts%fourr333333, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data36 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 6, C_SIZE_T))
  call c_f_pointer(ffts%data36, ffts%fourr36, [2*(nbox3/2+1),nbox2,nbox1,3,6])
  call c_f_pointer(ffts%data36, ffts%fourc36, [nbox3/2+1,nbox2,nbox1,3,6])
  ffts%plan36 = fftw_plan_many_dft_r2c(3, dimtrue, 3*6,  &
                  ffts%fourr36, dimr, 1, sizer,  &
                  ffts%fourc36, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan36 = fftw_plan_many_dft_c2r(3, dimtrue, 3*6, &
                  ffts%fourc36, dimc, 1, sizec,  &
                  ffts%fourr36, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data310 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 10, C_SIZE_T))
  call c_f_pointer(ffts%data310, ffts%fourr310, [2*(nbox3/2+1),nbox2,nbox1,3,10])
  call c_f_pointer(ffts%data310, ffts%fourc310, [nbox3/2+1,nbox2,nbox1,3,10])
  ffts%plan310 = fftw_plan_many_dft_r2c(3, dimtrue, 3*10,  &
                  ffts%fourr310, dimr, 1, sizer,  &
                  ffts%fourc310, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan310 = fftw_plan_many_dft_c2r(3, dimtrue, 3*10, &
                  ffts%fourc310, dimc, 1, sizec,  &
                  ffts%fourr310, dimr, 1, sizer, FFTW_ESTIMATE)

  ffts%data315 = fftw_alloc_complex(int((nbox3/2+1) * nbox2 * nbox1 * 3 * 15, C_SIZE_T))
  call c_f_pointer(ffts%data315, ffts%fourr315, [2*(nbox3/2+1),nbox2,nbox1,3,15])
  call c_f_pointer(ffts%data315, ffts%fourc315, [nbox3/2+1,nbox2,nbox1,3,15])
  ffts%plan315 = fftw_plan_many_dft_r2c(3, dimtrue, 3*15,  &
                  ffts%fourr315, dimr, 1, sizer,  &
                  ffts%fourc315, dimc, 1, sizec, FFTW_ESTIMATE)
  ffts%iplan315 = fftw_plan_many_dft_c2r(3, dimtrue, 3*15, &
                  ffts%fourc315, dimc, 1, sizec,  &
                  ffts%fourr315, dimr, 1, sizer, FFTW_ESTIMATE)

end

subroutine test_3D

  real(C_DOUBLE), pointer :: rarr(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: carr(:,:,:)
  type(C_PTR) :: plan, data, iplan
  integer :: L = 10, M = 12, N = 14
  integer :: i, j, k

  data = fftw_alloc_complex(int((L/2+1) * M * N, C_SIZE_T))
  call c_f_pointer(data, rarr, [2*(L/2+1),M,N])
  call c_f_pointer(data, carr, [L/2+1,M,N])
  plan = fftw_plan_dft_r2c_3d(N,M,L, rarr,carr, FFTW_ESTIMATE)
  iplan = fftw_plan_dft_c2r_3d(N,M,L, carr,rarr, FFTW_ESTIMATE)

  write(*,*) shape(rarr)
  read(*,*)
  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        rarr(i,j,k) = float(i**2*j+k**2)
        write(*,*) rarr(i,j,k)
      enddo
    enddo
  enddo
  do i = L+1,2*(L/2+1)
    do j = 1,M
      do k = 1,N                         
        rarr(i,j,k) = 0.0
      enddo
    enddo
  enddo
 
  call fftw_execute_dft_r2c(plan, rarr, carr)

  do i = 1,(L/2+1)
    do j = 1,M
      do k = 1,N                         
        write(*,*) carr(i,j,k)
      enddo
    enddo
  enddo

  call fftw_execute_dft_c2r(iplan, carr, rarr)

  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        write(*,*) rarr(i,j,k)/float(L*M*N)
      enddo
    enddo
  enddo

  call fftw_destroy_plan(plan)
  call fftw_destroy_plan(iplan)
  call fftw_free(data)

end
  
subroutine test_nD

  real(C_DOUBLE), pointer :: rarr(:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: carr(:,:,:,:)
  type(C_PTR) :: plan, data, iplan
  integer :: L = 11, M = 13, N = 15
  integer :: i, j, k, ii

  data = fftw_alloc_complex(int((L/2+1) * M * N * 3, C_SIZE_T))
  call c_f_pointer(data, rarr, [2*(L/2+1),M,N,3])
  call c_f_pointer(data, carr, [L/2+1,M,N,3])
  plan = fftw_plan_dft_r2c(4,[3,N,M,L], rarr,carr, FFTW_ESTIMATE)
  iplan = fftw_plan_dft_c2r(4,[3,N,M,L], carr,rarr, FFTW_ESTIMATE)

  write(*,*) shape(rarr)
  read(*,*)
  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        do ii = 1,3                         
          rarr(i,j,k,ii) = float(i**2*j+k**2-ii*(i+j+k))
        enddo
        write(*,*) rarr(i,j,k,:)
      enddo
    enddo
  enddo
  do i = L+1,2*(L/2+1)
    do j = 1,M
      do k = 1,N                         
        rarr(i,j,k,:) = 0.0
      enddo
    enddo
  enddo
  read(*,*)
 
  call fftw_execute_dft_r2c(plan, rarr, carr)

  do i = 1,(L/2+1)
    do j = 1,M
      do k = 1,N                         
        write(*,*) carr(i,j,k,:)
      enddo
    enddo
  enddo

  call fftw_execute_dft_c2r(iplan, carr, rarr)

  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        write(*,*) rarr(i,j,k,:)/float(L*M*N*3)
      enddo
    enddo
  enddo
  read(*,*)

  call fftw_destroy_plan(plan)
  call fftw_destroy_plan(iplan)
  call fftw_free(data)

end


subroutine test_nD_many

  real(C_DOUBLE), pointer :: rarr(:,:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: carr(:,:,:,:)
  type(C_PTR) :: plan, data, iplan
  integer :: L = 35, M = 23, N = 19
  integer :: i, j, k, ii

  data = fftw_alloc_complex(int((L/2+1) * M * N * 3, C_SIZE_T))
  call c_f_pointer(data, rarr, [2*(L/2+1),M,N,3])
  call c_f_pointer(data, carr, [L/2+1,M,N,3])
  ! plan = fftw_plan_dft_r2c(4,[3,N,M,L], rarr,carr, FFTW_ESTIMATE)
  ! iplan = fftw_plan_dft_c2r(4,[3,N,M,L], carr,rarr, FFTW_ESTIMATE)

  plan = fftw_plan_many_dft_r2c(3, [N,M,L], 3,  &
                  rarr, [N,M,2*(L/2+1)], 1, N*M*2*(L/2+1),  &
                  carr, [N,M,(L/2+1)], 1, N*M*(L/2+1), FFTW_ESTIMATE)
  iplan = fftw_plan_many_dft_c2r(3, [N,M,L], 3, &
                  carr, [N,M,(L/2+1)], 1, N*M*(L/2+1),  &
                  rarr, [N,M,2*(L/2+1)], 1, N*M*2*(L/2+1), FFTW_ESTIMATE)

  write(*,*) shape(rarr)
  read(*,*)
  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        do ii = 1,3                         
          rarr(i,j,k,ii) = float(i**2*j+k**2-ii*(i+j+k))
        enddo
        write(*,*) rarr(i,j,k,:)
      enddo
    enddo
  enddo
  do i = L+1,2*(L/2+1)
    do j = 1,M
      do k = 1,N                         
        rarr(i,j,k,:) = 0.0
      enddo
    enddo
  enddo
  read(*,*)
 
  call fftw_execute_dft_r2c(plan, rarr, carr)

  do i = 1,1! (L/2+1)
    do j = 1,1! M
      do k = 1,N                         
        write(*,*) carr(i,j,k,:)
      enddo
    enddo
  enddo

  call fftw_execute_dft_c2r(iplan, carr, rarr)

  do i = 1,L
    do j = 1,M
      do k = 1,N                         
        write(*,*) rarr(i,j,k,:)/float(L*M*N)
      enddo
    enddo
  enddo
  read(*,*)

  call fftw_destroy_plan(plan)
  call fftw_destroy_plan(iplan)
  call fftw_free(data)

end

subroutine test_fft_arr(ffts)

  type (fftw_type) :: ffts
  integer :: i, j, k, ii

  do i = 1,ffts%nbox3
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                    
        do ii = 1,3                         
          ffts%fourr3(i,j,k,ii) = float(i**2*j+k**2-ii*(i+j+k))
        enddo
        write(*,*) ffts%fourr3(i,j,k,:)
      enddo
    enddo
  enddo

  call fftw_execute_dft_r2c(ffts%plan3, ffts%fourr3, ffts%fourc3)

  do i = 1,1! ffts%nbox3cmplx
    do j = 1,1! ffts%nbox2
      do k = 1,ffts%nbox1                    
        write(*,*) ffts%fourc3(i,j,k,:)
      enddo
    enddo
  enddo

  call fftw_execute_dft_c2r(ffts%iplan3, ffts%fourc3, ffts%fourr3)
  
  do i = 1,ffts%nbox3
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                            
        write(*,*) ffts%fourr3(i,j,k,:)*ffts%normfact
      enddo
    enddo
  enddo

end

subroutine test_conv_arr(ffts)

  type (fftw_type) :: ffts
  integer :: i, j, k, ii

  do i = 1,ffts%nbox3
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                    
        do ii = 1,3                         
          ffts%fourr3(i,j,k,ii) = float(i**2*j+k**2-ii*(i+j+k))
        enddo
        ffts%fourr(i,j,k) = 1.0/(float(i**2+j**2+k**2))
      enddo
    enddo
  enddo

  call fftw_execute_dft_r2c(ffts%plan3, ffts%fourr3, ffts%fourc3)
  call fftw_execute_dft_r2c(ffts%plan, ffts%fourr, ffts%fourc)

  do i = 1,ffts%nbox3cmplx
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                    
        ffts%fourc3(i,j,k,:) = ffts%fourc3(i,j,k,:)*ffts%fourc(i,j,k)
      enddo
    enddo
  enddo

  do i = 1,ffts%nbox3cmplx
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                    
        write(777,*) i,j,k,ffts%fourc3(i,j,k,:)
      enddo
    enddo
  enddo

  call fftw_execute_dft_c2r(ffts%iplan3, ffts%fourc3, ffts%fourr3)

  do i = 1,ffts%nbox3
    do j = 1,ffts%nbox2
      do k = 1,ffts%nbox1                            
        write(888,*) i,j,k,ffts%fourr3(i,j,k,:)*ffts%normfact
      enddo
    enddo
  enddo

end

end module fourier_functions_r2c
