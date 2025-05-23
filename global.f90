module global
  !$ use omp_lib
  implicit none
  double precision, parameter  :: pi = 4.d0*datan(1.d0)
  double complex, parameter  :: ximag = (0.d0, 1.d0)
  double precision, parameter , dimension (3,3) :: id3 = reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),shape(id3))
  double precision, parameter , dimension (6,6) :: id6 = reshape((/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 1.0, 0.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0/),shape(id6))
  double precision, parameter , dimension (3,3,3) :: e = reshape((/ 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, &
    0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/),shape(e))
  integer, parameter  :: nphmx = 3
  integer, parameter  :: nmodmx = 3  
  integer, parameter  :: nsysmx = 24
  integer, parameter  :: ntwmmx = 0
  integer, parameter  :: ngrid_ind_mx = 10000
  integer :: npts1 = 32
  integer :: npts2 = 32
  integer :: npts3 = 32
  integer :: nbox1 = 0
  integer :: nbox2 = 0
  integer :: nbox3 = 0
  integer :: nbox = 0
  integer :: nel, nnode
  integer :: ng
  integer :: nthreads = 1
  integer :: ifull = 0
  integer :: iload = 0
  integer :: iJacobi = 0
  integer :: order = 0
  integer, parameter  :: igamma = 1
  integer, parameter  :: ngmx = 100

  integer :: iter, nclose_inter_el_ids_mx
  double precision, allocatable :: plot_el_scalar(:)
  double precision :: C(3,3,3,3), c066(6,6), s066(6,6), c0(3,3,3,3), s0(3,3,3,3)
  double precision :: Eapp(3,3), vtot, polaravg(3,3), err, strainavg(3,3), xc0, step
  double precision, allocatable :: Goper(:,:,:,:,:,:,:)
  double precision, allocatable :: Goper_box(:,:,:,:,:,:,:)
  double precision, allocatable :: Goperreal(:,:,:,:,:,:,:)
  double precision, allocatable :: Goperreal_approx(:,:,:,:,:,:,:)
  double precision, allocatable :: Goperreal_box(:,:,:,:,:,:,:)
  double precision, allocatable :: Gmat(:,:,:,:,:,:)
  double precision, allocatable :: Gamma_box(:,:,:,:,:,:,:)
  double precision, allocatable :: dGamma_dX_box(:,:,:,:,:,:,:,:)
  double precision, allocatable :: d2Gamma_dX2_box(:,:,:,:,:,:,:,:,:)
  double complex, allocatable :: Gamma_box_hat(:,:,:,:,:,:,:)
  double complex, allocatable :: dGamma_dX_box_hat(:,:,:,:,:,:,:,:)
  double complex, allocatable :: d2Gamma_dX2_box_hat(:,:,:,:,:,:,:,:,:)
  integer, allocatable :: el_id(:,:,:)
  logical :: FFT = .true.

contains

  subroutine allocate_globals

    allocate(Goper(3,3,3,3,npts1,npts2,npts3))
    allocate(Goper_box(3,3,3,3,nbox1,nbox2,nbox3))
    allocate(Goperreal(3,3,3,3,npts1,npts2,npts3))
    allocate(Goperreal_approx(3,3,3,3,npts1,npts2,npts3))
    allocate(Goperreal_box(3,3,3,3,nbox1,nbox2,nbox3))
    ! if (ifull == 1) allocate(Gmat(3,3,3,3,nel,nel))
    allocate(Gamma_box(3,3,3,3,nbox1,nbox2,nbox3))
    allocate(dGamma_dX_box(3,3,3,3,3,nbox1,nbox2,nbox3))
    allocate(d2Gamma_dX2_box(3,3,3,3,3,3,nbox1,nbox2,nbox3))
    allocate(Gamma_box_hat(nbox3,nbox2,nbox1,3,3,3,3))
    allocate(dGamma_dX_box_hat(nbox3,nbox2,nbox1,3,3,3,3,3))
    allocate(d2Gamma_dX2_box_hat(nbox3,nbox2,nbox1,3,3,3,3,3,3))
    allocate(el_id(npts1,npts2,npts3))
    Goper = 0
    Goper_box = 0
    Goperreal = 0
    Goperreal_approx = 0
    Goperreal_box = 0
    ! if (ifull == 1) Gmat = 0
    Gamma_box = 0
    dGamma_dX_box = 0
    d2Gamma_dX2_box = 0
    el_id = 0

  end

  subroutine setup_openmp
    !$ use omp_lib
    implicit none
    !$ if (nthreads == 0) then
    !$   nthreads=omp_get_max_threads()
    !!$   nthreads=int(nthreads/2)
    !$ endif
    !$ call omp_set_num_threads(nthreads)
    !$ write(*,'(A, I0, A)') 'Running OMP parallel on ', nthreads, ' threads'
  end subroutine setup_openmp

end module global