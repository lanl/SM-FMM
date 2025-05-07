module global
  !$ use omp_lib
  implicit none
  double precision, parameter  :: pi = 4.d0*datan(1.d0)
  double complex, parameter  :: ximag = (0.d0, 1.d0)
  double precision, parameter  :: sqrt2 = sqrt(2.0)
  double precision, parameter , dimension (3,3) :: id3 = reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),shape(id3))
  double precision, parameter , dimension (6,6) :: id6 = reshape((/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 1.0, 0.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
                                                                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0/),shape(id6))
  double precision, parameter , dimension (3,3,3) :: e = reshape((/ 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, &
    0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/),shape(e))
  integer, parameter, dimension (3,3) :: tens2Mandel = reshape((/1,6,5,6,2,4,5,4,3/),shape(tens2Mandel))
  integer, parameter  :: nphmx = 3
  integer, parameter  :: nmodmx = 3  
  integer, parameter  :: nsysmx = 24
  integer, parameter  :: ntwmmx = 0
  integer, parameter  :: ngrid_ind_mx = 100000
  integer, parameter  :: nbound_mx = 100000
  integer, parameter  :: nmx_shared_faces_per_edge = 50
  integer, parameter  :: nmx_shared_elements_per_edge = 50
  integer :: npts1 = 32
  integer :: npts2 = 32
  integer :: npts3 = 32
  integer :: nbox1 = 0
  integer :: nbox2 = 0
  integer :: nbox3 = 0
  integer :: nbox = 0
  integer :: nel, nnode, nfc, nedg
  integer :: ng
  integer :: nthreads = 1
  integer :: ifull = 0
  integer :: iload = 0
  integer :: iJacobi = 0
  integer :: order = 0
  integer :: nel_distorted = 0
  integer, parameter  :: igamma = 1
  integer, parameter  :: ngmx = 100
  integer, parameter  :: nintg_pt = 4
  integer, parameter  :: nintg_pt_quad_1d = 6!  3 or 6
  integer, parameter  :: nintg_pt_quad = 36! 9 or 36
  integer :: el_face_nodes(3,4), face_edge_nodes(2,3), face_edge_nodes_remaining(3)
  double precision :: aspect_mx = 100.0
  double precision :: dx_min = 1.0
  double precision :: tolerance = 1.0
  double precision :: tolerance_out = 1.0
  integer :: itmx = 100, itmin = 3
  integer :: max_near_edges_in_box = 0
  integer :: solution_method = 2
  integer :: itout_mx = 10, iter_out
  double precision :: min_dist_node = 0.1, min_dist_node_tmp = 0.1
  double precision :: coef_reg(4), err_bc, dx_th = 2.0, err_polar
  double precision :: G_coef1, G_coef2

  integer :: iter, nclose_inter_el_ids_mx
  double precision, allocatable :: plot_el_scalar(:)
  double precision :: C(3,3,3,3), c066(6,6), s066(6,6), c0(3,3,3,3), s0(3,3,3,3)
  double precision :: lambda, mu, k_penalty, k_penalty_incr
  double precision :: Eapp(3,3), vtot, polaravg(6), err, strainavg(6), xc0, step, stressavg(6), Sapp(3,3)
  double precision :: Edotapp(3,3)
  double precision, allocatable :: Gper(:,:,:,:,:)
  double precision, allocatable :: G_box(:,:,:,:,:)
  double precision, allocatable :: dG_dX_box(:,:,:,:,:,:)
  double precision, allocatable :: d2G_dX2_box(:,:,:,:,:,:,:)
  double precision, allocatable :: d3G_dX3_box(:,:,:,:,:,:,:,:)
  double precision, allocatable :: d4G_dX4_box(:,:,:,:,:,:,:,:,:)
  double complex, allocatable :: G_box_hat(:,:,:,:,:)
  double complex, allocatable :: dG_dX_box_hat(:,:,:,:,:,:)
  double complex, allocatable :: d2G_dX2_box_hat(:,:,:,:,:,:,:)
  double complex, allocatable :: d3G_dX3_box_hat(:,:,:,:,:,:,:,:)
  double complex, allocatable :: d4G_dX4_box_hat(:,:,:,:,:,:,:,:,:)
  integer, allocatable :: el_id(:,:,:), grid_ib(:,:,:)
  double precision, allocatable :: u_edge(:,:)
  logical :: FFT = .true.

  integer, allocatable :: el_distorted(:)

  logical :: periodic = .false.
  logical :: precalc_G = .false.
  logical :: plasticity = .false.
  
  integer :: ind2d21d(3,3), ind3d21d(3,3,3), ind4d21d(3,3,3,3)
  integer :: ind1d22d(2,6), ind1d23d(3,10), ind1d24d(4,15)
  double precision :: arr1d_2d_fact(6), arr1d_3d_fact(10), arr1d_4d_fact(15)

  double precision, allocatable :: slip_normal(:,:), slip_dir(:,:), Schmid(:,:,:)
  integer :: nsys
  double precision :: gamdot0 = 1.0, tauc
  integer :: rate_exp = 10
  double precision :: time_inc
  integer :: nincr, incr, id_gas, nit_non_loc

contains

  subroutine allocate_globals

    integer :: dim1, dim2, dim3

    if (periodic) then
      dim1 = nbox1
      dim2 = nbox2
      dim3 = nbox3
    else
      dim1 = nbox1*2 - 1
      dim2 = nbox2*2 - 1
      dim3 = nbox3*2 - 1
    endif

    allocate(Gper(3,3,npts1,npts2,npts3))
    allocate(G_box(3,3,dim1,dim2,dim3))
    allocate(dG_dX_box(3,3,3,dim1,dim2,dim3))
    allocate(d2G_dX2_box(3,3,3,3,dim1,dim2,dim3))
    allocate(d3G_dX3_box(3,3,3,3,3,dim1,dim2,dim3))
    allocate(d4G_dX4_box(3,3,3,3,3,3,dim1,dim2,dim3))
    allocate(G_box_hat(3,3,dim1,dim2,dim3/2+1))
    allocate(dG_dX_box_hat(3,3,3,dim1,dim2,dim3/2+1))
    allocate(d2G_dX2_box_hat(3,3,3,3,dim1,dim2,dim3/2+1))
    allocate(d3G_dX3_box_hat(3,3,3,3,3,dim1,dim2,dim3/2+1))
    allocate(d4G_dX4_box_hat(3,3,3,3,3,3,dim1,dim2,dim3/2+1))
    allocate(el_id(npts1,npts2,npts3))
    allocate(grid_ib(nbox1,nbox2,nbox3))
    G_box = 0
    dG_dX_box = 0
    d2G_dX2_box = 0
    d3G_dX3_box = 0
    d4G_dX4_box = 0
    G_box_hat = 0
    dG_dX_box_hat = 0
    d2G_dX2_box_hat = 0
    d3G_dX3_box_hat = 0
    d4G_dX4_box_hat = 0
    el_id = 0
    grid_ib = 0

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