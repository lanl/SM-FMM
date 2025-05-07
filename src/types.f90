module types
  ! use global
  implicit none

type intg_pt_type
  double precision :: Xhat(2)
  double precision :: Xhatn
  double precision :: w
end type intg_pt_type

type edge_type

  double precision :: xc(3) = 0
  double precision :: x(3,2) = 0
  integer, allocatable :: face(:)
  integer, allocatable :: face_close(:)
  integer, allocatable :: face_edge(:)
  double precision, allocatable :: face_G(:,:,:)
  double precision, allocatable :: face_close_G(:,:,:)
  integer, allocatable :: neigh_edges(:)
  integer, allocatable :: element(:)
  integer, allocatable :: box_id(:)
  integer, allocatable :: box_id_loc(:)
  integer :: box_id_unique
  integer :: nface = 0
  integer :: nel = 0
  double precision :: f(3) = 0
  double precision :: u(3) = 0
  double precision, allocatable :: dx(:,:)
  double precision, allocatable :: dxdx(:,:,:)
  double precision, allocatable :: G(:,:,:,:)

end type edge_type

type face_type

  double precision :: x(3) = 0
  integer :: el(2) = 0
  integer :: el_face(2) = 0
  integer :: nel = 0
  logical :: boundary = .false.
  integer :: edge(3) = 0
  double precision :: normal(3) = 0 ! outward to el 1
  double precision :: F(3) = 0
  double precision :: Fbc(3) = 0
  double precision :: area
  double precision :: u(3) = 0
  integer :: box_id(2)
  integer :: box_id_loc(2)
  integer :: nbox

end type face_type

type element_type

   double precision :: v = 0
   double precision :: x(3) = 0
   double precision :: dx(3) = 0
   double precision :: dxdx(3,3) = 0
   integer, allocatable :: grid_ind(:,:)
   double precision, allocatable :: x_grid(:,:)
   integer, allocatable :: close_inter_el_ids(:)
   integer :: nel_close = 0
   integer :: ngrid_ind = 0
   integer :: grain_id = 0
   integer :: box_id = 0
   integer :: face(4) = 0
   double precision :: face_normal_sgn(4) = 0
   double precision, allocatable :: Gamma(:,:,:,:,:)

   double precision :: Gamma_self(3,3,3,3) = 0
   double precision :: I_GdC_inv(3,3,3,3) = 0
   double precision :: I_C0S_inv(6,6) = 0
   double precision :: C(6,6) = 0
   double precision :: S(6,6) = 0
   double precision :: stress(6) = 0
   double precision :: stresst(6) = 0
   double precision :: strain(6) = 0
   double precision :: straint(6) = 0
   double precision :: disgrad(3,3) = 0
   double precision :: disgradt(3,3) = 0
   double precision :: strainold(6) = 0
   double precision :: polar(6) = 0
   double precision :: polarold(6) = 0
   double precision :: Q(3,3) = 0
   double precision :: eigen(6) = 0
   double precision :: strain_pl(6) = 0
   double precision :: edotp(6) = 0
   double precision :: strain_pl_vm = 0
   double precision :: strain_hyd_non_loc = 0
   double precision, allocatable :: Schmid(:,:)
   double precision, allocatable :: tauc(:)
   double precision, allocatable :: gamdot(:)

   double precision :: xnode (3,4) = 0
   integer :: ind_node(4) = 0

   ! double precision :: xintpt(3,5) = 0
   ! double precision :: w(5) = 0
   double precision :: xintpt(3,4) = 0
   double precision :: w(4) = 0

   double precision, allocatable :: x_bound(:,:)
   double precision, allocatable :: norm_bound(:,:)
   integer :: np_bound = 0

   integer, allocatable :: iel_interp(:)
   double precision, allocatable :: eig_dir(:,:), eval(:)
   double precision :: aspect

end type element_type

type box_type

   integer, allocatable :: element_ids(:)
   integer, allocatable :: edge_ids(:)
   integer, allocatable :: face_ids(:)
   integer, allocatable :: edge_unique_ids(:)
   integer, allocatable :: near_edge_ids(:)
   integer, allocatable :: near_edge_unique_ids(:)
   integer, allocatable :: near_bound_face_ids(:)
   double precision, allocatable :: edge_near_factor(:)
   double precision, allocatable :: face_F(:,:), edge_F(:,:), near_edge_f(:,:)
   integer :: nel = 0, nedg = 0, nnear = 0, nfar = 0, nnear_edge = 0, nnear_bound_face = 0
   integer :: ib_grid(3) = 0, nfc = 0
   double precision :: x(3) = 0
   double precision :: xstart(3) = 0
   double precision :: xend(3) = 0
   double precision :: edge_lngth(3) = 0
   integer, allocatable :: near_box_id(:)
   integer, allocatable :: far_box_id(:)
   integer :: ipstart(3), ipend(3)
   integer :: ngrid_ind = 0

   double precision, allocatable :: G(:,:,:)
   double precision, allocatable :: dG_dX(:,:,:,:)
   double precision, allocatable :: d2G_dX2(:,:,:,:,:)
   double precision, allocatable :: d3G_dX3(:,:,:,:,:,:)
   double precision, allocatable :: d4G_dX4(:,:,:,:,:,:,:)
   double precision :: outgoing0f(3)
   double precision :: outgoing1f(3,3)
   double precision :: outgoing2f(3,3,3)
   double precision :: outgoing3f(3,3,3,3)
   double precision :: outgoing4f(3,3,3,3,3)
   double precision :: incoming0f(3)
   double precision :: incoming1f(3,3)
   double precision :: incoming2f(3,3,3)
   double precision :: incoming3f(3,3,3,3)
   double precision :: incoming4f(3,3,3,3,3)

   double precision :: outgoing0f_1d(3)
   double precision :: outgoing1f_1d(3,3)
   double precision :: outgoing2f_1d(3,6)
   double precision :: outgoing3f_1d(3,10)
   double precision :: outgoing4f_1d(3,15)
   double precision :: incoming0f_1d(3)
   double precision :: incoming1f_1d(3,3)
   double precision :: incoming2f_1d(3,6)
   double precision :: incoming3f_1d(3,10)
   double precision :: incoming4f_1d(3,15)
  
end type box_type

type grain_type

   integer :: id, phase
   double precision :: Q(3,3), C66(6,6), C(3,3,3,3), eul(3), S(3,3,3,3), S66(6,6)
   double precision :: eigenstrain(3,3) = 0.0, strain(6) = 0.0, v = 0.0

end type grain_type

end
