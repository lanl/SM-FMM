module types
  use global
  implicit none

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
   double precision, allocatable :: Gamma(:,:,:,:,:)

   double precision :: Gamma_self(3,3,3,3) = 0
   double precision :: I_GdC_inv(3,3,3,3) = 0
   double precision :: C(3,3,3,3) = 0
   double precision :: stress(3,3) = 0
   double precision :: strain(3,3) = 0
   double precision :: strainold(3,3) = 0
   double precision :: polar(3,3) = 0
   double precision :: Q(3,3) = 0

   double precision :: xnode (3,4) = 0
   integer :: ind_node(4) = 0

   ! double precision :: xintpt(3,5) = 0
   ! double precision :: w(5) = 0
   double precision :: xintpt(3,4) = 0
   double precision :: w(4) = 0


end type element_type

type box_type

   integer, allocatable :: element_ids(:)
   integer :: nel = 0, nnear = 0, nfar = 0
   integer :: ib_grid(3) = 0
   double precision :: x(3) = 0
   double precision :: xstart(3) = 0
   double precision :: xend(3) = 0
   double precision :: edge_lngth(3) = 0
   integer, allocatable :: near_box_id(:)
   integer, allocatable :: far_box_id(:)
   integer :: ipstart(3), ipend(3)
   integer :: ngrid_ind = 0
   double precision, allocatable :: Gamma(:,:,:,:,:)
   double precision, allocatable :: dGamma_dX(:,:,:,:,:,:)
   double precision, allocatable :: d2Gamma_dX2(:,:,:,:,:,:,:)
   double precision :: outgoing0(3,3)
   double precision :: outgoing1(3,3,3)
   double precision :: outgoing2(3,3,3,3)
   double precision :: incoming0(3,3)
   double precision :: incoming1(3,3,3)
   double precision :: incoming2(3,3,3,3)
  
end type box_type


type grain_type

   integer :: id, phase
   double precision :: Q(3,3), C66(6,6), C(3,3,3,3), eul(3)

end type grain_type
  
public :: init

contains

  subroutine init(boxes)

    implicit none
    type (box_type), intent(inout) :: boxes


  end subroutine init


end
