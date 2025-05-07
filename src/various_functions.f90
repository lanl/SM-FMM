module various_functions

  implicit none

contains

  subroutine grains_to_elements(grains, elements)
    use global
    use types
    use tensor_functions

    type (grain_type), allocatable, intent (in) :: grains(:)
    type (element_type), allocatable, intent (inout) :: elements(:)

    integer :: iel, id

    do iel = 1, nel

      id = elements(iel)%grain_id
      elements(iel)%C = mandel_t4_to_t2(grains(id)%C)
      elements(iel)%S = mandel_t4_to_t2(grains(id)%S)
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
    vtot = 0.0
    do iel = 1, nel

      elements(iel)%strain = mandel_t2_to_t1(Eapp)
      elements(iel)%stress = matmul(elements(iel)%C, elements(iel)%strain - elements(iel)%eigen)
      elements(iel)%polar = elements(iel)%stress - elements(iel)%stresst - &
       matmul(c066, elements(iel)%strain - elements(iel)%straint)
      polaravg = polaravg + elements(iel)%polar*elements(iel)%v
      vtot = vtot + elements(iel)%v

    enddo
    polaravg = polaravg/vtot

    do iel = 1,nel
      elements(iel)%polar = elements(iel)%polar - polaravg
    enddo
 
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

  double precision :: xmin(3), xmax(3), x(3)
  double precision, allocatable :: x_grid(:,:)
  integer, allocatable :: grid_ind(:,:)
  integer :: iel, i, gridmin(3), gridmax(3), ip1, ip2, ip3, ic

  allocate(grid_ind(3,ngrid_ind_mx))
  allocate(x_grid(3,ngrid_ind_mx))

  vtot = 0.0
  do iel = 1, nel

    do i = 1,3
      xmin(i) = minval(elements(iel)%xnode(i,:))
      xmax(i) = maxval(elements(iel)%xnode(i,:))
    enddo

    gridmin = floor(xmin)
    gridmax = ceiling(xmax)

    ic = 0
    do ip1 = gridmin(1),gridmax(1)
      do ip2 = gridmin(2),gridmax(2)
        do ip3 = gridmin(3),gridmax(3)
          x(:) = float([ip1,ip2,ip3])
          if (point_in_tetra(x,elements(iel)%xnode)) then

            ic = ic + 1
            x_grid(:,ic) = x
            grid_ind(1,ic) = enforce_periodic(ip1, 1, npts1)
            grid_ind(2,ic) = enforce_periodic(ip2, 1, npts2)
            grid_ind(3,ic) = enforce_periodic(ip3, 1, npts3)

            if (el_id(grid_ind(1,ic),grid_ind(2,ic),grid_ind(3,ic)) .ne. 0) then 
               write(*,*) grid_ind(:,ic)
               write(*,*) el_id(grid_ind(1,ic),grid_ind(2,ic),grid_ind(3,ic)), iel
               stop
            endif

            el_id(grid_ind(1,ic),grid_ind(2,ic),grid_ind(3,ic)) = iel

          endif
        enddo
      enddo
    enddo

    ! if (ic == 0) write(*,*) 'element has no points', iel
    elements(iel)%ngrid_ind = ic
    allocate(elements(iel)%grid_ind(3,ic))
    allocate(elements(iel)%x_grid(3,ic))
    elements(iel)%grid_ind(:,1:ic) = grid_ind(:,1:ic)
    elements(iel)%x_grid(:,1:ic) = x_grid(:,1:ic)

    ! elements(iel)%v = float(elements(iel)%ngrid_ind)
    vtot = vtot + elements(iel)%v

  enddo

  deallocate(grid_ind)
  deallocate(x_grid)

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
     grid_ib(ib1,ib2,ib3) = ib

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
           xc(3).ge.boxes(ib)%xstart(3).and.xc(3).lt.boxes(ib)%xend(3)) then
           ! xc(3).ge.boxes(ib)%xstart(3).and.xc(3).lt.boxes(ib)%xend(3).and.elements(iel)%v > 0.0) then
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

subroutine edges_faces_to_boxes(edges,faces,elements,boxes)
  use global
  use types
  use tensor_functions

  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: iedg, ifc, iel, i, j, ic, nun, ib, ic_near, ib_near, iel1, iel2, ic_box_near, ic1, ic2, ictot
  integer, allocatable :: all_el_ids(:), unique_el_ids(:)
  integer, allocatable :: all_box_ids(:), unique_box_ids(:)
  integer, allocatable :: edge_ids(:), edge_ids_near(:), face_ids(:)
  integer, allocatable :: face_ids_near_bound(:)
  integer, allocatable :: edge_box_id_unique(:)
  double precision, allocatable :: edge_near_factor(:)
  double precision :: xc(3)
  logical :: flag
  logical, allocatable :: edges_flag(:)

  allocate(all_el_ids(nmx_shared_elements_per_edge))
  allocate(unique_el_ids(nmx_shared_elements_per_edge))
  allocate(all_box_ids(nmx_shared_elements_per_edge))
  allocate(unique_box_ids(nmx_shared_elements_per_edge))
  allocate(edges_flag(nedg))
  allocate(edge_box_id_unique(nedg))

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,ic,i,ifc,iel,all_el_ids,unique_el_ids,nun,all_box_ids,unique_box_ids) &
  !$OMP SHARED(nedg,edges,faces,elements) 
  !$OMP DO 
  do iedg = 1,nedg

    ic = 0
    do i = 1,edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      do j = 1,faces(ifc)%nel
        iel = faces(ifc)%el(j)
        ic = ic + 1
        if (ic > nmx_shared_elements_per_edge) stop 'ic > nmx_shared_elements_per_edge'
        all_el_ids(ic) = iel
      enddo
    enddo

    call uniqueInt(all_el_ids, ic, unique_el_ids, nun, nmx_shared_elements_per_edge)
    allocate(edges(iedg)%element(nun))
    edges(iedg)%element = unique_el_ids(1:nun)
    edges(iedg)%nel = nun

    do i = 1,edges(iedg)%nel
      iel = edges(iedg)%element(i)
      all_box_ids(i) = elements(iel)%box_id
    enddo
    call uniqueInt(all_box_ids, edges(iedg)%nel, unique_box_ids, nun, nmx_shared_elements_per_edge)
    allocate(edges(iedg)%box_id(nun))
    edges(iedg)%box_id = unique_box_ids(1:nun)
    allocate(edges(iedg)%box_id_loc(nun))

    ! write(*,*) iedg
    ! write(*,*) edges(iedg)%element
    ! write(*,*) edges(iedg)%box_id

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(all_el_ids)
  deallocate(unique_el_ids)
  deallocate(all_box_ids)
  deallocate(unique_box_ids)

  max_near_edges_in_box = 0
  allocate(edge_ids(nedg))
  allocate(edge_ids_near(nedg))
  allocate(edge_near_factor(nedg))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ic,ic_near,iedg,flag,ic_box_near,ib_near) &
  !$OMP PRIVATE(edge_ids,edge_ids_near,edge_near_factor) &
  !$OMP SHARED(nbox,nedg,edges,boxes)  &
  !$OMP REDUCTION(max:max_near_edges_in_box)
  !$OMP DO 
  do ib = 1, nbox

    ic = 0
    ic_near = 0
    do iedg = 1,nedg

      flag = .false.
      ic_box_near = 0
      do i = 1,size(edges(iedg)%box_id) ! loop over all boxes this edge belons to

        if (edges(iedg)%box_id(i) == ib) then 
          ic = ic + 1
          edge_ids(ic) = iedg
          edges(iedg)%box_id_loc(i) = ic
        endif

        ! if (.not.flag) then
          do j = 1,boxes(ib)%nnear ! loop over all near boxes
            ib_near = boxes(ib)%near_box_id(j)
            if (edges(iedg)%box_id(i) == ib_near) then
              ic_box_near = ic_box_near + 1 ! count to how many near boxes this edge belongs to
              if (.not.flag) then
                flag = .true.
                ic_near = ic_near + 1
                edge_ids_near(ic_near) = iedg
              endif
            endif
          enddo
        ! endif

      enddo
      if (flag) edge_near_factor(ic_near) = float(ic_box_near)/float(size(edges(iedg)%box_id)) ! factor for force (if fully in near domain 1.0, less than 1.0 if at boundary)

    enddo
    boxes(ib)%nedg = ic
    allocate(boxes(ib)%edge_ids(ic))
    boxes(ib)%edge_ids(1:ic) = edge_ids(1:ic)
    allocate(boxes(ib)%edge_f(3,ic))
    boxes(ib)%nnear_edge = ic_near
    allocate(boxes(ib)%near_edge_ids(ic_near))
    boxes(ib)%near_edge_ids(1:ic_near) = edge_ids_near(1:ic_near)
    allocate(boxes(ib)%edge_near_factor(ic_near))
    boxes(ib)%edge_near_factor(1:ic_near) = edge_near_factor(1:ic_near)
    allocate(boxes(ib)%near_edge_f(3,ic_near))
    boxes(ib)%near_edge_f = 0.0

    max_near_edges_in_box = maxval([max_near_edges_in_box,ic_near])

    ! write(*,*) ib
    ! write(*,*) boxes(ib)%edge_ids
    ! write(*,*) boxes(ib)%edge_near_factor
    ! read(*,*)

  enddo 
  !$OMP END DO
  !$OMP END PARALLEL
  write(*,*)'max_near_edges_in_box = ', max_near_edges_in_box

  ! unique
  ictot = 0
  edges_flag = .false.
  do ib = 1, nbox

    ic = 0
    do iedg = 1,nedg
      xc = edges(iedg)%xc
      if (.not.edges_flag(iedg)) then
        if (xc(1).ge.boxes(ib)%xstart(1)-1.0e-7.and.xc(1).lt.boxes(ib)%xend(1)+1.0e-7.and. &
            xc(2).ge.boxes(ib)%xstart(2)-1.0e-7.and.xc(2).lt.boxes(ib)%xend(2)+1.0e-7.and. &
            xc(3).ge.boxes(ib)%xstart(3)-1.0e-7.and.xc(3).lt.boxes(ib)%xend(3)+1.0e-7) then
          ic = ic + 1
          ! if (edges_flag(iedg)) then
          !   write(*,*) 'already assigned to a box' 
          !   write(*,*) edges(iedg)%xc
          !   write(*,*) edge_box_id_unique(iedg)
          !   write(*,*) boxes(edge_box_id_unique(iedg))%xstart
          !   write(*,*) boxes(edge_box_id_unique(iedg))%xend
          !   write(*,*) ib
          !   write(*,*) boxes(ib)%xstart
          !   write(*,*) boxes(ib)%xend
          !   read(*,*)
          ! endif
          edges_flag(iedg) = .true.
          edge_box_id_unique(iedg) = ib
          edges(iedg)%box_id_unique = ib
          ictot = ictot + 1
          edge_ids(ic) = iedg
        endif
      endif
    enddo

    allocate(boxes(ib)%edge_unique_ids(ic))
    boxes(ib)%edge_unique_ids = edge_ids(1:ic)

  enddo 
  write(*,*) ictot, nedg
  if (ictot .ne. nedg) stop 'ictot .ne. nedg in edges_faces_to_boxes'

  do ib = 1, nbox

    ic1 = 0
    ic2 = 0
    do i = 1, boxes(ib)%nnear
      ib_near = boxes(ib)%near_box_id(i)
      ic1 = ic2 + 1
      ic2 = ic2 + size(boxes(ib_near)%edge_unique_ids,1)
      edge_ids_near(ic1:ic2) = boxes(ib_near)%edge_unique_ids
    enddo
    allocate(boxes(ib)%near_edge_unique_ids(ic2))
    boxes(ib)%near_edge_unique_ids = edge_ids_near(1:ic2)

  enddo 

  ! ib = 1
  ! do i = 1,size(boxes(ib)%near_edge_ids,1)
  !   if (i <= size(boxes(ib)%near_edge_unique_ids,1)) then
  !     write(*,*)boxes(ib)%near_edge_ids(i), boxes(ib)%near_edge_unique_ids(i)
  !   else
  !     write(*,*)boxes(ib)%near_edge_ids(i)
  !   endif
  ! enddo

  deallocate(edge_ids)
  deallocate(edge_ids_near)
  deallocate(edge_near_factor)
  deallocate(edges_flag)
  deallocate(edge_box_id_unique)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,iel,iel1,iel2) &
  !$OMP SHARED(nfc,faces,elements) 
  !$OMP DO 
  do ifc = 1,nfc
 
   if (faces(ifc)%nel == 1) then
     iel = faces(ifc)%el(1)
     faces(ifc)%nbox = 1
     faces(ifc)%box_id(1) = elements(iel)%box_id
   elseif (faces(ifc)%nel == 2) then
     iel1 = faces(ifc)%el(1)
     iel2 = faces(ifc)%el(2)
     if (elements(iel1)%box_id == elements(iel2)%box_id) then
       faces(ifc)%nbox = 1
       faces(ifc)%box_id(1) = elements(iel1)%box_id
     else
       faces(ifc)%nbox = 2
       faces(ifc)%box_id(1) = elements(iel1)%box_id
       faces(ifc)%box_id(2) = elements(iel2)%box_id
     endif
   endif
 
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  allocate(face_ids_near_bound(nfc))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ic,ifc,ic_near,ib_near,face_ids_near_bound) &
  !$OMP SHARED(nbox,nfc,faces,boxes) 
  !$OMP DO 
  do ib = 1, nbox
    ic = 0
    do ifc = 1,nfc
  
      if (faces(ifc)%nbox == 2) then

        ic_near = 0
        do j = 1,boxes(ib)%nnear
          ib_near = boxes(ib)%near_box_id(j)
          if (faces(ifc)%box_id(1) == ib_near) then
            ic_near = ic_near + 1
            exit
          endif
        enddo

        do j = 1,boxes(ib)%nnear
          ib_near = boxes(ib)%near_box_id(j)
          if (faces(ifc)%box_id(2) == ib_near) then
            ic_near = ic_near + 1
            exit
          endif
        enddo

        if (ic_near == 1) then
          ic = ic + 1
          face_ids_near_bound(ic) = ifc
        endif

      endif
  
    enddo
    boxes(ib)%nnear_bound_face = ic
    allocate(boxes(ib)%near_bound_face_ids(ic))
    boxes(ib)%near_bound_face_ids(1:ic) = face_ids_near_bound(1:ic)

  enddo 
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(face_ids_near_bound)

  allocate(face_ids(nfc))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ic,ifc,ic_near,ib_near,face_ids,i) &
  !$OMP SHARED(nbox,nfc,faces,boxes) 
  !$OMP DO 
  do ib = 1, nbox
    ic = 0
    do ifc = 1,nfc
  
      do i = 1, faces(ifc)%nbox
        if (faces(ifc)%box_id(i) == ib) then
          ic = ic + 1
          face_ids(ic) = ifc
          faces(ifc)%box_id_loc(i) = ic
        endif
      enddo
  
    enddo
    boxes(ib)%nfc = ic
    allocate(boxes(ib)%face_ids(ic))
    boxes(ib)%face_ids(1:ic) = face_ids(1:ic)
    allocate(boxes(ib)%face_F(3,ic))

  enddo 
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(face_ids)

  do iedg = 1,nedg

    do i = 1, size(edges(iedg)%box_id)
      ib = edges(iedg)%box_id(i)
      ! write(*,*) iedg, boxes(ib)%edge_ids(edges(iedg)%box_id_loc(i))
      if (iedg .ne. boxes(ib)%edge_ids(edges(iedg)%box_id_loc(i))) then
        read(*,*)
      endif

    enddo

  enddo
  
  do ifc = 1,nfc

    do i = 1, faces(ifc)%nbox
      ib = faces(ifc)%box_id(i)
      ! write(*,*) ifc, boxes(ib)%face_ids(faces(ifc)%box_id_loc(i))
      if (ifc .ne. boxes(ib)%face_ids(faces(ifc)%box_id_loc(i))) then
        read(*,*)
      endif

    enddo

  enddo

  ! ib = 1
  ! do i = 1,boxes(ib)%nfc
  !   ifc = boxes(ib)%face_ids(i)
  !   write(*,*) ifc, faces(ifc)%box_id
  ! enddo
  ! stop

  ! ib = 10
  ! do i = 1,boxes(ib)%nnear_bound_face
  !   ifc = boxes(ib)%near_bound_face_ids(i)
  !   write(*,*) faces(ifc)%box_id(1:2)
  !   write(*,*) boxes(ib)%near_box_id
  ! enddo

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

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,icnear,icfar,ib2,dib,boxnear,boxfar) &
  !$OMP SHARED(nbox,boxes,periodic,nbox1,nbox2,nbox3) 
  !$OMP DO 
  do ib1 = 1,nbox

    icnear = 0
    icfar = 0
    do ib2 = 1,nbox

      dib = abs(boxes(ib1)%ib_grid - boxes(ib2)%ib_grid)
      if (dib(1).eq.nbox1 - 1 .and. periodic) dib(1) = 1
      if (dib(2).eq.nbox2 - 1 .and. periodic) dib(2) = 1
      if (dib(3).eq.nbox3 - 1 .and. periodic) dib(3) = 1

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
  !$OMP END DO
  !$OMP END PARALLEL

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
      if (ic > 0) then
        allocate(elements(iel1)%close_inter_el_ids(ic))
        elements(iel1)%close_inter_el_ids(1:ic) = close_inter_el_ids(1:ic)
      endif
    enddo
  enddo

  deallocate(close_inter_el_ids)

  write(*,*) 'nclose_inter_el_ids_mx', nclose_inter_el_ids_mx

end

subroutine box_far_inter_G(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dx(3), dir1_inc(3), dir2_inc(3), dir3_inc(3), inc
  integer :: ib, i, ibfar, ind(3), ic, ip1, ip2, ip3, ip1far, ip2far, ip3far

  do ib = 1,nbox
    allocate(boxes(ib)%G(3,3,boxes(ib)%nfar))
    if (order > 0) allocate(boxes(ib)%dG_dX(3,3,3,boxes(ib)%nfar))
    if (order > 1) allocate(boxes(ib)%d2G_dX2(3,3,3,3,boxes(ib)%nfar))
    if (order > 2) allocate(boxes(ib)%d3G_dX3(3,3,3,3,3,boxes(ib)%nfar))
    if (order > 3) allocate(boxes(ib)%d4G_dX4(3,3,3,3,3,3,boxes(ib)%nfar))
  enddo

  ! far field interactions
  ! inc = float(nint(boxes(1)%edge_lngth(1)/2.0))/1.0! 1.0
  inc = 1.0e-3
  if (periodic) inc = 1.0
  write(*,*) 'inc for derivatives of Gamma is ', inc
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,i, ibfar, dx, ind) &
  !$OMP SHARED(nbox, boxes, inc, dir1_inc, dir2_inc, dir3_inc,order, lambda, mu, periodic)
  !$OMP DO 
  do ib = 1,nbox
    ! write(*,*) ib
    do i = 1,boxes(ib)%nfar

      ibfar = boxes(ib)%far_box_id(i)
      ! write(*,*) ib,ibfar
      dx = - boxes(ibfar)%x + boxes(ib)%x
        
      boxes(ib)%G(:,:,i) = Green_function(dx, lambda, mu)
  
      if (order > 0) then
        boxes(ib)%dG_dX(:,:,1,i) = (Green_function(dx + dir1_inc, lambda, mu) - &
          Green_function(dx - dir1_inc, lambda, mu))/(2.0*inc)
        boxes(ib)%dG_dX(:,:,2,i) = (Green_function(dx + dir2_inc, lambda, mu) - &
          Green_function(dx - dir2_inc, lambda, mu))/(2.0*inc)
        boxes(ib)%dG_dX(:,:,3,i) = (Green_function(dx + dir3_inc, lambda, mu) - &
          Green_function(dx - dir3_inc, lambda, mu))/(2.0*inc)
      endif
 
      if (order > 1) then
!         boxes(ib)%d2G_dX2(:,:,1,1,i) = (Green_function(dx + dir1_inc, lambda, mu) &
!          - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir1_inc, lambda, mu))/(inc**2)
!         boxes(ib)%d2G_dX2(:,:,2,2,i) = (Green_function(dx + dir2_inc, lambda, mu) &
!          - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir2_inc, lambda, mu))/(inc**2)
!         boxes(ib)%d2G_dX2(:,:,3,3,i) = (Green_function(dx + dir3_inc, lambda, mu) &
!          - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir3_inc, lambda, mu))/(inc**2)
!   
!         boxes(ib)%d2G_dX2(:,:,1,2,i) = &
!          (Green_function(dx + dir1_inc + dir2_inc, lambda, mu) - Green_function(dx + dir1_inc - dir2_inc, lambda, mu) - &
!           Green_function(dx - dir1_inc + dir2_inc, lambda, mu) + Green_function(dx - dir1_inc - dir2_inc, lambda, mu))/(4.0*inc**2)
!   
!         boxes(ib)%d2G_dX2(:,:,1,3,i) = &
!          (Green_function(dx + dir1_inc + dir3_inc, lambda, mu) - Green_function(dx + dir1_inc - dir3_inc, lambda, mu) - &
!           Green_function(dx - dir1_inc + dir3_inc, lambda, mu) + Green_function(dx - dir1_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
!         
!         boxes(ib)%d2G_dX2(:,:,2,3,i) = &
!          (Green_function(dx + dir2_inc + dir3_inc, lambda, mu) - Green_function(dx + dir2_inc - dir3_inc, lambda, mu) - &
!           Green_function(dx - dir2_inc + dir3_inc, lambda, mu) + Green_function(dx - dir2_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
!   
!         boxes(ib)%d2G_dX2(:,:,2,1,i) = boxes(ib)%d2G_dX2(:,:,1,2,i)
!         boxes(ib)%d2G_dX2(:,:,3,1,i) = boxes(ib)%d2G_dX2(:,:,1,3,i)
!         boxes(ib)%d2G_dX2(:,:,3,2,i) = boxes(ib)%d2G_dX2(:,:,2,3,i)
        boxes(ib)%d2G_dX2(:,:,:,:,i) = Green_operator(dx, lambda, mu)
      endif

      if (order > 2) then
        boxes(ib)%d3G_dX3(:,:,:,:,1,i) = (Green_operator(dx + dir1_inc, lambda, mu) - &
          Green_operator(dx - dir1_inc, lambda, mu))/(2.0*inc)
        boxes(ib)%d3G_dX3(:,:,:,:,2,i) = (Green_operator(dx + dir2_inc, lambda, mu) - &
          Green_operator(dx - dir2_inc, lambda, mu))/(2.0*inc)
        boxes(ib)%d3G_dX3(:,:,:,:,3,i) = (Green_operator(dx + dir3_inc, lambda, mu) - &
          Green_operator(dx - dir3_inc, lambda, mu))/(2.0*inc)
      endif

      if (order > 3) then
        boxes(ib)%d4G_dX4(:,:,:,:,1,1,i) = (Green_operator(dx + dir1_inc, lambda, mu) &
         - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir1_inc, lambda, mu))/(inc**2)
        boxes(ib)%d4G_dX4(:,:,:,:,2,2,i) = (Green_operator(dx + dir2_inc, lambda, mu) &
         - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir2_inc, lambda, mu))/(inc**2)
        boxes(ib)%d4G_dX4(:,:,:,:,3,3,i) = (Green_operator(dx + dir3_inc, lambda, mu) &
         - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir3_inc, lambda, mu))/(inc**2)
   
        boxes(ib)%d4G_dX4(:,:,:,:,1,2,i) = &
         (Green_operator(dx + dir1_inc + dir2_inc, lambda, mu) - Green_operator(dx + dir1_inc - dir2_inc, lambda, mu) - &
          Green_operator(dx - dir1_inc + dir2_inc, lambda, mu) + Green_operator(dx - dir1_inc - dir2_inc, lambda, mu))/(4.0*inc**2)
   
        boxes(ib)%d4G_dX4(:,:,:,:,1,3,i) = &
         (Green_operator(dx + dir1_inc + dir3_inc, lambda, mu) - Green_operator(dx + dir1_inc - dir3_inc, lambda, mu) - &
          Green_operator(dx - dir1_inc + dir3_inc, lambda, mu) + Green_operator(dx - dir1_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
        
        boxes(ib)%d4G_dX4(:,:,:,:,2,3,i) = &
         (Green_operator(dx + dir2_inc + dir3_inc, lambda, mu) - Green_operator(dx + dir2_inc - dir3_inc, lambda, mu) - &
          Green_operator(dx - dir2_inc + dir3_inc, lambda, mu) + Green_operator(dx - dir2_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
   
        boxes(ib)%d4G_dX4(:,:,:,:,2,1,i) = boxes(ib)%d4G_dX4(:,:,:,:,1,2,i)
        boxes(ib)%d4G_dX4(:,:,:,:,3,1,i) = boxes(ib)%d4G_dX4(:,:,:,:,1,3,i)
        boxes(ib)%d4G_dX4(:,:,:,:,3,2,i) = boxes(ib)%d4G_dX4(:,:,:,:,2,3,i)
      endif

      ! boxes(ib)%dG_dX = 0.0
      ! boxes(ib)%d2G_dX2 = 0.0
  
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine box_far_inter_G_FFT(boxes, ffts)
  use global
  use types
  use tensor_functions
  use fourier_functions_r2c

  type (box_type), allocatable, intent (inout) :: boxes(:)
  type (fftw_type), intent (inout) :: ffts

  double precision :: dx(3), dir1_inc(3), dir2_inc(3), dir3_inc(3), inc, xorigin(3), edge_lngth(3), dum
  double precision :: Goper3333(3,3,3,3)
  integer :: ib, i, ibfar, ind(3), j, k, l, m, n, ib1, ib2, ib3, ibfartmp
  integer :: ib1shift, ib2shift, ib3shift

  inc = 1.0e-3
  write(*,*) 'inc for derivatives of Gamma is ', inc
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc

  edge_lngth(1) = (float(npts1)/float(nbox1)) ! length of each box edge (domain 0.5 to npts + 0.5)
  edge_lngth(2) = (float(npts2)/float(nbox2)) 
  edge_lngth(3) = (float(npts3)/float(nbox3)) 
  xorigin = 0.5 + edge_lngth*0.5 ! center of box (1,1,1) of true domain
  do ib3 = 1,nbox3*2 - 1
    do ib2 = 1,nbox2*2 - 1
      do ib1 = 1,nbox1*2 - 1

        ! ib1shift = ib1 - nbox1
        ! ib2shift = ib2 - nbox2
        ! ib3shift = ib3 - nbox3
        ! dx(1) = - float(ib1shift)*edge_lngth(1)
        ! dx(2) = - float(ib2shift)*edge_lngth(2)
        ! dx(3) = - float(ib3shift)*edge_lngth(3)
        ib1shift = ib1 - 1
        ib2shift = ib2 - 1
        ib3shift = ib3 - 1
        if (ib1 > nbox1) ib1shift = ib1shift - (nbox1*2 - 1)
        if (ib2 > nbox2) ib2shift = ib2shift - (nbox2*2 - 1)
        if (ib3 > nbox3) ib3shift = ib3shift - (nbox3*2 - 1)
        ! write(*,*) ib1, ib1shift
        ! read(*,*)
        dx(1) = float(ib1shift)*edge_lngth(1)
        dx(2) = float(ib2shift)*edge_lngth(2)
        dx(3) = float(ib3shift)*edge_lngth(3)

        if (abs(ib1shift) > 1 .or. abs(ib2shift) > 1 .or. abs(ib3shift) > 1) then
          G_box(:,:,ib1,ib2,ib3) = Green_function(dx, lambda, mu)

          if (order > 0) then
            dG_dX_box(:,:,1,ib1,ib2,ib3) = (Green_function(dx + dir1_inc, lambda, mu) - &
              Green_function(dx - dir1_inc, lambda, mu))/(2.0*inc)
            dG_dX_box(:,:,2,ib1,ib2,ib3) = (Green_function(dx + dir2_inc, lambda, mu) - &
              Green_function(dx - dir2_inc, lambda, mu))/(2.0*inc)
            dG_dX_box(:,:,3,ib1,ib2,ib3) = (Green_function(dx + dir3_inc, lambda, mu) - &
              Green_function(dx - dir3_inc, lambda, mu))/(2.0*inc)
          endif

          if (order > 1) then
!            d2G_dX2_box(:,:,1,1,ib1,ib2,ib3) = (Green_function(dx + dir1_inc, lambda, mu) &
!             - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir1_inc, lambda, mu))/(inc**2)
!            d2G_dX2_box(:,:,2,2,ib1,ib2,ib3) = (Green_function(dx + dir2_inc, lambda, mu) &
!             - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir2_inc, lambda, mu))/(inc**2)
!            d2G_dX2_box(:,:,3,3,ib1,ib2,ib3) = (Green_function(dx + dir3_inc, lambda, mu) &
!             - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir3_inc, lambda, mu))/(inc**2)
!      
!            d2G_dX2_box(:,:,1,2,ib1,ib2,ib3) = &
!             (Green_function(dx + dir1_inc + dir2_inc, lambda, mu) - &
!              Green_function(dx + dir1_inc - dir2_inc, lambda, mu) - &
!              Green_function(dx - dir1_inc + dir2_inc, lambda, mu) + &
!              Green_function(dx - dir1_inc - dir2_inc, lambda, mu))/(4.0*inc**2)
!      
!            d2G_dX2_box(:,:,1,3,ib1,ib2,ib3) = &
!             (Green_function(dx + dir1_inc + dir3_inc, lambda, mu) - &
!              Green_function(dx + dir1_inc - dir3_inc, lambda, mu) - &
!              Green_function(dx - dir1_inc + dir3_inc, lambda, mu) + &
!              Green_function(dx - dir1_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
!            
!            d2G_dX2_box(:,:,2,3,ib1,ib2,ib3) = &
!             (Green_function(dx + dir2_inc + dir3_inc, lambda, mu) - &
!              Green_function(dx + dir2_inc - dir3_inc, lambda, mu) - &
!              Green_function(dx - dir2_inc + dir3_inc, lambda, mu) + &
!              Green_function(dx - dir2_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
!      
!            d2G_dX2_box(:,:,2,1,ib1,ib2,ib3) = d2G_dX2_box(:,:,1,2,ib1,ib2,ib3)
!            d2G_dX2_box(:,:,3,1,ib1,ib2,ib3) = d2G_dX2_box(:,:,1,3,ib1,ib2,ib3)
!            d2G_dX2_box(:,:,3,2,ib1,ib2,ib3) = d2G_dX2_box(:,:,2,3,ib1,ib2,ib3)

            ! Goper3333 = Green_operator(dx, lambda, mu)
            d2G_dX2_box(:,:,:,:,ib1,ib2,ib3) =  Green_operator(dx, lambda, mu)
            ! do i = 1,3
            !   do j = 1,3
            !     do k = 1,3
            !       do l = 1,3
            !         write(*,*) Goper3333(i,j,k,l), d2G_dX2_box(i,j,k,l,ib1,ib2,ib3)
            !       enddo
            !     enddo
            !   enddo
            ! enddo
            ! stop

          endif

          if (order > 2) then
            d3G_dX3_box(:,:,:,:,1,ib1,ib2,ib3) = (Green_operator(dx + dir1_inc, lambda, mu) - &
              Green_operator(dx - dir1_inc, lambda, mu))/(2.0*inc)
            d3G_dX3_box(:,:,:,:,2,ib1,ib2,ib3) = (Green_operator(dx + dir2_inc, lambda, mu) - &
              Green_operator(dx - dir2_inc, lambda, mu))/(2.0*inc)
            d3G_dX3_box(:,:,:,:,3,ib1,ib2,ib3) = (Green_operator(dx + dir3_inc, lambda, mu) - &
              Green_operator(dx - dir3_inc, lambda, mu))/(2.0*inc)
          endif

          if (order > 3) then
            d4G_dX4_box(:,:,:,:,1,1,ib1,ib2,ib3) = (Green_operator(dx + dir1_inc, lambda, mu) &
             - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir1_inc, lambda, mu))/(inc**2)
            d4G_dX4_box(:,:,:,:,2,2,ib1,ib2,ib3) = (Green_operator(dx + dir2_inc, lambda, mu) &
             - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir2_inc, lambda, mu))/(inc**2)
            d4G_dX4_box(:,:,:,:,3,3,ib1,ib2,ib3) = (Green_operator(dx + dir3_inc, lambda, mu) &
             - 2.0*Green_operator(dx, lambda, mu) + Green_operator(dx - dir3_inc, lambda, mu))/(inc**2)
       
            d4G_dX4_box(:,:,:,:,1,2,ib1,ib2,ib3) = &
             (Green_operator(dx + dir1_inc + dir2_inc, lambda, mu) - Green_operator(dx + dir1_inc - dir2_inc, lambda, mu) - &
              Green_operator(dx - dir1_inc + dir2_inc, lambda, mu) + Green_operator(dx - dir1_inc - dir2_inc, lambda, mu)) &
              /(4.0*inc**2)
       
            d4G_dX4_box(:,:,:,:,1,3,ib1,ib2,ib3) = &
             (Green_operator(dx + dir1_inc + dir3_inc, lambda, mu) - Green_operator(dx + dir1_inc - dir3_inc, lambda, mu) - &
              Green_operator(dx - dir1_inc + dir3_inc, lambda, mu) + Green_operator(dx - dir1_inc - dir3_inc, lambda, mu)) &
              /(4.0*inc**2)
            
            d4G_dX4_box(:,:,:,:,2,3,ib1,ib2,ib3) = &
             (Green_operator(dx + dir2_inc + dir3_inc, lambda, mu) - Green_operator(dx + dir2_inc - dir3_inc, lambda, mu) - &
              Green_operator(dx - dir2_inc + dir3_inc, lambda, mu) + Green_operator(dx - dir2_inc - dir3_inc, lambda, mu)) &
              /(4.0*inc**2)
       
            d4G_dX4_box(:,:,:,:,2,1,ib1,ib2,ib3) = d4G_dX4_box(:,:,:,:,1,2,ib1,ib2,ib3)
            d4G_dX4_box(:,:,:,:,3,1,ib1,ib2,ib3) = d4G_dX4_box(:,:,:,:,1,3,ib1,ib2,ib3)
            d4G_dX4_box(:,:,:,:,3,2,ib1,ib2,ib3) = d4G_dX4_box(:,:,:,:,2,3,ib1,ib2,ib3)
          endif

        else
          G_box(:,:,ib1,ib2,ib3) = 0.0
          if (order > 0) then
            dG_dX_box(:,:,:,ib1,ib2,ib3) = 0.0
          endif
          if (order > 1) then
            d2G_dX2_box(:,:,:,:,ib1,ib2,ib3) = 0.0
          endif
          if (order > 2) then
            d3G_dX3_box(:,:,:,:,:,ib1,ib2,ib3) = 0.0
          endif
          if (order > 2) then
            d4G_dX4_box(:,:,:,:,:,:,ib1,ib2,ib3) = 0.0
          endif
          ! write(*,*)ib1,ib2,ib3,ib1shift,ib2shift,ib3shift
        endif
        ! write(*,*)ib1,ib2,ib3,dx
        ! write(*,*)ib1,ib2,ib3,norm2(dG_dX_box(:,:,:,ib1,ib2,ib3))

      enddo
    enddo
  enddo
  ! stop

  do ib3 = 1,nbox3*2 - 1
    do ib2 = 1,nbox2*2 - 1
      do ib1 = 1,nbox1*2 - 1
        do i = 1,3
          do j = 1,3
            ffts%fourr33(ib3,ib2,ib1,i,j) = G_box(i,j,ib1,ib2,ib3)
            if (order > 0) then
              do k = 1,3
                ffts%fourr333(ib3,ib2,ib1,i,j,k) = dG_dX_box(i,j,k,ib1,ib2,ib3)
                if (order > 1) then
                  do l = 1,3
                    ffts%fourr3333(ib3,ib2,ib1,i,j,k,l) = d2G_dX2_box(i,j,k,l,ib1,ib2,ib3)
                    if (order > 2) then
                      do m = 1,3
                        ffts%fourr33333(ib3,ib2,ib1,i,j,k,l,m) = d3G_dX3_box(i,j,k,l,m,ib1,ib2,ib3)
                        if (order > 3) then
                          do n = 1,3
                            ffts%fourr333333(ib3,ib2,ib1,i,j,k,l,m,n) = d4G_dX4_box(i,j,k,l,m,n,ib1,ib2,ib3)
                          enddo
                        endif
                      enddo
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
        enddo
        ! write(*,*) ib1,ib2,ib3
        ! write(*,*)norm2(G_box(:,:,ib1,ib2,ib3)),norm2(dG_dX_box(:,:,:,ib1,ib2,ib3)), &
        ! norm2(d2G_dX2_box(:,:,:,:,ib1,ib2,ib3)),norm2(d3G_dX3_box(:,:,:,:,:,ib1,ib2,ib3))
      enddo
    enddo
  enddo

  call fftw_execute_dft_r2c(ffts%plan33, ffts%fourr33, ffts%fourc33)
  if (order > 0) call fftw_execute_dft_r2c(ffts%plan333, ffts%fourr333, ffts%fourc333)
  if (order > 1) call fftw_execute_dft_r2c(ffts%plan3333, ffts%fourr3333, ffts%fourc3333)
  if (order > 2) call fftw_execute_dft_r2c(ffts%plan33333, ffts%fourr33333, ffts%fourc33333)
  if (order > 3) call fftw_execute_dft_r2c(ffts%plan333333, ffts%fourr333333, ffts%fourc333333)

  G_box_hat = 0.0
  do ib3 = 1,ffts%nbox3cmplx
    do ib2 = 1,ffts%nbox2
      do ib1 = 1,ffts%nbox1
        G_box_hat(:,:,ib1,ib2,ib3) = ffts%fourc33(ib3,ib2,ib1,:,:)
        if (order > 0) dG_dX_box_hat(:,:,:,ib1,ib2,ib3) = ffts%fourc333(ib3,ib2,ib1,:,:,:)
        if (order > 1) d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3) = ffts%fourc3333(ib3,ib2,ib1,:,:,:,:)
        if (order > 2) d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3) = ffts%fourc33333(ib3,ib2,ib1,:,:,:,:,:)
        if (order > 3) d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3) = ffts%fourc333333(ib3,ib2,ib1,:,:,:,:,:,:)
      enddo
    enddo
  enddo

  ! stop

  ! ! test FFT
  ! call ifft_tensor33box(iplan_advanced33box, fourgrid33box)
  ! if (order > 0) call ifft_tensor333box(iplan_advanced333box, fourgrid333box)
  ! if (order > 1) call ifft_tensor3333box(iplan_advanced3333box, fourgrid3333box)
  ! dum = 0.0
  ! do ib3 = 1,nbox3*2 - 1
  !   do ib2 = 1,nbox2*2 - 1
  !     do ib1 = 1,nbox1*2 - 1
  !       do i = 1,3
  !         do j = 1,3
  !           dum = dum + abs(real(fourgrid33box(ib3,ib2,ib1,j,i)) - G_box(i,j,ib1,ib2,ib3))
  !           if (order > 0) then
  !             do k = 1,3
  !               dum = dum + abs(real(fourgrid333box(ib3,ib2,ib1,k,j,i)) - dG_dX_box(i,j,k,ib1,ib2,ib3))
  !               if (order > 1) then
  !                 do l = 1,3
  !                   dum = dum + abs(real(fourgrid3333box(ib3,ib2,ib1,l,k,j,i)) - d2G_dX2_box(i,j,k,l,ib1,ib2,ib3))
  !                 enddo
  !               endif
  !             enddo
  !           endif
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! write(*,*) dum
  ! stop

end

subroutine outgoing_expansion_f(edges,faces,elements,boxes)
  use global
  use types
  use tensor_functions

  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: wedg = 1.0/3.0
  double precision :: dx(3), fact, dxdx(3,3), dxdxdx(3,3,3), f(3), dxdxdxdx(3,3,3,3), polar(3,3)
  integer :: ib, iedg_in_box, i, j, iedg, ind, k, ipt
  integer :: iel_in_box, iel, l, m, ifc_in_box, iedg_in_face, ifc

  ! outgoing expansions
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,iedg_in_box,iedg,i,j,fact,dx,dxdx,ind,dxdxdx,ifc,iel,dxdxdxdx,f,polar) &
  !$OMP SHARED(boxes,edges,nbox,order,periodic,elements,faces,wedg)
  !$OMP DO 
  do ib = 1,nbox

    boxes(ib)%outgoing0f = 0.0
    boxes(ib)%outgoing1f = 0.0
    boxes(ib)%outgoing2f = 0.0
    boxes(ib)%outgoing3f = 0.0
    boxes(ib)%outgoing4f = 0.0

    do ifc_in_box = 1,boxes(ib)%nfc
      ifc = boxes(ib)%face_ids(ifc_in_box)
      if (faces(ifc)%nel == 1) then
        iel = faces(ifc)%el(1)
        do iedg_in_face = 1,3
          iedg = faces(ifc)%edge(iedg_in_face)
          dx = edges(iedg)%xc - boxes(ib)%x
          do i = 1,3
            do j = 1,3
              dxdx(i,j) = dx(i)*dx(j)
              do k = 1,3
                dxdxdx(i,j,k) = dxdx(i,j)*dx(k)
                do l = 1,3
                  dxdxdxdx(i,j,k,l) = dxdxdx(i,j,k)*dx(l)
                enddo
              enddo
            enddo
          enddo

          polar = mandel_t1_to_t2(elements(iel)%polar)
          f = (matmul(polar,faces(ifc)%normal)*faces(ifc)%area + faces(ifc)%Fbc)*wedg
          boxes(ib)%outgoing0f = boxes(ib)%outgoing0f + f
          if (order > 0) then
            do i = 1,3
              boxes(ib)%outgoing1f(:,i) = boxes(ib)%outgoing1f(:,i) +  f*dx(i)
              if (order > 1) then
                do j = 1,3
                  boxes(ib)%outgoing2f(:,i,j) = boxes(ib)%outgoing2f(:,i,j) + f*dxdx(i,j)
                  if (order > 2) then
                    do k = 1,3
                      boxes(ib)%outgoing3f(:,i,j,k) = boxes(ib)%outgoing3f(:,i,j,k) + f*dxdxdx(i,j,k)
                      if (order > 3) then
                        do l = 1,3
                          boxes(ib)%outgoing4f(:,i,j,k,l) = boxes(ib)%outgoing4f(:,i,j,k,l) + f*dxdxdxdx(i,j,k,l)
                        enddo
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo

      endif
    enddo

!    do iedg_in_box = 1,boxes(ib)%nedg
!    ! do iedg_in_box = 1,size(boxes(ib)%edge_unique_ids) ! unique
!      iedg = boxes(ib)%edge_ids(iedg_in_box)
!      ! iedg = boxes(ib)%edge_unique_ids(iedg_in_box)
!      ! fact = 1.0/float(size(edges(iedg)%box_id,1))
!      ! dx = edges(iedg)%xc - boxes(ib)%x
!      ! if (periodic) dx = per_dx(dx)
!      ind = findloc(edges(iedg)%box_id, ib, 1)
!      ! if (ind == 0) then
!      !   write(*,*) 'ind == 0'
!      !   write(*,*) 'edges(iedg)%box_id', ib
!      ! endif
!      ! dx = edges(iedg)%dx(:,ind)
!      ! dxdx = edges(iedg)%dxdx(:,:,ind)
!      dx = edges(iedg)%xc - boxes(ib)%x
!      do i = 1,3
!        do j = 1,3
!          dxdx(i,j) = dx(i)*dx(j)
!          do k = 1,3
!            dxdxdx(i,j,k) = dxdx(i,j)*dx(k)
!          enddo
!        enddo
!      enddo
!      ! boxes(ib)%outgoing0f = boxes(ib)%outgoing0f + edges(iedg)%f! *fact
!      boxes(ib)%outgoing0f = boxes(ib)%outgoing0f + boxes(ib)%edge_f(:,iedg_in_box)
!      if (order > 0) then
!        do i = 1,3
!          ! boxes(ib)%outgoing1f(:,i) = boxes(ib)%outgoing1f(:,i) + edges(iedg)%f*dx(i)! *fact
!          boxes(ib)%outgoing1f(:,i) = boxes(ib)%outgoing1f(:,i) +  boxes(ib)%edge_f(:,iedg_in_box)*dx(i)! *fact
!          if (order > 1) then
!            do j = 1,3
!              ! boxes(ib)%outgoing2f(:,i,j) = boxes(ib)%outgoing2f(:,i,j) + edges(iedg)%f*dxdx(i,j)! *fact
!              boxes(ib)%outgoing2f(:,i,j) = boxes(ib)%outgoing2f(:,i,j) + boxes(ib)%edge_f(:,iedg_in_box)*dxdx(i,j)! *fact
!              if (order > 2) then
!                do k = 1,3
!                  ! boxes(ib)%outgoing3f(:,i,j,k) = boxes(ib)%outgoing3f(:,i,j,k) + edges(iedg)%f*dxdxdx(i,j,k)! *fact
!                  boxes(ib)%outgoing3f(:,i,j,k) = boxes(ib)%outgoing3f(:,i,j,k) + boxes(ib)%edge_f(:,iedg_in_box)*dxdxdx(i,j,k)! *fact
!                enddo
!              endif
!            enddo
!          endif
!        enddo
!      endif
!    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,iel_in_box,iel,i,j,k,l,dx,ipt,m,polar) &
  !$OMP SHARED(boxes,elements,nbox,order)
  !$OMP DO 
  do ib = 1,nbox

    ! boxes(ib)%outgoing0 = 0.0
    ! boxes(ib)%outgoing1 = 0.0
    ! boxes(ib)%outgoing2 = 0.0
    do iel_in_box = 1,size(boxes(ib)%element_ids) ! unique
      iel = boxes(ib)%element_ids(iel_in_box)
      polar = mandel_t1_to_t2(elements(iel)%polar)

      do ipt = 1,4

      dx = elements(iel)%xintpt(:,ipt) - boxes(ib)%x
      if (order > 0) then
        boxes(ib)%outgoing1f = boxes(ib)%outgoing1f - polar*elements(iel)%v*elements(iel)%w(ipt) ! *fact
        if (order > 1) then
          do i = 1,3
            do j = 1,3
              do k = 1,3
                boxes(ib)%outgoing2f(i,j,k) = boxes(ib)%outgoing2f(i,j,k) - 2.0*polar(i,j)*dx(k)*elements(iel)%v* &
                  elements(iel)%w(ipt)
                if (order > 2) then
                  do l = 1,3
                    boxes(ib)%outgoing3f(i,j,k,l) = boxes(ib)%outgoing3f(i,j,k,l) - &
                     3.0*polar(i,j)*dx(k)*dx(l)*elements(iel)%v*elements(iel)%w(ipt)
                     if (order > 3) then
                       do m = 1,3
                         boxes(ib)%outgoing4f(i,j,k,l,m) = boxes(ib)%outgoing4f(i,j,k,l,m) - &
                          4.0*polar(i,j)*dx(k)*dx(l)*dx(m)*elements(iel)%v*elements(iel)%w(ipt)
                       enddo
                     endif
                  enddo
                endif
              enddo
            enddo
          enddo
        endif
      endif

      enddo

    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,i,j,k,l) &
  !$OMP SHARED(nbox,boxes,order,ind2d21d,ind3d21d,ind4d21d)
  !$OMP DO 
  do ib = 1,nbox
    boxes(ib)%outgoing0f_1d = boxes(ib)%outgoing0f
    if (order > 0) boxes(ib)%outgoing1f_1d = boxes(ib)%outgoing1f

    if (order > 1) then
      boxes(ib)%outgoing2f_1d = 0.0
      do i = 1,3
        do j = 1,3
          boxes(ib)%outgoing2f_1d(:,ind2d21d(i,j)) = boxes(ib)%outgoing2f_1d(:,ind2d21d(i,j)) + &
            boxes(ib)%outgoing2f(:,i,j)
        enddo
      enddo
    endif

    if (order > 2) then
      boxes(ib)%outgoing3f_1d = 0.0
      do i = 1,3
        do j = 1,3
          do k = 1,3
            boxes(ib)%outgoing3f_1d(:,ind3d21d(i,j,k)) =  &
              boxes(ib)%outgoing3f_1d(:,ind3d21d(i,j,k)) + boxes(ib)%outgoing3f(:,i,j,k)
          enddo
        enddo
      enddo
    endif

    if (order > 3) then
      boxes(ib)%outgoing4f_1d = 0.0
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              boxes(ib)%outgoing4f_1d(:,ind4d21d(i,j,k,l)) =  &
                boxes(ib)%outgoing4f_1d(:,ind4d21d(i,j,k,l)) + boxes(ib)%outgoing4f(:,i,j,k,l)
            enddo
          enddo
        enddo
      enddo
    endif

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine incoming_expansion_f(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: dum
  integer :: ib, ibfar, i, ii, jj, kk

  ! incoming expansions
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ibfar,i) &
  !$OMP SHARED(nbox,boxes,order)
  !$OMP DO 
  do ib = 1,nbox

    boxes(ib)%incoming0f_1d = 0.0
    boxes(ib)%incoming1f_1d = 0.0
    boxes(ib)%incoming2f_1d = 0.0
    boxes(ib)%incoming3f_1d = 0.0
    boxes(ib)%incoming4f_1d = 0.0
    do i = 1,boxes(ib)%nfar
      ibfar = boxes(ib)%far_box_id(i)

      if (order == 0) then

        boxes(ib)%incoming0f_1d = boxes(ib)%incoming0f_1d + &
          matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f_1d)
  
      elseif (order == 1) then

        boxes(ib)%incoming0f_1d = boxes(ib)%incoming0f_1d + &
          matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f_1d) -&
          TijkTjk_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f_1d)
  
        boxes(ib)%incoming1f_1d = boxes(ib)%incoming1f_1d + &
         TijkTj_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f_1d)

      elseif (order == 2) then

        boxes(ib)%incoming0f_1d = boxes(ib)%incoming0f_1d + &
          matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f_1d) -&
          TijkTjk_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.5*TijklTjkl_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f_1d)

!         write(*,*) TijkTjk(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f)
!         write(*,*) TijkTjk_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f_1d)
!         read(*,*)
! 
!         write(*,*) 0.5*TijklTjkl(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f)
!         write(*,*) 0.5*TijklTjkl_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f_1d)
!         read(*,*)

        boxes(ib)%incoming1f_1d = boxes(ib)%incoming1f_1d + &
          TijkTj_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          TijklTjk_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f_1d)

!        write(*,*) TijkTj(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f)
!        write(*,*) TijkTj_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f_1d)
!        read(*,*)
!
!        write(*,*) TijklTjk(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f)
!        write(*,*) TijklTjk_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f_1d)
!        read(*,*)
  
        boxes(ib)%incoming2f_1d = boxes(ib)%incoming2f_1d + &
          0.5*TijklTj_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f_1d)

!        write(*,*) 0.5*TijklTj(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f)
!        write(*,*)
!        write(*,*) 0.5*TijklTj_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f_1d)
!        stop
  
      elseif (order == 3) then

        boxes(ib)%incoming0f_1d = boxes(ib)%incoming0f_1d + &
          matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f_1d) -&
          TijkTjk_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.5*TijklTjkl_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f_1d) -&
          0.16666666666666*TijklmTjklm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing3f_1d)

        boxes(ib)%incoming1f_1d = boxes(ib)%incoming1f_1d + &
          TijkTj_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          TijklTjk_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.5*TijklmTjlm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing2f_1d)

        boxes(ib)%incoming2f_1d = boxes(ib)%incoming2f_1d + &
          0.5*TijklTj_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          0.5*TijklmTjm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing1f_1d)

        boxes(ib)%incoming3f_1d = boxes(ib)%incoming3f_1d + &
          0.16666666666666*TijklmTj_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing0f_1d)

      elseif (order == 4) then

        ! stop 'zeroth order term needed for b.c. case'
        boxes(ib)%incoming0f_1d = boxes(ib)%incoming0f_1d + &
          matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f_1d) -&
          TijkTjk_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.5*TijklTjkl_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f_1d) -&
          0.16666666666666*TijklmTjklm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing3f_1d) + &
          1.0/24.0*TijklmnTjklmn_1d(boxes(ib)%d4G_dX4(:,:,:,:,:,:,i),boxes(ibfar)%outgoing4f_1d)

        boxes(ib)%incoming1f_1d = boxes(ib)%incoming1f_1d + &
          TijkTj_1d(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          TijklTjk_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.5*TijklmTjlm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing2f_1d) - &
          0.16666666666666*TijklmnTjlmn_1d(boxes(ib)%d4G_dX4(:,:,:,:,:,:,i),boxes(ibfar)%outgoing3f_1d)

        boxes(ib)%incoming2f_1d = boxes(ib)%incoming2f_1d + &
          0.5*TijklTj_1d(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          0.5*TijklmTjm_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing1f_1d) + &
          0.25*TijklmnTjmn_1d(boxes(ib)%d4G_dX4(:,:,:,:,:,:,i),boxes(ibfar)%outgoing2f_1d)

        boxes(ib)%incoming3f_1d = boxes(ib)%incoming3f_1d + &
          0.16666666666666*TijklmTj_1d(boxes(ib)%d3G_dX3(:,:,:,:,:,i),boxes(ibfar)%outgoing0f_1d) - &
          0.16666666666666*TijklmnTjn_1d(boxes(ib)%d4G_dX4(:,:,:,:,:,:,i),boxes(ibfar)%outgoing1f_1d)

        boxes(ib)%incoming4f_1d = boxes(ib)%incoming4f_1d + &
          1.0/24.0*TijklmnTj_1d(boxes(ib)%d4G_dX4(:,:,:,:,:,:,i),boxes(ibfar)%outgoing0f_1d)

      endif

    enddo
    ! write(*,*) boxes(ib)%incoming0f
    ! write(*,*) boxes(ib)%incoming1f
    ! write(*,*) boxes(ib)%incoming2f
    ! write(*,*) boxes(ib)%incoming3f
    ! ! write(*,*) boxes(ib)%incoming4f
    ! read(*,*)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! write(*,*) boxes(1)%incoming1f
  ! stop

end

subroutine incoming_expansion_f_per(boxes)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: ib, ibfar, i

  ! incoming expansions
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ibfar,i) &
  !$OMP SHARED(nbox,boxes,order)
  !$OMP DO 
  do ib = 1,nbox

    boxes(ib)%incoming0f = 0.0
    boxes(ib)%incoming1f = 0.0
    boxes(ib)%incoming2f = 0.0
    do i = 1,boxes(ib)%nfar
      ibfar = boxes(ib)%far_box_id(i)

      if (order == 0) then

        ! boxes(ib)%incoming0f = boxes(ib)%incoming0f + &
        !   matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f)
  
      elseif (order == 1) then

        ! boxes(ib)%incoming0f = boxes(ib)%incoming0f + &
        !   matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f) -&
        !   TijkTjk(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f)
  
        boxes(ib)%incoming1f = boxes(ib)%incoming1f + &
         TijkTj(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f)

      elseif (order == 2) then

        ! boxes(ib)%incoming0f = boxes(ib)%incoming0f + &
        !   matmul(boxes(ib)%G(:,:,i),boxes(ibfar)%outgoing0f) -&
        !   TijkTjk(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing1f) + &
        !   0.5*TijklTjkl(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing2f)
  
        boxes(ib)%incoming1f = boxes(ib)%incoming1f + &
         TijkTj(boxes(ib)%dG_dX(:,:,:,i),boxes(ibfar)%outgoing0f) - &
         TijklTjk(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing1f)
  
        boxes(ib)%incoming2f = boxes(ib)%incoming2f + &
          0.5*TijklTj(boxes(ib)%d2G_dX2(:,:,:,:,i),boxes(ibfar)%outgoing0f)
  
      endif

    enddo
    ! write(*,*) boxes(ib)%incoming0f
    ! write(*,*) boxes(ib)%incoming1f
    ! write(*,*) boxes(ib)%incoming2f
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! write(*,*) boxes(1)%incoming1f
  ! stop

end

subroutine incoming_expansion_f_FFT_1d(boxes, ffts)
  use global
  use types
  use tensor_functions
  use fourier_functions_r2c

  type (box_type), allocatable, intent (inout) :: boxes(:)
  type (fftw_type), intent (inout) :: ffts

  integer :: ib, i, ib1, ib2, ib3, j, k, l, m
  double precision :: dum, aux3(3), aux33(3,3), t1, t2, dt
  double complex :: aux33c(3,3), dumc, aux333c(3,3,3), outgoing0(3), outgoing1(3,3)
  double complex :: outgoing2(3,6), aux3333c(3,3,3,3), aux3c(3), outgoing3(3,10)
  double complex :: outgoing4(3,15), aux33333c(3,3,3,3,3)
  double complex :: aux36c(3,6), aux310c(3,10), aux315c(3,15)

  !$ t1 = omp_get_wtime()
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3) &
  !$OMP SHARED(order,ffts)
  !$OMP DO
  do ib3 = 1,ffts%nbox3
    do ib2 = 1,ffts%nbox2
      do ib1 = 1,ffts%nbox1
        ffts%fourr3(ib3,ib2,ib1,:) = 0
        if (order > 0) ffts%fourr33(ib3,ib2,ib1,:,:) = 0
        if (order > 1) ffts%fourr36(ib3,ib2,ib1,:,:) = 0
        if (order > 2) ffts%fourr310(ib3,ib2,ib1,:,:) = 0
        if (order > 3) ffts%fourr315(ib3,ib2,ib1,:,:) = 0
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    nulify fourgrid:',t2 - t1

  !$ t1 = omp_get_wtime()
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ib1,ib2,ib3,i,j,k) &
  !$OMP SHARED(nbox,boxes,ffts,nbox1,nbox2,nbox3,order)
  !$OMP DO
  do ib = 1,nbox
    ib1 = boxes(ib)%ib_grid(1) + nbox1 - 1
    ib2 = boxes(ib)%ib_grid(2) + nbox2 - 1
    ib3 = boxes(ib)%ib_grid(3) + nbox3 - 1
    
    ! write(*,*)ib, boxes(ib)%ib_grid(:), ib1, ib2, ib3
    do i = 1,3
      ffts%fourr3(ib3,ib2,ib1,i) = boxes(ib)%outgoing0f_1d(i)
      if (order > 0) then
        do j = 1,3
          ffts%fourr33(ib3,ib2,ib1,i,j) = boxes(ib)%outgoing1f_1d(i,j)
        enddo
      endif
      if (order > 1) then
        do j = 1,6
          ffts%fourr36(ib3,ib2,ib1,i,j) = boxes(ib)%outgoing2f_1d(i,j)
        enddo
      endif
      if (order > 2) then
        do j = 1,10
          ffts%fourr310(ib3,ib2,ib1,i,j) = boxes(ib)%outgoing3f_1d(i,j)
        enddo
      endif
      if (order > 3) then
        do j = 1,15
          ffts%fourr315(ib3,ib2,ib1,i,j) = boxes(ib)%outgoing4f_1d(i,j)
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    outgoing to fourgrid:',t2 - t1

  !$ t1 = omp_get_wtime()
  call fftw_execute_dft_r2c(ffts%plan3, ffts%fourr3, ffts%fourc3)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    FFT3:',t2 - t1
  !$ if (order > 0) t1 = omp_get_wtime()
  if (order > 0) call fftw_execute_dft_r2c(ffts%plan33, ffts%fourr33, ffts%fourc33)
  !$ if (order > 0) t2 = omp_get_wtime()
  !$ if (order > 0) write(*,*) '    FFT33:',t2 - t1
  !$ if (order > 1) t1 = omp_get_wtime()
  if (order > 1) call fftw_execute_dft_r2c(ffts%plan36, ffts%fourr36, ffts%fourc36)
  !$ if (order > 1) t2 = omp_get_wtime()
  !$ if (order > 1) write(*,*) '    FFT36:',t2 - t1
  !$ if (order > 2) t1 = omp_get_wtime()
  if (order > 2) call fftw_execute_dft_r2c(ffts%plan310, ffts%fourr310, ffts%fourc310)
  !$ if (order > 2) t2 = omp_get_wtime()
  !$ if (order > 2) write(*,*) '    FFT310:',t2 - t1
  !$ if (order > 3) t2 = omp_get_wtime()
  if (order > 3) call fftw_execute_dft_r2c(ffts%plan315, ffts%fourr315, ffts%fourc315)
  !$ if (order > 3) t2 = omp_get_wtime()
  !$ if (order > 3) write(*,*) '    FFT315:',t2 - t1

!   ! test fft
!   call ifft_tensor36box(iplan_advanced36box, fourgrid36box)
!   dum = 0.0
!   do ib = 1,nbox
!     ib1 = boxes(ib)%ib_grid(1) + nbox1 - 1
!     ib2 = boxes(ib)%ib_grid(2) + nbox2 - 1
!     ib3 = boxes(ib)%ib_grid(3) + nbox3 - 1
! 
!     do i = 1,3
!       do j = 1,6
!         dum = dum + abs(real(fourgrid36box(ib3,ib2,ib1,i,j)) - boxes(ib)%outgoing2f_1d(i,j))
!       enddo
!     enddo
!   enddo
!   write(*,*) dum
!   stop

  !$ t1 = omp_get_wtime()
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ib1,ib2,ib3,outgoing0,outgoing1,outgoing2,outgoing3,aux3c,aux33c,aux36c,aux310c) &
  !$OMP PRIVATE(aux315c,outgoing4) &
  !$OMP SHARED(ffts,G_box_hat,dG_dX_box_hat,d2G_dX2_box_hat,d3G_dX3_box_hat,d4G_dX4_box_hat,order) 
  !$OMP DO
  do ib3 = 1,ffts%nbox3cmplx
    do ib2 = 1,ffts%nbox2
      do ib1 = 1,ffts%nbox1

        if (order == 0) then

          outgoing0 = ffts%fourc3(ib3,ib2,ib1,:)

          aux3c = matmul(G_box_hat(:,:,ib1,ib2,ib3),outgoing0)
          ffts%fourc3(ib3,ib2,ib1,:) = aux3c

        elseif (order == 1) then

          outgoing0 = ffts%fourc3(ib3,ib2,ib1,:)
          outgoing1 = ffts%fourc33(ib3,ib2,ib1,:,:)

          aux3c = matmul(G_box_hat(:,:,ib1,ib2,ib3),outgoing0) - &
            TijkTjk_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing1)
          ffts%fourc3(ib3,ib2,ib1,:) = aux3c

          aux33c = TijkTj_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing0)
          ffts%fourc33(ib3,ib2,ib1,:,:) = aux33c

        elseif (order == 2) then

          outgoing0 = ffts%fourc3(ib3,ib2,ib1,:)
          outgoing1 = ffts%fourc33(ib3,ib2,ib1,:,:)
          outgoing2 = ffts%fourc36(ib3,ib2,ib1,:,:)

          aux3c = matmul(G_box_hat(:,:,ib1,ib2,ib3),outgoing0) - &
            TijkTjk_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.5*TijklTjkl_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing2)
          ffts%fourc3(ib3,ib2,ib1,:) = aux3c

          aux33c = TijkTj_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing0) - &
          TijklTjk_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing1)
          ffts%fourc33(ib3,ib2,ib1,:,:) = aux33c

          aux36c = 0.5*TijklTj_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing0)
          ffts%fourc36(ib3,ib2,ib1,:,:) = aux36c

        elseif (order == 3) then

          outgoing0 = ffts%fourc3(ib3,ib2,ib1,:)
          outgoing1 = ffts%fourc33(ib3,ib2,ib1,:,:)
          outgoing2 = ffts%fourc36(ib3,ib2,ib1,:,:)
          outgoing3 = ffts%fourc310(ib3,ib2,ib1,:,:)
  
          aux3c = matmul(G_box_hat(:,:,ib1,ib2,ib3),outgoing0) - &
            TijkTjk_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.5*TijklTjkl_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing2) - &
            0.16666666666666*TijklmTjklm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing3)
          ffts%fourc3(ib3,ib2,ib1,:) = aux3c

          aux33c = TijkTj_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing0) - &
            TijklTjk_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.5*TijklmTjlm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing2)
          ffts%fourc33(ib3,ib2,ib1,:,:) = aux33c

          aux36c = 0.5*TijklTj_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing0) - &
            0.5*TijklmTjm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing1)
          ffts%fourc36(ib3,ib2,ib1,:,:) = aux36c

          aux310c = & 
            0.16666666666666*TijklmTj_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing0)
          ffts%fourc310(ib3,ib2,ib1,:,:) = aux310c

        elseif (order == 4) then

          outgoing0 = ffts%fourc3(ib3,ib2,ib1,:)
          outgoing1 = ffts%fourc33(ib3,ib2,ib1,:,:)
          outgoing2 = ffts%fourc36(ib3,ib2,ib1,:,:)
          outgoing3 = ffts%fourc310(ib3,ib2,ib1,:,:)
          outgoing4 = ffts%fourc315(ib3,ib2,ib1,:,:)
  
          aux3c = matmul(G_box_hat(:,:,ib1,ib2,ib3),outgoing0) - &
            TijkTjk_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.5*TijklTjkl_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing2) - &
            0.16666666666666*TijklmTjklm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing3) + &
            1.0/24.0*TijklmnTjklmn_1d_cmplx(d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3),outgoing4)
          ffts%fourc3(ib3,ib2,ib1,:) = aux3c

          aux33c = TijkTj_1d_cmplx(dG_dX_box_hat(:,:,:,ib1,ib2,ib3),outgoing0) - &
            TijklTjk_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.5*TijklmTjlm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing2) - &
            0.16666666666666*TijklmnTjlmn_1d_cmplx(d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3),outgoing3)
          ffts%fourc33(ib3,ib2,ib1,:,:) = aux33c

          aux36c = 0.5*TijklTj_1d_cmplx(d2G_dX2_box_hat(:,:,:,:,ib1,ib2,ib3),outgoing0) - &
            0.5*TijklmTjm_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing1) + &
            0.25*TijklmnTjmn_1d_cmplx(d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3),outgoing2)
          ffts%fourc36(ib3,ib2,ib1,:,:) = aux36c

          aux310c = & 
            0.16666666666666*TijklmTj_1d_cmplx(d3G_dX3_box_hat(:,:,:,:,:,ib1,ib2,ib3),outgoing0) - &
            0.16666666666666*TijklmnTjn_1d_cmplx(d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3),outgoing1)
          ffts%fourc310(ib3,ib2,ib1,:,:) = aux310c

          aux315c = & 
            1.0/24.0*TijklmnTj_1d_cmplx(d4G_dX4_box_hat(:,:,:,:,:,:,ib1,ib2,ib3),outgoing0)
          ffts%fourc315(ib3,ib2,ib1,:,:) = aux315c

        endif

      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    convolution:',t2 - t1

  !$ t1 = omp_get_wtime()
  call fftw_execute_dft_c2r(ffts%iplan3, ffts%fourc3, ffts%fourr3)
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    IFFT3:',t2 - t1
  !$ if (order > 0) t1 = omp_get_wtime()
  if (order > 0) call fftw_execute_dft_c2r(ffts%iplan33, ffts%fourc33, ffts%fourr33)
  !$ if (order > 0) t2 = omp_get_wtime()
  !$ if (order > 0) write(*,*) '    IFFT33:',t2 - t1
  !$ if (order > 1) t1 = omp_get_wtime()
  if (order > 1) call fftw_execute_dft_c2r(ffts%iplan36, ffts%fourc36, ffts%fourr36)
    !$ if (order > 1) t2 = omp_get_wtime()
  !$ if (order > 1) write(*,*) '    IFFT36:',t2 - t1
  !$ if (order > 2) t1 = omp_get_wtime()
  if (order > 2) call fftw_execute_dft_c2r(ffts%iplan310, ffts%fourc310, ffts%fourr310)
  !$ if (order > 2) t2 = omp_get_wtime()
  !$ if (order > 2) write(*,*) '    IFFT310:',t2 - t1
  !$ if (order > 3) t2 = omp_get_wtime()
  if (order > 3) call fftw_execute_dft_c2r(ffts%iplan315, ffts%fourc315, ffts%fourr315)
  !$ if (order > 3) t2 = omp_get_wtime()
  !$ if (order > 3) write(*,*) '    IFFT315:',t2 - t1

  ! write(*,*) real(fourgrid33box(1,1,1,:,:))
  ! do ib3 = 1,nbox3*2 - 1
  !   do ib2 = 1,nbox2*2 - 1
  !     do ib1 = 1,nbox1*2 - 1
  !       if(abs(real(fourgrid33box(ib3,ib2,ib1,1,1)) - 3.0598907132652592E-005)<1.0e-12)then
  !         write(*,*) ib1,ib2,ib3,fourgrid3box(ib3,ib2,ib1,:)
  !         stop
  !       endif
  !     enddo
  !   enddo
  ! enddo
  ! stop


  !$ t1 = omp_get_wtime()
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib,ib1,ib2,ib3,i,j,k) &
  !$OMP SHARED(nbox,boxes,ffts,nbox1,nbox2,nbox3,order)
  !$OMP DO
  do ib = 1,nbox
    ib1 = boxes(ib)%ib_grid(1) + nbox1 - 1
    ib2 = boxes(ib)%ib_grid(2) + nbox2 - 1
    ib3 = boxes(ib)%ib_grid(3) + nbox3 - 1
    ! write(*,*) ib
    ! write(*,*) boxes(ib)%incoming0(:,:)
    do i = 1,3
      boxes(ib)%incoming0f_1d(i) = ffts%fourr3(ib3,ib2,ib1,i)*ffts%normfact
      if (order > 0) then
        do j = 1,3
          boxes(ib)%incoming1f_1d(i,j) = ffts%fourr33(ib3,ib2,ib1,i,j)*ffts%normfact
        enddo
      endif
      if (order > 1) then
        do j = 1,6
          boxes(ib)%incoming2f_1d(i,j) = ffts%fourr36(ib3,ib2,ib1,i,j)*ffts%normfact
        enddo
      endif
      if (order > 2) then
        do j = 1,10
          boxes(ib)%incoming3f_1d(i,j) = ffts%fourr310(ib3,ib2,ib1,i,j)*ffts%normfact
        enddo
      endif
      if (order > 3) then
        do j = 1,15
          boxes(ib)%incoming4f_1d(i,j) = ffts%fourr315(ib3,ib2,ib1,i,j)*ffts%normfact
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !$ t2 = omp_get_wtime()
  !$ write(*,*) '    fourgrid to incoming:',t2 - t1

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

function get_xi(x) result(xi)
  use tensor_functions

  double precision, intent(in) :: x(3)
  double precision :: xi(3)

  xi = (x - float(floor(x)))*2.0 - 1.0

end

function get_grid_ind(x,npts1,npts2,npts3) result(ind8)
  use tensor_functions

  double precision, intent(in) :: x(3)
  integer, intent(in) :: npts1,npts2,npts3
  integer :: ind8(8,3), ind1(3), i

  ! coordinates
  ind1 = floor(x)
  ind8(1,:) = ind1           ! 1
  ind8(2,:) = ind1 + [1,0,0] ! 2
  ind8(3,:) = ind1 + [1,1,0] ! 3
  ind8(4,:) = ind1 + [0,1,0] ! 4
  ind8(5,:) = ind1 + [0,0,1] ! 5
  ind8(6,:) = ind1 + [1,0,1] ! 6
  ind8(7,:) = ind1 + [1,1,1] ! 7
  ind8(8,:) = ind1 + [0,1,1] ! 8

  ind8 = ind8 + 1 ! shift for 1 for indices (x=0 corresponds i,j,k=1)

  do i = 1,8
    ind8(i,1) = enforce_periodic(ind8(i,1), 1, npts1)
    ind8(i,2) = enforce_periodic(ind8(i,2), 1, npts2)
    ind8(i,3) = enforce_periodic(ind8(i,3), 1, npts3)
  enddo

end


function shape_funct8(xi) result(phi)
  double precision, intent(in) :: xi(3)
  double precision :: phi(8), g, h, r

  g = xi(1)
  h = xi(2)
  r = xi(3)
  phi(1) = 0.125*(1.0 - g)*(1.0 - h)*(1.0 - r) ! abaqus
  phi(2) = 0.125*(1.0 + g)*(1.0 - h)*(1.0 - r)
  phi(3) = 0.125*(1.0 + g)*(1.0 + h)*(1.0 - r)
  phi(4) = 0.125*(1.0 - g)*(1.0 + h)*(1.0 - r)
  phi(5) = 0.125*(1.0 - g)*(1.0 - h)*(1.0 + r)
  phi(6) = 0.125*(1.0 + g)*(1.0 - h)*(1.0 + r)
  phi(7) = 0.125*(1.0 + g)*(1.0 + h)*(1.0 + r)
  phi(8) = 0.125*(1.0 - g)*(1.0 + h)*(1.0 + r)

end

function triangle_area(x) result(a)
  use tensor_functions

  double precision, intent (in) :: x(3,3)
  double precision :: v1(3), v2(3), a

  v1 = - x(:,1) + x(:,2)
  v2 = - x(:,1) + x(:,3)
  a = 0.5*norm2(cross_prod(v1,v2))

end

subroutine create_faces(elements, faces)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (out) :: faces(:)
  type (face_type), allocatable :: all_faces(:)

  integer :: in1, in2, in3, iel, ic, i, j, ifc, ind1, ind2, ifc1, ifc2, ind(3)
  integer :: nel_per_boxind = 10, n_boxind, nb1, nb2, nb3, ib1, ib2, ib3, nind_boxind
  integer, allocatable :: unique_faces_ind(:,:),ind_duplicate(:), boxind(:,:,:,:), ind_duplicate_box(:,:,:,:)
  integer, allocatable :: nboxind(:,:,:), unique_faces_ind_box(:,:,:,:,:), unique_faces_to_1d(:,:,:,:)
  double precision :: tol = 1.0e-5, x(3), dx1(3), dx2(3), dx(3), dir_range(3), x1(3), x2(3), x3(3)
  double precision :: xnode_face(3,3), boxind_edge_length, dxboxind(3)
  logical :: new

  allocate(all_faces(nel*4))
  allocate(unique_faces_ind(2,nel*4))

  el_face_nodes(:,1) = [1,2,3]
  el_face_nodes(:,2) = [1,2,4]
  el_face_nodes(:,3) = [1,3,4]
  el_face_nodes(:,4) = [2,3,4]

  dir_range(1) = float(npts1)
  dir_range(2) = float(npts2)
  dir_range(3) = float(npts3)

  n_boxind = nint(float(nel)/float(nel_per_boxind))
  boxind_edge_length = (dir_range(1)*dir_range(2)*dir_range(3)/float(n_boxind))**(0.333333)
  nb1 = nint(dir_range(1)/boxind_edge_length)
  nb2 = nint(dir_range(2)/boxind_edge_length)
  nb3 = nint(dir_range(3)/boxind_edge_length)
  dxboxind(1) = dir_range(1)/float(nb1)
  dxboxind(2) = dir_range(2)/float(nb2)
  dxboxind(3) = dir_range(3)/float(nb3)
  nind_boxind = nel_per_boxind*4*100
  allocate(boxind(nind_boxind,nb1,nb2,nb3))
  allocate(nboxind(nb1,nb2,nb3))
  allocate(ind_duplicate_box(nind_boxind,nb1,nb2,nb3))
  allocate(unique_faces_ind_box(2,nind_boxind,nb1,nb2,nb3))
  allocate(unique_faces_to_1d(2,nb1,nb2,nb3))

  write(*,*) '   all faces'
  ifc = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,ifc,x1,x2,x3,i) &
  !$OMP SHARED(nel,elements,el_face_nodes,all_faces,tol,dir_range)
  !$OMP DO 
  do iel = 1,nel
    ifc = (iel-1)*4 + 1
    x1 = elements(iel)%xnode(:,el_face_nodes(1,1))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,1))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,1))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 1
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 2
    x1 = elements(iel)%xnode(:,el_face_nodes(1,2))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,2))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,2))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 2
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 3
    x1 = elements(iel)%xnode(:,el_face_nodes(1,3))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,3))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,3))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 3
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 4
    x1 = elements(iel)%xnode(:,el_face_nodes(1,4))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,4))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,4))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 4
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,ic,ifc,ind) &
  !$OMP SHARED(nb1,nb2,nb3,nel,all_faces,boxind,nboxind,nind_boxind,dxboxind)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1

    ic = 0
    do ifc = 1,nel*4
      ind = ceiling(all_faces(ifc)%x/dxboxind)
      if (ind(1) == 0) ind(1) = 1
      if (ind(1) == nb1 + 1) ind(1) = nb1
      if (ind(2) == 0) ind(2) = 1
      if (ind(2) == nb2 + 1) ind(2) = nb2
      if (ind(3) == 0) ind(3) = 1
      if (ind(3) == nb3 + 1) ind(3) = nb3
      if (ind(1) == ib1 .and. ind(2) == ib2 .and. ind(3) == ib3) then
        ic = ic + 1
        if (ic > nind_boxind) then
          write(*,*) ic, nind_boxind
          stop 'ic > nind_boxind'
        endif
        boxind(ic,ib1,ib2,ib3) = ifc
      endif
    enddo
    nboxind(ib1,ib2,ib3) = ic

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  if (sum(nboxind) .ne. 4*nel) stop 'sum(nboxind) .ne. 4*nel'
  write(*,*) sum(nboxind), 4*nel

!  ic = 1
!  unique_faces_ind = 0
!  unique_faces_ind(1,ic) = 1
!  do ifc = 2,nel*4
!
!    x = all_faces(ifc)%x
!    new = .true.
!    if (.not.all_faces(ifc)%boundary) then ! if this is boundary face treat it as new
!      do i = 1,ic
!        if (unique_faces_ind(2,i) == 0) then ! check if it is the same only for faces that do not have duplicate assigned
!          ind1 = unique_faces_ind(1,i)
!          ! if (norm2(x - all_faces(ind1)%x) <= tol) then 
!          dx = abs(x - all_faces(ind1)%x)
!          if (.not.(all_faces(ind1)%boundary).and. &
!              (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol)) then ! duplicate, skip comparing with stored boundary faces
!            unique_faces_ind(2,i) = ifc
!            new = .false.
!            exit
!          endif
!        endif
!      enddo
!    endif
!
!    if (new) then
!
!      ic = ic + 1
!      unique_faces_ind(1,ic) = ifc
!
!    endif
!
!  enddo

  write(*,*) '   duplicates'
  ind_duplicate_box = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,i,ifc1,x,j,ifc2,dx) &
  !$OMP SHARED(nb1,nb2,nb3,nboxind,boxind,all_faces,tol,ind_duplicate_box)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1

    do i = 1, nboxind(ib1,ib2,ib3)
      ifc1 = boxind(i,ib1,ib2,ib3)
      x = all_faces(ifc1)%x
      if (.not.all_faces(ifc1)%boundary) then
        do j = 1, nboxind(ib1,ib2,ib3)
          ifc2 = boxind(j,ib1,ib2,ib3)
          dx = abs(x - all_faces(ifc2)%x)
          if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol .and. ifc1 .ne. ifc2) then ! duplicate, skip comparing with stored boundary faces
            ind_duplicate_box(i,ib1,ib2,ib3) = j
            exit
          endif
        enddo
      endif
    enddo

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   unique faces'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,ic,i,ifc,j) &
  !$OMP SHARED(nb1,nb2,nb3,boxind,nboxind,ind_duplicate_box,unique_faces_ind_box)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    ic = 0
    do i = 1, nboxind(ib1,ib2,ib3)
      ifc = boxind(i,ib1,ib2,ib3)
      if (ind_duplicate_box(i,ib1,ib2,ib3) == 0) then
        ic = ic + 1
        unique_faces_ind_box(1,ic,ib1,ib2,ib3) = ifc
        unique_faces_ind_box(2,ic,ib1,ib2,ib3) = 0
      else
        j = ind_duplicate_box(i,ib1,ib2,ib3)
        if (j > i) then
          ic = ic + 1
          unique_faces_ind_box(1,ic,ib1,ib2,ib3) = ifc
          unique_faces_ind_box(2,ic,ib1,ib2,ib3) = boxind(j,ib1,ib2,ib3)
        endif
      endif
    enddo
    nboxind(ib1,ib2,ib3) = ic ! overwrite

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  i = 1
  j = 1
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    j = j + nboxind(ib1,ib2,ib3) - 1
    unique_faces_to_1d(:,ib1,ib2,ib3) = [i, j]
    j = j + 1
    i = j

    ! write(*,*) unique_faces_to_1d(:,ib1,ib2,ib3), nboxind(ib1,ib2,ib3)

  enddo
  enddo
  enddo

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,i,j) &
  !$OMP SHARED(nb1,nb2,nb3,unique_faces_to_1d,unique_faces_ind,unique_faces_ind_box,nboxind)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    i = unique_faces_to_1d(1,ib1,ib2,ib3)
    j = unique_faces_to_1d(2,ib1,ib2,ib3)
    unique_faces_ind(:,i:j) = unique_faces_ind_box(:,1:nboxind(ib1,ib2,ib3),ib1,ib2,ib3)

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  ic = unique_faces_to_1d(2,nb1,nb2,nb3)
  write(*,*) ic

!   stop

!   write(*,*) '   duplicates'
!   allocate(ind_duplicate(nel*4))
!   ind_duplicate = 0
!   !$OMP PARALLEL DEFAULT(NONE) &
!   !$OMP PRIVATE(ifc1,x,ifc2,dx) &
!   !$OMP SHARED(nel,all_faces,ind_duplicate,tol)
!   !$OMP DO 
!   do ifc1 = 1,nel*4
!     x = all_faces(ifc1)%x
!     if (.not.all_faces(ifc1)%boundary) then ! if this is boundary face treat it as new
!       do ifc2 = 1,nel*4
!         dx = abs(x - all_faces(ifc2)%x)
!         if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol .and. ifc1 .ne. ifc2) then ! duplicate, skip comparing with stored boundary faces
!           ind_duplicate(ifc1) = ifc2
!           exit
!         endif
!       enddo
!     endif
!   enddo
!   !$OMP END DO
!   !$OMP END PARALLEL
! 
!   write(*,*) '   unique faces'
!   ic = 0
!   unique_faces_ind = 0
!   do ifc = 1,nel*4
!     if (mod(ifc,10000) == 0) write(*,*) float(ifc)/float(nel*4)
!     if (ind_duplicate(ifc) == 0) then
!       ic = ic + 1
!       unique_faces_ind(1,ic) = ifc
!     else
!       if (ind_duplicate(ifc) > ifc) then
!         ic = ic + 1
!         unique_faces_ind(1,ic) = ifc
!         unique_faces_ind(2,ic) = ind_duplicate(ifc)
!       endif
!     endif
!   enddo
! 
!   deallocate(ind_duplicate)

  write(*,*) '   create faces'
  allocate(faces(ic))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,ind1,iel,x1,x2,x3,xnode_face,i,ind2) &
  !$OMP SHARED(ic,faces,all_faces,elements,unique_faces_ind)
  !$OMP DO 
  do ifc = 1,ic
    ind1 = unique_faces_ind(1,ifc)
    faces(ifc)%x = all_faces(ind1)%x
    faces(ifc)%el(1) = all_faces(ind1)%el(1)
    faces(ifc)%el_face(1) = all_faces(ind1)%el_face(1)
    faces(ifc)%boundary = all_faces(ind1)%boundary
    faces(ifc)%normal = all_faces(ind1)%normal

    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
    xnode_face(:,1) = x1
    xnode_face(:,2) = x2
    xnode_face(:,3) = x3
    faces(ifc)%area = triangle_area(xnode_face)

    iel = faces(ifc)%el(1)
    i = faces(ifc)%el_face(1)
    if (elements(iel)%face(i) .ne. 0) stop 'face already assignd to this element place'
    elements(iel)%face(i) = ifc
    elements(iel)%face_normal_sgn(i) = 1.0 ! first is positive

    ind2 = unique_faces_ind(2,ifc)
    if (ind2 .ne. 0) then
      faces(ifc)%el(2) = all_faces(ind2)%el(1)
      faces(ifc)%el_face(2) = all_faces(ind2)%el_face(1)
      faces(ifc)%nel = 2

      iel = faces(ifc)%el(2)
      i = faces(ifc)%el_face(2)
      if (elements(iel)%face(i) .ne. 0) stop 'face already assignd to this element place'
      elements(iel)%face(i) = ifc
      elements(iel)%face_normal_sgn(i) = -1.0 ! second is negative

    else
      faces(ifc)%el(2) = 0
      faces(ifc)%nel = 1
      ! write(*,*) ifc, ind1
    endif

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  nfc = ic
  ! write(*,*)ic, nel

  ! verify
  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,i,iel,x,dx) &
  !$OMP SHARED(ic,faces,tol,elements,el_face_nodes)
  !$OMP DO 
  do ifc = 1,ic
    if (faces(ifc)%nel == 2) then

      i = faces(ifc)%el_face(1)
      iel = faces(ifc)%el(1)
      x = (elements(iel)%xnode(:,el_face_nodes(1,i)) + elements(iel)%xnode(:,el_face_nodes(2,i)) + &
       elements(iel)%xnode(:,el_face_nodes(3,i)))/3.0

      ! call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
      ! if (norm2( (x1 + x2 + x3)/3.0 - x) > tol) then
      !   stop 'problem'
      ! endif

      dx = abs(x - faces(ifc)%x)
      if (.not.((dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol)) ) then ! duplicate
        stop 'error in face'
      endif

    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,i) &
  !$OMP SHARED(nel,elements)
  !$OMP DO 
  do iel = 1,nel
    do i = 1,4
      if (elements(iel)%face(i) == 0) stop 'element without all faces'
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(all_faces)
  deallocate(unique_faces_ind)

  deallocate(boxind)
  deallocate(nboxind)
  deallocate(ind_duplicate_box)
  deallocate(unique_faces_ind_box)
  deallocate(unique_faces_to_1d)

  ! stop

end

subroutine create_edges(elements, faces, edges)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (edge_type), allocatable, intent (out) :: edges(:)
  type (edge_type), allocatable :: all_edges(:)

  integer :: iel, ifc, ind1, ind2, iedg, ic, i, j, ind, iedg1, iedg2, ind3d(3), k, l
  integer :: nel_per_boxind = 10, n_boxind, nb1, nb2, nb3, ib1, ib2, ib3, nind_boxind
  integer, allocatable :: unique_edges_ind(:,:), nduplicates(:), ind_duplicate(:,:), nduplicates_all(:)
  integer, allocatable :: nboxind(:,:,:), unique_edges_ind_box(:,:,:,:,:), unique_edges_to_1d(:,:,:,:)
  integer, allocatable :: ind_duplicate_box(:,:,:,:,:), boxind(:,:,:,:), nduplicate_box(:,:,:,:), nduplicate_all_box(:,:,:,:)
  double precision :: tol = 1.0e-5, x(3), dx1(3), dx2(3), dx(3),  x1(3), x2(3), x3(3)
  double precision :: xnode_face(3,3), boxind_edge_length,dxboxind(3)
  logical :: new

  allocate(all_edges(nfc*3))
  allocate(unique_edges_ind(nmx_shared_faces_per_edge,nfc*3))
  allocate(nduplicates(nfc*3))

  n_boxind = nint(float(nel)/float(nel_per_boxind))
  boxind_edge_length = (float(npts1*npts2*npts3)/float(n_boxind))**(0.333333)
  nb1 = nint(float(npts1)/boxind_edge_length)
  nb2 = nint(float(npts2)/boxind_edge_length)
  nb3 = nint(float(npts3)/boxind_edge_length)
  dxboxind(1) = float(npts1)/float(nb1)
  dxboxind(2) = float(npts2)/float(nb2)
  dxboxind(3) = float(npts3)/float(nb3)
  nind_boxind = nel_per_boxind*3*100
  allocate(boxind(nind_boxind,nb1,nb2,nb3))
  allocate(nboxind(nb1,nb2,nb3))
  allocate(ind_duplicate_box(nmx_shared_faces_per_edge,nind_boxind,nb1,nb2,nb3))
  allocate(unique_edges_ind_box(nmx_shared_faces_per_edge,nind_boxind,nb1,nb2,nb3))
  allocate(unique_edges_to_1d(nmx_shared_faces_per_edge,nb1,nb2,nb3))
  allocate(nduplicate_all_box(nind_boxind,nb1,nb2,nb3))
  allocate(nduplicate_box(nind_boxind,nb1,nb2,nb3))

  face_edge_nodes(:,1) = [1,2]
  face_edge_nodes(:,2) = [2,3]
  face_edge_nodes(:,3) = [3,1]

  face_edge_nodes_remaining(1) = 3
  face_edge_nodes_remaining(2) = 1
  face_edge_nodes_remaining(3) = 2

  write(*,*) '   all edges'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,iedg,iel,xnode_face,x1,x2) &
  !$OMP SHARED(nfc,faces,elements,face_edge_nodes,all_edges)
  !$OMP DO 
  do ifc = 1,nfc

    iedg = (ifc-1)*3 + 1
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,1))
    x2 = xnode_face(:,face_edge_nodes(2,1))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 1

    iedg = (ifc-1)*3 + 2
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,2))
    x2 = xnode_face(:,face_edge_nodes(2,2))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 2

    iedg = (ifc-1)*3 + 3
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,3))
    x2 = xnode_face(:,face_edge_nodes(2,3))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 3

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,ic,iedg,ind3d) &
  !$OMP SHARED(nb1,nb2,nb3,nfc,all_edges,dxboxind,nind_boxind,boxind,nboxind)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1

    ic = 0
    do iedg = 1,nfc*3
      ind3d = ceiling(all_edges(iedg)%xc/dxboxind)
      if (ind3d(1) == 0) ind3d(1) = 1
      if (ind3d(1) == nb1 + 1) ind3d(1) = nb1
      if (ind3d(2) == 0) ind3d(2) = 1
      if (ind3d(2) == nb2 + 1) ind3d(2) = nb2
      if (ind3d(3) == 0) ind3d(3) = 1
      if (ind3d(3) == nb3 + 1) ind3d(3) = nb3
      if (ind3d(1) == ib1 .and. ind3d(2) == ib2 .and. ind3d(3) == ib3) then
        ic = ic + 1
        if (ic > nind_boxind) then
          write(*,*) ic, nind_boxind
          stop 'ic > nind_boxind'
        endif
        boxind(ic,ib1,ib2,ib3) = iedg
      endif
    enddo
    nboxind(ib1,ib2,ib3) = ic

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  ! write(*,*) sum(nboxind), 3*nfc


!  ic = 1
!  unique_edges_ind = 0
!  nduplicates = 0
!  unique_edges_ind(1,ic) = 1
!  nduplicates(ic) = 1
!  do iedg = 2,nfc*3
!
!    x = all_edges(iedg)%xc
!    new = .true.
!    do i = 1,ic
!      do j = 1,nduplicates(i)
!        ind1 = unique_edges_ind(j,i)
!        dx = abs(x - all_edges(ind1)%xc)
!        if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol) then ! duplicate
!          nduplicates(i) = nduplicates(i) + 1
!          if (nduplicates(i) > nmx_shared_faces_per_edge) stop 'nduplicates(i) > nmx_shared_faces_per_edge'
!          unique_edges_ind(nduplicates(i),i) = iedg
!          new = .false.
!          exit
!        endif
!      enddo
!    enddo
!
!    if (new) then
!
!      ic = ic + 1
!      unique_edges_ind(1,ic) = iedg
!      nduplicates(ic) = 1
!
!    endif
!
!  enddo

  write(*,*) '   duplicates'
  ind_duplicate_box = 0
  nduplicate_all_box = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,i,iedg1,x,ic,j,iedg2,dx) &
  !$OMP SHARED(nb1,nb2,nb3,nboxind,boxind,all_edges,tol,ind_duplicate_box,nduplicate_all_box)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1

    do i = 1, nboxind(ib1,ib2,ib3)
      iedg1 = boxind(i,ib1,ib2,ib3)
      x = all_edges(iedg1)%xc

      ic = 0
      do j = 1, nboxind(ib1,ib2,ib3)
        iedg2 = boxind(j,ib1,ib2,ib3)
        dx = abs(x - all_edges(iedg2)%xc)
        if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol) then 
          ic = ic + 1
          ind_duplicate_box(ic,i,ib1,ib2,ib3) = j
          nduplicate_all_box(i,ib1,ib2,ib3) = nduplicate_all_box(i,ib1,ib2,ib3) + 1
        endif
      enddo

    enddo

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   unique edges'
  unique_edges_ind_box = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,ic,i,j,k,l) &
  !$OMP SHARED(nb1,nb2,nb3,nboxind,ind_duplicate_box,nduplicate_all_box) &
  !$OMP SHARED(unique_edges_ind_box,boxind,nduplicate_box)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    ic = 0
    do i = 1, nboxind(ib1,ib2,ib3)

      j = minval(ind_duplicate_box(:,i,ib1,ib2,ib3), mask = ind_duplicate_box(:,i,ib1,ib2,ib3) > 0)
      if (i == j) then
        ic = ic + 1
        do k = 1, nduplicate_all_box(i,ib1,ib2,ib3)
          l = ind_duplicate_box(k,i,ib1,ib2,ib3)
          unique_edges_ind_box(k,ic,ib1,ib2,ib3) = boxind(l,ib1,ib2,ib3)
        enddo
        nduplicate_box(ic,ib1,ib2,ib3) = nduplicate_all_box(i,ib1,ib2,ib3)
      endif
      
    enddo
    nboxind(ib1,ib2,ib3) = ic ! overwrite

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  i = 1
  j = 1
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    j = j + nboxind(ib1,ib2,ib3) - 1
    unique_edges_to_1d(:,ib1,ib2,ib3) = [i, j]
    j = j + 1
    i = j

    ! write(*,*) unique_edges_to_1d(:,ib1,ib2,ib3), nboxind(ib1,ib2,ib3)

  enddo
  enddo
  enddo

  unique_edges_ind = 0
  nduplicates = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1,ib2,ib3,i,j) &
  !$OMP SHARED(nb1,nb2,nb3,unique_edges_to_1d,unique_edges_ind) &
  !$OMP SHARED(unique_edges_ind_box,nduplicates,nduplicate_box,nboxind)
  !$OMP DO COLLAPSE(3)
  do ib3 = 1,nb3
  do ib2 = 1,nb2
  do ib1 = 1,nb1
    
    i = unique_edges_to_1d(1,ib1,ib2,ib3)
    j = unique_edges_to_1d(2,ib1,ib2,ib3)
    unique_edges_ind(:,i:j) = unique_edges_ind_box(:,1:nboxind(ib1,ib2,ib3),ib1,ib2,ib3)
    nduplicates(i:j) = nduplicate_box(1:nboxind(ib1,ib2,ib3),ib1,ib2,ib3)

  enddo
  enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  ic = unique_edges_to_1d(2,nb1,nb2,nb3)
  write(*,*) ic

!   write(*,*) '   duplicates'
!   allocate(ind_duplicate(nmx_shared_faces_per_edge,nfc*3))
!   allocate(nduplicates_all(nfc*3))
!   ind_duplicate = 0
!   nduplicates_all = 0
!   write(*,*)nfc*3
!   !$OMP PARALLEL DEFAULT(NONE) &
!   !$OMP PRIVATE(iedg1,x,ic,iedg2,dx) &
!   !$OMP SHARED(all_edges,nfc,ind_duplicate,nduplicates_all,tol)
!   !$OMP DO 
!   do iedg1 = 1,nfc*3
!     x = all_edges(iedg1)%xc
!     ic = 0
!     do iedg2 = 1,nfc*3
!       dx = abs(x - all_edges(iedg2)%xc)
!       if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol) then ! duplicate
!         ic = ic + 1
!         ind_duplicate(ic,iedg1) = iedg2
!         nduplicates_all(iedg1) = nduplicates_all(iedg1) + 1
!       endif
!     enddo
!     ! write(*,*)iedg1, ind_duplicate(:,iedg1)
!   enddo
!   !$OMP END DO
!   !$OMP END PARALLEL
! 
!   write(*,*) '   unique edges'
!   unique_edges_ind = 0
!   nduplicates = 0
!   ic = 0
!   do iedg1 = 1,nfc*3
! 
!     if (mod(iedg1,10000) == 0) write(*,*) float(iedg1)/float(nfc*3)
! 
!     i = minval(ind_duplicate(:,iedg1), mask = ind_duplicate(:,iedg1) > 0)
!     if (iedg1 == i) then
!       ic = ic + 1
!       unique_edges_ind(:,ic) = ind_duplicate(:,iedg1)
!       nduplicates(ic) = nduplicates_all(iedg1)
!     endif
! 
!   enddo
! 
!   deallocate(ind_duplicate)
!   deallocate(nduplicates_all)
!   ! stop

  ! do iedg = 1,ic
  !   write(*,*) iedg, nduplicates(iedg), unique_edges_ind(1:nduplicates(iedg),iedg)
  !   do i = 1,nduplicates(iedg)
  !     write(*,*) all_edges(unique_edges_ind(i,iedg))%xc
  !   enddo
  ! enddo

  write(*,*) '   create edges'
  allocate(edges(ic))
  allocate(u_edge(3,ic))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,ind1,i,ind,ifc) &
  !$OMP SHARED(ic,unique_edges_ind,edges,all_edges,nduplicates,faces)
  !$OMP DO 
  do iedg = 1,ic
    ind1 = unique_edges_ind(1,iedg)
    
    edges(iedg)%x = all_edges(ind1)%x
    edges(iedg)%xc = all_edges(ind1)%xc
    edges(iedg)%nface = nduplicates(iedg)

    allocate(edges(iedg)%face(edges(iedg)%nface))
    allocate(edges(iedg)%face_edge(edges(iedg)%nface))
    allocate(edges(iedg)%face_G(3,3,edges(iedg)%nface))

    do i = 1,nduplicates(iedg)
      ind = unique_edges_ind(i,iedg)
      edges(iedg)%face(i) = all_edges(ind)%face(1)
      edges(iedg)%face_edge(i) = all_edges(ind)%face_edge(1)
    enddo

    do i = 1,edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      if (faces(ifc)%edge(edges(iedg)%face_edge(i)) == 0) then
        faces(ifc)%edge(edges(iedg)%face_edge(i)) = iedg
      else
        stop 'edge already placed in this face'
      endif
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  nedg = ic

  write(*,*) '   close neighboring edges'
  ! close neighboring edges (belong to faces containing an edge)
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,ic,i,ifc,j) &
  !$OMP SHARED(nedg,edges,faces)
  !$OMP DO 
  do iedg = 1,nedg

    allocate(edges(iedg)%neigh_edges(edges(iedg)%nface*2)) ! 2 times number of faces this edge belongs to
    ic = 0
    do i = 1,edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      do j = 1,3
        if (faces(ifc)%edge(j) .ne. iedg) then
          ic = ic + 1
          edges(iedg)%neigh_edges(ic) = faces(ifc)%edge(j)
        endif
      enddo
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

!  iedg = 77
!  write(*,*) edges(iedg)%x(:,1)
!  write(*,*) edges(iedg)%x(:,2)
!  do i = 1,size(edges(iedg)%neigh_edges)
!    j = edges(iedg)%neigh_edges(i)
!    write(*,*) j
!    write(*,*) edges(j)%x(:,1)
!    write(*,*) edges(j)%x(:,2) ! must share one node with other
!  enddo

!  do iedg = 1,nedg
!    
!    write(*,*) iedg
!    do i = 1,edges(iedg)%nface
!      ifc = edges(iedg)%face(i)
!      iel = faces(ifc)%el(1)
!      call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
!      write(*,*) norm2(x1), norm2(x2), norm2(x3)
!    enddo
!
!  enddo

  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,i) &
  !$OMP SHARED(nfc,faces)
  !$OMP DO 
  do ifc = 1,nfc
    do i = 1,3
      if (faces(ifc)%edge(i) == 0) then
        stop 'faces(ifc)%edge(i) == 0'
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(all_edges)
  deallocate(unique_edges_ind)
  deallocate(nduplicates)

  deallocate(boxind)
  deallocate(nboxind)
  deallocate(ind_duplicate_box)
  deallocate(unique_edges_ind_box)
  deallocate(unique_edges_to_1d)
  deallocate(nduplicate_all_box)
  deallocate(nduplicate_box)

end

subroutine get_face_nodes(face, element, x1, x2, x3)
  use global, only : el_face_nodes
  use types
  use tensor_functions

  type (face_type), intent (in) :: face
  type (element_type), intent (in) :: element
  double precision, intent(out) :: x1(3), x2(3), x3(3)
  integer :: i

  i = face%el_face(1)
  x1 = element%xnode(:,el_face_nodes(1,i))
  x2 = element%xnode(:,el_face_nodes(2,i))
  x3 = element%xnode(:,el_face_nodes(3,i))

end

subroutine initialize_integration_points(intg_pt)
  use global, only : nintg_pt
  use types
  type (intg_pt_type), intent (inout) :: intg_pt(nintg_pt)

  intg_pt(1)%Xhat(1) = 0.16385495
  intg_pt(1)%Xhat(2) = 0.04756957
  intg_pt(1)%Xhatn = norm2(intg_pt(1)%Xhat)
  intg_pt(1)%w = 0.31161231

  intg_pt(2)%Xhat(1) = 0.61114353
  intg_pt(2)%Xhat(2) = 0.17753138
  intg_pt(2)%Xhatn = norm2(intg_pt(2)%Xhat)
  intg_pt(2)%w = 0.31161231

  intg_pt(3)%Xhat(1) = 0.04756957
  intg_pt(3)%Xhat(2) = 0.16385495
  intg_pt(3)%Xhatn = norm2(intg_pt(3)%Xhat)
  intg_pt(3)%w = 0.31161231

  intg_pt(4)%Xhat(1) = 0.17753138
  intg_pt(4)%Xhat(2) = 0.61114353
  intg_pt(4)%Xhatn = norm2(intg_pt(4)%Xhat)
  intg_pt(4)%w = 0.31161231

end

subroutine initialize_integration_points_quad(intg_pt_quad)
  use global, only : nintg_pt_quad, nintg_pt_quad_1d
  use types
  type (intg_pt_type), intent (inout) :: intg_pt_quad(nintg_pt_quad)
  double precision :: a(nintg_pt_quad_1d), w(nintg_pt_quad_1d), wtot
  integer :: i, j, ic

  ! order 3
  if (nintg_pt_quad_1d == 3) then
    a(1:3) = [-sqrt(0.6),0.0,sqrt(0.6)]
    w(1:3) = [5.0/9.0,8.0/9.0,5.0/9.0]
    ic = 0
    wtot = 0.0
    do i = 1,nintg_pt_quad_1d
      do j = 1,nintg_pt_quad_1d
        ic = ic + 1
        intg_pt_quad(ic)%Xhat(1) = a(i)
        intg_pt_quad(ic)%Xhat(2) = a(j)
        intg_pt_quad(ic)%Xhatn = norm2(intg_pt_quad(ic)%Xhat)
        intg_pt_quad(ic)%w = w(i)*w(j)
        wtot = wtot + intg_pt_quad(ic)%w 
        ! write(*,*) intg_pt_quad(ic)%Xhat
        ! write(*,*) intg_pt_quad(ic)%w 
      enddo
    enddo
    ! write(*,*) wtot

  ! order 6
  elseif (nintg_pt_quad_1d == 6) then

    a(1:6) = &
        [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, &
          0.2386191860831969,  0.6612093864662645,  0.9324695142031521]
    w(1:6) = &
        [ 0.1713244923791704,  0.3607615730481386,  0.4679139345726910, &
          0.4679139345726910,  0.3607615730481386,  0.1713244923791704]
  
    ic = 0
    wtot = 0.0
    do i = 1,nintg_pt_quad_1d
      do j = 1,nintg_pt_quad_1d
        ic = ic + 1
        intg_pt_quad(ic)%Xhat(1) = a(i)
        intg_pt_quad(ic)%Xhat(2) = a(j)
        intg_pt_quad(ic)%Xhatn = norm2(intg_pt_quad(ic)%Xhat)
        intg_pt_quad(ic)%w = w(i)*w(j)
        wtot = wtot + intg_pt_quad(ic)%w 
        ! write(*,*) intg_pt_quad(ic)%Xhat
        ! write(*,*) intg_pt_quad(ic)%w 
      enddo
    enddo

  endif

end

subroutine intg_pt_to_current(intg_pt, x1, x2, x3, xintpt)
  use global, only : nintg_pt
  use types
  type (intg_pt_type), intent (in) :: intg_pt(nintg_pt)
  double precision, intent (in) :: x1(3), x2(3), x3(3)
  double precision, intent (out) :: xintpt(3,nintg_pt)
  double precision :: h(3)
  integer :: ipt
  
  do ipt = 1,nintg_pt

    h(1) = 1.0 - intg_pt(ipt)%Xhat(1) - intg_pt(ipt)%Xhat(2)
    h(2) = intg_pt(ipt)%Xhat(1)
    h(3) = intg_pt(ipt)%Xhat(2)

    xintpt(:,ipt) = h(1)*x1 + h(2)*x2 + h(3)*x3

  enddo

end

subroutine jacobian(x1, x2, x3, jacob)
  use tensor_functions, only : determinant33, cross_prod
  double precision, intent (in) :: x1(3), x2(3), x3(3)
  double precision, intent (out) :: jacob
  double precision :: dh_dXhat(3,2), J22(2,2), n(3), e1(3), e2(3), e3(3), Qg2l(3,3)
  double precision :: x1loc(3), x2loc(3), x3loc(3)
  integer :: i, j
  
  dh_dXhat(1,:) = [-1.0, -1.0]
  dh_dXhat(2,:) = [ 1.0,  0.0]
  dh_dXhat(3,:) = [ 0.0,  1.0]
  
  n = cross_prod(x2 - x1, x3 - x1)
  n = n/norm2(n)

  e3 = n
  e1 = x2 - x1
  e1 = e1/norm2(e1)
  e2 = cross_prod(e3, e1)

  Qg2l(1,:) = e1
  Qg2l(2,:) = e2
  Qg2l(3,:) = e3

  x1loc = matmul(Qg2l, x1)
  x2loc = matmul(Qg2l, x2)
  x3loc = matmul(Qg2l, x3)

  if (.not. (abs(x1loc(3) - x2loc(3)) < 1.0e-7 .and. abs(x1loc(3) - x3loc(3)) < 1.0e-7) ) then
    write(*,*) abs(x1loc(3) - x2loc(3))
    write(*,*) abs(x1loc(3) - x3loc(3))
    stop 'third coordinat not the same in jacobian'
  endif

  do i = 1,2
    do j = 1,2
      J22(i,j) = x1loc(i)*dh_dXhat(1,j) + x2loc(i)*dh_dXhat(2,j) + x3loc(i)*dh_dXhat(3,j)
    enddo
  enddo
  jacob = abs(J22(1,1)*J22(2,2) - J22(1,2)*J22(2,1))

end

subroutine get_lambda_mu(c0, lambda, mu)
  double precision, intent(in) :: c0(3,3,3,3)
  double precision, intent(out) :: lambda, mu

  lambda = c0(1,1,2,2)
  mu = c0(2,3,2,3)

end

function force_isotropic(c66) result(c66iso)
  use global, only : id6
  double precision, intent(in) :: c66(6,6)
  double precision :: dum, dum1, c66iso(6,6)
  integer :: i

  dum = 0.0
  do i = 1,5
    dum = dum + c66(i,i)
  enddo
  dum = dum/5.0
  dum1 = c66(6,6)
  c66iso = 0.0
  c66iso = id6*dum
  c66iso(6,6) = dum1

end

function Green_function(dx, lambda, mu) result(G)
  use global, only : id3, pi, periodic, npts1, npts2, npts3, dx_min, coef_reg
  double precision, intent(in) :: dx(3), lambda, mu
  double precision :: G (3,3), GxR(3,3), dxp(3), r, one_over_r_reg
  integer :: i,j,k
  
  if (periodic) then
    G = Green_function_per_interp(dx)
  else
    ! GxR = Green_function_x_R(dx, lambda, mu)
    ! G = GxR/maxval([norm2(dx),dx_min])
    ! ! G = GxR/maxval([norm2(dx),sqrt(sqrt(norm2(dx)))])

    GxR = Green_function_x_R(dx, lambda, mu)
    r = norm2(dx)
    if (r >= dx_min) then
      G = GxR/r
    else
      if (r >= 1.0e-2) then
        one_over_r_reg = coef_reg(1)*r**3 + coef_reg(2)*r**2 + coef_reg(3)*r + coef_reg(4)
        G = GxR*one_over_r_reg
      else
        G = id3*coef_reg(4)*((lambda + 3.0*mu) + (lambda + mu)/3.0)/(8.0*pi*mu*(lambda + 2.0*mu))
      endif
    endif

  endif

end

function Green_function_x_R(dx, lambda, mu) result(GxR)
  use global, only : id3, pi
  double precision, intent(in) :: dx(3), lambda, mu
  double precision :: GxR(3,3), ee(3,3)
  double precision :: R, e(3)
  integer :: i,j

  R = norm2(dx)
  e = dx/R

  do i = 1,3
    do j = 1,3
      ee(i,j) = e(i)*e(j)
    enddo
  enddo

  GxR = ((lambda + 3.0*mu)*id3 + (lambda + mu)*ee)/(8.0*pi*mu*(lambda + 2.0*mu))

end

function Green_function_fast(dx) result(G)
  use global, only : id3, pi, G_coef1, G_coef2
  double precision, intent(in) :: dx(3)
  double precision :: G(3,3), ee(3,3)
  double precision :: R, e(3)
  integer :: i,j

  R = norm2(dx)
  e = dx/R
  G(1,1) = (G_coef1 + G_coef2*e(1)*e(1))/R
  G(2,2) = (G_coef1 + G_coef2*e(2)*e(2))/R
  G(3,3) = (G_coef1 + G_coef2*e(3)*e(3))/R
  G(1,2) = G_coef2*e(1)*e(2)/R
  G(1,3) = G_coef2*e(1)*e(3)/R
  G(2,3) = G_coef2*e(2)*e(3)/R
  G(2,1) = G(1,2)
  G(3,1) = G(1,3)
  G(3,2) = G(2,3)

end

function Green_function_fast_mod(dx) result(G)
  use global, only : id3, pi, G_coef1, G_coef2
  double precision, intent(in) :: dx(3)
  double precision :: G(3,3), ee(3,3)
  double precision :: R, e(3), Rmod
  integer :: i,j

  R = norm2(dx)
  Rmod = R + (0.5 - sign(0.5,R - 1.0e-10))*1.0e40
  e = dx/Rmod
  G(1,1) = (G_coef1 + G_coef2*e(1)*e(1))/Rmod
  G(2,2) = (G_coef1 + G_coef2*e(2)*e(2))/Rmod
  G(3,3) = (G_coef1 + G_coef2*e(3)*e(3))/Rmod
  G(1,2) = G_coef2*e(1)*e(2)/Rmod
  G(1,3) = G_coef2*e(1)*e(3)/Rmod
  G(2,3) = G_coef2*e(2)*e(3)/Rmod
  G(2,1) = G(1,2)
  G(3,1) = G(1,3)
  G(3,2) = G(2,3)

end

subroutine test_intgpt(intg_pt)
  use global
  use types
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)

  double precision :: x1(3), x2(3), x3(3), aux33(3,3), dum
  double precision :: xintpt(3,nintg_pt), fintpt(nintg_pt), jacob
  integer :: i 

  ! three points on a plane, rotated and translated in space
  x1 = 0.0
  x2 = [0.0,1.0,0.0]*0.25
  x3 = [0.0,1.0,1.0]*0.25
  aux33(1,:) = [1.0, 0.0, 0.0]
  aux33(2,:) = [0.0, cos(34.0*pi/180.0), -sin(34.0*pi/180.0)]
  aux33(3,:) = [0.0, sin(34.0*pi/180.0), cos(34.0*pi/180.0)]
  x1 = matmul(aux33, x1)
  x2 = matmul(aux33, x2)
  x3 = matmul(aux33, x3)
  aux33(1,:) = [cos(65.0*pi/180.0), 0.0, -sin(65.0*pi/180.0)]
  aux33(2,:) = [0.0, 1.0, 0.0]
  aux33(3,:) = [sin(65.0*pi/180.0), 0.0, cos(65.0*pi/180.0)]
  x1 = matmul(aux33, x1)
  x2 = matmul(aux33, x2)
  x3 = matmul(aux33, x3)
  aux33(1,:) = [cos(22.0*pi/180.0), -sin(22.0*pi/180.0), 0.0]
  aux33(2,:) = [sin(22.0*pi/180.0), cos(22.0*pi/180.0), 0.0]
  aux33(3,:) = [0.0, 0.0, 1.0]
  x1 = matmul(aux33, x1)
  x2 = matmul(aux33, x2)
  x3 = matmul(aux33, x3)

  x1 = x1 + [7.0,6.0,-4.5]
  x2 = x2 + [7.0,6.0,-4.5]
  x3 = x3 + [7.0,6.0,-4.5]

  ! integrate 1/R over triangle (with singularity at node 1)
  call intg_pt_to_current(intg_pt, x1, x2, x3, xintpt)
  call jacobian(x1, x2, x3, jacob)
  fintpt = 1.0
  dum = 0.0
  do i = 1, nintg_pt
    dum = dum + norm2(intg_pt(i)%Xhat)/norm2(xintpt(:,i) - x1)*fintpt(i)*jacob*intg_pt(i)%w
  enddo
  
  write(*,*) dum
end

function integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu) result(Gintg) 
  use types
  use global, only : nintg_pt
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  double precision, intent(in) :: x1(3), x2(3), x3(3), lambda, mu
  double precision :: Gintg(3,3), dxintpt(3,nintg_pt), GxR(3,3), dx(3), dx1(3), dx2(3), dx3(3)
  double precision ::jacob
  integer :: ip

  ! shift so that node 1 is origin
  ! dx1 = x1 - x1
  ! dx2 = x2 - x1
  ! dx3 = x3 - x1
  dx1 = x1 - x1 ! x1 is singular point
  dx2 = x1 - x2 
  dx3 = x1 - x3 

  ! integration points in triangle
  call intg_pt_to_current(intg_pt, dx1, dx2, dx3, dxintpt)

  ! Jacobian of transformation
  call jacobian(x1, x2, x3, jacob)

  ! integral as weighted sum
  Gintg = 0.0
  do ip = 1, nintg_pt
    GxR = Green_function_x_R(dxintpt(:,ip), lambda, mu)
    Gintg = Gintg + intg_pt(ip)%Xhatn/norm2(dxintpt(:,ip))*GxR*jacob*intg_pt(ip)%w
  enddo

end

subroutine test_Green_function_integration(intg_pt)
  use global
  use types
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)

  double precision :: x1(3), x2(3), x3(3), G_Gauss(3,3), Gintg(3,3)
  double precision :: xintpt(3,nintg_pt), fintpt(nintg_pt), jacob, x1pt, x2pt, x3pt
  double precision :: dx2, dx3, x(3)
  integer :: i, ip1, ip2, ip3, np1, np2, np3

  x1 = 0.0
  x2 = [0.0,1.0,0.0]*0.25
  x3 = [0.0,1.0,1.0]*0.25

  G_Gauss = integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)
  G_Gauss = integrate_Green_function(x1, x3, x2, intg_pt, lambda, mu) ! flip nodes, should give same
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)

  np2 = 1000
  np3 = 1000
  dx2 = 0.25/float(np2)
  dx3 = 0.25/float(np3)
  x2pt = 0.0
  x3pt = 0.0
  x1pt = 0.0
  Gintg = 0.0
  do ip2 = 1,np2
    x2pt = float(ip2)*dx2 - dx2*0.5
    do ip3 = 1,ip2
      x3pt = float(ip3)*dx3 - dx3*0.5
      x(1) = x1pt
      x(2) = x2pt
      x(3) = x3pt
      ! write(*,*)x2pt,x3pt
      Gintg = Gintg + Green_function(-x, lambda, mu)*dx2*dx3
    enddo
  enddo
  write(*,'(3E16.8)') (Gintg(i,:),i=1,3)

  ! stop

end


subroutine integrate_Green_function_edges(elements, faces, edges, intg_pt, intg_pt_quad)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1face(3), x2face(3), x3face(3), xnode_face(3,3)
  double precision :: x1(3), x2(3), x3(3), dx(3), xedge(3), xshift(3)
  integer :: iedg, i ,ifc, iel, j, face_edge

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,i,ifc,iel,x1face,x2face,x3face,xnode_face,face_edge,xedge,xshift) &
  !$OMP PRIVATE(x1,x2,x3) &
  !$OMP SHARED(nedg,edges,faces,elements,periodic,intg_pt_quad,lambda,mu) &
  !$OMP SHARED(face_edge_nodes,face_edge_nodes_remaining)
  !$OMP DO 
  do iedg = 1, nedg

    do i = 1, edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      iel = faces(ifc)%el(1)
      call get_face_nodes(faces(ifc),elements(iel),x1face, x2face, x3face)
      xnode_face(:,1) = x1face
      xnode_face(:,2) = x2face
      xnode_face(:,3) = x3face
      face_edge = edges(iedg)%face_edge(i)

      if (periodic) then
        xedge = 0.5*(xnode_face(:,face_edge_nodes(1,face_edge)) + xnode_face(:,face_edge_nodes(2,face_edge)))
        xshift = - xedge + edges(iedg)%xc
        xnode_face(:,1) = xnode_face(:,1) + xshift
        xnode_face(:,2) = xnode_face(:,2) + xshift
        xnode_face(:,3) = xnode_face(:,3) + xshift
        ! if (norm2(xshift) > 1.0e-5) then
        !   write(*,*) iedg, i, xshift
        !   write(*,*) edges(iedg)%xc
        !   write(*,*) xedge
        !   write(*,*) xnode_face(:,1)
        !   write(*,*) xnode_face(:,2)
        !   write(*,*) xnode_face(:,3)
        ! endif
      endif

      ! two triangles
      x1 = edges(iedg)%xc ! node with singularity (x1) is the center of edge
      x2 = xnode_face(:,face_edge_nodes_remaining(face_edge)) ! x2 is the opposite node on the face

      ! 1
      x3 = xnode_face(:,face_edge_nodes(1,face_edge))
      if (norm2((x2 + x3)*0.5 - x1) < 1.0e-7) stop 'nodes not properly assigned' ! if x2 is the opposite node on the face this cannot be true
      ! edges(iedg)%face_G(:,:,i) = integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)
      edges(iedg)%face_G(:,:,i) = integrate_Green_function_quad(x1, x2, x3, intg_pt_quad, lambda, mu)

      ! 2
      x3 = xnode_face(:,face_edge_nodes(2,face_edge))
      if (norm2((x2 + x3)*0.5 - x1) < 1.0e-7) stop 'nodes not properly assigned'
      ! edges(iedg)%face_G(:,:,i) = edges(iedg)%face_G(:,:,i) + integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)
      edges(iedg)%face_G(:,:,i) = edges(iedg)%face_G(:,:,i) + integrate_Green_function_quad(x1, x2, x3, intg_pt_quad, lambda, mu)
      edges(iedg)%face_G(:,:,i) = edges(iedg)%face_G(:,:,i)/faces(ifc)%area

!       write(*,*) iedg, i, (edges(iedg)%face_G(:,:,i))
! 
!       ! test
!       dx = edges(iedg)%xc - faces(ifc)%x
!       edges(iedg)%face_G(:,:,i) = Green_function(dx, lambda, mu)
! 
!       write(*,*) iedg, i, (edges(iedg)%face_G(:,:,i))
      
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! iedg = 1522
  ! do i = 1, edges(iedg)%nface
  !   write(*,*) i
  !   write(*,'(3E16.8)') (edges(iedg)%face_G(j,:,i), j = 1,3)
  ! enddo

end

subroutine calc_F_faces(elements, faces)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)

  double precision :: dpolar(3,3)
  integer :: iel1, iel2, ifc

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,iel1,iel2,dpolar) &
  !$OMP SHARED(nfc,elements,faces,periodic) 
  !$OMP DO 
  do ifc = 1, nfc

    iel1 = faces(ifc)%el(1)
 
    ! make sure it is outward normal to iel1
    if (dot_product(- elements(iel1)%x + faces(ifc)%x, faces(ifc)%normal) < 0.0 .and. (.not. periodic)) & ! will not hold for boundary faces in periodic case
      stop 'not outward normal w.r.t. el(1)'

    if (faces(ifc)%nel == 1) then
      ! faces(ifc)%F = matmul(elements(iel1)%polar, faces(ifc)%normal)*faces(ifc)%area
      ! faces(ifc)%F = faces(ifc)%F + matmul(elements(iel1)%stress - Sapp, faces(ifc)%normal)*faces(ifc)%area*0.1
      ! faces(ifc)%F = matmul(elements(iel1)%stress - Sapp, faces(ifc)%normal)*faces(ifc)%area*k_penalty + faces(ifc)%Fbc
      ! faces(ifc)%F = matmul(elements(iel1)%polar, faces(ifc)%normal)*faces(ifc)%area + faces(ifc)%Fbc
      ! faces(ifc)%F = matmul(elements(iel1)%polar, faces(ifc)%normal)*faces(ifc)%area + faces(ifc)%u*k_penalty + faces(ifc)%Fbc
      ! faces(ifc)%F = faces(ifc)%u*100.0 + faces(ifc)%Fbc
      faces(ifc)%F = - faces(ifc)%Fbc
      ! faces(ifc)%F = - faces(ifc)%u*k_penalty
      ! faces(ifc)%F = 0.0
      ! faces(ifc)%F = matmul(Sapp, faces(ifc)%normal)*faces(ifc)%area
    elseif (faces(ifc)%nel == 2) then
      iel2 = faces(ifc)%el(2)
      dpolar = mandel_t1_to_t2(elements(iel1)%polar - elements(iel2)%polar)
      faces(ifc)%F = matmul(dpolar, faces(ifc)%normal)*faces(ifc)%area ! minus for iel2 because normal is outward with respect to iel1
      ! make sure it is outward normal to iel1
      if (dot_product(- elements(iel2)%x + faces(ifc)%x, faces(ifc)%normal) > 0.0 .and. (.not. periodic)) then ! will not hold for boundary faces in periodic case
        ! write(*,*) elements(iel1)%xnode(:,1)
        ! write(*,*) elements(iel1)%xnode(:,2)
        ! write(*,*) elements(iel1)%xnode(:,3)
        ! write(*,*) elements(iel1)%xnode(:,4)
        ! write(*,*) elements(iel2)%xnode(:,1)
        ! write(*,*) elements(iel2)%xnode(:,2)
        ! write(*,*) elements(iel2)%xnode(:,3)
        ! write(*,*) elements(iel2)%xnode(:,4)
        ! write(*,*) dot_product(- elements(iel2)%x + faces(ifc)%x, faces(ifc)%normal) 
        ! write(*,*) iel2
        ! write(*,*) faces(ifc)%normal
        ! write(*,*) elements(iel2)%x
        ! write(*,*) - elements(iel2)%x + faces(ifc)%x
        stop 'outward normal w.r.t. el(2)'

      endif
    endif
    ! write(*,*) ifc, faces(ifc)%F

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end


subroutine calc_F_faces_box(elements, faces, boxes)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: fact1, fact2, dFn_sum, dpolar(3,3)
  double precision, allocatable :: face_F(:,:)
  integer :: iel1, iel2, ifc, ib, i


  allocate(face_F(3,nfc))
  face_F = 0.0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel1, iel2, ifc, ib, i, fact1, fact2, dpolar) &
  !$OMP SHARED(face_F,elements, faces, boxes, nbox, periodic) 
  !$OMP DO 
  do ib = 1, nbox

    do i = 1,boxes(ib)%nfc
      ifc = boxes(ib)%face_ids(i)

      if (faces(ifc)%boundary) then

        boxes(ib)%face_F(:,i) = - faces(ifc)%Fbc

      else
      
        iel1 = faces(ifc)%el(1)
        if (elements(iel1)%box_id == ib) then
          fact1 = 1.0
        else
          fact1 = 0.0
        endif
  
        if (faces(ifc)%nel == 1) then
          boxes(ib)%face_F(:,i) = - faces(ifc)%Fbc
        elseif (faces(ifc)%nel == 2) then
          iel2 = faces(ifc)%el(2)
          if (elements(iel2)%box_id == ib) then
            fact2 = 1.0
          else
            fact2 = 0.0
          endif
          dpolar = mandel_t1_to_t2(elements(iel1)%polar*fact1 - elements(iel2)%polar*fact2)
          boxes(ib)%face_F(:,i) = matmul(dpolar, faces(ifc)%normal)*faces(ifc)%area ! minus for iel2 because normal is outward with respect to iel1
          
          ! make sure it is outward normal to iel1
          if (dot_product(- elements(iel2)%x + faces(ifc)%x, faces(ifc)%normal) > 0.0 .and. (.not. periodic)) then ! will not hold for boundary faces in periodic case
            stop 'outward normal w.r.t. el(2)'
          endif
        endif

      endif

      ! test
      face_F(:,ifc) = face_F(:,ifc) + boxes(ib)%face_F(:,i)
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ! verify
  ! dFn_sum = 0.0
  ! do ifc = 1,nfc
  !   dFn_sum = dFn_sum + norm2(face_F(:,ifc) - faces(ifc)%F)
  ! enddo
  ! deallocate(face_F)
  ! write(*,*) 'dFn_sum', dFn_sum
  ! ! stop

end

subroutine calc_f_edges(faces, edges)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  
  double precision :: f(3), wedg = 1.0/3.0, favg(3)
  integer :: iedg, ifc, i

  ! favg = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,f,i,ifc) &
  !$OMP SHARED(nedg,edges,faces,wedg) 
  !$OMP DO 
  do iedg = 1, nedg
    ! write(*,*) iedg

    f = 0.0
    do i = 1, edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      f = f + wedg*faces(ifc)%F
      ! write(*,*) ifc, faces(ifc)%F
    enddo
    edges(iedg)%f = f
    ! favg = favg + f
    ! write(*,*) edges(iedg)%f
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
!  favg = favg/float(nedg)
!  write(*,*) favg
!  
!  do iedg = 1, nedg
!    edges(iedg)%f = edges(iedg)%f - favg
!  enddo

end


subroutine calc_f_edges_box(faces, edges, boxes)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (in) :: edges(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)
  
  double precision :: f(3), wedg = 1.0/3.0, favg(3), dfn_sum
  double precision, allocatable :: edge_f(:,:)
  integer :: iedg, ifc, i, j, box_id_loc, ib

  allocate(edge_f(3,nedg))
  edge_f = 0.0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg, ifc, i, j, box_id_loc, ib, f) &
  !$OMP SHARED(faces, edges, boxes, wedg, edge_f, nbox) 
  !$OMP DO 
  do ib = 1, nbox

    do i = 1,boxes(ib)%nedg
      iedg = boxes(ib)%edge_ids(i)

      f = 0.0
      do j = 1, edges(iedg)%nface
        ifc = edges(iedg)%face(j)
        if (faces(ifc)%box_id(1) == ib) then
          box_id_loc = faces(ifc)%box_id_loc(1)
          ! write(*,*) ifc, boxes(ib)%face_ids(box_id_loc)
          f = f + wedg*boxes(ib)%face_F(:,box_id_loc)
        elseif (faces(ifc)%nbox == 2 .and. faces(ifc)%box_id(2) == ib) then
          box_id_loc = faces(ifc)%box_id_loc(2)
          ! write(*,*) ifc, boxes(ib)%face_ids(box_id_loc), box_id_loc
          f = f + wedg*boxes(ib)%face_F(:,box_id_loc)
        endif
      enddo
      boxes(ib)%edge_f(:,i) = f
      edge_f(:,iedg) = edge_f(:,iedg) + boxes(ib)%edge_f(:,i)

    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! ! verify
  ! dfn_sum = 0.0
  ! do iedg = 1,nedg
  !   dfn_sum = dfn_sum + norm2(edge_f(:,iedg) - edges(iedg)%f)
  !   ! if (norm2(edge_f(:,iedg) - edges(iedg)%f) > 0.0) write(*,*)edge_f(:,iedg), edges(iedg)%f
  ! enddo
  ! deallocate(edge_f)
  ! write(*,*) 'dfn_sum', dfn_sum
  !  ! stop

end

subroutine convolution_edges(faces, edges)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)

  double precision :: uacc(3), G(3,3), dx(3), wedg = 1.0/3.0, f(3), Fface(3)
  integer :: iedg1, iedg2, i, j, ifc, ib1, iel1, iel2 

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,uacc,iedg2,dx,G,i,j,ifc,f) &
  !$OMP SHARED(nedg,edges,lambda,mu,faces,wedg,periodic) 
  !$OMP DO 
  do iedg1 = 1,nedg

    ! sum over all
    uacc = 0.0
    do iedg2 = 1,nedg

      if (iedg1 .ne. iedg2) then
        dx = edges(iedg1)%xc - edges(iedg2)%xc
        G = Green_function(dx, lambda, mu)
        uacc = uacc + matmul(G, edges(iedg2)%f)
      endif

    enddo

    do i = 1, edges(iedg1)%nface ! loop over faces containing this edge
      ifc = edges(iedg1)%face(i)

      ! remove effects of faces containing this edge
      do j = 1, 3 ! loop over edges in face
        iedg2 = faces(ifc)%edge(j)
        if (iedg2 .ne. iedg1) then
          dx = edges(iedg1)%xc - edges(iedg2)%xc
          G = Green_function(dx, lambda, mu)
          f = wedg*faces(ifc)%F ! contribution to edge force from this face
          uacc = uacc - matmul(G, f)
        endif
      enddo

      ! add true effect of faces containing this edge
      G = edges(iedg1)%face_G(:,:,i)
      uacc = uacc + matmul(G, faces(ifc)%F)

    enddo

    if (allocated(edges(iedg1)%face_close)) then
      
      do i = 1, size(edges(iedg1)%face_close) ! loop over close faces
        ifc = edges(iedg1)%face_close(i)
  
        ! remove effects of faces containing this edge
        do j = 1, 3 ! loop over edges in face
          iedg2 = faces(ifc)%edge(j)
          if (iedg2 .ne. iedg1) then
            dx = edges(iedg1)%xc - edges(iedg2)%xc
            G = Green_function(dx, lambda, mu)
            f = wedg*faces(ifc)%F ! contribution to edge force from this face
            uacc = uacc - matmul(G, f)
          endif
        enddo
  
        ! add true effect of faces containing this edge
        G = edges(iedg1)%face_close_G(:,:,i)
        uacc = uacc + matmul(G, faces(ifc)%F)
  
      enddo

    endif

    edges(iedg1)%u = - uacc
    ! write(*,*) iedg1, edges(iedg1)%u

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine convolution_edges_box(elements, faces, edges, boxes)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (in) :: faces(:)
  type (element_type), allocatable, intent (in) :: elements(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  double precision :: uacc(3), G(3,3), dx(3), wedg = 1.0/3.0, f(3), Fface(3)
  double precision :: dun_sum, dG_dX(3,3,3), polar(3,3)
  integer :: iedg1, iedg2, i, j, ifc, ib1, iel1, iel2, iedg_in_box, ib2, ib_ind_near
  integer :: ipt, ib_ind_far, iel_in_box, ifc2

  dun_sum = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(uacc, G, dx, f, Fface) &
  !$OMP PRIVATE(iedg1, iedg2, i, j, ifc, ib1, iel1, iel2, iedg_in_box, ib2, ib_ind_near) &
  !$OMP PRIVATE(ipt, ib_ind_far, iel_in_box, ifc2, dG_dX, polar) &
  !$OMP SHARED(nedg,faces,elements,edges,boxes,wedg,mu,lambda) 
  !$OMP DO 
  do iedg1 = 1,nedg

    ib1 = edges(iedg1)%box_id_unique

    ! sum over all
    uacc = 0.0
    ! do ib2 = 1, nbox
    do ib_ind_near = 1, boxes(ib1)%nnear
      ib2 = boxes(ib1)%near_box_id(ib_ind_near)
      do iedg_in_box = 1,boxes(ib2)%nedg
        iedg2 = boxes(ib2)%edge_ids(iedg_in_box)

        if (iedg1 .ne. iedg2) then
          dx = edges(iedg1)%xc - edges(iedg2)%xc
          ! G = Green_function(dx, lambda, mu)
          G = Green_function_fast(dx)
          uacc = uacc + matmul(G, boxes(ib2)%edge_f(:,iedg_in_box))
        endif

      enddo
    enddo

    do i = 1, edges(iedg1)%nface ! loop over faces containing this edge
      ifc = edges(iedg1)%face(i)

      ! remove effects of faces containing this edge
      do j = 1, 3 ! loop over edges in face
        iedg2 = faces(ifc)%edge(j)
        if (iedg2 .ne. iedg1) then
          dx = edges(iedg1)%xc - edges(iedg2)%xc
          G = Green_function(dx, lambda, mu)
          f = wedg*faces(ifc)%F ! contribution to edge force from this face
          uacc = uacc - matmul(G, f)
        endif
      enddo

      ! add true effect of faces containing this edge
      G = edges(iedg1)%face_G(:,:,i)
      uacc = uacc + matmul(G, faces(ifc)%F)

    enddo

    if (allocated(edges(iedg1)%face_close)) then
      
      do i = 1, size(edges(iedg1)%face_close) ! loop over close faces
        ifc = edges(iedg1)%face_close(i)
  
        ! remove effects of faces containing this edge
        do j = 1, 3 ! loop over edges in face
          iedg2 = faces(ifc)%edge(j)
          if (iedg2 .ne. iedg1) then
            dx = edges(iedg1)%xc - edges(iedg2)%xc
            G = Green_function(dx, lambda, mu)
            f = wedg*faces(ifc)%F ! contribution to edge force from this face
            uacc = uacc - matmul(G, f)
          endif
        enddo
  
        ! add true effect of faces containing this edge
        G = edges(iedg1)%face_close_G(:,:,i)
        uacc = uacc + matmul(G, faces(ifc)%F)
  
      enddo

    endif
    uacc = - uacc

    ! distant interaction sum
    do ib_ind_far = 1, boxes(ib1)%nfar
      ib2 = boxes(ib1)%far_box_id(ib_ind_far)
      do iel_in_box = 1,boxes(ib2)%nel
        iel2 = boxes(ib2)%element_ids(iel_in_box)

        dG_dX = 0.0
        do ipt = 1,4
          dx = edges(iedg1)%xc - elements(iel2)%xintpt(:,ipt)
          dG_dX = dG_dX + gradGreen_function_num(dx, lambda, mu)*elements(iel2)%w(ipt)
        enddo
        polar = mandel_t1_to_t2(elements(iel2)%polar)
        uacc = uacc + TijkTjk(dG_dX, polar*elements(iel2)%v)

        do i = 1,4
          ifc2 = elements(iel2)%face(i)
          if (faces(ifc2)%nel == 1) then
            G = 0.0
            do j = 1,3
              iedg2 = faces(ifc2)%edge(j)
              dx = edges(iedg1)%xc - edges(iedg2)%xc
              G = G + Green_function(dx, lambda, mu)*wedg
            enddo
            F = matmul(polar,faces(ifc2)%normal)*faces(ifc2)%area + faces(ifc2)%Fbc
            uacc = uacc + matmul(G, F)
          endif
        enddo

      enddo
    enddo

    edges(iedg1)%u = uacc
    ! dun_sum = dun_sum + norm2(edges(iedg1)%u + uacc)
    
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! write(*,*) 'dun_sum', dun_sum
  ! stop

end

subroutine calc_u_faces(edges, faces)
  use global
  use types
  use tensor_functions

  type (edge_type), allocatable, intent (in) :: edges(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  
  double precision :: wedg = 1.0/3.0, uacc(3)
  integer :: iedg, ifc, i

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,uacc,i,iedg) &
  !$OMP SHARED(nfc,wedg,faces,edges) 
  !$OMP DO 
  do ifc = 1, nfc
    uacc = 0.0
    do i = 1,3
      iedg = faces(ifc)%edge(i)
      uacc = uacc + edges(iedg)%u*faces(ifc)%area*wedg
    enddo
    faces(ifc)%u = uacc
    ! write(*,*) ifc, faces(ifc)%u
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

subroutine calc_disgrad_elements(faces, elements, boxes)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (in) :: faces(:)
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (box_type), allocatable, intent (in) :: boxes(:)
  
  double precision :: u(3), n(3), disgrad(3,3), sgn, Gamma_self_polar(3,3), polar(3,3)
  integer :: iel, ifc, i, j, k, ib

  ! strainavg = 0.0
  ! vtot = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,disgrad,i,ifc,sgn,u,n,Gamma_self_polar,polar) &
  !$OMP SHARED(nel,elements,faces,periodic,iJacobi,Eapp,polaravg,plasticity,solution_method,Edotapp,time_inc) &
  !$OMP REDUCTION(+:strainavg,vtot)
  !$OMP DO 
  do iel = 1, nel

    disgrad = 0.0
    do i = 1,4
      ifc = elements(iel)%face(i)
      sgn = elements(iel)%face_normal_sgn(i)
      u = faces(ifc)%u
      n = faces(ifc)%normal*sgn
      do j = 1,3
        do k = 1,3
          disgrad(j,k) = disgrad(j,k) + u(j)*n(k)
        enddo
      enddo

      if (dot_product(- elements(iel)%x + faces(ifc)%x, n) < 0.0 .and. (.not. periodic)) & ! will not hold for boundary faces in periodic case
        stop 'not outward normal w.r.t. element'
      
    enddo
    elements(iel)%disgrad = elements(iel)%disgradt + disgrad/elements(iel)%v + Edotapp*time_inc
    if (iJacobi == 1 .and. (.not.plasticity) .and. solution_method == 1) then
      polar = mandel_t1_to_t2(elements(iel)%polar + polaravg)
      Gamma_self_polar = TijklTkl(elements(iel)%Gamma_self, polar)
      elements(iel)%disgrad = TijklTkl(elements(iel)%I_GdC_inv,elements(iel)%disgrad - Gamma_self_polar)
    endif
    elements(iel)%strain = mandel_t2_to_t1(sym33(elements(iel)%disgrad))

    strainavg = strainavg + elements(iel)%strain*elements(iel)%v
    vtot = vtot + elements(iel)%v

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  strainavg = strainavg/vtot

  ! call correct_distorted(elements)

  strainavg = 0.0
  vtot = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(nel,elements) &
  !$OMP REDUCTION(+:strainavg,vtot)
  !$OMP DO 
  do iel = 1, nel

    strainavg = strainavg + elements(iel)%strain*elements(iel)%v
    vtot = vtot + elements(iel)%v
    ! write(*,*) iel, elements(iel)%strain

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  strainavg = strainavg/vtot

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(nel,elements,strainavg,Eapp)
  !$OMP DO 
  do iel = 1, nel
    elements(iel)%strain = elements(iel)%strain - strainavg + mandel_t2_to_t1(Eapp)
    elements(iel)%disgrad = elements(iel)%disgrad - mandel_t1_to_t2(strainavg) + Eapp
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

function Green_operator(x, lambda, mu) result(Goper)
  use global, only : id3, pi
  double precision, intent(in) :: x(3), lambda, mu
  double precision :: Goper(3,3,3,3), R
  integer :: i, j, k, l

  R = norm2(x)

  Goper=0.0
  do i = 1,3
  do j = 1,3
  do k = 1,3
  do l = 1,3
  Goper(i,j,k,l) = &
   -1.0/8.0/pi/mu*(id3(k,l)*R**(-3)-3*x(k)*x(l)*R**(-5))* &
   ((lambda+3.0*mu)/(lambda+2.0*mu)*id3(i,j)+ &
   (lambda+mu)/(lambda+2.0*mu)*x(i)*x(j)*R**(-2))- &
   1.0/8.0/pi/mu*(lambda+mu)/(lambda+2.0*mu)*x(k)*R**(-3)* &
   (id3(i,l)*x(j)*R**(-2)+x(i)*id3(j,l)*R**(-2) &
   -2.0*x(i)*x(j)*x(l)*R**(-4))+ &
   1.0/8.0/pi/mu*(lambda+mu)/(lambda+2.0*mu)* &
   (id3(i,k)*id3(j,l)*R**(-3)-3.0*id3(i,k)*x(j)*x(l)*R**(-5)+ &
   id3(j,k)*id3(i,l)*R**(-3)-3.0*id3(j,k)*x(i)*x(l)*R**(-5) &
   -2.0*id3(i,l)*x(j)*x(k)*R**(-5)- &
   2.0*x(i)*id3(j,l)*x(k)*R**(-5)- &
   2.0*x(i)*x(j)*id3(k,l)*R**(-5)+ &
   10.0*x(i)*x(j)*x(k)*x(l)*R**(-7))
  enddo
  enddo
  enddo
  enddo


end

function Green_operator_num(dx, lambda, mu) result(Goper)
  use global, only : id3, pi
  double precision, intent(in) :: dx(3), lambda, mu
  double precision :: Goper(3,3,3,3), inc, dir1_inc(3), dir2_inc(3), dir3_inc(3)
  integer :: i, j, k, l

  inc = 1.0e-4
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc

  Goper(:,:,1,1) = (Green_function(dx + dir1_inc, lambda, mu) &
  - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir1_inc, lambda, mu))/(inc**2)
  Goper(:,:,2,2) = (Green_function(dx + dir2_inc, lambda, mu) &
  - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir2_inc, lambda, mu))/(inc**2)
  Goper(:,:,3,3) = (Green_function(dx + dir3_inc, lambda, mu) &
  - 2.0*Green_function(dx, lambda, mu) + Green_function(dx - dir3_inc, lambda, mu))/(inc**2)
  
  Goper(:,:,1,2) = &
  (Green_function(dx + dir1_inc + dir2_inc, lambda, mu) - Green_function(dx + dir1_inc - dir2_inc, lambda, mu) - &
   Green_function(dx - dir1_inc + dir2_inc, lambda, mu) + Green_function(dx - dir1_inc - dir2_inc, lambda, mu))/(4.0*inc**2)
  
  Goper(:,:,1,3) = &
  (Green_function(dx + dir1_inc + dir3_inc, lambda, mu) - Green_function(dx + dir1_inc - dir3_inc, lambda, mu) - &
   Green_function(dx - dir1_inc + dir3_inc, lambda, mu) + Green_function(dx - dir1_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
  
  Goper(:,:,2,3) = &
  (Green_function(dx + dir2_inc + dir3_inc, lambda, mu) - Green_function(dx + dir2_inc - dir3_inc, lambda, mu) - &
   Green_function(dx - dir2_inc + dir3_inc, lambda, mu) + Green_function(dx - dir2_inc - dir3_inc, lambda, mu))/(4.0*inc**2)
  
  Goper(:,:,2,1) = Goper(:,:,1,2)
  Goper(:,:,3,1) = Goper(:,:,1,3)
  Goper(:,:,3,2) = Goper(:,:,2,3)

end

function gradGreen_function_num(dx, lambda, mu) result(dG_dX)
  use global, only : id3, pi
  double precision, intent(in) :: dx(3), lambda, mu
  double precision :: dG_dX(3,3,3), inc, dir1_inc(3), dir2_inc(3), dir3_inc(3)
  integer :: i, j, k, l

  inc = 1.0e-4
  dir1_inc = [1.0,0.0,0.0]*inc
  dir2_inc = [0.0,1.0,0.0]*inc
  dir3_inc = [0.0,0.0,1.0]*inc

  dG_dX(:,:,1) = (Green_function(dx + dir1_inc, lambda, mu) - &
    Green_function(dx - dir1_inc, lambda, mu))/(2.0*inc)
  dG_dX(:,:,2) = (Green_function(dx + dir2_inc, lambda, mu) - &
    Green_function(dx - dir2_inc, lambda, mu))/(2.0*inc)
  dG_dX(:,:,3) = (Green_function(dx + dir3_inc, lambda, mu) - &
    Green_function(dx - dir3_inc, lambda, mu))/(2.0*inc)

end

subroutine uniqueInt(a, n, b, nun, ndim)

  implicit none
  integer, intent(in) :: n, ndim
  integer, dimension(ndim), intent(in) :: a
  integer, dimension(ndim), intent(out) :: b
  integer, intent(out) :: nun
  integer :: min_val, max_val
  integer :: i

  i = 0
  b = 0
  min_val = minval(a(1:n)) - 1
  max_val = maxval(a(1:n))
  do while (min_val<max_val)
      i = i + 1
      min_val = minval(a(1:n), mask = a(1:n) > min_val)
      b(i) = min_val
  enddo
  nun = i

end

subroutine convolution_all(boxes, elements, faces, edges)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (in) :: boxes(:)
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)

  double precision :: uacc(3), G(3,3), dx(3), wedg = 1.0/3.0, f(3), Fface(3), u(3), n(3), disgrad(3,3), sgn, dxdx(3,3), fact
  double precision :: velunited, Gamma_self_polar(3,3), dxdxdx(3,3,3),dxdxdxdx(3,3,3,3),polar(3,3)
  integer :: iedg1, iedg2, i, j, ifc, ib, iel1, iel2, iedg, l, k, iel, ii, jj, ind, kk, ll, tedg, tf1, tf2, texp, t1, t2

  tedg = 0
  tf1 = 0
  tf2 = 0
  texp = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(uacc,G,dx,f,Fface,u,n,disgrad,sgn,dxdx,fact,ind,t1,t2,polar) &
  !$OMP PRIVATE(iedg1,iedg2,i,j,ifc,ib,iel1,iel2,iedg,l,k,iel,ii,jj,u_edge,Gamma_self_polar,dxdxdx,dxdxdxdx,ll) &
  !$OMP SHARED(boxes,lambda,mu,edges,faces,elements,wedg,order,nbox,periodic,iJacobi,Eapp,polaravg,precalc_G) &
  !$OMP SHARED(plasticity,solution_method) &
  !$OMP REDUCTION(+:tedg,tf1,tf2,texp)
  !$OMP DO
  do ib = 1,nbox
  
    ! loop over edges in box
    ! !$OMP PARALLEL DEFAULT(NONE) &
    ! !$OMP PRIVATE(j,iedg1,uacc,i,iedg2,fact,dx,G,ifc,f,dxdx) &
    ! !$OMP SHARED(boxes,ib,lambda,mu,edges,faces,wedg,order,u_edge)
    ! !$OMP DO
    do j = 1,boxes(ib)%nedg ! these are all including shared

      t1 = time()
      iedg1 = boxes(ib)%edge_ids(j)
      ind = findloc(edges(iedg1)%box_id, ib, 1)
      uacc = 0.0
      ! do i = 1,boxes(ib)%nnear_edge
      do i = 1,size(boxes(ib)%near_edge_unique_ids,1) ! these are unique
        ! iedg2 = boxes(ib)%near_edge_ids(i)
        iedg2 = boxes(ib)%near_edge_unique_ids(i)
        ! fact = boxes(ib)%edge_near_factor(i)
        
        ! full convolution over near edges
        if (iedg1 .ne. iedg2) then
          if (.not.(periodic.and.precalc_G)) then
            dx = edges(iedg1)%xc - edges(iedg2)%xc
            G = Green_function(dx, lambda, mu)
          else
            G = edges(iedg1)%G(:,:,ind,i)
          endif
          ! uacc = uacc + matmul(G, edges(iedg2)%f*fact)
          uacc = uacc + matmul(G, edges(iedg2)%f)
        endif
    
      enddo
      t2 = time()
      tedg = tedg + t2 - t1
  
      t1 = time()
      do i = 1, edges(iedg1)%nface ! loop over faces containing this edge
        ifc = edges(iedg1)%face(i)
    
        ! remove effects of faces containing this edge
        do k = 1, 3 ! loop over edges in face
          iedg2 = faces(ifc)%edge(k)
          if (iedg2 .ne. iedg1) then
            dx = edges(iedg1)%xc - edges(iedg2)%xc
            G = Green_function(dx, lambda, mu)
            f = wedg*faces(ifc)%F ! contribution to edge force from this face
            uacc = uacc - matmul(G, f)
          endif
        enddo
    
        ! add true effect of faces containing this edge
        G = edges(iedg1)%face_G(:,:,i)
        uacc = uacc + matmul(G, faces(ifc)%F)
    
      enddo
      t2 = time()
      tf1 = tf1 + t2 - t1

      t1 = time()
      if (allocated(edges(iedg1)%face_close)) then
        do i = 1, size(edges(iedg1)%face_close) ! loop over close faces 
          ifc = edges(iedg1)%face_close(i)
      
          ! remove effects of close faces
          do k = 1, 3 ! loop over edges in face
            iedg2 = faces(ifc)%edge(k)
            if (iedg2 .ne. iedg1) then
              dx = edges(iedg1)%xc - edges(iedg2)%xc
              G = Green_function(dx, lambda, mu)
              f = wedg*faces(ifc)%F ! contribution to edge force from this face
              uacc = uacc - matmul(G, f)
            endif
          enddo
      
          ! add true effect of faces containing this edge
          G = edges(iedg1)%face_close_G(:,:,i)
          uacc = uacc + matmul(G, faces(ifc)%F)
      
        enddo
      endif
      t2 = time()
      tf2 = tf2 + t2 - t1

      t1 = time()
      ! add effects of far field
      dx = edges(iedg1)%xc - boxes(ib)%x
      if (periodic) dx = per_dx(dx)
      do ii = 1,3
        do jj = 1,3
          dxdx(ii,jj) = dx(ii)*dx(jj)
          do kk = 1,3
            dxdxdx(ii,jj,kk) = dxdx(ii,jj)*dx(kk)
            do ll = 1,3
              dxdxdxdx(ii,jj,kk,ll) = dxdxdx(ii,jj,kk)*dx(ll)
            enddo
          enddo
        enddo
      enddo
      ! dx = edges(iedg1)%dx(:,ind)
      ! dxdx = edges(iedg1)%dxdx(:,:,ind)

      if (order == 0) then
        uacc = uacc + boxes(ib)%incoming0f
      elseif (order == 1) then
        uacc = uacc + boxes(ib)%incoming0f + matmul(boxes(ib)%incoming1f, dx)
      elseif (order == 2) then
        uacc = uacc + boxes(ib)%incoming0f + matmul(boxes(ib)%incoming1f, dx) + &
          TijkTjk(boxes(ib)%incoming2f, dxdx)
      elseif (order == 3) then
        uacc = uacc + boxes(ib)%incoming0f + matmul(boxes(ib)%incoming1f, dx) + &
          TijkTjk(boxes(ib)%incoming2f, dxdx) + &
          TijklTjkl(boxes(ib)%incoming3f, dxdxdx)
      elseif (order == 4) then
        uacc = uacc + boxes(ib)%incoming0f + matmul(boxes(ib)%incoming1f, dx) + &
          TijkTjk(boxes(ib)%incoming2f, dxdx) + &
          TijklTjkl(boxes(ib)%incoming3f, dxdxdx) + &
          TijklmTjklm(boxes(ib)%incoming4f, dxdxdxdx)
      endif
      t2 = time()
      texp = texp + t2 - t1

      edges(iedg1)%u = - uacc
      u_edge(:,iedg1) = - uacc

    enddo
    ! !$OMP END DO
    ! !$OMP END PARALLEL

    ! loop over elements in box
    ! !$OMP PARALLEL DEFAULT(NONE) &
    ! !$OMP PRIVATE(j,iel,disgrad,i,ifc,sgn,uacc,k,iedg,u,n,l) &
    ! !$OMP SHARED(boxes,ib,elements,edges,faces,wedg,u_edge)
    ! !$OMP DO
    do j = 1,boxes(ib)%nel

      iel = boxes(ib)%element_ids(j)

      disgrad = 0.0
      do i = 1,4
        ifc = elements(iel)%face(i)
        sgn = elements(iel)%face_normal_sgn(i)

        uacc = 0.0
        do k = 1,3
          iedg = faces(ifc)%edge(k)
          ! uacc = uacc + edges(iedg)%u*faces(ifc)%area*wedg
          uacc = uacc + u_edge(:,iedg)*faces(ifc)%area*wedg
        enddo
        ! if (faces(ifc)%boundary) uacc = 0.0
        faces(ifc)%u = uacc

        u = uacc ! faces(ifc)%u
        n = faces(ifc)%normal*sgn
        do l = 1,3
          do k = 1,3
            disgrad(l,k) = disgrad(l,k) + u(l)*n(k)
          enddo
        enddo
  
        if (dot_product(- elements(iel)%x + faces(ifc)%x, n) < 0.0 .and. (.not. periodic)) & ! will not hold for boundary faces in periodic case
          stop 'not outward normal w.r.t. element'
        
      enddo
      elements(iel)%disgrad = disgrad/elements(iel)%v ! + Eapp
      ! disgrad = elements(iel)%disgrad
      if (iJacobi == 1 .and. (.not.plasticity) .and. solution_method == 1) then
        polar = mandel_t1_to_t2(elements(iel)%polar + polaravg)
        Gamma_self_polar = TijklTkl(elements(iel)%Gamma_self, polar)
        elements(iel)%disgrad = TijklTkl(elements(iel)%I_GdC_inv,elements(iel)%disgrad - Gamma_self_polar)
      endif
      elements(iel)%strain = mandel_t2_to_t1(sym33(elements(iel)%disgrad))
      ! write(*,*)norm2(disgrad - elements(iel)%disgrad)

    enddo
    ! !$OMP END DO
    ! !$OMP END PARALLEL

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! call correct_distorted(elements)

  strainavg = 0.0
  vtot = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(nel,elements) &
  !$OMP REDUCTION(+:strainavg,vtot)
  !$OMP DO 
  do iel = 1, nel

    strainavg = strainavg + elements(iel)%strain*elements(iel)%v
    vtot = vtot + elements(iel)%v
    ! write(*,*) iel, elements(iel)%strain

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  strainavg = strainavg/vtot

  ! !$OMP PARALLEL DEFAULT(NONE) &
  ! !$OMP PRIVATE(iel) &
  ! !$OMP SHARED(nel,elements,strainavg,Eapp)
  ! !$OMP DO 
  ! do iel = 1, nel
  !   elements(iel)%strain = elements(iel)%strain - strainavg + Eapp
  ! enddo
  ! !$OMP END DO
  ! !$OMP END PARALLEL

  t1 = tedg + tf1 + tf2 + texp
  ! write(*,*) tedg/t1, tf1/t1, tf2/t1, texp/t1

end

subroutine convolution_all_box(boxes, elements, faces, edges)
  use global
  use types
  use tensor_functions

  type (box_type), allocatable, intent (in) :: boxes(:)
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)

  double precision :: uacc(3), G(3,3), dx(3), wedg = 1.0/3.0, f(3), Fface(3), u(3), n(3), disgrad(3,3), sgn, dxdx(3,3), fact
  double precision :: velunited, Gamma_self_polar(3,3), dxdxdx(3,3,3),dxdxdxdx(3,3,3,3)
  double precision :: dxdx_1d(6), dxdxdx_1d(10), dxdxdxdx_1d(15), polar(3,3), Eapp6(6)
  integer :: iedg1, iedg2, i, j, ifc, ib1, ib2, iel1, iel2, iedg, l, k, iel, ii, jj, ind, kk, ib, ll, tedg, tf1, tf2, texp, t1, t2
  integer :: ind1, ind2, ind3

  ! tedg = 0
  ! tf1 = 0
  ! tf2 = 0
  ! texp = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(uacc,G,dx,f,Fface,u,n,disgrad,sgn,dxdx,fact,ind,dxdxdxdx,ll,t1,t2) &
  !$OMP PRIVATE(iedg1,iedg2,i,j,ifc,ib1,iel1,iel2,iedg,l,k,iel,ii,jj,u_edge,Gamma_self_polar,dxdxdx,ib2) &
  !$OMP PRIVATE(dxdx_1d, dxdxdx_1d, dxdxdxdx_1d,ind1,ind2,ind3,polar) &
  !$OMP SHARED(boxes,lambda,mu,edges,faces,elements,wedg,order,nbox,periodic,iJacobi,Eapp,polaravg,precalc_G) &
  !$OMP SHARED(ind2d21d,ind3d21d,ind4d21d,arr1d_2d_fact,arr1d_3d_fact,arr1d_4d_fact,plasticity,solution_method) &
  !$OMP SHARED(time_inc,Edotapp) &
  !$OMP REDUCTION(+:tedg,tf1,tf2,texp)
  !$OMP DO
  do ib1 = 1,nbox
  
    ! loop over edges in box
    ! !$OMP PARALLEL DEFAULT(NONE) &
    ! !$OMP PRIVATE(j,iedg1,uacc,i,iedg2,fact,dx,G,ifc,f,dxdx) &
    ! !$OMP SHARED(boxes,ib1,lambda,mu,edges,faces,wedg,order,u_edge)
    ! !$OMP DO
    do j = 1,boxes(ib1)%nedg ! these are all including shared

      ! t1 = time()
      iedg1 = boxes(ib1)%edge_ids(j)
      uacc = 0.0

      ! old approach using repeated edges
!      do i = 1,boxes(ib1)%nnear
!        ib2 = boxes(ib1)%near_box_id(i)
!
!        do k = 1,boxes(ib2)%nedg
!          iedg2 = boxes(ib2)%edge_ids(k)
!        
!          ! full convolution over near edges
!          if (iedg1 .ne. iedg2) then
!            dx = edges(iedg1)%xc - edges(iedg2)%xc
!            ! G = Green_function(dx, lambda, mu)
!            G = Green_function_fast(dx)
!            uacc = uacc + matmul(G, boxes(ib2)%edge_f(:,k))
!          endif
!        enddo
!    
!      enddo

      ! new approach, no repetition
      do i = 1,size(boxes(ib1)%near_edge_ids)

        iedg2 = boxes(ib1)%near_edge_ids(i)
      
        ! full convolution over near edges
        if (iedg1 .ne. iedg2) then
          dx = edges(iedg1)%xc - edges(iedg2)%xc
          ! G = Green_function(dx, lambda, mu)
          G = Green_function_fast(dx)
          uacc = uacc + matmul(G, boxes(ib1)%near_edge_f(:,i))
        endif

      enddo

      ! t2 = time()
      ! tedg = tedg + t2 - t1
  
      ! t1 = time()
      do i = 1, edges(iedg1)%nface ! loop over faces containing this edge
        ifc = edges(iedg1)%face(i)
    
        ! remove effects of faces containing this edge
        do k = 1, 3 ! loop over edges in face
          iedg2 = faces(ifc)%edge(k)
          if (iedg2 .ne. iedg1) then
            dx = edges(iedg1)%xc - edges(iedg2)%xc
            G = Green_function(dx, lambda, mu)
            f = wedg*faces(ifc)%F ! contribution to edge force from this face
            uacc = uacc - matmul(G, f)
          endif
        enddo
    
        ! add true effect of faces containing this edge
        G = edges(iedg1)%face_G(:,:,i)
        uacc = uacc + matmul(G, faces(ifc)%F)
    
      enddo
      ! t2 = time()
      ! tf1 = tf1 + t2 - t1

      ! t1 = time()
      if (allocated(edges(iedg1)%face_close)) then
        do i = 1, size(edges(iedg1)%face_close) ! loop over close faces 
          ifc = edges(iedg1)%face_close(i)
      
          ! remove effects of close faces
          do k = 1, 3 ! loop over edges in face
            iedg2 = faces(ifc)%edge(k)
            if (iedg2 .ne. iedg1) then
              dx = edges(iedg1)%xc - edges(iedg2)%xc
              G = Green_function(dx, lambda, mu)
              f = wedg*faces(ifc)%F ! contribution to edge force from this face
              uacc = uacc - matmul(G, f)
            endif
          enddo
      
          ! add true effect of faces containing this edge
          G = edges(iedg1)%face_close_G(:,:,i)
          uacc = uacc + matmul(G, faces(ifc)%F)
      
        enddo
      endif
      uacc = - uacc
      ! t2 = time()
      ! tf2 = tf2 + t2 - t1

      ! t1 = time()
      ! add effects of far field
      ! dx = edges(iedg1)%xc - boxes(ib1)%x
      ! if (periodic) dx = per_dx(dx)
      ! do ii = 1,3
      !   do jj = 1,3
      !     dxdx(ii,jj) = dx(ii)*dx(jj)
      !     do kk = 1,3
      !       dxdxdx(ii,jj,kk) = dxdx(ii,jj)*dx(kk)
      !       do ll = 1,3
      !         dxdxdxdx(ii,jj,kk,ll) = dxdxdx(ii,jj,kk)*dx(ll)
      !       enddo
      !     enddo
      !   enddo
      ! enddo
      ! dx = edges(iedg1)%dx(:,ind)
      ! dxdx = edges(iedg1)%dxdx(:,:,ind)

      if (order == 0) then
        uacc = uacc + boxes(ib1)%incoming0f_1d
      elseif (order == 1) then
        dx = edges(iedg1)%xc - boxes(ib1)%x
        uacc = uacc + boxes(ib1)%incoming0f_1d + matmul(boxes(ib1)%incoming1f_1d, dx)
      elseif (order == 2) then
        dx = edges(iedg1)%xc - boxes(ib1)%x
        do ii = 1,3
          do jj = ii,3
            ind1 = ind2d21d(ii,jj)
            dxdx_1d(ind1) = dx(ii)*dx(jj)*arr1d_2d_fact(ind1)
          enddo
        enddo
        uacc = uacc + boxes(ib1)%incoming0f_1d + matmul(boxes(ib1)%incoming1f_1d, dx) + &
          matmul(boxes(ib1)%incoming2f_1d, dxdx_1d)
      elseif (order == 3) then
        dx = edges(iedg1)%xc - boxes(ib1)%x
        do ii = 1,3
          do jj = ii,3
            ind1 = ind2d21d(ii,jj)
            dxdx_1d(ind1) = dx(ii)*dx(jj)*arr1d_2d_fact(ind1)
            do kk = jj,3
              ind2 = ind3d21d(ii,jj,kk)
              dxdxdx_1d(ind2) = dx(ii)*dx(jj)*dx(kk)*arr1d_3d_fact(ind2)
            enddo
          enddo
        enddo
        uacc = uacc + boxes(ib1)%incoming0f_1d + matmul(boxes(ib1)%incoming1f_1d, dx) + &
          matmul(boxes(ib1)%incoming2f_1d, dxdx_1d) + &
          matmul(boxes(ib1)%incoming3f_1d, dxdxdx_1d)
      elseif (order == 4) then
        dx = edges(iedg1)%xc - boxes(ib1)%x
        do ii = 1,3
          do jj = ii,3
            ind1 = ind2d21d(ii,jj)
            dxdx_1d(ind1) = dx(ii)*dx(jj)*arr1d_2d_fact(ind1)
            do kk = jj,3
              ind2 = ind3d21d(ii,jj,kk)
              dxdxdx_1d(ind2) = dx(ii)*dx(jj)*dx(kk)*arr1d_3d_fact(ind2)
              do ll = kk,3
                ind3 = ind4d21d(ii,jj,kk,ll)
                dxdxdxdx_1d(ind3) = dx(ii)*dx(jj)*dx(kk)*dx(ll)*arr1d_4d_fact(ind3)
              enddo
            enddo
          enddo
        enddo
        uacc = uacc + boxes(ib1)%incoming0f_1d + matmul(boxes(ib1)%incoming1f_1d, dx) + &
          matmul(boxes(ib1)%incoming2f_1d, dxdx_1d) + &
          matmul(boxes(ib1)%incoming3f_1d, dxdxdx_1d) + &
          matmul(boxes(ib1)%incoming4f_1d, dxdxdxdx_1d)
      endif
      ! t2 = time()
      ! texp = texp + t2 - t1

      edges(iedg1)%u = uacc
      u_edge(:,iedg1) = uacc
      ! write(*,*)iedg1, uacc
      ! stop

    enddo
    ! !$OMP END DO
    ! !$OMP END PARALLEL

    ! loop over elements in box
    ! !$OMP PARALLEL DEFAULT(NONE) &
    ! !$OMP PRIVATE(j,iel,disgrad,i,ifc,sgn,uacc,k,iedg,u,n,l) &
    ! !$OMP SHARED(boxes,ib,elements,edges,faces,wedg,u_edge)
    ! !$OMP DO
    do j = 1,boxes(ib1)%nel

      iel = boxes(ib1)%element_ids(j)

      disgrad = 0.0
      do i = 1,4
        ifc = elements(iel)%face(i)
        sgn = elements(iel)%face_normal_sgn(i)

        uacc = 0.0
        do k = 1,3
          iedg = faces(ifc)%edge(k)
          ! uacc = uacc + edges(iedg)%u*faces(ifc)%area*wedg
          uacc = uacc + u_edge(:,iedg)*faces(ifc)%area*wedg
        enddo
        ! if (faces(ifc)%boundary) uacc = 0.0
        faces(ifc)%u = uacc

        u = uacc ! faces(ifc)%u
        n = faces(ifc)%normal*sgn
        do l = 1,3
          do k = 1,3
            disgrad(l,k) = disgrad(l,k) + u(l)*n(k)
          enddo
        enddo
  
        if (dot_product(- elements(iel)%x + faces(ifc)%x, n) < 0.0 .and. (.not. periodic)) & ! will not hold for boundary faces in periodic case
          stop 'not outward normal w.r.t. element'
        
      enddo
      elements(iel)%disgrad = elements(iel)%disgradt + disgrad/elements(iel)%v + Edotapp*time_inc
      ! disgrad = elements(iel)%disgrad
      if (iJacobi == 1 .and. (.not.plasticity) .and. solution_method == 1) then
        polar = mandel_t1_to_t2(elements(iel)%polar + polaravg)
        Gamma_self_polar = TijklTkl(elements(iel)%Gamma_self, polar)
        elements(iel)%disgrad = TijklTkl(elements(iel)%I_GdC_inv,elements(iel)%disgrad - Gamma_self_polar)
      endif
      elements(iel)%strain = mandel_t2_to_t1(sym33(elements(iel)%disgrad))
      ! write(*,*)norm2(disgrad - elements(iel)%disgrad)

    enddo
    ! !$OMP END DO
    ! !$OMP END PARALLEL

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! call correct_distorted(elements)

  strainavg = 0.0
  vtot = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(nel,elements) &
  !$OMP REDUCTION(+:strainavg,vtot)
  !$OMP DO 
  do iel = 1, nel

    strainavg = strainavg + elements(iel)%strain*elements(iel)%v
    vtot = vtot + elements(iel)%v
    ! write(*,*) iel, elements(iel)%strain

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  strainavg = strainavg/vtot

  Eapp6 = mandel_t2_to_t1(Eapp)
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(nel,elements,strainavg,Eapp6)
  !$OMP DO 
  do iel = 1, nel
    elements(iel)%strain = elements(iel)%strain - strainavg + Eapp6
    elements(iel)%disgrad = elements(iel)%disgrad + mandel_t1_to_t2(- strainavg + Eapp6)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ! t1 = tedg + tf1 + tf2 + texp
  ! write(*,*) tedg, tf1, tf2, texp

end

subroutine element_aspect_ratio(elements, faces)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)

  integer :: iel, info, inode, i, j, ifc, iel_attach, iel_neigh, iel_neigh_arr(4), ic, grain_id_arr(4)
  double precision :: work(64), eig_dir(3,3), eval(3), dx(3), face_area(4), eval33(3,3), dxdx(3,3)
  double precision, allocatable :: aspect(:)

  allocate(aspect(nel))
  do iel = 1,nel

    dxdx = 0.0
    do inode = 1,4
      dx = elements(iel)%xnode(:,inode) - elements(iel)%x
      do i = 1,3
        do j = 1,3
          dxdx(i,j) = dxdx(i,j) + dx(i)*dx(j)
        enddo
      enddo
    enddo
    dxdx = dxdx/4.0
    eig_dir = dxdx
    call dsyev('V','U',3,eig_dir,3,eval,work,64,info) ! LAPACK Symmetric eigenvalues & vectors

    ! ! test
    ! eval33 = 0.0
    ! eval33(1,1) = eval(1)
    ! eval33(2,2) = eval(2)
    ! eval33(3,3) = eval(3)
    ! write(*,*) norm2(matmul(matmul(eig_dir,eval33),transpose(eig_dir)) - dxdx)
    ! write(*,*) maxloc(eval)
    
    aspect(iel) = maxval(eval)/minval(eval)
    allocate(elements(iel)%eval(3))
    allocate(elements(iel)%eig_dir(3,3))
    elements(iel)%eval = eval
    elements(iel)%eig_dir = eig_dir
    elements(iel)%aspect = aspect(iel)

    if (aspect(iel) > aspect_mx) then

      write(*,*) iel, aspect(iel)

    endif
  
  enddo
! 
!   allocate(el_distorted(nel))
!   nel_distorted = 0
!   do iel = 1,nel
!     if (aspect(iel) > aspect_mx) then
!       
!       nel_distorted = nel_distorted + 1
!       el_distorted(nel_distorted) = iel
!       ic = 0
!       do i = 1,4 ! loop over faces
!         ifc = elements(iel)%face(i)
!         if (.not.faces(ifc)%boundary) then ! not a boundary face
! 
!           if (faces(ifc)%el(1) == iel) then
!             iel_neigh = faces(ifc)%el(2) ! neighboring element
!           elseif (faces(ifc)%el(2) == iel) then
!             iel_neigh = faces(ifc)%el(1)
!           else
!             stop 'error in faces(ifc)%el'
!           endif
! 
!           ! if (aspect(iel_neigh) < aspect_mx .and. elements(iel_neigh)%grain_id == elements(iel)%grain_id) then ! only if neighbor is not distorted
!           if (aspect(iel_neigh) < aspect_mx) then ! only if neighbor is not distorted
!             ic = ic + 1
!             face_area(ic) = faces(ifc)%area
!             iel_neigh_arr(ic) = iel_neigh
!             grain_id_arr(ic) = elements(iel_neigh)%grain_id
!           endif
! 
!         endif
! 
!       enddo
!       if (ic == 0) stop 'no good elements to use'
! 
! !      if (ic > 1) then
! !
! !        i = maxloc(face_area, 1, mask = grain_id_arr == elements(iel)%grain_id)
! !        write(*,*) iel
! !        write(*,*) elements(iel)%grain_id
! !        write(*,*) face_area
! !        write(*,*) grain_id_arr
! !        write(*,*) i
! !        read(*,*)       
! !        if (i > 0) then ! force only one
! !          ic = 1
! !          face_area(1) = face_area(i)
! !          iel_neigh_arr(ic) = iel_neigh_arr(i)
! !          grain_id_arr(ic) = grain_id_arr(i)
! !        endif
! !
! !      endif
! 
!       allocate(elements(iel)%iel_interp(ic))
!       elements(iel)%iel_interp = iel_neigh_arr(1:ic)
! 
!     endif
!   enddo

  deallocate(aspect)
  ! stop

  ! ! test
  ! iel = 6837
  ! nel_distorted = 1
  ! el_distorted(1) = iel
  ! ic = 1
  ! allocate(elements(iel)%iel_interp(ic))
  ! elements(iel)%iel_interp = 6702

end

subroutine correct_distorted(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: i, j, k, iel1, iel2
  double precision :: velunited, evalm, v(3), disgrad_v(3)

  do i = 1,nel_distorted
    iel1 = el_distorted(i)

    elements(iel1)%disgrad = elements(iel1)%disgrad*elements(iel1)%v
    velunited = elements(iel1)%v
    do j = 1,size(elements(iel1)%iel_interp)
      iel2 = elements(iel1)%iel_interp(j)
      elements(iel1)%disgrad = elements(iel1)%disgrad + elements(iel2)%disgrad*elements(iel2)%v
      velunited = velunited + elements(iel2)%v
    enddo
    elements(iel1)%disgrad = elements(iel1)%disgrad/velunited
    elements(iel1)%strain = mandel_t2_to_t1(sym33(elements(iel1)%disgrad))
    do j = 1,size(elements(iel1)%iel_interp)
      iel2 = elements(iel1)%iel_interp(j)
      elements(iel2)%disgrad = elements(iel1)%disgrad
      elements(iel2)%strain = elements(iel1)%strain
    enddo

!     write(*,*) iel1, norm2(elements(iel1)%disgrad)
!     v = elements(iel1)%eig_dir(:,1)
!     disgrad_v = matmul(elements(iel1)%disgrad, v)
!     write(*,*) matmul(elements(iel1)%disgrad, elements(iel1)%eig_dir(:,1))
!     write(*,*) matmul(elements(iel1)%disgrad, elements(iel1)%eig_dir(:,2))
!     write(*,*) matmul(elements(iel1)%disgrad, elements(iel1)%eig_dir(:,3))
!     do j = 1,3
!       do k = 1,3
!         elements(iel1)%disgrad(j,k) = elements(iel1)%disgrad(j,k) - disgrad_v(j)*v(k)
!       enddo
!     enddo
!     !  write(*,*) iel1, norm2(elements(iel1)%disgrad)
! 
! 
!     ! if (elements(iel1)%eval(3)/elements(iel1)%eval(2) > aspect_mx) then
!     !   v = elements(iel1)%eig_dir(:,2)
!     !   disgrad_v = matmul(elements(iel1)%disgrad, v)
!     !   do j = 1,3
!     !     do k = 1,3
!     !       elements(iel1)%disgrad(j,k) = elements(iel1)%disgrad(j,k) - disgrad_v(j)*v(k)
!     !     enddo
!     !   enddo
!     ! endif

  enddo

end

function per_dx(dx) result(dx_per)
  use global, only : npts1, npts2, npts3
  double precision, intent(in) :: dx(3)
  double precision :: dx_per(3)

  dx_per = dx
  if (abs(dx_per(1) - float(npts1)) < abs(dx_per(1))) then 
    dx_per(1) = dx_per(1) - float(npts1)
  elseif (abs(dx_per(1) + float(npts1)) < abs(dx_per(1))) then 
    dx_per(1) = dx_per(1) + float(npts1)
  endif
  if (abs(dx_per(2) - float(npts2)) < abs(dx_per(2))) then
    dx_per(2) = dx_per(2) - float(npts2)
  elseif (abs(dx_per(2) + float(npts2)) < abs(dx_per(2))) then
    dx_per(2) = dx_per(2) + float(npts2)
  endif
  if (abs(dx_per(3) - float(npts3)) < abs(dx_per(3))) then
    dx_per(3) = dx_per(3) - float(npts3)
  elseif (abs(dx_per(3) + float(npts3)) < abs(dx_per(3))) then
    dx_per(3) = dx_per(3) + float(npts3)
  endif
  
end


subroutine create_faces_per(elements, faces)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (out) :: faces(:)
  type (face_type), allocatable :: all_faces(:)

  integer :: in1, in2, in3, iel, ic, i, ifc, ind1, ind2, ifc1, ifc2
  integer, allocatable :: unique_faces_ind(:,:),ind_duplicate(:)
  double precision :: tol = 1.0e-5, x(3), dx1(3), dx2(3), dx(3), dir_range(3), x1(3), x2(3), x3(3)
  double precision :: xnode_face(3,3)
  logical :: new, flag

  allocate(all_faces(nel*4))
  allocate(unique_faces_ind(2,nel*4))

  el_face_nodes(:,1) = [1,2,3]
  el_face_nodes(:,2) = [1,2,4]
  el_face_nodes(:,3) = [1,3,4]
  el_face_nodes(:,4) = [2,3,4]

  dir_range(1) = float(npts1)
  dir_range(2) = float(npts2)
  dir_range(3) = float(npts3)

  write(*,*) '   all faces'
  ifc = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,ifc,x1,x2,x3,i) &
  !$OMP SHARED(nel,elements,el_face_nodes,all_faces,tol,dir_range)
  !$OMP DO 
  do iel = 1,nel
    ifc = (iel-1)*4 + 1
    x1 = elements(iel)%xnode(:,el_face_nodes(1,1))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,1))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,1))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 1
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        ! all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 2
    x1 = elements(iel)%xnode(:,el_face_nodes(1,2))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,2))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,2))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 2
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        ! all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 3
    x1 = elements(iel)%xnode(:,el_face_nodes(1,3))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,3))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,3))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 3
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        ! all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif

    ifc = (iel-1)*4 + 4
    x1 = elements(iel)%xnode(:,el_face_nodes(1,4))
    x2 = elements(iel)%xnode(:,el_face_nodes(2,4))
    x3 = elements(iel)%xnode(:,el_face_nodes(3,4))
    all_faces(ifc)%x = (x1 + x2 + x3)/3.0
    all_faces(ifc)%el(1) = iel
    all_faces(ifc)%el_face(1) = 4
    all_faces(ifc)%nel = 1
    all_faces(ifc)%boundary = .false.
    do i = 1,3
      if ( (abs(x1(i) - 0.5) <= tol .and. abs(x2(i) - 0.5) <= tol .and. abs(x3(i) - 0.5) <= tol) .or. &
           (abs(x1(i) - (dir_range(i)+0.5)) <= tol .and. abs(x2(i) - (dir_range(i)+0.5)) <= tol .and. &
            abs(x3(i) - (dir_range(i)+0.5)) <= tol) ) then
        ! all_faces(ifc)%boundary = .true.
      endif
    enddo
    all_faces(ifc)%normal = cross_prod(x2 - x1, x3 - x1)
    all_faces(ifc)%normal = all_faces(ifc)%normal/norm2(all_faces(ifc)%normal)
    if (dot_product(- elements(iel)%x + all_faces(ifc)%x, all_faces(ifc)%normal) < 0.0) then
      all_faces(ifc)%normal = - all_faces(ifc)%normal ! make sure it is outward normal
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
!  ic = 1
!  unique_faces_ind = 0
!  unique_faces_ind(1,ic) = 1
!  do ifc = 2,nel*4
!
!    x = all_faces(ifc)%x
!    new = .true.
!    if (.not.all_faces(ifc)%boundary) then ! if this is boundary face treat it as new
!      do i = 1,ic
!        if (unique_faces_ind(2,i) == 0) then ! check if it is the same only for faces that do not have duplicate assigned
!          ind1 = unique_faces_ind(1,i)
!          ! if (norm2(x - all_faces(ind1)%x) <= tol) then 
!          dx = abs(x - all_faces(ind1)%x)
!          if (.not.(all_faces(ind1)%boundary).and. &
!              (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol)) then ! duplicate, skip comparing with stored boundary faces
!            unique_faces_ind(2,i) = ifc
!            new = .false.
!            exit
!          endif
!        endif
!      enddo
!    endif
!
!    if (new) then
!
!      ic = ic + 1
!      unique_faces_ind(1,ic) = ifc
!
!    endif
!
!  enddo

  write(*,*) '   duplicates'
  allocate(ind_duplicate(nel*4))
  ind_duplicate = 0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc1,x,ifc2,dx) &
  !$OMP SHARED(nel,all_faces,ind_duplicate,tol,flag)
  !$OMP DO 
  do ifc1 = 1,nel*4
    x = all_faces(ifc1)%x
    ! if (.not.all_faces(ifc1)%boundary) then ! if this is boundary face treat it as new
      do ifc2 = 1,nel*4
        ! dx = abs(x - all_faces(ifc2)%x)
        dx = x - all_faces(ifc2)%x
        dx = abs(per_dx(dx))
        flag = .true.
        if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol .and. ifc1 .ne. ifc2) then ! duplicate, skip comparing with stored boundary faces
          ind_duplicate(ifc1) = ifc2
          flag = .false.
          exit
        endif
      enddo
      if (flag) write(*,*) ifc1, all_faces(ifc1)%el(1), x
      if (flag) stop
    ! endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   unique faces'
  ic = 0
  unique_faces_ind = 0
  do ifc = 1,nel*4
    if (mod(ifc,10000) == 0) write(*,*) float(ifc)/float(nel*4)
    if (ind_duplicate(ifc) == 0) then
      ic = ic + 1
      unique_faces_ind(1,ic) = ifc
    else
      if (ind_duplicate(ifc) > ifc) then
        ic = ic + 1
        unique_faces_ind(1,ic) = ifc
        unique_faces_ind(2,ic) = ind_duplicate(ifc)
      endif
    endif
  enddo

  deallocate(ind_duplicate)

  write(*,*) '   create faces'
  allocate(faces(ic))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,ind1,iel,x1,x2,x3,xnode_face,i,ind2) &
  !$OMP SHARED(ic,faces,all_faces,elements,unique_faces_ind)
  !$OMP DO 
  do ifc = 1,ic
    ind1 = unique_faces_ind(1,ifc)
    faces(ifc)%x = all_faces(ind1)%x
    faces(ifc)%el(1) = all_faces(ind1)%el(1)
    faces(ifc)%el_face(1) = all_faces(ind1)%el_face(1)
    faces(ifc)%boundary = all_faces(ind1)%boundary
    faces(ifc)%normal = all_faces(ind1)%normal

    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
    xnode_face(:,1) = x1
    xnode_face(:,2) = x2
    xnode_face(:,3) = x3
    faces(ifc)%area = triangle_area(xnode_face)

    iel = faces(ifc)%el(1)
    i = faces(ifc)%el_face(1)
    if (elements(iel)%face(i) .ne. 0) stop 'face already assignd to this element place'
    elements(iel)%face(i) = ifc
    elements(iel)%face_normal_sgn(i) = 1.0 ! first is positive

    ind2 = unique_faces_ind(2,ifc)
    if (ind2 .ne. 0) then
      faces(ifc)%el(2) = all_faces(ind2)%el(1)
      faces(ifc)%el_face(2) = all_faces(ind2)%el_face(1)
      faces(ifc)%nel = 2

      iel = faces(ifc)%el(2)
      i = faces(ifc)%el_face(2)
      if (elements(iel)%face(i) .ne. 0) stop 'face already assignd to this element place'
      elements(iel)%face(i) = ifc
      elements(iel)%face_normal_sgn(i) = -1.0 ! second is negative

    else
      faces(ifc)%el(2) = 0
      faces(ifc)%nel = 1
      ! write(*,*) ifc, ind1
    endif

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  nfc = ic
  write(*,*)ic, nel

  ! verify
  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,i,iel,x,dx) &
  !$OMP SHARED(ic,faces,tol,elements,el_face_nodes)
  !$OMP DO 
  do ifc = 1,ic
    if (faces(ifc)%nel == 2) then

      i = faces(ifc)%el_face(1)
      iel = faces(ifc)%el(1)
      x = (elements(iel)%xnode(:,el_face_nodes(1,i)) + elements(iel)%xnode(:,el_face_nodes(2,i)) + &
       elements(iel)%xnode(:,el_face_nodes(3,i)))/3.0

      ! call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
      ! if (norm2( (x1 + x2 + x3)/3.0 - x) > tol) then
      !   stop 'problem'
      ! endif

      dx = abs(x - faces(ifc)%x)
      if (.not.((dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol)) ) then ! duplicate
        stop 'error in face'
      endif

    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,i) &
  !$OMP SHARED(nel,elements)
  !$OMP DO 
  do iel = 1,nel
    do i = 1,4
      if (elements(iel)%face(i) == 0) stop 'element without all faces'
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(all_faces)
  deallocate(unique_faces_ind)

  ! stop

end


subroutine create_edges_per(elements, faces, edges)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (edge_type), allocatable, intent (out) :: edges(:)
  type (edge_type), allocatable :: all_edges(:)

  integer :: iel, ifc, ind1, ind2, iedg, ic, i, j, ind, iedg1, iedg2
  integer, allocatable :: unique_edges_ind(:,:), nduplicates(:), ind_duplicate(:,:), nduplicates_all(:)
  double precision :: tol = 1.0e-5, x(3), dx1(3), dx2(3), dx(3),  x1(3), x2(3), x3(3)
  double precision :: xnode_face(3,3)
  logical :: new

  allocate(all_edges(nfc*3))
  allocate(unique_edges_ind(nmx_shared_faces_per_edge,nfc*3))
  allocate(nduplicates(nfc*3))

  face_edge_nodes(:,1) = [1,2]
  face_edge_nodes(:,2) = [2,3]
  face_edge_nodes(:,3) = [3,1]

  face_edge_nodes_remaining(1) = 3
  face_edge_nodes_remaining(2) = 1
  face_edge_nodes_remaining(3) = 2

  write(*,*) '   all edges'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,iedg,iel,xnode_face,x1,x2) &
  !$OMP SHARED(nfc,faces,elements,face_edge_nodes,all_edges)
  !$OMP DO 
  do ifc = 1,nfc

    iedg = (ifc-1)*3 + 1
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,1))
    x2 = xnode_face(:,face_edge_nodes(2,1))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 1

    iedg = (ifc-1)*3 + 2
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,2))
    x2 = xnode_face(:,face_edge_nodes(2,2))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 2

    iedg = (ifc-1)*3 + 3
    iel = faces(ifc)%el(1)
    call get_face_nodes(faces(ifc),elements(iel),xnode_face(:,1), xnode_face(:,2), xnode_face(:,3))
    x1 = xnode_face(:,face_edge_nodes(1,3))
    x2 = xnode_face(:,face_edge_nodes(2,3))
    all_edges(iedg)%x(:,1) = x1
    all_edges(iedg)%x(:,2) = x2
    all_edges(iedg)%xc = (x1 + x2)*0.5
    allocate(all_edges(iedg)%face(1))
    allocate(all_edges(iedg)%face_edge(1))
    all_edges(iedg)%face(1) = ifc
    all_edges(iedg)%face_edge(1) = 3

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

!  ic = 1
!  unique_edges_ind = 0
!  nduplicates = 0
!  unique_edges_ind(1,ic) = 1
!  nduplicates(ic) = 1
!  do iedg = 2,nfc*3
!
!    x = all_edges(iedg)%xc
!    new = .true.
!    do i = 1,ic
!      do j = 1,nduplicates(i)
!        ind1 = unique_edges_ind(j,i)
!        dx = abs(x - all_edges(ind1)%xc)
!        if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol) then ! duplicate
!          nduplicates(i) = nduplicates(i) + 1
!          if (nduplicates(i) > nmx_shared_faces_per_edge) stop 'nduplicates(i) > nmx_shared_faces_per_edge'
!          unique_edges_ind(nduplicates(i),i) = iedg
!          new = .false.
!          exit
!        endif
!      enddo
!    enddo
!
!    if (new) then
!
!      ic = ic + 1
!      unique_edges_ind(1,ic) = iedg
!      nduplicates(ic) = 1
!
!    endif
!
!  enddo

  write(*,*) '   duplicates'
  allocate(ind_duplicate(nmx_shared_faces_per_edge,nfc*3))
  allocate(nduplicates_all(nfc*3))
  ind_duplicate = 0
  nduplicates_all = 0
  write(*,*)nfc*3
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,x,ic,iedg2,dx) &
  !$OMP SHARED(all_edges,nfc,ind_duplicate,nduplicates_all,tol)
  !$OMP DO 
  do iedg1 = 1,nfc*3
    x = all_edges(iedg1)%xc
    ic = 0
    do iedg2 = 1,nfc*3
      ! dx = abs(x - all_edges(iedg2)%xc)
      dx = x - all_edges(iedg2)%xc
      dx = abs(per_dx(dx))
      if (dx(1) <= tol .and. dx(2) <= tol .and. dx(3) <= tol) then ! duplicate
        ic = ic + 1
        ind_duplicate(ic,iedg1) = iedg2
        nduplicates_all(iedg1) = nduplicates_all(iedg1) + 1
      endif
    enddo
    ! write(*,*)iedg1, ind_duplicate(:,iedg1)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  write(*,*) '   unique edges'
  unique_edges_ind = 0
  nduplicates = 0
  ic = 0
  do iedg1 = 1,nfc*3

    if (mod(iedg1,10000) == 0) write(*,*) float(iedg1)/float(nfc*3)

    i = minval(ind_duplicate(:,iedg1), mask = ind_duplicate(:,iedg1) > 0)
    if (iedg1 == i) then
      ic = ic + 1
      unique_edges_ind(:,ic) = ind_duplicate(:,iedg1)
      nduplicates(ic) = nduplicates_all(iedg1)
    endif

  enddo

  deallocate(ind_duplicate)
  deallocate(nduplicates_all)
  ! stop

  ! do iedg = 1,ic
  !   write(*,*) iedg, nduplicates(iedg), unique_edges_ind(1:nduplicates(iedg),iedg)
  !   do i = 1,nduplicates(iedg)
  !     write(*,*) all_edges(unique_edges_ind(i,iedg))%xc
  !   enddo
  ! enddo

  write(*,*) '   create edges'
  allocate(edges(ic))
  allocate(u_edge(3,ic))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,ind1,i,ind,ifc) &
  !$OMP SHARED(ic,unique_edges_ind,edges,all_edges,nduplicates,faces)
  !$OMP DO 
  do iedg = 1,ic
    ind1 = unique_edges_ind(1,iedg)
    
    edges(iedg)%x = all_edges(ind1)%x
    edges(iedg)%xc = all_edges(ind1)%xc
    edges(iedg)%nface = nduplicates(iedg)

    allocate(edges(iedg)%face(edges(iedg)%nface))
    allocate(edges(iedg)%face_edge(edges(iedg)%nface))
    allocate(edges(iedg)%face_G(3,3,edges(iedg)%nface))

    do i = 1,nduplicates(iedg)
      ind = unique_edges_ind(i,iedg)
      edges(iedg)%face(i) = all_edges(ind)%face(1)
      edges(iedg)%face_edge(i) = all_edges(ind)%face_edge(1)
    enddo

    do i = 1,edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      if (faces(ifc)%edge(edges(iedg)%face_edge(i)) == 0) then
        faces(ifc)%edge(edges(iedg)%face_edge(i)) = iedg
      else
        stop 'edge already placed in this face'
      endif
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  nedg = ic

  write(*,*) '   close neighboring edges'
  ! close neighboring edges (belong to faces containing an edge)
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg,ic,i,ifc,j) &
  !$OMP SHARED(nedg,edges,faces)
  !$OMP DO 
  do iedg = 1,nedg

    allocate(edges(iedg)%neigh_edges(edges(iedg)%nface*2)) ! 2 times number of faces this edge belongs to
    ic = 0
    do i = 1,edges(iedg)%nface
      ifc = edges(iedg)%face(i)
      do j = 1,3
        if (faces(ifc)%edge(j) .ne. iedg) then
          ic = ic + 1
          edges(iedg)%neigh_edges(ic) = faces(ifc)%edge(j)
        endif
      enddo
    enddo

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

!  iedg = 77
!  write(*,*) edges(iedg)%x(:,1)
!  write(*,*) edges(iedg)%x(:,2)
!  do i = 1,size(edges(iedg)%neigh_edges)
!    j = edges(iedg)%neigh_edges(i)
!    write(*,*) j
!    write(*,*) edges(j)%x(:,1)
!    write(*,*) edges(j)%x(:,2) ! must share one node with other
!  enddo

!  do iedg = 1,nedg
!    
!    write(*,*) iedg
!    do i = 1,edges(iedg)%nface
!      ifc = edges(iedg)%face(i)
!      iel = faces(ifc)%el(1)
!      call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
!      write(*,*) norm2(x1), norm2(x2), norm2(x3)
!    enddo
!
!  enddo

  write(*,*) '   verify'
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ifc,i) &
  !$OMP SHARED(nfc,faces)
  !$OMP DO 
  do ifc = 1,nfc
    do i = 1,3
      if (faces(ifc)%edge(i) == 0) then
        stop 'faces(ifc)%edge(i) == 0'
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(all_edges)
  deallocate(unique_edges_ind)
  deallocate(nduplicates)

end


function Green_function_per_interp(x) result(Gper_interp)
  use global, only : Gper, npts1, npts2, npts3
  double precision, intent(in) :: x(3)
  double precision :: Gper_interp(3,3), xi(3), phi(8)
  integer :: ind8(8,3), i

  xi = get_xi(x)
  ind8 = get_grid_ind(x, npts1, npts2, npts3)
  phi = shape_funct8(xi)

  Gper_interp = 0.0
  do i = 1,8
    Gper_interp = Gper_interp + Gper(:,:,ind8(i,1),ind8(i,2),ind8(i,3))*phi(i)
  enddo

end

function Green_function_per_round(dx) result(G)
  use tensor_functions
  use global, only : Gper, npts1, npts2, npts3

  double precision, intent(in) :: dx(3)
  double precision :: G(3,3)
  integer :: i, j, k

  i = enforce_periodic(nint(dx(1)) + 1, 1, npts1)
  j = enforce_periodic(nint(dx(2)) + 1, 1, npts2)
  k = enforce_periodic(nint(dx(3)) + 1, 1, npts3)
  G = Gper(:,:,i,j,k) 

end

function Green_function_diff_round(dx, lambda, mu) result(dG)
  use tensor_functions
  use global, only : Gper, npts1, npts2, npts3

  double precision, intent(in) :: dx(3), lambda, mu
  double precision :: G(3,3), Gper_loc(3,3), dG(3,3), dxper(3), dx_round(3)
  integer :: i, j, k

  dx_round(1) = float(nint(dx(1)))
  dx_round(2) = float(nint(dx(2)))
  dx_round(3) = float(nint(dx(3)))
  if (norm2(dx) > 0.0) then
    dxper = per_dx(dx)
    G = Green_function(dxper, lambda, mu)
  else
    G = 0.0
    ! G = Green_function_per_round(dx)
  endif
  Gper_loc = Green_function_per_round(dx)
  dG = Gper_loc - G

end

subroutine self_interaction_Gamma(elements, faces, edges, intg_pt)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (in) :: edges(:)
  type (intg_pt_type), intent (in) :: intg_pt(nintg_pt)

  integer :: iel, i, j, k, l, ifc1, ifc2, iedg1, iedg2, ind, i1, i2, i3, i4, iel2, side, nproj
  double precision :: G_faceIJ(3,3,4,4), G(3,3), wedg = 1.0/3.0, dx(3), n1(3), n2(3), sgn1, sgn2
  double precision :: Gamma_self(3,3,3,3), x1(3), x2(3), x3(3), xface(3), dist_node(3), xproj(3)

  do iel = 1,nel
  ! write(*,*) 'iel', iel
  elements(iel)%Gamma_self = 0.0
  G_faceIJ = 0.0
  do i = 1,4
    ifc1 = elements(iel)%face(i)
    do j = 1,4
      ifc2 = elements(iel)%face(j)

      G = 0.0
      if (ifc1 == ifc2) then ! self interaction
        do k = 1,3
          iedg1 = faces(ifc1)%edge(k)
          ind = findloc(edges(iedg1)%face, ifc1, 1)
          if (ind == 0) then
            write(*,*) 'ind=0'
            write(*,*) edges(iedg1)%face, ifc1
          endif
          G = G - edges(iedg1)%face_G(:,:,ind)*faces(ifc1)%area*(wedg*faces(ifc1)%area)
        enddo
      else
        do k = 1,3
          iedg1 = faces(ifc1)%edge(k)
          ind = findloc(edges(iedg1)%face, ifc2, 1) 
          if (ind == 0) then ! this edge (iedg1) does not belong to ifc2

            if (allocated(edges(iedg1)%face_close)) then
              ind = findloc(edges(iedg1)%face_close, ifc2, 1) 
            else
              ind = 0
            endif
            if (ind == 0) then ! this face (ifc2) does not belong to close faces of iedg1, regular second oreder integration 
              do l = 1,3
                iedg2 = faces(ifc2)%edge(l)
                dx = edges(iedg1)%xc - edges(iedg2)%xc
                G = G - Green_function(dx, lambda, mu)*wedg*faces(ifc1)%area*wedg*faces(ifc2)%area
              enddo
            else ! this face belongs to iedg1 close faces, use the integrated expression
              G = G - edges(iedg1)%face_close_G(:,:,ind)*faces(ifc2)%area*(wedg*faces(ifc1)%area)
            endif
          else ! this edge belongs to ifc2, use the integrated expression
            G = G - edges(iedg1)%face_G(:,:,ind)*faces(ifc2)%area*(wedg*faces(ifc1)%area)
          endif
        enddo
      endif
      G_faceIJ(:,:,i,j) = G

      sgn1 = elements(iel)%face_normal_sgn(i)
      n1 = faces(ifc1)%normal*sgn1
      sgn2 = elements(iel)%face_normal_sgn(j)
      n2 = faces(ifc2)%normal*sgn2

      do i1 = 1,3
        do i2 = 1,3
          do i3 = 1,3
            do i4 = 1,3
              elements(iel)%Gamma_self(i1,i2,i3,i4) = elements(iel)%Gamma_self(i1,i2,i3,i4) + &
               G(i1,i3)*n1(i2)*n2(i4)/elements(iel)%V
            enddo
          enddo
        enddo
      enddo

    enddo
  enddo

!  ! using different integration points
!  Gamma_self = 0.0
!  do i = 1,4
!    ifc1 = elements(iel)%face(i)
!    do j = 1,4
!      ifc2 = elements(iel)%face(j)
!
!      ! decide how to integrate
!      G = 0.0
!      if (ifc1 == ifc2) then ! self-interaction of a triangle
!        do k = 1,3
!          iedg1 = faces(ifc1)%edge(k)
!          ind = findloc(edges(iedg1)%face, ifc1, 1)
!          G = G - edges(iedg1)%face_G(:,:,ind)*faces(ifc1)%area*(wedg*faces(ifc1)%area)
!        enddo
!      else
!
!        iel2 = faces(ifc2)%el(1)
!        call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
!
!        do k = 1,3 ! loop over edges of ifc1 
!          iedg1 = faces(ifc1)%edge(k)
!          
!          xface = point2plane(faces(ifc2)%x, faces(ifc2)%normal, edges(iedg1)%xc) ! project edge center to face ifc2
!          
!          dist_node = [norm2(xface - x1), norm2(xface - x2), norm2(xface - x3)] ! distance to face nodes
!          if (minval(dist_node) <= min_dist_node) then ! edge center projection close to one node
!
!            ind = minloc(dist_node,1)
!            if (ind == 1) then
!              G = G - integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!            elseif (ind == 2) then
!              G = G - integrate_Green_function(x2, x3, x1, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!            elseif (ind == 3) then
!              G = G - integrate_Green_function(x3, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!            endif
!
!          else ! edge center projection not close to node
!
!            if (point_in_triangle(x1, x2, x3, faces(ifc2)%area, xface)) then ! edge center projection in triangle
!
!              call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj)
!
!              if (nproj == 0) then! projection too close to nodes, use closest node
!
!                ! use closest node
!                dist_node = [norm2(xface - x1), norm2(xface - x2), norm2(xface - x3)] ! distance to face nodes
!                ind = minloc(dist_node,1)
!                if (ind == 1) then
!                  G = G - integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (ind == 2) then
!                  G = G - integrate_Green_function(x2, x3, x1, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (ind == 3) then
!                  G = G - integrate_Green_function(x3, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                endif
!
!              else
!                if (norm2(xproj - xface) < min_dist_node) then ! almost on side (use two triangles)
!
!                  if (side == 1) then ! on side x1-x2
!                    G = G - integrate_Green_function(xproj, x1, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                    G = G - integrate_Green_function(xproj, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  elseif (side == 2) then ! on side x2-x3
!                    G = G - integrate_Green_function(xproj, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                    G = G - integrate_Green_function(xproj, x1, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  elseif (side == 3) then ! on side x3-x1
!                    G = G - integrate_Green_function(xproj, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                    G = G - integrate_Green_function(xproj, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  endif
!  
!                else ! inside, use three triangles
!   
!                  G = G - integrate_Green_function(xface, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  G = G - integrate_Green_function(xface, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  G = G - integrate_Green_function(xface, x3, x1, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!  
!                endif
!              endif
!
!            else ! point is outside
!
!              call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj)
!
!              if (nproj > 0) then ! project to side and use two triangles
!
!                if (side == 1) then ! on side x1-x2
!                  G = G - integrate_Green_function(xproj, x1, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  G = G - integrate_Green_function(xproj, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (side == 2) then ! on side x2-x3
!                  G = G - integrate_Green_function(xproj, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  G = G - integrate_Green_function(xproj, x1, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (side == 3) then ! on side x3-x1
!                  G = G - integrate_Green_function(xproj, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                  G = G - integrate_Green_function(xproj, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                endif
!
!              elseif (nproj == 0) then ! projection too close to nodes, use closest node
!
!                ! use closest node
!                dist_node = [norm2(xface - x1), norm2(xface - x2), norm2(xface - x3)] ! distance to face nodes
!                ind = minloc(dist_node,1)
!                if (ind == 1) then
!                  G = G - integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (ind == 2) then
!                  G = G - integrate_Green_function(x2, x3, x1, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                elseif (ind == 3) then
!                  G = G - integrate_Green_function(x3, x1, x2, intg_pt, lambda, mu)*(wedg*faces(ifc1)%area) ! already includes area of ifc2
!                endif
!
!              endif
!
!            endif
!
!          endif
!
!        enddo
!
!      endif
!
!      sgn1 = elements(iel)%face_normal_sgn(i)
!      n1 = faces(ifc1)%normal*sgn1
!      sgn2 = elements(iel)%face_normal_sgn(j)
!      n2 = faces(ifc2)%normal*sgn2
!
!      do i1 = 1,3
!        do i2 = 1,3
!          do i3 = 1,3
!            do i4 = 1,3
!              Gamma_self(i1,i2,i3,i4) = Gamma_self(i1,i2,i3,i4) + &
!               G(i1,i3)*n1(i2)*n2(i4)/elements(iel)%V
!            enddo
!          enddo
!        enddo
!      enddo
!
!    enddo
!  enddo

  ! if (iel == 102) then
  !   write(*,*)'elements(iel)%Gamma_self', norm2(elements(iel)%Gamma_self)
  !   write(*,*)'Gamma_self', norm2(Gamma_self)
  !   stop
  ! endif

  enddo

end

subroutine calc_I_GdC_inv(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: id4(3,3,3,3), I_GdC(3,3,3,3), I_GdC_inv99(9,9)
  double precision :: I_GdC_inv3333(3,3,3,3), Cloc(3,3,3,3), Gamma_self_sym(3,3,3,3)
  double precision :: Gamma_self_sym66(6,6)
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

    Cloc = mandel_t2_to_t4(elements(iel)%c)
    I_GdC = id4 - TijklTklmn(elements(iel)%Gamma_self, Cloc - C0)
    I_GdC_inv99 = tensor4_to_tensor2(I_GdC)
    call lu_inverse(I_GdC_inv99,9)
    I_GdC_inv3333 = tensor2_to_tensor4(I_GdC_inv99)
    elements(iel)%I_GdC_inv = I_GdC_inv3333

    ! Gamma_self_sym = sym3333(elements(iel)%Gamma_self)
    ! Gamma_self_sym66 = mandel_t4_to_t2(Gamma_self_sym)

    ! elements(iel)%I_GdC_inv = id6 + matmul(Gamma_self_sym66, elements(iel)%C - C066)
    ! call lu_inverse(elements(iel)%I_GdC_inv,6) 

  enddo

end

subroutine edge_moments(edges, boxes)
  use global
  use types
  use tensor_functions

  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: iedg, i, ib, nb, j, k

  do iedg = 1,nedg

    nb = size(edges(iedg)%box_id)
    allocate(edges(iedg)%dx(3,nb))
    allocate(edges(iedg)%dxdx(3,3,nb))
    do i = 1,nb ! for each box edge belongs to
      ib = edges(iedg)%box_id(i)
      edges(iedg)%dx(:,i) = edges(iedg)%xc - boxes(ib)%x
      if (periodic) edges(iedg)%dx(:,i) = per_dx(edges(iedg)%dx(:,i))
      do j = 1,3
        do k = 1,3
          edges(iedg)%dxdx(j,k,i) = edges(iedg)%dx(j,i)*edges(iedg)%dx(k,i)
        enddo
      enddo

    enddo


  enddo

end

subroutine calc_G(edges, boxes)
  use global
  use types
  use tensor_functions

  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: iedg1, i, ib, iedg2, nb_edg, ind, j
  double precision :: dx(3), G(3,3)

  do iedg1 = 1, nedg
    nb_edg = size(edges(iedg1)%box_id)
    allocate(edges(iedg1)%G(3,3,nb_edg,max_near_edges_in_box))
  enddo

  do ib = 1,nbox
    do j = 1,boxes(ib)%nedg

      iedg1 = boxes(ib)%edge_ids(j)
      ind = findloc(edges(iedg1)%box_id, ib, 1)

      do i = 1, boxes(ib)%nnear_edge
        iedg2 = boxes(ib)%near_edge_ids(i)

        if (iedg1 .ne. iedg2) then
          dx = edges(iedg1)%xc - edges(iedg2)%xc
          G = Green_function(dx, lambda, mu)
          edges(iedg1)%G(:,:,ind,i) = G
        else
          edges(iedg1)%G(:,:,ind,i) = 0.0
        endif
    
      enddo
    enddo
  enddo

end

subroutine calc_I_C0S_inv(elements)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: id4(3,3,3,3), I_C0S_inv(3,3,3,3), I_C0S_inv99(9,9), S(3,3,3,3)
  integer :: i, iel, j, k, l

  ! do i = 1,3
  !   do j = 1,3
  !     do k = 1,3
  !       do l = 1,3
  !         id4(i,j,k,l) = id3(i,k)*id3(j,l)
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  do iel = 1,nel

    ! S = mandel_t2_to_t4(elements(iel)%S)
    ! I_C0S_inv = id4 + TijklTklmn(C0, S)
    ! I_C0S_inv99 = tensor4_to_tensor2(I_C0S_inv)
    ! call lu_inverse(I_C0S_inv99,9)
    ! I_C0S_inv = tensor2_to_tensor4(I_C0S_inv99)
    ! elements(iel)%I_C0S_inv = mandel_t4_to_t2(I_C0S_inv)

    elements(iel)%I_C0S_inv = id6 + matmul(C066, elements(iel)%S)
    call lu_inverse(elements(iel)%I_C0S_inv,6) 

  enddo

end

function point2plane(x_in_plane, normal, x) result(xplane)
  use types
  double precision, intent(in) :: x_in_plane(3), normal(3), x(3)
  double precision :: xplane(3), dx(3), dist

  dx = - x + x_in_plane ! vector from point to plane point
  dist = dot_product(dx, normal) ! normal distance to plane
  xplane = x + normal*dist
  
  ! ! test
  ! dx = - xplane + x_in_plane ! vector from point to face center
  ! dist = dot_product(dx, normal) ! normal distance to plane
  ! if (abs(dist) > 1.0e-9) stop 'point2face_plane'

end

function point_in_triangle(x1, x2, x3, area, x) result(in_triangle)
  use types
  double precision, intent(in) :: x1(3), x2(3), x3(3), area, x(3)
  double precision :: a1, a2, a3, x_vert(3,3)
  logical :: in_triangle

  x_vert(:,1) = x1
  x_vert(:,2) = x2
  x_vert(:,3) = x
  a1 = triangle_area(x_vert)

  x_vert(:,1) = x1
  x_vert(:,2) = x
  x_vert(:,3) = x3
  a2 = triangle_area(x_vert)

  x_vert(:,1) = x
  x_vert(:,2) = x2
  x_vert(:,3) = x3
  a3 = triangle_area(x_vert)

  if (abs((a1 + a2 + a3) - area) < 1.0e-9) then
    in_triangle = .true.
  else
    in_triangle = .false.
  endif

end

subroutine project_to_triangle(x1, x2, x3, x, side, xproj, nproj) ! ind tells us which side (1 for 1-2, 2 for 2-3, 3 for 3-1)
  use types
  use global, only : min_dist_node
  double precision, intent(in) :: x1(3), x2(3), x3(3), x(3)
  double precision, intent(out) :: xproj(3)
  integer, intent(out) :: side, nproj
  double precision :: v1(3), v2(3), v3(3), e1, e2, e3, proj, v(3), xproj_arr(3,3), dist_arr(3)
  integer :: ic, side_arr(3), ind
  logical :: in_triangle

  v1 = x2 - x1 ! from x1 to x2
  v2 = x3 - x2
  v3 = x1 - x3
  e1 = norm2(v1)
  e2 = norm2(v2)
  e3 = norm2(v3)
  v1 = v1/e1
  v2 = v2/e2
  v3 = v3/e3

  ic = 0
  v = x - x1 ! from x1 to x
  proj = dot_product(v, v1)
  if (proj > min_dist_node .and. proj < e1 - min_dist_node) then
    ic = ic + 1
    xproj_arr(:,ic) = x1 + v1*proj
    dist_arr(ic) = norm2(xproj_arr(:,ic) - x)
    side_arr(ic) = 1
  endif

  v = x - x2 ! from x2 to x
  proj = dot_product(v, v2)
  if (proj > min_dist_node .and. proj < e2 - min_dist_node) then
    ic = ic + 1
    xproj_arr(:,ic) = x2 + v2*proj
    dist_arr(ic) = norm2(xproj_arr(:,ic) - x)
    side_arr(ic) = 2
  endif

  v = x - x3 ! from x3 to x
  proj = dot_product(v, v3)
  if (proj > min_dist_node .and. proj < e3 - min_dist_node) then
    ic = ic + 1
    xproj_arr(:,ic) = x3 + v3*proj
    dist_arr(ic) = norm2(xproj_arr(:,ic) - x)
    side_arr(ic) = 3
  endif

  if (ic == 0) then
    ind = minloc([norm2(x - x1), norm2(x - x2), norm2(x - x3)],1)
    nproj = 1
    if (ind == 1) then
      xproj = x1
      side = 1
    elseif (ind == 2) then
      xproj = x2
      side = 2
    elseif (ind == 3) then
      xproj = x3
      side = 3
    endif
  else
    ind = minloc(dist_arr(1:ic),1)
    xproj = xproj_arr(:,ind)
    side = side_arr(ind)
    nproj = ic
  endif

end

subroutine test_geometry(faces, elements)
  use global
  use types
  use tensor_functions

  type (face_type), allocatable, intent (inout) :: faces(:)
  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: x(3), xplane(3), x1(3), x2(3), x3(3), xproj(3)
  integer :: ifc, iel, ind, nproj
  logical :: in_triangle

  ifc = 44
  iel = faces(ifc)%el(1)
  call get_face_nodes(faces(ifc), elements(iel), x1, x2, x3)
  write(*,*) 'x1', x1
  write(*,*) 'x2', x2
  write(*,*) 'x3', x3

  x = [5.0,-13.0,155.0]
  xplane = point2plane(faces(ifc)%x, faces(ifc)%normal, x)
  
  write(*,*) xplane

  in_triangle = point_in_triangle(x1, x2, x3, faces(ifc)%area, xplane)
  write(*,*) in_triangle

  in_triangle = point_in_triangle(x1, x2, x3, faces(ifc)%area, faces(ifc)%x)
  write(*,*) in_triangle

  call project_to_triangle(x1, x2, x3, xplane, ind, xproj, nproj)
  write(*,*) xproj

  call project_to_triangle(x1, x2, x3, faces(ifc)%x, ind, xproj, nproj)
  write(*,*) xproj

end


subroutine integrate_Green_function_edges_close(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (in) :: boxes(:)
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1face(3), x2face(3), x3face(3), xnode_face(3,3), G(3,3), xproj(3)
  double precision :: x1(3), x2(3), x3(3), dx(3), xedge(3), xshift(3), xface(3), dist_node(3)
  double precision :: x1quad(3), x2quad(3), x3quad(3), x4quad(3)
  double precision :: alpha1 = 0.0597158717, beta1 = 0.4701420641
  double precision :: alpha2 = 0.7974269853, beta2 = 0.1012865073
  integer :: iedg1, i ,ifc2, iel2, j, face_edge, ib1, iedg2, ib2, ic, ind, nproj, side, intgpt, ib
  integer :: ib_in_near, ibnear, ifc_in_box
  logical :: in_triangle

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,ic,ifc2,dx,xface,iel2,in_triangle,x1,x2,x3,xproj,nproj,side,ib,ib_in_near,ibnear,ifc_in_box) &
  !$OMP SHARED(nedg,nfc,edges,faces,elements,dx_th,boxes)
  !$OMP DO 
  do iedg1 = 1, nedg

    ib = edges(iedg1)%box_id_unique

    ic = 0
    ! write(*,*)'iedg1 :', iedg1
    do ib_in_near = 1, boxes(ib)%nnear
      ibnear = boxes(ib)%near_box_id(ib_in_near)

      do ifc_in_box = 1, boxes(ibnear)%nfc
        ifc2 = boxes(ibnear)%face_ids(ifc_in_box)

!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2)) then ! skip faces that are shared by two boxes if they are not in box containing the edge
        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2 .and. &
          (faces(ifc2)%box_id(1) == ib .or. faces(ifc2)%box_id(2) == ib))) then

          dx = edges(iedg1)%xc - faces(ifc2)%x
          xface = point2plane(faces(ifc2)%x, faces(ifc2)%normal, edges(iedg1)%xc) ! project edge center to face ifc2
          iel2 = faces(ifc2)%el(1)
          call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
          call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj) ! project the in-plane point to triangle
          in_triangle = point_in_triangle(x1, x2, x3, faces(ifc2)%area, xface) 
          if (((in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th) .or. &
              (.not. in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th .and. norm2(xface - xproj) < dx_th*0.75)) &
             .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th) then
            ic = ic + 1
          endif

        endif
  
      enddo
    enddo
    if (ic > 0) then
      ! write(*,*)'iedg1, ic :', iedg1, ic
      allocate(edges(iedg1)%face_close(ic))
      allocate(edges(iedg1)%face_close_G(3,3,ic))
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,ic,ifc2,dx,xface,iel2,in_triangle,x1,x2,x3,xproj,nproj,side) &
  !$OMP PRIVATE(G,x1quad,x2quad,x3quad,x4quad,ib,ib_in_near,ibnear,ifc_in_box) &
  !$OMP SHARED(nedg,nfc,edges,faces,elements,dx_th,intg_pt_quad,lambda,mu,boxes)
  !$OMP DO 
  do iedg1 = 1, nedg

    ib = edges(iedg1)%box_id_unique

    ic = 0
    ! write(*,*)'iedg1 :', iedg1
    do ib_in_near = 1, boxes(ib)%nnear
      ibnear = boxes(ib)%near_box_id(ib_in_near)

      do ifc_in_box = 1, boxes(ibnear)%nfc
        ifc2 = boxes(ibnear)%face_ids(ifc_in_box)

!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2)) then ! skip faces that are shared by two boxes if they are not in box containing the edge
        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2 .and. &
          (faces(ifc2)%box_id(1) == ib .or. faces(ifc2)%box_id(2) == ib))) then
    
          dx = edges(iedg1)%xc - faces(ifc2)%x
          xface = point2plane(faces(ifc2)%x, faces(ifc2)%normal, edges(iedg1)%xc) ! project edge center to face ifc2
          iel2 = faces(ifc2)%el(1)
          call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
          call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj) ! project the in-plane point to triangle
          in_triangle = point_in_triangle(x1, x2, x3, faces(ifc2)%area, xface) 
          if (((in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th) .or. &
              (.not. in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th .and. norm2(xface - xproj) < dx_th*0.75)) &
             .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th) then
            ! write(*,*)'ifc2 :', ifc2
            if (findloc(edges(iedg1)%face_close(1:ic), ifc2, 1) .ne. 0) &
             stop 'findloc(edges(iedg1)%face_close(1:ic), ifc2, 1) .ne. 0'
            ic = ic + 1
            edges(iedg1)%face_close(ic) = ifc2
    
            ! nearly singular
            iel2 = faces(ifc2)%el(1)
            call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
            ! write(33,'(A)') 'triangle_x = np.empty([3,4])'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,0] = [',x1(1),',',x1(2),',',x1(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,1] = [',x2(1),',',x2(2),',',x2(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,2] = [',x3(1),',',x3(2),',',x3(3),']'
            ! write(33,'(A)') 'triangle_x[:,3] = triangle_x[:,0]'
            ! write(33,'(A)') 'ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])'
            G = 0.0
            x1quad = x1
            x2quad = (x1 + x2)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x1 + x3)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'x0 = np.array([',edges(iedg1)%xc(1),',',edges(iedg1)%xc(2),',',edges(iedg1)%xc(3),'])'
            ! write(33,'(A)') 'ax.scatter(x0[0], x0[1], x0[2])'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad1_x[:,4] = quad1_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad1_x[0,:], quad1_x[1,:], quad1_x[2,:])'
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            x1quad = x2
            x2quad = (x2 + x3)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x2 + x1)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad2_x[:,4] = quad2_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad2_x[0,:], quad2_x[1,:], quad2_x[2,:])'
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            x1quad = x3
            x2quad = (x3 + x1)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x3 + x2)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad3_x[:,4] = quad3_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad3_x[0,:], quad3_x[1,:], quad3_x[2,:])'
            ! write(33,'(A)') 'plt.show()'
            ! close(33)
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            G = G/faces(ifc2)%area
            edges(iedg1)%face_close_G(:,:,ic) = G

          endif
  
        endif
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end


subroutine integrate_Green_function_edges_close_test(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (in) :: boxes(:)
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1face(3), x2face(3), x3face(3), xnode_face(3,3), G(3,3), xproj(3)
  double precision :: x1(3), x2(3), x3(3), dx(3), xedge(3), xshift(3), xface(3), dist_node(3)
  double precision :: x1quad(3), x2quad(3), x3quad(3), x4quad(3)
  double precision :: alpha1 = 0.0597158717, beta1 = 0.4701420641
  double precision :: alpha2 = 0.7974269853, beta2 = 0.1012865073
  integer, allocatable :: edge_face_close(:)
  integer :: iedg1, i ,ifc2, iel2, j, face_edge, ib1, iedg2, ib2, ic, ind, nproj, side, intgpt, ib
  integer :: ib_in_near, ibnear, ifc_in_box
  logical :: in_triangle

  allocate(edge_face_close(10000))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,ic,ifc2,dx,xface,iel2,in_triangle,x1,x2,x3,xproj,nproj,side,ib,ib_in_near,ibnear,ifc_in_box) &
  !$OMP PRIVATE(edge_face_close) &
  !$OMP SHARED(nedg,nfc,edges,faces,elements,dx_th,boxes)
  !$OMP DO 
  do iedg1 = 1, nedg

    ib = edges(iedg1)%box_id_unique

    ic = 0
    ! write(*,*)'iedg1 :', iedg1
    do ib_in_near = 1, boxes(ib)%nnear
      ibnear = boxes(ib)%near_box_id(ib_in_near)

      do ifc_in_box = 1, boxes(ibnear)%nfc
        ifc2 = boxes(ibnear)%face_ids(ifc_in_box)

!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2)) then ! skip faces that are shared by two boxes if they are not in box containing the edge
!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2 .and. &
!          (faces(ifc2)%box_id(1) == ib .or. faces(ifc2)%box_id(2) == ib))) then
        if (findloc(edge_face_close(1:ic), ifc2, 1) == 0) then

          dx = edges(iedg1)%xc - faces(ifc2)%x
          xface = point2plane(faces(ifc2)%x, faces(ifc2)%normal, edges(iedg1)%xc) ! project edge center to face ifc2
          iel2 = faces(ifc2)%el(1)
          call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
          call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj) ! project the in-plane point to triangle
          in_triangle = point_in_triangle(x1, x2, x3, faces(ifc2)%area, xface) 
          if (((in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th) .or. &
              (.not. in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th .and. norm2(xface - xproj) < dx_th*0.75)) &
             .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
          ! if (norm2(dx) <= dx_th) then
            ic = ic + 1
            edge_face_close(ic) = ifc2
          endif

        endif
  
      enddo
    enddo
    if (ic > 0) then
      ! write(*,*)'iedg1, ic :', iedg1, ic
      allocate(edges(iedg1)%face_close(ic))
      edges(iedg1)%face_close = edge_face_close(1:ic)
      allocate(edges(iedg1)%face_close_G(3,3,ic))
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iedg1,ic,ifc2,dx,xface,iel2,in_triangle,x1,x2,x3,xproj,nproj,side) &
  !$OMP PRIVATE(G,x1quad,x2quad,x3quad,x4quad,ib,ib_in_near,ibnear,ifc_in_box) &
  !$OMP SHARED(nedg,nfc,edges,faces,elements,dx_th,intg_pt_quad,lambda,mu,boxes)
  !$OMP DO 
  do iedg1 = 1, nedg

    if (allocated(edges(iedg1)%face_close)) then

    ib = edges(iedg1)%box_id_unique

    ic = 0
!    ! write(*,*)'iedg1 :', iedg1
!    do ib_in_near = 1, boxes(ib)%nnear
!      ibnear = boxes(ib)%near_box_id(ib_in_near)
!
!      do ifc_in_box = 1, boxes(ibnear)%nfc
!        ifc2 = boxes(ibnear)%face_ids(ifc_in_box)
    do i = 1,size(edges(iedg1)%face_close)
      ifc2 = edges(iedg1)%face_close(i)

!!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2)) then ! skip faces that are shared by two boxes if they are not in box containing the edge
!        if (.not. (ibnear .ne. ib .and. faces(ifc2)%nbox == 2 .and. &
!          (faces(ifc2)%box_id(1) == ib .or. faces(ifc2)%box_id(2) == ib))) then
!    
!          dx = edges(iedg1)%xc - faces(ifc2)%x
!          xface = point2plane(faces(ifc2)%x, faces(ifc2)%normal, edges(iedg1)%xc) ! project edge center to face ifc2
!          iel2 = faces(ifc2)%el(1)
!          call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
!          call project_to_triangle(x1, x2, x3, xface, side, xproj, nproj) ! project the in-plane point to triangle
!          in_triangle = point_in_triangle(x1, x2, x3, faces(ifc2)%area, xface) 
!          if (((in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th) .or. &
!              (.not. in_triangle .and. norm2(xface - edges(iedg1)%xc) < dx_th .and. norm2(xface - xproj) < dx_th*0.75)) &
!             .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
!          ! if (norm2(dx) <= dx_th .and. findloc(edges(iedg1)%face, ifc2, 1) == 0) then
!          ! if (norm2(dx) <= dx_th) then
            ! write(*,*)'ifc2 :', ifc2
            ic = ic + 1
            ! edges(iedg1)%face_close(ic) = ifc2
    
            ! nearly singular
            iel2 = faces(ifc2)%el(1)
            call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
            ! write(33,'(A)') 'triangle_x = np.empty([3,4])'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,0] = [',x1(1),',',x1(2),',',x1(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,1] = [',x2(1),',',x2(2),',',x2(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,2] = [',x3(1),',',x3(2),',',x3(3),']'
            ! write(33,'(A)') 'triangle_x[:,3] = triangle_x[:,0]'
            ! write(33,'(A)') 'ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])'
            G = 0.0
            x1quad = x1
            x2quad = (x1 + x2)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x1 + x3)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'x0 = np.array([',edges(iedg1)%xc(1),',',edges(iedg1)%xc(2),',',edges(iedg1)%xc(3),'])'
            ! write(33,'(A)') 'ax.scatter(x0[0], x0[1], x0[2])'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad1_x[:,4] = quad1_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad1_x[0,:], quad1_x[1,:], quad1_x[2,:])'
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            x1quad = x2
            x2quad = (x2 + x3)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x2 + x1)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad2_x[:,4] = quad2_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad2_x[0,:], quad2_x[1,:], quad2_x[2,:])'
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            x1quad = x3
            x2quad = (x3 + x1)*0.5
            x3quad = faces(ifc2)%x
            x4quad = (x3 + x2)*0.5
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
            ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
            ! write(33,'(A)') 'quad3_x[:,4] = quad3_x[:,0]'
            ! write(33,'(A)') 'ax.plot(quad3_x[0,:], quad3_x[1,:], quad3_x[2,:])'
            ! write(33,'(A)') 'plt.show()'
            ! close(33)
            G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
             x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
            G = G/faces(ifc2)%area
            edges(iedg1)%face_close_G(:,:,ic) = G

!          endif
!  
!        endif
!      enddo
    enddo

    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(edge_face_close)

end

subroutine integrate_Green_function_edges_dist_el(elements, faces, edges, boxes, intg_pt, intg_pt_quad)
  use global
  use types
  use tensor_functions

  type (element_type), allocatable, intent (in) :: elements(:)
  type (face_type), allocatable, intent (in) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)
  type (box_type), allocatable, intent (in) :: boxes(:)
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1face(3), x2face(3), x3face(3), xnode_face(3,3), G(3,3), xproj(3)
  double precision :: x1(3), x2(3), x3(3), dx(3), xedge(3), xshift(3), xface(3), dist_node(3)
  double precision :: x1quad(3), x2quad(3), x3quad(3), x4quad(3)
  double precision :: alpha1 = 0.0597158717, beta1 = 0.4701420641
  double precision :: alpha2 = 0.7974269853, beta2 = 0.1012865073
  integer, allocatable :: edge_face_close(:,:), nedge_face_close(:)
  integer :: iedg1, i ,ifc2, iel2, j, face_edge, ib1, iedg2, ib2, ic, ind, nproj, side, intgpt, ib
  integer :: ib_in_near, ibnear, ifc_in_box, iel, ifc1, k, idum
  logical :: in_triangle

  allocate(edge_face_close(200,nedg))
  allocate(nedge_face_close(nedg))

  nedge_face_close = 0
  do iel = 1, nel

    if (elements(iel)%aspect > aspect_mx) then

      ! write(*,*) 'iel', iel
      do i = 1,4
        ifc1 = elements(iel)%face(i)
        ! write(*,*) '  ifc1', ifc1
        ! write(*,*) '  edges in ifc1', faces(ifc1)%edge
        do j = 1,4
          ifc2 = elements(iel)%face(j)
          do k = 1,3
            iedg2 = faces(ifc2)%edge(k)
            ! write(*,*) '    iedg2', iedg2
            idum = nedge_face_close(iedg2)
            if (findloc(faces(ifc1)%edge, iedg2, 1) == 0 .and. & ! iedg2 does not belong to ifc1
                findloc(edge_face_close(1:idum, iedg2), ifc1, 1) == 0) then ! ifc1 was not previously added to iedg2 array
              nedge_face_close(iedg2) = nedge_face_close(iedg2) + 1
              edge_face_close(nedge_face_close(iedg2), iedg2) = ifc1
              ! write(*,*) '    in'
            endif
          enddo
        enddo
      enddo
    endif

  enddo

  do iedg1 = 1, nedg

    idum = nedge_face_close(iedg1)
    if (idum > 0) then
      allocate(edges(iedg1)%face_close(idum))
      allocate(edges(iedg1)%face_close_G(3,3,idum))
      edges(iedg1)%face_close(1:idum) = edge_face_close(1:idum,iedg1)

      do i = 1,idum
        ifc2 = edges(iedg1)%face_close(i)

        ! nearly singular
        iel2 = faces(ifc2)%el(1)
        call get_face_nodes(faces(ifc2), elements(iel2), x1, x2, x3) ! get nodes
        ! write(33,'(A)') 'triangle_x = np.empty([3,4])'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,0] = [',x1(1),',',x1(2),',',x1(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,1] = [',x2(1),',',x2(2),',',x2(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,2] = [',x3(1),',',x3(2),',',x3(3),']'
        ! write(33,'(A)') 'triangle_x[:,3] = triangle_x[:,0]'
        ! write(33,'(A)') 'ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])'
        G = 0.0
        x1quad = x1
        x2quad = (x1 + x2)*0.5
        x3quad = faces(ifc2)%x
        x4quad = (x1 + x3)*0.5
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'x0 = np.array([',edges(iedg1)%xc(1),',',edges(iedg1)%xc(2),',',edges(iedg1)%xc(3),'])'
        ! write(33,'(A)') 'ax.scatter(x0[0], x0[1], x0[2])'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad1_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
        ! write(33,'(A)') 'quad1_x[:,4] = quad1_x[:,0]'
        ! write(33,'(A)') 'ax.plot(quad1_x[0,:], quad1_x[1,:], quad1_x[2,:])'
        G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
         x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
        x1quad = x2
        x2quad = (x2 + x3)*0.5
        x3quad = faces(ifc2)%x
        x4quad = (x2 + x1)*0.5
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad2_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
        ! write(33,'(A)') 'quad2_x[:,4] = quad2_x[:,0]'
        ! write(33,'(A)') 'ax.plot(quad2_x[0,:], quad2_x[1,:], quad2_x[2,:])'
        G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
         x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
        x1quad = x3
        x2quad = (x3 + x1)*0.5
        x3quad = faces(ifc2)%x
        x4quad = (x3 + x2)*0.5
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,0] = [',x1quad(1),',',x1quad(2),',',x1quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,1] = [',x2quad(1),',',x2quad(2),',',x2quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,2] = [',x3quad(1),',',x3quad(2),',',x3quad(3),']'
        ! write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad3_x[:,3] = [',x4quad(1),',',x4quad(2),',',x4quad(3),']'
        ! write(33,'(A)') 'quad3_x[:,4] = quad3_x[:,0]'
        ! write(33,'(A)') 'ax.plot(quad3_x[0,:], quad3_x[1,:], quad3_x[2,:])'
        ! write(33,'(A)') 'plt.show()'
        ! close(33)
        G = G + integrate_Green_function_nearly_sing_quad(edges(iedg1)%xc, &
         x1quad, x2quad, x3quad, x4quad, intg_pt_quad, lambda, mu)
        G = G/faces(ifc2)%area
        edges(iedg1)%face_close_G(:,:,i) = G
      enddo

    endif
  enddo

  deallocate(edge_face_close)
  deallocate(nedge_face_close)

end

function calc_coef(d, f_0) result(coef)
  use tensor_functions
  double precision, intent(in) :: d, f_0
  double precision :: coef(4), mat(2,2), vec(2), f_d, dfdr_0, dfdr_d

  f_d = 1.0/d
  dfdr_0 = 0.0
  dfdr_d = -d**(-2)

  coef(4) = f_0
  coef(3) = 0.0

  mat(1,:) = [3.0*d**2, 2.0*d]
  mat(2,:) = [d**3, d**2]
  vec = [dfdr_d,f_d - coef(4)]

  call lu_inverse(mat,2)
  coef(1:2) = matmul(mat,vec)

end

subroutine test_Green_function_integration_quad(intg_pt, intg_pt_quad)
  use global
  use types
  type (intg_pt_type), intent(in) :: intg_pt(nintg_pt)
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1(3), x2(3), x3(3), G_Gauss(3,3), Gintg(3,3), dx(3)
  double precision :: xintpt(3,nintg_pt), fintpt(nintg_pt), jacob, x1pt, x2pt, x3pt, x123(3,3)
  double precision :: dx2, dx3, x(3), xhat(2), xbar(2), jacob_LW, jacob_ref2triang, wtot, w, oneoverx
  integer :: i, ip1, ip2, ip3, np1, np2, np3, iintg
  character (len=10) :: str

  open (unit = 33, file = 'plot_integration_pt.py', status = 'unknown')
  
  write(33,'(A)') 'import matplotlib.pyplot as plt'
  write(33,'(A)') 'import numpy as np'
  write(33,'(A)') 'import os'
  write(33,'(A)') 'import re'
  write(33,'(A)')
  write(str,'(I2)') nintg_pt_quad
  write(33,'(A)') 'xintgpt = np.empty([3,'//trim(str)//'])'

  x1 = 0.0
  x2 = [0.0,1.0,0.0]*0.25
  x3 = [0.0,1.0,1.0]*0.25

  write(*,*) reference_to_triangle([0.0, 0.0], x1, x2, x3)
  write(*,*) reference_to_triangle([2.0, 0.0], x1, x2, x3)
  write(*,*) reference_to_triangle([2.0, 2.0], x1, x2, x3)

  oneoverx = 0.0
  wtot = 0.0
  do iintg = 1, nintg_pt_quad

    xbar = intg_pt_quad(iintg)%Xhat
    xhat = Lachat_Watson_transform(xbar)
    x = reference_to_triangle(xhat, x1, x2, x3)

    jacob_LW = Lachat_Watson_transform_jacob(xbar)
    jacob_ref2triang = reference_to_triangle_jacob(xhat, x1, x2, x3)
    w = intg_pt_quad(iintg)%w*abs(jacob_LW)*abs(jacob_ref2triang)
    wtot = wtot + w

    dx = x1 - x 
    oneoverx = oneoverx + 1.0/norm2(dx)*w

    write(33,'(A,I2,A,E16.8,A,E16.8,A,E16.8,A)')'xintgpt[:,',iintg-1,'] = np.array([',x(1),',',x(2),',',x(3),'])'

  enddo

  x123(:,1) = x1
  x123(:,2) = X2
  x123(:,3) = x3
  write(*,*) wtot, triangle_area(x123)
  write(*,*) oneoverx

  write(33,'(A)') 'triangle_x = np.empty([3,4])'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,0] = [',x1(1),',',x1(2),',',x1(3),']'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,1] = [',x2(1),',',x2(2),',',x2(3),']'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'triangle_x[:,2] = [',x3(1),',',x3(2),',',x3(3),']'
  write(33,'(A)') 'triangle_x[:,3] = triangle_x[:,0]'

  write(33,'(A)') 'fig = plt.figure()'
  write(33,'(A)') 'ax = fig.add_subplot(111, projection=''3d'')'
  write(33,'(A)') 'ax.plot(triangle_x[0,:], triangle_x[1,:], triangle_x[2,:])'
  write(33,'(A)') 'ax.scatter(xintgpt[0,:], xintgpt[1,:], xintgpt[2,:])'
  write(33,'(A)') 'ax.axis(''equal'')'
  write(33,'(A)') 'plt.show()'
  close(33)

  write(*,*) 'G (Pina et al. 1981)'
  G_Gauss = integrate_Green_function(x1, x2, x3, intg_pt, lambda, mu)
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)
  G_Gauss = integrate_Green_function(x1, x3, x2, intg_pt, lambda, mu) ! flip nodes, should give same
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)
  
  write(*,*) 'G (Lachat-Watson transform)'
  G_Gauss = integrate_Green_function_quad(x1, x2, x3, intg_pt_quad, lambda, mu)
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)
  G_Gauss = integrate_Green_function_quad(x1, x3, x2, intg_pt_quad, lambda, mu) ! flip nodes, should give same
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)

  np2 = 1000
  np3 = 1000
  dx2 = 0.25/float(np2)
  dx3 = 0.25/float(np3)
  x2pt = 0.0
  x3pt = 0.0
  x1pt = 0.0
  Gintg = 0.0
  do ip2 = 1,np2
    x2pt = float(ip2)*dx2 - dx2*0.5
    do ip3 = 1,ip2
      x3pt = float(ip3)*dx3 - dx3*0.5
      x(1) = x1pt
      x(2) = x2pt
      x(3) = x3pt
      ! write(*,*)x2pt,x3pt
      Gintg = Gintg + Green_function(-x, lambda, mu)*dx2*dx3
    enddo
  enddo
  write(*,*) 'G (regular grid)'
  write(*,'(3E16.8)') (Gintg(i,:),i=1,3)

  ! stop

end

function Lachat_Watson_transform(Xbar) result(Xhat)
  double precision, intent(in) :: Xbar(2)
  double precision :: Xhat(2)

  Xhat(1) = Xbar(1) + 1.0
  Xhat(2) = 0.5*(Xbar(1) + 1.0)*(Xbar(2) + 1.0)

end

function Lachat_Watson_transform_jacob(Xbar) result(jacob)
  double precision, intent(in) :: Xbar(2)
  double precision :: jacob

  jacob = 0.5*(Xbar(1) + 1.0)

end

function reference_to_triangle(Xhat, X1, X2, X3) result(X)
  double precision, intent(in) :: Xhat(2), X1(3), X2(3), X3(3)
  double precision :: X(3), Xa(3), Xb(3), Xc(3)

  Xa = X1
  Xb = 0.5*(X2 - X1)
  Xc = 0.5*(X3 - X2)
  X = Xa + Xb*Xhat(1) + Xc*Xhat(2)

end

function reference_to_triangle_jacob(Xhat, X1, X2, X3) result(jacob)
  use tensor_functions
  double precision, intent(in) :: Xhat(2), X1(3), X2(3), X3(3)
  double precision :: n(3), e1(3), e2(3), e3(3), Qg2l(3,3), X1loc(3), X2loc(3)
  double precision :: X3loc(3), Xaloc(3), Xbloc(3), Xcloc(3), J22(2,2), jacob

  n = cross_prod(X2 - X1, X3 - X1)
  n = n/norm2(n)

  e3 = n
  e1 = x2 - x1
  e1 = e1/norm2(e1)
  e2 = cross_prod(e3, e1)

  Qg2l(1,:) = e1
  Qg2l(2,:) = e2
  Qg2l(3,:) = e3

  X1loc = matmul(Qg2l, X1)
  X2loc = matmul(Qg2l, X2)
  X3loc = matmul(Qg2l, X3)  

  Xaloc = X1loc
  Xbloc = 0.5*(X2loc - X1loc)
  Xcloc = 0.5*(X3loc - X2loc)

  if (.not. (abs(X1loc(3) - X2loc(3)) < 1.0e-7 .and. abs(X1loc(3) - X3loc(3)) < 1.0e-7) ) then
    write(*,*) abs(X1loc(3) - X2loc(3))
    write(*,*) abs(X1loc(3) - X3loc(3))
    stop 'third coordinat not the same in jacobian'
  endif

  J22(:,1) = Xbloc(1:2)
  J22(:,2) = Xcloc(1:2)
  jacob = abs(J22(1,1)*J22(2,2) - J22(1,2)*J22(2,1))

end

function integrate_Green_function_quad(x1, x2, x3, intg_pt_quad, lambda, mu) result(Gintg) 
  use types
  use global, only : nintg_pt_quad
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)
  double precision, intent(in) :: x1(3), x2(3), x3(3), lambda, mu
  double precision :: Gintg(3,3), GxR(3,3), dx(3), x123(3,3)
  double precision :: jacob, wtot, xbar(2), xhat(2), x(3), jacob_LW, jacob_ref2triang, w
  integer :: ip

  Gintg = 0.0
  wtot = 0.0
  do ip = 1, nintg_pt_quad

    xbar = intg_pt_quad(ip)%Xhat
    xhat = Lachat_Watson_transform(xbar)
    x = reference_to_triangle(xhat, x1, x2, x3)

    jacob_LW = Lachat_Watson_transform_jacob(xbar)
    jacob_ref2triang = reference_to_triangle_jacob(xhat, x1, x2, x3)
    w = intg_pt_quad(ip)%w*abs(jacob_LW)*abs(jacob_ref2triang)
    wtot = wtot + w

    dx = x1 - x 
    GxR = Green_function_x_R(dx, lambda, mu)
    Gintg = Gintg + GxR/norm2(dx)*w

  enddo

  x123(:,1) = x1
  x123(:,2) = X2
  x123(:,3) = x3
  if (abs(wtot - triangle_area(x123)) > 1.0e-7) then
    stop 'abs(wtot - triangle_area(x123)) > 1.0e-7'
  endif

end

function sinh_transform(Xbar, X0hat, b, mu, eta) result(Xhat)
  double precision, intent(in) :: Xbar(2), X0hat(2), b, mu(2), eta(2)
  double precision :: Xhat(2)

  Xhat(1) = X0hat(1) + b*sinh(mu(1)*Xbar(1) - eta(1))
  Xhat(2) = X0hat(2) + b*sinh(mu(2)*Xbar(2) - eta(2))

end

function sinh_transform_jacob(Xbar, X0hat, b, mu, eta) result(jacob)
  double precision, intent(in) :: Xbar(2), X0hat(2), b, mu(2), eta(2)
  double precision :: jacob

  jacob = abs((b*cosh(mu(1)*Xbar(1) - eta(1))*mu(1))*(b*cosh(mu(2)*Xbar(2) - eta(2))*mu(2)))

end

function reference_to_quad_3D(Xhat, X1, X2, X3, X4) result(X)
  double precision, intent(in) :: Xhat(2), X1(3), X2(3), X3(3), X4(3)
  double precision :: X(3), Xa(3), Xb(3), Xc(3), Xd(3)

  call calc_Xabcd(X1, X2, X3, X4, Xa, Xb, Xc, Xd)

  X = Xa + Xb*Xhat(1) + Xc*Xhat(2) + Xd*Xhat(1)*Xhat(2)

end

function reference_to_quad_2D(Xhat, X1, X2, X3, X4) result(X)
  double precision, intent(in) :: Xhat(2), X1(2), X2(2), X3(2), X4(2)
  double precision :: X(2), Xa(2), Xb(2), Xc(2), Xd(2)

  call calc_Xabcd(X1, X2, X3, X4, Xa, Xb, Xc, Xd)

  X = Xa + Xb*Xhat(1) + Xc*Xhat(2) + Xd*Xhat(1)*Xhat(2)

end

function reference_to_quad_2D_jacob(Xhat, X1, X2, X3, X4) result(jacob)
  double precision, intent(in) :: Xhat(2), X1(2), X2(2), X3(2), X4(2)
  double precision :: X(2), Xa(2), Xb(2), Xc(2), Xd(2), J22(2,2), jacob

  call calc_Xabcd(X1, X2, X3, X4, Xa, Xb, Xc, Xd)

  J22(:,1) = Xb + Xd*Xhat(2)
  J22(:,2) = Xc + Xd*Xhat(1)

  jacob = abs(J22(1,1)*J22(2,2) - J22(1,2)*J22(2,1))

end

function quad_to_reference_2D(X, X1, X2, X3, X4) result(Xhat)
  use tensor_functions
  double precision, intent(in) :: X(2), X1(2), X2(2), X3(2), X4(2)
  double precision :: Xhat(2), Xa(2), Xb(2), Xc(2), Xd(2), dRdXhat_inv(2,2), R(2), Rn, dRdXhat(2,2)
  double precision :: tol = 1.0e-5, Roldn, Xhatold(2), dXhat(2), e1(2), e2(2), d1, d2, Xhat_init(2)
  integer :: iter, i 

  call calc_Xabcd(X1, X2, X3, X4, Xa, Xb, Xc, Xd)

  Rn = 1.0
  ! Xhat = X
  ! Xhat = [0.0,0.0]
  e1 = -0.25*(X1 + X2 + X3 + X4) + 0.5*(X2 + X3)
  d1 = norm2(e1)
  e1 = e1/d1
  Xhat_init(1) = dot_product(X,e1)/d1
  e2 = -0.25*(X1 + X2 + X3 + X4) + 0.5*(X3 + X4)
  d2 = norm2(e2)
  e2 = e2/d2
  Xhat_init(2) = dot_product(X,e2)/d2
  Xhat = Xhat_init
  ! write(*,*) 'X:', X
  ! write(*,*) 'Xhat (initial):', Xhat
  do iter = 1, 10

    R = reference_to_quad_2D(Xhat, X1, X2, X3, X4) - X
    Roldn = norm2(R)
    ! write(*,*) 'Rn', Roldn
    if (Roldn < tol) then
      ! write(*,*) iter, Roldn
      return
    endif

    dRdXhat(:,1) = Xb + Xd*Xhat(2)
    dRdXhat(:,2) = Xc + Xd*Xhat(1)


    ! dRdXhat_inv = dRdXhat
    ! call lu_inverse(dRdXhat_inv,2)
    dRdXhat_inv(1,1) = dRdXhat(2,2)
    dRdXhat_inv(2,2) = dRdXhat(1,1)
    dRdXhat_inv(1,2) = -dRdXhat(1,2)
    dRdXhat_inv(2,1) = -dRdXhat(2,1)
    if (dRdXhat(1,1)*dRdXhat(2,2) - dRdXhat(1,2)*dRdXhat(2,1) < 1.0e-15) exit
    dRdXhat_inv = dRdXhat_inv/(dRdXhat(1,1)*dRdXhat(2,2) - dRdXhat(1,2)*dRdXhat(2,1))

    dXhat = - matmul(dRdXhat_inv, R)

    ! solve eq. sys.
    ! call lu_eqsystem(dRdXhat,R,2)
    ! dXhat = - R

    Xhatold = Xhat
    do i = 1,10
      Xhat = Xhatold + dXhat
      R = reference_to_quad_2D(Xhat, X1, X2, X3, X4) - X
      Rn = norm2(R)
      ! write(*,*) ' i, dXhat', i, dXhat
      ! write(*,*) ' i, Rn', i, Rn
      if (Rn < Roldn) then
        exit
      else
        dXhat = dXhat*0.5
      endif
    enddo

    if (Rn < tol) then
      ! write(*,*) 'Xhat solution', Xhat
      ! write(*,*) iter, Rn
      return
    endif

  enddo
  ! write(*,*) iter, Rn
  ! Xhat = Xhat_init
  Xhat = -100.0
  ! write(*,*) 'quad_to_reference_2D not converged'

end

subroutine calc_Xabcd(X1, X2, X3, X4, Xa, Xb, Xc, Xd)
  double precision, intent(in) :: X1(:), X2(:), X3(:), X4(:)
  double precision, dimension(size(X1)), intent(out) :: Xa(:), Xb(:), Xc(:), Xd(:)

  Xa = 0.25*(X1 + X2 + X3 + X4)
  Xb = 0.25*(-X1 + X2 + X3 - X4)
  Xc = 0.25*(-X1 - X2 + X3 + X4)
  Xd = 0.25*(X1 - X2 + X3 - X4)

end

subroutine transformation_mat_to_local_quad(X1, X2, X3, X4, Xcenter, Qg2l)
  use tensor_functions
  double precision, intent(in) :: X1(3), X2(3), X3(3), X4(3)
  double precision, intent(out) :: Xcenter(3), Qg2l(3,3)
  double precision :: e1(3), e2(3), e3(3), n(3)

  Xcenter = 0.25*(X1 + X2 + X3 + X4)

  n = cross_prod(x2 - x1, x4 - x1)
  n = n/norm2(n)

  e3 = n
  e1 = x2 - x1
  e1 = e1/norm2(e1)
  e2 = cross_prod(e3, e1)
  Qg2l(1,:) = e1
  Qg2l(2,:) = e2
  Qg2l(3,:) = e3

end

function transform_to_local_quad(X, Xcenter, Qg2l) result(Xloc)
  ! transformation to local frame attached to quad center
  double precision, intent(in) :: X(3), Xcenter(3), Qg2l(3,3)
  double precision :: Xloc(2), Xloc3(3)

  Xloc3 = matmul(Qg2l, X - Xcenter)
  Xloc = Xloc3(1:2)
  if (abs(Xloc3(3)) > 1.0e-7) stop 'abs(Xloc3(3)) > 1.0e-7'

end

function calc_mu(X0hat, b) result(mu)
  double precision, intent(in) :: X0hat(2), b
  double precision :: mu(2)

  mu(1) = 0.5*(asinh((1.0 + X0hat(1))/b) + asinh((1.0 - X0hat(1))/b))
  mu(2) = 0.5*(asinh((1.0 + X0hat(2))/b) + asinh((1.0 - X0hat(2))/b))

end

function calc_eta(X0hat, b) result(eta)
  double precision, intent(in) :: X0hat(2), b
  double precision :: eta(2)

  eta(1) = 0.5*(asinh((1.0 + X0hat(1))/b) - asinh((1.0 - X0hat(1))/b))
  eta(2) = 0.5*(asinh((1.0 + X0hat(2))/b) - asinh((1.0 - X0hat(2))/b))

end

subroutine test_Green_function_integration_nearly_sing(intg_pt_quad)
  use global
  use types
  use tensor_functions
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)

  double precision :: x1(3), x2(3), x3(3), x4(3), G_Gauss(3,3), Gintg(3,3), dx(3)
  double precision :: xintpt(3,nintg_pt), fintpt(nintg_pt), jacob, x1pt, x2pt, x3pt, x123(3,3)
  double precision :: normal(3), x0proj(3), x0(3), xcenter(3), x0proj_loc(2), Qg2l(3,3), x0hat(2)
  double precision :: b_param, eta_param(2), mu_param(2), xbar(2), xhat(2), x(3), wtot, w, area
  double precision :: x1loc(2), x2loc(2), x3loc(2), x4loc(2), jacob_sinh, jacob_ref2quad, oneoverx
  double precision :: x1start, x1end, x2start, x2end, area1, GxR(3,3), dx1, dx2, dx_1d
  integer :: i, ip1, ip2, ip3, np1, np2, np3, iintg, ip
  character (len=10) :: str

  x1 = [1.0,5.0,-3.0]
  x2 = [10.0,4.0,3.0]
  x3 = [-3.0,4.0,2.0]
  x4 = x1 + (x3 - x2 - 0.3*(x2-x1))

  ! x1 = [0.0,0.0,0.0]
  ! x2 = [5.0,0.0,0.0]
  ! x3 = [5.0,5.0,0.0]
  ! x4 = [0.0,5.0,0.0]
  
  x0 = [4.0,4.0,0.05]

  ! project the point
  normal = cross_prod(x2 - x1, x4 - x1)
  normal = normal/norm2(normal)
  x0proj = point2plane(x1, normal, x0)

  ! parameters for sinh transform
  call transformation_mat_to_local_quad(x1, x2, x3, x4, xcenter, Qg2l)
  x0proj_loc = transform_to_local_quad(x0proj, xcenter, Qg2l)
  x1loc = transform_to_local_quad(x1, xcenter, Qg2l)
  x2loc = transform_to_local_quad(x2, xcenter, Qg2l)
  x3loc = transform_to_local_quad(x3, xcenter, Qg2l)
  x4loc = transform_to_local_quad(x4, xcenter, Qg2l)
  x0hat = quad_to_reference_2D(x0proj_loc, x1loc, x2loc, x3loc, x4loc)
  b_param = norm2(x0proj - x0)
  mu_param = calc_mu(x0hat, b_param)
  eta_param = calc_eta(x0hat, b_param)

  ! ordering of nodes
  ! write(*,*)x1loc
  ! write(*,*)x2loc
  ! write(*,*)x3loc
  ! write(*,*)x4loc
  ! write(*,*)reference_to_quad_2D([-1.0,-1.0],x1loc, x2loc, x3loc, x4loc)
  ! write(*,*)reference_to_quad_2D([+1.0,-1.0],x1loc, x2loc, x3loc, x4loc)
  ! write(*,*)reference_to_quad_2D([+1.0,+1.0],x1loc, x2loc, x3loc, x4loc)
  ! write(*,*)reference_to_quad_2D([-1.0,+1.0],x1loc, x2loc, x3loc, x4loc)

  write(*,*) 'x0proj', x0proj
  write(*,*) 'x0proj_loc', x0proj_loc
  write(*,*) 'x0hat', x0hat
  write(*,*) 'b_param', b_param
  write(*,*) 'mu_param', mu_param
  write(*,*) 'eta_param', eta_param

  open (unit = 33, file = 'plot_quad.py', status = 'unknown')
  write(33,'(A)') 'import matplotlib.pyplot as plt'
  write(33,'(A)') 'import numpy as np'
  write(33,'(A)') 'import os'
  write(33,'(A)') 'import re'
  write(33,'(A)')
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'x0 = np.array([',x0(1),',',x0(2),',',x0(3),'])'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'x0proj = np.array([',x0proj(1),',',x0proj(2),',',x0proj(3),'])'
  write(33,'(A)') 'quad_x = np.empty([3,5])'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad_x[:,0] = [',x1(1),',',x1(2),',',x1(3),']'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad_x[:,1] = [',x2(1),',',x2(2),',',x2(3),']'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad_x[:,2] = [',x3(1),',',x3(2),',',x3(3),']'
  write(33,'(A,E16.8,A,E16.8,A,E16.8,A)') 'quad_x[:,3] = [',x4(1),',',x4(2),',',x4(3),']'
  write(33,'(A)') 'quad_x[:,4] = quad_x[:,0]'
  write(33,'(A,I2,A)') 'xintgpt = np.empty([3,',nintg_pt_quad,'])'

  G_Gauss = 0.0
  wtot = 0.0
  oneoverx = 0.0
  do ip = 1, nintg_pt_quad

    xbar = intg_pt_quad(ip)%Xhat
    xhat = sinh_transform(xbar, x0hat, b_param, mu_param, eta_param)
    x = reference_to_quad_3D(xhat, x1, x2, x3, x4)

    write(33,'(A,I2,A,E16.8,A,E16.8,A,E16.8,A)')'xintgpt[:,',ip-1,'] = np.array([',x(1),',',x(2),',',x(3),'])'

    jacob_sinh = sinh_transform_jacob(xbar, x0hat, b_param, mu_param, eta_param)
    jacob_ref2quad = reference_to_quad_2D_jacob(xhat, x1loc, x2loc, x3loc, x4loc)
    w = intg_pt_quad(ip)%w*abs(jacob_sinh)*abs(jacob_ref2quad)
    wtot = wtot + w

    ! write(*,*) 'jacob_ref2quad',jacob_ref2quad
    ! write(*,*) 'jacob_sinh',jacob_sinh
    ! write(*,*) 'w',w

    dx = x0 - x 
    oneoverx = oneoverx + 1.0/norm2(dx)*w
    GxR = Green_function_x_R(dx, lambda, mu)
    G_Gauss = G_Gauss + GxR/norm2(dx)*w

  enddo
  x123(:,1) = x1
  x123(:,2) = X2
  x123(:,3) = x3
  area = triangle_area(x123)
  x123(:,1) = x1
  x123(:,2) = x3
  x123(:,3) = x4
  area = area + triangle_area(x123)
  write(*,*) wtot, area
  write(*,*) 'oneoverx',oneoverx

  write(33,'(A)') 'fig = plt.figure()'
  write(33,'(A)') 'ax = fig.add_subplot(111, projection=''3d'')'
  write(33,'(A)') 'ax.plot(quad_x[0,:], quad_x[1,:], quad_x[2,:])'
  write(33,'(A)') 'ax.scatter(x0[0], x0[1], x0[2])'
  write(33,'(A)') 'ax.scatter(x0proj[0], x0proj[1], x0proj[2])'
  write(33,'(A)') 'ax.scatter(xintgpt[0], xintgpt[1], xintgpt[2])'

  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)
  write(*,*) 'G (nearly singular sinh transform)'
  G_Gauss = integrate_Green_function_nearly_sing_quad(x0, x1, x2, x3, x4, intg_pt_quad, lambda, mu)
  write(*,'(3E16.8)') (G_Gauss(i,:),i=1,3)

  ! in local frame
  write(33,'(A)') 'xintgpt_reg = np.array(['
  dx_1d = 0.05
  x1start = minval([x1loc(1),x2loc(1),x3loc(1),x4loc(1)]) - dx_1d
  x1end = maxval([x1loc(1),x2loc(1),x3loc(1),x4loc(1)]) + dx_1d
  x2start = minval([x1loc(2),x2loc(2),x3loc(2),x4loc(2)]) - dx_1d
  x2end = maxval([x1loc(2),x2loc(2),x3loc(2),x4loc(2)]) + dx_1d
  np1 = nint((x1end-x1start)/dx_1d)
  np2 = nint((x2end-x2start)/dx_1d)
  write(*,*)'x1start,x1end',x1start,x1end
  write(*,*)'x2start,x2end',x2start,x2end
  dx1 = (x1end-x1start)/float(np1)
  dx2 = (x2end-x2start)/float(np2)
  x1pt = x1start
  Gintg = 0.0
  do ip1 = 1,np1
    x1pt = x1pt + dx1
    x2pt = x2start
    do ip2 = 1,np2
      x2pt = x2pt + dx2
      ! write(*,*) x1pt, x2pt
      x(1) = x1pt
      x(2) = x2pt
      x(3) = 0.0

      x = matmul(transpose(Qg2l), x) + xcenter

      area1 = 0.0
      x123(:,1) = x
      x123(:,2) = x1
      x123(:,3) = x2
      area1 = area1 + triangle_area(x123)
      x123(:,1) = x
      x123(:,2) = x2
      x123(:,3) = x3
      area1 = area1 + triangle_area(x123)
      x123(:,1) = x
      x123(:,2) = x3
      x123(:,3) = x4
      area1 = area1 + triangle_area(x123)
      x123(:,1) = x
      x123(:,2) = x4
      x123(:,3) = x1
      area1 = area1 + triangle_area(x123)

      if (abs(area1 - area) < 1.0e-5) then

        ! write(*,*)x
        write(33,'(A,E16.8,A,E16.8,A,E16.8,A)')'[',x(1),',',x(2),',',x(3),'],'

        dx = x0 - x 
        Gintg = Gintg + Green_function(dx, lambda, mu)*dx1*dx2
      endif
    enddo
  enddo
  write(33,'(A)') '])'
  write(33,'(A)') 'ax.scatter(xintgpt_reg[:,0], xintgpt_reg[:,1], xintgpt_reg[:,2])'
  write(33,'(A)') 'plt.show()'
  write(33,'(A)') 'ax.axis(''equal'')'
  close(33)

  write(*,*) 'G (regular grid)'
  write(*,'(3E16.8)') (Gintg(i,:),i=1,3)

end

function integrate_Green_function_nearly_sing_quad(x0, x1, x2, x3, x4, intg_pt_quad, lambda, mu) result(Gintg) 
  use types
  use tensor_functions
  use global, only : nintg_pt_quad
  type (intg_pt_type), intent(in) :: intg_pt_quad(nintg_pt_quad)
  double precision, intent(in) :: x0(3), x1(3), x2(3), x3(3), x4(3), lambda, mu
  double precision :: normal(3), x0proj(3), xcenter(3), Qg2l(3,3), x0proj_loc(2) 
  double precision :: x1loc(2), x2loc(2), x3loc(2), x4loc(2), x0hat(2)
  double precision :: b_param, eta_param(2), mu_param(2), xbar(2), xhat(2), x(3), wtot, w, area
  double precision :: Gintg(3,3), jacob_sinh, jacob_ref2quad, dx(3), GxR(3,3), x123(3,3)
  integer :: ip
  logical :: sinh_trans

  ! project the point
  normal = cross_prod(x2 - x1, x4 - x1)
  normal = normal/norm2(normal)
  x0proj = point2plane(x1, normal, x0)

  ! parameters for sinh transform
  call transformation_mat_to_local_quad(x1, x2, x3, x4, xcenter, Qg2l)
  x0proj_loc = transform_to_local_quad(x0proj, xcenter, Qg2l)
  x1loc = transform_to_local_quad(x1, xcenter, Qg2l)
  x2loc = transform_to_local_quad(x2, xcenter, Qg2l)
  x3loc = transform_to_local_quad(x3, xcenter, Qg2l)
  x4loc = transform_to_local_quad(x4, xcenter, Qg2l)
  x0hat = quad_to_reference_2D(x0proj_loc, x1loc, x2loc, x3loc, x4loc)
  if (abs(x0hat(1) - (-100.0)) <= 1.0e-7 .and. abs(x0hat(2) - (-100.0)) <= 1.0e-7) then
    sinh_trans = .false.
  else
    sinh_trans = .true.
  endif
  b_param = norm2(x0proj - x0)
  if (b_param <= 1.0e-7) then
    ! write(*,*) x0proj, x0
    ! write(*,*) x1
    ! write(*,*) x2
    ! write(*,*) x3
    ! write(*,*) x4
    ! stop
    sinh_trans = .false.
  else
    sinh_trans = .true.
  endif
  if (sinh_trans) mu_param = calc_mu(x0hat, b_param)
  if (sinh_trans) eta_param = calc_eta(x0hat, b_param)

  Gintg = 0.0
  wtot = 0.0
  do ip = 1, nintg_pt_quad

    xbar = intg_pt_quad(ip)%Xhat
    if (sinh_trans) then
      xhat = sinh_transform(xbar, x0hat, b_param, mu_param, eta_param)
    else
      xhat = xbar
    endif
    x = reference_to_quad_3D(xhat, x1, x2, x3, x4)

    if (sinh_trans) then
      jacob_sinh = sinh_transform_jacob(xbar, x0hat, b_param, mu_param, eta_param)
    else
      jacob_sinh = 1.0
    endif
    jacob_ref2quad = reference_to_quad_2D_jacob(xhat, x1loc, x2loc, x3loc, x4loc)
    w = intg_pt_quad(ip)%w*abs(jacob_sinh)*abs(jacob_ref2quad)
    wtot = wtot + w

    dx = x0 - x 
    GxR = Green_function_x_R(dx, lambda, mu)
    Gintg = Gintg + GxR/norm2(dx)*w

  enddo

  ! x123(:,1) = x1
  ! x123(:,2) = X2
  ! x123(:,3) = x3
  ! if (abs(wtot - triangle_area(x123)) > 1.0e-2) then
  !   stop 'abs(wtot - triangle_area(x123)) > 1.0e-2'
  ! endif

end

subroutine analytical_sol_thermal_exp(elements,grains)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (grain_type), allocatable, intent (inout) :: grains(:)

  double precision :: xc(3), x(3), R, theta, phi, T(3,3), nu
  double precision :: strain_radial, strain_tang, a, V, eigen, strain_calc(6)
  double precision :: strain_loc(3,3), err_strain, avg_strain
  integer :: iel, igrain_incl

  igrain_incl = 2
  eigen = grains(igrain_incl)%eigenstrain(1,1)
  V = 0.0
  do iel = 1,nel
    if (norm2(elements(iel)%eigen) > 0.0) V = V + elements(iel)%v
  enddo
  ! v = 4.0/3.0*pi*a^3
  a = (3.0/4.0*v/pi)**0.333333333333333
  write(*,*) 'calc a', a
  ! write(*,*) 'input a'
  ! read(*,*) a

  xc(1) = 0.5 + float(npts1)*0.5
  xc(2) = 0.5 + float(npts2)*0.5
  xc(3) = 0.5 + float(npts3)*0.5
  nu = lambda/2.0/(lambda + mu)
  err_strain = 0.0
  avg_strain = 0.0
  do iel = 1,nel

    ! Vectors and Tensor Operations in Polar Coordinates, Continuum Mechanics - Polar Coordinates, Alan Bower
    ! to spherical
    x = elements(iel)%x - xc
    R = norm2(x)
    theta = acos(x(3)/R)
    phi = atan(x(2)/x(1))

    ! transformation matrix
    T(1,1) = sin(theta)*cos(phi)
    T(2,1) = sin(theta)*sin(phi)
    T(3,1) = cos(theta)
    T(1,2) = cos(theta)*cos(phi)
    T(2,2) = cos(theta)*sin(phi)
    T(3,2) = -sin(theta)
    T(1,3) = -sin(phi)
    T(2,3) = cos(phi)
    T(3,3) = 0.0

    ! strain (Mura, Eqs. (11.44-11.45))
    if (elements(iel)%grain_id == igrain_incl) then
      strain_radial = (1.0 + nu)/3.0/(1.0 - nu)*eigen
      strain_tang = (1.0 + nu)/3.0/(1.0 - nu)*eigen
    else
      strain_radial = -2.0/3.0*(1 + nu)/(1.0 - nu)*(a/R)**3*eigen
      strain_tang = 1.0/3.0*(1 + nu)/(1.0 - nu)*(a/R)**3*eigen
    endif

    ! write(*,*) 'x', x
    ! write(*,*) 'R,theta,phi', R,theta,phi
    ! write(*,*) 'strain_radial', strain_radial
    ! write(*,*) 'strain_tang', strain_tang
    ! write(*,*) 'nu', nu

    ! Applied Mechanics of Solids (A.F. Bower) Chapter 4
    ! 4.1.2 Simplified equations for spherically symmetric linear elasticity problems
    strain_loc = 0.0
    strain_loc(1,1) = strain_radial
    strain_loc(2,2) = strain_tang
    strain_loc(3,3) = strain_tang

    strain_calc = elements(iel)%strain
    elements(iel)%strain = mandel_t2_to_t1(matmul(T, matmul(strain_loc, transpose(T))))
    elements(iel)%stress = matmul(elements(iel)%C, elements(iel)%strain - elements(iel)%eigen)
    elements(iel)%polar = 0.0
    elements(iel)%polarold = 0.0

    err_strain = err_strain + norm2(strain_calc - elements(iel)%strain)
    avg_strain = avg_strain + norm2(elements(iel)%strain)
  enddo
  err_strain = err_strain/float(nel)
  avg_strain = avg_strain/float(nel)
  write(*,*) 'strain difference', err_strain, err_strain/avg_strain

end

subroutine total_load(boxes, faces, edges)
  use global
  use types
  use tensor_functions
  type (box_type), allocatable, intent (inout) :: boxes(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)

  double precision :: tot, conv_edg, conv_fc1, conv_fc2
  integer :: ib1, j, iedg1, i, ib2, k, iedg2, ifc

  conv_edg = 0.0
  conv_fc1 = 0.0
  conv_fc2 = 0.0

  do ib1 = 1,nbox

    do j = 1,boxes(ib1)%nedg ! these are all including shared

      iedg1 = boxes(ib1)%edge_ids(j)
      do i = 1,boxes(ib1)%nnear
        ib2 = boxes(ib1)%near_box_id(i)

        do k = 1,boxes(ib2)%nedg
          iedg2 = boxes(ib2)%edge_ids(k)
        
          ! full convolution over near edges
          if (iedg1 .ne. iedg2) then
            conv_edg = conv_edg + 1.0
          endif
        enddo
    
      enddo

      do i = 1, edges(iedg1)%nface ! loop over faces containing this edge
        ifc = edges(iedg1)%face(i)
    
        ! remove effects of faces containing this edge
        do k = 1, 3 ! loop over edges in face
          iedg2 = faces(ifc)%edge(k)
          if (iedg2 .ne. iedg1) then
            conv_fc1 = conv_fc1 + 1
          endif
        enddo
    
        ! add true effect of faces containing this edge
        conv_fc1 = conv_fc1 + 1
    
      enddo

      if (allocated(edges(iedg1)%face_close)) then
        do i = 1, size(edges(iedg1)%face_close) ! loop over close faces 
          ifc = edges(iedg1)%face_close(i)
      
          ! remove effects of close faces
          do k = 1, 3 ! loop over edges in face
            iedg2 = faces(ifc)%edge(k)
            if (iedg2 .ne. iedg1) then
              conv_fc2 = conv_fc2 + 1.0
            endif
          enddo
      
          ! add true effect of faces containing this edge
          conv_fc2 = conv_fc2 + 1.0
      
        enddo
      endif
    enddo

  enddo

  tot = (conv_edg + conv_fc1 + conv_fc2)
  write(*,*) conv_edg, conv_fc1, conv_fc2
  write(*,*) conv_edg/tot,conv_fc1/tot, conv_fc2/tot

end

subroutine calc_near_f_box(boxes, edges)
  use global
  use types
  use tensor_functions
  type (box_type), allocatable, intent (inout) :: boxes(:)
  type (edge_type), allocatable, intent (inout) :: edges(:)

  integer :: ib1, iedg_near, i, iedg2, ind, ib2
  double precision :: f(3)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(ib1, iedg_near, i, iedg2, ind, ib2, f) &
  !$OMP SHARED(nbox, boxes, edges)
  !$OMP DO 
  do ib1 = 1,nbox
    do iedg_near = 1, size(boxes(ib1)%near_edge_ids)
      iedg2 = boxes(ib1)%near_edge_ids(iedg_near)

      f = 0.0
      do i = 1, size(edges(iedg2)%box_id)
        ib2 = edges(iedg2)%box_id(i)
        if (findloc(boxes(ib1)%near_box_id, ib2, 1) .ne. 0) then ! in near boxes
          ind = edges(iedg2)%box_id_loc(i)
          f = f + boxes(ib2)%edge_f(:,ind)
        endif
      enddo
      boxes(ib1)%near_edge_f(:,iedg_near) = f
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
end

subroutine init_map_outgoing_exp_arr_to_1d
  use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
  use global, only : arr1d_2d_fact, arr1d_3d_fact, arr1d_4d_fact

  implicit none

  integer :: i, j, k, l, kk, ic, dum2(2), dum3(3), dum4(4), n1n2n3dum(3)
  integer, allocatable :: ind1d(:,:), n1n2n3(:,:)
  logical :: flag

  write(*,*) 'ind2d21d'
  ! ind2d
  ! 1d array
  ic = 0
  do i = 1,3
    do j = i,3
      ic = ic + 1
    enddo
  enddo
  allocate(ind1d(2,ic))
  allocate(n1n2n3(3,ic))
  ic = 0
  do i = 1,3
    do j = i,3
      ic = ic + 1
      ind1d(:,ic) = [i,j] ! indices in 1d array
      n1n2n3(:,ic) = 0
      do kk = 1,2
        if (ind1d(kk,ic) == 1) n1n2n3(1,ic) = n1n2n3(1,ic) + 1 ! number of 1s
        if (ind1d(kk,ic) == 2) n1n2n3(2,ic) = n1n2n3(2,ic) + 1 ! number of 2s
        if (ind1d(kk,ic) == 3) n1n2n3(3,ic) = n1n2n3(3,ic) + 1 ! number of 3s
      enddo
      write(*,*) 'ind1d', ind1d(:,ic)
      write(*,*) 'n1n2n3', n1n2n3(:,ic)
    enddo
  enddo

  do i = 1,3
    do j = 1,3
      dum2 = [i,j]
      n1n2n3dum = 0
      do kk = 1,2
        if (dum2(kk) == 1) n1n2n3dum(1) = n1n2n3dum(1) + 1 ! number of 1s
        if (dum2(kk) == 2) n1n2n3dum(2) = n1n2n3dum(2) + 1 ! number of 2s
        if (dum2(kk) == 3) n1n2n3dum(3) = n1n2n3dum(3) + 1 ! number of 3s
      enddo
      flag = .false.
      do kk = 1, ic
        if(n1n2n3dum(1) == n1n2n3(1,kk) .and. n1n2n3dum(2) == n1n2n3(2,kk) .and. n1n2n3dum(3) == n1n2n3(3,kk)) then
          ind2d21d(i,j) = kk
          flag = .true.
          exit
        endif
      enddo
      if (.not.flag) stop 'location not found'
      write(*,*) i,j,ind1d(:,ind2d21d(i,j)), ind2d21d(i,j)
    enddo
  enddo
  ind1d22d = ind1d
  deallocate(ind1d)
  deallocate(n1n2n3)

  write(*,*) 'ind3d21d'
  ! ind3d21d
  ! 1d array
  ic = 0
  do i = 1,3
    do j = i,3
      do k = j,3
        ic = ic + 1
      enddo
    enddo
  enddo
  allocate(ind1d(3,ic))
  allocate(n1n2n3(3,ic))
  ic = 0
  do i = 1,3
    do j = i,3
      do k = j,3
        ic = ic + 1
        ind1d(:,ic) = [i,j,k] ! indices in 1d array
        n1n2n3(:,ic) = 0
        do kk = 1,3
          if (ind1d(kk,ic) == 1) n1n2n3(1,ic) = n1n2n3(1,ic) + 1 ! number of 1s
          if (ind1d(kk,ic) == 2) n1n2n3(2,ic) = n1n2n3(2,ic) + 1 ! number of 2s
          if (ind1d(kk,ic) == 3) n1n2n3(3,ic) = n1n2n3(3,ic) + 1 ! number of 3s
        enddo
        write(*,*) 'ind1d', ind1d(:,ic)
        write(*,*) 'n1n2n3', n1n2n3(:,ic)
      enddo
    enddo
  enddo

  do i = 1,3
    do j = 1,3
      do k = 1,3
        dum3 = [i,j,k]
        n1n2n3dum = 0
        do kk = 1,3
          if (dum3(kk) == 1) n1n2n3dum(1) = n1n2n3dum(1) + 1 ! number of 1s
          if (dum3(kk) == 2) n1n2n3dum(2) = n1n2n3dum(2) + 1 ! number of 2s
          if (dum3(kk) == 3) n1n2n3dum(3) = n1n2n3dum(3) + 1 ! number of 3s
        enddo
        flag = .false.
        do kk = 1, ic
          if(n1n2n3dum(1) == n1n2n3(1,kk) .and. n1n2n3dum(2) == n1n2n3(2,kk) .and. n1n2n3dum(3) == n1n2n3(3,kk)) then
            ind3d21d(i,j,k) = kk
            flag = .true.
            exit
          endif
        enddo
        if (.not.flag) stop 'location not found'
        write(*,*) i,j,k,ind1d(:,ind3d21d(i,j,k)), ind3d21d(i,j,k)
      enddo
    enddo
  enddo
  ind1d23d = ind1d
  deallocate(ind1d)
  deallocate(n1n2n3)

  write(*,*) 'ind4d21d'
  ! ind4d21d
  ! 1d array
  ic = 0
  do i = 1,3
    do j = i,3
      do k = j,3
        do l = k,3
          ic = ic + 1
        enddo
      enddo
    enddo
  enddo
  allocate(ind1d(4,ic))
  allocate(n1n2n3(3,ic))
  ic = 0
  do i = 1,3
    do j = i,3
      do k = j,3
        do l = k,3
          ic = ic + 1
          ind1d(:,ic) = [i,j,k,l] ! indices in 1d array
          n1n2n3(:,ic) = 0
          do kk = 1,4
            if (ind1d(kk,ic) == 1) n1n2n3(1,ic) = n1n2n3(1,ic) + 1 ! number of 1s
            if (ind1d(kk,ic) == 2) n1n2n3(2,ic) = n1n2n3(2,ic) + 1 ! number of 2s
            if (ind1d(kk,ic) == 3) n1n2n3(3,ic) = n1n2n3(3,ic) + 1 ! number of 3s
          enddo
          write(*,*) 'ind1d', ind1d(:,ic)
          write(*,*) 'n1n2n3', n1n2n3(:,ic)
        enddo
      enddo
    enddo
  enddo

  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          dum4 = [i,j,k,l]
          n1n2n3dum = 0
          do kk = 1,4
            if (dum4(kk) == 1) n1n2n3dum(1) = n1n2n3dum(1) + 1 ! number of 1s
            if (dum4(kk) == 2) n1n2n3dum(2) = n1n2n3dum(2) + 1 ! number of 2s
            if (dum4(kk) == 3) n1n2n3dum(3) = n1n2n3dum(3) + 1 ! number of 3s
          enddo
          flag = .false.
          do kk = 1, ic
            if(n1n2n3dum(1) == n1n2n3(1,kk) .and. n1n2n3dum(2) == n1n2n3(2,kk) .and. n1n2n3dum(3) == n1n2n3(3,kk)) then
              ind4d21d(i,j,k,l) = kk
              flag = .true.
              exit
            endif
          enddo
          if (.not.flag) stop 'location not found'
          write(*,*) i,j,k,l,ind1d(:,ind4d21d(i,j,k,l)), ind4d21d(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  ind1d24d = ind1d
  deallocate(ind1d)
  deallocate(n1n2n3)

  ! factors
  arr1d_2d_fact = 0.0
  arr1d_3d_fact = 0.0
  arr1d_4d_fact = 0.0
  do i = 1,3
    do j = 1,3
      arr1d_2d_fact(ind2d21d(i,j)) = arr1d_2d_fact(ind2d21d(i,j)) + 1.0
      do k = 1,3
        arr1d_3d_fact(ind3d21d(i,j,k)) = arr1d_3d_fact(ind3d21d(i,j,k)) + 1.0
        do l = 1,3
          arr1d_4d_fact(ind4d21d(i,j,k,l)) = arr1d_4d_fact(ind4d21d(i,j,k,l)) + 1.0
        enddo
      enddo
    enddo
  enddo  
  write(*,*) arr1d_2d_fact

end

subroutine test_Mandel
  use tensor_functions
  implicit none
  double precision t33(3,3), t6(6), t3333(3,3,3,3), t66(6,6)
  double precision t33tmp(3,3), t3333tmp(3,3,3,3)
  integer :: i,j,k,l,ic

  ic = 0
  do i = 1,3
    do j = 1,3
      ic = ic + 1
      t33(i,j) = ic
    enddo
  enddo
  t33 = sym33(t33)
  do i = 1,3
    write(*,*)t33(i,:)
  enddo
  t6 = mandel_t2_to_t1(t33)
  write(*,*)t6
  t33tmp = mandel_t1_to_t2(t6)
  do i = 1,3
    write(*,*)t33tmp(i,:)
  enddo
  write(*,*)norm2(t33-t33tmp)

  ic = 0
  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          ic = ic + 1
          t3333(i,j,k,l) = ic
          write(*,*) i,j,k,l,t3333(i,j,k,l)
        enddo
      enddo
    enddo
  enddo

  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          t3333tmp(i,j,k,l) = (t3333(i,j,k,l) + t3333(j,i,k,l) + t3333(i,j,l,k) + t3333(j,i,l,k))/4.0
          write(*,*) i,j,k,l,t3333tmp(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  t3333 = t3333tmp
  t66 = mandel_t4_to_t2(t3333)
  do i = 1,6
    write(*,*)t66(i,:)
  enddo
  t3333tmp = mandel_t2_to_t4(t66)
  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          write(*,*) i, j, k, l, t3333tmp(i,j,k,l), t3333(i,j,k,l)
        enddo
      enddo
    enddo
  enddo
  write(*,*)norm2(t3333-t3333tmp)

end

subroutine calc_stress_polar_elastic(elements)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: iel1
  double precision :: strain_from_stress(6)

  ! stress and polarization
  polaravg = 0.0
  stressavg = 0.0
  vtot = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel1,strain_from_stress) &
  !$OMP SHARED(elements,nel,solution_method,C066) &
  !$OMP REDUCTION(+:polaravg,stressavg,vtot)
  !$OMP DO 
  do iel1 = 1,nel
    elements(iel1)%polarold = elements(iel1)%polar
    if (solution_method == 1) then
      elements(iel1)%stress = matmul(elements(iel1)%C, elements(iel1)%strain - elements(iel1)%eigen)
      elements(iel1)%polar = elements(iel1)%stress - matmul(c066, elements(iel1)%strain - elements(iel1)%straint)
    else
      elements(iel1)%stress = matmul(elements(iel1)%I_C0S_inv,elements(iel1)%stress + &
        matmul(C066, elements(iel1)%strain - elements(iel1)%eigen))
      strain_from_stress = matmul(elements(iel1)%S, elements(iel1)%stress) + elements(iel1)%eigen
      elements(iel1)%polar = elements(iel1)%stress - matmul(c066, strain_from_stress)
    endif
    
    polaravg = polaravg + elements(iel1)%polar*elements(iel1)%v
    stressavg = stressavg + elements(iel1)%stress*elements(iel1)%v
    vtot = vtot + elements(iel1)%v
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  polaravg = polaravg/vtot
  stressavg = stressavg/vtot
  err_polar = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel1) &
  !$OMP SHARED(elements,nel,polaravg) &
  !$OMP REDUCTION(+:err_polar)
  !$OMP DO 
  do iel1 = 1,nel
    elements(iel1)%polar = elements(iel1)%polar - polaravg
    err_polar = err_polar + sum((elements(iel1)%polarold - elements(iel1)%polar)**2)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  err_polar = sqrt(err_polar)

end

subroutine init_slip_fcc
  use global, only : slip_normal, slip_dir, Schmid, nsys

  integer :: i, j, k

  nsys = 12
  allocate(slip_normal(3,nsys))
  allocate(slip_dir(3,nsys))
  allocate(Schmid(3,3,nsys))

  slip_normal(:,1) = [1, 1,-1]
  slip_normal(:,2) = [1, 1,-1]
  slip_normal(:,3) = [1, 1,-1]
  slip_normal(:,4) = [1,-1,-1]
  slip_normal(:,5) = [1,-1,-1]
  slip_normal(:,6) = [1,-1,-1]
  slip_normal(:,7) = [1,-1, 1]
  slip_normal(:,8) = [1,-1, 1]
  slip_normal(:,9) = [1,-1, 1]
  slip_normal(:,10) = [1, 1, 1]
  slip_normal(:,11) = [1, 1, 1]
  slip_normal(:,12) = [1, 1, 1]
  slip_dir(:,1) = [0, 1, 1]
  slip_dir(:,2) = [1, 0, 1]
  slip_dir(:,3) = [1,-1, 0]
  slip_dir(:,4) = [0, 1,-1]
  slip_dir(:,5) = [1, 0, 1]
  slip_dir(:,6) = [1, 1, 0]
  slip_dir(:,7) = [0, 1, 1]
  slip_dir(:,8) = [1, 0,-1]
  slip_dir(:,9) = [1, 1, 0]
  slip_dir(:,10) = [0, 1,-1]
  slip_dir(:,11) = [1, 0,-1]
  slip_dir(:,12) = [1,-1, 0]

  do i = 1,nsys
    slip_normal(:,i) = slip_normal(:,i)/norm2(slip_normal(:,i))
    slip_dir(:,i) = slip_dir(:,i)/norm2(slip_dir(:,i))
    do j = 1,3
      do k = 1,3
        Schmid(j,k,i) = slip_dir(j,i)*slip_normal(k,i)
      enddo
    enddo
    write(*,*) slip_normal(:,i), slip_dir(:,i)
    ! write(*,*) Schmid(:,:,i)
  enddo

end

subroutine allocate_slip_arr(elements)
  use types
  use global, only : nsys, nel
  type (element_type), allocatable, intent (inout) :: elements(:)
  integer :: iel, i

  do iel = 1,nel
    allocate(elements(iel)%Schmid(6,nsys))
    allocate(elements(iel)%tauc(nsys))
    allocate(elements(iel)%gamdot(nsys))
    elements(iel)%Schmid = 0.0
    elements(iel)%tauc = 0.0
    elements(iel)%gamdot = 0.0
  enddo

end

subroutine calc_Schmid(elements)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: iel, i
  double precision :: m(3,3)
  
  do iel = 1,nel

    do i = 1,nsys
      m = matmul(matmul(elements(iel)%Q, Schmid(:,:,i)), transpose(elements(iel)%Q))
      m = sym33(m)
      elements(iel)%Schmid(:,i) = mandel_t2_to_t1(m)
    enddo

  enddo

end

subroutine assign_tauc(elements)
  use global, only : tauc, nel
  use types
  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: iel
  
  do iel = 1,nel
    elements(iel)%tauc = tauc
  enddo

end

subroutine calc_plastic_rate(stress, Schmid, tauc, gamdot, edotp, dedotp_dstress)
  use global, only : nsys, gamdot0, rate_exp
  double precision, intent(in) :: stress(6), Schmid(6,nsys), tauc(nsys)
  double precision, intent(out) :: gamdot(nsys), edotp(6), dedotp_dstress(6,6)
  double precision :: tmp1(nsys), tmp2(nsys), rss(nsys), tmp3(nsys,6), stress_dev(6), pressure
  integer :: isys

  ! stress_dev = stress
  ! pressure = sum(stress_dev(1:3))/3.0
  ! stress_dev(1:3) = stress_dev(1:3) - pressure
  ! rss = matmul(stress_dev,Schmid)/tauc
  
  rss = matmul(stress,Schmid)/tauc

  tmp1 = gamdot0*abs(rss)**(rate_exp-1)
  tmp2 = rate_exp*tmp1/tauc
  gamdot = tmp1*rss
  
  edotp =  matmul(Schmid,gamdot)
  do isys = 1,nsys
    tmp3(isys,:) = tmp2(isys)*Schmid(:,isys)
  enddo
  dedotp_dstress = matmul(Schmid,tmp3)

end

subroutine test_plastic_rate(elements)
  use tensor_functions
  use global, only : nsys
  use types
  type (element_type), allocatable, intent (inout) :: elements(:)

  double precision :: stress(6), edotp(6), dedotp_dstress(6,6), gamdot(nsys)
  integer :: iel

  stress = [100.0, 100.0, - 200.0, 0.0, 0.0, 0.0]
  iel = 7871
  
  write(*,*) elements(iel)%grain_id
  write(*,*) mandel_t1_to_t2(elements(iel)%Schmid(:,1))
  write(*,*) elements(iel)%Q
  write(*,*) elements(iel)%tauc
  call calc_plastic_rate(stress, elements(iel)%Schmid, elements(iel)%tauc, gamdot, edotp, dedotp_dstress)
  write(*,*) mandel_t1_to_t2(edotp)
  write(*,*) mandel_t2_to_t4(dedotp_dstress)
 
end

subroutine calc_stress_polar_plastic_AL(elements)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: iel, iterNR, itmaxal = 100, icut
  double precision :: strain_from_stress(6), erral, res_prev_iter_n, resn, stress(6), stressold(6), edotp(6)
  double precision :: dedotp_dstress(6,6), strain_el(6), strain_pl(6), Schmid_loc(6,nsys), tauc_loc(nsys), stresst(6)
  double precision :: gamdot_loc(nsys), strain_pl_old(6), S_loc(6,6), eigen(6), strain(6), res(6), stress_prev_iter(6)
  double precision :: jacobian(6,6), stressn, erre, ddgnorm, errs, Gamma_self(6,6), polar(6), polaravgold(6), straint(6)
  double precision :: tolerance_AL = 1.0e-7

  ! stress and polarization
  polaravgold = polaravg
  polaravg = 0.0
  stressavg = 0.0
  vtot = 0.0
  erre = 0.0
  errs = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(strain_from_stress, erral, res_prev_iter_n, resn, stress, stressold, edotp) &
  !$OMP PRIVATE(dedotp_dstress, strain_el, strain_pl, Schmid_loc, tauc_loc) &
  !$OMP PRIVATE(gamdot_loc, strain_pl_old, S_loc, eigen, strain, res, stress_prev_iter) &
  !$OMP PRIVATE(jacobian, stressn, iel, iterNR, icut, Gamma_self, polar, straint, stresst) &
  !$OMP SHARED(tolerance_AL, elements, nel, itmaxal, time_inc, c066, iJacobi, polaravgold, id_gas) &
  !$OMP REDUCTION(+:polaravg, stressavg, vtot, erre, errs)
  !$OMP DO 
  do iel = 1,nel

    elements(iel)%polarold = elements(iel)%polar

    if (elements(iel)%grain_id .ne. id_gas) then

      Schmid_loc = elements(iel)%Schmid
      tauc_loc = elements(iel)%tauc
      S_loc = elements(iel)%S
      eigen = elements(iel)%eigen
      strain = elements(iel)%strain
      straint = elements(iel)%straint
      stresst = elements(iel)%stresst
      if (iJacobi == 1) Gamma_self = mandel_t4_to_t2(sym3333(elements(iel)%Gamma_self))
      if (iJacobi == 1) polar = elements(iel)%polar + polaravgold
  
      stress = elements(iel)%stress
      stressold = elements(iel)%stress
      strain_pl_old = elements(iel)%strain_pl
      stressn = norm2(stress)
  
      iterNR = 0
      erral = 2.0*tolerance_AL
      res_prev_iter_n = 1.0
  
      ! inner loop over stress
      do while(iterNR.lt.itmaxal.and.erral.gt.tolerance_AL)
  
        iterNR = iterNR + 1
  
        ! loop to find step length
        icut = 0
        resn = 2.0*res_prev_iter_n
        do while(resn.gt.res_prev_iter_n .and. icut.lt.10)
  
          ! constitutive response
          call calc_plastic_rate(stress, Schmid_loc, tauc_loc, gamdot_loc, edotp, dedotp_dstress)
          strain_from_stress = matmul(S_loc, stress) + strain_pl_old + edotp*time_inc + eigen
  
          ! residual
          res = stress - stressold - matmul(c066, strain - strain_from_stress)
  
          if (iJacobi == 1) res = res + matmul(c066, matmul(Gamma_self, -polar + (stress - stresst) - &
           matmul(C066,strain_from_stress - straint)))
  
          resn = norm2(res)
          if (iterNR.eq.1) res_prev_iter_n = resn*2.0
  
          if (resn.gt.res_prev_iter_n) stress = stress_prev_iter + (stress - stress_prev_iter)*0.5
  
          icut = icut + 1
  
        enddo ! end do while step length
        res_prev_iter_n = resn
        stress_prev_iter = stress
  
        ! Build Jacobian
        jacobian = id6 + matmul(c066, S_loc + dedotp_dstress*time_inc)
        if (iJacobi == 1) jacobian = jacobian + &
         matmul(c066, matmul(Gamma_self, id6 - matmul(C066, S_loc + dedotp_dstress*time_inc)))
  
        ! Calculate new stress by solving the system -[J][delt_sg] = [R]
        call lu_eqsystem(jacobian, res, 6)
        stress = stress - res
  
        ! errors
        erral = sqrt(sum((stress - stress_prev_iter)**2))
        if (stressn > 0.0) erral = erral/stressn
        ! write(*,*) stress
        ! write(*,*) resn
  
      enddo ! enddo for NR

      if(iterNR.ge.itmaxal) write(*,*) 'unconverged'
  
      ! ! test residual
      ! call calc_plastic_rate(stress, Schmid_loc, tauc_loc, gamdot_loc, edotp, dedotp_dstress)
      ! strain_from_stress = matmul(S_loc, stress) + strain_pl_old + edotp*time_inc + eigen
      ! write(*,*) stress - stressold - matmul(c066, strain - strain_from_stress)
      ! stop

    else
    ! if (elements(iel)%grain_id .eq. id_gas) then

      stress = 0.0
      gamdot_loc = 0.0
      edotp = 0.0
      strain = 0.0
      strain_from_stress = 0.0
      stressold = 0.0

    endif

    elements(iel)%stress = stress
    elements(iel)%gamdot = gamdot_loc
    elements(iel)%edotp = edotp
    elements(iel)%polar = (stress - stresst) - matmul(c066, strain_from_stress - elements(iel)%straint)
 
    erre = erre + sum((strain - strain_from_stress)**2)*elements(iel)%v
    errs = errs + sum((stress - stressold)**2)*elements(iel)%v

    polaravg = polaravg + elements(iel)%polar*elements(iel)%v
    stressavg = stressavg + elements(iel)%stress*elements(iel)%v
    vtot = vtot + elements(iel)%v
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  polaravg = polaravg/vtot
  stressavg = stressavg/vtot
  erre = sqrt(erre/vtot)/norm2(Eapp)
  errs = sqrt(errs/vtot)/norm2(stressavg)
  write(*,*) 'erre:', erre
  write(*,*) 'errs:', errs
  err_polar = 0.0
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel) &
  !$OMP SHARED(elements,nel,polaravg) &
  !$OMP REDUCTION(+:err_polar)
  !$OMP DO 
  do iel = 1,nel
    elements(iel)%polar = elements(iel)%polar - polaravg
    err_polar = err_polar + sum((elements(iel)%polarold - elements(iel)%polar)**2)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  err_polar = sqrt(err_polar)
  err = maxval([erre,errs,err])

end

subroutine update_state(elements)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)

  integer :: iel
  double precision :: dstrain(6), ddisgrad(3,3), dstress(6)

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(iel,dstrain,ddisgrad,dstress) &
  !$OMP SHARED(elements,nel,time_inc,incr,nincr)
  !$OMP DO 
  do iel = 1,nel
    elements(iel)%strain_pl = elements(iel)%strain_pl + elements(iel)%edotp*time_inc
    elements(iel)%strain_pl_vm = elements(iel)%strain_pl_vm + strain_vm(elements(iel)%edotp)*time_inc
    dstrain = elements(iel)%strain - elements(iel)%straint
    elements(iel)%straint = elements(iel)%strain
    if (incr < nincr) elements(iel)%strain = elements(iel)%strain + dstrain
    ddisgrad = elements(iel)%disgrad - elements(iel)%disgradt
    elements(iel)%disgradt = elements(iel)%disgrad
    if (incr < nincr) elements(iel)%disgrad = elements(iel)%disgrad + ddisgrad
    dstress = elements(iel)%stress - elements(iel)%stresst
    elements(iel)%stresst = elements(iel)%stress
    if (incr < nincr) elements(iel)%stress = elements(iel)%stress + dstress

  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end


subroutine non_local_hyd_correction(elements, faces, grains, boxes)
  use global
  use types
  use tensor_functions
  type (element_type), allocatable, intent (inout) :: elements(:)
  type (face_type), allocatable, intent (inout) :: faces(:)
  type (grain_type), allocatable, intent (inout) :: grains(:)
  type (box_type), allocatable, intent (inout) :: boxes(:)

  integer :: iel, itnonloc, i, ifc, ielngh, igr, ib, j
  double precision :: dum, dum1, strain_bx(6)

  do itnonloc = 1,nit_non_loc
    write(*,*) 'non local regularization iter:', itnonloc
    do iel = 1,nel
      dum = 0.0
      dum1 = 0.0
      do i = 1,4
        ifc = elements(iel)%face(i)
        if (faces(ifc)%nel == 2) then
          ielngh = faces(ifc)%el(1)
          if (ielngh == iel) ielngh = faces(ifc)%el(2)
          dum = dum + (elements(ielngh)%strain(1) - elements(ielngh)%straint(1) + &
           elements(ielngh)%strain(2) - elements(ielngh)%straint(2) + &
           elements(ielngh)%strain(3) - elements(ielngh)%straint(3))/3.0*&
           elements(ielngh)%v
          dum1 = dum1 + elements(ielngh)%v
        endif
      enddo
      dum = dum + (elements(iel)%strain(1) - elements(iel)%straint(1) + &
       elements(iel)%strain(2) - elements(iel)%straint(2) + &
       elements(iel)%strain(3) - elements(iel)%straint(3))/3.0*&
       elements(iel)%v
      dum1 = dum1 + elements(iel)%v

      elements(iel)%strain_hyd_non_loc = dum/dum1
    enddo

    strainavg = 0.0
    do iel = 1,nel
      dum = (elements(iel)%strain(1) - elements(iel)%straint(1) + &
             elements(iel)%strain(2) - elements(iel)%straint(2) + &
             elements(iel)%strain(3) - elements(iel)%straint(3))/3.0
      do i = 1,3
        elements(iel)%strain(i) = elements(iel)%strain(i) - dum + elements(iel)%strain_hyd_non_loc
        elements(iel)%disgrad(i,i) = elements(iel)%disgrad(i,i) - dum + elements(iel)%strain_hyd_non_loc
      enddo
      strainavg = strainavg + elements(iel)%strain*elements(iel)%v
      ! write(*,*)igr,sum(elements(iel)%strain(1:3))
    enddo
    strainavg = strainavg/vtot
  enddo
  ! stop

!  do itnonloc = 1,nit_non_loc
!    write(*,*) 'non local regularization iter:', itnonloc
!    do iel = 1,nel
!      dum = 0.0
!      dum1 = 0.0
!      do i = 1,4
!        ifc = elements(iel)%face(i)
!        if (faces(ifc)%nel == 2) then
!          ielngh = faces(ifc)%el(1)
!          if (ielngh == iel) ielngh = faces(ifc)%el(2)
!          dum = dum + ((elements(ielngh)%strain(1) - elements(ielngh)%straint(1) + 1.0)* &
!                       (elements(ielngh)%strain(2) - elements(ielngh)%straint(2) + 1.0)* &
!                       (elements(ielngh)%strain(3) - elements(ielngh)%straint(3) + 1.0))*&
!                      elements(ielngh)%v
!          dum1 = dum1 + elements(ielngh)%v
!        endif
!      enddo
!      dum = dum + ((elements(iel)%strain(1) - elements(iel)%straint(1) + 1.0)* &
!                   (elements(iel)%strain(2) - elements(iel)%straint(2) + 1.0)* &
!                   (elements(iel)%strain(3) - elements(iel)%straint(3) + 1.0))*&
!                  elements(iel)%v
!      dum1 = dum1 + elements(iel)%v
!
!      elements(iel)%strain_hyd_non_loc = dum/dum1
!    enddo
!
!    strainavg = 0.0
!    do iel = 1,nel
!      dum = ((elements(iel)%strain(1) - elements(iel)%straint(1) + 1.0)* &
!             (elements(iel)%strain(2) - elements(iel)%straint(2) + 1.0)* &
!             (elements(iel)%strain(3) - elements(iel)%straint(3) + 1.0))
!      ! do i = 1,3
!      !   elements(iel)%strain(i) = elements(iel)%strain(i) - dum + elements(iel)%strain_hyd_non_loc
!      !   elements(iel)%disgrad(i,i) = elements(iel)%disgrad(i,i) - dum + elements(iel)%strain_hyd_non_loc
!      ! enddo
!      elements(iel)%strain = elements(iel)%straint + &
!       (elements(iel)%strain_hyd_non_loc/dum)**(1.0/3.0)*(elements(iel)%strain - elements(iel)%straint + &
!       [1.0,1.0,1.0,0.0,0.0,0.0]) - [1.0,1.0,1.0,0.0,0.0,0.0]
!      elements(iel)%disgrad = elements(iel)%disgradt + &
!       (elements(iel)%strain_hyd_non_loc/dum)**(1.0/3.0)*(elements(iel)%disgrad - elements(iel)%disgradt + id3) - id3
!      strainavg = strainavg + elements(iel)%strain*elements(iel)%v
!      ! write(*,*)igr,sum(elements(iel)%strain(1:3))
!    enddo
!    strainavg = strainavg/vtot
!  enddo
!  ! stop

!   do igr = 1,ng
!     grains(igr)%strain = 0.0
!     grains(igr)%v = 0.0
!   enddo
!   do iel = 1,nel
!     igr = elements(iel)%grain_id
!     grains(igr)%strain = grains(igr)%strain + elements(iel)%strain*elements(iel)%v
!     grains(igr)%v = grains(igr)%v + elements(iel)%v
!   enddo
!   do igr = 1,ng
!     grains(igr)%strain = grains(igr)%strain/grains(igr)%v
!   enddo
!   strainavg = 0.0
!   do iel = 1,nel
!     igr = elements(iel)%grain_id
!     dum = - (elements(iel)%strain(1) + elements(iel)%strain(2) + elements(iel)%strain(3))/3.0 + &
!     (grains(igr)%strain(1) + grains(igr)%strain(2) + grains(igr)%strain(3))/3.0
!     do i = 1,3
!       elements(iel)%strain(i) = elements(iel)%strain(i) + dum
!     enddo
!     strainavg = strainavg + elements(iel)%strain*elements(iel)%v
!     ! write(*,*)igr,sum(elements(iel)%strain(1:3))
!   enddo
!   strainavg = strainavg/vtot

!  strainavg = 0.0
!  do ib = 1,nbox
!
!    ! write(*,*) ib
!
!    do igr = 1,ng
!      grains(igr)%strain = 0.0
!      grains(igr)%v = 0.0
!    enddo
!    do i = 1,boxes(ib)%nel
!      iel = boxes(ib)%element_ids(i)
!      igr = elements(iel)%grain_id
!      grains(igr)%strain = grains(igr)%strain + elements(iel)%strain*elements(iel)%v
!      grains(igr)%v = grains(igr)%v + elements(iel)%v
!    enddo
!    do igr = 1,ng
!      if (grains(igr)%v > 0) grains(igr)%strain = grains(igr)%strain/grains(igr)%v
!      ! if (grains(igr)%v > 0) write(*,*) igr, grains(igr)%strain
!    enddo
!
!    do i = 1,boxes(ib)%nel
!      iel = boxes(ib)%element_ids(i)
!      igr = elements(iel)%grain_id
!      dum = - (elements(iel)%strain(1) + elements(iel)%strain(2) + elements(iel)%strain(3))/3.0 + &
!      (grains(igr)%strain(1) + grains(igr)%strain(2) + grains(igr)%strain(3))/3.0
!      do j = 1,3
!        elements(iel)%strain(j) = elements(iel)%strain(j) + dum
!      enddo
!      ! write(*,*) iel, igr, sum(elements(iel)%strain(1:3))
!      strainavg = strainavg + elements(iel)%strain*elements(iel)%v
!    enddo
!
!  enddo
!  strainavg = strainavg/vtot

!  strainavg = 0.0
!  do ib = 1,nbox
!
!    ! write(*,*) ib
!
!    strain_bx = 0.0
!    dum = 0.0
!    do i = 1,boxes(ib)%nel
!      iel = boxes(ib)%element_ids(i)
!      strain_bx = strain_bx + elements(iel)%strain*elements(iel)%v
!      dum = dum + elements(iel)%v
!    enddo
!    strain_bx = strain_bx/dum
!
!    do i = 1,boxes(ib)%nel
!      iel = boxes(ib)%element_ids(i)
!      dum = - (elements(iel)%strain(1) + elements(iel)%strain(2) + elements(iel)%strain(3))/3.0 + &
!      (strain_bx(1) + strain_bx(2) + strain_bx(3))/3.0
!      do j = 1,3
!        elements(iel)%strain(j) = elements(iel)%strain(j) + dum
!      enddo
!      ! write(*,*) iel, igr, sum(elements(iel)%strain(1:3))
!      strainavg = strainavg + elements(iel)%strain*elements(iel)%v
!    enddo
!
!  enddo
!  strainavg = strainavg/vtot

end


! subroutine non_local_hyd_correction(elements, faces)
!   use global
!   use types
!   use tensor_functions
!   type (element_type), allocatable, intent (inout) :: elements(:)
!   type (face_type), allocatable, intent (inout) :: faces(:)
! 
!   integer :: iel, itnonloc, i, ifc, ielngh
!   double precision :: dum, dum1
! 
!   do itnonloc = 1,nit_non_loc
!     write(*,*) 'non local regularization iter:', itnonloc
!     do iel = 1,nel
!       dum = 0.0
!       dum1 = 0.0
!       do i = 1,4
!         ifc = elements(iel)%face(i)
!         if (faces(ifc)%nel == 2) then
!           ielngh = faces(ifc)%el(1)
!           if (ielngh == iel) ielngh = faces(ifc)%el(2)
!           dum = dum + (elements(ielngh)%strain(1) + elements(ielngh)%strain(2) + elements(ielngh)%strain(3))/3.0*&
!            elements(ielngh)%v
!           dum1 = dum1 + elements(ielngh)%v
!         endif
!       enddo
!       dum = dum + (elements(iel)%strain(1) + elements(iel)%strain(2) + elements(iel)%strain(3))/3.0*&
!        elements(iel)%v
!       dum1 = dum1 + elements(iel)%v
! 
!       elements(iel)%strain_hyd_non_loc = dum/dum1
!     enddo
! 
!     strainavg = 0.0
!     do iel = 1,nel
!       dum = - (elements(iel)%strain(1) + elements(iel)%strain(2) + elements(iel)%strain(3))/3.0 + &
!        elements(iel)%strain_hyd_non_loc
!       do i = 1,3
!         elements(iel)%strain(i) = elements(iel)%strain(i) + dum
!         elements(iel)%disgrad(i,i) = elements(iel)%disgrad(i,i) + dum
!       enddo
!       strainavg = strainavg + elements(iel)%strain*elements(iel)%v
!       ! write(*,*)igr,sum(elements(iel)%strain(1:3))
!     enddo
!     strainavg = strainavg/vtot
!   enddo
!   ! stop
! 
! end

end module various_functions
