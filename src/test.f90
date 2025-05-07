  pure function TijnTjn_1d(Tijn,Tjn) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijn(3,3,3), Tjn(3,6)
    double precision :: Ti(3), dum
    integer :: i,j,k
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
        dum = dum + Tijk(i,j,k)*Tjk(j,k)
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijkTjk_1d

  pure function TijkTj_1d(Tijk,Tj) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijk(3,3,3), Tj(3)
    double precision :: Tik(3,3), dum
    integer :: i,j,k
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijk(i,j,k)*Tj(j)
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijkTj_1d

  pure function TijklTjkl_1d(Tijkl,Tjkl) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tjkl(3,10)
    double precision :: Ti(3), dum
    integer :: i,j,k,l
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = k,3
        ind2 = ind2d21d(k,l)
        dum = dum + Tijkl(i,j,k,l)*Tjkl(j,ind2)
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklTjkl_1d

  pure function TijklTjk_1d(Tijkl,Tjk) result (Til)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tjk(3,3)
    double precision :: Til(3,3), dum
    integer :: i,j,k,l
    integer :: ind1,ind2

    do i = 1,3
    do l = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tjk(j,k)
      enddo
      enddo
      Til(i,l) = dum
    enddo
    enddo

    return
  end function TijklTjk_1d

  pure function TijklTj_1d(Tijkl,Tj) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tj(3)
    double precision :: Tikl(3,10), dum
    integer :: i,j,k,l
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tj(j)
      enddo
      ind1 = ind2d21d(k,l)
      Tikl(i,ind1) = dum
    enddo
    enddo
    enddo

    return
  end function TijklTj_1d

