  pure function TijkTjk_1d(Tijk,Tjk) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijk(3,3,3), Tjk(3,3)
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
    double precision, intent(in) :: Tijkl(3,3,3,3), Tjkl(3,6)
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
    double precision :: Tikl(3,6), dum
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

  pure function TijklmTjklm_1d(Tijklm,Tjklm) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjklm(3,10)
    double precision :: Ti(3), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = k,3
      do m = l,3
        ind2 = ind3d21d(k,l,m)
        dum = dum + Tijklm(i,j,k,l,m)*Tjklm(j,ind2)
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmTjklm_1d

  pure function TijklmTjlm_1d(Tijklm,Tjlm) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjlm(3,6)
    double precision :: Tik(3,3), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = l,3
        ind2 = ind2d21d(l,m)
        dum = dum + Tijklm(i,j,k,l,m)*Tjlm(j,ind2)
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmTjlm_1d

  pure function TijklmTjm_1d(Tijklm,Tjm) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjm(3,3)
    double precision :: Tikl(3,6), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tjm(j,m)
      enddo
      enddo
      ind1 = ind2d21d(k,l)
      Tikl(i,ind1) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmTjm_1d

  pure function TijklmTj_1d(Tijklm,Tj) result (Tiklm)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tj(3)
    double precision :: Tiklm(3,10), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tj(j)
      enddo
      ind1 = ind3d21d(k,l,m)
      Tiklm(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmTj_1d

  pure function TijklmnTjklmn_1d(Tijklmn,Tjklmn) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjklmn(3,15)
    double precision :: Ti(3), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = k,3
      do m = l,3
      do n = m,3
        ind2 = ind4d21d(k,l,m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjklmn(j,ind2)
      enddo
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmnTjklmn_1d

  pure function TijklmnTjlmn_1d(Tijklmn,Tjlmn) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjlmn(3,10)
    double precision :: Tik(3,3), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = l,3
      do n = m,3
        ind2 = ind3d21d(l,m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjlmn(j,ind2)
      enddo
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmnTjlmn_1d

  pure function TijklmnTjmn_1d(Tijklmn,Tjmn) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjmn(3,6)
    double precision :: Tikl(3,6), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
      do n = m,3
        ind2 = ind2d21d(m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjmn(j,ind2)
      enddo
      enddo
      enddo
      ind1 = ind2d21d(k,l)
      Tikl(i,ind1) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmnTjmn_1d

  pure function TijklmnTjn_1d(Tijklmn,Tjn) result (Tiklm)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjn(3,3)
    double precision :: Tiklm(3,10), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
      dum = 0.0
      do j = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjn(j,n)
      enddo
      enddo
      ind1 = ind3d21d(k,l,m)
      Tiklm(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTjn_1d

  pure function TijklmnTj_1d(Tijklmn,Tj) result (Tiklmn)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tj(3)
    double precision :: Tiklmn(3,15), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
    do n = m,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tj(j)
      enddo
      ind1 = ind4d21d(k,l,m,n)
      Tiklmn(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTj_1d

  pure function TijkTjk_1d_cmplx(Tijk,Tjk) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijk(3,3,3), Tjk(3,3)
    double complex :: Ti(3), dum
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
  end function TijkTjk_1d_cmplx

  pure function TijkTj_1d_cmplx(Tijk,Tj) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijk(3,3,3), Tj(3)
    double complex :: Tik(3,3), dum
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
  end function TijkTj_1d_cmplx

  pure function TijklTjkl_1d_cmplx(Tijkl,Tjkl) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijkl(3,3,3,3), Tjkl(3,6)
    double complex :: Ti(3), dum
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
  end function TijklTjkl_1d_cmplx

  pure function TijklTjk_1d_cmplx(Tijkl,Tjk) result (Til)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijkl(3,3,3,3), Tjk(3,3)
    double complex :: Til(3,3), dum
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
  end function TijklTjk_1d_cmplx

  pure function TijklTj_1d_cmplx(Tijkl,Tj) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijkl(3,3,3,3), Tj(3)
    double complex :: Tikl(3,6), dum
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
  end function TijklTj_1d_cmplx

  pure function TijklmTjklm_1d_cmplx(Tijklm,Tjklm) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklm(3,3,3,3,3), Tjklm(3,10)
    double complex :: Ti(3), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = k,3
      do m = l,3
        ind2 = ind3d21d(k,l,m)
        dum = dum + Tijklm(i,j,k,l,m)*Tjklm(j,ind2)
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmTjklm_1d_cmplx

  pure function TijklmTjlm_1d_cmplx(Tijklm,Tjlm) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklm(3,3,3,3,3), Tjlm(3,6)
    double complex :: Tik(3,3), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = l,3
        ind2 = ind2d21d(l,m)
        dum = dum + Tijklm(i,j,k,l,m)*Tjlm(j,ind2)
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmTjlm_1d_cmplx

  pure function TijklmTjm_1d_cmplx(Tijklm,Tjm) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklm(3,3,3,3,3), Tjm(3,3)
    double complex :: Tikl(3,6), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tjm(j,m)
      enddo
      enddo
      ind1 = ind2d21d(k,l)
      Tikl(i,ind1) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmTjm_1d_cmplx

  pure function TijklmTj_1d_cmplx(Tijklm,Tj) result (Tiklm)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklm(3,3,3,3,3), Tj(3)
    double complex :: Tiklm(3,10), dum
    integer :: i,j,k,l,m
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tj(j)
      enddo
      ind1 = ind3d21d(k,l,m)
      Tiklm(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmTj_1d_cmplx

  pure function TijklmnTjklmn_1d_cmplx(Tijklmn,Tjklmn) result (Ti)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjklmn(3,15)
    double complex :: Ti(3), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = k,3
      do m = l,3
      do n = m,3
        ind2 = ind4d21d(k,l,m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjklmn(j,ind2)
      enddo
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmnTjklmn_1d_cmplx

  pure function TijklmnTjlmn_1d_cmplx(Tijklmn,Tjlmn) result (Tik)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjlmn(3,10)
    double complex :: Tik(3,3), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = l,3
      do n = m,3
        ind2 = ind3d21d(l,m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjlmn(j,ind2)
      enddo
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmnTjlmn_1d_cmplx

  pure function TijklmnTjmn_1d_cmplx(Tijklmn,Tjmn) result (Tikl)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjmn(3,6)
    double complex :: Tikl(3,6), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
      do n = m,3
        ind2 = ind2d21d(m,n)
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjmn(j,ind2)
      enddo
      enddo
      enddo
      ind1 = ind2d21d(k,l)
      Tikl(i,ind1) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmnTjmn_1d_cmplx

  pure function TijklmnTjn_1d_cmplx(Tijklmn,Tjn) result (Tiklm)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjn(3,3)
    double complex :: Tiklm(3,10), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
      dum = 0.0
      do j = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjn(j,n)
      enddo
      enddo
      ind1 = ind3d21d(k,l,m)
      Tiklm(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTjn_1d_cmplx

  pure function TijklmnTj_1d_cmplx(Tijklmn,Tj) result (Tiklmn)
    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
    implicit none
    double complex, intent(in) :: Tijklmn(3,3,3,3,3,3), Tj(3)
    double complex :: Tiklmn(3,15), dum
    integer :: i,j,k,l,m,n
    integer :: ind1,ind2

    do i = 1,3
    do k = 1,3
    do l = k,3
    do m = l,3
    do n = m,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tj(j)
      enddo
      ind1 = ind4d21d(k,l,m,n)
      Tiklmn(i,ind1) = dum
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTj_1d_cmplx

