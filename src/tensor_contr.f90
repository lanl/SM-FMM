  pure function TijklmTklm(Tijklm,Tklm) result (Tij)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tklm(3,3,3)
    double precision :: Tij(3,3), dum
    integer :: i,j,k,l,m

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tklm(k,l,m)
      enddo
      enddo
      enddo
      Tij(i,j) = dum
    enddo
    enddo

    return
  end function TijklmTklm

  pure function TmlkjiTmlk_cmplx(Tmlkji,Tmlk) result (Tji)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tmlk(3,3,3)
    double complex :: Tji(3,3), dum
    integer :: i,j,k,l,m

    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tmlk(m,l,k)
      enddo
      enddo
      enddo
      Tji(j,i) = dum
    enddo
    enddo

    return
  end function TmlkjiTmlk_cmplx

  pure function TijklmTkl(Tijklm,Tkl) result (Tijm)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tkl(3,3)
    double precision :: Tijm(3,3,3), dum
    integer :: i,j,k,l,m

    do i = 1,3
    do j = 1,3
    do m = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tkl(k,l)
      enddo
      enddo
      Tijm(i,j,m) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmTkl

  pure function TijkTjk(Tijk,Tjk) result (Ti)
    implicit none
    double precision, intent(in) :: Tijk(3,3,3), Tjk(3,3)
    double precision :: Ti(3), dum
    integer :: i,j,k

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
  end function TijkTjk

  pure function TijklTjkl(Tijkl,Tjkl) result (Ti)
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tjkl(3,3,3)
    double precision :: Ti(3), dum
    integer :: i,j,k,l

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tjkl(j,k,l)
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklTjkl

  pure function TijkTj(Tijk,Tj) result (Tik)
    implicit none
    double precision, intent(in) :: Tijk(3,3,3), Tj(3)
    double precision :: Tik(3,3), dum
    integer :: i,j,k

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
  end function TijkTj

  pure function TkjiTj_cmplx(Tkji,Tj) result (Tki)
    implicit none
    double complex, intent(in) :: Tkji(3,3,3), Tj(3)
    double complex :: Tki(3,3), dum
    integer :: i,j,k

    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tkji(k,j,i)*Tj(j)
      enddo
      Tki(k,i) = dum
    enddo
    enddo

    return
  end function TkjiTj_cmplx

  pure function TijklTjk(Tijkl,Tjk) result (Til)
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tjk(3,3)
    double precision :: Til(3,3), dum
    integer :: i,j,k,l

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
  end function TijklTjk

  pure function TlkjiTkj_cmplx(Tlkji,Tkj) result (Tli)
    implicit none
    double complex, intent(in) :: Tlkji(3,3,3,3), Tkj(3,3)
    double complex :: Tli(3,3), dum
    integer :: i,j,k,l

    do l = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
        dum = dum + Tlkji(l,k,j,i)*Tkj(k,j)
      enddo
      enddo
      Tli(l,i) = dum
    enddo
    enddo

    return
  end function TlkjiTkj_cmplx

  pure function TijklTj(Tijkl,Tj) result (Tikl)
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tj(3)
    double precision :: Tikl(3,3,3), dum
    integer :: i,j,k,l

    do i = 1,3
    do k = 1,3
    do l = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tj(j)
      enddo
      Tikl(i,k,l) = dum
    enddo
    enddo
    enddo

    return
  end function TijklTj

  pure function TlkjiTj_cmplx(Tlkji,Tj) result (Tlki)
    implicit none
    double complex, intent(in) :: Tlkji(3,3,3,3), Tj(3)
    double complex :: Tlki(3,3,3), dum
    integer :: i,j,k,l

    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tlkji(l,k,j,i)*Tj(j)
      enddo
      Tlki(l,k,i) = dum
    enddo
    enddo
    enddo

    return
  end function TlkjiTj_cmplx

  pure function TmlkjiTlk_cmplx(Tmlkji,Tlk) result (Tmji)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tlk(3,3)
    double complex :: Tmji(3,3,3), dum
    integer :: i,j,k,l,m

    do m = 1,3
    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tlk(l,k)
      enddo
      enddo
      Tmji(m,j,i) = dum
    enddo
    enddo
    enddo

    return
  end function TmlkjiTlk_cmplx

  pure function TijkTk(Tijk,Tk) result (Tij)
    implicit none
    double precision, intent(in) :: Tijk(3,3,3), Tk(3)
    double precision :: Tij(3,3), dum
    integer :: i,j,k

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
        dum = dum + Tijk(i,j,k)*Tk(k)
      enddo
      Tij(i,j) = dum
    enddo
    enddo

    return
  end function TijkTk

  pure function TijklmnTkl(Tijklmn,Tkl) result (Tijmn)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tkl(3,3)
    double precision :: Tijmn(3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do j = 1,3
    do m = 1,3
    do n = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tkl(k,l)
      enddo
      enddo
      Tijmn(i,j,m,n) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTkl

  pure function TnmlkjiTlk_cmplx(Tnmlkji,Tlk) result (Tnmji)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tlk(3,3)
    double complex :: Tnmji(3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do n = 1,3
    do m = 1,3
    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tlk(l,k)
      enddo
      enddo
      Tnmji(n,m,j,i) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTlk_cmplx

  pure function TijklmnTklm(Tijklmn,Tklm) result (Tijn)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tklm(3,3,3)
    double precision :: Tijn(3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do j = 1,3
    do n = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tklm(k,l,m)
      enddo
      enddo
      enddo
      Tijn(i,j,n) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmnTklm

  pure function TnmlkjiTmlk_cmplx(Tnmlkji,Tmlk) result (Tnji)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tmlk(3,3,3)
    double complex :: Tnji(3,3,3), dum
    integer :: i,j,k,l,m,n

    do n = 1,3
    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tmlk(m,l,k)
      enddo
      enddo
      enddo
      Tnji(n,j,i) = dum
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTmlk_cmplx

  pure function TijklmnTkln(Tijklmn,Tkln) result (Tijm)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tkln(3,3,3)
    double precision :: Tijm(3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do j = 1,3
    do m = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tkln(k,l,n)
      enddo
      enddo
      enddo
      Tijm(i,j,m) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmnTkln

  pure function TnmlkjiTnlk_cmplx(Tnmlkji,Tnlk) result (Tmji)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnlk(3,3,3)
    double complex :: Tmji(3,3,3), dum
    integer :: i,j,k,l,m,n

    do m = 1,3
    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnlk(n,l,k)
      enddo
      enddo
      enddo
      Tmji(m,j,i) = dum
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTnlk_cmplx

  pure function TijklmnTklmn(Tijklmn,Tklmn) result (Tij)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tklmn(3,3,3,3)
    double precision :: Tij(3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tklmn(k,l,m,n)
      enddo
      enddo
      enddo
      enddo
      Tij(i,j) = dum
    enddo
    enddo

    return
  end function TijklmnTklmn

  pure function TnmlkjiTnmlk_cmplx(Tnmlkji,Tnmlk) result (Tji)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnmlk(3,3,3,3)
    double complex :: Tji(3,3), dum
    integer :: i,j,k,l,m,n

    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnmlk(n,m,l,k)
      enddo
      enddo
      enddo
      enddo
      Tji(j,i) = dum
    enddo
    enddo

    return
  end function TnmlkjiTnmlk_cmplx

  pure function TijklTkl(Tijkl,Tkl) result (Tij)
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tkl(3,3)
    double precision :: Tij(3,3), dum
    integer :: i,j,k,l

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tkl(k,l)
      enddo
      enddo
      Tij(i,j) = dum
    enddo
    enddo

    return
  end function TijklTkl

  pure function TlkjiTlk_cmplx(Tlkji,Tlk) result (Tji)
    implicit none
    double complex, intent(in) :: Tlkji(3,3,3,3), Tlk(3,3)
    double complex :: Tji(3,3), dum
    integer :: i,j,k,l

    do j = 1,3
    do i = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tlkji(l,k,j,i)*Tlk(l,k)
      enddo
      enddo
      Tji(j,i) = dum
    enddo
    enddo

    return
  end function TlkjiTlk_cmplx

  pure function TijklTklmn(Tijkl,Tklmn) result (Tijmn)
    implicit none
    double precision, intent(in) :: Tijkl(3,3,3,3), Tklmn(3,3,3,3)
    double precision :: Tijmn(3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do j = 1,3
    do m = 1,3
    do n = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + Tijkl(i,j,k,l)*Tklmn(k,l,m,n)
      enddo
      enddo
      Tijmn(i,j,m,n) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklTklmn

  pure function TijklmTjklm(Tijklm,Tjklm) result (Ti)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjklm(3,3,3,3)
    double precision :: Ti(3), dum
    integer :: i,j,k,l,m

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tjklm(j,k,l,m)
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmTjklm

  pure function TijklmnTjklmn(Tijklmn,Tjklmn) result (Ti)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjklmn(3,3,3,3,3)
    double precision :: Ti(3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjklmn(j,k,l,m,n)
      enddo
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijklmnTjklmn

  pure function TijklmTjlm(Tijklm,Tjlm) result (Tik)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjlm(3,3,3)
    double precision :: Tik(3,3), dum
    integer :: i,j,k,l,m

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tjlm(j,l,m)
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmTjlm

  pure function TmlkjiTmlj_cmplx(Tmlkji,Tmlj) result (Tki)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tmlj(3,3,3)
    double complex :: Tki(3,3), dum
    integer :: i,j,k,l,m

    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tmlj(m,l,j)
      enddo
      enddo
      enddo
      Tki(k,i) = dum
    enddo
    enddo

    return
  end function TmlkjiTmlj_cmplx

  pure function TijklmTjm(Tijklm,Tjm) result (Tikl)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tjm(3,3)
    double precision :: Tikl(3,3,3), dum
    integer :: i,j,k,l,m

    do i = 1,3
    do k = 1,3
    do l = 1,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tjm(j,m)
      enddo
      enddo
      Tikl(i,k,l) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmTjm

  pure function TmlkjiTmj_cmplx(Tmlkji,Tmj) result (Tlki)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tmj(3,3)
    double complex :: Tlki(3,3,3), dum
    integer :: i,j,k,l,m

    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tmj(m,j)
      enddo
      enddo
      Tlki(l,k,i) = dum
    enddo
    enddo
    enddo

    return
  end function TmlkjiTmj_cmplx

  pure function TijklmTj(Tijklm,Tj) result (Tiklm)
    implicit none
    double precision, intent(in) :: Tijklm(3,3,3,3,3), Tj(3)
    double precision :: Tiklm(3,3,3,3), dum
    integer :: i,j,k,l,m

    do i = 1,3
    do k = 1,3
    do l = 1,3
    do m = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklm(i,j,k,l,m)*Tj(j)
      enddo
      Tiklm(i,k,l,m) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmTj

  pure function TmlkjiTj_cmplx(Tmlkji,Tj) result (Tmlki)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tj(3)
    double complex :: Tmlki(3,3,3,3), dum
    integer :: i,j,k,l,m

    do m = 1,3
    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tj(j)
      enddo
      Tmlki(m,l,k,i) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TmlkjiTj_cmplx

  pure function TijklmnTjlmn(Tijklmn,Tjlmn) result (Tik)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjlmn(3,3,3,3)
    double precision :: Tik(3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do k = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjlmn(j,l,m,n)
      enddo
      enddo
      enddo
      enddo
      Tik(i,k) = dum
    enddo
    enddo

    return
  end function TijklmnTjlmn

  pure function TijklmnTjmn(Tijklmn,Tjmn) result (Tikl)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjmn(3,3,3)
    double precision :: Tikl(3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do k = 1,3
    do l = 1,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjmn(j,m,n)
      enddo
      enddo
      enddo
      Tikl(i,k,l) = dum
    enddo
    enddo
    enddo

    return
  end function TijklmnTjmn

  pure function TijklmnTjn(Tijklmn,Tjn) result (Tiklm)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tjn(3,3)
    double precision :: Tiklm(3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do k = 1,3
    do l = 1,3
    do m = 1,3
      dum = 0.0
      do j = 1,3
      do n = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tjn(j,n)
      enddo
      enddo
      Tiklm(i,k,l,m) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTjn

  pure function TijklmnTj(Tijklmn,Tj) result (Tiklmn)
    implicit none
    double precision, intent(in) :: Tijklmn(3,3,3,3,3,3), Tj(3)
    double precision :: Tiklmn(3,3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
    do k = 1,3
    do l = 1,3
    do m = 1,3
    do n = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tijklmn(i,j,k,l,m,n)*Tj(j)
      enddo
      Tiklmn(i,k,l,m,n) = dum
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end function TijklmnTj

  pure function TkjiTkj_cmplx(Tkji,Tkj) result (Ti)
    implicit none
    double complex, intent(in) :: Tkji(3,3,3), Tkj(3,3)
    double complex :: Ti(3), dum
    integer :: i,j,k

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
        dum = dum + Tkji(k,j,i)*Tkj(k,j)
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TkjiTkj_cmplx

  pure function TlkjiTlkj_cmplx(Tlkji,Tlkj) result (Ti)
    implicit none
    double complex, intent(in) :: Tlkji(3,3,3,3), Tlkj(3,3,3)
    double complex :: Ti(3), dum
    integer :: i,j,k,l

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
        dum = dum + Tlkji(l,k,j,i)*Tlkj(l,k,j)
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TlkjiTlkj_cmplx

  pure function TmlkjiTmlkj_cmplx(Tmlkji,Tmlkj) result (Ti)
    implicit none
    double complex, intent(in) :: Tmlkji(3,3,3,3,3), Tmlkj(3,3,3,3)
    double complex :: Ti(3), dum
    integer :: i,j,k,l,m

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + Tmlkji(m,l,k,j,i)*Tmlkj(m,l,k,j)
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TmlkjiTmlkj_cmplx

  pure function TnmlkjiTnmlkj_cmplx(Tnmlkji,Tnmlkj) result (Ti)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnmlkj(3,3,3,3,3)
    double complex :: Ti(3), dum
    integer :: i,j,k,l,m,n

    do i = 1,3
      dum = 0.0
      do j = 1,3
      do k = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnmlkj(n,m,l,k,j)
      enddo
      enddo
      enddo
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TnmlkjiTnmlkj_cmplx

  pure function TnmlkjiTnmlj_cmplx(Tnmlkji,Tnmlj) result (Tki)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnmlj(3,3,3,3)
    double complex :: Tki(3,3), dum
    integer :: i,j,k,l,m,n

    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do l = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnmlj(n,m,l,j)
      enddo
      enddo
      enddo
      enddo
      Tki(k,i) = dum
    enddo
    enddo

    return
  end function TnmlkjiTnmlj_cmplx

  pure function TnmlkjiTnmj_cmplx(Tnmlkji,Tnmj) result (Tlki)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnmj(3,3,3)
    double complex :: Tlki(3,3,3), dum
    integer :: i,j,k,l,m,n

    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do m = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnmj(n,m,j)
      enddo
      enddo
      enddo
      Tlki(l,k,i) = dum
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTnmj_cmplx

  pure function TnmlkjiTnj_cmplx(Tnmlkji,Tnj) result (Tmlki)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tnj(3,3)
    double complex :: Tmlki(3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do m = 1,3
    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
      do n = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tnj(n,j)
      enddo
      enddo
      Tmlki(m,l,k,i) = dum
    enddo
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTnj_cmplx

  pure function TnmlkjiTj_cmplx(Tnmlkji,Tj) result (Tnmlki)
    implicit none
    double complex, intent(in) :: Tnmlkji(3,3,3,3,3,3), Tj(3)
    double complex :: Tnmlki(3,3,3,3,3), dum
    integer :: i,j,k,l,m,n

    do n = 1,3
    do m = 1,3
    do l = 1,3
    do k = 1,3
    do i = 1,3
      dum = 0.0
      do j = 1,3
        dum = dum + Tnmlkji(n,m,l,k,j,i)*Tj(j)
      enddo
      Tnmlki(n,m,l,k,i) = dum
    enddo
    enddo
    enddo
    enddo
    enddo

    return
  end function TnmlkjiTj_cmplx

