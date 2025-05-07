module tensor_functions
  implicit none

contains
  ! ************************************************************************
  !     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
  !
  !     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
  !     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
  !     (modif. 10/MAY/01 - KDIM version - R.L.)
  !
  !     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
  !     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
  !     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
  !             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
  !     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
  !     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
  !             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
  !     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
  ! **************************************************************************
  SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)
    implicit none

    integer, intent(in) :: iopt, kdim
    double precision, parameter :: RSQ2=sqrt(2.)/2.
    double precision, parameter :: NRSQ2=-sqrt(2.)/2.
    double precision, parameter :: RSQ3=sqrt(3.)/3.
    double precision, parameter :: RSQ6=sqrt(6.)/6.
    double precision, parameter :: NRSQ6_2=-sqrt(6.)/3.
    double precision, parameter :: RSQ23=sqrt(2.)*sqrt(3.)/3.
    double precision :: CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)

    ! *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
    IF(IOPT.EQ.1) THEN
      IF(KDIM==5) THEN
        C2(1,1)=-RSQ2*CE2(1)-RSQ6*CE2(2); C2(1,2)=RSQ2*CE2(5);             C2(1,3)=RSQ2*CE2(4)
        C2(2,1)=RSQ2*CE2(5);              C2(2,2)=RSQ2*CE2(1)-RSQ6*CE2(2); C2(2,3)=RSQ2*CE2(3)
        C2(3,1)=RSQ2*CE2(4);              C2(3,2)=RSQ2*CE2(3);             C2(3,3)=RSQ23*CE2(2)
      ELSE IF(KDIM==6) THEN
        C2(1,1)=RSQ3*CE2(6)-RSQ2*CE2(1)-RSQ6*CE2(2); C2(1,2)=RSQ2*CE2(5);                         C2(1,3)=RSQ2*CE2(4)
        C2(2,1)=RSQ2*CE2(5);                         C2(2,2)=RSQ3*CE2(6)+RSQ2*CE2(1)-RSQ6*CE2(2); C2(2,3)=RSQ2*CE2(3)
        C2(3,1)=RSQ2*CE2(4);                         C2(3,2)=RSQ2*CE2(3);                         C2(3,3)=RSQ3*CE2(6)+RSQ23*CE2(2)
      ENDIF
    ENDIF

    ! *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
    IF(IOPT.EQ.2) THEN
      CE2(1) =  RSQ2*C2(2,2) - RSQ2*C2(1,1)
      CE2(2) = RSQ23*C2(3,3) - RSQ6*C2(2,2) - RSQ6*C2(1,1)
      CE2(3) =  RSQ2*C2(2,3) + RSQ2*C2(3,2)
      CE2(4) =  RSQ2*C2(1,3) + RSQ2*C2(3,1)
      CE2(5) =  RSQ2*C2(1,2) + RSQ2*C2(2,1)
      if (KDIM == 6) then
        CE2(6) = RSQ3*C2(1,1) + RSQ3*C2(2,2) + RSQ3*C2(3,3)
      endif
    ENDIF

    ! *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
    IF(IOPT.EQ.3) THEN
      C4(1,1,1,1)=CE4(1,1)/2. + CE4(2,2)/6. + CE4(6,6)/3. + RSQ3/2*CE4(1,2) - RSQ6*CE4(1,6) + RSQ3/2*CE4(2,1) &
                - RSQ2/3*CE4(2,6) - RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(1,1,1,2)=RSQ6*CE4(6,5) - RSQ3/2*CE4(2,5) - CE4(1,5)/2.
      C4(1,1,1,3)=RSQ6*CE4(6,4) - RSQ3/2*CE4(2,4) - CE4(1,4)/2.
      C4(1,1,2,1)=RSQ6*CE4(6,5) - RSQ3/2*CE4(2,5) - CE4(1,5)/2.
      C4(1,1,2,2)=CE4(2,2)/6. - CE4(1,1)/2. + CE4(6,6)/3. + RSQ3/2*CE4(1,2) - RSQ6*CE4(1,6) - RSQ3/2*CE4(2,1) &
                - RSQ2/3*CE4(2,6) + RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(1,1,2,3)=RSQ6*CE4(6,3) - RSQ3/2*CE4(2,3) - CE4(1,3)/2.
      C4(1,1,3,1)=RSQ6*CE4(6,4) - RSQ3/2*CE4(2,4) - CE4(1,4)/2.
      C4(1,1,3,2)=RSQ6*CE4(6,3) - RSQ3/2*CE4(2,3) - CE4(1,3)/2.
      C4(1,1,3,3)=CE4(6,6)/3. - 2*CE4(2,2)/6. - RSQ3*CE4(1,2) - RSQ6*CE4(1,6) - RSQ2/3*CE4(2,6) + 2*RSQ2/3*CE4(6,2)
      C4(1,2,1,1)=RSQ6*CE4(5,6) - RSQ3/2*CE4(5,2) - CE4(5,1)/2.
      C4(1,2,1,2)=CE4(5,5)/2.
      C4(1,2,1,3)=CE4(5,4)/2.
      C4(1,2,2,1)=CE4(5,5)/2.
      C4(1,2,2,2)=CE4(5,1)/2. - RSQ3/2*CE4(5,2) + RSQ6*CE4(5,6)
      C4(1,2,2,3)=CE4(5,3)/2.
      C4(1,2,3,1)=CE4(5,4)/2.
      C4(1,2,3,2)=CE4(5,3)/2.
      C4(1,2,3,3)=RSQ3*CE4(5,2) + RSQ6*CE4(5,6)
      C4(1,3,1,1)=RSQ6*CE4(4,6) - RSQ3/2*CE4(4,2) - CE4(4,1)/2.
      C4(1,3,1,2)=CE4(4,5)/2.
      C4(1,3,1,3)=CE4(4,4)/2.
      C4(1,3,2,1)=CE4(4,5)/2.
      C4(1,3,2,2)=CE4(4,1)/2. - RSQ3/2*CE4(4,2) + RSQ6*CE4(4,6)
      C4(1,3,2,3)=CE4(4,3)/2.
      C4(1,3,3,1)=CE4(4,4)/2.
      C4(1,3,3,2)=CE4(4,3)/2.
      C4(1,3,3,3)=RSQ3*CE4(4,2) + RSQ6*CE4(4,6)
      C4(2,1,1,1)=RSQ6*CE4(5,6) - RSQ3/2*CE4(5,2) - CE4(5,1)/2.
      C4(2,1,1,2)=CE4(5,5)/2.
      C4(2,1,1,3)=CE4(5,4)/2.
      C4(2,1,2,1)=CE4(5,5)/2.
      C4(2,1,2,2)=CE4(5,1)/2. - RSQ3/2*CE4(5,2) + RSQ6*CE4(5,6)
      C4(2,1,2,3)=CE4(5,3)/2.
      C4(2,1,3,1)=CE4(5,4)/2.
      C4(2,1,3,2)=CE4(5,3)/2.
      C4(2,1,3,3)=RSQ3*CE4(5,2) + RSQ6*CE4(5,6)
      C4(2,2,1,1)=CE4(2,2)/6. - CE4(1,1)/2. + CE4(6,6)/3. - RSQ3/2*CE4(1,2) + RSQ6*CE4(1,6) + RSQ3/2*CE4(2,1) &
                - RSQ2/3*CE4(2,6) - RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(2,2,1,2)=CE4(1,5)/2. - RSQ3/2*CE4(2,5) + RSQ6*CE4(6,5)
      C4(2,2,1,3)=CE4(1,4)/2. - RSQ3/2*CE4(2,4) + RSQ6*CE4(6,4)
      C4(2,2,2,1)=CE4(1,5)/2. - RSQ3/2*CE4(2,5) + RSQ6*CE4(6,5)
      C4(2,2,2,2)=CE4(1,1)/2. + CE4(2,2)/6. + CE4(6,6)/3. - RSQ3/2*CE4(1,2) + RSQ6*CE4(1,6) - RSQ3/2*CE4(2,1) & 
                - RSQ2/3*CE4(2,6) + RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(2,2,2,3)=CE4(1,3)/2. - RSQ3/2*CE4(2,3) + RSQ6*CE4(6,3)
      C4(2,2,3,1)=CE4(1,4)/2. - RSQ3/2*CE4(2,4) + RSQ6*CE4(6,4)
      C4(2,2,3,2)=CE4(1,3)/2. - RSQ3/2*CE4(2,3) + RSQ6*CE4(6,3)
      C4(2,2,3,3)=CE4(6,6)/3. - 2*CE4(2,2)/6. + RSQ3*CE4(1,2) + RSQ6*CE4(1,6) - RSQ2/3*CE4(2,6) + 2*RSQ2/3*CE4(6,2)
      C4(2,3,1,1)=RSQ6*CE4(3,6) - RSQ3/2*CE4(3,2) - CE4(3,1)/2.
      C4(2,3,1,2)=CE4(3,5)/2.
      C4(2,3,1,3)=CE4(3,4)/2.
      C4(2,3,2,1)=CE4(3,5)/2.
      C4(2,3,2,2)=CE4(3,1)/2. - RSQ3/2*CE4(3,2) + RSQ6*CE4(3,6)
      C4(2,3,2,3)=CE4(3,3)/2.
      C4(2,3,3,1)=CE4(3,4)/2.
      C4(2,3,3,2)=CE4(3,3)/2.
      C4(2,3,3,3)=RSQ3*CE4(3,2) + RSQ6*CE4(3,6)
      C4(3,1,1,1)=RSQ6*CE4(4,6) - RSQ3/2*CE4(4,2) - CE4(4,1)/2.
      C4(3,1,1,2)=CE4(4,5)/2.
      C4(3,1,1,3)=CE4(4,4)/2.
      C4(3,1,2,1)=CE4(4,5)/2.
      C4(3,1,2,2)=CE4(4,1)/2. - RSQ3/2*CE4(4,2) + RSQ6*CE4(4,6)
      C4(3,1,2,3)=CE4(4,3)/2.
      C4(3,1,3,1)=CE4(4,4)/2.
      C4(3,1,3,2)=CE4(4,3)/2.
      C4(3,1,3,3)=RSQ3*CE4(4,2) + RSQ6*CE4(4,6)
      C4(3,2,1,1)=RSQ6*CE4(3,6) - RSQ3/2*CE4(3,2) - CE4(3,1)/2.
      C4(3,2,1,2)=CE4(3,5)/2.
      C4(3,2,1,3)=CE4(3,4)/2.
      C4(3,2,2,1)=CE4(3,5)/2.
      C4(3,2,2,2)=CE4(3,1)/2. - RSQ3/2*CE4(3,2) + RSQ6*CE4(3,6)
      C4(3,2,2,3)=CE4(3,3)/2.
      C4(3,2,3,1)=CE4(3,4)/2.
      C4(3,2,3,2)=CE4(3,3)/2.
      C4(3,2,3,3)=RSQ3*CE4(3,2) + RSQ6*CE4(3,6)
      C4(3,3,1,1)=CE4(6,6)/3. - 2*CE4(2,2)/6. - RSQ3*CE4(2,1) + 2*RSQ2/3*CE4(2,6) - RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(3,3,1,2)=RSQ3*CE4(2,5) + RSQ6*CE4(6,5)
      C4(3,3,1,3)=RSQ3*CE4(2,4) + RSQ6*CE4(6,4)
      C4(3,3,2,1)=RSQ3*CE4(2,5) + RSQ6*CE4(6,5)
      C4(3,3,2,2)=CE4(6,6)/3. - 2*CE4(2,2)/6. + RSQ3*CE4(2,1) + 2*RSQ2/3*CE4(2,6) + RSQ6*CE4(6,1) - RSQ2/3*CE4(6,2)
      C4(3,3,2,3)=RSQ3*CE4(2,3) + RSQ6*CE4(6,3)
      C4(3,3,3,1)=RSQ3*CE4(2,4) + RSQ6*CE4(6,4)
      C4(3,3,3,2)=RSQ3*CE4(2,3) + RSQ6*CE4(6,3)
      C4(3,3,3,3)=4*CE4(2,2)/6. + CE4(6,6)/3. + 2*RSQ2/3*CE4(2,6) + 2*RSQ2/3*CE4(6,2)
    ENDIF

    ! *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
    IF(IOPT.EQ.4) THEN
      CE4(1,1)=C4(1,1,1,1)/2. - C4(1,1,2,2)/2. - C4(2,2,1,1)/2. + C4(2,2,2,2)/2.
      CE4(1,2)=RSQ3/2*C4(1,1,1,1) + RSQ3/2*C4(1,1,2,2) - RSQ3*C4(1,1,3,3) - RSQ3/2*C4(2,2,1,1) - RSQ3/2*C4(2,2,2,2) &
            + RSQ3*C4(2,2,3,3)
      CE4(1,3)=C4(2,2,2,3)/2. - C4(1,1,3,2)/2. - C4(1,1,2,3)/2. + C4(2,2,3,2)/2.
      CE4(1,4)=C4(2,2,1,3)/2. - C4(1,1,3,1)/2. - C4(1,1,1,3)/2. + C4(2,2,3,1)/2.
      CE4(1,5)=C4(2,2,1,2)/2. - C4(1,1,2,1)/2. - C4(1,1,1,2)/2. + C4(2,2,2,1)/2.
      CE4(1,6)=RSQ6*C4(2,2,1,1) - RSQ6*C4(1,1,2,2) - RSQ6*C4(1,1,3,3) - RSQ6*C4(1,1,1,1) + RSQ6*C4(2,2,2,2) + RSQ6*C4(2,2,3,3)
      CE4(2,1)=RSQ3/2*C4(1,1,1,1) - RSQ3/2*C4(1,1,2,2) + RSQ3/2*C4(2,2,1,1) - RSQ3/2*C4(2,2,2,2) - RSQ3*C4(3,3,1,1) &
            + RSQ3*C4(3,3,2,2)
      CE4(2,2)=C4(1,1,1,1)/6. + C4(1,1,2,2)/6. - 2*C4(1,1,3,3)/6. + C4(2,2,1,1)/6. + C4(2,2,2,2)/6. - 2*C4(2,2,3,3)/6. &
            - 2*C4(3,3,1,1)/6. - 2*C4(3,3,2,2)/6. + 4*C4(3,3,3,3)/6.
      CE4(2,3)=RSQ3*C4(3,3,2,3) - RSQ3/2*C4(1,1,3,2) - RSQ3/2*C4(2,2,2,3) - RSQ3/2*C4(2,2,3,2) - RSQ3/2*C4(1,1,2,3) & 
            + RSQ3*C4(3,3,3,2)
      CE4(2,4)=RSQ3*C4(3,3,1,3) - RSQ3/2*C4(1,1,3,1) - RSQ3/2*C4(2,2,1,3) - RSQ3/2*C4(2,2,3,1) - RSQ3/2*C4(1,1,1,3) & 
            + RSQ3*C4(3,3,3,1)
      CE4(2,5)=RSQ3*C4(3,3,1,2) - RSQ3/2*C4(1,1,2,1) - RSQ3/2*C4(2,2,1,2) - RSQ3/2*C4(2,2,2,1) - RSQ3/2*C4(1,1,1,2) & 
            + RSQ3*C4(3,3,2,1)
      CE4(2,6)=2*RSQ2/3*C4(3,3,1,1) - RSQ2/3*C4(1,1,2,2) - RSQ2/3*C4(1,1,3,3) - RSQ2/3*C4(2,2,1,1) - RSQ2/3*C4(2,2,2,2) &
            - RSQ2/3*C4(2,2,3,3) - RSQ2/3*C4(1,1,1,1) + 2*RSQ2/3*C4(3,3,2,2) + 2*RSQ2/3*C4(3,3,3,3)
      CE4(3,1)=C4(2,3,2,2)/2. - C4(2,3,1,1)/2. - C4(3,2,1,1)/2. + C4(3,2,2,2)/2.
      CE4(3,2)=RSQ3*C4(2,3,3,3) - RSQ3/2*C4(2,3,2,2) - RSQ3/2*C4(2,3,1,1) - RSQ3/2*C4(3,2,1,1) - RSQ3/2*C4(3,2,2,2) &
            + RSQ3*C4(3,2,3,3)
      CE4(3,3)=C4(2,3,2,3)/2. + C4(2,3,3,2)/2. + C4(3,2,2,3)/2. + C4(3,2,3,2)/2.
      CE4(3,4)=C4(2,3,1,3)/2. + C4(2,3,3,1)/2. + C4(3,2,1,3)/2. + C4(3,2,3,1)/2.
      CE4(3,5)=C4(2,3,1,2)/2. + C4(2,3,2,1)/2. + C4(3,2,1,2)/2. + C4(3,2,2,1)/2.
      CE4(3,6)=RSQ6*C4(2,3,1,1) + RSQ6*C4(2,3,2,2) + RSQ6*C4(2,3,3,3) + RSQ6*C4(3,2,1,1) + RSQ6*C4(3,2,2,2) + RSQ6*C4(3,2,3,3)
      CE4(4,1)=C4(1,3,2,2)/2. - C4(1,3,1,1)/2. - C4(3,1,1,1)/2. + C4(3,1,2,2)/2.
      CE4(4,2)=RSQ3*C4(1,3,3,3) - RSQ3/2*C4(1,3,2,2) - RSQ3/2*C4(1,3,1,1) - RSQ3/2*C4(3,1,1,1) - RSQ3/2*C4(3,1,2,2) &
            + RSQ3*C4(3,1,3,3)
      CE4(4,3)=C4(1,3,2,3)/2. + C4(1,3,3,2)/2. + C4(3,1,2,3)/2. + C4(3,1,3,2)/2.
      CE4(4,4)=C4(1,3,1,3)/2. + C4(1,3,3,1)/2. + C4(3,1,1,3)/2. + C4(3,1,3,1)/2.
      CE4(4,5)=C4(1,3,1,2)/2. + C4(1,3,2,1)/2. + C4(3,1,1,2)/2. + C4(3,1,2,1)/2.
      CE4(4,6)=RSQ6*C4(1,3,1,1) + RSQ6*C4(1,3,2,2) + RSQ6*C4(1,3,3,3) + RSQ6*C4(3,1,1,1) + RSQ6*C4(3,1,2,2) + RSQ6*C4(3,1,3,3)
      CE4(5,1)=C4(1,2,2,2)/2. - C4(1,2,1,1)/2. - C4(2,1,1,1)/2. + C4(2,1,2,2)/2.
      CE4(5,2)=RSQ3*C4(1,2,3,3) - RSQ3/2*C4(1,2,2,2) - RSQ3/2*C4(1,2,1,1) - RSQ3/2*C4(2,1,1,1) - RSQ3/2*C4(2,1,2,2) &
            + RSQ3*C4(2,1,3,3)
      CE4(5,3)=C4(1,2,2,3)/2. + C4(1,2,3,2)/2. + C4(2,1,2,3)/2. + C4(2,1,3,2)/2.
      CE4(5,4)=C4(1,2,1,3)/2. + C4(1,2,3,1)/2. + C4(2,1,1,3)/2. + C4(2,1,3,1)/2.
      CE4(5,5)=C4(1,2,1,2)/2. + C4(1,2,2,1)/2. + C4(2,1,1,2)/2. + C4(2,1,2,1)/2.
      CE4(5,6)=RSQ6*C4(1,2,1,1) + RSQ6*C4(1,2,2,2) + RSQ6*C4(1,2,3,3) + RSQ6*C4(2,1,1,1) + RSQ6*C4(2,1,2,2) + RSQ6*C4(2,1,3,3)
      CE4(6,1)=RSQ6*C4(1,1,2,2) - RSQ6*C4(1,1,1,1) - RSQ6*C4(2,2,1,1) + RSQ6*C4(2,2,2,2) - RSQ6*C4(3,3,1,1) + RSQ6*C4(3,3,2,2)
      CE4(6,2)=2*RSQ2/3*C4(1,1,3,3) - RSQ2/3*C4(1,1,2,2) - RSQ2/3*C4(1,1,1,1) - RSQ2/3*C4(2,2,1,1) - RSQ2/3*C4(2,2,2,2) &
            + 2*RSQ2/3*C4(2,2,3,3) - RSQ2/3*C4(3,3,1,1) - RSQ2/3*C4(3,3,2,2) + 2*RSQ2/3*C4(3,3,3,3)
      CE4(6,3)=RSQ6*C4(1,1,2,3) + RSQ6*C4(1,1,3,2) + RSQ6*C4(2,2,2,3) + RSQ6*C4(2,2,3,2) + RSQ6*C4(3,3,2,3) + RSQ6*C4(3,3,3,2)
      CE4(6,4)=RSQ6*C4(1,1,1,3) + RSQ6*C4(1,1,3,1) + RSQ6*C4(2,2,1,3) + RSQ6*C4(2,2,3,1) + RSQ6*C4(3,3,1,3) + RSQ6*C4(3,3,3,1)
      CE4(6,5)=RSQ6*C4(1,1,1,2) + RSQ6*C4(1,1,2,1) + RSQ6*C4(2,2,1,2) + RSQ6*C4(2,2,2,1) + RSQ6*C4(3,3,1,2) + RSQ6*C4(3,3,2,1)
      CE4(6,6)=C4(1,1,1,1)/3. + C4(1,1,2,2)/3. + C4(1,1,3,3)/3. + C4(2,2,1,1)/3. + C4(2,2,2,2)/3. + C4(2,2,3,3)/3. &
            + C4(3,3,1,1)/3. + C4(3,3,2,2)/3. + C4(3,3,3,3)/3.
    ENDIF

    RETURN
  END SUBROUTINE CHG_BASIS

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !     SUBROUTINE VOIGT_VPSC   ---->   VERSION OF 09/02/98
  !
  !     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
  !     AND VICEVERSA IF IOPT=2.
  !     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
  !     AND VICEVERSA IF IOPT=4.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE VOIGT_VPSC(T1,T2,C2,C4,IOPT)
    implicit none

    integer, intent(in) :: iopt
    double precision :: T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
    integer :: IJV(6,2), i, j, n, m, i1, i2, j1, j2
    data ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

    IF(IOPT.EQ.1) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        T2(I1,I2)=T1(I)
        T2(I2,I1)=T1(I)
      ENDDO
    ENDIF

    IF(IOPT.EQ.2) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        T1(I)=T2(I1,I2)
      ENDDO
    ENDIF

    IF (IOPT.EQ.3) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        DO J=1,6
          J1=IJV(J,1)
          J2=IJV(J,2)
          C4(I1,I2,J1,J2)=C2(I,J)
          C4(I2,I1,J1,J2)=C2(I,J)
          C4(I1,I2,J2,J1)=C2(I,J)
          C4(I2,I1,J2,J1)=C2(I,J)
        ENDDO
      ENDDO
    ENDIF

    IF(IOPT.EQ.4) THEN
      DO I=1,6
        I1=IJV(I,1)
        I2=IJV(I,2)
        DO J=1,6
          J1=IJV(J,1)
          J2=IJV(J,2)
          C2(I,J)=C4(I1,I2,J1,J2)
        ENDDO
      ENDDO
    ENDIF

    RETURN
  END SUBROUTINE VOIGT_VPSC

  ! *************************************************************************
  ! CALCULATES THE VON MISES EQUIVALENT OF A NON-SYMMETRIC, NON-TRACELESS STRESS TENSOR
  ! *************************************************************************
  pure function vm_stress(dtensor)
    implicit none
    double precision, intent(in) :: dtensor(3,3)
    double precision :: h,vm_stress

    h = (dtensor(1,1)+dtensor(2,2)+dtensor(3,3))/3.
    vm_stress = sqrt(3./2. * ((dtensor(1,1) - h)**2 + (dtensor(2,2) - h)**2 + (dtensor(3,3) - h)**2 + & 
                               dtensor(1,2)**2      +  dtensor(1,3)**2      +  dtensor(2,3)**2     + &
                               dtensor(2,1)**2      +  dtensor(3,1)**2      +  dtensor(3,2)**2))
  end

  ! *************************************************************************
  ! CALCULATES THE VON MISES EQUIVALENT OF A NON-SYMMETRIC, NON-TRACELESS STRAIN TENSOR
  ! *************************************************************************
  pure function vm_strain(dtensor)
    implicit none
    double precision, intent(in) :: dtensor(3,3)
    double precision :: h,vm_strain

    h = (dtensor(1,1)+dtensor(2,2)+dtensor(3,3))/3.
    vm_strain = sqrt(2./3. * ((dtensor(1,1) - h)**2 + (dtensor(2,2) - h)**2 + (dtensor(3,3) - h)**2 + & 
                               dtensor(1,2)**2      +  dtensor(1,3)**2      +  dtensor(2,3)**2     + &
                               dtensor(2,1)**2      +  dtensor(3,1)**2      +  dtensor(3,2)**2))
  end

  ! *************************************************************************
  ! Calculates the principle stress and principle stress directions (eigenvalues & eigenvectors)
  ! Values are sorted from smallest to largest 
  ! *************************************************************************
  subroutine prin_stress(stensor,pstress,pstressdirs)
    implicit none
    double precision, intent(in) :: stensor(3,3)
    double precision, intent(out) :: pstress(3),pstressdirs(3,3)
    double precision :: work(64)
    integer :: info 

    pstressdirs=stensor
    call dsyev('V','U',3,pstressdirs,3,pstress,work,64,info) ! LAPACK Symmetric eigenvalues & vectors

    return
  end

  ! *************************************************************************
  ! Inverts a matrix using LU decomposition implemented in LAPACK
  ! a -> matrix to invert and inverted matrix
  ! n -> dimension of a
  ! *************************************************************************
  subroutine lu_inverse(a,n)
    implicit none
    integer :: n
    double precision :: a(n,n)
    double precision :: work(64)
    integer :: ipivot(n), info 

    call dgetrf(n,n,a,n,ipivot,info)
    if (info/=0) stop 'Singular matrix in lu_inverse'

    call dgetri(n,a,n,ipivot,work,64,info)
    if (info/=0) stop 'Failed inverse in lu_inverse'

    return 
  end subroutine lu_inverse

  ! *************************************************************************
  ! Solves a system of linear equations using LU decomposition implemented in LAPACK
  ! ax = b where b is replaced by x on exit. 
  ! n -> dimensions of a
  ! *************************************************************************
  subroutine lu_eqsystem(a,b,n)
    implicit none
    integer :: n
    double precision :: a(n,n),b(n)
    integer :: ipivot(n), info 

    call dgesv(n,1,a,n,ipivot,b,n,info) ! system of equations solver
    if (info/=0) stop 'Failed system solve in lu_eqsystem'

    return 
  end subroutine lu_eqsystem

  ! *************************************************************************
  ! Inverts a matrix using LU decomposition implemented in LAPACK
  ! F -> matrix to be decomposed
  ! R -> rotation matrix
  ! U -> right stretch tensor
  ! *************************************************************************
  subroutine polar_dcmp(F,R,U)
    implicit none

    double precision :: F(3,3),R(3,3),U(3,3),aux33(3,3)
    double precision :: Usvd(3,3),W3svd(3),VTsvd(3,3),Wsvd(3,3)
    double precision :: work201(201)
    integer :: info, i

    aux33=F

    call dgesvd('A','A',3,3,aux33,3,W3svd,Usvd,3,VTsvd,3,work201,201,info)
    if (info/=0) stop 'Error in LAPACK dgesvd'

    Wsvd=0.
    do i=1,3
      Wsvd(i,i)=W3svd(i)
    enddo

    R=matmul(Usvd,VTsvd)
    U=matmul(transpose(VTsvd),matmul(Wsvd,VTsvd))

  end subroutine polar_dcmp

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
  !
  !     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NROWSxNCOLS MATRICES
  !     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
  !     OF BOTH DATA.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  FUNCTION tmismatch (v1,v2,nrows,ncols)
    implicit none

    integer :: nrows, ncols
    double precision :: v1(36),v2(36)
    double precision :: v_dif(36),v_ave(36)
    double precision :: tmismatch
    integer i

    do i=1,nrows*ncols
      v_dif(i)=v1(i)-v2(i)
      v_ave(i)=0.5d0*(v1(i)+v2(i))
    enddo
    tmismatch=tnorm(v_dif,nrows,ncols)/tnorm(v_ave,nrows,ncols)

    RETURN
  END

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  !     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
  !
  !     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  FUNCTION tnorm(v,nrows,ncols)
    implicit none
    double precision :: v(36)
    integer :: nrows, ncols, i
    double precision :: tnorm

    tnorm=0.d0
    do i=1,nrows*ncols
      tnorm=tnorm+v(i)*v(i)
    enddo
    tnorm=sqrt(tnorm)

    RETURN
  END

  pure function determinant33(matrix)
    implicit none
    double precision, intent(in) :: matrix(3, 3)
    double precision :: determinant33

    determinant33 = &
       matrix(1,1)*matrix(2,2)*matrix(3,3) &
      +matrix(1,2)*matrix(2,3)*matrix(3,1) &
      +matrix(1,3)*matrix(2,1)*matrix(3,2) &
      -matrix(3,1)*matrix(2,2)*matrix(1,3) &
      -matrix(3,2)*matrix(2,3)*matrix(1,1) &
      -matrix(3,3)*matrix(2,1)*matrix(1,2)

    return
  end function determinant33

  pure function sym33(matrix)
    implicit none
    double precision, intent(in) :: matrix(3,3)
    double precision :: sym33(3,3)

    sym33(1,1)=matrix(1,1)
    sym33(2,2)=matrix(2,2)
    sym33(3,3)=matrix(3,3)
    sym33(1,2)=0.5*(matrix(1,2)+matrix(2,1))
    sym33(1,3)=0.5*(matrix(1,3)+matrix(3,1))
    sym33(2,3)=0.5*(matrix(2,3)+matrix(3,2))
    sym33(2,1)=sym33(1,2)
    sym33(3,1)=sym33(1,3)
    sym33(3,2)=sym33(2,3)

    return
  end function sym33

  pure function antisym33(matrix)
    implicit none
    double precision, intent(in) :: matrix(3,3)
    double precision :: antisym33(3,3)

    antisym33(1,1)=0.
    antisym33(2,2)=0.
    antisym33(3,3)=0.
    antisym33(1,2)=0.5*(matrix(1,2)-matrix(2,1))
    antisym33(1,3)=0.5*(matrix(1,3)-matrix(3,1))
    antisym33(2,3)=0.5*(matrix(2,3)-matrix(3,2))
    antisym33(2,1)=-antisym33(1,2)
    antisym33(3,1)=-antisym33(1,3)
    antisym33(3,2)=-antisym33(2,3)

    return
  end function antisym33

  pure function cross_prod(a,b) result (axb)
    implicit none
    double precision, intent(in) :: a(3), b(3)
    double precision :: axb(3)

    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

    return
  end function cross_prod

  pure function levi_civita(i,j,k)
    integer, intent(in) :: i,j,k
    integer :: levi_civita

    if ((i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)) then
      levi_civita = 1
    else if ((i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)) then
      levi_civita = -1
    else
      levi_civita = 0
    endif
    
  end function levi_civita

  subroutine euler(iopt,ph,th,tm,a)
    use global, only : pi
    implicit none
  
    integer :: iopt
    double precision :: ph, th, tm
    double precision :: a(3,3)
    double precision :: cph, cth, ctm, sph, sth, stm
  !     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
  !     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
  !     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
  !     ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.
    if(iopt.eq.1) then
      th=acos(a(3,3))
      if(abs(a(3,3)).ge.0.9999) then
        tm=0.
        ph=atan2(a(1,2),a(1,1))
      else
        sth=sin(th)
        tm=atan2(a(1,3)/sth,a(2,3)/sth)
        ph=atan2(a(3,1)/sth,-a(3,2)/sth)
      endif
      th=th*180./pi
      ph=ph*180./pi
      tm=tm*180./pi
    else if(iopt.eq.2) then
      sph=sin(ph)
      cph=cos(ph)
      sth=sin(th)
      cth=cos(th)
      stm=sin(tm)
      ctm=cos(tm)
      a(1,1)=ctm*cph-sph*stm*cth
      a(2,1)=-stm*cph-sph*ctm*cth
      a(3,1)=sph*sth
      a(1,2)=ctm*sph+cph*stm*cth
      a(2,2)=-sph*stm+cph*ctm*cth
      a(3,2)=-sth*cph
      a(1,3)=sth*stm
      a(2,3)=ctm*sth
      a(3,3)=cth
    endif
  
    return
  end

  pure function T4T2_cont(T4, T2)
    implicit none
    double precision, intent(in) :: T4(3,3,3,3), T2(3,3)
    double precision :: T4T2_cont(3,3), dum
    integer :: i, j, k, l

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
        dum = dum + T4(i,j,k,l)*T2(k,l)
      enddo
      enddo
      T4T2_cont(i,j) = dum
    enddo
    enddo

    return
  end function T4T2_cont

  pure function T5T3_cont(T5, T3)
    implicit none
    double precision, intent(in) :: T5(3,3,3,3,3), T3(3,3,3)
    double precision :: T5T3_cont(3,3), dum
    integer :: i, j, k, l, m

    do i = 1,3
    do j = 1,3
      dum = 0.0
      do k = 1,3
      do l = 1,3
      do m = 1,3
        dum = dum + T5(i,j,k,l,m)*T3(k,l,m)
      enddo
      enddo
      enddo
      T5T3_cont(i,j) = dum
    enddo
    enddo

    return
  end function T5T3_cont

  pure function enforce_periodic(i, minval, maxval) result (iper)
    integer, intent(in) :: i, minval, maxval
    integer :: iper

    if (i < minval) then
      iper = i + maxval
      ! if (iper < minval) iper = iper + maxval
    elseif (i > maxval) then
      iper = i - maxval
      ! if (iper > maxval) iper = iper - maxval
    else
      iper = i
    endif
  
  end

  pure function TijkTjk_var_dim(Tijk,Tjk) result (Ti)
    implicit none
    double precision, intent(in) :: Tijk(:,:,:), Tjk(:,:)
    double precision :: Ti(3), dum
    integer :: i,j,k

    do i = 1,3
      dum = 0.0
      do j = 1,size(Tijk,2)
      do k = 1,size(Tijk,3)
        dum = dum + Tijk(i,j,k)*Tjk(j,k)
      enddo
      enddo
      Ti(i) = dum
    enddo

    return
  end function TijkTjk_var_dim

!  pure function TijklTjkl_1d(Tijkl,Tjkl) result (Ti)
!    use global, only : ind2d21d, ind3d21d, ind4d21d, ind1d22d, ind1d23d, ind1d24d
!    implicit none
!    double precision, intent(in) :: Tijkl(3,3,3,3), Tjkl(3,6)
!    double precision :: Ti(3), dum
!    integer :: i,j,k,l,ind2
!
!    do i = 1,3
!      dum = 0.0
!      do j = 1,3
!      do k = 1,3
!      do l = k,3
!        ind2 = ind2d21d(k,l)
!        dum = dum + Tijkl(i,j,k,l)*Tjkl(j,ind2)
!      enddo
!      enddo
!      enddo
!      Ti(i) = dum
!    enddo
!
!    return
!  end function TijklTjkl_1d

  include 'tensor_contr.f90'
  include 'tensor_contr_1d.f90'

  pure function mandel_t2_to_t1(t2) result(t1)
    use global, only : sqrt2
    implicit none
    double precision, intent(in) :: t2(3,3)
    double precision :: t1(6)
    
    t1(1) = t2(1,1)
    t1(2) = t2(2,2)
    t1(3) = t2(3,3)
    t1(4) = t2(2,3)*sqrt2
    t1(5) = t2(1,3)*sqrt2
    t1(6) = t2(1,2)*sqrt2

  end

  pure function mandel_t1_to_t2(t1) result(t2)
    use global, only : sqrt2
    implicit none
    double precision, intent(in) :: t1(6)
    double precision :: t2(3,3)
    
    t2(1,1) = t1(1)
    t2(2,2) = t1(2)
    t2(3,3) = t1(3)
    t2(2,3) = t1(4)/sqrt2
    t2(1,3) = t1(5)/sqrt2
    t2(1,2) = t1(6)/sqrt2
    t2(3,2) = t2(2,3)
    t2(3,1) = t2(1,3)
    t2(2,1) = t2(1,2)

  end

  pure function mandel_t4_to_t2(t4) result(t2)
    use global, only : sqrt2
    implicit none
    double precision, intent(in) :: t4(3,3,3,3)
    double precision :: t2(6,6)
    
    t2(1,:) = [t4(1,1,1,1),t4(1,1,2,2),t4(1,1,3,3),t4(1,1,2,3)*sqrt2,t4(1,1,1,3)*sqrt2,t4(1,1,1,2)*sqrt2]
    t2(2,:) = [t4(2,2,1,1),t4(2,2,2,2),t4(2,2,3,3),t4(2,2,2,3)*sqrt2,t4(2,2,1,3)*sqrt2,t4(2,2,1,2)*sqrt2]
    t2(3,:) = [t4(3,3,1,1),t4(3,3,2,2),t4(3,3,3,3),t4(3,3,2,3)*sqrt2,t4(3,3,1,3)*sqrt2,t4(3,3,1,2)*sqrt2]
    t2(4,:) = [t4(2,3,1,1)*sqrt2,t4(2,3,2,2)*sqrt2,t4(2,3,3,3)*sqrt2,t4(2,3,2,3)*2.0,t4(2,3,1,3)*2.0,t4(2,3,1,2)*2.0]
    t2(5,:) = [t4(1,3,1,1)*sqrt2,t4(1,3,2,2)*sqrt2,t4(1,3,3,3)*sqrt2,t4(1,3,2,3)*2.0,t4(1,3,1,3)*2.0,t4(1,3,1,2)*2.0]
    t2(6,:) = [t4(1,2,1,1)*sqrt2,t4(1,2,2,2)*sqrt2,t4(1,2,3,3)*sqrt2,t4(1,2,2,3)*2.0,t4(1,2,1,3)*2.0,t4(1,2,1,2)*2.0]

  end

  pure function mandel_t2_to_t4(t2) result(t4)
    use global, only : sqrt2
    implicit none
    double precision, intent(in) :: t2(6,6)
    double precision :: t4(3,3,3,3)

    t4(1,1,1,1) = t2(1,1)
    t4(1,1,2,2) = t2(1,2)
    t4(1,1,3,3) = t2(1,3)
    t4(1,1,2,3) = t2(1,4)/sqrt2
    t4(1,1,1,3) = t2(1,5)/sqrt2
    t4(1,1,1,2) = t2(1,6)/sqrt2
    t4(1,1,3,2) = t4(1,1,2,3)
    t4(1,1,3,1) = t4(1,1,1,3)
    t4(1,1,2,1) = t4(1,1,1,2)

    t4(2,2,1,1) = t2(2,1)
    t4(2,2,2,2) = t2(2,2)
    t4(2,2,3,3) = t2(2,3)
    t4(2,2,2,3) = t2(2,4)/sqrt2
    t4(2,2,1,3) = t2(2,5)/sqrt2
    t4(2,2,1,2) = t2(2,6)/sqrt2
    t4(2,2,3,2) = t4(2,2,2,3)
    t4(2,2,3,1) = t4(2,2,1,3)
    t4(2,2,2,1) = t4(2,2,1,2)

    t4(3,3,1,1) = t2(3,1)
    t4(3,3,2,2) = t2(3,2)
    t4(3,3,3,3) = t2(3,3)
    t4(3,3,2,3) = t2(3,4)/sqrt2
    t4(3,3,1,3) = t2(3,5)/sqrt2
    t4(3,3,1,2) = t2(3,6)/sqrt2
    t4(3,3,3,2) = t4(3,3,2,3)
    t4(3,3,3,1) = t4(3,3,1,3)
    t4(3,3,2,1) = t4(3,3,1,2)

    t4(2,3,1,1) = t2(4,1)/sqrt2
    t4(2,3,2,2) = t2(4,2)/sqrt2
    t4(2,3,3,3) = t2(4,3)/sqrt2
    t4(2,3,2,3) = t2(4,4)/2.0
    t4(2,3,1,3) = t2(4,5)/2.0
    t4(2,3,1,2) = t2(4,6)/2.0
    t4(2,3,3,2) = t4(2,3,2,3)
    t4(2,3,3,1) = t4(2,3,1,3)
    t4(2,3,2,1) = t4(2,3,1,2)

    t4(3,2,1,1) = t4(2,3,1,1)
    t4(3,2,2,2) = t4(2,3,2,2)
    t4(3,2,3,3) = t4(2,3,3,3)
    t4(3,2,2,3) = t4(2,3,2,3)
    t4(3,2,1,3) = t4(2,3,1,3)
    t4(3,2,1,2) = t4(2,3,1,2)
    t4(3,2,3,2) = t4(2,3,3,2)
    t4(3,2,3,1) = t4(2,3,3,1)
    t4(3,2,2,1) = t4(2,3,2,1)

    t4(1,3,1,1) = t2(5,1)/sqrt2
    t4(1,3,2,2) = t2(5,2)/sqrt2
    t4(1,3,3,3) = t2(5,3)/sqrt2
    t4(1,3,2,3) = t2(5,4)/2.0
    t4(1,3,1,3) = t2(5,5)/2.0
    t4(1,3,1,2) = t2(5,6)/2.0
    t4(1,3,3,2) = t4(1,3,2,3)
    t4(1,3,3,1) = t4(1,3,1,3)
    t4(1,3,2,1) = t4(1,3,1,2)

    t4(3,1,1,1) = t4(1,3,1,1)
    t4(3,1,2,2) = t4(1,3,2,2)
    t4(3,1,3,3) = t4(1,3,3,3)
    t4(3,1,2,3) = t4(1,3,2,3)
    t4(3,1,1,3) = t4(1,3,1,3)
    t4(3,1,1,2) = t4(1,3,1,2)
    t4(3,1,3,2) = t4(1,3,3,2)
    t4(3,1,3,1) = t4(1,3,3,1)
    t4(3,1,2,1) = t4(1,3,2,1)

    t4(1,2,1,1) = t2(6,1)/sqrt2
    t4(1,2,2,2) = t2(6,2)/sqrt2
    t4(1,2,3,3) = t2(6,3)/sqrt2
    t4(1,2,2,3) = t2(6,4)/2.0
    t4(1,2,1,3) = t2(6,5)/2.0
    t4(1,2,1,2) = t2(6,6)/2.0
    t4(1,2,3,2) = t4(1,2,2,3)
    t4(1,2,3,1) = t4(1,2,1,3)
    t4(1,2,2,1) = t4(1,2,1,2)

    t4(2,1,1,1) = t4(1,2,1,1)
    t4(2,1,2,2) = t4(1,2,2,2)
    t4(2,1,3,3) = t4(1,2,3,3)
    t4(2,1,2,3) = t4(1,2,2,3)
    t4(2,1,1,3) = t4(1,2,1,3)
    t4(2,1,1,2) = t4(1,2,1,2)
    t4(2,1,3,2) = t4(1,2,3,2)
    t4(2,1,3,1) = t4(1,2,3,1)
    t4(2,1,2,1) = t4(1,2,2,1)

  end

  pure function stress_vm(stress) result(stress_vm_out)
    double precision, intent(in) :: stress(6)
    double precision :: stress_vm_out

    stress_vm_out = &
     sqrt(0.5*((stress(1) - stress(2))**2 + (stress(2) - stress(3))**2 + (stress(3) - stress(1))**2 + &
     3.0*(stress(4)**2 + stress(5)**2 + stress(6)**2)))

  end

  pure function strain_vm(strain) result(strain_vm_out)
    double precision, intent(in) :: strain(6)
    double precision :: strain_vm_out

    strain_vm_out = &
     sqrt(2.0/3.0*(((strain(1) - strain(2))**2 + (strain(2) - strain(3))**2 + (strain(3) - strain(1))**2)/3.0 + &
     strain(4)**2 + strain(5)**2 + strain(6)**2))

  end

  pure function sym3333(t3333) result(t3333sym)
    double precision, intent(in) :: t3333(3,3,3,3)
    double precision :: t3333sym(3,3,3,3)
    integer :: i,j,k,l

    do i = 1,3
      do j = 1,3
        do k = 1,3
          do l = 1,3
            t3333sym(i,j,k,l) = (t3333(i,j,k,l) + t3333(j,i,k,l) + t3333(i,j,l,k) + t3333(j,i,l,k))/4.0
          enddo
        enddo
      enddo
    enddo

  end

end module tensor_functions
