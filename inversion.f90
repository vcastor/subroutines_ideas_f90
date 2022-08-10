PROGRAM inverstion

!***********************************************************************
!----      This program compute the inverse matrix by Gaussian      ----
!----                         ~VCastor 2020                         ----
!***********************************************************************

IMPLICIT NONE
INTEGER                                   :: i, j, k, n, l, option, aux3
REAL(KIND=8)                              :: aux, aux2
REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: V
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: M, S
CHARACTER(LEN=50)                         :: archivo

WRITE(*,*) "This program needs a matrix you have a with the data or &
            you will write down now?"
WRITE(*,*) "Write down : 0"
WRITE(*,*) "Data file : 1"
WRITE(*,FMT='(A)',ADVANCE='NO') 'Select one option: '
READ(*,*) option

!-----------------------------------------------------------------------
!---- How the data will be read?
SELECTCASE(option)

  CASE(0)
    WRITE(*,*) "Write down the dimension of the matrix"
    WRITE(*,FMT='(A)',ADVANCE='NO') "Dimension: "
    READ(*,*) n
    ALLOCATE(M(n,n),V(2*n))

    WRITE(*,*) 'Element of the matrix'
    DO j = 1, n
      DO i = 1, n
        WRITE(*,FMT='(A1,I1,A1,I1,A3)',ADVANCE='NO') '(',i,',',j,'): '
        READ(*,*) M(i,j)
      ENDDO
    ENDDO

  CASE(1)
    WRITE(*,FMT='(A)',ADVANCE='NO') "File name of matrix data: "
    READ(*,*) archivo
    OPEN(11,FILE=archivo)
    READ(11,*) n
    ALLOCATE(M(n,n),V(2*n))
    DO i = 1, n
      READ(11,*) M(i,:)
    ENDDO

ENDSELECT

!-----------------------------------------------------------------------
!---- Write the system
DO i = 1, n
  WRITE(*,*) (M(i,j), j=1,n)
ENDDO
WRITE(*,*) '**********************************'

!-----------------------------------------------------------------------
!---- Now, compute the inverse

!---- First, add the identity
ALLOCATE(S(n,2*n))
S(:,:) = 0.
DO i = 1, n
  aux3 = n+i                   !"diagonal" on the right side
  S(i,aux3) = 1.
ENDDO
DO i = 1, n
  DO j = 1, n
    S(i,j) = M(i,j)            !data in the left side
  ENDDO
ENDDO

!---- The matrix with the identity
DO i = 1, n
  WRITE(*,*) (S(i,j), j=1,2*n)
ENDDO

!---- 2nd up diagonal
DO i = 1, n-1
  aux = S(i,i)                 !diagonal elemnt
  DO j = 1, 2*n
    V(j) = S(i,j)/aux          !divided by the diagonal elemnt
  ENDDO
  DO k = i+1, n
    aux2 = S(k,i)
    DO j = 1, 2*n              !Actualize the matrix
      S(k,j) = S(k,j) - V(j)*aux2
    ENDDO
  ENDDO
ENDDO

WRITE(*,*) '***************************************'
DO i = 1, n
  WRITE(*,*) (S(i,j), j=1,2*n)
ENDDO

!---- 3rd diagonal, zeros in the up part of matrix
DO i = 1, n-1
  DO l = i+1, n
    DO k = l, n
      aux = S(k,k)
      DO j = 1, 2*n
        V(j) = S(k,j)/aux
      ENDDO
      aux2 = S(i,k)
      DO j = 1, 2*n
        S(i,j) = S(i,j) - V(j)*aux2
      ENDDO
    ENDDO
  ENDDO
ENDDO

WRITE(*,*) '***************************************'
DO i = 1, n
    WRITE(*,*) (S(i,j), j=1,2*n)
ENDDO

!---- 4th diagonal equals to one
DO i = 1, n
  aux = S(i,i)
  DO j = 1, 2*n
    S(i,j) = S(i,j)/aux
ENDDO
ENDDO

WRITE(*,*) '***************************************'
WRITE(*,*) '*****The inverse in the left side******'
DO i = 1, n
  WRITE(*,*) (S(i,j), j=1,2*n)
ENDDO

!---- only the inverse

DO i = 1, n
  DO j = 1, n
    aux3 = n+j
    M(i,j) = S(i,aux3)
  ENDDO
ENDDO

WRITE(*,*) '***************************************'
WRITE(*,*) '**************The inverse**************'
DO i = 1, n
  WRITE(*,*) (M(i,j), j=1,n)
ENDDO

ENDPROGRAM inverstion 
