PROGRAM determinant

!***********************************************************************
!----       This program compute the determinant by LU matrices     ----
!----                           VCastor 2020                        ----
!-----------------------------------------------------------------------
!***********************************************************************

IMPLICIT NONE
INTEGER                                   :: i, j, k, n, option, aux3
REAL(KIND=8)                              :: aux, aux2, up, low, det
REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: V
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: M, L, U, P
CHARACTER(LEN=50)                         :: archivo

WRITE(*,*) "This program needs a matrix, you have a file &
            with the data or you will write down now?"
WRITE(*,*) "Write down : 0"
WRITE(*,*) "Data file : 1"
WRITE(*,FMT='(A)',ADVANCE='NO') 'Select one option: '
READ(*,*) option

!---- How the data will be read
SELECTCASE(option)

  CASE(0)
    WRITE(*,*) "Write down the dimension of the matrix"
    WRITE(*,FMT='(A)',ADVANCE='NO') "Dimension: "
    READ(*,*) n
    ALLOCATE(M(n,n),V(n),L(n,n),U(n,n))

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
      ALLOCATE(M(n,n),V(n),L(n,n),U(n,n))
      DO i = 1, n
        READ(11,*) M(i,:)
      ENDDO
    CLOSE(11)

ENDSELECT

!---- Write the system
WRITE(*,*) '***************************************'
WRITE(*,*) '*******The original matrix is: ********'
DO i = 1, n
    WRITE(*,*) (M(i,j), j=1,n)
ENDDO

!-----------------------------------------------------------------------
!---- LU method

!---- Set the matrices
L(:,:) = 0.d0
U(:,:) = M(:,:)
DO i = 1, n
  aux = U(i,i)                      !diagonal elemnt
  DO j = 1, n
    V(j)   = U(i,j)/aux             !divided by the diagonal elemnt
    U(i,j) = V(j)                   !U actualizaded
    IF (i .EQ. j) THEN
      L(i,i) = aux                  !diagonal for L
    ENDIF
  ENDDO
  DO k = i+1, n
    aux2 = U(k,i)
    DO j = 1, n
      U(k,j) = U(k,j) - V(j)*aux2
    ENDDO
    L(k,i) = aux2
  ENDDO
ENDDO

WRITE(*,*) '***************************************'
WRITE(*,*) '***********U matrix is : **************'
DO i = 1, n
  WRITE(*,*) (U(i,j), j=1,n)
ENDDO

WRITE(*,*) '***************************************'
WRITE(*,*) '***********L matrix is : **************'
DO i = 1, n
  WRITE(*,*) (L(i,j), j=1,n)
ENDDO

ALLOCATE(P(n,n))
P = MATMUL(L,U)
WRITE(*,*) '***************************************'
WRITE(*,*) '***********U \times L is: *************'
DO i = 1, n
  WRITE(*,*) (P(i,j), j=1,n)
ENDDO

!-----------------------------------------------------------------------
!---- determinants

!---- for U
up = 1.d0 low = 1.d0
DO i = 1, n
  up  = up*U(i,i)
  low = low*L(i,i)
ENDDO
!---- finally
det = up*low

WRITE(*,*) '***************************************'
WRITE(*,*) '   The determiant of the matrix is: '
WRITE(*,*) det
WRITE(*,*) '***************************************'


ENDPROGRAM determinant
