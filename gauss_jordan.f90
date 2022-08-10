PROGRAM gauss_jordan

!***********************************************************************
!----    This program is the Gauss Jordan method to diagonalizate   ----
!----    matrices                                                   ----
!----                        ~VCastor 2020                          ----
!***********************************************************************

IMPLICIT NONE
INTEGER                           :: i, j, k, r, c, l, option
REAL                              :: aux, aux2
REAL, DIMENSION(:), ALLOCATABLE   :: V
REAL, DIMENSION(:,:), ALLOCATABLE :: M
CHARACTER(LEN=50)                 :: archivo

WRITE(*,*) "The equation system needs to be in form as a matrix."
WRITE(*,*) "This program needs the matrix, you have a file &
            with the data or you will write down now?"
WRITE(*,*) "Write down : 0"
WRITE(*,*) "Data file : 1"
WRITE(*,FMT='(A)',ADVANCE='NO') 'Select one option: '
READ(*,*) option

!-----------------------------------------------------------------------
!---- How the data will be read
SELECTCASE(option)
  CASE(0)
    WRITE(*,*) "Write down the number of variables of the system, &
                the matrix will have that number of rows"
    WRITE(*,FMT='(A)',ADVANCE='NO') "Variables: "
    READ(*,*) r
    c = r + 1                  !columns
    ALLOCATE(M(r,c),V(c))      !dimension of the matrix

    WRITE(*,*) 'Element of the matrix'
    DO j = 1, c
      DO i = 1, r
        WRITE(*,FMT='(A1,I1,A1,I1,A3)',ADVANCE='NO') '(',i,',',j,'): '
        READ(*,*) M(i,j)
      ENDDO
    ENDDO

  CASE(1)
    WRITE(*,*) "The file would have the number of rows and columsn &
                in the first line, then the matrix"
    WRITE(*,FMT='(A)',ADVANCE='NO') "File name of matrix data: "
    READ(*,*) archivo
    OPEN(11,FILE=archivo)
      READ(11,*) r, c
      ALLOCATE(M(r,c),V(c))
      DO i = 1, r
        READ(11,*) M(i,:)
      ENDDO
    CLOSE(11)
ENDSELECT

!---- Write the system
WRITE(*,*) '**************************************'
WRITE(*,*) "Ok, the system is:"
DO i = 1, r
  WRITE(*,*) (M(i,j), j=1,c)
ENDDO

!---------------------------------------------------------
!---- Now, here we go with the Gauss Jordan method

!---- First up diagonal
DO i = 1, r-1           !for all rows except the las one
  aux = M(i,i)          !diagonal elemnt
  DO j = 1, c           !all elements in the row
    V(j) = M(i,j)/aux   !divided by the diagonal elemnt
  ENDDO
  DO k = i+1, r         !only with the rows under itself
    aux2 = M(k,i)
    DO j = 1, c         !Actualize the matrix j-times
      M(k,j) = M(k,j) - V(j)*aux2
    ENDDO
  ENDDO
ENDDO

WRITE(*,*) '**************************************'
WRITE(*,*) '*************Firsts steps*************'
DO i = 1, r
  WRITE(*,*) (M(i,j), j=1,c)
ENDDO

!---- 2nd zeros in up diagonal
DO i = 1, r-1                      !for all rowa except the last one
  DO l = i+1, r                    !only with the rows under itself
    aux = M(l,l)
    DO k = 1, c                    !all elemnts in the row
      V(k) = M(l,k)/aux
    ENDDO
    aux2 = M(i,l)
    DO j = 1, c                    !all elements in the row
      M(i,j) = M(i,j) - V(j)*aux2
    ENDDO
  ENDDO
ENDDO

!---- 3rd diagonal equals to one
DO i = 1, r
  aux = M(i,i)
  DO j = 1, c
    M(i,j) = M(i,j)/aux
  ENDDO
ENDDO

WRITE(*,*) '***************************************'
WRITE(*,*) 'The matrix with values of the variables'
DO i = 1, r
  WRITE(*,*) (M(i,j), j=1,c)
ENDDO

ENDPROGRAM gauss_jordan
