!!
!!                 data_handling_mod.f90
!!
!!
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!

MODULE data_handling
  use parameters
  !use math_functions
  implicit none
  
  !define public subroutines and functions
  public :: getnlines,readfile,printarray,printMatrix,getradiativeR
  public :: getsmallestfile,string, getj_l0,runGYRE,double_Npoints
  
  contains

function string(g) result(str)
    !   "Convert an integer to string."
      implicit none
        integer, intent(in) :: g
        character(len=20) :: str
        write (str, *) g
        str = adjustl(str)
end function string

function getnlines(filename) result(nlines)
  implicit NONE
  integer :: nlines, io
  character(len=*), intent(in) :: filename
  
  nlines = 0 
  open(101, action='read',file=trim(filename))
  DO
    READ(101,*,iostat=io)
    IF (io/=0) EXIT
    nlines = nlines + 1
  END DO
  close(101)
  nlines = nlines - 6

end function getnlines

subroutine readfile(filename,filearray)
  implicit NONE
  real (DP), intent(inout):: filearray(:,:)
  character(len=*), intent(in) :: filename
  character(len=50) :: dummy
  
  open(100, action='read',file=trim(filename))
    read(100,*) dummy
    read(100,*) dummy
    read(100,*) dummy
    read(100,*) dummy
    read(100,*) dummy
    read(100,*) filearray
  close(100)

end subroutine

function getj_l0(sumarray,freq) result(jvalue)
  implicit NONE
  integer:: m,jvalue,i
  real (DP), intent(in):: sumarray(:,:),freq
  real (DP) :: d, dnew

  m = size(sumarray,dim=2)

  jvalue = 1
  d = 100.
  do i=1,m
    if (int(sumarray(8,i)) == 0) then
      dnew = abs(sumarray(5,i)- freq)
      if (dnew < d) then
        jvalue = int(sumarray(7,i))
        d=dnew
      end if
    end if
  end do
end function getj_l0


subroutine printarray(array, m,name)
  implicit none
  real (DP), intent(in) :: array(m)
  integer, intent(in) :: m,name
  integer :: j
      do j=1,m
        write (name,*) array(j)
      end do
end subroutine printarray


subroutine printMatrix(array, p, m,name)
  implicit none
  real (DP), intent(in) :: array(p,m)
  integer, intent(in) :: p,m,name
  integer :: j
  !open(50, action='write',file='output/Matrix.txt')
      do j=1,m
        write (name,*) array(:,j)
      end do
  !close(50)50
end subroutine printMatrix


subroutine getsmallestfile(sumfile,jarray,jfinal,nlines,model)
  implicit none
  real (DP), intent(in) :: sumfile(:,:)
  integer :: lvalue, nnew,m,i
  integer, intent(in) :: jarray(:),model
  integer, intent(out) :: jfinal, nlines

  m = size(jarray,dim=1)

  nlines = 1e9
  do i=1,m
    lvalue = int(sumfile(8,jarray(i)))
    nnew = getnlines('input/model_'//trim(string(model))//'/detail.l'//trim(string(lvalue))//'.j'//trim(string(jarray(i)))//'.txt') !number of lines
    if (nnew < nlines) then
      nlines = nnew
      jfinal = jarray(i)
    end if

  end do

end subroutine getsmallestfile


subroutine getradiativeR(ri,rf,filearray,value)
  implicit none
  integer, intent(inout) :: ri, rf
  real(DP), intent(in) :: filearray(:,:)
  integer :: N,i, counti
  integer, intent(in) :: value !! defines if we use the GYRE file(0) or the MESA profile(1)

  N = size(filearray,dim=2)

  counti=0
  if (value ==0) then !! for GYRE file
    do i=1,N
      ! for RC star !filearray(2,i)/R_star >0.001
      !if (filearray(9,i) >0. .AND. counti==0 .AND. filearray(2,i)/R_star < 0.6 .AND. filearray(2,i)/R_star >0.001) then
      if (filearray(9,i) >0d0 .AND. counti==0 .AND. filearray(2,i)/R_star < 0.6d0) then
        counti=1
        ri = i!filearray(2,i)
      end if
      ! for RC star !filearray(2,i)/R_star >0.1
      !if (filearray(9,i) <0. .AND. counti==1 .AND. filearray(2,i)/R_star < 0.6 .AND. filearray(2,i)/R_star >0.1) then
      if (filearray(9,i) <0d0 .AND. counti==1 .AND. filearray(2,i)/R_star < 0.6d0) then
        rf = i-1!filearray(2,i-1)
        counti = 2
      end if
    end do
  else if (value == 1) then !! for MESA file
    do i=1,N
     if (filearray(31,i) >0d0 .AND. counti==0 .AND. (10**filearray(3,i))*R_sun/R_star < 0.6d0) then
        counti=1
        ri = i
      end if
      if (filearray(31,i) <0d0 .AND. counti==1 .AND. (10**filearray(3,i))*R_sun/R_star < 0.6d0) then
        rf = i-1
        counti = 2
      end if
    end do
  end if

end subroutine getradiativeR




subroutine runGYRE(model)
  implicit none
  integer, intent(in) :: model
  integer :: status,status0
  character(35) :: COMMAND,makedirectory


  !!execute_command_line run GYRE
  makedirectory = 'mkdir -p input/output'
  status0 = SYSTEM(trim(makedirectory))
  COMMAND = '$GYRE_DIR/bin/gyre gyre.in'
  status = SYSTEM(COMMAND)

  makedirectory = 'mv input/output input/model_'//trim(string(model))
  status = SYSTEM(trim(makedirectory))
  if ( status0 .ne. 0 .and. status .ne. 0 ) stop 'system: error'

end subroutine runGYRE


!not changed since GYRE new version
subroutine double_Npoints(file,double_file,N)
  implicit none
  real(DP), intent(in) :: file(:,:)
  real(DP), intent(inout) :: double_file(:,:)
  integer, intent(in) :: N 
  integer :: i,j

  !double number of mass points
  do i=1,N
    double_file(3,2*i -1) = file(3,i) !odd numbers are the same
    if (i/=N) then
      double_file(3,2*i) = (file(3,i) + file(3,i+1))/2.
    end if
  end do

  !only interpolate in radius

  do j=2,19
    if(j/=3) then
      !double_file(j,:) = interpolate(file(3,:),file(j,:),N,double_file(3,:),2*N-1)
    end if
  end do

end subroutine double_Npoints


END MODULE data_handling