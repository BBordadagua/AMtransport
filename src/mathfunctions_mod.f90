!!
!!                 data_handling_mod.f90
!!
!!
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!

MODULE math_functions
  use parameters
  use data_handling
  implicit none
  
  !define public subroutines and functions
  public :: spline,splint,interpolate
  public :: ludcmp,lubksb,inverse_matrix,tridag,invert_tridagmatrix
  public :: bksub,difeq,pinvs,red,solvde !solve Relaxation scheme
  public :: multiply_matrix_vector_triag
  
  contains

!cubic spline Numerical recipies book FORTRAN 77

!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f (xi), with
!x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the inter-
!polating function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!length n which contains the second derivatives of the interpolating function at the tabulated
!points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set
!the corresponding boundary condition for a natural spline, with zero second derivative on
!that boundary.

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
  REAL(DP), INTENT(IN) :: yp1,ypn
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: y2
  INTEGER :: n
  INTEGER :: i,k
  REAL(DP) :: p,qn,sig,un
  REAL(DP), DIMENSION(n) :: u

  if (yp1 > 0.99e30) then !The lower boundary condition is set either to be “natural”
    y2(1)=0.
    u(1)=0.
  else !or else to have a specified first derivative.
    y2(1)=-0.5
    u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  end if

  do i=2,n-1 !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do

  if (ypn > 0.99e30) then !The upper boundary condition is set either to be “natural”
    qn=0.
    un=0.
  else !or else to have a specified first derivative.
    qn=0.5
    un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  end if

  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

  do k=n-1,1,-1 !This is the backsubstitution loop of the tridiagonal algorithm.
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do

END SUBROUTINE spline

!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
!xai’s in order), and given the array y2a(1:n), which is the output from spline above,
!and given a value of x, this routine returns a cubic-spline interpolated value y.
function splint(xa,ya,y2a,n,x) result(y)
  IMPLICIT NONE
  REAL(DP), INTENT(IN), DIMENSION(:) :: xa,ya,y2a
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: y
  INTEGER, INTENT(IN) :: n
  INTEGER(kind = 4) :: khi,klo, k
  REAL(DP) :: a,b,h

  klo=1 

  !We will find the right place in the table by means of bisection.
  !This is optimal if sequential calls to this routine are at random
  !values of x. If sequential calls are in order, and closely
  !spaced, one would do better to store previous values of
  !klo and khi and test if they remain appropriate on the next call.

  khi=n
  1 if (khi-klo>1) then !value(i)
      k=(khi+klo)/2
      if(xa(k)>x)then
        khi=k
      else
        klo=k
      end if
  go to 1
  end if !klo and khi now bracket the input value of x.

  h=xa(khi)-xa(klo) !read smallest file
  
  if (h == 0.) then
    !print*,xa(khi),khi,xa(klo),klo
    !write(9000,*) xa
    stop
  end if
  !call nrerror('bad xa input in splint') !pause 'bad xa input in splint' !The xa’s must be distinct.

  a=(xa(khi)-x)/h                             !Cubic spline polynomial is now evaluated.
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+ ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

END function splint

function interpolate(xa,ya,n,x,m) result(y)
  implicit none
  REAL(DP), INTENT(IN) :: xa(:), ya(:)
  REAL(DP), INTENT(IN) :: x(:)
  REAL(DP) :: y(m)
  REAL(DP), allocatable :: y2a(:)
  INTEGER, INTENT(IN) :: n,m
  INTEGER :: j

  allocate(y2a(n))

  call spline(xa(:),ya(:),n,real(0.,16),real(0.,16),y2a(:))
  do j=1,m
    y(j) = splint(xa(:),ya(:),y2a(:),n,x(j))
  end do

  deallocate(y2a)

end function


!Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces 
!it by the LU decomposition of a rowwise permutation of itself. a and n are input. 
!a is output, arranged as in equation (2.3.14) above; indx(1:n) is an output vector 
!that records the row permutation effected by the partial pivoting; d is output as +-1
! depending on whether the number of row interchanges was even or odd, respectively. 
!This routine is used in combination with lubksb to solve linear equations or invert a matrix.
SUBROUTINE ludcmp(a,n,np,indx,d)
  IMPLICIT NONE
  INTEGER, intent(in) :: n,np
  INTEGER, intent(inout) :: indx(n)
  REAL (DP), intent(inout) :: a(np,np),d 
  INTEGER :: i,imax,j,k!,NMAX 
  REAL (DP) :: TINY
  PARAMETER (TINY=1.0e-20)  !Largest expected n, and a small number.
  !PARAMETER (NMAX=10000,TINY=1.0e-20)  !Largest expected n, and a small number.
  REAL (DP) :: aamax,dum,sum,vv(n)!vv(NMAX) !vv stores the implicit scaling of each row. 

  imax = 0 !! changed

  d=1. !No row interchanges yet.
  do i=1,n !Loop over rows to get the implicit scaling information. 
    aamax=0. 
    do j=1,n 
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j)) 
    enddo
    if (aamax.eq.0.) then!pause 'singular matrix in ludcmp' !No nonzero largest element. 
      stop
    end if
    vv(i)=1./aamax !Save the scaling. 
  enddo
  do j=1,n !This is the loop over columns of Crout's method. 
    do i=1,j-1 !This is equation (2.3.12) except for i = j. 
      sum=a(i,j) 
      do k=1,i-1 
        sum=sum-a(i,k)*a(k,j) 
      enddo 
      a(i,j)=sum 
    enddo 
    aamax=0. !Initialize for the search for largest pivot element. 
    do i=j,n !This is i = j of equation (2.3.12) and i = j +1:::N of equation (2.3.13). 
      sum=a(i,j) 
      do k=1,j-1 
        sum=sum-a(i,k)*a(k,j) 
      enddo 
      a(i,j)=sum 
      dum=vv(i)*abs(sum) !Figure of merit for the pivot. 
      if (dum.ge.aamax) then !Is it better than the best so far? 
        imax=i 
        aamax=dum 
      endif 
    enddo
    if (j.ne.imax) then !Do we need to interchange rows? 
      do k=1,n !Yes, do so... 
        dum=a(imax,k) 
        a(imax,k)=a(j,k) 
        a(j,k)=dum 
      enddo 
      d=-d           !...and change the parity of d. 
      vv(imax)=vv(j) !Also interchange the scale factor.
    endif
    indx(j)=imax 
    if(a(j,j).eq.0.)a(j,j)=TINY
    !If the pivot element is zero the matrix is singular 
    !(at least to the precision of the algorithm). 
    !For some applications on singular matrices, it is desirable to substitute 
    !TINY for zero.

    if(j.ne.n)then !Now, Finally, divide by the pivot element. 
      dum=1./a(j,j) 
      do i=j+1,n 
        a(i,j)=a(i,j)*dum 
      enddo 
    endif 
  enddo !Go back for the next column in the reduction. 
  return

END SUBROUTINE ludcmp

!Solves the set of n linear equations A * X = B.Hereais input, not as the matrix A
!but rather as its LU decomposition, determined by the routine ludcmp. indx is 
!input as the permutation vector returned by ludcmp. b(1:n) is input as the 
!right-hand side vector B, and returns with the solution vector X. a, n, np,and indx
!are not modified by this routine and can be left in place for successive calls 
!with different right-hand sides b. This routine takes into account the possibility
!that b will begin with many zero elements, so it is eficient for use in matrix inversion.
SUBROUTINE lubksb(a,n,np,indx,b) 
  IMPLICIT NONE
  INTEGER, intent(in) :: n,np,indx(n) 
  REAL (DP), intent(in) :: a(np,np)
  REAL (DP), intent(inout) :: b(n)
  INTEGER :: i,ii,j,ll 
  REAL (DP) :: sum
  
  !When ii is set to a positive value, it will become the index of the first 
  !nonvanishing element of b.Wenowdo the forward substitution, equation (2.3.6). 
  !The only new wrinkle is to unscramble the permutation as we go
  ii=0
  do i=1,n 
    ll=indx(i) 
    sum=b(ll) 
    b(ll)=b(i) 
    if (ii.ne.0)then 
      do j=ii,i-1 
        sum=sum-a(i,j)*b(j) 
      enddo 
    else if (sum.ne.0.) then
      ii=i !A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
    endif
    b(i)=sum 
  enddo 
  do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7). 
    sum=b(i) 
    do j=i+1,n 
      sum=sum-a(i,j)*b(j) 
    enddo 
    b(i)=sum/a(i,i) !Store a component of the solution vector X. 
  enddo 
  return

END SUBROUTINE lubksb

!Solves for a vector u(1:n) of length n the tridiagonal linear set given 
!by equation (2.4.1). a(1:n), b(1:n), c(1:n),andr(1:n) are input vectors 
!and are not modified. 
SUBROUTINE tridag(a,b,c,r,u,n) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: J!,NMAX
  REAL (DP), INTENT(IN) :: a(:),b(:),c(:),r(n)
  REAL (DP), INTENT(INOUT) :: u(n)
  !PARAMETER (NMAX=5000) !Parameter: NMAX is the maximum expected value of n. 
  REAL (DP) :: bet,gam(n)!,gam(NMAX) !One vector of workspace, gam is needed. 
  
  if(b(1).eq.0.) then
    print*, 'tridag: rewrite equations' 
    stop
  end if
  !If this happens then you should rewrite your equations as a set of order 
  !N − 1,with u2 trivially eliminated. 
  
  bet=b(1) 
  u(1)=r(1)/bet 
  do j=2,n !Decomposition and forward substitution. 
    gam(j)=c(j-1)/bet 
    bet=b(j)-a(j)*gam(j) 
    if(bet.eq.0.) then
      print*, 'tridag failed' !Algorithm fails; see below. 
      stop
    end if
    u(j)=(r(j)-a(j)*u(j-1))/bet 
  enddo 
  do j=n-1,1,-1 !Backsubstitution. 
    u(j)=u(j)-gam(j+1)*u(j+1) 
  enddo 
  return 
END SUBROUTINE tridag

!matrix*u = r , solution is u
SUBROUTINE invert_tridagmatrix(matrix,n,u,r)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL (DP), INTENT(IN) :: r(:),matrix(:,:) !matrix(n,n)
  REAL (DP), INTENT(INOUT) :: u(:)
  REAL (DP), allocatable :: a(:),b(:),c(:)
  INTEGER :: i

  allocate(a(n),b(n),c(n))
  !do i=1,n
  !  b(i) = matrix(i,i)
  !  if (i/=1) then
  !    a(i) = matrix(i-1,i)
  !  end if
  !  if (i/=n) then
  !    c(i) = matrix(i+1,i)
  !  end if
  !end do

  do i=1,n
    b(i) = matrix(2,i)
    if (i/=1) then
      a(i) = matrix(1,i)
    end if
    if (i/=n) then
      c(i) = matrix(3,i)
    end if
  end do

  call tridag(a,b,c,r,u,n) 

  deallocate(a,b,c)


END SUBROUTINE invert_tridagmatrix


!!using LU Decomposition
SUBROUTINE inverse_matrix(a,n,np,y)
  IMPLICIT NONE
  INTEGER, allocatable :: indx(:)
  INTEGER, intent(in):: n,np
  REAL (DP), intent(inout) :: a(np,np),y(np,np)
  REAL (DP) :: d
  INTEGER i,j
  !...

  allocate(indx(np))

  do i=1,n !Set up identity matrix. 
    do j=1,n 
      y(i,j)=0. 
    enddo 
    y(i,i)=1. 
  enddo 
  call ludcmp(a,n,np,indx,d) !Decompose the matrix just once. 
  do j=1,n !Find inverse by columns. 
    call lubksb(a,n,np,indx,y(1,j)) !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the address of the jth column of y. 
  enddo

  deallocate(indx)

  END SUBROUTINE inverse_matrix


  !USES bksub,difeq,pinvs,red
  !Driver routine for solution of two point boundary value problems by relaxation. 
  !itmax is the maximum number of iterations. conv is the convergence criterion 
  !(see text). slowc controls the fraction of corrections actually used after each 
  !iteration. scalv(1:nyj) contains typical sizes for each dependent variable, used 
  !to weight errors. indexv(1:nyj) lists the column ordering of variables used to 
  !construct the matrix s of derivatives. (The nb boundary conditions at the first 
  !mesh point must contain some dependence on the first nb variables listed in indexv.)
  !The problem involves ne equations for ne adjustable dependent variables at each 
  !point. At the first mesh point there are nb boundary conditions. There are a total
  !of m mesh points. y(1:nyj,1:nyk) is the two-dimensional array that contains the
  !initial guess for all the dependent variables at each mesh point. On each 
  !iteration, it is updated by the calculated correction. The arrays 
  !c(1:nci,1:ncj,1:nck), s(1:nsi,1:nsj) supply dummy storage used by the relaxation 
  !code; the minimum dimensions must satisfy: nci=ne, ncj=ne-nb+1, nck=m+1, nsi=ne, nsj=2*ne+1.
  SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m, &
    & y,c,s,MESAfile,oldMESAfile,delta_t,nu_diff)
    implicit none 
    real (DP), intent(in) :: MESAfile(:,:), delta_t, oldMESAfile(:,:),nu_diff
    real (DP), intent(in) :: conv,slowc,scalv(:)
    real (DP), intent(inout) :: c(:,:,:),s(:,:),y(:,:)
    integer, intent(in) :: itmax,m,ne,nb,indexv(:)
    integer :: NMAX 
    parameter (NMAX=10) !Largest expected value of ne.
    integer :: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8, &
    & j9,jc1,jcf,jv,k,k1,k2,km,kp,nvars,kmax(NMAX)
    real (DP) :: err,errj,fac,vmax,vz,ermax(NMAX)


    k1=1 !Set up row and column markers. 
    k2=m 
    nvars=ne*m 
    j1=1 
    j2=nb 
    j3=nb+1 
    j4=ne 
    j5=j4+j1 
    j6=j4+j2 
    j7=j4+j3 
    j8=j4+j4 
    j9=j8+j1 
    ic1=1 
    ic2=ne-nb 
    ic3=ic2+1 
    ic4=ne 
    jc1=1 
    jcf=ic3
    do it=1,itmax 
      !Primary iteration loop. 
      k=k1 !Boundary conditions at first point. 
      call difeq(k,k1,k2,j9,indexv,s,y,MESAfile,oldMESAfile,delta_t,nu_diff) 
      call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s) 
      do k=k1+1,k2 
        !Finite difference equations at all point pairs. 
        kp=k-1 
        call difeq(k,k1,k2,j9,indexv,s,y,MESAfile,oldMESAfile,delta_t,nu_diff) 
        call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s) 
        call pinvs(ic1,ic4,j3,j9,jc1,k,c,s) 
      enddo 
      k=k2+1 !Final boundary conditions. 
      call difeq(k,k1,k2,j9,indexv,s,y,MESAfile,oldMESAfile,delta_t,nu_diff) 
      call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s) 
      call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
      call bksub(ne,nb,jcf,k1,k2,c) !Backsubstitution. 
      err=0. 
      do j=1,ne !Convergence check, accumulate average error. 
        jv=indexv(j) 
        errj=0. 
        km=0 !! this can't be zero
        vmax=0. 
        do k=k1,k2 !Find point with largest error, for each dependent variable. 
          vz=abs(c(jv,1,k)) 
          if(vz.gt.vmax) then 
            vmax=vz 
            km=k 
          endif 
          errj=errj+vz 
        enddo 
        err=err+errj/scalv(j) !Note weighting for each dependent variable. 
        ermax(j)=c(jv,1,km)/scalv(j) 
        kmax(j)=km 
      enddo 
      err=err/nvars 
      fac=slowc/max(slowc,err) !Reduce correction applied when error is large. 
      do j=1,ne !Apply corrections. 
        jv=indexv(j) 
        do k=k1,k2 
          y(j,k)=y(j,k)-fac*c(jv,1,k) 
        enddo
      enddo 
      !write(*,100) it,err,fac 
      write(*,100) 'it',it,'   err',err,'   fac',fac 
      !Summary of corrections for this step. Point with largest error for each 
      !variable can be monitored by writing out kmax and ermax. 
      if(err.lt.conv) return 
    enddo
    !100 format(1x,i4,2f12.6)
    100 format(A3,i4,A6,ES14.4,A6,ES14.3)

    print*, 'itmax exceeded in solvde' !Convergence failed. 
    stop
    
    return

  END SUBROUTINE solvde

  !Backsubstitution, used internally by solvde. 
  SUBROUTINE bksub(ne,nb,jf,k1,k2,c)
    IMPLICIT NONE 
    integer, intent(in) :: ne,nb,k1,k2,jf
    integer :: i,im,j,k,kp,nbf
    real (DP), intent(inout) :: c(:,:,:) !c(nci,ncj,nck) 
    real (DP) :: xx
    
    nbf=ne-nb 
    im=1 
    do k=k2,k1,-1 
    !Use recurrence relations to eliminate remaining dependences. 
      if (k.eq.k1) im=nbf+1 !Special handling of first point. 
      kp=k+1 
      do j=1,nbf 
        xx=c(j,jf,kp) 
        do i=im,ne 
          c(i,jf,k)=c(i,jf,k)-c(i,j,k)*xx 
        enddo
      enddo 
    enddo 
    do k=k1,k2 !Reorder corrections to be in column 1. 
      kp=k+1 
      do i=1,nb 
        c(i,1,k)=c(i+nbf,jf,k) 
      enddo 
      do i=1,nbf 
        c(i+nb,1,k)=c(i,jf,kp) 
      enddo
    enddo 
    return 
  END SUBROUTINE bksub

  !Diagonalize the square subsection of the s matrix, and store 
  !the recursion coefficients in c; used internally by solvde. 
  SUBROUTINE pinvs(ie1,ie2,je1,jsf,jc1,k,c,s)
    IMPLICIT NONE
    integer, intent(in) :: ie1,ie2,je1,jsf,jc1,k
    real (DP), intent(inout) :: c(:,:,:),s(:,:)!c(nci,ncj,nck),s(nsi,nsj)
    INTEGER :: NMAX 
    PARAMETER (NMAX=10) 
    INTEGER :: i,icoff,id,ipiv,irow,j,jcoff,je2,jp,jpiv,js1,indxr(NMAX) 
    REAL (DP) :: big,dum,piv,pivinv,pscl(NMAX)

    ipiv=0  !! added
    jpiv=0  !! added
    jp = 0  !! added

    je2=je1+ie2-ie1 
    js1=je2+1 
    do i=ie1,ie2 
      !Implicit pivoting, as in x2.1. 
      big=0. 
      do j=je1,je2 
        if(abs(s(i,j)).gt.big) big=abs(s(i,j)) 
      enddo 
      if(big.eq.0.) then
        print*, 'singular matrix, row all 0 in pinvs' 
        stop
      end if
      pscl(i)=1./big 
      indxr(i)=0 
    enddo  
    do id=ie1,ie2 
      piv=0. 
      do i=ie1,ie2 
        !Find pivot element. 
        if(indxr(i).eq.0) then 
          big=0. 
          do j=je1,je2 
            if(abs(s(i,j)).gt.big) then 
              jp=j 
              big=abs(s(i,j)) 
            endif 
          enddo 
          if(big*pscl(i).gt.piv) then 
            ipiv=i 
            jpiv=jp 
            piv=big*pscl(i) 
          endif 
        endif 
      enddo 
      if(s(ipiv,jpiv).eq.0.) then
        print*, 'singular matrix in pinvs' 
        stop
      end if
      indxr(ipiv)=jpiv !In place reduction. Save column ordering. 
      pivinv=1./s(ipiv,jpiv) 
      do j=je1,jsf 
        !Normalize pivot row. 
        s(ipiv,j)=s(ipiv,j)*pivinv 
      enddo 
      s(ipiv,jpiv)=1. 
      do i=ie1,ie2 
        !Reduce nonpivot elements in column. 
        if(indxr(i).ne.jpiv) then 
          if(s(i,jpiv).ne.0.) then 
            dum=s(i,jpiv) 
            do j=je1,jsf 
              s(i,j)=s(i,j)-dum*s(ipiv,j) 
            enddo  
            s(i,jpiv)=0. 
          endif 
        endif 
      enddo 
    enddo 
    jcoff=jc1-js1 !Sort and store unreduced coefficients. 
    icoff=ie1-je1 
    do i=ie1,ie2 
      irow=indxr(i)+icoff 
      do j=js1,jsf 
        c(irow,j+jcoff,k)=s(i,j) 
      enddo
    enddo
    write(1000,*) c(:,1,1)
    return 
  END SUBROUTINE pinvs


  !Reduce columns jz1-jz2 of the s matrix, using previous results as 
  !stored in the c matrix. Only columns jm1-jm2,jmf are affected by 
  !the prior results. red is used internally by solvde.
  SUBROUTINE red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,s)
    IMPLICIT NONE 
    integer, intent(in) :: ic1,iz1,iz2,jc1,jcf,jm1,jm2,jmf,jz1,jz2,kc
    real (DP), intent(in) :: c(:,:,:)
    real (DP), intent(inout) :: s(:,:)
    integer :: i,ic,j,l,loff 
    real (DP) :: vx 
     
    loff=jc1-jm1 
    ic=ic1 
     
    do j=jz1,jz2 
      !Loop over columns to be zeroed. 
      do l=jm1,jm2 
        !Loop over columns altered. 
        vx=c(ic,l+loff,kc) 
        do i=iz1,iz2 
          !Loop over rows. 
          s(i,l)=s(i,l)-s(i,j)*vx 
        enddo 
      enddo 
      vx=c(ic,jcf,kc) 
      do i=iz1,iz2 
        !Plus final element. 
        s(i,jmf)=s(i,jmf)-s(i,j)*vx 
      enddo 
      ic=ic+1 
    enddo 
    return 
  END SUBROUTINE red


  !! Rewrite for the specific problem !!
  !The only information returned from difeq to solvde is the matrix of derivatives 
  !s(i,j); all other arguments are input to difeq and should not be altered. 
  !k indicates the current mesh point, or block number. k1,k2 label the first and 
  !last point in the mesh. If k=k1 or k>k2, the block involves the boundary conditions
  !at the first or final points; otherwise the block acts on FDEs coupling variables
  !at points k-1, k.
  !Returns matrix s(i,j) for solvde. 
  SUBROUTINE difeq(k,k1,k2,jsf,indexv,s,y,MESAfile,oldMESAfile,delta_t,nu_diff)
    IMPLICIT NONE 
    REAL (DP), intent(in) :: MESAfile(:,:), delta_t, oldMESAfile(:,:),nu_diff,y(:,:)
    INTEGER, intent(in) :: jsf,k,k1,k2,indexv(:)
    REAL (DP), intent(inout) :: s(:,:)
    !INTEGER :: mm,n,i
    real (DP) :: delta_nu,radius,rho,nu,radius_old,omega_old,g_surf!,r_dot
    real (DP) :: C1,C2,C3,C4,C5,C6
     

    if(k.eq.k1) then !Boundary condition at first point. 
      delta_nu   = MESAfile(7,k) !!! caution division by zero!!!
      radius     = MESAfile(2,k)
      rho        = MESAfile(3,k)
      g_surf     = MESAfile(6,k)
      nu         = MESAfile(7,k) 
      radius_old = oldMESAfile(2,k)
      omega_old  = oldMESAfile(4,k)
      C6 = (4.*PI*MESAfile(6,k)*MESAfile(3,k)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      
    !else if(k.gt.k2) then !Boundary conditions at last point. 
    else !Interior point. 
      delta_nu   = MESAfile(7,k) - MESAfile(7,k-1)
      radius     = (MESAfile(2,k) + MESAfile(2,k-1))/2.
      rho        = (MESAfile(3,k) + MESAfile(3,k-1))/2.
      g_surf     = (MESAfile(6,k) + MESAfile(6,k-1))/2.
      nu         = (MESAfile(7,k) + MESAfile(7,k-1))/2.
      radius_old = (oldMESAfile(2,k) + oldMESAfile(2,k-1))/2.
      omega_old  = (oldMESAfile(4,k) + oldMESAfile(4,k-1))/2.
      C6 = (4.*PI*MESAfile(6,k-1)*MESAfile(3,k-1)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
    end if
    C1 = radius**2 /delta_t
    C2 = radius_old**2 *omega_old /delta_t
    C3 = 0. !MESAfile(4,i)/(rho*R_sun**2)  !Jmodes
    C4 = (16.*PI*R_sun**4 *radius**4 *rho)/(9.*M_sun*g_surf*sqrt(nu)*delta_nu)
    C5 = (4.*PI*MESAfile(6,k)*MESAfile(3,k)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
    
    if(k.eq.k1) then !Boundary condition at first point. 
      s(1,2+indexv(1))=0.
      s(1,2+indexv(2))=1.
      s(1,jsf)=y(2,k1)
      !call printMatrix(s,2,5,1000)
      print*,s(1,jsf), jsf
    else if(k.gt.k2) then !Boundary conditions at last point. 
      s(2,2+indexv(1))=0.
      s(2,2+indexv(2))=1.
      s(2,jsf)=y(2,k2)
      call printMatrix(s,2,5,2000)
    else !Interior point. 
      s(1,indexv(1))=0.5d0*C1  !E1
      s(1,indexv(2))=C6
      s(1,2+indexv(1))= 0.5d0*C1
      s(1,2+indexv(2))=-C5

      s(2,indexv(1))=C4    !E2
      s(2,indexv(2))=0.5d0 
      s(2,2+indexv(1))=-C4
      s(2,2+indexv(2))=0.5d0

      s(1,jsf)=0.5d0*C1*(y(1,k)+y(1,k-1))-C2-C5*y(2,k)+C6*y(2,k-1)&
      &-0.5d0*C3*(y(1,k)+y(1,k-1)) !E1
      s(2,jsf)=0.5d0*(y(2,k)+y(2,k-1))-C4*(y(1,k)-y(1,k-1)) !E2
      call printMatrix(s,2,5,3000)
    endif 

    
    return 

  END SUBROUTINE difeq

  !! multiplication of matrix and vector v with a tridiagonal matrix
  !! input are 3 arrays that define the diagonal d, below diagonal b and above diagonal a of matrix
  function multiply_matrix_vector_triag(matrix_triag,v,m) result(r)
    implicit none
    real (DP), intent(in) :: v(:),matrix_triag(:,:)
    integer, intent(in) :: m
    real (DP) :: r(m)
    real (DP), allocatable :: d(:),a(:),b(:)
    integer :: i

    allocate(d(m),a(m),b(m))
    b(:) = matrix_triag(1,:)
    d(:) = matrix_triag(2,:)
    a(:) = matrix_triag(3,:)

    !! boundary conditions
    r(1) = d(1)*v(1) + a(1)*v(2)
    r(m) = b(m)*v(m-1) + d(m)*v(m)

    do i=2,m-1
      r(i) = b(i)*v(i-1) + d(i)*v(i) + a(i)*v(i+1)
    end do

    deallocate(a,b,d)

  end function multiply_matrix_vector_triag



END MODULE math_functions