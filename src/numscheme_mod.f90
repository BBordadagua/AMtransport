Module numerical_scheme
  use parameters
  use data_handling
  use math_functions
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: rotation_profile_implicit,rotation_profile_explicit,implicit_test
  public :: create_matrix_diffusionCN,create_matrix_diffusionCN_mod,CrankNicolson
  contains

  !!Backward Euler
  subroutine rotation_profile_implicit(F_total,MESAold,omega,model,m)
    implicit none
    real(DP) :: r_n_i1,r_n_i,omega_n_i(m),r_dot,Jacobian(m),f(m),rho_n_i,Tm(m),Tc(m),AMi,AMf
    real (DP), intent(inout) :: omega(:), F_total(:,:),MESAold(:,:)
    real(DP), allocatable :: new_F_total(:,:),MESAfile(:,:),MESAold_int(:,:)
    integer :: nsize,j,ri,rf,mold
    integer, intent(in) :: model,m
    character(len=18) :: dummy

    !allocate arrays
    allocate(MESAfile(19,m))
    open(200, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(200,*) dummy
      read(200,*) MESAfile
    close(200)

    mold = size(MESAold,dim=2) 
    nsize = size(F_total,dim=2)
    allocate(new_F_total(2,m))
    allocate(MESAold_int(4,m))

    !Interpolation of F_total ------------------------------------------------------------
    F_total(2,1) = 0. !making sure it's not nan
    new_F_total(2,:) = interpolate(F_total(3,:),F_total(2,:),nsize,MESAfile(3,:),m)
    new_F_total(1,:) = MESAfile(3,:)  !mass
    
    !interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),m)
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),m) !rho
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),m)
    MESAold_int(1,:) = MESAfile(3,:)
    
    !compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)

    omega(:) = MESAold_int(4,:)
    omega_n_i(:) = MESAold_int(4,:)

    open(250, action='write',file = 'output/rotj50_'//trim(string(model))//'.txt')
    !do i=1,100 !number of iterations of the implicit method
      AMi = 0.
      AMf = 0.

      open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
      do j=1,m
        r_n_i1 = MESAfile(2,j)
        r_n_i = MESAold_int(2,j)
        r_dot = (r_n_i1 - r_n_i)/dt
        rho_n_i = MESAold_int(3,j)
        !omega(j) = MESAold_int(4,j)

        f(j) = - 2.*r_dot/r_n_i + new_F_total(2,j)/(r_n_i**2 *rho_n_i*omega(j))
      
        !f(j) = - 2.*r_dot/r_n_i + new_F_total(3,j)

        Jacobian(j) = 1./dt - f(j)
        Jacobian(j) = f(j)/Jacobian(j)

        !compute timescales
        Tm(j) = abs(rho_n_i*r_n_i**2 *omega(j) /new_F_total(2,j))
        Tc(j) = abs(-r_n_i/r_dot)


        !compute omega
        if (j<ri .OR. j>rf) then
          omega(j) = omega_n_i(j)
        else
          omega(j) = omega_n_i(j) + Jacobian(j)*omega(j)
        end if


        !r_n_i1 = MESAfile(2,j)

        if (j<ri .OR. j>rf) then
        else
          AMi = AMi + MESAold_int(2,j)**2 * omega_n_i(j)
          AMf = AMf + r_n_i1**2 *omega(j)
        end if
      
        write(203,*) r_n_i1/MESAfile(2,rf),omega_n_i(j),Jacobian(j),-f(j),&
        &omega(j),Tm(j),Tc(j),dt
      
        if (j==5) then
          write(250,*) r_n_i1/MESAfile(2,rf),omega(j)
        end if

        omega_n_i(j) = omega(j)
        
      end do
      close(203)
    !end do
    close(250)


    print*, 'delta AM radiative zones = ', AMf - AMi

    deallocate(new_F_total,MESAfile,MESAold_int)

  end subroutine rotation_profile_implicit

  subroutine implicit_test(k,h,x0ini,y0ini)
    implicit none
    !dot_product(vector_a, vector_b) 	This function returns a scalar product of two input vectors, which must have the same length.
    !matmul (matrix_a, matrix_b) 	It returns the matrix product of two matrices
    integer ::  maxSteps,i
    real (DP) :: matrix(2,2), vector_x0(2), dx(2), x_new(2), x0,y0
    real (DP), intent(in) :: k, h,x0ini,y0ini

    maxSteps = 10

    matrix(1,1) = h/(1.+h)
    matrix(1,2) = 0.
    matrix(2,1) = 0.
    matrix(2,2) = h/(1.+k*h)

    x0 = x0ini
    y0 = y0ini

    do i=1,maxSteps
      print*, 'x0 = ', x0, ' y0 = ', y0
      vector_x0(1) = -x0
      vector_x0(2) = -k*y0

      x_new(1) = x0
      x_new(2) = y0

      dx = matmul(matrix,vector_x0) !!shoudl matrix be transposed???

      x_new = x_new + dx

      x0 = x_new(1)
      y0 = x_new(2)
    end do


  end subroutine implicit_test

  subroutine rotation_profile_explicit(MESAold,omega,model,m)
    implicit none
    real (DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,dev_r,dev_t!,Tm,Tc,J_dot,r_dot
    real (DP) :: AMi,AMf, u(m),nu,delta_r
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:)
    integer :: j,ri,rf,mold
    integer, intent(in) :: model,m
    character(len=18) :: dummy

    !allocate arrays
    allocate(MESAfile(19,m))
    open(200, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(200,*) dummy
      read(200,*) MESAfile
    close(200)

    mold = size(MESAold,dim=2) 
    allocate(MESAold_int(4,m))

    
    !interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),m)
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),m)
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),m)
    MESAold_int(1,:) = MESAfile(3,:)
  
    !compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')

    AMi = 0.
    AMf = 0.
    do j=1,m
      r_n_i1 = MESAfile(2,j)
      r_n_i = MESAold_int(2,j)
      rho_n_i = MESAold_int(3,j)
      omega_n_i = MESAold_int(4,j)

      !r_dot = (r_n_i1 - r_n_i)/dt
      !dev_r = 2.*r_dot*omega_n_i*dt/(r_n_i)
      !dev_t = omega_n_i

      !------------------------------------------------------------------
      !!J dot
      !A = 1e10/(3.*dt)
      !J_dot = -(3.*A)*dt/(rho_n_i*r_n_i**2)!new_F_total(2,j)*dt/(rho_n_i*r_n_i**2)
      !!flux of AM
      !A = 1e10/(3.*dt)
      !J_dot = -(dt/(rho_n_i*r_n_i**4))*(A*MESAold_int(2,j+1)**3 - A*r_n_i**3)/(MESAold_int(2,j+1)-r_n_i)
      
      !if (j<ri .OR. j>rf) then
      !  omega(j) = MESAold_int(4,j)
      !else
      !  omega(j) = J_dot - dev_r + dev_t
      !end if

      !compute timescales
      !Tm = abs(rho_n_i*r_n_i**2 *omega_n_i /new_F_total(2,j))
      !Tc = abs(-r_n_i/r_dot)

      !write(203,*) r_n_i1/MESAfile(2,rf),abs(dev_t),abs(J_dot),&
      !&abs(dev_r),omega(j), Tm,Tc,dt!,&
      !!&r_n_i,r_dot,MESAold_int(1,j)
      !-------------------------------------------------------------------

      !---------------- Diffusion eq. ---------------------------------
      if (j == 1) then
        u(j) = omega_n_i*r_n_i**2
        omega(j) = omega_n_i
        delta_r = 0.
      else if (j == m) then
        u(j) = omega_n_i*r_n_i**2
        omega(j) = u(j)/(r_n_i1**2)
        delta_r = 0.
      else
        nu = 1.!1.e14
        delta_r = MESAold_int(2,j+1) - r_n_i
        dev_r = dt*nu*(MESAold_int(3,j+1)*MESAold_int(2,j+1)**4 - rho_n_i*r_n_i**4)*&
        &(MESAold_int(4,j-1)-omega_n_i)/((delta_r**2)*rho_n_i*r_n_i**2) 
        dev_t = dt*nu*(r_n_i**2)*(MESAold_int(4,j+1)-2.*omega_n_i+MESAold_int(4,j-1))/(delta_r**2)

        u(j) = omega_n_i*r_n_i**2 + dev_r + dev_t

        omega(j) = u(j)/(r_n_i1**2)

        write(1000,*)  omega_n_i*r_n_i**2 ,dev_r, dev_t, dt*nu/(delta_r**2)
      end if

      write(203,*) r_n_i1,u(j),omega_n_i*r_n_i**2,0.,omega(j),omega_n_i,dev_r,dev_t,dt,(delta_r**2)/(2.*nu)
      
      !if (j<ri .OR. j>rf) then
      !else
        AMi = AMi + r_n_i**2 * omega_n_i
        AMf = AMf + r_n_i1**2 *omega(j)
      !end if
    end do


    close(203)
    print*, 'delta AM radiative zones = ', AMf - AMi
    !write(500,*)  model, AMi,AMf, AMf - AMi
    
    deallocate(MESAfile,MESAold_int)

  end subroutine rotation_profile_explicit



  subroutine create_matrix_diffusionCN(MESAfile,m,delta_t,nu_diff,matrix_A,matrix_B)
    implicit none
    real (DP), intent(in) :: MESAfile(:,:), delta_t,nu_diff
    real (DP), intent(inout) :: matrix_A(:,:),matrix_B(:,:)
    integer, intent(in) :: m 
    real (DP), allocatable :: lambda(:),alpha(:),beta(:)
    real (DP) :: delta_r,delta_rho!,r_dot
    integer :: i

    allocate(lambda(m),alpha(m),beta(m))

    do i=1,m
     if (i==m) then
        delta_r = MESAfile(2,i) - MESAfile(2,i-1)
        delta_rho = (log(MESAfile(3,i)) - log(MESAfile(3,i-1)))/delta_r
      else
        delta_r = MESAfile(2,i+1) - MESAfile(2,i)
        delta_rho = (log(MESAfile(3,i+1)) - log(MESAfile(3,i)))/delta_r
      end if
      !r_dot = MESAfile(2,i) - MESAold_int(2,i)
      lambda(i) = nu_diff*delta_t/(2.*(delta_r**2)) 
      alpha(i) = (nu_diff*delta_t/(4.*delta_r))*(delta_rho+ 4./MESAfile(2,i))
      beta(i)  = 0.!r_dot/MESAold_int(2,i)
    end do

    !!fill matrix_A
    matrix_A = 0.
    do i=2,m-1
      matrix_A(i,i)   = 1. + 2.*lambda(i) + beta(i)
      matrix_A(i-1,i) = -(lambda(i) - alpha(i))
      matrix_A(i+1,i) = -(lambda(i) + alpha(i))
    end do
    !!Boundary conditions
    matrix_A(2,1)   = -(lambda(1) + alpha(1))
    matrix_A(1,1)   = 1. + lambda(1) + beta(1) + alpha(1)
    matrix_A(m-1,m) = -(lambda(m) - alpha(m))
    matrix_A(m,m)   = 1. + lambda(m) + beta(m) - alpha(m)

    !!fill matrix_B
    matrix_B = 0.
    do i=2,m-1
      matrix_B(i,i)   = 1. - 2.*lambda(i) - beta(i)
      matrix_B(i-1,i) = lambda(i) - alpha(i)
      matrix_B(i+1,i) = lambda(i) + alpha(i)
    end do
    !!Boundary conditions
    matrix_B(2,1)   = lambda(1) + alpha(1)
    matrix_B(1,1)   = 1. - lambda(1) - beta(1) - alpha(1)
    matrix_B(m-1,m) = lambda(m) - alpha(m)
    matrix_B(m,m)   = 1. - lambda(m) - beta(m) + alpha(m)
    !print*,matrix_A(m,m) , matrix_A(m-1,m), matrix_A(1,1),matrix_A(2,1)


    deallocate(lambda,alpha,beta)

  end subroutine create_matrix_diffusionCN


  subroutine create_matrix_diffusionCN_mod(MESAfile,m,delta_t,nu_diff,matrix_A,matrix_B)
    implicit none
    real (DP), intent(in) :: MESAfile(:,:), delta_t,nu_diff
    real (DP), intent(inout) :: matrix_A(:,:),matrix_B(:,:)
    integer, intent(in) :: m 
    real (DP), allocatable :: lambda(:),beta(:)
    integer :: i

    allocate(lambda(m),beta(m))

    do i=1,m
      if (i == 1) then
        lambda(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i+1)*MESAfile(2,i+1)**4)&
          &/((MESAfile(2,i+1)-MESAfile(2,i))*(MESAfile(2,i+1)-MESAfile(2,i)))

        beta(i) = 0.
      else if (i == m) then
        lambda(i) = 0.

        beta(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i-1)*MESAfile(2,i-1)**4)&
          &/((MESAfile(2,i)-MESAfile(2,i-1))*(MESAfile(2,i)-MESAfile(2,i-1)))
      else
        lambda(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i+1)*MESAfile(2,i+1)**4)&
          &/((MESAfile(2,i+1)-MESAfile(2,i))*(MESAfile(2,i+1)-MESAfile(2,i-1)))

        beta(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i-1)*MESAfile(2,i-1)**4)&
          &/((MESAfile(2,i)-MESAfile(2,i-1))*(MESAfile(2,i+1)-MESAfile(2,i-1)))
      end if
    end do

    !!fill matrix_A
    matrix_A = 0.
    do i=2,m-1
      matrix_A(i,i)   = 1. + lambda(i) + beta(i)
      matrix_A(i-1,i) = -beta(i)
      matrix_A(i+1,i) = -lambda(i)
    end do
    !!Boundary conditions
    matrix_A(2,1)   = -lambda(1)
    matrix_A(1,1)   = 1. + lambda(1)
    matrix_A(m-1,m) = -beta(m)
    matrix_A(m,m)   = 1. + beta(m) 

    !!fill matrix_B
    matrix_B = 0.
    do i=2,m-1
      matrix_B(i,i)   = 1. - lambda(i) - beta(i)
      matrix_B(i-1,i) = beta(i)
      matrix_B(i+1,i) = lambda(i)
    end do
    !!Boundary conditions
    matrix_B(2,1)   = lambda(1) 
    matrix_B(1,1)   = 1. - lambda(1)
    matrix_B(m-1,m) = beta(m) 
    matrix_B(m,m)   = 1. - beta(m) 
    
    !print*,matrix_A(m,m) , matrix_A(m-1,m), matrix_A(1,1),matrix_A(2,1),&
    !&matrix_B(m,m) , matrix_B(m-1,m), matrix_B(1,1),matrix_B(2,1)

    deallocate(lambda,beta)

  end subroutine create_matrix_diffusionCN_mod


  function CrankNicolson(MESAfile,m,iter,delta_t,nu_diff) result(omega)
    implicit none
    real (DP), intent(in) :: delta_t, MESAfile(:,:),nu_diff
    integer, intent(in) :: iter,m
    real (DP), allocatable :: matrix_A(:,:), matrix_B(:,:), U(:),V(:)
    real (DP) :: omega(m)
    integer :: j

    !!allocate matrix and vector
    allocate(matrix_A(m,m),matrix_B(m,m))
    allocate(U(m),V(m))

    !! 'Finish creating matrix'
    call create_matrix_diffusionCN(MESAfile,m,delta_t,nu_diff,matrix_A,matrix_B)
    !call create_matrix_diffusionCN_mod(MESAfile,m,delta_t,nu_diff,matrix_A,matrix_B)
    !call printMatrix(matrix_A, m, m,1000)
   
    V = MESAfile(4,:) !omega
    do j=1,iter    !!subtimesteps
      !!initial conditions and boundary conditions
      !V(1) = V(2)
      !V(m) = V(m-1)

      !!'Finish multiplying matrix'
      U = 0.
      U = matmul(transpose(matrix_B),V)
      !print*, U(1),U(m),V(1),V(m)

      !!Tridiagonal matrix
      omega = 0.
      !call printMatrix(U, 1, m,2000)
      call invert_tridagmatrix(matrix_A,m,omega,U)
      !print*, 'Finish inverting matrix'

      V = omega
    end do

    deallocate(matrix_A,matrix_B,U,V)

  end function CrankNicolson
  

end module numerical_scheme
