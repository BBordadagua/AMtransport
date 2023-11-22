Module numerical_scheme
  use parameters
  use data_handling
  use math_functions
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: rotation_profile_implicit,rotation_profile_explicit,implicit_test
  public :: create_matrix_diffusionCN,create_matrix_diffusionCN_mod,create_matrix_Henyey
  public :: CrankNicolson, create_vector_Henyey
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
      matrix_A(1,i) = -(lambda(i) - alpha(i))
      matrix_A(2,i) = 1. + 2.*lambda(i) + beta(i) !diagonal
      matrix_A(3,i) = -(lambda(i) + alpha(i))
    end do
    !!Boundary conditions
    matrix_A(1,1) = 0. !b
    matrix_A(2,1) = 1. + lambda(1) + beta(1) + alpha(1) !diagonal
    matrix_A(3,1) = -(lambda(1) + alpha(1)) !a
    
    matrix_A(1,m) = -(lambda(m) - alpha(m)) !b
    matrix_A(2,m) = 1. + lambda(m) + beta(m) - alpha(m) !diagonal
    matrix_A(3,m) = 0. !a

    !!fill matrix_B
    matrix_B = 0.
    do i=2,m-1
      matrix_B(1,i) = lambda(i) - alpha(i)
      matrix_B(2,i) = 1. - 2.*lambda(i) - beta(i)
      matrix_B(3,i) = lambda(i) + alpha(i)
    end do
    !!Boundary conditions
    matrix_B(1,1) = 0.
    matrix_B(2,1) = 1. - lambda(1) - beta(1) - alpha(1)
    matrix_B(3,1) = lambda(1) + alpha(1)
    
    matrix_B(1,m) = lambda(m) - alpha(m)
    matrix_B(2,m) = 1. - lambda(m) - beta(m) + alpha(m)
    matrix_B(3,m) = 0.
    !print*,matrix_A(m,m) , matrix_A(m-1,m), matrix_A(1,1),matrix_A(2,1)


    deallocate(lambda,alpha,beta)

  end subroutine create_matrix_diffusionCN


  subroutine create_matrix_diffusionCN_mod(MESAfile,m,delta_t,nu_diff,rdot,matrix_A,matrix_B)
    implicit none
    real (DP), intent(in) :: MESAfile(:,:), delta_t,nu_diff,rdot(:)
    real (DP), intent(inout) :: matrix_A(:,:),matrix_B(:,:)
    integer, intent(in) :: m 
    real (DP), allocatable :: lambda(:),beta(:),alpha(:),epsilon(:),gamma(:),zeta(:)
    real (DP), allocatable :: a1(:), a2(:)
    integer :: i
    real (DP) :: const

    const = 1d10!1d-10

    allocate(lambda(m),beta(m),alpha(m),epsilon(m),gamma(m),zeta(m),a1(m),a2(m))

    do i=1,m
      if (i == 1) then
        !! diffusion terms
        lambda(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i+1)*MESAfile(2,i+1)**4)&
          &/((MESAfile(2,i+1)-MESAfile(2,i))*(MESAfile(2,i+1)-MESAfile(2,i)))

        beta(i) = 0.

        !! mixed modes term jdot
        alpha(i) = (MESAfile(5,i)*delta_t)/(MESAfile(3,i)*MESAfile(2,i)**2)

        !! contraction and expasion of the star
        epsilon(i) = rdot(i)*delta_t

        !! toy model flux
        !gamma(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i+1)+MESAfile(2,i))**2)/(MESAfile(2,i+1)-MESAfile(2,i))
        !zeta(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i)+MESAfile(2,i))**2)/(MESAfile(2,i+1)-MESAfile(2,i))

        !! j dot mid point
        a1(i) = 4.*(MESAfile(5,i)*delta_t)/((MESAfile(3,i)+MESAfile(3,i))*&
        &(MESAfile(2,i)+MESAfile(2,i))**2)
        a2(i) = 4.*(MESAfile(5,i)*delta_t)/((MESAfile(3,i)+MESAfile(3,i))*&
        &(MESAfile(2,i)+MESAfile(2,i))**2)

        a1(i) = 0.!0.5*(MESAfile(5,i-1)*delta_t)/((MESAfile(3,i-1))*&
        !&(MESAfile(2,i-1))**2)
        a2(i) = 0.5*(MESAfile(5,i)*delta_t)/((MESAfile(3,i))*&
        &(MESAfile(2,i))**2)

      else if (i == m) then
        lambda(i) = 0.

        beta(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i-1)*MESAfile(2,i-1)**4)&
          &/((MESAfile(2,i)-MESAfile(2,i-1))*(MESAfile(2,i)-MESAfile(2,i-1)))

        alpha(i) = (MESAfile(5,i)*delta_t)/(MESAfile(3,i)*MESAfile(2,i)**2)

        epsilon(i) = rdot(i)*delta_t

        !! toy model flux
        !gamma(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i)+MESAfile(2,i))**2)/(MESAfile(2,i)-MESAfile(2,i-1))
        !zeta(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i)+MESAfile(2,i-1))**2)/(MESAfile(2,i)-MESAfile(2,i-1))

        !! j dot mid point
        a1(i) = 4.*(MESAfile(5,i-1)*delta_t)/((MESAfile(3,i)+MESAfile(3,i-1))*&
        &(MESAfile(2,i)+MESAfile(2,i-1))**2)
        a2(i) = 4.*(MESAfile(5,i)*delta_t)/((MESAfile(3,i)+MESAfile(3,i-1))*&
        &(MESAfile(2,i)+MESAfile(2,i-1))**2)

        a1(i) = 0.5*(MESAfile(5,i-1)*delta_t)/((MESAfile(3,i-1))*&
        &(MESAfile(2,i-1))**2)
        a2(i) = 0.5*(MESAfile(5,i)*delta_t)/((MESAfile(3,i))*&
        &(MESAfile(2,i))**2)

      else
        lambda(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i+1)*MESAfile(2,i+1)**4)&
          &/((MESAfile(2,i+1)-MESAfile(2,i))*(MESAfile(2,i+1)-MESAfile(2,i-1)))

        beta(i) = (nu_diff*delta_t/(2.*MESAfile(3,i)*MESAfile(2,i)**4))*&
          &(MESAfile(3,i)*MESAfile(2,i)**4 + MESAfile(3,i-1)*MESAfile(2,i-1)**4)&
          &/((MESAfile(2,i)-MESAfile(2,i-1))*(MESAfile(2,i+1)-MESAfile(2,i-1)))

        alpha(i) = (MESAfile(5,i)*delta_t)/(MESAfile(3,i)*MESAfile(2,i)**2)

        epsilon(i) = rdot(i)*delta_t

        !! toy model flux
        !gamma(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i+1)+MESAfile(2,i))**2)/(MESAfile(2,i+1)-MESAfile(2,i-1))
        !zeta(i) = const*(delta_t/(8.*MESAfile(3,i)*MESAfile(2,i)**4))*&
        !&((MESAfile(2,i)+MESAfile(2,i-1))**2)/(MESAfile(2,i+1)-MESAfile(2,i-1))

        !! j dot mid point
        a1(i) = 4.*(MESAfile(5,i-1)*delta_t)/((MESAfile(3,i)+MESAfile(3,i-1))*&
        &(MESAfile(2,i)+MESAfile(2,i-1))**2)
        a2(i) = 4.*(MESAfile(5,i)*delta_t)/((MESAfile(3,i)+MESAfile(3,i-1))*&
        &(MESAfile(2,i)+MESAfile(2,i-1))**2)

        a1(i) = 0.5*(MESAfile(5,i-1)*delta_t)/((MESAfile(3,i-1))*&
        &(MESAfile(2,i-1))**2)
        a2(i) = 0.5*(MESAfile(5,i)*delta_t)/((MESAfile(3,i))*&
        &(MESAfile(2,i))**2)

      end if
    end do

    !! diffusion term
    !lambda = 0.
    !beta = 0.
    !! contraction and expansion term
    !epsilon = 0.
    !! jdot term
    alpha = 0.
    !! toy model flux term
    gamma = 0.
    zeta = 0.
    !! toy model integrated flux term
    a1 = 0.
    a2 = 0.

    !!fill matrix_A
    matrix_A = 0.
    do i=2,m-1
      matrix_A(1,i)= -beta(i)+zeta(i)-a1(i)
      matrix_A(2,i)= 1.+lambda(i)+beta(i)+epsilon(i)-alpha(i)-gamma(i)+zeta(i)-a2(i)
      matrix_A(3,i)= -lambda(i)-gamma(i)
    end do
    !!Boundary conditions
    matrix_A(1,1) = 0.
    matrix_A(2,1) = 1.+lambda(1)+epsilon(1)-alpha(1)-gamma(1)+2.*zeta(1)-a2(1)-a1(1)
    matrix_A(3,1) = -lambda(1)-gamma(1)

    matrix_A(1,m) = -beta(m)+zeta(m)-a1(m)
    matrix_A(2,m) = 1.+beta(m)+epsilon(m)-alpha(m)-2.*gamma(m)+zeta(m)-a2(m)
    matrix_A(3,m) = 0.

    !!fill matrix_B
    matrix_B = 0.
    do i=2,m-1
      matrix_B(1,i) = beta(i)-zeta(i)+a1(i)
      matrix_B(2,i) = 1.-lambda(i)-beta(i)-epsilon(i)+alpha(i)+gamma(i)-zeta(i)+a2(i)
      matrix_B(3,i) = lambda(i)+gamma(i)
    end do
    !!Boundary conditions
    matrix_B(1,1) = 0.
    matrix_B(2,1) = 1.-lambda(1)-epsilon(1)+alpha(1)+gamma(1)-2.*zeta(1)+a2(1)+a1(1)
    matrix_B(3,1) = lambda(1) +gamma(1)

    matrix_B(1,m) = beta(m) -zeta(m)+a1(m)
    matrix_B(2,m) = 1.-beta(m)-epsilon(m)+alpha(m)+2.*gamma(m)-zeta(m)+a2(m)
    matrix_B(3,m) = 0.

    !call printMatrix(matrix_A, 3, m,1000)

    deallocate(lambda,beta,alpha,epsilon,gamma,zeta,a1,a2)

  end subroutine create_matrix_diffusionCN_mod


  function CrankNicolson(MESAfile,m,iter,delta_t,nu_diff,rdot) result(omega)
    implicit none
    real (DP), intent(in) :: delta_t, MESAfile(:,:),nu_diff,rdot(:)
    integer, intent(in) :: iter,m
    real (DP), allocatable :: U(:),V(:),mA_triag(:,:),mB_triag(:,:)
    real (DP) :: omega(m)
    integer :: j

    !!allocate matrix and vector
    allocate(mA_triag(3,m),mB_triag(3,m))
    allocate(U(m),V(m))

    !! Create matrix
    !call create_matrix_diffusionCN(MESAfile,m,delta_t,nu_diff,mA_triag,mB_triag)
    call create_matrix_diffusionCN_mod(MESAfile,m,delta_t,nu_diff,rdot,mA_triag,mB_triag)
   
    V = MESAfile(4,:) !omega
    do j=1,iter    !!subtimesteps

      !! multiplying matrix
      U = 0.
      !U = matmul(transpose(matrix_B),V)
      U = multiply_matrix_vector_triag(mB_triag,V,m)

      !!invert Tridiagonal matrix
      omega = 0.
      call invert_tridagmatrix(mA_triag,m,omega,U)

      V = omega
    end do

    deallocate(U,V,mA_triag,mB_triag)

  end function CrankNicolson



  !! ----------------------------Henyey Scheme ----------------------------------
  
  subroutine create_matrix_Henyey(MESAfile,oldMESAfile,m,delta_t,nu_diff,matrix)
    implicit none
    real (DP), intent(in) :: MESAfile(:,:), delta_t,nu_diff, oldMESAfile(:,:)
    real (DP), intent(inout) :: matrix(:,:)
    integer, intent(in) :: m 
    real (DP), allocatable :: S1(:,:),S2(:,:),C1(:),C2(:),C3(:),C4(:),C5(:),C6(:)
    real (DP) :: delta_nu,radius,rho,nu,radius_old,omega_old,g_surf!,r_dot
    integer :: i

    allocate(S1(4,m),S2(4,m),C1(m),C2(m),C3(m),C4(m),C5(m),C6(m))

    !MESAfile(1,:) = mass
    !MESAfile(2,:) = radius
    !MESAfile(3,:) = rho
    !MESAfile(4,:) = omega
    !MESAfile(5,:) = Jmodes
    !MESAfile(6,:) = g
    !MESAfile(7,:) = nu

    do i=1,m
      !! boundary condition (domega/dnu = 0 in the boundary)
      if (i==1) then 
        delta_nu   = MESAfile(7,i) !!! caution division by zero!!!
        radius     = MESAfile(2,i)
        rho        = MESAfile(3,i)
        g_surf     = MESAfile(6,i)
        nu         = MESAfile(7,i) 
        radius_old = oldMESAfile(2,i)
        omega_old  = oldMESAfile(4,i)
        C6(i) = (4.*PI*MESAfile(6,i)*MESAfile(3,i)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      else
        delta_nu   = MESAfile(7,i) - MESAfile(7,i-1)
        radius     = (MESAfile(2,i) + MESAfile(2,i-1))/2.
        rho        = (MESAfile(3,i) + MESAfile(3,i-1))/2.
        g_surf     = (MESAfile(6,i) + MESAfile(6,i-1))/2.
        nu         = (MESAfile(7,i) + MESAfile(7,i-1))/2.
        radius_old = (oldMESAfile(2,i) + oldMESAfile(2,i-1))/2.
        omega_old  = (oldMESAfile(4,i) + oldMESAfile(4,i-1))/2.

        C6(i) = (4.*PI*MESAfile(6,i-1)*MESAfile(3,i-1)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      end if

      C1(i) = radius**2 /delta_t
      C2(i) = radius_old**2 *omega_old /delta_t
      C3(i) = 0. !MESAfile(4,i)/(rho*R_sun**2)  !Jmodes
      C4(i) = (16.*PI*R_sun**4 *radius**4 *rho)/(9.*M_sun*g_surf*sqrt(nu)*delta_nu)
      C5(i) = (4.*PI*MESAfile(6,i)*MESAfile(3,i)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)

    end do

      S1(1,:) = C1/2.
      S1(2,:) = C6
      S1(3,:) = C1/2.
      S1(4,:) = -C5

      S2(1,:) = C4
      S2(2,:) = 0.5d0
      S2(3,:) = -C4
      S2(4,:) = 0.5d0

    !!fill matrix_A
    matrix = 0.
    do i=3,2*m-1,2
      matrix(1,i) = S2(1,i)
      matrix(2,i) = S2(2,i)
      matrix(3,i) = S2(3,i)
      matrix(4,i) = S2(4,i)
      matrix(5,i) = 0.
    end do
    do i=1,2*m-1,2
      matrix(i,i+1)   = S1(1,i)
      matrix(i+1,i+1) = S1(2,i)
      matrix(i+2,i+1) = S1(3,i)
      matrix(i+3,i+1) = S1(4,i)

      matrix(i,i+2)   = S2(1,i)
      matrix(i+1,i+2) = S2(2,i)
      matrix(i+2,i+2) = S2(3,i)
      matrix(i+3,i+2) = S2(4,i)
    end do

    !!Boundary conditions  
    !!!! missing mixed modes!!!
    matrix(1,1)   = 0. !b
    matrix(2,1)   = 0. !c
    matrix(3,1)   = 0. !diagonal
    matrix(4,1)   = 1. !e
    matrix(5,1)   = 0. !f

    matrix(1,2*m+2) = 0.
    matrix(2,2*m+2) = 0.
    matrix(3,2*m+2) = 1.  !diagonal
    matrix(4,2*m+2) = 0.
    matrix(5,2*m+2) = 0.


    deallocate(S1,S2,C1,C2,C3,C4,C5,C6)

  end subroutine create_matrix_Henyey

  subroutine create_vector_Henyey(MESAfile,oldMESAfile,m,delta_t,nu_diff,vector_E,omega)
    implicit none
    real (DP), intent(in) :: MESAfile(:,:), delta_t,nu_diff, oldMESAfile(:,:),omega(:)
    real (DP), intent(inout) :: vector_E(:)
    integer, intent(in) :: m 
    real (DP), allocatable :: C1(:),C2(:),C3(:),C4(:),C5(:),C6(:),y1(:),y1k(:),y2(:),y2k(:)
    real (DP) :: delta_nu,radius,rho,nu,radius_old,omega_old,g_surf!,r_dot
    integer :: i

    allocate(C1(m),C2(m),C3(m),C4(m),C5(m),C6(m),y1(m),y1k(m),y2(m),y2k(m))

    !MESAfile(1,:) = mass
    !MESAfile(2,:) = radius
    !MESAfile(3,:) = rho
    !MESAfile(4,:) = omega
    !MESAfile(5,:) = Jmodes
    !MESAfile(6,:) = g
    !MESAfile(7,:) = nu

    do i=1,m
      !! boundary condition (domega/dnu = 0 in the boundary)
      if (i==1) then 
        delta_nu   = MESAfile(7,i) !!! caution division by zero!!!
        radius     = MESAfile(2,i)
        rho        = MESAfile(3,i)
        g_surf     = MESAfile(6,i)
        nu         = MESAfile(7,i) 
        radius_old = oldMESAfile(2,i)
        omega_old  = oldMESAfile(4,i)
        C6(i) = (4.*PI*MESAfile(6,i)*MESAfile(3,i)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      else
        delta_nu   = MESAfile(7,i) - MESAfile(7,i-1)
        radius     = (MESAfile(2,i) + MESAfile(2,i-1))/2.
        rho        = (MESAfile(3,i) + MESAfile(3,i-1))/2.
        g_surf     = (MESAfile(6,i) + MESAfile(6,i-1))/2.
        nu         = (MESAfile(7,i) + MESAfile(7,i-1))/2.
        radius_old = (oldMESAfile(2,i) + oldMESAfile(2,i-1))/2.
        omega_old  = (oldMESAfile(4,i) + oldMESAfile(4,i-1))/2.

        C6(i) = (4.*PI*MESAfile(6,i-1)*MESAfile(3,i-1)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      end if

      C1(i) = radius**2 /delta_t
      C2(i) = radius_old**2 *omega_old /delta_t
      C3(i) = 0. !MESAfile(4,i)/(rho*R_sun**2)  !Jmodes
      C4(i) = (16.*PI*R_sun**4 *radius**4 *rho)/(9.*M_sun*g_surf*sqrt(nu)*delta_nu)
      C5(i) = (4.*PI*MESAfile(6,i)*MESAfile(3,i)*nu_diff)/(M_sun*sqrt(nu)*delta_nu)
      
    end do

    do i=1,m
      y1k(i) =  omega(i)
    end do
    do i = 1,m
      if(i == 1) then
        y1(i) =  y1k(i) !k-1/2
        y2k(i) = C4(i+1)*(y1k(i+1)-y1k(i))/2. !k-1/2 + k+1/2
        y2(i) = 0. !k-1/2
      else
        y1(i) =  (y1k(i) + y1k(i-1))/2. !k-1/2
        y2k(i) =  C4(i)*(y1k(i)-y1k(i-1))/2. + C4(i+1)*(y1k(i+1)-y1k(i))/2. !k-1/2 + k+1/2
        y2(i) = C4(i)*(y1k(i) - y1k(i-1))/2. !k-1/2
      end if
    end do

    !! fill vector E, size 2 boundary cond + 2N
    do i=2,2*m+1,2
      vector_E(i) = -(C1(i)*y1(i)-C2(i)-C5(i)*y2k(i)+C6(i)*y2k(i-1)-C3(i)) !!k-1/2
      vector_E(i+1) = -(y2(i)-C4(i)*(y1k(i) - y1k(i-1))) !!k-1/2
    end do
    !!Boundary conditions
    vector_E(1) = -y2(1)
    vector_E(2*m+2) = -y2(m)

    deallocate(C1,C2,C3,C4,C5,C6,y1,y1k,y2,y2k)

  end subroutine create_vector_Henyey

end module numerical_scheme