Module rotationprofile_mod
  use parameters
  use data_handling
  use math_functions
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: get_modelN,rot_profile_conserveAM,rotation_profile_implicit,rotation_profile_explicit!,integration_coeff
  public :: iniprofile, rotation_profile, insertcol_MESAfile,get_r_and_w,implicit_test,rotation_profile_Henyey
  public :: rotation_profile_CN_radiative
  contains

  function get_modelN(MESAhistory) result(model)
    implicit none
    real(DP), intent(in) :: MESAhistory(:,:)
    real(DP) ::logL,bumpL,d
    integer :: m, i,model,start
 
    start = 500  !avoid MS where luminosity decreases
    m = size(MESAhistory,dim=2)
    model = start
    bumpL=0.

    logL = MESAhistory(38,start-1)
    do i=start,m
      if(logL> MESAhistory(38,i)) then
        bumpL = MESAhistory(38,i)
        EXIT
      end if
      logL = MESAhistory(38,i)
    end do

    d=100.
    do i=start,m
      if (abs(MESAhistory(38,i) - 0.85*bumpL) <d) then
        model = i
        d = abs(abs(MESAhistory(38,i) - 0.85*bumpL))
      end if
    end do

    print*,'Model = ',model


  end function get_modelN

  !find the radius and width of the burning shell
  subroutine get_r_and_w(MESAfile,r_t,w)
    implicit none
    real(DP), intent(in) :: MESAfile(:,:)
    real(DP), intent(inout) :: r_t,w
    real(DP) :: max,r1,r2,d
    integer :: i,m,j

    m = size(MESAfile,dim=2)

    r_t = 0.
    w=0.
    r2=0.
    r1=0.
    max = 0.
    j=1
    do i=1,m
      !! for RC star
      !if(MESAfile(1,i)/R_star > 0.003) then
        if(MESAfile(15,i) > max) then
          max = MESAfile(15,i)
          r_t = MESAfile(1,i)
          j = i
        end if
      !end if
    end do
    r_t = r_t/R_star
    
    max = max/3.
    d=10000.
    do i=j,1,-1
      if (abs(MESAfile(15,i) - max) <d) then
        r1 = MESAfile(1,i)
        d = abs(MESAfile(15,i) - max)
      end if
    end do
    d=10000.
    do i=j,m
      if (abs(MESAfile(15,i) - max) <d) then
        r2 = MESAfile(1,i)
        d = abs(MESAfile(15,i) - max)
      end if
    end do

    w=(r2-r1)/R_star

  end subroutine get_r_and_w

  function iniprofile(MESAfile,m,model) result(rot)
    implicit none
    integer, intent(in) :: m,model
    real(DP) :: Omega_s,Omega_c,r_t,w, rot(m), r(m)
    real(DP), intent(in) :: MESAfile(:,:)
    !end of MS constant profile
    !omega = 400.e-9 * 2*PI

    !!Belkacem et al. (2015b) table 2. values

    call get_r_and_w(MESAfile,r_t,w)

    Omega_s = 0.
    Omega_c = 0.


    if (model < 400) then
      !!SG model
      Omega_s = 200*2.*PI*1e-9
      Omega_c = 600.*2.*PI*1e-9
    else if (model > 400 .AND. model < 500) then
      !!Early RG model
      Omega_s = 55.*2.*PI*1e-9
      Omega_c = 1000.*2.*PI*1e-9
    else if (model > 500) then
      !!RG model
      Omega_s = 8.1*2.*PI*1e-9
      Omega_c = 500.*2.*PI*1e-9
    end if

    !!RC model  !Deheuvels et al. (2015) Table 6
    !Omega_s = 70.*2.*PI*1e-9
    !Omega_c = 160.*2.*PI*1e-9


    r = MESAfile(1,:)/R_star
    rot = Omega_s + (Omega_c-Omega_s)*(1.+ERF((r_t-r)/w))/2.
    
  end function iniprofile


  subroutine insertcol_MESAfile(model,number,omega)
    implicit none
    real (DP), allocatable :: filearray(:,:)
    real(DP), intent(in) :: omega(:)
    integer, allocatable :: a1(:)
    real(DP) :: array(5)
    integer, intent(in) :: model
    integer :: N,i,number
  
    N = getnlines('MESA/profile'//trim(string(model))//'.data.GYRE')
    N=N+5
    allocate(filearray(18,N),a1(N))
    
    open(200, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(200,*) array
      do i=1,N
        read(200,*) a1(i), filearray(:,i)
      end do
    close(200)

    open(205, action='write',file = 'MESA/profile.data.GYRE') !changed!!!
    write(205,*) int(array(1)), array(2), array(3), array(4), int(array(5))
      do i=1,N
        write(205,*) a1(i), filearray(:,i)
      end do
    close(205)
  
    !insert rot profile in new Mesa profile
    if (number == 1) then
      filearray(18,:) = iniprofile(filearray,N,model)
    else if (number ==3) then
      filearray(18,:) = omega(:)
    end if
  
    if (number == 1 .OR. number == 3) then
      !open(201, action='write',file = 'MESA/profile.data.GYRE') !changed!!!
      open(201, action='write',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      write(201,*) int(array(1)), array(2), array(3), array(4), int(array(5))
      do i=1,N
        write(201,*) a1(i), filearray(:,i)
      end do
      close(201)
    end if


    !open(202, action='write',file = 'MESA/new_rotprofile/profile'//trim(string(model))//'.data.GYRE')
    !  write(202,*) int(array(1)), array(2), array(3), array(4), int(array(5))
    !  do i=1,N
     !   write(202,*) a1(i), filearray(:,i)
     ! end do
    !close(202)
  
    
    deallocate(filearray,a1)
  
  end subroutine insertcol_MESAfile
 

  !!Backward Euler
  subroutine rotation_profile_implicit(F_total,MESAold,omega,model,m)
    implicit none
    real(DP) :: r_n_i1,r_n_i,omega_n_i(m),r_dot,Jacobian(m),f(m),rho_n_i,Tm(m),Tc(m),AMi,AMf
    real (DP), intent(inout) :: omega(:), F_total(:,:),MESAold(:,:)
    real(DP), allocatable :: new_F_total(:,:),MESAfile(:,:),MESAold_int(:,:)
    integer :: nsize,j,ri,rf,mold,i
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

      dx = matmul(matrix,vector_x0)

      x_new = x_new + dx

      x0 = x_new(1)
      y0 = x_new(2)
    end do


  end subroutine implicit_test

  subroutine rotation_profile_explicit(MESAold,omega,model,m)
    implicit none
    real (DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,r_dot,J_dot,dev_r,dev_t,Tm,Tc,A
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

  subroutine rotation_profile(F_total,MESAold,omega,model,m)
    implicit none
    real(DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,r_dot,J_dot,dev_r,dev_t,Tm,Tc
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
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),m)
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),m)
    MESAold_int(1,:) = MESAfile(3,:)
  
    !compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')

    do j=1,m
      r_n_i1 = MESAfile(2,j)
      r_n_i = MESAold_int(2,j)
      rho_n_i = MESAold_int(3,j)
      omega_n_i = MESAold_int(4,j)

      r_dot = (r_n_i1 - r_n_i)/dt
      dev_r = 2.*r_dot*omega_n_i*dt/(r_n_i)
      dev_t = omega_n_i
      J_dot = new_F_total(2,j)*dt/(rho_n_i*r_n_i**2)
      
      if (j<ri .OR. j>rf) then
        omega(j) = MESAold_int(4,j)
      else
        omega(j) = J_dot - dev_r + dev_t
      end if

      !compute timescales
      Tm = abs(rho_n_i*r_n_i**2 *omega_n_i /new_F_total(2,j))
      Tc = abs(-r_n_i/r_dot)

      write(203,*) r_n_i1/MESAfile(2,rf),abs(dev_t),abs(J_dot),&
      &abs(dev_r),omega(j), Tm,Tc,dt,&
      &r_n_i,r_dot,MESAold_int(1,j)
    end do


    close(203)
    
    deallocate(new_F_total,MESAfile,MESAold_int)

  end subroutine rotation_profile

  subroutine rot_profile_conserveAM(MESAold,omega,model,m)
    implicit none
    real(DP) :: r_n_i1,r_n_i,omega_n_i
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real(DP), allocatable :: MESAfile(:,:),MESAold_int(:,:)
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
    do j=1,m
      r_n_i1 = MESAfile(2,j)
      r_n_i = MESAold_int(2,j)
      omega_n_i = MESAold_int(4,j)
      
      if (j<ri .OR. j>rf) then
        omega(j) = MESAold_int(4,j)
      else
        omega(j) = r_n_i1*omega_n_i/(3*r_n_i1-2*r_n_i)
      end if

      write(203,*) r_n_i1/MESAfile(2,rf),MESAold_int(1,j)/M_star,&
      &0.,0., omega(j)

    end do
    close(203)
    
    deallocate(MESAfile,MESAold_int)

  end subroutine rot_profile_conserveAM


  !function integration_coeff(MESAfile,F_total,k) result(C) !array with 33 entries
    !implicit none
    !real (DP), intent(in) :: F_total(:,:),MESAfile(:,:)
    !integer, intent(in) :: k
    !real (DP) :: r,r_k,r_k_1,r_t,omega_t,rho,rho_k,rho_k_1,nu,deltanu,g_small,g_k,c_p,T,p
    !real (DP) :: C(33)

    !C(1) = r**2/dt
    !C(2) = r_t**2 * omega_t /dt
    !C(3) = (8.*PI/15.)*(R_sun**2 /M_sun)*(rho_k*r_k**4)/(sqrt(nu)*deltanu)
    !C(4) = (8.*PI/15.)*(R_sun**2 /M_sun)*(rho_k_1*r_k_1**4)/(sqrt(nu)*deltanu)
    !C(5) = 0.
    !C(6) = 0.
    !C(7) = (M_sun/L_sun)*(nu**(3./2.) *g_small *c_p*rho*T)/(p*L_star)
  !end function


  !!Crank-Nicolson
  subroutine rotation_profile_CN(MESAold,omega,model,m)
    implicit none
    real (DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,r_dot,J_dot,dev_r,dev_t,Tm,Tc,A
    real (DP) :: AMi,AMf,nu,delta_r,delta_rho, radius(m)
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:)
    real (DP), allocatable :: matrix_A(:,:), matrix_B(:,:), matrix_Ainv(:,:), U(:),V(:)
    real (DP), allocatable :: gamma(:), alpha(:), beta(:)
    integer :: j,ri,rf,mold,i
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

    !!create uniform grid in radius
    !delta_r = (MESAold_int(2,m) - MESAold_int(2,1))/m
    !radius(1) = MESAold_int(2,1)
    !do i=2,m
    !  radius(i) = radius(i-1) + delta_r
    !end do
  
    !compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)

    !! compute gamma, alpha and beta coefficients
    allocate(gamma(m),alpha(m),beta(m),U(m),V(m))
    nu = 1.e12!1. !!Difusion coefficient
    do i=2,m
      if (i==m) then
        !delta_r = MESAold_int(2,i) - MESAold_int(2,i-1)
        delta_rho = (MESAold_int(3,i) - MESAold_int(3,i-1))/delta_r
      else
        !delta_r = MESAold_int(2,i+1) - MESAold_int(2,i)
        delta_rho = (MESAold_int(3,i+1) - MESAold_int(3,i))/delta_r
      end if
      r_dot = MESAfile(2,i) - MESAold_int(2,i)
      
      gamma(i) = nu*dt/(2.*(delta_r**2))
      alpha(i) = 0.!(nu*dt/(4.*delta_r))*(delta_rho/MESAold_int(3,i) + 4./MESAold_int(2,i))
      beta(i)  = 0.!r_dot/MESAold_int(2,i)

    end do

    !!initial conditions and boundary conditions
    V = MESAold_int(4,:)
    V(1) = V(2)
    V(m) = V(m-1)

    !!allocate matrix and vector
    allocate(matrix_A(m,m),matrix_B(m,m),matrix_Ainv(m,m))

    !!fill matrix_A
    matrix_A = 0.
    do i=2,m-1
      matrix_A(i,i)   = 1. + 2.*gamma(i) + beta(i)
      matrix_A(i-1,i) = -(gamma(i) - alpha(i))
      matrix_A(i+1,i) = -(gamma(i) + alpha(i))
    end do
    !!Boundary conditions
    matrix_A(2,1)   = -(gamma(2) + alpha(2))
    matrix_A(1,1)   = -matrix_A(2,1)
    matrix_A(m-1,m) = -(gamma(m-1) - alpha(m-1))
    matrix_A(m,m)   = 1. + gamma(m) + beta(m) - alpha(m)

    !!fill matrix_B
    matrix_B = 0.
    do i=2,m-1
      matrix_B(i,i)   = 1. - 2.*gamma(i) - beta(i)
      matrix_B(i-1,i) = gamma(i) - alpha(i)
      matrix_B(i+1,i) = gamma(i) + alpha(i)
    end do
    !!Boundary conditions
    matrix_B(2,1)   = gamma(2) + alpha(2)
    matrix_B(1,1)   = -matrix_B(2,1)
    matrix_B(m-1,m) = gamma(m-1) - alpha(m-1)
    matrix_B(m,m)   = 1. - gamma(m) - beta(m) + alpha(m)
    !print*,matrix_A(1,1) , matrix_A(2,1), matrix_B(1,1),matrix_B(2,1)

    print*, 'Finish creating matrix'
    call printMatrix(matrix_A, m, m,1000)

    U = matmul(matrix_B,V)
    print*, 'Finish multiplying matrix'

    !! LU decompostion
    !V=U
    !call inverse_matrix(matrix_A,m,m,matrix_Ainv)
    !print*, 'Finish inverting matrix'
    !U = matmul(matrix_Ainv,V)
    !print*, 'Finish numerical method'
    !omega = U

    !!Tridiagonal matrix
    omega = 0.
    !call printMatrix(U, 1, m,2000)
    call invert_tridagmatrix(matrix_A,m,omega,U)
    print*, 'Finish inverting matrix'
  

    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    AMi = 0.
    AMf = 0.
    do j=1,m-1
      write(203,*) radius(j),0.,0.,0.,omega(j),MESAold_int(4,j),MESAold_int(3,j),&
      &(MESAfile(2,j) - MESAold_int(2,j))/dt,MESAold_int(2,j+1) - MESAold_int(2,j),&
      &MESAold_int(1,j+1) - MESAold_int(1,j)  
      AMi = AMi + MESAold_int(2,j)**2 * V(j)!omega_n_i
      AMf = AMf + MESAfile(2,j)**2 *omega(j)
    end do
    close(203)

    print*, 'delta AM radiative zones = ', AMf - AMi
    
    deallocate(matrix_A,matrix_B,matrix_Ainv,U,V,gamma,alpha,beta)
    deallocate(MESAfile,MESAold_int)

  end subroutine rotation_profile_CN


  !!Crank-Nicolson
  subroutine rotation_profile_CN_radiative(MESAold,omega,model,m)
    implicit none
    real (DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,r_dot,J_dot,dev_r,dev_t,Tm,Tc,A
    real (DP) :: AMi,AMf,nu,delta_r,delta_rho,dm,rmin
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:),omega_old(:),radius_old(:)
    real (DP), allocatable :: matrix_A(:,:), matrix_B(:,:), matrix_Ainv(:,:), U(:),V(:)
    real (DP), allocatable :: gamma(:), alpha(:), beta(:), omega_rad(:), radius(:),rho_old(:),rho(:)
    integer :: j,ri,rf,mold,i,mrad, muni_r
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


    !! create uniform grid in radius only in radiative part - simple version
    !delta_r = (MESAold_int(2,m) - MESAold_int(2,1))/m
    !radius(1) = MESAold_int(2,1)
    !do i=2,m
    !  radius(i) = radius(i-1) + delta_r
    !end do
    !! --------------------------------------------------------------------

    !! compute radiative region radius -------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    mrad = rf - 100
    allocate(omega_old(mrad),radius_old(mrad),rho_old(mrad))
    do i=1,mrad
      omega_old(i) = MESAold_int(4,i)
      radius_old(i) = MESAold_int(2,i)
      rho_old(i) = MESAold_int(3,i)
    end do

    !! create uniform grid in radius only in radiative part
    rmin = 1.e20
    do i=1,mrad-1
      delta_r = MESAold_int(2,i+1) - MESAold_int(2,i)
      if (delta_r < rmin) then
        rmin = delta_r
      end if
    end do
    delta_r = 20.*rmin

    muni_r = int((MESAold_int(2,mrad) - MESAold_int(2,1))/delta_r)
    print*,'muni_r ',muni_r, 'deltar = ',delta_r
    allocate(radius(muni_r),omega_rad(muni_r),rho(muni_r))
    radius(1) = MESAold_int(2,1)
    do i=2,muni_r
      radius(i) = radius(i-1) + delta_r
    end do

    !interpolation of for uniform grid ---------------------------------------------------------
    omega_rad = interpolate(radius_old,omega_old,mrad,radius,muni_r)
    rho = interpolate(radius_old,rho_old,mrad,radius,muni_r)
    mrad = muni_r

    !! compute rotation profile --------------------------------------------------------------
    
    !! compute gamma, alpha and beta coefficients
    allocate(gamma(mrad),alpha(mrad),beta(mrad),U(mrad),V(mrad))
    !!allocate matrix and vector
    allocate(matrix_A(mrad,mrad),matrix_B(mrad,mrad),matrix_Ainv(mrad,mrad))
    !do i=1,mrad
    !  omega_rad(i) = MESAold_int(4,i)
    !end do

    dt = dt*1.e-2
    nu = 1.e5!1.e12!1. !!Difusion coefficient
    do i=1,mrad
      if (i==1) then
        alpha(i) = 0.!(nu*dt/(4.*delta_r))*(delta_rho+ 4./radius(i))
      else if (i==mrad) then
        !delta_r = MESAold_int(2,i) - MESAold_int(2,i-1)
        !dm = MESAold_int(1,i) - MESAold_int(1,i-1)
        !delta_rho = (MESAold_int(3,i) - MESAold_int(3,i-1))/delta_r
        delta_rho = (log(rho(i)) - log(rho(i-1)))/delta_r
        alpha(i) = (nu*dt/(4.*delta_r))*(delta_rho)!+ 4./radius(i))
      else
        !delta_r = MESAold_int(2,i+1) - MESAold_int(2,i)
        !dm = MESAold_int(1,i+1) - MESAold_int(1,i)
        !delta_rho = (MESAold_int(3,i+1) - MESAold_int(3,i))/delta_r
        delta_rho = (log(rho(i+1)) - log(rho(i)))/delta_r
        alpha(i) = (nu*dt/(4.*delta_r))*(delta_rho)!+ 4./radius(i))
      end if
      !r_dot = MESAfile(2,i) - MESAold_int(2,i)
      gamma(i) = nu*dt/(2.*(delta_r**2)) 
      
      !alpha(i) = (nu*dt/(4.*delta_r))*(delta_rho/MESAold_int(3,i) + 4./MESAold_int(2,i))
      beta(i)  = 0.!r_dot/MESAold_int(2,i)
      !alpha(i) = 0.
    end do

    !!fill matrix_A
    matrix_A = 0.
    do i=2,mrad-1
      matrix_A(i,i)   = 1. + 2.*gamma(i) + beta(i)
      matrix_A(i-1,i) = -(gamma(i) - alpha(i))
      matrix_A(i+1,i) = -(gamma(i) + alpha(i))
    end do
    !!Boundary conditions
    matrix_A(2,1)   = -(gamma(2) + alpha(2))
    matrix_A(1,1)   = 1. + gamma(1) + beta(1) + alpha(1)!-matrix_A(2,1)
    matrix_A(mrad-1,mrad) = -(gamma(mrad-1) - alpha(mrad-1))
    matrix_A(mrad,mrad)   = 1. + gamma(mrad) + beta(mrad) - alpha(mrad)

    !!fill matrix_B
    matrix_B = 0.
    do i=2,mrad-1
      matrix_B(i,i)   = 1. - 2.*gamma(i) - beta(i)
      matrix_B(i-1,i) = gamma(i) - alpha(i)
      matrix_B(i+1,i) = gamma(i) + alpha(i)
    end do
    !!Boundary conditions
    matrix_B(2,1)   = gamma(2) + alpha(2)
    matrix_B(1,1)   = 1. - gamma(1) - beta(1) - alpha(1)!-matrix_B(2,1)
    matrix_B(mrad-1,mrad) = gamma(mrad-1) - alpha(mrad-1)
    matrix_B(mrad,mrad)   = 1. - gamma(mrad) - beta(mrad) + alpha(mrad)
    print*,matrix_A(mrad,mrad) , matrix_A(mrad-1,mrad), matrix_A(1,1),matrix_A(2,1),&
    &gamma(1), gamma(2) , alpha(2), alpha(mrad-1), alpha(mrad)

    !print*, 'Finish creating matrix'
    !call printMatrix(matrix_A, mrad, mrad,1000)
   
    do j=1,100
    !!initial conditions and boundary conditions
    V = omega_rad
    V(1) = V(2)
    V(mrad) = V(mrad-1)

    U = 0.
    U = matmul(matrix_B,V)
    
    !print*, 'Finish multiplying matrix'

    !! LU decompostion
    !V=U
    !call inverse_matrix(matrix_A,m,m,matrix_Ainv)
    !print*, 'Finish inverting matrix'
    !U = matmul(matrix_Ainv,V)
    !print*, 'Finish numerical method'
    !omega = U

    !!Tridiagonal matrix
    omega_rad = 0.
    !call printMatrix(U, 1, m,2000)
    call invert_tridagmatrix(matrix_A,mrad,omega_rad,U)
    
    print*, 'Finish inverting matrix'


  end do

  !interpolation of to the original grid ---------------------------------------------------------
  omega_old = interpolate(radius,omega_rad,mrad,radius_old,rf - 100)
  print*,'here3'

  !convective envelope stays constant omega
  mrad = rf - 100
  do i=1,m
    if (i > mrad) then
      omega(i) = omega_old(mrad)
    else
      omega(i) = omega_old(i)
    end if
  end do


    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    AMi = 0.
    AMf = 0.
    do j=1,m-1
      write(203,*) MESAfile(2,j),MESAold_int(2,j),0.,0.,omega(j),MESAold_int(4,j),MESAold_int(3,j),&
      &(MESAfile(2,j) - MESAold_int(2,j))/dt,MESAold_int(2,j+1) - MESAold_int(2,j),&
      &MESAold_int(1,j+1) - MESAold_int(1,j)  
      AMi = AMi + MESAold_int(2,j)**2 * MESAold_int(4,j)
      AMf = AMf + MESAfile(2,j)**2 *omega(j)
    end do
    close(203)

    print*, 'delta AM radiative zones = ', AMf - AMi,' dt = ',dt/(365.*24.*60.*60.),'dr = ',delta_r
    
    deallocate(matrix_A,matrix_B,matrix_Ainv,U,V,gamma,alpha,beta)
    deallocate(MESAfile,MESAold_int,radius,omega_rad,omega_old,radius_old)

  end subroutine rotation_profile_CN_radiative


  !!Henyey Scheme
  subroutine rotation_profile_Henyey(MESAold,omega,model,m)
    implicit none
    !real (DP) :: r_n_i1,r_n_i,rho_n_i,omega_n_i,r_dot,J_dot,dev_r,dev_t,Tm,Tc,A
    !real (DP) :: AMi,AMf,nu,delta_r,delta_rho
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:)
    integer :: j,ri,rf,mold,i
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


    deallocate(MESAfile,MESAold_int)
  end subroutine rotation_profile_Henyey

end module rotationprofile_mod
