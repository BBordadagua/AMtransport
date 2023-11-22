Module rotationprofile_mod
  use parameters
  use data_handling
  use math_functions
  use numerical_scheme
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: insertcol_MESAfile, get_r_and_w, iniprofile,get_modelN
  public :: rotation_profile,rotation_profile_Henyey,rot_profile_conserveAM
  public :: rotation_profile_CN_radiative,rotation_profile_CN_nonunigrid
  public :: create_uniform_grid_size,create_uniform_grid
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
  subroutine rotation_profile_CN_nonunigrid(MESAold,omega,model,m,F_total)
    implicit none
    real (DP) :: AMi,AMf, delta_t,nu_diff,o1,o2 !, delta_r
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), intent(in) :: F_total(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:),MESA_unigrid(:,:)
    real (DP), allocatable ::  omega_rad(:), r_dot(:)!, radius(:),omega_old(:)
    integer :: j,mold,i, m_nonuni,iter,ri,rf !,mrad,
    integer, intent(in) :: model,m
    character(len=18) :: dummy

    !allocate arrays
    allocate(MESAfile(19,m))
    open(200, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(200,*) dummy
      read(200,*) MESAfile
    close(200)

    mold = size(MESAold,dim=2) 
    allocate(MESAold_int(5,m))

    
    !! interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),m)
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),m)
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),m)
    MESAold_int(1,:) = MESAfile(3,:)
    MESAold_int(5,:) = interpolate(F_total(3,:),F_total(2,:),size(F_total,dim=2),MESAfile(3,:),m)

    
    !! take r=0 from the arrays
    call getradiativeR(ri,rf,MESAfile)
    m_nonuni = m-1!rf-100 !! how many points before convective boundary
    allocate(MESA_unigrid(5,m_nonuni))
    do i=1,m_nonuni !skiping r=0
      MESA_unigrid(1,i) = MESAold_int(1,i+1) !mass
      MESA_unigrid(2,i) = MESAold_int(2,i+1) !radius
      MESA_unigrid(3,i) = MESAold_int(3,i+1) !rho
      MESA_unigrid(4,i) = MESAold_int(4,i+1) !omega
      if (i>rf) then
        MESA_unigrid(5,i) = 0.
      else
        MESA_unigrid(5,i) = MESAold_int(5,i+1) !jdot
      end if
    end do

    open(300, action='write',file='output/total_flux_'//trim(string(model))//'.txt')
    call printMatrix(MESA_unigrid, 5, m_nonuni,300)
    close(300) 

    !! compute rdot
    allocate(r_dot(m_nonuni))
    do i=1,m_nonuni
      r_dot(i) = (MESAfile(2,i+1) - MESA_unigrid(2,i) )/(MESA_unigrid(2,i)*dt)
    end do
    
  !! compute rotation profile --------------------------------------------------------------
  nu_diff = 1d5!1d7 !1.1d3 !min for stability
  iter = 1000
  delta_t = 0.001*dt
  omega_rad = CrankNicolson(MESA_unigrid,m_nonuni,iter,delta_t,nu_diff,r_dot)

  !convective envelope stays constant omega
  omega(1) = omega_rad(1)
  do i=2,m
    if (i > m_nonuni) then
      omega(i) = omega_rad(m_nonuni)
    else
      omega(i) = omega_rad(i-1)
    end if
  end do


    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    
    do j=1,m-1
      write(203,*) MESAfile(2,j),MESAold_int(2,j),0.,0.,omega(j),MESAold_int(4,j),MESAold_int(3,j),&
      &(MESAfile(2,j) - MESAold_int(2,j))/dt,MESAold_int(2,j+1) - MESAold_int(2,j),&
      &MESAold_int(1,j+1) - MESAold_int(1,j)  
    end do
    close(203)

    !! checking if AM is conserved
    AMi = 0.
    AMf = 0.
    o1 = 0.
    o2 = 0.
    do j=2,m_nonuni
      AMi = AMi + MESAold_int(2,j)**2 * MESAold_int(4,j)
      AMf = AMf + MESAfile(2,j)**2 *omega(j)
      o1 = o1 + MESAold_int(4,j)
      o2 = o2 + omega(j)
      !write(203,*) MESA_unigrid(2,j),0,0,0,omega_rad(j)
    end do

    print*, 'delta AM  = ', (AMf - AMi)/AMi, o1-o2
    print*,'AMi = ',AMi
    print*,'AMf = ',AMf
    print*, ' dt[s] = ',delta_t,' dt[yr] = ',delta_t/(365.*24.*60.*60.)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    
    deallocate(MESAfile,MESAold_int,MESA_unigrid,omega_rad,r_dot)

  end subroutine rotation_profile_CN_nonunigrid


  !!Crank-Nicolson
  subroutine rotation_profile_CN_radiative(MESAold,omega,model,m)
    implicit none
    real (DP) :: AMi,AMf, delta_r, delta_t,nu_diff, o1, o2
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:),MESA_unigrid(:,:)
    real (DP), allocatable ::  omega_rad(:), radius(:),omega_old(:),rdot(:)
    integer :: j,ri,rf,mold,i,mrad, m_nonuni,iter
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


    !! create uniform grid
    call getradiativeR(ri,rf,MESAfile)
    m_nonuni = m-1!rf!-100 !! how many points before convective boundary
    mrad = create_uniform_grid_size(MESAold_int,MESAfile,m_nonuni)
    allocate(MESA_unigrid(4,mrad))
    call create_uniform_grid(MESAold_int,MESAfile,MESA_unigrid,m_nonuni,delta_r)

    !! compute rdot
    allocate(rdot(m_nonuni))
    do i=1,m_nonuni
      rdot(i) = (MESAfile(2,i+1) - MESA_unigrid(2,i) )/MESA_unigrid(2,i)
    end do
    

  !! compute rotation profile --------------------------------------------------------------
  nu_diff = 1d7 !1.1d3 !min for stability
  iter = 100
  delta_t = 0.01*dt
  omega_rad = CrankNicolson(MESA_unigrid,mrad,iter,delta_t,nu_diff,rdot)

  !interpolation of to the original non uniform grid ---------------------------------------------------------
  allocate(radius(m_nonuni),omega_old(m_nonuni))
  do i=1,m_nonuni !skiping r=0
    radius(i) = MESAold_int(2,i+1)
  end do
  omega_old = interpolate(MESA_unigrid(2,:),omega_rad,mrad,radius,m_nonuni)
    !! r=0 is not being accounted !!!!!!!!

  !convective envelope stays constant omega
  omega(1) = omega_old(1)
  do i=2,m
    if (i > m_nonuni) then
      omega(i) = omega_old(m_nonuni)
    else
      omega(i) = omega_old(i-1)
    end if
  end do


    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    
    do j=1,m-1
      write(203,*) MESAfile(2,j),MESAold_int(2,j),0.,0.,omega(j),MESAold_int(4,j),MESAold_int(3,j),&
      &(MESAfile(2,j) - MESAold_int(2,j))/dt,MESAold_int(2,j+1) - MESAold_int(2,j),&
      &MESAold_int(1,j+1) - MESAold_int(1,j)  
    end do
    close(203)

    !! checking if AM is conserved
    AMi = 0.
    AMf = 0.
    o1 = 0.
    o2 = 0.
    do j=1,m_nonuni
      AMi = AMi + MESAold_int(2,j)**2 * MESAold_int(4,j)
      AMf = AMf + MESAold_int(2,j)**2 *omega(j)
      o1 = o1 + MESAold_int(4,j)
      o2 = o2 + omega(j)
      !write(203,*) MESA_unigrid(2,j),0,0,0,omega_rad(j)
    end do
    !print*, 'm_nonuniform grid = ',m_nonuni

    print*, 'delta AM radiative zones = ', (AMf - AMi)/AMi, o2-o1
    print*, 'AMinitial = ', AMi, ' AMfin = ', AMf 
    print*, ' dt[s] = ',delta_t,' dt[yr] = ',delta_t/(365.*24.*60.*60.)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    
    deallocate(omega_old,radius,omega_rad)
    deallocate(MESAfile,MESAold_int,MESA_unigrid,rdot)

  end subroutine rotation_profile_CN_radiative



  !! Henyey Scheme -----------------------------------------------------------------------------------
  subroutine rotation_profile_Henyey(MESAold,omega,model,mmesa,F_total)
    implicit none
    real (DP) :: delta_t,nu_diff
    real (DP), intent(in) :: F_total(:,:)
    real (DP), intent(inout) :: omega(:),MESAold(:,:)
    real (DP), allocatable :: MESAfile(:,:),MESAold_int(:,:)
    integer :: ri,rf,mold,j
    integer, intent(in) :: model,mmesa
    character(len=18) :: dummy
    !! relaxation scheme variables
    integer :: ne,m,nb,nci,ncj,nck,nsi,nsj,nyj,nyk,itmax
    integer, allocatable :: indexv(:)
    real (DP) :: conv, slowc,AMi,AMf
    real (DP), allocatable :: c(:,:,:),s(:,:),scalv(:),y(:,:)

    !allocate arrays
    mold = size(MESAold,dim=2)
    allocate(MESAfile(19,mmesa),MESAold_int(7,mmesa))
    open(200, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(200,*) dummy
      read(200,*) MESAfile
    close(200)

    !interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),mmesa) !radius
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),mmesa) !rho
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),mmesa) !omega
    MESAold_int(1,:) = MESAfile(3,:) !mass
    MESAold_int(5,:) = interpolate(F_total(3,:),F_total(2,:),size(F_total,dim=2),MESAfile(3,:),mmesa)
    MESAold_int(5,1) = 0. !! excluding the infinit point
    MESAold_int(6,:) = 1.!g
    MESAold_int(7,:) = (MESAfile(3,:)/M_sun)**(2./3.) !nu

    !compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    omega = 0.
    
    !! relaxation method ------------------------------------------------------------------
    m=rf !! total of mesh points
    ne = 2 !! The problem involves ne equations for ne adjustable dependent variables at each point
    nyj = ne !! number of dependent variables
    nyk = m !! total of mesh points

    !! The arrays c(1:nci,1:ncj,1:nck), s(1:nsi,1:nsj) supply dummy storage used by the relaxation code
    !! the minimum dimensions must satisfy: nci=ne, ncj=ne-nb+1, nck=m+1, nsi=ne, nsj=2*ne+1.
    nb=1  !! The nb boundary conditions at the first mesh point
    nci=ne
    ncj=ne-nb+1
    nck=m+1
    nsi=ne
    nsj=2*ne+1

    !! allocate arrays
    allocate(indexv(ne),c(nci,ncj,nck),s(nsi,nsj),scalv(ne),y(ne,m))
    
    itmax=1000!100   !! maximum number of iterations
    conv=5.e-6  !! the convergence criterion
    slowc=1.    !! controls the fraction of corrections actually used after each iteration.
    nu_diff = 1d3!1d5 !! difusion coefficient
    delta_t = 0.001*dt !0.01*dt !! timestep
    
    !! indexv(j) describes which column of s(i,j) the variable y(j) has been put in
    !! The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv
    indexv(1)=2 
    indexv(2)=1
    
    !! scalv(1:nyj) contains typical sizes for each dependent variable, used to weight errors.
    scalv(1) = 1.     !abs(anorm) !!!
    scalv(2) = 1.     !max(abs(anorm),y(2,M)) !!!change!!!

    !! variables - initial guesses
    do j=1,m
     y(1,j)=MESAold_int(4,j)
     y(2,j)=(2./sqrt(PI))*exp(-MESAold_int(2,j)**2)
    end do
    y(2,1)=0d0
    y(2,m)=0d0


    call solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m,y,c,s,&
    &MESAfile,MESAold_int,delta_t,nu_diff)

    do j=1,m
      omega(j) = y(1,j)
    end do
    

    !! print final rotation profile for this model
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    do j=1,m-1
      write(203,*) MESAfile(2,j),MESAold_int(2,j),0.,0.,omega(j),MESAold_int(4,j) 
    end do
    close(203)

    !! checking if AM is conserved
    AMi = 0.
    AMf = 0.
    do j=1,m
      AMi = AMi + MESAold_int(2,j)**2 * MESAold_int(4,j)
      AMf = AMf + MESAold_int(2,j)**2 *omega(j)
    end do
  
    write(*,555) 'delta AM rad zones', (AMf - AMi)/AMi
    write(*,556) ' AMini', AMi, '   AMfin', AMf 
    write(*,556) ' dt[s]',delta_t,'  dt[yr]',delta_t/(365.*24.*60.*60.)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    555 format(A19, ES14.3)
    556 format(A6, ES14.3,A8,ES14.3)

    deallocate(MESAfile,MESAold_int,c,s,scalv,indexv,y)
  end subroutine rotation_profile_Henyey



  subroutine create_uniform_grid(MESAold_int,MESAfile,MESA_unigrid,mrad,delta_r)
    implicit none
    integer :: ri,rf,i,muni_r
    integer, intent(in) :: mrad
    real (DP), intent(in) :: MESAold_int(:,:),MESAfile(:,:)
    real (DP), intent(inout) :: MESA_unigrid(:,:), delta_r
    real (DP), allocatable :: omega_old(:),radius_old(:),rho_old(:),mass_old(:),radius(:)
    real (DP) :: rmin


    !! compute radiative region radius -------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    !mrad = rf - 100
    allocate(omega_old(mrad),radius_old(mrad),rho_old(mrad),mass_old(mrad))
    do i=1,mrad !skiping r=0
      mass_old(i)   = MESAold_int(1,i+1)
      radius_old(i) = MESAold_int(2,i+1)
      rho_old(i)    = MESAold_int(3,i+1)
      omega_old(i)  = MESAold_int(4,i+1)
    end do

    !! create uniform grid in radius only in radiative part (ignoring r=0)
    rmin = 1.e20
    do i=2,mrad
      delta_r = MESAold_int(2,i+1) - MESAold_int(2,i)
      if (delta_r < rmin) then
        rmin = delta_r
      end if
    end do
    delta_r = rmin

    muni_r = int((MESAold_int(2,mrad+1) - MESAold_int(2,2))/delta_r)
    if (muni_r > 8000) then
      muni_r = 8000
      delta_r = int((MESAold_int(2,mrad+1) - MESAold_int(2,2))/muni_r)
    end if
    print*,'muni_r ',muni_r, 'deltar = ',delta_r
    !allocate(radius(muni_r),omega_rad(muni_r),rho(muni_r))

    !! create new equally spaced radius
    allocate(radius(muni_r))
    radius(1) = MESAold_int(2,2) !ignoring r=0
    do i=2,muni_r
      radius(i) = radius(i-1) + delta_r
    end do

    !interpolation of for uniform grid ---------------------------------------------------------
    MESA_unigrid(2,:) = radius
    MESA_unigrid(3,:) = interpolate(radius_old,rho_old,mrad,radius,muni_r) !rho
    MESA_unigrid(4,:) = interpolate(radius_old,omega_old,mrad,radius,muni_r) !omega
    MESA_unigrid(1,:) = interpolate(radius_old,mass_old,mrad,radius,muni_r) !mass

    deallocate(radius_old,radius,omega_old,rho_old,mass_old)

  end subroutine create_uniform_grid

  function create_uniform_grid_size(MESAold_int,MESAfile,mrad) result(muni_r)
    implicit none
    integer :: ri,rf,i, muni_r
    integer, intent(in) :: mrad
    real (DP), intent(in) :: MESAold_int(:,:),MESAfile(:,:)
    real (DP) :: delta_r,rmin

    !! compute radiative region radius -------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    !mrad = rf - 100
    !! create uniform grid in radius only in radiative part (ignoring r=0)
    rmin = 1.e20
    do i=2,mrad
      delta_r = MESAold_int(2,i+1) - MESAold_int(2,i)
      if (delta_r < rmin) then
        rmin = delta_r
      end if
    end do
    delta_r = rmin

    muni_r = int((MESAold_int(2,mrad+1) - MESAold_int(2,2))/delta_r)
    if (muni_r > 8000) then
      muni_r = 8000
      delta_r = int((MESAold_int(2,mrad+1) - MESAold_int(2,2))/muni_r)
    end if

  end function create_uniform_grid_size





  

end module rotationprofile_mod
