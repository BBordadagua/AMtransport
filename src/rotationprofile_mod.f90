Module rotationprofile_mod
  use parameters
  use data_handling
  use math_functions
  use numerical_scheme
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: insertcol_MESAfile, get_r_and_w, iniprofile,get_modelN
  public :: rotation_profile,rotation_profile_Henyey,rot_profile_conserveAM
  public :: rotation_profile_CN
  public :: getG
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

    Omega_s = 0d0
    Omega_c = 0d0


    if (model < 400) then
      !!SG model
      Omega_s = 200d0*2d0*PI*1.e-9
      Omega_c = 600d0*2d0*PI*1.e-9
    else if (model > 400 .AND. model < 500) then
      !!Early RG model
      Omega_s = 55d0*2d0*PI*1.e-9
      Omega_c = 1000d0*2d0*PI*1.e-9
    else if (model > 500) then
      !!RG model
      Omega_s = 8.1d0*2d0*PI*1.e-9
      Omega_c = 500d0*2d0*PI*1.e-9
    end if

    !!RC model  !Deheuvels et al. (2015) Table 6
    !Omega_s = 70.*2.*PI*1e-9
    !Omega_c = 160.*2.*PI*1e-9


    r = MESAfile(1,:)/R_star
    rot = Omega_s + (Omega_c-Omega_s)*(1d0+ERF((r_t-r)/w))/2d0
    
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


  !! ---------------------------------------------------------------------------
  !! ---------------------- Crank-Nicolson at next model non-uniform grid ---------------------------
  !! ---------------------------------------------------------------------------
 
  subroutine rotation_profile_CN(MESAold,MESAfile,omega,model,F_total)
    implicit none
    real (DP) :: AMi,AMf, delta_t,nu_diff,o1,o2 !, delta_r
    real (DP), intent(inout) :: omega(:)
    real (DP), intent(in) :: F_total(:,:),MESAfile(:,:),MESAold(:,:)
    real (DP), allocatable :: MESAold_int(:,:),Mo(:,:),Mn(:,:)
    real (DP), allocatable ::  omega_rad(:), r_dot(:)!, radius(:),omega_old(:)
    integer :: j,mold,i, mcn,iter,ri,rf,m !,mrad,
    integer, intent(in) :: model

    !allocate arrays
    mold = size(MESAold,dim=2) 
    m = size(MESAfile,dim=2) 
    allocate(MESAold_int(5,m))

    
    !! interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),m)
    MESAold_int(3,:) = interpolate(MESAold(3,:),MESAold(7,:),mold,MESAfile(3,:),m)
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),m)
    MESAold_int(1,:) = MESAfile(3,:)
    MESAold_int(5,:) = F_total(2,:)!interpolate(F_total(3,:),F_total(2,:),size(F_total,dim=2),MESAfile(3,:),m)

    
    !! include r=0 in the arrays
    call getradiativeR(ri,rf,MESAfile)
    mcn = m!rf-100 !! how many points before convective boundary
    allocate(Mo(5,mcn),Mn(5,mcn),r_dot(mcn))
      Mo(1,:) = MESAold_int(1,:) !mass
      Mo(2,:) = MESAold_int(2,:) !radius
      Mo(3,:) = MESAold_int(3,:) !rho
      Mo(4,:) = MESAold_int(4,:) !omega

      Mn(1,:) = MESAfile(3,:) !mass
      Mn(2,:) = MESAfile(2,:) !radius
      Mn(3,:) = MESAfile(7,:) !rho
      Mn(4,:) = MESAold_int(4,:) !omega

      do i=1,mcn
        !! compute rdot
        r_dot(i) = (Mn(2,i) - Mo(2,i) )/dt
        if (i>rf) then
          Mo(5,i) = 0d0
          Mn(5,i) = 0d0
        else
          Mo(5,i) = MESAold_int(5,i) !jdot
          Mn(5,i) = MESAold_int(5,i) !jdot
        end if
      end do

  !! compute rotation profile --------------------------------------------------------------
  nu_diff = 1d7!1d7 !1.1d3 !min for stability
  iter = 100
  delta_t = 0.01d0*dt
  omega_rad = CrankNicolson(Mn,mcn,iter,delta_t,nu_diff,r_dot)

  !convective envelope stays constant omega
  omega = omega_rad
  
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    write(203,*) MESAfile(2,1)/R_sun,MESAold_int(2,1)/R_sun,0.,0.,omega(1),MESAold_int(4,1),&
      &0d0
    do j=2,m
      write(203,*) MESAfile(2,j)/R_sun,MESAold_int(2,j)/R_sun,0.,0.,omega(j),MESAold_int(4,j),&
      &(omega(j)-omega(j-1))/(MESAfile(2,j)/R_sun - MESAfile(2,j-1)/R_sun)
      !MESAold_int(3,j),&
      !&(MESAfile(2,j) - MESAold_int(2,j))/dt,MESAold_int(2,j+1) - MESAold_int(2,j),&
    end do
    close(203)

    !! checking if AM is conserved
    AMi = 0.
    AMf = 0.
    o1 = 0.
    o2 = 0.
    do j=2,mcn
      AMi = AMi + MESAold_int(2,j)**2 * MESAold_int(4,j)
      AMf = AMf + MESAfile(2,j)**2 *omega(j)
      o1 = o1 + MESAold_int(4,j)
      o2 = o2 + omega(j)
      !write(203,*) MESA_unigrid(2,j),0,0,0,omega_rad(j)
    end do

    print*, 'delta AM  = ', (AMf - AMi)/AMi, o1-o2
    print*,'AMi = ',AMi
    print*,'AMf = ',AMf
    print*, ' dt[s] = ',delta_t,' dt[yr] = ',delta_t/(365d0*24d0*60d0*60d0)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    
    deallocate(MESAold_int,Mo,Mn,omega_rad,r_dot)

  end subroutine rotation_profile_CN

  !! ---------------------------------------------------------------------------
  !! --------------------------------- Henyey Scheme ---------------------------
  !! ---------------------------------------------------------------------------
  subroutine rotation_profile_Henyey(MESAold,MESAfile,omega,model,F_total)
    implicit none
    real (DP) :: delta_t,nu_diff,C4
    real (DP), intent(in) :: F_total(:,:),MESAfile(:,:),MESAold(:,:)
    real (DP), intent(inout) :: omega(:)
    real (DP), allocatable :: MESAold_int(:,:),Mn(:,:),Mo(:,:)
    integer :: ri,rf,mold,j,i,niter,mmesa
    integer, intent(in) :: model
    !! relaxation scheme variables
    integer :: ne,m,nb,nci,ncj,nck,nsi,nsj,nyj,nyk,itmax
    integer, allocatable :: indexv(:)
    real (DP) :: conv, slowc,AMi,AMf
    real (DP), allocatable :: c(:,:,:),s(:,:),scalv(:),y(:,:)

    !allocate arrays
    mold = size(MESAold,dim=2)
    mmesa = size(MESAfile,dim=2)
    allocate(MESAold_int(5,mmesa))

    !! interpolation of MESA files ---------------------------------------------------------
    MESAold_int(2,:) = interpolate(MESAold(3,:),MESAold(2,:),mold,MESAfile(3,:),mmesa) !radius
    MESAold_int(4,:) = interpolate(MESAold(3,:),MESAold(19,:),mold,MESAfile(3,:),mmesa) !omega
    !MESAold_int(1,:) = MESAfile(3,:) !mass
    MESAold_int(3,:) = F_total(2,:)!interpolate(F_total(3,:),F_total(2,:),size(F_total,dim=2),MESAfile(3,:),mmesa) !jdot
    MESAold_int(5,:) = 1d0!getG(model,mmesa) !g

    allocate(Mn(7,mmesa),Mo(4,mmesa))
    
      Mn(1,:) = MESAfile(3,:) !mass
      Mn(2,:) = MESAfile(2,:)/R_sun !radius
      Mn(3,:) = MESAfile(7,:) !rho
      Mn(4,:) = MESAfile(19,:) !omega
      Mn(6,:) = MESAold_int(5,:) !1d0 !interpolate(MESAold(3,:),getG(model-1,mold),mold,MESAfile(3,:),mmesa) !g
      Mn(7,:) = (MESAfile(3,:)/M_sun)**(2d0/3d0) !nu
    
      Mo(1,:) = MESAfile(3,:) !mass
      Mo(2,:) = MESAold_int(2,:)/R_sun !radius
      Mo(4,:) = MESAold_int(4,:) !omega

      call getradiativeR(ri,rf,MESAfile)
      do i=1,mmesa
        if (i>rf) then
          Mn(5,i) = 0d0
        else
          Mn(5,i) = MESAold_int(3,i) !jdot
        end if
      end do
    


    !! compute rotation profile --------------------------------------------------------------
    ri=0.
    rf=0.
    call getradiativeR(ri,rf,MESAfile)
    omega = 0.
    
    !! relaxation method ------------------------------------------------------------------
    m = mmesa!mmesa-1 !rf !! radiative part !! total of mesh points
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
    conv=1d-10!!Cesam2k20 !5d-6  !! the convergence criterion
    slowc=1d0    !! controls the fraction of corrections actually used after each iteration.
    nu_diff = 1d7!1d7 !! difusion coefficient
    delta_t = 0.001d0*dt!0.001d0*dt  !! timestep
    niter = 1000
    
    !! indexv(j) describes which column of s(i,j) the variable y(j) has been put in
    !! The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv
    indexv(1)=2 
    indexv(2)=1

    !! variables - initial guesses !!! change slightly so it doesnt cancel and get singular matrix
    y(1,1)=Mo(4,1)
    y(1,m)=Mo(4,m)
    y(2,1)=0d0
    y(2,m)=0d0
    do j=2,m-1
      y(1,j)=Mo(4,j) !! omega

      !C4 = (16d0*PI*(R_sun**4) *((0.5d0*(Mn(2,j)+Mn(2,j-1)))**4)&
      !& *(0.5d0*(Mn(3,j)+Mn(3,j-1))))&
      !&/(9d0*M_sun*(0.5d0*(Mn(6,j) + Mn(6,j-1)))*&
      !&dsqrt(0.5d0*(Mn(7,j)+Mn(7,j-1)))*(Mn(7,j)-Mn(7,j-1)))

      !C4 = 1d0/(Mn(7,j)-Mn(7,j-1))
      !y(2,j)=C4*(y(1,j)-y(1,j-1)) !! derivative of omega
      y(2,j)=(y(1,j)-y(1,j-1))/(dsqrt(0.5d0*(Mn(7,j)+Mn(7,j-1)))*(Mn(7,j)-Mn(7,j-1)))

      write(555,*) Mn(2,j), y(1,j), y(2,j),y(1,j)-y(1,j-1),C4
      
    end do
   
    do i=1,niter
      if (i/=1) then
        Mo(4,:) = y(1,:) !omega
        y(2,1) = 0d0
        y(2,m) = 0d0
        do j=2,m-1
          y(2,j)=(y(1,j)-y(1,j-1))/(dsqrt(0.5d0*(Mn(7,j)+Mn(7,j-1)))*(Mn(7,j)-Mn(7,j-1)))
        end do
      end if
      !! solve differential equations using relaxation scheme
      call solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m,y,c,s,&
      &Mn,Mo,delta_t,nu_diff)
    end do

    omega = y(1,:)

    !! print final rotation profile for this model
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    write(203,*) 0d0,0d0,0.,0.,omega(1),0d0 ,y(2,1),0d0
    do j=2,m
      write(203,*) Mn(2,j),Mo(2,j),0.,0.,omega(j),Mo(4,j) ,y(2,j),&
      &(omega(j)-omega(j-1))/(Mn(2,j)-Mn(2,j-1))
      !&(omega(j)-omega(j-1))/(dsqrt(0.5d0*(Mn(7,j)+Mn(7,j-1)))*(Mn(7,j)-Mn(7,j-1)))
    end do
    close(203)

    !! checking if AM is conserved
    AMi = 0d0
    AMf = 0d0
    do j=1,m
      AMi = AMi + ((Mo(2,j)*R_sun)**2) * Mo(4,j)
      AMf = AMf + ((Mn(2,j)*R_sun)**2) *omega(j)
    end do
  
    write(*,555) 'delta AM rad zones', (AMf - AMi)/AMi
    write(*,556) ' AMini', AMi, '   AMfin', AMf 
    write(*,556) ' dt[s]',delta_t,'  dt[yr]',delta_t/(365d0*24d0*60d0*60d0)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    555 format(A19, ES14.3)
    556 format(A6, ES14.3,A8,ES14.3)

    deallocate(MESAold_int,c,s,scalv,indexv,y,Mn,Mo)
  end subroutine rotation_profile_Henyey


  function getG(model,m) result(G)
    implicit none
    real(DP), allocatable :: MESAprofile(:,:)
    real(DP) :: G(m)
    integer, intent(in) :: model,m

    !m = getnlines('MESA/profile'//trim(string(model))//'.data')
    allocate(MESAprofile(37,m))

    call readfile('MESA/profile'//trim(string(model))//'.data',MESAprofile)

    G = MESAprofile(12,:)
    deallocate(MESAprofile)
  
  end function getG

  



  

end module rotationprofile_mod
