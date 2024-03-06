Module rotationprofile_mod
  use parameters
  use data_handling
  use math_functions
  use numerical_scheme
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: insertcol_GYREfile, get_r_and_w, iniprofile,get_modelN
  public :: rotation_profile_Henyey, rotation_profile_CN_new, rotation_profile_relaxationscheme_tot
  public :: get_MESAfiles, check_AMconservation,set_initialguess,get_MESAfiles_tot
  contains

  function get_modelN(MESAhistory) result(model)
    implicit none
    real(DP), intent(in) :: MESAhistory(:,:)
    real(DP) ::logL,bumpL,d
    integer :: m, i,model,start
 
    start = 500  !avoid MS where luminosity decreases
    m = size(MESAhistory,dim=2)
    model = start
    bumpL=0d0

    logL = MESAhistory(38,start-1)
    do i=start,m
      if(logL> MESAhistory(38,i)) then
        bumpL = MESAhistory(38,i)
        EXIT
      end if
      logL = MESAhistory(38,i)
    end do

    d=100d0
    do i=start,m
      if (abs(MESAhistory(38,i) - 0.85d0*bumpL) <d) then
        model = i
        d = abs(abs(MESAhistory(38,i) - 0.85d0*bumpL))
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

    r_t = 0d0
    w=0d0
    r2=0d0
    r1=0d0
    max = 0d0
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
    
    max = max/3d0
    d=10000d0
    do i=j,1,-1
      if (abs(MESAfile(15,i) - max) <d) then
        r1 = MESAfile(1,i)
        d = abs(MESAfile(15,i) - max)
      end if
    end do
    d=10000d0
    do i=j,m
      if (abs(MESAfile(15,i) - max) <d) then
        r2 = MESAfile(1,i)
        d = abs(MESAfile(15,i) - max)
      end if
    end do

    w=(r2-r1)/R_star

  end subroutine get_r_and_w

  function iniprofile(MESAfile,m,model,value) result(rot)
    implicit none
    integer, intent(in) :: m,model,value
    real(DP) :: Omega_s,Omega_c,r_t,w, rot(m), r(m), array(5)
    real(DP), intent(inout) :: MESAfile(:,:)
    integer, allocatable :: a1(:)
    integer :: i

    
    if (value == 1) then !! read file from txt
      allocate(a1(m))
      open(2000, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
      read(2000,*) array
      do i=1,m
        read(2000,*) a1(i), MESAfile(:,i)
      end do
      close(2000)
      deallocate(a1)
    end if
    

    !end of MS constant profile
    !omega = 400.e-9 * 2*PI

    !!Belkacem et al. (2015b) table 2. values

    call get_r_and_w(MESAfile,r_t,w)

    Omega_s = 0d0
    Omega_c = 0d0


    if (model < 400) then
      !!SG model
      Omega_s = 200d0*2d0*PI*1d-9
      Omega_c = 600d0*2d0*PI*1d-9
    else if (model > 400 .AND. model < 500) then
      !!Early RG model
      Omega_s = 55d0*2d0*PI*1d-9
      Omega_c = 1000d0*2d0*PI*1d-9
    else if (model > 500) then
      !!RG model
      Omega_s = 8.1d0*2d0*PI*1d-9
      Omega_c = 500d0*2d0*PI*1d-9
    end if

    !!RC model  !Deheuvels et al. (2015) Table 6
    !Omega_s = 70.*2.*PI*1e-9
    !Omega_c = 160.*2.*PI*1e-9


    r = MESAfile(1,:)/R_star
    rot = Omega_s + (Omega_c-Omega_s)*(1d0+ERF((r_t-r)/w))/2d0
    
  end function iniprofile


  subroutine insertcol_GYREfile(model,number,omega)
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
      filearray(18,:) = iniprofile(filearray,N,model,0)
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

    deallocate(filearray,a1)
  
  end subroutine insertcol_GYREfile


  !! ---------------------------------------------------------------------------
  !! ---------------------- Crank-Nicolson at next model non-uniform grid ---------------------------
  !! ---------------------------------------------------------------------------
 
  subroutine rotation_profile_CN_new(omega,model,dmodel,jmodes,omega_old)
    implicit none
    real (DP) :: delta_t
    real (DP), intent(inout) :: omega(:)
    real (DP), intent(in) :: jmodes(:,:),omega_old(:)
    real (DP), allocatable :: Mo(:,:),Mn(:,:)
    integer :: j,iter,m
    integer, intent(in) :: model,dmodel

    !allocate and fill arrays
    m = getnlines('MESA/profile'//trim(string(model+dmodel))//'.data')
    allocate(Mn(35,m),Mo(6,m))
    call get_MESAfiles(Mn,Mo,model,dmodel,omega_old,jmodes)

    !! print initial rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/inirot_'//trim(string(model))//'.txt')
    write(203,*) 0d0,0d0,0.,0.,0.,Mo(4,1)
    do j=2,m
      write(203,*) 10**Mo(2,j),0.,0.,0.,0.,Mo(4,j)
    end do
    close(203)

    !! compute rotation profile --------------------------------------------------------------
    iter = 100
    delta_t = dt/iter
    omega = CrankNicolson_new(Mn,m,iter,delta_t,Mo)
  
    open(203, action='write',file = 'output/rot_'//trim(string(model+dmodel))//'.txt')
    write(203,*) 0.,0.,0.,0.,omega(1),Mo(4,1)!,0d0
    do j=2,m
      write(203,*) 10**Mn(3,j),10**Mo(2,j),0.,0.,omega(j),Mo(4,j)!,&
      !&(omega(j)-omega(j-1))/(10**Mn(3,j) - 10**Mn(3,j-1))
    end do
    close(203)

    !! checking if AM is conserved
    call check_AMconservation(Mo,Mn,omega,delta_t)
    
    deallocate(Mo,Mn)

  end subroutine rotation_profile_CN_new

  subroutine get_MESAfiles(Mn,Mo,model,dmodel,omega,jmodes)
    implicit none
    real (DP), intent(inout) :: Mn(:,:), Mo(:,:)
    real (DP), allocatable :: Mo_full(:,:),Mn_original(:,:),Mo_original(:,:)
    real (DP), intent(in) :: omega(:), jmodes(:,:)
    integer, intent(in) :: model, dmodel
    integer :: mold,mnew,i,ri,rf

    !! MESA profile from previous timestep
    mold = getnlines('MESA/profile'//trim(string(model))//'.data')
    allocate(Mo_full(size(Mn,dim=1),mold),Mo_original(size(Mn,dim=1),mold))
    call readfile('MESA/profile'//trim(string(model))//'.data',Mo_original)

    !! MESA profile from current timestep
    mnew = size(Mn,dim=2)
    allocate(Mn_original(size(Mn,dim=1),mnew))
    call readfile('MESA/profile'//trim(string(model+dmodel))//'.data',Mn_original)

    !! invert MESA files to start in r=0 
    do i=1,mnew
      Mn(:,i) =  Mn_original(:,mnew+1-i)
    end do
    do i=1,mold
      Mo_full(:,i) =  Mo_original(:,mold+1-i)
    end do
    
    !! interpolation of MESA files
    Mo(1,:) = Mn(2,:) !mass !(M_sun)
    Mo(2,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(3,:),mold,Mn(2,:),mnew)) !log radius
    Mo(3,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(5,:),mold,Mn(2,:),mnew)) !log rho
    Mo(4,:) = interpolate(Mo_full(2,:),omega,mold,Mn(2,:),mnew) !omega
    Mo(5,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(10,:),mold,Mn(2,:),mnew)) !log gsurf
    Mo(6,:) = 5d0!interpolate(jmodes(3,:)/M_sun,jmodes(2,:),size(jmodes,dim=2),Mn(2,:),mnew) !jmodes
    !! add N2
    call getradiativeR(ri,rf,Mn,1)
    Mo(6,1) = 0d0
      do i=1,mnew
        if(i<ri .or. i>rf) then
          Mo(6,i) = 0d0 !! mixed modes not valid in convection zones
        end if
      end do
     
    deallocate(Mo_full,Mn_original,Mo_original)

  end subroutine get_MESAfiles

  subroutine get_MESAfiles_tot(Mn,Mo,model,dmodel,y,jmodes)
    implicit none
    real (DP), intent(inout) :: Mn(:,:), Mo(:,:)
    real (DP), allocatable :: Mo_full(:,:),Mn_original(:,:),Mo_original(:,:)
    real (DP), intent(in) :: y(:,:), jmodes(:,:)
    integer, intent(in) :: model, dmodel
    integer :: mold,mnew,i,ri,rf

    !! MESA profile from previous timestep
    mold = getnlines('MESA/profile'//trim(string(model))//'.data')
    allocate(Mo_full(size(Mn,dim=1),mold),Mo_original(size(Mn,dim=1),mold))
    call readfile('MESA/profile'//trim(string(model))//'.data',Mo_original)

    !! MESA profile from current timestep
    mnew = size(Mn,dim=2)
    allocate(Mn_original(size(Mn,dim=1),mnew))
    call readfile('MESA/profile'//trim(string(model+dmodel))//'.data',Mn_original)

    !! invert MESA files to start in r=0 
    do i=1,mnew
      Mn(:,i) =  Mn_original(:,mnew+1-i)
    end do
    do i=1,mold
      Mo_full(:,i) =  Mo_original(:,mold+1-i)
    end do
    print*, 'here 1'
    
    !! interpolation of MESA files
    Mo(1,:) = Mn(2,:) !mass !(M_sun)
    Mo(2,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(3,:),mold,Mn(2,:),mnew)) !log radius
    Mo(3,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(5,:),mold,Mn(2,:),mnew)) !log rho
    Mo(4,:) = interpolate(Mo_full(2,:),y(1,:),mold,Mn(2,:),mnew) !omega y1
    Mo(5,:) = dlog10(interpolate(Mo_full(2,:),10**Mo_full(10,:),mold,Mn(2,:),mnew)) !log gsurf
    Mo(6,:) = 5d0 !interpolate(jmodes(3,:)/M_sun,jmodes(2,:),size(jmodes,dim=2),Mn(2,:),mnew) !jmodes
    Mo(7,:) = interpolate(Mo_full(2,:),y(2,:),mold,Mn(2,:),mnew) !y2
    Mo(8,:) = interpolate(Mo_full(2,:),y(3,:),mold,Mn(2,:),mnew) !y3
    Mo(9,:) = interpolate(Mo_full(2,:),y(4,:),mold,Mn(2,:),mnew) !y4
    Mo(10,:) = interpolate(Mo_full(2,:),y(5,:),mold,Mn(2,:),mnew) !y5
    !! add N2
    call getradiativeR(ri,rf,Mn,1)
    Mo(6,1) = 0d0
      do i=1,mnew
        if(i<ri .or. i>rf) then
          Mo(6,i) = 0d0 !! mixed modes not valid in convection zones
        end if
      end do
      print*, 'here 2'
     
    deallocate(Mo_full,Mn_original,Mo_original)

  end subroutine get_MESAfiles_tot

  !! ---------------------------------------------------------------------------
  !! --------------------------------- Henyey Scheme ---------------------------
  !! ---------------------------------------------------------------------------
  subroutine rotation_profile_Henyey(omega,model,dmodel,jmodes,omega_old)
    implicit none
    real (DP) :: delta_t
    real (DP), intent(in) :: jmodes(:,:),omega_old(:)
    real (DP), intent(inout) :: omega(:)
    real (DP), allocatable :: Mn(:,:),Mo(:,:)
    integer :: j,i,niter,mmesa
    integer, intent(in) :: model,dmodel
    !! relaxation scheme variables
    integer :: ne,m,nb,nci,ncj,nck,nsi,nsj,nyj,nyk,itmax,k1,k2
    integer, allocatable :: indexv(:)
    real (DP) :: conv, slowc
    real (DP), allocatable :: c(:,:,:),s(:,:),scalv(:),y(:,:), old_radius(:),dr(:)
    real (DP), allocatable :: cte(:,:),var(:,:)

    
    !allocate and fill arrays
    mmesa = getnlines('MESA/profile'//trim(string(model+dmodel))//'.data')
    allocate(Mn(35,mmesa),Mo(6,mmesa))
    call get_MESAfiles(Mn,Mo,model,dmodel,omega_old,jmodes)

    Mo(4,:) = 3e-5!3d-7 !! omega old constant !! changed
    

    !! print initial rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/inirot_'//trim(string(model))//'.txt')
    write(203,*) 0.,0.,0.,0.,0.,Mo(4,1) 
    do j=1,mmesa
      write(203,*) 10**Mo(2,j),0.,0.,0.,0.,Mo(4,j) !,y(2,j),&
      !&(omega(j)-omega(j-1))/(var(3,j)-var(3,j-1))
    end do
    close(203)
    
    !! defining the begining and ending point of integration
    call getradiativeR(k1,k2,Mn,1)
    print*,k1,k2
    !k1 = 1 !! first point
    !k2 = mmesa !! last point
    !! relaxation method ------------------------------------------------------------------
    m = k2 !! radiative part !! total of mesh points
    ne = 2 !! The problem involves ne equations for ne adjustable dependent variables at each point
    nyj = ne !! number of dependent variables
    nyk = m !! total of mesh points

    itmax=1000!100   !! maximum number of iterations
    conv=1d-10!!Cesam2k20 !5d-6  !! the convergence criterion
    slowc=1d0    !! controls the fraction of corrections actually used after each iteration.
    niter = 10 !! iterations for smaller timesteps
    delta_t = dt/niter !! timestep

    !! The arrays c(1:nci,1:ncj,1:nck), s(1:nsi,1:nsj) supply dummy storage used by the relaxation code
    !! the minimum dimensions must satisfy: nci=ne, ncj=ne-nb+1, nck=m+1, nsi=ne, nsj=2*ne+1.
    nb=1  !! The nb boundary conditions at the first mesh point
    nci=ne
    ncj=ne-nb+1
    nck=m+1
    nsi=ne
    nsj=2*ne+1

    !! allocate arrays
    allocate(indexv(ne),c(nci,ncj,nck),s(nsi,nsj),scalv(ne),y(ne,m),dr(m),old_radius(m))

    !! indexv(j) describes which column of s(i,j) the variable y(j) has been put in
    !! The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv
    indexv(1)=2 
    indexv(2)=1
  
    !! compute constants
    allocate(var(50,m),cte(60,m)) !! number of constants
    call variables(k1,k2,Mn,Mo,delta_t,var)
    dr = (var(3,:) - var(8,:))/niter !! divide the radius into smaller steps according to smaller timesteps
    old_radius = var(8,:)
    var(3,:) = old_radius + dr
  
    call set_initialguess(k1,k2,cte,var,y,0)
   
    !! iterations for smaller timesteps
    do i=1,niter 
      if (i/=1) then
        !! previous omega result is the new omega initial
        var(9,:) = y(1,:) !omega
        !y(2,k1) = 0d0
        !y(2,k2) = 0d0
        !do j=k1+1,k2-1
        !  y(2,j)=(y(1,j)-y(1,j-1))/cte(1,j)
        !end do

        var(8,:) = var(3,:)
        var(3,:) = old_radius + i*dr
      end if
      
      !call constants(1,m,cte,var,y) !! update constants
      !! solve differential equations using relaxation scheme
      call solvde(itmax,conv,slowc,scalv,indexv,ne,nb,k2,y,c,s,cte,var,k1,k2,0)
    end do

    do j=k1,mmesa
      if (j<k1) then
        omega(j) = 0d0!y(1,k1)*var(3,k1)*var(3,k1)/((10**Mn(3,j))*(10**Mn(3,j)))
      else if (j>k2) then
        omega(j) = y(1,k2)!y(1,k2)*var(3,k2)*var(3,k2)/((10**Mn(3,j))*(10**Mn(3,j)))
      else 
        omega(j) = y(1,j)
      end if
    end do

    !! print final rotation profile for this model ------------------------------------------
    open(210, action='write',file = 'output/rot_'//trim(string(model+dmodel))//'.txt')
    write(210,*) 0.,0.,0.,0.,omega(k1),Mo(4,k1)! ,y(2,k1),0d0
    do j=k1+1,mmesa
      write(210,*) 10**Mn(3,j),10**Mo(2,j),0.,0.,omega(j),Mo(4,j) !,y(2,j),&
      !&(omega(j)-omega(j-1))/(var(3,j)-var(3,j-1))
    end do
    close(210)

    call check_AMconservation(Mo,Mn,omega,delta_t)

    deallocate(c,s,scalv,indexv,y,Mn,Mo,dr,old_radius)
    deallocate(cte,var)
  end subroutine rotation_profile_Henyey


  subroutine rotation_profile_relaxationscheme_tot(y_ini,model,dmodel,jmodes,y_old,count)
    implicit none
    real (DP) :: delta_t
    real (DP), intent(in) :: jmodes(:,:),y_old(:,:)
    real (DP), intent(inout) :: y_ini(:,:)
    real (DP), allocatable :: Mn(:,:),Mo(:,:)
    integer :: j,i,niter,mmesa,k
    integer, intent(in) :: model,dmodel
    !! relaxation scheme variables
    integer :: ne,m,nb,nci,ncj,nck,nsi,nsj,nyj,nyk,itmax,k1,k2
    integer, allocatable :: indexv(:)
    real (DP) :: conv, slowc
    real (DP), allocatable :: c(:,:,:),s(:,:),scalv(:),y(:,:), old_radius(:),dr(:)
    real (DP), allocatable :: cte(:,:),var(:,:)
    integer, intent(inout) :: count
    
    !allocate and fill arrays
    mmesa = getnlines('MESA/profile'//trim(string(model+dmodel))//'.data')
    allocate(Mn(35,mmesa),Mo(10,mmesa))!,Mo(6,mmesa))
    call get_MESAfiles_tot(Mn,Mo,model,dmodel,y_old,jmodes)
    
    !Mo(4,:) = 3e-5!3d-7 !! omega old constant !! changed

    !! print initial rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/inirot_'//trim(string(model))//'.txt')
    write(203,*) 0.,0.,0.,0.,0.,Mo(4,1) 
    do j=1,mmesa
      write(203,*) 10**Mo(2,j),0.,0.,0.,0.,Mo(4,j) !,y(2,j),&
      !&(omega(j)-omega(j-1))/(var(3,j)-var(3,j-1))
    end do
    close(203)
    
    
    !! defining the begining and ending point of integration
    call getradiativeR(k1,k2,Mn,1) 
    !k1 = 1 !! first point
    !k2 = mmesa !! last point
    print*,'k1=',k1,'k2=',k2 ,10**Mn(3,k2)

    !! relaxation method ------------------------------------------------------------------
    m = k2 !! radiative part !! total of mesh points
    ne = 5 !! The problem involves ne equations for ne adjustable dependent variables at each point
    nb = 2  !! The nb boundary conditions at the first mesh point
    nyj = ne !! number of dependent variables
    nyk = m !! total of mesh points

    itmax=1000!100   !! maximum number of iterations
    conv=1d-10!!Cesam2k20 !5d-6  !! the convergence criterion
    slowc=1d0    !! controls the fraction of corrections actually used after each iteration.
    niter = 1!10 !! iterations for smaller timesteps
    delta_t = dt/niter !! timestep

    !! The arrays c(1:nci,1:ncj,1:nck), s(1:nsi,1:nsj) supply dummy storage used by the relaxation code
    !! the minimum dimensions must satisfy: nci=ne, ncj=ne-nb+1, nck=m+1, nsi=ne, nsj=2*ne+1.
    nci=ne
    ncj=ne-nb+1
    nck=m+1
    nsi=ne
    nsj=2*ne+1

    !! allocate arrays
    allocate(indexv(ne),c(nci,ncj,nck),s(nsi,nsj),scalv(ne),y(ne,m),dr(m),old_radius(m))

    !! indexv(j) describes which column of s(i,j) the variable y(j) has been put in
    !! The nb boundary conditions at the first mesh point must contain some dependence on the first nb variables listed in indexv
    indexv(1)=3 !!CHANGED
    indexv(2)=1 !!CHANGED
    indexv(3)=2 !!CHANGED
    indexv(4)=4
    indexv(5)=5
    
    !! compute constants
    allocate(var(50,m),cte(60,m)) !! number of constants
    call variables(k1,k2,Mn,Mo,delta_t,var)
    dr = (var(3,:) - var(8,:))/niter !! divide the radius into smaller steps according to smaller timesteps
    old_radius = var(8,:)
    var(3,:) = old_radius + dr
  
    if (count ==0) then
      call set_initialguess(k1,k2,cte,var,y,1)
      count = 1
    else 
      y(1,:) = Mo(4,:) !y1
      y(2,:) = Mo(7,:) !y2
      y(3,:) = Mo(8,:) !y3
      y(4,:) = Mo(9,:) !y4
      y(5,:) = Mo(10,:) !y5
      var(31,:) = y(5,:)   !! lambda old
      var(19,k1) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*((var(4,k1)*var(3,k1)**4)/var(6,k1))*y(2,k1) !! theta_old
      var(19,k1+1:k2) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*&
                      & ((0.5d0*(var(4,k1+1:k2)+var(4,k1:k2-1))*(0.5d0*(var(3,k1+1:k2)+var(3,k1:k2-1)))**4)&
                      &/(0.5d0*(var(6,k1+1:k2)+var(6,k1:k2-1))))*y(2,k1+1:k2) !! theta_old 
      var(20,:) = y(1,1) !! omega(1)
      call constants_tot(k1,k2,cte,var,y)
    end if
    

    open(1002, action='write',file='var.txt')
    open(1003, action='write',file='cte.txt')
    do i=k1,k2
      write(1002,*) var(:,i)
      write(1003,*) cte(:,i)
    end do
    close(1002)
    close(1003)
    
    !! iterations for smaller timesteps
    do i=1,niter 
      if (i/=1) then
        !! previous omega result is the new omega initial
        var(9,:) = y(1,:) !omega
        var(31,:) =  y(5,:)!! lambda old
        var(19,:) =  cte(41,:)*y(2,:)!! theta_old
        var(20,:) = y(1,1)

        var(8,:) = var(3,:)
        var(3,:) = old_radius + i*dr
      end if
      
      !call constants(1,m,cte,var,y) !! update constants
      !! solve differential equations using relaxation scheme
      call solvde(itmax,conv,slowc,scalv,indexv,ne,nb,k2,y,c,s,cte,var,k1,k2,1)
    end do

    do j=k1,mmesa
      if (j<k1) then
        y_ini(1,j) = 0d0!y(1,k1)*var(3,k1)*var(3,k1)/((10**Mn(3,j))*(10**Mn(3,j)))
        y_ini(:,j) = 0d0
      else if (j>k2) then
        y_ini(1,j) = y(1,k2)!y(1,k2)*var(3,k2)*var(3,k2)/((10**Mn(3,j))*(10**Mn(3,j)))
        y_ini(2,j) = 0d0
        y_ini(3,j) = 0d0
        y_ini(4,j) = 0d0
        y_ini(5,j) = 0d0
      else 
        y_ini(:,j) = y(:,j)
      end if
    end do

    !! print final rotation profile for this model ------------------------------------------
    open(210, action='write',file = 'output/rot_'//trim(string(model+dmodel))//'.txt')
    !write(210,*) 0.,0.,0.,0.,y_ini(1,k1),Mo(4,k1)! ,y(2,k1),0d0
    write(210,*) var(2,k1),var(3,k1), y(1,k1),y(2,k1),y(3,k1),y(4,k1),y(5,k1),cte(44,k1),cte(45,k1),cte(46,k1)
    do k=k1+1,k2!mmesa
      !write(210,*) 10**Mn(3,j),10**Mo(2,j),0.,0.,y_ini(1,j),Mo(4,j) !,y(2,j),&
      !&(omega(j)-omega(j-1))/(var(3,j)-var(3,j-1))
      write(210,*) var(2,k),var(3,k), y(1,k),y(2,k),y(3,k),y(4,k),y(5,k),cte(44,k),cte(45,k),cte(46,k)
    end do
    close(210)

    call check_AMconservation(Mo,Mn,y_ini(1,:),delta_t)

    deallocate(c,s,scalv,indexv,y,Mn,Mo,dr,old_radius)
    deallocate(cte,var)
  end subroutine rotation_profile_relaxationscheme_tot


  subroutine set_initialguess(k1,k2,cte,var,y,value)
    implicit none
    real (DP), intent(inout) :: cte(:,:),y(:,:),var(:,:)
    integer, intent(in) :: k1,k2,value
    integer :: k
    real (DP), allocatable :: y1(:),y2(:),y4(:),y5(:),psi,psik,psik_1
    real (DP) :: aux1, aux2, aux3, aux4,aux5

    !! variables - initial guesses
    open(1001, action='write',file='initial_guess.txt')
    !! omega
    y(1,:)= var(9,:) !!changed
    !var(40,:) = var(9,k2)*var(8,k2)*var(8,k2) !! constant AM at the convective zone of previous timestep

    !! Theta = 1/nu12 *domega/dnu
    y(2,k1)=0d0 !!boundary condition
    y(2,k2)=0d0 !!boundary condition
    do k=k1+1,k2-1 
      y(2,k)=(y(1,k)-y(1,k-1))/((var(1,k)-var(1,k-1))*dsqrt(var(1,k)))
    end do

    if (value == 0) then
      call constants(k1,k2,cte,var,y)
     do k=k1,k2
      write(1001,*) var(2,k),var(3,k), y(1,k),y(2,k)
     end do
    else
      !! Lambda
      y(5,:) = 0d0

      allocate(y1(size(y,dim=2)),y2(size(y,dim=2)),y4(size(y,dim=2)),y5(size(y,dim=2)))
      y1(k1) = y(1,k1) !boundary
      y2(k1) = y(2,k1) !boundary
      y5(k1) = y(5,k1) !boundary
      do k=k1+1,k2
        y1(k) = 0.5d0*(y(1,k)+y(1,k-1))
        y2(k) = 0.5d0*(y(2,k)+y(2,k-1))
        y5(k) = 0.5d0*(y(5,k)+y(5,k-1))

        var(19,k) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*&
                      & ((0.5d0*(var(4,k)+var(4,k-1))*(0.5d0*(var(3,k)+var(3,k-1)))**4)&
                      &/(0.5d0*(var(6,k)+var(6,k-1))))*y(2,k) !! theta_old
      end do
      var(19,k1) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*((var(4,k1)*var(3,k1)**4)/var(6,k1))*y(2,k1) !! theta_old
      var(20,:) = y(1,1) !! omega(1)
      var(31,:) = y(5,:)   !! lambda old
      call constants_tot(k1,k2,cte,var,y)

      
      do k=k1+1,k2
        !psi = cte(35,k)*y1(k)*cte(41,k)*y2(k) !- cte(32,k)*y5(k)
        !psik = cte(33,k)*y(1,k)*cte(42,k)*y(2,k) !- cte(30,k)*y(5,k)
        !psik_1 = cte(34,k)*y(1,k-1)*cte(43,k)*y(2,k-1)!-cte(31,k)*y(5,k-1)

        !! A
        !y(4,k1) = 0d0 !!boundary condition !!not sure
        !y(4,k) = cte(27,k)*(psik-psik_1)/cte(1,k) &
        !  & + cte(28,k)*psi/cte(1,k) &
        !  !- cte(29,k)*y5(k)/cte(1,k) &
        !  & - y1(k)*cte(41,k)*y2(k)

        !y4(k) = 0.5d0*(y(4,k)+y(4,k-1))

        y(4,:) = 0d0
        y4 = 0d0
        

        !! U2
        !aux1 = 1d0-cte(40,k)*y1(k)*y1(k)
        !!aux2 = cte(23,k)+cte(24,k)*aux1
        !aux4 = cte(13,k)*(-1d0+cte(14,k)+2d0*cte(15,k)*y1(k)*y1(k))
        !y(3,k1) = 0d0!!boundary condition
        !y(3,k) = (cte(12,k)*aux1 &
         !     & - aux4*y1(k)*y1(k) + cte(13,k)*cte(15,k)*y1(k)**4 &
              !& + cte(16,k)*(y(4,k)-y(4,k-1)) &
              !& + aux2*cte(41,k)*y1(k)*y2(k) &
              !& + cte(20,k)*y4(k) &
              !& + cte(25,k)*y5(k) &
              !& + cte(26,k)*psi &
         !     &)/(cte(11,k)*aux1)

        !! Cesam initial condition
        y(3,k) = 0d0
        aux1 = (8d0*L_sun*R_sun**5 )/(3d0*G**2 *M_sun**3)
        aux2 = (var(12,k)*var(13,k)*var(3,k)**5)/((var(2,k)**3)*var(10,k)*var(4,k)*var(11,k))
        aux3 = 1d0 - var(33,k)
        aux4 = var(14,k) - var(15,k)
        aux5 = (3d0*L_sun*var(34,k)*var(13,k)*var(12,k))/(16d0*PI*a_raddensity*c_speed*M_sun*G*var(2,k)*var(11,k)**4) !! gradrad
        if (aux5<var(14,k)) y(3,k) = -(y(1,1)**2)*aux1*aux2*aux3/aux4 !! negative was added!!
        y(3,k2) = 0d0 !!boundary condition
        
        
      if (k==k1+1) then
        write(1001,*) var(2,k1),var(3,k1), y(1,k1),y(2,k1),y(3,k1),y(4,k1),y(5,k1), &
        & cte(11,k1),&
        & 1d0/(cte(11,k1)*(1d0-cte(40,k1)*y1(k1)*y1(k1)))
      end if

      write(1001,*) var(2,k),var(3,k), y(1,k),y(2,k),y(3,k),y(4,k),y(5,k), &
      cte(11,k), 1d0/(cte(11,k)*(1d0-cte(40,k)*y1(k)*y1(k)))
      
      end do

      deallocate(y1,y2,y4,y5)
    end if

    close(1001)
  end subroutine set_initialguess

  subroutine check_AMconservation(Mo,Mn,omega,delta_t)
    implicit none
    real (DP) :: AMi,AMf
    real (DP), intent(in) :: Mo(:,:), Mn(:,:), omega(:), delta_t
    integer :: j

    !! checking if AM is conserved
    AMi = 0d0
    AMf = 0d0
    do j=1,size(Mn,dim=2)
      AMi = AMi + ((Mo(2,j)*R_sun)**2) * Mo(4,j)
      AMf = AMf + ((Mn(3,j)*R_sun)**2) *omega(j)
    end do
  
    write(*,555) 'delta AM rad zones', (AMf - AMi)/AMi
    write(*,556) ' AMini', AMi, '   AMfin', AMf 
    write(*,556) ' dt[s]',delta_t,'  dt[yr]',delta_t/(365d0*24d0*60d0*60d0)
    !print*,'dt-dtmin (stable<0) = ',delta_t - (delta_r**2)/(2.*nu_diff)
    print*,'  '
    555 format(A19, ES14.3)
    556 format(A6, ES14.3,A8,ES14.3)
  end subroutine check_AMconservation

end module rotationprofile_mod
