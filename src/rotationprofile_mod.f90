Module rotationprofile_mod
  use parameters
  use data_handling
  use math_functions
  use numerical_scheme
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: insertcol_GYREfile, get_r_and_w, iniprofile,get_modelN
  public :: rotation_profile_Henyey, rotation_profile_CN_new
  public :: get_MESAfiles, check_AMconservation,set_initialguess
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
    allocate(Mn(33,m),Mo(6,m))
    call get_MESAfiles(Mn,Mo,model,dmodel,omega_old,jmodes)

    !! print initial rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    write(203,*) 0d0,0d0,0.,0.,Mo(4,1)
    do j=2,m
      write(203,*) 10**Mo(2,j),0.,0.,0.,Mo(4,j)
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
    Mo(6,:) = interpolate(jmodes(3,:)/M_sun,jmodes(2,:),size(jmodes,dim=2),Mn(2,:),mnew) !jmodes
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
    allocate(Mn(33,mmesa),Mo(6,mmesa))
    call get_MESAfiles(Mn,Mo,model,dmodel,omega_old,jmodes)
    
    !! relaxation method ------------------------------------------------------------------
    m = mmesa!mmesa-1 !rf !! radiative part !! total of mesh points
    ne = 2 !! The problem involves ne equations for ne adjustable dependent variables at each point
    nyj = ne !! number of dependent variables
    nyk = m !! total of mesh points
    k1 = 1 !! first point
    k2 = m !! last point

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
    !indexv(3)=3
    !indexv(4)=4
    !indexv(5)=5
    
    !! compute constants
    allocate(var(50,m),cte(55,m)) !! number of constants
    call variables(Mn,Mo,delta_t,var,y)
    dr = (var(3,:) - var(8,:))/niter !! divide the radius into smaller steps according to smaller timesteps
    old_radius = var(8,:)
    var(3,:) = old_radius + dr
  
    call set_initialguess(k1,k2,cte,var,y,model)
   
    !! iterations for smaller timesteps
    do i=1,niter 
      if (i/=1) then
        var(8,:) = var(3,:)
        var(3,:) = old_radius + i*dr
        !! previous omega result is the new omega initial
        var(9,:) = y(1,:) !omega
        y(2,k1) = 0d0
        y(2,k2) = 0d0
        do j=k1+1,k2-1
          y(2,j)=(y(1,j)-y(1,j-1))/cte(1,j)
        end do
      end if
      
      var(31,:) = 0d0    !! lambda old
      var(19,:) = 0d0 !(2d0*var(3,:)*var(3,:)*var(9,:))/(3d0*var(6,:))    !! theta_old 
      !call constants(1,m,cte,var,y) !! update constants
      !! solve differential equations using relaxation scheme
      call solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m,y,c,s,cte,var)
    end do

    omega = y(1,:)

    !! print final rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/rot_'//trim(string(model+dmodel))//'.txt')
    write(203,*) 0d0,0d0,0.,0.,omega(k1),var(9,k1) ,y(2,k1),0d0
    do j=k1+1,k2
      write(203,*) var(3,j),var(8,j),0.,0.,omega(j),var(9,j) ,y(2,j),&
      &(omega(j)-omega(j-1))/(var(3,j)-var(3,j-1))
    end do
    close(203)

    call check_AMconservation(Mo,Mn,omega,delta_t)

    deallocate(c,s,scalv,indexv,y,Mn,Mo,dr,old_radius)
    deallocate(cte,var)
  end subroutine rotation_profile_Henyey

  subroutine set_initialguess(k1,k2,cte,var,y,model)
    implicit none
    real (DP), intent(inout) :: cte(:,:),y(:,:)
    real (DP), intent(in) :: var(:,:)
    integer, intent(in) :: k1,k2,model
    integer :: j
    !real (DP),allocatable :: y1(:),y2(:),y3(:),y4(:),y5(:)

    !allocate(y1(size(y,dim=2)),y2(size(y,dim=2)),y3(size(y,dim=2)),y4(size(y,dim=2)),y5(size(y,dim=2)))

    call constants(k1,k2,cte,var,y)

    !! variables - initial guesses 
    !!! change slightly so it doesnt cancel and get singular matrix
    
    !! omega
    y(1,:)=var(9,:) 

    !! Theta = 1/nu12 *domega/dnu
    y(2,k1)=0d0
    y(2,k2)=0d0
    do j=k1+1,k2-1 
      y(2,j)=(y(1,j)-y(1,j-1))/cte(1,j) 
    end do

    !! U2
    !y(3,:) = 0d0 

    !! A
    !y(4,:) = 0d0 

    !! Lambda
    !y(5,:) = 0d0 

    !deallocate(y1,y2,y3,y4,y5)

    !! print initial rotation profile for this model ------------------------------------------
    open(203, action='write',file = 'output/rot_'//trim(string(model))//'.txt')
    write(203,*) 0d0,0d0,0.,0.,var(9,k1)
    do j=k1+1,k2
      write(203,*) var(8,j),0.,0.,0.,var(9,j)
    end do
    close(203)
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
