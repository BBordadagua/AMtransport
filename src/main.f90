!!
!!                 mixedmodesflux.f90
!!    Explanation:
!!     
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!


program main
    use parameters
    use coefficient_equations
    use data_handling
    use math_functions
    use rotationprofile_mod
    use numerical_scheme
    use compute_mixedmodesflux

  IMPLICIT NONE
  
  !variable declaration
  integer :: NMESA,model,dmodel,mmodes,count
  integer :: status,modelini,modelend,version,mprofile,GYRE_var
  real (DP), allocatable :: MESAhistory(:,:),GYREfile(:,:)
  real (DP), allocatable :: y(:,:),F_total(:,:),y_old(:,:)
  real (DP) ::  dummyvector(5)

  !! user settings ----------------------------------------------------------------
  print*, 'Give initial, final model number, and delta model (e.g., 0 500 10):'
  read*, modelini, modelend, dmodel

  !print*, 'Use simplified eq. (0) or full eq. (1):'
  !read*, version
  version = 0
  if (version /= 0 .AND. version /= 1) then
    print*, 'select only 0 or 1'
    stop
  end if

  GYRE_var = 1
  !print*, 'Run with GYRE (0), without GYRE (1) or only GYRE (2):'
  !read*, GYRE_var
  if (GYRE_var ==0) then
    !status = SYSTEM('rm -r input/*')
  else if (GYRE_var /= 0 .AND. GYRE_var /= 1 .AND. GYRE_var /= 2) then
    print*, 'select only 0, 1 or 2'
    stop
  end if

  !! only run GYRE
  if (GYRE_var ==2) then
    status = SYSTEM('rm -r input/*')
    do model = modelini,modelend,dmodel
      call insertcol_GYREfile(model,0,dummyvector)
      call runGYRE(model)
    end do
    stop
  end if


  !! read MESA history file -----------------------------------------------------
  NMESA = getnlines('MESA/history.data')
  allocate(MESAhistory(59,NMESA), STAT = status)
  IF(status/=0) STOP
  call readfile('MESA/history.data',MESAhistory)
  !! get models with the same position in the HR diagram
  !modelini = get_modelN(MESAhistory) 
  !modelend = modelini+4!1
  !dmodel = 1

  !! set general variables of the model
  call set_parameters(MESAhistory,modelini)

  mprofile = 5+getnlines('MESA/profile'//trim(string(modelini))//'.data.GYRE')
  allocate(y_old(5,mprofile),GYREfile(18,mprofile))
  y_old = 0d0
  y_old(1,:) = 3e-5!iniprofile(GYREfile,mprofile,modelini,1)!! compute initial rotation profile
  count = 0
  
  
  !! loop to get data of each model ----------------------------------------------
  do model = modelini,modelend,dmodel

    if (GYRE_var == 0) then
      !! generate an initial rotation profile
      !call insertcol_GYREfile(modelini,1,dummyvector) !! equations do not have omega
      !call runGYRE(model) !! run GYRE
      call compute_flux(model,dmodel) !! compute mixed modes flux
    else
      mprofile = getnlines('MESA/profile'//trim(string(model+dmodel))//'.data')
      allocate(y(5,mprofile))
      
      !! get jmodes from file
      !mmodes = 6+getnlines('mixed_modes/total_flux_'//trim(string(model+dmodel))//'.txt')
      !allocate(F_total(3,mmodes))
      !open(2000, action='read',file = 'mixed_modes/total_flux_'//trim(string(model+dmodel))//'.txt')
      !read(2000,*) F_total
      !close(2000)
      allocate(F_total(3,100)) !! For MS stars
      F_total = 1d0 !! For MS stars
      

      call set_parameters(MESAhistory,model+dmodel)
      
      !! compute new rotation profile
      call rotation_profile_relaxationscheme_tot(y,model,dmodel,F_total,y_old,count)
      !call rotation_profile_Henyey(y(1,:),model,dmodel,F_total,y_old(1,:))
      !call rotation_profile_CN_new(y(1,:),model,dmodel,F_total,y_old(1,:))
      
      deallocate(y_old)
      allocate(y_old(5,mprofile))
      y_old = y
      deallocate(F_total,y)
      
    end if
  end do

  if (GYRE_var /= 0) deallocate(y_old)

  deallocate(MESAhistory,GYREfile)

end program main
