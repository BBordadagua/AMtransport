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
    use compute_modeamplitude
    use data_handling
    use math_functions
    !use meridionalcirculation_mod
    use rotationprofile_mod
    use numerical_scheme
    use compute_mixedmodesflux

  IMPLICIT NONE
  
  !variable declaration
  integer :: NMESA,model,dmodel
  integer :: status,modelini,modelend,version,mprofile,mnew,GYRE_var !,ri,rf
  real (DP), allocatable :: MESAhistory(:,:),MESAprofile(:,:),omega(:),F_total(:,:),MESAnewprofile(:,:)
  real (DP) ::  dummyvector(5)
  character(35) :: dummy

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
    status = SYSTEM('rm -r input/*')
  else if (GYRE_var /= 0 .AND. GYRE_var /= 1 .AND. GYRE_var /= 2) then
    print*, 'select only 0, 1 or 2'
    stop
  end if

  !! only run GYRE
  if (GYRE_var ==2) then
    status = SYSTEM('rm -r input/*')
    do model = modelini,modelend,dmodel
      call insertcol_MESAfile(model,0,dummyvector)
      call runGYRE(model)
    end do
    stop
  end if


  !! read MESA history file -----------------------------------------------------
  NMESA = getnlines('MESA/history.data')
  allocate(MESAhistory(67,NMESA), STAT = status)
  IF(status/=0) STOP
  call readfile('MESA/history.data',MESAhistory)
  !! get models with the same position in the HR diagram
  !modelini = get_modelN(MESAhistory) 
  !modelend = modelini+4!1
  !dmodel = 1

  !! insert ini rotation profile in MESA file
  call set_variables(MESAhistory,modelini)
  call insertcol_MESAfile(modelini,1,dummyvector)

  !! read MESAfirst current profile
  mprofile = 5+getnlines('MESA/profile'//trim(string(modelini))//'.data.GYRE')
  allocate(MESAprofile(19,mprofile))
  open(301, action='read',file = 'MESA/profile'//trim(string(modelini))//'.data.GYRE')
  read(301,*) dummy
  read(301,*) MESAprofile
  close(301)

  !! set general variables of the model
  call set_variables(MESAhistory,modelini)


  !! loop to get data of each model ----------------------------------------------
  do model = modelini,modelend,dmodel

    !! read next mesa profile
    mnew = 5+getnlines('MESA/profile'//trim(string(model+dmodel))//'.data.GYRE')
    allocate(MESAnewprofile(19,mnew),omega(mnew))
    open(200, action='read',file = 'MESA/profile'//trim(string(model+dmodel))//'.data.GYRE')
    read(200,*) dummy
    read(200,*) MESAnewprofile
    close(200)

    !! compute mixed modes flux
    allocate(F_total(3,mnew))
    if (GYRE_var == 0) then
      call runGYRE(model) !! run GYRE
      F_total = 0d0
      call compute_flux(model,dmodel,MESAnewprofile,F_total)
    else
      open(2000, action='read',file = 'mixed_modes/total_flux_'//trim(string(model+dmodel))//'.txt')
      read(2000,*) F_total
      close(2000)
    end if 
    

    call set_variables(MESAhistory,model+dmodel)

    !! compute new rotation profile
    !omega = 0d0
    !call rotation_profile_CN(MESAprofile,MESAnewprofile,omega,model+dmodel,F_total)
    call rotation_profile_Henyey(MESAprofile,MESAnewprofile,omega,model+dmodel,F_total)

    !! insert new rotation profile in the next MESA profile
    call insertcol_MESAfile(model+dmodel,3,omega)
    
    !free the memory
    deallocate(MESAprofile)
    allocate(MESAprofile(19,mnew))
    MESAprofile = MESAnewprofile !! current model is the previous next model
    MESAprofile(19,:) = omega
    mprofile = mnew

    deallocate(omega,F_total,MESAnewprofile)
    
  end do

  deallocate(MESAprofile)
  deallocate(MESAhistory)

end program main
