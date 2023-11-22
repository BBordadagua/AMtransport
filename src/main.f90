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
  integer :: NMESA,model,dmodel,i
  integer :: status,modelini,modelend,version,mprofile,ri,rf,mnew,GYRE_var
  real (DP), allocatable :: MESAhistory(:,:), MESAprofile(:,:),omega(:),F_total(:,:)
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

  !! loop to get data of each model ----------------------------------------------
  do model = modelini,modelend,dmodel

    !! read MESA profile
    mprofile = getnlines('MESA/profile'//trim(string(model))//'.data.GYRE')
    allocate(MESAprofile(19,mprofile+5))
    open(301, action='read',file = 'MESA/profile'//trim(string(model))//'.data.GYRE')
    read(301,*) dummy
    read(301,*) MESAprofile
    close(301)

    !! set general variables of the model
    call set_variables(MESAhistory,modelini)

    !! obtain radius of radiative region from MESA profile N2
    call getradiativeR(ri,rf,MESAprofile)

    !! run GYRE
    if (GYRE_var == 0) then
      call runGYRE(model)
    end if 

    !! compute mixed modes flux
    allocate(F_total(3,size(MESAprofile,dim=2)))
    F_total = 0.
    do i=1,size(MESAprofile,dim=2)
      F_total(:,i) = i
    end do
    !call compute_flux(model,MESAprofile,F_total)
    
    !! compute new rotation profile
    mnew = 5+getnlines('MESA/profile'//trim(string(model+dmodel))//'.data.GYRE')
    allocate(omega(mnew))
    !call rot_profile_conserveAM(MESAprofile,omega,model+dmodel,mnew)
    !!call rotation_profile(F_total,MESAprofile,omega,model+dmodel,mnew)
    !call rotation_profile_implicit(F_total,MESAprofile,omega,model+dmodel,mnew)
    !call rotation_profile_explicit(MESAprofile,omega,model+dmodel,mnew)
    !call rotation_profile_CN_radiative(MESAprofile,omega,model+dmodel,mnew)
    !call rotation_profile_CN_nonunigrid(MESAprofile,omega,model+dmodel,mnew,F_total)
    call rotation_profile_Henyey(MESAprofile,omega,model+dmodel,mnew,F_total)


    !! insert new rotation profile in the next MESA profile
    call insertcol_MESAfile(model+dmodel,3,omega)
    
    !free the memory
    deallocate(MESAprofile,omega,F_total)
    
  end do

  deallocate(MESAhistory)

end program main
