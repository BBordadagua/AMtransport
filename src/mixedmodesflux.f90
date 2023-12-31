Module compute_mixedmodesflux
    use parameters
    use coefficient_equations
    use compute_modeamplitude
    use data_handling
    use math_functions
    !use meridionalcirculation_mod
    !use rotationprofile_mod
    !use numerical_scheme

  IMPLICIT NONE


  public :: compute_flux, sum_flux

  contains
  
  subroutine compute_flux(model,dmodel,MESAprofile,Flux)
    implicit none
    integer :: Nsum,status,count,nsize,i,ncolumns,jsmallest,nlines,lsmallest,ri,rf
    integer, allocatable :: jarray(:)
    integer, intent(in) :: model,dmodel
    real (DP), allocatable :: summaryfile(:,:), smalldetailfile(:,:),F_total(:,:)
    real (DP), intent(inout) :: Flux(:,:)
    real (DP), intent(in) :: MESAprofile(:,:)


    !! allocate the GYRE summary file array
    Nsum = getnlines('input/model_'//trim(string(model))//'/summary.txt')
    allocate(summaryfile(14,Nsum), STAT = status)
    IF(status/=0) STOP
    call readfile('input/model_'//trim(string(model))//'/summary.txt',summaryfile)

    !! use all the files generated by GYRE except l=0
    count = 0
    do i=1,Nsum
      if (int(summaryfile(8,i)) == 0) then
        count = count + 1
      end if
    end do
    nsize = Nsum - count

    !! allocate array with the order of all GYRE files 
    allocate(jarray(nsize), STAT = status)
    IF(status/=0) STOP
    do i=1,nsize
      jarray(i) = int(summaryfile(7,i+count))
    end do

    !! get GYRE smallest file
    ncolumns = 19 !number of columns of GYRE detail file
    call getsmallestfile(summaryfile,jarray,jsmallest, nlines,model) ! determine smallest file
    
    !! create flux array and smallest detail file
    allocate(F_total(3,nlines),smalldetailfile(ncolumns,nlines),STAT=status)
    IF(status/=0) STOP
    lsmallest = int(summaryfile(8,jsmallest))
    call readfile('input/model_'//trim(string(model))//'/detail.l'//trim(string(lsmallest))//&
    &'.j'//trim(string(jsmallest))//'.txt',smalldetailfile)

    
    !! compute F_total
    F_total(2,:) = 0.
    call sum_flux(nsize,model,jarray,smalldetailfile,nlines,ncolumns,F_total,MESAprofile,summaryfile)

    call getradiativeR(ri,rf,MESAprofile)
    !F_total(1,:) = smalldetailfile(15,:)*R_star/MESAprofile(2,rf) !radius normalized by base of convection zone
    !F_total(3,:) = smalldetailfile(3,:) !mass

    !! save output of AM by mixed modes
    !open(300, action='write',file='output/total_flux_'//trim(string(model))//'.txt')
    !call printMatrix(F_total, 3, nlines,300)
    !close(300) 

    !F_total(1,:) = smalldetailfile(15,:)*R_star*MESAprofile(2,rf)
    F_total(1,:) = smalldetailfile(15,:)*R_star !! radius in cm
    F_total(3,:) = smalldetailfile(3,:) !mass in g

    !! interpolate mixed mode flux to have the same size as MESA profile
    Flux(1,:) = interpolate(F_total(3,:),F_total(1,:),size(F_total,dim=2),MESAprofile(3,:),size(MESAprofile,dim=2)) !! radius
    Flux(2,:) = interpolate(F_total(3,:),F_total(2,:),size(F_total,dim=2),MESAprofile(3,:),size(MESAprofile,dim=2)) !! jdot
    Flux(3,:) = MESAprofile(3,:) !! mass

    !save output
    open(300, action='write',file='mixed_modes/total_flux_'//trim(string(model+dmodel))//'.txt')
    call printMatrix(Flux, 3, size(MESAprofile,dim=2),300)
    close(300) 
    
    deallocate(F_total,smalldetailfile,jarray,summaryfile)

  end subroutine compute_flux



  subroutine sum_flux(nsize,model,jarray,smalldetailfile,nlines,ncolumns,F_total,MESAprofile,summaryfile)
    implicit none
    real (DP), intent(inout) :: F_total(:,:)
    real (DP), intent(in) :: smalldetailfile(:,:), MESAprofile(:,:),summaryfile(:,:)
    real (DP), allocatable :: detailfile(:,:), coeff(:)
    real (DP) :: amp_ml
    integer, intent(in) :: nsize,model,jarray(:),nlines,ncolumns
    integer :: lvalue,m,i,status

    !open(50,action='write',file='output/amp_'//trim(string(model))//'.txt')

    do i=1,nsize

      !! get detailfile
      lvalue = int(summaryfile(8,jarray(i)))
      m=getnlines('input/model_'//trim(string(model))//'/detail.l'//trim(string(lvalue))//'.j'//trim(string(jarray(i)))//'.txt')
      allocate(detailfile(ncolumns,m),coeff(m), STAT = status)
      IF(status/=0) STOP
      call readfile('input/model_'//trim(string(model))//'/detail.l'//trim(string(lvalue))//&
      &'.j'//trim(string(jarray(i)))//'.txt',detailfile)
        
      !! compute coefficient of the mixed modes
      coeff = 0.
      !! full equation
      !call calculate_gwaves(detailfile,m,summaryfile,jarray(i),coeff,MESAprofile)
      !call calculate_coefficients(detailfile,m,summaryfile,jarray(i),coeff,MESAprofile)

      !! simplified equation
      !call calculate_gwaves_simplified(detailfile,m,summaryfile,jarray(i),coeff,MESAprofile)
      call calculate_simplifiedeq(detailfile,summaryfile,jarray(i),coeff,MESAprofile)

      !! compute amplitude of the modes
      amp_ml = 0.
      call compute_amplitude(detailfile,summaryfile,jarray(i),amp_ml)

      !! sum flux over all modes
      coeff(1)=0.
      F_total(2,:) = F_total(2,:) + interpolate(detailfile(15,:),amp_ml*coeff,m,smalldetailfile(15,:),nlines)
        
      deallocate(detailfile,coeff)
    end do
    !close(50)

  end subroutine sum_flux


      

      !derive gwaves to give jdot
      !F_total(1,:) = smalldetailfile(15,:)*R_star
      !do i=1,nlines-1
      !  F_total(3,i) = - derive(F_total(2,i),F_total(2,i+1),F_total(1,i+1) - F_total(1,i))/(2.*F_total(1,i)**2) 
      !end do
      !F_total(2,:) = F_total(3,:)

      !integrate jdot to give gwaves
      !F_total(1,:) = smalldetailfile(15,:)*R_star
      !do i=2,nlines
      !  F_total(3,i) = -(F_total(2,i-1)*F_total(1,i-1)**2 + F_total(2,i)*F_total(1,i)**2)*(F_total(1,i) - F_total(1,i-1))
      !end do
      !F_total(2,:) = F_total(3,:)

      !Add jdot to the next shell
      !do i=4,nlines
      !  F_total(3,i) = F_total(2,i) + abs(F_total(2,i-1))
      !end do
      !do i=4,nlines
      !  F_total(2,i) = F_total(3,i)
      !end do

       


    end module compute_mixedmodesflux
