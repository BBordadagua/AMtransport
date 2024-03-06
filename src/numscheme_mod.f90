Module numerical_scheme
  use parameters
  use data_handling
  use math_functions
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: variables
  public :: CrankNicolson_new,create_matrix_diffusionCN_new,get_CNcoeff
  public :: difeq,difeq_total,constants,solvde,constants_tot,smoothing_func
  contains

  function smoothing_func(y,k1,k2) result(s)
    implicit none
    real (DP), intent(in) :: y(:)
    integer, intent(in) :: k1,k2
    real (DP) :: s(k2)
    integer :: k
    s=y
    do k=k1+3,k2-3
      s(k) = y(k-2) + 2d0*y(k-1) + 3d0*y(k) + 2d0*y(k+1) + y(k+2)
      s(k) = s(k)/9d0
    end do
  end function

  subroutine variables(k1,k2,Mn,Mo,delta_t,var)
    implicit none
    real (DP), intent(inout) :: var(:,:)
    real (DP), intent(in) :: Mn(:,:), delta_t, Mo(:,:)
    integer, intent(in) :: k1,k2
    integer :: k,k1o,k2o,m
    real (DP) :: nu1,nu2,dnu

    !! variables computed using derivatives !!only for interior points
    !var(33,k1) = 0d0 !check
    var(37,k1) = 0d0 !check

    do k=k1,k2
      var(34,k) = Mn(26,k) !!kappa oppacity
      var(35,k) = Mn(13,k) !! mu

      var(1,k) = Mn(2,k)**(2d0/3d0) !! nu 
      var(2,k) = Mn(2,k) !! mass (M_sun)
      var(3,k) = 10**Mn(3,k) !! radius (R_sun)
      var(4,k) = 10**Mn(5,k) !! rho
      var(5,k) = Mo(6,k) !! jmodes
      var(6,k) = 10**Mn(10,k) !! g surf
      var(7,k) = 1d7 !! Dv incomplete
      var(8,k) = 10**Mo(2,k) !! radius old
      var(9,k) = Mo(4,k) !! omega old
      var(10,k) = Mn(15,k)    !! cp
      var(11,k) = 10**Mn(4,k)    !! temperature
      var(12,k) = 10**Mn(6,k)    !! pressure
      var(13,k) = Mn(29,k)   !! luminosity (L_sun)
      var(14,k) = Mn(14,k)    !! nabla_ad
      var(15,k) = Mn(30,k)    !! nabla !!not sure!!
      var(17,k) = 1d0 !! phi !(dlog(var(4,k))-dlog(var(4,k-1)))/(dlog(var(35,k))-dlog(var(35,k-1)))   !constant pressure and T
      var(16,k) = Mn(35,k)*Mn(12,k)*R_sun/(var(6,k)*var(17,k)) !! nabla_mu !(dlog(var(35,k))-dlog(var(35,k-1)))/(dlog(var(12,k))-dlog(var(12,k-1)))  
      var(18,k) = var(11,k)*var(4,k)*var(10,k)*var(14,k)/var(12,k) !! delta (eq. 4.21 Kippenhahn)
      !var(18,k) = (var(12,k)+a_raddensity*(var(11,k)**4))/(var(12,k)-a_raddensity*(var(11,k)**4)/3d0)!! delta !! (4-3*beta)/beta , beta=pgas_div_p
      !var(19,k) = (2d0*var(3,k)*var(3,k)*var(9,k))/(3d0*var(6,k))    !! theta_old 
      !var(20,k) = y(1,1)    !! omega_0  !!should be omega actual !! change
      var(21,k) = 10**Mn(5,1)    !! rho_0
      var(22,k) = 1.5d-6 !! beta for Dh (Richard & Zahn 1999) !!accouting for the dividion by 10
      var(23,k) = Mn(11,k)    !! HT !!not sure!!
      var(24,k) = (4d0*a_raddensity*c_speed*var(11,k)**3)/(3d0*var(4,k)*var(4,k)*var(10,k)*var(34,k))    !! grand K thermal conductivity
      var(25,k) = Mn(20,k)    !! eps_nuc 
      !if (var(25,k)>0d0) 
      var(26,k) = Mn(22,k)    !! eps_t 
      !if (var(25,k)>0d0) 
      var(27,k) = var(17,k)*Mn(21,k)    !! eps_mu !phi_e*ro/eps_nuc*depsro Cesam
      var(28,k) = 3d0-Mn(28,k)/var(34,k)+var(18,k)*(1d0+Mn(27,k)/var(34,k))  !! chi_t 
      var(29,k) = -var(17,k)*(1d0+Mn(27,k)/var(34,k)) !! chi_mu
      var(30,k) = Mn(12,k)*R_sun !R_sun**2 *var(3,k)**2 *var(12,k)/(G*M_sun*var(4,k)*Mn(2,k))     !! Hp
      !var(31,k) = 0d0    !! lambda old
      var(32,k) = delta_t !! delta_t
      if (k/=k1) var(33,k) =  (dlog(var(13,k))-dlog(var(13,k-1)))&
        &/(dlog(var(2,k))-dlog(var(2,k-1)))!! dln L/dln m  !! according to page 146 cesam is correct
      var(33,k1) = var(33,k1+1) !! shouldn't it be zero?
        
      !if (k/=k1) var(36,k) = (y(3,k)-y(3,k-1))/((var(1,k)-var(1,k-1))) !! for Dh
      if (k/=k1) var(37,k) = (dlog(var(4,k))-dlog(var(4,k-1)))/((var(1,k)-var(1,k-1))) !! for Dh
      var(37,k1) = var(37,k1+1) !! not sure!!

      var(38,k) = Mn(35,k) !! brunt_N2_composition_term !N2_mu
      var(39,k) = Mn(34,k) !! brunt_N2_structure_term  !N2_T
      var(40,k) = Mn(33,k) !! Brunt N2

      !var(40,:) = 0d0!! constant AM at the convective zone of previous timestep
      !call getradiativeR(k1o,k2o,Mo,1) !!change
      k1o = k1
      k2o = k2
      m = size(Mn,dim=2)
      var(41,:) = 0.5d0*((10**Mn(3,m))**2 *(Mn(2,m)**(1d0/3d0)) +(10**Mn(3,k2))**2 *(Mn(2,k2)**(1d0/3d0)))*&
      & ((Mn(2,m)**(2d0/3d0))-(Mn(2,k2)**(2d0/3d0)))!! Itop
      var(42,:) = 0.5d0*(Mo(4,m)*(10**Mo(2,m))**2 *(Mn(2,m)**(1d0/3d0)) + Mo(4,k2o)*(10**Mo(2,k2o))**2 *(Mn(2,k2o)**(1d0/3d0)))*&
      & ((Mn(2,m)**(2d0/3d0))-(Mn(2,k2o)**(2d0/3d0)))!! Itop_t
      var(43,:) =  0.5d0*((10**Mn(3,k1o))**2 *(Mn(2,k1o)**(1d0/3d0)) +(10**Mn(3,1))**2 *(Mn(2,1)**(1d0/3d0)))*&
      & ((Mn(2,k1o)**(2d0/3d0))-(Mn(2,1)**(2d0/3d0)))!! Ibot
      var(44,:) = 0.5d0*(Mo(4,k1o)*(10**Mo(2,k1o))**2 *(Mn(2,k1o)**(1d0/3d0)) + Mo(4,1)*(10**Mo(2,1))**2 *(Mn(2,1)**(1d0/3d0)))*&
      & ((Mn(2,k1o)**(2d0/3d0))-(Mn(2,1)**(2d0/3d0)))!! Ibot_t

    end do

    !! compute integration boundary conditions I_top
    var(41,:) = 0d0
    var(42,:) = 0d0
    m = size(Mn,dim=2)
    do k=k2,m
      if (k == k2) then
        nu1 = Mn(2,k2)**(2d0/3d0)
        nu2 = Mn(2,k2+1)**(2d0/3d0)
      else if (k == m) then
        nu1 = Mn(2,m-1)**(2d0/3d0)
        nu2 = Mn(2,m)**(2d0/3d0)
      else
        nu1 = Mn(2,k-1)**(2d0/3d0)
        nu2 = Mn(2,k+1)**(2d0/3d0)
      end if
      dnu = 0.5d0*(nu2-nu1)*dsqrt(Mn(2,k))
      var(41,:) = var(41,:) + dnu*(10**Mn(3,k))*(10**Mn(3,k))
      var(42,:) = var(42,:) + dnu*(10**Mo(2,k))*(10**Mo(2,k))*Mo(4,k) 
    end do

  end subroutine variables

  !! --------------------------------------------------------------------------------------------------
  !! Relaxation scheme --------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------------------

  subroutine constants(k1,k2,cte,var_old,y)
    implicit none
    real (DP), intent(inout) :: cte(:,:)
    real (DP), intent(in) :: var_old(:,:),y(:,:)
    integer, intent(in) :: k1,k2
    integer :: k
    real (DP),allocatable :: nu12(:),nu32(:),delta_nu(:),var_k(:,:),var(:,:)
    real (DP),allocatable :: y1(:),y2(:),y3(:)!,y4(:),y5(:)

    allocate(nu12(size(var_old,dim=2)),nu32(size(var_old,dim=2)),&
    &delta_nu(size(var_old,dim=2)),var_k(size(var_old,dim=1),size(var_old,dim=2)),&
    &var(size(var_old,dim=1),size(var_old,dim=2)),&
    &y1(size(y,dim=2)),y2(size(y,dim=2)),y3(size(y,dim=2)))!,y4(size(y,dim=2)),y5(size(y,dim=2)))
    
    var_k = var_old
    var(:,k1) = var_old(:,k1)
    do k=k1+1,k2
      y1(k) = 0.5d0*(y(1,k)+y(1,k-1))
      y2(k) = 0.5d0*(y(2,k)+y(2,k-1))
      delta_nu(k) = var_k(1,k) - var_k(1,k-1)
      !var(:,k) = 0.5d0*(var_old(:,k)+var_old(:,k-1))
      var(:,k) = 0.5d0*(var_k(:,k)+var_k(:,k-1))
    end do
    delta_nu(k1) = 0d0
    nu12 = dsqrt(var(1,:))
    nu32 = var(1,:)**(3d0/2d0)
    cte(46,:) = 1d7

    do k=k1,k2
      cte(1,k) = delta_nu(k)*nu12(k)
      cte(2,k) = (cte(1,k)*var(3,k)*var(3,k))/var_k(32,k)  !C1 cesam
      cte(3,k) = (cte(1,k)*var(8,k)*var(8,k)*var(9,k)) /var_k(32,k)  !C2 cesam
      cte(4,k) = (64d0*PI*PI*(R_sun**4))/(9d0*M_sun*M_sun) !! diffusion !C18 cesam
      cte(5,k) = cte(4,k)*(cte(46,k)*(var_k(3,k)**6)*(var_k(4,k)**2))  !! diffusion !C5 cesam
      if (k/=k1) cte(6,k) = cte(4,k)*(cte(46,k-1)*(var_k(3,k-1)**6)*(var_k(4,k-1)**2)) !! diffusion !C6 cesam
      cte(7,k) = cte(1,k)*var_k(5,k)/(R_sun*R_sun*var_k(4,k)) !! mixed modes !C33 cesam
      !!cte(7,k) = cte(1,k)*var(5,k)/(R_sun*R_sun*var(4,k)) !! mixed modes !C33 cesam
      if (k/=k1) cte(8,k) = cte(1,k)*var_k(5,k-1)/(R_sun*R_sun*var_k(4,k-1)) !! mixed modes 
      
    end do

    
    deallocate(nu12,nu32,delta_nu,var_k,var,y1,y2,y3)!,y4,y5)
  end subroutine constants

  subroutine constants_tot(k1,k2,cte,var_old,y)
    implicit none
    real (DP), intent(inout) :: cte(:,:)
    real (DP), intent(in) :: y(:,:) !,var_old(:,:)
    real (DP), intent(inout) :: var_old(:,:) !!for debuging only
    integer, intent(in) :: k1,k2
    integer :: k
    real (DP),allocatable :: nu12(:),nu32(:),delta_nu(:),var_k(:,:),var(:,:)
    real (DP),allocatable :: y1(:),y2(:),y3(:)!,y4(:),y5(:)
    real (DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8,aux9,aux10,aux11,nu_kin

    allocate(nu12(size(var_old,dim=2)),nu32(size(var_old,dim=2)),&
    &delta_nu(size(var_old,dim=2)),var_k(size(var_old,dim=1),size(var_old,dim=2)),&
    &var(size(var_old,dim=1),size(var_old,dim=2)),&
    &y1(size(y,dim=2)),y2(size(y,dim=2)),y3(size(y,dim=2)))!,y4(size(y,dim=2)),y5(size(y,dim=2)))
    
    
    var_old(36,k1) = 0d0
    !var_k(36,k1+1:k2) = (y(3,k1+1:k2)-y(3,k1:k2-1))/((var_k(1,k1+1:k2)-var_k(1,k1:k2-1))) !! for Dh
    do k=k1+1,k2
      var_old(36,k) = (y(3,k)-y(3,k-1))/((var_old(1,k)-var_old(1,k-1))) !! for Dh
    end do

    var_k = var_old

    !var(:,k1) = var_old(:,k1) !!old version
    var(:,k1) = var_k(:,k1)
    y1(k1) = y(1,k1)
    y2(k1) = y(2,k1)
    y3(k1) = y(3,k1)
    do k=k1+1,k2
      y1(k) = 0.5d0*(y(1,k)+y(1,k-1))
      y2(k) = 0.5d0*(y(2,k)+y(2,k-1))
      y3(k) = 0.5d0*(y(3,k)+y(3,k-1))
    
      delta_nu(k) = var_k(1,k) - var_k(1,k-1)
      !var(:,k) = 0.5d0*(var_old(:,k)+var_old(:,k-1)) !!old version
      var(:,k) = 0.5d0*(var_k(:,k)+var_k(:,k-1))
    end do
    delta_nu(k1) = 0d0
    nu12 = dsqrt(var(1,:))
    nu32 = var(1,:)**(3d0/2d0)

    do k=k1,k2
      !! computing Dh(k)
      aux1 = (var(22,k)*8d0*PI*R_sun**6 /(9d0*M_sun))*var_k(4,k)*var_k(3,k)**6
      aux2 = var(36,k)/dsqrt(var_k(1,k))
      aux3 = var(37,k)/dsqrt(var_k(1,k))
      aux4 = var(22,k)*R_sun**4 *var_k(3,k)**4 /(3d0*dsqrt(var_k(1,k))*R_sun*var_k(3,k))

      cte(44,k) = aux1*y(1,k)*(aux2 + aux3*y(3,k) - (3d0/2d0)*y(2,k)*y(3,k)/y(1,k)) - aux4*y(1,k)*y(3,k) !!Dh^2(k)
      cte(44,k) = sqrt(abs(cte(44,k))) !!Dh(k)

      !! computing Dh(k) CESAM version
      aux1 = ((var(22,k)*8d0*PI*R_sun**6)/(9d0*M_sun*dsqrt(var_k(1,k))))*var_k(4,k)*var_k(3,k)**6 !!C1
      aux2 = ((var(22,k)*8d0*PI*R_sun**6)/(9d0*M_sun*dsqrt(var_k(1,k))))*var_k(4,k)*var_k(3,k)**6 *var(37,k) &
      & - (var(22,k)*R_sun**3 * var_k(3,k)**3) /3d0 !!C2
      aux3 = (var(22,k)*4d0*PI*R_sun**6 *var_k(4,k)*var_k(3,k)**6)/(3d0*M_sun*dsqrt(var_k(1,k))) !!C3
      cte(44,k) = aux1*y(1,k)*var(36,k) + aux2*y(1,k)*y(3,k) - aux3*y(2,k)*y(3,k)*dsqrt(var_k(1,k)) 
      cte(44,k) = sqrt(abs(cte(44,k)))

      !! computing Dh middle point
      aux5 = (var(22,k)*8d0*PI*R_sun**6 /(9d0*M_sun))*var(4,k)*var(3,k)**6
      aux6 = var_k(36,k)/nu12(k)
      aux7 = var_k(37,k)/nu12(k)
      aux8 = var(22,k)*R_sun**4 *var(3,k)**4 /(3d0*nu12(k)*R_sun*var(3,k))

      cte(45,k) = aux5*y1(k)*(aux6 + aux7*y3(k) - (3d0/2d0)*y2(k)*y3(k)/y1(k)) - aux8 *y1(k)*y3(k) !!Dh^2 middle point
      cte(45,k) = sqrt(abs(cte(45,k))) !!Dh middle point

      !! computing Dv
      aux9 = (1d0/6d0)*(64d0*PI*PI*R_sun**6 *var_k(3,k)**6 *var_k(4,k)**2)/(9d0*M_sun*M_sun)
      aux10 = var_k(38,k)+var_k(39,k)
      if (k/=k1) cte(46,k) = aux9*y(2,k)*y(2,k)*&
                & (cte(44,k))*((cte(44,k))+var_k(24,k))/((cte(44,k))*aux10 + (var_k(24,k)*var_k(38,k))) !!Dv
                

      !! using a power law
      !cte(44,k) = (var_k(3,k))**(3d0/2d0)*1d7
      !cte(45,k) = (var(3,k))**(3d0/2d0)*1d7
      !cte(46,k) = (var(3,k))**(3d0/2d0)*1d4

      !cte(44,k) = 5d4 * var_k(3,k)*R_sun*abs(y(3,k)) !Formalisme simplified Castro, Vauclair & Richard
      !cte(45,k) = 5d4 * var(3,k)*R_sun*abs(y3(k))
      !cte(46,k) = (6d0-1d0/(30d0*5d4) )* var_k(3,k)*R_sun*abs(y(3,k))

      !write(607,*) var(3,k),cte(45,k),cte(44,k),cte(46,k)
      !cte(46,k) = 1d-7!1d7 !!Dv

      nu_kin = 7d0
      if (cte(44,k) <= nu_kin .and. (var_k(3,k) < 0.5d0) ) then
        cte(44,k) = nu_kin
        cte(45,k) = nu_kin
        cte(46,k) = nu_kin
      else
        !! computing Dv(k) CESAM version
        aux9 = (1d0/6d0)*(64d0*PI*PI*R_sun**6 *var_k(3,k)**6 *var_k(4,k)**2)/(9d0*M_sun*M_sun)
        aux10 = 1d0 + var_k(24,k)/cte(44,k)
        aux11 = var_k(39,k) + var(38,k)*aux10
        cte(46,k) = aux9*cte(44,k)*aux10*y(2,k)*y(2,k)/aux11
        

        if (cte(46,k) > nu_kin) then
          cte(46,k) = cte(46,k) + nu_kin
        else
          cte(46,k) = nu_kin
        end if
      end if
      
    end do

    !! smoothing Dh Dv
    !var(44,:) = smoothing_func(var(44,:),k1,k2)
    !var(46,:) = smoothing_func(var(46,:),k1,k2)

    do k=k1,k2 
      cte(1,k) = delta_nu(k)*nu12(k)
      cte(2,k) = (cte(1,k)*var(3,k)*var(3,k))/var_k(32,k)  !C1 cesam
      cte(3,k) = (cte(1,k)*var(8,k)*var(8,k)*var(9,k)) /var_k(32,k)  !C2 cesam
      cte(4,k) = (64d0*PI*PI*(R_sun**4))/(9d0*M_sun*M_sun) !! diffusion !C18 cesam
      cte(5,k) = cte(4,k)*(cte(46,k)*(var_k(3,k)**6)*(var_k(4,k)**2))  !! diffusion !C5 cesam
      if (k/=k1) cte(6,k) = cte(4,k)*(cte(46,k-1)*(var_k(3,k-1)**6)*(var_k(4,k-1)**2)) !! diffusion !C6 cesam
      cte(7,k) = cte(1,k)*var_k(5,k)/(R_sun*R_sun*var_k(4,k)) !! mixed modes !C33 cesam
      !!cte(7,k) = cte(1,k)*var(5,k)/(R_sun*R_sun*var(4,k)) !! mixed modes !C33 cesam
      if (k/=k1) cte(8,k) = cte(1,k)*var_k(5,k-1)/(R_sun*R_sun*var_k(4,k-1)) !! mixed modes 
      cte(9,k) = ((8d0*PI*R_sun*R_sun)/(15d0*M_sun))*var_k(4,k)*var_k(3,k)**4 !C3 cesam
      if (k/=k1) cte(10,k) = ((8d0*PI*R_sun*R_sun)/(15d0*M_sun))*var_k(4,k-1)*var_k(3,k-1)**4 !C4 cesam
      !cte(11,k) = (M_sun/L_sun)*(nu32(k)*var(6,k)*var(10,k)*var(4,k)*var(11,k)/(var(12,k)*var(13,k)))*&
      !        &(var(14,k)-var(15,k)+var(17,k)*var(16,k)/var(18,k)) !C7 cesam
      cte(11,k) = (M_sun/L_sun)*(nu32(k)*var(6,k)*var(10,k)*var(4,k)*var(11,k)/(var(12,k)*var(13,k)))*&
              &abs(var(14,k)-var(15,k)+var(17,k)*var(16,k)/var(18,k)) !C7 cesam !modified
      cte(12,k) = -(M_sun/(L_sun*var(32,k)))*(nu32(k)*var(10,k)*var(11,k)*var(9,k)*var(19,k)/(var(13,k)*var(18,k))) !C8 cesam
      cte(13,k) = 8d0*R_sun*var(3,k)/(3d0*var(6,k)) !C9 cesam
      cte(14,k) = var(33,k) - var(20,k)*var(20,k)/(2d0*PI*G*var(21,k)) !C10 cesam
      cte(15,k) = 1d0/(2d0*PI*G*var(4,k)) !C11 cesam
      if (k/=k1) cte(16,k) = (2d0*var(1,k))/(3d0*delta_nu(k)) !C12 cesam
      cte(17,k) = (3d0*M_sun/(2d0*PI*R_sun**4))*(nu32(k)*cte(45,k)*var(23,k)/(var(4,k)*var(24,k)*var(18,k)*var(3,k)**4))!C1_star cesam
      cte(18,k) = (M_sun/(2d0*PI*R_sun**3))*(nu32(k)/(var(4,k)*var(3,k)**3)) !C2_star cesam
      cte(19,k) = (M_sun/L_sun)*(var(25,k)*nu32(k)/var(13,k)) !C3_star cesam
      cte(20,k) = var(33,k) !C4_star cesam !C15 cesam
      cte(21,k) = (3d0*M_sun/(2d0*PI*R_sun**4))*(nu32(k)*var(23,k)/(var(4,k)*var(3,k)**4)) !C5_star cesam
      cte(22,k) = (M_sun/L_sun)*(var(25,k)*var(26,k)*nu32(k)/var(13,k)) !C6_star cesam
      cte(23,k) = cte(17,k)-cte(18,k)+cte(19,k)-cte(20,k) !C13 cesam
      cte(24,k) = (M_sun/(L_sun*var(32,k)))*(nu32(k)*var(10,k)*var(11,k)/(var(13,k)*var(18,k))) !C14 cesam
      cte(25,k) = (M_sun/L_sun)*(var(25,k)*var(27,k)*nu32(k)/var(13,k)) !C16 cesam
      cte(26,k) = cte(21,k) - cte(22,k) ! C17 cesam
      cte(27,k) = -(8d0*PI*R_sun*R_sun/(3d0*M_sun))*(var(3,k)*var(3,k)*var(4,k)*var(23,k)) !C19 cesam
      cte(28,k) = cte(1,k)*(1d0+var(28,k)) !C20 cesam
      cte(29,k) = cte(1,k)*var(29,k) !C21 cesam
      cte(30,k) = var_k(17,k)/var_k(18,k) !C22 cesam
      if (k/=k1) cte(31,k) = var_k(17,k-1)/var_k(18,k-1) !C23 cesam
      cte(32,k) = var(17,k)/var(18,k) !C24 cesam
      cte(33,k) = 1d0/var_k(18,k) !C25 cesam
      if (k/=k1) cte(34,k) = 1d0/var_k(18,k-1) !C26 cesam
      cte(35,k) = 1d0/var(18,k) !C27 cesam
      if (k/=k1) cte(36,k) = var_k(16,k-1)/var_k(30,k-1) !C28 cesam
      cte(37,k) = 1d0/var(32,k) !C29 cesam
      if (k/=k1) cte(38,k) = (6d0/(R_sun*R_sun))*((cte(44,k-1))/(var_k(3,k-1)*var_k(3,k-1))) !C30 cesam
      if (k/=k1) cte(39,k) = var(31,k-1)/var(32,k-1) !C31 cesam
      cte(40,k) = (2d0*R_sun**3/(3d0*G*M_sun))*(var(3,k)**3/nu32(k)) !C32 cesam
      cte(41,k) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*((var(4,k)*var(3,k)**4)/var(6,k)) 
      cte(42,k) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*((var_k(4,k)*var_k(3,k)**4)/var_k(6,k))
      if (k/=k1) cte(43,k) = ((16d0*PI*R_sun**4)/(9d0*M_sun))*((var_k(4,k-1)*var_k(3,k-1)**4)/var_k(6,k-1))
      
    end do

    
    !! constants for boundary conditions !! solid body rotation in convective zones
    cte(47,:) = var_k(41,:) !! Itop cesam
    cte(48,:) = var_k(42,:) !! Itop_t cesam !!change
    cte(49,:) = 1d0/var_k(32,:) !! K1 cesam
    cte(50,:) = (8d0*PI*R_sun*R_sun/(15d0*M_sun))*var(4,:)*var(3,:)**4 !! K2 cesam
    cte(51,:) = 0d0 !! K3 cesam
    cte(52,:) = 0d0 !! K4 cesam
    cte(53,:) = 0.5d0*(var(5,k2)*y1(k2)*nu12(k2)/(R_sun*R_sun*var(4,k2)) &
    & + var(5,k1)*y1(k1)*nu12(k1)/(R_sun*R_sun*var(4,k1)))*(var(1,k2)-var(1,k1)) !! K5 jmodes omega int
    cte(54,:) = 0.5d0*(var(5,k2)*nu12(k2)/(R_sun*R_sun*var(4,k2)) &
    & + var(5,k1)*nu12(k1)/(R_sun*R_sun*var(4,k1)))*(var(1,k2)-var(1,k1)) !! K6 jmodes_tot
    cte(55,:) = var_k(43,:)!! Ibot cesam
    cte(56,:) = var_k(44,:)!! I_bot_t cesam


    deallocate(nu12,nu32,delta_nu,var_k,var,y1,y2,y3)!,y4,y5)
  end subroutine constants_tot


  !! Rewrite for the specific problem !!
  !The only information returned from difeq to solvde is the matrix of derivatives 
  !s(i,j); all other arguments are input to difeq and should not be altered. 
  !k indicates the current mesh point, or block number. k1,k2 label the first and 
  !last point in the mesh. If k=k1 or k>k2, the block involves the boundary conditions
  !at the first or final points; otherwise the block acts on FDEs coupling variables
  !at points k-1, k.
  !Returns matrix s(i,j) for solvde. 
  SUBROUTINE difeq(k,k1,k2,jsf,indexv,s,y,cte)
    IMPLICIT NONE 
    REAL (DP), intent(in) :: cte(:,:),y(:,:)
    INTEGER, intent(in) :: jsf,k,k1,k2,indexv(:)
    REAL (DP), intent(inout) :: s(:,:)
    
    !s=0d0 !! added makes matrix singular
    
    if(k.eq.k1) then !Boundary condition at first point. 
      s(1,2+indexv(1))=0d0 
      s(1,2+indexv(2))=1d0
      s(1,jsf)=y(2,k1)

      s(2,2+indexv(1))=0d0
      s(2,2+indexv(2))=1d0
      s(2,jsf)=y(2,k1)
    
      !call printMatrix(s,2,5,5000)
    else if(k.gt.k2) then !Boundary conditions at last point. 
      s(1,2+indexv(1))=0d0 
      s(1,2+indexv(2))=1d0
      s(1,jsf)=y(2,k2)

      s(2,2+indexv(1))=0d0 
      s(2,2+indexv(2))=1d0
      s(2,jsf)=y(2,k2)
  
      !call printMatrix(s,2,5,2000)
    else !Interior point. 
      
      s(1,indexv(1))=0.5d0*(cte(2,k))!-cte(8,k))   !E1
      s(1,indexv(2))=cte(6,k)
      s(1,2+indexv(1))= 0.5d0*(cte(2,k))!-cte(7,k))
      s(1,2+indexv(2))=-cte(5,k)

      s(2,indexv(1))=1d0     !E2
      s(2,indexv(2))=0.5d0*cte(1,k)
      s(2,2+indexv(1))=-1d0
      s(2,2+indexv(2))=0.5d0*cte(1,k)

      !! expression with rdot linear interpolation
      s(1,jsf)=0.5d0*cte(2,k)*(y(1,k)+y(1,k-1))-cte(3,k)&
              & -cte(5,k)*y(2,k)+cte(6,k)*y(2,k-1) !& !! diffusion
              !& -0.5d0*cte(7,k)*y(1,k) - 0.5d0*cte(8,k)*y(1,k-1) !!mixed modes
      s(2,jsf)=0.5d0*cte(1,k)*(y(2,k)+y(2,k-1))-(y(1,k)-y(1,k-1)) !E2

      !call printMatrix(s,2,5,3000)
    endif 

    return 

  END SUBROUTINE difeq

  !! prescription of Talon & Zahn 1997 (Cesam2k20)
  SUBROUTINE difeq_total(k,k1,k2,jsf,indexv,s,y,cte)
    IMPLICIT NONE 
    REAL (DP), intent(in) :: cte(:,:),y(:,:)
    INTEGER, intent(in) :: jsf,k,k1,k2,indexv(:)
    REAL (DP), intent(inout) :: s(:,:)
    real (DP) :: y1,y2,y3,y4,y5,psi,psik,psik_1
    real (DP) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    
    !s=0d0 !! added makes matrix singular
    
    if(k.eq.k1) then !Boundary condition at first point. 
      !! considering only radiative center
      s(4,5+indexv(1))=0d0 
      s(4,5+indexv(2))=1d0  !! should be a factor constant??
      s(4,5+indexv(3))=0d0 
      s(4,5+indexv(4))=0d0 
      s(4,5+indexv(5))=0d0 
      s(4,jsf)=y(2,k1) 

      s(5,5+indexv(1))=0d0 
      s(5,5+indexv(2))=0d0 
      s(5,5+indexv(3))=1d0
      s(5,5+indexv(4))=0d0 
      s(5,5+indexv(5))=0d0 
      s(5,jsf)=y(3,k1)
    
       !call printMatrix(s,5,11,7000)
    else if(k.gt.k2) then !Boundary conditions at last point. 
      s(1,5+indexv(1))=cte(49,k2)*cte(47,k2) + cte(50,k2)*y(3,k2)-cte(52,k2)!-cte(54,k2)!mixedmodes
      !s(1,5+indexv(1))=cte(49,k)*cte(47,k2) + cte(50,k)*y(3,k2)-cte(52,k2)!-cte(54,k2)!mixedmodes
      s(1,5+indexv(2))=0d0 
      s(1,5+indexv(3))=cte(50,k2)*y(1,k2) 
      s(1,5+indexv(4))=0d0 
      s(1,5+indexv(5))=0d0 
      s(1,jsf)=cte(49,k2)*(cte(47,k2)*y(1,k2)-cte(48,k2))+cte(50,k2)*y(1,k2)*y(3,k2)-cte(51,k2)!+cte(53,k2)!mixedmodes
      !s(1,jsf)=cte(49,k)*(cte(47,k2)*y(1,k2)-cte(48,k2))+cte(50,k)*y(1,k2)*y(3,k2)-cte(51,k2)!+cte(53,k2)!mixedmodes

      s(2,5+indexv(1))=0d0 
      s(2,5+indexv(2))=1d0 !! should be a factor constant??
      s(2,5+indexv(3))=0d0 
      s(2,5+indexv(4))=0d0 
      s(2,5+indexv(5))=0d0 
      s(2,jsf)=y(2,k2) 

      s(3,5+indexv(1))=0d0 
      s(3,5+indexv(2))=0d0 
      s(3,5+indexv(3))=0d0 
      s(3,5+indexv(4))=0d0 
      s(3,5+indexv(5))=1d0 
      s(3,jsf)=y(5,k2)
  
      !call printMatrix(s,5,11,11000)
    else !! Interior points

      y1 = 0.5d0*(y(1,k)+y(1,k-1))
      y2 = 0.5d0*(y(2,k)+y(2,k-1))
      y3 = 0.5d0*(y(3,k)+y(3,k-1))
      y4 = 0.5d0*(y(4,k)+y(4,k-1))
      y5 = 0.5d0*(y(5,k)+y(5,k-1))
      psi = cte(35,k)*y1*cte(41,k)*y2 - cte(32,k)*y5
      psik = cte(33,k)*y(1,k)*cte(42,k)*y(2,k) - cte(30,k)*y(5,k)
      psik_1 = cte(34,k)*y(1,k-1)*cte(43,k)*y(2,k-1) - cte(31,k)*y(5,k-1)
      
      !! Omega y(1,:) ------------------------------------------------------------------------------
      !! derivatives dE/dy(k-1)
      s(1,indexv(1))=0.5d0*cte(2,k)+cte(10,k)*y(3,k-1) !-0.5d0*cte(8,k) !mixedmodes
      s(1,indexv(2))=cte(6,k)
      s(1,indexv(3))=cte(10,k)*y(1,k-1)
      s(1,indexv(4))=0d0
      s(1,indexv(5))=0d0
      !! derivatives dE/dy(k)
      s(1,5+indexv(1))=0.5d0*cte(2,k)-cte(9,k)*y(3,k) !-0.5d0*cte(7,k) !mixedmodes
      s(1,5+indexv(2))=-cte(5,k)
      s(1,5+indexv(3))=-cte(9,k)*y(1,k)
      s(1,5+indexv(4))=0d0
      s(1,5+indexv(5))=0d0
      !! equation E1
      s(1,jsf) = cte(2,k)*y1 - cte(3,k)&
        &-cte(9,k)*y(1,k)*y(3,k) + cte(10,k)*y(1,k-1)*y(3,k-1)& !! meridional circ
        &-cte(5,k)*y(2,k) + cte(6,k)*y(2,k-1)!&                  !! diffusion
        !& -0.5d0*cte(7,k)*y(1,k) - 0.5d0*cte(8,k)*y(1,k-1)      !! mixed modes
        !&-cte(7,k)*y1                                           !! mixed modes
      
      !! Theta = domega/dr y(2,:) ------------------------------------------------------------
      !! derivatives dE/dy(k-1)
      s(2,indexv(1))=1d0
      s(2,indexv(2))=0.5d0*cte(1,k)
      s(2,indexv(3))=0d0
      s(2,indexv(4))=0d0
      s(2,indexv(5))=0d0
      !! derivatives dE/dy(k)
      s(2,5+indexv(1))=-1d0
      s(2,5+indexv(2))=0.5d0*cte(1,k)
      s(2,5+indexv(3))=0d0
      s(2,5+indexv(4))=0d0
      s(2,5+indexv(5))=0d0
      !! equation E2
      s(2,jsf) = cte(1,k)*y2 - (y(1,k)-y(1,k-1))

      !! U2 y(3,:) ------------------------------------------------------------------
      aux1 = 1d0-cte(40,k)*y1*y1
      aux2 = cte(23,k)+cte(24,k)*aux1
      aux3 = -0.5d0*(aux2+cte(26,k)*cte(35,k))*cte(41,k)
      aux4 = cte(13,k)*(-1d0+cte(14,k)+2d0*cte(15,k)*y1*y1)
      !! derivatives dE/dy(k-1)
      s(3,indexv(1))= aux4*y1 + aux3*y2 &
                    & - cte(40,k)*y1*(cte(11,k)*y3-cte(12,k)-cte(24,k)*cte(41,k)*y1*y2)  
      s(3,indexv(2))=aux3*y1
      s(3,indexv(3))=0.5d0*cte(11,k)*aux1
      s(3,indexv(4))=cte(16,k) - 0.5d0*cte(20,k)
      s(3,indexv(5))=0.5d0*(cte(26,k)*cte(32,k)-cte(25,k))
      !! derivatives dE/dy(k)
      s(3,5+indexv(1))=s(3,indexv(1))
      s(3,5+indexv(2))=s(3,indexv(2))
      s(3,5+indexv(3))=s(3,indexv(3))
      s(3,5+indexv(4))=s(3,indexv(4)) - 2d0*cte(16,k)
      s(3,5+indexv(5))=s(3,indexv(5))
      !! equation E3
      s(3,jsf) = cte(11,k)*aux1*y3 &
        & - cte(12,k)*aux1 &
        & + aux4*y1*y1 - cte(13,k)*cte(15,k)*y1**4 &
        & - cte(16,k)*(y(4,k)-y(4,k-1)) &
        & - aux2*cte(41,k)*y1*y2 &
        & - cte(20,k)*y4 &
        & - cte(25,k)*y5 &
        & - cte(26,k)*psi

      !! Upsilon y(4,:) ----------------------------------------------------------------
      aux5 = cte(27,k)*cte(34,k)*cte(43,k)
      aux6 = 0.5d0*(cte(1,k)-cte(28,k)*cte(35,k))*cte(41,k)
      aux7 = 0.5d0*(cte(29,k)+cte(28,k)*cte(32,k))
      aux8 = cte(27,k)*cte(33,k)*cte(42,k)
      !! derivatives dE/dy(k-1)
      s(4,indexv(1))=aux5*y(2,k-1)+aux6*y2
      s(4,indexv(2))=aux5*y(1,k-1)+aux6*y1
      s(4,indexv(3))=0d0
      s(4,indexv(4))=0.5d0*cte(1,k)
      s(4,indexv(5))=aux7 - cte(27,k)*cte(31,k)
      !! derivatives dE/dy(k)
      s(4,5+indexv(1))=-aux8*y(2,k)+aux6*y2
      s(4,5+indexv(2))=-aux8*y(1,k)+aux6*y1
      s(4,5+indexv(3))=0d0
      s(4,5+indexv(4))=0.5d0*cte(1,k)
      s(4,5+indexv(5))=aux7 + cte(27,k)*cte(30,k)
      !! equation E4
      s(4,jsf) = cte(1,k)*y4 &
        & - cte(27,k)*(psik-psik_1) &
        & - cte(28,k)*psi &
        & + cte(29,k)*y5 &
        & + cte(1,k)*y1*cte(41,k)*y2

      !! Lambda y(5,:) ---------------------------------------------------------------------
      !! derivatives dE/dy(k-1)
      s(5,indexv(1))=0d0 
      s(5,indexv(2))=0d0
      s(5,indexv(3))=-cte(36,k)
      s(5,indexv(4))=0d0
      s(5,indexv(5))=cte(37,k)+cte(38,k)
      !! derivatives dE/dy(k)
      s(5,5+indexv(1))=0d0
      s(5,5+indexv(2))=0d0
      s(5,5+indexv(3))=0d0
      s(5,5+indexv(4))=0d0
      s(5,5+indexv(5))=0d0
      !! equation E5
      s(5,jsf) = (cte(37,k)+cte(38,k))*y(5,k-1) &
        & - cte(39,k) &
        & - cte(36,k)*y(3,k-1)

      !call printMatrix(s,5,11,9000)
    endif 

    return 

  END SUBROUTINE difeq_total

  !USES bksub,difeq,pinvs,red
  !Driver routine for solution of two point boundary value problems by relaxation. 
  !itmax is the maximum number of iterations. conv is the convergence criterion 
  !(see text). slowc controls the fraction of corrections actually used after each 
  !iteration. scalv(1:nyj) contains typical sizes for each dependent variable, used 
  !to weight errors. indexv(1:nyj) lists the column ordering of variables used to 
  !construct the matrix s of derivatives. (The nb boundary conditions at the first 
  !mesh point must contain some dependence on the first nb variables listed in indexv.)
  !The problem involves ne equations for ne adjustable dependent variables at each 
  !point. At the first mesh point there are nb boundary conditions. There are a total
  !of m mesh points. y(1:nyj,1:nyk) is the two-dimensional array that contains the
  !initial guess for all the dependent variables at each mesh point. On each 
  !iteration, it is updated by the calculated correction. The arrays 
  !c(1:nci,1:ncj,1:nck), s(1:nsi,1:nsj) supply dummy storage used by the relaxation 
  !code; the minimum dimensions must satisfy: nci=ne, ncj=ne-nb+1, nck=m+1, nsi=ne, nsj=2*ne+1.
  SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m, &
    & y,c,s,cte,var,k1,k2,value)
    implicit none 
    !real (DP), intent(in) :: var(:,:) !for debuging only
    real (DP), intent(in) :: conv,slowc !! changed
    real (DP), intent(inout) :: c(:,:,:),s(:,:),y(:,:),scalv(:),cte(:,:),var(:,:) !for debuging only
    integer, intent(in) :: itmax,m,ne,nb,indexv(:),k1,k2,value
    integer :: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8, &
    & j9,jc1,jcf,jv,k,km,kp,nvars,kmax(ne)
    real (DP) :: err,errj,fac,vmax,vz,ermax(ne)
    integer :: i

    s=0d0
    c=0d0

    !k1=1 !Set up row and column markers. 
    !k2=m 
    nvars=ne*m 
    j1=1 
    j2=nb 
    j3=nb+1 
    j4=ne 
    j5=j4+j1 
    j6=j4+j2 
    j7=j4+j3 
    j8=j4+j4 
    j9=j8+j1 
    ic1=1 
    ic2=ne-nb 
    ic3=ic2+1 
    ic4=ne 
    jc1=1 
    jcf=ic3

    if (value == 0) then
      call constants(k1,k2,cte,var,y)
    else 
      call constants_tot(k1,k2,cte,var,y)
    end if
    open(15, action='write',file='var2.txt')
    open(16, action='write',file='cte2.txt')
    do i=k1,k2
      write(15,*) var(:,i)
      write(16,*) cte(:,i)
    end do
    close(15)
    close(16)

    
    do it=1,itmax 
      !if (value == 0) then
      !  call constants(k1,k2,cte,var,y)
      !else 
      !  call constants_tot(k1,k2,cte,var,y)
      !end if
      
      open(1004, action='write',file='output/it'//trim(string(it))//'.txt')
      do i=k1,k2
        if (value == 0) write(1004,*) var(3,i),y(1,i),y(2,i)
        if (value /= 0) write(1004,*) var(3,i),y(1,i),y(2,i),y(3,i),y(4,i),y(5,i),cte(44,i),cte(45,i),cte(46,i)
      end do
      close(1004)

      !Primary iteration loop. 
      k=k1 !Boundary conditions at first point. 
      if (value == 0) then
        call difeq(k,k1,k2,j9,indexv,s,y,cte) 
      else 
        call difeq_total(k,k1,k2,j9,indexv,s,y,cte) 
      end if
      call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s) 
      
      do k=k1+1,k2 
        !Finite difference equations at all point pairs. 
        kp=k-1 
        if (value == 0) then
          call difeq(k,k1,k2,j9,indexv,s,y,cte) 
        else
          call difeq_total(k,k1,k2,j9,indexv,s,y,cte) 
        end if
        call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s) 
        call pinvs(ic1,ic4,j3,j9,jc1,k,c,s) 
      enddo 

      k=k2+1 !Final boundary conditions.
      if (value == 0) then 
        call difeq(k,k1,k2,j9,indexv,s,y,cte) 
      else
        call difeq_total(k,k1,k2,j9,indexv,s,y,cte)
      end if
      call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s) 
      call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
      call bksub(ne,nb,jcf,k1,k2,c) !Backsubstitution. 
      err=0d0

      do j=1,ne !Convergence check, accumulate average error. 

      !! scalv(1:nyj) contains typical sizes for each dependent variable, used to weight errors.
        !! cesam2k20 version
        !scalv = maxval(dabs(y(1:nrot_eq,:)),2) !! Cesam2k20
        scalv(j) = maxval(dabs(y(j,:))) !! each variable has different scaling
        if (scalv(j) < 1.0d-20) scalv = 1.0d0 !! Cesam2k20

        jv=indexv(j) 
        errj=0d0
        km=0 !! this can't be zero 
        vmax=0d0 
        do k=k1,k2 !Find point with largest error, for each dependent variable. 
          vz=dabs(c(jv,1,k)) 
          if(vz.gt.vmax) then 
            vmax=vz 
            km=k 
          endif 
          errj=errj+vz 
        enddo 
        !err=err+errj/scalv(j) !Note weighting for each dependent variable. 
        !ermax(j)=c(jv,1,km)/scalv(j) 
        ermax(j)=dabs(c(jv,1,km))/scalv(j) !! Cesam2k20!
        kmax(j)=km 
      enddo 
      
      err = maxval(ermax) !! Cesam2k20
      
      !err=err/nvars 
      !fac=slowc/max(slowc,err) !Reduce correction applied when error is large.
      if(err > 1.0d2) then !! cesam2k20
        fac = 10.0d0*slowc/max(slowc,err)
      else
        fac = slowc/max(slowc,err)
      endif
      
      do j=1,ne !Apply corrections. 
        jv=indexv(j) 
        do k=k1,k2 
          !if (jv == 1) write(9000,*) y(2,k),c(1,1,k), fac, err, scalv(j)
          !if (jv == 2) write(8000,*) y(1,k),c(2,1,k), fac, err, scalv(j)
          y(j,k)=y(j,k)-fac*c(jv,1,k) 
          
        enddo
        !if (jv == 1) write(9000,*) '  '
        !if (jv == 2) write(8000,*) '  '
      enddo 
      write(*,100) 'it',it,'   err',err,'   fac',fac
      !! Summary of corrections for this step. Point with largest error for each 
      !! variable can be monitored by writing out kmax and ermax. 
      print*, 'jv ', 1, 'kmax ',kmax(1), 'ermax ', ermax(1)
      print*, 'jv ', 2, 'kmax ',kmax(2), 'ermax ', ermax(2)
      if (value /= 0) print*, 'jv ', 3, 'kmax ',kmax(3), 'ermax ', ermax(3)
      if (value /= 0) print*, 'jv ', 4, 'kmax ',kmax(4), 'ermax ', ermax(4)
      if (value /= 0) print*, 'jv ', 5, 'kmax ',kmax(5), 'ermax ', ermax(5)
      if(err.lt.conv) return 
    enddo
    !100 format(1x,i4,2f12.6)
    100 format(A3,i4,A6,ES14.4,A6,ES14.3)


    print*, 'itmax exceeded in solvde' !Convergence failed. 
    stop

    
    return

  END SUBROUTINE solvde


!! --------------------------------------------------------------------------------------------------
!! Crank-Nicolson scheme ----------------------------------------------------------------------------
!! --------------------------------------------------------------------------------------------------

  function CrankNicolson_new(Mn,m,iter,delta_t,Mo) result(omega)
    implicit none
    real (DP), intent(in) :: delta_t
    integer, intent(in) :: iter,m
    real (DP), intent(inout) :: Mn(:,:),Mo(:,:)
    real (DP), allocatable :: U(:),V(:),mA_triag(:,:),mB_triag(:,:),dr(:),old_radius(:),omega_old(:)
    real (DP) :: omega(m)
    integer :: j,k1,k2
    real (DP), allocatable :: var(:,:)

    !! defining the begining and ending point of integration
    call getradiativeR(k1,k2,Mn,1)
    !k1 = 1 !ri
    !k2 = m !rf

    !!allocate matrix and vector
    allocate(mA_triag(3,k2),mB_triag(3,k2))
    allocate(U(k2),V(k2),dr(k2),old_radius(k2),omega_old(k2))

    !! compute constants
    allocate(var(50,k2)) !! number of constants
    call variables(k1,k2,Mn,Mo,delta_t,var)

    dr = (var(3,:) - var(8,:))/iter !! divide the radius into smaller steps according to smaller timesteps
    old_radius = var(8,:)
    var(3,:) = old_radius + dr !! radius current model

    do j=1,k2
      omega_old(j) = var(9,j)
    end do
    V = omega_old !omega old model
    
    
    do j=1,iter    !!subtimesteps
      if (j/=1) then
        var(8,:) = var(3,:)
        var(3,:) = old_radius + j*dr
      end if
      
      !! Create matrix
      call create_matrix_diffusionCN_new(k2,mA_triag,mB_triag,var)

      !! multiplying matrix
      U = 0d0
      !U = matmul(transpose(matrix_B),V)
      U = multiply_matrix_vector_triag(mB_triag,V,k2)

      !!invert Tridiagonal matrix
      omega_old = 0d0
      call invert_tridagmatrix(mA_triag,k2,omega_old,U)

      V = omega_old
    end do

    do j=k1,m
      if (j>k2) then
        omega(j) = omega_old(k2)
      else 
        omega(j) = omega_old(j)
      end if
    end do

    deallocate(U,V,mA_triag,mB_triag,dr,old_radius,var,omega_old)

  end function CrankNicolson_new

  subroutine get_CNcoeff(var,lambda,beta,alpha,epsilon,C1,a1,a2)
    implicit none
    real (DP), intent(in) :: var(:,:)
    real (DP) :: rho,radius,delta_t,nu_diff
    integer :: i,m
    real (DP), intent(inout) :: lambda(:),beta(:),alpha(:),epsilon(:),C1(:)
    real (DP), intent(inout) :: a1(:), a2(:)

    m = size(var,dim=2)
    
    do i=1,m
      rho = var(4,i)
      radius = var(3,i)*R_sun
      delta_t = var(32,i)
      nu_diff = var(7,i)
      if (i == 1) then
        !! time derivative term
        C1(i) = rho*(radius**4)

        !! diffusion terms
        lambda(i) = (nu_diff*delta_t/2d0)*&
          &(rho*radius**4 + var(4,i+1)*(var(3,i+1)*R_sun)**4)&
          &/((var(3,i+1)*R_sun-radius)*(var(3,i+1)*R_sun-radius))

        beta(i) = 0d0

        !! mixed modes term jdot
        alpha(i) = radius*radius*var(5,i)*delta_t

        !! contraction and expasion of the star
        epsilon(i) = R_sun*(var(3,i)-var(8,i))*rho*(radius*radius*radius)

        !! j dot mid point
        !a1(i) = 4d0*(var(5,i)*delta_t)
        !a2(i) = 4d0*(var(5,i)*delta_t)

        a1(i) = 0d0!0.5d0*(var(5,i-1)*delta_t)/((var(4,i-1))*&
        !&(var(3,i-1)*R_sun)**2)
        a2(i) = 0.5d0*var(5,i)*delta_t*radius*radius

      else if (i == m) then

        C1(i) = rho*(radius**4)

        lambda(i) = 0d0

        beta(i) = (nu_diff*delta_t/2d0)*&
          &(rho*radius**4 + var(4,i-1)*(var(3,i-1)*R_sun)**4)&
          &/((radius-var(3,i-1)*R_sun)*(radius-var(3,i-1)*R_sun))

        alpha(i) = radius*radius*(var(5,i)*delta_t)

        epsilon(i) = (radius-var(8,i)*R_sun)*rho*(radius*radius*radius)

        !! j dot mid point
        !a1(i) = 4d0*(var(5,i-1)*delta_t)
        !a2(i) = 4d0*(var(5,i)*delta_t)

        a1(i) = radius*radius*0.5d0*(var(5,i-1)*delta_t)
        a2(i) = radius*radius*0.5d0*(var(5,i)*delta_t)

      else

        C1(i) = rho*(radius**4)

        lambda(i) = (nu_diff*delta_t/2d0)*&
          &(rho*radius**4 + var(4,i+1)*(var(3,i+1)*R_sun)**4)&
          &/((var(3,i+1)*R_sun-radius)*(var(3,i+1)*R_sun-var(3,i-1)*R_sun))

        beta(i) = (nu_diff*delta_t/2d0)*&
          &(rho*radius**4 + var(4,i-1)*(var(3,i-1)*R_sun)**4)&
          &/((radius-var(3,i-1)*R_sun)*(var(3,i+1)*R_sun-var(3,i-1)*R_sun))

        alpha(i) = radius*radius*(var(5,i)*delta_t)

        epsilon(i) = (radius-var(8,i)*R_sun)*rho*(radius*radius*radius)

        !! j dot mid point
        !a1(i) = 4d0*(var(5,i-1)*delta_t)
        !a2(i) = 4d0*(var(5,i)*delta_t)

        a1(i) = radius*radius*0.5d0*(var(5,i-1)*delta_t)
        a2(i) = radius*radius*0.5d0*(var(5,i)*delta_t)

      end if
    end do
  end subroutine get_CNcoeff

  subroutine create_matrix_diffusionCN_new(m,matrix_A,matrix_B,var)
    implicit none
    real (DP), intent(in) :: var(:,:)
    real (DP), intent(inout) :: matrix_A(:,:),matrix_B(:,:)
    integer, intent(in) :: m 
    real (DP), allocatable :: lambda(:),beta(:),alpha(:),epsilon(:),C1(:)
    real (DP), allocatable :: a1(:), a2(:)
    integer :: i
    

    allocate(lambda(m),beta(m),alpha(m),epsilon(m),a1(m),a2(m),C1(m))

    call get_CNcoeff(var,lambda,beta,alpha,epsilon,C1,a1,a2)

    !! diffusion term
    !lambda = 0.
    !beta = 0.
    !! contraction and expansion term
    !epsilon = 0.
    !! jdot term
    alpha = 0d0
    !! j dot mid point term
    a1 = 0d0
    a2 = 0d0

    !!fill matrix_A
    matrix_A = 0d0
    do i=2,m-1
      matrix_A(1,i)= -beta(i)-a1(i)
      matrix_A(2,i)= C1(i)+lambda(i)+beta(i)+epsilon(i)-alpha(i)-a2(i)
      matrix_A(3,i)= -lambda(i)
    end do
    !!Boundary conditions
    matrix_A(1,1) = 0d0
    matrix_A(2,1) = C1(1)+lambda(1)+epsilon(1)-alpha(1)-a2(1)-a1(1)
    matrix_A(3,1) = -lambda(1)

    matrix_A(1,m) = -beta(m)-a1(m)
    matrix_A(2,m) = C1(m)+beta(m)+epsilon(m)-alpha(m)-a2(m)
    matrix_A(3,m) = 0d0

    !!fill matrix_B
    matrix_B = 0d0
    do i=2,m-1
      matrix_B(1,i) = beta(i)+a1(i)
      matrix_B(2,i) = C1(i)-lambda(i)-beta(i)-epsilon(i)+alpha(i)+a2(i)
      matrix_B(3,i) = lambda(i)
    end do
    !!Boundary conditions
    matrix_B(1,1) = 0d0
    matrix_B(2,1) = C1(1)-lambda(1)-epsilon(1)+alpha(1)+a2(1)+a1(1)
    matrix_B(3,1) = lambda(1)

    matrix_B(1,m) = beta(m)+a1(m)
    matrix_B(2,m) = C1(m)-beta(m)-epsilon(m)+alpha(m)+a2(m)
    matrix_B(3,m) = 0d0

    !call printMatrix(matrix_A, 3, m,1000)

    deallocate(lambda,beta,alpha,epsilon,a1,a2,C1)

  end subroutine create_matrix_diffusionCN_new

end module numerical_scheme