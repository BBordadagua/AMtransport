Module numerical_scheme
  use parameters
  use data_handling
  use math_functions
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: variables
  public :: CrankNicolson_new,create_matrix_diffusionCN_new,get_CNcoeff
  public :: difeq,difeq_total,constants,solvde
  contains

  subroutine variables(Mn,Mo,delta_t,var,y)
    implicit none
    real (DP), intent(inout) :: var(:,:)
    real (DP), intent(in) :: Mn(:,:), delta_t, Mo(:,:),y(:,:)
    integer :: k

    var(1,:) = Mn(2,:)**(2d0/3d0) !! nu 
    var(2,:) = Mn(2,:) !! mass (M_sun)
    var(3,:) = 10**Mn(3,:) !! radius (R_sun)
    var(4,:) = 10**Mn(5,:) !! rho
    var(5,:) = Mo(6,:) !! jmodes
    var(6,:) = 10**Mn(10,:) !! g surf
    var(7,:) = 1d7 !! Dv incomplete
    var(8,:) = 10**Mo(2,:) !! radius old
    var(9,:) = Mo(4,:) !! omega old
    var(10,:) = Mn(15,:)    !! cp
    var(11,:) = 10**Mn(4,:)    !! temperature
    var(12,:) = 10**Mn(6,:)    !! pressure
    var(13,:) = Mn(27,:)    !! luminosity (L_sun)
    var(14,:) = Mn(14,:)    !! nabla_ad
    var(15,:) = Mn(28,:)    !! nabla !!not sure!!

    var(34,:) = Mn(25,:) !!kappa oppacity
    var(35,:) = Mn(13,:) !! mu

    var(18,:) = -Mn(17,:)/Mn(16,:)   !! delta 
    var(19,:) = (2d0*var(3,:)*var(3,:)*var(9,:))/(3d0*var(6,:))    !! theta_old 
    !var(20,:) = y(1,1)    !! omega_0  !!should be omega actual !! change
    var(21,:) = 10**Mn(5,1)    !! rho_0
    var(22,:) = 1.5d-6 !! beta for Dh (Richard & Zahn 1999) !!accouting for the dividion by 10
    var(23,:) = Mn(11,:)    !! HT 
    var(24,:) = (4d0*7.5646d-15*c_speed*var(11,:)**3)/(3d0*var(4,:)*var(4,:)*var(10,:)*var(34,:))    !! grand K thermal conductivity
    var(25,:) = Mn(20,:)    !! eps_nuc 
    var(26,:) = Mn(21,:)    !! eps_t  !!not sure!!
  
    var(28,:) = 3d0-Mn(26,:)+var(18,:)  !! chi_t !!not sure!!
    
    var(30,:) = Mn(12,:)    !! HP
    !var(31,:) = 0d0    !! lambda old
    var(32,:) = delta_t!! delta_t

    var(38,:) = Mn(32,:) !! N2_t
    var(39,:) = Mn(33,:) !! N2_mu

    !! variables computed using derivatives !!only for interior points
    var(16,1) = 0d0 !check
    var(17,1) = 0d0 !check
    var(27,1) = 0d0 !check
    var(29,1) = 0d0 !check
    var(33,1) = 0d0 !check
    var(34,1) = 0d0 !check
    var(35,1) = 0d0 !check
    do k=2,size(var,dim=2)
      var(16,k) = (dlog(var(35,k))-dlog(var(35,k-1)))/(dlog(var(12,k))-dlog(var(12,k-1)))     !! nabla_mu
      var(17,k) = (dlog(var(4,k))-dlog(var(4,k-1)))/(dlog(var(35,k))-dlog(var(35,k-1)))    !! phi   !constant pressure and T
      var(27,k) = (dlog(var(20,k))-dlog(var(20,k-1)))/(dlog(var(35,k))-dlog(var(35,k-1)))   !! eps_mu !!not sure because P and T constant
      var(29,k) =  -(dlog(var(34,k))-dlog(var(34,k-1)))/(dlog(var(35,k))-dlog(var(35,k-1))) -var(17,k)  !! chi_mu
      var(33,k) =  (dlog(var(13,k))-dlog(var(13,k-1)))&
      &/(dlog(var(2,k))-dlog(var(2,k-1)))!! dln L/dln m !!not sure if should be in Msun and Lsun
      !var(36,k) = (y(3,k)-y(3,k-1))/((var(1,k)-var(1,k-1))) !! for Dh
      var(37,k) = (dlog(var(4,k))-dlog(var(4,k-1)))/((var(1,k)-var(1,k-1))) !! for Dh
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
    real (DP),allocatable :: y1(:),y2(:),y3(:),y4(:),y5(:)
    real (DP) :: nu_star

    allocate(nu12(size(var_old,dim=2)),nu32(size(var_old,dim=2)),&
    &delta_nu(size(var_old,dim=2)),var_k(size(var_old,dim=1),size(var_old,dim=2)),&
    &var(size(var_old,dim=1),size(var_old,dim=2)),&
    &y1(size(y,dim=2)),y2(size(y,dim=2)),y3(size(y,dim=2)),y4(size(y,dim=2)),y5(size(y,dim=2)))
    
    var_k = var_old
    do k=k1+1,k2
      delta_nu(k) = var_k(1,k) - var_k(1,k-1)
      var(:,k) = 0.5d0*(var_old(:,k)+var_old(:,k-1))
    end do
    delta_nu(k1) = 0d0
    nu12 = dsqrt(var(1,:))
    nu32 = var(1,:)**(3d0/2d0)

    !! problem dividing by delta_nu
    !! problem when k-1
        
    cte(1,:) = delta_nu*nu12
    cte(2,:) = (cte(1,:)*var(3,:)*var(3,:))/var_k(32,:)  !C1 cesam
    cte(3,:) = (cte(1,:)*var(8,:)*var(8,:)*var(9,:)) /var_k(32,:)  !C2 cesam
    cte(4,:) = (64d0*PI*PI*(R_sun**4))/(9d0*M_sun*M_sun) !! diffusion !C18 cesam
    
    cte(7,:) = cte(1,:)*var(5,:)/(R_sun*R_sun*var(4,:)) !! mixed modes !C33 cesam
    cte(8,:) = (8d0*PI*R_sun*R_sun)/(15d0*M_sun)
    cte(9,:) = cte(8,:)*var_k(4,:)*var_k(3,:)**4 !C3 cesam

    cte(11,:) = (M_sun/L_sun)*(nu32*var(6,:)*var(10,:)*var(4,:)*var(11,:)/(var(12,:)*var(13,:)))*&
              &(var(14,:)-var(15,:)+var(17,:)*var(16,:)/var(18,:)) !C7 cesam
    cte(12,:) = -(M_sun/(L_sun*var(32,:)))*(nu32*var(10,:)*var(11,:)*var(9,:)*var(19,:)/(var(13,:)*var(18,:))) !C8 cesam
    cte(13,:) = 8d0*R_sun*var(3,:)/(3d0*var(6,:)) !C9 cesam
    cte(14,:) = var(33,:) - var(20,:)*var(20,:)/(2d0*PI*G*var(21,:)) !C10 cesam
    cte(15,:) = 1d0/(2d0*PI*G*var(4,:)) !C11 cesam
    cte(16,:) = (2d0*var(1,:))/(3d0*delta_nu) !C12 cesam

    cte(18,:) = (M_sun/(2d0*PI*R_sun**3))*(nu32/(var(4,:)*var(3,:)**3)) !C2_star cesam
    cte(19,:) = (M_sun/L_sun)*(var(25,:)*nu32/var(13,:)) !C3_star cesam
    cte(20,:) = var(33,:) !C4_star cesam !C15 cesam
    cte(21,:) = (3d0*M_sun/(2d0*PI*R_sun**4))*(nu32*var(23,:)/(var(4,:)*var(3,:)**4)) !C5_star cesam
    cte(22,:) = (M_sun/L_sun)*(var(25,:)*var(26,:)*nu32/var(13,:)) !C6_star cesam
    cte(23,:) = cte(17,:)-cte(18,:)+cte(19,:)-cte(20,:) !C13 cesam
    cte(24,:) = (M_sun/(L_sun*var(32,:)))*(nu32*var(10,:)*var(11,:)/(var(13,:)*var(18,:))) !C14 cesam
    cte(25,:) = (M_sun/L_sun)*(var(25,:)*var(27,:)*nu32/var(13,:)) !C16 cesam
    cte(26,:) = cte(21,:) - cte(22,:) ! C17 cesam
    cte(27,:) = -(8d0*PI*R_sun*R_sun/(3d0*M_sun))*(var(3,:)*var(3,:)*var(4,:)*var(23,:)) !C19 cesam
    cte(28,:) = cte(1,:)*(1d0+var(28,:)) !C20 cesam
    cte(29,:) = cte(1,:)*var(29,:) !C21 cesam
    cte(30,:) = var_k(17,:)/var_k(18,:) !C22 cesam

    cte(32,:) = var(17,:)/var(18,:) !C24 cesam
    cte(33,:) = 1d0/var_k(18,:) !C25 cesam
    
    cte(35,:) = 1d0/var(18,:) !C27 cesam
    
    cte(37,:) = 1d0/var(32,:) !C29 cesam
    
    cte(40,:) = (2d0*R_sun**3/(3d0*G*M_sun))*(var(3,:)**3/nu32) !C32 cesam
    cte(41,:) = ((16d0*PI*R_sun**4)/9d0*M_sun)*((var(4,:)*var(3,:)**4)/var(6,:)) 
    cte(42,:) = ((16d0*PI*R_sun**4)/9d0*M_sun)*((var_k(4,:)*var_k(3,:)**4)/var_k(6,:))

    cte(45,:) = (var(22,:)*8d0*PI*R_sun**6 /(9d0*M_sun))*(var(4,:)*var(3,:)**6 *y1)*&
    & (var_k(34,:)/nu12 + y3*var_k(35,:)/nu12-3d0*y2*y3/(2d0*y1))&
    & - (var(22,:)*R_sun**4 *var(3,:)**4 *y1/3d0)*(y3/(nu12*R_sun*var(3,:))) !!Dh middle point

    cte(17,:) = (3d0*M_sun/(2d0*PI*R_sun**4))*(nu32*cte(45,:)*var(23,:)/(var(4,:)*var(24,:)*var(18,:)*var(3,:)**4))!C1_star cesam


    do k=k1+1,k2
      y1 = 0.5d0*(y(1,k)+y(1,k-1))
      y2 = 0.5d0*(y(2,k)+y(2,k-1))
      !y3 = 0.5d0*(y(3,k)+y(3,k-1))
      
      cte(10,k) = cte(8,k)*var_k(4,k-1)*var_k(3,k-1)**4 !C4 cesam
      cte(31,k) = var_k(17,k-1)/var_k(18,k-1) !C23 cesam
      cte(34,k) = 1d0/var_k(18,k-1) !C26 cesam
      cte(36,k) = var_k(16,k-1)/var_k(30,k-1) !C28 cesam
      
      cte(39,k) = var(31,k-1)/var(32,k-1) !C31 cesam
      cte(43,k) = ((16d0*PI*R_sun**4)/9d0*M_sun)*((var_k(4,k-1)*var_k(3,k-1)**4)/var_k(6,k-1))
      !cte(44,k-1) = (var(22,k-1)*8d0*PI*R_sun**6 /(9d0*M_sun))*(var_k(4,k-1)*var_k(3,k-1)**6 *y(1,k-1))*&
      !& (var(36,k-1)/dsqrt(var_k(1,k-1)) + y(3,k-1)*var(37,k-1)/dsqrt(var_k(1,k-1))-3d0*y(2,k-1)*y(3,k-1)/(2d0*y(1,k-1)))&
      !& - (var(22,k-1)*R_sun**4 *var_k(3,k-1)**4 *y(1,k-1)/3d0)*(y(3,k-1)/(dsqrt(var_k(1,k-1))*R_sun*var_k(3,k-1))) !!Dh(k-1)
      cte(38,k) = (6d0/(R_sun*R_sun))*(dsqrt(cte(44,k-1))/(var_k(3,k-1)*var_k(3,k-1))) !C30 cesam

    end do
    !cte(44,k2) = (var(22,k2)*8d0*PI*R_sun**6 /(9d0*M_sun))*(var_k(4,k2)*var_k(3,k2)**6 *y(1,k2))*&
    !  & (var(36,k2)/dsqrt(var_k(1,k2)) + y(3,k2)*var(37,k2)/dsqrt(var_k(1,k2))-3d0*y(2,k2)*y(3,k2)/(2d0*y(1,k2)))&
    !  & - (var(22,k2)*R_sun**4 *var_k(3,k2)**4 *y(1,k2)/3d0)*(y(3,k2)/(dsqrt(var_k(1,k2))*R_sun*var_k(3,k2))) !!Dh(k-1)

    !var(7,:) = (1d0/6d0)*(64d0*PI*PI*R_sun**6 *var(3,:)**6 *var(4,:)**2 *y(2,:)*y(2,:))/(9d0*M_sun*M_sun) !! Dv incomplete
    !cte(46,:) = var_k(7,:)*dsqrt(cte(44,:))*(1d0+var_k(24,:)/dsqrt(cte(44,:)))/&
    !&(var_k(38,:)+var_k(39,:)*(1d0+var_k(24,:)/dsqrt(cte(44,:)))) !!Dv
    cte(46,:) = 1d7 !!Dv
    cte(5,:) = cte(4,:)*(cte(46,:)*(var_k(3,:)**6)*(var_k(4,:)**2))  !! diffusion !C5 cesam

    do k=k1+1,k2
      cte(6,k) = cte(4,k)*(cte(46,k-1)*(var_k(3,k-1)**6)*(var_k(4,k-1)**2)) !! diffusion !C6 cesam
    end do

    !! constants for boundary conditions !! local conservation of AM in CZs
    nu_star = (M_star/M_sun)**(2d0/3d0) !! nu_star cesam
    cte(47,k2) = 2d0*var(3,k2)*var(3,k2)*((M_star/M_sun)-nu32(k2))/3d0 !! Itop cesam
    cte(48,:) = 2d0*var(3,k2)*var(3,k2)*((M_star/M_sun)-nu32(k2))/3d0 !! Itop_t cesam
    cte(49,:) = 1d0/var_k(32,:) !! K1 cesam
    cte(50,:) = (8d0*PI*R_sun*R_sun/(15d0*M_sun))*var(4,:)*var(3,:)**4 !! K2 cesam
    cte(51,:) = 0d0 !! K3 cesam
    cte(52,:) = 0d0 !! K4 cesam
    cte(53,:) = 0d0 !! K5 cesam


    deallocate(nu12,nu32,delta_nu,var_k,var,y1,y2,y3,y4,y5)
  end subroutine constants


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
    
      !call printMatrix(s,2,5,7000)
    else if(k.gt.k2) then !Boundary conditions at last point. 
      s(1,2+indexv(1))=0d0 
      s(1,2+indexv(2))=1d0
      s(1,jsf)=y(2,k2)

      s(2,2+indexv(1))=0d0 
      s(2,2+indexv(2))=1d0
      s(2,jsf)=y(2,k2)
  
      !call printMatrix(s,2,5,2000)
    else !Interior point. 
      
      s(1,indexv(1))=0.5d0*(cte(2,k)-cte(7,k))   !E1
      s(1,indexv(2))=cte(6,k)
      s(1,2+indexv(1))= 0.5d0*(cte(2,k)-cte(7,k))
      s(1,2+indexv(2))=-cte(5,k)

      s(2,indexv(1))=1d0     !E2
      s(2,indexv(2))=0.5d0*cte(1,k)
      s(2,2+indexv(1))=-1d0
      s(2,2+indexv(2))=0.5d0*cte(1,k)

      !! expression with rdot linear interpolation
      s(1,jsf)=0.5d0*cte(2,k)*(y(1,k)+y(1,k-1))-cte(3,k)-cte(5,k)*y(2,k)+cte(6,k)*y(2,k-1)&
      &-0.5d0*cte(7,k)*(y(1,k)+y(1,k-1))
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
    
    !s=0d0 !! added makes matrix singular
    
    if(k.eq.k1) then !Boundary condition at first point. 
      !! considering only radiative center
      s(1,2+indexv(1))=0d0 
      s(1,2+indexv(2))=1d0  !! should be a factor constant??
      s(1,2+indexv(3))=0d0 
      s(1,2+indexv(4))=0d0 
      s(1,2+indexv(5))=0d0 
      s(1,jsf)=y(2,k1) 

      s(2,2+indexv(1))=0d0 
      s(2,2+indexv(2))=0d0 
      s(2,2+indexv(3))=1d0 
      s(2,2+indexv(4))=0d0 
      s(2,2+indexv(5))=0d0 
      s(2,jsf)=y(3,k1)
    
      !call printMatrix(s,2,5,7000)
    else if(k.gt.k2) then !Boundary conditions at last point. 
      s(1,2+indexv(1))=cte(49,k2)*cte(47,k2) + cte(50,k2)*y(3,k2)-cte(52,k2)
      s(1,2+indexv(2))=0d0 
      s(1,2+indexv(3))=cte(50,k2)*y(1,k2) 
      s(1,2+indexv(4))=0d0 
      s(1,2+indexv(5))=0d0 
      s(1,jsf)=cte(49,k2)*(cte(47,k2)*y(1,k2)-cte(48,k2)) + cte(50,k2)*y(1,k2)*y(3,k2)-cte(51,k2)

      s(2,2+indexv(1))=0d0 
      s(2,2+indexv(2))=1d0 !! should be a factor constant??
      s(2,2+indexv(3))=0d0 
      s(2,2+indexv(4))=0d0 
      s(2,2+indexv(5))=0d0 
      s(2,jsf)=y(2,k2) 

      s(3,2+indexv(1))=0d0 
      s(3,2+indexv(2))=0d0 
      s(3,2+indexv(3))=0d0 
      s(3,2+indexv(4))=0d0 
      s(3,2+indexv(5))=1d0 
      s(3,jsf)=y(5,k2)
  
      !call printMatrix(s,2,5,2000)
    else !! Interior points

      y1 = 0.5d0*(y(1,k)+y(1,k-1))
      y2 = 0.5d0*(y(2,k)+y(2,k-1))
      y3 = 0.5d0*(y(3,k)+y(3,k-1))
      y4 = 0.5d0*(y(4,k)+y(4,k-1))
      y5 = 0.5d0*(y(5,k)+y(5,k-1))
      psi = cte(35,k)*y1*cte(41,k)*y2 - cte(32,k)*y5
      psik = cte(33,k)*y(1,k)*cte(42,k)*y(2,k) - cte(30,k)*y(5,k)
      psik_1 = cte(34,k)*y(1,k-1)*cte(43,k)*y(2,k-1)-cte(31,k)*y(5,k-1)
      
      !! Omega y(1,:)
      !! derivatives dE/dy(k-1)
      s(1,indexv(1))=0.5d0*(cte(2,k)-cte(7,k))+cte(10,k)*y(3,k-1)
      s(1,indexv(2))=cte(6,k)
      s(1,indexv(3))=cte(10,k)*y(1,k-1)
      s(1,indexv(4))=0d0
      s(1,indexv(5))=0d0
      !! derivatives dE/dy(k)
      s(1,2+indexv(1))=0.5d0*(cte(2,k)-cte(7,k))-cte(9,k)*y(3,k)
      s(1,2+indexv(2))=-cte(5,k)
      s(1,2+indexv(3))=-cte(9,k)*y(1,k)
      s(1,2+indexv(4))=0d0
      s(1,2+indexv(5))=0d0
      !! equation E1
      s(1,jsf) = cte(2,k)*y1 - cte(3,k)&
        &-cte(9,k)*y(1,k)*y(3,k) + cte(10,k)*y(1,k-1)*y(3,k-1)& !! meridional circ
        &-cte(5,k)*y(2,k) + cte(6,k)*y(2,k-1)&                  !! diffusion
        &-cte(7,k)*y1                                           !! mixed modes
      
      !! Theta = domega/dr y(2,:)
      !! derivatives dE/dy(k-1)
      s(2,indexv(1))=1d0
      s(2,indexv(2))=0.5d0*cte(1,k)
      s(2,indexv(3))=0d0
      s(2,indexv(4))=0d0
      s(2,indexv(5))=0d0
      !! derivatives dE/dy(k)
      s(2,2+indexv(1))=-1d0
      s(2,2+indexv(2))=0.5d0*cte(1,k)
      s(2,2+indexv(3))=0d0
      s(2,2+indexv(4))=0d0
      s(2,2+indexv(5))=0d0
      !! equation E2
      s(2,jsf) = cte(1,k)*y2 - (y(1,k)-y(1,k-1))

      !! U2 y(3,:)
      !! derivatives dE/dy(k-1)
      s(3,indexv(1))=cte(13,k)*y1*(2d0*cte(15,k)*y1*y1-1d0+cte(14,k)) &
                    & -0.5d0*(cte(23,k)+cte(24,k)*(1d0-cte(40,k)*y1*y1)+cte(26,k)*cte(35,k))*cte(41,k)*y2 &
                    & - cte(40,k)*y1*(cte(11,k)*y3-cte(12,k)-cte(24,k)*cte(41,k)*y1*y2)  
      s(3,indexv(2))=-0.5d0*(cte(23,k)+cte(24,k)*(1d0-cte(40,k)*y1*y1)+cte(26,k)*cte(35,k))*cte(41,k)*y1
      s(3,indexv(3))=0.5d0*cte(11,k)*(1d0-cte(40,k)*y1*y1)
      s(3,indexv(4))=cte(16,k) - 0.5d0*cte(20,k)
      s(3,indexv(5))=0.5d0*(cte(26,k)*cte(32,k)-cte(25,k))
      !! derivatives dE/dy(k)
      s(3,2+indexv(1))=s(3,indexv(1))
      s(3,2+indexv(2))=s(3,indexv(2))
      s(3,2+indexv(3))=s(3,indexv(3))
      s(3,2+indexv(4))=s(3,indexv(4)) - 2d0*cte(16,k)
      s(3,2+indexv(5))=s(3,indexv(5))
      !! equation E3
      s(3,jsf) = cte(11,k)*(1d0-cte(40,k)*y1*y1)*y3 &
        & - cte(12,k)*(1d0-cte(40,k)*y1*y1) &
        & - cte(13,k)*y1*y1*(1d0-cte(14,k)-cte(15,k)*y1*y1) &
        & - cte(16,k)*(y(4,k)-y(4,k-1)) &
        & - (cte(23,k)+cte(25,k)*(1d0-cte(40,k)*y1*y1))*y1*cte(41,k)*y2 &
        & - cte(20,k)*y4 - cte(25,k)*y5 - cte(26,k)*psi

      !! A y(4,:)
      !! derivatives dE/dy(k-1)
      s(4,indexv(1))=cte(27,k)*cte(34,k)*cte(43,k)*y(2,k-1)+0.5d0*(1d0-cte(28,k)*cte(35,k))*cte(41,k)*y2
      s(4,indexv(2))=cte(27,k)*cte(34,k)*cte(43,k)*y(1,k-1)+0.5d0*(1d0-cte(28,k)*cte(35,k))*cte(41,k)*y1
      s(4,indexv(3))=0d0
      s(4,indexv(4))=0.5d0
      s(4,indexv(5))=0.5d0*(cte(29,k)+cte(28,k)*cte(32,k)) - cte(27,k)*cte(31,k)
      !! derivatives dE/dy(k)
      s(4,2+indexv(1))=-cte(27,k)*cte(33,k)*cte(42,k)*y(2,k)+0.5d0*(1d0-cte(28,k)*cte(35,k))*cte(41,k)*y2
      s(4,2+indexv(2))=-cte(27,k)*cte(33,k)*cte(42,k)*y(1,k)+0.5d0*(1d0-cte(28,k)*cte(35,k))*cte(41,k)*y1
      s(4,2+indexv(3))=0d0
      s(4,2+indexv(4))=0.5d0
      s(4,2+indexv(5))=0.5d0*(cte(29,k)+cte(28,k)*cte(32,k)) + cte(27,k)*cte(30,k)
      !! equation E4
      s(4,jsf) = cte(1,k)*y4 - cte(27,k)*(psik-psik_1) &
        & - cte(28,k)*psi + cte(29,k)*y5 &
        & + cte(1,k)*y1*cte(41,k)*y2

      !! Lambda y(5,:)
      !! derivatives dE/dy(k-1)
      s(5,indexv(1))=0d0 
      s(5,indexv(2))=0d0
      s(5,indexv(3))=-cte(36,k)
      s(5,indexv(4))=0d0
      s(5,indexv(5))=(cte(37,k)+cte(38,k))
      !! derivatives dE/dy(k)
      s(5,2+indexv(1))=0d0
      s(5,2+indexv(2))=0d0
      s(5,2+indexv(3))=0d0
      s(5,2+indexv(4))=0d0
      s(5,2+indexv(5))=0d0
      !! equation E5
      s(5,jsf) = (cte(37,k)+cte(38,k))*y(5,k-1) &
        & - cte(39,k) - cte(36,k)*y(3,k-1)

      !call printMatrix(s,2,5,3000)
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
    & y,c,s,cte,var)
    implicit none 
    real (DP), intent(in) :: var(:,:)
    real (DP), intent(in) :: conv,slowc !! changed
    real (DP), intent(inout) :: c(:,:,:),s(:,:),y(:,:),scalv(:),cte(:,:)
    integer, intent(in) :: itmax,m,ne,nb,indexv(:)
    integer :: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8, &
    & j9,jc1,jcf,jv,k,k1,k2,km,kp,nvars,kmax(ne)
    real (DP) :: err,errj,fac,vmax,vz,ermax(ne)

    s=0d0
    c=0d0


    k1=1 !Set up row and column markers. 
    k2=m 
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
    
    do it=1,itmax 
      call constants(k1,k2,cte,var,y)
      !Primary iteration loop. 
      k=k1 !Boundary conditions at first point. 
      call difeq(k,k1,k2,j9,indexv,s,y,cte) 
      call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s) 
      do k=k1+1,k2 
        !Finite difference equations at all point pairs. 
        kp=k-1 
        call difeq(k,k1,k2,j9,indexv,s,y,cte) 
        call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s) 
        call pinvs(ic1,ic4,j3,j9,jc1,k,c,s) 
      enddo 
      k=k2+1 !Final boundary conditions. 
      call difeq(k,k1,k2,j9,indexv,s,y,cte) 
      call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s) 
      call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
      call bksub(ne,nb,jcf,k1,k2,c) !Backsubstitution. 
      err=0d0


      do j=1,ne !Convergence check, accumulate average error. 

      !! scalv(1:nyj) contains typical sizes for each dependent variable, used to weight errors.
        !! cesam2k20 version
        !scalv = maxval(dabs(y(1:nrot_eq,:)),2) !! Cesam2k20
        scalv(j) = maxval(dabs(y(j,:))) !! each variable has different scaling
        !if (scalv(j) < 1.0d-20) scalv = 1.0d0

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
      !write(*,100) 'it',it,'   err',err,'   fac',fac
      !! Summary of corrections for this step. Point with largest error for each 
      !! variable can be monitored by writing out kmax and ermax. 
      !print*, 'jv ', 1, 'kmax ',kmax(1), 'ermax ', ermax(1)
      !print*, 'jv ', 2, 'kmax ',kmax(2), 'ermax ', ermax(2)
      if(err.lt.conv) return 
    enddo
    !100 format(1x,i4,2f12.6)
    !100 format(A3,i4,A6,ES14.4,A6,ES14.3)


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
    real (DP), allocatable :: U(:),V(:),mA_triag(:,:),mB_triag(:,:),dr(:),old_radius(:)
    real (DP) :: omega(m)
    integer :: j
    real (DP), allocatable :: var(:,:)

    !!allocate matrix and vector
    allocate(mA_triag(3,m),mB_triag(3,m))
    allocate(U(m),V(m),dr(m),old_radius(m))

    !! compute constants
    allocate(var(50,m)) !! number of constants
    call variables(Mn,Mo,delta_t,var,var)

    dr = (var(3,:) - var(8,:))/iter !! divide the radius into smaller steps according to smaller timesteps
    V = var(9,:) !omega old model
    old_radius = var(8,:)
    var(3,:) = old_radius + dr !! radius current model
    
    do j=1,iter    !!subtimesteps
      if (j/=1) then
        var(8,:) = var(3,:)
        var(3,:) = old_radius + j*dr
      end if
      
      !! Create matrix
      call create_matrix_diffusionCN_new(m,mA_triag,mB_triag,var)

      !! multiplying matrix
      U = 0d0
      !U = matmul(transpose(matrix_B),V)
      U = multiply_matrix_vector_triag(mB_triag,V,m)

      !!invert Tridiagonal matrix
      omega = 0d0
      call invert_tridagmatrix(mA_triag,m,omega,U)

      V = omega
    end do

    deallocate(U,V,mA_triag,mB_triag,dr,old_radius,var)

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