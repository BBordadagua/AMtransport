!!
!!                 coefficients_mod.f90
!!
!!  This program contains the module 'coefficient_equations' with functions
!!      and subroutines that calculate the coeffcients A,B,C,D
!!    of the flux equation.
!!
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!

MODULE coefficient_equations
  use parameters
  use data_handling
  use math_functions
  implicit none
  
  !define public subroutines and functions
  public :: calculate_coefficients
  public :: derive,compute_J,compute_K0,compute_K1,compute_K2,derive_2ndorder,calculate_gwaves,calculate_gwaves_simplified
  public :: compute_zeta0,compute_zeta1,compute_zeta2,compute_zeta3,computeA,computeG,compute_cs2
  public :: calculate_simplifiedeq,compute_cp,alpha,computeKr2,getL,getOmega
  public :: derive_moreprecision
  contains

  subroutine calculate_gwaves(array,m,summaryfile,k,coeff,MESAprofile)
    implicit none
    real (DP), intent(in) :: array(:,:), summaryfile(:,:),MESAprofile(:,:)
    real (DP), intent(inout) :: coeff(:)
    integer, intent(in) :: m,k
    real (DP) ,dimension(m) :: dr,r, omega,Lr, firstterm,secondterm,thirdterm,alph
    real (DP) ,dimension(m) :: zeta0,zeta1,zeta2,zeta3,kr2,A,G_func, um,dois,tres,quatro
    real (DP) ,dimension(m) :: x,xi_r,Gamma_1,omega_rot,P,rho,V,V_2,As
    real (DP) :: omega_R, K0,K1,K2, l, m_mode,freq, dev1tot
    integer :: i

    l      = summaryfile(8,k)
    m_mode = summaryfile(9,k)
    freq   = summaryfile(5,k)
    omega_R = freq*10**(-6.)*2.*PI !rad/s
    Lr = getL(MESAprofile,array,m)
    omega_rot = getOmega(MESAprofile,array,m)

    if ( k == 369 ) then
      open(20,action='write',file='output/detail.l'//trim(string(int(l)))//'.j'//trim(string(k))//'.txt')
    end if
    
    do i=1,m
    !put the variables into arrays
      x(i)         = array(15,i)
      xi_r(i)      = array(18,i)
      V_2(i)       = array(7,i)
      As(i)        = array(1,i)
      Gamma_1(i)   = array(2,i)
      P(i)         = array(5,i)
      rho(i)       = array(14,i)

    !put the functions into arrays
      r(i) = array(15,i)*R_star
      kr2(i) = computeKr2(array,i,l,freq)
      K0 = compute_K0(m_mode,l)
      K1 = compute_K1(m_mode,l)
      K2 = compute_K1(m_mode,l)
      A(i) = computeA(array,i)
      zeta0(i) = compute_zeta0(array,i,Lr(i))
      zeta1(i) = compute_zeta1(array,i,l)
      zeta2(i) = compute_zeta2(array,i,l,freq,Lr(i))
      zeta3(i) = compute_zeta3(array,i,l)
      G_func(i) = computeG(array,i,freq,l)
      alph(i) = alpha(array,i,Lr(i))

      V(i) = x(i)**2 * V_2(i)

    end do

    do i=1,m-1
      dr(i) = r(i+1) - r(i)
      dev1tot = 2.*r(i)*omega_rot(i) + (r(i)**2)*derive(omega_rot(i),omega_rot(i+1),dr(i))

      !positive m
      omega(i) = omega_R + m_mode*omega_rot(i)
      firstterm(i) = (r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i)/omega(i))*(&
      & - m_mode*l*(l+1.)*G_func(i)*Gamma_1(i)*P(i)*omega_R/omega(i) &
      & - m_mode*Gamma_1(i)*P(i) &
      & - 2.*rho(i)*omega_R*G_func(i)*r(i)**2 * K2 *omega_rot(i))!*(1./(8.*PI))

      secondterm(i) = -(m_mode*r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i))*(&
      & l*(l+1.)*Gamma_1(i)*P(i)*G_func(i)/omega(i) &
      & + Gamma_1(i)*P(i)/omega_R &
      & - omega(i)*rho(i)*G_func(i)*r(i)**2)!*(1./(8.*PI))

      thirdterm(i) = (zeta2(i)*((xi_r(i))**2)*r(i)**2)*(&
      !& -G_func(i)*l*(l+1)*m_mode*omega_R*Gamma_1(i)*P(i)*zeta3(i)/omega(i) &
      !& + m_mode*omega_R*Gamma_1(i)*P(i)*(-zeta1(i) - 2./r(i) - A(i))/omega(i) &
      & - rho(i)*omega_R*m_mode*omega(i)*r(i)**2 *G_func(i)*zeta3(i) &
      & - rho(i)*omega_R*K1*dev1tot)!*(1./(8.*PI))

      !negative m
      omega(i) = omega_R - m_mode*omega_rot(i)
      firstterm(i) = firstterm(i) + (r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i)/omega(i))*(&
      & + m_mode*l*(l+1.)*G_func(i)*Gamma_1(i)*P(i)*omega_R/omega(i) &
      & + m_mode*Gamma_1(i)*P(i) &
      & - 2.*rho(i)*omega_R*G_func(i)*r(i)**2 * K2 *omega_rot(i))!*(1./(8.*PI))

      secondterm(i) = secondterm(i) -(-m_mode*r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i))*(&
      & l*(l+1.)*Gamma_1(i)*P(i)*G_func(i)/omega(i) &
      & + Gamma_1(i)*P(i)/omega_R &
      & - omega(i)*rho(i)*G_func(i)*r(i)**2)!*(1./(8.*PI))

      thirdterm(i) = thirdterm(i) + (zeta2(i)*((xi_r(i))**2)*r(i)**2)*(&
      !& +G_func(i)*l*(l+1)*m_mode*omega_R*Gamma_1(i)*P(i)*zeta3(i)/omega(i) &
      !& - m_mode*omega_R*Gamma_1(i)*P(i)*(-zeta1(i) - 2./r(i) - A(i))/omega(i) &
      & + rho(i)*omega_R*m_mode*omega(i)*r(i)**2 *G_func(i)*zeta3(i) &
      & - rho(i)*omega_R*K1*dev1tot)!*(1./(8.*PI))

      !simplified eq.
      !thirdterm(i) = 2.*(m_mode**2)*omega_rot(i)*rho(i)*kr2(i)*&
      !&((xi_r(i))**2)*(r(i)**4)*G_func(i)*alph(i)*zeta3(i)


      coeff(i) = (thirdterm(i) + secondterm(i) +firstterm(i)) 

      omega(i) = omega_R + m_mode*omega_rot(i)
      um(i) = (r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i))*(&
      & - m_mode*l*(l+1.)*G_func(i) )
      !omega(i) = omega_R - m_mode*omega_rot(i)
      !um(i) = um(i) + (r(i)**2 *kr2(i)*((xi_r(i))**2)*zeta0(i))*(&
      !& + m_mode*l*(l+1.)*G_func(i) )

      omega(i) = omega_R + m_mode*omega_rot(i)
      dois(i) =  - Gamma_1(i) 
      !omega(i) = omega_R - m_mode*omega_rot(i)
      !dois(i) = dois(i)  + Gamma_1(i) 

      omega(i) = omega_R + m_mode*omega_rot(i)
      tres(i) = - P(i) 
      !omega(i) = omega_R - m_mode*omega_rot(i)
      !tres(i) = tres(i) + P(i) 

      omega(i) = omega_R + m_mode*omega_rot(i)
      quatro(i) =  - 1./(omega(i))
      omega(i) = omega_R - m_mode*omega_rot(i)
      quatro(i) = quatro(i) +1./(omega(i))





      !! ----------------------------------------------------------

    end do

    do i=1,m-1
   
      if ( k == 369 ) then
        write(20,*) r(i), coeff(i),firstterm(i), secondterm(i), thirdterm(i),&
        &um(i),dois(i),tres(i),quatro(i)
       
      end if
    end do
    close(20)
    
  end subroutine calculate_gwaves

  subroutine calculate_gwaves_simplified(array,m,summaryfile,k,coeff,MESAprofile)
    implicit none
    real (DP), intent(in) :: array(:,:), summaryfile(:,:),MESAprofile(:,:)
    real (DP), intent(inout) :: coeff(:)
    integer, intent(in) :: m,k
    real (DP) ,dimension(m) :: r, omega,Lr,As,c_1,N2
    real (DP) ,dimension(m) :: kr2,A
    real (DP) ,dimension(m) :: x,xi_r,Gamma_1,omega_rot,Pr,rho,V_2,cs2,alph,V
    real (DP) :: omega_R, l, m_mode,freq
    integer :: i

    l      = summaryfile(8,k)
    m_mode = summaryfile(9,k)
    freq   = summaryfile(5,k)
    omega_R = freq*10**(-6.)*2.*PI !rad/s
    Lr = getL(MESAprofile,array,m)
    omega_rot = getOmega(MESAprofile,array,m)

    if ( k == 369 ) then
      open(20,action='write',file='output/detail.l'//trim(string(int(l)))//'.j'//trim(string(k))//'.txt')
    end if
    
    do i=1,m
    !put the variables into arrays
      x(i)         = array(15,i)
      xi_r(i)      = array(18,i)
      As(i)        = array(1,i)
      c_1(i)       = array(8,i)
      Gamma_1(i)   = array(2,i)
      Pr(i)         = array(5,i)
      rho(i)       = array(14,i)
      V_2(i)       = array(7,i)
      !nabla(i)     = array(12,i)
      !nabla_ad(i)  = array(13,i)

    !put the functions into arrays
      r(i) = array(15,i)*R_star
      omega(i) = omega_R + m_mode*omega_rot(i)
      kr2(i) = computeKr2(array,i,l,freq)
      !K0 = compute_K0(m_mode,l)
      !K1 = compute_K1(m_mode,l)
      !K2 = compute_K1(m_mode,l)
      A(i) = computeA(array,i)
      !zeta0(i) = compute_zeta0(array,i,Lr(i))
      !zeta1(i) = compute_zeta1(array,i,l)
      !zeta2(i) = compute_zeta2(array,i,l,freq,Lr(i))
      !zeta3(i) = compute_zeta3(array,i,l)
      !G_func(i) = computeG(array,i,freq,l)

      !ds_dr(i) = - (compute_cp(array,i) *(nabla(i) - nabla_ad(i)) * V_2(i) *r(i))/(R_star**2)
      alph(i) = alpha(array,i,Lr(i))
      cs2(i) = compute_cs2(array,i)
      V(i) = -x(i)**2 * V_2(i)
      N2(i) = (G*M_star/(R_star**3))*(As(i)/c_1(i))

    end do

    do i=1,m

      coeff(i) = -(abs(m_mode**2)/(omega_R**2))*kr2(i)*abs(xi_r(i))**2 &
      &*omega_rot(i)*alph(i)*r(i)**2 *rho(i)*N2(i)/A(i)!*Pr(i)*V(i)/r(i)

      !coeff(i) = (abs(m_mode**2)/(omega_R**2))*kr2(i)*abs(xi_r(i))**2 &
      !&*omega_rot(i)*alph(i)*r(i)**2 *(A(i) - zeta3(i))*cs2(i)*rho(i)!*Gamma_1(i)*P(i)

      if ( k == 369 ) then
        write(20,*) r(i), coeff(i),rho(i)*N2(i)/A(i),&
        &Pr(i)*V(i)/r(i)
        !write(20,*) r(i), coeff(i),&
        !&(abs(m_mode**2)/(omega_R**2))*kr2(i)*abs(xi_r(i))**2 *omega_rot(i)*alph(i)*r(i)**2 *A(i)*cs2(i)*rho(i),&
        !&-(abs(m_mode**2)/(omega_R**2))*kr2(i)*abs(xi_r(i))**2 *omega_rot(i)*alph(i)*r(i)**2 *zeta3(i)**cs2(i)*rho(i)
      end if

    end do
    close(20)
    
  end subroutine calculate_gwaves_simplified

  subroutine calculate_coefficients(array,m,summaryfile,k,coeff,MESAprofile)
    implicit none
    real (DP), intent(in) :: array(:,:), summaryfile(:,:),MESAprofile(:,:)
    real (DP), intent(inout) :: coeff(:)
    integer, intent(in) :: m,k
    real (DP) ,dimension(m) :: A_ml , B_ml, C_ml, D_ml, dr,r, omega,cs2,Lr
    real (DP) ,dimension(m) :: zeta0,zeta1,zeta2,zeta3,alph,kr2,A,G_func,N2
    real (DP) ,dimension(m) :: x,xi_r,xi_h,V_2,As,c_1,Gamma_1,nabla,nabla_ad,delta,c_thk,omega_rot,P,rho,T
    real (DP) :: omega_R, G1, K0,K1,K2, l, m_mode,freq
    real (DP) :: dev01,dev02,devtot,devtot1,dzeta0,dzeta1,dzeta2,dzeta3,dev2,dev1,dev2tot,dev1tot
    integer :: i

    dev01 = 0.
    dev02 = 0.
    devtot = 0.
    devtot1 = 0.
    dzeta0=0.
    dzeta1=0.
    dzeta2=0.
    dzeta3=0.
    dev1 = 0.
    dev2 = 0.
    dev1tot = 0.
    dev2tot = 0.

    l      = summaryfile(8,k)
    m_mode = summaryfile(9,k)
    freq   = summaryfile(5,k)
    omega_R = freq*10**(-6.)*2.*PI !rad/s
    Lr = getL(MESAprofile,array,m)
    omega_rot = getOmega(MESAprofile,array,m)

    !open(20,action='write',file='output/detail.l'//trim(string(int(l)))//'.j'//trim(string(k))//'.txt')

    
    do i=1,m
    !put the variables into arrays
      x(i)         = array(15,i)
      xi_r(i)      = array(18,i)!*R_star/abs(array(2,m))
      xi_h(i)      = array(16,i)!*R_star/abs(array(2,m)) 
      V_2(i)       = array(7,i)
      As(i)        = array(1,i)
      c_1(i)       = array(8,i)
      Gamma_1(i)   = array(2,i)
      nabla(i)     = array(12,i)
      nabla_ad(i)  = array(13,i)
      delta(i)     = array(10,i)
      c_thk(i)     = array(9,i)
      !omega_rot(i) = array(4,i)*sqrt((G*M_star)/(R_star**3))
      P(i)         = array(5,i)
      rho(i)       = array(14,i)
      T(i)         = array(6,i)

    !put the functions into arrays
      r(i) = array(15,i)*R_star
      omega(i) = omega_R + m_mode*omega_rot(i)
      kr2(i) = computeKr2(array,i,l,freq)
      alph(i) = alpha(array,i,Lr(i))
      A(i) = computeA(array,i)
      K0 = compute_K0(m_mode,l)
      K1 = compute_K1(m_mode,l)
      K2 = compute_K1(m_mode,l)
      zeta0(i) = compute_zeta0(array,i,Lr(i))
      zeta1(i) = compute_zeta1(array,i,l)
      zeta2(i) = compute_zeta2(array,i,l,freq,Lr(i))
      zeta3(i) = compute_zeta3(array,i,l)
      G_func(i) = computeG(array,i,freq,l)
      cs2(i) = compute_cs2(array,i)
      N2(i) = (G*M_star/(R_star**3))*(As(i)/c_1(i)) !(rad/s)^2

    end do

    !compute A_ml
    do i=2,m-2
      K1 = compute_K1(m_mode,l)
      A_ml(i) = (1./(8.*PI)) * rho(i) * alph(i) * kr2(i) * ((xi_r(i))**2) *(K1 + compute_K1(-m_mode,l))
    end do


    !compute B_ml
    do i=2,m-2
      
      dr(i) = r(i+1) - r(i)
     
      dev01 = derive(alph(i), alph(i+1),dr(i)) 
      dev02 = derive(kr2(i), kr2(i+1),dr(i))
      devtot = 2./r(i) + (1./alph(i))*dev01 + (1./kr2(i))*dev02
      
      !positive m
      omega(i) = omega_R + m_mode*omega_rot(i)
      K0 = compute_K0(m_mode,l)
      K1 = compute_K1(m_mode,l)
      B_ml(i) = (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI) * (((alph(i)*K1)/2.)*(A(i) + 2.*zeta1(i) + devtot) + &
              (1./(2.*cs2(i)**2))*(r(i)**2 * omega_R * omega(i) * G_func(i) * zeta0(i) * K1) +&
                ((m_mode**2) * alph(i) * G_func(i) * zeta3(i))/2. - &
               (omega_R/omega(i)) * G_func(i) * zeta0(i) * ((m_mode**2) + K0/2.))

      !negative m
      omega(i) = omega_R - m_mode*omega_rot(i)
      K0 = compute_K0(-m_mode,l)
      K1 = compute_K1(-m_mode,l)
      B_ml(i) = B_ml(i) + (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI) * (((alph(i)*K1)/2.)*(A(i) + 2.*zeta1(i) + devtot) + &
              (1./(2.*cs2(i)**2))*(r(i)**2 * omega_R * omega(i) * G_func(i) * zeta0(i) * K1) +&
                ((m_mode**2) * alph(i) * G_func(i) * zeta3(i))/2. - &
               (omega_R/omega(i)) * G_func(i) * zeta0(i) * ((m_mode**2) + K0/2.))         
               
    end do


    !compute C_ml
    do i=2,m-2
      dzeta0 = derive(zeta0(i), zeta0(i+1),dr(i))
      dev01 = derive(G_func(i), G_func(i+1),dr(i))
      dev02 = derive(kr2(i), kr2(i+1),dr(i))
      devtot = 1./r(i) + (1./G_func(i))*dev01 + (1./kr2(i))*dev02 + (1./zeta0(i))*dzeta0

      !positive m
      omega(i) = omega_R + m_mode*omega_rot(i)
      K2 = compute_K1(m_mode,l)
      C_ml(i) = (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI)*&
      ( (-omega_R/omega(i))*(1./Gamma_1(i))*((-V_2(i)*x(i)**2)/r(i))*G_func(i)*K2*(zeta0(i) - alph(i)*zeta3(i))&
      + (omega_R/omega(i))*G_func(i)*zeta0(i)*K2 *(((A(i)*r(i) - 1.)/r(i)) + devtot + zeta1(i) ) )

      !negative m
      omega(i) = omega_R - m_mode*omega_rot(i)
      K2 = compute_K1(-m_mode,l)
      C_ml(i) = C_ml(i) + (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI)*&
      ( (-omega_R/omega(i))*(1./Gamma_1(i))*((-V_2(i)*x(i)**2)/r(i))*G_func(i)*K2*(zeta0(i) - alph(i)*zeta3(i))&
      + (omega_R/omega(i))*G_func(i)*zeta0(i)*K2 *(((A(i)*r(i) - 1.)/r(i)) + devtot + zeta1(i) ) )
    end do

    if ( k == 369 ) then
      open(20,action='write',file='output/detail.l'//trim(string(int(l)))//'.j'//trim(string(k))//'.txt')
    end if

    !compute D_ml
    do i=2,m-2
      
      G1 = -(1./(l*(l+1.)))
      dev02 = derive(kr2(i), kr2(i+1),dr(i))
      dzeta0 = derive(zeta0(i), zeta0(i+1),dr(i))
      dzeta1 = l*(l+1.)*(-(xi_h(i))/(r(i)**2 * xi_r(i)) + derive(xi_h(i),xi_h(i+1),dr(i))/(r(i)*xi_r(i)) -&
               (xi_h(i)*derive(xi_r(i),xi_r(i+1),dr(i)))/(r(i)*xi_r(i)**2))
      dzeta2 = derive(zeta2(i), zeta2(i+1),dr(i))
      dzeta3 = -2./r(i)**2 + dzeta1 - &
      (1./R_star**2)*derive((r(i) *V_2(i))/(gamma_1(i)), (r(i+1) *V_2(i+1))/(gamma_1(i+1)),dr(i))

      devtot = zeta0(i)/r(i) + (zeta0(i)/kr2(i))*dev02 + dzeta0
      devtot1 = zeta3(i)/r(i) + dzeta3

      !positive m
      omega(i) = omega_R + m_mode*omega_rot(i)
      D_ml(i) = m_mode*omega(i)*(&
      (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI))*((alph(i)/2.)*(1.-N2(i)/(omega(i)**2)) &
      -((r(i)**2 *G1)/2.)*(zeta0(i)*A(i) + zeta0(i)/r(i) + zeta0(i)*zeta1(i) + devtot &
      -(-zeta0(i)*r(i) *V_2(i))/(gamma_1(i)*R_star**2)) + ((r(i)**2 *G1*alph(i))/2.)*&
      (zeta3(i)*A(i) + (3.*zeta3(i))/r(i) + 2.*zeta1(i)*zeta3(i) +&
      devtot1 + (zeta3(i)/zeta2(i))*dzeta2) )

      !negative m
      omega(i) = omega_R - m_mode*omega_rot(i)
      D_ml(i) = D_ml(i) - m_mode*omega(i)*(&
      (rho(i)*kr2(i)*xi_r(i)**2)/(4.*PI))*((alph(i)/2.)*(1.-N2(i)/(omega(i)**2)) &
      -((r(i)**2 *G1)/2.)*(zeta0(i)*A(i) + zeta0(i)/r(i) + zeta0(i)*zeta1(i) + devtot &
      -(-zeta0(i)*r(i) *V_2(i))/(gamma_1(i)*R_star**2)) + ((r(i)**2 *G1*alph(i))/2.)*&
      (zeta3(i)*A(i) + (3.*zeta3(i))/r(i) + 2.*zeta1(i)*zeta3(i) +&
      devtot1 + (zeta3(i)/zeta2(i))*dzeta2) )
      

      D_ml(i) = 0.
      D_ml(i) = m_mode*((alph(i)*rho(i)*kr2(i)*xi_r(i)**2)/2.)*(&
      -N2(i)/(omega_R + m_mode*omega_rot(i))&
      &+N2(i)/(omega_R - m_mode*omega_rot(i)) )!&
      !&+ (omega_R + m_mode*omega_rot(i))*(&
      !  &((r(i)**2 *G1))*( 2.*zeta1(i)*zeta3(i) +devtot1) )&
      !  &-(omega_R - m_mode*omega_rot(i))*(&
      !  &((r(i)**2 *G1))*( 2.*zeta1(i)*zeta3(i) +devtot1 ) ))

      if ( k == 369 ) then
        write(20,*) r(i),D_ml(i),&
        &m_mode*((alph(i)*rho(i)*kr2(i)*xi_r(i)**2)/(2.))*(&
        &-N2(i)/(omega_R + m_mode*omega_rot(i)) + N2(i)/(omega_R - m_mode*omega_rot(i))),&
        &((m_mode*rho(i)*kr2(i)*xi_r(i)**2))*((r(i)**2 *G1*alph(i))/2.)*(2.*zeta1(i)*zeta3(i))*(&
        &(omega_R + m_mode*omega_rot(i))-(omega_R - m_mode*omega_rot(i))),&
        &((m_mode*rho(i)*kr2(i)*xi_r(i)**2))*((r(i)**2 *G1*alph(i))/2.)*(devtot1)*(&
        &(omega_R + m_mode*omega_rot(i))-(omega_R - m_mode*omega_rot(i)))

     
      end if

    end do
    close(20)


  !compute total AM transport
    do i=2,m-2
      dev2 = derive_2ndorder(omega_rot(i),omega_rot(i+1),omega_rot(i+2),dr(i))
      dev1 = derive(omega_rot(i),omega_rot(i+1),dr(i))
      dev2tot = 2.*omega_rot(i) + 4.*r(i)*dev1 + (r(i)**2)*dev2
      dev1tot = 2.*r(i)*omega_rot(i) + (r(i)**2)*dev1

      !coeff(i) = A_ml(i)*dev2tot + B_ml(i)*dev1tot + C_ml(i)*(r(i)**2 *omega_rot(i)) + D_ml(i)

      coeff(i) = D_ml(i)
      !coeff(i) = -(m_mode/2.)*((alph(i)*rho(i)*kr2(i)*N2(i)*xi_r(i)**2))*(&
      !&1./(omega_R + m_mode*omega_rot(i)) - 1./(omega_R - m_mode*omega_rot(i)))

     !write(20,*) x(i), A_ml(i)*dev2tot, B_ml(i)*dev1tot, C_ml(i)*(r(i)**2 *omega_rot(i)), m_mode*omega(i)*D_ml(i)
    end do
    
    !close(20)
    
  end subroutine calculate_coefficients



  subroutine calculate_simplifiedeq(array,summaryfile,j,coeff,MESAprofile)
    implicit none
    real (DP), intent(in) :: array(:,:),summaryfile(:,:),MESAprofile(:,:)
    real (DP), intent(inout) :: coeff(:)
    real (DP) :: r, alph,kr2,N2,xi_r,As,c_1,rho,omega_R,l,m_mode,freq
    integer, intent(in) :: j
    integer :: i,m
    real (DP) ::V_2,nabla,nabla_ad,T,ds_dr
    real (DP), allocatable :: Lr(:),omega_rot(:)

    m = size(array,dim=2)
    allocate(Lr(m),omega_rot(m))
 
    l      = summaryfile(8,j)
    m_mode = summaryfile(9,j)
    freq   = summaryfile(5,j)
    omega_R = freq*10**(-6.)*2.*PI !rad/s
    Lr = getL(MESAprofile,array,m)
    omega_rot = getOmega(MESAprofile,array,m)
    
    if ( j == 369 ) then
      open(305,action='write',file='output/detail.l'//trim(string(int(l)))//'.j'//trim(string(j))//'.txt')
    end if
    do i=1,m
    !put the variables into variables
      xi_r      = array(18,i)
      V_2       = array(7,i)
      As        = array(1,i)
      c_1       = array(8,i)
      nabla     = array(12,i)
      nabla_ad  = array(13,i)
      rho       = array(14,i)
      T         = array(6,i)
      !omega_rot = array(4,i)*sqrt((G*M_star)/(R_star**3))

    !put the functions into variables
      r = array(15,i)*R_star
      ds_dr = - (compute_cp(array,i) *(nabla - nabla_ad) * V_2 *r)/(R_star**2)
      kr2 = computeKr2(array,i,l,freq)
      alph = alpha(array,i,Lr(i))
      N2 = (G*M_star/(R_star**3))*(As/c_1) !(1/s)^2

      !S_l2 = ((compute_cs2(array,i)* l*(l+1.))/r**2 )!*(G*M_star/(R_star**3))


      !coeff(i) = (2.*abs(m_mode**2)/(omega_R**2)) *rho*kr2*abs(xi_r)**2 *omega_rot(i)*alph*N2
      coeff(i) = 2.*(abs(m_mode**2)/(omega_R**2)) *rho*kr2*abs(xi_r)**2 *alph*N2!*omega_rot(i)
    
      if (j == 369) then
        write (305,*) r, coeff(i), omega_rot(i),&
        (2.*abs(m_mode**2)/(omega_R**2)),rho,kr2,abs(xi_r)**2,alph,N2!,compute_cs2(array,i)
        !&L_star,T,nabla,nabla_ad,compute_cp(array,i),V_2,((L_star)/(4.*PI*r**2 *rho*T)),(nabla_ad/nabla -1.),1./ds_dr
        !&(omega_R**2/compute_cs2(array,i)),(1. - S_l2/(omega_R**2)),(1. - N2/(omega_R**2)),S_l2,omega_R**2
      end if

    end do
    !close(20)
    close(305)


    deallocate(Lr,omega_rot)

  end subroutine calculate_simplifiedeq



  function getL(MESAprofile,detailfile,m) result(Lr)
    implicit none
    real(DP), intent(in) :: MESAprofile(:,:),detailfile(:,:)
    real(DP) :: Lr(m)
    integer, intent(in) :: m
    integer :: N

    N = size(MESAprofile,dim=2)

    Lr = interpolate(MESAprofile(2,:),MESAprofile(4,:),N,detailfile(15,:)*R_star,m)
  
  end function getL


  function getOmega(MESAprofile,detailfile,m) result(omega)
    implicit none
    real(DP), intent(in) :: MESAprofile(:,:),detailfile(:,:)
    real(DP) :: omega(m)
    integer, intent(in) :: m
    integer :: N

    N = size(MESAprofile,dim=2)

    omega = interpolate(MESAprofile(2,:),MESAprofile(19,:),N,detailfile(15,:)*R_star,m)
    
  end function getOmega

  function compute_cp(array,j) result(cp)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in) :: array(:,:)
    real (DP) :: rho,T,c_thk, cp
    
    c_thk = array(9,j)
    rho   = array(14,j)
    T     = array(6,j)

    cp = ((c_thk * L_star)/(4. * PI * R_star**3 *T*rho)) * (sqrt((R_star**3)/(G*M_star)))
  end function compute_cp


  function alpha(array,j,Lr) result(a)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in):: array(:,:),Lr
    real (DP) :: x,rho,T,nabla_ad,nabla,ds_dr,r,V_2, a

    x        = array(15,j)
    V_2      = array(7,j)
    nabla    = array(12,j)
    nabla_ad = array(13,j)
    rho      = array(14,j)
    T        = array(6,j)
    
    r = x*R_star

    !Kevin Belkacem code
    !ds_dr = cp*(dlnT_dr-gradad*dlnP_dr)

    ds_dr = - (compute_cp(array,j) *(nabla - nabla_ad) * V_2 *r)/(R_star**2)
    
    a =  -(Lr/(4.*PI*r**2 *rho*T)) * ((nabla_ad/nabla) - 1.) /ds_dr
    
  end function alpha

  function computeKr2(array,j,l,freq) result(k)
    implicit none
    integer, intent(in) :: j
    real (DP) :: As,c_1,x,r,omega2,N2,k,S_l2,cs2
    real(DP), intent(in) :: l, freq, array(:,:)
 
    x   = array(15,j)
    As  = array(1,j)
    c_1 = array(8,j)

    r = x*R_star
    omega2 = (freq*10**(-6.)*2.*PI)**2
    cs2 = compute_cs2(array,j)
    N2 = (As/c_1) *(G*M_star/(R_star**3))!rad/s
    S_l2 = ((cs2* l*(l+1.))/r**2 )!*(G*M_star/(R_star**3))
   
    k= omega2/cs2 *(1. - S_l2/omega2)*(1. - N2/omega2)
    
  end function computeKr2

  function compute_cs2(array,j) result(cs)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in) :: array(:,:)
    real(DP) :: gamma_1, V_2, c_1, cs

    V_2     = array(7,j)
    c_1     = array(8,j)
    gamma_1 = array(2,j)
  
    cs = (gamma_1/(V_2*c_1))*((G*M_star)/R_star)
  end function compute_cs2

  function derive(f,f1,dr) result(df)
    implicit none
    real (DP), intent(in) :: f, f1, dr
    real (DP):: df

    df = (f1-f)/dr
  end function derive

  function derive_2ndorder(f,f1,f2,dr) result(df)
    implicit none
    real (DP), intent(in) :: f, f1, f2, dr
    real (DP) :: df

    df = (f2 - 2.*f1 + f)/(dr**2)
  end function derive_2ndorder

  function derive_moreprecision(f,nlines) result(df)
    implicit none
    real(DP), intent(in) :: f(:,:)
    integer, intent(in) :: nlines
    integer :: i
    real (DP) :: df(nlines),dr

    !O(∆x2) forward difference approximation
    do i=1,2
      dr = f(1,i+1) - f(1,i)
      df(i) = (-3.*f(2,i) + 4.*f(2,i+1) -f(2,i+2))/(2.*dr)
      df(i) = -df(i)/(2.*f(1,i)**2)
    end do

    !O(∆x2) backward difference approximation
    do i=nlines-1,nlines
      dr = f(1,i) - f(1,i-1)
      df(i) = (3.*f(2,i) - 4.*f(2,i-1) + f(2,i-2))/(2.*dr)
      df(i) = -df(i)/(2.*f(1,i)**2)
    end do

    !O(∆x4) centered difference approximation
    do i=3,nlines-2
      dr = f(1,i+1) - f(1,i)
      df(i) = (-f(2,i+2) + 8.*f(2,i+1) - 8.*f(2,i-1) + f(2,i-2))/(12.*dr)
      df(i) = -df(i)/(2.*f(1,i)**2)
    end do

  end function

  function compute_J(m,l) result(J)
    implicit none
    real(DP), intent(in) :: m,l
    real (DP) :: J

    if (l>abs(m)) then
      J = sqrt(((l**2 - m**2))/(4.* (l**2) - 1.))
    else
      J = 0.
    end if

  end function compute_J

  function compute_K0(m,l) result(k)
    implicit none
    real (DP), intent(in) :: m,l
    real (DP) :: k

    k = (l**2) * (compute_J(m,l+1.))**2 + ((l+1.)**2)*((compute_J(m,l))**2)

  end function compute_K0

  function compute_K1(m,l) result(k)
    implicit none
    real (DP), intent(in) :: m,l
    real (DP):: k

    k = 1. -(compute_J(m,l+1.))**2 - ((compute_J(m,l))**2)

  end function compute_K1

  function compute_K2(m,l) result(k)
    implicit none
    real (DP), intent(in) :: m,l
    real (DP) :: k

    k = l*(compute_J(m,l+1.))**2 - (l+1)*((compute_J(m,l))**2)

  end function compute_K2

  function compute_zeta0(array,j,Lr) result(z)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in) :: array(:,:),Lr
    real (DP) :: x,rho,T,nabla_ad,nabla,r,rhoT,c_thk,cp, z

    x        = array(15,j)
    nabla    = array(12,j)
    nabla_ad = array(13,j)
    rhoT     = -array(10,j)
    c_thk    = array(9,j)
    rho      = array(14,j)
    T        = array(6,j)

    r = x*R_star
    cp = compute_cp(array,j)

    z = - (Lr/(4.*PI*r**2 * rho * T))*((nabla_ad/nabla) - 1d0)*(rhoT/cp)

  end function compute_zeta0

  function compute_zeta1(array,j,l) result(z)
    implicit none
    integer, intent(in) :: j
    real(DP), intent(in) :: l, array(:,:)
    real (DP) :: x,xi_r,xi_h,r, z
    integer :: m
    m = size(array,dim=2)

    x    = array(15,j)
    xi_r = array(18,j)*R_star!/abs(array(18,m))
    xi_h = array(16,j)*R_star!/abs(array(18,m))

    r = x*R_star

    z = ((l*(l+1.))/r) * (xi_h/xi_r)
  end function compute_zeta1

  function compute_zeta2(array,j,l,freq,Lr) result(z)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in) :: l, freq, array(:,:),Lr
    real (DP) :: omega, z
    
    omega = freq*10**(-6.)*2.*PI

    z = (alpha(array,j,Lr)* computeKr2(array,j,l,freq)) /omega
  end function compute_zeta2

  function compute_zeta3(array,j,l) result(z)
    implicit none
    integer, intent(in) :: j
    real(DP), intent(in) :: l, array(:,:)
    real (DP) :: r,x,gamma_1,V2,V, z

    x       = array(15,j)
    V2      = array(7,j)
    gamma_1 = array(2,j)
    
    r = x*R_star
    V = x**2 * V2

    z = 2./r + compute_zeta1(array,j,l) - (1./gamma_1)*V/r
  end function compute_zeta3

  function computeA(array,j) result(a)
    implicit none
    integer, intent(in) :: j
    real (DP), intent(in) :: array(:,:)
    real (DP) :: x,A_star,r,a

    x      = array(15,j)
    A_star = array(1,j)

    r = x*R_star
    a = -A_star/r
  end function computeA


  function computeG(array,j,freq,l) result(g)
    implicit none
    integer, intent(in) :: j
    real(DP), intent(in) :: freq, l, array(:,:)
    real (DP) :: x,omega,S_l,r,cs2, g

    x = array(15,j)

    r = x*R_star
    cs2 = compute_cs2(array,j)
    S_l= (((l*(l+1.)))*cs2) / r**2
    omega = freq*10**(-6.)*2.*PI
    
    g = ((cs2)/(r**2 * omega**2))*(1./(1. - S_l/(omega**2)))
  end function computeG


  END MODULE coefficient_equations