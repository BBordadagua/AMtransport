                                                                                                                                                                                          !!
!!                 parameters.f90
!!
!!     This program contains parameters to be used
!!           on the other two programs.
!!
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!

MODULE parameters
  IMPLICIT NONE

  public :: set_parameters!,getlmax

  ! numerical precision
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP  = KIND(1.0)
  INTEGER, PARAMETER :: DP  = KIND(1.d0)

  !Constants
  REAL (DP) :: G = 6.67259000d-08 !!cm^3 g^-1 s^-2
  REAL (DP) :: PI = 3.14159265358979d0
  REAL (DP) :: c_speed = 2.99792458d10 !!cm/s
  REAL (DP) :: a_raddensity = 7.5646d-15 !!erg cm^-3 K^-4 !!radiation density constant

  !Solar values
  REAl (DP) :: R_sun = 6.95700d10 !cm
  REAl (DP) :: M_sun = 1.98847d33 !g
  REAl (DP) :: L_sun = 3.846d33   !erg/s
  REAl (DP) :: T_eff_sun = 5777d0   !K
  REAl (DP) :: nu_max_sun = 3050d0  !μHz
  REAl (DP) :: delta_nu_sun = 135d0 !μHz
  REAL (DP) :: A_numax_sun = 20d0!18.7 !+- 0.7 cm/s Kjeldsen et al.2008
  !Vsun_max = 20.  # solar maximum at nu_max [cm/s] !Kevin Belkacem code

  !General properties of the model
  REAL (DP) :: M_star         !g
  REAL (DP) :: R_star         !cm
  REAL (DP) :: L_star         !erg/s
  REAL (DP) :: T_eff_star     !K 
  REAL (DP) :: nu_max_star    !μHz 
  REAL (DP) :: delta_nu_star  !μHz 
  REAl (DP) :: dt             !s

  !Variables dependent on the evolutionary stage
  !REAL (DP) :: delta_env = 70 !M1  !30 !μHz M2 Belkacem 2015
  INTEGER :: lmax = 2 !input parameter

  contains

  subroutine set_parameters(MESAfile,model) !!MESA version 23.05.1
    implicit none
    integer, intent(in) :: model
    real (DP), intent(in)::  MESAfile(:,:)

    M_star        = MESAfile(5,model)*M_sun
    T_eff_star    = 10**MESAfile(37,model)
    L_star        = (10**MESAfile(38,model))*L_sun
    R_star        = (10**MESAfile(39,model))*R_sun
    nu_max_star   = MESAfile(57,model)
    delta_nu_star = MESAfile(56,model)
    dt            = (10**MESAfile(4,model))*365d0*24d0*60d0*60d0 !s

    !lmax = getlmax(summaryfile,MESAfile,model)

  end subroutine

    
  !function getlmax(summaryfile,MESAfile,model) result(l_max)
  !  implicit none
  !  integer, intent(in) :: model
  !  real (DP), intent(in):: summaryfile(:,:), MESAfile(:,:)
  !  integer :: l_max

  !  l_max = 2

  !end function getlmax

END MODULE parameters
