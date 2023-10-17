Module meridionalcirculation_mod
  use parameters
  use data_handling
  IMPLICIT NONE
  
  !define public subroutines and functions
  public :: 

  contains

  subroutine compute_U2(MESAprofile,m,U2)
    implicit none
    real (DP), intent(in) :: MESAprofile(:,:)
    real (DP), intent(inout) :: U2(:)
    real (DP) ,dimension(m) :: P,cp,rho,T,g,nabla_ad,nabla,phi,delta,gradmu,gs,Lr,Eomega,Emu,theta
    integer :: i

  end subroutine compute_U2

  subroutine compute_Eomega(MESAprofile)
    implicit none
  end subroutine compute_Eomega

  subroutine compute_Emu(MESAprofile)
    implicit none
  end subroutine compute_Emu

  subroutine compute_theta(MESAprofile)
    implicit none
  end subroutine compute_theta

  !stationary case: eq. 4.40 Maeder & Zahn 1998
  subroutine compute_lambda(MESAprofile) 
    implicit none
  end subroutine compute_lambda

  subroutine compute_delta(MESAprofile) 
    implicit none
  end subroutine compute_delta

  subroutine compute_phi(MESAprofile) 
    implicit none
  end subroutine compute_phi

  

end module meridionalcirculation_mod
