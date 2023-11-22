!!
!!                 modeamplitude_mod.f90
!!
!!  This program contains the module 'compute_modeamplitude' with functions
!!      and subroutines that calculate the mode amplitudes.
!!
!!  compile with: > make
!!  usage:        > ./program.exec
!!  clean exec:   > make clean
!!

MODULE compute_modeamplitude
  use parameters
  use data_handling
  use math_functions
  implicit none
  
  !define public subroutines and functions
  public :: computel0,compute_amplitude
  
  contains

  subroutine compute_amplitude(array,summaryfile,j,amp_ml)
    implicit none
    real (DP), intent(in) :: array(:,:), summaryfile(:,:)
    real (DP), intent(inout) :: amp_ml
    real (DP) ::  V2,omega_R, axi2,freq,Enorm0,Enorml,inertia_ratio
    integer, intent(in) :: j
    integer :: m

    freq = summaryfile(5,j)
    Enorm0 = summaryfile(1,getj_l0(summaryfile,freq))
    Enorml = summaryfile(1,j)
    inertia_ratio = Enorm0/Enorml

    m = size(array,dim=2)
    V2 = computel0(freq)*inertia_ratio
    omega_R = freq*10**(-6.)*2.*PI

    axi2 = 2.*V2/(omega_R**2)
    amp_ml = axi2/(abs(array(18,m))**2)

    !# computation of raidative damping    !Kevin Belkacem code
    !omega = 2.*pi*final
    !eta_rad=zeros((len(omega)))
    !eta_sun = 0.9*1e-6*pi
    !for i in range(0,len(omega),1):
	    !eta_rad[i] = (ell*(ell+1.))**(1.5)/8./pi/omega[i]**3/theta_g[i]
	    !tmp_eta = (gradad/grad-1.)*gradad*BV*g*L/P/r**5
	    !eta_rad[i] = eta_rad[i] * trapz(tmp_eta[1:i_cut],r[1:i_cut])
	    !eta_p = eta_sun * (T_eff/Tsun)**(12.)


    !! print amplitude info
    !write(50,*) summaryfile(8,j),freq, sqrt(V2),summaryfile(5,getj_l0(summaryfile,freq)),&
    !& sqrt(computel0(summaryfile(5,getj_l0(summaryfile,freq)))), amp_ml,m!,abs(array(2,5))
    

  end subroutine


  function computel0(freq) result(V_final2)
    implicit none
    real (DP), intent(in) :: freq
    real (DP) :: V_final2,Vmax, sigma , delta_env

    Vmax = A_numax_sun * (T_eff_star/T_eff_sun)**(-0.77) * (nu_max_star/nu_max_sun)**(-1.15) * (delta_nu_star/delta_nu_sun)**(0.65)
    !Vmax = A_numax_sun * (T_eff_star/T_eff_sun)**(-0.77) * (nu_max_star/nu_max_sun)**(-1.15) * (delta_nu_star/delta_nu_sun)**(0.55)
    !Vmax = 300.

    !computation of radial mode amplitudes (Gaussian enveloppe)
    delta_env = (2.37+0.77*M_star/M_sun) * delta_nu_star   !! Eq.(4) Mosser et al. (2012) !RGB
    !delta_env = 30.
    sigma = delta_env/(2.*sqrt(2.*log(2.)))


    V_final2 = Vmax**2 * exp(-((freq-nu_max_star)**(2)) / (2.*sigma**2))

  end function computel0



  END MODULE compute_modeamplitude