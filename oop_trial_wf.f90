module trial_wf_pot
  implicit real(8) (a-h,o-z)
  integer, parameter :: dim = 24 !number of dimensions
  !real(8), dimension(dim) :: k_pot = 1.2 !parameters for the harmonic potential
  real(8), dimension(dim) :: k_pot = (/7.4752000e-01,&
     6.0555556e-01,&
     7.8544000e-01,&
     8.1333333e-01,&
     6.1551020e-01,&
     7.8544000e-01,&
     5.0161111e-01,&
     5.8110727e-02,&
     8.9668639e-02,&
     9.7147929e-02,&
     3.5093750e-01,&
     8.6785714e-02,&
     9.7147929e-02,&
     3.2971429e-01,&
     8.6785714e-02,&
     5.0511111e-01,&
     5.8110727e-02,&
     8.9668639e-02,&
     9.7147929e-02,&
     3.5093750e-01,&
     8.9680473e-02,&
     9.7147929e-02,&
     3.2971429e-01,&
     8.9680473e-02/)
  real(8), dimension(dim) :: masses = (/12.011, 12.011, 12.011, 12.011, 12.011, 12.011,&
  1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007, 1.007,&
   1.007, 1.007, 1.007, 1.007, 1.007, 1.007/)
  !real(8), dimension(dim) :: masses = 1.0
  !real(8), dimension(dim) :: sigma2_shift =(/0.9143, -0.0292, 0.6006, -0.7162, &
  !-0.1565, 0.8315, 0.5844, 0.9190, 0.3115, -0.9286, 0.6983, 0.8680/)
  real(8), dimension(dim) :: mu0 = 0.0
  real(8), dimension(dim) :: sigma2 = 1.0 !sigma^2
  !real(8), dimension(dim) :: mu_shift =(/0.3575, 0.5155, 0.4863, -0.2155,&
  !0.3110, -0.6576, 0.4121, -0.9363, -0.4462, -0.9077, -0.8057, 0.6469/)
  !real(8) :: shift_mult = 0.001

contains

!calculate the drift velocity grad(psi)/psi-------------------------------------
function drift(rxyz) result(v)
  implicit real(8) (a-h,o-z)
  real(8), dimension(3,8), intent(in) :: rxyz
  real(8), dimension(dim) :: v, r, mu

  sigma2 = 1/sqrt(k_pot*masses)
  !sigma2 = sigma2 + shift_mult*sigma2_shift
  !mu = mu0 + shift_mult*mu_shift
  mu = mu0

  r = reshape(rxyz,shape(r))

  do i=1,dim
    !v(i) = (psi(rxyz+h,k_pot,dim) - psi(rxyz-h,k_pot,dim))/(2*h_step)/psi(rxyz,k_pot,dim)
    v(i) = -(r(i)-mu(i))/sigma2(i)
  end do

end function drift

!calculate the local energy psi*H/psi-------------------------------------------
real(8) function el(rxyz)
  implicit real(8) (a-h,o-z)
  real(8), dimension(dim), intent(in) :: rxyz
  real(8) ,dimension(dim) :: r, mu

  sigma2 = 1/sqrt(k_pot*masses)
  !sigma2 = sigma2 + shift_mult*sigma2_shift
  !mu = mu0 + shift_mult*mu_shift
  mu = mu0

  r = reshape(rxyz,shape(r))

  el = 0
  do i=1,dim
    el = el + 0.5*k_pot(i)*(r(i)-mu(i))**2 - 0.5/sigma2(i)*((r(i)-mu(i))**2/sigma2(i) - 1)/masses(i)
    !el = el + 0.5*k_pot(i)*rxyz(i)**2
  end do

end function el

subroutine calc_e0_theory(e0_theory)
  implicit real(8) (a-h,o-z)

  e0_theory = 0
  do i=1,dim
    e0_theory = e0_theory + 0.5*sqrt(k_pot(i)/masses(i))
  end do

end subroutine calc_e0_theory

end module trial_wf_pot
