module trial_wf_pot
  implicit real(8) (a-h,o-z)
  integer, parameter :: dim = 12 !number of dimensions
  !real(8), dimension(dim) :: k_pot = 1.2 !parameters for the harmonic potential
  real(8), dimension(dim) :: k_pot = (/2.0, 1.1, 0.4, 1.8, 3.2, 1.0, 3.4, 1.4, 1.1, &
  0.6, 1.9, 0.6/)
  real(8), dimension(dim) :: sigma2_shift =(/0.9143, -0.0292, 0.6006, -0.7162, &
  -0.1565, 0.8315, 0.5844, 0.9190, 0.3115, -0.9286, 0.6983, 0.8680/)
  real(8), dimension(dim) :: mu0 = -12.0
  real(8), dimension(dim) :: sigma2 = 1.0 !sigma^2
  real(8), dimension(dim) :: mu_shift =(/0.3575, 0.5155, 0.4863, -0.2155,&
  0.3110, -0.6576, 0.4121, -0.9363, -0.4462, -0.9077, -0.8057, 0.6469/)
  real(8) :: shift_mult = 0.001

contains

!calculate the drift velocity grad(psi)/psi-------------------------------------
function drift(rxyz) result(v)
  implicit real(8) (a-h,o-z)
  real(8), dimension(3,8), intent(in) :: rxyz
  real(8), dimension(dim) :: v, r, mu

  sigma2 = 1/sqrt(k_pot)
  sigma2 = sigma2 + shift_mult*sigma2_shift
  mu = mu0 + shift_mult*mu_shift

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

  sigma2 = 1/sqrt(k_pot)
  sigma2 = sigma2 + shift_mult*sigma2_shift
  mu = mu0 + shift_mult*mu_shift

  r = reshape(rxyz,shape(r))

  el = 0
  do i=1,dim
    el = el + 0.5*k_pot(i)*(r(i)-mu(i))**2 - 0.5/sigma2(i)*((r(i)-mu(i))**2/sigma2(i) - 1)
    !el = el + 0.5*k_pot(i)*rxyz(i)**2
  end do

end function el

subroutine calc_e0_theory(e0_theory)
  implicit real(8) (a-h,o-z)

  e0_theory = 0
  do i=1,dim
    e0_theory = e0_theory + 0.5*sqrt(k_pot(i))
  end do

end subroutine calc_e0_theory

end module trial_wf_pot
