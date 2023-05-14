include 'mkl_rci.f90'

program main
  
  ! 3D model, rate independent, mixed stress-strain control
  
  use subs
  use diag3x3
  use MKL_RCI
  use MKL_RCI_type
  use data_types
  use csv_file

  implicit none

  external resid_coupled
  external resid_yield


  ! variable declaration

  REAL :: a1i, a2i, a1s, a2s, ax, b, bi, bs, c11A, c12A, c13A, c33A, c44A, c11M, c12M, c13M, c33M, c44M, cpA, cpM, delta_i, &
  dlamcr, dlamcr0, dmc, dmc0, e, EA, EM, fll, flp, fmp, Gc_temp, kc, ks, L, lam2use, Lamlam2use, lambda_i, mc, ms, nll, nlp, nmp, &
  nuA, nuM, Q, r1, r2, rs, s, s1, s2, s3, s4, T0, TAf, TAs, Ti, TMf, TMs, tol0, tor, xi_r, xi_r_i

  REAL, DIMENSION(6) :: eps_stop, temp1, temp2, tempDm, v, v1, v2, xi1, xi2, xi_a_i, xi_p_i

  REAL, DIMENSION(12) :: xi0, xi0_y

  REAL, DIMENSION(22) :: x0

  REAL, DIMENSION(32) :: Tdiff

  REAL, DIMENSION(3,3) :: dGcxi_temp, dGceta_temp, epsdiff01, epsdiff02, epsdiff03, epsdiff04, t, t1, t2, t3, t4, test_d_m, &
  testDm, xi_a, xi_p

  REAL, DIMENSION(6,6) :: C, C_a, CA, Ci, CM, Cs, M, Minv, m_, m_inv

  REAL, DIMENSION(:), ALLOCATABLE :: alph_m, d_lam, del_lam, d_m_, d_r, eps13, eps22, eps23, eps33, eps_eq, epse13, epse22, &
  epse23, epse33, epsexp_eq, epsm13, epsm22, epsm23, epsm33, epsp13, epsp22, epsp23, epsp33, Gc, gi, Gs, gs_, lam, Lam_i, Lam_lam, &
  r, sig13, sig22, sig23, sig33, sig_eq, sigd13, sigd22, sigd23, sigd33, sigexp_eq, theta, theta0, yieldcheck, init

  REAL, DIMENSION(:,:), ALLOCATABLE :: data, fjac_coupled, fjac_resid, xi, xi_y

  REAL, DIMENSION(:,:,:),ALLOCATABLE :: a, d_a, D_m, d_eps_p, dgi_, dGceta, dGcxi, dGi, dGs, eps, eps_e, eps_m, eps_p, epsdiff, &
  sig, sig_dev

  INTEGER :: eqval, ii, iter, iter1, iter2, jj, lf, lx, n_inc, n_seg, n_step, RCI_request, res, st_cr, successful, zz, tt

  INTEGER, DIMENSION(6) :: info
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: flag

  LOGICAL :: cmv, conflag, dmsigninv

  CHARACTER(len = 5) :: fn1, fn2

  CHARACTER(len = 13) :: filename

  TYPE(HANDLE_TR) :: handle

  TYPE(resid_yield_args) :: args_yield

  TYPE(resid_coupled_args) :: args_coupled
  
  TYPE(cell), DIMENSION(:), ALLOCATABLE :: fval, x

  
  ! strain loading
  ax = 0.0100
  tor = 0.0000


  ! material parameters


  !! eval-based initiation surface parameters

  a1i = 0.3       ! 0 = tension-compression symmetry, max +/-1
  a2i = 2         ! recommended to be set between 2 and 5, max of 10, min of 2
  bi = 0.06


  !! anisotropy tensor for initiation

  Ci = Ci*0
  Ci(1,1) = 1
  Ci(2,2) = 1
  Ci(3,3) = 1
  Ci(4,4) = 1
  Ci(5,5) = 1
  Ci(6,6) = 1
  

  !! eval-based saturation surface parameters

  a1s = 0.3       ! 0 = tension-compression symmetry, max +/-1
  a2s = 2         ! recommended to be set between 2 and 5, max of 10, min of 2
  bs = 0.03       ! should be less than bi in most cases


  !! anisotropy tensor for saturation

  Cs = Cs*0
  Cs(1,1) = 1
  Cs(2,2) = 1
  Cs(3,3) = 1
  Cs(4,4) = 1
  Cs(5,5) = 1
  Cs(6,6) = 1

  
  !! plasticity parameters

  b = 1           ! isotropic hardening coefficient (s^-1)
  Q = 2           ! isotropic hardening coefficient (MPa*s)
  e = 1.75        ! isotropic hardening parameter (MPa^2)
  
  
  !! kinematic hardening coefficients

  C_a(:,:) = 0
  m_(:,:) = 0
  M(:,:) = 0

  do jj = 1, 6
    C_a(jj,jj) = 150e3
    m_(jj,jj) = 80e3
    M(jj,jj) = 120e3
    m_inv(jj,jj) = 1/80e3
    Minv(jj,jj) = 1/120e3
  end do


  !! material properties

  L = 107         ! latent heat of transformation (MJ/m^3)
  dmc0 = 50       ! martensite reorientation onset stress (MPa)

  TMs = 221.4       ! martensite start temperature (K)
  TMf = 185.6       ! martensite finish temperature (K)
  cpM = 2.79        ! martensite heat capacity (MJ/m^3*K)
  EM = 35e3         ! martensite Young's modulus (GPa)
  nuM = 0.45        ! martensite Poisson's ratio

  c12M = EM*nuM/((1 + nuM)*(1 - 2*nuM))       ! martensite lambda
  c44M = 2*EM/(2*(1+nuM))                     ! martensite 2*mu
  c11M = c12M + c44M                          ! martensite lambda + 2*mu
  c33M = c11M                                 ! martensite lambda + 2*mu
  c13M = c12M                                 ! martensite lambda

  CM = CM*0
  CM(1,1) = c11M
  CM(2,1) = c12M
  CM(3,1) = c13M
  CM(1,2) = c12M
  CM(2,2) = c11M
  CM(3,2) = c13M
  CM(1,3) = c13M
  CM(2,3) = c13M
  CM(3,3) = c33M
  CM(4,4) = c44M
  CM(5,5) = c44M
  CM(6,6) = c44M

  TAs = 266.6       ! austenite start temperature (K)
  TAf = 291.1       ! austenite finish temperature (K)
  cpA = 2.79        ! austenite heat capacity (MJ/m^3*K)
  EA = 55e3         ! austenite Young's modulus (GPa)
  nuA = 0.45        ! austenite Poisson's ratio

  c12A = EA*nuA/((1 + nuA)*(1 - 2*nuA))       ! austenite lambda
  c44A = 2*EA/(2*(1+nuA))                     ! austenite 2*mu
  c11A = c12A + c44A                          ! austenite lambda + 2*mu
  c33A = c11A                                 ! austenite lambda + 2*mu
  c13A = c12A                                 ! austenite lambda

  CA = CA*0
  CA(1,1) = c11A
  CA(2,1) = c12A
  CA(3,1) = c13A
  CA(1,2) = c12A
  CA(2,2) = c11A
  CA(3,2) = c13A
  CA(1,3) = c13A
  CA(2,3) = c13A
  CA(3,3) = c33A
  CA(4,4) = c44A
  CA(5,5) = c44A
  CA(6,6) = c44A


  ! model calibration parameters

  ks = 4          ! saturation slope
  ms = 3.5        ! saturation exponent
  kc = 5e3        ! coupled transformation-plasticity hardening parameter
  mc = 0.5        ! coupled transformation-plasticity hardening exponent
  flp = 1         ! friction (isotropic transformation hardening) coefficient
  nlp = 0.5       ! friction (isotropic transformation hardening) exponent
  fmp = 0         ! friction (isotropic martensite strain hardening) coefficient
  nmp = 0.5       ! friction (isotropic martensite strain hardening) exponent
  fll = 3         ! phenomenological friction fitting coefficient
  nll = 1         ! phenomenological friction fitting exponent


  !! calculated parameters

  T0 = (TAf + TMs)/2                ! equilibrium temperature
  ! dlamcr0 = L/T0*(T0 - TMs)        ! transformation critical driving force
  dlamcr0 = 7.28


  !! temperature profile

  Ti = 295
  Tdiff = Tdiff*0


  ! load profile

  ! epsdiff01 = epsdiff01*0
  ! epsdiff01(3,3) = ax
  ! epsdiff02 = epsdiff02*0
  ! epsdiff02(3,2) = tor
  ! epsdiff02(2,3) = tor
  ! epsdiff03 = epsdiff04*0
  ! epsdiff03(3,2) = -tor
  ! epsdiff03(2,3) = -tor
  ! epsdiff04 = epsdiff04*0
  ! epsdiff04(3,3) = -ax

  ! ALLOCATE(epsdiff(3,3,4))
  ! epsdiff(:,:,1) = epsdiff01
  ! epsdiff(:,:,2) = epsdiff02
  ! epsdiff(:,:,3) = epsdiff03
  ! epsdiff(:,:,4) = epsdiff04

  epsdiff01 = epsdiff01*0
  epsdiff01(3,3) = ax
  epsdiff02 = epsdiff02*0
  epsdiff02(3,3) = -ax

  ALLOCATE(epsdiff(3,3,2))
  epsdiff(:,:,1) = epsdiff01
  epsdiff(:,:,2) = epsdiff02



  ! simulation parameters

  cmv = .TRUE.
  n_seg = 2
  n_step = 50
  n_inc = n_seg*n_step
  tol0 = 1e-6
  iter1 = 1000
  iter2 = 100
  rs = 1
  eps_stop = tol0


  ! build loading path

  
  !! allocate space

  ALLOCATE(eps(3,3,n_inc))
  ALLOCATE(theta0(n_inc))


  !! set initial conditions

  eps(:,:,1) = 0;
  theta0(1) = Ti


  !! apply ramp loading

  do jj = 1, n_seg
    if (jj > 1) then
      eps(:,:,n_step*(jj - 1) + 1) = eps(:,:,n_step*(jj - 1))
      theta0(n_step*(jj - 1) + 1) = theta0(n_step*(jj - 1))
    end if
    do ii = 2, n_step
      eps(:,:,n_step*(jj - 1) + ii) = eps(:,:,n_step*(jj - 1) + 1) + epsdiff(:,:,jj)/(n_step - 1)*(ii - 1)
      theta0(n_step*(jj - 1) + ii) = theta0(n_step*(jj - 1) + 1) + Tdiff(jj)/(n_step - 1)*(ii - 1)
    end do
  end do


  ! allocate space for variables


  !! tensors

  ALLOCATE(sig(3,3,n_inc))            ! stress
  ALLOCATE(sig_dev(3,3,n_inc))        ! deviatoric stress
  ALLOCATE(eps_e(3,3,n_inc))          ! elastic strain
  ALLOCATE(eps_m(3,3,n_inc))          ! martensite strain
  ALLOCATE(eps_p(3,3,n_inc))          ! plastic strain
  ALLOCATE(d_eps_p(3,3,n_inc))        ! plastic strain driving force
  ALLOCATE(a(3,3,n_inc))              ! backstrain
  ALLOCATE(d_a(3,3,n_inc))            ! backstrain driving force
  ALLOCATE(xi_y(12,n_inc))            ! yield trial variables
  ALLOCATE(xi(12,n_inc))              ! plastic trial variables
  ALLOCATE(dgi_(3,3,n_inc))           ! initiation surface derivatives
  ALLOCATE(dGi(3,3,n_inc))            ! initiation surface derivatives
  ALLOCATE(dGs(3,3,n_inc))            ! saturation surface derivatives
  ALLOCATE(dGcxi(3,3,n_inc))          ! coupling surface derivatives
  ALLOCATE(dGceta(3,3,n_inc))         ! coupling surface derivatives
  ALLOCATE(D_m(3,3,n_inc))            ! martensite strain increment direction


  !! scalars

  ALLOCATE(lam(n_inc))                ! volume fraction
  ALLOCATE(del_lam(n_inc))            ! volume fraction increment
  ALLOCATE(theta(n_inc))              ! temperature
  ALLOCATE(d_lam(n_inc))              ! volume fraction driving force
  ALLOCATE(d_m_(n_inc))               ! martensite strain driving force
  ALLOCATE(r(n_inc))                  ! cumulative plastic strain magnitude
  ALLOCATE(d_r(n_inc))                ! cumulative plastic strain driving force
  ALLOCATE(Gi(n_inc))                 ! initiation surface energy
  ALLOCATE(gs_(n_inc))                ! saturation surface energy
  ALLOCATE(Gs(n_inc))                 ! saturation surface energy
  ALLOCATE(Gc(n_inc))                 ! coupling surface energy
  ALLOCATE(Lam_i(n_inc))              ! initiation multiplier
  ALLOCATE(Lam_lam(n_inc))            ! volume fraction multiplier
  ALLOCATE(alph_m(n_inc))             ! martensite strain increment magnitude
  ALLOCATE(yieldcheck(n_inc))         ! plastic yield predictor check
  ALLOCATE(flag(6,n_inc))             ! flags for selecting equations


  !! cells

  ALLOCATE(x(n_inc))
  ALLOCATE(fval(n_inc))


  !! plotting variables

  ALLOCATE(sig13(n_inc))
  ALLOCATE(eps13(n_inc))
  ALLOCATE(epsm13(n_inc))
  ALLOCATE(epsp13(n_inc))
  ALLOCATE(epse13(n_inc))
  ALLOCATE(sigd13(n_inc))
  ALLOCATE(sig22(n_inc))
  ALLOCATE(eps22(n_inc))
  ALLOCATE(epsm22(n_inc))
  ALLOCATE(epsp22(n_inc))
  ALLOCATE(epse22(n_inc))
  ALLOCATE(sigd22(n_inc))
  ALLOCATE(sig23(n_inc))
  ALLOCATE(eps23(n_inc))
  ALLOCATE(epsm23(n_inc))
  ALLOCATE(epsp23(n_inc))
  ALLOCATE(epse23(n_inc))
  ALLOCATE(sigd23(n_inc))
  ALLOCATE(sig33(n_inc))
  ALLOCATE(eps33(n_inc))
  ALLOCATE(epsm33(n_inc))
  ALLOCATE(epsp33(n_inc))
  ALLOCATE(epse33(n_inc))
  ALLOCATE(sigd33(n_inc))
  ALLOCATE(sig_eq(n_inc))
  ALLOCATE(eps_eq(n_inc))
  ALLOCATE(sigexp_eq(n_inc))
  ALLOCATE(epsexp_eq(n_inc)) 


  ! initial values


  !! tensors

  sig(:,:,1) = 0
  sig_dev(:,:,1) = 0
  eps_e(:,:,1) = 0
  eps_m(:,:,1) = 0
  eps_p(:,:,1) = 0
  d_eps_p(:,:,1) = 0
  a(:,:,1) = 0
  d_a(:,:,1) = 0
  xi_y(:,1) = 0
  xi(:,1) = 0
  dgi_(:,:,1) = 0
  dGi(:,:,1) = 0
  dGs(:,:,1) = 0
  dGcxi(:,:,1) = 0
  dGceta(:,:,1) = 0
  D_m(:,:,1) = 0
  


  !! scalars
  
  del_lam(1) = 0
  d_lam(1) = 0
  d_m_(1) = 0
  r(1) = 0
  d_r(1) = 0
  Gi(1) = -1
  gs_(1) = -1
  Gs(1) = -1
  Gc(1) = -1
  Lam_i(1) = 0
  Lam_lam(1) = 0
  alph_m(1) = 0
  yieldcheck(1) = 0
  flag(:,1) = 0
  xi0_y(:) = 0
  xi0(:) = 0
  x0(:) = 0
  
  dmsigninv = .FALSE.
  

  !! plotting variables

  sig13(1) = 0
  eps13(1) = 0
  epsm13(1) = 0
  epsp13(1) = 0
  epse13(1) = 0
  sigd13(1) = 0
  sig22(1) = 0
  eps22(1) = 0
  epsm22(1) = 0
  epsp22(1) = 0
  epse22(1) = 0
  sigd22(1) = 0
  sig23(1) = 0
  eps23(1) = 0
  epsm23(1) = 0
  epsp23(1) = 0
  epse23(1) = 0
  sigd23(1) = 0
  sig33(1) = 0
  eps33(1) = 0
  epsm33(1) = 0
  epsp33(1) = 0
  epse33(1) = 0
  sigd33(1) = 0
  sig_eq(1) = 0
  eps_eq(1) = 0
  sigexp_eq(1) = 0
  epsexp_eq(1) = 0

  
  !! define non-zero initial conditions

  theta(1) = theta0(1)


  !! initial volume fraction

  if (theta0(1) < TMs .and. theta0(1) > TMf) then
    lam(1) = 1 - (theta0(1) - TMf)/(TMs - TMf)
  else if ( theta0(1) < TMf ) then
    lam(1) = 1
  end if


  !! volume fraction driving force at initial conditions

  d_lam(1) = (cpM - cpA)*theta(1)*LOG(theta(1)/T0) - L*(theta(1) - T0)/T0 - fll*lam(1)**nll


  !! set initial transformation Lagrange multiplier if non-zero

  if ((d_lam(1) - dlamcr0 > tol0 .and. 1 - lam(1) < tol0) .or. (d_lam(1) + dlamcr0 < -tol0 .and. lam(1) < tol0)) then
    Lam_lam(1) = (cpM - cpA)*theta(1)*LOG(theta(1)/T0)  - L*(theta(1) - T0)/T0 - fll*lam(1)**nll - dlamcr0
  end if


  ! run the simulation

  do ii = 2,6
    
    !! calculate rule of mixture quantities

    C = lam(ii - 1)*CM + (1 - lam(ii - 1))*CA


    !! isothermal

    theta(ii) = theta0(ii)


    !! calculate strains along non-loading axes to enforce sig_11 = sig_22 = sig_12 = sig_13 = 0
    
    CALL axialtorstrain(eps(3,2,ii), eps(3,3,ii), C, lam(ii - 1), eps_m(:,:,ii - 1), eps_p(:,:,ii - 1), cmv, eps(:,:,ii))


    !! predicted elastic strain

    eps_e(:,:,ii) = eps(:,:,ii) - lam(ii - 1)*eps_m(:,:,ii - 1) - eps_p(:,:,ii - 1)


    !! predicted stress

    CALL mat2voigt(eps_e(:,:,ii), v1)
    CALL voigt2mat(MATMUL(C, v1), sig(:,:,ii))
    CALL dev(sig(:,:,ii), sig_dev(:,:,ii))
    dmc = dmc0*SQRT(2./3.) + fmp*r(ii - 1)**nmp


    !! martensite strain driving force prediction

    test_d_m = sig_dev(:,:,ii) - dGs(:,:,ii - 1) - dGcxi(:,:,ii - 1)
    
    CALL mag(test_d_m, d_m_(ii))

    CALL non0mag(test_d_m - Lam_i(ii - 1)*dgi_(:,:,ii - 1), s)

    testDm = (test_d_m - Lam_i(ii - 1)*dgi_(:,:,ii - 1))/s

    if (d_m_(ii) - dmc > tol0) then
      CALL mag(testDm - D_m(:,:,ii - 1), s)
      if (s < tol0 .or. s == 2) then
        flag(3,ii) = 1
      else
        flag(1:3,ii) = 1
      end if
    end if

    CALL scalar22(dmc*testDm, test_d_m, s1)
    CALL scalar22(testDm, D_m(:,:,ii - 1), s2)

    if ((s1 < 0 .or. s2 < 0) .and. .not. dmsigninv) then
      dmsigninv = .TRUE.
    end if


    !! reset eps_m if necessary

    CALL mag(eps_m(:,:,ii - 1), s)

    if (theta(ii) > TAf .and. lam(ii - 1) < tol0 .and. (dmsigninv .or. d_m_(ii) - dmc < tol0) .and. s > 0) then
      eps_m(:,:,ii - 1) = 0
      alph_m(ii - 1) = 0
      D_m(:,:,ii - 1) = 0
      Lam_i(ii - 1) = 0
      CALL non0mag(test_d_m, s)
      testDm = test_d_m/s
      dmsigninv = .FALSE.
      if (d_m_(ii) - dmc > tol0) then
        flag(1:3,ii) = 1
      end if   
    end if


    !! volume fraction driving force prediction

    CALL scalar22(sig(:,:,ii), eps_m(:,:,ii - 1), s1)
    CALL mat2voigt(eps_e(:,:,ii), v1)
    CALL voigt2mat(MATMUL(CM - CA, v1), t1)
    CALL scalar22(eps_m(:,:,ii - 1), t1, s2)
    CALL scalar22(eps_m(:,:,ii - 1), dGs(:,:,ii - 1), s3)
    CALL scalar22(eps_m(:,:,ii - 1), dGcxi(:,:,ii - 1), s4)
    d_lam(ii) = s1 - s2/2 + (cpM - cpA)*theta(ii)*LOG(theta(ii)/T0) - L*(theta(ii) - T0)/T0 - s3 - fll*lam(ii - 1)**nll - s4

    dlamcr = dlamcr0 + flp*r(ii - 1)**nlp

    CALL mag(eps_m(:,:,ii - 1), s)

    if (d_lam(ii) - dlamcr > tol0 .and. lam(ii - 1) < 1 .and. s > tol0) then
      flag(4,ii) = 1
      if (1 - lam(ii - 1) < tol0) then
        flag(5,ii) = 1
      end if
    else if (d_lam(ii) + dlamcr < - tol0 .and. lam(ii - 1) > 0 .and. s > tol0) then
      dlamcr = -dlamcr
      flag(4,ii) = 1
      if (lam(ii - 1) < tol0) then
        flag(5,ii) = 1
      end if
    end if

    if ((d_lam(ii) - ABS(dlamcr) > tol0 .and. 1 - lam(ii - 1) < tol0) .or. (d_lam(ii) + ABS(dlamcr) < -tol0 .and. lam(ii - 1) < &
    tol0)) then
      flag(5,ii) = 1
    end if


    !! eps_p driving force prediction

    d_eps_p(:,:,ii) = sig_dev(:,:,ii) - dGceta(:,:,ii - 1)


    !! a driving force prediction

    CALL mat2voigt(a(:,:,ii - 1), v1)
    CALL voigt2mat(MATMUL(C_a, v1), t1)

    d_a(:,:,ii) = -2./3.*t1


    !! r driving force prediction

    d_r(ii) = -b*Q*r(ii - 1)


    !! calculate xi_p, xi_a, xi_r for plastic yield check

    CALL mag(eps_p(:,:,ii - 1), s)

    if (s < tol0) then
      CALL mat2voigt(sig_dev(:,:,ii), v1)
      CALL voigt2mat(MATMUL(m_inv + Minv, v1), t1)
      CALL scalar22(t1, sig_dev(:,:,ii), s1)
      lambda_i = 1/SQRT(8.)*SQRT(s1)
      xi_p = 0.25/lambda_i*t1
      CALL mat2voigt(xi_p, xi1)
      xi_a = 0.25/lambda_i*t1
      CALL mat2voigt(xi_a, xi2)
      xi_y(:,ii) = (/xi1(2:6), xi2(2:6), xi_r, lambda_i/)
    else
      xi0_y = xi_y(:,ii - 1)


      !! evaluate elastic predictor function

      
      !!! arguments for resid_yield function

      args_yield%d_eps_p = d_eps_p(:,:,ii)
      args_yield%d_a = d_a(:,:,ii)
      args_yield%d_r = d_r(ii)
      args_yield%M = M
      args_yield%m_ = m_
      args_yield%b = b
      args_yield%e = e

      lf = 12
      lx = 12

      ALLOCATE(fjac_resid(lf,lx))
      fjac_resid(:,:) = 0.


      !!! call solver

      if (strnlsp_init(handle, lx, lf, xi0_y, eps_stop, iter1, iter2, rs) /= TR_SUCCESS) then
        stop '| error in strnlsp_init 631'
      end if

      if (strnlsp_check(handle, lx, lf, fjac_resid, xi_y(:,ii), eps_stop, info) /= TR_SUCCESS) then
        stop '| error in strnlsp_check 635'
      end if  

      RCI_request = 0
      successful = 0

      do while (successful == 0)

        if (strnlsp_solve(handle, xi_y(:,ii), fjac_resid, RCI_request) /= TR_SUCCESS) then  
          stop '| error in strnlsp_solve 644'
        end if

        if (RCI_request == -1 .or. RCI_request == -2 .or. RCI_request == -3 .or. RCI_request == -4 .or. RCI_request == -5 .or. &
        RCI_request == -6) then
          successful = 1
        end if

        if (RCI_request == 1) then
          CALL resid_yield(lf, lx, xi0_y, xi_y(:,ii), args_yield)
        end if

        if (RCI_request == 2) then
          if (sjacobix(resid_yield, lf, lx, fjac_resid, xi0_y, tol0, %VAL(LOC(args_yield))) /= TR_SUCCESS) then
            stop '| error in sjacobix 658'
          end if
        end if
      end do

      if (strnlsp_get(handle, iter, st_cr, r1, r2) /= TR_SUCCESS) then
        stop '| error in strnlsp_get 664'
      end if

      if (strnlsp_delete(handle) /= TR_SUCCESS) then
        stop '| error in strnlsp_delete 668'
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!! if lambda is negative, it is non-physical, need to reset the inital guesses and solve again

      if (xi_y(12,ii) < -tol0) then

        CALL mat2voigt(sig_dev(:,:,ii), v)
        CALL voigt2mat(MATMUL(m_inv + Minv, v), t1)
        CALL voigt2mat(MATMUL(m_inv - Minv, v), t2)
        CALL scalar22(t1, sig_dev(:,:,ii), s)

        lambda_i = 1/SQRT(8.)*SQRT(s)
        CALL mat2voigt(0.25/lambda_i*t1, xi_p_i)
        CALL mat2voigt(0.25/lambda_i*t2, xi_a_i)
        xi_r_i = b
      
        xi0_y = (/xi_p_i(2:6), xi_a_i(2:6), xi_r_i, lambda_i/)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        fjac_resid(:,:) = 0
        
        if (strnlsp_init(handle, lx, lf, xi0_y, eps_stop, iter1, iter2, rs) /= TR_SUCCESS) then
          stop '| error in strnlsp_init 697'
        end if

        if (strnlsp_check(handle, lx, lf, fjac_resid, xi_y(:,ii), eps_stop, info) /= TR_SUCCESS) then
          stop '| error in strnlsp_check 701'
        end if  

        RCI_request = 0
        successful = 0

        do while (successful == 0)

          if (strnlsp_solve(handle, xi_y(:,ii), fjac_resid, RCI_request) /= TR_SUCCESS) then  
            stop '| error in strnlsp_solve 710'
          end if

          if (RCI_request == -1 .or. RCI_request == -2 .or. RCI_request == -3 .or. RCI_request == -4 .or. RCI_request == -5 .or. &
          RCI_request == -6) then
            successful = 1
          end if

          if (RCI_request == 1) then
            CALL resid_yield(lf, lx, xi0_y, xi_y(:,ii), args_yield)
          end if

          if (RCI_request == 2) then
            if (sjacobix(resid_yield, lf, lx, fjac_resid, xi0_y, tol0, %VAL(LOC(args_yield))) /= TR_SUCCESS) then
              stop '| error in sjacobix 724'
            end if
          end if
        end do

        if (strnlsp_get(handle, iter, st_cr, r1, r2) /= TR_SUCCESS) then
          stop '| error in strnlsp_get 730'
        end if

        if (strnlsp_delete(handle) /= TR_SUCCESS) then
          stop '| error in strnlsp_delete 734'
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      end if
      
      CALL voigt2mat((/-(xi_y(1,ii) + xi_y(2,ii)), xi_y(1:5,ii)/), xi_p)
      CALL voigt2mat((/-(xi_y(6,ii) + xi_y(7,ii)), xi_y(6:10,ii)/), xi_a)
      xi_r = xi_y(11,ii)

    end if


    !! check for yielding

    CALL scalar22(xi_p, d_eps_p(:,:,ii), s1)
    CALL scalar22(xi_a, d_a(:,:,ii), s2)
    yieldcheck(ii) = s1 + s2 + xi_r*d_r(ii)
    eqval = 0

    if (yieldcheck(ii) > 1 .or. yieldcheck(ii) < -tol0) then

      flag(6,ii) = 1

      if (yieldcheck(ii) > 1) then
        eqval = 1
      end if

      xi0 = xi(:,ii - 1)

      if (norm2(xi(:,ii - 1)) < tol0) then
        CALL mat2voigt(xi_p, v1)
        CALL voigt2mat(MATMUL(C, v1), t1)
        CALL dev(t1, t2)
        CALL scalar22(t2, xi_p, s1)
        CALL mat2voigt(xi_a, v2)
        CALL voigt2mat(MATMUL(C_a, v2), t3)
        CALL dev(t3, t4)
        CALL scalar22(t4, xi_a, s2)
        delta_i = (yieldcheck(ii) - eqval)/s1 + 2./3.*s2 + b*Q*xi_r**2
        CALL mat2voigt(xi_p, temp1)
        CALL mat2voigt(xi_a, temp2)
        xi0 = (/temp1(2:6), temp2(2:6), xi_r, delta_i/)
      else if (xi0(12) < tol0) then
        xi0 = xi_y(:,ii)
         CALL mat2voigt(xi_p, v1)
        CALL voigt2mat(MATMUL(C, v1), t1)
        CALL dev(t1, t2)
        CALL scalar22(t2, xi_p, s1)
        CALL mat2voigt(xi_a, v2)
        CALL voigt2mat(MATMUL(C_a, v2), t3)
        CALL dev(t3, t4)
        CALL scalar22(t4, xi_a, s2)
        xi0(12) = (yieldcheck(ii) - eqval)/s1 + 2./3.*s2  + b*Q*xi_r**2
      end if

      CALL mag(eps_p(:,:,ii - 1), s)

      if (flag(4,ii) == 0 .and. lam(ii - 1) < tol0 .and. s < tol0) then
        CALL surf_coupled(lam(ii - 1)*eps_m(:,:,ii - 1), xi0(12)*xi_p, kc, mc, Gc_temp, dGcxi_temp, dGceta_temp)
        CALL non0mag(dGcxi_temp, s1)
        testDm = -dGcxi_temp/s1
        CALL non0mag(D_m(:,:,ii - 1) - testDm, s2)
        testDm = (D_m(:,:,ii - 1) - testDm)/s2
      end if

    end if

    !! correct for inelastic evolution

    if (NORM2(REAL(flag(:,ii))) > 0) then

      conflag = .TRUE.
      lam2use = lam(ii - 1)
      Lamlam2use = Lam_lam(ii - 1)

      do while (conflag .and. NORM2(REAL(flag(:,ii))) > 0)
        
        !!! build initial guess vector
        
        zz = 1

        if (flag(1,ii) == 1) then
          x0(zz) = alph_m(ii - 1)
          if (ABS(x0(zz)) < tol0 .or. ABS(gi(ii - 1)) < tol0) then
            x0(zz) = 1e-4
          end if
          zz = zz + 1
        end if

        if (flag(2,ii) == 1) then
          CALL mat2voigt(testDm, tempDm)
          x0(zz:zz + 4) = tempDm(2:6)
          zz = zz + 5
        end if

        if (flag(3,ii) == 1) then
          x0(zz) = Lam_i(ii - 1)
          if (ABS(x0(zz)) < tol0) then
            x0(zz) = 1
          end if
          zz = zz + 1
        end if
        
        if (flag(4,ii) == 1) then
          if (lam(ii - 1) + del_lam(ii - 1) < 0) then
            x0(zz) = (tol0 - lam(ii - 1))/2
          else if (lam(ii - 1) + del_lam(ii - 1) > 1) then
            x0(zz) = (1 - lam(ii - 1))/2
          else if (del_lam(ii - 1) == 0) then
            x0(zz) = SIGN(.0001, dlamcr)
          else
            x0(zz) = del_lam(ii - 1)
          end if
          zz = zz + 1
        end if
      
        if (flag(5,ii) == 1) then
          x0(zz) = Lamlam2use
          zz = zz + 1
        end if

        if (flag(6,ii) == 1) then
          x0(zz:zz+11) = xi0
          zz = zz + 12

        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!! arguments for resid_coupled function
        args_coupled%eps = eps(:,:,ii)
        args_coupled%theta = theta(ii)
        args_coupled%eps_p_pr = eps_p(:,:,ii - 1)
        args_coupled%r_pr = r(ii - 1)
        args_coupled%a_pr = a(:,:,ii - 1)
        args_coupled%eps_m_pr = eps_m(:,:,ii - 1)
        args_coupled%lam_pr = lam2use
        args_coupled%bi = bi
        args_coupled%a1i = a1i
        args_coupled%a2i = a2i
        args_coupled%Ci = Ci
        args_coupled%bs = bs
        args_coupled%ks = ks
        args_coupled%ms = ms
        args_coupled%a1s = a1s
        args_coupled%a2s = a2s
        args_coupled%Cs = Cs
        args_coupled%kc = kc
        args_coupled%mc = mc
        args_coupled%CA = CA
        args_coupled%CM = CM
        args_coupled%cpA = cpA
        args_coupled%cpM = cpM
        args_coupled%T0 = T0
        args_coupled%L = L
        args_coupled%dmc0 = dmc0
        args_coupled%dlamcr0 = dlamcr0
        args_coupled%dlamcr = dlamcr
        args_coupled%flag = flag(:,ii)
        args_coupled%cmv = cmv
        args_coupled%D_m = D_m(:,:,ii - 1)
        args_coupled%flp = flp
        args_coupled%fmp = fmp
        args_coupled%fll = fll
        args_coupled%nlp = nlp
        args_coupled%nmp = nmp
        args_coupled%nll = nll
        args_coupled%C_a = C_a
        args_coupled%M = M
        args_coupled%m_ = m_
        args_coupled%b = b
        args_coupled%Q = Q
        args_coupled%e = e
        args_coupled%eqval = eqval
        args_coupled%xi0 = xi0

        
        ! print*,ii
        ! print*,''
        ! print*,args_coupled%eps
        ! print*,''
        ! print*,args_coupled%theta
        ! print*,''
        ! print*,args_coupled%eps_p_pr
        ! print*,''
        ! print*,args_coupled%r_pr
        ! print*,''
        ! print*,args_coupled%a_pr
        ! print*,''
        ! print*,args_coupled%eps_m_pr
        ! print*,''
        ! print*,args_coupled%lam_pr
        ! print*,''
        ! print*,args_coupled%bi
        ! print*,''
        ! print*,args_coupled%a1i
        ! print*,''
        ! print*,args_coupled%a2i
        ! print*,''
        ! print*,args_coupled%Ci
        ! print*,''
        ! print*,args_coupled%bs
        ! print*,''
        ! print*,args_coupled%ks
        ! print*,''
        ! print*,args_coupled%ms
        ! print*,''
        ! print*, args_coupled%a1s
        ! print*,''
        ! print*,args_coupled%a2s
        ! print*,''
        ! print*,args_coupled%Cs
        ! print*,''
        ! print*,args_coupled%kc
        ! print*,''
        ! print*, args_coupled%mc
        ! print*,''
        ! print*,args_coupled%CA
        ! print*,''
        ! print*,args_coupled%CM
        ! print*,''
        ! print*,args_coupled%cpA 
        ! print*,''
        ! print*,args_coupled%cpM
        ! print*,''
        ! print*,args_coupled%T0
        ! print*,''
        ! print*,args_coupled%L
        ! print*,''
        ! print*,args_coupled%dmc0
        ! print*,''
        ! print*,args_coupled%dlamcr0
        ! print*,''
        ! print*,args_coupled%dlamcr
        ! print*,''
        ! print*,args_coupled%flag
        ! print*,''
        ! print*,args_coupled%cmv
        ! print*,''
        ! print*,args_coupled%D_m
        ! print*,''
        ! print*,args_coupled%flp
        ! print*,''
        ! print*,args_coupled%fmp
        ! print*,''
        ! print*,args_coupled%fll
        ! print*,''
        ! print*,args_coupled%nlp
        ! print*,''
        ! print*,args_coupled%nmp
        ! print*,''
        ! print*,args_coupled%nll
        ! print*,''
        ! print*,args_coupled%C_a
        ! print*,''
        ! print*,args_coupled%M
        ! print*,''
        ! print*,args_coupled%m_
        ! print*,''
        ! print*,args_coupled%b
        ! print*,''
        ! print*,args_coupled%Q
        ! print*,''
        ! print*,args_coupled%e
        ! print*,''
        ! print*,args_coupled%eqval
        ! print*,''
        ! print*,args_coupled%xi0
        ! print*,''
        ! print*,'-----'

        
        lx = zz - 1
        lf = lx

        if (ALLOCATED(fjac_coupled)) then
          DEALLOCATE(fjac_coupled)
        end if
        ALLOCATE(fjac_coupled(lf,lx))
        fjac_coupled(:,:) = 0
        
        if (ALLOCATED(x(ii)%component)) then
          DEALLOCATE(x(ii)%component)
        end if
        ALLOCATE(x(ii)%component(lx))

        if (ALLOCATED(fval(ii)%component)) then
          DEALLOCATE(fval(ii)%component)
        end if
        ALLOCATE(fval(ii)%component(lf))

        x(ii)%component = x0(1:zz - 1)
        fval(ii)%component = 0

        if (strnlsp_init(handle, lx, lf, x(ii)%component, eps_stop, iter1, iter2, rs) /= TR_SUCCESS) then
          stop '| error in strnlsp_init 929'
        end if
        
        if (strnlsp_check(handle, lx, lf, fjac_coupled, fval(ii)%component, eps_stop, info) /= TR_SUCCESS) then
          stop '| error in strnlsp_check 933'
        end if  
        
        RCI_request = 0
        successful = 0

        do while (successful == 0)

          if (strnlsp_solve(handle, fval(ii)%component, fjac_coupled, RCI_request) /= TR_SUCCESS) then  
            stop '| error in strnlsp_solve 942'
          end if
          
          if (RCI_request == -1 .or. RCI_request == -2 .or. RCI_request == -3 .or. RCI_request == -4 .or. RCI_request == -5 .or. &
          RCI_request == -6) then
            successful = 1
          end if

          if (RCI_request == 1) then
            CALL resid_coupled(lx, lf, x(ii)%component, fval(ii)%component, args_coupled)
          end if

          if (RCI_request == 2) then
            if (sjacobix(resid_coupled, lx, lf, fjac_coupled, x(ii)%component, tol0, %VAL(LOC(args_coupled))) /= TR_SUCCESS) then
              stop '| error in sjacobix 956'
            end if
          end if
        end do

        if (strnlsp_get(handle, iter, st_cr, r1, r2) /= TR_SUCCESS) then
          stop '| error in strnlsp_get 962'
        end if
        
        if (strnlsp_delete(handle) /= TR_SUCCESS) then
          stop '| error in strnlsp_delete 966'
        end if

        print*,ii
        print*,''
        print*,'x'
        print*,x(ii)%component
        print*,''
        print*,'fval'
        print*,fval(ii)%component
        print*,''
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !! update inelastic variables

        zz = 1

        if (flag(1,ii) == 1) then
          alph_m(ii) = x(ii)%component(zz)
          zz = zz + 1
        end if

        CALL mag(testDm - D_m(:,:,ii - 1),s)
        if (flag(2,ii) == 1) then
          CALL voigt2mat((/-(x(ii)%component(zz) + x(ii)%component(zz + 1)), x(ii)%component(zz:zz + 4)/), D_m(:,:,ii))
          zz = zz + 5
        else if (s == 2) then
          D_m(:,:,ii) = testDm
        else
          D_m(:,:,ii) = D_m(:,:,ii - 1)
        end if
        

        if (flag(3,ii) == 1) then
          Lam_i(ii) = x(ii)%component(zz)
          zz = zz + 1
        end if

        if (flag(4,ii) == 1) then
          del_lam(ii) = x(ii)%component(zz)
          zz = zz + 1
        end if

        if (flag(5,ii) == 1) then
          Lam_lam(ii) = x(ii)%component(zz)
          zz = zz + 1
        end if

        if (flag(6,ii) == 1) then
          xi(:,ii) = (/xi(1:11,ii - 1), 0./)
          CALL voigt2mat((/-(xi(1,ii) + xi(2,ii)), xi(1:5,ii)/), t)
          eps_p(:,:,ii) = eps_p(:,:,ii - 1) + xi(12,ii)*t
          CALL voigt2mat((/-(xi(6,ii) + xi(7,ii)), xi(6:10,ii)/), t)
          a(:,:,ii) = a(:,:,ii - 1) + xi(12,ii)*t
          r(ii) = r(ii - 1) + xi(12,ii)*xi(11,ii)
        else
          xi(:,ii) = (/xi(1:11,ii), 0./)
          eps_p(:,:,ii) = eps_p(:,:,ii - 1)
          a(:,:,ii) = a(:,:,ii - 1)
          r(ii) = r(ii - 1)
        end if

        !! update variables and surfaces
        
        lam(ii) = lam(ii - 1) + del_lam(ii)

        eps_m(:,:,ii) = eps_m(:,:,ii - 1) + alph_m(ii)*D_m(:,:,ii)

        C = lam(ii)*CM + (1 - lam(ii))*CA
        CALL axialtorstrain(eps(3,2,ii), eps(3,3,ii), C, lam(ii), eps_m(:,:,ii), eps_p(:,:,ii), cmv, eps(:,:,ii))
        eps_e(:,:,ii) = eps(:,:,ii) - lam(ii)*eps_m(:,:,ii) - eps_p(:,:,ii)
        CALL mat2voigt(eps_e(:,:,ii), v)
        CALL voigt2mat(MATMUL(C, v), sig(:,:,ii))
        CALL dev(sig(:,:,ii), sig_dev(:,:,ii))

        CALL mat2voigt(eps_m(:,:,ii), v)
        CALL voigt2mat(MATMUL(Ci, v), t)
        CALL surf_init(t, a1i, a2i, bi, gi(ii), dgi_(:,:,ii))

        CALL voigt2mat(MATMUL(Cs, v), t)
        CALL surf_sat(lam(ii)*t, a1s, a2s, bs, ks, ms, gs_(ii), Gs(ii), dGs(:,:,ii))

        CALL surf_coupled(lam(ii)*eps_m(:,:,ii), eps_p(:,:,ii), kc, mc, Gc(ii), dGcxi(:,:,ii), dGceta(:,:,ii))

        CALL mag(sig_dev(:,:,ii) - Lam_i(ii)*dgi_(:,:,ii) - dGs(:,:,ii) - dGcxi(:,:,ii), d_m_(ii))

        
        !! check to ensure hard constraints were not violated
        
        conflag = .FALSE.

        if (alph_m(ii) < - tol0 .and. flag(1,ii) == 1) then
          alph_m(ii) = 0
          flag(1:2,ii) = 0
          conflag = .TRUE.
        end if

        if (d_m_(ii) - dmc > tol0 .and. flag(3,ii) == 0) then
          flag(1:3,ii) = 1
          conflag = .TRUE.
        end if

        if (Lam_i(ii) < -tol0) then
          Lam_i(ii) = 0
          flag(3,ii) = 1
          conflag = .TRUE.
        end if

        if (lam(ii) > 1) then
          del_lam(ii) = 1 - lam(ii - 1)
          lam(ii) = 1
          lam2use = lam(ii)
          Lamlam2use = 1
          flag(4,ii) = 0
          flag(5,ii) = 1
          conflag = .TRUE.
        end if

        if (lam(ii) < 0) then
          del_lam(ii) = -lam(ii - 1)
          lam(ii) = 0
          lam2use = lam(ii)
          Lamlam2use = 0
          flag(4,ii) = 1
          flag(5,ii) = 0
          conflag = .TRUE.
        end if

        if (Lam_lam(ii)*del_lam(ii) < 0) then
          Lam_lam(ii) = 0
          Lamlam2use = 0
          flag(4,ii) = 1
          flag(5,ii) = 0
          conflag = .TRUE.
        end if

      end do
      
      
      !! update volume fraction driving force

      CALL scalar22(sig(:,:,ii), eps_m(:,:,ii - 1), s1)
      CALL mat2voigt(eps_e(:,:,ii), v)
      CALL voigt2mat(MATMUL(CM - CA, v), t)
      CALL scalar22(eps_e(:,:,ii), t, s2)
      CALL scalar22(eps_m(:,:,ii), dGs(:,:,ii), s3)
      CALL scalar22(eps_m(:,:,ii), dGcxi(:,:,ii), s4)

      d_lam(ii) = s1 - s2/2 + (cpM - cpA)*theta(ii)*LOG(theta(ii)/T0) - L*(theta(ii) - T0)/T0 - s3 - fll*lam(ii)**nll - s4
      

      !! update plastic driving forces

      d_eps_p(:,:,ii) = sig_dev(:,:,ii) - dGceta(:,:,ii)

      CALL mat2voigt(a(:,:,ii), v)
      CALL voigt2mat(MATMUL(C_a, v), t)
      d_a(:,:,ii) = -2./3.*t

      d_r(ii) = -b*Q*r(ii)

    else

      D_m(:,:,ii) = D_m(:,:,ii - 1)
      eps_m(:,:,ii) = eps_m(:,:,ii - 1)
      Lam_i(ii) = Lam_i(ii - 1)
      lam(ii) = lam(ii - 1)
      Lam_lam(ii) = Lam_lam(ii - 1)
      eps_p(:,:,ii) = eps_p(:,:,ii - 1)
      a(:,:,ii) = a(:,:,ii - 1)
      r(ii) = r(ii - 1)
      xi(:,ii) = (/xi(1:11,ii - 1), 0./)

      CALL mat2voigt(eps_m(:,:,ii), v)
      CALL voigt2mat(MATMUL(Ci, v), t)
      CALL surf_init(t, a1i, a2i, bi, gi(ii), dgi_(:,:,ii))

      CALL voigt2mat(MATMUL(Cs, v), t) 
      CALL surf_sat(lam(ii)*t, a1s, a2s, bs, ks, ms, gs_(ii), Gs(ii), dGs(:,:,ii))

      CALL surf_coupled(lam(ii)*eps_m(:,:,ii), eps_p(:,:,ii), kc, mc, Gc(ii), dGcxi(:,:,ii), dGceta(:,:,ii))

      CALL mag(sig_dev(:,:,ii) - Lam_i(ii)*dgi(:,:,ii) - dGs(:,:,ii) - dGcxi(:,:,ii), d_m_(ii))

    end if
    

  !! update plotting varaibles

  eps13(ii) = eps(3,1,ii)
  epsm13(ii) = eps_m(3,1,ii)
  epsp13(ii) = eps_p(3,1,ii)
  epse13(ii) = eps_e(3,1,ii)
  sig13(ii) = sig(3,1,ii)
  sigd13(ii) = sig_dev(3,1,ii)
  
  eps22(ii) = eps(2,2,ii)
  epsm22(ii) = eps_m(2,2,ii)
  epsp22(ii) = eps_p(2,2,ii)
  epse22(ii) = eps_e(2,2,ii)
  sig22(ii) = sig(2,2,ii)
  sigd22(ii) = sig_dev(2,2,ii)

  eps23(ii) = eps(3,2,ii)
  epsm23(ii) = eps_m(3,2,ii)
  epsp23(ii) = eps_p(3,2,ii)
  epse23(ii) = eps_e(3,2,ii)
  sig23(ii) = sig(3,2,ii)
  sigd23(ii) = sig_dev(3,2,ii)

  eps33(ii) = eps(3,3,ii)
  epsm33(ii) = eps_m(3,3,ii)
  epsp33(ii) = eps_p(3,3,ii)
  epse33(ii) = eps_e(3,3,ii)
  sig33(ii) = sig(3,3,ii)
  sigd33(ii) = sig_dev(3,3,ii)

  CALL scalar22(sig_dev(:,:,ii), sig_dev(:,:,ii), s)
  sig_eq(ii) = SQRT(3./2.*s)
  CALL dev(eps(:,:,ii), t)
  CALL scalar22(t, t, s)
  eps_eq(ii) = SQRT(2./3.*s)
  sigexp_eq(ii) = SQRT(sig33(ii)**2 + 3*sig23(ii)**2)
  epsexp_eq(ii) = SQRT(eps33(ii)**2 + 4./3.*eps23(ii)**2)

  end do
  

  !! export results to csv file

  write(fn1,'(f5.4)') ax
  write(fn2, '(f5.4)') tor
  filename = fn1(2:5)//'-'//fn2(2:5)//'.csv'

  ALLOCATE(data(8,n_inc))
  data(1,:) = eps_eq
  data(2,:) = lam
  data(3,:) = eps33
  data(4,:) = sig33
  data(5,:) = eps23
  data(6,:) = sig23
  data(7,:) = ABS(d_lam)
  data(8,:) = d_m_

  OPEN(unit = 1,file = filename,status = 'unknown')
  CALL csv_write(1,data)

end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine resid_coupled(lx, lf, x, f, args_coupled)

  use subs
  use data_types

  implicit none
  
  INTEGER :: lf, lx
  REAL :: x(lx), f(lf)
  TYPE(resid_coupled_args) :: args_coupled 
  
  REAL :: alph_m, C(6,6), d_a(3,3), d_eps_p(3,3), d_m_c(3,3), del_lam, delta, dGceta(3,3), dGcxi(3,3), dgi(3,3), dgs(3,3), &
  dh_dxi_a(3,3), dh_dxi_p(3,3), dh_dxi_r, dlamcr_c, dlamsign, dmc, eps(3,3), eps_e(3,3), Gc, gi, Gs, gs_, Lam_i, Lam_lam,lcv, &
  s, s1, s2, s3, s4, sig(3,3), sigdev(3,3), t(3,3), t1(3,3), t2(3,3), v(6), v1(6), v2(6), xi_a(3,3), xi_dot_grad_h, xi_p(3,3), &
  xi_r
  INTEGER :: zz

  D_m_c = args_coupled%D_m
  dlamcr_c = args_coupled%dlamcr


  ! reduce unknowns as necessary

  zz = 1

  if (args_coupled%flag(1) == 0) then
    alph_m = 0
  else
    alph_m = x(zz)
    zz = zz + 1
  end if

  if (args_coupled%flag(2) == 1) then
    CALL voigt2mat((/-(x(zz) + x(zz + 1)), x(zz:zz + 4)/), D_m_c)
    zz = zz + 5
  end if

  if (args_coupled%flag(3) == 0) then
    Lam_i = 0
  else
    Lam_i = x(zz)
    zz = zz + 1
  end if

  if (args_coupled%flag(4) == 0) then
    del_lam = 0
  else
    del_lam = x(zz)
    zz = zz + 1
  end if

  if (args_coupled%flag(5) == 0) then
    Lam_lam = 0
  else
    Lam_lam = x(zz)
    zz = zz + 1
  end if

  if (args_coupled%flag(6) == 0) then
    CALL voigt2mat((/-(args_coupled%xi0(1) + args_coupled%xi0(2)), args_coupled%xi0(1:5)/), xi_p)
    CALL voigt2mat((/-(args_coupled%xi0(6) + args_coupled%xi0(7)), args_coupled%xi0(6:10)/), xi_a)
    xi_r = args_coupled%xi0(11)
    delta = 0
  else
    CALL voigt2mat((/-(x(zz) + x(zz + 1)), x(zz:zz + 4)/), xi_p)
    CALL voigt2mat((/-(x(zz + 5) + x(zz + 6)), x(zz + 5:zz + 9)/), xi_a)
    xi_r = x(zz + 10)
    delta = x(zz + 11)
  end if


  ! make preliminary calculations

  C = (args_coupled%lam_pr + del_lam)*args_coupled%CM + (1 - (args_coupled%lam_pr + del_lam))*args_coupled%CA
  
  CALL axialtorstrain(args_coupled%eps(3,2), args_coupled%eps(3,3), C, args_coupled%lam_pr + del_lam, args_coupled%eps_m_pr &
  + alph_m*D_m_c, args_coupled%eps_p_pr + delta*xi_p, args_coupled%cmv, eps)
  
  eps_e = eps - (args_coupled%lam_pr + del_lam)*(args_coupled%eps_m_pr + alph_m*D_m_c) - (args_coupled%eps_p_pr + delta*xi_p)
  
  CALL mat2voigt(eps_e, v)
  CALL voigt2mat(MATMUL(C, v), sig)
  CALL dev(sig, sigdev)

  dlamsign = SIGN(1., dlamcr_c)
  dlamcr_c = dlamsign*(args_coupled%dlamcr0 + args_coupled%flp*(args_coupled%r_pr + delta*xi_r)**args_coupled%nlp)

  dmc = args_coupled%dmc0*SQRT(2./3.) + args_coupled%fmp*(args_coupled%r_pr + delta*xi_r)**args_coupled%nmp


  ! transformation surfaces
  
  CALL mat2voigt(args_coupled%eps_m_pr + alph_m*D_m_c, v)
  CALL voigt2mat(MATMUL(args_coupled%Ci, v), t)
  CALL surf_init(t, args_coupled%a1i, args_coupled%a2i, args_coupled%bi, gi, dgi)
  CALL voigt2mat(MATMUL(args_coupled%Cs, v), t)
  CALL surf_sat((args_coupled%lam_pr + del_lam)*t, args_coupled%a1s, args_coupled%a2s, args_coupled%bs, args_coupled%ks, &
  args_coupled%ms, gs_, Gs, dGs)

  CALL surf_coupled((args_coupled%lam_pr + del_lam)*(args_coupled%eps_m_pr + alph_m*D_m_c), args_coupled%eps_p_pr + delta*xi_p, &
  args_coupled%kc, args_coupled%mc, Gc, dGcxi, dGceta)

  if (dlamsign == 1) then
    lcv = 1
  else
    lcv = 0
  end if


  ! plastic variables and updates

  CALL mat2voigt(args_coupled%a_pr + delta*xi_a, v)
  CALL voigt2mat(MATMUL(args_coupled%C_a, v), t)
  d_a = -2./3.*t

  d_eps_p = sigdev - dGceta


  ! grad(h) wrt del_U_p

  CALL mat2voigt(xi_p + xi_a, v1)
  CALL voigt2mat(MATMUL(args_coupled%m_, v1), t1)
  CALL mat2voigt(xi_p - xi_a, v2)
  CALL voigt2mat(MATMUL(args_coupled%M, v2), t2)

  dh_dxi_p = t1 + t2
  dh_dxi_a = t1 - t2
  dh_dxi_r = args_coupled%e*(xi_r - args_coupled%b)


  ! scalar product of plastic strain increments with grad_h

  CALL scalar22(xi_p, dh_dxi_p, s1)
  CALL scalar22(xi_a, dh_dxi_a, s2)
  xi_dot_grad_h = s1 + s2 + xi_r*dh_dxi_r


  ! formulate residuals

  zz = 1

  if (args_coupled%flag(1) == 1) then
    if (args_coupled%flag(2) == 0 .and. args_coupled%lam_pr + del_lam /= 0) then
      CALL mag(sigdev - dGs - dGcxi - Lam_i*dgi, s)
      f(zz) = dmc - s
    else if (args_coupled%flag(3) == 1 .or. (args_coupled%flag(2) == 0 .and. args_coupled%lam_pr + del_lam == 0)) then
      f(zz) = gi 
    else
      CALL scalar22(D_m_c,D_m_c,s)
      f(zz) = s - 1
    end if
    zz = zz + 1
  end if

  if (args_coupled%flag(2) == 1) then
    f(zz) = dmc*D_m_c(2,2) - sigdev(2,2) + Lam_i*dgi(2,2) + dGs(2,2) + dGcxi(2,2)
    f(zz + 1) = dmc*D_m_c(3,3) - sigdev(3,3) + Lam_i*dgi(3,3) + dGs(3,3) + dGcxi(3,3) 
    f(zz + 2) = dmc*D_m_c(2,1) - sigdev(2,1) + Lam_i*dgi(2,1) + dGs(2,1) + dGcxi(2,1)
    f(zz + 3) = dmc*D_m_c(3,1) - sigdev(3,1) + Lam_i*dgi(3,1) + dGs(3,1) + dGcxi(3,1)
    f(zz + 4) = dmc*D_m_c(3,2) - sigdev(3,2) + Lam_i*dgi(3,2) + dGs(3,2) + dGcxi(3,2)
    zz = zz + 5
  end if

  if (args_coupled%flag(3) == 1) then
    if (args_coupled%flag(1) == 1) then
      CALL scalar22(D_m_c, D_m_c, s)
      f(zz) = s - 1
    else
      CALL mag(sigdev - dGs - dGcxi - Lam_i*dgi, s)
      f(zz) = dmc - s
    end if
    zz = zz + 1
  end if

  if (args_coupled%flag(4) == 1) then
    CALL scalar22(sig, args_coupled%eps_m_pr + alph_m*D_m_c, s1)
    CALL mat2voigt(eps_e, v)
    CALL voigt2mat(MATMUL(args_coupled%CM - args_coupled%CA, v), t)
    CALL scalar22(eps_e, t, s2)
    CALL scalar22(args_coupled%eps_m_pr + alph_m*D_m_c, dGs, s3)
    CALL scalar22(args_coupled%eps_m_pr + alph_m*D_m_c, dGcxi, s4)
    f(zz) = - s1 + 0.5*s2 - (args_coupled%cpM - args_coupled%cpA)*args_coupled%theta*LOG(args_coupled%theta/args_coupled%T0) &
    + args_coupled%L*(args_coupled%theta - args_coupled%T0)/args_coupled%T0 + args_coupled%fll*(args_coupled%lam_pr &
    + del_lam)**args_coupled%nll + s3 + s4 + Lam_lam + dlamcr_c
    zz = zz + 1
  end if

  if (args_coupled%flag(5) == 1 .and. args_coupled%flag(4) == 1) then
    f(zz) = args_coupled%lam_pr + del_lam - lcv
    zz = zz + 1
  else if (args_coupled%flag(5) == 1) then
    CALL scalar22(sig, args_coupled%eps_m_pr + alph_m*D_m_c, s1)
    CALL mat2voigt(eps_e, v)
    CALL voigt2mat(MATMUL(args_coupled%CM - args_coupled%CA, v), t)
    CALL scalar22(eps_e, t, s2)
    CALL scalar22(args_coupled%eps_m_pr + alph_m*D_m_c, dGs, s3)
    CALL scalar22(args_coupled%eps_m_pr + alph_m*D_m_c, dGcxi, s4)
    f(zz) = - s1 + 0.5*s2 - (args_coupled%cpM - args_coupled%cpA)*args_coupled%theta*LOG(args_coupled%theta/args_coupled%T0) &
    + args_coupled%L*(args_coupled%theta - args_coupled%T0)/args_coupled%T0 + args_coupled%fll*(args_coupled%lam_pr &
    + del_lam)**args_coupled%nll + s3 + s4 + Lam_lam + dlamcr_c
    zz = zz + 1
  end if

  if (args_coupled%flag(6) == 1) then
    f(zz) = args_coupled%eqval*dh_dxi_p(2,2) - xi_dot_grad_h*d_eps_p(2,2)
    f(zz + 1) = args_coupled%eqval*dh_dxi_p(3,3) - xi_dot_grad_h*d_eps_p(3,3)
    f(zz + 2) = args_coupled%eqval*dh_dxi_p(3,2) - xi_dot_grad_h*d_eps_p(3,2)
    f(zz + 3) = args_coupled%eqval*dh_dxi_p(3,1) - xi_dot_grad_h*d_eps_p(3,1)
    f(zz + 4) = args_coupled%eqval*dh_dxi_p(2,1) - xi_dot_grad_h*d_eps_p(2,1)
    f(zz + 5) = args_coupled%eqval*dh_dxi_a(2,2) - xi_dot_grad_h*d_a(2,2)
    f(zz + 6) = args_coupled%eqval*dh_dxi_a(3,3) - xi_dot_grad_h*d_a(3,3)
    f(zz + 7) = args_coupled%eqval*dh_dxi_a(3,2) - xi_dot_grad_h*d_a(3,2)
    f(zz + 8) = args_coupled%eqval*dh_dxi_a(3,1) - xi_dot_grad_h*d_a(3,1)
    f(zz + 9) = args_coupled%eqval*dh_dxi_a(2,1) - xi_dot_grad_h*d_a(2,1)
    f(zz + 10) = args_coupled%eqval*dh_dxi_r + xi_dot_grad_h*args_coupled%b*args_coupled%Q*(args_coupled%r_pr + delta*xi_r)
    CALL mat2voigt(xi_p + xi_a, v1)
    CALL voigt2mat(MATMUL(args_coupled%m_, v1), t1)
    CALL scalar22(xi_p + xi_a, t1, s1)
    CALL mat2voigt(xi_p - xi_a, v2)
    CALL voigt2mat(MATMUL(args_coupled%M, v2), t2)
    CALL scalar22(xi_p - xi_a, t2, s2)
    f(zz + 11) = 2 - (s1 + s2 + args_coupled%e*(xi_r - args_coupled%b)**2)
  end if

end subroutine resid_coupled


subroutine resid_yield(lx, lf, x, f, args_yield)

  use subs
  use data_types

  implicit none

  INTEGER, INTENT(IN) :: lf, lx
  REAL, INTENT(IN) :: x(lx)
  REAL, INTENT(OUT) :: f(lf)
  TYPE(resid_yield_args), INTENT(IN) :: args_yield
  
  REAL :: dh_dxi1(3,3), dh_dxi2(3,3), dh_dxi3, lambda, s1, s2, t1(3,3), t2(3,3), v1(6), v2(6), xi1(3,3), xi2(3,3), xi3

  CALL voigt2mat((/-(x(1) + x(2)), x(1:5)/), xi1)
  CALL voigt2mat((/-(x(6) + x(7)), x(6:10)/), xi2)
  xi3 = x(11)
  lambda = x(12)

  CALL mat2voigt(xi1 + xi2, v1)
  CALL mat2voigt(xi1 - xi2, v2)
  CALL voigt2mat(MATMUL(args_yield%m_, v1), t1)
  CALL voigt2mat(MATMUL(args_yield%M, v2), t2)

  dh_dxi1 = t1 + t2
  dh_dxi2 = t1 - t2
  dh_dxi3 = args_yield%e*(xi3 - args_yield%b)

  f(1) = lambda*dh_dxi1(2,2) - args_yield%d_eps_p(2,2)
  f(2) = lambda*dh_dxi1(3,3) - args_yield%d_eps_p(3,3)
  f(3) = lambda*dh_dxi1(3,2) - args_yield%d_eps_p(3,2)
  f(4) = lambda*dh_dxi1(3,1) - args_yield%d_eps_p(3,1)
  f(5) = lambda*dh_dxi1(2,1) - args_yield%d_eps_p(2,1)
  f(6) = lambda*dh_dxi2(2,2) - args_yield%d_a(2,2)
  f(7) = lambda*dh_dxi2(3,3) - args_yield%d_a(3,3)
  f(8) = lambda*dh_dxi2(3,2) - args_yield%d_a(3,2)
  f(9) = lambda*dh_dxi2(3,1) - args_yield%d_a(3,1)
  f(10) = lambda*dh_dxi2(2,1) - args_yield%d_a(2,1)
  f(11) = lambda*dh_dxi3 - args_yield%d_r
  CALL scalar22(xi1 + xi2, t1, s1)
  CALL scalar22(xi1 - xi2, t2, s2)
  f(12) = 2 - s1 + s2 + args_yield%e*(xi3 - args_yield%b)**2

end subroutine resid_yield