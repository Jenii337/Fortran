module subs

use diag3x3

contains


  subroutine axialtorstrain(eps23, eps33, C, lam, eps_m, eps_p, cmv, eps)

    implicit none

    REAL, INTENT(IN)  :: eps23, eps33, C(6,6), lam, eps_m(3,3), eps_p(3,3)
    LOGICAL, INTENT(IN) :: cmv
    REAL, INTENT(OUT) :: eps(3,3)
    REAL :: eps11, eps12, eps13, eps22, eps33_c

    if (.not. cmv) then
      eps33_c = lam*eps_m(3,3) + eps_p(3,3)
    else
      eps33_c = eps33
    end if

    eps11 = (C(1,2)*C(1,3) - C(1,1)*C(1,3))/(C(1,1)**2 - C(1,2)**2)*(eps33_c - lam*eps_m(3,3) - eps_p(3,3)) + lam*eps_m(1,1) &
    + eps_p(1,1)

    eps22 = (C(1,2)*C(1,3) - C(1,1)*C(1,3))/(C(1,1)**2 - C(1,2)**2)*(eps33_c - lam*eps_m(3,3) - eps_p(3,3)) + lam*eps_m(2,2) &
    + eps_p(2,2)

    eps12 = lam*eps_m(1,2) + eps_p(1,2)

    eps13 = lam*eps_m(1,3) + eps_p(1,3)

    eps(1,1) = eps11
    eps(1,2) = eps12
    eps(1,3) = eps13
    eps(2,1) = eps12
    eps(2,2) = eps22
    eps(2,3) = eps23
    eps(3,1) = eps13
    eps(3,2) = eps23
    eps(3,3) = eps33_c

  end subroutine axialtorstrain
  

  subroutine dev(tens2, devtens2)

  implicit none
  
  REAL, INTENT(IN) :: tens2(3,3)
  REAL, INTENT(OUT) :: devtens2(3,3)
  INTEGER :: ii

  devtens2 = tens2

  do ii = 1, 3
    devtens2(ii,ii) = tens2(ii,ii) - 1./3.*(tens2(1,1) + tens2(2,2)+ tens2(3,3))
  end do

  end subroutine dev


  subroutine dyad11(vecin1, vecin2, tens2out)

    implicit none
    
    REAL, INTENT(IN) :: vecin1(3), vecin2(3)
    REAL, INTENT(OUT) :: tens2out(3,3)
    INTEGER :: ii, jj

    tens2out (:,:) = 0
    do ii = 1, 3
      do jj = 1, 3
        tens2out(jj,ii) = vecin1(ii)*vecin2(jj)
      end do
    end do

  end subroutine dyad11


  subroutine mag(tensor, scalout)

    implicit none
    
    REAL, INTENT(IN) :: tensor(3,3)
    REAL, INTENT(OUT) :: scalout

    call scalar22(tensor, tensor, scalout)
    scalout = SQRT(scalout)

  end subroutine mag


  subroutine mat2voigt(matin, vecout)

    implicit none
    
    REAL, INTENT(IN) :: matin(3,3)
    REAL, INTENT(OUT) :: vecout(6,1)

    vecout(1,1) = matin(1,1)
    vecout(2,1) = matin(2,2)
    vecout(3,1) = matin(3,3)
    vecout(4,1) = matin(3,2)
    vecout(5,1) = matin(3,1)
    vecout(6,1) = matin(2,1)

  end subroutine mat2voigt


  subroutine non0mag(arg, scalar)

    implicit none
    
    REAL, INTENT(IN) :: arg(3,3)
    REAL, INTENT(OUT) :: scalar

    CALL mag(arg, scalar)

    if (scalar == 0) then
      scalar = 1
    end if

  end subroutine non0mag


  subroutine non0norm(arg, scalar)

    implicit none
    
    REAL, INTENT(IN) :: arg
    REAL, INTENT(OUT) :: scalar

    scalar = ABS(arg)

    if (scalar == 0) then
      scalar = 1
    end if

  end subroutine non0norm


  subroutine pvalsderivs(dev2tensin, lam1, lam2, lam3, dlam1, dlam2, dlam3)

    implicit none
    
    REAL, INTENT(IN) :: dev2tensin(3,3)
    REAL, INTENT(OUT) :: lam1, lam2, lam3, dlam1(3,3), dlam2(3,3), dlam3(3,3)
    REAL :: w(3), q(3,3), tens(3,3)

    tens = dev2tensin

    CALL eigvec3x3(tens, w, q)

    lam1 = w(1)
    lam2 = w(2)
    lam3 = w(3)

    CALL dyad11(q(1,:), q(1,:), dlam1)
    CALL dyad11(q(2,:), q(2,:), dlam2)
    CALL dyad11(q(3,:), q(3,:), dlam3)

  end subroutine pvalsderivs


  subroutine scalar22(tens2in1, tens2in2, scalout)

    implicit none
    
    REAL, INTENT(IN) :: tens2in1(3,3), tens2in2(3,3)
    REAL, INTENT(OUT) :: scalout
    INTEGER :: ii, jj

    scalout = 0;

    do ii = 1, 3
      do jj = 1, 3
        scalout = scalout + tens2in1(ii,jj)*tens2in2(ii,jj)
      end do
    end do

  end subroutine scalar22


  subroutine surf_coupled(xi, eta, k, m, G, dGxi, dGeta)
  
    implicit none
    
    REAL, INTENT(IN) :: xi(3,3), eta(3,3), k, m
    REAL, INTENT(OUT) :: G, dGxi(3,3), dGeta(3,3)
    REAL :: s1, s2

    CALL mag(eta, s1)
    CALL scalar22(xi, eta, s2)

    G = -k*s1**m*s2

    if (s1 == 0) then 
      dGxi = dGxi*0
      dGeta = dGeta*0
    else
      dGxi = -k*s1**m*eta
      dGeta = -k*(xi*s1**m + s2*m*s1**(m-2)*eta)
    end if

  end subroutine surf_coupled
  

  subroutine surf_init(xi, a1, a2, b, g, dg)

    implicit none
    
    REAL, INTENT(IN) :: xi(3,3), a1, a2, b
    REAL, INTENT(OUT) :: g, dg(3,3)
    REAL :: dg_c(3,3), dlam1(3,3), dlam2(3,3), dlam3(3,3), lam1, lam2, lam3, s1, s2, s3

    CALL pvalsderivs(xi, lam1, lam2, lam3, dlam1, dlam2, dlam3)

    g = ((ABS(lam1) - a1*lam1)**a2 + (ABS(lam2) - a1*lam2)**a2 + (ABS(lam3) - a1*lam3)**a2 - b**a2)/b**a2

    CALL non0norm(lam1, s1)
    CALL non0norm(lam2, s2)
    CALL non0norm(lam3, s3)

    dg_c = (a2*(ABS(lam1) - a1*lam1)**(a2-1)*(lam1/s1 - a1)*dlam1 + a2*(ABS(lam2) - a1*lam2)**(a2-1)*(lam2/s2 - a1)*dlam2 &
    + a2*(ABS(lam3) - a1*lam3)**(a2-1)*(lam3/s3 - a1)*dlam3)/b**a2
    
    CALL dev(dg_c, dg)
    

  end subroutine surf_init
  

  subroutine surf_sat(xi, a1, a2, b, k, m, g_, G, dG)

    implicit none
    
    REAL, INTENT(IN) :: xi(3,3), a1, a2, b, k, m
    REAL, INTENT(OUT) :: g_, G, dG(3,3)
    REAL :: dg_(3,3), dlam1(3,3), dlam2(3,3), dlam3(3,3), lam1, lam2, lam3, gp, s ,s1, s2, s3, t(3,3)

    CALL pvalsderivs(xi, lam1, lam2, lam3, dlam1, dlam2, dlam3)

    g_ = ((ABS(lam1) - a1*lam1)**a2 + (ABS(lam2) - a1*lam2)**a2 + (ABS(lam3) - a1*lam3)**a2 - b**a2)/b**a2

    if (g_ > 1) then
      gp = LOG(g)
    else
      gp = 0
    end if

    G = k*gp**m

    CALL non0norm(lam1, s1)
    CALL non0norm(lam2, s2)
    CALL non0norm(lam3, s3)

    dg_ = (a2*(ABS(lam1) - a1*lam1)**(a2 - 1)*(lam1/s1 - a1)*dlam1 + a2*(ABS(lam2) - a1*lam2)**(a2 - 1)*(lam2/s2 - a1)*dlam2 &
    + a2*(ABS(lam3) - a1*lam3)**(a2 - 1)*(lam3/s3 - a1)*dlam3)/b**a2

    if (g > 1) then
      CALL non0norm(LOG(g_), s)
      CALL dev((m*k*gp**(m - 1))*dg_/2/g_*(LOG(g_)/s + 1), dG)
    else
      CALL non0norm(g_, s)
      CALL dev(dg_, t)
      dG = (m*k*gp**(m - 1)/2*(g_/s + 1))*t
    end if

  end subroutine surf_sat
  

  subroutine voigt2mat(vecin, matout)

    implicit none
    
    REAL, INTENT(IN) :: vecin(6)
    REAL, INTENT(OUT) :: matout(3,3)

    matout(1,1) = vecin(1)
    matout(2,2) = vecin(2)
    matout(3,3) = vecin(3)
    matout(3,2) = vecin(4)
    matout(3,1) = vecin(5)
    matout(2,1) = vecin(6)
    matout(1,2) = matout(2,1)
    matout(1,3) = matout(3,1)
    matout(2,3) = matout(3,2)

  end subroutine voigt2mat



end module subs



