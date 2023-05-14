module data_types

type, public :: resid_yield_args
    REAL :: d_eps_p(3,3), d_a(3,3), d_r, M(6,6), m_(6,6), b, e
end type resid_yield_args

type, public :: resid_coupled_args
    REAL :: eps(3,3), theta, eps_p_pr(3,3), r_pr, a_pr(3,3), eps_m_pr(3,3), lam_pr, bi, a1i, a2i, Ci(6,6), bs, ks, ms, a1s, a2s, &
    Cs(6,6), kc, mc, CA(6,6), CM(6,6), cpA, cpM, T0, L, dmc0, dlamcr0, dlamcr, D_m(3,3), flp, fmp, fll, nlp, nmp, nll, &
    C_a(6,6), M(6,6), m_(6,6), b, Q, e, eqval, xi0(12)
    INTEGER :: flag(6)
    LOGICAL :: cmv
end type resid_coupled_args

type, public :: cell
    REAL, DIMENSION(:), ALLOCATABLE :: component
end type cell

end module