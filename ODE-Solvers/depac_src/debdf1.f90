MODULE DEBDF1
  !USE service, ONLY : SP
    use mod_kinds, only : i4,sp
    use omp_lib
  IMPLICIT NONE
  REAL(sp) :: told_com, conit_com, crate_com, hold_com, rc_com, rmax_com, el0_com, &
    h_com, hmin_com, hmxi_com, hu_com, tn_com, uround_com
  !$OMP THREADPRIVATE(told_com, conit_com, crate_com, hold_com, rc_com, rmax_com, el0_com, &
  !$OMP& h_com, hmin_com, hmxi_com, hu_com, tn_com, uround_com)
  REAL(sp) :: el_com(13), elco_com(13,12), tesco_com(3,12)
  !$OMP THREADPRIVATE(el_com,elco_com,tesco_com)
  INTEGER(i4) :: iquit_com, init_com, iyh_com, iewt_com, iacor_com, isavf_com, iwm_com, &
    ksteps_com, ialth_com, ipup_com, lmax_com, meo_com, nqnyh_com, nstepj_com, &
    ibegin_com, itol_com, iinteg_com, itstop_com, ijac_com, iband_com, ier_com, &
    jstart_com, kflag_com, l_com, meth_com, miter_com, maxord_com, n_com, nq_com, &
    nst_com, nfe_com, nje_com, nqu_com
  !$OMP THREADPRIVATE(iquit_com, init_com, iyh_com, iewt_com, iacor_com, isavf_com, iwm_com, &
  !$OMP&  ksteps_com, ialth_com, ipup_com, lmax_com, meo_com, nqnyh_com, nstepj_com, &
  !$OMP&  ibegin_com, itol_com, iinteg_com, itstop_com, ijac_com, iband_com, ier_com, &
  !$OMP&  jstart_com, kflag_com, l_com, meth_com, miter_com, maxord_com, n_com, nq_com, &
  !$OMP&  nst_com, nfe_com, nje_com, nqu_com)
END MODULE DEBDF1
