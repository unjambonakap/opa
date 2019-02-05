/************************
 *  Extra utils for NR_FP
 ************************/

#ifndef FPLLL_NR_FP_MISC_H
#define FPLLL_NR_FP_MISC_H

FPLLL_BEGIN_NAMESPACE


/* -------------------------
 *   set_z (Z_NR --> FP_NR)
 * ------------------------- */


/* set_z (to double) */
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}
#endif

template<> template<>
inline void FP_NR<double>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = a.get_d();
}


/* set_z (to long double) */
#ifdef FPLLL_WITH_LONG_DOUBLE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}
#endif

template<> template<>
inline void FP_NR<long double>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = a.get_ld();
}

#endif


/* set_z (to dpe_t) */
#ifdef FPLLL_WITH_DPE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<long>& a, mp_rnd_t /*rnd*/) {
  dpe_set_si(data, a.get_data());
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  dpe_set_d(data, a.get_data());
}
#endif

template<> template<>
inline void FP_NR<dpe_t>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  dpe_set_z(data, const_cast<mpz_t&>(a.get_data()));
}

#endif


/* set_z (to dd_real) */
#ifdef FPLLL_WITH_QD

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dd_real>::set_z(const Z_NR<long>& a, mp_rnd_t) {
  data = a.get_data();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dd_real>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_data();
}
#endif

template<> template<>
inline void FP_NR<dd_real>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = mpz_get_d(a.get_data());
}

#endif


/* set_z (to qd_real) */
#ifdef FPLLL_WITH_QD

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<qd_real>::set_z(const Z_NR<long>& a, mp_rnd_t) {
  data = a.get_data();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<qd_real>::set_z(const Z_NR<double>& a, mp_rnd_t /*rnd*/) {
  data = a.get_data();
}
#endif

template<> template<>
inline void FP_NR<qd_real>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t /*rnd*/) {
  data = mpz_get_d(a.get_data());
}

#endif


/* set_z (to mpfr_t) */
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<long>& a, mp_rnd_t rnd) {
  mpfr_set_si(data, a.get_data(), rnd);
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<double>& a, mp_rnd_t rnd) {
  mpfr_set_d(data, a.get_data(), rnd);
}
#endif

template<> template<>
inline void FP_NR<mpfr_t>::set_z(const Z_NR<mpz_t>& a, mp_rnd_t rnd) {
  mpfr_set_z(data, a.get_data(), rnd);
}


/* -----------------------------
 *   get_z_exp (FP_NR --> Z_NR)
 * ----------------------------- */


/* get_z_exp_we (double --> Z_NR) */
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<long>& a, long& expo, long expo_add) const {
  expo = 0;
  a = static_cast<long>(ldexp(data, expo_add));
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<double>& a, long& expo, long expo_add) const {
  expo = 0;
  a.get_data() = trunc(ldexp(data, expo_add));
}
#endif

template<> template<>
inline void FP_NR<double>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expo_add) const {
  expo = max(exponent() + expo_add - numeric_limits<double>::digits, 0L);
  mpz_set_d(a.get_data(), ldexp(data, expo_add - expo));
}

template<> template<class Z>
inline void FP_NR<double>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}


/* get_z_exp_we (long double --> Z_NR) */
#ifdef FPLLL_WITH_LONG_DOUBLE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<long>& a, long& expo, long expo_add) const {
  expo = 0;
  a = static_cast<long>(ldexpl(data, expo_add));
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<double>& a, long& expo, long expo_add) const {
  expo = 0;
  a.get_data() = trunc(static_cast<double>(ldexpl(data, expo_add)));
}
#endif

template<> template<>
inline void FP_NR<long double>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expo_add) const {
  expo = max(exponent() + expo_add - numeric_limits<long double>::digits, 0L);
  /* If expo > 0, then
     expo_add - expo = numeric_limits<long double>::digits - exponent()
     which implies that ldexpl(data, expo_add - expo) is an integer */
  LDConvHelper::mpz_set_ld(a.get_data(), trunc(ldexpl(data, expo_add - expo)));
}

template<> template<class Z>
inline void FP_NR<long double>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}

#endif


/* get_z_exp_we (dpe_t --> Z_NR) */
#ifdef FPLLL_WITH_DPE

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<long>& a, long& expo) const {
  expo = 0;
  a = get_si();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<double>& a, long& expo) const {
  expo = max(DPE_EXP(data) - numeric_limits<double>::digits, 0);
  a.get_data() = trunc(ldexp(DPE_MANT(data), DPE_EXP(data) - expo));
}
#endif

template<> template<>
inline void FP_NR<dpe_t>::get_z_exp(Z_NR<mpz_t>& a, long& expo) const {
  expo = max(DPE_EXP(data) - DPE_BITSIZE, 0);
  mpz_set_d(a.get_data(), trunc(ldexp(DPE_MANT(data), DPE_EXP(data) - expo)));
}

template<> template<class Z>
inline void FP_NR<dpe_t>::get_z_exp_we(Z_NR<Z>& a, long& expo, long /*expo_add*/) const {
  return get_z_exp(a, expo);
}

#endif


/* get_z_exp and get_z_exp_we (dd_real and qd_real --> Z_NR) */
#ifdef FPLLL_WITH_QD

template<> template<>
inline void FP_NR<dd_real>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expo_add) const {
  // double-double has almost the same exp as double
  expo = max(exponent() + expo_add - numeric_limits<double>::digits, 0L);
  // losing precision here
  mpz_set_d(a.get_data(),::to_double(::ldexp(data, expo_add - expo)));
}
template<> template<>
inline void FP_NR<qd_real>::get_z_exp_we(Z_NR<mpz_t>& a, long& expo, long expo_add) const {
  expo = max(exponent() + expo_add - numeric_limits<double>::digits, 0L);
  mpz_set_d(a.get_data(),::to_double(::ldexp(data, expo_add - expo)));
}

#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<dd_real>::get_z_exp_we(Z_NR<long>& a, long& expo, long expo_add) const {
  expo = 0;
  a = ::to_int(::ldexp(data, expo_add));
}
template<> template<>
inline void FP_NR<qd_real>::get_z_exp_we(Z_NR<long>& a, long& expo, long expo_add) const {
  expo = 0;
  a = ::to_int(::ldexp(data, expo_add));
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<dd_real>::get_z_exp_we(Z_NR<double>& a, long& expo, long expo_add) const {
  expo = 0;
  a = get_si();
}
template<> template<>
inline void FP_NR<qd_real>::get_z_exp_we(Z_NR<double>& a, long& expo, long expo_add) const {
  expo = 0;
  a = get_si();
}
#endif

template<> template<class Z>
inline void FP_NR<dd_real>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}
template<> template<class Z>
inline void FP_NR<qd_real>::get_z_exp(Z_NR<Z>& a, long& expo) const {
  return get_z_exp_we(a, expo, 0);
}

#endif


/* get_z_exp_we (from mpfr_t) */
#ifdef FPLLL_WITH_ZLONG
template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<long>& a, long& expo) const {
  expo = 0;
  a = get_si();
}
#endif

#ifdef FPLLL_WITH_ZDOUBLE
template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<double>& a, long& expo) const {
  expo = 0;
  a.get_data() = trunc(mpfr_get_d(data, GMP_RNDZ));
}
#endif

template<> template<>
inline void FP_NR<mpfr_t>::get_z_exp(Z_NR<mpz_t>& a, long& expo) const {
  expo = mpfr_get_z_exp(a.get_data(), data);
  if (expo < 0) {
    a.mul_2si(a, expo);
    expo = 0;
  }
}

template<> template<class Z>
inline void FP_NR<mpfr_t>::get_z_exp_we(Z_NR<Z>& a, long& expo, long /*expo_add*/) const {
  return get_z_exp(a, expo);
}


FPLLL_END_NAMESPACE

#endif
