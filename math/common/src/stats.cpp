#include <opa/math/common/stats.h>

OPA_NAMESPACE(opa, math, common)

bignum nchoosek(int n, int k) {
  bignum res = 1;
  if (k < 0 || k > n)
    return 0;
  REP (i, k)
    res = res * (n - i);
  REP (i, k)
    res = res / (i + 1);
  return res;
}

Float erf_binomial(const Float &prob, const Float &nsamples,
                   const Float &obs_below) {
  Float var = prob * (Float(1) - prob);
  Float mean = prob;
  return FloatUtil::erf<Float>((obs_below - mean * nsamples) /
                               std::sqrt(var * nsamples)) /
           2 +
         0.5;
}

Float erf_binomial_small(const Float &prob, int nsamples, int obs_below) {
  Float res = 0;
  REP (i, obs_below + 1) {
    res += nchoosek(nsamples, i).getfloat() * prob.pow(i) *
           (Float(1) - prob).pow(nsamples - i);
  }
  return res;
}

Float erf_at_conf(const Float &prob, const Float &conf) { OPA_CHECK0(false); }

OPA_NAMESPACE_END(opa, math, common)
