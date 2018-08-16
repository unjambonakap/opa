#pragma once

#include <opa_common.h>
#include <opa/dsp/sym/base.h>
#include <opa/dsp/sym/dsp_base.h>

OPA_NAMESPACE(opa, dsp)

template <class T>
CTFunc<T> cdma_func(CDMAEncoder &encoder, int nsym, const T &chip_freq,
                    CTFunc<T> *debug, CTFunc<T> *co_func) {
  CTFunc<T> res = ConstFunction<T>(0);

  T chip_period = T(1) / chip_freq;
  T delay = 0;
  T rolloff = 0.5;
  std::vector<std::pair<T, Complex<T> > > tb;
  while (encoder.data_provider()->sym_count() < nsym) {
    double sym = encoder.get_mapped(encoder.get_next());
    tb.pb(MP(delay, Complex<T>(sym)));
    auto cur = compose(RRC<T>(chip_period, rolloff) * ConstFunction<T>(sym),
                       DelayFunction<T>(delay));
    res = res + cur;
    delay += chip_period;
  }
  if (debug) *debug = StepFunction(tb);
  if (co_func) {
    *co_func = res * ExpFunc<T>(1, chip_freq / 11, 0);
  }
  return res;
}

OPA_NAMESPACE_END(opa, dsp)
