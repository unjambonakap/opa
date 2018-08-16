#pragma once

#include <opa_common.h>
#include <opa/utils/csv.h>
#include <opa/dsp/sym/base.h>
#include <opa/math/common/FFT.h>
#include <opa/math/common/float.h>
#include <opa/math/common/RealField.h>

OPA_NAMESPACE(opa, dsp)

template <typename T> using TFunc = FunctionSptr<T, T>;
template <typename T> using CTFunc = FunctionSptr<Complex<T>, T>;

template <typename T> Complex<T> ang_to_complex(const T &ang) {
  return Complex<T>(cos(ang), sin(ang));
}

template <typename R, typename A> class IntegralComputer : public Function<R> {
public:
  struct IntegralParams {
    A from;
    A to;
    int n;
    FunctionSptr<R, A> a;
  };

  virtual void init(const IntegralParams &params) { m_params = params; }

  virtual R eval() const override {
    A cur = m_params.from;
    A step = (m_params.to - m_params.from) / (m_params.n - 1);
    R res = R(0);
    REP (i, m_params.n) {
      res += m_params.a->eval(cur);
      cur += step;
    }
    return res / R(m_params.n);
  }

private:
  IntegralParams m_params;
};

/*
template <typename T>
std::shared_ptr<typename std::remove_pointer<
  typename std::result_of<decltype(&T::base)(T)>::type>::type>
to_base(std::shared_ptr<T> a) {
  return std::static_pointer_cast<
    T, typename std::remove_pointer<
         typename std::result_of<decltype(&T::base)(T)>::type>::type>(a);
}
*/

template <typename T> FunctionSptr<T, T> DelayFunction(const T &t) {
  return FuncFromFunc<T, T>::Make([t](const T &a) { return a - t; });
}
template <typename T> FunctionSptr<T, T> RevFunction() {
  return FuncFromFunc<T, T>::Make([](const T &a) { return -a; });
}

template <typename T> CTFunc<T> ExpFunc(const T &a, const T &f, const T &phi) {
  T w = OPA_MATH::FloatUtil::get_pi<T>() * 2 * f;
  return FuncFromFunc<Complex<T>, T>::Make([a, w, phi](const T &t) {
    return a * OPA_MATH::FloatUtil::iexp(w * t + phi);
  });
}

template <typename T> CTFunc<T> ConstFunction(const T &c) {
  return FuncFromFunc<Complex<T>, T>::Make(
    [c](const T &t) { return Complex<T>(c); });
}

template <typename T> CTFunc<T> WhiteNoise(const T &e, const T &n) {
  return FuncFromFunc<Complex<T>, T>::Make([e, n](const T &t) {
    return Complex<T>(OPA_MATH::FloatUtil::get_gaussian<T>(e, n),
                      OPA_MATH::FloatUtil::get_gaussian<T>(e, n));
  });
}

template <typename T>
CTFunc<T> StepFunction(const std::vector<std::pair<T, Complex<T> > > &tb) {
  auto tmp = tb;
  sort(ALL(tmp));
  return FuncFromFunc<Complex<T>, T>::Make([tmp](const T &t) {
    int pos = std::lower_bound(ALL(tmp), MP(t, Complex<T>(0))) - tmp.begin();
    if (pos == tmp.size())
      return tmp.back().ND;
    if (pos == 0)
      return tmp[0].ND;

    pos -= (t < (tmp[pos].ST + tmp[pos - 1].ST) / 2);
    return tmp[pos].ND;
  });
}

template <typename T> CTFunc<T> RRC(const T &ts, const T &rolloff) {
  return FuncFromFunc<Complex<T>, T>::Make([ts, rolloff](const T &t) {
    T res;
    T pi = OPA_MATH::FloatUtil::get_pi<T>();
    if (OPA_MATH::FloatUtil::eq<T>(t, 0, 1e-9)) {
      res = (T(1) - rolloff + T(4) * rolloff / pi) / sqrt(ts);
    } else {
      T tmp = std::abs((t * rolloff * 4 / ts));
      if (!OPA_MATH::FloatUtil::eq<T>(rolloff, 0, 1e-9) &&
          OPA_MATH::FloatUtil::eq<T>(tmp, 1, 1e-9)) {
        tmp = pi / rolloff / 4;
        res = (T(1) + T(2) / pi) * sin(tmp) + (T(1) - T(2) / pi) * cos(tmp);
        res *= rolloff / sqrt(rolloff * 2);
      } else {
        res = (sin(pi * t / ts * (T(1) - rolloff)) +
               rolloff * 4 * t / ts * cos(pi * t / ts * (T(1) + rolloff))) /
              sqrt(ts) /
              (t / ts * pi * (T(1) - pow((rolloff * 4 * t / ts), 2)));
      }
    }
    return Complex<T>(res);
  });
}

template <typename T>
std::vector<Complex<T> > fast_correl(const CTFunc<T> &a, const CTFunc<T> &b,
                                     const Range<T> &range) {
  int npw = log2_high_bit(range.tb().size() - 1) + 1;
  int n = 1 << npw;

  OPA_MATH::FFT2<Complex<T> > fft;
  OPA_MATH::RealField<Complex<T> > r;
  fft.init(&r, npw, OPA_MATH::FloatUtil::iexp<T>(
                      (OPA_MATH::FloatUtil::get_pi<T>() * 2 / n)));
  std::vector<Complex<T> > d1, d2;
  for (auto x : range.tb()) {
    d1.pb(a->eval(x));
    d2.pb(std::conj(b->eval(x)));
  }
  OPA_DISP0(d1.size(), n, npw, range.tb().size());
  auto r1 = fft.fft(d1);
  auto r2 = fft.fft(d2);

  std::vector<Complex<T> > t1;
  REP (i, r1.size())
    t1.pb(r1[i] * std::conj(r2[i]));

  auto res = fft.ifft(t1);
  return res;
}

template <typename T>
CTFunc<T> Correlation(const CTFunc<T> &a, const CTFunc<T> &b, const T &window,
                      int npoints) {
  npoints |= 1;
  return FuncFromFunc<Complex<T>, T>::Make([a, b, window, npoints](const T &t) {
    Complex<T> res;
    Range<T> range;
    range.build_range(-window / 2, window / 2, npoints);
    for (auto &pos : range.tb()) {
      res += a->eval(pos) * std::conj(b->eval(t + pos));
    }
    return res;
  });
}

template <typename T, typename Stream>
opa::utils::CsvRecordWriterSptr<std::pair<T, Complex<T> >, Stream>
complex_dumper() {
  auto x =
    new opa::utils::CsvRecordWriterFromFunc<std::pair<T, Complex<T> >, Stream>(
      { "t", "x", "y", "ang", "mag" },
      [](Stream &s, const std::pair<T, Complex<T> > &a) {
        Complex<T> res = a.ND;
        s << opa::utils::Join(",", a.ST, res.real(), res.imag(), arg(res),
                              std::abs(res));
      });
  return opa::utils::CsvRecordWriterSptr<std::pair<T, Complex<T> >, Stream>(x);
}
template <typename T, typename Stream>
opa::utils::CsvRecordWriterSptr<T, Stream> complex_dumper2() {
  auto x = new opa::utils::CsvRecordWriterFromFunc<T, Stream>(
    { "x", "y", "ang", "mag" }, [](Stream &s, const T &res) {
      s << opa::utils::Join(",", res.real(), res.imag(), arg(res),
                            std::abs(res));
    });
  return opa::utils::CsvRecordWriterSptr<T, Stream>(x);
}

template <typename R, typename T, typename Stream>
void to_csv(Stream &s, const std::vector<FunctionSptr<R, T> > &funcs,
            const Range<T> &r) {
  typedef std::pair<T, R> DataType;
  auto record_writer =
    std::make_shared<opa::utils::AggregateCsvRecordWriter<DataType, Stream> >();
  record_writer->init(complex_dumper<T, Stream>(), funcs.size());

  auto writer = opa::utils::CsvWriter<std::vector<DataType>, Stream>();
  writer.init(record_writer, &s);

  for (auto &x : r.tb()) {
    std::vector<DataType> data;
    for (auto &f : funcs)
      data.pb(MP(x, f->eval(x)));
    writer.add(data);
  }
}

template <typename T, typename Stream>
void to_csv(Stream &s, const std::vector<T> &tb) {
  auto record_writer = complex_dumper2<T, Stream>();
  auto writer = opa::utils::CsvWriter<T, Stream>();
  writer.init(record_writer, &s);
  for (auto &x : tb)
    writer.add(x);
}

OPA_NAMESPACE_END(opa, dsp)
