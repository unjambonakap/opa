#include <opa/dsp/gps/l1ca.h>
#include <opa/dsp/reader.h>
#include <opa/dsp/sym/dsp_base.h>
#include <opa/math/common/FFT.h>
#include <opa/utils/clone.h>
#include <opa_common.h>

DEFINE_string(filename, "", "");
DEFINE_bool(costas, false, "");
DEFINE_bool(filter, false, "");
DEFINE_string(out_dir, "", "");

using namespace std;
using namespace opa::dsp;
using namespace opa::utils;
using namespace opa::dsp::gps;
using namespace opa::math::common;
typedef Complex<s8> T;
typedef complex<double> Type;
double samp_f = 8738133.333;

double angle_diff(double diff) {
  while (diff < -pi)
    diff += 2 * pi;
  while (diff >= pi)
    diff -= 2 * pi;
  return diff;
}

template <class T> class MemUnit {
public:
  void init(int n_) {
    this->n = max(n_, 1);
    pos = 0;
    tb.resize(n, 0);
  }

  const T &get(int p) const { return tb[(pos + n - p) % n]; }
  void set(const T &v) {
    pos = (pos + 1) % n;
    tb[pos] = v;
  }

  int n;
  int pos;
  vector<T> tb;
};

template <class T>
class Filter : public opa::utils::Initable,
               public opa::utils::Clonable<Filter<T> > {
public:
  virtual ~Filter() {}
  virtual T push_internal(const T &v) = 0;
  virtual void init(int mem_in, int mem_out) {
    opa::utils::Initable::init();
    in.init(mem_in);
    out.init(mem_out);
  }

  virtual T push(const T &v) final {
    in.set(v);
    T res = push_internal(v);
    out.set(res);
    return res;
  }

  inline T get_in(int id) const { return in.get(id); }
  inline T get_out(int id) const { return out.get(id); }

  MemUnit<T> in, out;
};
OPA_DECL_SPTR_TMPL(Filter, FilterSptr);

template <class T> class FIRFilter : public Filter<T> {
public:
  static FIRFilter<T> Delay(int n) {
    FIRFilter<T> filter;
    std::vector<T> coeffs(n, 0);
    coeffs.back() = 1;
    filter.init(coeffs);
    return filter;
  }

  void init(const std::vector<T> &coeffs_) {
    Filter<T>::init(coeffs_.size(), 0);
    this->coeffs = coeffs_;
    n = coeffs.size();
  }

  virtual T push_internal(const T &v) override {
    T res = 0;
    REP (i, n)
      res += this->get_in(i) * T(coeffs[i]);
    return res;
  }
  OPA_DECL_CLONE(Filter<T>, FIRFilter<T>);

  int n;
  std::vector<T> coeffs;
};

template <class T> class SecondOrderIIRFilter : public Filter<T> {
public:
  struct Params {
    std::vector<double> coeff_line;
    double scale;
  };

  virtual void init(const Params &params) {
    Filter<T>::init(0, 0);
    this->m_params = params;
    m_mem.init(2);
  }

  virtual T push_internal(const T &v) override {
    T res = 0;
    const auto &c = m_params.coeff_line;
    T tmp =
      v * m_params.scale * c[3] - m_mem.get(0) * c[4] - m_mem.get(1) * c[5];
    res = tmp * c[0] + m_mem.get(0) * c[1] + m_mem.get(1) * c[2];
    m_mem.set(tmp);
    return res;
  }
  OPA_DECL_CLONE(Filter<T>, SecondOrderIIRFilter<T>);

private:
  Params m_params;
  MemUnit<T> m_mem;
};

template <class T> class IIRFilter : public Filter<T> {
public:
  struct Params {
    std::vector<typename SecondOrderIIRFilter<T>::Params> lines;
  };

  virtual void init(const Params &params) {
    Filter<T>::init(0, 0);
    REP (i, params.lines.size()) {
      m_filters.emplace_back();
      m_filters.back().init(params.lines[i]);
    }
  }

  virtual T push_internal(const T &v) override {
    T res = v;
    for (auto &x : m_filters) {
      res = x.push(res);
    }
    return res;
  }

  OPA_DECL_CLONE(Filter<T>, IIRFilter<T>);

  Params m_params;
  vector<SecondOrderIIRFilter<T> > m_filters;
};

template <class T> class IIRFilter2 : public Filter<T> {
public:
  struct Params {
    std::vector<T> fir_coeffs;
    std::vector<T> iir_coeffs;
  };

  virtual void init(const Params &params) {
    m_params = params;
    Filter<T>::init(m_params.fir_coeffs.size(), m_params.iir_coeffs.size());
  }

  virtual T push_internal(const T &v) override {
    T res = 0;
    REP (i, m_params.fir_coeffs.size())
      res += this->get_in(i) * m_params.fir_coeffs[i];

    FOR (i, 1, m_params.iir_coeffs.size())
      res -= this->get_out(i - 1) * m_params.iir_coeffs[i];
    res /= m_params.iir_coeffs[0];
    return res;
  }

  OPA_DECL_CLONE(Filter<T>, IIRFilter2<T>);

  Params m_params;
};

class CostasLoop {
public:
  void init(double bw, bool mode) {
    this->bw = bw;
    damp = sqrt(2.) / 2.;
    // integrator filter (delay is in f1)
    f2.init({ { 1 }, { 1, -1 } });

    double alpha = 2 * 2;
    double beta = 2 * bw * damp * 2;
    double wo2 = bw * bw;
    double norm = alpha + beta;
    norm = 1;
    alpha /= norm;
    beta /= norm;
    wo2 /= norm;
    f1.init({ { 0, wo2, 2 * wo2, wo2 }, { beta + alpha, beta - alpha } });

    float loop_bw = bw;
    float denom = (1.0 + 2.0 * damp * loop_bw + loop_bw * loop_bw);
    d_alpha = (4 * damp * loop_bw) / denom;
    d_beta = (4 * loop_bw * loop_bw) / denom;
    m_mode = mode;

    {
      int nn = 10;
      vector<double> tb(nn + 1, 0);
      double fx = 1. / nn;
      tb[0] = fx;
      tb[nn] = -fx;

      derivative_filter.init({ tb });
    }
  }

  complex<double> push(complex<double> v) {
    err = angle_diff(2 * (std::arg(v) - cur_ang));

    double y;
    if (m_mode) {
      // Type tmp = v * FloatUtil::iexp<double>(-d_phase);
      // tmp /= max(1., std::abs(v));
      // err = tmp.real() * tmp.imag();
      y = d_phase;
      d_freq = d_freq + d_beta * err;
      d_phase = d_phase + d_freq + d_alpha * err;
    } else {
      y = f2.push(f1.push(err));
    }
    cur_ang = y;
    cur_w = derivative_filter.push(y);

    return v * FloatUtil::iexp<double>(-y);
  }

  bool m_mode;
  double d_freq = 0;
  double d_phase = 0;
  double d_alpha = 0;
  double d_beta = 0;

  double cur_w;
  double cur_ang;
  double err;

  double bw;
  double damp;
  FIRFilter<double> derivative_filter;
  IIRFilter2<double> f1, f2;
};

class ControlLoop {
public:
  void init(double bw) {
    this->bw = bw;
    damp = sqrt(2.) / 2.;
    // integrator filter (delay is in f1)

    float loop_bw = bw;
    float denom = (1.0 + 2.0 * damp * loop_bw + loop_bw * loop_bw);
    d_alpha = (4 * damp * loop_bw) / denom;
    d_beta = (4 * loop_bw * loop_bw) / denom;
  }

  double push(double err) {

    double y;
    y = d_phase;
    d_freq = d_freq + d_beta * err;
    d_phase = d_phase + d_freq + d_alpha * err;

    return y;
  }

  double d_freq = 0;
  double d_phase = 0;
  double d_alpha = 0;
  double d_beta = 0;
  double bw;
  double damp;
};

class GrcCostatsLoop {
public:
  void init(double bw) { m_loop.init(bw); }
  complex<double> push(const complex<double> &x) {
    complex<double> cur = FloatUtil::iexp<double>(-m_loop.d_phase) * x;
    double err = cur.real() * cur.imag();
    double mag = std::abs(cur);
    mag = max(1., mag);
    err /= mag;
    err = angle_diff(2 * (std::arg(cur)));
    m_loop.push(err);
    return cur;
  }

  ControlLoop m_loop;
};

FIRFilter<Type> lp_filter1;
IIRFilter<Type> iir_filter1;
void init_data() {

  vector<Type> coeffs = {
    0.0002904875165031577, 0.0006182507430321314, 0.0011320780290646976,
    0.0017400524253564794, 0.00229736517390444,   0.0025749917628237477,
    0.0022818979784289565, 0.001115521863788132,  -0.001163233829050091,
    -0.004640490054843419, -0.009163929475687577, -0.01428355143267104,
    -0.019238727485646673, -0.02300834519958391,  -0.024427693674656338,
    -0.022361304871480067, -0.015909662544951266, -0.004611596960447465,
    0.011395030372279407,  0.031294141145370005,  0.05361749479107545,
    0.07638390774594916,   0.09733756558534401,   0.11425045179170212,
    0.12524062952886045,   0.12905119980320923,   0.12524062952886045,
    0.11425045179170212,   0.09733756558534401,   0.07638390774594916,
    0.05361749479107545,   0.031294141145370005,  0.011395030372279407,
    -0.004611596960447465, -0.015909662544951266, -0.022361304871480067,
    -0.024427693674656338, -0.02300834519958391,  -0.019238727485646673,
    -0.01428355143267104,  -0.009163929475687577, -0.004640490054843419,
    -0.001163233829050091, 0.001115521863788132,  0.0022818979784289565,
    0.0025749917628237477, 0.00229736517390444,   0.0017400524253564794,
    0.0011320780290646976, 0.0006182507430321314, 0.0002904875165031577
  };

  lp_filter1.init(coeffs);
  iir_filter1.init(
    { { { { 1.0, 2.0, 1.0, 1.0, -1.9968534854432143, 0.9968633395607119 },
          2.4635293743901034e-06 },
        { { 1.0, 1.0, 0.0, 1.0, -0.9968633318334379, 0.0 },
          0.0015683340832809845 } } });
}

void test_costas(const std::vector<T> &samples, double w0) {
  auto filename = FLAGS_out_dir + "/costas.out";
  auto writer = SimpleFileWriter(filename, { "x", "y", "w", "ang", "res_ang",
                                             "ox", "oy", "fx", "fy", "err" });
  CostasLoop cl;
  w0 += 0;
  cl.init(0.0628, true);
  int n = samples.size();
  n = min(n, int(1e5));

  ReaderC8 reader("/tmp/data.out");
  vector<Type> check_tb;
  Type tmpv;
  while (reader.read(tmpv)) {
    check_tb.pb(tmpv);
  }

  REP (i, n) {
    complex<double> cur(samples[i].real(), samples[i].imag());
    // complex<double> tmp = FloatUtil::iexp(-0.3 * i);
    // cur = cur * tmp;
    complex<double> res = cl.push(cur);
    // res = cur;
    complex<double> fv = FloatUtil::iexp(cl.cur_ang);

    OPA_DISP0(cur, res, check_tb[i]);
    if (i >= check_tb.size())
      break;

    if (0) {
      writer->push(res.real(), res.imag(), cl.cur_w, std::arg(res),
                   std::arg(cur), cur.real(), cur.imag(), fv.real(), fv.imag(),
                   cl.err);
    }
  }
}

vector<T> simulated_carrier(double w0, int n = -1) {
  double ww0 = 0.001;
  double aww0 = 0.001;
  double noise_0 = 5;
  double a_w = 70;
  vector<T> res;

  double ang = 0;
  if (n == -1)
    n = 1e5;
  REP (i, n) {
    T cur;
    double nx = 1 * FloatUtil::get_gaussian<double>(0, 5);
    double ny = 1 * FloatUtil::get_gaussian<double>(0, 5);
    cur = T(a_w * cos(ang) + nx, a_w * sin(ang) + ny);
    ang += w0 + aww0 * cos(ww0 * i);
    res.pb(cur);
  }
  return res;
}

void test_filter(Filter<Type> *f) {
  FFT2<Type> fft;
  RealField<Type> r;
  int npw = 12;
  int n = 1 << npw;
  fft.init(&r, npw,
           FloatUtil::iexp<double>((FloatUtil::get_pi<double>() * 2 / n)));
  vector<Type> tb(n, 1);

  auto res = fft.ifft(tb);
  auto filename = FLAGS_out_dir + "/filter.out";

  {
    auto tmp = simulated_carrier(0.001, n);
    REP (i, n)
      res[i] = Type(tmp[i].real(), tmp[i].imag());
  }

  vector<Type> tsf(n);
  REP (i, n)
    tsf[i] = f->push(res[i]);
  // tsf = fft.fft(tsf);

  auto writer = SimpleFileWriter(filename, { "w", "x", "y", "ox", "oy" });
  REP (i, n) {
    writer->push(std::abs(tsf[i]), tsf[i].real(), tsf[i].imag(), res[i].real(),
                 res[i].imag());
  }
}

template <class T> class DataStream {
public:
  virtual T get() { OPA_CHECK0(false); }
  virtual void push(const T &a) { OPA_CHECK0(false); }
  virtual bool has_more() const { OPA_CHECK0(false); }
};

template <class T> class FileSourceDataStream : public DataStream<T> {
public:
  FileSourceDataStream(StringRef filename) {
    m_reader.init(filename);
    next = get();
  }

  virtual bool has_more() const { return !m_reader.is_eof(); }
  virtual T get() {
    T ret = next;
    m_reader.read<T>(next);
    return ret;
  }
  virtual void push(const T &a) { OPA_CHECK0(false); }

private:
  T next;
  ReaderC8 m_reader;
};

template <class T> class FileSinkDataStream : public DataStream<T> {
public:
  FileSinkDataStream(StringRef filename) { m_writer.init(filename); }

  virtual T get() { OPA_CHECK0(false); }
  virtual void push(const T &a) { m_writer.write(a); }

private:
  Writer m_writer;
};

template <class T> class DequeDataStream : public DataStream<T> {
public:
  virtual bool has_more() const { return m_q.size() > 0; }
  virtual T get() {
    T res = m_q.front();
    m_q.pop_front();
    return res;
  }
  virtual void push(const T &a) { m_q.push_back(a); }

private:
  std::deque<T> m_q;
};

template <typename A, typename B>
void conv(DataStream<complex<A> > &in, DataStream<complex<B> > &out) {
  int i = 0;
  while (in.has_more()) {
    OPA_DISP0(i++);
    complex<A> tmp = in.get();
    out.push(complex<B>(tmp.real(), tmp.imag()) / B(128.));
  }
}

template <class T>
void filter_data(DataStream<complex<T> > &in, DataStream<complex<T> > &out) {
  GrcCostatsLoop loop;
  loop.init(0.0628);
  int i = 0;
  while (in.has_more()) {
    auto e = in.get();
    complex<T> tmp(e.real(), e.imag());
    auto res = loop.push(tmp);
    OPA_DISP0(res);
    out.push(complex<T>(res.real(), res.imag()));
  }
}

void test_tracking() {
  init_data();
  OPA_CHECK0(!FLAGS_filename.empty());
  OPA_CHECK0(!FLAGS_out_dir.empty());
  auto sats = Satellites::get();
  ReaderC8 reader(FLAGS_filename);

  if (0) {
    REP (i, sats->nsats()) {
      auto e = sats->get(i);
      Writer writer(stdsprintf("/tmp/sat_%d.bin", i));
      REP (i, PRN_PERIOD) { writer.write<u8>(e.get_l1(i)); }
    }
    return;
  }

  if (1) {
    FileSourceDataStream<complex<s8> > fin(FLAGS_filename);
    DequeDataStream<complex<float> > fin_2;
    FileSinkDataStream<complex<float> > fout(FLAGS_out_dir + "/filtered.out");

    conv(fin, fin_2);
    filter_data(fin_2, fout);
    return;
  }

  int pn_period = 1023;
  double chip_rate = pn_period * 1000;
  double bit_spacing = samp_f / chip_rate;
  int n = pn_period * 40 * bit_spacing;

  vector<T> tb(n);
  REP (i, n)
    reader.read1(tb[i]);

  if (FLAGS_filter) {
    if (1) {
      test_filter(&iir_filter1);
    } else {
      test_filter(&lp_filter1);
    }
    return;
  }

  if (FLAGS_costas) {
    double w0 = 0.003;
    // tb = simulated_carrier(w0);
    test_costas(tb, w0);
    return;
  }

  int sz_pn_off = pn_period * 2;
  int sz_carrier_off = 100;
  // sz_pn_off = 100;

  OPA_DISP0(n, sz_pn_off, sz_carrier_off);
  double carrier_off = 420000 * 0;

  double carrier_off_lb = -1000, carrier_off_ub = 1000;
  double pn_off_lb = 0, pn_off_ub = pn_period * bit_spacing;
  double carrier_delta_off = 500;
  sz_carrier_off = 100;
  carrier_off_lb = carrier_off - sz_carrier_off * carrier_delta_off / 2;
  carrier_off_ub = carrier_off_lb + carrier_delta_off * sz_carrier_off;

  // number of samples allowed to have carrier rotation of 10%
  int correl_nsamp =
    0.1 / (carrier_delta_off / 2 / samp_f) / pn_period / bit_spacing;
  OPA_DISP("will drift by ",
           carrier_delta_off / 2 / samp_f * pn_period * bit_spacing);
  correl_nsamp = max(correl_nsamp, 1);
  OPA_DISP0(correl_nsamp);

  auto range_carrier_offset =
    Range<double>::Build_range(carrier_off_lb, carrier_off_ub, sz_carrier_off);
  auto range_pn_offset =
    Range<double>::Build_range(pn_off_lb, pn_off_ub, sz_pn_off);

  // sat_id, carrier_off, pn_pos, correl
  typedef std::tuple<int, double, double, double> OutType;
  std::ofstream of(FLAGS_out_dir + "/tracking_res.csv", std::ofstream::binary);
  CsvTupleRecordWriter<std::ofstream, int, double, double, double>
    record_writer({ "sat_id", "carrier_off", "pn_off", "correl" });
  CsvWriter<OutType, std::ofstream> writer(record_writer.get_dummy_sptr(), &of);

  REP (sat_id, sats->nsats()) {
    if (sat_id < 26 || sat_id > 26)
      continue;
    auto sat = sats->get(sat_id);
    for (auto &v1 : range_carrier_offset.tb()) {
      double w = FloatUtil::get_pi<double>() * 2 * v1 / samp_f;
      for (auto &v2 : range_pn_offset.tb()) {

        if (0) {
          double pos = v2;
          double metric = 0;
          int cnt = 0;
          REP (kk, correl_nsamp) {
            complex<double> cur = 0;
            REP (j, pn_period) {
              int id1 = round(pos);
              OPA_CHECK0(id1 < n);
              auto off = FloatUtil::iexp<double>(w * pos);
              off = off * complex<double>(tb[id1].real(), tb[id1].imag());
              int v = 2 * sat.get_l1(j) - 1;
              cur += complex<double>(v * off.real(), v * off.imag());
              pos += bit_spacing;
            }
            metric += abs(cur);
            ++cnt;
          }
          metric /= cnt;
          writer.add(OutType(sat_id, v1, v2, metric));
        } else {

          complex<double> correl = 0;
          int nx = correl_nsamp * pn_period * bit_spacing;
          REP (i, nx) {
            double pos = v2 + i;
            int sgn = 2 * sat.get_l1(i / pn_period / bit_spacing) - 1;
            complex<double> cur(tb[pos].real(), tb[pos].imag());
            complex<double> expected =
              FloatUtil::iexp<double>(w * pos) * complex<double>(sgn, 0);
            correl += (cur * expected);
          }
          correl /= nx;
          writer.add(OutType(sat_id, v1, v2, abs(correl)));
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  opa::init::opa_init(argc, argv);

  test_tracking();
  return 0;
}
