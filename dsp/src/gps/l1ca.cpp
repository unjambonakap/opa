#include <opa/dsp/gps/l1ca.h>
#include <opa/utils/csv.h>
#include <opa/dsp/resources.h>
#include <opa/crypto/lfsr.h>
#include <opa/math/common/Types.h>

using namespace std;
using namespace opa::utils;
using namespace opa::crypto;
using namespace opa::math::common;

OPA_NAMESPACE(opa, dsp, gps)

Satellites::Satellites() {
  vector<tuple<int, int, int, std::string> > sat_data;
  {
    istringstream iss(l1ca_dat);
    CsvReader<istringstream> reader(&iss);
    CsvFieldReader<int, int, int, std::string> field_reader(&reader);
    sat_data = field_reader.get_all();
  }

  vector<u32> seq1;
  vector<u32> seq2;
  {
    LFSR<u32> l1, l2;
    Poly<u32> state = PR_GF2.import({ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
    Poly<u32> p1 = PR_GF2.import({ 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, false);
    Poly<u32> p2 = PR_GF2.import({ 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1 }, false);
    l1.init(&GF2, state, p1);
    l2.init(&GF2, state, p2);
    REP (i, PRN_PERIOD) {
      seq1.pb(l1.get_next_non_galois_gps_style());
      seq2.pb(l2.get_next_non_galois_gps_style());
    }
  }

  REP (i, sat_data.size()) {
    Satellite sat;
    vector<u32> seqi;
    int delay = std::get<1>(sat_data[i]);
    int need = bignum(std::get<3>(sat_data[i]), 8).getu32();

    delay = PRN_PERIOD - delay;
    REP (t, PRN_PERIOD) { seqi.pb(seq1[t] ^ seq2[(delay + t) % PRN_PERIOD]); }

    if (0) {
      REP (j, 10)
        printf("%d", seqi[j]);
      OPA_DISP(" ", std::get<3>(sat_data[i]), need);
    }

    REP (j, 10)
      OPA_CHECK0((need >> j & 1) == seqi[9 - j]);

    sat.init({ seqi });
    m_sats.pb(sat);
  }
}

void Satellite::init(const Params &params) {
  opa::utils::Initable::init();
  m_params = params;
}

int Satellite::get_l1(int x) const {
  x = (x % PRN_PERIOD + PRN_PERIOD) % PRN_PERIOD;
  return m_params.l1ca_seq[x];
}

OPA_NAMESPACE_END(opa, dsp, gps)
