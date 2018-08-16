#include <gsl/gsl_poly.h>
#include <opa/crypto/la/algo.h>
#include <opa/crypto/stream/target.h>
#include <opa/crypto/stream/plan.h>
#include <opa/crypto/stream/solver.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/math/common/stats.h>
#include <opa/utils/csv.h>
#include <opa/utils/misc.h>
#include <opa_common.h>
#include <opa/math/common/algebra/basis.h>
#include <opa/math/common/Matrix.h>
#include <opa/threading/runner.h>

using namespace std;
using namespace opa::crypto::stream;
using namespace opa::crypto::la;
using namespace opa::crypto;
using namespace opa::math::common;
using namespace opa::utils;
DEFINE_string(obs_file, "", "");
DEFINE_string(out_file, "", "");
DEFINE_string(coeffs_file, "", "");
DEFINE_string(check_key, "", "");

void do_solve(const SboxLfsr &sbox_lfsr, const ObservedData &obs_data,
              const BitVec &ans, bool check) {
  SboxLfsrSolver solver;
  SboxLfsrSolver::Params params = { sbox_lfsr, ans, check };
  solver.init(params);
  set<int> known_data;
  for (auto &x : obs_data) {
    known_data.insert(x.ST);
  }

  solver.load_plan_from_data(known_data);
  BitVec key = solver.solve(obs_data);
  OPA_DISP("Got key >> ", key.str());

  if (!FLAGS_out_file.empty()) {
    std::ofstream ofs(FLAGS_out_file);
    ofs << key.str();
  }
}

void load_sbox_lfsr(SboxLfsr &lfsr) {
  SboxDesc desc;
  int NI, NO;
  int nlfsr;

  vector<vi> lfsr_coeffs;
  {
    std::ifstream ifs(FLAGS_coeffs_file);
    ifs >> NI >> NO >> nlfsr;
    vector<u8> sbox_coeffs;
    REP (i, 1 << NI) {
      int a;
      ifs >> a;
      sbox_coeffs.pb(a);
    }
    desc.init(NO, sbox_coeffs);

    lfsr_coeffs.resize(nlfsr);
    REP (i, nlfsr) {
      int nx;
      ifs >> nx;
      auto &cur = lfsr_coeffs[i];
      cur.resize(nx);
      REP (j, nx)
        ifs >> cur[j];
    }
    OPA_DISP0(lfsr_coeffs);
  }

  SboxLfsr::Params params;
  params.sbox = desc;

  REP (i, nlfsr) {
    LFSR_GF2_small lfsr;
    u64 poly = 0;
    for (auto e : lfsr_coeffs[i])
      poly |= 1ull << e;

    lfsr.init(poly, lfsr_coeffs[i].back(), 0);
    params.lfsrs.emplace_back(lfsr);
  }
  lfsr.init(params);
}

void load_observed_data(ObservedData &observed_data) {

  std::ifstream ifs(FLAGS_obs_file);
  CsvReader<ifstream> reader(&ifs, false /* headers */);
  CsvFieldReader<int, int> field_reader(&reader);
  auto rels = field_reader.get_all();

  puts("read all");
  for (auto &e : rels) {
    observed_data[std::get<0>(e)] = std::get<1>(e);
  }
}

BitVec get_key_from_str(const string &tmp) {
  BitVec ans_input(tmp.size());
  REP (i, tmp.size())
    ans_input.set(i, tmp[i] == '1');
  return ans_input;
}

void test1() {

  constexpr int NI = 6;
  constexpr int NO = 4;
  double P0 = 0.70;
  int nweak = 2;
  vector<int> degs={10,11,9,12,13,14};

  int nlfsr = degs.size();
  SboxDesc desc;
  if (1) {
    desc.init(NO, {
                    7, 6, 5, 10, 8, 1, 12, 13, 6, 11, 15, 11, 1, 6, 2, 7, 0, 2,
                    8, 12, 3, 2, 15, 0, 1, 15, 9, 7, 13, 6, 7, 5, 9, 11, 3, 3,
                    12, 12, 5, 10, 14, 14, 1, 4, 13, 3, 5, 10, 4, 9, 11, 15, 10,
                    14, 8, 13, 14, 2, 4, 0, 0, 4, 9, 8,
                  });
  } else {
    int id = 0;
    while (true) {
      desc = SboxDesc::RandWeak(NI, NO, nweak, P0);
      rng.seed(id++);
      desc = SboxDesc::Rand(NI, NO);
      SboxBlock blk;
      blk.init(desc);
      auto rels = blk.get_relations().tb;
      reverse(ALL(rels));
      bool ok = rels[0].cost.prob() > 0.65;
      if (ok) {
        REP (i, 1 << NI) { printf("%d,", desc.get(i)); }
        puts("");

        puts("GOT ");
        OPA_DISP0(id - 1, rels);
        rng.seed(id - 1);
        OPA_DISP0(rng());
        OPA_DISP0(rng());
      }
    }
  }

  SboxLfsr lfsr;
  SboxLfsr::Params params;
  params.sbox = desc;

  REP (i, nlfsr) {
    auto lfsr = LFSR<u32>();
    lfsr.init_rand(degs[i], &GF2);
    lfsr.get_poly().disp();
    params.lfsrs.pb(LFSR_GF2_small::from_lfsr(lfsr));
  }

  lfsr.init(params);
  BitVec ans_input;
  ans_input = lfsr.get_state();
  OPA_DISP("GOT ANS >> ", ans_input.str());

  ObservedData observed_data;
  int end = 30000 + 10;

  REP (i, end) {
    int cur = lfsr.get();
    observed_data[i] = cur;
  }

  {
    vi filter;
    REP (i, end / NO)
      filter.pb(i);
    random_shuffle(ALL(filter));

    ObservedData obs2 = observed_data;
    observed_data.clear();
    int nkeep = end * 9 / 10;
    // keep ~ nkeep relations
    filter.resize(nkeep / NO);
    for (auto i : filter) {
      REP (j, NO)
        observed_data[i * NO + j] = obs2[i * NO + j];
    }
  }

  do_solve(lfsr, observed_data, ans_input, true);
}

void solve_real() {
  OPA_CHECK0(FLAGS_obs_file.size() > 0);
  SboxLfsr lfsr;
  load_sbox_lfsr(lfsr);

  ObservedData observed_data;
  load_observed_data(observed_data);

  auto ans = get_key_from_str(FLAGS_check_key);

  do_solve(lfsr, observed_data, ans, !FLAGS_check_key.empty());
}

int main(int argc, char **argv) {
  opa::init::opa_init(argc, argv);
  if (1) {
    opa::threading::Runner::EasySetup();
  }

  if (1) {
    test1();
  } else {
    solve_real();
  }

  return 0;
}
