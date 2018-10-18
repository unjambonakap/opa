#include "jit.h"
#include <opa/crypto/la/algo.h>
#include <opa/crypto/la/base.h>
#include <opa/crypto/la/context.h>
#include <opa/crypto/la/relation_driver.h>
#include <opa/math/common/Types.h>
#include <opa_common.h>

using namespace std;

OPA_NAMESPACE(opa, crypto, la)

void ExactRelations::init(int input_size, int output_size,
                          const vector<Relation> &rels) {
  std::vector<int> filtered;
  REP (i, rels.size()) {
    if (!rels[i].cost.is_exact())
      continue;
    filtered.pb(i);
  }

  m_mat.initialize(&opa::math::common::GF2, filtered.size(),
                   output_size + input_size);
  this->input_size = input_size;
  this->output_size = output_size;

  REP (i, filtered.size()) {
    auto &r = rels[filtered[i]];
    RangeCoverage row;
    row.add(r.out).add(r.in.shift(output_size));
    m_mat.setRow(i, row.to_binary_vec(m_mat.getM()));
  }
  m_mat.row_echelon(-1, output_size);
}

bool ExactRelations::try_rel(const Basis &out, Basis &res) const {
  out.check_init();
  opa::math::common::Matrix<u32> m2;
  m2.initialize(out.mat().get_ring(), out.dim(), input_size + output_size);
  m2.set_submatrix(out.mat().get_submatrix(0, 0, out.dim()));

  std::vector<std::vector<u32> > tb;
  if (!m_mat.reduce(m2, output_size))
    return false;
  res.mat().affect(m2.get_mutable(0, output_size));
  res.reduce();
  return true;
}

CipherBlock::Params::Params() {
  m_rel_driver = SPTR(SimpleRelDriver)(new SimpleRelDriver(0.5, 0.5e3));
}

CipherBlock::Params &CipherBlock::Params::rel_driver(RelationDriver *driver) {
  m_rel_driver.reset(driver);
  return *this;
}

void CipherBlock::get_relations_dumb(Relations &rels) const {
  vector<u32> evals(1 << input_size());
  do_evals(&evals[0]);
  FOR (i, 1, 1 << input_size())
    FOR (j, 1, 1 << output_size()) {
      int cnt = 0;
      REP (k, 1 << input_size())
        cnt += opa::utils::dot(i, k) ^ opa::utils::dot((u32)j, evals[k]);
      Relation::RelationCost cost;
      cost.from_count(cnt, 1 << input_size());
      OPA_DISP0("REL >> ", i, j, cnt, 1<<input_size(), cost.bias());
      if (!driver()->should_add_rel(cost))
        continue;

      Relation r;
      r.in.from(i);
      r.out.from(j);
      r.c = cnt >= (1 << (input_size() - 1));
      r.cost = cost;
      driver()->add_rel(r);
    }
  driver()->set_rels(rels);
}

void CipherBlock::get_relations_walsh_diff(Relations &rels) const {

  int n = input_size();
  int m = output_size();
  vector<u32> evals(1 << n);
  do_evals(&evals[0]);

  opa::utils::ScopedHeapAlloc<u64> tb(1 << (n + m));
  memset(tb.buf(), 0, tb.bytesize());
  REP (i, 1 << n) { ++tb[i | (evals[i] << n)]; }
  do_walsh<u64>(tb.buf(), n + m);

  REP (i, tb.size())
    tb[i] = tb[i] * tb[i];
  do_iwalsh(tb.buf(), n + m);

  REP (i, 1 << (n + m)) {
    int in = i & (1 << n) - 1;
    int out = i >> n;
    if (!in || !out)
      continue;
    Relation::RelationCost cost;
    int v = tb[i] >> (n + m);
    cost.from_walsh_count(v, 1 << n);

    if (0) {
      int cnt_check = 0;
      REP (i, 1 << n) { cnt_check += (evals[i] ^ evals[i ^ in]) == out; }
      OPA_DISP0(in, out, v, cnt_check, 1 << n);
    }

    if (!driver()->should_add_rel(cost))
      continue;
    Relation rel;
    rel.in = RangeCoverage::From(in, n);
    rel.out = RangeCoverage::From(out, m);
    rel.c = 0;
    rel.cost = cost;
    driver()->add_rel(rel);
  }
  driver()->set_rels(rels);
}

void CipherBlock::get_relations_walsh(Relations &rels) const {
  OPA_CHECK0(input_size() <= 16);
  opa::utils::ScopedHeapAlloc<int> tb(1 << input_size());

  vector<u32> evals(1 << input_size());
  do_evals(&evals[0]);

  FOR (om, 1, 1 << output_size()) {
    REP (i, 1 << input_size()) {
      int x = opa::utils::dot((u32)evals[i], (u32)om);
      tb[i] = 2 * x - 1;
    }
    do_walsh(tb.buf(), input_size());

    FOR (i, 1, 1 << input_size()) {
      Relation::RelationCost cost;
      cost.from_walsh_count(tb[i], 1 << input_size());
      OPA_DISP0(cost.bias());
      if (!driver()->should_add_rel(cost))
        continue;
      Relation rel;
      rel.in = RangeCoverage::From(i, input_size());
      rel.out = RangeCoverage::From(om, output_size());
      rel.c = tb[i] >= 0;
      rel.cost = cost;
      driver()->add_rel(rel);
    }
  }
  driver()->set_rels(rels);
}

void CipherBlock::do_get_pre(const Basis &out, Basis &res) const {
  const auto &exact_rels = get_exact_relations();
  res.init(input_size());
  if (exact_rels.try_rel(out, res)) {
  } else {
    do_get_pre_fail(out, res);
  }
  res.reduce();
}

std::vector<RangeCoverage>
CipherBlock::get_pre(const RangeCoverage &out) const {
  OPA_CHECK0(out.size() > 0);
  Basis o, tmp;
  o.init(output_size());
  o.add(out);
  OPA_CHECK0(o.dim() > 0);

  do_get_pre(o, tmp);
  std::vector<RangeCoverage> res;
  REP (i, tmp.dim())
    res.pb(RangeCoverage::From_binary_vec(tmp.mat().get_row(i).tovec()));
  return res;
}

void CipherBlock::fini() {
  if (m_func)
    context()->runtime().release((void *)m_func);
  Initable::fini();
}

BitVec CipherBlock::evaluate_jit(const BitVec &iv, const BitVec &kv) const {

  std::unique_ptr<u8> key_store;
  int ni = raw_input_size() / call().input_base_size;
  FunctionCaller caller;
  caller.func_ptr = get_jit_func();

  REP (i, ni) {
    caller.args.pb(
      iv.to<uintptr_t>(i * call().input_base_size, call().input_base_size));
  }

  if (key_size() > 0) {
    int nk = key_size() / call().key_base_size;
    if (nk > 1) {
      key_store.reset(new u8[key_size() / 8]);
      OPA_DISP("keystore > >", key_size(), uintptr_t(key_store.get()));
      int kb = call().key_base_size;
      REP (i, nk) {
        if (kb == 8) {
          *((u8 *)key_store.get() + i) = kv.to<u8>(i * kb, kb);
        } else if (kb == 16) {
          *((u16 *)key_store.get() + i) = kv.to<u16>(i * kb, kb);
        } else {
          OPA_CHECK0(0);
        }
      }
      caller.args.pb((uintptr_t)key_store.get());
    } else {
      caller.args.pb(kv.to<uintptr_t>());
    }
  }
  u64 res = (u64)caller.call();

  return BitVec::From(res, output_size());
}

void CipherBlock::create_jit_builder(JitBuilder &builder,
                                     JitContext &jc) const {
  OPA_CHECK0(jc.c != 0);
  JitBuilder res(jc.c);
  jc.vars.emplace_back(
    new asmjit::X86GpVar(*jc.c, context()->jit().sizeToType(output_size())));
  res.output_var = jc.vars.back().get();

  int var_size = call().input_base_size;
  OPA_CHECK0(raw_input_size() % var_size == 0);

  int nvars = raw_input_size() / var_size;
  REP (i, nvars) {
    jc.vars.emplace_back(
      new asmjit::X86GpVar(*jc.c, context()->jit().sizeToType(var_size)));
    res.input_vars.pb(jc.vars.back().get());
  }

  res.key_var = nullptr;
  if (key_size() != 0) {
    OPA_DISP("CREATE JIT >> ", call().key_base_size, key_size());
    if (call().key_base_size < key_size()) {
      OPA_CHECK0(key_size() % call().key_base_size == 0);
      jc.vars.emplace_back(
        new asmjit::X86GpVar(*jc.c, asmjit::kVarTypeUIntPtr));
    } else {
      OPA_CHECK0(call().key_base_size == key_size());
      jc.vars.emplace_back(
        new asmjit::X86GpVar(*jc.c, context()->jit().sizeToType(key_size())));
    }
    res.key_var = jc.vars.back().get();
  }

  builder = res;
}

uintptr_t CipherBlock::get_jit_func() const {
  if (!m_func) {

    asmjit::FileLogger logger(stdout);
    int var_size = call().input_base_size;
    int nvars = input_size() / var_size;
    asmjit::X86Compiler c(&context()->runtime());
    c.setLogger(&logger);
    JitContext ctx;
    ctx.c = &c;
    JitBuilder builder;
    create_jit_builder(builder, ctx);

    asmjit::FuncBuilderX func_builder;

    OPA_DISP("has input vars >> ", builder.input_vars.size());
    for (auto i : builder.input_vars)
      func_builder.addArg(i->getVarType());
    if (builder.key_var)
      func_builder.addArg(builder.key_var->getVarType());
    func_builder.setRet(builder.output_var->getVarType());

    c.addFunc(asmjit::kFuncConvHost, func_builder);
    int pos = 0;
    for (auto i : builder.input_vars)
      c.setArg(pos++, *i);
    if (builder.key_var)
      c.setArg(pos++, *builder.key_var);

    setup_jit(builder);
    c.ret(*builder.output_var);
    OPA_DISP("FINALIZING ", desc());
    c.endFunc();

    m_func = (uintptr_t)c.make();
    OPA_DISP("func at >> ", m_func);
  }
  OPA_CHECK0(m_func != 0);
  return m_func;
}

BitVec CipherBlock::fast_eval_fe(const BitVec &input, const BitVec &kv) const {
  BitVec in(input_size());
  in.copy(input);
  in.copy(kv, raw_input_size());
  return fast_eval_fe(in);
}

BitVec CipherBlock::fast_eval_fe(const BitVec &iv) const {
  OPA_CHECK0(fast_eval());
  int s = call().n_input_args + call().n_key_args;
  u64 res;

  if (call().key_array) {
    u64 a = iv.to<u64>(0, raw_input_size());
    u64 tb[call().n_key_args];
    REP (i, call().n_key_args) {
      tb[i] = iv.to<u64>(raw_input_size() + call().key_base_size * i,
                         call().key_base_size);
    }
    res = fast_eval_tb(a, tb);
  } else if (s == 1) {
    res = fast_eval1(iv.to<u64>());
  } else if (s == 2) {
    res = fast_eval2(iv.to<u64>(0, call().input_base_size),
                     iv.to<u64>(call().input_base_size));
  } else
    OPA_CHECK(0, s, desc());
  return BitVec::From(res, output_size());
}

void CipherBlock::do_evals(u32 *tb) const {
  if (fast_eval()) {
    // OPA_DISP0(call().input_array, call().key_array, call().key_size);
    OPA_CHECK0(!call().input_array && !call().key_array);

    int s = call().n_input_args + call().n_key_args;
    // OPA_DISP0(call().n_input_args, call().n_key_args);
    if (s == 1) {
      REP (i, 1 << input_size()) { tb[i] = fast_eval1(i); }

    } else if (s == 2) {
      int s0 = call().input_base_size;
      int s1 = call().n_input_args == 1 ? call().key_base_size
                                        : call().input_base_size;
      int M0 = (1 << s0) - 1;

      REP (i, 1 << input_size()) { tb[i] = fast_eval2(i & M0, i >> s0); }
    } else
      OPA_CHECK0(0);

  } else {
    REP (i, 1 << input_size()) {
      tb[i] = evaluate(BitVec::From(i, input_size())).to<u32>();
    }
  }
}

void CipherBlock::test_pre(int nr, const RangeCoverage &r) {
  auto res = get_pre(r);
  BitVec kv = BitVec::rand(key_size());
  auto data = gen_rand(nr, kv);
  for (auto &e : res) {
    e.filter_greater(raw_input_size() - 1);
    OPA_DISP0(e);
  }

  OPA_CHECK0(res.size() <= 32);

  std::map<u32, int> seen;
  for (auto &x : data.x) {
    u32 a = 0;
    for (auto &e : res) {
      a = a << 1 | x.ST.dot(e);
    }
    u32 b = x.ND.dot(r);
    if (seen.count(a))
      OPA_CHECK_EQ(seen[a], b, a);
    seen[a] = b;
  }
  OPA_DISP("TEST PRE OK", seen.size(), res.size(), nr);
}

void CipherBlock::test_rel(int nr, const RangeCoverage &irel,
                           const RangeCoverage &orel, bool diff) const {
  if (1 || input_size() - key_size() > 16) {
    int cnt = 0;
    BitVec kv = BitVec::rand(key_size());
    auto data = gen_rand(nr, kv);
    if (diff) {
      for (auto &x : data.x) {
        BitVec tmp(input_size());
        tmp.copy(x.first);
        tmp.copy(kv, raw_input_size());

        BitVec tmp2 = tmp;
        tmp2.sxorz(irel);
        BitVec res2 = fast_eval_fe(tmp2);
        OPA_DISP0(tmp, tmp2, irel, x.second, res2, orel);

        res2.sxorz(x.second);
        res2.sxorz(orel);
        cnt += res2.isz();
      }

    } else {
      for (auto &x : data.x) {
        BitVec tmp(input_size());
        tmp.copy(x.first);
        tmp.copy(kv, raw_input_size());
        cnt += tmp.dot(irel) ^ x.second.dot(orel);
      }
    }
    OPA_DISP0(cnt, nr, 1. * cnt / nr);
  }
}

BitVec CipherBlock::evaluate(const BitVec &iv_) const {
  OPA_CHECK0(input_size() == iv_.size());
  BitVec iv = iv_.extract(0, raw_input_size());
  BitVec kv = iv_.extract(raw_input_size());

  BitVec res;
  res.init(output_size());

  do_evaluate(iv, kv, res);
  OPA_CHECK0(res.size() == output_size());
  return res;
}

double CipherBlock::test_outrel(const BitVec &kv, const OutRel &orel) const {
  OPA_CHECK0(kv.size() == key_size());
  if (orel.in_diff.size() > 0) {
    int cnt = 0;
    int pos = 0;
    for (auto &x : orel.in_diff) {
      ++pos;
      BitVec v = evaluate2(BitVec::From(x.b1, raw_input_size()), kv);
      BitVec v2 = evaluate2(BitVec::From(x.b2, raw_input_size()), kv);
      cnt += (v.to<u64>() ^ v2.to<u64>()) == x.diff_res;
    }
    double tmp = 1. * cnt / orel.in_diff.size();
    return tmp;

  } else {
    int cnt = 0;
    int pos = 0;
    for (auto &x : orel.in) {
      ++pos;
      BitVec v = evaluate2(BitVec::From(x.first, raw_input_size()), kv);
      BitVec tmp = BitVec::From(x.first, input_size());
      tmp.copy(kv, raw_input_size());

      cnt += v.dot(BitVec::From(orel.out, output_size())) ^ x.second;
    }

    // OPA_DISP("GOT ", cnt, orel.in.size());
    double tmp = 1. * cnt / orel.in.size();
    tmp = fabs(tmp - 0.5) * 2;
    return tmp;
    // OPA_DISP0(cnt, orel.in.size(), tmp, fabs(tmp - 0.5) * 2);
  }
}

CipherData CipherBlock::gen_rand(int nr, const BitVec &kv) const {
  CipherData data;
  REP (i, nr) {
    BitVec in = BitVec::rand(raw_input_size());
    BitVec out = fast_eval_fe(in, kv);
    data.x.pb(MP(in, out));
  }
  return data;
}

void CipherBlock::gen_rand_diff(const BitVec &diff, int nr, const BitVec &kv,
                                CipherData *out_cipher_data) const {
  REP (i, nr) {
    BitVec in1 = BitVec::rand(raw_input_size());
    BitVec in2 = in1.xorz(diff);
    BitVec out1 = fast_eval_fe(in1, kv);
    BitVec out2 = fast_eval_fe(in2, kv);
    out_cipher_data->x_diff.emplace_back(MP(in1, out1), MP(in2, out2));
  }
}

bool CipherBlock::verify_data(const CipherData &data, const BitVec &kv) const {
  for (const auto &e : data.x) {
    BitVec out = fast_eval_fe(e.ST, kv);
    if (out != e.second)
      return false;
  }

  for (const auto &e : data.x_diff) {
    BitVec o1 = fast_eval_fe(e.ST.ST, kv);
    if (o1 != e.ST.ND) {
      return false;
    }

    BitVec o2 = fast_eval_fe(e.ND.ST, kv);
    if (o2 != e.ND.ND) {
      return false;
    }
  }

  return true;
}

OPA_NAMESPACE_END(opa, crypto, la)
