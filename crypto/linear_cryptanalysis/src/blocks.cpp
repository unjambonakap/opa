#include <opa_common.h>

#include "jit.h"
#include <opa/crypto/la/blocks.h>

using namespace opa::math::common;
using namespace std;
using namespace opa::utils;

OPA_NAMESPACE(opa, crypto, la)

SboxDesc SboxDesc::Rand(int in_bitsize, int out_bitsize) {
  SboxDesc res;
  int ni = 1 << in_bitsize;
  int no = 1 << out_bitsize;
  std::vector<u8> tmp(ni);
  REP (i, ni)
    tmp[i] = i % no;
  random_shuffle(ALL(tmp));
  // REP (i, ni)
  //  tmp[i] = opa::math::common::rng() % no;
  res.init(out_bitsize, tmp);
  return res;
}

int update_rel(vector<u8> &tb, set<u8> &can, int rel, int ni, int ob,
               double prob, bool dry_run) {
  vector<u8> can2;
  int good = 0;
  REP (k, ni) {
    int v = dot(k, rel);
    v ^= (tb[k] >> ob) & 1;
    good += v == 0;
    if (v == 1 && can.count(k))
      can2.pb(k);
  }
  int v1 = good + can2.size();
  int v2 = ni - good + (can.size() - can2.size());
  if (dry_run) {
    return max(v1, v2);
  }

  double need = ni * prob;

  if (v1 < v2 || (0 && v1 == v2 && good > ni - good)) {
    for (auto x : can2)
      can.erase(x);
    good = ni - good;
  } else {
    can = set<u8>(ALL(can2));
  }
  vector<u8> tmp(ALL(can));

  while (need > good && tmp.size()) {
    int e = rng() % tmp.size();
    int x = tmp[e];
    tmp[e] = tmp.back();
    tmp.pop_back();
    can.erase(x);
    ++good;
    tb[x] ^= 1 << ob;
  }

  return 0;
}

SboxDesc SboxDesc::RandWeak(int in_bitsize, int out_bitsize, int nweak,
                            double mx_v) {
  OPA_CHECK0(nweak <= out_bitsize);
  int ni = 1 << in_bitsize;
  int no = 1 << out_bitsize;
  SboxDesc res = SboxDesc::Rand(in_bitsize, out_bitsize);

  REP (ob, nweak) {
    vector<double> targets;
    set<int> rels;
    // blacklist
    rels.insert(0);
    double prob = mx_v;
    set<u8> can;
    REP (k, ni)
      can.insert(k);

    REP (j, 1) {

      int best = -1, bestv = -1;
      int ntry = j == 0 ? 1 : 10;
      ntry = 1;
      REP (tryid, ntry) {
        int rel;
        while (true) {
          rel = rng() % (1 << in_bitsize);
          if (rels.count(rel))
            continue;
          break;
        }
        int nv = update_rel(res.mp, can, rel, ni, ob, prob, true);
        if (nv > bestv)
          bestv = nv, best = rel;
      }
      update_rel(res.mp, can, best, ni, ob, prob, false);
      int cnt = 0;
      REP (k, ni) {
        int x = dot(best, k);
        x ^= (res.mp[k] >> ob & 1);
        cnt += x;
      }
      // OPA_DISP0(cnt, 1. * cnt / ni, prob, best);
      rels.insert(best);
      targets.pb(prob);
      prob *= 0.90;
    }
  }

  if (0)
    while (true) {
      Matrix<u32> mix_mat = Matrix<u32>::rand(&GF2, out_bitsize, out_bitsize);
      if (mix_mat.get_det() != 0) {
        REP (k, ni) {
          Matrix<u32> v(&GF2, out_bitsize, 1);
          REP (ob, out_bitsize)
            v(ob) = res.get(k) >> ob & 1;
          v = mix_mat * v;
          res.mp[k] = 0;
          REP (ob, out_bitsize)
            res.mp[k] |= v(ob) << ob;
        }
        break;
      }
    }
  return res;
}

void ExpansionBlock::do_get_relations_diff(Relations &rels) const {
  vector<vi> from(input_size());
  REP (i, output_size()) { from[m_ed.from(i)].pb(i); }
  REP (i, input_size()) {
    BitVec v(input_size() + output_size());
    v.toggle(i);
    for (auto &x : from[i])
      v.toggle(input_size() + x);
    rels.add_basis_vec(v);
  }
}

void ExpansionBlock::do_get_relations(Relations &rels) const {
  REP (i, m_ed.size()) {
    Relation r;

    r.in = { m_ed.from(i) };
    r.out = { i };
    r.cost.set_bias(1);
    rels.add_rel(r);
  }
}

u64 PermutationDesc::get2(u64 v) const {

  if (size() <= 16)
    return ptr()[v];
  u64 res = 0;
  REP (i, size())
    res |= u64(v >> i & 1) << get(i);
  return res;
}

void AddBlock::do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const {
  u64 a = iv.to<u64>(0, input_size() / 2);
  u64 b = iv.to<u64>(input_size() / 2, input_size() / 2);
  ov.from(a + b);
}

void AddBlock::do_get_relations_diff(Relations &rels) const {
  get_relations_walsh_diff(rels);
}

void AddBlock::do_get_relations(Relations &rels) const {
  get_relations_walsh(rels); // i am not a smart man
}

void AddBlock::do_get_pre_fail(const Basis &out, Basis &res) const {
  int mx = -1;
  REP (i, out.dim()) {
    REP (j, out.sz())
      if (out.mat().get(i, j))
        mx = std::max(mx, j);
  }
  REP (i, mx + 1) {
    res.add(RangeCoverage().add(i));
    res.add(RangeCoverage().add(sz + i));
  }
}

void AddBlock::setup_jit(const JitBuilder &builder) const {
  builder.c->mov(*builder.output_var, *builder.input_vars[0]);
  builder.c->add(*builder.output_var, *builder.input_vars[1]);
}

void AddConstBlock::do_evaluate(const BitVec &iv, const BitVec &kv,
                                BitVec &ov) const {
  u64 a = iv.to<u64>(0, input_size());
  ov.from(a + m_v);
}

void AddConstBlock::do_get_relations_diff(Relations &rels) const {
  get_relations_walsh_diff(rels);
}
void AddConstBlock::do_get_relations(Relations &rels) const {
  get_relations_walsh(rels); // i am not a smart man
}

void AddConstBlock::do_get_pre_fail(const Basis &out, Basis &res) const {
  int mx = -1;
  REP (i, out.dim()) {
    REP (j, out.sz())
      if (out.mat().get(i, j))
        mx = std::max(mx, j);
  }
  REP (i, mx + 1) { res.add(RangeCoverage().add(i)); }
}

void AddConstBlock::setup_jit(const JitBuilder &builder) const {
  builder.c->mov(*builder.output_var, *builder.input_vars[0]);
  builder.c->add(*builder.output_var, asmjit::Imm(m_v));
}

void XorBlock::do_get_relations_diff(Relations &rels) const {
  REP (i, m_blk_size) {
    REP (j, 3) {
      BitVec v(input_size()+output_size());
      if (j != 0)
        v.toggle(i);
      if (j != 1)
        v.toggle(m_blk_size+i);
      if (j != 2)
        v.toggle(input_size()+i);
      rels.add_basis_vec(v);
    }
  }
}

void XorBlock::do_get_relations(Relations &rels) const {
  REP (i, m_blk_size) {
    Relation r;
    r.in.add(i);
    r.in.add(m_blk_size + i);
    r.out.add(i);
    r.cost.set_bias(1);
    rels.add_rel(r);
  }
}

void XorBlock::setup_jit(const JitBuilder &builder) const {
  builder.c->mov(*builder.output_var, *builder.input_vars[0]);
  builder.c->xor_(*builder.output_var, *builder.input_vars[1]);
}

void PermutationBlock::do_get_relations_diff(Relations &rels) const {
  REP (i, input_size()) {
    BitVec v(input_size()+output_size());
    v.toggle(i);
    v.toggle(input_size()+m_perm.get(i));
    rels.add_basis_vec(v);
  }
  return do_get_relations(rels);
}

void PermutationBlock::do_get_relations(Relations &rels) const {
  REP (i, m_perm.size()) {
    Relation r;
    r.in.add(i);
    r.out.add(m_perm.get(i));
    r.cost.set_bias(1);
    rels.add_rel(r);
  }
}

void PermutationBlock::setup_jit(const JitBuilder &builder) const {
  /*
  asmjit::X86GpVar tmp(*builder.c, asmjit::kVarTypeUIntPtr);
  asmjit::X86GpVar tmp2(*builder.c, asmjit::kVarTypeUIntPtr);
  builder.c->mov(tmp, asmjit::Imm((u64)m_perm.ptr()));
  builder.c->movzx(tmp2, *builder.input_vars[0]);
  int shift = 1;
  builder.c->mov(*builder.output_var, asmjit::x86::ptr(tmp, tmp2, shift));
  */
}

void SboxBlock::setup_jit(const JitBuilder &builder) const {
  /*
  asmjit::X86GpVar tmp(*builder.c, asmjit::kVarTypeUIntPtr);
  asmjit::X86GpVar tmp2(*builder.c, asmjit::kVarTypeUIntPtr);
  builder.c->mov(tmp, asmjit::Imm((u64)sbox.ptr()));
  builder.c->movzx(tmp2, *builder.input_vars[0]);
  builder.c->mov(*builder.output_var, asmjit::x86::ptr(tmp, tmp2));
  */
}

void SboxBlock::do_get_relations_diff(Relations &rels) const {
  get_relations_walsh_diff(rels);
}

void SboxBlock::do_get_relations(Relations &rels) const {
  get_relations_walsh(rels);
}

OPA_NAMESPACE_END(opa, crypto, la)
