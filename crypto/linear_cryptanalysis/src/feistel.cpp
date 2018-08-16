#include <opa/crypto/la/feistel.h>
#include <opa/crypto/la/blocks.h>
#include <opa/math/common/Utils.h>

using namespace opa::OR;

#include <opa/threading/mapper.h>

// OPA_NAMESPACE(opa, threading)
// OPA_CLOUDY_JOB_IMPL_TMPL(MapJob, opa::utils::DoubleRes, int,
//                         SPTR(const opa::crypto::la::FeistelBlock),
//                         opa::crypto::la::RelState)
// OPA_NAMESPACE_END(opa, threading)
//

OPA_NAMESPACE(opa, crypto, la)

// template <>
// opa::threading::JobId opa::threading::MapJob<
//  opa::utils::DoubleRes, int, SPTR(const opa::crypto::la::FeistelBlock),
//  opa::crypto::la::RelState>::StaticJobId = -1;
// template <>
// std::string opa::threading::MapJob<opa::utils::DoubleRes, int,
//                                   SPTR(const opa::crypto::la::FeistelBlock),
//                                   opa::crypto::la::RelState>::JobName = "";
//
// OPA_CLOUDY_DECL_MAPPER(
//  test_r2_mapper,
//  std::function<opa::utils::DoubleRes(int, SPTR(const FeistelBlock),
//  RelState)>(
//    [](int nr, SPTR(const FeistelBlock) blk, RelState r) {
//      return opa::utils::DoubleRes(blk->test_r2_internal(nr, r));
//    }),
//  opa::utils::DoubleRes, int, SPTR(const FeistelBlock), RelState);

// OPA_REG_BLOCK(FeistelBlock);

void FeistelBlock::u64_to_lr(u64 lr, u32 &l, u32 &r) const {
  l = lr & (1ull << blk_size()) - 1;
  r = lr >> blk_size();
}

u64 FeistelBlock::lr_to_u64(u32 l, u32 r) const {
  return u64(r) << blk_size() | l;
}

u64 FeistelBlock::fast_eval_tb(u64 a, const u64 *b) const {
  u32 l, r;
  u64_to_lr(a, l, r);

  if (mid_xor()) {
    l ^= b[nround()];
    r ^= b[nround() + 1];
  }
  // std::swap(l, r);

  // OPA_DISP0(l, r, b[nround], b[nround+1]);
  REP (i, nround()) {
    u32 tmp = round()->fast_eval2(r, b[i]);
    u32 nr = l ^ tmp;
    l = r;
    r = nr;
    // OPA_DISP0(i, l, r, tmp, b[i]);
  }
  std::swap(l, r);
  return lr_to_u64(l, r);
}

void FeistelBlock::init(const Params &params) {
  m_params = params;
  params.round->check_init();
  // OPA_DISP("INIT FEISTELBLock ", params.nround,
  // uintptr_t(params.round.get()),
  //         params.mid_xor);
  // OPA_DISP("Feistel round func: ", round()->raw_input_size(),
  //         round()->input_size(), round()->output_size());

  m_blk_size = round()->raw_input_size();
  OPA_DISP0(m_blk_size, round()->raw_input_size(), round()->key_size(),
            round()->output_size());
  OPA_CHECK_EQ0(round()->output_size(), m_blk_size);

  int nkey = nround() * round()->key_size();
  if (mid_xor())
    nkey += 2 * m_blk_size;

  m_params.call(CallInfo(m_blk_size * 2, round()->key_size(), true, true))
    .input_size(m_blk_size * 2)
    .output_size(m_blk_size * 2)
    .key_size((nkey))
    .fast_eval(true);
  CipherGraph::init(m_params);
}

void FeistelBlock::init_graph() {
  round()->check_init();

  xorblock = context()->instanciate<XorBlock>();
  xorblock->init(m_blk_size);
  idb = context()->instanciate<IdBlock>();
  idb->init(m_blk_size);

  CipherNode *l = add_node(idb.get());
  CipherNode *r = add_node(idb.get());

  smart_plug({ PlugDesc(input, 0, m_blk_size) }, l);
  smart_plug({ PlugDesc(input, m_blk_size, m_blk_size) }, r);
  if (mid_xor()) {
    CipherNode *xl = add_node(xorblock.get());
    CipherNode *xr = add_node(xorblock.get());
    int basepos = nround() * round()->key_size();
    smart_plug({ l, PlugDesc(key, basepos, m_blk_size) }, xl);
    smart_plug({ r, PlugDesc(key, basepos + m_blk_size, m_blk_size) }, xr);
    l = xl;
    r = xr;
  }

  REP (i, nround()) {

    CipherNode *nr = add_node(xorblock.get());
    CipherNode *tmp_round = add_node(round().get());

    smart_plug(
      { r, PlugDesc(key, i * round()->key_size(), round()->key_size()) },
      tmp_round);
    smart_plug({ l, tmp_round }, nr);

    l = r;
    r = nr;
  }

  std::swap(l, r); // undo swap

  add_range_edges(l, output, 0, 0, m_blk_size);
  add_range_edges(r, output, 0, m_blk_size, m_blk_size);

  desc() =
    opa::utils::SPrintf("FeistelBlock(%d, %d), graph_desc=%s", output_size(),
                        nround(), build_graph_desc().c_str());
}

void FeistelBlock::do_get_relations(Relations &rels) const {
  const Relations &brel = round()->get_relations();
  OPA_CHECK0(0);
}

BitVec FeistelBlock::get_round_key(const BitVec &kv, int r) const {
  return kv.extract(r * round()->key_size(), round()->key_size());
}

void FeistelBlock::undo_last_round_one(const CipherData::IOPair &pair,
                                       const BitVec &kv,
                                       CipherData::IOPair *output_pair) const {
  u32 l, r;
  this->u64_to_lr(pair.second.to<u64>(), l, r);
  l ^= round()->fast_eval2(r, kv.to<u64>());
  std::swap(l, r);
  *output_pair = MP(pair.first, BitVec::From(lr_to_u64(l, r), output_size()));
}

void FeistelBlock::undo_last_round(const CipherData &data, const BitVec &kv,
                                   CipherData *output_data) const {
  for (auto &e : data.x_diff) {
    output_data->x_diff.emplace_back();
    this->undo_last_round_one(e.first, kv, &output_data->x_diff.back().first);
    this->undo_last_round_one(e.second, kv, &output_data->x_diff.back().second);
  }
  for (auto &e : data.x) {
    output_data->x.emplace_back();
    this->undo_last_round_one(e, kv, &output_data->x.back());
  }
}

RangeCoverage FeistelBlock::range_lr_to_one(const RangeCoverage &l,
                                            const RangeCoverage &r) const {
  return l + r.shift(blk_size());
}

BitVec FeistelBlock::range_lr_to_bitvec(const RangeCoverage &l,
                                        const RangeCoverage &r) const {
  return BitVec::FromRange(range_lr_to_one(l, r), raw_input_size());
}

OPA_NAMESPACE_END(opa, crypto, la)
