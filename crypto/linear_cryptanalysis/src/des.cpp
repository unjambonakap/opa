#include <opa_common.h>
#include <opa/crypto/la/graph.h>
#include <opa/crypto/la/des.h>

using namespace std;
OPA_NAMESPACE(opa, crypto, la)

OPA_REG_BLOCK(DESBlock)

vector<u8> tsf_pydes_sbox(const vector<u8> &tb) {
  vector<u8> res(1 << 6);
  PermutationDesc din;
  PermutationDesc dout;
  din.init({ 4, 3, 2, 1, 5, 0 });
  dout.init({ 3, 2, 1, 0 });

  REP (i, 1 << 6) { res[din.get2(i)] = dout.get2(tb[i]); }
  return res;
}

void DESBlock::init() {
  /* tsf code from dumb pydes code
                Bn[pos] = (v & 8) >> 3
                Bn[pos + 1] = (v & 4) >> 2
                Bn[pos + 2] = (v & 2) >> 1
                Bn[pos + 3] = v & 1
                m = (B[j][0] << 1) + B[j][5]
                n = (B[j][1] << 3) + (B[j][2] << 2) + (B[j][3] << 1) + B[j][4]
                */

  m_sbox[0].init(4, {
                      14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7, 0,
                      15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8, 4, 1,
                      14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0, 15, 12, 8,
                      2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13,
                    });
  m_sbox[1].init(4, { 15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10, 3,
                      13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5, 0, 14,
                      7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15, 13, 8, 10,
                      1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9 });
  m_sbox[2].init(4, { 10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8, 13,
                      7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1, 13, 6,
                      4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7, 1, 10, 13,
                      0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12 });
  m_sbox[3].init(4, { 7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15, 13,
                      8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9, 10, 6,
                      9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4, 3, 15, 0,
                      6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14 });
  m_sbox[4].init(4, { 2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9, 14,
                      11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6, 4, 2, 1,
                      11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14, 11, 8, 12, 7,
                      1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3 });
  m_sbox[5].init(4, { 12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11, 10,
                      15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8, 9, 14,
                      15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6, 4, 3, 2,
                      12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13 });
  m_sbox[6].init(4, { 4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1, 13,
                      0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6, 1, 4,
                      11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2, 6, 11, 13,
                      8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12 });
  m_sbox[7].init(4, { 13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7, 1,
                      15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2, 7, 11,
                      4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8, 2, 1, 14,
                      7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11 });
  m_ed.init({ 31, 0, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 7, 8, 9, 10, 11, 12, 11, 12,
              13, 14, 15, 16, 15, 16, 17, 18, 19, 20, 19, 20, 21, 22, 23, 24,
              23, 24, 25, 26, 27, 28, 27, 28, 29, 30, 31, 0 });

  REP (i, 8)
    m_sbox[i].mp = tsf_pydes_sbox(m_sbox[i].mp);

  m_perm.init({ 15, 6, 19, 20, 28, 11, 27, 16, 0, 14, 22, 25, 4, 17, 30, 9, 1,
                7, 23, 13, 31, 26, 2, 8, 18, 12, 29, 5, 21, 10, 3, 24 },
              true);

  CipherGraph::init(Params()
                      .call(CallInfo(blk_size, keysize, true))
                      .input_size(blk_size)
                      .output_size(blk_size)
                      .key_size(keysize)
                      .fast_eval(true));
  desc() = opa::utils::SPrintf("DesBlock()");
}

void DESBlock::init_graph() {

  m_eb = context()->instanciate<ExpansionBlock>();
  m_eb->init(raw_input_size(), m_ed);
  auto expn = add_node(m_eb.get());
  smart_plug({ input }, expn);
  std::vector<PlugDesc> nodes;

  m_xor = context()->instanciate<XorBlock>();
  m_xor->init(key_size());
  auto xorn = add_node(m_xor.get());
  smart_plug({ expn, key }, xorn);

  REP (i, 8) {
    m_sboxb[i] = context()->instanciate<SboxBlock>();
    m_sboxb[i]->init(m_sbox[i]);
    auto node = add_node(m_sboxb[i].get());
    smart_plug({ PlugDesc(xorn, i * 6, 6) }, node);
    nodes.emplace_back(node);
  }

  m_permb = context()->instanciate<PermutationBlock>();
  m_permb->init(m_perm, false);
  auto permn = add_node(m_permb.get());
  smart_plug(nodes, permn);
  smart_plug({ permn }, output);
}

u64 DESBlock::fast_eval2(u64 a, u64 b) const {
  a = m_ed.get(a);
  a ^= b;
  u64 res = 0;
  REP (i, 8) { res |= u64(m_sbox[i].get(a >> 6 * i & 0x3f)) << 4 * i; }
  res = m_perm.get2(res);
  return res;
}

OPA_NAMESPACE_END(opa, crypto, la)
