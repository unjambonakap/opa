#pragma once
#include <opa_common.h>
#include <opa/crypto/la/graph.h>
#include <opa/crypto/la/relation_driver.h>
#include <opa/crypto/la/blocks.h>

OPA_NAMESPACE_DECL3(opa, crypto, la)

template <int N> class FealGBlockN : public CipherGraph {
public:
  static constexpr int blk_size = N;
  virtual void init(int mod) {
    this->mod = mod;
    OPA_DISP("GOGO ", blk_size, mod);
    CipherGraph::init(Params()
                        .call(CallInfo(blk_size, 0, true))
                        .input_size(2 * blk_size)
                        .output_size(blk_size)
                        .fast_eval(true));
    auto dr = ((SimpleRelDriver *)driver());
    dr->fini();
    dr->init(0.5, 100);
  }
  virtual void init_graph() override;

  virtual u64 fast_eval2(u64 a, u64 b) const override {
    return opa::utils::rotl<N>((a + b + mod) & ((1 << blk_size) - 1), 3);
  }

  OPA_TGEN_CIPHERBLOCK(int)
  UPTR(AddBlock) add1;
  UPTR(AddConstBlock) add2;
  UPTR(RotationBlock) rotb;
  OPA_ACCESSOR_R(int, mod, params);

  int mod;
};

typedef FealGBlockN<8> FealGBlock;

template <int N> class FealFBoxBlockN : public CipherGraph {
public:
  static constexpr int blk_size = N;
  static constexpr int base_size = N / 4;
  virtual void init() override {
    CipherGraph::init(Params()
                        .call(CallInfo(blk_size, blk_size, true))
                        .input_size(blk_size)
                        .output_size(blk_size)
                        .key_size(blk_size)
                        .fast_eval(true));
    auto dr = ((SimpleRelDriver *)driver());
    dr->fini();
    dr->init(0.5, 100);
  }

  virtual void init_graph() override;

  virtual u64 fast_eval2(u64 a, u64 b) const override;
  OPA_TGEN_CIPHERBLOCK0()

  UPTR(XorBlock) xor8;
  UPTR(XorBlock) xor32;
  UPTR(FealGBlockN<base_size>) g0;
  UPTR(FealGBlockN<base_size>) g1;
};

typedef FealFBoxBlockN<32> FealFBoxBlock;

union u32_un {
  u32 v32;
  u16 v16[2];
  u8 v8[4];
};

template <int N> void FealGBlockN<N>::init_graph() {
  add1 = context()->template instanciate<AddBlock>();
  add1->init(blk_size);

  CipherNode *cur;
  CipherNode *tmp;

  cur = add_node(add1.get());
  smart_plug({ input }, cur);

  if (mod) {
    add2 = context()->template instanciate<AddConstBlock>();
    add2->init(mod, blk_size);
    tmp = add_node(add2.get());
    smart_plug({ cur }, tmp);
    cur = tmp;
  }

  rotb = context()->template instanciate<RotationBlock>();
  rotb->init(blk_size, 3);
  tmp = add_node(rotb.get());

  smart_plug({ cur }, tmp);
  smart_plug({ tmp }, output);
  desc() = opa::utils::SPrintf("FealG<%d>(%d)", N, mod);
}

template <int N> void FealFBoxBlockN<N>::init_graph() {

  CipherNode *cur;
  CipherNode *tmp;

  xor32 = context()->template instanciate<XorBlock>();
  xor32->init(blk_size);

  xor8 = context()->template instanciate<XorBlock>();
  xor8->init(base_size);

  g0 = context()->template instanciate<FealGBlockN<base_size> >();
  g0->init(0);
  g1 = context()->template instanciate<FealGBlockN<base_size> >();
  g1->init(1);

  cur = add_node(xor32.get());
  smart_plug({ input, key }, cur);

  auto t0 = add_node(xor8.get());
  smart_plug({ PlugDesc(cur, base_size * 2, base_size),
               PlugDesc(cur, base_size * 3, base_size) },
             t0);

  auto t1 = add_node(xor8.get());
  smart_plug(
    { PlugDesc(cur, 0, base_size), PlugDesc(cur, base_size, base_size) }, t1);

  auto y0 = add_node(g0.get());
  auto y1 = add_node(g1.get());
  auto y2 = add_node(g0.get());
  auto y3 = add_node(g1.get());

  smart_plug({ t1, t0 }, y1);
  smart_plug({ PlugDesc(cur, 0, base_size), y1 }, y0);
  smart_plug({ t0, y1 }, y2);
  smart_plug({ PlugDesc(cur, base_size * 3, base_size), y2 }, y3);

  smart_plug({ y3, y2, y1, y0 }, output);
  desc() = opa::utils::SPrintf("FealF");
}

template <int N> u64 FealFBoxBlockN<N>::fast_eval2(u64 a, u64 b) const {
  a ^= b;
  u32_un x;
  int mask = (1 << base_size) - 1;
  x.v8[0] = a & mask;
  x.v8[1] = a >> base_size & mask;
  x.v8[2] = a >> 2 * base_size & mask;
  x.v8[3] = a >> 3 * base_size & mask;
  // x.v32 = a;
  u8 t0 = x.v8[2] ^ x.v8[3];
  u8 y1 = g1->fast_eval2(x.v8[0] ^ x.v8[1], t0);
  u8 y0 = g0->fast_eval2(x.v8[0], y1);
  u8 y2 = g0->fast_eval2(t0, y1);
  u8 y3 = g1->fast_eval2(x.v8[3], y2);

  // def gBox(a,b,mode):
  //    return rot3((a+b+mode)%256)
  // def fBox(plain):
  //    t0 = (plain[2] ^ plain[3])
  //    y1 = gBox(plain[0] ^ plain[1], t0, 1)
  //    y0 = gBox(plain[0], y1, 0)
  //    y2 = gBox(t0, y1, 0)
  //    y3 = gBox(plain[3], y2, 1)

  u32_un res;
  res.v32 = y3 | y2 << base_size | y1 << 2 * base_size | y0 << 3 * base_size;
  // res.v8[0] = y3;
  // res.v8[1] = y2;
  // res.v8[2] = y1;
  // res.v8[3] = y0;
  return res.v32;
}

template <> u64 FealFBoxBlock::fast_eval2(u64 a, u64 b) const {
  a ^= b;
  u32_un x;
  x.v32 = a;
  u8 t0 = x.v8[2] ^ x.v8[3];
  u8 y1 = g1->fast_eval2(x.v8[0] ^ x.v8[1], t0);
  u8 y0 = g0->fast_eval2(x.v8[0], y1);
  u8 y2 = g0->fast_eval2(t0, y1);
  u8 y3 = g1->fast_eval2(x.v8[3], y2);

  // def gBox(a,b,mode):
  //    return rot3((a+b+mode)%256)
  // def fBox(plain):
  //    t0 = (plain[2] ^ plain[3])
  //    y1 = gBox(plain[0] ^ plain[1], t0, 1)
  //    y0 = gBox(plain[0], y1, 0)
  //    y2 = gBox(t0, y1, 0)
  //    y3 = gBox(plain[3], y2, 1)

  u32_un res;
  res.v8[0] = y3;
  res.v8[1] = y2;
  res.v8[2] = y1;
  res.v8[3] = y0;
  return res.v32;
}

OPA_NAMESPACE_DECL3_END
