#include <opa_common.h>
#include <opa/crypto/la/blocks.h>
#include <opa/crypto/la/relation_driver.h>
#include <opa/math/common/GF_q.h>
#include <opa/math/common/Types.h>

using namespace std;
using namespace opa::crypto::la;
using namespace opa::math::common;

class GF2Block : public CipherBlock {
public:
  virtual void init(u32 deg) {
    CipherBlock::init(Params()
                        .call(CallInfo(deg, 0, true))
                        .input_size(deg)
                        .output_size(deg)
                        .fast_eval(true)
                        .rel_driver(new SimpleRelDriver(0.001, 100)));
    desc() = opa::utils::SPrintf("Gf2block(sz=%d)", deg);
    this->field.init(&GF2, deg);
  }

  virtual void do_evaluate(const Value &iv, const Value &kv,
                           Value &ov) const override {
    OPA_CHECK0(false);
  }

  virtual void do_get_relations(std::vector<Relation> &rels) const override {
    return get_relations_walsh(rels);
  }

  virtual void do_get_pre_fail(const Basis &out, Basis &res) const override {
    OPA_CHECK0(false);
  }

  virtual void setup_jit(const JitBuilder &builder) const override {
    OPA_CHECK0(false);
  }

  virtual u64 fast_eval1(u64 a) const override {
    auto aa = field.import_base(a);
    return field.to_base(aa.powm(5)).getu32();
  }

  GF_q<u32> field;
};

int main(int argc, char **argv) {
  opa::init::opa_init(argc, argv);

  GF2Block blk;
  blk.init(5);
  auto rels = blk.get_relations();
  OPA_DISP0(rels.size());

  for (auto &rel : rels) {
    OPA_DISP0(rel);
  }

  return 0;
}
