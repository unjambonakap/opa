#include <opa/crypto/la/tester.h>
#include <opa/crypto/la/base.h>

OPA_NAMESPACE(opa, crypto, la)

void test_fast_eval(int nr, const CipherBlock *block) {
  REP (i, nr) {
    BitVec ov1, ov2;
    BitVec iv = BitVec::rand(block->input_size());
    OPA_DISP("START TEST >> ", iv);
    ov1 = block->evaluate(iv);
    ov2 = block->fast_eval_fe(iv);
    OPA_CHECK_EQ(ov1, ov2, iv);
  }
}

OPA_NAMESPACE_END(opa, crypto, la)
