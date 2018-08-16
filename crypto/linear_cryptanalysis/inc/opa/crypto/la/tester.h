#pragma once
#include <opa_common.h>

OPA_NAMESPACE(opa, crypto, la)

class CipherBlock;

void test_fast_eval(int nr, const CipherBlock *block);
OPA_NAMESPACE_END(opa, crypto, la)
