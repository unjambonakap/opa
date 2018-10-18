%module opa_crypto_linear_cryptanalysis_swig

%include "math_common_swig_base.i"

%{
#include "opa/crypto/la/context.h"
#include "opa/crypto/la/base.h"
#include "opa/crypto/la/blocks.h"

using namespace opa::crypto::la;
%}

namespace std {
  %template(vu8) std::vector<u8>;
}


%include "opa/crypto/la/context.h"
%include "opa/crypto/la/base.h"
%include "opa/crypto/la/blocks.h"



%apply std::vector<u8> {std::string* };
