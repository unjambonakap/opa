%module opa_crypto_swig

%include "math_common_swig_base.i"

%apply (const u8 *src, int n) { (const u8 *charset, int n) }
%apply (const u8 *src, int n) { (const u8 *data, int n) }
%apply (const u8 *src) { (const u8 *data) }

%{
#include "opa/utils/serialize.h"

#include "opa/crypto/base.h"
#include "opa/crypto/cracker.h"
#include "opa/math/common/Ring.h"
#include "opa/math/common/Field.h"
#include "opa/math/common/GF_p.h"
#include "opa/math/common/bignum.h"
#include "opa/crypto/zip.h"
#include "opa/crypto/hash.h"
#include "opa/crypto/misc.h"
#include "opa/crypto/cracker_job.h"
#include "opa/crypto/cracker_misc.h"
#include "opa/crypto/lfsr.h"
#include "opa/crypto/lfsr_small.h"
#include "opa/crypto/aes.h"
#include "opa/crypto/padding.h"

using namespace opa::crypto;
using namespace opa::utils;
using namespace opa::crypto::cracker;
%}


%shared_ptr(opa::utils::ProtobufParams);
%shared_ptr(opa::utils::BaseStorable);
%shared_ptr(opa::utils::BaseStorage);
%shared_ptr(opa::utils::StrStrMapperFunc);
%shared_ptr(opa::crypto::cracker::ShardParams)
%shared_ptr(opa::crypto::cracker::Res)
%shared_ptr(opa::crypto::Md5)
%shared_ptr(opa::crypto::Sha1)
%shared_ptr(opa::crypto::Sha256)
%shared_ptr(opa::crypto::Sha512)
%shared_ptr(opa::crypto::Hash)
%shared_ptr(opa::crypto::cracker::PatternChecker)
%shared_ptr(opa::crypto::cracker::PatternCheckerParams)
%shared_ptr(opa::crypto::cracker::MultiplePatternChecker)
%shared_ptr(opa::crypto::cracker::MultiplePatternCheckerParams)
%shared_ptr(opa::crypto::cracker::CrackerParams)
%shared_ptr(opa::crypto::cracker::MapperAndCheckerParams)
%shared_ptr(opa::crypto::cracker::MapperAndChecker)
%shared_ptr(opa::crypto::cracker::CrackerChecker)
%shared_ptr(opa::crypto::cracker::Cracker)
%shared_ptr(opa::crypto::cracker::Pattern)
%shared_ptr(opa::utils::DoubleRes)
%shared_ptr(opa::utils::MapperFunc<opa::stolen::StringRef, std::string>)
%shared_ptr(opa::crypto::LFSR_GF2_small)
%shared_ptr(opa::crypto::LFSR<u32>)


%include "opa/utils/serialize.h"
%include "opa/crypto/base.h"

namespace std{
  %template(PatternVector) std::vector<opa::crypto::cracker::Pattern>;
  %template(PatternCheckerParamsVector) std::vector<opa::crypto::cracker::PatternCheckerParams>;
}

%apply const ::std::vector<opa::crypto::cracker::Pattern>& {::std::vector<opa::crypto::cracker::Pattern>* };



%rename("%(regex:/^server_.*/$ignore/)s") "";
%rename("%(regex:/^worker_.*/$ignore/)s") "";
%rename("%(regex:/^worker_.*/$ignore/)s") "";
%ignore "Build_StatusMsg";
//%ignore "opa::crypto::LFSR<T>::rand";
%include "opa/threading/data.h"
%include "opa/threading/job.h"
%include "opa/threading/auto_job.h"
%include "opa/threading/dispatcher.h"
%include "opa/threading/runner.h"

namespace opa{
namespace threading{
//%template(JAMBON) AutoVectorJob<opa::crypto::cracker::ShardParams, opa::crypto::cracker::Res>;
}
}

namespace opa{
namespace utils{
%template(StrStrMapperFuncSwig) MapperFunc<opa::stolen::StringRef, std::string>;
}
}


#%template(non_owned_sptr_MapperAndCheker) opa::crypto::cracker::CrackerChecker::non_owned_sptr<opa::crypto::cracker::MapperAndChecker>;
#%template(non_owned_sptr_Md5) opa::crypto::cracker::CrackerChecker::non_owned_sptr<opa::crypto::Md5>;


%include "opa/stolen/StringRef.h"
%include "opa/crypto/cracker.h"
%include "opa/crypto/zip.h"
%include "opa/crypto/hash.h"
%include "opa/crypto/misc.h"
%include "opa/crypto/cracker_job.h"
%include "opa/crypto/cracker_misc.h"
%include "opa/crypto/lfsr.h"
%include "opa/crypto/aes.h"
%include "opa/crypto/padding.h"

namespace opa{
namespace crypto{
%template(LFSR_u32) opa::crypto::LFSR<u32>;
}
}

%include "opa/crypto/lfsr_small.h"


