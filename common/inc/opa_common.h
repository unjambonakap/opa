#pragma once

#include <opa/predef.h>
#include <opa_common_base.h>
#include <opa/utils/stacktrace.h>
#if OPA_SWIG
#define DECLARE_int32(...)
#endif

#include <gflags/gflags.h>
#include <glib/core/stringpiece.h>
#include <glib/gtl/map_util.h>
#include <glog/logging.h>
#include <opa/stolen/StringRef.h>

namespace opa {
using glib::StringPiece;
}

OPA_NM_INIT

void opa_init(int argc, char **argv);
class OpaInitRegister {
public:
  typedef void (*InitFunc)();
  OpaInitRegister(InitFunc init_func, const std::string &name);
};

OPA_NM_INIT_END

#define OPA_REGISTER_INIT(name, func)                                          \
  static ::opa::init::OpaInitRegister RegisterCl_##name(func, #name);
