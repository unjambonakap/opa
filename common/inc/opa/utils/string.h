#pragma once

#include "opa/utils/string_base.h"
#include <opa/stolen/StringRef.h>
#include <opa_common.h>
#include <absl/strings/string_view.h>
#include <opa_inc.h>
#include <string>
namespace opa {
using absl::string_view;
}


OPA_NAMESPACE_DECL2(opa, utils)

std::string StringPrintf(string_view fmt, ...);

OPA_NAMESPACE_DECL2_END
