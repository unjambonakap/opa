#pragma once

#include "opa/utils/string_base.h"
#include <opa/stolen/StringRef.h>
#include <opa_common.h>
#include <opa_inc.h>
#include <string>
#include <glib/core/stringpiece.h>
namespace opa {
using glib::StringPiece;
}


OPA_NAMESPACE_DECL2(opa, utils)

std::string StringPrintf(StringPiece fmt, ...);

OPA_NAMESPACE_DECL2_END
