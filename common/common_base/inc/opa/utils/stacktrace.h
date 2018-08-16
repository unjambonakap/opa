// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

#pragma once

#include <cstdlib>
#include <stdio.h>
#include <execinfo.h>
#include <cxxabi.h>
#include <string>
#include <opa_common_base.h>

OPA_NAMESPACE_DECL2(opa, utils)

/** Print a demangled stack backtrace of the caller function to FILE* out. */
std::string get_stacktrace(unsigned int max_frames = 63);

OPA_NAMESPACE_DECL2_END

