#include <plll/config.hpp>

#undef DEBUG // Otherwise, dlalloc contains a function (do_check_malloc_state) which is never used
#include "dlalloc.c"
