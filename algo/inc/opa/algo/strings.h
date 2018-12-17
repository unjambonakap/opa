#include <opa/algo/base.h>

OPA_NAMESPACE_DECL2(opa, algo)

std::vector<int> kmp_matches(const glib::StringPiece &hay,
                             const glib::StringPiece &needle);

OPA_NAMESPACE_DECL2_END
