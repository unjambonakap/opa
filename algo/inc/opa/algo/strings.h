#include <opa/algo/base.h>
#include <string_view>

OPA_NAMESPACE_DECL2(opa, algo)

std::vector<int> kmp_matches(const std::string_view &hay,
                             const std::string_view &needle);

OPA_NAMESPACE_DECL2_END
