#include <opa/math/common/Poly.h>
#include <opa/utils/string.h>

OPA_NAMESPACE_DECL3(opa, math, common)

std::string get_poly_letter(int id){
    assert(id<26);
    return std::string(1, 'A'+id);
}

OPA_NAMESPACE_DECL3_END


