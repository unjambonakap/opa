#include <opa_common.h>

OPA_NAMESPACE_DECL2(opa, crypto)

class HillCipher {
  public:
    virtual ~HillCipher() {}
    virtual bool setKey(const std::vector<u32> &m) = 0;
    virtual std::vector<u32> decrypt(const std::vector<u32> &cipher) = 0;
    virtual std::vector<u32> encrypt(const std::vector<u32> &plain) = 0;
    static HillCipher *get(u32 n, u32 mod);
};

OPA_NAMESPACE_DECL2_END
