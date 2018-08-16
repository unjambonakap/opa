#pragma once

#include <opa/math/common/Field.h>

OPA_NAMESPACE_DECL3(opa, math, common)

class GF_qZech : public Field<u32> {
  private:
    u32 m_qm1;
    u32 *zechslog;
    u32 *lst_u32;

  public:
    GF_qZech(u32 size, u32 car, const std::vector<u32> &logTable,
             const std::vector<u32> &u32Table);
    ~GF_qZech();

    virtual u32 mul(const u32 &a, const u32 &b) const;
    virtual u32 add(const u32 &a, const u32 &b) const;
    virtual u32 inv(const u32 &a) const;
    virtual u32 neg(const u32 &a) const;

    virtual bool isZ(const u32 &a) const;
    virtual bool isE(const u32 &a) const;
    virtual u32 import(const u32 &a) const;
    virtual u32 getPrimitiveElem() const { return 2; }

    virtual u32 getZ() const;
    virtual u32 getE() const;
    virtual u32 importu32(const u32 &a) const;
    virtual bignum getOrderOf(const u32 &a) const;
};

OPA_NAMESPACE_DECL3_END
