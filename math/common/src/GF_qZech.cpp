#include "common/GF_qZech.h"
#include "common/Utils.h"
#include "common/PolyRing.h"

OPA_NAMESPACE_DECL3(opa, math, common)

GF_qZech::GF_qZech(u32 size, u32 car, const std::vector<u32> &logTable,
                   const std::vector<u32> &u32Table)
    : Field<u32>(size, car) {

    activateU32();

    m_qm1 = size - 1;
    zechslog = (u32 *)malloc((m_qm1 + 1) * sizeof(u32));
    lst_u32 = (u32 *)malloc(getCar().getu32() * sizeof(u32));

    for (u32 i = 0; i < logTable.size(); ++i)
        zechslog[i] = logTable[i];
    for (u32 i = 0; i < u32Table.size(); ++i)
        lst_u32[i] = u32Table[i];
}

GF_qZech::~GF_qZech() {

    free(zechslog);
    free(lst_u32);
}

u32 GF_qZech::mul(const u32 &a, const u32 &b) const {
    if (a == 0 || b == 0)
        return 0;
    return (a + b - 2) % m_qm1 + 1;
}

u32 GF_qZech::add(const u32 &a, const u32 &b) const {
    if (a == 0)
        return b;
    if (b == 0)
        return a;
    u32 c = (a - b + m_qm1) % m_qm1 + 1;
    return mul(b, zechslog[c]);
}

u32 GF_qZech::inv(const u32 &a) const {
    assert(a != 0);
    return m_qm1 - a + 2;
}

u32 GF_qZech::neg(const u32 &a) const {
    if (a == 0)
        return 0;
    if (getCar() == 2)
        return a;
    return mul(1 + m_qm1 / 2, a);
}

bool GF_qZech::isZ(const u32 &a) const { return a == 0; }

bool GF_qZech::isE(const u32 &a) const { return a == 1; }

u32 GF_qZech::getZ() const { return 0; }

u32 GF_qZech::getE() const { return 1; }

bignum GF_qZech::getOrderOf(const u32 &a) const {
    assert(a != 0);
    return m_qm1 / gcd(m_qm1, a - 1);
}

u32 GF_qZech::importu32(const u32 &a) const {
    return lst_u32[a % getCar().getu32()];
}

u32 GF_qZech::import(const u32 &a) const { return a % getSizeU32(); }

OPA_NAMESPACE_DECL3_END
