#include <opa/math/common/base/value.h>

#include <opa/math/common/base/range_coverage.h>

OPA_NAMESPACE(opa, math, common)

Value Value::FromBytes(const RangeCoverage &range, opa::stolen::StringRef s) {
  Value res(range.all().size());
  for (auto i : range.all()) res.setb(i, s[i / 8] >> (i & 7) & 1);
  return res;
}
void Value::sxor(const RangeCoverage &x) {
  for (auto i : x.all()) toggleb(i);
}

int Value::dot(const RangeCoverage &x) const {
  int res = 0;
  for (auto i : x.all()) res ^= getb(i);
  return res;
}
Value Value::and_(const RangeCoverage &x) const {
  Value res(size());
  REP(i, size())
  if (x.in(i)) res.setb(i, getb(i));

  return res;
}

OPA_NAMESPACE_END(opa, math, common)
