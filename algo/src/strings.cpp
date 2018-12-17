#include <opa/algo/strings.h>

OPA_NAMESPACE_DECL2(opa, algo)
std::vector<int> kmp_matches(const glib::StringPiece &hay,
                             const glib::StringPiece &needle) {

  std::vector<int> res;

  std::vector<int> fail(hay.size() + 2);
  fail[0] = -1;
  fail[1] = 0;
  int p = 0;
  FOR (i, 1, needle.size()) {
    while (p != -1 && needle[p] != needle[i]) p = fail[p];
    ++p;
    fail[i + 1] = fail[p];
  }

  p = 0;
  int nn = needle.size();
  REP(i, hay.size()){
    while (p != -1 && needle[p] != hay[i]) p = fail[p];
    ++p;
    if (p == nn) res.push_back(i+1-nn), p = fail[p];
  }

  return res;
}

OPA_NAMESPACE_DECL2_END

