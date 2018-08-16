#include "cracker_job.h"

OPA_NM_CRYPTO_CRACKER
std::string CHARSET_UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
std::string CHARSET_LOWER = "abcdefghijklmnopqrstuvwxyz";
std::string CHARSET_ALPHA =
  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
std::string CHARSET_NUM = "0123456789";
std::string CHARSET_ALPHANUM =
  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
std::string CHARSET_ALL;

void init_charset() {
  CHARSET_ALL.resize(256);
  REP (i, 256)
    CHARSET_ALL[i] = i;
}

OPA_REGISTER_INIT(init_charset, init_charset);
OPA_CLOUDY_REGISTER_BASE(Cracker);
OPA_CLOUDY_JOB_IMPL(Cracker);

OPA_NM_CRYPTO_CRACKER_END
