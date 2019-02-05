#pragma once

#include <opa/crypto/base.h>

OPA_NM_CRYPTO

std::string pkcs7(opa::stolen::StringRef s, int sz);
std::string pkcs1(opa::stolen::StringRef s, int sz);
std::string zeropad(opa::stolen::StringRef s, int sz);
std::pair<std::string, bool> rpkcs7(opa::stolen::StringRef s, int bs);
std::pair<std::string, bool> rpkcs1(opa::stolen::StringRef s);

OPA_NM_CRYPTO_END
