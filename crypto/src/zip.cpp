#include "zip.h"
#include <opa/utils/string.h>

using namespace std;

OPA_NM_CRYPTO_CRACKER

static u32 crc_32_table[256];

/*
 * Fill the CRC table, if not done already.
 */
static void make_crc_tab() {
  u32 s, t, v;
  static bool done = false;

  if (done)
    return;
  for (t = 0; t < 256; t++) {
    v = t;
    for (s = 0; s < 8; s++)
      v = (v >> 1) ^ ((v & 1) * (u32)0xedb88320L);
    crc_32_table[t] = v;
  }
  done = true;
}

#define CRC32(c, b) (crc_32_table[((int)(c) ^ (b)) & 0xff] ^ ((c) >> 8))

/*
 * Return the next byte in the pseudo-random sequence.
 */
#define DECRYPT_BYTE_ZIP(keys, t)                                              \
  {                                                                            \
    u16 temp = (u16)keys[2] | 2;                                               \
    t = (((temp * (temp ^ 1U)) >> 8) & 0xff);                                  \
  }

/*
 * Update the encryption keys with the next byte of plain text.
 */
#define UPDATE_KEYS_ZIP(keys, c)                                               \
  {                                                                            \
    keys[0] = CRC32(keys[0], (c));                                             \
    keys[1] += keys[0] & 0xff;                                                 \
    keys[1] = keys[1] * 134775813L + 1;                                        \
    keys[2] = CRC32(keys[2], (int)(keys[1] >> 24));                            \
  }

void ZipState::reset(const u8 *buf, int n) {
  make_crc_tab();
  m_keys[0] = 305419896L;
  m_keys[1] = 591751049L;
  m_keys[2] = 878082192L;
  REP (i, n)
    UPDATE_KEYS_ZIP(m_keys, buf[i]);
}

void ZipState::encrypt(const u8 *src, u8 *dest, int n) {
  u16 t;

  REP (i, n) {
    DECRYPT_BYTE_ZIP(m_keys, t);
    UPDATE_KEYS_ZIP(m_keys, src[i]);
    dest[i] = t ^ src[i];
  }
}

void ZipState::decrypt(const u8 *src, u8 *dest, int n) {
  u16 t;
  REP (i, n) {
    DECRYPT_BYTE_ZIP(m_keys, t);
    dest[i] = src[i] ^ t;
    UPDATE_KEYS_ZIP(m_keys, dest[i]);
  }
}

OPA_NM_CRYPTO_CRACKER_END
