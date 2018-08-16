#include <stdint.h>
#include <string.h>

typedef uint8_t u8;
typedef int32_t s32;
typedef int64_t limb;

int curve25519_donna(const u8 *secret, const u8 *basepoint, u8 *mypublic);
