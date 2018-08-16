#ifndef _H_LIB_STRUCT_RING_BUFFER_H_
#define _H_LIB_STRUCT_RING_BUFFER_H_

#include <stdint.h>
#include <stdbool.h>
#include <opa_inc.h>

typedef struct RingBuffer {
  volatile u16 pos;
  volatile u16 right;
  volatile u16 n;
  volatile u16 sz;
  u8 *buf;

} RingBuffer;

void RingBuffer_init(RingBuffer *this, u8 *buf, u16 sz);

bool RingBuffer_is_full(RingBuffer *this);
bool RingBuffer_is_empty(RingBuffer *this);
u16 RingBuffer_free(RingBuffer *this);
u16 RingBuffer_used(RingBuffer *this);
u16 RingBuffer_size(RingBuffer *this);
u16 RingBuffer_capacity(RingBuffer *this);
void RingBuffer_advance(RingBuffer *this, u16 count);

u16 RingBuffer_remap(RingBuffer *this, u16 pos);

bool RingBuffer_push(RingBuffer *this, u8 v);
u8 RingBuffer_top(RingBuffer *this);

u8 RingBuffer_pop(RingBuffer *this);
u8 RingBuffer_get(RingBuffer *this, u16 pos);

void RingBuffer_commit_write(RingBuffer *this, u16 count);
void RingBuffer_commit_read(RingBuffer *this, u16 count);
u8 *RingBuffer_get_write_buf(RingBuffer *this);
u16 RingBuffer_get_write_size(RingBuffer *this);
u8 *RingBuffer_get_read_buf(RingBuffer *this);
u16 RingBuffer_get_read_size(RingBuffer *this);
u16 RingBuffer_fill(RingBuffer *this, const u8 *buf, u32 len);
u16 RingBuffer_read_to_buf(RingBuffer *this, u8 *buf, u32 len);

#endif // _H_LIB_STRUCT_RING_BUFFER_H_
