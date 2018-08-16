#include <opa/embed/c_struct/ring_buffer.h>

#include <assert.h>
#include <string.h>
#define CHECK_BOUNDS assert(this->n >= 0 && this->n <= this->sz);

void RingBuffer_init(RingBuffer *this, u8 *buf, u16 sz) {
  this->buf = buf;
  this->sz = sz;
  this->n = 0;
  this->pos = 0;
  this->right = 0;
}

bool RingBuffer_is_full(RingBuffer *this) { return this->n == this->sz; }
bool RingBuffer_is_empty(RingBuffer *this) { return this->n == 0; }
u16 RingBuffer_free(RingBuffer *this) { return this->sz - this->n; }
u16 RingBuffer_used(RingBuffer *this) { return this->n; }

u16 RingBuffer_remap(RingBuffer *this, u16 pos) {
  while (pos >= this->sz) pos -= this->sz;
  return pos;
}

bool RingBuffer_push(RingBuffer *this, u8 v) {
  if (RingBuffer_is_full(this)) return false;
  this->buf[this->right] = v;
  RingBuffer_commit_write(this, 1);
  return true;
}
u8 RingBuffer_top(RingBuffer *this) { return this->buf[this->pos]; }

u8 RingBuffer_pop(RingBuffer *this) {
  u8 res;
  if (RingBuffer_is_empty(this)) return 0xff;
  res = RingBuffer_top(this);
  RingBuffer_commit_read(this, 1);
  return res;
}

u16 RingBuffer_size(RingBuffer *this) { return this->n; }
u16 RingBuffer_capacity(RingBuffer *this) { return this->sz; }
u8 RingBuffer_get(RingBuffer *this, u16 pos) {
  return this->buf[RingBuffer_remap(this, this->pos + pos)];
}

void RingBuffer_commit_write(RingBuffer *this, u16 count) {
  this->n += count;
  CHECK_BOUNDS;
  this->right = RingBuffer_remap(this, this->right + count);
}

void RingBuffer_commit_read(RingBuffer *this, u16 count) {
  this->n -= count;
  CHECK_BOUNDS;
  this->pos = RingBuffer_remap(this, this->pos + count);
}

u8 *RingBuffer_get_write_buf(RingBuffer *this) {
  return this->buf + this->right;
}

u16 RingBuffer_get_write_size(RingBuffer *this) {
  return MIN(this->sz - this->n, this->sz - this->right);
}

u8 *RingBuffer_get_read_buf(RingBuffer *this) { return this->buf + this->pos; }

u16 RingBuffer_get_read_size(RingBuffer *this) {
  return MIN(this->n, this->sz - this->pos);
}

u16 RingBuffer_fill(RingBuffer *this, const u8 *buf, u32 len) {
  u16 avail_len = RingBuffer_get_write_size(this);
  avail_len = MIN(avail_len, len);
  memcpy(RingBuffer_get_write_buf(this), buf, avail_len);
  RingBuffer_commit_write(this, avail_len);
  return avail_len;
}

u16 RingBuffer_read_to_buf(RingBuffer *this, u8 *buf, u32 len) {
  const u8 *read_buf = RingBuffer_get_read_buf(this);
  u16 read_size = RingBuffer_get_read_size(this);
  read_size = MIN(read_size, len);
  memcpy(buf, read_buf, read_size);
  RingBuffer_commit_read(this, read_size);
  return read_size;
}
