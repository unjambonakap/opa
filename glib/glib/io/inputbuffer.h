/* Copyright 2015 The TensorFlow Authors. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#ifndef TENSORFLOW_LIB_IO_INPUTBUFFER_H_
#define TENSORFLOW_LIB_IO_INPUTBUFFER_H_

#include <string>
#include "glib/core/coding.h"
#include "glib/core/status.h"
#include "glib/platform/env.h"
#include "glib/platform/macros.h"
#include "glib/platform/types.h"

namespace glib {
namespace io {

// An InputBuffer provides a buffer on top of a RandomAccessFile.
// A given instance of an InputBuffer is NOT safe for concurrent use
// by multiple threads
class InputBuffer {
 public:
  // Create an InputBuffer for "file" with a buffer size of
  // "buffer_bytes" bytes.  'file' must outlive *this.
  InputBuffer(RandomAccessFile* file, size_t buffer_bytes);
  ~InputBuffer();

  // Read one text line of data into "*result" until end-of-file or a
  // \n is read.  (The \n is not included in the result.)  Overwrites
  // any existing data in *result.
  //
  // If successful, returns OK.  If we are already at the end of the
  // file, we return an OUT_OF_RANGE error.  Otherwise, we return
  // some other non-OK status.
  Status ReadLine(string* result);

  // Reads bytes_to_read bytes into *result, overwriting *result.
  //
  // If successful, returns OK.  If we there are not enough bytes to
  // read before the end of the file, we return an OUT_OF_RANGE error.
  // Otherwise, we return some other non-OK status.
  Status ReadNBytes(int64 bytes_to_read, string* result);

  // An overload that writes to char*.  Caller must ensure result[0,
  // bytes_to_read) is valid to be overwritten.  Returns OK iff "*bytes_read ==
  // bytes_to_read".
  Status ReadNBytes(int64 bytes_to_read, char* result, size_t* bytes_read);

  // Reads a single varint32.
  Status ReadVarint32(uint32* result);

  // Like ReadNBytes() without returning the bytes read.
  Status SkipNBytes(int64 bytes_to_skip);

  // Seek to this offset within the file.
  //
  // If we seek to somewhere within our pre-buffered data, we will re-use what
  // data we can.  Otherwise, Seek() throws out the current buffer and the next
  // read will trigger a File::Read().
  Status Seek(int64 position);

  // Returns the position in the file.
  int64 Tell() const { return file_pos_ - (limit_ - pos_); }

  // Returns the underlying RandomAccessFile.
  RandomAccessFile* file() const { return file_; }

 private:
  Status FillBuffer();

  // Internal slow-path routine used by ReadVarint32().
  Status ReadVarint32Fallback(uint32* result);

  RandomAccessFile* file_;  // Not owned
  int64 file_pos_;          // Next position to read from in "file_"
  size_t size_;             // Size of "buf_"
  char* buf_;               // The buffer itself
  // [pos_,limit_) hold the "limit_ - pos_" bytes just before "file_pos_"
  char* pos_;    // Current position in "buf"
  char* limit_;  // Just past end of valid data in "buf"

  TF_DISALLOW_COPY_AND_ASSIGN(InputBuffer);
};

// Implementation details.

// Inlined for performance.
inline Status InputBuffer::ReadVarint32(uint32* result) {
  if (pos_ + core::kMaxVarint32Bytes <= limit_) {
    // Fast path: directly parse from buffered data.
    // Reads strictly from the range [pos_, limit_).
    const char* offset = core::GetVarint32Ptr(pos_, limit_, result);
    if (offset == nullptr) return errors::OutOfRange("Parsed past limit.");
    pos_ = const_cast<char*>(offset);
    return Status::OK();
  } else {
    return ReadVarint32Fallback(result);
  }
}

}  // namespace io
}  // namespace glib

#endif  // TENSORFLOW_LIB_IO_INPUTBUFFER_H_
