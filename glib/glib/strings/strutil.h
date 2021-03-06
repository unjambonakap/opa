// Protocol Buffers - Google's data interchange format
// Copyright 2008 Google Inc.  All rights reserved.
// https://developers.google.com/protocol-buffers/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// from google3/strings/strutil.h

#ifndef GLIB_PROTOBUF_STUBS_STRUTIL_H__
#define GLIB_PROTOBUF_STUBS_STRUTIL_H__

#define LIBPROTOBUF_EXPORT
#include <stdlib.h>
#include <vector>
#include "glib/core/stringpiece.h"
#include "glib/strings/numbers.h"
#include "glib/platform/macros.h"
#include "glib/platform/types.h"

namespace Eigen {
struct half;
}

namespace glib {

// ===================================================================
// from google3/base/scoped_ptr.h

namespace internal {

//  This is an implementation designed to match the anticipated future TR2
//  implementation of the scoped_ptr class, and its closely-related brethren,
//  scoped_array, scoped_ptr_malloc, and make_scoped_ptr.

template <class C> class scoped_ptr;
template <class C> class scoped_array;

// A scoped_ptr<T> is like a T*, except that the destructor of scoped_ptr<T>
// automatically deletes the pointer it holds (if any).
// That is, scoped_ptr<T> owns the T object that it points to.
// Like a T*, a scoped_ptr<T> may hold either NULL or a pointer to a T object.
//
// The size of a scoped_ptr is small:
// sizeof(scoped_ptr<C>) == sizeof(C*)
template <class C>
class scoped_ptr {
 public:

  // The element type
  typedef C element_type;

  // Constructor.  Defaults to initializing with NULL.
  // There is no way to create an uninitialized scoped_ptr.
  // The input parameter must be allocated with new.
  explicit scoped_ptr(C* p = NULL) : ptr_(p) { }

  // Destructor.  If there is a C object, delete it.
  // We don't need to test ptr_ == NULL because C++ does that for us.
  ~scoped_ptr() {
    enum { type_must_be_complete = sizeof(C) };
    delete ptr_;
  }

  // Reset.  Deletes the current owned object, if any.
  // Then takes ownership of a new object, if given.
  // this->reset(this->get()) works.
  void reset(C* p = NULL) {
    if (p != ptr_) {
      enum { type_must_be_complete = sizeof(C) };
      delete ptr_;
      ptr_ = p;
    }
  }

  // Accessors to get the owned object.
  // operator* and operator-> will assert() if there is no current object.
  C& operator*() const {
    assert(ptr_ != NULL);
    return *ptr_;
  }
  C* operator->() const  {
    assert(ptr_ != NULL);
    return ptr_;
  }
  C* get() const { return ptr_; }

  // Comparison operators.
  // These return whether two scoped_ptr refer to the same object, not just to
  // two different but equal objects.
  bool operator==(C* p) const { return ptr_ == p; }
  bool operator!=(C* p) const { return ptr_ != p; }

  // Swap two scoped pointers.
  void swap(scoped_ptr& p2) {
    C* tmp = ptr_;
    ptr_ = p2.ptr_;
    p2.ptr_ = tmp;
  }

  // Release a pointer.
  // The return value is the current pointer held by this object.
  // If this object holds a NULL pointer, the return value is NULL.
  // After this operation, this object will hold a NULL pointer,
  // and will not own the object any more.
  C* release() {
    C* retVal = ptr_;
    ptr_ = NULL;
    return retVal;
  }

 private:
  C* ptr_;

  // Forbid comparison of scoped_ptr types.  If C2 != C, it totally doesn't
  // make sense, and if C2 == C, it still doesn't make sense because you should
  // never have the same object owned by two different scoped_ptrs.
  template <class C2> bool operator==(scoped_ptr<C2> const& p2) const;
  template <class C2> bool operator!=(scoped_ptr<C2> const& p2) const;

  // Disallow evil constructors
  scoped_ptr(const scoped_ptr&);
  void operator=(const scoped_ptr&);
};

// scoped_array<C> is like scoped_ptr<C>, except that the caller must allocate
// with new [] and the destructor deletes objects with delete [].
//
// As with scoped_ptr<C>, a scoped_array<C> either points to an object
// or is NULL.  A scoped_array<C> owns the object that it points to.
//
// Size: sizeof(scoped_array<C>) == sizeof(C*)
template <class C>
class scoped_array {
 public:

  // The element type
  typedef C element_type;

  // Constructor.  Defaults to initializing with NULL.
  // There is no way to create an uninitialized scoped_array.
  // The input parameter must be allocated with new [].
  explicit scoped_array(C* p = NULL) : array_(p) { }

  // Destructor.  If there is a C object, delete it.
  // We don't need to test ptr_ == NULL because C++ does that for us.
  ~scoped_array() {
    enum { type_must_be_complete = sizeof(C) };
    delete[] array_;
  }

  // Reset.  Deletes the current owned object, if any.
  // Then takes ownership of a new object, if given.
  // this->reset(this->get()) works.
  void reset(C* p = NULL) {
    if (p != array_) {
      enum { type_must_be_complete = sizeof(C) };
      delete[] array_;
      array_ = p;
    }
  }

  // Get one element of the current object.
  // Will assert() if there is no current object, or index i is negative.
  C& operator[](std::ptrdiff_t i) const {
    assert(i >= 0);
    assert(array_ != NULL);
    return array_[i];
  }

  // Get a pointer to the zeroth element of the current object.
  // If there is no current object, return NULL.
  C* get() const {
    return array_;
  }

  // Comparison operators.
  // These return whether two scoped_array refer to the same object, not just to
  // two different but equal objects.
  bool operator==(C* p) const { return array_ == p; }
  bool operator!=(C* p) const { return array_ != p; }

  // Swap two scoped arrays.
  void swap(scoped_array& p2) {
    C* tmp = array_;
    array_ = p2.array_;
    p2.array_ = tmp;
  }

  // Release an array.
  // The return value is the current pointer held by this object.
  // If this object holds a NULL pointer, the return value is NULL.
  // After this operation, this object will hold a NULL pointer,
  // and will not own the object any more.
  C* release() {
    C* retVal = array_;
    array_ = NULL;
    return retVal;
  }

 private:
  C* array_;

  // Forbid comparison of different scoped_array types.
  template <class C2> bool operator==(scoped_array<C2> const& p2) const;
  template <class C2> bool operator!=(scoped_array<C2> const& p2) const;

  // Disallow evil constructors
  scoped_array(const scoped_array&);
  void operator=(const scoped_array&);
};

}  // namespace internal

// We made these internal so that they would show up as such in the docs,
// but we don't want to stick "internal::" in front of them everywhere.
using internal::scoped_ptr;
using internal::scoped_array;

#ifdef _MSC_VER
#define strtoll  _strtoi64
#define strtoull _strtoui64
#elif defined(__DECCXX) && defined(__osf__)
// HP C++ on Tru64 does not have strtoll, but strtol is already 64-bit.
#define strtoll strtol
#define strtoull strtoul
#endif

// ----------------------------------------------------------------------
// ascii_isalnum()
//    Check if an ASCII character is alphanumeric.  We can't use ctype's
//    isalnum() because it is affected by locale.  This function is applied
//    to identifiers in the protocol buffer language, not to natural-language
//    strings, so locale should not be taken into account.
// ascii_isdigit()
//    Like above, but only accepts digits.
// ascii_isspace()
//    Check if the character is a space character.
// ----------------------------------------------------------------------

inline bool ascii_isalnum(char c) {
  return ('a' <= c && c <= 'z') ||
         ('A' <= c && c <= 'Z') ||
         ('0' <= c && c <= '9');
}

inline bool ascii_isdigit(char c) {
  return ('0' <= c && c <= '9');
}

inline bool ascii_isspace(char c) {
  return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' ||
      c == '\r';
}

inline bool ascii_isupper(char c) {
  return c >= 'A' && c <= 'Z';
}

inline bool ascii_islower(char c) {
  return c >= 'a' && c <= 'z';
}

inline char ascii_toupper(char c) {
  return ascii_islower(c) ? c - ('a' - 'A') : c;
}

inline char ascii_tolower(char c) {
  return ascii_isupper(c) ? c + ('a' - 'A') : c;
}

inline int hex_digit_to_int(char c) {
  /* Assume ASCII. */
  int x = static_cast<unsigned char>(c);
  if (x > '9') {
    x += 9;
  }
  return x & 0xf;
}

// ----------------------------------------------------------------------
// HasPrefixString()
//    Check if a string begins with a given prefix.
// StripPrefixString()
//    Given a string and a putative prefix, returns the string minus the
//    prefix string if the prefix matches, otherwise the original
//    string.
// ----------------------------------------------------------------------
inline bool HasPrefixString(const std::string& str,
                            const std::string& prefix) {
  return str.size() >= prefix.size() &&
         str.compare(0, prefix.size(), prefix) == 0;
}

inline string StripPrefixString(const std::string& str, const std::string& prefix) {
  if (HasPrefixString(str, prefix)) {
    return str.substr(prefix.size());
  } else {
    return str;
  }
}

// ----------------------------------------------------------------------
// HasSuffixString()
//    Return true if str ends in suffix.
// StripSuffixString()
//    Given a string and a putative suffix, returns the string minus the
//    suffix string if the suffix matches, otherwise the original
//    string.
// ----------------------------------------------------------------------
inline bool HasSuffixString(const std::string& str,
                            const std::string& suffix) {
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline string StripSuffixString(const std::string& str, const std::string& suffix) {
  if (HasSuffixString(str, suffix)) {
    return str.substr(0, str.size() - suffix.size());
  } else {
    return str;
  }
}

// ----------------------------------------------------------------------
// StripString
//    Replaces any occurrence of the character 'remove' (or the characters
//    in 'remove') with the character 'replacewith'.
//    Good for keeping html characters or protocol characters (\t) out
//    of places where they might cause a problem.
// StripWhitespace
//    Removes whitespaces from both ends of the given string.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT void StripString(string* s, const char* remove,
                                    char replacewith);

LIBPROTOBUF_EXPORT void StripWhitespace(string* s);


// ----------------------------------------------------------------------
// LowerString()
// UpperString()
// ToUpper()
//    Convert the characters in "s" to lowercase or uppercase.  ASCII-only:
//    these functions intentionally ignore locale because they are applied to
//    identifiers used in the Protocol Buffer language, not to natural-language
//    strings.
// ----------------------------------------------------------------------

inline void LowerString(std::string * s) {
  std::string::iterator end = s->end();
  for (std::string::iterator i = s->begin(); i != end; ++i) {
    // tolower() changes based on locale.  We don't want this!
    if ('A' <= *i && *i <= 'Z') *i += 'a' - 'A';
  }
}

inline void UpperString(std::string * s) {
  std::string::iterator end = s->end();
  for (std::string::iterator i = s->begin(); i != end; ++i) {
    // toupper() changes based on locale.  We don't want this!
    if ('a' <= *i && *i <= 'z') *i += 'A' - 'a';
  }
}

inline std::string ToUpper(const std::string& s) {
  std::string out = s;
  UpperString(&out);
  return out;
}

// ----------------------------------------------------------------------
// StringReplace()
//    Give me a string and two patterns "old" and "new", and I replace
//    the first instance of "old" in the string with "new", if it
//    exists.  RETURN a new string, regardless of whether the replacement
//    happened or not.
// ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT string StringReplace(const std::string& s, const std::string& oldsub,
                                        const std::string& newsub, bool replace_all);

// ----------------------------------------------------------------------
// SplitStringUsing()
//    Split a string using a character delimiter. Append the components
//    to 'result'.  If there are consecutive delimiters, this function skips
//    over all of them.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT void SplitStringUsing(const std::string& full, const char* delim,
                                         std::vector<std::string>* res);

// Split a string using one or more byte delimiters, presented
// as a nul-terminated c string. Append the components to 'result'.
// If there are consecutive delimiters, this function will return
// corresponding empty strings.  If you want to drop the empty
// strings, try SplitStringUsing().
//
// If "full" is the empty string, yields an empty string as the only value.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT void SplitStringAllowEmpty(const std::string& full,
                                              const char* delim,
                                              std::vector<std::string>* result);

// ----------------------------------------------------------------------
// Split()
//    Split a string using a character delimiter.
// ----------------------------------------------------------------------
inline std::vector<std::string> Split(
    const std::string& full, const char* delim, bool skip_empty = true) {
  std::vector<std::string> result;
  if (skip_empty) {
    SplitStringUsing(full, delim, &result);
  } else {
    SplitStringAllowEmpty(full, delim, &result);
  }
  return result;
}

// ----------------------------------------------------------------------
// JoinStrings()
//    These methods concatenate a std::vector of strings into a C++ string, using
//    the C-string "delim" as a separator between components. There are two
//    flavors of the function, one flavor returns the concatenated string,
//    another takes a pointer to the target string. In the latter case the
//    target string is cleared and overwritten.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT void JoinStrings(const std::vector<std::string>& components,
                                    const char* delim, std::string* result);

inline std::string JoinStrings(const std::vector<std::string>& components,
                          const char* delim) {
  std::string result;
  JoinStrings(components, delim, &result);
  return result;
}

// ----------------------------------------------------------------------
// UnescapeCEscapeSequences()
//    Copies "source" to "dest", rewriting C-style escape sequences
//    -- '\n', '\r', '\\', '\ooo', etc -- to their ASCII
//    equivalents.  "dest" must be sufficiently large to hold all
//    the characters in the rewritten string (i.e. at least as large
//    as strlen(source) + 1 should be safe, since the replacements
//    are always shorter than the original escaped sequences).  It's
//    safe for source and dest to be the same.  RETURNS the length
//    of dest.
//
//    It allows hex sequences \xhh, or generally \xhhhhh with an
//    arbitrary number of hex digits, but all of them together must
//    specify a value of a single byte (e.g. \x0045 is equivalent
//    to \x45, and \x1234 is erroneous).
//
//    It also allows escape sequences of the form \uhhhh (exactly four
//    hex digits, upper or lower case) or \Uhhhhhhhh (exactly eight
//    hex digits, upper or lower case) to specify a Unicode code
//    point. The dest array will contain the UTF8-encoded version of
//    that code-point (e.g., if source contains \u2019, then dest will
//    contain the three bytes 0xE2, 0x80, and 0x99).
//
//    Errors: In the first form of the call, errors are reported with
//    LOG(ERROR). The same is true for the second form of the call if
//    the pointer to the string std::vector is NULL; otherwise, error
//    messages are stored in the std::vector. In either case, the effect on
//    the dest array is not defined, but rest of the source will be
//    processed.
//    ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT int UnescapeCEscapeSequences(const char* source, char* dest);
LIBPROTOBUF_EXPORT int UnescapeCEscapeSequences(const char* source, char* dest,
                                                std::vector<std::string> *errors);

// ----------------------------------------------------------------------
// UnescapeCEscapeString()
//    This does the same thing as UnescapeCEscapeSequences, but creates
//    a new string. The caller does not need to worry about allocating
//    a dest buffer. This should be used for non performance critical
//    tasks such as printing debug messages. It is safe for src and dest
//    to be the same.
//
//    The second call stores its errors in a supplied string std::vector.
//    If the string std::vector pointer is NULL, it reports the errors with LOG().
//
//    In the first and second calls, the length of dest is returned. In the
//    the third call, the new string is returned.
// ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT int UnescapeCEscapeString(const std::string& src, string* dest);
LIBPROTOBUF_EXPORT int UnescapeCEscapeString(const std::string& src, string* dest,
                                             std::vector<std::string> *errors);
LIBPROTOBUF_EXPORT string UnescapeCEscapeString(const std::string& src);

// ----------------------------------------------------------------------
// CEscapeString()
//    Copies 'src' to 'dest', escaping dangerous characters using
//    C-style escape sequences. This is very useful for preparing query
//    flags. 'src' and 'dest' should not overlap.
//    Returns the number of bytes written to 'dest' (not including the \0)
//    or -1 if there was insufficient space.
//
//    Currently only \n, \r, \t, ", ', \ and !isprint() chars are escaped.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int CEscapeString(const char* src, int src_len,
                                     char* dest, int dest_len);

// ----------------------------------------------------------------------
// CEscape()
//    More convenient form of CEscapeString: returns result as a "string".
//    This version is slower than CEscapeString() because it does more
//    allocation.  However, it is much more convenient to use in
//    non-speed-critical code like logging messages etc.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT string CEscape(const std::string& src);

namespace strings {
// Like CEscape() but does not escape bytes with the upper bit set.
LIBPROTOBUF_EXPORT string Utf8SafeCEscape(const std::string& src);

// Like CEscape() but uses hex (\x) escapes instead of octals.
LIBPROTOBUF_EXPORT string CHexEscape(const std::string& src);
}  // namespace strings

// ----------------------------------------------------------------------
// strto32()
// strtou32()
// strto64()
// strtou64()
//    Architecture-neutral plug compatible replacements for strtol() and
//    strtoul().  Long's have different lengths on ILP-32 and LP-64
//    platforms, so using these is safer, from the point of view of
//    overflow behavior, than using the standard libc functions.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int32 strto32_adaptor(const char *nptr, char **endptr,
                                         int base);
LIBPROTOBUF_EXPORT uint32 strtou32_adaptor(const char *nptr, char **endptr,
                                           int base);

inline int32 strto32(const char *nptr, char **endptr, int base) {
  if (sizeof(int32) == sizeof(long))
    return strtol(nptr, endptr, base);
  else
    return strto32_adaptor(nptr, endptr, base);
}

inline uint32 strtou32(const char *nptr, char **endptr, int base) {
  if (sizeof(uint32) == sizeof(unsigned long))
    return strtoul(nptr, endptr, base);
  else
    return strtou32_adaptor(nptr, endptr, base);
}

// For now, long long is 64-bit on all the platforms we care about, so these
// functions can simply pass the call to strto[u]ll.
inline int64 strto64(const char *nptr, char **endptr, int base) {
  //GOOGLE_COMPILE_ASSERT(sizeof(int64) == sizeof(long long),
  //                      sizeof_int64_is_not_sizeof_long_long);
  return strtoll(nptr, endptr, base);
}

inline uint64 strtou64(const char *nptr, char **endptr, int base) {
  //GOOGLE_COMPILE_ASSERT(sizeof(uint64) == sizeof(unsigned long long),
  //                      sizeof_uint64_is_not_sizeof_long_long);
  return strtoull(nptr, endptr, base);
}

// ----------------------------------------------------------------------
// safe_strtob()
// safe_strto32()
// safe_strtou32()
// safe_strto64()
// safe_strtou64()
// safe_strtof()
// safe_strtod()
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT bool safe_strtob(StringPiece str, bool* value);

LIBPROTOBUF_EXPORT bool safe_strto32(const std::string& str, int32* value);
LIBPROTOBUF_EXPORT bool safe_strtou32(const std::string& str, uint32* value);
inline bool safe_strto32(const char* str, int32* value) {
  return safe_strto32(string(str), value);
}
inline bool safe_strto32(StringPiece str, int32* value) {
  return safe_strto32(str.ToString(), value);
}
inline bool safe_strtou32(const char* str, uint32* value) {
  return safe_strtou32(string(str), value);
}
inline bool safe_strtou32(StringPiece str, uint32* value) {
  return safe_strtou32(str.ToString(), value);
}

LIBPROTOBUF_EXPORT bool safe_strto64(const std::string& str, int64* value);
LIBPROTOBUF_EXPORT bool safe_strtou64(const std::string& str, uint64* value);
inline bool safe_strto64(const char* str, int64* value) {
  return safe_strto64(string(str), value);
}
inline bool safe_strto64(StringPiece str, int64* value) {
  return safe_strto64(str.ToString(), value);
}
inline bool safe_strtou64(const char* str, uint64* value) {
  return safe_strtou64(string(str), value);
}
inline bool safe_strtou64(StringPiece str, uint64* value) {
  return safe_strtou64(str.ToString(), value);
}

LIBPROTOBUF_EXPORT bool safe_strtof(const char* str, float* value);
LIBPROTOBUF_EXPORT bool safe_strtod(const char* str, double* value);
inline bool safe_strtof(const std::string& str, float* value) {
  return safe_strtof(str.c_str(), value);
}
inline bool safe_strtod(const std::string& str, double* value) {
  return safe_strtod(str.c_str(), value);
}
inline bool safe_strtof(StringPiece str, float* value) {
  return safe_strtof(str.ToString(), value);
}
inline bool safe_strtod(StringPiece str, double* value) {
  return safe_strtod(str.ToString(), value);
}

// ----------------------------------------------------------------------
// FastIntToBuffer()
// FastHexToBuffer()
// FastHex64ToBuffer()
// FastHex32ToBuffer()
// FastTimeToBuffer()
//    These are intended for speed.  FastIntToBuffer() assumes the
//    integer is non-negative.  FastHexToBuffer() puts output in
//    hex rather than decimal.  FastTimeToBuffer() puts the output
//    into RFC822 format.
//
//    FastHex64ToBuffer() puts a 64-bit unsigned value in hex-format,
//    padded to exactly 16 bytes (plus one byte for '\0')
//
//    FastHex32ToBuffer() puts a 32-bit unsigned value in hex-format,
//    padded to exactly 8 bytes (plus one byte for '\0')
//
//       All functions take the output buffer as an arg.
//    They all return a pointer to the beginning of the output,
//    which may not be the beginning of the input buffer.
// ----------------------------------------------------------------------

// Suggested buffer size for FastToBuffer functions.  Also works with
// DoubleToBuffer() and FloatToBuffer().
static const int kFastToBufferSize = 32;

LIBPROTOBUF_EXPORT char* FastInt32ToBuffer(int32 i, char* buffer);
LIBPROTOBUF_EXPORT char* FastInt64ToBuffer(int64 i, char* buffer);
char* FastUInt32ToBuffer(uint32 i, char* buffer);  // inline below
char* FastUInt64ToBuffer(uint64 i, char* buffer);  // inline below
LIBPROTOBUF_EXPORT char* FastHexToBuffer(int i, char* buffer);
LIBPROTOBUF_EXPORT char* FastHex64ToBuffer(uint64 i, char* buffer);
LIBPROTOBUF_EXPORT char* FastHex32ToBuffer(uint32 i, char* buffer);

// at least 22 bytes long
inline char* FastIntToBuffer(int i, char* buffer) {
  return (sizeof(i) == 4 ?
          FastInt32ToBuffer(i, buffer) : FastInt64ToBuffer(i, buffer));
}
inline char* FastUIntToBuffer(unsigned int i, char* buffer) {
  return (sizeof(i) == 4 ?
          FastUInt32ToBuffer(i, buffer) : FastUInt64ToBuffer(i, buffer));
}
inline char* FastLongToBuffer(long i, char* buffer) {
  return (sizeof(i) == 4 ?
          FastInt32ToBuffer(i, buffer) : FastInt64ToBuffer(i, buffer));
}
inline char* FastULongToBuffer(unsigned long i, char* buffer) {
  return (sizeof(i) == 4 ?
          FastUInt32ToBuffer(i, buffer) : FastUInt64ToBuffer(i, buffer));
}

// ----------------------------------------------------------------------
// FastInt32ToBufferLeft()
// FastUInt32ToBufferLeft()
// FastInt64ToBufferLeft()
// FastUInt64ToBufferLeft()
//
// Like the Fast*ToBuffer() functions above, these are intended for speed.
// Unlike the Fast*ToBuffer() functions, however, these functions write
// their output to the beginning of the buffer (hence the name, as the
// output is left-aligned).  The caller is responsible for ensuring that
// the buffer has enough space to hold the output.
//
// Returns a pointer to the end of the string (i.e. the null character
// terminating the string).
// ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT char* FastInt32ToBufferLeft(int32 i, char* buffer);
LIBPROTOBUF_EXPORT char* FastUInt32ToBufferLeft(uint32 i, char* buffer);
LIBPROTOBUF_EXPORT char* FastInt64ToBufferLeft(int64 i, char* buffer);
LIBPROTOBUF_EXPORT char* FastUInt64ToBufferLeft(uint64 i, char* buffer);

// Just define these in terms of the above.
inline char* FastUInt32ToBuffer(uint32 i, char* buffer) {
  FastUInt32ToBufferLeft(i, buffer);
  return buffer;
}
inline char* FastUInt64ToBuffer(uint64 i, char* buffer) {
  FastUInt64ToBufferLeft(i, buffer);
  return buffer;
}

inline string SimpleBtoa(bool value) {
  return value ? "true" : "false";
}

// ----------------------------------------------------------------------
// SimpleItoa()
//    Description: converts an integer to a string.
//
//    Return value: string
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT string SimpleItoa(int i);
LIBPROTOBUF_EXPORT string SimpleItoa(unsigned int i);
LIBPROTOBUF_EXPORT string SimpleItoa(long i);
LIBPROTOBUF_EXPORT string SimpleItoa(unsigned long i);
LIBPROTOBUF_EXPORT string SimpleItoa(long long i);
LIBPROTOBUF_EXPORT string SimpleItoa(unsigned long long i);

// ----------------------------------------------------------------------
// SimpleDtoa()
// SimpleFtoa()
// DoubleToBuffer()
// FloatToBuffer()
//    Description: converts a double or float to a string which, if
//    passed to NoLocaleStrtod(), will produce the exact same original double
//    (except in case of NaN; all NaNs are considered the same value).
//    We try to keep the string short but it's not guaranteed to be as
//    short as possible.
//
//    DoubleToBuffer() and FloatToBuffer() write the text to the given
//    buffer and return it.  The buffer must be at least
//    kDoubleToBufferSize bytes for doubles and kFloatToBufferSize
//    bytes for floats.  kFastToBufferSize is also guaranteed to be large
//    enough to hold either.
//
//    Return value: string
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT string SimpleDtoa(double value);
LIBPROTOBUF_EXPORT string SimpleFtoa(float value);

LIBPROTOBUF_EXPORT char* DoubleToBuffer(double i, char* buffer);
LIBPROTOBUF_EXPORT char* FloatToBuffer(float i, char* buffer);

// In practice, doubles should never need more than 24 bytes and floats
// should never need more than 14 (including null terminators), but we
// overestimate to be safe.
static const int kDoubleToBufferSize = 32;
static const int kFloatToBufferSize = 24;

namespace strings {

enum PadSpec {
  NO_PAD = 1,
  ZERO_PAD_2,
  ZERO_PAD_3,
  ZERO_PAD_4,
  ZERO_PAD_5,
  ZERO_PAD_6,
  ZERO_PAD_7,
  ZERO_PAD_8,
  ZERO_PAD_9,
  ZERO_PAD_10,
  ZERO_PAD_11,
  ZERO_PAD_12,
  ZERO_PAD_13,
  ZERO_PAD_14,
  ZERO_PAD_15,
  ZERO_PAD_16,
};

struct Hex {
  uint64 value;
  enum PadSpec spec;
  template <class Int>
  explicit Hex(Int v, PadSpec s = NO_PAD)
      : spec(s) {
    // Prevent sign-extension by casting integers to
    // their unsigned counterparts.
#ifdef LANG_CXX11
    static_assert(
        sizeof(v) == 1 || sizeof(v) == 2 || sizeof(v) == 4 || sizeof(v) == 8,
        "Unknown integer type");
#endif
    value = sizeof(v) == 1 ? static_cast<uint8>(v)
          : sizeof(v) == 2 ? static_cast<uint16>(v)
          : sizeof(v) == 4 ? static_cast<uint32>(v)
          : static_cast<uint64>(v);
  }
};

class AlphaNum {
 public:
  // No bool ctor -- bools convert to an integral type.
  // A bool ctor would also convert incoming pointers (bletch).

  AlphaNum(int i32)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastInt32ToBufferLeft(i32, digits_) - &digits_[0]) {}
  AlphaNum(unsigned int u32)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastUInt32ToBufferLeft(u32, digits_) - &digits_[0]) {}
  AlphaNum(long x)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastInt64ToBufferLeft(x, digits_) - &digits_[0]) {}
  AlphaNum(unsigned long x)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastUInt64ToBufferLeft(x, digits_) - &digits_[0]) {}
  AlphaNum(long long int i64)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastInt64ToBufferLeft(i64, digits_) - &digits_[0]) {}
  AlphaNum(unsigned long long int u64)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastUInt64ToBufferLeft(u64, digits_) - &digits_[0]) {}

  AlphaNum(float f)  // NOLINT(runtime/explicit)
      : piece_(digits_, strlen(FloatToBuffer(f, digits_))) {}
  AlphaNum(double f)  // NOLINT(runtime/explicit)
      : piece_(digits_, strlen(DoubleToBuffer(f, digits_))) {}

  AlphaNum(const Eigen::half &f);  // NOLINT(runtime/explicit)
  AlphaNum(Hex hex);  // NOLINT(runtime/explicit)

  AlphaNum(const char *c_str) : piece_(c_str) {}   // NOLINT(runtime/explicit)
  AlphaNum(const StringPiece &pc) : piece_(pc) {}  // NOLINT(runtime/explicit)
  AlphaNum(const glib::string &str)          // NOLINT(runtime/explicit)
      : piece_(str) {}

  StringPiece::size_type size() const { return piece_.size(); }
  const char *data() const { return piece_.data(); }
  StringPiece Piece() const { return piece_; }

 private:
  StringPiece piece_;
  char digits_[kFastToBufferSize];

  // Use ":" not ':'
  AlphaNum(char c);  // NOLINT(runtime/explicit)

  TF_DISALLOW_COPY_AND_ASSIGN(AlphaNum);
};

extern AlphaNum gEmptyAlphaNum;

}  // namespace strings

using strings::AlphaNum;

// ----------------------------------------------------------------------
// StrCat()
//    This merges the given strings or numbers, with no delimiter.  This
//    is designed to be the fastest possible way to construct a string out
//    of a mix of raw C strings, strings, bool values,
//    and numeric values.
//
//    Don't use this for user-visible strings.  The localization process
//    works poorly on strings built up out of fragments.
//
//    For clarity and performance, don't use StrCat when appending to a
//    string.  In particular, avoid using any of these (anti-)patterns:
//      str.append(StrCat(...)
//      str += StrCat(...)
//      str = StrCat(str, ...)
//    where the last is the worse, with the potential to change a loop
//    from a linear time operation with O(1) dynamic allocations into a
//    quadratic time operation with O(n) dynamic allocations.  StrAppend
//    is a better choice than any of the above, subject to the restriction
//    of StrAppend(&str, a, b, c, ...) that none of the a, b, c, ... may
//    be a reference into str.
// ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d,
                                 const AlphaNum& e);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d,
                                 const AlphaNum& e, const AlphaNum& f);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d,
                                 const AlphaNum& e, const AlphaNum& f,
                                 const AlphaNum& g);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d,
                                 const AlphaNum& e, const AlphaNum& f,
                                 const AlphaNum& g, const AlphaNum& h);
LIBPROTOBUF_EXPORT string StrCat(const AlphaNum& a, const AlphaNum& b,
                                 const AlphaNum& c, const AlphaNum& d,
                                 const AlphaNum& e, const AlphaNum& f,
                                 const AlphaNum& g, const AlphaNum& h,
                                 const AlphaNum& i);

inline string StrCat(const AlphaNum& a) { return string(a.data(), a.size()); }

// ----------------------------------------------------------------------
// StrAppend()
//    Same as above, but adds the output to the given string.
//    WARNING: For speed, StrAppend does not try to check each of its input
//    arguments to be sure that they are not a subset of the string being
//    appended to.  That is, while this will work:
//
//    string s = "foo";
//    s += s;
//
//    This will not (necessarily) work:
//
//    string s = "foo";
//    StrAppend(&s, s);
//
//    Note: while StrCat supports appending up to 9 arguments, StrAppend
//    is currently limited to 4.  That's rarely an issue except when
//    automatically transforming StrCat to StrAppend, and can easily be
//    worked around as consecutive calls to StrAppend are quite efficient.
// ----------------------------------------------------------------------

LIBPROTOBUF_EXPORT void StrAppend(string* dest, const AlphaNum& a);
LIBPROTOBUF_EXPORT void StrAppend(string* dest, const AlphaNum& a,
                                  const AlphaNum& b);
LIBPROTOBUF_EXPORT void StrAppend(string* dest, const AlphaNum& a,
                                  const AlphaNum& b, const AlphaNum& c);
LIBPROTOBUF_EXPORT void StrAppend(string* dest, const AlphaNum& a,
                                  const AlphaNum& b, const AlphaNum& c,
                                  const AlphaNum& d);

// ----------------------------------------------------------------------
// Join()
//    These methods concatenate a range of components into a C++ string, using
//    the C-string "delim" as a separator between components.
// ----------------------------------------------------------------------
template <typename Iterator>
void Join(Iterator start, Iterator end,
          const char* delim, string* result) {
  for (Iterator it = start; it != end; ++it) {
    if (it != start) {
      result->append(delim);
    }
    StrAppend(result, *it);
  }
}

template <typename Range>
string Join(const Range& components,
            const char* delim) {
  string result;
  Join(components.begin(), components.end(), delim, &result);
  return result;
}

// ----------------------------------------------------------------------
// ToHex()
//    Return a lower-case hex string representation of the given integer.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT string ToHex(uint64 num);

// ----------------------------------------------------------------------
// GlobalReplaceSubstring()
//    Replaces all instances of a substring in a string.  Does nothing
//    if 'substring' is empty.  Returns the number of replacements.
//
//    NOTE: The string pieces must not overlap s.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int GlobalReplaceSubstring(const std::string& substring,
                                              const std::string& replacement,
                                              string* s);

// ----------------------------------------------------------------------
// Base64Unescape()
//    Converts "src" which is encoded in Base64 to its binary equivalent and
//    writes it to "dest". If src contains invalid characters, dest is cleared
//    and the function returns false. Returns true on success.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT bool Base64Unescape(StringPiece src, string* dest);

// ----------------------------------------------------------------------
// WebSafeBase64Unescape()
//    This is a variation of Base64Unescape which uses '-' instead of '+', and
//    '_' instead of '/'. src is not null terminated, instead specify len. I
//    recommend that slen<szdest, but we honor szdest anyway.
//    RETURNS the length of dest, or -1 if src contains invalid chars.

//    The variation that stores into a string clears the string first, and
//    returns false (with dest empty) if src contains invalid chars; for
//    this version src and dest must be different strings.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int WebSafeBase64Unescape(const char* src, int slen,
                                             char* dest, int szdest);
LIBPROTOBUF_EXPORT bool WebSafeBase64Unescape(StringPiece src, string* dest);

// Return the length to use for the output buffer given to the base64 escape
// routines. Make sure to use the same value for do_padding in both.
// This function may return incorrect results if given input_len values that
// are extremely high, which should happen rarely.
LIBPROTOBUF_EXPORT int CalculateBase64EscapedLen(int input_len,
                                                 bool do_padding);
// Use this version when calling Base64Escape without a do_padding arg.
LIBPROTOBUF_EXPORT int CalculateBase64EscapedLen(int input_len);

// ----------------------------------------------------------------------
// Base64Escape()
// WebSafeBase64Escape()
//    Encode "src" to "dest" using base64 encoding.
//    src is not null terminated, instead specify len.
//    'dest' should have at least CalculateBase64EscapedLen() length.
//    RETURNS the length of dest.
//    The WebSafe variation use '-' instead of '+' and '_' instead of '/'
//    so that we can place the out in the URL or cookies without having
//    to escape them.  It also has an extra parameter "do_padding",
//    which when set to false will prevent padding with "=".
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int Base64Escape(const unsigned char* src, int slen,
                                    char* dest, int szdest);
LIBPROTOBUF_EXPORT int WebSafeBase64Escape(
    const unsigned char* src, int slen, char* dest,
    int szdest, bool do_padding);
// Encode src into dest with padding.
LIBPROTOBUF_EXPORT void Base64Escape(StringPiece src, string* dest);
// Encode src into dest web-safely without padding.
LIBPROTOBUF_EXPORT void WebSafeBase64Escape(StringPiece src, string* dest);
// Encode src into dest web-safely with padding.
LIBPROTOBUF_EXPORT void WebSafeBase64EscapeWithPadding(StringPiece src,
                                                       string* dest);

LIBPROTOBUF_EXPORT void Base64Escape(const unsigned char* src, int szsrc,
                                     string* dest, bool do_padding);
LIBPROTOBUF_EXPORT void WebSafeBase64Escape(const unsigned char* src, int szsrc,
                                            string* dest, bool do_padding);

static const int UTFmax = 4;
// ----------------------------------------------------------------------
// EncodeAsUTF8Char()
//  Helper to append a Unicode code point to a string as UTF8, without bringing
//  in any external dependencies. The output buffer must be as least 4 bytes
//  large.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int EncodeAsUTF8Char(uint32 code_point, char* output);

// ----------------------------------------------------------------------
// UTF8FirstLetterNumBytes()
//   Length of the first UTF-8 character.
// ----------------------------------------------------------------------
LIBPROTOBUF_EXPORT int UTF8FirstLetterNumBytes(const char* src, int len);

}  // namespace google

#endif  // GOOGLE_PROTOBUF_STUBS_STRUTIL_H__
