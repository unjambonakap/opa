%include "stdint.i"
%include "std_string.i"
%include "std_array.i"
%include "typemaps.i"
%include "std_vector.i"
%include "std_set.i"
%include "std_map.i"
%include "std_shared_ptr.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%{
#include <stdint.h>
#include <opa_inc.h>
#include <opa_callback.h>
#include <opa_common_base.h>
#include <opa_common.h>
#include <glib/core/stringpiece.h>
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply (float* IN_ARRAY1, int DIM1) {(const float* data, int n)};
%apply (double* IN_ARRAY1, int DIM1) {(const double* data, int n)};

%feature("director") opa::OpaCallback;

%include "opa_callback.h"
%include "opa_inc.h"
%include "opa_common_base.h"
%include "opa_common.h"
#include "glib/core/stringpiece.h"
%ignore DECLARE_;

namespace std {
  %template(si) std::set<int>;
  %template(mii) std::map<int, int>;
  %template(vi) std::vector<int>;
  %template(vvi) std::vector<std::vector<int>>;
  %template(vd) std::vector<double>;
  %template(v_u32) std::vector<u32>;
  %template(vf) std::vector<float>;
  %template(vs) std::vector<std::string>;
  %template(vvd) std::vector<std::vector<double>>;
  %template(vvf) std::vector<std::vector<float>>;
}


%typemap(in) u16 {
  $1=PyLong_AsUnsignedLong($input);
}

%typemap(in) u64 {
  $1=PyLong_AsUnsignedLongLong($input);
}
%{
#include "swig_common.h"
%}

%typemap(in) glib::StringPiece (swig_tsf_data data){
  $1=data.get_piece($input);
  if (!data.ok) SWIG_fail;
}
%typemap(in) const glib::StringPiece& (swig_tsf_data data){
  $1=data.get_piece_ptr($input);
  if (!data.ok) SWIG_fail;
}
%typemap(in) const opa::stolen::StringRef& (swig_tsf_data data){
  $1=data.get_ref_ptr($input);
  if (!data.ok) SWIG_fail;
}
%typemap(in) opa::stolen::StringRef (swig_tsf_data data){
  $1=data.get_ref($input);
  if (!data.ok) SWIG_fail;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_SWIGOBJECT) glib::StringPiece {$1 = 1; }
%typemap(typecheck, precedence=SWIG_TYPECHECK_SWIGOBJECT) opa::stolen::StringRef{$1 = 1; }



%typemap(out) std::string {
  $result=swig_helper::convert($1);
}
%typemap(out) std::pair<std::string, bool> {
  $result=swig_helper::convert($1);
}

%typemap(in) const std::string& (swig_tsf_data data){
  $1=data.get_str($input);
  if (!data.ok) SWIG_fail;
}
%typemap(freearg) const std::string& { 
}

%typemap(in, numinputs=1) (const u8 *src, u8 *dest, int n) (swig_tsf_data data) { 
  $1=data.get_u8($input);
  $3=data.n;
  if (!data.ok) SWIG_fail;
  $2=(u8*)malloc($3+1);
}

%typemap(argout) (const u8 *src, u8 *dest, int n) {
  %append_output(PyBytes_FromStringAndSize((const char*)$2, $3)); 
}

%typemap(freearg) (const u8 *src, u8 *dest, int n) {
  free($2);
}

%typemap(in, numinputs=1) (const u8 *src, int n) (swig_tsf_data data){
  $1=data.get_u8($input);
  $2=data.n;
  if (!data.ok) SWIG_fail;
}
%typemap(freearg) (const u8 *src, int n) {}

%typemap(in, numinputs=1) (const u8 *src) (swig_tsf_data data){
  $1=data.get_u8($input);
  if (!data.ok) SWIG_fail;
}
%typemap(freearg) (const u8 *src) {}


%typemap(in, numinputs=0) (u64 &res) (u64 tmp) {
  $1=&tmp;
}

%typemap(argout) (u64 &res) {
  %append_output(PyLong_FromUnsignedLongLong(*$1));
}

%typemap(in, numinputs=0) (std::string &out) (swig_tsf_data data) { 
  $1 = new std::string;
}
%typemap(argout) (std::string &out) {
  %append_output(PyBytes_FromStringAndSize((const char*)$1->data(), $1->size())); 
}

%typemap(freearg) (std::string &out) {
  delete $1;
}

%apply const std::string& {std::string* };
%apply const std::vector<int>& {std::vector<int>* };
%apply const std::vector<double>& {std::vector<double>* };
%apply const std::vector<std::string>& {std::vector<std::string>* };

