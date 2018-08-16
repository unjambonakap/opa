%module(directors="1") opa_wrapper_swig

%include "std_string.i"
%include "typemaps.i"
%include "opa.i"


%{
extern "C" {
#include "opa/wrapper/curve25519-donna.h"
}
%}

%typemap(in, numinputs=0) (u8 *mypublic) (swig_tsf_data data) { 
  $1=(u8*)malloc(32);
}

%typemap(argout) (u8 *mypublic) {
  %append_output(PyBytes_FromStringAndSize((const char*)$1, 32)); 
}

%typemap(freearg) (u8 *mypublic) {
  free($1);
}

%typemap(in, numinputs=1) (const u8*) (swig_tsf_data data){
  $1=data.get_u8($input);
  if (!data.ok) SWIG_fail;
}
%typemap(freearg) (const u8 *) {}


extern "C" {
%include "opa/wrapper/curve25519-donna.h"
}
