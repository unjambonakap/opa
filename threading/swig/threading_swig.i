%module opa_threading_swig

%include "opa.i"

%apply (const u8 *src, int n) { (const u8 *charset, int n) }
%apply (const u8 *src, int n) { (const u8 *data, int n) }
%apply (const u8 *src) { (const u8 *data) }


%{
#include "opa/threading/job.h"
#include "opa/threading/runner.h"
#include "opa/threading/dispatcher.h"

using namespace opa::threading;
using namespace opa::utils;
%}

%include "opa/threading/job.h"
%include "opa/threading/dispatcher.h"
%include "opa/threading/runner.h"
