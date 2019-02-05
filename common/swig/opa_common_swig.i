%module opa_common_swig


%include "opa.i"


%{
#include "opa/common_swig.h"

using namespace opa::init;
#include "opa/utils/misc.h"
#include "opa/utils/DataStruct.h"
%}


%shared_ptr(opa::utils::Base);
%shared_ptr(opa::utils::Initable);

%include "opa/common_swig.h"
%include "opa/utils/misc.h"
%include "opa/utils/DataStruct.h"
