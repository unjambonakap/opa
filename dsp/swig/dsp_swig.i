%module opa_dsp_swig

%include "math_common_swig_base.i"


%{
#include "opa/utils/serialize.h"

#include "opa/dsp/gps/l1ca.h"

using namespace opa::dsp;
using namespace opa::dsp::gps;
using namespace opa::utils;
%}


%include "opa/utils/serialize.h"
%include "opa/dsp/gps/l1ca.h"


