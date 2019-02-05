%module opa_or_swig

%include "opa.i"


%{
#include "opa/or/grid_search.h"
#include "opa/or/adaptative_search.h"

using namespace opa::OR;
using namespace opa::utils;
%}


%rename(DspBinaryMatcherV) opa::OR::DspBinaryMatcher::DspBinaryMatcher(const std::vector<float> &data);

%include "opa/or/or_common.h"
%include "opa/or/grid_search.h"
%include "opa/or/adaptative_search.h"

