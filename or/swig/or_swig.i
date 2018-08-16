%module opa_or_swig

%include "opa.i"


%{
#define SWIG_FILE_WITH_INIT
#include "opa/or/grid_search.h"
#include "opa/or/adaptative_search.h"

using namespace opa::OR;
using namespace opa::utils;
%}

%include "numpy.i"
%init %{
import_array();
%}

%rename(DspBinaryMatcherV) opa::OR::DspBinaryMatcher::DspBinaryMatcher(const std::vector<float> &data);

%apply (float* IN_ARRAY1, int DIM1) {(const float* data, int n)};
%apply (double* IN_ARRAY1, int DIM1) {(const double* data, int n)};
%include "opa/or/or_common.h"
%include "opa/or/grid_search.h"
%include "opa/or/adaptative_search.h"

