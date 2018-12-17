%module opa_algo_swig

%include "opa.i"


%{
#include "opa/common_swig.h"
#include "opa/algo/graph.h"
#include "opa/algo/base.h"
#include "opa/algo/strings.h"

using namespace opa::init;
%}
namespace opa{
namespace algo{
}
}

%shared_ptr(opa::algo::FastGraph);
%ignore opa::algo::graph::rmp;
%include "opa/common_swig.h"
%include "opa/algo/graph.h"
%include "opa/algo/base.h"
%include "opa/algo/strings.h"

namespace std {
  %template(vec_graph) std::vector<std::shared_ptr<opa::algo::FastGraph>>;
}
