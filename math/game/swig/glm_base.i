%{
#define GLM_FORCE_SWIZZLE
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/glm.hpp>
#include <opa/math/game/conf.h>
#include <opa/math/game/quat.h>
%}


%{
#include "helper.h"
%}


%define VEC_DEF(type, name, N)
%typemap(in) name(SwigSeqLoader loader) {
  std::vector<type> tb = loader.load<type>($input, N);
  if (!loader.ok) {
    SWIG_fail;
  }
  $1 = vec_to_glm<N,type>(tb);
}


%typemap(in) name& (SwigSeqLoader loader, name tmp) {
  std::vector<type> tb = loader.load<type>($input, N);
  if (!loader.ok) {
    SWIG_fail;
  }
  tmp =  vec_to_glm<N,type>(tb);
  $1 = &tmp;
}

%typemap(out) name (SwigVectorHelper helper), name& (SwigVectorHelper helper) {
  helper.add($1);
  $result = helper.obj;
}

%typemap(in, numinputs=0) name& res (name tmp) { 
  $1=&tmp;
}

%typemap(argout) name &res (SwigVectorHelper helper) {
  %append_output(helper.add(*$1).obj);
}
%enddef

%define MAT_DEF(type, name, N)
%typemap(in) name(SwigSeqLoader loader) {
  std::vector<type> tb = loader.load<type>($input, N);
  if (!loader.ok) {
    SWIG_fail;
  }
  $1 = mat_to_glm<N,type>(tb);
}


%typemap(in) name& (SwigSeqLoader loader, name tmp) {
  std::vector<type> tb = loader.load<type>($input, N);
  if (!loader.ok) {
    SWIG_fail;
  }
  tmp =  mat_to_glm<N,type>(tb);
  $1 = &tmp;
}

%typemap(out) name (SwigVectorHelper helper), name& (SwigVectorHelper helper) {
  helper.add($1);
  $result = helper.obj;
}

%typemap(in, numinputs=0) name& res (name tmp) { 
  $1=&tmp;
}

%typemap(argout) name &res (SwigVectorHelper helper) {
  %append_output(helper.add(*$1).obj);
}
%enddef



%typemap(typecheck) std::string_view {}

%typemap(typecheck, precedence=SWIG_TYPECHECK_SWIGOBJECT) glm::vec3 {$1 = SwigSeqLoader().can_convert<$1_ltype>($input); }
%typemap(typecheck, precedence=SWIG_TYPECHECK_CHAR_PTR) opa::Pos&, const opa::Pos&, glm::vec3&, const glm::vec3& {$1 = SwigSeqLoader().can_convert<$*1_ltype>($input);}

%typemap(in) glm::quat(SwigSeqLoader loader) {
  std::vector<double> tb = loader.load<double>($input, 4);
  if (!loader.ok) {
    SWIG_fail;
  }
  $1 = opa::math::game::quat_from_vec4(vec_to_glm<4,double>(tb));
}


%typemap(in) glm::quat& (SwigSeqLoader loader, glm::quat tmp) {
  std::vector<double> tb = loader.load<double>($input, 4);
  if (!loader.ok) {
    SWIG_fail;
  }
  tmp = opa::math::game::quat_from_vec4(vec_to_glm<4,double>(tb));
  $1 = &tmp;
}

%typemap(out) glm::quat (SwigVectorHelper helper), glm::quat& (SwigVectorHelper helper) {
  helper.add($1);
  $result = helper.obj;
}



VEC_DEF(int, glm::ivec2, 2);
VEC_DEF(double, glm::vec2, 2);
VEC_DEF(int, glm::ivec3, 3);
VEC_DEF(double, glm::vec3, 3);
VEC_DEF(int, glm::ivec4, 4);
VEC_DEF(double, glm::vec4, 4);

MAT_DEF(float, glm::mat2, 2);
MAT_DEF(float, glm::mat3, 3);
MAT_DEF(float, glm::mat4, 4);


%typemap(memberin) glm::vec2, glm::vec2 &, opa::Pos2, opa::Pos2& {
  $1 = $input;
}
%typemap(memberin) glm::vec4, glm::vec4 &, opa::Pos4, opa::Pos4& {
  $1 = $input;
}
%typemap(memberin) glm::vec3, glm::vec3 &, opa::Pos, opa::Pos& {
  $1 = $input;
}

%typemap(memberin) glm::ivec2, glm::ivec2 &, opa::IPos2, opa::IPos2& {
  $1 = $input;
}

%ignore opa::math::game::SphereManifold;
%ignore opa::math::game::RotManifold;

%include "glm/glm.hpp"
%include "opa/math/game/conf.h"


%typemap(typecheck) opa::PointVec {$1 = SwigSeqLoader().can_convert<$1_ltype>($input);}
%typemap(typecheck) opa::Point2Vec {$1 = SwigSeqLoader().can_convert<$1_ltype>($input);}

%typemap(typecheck) std::vector<opa::Pos>&, const std::vector<opa::Pos>&, opa::PointVec&, const opa::PointVec& {$1 = SwigSeqLoader().can_convert<$*1_ltype>($input);}
%typemap(out) std::vector<opa::Pos2> (SwigVectorHelper helper), opa::Point2Vec (SwigVectorHelper helper){
  REP(i, $1.size()){
    helper.add(SwigVectorHelper().add((*(&$1))[i]).obj);
  }
  $result = helper.obj;
}

%typemap(out) std::vector<opa::Pos4> (SwigVectorHelper helper) {
  REP(i, $1.size()){
    helper.add(SwigVectorHelper().add((*(&$1))[i]).obj);
  }
  $result = helper.obj;
}

%typemap(out) std::vector<opa::Pos> (SwigVectorHelper helper) {
  REP(i, $1.size()){
    helper.add(SwigVectorHelper().add((*(&$1))[i]).obj);
  }
  $result = helper.obj;
}

%typemap(out) std::vector<opa::Pos2> (SwigVectorHelper helper) {
  REP(i, $1.size()){
    helper.add(SwigVectorHelper().add((*(&$1))[i]).obj);
  }
  $result = helper.obj;
}


%typemap(in) const opa::PointVec& (SwigSeqLoader loader, std::vector<opa::Pos> tmp), const std::vector<opa::Pos>& (SwigSeqLoader loader, std::vector<opa::Pos> tmp) {
  tmp = loader.load<opa::Pos>($input);
  if (!loader.ok) {
    SWIG_fail;
  }
  $1 = &tmp;
}

%typemap(in) const opa::Point2Vec& (SwigSeqLoader loader, std::vector<opa::Pos2> tmp), const std::vector<opa::Pos2>& (SwigSeqLoader loader, std::vector<opa::Pos2> tmp) {
  tmp = loader.load<opa::Pos2>($input);
  if (!loader.ok) {
    SWIG_fail;
  }
  $1 = &tmp;
}

%typemap(freearg) const opa::PointVec&, const std::vector<opa::Pos>&  { }
%typemap(freearg) const opa::Point2Vec&, const std::vector<opa::Pos2>&  { }

%typemap(in) const opa::Tr2D& (SwigSeqLoader loader, opa::Tr2D tmp) {
  tmp = loader.load_obj<opa::Tr2D>($input);
  $1 = &tmp;
}

%typemap(in) opa::Tr2D, (SwigSeqLoader loader) {
  $1 = loader.load_obj<opa::Tr2D>($input);
}

%apply const opa::Tr2D& { const std::array<glm::vec2, 3>& }
%apply opa::Tr2D { std::array<glm::vec2, 3> }

%typemap(freearg) const opa::Tr2D& (SwigSeqLoader loader, opa::Tr2D tmp) {}

%typemap(out) opa::Tr2D (SwigVectorHelper helper) {
  helper.add($1);
  $result = helper.obj;
}


namespace std {
  %template(vec_pos) std::vector<opa::Pos>;
  %template(vec_vec_pos) std::vector<std::vector<glm::vec3>>;
  %template(vec_pos2) std::vector<opa::Pos2>;
  %template(vec_pos4) std::vector<opa::Pos4>;
  %template(tr2d) std::array<opa::Pos2, 3>;
}

%include "opa/math/game/quat.h"
