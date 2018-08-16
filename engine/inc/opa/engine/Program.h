#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class Program : public opa::utils::Initable {
public:
  Program() {}

  virtual void init();
  virtual void fini();

  void add_shader(GLenum shader_type, const char *content);
  void setup();
  void register_attribute(Attribute *attr);
  void bind();

  OPA_ACCESSOR_R(GLuint, m_id, id)

private:
  std::vector<Attribute *> m_attrs;
  GLuint m_id;
};

OPA_NAMESPACE_DECL2_END
