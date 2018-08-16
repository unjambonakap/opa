#include "Program.h"

using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)
void Program::init() {
  opa::utils::Initable::init();
  m_id = glCreateProgram();
}

void Program::setup() {

  GLint result;
  glLinkProgram(m_id);
  glGetProgramiv(m_id, GL_LINK_STATUS, &result);
  OPA_CHECK(result == GL_TRUE, "Failed to link");

  bind();
  for (auto attr : m_attrs) attr->init(this);
}

void Program::register_attribute(Attribute *attr) { m_attrs.pb(attr); }

void Program::bind() { GL_CHECK(glUseProgram(m_id)); }

void Program::fini() {
  m_attrs.clear();
  glDeleteProgram(m_id);
  opa::utils::Initable::fini();
}

void Program::add_shader(GLenum shader_type, const char *content) {
  check_init();

  GLint res;
  GLuint shader = glCreateShader(shader_type);
  GL_CHECK(glShaderSource(shader, 1, (const char **)&content, NULL));
  GL_CHECK(glCompileShader(shader));

  glGetShaderiv(shader, GL_COMPILE_STATUS, &res);
  if (res == GL_FALSE) {
    int len;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
    fprintf(stderr, "error in vertex shader %s\n", content);
    GLboolean supported;

    glGetBooleanv(GL_SHADER_COMPILER, &supported);
    printf("supported >> %d\n", supported);
    if (len > 1) {
      char *tmp = (char *)malloc(len);
      glGetInfoLogARB(shader, len, &len, tmp);
      cerr << tmp << endl;
      free(tmp);
    }
    glDeleteShader(shader);
    OPA_CHECK0(false);
  }

  GL_CHECK(glAttachShader(m_id, shader));
}

OPA_NAMESPACE_DECL2_END
