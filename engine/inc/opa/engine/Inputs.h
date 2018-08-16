#pragma once

#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

enum InputSource {
    Key,
    Mouse,
    Button,
    Joystick,
    Focus,
};

#if !OPA_SWIG
#define OPA_INPUT_DECL(name, glfw_val, source) name,
enum InputEnum {
#include <opa/engine/Inputs_vals.h>
    ENUM_END,
};
#undef OPA_INPUT_DECL
#endif

class Input {

  public:
    Input();
    Input(InputEnum type, int glfw_val, const std::string &name,
          InputSource source);
    static Input *from_type(InputEnum type);
    static Input *from_glfw_key(int key);
    static Input *from_glfw_button(int key);
    static Input *from_string(const std::string &str);

    OPA_ACCESSOR_R(std::string, m_name, name)
    OPA_ACCESSOR_R(InputEnum, m_type, type)
    OPA_ACCESSOR_R(InputSource, m_source, source)
  private:
    InputEnum m_type;

    InputSource m_source;
    int m_glfw_val;
    std::string m_name;
};

OPA_NAMESPACE_DECL2_END
