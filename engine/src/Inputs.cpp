#include "Inputs.h"

OPA_NAMESPACE_DECL2(opa, engine)

static Input s_input_data[InputEnum::ENUM_END];
static std::map<int, Input *> s_key_map;
static std::map<int, Input *> s_button_map;
static std::map<int, Input *> s_joy_map;

Input::Input() {}

Input::Input(InputEnum type, int glfw_val, const std::string &name, InputSource source) {
    m_type = type;
    m_glfw_val = glfw_val;
    m_name = name;
    m_source = source;
}

void handleNewSource(InputEnum type, int glfw_val, const std::string &name, InputSource source) {
    switch (source) {
        OPA_CASE(InputSource::Key, s_key_map[glfw_val] = &s_input_data[type];)
        OPA_CASE(InputSource::Button, s_button_map[glfw_val] = &s_input_data[type];)
        OPA_CASE(InputSource::Joystick, s_joy_map[glfw_val] = &s_input_data[type];)
    default:
        break;
    }

    s_input_data[type] = Input(type, glfw_val, name, source);
}

static void InputInit() {
#define OPA_INPUT_DECL(name, glfw_val, source)                                                     \
    handleNewSource(InputEnum::name, glfw_val, #name, InputSource::source);
#include "Inputs_vals.h"
#undef OPA_INPUT_DECL
}
OPA_REGISTER_INIT(input_init, InputInit);

Input *Input::from_type(InputEnum type) { return &s_input_data[type]; }
Input *Input::from_glfw_key(int key) { return s_key_map[key]; }
Input *Input::from_glfw_button(int key) { return s_button_map[key]; }
Input *Input::from_string(const std::string &str) { return nullptr; }

OPA_NAMESPACE_DECL2_END
