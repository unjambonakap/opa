#include "InputHandler.h"
#include "Renderer.h"

OPA_NAMESPACE_DECL2(opa, engine)

InputHandler::InputHandler() { m_game = 0; }

static std::vector<InputActionType> s_glfw_to_action;
void actions_init() {
  s_glfw_to_action.resize(3);
  s_glfw_to_action[GLFW_PRESS] = InputActionType::PRESS;
  s_glfw_to_action[GLFW_RELEASE] = InputActionType::RELEASE;
  s_glfw_to_action[GLFW_REPEAT] = InputActionType::REPEAT;
}
OPA_REGISTER_INIT(ActionsInit, actions_init)

InputActionType ActionHelper::action_type_from_glfw(int action) {
  return s_glfw_to_action[action];
}

void InputTracker::reset() {
  m_tracks.resize(InputEnum::ENUM_END);
  for (auto &x : m_tracks) x.reset();
}

void ModifierHelper::set(ModifierType type, bool val) {
  switch (type) {
    OPA_CASE(ModifierType::Shift, shift = val;)
    OPA_CASE(ModifierType::Ctrl, ctrl = val;)
    OPA_CASE(ModifierType::Alt, alt = val;)
  default:
    OPA_CHECK0(false);
  }
}

void ModifierHelper::load_from_glfw(int mods) {
  shift = bool(mods & GLFW_MOD_SHIFT);
  ctrl = bool(mods & GLFW_MOD_CONTROL);
  alt = bool(mods & GLFW_MOD_ALT);
}

void InputHandler::init(Game *game) {
  opa::utils::Initable::init();
  m_game = game;
  renderer()->check_init();

  m_tracker.reset();
  m_cur_run = 0;
  if (!renderer()->offscreen()) {

    glfwSetKeyCallback(renderer()->window(), InputHandler::Glfw_key_cb_entry);
    glfwSetMouseButtonCallback(renderer()->window(),
                               InputHandler::Glfw_mouse_key_cb_entry);
    glfwSetCursorPosCallback(renderer()->window(),
                             InputHandler::Glfw_mouse_cb_entry);
    glfwSetScrollCallback(renderer()->window(),
                          InputHandler::Glfw_scroll_cb_entry);
    glfwSetCursorEnterCallback(renderer()->window(),
                               InputHandler::Glfw_focus_cb_entry);
  }

  InputCallerPtr i_mod_cb =
    std::make_shared<InputCaller>([this](InputHandler *) {
      this->impl_mod_cb();
      return true;
    });

  FE_LIST(InputEnum, x, input_cbs()[x].add(i_mod_cb),
          { InputEnum::Key_LEFT_CONTROL, InputEnum::Key_RIGHT_CONTROL,
            InputEnum::Key_LEFT_SHIFT, InputEnum::Key_RIGHT_SHIFT,
            InputEnum::Key_LEFT_ALT });

  this->m_cb_caller.add(InputMark::InputMark_Spec,
                        std::make_shared<InputCaller>([&](InputHandler *) {
                          bool found = false;
                          for (auto cb : m_cbs[m_input->type()].it()) {
                            (*cb)(this);
                            found = true;
                          }
                          return found;
                        }));
}

void InputHandler::enable_cursor(bool enable) {
  if (!renderer()->offscreen()) {
    glfwSetInputMode(renderer()->window(), GLFW_CURSOR,
                     enable ? GLFW_CURSOR_NORMAL : GLFW_CURSOR_DISABLED);
  }
}

void InputHandler::set_cursor(const glm::vec2 &pos) {
  if (!renderer()->offscreen()) {
    m_mouse.set(glm::vec2(pos));
    glfwSetCursorPos(renderer()->window(), pos.x, pos.y);
  }
}

void InputHandler::fini() {}

void InputHandler::do_step() {
  ++m_cur_run;
  m_tracker.set_last_run(m_cur_run);

  if (!renderer()->offscreen()) {
    m_mouse.set(m_mouse.cur());
    m_scroll = glm::vec2();
    glfwPollEvents();

    double xpos, ypos;
    glfwGetCursorPos(renderer()->window(), &xpos, &ypos);
    m_mouse.set(glm::vec2(xpos, ypos));
  }
}
void InputHandler::update() {}

void InputHandler::handle_cb_for_input() {
  if (m_action.pressed()) {
    m_tracker.get(m_input->type()).pressed = true;
    m_tracker.get(m_input->type()).last_press = m_cur_run;
  } else if (m_action.released())
    m_tracker.get(m_input->type()).pressed = false;
  else if (m_action.repeat())
    m_tracker.get(m_input->type()).pressed = false;

  this->m_cb_caller.do_loop(this);
}

void InputHandler::key_cb(int key, int scancode, int action, int mods) {
  m_input = Input::from_glfw_key(key);
  m_action.from_glfw(action);
  m_mods.load_from_glfw(mods);

  handle_cb_for_input();
}

void InputHandler::mouse_key_cb(int button, int action, int mods) {
  m_input = Input::from_glfw_button(button);
  m_action.from_glfw(action);
  m_mods.load_from_glfw(mods);
  handle_cb_for_input();
}

void InputHandler::mouse_cb(double x, double y) {
  m_input = Input::from_type(InputEnum::Mouse_Mov);
  handle_cb_for_input();
}

void InputHandler::scroll_cb(double x, double y) {
  m_input = Input::from_type(InputEnum::Mouse_Scroll);
  m_scroll.x = x;
  m_scroll.y = y;
  handle_cb_for_input();
}

void InputHandler::focus_cb(bool enter) {
  m_input = Input::from_type(InputEnum::Focus_Enter);
  m_action.set(enter ? InputActionType::PRESS : InputActionType::RELEASE);
  m_mouse.reset();

  double xpos, ypos;
  handle_cb_for_input();
}

void InputHandler::center_cursor() {
  set_cursor(glm::vec2(renderer()->w() / 2., renderer()->h() / 2.));
}

void InputHandler::impl_mod_cb() {
  switch (m_input->type()) {
    OPA_CASE(Key_LEFT_CONTROL,
             m_mods.set(ModifierType::Ctrl, m_action.pressed());)
    OPA_CASE(Key_RIGHT_CONTROL,
             m_mods.set(ModifierType::Ctrl, m_action.pressed());)
    OPA_CASE(Key_LEFT_SHIFT,
             m_mods.set(ModifierType::Shift, m_action.pressed());)
    OPA_CASE(Key_RIGHT_SHIFT,
             m_mods.set(ModifierType::Shift, m_action.pressed());)
    OPA_CASE(Key_LEFT_ALT, m_mods.set(ModifierType::Alt, m_action.pressed());)
  default:
    break;
  }
}

OPA_NAMESPACE_DECL2_END
