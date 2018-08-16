#pragma once

#include <opa/engine/conf.h>

#include <GLFW/glfw3.h>
#include <opa/engine/Game.h>
#include <opa/engine/Inputs.h>

OPA_NAMESPACE_DECL2(opa, engine)

class InputHandler;


typedef u64 RunId;

enum InputActionType {
  PRESS,
  RELEASE,
  REPEAT,
};

enum ModifierType { Ctrl, Alt, Shift };

class ActionHelper {
public:
  ActionHelper(InputActionType a) { set(a); }
  ActionHelper() {}
  bool pressed() const { return m_action == InputActionType::PRESS; }
  bool repeat() const { return m_action == InputActionType::REPEAT; }
  bool released() const { return m_action == InputActionType::RELEASE; }

  void set(InputActionType a) { m_action = a; }
  InputActionType get() const { return m_action; }

  void from_glfw(int action) { set(action_type_from_glfw(action)); }
  static InputActionType action_type_from_glfw(int action);

private:
  InputActionType m_action;
};

class ModifierHelper {

public:
  ModifierHelper() { ctrl = shift = alt = 0; }
  void load_from_glfw(int mods);
  void set(ModifierType type, bool val);

  bool ctrl;
  bool shift;
  bool alt;
};

class MouseHelper {
public:
  MouseHelper() { reset(); }
  void reset() { m_internal = 2; }
  void set(const glm::vec2 &pos) {
    m_last = m_cur;
    m_cur = pos;
    m_diff = m_cur - m_last;
    m_internal >>= 1;
  }
  bool has() const { return m_internal == 0; }

  OPA_ACCESSOR_R(glm::vec2, m_last, last)
  OPA_ACCESSOR_R(glm::vec2, m_cur, cur)
  OPA_ACCESSOR_R(glm::vec2, m_diff, diff)

private:
  glm::vec2 m_last;
  glm::vec2 m_cur;
  glm::vec2 m_diff;

  int m_internal;
};

class InputTracker {
public:
  class InputEntry {
  public:
    bool pressed;
    Time time_change;
    RunId last_press;
    void reset() {
      pressed = false;
      last_press = 0;
    }
  };
  void reset();
  bool just_pressed(InputEnum key) const {
    return get(key).last_press == m_last_run;
  }

  InputEntry &get(InputEnum key) { return m_tracks[key]; }
  const InputEntry &get(InputEnum key) const { return m_tracks[key]; }

  void set_last_run(RunId run) { m_last_run = run; }

private:
  RunId m_last_run;

  std::vector<InputEntry> m_tracks;
};

class CallbackRegisterer {
public:
  void add(InputCallerPtr cb) { m_cbs.insert(cb); }
  void remove(InputCallerPtr cb) { m_cbs.erase(cb); }
  OPA_ACCESSOR_R(std::set<InputCallerPtr>, m_cbs, it)
private:
  std::set<InputCallerPtr> m_cbs;
};

class InputHandler : public opa::utils::Initable {
public:
  virtual void init(Game *game);
  virtual void fini();

  InputHandler();

  void do_step();

  OPA_ACCESSOR_PTR(Renderer, m_game->renderer(), renderer)

  OPA_ACCESSOR_PTR(Input, m_input, input);

  OPA_ACCESSOR_R(ActionHelper, m_action, action);
  OPA_ACCESSOR_R(MouseHelper, m_mouse, mouse);
  OPA_ACCESSOR_R(glm::vec2, m_scroll, scroll);
  OPA_ACCESSOR_R(ModifierHelper, m_mods, mods);

  OPA_ACCESSOR_R(InputTracker, m_tracker, tracker);
  OPA_ACCESSOR(RankCaller<InputCallerPtr>, m_cb_caller, cb_caller);

  OPA_ACCESSOR_PTR(CallbackRegisterer, m_cbs, input_cbs);

  void enable_cursor(bool enable);
  void set_cursor(const glm::vec2 &pos);
  void center_cursor();

private:
  static InputHandler *GetInputHandler() {
    return Game::GetInstance()->input_handler();
  }

  static void Glfw_key_cb_entry(GLFWwindow *window, int key, int scancode,
                                int action, int mods) {
    return GetInputHandler()->key_cb(key, scancode, action, mods);
  }
  static void Glfw_mouse_key_cb_entry(GLFWwindow *window, int button,
                                      int action, int mods) {
    return GetInputHandler()->mouse_key_cb(button, action, mods);
  }
  static void Glfw_mouse_cb_entry(GLFWwindow *window, double x, double y) {
    return GetInputHandler()->mouse_cb(x, y);
  }
  static void Glfw_scroll_cb_entry(GLFWwindow *window, double x, double y) {
    return GetInputHandler()->scroll_cb(x, y);
  }
  static void Glfw_focus_cb_entry(GLFWwindow *window, int entered) {
    return GetInputHandler()->focus_cb(entered);
  }

  // called from static
  void key_cb(int key, int scancode, int action, int mods);
  void mouse_key_cb(int button, int action, int mods);
  void mouse_cb(double x, double y);
  void scroll_cb(double x, double y);
  void focus_cb(bool enter);

  void update();
  void handle_cb_for_input();
  void impl_mod_cb();

  Game *m_game;

  Input *m_input;
  MouseHelper m_mouse;
  glm::vec2 m_scroll;

  ActionHelper m_action;
  ModifierHelper m_mods;
  InputTracker m_tracker;
  RankCaller<InputCallerPtr> m_cb_caller;

  CallbackRegisterer m_cbs[InputEnum::ENUM_END];
  RunId m_cur_run;
};
OPA_DECL_SPTR(InputHandler, InputHandlerPtr)

OPA_NAMESPACE_DECL2_END
