#pragma once

#include <opa/engine/Camera.h>
#include <opa/engine/Game.h>
#include <opa/engine/InputHandler.h>
#include <opa/engine/conf.h>

class GLFWwindow;

OPA_NAMESPACE_DECL2(opa, engine)
class FrameBufferWrapper;
OPA_DECL_SPTR(FrameBufferWrapper, FrameBufferWrapperPtr)

struct WindowParameters {
  WindowParameters(std::string _title = "", int _w = -1, int _h = -1);
  std::string title;
  int w;
  int h;
};

class Renderer : public opa::utils::Initable {
public:
  Renderer();
  virtual void init(Game *game, const WindowParameters &window_params);
  virtual void fini();
  void init_offscreen();

  void render();
  ~Renderer();
  bool should_quit();

  OPA_ACCESSOR_PTR(GLFWwindow, m_window, window)
  OPA_ACCESSOR_PTR(Game, m_game, game)
  OPA_ACCESSOR_PTR(Scene, m_game->scene(), scene)
  OPA_ACCESSOR_PTR(Camera, m_game->camera(), camera)

  OPA_ACCESSOR_R(u32, m_w, w);
  OPA_ACCESSOR_R(u32, m_h, h);
  OPA_ACCESSOR_R(bool, m_offscreen, offscreen);

  void toggle_wireframe();
  void get_fbo(Image *img) const;

private:
  virtual void init() {}
  bool m_is_wireframe;
  bool m_offscreen;

  FrameBufferWrapperPtr m_fb;
  Game *m_game;
  GLFWwindow *m_window;
  void* m_offscreen_display;
  u32 m_w, m_h;
};
OPA_DECL_SPTR(Renderer, RendererPtr)

OPA_NAMESPACE_DECL2_END
