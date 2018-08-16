#include "Renderer.h"
#include "Camera.h"
#include "Scene.h"
#include <EGL/egl.h>
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glxew.h>
#include <GLFW/glfw3.h>

DEFINE_int32(width, 800, "");
DEFINE_int32(height, 600, "");
DEFINE_string(title, "renderer", "");
DEFINE_bool(render, true, "");
DEFINE_bool(render_to_texture, false, "");
DEFINE_bool(offscreen, false, "");
DEFINE_int32(render_nfbo, 3, "");
DEFINE_int32(opengl_major, 4, "");
DEFINE_int32(opengl_minor, 5, "");

DEFINE_bool(window, true, "");

using namespace opa::math::game;
OPA_NAMESPACE_DECL2(opa, engine)

class FrameBufferWrapper {

public:
  void init(int w, int h, int n) {

    this->n = n;
    fbs.resize(n);
    texs.resize(n);
    depth_bufs.resize(n);

    glGenFramebuffers(n, fbs.data());
    glGenTextures(n, texs.data());
    glGenRenderbuffers(n, depth_bufs.data());

    REP (i, n) {
      glBindFramebuffer(GL_FRAMEBUFFER, fbs[i]);

      // "Bind" the newly created texture : all future texture functions will
      // modify this texture
      glBindTexture(GL_TEXTURE_2D, texs[i]);

      // Give an empty image to OpenGL ( the last "0" )
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE,
                   0);

      // Poor filtering. Needed !
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      glBindRenderbuffer(GL_RENDERBUFFER, depth_bufs[i]);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, 1024, 768);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                GL_RENDERBUFFER, depth_bufs[i]);

      glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texs[i], 0);

      // Set the list of draw buffers.
      GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
      glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

      OPA_CHECK0(glCheckFramebufferStatus(GL_FRAMEBUFFER) ==
                 GL_FRAMEBUFFER_COMPLETE);
    }
    cur = 0;
  }

  void fini() {}

  GLuint get_fb() {
    cur = (cur + 1) % n;
    return fbs[cur];
  }

  int n;

  int cur;
  std::vector<GLuint> fbs;
  std::vector<GLuint> texs;
  std::vector<GLuint> depth_bufs;
};

void error_cb(int error, const char *desc) {
  fprintf(stderr, "Glfw error code=%d, desc=%s\n", error, desc);
}
void GlobalInit() {
  if (FLAGS_render && !FLAGS_offscreen) {
    OPA_CHECK0(glfwInit() == GL_TRUE);
    glfwSetErrorCallback(error_cb);
  }
}

WindowParameters::WindowParameters(std::string _title, int _w, int _h) {
  title = _title;
  w = _w;
  h = _h;

  if (w == -1) w = FLAGS_width;
  if (h == -1) h = FLAGS_height;
  if (title.size() == 0) title = FLAGS_title;
}

Renderer::Renderer() {
  m_window = 0;
  m_is_wireframe = false;
}

Renderer::~Renderer() {}
void Renderer::toggle_wireframe() {
  if (!m_is_wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  m_is_wireframe ^= 1;
}

void Renderer::init_offscreen() {
  const EGLint configAttribs[] = { EGL_SURFACE_TYPE,
                                   EGL_PBUFFER_BIT,
                                   EGL_BLUE_SIZE,
                                   8,
                                   EGL_GREEN_SIZE,
                                   8,
                                   EGL_RED_SIZE,
                                   8,
                                   EGL_DEPTH_SIZE,
                                   8,
                                   EGL_CONTEXT_MAJOR_VERSION,
                                   FLAGS_opengl_major,
                                   EGL_CONTEXT_MINOR_VERSION,
                                   FLAGS_opengl_minor,
                                   EGL_RENDERABLE_TYPE,
                                   EGL_OPENGL_BIT,
                                   EGL_NONE };

  static const EGLint pbufferAttribs[] = {
    EGL_WIDTH, (EGLint)m_w, EGL_HEIGHT, (EGLint)m_h, EGL_NONE,
  };

  // 1. Initialize EGL
  m_offscreen_display = eglGetDisplay(EGL_DEFAULT_DISPLAY);

  EGLint major, minor;

  OPA_CHECK0(eglInitialize(m_offscreen_display, &major, &minor));
  OPA_DISP("On ", major, minor);

  // 2. Select an appropriate configuration
  EGLint numConfigs;
  EGLConfig eglCfg;

  eglChooseConfig(m_offscreen_display, configAttribs, &eglCfg, 1, &numConfigs);

  // 3. Create a surface
  EGLSurface eglSurf =
    eglCreatePbufferSurface(m_offscreen_display, eglCfg, pbufferAttribs);

  // 4. Bind the API
  eglBindAPI(EGL_OPENGL_API);

  // 5. Create a context and make it current
  EGLContext eglCtx =
    eglCreateContext(m_offscreen_display, eglCfg, EGL_NO_CONTEXT, NULL);

  eglMakeCurrent(m_offscreen_display, eglSurf, eglSurf, eglCtx);
}

void Renderer::init(Game *game, const WindowParameters &window_params) {
  opa::utils::Initable::init();
  m_game = game;
  m_w = window_params.w;
  m_h = window_params.h;
  m_offscreen = FLAGS_offscreen;

  if (FLAGS_offscreen) {
    init_offscreen();

  } else {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, FLAGS_opengl_major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, FLAGS_opengl_minor);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if (!FLAGS_window) {
      glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    }

    m_window = glfwCreateWindow(window_params.w, window_params.h,
                                window_params.title.c_str(), NULL, NULL);
    OPA_CHECK(m_window, "failed to create window");

    glfwMakeContextCurrent(m_window);
    glfwSetWindowPos(m_window, 0, 0);
    OPA_CHECK0(glfwExtensionSupported("GL_ARB_pixel_buffer_object"));
  }

  camera()->update_viewport(m_w, m_h);
  camera()->setup();

  GLenum status = glewInit();
  OPA_CHECK(status == GLEW_OK, "Error: %s\n", glewGetErrorString(status));
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  GL_CHECK(;);

  if (0) {

    GLint n;
    glGetIntegerv(GL_NUM_EXTENSIONS, &n);

    OPA_DISP0(n);
    for (GLint i = 0; i < n; i++) {
      const char *extension = (const char *)glGetStringi(GL_EXTENSIONS, i);
      OPA_DISP0("Extension #%d: %s", i, extension);
    }
  }

  if (FLAGS_render_to_texture) {
    m_fb.reset(new FrameBufferWrapper);
    m_fb->init(m_w, m_h, FLAGS_render_nfbo);
  }
}
void Renderer::fini() {
  if (FLAGS_offscreen) {
    eglTerminate(m_offscreen_display);

  } else {
    if (m_window) {
      glfwDestroyWindow(m_window);
      m_window = 0;
    }
  }
  opa::utils::Initable::fini();
}

void Renderer::render() {
  check_init();

  int cur_w, cur_h;

  glfwGetFramebufferSize(m_window, &cur_w, &cur_h);

  if (m_fb) {
    glBindFramebuffer(GL_FRAMEBUFFER, m_fb->get_fb());
  } else {
    glDrawBuffer(GL_BACK);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  glViewport(0, 0, m_w, m_h);

  glClearColor(0.0, 0.0, 0.5, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  scene()->draw();

  if (FLAGS_render_to_texture) {
    glFinish();
  } else {
    glfwSwapBuffers(m_window);
  }

  if (m_w != cur_w || m_h != cur_h) camera()->update_viewport(cur_w, cur_h);

  if (m_fb) {
    //glBindFramebuffer(GL_FRAMEBUFFER, m_fb->fb);
  } else {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glReadBuffer(GL_FRONT);
  }
  GL_CHECK(;);

  m_w = cur_w;
  m_h = cur_h;
}

void Renderer::get_fbo(Image *img) const {
  glReadPixels(0, 0, img->w(), img->h(), img->get_gl_enum(), GL_UNSIGNED_BYTE,
               img->buf());
}

bool Renderer::should_quit() { return glfwWindowShouldClose(m_window); }

OPA_REGISTER_INIT(renderer_init, GlobalInit);

OPA_NAMESPACE_DECL2_END
