#pragma once

#include <experimental/optional>
#include <opa/engine/conf.h>
#include <opa/engine/decl.h>

#include <GLFW/glfw3.h>

OPA_NAMESPACE_DECL2(opa, engine)

typedef u32 ObjId;

template <typename RetVal, typename... Args>
class Updater0 : public opa::utils::Base {
public:
  typedef std::function<RetVal(Args... args)> CallerFunc;

  virtual ~Updater0() {}
  Updater0() {}
  Updater0(CallerFunc func) { this->func = func; }

  template <class U> Updater0(U *v) {
    this->func = [v](Args... args) { return v->swig_call(args...); };
  }

  virtual RetVal update(Args... args) {
    if (func) {
      return (*func)(args...);
    } else {
      return this->_update(args...);
    }
  };

  RetVal operator()(Args... args) { return update(args...); }

  virtual RetVal _update(Args... args) { return RetVal(); } // to override

  std::experimental::optional<CallerFunc> func;
};

class Updater : public Updater0<void> {
public:
  using Updater0<void>::Updater0;
  virtual void _update() override {}
};
OPA_DECL_SPTR(Updater, UpdaterPtr);

class InputCaller : public Updater0<bool, InputHandler *> {
public:
  using Updater0<bool, InputHandler *>::Updater0;

  virtual bool _update(InputHandler *handler) override { return false; }
};
OPA_DECL_SPTR(InputCaller, InputCallerPtr);

class Updater_Swig {
public:
  virtual ~Updater_Swig() {}
  virtual void update();

  void swig_call();

  UpdaterPtr to_updater();
};
OPA_DECL_SPTR(Updater_Swig, Updater_SwigPtr);

class InputCaller_Swig {
public:
  virtual ~InputCaller_Swig() {}
  virtual bool update(InputHandler *handler);
  InputCallerPtr to_input_caller();
  bool swig_call(InputHandler *handler);
};
OPA_DECL_SPTR(InputCaller_Swig, InputCaller_SwigPtr);

enum UpdateMark : int {
  StartLoop = 0,
  InputUpdate = 100,
  CameraUpdate = 150,
  Draw = 200,
  PostDraw = 201,
  EndLoop = 300,
};

enum InputMark : int {
  InputMark_Spec = 0,
  InputMark_Generaal = 100,
};

template <class T> class RankCaller {
public:
  typedef int RankId;
  void add(RankId id, T func) { m_funcs[id].insert(func); }
  void remove(RankId id, T func) { m_funcs[id].erase(func); }

  template <class U> void do_loop0() {
    for (const auto &lst : m_funcs)
      for (const auto &func : lst.ND) (*func)();
  }

  template <class U> bool do_loop(U x) {
    for (const auto &lst : m_funcs)
      for (const auto &func : lst.ND)
        if ((*func)(x)) return true;
    return false;
  }

  template <class U> bool do_loop_one(int rank, U x) {
    for (const auto &func : glib::gtl::FindPtrOrNull(m_funcs, x))
      if ((*func)(x)) return true;
    return false;
  }

  std::map<RankId, std::set<T> > m_funcs;
};

using UpdateManager = RankCaller<UpdaterPtr>;

template <typename T> struct TypeName {
  static std::string Get() { return typeid(T).name(); }
};

class ObjectStore {
public:
  template <class T> void add_singleton(SPTR(T) ptr) {
    m_singletons[TypeName<T>::Get()] =
      std::static_pointer_cast<opa::utils::Base, T>(ptr);
  }
  template <class T> SPTR(T) get() const {
    std::string name = TypeName<T>::Get();
    OPA_CHECK0(m_singletons.count(name));
    return std::static_pointer_cast<T, opa::utils::Base>(
      m_singletons.find(name)->ND);
  }

private:
  std::unordered_map<std::string, opa::utils::BasePtr> m_singletons;
};

class Game : public opa::utils::Initable {
public:
  Game();
  virtual void init();

  void run();
  static Game *GetInstance();

  OPA_ACCESSOR_PTR_DECL(Game, InputHandler, m_input_handler.get(),
                        input_handler);
  OPA_ACCESSOR_PTR_DECL(Game, Camera, m_camera.get(), camera);
  OPA_ACCESSOR_PTR_DECL(Game, Renderer, m_renderer.get(), renderer);
  OPA_ACCESSOR_PTR_DECL(Game, Scene, m_scene.get(), scene);
  OPA_ACCESSOR_PTR_DECL(Game, SceneTreeHelper, m_tree.get(), tree);
  OPA_ACCESSOR_PTR_DECL(Game, DataCollector, m_data_collector.get(), data_collector);
  OPA_ACCESSOR_R(ObjectStore, m_store, store);

  OPA_ACCESSOR_R(Time, m_loop_time, time)

  void register_obj(SceneObjPtr obj, SceneObj *par);
  void remove_obj(SceneObj *obj);
  void check_register(SceneObj *obj);
  SceneObj *get_obj(ObjId id);
  OPA_ACCESSOR(UpdateManager, m_updater, updater);
  OPA_ACCESSOR_R(ColorPalette, m_palette, palette);

private:
  InputHandlerPtr m_input_handler;
  CameraPtr m_camera;
  RendererPtr m_renderer;
  ScenePtr m_scene;
  ServiceManagerPtr m_services;
  SceneTreeHelperPtr m_tree;
  UpdateManager m_updater;
  Time m_loop_time;
  std::unordered_map<ObjId, SceneObjPtr> m_obj_store;
  RenderStreamerPtr m_render_streamer;
  DataCollectorPtr m_data_collector;

  ObjectStore m_store;
  ColorPalette m_palette;
  ObjId m_cur;
};
OPA_DECL_SPTR(Game, GamePtr)

OPA_NAMESPACE_DECL2_END
