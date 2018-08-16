#include "Game.h"

#include "Camera.h"
#include "CameraObserver.h"
#include "InputHandler.h"
#include "PosUpdater.h"
#include "RenderStreamer.h"
#include "Renderer.h"
#include "RotatorService.h"
#include "Scene.h"
#include "SceneTreeHelper.h"
#include "Service.h"
#include "Sphere.h"
#include "data_collector.h"
#include "graph/graph.h"
#include "primitives/Base.h"
#include "procgen/Texture.h"

#include "mesh/MeshDrawer.h"
#include "mesh/MeshIsec.h"
#include "mesh/mesh.h"
#include <opa/math/common/Utils.h>
#include <opa/math/game/mesh.h>
#include <opa/math/game/mesh_util.h>

DEFINE_bool(game_default_objs, false, "");
using namespace opa::math::game;
using namespace std;

OPA_NAMESPACE_DECL2(opa, engine)

OPA_ACCESSOR_PTR_IMPL(Game, InputHandler, m_input_handler.get(), input_handler)
OPA_ACCESSOR_PTR_IMPL(Game, Camera, m_camera.get(), camera)
OPA_ACCESSOR_PTR_IMPL(Game, Renderer, m_renderer.get(), renderer)
OPA_ACCESSOR_PTR_IMPL(Game, Scene, m_scene.get(), scene)
OPA_ACCESSOR_PTR_IMPL(Game, SceneTreeHelper, m_tree.get(), tree)
OPA_ACCESSOR_PTR_IMPL(Game, DataCollector, m_data_collector.get(), data_collector)

Game *s_game_instance = 0;
Game *Game::GetInstance() { return s_game_instance; }

UpdaterPtr Updater_Swig::to_updater() {
  OPA_DISP("Make ", this);
  return std::make_shared<Updater>((Updater_Swig *)this);
}
InputCallerPtr InputCaller_Swig::to_input_caller() {
  OPA_DISP("Make ", this);
  return std::make_shared<InputCaller>((InputCaller_Swig *)this);
}

void Updater_Swig::swig_call() { return this->update(); }
bool InputCaller_Swig::swig_call(InputHandler *handler) {
  return update(handler);
}
void Updater_Swig::update() { OPA_CHECK0(false); }
bool InputCaller_Swig::update(InputHandler *handler) {
  OPA_CHECK0(false);
  return false;
}

Game::Game() {
  OPA_CHECK0(GetInstance() == 0);
  s_game_instance = this;
  m_cur = 0;
}

void Game::init() {
  opa::utils::Initable::init();

  m_renderer.reset(new Renderer);
  m_input_handler.reset(new InputHandler);
  m_scene.reset(new Scene);
  m_services.reset(new ServiceManager);
  m_tree.reset(new SceneTreeHelper);
  m_render_streamer.reset(new RenderStreamer);
  m_data_collector.reset(new DataCollector);

  WindowParameters params;
  m_camera.reset(new Camera);
  m_renderer->init(this, params);

  m_input_handler->init(this);
  m_palette.resize_rand(100);

  m_scene->init(this);
  m_services->init(this);

  m_data_collector->init(this);

  m_tree->init();
  register_obj(m_scene, nullptr);

  m_render_streamer->init(m_renderer.get());

  {
    SPTR(PosUpdater) updater(new PosUpdater);
    updater->init(this, UpdateMark(UpdateMark::CameraUpdate + 1));
    m_store.add_singleton<PosUpdater>(updater);
  }

  updater().add(UpdateMark::InputUpdate, UpdaterPtr(new Updater([this]() {
                  m_input_handler->do_step();
                })));
  updater().add(UpdateMark::Draw,
                UpdaterPtr(new Updater([this]() { m_renderer->render(); })));

  updater().add(UpdateMark::CameraUpdate, UpdaterPtr(new Updater([this]() {
                  PosUpdaterPtr ptr = m_store.get<PosUpdater>();
                  ptr->push(m_camera.get());
                })));

  updater().add(UpdateMark::PostDraw, UpdaterPtr(new Updater([this]() {
                  m_render_streamer->update();
                })));

  updater().add(UpdateMark::EndLoop, UpdaterPtr(new Updater([this]() {
                  m_data_collector->update();
                })));

  InputCallerPtr wireframe_cb =
    std::make_shared<InputCaller>([this](InputHandler *handler) {
      if (handler->action().pressed()) renderer()->toggle_wireframe();
      return true;
    });
  m_input_handler->input_cbs()[InputEnum::Key_W].add(wireframe_cb);

  m_input_handler->input_cbs()[InputEnum::Key_I].add(
    std::make_shared<InputCaller>([this](InputHandler *handler) {
      if (!handler->action().pressed()) return false;
      IntersectionResult res;
      Ray ray = camera()->get_click_dir(input_handler()->mouse().cur());

      if (scene()->isec().find_intersection(ray, res, true)) {
        std::vector<SceneObj *> objs =
          tree()->list_children(tree()->get_non_weak(res.obj));
        objs.push_back(res.obj);

        std::vector<opa::math::game::MeshData> meshes;
        OPA_DISP0(res.obj->pos_obj().pos());
        OPA_DISP0(res.obj->pos()->get_self_pos(scene()));
        for (auto &e : objs) {
          if (e->status().drawable) {
            meshes.push_back(
              { e->get_mesh(),
                e->pos()->get_mat_to(glm::mat4(1.), camera()->obj().get()) });
          }
        }
        OPA_DISPL("now for the stuff ", objs.size(), meshes.size());
        for (auto &pt : opa::math::game::to_point_cloud(meshes)) {
          OPA_DISP0(ApplyMat(camera()->get_draw_data().proj(), pt), pt);
        }
      }
      return true;
    }));

  m_services->enable_service(CameraObserver::KEY);
  m_services->enable_service(RotatorService::KEY);

  SceneObjPtr cam_obj(new SceneObj());
  cam_obj->status().drawable = 0;
  cam_obj->init(this);
  cam_obj->pos_obj().pos() = glm::vec3(16, 16, 16);
  cam_obj->pos_obj().rot() = quat_look_at_safe(-vec_y - vec_x - vec_z, vec_z);
  m_camera->obj() = cam_obj;
  register_obj(cam_obj, m_scene.get());

  if (FLAGS_game_default_objs) {
    if (0) {
      auto sphere = std::make_shared<Sphere>();
      sphere->init(this, 1);
      sphere->build(1000);
      sphere->pos_obj().pos() = glm::vec3(0, 0, 0);
      sphere->pos_obj().rot() = quat_look_at(vec_x + vec_y, vec_z);
      register_obj(sphere, m_scene.get());
    }

    if (0) {
      auto square = std::make_shared<Square>();
      float side = 3;
      square->init(this, side, m_palette.get_rand(), {});
      square->pos_obj().pos() = glm::vec3(-side / 2, -side / 2, -3);
      square->pos_obj().rot() = quat_look_at(vec_x, vec_z);
      register_obj(square, m_scene.get());
    }

    if (1) {
      auto cube = std::make_shared<Cube>();
      std::vector<TexturePtr> textures;
      REP (i, 6)
        textures.pb(m_palette.get_rand());

      register_obj(cube, m_scene.get());
      cube->init(this, 2, textures.data());
      cube->pos_obj().pos() = glm::vec3(-5, 0, 0);
    }

    if (2) {
      auto icosahedron = std::make_shared<Icosahedron>();
      std::vector<TexturePtr> textures;
      REP (i, 20)
        textures.pb(m_palette.get_rand());

      icosahedron->init(this, 2, textures.data());
      icosahedron->pos_obj().pos() = glm::vec3(5, 0, 0);
      register_obj(icosahedron, m_scene.get());
    }

    if (0) {
      auto obj = std::make_shared<SceneObj>();
      obj->init(this);
      std::vector<math::game::VertexInfoData> cloud;
      int np = 20;
      REP (i, np)
        cloud.pb(math::game::VertexInfoData(i, vec_rand_uni() * 5));

      auto polyhedra = std::make_shared<DisplayMesh>();
      math::game::MeshBuilder::PolyhedraFromCloud(polyhedra.get(), cloud);

      polyhedra->init();
      auto tex = new MultiColorTexture();
      int nc = 1000;
      tex->init(nc);
      polyhedra->set_texture(TexturePtr(tex));

      int pos = 0;
      for (auto face : polyhedra->list_faces()) {
        for (auto eid : math::game::Mesh::Face_EdgeWalker(*polyhedra, face)) {
          math::game::Edge &e = polyhedra->get_edge(eid);
          math::game::VertexAttribute &attr =
            polyhedra->get_attr(e.vertex_attr);
          attr.uv = tex->get(pos);
        }
        pos = (pos + 1) % nc;
      }

      auto drawer = new mesh::MeshDrawer;
      drawer->init(mesh::MeshDrawer::Params(polyhedra));
      auto isec = new mesh::MeshIsec;
      isec->init(polyhedra);

      obj->drawer().register_drawer(DrawerPtr(drawer));
      obj->isec().register_isec(IsecPtr(isec));
      register_obj(obj, m_scene.get());

      updater().add(
        UpdateMark(UpdateMark::Draw - 1),
        UpdaterPtr(new Updater([polyhedra, tex]() {
          return;

          puts("ISNERTING");
          usleep(1e6);
          polyhedra->orig_map().clear();
          math::game::MeshBuilder::ConvexPolyhedra_Extend(
            polyhedra.get(),
            { math::game::VertexInfoData(
              0, vec_rand_uni() *
                   5) /*, math::game::VertexInfoData(1, vec_rand_uni()*5)*/ });
          if (polyhedra->orig_map().size()) {
            VertexId vid = polyhedra->orig_map().begin()->ST;

            for (auto face :
                 math::game::Mesh::Vertex_FaceWalker(*polyhedra, vid)) {
              int id = opa::math::common::rng64() % tex->size();
              for (auto eid :
                   math::game::Mesh::Face_EdgeWalker(*polyhedra, face)) {
                math::game::Edge &e = polyhedra->get_edge(eid);
                math::game::VertexAttribute &attr =
                  polyhedra->get_attr(e.vertex_attr);
                attr.uv = tex->get(id);
              }
            }
          }
        })));
    }
  }
}

void Game::register_obj(SceneObjPtr obj, SceneObj *par) {
  obj->id() = m_cur++;
  m_tree->add(obj.get(), par);
  m_obj_store[obj->id()] = obj;
}

void Game::check_register(SceneObj *obj) {
  OPA_CHECK(m_obj_store.count(obj->id()) > 0, "forgot to register, obj=%p", obj,
            obj->id());
}

void Game::remove_obj(SceneObj *obj) {
  OPA_DISP("UNREGISTERING ", obj->id());
  auto children = m_tree->list_children(obj);
  m_tree->remove_from_map_only(children);
  m_tree->remove(obj);
  for (auto &child : children) m_obj_store.erase(child->id());
}

SceneObj *Game::get_obj(ObjId id) { return m_obj_store.find(id)->second.get(); }

void Game::run() {
  check_init();

  while (!renderer()->should_quit()) {
    // getc(stdin);
    m_loop_time.us() = glfwGetTime() * 1e6;

    updater().do_loop0<void>();
  }
}

OPA_NAMESPACE_DECL2_END
