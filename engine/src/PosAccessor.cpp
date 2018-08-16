#include "PosAccessor.h"
#include "Scene.h"
#include "SceneTreeHelper.h"

using namespace opa::math::game;

OPA_NAMESPACE_DECL2(opa, engine)

PosAccessor::PosAccessor() {}
void PosAccessor::init(SceneObj *obj) {
  opa::utils::Initable::init();
  m_obj = obj;
}

const SceneObj *PosAccessor::normalize(const SceneObj *base) const {
  if (base == 0)
    return obj()->game()->scene();
  return base;
}

void PosAccessor::do_mov(const glm::vec3 &dp, const SceneObj *base) {
  obj()->pos_obj().mov_rel(get_dir_from(dp, base));
}

void PosAccessor::do_rot(const Rot &v, const SceneObj *base) {
  obj()->pos_obj().rot_rel(get_rot_from(v, base));
}

Dir PosAccessor::get_front(const SceneObj *base) const {
  return get_dir_to(vec_x, base);
}
Dir PosAccessor::get_left(const SceneObj *base) const {
  return get_dir_to(vec_y, base);
}
Dir PosAccessor::get_up(const SceneObj *base) const {
  return get_dir_to(vec_z, base);
}

Rot PosAccessor::get_rot_from(const Rot &v, const SceneObj *base) {
  return rot_apply_path(base, obj(), v);
}
Rot PosAccessor::get_rot_to(const Rot &v, const SceneObj *base) {
  return rot_apply_path(obj(), base, v);
}

Dir PosAccessor::get_dir_to(const Dir &dir, const SceneObj *base) const {
  return dir_apply_path(obj(), base, dir);
}
Dir PosAccessor::get_dir_from(const Dir &dir, const SceneObj *base) const {
  return dir_apply_path(base, obj(), dir);
}

Pos PosAccessor::get_pos_to(const Pos &pos, const SceneObj *base) const {
  return pos_apply_path(obj(), base, pos);
}
Pos PosAccessor::get_pos_from(const Pos &pos, const SceneObj *base) const {
  return pos_apply_path(base, obj(), pos);
}

glm::mat4 PosAccessor::get_mat_to(const glm::mat4 &mat,
                                  const SceneObj *base) const {
  return mat_apply_path(obj(), base, mat);
}
glm::mat4 PosAccessor::get_mat_from(const glm::mat4 &mat,
                                    const SceneObj *base) const {
  return mat_apply_path(base, obj(), mat);
}

Pos PosAccessor::get_self_pos(const SceneObj *base) const {
  return get_pos_to(Pos(), base);
}

#define OPA_FE_PATH(from, to, op_up, op_down)                                  \
  {                                                                            \
    from = normalize(from);                                                    \
    to = normalize(to);                                                        \
    SceneTreeHelper::TreePath path =                                           \
      obj()->game()->tree()->get_path(from, to);                               \
    for (auto &elem : path) {                                                  \
      if (elem.d == SceneTreeHelper::TreeDir::Up) {                            \
        op_up;                                                                 \
      } else {                                                                 \
        op_down;                                                               \
      }                                                                        \
    }                                                                          \
  }

Rot PosAccessor::rot_apply_path(const SceneObj *from, const SceneObj *to,
                                const Rot &v) const {
  auto cur = v;
  OPA_FE_PATH(from, to, cur = elem.x->obj->pos_obj().rot_to_world() * cur,
              cur = elem.x->obj->pos_obj().rot_to_local() * cur);
  return cur;
}

Pos PosAccessor::pos_apply_path(const SceneObj *from, const SceneObj *to,
                                const Pos &v) const {
  auto cur = v;
  //std::cout << "CUR POS <<" << v << std::endl;
              //std::cout << "CUR POS <<" << cur << std::endl;
  OPA_FE_PATH(from, to, cur = elem.x->obj->pos_obj().get_world_pos(cur);
              , cur = elem.x->obj->pos_obj().get_local_pos(cur));
  return cur;
}
Dir PosAccessor::dir_apply_path(const SceneObj *from, const SceneObj *to,
                                const Dir &v) const {

  auto cur = v;
  OPA_FE_PATH(from, to, cur = elem.x->obj->pos_obj().get_world_dir(cur),
              cur = elem.x->obj->pos_obj().get_local_dir(cur));
  return cur;
}

glm::mat4 PosAccessor::mat_apply_path(const SceneObj *from, const SceneObj *to,
                                      const glm::mat4 &v) const {
  auto cur = v;
  OPA_FE_PATH(from, to, cur = elem.x->obj->pos_obj().mat_to_world() * cur,
              cur = elem.x->obj->pos_obj().mat_to_local() * cur);
  return cur;
}

OPA_NAMESPACE_DECL2_END
