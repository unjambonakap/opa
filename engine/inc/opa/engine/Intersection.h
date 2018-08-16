#pragma once

#include <opa/engine/Game.h>
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)

class TriangleContainer;
class SceneObj;
OPA_DECL_SPTR(TriangleContainer, TriangleContainerPtr)

struct IntersectionResult {
  bool found;
  Pos pos;
  float dist;
  SceneObj *obj;
  std::string debug;
  IntersectionResult() { reset(); }

  void reset() {
    found = false;
    dist = 1e100;
    obj = nullptr;
  }
  bool update(const IntersectionResult &cnd) {
    if (cnd.found && cnd.dist < dist) {
      *this = cnd;
      return true;
    }
    return false;
  }

  void update_success(const Pos &pos, const Ray &ray, SceneObj *obj) {
    found = true;
    this->pos = pos;
    dist = glm::length(pos - ray.pos);
    this->obj = obj;
  }

  Pos get_world_pos();
};
OPA_DECL_SPTR(IntersectionResult, IntersectionResultPtr)

class Isec : public opa::utils::Initable {
public:
  virtual bool find_intersection(const Ray &ray, IntersectionResult &res) = 0;
  void update_intersection(const Ray &ray, IntersectionResult &cur);

  void set_obj(SceneObj *obj) { m_obj = obj; }
  OPA_ACCESSOR_PTR(SceneObj, m_obj, obj);

protected:
  SceneObj *m_obj;
};
OPA_DECL_SPTR(Isec, IsecPtr);

class IntersectionHandler : public opa::utils::Initable {
public:
  virtual void init(SceneObj *obj) {
    opa::utils::Initable::init();
    m_obj = obj;
  }

  virtual bool find_intersection(const Ray &ray, IntersectionResult &res,
                                 bool from_cam = false);
  void update_intersection(const Ray &ray, IntersectionResult &cur);
  void register_isec(IsecPtr isec);

  OPA_ACCESSOR_PTR(SceneObj, m_obj, obj);

protected:
private:
  virtual void init() override {}
  SceneObj *m_obj;
  std::set<IsecPtr> m_isec;
};
OPA_DECL_SPTR(IntersectionHandler, IntersectionHandlerPtr);

struct SphereIsecParams {
  float r;
  SphereIsecParams() {}
  SphereIsecParams(float r) { this->r = r; };
};

class SphereIsec : public Isec {
public:
  virtual void init(const SphereIsecParams &params) {
    Isec::init();
    this->m_params = params;
  }

  virtual bool find_intersection(const Ray &ray,
                                 IntersectionResult &res) override;

private:
  SphereIsecParams m_params;
};
OPA_DECL_SPTR(SphereIsec, SphereIsecPtr);

struct TriangleContainerIsecParams {
  TriangleContainer *tc;
  TriangleContainerIsecParams() {}
  TriangleContainerIsecParams(TriangleContainer *tc) { this->tc = tc; }
};
class TriangleContainerIsec : public Isec {
public:
  virtual bool find_intersection(const Ray &ray,
                                 IntersectionResult &res) override;
  virtual void init(const TriangleContainerIsecParams &params) {
    Isec::init();
    this->m_params = params;
  }

private:
  TriangleContainerIsecParams m_params;
};
OPA_DECL_SPTR(TriangleContainerIsec, TriangleContainerIsecPtr);

OPA_NAMESPACE_DECL2_END
