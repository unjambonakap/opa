#pragma once

#include <opa/engine/Game.h>
#include <opa/engine/MipMapBuilder.h>
#include <opa/engine/PosUpdater.h>
#include <opa/engine/Scene.h>
#include <opa/engine/conf.h>
#include <opa/engine/proto/common.pb.h>

class GDALDataset;

OPA_NAMESPACE_DECL2(opa, engine)

typedef std::function<Pos(Pos2)> Sampler2DFunc;
OPA_DECL_SPTR(Sampler2DFunc, Sampler2DFuncPtr);

class BoundingBox2D {
public:
  BoundingBox2D &from(const std::vector<Pos2> &tb);

  BoundingBox2D get_intersect(const BoundingBox2D &x) const;

  BoundingBox2D &intersect(const BoundingBox2D &x);
  float area();

  Pos2 low;
  Pos2 high;
};

// useful?
class CurvatureSmoother {};

struct AdaptatifMmqtData {
  std::vector<utils::IdType> strip_ids;
  int count = 0;
};

class AdaptatifMmqt : public opa::utils::Initable,
                      public Mmqt<AdaptatifMmqtData> {
public:
  struct Params {
    Params(SceneObj *obj, TriangleContainer *container, Sampler2DFuncPtr func,
           ResourceManager rsc, TexturePtr default_tex) {
      this->obj = obj;
      this->container = container;
      this->func = func;
      this->rsc = rsc;
      this->default_tex = default_tex;
    }
    Params() {}

    SceneObj *obj;
    TriangleContainer *container;
    Sampler2DFuncPtr func;
    ResourceManager rsc;
    TexturePtr default_tex;
  };
  struct UpdateData {
    UpdateData(const Mat4 *proj, const Pos2 &viewport, const Pos &front_vec) {
      this->proj = proj;
      this->viewport = viewport;
      this->front_vec = front_vec;
    }

    Pos front_vec;
    const Mat4 *proj;
    IPos2 viewport;
  };

  virtual void init(const Params &params);
  void update(const UpdateData &data);
  void clear(MmqtNodeId id);
  void clear_node(MmqtNodeId id);
  void activate_node(MmqtNodeId id);

  void update_node(MmqtNodeId id);

private:
  void update(const UpdateData &data, MmqtNodeId node_id);
  Params m_params;
};

class Sampler2D : public opa::utils::Initable {
public:
  virtual void init(const Sampler2DFunc &func, TriangleContainer *container,
                    TexturePtr texture);
  void sample(Pos2 low, Pos2 high, int xdiv, int ydiv);
  void sample(Pos2 low, Pos2 high, int ndiv, float area);

  std::vector<utils::IdType> ids;

private:
  virtual void fini() override {}
  virtual void init() override {}
  Sampler2DFunc m_func;
  TriangleContainer *m_container;
  TexturePtr m_texture;
};

class GdalLoader : public opa::utils::Initable {
public:
  void init(const std::string &filename);
  virtual ~GdalLoader();
  std::string info();
  virtual void fini() override;
  ImagePtr get_img(int w, int h) const;

  struct SpheroidParams {
    double a1, a2;
  };
  struct Gdal_SphereData {
    double xo, yo;
    double dx, dy;
    int nx, ny;
  };
  void setup_proto(proto::GdalData &data) const;
  OPA_ACCESSOR_PTR(GDALDataset, m_ds, ds);

private:
  virtual void init() override {}
  Gdal_SphereData m_priv;
  SpheroidParams m_spheroid;
  GDALDataset *m_ds;
};
OPA_DECL_SPTR(GdalLoader, GdalLoaderPtr);

class GdalBase : public SceneObj {
public:
  class GdalPosSub : public PosUpdaterSub {
  public:
    GdalPosSub(GdalBase *base) { m_base = base; }

    virtual void notify(Camera *cam) override;

  private:
    GdalBase *m_base;
  };

  virtual void init(Game *game, float scale, const std::string &info_file);

private:
  virtual void init(Game *game) override {}
  std::string m_base_path;
  float m_scale;
  AdaptatifMmqt m_mmqt;
  TexturePtr m_tex;
  proto::GdalData m_data;
  TriangleContainer m_container;
};
OPA_DECL_SPTR(GdalBase, GdalBasePtr);

class GdalObj : public SceneObj {
public:
  virtual void init(Game *game, const std::vector<std::string> &lst);

private:
  virtual void init(Game *game) override {}
  virtual void fini() override;

  std::vector<GdalBasePtr> m_entries;
};

OPA_NAMESPACE_DECL2_END
