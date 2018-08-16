#include "gdal_obj.h"

#include <Poco/Path.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <opa/engine/Intersection.h>
#include <opa/engine/MipMapBuilder.h>
#include <opa/engine/PosUpdater.h>
#include <opa/engine/proto/common.pb.h>
#include <opa/utils/string.h>

using namespace opa::math::game;
using namespace opa::utils;
using namespace std;
const float eps = 1e-5;

OPA_NAMESPACE_DECL2(opa, engine)

BoundingBox2D ScreenBB;
void init_gdal() {
  // GDALAllRegister();
  {
    vector<Pos2> tb;
    REP (i, 4)
      tb.pb(2 * Pos2(i / 2, i % 2) - Pos2(1, 1));
    ScreenBB.from(tb);
  }
}
OPA_REGISTER_INIT(gdal, init_gdal);

BoundingBox2D &BoundingBox2D::from(const std::vector<Pos2> &tb) {
  low.x = std::min_element(ALL(tb), GlmVecXCmp<Pos2>())->x;
  low.y = std::min_element(ALL(tb), GlmVecYCmp<Pos2>())->y;
  high.x = std::max_element(ALL(tb), GlmVecXCmp<Pos2>())->x;
  high.y = std::max_element(ALL(tb), GlmVecYCmp<Pos2>())->y;
  return *this;
}

BoundingBox2D BoundingBox2D::get_intersect(const BoundingBox2D &x) const {
  BoundingBox2D a = *this;
  a.intersect(x);
  return a;
}

BoundingBox2D &BoundingBox2D::intersect(const BoundingBox2D &x) {
  low.x = max(low.x, x.low.x);
  low.y = max(low.y, x.low.y);
  high.x = min(high.x, x.high.x);
  high.y = min(high.y, x.high.y);
  if (low.x > high.x || low.y > high.y) low = Pos2(), high = Pos2();
  return *this;
}
float BoundingBox2D::area() {
  Pos2 v = high - low;
  return v.x * v.y;
}

void AdaptatifMmqt::init(const Params &params) { m_params = params; }

void AdaptatifMmqt::update(const UpdateData &data) {
  // OPA_DISP("got >> ", get(root()).data.count);
  update(data, root());
}

// 3 types of actions:
// 1: need more precision
// 2: use current ndoe
// 3: out of screen -> kill active
void AdaptatifMmqt::update(const UpdateData &data, MmqtNodeId id) {
  if (id == InvalidId) return;

  auto &node = this->get(id);
  vector<Pos> tb;
  // TODO: put it in new function
  Pos2 start(node.xo, node.yo);
  Pos2 dp(1. * node.w / node.img_w, 1. * node.h / node.img_h);
  float precision = dp.x * dp.y;

  REP (i, 4)
    tb.pb(m_params.func->operator()(start +
                                    Pos2(node.w, node.h) * Pos2(i / 2, i % 2)));
  // central -> normal is sum
  Pos normal;
  REP (i, 4)
    normal += tb[i];
  normal = glm::normalize(normal);

  vector<Pos2> proj_tb;
  for (auto &x : tb) {
    auto tmp = *data.proj * Pos4(x, 1);
    proj_tb.pb(tmp.xy() / tmp.w);
  }

  BoundingBox2D bb;
  bb.from(proj_tb);
  float tot_area = bb.area();
  float area_screen = bb.get_intersect(ScreenBB).area();
  float tot_screen = ScreenBB.area(); // 4
  Pos2 img_dim(node.img_w, node.img_h);
  img_dim /= data.viewport;
  float img_ratio = img_dim.x * img_dim.y;
  float rec_ratio = (tot_area / tot_screen) / img_ratio;

  bool is_active = node.data.strip_ids.size() > 0;
  bool want_active = false;
  bool clear_children = false;
  bool is_leaf = node.tb[0] == InvalidId;

  bool drop = OPA_FLOAT_EQ(area_screen, 0, eps);
  drop |= glm::dot(normal, data.front_vec) > 0;

  // OPA_DISP("update >> ", id, rec_ratio, area_screen, node.data.count,
  // proj_tb);
  if (drop) { // case 3
    // OPA_DISP0("Clear children ", area_screen, proj_tb);
    clear_children = true;
  } else if (is_leaf || rec_ratio < 4) { // case 1
    want_active = true;
    clear_children = true;
  } else { // case 2
  }

  if (want_active && is_active) return;
  if (!want_active)
    clear_node(id);
  else {
    activate_node(id);
  }

  if (clear_children) {
    clear(id);
  } else {
    REP (i, 4)
      update(data, node.tb[i]);
    update_node(id);
  }
}
void AdaptatifMmqt::activate_node(MmqtNodeId id) {
  auto &node = this->get(id);

  TexturePtr tex;

  if (m_params.rsc.exists(node.resource)) {
    Poco::File image_file = m_params.rsc.open(node.resource).path();
    std::ifstream ifs(image_file.path(), std::ifstream::binary);
    Image img;
    img.init(node.img_w, node.img_h, Image::ImageFormat::Grey);
    img.load(ifs);
    tex.reset(new Texture);
    tex->init(img);
  } else {
    tex = m_params.default_tex;
  }

  Sampler2D sampler;
  Pos2 start(node.xo, node.yo);
  Pos2 end = start + Pos2(node.w, node.h);
  sampler.init(*m_params.func, m_params.container, tex);
  sampler.sample(start, end, 2, 2);
  node.data.strip_ids = sampler.ids;
  node.data.count++;
}

void AdaptatifMmqt::clear_node(MmqtNodeId id) {
  if (id == InvalidId) return;
  auto &node = this->get(id);
  if (node.data.strip_ids.size() > 0) {
    // OPA_DISP("remove ", node.data.strip_ids.size());
    for (auto x : node.data.strip_ids) m_params.container->strips().remove(x);
    node.data.strip_ids.clear();
    --node.data.count;
  }
}

void AdaptatifMmqt::clear(MmqtNodeId id) {
  if (id == InvalidId) return;
  auto &node = this->get(id);
  if (!node.data.count) return;
  // does not clear self, only children
  REP (i, 4) {
    clear_node(node.tb[i]);
    clear(node.tb[i]);
  }
  update_node(id);
}

void AdaptatifMmqt::update_node(MmqtNodeId id) {
  if (id == InvalidId) return;
  auto &node = this->get(id);

  node.data.count = node.data.strip_ids.size() > 0;
  REP (i, 4) {
    if (node.tb[i] == InvalidId) continue;
    node.data.count += this->get(node.tb[i]).data.count;
  }
}

void Sampler2D::init(const Sampler2DFunc &func, TriangleContainer *container,
                     TexturePtr texture) {
  opa::utils::Initable::init();
  m_func = func;
  m_container = container;
  m_texture = texture;
}

void Sampler2D::sample(Pos2 low, Pos2 high, int xdiv, int ydiv) {
  Pos2 diff = high - low;
  diff.x /= xdiv;
  diff.y /= ydiv;

  vector<vector<pair<Pos, Pos2> > > strips;
  REP (j, ydiv + 1) {
    strips.pb(vector<pair<Pos, Pos2> >());
    REP (i, xdiv + 1) {
      Pos2 cur = low + diff * Pos2(i, j);
      Pos p = m_func(cur);
      strips.back().pb(MP(p, Pos2(i, j) / Pos2(xdiv, ydiv)));
    }

    if (j > 0) {
      auto &strip = m_container->strips().get_new2();
      ids.pb(strip.id);
      strip.init(2 * xdiv, m_texture);
      REP (i, xdiv + 1) {
        strip.data()[i * 2].pos = strips[j - 1][i].ST;
        strip.data()[i * 2].texpos = strips[j - 1][i].ND;
        strip.data()[i * 2 + 1].pos = strips[j][i].ST;
        strip.data()[i * 2 + 1].texpos = strips[j][i].ND;
      }
      strip.refresh();
    }
  }

  /*
     //area computation
  vector<Pos2> tb;
  REP(j, 4) tb.pb(low + diff * Pos2(j / 2, j % 2));

  float area = 0;
  REP(i, xdiv) REP(j, ydiv) {
      vector<Pos> now;
      Pos2 base = Pos(i, j) * diff;
      REP(j, 4) now.pb(m_func(tb[j] + base));

      float sum = 0;
      REP(j, 1) sum += = tr_area(now[1] - now[j * 3], now[2] - now[j * 3]);
      area += sum;
  }
  */

  // sample(low, high, xdiv * ydiv, area);
}

void Sampler2D::sample(Pos2 low, Pos2 high, int ndiv, float area) {
  OPA_ABORT0(false);
}

ImagePtr GdalLoader::get_img(int w, int h) const {
  Image *res = new Image();
  res->init(w, h, Image::ImageFormat::Grey);

  GDALRasterBand *x = m_ds->GetRasterBand(1);
  opa::utils::ScopedHeapAlloc<u16> buf(w * h * 2);
  OPA_CHECK0(!x->RasterIO(GF_Read, 0, 0, x->GetXSize(), x->GetYSize(), buf.buf(), w,
                    h, GDT_UInt16, 0, 0));
  REP (i, w)
    REP (j, h) {
      res->set_col(i, j, Color(min((u16)255, buf[j * w + i]), 0, 0));
    }
  return ImagePtr(res);
}

void GdalLoader::init(const std::string &filename) {
  opa::utils::Initable::init();
  m_ds = (GDALDataset *)GDALOpen(filename.c_str(), GA_ReadOnly);
  OPA_CHECK(m_ds, "Can't open file %s", filename.c_str());

  OGRSpatialReference ref(m_ds->GetProjectionRef());

  m_spheroid.a1 = ref.GetSemiMajor(0);
  m_spheroid.a2 = ref.GetSemiMinor(0);

  double tsf[6];
  OPA_CHECK0(m_ds->GetGeoTransform(tsf) == CE_None);
  m_priv.xo = deg_to_rad(tsf[0]);
  m_priv.yo = deg_to_rad(tsf[3]);
  m_priv.dx = deg_to_rad(tsf[1]);
  m_priv.dy = deg_to_rad(tsf[5]);
  m_priv.nx = m_ds->GetRasterXSize();
  m_priv.ny = m_ds->GetRasterYSize();
}

GdalLoader::~GdalLoader() { fini(); }
void GdalLoader::fini() {
  GDALClose((GDALDatasetH)m_ds);
  opa::utils::Initable::fini();
}

void GdalLoader::setup_proto(proto::GdalData &data) const {
  data.set_x(m_priv.xo);
  data.set_y(m_priv.yo);
  data.set_dx(m_priv.dx);
  data.set_dy(m_priv.dy);
  data.set_nx(m_priv.nx);
  data.set_ny(m_priv.ny);
  data.set_a1(m_spheroid.a1);
  data.set_a2(m_spheroid.a2);
}

std::string GdalLoader::info() {

  std::stringstream ss;
  double adfGeoTransform[6];
  ss << stdsprintf("Driver: %s/%s\n", m_ds->GetDriver()->GetDescription(),
                   m_ds->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));

  ss << stdsprintf("Size is %dx%dx%d\n", m_ds->GetRasterXSize(),
                   m_ds->GetRasterYSize(), m_ds->GetRasterCount());

  if (m_ds->GetProjectionRef() != NULL)
    ss << stdsprintf("Projection is `%s'\n", m_ds->GetProjectionRef());

  if (m_ds->GetGeoTransform(adfGeoTransform) == CE_None) {
    ss << stdsprintf("Origin = (%.6f,%.6f)\n", adfGeoTransform[0],
                     adfGeoTransform[3]);

    ss << stdsprintf("Pixel Size = (%.6f,%.6f)\n", adfGeoTransform[1],
                     adfGeoTransform[5]);
    ss << stdsprintf("Pixel Size 2= (%.6f,%.6f)\n", adfGeoTransform[2],
                     adfGeoTransform[4]);
  }
  return ss.str();
}

void GdalBase::GdalPosSub::notify(Camera *cam) {
  DrawData data = cam->get_draw_data();
  Mat4 mat = data.proj() * m_base->pos()->get_mat_to(glm::mat4(1.), data.base());
  Mat4 imat = glm::inverse(mat);
  Pos4 tmp = imat * Pos4(vec_z, 1);
  Pos front_vec = glm::normalize(tmp.xyz() / tmp.w);
  m_base->m_mmqt.update(
    AdaptatifMmqt::UpdateData(&mat, cam->viewport(), front_vec));
}

void GdalBase::init(Game *game, float scale, const std::string &info_file) {
  SceneObj::init(game);
  auto path = Poco::Path(info_file);
  path.setFileName("");
  m_base_path = path.toString();

  std::ifstream ifs(info_file, std::ifstream::binary);
  OPA_CHECK0(m_data.ParseFromIstream(&ifs));
  m_mmqt.load(m_data.mmqt());

  m_data.set_a1(m_data.a1() * scale);
  m_data.set_a2(m_data.a2() * scale);

  OPA_DISP0(OPA_FIELDS(m_data, a1(), a2(), dx(), dy()));
  Sampler2DFuncPtr func(new Sampler2DFunc([this](Pos2 c) {
    c = Pos2(m_data.x(), m_data.y()) + c * Pos2(m_data.dx(), m_data.dy());
    double v1 = m_data.a1() * cos(c.y), v2 = m_data.a2() * sin(c.y);
    return Pos(v1 * cos(c.x), v1 * sin(c.x), v2);
  }));
  ResourceManager rsc(ResourceManager::Params(MM_RSC_PREFIX, "", m_base_path));
  TexturePtr default_tex = game->palette().get_rand();
  m_mmqt.init(
    AdaptatifMmqt::Params(this, &m_container, func, rsc, default_tex));

  {
    auto x = new TriangleContainerDrawer;
    x->init(TriangleContainerDrawer::Params(&m_container));
    drawer().register_drawer(DrawerPtr(x));
  }

  {
    auto x = new TriangleContainerIsec;
    x->init(TriangleContainerIsecParams(&m_container));
    isec().register_isec(IsecPtr(x));
  }
  PosUpdaterPtr ptr = game->store().get<PosUpdater>();
  ptr->add_sub(PosUpdaterSubPtr(new GdalPosSub(this)));
}

void GdalObj::fini() {
  SceneObj::fini();
  m_entries.clear();
}

void GdalObj::init(Game *game, const std::vector<std::string> &lst) {
  SceneObj::init(game);

  for (int i = 0; i < lst.size(); ++i) {
    auto cur = new GdalBase;
    m_entries.pb(GdalBasePtr(cur));
    cur->init(game, 1. / 2000000, lst[i]);
    cur->fixed() = true;
    game->register_obj(SceneObjPtr(cur), this);
  }

  /*

  GDALRasterBand *poBand;
  int nBlockXSize, nBlockYSize;
  int bGotMin, bGotMax;
  double adfMinMax[2];

  poBand = data->GetRasterBand(1);
  printf(">> %d %d\n", poBand->GetXSize(), poBand->GetYSize());
  poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);
  printf("Block=%dx%d Type=%s, ColorInterp=%s\n", nBlockXSize, nBlockYSize,
         GDALGetDataTypeName(poBand->GetRasterDataType()),
         GDALGetColorInterpretationName(poBand->GetColorInterpretation()));

  adfMinMax[0] = poBand->GetMinimum(&bGotMin);
  adfMinMax[1] = poBand->GetMaximum(&bGotMax);
  if (!(bGotMin && bGotMax))
      GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);

  printf("Min=%.3fd, Max=%.3f\n", adfMinMax[0], adfMinMax[1]);

  if (poBand->GetOverviewCount() > 0)
      printf("Band has %d overviews.\n", poBand->GetOverviewCount());

  if (poBand->GetColorTable() != NULL)
      printf("Band has a color table with %d entries.\n",
             poBand->GetColorTable()->GetColorEntryCount());

  */
}

OPA_NAMESPACE_DECL2_END
