#include "MipMapBuilder.h"

using namespace opa::utils;
using namespace std;

const int bound = 128;
const int bound2 = 256;
OPA_NAMESPACE_DECL2(opa, engine)

std::string RandomGenerator::get_hex_str(int byte_len) {
  std::string buf(byte_len, 0);
  OPA_CHECK0(read(m_fd, (char *)buf.data(), byte_len) == byte_len);
  return b2h(buf);
}

RandomGenerator::RandomGenerator() { m_fd = open("/dev/urandom", O_RDONLY); }
RandomGenerator::~RandomGenerator() { close(m_fd); }

ResourceManager::ResourceManager() {}
ResourceManager::ResourceManager(const Params &params) { init(params); }

bool ResourceManager::exists(const ResourceId &id) const {
  return get_file(id).exists();
}

void ResourceManager::init(const Params &params) { m_params = params; }

Poco::File ResourceManager::get_file(const ResourceId &id) const {
  std::string filename = stdsprintf("%s%s%s", m_params.prefix.c_str(),
                                    id.c_str(), m_params.suffix.c_str());
  return Poco::File(Poco::Path(m_params.path, filename));
}
Poco::File ResourceManager::open(const ResourceId &id) {
  Poco::File file = get_file(id);
  OPA_CHECK0(file.canRead() && file.isFile());
  return file;
}

std::pair<Poco::File, ResourceManager::ResourceId> ResourceManager::create() {

  RandomGenerator rng;
  while (true) {
    ResourceId id = rng.get_hex_str(m_params.rng_len);
    Poco::File file = get_file(id);
    if (!file.createFile())
      continue;
    return MP(Poco::File(file), id);
  }
}

void MipMapBuilder::init(const Image *image, const std::string &dest_path) {
  m_image = image;
  m_rsc.init(ResourceManager::Params(MM_RSC_PREFIX, "", dest_path));
  m_mmqt.reset(new Mmqt<int>);
}

SPTR(Mmqt<int>) MipMapBuilder::build() {
  m_mmqt->root() = m_mmqt->pool().get_new();
  setup_node(m_mmqt->root(), 0, 0, m_image->w(), m_image->h());
  return m_mmqt;
}

void MipMapBuilder::get_proto(proto::Mmqt &data) const { m_mmqt->store(data); }

ImagePtr MipMapBuilder::setup_node(MmqtNodeId id, int xi, int yi, int w,
                                   int h) {
  std::vector<MmqtNodeId> children;
  bool last = false;
  if (w * h < get_lim())
    last = true;

  if (!last)
    children = m_mmqt->pool().get_new_vec(4);

  {
    Mmqt<int>::Node &node = m_mmqt->pool().get(id);
    node.xo = xi;
    node.yo = yi;
    node.w = w;
    node.h = h;
    REP (i, 4)
      node.tb[i] = InvalidId;
    REP (i, children.size())
      node.tb[i] = children[i];
  }

  ImagePtr res_img;

  if (!last) {
    int w2 = w / 2;
    int h2 = h / 2;
    int wl[] = { w2, w - w2 };
    int hl[] = { h2, h - h2 };
    vector<ImagePtr> imgs;
    ImageResampler::Params params;
    REP (i, 4) {
      int curx = xi + (i % 2) * w2;
      int cury = yi + (i / 2) * h2;
      int curw = wl[i % 2];
      int curh = hl[i / 2];
      imgs.pb(setup_node(children[i], curx, cury, curw, curh));
    }
    w2 = imgs[0]->w();
    h2 = imgs[0]->h();
    REP (i, 4) {
      params.entries.pb(
        ImageResampler::ImageEntry(imgs[i].get(), (i % 2) * w2, (i / 2) * h2));
    }

    int dw = min(bound, 2 * w2);
    int dh = min(bound, 2 * h2);
    ImageResampler resampler;
    res_img = resampler.resample(params, dw, dh);

  } else {
    ImagePtr tmp = m_image->extract(xi, yi, w, h);
    ImageResampler::Params params(
      { ImageResampler::ImageEntry(tmp.get(), 0, 0) });
    ImageResampler resampler;
    res_img = resampler.resample(params, min(bound, w), min(bound, h));
  }

  {
    Mmqt<int>::Node &node = m_mmqt->pool().get(id);
    node.img_w = res_img->w();
    node.img_h = res_img->h();
  }

  auto data = m_rsc.create();
  m_mmqt->pool().get(id).resource = data.ND;
  std::ofstream os(data.ST.path(), std::ofstream::binary);
  res_img->dump(os);

  return res_img;
}

int MipMapBuilder::get_lim() const { return bound2 * bound2; }

OPA_NAMESPACE_DECL2_END
