#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/proto/common.pb.h>

OPA_NAMESPACE_DECL2(opa, engine)

static const std::string MM_RSC_PREFIX = "oparsc_";

class RandomGenerator {
  public:
    RandomGenerator();
    ~RandomGenerator();
    // len even
    std::string get_hex_str(int byte_len);

  private:
    int m_fd;
};

class ResourceManager : public opa::utils::Initable {
  public:
    typedef std::string ResourceId;
    struct Params {
        std::string prefix;
        std::string suffix;
        std::string path;
        int rng_len;
        Params() {}
        Params(const std::string &prefix, const std::string &suffix, const std::string &path) {
            this->prefix = prefix;
            this->suffix = suffix;
            this->path = path;
            this->rng_len = 8;
        }
    };

    ResourceManager();
    ResourceManager(const Params &params);

    virtual void init(const Params &params);
    bool exists(const ResourceId &id) const;

    Poco::File open(const ResourceId &id);
    std::pair<Poco::File, ResourceId> create();

  private:
    virtual void init() override {}
    Poco::File get_file(const ResourceId &id) const;
    int m_rng_len;
    Params m_params;
};

typedef opa::utils::IdType MmqtNodeId;
template <class T> class Mmqt {
  public:
    struct Node : public opa::utils::IdObj {
        MmqtNodeId tb[4];
        int xo, yo;
        int w, h;
        int img_w;
        int img_h;
        ResourceManager::ResourceId resource;
        T data;

        void init() { REP(i, 4) tb[i] = utils::InvalidId; }

        void load(const proto::MmqtNode &data) {
            id = data.id();
            xo = data.xo();
            yo = data.yo();
            w = data.w();
            h = data.h();
            img_w = data.img_w();
            img_h = data.img_h();
            resource = data.resource();
            REP(i, 4) tb[i] = utils::InvalidId;
            REP(i, data.children().size())
            tb[i] = data.children(i);
        }

        void store(proto::MmqtNode &data) const {
            data.set_xo(xo);
            data.set_yo(yo);
            data.set_id(id);
            data.set_w(w);
            data.set_h(h);
            data.set_img_w(img_w);
            data.set_img_h(img_h);
            data.set_resource(resource);
            REP(i, 4) data.add_children(tb[i]);
        }
    };

    void load(const proto::Mmqt &data) {
        {
            typename opa::utils::ObjectPool<Mmqt<T>::Node>::ObjectPoolLoader loader(m_pool);
            for (auto &x : data.nodes())
                loader.load(x.id());
        }
        m_root = data.root();

        for (auto &x : data.nodes()) {
            Node &node = m_pool.get(x.id());
            node.load(x);
        }
        printf(">> %d\n", data.nodes_size());
    }

    void store(proto::Mmqt &data) const {
        data.set_root(m_root);
        for (auto &x : m_pool.used())
            m_pool.get(x).store(*data.add_nodes());
    }

    OPA_ACCESSOR(opa::utils::ObjectPool<Node>, m_pool, pool);
    OPA_ACCESSOR(MmqtNodeId, m_root, root);
    Node &get(MmqtNodeId id) { return m_pool.get(id); }

  private:
    MmqtNodeId m_root;
    opa::utils::ObjectPool<Node> m_pool;
};

class MipMapBuilder : public opa::utils::Initable {
  public:
    virtual void init(const Image *image, const std::string &dest_path);
    SPTR(Mmqt<int>) build();
    void get_proto(proto::Mmqt &data) const;

  private:
    int get_lim() const;
    ImagePtr setup_node(MmqtNodeId id, int xi, int yi, int w, int h);
    virtual void init() override {}

    ResourceManager m_rsc;
    const Image *m_image;
    std::string m_dest_path;
    SPTR(Mmqt<int>) m_mmqt;
};

OPA_NAMESPACE_DECL2_END
