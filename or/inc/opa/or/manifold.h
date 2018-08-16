#pragma once

#include <opa/or/grid_search.h>
#include <opa/or/or_common.h>

OPA_NAMESPACE_OR

template <class PointType_> class Manifold {
public:
  typedef PointType_ PointType;
  typedef std::function<double(const PointType &, const PointType &)>
    DistFuncType;

  template <typename Data> struct RequestBase {
    PointType pt; // point to sample neighborhood from
    bool is_grid;
    Data data;
  };

  virtual ~Manifold() {}

  virtual double compute_dist(const PointType &a, const PointType &b) const = 0;
};

class Manifold1D : public Manifold<double> {
public:
  struct Request1D {
    double region;
    int count;
  };
  typedef RequestBase<Request1D> Request;
  virtual std::vector<PointType> sample(const Request &req) const = 0;
  virtual double average_dist(const Request &req) const = 0;
};

class C1Manifold : public Manifold1D {
public:
  C1Manifold(double scale = 1.) : m_scale(scale) {}
  virtual std::vector<PointType> sample(const Request &req) const override {
    if (req.is_grid) {
      return utils::Range<double>::Build_range(0, m_scale, req.data.count).tb();

    } else {
      return utils::Range<double>::Build_range(req.pt - req.data.region,
                                               req.pt + req.data.region,
                                               req.data.count)
        .tb();
    }
  }

  virtual double average_dist(const Request &req) const override {
    double interval = 1.;
    if (!req.is_grid) {
      interval = req.data.region * 2;
    }
    return interval / (req.data.count - 1);
  }

  virtual double compute_dist(const PointType &a,
                              const PointType &b) const override {
    PointType res = std::abs(a - b);
    if (res > m_scale / 2)
      res -= m_scale / 2;
    return res;
  }
  double m_scale;
};

class RManifold : public Manifold1D {
public:
  RManifold(double scale = 1.) : m_scale(scale) {}
  std::vector<PointType> sample(const Request &req) const override {
    if (req.is_grid) {
      return utils::Range<double>::Build_range(0, m_scale, req.data.count).tb();

    } else {
      return utils::Range<double>::Build_range(
               std::max<double>(0, req.pt - req.data.region),
               std::min<double>(m_scale, req.pt + req.data.region),
               req.data.count)
        .tb();
    }
  }

  virtual double average_dist(const Request &req) const override {
    double interval = 1.;
    if (!req.is_grid) {
      interval = req.data.region * 2;
    }
    return interval / (req.data.count - 1);
  }

  virtual double compute_dist(const PointType &a,
                              const PointType &b) const override {
    return std::abs(a - b);
  }
  double m_scale;
};

template <class PointType> class NManifold : public Manifold<PointType> {
public:
  struct DimReq {
    Manifold1D::Request1D req;
    int dim_id;
  };

  struct NReq {
    std::vector<DimReq> dims;
  };

  typedef typename Manifold<PointType>::template RequestBase<NReq> Request;

  NManifold(const std::vector<const Manifold1D *> &submanifolds)
      : m_submanifolds(submanifolds) {}

  std::vector<PointType> sample(const Request &req) const {
    std::vector<std::vector<double> > ranges_1d;
    for (auto &req1 : req.data.dims) {
      int dim = req1.dim_id;
      const auto &cur = get_sub(dim);
      ranges_1d.push_back(
        cur.sample(Manifold1D::Request{ req.pt[dim], req.is_grid, req1.req }));
    }

    std::vector<PointType> res;
    CrossProdGen<double> cpg(ranges_1d, [&](const std::vector<double> &state) {
      PointType cur = req.pt;
      REP (i, req.data.dims.size()) {
        cur[req.data.dims[i].dim_id] = state[i];
      }
      res.push_back(cur);
    });
    return res;
  }

  const Manifold1D &get_sub(int pos) const { return *m_submanifolds[pos]; }
  int ndim() const { return m_submanifolds.size(); }

  virtual std::vector<double> average_dist(const Request &req) const {
    std::vector<double> res(ndim());
    for (auto &req1 : req.data.dims) {
      res[req1.dim_id] = get_sub(req1.dim_id).average_dist(
        Manifold1D::Request{ 0., req.is_grid, req1.req });
    }
    return res;
  }

  std::vector<const Manifold1D *> m_submanifolds;
};

OPA_NAMESPACE_OR_END
