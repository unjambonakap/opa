#pragma once
#include <opa/engine/conf.h>

OPA_NAMESPACE_DECL2(opa, engine)
class DrawData;

struct PointData {
  Pos pos;
  Col col;
};

class PointCloud : public Buffer, public opa::utils::IdObj {
public:
  PointCloud();

  virtual void init(int n);
  virtual void fini();
  void refresh();

  void draw(const DrawData &draw_data);

  OPA_ACCESSOR_R(u32, m_n, n)

  PointData &get_pt(int i) { return m_data[i]; }
  OPA_ACCESSOR(double, m_r_fact, r_fact);

private:
  virtual void init() {}

  double m_r_fact = 0.01;
  u32 m_n;
  std::vector<PointData> m_data;
  Program *m_program;
  GLuint vao;
};
OPA_DECL_SPTR(PointCloud, PointCloudPtr)

OPA_NAMESPACE_DECL2_END
