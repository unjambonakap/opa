#pragma once

#include <opa/math/game/conf.h>
#include <opa/math/game/mesh.h>
#include <opa/math/game/mesh_util.h>
#include <opa/math/game/proto/common.pb.h>
#include <opa/utils/files.h>

OPA_NAMESPACE_DECL3(opa, math, game)
enum ExportType {
  STL = 1,
  OWN = 2,
};

inline void PosToProto(const Pos &pos, proto::Pos *out_pos) {
  out_pos->set_x(pos.x);
  out_pos->set_y(pos.y);
  out_pos->set_z(pos.z);
}

class Exporter {
public:
  ExportType m_type;
  Exporter(ExportType type = ExportType::OWN) : m_type(type) {}

  template <class T> proto::Mesh *to_mesh(const T &mesh, proto::Mesh *out_mesh);
  template <class T>
  proto::Mesh *to_mesh(const T *mesh, proto::Mesh *out_mesh) {
    to_mesh(*mesh, out_mesh);
    return out_mesh;
  }

  template <class T>
  proto::Mesh *to_mesh(const SPTR(T) & mesh, proto::Mesh *out_mesh) {
    to_mesh(*mesh, out_mesh);
    return out_mesh;
  }
  template <class T>
  proto::Mesh *to_mesh(const UPTR(T) & mesh, proto::Mesh *out_mesh) {
    to_mesh(*mesh, out_mesh);
    return out_mesh;
  }

  template <class T>
  proto::MeshList *to_mesh_list(const std::vector<T> &tb,
                                proto::MeshList *out_meshes,
                                const std::string &prefix = "") {
    utils::FilenameSharder sharder;
    if (prefix.size() != 0) {
      sharder.set_pattern(prefix + "_{}").set_number_elem(tb.size()).build();
    }

    for (auto &e : tb) {
      proto::Mesh *out_mesh = out_meshes->add_mesh();
      to_mesh(e, out_mesh);
      if (sharder.valid()) {
        out_mesh->set_name(sharder.get());
      }
    }
    return out_meshes;
  }
};

template <>
inline proto::Mesh *Exporter::to_mesh(const FaceCollection &tc,
                                      proto::Mesh *out_mesh) {
  if (m_type == ExportType::STL) {
    tc.to_stl_buf(out_mesh->mutable_stl_content());
  } else {
    proto::MeshFace *out_face = out_mesh->add_face();
    for (auto &face : tc.faces()) {
      for (auto &pt : face) PosToProto(pt, out_face->add_vertex());
    }
  }
  return out_mesh;
}

template <>
inline proto::Mesh *Exporter::to_mesh(const Mesh &mesh, proto::Mesh *out_mesh) {
  return to_mesh(*mesh.to_faces(), out_mesh);
}
OPA_NAMESPACE_DECL3_END
