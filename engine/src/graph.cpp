#include "graph/graph.h"
#include "resources.h"

const float eps = 1e-7;
using namespace std;
using namespace opa::math::game;

OPA_NAMESPACE_DECL3(opa, engine, graph)

Program point_prog;
Program line_prog;

Attribute point_unif_cam_left(Attribute::UNIFORM, "cam_left");
Attribute point_unif_cam_up(Attribute::UNIFORM, "cam_up");
Attribute point_unif_view(Attribute::UNIFORM, "view");
Attribute point_attr_sq_pos(Attribute::ATTRIB, "sq_pos");
Attribute point_attr_pos(Attribute::ATTRIB, "pos");
Attribute point_attr_color(Attribute::ATTRIB, "color");

Program *get_program_point() {
  if (!point_prog.is_init()) {
    point_prog.init();
    point_prog.add_shader(GL_FRAGMENT_SHADER, point_fs);
    point_prog.add_shader(GL_VERTEX_SHADER, point_vs);

    point_prog.register_attribute(&point_attr_pos);
    point_prog.register_attribute(&point_attr_color);
    point_prog.register_attribute(&point_attr_sq_pos);
    point_prog.register_attribute(&point_unif_cam_left);
    point_prog.register_attribute(&point_unif_cam_up);
    point_prog.register_attribute(&point_unif_view);
    point_prog.setup();
  }
  return &point_prog;
}

Attribute line_attr_sq_pos(Attribute::ATTRIB, "sq_pos");
Attribute line_attr_pos1(Attribute::ATTRIB, "pos1");
Attribute line_attr_pos2(Attribute::ATTRIB, "pos2");
Attribute line_attr_color(Attribute::ATTRIB, "color");
Attribute line_unif_cam_left(Attribute::UNIFORM, "cam_left");
Attribute line_unif_cam_up(Attribute::UNIFORM, "cam_up");
Attribute line_unif_view(Attribute::UNIFORM, "view");

Program *get_program_line() {
  if (!line_prog.is_init()) {
    line_prog.init();
    line_prog.add_shader(GL_FRAGMENT_SHADER, line_fs);
    line_prog.add_shader(GL_VERTEX_SHADER, line_vs);

    line_prog.register_attribute(&line_attr_pos1);
    line_prog.register_attribute(&line_attr_pos2);
    line_prog.register_attribute(&line_attr_color);
    line_prog.register_attribute(&line_attr_sq_pos);
    line_prog.register_attribute(&line_unif_cam_left);
    line_prog.register_attribute(&line_unif_cam_up);
    line_prog.register_attribute(&line_unif_view);
    line_prog.setup();
  }
  return &line_prog;
}

OPA_NAMESPACE_DECL3_END
