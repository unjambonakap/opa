#include "mesh/MeshDrawer.h"
#include "resources.h"

OPA_NAMESPACE_DECL3(opa, engine, mesh)

Attribute s1_attrib_pos(Attribute::ATTRIB, "pos");
Attribute s1_attrib_texpos(Attribute::ATTRIB, "texpos");
Attribute s1_attrib_view(Attribute::UNIFORM, "view");
Attribute s1_attrib_texture(Attribute::UNIFORM, "texture1");

Program s1_prog;

Program *get_program_s1() {
    if (!s1_prog.is_init()) {
        s1_prog.init();
        s1_prog.add_shader(GL_FRAGMENT_SHADER, s1_fs);
        s1_prog.add_shader(GL_VERTEX_SHADER, s1_vs);

        s1_prog.register_attribute(&s1_attrib_pos);
        s1_prog.register_attribute(&s1_attrib_texpos);
        s1_prog.register_attribute(&s1_attrib_view);
        s1_prog.register_attribute(&s1_attrib_texture);
        s1_prog.setup();
    }
    return &s1_prog;
}

void MeshDrawer::init(const Params &params) {
    Drawer::init();

    m_params = params;
    m_program = get_program_s1();
}

void MeshDrawer::draw(const DrawData &data) {
    m_params.mesh->prepare_draw();
    check_init();

    m_program->bind();
    ScopedAttribute a1(s1_attrib_pos);
    ScopedAttribute a2(s1_attrib_texpos);
    ScopedAttribute a3(s1_attrib_view);
    ScopedAttribute a4(s1_attrib_texture);
    Buffer &attr_buf = m_params.mesh->buffer();

    GL_CHECK(glUniformMatrix4fv(s1_attrib_view.id(), 1, GL_FALSE, glm::value_ptr(data.proj())));

    glActiveTexture(GL_TEXTURE0);
    m_params.mesh->texture()->bind();
    GL_CHECK(glUniform1i(s1_attrib_texture.id(), 0));


    int id_attr = 0;
    attr_buf.vertex_bind(id_attr, 0, sizeof(opa::math::game::VertexAttribute));

    {
        ScopedBuffer aa(attr_buf);
        glVertexAttribFormat(s1_attrib_texpos.id(), 2, GL_FLOAT, GL_FALSE,
                             offsetof(opa::math::game::VertexAttribute, uv));
        glVertexAttribBinding(s1_attrib_texpos.id(), id_attr);

        glVertexAttribFormat(s1_attrib_pos.id(), 3, GL_FLOAT, GL_FALSE,
                             offsetof(opa::math::game::VertexAttribute, pos));
        glVertexAttribBinding(s1_attrib_pos.id(), id_attr);
    }

    m_params.mesh->vbo().draw();
}

OPA_NAMESPACE_DECL3_END
