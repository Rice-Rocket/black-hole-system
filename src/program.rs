#[derive(Copy, Clone)]
pub struct Vertex {
    position: [f32; 2],
}
implement_vertex!(Vertex, position);


pub fn load_program(display: &glium::Display) -> (glium::VertexBuffer<Vertex>, glium::index::NoIndices, glium::program::Program) {
    let vertex_buffer = glium::VertexBuffer::new(display, &[
        Vertex { position: [-1.0, -1.0] },
        Vertex { position: [3.0, -1.0] },
        Vertex { position: [-1.0, 3.0] }
    ]).unwrap();
    let indices = glium::index::NoIndices(glium::index::PrimitiveType::TriangleStrip);

    let vertex_shader_src = include_str!("shaders/vertex.glsl");
    let fragment_shader_src = include_str!("shaders/fragment.glsl");

    let program = glium::Program::from_source(display, vertex_shader_src, fragment_shader_src, None).unwrap();
    return (vertex_buffer, indices, program);
}