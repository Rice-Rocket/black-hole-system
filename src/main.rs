#[allow(unused_variables)]
#[macro_use]
extern crate glium;
use glium::glutin;
use glium::Surface;
use std::io::Cursor;


use program::*;
use camera::*;
#[path = "program.rs"] mod program;
#[path = "camera.rs"] mod camera;


fn set_camera() -> Camera {
    return Camera::new(
        point3(3., 0., 2.),
        vec3(0., 0., -1.),
        0.0, 1.0
    );
}

fn input(cam: &mut Camera, held_keys: &[bool; 255], prev_keys: &[bool; 255]) {
    let movement_speed = 0.04;
    let rotate_speed = 0.03;

    if held_keys[glutin::event::VirtualKeyCode::A as usize] {
        cam.move_x(movement_speed); }
    if held_keys[glutin::event::VirtualKeyCode::D as usize] {
        cam.move_x(-movement_speed); }
    if held_keys[glutin::event::VirtualKeyCode::W as usize] {
        cam.move_z(movement_speed); }
    if held_keys[glutin::event::VirtualKeyCode::S as usize] {
        cam.move_z(-movement_speed); }
    if held_keys[glutin::event::VirtualKeyCode::Q as usize] {
        cam.move_y(movement_speed); }
    if held_keys[glutin::event::VirtualKeyCode::E as usize] {
        cam.move_y(-movement_speed); }

    if held_keys[glutin::event::VirtualKeyCode::Left as usize] {
        cam.rotate_x(rotate_speed); }
    if held_keys[glutin::event::VirtualKeyCode::Right as usize] {
        cam.rotate_x(-rotate_speed); }
    if held_keys[glutin::event::VirtualKeyCode::Up as usize] {
        cam.rotate_y(rotate_speed); }
    if held_keys[glutin::event::VirtualKeyCode::Down as usize] {
        cam.rotate_y(-rotate_speed); }

}



fn main() {
    let event_loop = glutin::event_loop::EventLoop::new();
    let wb = glutin::window::WindowBuilder::new();
    let cb = glutin::ContextBuilder::new();
    let display = glium::Display::new(wb, cb, &event_loop).unwrap();
    let (vertex_buffer, indices, program) = load_program(&display);

    let mut cam = set_camera();
    let start_time = std::time::SystemTime::now();
    let mut time = 0f32;
    let mut held_keys = [false; 255];
    let mut prev_keys = [false; 255];
    let mut mouse = [0f32; 4];

    let starmap = image::load(Cursor::new(&include_bytes!("../resources/starmap.png")), 
                                                                image::ImageFormat::Png).unwrap().to_rgb8();
    let starmap_dims = starmap.dimensions();
    let starmap = glium::texture::RawImage2d::from_raw_rgb_reversed(&starmap.into_raw(), starmap_dims);
    let starmap_texture = glium::texture::SrgbTexture2d::new(&display, starmap).unwrap();

    let rgba_noise = image::load(Cursor::new(&include_bytes!("../resources/rgba_noise.png")),
                                                                    image::ImageFormat::Png).unwrap().to_rgba32f();
    let rgba_noise_dims = rgba_noise.dimensions();
    let rgba_noise = glium::texture::RawImage2d::from_raw_rgba_reversed(&rgba_noise.into_raw(), rgba_noise_dims);
    let rgba_noise_tex = glium::texture::SrgbTexture2d::new(&display, rgba_noise).unwrap();

    let organic_texture = image::load(Cursor::new(&include_bytes!("../resources/organic_tex.jpeg")),
                                                                    image::ImageFormat::Jpeg).unwrap().to_rgba32f();
    let organic_texture_dims = organic_texture.dimensions();
    let organic_texture = glium::texture::RawImage2d::from_raw_rgba_reversed(&organic_texture.into_raw(), organic_texture_dims);
    let organic_texture_tex = glium::texture::SrgbTexture2d::new(&display, organic_texture).unwrap();
    event_loop.run(move |ev, _, control_flow| {
        match ev {
            glutin::event::Event::WindowEvent { event, .. } => { match event {
                glutin::event::WindowEvent::CloseRequested => {
                    *control_flow = glutin::event_loop::ControlFlow::Exit;
                    return;
                },
                glutin::event::WindowEvent::KeyboardInput {
                    input: glutin::event::KeyboardInput {
                        virtual_keycode: Some(keycode),
                        state,
                        .. }, .. } => { match state {
                        glutin::event::ElementState::Pressed => {
                            held_keys[keycode as usize] = true;
                        },
                        glutin::event::ElementState::Released => {
                            held_keys[keycode as usize] = false;
                        }
                    }
                },
                glutin::event::WindowEvent::MouseInput { button, state, .. } => { match state {
                        glutin::event::ElementState::Pressed => {
                            match button {
                                glutin::event::MouseButton::Left => { mouse[2] = 1.0; },
                                glutin::event::MouseButton::Right => { mouse[3] = 1.0; },
                                _ => ()
                            }
                        },
                        glutin::event::ElementState::Released => {
                            match button {
                                glutin::event::MouseButton::Left => { mouse[2] = 0.0; },
                                glutin::event::MouseButton::Right => { mouse[3] = 0.0; },
                                _ => ()
                            }
                        }
                    }
                },
                glutin::event::WindowEvent::CursorMoved { position, .. } => {
                    mouse[0] = position.x as f32;
                    mouse[1] = position.y as f32;
                }
                _ => return,
                };
            },
            glutin::event::Event::MainEventsCleared => {
                input(&mut cam, &held_keys, &prev_keys);
            },
            glutin::event::Event::NewEvents(cause) => { match cause {
                    glutin::event::StartCause::ResumeTimeReached { .. } => {
                        display.gl_window().window().request_redraw();
                    },
                    glutin::event::StartCause::Init => (),
                    _ => return,
                }
                prev_keys.copy_from_slice(&held_keys);
            },
            _ => (),
        }

        let next_frame_time = std::time::Instant::now() + 
        std::time::Duration::from_nanos(16_666_667);
        *control_flow = glutin::event_loop::ControlFlow::WaitUntil(next_frame_time);
        let time_now = std::time::SystemTime::now();
        time = time_now.duration_since(start_time).unwrap().as_secs_f32();

        let mut target = display.draw();
        target.clear_color(0.0, 0.0, 0.0, 0.0);
        target.draw(&vertex_buffer, &indices, &program, &uniform! {
            time: time, 
            resolution: [display.get_framebuffer_dimensions().0 as f32, display.get_framebuffer_dimensions().1 as f32],
            mouse: mouse,

            camera: cam.as_data(),
            camera_origin: cam.origin.to_tuple(),
            camera_focal_length: cam.focal_length,

            starmap: &starmap_texture,
            rgba_noise: &rgba_noise_tex,
            organic_tex: &organic_texture_tex,
        }, &Default::default()).unwrap();
        target.finish().unwrap();
    });


}