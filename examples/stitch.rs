
use tri_mesh::prelude::*;
use tri_mesh::prelude::Vec3 as Vec3;
use tri_mesh::prelude::vec3 as vec3;
use tri_mesh::prelude::vec4 as vec4;

/// Loads the meshes and scale/translate it.
fn on_startup(scene_center: &Vec3, scene_radius: f64) -> (Mesh, Mesh)
{
    let mut mesh = MeshBuilder::new().with_obj(include_str!("assets/bunny.obj").to_string()).build().unwrap();
    transform(&mut mesh, scene_center, scene_radius);
    let mut other_mesh = MeshBuilder::new().with_obj(include_str!("assets/blob.obj").to_string()).build().unwrap();
    transform(&mut other_mesh, scene_center, 0.5 * scene_radius);
    (mesh, other_mesh)
}

/// Translates the mesh to the scene center and scales it such that the size of the biggest side of the bounding box is half a scene radius
fn transform(mesh: &mut Mesh, scene_center: &Vec3, scene_radius: f64)
{
    let (min, max) = mesh.extreme_coordinates();
    mesh.translate(-0.5 * (max + min)); // Translate such that the mesh center is in origo.
    let size = max - min;
    mesh.scale(0.5 * scene_radius / size.x.max(size.y).max(size.z)); // Scale the mesh such that the size of the biggest side of the bounding box is half a scene radius
    mesh.translate(*scene_center); // Translate the mesh to the scene center
}

/// When the user clicks, we see if the model is hit. If it is, we
/// * translate the other mesh to the picked position,
/// * splits the meshes at their intersection and
/// * merges the sub meshes into new combinations.
fn on_click(mesh: &mut Mesh, other_mesh: &mut Mesh, ray_start_point: &Vec3, ray_direction: &Vec3) -> Option<Vec<Mesh>>
{
    if let Some(Intersection::Point {point, ..}) = mesh.ray_intersection(ray_start_point, ray_direction) {
        // Translate
        other_mesh.translate(point - other_mesh.axis_aligned_bounding_box_center());

        // Split at intersection
        let (meshes1, meshes2) = mesh.split_at_intersection(other_mesh);

        // Merge sub meshes
        let mut result_meshes = Vec::new();
        for mesh1 in meshes1.iter() {
            for mesh2 in meshes2.iter() {
                let mut result = mesh1.clone();
                if result.merge_with(mesh2).is_ok()
                {
                    result_meshes.push(result);
                }
            }
        }
        if result_meshes.len() == 0 { None } else {Some(result_meshes)}
    }
    else {None}
}

///
/// Above: Everything related to tri-mesh
/// Below: Visualisation of the mesh, event handling and so on
///
use dust::*;
use dust::objects::*;
use dust::window::{event::*, Window};

fn main()
{
    let scene_radius = 10.0;
    let scene_center = dust::vec3(0.0, 5.0, 0.0);
    let (mut mesh, mut other_mesh) = on_startup(&vec3(scene_center.x as f64, scene_center.y as f64, scene_center.z as f64), scene_radius as f64);

    let mut window = Window::new_default("Stitch tool").unwrap();
    let (framebuffer_width, framebuffer_height) = window.framebuffer_size();
    let window_size = window.size();
    let gl = window.gl();

    // Renderer
    let renderer = DeferredPipeline::new(&gl, framebuffer_width, framebuffer_height, true, vec4(0.8, 0.8, 0.8, 1.0)).unwrap();

    // Camera
    let mut camera = camera::PerspectiveCamera::new(scene_center + scene_radius * vec3(1.0, 1.0, 1.0).normalize(), scene_center,
                                                    vec3(0.0, 1.0, 0.0),degrees(45.0), framebuffer_width as f32 / framebuffer_height as f32, 0.1, 1000.0);

    // Objects
    let positions: Vec<f32> = mesh.positions_buffer().iter().map(|v| *v as f32).collect();
    let normals: Vec<f32> = mesh.normals_buffer().iter().map(|v| *v as f32).collect();

    let mut wireframe_model = Wireframe::new(&gl, &mesh.indices_buffer(), &positions, 0.02);
    wireframe_model.set_parameters(0.8, 0.2, 5.0);
    wireframe_model.set_color(&vec3(0.9, 0.2, 0.2));

    let mut model = ShadedMesh::new(&gl, &mesh.indices_buffer(), &att!["position" => (positions, 3), "normal" => (normals, 3)]).unwrap();
    model.color = vec3(0.8, 0.8, 0.8);

    let other_positions: Vec<f32> = other_mesh.positions_buffer().iter().map(|v| *v as f32).collect();
    let other_normals: Vec<f32> = other_mesh.normals_buffer().iter().map(|v| *v as f32).collect();
    let mut other_model = ShadedMesh::new(&gl, &other_mesh.indices_buffer(), &att!["position" => (other_positions, 3), "normal" => (other_normals, 3)]).unwrap();
    other_model.color = vec3(0.8, 0.8, 0.8);

    let plane_positions: Vec<f32> = vec![
        -1.0, 0.0, -1.0,
        1.0, 0.0, -1.0,
        1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0
    ];
    let plane_normals: Vec<f32> = vec![
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0
    ];
    let plane_indices: Vec<u32> = vec![
        0, 2, 1,
        0, 3, 2,
    ];
    let mut plane = ShadedMesh::new(&gl, &plane_indices, &att!["position" => (plane_positions, 3), "normal" => (plane_normals, 3)]).unwrap();
    plane.diffuse_intensity = 0.2;
    plane.specular_intensity = 0.4;
    plane.specular_power = 20.0;

    let mut ambient_light = light::AmbientLight::new();
    ambient_light.base.intensity = 0.4;

    let mut dir = dust::vec3(-1.0, -1.0, -1.0).normalize();
    let mut light1 = light::SpotLight::new(scene_center - 2.0 * scene_radius * dir, dir);
    light1.enable_shadows(&gl, scene_radius * 4.0).unwrap();
    light1.base.intensity = 0.75;

    dir = dust::vec3(-1.0, -1.0, 1.0).normalize();
    let mut light2 = light::SpotLight::new(scene_center - 2.0 * scene_radius * dir, dir);
    light2.enable_shadows(&gl, scene_radius * 4.0).unwrap();
    light2.base.intensity = 0.75;

    dir = dust::vec3(1.0, -1.0, 1.0).normalize();
    let mut light3 = light::SpotLight::new(scene_center - 2.0 * scene_radius * dir, dir);
    light3.enable_shadows(&gl, scene_radius * 4.0).unwrap();
    light3.base.intensity = 0.75;

    dir = dust::vec3(1.0, -1.0, -1.0).normalize();
    let mut light4 = light::SpotLight::new(scene_center - 2.0 * scene_radius * dir, dir);
    light4.enable_shadows(&gl, scene_radius * 4.0).unwrap();
    light4.base.intensity = 0.75;

    let mut camera_handler = camerahandler::CameraHandler::new(camerahandler::CameraState::SPHERICAL);

    let mut chosen = 0;
    let mut results: Option<Vec<(Mesh, ShadedMesh, Wireframe)>> = None;
    // main loop
    window.render_loop(move |events, _elapsed_time|
    {
        for event in events {
            match event {
                Event::Key {state, kind} => {
                    if kind == "Tab" && *state == State::Pressed
                    {
                        camera_handler.next_state();
                    }
                },
                Event::MouseClick {state, button, position} => {
                    if *button == MouseButton::Left
                    {
                        if *state == State::Pressed
                        {
                            let (x, y) = (position.0 / window_size.0 as f64, position.1 / window_size.1 as f64);
                            let p = camera.position();
                            let ray_start_point = vec3(p.x as f64, p.y as f64, p.z as f64);
                            let dir = camera.view_direction_at((x, y));
                            let ray_direction = vec3(dir.x as f64, dir.y as f64, dir.z as f64);

                            if let Some(ref result) = results {
                                if let Some(Intersection::Point {..}) = result[chosen].0.ray_intersection(&ray_start_point, &ray_direction) {
                                    chosen = (chosen + 1) % result.len();
                                }
                                else {
                                    camera_handler.start_rotation();
                                }
                            }
                            else {
                                if let Some(mut meshes) = on_click(&mut mesh, &mut other_mesh, &ray_start_point, &ray_direction) {
                                    let mut result = Vec::new();
                                    for mesh in meshes.drain(..) {
                                        let positions: Vec<f32> = mesh.positions_buffer().iter().map(|v| *v as f32).collect();
                                        let normals: Vec<f32> = mesh.normals_buffer().iter().map(|v| *v as f32).collect();

                                        let mut wireframe_model = Wireframe::new(&gl, &mesh.indices_buffer(), &positions, 0.02);
                                        wireframe_model.set_parameters(0.8, 0.2, 5.0);
                                        wireframe_model.set_color(&vec3(0.9, 0.2, 0.2));

                                        let mut model = ShadedMesh::new(&gl, &mesh.indices_buffer(), &att!["position" => (positions, 3), "normal" => (normals, 3)]).unwrap();
                                        model.color = vec3(0.8, 0.8, 0.8);

                                        result.push((mesh, model, wireframe_model));
                                    }
                                    results = Some(result);

                                }
                                else {
                                    camera_handler.start_rotation();
                                }
                            }
                        }
                        else {
                            camera_handler.end_rotation()
                        }
                    }
                    else if *button == MouseButton::Right && *state == State::Pressed
                    {
                        results = None;
                        chosen = 0;
                    }
                },
                Event::MouseWheel {delta} => {
                    camera_handler.zoom(&mut camera, *delta as f32);
                },
                Event::MouseMotion {delta} => {
                    camera_handler.rotate(&mut camera, delta.0 as f32, delta.1 as f32);
                }
            }
        }

        // Draw
        let render_scene = |camera: &Camera| {
            let model_matrix = dust::Mat4::identity();
            if let Some(ref result) = results
            {
                result[chosen].1.render(&model_matrix, camera);
                result[chosen].2.render(camera);
            }
            else {
                model.render(&model_matrix, camera);
                wireframe_model.render(camera);

                let camera_pos = camera.position();
                let trans = dust::Mat4::from_translation(0.5 * scene_radius * vec3(-camera_pos.z, 0.0, camera_pos.x).normalize() - vec3(scene_center.x, 0.0, scene_center.z));
                other_model.render(&trans, camera);
            }
        };

        // Shadow pass
        light1.shadow_cast_begin().unwrap();
        render_scene(light1.shadow_camera().unwrap());

        light2.shadow_cast_begin().unwrap();
        render_scene(light2.shadow_camera().unwrap());

        light3.shadow_cast_begin().unwrap();
        render_scene(light3.shadow_camera().unwrap());

        light4.shadow_cast_begin().unwrap();
        render_scene(light4.shadow_camera().unwrap());

        // Geometry pass
        renderer.geometry_pass_begin().unwrap();
        render_scene(&camera);
        plane.render(&dust::Mat4::from_scale(100.0), &camera);

        // Light pass
        renderer.light_pass_begin(&camera).unwrap();
        renderer.shine_ambient_light(&ambient_light).unwrap();
        renderer.shine_spot_light(&light1).unwrap();
        renderer.shine_spot_light(&light2).unwrap();
        renderer.shine_spot_light(&light3).unwrap();
        renderer.shine_spot_light(&light4).unwrap();

        renderer.copy_to_screen().unwrap();
    }).unwrap();
}