#[path = "vec3.rs"] mod vec3;
pub use vec3::*;



#[derive(Clone)]
pub struct Camera {
    pub origin: Point3,
    pub p: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub w: Vec3,
    pub focal_length: f32
}

impl Camera {
    pub fn new(origin: Point3, target: Point3, roll: f32, focal_length: f32) -> Self {
        let w = (target - origin).normalize();
        let p = vec3(roll.sin(), roll.cos(), 0.0);
        let u = w.cross(p).normalize();
        let v = u.cross(w);
        Self {
            origin: origin,
            p: p,
            u, v, w,
            focal_length
        }
    }
    pub fn rotate_x(&mut self, delta: f32) {
        self.w = (self.w + self.u * -delta).normalize();
        self.u = self.w.cross(self.v).normalize();
    }
    pub fn rotate_y(&mut self, delta: f32) {
        self.w = (self.w + self.v * delta).normalize();
        self.v = self.u.cross(self.w);
    }
    pub fn move_x(&mut self, delta: f32) {
        let left = Vec3::new(0., 1., 0.).cross(self.w).normalize() * delta;
        self.origin = self.origin + left;
    }
    pub fn move_y(&mut self, delta: f32) {
        self.origin = self.origin + Vec3::new(0., delta, 0.);
    }
    pub fn move_z(&mut self, delta: f32) {
        self.origin = self.origin + self.w * delta;
    }
    pub fn as_data(&self) -> [[f32; 3]; 3] {
        [
            [self.u.x, self.u.y, self.u.z],
            [self.v.x, self.v.y, self.v.z],
            [self.w.x, self.w.y, self.w.z],
        ]
    }
}