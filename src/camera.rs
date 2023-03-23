#[path = "primitive.rs"] mod primitive;
pub use primitive::*;



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
        let p = Vec3::new(roll.sin(), roll.cos(), 0.0);
        let u = w.cross(p).normalize();
        let v = u.cross(w).normalize();
        Self {
            origin: origin,
            p: p,
            u, v, w,
            focal_length
        }
    }
    pub fn rotate_x(&mut self, delta: f32) {
        self.w = (self.w + self.u * -delta).normalize();
        self.u = self.w.cross(self.p).normalize();
        self.v = self.u.cross(self.w).normalize();
    }
    pub fn rotate_y(&mut self, delta: f32) {
        self.w = (self.w + self.v * -delta).normalize();
        self.u = self.w.cross(self.p).normalize();
        self.v = self.u.cross(self.w).normalize();
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
            [self.u.x, self.v.x, self.w.x],
            [self.u.y, self.v.y, self.w.y],
            [self.u.z, self.v.z, self.w.z],
        ]
    }
}