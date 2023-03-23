#[allow(dead_code)]
use std::ops::{Add, Sub, Mul, Div, Neg, Index, IndexMut};
#[path = "utils.rs"] mod utils;
pub use utils::*;




#[derive(Clone, Copy, PartialEq)]
pub struct Vec2 {
    pub x: f32,
    pub y: f32,
}

impl Add for Vec2 {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Vec2 {
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<f32> for Vec2 {
    type Output = Self;
    fn mul(self, scaler: f32) -> Self::Output {
        Self {
            x: self.x * scaler,
            y: self.y * scaler,
        }
    }
}

impl Mul for Vec2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        Self {
            x: self.x * other.x,
            y: self.y * other.y,
        }
    }
}

impl Div<f32> for Vec2 {
    type Output = Self;
    fn div(self, scaler: f32) -> Self::Output {
        Self {
            x: self.x / scaler,
            y: self.y / scaler,
        }
    }
}

impl Div for Vec2 {
    type Output = Self;
    fn div(self, other: Self) -> Self::Output {
        Self {
            x: self.x / other.x,
            y: self.y / other.y,
        }
    }
}

impl Neg for Vec2 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self * -1.0
    }
}

impl Index<usize> for Vec2 {
    type Output = f32;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Cannot index into Vec2 with {}", index)
        }
    }
}

impl IndexMut<usize> for Vec2 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("Cannot index into Vec2 with {}", index)
        }
    }
}

impl Index<Axis> for Vec2 {
    type Output = f32;
    fn index(&self, index: Axis) -> &Self::Output {
        match index {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            Axis::Z => &self.x
        }
    }
}

impl IndexMut<Axis> for Vec2 {
    fn index_mut(&mut self, index: Axis) -> &mut Self::Output {
        match index {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            Axis::Z => &mut self.x
        }
    }
}

impl Vec2 {
    pub fn new(x: f32, y: f32) -> Self {
        Self {
            x: x,
            y: y,
        }
    }
    pub fn origin() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
        }
    }
    pub fn set_to(&mut self, other: Self) {
        self.x = other.x;
        self.y = other.y;
    }
    pub fn near_zero(&self) -> bool {
        let s = 0.00000001;
        (self.x.abs() < s) && (self.y.abs() < s)
    }
    pub fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y
    }
    pub fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }
    pub fn normalize(&self) -> Self {
        self.clone() / self.length()
    }
    pub fn dot(&self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y
    }
    pub fn to_tuple(&self) -> (f32, f32) {
        (self.x, self.y)
    }
}





#[derive(Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32
}

impl Add<Vec3> for Vec3 {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z
        }
    }
}

impl Add<f32> for Vec3 {
    type Output = Self;
    fn add(self, other: f32) -> Self::Output {
        Self {
            x: self.x + other,
            y: self.y + other,
            z: self.z + other,
        }
    }
}

impl Sub<Vec3> for Vec3 {
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z
        }
    }
}

impl Sub<f32> for Vec3 {
    type Output = Self;
    fn sub(self, other: f32) -> Self::Output {
        Self {
            x: self.x - other,
            y: self.y - other,
            z: self.z - other
        }
    }
}

impl Mul<f32> for Vec3 {
    type Output = Self;
    fn mul(self, scaler: f32) -> Self::Output {
        Self {
            x: self.x * scaler,
            y: self.y * scaler,
            z: self.z * scaler
        }
    }
}

impl Mul for Vec3 {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        Self {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z
        }
    }
}

impl Div<f32> for Vec3 {
    type Output = Self;
    fn div(self, scaler: f32) -> Self::Output {
        Self {
            x: self.x / scaler,
            y: self.y / scaler,
            z: self.z / scaler
        }
    }
}

impl Div for Vec3 {
    type Output = Self;
    fn div(self, other: Self) -> Self::Output {
        Self {
            x: self.x / other.x,
            y: self.y / other.y,
            z: self.z / other.z
        }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self * -1.0
    }
}

impl Index<usize> for Vec3 {
    type Output = f32;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Cannot index into Vec3 with {}", index)
        }
    }
}

impl IndexMut<usize> for Vec3 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Cannot index into Vec3 with {}", index)
        }
    }
}

impl Index<Axis> for Vec3 {
    type Output = f32;
    fn index(&self, index: Axis) -> &Self::Output {
        match index {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            Axis::Z => &self.z
        }
    }
}

impl IndexMut<Axis> for Vec3 {
    fn index_mut(&mut self, index: Axis) -> &mut Self::Output {
        match index {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            Axis::Z => &mut self.z,
        }
    }
}

impl Vec3 {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            x: x,
            y: y,
            z: z
        }
    }
    pub fn origin() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0
        }
    }
    pub fn set_to(&mut self, other: Self) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }
    pub fn near_zero(&self) -> bool {
        let s = 0.00000001;
        (self.x.abs() < s) && (self.y.abs() < s) && (self.z.abs() < s)
    }
    pub fn length_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> f32 {
        self.length_squared().sqrt()
    }
    pub fn normalize(&self) -> Self {
        self.clone() / self.length()
    }
    pub fn abs(&self) -> Self {
        Self {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
    pub fn max(&self, cap: f32) -> Self {
        Self {
            x: self.x.max(cap),
            y: self.y.max(cap),
            z: self.z.max(cap),
        }
    }
    pub fn min(&self, cap: f32) -> Self {
        Self {
            x: self.x.min(cap),
            y: self.y.min(cap),
            z: self.z.min(cap)
        }
    }
    pub fn clamp(&self, min: f32, max: f32) -> Self {
        self.max(min).min(max)
    }
    pub fn lerp(&self, other: Self, t: f32) -> Self {
        Self {
            x: self.x + (other.x - self.x) * t,
            y: self.y + (other.y - self.y) * t,
            z: self.z + (other.z - self.z) * t,
        }
    }
    pub fn dot(&self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(&self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x
        }
    }
    pub fn to_tuple(&self) -> [f32; 3] {
        [self.x, self.y, self.z]
    }
}

pub fn vec3(x: f32, y: f32, z: f32) -> Vec3 {
    Vec3 {
        x: x,
        y: y,
        z: z
    }
}
pub fn vec2(x: f32, y: f32) -> Vec2 {
    Vec2 {
        x: x,
        y: y
    }
}

pub fn point3(x: f32, y: f32, z: f32) -> Vec3 {
    vec3(x, y, z)
}

pub fn rgb(r: f32, g: f32, b: f32) -> Vec3 {
    vec3(r, g, b)
}


pub type Rgb = Vec3;
pub type Point3 = Vec3;




#[derive(Clone)]
pub struct ONB {
    pub u: Vec3,
    pub v: Vec3,
    pub w: Vec3
}

impl ONB {
    pub fn build_from_w(n: Vec3) -> Self {
        let w = n.normalize();
        let a = if w.x.abs() > 0.9 { Vec3::new(0., 1., 0.) } else { Vec3::new(1., 0., 0.) };
        let v = w.cross(a).normalize();
        let u = w.cross(v);
        Self {
            u: u,
            v: v,
            w: w
        }
    }
    pub fn local(&self, a: Vec3) -> Vec3 {
        self.u * a.x + self.v * a.y + self.w * a.z
    }
}

