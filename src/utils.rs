#[allow(dead_code)]
#[path = "constants.rs"] mod constants;
pub use constants::*;
pub use std::f32::consts::{PI, TAU};


#[derive(Clone, Copy)]
pub enum Axis {
    X,
    Y,
    Z
}

impl Axis {
    pub fn as_int(&self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2
        }
    }
    pub fn perpendicular(&self, other: &Self) -> Self {
        match self.as_int() + other.as_int() {
            1 => Self::Z,
            2 => Self::Y,
            _ => Self::X
        }
    }
    pub fn others(&self) -> (Self, Self) {
        match self {
            Self::X => (Self::Y, Self::Z),
            Self::Y => (Self::X, Self::Z),
            Self::Z => (Self::X, Self::Y)
        }
    }
}


pub fn to_radians(degrees: f32) -> f32 {
    degrees * PI / 180.0
}

pub fn clamp(x: f32, min: f32, max: f32) -> f32{
    if x < min {
        return min;
    }
    if x > max {
        return max;
    }
    return x;
}