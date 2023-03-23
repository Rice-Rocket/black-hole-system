#[allow(dead_code, unused)]
#[path = "vec3.rs"] mod vec3;

pub use vec3::*;




#[derive(PartialEq, Clone)]
pub enum Primitive {
    Sphere(f32, Rgb), // Radius
    Plane(Vec3, f32, Rgb), // Plane Normal, Distance Along Normal
    Cuboid(Vec3, f32, Rgb), // Dimensions, Rounding Value
    BoxFrame(Vec3, f32, Rgb), // Dimensions, Edge Thickness
    Torus(f32, f32, Rgb), // Inner Radius, Outer Radius
    Horseshoe(f32, f32, f32, f32, Rgb), // Inner Radius, Outer Radius, Cap A, Cap B
    Link(f32, f32, f32, Rgb), // Length, Radius 1, Radius 2
    Cone(f32, f32, Rgb), // Angle, Height
    HexagonalPrism(f32, f32, Rgb), // Radius, Height
    TriangularPrism(f32, f32, Rgb), // Radius, Height
    Capsule(f32, f32, Rgb), // Height, Radius
    CappedCylinder(f32, f32, Rgb), // Height, Radius
    RoundedCylinder(f32, f32, f32, Rgb), // Radius Max, Radius Min, Height
    CappedCone(f32, f32, f32, Rgb), // Height, Radius 1, Radius 2
    SolidAngle(f32, f32, Rgb), // Angle, Rounding
    CutSphere(f32, f32, Rgb), // Radius, Height
    CutHollowSphere(f32, f32, f32, Rgb), // Radius, Height, Thickness
    DeathStar(f32, f32, f32, Rgb), // Radius 1, Radius 2, d
    RoundCone(f32, f32, f32, Rgb), // Radius 1, Radius 2, Height
    Ellipsoid(Vec3, Rgb), // Radii
    Rhombus(f32, f32, f32, f32, Rgb), // la, lb, Height, ra
    Octahedron(f32, Rgb), // Side length
    Pyramid(f32, Rgb), // Height
    Triangle(Point3, Point3, Point3, Rgb), // Point 1, Point 2, Point 3
}

impl Primitive {
    pub fn id(&self) -> f32 {
        match &self {
            Self::Sphere(..) => 1.,
            Self::Plane(..) => 2.,
            Self::Cuboid(..) => 3.,
            Self::BoxFrame(..) => 4.,
            Self::Torus(..) => 5.,
            Self::Horseshoe(..) => 6.,
            Self::Link(..) => 7.,
            Self::Cone(..) => 8.,
            Self::HexagonalPrism(..) => 9.,
            Self::TriangularPrism(..) => 10.,
            Self::Capsule(..) => 11.,
            Self::CappedCylinder(..) => 12.,
            Self::RoundedCylinder(..) => 13.,
            Self::CappedCone(..) => 14.,
            Self::SolidAngle(..) => 15.,
            Self::CutSphere(..) => 16.,
            Self::CutHollowSphere(..) => 17.,
            Self::DeathStar(..) => 18.,
            Self::RoundCone(..) => 19.,
            Self::Ellipsoid(..) => 20.,
            Self::Rhombus(..) => 21.,
            Self::Octahedron(..) => 22.,
            Self::Pyramid(..) => 23.,
            Self::Triangle(..) => 24.,
        }
    }
    pub fn as_str(&self) -> String {
        match &self {
            Self::Sphere(..) => "Sphere".to_string(),
            Self::Plane(..) => "Plane".to_string(),
            Self::Cuboid(..) => "Cuboid".to_string(),
            Self::BoxFrame(..) => "Box Frame".to_string(),
            Self::Torus(..) => "Torus".to_string(),
            Self::Horseshoe(..) => "Horseshoe".to_string(),
            Self::Link(..) => "Link".to_string(),
            Self::Cone(..) => "Cone".to_string(),
            Self::HexagonalPrism(..) => "Hexagonal Prism".to_string(),
            Self::TriangularPrism(..) => "Triangular Prism".to_string(),
            Self::Capsule(..) => "Capsule".to_string(),
            Self::CappedCylinder(..) => "Capped Cylinder".to_string(),
            Self::RoundedCylinder(..) => "Rounded Cylinder".to_string(),
            Self::CappedCone(..) => "Capped Cone".to_string(),
            Self::SolidAngle(..) => "Solid Angle".to_string(),
            Self::CutSphere(..) => "Cut Sphere".to_string(),
            Self::CutHollowSphere(..) => "Cut Hollow Sphere".to_string(),
            Self::DeathStar(..) => "Death Star".to_string(),
            Self::RoundCone(..) => "Round Cone".to_string(),
            Self::Ellipsoid(..) => "Ellipsoid".to_string(),
            Self::Rhombus(..) => "Rhombus".to_string(),
            Self::Octahedron(..) => "Octahedron".to_string(),
            Self::Pyramid(..) => "Pyramid".to_string(),
            Self::Triangle(..) => "Triangle".to_string(),
        }
    }
    pub fn as_data(&self, op_group: f32) -> [[f32; 4]; 4] {
        match &self {
            Self::Sphere(rad, color) => [
                [self.id(), *rad, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Plane(normal, k, color) => [
                [self.id(), *k, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [normal.x, normal.y, normal.z, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Cuboid(dims, rounding, color) => [
                [self.id(), *rounding, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [dims.x, dims.y, dims.z, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::BoxFrame(dims, edge, color) => [
                [self.id(), *edge, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [dims.x, dims.y, dims.z, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Torus(in_rad, out_rad, color) => [
                [self.id(), *in_rad, *out_rad, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Horseshoe(in_rad, out_rad, capa, capb, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*in_rad, *out_rad, *capa, *capb],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Link(length, rad1, rad2, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*length, *rad1, *rad2, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Cone(angle, height, color) => [
                [self.id(), *angle, *height, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::HexagonalPrism(radius, height, color) => [
                [self.id(), *radius, *height, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::TriangularPrism(radius, height, color) => [
                [self.id(), *radius, *height, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Capsule(height, radius, color) => [
                [self.id(), *height, *radius, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::CappedCylinder(height, radius, color) => [
                [self.id(), *height, *radius, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::RoundedCylinder(rad1, rad2, height, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*rad1, *rad2, *height, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::CappedCone(height, rad1, rad2, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*height, *rad1, *rad2, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::SolidAngle(angle, rounding, color) => [
                [self.id(), *angle, *rounding, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::CutSphere(radius, height, color) => [
                [self.id(), *radius, *height, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::CutHollowSphere(radius, height, thickness, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*radius, *height, *thickness, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::DeathStar(radius1, radius2, d, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*radius1, *radius2, *d, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::RoundCone(rad1, rad2, height, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*rad1, *rad2, *height, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Ellipsoid(radii, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [radii.x, radii.y, radii.z, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Rhombus(la, lb, height, ra, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [*la, *lb, *height, *ra],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Octahedron(side, color) => [
                [self.id(), *side, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Pyramid(height, color) => [
                [self.id(), *height, 0.0, op_group],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [color.x, color.y, color.z, 0.0]],
            Self::Triangle(p1, p2, p3, color) => [
                [self.id(), 0.0, 0.0, op_group],
                [p1.x, p1.y, p1.z, p3.x],
                [p2.x, p2.y, p2.z, p3.y],
                [color.x, color.y, color.z, p3.z]],
        }
    }
    // pub fn center(&self) -> Point3 {
    //     match &self {

    //     }
    // }
}


#[derive(Clone, Copy)]
pub enum BooleanOpType {
    Union,
    Intersect,
    Subtract,
    SmoothUnion(f32),
    SmoothIntersect(f32),
}



pub struct BooleanOp {
    operator: BooleanOpType,
    pub obj_uids: Vec<usize>,
    pub uid: usize
}

// * Csgs have indices to two primitives. csgs just dictate how those two act together. 
// * in the end, there will be a list of csgs that will be passed to the shader
impl BooleanOp {
    pub fn new(operator: BooleanOpType, objs: Vec<usize>) -> Self {
        Self {
            operator: operator,
            obj_uids: objs,
            uid: 0
        }
    }
    pub fn get(&self) -> [f32; 2] {
        match self.operator {
            BooleanOpType::Union => [1., 0.0],
            BooleanOpType::Intersect => [2., 0.0],
            BooleanOpType::Subtract => [3., 0.0],
            BooleanOpType::SmoothUnion(k) => [4., k],
            BooleanOpType::SmoothIntersect(k) => [5., k],
        }
    }
}



#[derive(Clone, Copy)]
pub struct Transform {
    pub translate: Vec3,
    pub rotate: Vec3,
    pub scale: Vec3,
    // pub shear: Vec3,
}

impl Transform {
    pub fn new(translate: Vec3, rotate: Vec3, scale: Vec3) -> Self {
        Self {
            translate, 
            rotate, 
            scale
        }
    }
    pub fn translation(translate: Vec3) -> Self {
        Self {
            translate,
            rotate: vec3(0., 0., 0.),
            scale: vec3(1., 1., 1.)
        }
    }
    pub fn rotation(rotate: Vec3) -> Self {
        Self {
            translate: vec3(0., 0., 0.),
            rotate,
            scale: vec3(1., 1., 1.)
        }
    }
    pub fn scaler(scale: Vec3) -> Self {
        Self {
            translate: vec3(0., 0., 0.),
            rotate: vec3(0., 0., 0.),
            scale: scale
        }
    }
    pub fn transrot(translate: Vec3, rotate: Vec3) -> Self {
        Self {
            translate,
            rotate,
            scale: vec3(1., 1., 1.)
        }
    }
    pub fn transscale(translate: Vec3, scale: Vec3) -> Self {
        Self {
            translate,
            rotate: vec3(0., 0., 0.),
            scale
        }
    }
    pub fn rotscale(rotate: Vec3, scale: Vec3) -> Self {
        Self {
            translate: vec3(0., 0., 0.),
            rotate,
            scale
        }
    }
    pub fn none() -> Self {
        Self {
            translate: vec3(0., 0., 0.),
            rotate: vec3(0., 0., 0.),
            scale: vec3(1., 1., 1.)
        }
    }
    pub fn get_data(&self) -> [[f32; 4]; 4] {
        [
            [self.translate.x, self.translate.y, self.translate.z, 0.0],
            [self.rotate.x, self.rotate.y, self.rotate.z, 0.0],
            [self.scale.x, self.scale.y, self.scale.z, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    }
}