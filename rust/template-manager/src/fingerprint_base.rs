#![allow(dead_code)] // Muito codigo nao acessado neste momento, supress warnings in this file
use serde::{Serialize, Deserialize};

use std::ops;
//use std::vec;

pub const PI: f32 = std::f32::consts::PI;
pub const TWO_PI: f32 = 2.0 * std::f32::consts::PI;
pub const RAD_TO_DEG: f32 = 180.0 / std::f32::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Point {
    pub x: i32,
    pub y: i32
}

impl Point {
    pub fn compute_angle_degrees(&self, p: &Point) -> i16 {
        compute_angle_degrees((p.x - self.x) as i16, (p.y - self.y) as i16)
    }
}

impl ops::Sub<Point> for Point {
    type Output = Point;

    fn sub(self, _p: Point) -> Point {
        Point {
            x: self.x - _p.x,
            y: self.y - _p.y,
        }
    }
}

//TODO: implement std::cmp::PartialOrd::gt/lt

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Minutia {
    pub position: Point,
    pub direction: i16,
    pub index: u8,
}

impl Minutia {
    pub fn arc(&self) -> f32 {
        2.0 * self.direction as f32 * PI / 180.0 // no original é 2.0 *
    }
}

pub fn compute_angle_diff_180(angle1: i16, angle2: i16) -> i16 {
    let diff = i16::abs_diff(angle1, angle2);
    (if diff < 180 { diff } else { 360 - diff}) as i16
}

pub fn compute_angle_diff_360(angle1: i16, angle2: i16) -> i16 {
    let diff: i16 = angle1 - angle2;
    if diff < 0 { diff + 360 } else { diff }
}

pub fn compute_squared_euclidean_distance(p1: &Point, p2: &Point) -> i32 {
    let dx = (p1.x - p2.x) as f32;
    let dy = (p1.y - p2.y) as f32;

    f32::sqrt(dx * dx + dy * dy) as i32
}

pub fn compute_angle_radians(dx: f32, dy: f32) -> f32 {
    let mut angle = f32::atan2(dy, dx);
    if angle < 0.0 {
        angle += TWO_PI
    }
    angle
}

pub fn compute_angle_degrees(dx: i16, dy: i16) -> i16 {
    (compute_angle_radians(dx as f32, dy as f32) * RAD_TO_DEG + 0.5) as i16
}

pub fn compute_adjacency_matrix(minutiae: &Vec<Minutia>) -> Vec<i32> {
    let len = minutiae.len();
    let mut adjacency = vec![0; len];
    
    for i in 0..len {
        for j in (i+1)..len {
            let dist = compute_squared_euclidean_distance(&minutiae[i].position, &minutiae[j].position);
            adjacency[i * len + j] = dist;
            adjacency[j * len + i] = dist;
        }
    }

    adjacency
}

// // TODO: read template from binary iso
// pub fn read_record_reader() -> i32 {
//     0
// }

// Tests

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angle_in_degrees() {
        assert_eq!((f32::atan2(2.0, 2.0) * RAD_TO_DEG + 0.5) as i16, 45);
    }

    #[test]
    fn subtract_points() {
        let p1 = Point{x:3, y:5};
        let p2 = Point{x:1, y:3};
        let p3 = p1 - p2;

        assert_eq!(p3.x, 2);
        assert_eq!(p3.y, 2);
    }

    #[test]
    fn compare_points() {
        let p1 = Point{x:2, y:2};
        let p2 = Point{x:2, y:2};

        assert_eq!(p1 == p2, true);
    }

    #[test]
    fn angle_between_points() {
        let p1 = Point{x:0, y:0};
        let p2 = Point{x:2, y:2};

        assert_eq!(p1.compute_angle_degrees(&p2), 45);
    }

    #[test]
    fn compare_minutia() {
        let m1: Minutia = Minutia { position: Point{x: 2, y: 2}, direction: 100, index: 1 };
        let m2: Minutia = Minutia { position: Point{x: 2, y: 2}, direction: 100, index: 1 };

        assert_eq!(m1 == m2, true);
    }
    
    #[test]
    fn t_compute_angle_diff_180() {
        let d1 = compute_angle_diff_180(150, 250);
        let d2 = compute_angle_diff_180(50, 250);

        assert_eq!(d1, 100);
        assert_eq!(d2, 160);
    }

    #[test]
    fn t_compute_squared_euclidean_distance() {
        let p1 = Point{x:135, y:225};
        let p2 = Point{x:144, y:281};
    
        let d = compute_squared_euclidean_distance(&p1, &p2);
        assert_eq!(d, 56);
    }

}