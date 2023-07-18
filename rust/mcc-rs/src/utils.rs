//#![allow(dead_code, unused_imports)] // Muito codigo nao acessado neste momento, supress warnings in this file

use crate::constants::*;
use template_manager::fingerprint_base::{Minutia, PI, TWO_PI};

use libm::{sqrtf, expf, erf, sqrt};
use geo::{LineString, Polygon, ConvexHull, EuclideanDistance};


/// Calculate angle associated with a cell
/// (1) on paper
#[inline(always)]
pub(crate) fn cell_angle(k: i32) -> f32 {
    - PI + (k as f32 - 0.5) * DELTA_D
}

/// Calculate cell-center based on minutiae and pre-calculated sin & cos
/// (2) on paper
#[inline(always)] //always? never?
pub(crate) fn cell_center(m: &Minutia, i: usize, j: usize, sin: f32, cos: f32) -> (i32, i32) {

    let i_term = i as f32 - (NS + 1.0) / 2.0;
    let j_term = j as f32 - (NS + 1.0) / 2.0;

    let x: f32 = 0.5 + m.position.x as f32 + DELTA_S * (cos * i_term + sin * j_term);
    let y: f32 = 0.5 + m.position.y as f32 + DELTA_S * (-sin * i_term + cos * j_term);

    (x as i32, y as i32)
}


/// euclidean distance
#[inline(always)]
pub(crate) fn _ed(pi: &Minutia, pj: &Minutia) -> f32 {
    let dx: f32 = (pi.position.x - pj.position.x) as f32;
    let dy: f32 = (pi.position.y - pj.position.y) as f32;
    
    sqrtf(dx*dx+dy*dy)
}

/// euclidean distance to point
#[inline(always)]
pub(crate) fn ed_to_point(m: &Minutia, x: i32, y: i32) -> f32 {
    let dx: f32 = (m.position.x - x) as f32;
    let dy: f32 = (m.position.y - y) as f32;
    
    sqrtf(dx*dx+dy*dy)
}

/// Parametrized sigmoid
/// (5) on paper
#[inline(always)]
pub(crate) fn sigmoid(v: f32, u: f32, t: f32) -> f32 {
    1.0 / (1.0 + expf(-t * (v - u)))
}

/// Gaussian used in distance contribution calculation
/// (7) on paper
#[inline(always)]
pub(crate) fn gaussian(t: f32) -> f32 {
    (1.0 / (SIGMA_S * SQRT_TWO_PI)) * expf(-(t*t) / (2.0 * SIGMA_S * SIGMA_S))
}

/// Angle difference
/// (9) on paper
#[inline(always)]
pub(crate) fn ang_diff(t1: f32, t2: f32) -> f32 {
    let diff = t1 - t2;
    if (-PI..PI).contains(&diff) { 
        diff
    }
    else if diff < -PI  {
       diff + TWO_PI
    }
    else {
        diff - TWO_PI
    }
}

/// Calculate de convex hull of a set of minutiae
pub(crate) fn convex_hull(minutiae: &[Minutia]) -> Polygon {
    let mut coords = Vec::new();
    //let mut ls = LineString::new();
    for m in minutiae {
        coords.push((m.position.x as f64, m.position.y as f64));
    }

    let poly = Polygon::new(LineString::from(coords), vec![]);
    
    poly.convex_hull()
}

/// checks if a point is inside of a convex hull
/// with an offset omega
#[inline(always)]
pub(crate) fn is_point_inside_chull(chull: &Polygon, px: i32, py: i32, omega: f64) -> bool {
    let p: geo::Point = (px as f64, py as f64).into();

    if chull.euclidean_distance(&p) > omega {
        return false
    }

    true
}

/// Area under curve, calculated using wolfram
/// https://www.wolframalpha.com/input?i2d=true&i=Divide%5B1%2Csigma+*+Sqrt%5B2+*+pi%5D%5D+*+Integrate%5BPower%5Be%2CDivide%5B-Power%5Bt%2C2%5D%2C2*Power%5Bsigma%2C2%5D%5D%5D%2C%7Bt%2Cx+-+Divide%5BD%2C2%5D%2Cx+%2B+Divide%5BD%2C2%5D%7D%5D
/// (11) on paper
#[inline(always)]
pub(crate) fn area_under_gaussian(v: f64) -> f64 {
    0.5 * (erf((0.5 * DELTA_D as f64 + v) / (sqrt(2.0) * SIGMA_D as f64)) - erf((v - 0.5 * DELTA_D as f64) / (sqrt(2.0) * SIGMA_D as f64)))
}

/// Formula to linearize the idxs
/// slightly different from the paper as my arrays start at 0
/// (13) on paper
#[inline(always)]
pub(crate) fn linearize_idxs(i: usize, j: usize, k: usize) -> usize {
    (k) * NS as usize * NS as usize + (j) * NS as usize + i + 1
}