use template_manager::fingerprint_base::PI;

use bitvec::prelude::*;
//use::ordered_float::OrderedFloat;

pub(crate) const CYLINDER_RADIUS: f32 = 70.0;
pub(crate) const NS: f32 = 16.0; // 8.0
pub(crate) const ND: f32 = 6.0;
pub(crate) const DELTA_S: f32 = 2.0 * CYLINDER_RADIUS / NS;
pub(crate) const DELTA_D: f32 = 2.0 * PI / ND;

pub(crate) const SIGMA_S: f32 = 28.0 / 3.0;
pub(crate) const SIGMA_D: f32 = 2.0 * PI / 9.0;
pub(crate) const SQRT_TWO_PI: f32 = 2.5066283;

pub(crate) const OMEGA: f64 = 50.0;

pub(crate) const MIN_VALID_CELLS: f32 = 0.75;  // % minimo de celulas para um cilindro ser valido
pub(crate) const MIN_CONTRIBUTING_MINUTIAE: i32 = 2;        // numero minimo de minucias contribuintes para um cilindro valido
pub(crate) const MAX_ANG_DIFF: f32 = PI / 2.0; // angulo maximo na orientação da minucia/cilindro

pub(crate) const NP_MIN: f32 = 4.0;
pub(crate) const NP_MAX: f32 = 12.0;

pub const CYLINDER_SIZE: usize = (NS * NS * ND) as usize;
// pub type LinearizedCylinder = [f32; CYLINDER_SIZE];

// bitwise params and types
pub(crate) const MICRO_PSI: f32 = 0.01;
pub type LinearizedBitCylinder = BitArray<[u32; CYLINDER_SIZE / 32]>;
//pub type LinearizedBitCylinder = BitArray<[u64; CYLINDER_SIZE / 64]>;

// // relaxation parameters
// const WR: OrderedFloat<f32> = OrderedFloat(0.5);
// const RELAXING_ITERATIONS: usize = 5;