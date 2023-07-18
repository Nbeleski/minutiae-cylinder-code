#![allow(dead_code, unused_imports)] // Muito codigo nao acessado neste momento, supress warnings in this file

use crate::constants::*;
use crate::utils::*;

use ndarray::Order;
use template_manager::fingerprint_base::{Point, Minutia, PI, TWO_PI};

extern crate ndarray;
use ndarray::{Array2, Array3};
use libm::{sqrtf, cosf, sinf, expf, erf, sqrt, atan2f};
use geo::{LineString, Polygon, ConvexHull, EuclideanDistance};
use::ordered_float::OrderedFloat;
use pathfinding::prelude::{kuhn_munkres, Matrix};
use bitvec::prelude::*;

/// Bit-implementation Cylinder set for a set of minutiae
pub struct BitCylinderSet {
    pub cylinders:  Vec<LinearizedBitCylinder>,
    validities: Vec<LinearizedBitCylinder>,
    minutiae_reference: Vec<u8>,
    // angle?
}

/// a Cylinder set representation of a minutiae set
/// where each cylinder is linearized...
/// the object keeps a reference to the original minutiae index
/// as it is needed to check orientation during matching
/// (or coordinates during relaxation)
impl BitCylinderSet {
    pub fn new(minutiae: &[Minutia]) -> Self {

        let mut cylinders = Vec::with_capacity(minutiae.len());
        let mut validities = Vec::with_capacity(minutiae.len());
        let mut minutiae_reference = Vec::with_capacity(minutiae.len());

        // convex hull used to check cell validity
        let chull = convex_hull(minutiae);

        // calculate all phi vals
        let phi_vals: Vec<f32> = (1..ND as i32 + 1).map(cell_angle).collect();

        // for each minutiae in the set
        for (idx_probe, m) in minutiae.iter().enumerate() {

            // 3D cylinder
            let mut cylinder = Array3::<f32>::zeros((NS as usize, NS as usize, ND as usize));

            let mut tmp_cylinder = vec![0.0; CYLINDER_SIZE];
            let mut cyl_values = bitarr![u32, Lsb0; 0; CYLINDER_SIZE];
            let mut cyl_validity = !bitarr![u32, Lsb0; 0; CYLINDER_SIZE];
            // let mut cyl_values = bitarr![u64, Lsb0; 0; CYLINDER_SIZE];
            // let mut cyl_validity = !bitarr![u64, Lsb0; 0; CYLINDER_SIZE];
            
            // precalc sin & cos to avoid repeated calls
            let sin = sinf(m.arc());
            let cos = cosf(m.arc());

            // invalidating checks
            //let mut contributing_minutiae = 0;
            let mut contributing_minutiae_vec = vec![0; minutiae.len()];
            let mut invalid_cells: f32 = 0.0;

            // for each cell
            for i in 0..NS as usize {
                for j in 0..NS as usize {

                    // cell center
                    let (cx, cy) = cell_center(m, i+1, j+1, sin ,cos);

                    
                    // if cell center is farther than the radius or outside de convexhull, invalidate cell
                    if ed_to_point(m, cx, cy) > CYLINDER_RADIUS || !is_point_inside_chull(&chull,cx, cy, OMEGA) {
                        
                        // invalidate the cell on all K levels
                        for k_val in 0..ND as usize {
                            cyl_validity.set(linearize_idxs(i, j, k_val)-1, false);
                            cylinder[[i, j, k_val]] = -1.0;
                        }

                        invalid_cells += 1.0;
                        continue;
                    }

                    // given a cell, we need to calculate 
                    // the contributions of every other minutiae mt
                    // caso contrário, calculamos o valor checamos a contribuicao de cada outra minucia
                    for (idx_cand, mt) in minutiae.iter().enumerate() {

                        // skip itself
                        if idx_probe == idx_cand {
                            continue;
                        }

                        // distance mt -> cell
                        let dist = ed_to_point(mt, cx, cy);

                        // only the closer minutiae contribute
                        if dist > 3.0 * SIGMA_S {
                            continue;
                        }

                        let spatial_contribution = gaussian(dist);
                        //contributing_minutiae += 1;
                        contributing_minutiae_vec[idx_cand] = 1;

                        // Lastly, for each K level on the cylinder
                        for k in 0..ND as usize {
                        
                            // pega o d_phi pre calculado
                            let d_phi = phi_vals[k];
    
                            let ang_diff = ang_diff(d_phi, ang_diff(m.arc(), mt.arc()));
                            let directional_contribution = area_under_gaussian(ang_diff as f64) as f32;
                            
                            tmp_cylinder[linearize_idxs(i, j, k)-1] += spatial_contribution * directional_contribution;
                            cylinder[[i as usize, j as usize, k as usize]] += spatial_contribution * directional_contribution;
                        } // k
                    } // mt
                } // j
            } // i

            //println!("{} {}", contributing_minutiae, contributing_minutiae_2.iter().sum::<i32>());

            // if cylinder area is too small or too few minutiae contribute
            // we do not add it to the final cylinder set
            if 1.0 - invalid_cells / (NS * NS) < MIN_VALID_CELLS || contributing_minutiae_vec.iter().sum::<i32>() < MIN_CONTRIBUTING_MINUTIAE {
                continue;
            }

            //println!("{:?}", tmp_cylinder.len());

            // get the linearized bit cylinder
            for i in 0..NS as usize {
                for j in 0..NS as usize {
                    for k in 0..ND as usize {
                        if cylinder[[i as usize, j as usize, k as usize]] > MICRO_PSI {
                            cyl_values.set(linearize_idxs(i, j, k)-1, true);
                        }
                    }
                }
            }

            cylinders.push(cyl_values);
            validities.push(cyl_validity);
            minutiae_reference.push(idx_probe as u8);

        } // m

        BitCylinderSet { 
            cylinders,
            validities,
            minutiae_reference,
        }
    }

    pub fn len(&self) -> usize {
        self.cylinders.len()
    }
}

//
// --- Matching -----------------------------------------------------------------------------------
//


#[inline]
pub fn bitwise_sum(b: &LinearizedBitCylinder) -> f32 {
    b.count_ones() as f32
}

#[inline]
fn bitwise_norm(b: &LinearizedBitCylinder) -> f32 {
    sqrtf(bitwise_sum(b))
}

/// calculate all similarities 
pub fn local_matching(c1: &BitCylinderSet, m1: &[Minutia], c2: &BitCylinderSet, m2: &[Minutia]) -> Vec<OrderedFloat<f32>> {

    let mut similarities = Vec::with_capacity(c1.len() * c2.len()); 

    // allocate memory once, afterall they all have fixed size
    let mut validity_ab: LinearizedBitCylinder;
    let mut c_ab: LinearizedBitCylinder;
    let mut c_ba: LinearizedBitCylinder;
    let mut xor: LinearizedBitCylinder;

    for i in 0..c1.len() { //ca = c1[i]
        for j in 0..c2.len() { //cb = c2[j]

            // CONDITION 1: angle difference
            if ang_diff(m1[c1.minutiae_reference[i] as usize].arc(), m2[c2.minutiae_reference[j] as usize].arc()) > MAX_ANG_DIFF {
                similarities.push(OrderedFloat(0.0));
                continue;
            }

            // bitwise operations
            validity_ab = c1.validities[i] & c2.validities[j];
            c_ab = c1.cylinders[i] & validity_ab;
            c_ba = c2.cylinders[j] & validity_ab;

            // CONDITION 2: number of matchable cells in cylinder
            //if OrderedFloat(bitwise_sum(&validity_ab) / CYLINDER_SIZE as f32) < OrderedFloat(0.6) {
            if 100 * bitwise_sum(&validity_ab) as usize / CYLINDER_SIZE < 60 {
                similarities.push(OrderedFloat(0.0));
                continue;
            }

            let norm_c_ab = bitwise_norm(&c_ab);
            let norm_c_ba = bitwise_norm(&c_ba);

            // CONDITION 3: norm sum != 0
            if norm_c_ab + norm_c_ba == 0.0 {
                similarities.push(OrderedFloat(0.0));
                continue;
            }

            xor = c_ab ^ c_ba;
            let lambda = 1.0 - ( bitwise_norm(&xor) / (norm_c_ab + norm_c_ba));
            similarities.push(OrderedFloat(lambda));
        }
    }

    similarities
}

pub fn local_matching_lsa(c1: &BitCylinderSet, m1: &[Minutia], c2: &BitCylinderSet, m2: &[Minutia]) -> Vec<OrderedFloat<f32>> {
    
    // The kuhn-munkres solver requires the smaller set first,
    // so before applying the local matching we get the correct
    // order by size
    let smaller_cyl_ref: &BitCylinderSet;
    let smaller_min_ref: &[Minutia];
    let larger_cyl_ref: &BitCylinderSet;
    let larger_min_ref: &[Minutia];

    if c1.len() < c2.len() {
        smaller_cyl_ref = c1;
        smaller_min_ref = m1;
        larger_cyl_ref = c2;
        larger_min_ref = m2;
    }
    else {
        larger_cyl_ref = c1;
        larger_min_ref = m1;
        smaller_cyl_ref = c2;
        smaller_min_ref = m2;
    }

    let similarities = local_matching(smaller_cyl_ref, smaller_min_ref, larger_cyl_ref, larger_min_ref);
    let similarity_matrix = Matrix::from_vec(smaller_cyl_ref.len(), larger_cyl_ref.len(), similarities).unwrap();

    // solving the assignment problem
    let (_, mapping) = kuhn_munkres(&similarity_matrix);

    // save the top similarities found using the kuhn-munkres algorithm
    let mut top_similarities = Vec::new();
    for (idx_c1, idx_c2) in mapping.iter().enumerate() {
        top_similarities.push(similarity_matrix[(idx_c1, *idx_c2)]);
    }

    top_similarities
}

/// calculate np - np best similarities
pub fn calc_np(c1: &BitCylinderSet, c2: &BitCylinderSet) -> usize {
    (NP_MIN + (NP_MAX - NP_MIN) * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as usize
}

/// global score
/// (23) on paper
pub fn global_matching(similarities: &mut [OrderedFloat<f32>], np: usize) -> f32 {

    similarities.sort_by_key(|v| std::cmp::Reverse(*v));
    let sum: OrderedFloat<f32> = similarities[0..np].iter().sum(); //similarities.truncate(np as usize);
    
    sum.0 / np as f32
}

//
// --- Tests -------------------------------------------------------------------------------------- 
//
#[cfg(test)]
mod tests {
    use geo::EuclideanDistance;

    use super::*;

    // #[test]
    // fn t_cell_center() {
    //     let m1 = Minutia {
    //         position: Point { x: 100, y: 100 },
    //         direction: 0,
    //         index: 0,
    //     };

    //     assert_eq!(cell_center(&m1, 1, 1), (37, 37));
    //     assert_eq!(cell_center(&m1, 10, 10), (163, 163));
    // }
}

