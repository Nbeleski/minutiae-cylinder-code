mod constants;
mod utils;
mod mcc;
pub mod mcc_bit;


use template_manager::fingerprint_base as fpb;
use template_manager::fingerprint_base::Minutia;
use mcc_bit::BitCylinderSet;

// MCC16b

pub fn build_cylinders(probe_minutiae: &[fpb::Minutia]) -> BitCylinderSet {
    BitCylinderSet::new(probe_minutiae)
}

pub fn match_cylinder_to_template(probe_cylinderset: &BitCylinderSet, probe_minutiae: &[Minutia], reference_minutiae: &[Minutia]) -> f32 {
    let reference_cylinderset = BitCylinderSet::new(reference_minutiae);

    let mut km_similarities = mcc_bit::local_matching_lsa(probe_cylinderset, probe_minutiae,
                                                                                &reference_cylinderset, reference_minutiae);

    let np = mcc_bit::calc_np(probe_cylinderset, &reference_cylinderset);
    mcc_bit::global_matching(&mut km_similarities, np)
}

pub fn match_cylinder_to_cylinder(probe_cylinderset: &BitCylinderSet, probe_minutiae: &[Minutia], reference_cylinderset: &BitCylinderSet, reference_minutiae: &[Minutia]) -> f32 {
    let mut km_similarities = mcc_bit::local_matching_lsa(probe_cylinderset, probe_minutiae,
                                                                                reference_cylinderset, reference_minutiae);

    let np = mcc_bit::calc_np(probe_cylinderset, reference_cylinderset);
    mcc_bit::global_matching(&mut km_similarities, np)
}

// ---

// --- Tests ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::mcc::{match_bitwise_cylinder_sets_lss, match_bitwise_cylinder_sets_lss_r, bitwise_norm, CYLINDER_SIZE};
    use template_manager::template::Template;
    use bitvec::prelude::*;

    use super::*;


    #[test]
    fn t_bitvec() {
        let mut b1 = bitarr![u32, Lsb0; 0; CYLINDER_SIZE];
        //let mut b2 = bitarr![u32, Lsb0; 0; CYLINDER_SIZE];
        
        b1.set(0, true);
        b1.set(1, true);
        b1.set(2, true);
        
        println!("bitwise_norm: {}", bitwise_norm(&b1));
        println!("mem::size_of_val: {}", std::mem::size_of_val(&b1)); // 192b
        //println!("b1: {:?}", b1);
        
    }

    #[test]
    fn t_1x1_bitwise_cmp() {

        let buffer_t1 = std::fs::read("templates/1023_1_1.incits").unwrap();
        let t1 = Template::new_from_incits_buffer(1, 1, 1, &buffer_t1, 0, 0).unwrap();
        
        let buffer_t2 = std::fs::read("templates/1023_2_1.incits").unwrap();
        let t2 = Template::new_from_incits_buffer(1, 1, 1, &buffer_t2, 0, 0).unwrap();

        let buffer_t3 = std::fs::read("templates/1072_5_1.incits").unwrap();
        let t3 = Template::new_from_incits_buffer(1, 1, 1, &buffer_t3, 0, 0).unwrap();

        let c1 = mcc_bit::BitCylinderSet::new(&t1.minutiae);
        let c2 = mcc_bit::BitCylinderSet::new(&t2.minutiae);
        let c3 = mcc_bit::BitCylinderSet::new(&t3.minutiae);

        println!("t1 -> num_minutia: {}\tnum_valid_cyl: {}", t1.minutiae.len(), c1.cylinders.len());
        println!("t2 -> num_minutia: {}\tnum_valid_cyl: {}", t2.minutiae.len(), c2.cylinders.len());
        println!("t3 -> num_minutia: {}\tnum_valid_cyl: {}", t3.minutiae.len(), c3.cylinders.len());

        let mut similarities = mcc_bit::local_matching(&c1, &t1.minutiae, &c2, &t2.minutiae);
        let mut km_similarities = mcc_bit::local_matching_lsa(&c1, &t1.minutiae, &c2, &t2.minutiae);
        let np = mcc_bit::calc_np(&c1, &c2);
        let score_lss = mcc_bit::global_matching(&mut similarities, np);
        let score_lsa = mcc_bit::global_matching(&mut km_similarities, np);

        println!("genuine score lss: {}", score_lss);
        println!("genuine score lsa: {}", score_lsa);

        similarities = mcc_bit::local_matching(&c1, &t1.minutiae, &c3, &t3.minutiae);
        km_similarities = mcc_bit::local_matching_lsa(&c1, &t1.minutiae, &c3, &t3.minutiae);
        let np = mcc_bit::calc_np(&c1, &c3);
        let score_lss = mcc_bit::global_matching(&mut similarities, np);
        let score_lsa = mcc_bit::global_matching(&mut km_similarities, np);

        println!("imposter score lss: {}", score_lss);
        println!("imposter score lsa: {}", score_lsa);

    }

    #[test]
    fn t_1x1_bitwise() {

        let buffer_t1 = std::fs::read("templates/1023_1_1.incits").unwrap();
        let t1 = Template::new_from_incits_buffer(1, 1, 1, &buffer_t1, 0, 0).unwrap();
        
        let buffer_t2 = std::fs::read("templates/1023_2_1.incits").unwrap();
        let t2 = Template::new_from_incits_buffer(1, 1, 1, &buffer_t2, 0, 0).unwrap();

        let (c1, v1, o1, i1) = mcc::calc_bitwise_cylinder_code(&t1.minutiae);
        let (c2, v2, o2, i2) = mcc::calc_bitwise_cylinder_code(&t2.minutiae);

        println!("t1 -> num_minutia: {}\tnum_valid_cyl: {}", t1.minutiae.len(), c1.len());
        println!("t2 -> num_minutia: {}\tnum_valid_cyl: {}", t2.minutiae.len(), c2.len());

        println!("score LSS:   {}", match_bitwise_cylinder_sets_lss(&c1, &v1, &o1, &c2, &v2, &o2));
        println!("score LSS-R: {}", match_bitwise_cylinder_sets_lss_r(&c1, &v1, &t1.minutiae, &i1, &c2, &v2, &t2.minutiae, &i2));
    }
}