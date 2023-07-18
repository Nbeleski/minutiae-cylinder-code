#![allow(dead_code, unused_imports)] // Muito codigo nao acessado neste momento, supress warnings in this file

use template_manager::fingerprint_base::{Point, Minutia, PI, TWO_PI};

extern crate ndarray;
use ndarray::{Array2, Array3};
use libm::{sqrtf, cosf, sinf, expf, erf, sqrt, atan2f};
use geo::{LineString, Polygon, ConvexHull, EuclideanDistance};
use::ordered_float::OrderedFloat;
use pathfinding::prelude::{kuhn_munkres, Matrix};
use bitvec::prelude::*;

const CYLINDER_RADIUS: f32 = 70.0;
pub const NS: f32 = 16.0;
pub const ND: f32 = 6.0;
const DELTA_S: f32 = 2.0 * CYLINDER_RADIUS / NS;
const DELTA_D: f32 = 2.0 * PI / ND;

const SIGMA_S: f32 = 28.0 / 3.0;
const SIGMA_D: f32 = 2.0 * PI / 9.0;
const THREE_SIGMA_S: f32 = 3.0 * SIGMA_S;
const SQRT_TWO_PI: f32 = 2.5066283;

const OMEGA: f64 = 50.0;

const MIN_VALID_CELLS: f32 = 0.75;  // % minimo de celulas para um cilindro ser valido
const MIN_MINUTIAE: i32 = 2;        // numero minimo de minucias contribuintes para um cilindro valido
const MAX_ANG_DIFF: f32 = PI / 2.0; // angulo maximo na orientação da minucia/cilindro

pub const CYLINDER_SIZE: usize = (NS * NS * ND) as usize;
pub type LinearizedCylinder = [f32; 1 + (NS * NS * ND) as usize];

// bitwise params and types
const MICRO_PSI: f32 = 0.01;
pub type BitWiseCylinder = BitArray<[u32; CYLINDER_SIZE / 32]>; // CYLINDER_SIZE + 1 ?

// relaxation parameters
const WR: OrderedFloat<f32> = OrderedFloat(0.5);
const RELAXING_ITERATIONS: usize = 5;

//
// 1 <= i, j <= NS 
//
#[inline(always)] //always? never?
fn cell_center(m: &Minutia, i: usize, j: usize) -> (i32, i32) {

    // TODO: precalcular sin e cos e passar como argumentos...
    let cos = cosf(m.arc());
    let sin = sinf(m.arc());
    let i_term = i as f32 - (NS + 1.0) / 2.0;
    let j_term = j as f32 - (NS + 1.0) / 2.0;

    let x: f32 = 0.5 + m.position.x as f32 + DELTA_S * (cos * i_term + sin * j_term);
    let y: f32 = 0.5 + m.position.y as f32 + DELTA_S * (-sin * i_term + cos * j_term);

    (x as i32, y as i32)
}

// Mesmo que acima, porem com sin e cos precalculados
#[inline(always)] //always? never?
fn cell_center_2(m: &Minutia, i: usize, j: usize, sin: f32, cos: f32) -> (i32, i32) {

    let i_term = i as f32 - (NS + 1.0) / 2.0;
    let j_term = j as f32 - (NS + 1.0) / 2.0;

    let x: f32 = 0.5 + m.position.x as f32 + DELTA_S * (cos * i_term + sin * j_term);
    let y: f32 = 0.5 + m.position.y as f32 + DELTA_S * (-sin * i_term + cos * j_term);

    (x as i32, y as i32)
}

//
// Angulo associado a uma cell
//
#[inline(always)]
fn cell_angle(k: i32) -> f32 {
    - PI + (k as f32 - 0.5) * DELTA_D
}


#[inline(always)]
fn ed_to_point(m: &Minutia, x: i32, y: i32) -> f32 {
    let dx: f32 = (m.position.x - x) as f32;
    let dy: f32 = (m.position.y - y) as f32;
    
    sqrtf(dx*dx+dy*dy)
}

#[inline(always)]
pub fn ed(pi: &Minutia, pj: &Minutia) -> f32 {
    let dx: f32 = (pi.position.x - pj.position.x) as f32;
    let dy: f32 = (pi.position.y - pj.position.y) as f32;
    
    sqrtf(dx*dx+dy*dy)
}

#[inline(always)]
fn gaussian(t: f32) -> f32 {
    (1.0 / (SIGMA_S * SQRT_TWO_PI)) * expf(-(t*t) / (2.0 * SIGMA_S * SIGMA_S))
}

#[inline(always)]
fn ang_diff(t1: f32, t2: f32) -> f32 {
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

// resolução daquela integral maluca usando wolfram
// https://www.wolframalpha.com/input?i2d=true&i=Divide%5B1%2Csigma+*+Sqrt%5B2+*+pi%5D%5D+*+Integrate%5BPower%5Be%2CDivide%5B-Power%5Bt%2C2%5D%2C2*Power%5Bsigma%2C2%5D%5D%5D%2C%7Bt%2Cx+-+Divide%5BD%2C2%5D%2Cx+%2B+Divide%5BD%2C2%5D%7D%5D
#[inline(always)]
fn area_under_gaussian(v: f64) -> f64 {
    0.5 * (erf((0.5 * DELTA_D as f64 + v) / (sqrt(2.0) * SIGMA_D as f64)) - erf((v - 0.5 * DELTA_D as f64) / (sqrt(2.0) * SIGMA_D as f64)))
}

/// Retorna o fecho convexo de um conjunto de minucias
fn convex_hull(minutiae: &[Minutia]) -> Polygon {
    let mut coords = Vec::new();
    //let mut ls = LineString::new();
    for m in minutiae {
        coords.push((m.position.x as f64, m.position.y as f64));
    }

    let poly = Polygon::new(LineString::from(coords), vec![]);
    
    poly.convex_hull()
}

/// Calcula a distancia de um ponto até um poligono com um offset de aceitação
/// usado para verificar se um centro de celula está (quase) dentro do fecho de minucias
#[inline(always)]
fn is_point_inside_chull(chull: &Polygon, px: i32, py: i32, omega: f64) -> bool {
    let p: geo::Point = (px as f64, py as f64).into();

    if chull.euclidean_distance(&p) > omega {
        return false
    }

    true
}

#[inline(always)]
pub fn sigmoid(v: f32, u: f32, t: f32) -> f32 {
    1.0 / (1.0 + expf(-t * (v - u)))
}


/// Calcula o cylinder set dado um conjunto de minucias
/// Cada cilindro da saída tem tamanho CYLINDER_SIZE + 1 pois a orientação
/// é anexada ao array (facilita na hora de fazer o matching...)
pub fn calc_cylinder_code(minutiae: &[Minutia]) -> (Vec<[f32; 1 + CYLINDER_SIZE]>, Vec<usize>) {

    // output
    let mut cylinder_set = Vec::with_capacity(minutiae.len());
    let mut minutiae_index = Vec::new();

    // calcula o fecho convexo
    let chull = convex_hull(minutiae);

    // precalcula todos os valores de d_phi (pois nao mudam)
    let phi_vals: Vec<f32> = (1..ND as i32 + 1).map(cell_angle).collect();

    // para cada minucia do template
    for (idx_probe, m) in minutiae.iter().enumerate() {
    
        // C_m no artigo, é o cilindro de uma minucia
        let mut cylinder = Array3::<f32>::zeros((NS as usize, NS as usize, ND as usize));
        
        let mut contributing_minutiae = 0;
        
        let mut invalid_cells: f32 = 0.0;

        // Para cada celula do cilindro da minucia m
        for i in 0..NS as usize {
            for j in 0..NS as usize {

                // pega o centro da celula
                let (cx, cy) = cell_center(m, (i+1).try_into().unwrap(), (j+1).try_into().unwrap());

                // se o centro da celula esta fora do cilindro
                // ou fora do convexhull de minucias, 
                if ed_to_point(m, cx, cy) > CYLINDER_RADIUS || !is_point_inside_chull(&chull,cx, cy, OMEGA) {
                    
                    // seta todos os Ks para esse I e J como invalidos
                    for k_val in 0..ND as usize {
                        cylinder[[i, j, k_val]] = -1.0;
                        
                    }

                    // incrementa o numero de cells invalidas
                    invalid_cells += 1.0;
                    
                    continue;
                }

                
                // caso contrário, calculamos o valor checamos a contribuicao de cada outra minucia
                for (idx_cand, mt) in minutiae.iter().enumerate() {

                    // nao precisa calcular a influencia dela mesma no cilindro
                    if idx_probe == idx_cand {
                        continue;
                    }

                    // distancia da celula até a minucia candidata
                    let dist = ed_to_point(mt, cx, cy);

                    // Se a celula esta muito longe da minucia, continue
                    if dist > 3.0 * SIGMA_S {
                        continue
                    }

                    // contrib
                    let spatial_contribution = gaussian(dist);

                    // incrementa o numero de contribs
                    contributing_minutiae += 1;

                    for k in 0..ND as usize {
                        
                        // pega o d_phi pre calculado
                        let d_phi = phi_vals[k];

                        let ang_diff = ang_diff(d_phi, ang_diff(m.arc(), mt.arc()));
                        let directional_contribution = area_under_gaussian(ang_diff as f64) as f32;

                        //cylinder[[i as usize, j as usize, k as usize]] += sigmoid(spatial_contribution * directional_contribution, 0.01, 400.0);
                        cylinder[[i as usize, j as usize, k as usize]] += spatial_contribution * directional_contribution;
                    } // k
                } // mt
            } // j
        } // i

        // Cilindro invalido?
        //println!("invalid? {} {}", 1.0 - invalid_cells / (NS * NS), contributing_minutiae);
        if 1.0 - invalid_cells / (NS * NS) < MIN_VALID_CELLS || contributing_minutiae < 2 {
            continue;
        }

        // lineariza a resposta
        // TODO: como resolver isso direto no loop principal
        let mut linearized_cyl: LinearizedCylinder = [0.0; 1 + CYLINDER_SIZE];
        for i in 0..NS as i32 {
            for j in 0..NS as i32 {
                for k in 0..ND as i32 {
                    let idx = linearize_idxs(i, j, k) as usize;
                    let val = cylinder[(i as usize, j as usize, k as usize)];
                    if val != -1.0 {
                        linearized_cyl[idx-1] = sigmoid(val, 0.01, 400.0);
                    }
                    else {
                        linearized_cyl[idx-1] = val;
                    }
                }
            }
        }

        // coloca o angulo da minucia na ultima posição do cilindro... //TODO: uma maneira mais descritiva de fazer isso?
        linearized_cyl[linearized_cyl.len() - 1] = 2.0 * m.direction as f32 * PI / 180.0;

        cylinder_set.push(linearized_cyl);
        minutiae_index.push(idx_probe);

    } // m

    (cylinder_set, minutiae_index)
}

/// Calcula o cylinder set dado um conjunto de minucias
/// Cada cilindro da saída tem tamanho CYLINDER_SIZE + 1 pois a orientação
/// é anexada ao array (facilita na hora de fazer o matching...)
pub fn calc_bitwise_cylinder_code(minutiae: &[Minutia]) -> (Vec<BitWiseCylinder>, Vec<BitWiseCylinder>, Vec<f32>, Vec<usize>) {

    // output
    // No caso do bitwise:
    // um vec com os valores em bit para cada celula de cada cilindro
    // um vec com a validade de cada celula de cada cilindro
    // um vec com a orientação de cada cilindro
    let mut cylinder_set_values = Vec::with_capacity(minutiae.len());
    let mut cylinder_set_validities = Vec::with_capacity(minutiae.len());
    let mut cylinder_orientations = Vec::with_capacity(minutiae.len());
    let mut minutiae_index = Vec::new();

    // calcula o fecho convexo
    let chull = convex_hull(minutiae);

    // precalcula todos os valores de d_phi (pois nao mudam) //TODO: gerar compile time, se possivel
    let phi_vals: Vec<f32> = (1..ND as i32 + 1).map(cell_angle).collect();

    // para cada minucia do template
    for (idx_probe, m) in minutiae.iter().enumerate() {
    
        // C_m no artigo, é o cilindro de uma minucia
        let mut cylinder = Array3::<f32>::zeros((NS as usize, NS as usize, ND as usize));
        
        let mut contributing_minutiae = 0;
        
        let mut invalid_cells: f32 = 0.0;

        let sin = sinf(m.arc());
        let cos = cosf(m.arc());

        // Para cada celula do cilindro da minucia m
        for i in 0..NS as usize {
            for j in 0..NS as usize {

                // pega o centro da celula
                //let (cx, cy) = cell_center(m, (i+1).try_into().unwrap(), (j+1).try_into().unwrap());
                let (cx, cy) = cell_center_2(m, (i+1).try_into().unwrap(), (j+1).try_into().unwrap(), sin, cos);

                // se o centro da celula esta fora do cilindro
                // ou fora do convexhull de minucias, 
                if ed_to_point(m, cx, cy) > CYLINDER_RADIUS || !is_point_inside_chull(&chull,cx, cy, OMEGA) {
                    
                    // seta todos os Ks para esse I e J como invalidos
                    for k_val in 0..ND as usize {
                        cylinder[[i, j, k_val]] = -1.0;
                        
                    }

                    // incrementa o numero de cells invalidas
                    invalid_cells += 1.0;
                    
                    continue;
                }

                
                // caso contrário, calculamos o valor checamos a contribuicao de cada outra minucia
                for (idx_cand, mt) in minutiae.iter().enumerate() {

                    // nao precisa calcular a influencia dela mesma no cilindro
                    if idx_probe == idx_cand {
                        continue;
                    }

                    // distancia da celula a minucia candidata
                    let dist = ed_to_point(mt, cx, cy);

                    // Se a celula esta muito longe da minucia, continue
                    if dist > 3.0 * SIGMA_S {
                        continue
                    }

                    // contrib
                    let spatial_contribution = gaussian(dist);

                    // incrementa o numero de contribs
                    contributing_minutiae += 1;

                    for k in 0..ND as usize {
                        
                        // pega o d_phi pre calculado
                        let d_phi = phi_vals[k];

                        let ang_diff = ang_diff(d_phi, ang_diff(m.arc(), mt.arc()));
                        let directional_contribution = area_under_gaussian(ang_diff as f64) as f32;

                        cylinder[[i as usize, j as usize, k as usize]] += spatial_contribution * directional_contribution;
                    } // k
                } // mt
            } // j
        } // i

        // Cilindro invalido?
        //println!("invalid? {} {}", 1.0 - invalid_cells / (NS * NS), contributing_minutiae);
        if 1.0 - invalid_cells / (NS * NS) < MIN_VALID_CELLS || contributing_minutiae < 2 {
            continue;
        }

        // Cria um array de bits, o tamanho aqui é o cilindro mais o espaço pra um f32
        // por exemplo 600 + 32 vai gerar um array de tamanho 640 (multiplo de 32)
        // o outro array é a validade de cada cilindro, inicia-se em 1 e marca com 0 cada invalido (-1.0)
        let mut linearized_bit_cyl_values = bitarr![u32, Lsb0; 0; CYLINDER_SIZE];
        let mut linearized_bit_cyl_validity = !bitarr![u32, Lsb0; 0; CYLINDER_SIZE];

        // lineariza a resposta
        // TODO: como resolver isso direto no loop principal
        for i in 0..NS as i32 {
            for j in 0..NS as i32 {
                for k in 0..ND as i32 {
                    if cylinder[(i as usize, j as usize, k as usize)] == -1.0 {
                        let idx = linearize_idxs(i, j, k) as usize;
                        linearized_bit_cyl_validity.set(idx-1, false);
                    }
                    else if cylinder[(i as usize, j as usize, k as usize)] > MICRO_PSI {
                        let idx = linearize_idxs(i, j, k) as usize;
                        linearized_bit_cyl_values.set(idx-1, true);
                    }
                }
            }
        }

        // TODO: adiciona o angulo nos ultimos bytes
        //linearized_cyl[linearized_cyl.len() - 1] = 2.0 * m.direction as f32 * PI / 180.0;

        cylinder_set_values.push(linearized_bit_cyl_values);
        cylinder_set_validities.push(linearized_bit_cyl_validity);
        cylinder_orientations.push(2.0 * m.direction as f32 * PI / 180.0);
        minutiae_index.push(idx_probe);

    } // m

    (cylinder_set_values, cylinder_set_validities, cylinder_orientations, minutiae_index)
}

// --- matching utils ---

// funcoes de relaxameto

// relaxamento com relação a distancia
#[inline(always)]
fn relaxing_d1(a1: &Minutia, a2: &Minutia, b1: &Minutia, b2: &Minutia) -> f32 {
    sigmoid((ed(a1, a2) - ed(b1, b2)).abs(), 5.0, -8.0 / 5.0)
}

// relaxamento com relação ao angulo
#[inline(always)]
fn relaxing_d2(a1: &Minutia, a2: &Minutia, b1: &Minutia, b2: &Minutia) -> f32 {
    sigmoid(ang_diff(ang_diff(a1.arc(), a2.arc()), ang_diff(b1.arc(), b2.arc())).abs(), PI / 12.0, -30.0)
}

// relaxamento com relação ao angulo radial 
#[inline(always)]
fn relaxing_d3(a1: &Minutia, a2: &Minutia, b1: &Minutia, b2: &Minutia) -> f32 {
    let dra = ang_diff(a1.arc(), atan2f((-a2.position.y + a1.position.y) as f32, (a2.position.x - a1.position.x) as f32));
    let drb = ang_diff(b1.arc(), atan2f((-b2.position.y + b1.position.y) as f32, (b2.position.x - b1.position.x) as f32));
    sigmoid(ang_diff(dra, drb).abs(), PI / 12.0, -30.0)
}

fn weighting_factor(a1: &Minutia, a2: &Minutia, b1: &Minutia, b2: &Minutia) -> f32 {
    relaxing_d1(a1, a2, b1, b2) * relaxing_d2(a1, a2, b1, b2) * relaxing_d3(a1, a2, b1, b2)
}


// --- Matching bitwise ---

#[inline(always)]
pub fn bitwise_sum(b: &BitWiseCylinder) -> f32 {
    //TODO: implementacao naive
    let mut sum: f32 = 0.0;
    for bit in b {
        if *bit {
            sum += 1.0;
        }
    }

    sum
}

#[inline(always)]
pub fn bitwise_norm(b: &BitWiseCylinder) -> f32 {
    sqrtf(bitwise_sum(b))
}

fn calc_bitwise_similarities(
    c1: &[BitWiseCylinder], v1: &[BitWiseCylinder], o1: &[f32],
    c2: &[BitWiseCylinder], v2: &[BitWiseCylinder], o2: &[f32]) -> Vec<OrderedFloat<f32>> {

    let mut all_similarities = Vec::new();

    for (i, (ca, va)) in c1.iter().zip(v1.iter()).enumerate() {

        for (j, (cb, vb)) in c2.iter().zip(v2.iter()).enumerate() {

            // CONDICAO 1 - checa o angulo
            if (o1[i] - o2[j]).abs() > MAX_ANG_DIFF {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let vab = *va & *vb;
            let cab = *ca & vab;
            let cba = *cb & vab;

            // CONDICAO 2 - pelo menos 60% matchable
            if OrderedFloat(bitwise_sum(&vab) / CYLINDER_SIZE as f32) < OrderedFloat(0.6) {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            // CONDICAO 3 - != 0
            if bitwise_norm(&cab) + bitwise_norm(&cba) == 0.0 {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let xor = cab ^ cba;
            let similarity = 1.0 - (bitwise_norm(&xor)) / (bitwise_norm(&cab) + bitwise_norm(&cba));
            all_similarities.push(OrderedFloat(similarity));
        }
    }

    all_similarities
}

pub fn match_bitwise_cylinder_sets_lss(
    c1: &[BitWiseCylinder], v1: &[BitWiseCylinder], o1: &[f32],
    c2: &[BitWiseCylinder], v2: &[BitWiseCylinder], o2: &[f32]) -> f32 {

    // calcula-se todas as similaridades N:M
    let mut all_similarities = calc_bitwise_similarities(c1, v1, o1, c2, v2, o2);

    // numero de top pares usados no score global (24) no paper
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    // ordena por escore decrescente, trunca em np e acumula
    all_similarities.sort_by_key(|v| std::cmp::Reverse(*v));
    all_similarities.truncate(np as usize);
    let sum: OrderedFloat<f32> = all_similarities.iter().sum(); // (23) no paper
    
    sum.0 / np as f32
}

pub fn match_bitwise_cylinder_sets_lss_r(
    c1: &[BitWiseCylinder], v1: &[BitWiseCylinder], m1: &[Minutia], i1: &[usize],
    c2: &[BitWiseCylinder], v2: &[BitWiseCylinder], m2: &[Minutia], i2: &[usize]) -> f32 {

    let mut current_similarities = Vec::new();

    // guarda o indice das minucias que formam cada par
    let mut minutiae_pair_idxs = Vec::new();

    // calcula-se todas as similaridades N:M
    for (i, (ca, va)) in c1.iter().zip(v1.iter()).enumerate() {

        for (j, (cb, vb)) in c2.iter().zip(v2.iter()).enumerate() {

            minutiae_pair_idxs.push((i1[i], i2[j]));

            // CONDICAO 1 - checa o angulo
            if (m1[i].arc() - m2[j].arc()).abs() > MAX_ANG_DIFF {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let vab = *va & *vb;
            let cab = *ca & vab;
            let cba = *cb & vab;

            // CONDICAO 2 - pelo menos 60% matchable
            if OrderedFloat(bitwise_sum(&vab) / CYLINDER_SIZE as f32) < OrderedFloat(0.6) {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            // CONDICAO 3 - != 0
            if bitwise_norm(&cab) + bitwise_norm(&cba) == 0.0 {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let xor = cab ^ cba;
            let similarity = 1.0 - (bitwise_norm(&xor)) / (bitwise_norm(&cab) + bitwise_norm(&cba));
            current_similarities.push(OrderedFloat(similarity));
        }
    }

    // conjunto preliminar de similaridades
    let nr = core::cmp::min(c1.len(), c2.len()) as i32;

    // obtem a ordem de pares
    let mut ordered_idxs: Vec<usize> = (0..current_similarities.len()).collect();
    ordered_idxs.sort_by_key(|&i| std::cmp::Reverse(&current_similarities[i]));
    ordered_idxs.truncate(nr as usize);

    let original_similarities = current_similarities.clone();
    for _ in 0..RELAXING_ITERATIONS {

        // é preciso fazer uma copia a cada iteraçao para evitar alteraçoes enquanto o loop acontece
        let previous_similarities = current_similarities.clone();

        for t in 0..nr as usize {
            let mut sum: OrderedFloat<f32> = OrderedFloat(0.0);
            let (a_idx_1, b_idx_1) = minutiae_pair_idxs[ordered_idxs[t]]; // obtem o indice ordenado por similaridade original
            
            //println!("Para o par {} {} com similaridade {}:", a_idx_1, b_idx_1, previous_similarities[idx_vals[t]]);
            for k in 0..nr as usize {
                let (a_idx_2, b_idx_2) = minutiae_pair_idxs[ordered_idxs[k]];
                if k == t {
                    continue
                }

                //println!("   {a_idx_2} {b_idx_2} wf: {}\t  prev: {}\t avg:  {}", weighting_factor(&m1[a_idx_1], &m1[a_idx_2], &m2[b_idx_1], &m2[b_idx_2]), previous_similarities[idx_vals[k]], sum / OrderedFloat(nr as f32 - 1.0));
                sum += OrderedFloat(weighting_factor(&m1[a_idx_1], &m1[a_idx_2], &m2[b_idx_1], &m2[b_idx_2])) * previous_similarities[ordered_idxs[k]];
            }

            //println!("   sum: {} {}", current_similarities[idx_vals[t]], sum / OrderedFloat(nr as f32 - 1.0));
            current_similarities[ordered_idxs[t]] = WR * previous_similarities[ordered_idxs[t]] + (OrderedFloat(1.0) - WR) * sum / OrderedFloat(nr as f32 - 1.0);
        }
    }

    let mut efficiencies = Vec::new(); 
    for i in 0..nr as usize {

        if original_similarities[ordered_idxs[i]] != 0.0 {
            efficiencies.push(current_similarities[ordered_idxs[i]] / original_similarities[ordered_idxs[i]]);
        }
        else {
            efficiencies.push(OrderedFloat(0.0));
        }
        
    }

    // pega os np melhores eficiencias
    let mut final_idxs: Vec<usize> = (0..nr as usize).collect();
    final_idxs.sort_by_key(|&i| std::cmp::Reverse(&efficiencies[i]));
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    let mut sum: f32 = 0.0;
    for i in 0..np as usize {
        sum += *current_similarities[ordered_idxs[final_idxs[i]]];
    }
    
    sum / np as f32
}

// --- Matching f32 ---

/// Calcula todas as similaridades possíveis
fn calc_all_similarities(c1: &[LinearizedCylinder], c2: &[LinearizedCylinder]) -> Vec<OrderedFloat<f32>> {
    // onde guardamos todas as similaridades
    let mut all_similarities = Vec::new();

    for (i, ca) in c1.iter().enumerate() {

        // pula o ultimo valor de ca
        if i == ca.len() - 1 {
            continue;
        }

        for (j, cb) in c2.iter().enumerate() {

            // pula o ultimo valor de cb
            if j == cb.len() - 1 {
                continue;
            }

            // CONDICAO 1 - se a dif de angulo for grande, pula
            if (ca.last().unwrap() - cb.last().unwrap()).abs() > MAX_ANG_DIFF {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mut cab: [f32; CYLINDER_SIZE]     = [0.0; CYLINDER_SIZE];
            let mut cba: [f32; CYLINDER_SIZE]     = [0.0; CYLINDER_SIZE];
            let mut cab_cba: [f32; CYLINDER_SIZE] = [0.0; CYLINDER_SIZE];

            for t in 0..CYLINDER_SIZE {
                cab[t] = if ca[t] != -1.0 && cb[t] != -1.0 { ca[t] } else { 0.0 };
                cba[t] = if ca[t] != -1.0 && cb[t] != -1.0 { cb[t] } else { 0.0 };
                cab_cba[t] = cab[t] - cba[t];
            }

            // CONDICAO 2 - pelo menos 60% matchable
            if OrderedFloat(cba.iter().copied().filter(|x| *x > 0.0).count() as f32 / cba.len() as f32) < OrderedFloat(0.6) {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mod_cab: f32 = cab.iter().copied().map(|x| x * x).sum();
            let mod_cba: f32 = cba.iter().copied().map(|x| x * x).sum();

            // CONDICAO 3 - essa soma deve ser diferente de zero
            if mod_cab + mod_cba == 0.0 {
                all_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mod_cab_cba: f32 = cab_cba.iter().copied().map(|x| x * x).sum();

            let similarity = 1.0 - mod_cab_cba / (mod_cab + mod_cba);
            all_similarities.push(OrderedFloat(similarity));
        }
    }

    all_similarities
}

/// matching de cilindros simplificado
pub fn match_cylinder_sets_lss(c1: &[LinearizedCylinder], c2: &[LinearizedCylinder]) -> f32 {

    // calcula-se todas as similaridades N:M
    let mut all_similarities = calc_all_similarities(c1, c2);

    // numero de top pares usados no score global (24) no paper
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    // ordena por escore decrescente, trunca em np e acumula
    all_similarities.sort_by_key(|v| std::cmp::Reverse(*v));
    all_similarities.truncate(np as usize);
    let sum: OrderedFloat<f32> = all_similarities.iter().sum(); // (23) no paper
    
    sum.0 / np as f32
}

/// Matching de cilindros usando algoritmo hungaro
pub fn match_cylinder_sets_lsa(c1: &[LinearizedCylinder], c2: &[LinearizedCylinder]) -> f32 {

    // é preciso usar os cilindros ordenados por tamanho
    // pois o algoritmo hungaro requer isso
    let smaller_ref: &[LinearizedCylinder];
    let larger_ref: &[LinearizedCylinder];

    if c1.len() < c2.len() {
        smaller_ref = c1;
        larger_ref = c2;
    }
    else {
        larger_ref = c1;
        smaller_ref = c2;
    }

    // calcula-se todas as similaridades N:M
    let all_similarities = calc_all_similarities(smaller_ref, larger_ref);

    // cria-se a matrix para resolver o assingment problem
    let m1 = Matrix::from_vec(smaller_ref.len(), larger_ref.len(), all_similarities).unwrap();
    let (_, mapping) = kuhn_munkres(&m1); 

    // numero de top pares usados no score global (24) no paper
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    // guarda os pares usados no mapeamento gerado pelo algoritmo hungaro
    let mut top_similarities = Vec::new();
    for (idx_c1, idx_c2) in mapping.iter().enumerate() {
        top_similarities.push(m1[(idx_c1, *idx_c2)]);
    }

    // ordena por escore decrescente, trunca em np e acumula
    top_similarities.sort_by_key(|v| std::cmp::Reverse(*v));
    top_similarities.truncate(np as usize);
    let sum: OrderedFloat<f32> = top_similarities.iter().sum(); // (23) no paper
    
    sum.0 / np as f32
}

/// match de cilindros simplificado com relaxamento
pub fn match_cylinder_sets_lss_r(c1: &[LinearizedCylinder], m1: &[Minutia], i1: &[usize], c2: &[LinearizedCylinder], m2: &[Minutia], i2: &[usize]) -> f32 {

    // TODO: eu preciso saber qual minucia gerou qual cilindro

    // calcula-se todas as similaridades N:M
    let mut current_similarities = Vec::new();

    // guarda o indice das minucias que formam cada par
    let mut similarity_pair_idxs = Vec::new();

    for (i, ca) in c1.iter().enumerate() {

        // pula o ultimo valor de ca
        if i == ca.len() - 1 {
            continue;
        }

        for (j, cb) in c2.iter().enumerate() {

            // pula o ultimo valor de cb
            if j == cb.len() - 1 {
                continue;
            }

            similarity_pair_idxs.push((i1[i], i2[j]));

            // CONDICAO 1 - se a dif de angulo for grande, pula
            if (ca.last().unwrap() - cb.last().unwrap()).abs() > MAX_ANG_DIFF {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mut cab: [f32; CYLINDER_SIZE]     = [0.0; CYLINDER_SIZE];
            let mut cba: [f32; CYLINDER_SIZE]     = [0.0; CYLINDER_SIZE];
            let mut cab_cba: [f32; CYLINDER_SIZE] = [0.0; CYLINDER_SIZE];

            for t in 0..CYLINDER_SIZE {
                cab[t] = if ca[t] != -1.0 && cb[t] != -1.0 { ca[t] } else { 0.0 };
                cba[t] = if ca[t] != -1.0 && cb[t] != -1.0 { cb[t] } else { 0.0 };
                cab_cba[t] = cab[t] - cba[t];
            }

            // CONDICAO 2 - pelo menos 60% matchable
            if OrderedFloat(cba.iter().copied().filter(|x| *x > 0.0).count() as f32 / cba.len() as f32) < OrderedFloat(0.6) {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mod_cab: f32 = cab.iter().copied().map(|x| x * x).sum();
            let mod_cba: f32 = cba.iter().copied().map(|x| x * x).sum();

            // CONDICAO 3 - essa soma deve ser diferente de zero
            if mod_cab + mod_cba == 0.0 {
                current_similarities.push(OrderedFloat(0.0));
                continue;
            }

            let mod_cab_cba: f32 = cab_cba.iter().copied().map(|x| x * x).sum();

            let similarity = 1.0 - mod_cab_cba / (mod_cab + mod_cba);
            current_similarities.push(OrderedFloat(similarity));
        }
    }

    // conjunto preliminar de similaridades
    let nr = core::cmp::min(c1.len(), c2.len()) as i32;

    // TODO: eu preciso descobrir qual par de minucias resulta em qual similaridade...
    // para isso 
    // obtem a ordem de pares
    let mut idx_vals: Vec<usize> = (0..current_similarities.len()).collect();
    idx_vals.sort_by_key(|&i| std::cmp::Reverse(&current_similarities[i]));
    idx_vals.truncate(nr as usize);


    let original_similarities = current_similarities.clone();
    for _ in 0..RELAXING_ITERATIONS {

        // é preciso fazer uma copia a cada iteraçao para evitar alteraçoes enquanto o loop acontece
        let previous_similarities = current_similarities.clone();

        for t in 0..nr as usize {
            let mut sum: OrderedFloat<f32> = OrderedFloat(0.0);
            let (a_idx_1, b_idx_1) = similarity_pair_idxs[idx_vals[t]]; // obtem o indice ordenado por similaridade original
            
            for k in 0..nr as usize {
                let (a_idx_2, b_idx_2) = similarity_pair_idxs[idx_vals[k]];
                if k == t {
                    continue
                }

                sum += OrderedFloat(weighting_factor(&m1[a_idx_1], &m1[a_idx_2], &m2[b_idx_1], &m2[b_idx_2])) * previous_similarities[idx_vals[k]];
            }

            current_similarities[idx_vals[t]] = WR * previous_similarities[idx_vals[t]] + (OrderedFloat(1.0) - WR) * sum / OrderedFloat(nr as f32 - 1.0);
        }
    }

    let mut efficiencies = Vec::new(); 
    for i in 0..nr as usize {

        if original_similarities[idx_vals[i]] != 0.0 {
            efficiencies.push(current_similarities[idx_vals[i]] / original_similarities[idx_vals[i]]);
        }
        else {
            efficiencies.push(OrderedFloat(0.0));
        }
    }

    // for ef in &efficiencies {
    //     if *ef > OrderedFloat(1.0) {
    //         println!("Shouldn't this be larger than 1.0 sometimes?");
    //     }
    // }

    // pega os np melhores eficiencias
    let mut final_idxs: Vec<usize> = (0..nr as usize).collect();
    final_idxs.sort_by_key(|&i| std::cmp::Reverse(&efficiencies[i]));
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    let mut sum: f32 = 0.0;
    for i in 0..np as usize {
        sum += *current_similarities[idx_vals[final_idxs[i]]];
    }
    
    sum / np as f32
}

/// Matching de cilindros usando algoritmo hungaro
pub fn match_cylinder_sets_lsa_r(c1: &[LinearizedCylinder], m1: &[Minutia], i1: &[usize], c2: &[LinearizedCylinder], m2: &[Minutia], i2: &[usize]) -> f32 {

    // é preciso usar os cilindros ordenados por tamanho
    // pois o algoritmo hungaro requer isso
    let smaller_cyl: &[LinearizedCylinder];
    let smaller_min: &[Minutia];
    let smaller_idx: &[usize];
    let larger_cyl: &[LinearizedCylinder];
    let larger_min: &[Minutia];
    let larger_idx: &[usize];

    if c1.len() < c2.len() {
        smaller_cyl = c1;
        smaller_min = m1;
        smaller_idx = i1;
        larger_cyl = c2;
        larger_min = m2;
        larger_idx = i2;
    }
    else {
        smaller_cyl = c2;
        smaller_min = m2;
        smaller_idx = i2;
        larger_cyl = c1;
        larger_min = m1;
        larger_idx = i1;
    }

    // calcula-se todas as similaridades N:M
    let all_similarities = calc_all_similarities(smaller_cyl, larger_cyl);

    // cria-se a matrix para resolver o assingment problem
    let similarities_matrix = Matrix::from_vec(smaller_cyl.len(), larger_cyl.len(), all_similarities).unwrap();
    let (_, mapping) = kuhn_munkres(&similarities_matrix); 

    // conjunto preliminar de similaridades
    let nr = core::cmp::min(smaller_cyl.len(), larger_cyl.len()) as i32;

    // guarda os pares usados no mapeamento gerado pelo algoritmo hungaro
    let mut current_similarities = Vec::new();
    for (idx_c1, idx_c2) in mapping.iter().enumerate() {
        current_similarities.push(similarities_matrix[(idx_c1, *idx_c2)]);
    }

    let original_similarities = current_similarities.clone();
    for _ in 0..RELAXING_ITERATIONS {

        // é preciso fazer uma copia a cada iteraçao para evitar alteraçoes enquanto o loop acontece
        let previous_similarities = current_similarities.clone();

        for t in 0..nr as usize {
            let mut sum: OrderedFloat<f32> = OrderedFloat(0.0);
            let (a_idx_1, b_idx_1) = (smaller_idx[t], larger_idx[mapping[t]]); // obtem o indice ordenado por similaridade original
            
            for k in 0..nr as usize {
                let (a_idx_2, b_idx_2) = (smaller_idx[k], larger_idx[mapping[k]]);
                if k == t {
                    continue
                }

                sum += OrderedFloat(weighting_factor(&smaller_min[a_idx_1], &smaller_min[a_idx_2], &larger_min[b_idx_1], &larger_min[b_idx_2])) * previous_similarities[k];
            }

            current_similarities[t] = WR * previous_similarities[t] + (OrderedFloat(1.0) - WR) * sum / OrderedFloat(nr as f32 - 1.0);
        }
    }

    let mut efficiencies = Vec::new(); 
    for i in 0..nr as usize {

        if original_similarities[i] != 0.0 {
            efficiencies.push(current_similarities[i] / original_similarities[i]);
        }
        else {
            efficiencies.push(OrderedFloat(0.0));
        }
        
    }

    // pega os np melhores eficiencias
    let mut final_idxs: Vec<usize> = (0..nr as usize).collect();
    final_idxs.sort_by_key(|&i| std::cmp::Reverse(&efficiencies[i]));
    let np = (4.0 + 8.0 * sigmoid(std::cmp::min(c1.len(), c2.len()) as f32, 20.0, 2.0/5.0) + 0.5) as i32;

    let mut sum: f32 = 0.0;
    for i in 0..np as usize {
        sum += *current_similarities[final_idxs[i]];
    }
    
    sum / np as f32

}

#[inline(always)]
pub fn linearize_idxs(i: i32, j: i32, k: i32) -> i32 {
    (k) * NS as i32 * NS as i32 + (j) * NS as i32 + i + 1
}

// --- testes --- 


#[cfg(test)]
mod tests {
    use geo::EuclideanDistance;

    use super::*;

    #[test]
    fn t_cell_center() {
        let m1 = Minutia {
            position: Point { x: 100, y: 100 },
            direction: 0,
            index: 0,
        };

        assert_eq!(cell_center(&m1, 1, 1), (34, 34));
        assert_eq!(cell_center(&m1, 10, 10), (113, 113));
    }

    #[test]
    fn t_convex_hull() {
        let p1 = (1.0, 1.0);
        let p2 = (1.0, -1.0);
        let p3 = (-1.0, -1.0);
        let p4 = (-1.0, 1.0);

        let coords = vec![p1, p2, p3, p4];

        // fast but non-robust
        let poly = Polygon::new(LineString::from(coords), vec![]);
        let res = poly.convex_hull();

        let p5: geo::Point = (0.5, 0.5).into();
        let p6: geo::Point = (2.0, 2.0).into();

        assert_eq!(res.euclidean_distance(&p5), 0.0);
        assert_eq!(res.euclidean_distance(&p6), 1.4142135623730951);
    }
}