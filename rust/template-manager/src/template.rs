#![allow(dead_code)] // Muito codigo nao acessado neste momento, supress warnings in this file

use super::fingerprint_base as fpb;
use std::fmt;

use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct Template {
    pub(crate) nia: i32,
    pub(crate) seq: u8,
    pub(crate) fgp: u8,
    pub core: fpb::Point,
    pub minutiae: Vec<fpb::Minutia>,
}

impl Template {
    pub fn new(n: i32, s: u8, f: u8, minutiae: Vec<fpb::Minutia>, n1x: i32, n1y: i32) -> Template {
        Template {
            nia: n,
            seq: s,
            fgp: f,
            core: fpb::Point { x: n1x, y: n1y },
            minutiae,
        }
    }

    pub fn reg(&self) -> (i32, u8, u8) {
        (self.nia, self.seq, self.fgp)
    }

    pub fn new_from_csv(n: i32, s: u8, f: u8, filepath: &str) -> Template {
        let mut vec_minutiae: Vec<fpb::Minutia> = Vec::new();
    
        let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(filepath).unwrap();

        //let mut rdr = Reader::from_path(filepath).unwrap();
        for result in rdr.records() {
            let record = result.unwrap();
    
            vec_minutiae.push(fpb::Minutia {
                position: fpb::Point{
                    x: record[0].parse().unwrap(),
                    y: record[1].parse().unwrap(),
                },
                direction: record[2].parse::<i16>().unwrap(),
                index: 0
            });
        }
        
        Template::new(n, s, f, vec_minutiae, 0, 0)
    }

    pub fn new_from_incits_buffer(n: i32, s: i32, f: i32, templ: &[u8], n1x: i32, n1y: i32) -> Option<Template> {
        let decoded_template = decode_template(templ)?;
        Some(Template::new(n, s as u8, f as u8, decoded_template, n1x, n1y))
    }
}

impl fmt::Debug for Template {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Template")
        .field("nia", &self.nia)
        .field("seq", &self.seq)
        .field("fgp", &self.fgp)
        .field("core", &self.core)
        .field("num_minucias", &self.minutiae.len())
        .finish()
    }
}

//
// fn used to decode an incits378 template buffer
//
#[allow(unused_variables, unused_mut, unused_assignments)]
fn decode_template(buffer: &[u8]) -> Option<Vec<fpb::Minutia>> {
    let mut res = Vec::new();
    let mut pos: usize = 0;

    //Format Identifier
    if buffer[0] != b'F' || buffer[1] != b'M' || buffer[2] != b'R' || buffer[3] != 0 {
        return None
    }

    //Version of the standard
	if buffer[4] != b' ' || buffer[5] != b'2' || buffer[6] != b'0' || buffer[7] != 0 {
        return None
    }

    let mut length_of_record: i32 = 0;

    if buffer[8] == 0 && buffer[9] == 0 {
        length_of_record |= (buffer[10] as i32) << 24;
        length_of_record |= (buffer[11] as i32) << 16;
        length_of_record |= (buffer[12] as i32) << 8;
        length_of_record |= buffer[13] as i32;
        pos = 14;
    }
    else {
        length_of_record |= (buffer[8] as i32) << 8;
        length_of_record |= buffer[9] as i32;
        pos = 10;
    }

    if length_of_record == 32 {
        return None
    }

    pos += 6;

    //Image Size
    let mut image_width: i16 = (buffer[pos] as i16) << 8; 
    pos += 1;
    image_width |= buffer[pos] as i16; 
    pos += 1;

    let mut image_height: i16 = (buffer[pos] as i16) << 8; 
    pos += 1;
    image_height |= buffer[pos] as i16; 
    pos += 1;

    // Image resolution
    let mut horizontal_res: i16 = (buffer[pos] as i16) << 8;
    pos += 1;
    horizontal_res |= buffer[pos] as i16;
    pos += 1;
    horizontal_res = (horizontal_res as f32 * 2.54 + 0.5) as i16;

    let mut vertical_res: i16 = (buffer[pos] as i16) << 8;
    pos += 1;
    vertical_res |= buffer[pos] as i16;
    pos += 1;
    vertical_res = (vertical_res as f32 * 2.54 + 0.5) as i16;

    // Number of finger views and reserved byte
    pos += 2;

    let fgp: u8 = buffer[pos];
    pos += 1;

    let impression_type: u8 = buffer[pos] & 0x0F;
    pos += 1;

    let quality: u8 = buffer[pos];
    pos += 1;

    let num_minutiae = buffer[pos];
    pos += 1;
    res.reserve(num_minutiae as usize);


    for i in 0..num_minutiae as usize {
        let mtype = (buffer[pos] >> 6) & 0x03;

        let mut pos_x: i16 = ((buffer[pos] & 0x3F) as i16) << 8;
        pos_x |= buffer[pos + 1] as i16;

        let mut pos_y: i16 = ((buffer[pos + 2] & 0x3F) as i16) << 8;
        pos_y |= buffer[pos + 3] as i16;

        let ori: i16 = buffer[pos + 4] as i16;

        let quality = buffer[pos + 5];

        res.push(fpb::Minutia {
            position: fpb::Point{
                x: pos_x as i32,
                y: pos_y as i32,
            },
            direction: ori,
            index: 0
        });
        
        pos += 6;
    }

    Some(res)
}