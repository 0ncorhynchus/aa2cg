extern crate getopts;
extern crate rustc_serialize;
extern crate csv;
extern crate biost;
extern crate groio;

use getopts::Options;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::env;
use biost::Vector3d;
use groio::*;

#[derive(RustcDecodable)]
struct AACGItem {
    aa: String,
    cg: String,
    mass: f32,
}

type AACGMap = HashMap<String, Vec<(String, f32)>>;

fn read_mapping(filename: &str) -> csv::Result<AACGMap> {
    let mut reader = csv::Reader::from_file(filename)?
                                 .has_headers(true);
    let mut map = AACGMap::new();
    for record in reader.decode() {
        let item: AACGItem = record?;
        let list = map.entry(item.cg.trim().to_string()).or_insert(Vec::new());
        list.push((item.aa.trim().to_string(), item.mass));
    }
    Ok(map)
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options] CSV", program);
    print!("{}", opts.usage(&brief));
}

fn main() {
    let args: Vec<_> = env::args().collect();
    let program = &args[0];

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("i", "", "set input file name", "FILE");
    opts.optopt("o", "", "set output file name", "FILE");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            println!("{}", f.to_string());
            print_usage(program, opts);
            return;
        }
    };

    if matches.opt_present("h") {
        print_usage(program, opts);
        return;
    }

    let mapping = if !matches.free.is_empty() {
        &matches.free[0]
    } else {
        print_usage(program, opts);
        return;
    };

    let aacgmap = match read_mapping(mapping) {
        Ok(map) => map,
        Err(err) => {
            println!("{}", err);
            return;
        }
    };

    let input: Structure = match matches.opt_str("i") {
        Some(file) => {
            let mut file = File::open(&file).unwrap();
            let mut s = String::new();
            file.read_to_string(&mut s).unwrap();
            s.parse().unwrap()
        },
        None => {
            return;
        }
    };

    let mut atoms = Vec::new();
    let mut atom_number = 1;
    for residue in input.residues() {
        for (cg, aalist) in aacgmap.iter() {
            let mut total = 0.0_f32;
            let mut position = Vector3d::new(0.0, 0.0, 0.0);
            let mut velocity = Vector3d::new(0.0, 0.0, 0.0);

            for &(ref aa, ref mass) in aalist {
                total += *mass;
                let coord = residue.atoms.get(aa).unwrap();
                position.x += *mass * coord.position.x;
                position.y += *mass * coord.position.y;
                position.z += *mass * coord.position.z;
                velocity.x += *mass * coord.velocity.x;
                velocity.y += *mass * coord.velocity.y;
                velocity.z += *mass * coord.velocity.z;
            }

            position.x /= total;
            position.y /= total;
            position.z /= total;
            velocity.x /= total;
            velocity.y /= total;
            velocity.z /= total;

            atoms.push(Atom {
                res_number: residue.number,
                res_name: residue.name.clone(),
                atom_name: cg.to_string(),
                atom_number: atom_number,
                position: position,
                velocity: velocity,
            });
            atom_number += 1;
        }
    }
    println!("{}", Structure::new(input.title().to_string(),
                                  atoms,
                                  input.box_size().clone()));
}
