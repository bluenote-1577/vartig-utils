use std::fs::File;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::io::{prelude::*, BufReader};


fn get_range(token: &str) -> (usize, usize){
        let range = token.split(':').collect::<Vec<&str>>()[1].split('-').collect::<Vec<&str>>();
        let r1 = range[0].parse::<usize>().unwrap();
        let r2 = range[1].parse::<usize>().unwrap();
        return (r1,r2);
}

fn get_val(token: &str) -> f64{
        let val = token.split(':').collect::<Vec<&str>>()[1];
        let ret =  val.parse::<f64>().unwrap();
        return ret;
}

#[derive(Default, Debug, Clone, PartialEq)]
pub struct VartigAln{
    pub name1: String,
    pub name2: String,
    pub vtig1_len: usize,
    pub vtig2_len: usize,
    pub snp_range1: (usize, usize),
    pub snp_range2: (usize, usize),
    pub base_range1: (usize, usize),
    pub base_range2: (usize, usize),
    pub snp_identity: f64,
    pub gapless_base_identity: f64,
    pub cov1: f64,
    pub cov2: f64,
    pub diff: f64,
    pub same: f64,
}

#[derive(Default, Debug, Clone, PartialEq)]
pub struct Vartig {
    //closed interval
    pub name: String,
    pub index: usize,
    pub baserange: (usize, usize),
    pub snprange: (usize, usize),
    pub allele_vec: HashMap<usize,i8>,
    pub err: f64,
    pub cov: f64,
    pub hapq: f64
}

pub fn get_vartigs_from_file(file_name: &str) -> Vec<Vartig> {
    let mut toret = vec![];
    let file = File::open(file_name).unwrap();
    let reader = BufReader::new(file);
    let mut hapq_pass = false;
    let hapq_cutoff = 1.;

    for line in reader.lines() {
        let l = line.unwrap();
        if l.chars().nth(0).unwrap() == '>'{
            let spl = l.split_whitespace();
            let mut vartig = Vartig::default();
            let spl_vec = spl.collect::<Vec<&str>>();
            let col = spl_vec[0][1..].split('_').collect::<Vec<&str>>();



            let tmp = col[0..col.len()-1].join("_").to_string();
            let remove = tmp.split('/').collect::<Vec<&str>>();
            vartig.name = remove[remove.len()-1].to_string();
            vartig.snprange = get_range(spl_vec[1]);
            vartig.baserange = get_range(spl_vec[2]);
            vartig.cov = get_val(spl_vec[3]);
            vartig.err = get_val(spl_vec[4]);
            vartig.hapq  = get_val(spl_vec[5]);
            vartig.index = toret.len();
            if !vartig.err.is_nan() && vartig.hapq > hapq_cutoff{
                toret.push(vartig);
                hapq_pass = true;
            }
            else{
                hapq_pass = false;
            }
        }
        else if hapq_pass{
           let bytes = l.as_bytes(); 
           let last_vartig = toret.last_mut().unwrap();
           for (i,c) in bytes.iter().enumerate(){
               //question mark
               if *c == 63 {
//                   last_vartig.allele_vec.insert(last_vartig.snprange.0 + i,-1); 
               }
               else{
                   last_vartig.allele_vec.insert(last_vartig.snprange.0 + i, (c - 48) as i8);
               }
           }
        }
    }
    toret
}
