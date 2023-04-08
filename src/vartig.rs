use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{prelude::*, BufReader};

fn get_range(token: &str) -> (usize, usize) {
    let range = token.split(':').collect::<Vec<&str>>()[1]
        .split('-')
        .collect::<Vec<&str>>();
    let r1 = range[0].parse::<usize>().unwrap();
    let r2 = range[1].parse::<usize>().unwrap();
    return (r1, r2);
}

fn get_val(token: &str) -> f64 {
    let val = token.split(':').collect::<Vec<&str>>()[1];
    let ret = val.parse::<f64>().unwrap();
    return ret;
}

#[derive(Default, Debug, Clone, PartialEq)]
pub struct VartigAln {
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
    pub cov1: Option<f64>,
    pub cov2: Option<f64>,
    pub diff: f64,
    pub same: f64,
}

#[derive(Default, Debug, Clone, PartialEq)]
pub struct Vartig {
    //closed interval
    pub name: String,
    pub contig: Option<String>,
    pub index: usize,
    pub snprange: (usize, usize),
    pub allele_vec: HashMap<usize, i8>,
    pub baserange: (usize, usize),
    pub err: Option<f64>,
    pub cov: Option<f64>,
    pub hapq: Option<f64>,
}

pub fn get_vartigs_from_file(file_name: &str) -> Vec<Vartig> {
    let mut toret = vec![];
    let file = File::open(file_name).unwrap();
    let reader = BufReader::new(file);
    let mut hapq_pass = false;
    let hapq_cutoff = 1.;

    for line in reader.lines() {
        let l = line.unwrap();
        if l.chars().nth(0).unwrap() == '>' {
            let spl = l.split_whitespace();
            let mut vartig = Vartig::default();
            let spl_vec = spl.collect::<Vec<&str>>();
            let name = &spl_vec[0][1..];
            //            let tmp = col[0..col.len()-1].join("_").to_string();
            //            let remove = tmp.split('/').collect::<Vec<&str>>();
            vartig.name = name.to_string();

            let re_snprange = Regex::new(r"SNPRANGE:(\d+)-(\d+)").unwrap();
            let re_contig = Regex::new(r"[ \t]CONTIG:([^\s]+)").unwrap();
            let re_baserange = Regex::new(r"BASERANGE:(\d+)-(\d+)").unwrap();
            let re_cov = Regex::new(r":COV:(/\d+\.?\d*/)").unwrap();
            let re_err = Regex::new(r"ERR:(/\d+\.?\d*/)").unwrap();
            let re_hapq = Regex::new(r"HAPQ:(\d+)").unwrap();

            let sr_cap = re_snprange.captures(&l).unwrap();
            let contig_cap = re_contig.captures(&l);
            let br_cap = re_baserange.captures(&l);
            let cov_cap = re_cov.captures(&l);
            let err_cap = re_err.captures(&l);
            let hapq_cap = re_hapq.captures(&l);

            vartig.snprange = (
                sr_cap.get(1).unwrap().as_str().parse().unwrap(),
                sr_cap.get(2).unwrap().as_str().parse().unwrap(),
            );
            let br_cap = br_cap.unwrap();
            vartig.baserange = (
                    br_cap.get(1).unwrap().as_str().parse().unwrap(),
                    br_cap.get(2).unwrap().as_str().parse().unwrap(),
                );
            if contig_cap.is_none(){
                vartig.contig = None;
            }
            else{
                vartig.contig = Some(contig_cap.unwrap().get(1).unwrap().as_str().to_string());
            }
            if cov_cap.is_none() {
                vartig.cov = None;
            } else {
                vartig.cov = Some(cov_cap.unwrap().get(1).unwrap().as_str().parse().unwrap());
            }
            if err_cap.is_none() {
                vartig.err= None;
            } else {
                vartig.err= Some(err_cap.unwrap().get(1).unwrap().as_str().parse().unwrap());
            }
            if hapq_cap.is_none() {
                vartig.hapq = None;
            } else {
                vartig.hapq = Some(hapq_cap.unwrap().get(1).unwrap().as_str().parse().unwrap());
            }
            vartig.index = toret.len();
            toret.push(vartig)
        } else {
            let bytes = l.as_bytes();
            let last_vartig = toret.last_mut().unwrap();
            for (i, c) in bytes.iter().enumerate() {
                //question mark
                if *c == 63 {
                    //                   last_vartig.allele_vec.insert(last_vartig.snprange.0 + i,-1);
                } else {
                    last_vartig
                        .allele_vec
                        .insert(last_vartig.snprange.0 + i, (c - 48) as i8);
                }
            }
        }
    }
    toret
}
