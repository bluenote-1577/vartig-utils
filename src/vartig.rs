use std::collections::HashMap;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{prelude::*, BufReader};

fn get_range(token: &str) -> (usize, usize) {
    let range = token
        .split('-')
        .collect::<Vec<&str>>();
    let r1 = range[0].parse::<usize>().unwrap();
    let r2 = range[1].parse::<usize>().unwrap();
    return (r1,r2);
}

fn get_val(token: &str) -> f64 {
    let val = token;
    let ret = val.parse::<f64>().unwrap();
    return ret
}

fn get_pair(token: &str) -> (String, String){
    let spl = token.split(':').collect::<Vec<&str>>();
    let val = (spl[0].to_string(), spl[1].to_string());
    val
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

            let mut tag_dict:HashMap<_,_> = HashMap::default();
            for thing in spl_vec{
                if thing.contains(':'){
                    let (key, val) = get_pair(thing);
                    tag_dict.insert(key, val);
                }
            }
            
            let snprange = tag_dict.get("SNPRANGE").unwrap();
            vartig.snprange = get_range(snprange);
            vartig.baserange = get_range(tag_dict.get("BASERANGE").unwrap());
            let contig = tag_dict.get("CONTIG");
            let cov = tag_dict.get("COV");
            let err = tag_dict.get("ERR");
            let hapq = tag_dict.get("HAPQ");

            if contig.is_none(){
                vartig.contig = None;
            }
            else{
                vartig.contig = Some(contig.unwrap().clone());
            }

            if cov.is_none() {
                vartig.cov = None;
            } else {
                vartig.cov = Some(get_val(cov.unwrap()));
            }
            if err.is_none() {
                vartig.err= None;
            } else {
                vartig.err= Some(get_val(err.unwrap()));
            }
            if hapq.is_none() {
                vartig.hapq = None;
            } else {
                vartig.hapq = Some(get_val(hapq.unwrap()));
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
