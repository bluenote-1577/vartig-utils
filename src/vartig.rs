use std::collections::HashMap;
use rust_htslib::{bcf, bcf::Read as DUMMY_NAME2};
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

pub fn get_vartigs_from_file(file_name: &str, hapq_cutoff: usize, vcf: Option<&String>, other_vcf: Option<&String>) -> Vec<Vartig> {
    let mut toret = vec![];
    let file = File::open(file_name).unwrap();
    let reader = BufReader::new(file);
    let hapq_cutoff = hapq_cutoff as f64;
    let vcf_snppos_to_allele;
    let vcf_snppos_to_pos;
    if vcf.is_some(){
        (vcf_snppos_to_allele, vcf_snppos_to_pos) = get_vcf_profile(vcf.unwrap())
    }
    else{
        vcf_snppos_to_allele = HashMap::default();
        vcf_snppos_to_pos = HashMap::default();
    }

    let other_vcf_snppos_to_allele;
    let other_vcf_snppos_to_pos;

    if other_vcf.is_some(){
        (other_vcf_snppos_to_allele, other_vcf_snppos_to_pos) = get_vcf_profile(other_vcf.unwrap())
    }
    else{
        other_vcf_snppos_to_allele = HashMap::default();
        other_vcf_snppos_to_pos = HashMap::default();

    }

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
                if vartig.hapq.unwrap() < hapq_cutoff{
                    continue
                }
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
                    //byte to int by -48
                    if vcf.is_none(){
                        last_vartig
                            .allele_vec
                            .insert(last_vartig.snprange.0 + i, (c - 48) as i8);
                    }
                    else{
                        let ctg = last_vartig.contig.as_ref().unwrap();
                        let snp_count = last_vartig.snprange.0 + i;
                        let base = vcf_snppos_to_pos[ctg][&snp_count];
                        last_vartig
                            .allele_vec
                            .insert(base, vcf_snppos_to_allele[ctg][&base][(c - 48) as usize]);

                    }
                }
            }
        }
    }

    if toret.len() == 0{
        return toret
    }

    if other_vcf.is_some(){
        let mut bp_allele_vec = other_vcf_snppos_to_allele[&toret.first().unwrap().contig.as_ref().unwrap().clone()].iter().collect::<Vec<_>>();
        bp_allele_vec.sort();
        let bp_vec = bp_allele_vec.iter().map(|x| *x.0).collect::<Vec<usize>>();
        let allele_vec = bp_allele_vec.iter().map(|x| x.1[0]).collect::<Vec<i8>>();
        for vartig in toret.iter_mut(){
            merge_vcfs_vartig(vartig, &bp_vec, &allele_vec)
        }

    }

    return toret
}


fn merge_vcfs_vartig(vartig: &mut Vartig, bp_vec: &Vec<usize>, allele_vec: &Vec<i8>){
    let start = vartig.baserange.0;
    let end = vartig.baserange.1;

    let lind;
    let rind;
    match bp_vec.binary_search(&start){
        Ok(v)=>{lind = v}
        Err(v)=>{lind = v}
    }
    match bp_vec.binary_search(&end){
        Ok(v)=>{rind = v}
        Err(v)=>{rind = v}
    }
    //dbg!(lind,rind);

    for i in lind..rind{
        let bp = bp_vec[i];
        if !vartig.allele_vec.contains_key(&bp){
            vartig.allele_vec.insert(bp_vec[i], allele_vec[i]);
         }
    }

}


pub fn get_vcf_profile<'a>(vcf_file: &str) -> (HashMap<String,HashMap<usize, Vec<i8>>>, HashMap<String, HashMap<usize, usize>>){
    let mut vcf = match bcf::Reader::from_path(vcf_file) {
        Ok(vcf) => vcf,
        Err(_) =>{ panic!("rust_htslib had an error while reading the VCF file. Exiting.");},
    };
    let mut vcf_pos_allele_map = HashMap::default();
    let mut vcf_snp_to_pos_map = HashMap::default();
    let vcf_header = vcf.header().clone();

    let mut last_ref_chrom = String::default();
    let mut snp_count = 0;
    for rec in vcf.records() {
        snp_count += 1;
        let unr = rec.unwrap();
        let alleles = unr.alleles();
        let mut is_snp = true;

        let record_rid = unr.rid().unwrap();
        let contig_name =
            String::from_utf8(vcf_header.rid2name(record_rid).unwrap().to_vec()).unwrap();
        let pos = unr.pos();

        if last_ref_chrom != contig_name{
            last_ref_chrom = contig_name.clone();
            snp_count = 1;
        }
        let pos_allele_map = vcf_pos_allele_map
            .entry(contig_name.clone())
            .or_insert(HashMap::default());
        let snp_pos_map = vcf_snp_to_pos_map
            .entry(contig_name.clone())
            .or_insert(HashMap::default());
        snp_pos_map.insert(snp_count,pos as usize);

        let mut al_vec = vec![];

        for allele in alleles.iter() {
            if allele.len() > 1 {
                is_snp = false;
                break;
            }
            al_vec.push(allele[0] as i8);
        }

        if !is_snp {
            continue;
        }
        pos_allele_map.insert(pos as usize,al_vec);
    }

    return (vcf_pos_allele_map, vcf_snp_to_pos_map);
}

