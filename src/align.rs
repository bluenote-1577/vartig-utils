use crate::vartig::*;
use std::collections::HashMap;
use std::collections::HashSet;
use fxhash::FxHashMap;

fn same_diff(v1: &Vartig, v2: &Vartig) -> (f64, f64){
    let mut same = 0.;
    let mut diff = 0.;
    let (smaller,bigger) = if v1.allele_vec.len() < v2.allele_vec.len() { (&v1, &v2)} else{ (&v2, &v1)};
    for (pos,allele) in smaller.allele_vec.iter(){
        if bigger.allele_vec.contains_key(pos){
            if *allele == bigger.allele_vec[pos]{
                same += 1.;
            }
            else{
                diff += 1.;
            }
        }
    }
    return (same,diff)
}

fn hit_score(v1: &Vartig, v2: &Vartig) -> f64{
    let (same,diff) = same_diff(v1,v2);
    
    //let mult = v1.err + v2.err;
    let score = same * (0.20) - diff;
    //let score = same - diff;
    //let score = same;
    if score.is_nan(){
        dbg!(same, diff, v1.err, v2.err);
        panic!();
    }
    return score

}

fn overlap_percent(x1: usize, x2: usize, y1: usize, y2: usize) -> f64 {
    let intersect = usize::max(usize::min(x2 - y1 + 1, y2 - x1 + 1), 0);
    let min_length = usize::min(x2 - x1 + 1, y2-y1 + 1);
    let p = intersect as f64 / min_length as f64;
    if p > 1. {
        return 1.;
    //        println!("{},{},{},{}",x1,x2,y1,y2);
    //        panic!();
    } else {
        return p;
    }
}

fn overlap(v1: &Vartig, v2: &Vartig) -> bool{
    if v1.snprange.1 >= v2.snprange.0 &&
        v1.snprange.0 <= v2.snprange.1{
            if overlap_percent(v1.snprange.0, v1.snprange.1, v2.snprange.0, v2.snprange.1) > 0.1{
            return true;
            }
            else{
                return false;
            }
    }
    if v2.snprange.1 >= v1.snprange.0 &&
        v2.snprange.0 <= v1.snprange.1{
            if overlap_percent(v1.snprange.0, v1.snprange.1, v2.snprange.0, v2.snprange.1) > 0.1{
            return true;
            }
            else{
                return false;
            }
    }

    return false;

}

pub fn align_vartig(q_vartig: &[Vartig], r_vartig: &[Vartig]) -> Vec<VartigAln>{
    let mut all_align_hits = vec![];
    let mut pos_to_vartig_r : FxHashMap<usize, Vec<&Vartig>>= HashMap::default();
    
    for vartig in r_vartig.iter(){
        for pos in vartig.allele_vec.keys(){
            let pos_vartigs = pos_to_vartig_r.entry(*pos).or_insert(vec![]);
            pos_vartigs.push(vartig);
        }
    }

    for vartig in q_vartig.iter(){
        let mut set_of_hits : HashSet<usize> = HashSet::default();
        for pos in vartig.allele_vec.keys(){
            if pos_to_vartig_r.contains_key(&pos){
                for vartig_hit in pos_to_vartig_r[&pos].iter(){
                    set_of_hits.insert(vartig_hit.index);
                }
            }
        }
        let mut vartig_hits = vec![];
        for index in set_of_hits{
            if vartig.contig.is_none() || r_vartig[index].contig.is_none(){
                vartig_hits.push(&r_vartig[index]);
            }
            else if vartig.contig == r_vartig[index].contig{
                vartig_hits.push(&r_vartig[index]);
            }
        }
        vartig_hits.sort_by(|x,y| x.snprange.cmp(&y.snprange));
        //Dynamic programming with traceback here
        let alignment_hits = dp_align(vartig, &vartig_hits);
        for align_hit in alignment_hits{
            all_align_hits.push(align_hit);
        }
    }
    all_align_hits
}

fn dp_align<'a>(q_vartig: &Vartig, vartig_hits: &Vec<&'a Vartig>) -> Vec<VartigAln>{

    let mut scores = vec![];
    for i in 0..vartig_hits.len(){
        let hit_i = &vartig_hits[i];
        let score = hit_score(q_vartig,hit_i);
        scores.push(score);
    }
    let max_band = 100;
    let mut traceback :Vec<usize> = (0..scores.len()).collect();
    let mut best_scores = scores.clone();
    for i in 0..vartig_hits.len(){
        let mut curr_best_score = scores[i];
        let v1 = &vartig_hits[i];
        let min_ind = i64::max(0, i as i64 - max_band);
        for j in min_ind as usize ..i{
            let v2 = &vartig_hits[j];
            let ol = overlap(v1,v2);
            //println!("NAME: {}, {} {}-{}, {}-{}; {}",&q_vartig.name, ol, v1.snprange.0, v1.snprange.1, v2.snprange.0, v2.snprange.1, scores[j]);
            if !ol{
                if scores[i] + best_scores[j] > curr_best_score{
                    traceback[i] = j;
                    curr_best_score = scores[i] + best_scores[j]
                }
            }
        }
        best_scores[i] = curr_best_score;
    }


    let mut clone_best_scores = best_scores.clone();
    clone_best_scores.sort_by(|x,y| y.partial_cmp(&x).unwrap());

    let mut num_max_scores = 0;
    for i in 0..clone_best_scores.len(){
        if clone_best_scores[i] == clone_best_scores[0]{
            num_max_scores += 1;
        }
    }
    if num_max_scores > 1{
        num_max_scores = 1;
    }

    let mut best_indices : Vec<(usize,&f64)> = best_scores.iter().enumerate().collect::<Vec<(_,_)>>();
    best_indices.sort_by(|x,y| y.1.partial_cmp(&x.1).unwrap());

    let mut already_hit: HashSet<&String> = HashSet::default();
    let mut toret = vec![];

    for i in 0..num_max_scores{
        let index_of_max = Some(best_indices[i].0);
//        let index_of_max: Option<usize> = best_scores 
//            .iter()
//            .enumerate()
//            .max_by(|(_, a), (_, b)| a.total_cmp(b))
//            .map(|(index, _)| index);

        let mut align_hits = vec![];
        if index_of_max.is_some(){
            let mut index = index_of_max.unwrap();
            while traceback[index] != index{
                align_hits.push(vartig_hits[index]);
                index = traceback[index];
            }
            align_hits.push(vartig_hits[index]);
        }
        for hit in align_hits{
            if !already_hit.contains(&hit.name){
                toret.push(align_stats(q_vartig, hit));
                already_hit.insert(&hit.name);
            }
        }
    }

    return toret;
}

fn align_stats(v1: &Vartig, v2: &Vartig) -> VartigAln{
    let (same,diff) = same_diff(v1,v2);
    let snp_identity = same / (same + diff);
//    let ol_len;
//    if v2.baserange.1 > v1.baserange.0{
//        ol_len = v1.baserange.1 - v2.baserange.0 + 1;
//    }
//    else{
//        ol_len = v2.baserange.1 - v1.baserange.0 + 1;
//    }
//    let gapless_base_identity = 1.0 - diff / ol_len as f64;
    return VartigAln{
        name1: v1.name.clone(),
        name2: v2.name.clone(),
        vtig1_len: v1.allele_vec.len(),
        vtig2_len: v2.allele_vec.len(),
        snp_range1: v1.snprange,
        snp_range2: v2.snprange,
        base_range1: v1.baserange,
        base_range2: v2.baserange,
        snp_identity,
        gapless_base_identity: 0.,
        cov1: v1.cov,
        cov2: v2.cov,
        same,
        diff,
    }
}
