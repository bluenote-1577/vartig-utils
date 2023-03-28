use clap::{AppSettings, Arg, Command, SubCommand};
use std::collections::HashSet;
use vtig::align;
use vtig::vartig;

fn main() {
    let matches = Command::new("vartig")
        .version("0.1.0")
        .setting(AppSettings::ArgRequiredElseHelp)
        .about("TODO")
        .subcommand(
            SubCommand::with_name("map")
                .arg(
                    Arg::new("haplotigs1")
                        .value_name("HAPLOTIGS FILE")
                        .index(1)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("haplotigs2")
                        .value_name("HAPLOTIGS FILE")
                        .index(2)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("match cutoff")
                        .value_name("INT")
                        .short('m')
                        .takes_value(true),
                ),
        )
        .subcommand(
            SubCommand::with_name("dist")
                .arg(
                    Arg::new("haplotigs1")
                        .value_name("HAPLOTIGS FILE")
                        .index(1)
                        .takes_value(true),
                )
                .arg(
                    Arg::new("match cutoff")
                        .value_name("INT")
                        .short('m')
                        .takes_value(true),
                )
                .arg(
                    Arg::new("haplotigs2")
                        .value_name("HAPLOTIGS FILE")
                        .index(2)
                        .takes_value(true),
                ),
        )
        .get_matches();

    //let file_name = "haplotigs.fa";

    if matches.subcommand_name().unwrap() == "map" {
        let matches = matches.subcommand_matches("map").unwrap();
        let file_name1 = matches.value_of("haplotigs1").unwrap();
        let file_names2 = matches.values_of("haplotigs2").unwrap();
        let match_cutoff = matches
            .value_of("match cutoff")
            .unwrap_or("0")
            .parse::<f64>()
            .unwrap();

        let tig1 = vartig::get_vartigs_from_file(&file_name1);
        let mut tigs2 = vec![];
        for file_name in file_names2 {
            let tig2 = vartig::get_vartigs_from_file(&file_name);
            tigs2.extend(tig2);
        }

        let forward = align::align_vartig(&tig1, &tigs2);
        let backward = align::align_vartig(&tigs2, &tig1);

        let mut backward_match_set: HashSet<(&str, &str)> = HashSet::default();
        let mut good_matches = vec![];
        for aln in backward.iter() {
            backward_match_set.insert((&aln.name2, &aln.name1));
        }
        for aln in forward {
            if backward_match_set.contains(&(&aln.name1, &aln.name2)) {
                good_matches.push(aln);
            }
        }

        let mut used_matches: HashSet<String> = HashSet::default();

        for aln in good_matches {
            if !used_matches.contains(&aln.name1) {
                println!("#{}\tLEN:{}", aln.name1, aln.vtig1_len);
                used_matches.insert(aln.name1.clone());
            }
            if aln.same + aln.diff > match_cutoff {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    aln.name1,
                    aln.name2,
                    aln.snp_identity,
                    aln.same,
                    aln.diff,
                    aln.cov1,
                    aln.cov2,
                    aln.gapless_base_identity,
                    aln.snp_range1.0,
                    aln.snp_range1.1,
                    aln.snp_range2.0,
                    aln.snp_range2.1,
                    aln.base_range1.0,
                    aln.base_range1.1,
                    aln.base_range2.0,
                    aln.base_range2.1
                );
            }
        }
    }
    if matches.subcommand_name().unwrap() == "dist" {
        let len_cutoff = 2;
        let error_cutoff = 0.0;
        let match_cutoff = matches
            .value_of("match cutoff")
            .unwrap_or("0")
            .parse::<f64>()
            .unwrap();
        let matches = matches.subcommand_matches("dist").unwrap();
        let file_name1 = matches.value_of("haplotigs1").unwrap();
        let file_name2 = matches.value_of("haplotigs2").unwrap();
        let tig1 = vartig::get_vartigs_from_file(&file_name1);
        let tig2 = vartig::get_vartigs_from_file(&file_name2);
        let avg_cov_1 = tig1.iter().map(|x| x.cov).sum::<f64>() / (tig1.len() as f64);
        let avg_cov_2 = tig2.iter().map(|x| x.cov).sum::<f64>() / (tig2.len() as f64);

        let forward = align::align_vartig(&tig1, &tig2);
        let backward = align::align_vartig(&tig2, &tig1);

        let mut backward_match_set: HashSet<(&str, &str)> = HashSet::default();
        let mut good_matches = vec![];
        for aln in backward.iter() {
            backward_match_set.insert((&aln.name2, &aln.name1));
        }
        for aln in forward {
            if backward_match_set.contains(&(&aln.name1, &aln.name2)) {
                good_matches.push(aln);
            }
        }

        let tig1_pass = tig1
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .map(|x| x.name.clone())
            .collect::<HashSet<String>>();
        let tig2_pass = tig2
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .map(|x| x.name.clone())
            .collect::<HashSet<String>>();
        let total_alleles_1: usize = tig1
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .map(|x| x.allele_vec.len())
            .sum();
        let total_alleles_2: usize = tig2
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .map(|x| x.allele_vec.len())
            .sum();

        let mut cov_disc = 0.;
        let mut aligned_ctgs_1: HashSet<String> = HashSet::default();
        let mut aligned_ctgs_2: HashSet<String> = HashSet::default();
        let mut total_aligned = 0.;
        let mut errors = 0.;
        for aln in good_matches.iter() {
            if aln.same + aln.diff > match_cutoff {
                cov_disc += f64::abs((aln.cov1 / avg_cov_1 * avg_cov_2 / aln.cov2).log(2.));
                if tig1_pass.contains(&aln.name1) || tig2_pass.contains(&aln.name2) {
                    total_aligned += aln.same + aln.diff;
                    errors += aln.diff
                }
                if tig1_pass.contains(&aln.name1) {
                    aligned_ctgs_1.insert(aln.name1.clone());
                }
                if tig2_pass.contains(&aln.name2) {
                    aligned_ctgs_2.insert(aln.name2.clone());
                }
            }
        }

        let cutoff_tig1_len = tig1
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .collect::<Vec<&vartig::Vartig>>()
            .len();
        let cutoff_tig2_len = tig2
            .iter()
            .filter(|x| x.allele_vec.len() > len_cutoff)
            .collect::<Vec<&vartig::Vartig>>()
            .len();

        let cov_disc = cov_disc / (good_matches.len() as f64);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "err_rate",
            "cov_discrepancy",
            "aligned_alleles_1_frac",
            "aligned_alleles_2_frac",
            "total_allele_align",
            "total_allele_1",
            "total_allele_2"
        );
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            errors / total_aligned,
            cov_disc,
            total_aligned / total_alleles_1 as f64,
            total_aligned / total_alleles_2 as f64,
            total_aligned,
            total_alleles_1 as f64,
            total_alleles_2 as f64
        );
    }
}
