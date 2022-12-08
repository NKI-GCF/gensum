use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::str::FromStr;
use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

mod app;
mod gtf;

use app::{quantify_bam, GeneMap, QuantMethod, Strandness};

#[derive(Debug, Clone, Parser)]
#[command(author, version, about)]
pub struct Config {
    /// The .gtf reference transcriptome file.
    #[arg(short = 'f', long, value_name = "FILE", required = true)]
    gtf: PathBuf,

    /// The bam file to quantify.
    #[arg(short, long, value_name = "FILE", required = true)]
    bam: String,

    /// The output file
    #[arg(short, long, value_name = "FILE", required = false)]
    output: Option<PathBuf>,

    /// The minimum required mapping quality to include.
    #[arg(short = 'q', long, value_name = "0-255", default_value = "10", value_parser = clap::value_parser!(u8).range(..256))]
    mapq: u8,

    /// The RNA library strandness [F]orward, [R]everse or [U]nstranded.
    #[arg(value_enum, short, long, default_value = "U", value_parser = Strandness::from_str)]
    strandness: Strandness,

    /// The quantification method. 'union' counts all genes that overlap any
    /// part of the reads, 'strict' requires the read to map within the exon
    /// boundaries.
    #[arg(value_enum, short = 'm', long = "method", default_value = "union", verbatim_doc_comment)]
    method: QuantMethod,

    /// A comma seperated list of sequence types. Allowed are non overlapping:
    /// exon, gene, mRNA, cds, intron, polyA_sequence, polyA_site, five_prime_UTR
    /// and three_prime_UTR.
    #[arg(short = 'T', long, default_value = "exon", value_parser = parse_seq_types, verbatim_doc_comment)]
    seq_types: HashSet<String>,

    /// Also count read (pairs) marked as (optical) duplicate, default excludes
    /// duplicates. Requires a bam files processed with a markdups tool.
    #[arg(short = 'd', long, verbatim_doc_comment)]
    usedups: bool,

    /// Do not count paired-end reads that have only 1 mapped end.
    /// Default allows one mapped end. Only affects paired-end reads.
    #[arg(short = 'n', long, verbatim_doc_comment)]
    nosingletons: bool,

    /// Include supplementary reads. Secondary reads (alternate alignments) not
    /// included.
    #[arg(short = 's', long, verbatim_doc_comment)]
    use_supplementary: bool,
}

fn parse_seq_types(s: &str) -> Result<HashSet<String>, String> {
    let mut seq_types = HashSet::new();
    for part in s.split(',') {
        match part {
            "gene" => {
                if !seq_types.is_empty() {
                    return Err("\"gene\" overlaps other seq_types".to_string());
                }
            }
            "mRNA" => {
                if !seq_types.is_disjoint(
                    &["exon", "gene", "cds", "five_prime_UTR", "three_prime_UTR"]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                ) {
                    return Err("\"mRNA\" overlaps other requested seq_types".to_string());
                }
            }
            "cds" => {
                if !seq_types.is_disjoint(
                    &["exon", "gene", "mRNA"]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                ) {
                    return Err("\"mRNA\" overlaps other requested seq_types".to_string());
                }
            }
            "exon" | "intron" | "polyA_sequence" | "polyA_site" | "five_prime_UTR"
            | "three_prime_UTR" => {
                if seq_types.contains("gene") {
                    return Err("\"gene\" overlaps other seq_types".to_string());
                } else if seq_types.contains("mRNA") && part != "intron" && !part.contains("polyA")
                {
                    return Err("\"mRNA\" overlaps other seq_types".to_string());
                } else if seq_types.contains("cds") && part == "exon" {
                    return Err("\"cds\" overlaps with exons".to_string());
                }
            }
            _ => {
                return Err(r#"
Use one of "exon", "gene", "mRNA", "cds", "intron", "polyA_sequence", "polyA_site", "five_prime_UTR" and "three_prime_UTR"."#.to_string())
            }
        }
        if !seq_types.insert(part.to_string()) {
            eprintln!("duplicate seq-type {part}");
        }
    }
    Ok(seq_types)
}

fn main() -> Result<()> {
    let config = Config::parse();

    let gm = GeneMap::from_gtf(&config)?;

    let res = quantify_bam(&config, &gm)?;

    if let Some(f) = config.output {
        let o = File::create(f)?;
        res.write(o, &gm)?;
    } else {
        let stdout = io::stdout();
        let stdout = stdout.lock();
        res.write(stdout, &gm)?;
    }

    Ok(())
}
