use std::collections::HashSet;
use std::fs::File;
use std::io;

use anyhow::Result;
use clap::Parser;

mod app;
mod gtf;

use app::{quantify_bam, GeneMap, QuantMethod, Strandness};

#[derive(Debug, Clone, Parser)]
#[command(author, version, about)]
pub struct Config {
    #[arg(short = 'f', long, value_name = "FILE", required = true)]
    gtf: String,

    #[arg(short, long, value_name = "FILE", required = true)]
    bam: String,

    #[arg(short, long, value_name = "FILE", required = false)]
    output: Option<String>,

    #[arg(short = 'q', long, value_name = "0-255", default_value = "10", value_parser = clap::value_parser!(u8).range(..256))]
    mapq: u8,

    #[arg(value_enum, short, long, default_value = "U")]
    strandness: Strandness,

    #[arg(value_enum, short = 'm', long = "method", default_value = "union")]
    method: QuantMethod,

    #[arg(short = 'd', long)]
    usedups: bool,

    #[arg(short = 's', long)]
    use_supplementary: bool,

    #[arg(short = 'n', long)]
    nosingletons: bool,

    #[arg(short = 'T', long, default_value = "exon", value_parser = parse_seq_types)]
    seq_types: HashSet<String>,
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
                if !seq_types.insert(part.to_string()) {
                    eprintln!("duplicate seq-type {part}");
                }
            }
            _ => {
                return Err(format!(
                    r#"
    {part}? supported are --seq-types "exon", "gene", "mRNA", "cds", "intron", "polyA_sequence"
    "polyA_site", "five_prime_UTR" and "three_prime_UTR".
"#
                ))
            }
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
