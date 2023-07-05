use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use anyhow::{Result, Context};

mod gtf;
mod app;

use app::{GeneMap, QuantMethod, Strandness, quantify_bam};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None, max_term_width = 120)]
pub struct Args {
    /// The bam file to quantify
    #[clap(short, long, value_name = "FILE")]
    bam: PathBuf,

    /// The .gtf reference transcriptome file. This file may be (b)gzipped.
    #[clap(short, long, value_name = "FILE")]
    gtf: PathBuf,

    /// The output file (TXT), default: stdout
    #[clap(short, long, value_name = "FILE")]
    out: Option<PathBuf>,

    /// The quantification method, 'strict' or 'union'. 'union' counts all genes that overlap any
    /// part of the reads, 'strict' requires the read to map within the exon boundaries
    #[clap(long, short, default_value = "union")]
    method: QuantMethod,

    /// The RNA library strandness [F]orward, [R]everse or [U]nstranded
    #[clap(long, short, default_value = "U")]
    strandness: Strandness,

    /// The minimum required mapping quality required for a read to be counted
    #[clap(long, short = 'q', value_name = "0-255", default_value_t = 10)]
    mapq: u8,

    /// Also count read (pairs) marked as (optical) duplicate, default excludes duplicates.
    /// Requires a bam files processed with a markduplicates tool
    #[clap(long, short = 'd')]
    usedups: bool,

    /// Do not count paired-end reads that have only 1 mapped end (singletons). Default allows one
    /// mapped end.  Only affects paired-end reads.
    #[clap(long = "nosingle")]
    nosingletons: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let gm = GeneMap::from_gtf(&args.gtf)?;

    let res = quantify_bam(&args.bam, &args, &gm)?;

    if let Some(f) = args.out.as_ref() {
        let o = File::create(f)?;
        res.write(o, &gm)?;
    } else {
        let stdout = io::stdout();
        let stdout = stdout.lock();
        res.write(stdout, &gm)?;
    }

    Ok(())
}
