use std::fs::File;
use std::io;

use clap::{crate_authors, crate_version, App, Arg};
use anyhow::{Result, Context};

mod gtf;
mod app;

use app::{GeneMap, Config, quantify_bam};

fn main() -> Result<()> {
    let matches = App::new("gensum")
        .author(crate_authors!())
        .version(crate_version!())
        .arg(Arg::with_name("gtf")
                 .help("The .gtf reference transcriptome file")
                 .long("gtf")
                 .short("f")
                 .required(true)
                 .takes_value(true)
                 .value_name("FILE"))
        .arg(Arg::with_name("bam")
                 .help("The bam file to quantify")
                 .long("bam")
                 .short("b")
                 .required(true)
                 .takes_value(true)
                 .value_name("FILE"))
        .arg(Arg::with_name("output")
                 .help("The output file")
                 .long("output")
                 .short("o")
                 .required(false)
                 .takes_value(true)
                 .value_name("FILE"))
        .arg(Arg::with_name("mapq")
                 .help("The minimum required mapping quality to include")
                 .long("mapq")
                 .short("q")
                 .required(false)
                 .takes_value(true)
                 .default_value("10")
                 .value_name("0-255"))
        .arg(Arg::with_name("strandness")
                 .help("The RNA library strandness [F]orward, [R]everse or [U]nstranded")
                 .long("strandness")
                 .short("s")
                 .required(false)
                 .takes_value(true)
                 .possible_values(&["F", "R", "U"])
                 .default_value("U"))
        .arg(Arg::with_name("qmethod")
                 .help("The quantification method, 'strict' or 'union'. 'union' counts all \
                     genes that overlap any part of the reads, 'strict' requires the read to \
                     map within the exon boundaries.")
                 .long("method")
                 .short("m")
                 .required(false)
                 .takes_value(true)
                 .possible_values(&["union", "strict"])
                 .default_value("union"))
        .arg(Arg::with_name("usedups")
                 .help("Also count read (pairs) marked as (optical) duplicate, default \
                     excludes duplicates. Requires a bam files processed with a markdups tool")
                 .long("usedups")
                 .short("d")
                 .required(false)
                 .takes_value(false))
        .arg(Arg::with_name("nosingletons")
                 .help("Do not count paired-end reads that have only 1 mapped end. \
                     Default allows one mapped end. Only affects paired-end reads.")
                 .long("nosingle")
                 .short("n")
                 .required(false)
                 .takes_value(false))
        .get_matches();

    let usedups = matches.is_present("usedups");
    let nosingletons = matches.is_present("nosingletons");
    let mapq: u8 = matches
        .value_of("mapq")
        .unwrap()
        .parse()
        .context("Min mapq cutoff value should be 0-255")?;

    let strandness = matches.value_of("strandness").unwrap().parse().unwrap();
    let method = matches.value_of("qmethod").unwrap().parse().unwrap();

    let config = Config { 
        usedups,
        nosingletons,
        mapq,
        method,
        strandness
    };

    let r = matches.value_of_os("gtf").unwrap();
    let gm = GeneMap::from_gtf(r)?;

    let bam = matches.value_of_os("bam").unwrap();
    let res = quantify_bam(bam, config, &gm)?;

    if let Some(f) = matches.value_of_os("output") {
        let o = File::create(f)?;
        res.write(o, &gm)?;
    } else {
        let stdout = io::stdout();
        stdout.lock();
        res.write(stdout, &gm)?;
    }

    Ok(())
}
