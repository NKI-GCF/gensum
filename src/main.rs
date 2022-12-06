use std::fs::File;
use std::io;

use clap::Parser;
use anyhow::Result;

mod gtf;
mod app;

use app::{GeneMap, quantify_bam, Strandness, QuantMethod};

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
}

fn main() -> Result<()> {
    let config = Config::parse();

    let gm = GeneMap::from_gtf(config.gtf.as_str())?;

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
