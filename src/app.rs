use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Write, BufWriter};
use std::ops::Range;
use std::path::Path;
use std::cmp::{Ord, PartialOrd, Ordering};
use std::str::FromStr;


use anyhow::{anyhow, Result};
use indexmap::IndexSet;
use itoa;
use nclist::{NClist, Interval};
use rust_htslib::{bam, bam::Read, bam::record::Cigar};

use crate::gtf::{GtfReader, Strand};


#[derive(Debug, Copy, Clone)]
pub struct Config {
    pub usedups: bool,
    pub nosingletons: bool,
    pub mapq: u8,
    pub method: QuantMethod,
    pub strandness: Strandness
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum QuantMethod {
    Union,
    Strict
}

impl FromStr for QuantMethod {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "strict" => Ok(QuantMethod::Strict),
            "union" => Ok(QuantMethod::Union),
            _ => Err("Unknown quantification method")
        }
    }
}


#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Strandness {
    Forward,
    Reverse,
    Unstranded
}

impl Strandness {
    #[inline]
    fn matches_bam_record(self, r: &bam::Record, target: Strand) -> bool {
        if self == Strandness::Unstranded {
            return true;
        }
        
        // Asuming a FR library ( --->____<--- )
        // Determine if fragment is forward
        let fragment_forward = if r.is_paired() {
            (r.is_first_in_template() && !r.is_reverse()) ||
                (r.is_last_in_template() && r.is_reverse())
        } else {
            !r.is_reverse()
        };

        match (self, target) {
            (_, Strand::Unknown) => true,
            (Strandness::Forward, Strand::Forward) |
            (Strandness::Reverse, Strand::Reverse) => fragment_forward,
            (Strandness::Forward, Strand::Reverse) |
            (Strandness::Reverse, Strand::Forward) => !fragment_forward,
            (Strandness::Unstranded, _) => unreachable!(),
        }
    }
}

impl FromStr for Strandness {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "F" => Ok(Strandness::Forward),
            "R" => Ok(Strandness::Reverse),
            "U" => Ok(Strandness::Unstranded),
            _ => Err("Unknown strandness")
        }
    }
}

/// Exon is defined by it's coordinates and references a parent Gene
#[derive(Debug, Eq, PartialEq)]
struct Exon {
    id: usize,
    strand: Strand,
    range: Range<i64>,
}

impl Ord  for Exon {
    fn cmp(&self, other: &Self) -> Ordering {
        self.id.cmp(&other.id)
            .then(self.range.start.cmp(&other.range.start))
            .then(self.range.end.cmp(&other.range.end))
    }
}

impl PartialOrd  for Exon {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.id.cmp(&other.id)
            .then(self.range.start.cmp(&other.range.start))
            .then(self.range.end.cmp(&other.range.end)))
    }
}

impl Interval for Exon {
    type Coord = i64;
    fn start(&self) -> &Self::Coord {
        &self.range.start
    }

    fn end(&self) -> &Self::Coord {
        &self.range.end
    }
}

fn get_index_or_insert_owned(map: &mut IndexSet<String>, v: &str) -> usize {
    if !map.contains(v) {
        map.insert_full(v.to_owned()).0
    } else {
        map.get_full(v).unwrap().0
    }
}

pub struct GeneMap {
    genes: IndexSet<String>,
    seq_names: IndexSet<String>,
    intervals: Vec<NClist<Exon>>,
}

impl GeneMap {
    pub fn from_gtf<P: AsRef<Path>>(p: P) -> Result<GeneMap> {
        //open gtf
        let f = File::open(p)?;
        let mut reader = GtfReader::new(BufReader::new(f));
        
        let mut genes = IndexSet::new();
        let mut seq_names = IndexSet::new();
        let mut exons = Vec::new();


        //iterate records
        let mut gtfline = String::new();
        let mut n = 0;
        loop {
            n += 1;
            gtfline.clear();
            if !reader.advance_record()? {
                break;
            }

            if let Some(r) = reader.parse_exon()? {

                let gene_idx = get_index_or_insert_owned(&mut genes, r.id);
                let chr_idx = get_index_or_insert_owned(&mut seq_names, r.seq_name);

                 if r.end - r.start < 0 {
                     eprintln!("Yikes: {} {} {}", r.start, r.end, gtfline);
                     continue;
                 }

                if exons.len() == chr_idx {
                    exons.push(Vec::new());
                }

                // gtf exon coordinates are 1 based and closed end
                // bam files are 0 based, and nclist expects half open
                exons[chr_idx].push(Exon {id: gene_idx, strand: r.strand, range: r.start-1..r.end });
            }
        }

        //Create the NClists
        let mut numexons = 0;
        let mut numexonsdd = 0;
        let intervals = exons.into_iter().map(|mut v| {
            numexons += v.len();
            v.sort();
            v.dedup(); //deduplication saves around 50% because of comparable isoforms (and havanna entries)
            numexonsdd += v.len();
            NClist::from_vec(v)
                .map_err(|_| anyhow!("Cannot create interval search list, all ranges must be > 1"))
        }).collect::<Result<_, _>>()?;

        eprintln!("{} lines in GTF, parsed {} exons, {} unique geneid-exon ranges", n, numexons, numexonsdd);

        Ok(GeneMap { genes, seq_names, intervals })
    }

    #[inline]
    pub fn hit_name(&self, i: usize) -> Option<&String> {
        self.genes.get_index(i)
    }

}

#[derive(Eq, PartialEq)]
enum SegmentHit {
    Hit(usize),
    Nohit,
    Ambiguous
}

#[derive(Default)]
pub struct ReadMappings {
    qc_failed: usize,
    unmapped: usize,
    secondary: usize,
    duplicated: usize,
    ambiguous: usize,
    ambiguous_pair: usize,
    notingtf: usize,
    mapq: usize,
    nohit: usize,
    hit: Vec<usize>
}

impl ReadMappings {
    pub fn new(n: usize) -> ReadMappings {
        ReadMappings { hit: vec![0; n], ..Default::default() }
    }

    fn count_hit(&mut self, h: SegmentHit) {
        match h {
            SegmentHit::Nohit => self.nohit += 1,
            SegmentHit::Ambiguous => self.ambiguous += 1,
            SegmentHit::Hit(id) => self.hit[id] += 1,
        }
    }

    pub fn write<W: Write>(&self, o: W, genes: &GeneMap) -> Result<()> {

        let mut w = BufWriter::new(o);
        for (geneidx, &count) in self.hit.iter().enumerate() {
            w.write_all(genes.hit_name(geneidx).unwrap().as_bytes())?;
            w.write_all(&[b'\t'])?;
            itoa::write(&mut w, count)?;
            w.write_all(&[b'\n'])?;
        }

        writeln!(w, "qc_failed\t{}", self.qc_failed)?;
        writeln!(w, "unmapped\t{}", self.unmapped)?;
        writeln!(w, "secondary_alignments\t{}", self.secondary)?;
        writeln!(w, "marked_duplicated\t{}", self.duplicated)?;
        writeln!(w, "ambiguous\t{}", self.ambiguous)?;
        writeln!(w, "ambiguous_pair\t{}", self.ambiguous_pair)?;
        writeln!(w, "chr_not_in_gtf\t{}", self.notingtf)?;
        writeln!(w, "nohit\t{}", self.nohit)?;

        Ok(())
    }
}

pub fn quantify_bam<P: AsRef<Path>>(bam_file: P, config: Config, genemap: &GeneMap) -> Result<ReadMappings> {
    //open bam
    let mut bam = bam::Reader::from_path(bam_file)?;
    // test from command line show improve until 4 cpu's
    bam.set_threads(4)?;

    //intersect header chr list with rr
    let header = bam.header();
    let tid_map: Vec<_> = header.target_names().iter()
        .map(|v| String::from_utf8_lossy(v))
        .map(|name| genemap.seq_names.iter().position(|n| name == n.as_ref())).collect();

    //quantify
    let mut delayed = HashMap::new();
    let mut counts = ReadMappings::new(genemap.genes.len());


    for record in bam.records() {
        let record = record?;
            if record.is_unmapped() {
                counts.unmapped += 1;
                continue;
            }

            if record.is_quality_check_failed() {
                counts.qc_failed += 1;
            }
            if record.is_secondary() || record.is_supplementary() {
                counts.secondary += 1;
                continue;
            }

            if !config.usedups && record.is_duplicate() {
                counts.duplicated += 1;
                continue;
            }

            if record.mapq() < config.mapq {
                counts.mapq += 1;
                continue;
            }

            if let Some(ref_chr_id) = tid_map[record.tid() as usize] {
                let ref_chr_map = &genemap.intervals[ref_chr_id];
                if record.is_paired() {
                    if record.is_mate_unmapped() && !config.nosingletons {
                        counts.count_hit(map_segments(&record, ref_chr_map, config));
                    } else {
                        //is the mate on the same chromosome? if not than this read pair is ambiguous
                        if record.tid() != record.mtid() {
                            counts.ambiguous_pair += 1;
                        } else if let Some(mate) = delayed.remove(record.qname()) {
                            let m1 = map_segments(&record, ref_chr_map, config);
                            let m2 = map_segments(&mate, ref_chr_map, config);
                            if m1 == m2 {
                                counts.count_hit(map_segments(&record, ref_chr_map, config));
                            } else {
                                counts.ambiguous_pair += 1;
                            }
                        } else {
                            delayed.insert(record.qname().to_vec(), record);
                        }
                    }
                } else {
                    //Single-end read
                    counts.count_hit(map_segments(&record, ref_chr_map, config));
                }
            } else {
                // this chr was not in the gtf
                counts.notingtf += 1;
            }
    }
    Ok(counts)
}



fn map_segments(r: &bam::Record, map: &NClist<Exon>, config: Config) -> SegmentHit {
    //Store the first gene hit id
    let mut  target_id = None;
    let cigar = r.cigar();

    let strict = config.method == QuantMethod::Strict;
    let strandness = config.strandness;

    // use cigar line to filter the cigar to ranges (o) that lie on the genome
    for o in cigar.into_iter().scan(r.pos(), |pos, c| {
        match c {
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                *pos += *n as i64;
                Some(None)
            },
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                Some(None)
            },
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let r = *pos..*pos + *n as i64;
                *pos += *n as i64;
                Some(Some(r))
            }
        }
    }).filter_map(|c| c)
    {
        //match this segment's genomic region to exons and filter based on program configuration
        let exons =  map.overlaps(&o)
            .filter(|e| !strict || (o.start >= *e.start() && o.end <= *e.end()))
            .filter(|e| strandness.matches_bam_record(r, e.strand));

        // check that all overlapping exons map to the same gene
        let mut segment_ambiguous = false;
        let mut segment_id = None;

        for exon in exons {
            if let Some(id) = segment_id {
                if !strict && (id != exon.id) {
                    // in  union mode any part linking to a different gene makes it ambiguous
                    return SegmentHit::Ambiguous;
                } else if strict && id != exon.id {
                    // in strict mode ambigous segments can be recued if a unique mapping is 
                    // available from other segments
                    segment_ambiguous = true;
                    break;
                } 
            } else {
                segment_id = Some(exon.id)
            }
        }

        //strict requires al segments overlap the same gene
        if strict && segment_id.is_none() {
            return SegmentHit::Nohit;
        }

        if !segment_ambiguous {
            if target_id.is_some() && segment_id.is_some() && target_id != segment_id {
                return SegmentHit::Ambiguous
            }
            if target_id.is_none() && segment_id.is_some() {
                target_id = segment_id;
            }
        } 
    }

    if let Some(id) = target_id {
        SegmentHit::Hit(id)
    } else {
        SegmentHit::Nohit
    }
}



