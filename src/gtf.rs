use std::convert::TryFrom;
use std::io::{self, Read, BufRead, BufReader};

use anyhow::{Result, Context};
use atoi::atoi;

pub struct GtfReader<R> {
    reader: BufReader<R>,
}

impl<R: Read> GtfReader<R> {
    pub fn new(r: R) -> GtfReader<R> {
        let reader = BufReader::new(r);
        GtfReader { reader }
    }

    pub fn read_record(&mut self, record: &mut GtfRecord) -> io::Result<usize> {
        loop {
            let n = self.reader.read_until(b'\n', record.clear_buf_mut())?;
            if !record.is_comment() {
                break Ok(n);
            }
        }
    }
}

pub struct GtfRecord(Vec<u8>);

impl GtfRecord {
    pub fn new() -> GtfRecord {
        GtfRecord(Vec::new())
    }

    fn clear_buf_mut(&mut self) -> &mut Vec<u8> {
        self.0.clear();
        &mut self.0
    }

    pub fn is_comment(&self) -> bool {
        self.0.first() == Some(&b'#')
    }

    /// attempt to parse the current GTF record as an exon
    /// Returns None for any other type
    /// Fails when unable to parse or required attributes (gene_id)
    /// are not present
    pub fn parse_exon(&self) -> Result<Option<GtfExon>> {
        let mut s = self.0.split(|&b| b == b'\t');
        let seq_name = s.next()
            .ok_or_else(|| data_error(&self.0))
            .context("No seqname in gtf line")?;
        //skip source
        let seq_type = s.nth(1)
            .ok_or_else(|| data_error(&self.0))
            .context("No seqtype in gtf line")?;
        //eprintln!("type {}", seq_type);
        if seq_type == b"exon" {
            let start = s.next().and_then(atoi)
                .ok_or_else(|| data_error(&self.0))
                .context("Invalid start")?;
            let end = s.next().and_then(atoi)
                .ok_or_else(|| data_error(&self.0))
                .context("Invalid end")?;
            let strand = s.nth(1)
                .ok_or_else(|| "No strand")
                .and_then(Strand::try_from)
                .map_err(|_| data_error(&self.0))
                .context("Invalid strand")?;

            let attrs = s.nth(1).ok_or_else(|| data_error(&self.0)).context("No attributes")?;

            // split attrs on ';'
            // in the ensembl gtf the gene_id is the first entry so this is not
            // really necessary.
            let mut attr = attrs.split(|&b| b == b';');
            let id = attr.find(|s| s.starts_with(b"gene_id "))
                .map(|s| &s[9..s.len()-1])
                .ok_or_else(|| data_error(&self.0)).context("No gene_id in attributes")?;

            Ok(Some(GtfExon { seq_name, start, end, strand, id}))
        } else {
            Ok(None)
        }
    }
}

impl std::fmt::Display for GtfRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8_lossy(&self.0))
    }
}

fn data_error(s: &[u8]) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, String::from_utf8_lossy(s))
}

#[derive(Debug)]
pub struct GtfExon<'a> {
    pub seq_name: &'a [u8],
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub id: &'a [u8]
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown
}

impl TryFrom<&[u8]> for Strand {
    type Error = &'static str;
    fn try_from(s: &[u8]) -> Result<Self, Self::Error> {
        match s {
            b"+" => Ok(Strand::Forward),
            b"-" => Ok(Strand::Reverse),
            b"." => Ok(Strand::Unknown),
            _ => Err("GTF strand not +/-/.")
        }
    }
}


#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::*;

    const GTF:&str = r#"#!genome-build GRCh38.p12
6	ensembl_havana	gene	170554302	170572870	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
6	havana	transcript	170554302	170566957	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	exon	170554302	170554463	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "1"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001701648"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	exon	170556882	170557083	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "2"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001510679"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	CDS	170557030	170557083	.	+	0	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "2"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000400008"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
"#;


    #[test]
    fn read() {
        let mut reader = GtfReader::new(Cursor::new(GTF));
        let mut record = GtfRecord::new();

        //gene entry
        assert!(matches!(reader.read_record(&mut record), Ok(n) if n > 0));
        assert!(matches!(record.parse_exon(), Ok(None)));

        //transcript entry
        assert!(matches!(reader.read_record(&mut record), Ok(n) if n > 0));
        assert!(matches!(record.parse_exon(), Ok(None)));

        // two exons
        assert!(matches!(reader.read_record(&mut record), Ok(n) if n > 0));
        assert!(matches!(record.parse_exon(), Ok(Some(r)) if r.id == b"ENSG00000112592"));

        assert!(matches!(reader.read_record(&mut record), Ok(n) if n > 0));
        assert!(matches!(record.parse_exon(), Ok(Some(r)) if r.id == b"ENSG00000112592"));

        // and a CDS
        assert!(matches!(reader.read_record(&mut record), Ok(n) if n > 0));
        assert!(matches!(record.parse_exon(), Ok(None)));

        //EOF
        assert!(matches!(reader.read_record(&mut record), Ok(0)));
    }
}


