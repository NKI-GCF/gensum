use std::io::{self, BufRead};
use std::str::FromStr;

use anyhow::{Result, Context, anyhow};

pub struct GtfReader<B> {
    s: io::Split<B>,
    last_line: Vec<u8>,
}

impl<B: BufRead> GtfReader<B> {
    pub fn new(r: B) -> Result<GtfReader<B>> {
        let mut s = r.split(b'\n');
        while let Some(last_line) = s.next().transpose()? {
            if last_line[0] != b'#' {
                return Ok(GtfReader { s, last_line })
            }
        }
        Err(anyhow!("No or only comment lines in gtf"))
    }

    /// consume from the buffreader until a record line is read
    pub fn advance_record(&mut self) -> Result<bool> {
        if let Some(line) = self.s.next().transpose()? {
            self.last_line = line;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// attempt to parse the current GTF line as an exon
    /// Returns None for any other type
    /// Fails when unable to parse or required attributes (gene_id)
    /// are not present
    pub fn parse_exon(&self) -> Result<Option<GtfExon>> {
        let line = &self.last_line;
        let mut s = line.split(|&c| c == b'\t');
        let seq_name = s.next().ok_or_else(|| data_error(&line))?;
        //eprintln!("chr {}", seq_name);
        //skip source
        let seq_type = s.nth(1).ok_or_else(|| data_error(&line))?;
        //eprintln!("type {}", seq_type);
        if seq_type == b"exon" {
            let start = s.next()
                .map(|u| std::str::from_utf8(&u)).transpose()?
                .ok_or_else(|| data_error(&line))?
                .parse().context("Invalid start")?;
            let end = s.next()
                .map(|u| std::str::from_utf8(&u)).transpose()?
                .ok_or_else(|| data_error(&line))?
                .parse().context("Invalid end")?;
            let strand = s.nth(1)
                .map(|u| std::str::from_utf8(&u)).transpose()?
                .ok_or_else(|| data_error(&line))?
                .parse().map_err(|_| data_error(&line)).context("Invalid strand")?;

            let attrs = s.nth(1).ok_or_else(|| data_error(&line)).context("No attributes")?;
            //eprintln!("attr {}", attrs);

            // split attrs on ';'
            let mut attr = attrs.split(|&c| c == b';');
            let id = attr.find(|s| &s[0..8] == b"gene_id ")
                .map(|s| &s[9..(s.len()-1)])
                .ok_or_else(|| data_error(&line)).context("No gene_id in attributes")?;

            Ok(Some(GtfExon { seq_name, seq_type, start, end, strand, id}))
        } else {
            Ok(None)
        }
    }
}

fn data_error(s: &[u8]) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, std::str::from_utf8(&s).unwrap())
}

#[derive(Debug)]
pub struct GtfExon<'a> {
    pub seq_name: &'a [u8],
    pub seq_type: &'a [u8],
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

impl FromStr for Strand {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err("GTF strand not +/-/.")
        }
    }
}


#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::GtfReader;

    const GTF:&str = r#"#!genome-build GRCh38.p12
6	ensembl_havana	gene	170554302	170572870	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
6	havana	transcript	170554302	170566957	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	exon	170554302	170554463	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "1"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001701648"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	exon	170556882	170557083	.	+	.	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "2"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; exon_id "ENSE00001510679"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
6	havana	CDS	170557030	170557083	.	+	0	gene_id "ENSG00000112592"; gene_version "13"; transcript_id "ENST00000421512"; transcript_version "5"; exon_number "2"; gene_name "TBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TBP-203"; transcript_source "havana"; transcript_biotype "protein_coding"; protein_id "ENSP00000400008"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF"; transcript_support_level "1";
"#;


    #[test]
    fn read() {
        let mut reader = GtfReader::new(Cursor::new(GTF)).unwrap();

        //gene entry
        assert!(reader.parse_exon().unwrap().is_none());
        assert!(reader.advance_record().unwrap());

        //transcript entry
        assert!(reader.parse_exon().unwrap().is_none());
        assert!(reader.advance_record().unwrap());

        // two exons
        let r = reader.parse_exon().unwrap().unwrap();
        assert_eq!(r.id, b"ENSG00000112592");
        assert!(reader.advance_record().unwrap());

        let r = reader.parse_exon().unwrap().unwrap();
        assert_eq!(r.id, b"ENSG00000112592");
        assert!(reader.advance_record().unwrap());

        // and a CDS
        assert!(reader.parse_exon().unwrap().is_none());
        //EOF
        assert!(!reader.advance_record().unwrap());

    }
}


