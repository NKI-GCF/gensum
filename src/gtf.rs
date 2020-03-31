use std::io::{self, BufRead};
use std::str::FromStr;

use anyhow::{Result, Context};

pub struct GtfReader<B> {
    r: B,
    last_line: String,
}

impl<B: BufRead> GtfReader<B> {
    pub fn new(r: B) -> GtfReader<B> {
        GtfReader { r, last_line: String::new() }
    }

    /// consume from the buffreader until a record line is read
    pub fn advance_record(&mut self) -> io::Result<bool> {
        loop {
            self.last_line.clear();
            match self.r.read_line(&mut self.last_line)? {
                0 => return Ok(false),
                _ if self.last_line.starts_with('#') => {},
                _ => break
            }
        }
        Ok(true)
    }

    /// attempt to parse the current GTF line as an exon
    /// Returns None for any other type
    /// Fails when unable to parse or required attributes (gene_id)
    /// are not present
    pub fn parse_exon(&self) -> Result<Option<GtfExon>> {
        let line = self.last_line.as_str();
        let mut s = line.split('\t');
        let seq_name = s.next().ok_or_else(|| data_error(line))?;
        //eprintln!("chr {}", seq_name);
        //skip source
        let seq_type = s.nth(1).ok_or_else(|| data_error(line))?;
        //eprintln!("type {}", seq_type);
        if seq_type == "exon" {
            let start = s.next()
                .ok_or_else(|| data_error(line))?
                .parse().context("Invalid start")?;
            let end = s.next()
                .ok_or_else(|| data_error(line))?
                .parse().context("Invalid end")?;
            let strand = s.nth(1)
                .ok_or_else(|| data_error(line))?
                .parse().map_err(|_| data_error(line)).context("Invalid strand")?;

            let attrs = s.nth(1).ok_or_else(|| data_error(line)).context("No attributes")?;
            //eprintln!("attr {}", attrs);

            // split attrs on ';'
            let mut attr = attrs.split(';');
            let id = attr.find(|s| s.starts_with("gene_id "))
                .map(|s| s.split_at(9).1.trim_matches('"'))
                .ok_or_else(|| data_error(line)).context("No gene_id in attributes")?;

            Ok(Some(GtfExon { seq_name, seq_type, start, end, strand, id}))
        } else {
            Ok(None)
        }
    }
}

fn data_error(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}

#[derive(Debug)]
pub struct GtfExon<'a> {
    pub seq_name: &'a str,
    pub seq_type: &'a str,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub id: &'a str
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum Strand {
    Foward,
    Reverse,
    Unknown
}

impl FromStr for Strand {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Foward),
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
        let mut reader = GtfReader::new(Cursor::new(GTF));

        //gene entry
        assert!(reader.advance_record().unwrap());
        assert!(reader.parse_exon().unwrap().is_none());

        //transcript entry
        assert!(reader.advance_record().unwrap());
        assert!(reader.parse_exon().unwrap().is_none());

        // two exons
        assert!(reader.advance_record().unwrap());
        let r = reader.parse_exon().unwrap().unwrap();
        assert_eq!(r.id, "ENSG00000112592");

        assert!(reader.advance_record().unwrap());
        let r = reader.parse_exon().unwrap().unwrap();
        assert_eq!(r.id, "ENSG00000112592");

        // and a CDS
        assert!(reader.advance_record().unwrap());
        assert!(reader.parse_exon().unwrap().is_none());

        //EOF
        assert!(!reader.advance_record().unwrap());

    }
}


