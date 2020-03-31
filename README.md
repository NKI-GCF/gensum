# GenSum

## RNA-seq gene expression counter
Counts reads on genes. Fast. Inspired by htseq-count `GenSum` support two
quantification methods: union and strict, which are compatible with
`htseq-count` `union` and `intersection-strict`. `GenSum` and `htseq-count` produce
identical counts in the union mode, and minor differences in strict mode (1 in
35000 reads). `GeSum` is much much faster. A 15 million read `bam` file can be
quantified in around 5s while `htseq-count` takes about 11 minutes.


## Input
Required input is a GTF file that contains the exons and the gene-ids. GenSum
uses a very naive parser to extract only 'exon' features. The GTF exon features
are required to have an entry in the attributes that contains: `gene_id: "<the
gene id to count>"`. It it recommended to use the files generated by the ensembl
team at: http://ftp.ensembl.org/pub/current_gtf/

The second input is the .bam file created by an aligner. TopHat/HiSat2/STAR
should all work fine. Stranded libraries as well as paired end data are
supported. When using a stranded RNA library supply the library type using the
`--strandness` flag to restrict counting the correctly oriented reads.
Paired-end reads are expected by in oriented inwards (--->...<---). Sorting,
neither position nor name sorting is required, but paired end reads are store
until the mate is encountered which can affect memory usage.

## Options
```
USAGE:
    gensum [FLAGS] [OPTIONS] --bam <FILE> --gtf <FILE>

FLAGS:
    -h, --help        Prints help information
    -n, --nosingle    Do not count paired-end reads that have only 1 mapped end. Default allows one mapped end. Only
                      affects paired-end reads.
    -d, --usedups     Also count read (pairs) marked as (optical) duplicate, default excludes duplicates. Requires a bam
                      files processed with a markdups tool
    -V, --version     Prints version information

OPTIONS:
    -b, --bam <FILE>                 The bam file to quantify
    -f, --gtf <FILE>                 The .gtf reference transcriptome file
    -q, --mapq <0-255>               The minimum required mapping quality to include [default: 10]
    -o, --output <FILE>              The output file
    -m, --method <qmethod>           The quantification method, 'strict' or 'union'. 'union' counts all genes that
                                     overlap any part of the reads, 'strict' requires the read to map within the exon
                                     boundaries. [default: union]  [possible values: union, strict]
    -s, --strandness <strandness>    The RNA library strandness [F]orward, [R]everse or [U]nstranded [default: U]
                                     [possible values: F, R, U]
```

## Output
The output is a simple two column `<tab>` delimited file. The first column
contains the `gene_id` or a descriptive name for unassigned reads. The second
column the counts on that gene.

