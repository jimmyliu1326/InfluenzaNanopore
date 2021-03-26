# InfluenzaNanopore

## Description
The snakemake pipeline functions to construct the consensus sequence of the **full-length Influenza A genome** from Nanopore reads. The pipeline uses Centrifuge to determine the segment identity (segment 1-8) of each read from a single fastq file and subsequently performs multiple rounds of genome polishing steps per segment using medaka and racon.

## Usage
```
Required arguments:

-i|--input    Path to input samples.csv containing sample name and path to .fastq per line
-o|--output   Path to output directory, the final consensus sequence will be found under consensus/
--db          Path to Centrifuge database for taxonomic and segment classification

Optional arguments:

-t|--threads        Number of threads [Default = 32]
-s|--segment        Target specific Influenza A genomic segments for consensus calling with each segment number delimited by a comma (Example: -s 1,2,5,6)
-m|--model          Specify the flowcell chemistry used for Nanopore sequencing {Options: r9, r10} [Default = r9]
--notrim            Disable adaptor trimming by Porechop
--keep-tmp          Keep all temporary files
-h|--help           Display help message
```

Example command line for pipeline execution:
```
influenza_consensus.sh -i samples.csv -o /path/to/output --db /path/to/centrifuge/database
```

## Dependencies
* R >= 3.6
* medaka == 1.0.3
* racon >= 1.4.13
* centrifuge >= 1.0.3
* seqtk >= 1.3
* snakemake >= 5.30.1
* porechop >= 0.2.4
