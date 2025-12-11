# AnnoOddities

AnnoOddies is a Python utility for detecting, identifying and characterising oddities in genome annotations. It parses and integrates statistics from multiple tools - including [AGAT](https://github.com/NBISweden/AGAT), [Mikado](https://github.com/EI-CoreBioinformatics/mikado), and [GFFread](https://github.com/gpertea/gffread), to generate a harmonised set of extended metrics that help assess annotation quality and highlight potential anomalies within genome annotations.

In addition to summarising standard outputs, AnnoOddities computes a range of additional “oddity” measures, for example, counts of unusually large or small introns, exons, and untranslated regions; detection of transcripts lacking start or stop codons; identification of in-frame stop codons; estimation of canonical intron proportions; and more. These metrics help flag potential structural inconsistencies or systematic issues in gene models. The consolidated outputs include both summary statistics and detailed per-transcript results, exported in widely used machine- and human-readable formats (TSV, CSV, JSON, YAML, TOML) for flexible downstream analysis and visualisation.

Currently, AnnoOddities runs [AGAT](https://github.com/NBISweden/AGAT), [Mikado](https://github.com/EI-CoreBioinformatics/mikado), and [GFFread](https://github.com/gpertea/gffread) to generate and combine statistics into unified reports. Future developments will focus on extending the workflow to incorporate additional quality assessment tools such as BUSCO and OMArk, enabling integrated comparisons of annotation and genome completeness metrics.

This tool was developed as part of **[BioHackathon 2025](https://biohackathon-europe.org/), Project 23: Streamlining FAIR Metadata for Biodiversity Genome Annotations**.
Project details: https://github.com/elixir-europe/biohackathon-projects-2025/blob/main/23.md


## Installation

All installation methods below will install AnnoOddities along with its dependencies.

### Docker Installation
AnnoOddities can be installed with Docker. If you don't have Docker, please install [docker](https://docs.docker.com/get-docker/) first. Then you can pull the Docker image with AnnoOddities installed

```console
VERSION=0.1.0
docker run gemygk/annooddities:v${VERSION} annooddities -h
```

### Singularity Installation
AnnoOddities can be installed with Singularity. If you don't have Singularity, please install [singularity](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html#quick-installation-steps) first. Then you can pull the singularity image with AnnoOddities installed.

We can directly run AnnoOddities from the Singularity image hosted on DockerHub
```console
VERSION=0.1.0
singularity exec docker://gemygk/annooddities:v${VERSION} annooddities -h
```

Or, we can build and run a Singularity image, following the steps below:
```console
# Create a Singularity definition file, like below:

$ cat annooddities-0.1.0.def
bootstrap: docker
from: gemygk/annooddities:v0.1.0

# Build the Singularity image
$ sudo singularity build annooddities-0.1.0.sif annooddities-0.1.0.def

# Execute AnnoOddities from the Singularity image
$ singularity exec annooddities-0.1.0.sif annooddities -h
```

### Usage
```console
$ annooddities -h
usage: annooddities.py [-h] --genome_fasta GENOME_FASTA --gff3_file GFF3_FILE [--five_prime_utr_length FIVE_PRIME_UTR_LENGTH] [--three_prime_utr_length THREE_PRIME_UTR_LENGTH] [--five_utr_num FIVE_UTR_NUM] [--three_utr_num THREE_UTR_NUM]
                        [--min_intron_length MIN_INTRON_LENGTH] [--max_intron_length MAX_INTRON_LENGTH] [--min_exon_length MIN_EXON_LENGTH] [--max_exon_length MAX_EXON_LENGTH] [--selected_cds_fraction SELECTED_CDS_FRACTION]
                        [--canonical_intron_motifs CANONICAL_INTRON_MOTIFS] [--output_prefix OUTPUT_PREFIX] [--force] [--verbosity {debug,info,warning,error,critical}]

Find annotation oddities from GFF3 and genome FASTA files

optional arguments:
  -h, --help            show this help message and exit
  --genome_fasta GENOME_FASTA
                        Provide Genome FASTA file
  --gff3_file GFF3_FILE
                        Provide GFF3 file with transcript annotations
  --output_prefix OUTPUT_PREFIX
                        Provide sample prefix for the output table [default:output]
  --force               Force rerun even if output files exist
  --verbosity {debug,info,warning,error,critical}
                        Set logging verbosity level [default:info]

Oddity Thresholds:
  --five_prime_utr_length FIVE_PRIME_UTR_LENGTH
                        Threshold for 5' UTR length oddity [default:10000]
  --three_prime_utr_length THREE_PRIME_UTR_LENGTH
                        Threshold for 3' UTR length oddity [default:10000]
  --five_utr_num FIVE_UTR_NUM
                        Threshold for 5' UTR exon number oddity [default:5]
  --three_utr_num THREE_UTR_NUM
                        Threshold for 3' UTR exon number oddity [default:4]
  --min_intron_length MIN_INTRON_LENGTH
                        Threshold for minimum intron length oddity [default:5]
  --max_intron_length MAX_INTRON_LENGTH
                        Threshold for maximum intron length oddity.
                        Below are some guidelines:
                        - For fungi species, consider setting this to 1000 bp.
                        - For plant species, consider setting this to 10000 bp.
                        - For invertebrates species, consider setting this to 60000 bp.
                        - For vertebrates species, consider setting this to 120000 bp.
                        [default:120000]
  --min_exon_length MIN_EXON_LENGTH
                        Threshold for minimum exon length oddity [default:5]
  --max_exon_length MAX_EXON_LENGTH
                        Threshold for maximum exon length oddity [default:10000]
  --selected_cds_fraction SELECTED_CDS_FRACTION
                        Threshold for selected CDS fraction oddity.
                        This is the proportion of coding sequence to that of the transcript.
                        Values range from 0.0 to 1.0 [default:0.3]
  --canonical_intron_motifs CANONICAL_INTRON_MOTIFS
                        Comma-separated list of canonical intron motifs. [default:'GT-AG,GC-AG,AT-AC']
```

## Running AnnoOddities

### Example Command
To run AnnoOddities, use the command line interface with the required arguments for the genome FASTA file and GFF3 annotation file. For example:
```console
annooddities \
    --genome_fasta input_genome.fna \
    --gff3_file input_genome.gff
```

You can also specify optional parameters to customise the oddity detection thresholds. For example, to set custom thresholds for minimum exon length and maximum intron length (for plants), you can run:
```console
annooddities \
    --genome_fasta input_genome.fna \
    --gff3_file input_genome.gff \
    --min_exon_length 3 \
    --max_intron_length 10000
```

Or, if you want to restrict canonical intron motifs to only 'GT-AG', you can run:
```console
annooddities \
    --genome_fasta input_genome.fna \
    --gff3_file input_genome.gff \
    --canonical_intron_motifs 'GT-AG'
```

## Output

An example output directory structure:
```console
CMD: annooddities --genome_fasta input_genome.fna --gff3_file input_genome.gff

DIRECTORY: output_directory
├── input_genome.fna
├── input_genome.gff
├── input_genome.agat.log
├── output.agat_standardised.gff
├── output.agat_standardised.log
├── output.agat_standardised.agat.log
├── output.agat_sp_statistics.yaml
├── output.agat_sp_statistics.txt
├── output.agat_sp_statistics.log
├── output.agat_sp_statistics.json
├── output.agat_sp_statistics.toml
├── input_genome.fna.fai
├── output.gffread_table.log
├── output.gffread_table.tbl
├── output.mikado_tab_stats.log
├── output.mikado_tab_stats.tsv
├── output.mikado_stats.tsv
├── output.mikado_stats.yaml
├── output.mikado_summary_stats.tsv
├── output.mikado_stats.json
├── output.mikado_stats.toml
├── oddity_files
│   ├── output.AnnoOddities.has_inframe_stop.gff
│   ├── output.AnnoOddities.max_exon_length_gt_10000.gff
│   ├── output.AnnoOddities.min_exon_length_lte_5.gff
│   ├── output.AnnoOddities.exon_num_eq_1.gff
│   ├── output.AnnoOddities.exon_num_gt_1.gff
│   ├── output.AnnoOddities.is_fragment.gff
│   ├── output.AnnoOddities.not_has_start_codon.gff
│   ├── output.AnnoOddities.not_has_stop_codon.gff
│   ├── output.AnnoOddities.not_is_complete.gff
│   └── output.AnnoOddities.selected_cds_fraction_lte_0.3.gff
├── output.AnnoOddities.combined_statistics.yaml
├── output.AnnoOddities.combined_statistics.json
├── output.AnnoOddities.combined_statistics.toml
├── output.AnnoOddities.all_stats.tsv
├── output.AnnoOddities.gff
└── output.AnnoOddities.oddity_summary.txt

```

### AnnoOddities Summary
Provides both high-level summaries and more granular diagnostic reports to support manual review or automated pipelines.

Below is an example of the `AnnoOddities.oddity_summary.txt` output file:
```console
FILE: output.AnnoOddities.oddity_summary.txt
AnnoOddities                      output
exon_num == 1                     4520
exon_num > 1                      12900
five_utr_length > 10000           0
five_utr_num > 5                  0
three_utr_length > 10000          0
three_utr_num > 4                 0
not is_complete                   428
not has_start_codon               174
not has_stop_codon                280
is_fragment                       26
has_inframe_stop                  1
max_exon_length > 10000           8
max_intron_length > 120000        0
min_exon_length <= 5              308
0 < min_intron_length <= 5        0
selected_cds_fraction <= 0.3      158
canonical_intron_proportion != 1  0
only_non_canonical_splicing       0
suspicious_splicing               0
```

### AnnoOddities All Stats
Detailed per-transcript statistics (as detailed below) computed by AnnoOddities.

Below are the columns included in the `AnnoOddities.all_stats.tsv` output file:
```console
FILE: output.AnnoOddities.all_stats.tsv
1 transcript_id                             - transcript identifier
2 gene_id                                   - gene identifier
3 chromosome                                - chromosome name or number
4 start                                     - start position
5 end                                       - end position
6 region                                    - genomic region in the format chr:start-end
7 strand                                    - strand information
8 exon_num                                  - number of exons
9 exons                                     - list of exon positions
10 exon_lengths                             - list of exon lengths
11 max_exon_length                          - maximum exon length
12 min_exon_length                          - minimum exon length
13 total_exon_length                        - total length of all exons
14 cds_exon_num                             - number of CDS exons
15 cds_exons                                - list of CDS exon positions
16 cds_exon_lengths                         - list of CDS exon lengths
17 max_cds_length                           - maximum CDS exon length
18 min_cds_length                           - minimum CDS exon length
19 total_cds_length                         - total length of all CDS exons
20 cds_cdna_ratio                           - ratio of CDS length to cDNA length
21 intron_num                               - number of introns
22 introns                                  - list of intron positions
23 intron_lengths                           - list of intron lengths
24 max_intron_length                        - maximum intron length
25 min_intron_length                        - minimum intron length
26 total_intron_length                      - total length of all introns
27 five_utr_num                             - number of 5' UTRs
28 five_utr_length                          - total length of 5' UTRs
29 three_utr_num                            - number of 3' UTRs
30 three_utr_length                         - total length of 3' UTRs
31 has_start_codon                          - presence of start codon
32 has_stop_codon                           - presence of stop codon
33 is_complete                              - completeness of the transcript
34 is_fragment                              - whether the transcript is a fragment
35 has_inframe_stop                         - presence of in-frame stop codons
36 known_strand_junctions_dict              - dictionary of known strand junctions
37 known_strand_junctions_str               - string representation of known strand junctions
38 unknown_strand_junctions_dict            - dictionary of unknown strand junctions
39 unknown_strand_junctions_str             - string representation of unknown strand junctions
40 canonical_intron_count                   - count of canonical introns
41 non_canonical_intron_count               - count of non-canonical introns
42 suspicious_intron_count                  - count of suspicious introns
43 unknown_canonical_intron_count           - count of unknown canonical introns
44 unknown_non_canonical_intron_count       - count of unknown non-canonical introns
45 unknown_suspicious_intron_count          - count of unknown suspicious introns
46 unknown_predicted_strand                 - predicted strand for unknown introns
47 canonical_intron_proportion              - proportion of canonical introns
48 only_non_canonical_splicing              - whether only non-canonical splicing is present
49 suspicious_splicing                      - presence of suspicious splicing
50 matched_oddities                         - matched oddities
```

### AnnoOddities Combined Statistics
Consolidated statistics from AGAT, Mikado, GFFread, and AnnoOddities.
```console
FILE: output.AnnoOddities.combined_statistics.yaml
FILE: output.AnnoOddities.combined_statistics.json
FILE: output.AnnoOddities.combined_statistics.toml
```

### AnnoOddities Individual Oddity GFFs
Separate GFF3 files for each detected oddity, facilitating targeted review and analysis.
```console
DIRECTORY: oddity_files
FILE: output.AnnoOddities.has_inframe_stop.gff
FILE: output.AnnoOddities.max_exon_length_gt_10000.gff
FILE: output.AnnoOddities.min_exon_length_lte_5.gff
FILE: output.AnnoOddities.exon_num_eq_1.gff
FILE: output.AnnoOddities.exon_num_gt_1.gff
FILE: output.AnnoOddities.is_fragment.gff
FILE: output.AnnoOddities.not_has_start_codon.gff
FILE: output.AnnoOddities.not_has_stop_codon.gff
FILE: output.AnnoOddities.not_is_complete.gff
FILE: output.AnnoOddities.selected_cds_fraction_lte_0.3.gff
```

### AnnoOddities GFF File
A GFF3 file containing all transcripts annotated with their respective oddities. The oddities are included in the attributes column for easy identification, appeneded to the 'Note' field to the respective transcript entries.
```console
FILE: output.AnnoOddities.gff
```


## License

MIT
