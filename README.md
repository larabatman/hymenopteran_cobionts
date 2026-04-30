# hymenopteran_cobionts
A cross-species pipeline for mining cobionts in Hymenopteran genomes

## OVERVIEW

The project targets 271 hymenopteran species with publicly available Darwin Tree of Life (DToL) genome assemblies. The goal is to build a cobiont-species matrix from host-associated organisms (bacteria, fungi, protist) embedded in or co-sequenced with the host genome. This matrix serves as a covariate for genome evolution studies conducted by collaborating Sanger groups. 

The work is organized in two phases: 
1. Exploratory phase (5 test species): manual pipeline development and evaluation of orthogonal evidence layers for host-cobiont separation, documented below.
2. Production phase (271 species): batch runs using the sanger-tol/metagenomeassembly pipeline with additional binners, followed by RNA-seq evidence integration where data is available 

### Test species

| Species| Wolbachia evidence | Hi-C available | Sample type |
| ---- | ---- | ---- | ----- |
|Bombus pratorum | No | No | Head & Thorax adult female|
| Bombus terrestris | No | Yes | Head & Thorax adult female |
| Lasioglossum morio | No | No | Whole body adult female |
| Lasioglossum pauxillum| Yes | Yes | Whole body adult female |
| Nomada fabricana | Yes | Yes | Whole body adult female |

## REPOSITORY STRUCTURE
scripts/
species/        # curated species list (input)
data/           # large or derived datasets (gitignored)
logs/           # HPC logs (gitignored)

## SCRIPTS

### Phase 0: Dataset exploration, NCBI-ENA metadata linkage

The project involves 271 hymenopteran species. Manually checking, species by species, whether PacBio, Hi-C, or RNA-seq data exist, and whether they are coherent across studies or individuals is error-prone, irreproducible and unrealistic. This difficulty is accentuated by the fact that links between genome assemblies, sequencing runs, studies and biological individuals are inconsistently curated across NCBI and ENA. Moreover, assessing the biological context (sex, tissue of origin, developmental stage) introduces an additional layer of complexity.
Using assemblies as anchor points, data is discovered independently and reconnected when metadata agrees. BioSamples are then queried from NCBI to recover biological metadata, following parent-child relationships when necessary. All following scripts were executed through a wrapper script that submitted each task as independent SLURM job SBATCH directives.

#### 01a_map_gca_to_runs.py
As such, a Python script allows to retrieve raw sequencing dat aassociated with hymenopteran genome assemblies (01a_map_gca_to_runs.py). For each species listed in species/Hymenopteran_genomes.csv, it:

- Resolves assembly metadata from NCBI
- Queries ENA for sequencing runs using evidence hierarchy
- Classifies seuqencing technologies: PacBio, Hi-C, RNA-seq, others
- Emits a normalized TSV summary for downstream filtering

Progress and API warning are written to stderr, and TSV results are written to stdout. It is designed to run on HPC systems, with rate limiting, retry logic and enforced IPv4 networking.

The input is a CSV with at least 'species' and 'accession' as (GCA_*). 
The output is a TSV, with one row per assembly, containing:

- Assembly identifiers: taxid, BioProject, BioSample
- Evidence lebel used to retrieve runs
- Presence/ absence flags for PacBio, Hi-C and RNA-seq
- Lists of run accessions and Biosamples
- Detection of shared PacBio/ Hi-C BioSamples (same individual)
- Other sequencing strategies present

##### Output columns
| Column                        | Description                                                                                    |
| ----------------------------- | ---------------------------------------------------------------------------------------------- |
| `species`                     | Species name (from input CSV)                                                                  |
| `assembly_accession`          | Genome assembly accession (GCA_*)                                                              |
| `taxid`                       | NCBI taxonomic identifier resolved from the assembly                                           |
| `bioproject`                  | BioProject associated with the assembly (NCBI assembly metadata)                               |
| `biosample`                   | BioSample associated with the assembly (NCBI assembly metadata)                                |
| `run_mapping_level`           | Evidence level used to retrieve ENA runs (`BIOPROJECT`, `BIOSAMPLE`, `TAXID_FALLBACK`, `NONE`) |
| `has_pacbio`                  | Whether PacBio sequencing runs were detected                                                   |
| `has_hic`                     | Whether Hi-C sequencing runs were detected                                                     |
| `has_rnaseq`                  | Whether RNA-seq runs were detected                                                             |
| `pacbio_runs`                 | Comma-separated list of PacBio run accessions                                                  |
| `hic_runs`                    | Comma-separated list of Hi-C run accessions                                                    |
| `rnaseq_runs`                 | Comma-separated list of RNA-seq run accessions                                                 |
| `pacbio_biosamples`           | BioSample accessions associated with PacBio runs                                               |
| `hic_biosamples`              | BioSample accessions associated with Hi-C runs                                                 |
| `shared_pacbio_hic_biosample` | Boolean indicating whether PacBio and Hi-C data share at least one BioSample (same individual) |
| `other_sequencing_methods`    | Other ENA library strategies detected (excluding PacBio, Hi-C, RNA-seq)                        |

#### 01b_summarize_data_types.R
This script summarizes sequencing data availability across species using the ouput of the GCA to ENA inventory. As primary signals, it focuses on PacBio and Hi-C availability and uses BioSample sharing as a biological coherence indicator.

The input is data/data_inventory.tsv previously produced.
The output is a printed summary table as well as TSV files written to data/inventory_analysis for downstream inspection. 

#### 01c_map_biological_context.py
This script annotates ENA or assembly-linked BioSample accessions with biological metadata (sex, tissue. developmental stage, description) using NCBI BioSample as source. 
It is designed to work downstream of the sequencing inventory, to recover biological context at the individual level. 

The input is a TSV file with columns 'pacbio_biosamples' and 'hic_biosamples' required.
The output is a TSV with on row per input row and BioSample. Original columns are preserved and annotation fields are appended. 

##### Appended columns:
| Column             | Description                             |
| ------------------ | --------------------------------------- |
| `sample_accession` | BioSample accession being annotated     |
| `sex`              | Biological sex (or `UNKNOWN`)           |
| `tissue`           | Tissue / organism part (or `UNKNOWN`)   |
| `dev_stage`        | Developmental stage                     |
| `description`      | BioSample title / free-text description |

#### 01d_explore_biosample_metadata.R
This script performs an exploratory analysis of the BioSample-level annotations, reflecting the availability and structure of BioSample metadata and not necessarily the biological composition of the underlying genome assemblies.

The input is data/annotation/data_inventory_biosample_annotation.tsv.
The outputs are multiple TSV files written to data/annotation/exploration/

#### 01e_filter_sequencing_biological_context.R
This script constructs a biologically and technically restricted BioSample dataset for downstream anayses. It selects BioSamples that are technically coherent, sharing PacBio and Hi-C evidence from the same individual, and biologically constrained to whole-body or abdomen female tissues.
The input is data/annotation/data_inventory_biosample_annotation.tsv.
The output is data/annotation/restricted_dataset.tsv, a filtered BioSample-level dataset, as well as printed summary statistics. 

### Phase 1: Exploratory pipeline, host-cobiont separation (5 test species)

The main logic of this section is to separate insect (host) contigs from co-ocurring cobiont contigs embedded in the sequencing data. By collecting several orthogonal layers of evidence, a model of th ehost genome is inferred and substracted from the whole assembly. The remainder can be further analyzed for cobiont grenomes through binning, while the host genome can be mined for integrated cobionts. 

Scripts are organized as stages and accept a uniform interface: 

sbatch STAGE_scripts.sh ${SPECIES} ${ASM_MODE}

Where ${SPECIES} is the species name (e.g., Bombus_terrestris) and ${ASM_MODE} is either:
- hic when Hi.C data is available
- bp when Hi-C data is not available

#### Data retrieval 

To retrieve sequencing runs from ENA: 

for run in ERR8702822 ERR8702823 ERR9081702; do
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp" \
    | tail -n +2 | cut -f2 | tr ';' '\n' | sed 's#^#https://#' \
    | xargs -n1 wget -c
done

#### Stage A: Assembly and Quality Control

##### A1: Reads QC

###### A1a: HiFi read QC


Scripts:
A1a_hifi_QC.sh and A1_reads_diagnostics.R filter the raw sequencing HiFi reads prior to assembly and is meant to be launched in parallel with A1b for Hi-C QC when relevant. 

HiFi reads:
- Raw statistics are extracted with seqkit stats
- Length-filtered with fastplong: minimum 3 kbp, maximum 2X raw read N50 kbp (auto-estimated or manually overriden for species flagged in curation)
- Adapter trimming with fastplong
- Post-filter statistics and FastQC report are generated
- Raw and filtered read length distributions are extracted for the R diagnostics script 

R diagnostics: produces before and after read length histograms with bindiwth set at 0.5 kbp for cross-species comparability, and prints N50, read count, and summary statistics. 

Usage: 
sbatch A1a_hifi_QC.sh <species> [hifi_max_len]

Inputs:
- PacBio HiFi reads: reads/pacbio_hifi/${SPECIES}/*.fastq.gz

Outputs: 
- results/${SPECIES}_stages/read_qc/hifi_qc/hifi_raw_stats.tsv
- results/${SPECIES}_stages/read_qc/hifi_filtered/hifi_filtered.fastq.gz
- results/${SPECIES}_stages/read_qc/hifi_qc/hifi_filtered_stats.tsv
- results/${SPECIES}_stages/read_qc/hifi_qc/fastplong_report.{html,json}
- results/${SPECIES}_stages/read_qc/hifi_qc/hifi_length_distribution_before_after.png
- results/${SPECIES}_stages/read_qc/hifi_qc/hifi_length_cutoffs.tsv

Tools: fastplong, FastQC, seqkit, R (ggplot2)

###### A1b: Hi-C reads QC

Scripts: 
A1b_hic_qc.sh quality controls and filters Hi-C reads, meant to be launched in parallel with A1a. 

Hi-C reads:
- Pre-cleaning FastQC and seqkit stats
- Adapter trimming and quality filtering with fastp: minimal quality of called base Q20, minimal length of 50 bp, PE adapter detection 
- Post-cleaning FastQC and seqkit stats

Usage:
sbatch A1b_hic_qc.sh <species>

Inputs:
- Hi-C reads: reads/hic/${SPECIES}/hic_R1.fastq.gz, hic_R2.fastq.gz

Outputs: 
- reads/${SPECIES}/hic_clean/hic_R{1,2}.clean.fastq.gz (if Hi-C)
- results/${SPECIES}_stages/read_qc/hic_qc/fastè_hic_report.{html,json}

Tools: fastp, FastQC, seqkit

Note: When multiple Hi-C technical replicates exist, they are concatenated into a single R1/R2 pair using concatenate_hic_replicates.sh prior to this stage, since additional information increases contact density and reduces noise. 


##### A2: hifiasm assembly 

Scripts: 
A2_hifiasm.sh assembles filtered HiFi reads, optionally integrating Hi-C contact infromation for phasing. 
- ASM_MODE=bp: HiFi-only assembly, produces asm.bp.p_ctg.gfa
- ASM_MODE=hic: HiFi + Hi-C assembly, produces asm.hic.p_ctg.gfa

The primary GFA is converted to FASTA by extraction segment lines (with S records). Basic assembly statistics are generated with seqkit stats.

Usage: 
sbatch A2_hifiasm.sh <species> <asm_mode>

Inputs:
- Filtered HiFi: results/${SPECIES}_stages/read_qc/hifi_filtered/hifi_filtered.fastq.gz
- (Optional) Cleaned Hi-C: reads/${SPECIES}/hic_clean/hic_R{1,2}.clean.fastq.gz

Outputs: 
- assemblies/hifiasm/${SPECIES}/asm.{bp,hic}.p_ctg.{gfa,fasta}
- assemblies/hifiasm/${SPECIES}assembly_basic_stats.tsv

Tools: hifiasm, seqkit

##### A3: Read mapping validation

Scripts:
A3_read_mapping.sh maps filtered HiFi reads back to the assembly to assess mapping rate and per-contig coverage. 
1. Aligns with minimap2 (-ax map-hifi preset), pipes to samtools sort for coordinate-sorted BAM
2. Indexes the BAM with samtools index
3. Generates mapping summaries: samtools flagstat and samtools stats
4. Computes per-contig coverage with samtools coverage
5. Calculates length-weighted genome-wide mean depth

Usage:
sbatch A3_read_mapping.sh <species> <asm_mode>

Inputs:
- Assembly: assemblies/hifiasm/${SPECIES}/asm.{bp,hic}.p_ctg.fasta
- Filtered HiFi: results/${SPECIES}_stages/read_qc/hifi_filtered/hifi_filtered.fastq.gz

Outputs: 
- results/${SPECIES}_stages/assembly_qc/reads.bam + .bai
- results/${SPECIES}_stages/assembly_qc/mapping_flagstat.txt
- results/${SPECIES}_stages/assembly_qc/mapping_stats.txt
- results/${SPECIES}_stages/assembly_qc/coverage_per_contig.tsv
- results/${SPECIES}_stages/assembly_qc/coverage_summary.tsv

Tools: minimap2, samtools

##### A4: Assembly QC sumary

Scripts: 
A4_summarize_QC.sh and A4_summarize_QC.R parse stage A outputs into a standardized QC table. The R script merges assembly size, contig count, average contig length, mapping rate and mean depth into a species-level summary. A QC status flag is assigned (PASS or FAIL_low_mapping if mapping rate < 95%). It also updates a global summary file (results/global_qc_summary.tsv) across species, upserting by species name. 

Usage:

A4_summarize_QC.sh <species>

Output: 
- results/${SPECIES}_stages/assembly_qc/assembly_qc_status.tsv
- results/global_qc_summary.tsv (global, appended and updated)

| Column | Description|
| ----- | ----- |
| species| Species name |
| assembly_size | Total assembly length (bp) |
| n_contigs | Number of contigs |
| avg_len | Average contig length |
| mapped_pct| % reads mapped from flagstat |
| mean_depth | Length-weighted mean depth |
| status| PASS/FAIL_low_mapping |

Tools: R (dplyr, readr, stingr)

#### Stage B: Orthogonal evidence, coverage-based host backbone 

##### B1: Per-contig GC, length and coverage

Scripts: 
B1_gc_cov.sh computes per-contig sequence composition and coverage from the assembly and existing BAM. 
1. GC and length: seqkit fx2tab extracts contig name, length and GC fraction into gc_len.tsv
2. Coverage: samtools coverage on the BAM from A3 writes to coverage.tsv, which is parsed to retain contig name and mean depth to cov.tsv
3. Merge: an awk join on contig name produces the final gc_cov.tsv

Usage_
sbatch B1_gc_cov.sh <species> <asm_mode>

Inputs:
- Assembly: assemblies/hifiasm/${SPECIES}/asm.{bp,hic}.p_ctg.fasta
- BAM from A3: results/${SPECIES}_stages/assembly_qc/reads.bam

Output:
- results/${SPECIES}_stages/gc_cov/gc_cov.tsh

| Column | Description |
| ----- | ----- |
| contig | Contig ID (ptg*) |
| len | Contig length (bp) |
| gc | GC content (%) |
| mean_cov | Mean read depth |

Tools: seqkit, samtools 

##### B2: Coverage modeling and host backbone definition

Scripts: 
B2_identify_dominant_cov_mode.sh and B2_identify_dominant_cov_mode.R define the host backbone by modeling the coverage distribution with a length-weighted Median Absolute Deviation (MAD) approach. The rationale is that the majority of assembled DNA is host, such that the coverage median approximates host depth, and MAD provides a robust spread estimate. 

Method: 
1. Log-transform mean coverage (log10 with e = 1e-6 for zeros)
2. Compute length-weigthed median (host coverage center)
3. Compute length-weighted MAD with normal consistency constant (1.4826)
4. Classify contigs by absolute z-score:
- |z| < 2: host_like
- 2 ≤ |z| < 4: ambiguous
- |z| ≥ 4: coverage_outlier
- |z| > 6: additional is_extreme flag
5. Direction is annotated as high, low or neutral 

Usage: 
sbatch B1b_run_identify_dominant_cov_mode.sh <species>

Input:
- results/${SPECIES}_stages/gc_cov/gc_cov.tsv

Outputs:
- results/${SPECIES}_stages/host_backbone/coverage_classification.tsv is the full table with z-scores and classes
- results/${SPECIES}_stages/host_backbone/host_backbone.tsv are the host-like contigs only
- results/${SPECIES}_stages/host_backbone/coverage_backbone_summary.tsv
- results/${SPECIES}_stages/host_backbone/bloplot_covergae_class.png for a base GC-vs-coverage blobplot
- results/${SPECIES}_stages/host_bakcbone/coverage_histogram.png for coverage distirbution with MAD median line 

| Column | Description |
| ----- | ----- |
| species | Species name |
| n_total_contigs| Total contigs |
| host_median_logcov| Length-weighted median log10(coverage) |
| host_mad_logcov| length-weighted MAD of log10(coverage) |
| n_backbone | Number of host-like contigs |
| percent_backbone| % of assembly bp in host backbone |

Tools: R (dplyr, readr, matrixStats, ggplot2)

#### Stage C: Orthogonal evidence, BUSCO marker gene analysis

##### C1: BUSCO runs

Scripts:
C1_busco.sh runs BUSCO in genome mode against three lineage databases to detect conserved marker genes on each contig. The logic is that Hymenoptera and Arthropoda hist anchor host contigs while Bacteria hits flag potential cobionts. Contigs carrying both sets of marker are candidates for lateral gene transfer or host-integrated cobionts. 
Three lineage databases are used: 
- hymenoptera_odb10 as strong host anchors
- arthropoda_odb10 as supportive host anchors 
- bacteria_odb10 as anti-host or cobiont signal

Usage: 
sbatch C1_busco.sh <species> <asm_mode>

Input:
- Assembly: assemblies/hifiasm/${SPECIES}/asm.{bp,hic}.p_ctg.fasta

Outputs:
- Standard BUSCO output directories under results/${SPECIES}_stages/busco_anchor

Tools: BUSCO

##### C2: Unified BUSCO analysis

Scripts:
C2_run_busco_analysis.sh and C2_busco_analysis.R consolidate contig-level BUSCO collapse, anchor validation against the coverage backbone and all blobplot generation. 

1. Contig-level collapse: reads the three full_table.tsv files using a robust parser that handles broken TSV formatting. For each contig, computes total BUSCO hits, per-lineage counts and presence flags, status breakdown (Complete, Fragmented, Duplicated) per lineage, and functional annotations (BUSCO gene descriptions per lineage). The result is merged with the coverage classification table from B2.

2. Anchor validation: assigns anchor tiers to each contig based on BUSCO signal:
- strong_hymenoptera for contigs that have hymenoptera BUSCO hits
- supportive_arthtropoda for contigs that have arthropoda but not hymenoptera hits
- none for contigs without host markers
It also cross-tabulates anchor type and coverage class with bacterial signal and coverage class to assess whether BUSCO-positive host contigs fall within the expected host coverage range and whether bacterial BUSCO contigs are enriched among coverage outliers. 

3. Blobplots: all plots use log10(mean_cov) on y-axis and GC% on x-axis, with point size proportional to log10(contig_length). Filters applied: mean_cov > 0, len > 0 and GC between 13 and 80 %. 

| Plot | Description |
| ----- | ----- |
| bloplot_hym_coverage.png| Coverage class fille + hymeoptera BUSCO-positive contigs as ring overlay |
| blobplot_arth_coverage.png | Coverage class fill + arthropoda BUSCO-positive contigs as ring overlay |
| blobplot_bact_coverage.png | Coverage class fill + bacteria BUSCO-positive contigs as ring overlay |
| ${SPECIES}_blobplot_hym_gradient.png | Continuous gradient by hymenoptera BUSCO density (hits per Mb) |
| ${SPECIES}_blobplot_arth_gradient.png | Continuous gradient by arthropoda BUSCO density (hits per Mb) |
| ${SPECIES}_blobplot_bact_gradient.png | Continuous gradient by bacteria BUSCO density (hits per Mb) |

Usage: 
sbatch C2_run_busco_analysis.sh <species>

Inputs: 
- BUSCO full_table.tsv for each lineage, from C1
- results/${species}_stages/host_backbone/coverage_classification.tsv

Outputs: all in results/${SPECIES}_stages/busco_analysis/
- busco_gene_level.tsv: every BUSCO hit with contig, busco_id, status, lineage, description
- busco_per_contig_summary.tsv: contig-level counts, statuses, annotations, merged with coverage
- busco_mixed_hym_bact.tsv: subset of contigs with both hymenoptera and bacteria hits
- buscho_ancho_tiers.tsv: anchor type classification per contig
- busco_anchor_coverage_distribution.tsv: cross-tabulation of anchor type and coverage class
- busco_bacteria_coverage_distirbution.tsv: cross-tabulation of anchor type and coverage class
- busco_anchor_full_eval.tsv: all contigs with anchor type, bacterial signal and coverage feature
- blobplot_summary.tsv: contig counts with total, hym, arth, bact, mixed

Tools: R (dplyr, readr, stringr, ggplot2)

#### Stage D: Repetitive element annotation 

Repetitive elements are investigated to explain the striking heterogeneity observed in contig coverage and GC content. Repeat-rich regions can inflate coverage estimates and cofound host-cobiont separation. 

##### D1: EDTA

Scripts:
D1_edta.sh runs the EDTA pipeline for whole-genome TE annotation and produces a classified TE library and genome-wide GFF3 annotation. 

Usage: 
sbatch D1_edta.sh <species> <asm_mode>

Parameters: --species others, --sensitive 1, --anno 1

Outputs: results/${SPECIES}_stage/edta
- .TElib.fa: TE library
- .TEanno.gff3: annotated TE library
- summary statistics

Tools: EDTA

##### D2: RepeatMasker + TRF

Scripts:
D2_rm_trf.sh runs RepeatMasker with the EDTA-derived TE library for classified repeat masking and TRF for microsatellite and tandem repeat detection. These complement EDTA by capturing simple repeats and low-complexity regions not covered by the TE library. 
RepeatMasker is run with -xsmall sor softmasking, -gff and -no_is to skip bacterial IS check. 
TRF uses default parameters with match=2, mismatch=7, delta=7, PM=80, PI=10, minscore=50 and maxperiod=500.

Usage: 
sbatch D2_rm_trf.sh <species> <asm_mode>

Outputs: 
- results/${SPECIES}_stages/repeat_masking/repeatmasker/ RepeatMasker output (.out, .gff, masked FASTA)
- results/${SPECIES}_stages/repeat_masking/trf/ TRF .dat output

Tools: RepeatMasker, TRF

##### D3: Unified repeat analysis

Scripts: 
D3_run_repeat_analysis.sh and D3_repeat_analysis.R consolidate EDTA GFF3 parsing, RepeatMasker and TRF parsing, blobplot generation and summary tables. 

1. Parsing: reads the EDTA GFF3 for TE intervals, RepeatMasker .out for classified repeats (TEs, simple repeats, low-complexity, satellites), and TRF .dat for tandem repeats. All overlapping intervals are merged per-contig to avoid double counting. RepeatMasker hits are classified into broader categories: te, simple_repeat, low_xomplexity, satellite, rna.

2. Blobplots: all blobplots use the same axes and styling as stage C blobplots (log10(mean_cov) as y-axis, GC% as x-axis and size by log10(contig_length)). Gradient coloring by repeat fraction on 0-1 scale

| Plot | Description |
| ----- | ----- |
| ${SPECIES}_blobplot_te.png | TE fracion from EDTA |
| ${SPECIES}_blobplot_non_te_repeats.png | Non-Te repeats|
| ${SPECIES}_blobplot_total_repeats.png | Total repeat fraction from all sources

Usage:
sbatch D3_run_repeat_analysis.sh <species> <asm_mode>

Inputs
- EDTA GFF3: results/${SPECIES}_stages/edta/${èREFIX}.mod.EDTA.TEanno.gff3
- RepeatMasker .out: results/${SPECIES}_stages/repeat_masking/repeatmasker/${PREFIX}.out
- TRF .dat: results/${SPECIES}_stages/repeat_masking/trf/${PREFIX}.2.7.7.80.10.50.500.dat
- Assembly: assemblies/hifiasm/${SPECIES}/asm.{bp,hic}.p_ctg.fasta
- results/${SPECIES}_stages/host_backbone/coverage_classification.tsv from B1b

Outputs: all in results/${SPECIES}_stages/repeat_analysis/
- contig_repeat_coverage.tsv: per-contig bp and fractions for each category
- repeat_coverage_distirbution.tsv: cross-tabulation of repeat statues (high/low) and coverage class
- repeat_genome_summary.tsv: genome-wide repeat percentages by category

The contig_repeat_coverage.tsv columns include:

| Column | Description |
| ----- | ----- |
| contig| Contig ID |
| length | Contig length (bp) |
|  edta_te_bp | TE bp from EDTA GFF3 (merged intervals) |
| rm_te_bp | TE bp from RepeatMasker |
| simple_repeat_bp | Simple repeat bp (RepeatMasker) |
| low_complexity_bp | Low-complexity bp (RepeatMasker) |
| satellite_bp | Satellie bp (RepeatMasker) |
| tandem_repeat_bp | Tandem repeats bp (TRF) |
| non_te_repeat_bp | Non-TE repeat bp |
| total_repeat_bp | Total repeat bp, all sources merged |
| edta_te_fraction | TE fraction from EDTA |
| non_te_repeat_fraction | Non-TE repeat fraction |
| total_repeat_fraction | Total repeat fraction |
| flagged_repetitive| yes/no based on threshold (default at ≥50 %) |

Tools: R (dplyr, readr, ggplot2)

# TOOLS AND VERSIONS

| Tool | Version | Source | 
| Python |  3.9.5 | GCCcore-10.3.0 |
| R | 4.2.1 | foss-2021a |
| hifiasm | 0.16.1 |  GCCcore-10.3.0 |
| SeqKit | 2.6.1 | module |
| Kraken2 | 2.1.2 | gompi-2021a |
| minimap2 | 2.20 | GCCcore-10.3.0 |
| SAMtools| 1.13 | GCC-10.3.0 |
| BUSCO | 5.4.2 | foss-2021a
| FastQC | 0.11.9 | Java-11 |
| fastp | 0.23.4 | GCCC-10--3-0 |
| EDTA | 2.2 | Apptainer container |
| RepeatMasker | | 4.1.5 | foss-2021a |
| TRF | 4.09.1 | GCC-10.3.0 |
| Nextflow | | isolated environment, Java 17.0.6 |
 