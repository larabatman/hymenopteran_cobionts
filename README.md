# hymenopteran_cobionts
A cross-species pipeline for mining cobionts in Hymenopteran genomes

## REPOSITORY STRUCTURE
scripts/
species/        # curated species list (input)
data/           # large or derived datasets (gitignored)
logs/           # HPC logs (gitignored)

## SCRIPTS

### DATASET EXPLORATION: NCBI-ENA metadata linkage scripts

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

## TOOLS VERSION
module Python/3.9.5-GCCcore-10.3.0
module R/4.2.1-foss-2021a