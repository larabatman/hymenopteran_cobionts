# hymenopteran_cobionts
A cross-species pipeline for mining cobionts in Hymenopteran genomes

## REPOSITORY STRUCTURE
scripts/
species/        # curated species list (input)
data/           # large or derived datasets (gitignored)
logs/           # HPC logs (gitignored)

## SCRIPTS

### NCBI-ENA metadata linkage script
The project involves 272 hymenopteran species. Manually checking, species by species, whether PacBio, Hi-C, or RNA-seq data exist, and whther they are coherent across studies or individuals is error-prone, irreproducible and unrealistic. Furthermore, links between assemblies, reads and studies are inconsistently curated across NCBI ad ENA. Using assemblies as anchor poinrs, data is discovered independently and reconnected when metadata agrees.
As such, a Python script allows to retrieve raw sequencing dat aassociated with hymenopteran genome assemblies (ncbi_GCA_query.py). For each species listed in species/Hymenopteran_genomes.csv, it:
- Queries NCBI Assembly via GCA by resolving the assembly accession to NCBI taxonomic ID (taxid) and Assembly-associated BioProject
- Queries ENA independently via taxid to retrieve all raw seequencing runs for the species, while identifying data types: PacBio, Hi-C, RNA-seq. 
- Assesses the coherence across datasets by checking the overlap of BioProjects and BioSamples, to determine whether different data types are likely to originate from the same biological context (same individual or study)
- Outputs a TSV summary to stdout, with one row per species, allowing filtering and analysis downstream 

#### Usage 
python3 scripts/ncbi_GCA_query.py > species_raw_data_summary.tsv

Progress and API warning are written to stderr, and TSV results are written to stdout. It is designed to run on HPC systems, with rate limiting, retry logic and enforced IPv4 networking. 

#### Output columns
| Column                        | Description                                |
| ----------------------------- | ------------------------------------------ |
| `species`                     | Species name                               |
| `assembly_accession`          | GCA accession                              |
| `taxid`                       | NCBI taxonomic ID                          |
| `assembly_bioproject`         | BioProject associated with the assembly    |
| `has_pacbio`                  | PacBio data present                        |
| `has_hic`                     | Hi-C data present                          |
| `has_rnaseq`                  | RNA-seq data present                       |
| `pacbio_runs`                 | PacBio run accessions                      |
| `hic_runs`                    | Hi-C run accessions                        |
| `rnaseq_runs`                 | RNA-seq run accessions                     |
| `pacbio_projects`             | BioProjects for PacBio runs                |
| `hic_projects`                | BioProjects for Hi-C runs                  |
| `pacbio_biosamples`           | BioSamples for PacBio runs                 |
| `hic_biosamples`              | BioSamples for Hi-C runs                   |
| `shared_pacbio_hic_project`   | Shared BioProject between PacBio and Hi-C  |
| `shared_pacbio_hic_biosample` | Shared BioSample between PacBio and Hi-C   |
| `assembly_project_in_pacbio`  | Assembly BioProject appears in PacBio data |
| `assembly_project_in_hic`     | Assembly BioProject appears in Hi-C data   |

## TOOLS VERSION
module Python/3.9.5-GCCcore-10.3.0
module R/4.2.1-foss-2021a