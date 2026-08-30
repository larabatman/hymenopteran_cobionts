#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=metaasm_full_binners
#SBATCH --output=logs/metaasm_full_binners_%j.out
#SBATCH --error=logs/metaasm_full_binners_%j.err

# Pilot job for the NextFlow sanger-tol/metagenomeassembly pipeline 
# SLURM directives above define the allocation that runs NextFlow, and then it submits its own SLORM hobs with specific resources from the NextFlow config file defined below
# Usage:
# sbatch run_full_binner_hic.sh species


set -euo pipefail

SPECIES="$1"
# Run tag for multiple attempts on the same species
RUN_TAG="attempt_hic"
# Working directory
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
# Path to the fork of the pipeline
PIPELINE_DIR="${PROJECT_ROOT}/metagenomeassembly"

# Nextflow specificities: working directory that contains
# nf/ with the input YAML and nextflow config files that are created below
# work/ with the process directories for each task
# .nextfow/ internal state that allows the -resume
# .apptainer_cache for the downladed images
# logs/ with the pipeline reports as timelines and traces
RUN_DIR="$PROJECT_ROOT/runs/${SPECIES}/${RUN_TAG}"
# The publishing directory is separated from the work directory, and that is where the final results are 
OUTDIR="$PROJECT_ROOT/results/${SPECIES}/full_binners"

# Input directories: hifi fasta reads and hic.cram
HIFI_FASTA="${PROJECT_ROOT}/reads/processed/${SPECIES}/hifi_filtered.fasta.gz"
HIC_CRAM="${PROJECT_ROOT}/reads/processed/${SPECIES}/hic.cram"

# Use scratch and another for Apptainer
export TMPDIR="/scratch/$USER/tmp_${SLURM_JOB_ID}"
export APPTAINER_TMPDIR="/scratch/$USER/apptainer_${SLURM_JOB_ID}"

# Databases directory
# GTDB-Tk release r226
GTDB_DB="${PROJECT_ROOT}/databases/gtdb_r226/release226"
HMM_TIGRFAM="${GTDB_DB}/markers/tigrfam/tigrfam.hmm"
HMM_PFAM="${GTDB_DB}/markers/pfam/Pfam-A.hmm"
GTDB_BAC120_META="${PROJECT_ROOT}/databases/gtdb_r226/bac120_metadata_r226.tsv"
GTDB_AR53_META="${PROJECT_ROOT}/databases/gtdb_r226/ar53_metadata_r226.tsv"
# CheckM2 with UniRef100 proteins  
CHECKM2_DB="${PROJECT_ROOT}/databases/checkm2_102/CheckM2_database/uniref100.KO.1.dmnd"
# geNomad
GENOMAD_DB="${PROJECT_ROOT}/databases/genomad_db"

# Modules
module purge
module load Java/17.0.6

mkdir -p logs "$RUN_DIR"/{nf,logs,.nextflow,.apptainer_cache,work} "$OUTDIR" "$TMPDIR" "$APPTAINER_TMPDIR"

# Use scratch and another for Apptainer
export TMPDIR="/scratch/$USER/tmp_${SLURM_JOB_ID}"
export APPTAINER_TMPDIR="/scratch/$USER/apptainer_${SLURM_JOB_ID}"
# Nextflow environment:
# Define NXF_HOME as the place to store metadata to isolate them between runs
export NXF_HOME="$RUN_DIR/.nextflow"
# Define NXF_WORK as the working directory for each proces
export NXF_WORK="$RUN_DIR/work"
# Define NXF_CACHE_DIR 
export NXF_CACHE_DIR="$RUN_DIR/.nextflow/cache"
# Define NXF_OPTS with Java options for the NextFlow process, not the tools
export NXF_OPTS='-Xms4g -Xmx16g'
# Define APPTAINER_CACHEDIR for container images to share images between runs
export APPTAINER_CACHEDIR="$RUN_DIR/.apptainer_cache"
export NXF_SINGULARITY_CACHEDIR="$RUN_DIR/.apptainer_cache"

# Write input YAML: 
# id are the sample identifier to name output files, species here
# pacbio: fasta: with the list of HiFi reads
# hic: cram with the HiC cram file path en hic: enzmyes: with the restriction enzymes used in the Hic protocol
cat > "$RUN_DIR/nf/input.yaml" <<EOF
id: ${SPECIES}

pacbio:
  fasta:
    - ${HIFI_FASTA}

hic:
  cram:
    - ${HIC_CRAM}
  enzymes:
    - DpnII
    - HinfI
    - MseI
    - DdeI
EOF

# Write Nextflow config file 
# Generated each run and passed through NextFlow nextflow -c with the resource allocation for each process, and the executor settings
cat > "$RUN_DIR/nf/nextflow.config" <<EOF

process.executor = 'slurm'
process.queue    = 'pibu_el8'
workDir          = System.getenv('NXF_WORK')

process {
  scratch   = false
  maxForks  = 4
  errorStrategy = 'finish'

  beforeScript = '''
  set -e
  export TMPDIR=\$PWD/tmp
  mkdir -p \$TMPDIR
  '''

  withName: /.*METAMDBG_ASM/ {
    cpus   = 8
    memory = '100 GB'
    time   = '120h'
    // Pass thread count
    ext.args = '--threads 8'
  }

  withName: /.*PYRODIGAL/ {
    cpus   = 8
    memory = '32 GB'
    time   = '4h'
  }

  withName: /.*FASTXALIGN_MINIMAP2ALIGN/ {
    cpus   = 8
    memory = '60 GB'
    scratch = true
  }

  withName: /.*MINIMAP2_INDEX/ {
    cpus          = 4
    // Memory scales with assembly size, which varies across species: 80 GB * attempt gives 80, then 160, then 240 GB across three retries.
    memory        = { 80.GB * task.attempt }
    time          = '4h'
    // Setting a retry to not fail the pipeline immediately
    errorStrategy = 'retry'
    maxRetries    = 3
  }

  withName: /.*INFERNAL_CMSEARCH/ {
    cpus          = 6
    // Here memory also scales with assembly size: 40 GB * attempt gives 40, then 80, then 120 GB across retries.
    memory        = { 40.GB * task.attempt }
    time          = '8h'
    errorStrategy = 'retry'
    maxRetries    = 3
  }

  withName: /.*CHECKM2_PREDICT/ {
    cpus   = 4
    memory = '96 GB'
    time   = '24h'
    // Apptainer options: --writable-tmpfs to let CheckM2 write on top of the image, files live in memory
    // -B /scratch to map /scratch inside the container to write the temporary files 
    containerOptions = '--writable-tmpfs -B /scratch'
    // Passing arguments to CheckM2: --extension gz since the bins are gzipped fasta fules
    // --tmpdir to allow diamond to write intermediate files in scratch
    // --lowmem to stabilize the memory usage 
    ext.args = '--extension gz --tmpdir /scratch/checkm2_tmp --lowmem'
    beforeScript = '''
    set -e
    export TMPDIR=/scratch/\${USER}/\${SLURM_JOB_ID:-\$\$}
    mkdir -p \$TMPDIR
    mkdir -p /scratch/checkm2_tmp
    '''
    errorStrategy = 'retry'
    maxRetries    = 2
  }

  withName: /.*METATOR_PIPELINE/ {
    // -B /scratch metator also needs to write temporary files
    containerOptions = '-B /scratch'
    beforeScript = '''
    set -e
    export TMPDIR=/scratch/${USER}/metator_${SLURM_JOB_ID:-$$}
    mkdir -p \$TMPDIR
    '''
  }

  withName: /.*TRNASCAN_SE/ {
    cpus   = 6
    memory = '8 GB'
    time   = '4h'
    // tRNA-SCAN needed tmpdir to be specified in bash and apptainer 
    beforeScript = '''
    set -e
    export TMPDIR=\${PWD}/tmp
    export APPTAINERENV_TMPDIR=\${PWD}/tmp
    export SINGULARITYENV_TMPDIR=\${PWD}/tmp
    mkdir -p \${PWD}/tmp
    '''
    containerOptions = '--env TMPDIR=\${PWD}/tmp'
  }

}

singularity {
  pullTimeout = '1h'
  enabled     = true
  autoMounts  = true
  cacheDir    = "$RUN_DIR/.apptainer_cache"
}

report {
  enabled   = true
  overwrite = true
  file      = "$RUN_DIR/logs/report.html"
}

timeline {
  enabled   = true
  overwrite = true
  file      = "$RUN_DIR/logs/timeline.html"
}

trace {
  enabled   = true
  overwrite = true
  file      = "$RUN_DIR/logs/trace.txt"
}

EOF
# Run the pipeline on the fork 
# Move to running directpry to create all NetFlow files there 
cd "$RUN_DIR"
# nextflow run to execute the pipeline from a local for directory
nextflow run ${PIPELINE_DIR} \
  -profile singularity \
  -c "$RUN_DIR/nf/nextflow.config" \
  --input "$RUN_DIR/nf/input.yaml" \
  --outdir "$OUTDIR" \
  \
  --long_read_mapping_reads_per_chunk 1000000 \
  --hic_aligner bwamem2 \
  --hic_mapping_cram_bin_size 10000 \
  --hic_mapping_minq 10 \
  \
  --enable_metamdbg \
  --enable_rrna_prediction \
  --enable_genomad \
  --genomad_db "$GENOMAD_DB" \
  \
  --enable_binning \
  --enable_metabat2 \
  --enable_maxbin2 \
  --enable_bin3c \
  --enable_metator \
  --enable_semibin2 \
  --enable_lorbin \
  --minimum_contig_size 3000 \
  --minimum_hifi_perc_identity 99 \
  --minimum_bin_size 150000 \
  \
  --enable_bin_refinement \
  --enable_dastool \
  --enable_magscot \
  \
  --enable_binqc \
  --enable_checkm2 \
  --checkm2_db "$CHECKM2_DB" \
  --enable_trnascan_se \
  \
  --enable_taxonomy \
  --enable_gtdbtk \
  --gtdbtk_db "$GTDB_DB" \
  --hmm_gtdb_tigrfam "$HMM_TIGRFAM" \
  --hmm_gtdb_pfam "$HMM_PFAM" \
  --gtdb_bac120_metadata "$GTDB_BAC120_META" \
  --gtdb_ar53_metadata "$GTDB_AR53_META" \
  --gtdbtk_use_full_tree \
  --gtdbtk_min_completeness 50 \
  --gtdbtk_max_contamination 10 \
  -resume