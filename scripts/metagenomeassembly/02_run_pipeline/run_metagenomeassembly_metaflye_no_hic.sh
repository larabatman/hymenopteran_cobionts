#!/usr/bin/env bash
#SBATCH --partition=pibu_el8
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --job-name=metaasm_metaflye
#SBATCH --output=logs/metaasm_metaflye_%j.out
#SBATCH --error=logs/metaasm_metaflye_%j.err

# Pilot job for the NextFlow sanger-tol/metagenomeassembly pipeline 
# SLURM directives above define the allocation that runs NextFlow, and then it submits its own SLORM hobs with specific resources from the NextFlow config file defined below
# This is the variant to enter the pipeline with an existing metaFlye assembly, so no metaMDBG
# Moreover, this is the ariant for which there is no HiC data, so HiC parameters are disabled
# Usage:
# sbatch run_metagenomeassembly_metaflye.sh species

set -euo pipefail


SPECIES="${1}"

RUN_TAG="attempt_hic_metaflye"
PROJECT_ROOT="/data/projects/p2025-0083_mining_cobionts"
PIPELINE_DIR="${PROJECT_ROOT}/metagenomeassembly"
RUN_DIR="${PROJECT_ROOT}/runs/${SPECIES}/${RUN_TAG}"
OUTDIR="${PROJECT_ROOT}/results/${SPECIES}/metaflye_binners"
# Input files: CLR cleaned fasta
CLR_FASTA="${PROJECT_ROOT}/reads/processed/${SPECIES}/clr_filtered.fasta.gz"
# Input assembly: metaflye
METAFLYE_FASTA="${PROJECT_ROOT}/assemblies/${SPECIES}/metaflye/${SPECIES}_metaflye.fasta.gz"

GTDB_DB="${PROJECT_ROOT}/databases/gtdb_r226/release226"
HMM_TIGRFAM="${GTDB_DB}/markers/tigrfam/tigrfam.hmm"
HMM_PFAM="${GTDB_DB}/markers/pfam/Pfam-A.hmm"
GTDB_BAC120_META="${PROJECT_ROOT}/databases/gtdb_r226/bac120_metadata_r226.tsv"
GTDB_AR53_META="${PROJECT_ROOT}/databases/gtdb_r226/ar53_metadata_r226.tsv"
CHECKM2_DB="${PROJECT_ROOT}/databases/checkm2_102/CheckM2_database/uniref100.KO.1.dmnd"
GENOMAD_DB="${PROJECT_ROOT}/databases/genomad_db"

export TMPDIR="/scratch/${USER}/tmp_${SLURM_JOB_ID}"
export APPTAINER_TMPDIR="/scratch/${USER}/apptainer_${SLURM_JOB_ID}"

module purge
module load Java/17.0.6

mkdir -p logs "${RUN_DIR}"/{nf,logs,.nextflow,.apptainer_cache,work} "${OUTDIR}" "${TMPDIR}" "${APPTAINER_TMPDIR}"

export NXF_HOME="${RUN_DIR}/.nextflow"
export NXF_WORK="${RUN_DIR}/work"
export NXF_CACHE_DIR="${RUN_DIR}/.nextflow/cache"
export NXF_OPTS='-Xms4g -Xmx16g'
export APPTAINER_CACHEDIR="${RUN_DIR}/.apptainer_cache"
export NXF_SINGULARITY_CACHEDIR="${RUN_DIR}/.apptainer_cache"

# Write input YAML
# assembly: assembler metamdbg tag, but actually metaflye and assembly: fasta points to the metaflye assembly
cat > "${RUN_DIR}/nf/input.yaml" <<EOF
id: ${SAMPLE_ID}

assembly:
  fasta: ${METAFLYE_FASTA}
  assembler: metamdbg

pacbio:
  fasta:
    - ${CLR_FASTA}
EOF

# Write NextFlow config file
cat > "${RUN_DIR}/nf/nextflow.config" <<EOF

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
    ext.args = '--threads 8'
  }

  withName: /.*PYRODIGAL/ {
  cpus          = 8
  memory        = { 48.GB * task.attempt }   // 48 → 96 → 144 GB
  time          = '4h'
  errorStrategy = 'retry'
  maxRetries    = 2
  }

  withName: /.*FASTXALIGN_MINIMAP2ALIGN/ {
    cpus   = 8
    memory = '60 GB'
    scratch = true
  }

    withName: /.*MAXBIN2/ {
    errorStrategy = 'ignore'
  }

  withName: /.*MINIMAP2_INDEX/ {
    cpus          = 4
    memory        = { 80.GB * task.attempt }
    time          = '4h'
    errorStrategy = 'retry'
    maxRetries    = 3
  }

  withName: /.*INFERNAL_CMSEARCH/ {
    cpus          = 6
    memory        = { 40.GB * task.attempt }
    time          = '8h'
    errorStrategy = 'retry'
    maxRetries    = 3
  }

  withName: /.*CHECKM2_PREDICT/ {
    cpus   = 4
    memory = '96 GB'
    time   = '24h'
    containerOptions = '--writable-tmpfs -B /scratch'
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

# Run local nextflow pipeline:
# --enable_metamdbg is not passed, and --minimum_hifi_perc_identity is lowered to 85
# disabled hic parameters
cd "${RUN_DIR}"
nextflow run "${PIPELINE_DIR}" \
  -profile singularity \
  -c "${RUN_DIR}/nf/nextflow.config" \
  --input "${RUN_DIR}/nf/input.yaml" \
  --outdir "${OUTDIR}" \
  \
  --long_read_mapping_reads_per_chunk 1000000 \
  \
  --enable_rrna_prediction \
  --enable_genomad \
  --genomad_db "${GENOMAD_DB}" \
  \
  --enable_binning \
  --enable_metabat2 \
  --enable_maxbin2 \
  --enable_semibin2 \
  --enable_lorbin \
  --minimum_contig_size 2000 \
  --minimum_hifi_perc_identity 85 \
  --minimum_bin_size 150000 \
  \
  --enable_bin_refinement \
  --enable_dastool \
  --enable_magscot \
  \
  --enable_binqc \
  --enable_checkm2 \
  --checkm2_db "${CHECKM2_DB}" \
  --enable_trnascan_se \
  \
  --enable_taxonomy \
  --enable_gtdbtk \
  --gtdbtk_db "${GTDB_DB}" \
  --hmm_gtdb_tigrfam "${HMM_TIGRFAM}" \
  --hmm_gtdb_pfam "${HMM_PFAM}" \
  --gtdb_bac120_metadata "${GTDB_BAC120_META}" \
  --gtdb_ar53_metadata "${GTDB_AR53_META}" \
  --gtdbtk_use_full_tree \
  --gtdbtk_min_completeness 50 \
  --gtdbtk_max_contamination 10 \
  -resume