#!/bin/bash
# Welcome to the config!

###     User? Please edit me:   ###

## Essential
# General use
greedy_gpu="true" # Would you like your GPU to be fully utalised?
hardware_use="heavy" # 'super-light' uses 20%, 'light' uses 45%, and 'heavy' will use 90% of your computational resourses
amplicon_pre_set="515y-926r" # premade configs to change amplicon handeling
# Hardware
export all_cores="20" # All threads availible
export all_RAM="250" # All RAM (memory)
# Conda envs
export qiime="qiime2_amplicon_env"
export biopy="biopy38_env"
export medaka="medaka_env"
export python="python3.10"
# Where are we?
export w_d="/mnt/data/lj752/tools/Nanopore_amplicon_pipeline" # Directory where tool is stored

# Dorado personalisation (so you can copy input from 'ns dorado' prompts)
export defualt_kit="SQK-NBD114-24" # Dorado
export defualt_model="dna_r10.4.1_e8.2_400bps_sup@v4.3.0" # Dorado
export medaka_model="r1041_e82_400bps_sup_v4.3.0"
# Filtering by quality and chimeras
export phred="20" # Average Phred desired
export phred_t="20" # Trimming of ends (recommended slighlty lower than phred above)
export norm_factor="100000"








###     User? Dont edit me please!      ###

# Defualt tool information
export pipeline_name="Nanopore_amplicon_pipeline"
export pipeline_database_version="0.3 (03/24)"
export default_database_version="138.1"
export default_database_name='SILVA'
# Locations
#export default_database="${w_d}/bin/databases/SILVA_138.1_NR99.fasta"
#export default_classifier="${w_d}/bin/databases/SILVA_138.1_NR99.qza"
export default_16s_classifier="${w_d}/bin/databases/16s_SILVA_138.1_NR99.qza"
export default_18s_classifier="${w_d}/bin/databases/18s_SILVA_138.1_NR99.qza"
export fasta_18s_database="${w_d}/bin/databases/18s_SILVA_138.1_NR99.fasta"
export fasta_16s_database="${w_d}/bin/databases/16s_SILVA_138.1_NR99.fasta"
#export NCBI="/mnt/data/lj752/tools/databases/NCBI_nt_03-24/nt"
# Python
export progress_monitor="${w_d}/bin/scripts/progress_monitor.py"
#export blastn_con="${w_d}/bin/scripts/blastn_consensus.py"
#export abundance_sum="${w_d}/bin/scripts/abundance_sum.py"
#export merge_taxonomy_abundance="${w_d}/bin/scripts/merge_taxonomy_abundance.py"
#export trim_otu="${w_d}/bin/scripts/trim_otu.py"
#export qc_microbiome="${w_d}/bin/scripts/qc_microbiome.py"
export tax_to_fasta="${w_d}/bin/scripts/merge_taxonomy_to_fasta.py"
#export calculate_loss_percentage="${w_d}/bin/scripts/calculate_loss_percentage.py"
export normalise="${w_d}/bin/scripts/normalise.py"
export merge_16s_18s="${w_d}/bin/scripts/merge_16s_18s.py"
export trim_tsv="${w_d}/bin/scripts/trim_tsv.py"
