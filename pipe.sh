#!/bin/bash

# Load configuration variables from config.sh
source /mnt/data/lj752/tools/Nanopore_amplicon_pipeline/config.sh

# Load sub-configuration for specific pipeline parameters
source /mnt/data/lj752/tools/Nanopore_amplicon_pipeline/bin/subconfig/NoCHIMERA-MM-50_ST-1_Q-0.7_KKc-0.05_KKhg-4.sh

# Save the current directory path to a variable
current_dir=$(pwd)

# Generate a timestamp for the run
run_time=$(date +"%d-%m-%Y_%H-%M-%S")

# Set raw_data and id from command-line arguments
raw_data="$1"  # The first argument is the input file path
id="$2"        # The second argument is the sample ID
log="$3"       # The third argument is the log file name

# Append a separator to the log file
echo -e "///////////////////////////////////////////////////////////////////////////////////////////" >> "${log}"

# Create a directory for the sample ID
mkdir -p ./${id}

# Change to the new directory
cd ./${id}

# Total number of tasks, adjust this based on your actual tasks
total_tasks=12

# Define the progress file path
progress_file="${current_dir}/${id}/progress.txt"

# Count the number of reads in the raw data file
input_count=$(grep -c '^@' "$raw_data")

# Initialize progress file
echo "0,...and so it begins: ${input_count} raw reads" > "${progress_file}"

# Log the start message
echo "0,...and so it begins: ${input_count} raw reads" >> "${log}"

# Start the Python progress monitor in the background and pass the progress file path to it
python "${progress_monitor}" "${progress_file}" "${total_tasks}" "${id}" "${log}" &

# Save the PID of the background Python job to ensure we can wait for this specific job later
PYTHON_PID=$!

# Ensure exit cleans up after itself
trap 'kill ${PYTHON_PID} 2>/dev/null; rm -f ${progress_file}; mv ${log} ${current_dir}/${id}/logs/; exit' EXIT INT TERM

# Setup hardware configuration based on the variable hardware_use
if [ "${hardware_use}" = "light" ]; then
  source "${subconfig}/hardware-light.sh"
elif [ "${hardware_use}" = "heavy" ]; then
  source "${subconfig}/hardware-heavy.sh"
elif [ "${hardware_use}" = "super-light" ]; then
  source "${subconfig}/hardware-super-light.sh"
fi

# Setup GPU configuration if greedy_gpu is true
if [ "${greedy_gpu}" = "true" ]; then
  source "${subconfig}/GPU-greedy.sh"
fi

# Setup amplicon configuration based on amplicon_pre_set
if [ "${amplicon_pre_set}" = "515y-926r" ]; then
  source "${subconfig}/AMP_515y-926r.sh"
fi

# Create directories for different stages of the pipeline
mkdir -p ./prep/bin
mkdir -p ./prep/filter
mkdir -p ./prep/polishing
mkdir -p ./prep/RAW
mkdir -p ./PROK/fasta
mkdir -p ./EUK/fasta
mkdir -p ./merge/bin
mkdir -p ./logs

##### Filtration and trimming by quality ######
echo "1,Filtering and trimming ${input_count} reads to Phred ${phred} and < ${filtered_amplicon_length}" > "${progress_file}"
echo "1,Filtering and trimming ${input_count} reads to Phred ${phred} and < ${filtered_amplicon_length}" >> "${log}"

if [ ! -s "${raw_data}" ]; then
    echo "### Input file ${raw_data} is empty or does not exist ###" >> "${log}"
    exit 1
fi

fastp \
--in1 "${raw_data}" \
--out1 "./prep/filter/${id}_${phred}_Phred.fastq" \
--cut_front \
--cut_tail \
--cut_mean_quality "${phred_t}" \
--length_required 200 \
--average_qual "${phred}" \
--json "./prep/bin/${id}_fastp_report.json" \
--html "./prep/bin/${id}_fastp_report.html" \
>> "${log}" 2>&1 || { echo "### fastp failed ###" >> "${log}"; exit 1; }

if [ ! -s "./prep/filter/${id}_${phred}_Phred.fastq" ]; then
    echo "### Filtered FASTQ file is empty ###" >> "${log}"
    exit 1
fi

seq_count_raw=$(grep -c '^@' "./prep/filter/${id}_${phred}_Phred.fastq")
percentage_filt_retained=$(echo "scale=2; ${seq_count_raw} / ${input_count} * 100" | bc)

if [ "${seq_count_raw}" -eq 0 ]; then
    echo "### No reads retained after filtering ###" >> "${log}"
    exit 1
fi

##### Chimera removal ######
echo "2,Searching for Chimeras in ${seq_count_raw} reads (${percentage_filt_retained}% retained from filtration)" > "${progress_file}"
echo "2,Searching for Chimeras in ${seq_count_raw} reads (${percentage_filt_retained}% retained from filtration)" >> "${log}"

# Convert the filtered FASTQ file to FASTA
vsearch --fastq_filter "./prep/filter/${id}_${phred}_Phred.fastq" --fastaout "./prep/${id}_filtered.fasta" --fastq_qmax 50 >> "${log}" 2>&1 || { echo "### vsearch fastq_filter failed ###" >> "${log}"; exit 1; }

# Check if the filtered FASTA file is created
if [ ! -s "./prep/${id}_filtered.fasta" ]; then
    echo "### vsearch fastq_filter failed: filtered file is empty ###" >> "${log}"
    exit 1
fi

# Run vsearch to detect and remove chimeras
vsearch --uchime_ref "./prep/${id}_filtered.fasta" --db "${fasta_16s_database}" --nonchimeras "./prep/${id}_nonchimeric.fasta" --chimeras "./prep/${id}_chimeras.fasta" --threads 20 >> "${log}" 2>&1

# Check if non-chimeric file is created
if [ ! -s "./prep/${id}_nonchimeric.fasta" ]; then
    echo "### vsearch chimera detection failed: non-chimeric file is empty ###" >> "${log}"
    exit 1
fi

# Deduplicate sequences using vsearch
vsearch --derep_fulllength "./prep/${id}_nonchimeric.fasta" --output "./prep/${id}_nonchimeric_dedup.fasta" --sizeout >> "${log}" 2>&1 || { echo "### vsearch deduplication failed ###" >> "${log}"; exit 1; }

# QC
seq_count_nonchimeric=$(grep -c '^>' "./prep/${id}_nonchimeric_dedup.fasta")
percentage_nonchimeric_retained=$(echo "scale=2; ${seq_count_nonchimeric} / ${seq_count_raw} * 100" | bc)

echo "3,Isolating ${seq_count_nonchimeric} (partially dereplicated count) non-chimeric amplicons" > "${progress_file}"
echo "3,Isolating ${seq_count_nonchimeric} (partially dereplicated count) non-chimeric amplicons" >> "${log}"

# Mapping with minimap2
minimap2 -ax map-ont -t 18 "./prep/${id}_nonchimeric_dedup.fasta" "./prep/filter/${id}_${phred}_Phred.fastq" > "./prep/${id}_mapped.sam"
samtools view -Sb "./prep/${id}_mapped.sam" > "./prep/${id}_mapped.bam"
samtools sort "./prep/${id}_mapped.bam" -o "./prep/${id}_sorted.bam"
samtools index "./prep/${id}_sorted.bam"

# Check if the mapping step succeeded
if [ ! -s "./prep/${id}_mapped.sam" ]; then
    echo "### minimap2 mapping failed: SAM file is empty ###" >> "${log}"
    exit 1
fi
if [ ! -s "./prep/${id}_sorted.bam" ]; then
    echo "### samtools sorting failed: BAM file is empty ###" >> "${log}"
    exit 1
fi

# Extract mapped reads
samtools view -@ "${cores}" -b -F 4 "./prep/${id}_sorted.bam" | samtools fasta -@ "${cores}" > "./prep/${id}_filtered_mapped.fasta" 2>> "${log}"

seq_count_filtered=$(grep -c '^>' "./prep/${id}_filtered_mapped.fasta")
percentage_retained=$(echo "scale=5; ${seq_count_filtered} / ${input_count} * 100" | bc)

# Log counts after filtering
echo "Filtered read count after mapping: ${seq_count_filtered}" >> "${log}"

if [ "${seq_count_filtered}" -eq 0 ]; then
  echo "No reads retained after filtering. Exiting." >> "${log}"
  exit 1
fi

###### Binning and error correction ######
echo "4,Binning ${seq_count_filtered} reads for error correcting using ${default_database_name} ${default_database_version} (${percentage_retained}% retained from QC)" > "${progress_file}"
echo "4,Binning ${seq_count_filtered} reads for error correcting using ${default_database_name} ${default_database_version} (${percentage_retained}% retained from QC)" >> "${log}"

## Binning
# 16S binning
minimap2 -ax map-ont -t "${cores}" --secondary=no -s "${minimap_alighment_quality_filter}" -r "${gene_coverage_16s}" "${fasta_16s_database}" "./prep/${id}_filtered_mapped.fasta" > "./prep/${id}_16s_mapped.sam" 2>> "${log}"
samtools view -@ "${cores}" -bS "./prep/${id}_16s_mapped.sam" > "./prep/${id}_16s_mapped.bam" 2>> "${log}"
samtools sort -@ "${cores}" -m "${memory_per_core}G" "./prep/${id}_16s_mapped.bam" -o "./prep/${id}_16s_mapped_sorted.bam" 2>> "${log}"
samtools index "./prep/${id}_16s_mapped_sorted.bam" 2>> "${log}"
samtools view -@ "${cores}" -b -F 4 -F 2048 -q "${samtools_quality_filter}" "./prep/${id}_16s_mapped_sorted.bam" | samtools fasta -@ "${cores}" > "./prep/${id}_16s_bin.fasta" 2>> "${log}"

# 18S binning
minimap2 -ax map-ont -t "${cores}" --secondary=no -s "${minimap_alighment_quality_filter}" -r "${gene_coverage_18s}" "${fasta_18s_database}" "./prep/${id}_filtered_mapped.fasta" > "./prep/${id}_18s_mapped.sam" 2>> "${log}"
samtools view -@ "${cores}" -bS "./prep/${id}_18s_mapped.sam" > "./prep/${id}_18s_mapped.bam" 2>> "${log}"
samtools sort -@ "${cores}" -m "${memory_per_core}G" "./prep/${id}_18s_mapped.bam" -o "./prep/${id}_18s_mapped_sorted.bam" 2>> "${log}"
samtools index "./prep/${id}_18s_mapped_sorted.bam" 2>> "${log}"
samtools view -@ "${cores}" -b -F 4 -F 2048 -q "${samtools_quality_filter}" "./prep/${id}_18s_mapped_sorted.bam" | samtools fasta -@ "${cores}" > "./prep/${id}_18s_bin.fasta" 2>> "${log}"

# QC
binned_18s=$(grep -c '^>' "./prep/${id}_18s_bin.fasta")
binned_16s=$(grep -c '^>' "./prep/${id}_16s_bin.fasta")
per_18s_bin=$(echo "scale=5; (($binned_18s / ${seq_count_filtered}) * 100)" | bc)
per_16s_bin=$(echo "scale=5; (($binned_16s / ${seq_count_filtered}) * 100)" | bc)
binned_read_representation=$(echo "scale=5; (($per_16s_bin + $per_18s_bin) - 100)" | bc)
read_number_discrepancy=$(echo "scale=5; (($binned_16s + $binned_18s) - $seq_count_filtered)" | bc)

echo "###################################
16s/18s binning completion stats:
  16s (${per_16s_bin}%) ${binned_16s} reads
  18s (${per_18s_bin}%) ${binned_18s} reads
discrepancy = ${binned_read_representation}% (${read_number_discrepancy} reads)
       --used ${default_database_name} ${default_database_version}--" >> "${log}"

## Medaka iteration 1
# Log and update progress for Medaka iteration 1
echo "5,Medaka iteration 1 - using ${binned_16s} (${per_16s_bin}%) 16s and ${binned_18s} (${per_18s_bin}%) 18s reads" > "${progress_file}"
echo "5,Medaka iteration 1 - using ${binned_16s} (${per_16s_bin}%) 16s and ${binned_18s} (${per_18s_bin}%) 18s reads" >> "${log}"

# Activate Medaka environment
eval "$(conda shell.bash hook)"
conda activate "${medaka}" || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

# Run Medaka consensus for 16s
echo "Iteration 1 - 16s" >> "./prep/medaka_stdout_log.txt"
medaka_consensus -i "./prep/${id}_16s_bin.fasta" -d "${fasta_16s_database}" -o "./prep/${id}_16s_medaka_1" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g -x >> "./prep/medaka_stdout_log.txt" 2>&1

# Run Medaka consensus for 18s
echo "Iteration 1 - 18s" >> "./prep/medaka_stdout_log.txt"
medaka_consensus -i "./prep/${id}_18s_bin.fasta" -d "${fasta_18s_database}" -o "./prep/${id}_18s_medaka_1" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g -x >> "./prep/medaka_stdout_log.txt" 2>&1

# QC
seq_count_medaka_16s=$(grep -c '^>' "./prep/${id}_16s_medaka_1/consensus.fasta")
seq_count_medaka_18s=$(grep -c '^>' "./prep/${id}_18s_medaka_1/consensus.fasta")
percentage_retained_16s_medaka=$(echo "scale=5; ${seq_count_medaka_16s} / ${binned_16s} * 100" | bc)
percentage_retained_18s_medaka=$(echo "scale=5; ${seq_count_medaka_18s} / ${binned_18s} * 100" | bc)
conda deactivate

# Copy consensus files for further processing
cp "./prep/${id}_16s_medaka_1/consensus.fasta" "./prep/${id}_16s.fasta"
cp "./prep/${id}_18s_medaka_1/consensus.fasta" "./prep/${id}_18s.fasta"

# Count sequences in the final consensus files
seq_count_medaka_16s_2=$(grep -c '^>' "./prep/${id}_16s.fasta")
seq_count_medaka_18s_2=$(grep -c '^>' "./prep/${id}_18s.fasta")
total_amplicons_CON=$(echo "scale=2; (${seq_count_medaka_16s_2} + ${seq_count_medaka_18s_2})" | bc)

# Log final stats
echo "###################################
16s/18s error correction stats:
  16s ${binned_16s} reads reduced to ${seq_count_medaka_16s_2} amplicons (${percentage_retained_16s_medaka}%)
  18s ${binned_18s} reads reduced to ${seq_count_medaka_18s_2} amplicons (${percentage_retained_18s_medaka}%)
  totaling ${total_amplicons_CON} amplicons from ${seq_count_filtered} reads" >> "${log}"

###### Classification and abundance extraction ######
echo "7,Importing and classifying ${seq_count_medaka_16s_2} polished 16S amplicons in QIIME2 (${default_database_name} ${default_database_version})" > "${progress_file}"
echo "7,Importing and classifying ${seq_count_medaka_16s_2} polished 16S amplicons in QIIME2 (${default_database_name} ${default_database_version})" >> "${log}"
eval "$(conda shell.bash hook)"
conda activate "${biopy}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

# Import 16S reads to QIIME2
seqtk seq -F I "./prep/${id}_16s.fasta" > "./prep/${id}_16s.fastq" 2>> "${log}"
seqtk seq -F I "./prep/${id}_18s.fasta" > "./prep/${id}_18s.fastq" 2>> "${log}"
conda deactivate
eval "$(conda shell.bash hook)"
conda activate "${qiime}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

# Import and classify 16S amplicons
echo -e "sample-id\tabsolute-filepath" > "./prep/${id}_manifest_16s.tsv"
echo -e "${id}\t${current_dir}/${id}/prep/${id}_16s.fastq" >> "./prep/${id}_manifest_16s.tsv"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "./prep/${id}_manifest_16s.tsv" \
  --output-path "./prep/${id}_CON_qiime2_input_16s.qza" \
  --input-format SingleEndFastqManifestPhred33V2 \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: RAW data importing to qiime###" >> "${log}"; exit 1; }

qiime vsearch dereplicate-sequences \
  --i-sequences "./prep/${id}_CON_qiime2_input_16s.qza" \
  --o-dereplicated-table "./prep/${id}_CON_table-derep_16s.qza" \
  --o-dereplicated-sequences "./prep/${id}_CON_rep-seqs-derep_16s.qza" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: RAW data dereplication###" >> "${log}"; exit 1; }

qiime feature-classifier classify-sklearn \
  --i-classifier "${default_16s_classifier}" \
  --i-reads "./prep/${id}_CON_rep-seqs-derep_16s.qza" \
  --o-classification "./prep/${id}_taxonomy_16s.qza" \
  --p-confidence "${QIIME_confidence}" \
  --p-n-jobs "${Q_cores}" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: 16s taxonomic classification QIIME2###" >> "${log}"; exit 1; }

qiime tools export --input-path "./prep/${id}_taxonomy_16s.qza" --output-path "./prep/${id}_taxonomy_16s" >> "${log}" 2>&1
qiime tools export --input-path "./prep/${id}_CON_rep-seqs-derep_16s.qza" --output-path "./prep/${id}_CON_rep-seqs-derep_16s" >> "${log}" 2>&1

# Save taxonomy files as TSV
cp ./prep/${id}_taxonomy_16s/taxonomy.tsv ./prep/${id}_taxonomy_16s.tsv
cp ./prep/${id}_CON_rep-seqs-derep_16s/dna-sequences.fasta ./prep/${id}_CON_rep-seqs-derep_16s.fasta

percent_lost_16s=$(python "${tax_to_fasta}" "./prep/${id}_CON_rep-seqs-derep_16s/dna-sequences.fasta" "./prep/${id}_taxonomy_16s/taxonomy.tsv" "${fasta_16s_database}" "./prep/${id}_CON_taxonomy_rep_seq_16s.fasta" "7" "${log}" "${cores}")
sleep 5

# Import and classify 18S amplicons
echo "8,Importing and classifying ${seq_count_medaka_18s_2} polished 18S amplicons in QIIME2 (${default_database_name} ${default_database_version})" > "${progress_file}"
echo "8,Importing and classifying ${seq_count_medaka_18s_2} polished 18S amplicons in QIIME2 (${default_database_name} ${default_database_version})" >> "${log}"

echo -e "sample-id\tabsolute-filepath" > "./prep/${id}_manifest_18s.tsv"
echo -e "${id}\t${current_dir}/${id}/prep/${id}_18s.fastq" >> "./prep/${id}_manifest_18s.tsv"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "./prep/${id}_manifest_18s.tsv" \
  --output-path "./prep/${id}_CON_qiime2_input_18s.qza" \
  --input-format SingleEndFastqManifestPhred33V2 \
 >> "${log}" 2>&1 || { sleep 5; echo "###failed here: 18s CON data importing to qiime###" >> "${log}"; exit 1; }

qiime vsearch dereplicate-sequences \
  --i-sequences "./prep/${id}_CON_qiime2_input_18s.qza" \
  --o-dereplicated-table "./prep/${id}_CON_table-derep_18s.qza" \
  --o-dereplicated-sequences "./prep/${id}_CON_rep-seqs-derep_18s.qza" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: 18s CON data dereplication###" >> "${log}"; exit 1; }

qiime feature-classifier classify-sklearn \
  --i-classifier "${default_18s_classifier}" \
  --i-reads "./prep/${id}_CON_rep-seqs-derep_18s.qza" \
  --o-classification "./prep/${id}_taxonomy_18s.qza" \
  --p-confidence "${QIIME_confidence}" \
  --p-n-jobs "${Q_cores}" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: 18s taxonomic classification QIIME2###" >> "${log}"; exit 1; }

qiime tools export --input-path "./prep/${id}_taxonomy_18s.qza" --output-path "./prep/${id}_taxonomy_18s" >> "${log}" 2>&1
qiime tools export --input-path "./prep/${id}_CON_rep-seqs-derep_18s.qza" --output-path "./prep/${id}_CON_rep-seqs-derep_18s" >> "${log}" 2>&1

# Save taxonomy files as TSV
cp ./prep/${id}_taxonomy_18s/taxonomy.tsv ./prep/${id}_taxonomy_18s.tsv
cp ./prep/${id}_CON_rep-seqs-derep_18s/dna-sequences.fasta ./prep/${id}_CON_rep-seqs-derep_18s.fasta

percent_lost_18s=$(python "${tax_to_fasta}" "./prep/${id}_CON_rep-seqs-derep_18s/dna-sequences.fasta" "./prep/${id}_taxonomy_18s/taxonomy.tsv" "${fasta_18s_database}" "./prep/${id}_CON_taxonomy_rep_seq_18s.fasta" "11" "${log}" "${cores}")
conda deactivate

### Kraken2 and Bracken
# Log and update progress for Kraken2 and Bracken analysis
echo "9,RAW 16S reads being assigned to CON 16s (${percent_lost_16s}% lost; 0.00% ideal) reads in Kraken2" > "${progress_file}"
echo "9,RAW 16S reads being assigned to CON 16s (${percent_lost_16s}% lost; 0.00% ideal) reads in Kraken2" >> "${log}"

# Activate Biopython environment
eval "$(conda shell.bash hook)"
conda activate "${biopy}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

# Ensure the output directory exists
mkdir -p "./prep/${id}/"

## 16S
# Generate CON database
stdbuf -oL kraken2-build --add-to-library "./prep/${id}_CON_taxonomy_rep_seq_16s.fasta" --db "./prep/${id}_16s_${run_time}-CON" >> "${log}" 2>&1
stdbuf -oL kraken2-build --db "./prep/${id}_16s_${run_time}-CON" --download-taxonomy >> "${log}" 2>&1
stdbuf -oL kraken2-build --db "./prep/${id}_16s_${run_time}-CON" --kmer-len 15 --minimizer-len 10 --minimizer-spaces 2 --threads "${cores}" --build >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Kraken 16s database creation###" >> "${log}"; exit 1; }

# Verify the paths for Kraken2 outputs
echo "Verifying paths for Kraken2 outputs..." >> "${log}"
echo "Path for classification output: ./prep/${id}/kraken_classification_16s.txt" >> "${log}"
echo "Path for report output: ./prep/${id}/kraken_report_16s.txt" >> "${log}"

# Classify 16S RAW Data
stdbuf -oL kraken2 --db "./prep/${id}_16s_${run_time}-CON" \
  --output "./prep/${id}/kraken_classification_16s.txt" \
  --report "./prep/${id}/kraken_report_16s.txt" \
  --threads "${cores}" \
  --confidence "${KrKn_confidence}" \
  --minimum-hit-groups "${KrKn_hg}" \
  "./prep/${id}_16s_bin.fasta" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Kraken matching raw and 16s con data ###" >> "${log}"; exit 1; }

# Check if the Kraken2 report file for 16S exists
if [ ! -f "./prep/${id}/kraken_report_16s.txt" ]; then
    echo "###failed here: Kraken report 16s file not found ###" >> "${log}"
    exit 1
fi

# Recalculate 16s abundances
bracken-build -l "${len_16s}" -d "./prep/${id}_16s_${run_time}-CON" -t "${cores}" >> "${log}" 2>&1
bracken -d "./prep/${id}_16s_${run_time}-CON" -i "./prep/${id}/kraken_report_16s.txt" -o "./prep/${id}/braken_16s.tsv" -r "${len_16s}" -l S -t "${cores}" >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Bracken matching raw and 16s con data ###" >> "${log}"; exit 1; }

# Log and update progress for Kraken2 analysis for 18S
echo "10,RAW 18S reads being assigned to CON 18s (${percent_lost_18s}% lost; 0.00% ideal) reads in Kraken2" > "${progress_file}"
echo "10,RAW 18S reads being assigned to CON 18s (${percent_lost_18s}% lost; 0.00% ideal) reads in Kraken2" >> "${log}"

## 18S
# Generate CON database
stdbuf -oL kraken2-build --add-to-library "./prep/${id}_CON_taxonomy_rep_seq_18s.fasta" --db "./prep/${id}_18s_${run_time}-CON" >> "${log}" 2>&1
stdbuf -oL kraken2-build --db "./prep/${id}_18s_${run_time}-CON" --download-taxonomy >> "${log}" 2>&1
stdbuf -oL kraken2-build --db "./prep/${id}_18s_${run_time}-CON" --kmer-len 15 --minimizer-len 10 --minimizer-spaces 2 --threads "${cores}" --build >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Kraken 18s database creation###" >> "${log}"; exit 1; }

# Verify the paths for Kraken2 outputs
echo "Verifying paths for Kraken2 outputs..." >> "${log}"
echo "Path for classification output: ./prep/${id}/kraken_classification_18s.txt" >> "${log}"
echo "Path for report output: ./prep/${id}/kraken_report_18s.txt" >> "${log}"

# Classify 18S RAW Data
stdbuf -oL kraken2 --db "./prep/${id}_18s_${run_time}-CON" \
  --output "./prep/${id}/kraken_classification_18s.txt" \
  --report "./prep/${id}/kraken_report_18s.txt" \
  --threads "${cores}" \
  --confidence "${KrKn_confidence}" \
  --minimum-hit-groups "${KrKn_hg}" \
  "./prep/${id}_18s_bin.fasta" \
  >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Kraken matching raw and 18s con data ###" >> "${log}"; exit 1; }

# Check if the Kraken2 report file for 18S exists
if [ ! -f "./prep/${id}/kraken_report_18s.txt" ]; then
    echo "###failed here: Kraken report 18s file not found ###" >> "${log}"
    exit 1
fi

# Recalculate 18s abundances
bracken-build -l "${len_18s}" -d "./prep/${id}_18s_${run_time}-CON" -t "${cores}" >> "${log}" 2>&1
bracken -d "./prep/${id}_18s_${run_time}-CON" -i "./prep/${id}/kraken_report_18s.txt" -o "./prep/${id}/braken_18s.tsv" -r "${len_18s}" -l S -t "${cores}" >> "${log}" 2>&1 || { sleep 5; echo "###failed here: Bracken matching raw and 18s con data ###" >> "${log}"; exit 1; }

# Deactivate conda environment
conda deactivate

### Normalisation, bias correction and merge
# Log and update progress for normalization, bias correction, and merging data
echo "11,Normalising, scaling and merging data" > "${progress_file}"
echo "11,Normalising, scaling and merging data" >> "${log}"

# Activate Python environment
eval "$(conda shell.bash hook)"
conda activate "${python}"

# Remove excess from bracken outputs and re-name headers Python doesnt like bracken headers
python "${trim_tsv}" "${current_dir}/${id}${output_p}/braken_16s.tsv" "${current_dir}/${id}${in_m}/braken_16s_trimmed.tsv" "${log}"
python "${trim_tsv}" "${current_dir}/${id}${output_e}/braken_18s.tsv" "${current_dir}/${id}${in_m}/braken_18s_trimmed.tsv" "${log}"

# Bias correction
python "${bias_correction}" "${current_dir}/${id}${in_m}/braken_16s_trimmed.tsv" "${bias_factor_16s}" "${current_dir}/${id}${in_m}/braken_16s_bias_corrected.tsv" "${log}"
python "${bias_correction}" "${current_dir}/${id}${in_m}/braken_18s_trimmed.tsv" "${bias_factor_18s}" "${current_dir}/${id}${in_m}/braken_18s_bias_corrected.tsv" "${log}"

# Count total reads and normalize
python "${normalise}" "${current_dir}/${id}${in_m}/braken_16s_bias_corrected.tsv" "${norm_factor}" "${current_dir}/${id}${in_m}/16s_microbiome.tsv" "${log}"
python "${normalise}" "${current_dir}/${id}${in_m}/braken_18s_bias_corrected.tsv" "${norm_factor}" "${current_dir}/${id}${in_m}/18s_microbiome.tsv" "${log}"

# Merge microbiomes
python "${merge_16s_18s}" "${current_dir}/${id}${in_m}/18s_microbiome.tsv" "${current_dir}/${id}${in_m}/16s_microbiome.tsv" "${current_dir}/${id}${out_m}/${run_time}-${id}_microbiome.tsv" "${log}"

# Deactivate conda environment
conda deactivate

# Log completion message
echo "12,Done: see ${current_dir}/${id}/${run_time}-${id}_microbiome.tsv" > "${progress_file}"
echo "12,Done: see ${current_dir}/${id}/${run_time}-${id}_microbiome.tsv" >> "${log}"

# Wait for the Python script to finish before removing the progress file and exiting
sleep 2
wait "${PYTHON_PID}"
rm -f "${progress_file}"
