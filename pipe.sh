#!/bin/bash

# Resourses and variables
source config.sh
cd=$(pwd)
run_time=$(date +"%d-%m-%Y_%H-%M-%S")
# Set raw_data and id from command-line arguments
raw_data="$1"  # The first argument is the input file path
id="$2"        # The second argument is the sample ID
log="$3"	# The log name in
echo -e "///////////////////////////////////////////////////////////////////////////////////////////" >> "${log}"



###	Progress bar	###
#mkdir ./${id}
cd ./${id}


total_tasks=12  # Total number of tasks, adjust this based on your actual tasks
progress_file="${cd}/${id}/progress.txt" # Progress monitor
input_count=$(grep -c '^@' "$raw_data")

echo "0,...and so it begins: ${input_count} raw reads" > $progress_file
echo "0,...and so it begins: ${input_count} raw reads" >> ${log}
# Start the Python progress monitor in the background and pass the progress file path to it
python $progress_monitor "$progress_file" "$total_tasks" "$id" "${log}" &
# Save the PID of the background Python job to ensure we can wait for this specific job later
PYTHON_PID=$!
# Ensure exit cleans up after itself # Ensure cancleing the script doesnt leave redundant progress.txt script (also moved log)
trap 'kill $PYTHON_PID 2>/dev/null; rm -f $progress_file; mv $log $cd/$id/logs/; exit' EXIT INT TERM
# Setup
if [ "$hardware_use" = "light" ]; then
  source "$subconfig/hardware-light.sh"
elif [ "$hardware_use" = "heavy" ]; then
  source "$subconfig/hardware-heavy.sh"
elif [ "$hardware_use" = "super-light" ]; then
  source "$subconfig/hardware-super-light.sh"
fi
if [ "$greedy_gpu" = "true" ]; then
  source "$subconfig/GPU-greedy.sh"
fi

if [ "$amplicon_pre_set" = "515y-926r" ]; then
  source "$subconfig/AMP_515y-926r.sh"
fi

mkdir -p ./prep/bin
mkdir -p ./prep/filter
mkdir -p ./prep/polishing
mkdir -p ./prep/RAW
mkdir -p ./PROK/fasta
mkdir -p ./EUK/fasta
mkdir -p ./merge/bin
mkdir -p ./logs


######		Filtration and trimming by quality		######
#echo "1,Filtering and trimming ${input_count} reads to Phred ${phred} and < ${filtered_amplicon_length}" > $progress_file
#echo "1,Filtering and trimming ${input_count} reads to Phred ${phred} and < ${filtered_amplicon_length}" >> ${log}
#fastp \
#--in1 "${raw_data}" \
#--out1 ".${bb_f}/${id}_${phred}_Phred.fastq" \
#--cut_front \
#--cut_tail \
#--cut_mean_quality ${phred_t} \
#--length_required 200 \
#--average_qual ${phred} \
#--json ".${bb_out}/${id}_fastp_report.json" \
#--html ".${bb_out}/${id}_fastp_report.html" \
#>> ${log} 2>&1 || { echo "### fastp failed ###" >> "${log}"; exit 1; }

seq_count_raw=$(grep -c '^@' ".${bb_f}/${id}_${phred}_Phred.fastq")
percentage_filt_retained=$(echo "scale=2; $seq_count_raw / $input_count * 100" | bc)



#####		Chimera removal		######
#echo "2,Searching for Chimeras in ${seq_count_raw} reads (${percentage_filt_retained}% retained from filtration)" > $progress_file
#echo "2,Searching for Chimeras in ${seq_count_raw} reads (${percentage_filt_retained}% retained from filtration)" >> ${log}
#eval "$(conda shell.bash hook)"
#conda activate "${qiime}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

# Manifest
#echo -e "sample-id,absolute-filepath,direction" > manifest.csv
#echo -e "${id},${cd}/${id}${bb_f}/${id}_${phred}_Phred.fastq,forward" >> manifest.csv

# Import
#qiime tools import \
#  --type 'SampleData[SequencesWithQuality]' \
#  --input-path "manifest.csv" \
#  --output-path ".${bb_f}/${id}_${phred}_Phred.qza" \
#  --input-format SingleEndFastqManifestPhred33 \
#  >> ${log} 2>&1

# Dereplicate for computation simplicity
#qiime vsearch dereplicate-sequences \
#  --i-sequences ".${bb_f}/${id}_${phred}_Phred.qza" \
#  --o-dereplicated-table ".${bb_out}/table-derep.qza" \
#  --o-dereplicated-sequences ".${bb_out}/rep-seqs-derep.qza" \
#  >> ${log} 2>&1

# Chimera removal
#qiime vsearch uchime-ref \
#  --i-sequences ".${bb_out}/rep-seqs-derep.qza" \
#  --i-table ".${bb_out}/table-derep.qza" \
#  --i-reference-sequences "${default_classifier}" \
#  --o-chimeras ".${bb_out}/chimeras.qza" \
#  --o-nonchimeras ".${bb}/nonchimeras.qza" \
#  --o-stats ".${bb_out}/uchime-stats.qza" \
#  --quiet \
#  --p-threads ${cores} \
#  >> ${log} 2>&1

# Export non chimeric structures
#qiime tools export \
#  --input-path ".${bb}/nonchimeras.qza" \
#  --output-path ".${bb}/${id}_nonchimeric" \
#  >> ${log} 2>&1
#conda deactivate

# Filter raw data
#mv ".${bb}/${id}_nonchimeric/dna-sequences.fasta" ".${bb}/${id}_nonchimeric.fasta"
#rm -r ".${bb}/${id}_nonchimeric"
#eval "$(conda shell.bash hook)"
#conda activate "${biopy}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }
#seqtk seq -F I ".${bb}/${id}_nonchimeric.fasta" > ".${bb}/${id}_nonchimeric.fastq" 2>> "${log}"
#conda deactivate
#seq_count_nonchimeric=$(grep -c '^@' ".${bb}/${id}_nonchimeric.fastq")

#echo "3,Isolating ${seq_count_nonchimeric} (partially dereplicated count) non-chimeric amplicons" > $progress_file
#echo "3,Isolating ${seq_count_nonchimeric} (partially dereplicated count) non-chimeric amplicons" >> ${log}
# Minimap mapping non chimeric reads agianst filtered reads
#minimap2 -ax map-ont -t "${cores}" ".${bb}/${id}_nonchimeric.fastq" ".${bb_f}/${id}_${phred}_Phred.fastq" > ".${bb}/${id}_mapped.sam" 2>> "${log}"
# Extract non chimeric strcutures (This approach allows chimera detection before binning for optimal results, but without loss of abundance data through binning of depricated data)
#samtools view -@ "${cores}" -bS ".${bb}/${id}_mapped.sam" | samtools sort -@ "${cores}" -m "${memory_per_core}G" -o ".${bb}/${id}_mapped_sorted.bam"  2>> "${log}"
#samtools view -@ "${cores}" -b -F 4 ".${bb}/${id}_mapped_sorted.bam" | samtools fastq -@ "${cores}" > ".${bb}/${id}_filtered.fastq" 2>> "${log}"
# Count and present as fraction
seq_count_filtered=$(grep -c '^@' ".${bb}/${id}_filtered.fastq")
percentage_retained=$(echo "scale=5; $seq_count_filtered / $input_count * 100" | bc)



######		 Binning and error correction		######
echo "4,Binning ${seq_count_filtered} reads for error correcting using ${default_database_name} ${default_database_version} (${percentage_retained}% retained from QC)" > $progress_file
echo "4,Binning ${seq_count_filtered} reads for error correcting using ${default_database_name} ${default_database_version} (${percentage_retained}% retained from QC)" >> ${log}
## Binnning
# 16S
minimap2 -ax map-ont -t "${cores}" --secondary=no -s "$minimap_alighment_quality_filter" -r "$gene_coverage_16s" "${fasta_16s_database}" ".${bb}/${id}_filtered.fastq" > ".${bb_p}/${id}_16s_mapped.sam" 2>> "${log}"
samtools view -@ "${cores}" -bS ".${bb_p}/${id}_16s_mapped.sam" > ".${bb_p}/${id}_16s_mapped.bam" 2>> "${log}"
samtools sort -@ "${cores}" -m "${memory_per_core}G" ".${bb_p}/${id}_16s_mapped.bam" -o ".${bb_p}/${id}_16s_mapped_sorted.bam" 2>> "${log}"
samtools index ".${bb_p}/${id}_16s_mapped_sorted.bam" 2>> "${log}"
samtools view -@ "${cores}" -b -F 4 -F 2048 -q "$samtools_quality_filter" ".${bb_p}/${id}_16s_mapped_sorted.bam" | samtools fasta -@ "${cores}" > ".${output_s}/${id}_16s_bin.fasta" 2>> "${log}"
# 18S
minimap2 -ax map-ont -t "${cores}" --secondary=no -s "$minimap_alighment_quality_filter" -r "$gene_coverage_18s" "${fasta_18s_database}" ".${bb}/${id}_filtered.fastq" > ".${bb_p}/${id}_18s_mapped.sam" 2>> "${log}"
samtools view -@ "${cores}" -bS ".${bb_p}/${id}_18s_mapped.sam" > ".${bb_p}/${id}_18s_mapped.bam" 2>> "${log}"
samtools sort -@ "${cores}" -m "${memory_per_core}G" ".${bb_p}/${id}_18s_mapped.bam" -o ".${bb_p}/${id}_18s_mapped_sorted.bam" 2>> "${log}"
samtools index ".${bb_p}/${id}_18s_mapped_sorted.bam" 2>> "${log}"
samtools view -@ "${cores}" -b -F 4 -F 2048 -q "$samtools_quality_filter" ".${bb_p}/${id}_18s_mapped_sorted.bam" | samtools fasta -@ "${cores}" > ".${output_s}/${id}_18s_bin.fasta" 2>> "${log}"
# QC
binned_18s=$(grep -c '^>' ".${output_s}/${id}_18s_bin.fasta")
binned_16s=$(grep -c '^>' ".${output_s}/${id}_16s_bin.fasta")
per_18s_bin=$(echo "scale=5; (($binned_18s / $seq_count_filtered) * 100)" | bc)
per_16s_bin=$(echo "scale=5; (($binned_16s / $seq_count_filtered) * 100)" | bc)
binned_read_representation=$(echo "scale=5; (($per_16s_bin + $per_18s_bin) - 100)" | bc)
read_number_discrepancy=$(echo "scale=5; (($binned_16s + $binned_18s) - $seq_count_filtered)" | bc)

echo "###################################
16s/18s binnning completion stats:
  16s ($per_16s_bin%) $binned_16s reads
  18s ($per_18s_bin%) $binned_18s reads
discrepancy = $binned_read_representation% ($read_number_discrepancy reads)
       --used ${default_database_name} ${default_database_version}--" >> ${log}
echo "" > ".${bb_p}/medaka_stdout_log.txt"
## Medaka iteration 1
echo "5,Medaka iteration 1 - using $binned_16s ($per_16s_bin%) 16s and $binned_18s ($per_18s_bin%) 18s reads" > $progress_file
echo "5,Medaka iteration 1 - using $binned_16s ($per_16s_bin%) 16s and $binned_18s ($per_18s_bin%) 18s reads" >> ${log}
eval "$(conda shell.bash hook)"
conda activate "${medaka}" || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }
echo "Iteration 1 - 16s" >> ".${bb_p}/medaka_stdout_log.txt"
medaka_consensus -i ".${output_s}/${id}_16s_bin.fasta" -d "${fasta_16s_database}" -o ".${bb_p}/${id}_16s_medaka_1" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g -x >> ".${bb_p}/medaka_stdout_log.txt" 2>&1
echo "Iteration 1 - 18s" >> ".${bb_p}/medaka_stdout_log.txt"
medaka_consensus -i ".${output_s}/${id}_18s_bin.fasta" -d "${fasta_18s_database}" -o ".${bb_p}/${id}_18s_medaka_1" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g -x >> ".${bb_p}/medaka_stdout_log.txt" 2>&1
# QC
seq_count_medaka_16s=$(grep -c '^>' ".${bb_p}/${id}_16s_medaka_1/consensus.fasta")
seq_count_medaka_18s=$(grep -c '^>' ".${bb_p}/${id}_18s_medaka_1/consensus.fasta")
percentage_retained_16s_medaka=$(echo "scale=5; $seq_count_medaka_16s / $binned_16s * 100" | bc)
percentage_retained_18s_medaka=$(echo "scale=5; $seq_count_medaka_18s / $binned_18s * 100" | bc)

## Medaka iteration 2
echo "6,Medaka iteration 2 - $seq_count_medaka_16s ($percentage_retained_16s_medaka%) 16s and $seq_count_medaka_18s ($percentage_retained_18s_medaka%) 18s varients found" > $progress_file
echo "6,Medaka iteration 2 - $seq_count_medaka_16s ($percentage_retained_16s_medaka%) 16s and $seq_count_medaka_18s ($percentage_retained_18s_medaka%) 18s varients found" >> ${log}
echo "Iteration 2 - 16s" >> ".${bb_p}/medaka_stdout_log.txt"
medaka_consensus -i ".${bb_p}/${id}_16s_medaka_1/consensus.fasta" -d "${fasta_16s_database}" -o ".${bb_p}/${id}_16s_medaka_2" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g >> ".${bb_p}/medaka_stdout_log.txt" 2>&1
echo "Iteration 2 - 18s" >> ".${bb_p}/medaka_stdout_log.txt"
medaka_consensus -i ".${bb_p}/${id}_18s_medaka_1/consensus.fasta" -d "${fasta_18s_database}" -o ".${bb_p}/${id}_18s_medaka_2" -m "${medaka_model}" -t "${para_cores}" -b "${para_RAM}" -f -g >> ".${bb_p}/medaka_stdout_log.txt" 2>&1
# QC
conda deactivate
cp ".${bb_p}/${id}_16s_medaka_2/consensus.fasta" ".${input_p}/${id}.fasta"
cp ".${bb_p}/${id}_18s_medaka_2/consensus.fasta" ".${input_e}/${id}.fasta"
seq_count_medaka_16s_2=$(grep -c '^>' ".${input_p}/${id}.fasta")
seq_count_medaka_18s_2=$(grep -c '^>' ".${input_e}/${id}.fasta")
percentage_retained_16s_medaka_it2=$(echo "scale=5; $seq_count_medaka_16s_2 / $binned_16s * 100" | bc)
percentage_retained_18s_medaka_it2=$(echo "scale=5; $seq_count_medaka_18s_2 / $binned_18s * 100" | bc)
total_amplicons_CON=$(echo "scale=5; ($seq_count_medaka_16s_2 + $seq_count_medaka_18s_2)" | bc)
echo "###################################
16s/18s error correction stats:
  16s $binned_16s reads reduced to $seq_count_medaka_16s_2 amplicons ($percentage_retained_16s_medaka_it2%)
  18s $binned_18s reads reduced to $seq_count_medaka_18s_2 amplicons ($percentage_retained_18s_medaka_it2%)
  totaling $total_amplicons_CON amplicons from $seq_count_filtered reads" >> ${log}



######		Classification and abundance extraction		######
echo "7,Importing and classifying ${seq_count_medaka_16s_2} ($percentage_retained_16s_medaka_it2% retained) polished 16S amplicons in QIIME2 ($default_database_name $default_database_version)" > $progress_file
echo "7,Importing and classifying ${seq_count_medaka_16s_2} ($percentage_retained_16s_medaka_it2% retained) polished 16S amplicons in QIIME2 ($default_database_name $default_database_version)" >> ${log}
eval "$(conda shell.bash hook)"
conda activate "${biopy}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }
seqtk seq -F I ".${input_p}/${id}.fasta" > ".${input_p}/${id}.fastq" 2>> "${log}" # Convert fasta to fastq as QIIME only functions with sampledata (fastq) not featuredata (fasta)
seqtk seq -F I ".${input_e}/${id}.fasta" > ".${input_e}/${id}.fastq" 2>> "${log}" # Convert fasta to fastq as QIIME only functions with sampledata (fastq) not featuredata (fasta)
conda deactivate
eval "$(conda shell.bash hook)"
conda activate "${qiime}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }

## 16S
# Import
echo -e "sample-id\tabsolute-filepath" > ".${input_p}/manifest.tsv"
echo -e "${id}\t${cd}/${id}${input_p}/${id}.fastq" >> ".${input_p}/manifest.tsv"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ".${input_p}/manifest.tsv" \
  --output-path ".${output_p}/CON_qiime2_input.qza" \
  --input-format SingleEndFastqManifestPhred33V2 \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: RAW data importing to qiime###" >> "${log}"; exit 1; }

# CON data Dereplication
qiime vsearch dereplicate-sequences \
  --i-sequences ".${output_p}/CON_qiime2_input.qza" \
  --o-dereplicated-table ".${output_p}/CON_table-derep.qza" \
  --o-dereplicated-sequences ".${output_p}/CON_rep-seqs-derep.qza" \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: RAW data dereplication###" >> "${log}"; exit 1; }

# Run CON taxonomic classification for SILVA
qiime feature-classifier classify-sklearn \
  --i-classifier "${default_16s_classifier}" \
  --i-reads ".${output_p}/CON_rep-seqs-derep.qza" \
  --o-classification ".${output_p}/taxonomy.qza" \
  --p-confidence "$QIIME_confidence" \
  --p-n-jobs ${Q_cores} \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: 16s taxonomic classification QIIME2###" >> "${log}"; exit 1; }

# Export processed files
qiime tools export --input-path ".${output_p}/taxonomy.qza" --output-path ".${output_p}/taxonomy" >> ${log} 2>&1
qiime tools export --input-path ".${output_p}/CON_rep-seqs-derep.qza" --output-path ".${output_p}/CON_rep-seqs-derep" >> ${log} 2>&1
percent_lost_16s=$(python $tax_to_fasta ".${output_p}/CON_rep-seqs-derep/dna-sequences.fasta" ".${output_p}/taxonomy/taxonomy.tsv" "$fasta_16s_database" ".${output_p}/CON_taxonomy_rep_seq.fasta" "7" "$log" "$cores")
sleep 5

## 18S
echo "8,Importing and classifying ${seq_count_medaka_18s_2} ($percentage_retained_18s_medaka_it2% retained) polished 16S amplicons in QIIME2 ($default_database_name $default_database_version)" > $progress_file
echo "8,Importing and classifying ${seq_count_medaka_18s_2} ($percentage_retained_18s_medaka_it2% retained) polished 16S amplicons in QIIME2 ($default_database_name $default_database_version)" >> ${log}
echo -e "sample-id\tabsolute-filepath" > ".${input_e}/manifest.tsv"
echo -e "${id}\t${cd}/${id}${input_e}/${id}.fastq" >> ".${input_e}/manifest.tsv"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ".${input_e}/manifest.tsv" \
  --output-path ".${output_e}/CON_qiime2_input.qza" \
  --input-format SingleEndFastqManifestPhred33V2 \
 >> ${log} 2>&1 || { sleep 5; echo "###failed here: 18s CON data importing to qiime###" >> "${log}"; exit 1; }

# RAW data Dereplication
qiime vsearch dereplicate-sequences \
  --i-sequences ".${output_e}/CON_qiime2_input.qza" \
  --o-dereplicated-table ".${output_e}/CON_table-derep.qza" \
  --o-dereplicated-sequences ".${output_e}/CON_rep-seqs-derep.qza" \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: 18s CON data dereplication###" >> "${log}"; exit 1; }

# Run CON taxonomic classification for SILVA
qiime feature-classifier classify-sklearn \
  --i-classifier "${default_18s_classifier}" \
  --i-reads ".${output_e}/CON_rep-seqs-derep.qza" \
  --o-classification ".${output_e}/taxonomy.qza" \
  --p-confidence "$QIIME_confidence" \
  --p-n-jobs ${Q_cores} \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: 18s taxonomic classification QIIME2###" >> "${log}"; exit 1; }

# Export processed files
qiime tools export --input-path ".${output_e}/taxonomy.qza" --output-path ".${output_e}/taxonomy" >> ${log} 2>&1
qiime tools export --input-path ".${output_e}/CON_rep-seqs-derep.qza" --output-path ".${output_e}/CON_rep-seqs-derep" >> ${log} 2>&1
percent_lost_18s=$(python $tax_to_fasta ".${output_e}/CON_rep-seqs-derep/dna-sequences.fasta" ".${output_e}/taxonomy/taxonomy.tsv" "$fasta_18s_database" ".${output_e}/CON_taxonomy_rep_seq.fasta" "11" "$log" "$cores")
conda deactivate


### Kraken2 and Bracken
echo "9,RAW 16S reads being assigned to CON 16s (${percent_lost_16s}% lost; 0.00% ideal) reads in Kraken2" > $progress_file
echo "9,RAW 16S reads being assigned to CON 16s (${percent_lost_16s}% lost; 0.00% ideal) reads in Kraken2" >> ${log}
eval "$(conda shell.bash hook)"
conda activate "${biopy}"  || { echo "${er}ERROR:${in} Please ensure config.sh contains correct environment names${r} "; exit 1; }
## 16S
# Generate CON database
stdbuf -oL kraken2-build --add-to-library ".${output_p}/CON_taxonomy_rep_seq.fasta" --db ".${output_p}/${id}_16s_${run_time}-CON" >> ${log} 2>&1
stdbuf -oL kraken2-build --db ".${output_p}/${id}_16s_${run_time}-CON" --download-taxonomy >> ${log} 2>&1
stdbuf -oL kraken2-build --db ".${output_p}/${id}_16s_${run_time}-CON" --kmer-len 15 --minimizer-len 10 --minimizer-spaces 2 --threads ${cores} --build >> ${log} 2>&1 || { sleep 5; echo "###failed here: Kraken 16s database creation###" >> "${log}"; exit 1; }
# Classify 16S RAW Data
stdbuf -oL kraken2 --db ".${output_p}/${id}_16s_${run_time}-CON" \
  --output ".${output_p}/kraken_classification_16s.txt" \
  --report ".${output_p}/kraken_report_16s.txt" \
  --threads ${cores} \
  --confidence "$KrKn_confidence" \
  --minimum-hit-groups "$KrKn_hg" \
  ".${output_s}/${id}_16s_bin.fasta" \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: Kraken matching raw and 16s con data ###" >> "${log}"; exit 1; }
# Recalculate 16s abundances
bracken-build -l "$len_16s" -d ".${output_p}/${id}_16s_${run_time}-CON" -t ${cores} >> ${log} 2>&1
bracken -d ".${output_p}/${id}_16s_${run_time}-CON" -i ".${output_p}/kraken_report_16s.txt" -o ".${output_p}/braken_16s.tsv" -r "$len_16s" -l S -t ${cores} >> ${log} 2>&1 || { sleep 5; echo "###failed here: Bracken matching raw and 16s con data ###" >> "${log}"; exit 1; }

echo "10,RAW 18S reads being assigned to CON 18s (${percent_lost_18s}% lost; 0.00% ideal) reads in Kraken2" > $progress_file
echo "10,RAW 18S reads being assigned to CON 18s (${percent_lost_18s}% lost; 0.00% ideal) reads in Kraken2" >> ${log}
## 18S
# Generate CON database
stdbuf -oL kraken2-build --add-to-library ".${output_e}/CON_taxonomy_rep_seq.fasta" --db ".${output_e}/${id}_18s_${run_time}-CON" >> ${log} 2>&1
stdbuf -oL kraken2-build --db ".${output_e}/${id}_18s_${run_time}-CON" --download-taxonomy >> ${log} 2>&1
stdbuf -oL kraken2-build --db ".${output_e}/${id}_18s_${run_time}-CON" --kmer-len 15 --minimizer-len 10 --minimizer-spaces 2 --threads ${cores} --build >> ${log} 2>&1 || { sleep 5; echo "###failed here: Kraken 18s database creation###" >> "${log}"; exit 1; }
# Classify 18S RAW Data
stdbuf -oL kraken2 --db ".${output_e}/${id}_18s_${run_time}-CON" \
  --output ".${output_e}/kraken_classification_18s.txt" \
  --report ".${output_e}/kraken_report_18s.txt" \
  --threads ${cores} \
  --confidence "$KrKn_confidence" \
  --minimum-hit-groups "$KrKn_hg" \
  ".${output_s}/${id}_18s_bin.fasta" \
  >> ${log} 2>&1 || { sleep 5; echo "###failed here: Kraken matching raw and 18s con data ###" >> "${log}"; exit 1; }
# Recalculate 18s abundances
bracken-build -l "$len_18s" -d ".${output_e}/${id}_18s_${run_time}-CON" -t ${cores} >> ${log} 2>&1
bracken -d ".${output_e}/${id}_18s_${run_time}-CON" -i ".${output_e}/kraken_report_18s.txt" -o ".${output_e}/braken_18s.tsv" -r "$len_18s" -l S -t ${cores} >> ${log} 2>&1 || { sleep 5; echo "###failed here: Bracken matching raw and 16s con data ###" >> "${log}"; exit 1; }

conda deactivate


###	Normalisation, bias correction and merge
echo "11,Normalising, scaling and merging data" > $progress_file
echo "11,Normalising, scaling and merging data" >> ${log}
eval "$(conda shell.bash hook)"
conda activate "$python"
# Remove excess from bracken outputs and re-name headers (python doesnt like bracken headers)
python ${trim_tsv} "${cd}/${id}${output_p}/braken_16s.tsv" "${cd}/${id}${in_m}/braken_16s_trimmed.tsv" "${log}"
python ${trim_tsv} "${cd}/${id}${output_e}/braken_18s.tsv" "${cd}/${id}${in_m}/braken_18s_trimmed.tsv" "${log}"
# bias correction
python ${bias_correction} "${cd}/${id}${in_m}/braken_16s_trimmed.tsv" "$bias_factor_16s" "${cd}/${id}${in_m}/braken_16s_bias_corrected.tsv" "${log}"
python ${bias_correction} "${cd}/${id}${in_m}/braken_18s_trimmed.tsv" "$bias_factor_18s" "${cd}/${id}${in_m}/braken_18s_bias_corrected.tsv" "${log}"
# Count total reads and normalise
python ${normalise} "${cd}/${id}${in_m}/braken_16s_bias_corrected.tsv" "$norm_factor" "${cd}/${id}${in_m}/16s_microbiome.tsv" "${log}"
python ${normalise} "${cd}/${id}${in_m}/braken_18s_bias_corrected.tsv" "$norm_factor" "${cd}/${id}${in_m}/18s_microbiome.tsv" "${log}"
# Merge microbiomes
python ${merge_16s_18s} "${cd}/${id}${in_m}/18s_microbiome.tsv" "${cd}/${id}${in_m}/16s_microbiome.tsv" "${cd}/${id}${out_m}/${run_time}-${id}_microbiome.tsv" "${log}"
conda deactivate

echo "12,Done: see ${cd}/${id}/${run_time}-${id}_microbiome.tsv" > $progress_file 
echo "12,Done: see ${cd}/${id}/${run_time}-${id}_microbiome.tsv" >> ${log}
sleep 2

# Make sure to wait for the Python script to finish before removing the progress file and exiting
wait $PYTHON_PID
rm -f $progress_file


