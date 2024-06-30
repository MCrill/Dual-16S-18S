#!/bin/bash

# Find directories
source config.sh
sub="${w_d}/bin/scripts"
cd=$(pwd)
current_time=$(date '+%H-%M-%S_%d%b%Y')

###	 HELP	###
## General
# Display help regardless of bash syntax error in -help
function display_help() {
    echo "${in} Usage: nap <dorado/pipe/stat/update-database> [options]...${r} "
    echo "${in} Try: nap <tool name> -h or --help for more specific information ${r} "
}

if [ "$1" = "-help" ] || [ "$1" = "-h" ]; then
    display_help
    exit 0
fi
# DORADO
if [[ "$1" == "dorado" ]]; then
    if [[ "$2" == "-h" || "$2" == "--help" ]]; then
        echo "${er} Usage:${in} nap dorado (simply be in the directory containing your ./pod5/ ${r}"
        exit 0
    fi
fi
# PIPE
if [[ "$1" == "pipe" ]]; then
    if [[ "$2" == "-h" || "$2" == "--help" ]]; then
        echo "${er} Usage:${in} nap pipe <(barcode number:)01> <corresponding_sample_id_1> <02> <corresponding_sample_id_2>.... ${r}"
        exit 0
    fi
fi
# Update database
if [[ "$1" == "update-database" ]]; then
    if [[ "$2" == "-h" || "$2" == "--help" ]]; then
        echo "${er} Usage:${in} nap update-database; simply be in the database folder along with your database.fasta${r}"
        exit 0
    fi
fi

###	Data identification	###
TOOL_NAME=$1
shift

# Script 'nap dorado'
if [ "$TOOL_NAME" = "dorado" ]; then
  bash "${sub}/dorado.sh"
fi

# Script 'nap dorado-hq'
if [ "$TOOL_NAME" = "dorado-hq" ]; then
  bash "${sub}/dorado-strict.sh"
fi

# Script 'nap update-database'
if [ "$TOOL_NAME" = "update-database" ]; then
  bash "${sub}/update-database.sh"
fi

# Script to run 'nap pipe' for each sample given and export variable
if [ "$TOOL_NAME" = "pipe" ]; then
    # Check if there are sufficient arguments for sample processing
    if [ $# -lt 2 ]; then
        echo "${in}Format: nap pipe <01> <corresponding_sample_id_1> <02> <corresponding_sample_id_2> ... for up to 12 samples (please copy barcode file names) ${r} "
        display_help
        exit 1
    fi

    # Prepare log file
    log_file="${w_d}/bin/logs/$(date +'%Y%m%d_%H-%M-%S')_pipe_log.txt"
    echo "Date: $(date)" > "$log_file"
    echo -e "\nSample ID\tBarcode\tFile Location" >> "$log_file"
    # Run nap pipe for each file sequentially
    while [ $# -ge 2 ]; do
        barcode="$1"
        sample_id="$2"
        sample_file=$(find ${cd}/raw_data -name "*barcode${barcode}.fastq" -type f -print -quit)

        if [ -n "$sample_file" ]; then
            # Pass the log file path as the third argument to pipe.sh
            echo -e "${sample_id}\t${barcode}\t${sample_file}" >> "$log_file"
            bash "${sub}/pipe.sh" "${sample_file}" "${sample_id}" "${log_file}"
        else
            echo "${er}ERROR: ${in}File for barcode ${cd}/raw_data/*${barcode}.fastq not found ${r}"
            display_help
            exit 1
        fi

        shift 2
    done
    echo -e "${in}Log file created: ${log_file} ${r}"
fi

## PCoA python
if [ "$TOOL_NAME" = "PCoA" ]; then
  eval "$(conda shell.bash hook)"
  conda activate $python
  # Check for help option
  if [[ "$1" == "-h" || "$1" == "--help" ]]; then
      echo "Usage: nap PCoA <metric> <input_file1.tsv> <input_file2.tsv> ... or <./path/to/all/*.tsv>"
      echo "  metric           Metric to use for beta-diversity calculation (e.g., braycurtis, jaccard)."
      exit 0
  fi
  # Check minimum number of arguments
  if [ $# -lt 2 ]; then
      echo "Error: Not enough arguments provided."
      echo "Usage: nap PCoA <metric> <input_file1.tsv> <input_file2.tsv> ... or <./path/to/all/*.tsv>"
      echo "  metric           Metric to use for beta-diversity calculation (e.g., braycurtis, jaccard)."
      exit 1
  fi

  # Extract arguments
  METRIC="$1"
  shift 1

  # Validate files and construct the input files string
  INPUT_FILES=""
  FILE_COUNT=0
  while [ $# -gt 0 ]; do
      if [[ "$1" == *"*"* ]]; then
          # Expand wildcard and check files
          EXPANDED_FILES=$(find . -path "$1")
          for file in $EXPANDED_FILES; do
              if [ ! -f "$file" ]; then
                  echo "Error: File $file not found."
                  display_help
                  exit 1
              fi
              # Remove any leading './' in the file path
              file="${file#./}"
              INPUT_FILES+="${cd}/$file "
              ((FILE_COUNT++))
          done
          if [ $FILE_COUNT -eq 0 ]; then
              echo "Error: No files match your wildcard $1."
              display_help
              exit 1
          fi
      else
          # Check for a single file
          if [ ! -f "$1" ]; then
              echo "Error: File $1 not found."
              display_help
              exit 1
          fi
          # Remove any leading './' in the file path
          file="${1#./}"
          INPUT_FILES+="${cd}/$file "
          ((FILE_COUNT++))
      fi
      shift
  done

  echo "Detected $FILE_COUNT input files."
  echo "python "$PCoA" "${cd}/PCoA-${METRIC}_${current_time}.png" "$METRIC" -i "$INPUT_FILES""
  # Prepare the Python command line and execute it
  #python "$PCoA" "${cd}/PCoA-${METRIC}_${current_time}.png" "$METRIC" -i "$INPUT_FILES" 
  #xdg-open "${cd}/PCoA-${METRIC}_${current_time}.png"
fi






# DELETE LATER
# Script to run 'nap pipe' for each sample given and export variable
if [ "$TOOL_NAME" = "pipe-1" ]; then
    # Check if there are sufficient arguments for sample processing
    if [ $# -lt 2 ]; then
        echo "${in}Format: nap pipe <01> <corresponding_sample_id_1> <02> <corresponding_sample_id_2> ... for up to 12 samples (please copy barcode file names) ${r} "
        display_help
        exit 1
    fi

    # Prepare log file
    log_file="${w_d}/bin/logs/$(date +'%Y%m%d_%H-%M-%S')_pipe_log.txt"
    echo "Date: $(date)" > "$log_file"
    echo -e "\nSample ID\tBarcode\tFile Location" >> "$log_file"
    # Run nap pipe for each file sequentially
    while [ $# -ge 2 ]; do
        barcode="$1"
        sample_id="$2"
        sample_file=$(find ${cd}/raw_data -name "*barcode${barcode}.fastq" -type f -print -quit)

        if [ -n "$sample_file" ]; then
            # Pass the log file path as the third argument to pipe.sh
            echo -e "${sample_id}\t${barcode}\t${sample_file}" >> "$log_file"
            bash "${sub}/pipe-1.sh" "${sample_file}" "${sample_id}" "${log_file}"
        else
            echo "${er}ERROR: ${in}File for barcode ${cd}/raw_data/*${barcode}.fastq not found ${r}"
            display_help
            exit 1
        fi

        shift 2
    done
    echo -e "${in}Log file created: ${log_file} ${r}"

fi




