#!/bin/bash

source config.sh
# Easy edit!
prefix=$1
ex_type=$2
filter="mammalian_microbiome_inclusive"
if [ -z "$1" ]; then
    echo "${er}ERROR:${in} please provide prefix to match fasta (which must be in /databases/) ${r}"
    exit 0
fi

filter_path="${w_d}/bin/scripts/${filter}.sh"
input="${w_d}/bin/databases/${prefix}.fasta"
database_dir="${w_d}/bin/databases/"

if [ "$ex_type" = "0" ]; then
# GO (0 - manual)
echo "${in}Generating databases using filter: ${er}${filter}${r}, this should take between 10-30 mins"
# Go (1- build)
else
echo "${in}This script is currently formatted to function with ${er}${default_database_name} ${default_database_version}${in}, using filter preset: ${er}${filter}${r}"
echo "${in}Edit update-database.sh to modify this pre-set${r}
"
fi

# Intitate awk for 16s and 18s specific groups
bash "${filter_path}" "$input" "$prefix" "$database_dir"

# Count reads
total_16s=$(grep -c '^>' "${database_dir}16s_${prefix}.fasta")
total_18s=$(grep -c '^>' "${database_dir}18s_${prefix}.fasta")
total_retain=$((total_16s + total_18s))
silva_total=$(grep -c '^>' "${input}")
discarded_total=$((silva_total - total_retain))
percentage_retain=$(awk "BEGIN {printf \"%.5f\", (${total_retain}/${silva_total})*100}")

# Display the counts and percentage
if [ "$discarded_total" -gt 1 ]; then
    echo "${su}Success:${in} Database filtration complete, ${percentage_retain}% read retention"
    echo "${in} 16S = ${total_16s} amplicons ${r}"
    echo "${in} 18S = ${total_18s} amplicons ${r}"
else
    echo "${er}ERROR:${in} Database filtration failed"
    echo "${in} 16S = ${total_16s} ${r}"
    echo "${in} 18S = ${total_18s} ${r}"
    echo "${in} Discarded | total = $discarded_total | $silva_total ${r}"
    printf "${in}Press ${r}Enter${in} to continue regardless?${r}"
fi

# 16S + 18S = FILERED
cat "${database_dir}18s_${prefix}.fasta" "${database_dir}16s_${prefix}.fasta" > "${database_dir}filtered_${prefix}.fasta"

# Function to update the config file
update_config() {
    local var_name=$1
    local new_value=$2
    if grep -q "^export ${var_name}=" "../../config.sh"; then
        sed -i "s|^export ${var_name}=.*|export ${var_name}=\"${new_value}\"|" "../../config.sh"
    else
        echo "${er}ERROR:${in} failed to update config for ${var_name} ${r}"
    fi
}
# Update config with the new database paths
update_config "fasta_16s_database" "${database_dir}16s_${prefix}.fasta"
update_config "fasta_18s_database" "${database_dir}18s_${prefix}.fasta"
update_config "fasta_filtered_database" "${database_dir}filtered_${prefix}.fasta"
update_config "kraken_16s_database" "${database_dir}16s_${prefix}"
update_config "kraken_18s_database" "${database_dir}18s_${prefix}"

# Identify database version, and update config is needed
prefix_version=$(echo "$prefix" | grep -oP '\d+\.\d+')
if [ "$default_database_version" != "$prefix_version" ]; then
    update_config "default_database_version" "$prefix_version"
fi

# verify kraken database isnt being added too
check_and_delete() {
    local db_path=$1
    if [ -d "${db_path}" ]; then
        echo "${er}WARNING:${in} deleting previous database before continuing (avoiding concatination)${r}"
        rm -rf "${db_path}"    
    fi
}
# Kraken 16S
echo "   ${in}Kraken 16S database construction: ${r}"
check_and_delete "${database_dir}16s_${prefix}"
kraken2-build --add-to-library "${database_dir}16s_${prefix}.fasta" --db "${database_dir}16s_${prefix}"
kraken2-build --db "${database_dir}16s_${prefix}" --download-taxonomy
kraken2-build --build --db "${database_dir}16s_${prefix}" --kmer-len 35 --minimizer-len 31 --minimizer-spaces 7 --threads ${all_cores}
# Kraken 18S
echo "   ${in}Kraken 18S database construction: ${r}"
check_and_delete "${database_dir}18s_${prefix}"
kraken2-build --add-to-library "${database_dir}18s_${prefix}.fasta" --db "${database_dir}18s_${prefix}"
kraken2-build --db "${database_dir}18s_${prefix}" --download-taxonomy
kraken2-build --build --db "${database_dir}18s_${prefix}" --kmer-len 35 --minimizer-len 31 --minimizer-spaces 7 --threads ${all_cores} 
