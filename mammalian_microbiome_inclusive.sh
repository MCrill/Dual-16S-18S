#!/bin/bash
# This filtration approach: includes only organisoms which are known microbiome members (mammalian microbiome), as opposed to filtering out unknown (inclusive).

input=$1
prefix=$2
out_dir=$3

# Create the output file if it doesn't exist
echo > "${out_dir}16s_${prefix}.fasta"
echo > "${out_dir}18s_${prefix}.fasta"

# 16S filter
awk '
/^>/ {
    printit = 0
    if ($0 ~ / Bacteria;/ || $0 ~ / Archaea;/) {
        if ($0 ~ /;Firmicutes;/ || $0 ~ /;Bacteroidetes;/ || $0 ~ /;Actinobacteria;/ ||
            $0 ~ /;Proteobacteria;/ || $0 ~ /;Verrucomicrobia;/ || $0 ~ /;Tenericutes;/ ||
            $0 ~ /;Spirochaetes;/ || $0 ~ /;Cyanobacteria;/ || $0 ~ /;Euryarchaeota;/ ||
            $0 ~ /;Thaumarchaeota;/ || $0 ~ /;Methanosphaera;/ || $0 ~ /;Methanocorpusculum;/ ||
            $0 ~ /;Thermoplasma;/ || $0 ~ /;Sulfolobus;/) {
            printit = 1
        }
    }
}
{
    if (printit) {
        print > "'"${out_dir}16s_${prefix}.fasta"'"
    }
}
' "${input}"

# 18S filter
awk '
/^>/ {
    printit = 0
    if ($0 ~ / Eukaryota;/) {
        if ($0 ~ /;Ascomycota;/ || $0 ~ /;Basidiomycota;/ || $0 ~ /;Mucoromycota;/ ||
            $0 ~ /;Chytridiomycota;/ || $0 ~ /;Glomeromycota;/ || $0 ~ /;Amoebozoa;/ ||
            $0 ~ /;Apicomplexa;/ || $0 ~ /;Ciliates;/ || $0 ~ /;Excavata;/) {
            printit = 1
        }
    }
}
{
    if (printit) {
        print > "'"${out_dir}18s_${prefix}.fasta"'"
    }
}
' "${input}"