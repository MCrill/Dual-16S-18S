# Dual 16S & 18S rRNA Sequencing Pipeline

Bioinformatics Msc Dissertation Project for dual 16S and 18S rRNA Sequencing Pipeline

 This pipeline integrates various tools and scripts designed to process raw sequencing data, perform quality control, filter reads, identify chimeras, and classify sequences taxonomically.

The pipeline begins with the config.sh script, which sets up essential parameters and environment variables needed for the analysis. The mammalian_microbiome_inclusive.sh script filters input sequences to include only known microbiome members from mammalian sources. The update-database.sh script updates the database with filtered sequences, and the nap.sh script orchestrates the execution of different stages of the pipeline, including the filtering, chimera detection, and classification.

Key features of the pipeline include:

Configuration Flexibility: Easily customizable parameters for hardware utilization, quality thresholds, and specific amplicon handling.

Comprehensive Filtration: Scripts like mammalian_microbiome_inclusive.sh ensure only relevant microbiome sequences are included.

Automated Database Updates: update-database.sh automates the integration of new sequence data into the existing database.

End-to-End Analysis: The pipe.sh script finalizes the pipeline, performing thorough processing from raw data to classified output.

This pipeline is suitable for researchers looking to analyze large-scale microbiome data with a focus on both bacterial and eukaryotic components, ensuring high-quality, reproducible results. For detailed usage instructions and parameter adjustments, refer to the individual script documentation within the repository.
