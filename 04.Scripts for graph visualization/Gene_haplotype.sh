#!/bin/bash

# Extract variation information within the entire gene region
# Use bcftools to filter and extract variants on a specified chromosome and region
# Input file: compressed VCF (e.g., input_variants.vcf.gz)
# Output: compressed VCF file (e.g., extracted_region.vcf.gz)
echo "Extracting variants in the specified gene region..."
time bcftools view -r <chromosome>:<start_position>-<end_position> input_variants.vcf.gz -Oz -o extracted_region.vcf.gz

# Index the output VCF file to allow for quick access and queries
echo "Indexing the VCF file..."
bcftools index extracted_region.vcf.gz

# Plot a haplotype heatmap using a Python script
# Parameters explained:
# --vcffile: Path to the input VCF file
# --maf: Minor Allele Frequency threshold for filtering
# --querylist: File containing a list of samples to include in the analysis
# --ticklabelsize: Font size for tick labels on the heatmap
# --groupfile: File specifying population group information for samples
# --region: Specific region to plot, defined by start and end positions
# --outfile: Name of the output PDF file for the heatmap
echo "Generating haplotype heatmap..."
python3 /path/to/plot_vcfHeatmap.py \
  --vcffile extracted_region.vcf.gz \
  --maf <maf_threshold> \
  --querylist sample_list.txt \
  --ticklabelsize 1 \
  --groupfile population_info.tsv \
  --region <chromosome>:<start_position>-<end_position> \
  --outfile output_heatmap.pdf

echo "Script completed successfully!"
