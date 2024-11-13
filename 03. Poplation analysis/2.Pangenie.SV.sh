#Structural variation detection

# Activate the Conda environment for PanGenie
conda activate pangenie

# Use PanGenie-index to create an index file
# -v specifies the input VCF file containing known structural variants
# -r specifies the reference genome sequence file (in FASTA format)
# -t specifies the number of threads to use (10 in this case)
# -o specifies the output prefix for the index files (e.g., "Pig")
PanGenie-index -v ./path/Data/PanGenie/27.chr20.pang.SV.vcf \
               -r ./path/Data/PanGenie/Sscrofall.chr20.fa \
               -t 10 \
               -o Pig

# Perform structural variant genotyping using PanGenie
# Use process substitution <() to decompress the paired-end FASTQ files
# -i specifies the input paired-end reads (decompressed from .fastq.gz)
# -r specifies the reference genome sequence file
# -v specifies the input VCF file with known structural variants
# -o specifies the output prefix
# -s specifies the sample name
# -t and -j specify the number of threads to use (both set to 10)
time PanGenie -i <(zcat ./path/Data/PanGenie/sample1.clean_1.fastq.gz \
                   ./path/Data/PanGenie/sample1.clean_2.fastq.gz) \
              -r ./path/Data/PanGenie/Sscrofall.chr20.fa \
              -v ./path/Data/PanGenie/27.chr20.pang.SV.vcf \
              -o sample1 \
              -s sample1 \
              -t 10 \
              -j 10

#Using manta to detect SV in next-generation sequencing data

# Create a new Conda environment named "manta"
conda create -n manta

# Activate the "manta" environment
conda activate manta

# Install Manta from the Bioconda channel
conda install -c bioconda manta

# Generate configuration scripts for Manta in a loop
for i in `cat list`
do
  configManta.py \
  --bam ${i}.bam \  # Specify the input BAM file for each sample
  --referenceFasta ./data/pangenome/mapping/Reference/GCF_000003025.6_Sscrofa11.1_genomic.fna \  # Specify the reference genome in FASTA format
  --runDir ./${i}.manta  # Specify the output directory for Manta's run scripts
done




