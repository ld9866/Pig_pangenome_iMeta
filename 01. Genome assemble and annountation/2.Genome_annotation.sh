#Genome annotation


#1.Masking repetitive sequences
# Perform genome-wide repeat sequence analysis using RepeatModeler by building a repeat sequence database
BuildDatabase -name example_lib sample1.fasta
RepeatModeler -pa 10 -database example_lib -LTRStruct >& RepeatModeler.run.log

# Use RepeatMasker to annotate repeat sequences in the genome using a custom repeat library
RepeatMasker sample1.fasta -lib sample1.families.fa -pa 10 -poly -html -gff 1>log.o.txt 2>log.e.txt




#2.Use BRAKER for gene prediction and annotation, integrating evidence from RNA-Seq and protein homology data


#!/bin/bash

# Build the HISAT2 index for the reference genome
hisat2-build -p 4 genome.fa genome

# Align RNA-Seq paired-end reads to the reference genome using HISAT2
hisat2 -p 10 -x ./genome -1 Sample.R1.fastq -2 Sample.R2.fastq --rna-strandness RF --fr -S Sample.sam

# Convert the SAM file to a BAM file, sort, and index the BAM file using Samtools
samtools view -uS Sample.sam | samtools sort - -o Sample.sorted.bam
samtools index Sample.sorted.bam


# Run BRAKER1 for gene prediction using multiple RNA-Seq BAM files as evidence
braker.pl --genome=sample_genome.fasta.masked --species=sample_species \
          --bam sample1.bam,sample2.bam,sample3.bam,sample4.bam,sample5.bam,sample6.bam,sample7.bam,sample8.bam, \
          sample9.bam,sample10.bam,sample11.bam,sample12.bam,sample13.bam,sample14.bam,sample15.bam,sample16.bam, \
          sample17.bam,sample18.bam,sample19.bam,sample20.bam,sample21.bam,sample22.bam,sample23.bam,sample24.bam, \
          sample25.bam,sample26.bam \
          --softmasking --cores 20 --BAMTOOLS_PATH=/path/to/bamtools/


# Step 1: Generate high-quality CCS reads from subreads.bam using minimum read quality threshold of 0.9
ccs [movie].subreads.bam [movie].ccs.bam --min-rq 0.9

# Step 2: Perform transcript demultiplexing using Lima with Iso-Seq primers
lima --isoseq --dump-clips --no-pbi --peek-guess -j 24 ccs.bam primers.fasta demux.bam

# Step 3: Refine the consensus reads to improve accuracy and add PolyA tail filtering
isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta flnc.bam

# Step 4: Cluster the refined full-length non-concordant (FLNC) reads
isoseq3 cluster flnc.bam polished.bam --verbose --use-qvs

# Step 5: Convert the final FLNC BAM file to FASTQ format for further analysis
bamtools convert -format fastq -in flnc.bam > flnc.fastq




# Perform gene prediction using BRAKER with protein evidence and a soft-masked genome

# Run BRAKER2 on a sample genome with protein sequence data
braker.pl --genome=sample1_genome.fa --species=sample_species --prot_seq=pig_proteins.fa --softmasking --cores 20


# Run TSEBRA to merge BRAKER1 and BRAKER2 predictions using RNA-Seq and protein evidence
# Input files include Augustus hints, GFF files, and repeat modeler results

/home/user/Software/TSEBRA_longreads/TSEBRA-long_reads/bin/tsebra.py \
    -g sample1.augustus.hints.gtf,sample2.augustus.hints.gtf \  # Augustus hint files from RNA-Seq and protein data
    -e sample1.hintsfile.gff,sample2.hintsfile.gff \  # GFF files containing hints from RNA-Seq and protein evidence
    -l sample_repeat_modeler.gtf \  # Repeat modeler output with repeat and transcript information
    --cfg /home/user/Software/TSEBRA_longreads/TSEBRA-long_reads/config/long_reads.cfg \  # Configuration file for TSEBRA
    -o output_file.gtf  # Output file containing the merged gene predictions in GTF format


# Step 1: Extract CDS and protein sequences from the TSEBRA GTF file
# Using gffread to extract CDS sequences in FASTA format and corresponding protein sequences
gffread tsebra2.gtf -g genome.fasta.masked -x cds.fasta -y protein.fa

# Step 2: Perform sequence clustering to remove redundancy
# Using CD-HIT to cluster protein sequences based on 90% identity
cd-hit -i protein_0.5.fa -o protein_0.6.fa -c 0.9 -d 0 -T 8


# Define the input protein file (protein sequences after CD-HIT clustering)
INPUT_PROTEIN="protein.cdhit.fa"

# Define the output log file for all commands
LOG_FILE="log"

# Define paths to the databases (use variables to represent database paths)
NR_DB="/path/to/nr.dmnd"          # Nr database
SWISS_DB="/path/to/uniprot.dmnd"  # Swiss-Prot database
TREMDB_DB="/path/to/trembl.dmnd"  # TrEMBL database

# Define output files for each database match result
OUTPUT_NR="match.Nr.tsv"
OUTPUT_SWISS="match.Swiss.port.tsv"
OUTPUT_TREMDB="match.TrEMBL.tsv"

# Nr database annotation
diamond blastp --max-target-seqs 1 \
               --evalue 1e-5 \
               --max-hsps 1 \
               -d $NR_DB \
               --threads $THREADS \
               -q $INPUT_PROTEIN \
               -o $OUTPUT_NR \
               &> $LOG_FILE

# Swiss-Prot database annotation
diamond blastp --max-target-seqs 1 \
               --evalue 1e-5 \
               --max-hsps 1 \
               -d $SWISS_DB \
               --threads $THREADS \
               -q $INPUT_PROTEIN \
               -o $OUTPUT_SWISS \
               &>> $LOG_FILE

# TrEMBL database annotation
diamond blastp --max-target-seqs 1 \
               --evalue 1e-5 \
               --max-hsps 1 \
               -d $TREMDB_DB \
               --threads $THREADS \
               -q $INPUT_PROTEIN \
               -o $OUTPUT_TREMDB \
               &>> $LOG_FILE



# Run tRNAscan-SE for tRNA prediction
tRNAscan-SE \
  -E \                       # Enable eukaryotic tRNA detection mode
  -o tRNA.out \              # Specify output file for predicted tRNAs in tabular format
  -f rRNA.ss \               # Save tRNA secondary structure predictions to this file
  -m tRNA.stats \            # Output statistics about the tRNA predictions to this file
  -b tRNA.bed \              # Output tRNA locations in BED format for genome browsers
  -j tRNA.gff \              # Save tRNA annotation in GFF3 format
  -a tRNA.fa \               # Output tRNA sequences in FASTA format
  -l log \                   # Write the log file with details of the run
  --thread 10 \              # Use 10 threads for parallel processing
  sample1.fasta                # Input genome sequence file for the  breed
