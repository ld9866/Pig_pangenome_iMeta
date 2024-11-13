# Population genetic analysis

# Step 1: Convert VCF file to PLINK format
plink1.9 --vcf input_data.vcf.gz --out output_prefix  # Convert VCF to PLINK format (BED/BIM/FAM files)

# Step 2: Extract SNPs and remove indels
plink1.9 --bfile output_prefix --snps-only --out SNP_only --make-bed  # Keep only SNPs and create a new binary file

# Step 3: Convert PLINK binary files to MAP and PED format
plink1.9 --bfile SNP_only --recode --out SNP_data  # Recode binary files to MAP and PED format

# Step 4: Run Admixture analysis for K values from 1 to 5
for K in 1 2 3 4 5; do
    admixture --cv SNP_data.bed $K | tee log${K}.out  # Perform Admixture analysis with cross-validation and save logs
done

library(ggplot2)
data<-read.table("testacc.eigenvec",header=T)
ggplot(data,aes(x=pc1,y=pc2))+ geom_point()

(ggplot2)
pca <- read.table("pca_info.txt",sep = "\t",header = T)
ggplot(pca, aes(x=pca1,y=pca2)) +
        geom_point(aes(color=pop, shape=pop),size=1.5)+
        labs(x="PC1",y="PC2")+theme_bw()+theme(legend.title = element_blank())


# Loop through VCF files to submit jobs for calculating nucleotide diversity (pi)
for i in /path/to/your/vcf/files/*.raw.snp.vcf.gz
do 
    # Submit job using jsub with specified resources and output settings
    jsub -n 2 -M 100000000 -o $(basename $i .raw.snp.vcf.gz).o -e $(basename $i .raw.snp.vcf.gz).e -J $(basename $i .raw.snp.vcf.gz).pi \
         sh calpi_script.sh $(basename $i .raw.snp.vcf.gz)
done

# Loop through chromosomes 1 to 24 to create XPEHH scripts
for i in {1..24}
do 
    echo $i  # Print the chromosome number
    echo "#!/bin/bash" > $i.xpehh.sh  # Create a shell script for each chromosome
    # Extract samples for population A (Swamp)
    echo "vcftools --gzvcf /path/to/chr$i.imp.phase.vcf.gz --keep /path/to/sample_A_list --recode --out Swamp.chr$i" >> $i.xpehh.sh
    # Extract samples for population B (River)
    echo "vcftools --gzvcf /path/to/chr$i.imp.phase.vcf.gz --keep /path/to/sample_B_list --recode --out River.chr$i" >> $i.xpehh.sh
    # Add IDs to VCF files
    echo "perl add.id.vcf.pl Swamp.chr$i.recode.vcf Swamp.chr$i.vcf" >> $i.xpehh.sh
    echo "perl add.id.vcf.pl River.chr$i.recode.vcf River.chr$i.vcf" >> $i.xpehh.sh
    # Generate map file for Selscan
    echo "awk '{print \$1,\$3,\$2,\$2}' OFS=\"\t\" Swamp.chr$i.vcf > chr$i.map" >> $i.xpehh.sh
    # Run Selscan to calculate XP-EHH
    echo "/path/to/selscan --xpehh --vcf Swamp.chr$i.vcf --vcf-ref River.chr$i.vcf --map chr$i.map --cutoff 0.01 --out Swamp_River.$i.xpehh --threads 4" >> $i.xpehh.sh
    chmod 755 $i.xpehh.sh  # Make the script executable
done

# Normalize XP-EHH results
for i in {1..24}
do
    echo $i  # Print the chromosome number
    /path/to/selscan/bin/linux/norm --xpehh --files Swamp_River.$i.xpehh.xpehh.out --bp-win --winsize 200000
done

# Perform window-based processing of normalized results
for i in {1..24}
do
    echo $i  # Print the chromosome number
    # Create a shell script for window processing
    echo "#!/bin/bash" > $i.XPEHH.window.sh
    echo "perl XPEHH.window.pl /path/to/chr$i.fa.fai XL_CN.$i.xpehh.xpehh.out.norm 50000 20000 XL_CN.$i.50k.norm" >> $i.XPEHH.window.sh
    chmod 755 $i.XPEHH.window.sh  # Make the script executable
    # Submit the window processing job using jsub
    jsub -q jynodequeue -R "rusage[res=1]span[hosts=1]" -J $i.XPEHH.window -e $i.XPEHH.window.e -o $i.XPEHH.window.o -n 4 -M 40000000 ./$i.XPEHH.window.sh
done

# Compile results from window-based processing into a single file
for i in {1..24}
do 
    echo $i  # Print the chromosome number
    # Format output and append to a combined results file
    less GZB_CN.$i.xpehh.xpehh.out.norm.50kb.windows | awk -F '\t' '{print "'$i'""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' >> GZB_CN.xpehh_50K
done




# Calculate Fst using VCFtools with specified filters and window parameters
vcftools --gzvcf ./1.beagle.vcf.gz --weir-fst-pop high.33.txt --weir-fst-pop low.195.txt \
         --fst-window-size 50000 --fst-window-step 25000 --out chr1.elevation --max-missing 0.9 --maf 0.05

# Remove the header line from each .fst file
for i in *.fst; do sed -i '1d' $i; done 

# Concatenate all windowed Fst results into one file
cat *.Fst_cashmere.windowed.weir.fst > Fst_cashmere.windowed.weir.fst

# Add a header to the Fst result file
# CHROM: Chromosome, BIN_START: Start of window, BIN_END: End of window
# N_VARIANTS: Number of variants in window, WEIGHTED_FST: Weighted Fst, MEAN_FST: Mean Fst
echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST" > header.txt
cat header.txt Fst_cashmere.windowed.weir.fst > Fst_cashmere.windowed_with_header.weir.fst
rm header.txt  # Remove temporary header file





# Plot Manhattan plot using a custom Python script
# infile: Input file, chr-col: Chromosome column, loc-col: Location column
# val-col: Value column (Fst), outfile: Output file for the plot
python plot_Manhattan.py --infile FST_IMF.windowed.weir.fst --chr-col CHROM \
                         --loc-col BIN_START --val-col WEIGHTED_FST --outfile FST.Manhattan.pdf \
                         --xlabel pos --ylabel fst

# Run XPCLR for population differentiation analysis
xpclr --format txt --popA chr1.pop.high.geno --popB chr1.pop.low.geno \
      --map chr01.snp --out chr01.out --chr 1 -w1 0.005 200 2000 1

# Explanation of XPCLR parameters:
# -out: Output file
# -format: Input file format (options: vcf, hdf5, txt); details at: https://github.com/hardingnj/xpclr
# -popA: Sample names for population A
# -popB: Sample names for population B
# -w1: Calculation window size, number of SNPs in each window, and step size


# Build a custom SnpEff database for signal selection mutation site annotation
java -jar /path/to/snpEff/snpEff.jar build -gtf22 -c /path/to/snpEff.config -v AT_10

# Explanation of parameters:
# -gtf22: Use GTF format version 22 for annotation
# -c: Path to the SnpEff configuration file
# -v: Verbose mode for more detailed output
# AT_10: The name of the genome version or database being built

