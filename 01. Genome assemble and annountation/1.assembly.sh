##Genome assembly script

##1.Genome assembly using hifiasm

#install
conda create --name hifiasm
conda activate hifiasm
conda install -c bioconda hifiasm

READS=/path/to/your/hifi_reads.fasta.gz # HiFi reads path
OUTPUT_PREFIX=genome_assembly            # Output file prefix
THREADS=32                                # Number of threads to use

# run Hifiasm
hifiasm -o ${OUTPUT_PREFIX} -t ${THREADS} ${READS}

#### The FASTA file can be produced from GFA as follows
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa


#2.Hi-C mounted chromosome-level genome
# juicer

#build index
bwa index sample.hic.p_ctg.fa

#running
python ./juicer-master/misc/generate_site_positions.py MboI draft.genome sample1.hic.p_ctg.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' draft.genome_MboI.txt > draft.genome.chrom.sizes

awk -f Software/3d-dna-master/edit/edit-fasta-according-to-new-cprops.awk sample1.hic.p_ctg.rawchrom.cprops test.hic.p_ctg.fa
python ./juicer-master/misc/generate_site_positions.py MboI draft.genome sample1.fa
awk -f ./3d-dna-master/utils/wrap-fasta-sequence.awk sample1.fa > wrapped.fasta



#3.Use ragtag to mount chromosomes to pseudo-chromosome level
ragtag.py correct Ref.fasta sample.contig.fasta
ragtag.py scaffold Ref.fasta ragtag_output/ragtag.correct.fasta -o final_result
ragtag.py patch ./final_result/ragtag.scaffold.fasta sample.contig.fasta -t 16
ragtag.py merge -b hic.bam query.fasta out_*/*.agp

###Checks that the provided AGP file is in the correct format
ragtag.py agpcheck <asm1.agp> [<asm2.agp> ... <asmN.agp>]

# 4.Modify the nomenclature of each chromosome in the genome

### Export genome chromosome number
seqkit seq wrapped_HiC.fasta -n > HiC_scaffold.list #Print sequence ID full name
### Write a list to extract specific chromosomes
seqtk subseq wrapped_HiC.fasta list > wrapped_HiC.extract.fasta
### Modify the chromosome number of the unanchored chromosome to chr
seqtk rename wrapped_HiC.fasta chr > wrapped_HiC.rename.fasta

#5.Genome quality assessment

#QUAST-LG CGAL QUAST REAPR CEGMA GenomeRibbon LAP
./quast-5.1.0rc1/quast.py -o ./QUAST -t 20 sample1.fasta
busco -m genome -i sample1.fasta -o sample1 -l ./busco_downloads/lineages/vertebrata_odb10 -c 10
