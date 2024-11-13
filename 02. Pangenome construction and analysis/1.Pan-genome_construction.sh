#Pan-genome construction

# Run cactus-minigraph to construct a pan-genome graph using  Minigraph Cactus
# - Mount the working directory to /data inside the container
# - Use 20 cores for mapping, with specified memory and disk limits

docker run -v /path/to/your/working_directory:/data --rm -it \
  quay.io/comparative-genomics-toolkit/cactus:v2.4.2 cactus-minigraph \
  --mapCores 20 \                      # Number of cores for mapping
  --defaultMemory 40G \                # Default memory for job allocation
  --defaultCores 4 \                   # Default number of cores for tasks
  --defaultDisk 40G \                  # Default disk space allocation
  --maxMemory 200G \                   # Maximum memory limit
  --maxDisk 400G \                     # Maximum disk space limit
  /data/jobstore \                     # Path to the job store (for storing intermediate results)
  /data/fa.list \                      # List of input FASTA files
  primates.gfa.gz --reference Sscrofall  # Input GFA graph and the reference genome

# Run cactus-graphmap to align sequences to the constructed pan-genome graph
# - Map sequences and produce alignments in PAF format

docker run -v /path/to/your/working_directory:/data --rm -it \
  quay.io/comparative-genomics-toolkit/cactus:v2.4.2 cactus-graphmap \
  --mapCores 20 \                      # Number of cores for mapping
  --defaultMemory 40G \                # Default memory for job allocation
  --defaultCores 4 \                   # Default number of cores for tasks
  --defaultDisk 40G \                  # Default disk space allocation
  --maxMemory 200G \                   # Maximum memory limit
  --maxDisk 400G \                     # Maximum disk space limit
  /data/jobstore \                     # Path to the job store (for storing intermediate results)
  /data/fa.list \                      # List of input FASTA files
  primates.gfa.gz \                    # Input GFA graph
  primates.paf --reference Sscrofall \ # Output PAF file and reference genome
  --outputFasta primates.gfa.fa.gz     # Output FASTA file of aligned sequences

# Run cactus-align to generate a HAL format multiple genome alignment
# - Use the PAF alignments to build a pangenome alignment

docker run -v /path/to/your/working_directory:/data --rm -it \
  quay.io/comparative-genomics-toolkit/cactus:v2.4.2 cactus-align \
  /data/jobstore \                     # Path to the job store (for storing intermediate results)
  /data/fa.list \                      # List of input FASTA files
  primates.paf primates.hal \          # Input PAF file and output HAL file
  --pangenome \                        # Flag for pangenome mode
  --outVG \                            # Output as a variation graph
  --reference Sscrofall                # Reference genome for the alignment

# Run cactus-graphmap-join to join mapped graphs and perform final filtering
# - Output as a GFA variation graph with filtering options

docker run -v /path/to/your/working_directory:/data --rm -it \
  quay.io/comparative-genomics-toolkit/cactus:v2.4.2 cactus-graphmap-join \
  /data/jobstore3 \                    # Path to the job store for this stage
  --vg /data/primates.vg \             # Input VG file
  --defaultMemory 300G \               # Default memory for job allocation
  --defaultCores 10 \                  # Default number of cores for tasks
  --defaultDisk 200G \                 # Default disk space allocation
  --maxCores 16 \                      # Maximum number of cores
  --maxMemory 500G \                   # Maximum memory limit
  --maxDisk 500G \                     # Maximum disk space limit
  --logFile cactus.log2 \              # Log file for monitoring progress
  --outDir /data/primates-pg2 \        # Output directory for the results
  --outName primates-pg2 \             # Name of the output files
  --reference Sscrofall \              # Reference genome for mapping
  --vcf --giraffe clip filter \        # Output VCF format, use Giraffe mapper, clip/filter options
  --chrom-vg                           # Produce VG files for each chromosome
# Optional parameters for advanced filtering:
# --wlineSep "." --clipLength 10000 --clipNonMinigraph --vgClipOpts "-d 2 -m 1000"


