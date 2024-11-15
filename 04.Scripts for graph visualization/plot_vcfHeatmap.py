import os
import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def loadgroup(groupfile):
    """
    Load sample grouping information from a file.
    The file should contain two columns: [Sample ID, Group Name].
    """
    id2groups = {}
    with open(groupfile) as f:
        for line in f:
            tline = line.strip().split()
            id2groups[tline[0]] = tline[1]
    return id2groups

@click.command()
@click.option('--vcffile', help='Input genotype VCF file')
@click.option('--outlist', help='List of samples for the outgroup (used to determine the derived allele)', default=None)
@click.option('--querylist', help='List of samples to be plotted (one sample per line in a text file)')
@click.option('--region', help='Region to plot, e.g., 12:1000-2000', default=None)
@click.option('--maf', help='Minimum allele frequency threshold, default=0.2', type=float, default=0.2)
@click.option('--groupfile', help='Sample grouping file, two columns: [Sample ID, Group Name]')
@click.option('--outfile', help='Output image file')
@click.option('--figsize', nargs=2, type=float, help='Figure dimensions, default=(20, 5)', default=(20, 5))
@click.option('--ticklabelsize', help='Font size for tick labels, default=1', default=1, type=int)
@click.option('--dpi', help='Image resolution, default=500', default=500)
def main(vcffile, outlist, querylist, region, maf, groupfile, outfile, figsize, ticklabelsize, dpi):
    """
    Generate a haplotype heatmap from a VCF file.
    The most frequent allele in the outgroup is designated as the ancestral allele.
    If no outgroup is provided, the reference and alternate alleles from the VCF are used.
    """
    querysamples = [x.strip() for x in open(querylist)]
    if outlist:
        outsamples = [x.strip() for x in open(outlist)]
        vcf_outgroup = VCF(vcffile, gts012=True, samples=outsamples)
    vcf_query = VCF(vcffile, gts012=True, samples=querysamples)
    if len(querysamples) > len(vcf_query.samples):
        miss = set(querysamples) - set(vcf_query.samples)
        print(f'Query samples missing: {miss}')
        for ind in miss:
            querysamples.remove(ind)
    
    # Data storage for the heatmap
    df = []
    index = []
    
    if outlist:
        # Determine the ancestral allele using the outgroup
        for variant_outgroup, variant_query in zip(vcf_outgroup(region), vcf_query(region)):
            counts = np.bincount(variant_outgroup.gt_types)  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
            try:
                major_gt = np.argmax([counts[0], counts[2]])  # Compare counts of HOM_REF and HOM_ALT
            except IndexError:  # Handle cases with no HOM_ALT
                major_gt = 0
            arr = variant_query.gt_types
            if major_gt == 0:
                # If the major allele is ALT, switch REF and ALT
                arr[arr == 2] = -9
                arr[arr == 0] = 2
                arr[arr == -9] = 0
            df.append(arr.tolist())
            index.append(variant_query.POS)
    else:
        for variant_query in vcf_query(region):
            arr = variant_query.gt_types
            df.append(arr.tolist())
            index.append(variant_query.POS)

    # Create a DataFrame for heatmap plotting
    df = pd.DataFrame(df, columns=vcf_query.samples, index=index)
    df = df[querysamples]  # Reorder columns based on the query list
    df = df.replace(3, np.nan)  # Replace UNKNOWN values with NaN
    print(f'{os.path.basename(vcffile)} {region}:\n{df.shape}')

    # Rename columns to include group information
    id2groups = loadgroup(groupfile)
    df.columns = [f'{id2groups[x]}_{x}' for x in df.columns]

    # Filter by minor allele frequency (MAF)
    freqs = df.sum(axis=1).values / (df.count(axis=1).values * 2)
    df = df.loc[((1 - maf) >= freqs) & (freqs >= maf), :]
    print(f'Filtered by MAF({maf}):\n{df.shape}')

    # Plot the heatmap
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_facecolor("grey")
    sns.heatmap(df.T, yticklabels=1, cmap='OrRd', ax=ax)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(ticklabelsize)
    plt.savefig(outfile, dpi=dpi)
    plt.close()

if __name__ == '__main__':
    main()
