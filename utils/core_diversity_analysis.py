import argparse
import os
import pandas as pd
import unifrac

from skbio.stats.ordination import pcoa

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--biom', help='Path to feature table in biom format')
parser.add_argument('-t', '--tree', help='Path to phylogenetic tree in newick format')
parser.add_argument('-o', '--outdir', help='Path to directory for output files.')
args = parser.parse_args()

BIOM_FP = args.biom
TREE_FP = args.tree
OUTDIR = args.outdir

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

# unifrac distance matrix filepaths
UU_FP = os.path.join(OUTDIR, './unifrac_weighted_distance.tsv')
WU_FP = os.path.join(OUTDIR, './unifrac_unweighted_distance.tsv')

# pcoa coordinate filepaths
UU_PCOA_FP = os.path.join(OUTDIR, './unifrac_unweighted_pcoa.tsv')
WU_PCOA_FP = os.path.join(OUTDIR, './unifrac_weighted_pcoa.tsv')

# pcoa variance explained filepaths
UU_VAR_FP = os.path.join(OUTDIR, './unifrac_unweighted_variance.tsv')
WU_VAR_FP = os.path.join(OUTDIR, './unifrac_weighted_variance.tsv')

# faith's pd filepath
ALPHA_FP = os.path.join(OUTDIR, './faith_pd.tsv')

# main function
# calculates weighted and unweighted unifrac and faith pd
# spits out distance matrices, pcoa coordinates, proportion explained
# and alpha diversity metrics.
def calc_div():

    uu = unifrac.unweighted(table=BIOM_FP,
                            phylogeny=TREE_FP)

    wu = unifrac.weighted_unnormalized(table=BIOM_FP,
                                    phylogeny=TREE_FP)

    alpha = unifrac.faith_pd(biom_filename=BIOM_FP,
                            tree_filename=TREE_FP)

    # convert skbio distance matrix to pandas dataframe
    # and then write to file
    uu.to_data_frame().to_csv(UU_FP, sep='\t', index=False)
    wu.to_data_frame().to_csv(WU_FP, sep='\t', index=False)

    # calculate pcoa for weighted and unweighted 
    # unifrac distances
    uu_pcoa = pcoa(uu.data, number_of_dimensions=10)
    wu_pcoa = pcoa(wu.data, number_of_dimensions=10)

    # insert column at 0th index for sample id by
    # referencing the rownames (index) for the 
    # unifrac distances 
    uu_pcoa.samples.insert(loc=0, 
                       column='sampleid',
                       value=uu.ids)
    wu_pcoa.samples.insert(loc=0, 
                        column='sampleid',
                        value=wu.ids)
    
    # write the pcoa coordinates to a tsv
    uu_pcoa.samples.to_csv(UU_PCOA_FP, sep = '\t', index = False)
    wu_pcoa.samples.to_csv(WU_PCOA_FP, sep = '\t', index = False)

    # convert the pcoa variance explained into a dataframe
    uu_variance = pd.DataFrame(uu_pcoa.proportion_explained, columns=['variance'])
    wu_variance = pd.DataFrame(wu_pcoa.proportion_explained, columns=['variance'])

    # insert a column at the 0th index for the pcs
    # by referencing the rownames (index)
    uu_variance.insert(loc=0,
                       column='PC',
                       value=uu_variance.index)
    wu_variance.insert(loc=0,
                       column='PC',
                       value=wu_variance.index)

    # write the pcoa variance explained to a tsv
    uu_variance.to_csv(UU_VAR_FP, sep='\t', index=False)
    wu_variance.to_csv(WU_VAR_FP, sep='\t', index=False)

    alpha = pd.DataFrame(alpha, columns=['faith_pd'])
    alpha.insert(loc=0,
                column='sampleid',
                value=alpha.index)
    alpha.to_csv(ALPHA_FP, sep='\t', index=False)


if __name__ == '__main__':
    calc_div()