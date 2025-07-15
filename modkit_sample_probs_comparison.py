#!/usr/bin/python 

# modkit_sample_probs_comparison.py
# script to generate comparative line plots of methylation likelihood for multiple samples for benchmarking
# uses modkit sample-probs probabilities.tsv as input

import argparse
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import time
import psutil
import os
import re

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Compare modkit sample-probs methylation probabilities TSVs between samples by drawing counts or fractions vs. methylation likelihood for given base/modification. Uses all provided bases/modifications by default.")
    # argument for inputs - 1 or more inputs
    parser.add_argument("--input", required=True, nargs="*", help="Input methylation probabilities TSV file(s).")
    # argument for names - 0 or more names
    parser.add_argument("--names", required=False, nargs="+", help="Name(s) of input methylation probabilities TSV file(s).")
    # argument for output prefix
    parser.add_argument("--output_prefix", required=True, help="Prefix for output lineplots.")
    # argument for plot title
    parser.add_argument("--plot_title", required=False, default="Modkit sample-probs ML benchmark", help="Title for methylation likelihood lineplots.")
    # argument for plotting counts or fractions
    parser.add_argument("--dependent_variable",choices=["counts", "fractions"],default="counts",required=False,help="Dependent variable (y-axis) to use in lineplot (counts or fractions).")
    # argument for minimum methylation likelihood
    parser.add_argument("--min_ml", required=False, type=float, help="Minimum methylation likelihood to plot (between 0 and 1).")
    # argument for maximum methylation likelihood
    parser.add_argument("--max_ml", required=False, type=float, help="Maximum methylation likelihood to plot (between 0 and 1).")
    # return parsed arguments
    return parser.parse_args()
    
# find unique list of base modifications and relabel to modification type
def get_bases_modifications(meth_probabilities_df):
    # get unique set of 'code'/'primary_base' combinations
    base_mod_combos=meth_probabilities_df[['code','primary_base']].drop_duplicates(ignore_index=True)
    # relabel to modification type - e.g., 6mA, 5mC, 5hmC, etc.
    base_mod_combos.insert(2,'label',np.nan)
    # currently handles A, C, 6mA, 5mC, and 5hmC
    base_mod_combos.loc[(base_mod_combos['code']=='-') & (base_mod_combos['primary_base']=='A'),'label']='A'
    base_mod_combos.loc[(base_mod_combos['code']=='a') & (base_mod_combos['primary_base']=='A'),'label']='6mA'
    base_mod_combos.loc[(base_mod_combos['code']=='-') & (base_mod_combos['primary_base']=='C'),'label']='C'
    base_mod_combos.loc[(base_mod_combos['code']=='h') & (base_mod_combos['primary_base']=='C'),'label']='5hmC'
    base_mod_combos.loc[(base_mod_combos['code']=='m') & (base_mod_combos['primary_base']=='C'),'label']='5mC'
    # return base_mod_combos object
    return(base_mod_combos)

# make methylation likelihood plots for all base/modification combos    
def meth_likelihood_plot(base_mod_combos,concat_meth_probabilities_df,dependent_variable,plot_title,output_prefix,min_ml,max_ml):
    for idx in range(len(base_mod_combos)):
        # for debugging
        # print(idx)
        # get base mod subset df
        current_base_mod_subset_df=concat_meth_probabilities_df.loc[(concat_meth_probabilities_df['code']==base_mod_combos.loc[idx,'code']) & (concat_meth_probabilities_df['primary_base']==base_mod_combos.loc[idx,'primary_base'])]
        # for debugging
        # print(current_base_mod_subset_df)
        # filter by min and max ML
        # if (min_ml is not None):
        #     current_base_mod_subset_df=current_base_mod_subset_df[current_base_mod_subset_df['range_start']>min_ml]
        # if (max_ml is not None):
        #    current_base_mod_subset_df=current_base_mod_subset_df[current_base_mod_subset_df['range_start']<max_ml]
        # initialize figure
        fig, ax = plt.subplots()
        # set up lineplot - different configuration if plotting raw site counts or fractions
        if (dependent_variable == 'counts'):
            ax = sb.lineplot(x='range_start', y='count', hue='name', data=current_base_mod_subset_df)
            # use logarithmic y axis scale for counts
            ax.set_yscale('log')
            # set axis labels
            ax.set(xlabel=base_mod_combos.loc[idx,'label'] + " methylation likelihood",ylabel="Counts")
        elif (dependent_variable == 'fractions'):
            ax = sb.lineplot(x='range_start', y='frac', hue='name', data=current_base_mod_subset_df)
            # set axis labels
            ax.set(xlabel=base_mod_combos.loc[idx,'label'] + " methylation likelihood",ylabel="Fractions")
        # set x axis limits based on min and max ML (methylation likelihood) - tried doing ML filtering before
        ax.set_xlim(min_ml,max_ml)
        # set plot title
        ax.set(title=plot_title)
        # set legend title
        ax.get_legend().set_title("Input")
        # save figure - file name is (output_prefix)_(modification_name)_ML_lineplot.png
        fig.savefig(output_prefix + "_" + base_mod_combos.loc[idx,'label'] + "_ML_lineplot.png", format='png', dpi=300, bbox_inches='tight')
        # close figure
        fig.clf()
        
# main script subroutine
def main():
    # Parse the arguments
    args = parse_args()
    # throw error if no input file provided
    if args.input is None:
        quit('ERROR: No input file (-input) provided!')
    # throw error if no names provided if multiple input files provided
    if len(args.input)>1:
        if len(args.names)<=1:
            quit('ERROR: Multiple input files provided but not multiple names (-names).')
    # import methylation probabilities TSV files
    meth_probabilities_df_list=[0] * len(args.input)
    for idx, i in enumerate(args.input):
        meth_probabilities_df_list[idx]=pd.read_csv(i,sep='\t')
        # add name column to each imported data frame based on --names argument order
        meth_probabilities_df_list[idx]['name']=args.names[idx]
    # concatenate data frames into single data frame
    concat_meth_probabilities_df=pd.concat(meth_probabilities_df_list[:],ignore_index=True)
    # get base/modification list
    unique_base_mod_pairs=get_bases_modifications(concat_meth_probabilities_df)
    # make methylation likelihood plots
    meth_likelihood_plot(unique_base_mod_pairs,concat_meth_probabilities_df,args.dependent_variable,args.plot_title,args.output_prefix,args.min_ml,args.max_ml)
    
if __name__ == "__main__":
    main()