#!/usr/bin/python 

# CARDlongread_methylation_entropy_pairwise_comparison.py
# script to compare methylation entropies between two samples and further analyze relationships with respect to methylation differences and supporting coverage
# uses modkit entropy, modkit dmr, and DSS/bsseq DMR tab-delimited outputs as input for analysis and visualization

import argparse
import pandas as pd
# import polars as pl
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import time
import psutil
import os
import re

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Compare methylation entropies between two ONT sequenced samples and further analyze relationships with respect to methylation differences and supporting coverage.")
    # arguments for sample names
    parser.add_argument("--sample_name_1", required=True, help="Name of first sample.")
    parser.add_argument("--sample_name_2", required=True, help="Name of second sample.")
    # argument for entropy inputs
    # overall entropy
    parser.add_argument("--sample_1_bulk_entropy", required=True, help="Default ONT entropy input for sample 1 (e.g., 50 bp windows).")
    parser.add_argument("--sample_2_bulk_entropy", required=True, help="Default ONT entropy input for sample 2 (e.g., 50 bp windows).")
    # modkit dmr entropy 
    parser.add_argument("--sample_1_modkit_dmr_entropy", required=True, help="Modkit DMR ONT entropy input for sample 1 (entropy per DMR).")
    parser.add_argument("--sample_2_modkit_dmr_entropy", required=True, help="Modkit DMR ONT entropy input for sample 2 (entropy per DMR).")
    # DSS/bsseq unsmoothed dmr entropy
    parser.add_argument("--sample_1_dss_unsmoothed_dmr_entropy", required=True, help="Unsmoothed DSS/bsseq DMR ONT entropy input for sample 1 (entropy per DMR).")
    parser.add_argument("--sample_2_dss_unsmoothed_dmr_entropy", required=True, help="Unsmoothed DSS/bsseq DMR ONT entropy input for sample 2 (entropy per DMR).")
    # DSS/bsseq smoothed dmr entropy
    parser.add_argument("--sample_1_dss_smoothed_dmr_entropy", required=True, help="Smoothed DSS/bsseq DMR ONT entropy input for sample 1 (entropy per DMR).")
    parser.add_argument("--sample_2_dss_smoothed_dmr_entropy", required=True, help="Smoothed DSS/bsseq DMR ONT entropy input for sample 2 (entropy per DMR).")
    # modkit dmrs
    parser.add_argument("--modkit_dmr_segments", required=True, help="Modkit DMR segments for sample 1 vs. sample 2 (modkit dmr pair output).")
    # DSS/bsseq unsmoothed dmrs
    parser.add_argument("--dss_unsmoothed_dmrs", required=True, help="Unsmoothed DSS/bsseq DMRs for sample 1 vs. sample 2 (bsseq callDMR output on unsmoothed methylation inputs).")
    # DSS/bsseq smoothed dmrs
    parser.add_argument("--dss_smoothed_dmrs", required=True, help="Smoothed DSS/bsseq DMRs for sample 1 vs. sample 2 (bsseq callDMR output on smoothed inputs).")
    # argument for output prefix
    parser.add_argument("--output_prefix", required=True, help="Prefix for output lineplots.")
    # argument for plot title
    parser.add_argument("--plot_title", required=False, default="ONT pairwise entropy comparison", help="Title for each output plot.")
    # argument for read count cutoff
    parser.add_argument("--read_count_cutoff", required=False, default=500, help="Read count cutoff for read count vs. methylation entropy plot.")
    # return parsed arguments
    return parser.parse_args()

# subroutine to make per sample entropy distribution comparison histograms
# include bulk, modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_entropy_distribution_histogram(per_sample_entropy_df,sample_name,plot_title,output_prefix):
    # rename "name" column to "Region type"
    per_sample_entropy_df=per_sample_entropy_df.rename(columns={'name': 'Region type'})
    # initialize figure
    fig, ax = plt.subplots()
    # plot per sample entropy distribution histogram with dodge set to true (side by side bars per data type)
    ax = sb.histplot(data=per_sample_entropy_df,x="mean_entropy",hue="Region type",multiple="dodge",stat="proportion",common_norm=False,kde=True)
    # label x-axis
    ax.set(xlabel="Methylation entropy per region")
    # label y-axis
    ax.set(ylabel="Proportion of regions")
    # label legend
    # plt.legend(title="Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_entropy_distribution_histogram.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_entropy_distribution_histogram.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    
# subroutine to plot sample 1 vs. sample 2 pairwise entropy on scatterplot
# include bulk, modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def pairwise_entropy_scatterplot(both_samples_entropy_df,sample_name_1,sample_name_2,plot_title,output_prefix):
    # initialize figure
    fig, ax = plt.subplots()
    # plot scatterplot of sample 1 entropies against sample 2 entropies for common regions in each sample
    ax = sb.scatterplot(data=both_samples_entropy_df,x="mean_entropy_x",y="mean_entropy_y",hue="common_name")
    # label axes
    ax.set(xlabel=sample_name_1 + " methylation entropy per region")
    ax.set(ylabel=sample_name_2 + " methylation entropy per region")
    # label legend
    legend = ax.legend()
    legend.set_title("Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name_1)_v_(sample_name_2)_pairwise_entropy_scatterplot.png
    fig.savefig(output_prefix + "_" + sample_name_1 + "_v_" + sample_name_2 + "_pairwise_entropy_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    
# subroutine to make per sample entropy vs. supporting read count/read proportion scatterplots
# include bulk, modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_entropy_read_count_scatterplot(per_sample_entropy_df,sample_name,plot_title,output_prefix,read_count_cutoff):
    # initialize figure
    fig, ax = plt.subplots()
    # plot scatterplot of entropies vs. supporting read counts for all entropy types per sample
    ax = sb.scatterplot(data=per_sample_entropy_df,x="mean_num_reads",y="mean_entropy",hue="name")
    # label axes
    ax.set(xlabel="Number of reads supporting entropy call")
    ax.set(ylabel="Methylation entropy per region")
    # set x-limit
    ax.set_xlim(left=0,right=read_count_cutoff)
    # label legend
    legend = ax.legend()
    legend.set_title("Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_entropy_read_proportion_scatterplot.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_entropy_read_count_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()

# subroutine to make per sample entropy vs. meth changes plots
# include modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_entropy_methylation_changes_scatterplot(per_sample_entropy_methylation_df,sample_name,plot_title,output_prefix):
    # initialize figure
    fig, ax = plt.subplots()
    # plot entropy vs. respective methylation changes per DMR per sample
    ax = sb.scatterplot(data=per_sample_entropy_methylation_df,x="mean_entropy",y="effect_size",hue="name")
    # label axes
    ax.set(xlabel="Methylation entropy per DMR")
    ax.set(ylabel="Methylation change per DMR")
    # label legend
    legend = ax.legend()
    legend.set_title("Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_entropy_methylation_changes_scatterplot.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_entropy_methylation_changes_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()

# subroutine to make per sample entropy vs. DMR length plots
# include modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_entropy_DMR_length_scatterplot(per_sample_entropy_methylation_df,sample_name,plot_title,output_prefix):
    # initialize figure
    fig, ax = plt.subplots()
    # plot entropy vs. respective methylation changes per DMR per sample
    ax = sb.scatterplot(data=per_sample_entropy_methylation_df,x="DMR length",y="mean_entropy",hue="name")
    # label axes
    ax.set(xlabel="DMR length (bp)")
    ax.set(ylabel="Methylation entropy per DMR")
    # label legend
    legend = ax.legend()
    legend.set_title("Region type")
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_entropy_DMR_length_scatterplot.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_entropy_DMR_length_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()

# subroutine to make DMR change vs. DMR length plots
# include modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_DMR_change_DMR_length_scatterplot(per_sample_entropy_methylation_df,sample_name,plot_title,output_prefix):
    # initialize figure
    fig, ax = plt.subplots()
    # plot entropy vs. respective methylation changes per DMR per sample
    ax = sb.scatterplot(data=per_sample_entropy_methylation_df,x="DMR length",y="effect_size",hue="name")
    # label axes
    ax.set(xlabel="DMR length (bp)")
    ax.set(ylabel="Methylation change per DMR")
    # label legend
    legend = ax.legend()
    legend.set_title("Region type")
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_DMR_change_DMR_length_scatterplot.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_DMR_change_DMR_length_scatterplot.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()

# subroutine to make DMR length histogram
# include modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_DMR_length_distribution_histogram(per_sample_entropy_methylation_df,sample_name,plot_title,output_prefix):
    # rename "name" column to "Region type"
    per_sample_entropy_methylation_df=per_sample_entropy_methylation_df.rename(columns={'name': 'Region type'})
    # initialize figure
    fig, ax = plt.subplots()
    # plot per sample entropy distribution histogram with dodge set to true (side by side bars per data type)
    ax = sb.histplot(data=per_sample_entropy_methylation_df,x="DMR length",hue="Region type",multiple="dodge",stat="proportion",common_norm=False,kde=True)
    # label x-axis
    ax.set(xlabel="DMR length (bp)")
    # label y-axis
    ax.set(ylabel="Proportion of DMRs")
    # label legend
    # plt.legend(title="Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_DMR_length_distribution_histogram.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_DMR_length_distribution_histogram.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()

# subroutine to make DMR change histogram
# include modkit dmr pair, unsmoothed DSS DMR, and smoothed DSS DMR on one plot
def per_sample_DMR_change_distribution_histogram(per_sample_entropy_methylation_df,sample_name,plot_title,output_prefix):
    # rename "name" column to "Region type"
    per_sample_entropy_methylation_df=per_sample_entropy_methylation_df.rename(columns={'name': 'Region type'})
    # initialize figure
    fig, ax = plt.subplots()
    # plot per sample entropy distribution histogram with dodge set to true (side by side bars per data type)
    ax = sb.histplot(data=per_sample_entropy_methylation_df,x="effect_size",hue="Region type",multiple="dodge",stat="proportion",common_norm=False,kde=True)
    # label x-axis
    ax.set(xlabel="Methylation change per DMR")
    # label y-axis
    ax.set(ylabel="Proportion of DMRs")
    # label legend
    # plt.legend(title="Region type")
    # add title
    ax.set_title(plot_title)
    # save figure - file name is (output_prefix)_(sample_name)_per_sample_DMR_change_distribution_histogram.png
    fig.savefig(output_prefix + "_" + sample_name + "_per_sample_DMR_change_distribution_histogram.png", format='png', dpi=300, bbox_inches='tight')
    # close figure
    fig.clf()
    
# main script subroutine
def main():
    # Parse the arguments
    args = parse_args()
    # load entropy files
    sample_1_bulk_entropy_df=pd.read_csv(args.sample_1_bulk_entropy,sep="\t",header=None)
    sample_2_bulk_entropy_df=pd.read_csv(args.sample_2_bulk_entropy,sep="\t",header=None)
    sample_1_modkit_dmr_entropy_df=pd.read_csv(args.sample_1_modkit_dmr_entropy,sep="\t",header=None)
    sample_2_modkit_dmr_entropy_df=pd.read_csv(args.sample_2_modkit_dmr_entropy,sep="\t",header=None)
    sample_1_dss_unsmoothed_dmr_entropy_df=pd.read_csv(args.sample_1_dss_unsmoothed_dmr_entropy,sep="\t",header=None)
    sample_2_dss_unsmoothed_dmr_entropy_df=pd.read_csv(args.sample_2_dss_unsmoothed_dmr_entropy,sep="\t",header=None)
    sample_1_dss_smoothed_dmr_entropy_df=pd.read_csv(args.sample_1_dss_smoothed_dmr_entropy,sep="\t",header=None)
    sample_2_dss_smoothed_dmr_entropy_df=pd.read_csv(args.sample_2_dss_smoothed_dmr_entropy,sep="\t",header=None)
    # add headers (column names) where necessary
    sample_1_bulk_entropy_df.columns=['chrom','start','end','entropy','strand','num_reads']
    sample_2_bulk_entropy_df.columns=['chrom','start','end','entropy','strand','num_reads']
    sample_1_modkit_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    sample_2_modkit_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    sample_1_dss_unsmoothed_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    sample_1_dss_smoothed_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    sample_2_dss_unsmoothed_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    sample_2_dss_smoothed_dmr_entropy_df.columns=['chrom','start','end','region_name','mean_entropy','strand','median_entropy','min_entropy','max_entropy','mean_num_reads','min_num_reads','max_num_reads','successful_window_count','failed_window_count']
    # load DMRs
    modkit_dmr_segments_df=pd.read_csv(args.modkit_dmr_segments,sep="\t",header=None)
    dss_unsmoothed_dmr_df=pd.read_csv(args.dss_unsmoothed_dmrs,sep="\t")
    # rename columns
    dss_unsmoothed_dmr_df=dss_unsmoothed_dmr_df.rename(columns={'chr': 'chrom'})
    dss_smoothed_dmr_df=pd.read_csv(args.dss_smoothed_dmrs,sep="\t")
    # rename columns
    dss_smoothed_dmr_df=dss_smoothed_dmr_df.rename(columns={'chr': 'chrom'})
    # add headers (column names) where necessary
    modkit_dmr_segments_df.columns=['chrom','start','end','state-name','score','N-sites','sample_a_counts','sample_b_counts','sample_a_percents','sample_b_percents','sample_a_fraction_modified','sample_b_fraction_modified','effect_size','cohen_h','cohen_h_low','cohen_h_high']
    # add region name column for table join
    # note that you need to first convert start/end positions to integer, and then to string...note trailing .0 in smoothed dmr values
    modkit_dmr_segments_df['region_name']=modkit_dmr_segments_df['chrom'].astype(str)+":"+modkit_dmr_segments_df['start'].astype(int).astype(str)+"-"+modkit_dmr_segments_df['end'].astype(int).astype(str)
    dss_unsmoothed_dmr_df['region_name']=dss_unsmoothed_dmr_df['chrom'].astype(str)+":"+dss_unsmoothed_dmr_df['start'].astype(int).astype(str)+"-"+dss_unsmoothed_dmr_df['end'].astype(int).astype(str)
    dss_smoothed_dmr_df['region_name']=dss_smoothed_dmr_df['chrom'].astype(str)+":"+dss_smoothed_dmr_df['start'].astype(int).astype(str)+"-"+dss_smoothed_dmr_df['end'].astype(int).astype(str)
    # get dmr segments tagged different
    # modkit_dmr_segments_df=modkit_dmr_segments_df[modkit_dmr_segments_df['state-name']=='different']
    # combine entropies and DMRs appropriately
    # all sample 1
    sample_1_modkit_dmr_segments_entropy_df=pd.merge(sample_1_modkit_dmr_entropy_df,modkit_dmr_segments_df,on=['region_name'])
    sample_1_dss_unsmoothed_dmr_entropy_dmrs_df=pd.merge(sample_1_dss_unsmoothed_dmr_entropy_df,dss_unsmoothed_dmr_df,on=['region_name'])
    sample_1_dss_smoothed_dmr_entropy_dmrs_df=pd.merge(sample_1_dss_smoothed_dmr_entropy_df,dss_smoothed_dmr_df,on=['region_name'])
    # filter for modkit dmr segments that are truly different in methylation
    sample_1_modkit_dmr_segments_entropy_df=sample_1_modkit_dmr_segments_entropy_df[sample_1_modkit_dmr_segments_entropy_df['state-name']=='different']
    # get just necessary columns for plotting
    sample_1_modkit_dmr_segments_entropy_diff_df=sample_1_modkit_dmr_segments_entropy_df[['region_name','mean_entropy','N-sites','effect_size']].rename(columns={'N-sites': 'DMR length'})
    sample_1_dss_unsmoothed_dmr_entropy_dmrs_diff_df=sample_1_dss_unsmoothed_dmr_entropy_dmrs_df[['region_name','mean_entropy','length','diff.Methy']].rename(columns={'diff.Methy': 'effect_size', 'length': 'DMR length'})
    sample_1_dss_smoothed_dmr_entropy_dmrs_diff_df=sample_1_dss_smoothed_dmr_entropy_dmrs_df[['region_name','mean_entropy','length','diff.Methy']].rename(columns={'diff.Methy': 'effect_size', 'length': 'DMR length'})
    # name by region type
    sample_1_modkit_dmr_segments_entropy_diff_df.insert(loc=0, column='name', value=args.sample_name_1 + " modkit DMR segments")
    sample_1_dss_unsmoothed_dmr_entropy_dmrs_diff_df.insert(loc=0, column='name', value=args.sample_name_1 + " DSS unsmoothed DMRs")
    sample_1_dss_smoothed_dmr_entropy_dmrs_diff_df.insert(loc=0, column='name', value=args.sample_name_1 + " DSS smoothed DMRs")
    # concatenate sample 1 entropy/dmr tables
    sample_1_concat_dmr_entropy_table=pd.concat([sample_1_modkit_dmr_segments_entropy_diff_df,sample_1_dss_unsmoothed_dmr_entropy_dmrs_diff_df,sample_1_dss_smoothed_dmr_entropy_dmrs_diff_df],ignore_index=True)
    # debugging
    # print(sample_1_concat_dmr_entropy_table)
    # get name column data type
    # print(sample_1_concat_dmr_entropy_table['name'].dtype)
    # set name column to string data type
    sample_1_concat_dmr_entropy_table['name'] = sample_1_concat_dmr_entropy_table['name'].astype('string')
    # sample 1 just entropies after join
    sample_1_bulk_entropy_renamed_df=sample_1_bulk_entropy_df.rename(columns={'chrom': 'chrom_x', 'start': 'start_x', 'end': 'end_x', 'entropy': 'mean_entropy', 'num_reads': 'mean_num_reads'}).drop(columns=['strand'])
    sample_1_modkit_dmr_segments_entropy_only_df=sample_1_modkit_dmr_segments_entropy_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    sample_1_dss_unsmoothed_dmr_entropy_only_dmrs_df=sample_1_dss_unsmoothed_dmr_entropy_dmrs_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    sample_1_dss_smoothed_dmr_entropy_only_dmrs_df=sample_1_dss_smoothed_dmr_entropy_dmrs_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    # name entropy by region type
    sample_1_bulk_entropy_renamed_df.insert(loc=0, column='name', value=args.sample_name_1 + " genomic windows")
    # sample_1_bulk_entropy_renamed_df['name']=args.sample_name_1 + " genomic windows"
    sample_1_modkit_dmr_segments_entropy_only_df.insert(loc=0, column='name', value=args.sample_name_1 + " modkit DMR segments")
    # sample_1_modkit_dmr_segments_entropy_only_df.loc[:,'name']=args.sample_name_1 + " modkit DMR segments"
    sample_1_dss_unsmoothed_dmr_entropy_only_dmrs_df.insert(loc=0, column='name', value=args.sample_name_1 + " DSS unsmoothed DMRs")
    # sample_1_dss_unsmoothed_dmr_entropy_only_dmrs_df.loc[:,'name']=args.sample_name_1 + " DSS unsmoothed DMRs"
    sample_1_dss_smoothed_dmr_entropy_only_dmrs_df.insert(loc=0, column='name', value=args.sample_name_1 + " DSS smoothed DMRs")
    # sample_1_dss_smoothed_dmr_entropy_only_dmrs_df.loc[:,'name']=args.sample_name_1 + " DSS smoothed DMRs"
    # concat above tables
    sample_1_concat_entropy_table=pd.concat([sample_1_bulk_entropy_renamed_df,sample_1_modkit_dmr_segments_entropy_only_df,sample_1_dss_unsmoothed_dmr_entropy_only_dmrs_df,sample_1_dss_smoothed_dmr_entropy_only_dmrs_df],ignore_index=True)
    # set name column to string data type
    sample_1_concat_entropy_table['name'] = sample_1_concat_entropy_table['name'].astype('string')
    # all sample 2
    sample_2_modkit_dmr_segments_entropy_df=pd.merge(sample_2_modkit_dmr_entropy_df,modkit_dmr_segments_df,on=['region_name'])
    sample_2_dss_unsmoothed_dmr_entropy_dmrs_df=pd.merge(sample_2_dss_unsmoothed_dmr_entropy_df,dss_unsmoothed_dmr_df,on=['region_name'])
    sample_2_dss_smoothed_dmr_entropy_dmrs_df=pd.merge(sample_2_dss_smoothed_dmr_entropy_df,dss_smoothed_dmr_df,on=['region_name'])
    # filter for modkit dmr segments that are truly different in methylation
    sample_2_modkit_dmr_segments_entropy_df=sample_2_modkit_dmr_segments_entropy_df[sample_2_modkit_dmr_segments_entropy_df['state-name']=='different']
    # get just necessary columns for plotting
    sample_2_modkit_dmr_segments_entropy_diff_df=sample_2_modkit_dmr_segments_entropy_df[['region_name','mean_entropy','N-sites','effect_size']].rename(columns={'N-sites': 'DMR length'})
    sample_2_dss_unsmoothed_dmr_entropy_dmrs_diff_df=sample_2_dss_unsmoothed_dmr_entropy_dmrs_df[['region_name','mean_entropy','length','diff.Methy']].rename(columns={'diff.Methy': 'effect_size', 'length': 'DMR length'})
    sample_2_dss_smoothed_dmr_entropy_dmrs_diff_df=sample_2_dss_smoothed_dmr_entropy_dmrs_df[['region_name','mean_entropy','length','diff.Methy']].rename(columns={'diff.Methy': 'effect_size', 'length': 'DMR length'})
    # name by region type
    sample_2_modkit_dmr_segments_entropy_diff_df.insert(loc=0, column='name', value=args.sample_name_2 + " modkit DMR segments")
    sample_2_dss_unsmoothed_dmr_entropy_dmrs_diff_df.insert(loc=0, column='name', value=args.sample_name_2 + " DSS unsmoothed DMRs")
    sample_2_dss_smoothed_dmr_entropy_dmrs_diff_df.insert(loc=0, column='name', value=args.sample_name_2 + " DSS smoothed DMRs")
    # concatenate sample 1 entropy/dmr tables
    sample_2_concat_dmr_entropy_table=pd.concat([sample_2_modkit_dmr_segments_entropy_diff_df,sample_2_dss_unsmoothed_dmr_entropy_dmrs_diff_df,sample_2_dss_smoothed_dmr_entropy_dmrs_diff_df],ignore_index=True)
    # set name column to string data type
    sample_2_concat_dmr_entropy_table['name'] = sample_2_concat_dmr_entropy_table['name'].astype('string')
    # sample 2 just entropies after join
    sample_2_bulk_entropy_renamed_df=sample_2_bulk_entropy_df.rename(columns={'chrom': 'chrom_x', 'start': 'start_x', 'end': 'end_x', 'entropy': 'mean_entropy', 'num_reads': 'mean_num_reads'}).drop(columns=['strand'])
    sample_2_modkit_dmr_segments_entropy_only_df=sample_2_modkit_dmr_segments_entropy_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    sample_2_dss_unsmoothed_dmr_entropy_only_dmrs_df=sample_2_dss_unsmoothed_dmr_entropy_dmrs_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    sample_2_dss_smoothed_dmr_entropy_only_dmrs_df=sample_2_dss_smoothed_dmr_entropy_dmrs_df[['chrom_x','start_x','end_x','mean_entropy','mean_num_reads']]
    # name entropy by region type
    sample_2_bulk_entropy_renamed_df.insert(loc=0, column='name', value=args.sample_name_2 + " genomic windows")
    sample_2_modkit_dmr_segments_entropy_only_df.insert(loc=0, column='name', value=args.sample_name_2 + " modkit DMR segments")
    sample_2_dss_unsmoothed_dmr_entropy_only_dmrs_df.insert(loc=0, column='name', value=args.sample_name_2 + " DSS unsmoothed DMRs")
    sample_2_dss_smoothed_dmr_entropy_only_dmrs_df.insert(loc=0, column='name', value=args.sample_name_2 + " DSS smoothed DMRs")
    # concat above tables
    sample_2_concat_entropy_table=pd.concat([sample_2_bulk_entropy_renamed_df,sample_2_modkit_dmr_segments_entropy_only_df,sample_2_dss_unsmoothed_dmr_entropy_only_dmrs_df,sample_2_dss_smoothed_dmr_entropy_only_dmrs_df],ignore_index=True)
    # set name column to string data type
    sample_2_concat_entropy_table['name'] = sample_2_concat_entropy_table['name'].astype('string')
    # sample 1 + 2 combined
    samples_1_and_2_combined_concat_entropies=pd.merge(sample_1_concat_entropy_table,sample_2_concat_entropy_table,on=['chrom_x','start_x','end_x'])
    samples_1_and_2_combined_concat_entropies['common_name']=samples_1_and_2_combined_concat_entropies['name_x']
    # remove sample name 1 from common name column.
    pattern = re.escape(args.sample_name_1) + r'\s'
    samples_1_and_2_combined_concat_entropies['common_name']=samples_1_and_2_combined_concat_entropies['common_name'].str.replace(pattern, '', regex=True, n=1)
    samples_1_and_2_combined_concat_entropies['common_name']=samples_1_and_2_combined_concat_entropies['common_name'].str.replace('genomic windows', 'Genomic windows', regex=True, n=1)
    # make 'common_name' data type string
    samples_1_and_2_combined_concat_entropies['common_name']=samples_1_and_2_combined_concat_entropies['common_name'].astype('string')
    # proceed with plotting
    # concatenate together in single labeled data frame
    # make each output plot
    # entropy histogram for sample 1
    per_sample_entropy_distribution_histogram(sample_1_concat_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # entropy histogram for sample 2
    per_sample_entropy_distribution_histogram(sample_2_concat_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
    # pairwise entropy scatterplot for both samples
    pairwise_entropy_scatterplot(samples_1_and_2_combined_concat_entropies,args.sample_name_1,args.sample_name_2,args.plot_title,args.output_prefix)
    # entropy vs. read count scatterplot for sample 1
    per_sample_entropy_read_count_scatterplot(sample_1_concat_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix,args.read_count_cutoff)
    # entropy vs. read count scatterplot for sample 2
    per_sample_entropy_read_count_scatterplot(sample_2_concat_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix,args.read_count_cutoff)
    # entropy vs. methylation change scatterplot for sample 1
    per_sample_entropy_methylation_changes_scatterplot(sample_1_concat_dmr_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # entropy vs. methylation change scatterplot for sample 2
    per_sample_entropy_methylation_changes_scatterplot(sample_2_concat_dmr_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
    # entropy vs. dmr length for sample 1
    per_sample_entropy_DMR_length_scatterplot(sample_1_concat_dmr_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # entropy vs. dmr length for sample 2
    per_sample_entropy_DMR_length_scatterplot(sample_2_concat_dmr_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
    # DMR change vs. DMR length for sample 1
    per_sample_DMR_change_DMR_length_scatterplot(sample_1_concat_dmr_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # DMR change vs. DMR length for sample 2
    per_sample_DMR_change_DMR_length_scatterplot(sample_2_concat_dmr_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
    # DMR length histogram for sample 1
    per_sample_DMR_length_distribution_histogram(sample_1_concat_dmr_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # DMR length histogram for sample 2
    per_sample_DMR_length_distribution_histogram(sample_2_concat_dmr_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
    # DMR change histogram for sample 1
    per_sample_DMR_change_distribution_histogram(sample_1_concat_dmr_entropy_table,args.sample_name_1,args.plot_title,args.output_prefix)
    # DMR change histogram for sample 2
    per_sample_DMR_change_distribution_histogram(sample_2_concat_dmr_entropy_table,args.sample_name_2,args.plot_title,args.output_prefix)
if __name__ == "__main__":
    main()