# NIA CARD Long Read Sequencing Methylation Benchmarking
Oxford Nanopore Technologies (ONT) long read sequencing data not only sequences longer fragments of nucleic acid than short read Illumina sequencing but also identifies modifications in the native nucleic acid. Much like sequencing accuracy, methylation calling accuracy has steadily improved with chemistry, sampling rate, and basecalling model improvements. We thus provide scripts and guidelines for standard benchmarking of methylation calling from Oxford Nanopore (ONT) long read sequencing data to evaluate changes in performance.
## Benchmarking within modkit
The [modkit](https://github.com/nanoporetech/modkit) methylation modification toolkit ONT provides includes a subcommand (sample-probs) to benchmark methylation from 10,042 reads in an alignment (mapped BAM) by default.

We have run ```modkit sample-probs``` to generate methylation probabilities TSV files like so to generate both methylation likelihood thresholds for ```modkit pileup``` and to generate histograms (```--hist```) of methylation likelihoods for each detected modification:
```
modkit sample-probs --force --hist -t 64 sample.bam -o MODKIT/sample --prefix sample
```
## Usage
We have included a script to compare modkit sample-probs methylation benchmarks across different samples by drawing methylation likelihood lineplots.
```
usage: modkit_sample_probs_comparison.py [-h] --input [INPUT ...] [--names NAMES [NAMES ...]] --output_prefix OUTPUT_PREFIX [--plot_title PLOT_TITLE] [--dependent_variable {counts,fractions}] [--min_ml MIN_ML] [--max_ml MAX_ML]

Compare modkit sample-probs methylation probabilities TSVs between samples by drawing counts or fractions vs. methylation likelihood for given base/modification. Uses all provided bases/modifications by default.

optional arguments:
  -h, --help            show this help message and exit
  --input [INPUT ...]   Input methylation probabilities TSV file(s).
  --names NAMES [NAMES ...]
                        Name(s) of input methylation probabilities TSV file(s).
  --output_prefix OUTPUT_PREFIX
                        Prefix for output lineplots.
  --plot_title PLOT_TITLE
                        Title for methylation likelihood lineplots.
  --dependent_variable {counts,fractions}
                        Dependent variable (y-axis) to use in lineplot (counts or fractions).
  --min_ml MIN_ML       Minimum methylation likelihood to plot (between 0 and 1).
  --max_ml MAX_ML       Maximum methylation likelihood to plot (between 0 and 1).
```
We also added a script to visualize pairwise comparisons between methylation entropies of two samples for which differentially methylated regions (DMRs) were identified with several different methods (modkit dmr pair fine-grained, DSS/bsseq with smoothing, and DSS/bsseq without smoothing). This script was developed as part of the BCM HGSC hackathon [MethSmoothEval](https://github.com/collaborativebioinformatics/MethSmoothEval) project.
```
usage: CARDlongread_methylation_entropy_pairwise_comparison.py [-h] --sample_name_1 SAMPLE_NAME_1 --sample_name_2 SAMPLE_NAME_2 --sample_1_bulk_entropy SAMPLE_1_BULK_ENTROPY --sample_2_bulk_entropy SAMPLE_2_BULK_ENTROPY
                                                               --sample_1_modkit_dmr_entropy SAMPLE_1_MODKIT_DMR_ENTROPY --sample_2_modkit_dmr_entropy SAMPLE_2_MODKIT_DMR_ENTROPY --sample_1_dss_unsmoothed_dmr_entropy
                                                               SAMPLE_1_DSS_UNSMOOTHED_DMR_ENTROPY --sample_2_dss_unsmoothed_dmr_entropy SAMPLE_2_DSS_UNSMOOTHED_DMR_ENTROPY --sample_1_dss_smoothed_dmr_entropy SAMPLE_1_DSS_SMOOTHED_DMR_ENTROPY
                                                               --sample_2_dss_smoothed_dmr_entropy SAMPLE_2_DSS_SMOOTHED_DMR_ENTROPY --modkit_dmr_segments MODKIT_DMR_SEGMENTS --dss_unsmoothed_dmrs DSS_UNSMOOTHED_DMRS --dss_smoothed_dmrs
                                                               DSS_SMOOTHED_DMRS --output_prefix OUTPUT_PREFIX [--plot_title PLOT_TITLE] [--read_count_cutoff READ_COUNT_CUTOFF]

Compare methylation entropies between two ONT sequenced samples and further analyze relationships with respect to methylation differences and supporting coverage.

optional arguments:
  -h, --help            show this help message and exit
  --sample_name_1 SAMPLE_NAME_1
                        Name of first sample.
  --sample_name_2 SAMPLE_NAME_2
                        Name of second sample.
  --sample_1_bulk_entropy SAMPLE_1_BULK_ENTROPY
                        Default ONT entropy input for sample 1 (e.g., 50 bp windows).
  --sample_2_bulk_entropy SAMPLE_2_BULK_ENTROPY
                        Default ONT entropy input for sample 2 (e.g., 50 bp windows).
  --sample_1_modkit_dmr_entropy SAMPLE_1_MODKIT_DMR_ENTROPY
                        Modkit DMR ONT entropy input for sample 1 (entropy per DMR).
  --sample_2_modkit_dmr_entropy SAMPLE_2_MODKIT_DMR_ENTROPY
                        Modkit DMR ONT entropy input for sample 2 (entropy per DMR).
  --sample_1_dss_unsmoothed_dmr_entropy SAMPLE_1_DSS_UNSMOOTHED_DMR_ENTROPY
                        Unsmoothed DSS/bsseq DMR ONT entropy input for sample 1 (entropy per DMR).
  --sample_2_dss_unsmoothed_dmr_entropy SAMPLE_2_DSS_UNSMOOTHED_DMR_ENTROPY
                        Unsmoothed DSS/bsseq DMR ONT entropy input for sample 2 (entropy per DMR).
  --sample_1_dss_smoothed_dmr_entropy SAMPLE_1_DSS_SMOOTHED_DMR_ENTROPY
                        Smoothed DSS/bsseq DMR ONT entropy input for sample 1 (entropy per DMR).
  --sample_2_dss_smoothed_dmr_entropy SAMPLE_2_DSS_SMOOTHED_DMR_ENTROPY
                        Smoothed DSS/bsseq DMR ONT entropy input for sample 2 (entropy per DMR).
  --modkit_dmr_segments MODKIT_DMR_SEGMENTS
                        Modkit DMR segments for sample 1 vs. sample 2 (modkit dmr pair output).
  --dss_unsmoothed_dmrs DSS_UNSMOOTHED_DMRS
                        Unsmoothed DSS/bsseq DMRs for sample 1 vs. sample 2 (bsseq callDMR output on unsmoothed methylation inputs).
  --dss_smoothed_dmrs DSS_SMOOTHED_DMRS
                        Smoothed DSS/bsseq DMRs for sample 1 vs. sample 2 (bsseq callDMR output on smoothed inputs).
  --output_prefix OUTPUT_PREFIX
                        Prefix for output lineplots.
  --plot_title PLOT_TITLE
                        Title for each output plot.
  --read_count_cutoff READ_COUNT_CUTOFF
                        Read count cutoff for read count vs. methylation entropy plot.
```
## Example outputs
<img width="1705" height="1357" alt="image" src="https://github.com/user-attachments/assets/143227de-f269-4e7e-a4c5-ca3598708cec" />
