# NIA CARD Long Read Sequencing Methylation Benchmarking
Oxford Nanopore Technologies (ONT) long read sequencing data not only sequences longer fragments of nucleic acid than short read Illumina sequencing but also identifies modifications in the native nucleic acid. Much like sequencing accuracy, methylation calling accuracy has steadily improved with chemistry, sampling rate, and basecalling model improvements. We thus provide scripts and guidelines for standard benchmarking of methylation calling from Oxford Nanopore (ONT) long read sequencing data to evaluate changes in performance.
## Benchmarking within modkit
The [modkit](https://github.com/nanoporetech/modkit) methylation modification toolkit ONT provides includes a subcommand to benchmark methylation from 10,042 reads in an alignment (mapped BAM) by default.
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
