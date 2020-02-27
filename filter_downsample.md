# Filtering and downsampling in bamQC

## Overview

This document describes filtering and downsampling which may be applied to BAM files before bamQC metric computation.

## Filtering

We remove the following undesirable read types:
- Non-primary alignments: Important to avoid double-counting (or multi-counting) of reads. These are defined as having samtools flag 256 (secondary alignment) or 2048 (supplementary alignment).
- Unmapped reads. Defined as having samtools flag 4 (segment unmapped).
- Low-quality alignments: Reads with a low alignment quality score (default minimum = 30)

Totals of non-primary, unmapped, and low-quality reads are recorded in workflow output. The excluded reads are not used for subsequent QC metric computation.

## Downsampling

Very large BAM files are downsampled in order to enable faster QC. Downsampling consists of choosing a subset of reads at random. The downsampled file is used as input to the "slow" computationally intensive metrics. "Fast" metrics do not require downsampling, and are always computed on the full-sized filtered input.

Picard MarkDuplicates is considered a "slow" metric. For other metrics, see the [bam-qc-metrics documentation](https://github.com/oicr-gsi/bam-qc-metrics/blob/master/metrics.md).

### Threshold

If a BAM file has more reads than a threshold T, downsample to a subset containing T reads.

The default value of T is 100000. For example:
- A BAM file with 17000000 reads will be downsampled to 100000.
- A BAM file with 100001 reads will be downsampled to 100000.
- A BAM file with 99999 reads will not be downsampled.

So, if the original number of reads is 100000 or greater, downsampling is guaranteed not to reduce it below 100000.

Downsampling thresholds are always evaluated on the BAM file _after_ filtering.

### Methods

There are two ways of downsampling with samtools.

#### Random downsampling

This is done using `samtools view -s $FOO.$BAR`, where `$FOO` is a random seed and `$BAR` is a decimal expressing the probability of retaining a read. For example, to sample reads with a probability of 0.05 and random seed 99, use `samtools view -s 99.05`.

This method is very fast. Because it is probabilistic in nature, it does not sample an exact number of reads. Given 1 million reads and a sampling parameter of 0.05, it will iterate over the reads one by one, with a 0.05 probability of keeping each one. This will retain approximately 50000 reads, but the final total could equally well be 49993 or 50010.

##### Exact downsampling

We would like to sample an exact number of reads. The procedure recommended by the samtools developers is to use `samtools collate`, `awk` and `samtools sort`. This is rather slow and inefficient, making it impractical for very large inputs.

Unlike random downsampling, this method cannot take a random seed as input; given an input BAM file and output size, it will always select the same subset of the input.

### Predownsampling and downsampling

In order to combine the speed of random sampling with the precision of exact sampling, we implement downsampling as a two-stage process. For sufficiently large input, we first _predownsample_ randomly, getting an intermediate set somewhat larger than the desired output. Then, we _downsample_ exactly to get the final output set.

Predownsampling is applied if the original input is at least 2 times the size of the final output, _and_ the final output size is above an absolute minimum of 10000 reads. (For very small datasets, random sampling might accidentally produce an intermediate set smaller than the desired final size.) The predownsampling target size is 1.5 times the size of final output.

Predownsampling and downsampling are implemented as a single streamed command, so the predownsampled intermediate data is never written to disk.

### Examples

Let `N` be the size of the original input, and `T` be the target size for downsampled output.

- `N=25000, T=100000`. Because `N < T`, no downsampling of any kind takes place.
- `N=125000, T=100000`. Because `N < 2T`, no predownsampling takes place. The input is exactly downsampled to 100000 reads.
- `N=250000, T=100000`. Because `N > 2T` and `T > 10000` we apply predownsampling. The input is randomly predownsampled to an intermediate set of approximately 150000 reads. Then, the intermediate set is exactly downsampled to the final set of 100000 reads.
- `N=2500, T=1000`. Although `N > 2T`, we have `T < 10000`; so predownsampling is omitted, and we exactly downsample to 1000 reads.
