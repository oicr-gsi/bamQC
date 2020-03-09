# Filtering and downsampling in bamQC

## Overview

This document describes filtering and downsampling which may be applied to BAM files before bamQC metric computation.

_Filtering_ removes undesirable reads. _Downsampling_ obtains a subset of very large BAM files to enable faster processing. All downsampling operations, including checks of size thresholds, are applied after filtering. Read pairing is preserved; if read 1 of an aligned pair is downsampled, read 2 will be as well.

Sections are as follows:
- Filtering
- Downsampling: MarkDuplicates
- Downsampling: Other

## Filtering

We remove the following undesirable read types:
- Non-primary alignments: Important to avoid double-counting (or multi-counting) of reads. These are defined as having samtools flag 256 (secondary alignment) or 2048 (supplementary alignment).
- Unmapped reads. Defined as having samtools flag 4 (segment unmapped).
- Low-quality alignments: Reads with a low alignment quality score (default minimum = 30)

Totals of non-primary, unmapped, and low-quality reads are recorded in workflow output. The excluded reads are not used for subsequent QC metric computation.

## Downsampling: MarkDuplicates

MarkDuplicates is not amenable to random downsampling. This is because the number of duplicate pairs retained scales as the square of the downsampled fraction. (Recall that read pairing is always preserved by downsampling; a "duplicate pair" is a pair of paired reads.)

For example, suppose we have a BAM file with 250 million reads and a 10% duplicate rate, giving 25 million duplicates. We downsample to 1 million reads, retaining `1/250` of the original input. So the probability of downsampling *both* members of a duplicate is `(1/250)*(1/250) = 1/62500`, and the downsampled set has approximately 400 duplicate pairs. But because of the random nature of this process, the number sampled may be significantly larger or smaller. In empirical tests, naively scaling up from a sample of `1/62500` has given results at least 25% off the true value.

To obtain a more robust estimate, we use samtools to downsample reads aligned to a specific region of the human genome reference. Since both members of a duplicate should be aligned to the same locus, this method preserves duplicate pairs roughly in proportion to the fraction downsampled. (In fact there will be some variation as duplicates are not evenly distributed across the genome, but the result is still more robust than random downsampling.)

We apply downsampling for MarkDuplicates by choosing some or all of chromosome 1, which comprises about 8% of the human genome. The number of reads downsampled is not exact, as it depends on the coverage of our chosen region. Thresholds have been chosen to retain roughly between 1 million and 10 million reads. The lower end of this range is large enough for an informative sample of duplicate metrics; the upper end is small enough for MarkDuplicates to be run quickly without excessive memory usage.

### Downsampling MarkDuplicates thresholds

| Total input reads  | Downsample range         | Approx. fraction of genome |
| -------------------|--------------------------|----------------------------|
| <= 10<sup>7</sup>  | None                     | All                        |
| <= 10<sup>8</sup>  | chr1                     | 1/12                       |
| <= 10<sup>9</sup>  | chr1, 25 Mbase window    | 1/124                      |
| <= 10<sup>10</sup> | chr1, 2.5 Mbase window   | 1/1240                     |
| <= 10<sup>11</sup> | chr1, 0.25 Mbase window  | 1/12405                    |
| > 10<sup>11</sup>  | chr1, 0.025 Mbase window | 1/124049                   |

### TODO: Targeted sequencing

For targeted sequencing which excludes chromosome 1, the above method will not work. See the relevant [Github issue](https://github.com/oicr-gsi/bam-qc/issues/13) for details.

## Downsampling: Other

For "slow" computationally intensive metrics other than MarkDuplicates, we downsample by choosing a subset of reads at random. "Fast" metrics do not require downsampling, and are always computed on the full-sized filtered input.

The distinction between "fast" and "slow" metrics is covered in the [bam-qc-metrics documentation](https://github.com/oicr-gsi/bam-qc-metrics/blob/master/metrics.md).

### Threshold

If a BAM file has more reads than a threshold T, downsample to a subset containing T reads.

The default value of T is 100000. For example:
- A BAM file with 17000000 reads will be downsampled to 100000.
- A BAM file with 100001 reads will be downsampled to 100000.
- A BAM file with 99999 reads will not be downsampled.

So, if the original number of reads is 100000 or greater, downsampling is guaranteed not to reduce it below 100000.

Downsampling thresholds are always evaluated on the BAM file _after_ filtering.

### Methods

We use two ways of downsampling with samtools.

#### Random downsampling

This is done using `samtools view -s $FOO.$BAR`, where `$FOO` is a random seed and `$BAR` is a decimal expressing the probability of retaining a read. For example, to sample reads with a probability of 0.05 and random seed 99, use `samtools view -s 99.05`.

This method is very fast. Because it is probabilistic in nature, it does not sample an exact number of reads. Given 1 million reads and a sampling parameter of 0.05, it will iterate over the reads one by one, with a 0.05 probability of keeping each read. This will retain approximately 50000 reads, but the final total could plausibly be 49993 or 50010.

##### Exact downsampling

We would like to sample an exact number of reads. The procedure [recommended by the samtools developers](https://github.com/samtools/samtools/issues/931) is to use `samtools collate`, `awk` and `samtools sort`. This is rather slow and inefficient, making it impractical for very large inputs.

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
