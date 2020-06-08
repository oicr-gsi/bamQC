# Filtering and downsampling in bamQC

## Overview

This document describes filtering and downsampling which may be applied to BAM files before bamQC metric computation.

_Filtering_ removes undesirable reads. _Downsampling_ obtains a subset of larger BAM files to enable faster processing. All downsampling operations, including checks of size thresholds, are applied after filtering. Read pairing is preserved; if read 1 of an aligned pair is downsampled, read 2 will be as well.

Picard MarkDuplicates has a separate downsampling method from other metrics. The two downsampling operations are implemented as separate WDL tasks, which can be efficiently parallelized by Cromwell.

Downsampling is not applied when using the mosdepth tool to compute depth of coverage.

## Filtering

We remove the following undesirable read types:
- **Non-primary alignments**: Important to avoid double-counting (or multi-counting) of reads. Defined as having samtools flag 256 (secondary alignment) or 2048 (supplementary alignment).
- **Unmapped reads**: Defined as having samtools flag 4 (segment unmapped).
- **Low-quality alignments**: Reads with a low alignment quality score (default minimum = 30)

Totals of non-primary, unmapped, and low-quality reads are recorded in workflow output. The excluded reads are not used for subsequent QC metric computation.

## Random Downsampling

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

We combine two ways of downsampling with samtools.

#### Approximate downsampling

This is done using `samtools view -s $SEED.$PROB`, where `$SEED` is a random seed and `$PROB` is a decimal expressing the probability of retaining a read. For example, to sample reads with a probability of 0.05 and random seed 99, use `samtools view -s 99.05`.

This method is very fast. Because it is probabilistic in nature, it does not sample an exact number of reads. Given 1 million reads and a sampling parameter of 0.05, it will iterate over the reads one by one, with a 0.05 probability of keeping each read. This will retain approximately 50000 reads, but the final total could plausibly be 49993 or 50010.

##### Exact downsampling

We would like to sample an exact number of reads. The procedure [recommended by the samtools developers](https://github.com/samtools/samtools/issues/931) is to use `samtools collate`, `awk` and `samtools sort`. This is rather slow and inefficient, making it impractical for very large inputs.

Unlike random downsampling, this method cannot take a random seed as input; given an input BAM file and output size, it will always select the same subset of the input.

### Predownsampling and downsampling

In order to combine the speed of random sampling with the precision of exact sampling, we implement downsampling as a two-stage process. For sufficiently large input, we first _predownsample_ approximately, getting an set of intermediate size somewhat larger than the desired output. Then, we _downsample_ exactly to get the final output set.

Predownsampling is applied if the original input is at least 2 times the size of the final output, _and_ the final output size is above an absolute minimum of 10000 reads. (For very small datasets, random sampling might accidentally produce an intermediate set smaller than the desired final size.) The predownsampling target size is 1.5 times the size of final output.

Predownsampling and downsampling are implemented as a single streamed command, so the predownsampled intermediate data is never written to disk.

### Examples

Let `N` be the size of the original input, and `T` be the target size for downsampled output.

- `N=25000, T=100000`. Because `N < T`, no downsampling of any kind takes place.
- `N=125000, T=100000`. Because `N < 2T`, no predownsampling takes place. The input is exactly downsampled to 100000 reads.
- `N=250000, T=100000`. Because `N > 2T` and `T > 10000` we apply predownsampling. The input is randomly predownsampled to an intermediate set of approximately 150000 reads. Then, the intermediate set is exactly downsampled to the final set of 100000 reads.
- `N=2500, T=1000`. Although `N > 2T`, we have `T < 10000`; so predownsampling is omitted, and we exactly downsample to 1000 reads.

### WDL Parameters for Random Downsampling


| Parameter            | Description                                     |
|----------------------|-------------------------------------------------|
| `targetReads`       | Desired number of reads in downsampled set. Most likely parameter to need adjustment; the others are details of implementation. |
| `minReadsAbsolute` | Absolute minimum number of reads to allow pre-downsmapling. Defaults to 10000. |
| `minReadsRelative`   | Minimum size of original dataset, relative to `targetReads`, in order for pre-downsampling to be carried out. Defaults to 2. |
| `precision`           | Number of digits to retain in the fractional string for pre-downsampling. Defaults to 8. |
| `preDSMultiplier`    |Approximate size of pre-downsampled set (if any), relative to `targetReads`. Defaults to 1.5. |

## Regional Downsampling for MarkDuplicates

### Motivation

Random downsampling is not recommended for MarkDuplicates. This is because the number of duplicate pairs retained scales as the square of the randomly downsampled fraction. (Recall that read pairing is always preserved by downsampling; a "duplicate pair" is a pair of paired reads.)

For example, suppose we have a BAM file with 250 million reads and a 10% duplicate rate, giving 25 million duplicates. We downsample to 1 million reads, retaining `1/250` of the original input. So the probability of downsampling *both* members of a duplicate is `(1/250)*(1/250) = 1/62500`, and the downsampled set has approximately 400 duplicate pairs. But because of the random nature of this process, the number sampled may be significantly larger or smaller. In empirical tests, naively scaling up from a sample of `1/62500` has given results at least 25% off the true value.

To obtain a more robust estimate, we use samtools to downsample reads aligned to a specific region of the human genome. Since both members of a duplicate should be aligned to the same locus, this method preserves duplicate pairs roughly in proportion to the fraction downsampled. In fact there will be some variation as duplicates are not evenly distributed across the genome, but the result is still more robust than random downsampling.

The number of reads downsampled by this method cannot be exactly predicted, as it depends on coverage of the chosen sub-region.

### Parameters

Downsampling is applied to specific chromosomes; and for very large files, to specific intervals within each chromosome. The chromosome and interval settings may be customised if needed.

#### WDL Parameters for Regional Downsampling

| Parameter                                            | Description          |
|------------------------------------------------------|----------------------|
| `threshold` | Minimum number of reads for downsampling. Defaults to 10 million. |
| `chromosomes` | Array of chromosome names for downsampling. Defaults to `["chr12", "chr13"]` |
| `baseInterval` | Base width of interval for downsampling. Defaults to 15000. |
| `intervalStart`| Start of downsampling interval on each chromosome. Defaults to 100000. |
| `customRegions` | Custom downsampling regions. Format is string input to samtools, eg. `chr1:1-1000000 chr2:10001-20000`. <br/>Defaults to the empty string `""`. If set to a value other than `""`, this will override the `chromosomes`, `baseInterval`, and `intervalStart` parameters. |

#### Custom regions

If the `customRegions` parameter is in effect, downsampling is very simple. If the number of reads is greater than `threshold`, reads are downsampled to the intervals specified in `customRegions`; otherwise, no downsampling is done.

#### Interval construction

The following only applies if `customRegions` is *not* in effect.

Let `R` be the number of reads and `T` be the minimum threshold for downsampling.

- If `R <= T`, no downsampling takes place.
- If `T < R <= 10T`, the entire sequence of the chromosomes specified in `chromosomes` is used for downsampling.
- If `10T < R <= 100T`, the same interval is used within each chromosome. The interval begins at `intervalStart + 1` and ends at `intervalStart + baseInterval*1000`.
- For larger values of R, the interval width is scaled down as shown in the table.

| Input reads (R) vs. threshold (T)                 | Width wrt base (B)  | Downsampling |
| --------------------------------------------------|---------------------|--------------|
| <pre>R <= T</pre>                                 | None                | No           |
| <pre>T < R <= 10T</pre>                           | (Entire chromosome) | Yes          |
| <pre>10T < R <= 10<sup>2</sup>T</pre>             | 10<sup>3</sup>B     | Yes          |
| <pre>10<sup>2</sup>T < R <= 10<sup>3</sup>T</pre> | 10<sup>2</sup>B     | Yes          |
| <pre>10<sup>3</sup>T < R <= 10<sup>4</sup>T</pre> | 10B                 | Yes          |
| <pre>R > 10<sup>4</sup>T</pre>                    | B                   | Yes          |

With the default interval values `T=10000000` and `B=15000`, and default chromosomes, we have:

| Input reads (R)                                   | DS chromosomes | Interval width      |
| --------------------------------------------------|----------------|---------------------|
| <pre>R <= 10<sup>7</sup></pre>                    | None           | None                |
| <pre>10<sup>7</sup> < R <= 10<sup>8</sup></pre>   | chr12 & chr13  | (Entire chromosome) |
| <pre>10<sup>8</sup> < R <= 10<sup>9</sup></pre>   | chr12 & chr13  | 1.5x10<sup>7</sup>  |
| <pre>10<sup>9</sup> < R <= 10<sup>10</sup></pre>  | chr12 & chr13  | 1.5x10<sup>6</sup>  |
| <pre>10<sup>10</sup> < R <= 10<sup>11</sup></pre> | chr12 & chr13  | 1.5x10<sup>5</sup>  |
| <pre>R > 10<sup>11</sup></pre>                    | chr12 & chr13  | 1.5x10<sup>4</sup>  |

The default parameters were chosen so that, when downsampling a whole-genome BAM file, the downsampled set is roughly between 1 million and 10 million reads. Empirically, this size range has been large enough to give a reasonably accurate sample, but small enough for Picard MarkDuplicates to be tractable. (Note that the default chromosomes 12 and 13 together are approximately 8% of the human genome.)

#### Parameters for targeted libraries

By default we downsample reads aligned to chromosomes 12 and 13, because they appear in all targeted sequencing panels in use at OICR as of 2020-03-10. For an up-to-date list of panels, see the [interval-files repository](https://bitbucket.oicr.on.ca/projects/GSI/repos/interval-files/browse).

We also need to consider what fraction of the targets falls within the targeted interval. Consider three example scenarios:

1. *Too few reads*: None of the targeted regions fall within the downsampled chromosomes. The downsampled set contains no reads of interest, so analysis fails.
2. *About enough reads*: Approximately 10% of the target regions fall within the downsampled chromosomes. Downsampling behaves roughly as on a whole-genome sample.
3. *Too many reads*: All of the targeted regions fall within the downsampled chromosomes. The downsampled set is roughly 10 times larger than for a whole-genome sample. This may be acceptable; but could cause problems if the downsampled set is too large for Picard MarkDuplicates to be run efficiently.

We need to choose a suitable *depth* and *range* of downsampling for a targeted library. This is done respectively by setting the `threshold` parameter; and setting either `customRegions` or the chromosome and interval parameters. For example, the "too many reads" scenario above can be addressed by reducing the threshold for downsampling.

