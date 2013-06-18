#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;

use JSON::PP; # imports encode_json, decode_json, to_json and from_json

my $defaultSampleRate = 1000;
my $defaultInsertMax = 1500;
my $defaultBedFile = "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed";
my $defaultQualCut = 30;

# reads sam format from STDIN

sub usage
{
	print "\nUsage is samtools view aligned.sorted.bam | samStats.pl [options].\n";
	print "Bam file must be sorted by coordinates.\n";
	print "Options are as follows:\n";
	print "\t-s sample rate.  Defines how often to sample the reads (default $defaultSampleRate).\n";
	print "\t-i normal insert max.  Defines upper limit to what is considered a normal insert (default $defaultInsertMax).\n";
	print "\t-q mapping quality cut.  Reads that map with a quality worse than this will not be considered \"uniquely mapped\" (default $defaultQualCut).\n";
	print "\t-r target.bed.  Bed file containing targets to calculate coverage against (default $defaultBedFile).\n";
	print "\t\tNOTE: target.bed file MUST be sorted in the same order as the bam file.\n";
	print "\t-c enables % base covered reporting. Sets sample rate to 1 (overrides -s), and runs long.\n";
	print "\t-j additional JSON formatted data string FILEPATH!!. e.g. '{\"sample\":\"TLCR2C10n\",\"library\":\"TLCR_2C10_nn_n_PE_501_nn\",\"barcode\":\"TAGCTT\",\"instrument\":\"h802\",\"run name\":\"110616_SN802_0060_AC00FWACXX\",\"lane\":\"4\"}'\n";
	print "\t-g group id. Group ID provided by the investigator to guide sample grouping.\n";
	print "\t-d group id description. Text to describe the purpose of the Group ID. This string may be supplied by the investigator.\n";
	print "\t-n external name. The name the investigator uses to refer to the sample.\n";
	print "\t-w workflow name. Name of the workflow used to generate the bam file.\n";
	print "\t-v workflow version. Version of the workflow used to generate the bam file.\n";
	print "\t-t workflow timestamp. Workflow run creation timestamp (number of seconds since epoch).\n";
	print "\t-h displays this usage message.\n";

	die "\n@_\n\n";
}

my $sampleRate = $defaultSampleRate;
my $normalInsertMax = $defaultInsertMax;
my $bedFile = $defaultBedFile;
my $qualCut = $defaultQualCut;
my %jsonHash;
my $reportBasesCovered = 0;

my $group_id;
my $group_id_description;
my $external_name;
my $workflow_name;
my $workflow_version;
my $workflow_timestamp;

my $line;
my @fields;

my $opt_string = "s:i:l:r:j:q:c:g:d:n:w:v:t:h";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}
if (exists $opt{s})
{
	$sampleRate = $opt{s};
}
if (exists $opt{i})
{
	$normalInsertMax = $opt{i};
}
if (exists $opt{r})
{
	$bedFile = $opt{r};
}
if (exists $opt{q})
{
	$qualCut = $opt{q};
}

if (exists $opt{j})
{
	open(JSONFILE, "$opt{j}") or usage("Couldn't open $opt{j}\n");
	while($line = <JSONFILE>)
	{
		chomp $line;
		%jsonHash = %{ decode_json($line) };
	}
	close JSONFILE;
}

if (exists $opt{g}) {
	$group_id =  $opt{g};
}

if (exists $opt{d}) {
	$group_id_description = $opt{d};
}

if (exists $opt{n}) {
	$external_name = $opt{n};
}

if (exists $opt{w}) {
	$workflow_name = $opt{w};
}

if (exists $opt{v}) {
	$workflow_version = $opt{v};
}

if (exists $opt{t}) {
	$workflow_timestamp = $opt{t};
}


if (exists $opt{c})
{
	$reportBasesCovered = 1;
	$sampleRate = 1;
}

if ($reportBasesCovered ==  0)
{
	warn "Only sampling every $sampleRate reads.\n";
}
else
{
	warn "Sampling every read and calculating bases covered.\n";
}


my %bedStart;
my %bedStop;
my $bedPos = 0;
my $bedCurrentChr = "null";
my %bedHist;

my $mappedBases;
my $targetSize = 0;
my $numberOfTargets;
open (BEDFILE, $bedFile) or usage("Couldn't open target file: $bedFile.\n");

while ($line = <BEDFILE>)
{
	chomp $line;
	@fields = split(/\t/, $line);

	push (@{ $bedStart{$fields[0]} }, $fields[1]);
	push (@{ $bedStop{$fields[0]} }, $fields[2]);

	$targetSize += $fields[2] - $fields[1];

	$bedHist{"$fields[0]\t$fields[1]\t$fields[2]"} = 0;
}
close BEDFILE;

$numberOfTargets = scalar(keys %bedHist);
warn "Loaded " . $numberOfTargets . " targets.\n\n";


my @cigarOp;
my @mdOp;

my $timeToSampleCount = 0;

my $nonPrimaryReadCount = 0;
my $qualFailReadCount = 0;
my $mappedReadCount = 0;
my $unmappedReadCount = 0;

my $pairedReadCount = 0;
my $properPairReadCount = 0;
my $mateUnmappedReadCount = 0;

my $sampledCount = 0;
my $basesCount = 0;
my $mappedCount = 0;
my $alignedCount = 0;
my $softClipCount = 0;
my $hardClipCount = 0;
my $mismatchCount = 0;
my $insertCount = 0;
my $deletionCount = 0;

my $readsOnTarget = 0;

my $readsMissingMDtags = 0;

my $pairsMappedToDifferentChr = 0;
my $pairsMappedAbnormallyFar = 0;

my %normalInsertSizes;


my %startPointHash;
my $currentChr = "null";
my $currentStartPoint = "null";
my $currentReadsCount = 0;
my $readStart;
my $readsPerStartPoint = 1;   # wee bit sloppy, but saves a check at every start point
my $startPointCount = 0;


my $startPoint;
my $posOffset;
my $numStartPoints;
my $totalDepth;

my $readLength;
my $readLengthMean;
my %readLengthHist;

my $qualString;
my $cycle;

my %qualsByCycle;

my %alignedByCycle;
my %softClipByCycle;
my %hardClipByCycle;
my %mismatchByCycle;
my %insertionByCycle;
my %deletionByCycle;

my $firstOrSecondRead;

my %pairStartHash;

my %runningBaseCoverage;			# $runningBaseCoverage{chr}{pos}{startPoint}{count}
my %nonCollapsedCoverageHist;
my %collapsedCoverageHist;

my $onTarget;

my @sortedChars = qw(! " placeholderForPound $ % & ' \( \) * + placeholderForComma - . / 0 1 2 3 4 5 6 7 8 9 : ; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \\ ] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | } ~);
$sortedChars[2] = "#";
$sortedChars[11] = ",";


while ($line = <STDIN>)
{
	chomp $line;
	@fields = split(/\t/, $line);


	if ($fields[1] & 256)			# not primary alignment
	{
		$nonPrimaryReadCount++;
	}
	elsif ($fields[1] & 4)			# read unmapped
	{
		$unmappedReadCount++;
	}
	elsif ($fields[4] < $qualCut)			# below qual cut
	{
		$qualFailReadCount++;
	}
	else		# read is mapped (I hope)
	{
		if ($fields[1] & 1)
		{
			$pairedReadCount++;
			if ($fields[1] & 8)
			{
				$mateUnmappedReadCount++;
			}
		}
		if ($fields[1] & 2)
		{
			$properPairReadCount++;
		}

		$mappedReadCount++;

		if ("$fields[2]\t$fields[3]" ne $currentStartPoint)
		{
			for my $sp (keys %pairStartHash)
			{
				$startPointCount++;
				$readsPerStartPoint += ($pairStartHash{$sp} - $readsPerStartPoint) / $startPointCount;
				delete $pairStartHash{$sp};
			}
			
			$currentStartPoint = "$fields[2]\t$fields[3]";
		}

		$pairStartHash{"$fields[2]\t$fields[3]\t$fields[7]"}++;


		$timeToSampleCount++;
		if (($timeToSampleCount >= $sampleRate) or ($reportBasesCovered == 1))
		{
			$timeToSampleCount = 0;
	
	
			$sampledCount++;
	
#			if ($fields[2] ne $currentChr)
#			{
#				$currentChr = $fields[2];
#				warn "Processing $currentChr\n";
#			}
	
			@cigarOp = procCigar($fields[5]);
			@mdOp = ();
			foreach my $f (@fields)
			{
				if ($f =~ /MD:Z:(.*)/)
				{
					@mdOp = procMD($1);
				}
			}
			if (scalar(@mdOp) == 0)
			{
			#	warn "Couldn't find MD tag!\n";
			#	print $line . "\n";
				$readsMissingMDtags++;
			}
			
			if ($fields[1] & 16)		# read is mapped to reverse strand
			{
				@cigarOp = reverse(@cigarOp);
				@mdOp = reverse(@mdOp);
				$qualString = reverse($fields[10]);
			}
			else
			{
				$qualString = $fields[10];
			}

			$cycle = 1;
			$posOffset = 0;
			$mappedBases = 0;

			if ($fields[1] & 64)
			{
				$firstOrSecondRead = "R1";
			}
			elsif ($fields[1] & 128)
			{
				$firstOrSecondRead = "R2";
			}
			else
			{
				$firstOrSecondRead = "R?";
			}

			$onTarget = 0;
			if (exists $bedStart{$fields[2]})
			{
				if ($bedCurrentChr ne $fields[2])
				{
					$bedPos = 0;
					$bedCurrentChr = $fields[2];
				}

				while (($bedPos < scalar @{ $bedStop{$fields[2]} }) and ($fields[3] > $bedStop{$fields[2]}[$bedPos]))		# need to advance current bed target region
				{
					$bedPos++;
				}

				if (($bedPos < scalar @{ $bedStop{$fields[2]} }) and (($fields[3] <= $bedStop{$fields[2]}[$bedPos]) and (($fields[3] + $mappedBases) >= $bedStart{$fields[2]}[$bedPos])))
				{
					$readsOnTarget++;
					$onTarget = 1;
					$bedHist{"$fields[2]\t$bedStart{$fields[2]}[$bedPos]\t$bedStop{$fields[2]}[$bedPos]"} += $mappedBases;
				}
			}

			$startPoint = $fields[3];
			$readLength = 0;
			foreach my $c (@cigarOp)
			{
				if ($c =~ /(.*)M/)
				{
					$alignedCount += $1;
					for (my $i = 0; $i < $1; $i++)
					{
						$alignedByCycle{$firstOrSecondRead}{$cycle}++;
						$cycle++;
					}
					if ($reportBasesCovered == 1)
					{
						if ($onTarget == 1)
						{
							for (my $i = 0; $i < $1; $i++)
							{
								if ((($startPoint + $posOffset) >= $bedStart{$fields[2]}[$bedPos]) and (($startPoint + $posOffset) <= $bedStop{$fields[2]}[$bedPos]))
								{
									$runningBaseCoverage{$fields[2]}{$startPoint + $posOffset}{"$fields[2]\t$fields[3]\t$fields[7]"}++;
									$posOffset++;
								}
							}
						}
					}
					$readLength += $1;
					$mappedBases += $1;
				}
				elsif ($c =~ /(.*)H/)
				{
					$hardClipCount += $1;
					for (my $i = 0; $i < $1; $i++)
					{
						$hardClipByCycle{$firstOrSecondRead}{$cycle}++;
						$cycle++;
					}
					$readLength += $1;
				}
				elsif ($c =~ /(.*)S/)
				{
					$softClipCount += $1;
					for (my $i = 0; $i < $1; $i++)
					{
						$softClipByCycle{$firstOrSecondRead}{$cycle}++;
						$cycle++;
					}
					$readLength += $1;
				}
				elsif ($c =~ /(.*)I/)
				{
					$insertCount += $1;
					for (my $i = 0; $i < $1; $i++)
					{
						$insertionByCycle{$firstOrSecondRead}{$cycle}++;
						$cycle++;
					}
					$readLength += $1;
				}
				elsif ($c =~ /(.*)D/)
				{
					$deletionCount += $1;
					for (my $i = 0; $i < $1; $i++)
					{
						$deletionByCycle{$firstOrSecondRead}{$cycle}++;
					}
					if ($reportBasesCovered == 1)
					{
						if ($onTarget == 1)
						{
							for (my $i = 0; $i < $1; $i++)
							{
								if ((($startPoint + $posOffset) >= $bedStart{$fields[2]}[$bedPos]) and (($startPoint + $posOffset) <= $bedStop{$fields[2]}[$bedPos]))
								{
									$runningBaseCoverage{$fields[2]}{$startPoint + $posOffset}{"$fields[2]\t$fields[3]\t$fields[7]"}++;
									$posOffset++;
								}
							}
						}
					}

					$mappedBases += $1;
				}
				else
				{
					die "Can't handle CIGAR operation: $c\n";
				}
			}

			if ($reportBasesCovered == 1)
			{
				if ($onTarget == 1)
				{
				for my $chr (keys %runningBaseCoverage)
				{
					if ($chr ne $fields[2])
					{
						for my $p (keys %{ $runningBaseCoverage{$chr} })
						{
							$numStartPoints = scalar(keys %{ $runningBaseCoverage{$chr}{$p} });
							for (my $i = 1; $i <= $numStartPoints; $i++)
							{
								$collapsedCoverageHist{$i}++;
							}

							$totalDepth = 0;
							for my $sp (keys %{ $runningBaseCoverage{$chr}{$p} })
							{
								$totalDepth += $runningBaseCoverage{$chr}{$p}{$sp};
							}
							for (my $i = 1; $i <= $totalDepth; $i++)
							{
								$nonCollapsedCoverageHist{$i}++;
							}
						}
						delete $runningBaseCoverage{$chr};
					}
					else
					{
						for my $p (keys %{ $runningBaseCoverage{$chr} })
						{
							if ($p < $fields[3])
							{
								$numStartPoints = scalar(keys %{ $runningBaseCoverage{$chr}{$p} });
								for (my $i = 1; $i <= $numStartPoints; $i++)
								{
									$collapsedCoverageHist{$i}++;
								}

								$totalDepth = 0;
								for my $sp (keys %{ $runningBaseCoverage{$chr}{$p} })
								{
									$totalDepth += $runningBaseCoverage{$chr}{$p}{$sp};
								}
								for (my $i = 1; $i <= $totalDepth; $i++)
								{
									$nonCollapsedCoverageHist{$i}++;
								}
								delete $runningBaseCoverage{$chr}{$p};
							}
						}
					}

				}
				}
			}



			if (exists $readLengthHist{$firstOrSecondRead}{$readLength})
			{
				$readLengthHist{$firstOrSecondRead}{$readLength}++;
			}
			else
			{
				$readLengthHist{$firstOrSecondRead}{$readLength} = 1;
			}

			$cycle = 1;
			foreach my $md (@mdOp)
			{
				if ($md =~ /^([0-9]+)/)
				{
					# ignore matching bases
					$cycle += $1;
				}
				elsif ($md =~ /^(\^[A-Z]+)/)
				{
					# ignore deletions
				}
				elsif ($md =~ /^([A-Z]+)/)		# mismatch!
				{
					foreach my $i (split(//, $1))
					{
						$mismatchCount++;
						$mismatchByCycle{$firstOrSecondRead}{$cycle}++;
						$cycle++;
					}
				}
				else
				{
					die "Couldn't handle MD operation: $md\n";
				}
			}


			if ($fields[8] > 0)
			{
				if ($fields[6] eq "=")
				{
					if ($fields[8] < $normalInsertMax)
					{
						if (exists $normalInsertSizes{$fields[8]})
						{
							$normalInsertSizes{$fields[8]}++;
						}
						else
						{
							$normalInsertSizes{$fields[8]} = 1;
						}
					}
					else # abnormally long read
					{
						$pairsMappedAbnormallyFar++;
					}
				}
				else		# pair is mapped to a different chromosome
				{
					$pairsMappedToDifferentChr++;
				}
	
			}


			# grab base quality values

			$cycle = 1;
			for my $q (split(//, $qualString))
			{
				$qualsByCycle{$firstOrSecondRead}{$q}{$cycle}++;
				$cycle++;
			}
		}
	}
}

my $averageReadLength = 0;
if ($mappedReadCount > 0)
{
	$averageReadLength = (($alignedCount + $softClipCount + $hardClipCount + $insertCount) / $mappedReadCount) * $sampleRate;
}

my $meanInsert = meanHist(\%normalInsertSizes);
my $stdevInsert = stdevHist(\%normalInsertSizes);

my $maxInsert = 0;
my $minInsert = 99999;
for my $i (keys %normalInsertSizes)
{
	if ($i < $minInsert)
	{
		$minInsert = $i;
	}
	if ($i > $maxInsert)
	{
		$maxInsert = $i;
	}
}
for (my $i = $minInsert; $i <= $maxInsert; $i++)
{
	unless (exists $normalInsertSizes{$i})
	{
		$normalInsertSizes{$i} = 0;
	}
}

if ($reportBasesCovered == 1)
{
	# process any remaining coverages

	for my $chr (keys %runningBaseCoverage)
	{
		for my $p (keys %{ $runningBaseCoverage{$chr} })
		{
			$numStartPoints = scalar(keys %{ $runningBaseCoverage{$chr}{$p} });
			for (my $i = 1; $i <= $numStartPoints; $i++)
			{
				$collapsedCoverageHist{$i}++;
			}

			$totalDepth = 0;
			for my $sp (keys %{ $runningBaseCoverage{$chr}{$p} })
			{
				$totalDepth += $runningBaseCoverage{$chr}{$p}{$sp};
			}
			for (my $i = 1; $i <= $totalDepth; $i++)
			{
				$nonCollapsedCoverageHist{$i}++;
			}
		}
	}
}

my %qualHist;
my $qualSum;
my $qualCount;

my %qualLine;

my %maxReadLength;

for my $ends (qw(R1 R2 R?))
{
	$maxReadLength{$ends} = 0;
	for my $q (@sortedChars)
	{
		if (exists $qualsByCycle{$ends}{$q})
		{
			for my $c (keys %{ $qualsByCycle{$ends}{$q} })
			{
				if ($c > $maxReadLength{$ends})
				{
					$maxReadLength{$ends} = $c;
				}
			}
		}
	}
}

for my $ends (qw(R1 R2 R?))
{
	for (my $c = 1; $c <= $maxReadLength{$ends}; $c++)
	{
		unless (exists $alignedByCycle{$ends}{$c})
		{
			$alignedByCycle{$ends}{$c} = 0;
		}
		unless (exists $mismatchByCycle{$ends}{$c})
		{
			$mismatchByCycle{$ends}{$c} = 0;
		}
		unless (exists $insertionByCycle{$ends}{$c})
		{
			$insertionByCycle{$ends}{$c} = 0;
		}
		unless (exists $deletionByCycle{$ends}{$c})
		{
			$deletionByCycle{$ends}{$c} = 0;
		}
		unless (exists $softClipByCycle{$ends}{$c})
		{
			$softClipByCycle{$ends}{$c} = 0;
		}
		unless (exists $hardClipByCycle{$ends}{$c})
		{
			$hardClipByCycle{$ends}{$c} = 0;
		}
	}
}


for my $ends (qw(R1 R2 R?))
{
	for (my $c = 1; $c <= $maxReadLength{$ends}; $c++)
	{
		$qualSum = 0;
		$qualCount = 0;
		for my $q (@sortedChars)
		{
			if (exists $qualsByCycle{$ends}{$q}{$c})
			{
				$qualCount += $qualsByCycle{$ends}{$q}{$c};
				$qualSum += $qualsByCycle{$ends}{$q}{$c} * toPhred($q);
	
				if ($qualsByCycle{$ends}{$q}{$c} > 0)
				{
					if (exists $qualHist{$ends}{toPhred($q)})
					{
						$qualHist{$ends}{toPhred($q)} += $qualsByCycle{$ends}{$q}{$c};
					}
					else
					{
						$qualHist{$ends}{toPhred($q)} = $qualsByCycle{$ends}{$q}{$c};
					}
				}
	
			}
		}
	
		unless ($qualCount == 0)
		{
			$qualLine{$ends}{$c} = int($qualSum / $qualCount);
		}
		else
		{
			$qualLine{$ends}{$c} = 0;
		}
	}
}


for my $ends (qw(R1 R2 R?))
{
	for my $q (@sortedChars)
	{
		if (exists $qualHist{$ends}{toPhred($q)})
		{
			$qualHist{$ends}{toPhred($q)} = $qualHist{$ends}{toPhred($q)} * $sampleRate;
		}
	}

	for my $c (keys %{ $alignedByCycle{$ends} })
	{
		$alignedByCycle{$ends}{$c} = $alignedByCycle{$ends}{$c} * $sampleRate;
	}
	for my $c (keys %{ $mismatchByCycle{$ends} })
	{
		$mismatchByCycle{$ends}{$c} = $mismatchByCycle{$ends}{$c} * $sampleRate;
	}
	for my $c (keys %{ $insertionByCycle{$ends} })
	{
		$insertionByCycle{$ends}{$c} = $insertionByCycle{$ends}{$c} * $sampleRate;
	}
	for my $c (keys %{ $deletionByCycle{$ends} })
	{
		$deletionByCycle{$ends}{$c} = $deletionByCycle{$ends}{$c} * $sampleRate;
	}
	for my $c (keys %{ $softClipByCycle{$ends} })
	{
		$softClipByCycle{$ends}{$c} = $softClipByCycle{$ends}{$c} * $sampleRate;
	}
	for my $c (keys %{ $hardClipByCycle{$ends} })
	{
		$hardClipByCycle{$ends}{$c} = $hardClipByCycle{$ends}{$c} * $sampleRate;
	}

	for my $l (keys %{ $readLengthHist{$ends} })
	{
		$readLengthHist{$ends}{$l} = $readLengthHist{$ends}{$l} * $sampleRate;
	}
}


for my $i (keys %normalInsertSizes)
{
	$normalInsertSizes{$i} = $normalInsertSizes{$i} * $sampleRate;
}

my $lengthSum;
my $lengthCount;

my %averageReadLength;
for my $ends (qw(R1 R2 R?))
{
	$lengthSum = 0;
	$lengthCount = 0;
	for my $l (keys %{ $readLengthHist{$ends} })
	{
		$lengthSum += $readLengthHist{$ends}{$l} * $l;
		$lengthCount += $readLengthHist{$ends}{$l};
	}
	unless ($lengthCount == 0)
	{
		$averageReadLength{$ends} = $lengthSum / $lengthCount;
	}
	else
	{
		$averageReadLength{$ends} = 0;
	}
}

warn "\n\n";
my $totalReads = $nonPrimaryReadCount + $qualFailReadCount + $unmappedReadCount + $mappedReadCount;
warn "Total reads: $totalReads\n";
warn "Mapped reads: $mappedReadCount\n";
warn "Non primary reads: $nonPrimaryReadCount\n";
warn "MAPQ < $qualCut reads: $qualFailReadCount\n";
warn "Unmapped reads: $unmappedReadCount\n";

warn "Sampled reads: $sampledCount\n";
warn "Reads on target: $readsOnTarget\n\n";

warn "Reads missing MD tags!: $readsMissingMDtags\n\n";

warn "Aligned bases: $alignedCount\n";
warn "Soft clipped bases: $softClipCount\n";
warn "Hard clipped bases: $hardClipCount\n";
warn "Mismatched bases: $mismatchCount\n";
warn "Deleted base count: $deletionCount\n";
warn "Inserted base count: $insertCount\n";
warn "Average read length: $averageReadLength\n\n";

warn "Mean insert: $meanInsert\n";
warn "Stdev insert: $stdevInsert\n";
warn "Pairs with insert longer than $normalInsertMax: $pairsMappedAbnormallyFar\n";
warn "Pairs mapped to different chromosomes: $pairsMappedToDifferentChr\n\n";

warn "Reads per start point: $readsPerStartPoint\n\n";

#warn "Mismatches by cycle: ";
#printHist(\%mismatchByCycle);
#warn "\n";
#warn "Deletions by cycle: ";
#printHist(\%deletionByCycle);
#warn "\n";
#warn "Insertions by cycle: ";
#printHist(\%insertionByCycle);
#warn "\n";
#warn "Soft clips by cycle: ";
#printHist(\%softClipByCycle);
#warn "\n";
#warn "Hard clips by cycle: ";
#printHist(\%hardClipByCycle);
#warn "\n";
#warn "Average quality by cycle: ";
#printHist(\%qualHist);
#warn "\n";


my $coverage;
my $chr;
my $start;
my $end;
for my $reg (sort keys %bedHist)
{
	($chr, $start, $end) = split(/\t/, $reg);


	$coverage = int(($bedHist{$reg} / ($end - $start)) * $sampleRate);

	if ($coverage > 0)
	{
#		print "$chr:$start-$end: ${coverage}x\n";
	}
}

my $numEnds = "single end";
if ($properPairReadCount > 0)
{
	$numEnds = "paired end";
}

if ($group_id) {
	$jsonHash{"group id"} = $group_id;
}

if ($group_id_description) {
	$jsonHash{"group id description"} = $group_id_description;
}

if ($external_name) {
	$jsonHash{"external name"} = $external_name;
}

if ($workflow_name) {
	$jsonHash{"workflow name"} = $workflow_name;
}

if ($workflow_version) {
	$jsonHash{"workflow version"} = $workflow_version;
}

if ($workflow_timestamp) {
	$jsonHash{"workflow run creation timestamp"} = $workflow_timestamp;
}

$jsonHash{"number of ends"} = $numEnds;

$jsonHash{"average read length"} = $averageReadLength;
$jsonHash{"insert mean"} = $meanInsert;
$jsonHash{"insert stdev"} = $stdevInsert;

$jsonHash{"total reads"} = $mappedReadCount + $unmappedReadCount + $qualFailReadCount + $nonPrimaryReadCount;;
$jsonHash{"mapped reads"} = $mappedReadCount;
$jsonHash{"unmapped reads"} = $unmappedReadCount;
$jsonHash{"non primary reads"} = $nonPrimaryReadCount;
$jsonHash{"paired reads"} = $pairedReadCount;
$jsonHash{"properly paired reads"} = $properPairReadCount;
$jsonHash{"mate unmaped reads"} = $mateUnmappedReadCount;

$jsonHash{"qual fail reads"} = $qualFailReadCount;
$jsonHash{"qual cut"} = $qualCut;

$jsonHash{"aligned bases"} = ($alignedCount + $insertCount) * $sampleRate;
$jsonHash{"mismatch bases"} = $mismatchCount * $sampleRate;
$jsonHash{"inserted bases"} = $insertCount * $sampleRate;
$jsonHash{"deleted bases"} = $deletionCount * $sampleRate;
$jsonHash{"soft clip bases"} = $softClipCount * $sampleRate;
$jsonHash{"hard clip bases"} = $hardClipCount * $sampleRate;

$jsonHash{"reads per start point"} = $readsPerStartPoint;
$jsonHash{"reads on target"} = $readsOnTarget * $sampleRate;
$jsonHash{"target file"} = $bedFile;
$jsonHash{"target size"} = $targetSize;
$jsonHash{"number of targets"} = $numberOfTargets;

if ($reportBasesCovered == 1)
{
	$jsonHash{"non collapsed bases covered"} = \%nonCollapsedCoverageHist;
	$jsonHash{"collapsed bases covered"} = \%collapsedCoverageHist;
}

$jsonHash{"insert histogram"} = \%normalInsertSizes;

if (exists $qualHist{"R1"})
{
	$jsonHash{"read 1 quality histogram"} = $qualHist{"R1"};
	$jsonHash{"read 1 quality by cycle"} = $qualLine{"R1"};
	$jsonHash{"read 1 length histogram"} = $readLengthHist{"R1"};
	$jsonHash{"read 1 average length"} = $averageReadLength{"R1"};
	$jsonHash{"read 1 aligned by cycle"} = $alignedByCycle{"R1"};
	$jsonHash{"read 1 mismatch by cycle"} = $mismatchByCycle{"R1"};
	$jsonHash{"read 1 insertion by cycle"} = $insertionByCycle{"R1"};
	$jsonHash{"read 1 deletion by cycle"} = $deletionByCycle{"R1"};
	$jsonHash{"read 1 soft clip by cycle"} = $softClipByCycle{"R1"};
	$jsonHash{"read 1 hard clip by cycle"} = $hardClipByCycle{"R1"};
}
if (exists $qualHist{"R2"})
{
	$jsonHash{"read 2 quality histogram"} = $qualHist{"R2"};
	$jsonHash{"read 2 quality by cycle"} = $qualLine{"R2"};
	$jsonHash{"read 2 length histogram"} = $readLengthHist{"R2"};
	$jsonHash{"read 2 average length"} = $averageReadLength{"R2"};
	$jsonHash{"read 2 aligned by cycle"} = $alignedByCycle{"R2"};
	$jsonHash{"read 2 mismatch by cycle"} = $mismatchByCycle{"R2"};
	$jsonHash{"read 2 insertion by cycle"} = $insertionByCycle{"R2"};
	$jsonHash{"read 2 deletion by cycle"} = $deletionByCycle{"R2"};
	$jsonHash{"read 2 soft clip by cycle"} = $softClipByCycle{"R2"};
	$jsonHash{"read 2 hard clip by cycle"} = $hardClipByCycle{"R2"};
}
if (exists $qualHist{"R?"})
{
	$jsonHash{"read ? quality histogram"} = $qualHist{"R?"};
	$jsonHash{"read ? quality by cycle"} = $qualLine{"R?"};
	$jsonHash{"read ? length histogram"} = $readLengthHist{"R?"};
	$jsonHash{"read ? average length"} = $averageReadLength{"R?"};
	$jsonHash{"read ? aligned by cycle"} = $alignedByCycle{"R?"};
	$jsonHash{"read ? mismatch by cycle"} = $mismatchByCycle{"R?"};
	$jsonHash{"read ? insertion by cycle"} = $insertionByCycle{"R?"};
	$jsonHash{"read ? deletion by cycle"} = $deletionByCycle{"R?"};
	$jsonHash{"read ? soft clip by cycle"} = $softClipByCycle{"R?"};
	$jsonHash{"read ? hard clip by cycle"} = $hardClipByCycle{"R?"};
}


print encode_json(\%jsonHash);


sub initializeCycleHash
{
	my $hash = $_[0];
	my $max = $_[1];

	for (my $i = 1; $i < $max; $i++)
	{
		$hash->{$i} = 0;
	}

}


sub printHist
{
	my %hist = %{ $_[0] };

	for my $i (sort {$a <=> $b} keys %hist)
	{
#		warn " $hist{$i}";
	}
}

sub toPhred
{
	my $char = $_[0];
	my $ascii = ord($char);
	my $offset = 33;
	return $ascii - $offset;
}


sub mean
{
	my $val = $_[0];
	my $sum = 0;
	my $count = 0;
	for my $v (@{ $val })
	{
		$sum += $v;
		$count++;
	}
	if ($count > 0)
	{
		return $sum / $count;
	}
	else
	{
		return 0;
	}
}

sub stdev
{
	my $val = $_[0];
	my $mean = mean($val);
	my $squareDiff = 0;
	my $count = 0;
	for my $v (@{ $val })
	{
		$squareDiff += (($v - $mean)*($v - $mean));
		$count++;
	}
	if ($count > 0)
	{
		return sqrt($squareDiff / $count);
	}
	else
	{
		return 0;
	}
}

sub meanHist
{
	my %val = %{ $_[0] };
	my $sum = 0;
	my $count = 0;
	for my $v (keys %val)
	{
		$sum += ($v * $val{$v});
		$count += $val{$v};
	}

	if ($count > 0)
	{
		return $sum / $count;
	}
	else
	{
		return 0;
	}
}

sub stdevHist
{
	my %val = %{ $_[0] };
	my $mean = meanHist(\%val);
	my $squareDiff = 0;
	my $count = 0;
	for my $v (keys %val)
	{
		$squareDiff += ((($v - $mean)*($v - $mean)) * $val{$v});
		$count += $val{$v};
	}
	if ($count > 0)
	{
		return sqrt($squareDiff / $count);
	}
	else
	{
		return 0;
	}
}



sub procCigar
{
    my $cigar = $_[0];
    my @cigarOp;

    while ($cigar =~ /^([0-9]+[MIDNSHPX=]).*$/)
    {
        push (@cigarOp, $1);
        $cigar =~ s/$1//;
    }

    return @cigarOp;
}

sub procMD
{
    my $md = $_[0];
    my @mdOp;

    while ($md ne "")
    {
        if ($md =~ /^([0-9]+)/)
        {
            push(@mdOp, $1);
            $md =~ s/^$1//;
        }
        if ($md =~ /^([A-Z]+)/)
        {
            push(@mdOp, $1);
            $md =~ s/^$1//;
        }
        if ($md =~ /^(\^)([A-Z]+)/)
        {
            push(@mdOp, "^$2");
            $md =~ s/^\^$2//;
        }
    }

    return @mdOp;
}

sub findStart
{
        my @cigarOp = @{ $_[0] };
        my $start = $_[1];

        if ($cigarOp[0] =~ /(.*)S/)     # if first cigar operation is a soft clip, adjust start point
        {
                $start -= $1;
        }

        return $start;
}

sub findEnd
{
        my @cigarOp = @{ $_[0] };

        my $end = $_[1];

        if ($cigarOp[0] =~ /(.*)S/)     # if first cigar operation is a soft clip, adjust start point
        {
                $end -= $1;
        }

        foreach my $cig (@cigarOp)
        {
                if ($cig =~ /(.*)S/)
                {
                        $end += $1;
                }
                elsif ($cig =~ /(.*)M/)
                {
                        $end += $1;
                }
                elsif ($cig =~ /(.*)D/)
                {
                        $end += $1;
                }
        }

        return $end;
}

