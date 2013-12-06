#!/usr/bin/perl -w
$|++;
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use File::Copy;

######################################################################################################################################################
#
#	Description
#		This is a perl script to compare the multiple pileup perl storables generated wiggleToPerlStorable.
#
#	Input
#		--plsIndexListPath=		path; a file contains all paths of the index.hsh.pls, in format of: sample\tplsIndexPath\n;
#		--refFastaPath=			reference genome sequence, used for getting the name of the contig, as the reference to look for contig data in the pileup files;
#		--strandSpecific=		'yes' or 'no';[no] strand specific or not; if yes, will compare + and - seperately; if not, will merge + and -;
#		--refGffPath=			path of a GFF;
#		--compareMode=			"refBased" or "nrPair" or "allPair"; ["nrPair"] "refBased" will compare all others samples against a user-defined sample (from stdin); "nrPair" will compare non-redundant pairs only; "allPair" will compare all possible pairs, including the redudnant pairs like sample1.vs.sample2 and sample2.vs.sample1;
#		--outDir=				output directory; default = ./wiggleMerger/
#
#	Output
#		
#	Usage
#		./pileupPerlStorableComparer_v0.1.pl --plsIndexListPath=PlasmodiumTexVsPolyA.txt --refFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/Pf3D7_01_v3_JAN2012_withMitoPlstd.fa
#
#	History:
#		
#		V0.1
#			debut
#
######################################################################################################################################################
#	
#==========================================================Main body starts==========================================================================#
#----------Read parameters ----------#
my ($plsIndexListPath, $refFastaPath, $strandSpecific, $refGffPath, $compareMode, $outDir) = &readParameters();
printCMDLogOrFinishMessage("CMDLog");

#----------Read contig names
my $cntgLenHsh_ref = readContigLengthFromFasta();

my ($SSRngByCntgByWholeGeneHsh_ref, $ctgryByGeneHsh_ref, $strandHsh_ref, $XstrandRngByCntgByWholeGeneHsh_ref, $exonRngByGeneHsh_ref, $CDSRngByGeneHsh_ref, $intronRngByGeneHsh_ref, $junctStrByIntronIDHsh_ref, $intronRngXS_ref, $strandByIntronHsh_ref, $cntgByGeneHsh_ref, $strandByJunctStrHsh_ref, $descriptionByGeneHsh_ref, $exonLengthByGeneHsh_ref) = readGff($refGffPath);

#----------generate Extended Known Feature Rng
my $extendRange = 200;
my $extendedRngHsh_ref = generateExtendedKnownFeatureRng($XstrandRngByCntgByWholeGeneHsh_ref, $extendRange, $ctgryByGeneHsh_ref);

#----------Check pls index
my ($cntgCovPlsPathBySamByCntgHsh_ref, $comparisonInfoHsh_ref, $plsIndexPathHsh_ref) = getCntgCovPlsPaths($plsIndexListPath, $compareMode, $outDir);

#----decompress the storables, assuming they were compressed in gz (it's ok if even they are not compressed, wont crash)
compressOrDecompressAllStorable($cntgCovPlsPathBySamByCntgHsh_ref, "decompress");

#----Identify transfrags from all samples
my $transfragCovCutoff = 2;
my $minContinuousRegSize = 200;
my $maxGapTorlerance = 50;
my @strandToCompareAry = qw /+ -/;
@strandToCompareAry = qw /./ if $strandSpecific eq 'no';
my ($allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, $allTransfragInfoHsh_ref) = generateTransfragsForAllSamples(\@strandToCompareAry, $minContinuousRegSize, $maxGapTorlerance, $transfragCovCutoff, $plsIndexListPath, $outDir, $cntgLenHsh_ref, $plsIndexPathHsh_ref, $extendedRngHsh_ref, $strandHsh_ref);

#----overlay the transfrags of different samples
($allTransfragInfoHsh_ref, $comparisonInfoHsh_ref) = compareTransfragSamplePair($comparisonInfoHsh_ref, $cntgLenHsh_ref, $outDir, $allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, $allTransfragInfoHsh_ref);
my %allTransfragInfoHsh = %{$allTransfragInfoHsh_ref};
my %comparisonInfoHsh = %{$comparisonInfoHsh_ref};

#---print all transfrag info
printAllTransfragInfo($allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, \%allTransfragInfoHsh, $outDir);

#----overlay the coverage at different cutoffs
my $posCovCutoff = 5;#---hard coded cutoffs
$comparisonInfoHsh_ref = compareCntgCovPlsSamplePair($cntgCovPlsPathBySamByCntgHsh_ref, \%comparisonInfoHsh, $cntgLenHsh_ref, $posCovCutoff, $outDir, \@strandToCompareAry);
%comparisonInfoHsh = %{$comparisonInfoHsh_ref};

#----print comparisonInfoHsh
printComparisonInfoHsh(\%comparisonInfoHsh, $outDir);

#----compress all storables for better storage
compressOrDecompressAllStorable($cntgCovPlsPathBySamByCntgHsh_ref, "compress");

printCMDLogOrFinishMessage("finishMessage");

exit;

#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	my $dirPath = dirname(rel2abs($0));
	my $outDir = "$dirPath/pileupPerlStorableComparer/";
	my $strandSpecific = 'no';
	my $refGffPath = undef;
	my $refFastaPath = undef;
	my $compareMode = 'nrPair';

	foreach my $param (@ARGV) {
		if ($param =~ m/--plsIndexListPath=/) {$plsIndexListPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--refFastaPath=/) {$refFastaPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--strandSpecific=/) {$strandSpecific = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--refGffPath=/) {$refGffPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--compareMode=/) {$compareMode = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	#---check the files
	open (TEST, "$refFastaPath") || die "Can't open refFastaPath\n"; close TEST;
	open (TEST, "$refGffPath") || die "Can't open refGffPath\n"; close TEST;
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -pm 777 $outDir/storable";
	system "mkdir -pm 777 $outDir/GFF";
	
	return ($plsIndexListPath, $refFastaPath, $strandSpecific, $refGffPath, $compareMode, $outDir);
}
########################################################################## getCntgCovPlsPaths
sub getCntgCovPlsPaths {
	
	my ($plsIndexListPath, $compareMode, $userRef, $outDir) = @_;
	
	my %plsIndexPathHsh = ();
	my @sampleNameAry = ();
	my %cntgCovPlsPathBySamByCntgHsh = ();
	
	open (IDXLISTPATH, $plsIndexListPath);
	while (<IDXLISTPATH>) {
		chomp;
		next if $_ =~ m/^#/;
		my ($sample, $plsIndexPath) = split /\t/;
		my ($plsIndexName, $plsIndexDir, $plsIndexSuffix) = fileparse($plsIndexPath, qr/\.[^.]*/);
		die "wiggle file of $sample has to be .pls\n" if ($plsIndexSuffix ne ".pls");
		die "wiggle file of $sample doesn't exists\n" unless (-s $plsIndexPath);
		print "$sample .pls index checked\n";
		$plsIndexPathHsh{$sample} = $plsIndexPath;
		push @sampleNameAry, $sample;
		my %plsIndexHsh = %{retrieve($plsIndexPath)};
		my (undef, $cntgCovStroableDir, undef) = fileparse($plsIndexPath, qr/\.[^.]*/);
		foreach my $cntg (keys %plsIndexHsh) {
			my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
			${$cntgCovPlsPathBySamByCntgHsh{$sample}}{$cntg} = $cntgCovPlsPath;
		}
	}
	close IDXLISTPATH;

	#---generate pairs of comparison
	@sampleNameAry = sort {$a cmp $b} @sampleNameAry;
	my %comparisonInfoHsh = ();
	if ($compareMode eq 'nrPair') {
		foreach my $iRef(0..($#sampleNameAry - 1)) {
			foreach my $iQry(($iRef + 1)..$#sampleNameAry) {
				my $refSample = $sampleNameAry[$iRef];
				my $qrySample = $sampleNameAry[$iQry];
				my $comparisonTag = join ".vs.", ($refSample, $qrySample); 
				${$comparisonInfoHsh{$comparisonTag}}{'refSample'} = $refSample;
				${$comparisonInfoHsh{$comparisonTag}}{'qrySample'} = $qrySample;
			}
		}
	} elsif ($compareMode eq 'allPair') {
	
		foreach my $iRef(0..$#sampleNameAry) {
			foreach my $iQry(0..$#sampleNameAry) {
				my $refSample = $sampleNameAry[$iRef];
				my $qrySample = $sampleNameAry[$iQry];
				next if $refSample eq $qrySample;
				my $comparisonTag = join ".vs.", ($refSample, $qrySample); 
				${$comparisonInfoHsh{$comparisonTag}}{'refSample'} = $refSample;
				${$comparisonInfoHsh{$comparisonTag}}{'qrySample'} = $qrySample;
			}
		}
	
	} elsif($compareMode eq 'refBased') {
		
		print "Please enter the sample number to be used as the reference:\n";
		foreach my $i (0..$#sampleNameAry) {
			print "$i) $sampleNameAry[$i]\n";
		}
		my $iRef = -1;
		while (($iRef < 0) or ($iRef > $#sampleNameAry)) {
			print "your choice:";
			chomp ($iRef = <STDIN>);
		}
		
		my $refSample = $sampleNameAry[$iRef];
		foreach my $iQry(0..$#sampleNameAry) {
			my $qrySample = $sampleNameAry[$iQry];
			next if $refSample eq $qrySample;
			my $comparisonTag = join ".vs.", ($refSample, $qrySample); 
			${$comparisonInfoHsh{$comparisonTag}}{'refSample'} = $refSample;
			${$comparisonInfoHsh{$comparisonTag}}{'qrySample'} = $qrySample;
		}
	
	} else {
		
		die "compareMode has to be nrPair, allPair or refBased\n";
	}

	return \%cntgCovPlsPathBySamByCntgHsh, \%comparisonInfoHsh, \%plsIndexPathHsh;
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		my $dirPath = dirname(rel2abs($0));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($0, qr/\.[^.]*/);
		print CMDLOG "[".$runTime."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readContigLengthFromFasta
sub readContigLengthFromFasta {

	print "Reading $refFastaPath for contig names.\n";
	open (INFILE, $refFastaPath);
	my (%cntgLenHsh, $seqName, $length, $seq);

	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
			
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			$seqName = $curntLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} else {#--seq line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($curntLine);
		}
			
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seqName = $nextLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} elsif (eof(INFILE)) {#---this is the last line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($nextLine);
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;

	my $contigNum = keys %cntgLenHsh;

	print "Totally $contigNum contig names stored.\n";
	
	#---print to check
	#foreach my $cntg (sort {$a cmp $b} keys %cntgLenHsh) {print $cntg."\t".$cntgLenHsh{$cntg}."\n";}

	return (\%cntgLenHsh);
}
########################################################################## compressOrDecompressAllStorable
sub compressOrDecompressAllStorable {

	#compressOrDecompressAllStorable($cntgCovPlsPathBySamByCntgHsh_ref, $compressOrDecompress);

	my ($cntgCovPlsPathBySamByCntgHsh_ref, $compressOrDecompress) = @_;

	#---check for pigz
	my $compresser = "gzip";
	my $stdout = `pigz --version 2>&1`;
	if ($stdout !~ m/command not found/) {
		$compresser = 'pigz';
	}
	
	print "Try $compressOrDecompress storables using $compresser\n";

	my %cntgCovPlsPathBySamByCntgHsh = %{$cntgCovPlsPathBySamByCntgHsh_ref};

	foreach my $sample (keys %cntgCovPlsPathBySamByCntgHsh) {
		foreach my $cntg (keys %{$cntgCovPlsPathBySamByCntgHsh{$sample}}) {
			my $cntgCovPlsPath = ${$cntgCovPlsPathBySamByCntgHsh{$sample}}{$cntg};
			if (($compressOrDecompress eq "compress") and (-s $cntgCovPlsPath)) {
				system ("$compresser -f $cntgCovPlsPath");
			} elsif (($compressOrDecompress eq "decompress") and (-s $cntgCovPlsPath.".gz")) {
				system ("$compresser -df $cntgCovPlsPath.gz");
			} else {
				next;
			}
		}
	}
}
########################################################################## compareCntgCovPlsSamplePair
sub compareCntgCovPlsSamplePair {


	my ($cntgCovPlsPathBySamByCntgHsh_ref, $comparisonInfoHsh_ref, $cntgLenHsh_ref, $posCovCutoff, $outDir, $strandToCompareAry_ref) = @_;
	
	my %cntgCovPlsPathBySamByCntgHsh = %{$cntgCovPlsPathBySamByCntgHsh_ref};
	my %comparisonInfoHsh = %{$comparisonInfoHsh_ref};
	my @strandToCompareAry = @{$strandToCompareAry_ref};
	my %cntgLenHsh = %{$cntgLenHsh_ref};

	#---define overall hashes
	my %overlapHsh = ();
	my $cntgNum = keys %cntgLenHsh;
	
	my $genomeSize = 0;
	foreach my $cntg (sort keys %cntgLenHsh) {
		$genomeSize += $cntgLenHsh{$cntg};
	}
	
	my $totalPos = $genomeSize;
	
	#----loop through each pair of the comparison
	#----loop through each pair of the comparison
	my $compNum = keys %comparisonInfoHsh;
	
	foreach my $comparisonTag (keys %comparisonInfoHsh) {

		my $refSample = ${$comparisonInfoHsh{$comparisonTag}}{'refSample'};
		my $qrySample = ${$comparisonInfoHsh{$comparisonTag}}{'qrySample'};
		
		print "Start overlapping $refSample vs $qrySample\n";
		
		my $cntgProc = 0;
		my $posProc = 0;
		
		#----loop through each cutoff
		my $refSpfcCount = my $qrySpfcCount = my $bothCount = my $noneCount = 0;

		foreach my $cntg (sort {$a cmp $b} keys %cntgLenHsh) {
			$cntgProc++;
			my $qryCntgCovPlsPath = ${$cntgCovPlsPathBySamByCntgHsh{$qrySample}}{$cntg};
			system "gzip -df $qryCntgCovPlsPath" if -s "$qryCntgCovPlsPath.gz";
			my @qryCntgCovAry = @{retrieve($qryCntgCovPlsPath)};
			my $refCntgCovPlsPath = ${$cntgCovPlsPathBySamByCntgHsh{$refSample}}{$cntg};
			system "gzip -df $refCntgCovPlsPath" if -s "$refCntgCovPlsPath.gz";
			my @refCntgCovAry = @{retrieve($refCntgCovPlsPath)};
			my $cntgLen = $cntgLenHsh{$cntg};
			
			foreach my $index (0..$cntgLen) {
				
				$posProc++;
				if (not ($posProc % 10000)) {
					my $procPct = sprintf "%.2f", 100*$posProc/$totalPos;
					my $refPct = sprintf "%.2f", 100*$refSpfcCount/($posProc*@strandToCompareAry);
					my $qryPct = sprintf "%.2f", 100*$qrySpfcCount/($posProc*@strandToCompareAry);
					my $bothPct = sprintf "%.2f", 100*$bothCount/($posProc*@strandToCompareAry);
					my $nonePct = sprintf "%.2f", 100*$noneCount/($posProc*@strandToCompareAry);
					print "Processing $cntgProc of $cntgNum cntg at cov $posCovCutoff. $procPct% processed: $refSample:$refPct% $qrySample:$qryPct% both:$bothPct% none:$nonePct%      \r";
				}

				my %tmpRefCovHsh; #----a tmp hsh to store all the values
				my %tmpQryCovHsh; #----a tmp hsh to store all the values
				foreach my $strand ('+', '-', '.') {
					$tmpRefCovHsh{$strand} = 0;
					$tmpQryCovHsh{$strand} = 0;
				}
				
				#---define plusCov, minusCov and bothCov
				($tmpRefCovHsh{'+'}, $tmpRefCovHsh{'-'}) = split /,/, $refCntgCovAry[$index] if (defined $refCntgCovAry[$index]);
				($tmpQryCovHsh{'+'}, $tmpQryCovHsh{'-'}) = split /,/, $qryCntgCovAry[$index] if (defined $qryCntgCovAry[$index]);

				$tmpRefCovHsh{'.'} = $tmpRefCovHsh{'+'} + $tmpRefCovHsh{'-'};
				$tmpQryCovHsh{'.'} = $tmpQryCovHsh{'+'} + $tmpQryCovHsh{'-'};
				
				#---check if plusCov, minusCov or bothCov are sample specific
				foreach my $strand (@strandToCompareAry) {
					my $smplSpfc;

					if (($tmpRefCovHsh{$strand} >= $posCovCutoff) and ($tmpQryCovHsh{$strand} < $posCovCutoff)) {
						$smplSpfc = $refSample;
						$refSpfcCount++;
					} elsif (($tmpQryCovHsh{$strand} >= $posCovCutoff) and ($tmpRefCovHsh{$strand} < $posCovCutoff)) {
						$smplSpfc = $qrySample;
						$qrySpfcCount++;
					} elsif (($tmpRefCovHsh{$strand} >= $posCovCutoff) and ($tmpQryCovHsh{$strand} >= $posCovCutoff)) {
						$bothCount++;
						$smplSpfc = 'both';
					} elsif (($tmpRefCovHsh{$strand} < $posCovCutoff) and ($tmpQryCovHsh{$strand} < $posCovCutoff)) {
						$noneCount++;
						$smplSpfc = 'none';
					} else {
						die;
					}
					

				}#---end of foreach my $strand (@strandToCompareAry) {
			}#---end of foreach my $index (0..$#cntgLen) {
		}#---end of foreach my $cntg (sort keys %cntgLenHsh) {
		
		my $refPct = sprintf "%.2f", 100*$refSpfcCount/($totalPos*@strandToCompareAry);
		my $qryPct = sprintf "%.2f", 100*$qrySpfcCount/($totalPos*@strandToCompareAry);
		my $bothPct = sprintf "%.2f", 100*$bothCount/($totalPos*@strandToCompareAry);
		my $nonePct = sprintf "%.2f", 100*$noneCount/($totalPos*@strandToCompareAry);

		${${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}{'refPct'} = $refPct;
		${${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}{'qryPct'} = $qryPct;
		${${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}{'bothPct'} = $bothPct;
		${${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}{'nonePct'} = $nonePct;
		
		print "\n";
		
	}

	return \%comparisonInfoHsh;

}
########################################################################## 
sub detectContinousRegionFromPileupPls {

	#---detectContinousRegionFromPileupPls($minContinuousRegSize, $maxGapTorlerance, $plsIndexPath, $covCutoff);

	my ($strandToCompareAry_ref, $minContinuousRegSize, $maxGapTorlerance, $plsIndexPath, $covCutoff, $cntgLenHsh_ref) = @_;

	my %cntgLenHsh = %{$cntgLenHsh_ref};
	my @strandToCompareAry = @{$strandToCompareAry_ref};
	my $genomeSize = 0;
	foreach my $cntg (keys %cntgLenHsh) {
		$genomeSize += $cntgLenHsh{$cntg};
	}
	my $totalPosToProc = @strandToCompareAry*$genomeSize;
	
	my %continousRegionRngHsh = ();
	my %continousRegionStrandHsh = ();
	my %continousRegionCovPerNtHsh = ();
	
	my %plsIndexHsh = %{retrieve($plsIndexPath)};
	my (undef, $cntgCovStroableDir, undef) = fileparse($plsIndexPath, qr/\.[^.]*/);
	
	my %cntnusRegNumHsh = ();
	
	my $procPos = 0;

	my $totalContinusousRegNum = 0;

	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		system "gzip -df $cntgCovPlsPath" if -s "$cntgCovPlsPath.gz";
		my @cntgCovAry = @{retrieve($cntgCovPlsPath)};
		
		my %tmpStreakHsh = ();
		my %tmpCovSumHsh = ();
		foreach my $strand (@strandToCompareAry) {
			@{${$tmpStreakHsh{$strand}}{'hit'}} = ();
			@{${$tmpStreakHsh{$strand}}{'gap'}} = ();
			$tmpCovSumHsh{$strand} = 0;
		}
		
		foreach my $index (0..$#cntgCovAry) {
			my %tmpCovHsh = ();
			foreach my $strand ('+', '-') {
				$tmpCovHsh{$strand} = 0;
			}
			($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry[$index] if (defined $cntgCovAry[$index]);
			$tmpCovHsh{'.'} = $tmpCovHsh{'+'}+$tmpCovHsh{'-'};
			
			foreach my $strand (@strandToCompareAry) {
				$procPos++;
				
				if (not ($procPos % 10000)) {
					my $procPct = sprintf "%.2f", 100*$procPos/$totalPosToProc;
					print "scanning continuous region......$procPct% : $totalContinusousRegNum found        \r";
				}
				
				#----if inside hit region, tune down the dynamicCovCutoff to 1, make the "within hit gap" as 0;
				#my $dynamicCovCutoff = $covCutoff;
				# $dynamicCovCutoff = 1 if @{${$tmpStreakHsh{$strand}}{'hit'}};#---still testing
				
				if ($tmpCovHsh{$strand} >= $covCutoff) {#---in  hit region
					push @{${$tmpStreakHsh{$strand}}{'hit'}}, $index;
					$tmpCovSumHsh{$strand} += $tmpCovHsh{$strand};
					@{${$tmpStreakHsh{$strand}}{'gap'}} = ();
					
				} elsif (@{${$tmpStreakHsh{$strand}}{'hit'}}) {#----in gap region and hit is > 1
					
					push @{${$tmpStreakHsh{$strand}}{'gap'}}, $index; #---increase the gap
					
					if ((@{${$tmpStreakHsh{$strand}}{'gap'}} >= $maxGapTorlerance) or ($index == $#cntgCovAry)){#---larger than maximum gap torlerance or the last index

						if (@{${$tmpStreakHsh{$strand}}{'hit'}} >= $minContinuousRegSize) {#---larger than minimum continousus region size
							my $hitRegStart = ${${$tmpStreakHsh{$strand}}{'hit'}}[0] + 1;
							my $hitRegEnd = ${${$tmpStreakHsh{$strand}}{'hit'}}[-1] + 1;
							$cntnusRegNumHsh{$strand}++;
							my $cntnusRegID = join "_", ("reg", $totalContinusousRegNum);
							${$continousRegionRngHsh{$cntg}}{$cntnusRegID} = join ",", ($hitRegStart, $hitRegEnd);
							$continousRegionStrandHsh{$cntnusRegID} = $strand;
							$continousRegionCovPerNtHsh{$cntnusRegID} = sprintf "%.2f", $tmpCovSumHsh{$strand}/@{${$tmpStreakHsh{$strand}}{'hit'}};
							$totalContinusousRegNum++;
						}
						@{${$tmpStreakHsh{$strand}}{'hit'}} = (); #---reset the hit since gap is larger than torlerated
						$tmpCovSumHsh{$strand} = 0;
					}
				}
			}
		}
	}
	print "scanning continuous region......done                      \n";

	return \%continousRegionRngHsh, \%continousRegionStrandHsh, \%continousRegionCovPerNtHsh;
}
########################################################################## 
sub generateTransfragsForAllSamples {

	# generateTransfragsForAllSamples($minContinuousRegSize, $maxGapTorlerance, $covCutoff, $plsIndexListPath, $outDir, $cntgLenHsh_ref, $plsIndexPathHsh_ref, $extendedRngHsh_ref, $strandHsh_ref);
	my ($strandToCompareAry_ref, $minContinuousRegSize, $maxGapTorlerance, $covCutoff, $plsIndexListPath, $outDir, $cntgLenHsh_ref, $plsIndexPathHsh_ref, $extendedRngHsh_ref, $strandHsh_ref) = @_;

	my %plsIndexPathHsh = %{$plsIndexPathHsh_ref};

	my %allTransfragRngHsh = ();
	my %allTransfragStrandHsh = ();
	my %allTransfragInfoHsh = ();
	
	foreach my $sample (sort keys %plsIndexPathHsh) {
		my $plsIndexPath = $plsIndexPathHsh{$sample};
		print "Scanning transfrag for $sample\n";
		my ($transfragRngHsh_ref, $transfragStrandHsh_ref, $transfragCovPerNtHsh_ref) = detectContinousRegionFromPileupPls($strandToCompareAry_ref, $minContinuousRegSize, $maxGapTorlerance, $plsIndexPath, $covCutoff, $cntgLenHsh_ref);
		my %transfragRngHsh = %{$transfragRngHsh_ref};
		my %transfragStrandHsh = %{$transfragStrandHsh_ref};
		my %transfragCovPerNtHsh = %{$transfragCovPerNtHsh_ref};

		#---check whether the transfrag overlap with existing feature
		my ($transfragOvrlpAnnotatedHsh_ref, $annotatedTransfragStrandHsh_ref, $unannotatedTransfragStrandHsh_ref, $transfragAnnoFturHitHsh_ref) = checkOverlapWithExsistingFeature($transfragRngHsh_ref, $transfragStrandHsh_ref, $extendedRngHsh_ref, $strandHsh_ref);
		my %transfragOvrlpAnnotatedHsh = %{$transfragOvrlpAnnotatedHsh_ref};
		my %transfragAnnoFturHitHsh = %{$transfragAnnoFturHitHsh_ref};

		printGFF($transfragRngHsh_ref, $transfragStrandHsh_ref, $transfragCovPerNtHsh_ref, "$outDir/GFF/$sample.all.transfrag.gff");
		printGFF($transfragRngHsh_ref, $annotatedTransfragStrandHsh_ref, $transfragCovPerNtHsh_ref, "$outDir/GFF/$sample.anno.transfrag.gff");
		printGFF($transfragRngHsh_ref, $unannotatedTransfragStrandHsh_ref, $transfragCovPerNtHsh_ref, "$outDir/GFF/$sample.unanno.transfrag.gff");
		
		foreach my $cntg (keys %transfragRngHsh) {
			foreach my $transfragID (keys %{$transfragRngHsh{$cntg}}) {
				${${$allTransfragRngHsh{$sample}}{$cntg}}{$transfragID} = ${$transfragRngHsh{$cntg}}{$transfragID};
				${$allTransfragStrandHsh{$sample}}{$transfragID} = $transfragStrandHsh{$transfragID};
				${${$allTransfragInfoHsh{$sample}}{$transfragID}}{'ovrlpAnnotated'} = $transfragOvrlpAnnotatedHsh{$transfragID};
				${${$allTransfragInfoHsh{$sample}}{$transfragID}}{'covPerNt'} = $transfragCovPerNtHsh{$transfragID};
				${${$allTransfragInfoHsh{$sample}}{$transfragID}}{'annoFturHit'} = $transfragAnnoFturHitHsh{$transfragID};
			}
		}
	}
	
	return \%allTransfragRngHsh, \%allTransfragStrandHsh, \%allTransfragInfoHsh;
}
########################################################################## checkOverlapAndProximity
sub checkOverlapAndProximity {
#
#	dependence: printProgressScale, updateProgressBar
#
#					The 7 scenes of overlapping and proximity
#
#
# case 0: complete overlapp (($refStart == $qryStart) && ($refEnd == $qryEnd))
#
# case 1: overlapHead case 2: overlapTail	 case 3: cover		 case 4: within		case 5: prxmtyTail	 case 6: prxmtyHead
#
#ref |--------|		 |---------|	 |-------------|	 |-----|					|-----|				 			 	 |-------
#Qry 	<=========>	 <==========>		 <=========>	 <==========>						<==========>	 <==========>
#
# ($refStart<$qryStart)&&	 ($refStart>=$qryStart)&&	 ($refStart<$qryStart)&& ($refStart>$qryStart)&& ($refEnd<=$qryStart)&&		($refStart>$qryStart)&&
# ($refEnd>=$qryStart)&&	 ($refStart<=$qryEnd)&&	 ($refEnd>$qryEnd)	 ($refEnd<$qryEnd)			($refEnd<$qryEnd)			($refStart>=$qryEnd)
# ($refEnd<=$qryEnd)	 ($refEnd>$qryEnd)
#
	#---incoming variables
	my %refRngXSHsh = %{$_[0]};
	my %qryRngXSHsh = %{$_[1]};
	my %refstrandHsh = %{$_[2]};
	my %qrystrandHsh = %{$_[3]};
	my $checkPrxmty = $_[4]; #---yes or no
	my $reportExactMatch = $_[5]; #---yes or no


	my $refFturToProc = my $procRefFtur = 0;
	foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {
		foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {
			$refFturToProc++;
		}
	}

	
	#---outgoing variables
	my (%SSHitByRefHsh, %SSHitByQryHsh, %XSHitByRefHsh, %XSHitByQryHsh, %SSPrxmtyByRefHsh, %SSPrxmtyByQryHsh, %XSPrxmtyByRefHsh, %XSPrxmtyByQryHsh);
	
	
	foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {

		my (%tmpSSPrxmtyByRefHsh, %tmpSSPrxmtyByQryHsh, %tmpXSPrxmtyByRefHsh, %tmpXSPrxmtyByQryHsh);

		if ((exists $qryRngXSHsh{$cntg}) and (exists $refRngXSHsh{$cntg})) {#---if there are both ref and qry can both be found on cntg
			foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {#--- all ftur on the $strand of $cntg of refGff
				$procRefFtur++;
				print "Checking overlapping: $procRefFtur of $refFturToProc done          \r";
				
				my ($refStart, $refEnd) = split /,/, ${$refRngXSHsh{$cntg}}{$refFtur};

				foreach my $qryFtur (sort {$a cmp $b} keys %{$qryRngXSHsh{$cntg}}) {#--- all ftur on the $strand of $cntg of QryGtf

					my $samestrand = "no";
					if (($refstrandHsh{$refFtur} eq $qrystrandHsh{$qryFtur}) or ($refstrandHsh{$refFtur} eq '.') or ($qrystrandHsh{$qryFtur} eq '.')) {#---'.' for both strand
						$samestrand = "yes";
					}

					my ($qryStart, $qryEnd) = split /,/, ${$qryRngXSHsh{$cntg}}{$qryFtur};

					if (($refStart == $qryStart) && ($refEnd == $qryEnd)) {#---scene 0
						
						if ($reportExactMatch eq "yes") {
							${$XSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
							${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 0;

							if ($samestrand eq "yes") {
								${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 0;
								${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 0;
							}
						}

					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 1

						${$XSHitByRefHsh{$refFtur}}{$qryFtur} = 1;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 1;

						if ($samestrand eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 1;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 1;
						}

					} elsif (($refStart>=$qryStart)&&($refStart<=$qryEnd)&&($refEnd>=$qryEnd)) {#---scene 2

						${$XSHitByRefHsh{$refFtur}}{$qryFtur} = 2;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 2;

						if ($samestrand eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 2;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 2;
						}

					} elsif (($refStart<=$qryStart)&&($refEnd>=$qryEnd)) {#---scene 3

						${$XSHitByRefHsh{$refFtur}}{$qryFtur} = 3;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 3;

						if ($samestrand eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 3;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 3;
						}

					} elsif (($refStart>=$qryStart)&&($refEnd<=$qryEnd)) {#---scene 4

						${$XSHitByRefHsh{$refFtur}}{$qryFtur} = 4;
						${$XSHitByQryHsh{$qryFtur}}{$refFtur} = 4;

						if ($samestrand eq "yes") {
							${$SSHitByRefHsh{$refFtur}}{$qryFtur} = 4;
							${$SSHitByQryHsh{$qryFtur}}{$refFtur} = 4;
						}

					#------Proximity with ref's tail proximal to qry's head
					} elsif (($refEnd<=$qryStart)&&($refEnd<$qryEnd)) {#---scene 5 ---> ref Tail, qry Head

						if ($checkPrxmty eq "yes") {
							my $tmpPrmxty = $qryStart - $refEnd;
							${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}{$qryFtur} = $tmpPrmxty;
							${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur} = $tmpPrmxty;

							if ($samestrand eq "yes") {
								${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}{$qryFtur} = $tmpPrmxty;
								${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur} = $tmpPrmxty;
							}
						}

					#------Proximity with ref's head proximal to qry's tail
					} elsif (($refStart>=$qryEnd)&&($refStart>$qryStart)) {#---scene 6 ---> ref Head, qry Tail

						if ($checkPrxmty eq "yes") {
							my $tmpPrmxty = $refStart - $qryEnd;
							${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}{$qryFtur} = $tmpPrmxty;
							${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur} = $tmpPrmxty;

							if ($samestrand eq "yes") {
								${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}{$qryFtur} = $tmpPrmxty;
								${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur} = $tmpPrmxty;
							}
						}

					} else {#---BUG! possibly other scene?
						print "refStart=$refStart; refEnd=$refEnd; qryStart=$qryStart; qryEnd=$qryEnd\n";
						die "Unexpected overlapping scene between $refFtur and $qryFtur. It's a Bug. Program qutting.\n";
					}
				}
			}#---end of foreach my $refFtur (sort {$a cmp $b} keys %{$refRngXSHsh{$cntg}}) {#--- all ftur on the $strand of $cntg of refGff
		} #---end of if (exists $qryRngXSHsh{$cntg}) {

		#---find the closest proximity for all refs
		if ($checkPrxmty eq "yes") {

			#---for all ref based info
			foreach my $refFtur (keys %{$refRngXSHsh{$cntg}}) {

				#---in cases if the proximity are edges
				${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"});
				${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"});
				${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"});
				${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"});

				#---for all XS heads
				foreach my $qryFtur (sort {${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}{$a} <=> ${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}{$b}} keys %{${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}) {
					@{${$XSPrxmtyByRefHsh{$refFtur}}{"H"}} = (${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"H"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all XS tails
				foreach my $qryFtur (sort {${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}{$a} <=> ${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}{$b}} keys %{${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}) {
					@{${$XSPrxmtyByRefHsh{$refFtur}}{"T"}} = (${${$tmpXSPrxmtyByRefHsh{$refFtur}}{"T"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all SS heads
				foreach my $qryFtur (sort {${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}{$a} <=> ${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}{$b}} keys %{${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}) {
					@{${$SSPrxmtyByRefHsh{$refFtur}}{"H"}} = (${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"H"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

				#---for all SS tails
				foreach my $qryFtur (sort {${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}{$a} <=> ${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}{$b}} keys %{${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}) {
					@{${$SSPrxmtyByRefHsh{$refFtur}}{"T"}} = (${${$tmpSSPrxmtyByRefHsh{$refFtur}}{"T"}}{$qryFtur}, $qryFtur);
					last; #---sample the smallest only
				}

			}#---end of foreach my $refFtur (keys %refRngXSHsh)

			#---for all qry based info
			foreach my $qryFtur (keys %{$qryRngXSHsh{$cntg}}) {

				#---in cases if the proximity are edges
				${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"});
				${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"});
				${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"});
				${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{"edge"} = -999 if (not exists ${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"});

				#---for all XS heads
				foreach my $refFtur (sort {${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$a} <=> ${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$b}} keys %{${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}) {
					@{${$XSPrxmtyByQryHsh{$qryFtur}}{"H"}} = (${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all XS tails
				foreach my $refFtur (sort {${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$a} <=> ${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$b}} keys %{${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}) {
					@{${$XSPrxmtyByQryHsh{$qryFtur}}{"T"}} = (${${$tmpXSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all SS heads
				foreach my $refFtur (sort {${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$a} <=> ${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$b}} keys %{${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}) {
					@{${$SSPrxmtyByQryHsh{$qryFtur}}{"H"}} = (${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"H"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

				#---for all SS tails
				foreach my $refFtur (sort {${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$a} <=> ${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$b}} keys %{${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}) {
					@{${$SSPrxmtyByQryHsh{$qryFtur}}{"T"}} = (${${$tmpSSPrxmtyByQryHsh{$qryFtur}}{"T"}}{$refFtur}, $refFtur);
					last; #---sample the smallest only
				}

			}#---end of foreach my $qryFtur (keys %qryRngXSHsh) {
		}

	}#---end foreach my $cntg (sort {$a cmp $b} keys %refRngXSHsh) {

	print "Checking overlapping: All of $refFturToProc done                \n";

	return (\%SSHitByRefHsh, \%SSHitByQryHsh, \%XSHitByRefHsh, \%XSHitByQryHsh, \%SSPrxmtyByRefHsh, \%SSPrxmtyByQryHsh, \%XSPrxmtyByRefHsh, \%XSPrxmtyByQryHsh);
}
########################################################################## printTransfragGFF
sub printGFF {

	my %rngHsh = %{$_[0]}; 
	my %strandHsh = %{$_[1]};
	my %covPerNtHsh = %{$_[2]};
	my $outGFFPath = $_[3];

	#DS571232	ApiDB	gene	16744	17819	.	+	.	ID=EHI_074040;Name=EHI_074040;description=hypothetical+protein;size=1076;web_id=EHI_074040;locus_tag=EHI_074040;size=1076;Alias=3294.t00006,3294.m00006
	#DS571232	ApiDB	mRNA	16744	17819	.	+	.	ID=rna_EHI_074040-1;Name=EHI_074040-1;description=hypothetical+protein;size=1076;Parent=EHI_074040;Dbxref=ApiDB:EHI_074040,NCBI_gi:56470532,NCBI_gi:67475836,taxon:294381
	#DS571232	ApiDB	CDS		16744	16955	.	+	0	ID=cds_EHI_074040-1;Name=cds;description=.;size=212;Parent=rna_EHI_074040-1
	#DS571232	ApiDB	CDS		17003	17819	.	+	1	ID=cds_EHI_074040-1;Name=cds;description=.;size=817;Parent=rna_EHI_074040-1
	#DS571232	ApiDB	exon	16744	16955	.	+	.	ID=exon_EHI_074040-1;Name=exon;description=exon;size=212;Parent=rna_EHI_074040-1
	#DS571232	ApiDB	exon	17003	17819	.	+	.	ID=exon_EHI_074040-2;Name=exon;description=exon;size=817;Parent=rna_EHI_074040-1

	print "Printing $outGFFPath\n";

	open (GFFOUT, ">$outGFFPath");
	print GFFOUT "##gff-version	3\n";

	foreach my $cntg (sort {$a cmp $b} keys %rngHsh) {
		foreach my $ID (sort {$a cmp $b} keys %{$rngHsh{$cntg}}) {
			if (exists $strandHsh{$ID}) {
				my $covPerNt = 'unknown';
				$covPerNt = $covPerNtHsh{$ID} if exists $covPerNtHsh{$ID};
				my $strand = $strandHsh{$ID};
				my ($start, $end) = split /,/, ${$rngHsh{$cntg}}{$ID};
				my $length = $end - $start;
				print GFFOUT "$cntg"."\t"."BCP"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."."."\t"."$strand"."\t"."."."\t"."ID=$ID;Name=$ID.cov$covPerNt;description=$ID.cov$covPerNt;size=$length;\n";
				print GFFOUT "$cntg"."\t"."BCP"."\t"."mRNA"."\t"."$start"."\t"."$end"."\t"."."."\t"."$strand"."\t"."."."\t"."ID=rna_$ID\-1;Name=$ID.cov$covPerNt\-1;description=cov$covPerNt;size=$length;Parent=$ID;covPerNt=$covPerNt\n";
				print GFFOUT "$cntg"."\t"."BCP"."\t"."exon"."\t"."$start"."\t"."$end"."\t"."."."\t"."$strand"."\t"."."."\t"."ID=exon_$ID\-1;Name=exon;description=exon;size=$length;Parent=rna_$ID\-1;\n";
			}
		}
	}
	close GFFOUT;
}
########################################################################## readGff
sub readGff {

	#---the whole subrountine was inherited from pileupCounter_v0.8 so it looks a bit redundant. Will come back later to clean it up.

	#---variables to retun
	my (%strandHsh, %cntgByGeneHsh, %exonRngByGeneHsh, %exonNumByCntgHsh, %geneExonLenHsh, %geneCDSLenHsh, %ctgryReadCountHsh, %CDSRngByGeneHsh, %geneByCtgryHsh, %ctgryByGeneHsh, %descriptionByGeneHsh, %exonLengthByGeneHsh);
	my (%geneByRNAHsh, %CDSCountHsh, %exonCountHsh, %geneExonLocationHsh);
	my (%SSRngByCntgByWholeGeneHsh, %XstrandRngByCntgByWholeGeneHsh);

	#---read the gff
	my $gffPath = $_[0];

	my @gffPathSplt = split /\//, $gffPath;
	my $gffFileName = $gffPathSplt[-1];
	open (INFILE, $gffPath) || die "Cannot open $gffFileName";
	print "Reading $gffFileName.\n";
	while (my $theLine = <INFILE>) {
		chomp $theLine;
		if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
			my @theLineSplt = split (/\t/, $theLine);
			my $cntg = $theLineSplt[0];

			my $geneCategory = $theLineSplt[2];
			my $featureStart = $theLineSplt[3];
			my $featureEnd = $theLineSplt[4];
			my $strand = $theLineSplt[6];
			my $attribute = $theLineSplt[8];
			my @attributeSplt = split /;/, $attribute;
			my ($unqID, $parent, $description);
			foreach my $theAttribute (@attributeSplt) {
				if ($theAttribute =~ m/^ID=/) {$unqID = substr ($theAttribute, index ($theAttribute, "=")+1);}
				if ($theAttribute =~ m/^Parent=/) {$parent = substr ($theAttribute, index ($theAttribute, "=")+1);}
				if ($theAttribute =~ m/^description=/) {$description = substr ($theAttribute, index ($theAttribute, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene

				my $geneID = $unqID;
				$strandHsh{$geneID} = $strand;
				$description = "unknown" if (not defined $description);
				$descriptionByGeneHsh{$geneID} = $description;
				$cntgByGeneHsh{$geneID} = $cntg;
				${${${$SSRngByCntgByWholeGeneHsh{$cntg}}{$strand}}{$geneID}}{"start"} = $featureStart;
				${${${$SSRngByCntgByWholeGeneHsh{$cntg}}{$strand}}{$geneID}}{"end"} = $featureEnd;

				${${$XstrandRngByCntgByWholeGeneHsh{$cntg}}{$geneID}}{"start"} = $featureStart;
				${${$XstrandRngByCntgByWholeGeneHsh{$cntg}}{$geneID}}{"end"} = $featureEnd;

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes

				# The CDS is ignored at the moment, until it reaches the point that we are looking at UTRs
				#
				my $mRNAID = $parent;
				my $geneID = $geneByRNAHsh{$mRNAID};
				$CDSCountHsh{$geneID}++;
				my $CDSCount = $CDSCountHsh{$geneID};
				${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"} = $featureStart;
				${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} = $featureEnd;
			 	$geneCDSLenHsh{$geneID} = 0 if $CDSCount == 1; #---define the length hashfor the 1st time
			 	$geneCDSLenHsh{$geneID} += ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"end"} - ${${$CDSRngByGeneHsh{$geneID}}{$CDSCount}}{"start"};

			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later

				my $exonID = $unqID;
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh{$RNAID};
				my $locationTag = $cntg.":".$featureStart.":".$featureEnd;

				${$geneExonLocationHsh{$locationTag}}{$geneID}++;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"} = $featureStart;
				${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} = $featureEnd;
				my $exonNum = keys %{$exonRngByGeneHsh{$geneID}};
			 	$exonLengthByGeneHsh{$geneID} = 0 if $exonNum == 1; #---define the length hashfor the 1st time
			 	$exonLengthByGeneHsh{$geneID} += ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"end"} - ${${$exonRngByGeneHsh{$geneID}}{$exonID}}{"start"};

			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh{$RNAID} = $geneID;
				$ctgryByGeneHsh{$geneID} = $geneCategory;
				$geneByCtgryHsh{$geneCategory} = $geneID;

				if (not(exists $ctgryReadCountHsh{$geneCategory})) {#---initialize the $geneCategory for all category
					${$ctgryReadCountHsh{$geneCategory}}{"s"} = 0;
					${$ctgryReadCountHsh{$geneCategory}}{"a"} = 0;
				}
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close INFILE;

	my ($intronRngByGeneHsh_ref, $junctStrByIntronIDHsh_ref, $intronRngXS_ref, $strandByIntronHsh_ref, $strandByJunctStrHsh_ref) = getIntronFromExonRng(\%exonRngByGeneHsh, \%cntgByGeneHsh, \%strandHsh, "Getting the intron boundaries on Gff", "ref");

	my $filteredGeneNum = keys %strandHsh;

	print "Totally $filteredGeneNum gene have been stored.\n";

	return (\%SSRngByCntgByWholeGeneHsh, \%ctgryByGeneHsh, \%strandHsh, \%XstrandRngByCntgByWholeGeneHsh, \%exonRngByGeneHsh, \%CDSRngByGeneHsh, $intronRngByGeneHsh_ref, $junctStrByIntronIDHsh_ref, $intronRngXS_ref, $strandByIntronHsh_ref, \%cntgByGeneHsh, $strandByJunctStrHsh_ref, \%descriptionByGeneHsh, \%exonLengthByGeneHsh);

	#################################################### getIntronFromExonRng ############################################################################
	sub getIntronFromExonRng {

		#---incoming variables
		my %exonRngByGeneHsh = %{$_[0]};
		my %cntgByGeneHsh = %{$_[1]};
		my %strandHsh = %{$_[2]};
		my $strToPrint = $_[3];
		my $refOrNGS = $_[4];

		#---outgoing variables
		my (%intronRngByGeneHsh, %junctStrByIntronIDHsh, %XstrandRngByCntgByIntronID, %strandByIntronIDHsh, %strandByJunctStrHsh);

		#---on screen progress scale
		my $totalGeneNum = keys %exonRngByGeneHsh;

		#---go through each gene, see if there's intron, and store the ranges
		foreach my $geneID (sort {$a cmp $b} keys %exonRngByGeneHsh) {


			#---get the cntg and strand
			my $cntg = $cntgByGeneHsh{$geneID};
			my $strand = $strandHsh{$geneID};

			#---get all exon bounds
			my @tmpBoundAry;
			foreach my $exon (keys %{$exonRngByGeneHsh{$geneID}}) {
				push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"start"};
				push @tmpBoundAry, ${${$exonRngByGeneHsh{$geneID}}{$exon}}{"end"};
			}

			#---if more than one exon
			if (@tmpBoundAry > 2) {
				#---sort the bounds
				my @sortedTmpBoundAry = sort {$a <=> $b} @tmpBoundAry;
				#---get the intronBounds, go to all odd number indexes and take itself and the +1 index values, i.e. $i = 1, 3, 5 if there're 0,1,2,3,4,5,6,7
				my $intronNum = 0;
				for (my $i=1; $i<(@sortedTmpBoundAry-1); $i=$i+2) {
					$intronNum++;
					my $intronStart = $sortedTmpBoundAry[$i]+1;
					my $intronEnd = $sortedTmpBoundAry[$i+1]-1;
					my $intronID = $geneID.":".$intronNum;

					#---determine the junct is on ref, prominent or minor isoform
					my $prominentMinorRef;
					if ($refOrNGS eq "ref") {
						$prominentMinorRef = "R";
					} else {
						my @geneIDSplt = split /\./, $geneID;
						my $isoformNum = $geneIDSplt[-1];
						if ($isoformNum == 0) {#---prominent isoform
							$prominentMinorRef = "P";
						} else {
							$prominentMinorRef = "M";
						}
					}
					my $junctStr = join ":", ($cntg, $intronStart, $intronEnd);
					${${$intronRngByGeneHsh{$geneID}}{$intronNum}}{"start"} = $intronStart;
					${${$intronRngByGeneHsh{$geneID}}{$intronNum}}{"end"} = $intronEnd;
					${${$XstrandRngByCntgByIntronID{$cntg}}{$intronID}}{"start"} = $intronStart;
					${${$XstrandRngByCntgByIntronID{$cntg}}{$intronID}}{"end"} = $intronEnd;
					$strandByIntronIDHsh{$intronID} = $strand;
					$junctStrByIntronIDHsh{$intronID} = $junctStr;
					$strandByJunctStrHsh{$junctStr} = $strand;
				}
			}
		}

		return (\%intronRngByGeneHsh, \%junctStrByIntronIDHsh, \%XstrandRngByCntgByIntronID, \%strandByIntronIDHsh, \%strandByJunctStrHsh)
	}
}
########################################################################## generateExtendedKnownFeatureRng
sub generateExtendedKnownFeatureRng {
	
	#---my $extendedRngHsh_ref = generateExtendedKnownFeatureRng($XstrandRngByCntgByWholeGeneHsh_ref, $extendRange);
	
	my ($XstrandRngByCntgByWholeGeneHsh_ref, $extendRange, $ctgryByGeneHsh_ref) = @_;
	
	my %XstrandRngByCntgByWholeGeneHsh = %{$XstrandRngByCntgByWholeGeneHsh_ref};
	my %ctgryByGeneHsh = %{$ctgryByGeneHsh_ref};
	
	my %extendedRngHsh = ();
	
	foreach my $cntg (keys %XstrandRngByCntgByWholeGeneHsh) {
		foreach my $geneID (keys %{$XstrandRngByCntgByWholeGeneHsh{$cntg}}) {
			my $ctgry = $ctgryByGeneHsh{$geneID};
			if ($ctgry =~ m/mRNA|rRNA|tRNA|snRNA|snoRNA|pseudogenic_transcript/) {
				my $start = ${${$XstrandRngByCntgByWholeGeneHsh{$cntg}}{$geneID}}{"start"} - $extendRange;
				my $end = ${${$XstrandRngByCntgByWholeGeneHsh{$cntg}}{$geneID}}{"end"} + $extendRange;
				${$extendedRngHsh{$cntg}}{$geneID} = join ",", ($start, $end);
			}
		}
	}
	
	return \%extendedRngHsh;
}
########################################################################## checkOverlapWithExsistingFeature
sub checkOverlapWithExsistingFeature {
	
	#---my $transfragOvrlpAnnotatedHsh_ref = checkOverlapWithExsistingFeature($transfragRngHsh_ref, $transfragStrandHsh_ref, $extendedRngHsh_ref, $strandHsh_ref);

	my ($transfragRngHsh_ref, $transfragStrandHsh_ref, $extendedRngHsh_ref, $strandHsh_ref) = @_;
	
	my ($SSHitByTransfragHsh_ref, undef, undef, undef, undef, undef, undef, undef) = checkOverlapAndProximity($transfragRngHsh_ref, $extendedRngHsh_ref, $transfragStrandHsh_ref, $strandHsh_ref, 'no', 'no');
	
	my %SSHitByTransfragHsh = %{$SSHitByTransfragHsh_ref};
	
	my %transfragOvrlpAnnotatedHsh = ();
	
	my %transfragStrandHsh = %{$transfragStrandHsh_ref};
	my %annotatedTransfragStrandHsh = ();
	my %unannotatedTransfragStrandHsh = ();
	my %transfragAnnoFturHitHsh = ();
	
	my $annotatedCount = my $unannotatedCount = 0;
	
	foreach my $transfragID (keys %transfragStrandHsh) {
		my $ovrlpAnnotated;
		my @annoFturHitAry = qw /null/;
		if (exists $SSHitByTransfragHsh{$transfragID}) {
			@annoFturHitAry = ();
			foreach my $annoFturHit (keys %{$SSHitByTransfragHsh{$transfragID}}) {
				push @annoFturHitAry, $annoFturHit;
			}
			$ovrlpAnnotated = 'annotated' ;
			$annotatedCount++;
			$annotatedTransfragStrandHsh{$transfragID} = $transfragStrandHsh{$transfragID};
		} else {
			$ovrlpAnnotated = "unknown";
			$unannotatedCount++;
			$unannotatedTransfragStrandHsh{$transfragID} = $transfragStrandHsh{$transfragID};
		}
		$transfragOvrlpAnnotatedHsh{$transfragID} = $ovrlpAnnotated;
		$transfragAnnoFturHitHsh{$transfragID} = join ",", @annoFturHitAry;
		print "annotated:$annotatedCount vs unannotated:$unannotatedCount      \r";
	}
	
	print "annotated:$annotatedCount vs unannotated:$unannotatedCount      \n";
	
	return \%transfragOvrlpAnnotatedHsh, \%annotatedTransfragStrandHsh, \%unannotatedTransfragStrandHsh, \%transfragAnnoFturHitHsh;
}
########################################################################## compareCntgCovPlsSamplePair
sub compareTransfragSamplePair {

	#---compareTransfragSamplePair($comparisonInfoHsh_ref, $cntgLenHsh_ref, $outDir, $allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, $allTransfragInfoHsh_ref);
	my ($comparisonInfoHsh_ref, $cntgLenHsh_ref, $outDir, $allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, $allTransfragInfoHsh_ref) = @_;
	
	my %comparisonInfoHsh = %{$comparisonInfoHsh_ref};
	my %cntgLenHsh = %{$cntgLenHsh_ref};
	my %allTransfragRngHsh = %{$allTransfragRngHsh_ref};
	my %allTransfragStrandHsh = %{$allTransfragStrandHsh_ref};
	my %allTransfragInfoHsh = %{$allTransfragInfoHsh_ref};

	#----loop through each pair of the comparison
	my $compNum = keys %comparisonInfoHsh;
	
	foreach my $comparisonTag (keys %comparisonInfoHsh) {

		my $refSample = ${$comparisonInfoHsh{$comparisonTag}}{'refSample'};
		my $qrySample = ${$comparisonInfoHsh{$comparisonTag}}{'qrySample'};
		my %tmpTransfragTypeCountHsh = ();
		
		#----get the transfrag ranges
		my %refTransfragRngHsh = ();
		my %qryTransfragRngHsh = ();
		my %refTransfragStrandHsh = ();
		my %qryTransfragStrandHsh = ();
		my %refTransfragOvrlpAnnotatedHsh = ();
		my %qryTransfragOvrlpAnnotatedHsh = ();
		my %refTransfragCovPerNtHsh = ();
		my %qryTransfragCovPerNtHsh = ();
		
		foreach my $cntg (keys %cntgLenHsh) {
			if (exists ${$allTransfragRngHsh{$refSample}}{$cntg}) {
				foreach my $transfragID (keys %{${$allTransfragRngHsh{$refSample}}{$cntg}}) {
					${$refTransfragRngHsh{$cntg}}{$transfragID} = ${${$allTransfragRngHsh{$refSample}}{$cntg}}{$transfragID};
					$refTransfragStrandHsh{$transfragID} = ${$allTransfragStrandHsh{$refSample}}{$transfragID};
					$refTransfragOvrlpAnnotatedHsh{$transfragID} = ${${$allTransfragInfoHsh{$refSample}}{$transfragID}}{'ovrlpAnnotated'};
					$refTransfragCovPerNtHsh{$transfragID} = ${${$allTransfragInfoHsh{$refSample}}{$transfragID}}{'covPerNt'};
				}
			}
			
			if (exists ${$allTransfragRngHsh{$qrySample}}{$cntg}) {
				foreach my $transfragID (keys %{${$allTransfragRngHsh{$qrySample}}{$cntg}}) {
					${$qryTransfragRngHsh{$cntg}}{$transfragID} = ${${$allTransfragRngHsh{$qrySample}}{$cntg}}{$transfragID};
					$qryTransfragStrandHsh{$transfragID} = ${$allTransfragStrandHsh{$qrySample}}{$transfragID};
					$qryTransfragOvrlpAnnotatedHsh{$transfragID} = ${${$allTransfragInfoHsh{$qrySample}}{$transfragID}}{'ovrlpAnnotated'};
					$qryTransfragCovPerNtHsh{$transfragID} = ${${$allTransfragInfoHsh{$qrySample}}{$transfragID}}{'covPerNt'};
				}
			}
		}
		
		#--check the overlapping
		my ($SSHitByRefHsh_ref, $SSHitByQryHsh_ref, undef, undef, undef, undef, undef, undef) = checkOverlapAndProximity(\%refTransfragRngHsh, \%qryTransfragRngHsh, \%refTransfragStrandHsh, \%qryTransfragStrandHsh, 'no', 'yes');
		my %SSHitByRefHsh = %{$SSHitByRefHsh_ref};
		my %SSHitByQryHsh = %{$SSHitByQryHsh_ref};
		
		#----assign the overlapping of transfrag
		my %refCommonTransfragStrandHsh = ();
		my %refUniqueTransfragStrandHsh = ();
		my %qryCommonTransfragStrandHsh = ();
		my %qryUniqueTransfragStrandHsh = ();
		
		foreach my $cntg (keys %cntgLenHsh) {
			if (exists ${$allTransfragRngHsh{$refSample}}{$cntg}) {
				foreach my $transfragID (keys %{${$allTransfragRngHsh{$refSample}}{$cntg}}) {
					my $sampleSpecifc;
					if (exists $SSHitByRefHsh{$transfragID}) {
						$sampleSpecifc = 'common';
						$refCommonTransfragStrandHsh{$transfragID} = $refTransfragStrandHsh{$transfragID};
					} else {
						$sampleSpecifc = 'unique';
						$refUniqueTransfragStrandHsh{$transfragID} = $refTransfragStrandHsh{$transfragID};
					}
					
					my $ovrlpAnnotated = $refTransfragOvrlpAnnotatedHsh{$transfragID};
					${${$tmpTransfragTypeCountHsh{'ref'}}{$ovrlpAnnotated}}{$sampleSpecifc}++;
					${${$allTransfragInfoHsh{$refSample}}{$transfragID}}{"vs.$qrySample"} = $sampleSpecifc;
				}
			}
			
			if (exists ${$allTransfragRngHsh{$qrySample}}{$cntg}) {
				foreach my $transfragID (keys %{${$allTransfragRngHsh{$qrySample}}{$cntg}}) {
					my $sampleSpecifc;
					if (exists $SSHitByQryHsh{$transfragID}) {
						$sampleSpecifc = 'common';
						$qryCommonTransfragStrandHsh{$transfragID} = $qryTransfragStrandHsh{$transfragID};
					} else {
						$sampleSpecifc = 'unique';
						$qryUniqueTransfragStrandHsh{$transfragID} = $qryTransfragStrandHsh{$transfragID};
					}

					my $ovrlpAnnotated = $qryTransfragOvrlpAnnotatedHsh{$transfragID};
					${${$tmpTransfragTypeCountHsh{'qry'}}{$ovrlpAnnotated}}{$sampleSpecifc}++;
					${${$allTransfragInfoHsh{$qrySample}}{$transfragID}}{"vs.$refSample"} = $sampleSpecifc;
				}
			}
		}
		
		#---gather the transfrag type info
		my %totalCountHsh = ();
		$totalCountHsh{'ref'} = keys %refTransfragStrandHsh;
		$totalCountHsh{'qry'} = keys %qryTransfragStrandHsh;
		foreach my $refQry ('ref', 'qry') {
			foreach my $ovrlpAnnotated ('annotated', 'unknown') {
				foreach my $sampleSpecifc ('common', 'unique') {	
					my $count = 0;
					$count = ${${$tmpTransfragTypeCountHsh{$refQry}}{$ovrlpAnnotated}}{$sampleSpecifc} if exists ${${$tmpTransfragTypeCountHsh{$refQry}}{$ovrlpAnnotated}}{$sampleSpecifc};
					my $pct = sprintf "%.2f", 100*$count/$totalCountHsh{$refQry};
					my $tag = $refQry.".".$ovrlpAnnotated.".".$sampleSpecifc;
					my $value = $count."[$pct]";
					${${$comparisonInfoHsh{$comparisonTag}}{'transfragOvrlap'}}{$tag} = $value;
				}
			}
		}
		
		#----print the GFF
		printGFF(\%refTransfragRngHsh, \%refCommonTransfragStrandHsh, \%refTransfragCovPerNtHsh, "$outDir/GFF/$refSample.vs.$qrySample.common.transfrag.gff");
		printGFF(\%refTransfragRngHsh, \%refUniqueTransfragStrandHsh, \%refTransfragCovPerNtHsh, "$outDir/GFF/$refSample.vs.$qrySample.unique.transfrag.gff");
		printGFF(\%qryTransfragRngHsh, \%qryCommonTransfragStrandHsh, \%qryTransfragCovPerNtHsh, "$outDir/GFF/$qrySample.vs.$refSample.common.transfrag.gff");
		printGFF(\%qryTransfragRngHsh, \%qryUniqueTransfragStrandHsh, \%qryTransfragCovPerNtHsh, "$outDir/GFF/$qrySample.vs.$refSample.unique.transfrag.gff");
		
	}
	
	return \%allTransfragInfoHsh, \%comparisonInfoHsh;
}
########################################################################## printAllTransfragInfo
sub printAllTransfragInfo {
	
	my ($allTransfragRngHsh_ref, $allTransfragStrandHsh_ref, $allTransfragInfoHsh_ref, $outDir) = @_;

	print "Printing GFF info\n";
	my %allTransfragRngHsh = %{$allTransfragRngHsh_ref};
	my %allTransfragStrandHsh = %{$allTransfragStrandHsh_ref};
	my %allTransfragInfoHsh = %{$allTransfragInfoHsh_ref};

	foreach my $sample (keys %allTransfragInfoHsh) {
		open TRNFRGINFO, ">$outDir/$sample.GFF.info.log.tsv";
		
		#---print header
		foreach my $cntg (sort keys %{$allTransfragRngHsh{$sample}}) {
			foreach my $transfragID (sort keys %{${$allTransfragRngHsh{$sample}}{$cntg}}) {
				my @outputAry = ();
				push @outputAry, ('cntg', 'start', 'end', 'strand');
				foreach my $info (sort keys %{${$allTransfragInfoHsh{$sample}}{$transfragID}}) {
					push @outputAry, $info;
				}
				print TRNFRGINFO join "", ((join "\t", @outputAry), "\n");
				last;
			}
			last;
		}
		
		#---print content
		foreach my $cntg (sort keys %{$allTransfragRngHsh{$sample}}) {
			foreach my $transfragID (sort keys %{${$allTransfragRngHsh{$sample}}{$cntg}}) {
				my @outputAry = ();
				my ($start, $end) = split /,/, ${${$allTransfragRngHsh{$sample}}{$cntg}}{$transfragID};
				my $strand = ${$allTransfragStrandHsh{$sample}}{$transfragID};
				push @outputAry, ($cntg, $start, $end, $strand);
				foreach my $info (sort keys %{${$allTransfragInfoHsh{$sample}}{$transfragID}}) {
					push @outputAry, ${${$allTransfragInfoHsh{$sample}}{$transfragID}}{$info};
				}
				print TRNFRGINFO join "", ((join "\t", @outputAry), "\n");
			}
		}
	
		close TRNFRGINFO;
	}
	
}
########################################################################## printComparisonInfoHsh
sub printComparisonInfoHsh {
	
	#printComparisonInfoHsh(\%comparisonInfoHsh, $outDir);
	my ($comparisonInfoHsh_ref, $outDir) = @_;
	
	my $comparisonInfoHsh = %{$comparisonInfoHsh_ref};
	
	open COMLOG, ">$outDir/comparison.pairs.log.tsv";
	#----print header
	foreach my $comparisonTag (keys %comparisonInfoHsh){
		my @outputAry = ();
		push @outputAry, ('comparisonTag', 'refSample', 'qrySample');
		foreach my $posOvrlap (sort keys %{${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}){
			push @outputAry, $posOvrlap;
		}
		foreach my $transfragOvrlap (sort keys %{${$comparisonInfoHsh{$comparisonTag}}{'transfragOvrlap'}}){
			push @outputAry, $transfragOvrlap;
		}
		print COMLOG join "", ((join "\t", @outputAry), "\n");
		last;
	}

	foreach my $comparisonTag (keys %comparisonInfoHsh){
		my @outputAry = ();
		my $refSample = ${$comparisonInfoHsh{$comparisonTag}}{'refSample'};
		my $qrySample = ${$comparisonInfoHsh{$comparisonTag}}{'qrySample'};
		push @outputAry, ($comparisonTag, $refSample, $qrySample);
		foreach my $posOvrlap (sort keys %{${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}){
			push @outputAry, ${${$comparisonInfoHsh{$comparisonTag}}{'posOvrlap'}}{$posOvrlap};
		}
		foreach my $transfragOvrlap (sort keys %{${$comparisonInfoHsh{$comparisonTag}}{'transfragOvrlap'}}){
			push @outputAry, ${${$comparisonInfoHsh{$comparisonTag}}{'transfragOvrlap'}}{$transfragOvrlap};
		}

		print COMLOG join "", ((join "\t", @outputAry), "\n");
	}
	close COMLOG;

}
