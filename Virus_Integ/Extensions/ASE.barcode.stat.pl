#!/usr/bin/perl -w
use strict;
use File::Spec qw/catfile/;
use FindBin qw/$RealBin $RealScript/;
use Getopt::Long;
use List::Util qw[min max sum];
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr trible_run_for_success/;

my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);

$MODULE_NAME = 'ase_barc';
#----- version --------
$VERSION = "0.01";
$DATE = '2017-06-23';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


my ($Bam);
my ($PosToCheckFile);
my ($OutDir);
my ($samtools);
my ($readEdgeFilter) = (5);
my ($Help);

GetOptions(
	# files
	"-bam=s"	=>	\$Bam,
	"-pcf=s"	=>	\$PosToCheckFile,
	"-odir=s"	=>	\$OutDir,
	# tools
	"-samtl=s"	=>	\$samtools,
	# options
	"-rdeg=i"	=>	\$readEdgeFilter,
	# help
	"-h|help"	=>	\$Help
);

if(    !$Bam || !-e $Bam
	|| !$PosToCheckFile || !-e $PosToCheckFile
	|| !$OutDir || !-d $OutDir
	|| !$samtools || !-e $samtools
){
	warn_and_exit &return_HELP_INFO;
}

# read pos information
my %PosInfo;
&load_pos_to_check( { PosInfo_Href => \%PosInfo } );

# extract bam line for pos
&extract_bam_line( { PosInfo_Href => \%PosInfo } );

# analysis aligned reads for each position for barcodes
&stat_pos_barcode( { PosInfo_Href => \%PosInfo } );

# generate result
&generate_result( { PosInfo_Href => \%PosInfo } );

#--- return HELP_INFO
sub return_HELP_INFO{

 return "
	 Usage:   perl $RealScript <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -bam    [s]  Bam alignment file. <required>
	    -pcf    [s]  file stores the position to check. <required>
	                 format: chr  pos  ref_allele  alt_allele
	    -odir   [s]  folder to store results. <required>

	   # Tools #
	    -samtl  [s]  SamTools path. <required>

	   # Options #
	    -rdeg   [i]  Ignore base locates <= N bp away from read edges. [5]

	    -h|help      Display this help info.

 	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}

#--- read pos file ---
sub load_pos_to_check{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $PosInfo_Href = $Option_Href->{PosInfo_Href};

	open (PIF, Try_GZ_Read($PosToCheckFile)) || die "fail read PosToCheckFile: $!\n";
	while(<PIF>){
		my ($chr, $pos, $ref_base, $alt_base) = (split)[0..3];
		$PosInfo_Href->{$chr}->{$pos} = { ref_base=>uc($ref_base),
										  alt_base=>uc($alt_base),
										  out_pref=>"$chr.p$pos"     };
	}
	close PIF;

	# inform
	print `date`."[INFO]:\tLoad position to check ok.\n";
}

#--- extract bam line for each position ---
sub extract_bam_line{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $PosInfo_Href = $Option_Href->{PosInfo_Href};

	for my $chr (sort keys %$PosInfo_Href){
		for my $pos (sort {$a<=>$b} keys %{$PosInfo_Href->{$chr}}){
			my $outPref = $PosInfo_Href->{$chr}->{$pos}->{out_pref};
			my $subBam = File::Spec->catfile($OutDir, "$outPref.bam");
			# -F 0xD00 = 0x800 (supplementary) + 0x100 (secondary) + 0x400 (duplicate)
			my $command = "$samtools view -h -F 0xD00 $Bam $chr:$pos-$pos | $samtools view -b -o $subBam";
			trible_run_for_success($command, 'GetSubBam', {esdo_Nvb=>1});
			# record
			$PosInfo_Href->{$chr}->{$pos}->{sub_bam} = $subBam;
		}
	}

	# inform
	print `date`."[INFO]:\tExtract bam line ok.\n";
}

#--- read bam line for barcode support of each position ---
sub stat_pos_barcode{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $PosInfo_Href = $Option_Href->{PosInfo_Href};

	for my $chr (sort keys %$PosInfo_Href){
		for my $pos (sort {$a<=>$b} keys %{$PosInfo_Href->{$chr}}){
			my $subBam = $PosInfo_Href->{$chr}->{$pos}->{sub_bam};
			# read bam
			open (SBAM, "$samtools view -F 0xD00 $subBam |") || die "fail read subBam: $!\n";
			while(<SBAM>){
				chomp;
				my ($PID, $FLAG, $mSeq, $mPos, $CIGAR, $ReadSeq) = (split /\t+/)[0, 1, 2, 3, 5, 9];
				# filter when have other cigar operators
				next if( $CIGAR =~ /[^IDSM\d]/ ); 
				# get barcode
				## instance: E00492:112:HFVWVALXX:7:2212:13210:54506:BC:GTACTCCCACCGGTCA-1
				my ($barcode) = ( $PID =~ /:BC:([^:]+)$/ );
				next if( !defined($barcode) || $barcode =~ /NA/i );
				# store basic info
				my $readLen = length($ReadSeq);
				my $m_Orit = ($FLAG & 0x10) ? 'R' : 'F';
				# deal mutations from CIGAR
				my $i = 0;
				while(length($CIGAR)!=0){
					my ($len, $type) = ($CIGAR =~ /^(\d+)(\D)/);
					$CIGAR =~ s/^\d+\D//;
					# mapped with possible mismatch ONLY
					if($type eq 'M'){
						my @deal_seq = split //, substr($ReadSeq, 0, $len);
						if( $pos <= $mPos + $i + $len - 1 ){
							my $shift = $pos - $mPos - $i;
							# read edge control
							if(    $readLen - length($ReadSeq) + $shift + 1 <= $readEdgeFilter # 5 primer
								|| length($ReadSeq) - $shift <= $readEdgeFilter # 3 primer
							){
								next;
							}
							my $allele = uc($deal_seq[$shift]);
							if( $allele eq $PosInfo_Href->{$chr}->{$pos}->{ref_base} ){
								$PosInfo_Href->{$chr}->{$pos}->{ref_supp}->{$m_Orit}->{count} ++;
								$PosInfo_Href->{$chr}->{$pos}->{ref_supp}->{$m_Orit}->{barc}->{$barcode} = 1;
							}
							elsif( $allele eq $PosInfo_Href->{$chr}->{$pos}->{alt_base} ){
								$PosInfo_Href->{$chr}->{$pos}->{alt_supp}->{$m_Orit}->{count} ++;
								$PosInfo_Href->{$chr}->{$pos}->{alt_supp}->{$m_Orit}->{barc}->{$barcode} = 1;
							}
						}
						# update for next block
						$i += $len;
						$ReadSeq = substr($ReadSeq, $len);
					}
					elsif(   $type eq 'I'
						  || $type eq 'S'
					){
						# update ReadSeq for next block
						$ReadSeq = substr($ReadSeq, $len);
					}
					elsif($type eq 'D'){
						# update i for next block
						$i += $len;
					}
				}
			}
			close SBAM;
			# inform
			print `date`."[INFO]:\tAnalysis Bam of $chr:$pos ok.\n";
		}
	}
}

#--- read bam line for barcode support of each position ---
sub generate_result{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $PosInfo_Href = $Option_Href->{PosInfo_Href};

	for my $chr (sort keys %$PosInfo_Href){
		for my $pos (sort {$a<=>$b} keys %{$PosInfo_Href->{$chr}}){
			my $outPref = $PosInfo_Href->{$chr}->{$pos}->{out_pref};
			my $outStat = File::Spec->catfile($OutDir, "$outPref.barc_list");
			open (BARC, Try_GZ_Write($outStat)) || die "fail write barcode outStat: $!\n";
			for my $type ('ref', 'alt'){
				my $base_key = $type.'_base';
				my $supp_key = $type.'_supp';
				my $base = $PosInfo_Href->{$chr}->{$pos}->{$base_key};
				for my $m_Orit ('R', 'F'){
					print BARC "$_\t$type\t$m_Orit\n" for sort keys %{$PosInfo_Href->{$chr}->{$pos}->{$supp_key}->{$m_Orit}->{barc}};
				}
			}
			close BARC;
		}
	}

	# inform
	print `date`."[INFO]:\tGenerate barcode stat files ok.\n";
}
