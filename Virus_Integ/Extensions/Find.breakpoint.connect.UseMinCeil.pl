#!/usr/bin/perl -w
use strict;
use File::Spec qw/catfile/;
use FindBin qw/$RealBin $RealScript/;
use Getopt::Long;
use List::Util qw[min max sum];
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr trible_run_for_success/;

my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);

$MODULE_NAME = 'bK_connect';
#----- version --------
$VERSION = "0.03";
$DATE = '2017-06-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


my ($LocalHap_txt);
my ($BK_1,$BK_2);
my ($barcode_SpanStat_file);
my ($segment_length_file);
my ($max_allow_span_size) = (1E5);
my ($Help);

GetOptions(
	# files
	"-lhf=s"	=>	\$LocalHap_txt,
	"-bk1=s"	=>	\$BK_1,
	"-bk2=s"	=>	\$BK_2,
	"-bct=s"	=>	\$barcode_SpanStat_file,
	"-slf=s"	=>	\$segment_length_file,
	"-msp=i"	=>	\$max_allow_span_size,
	# help
	"-h|help"	=>	\$Help
);

if(    !$LocalHap_txt || !-e $LocalHap_txt
	|| !$BK_1
	|| !$BK_2
	|| !$barcode_SpanStat_file || !-e $barcode_SpanStat_file
	|| !$segment_length_file || !-e $segment_length_file
){
	warn_and_exit &return_HELP_INFO;
}

# read barcode span stat
my %BarcSpanStat;
my $all_count = 
&load_barcode_span( { BarcSpanStat_Href => \%BarcSpanStat } );

# load segment length
my %SegmentLen;
&load_segment_length( { SegmentLen_Href => \%SegmentLen } );

# read haplotype structure
my @HapTypeStr;
&load_haplotype_structure( { HapTypeStr_Aref => \@HapTypeStr } );

# find all connection
# calculate connections' length and probability
my %Distance;
&find_all_connection( { 
						HapTypeStr_Aref => \@HapTypeStr,
						Distance_Href => \%Distance,
						SegmentLen_Href => \%SegmentLen,
					} );

# calculate each distance probability
&output_distance_probability({
								BarcSpanStat_Href => \%BarcSpanStat,
								Distance_Href => \%Distance
							});


#--- return HELP_INFO
sub return_HELP_INFO{

 return "
	 Usage:   perl $RealScript <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -lhf    [s]  File stores localHaplotype structure. <required>
	    -bk1    [s]  Breakpoint NO.1 <required>
	    -bk2    [s]  Breakpoint NO.2 <required>
	    -bct    [s]  Stat file of barcode span size. <required>
	    -slf    [s]  Segment length file. <required>
	    -msp    [i]  The maximum allowed span size of barcode. [1E5]

	    -h|help      Display this help info.

 	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}

#--- read barcode span stat file---
sub load_barcode_span{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $BarcSpanStat_Href = $Option_Href->{BarcSpanStat_Href};

	my $all_count = 0;
	open (BCS, Try_GZ_Read($barcode_SpanStat_file)) || die "fail read barcode_SpanStat_file: $!\n";
	while(<BCS>){
		my ($count, $span) = (split)[0,1];
		next if($span == 151);
		$BarcSpanStat_Href->{$span} = $count;
		$all_count += $count;
	}
	close BCS;

	# inform
	warn `date`."[INFO]:\tLoad barcode span ok.\n";

	return $all_count;
}

#--- read haplotype structure file ---
sub load_haplotype_structure{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $HapTypeStr_Aref = $Option_Href->{HapTypeStr_Aref};

	my $hap_NO = 0;
	open (HLT, Try_GZ_Read($LocalHap_txt)) || die "fail read LocalHap_txt: $!\n";
	while(<HLT>){
		if(/^#/){
			($hap_NO) = (split)[1];
			$hap_NO --;
			next;
		}
		my ($string, $cn) = (split)[0,1];
		push @{$HapTypeStr_Aref->[$hap_NO]}, split //,$string for (1 .. $cn);
	}
	close HLT;

	# inform
	warn `date`."[INFO]:\tLoad haplotype structure ok.\n";
}

#--- read segment length ---
sub load_segment_length{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $SegmentLen_Href = $Option_Href->{SegmentLen_Href};

	open (SL, Try_GZ_Read($segment_length_file)) || die "fail read segment_length_file: $!\n";
	while(<SL>){
		my ($segid, $length) = (split)[0,1];
		$SegmentLen_Href->{$segid} = $length;
	}
	close SL;

	# inform
	warn `date`."[INFO]:\tLoad segment length ok.\n";
}

#--- find all connections of two breakpoints in local haplotype ---
sub find_all_connection{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $HapTypeStr_Aref = $Option_Href->{HapTypeStr_Aref};
	my $Distance_Href = $Option_Href->{Distance_Href};
	my $SegmentLen_Href = $Option_Href->{SegmentLen_Href};


	for my $str_Aref (@$HapTypeStr_Aref){
		my $SegAmount = scalar(@$str_Aref);
		my $Connect_Href = {};
		# find all neighbour connections
		for my $BK_Aref ( [$BK_1,$BK_2], [$BK_2,$BK_1] ){
			my $i_BK = $BK_Aref->[0];
			my $j_BK = $BK_Aref->[1];
			# find i_BK
			for my $i ( 0 .. $SegAmount-2 ){
				my $i_lf_seg = $str_Aref->[$i];
				my $i_rt_seg = $str_Aref->[$i+1];
				next if( "$i_lf_seg$i_rt_seg" ne $i_BK );
				# find nearest j_BK from upstream
				for my $j ( reverse 1 .. $i ){
					my $j_lf_seg = $str_Aref->[$j-1];
					my $j_rt_seg = $str_Aref->[$j];
					# last if( "$j_lf_seg$j_rt_seg" eq $i_BK ); # found itself
					next if( "$j_lf_seg$j_rt_seg" ne $j_BK ); # till found partner
					# record
					my ($lf_idx, $rt_idx) = ($j, $i);
					# get length
					my $conn_dist = 0;
					$conn_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( $lf_idx .. $rt_idx );
					if($conn_dist >= $max_allow_span_size){
						last;
					}
					else{
						$Connect_Href->{"$lf_idx,$rt_idx"} = 1;
					}
					# last; # stop when find the nearest j_BK
				}
				# find nearest j_BK from downstream
				for my $j ( $i+1 .. $SegAmount-2 ){
					my $j_lf_seg = $str_Aref->[$j];
					my $j_rt_seg = $str_Aref->[$j+1];
					# last if( "$j_lf_seg$j_rt_seg" eq $i_BK ); # found itself
					next if( "$j_lf_seg$j_rt_seg" ne $j_BK ); # till found partner
					# record
					my ($lf_idx, $rt_idx) = ($i+1, $j);
					# get length
					my $conn_dist = 0;
					$conn_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( $lf_idx .. $rt_idx );
					if($conn_dist >= $max_allow_span_size){
						last;
					}
					else{
						$Connect_Href->{"$lf_idx,$rt_idx"} = 1;
					}
					# last; # stop when find the nearest j_BK
				}
			}
		}
		# find the nearest edge (meet any one of BK_1 and BK_2)
		for my $conn_idx_pair (keys %$Connect_Href){
			my ($lf_idx, $rt_idx) = split /,/,$conn_idx_pair;
			# find nearest BK_1/2 from upstream
			my $upstream_dist = 0;
			for my $idx ( reverse 1 .. $lf_idx-1 ){
				my $lf_seg = $str_Aref->[$idx-1];
				my $rt_seg = $str_Aref->[$idx];
				next if( "$lf_seg$rt_seg" ne $BK_1 && "$lf_seg$rt_seg" ne $BK_2 );
				# record
				$upstream_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( $idx .. $lf_idx-1 );
				last; # stop when find the nearest BK
			}
			$upstream_dist ||= 3E10;
			# find nearest BK_1/2 from downstream
			my $downstream_dist = 0;
			for my $idx ( $rt_idx+1 .. $SegAmount-2 ){
				my $lf_seg = $str_Aref->[$idx];
				my $rt_seg = $str_Aref->[$idx+1];
				next if( "$lf_seg$rt_seg" ne $BK_1 && "$lf_seg$rt_seg" ne $BK_2 );
				# record
				$downstream_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( $rt_idx+1 .. $idx );
				last; # stop when find the nearest BK
			}
			$downstream_dist ||= 3E10;
			# get distance
			my $conn_dist = 0;
			$conn_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( $lf_idx .. $rt_idx );
			# compare
			my $min_nearest_dist = min($upstream_dist, $downstream_dist);
			if($min_nearest_dist == 3E10){
				$min_nearest_dist = 0;
				$min_nearest_dist += $SegmentLen_Href->{$str_Aref->[$_]} for ( 0 .. $SegAmount-1 );
				$min_nearest_dist -= $conn_dist;
			}
			my $ceil_dist = $conn_dist + $min_nearest_dist;
			# record
			$Distance_Href->{$conn_dist}->{$ceil_dist} ++;
		}
	}

	# inform
	warn `date`."[INFO]:\tFind all connections of two breakpoints ($BK_1,$BK_2) ok.\n";
}

#--- calculate each distance probability ---
sub output_distance_probability{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $BarcSpanStat_Href = $Option_Href->{BarcSpanStat_Href};
	my $Distance_Href = $Option_Href->{Distance_Href};

	my $merge_dist_prob = 0;
	for my $conn_dist (sort {$a<=>$b} keys %$Distance_Href){
		for my $ceil_dist (sort {$a<=>$b} keys %{$Distance_Href->{$conn_dist}}){
			# get the prob
			my $sum_count = 0;
			for my $span (sort {$a<=>$b} keys %$BarcSpanStat_Href){
				next if($span < $conn_dist);
				last if($span > $ceil_dist);
				$sum_count += $BarcSpanStat_Href->{$span};
			}
			my $prob = sprintf "%.5f", $sum_count/$all_count;
			last if( $prob == 0);
			my $CN = $Distance_Href->{$conn_dist}->{$ceil_dist};
			$merge_dist_prob += $CN * $prob;
			warn "$conn_dist\t$ceil_dist\t$CN\t$prob\n";
		}
	}

	print "$BK_1\t$BK_2\t$merge_dist_prob\n";
}
