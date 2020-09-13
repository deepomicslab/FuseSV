package FuseSV::Virus_Integ::EvaluateIns;

use strict;
use warnings;
use File::Spec qw/catfile/;
use List::Util qw[min max sum];
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr trible_run_for_success/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  evaluate_insert_size_for_all_ReadsGroup
			  get_average_value
			  get_trimmed_average
			  engineer_3times_SD_evaluation
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'EvaluateIns';
#----- version --------
$VERSION = "0.04";
$DATE = '2017-06-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						evaluate_insert_size_for_all_ReadsGroup
						load_insert_size_report
						SingleRG_engineer_3times_SD_evaluation
						get_average_value
						get_trimmed_average
						engineer_3times_SD_evaluation
						write_RG_info_report
						read_bam_for_proper_mapped_pe_reads
					 /;

#--- evaluate insert size from bam files OR read previous RG report ---
sub evaluate_insert_size_for_all_ReadsGroup{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	# basic
	my $bam_Href = $Option_Href->{bam_Href};
	my $source_bam_Href = $Option_Href->{source_bam_Href};
	my $out_prefix = $Option_Href->{out_prefix};
	my $Standard_Amount = $Option_Href->{Standard_Amount};
	my $min_pid_for_ins = $Option_Href->{min_pid_for_ins};
	# tools
	my $samtools = $Option_Href->{samtools};
	# actions bool
	my $SKIP_EVALUE_INS_bool = $Option_Href->{SKIP_EVALUE_INS};
	my $Allow_No_RG_ID_bool = $Option_Href->{Allow_No_RG_ID};
	my $Allow_SoftClip_For_InsEvalue_bool = $Option_Href->{Allow_SoftClip_For_InsEvalue};
	# global containers
	my $RG_Href = $Option_Href->{RG_Href};
	my $INS_Info_Href = $Option_Href->{INS_Info_Href};
	my $Rpref2Info_Href = $Option_Href->{Rpref2Info_Href};
	my $Bam2RG_ID_Href = $Option_Href->{Bam2RG_ID_Href};
	## about PCW
	my $PCW_Type_Href = $Option_Href->{PCW_Type_Href};
	my $Non_PCW_TYPE_Href = $Option_Href->{Non_PCW_TYPE_Href};
	## about the maximum value
	my $Max_Ins_Sref = $Option_Href->{Max_Ins_Sref};
	my $Max_SD_Sref = $Option_Href->{Max_SD_Sref};


	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	for my $tissue (sort keys %$bam_Href){

		# ins report folder
		my $ins_dir = File::Spec->catfile("$out_prefix-evaluate_ins", $tissue);

		if($SKIP_EVALUE_INS_bool){
			#-------------------------------#
			# load previous ins report file #
			#-------------------------------#
			# find ins report
			my $ins_report_glob = File::Spec->catfile($ins_dir, "RG_ID.*.ins.report");
			# read ins report
			for my $ins_report ( sort glob($ins_report_glob) ){
				&load_insert_size_report(
						{
							ins_report => $ins_report,
							tissue_type => $tissue,
							INS_Info_Href => $INS_Info_Href,
							Rpref2Info_Href => $Rpref2Info_Href,
							PCW_Type_Href => $PCW_Type_Href,
							Non_PCW_TYPE_Href => $Non_PCW_TYPE_Href
						}
					);
			}
			# alert
			stout_and_sterr `date`."Load All RG_IDs' information of tissue $tissue ok.\n";
		}
		else{
			# dir to store RG report
			`rm -rf $ins_dir`;
			`mkdir -p $ins_dir`;
			# to store ins count of each RG_ID
			my %RG_ID2INS_Stat;
			# overlap 3' end
			my %End3_Overlap;
			# RG_ID pid counts
			my %RG_ID2Stat_Amount;
			for my $bam (sort keys %{$bam_Href->{$tissue}}){
				$RG_ID2Stat_Amount{$_} = 0 for @{ $Bam2RG_ID_Href->{$tissue}->{$bam} };
			}
			# remained RG_ID
			my $Remain_RG_ID = scalar(keys %RG_ID2Stat_Amount);
			#---------------------------#
			# read bam file to load ins #
			#---------------------------#
			for my $useful_bam_Href ($bam_Href, $source_bam_Href){
				next unless( exists($useful_bam_Href->{$tissue}) );
				for my $bam (sort keys %{$useful_bam_Href->{$tissue}}){
					&read_bam_for_proper_mapped_pe_reads(
							{
								bam => $bam,
								tissue => $tissue,
								samtools => $samtools,
								Allow_No_RG_ID => ($Allow_No_RG_ID_bool || 0),
								Allow_SoftClip_For_InsEvalue => ($Allow_SoftClip_For_InsEvalue_bool || 0),
								RG_Href => $RG_Href,
								RG_ID2Stat_Amount_Href => \%RG_ID2Stat_Amount,
								Standard_Amount => $Standard_Amount,
								RG_ID2INS_Stat_Href => \%RG_ID2INS_Stat,
								Remain_RG_ID_Sref => \$Remain_RG_ID,
								End3_Overlap_Href => \%End3_Overlap
							}
						);
				}
				# all RG_IDs are ok for ins stat
				last if($Remain_RG_ID == 0);
			}
			#------------------------#
			#  Evaluate insert size  #
			#------------------------#
			my $RG_NO_count = 0;
			# each RG_ID
			for my $RG_ID (sort keys %RG_ID2INS_Stat){
				$RG_NO_count ++;
				## check the pid amount
				if($RG_ID2Stat_Amount{$RG_ID} < $min_pid_for_ins){
					stout_and_sterr `date`."RG $RG_ID used pid for ins stat is smaller than $min_pid_for_ins. So, ignore it.\n";
					next;
				}
				## get insertsize info
				my ( $ins_peak_Aref, $average_ins, $ins_sd ) = 
					&SingleRG_engineer_3times_SD_evaluation(
							{
								RG_ID2INS_Stat_Href => \%RG_ID2INS_Stat,
								RG_ID => $RG_ID,
								RG_ID2Stat_Amount_Href => \%RG_ID2Stat_Amount
							}
						);
				## overlap 3' end library ? 20% threshold.
				my $overlap3 = ((($End3_Overlap{$RG_ID}{1}||0) / $RG_ID2Stat_Amount{$RG_ID}) > 0.2)?'yes':'no';
				## record ins info for this RG_ID
				$INS_Info_Href->{$RG_ID} = {
										UsedPE => $RG_ID2Stat_Amount{$RG_ID},
										Tissue => $tissue,
										MaxRlen => $RG_Href->{$tissue}->{$RG_ID}->[2], # max_rlen_Aref
										PeakInsAmout => $ins_peak_Aref->[0],
										PeakIns => $ins_peak_Aref->[1],
										MeanIns => $average_ins,
										Ins_SD => $ins_sd,
										RidPref => $RG_Href->{$tissue}->{$RG_ID}->[0], # Readid_Same_Prefix
										LB_ID => $RG_Href->{$tissue}->{$RG_ID}->[3],
										Olp3end => $overlap3
									};
				## output ins report
				my $ins_report = File::Spec->catfile($ins_dir, "RG_ID.NO$RG_NO_count.ins.report");
				&write_RG_info_report(
						{
							ins_report => $ins_report,
							RG_ID => $RG_ID,
							RG_Info_Href => $INS_Info_Href
						}
					);
				## inform
				stout_and_sterr `date`."Evaluate ins of tissue $tissue RG_ID $RG_ID ok.\n";
			}
		}

		#--------------------------#
		#  gather ins information  #
		#--------------------------#
		for my $RG_ID (sort keys %$INS_Info_Href){

			## record the info to concerned read prefix
			$Rpref2Info_Href->{ $INS_Info_Href->{$RG_ID}->{RidPref} } = 
					[
						$INS_Info_Href->{$RG_ID}->{MaxRlen},
						$INS_Info_Href->{$RG_ID}->{MeanIns},
						int( $INS_Info_Href->{$RG_ID}->{Ins_SD} ),
						$INS_Info_Href->{$RG_ID}->{Olp3end},
						$INS_Info_Href->{$RG_ID}->{Tissue},
						$INS_Info_Href->{$RG_ID}->{LB_ID}
					];


			## PCW type for next bam reading both unmapped reads
			if( $INS_Info_Href->{$RG_ID}->{Olp3end} eq 'yes' ){
				$PCW_Type_Href->{$RG_ID} = 1;
			}
			else{
				$Non_PCW_TYPE_Href->{$RG_ID} = 1;
			}

			## store the max ins for following work
			($$Max_Ins_Sref,$$Max_SD_Sref) = ($$Max_Ins_Sref < $INS_Info_Href->{$RG_ID}->{MeanIns})
												  ? ($INS_Info_Href->{$RG_ID}->{MeanIns},$INS_Info_Href->{$RG_ID}->{Ins_SD})
												  : ($$Max_Ins_Sref,$$Max_SD_Sref);
		}

		#--------------------------------#
		# check all RG_ID of this tissue #
		#--------------------------------#
		for my $RG_ID (keys %{ $RG_Href->{$tissue} }){
			if(!exists($INS_Info_Href->{$RG_ID})){ # sometimes, rarely, some RG_ID cannot get enough information for INS_Info estimate.
				stout_and_sterr `date`."<WARN>\tRG_ID: $RG_ID of Tissue $tissue lacks enough information of estimation. All work will skip it.\n";
			}
		}
	}
}

#--- read bam for insertsize of properly mapped PE reads ---
sub read_bam_for_proper_mapped_pe_reads{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $bam = $Option_Href->{bam};
	my $samtools = $Option_Href->{samtools};
	my $Allow_No_RG_ID_bool = $Option_Href->{Allow_No_RG_ID};
	my $tissue = $Option_Href->{tissue};
	my $RG_Href = $Option_Href->{RG_Href};
	my $RG_ID2Stat_Amount_Href = $Option_Href->{RG_ID2Stat_Amount_Href};
	my $Standard_Amount = $Option_Href->{Standard_Amount};
	my $Allow_SoftClip_For_InsEvalue_bool = $Option_Href->{Allow_SoftClip_For_InsEvalue};
	my $RG_ID2INS_Stat_Href = $Option_Href->{RG_ID2INS_Stat_Href};
	my $Remain_RG_ID_Sref = $Option_Href->{Remain_RG_ID_Sref};
	my $End3_Overlap_Href = $Option_Href->{End3_Overlap_Href};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	# to store the pid
	my %RG_ID2Stat_Pid;
	# properly mapped(0x2) and non-dup(0x400), OLD: [awk '\$2~/P/ && \$2!~/d/' |]
	# -F 0x800, supplementary alignment
	# 0x100 + 0x800 + 0x400 = 0xD00
	open (BAM,"$samtools view -f 0x2 -F 0xD00 $bam |") || die"fail $bam: $!\n";
	while(<BAM>){
		# RG_ID
		my ($RG_ID) = (/RG\:Z\:(\S+)/);
		if( !$RG_ID ){
			warn_and_exit "Cannot get RG_ID from bam reads line $. of $tissue bam:\n$bam\n";
		}
		die "The RG_ID $RG_ID is not recognized\nGot from $tissue bam:\n$bam\n" if(!exists($RG_Href->{$tissue}->{$RG_ID}));
		## RG_ID info
		my ($rid_pref,$bam,$max_rlen_Aref) = @{$RG_Href->{$tissue}->{$RG_ID}};
		## if reach the Standard_Amount, no need to stat for this RG_ID
		next if($RG_ID2Stat_Amount_Href->{$RG_ID} == $Standard_Amount);
		# info of mapping
		my ($pid,$FLAG,$CIGAR,$TLEN,$rseq) = (split)[0,1,5,8,9];
		## add the customized readid/pid prefix
		$pid = $rid_pref.$pid;
		if(exists($RG_ID2Stat_Pid{$RG_ID}{$pid})){
			die "This Bam's properly mapped paired-end, the 'R' and 'r' is not right!\n$bam\n$_\n" unless($FLAG & 0x10);
			# if($CIGAR !~ /S|N|H/){ # good read
			if(     $CIGAR !~ /N|H/
			    &&  ( $CIGAR !~ /S/ || $Allow_SoftClip_For_InsEvalue_bool )
			   ){ # good read
			   	# this is the minus-orientation mapped end, and must not been large gap, such as soft-clipped
				# however, once allow the soft-clipped reads into, calculate the bilateral 'S' part
				# plus -> left; minus -> right
				if( $Allow_SoftClip_For_InsEvalue_bool && $CIGAR =~ /(\d+)S$/ ){
					# retrieve the right S part
					$RG_ID2Stat_Pid{$RG_ID}{$pid}->[0] += $1;
				}
				my ($ins,$up_rlen) = @{$RG_ID2Stat_Pid{$RG_ID}{$pid}};
				$RG_ID2INS_Stat_Href->{$RG_ID}->{$ins}++; # record
				$RG_ID2Stat_Amount_Href->{$RG_ID}++; # add 1 to the ok stat pid amount
				if($RG_ID2Stat_Amount_Href->{$RG_ID} == $Standard_Amount){ # if this RG_ID is ok for ins stat
					$$Remain_RG_ID_Sref --;
					delete $RG_ID2Stat_Pid{$RG_ID}; # release memory
				}
				### overlap at the 3' end ?
				my $overlap_sign = ($ins > ($up_rlen+length($rseq)))?0:1;
				$End3_Overlap_Href->{$RG_ID}->{$overlap_sign}++; # record number
				### the read length
				if($FLAG & 0x40){
					$$max_rlen_Aref[0] = max($$max_rlen_Aref[0] , length($rseq));
					$$max_rlen_Aref[1] = max($$max_rlen_Aref[1] , $up_rlen);
				}
				else{ # now-read is /2
					$$max_rlen_Aref[0] = max($$max_rlen_Aref[0] , $up_rlen);
					$$max_rlen_Aref[1] = max($$max_rlen_Aref[1] , length($rseq));
				}
			}
			## never use it, bye-bye
			delete $RG_ID2Stat_Pid{$RG_ID}{$pid};
		}
		elsif(     ( $FLAG & 0x20  && $CIGAR !~ /N|H/ )
			   &&  ( $CIGAR !~ /S/ || $Allow_SoftClip_For_InsEvalue_bool )
			  ){ # store the TLEN, it is the insertsize.
			# this is the plus-orientation mapped end, and must not been large gap, such as soft-clipped
			# however, once allow the soft-clipped reads into, calculate the bilateral 'S' part
			# plus -> left; minus -> right
			if( $Allow_SoftClip_For_InsEvalue_bool && $CIGAR =~ /^(\d+)S/ ){
				# retrieve the left S part
				$TLEN += $1;
			}
			$RG_ID2Stat_Pid{$RG_ID}{$pid} = [ $TLEN, length($rseq) ];
		}
		# all RG_IDs are ok for ins stat
		last if($$Remain_RG_ID_Sref == 0);
	}
	close BAM;
	# inform
	stout_and_sterr `date`."[INFO]\tExtract All RG_IDs' information of bam file ok:\n"
						  ."\t$bam\n";
}

#--- write ReadsGroup info report ---
sub write_RG_info_report{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $ins_report = $Option_Href->{ins_report};
	my $RG_Info_Href  = $Option_Href->{RG_Info_Href};
	my $RG_ID = $Option_Href->{RG_ID};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	# structure of report
	my $single_RG_ID_Info_Href = $RG_Info_Href->{$RG_ID};
	my @Structure = (
						['RG_ID',         $RG_ID],
						['#_of_used_PE',  $single_RG_ID_Info_Href->{UsedPE}],
						['Max-Rlen-End1', $single_RG_ID_Info_Href->{MaxRlen}->[0]],
						['Max-Rlen-End2', $single_RG_ID_Info_Href->{MaxRlen}->[1]],
						['Peak-ins',      $single_RG_ID_Info_Href->{PeakIns}],
						['Peak-Amount',   $single_RG_ID_Info_Href->{PeakInsAmout}],
						['Mean-ins',      $single_RG_ID_Info_Href->{MeanIns}],
						['Stat-SD-ins',   $single_RG_ID_Info_Href->{Ins_SD}],
						['Rid-samepref',  $single_RG_ID_Info_Href->{RidPref}],
						['Library-ID',    $single_RG_ID_Info_Href->{LB_ID}],
						['Overlap3end',   $single_RG_ID_Info_Href->{Olp3end}]
					);
	# output
	open (IR,Try_GZ_Write($ins_report)) || die"fail $ins_report: $!\n";
	printf IR "%-15s%s\n", ($_->[0].':'), $_->[1] for @Structure;
	close IR;
}

#--- engineering 3-times-SD filtering for one RG ID ---
sub SingleRG_engineer_3times_SD_evaluation{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $RG_ID2INS_Stat_Href = $Option_Href->{RG_ID2INS_Stat_Href};
	my $RG_ID = $Option_Href->{RG_ID};
	my $RG_ID2Stat_Amount_Href = $Option_Href->{RG_ID2Stat_Amount_Href};

	# system variants
	# my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	my $Dist_Href = &engineer_3times_SD_evaluation(
										{
											InsStat_Href => $RG_ID2INS_Stat_Href->{$RG_ID},
											stat_amount  => $RG_ID2Stat_Amount_Href->{$RG_ID}
										}
									  );

	# update
	$RG_ID2Stat_Amount_Href->{$RG_ID} = $Dist_Href->{stat_amount};

	# return
	return  (	$Dist_Href->{peak_value_Aref},
				$Dist_Href->{average_value},
				$Dist_Href->{value_sd}          );
}

#--- return average(mean) value ---
# default to use TMM (ratio=0), means to accept all records.
sub get_average_value{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $value_Aref = $Option_Href->{value_Aref};
	my $Eng_3SD_bool = $Option_Href->{Eng_3SD_bool} || 0;
	my $TrimRatio = $Option_Href->{TrimRatio} || 0;

	if( $Eng_3SD_bool ){
		my %Stat_value2count;
		$Stat_value2count{$_}++ for @$value_Aref;
		my $Dist_Href = &engineer_3times_SD_evaluation(
							{
								InsStat_Href => \%Stat_value2count,
								stat_amount  => scalar(@$value_Aref)
							}
						);
		return int($Dist_Href->{average_value});
	}
	else{
		return int( &get_trimmed_average(
						{
							array_Aref => $value_Aref,
							TrimRatio => $TrimRatio
						}
					)
			   );
	}
}

#--- TMM average method ---
sub get_trimmed_average{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $array_Aref = $Option_Href->{array_Aref};
	my $TrimRatio = $Option_Href->{TrimRatio} || 0;

	if( $TrimRatio == 0 ){
		return sum(@$array_Aref) / scalar(@$array_Aref);
	}
	else{
		my $count = scalar(@$array_Aref);
		my $thereshold = int(abs($TrimRatio) * $count);
		my @numbers = sort {$a<=>$b} @$array_Aref;
		# only tail operation
		if( $TrimRatio < 0 ){
			return sum(@numbers[$thereshold .. $#numbers]) / ($count-$thereshold);
		}
		# bilaterally operation
		else{
			return sum(@numbers[$thereshold .. $#numbers-$thereshold]) / ($count-2*$thereshold);
		}
	}
}

#--- engineering 3-times-SD filtering ---
sub engineer_3times_SD_evaluation{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $InsStat_Href = $Option_Href->{InsStat_Href};
	my $stat_amount = $Option_Href->{stat_amount};

	## get value distribution
	my @value_peak = (0,0);
	my $average_value = 0;
	my $value_sd = 0;
	my $median_value = 0;
	ENG_GET_INS_INFO:{ # engineering 3-times-SD filtering
		## peak and sum of value
		my $sum = 0;
		for my $value (sort {$a<=>$b} keys %$InsStat_Href){
			my $amount = $InsStat_Href->{$value};
			## the peak of value distribution
			@value_peak = ($amount,$value) if($amount > $value_peak[0]);
			## sum
			$sum += $value * $amount;
		}
		## average value
		$average_value = $sum / $stat_amount;
		## sd of value
		for my $value (sort {$a<=>$b} keys %$InsStat_Href){
			my $amount = $InsStat_Href->{$value};
			$value_sd += (($value - $average_value) ** 2) * $amount;
		}
		$value_sd = sqrt($value_sd / $stat_amount);
		## 3-times-SD filtering
		my $obvious_abnormal_number = 0;
		for my $OB_abnormal_value (grep abs($_-$average_value) > 3*$value_sd , keys %$InsStat_Href){
			$obvious_abnormal_number += $InsStat_Href->{$OB_abnormal_value};
			delete $InsStat_Href->{$OB_abnormal_value}; # discard the abnormals
		}
		if($obvious_abnormal_number){
			$stat_amount -= $obvious_abnormal_number; # change the number of pid-for-value-evalue
			redo ENG_GET_INS_INFO;
		}
		# median value
		my $count = 0;
		for my $value (sort {$a<=>$b} keys %$InsStat_Href){
			$count += $InsStat_Href->{$value};
			if($count >= $stat_amount/2){
				$median_value = $value;
				last;
			}
		}
		# refine average_value
		$average_value = int($average_value * 100) / 100;
		$value_sd = int($value_sd * 100) / 100;
	}

	# return
	return  {	peak_value_Aref => \@value_peak,
				median_value => $median_value,
				average_value => $average_value,
				value_sd => $value_sd,
				stat_amount => $stat_amount        };
}

#--- read insert size info from one report ---
sub load_insert_size_report{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $ins_report  = $Option_Href->{ins_report};
	my $tissue = $Option_Href->{tissue_type};
	my $INS_Info_Href = $Option_Href->{INS_Info_Href};
	my $Rpref2Info_Href = $Option_Href->{Rpref2Info_Href};
	# about the PCW type
	my $PCW_Type_Href = $Option_Href->{PCW_Type_Href};
	my $Non_PCW_TYPE_Href = $Option_Href->{Non_PCW_TYPE_Href};

	# system variants
	# my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	my ($RG_ID,$LB_ID);
	my ($average_ins,$ins_sd);
	my $Readid_Same_Prefix;
	my $overlap3;
	my ($MaxRlen1,$MaxRlen2);
	open (IR,Try_GZ_Read($ins_report))||die"fail $ins_report: $!\n";
	while(<IR>){
		my ($item, $value) = (split /[\s\:]+/)[0,1];
		if($item eq 'RG_ID'){
			$RG_ID = $value;
			$INS_Info_Href->{$RG_ID} = {};
			$INS_Info_Href->{$RG_ID}->{Tissue}  = $tissue;
		}
		# elsif($item eq 'Max-Rlen-User'){
		# 	$INS_Info_Href->{$RG_ID}{MaxRlen} = $value;
		# }
		elsif($item eq '#_of_used_PE'){
			$INS_Info_Href->{$RG_ID}->{UsedPE} = $value;
		}
		elsif($item eq 'Max-Rlen-End1'){
			$MaxRlen1 = $value;
		}
		elsif($item eq 'Max-Rlen-End2'){
			$MaxRlen2 = $value;
		}
		elsif($item eq 'Peak-ins'){
			$INS_Info_Href->{$RG_ID}->{PeakIns} = $value;
		}
		elsif($item eq 'Peak-Amount'){
			$INS_Info_Href->{$RG_ID}->{PeakInsAmout} = $value;
		}
		elsif($item eq 'Mean-ins'){
			$INS_Info_Href->{$RG_ID}->{MeanIns} = $value;
			$average_ins = $value;
		}
		elsif($item eq 'Stat-SD-ins'){
			$INS_Info_Href->{$RG_ID}->{Ins_SD} = $value;
			$ins_sd = $value;
		}
		elsif($item eq 'Rid-samepref'){
			$INS_Info_Href->{$RG_ID}->{RidPref} = $value;
			$Readid_Same_Prefix = $value;
		}
		elsif($item eq 'Library-ID'){
			$INS_Info_Href->{$RG_ID}->{LB_ID} = $value;
			$LB_ID = $value;
		}
		elsif($item eq 'Overlap3end'){
			$INS_Info_Href->{$RG_ID}->{Olp3end} = $value; # yes or no
			$overlap3 = $value;
		}
	}
	# max rlen
	$INS_Info_Href->{$RG_ID}->{MaxRlen} = [$MaxRlen1,$MaxRlen2];

	# inform
	stout_and_sterr `date`."Load Insert size info of tissue $tissue RG_ID $RG_ID ok.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
