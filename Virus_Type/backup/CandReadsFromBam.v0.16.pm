package FuseSV::Virus_Type::CandReadsFromBam;

use strict;
use warnings;
use File::Basename qw/basename dirname/;
use List::Util qw[min max sum];
use Getopt::Long;
use SOAPfuse::General_Operation qw/warn_and_exit trible_run_for_success stout_and_sterr/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Virus_Type::VirusTypeMain;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  extract_candidate_reads_from_bam
			  extract_fastq_from_nSortBam
			  read_host_fai
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'CandReadsFromBam';
#----- version --------
$VERSION = "0.16";
$DATE = '2018-07-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						extract_candidate_reads_from_bam
						load_RG_info_report
						read_host_fai
						get_RG_ID_prefix
						get_RG_ID_from_bam_line
						extract_fastq_from_nSortBam
						convert_nonhost_alignment_to_unmap
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} cand_reads <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -bam    [s]  the indexed-bam alignment file. <required>
	                  Assume to be in Sanger-quality format (column NO.11).
	    -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
	                  only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
	    -odir   [s]  directory to store outputs, must exist already. <required>
	    -host   [s]  the .fai file of only host genome reference. <required>
	                  such as pure human reference fai file, created by samtools faidx.

	   # Software Required #
	    -stl    [s]  SamTools, minimum version: 1.3. <required>

	   # Candidate Reads Option #
	    -tstag  [s]  add ONE uppercase prefix to read-id to distinguish reads from different tissues. [Disabled]
	    -minmq  [i]  PE reads of which anyone's mapQ lower than INT will be extract for virus analysis. [30]
	                  The properly mapping PE-reads (i.e., 0x2 in flag) will be all ignored.
	    -minrst [i]  only prepare tmp file for host-refseg larger than INT when extracting candidate reads. [3E8]
	                  Donot set it very small, because the more tmp-file handles, more buffer memory will be cost.
	    -minlen [i]  minimum candidate read length. [50]
	    -disdup      discard duplicated reads with SAM-flag 0x400. [Disabled]
	    -bartag [s]  TAG in bam to extract barcode info, such as 'BX:Z:' of 10X technology. [none]
	                  The extracted barcode will be concatenated with paired-end ID.

	    -h|help   Display this help info.

 	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
	FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool;
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-bam:s"	=> \$V_Href->{bam},
		"-id:s"		=> \$V_Href->{sample_ID},
		"-odir:s"	=> \$V_Href->{virusAnaFolder},
		"-host:s"	=> \$V_Href->{host_fai},
		# software
		"-stl:s"	=> \$V_Href->{samtools},
		# options
		## cand_reads for alignment on virus
		"-minmq:i"	=> \$V_Href->{min_mapQ},
		"-minrst:i"	=> \$V_Href->{min_refseg_size_for_tmp},
		"-minlen:i"	=> \$V_Href->{min_CandRead_length},
		"-disdup"	=> \$V_Href->{discard_DupReads},
		"-tstag:s"	=> \$V_Href->{tissue_tag},
		"-bartag:s"	=> \$V_Href->{barcode_SamTag},
		# "-tltag:s"	=> \$V_Href->{Tool_Tag}, # hidden option
		# help
		"-h|help"	=> \$V_Href->{HELP},
		# for debug
		"-debug"	=> \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || (!$V_Href->{bam} || !-e $V_Href->{bam})
			 || (!$V_Href->{sample_ID} || $V_Href->{sample_ID} =~ /[^\w\-]/)
			 || (!$V_Href->{virusAnaFolder} || !-d $V_Href->{virusAnaFolder})
			 || (!$V_Href->{host_fai} || !-e $V_Href->{host_fai})
			 || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
			 || (defined $V_Href->{tissue_tag} && $V_Href->{tissue_tag} !~ /^[A-Z]$/)
			);
}

#--- read bam file to collect the candidate reads in PE mode ---
# 1. non-host ref; 2. low_mapQ<$min_mapQ (host); 3. unmapped reads
sub extract_candidate_reads_from_bam{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Host_Refseg_Href = $V_Href->{Host_Refseg_Href};
	my $RG_Info_Href = $V_Href->{RG_Info_Href};

	$V_Href->{CandFqFolder} = GetPath( filekey => 'CandFqFolder' );
	$V_Href->{RG_info_report} = GetPath( filekey => 'RG_info_report' );

	#---------------------------#
	# mode 1. extract candidate #
	#  reads from given bam     #
	#---------------------------#

	# run mode check
	unless(    $V_Href->{run_mode} == 0 # run_whole mode must run this step
			|| $V_Href->{run_mode} == 1 # just this step
		  ){
		if(!-e $V_Href->{RG_info_report}){ # In other modes, this report will be loaded on at first
			warn_and_exit "<ERROR>\tCannot find report for each ReadGroup (RG):\n".
						  "\t $V_Href->{RG_info_report}\n".
						  "\tplease run command 'cand_reads'.\n";
		}
		else{
			&load_RG_info_report;
		}
		return;
	}

	# load host-refseg(s)
	&read_host_fai;

	# for memory saving, prepare to store anomalous pairs with other refseg
	$V_Href->{CandFq_TmpFolder} = GetPath( filekey => 'CandFq_TmpFolder' );
	`rm -rf $V_Href->{CandFq_TmpFolder}`;
	`mkdir -p $V_Href->{CandFq_TmpFolder}`;

	#---------------------------#
	# read bam header to set up #
	#---------------------------#
	my $refseg_NO = 0;
	my @bam_header;
	my $RG_ID_NO = 0;
	my @HostRefseg_to_tmp;
	open (BAM_H,"$V_Href->{samtools} view -H $V_Href->{bam} |") || die"fail samtools view header of $V_Href->{bam}: $!\n";
	while(<BAM_H>){
		if(/^\@SQ/){
			my ($refseg,$ref_len) = (/SN:(\S+)\sLN:(\d+)/);
			# only need host refseg
			next if(!exists($Host_Refseg_Href->{$refseg}));
			push @bam_header, $_; # store for next output
			# prepare tmp for this ref
			next if($ref_len < $V_Href->{min_refseg_size_for_tmp});
			push @HostRefseg_to_tmp, $refseg;
			$refseg_NO++;
			$Host_Refseg_Href->{$refseg}->{tmp_file} = GetPath( filekey => 'CandFq_RefSegTmpFile', refseg_NO => $refseg_NO );
			open ($Host_Refseg_Href->{$refseg}->{tmp_file_wH}, Try_GZ_Write($Host_Refseg_Href->{$refseg}->{tmp_file})) || die "fail write $Host_Refseg_Href->{$refseg}->{tmp_file}: $!\n";
			stout_and_sterr `date`."[INFO]\tOpen tmp file handle for host refseg: $refseg.\n";
		}
		elsif(/^\@RG/){
			push @bam_header, $_; # store for next output
			my ($RG_ID) = (/ID:(\S+)/);
			if(exists($RG_Info_Href->{$RG_ID})){ # avoid the duplicated RG(s)
				warn_and_exit "<ERROR>\tThe RG $RG_ID shows two times in your bam:\n\t $V_Href->{bam}\n";
			}
			else{
				$RG_ID_NO++;
				$RG_Info_Href->{$RG_ID}->{RG_NO} = $RG_ID_NO;
				$RG_Info_Href->{$RG_ID}->{RG_rID_pref} = $V_Href->{tissue_tag} . get_RG_ID_prefix( $RG_ID_NO );
				$RG_Info_Href->{$RG_ID}->{CandFqGz} = GetPath( filekey => 'CandFqGz', CandFq_RG_ID_NO => $RG_ID_NO );
				$RG_Info_Href->{$RG_ID}->{Cand_nSortBam} = GetPath( filekey => 'Cand_nSortBam', CandFq_RG_ID_NO => $RG_ID_NO );
				$RG_Info_Href->{$RG_ID}->{max_rlen} = [0,0];
			}
		}
		elsif(/^\@HD/){ # discard @PG and @CO
			push @bam_header, $_; # store for next output
		}
	}
	close BAM_H;

	# informs on the host-refseg tmp files
	if($refseg_NO){
		stout_and_sterr `date`."[INFO]\tOpen tmp file handles for $refseg_NO host refsegs ok.\n";
	}
	else{
		`rm -rf $V_Href->{CandFq_TmpFolder}`; # no use
	}

	#-------------------------------#
	# open file-handle(s) for RG(s) #
	#-------------------------------#
	for my $RG_ID (sort keys %$RG_Info_Href){
		my $RG_ID_NO = $RG_Info_Href->{$RG_ID}->{RG_NO};
		$RG_Info_Href->{$RG_ID}->{CandBam} = GetPath( filekey => 'CandBam', CandFq_RG_ID_NO => $RG_ID_NO );
		open ($RG_Info_Href->{$RG_ID}->{CandBam_wH},"| $V_Href->{samtools} view -b -S -o $RG_Info_Href->{$RG_ID}->{CandBam}") || die"fail samtools write $RG_Info_Href->{$RG_ID}->{CandBam}: $!\n";
		for my $header_line (@bam_header){
			if(    $header_line !~ /^\@RG/
				|| $header_line =~ /ID:$RG_ID\s/
			  ){
				print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $header_line;
			}
		}
	}

	# inform
	stout_and_sterr `date`."[INFO]\tstarts extracting candidate reads from bam.\n";

	#------------------------------#
	# load bam for candidate reads #
	#------------------------------#
	my %ProperMap;
	my %pid_mapQ;
	my $last_Mrefseg = '';
	#   0x800: to avoid supplementary alignment
	# + 0x100: to avoid secondary alignment
	# = 0x900
	open (BAM,"$V_Href->{samtools} view -F 0x900 $V_Href->{bam} |") || die"fail samtools view $V_Href->{bam}: $!\n";
	while(<BAM>){
		my $bam_line = $_;
		my @bam_line = (split /\t+/);
		my ($pid, $flag, $self_Mrefseg, $MapQ, $Cigar, $mate_Mrefseg, $r_seq) = @bam_line[0,1,2,4,5,6,9];
		# read length
		my $rlen = length($r_seq);

		#-----------------------------------#
		# pairs properly map to host-refseg #
		#-----------------------------------#
		# once has SoftClip, count into cand_reads
		if( exists($Host_Refseg_Href->{$self_Mrefseg}) && ($flag & 0x2) ){
			if( exists($ProperMap{$pid}) ){
				# have SoftClip
				## as the bam is coordinates sorted.
				## the stored end must be puls-mapped, it must has SoftClip at its left edge if it's one junc-read
				## the new end must be minus-mapped, it must has SoftClip at its right edge if it's one junc-read
				if( $ProperMap{$pid}->{cigar} =~ /^\d+S/ || $Cigar =~ /S$/ ){
					# RG info
		  			my $RG_ID = &get_RG_ID_from_bam_line($bam_line);
		  			# output to specific file handle
		  			print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $ProperMap{$pid}->{bamline};
					print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $bam_line;
					# update max_rlen
					for my $Aref ([$ProperMap{$pid}->{flag}, $ProperMap{$pid}->{rlen}], [$flag, $rlen]){
						my $end_NO = ($Aref->[0] & 0x40) ? 0 : 1;
						$RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO] = max( $RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO], $Aref->[1] );
					}
				}
				# else{ # future for InsertSize evaluation
				# }
				# sweep
				delete $ProperMap{$pid};
			}
			else{ # store it for quick retrieve later
				$ProperMap{$pid} = { bamline=>$bam_line, flag=>$flag, cigar=>$Cigar, rlen=>$rlen };
			}
			# skip next operations
			next;
		}

		#------------------------------------#
		# tmp file specific for stroing mapQ #
		#------------------------------------#
		# new Mrefseg! load its tmp file.
		if($last_Mrefseg ne $self_Mrefseg){
			if(exists($Host_Refseg_Href->{$self_Mrefseg})){
				# inform
				stout_and_sterr `date`."[INFO]\tdealing host-refseg $self_Mrefseg\n";
				# tmp file?
				if(exists($Host_Refseg_Href->{$self_Mrefseg}->{tmp_file_wH})){ # some host-refseg, may not have tmp file
					close $Host_Refseg_Href->{$self_Mrefseg}->{tmp_file_wH};
					open (TMP,Try_GZ_Read($Host_Refseg_Href->{$self_Mrefseg}->{tmp_file})) || die"fail read $Host_Refseg_Href->{$self_Mrefseg}->{tmp_file}: $!\n";
					while(<TMP>){
						my @ele = split /\t+/;
						chomp(my $tmp_pid = $ele[0]);
						if(scalar(@ele) == 1){
							$pid_mapQ{$tmp_pid} = ""; # empty as it has been outputed previously
						}
						else{
							$pid_mapQ{$tmp_pid} = join("\t",@ele); # SAM-info, not chomp
						}
					}
					close TMP;
					# sweep
					`rm -rf $Host_Refseg_Href->{$self_Mrefseg}->{tmp_file}`;
					delete $Host_Refseg_Href->{$self_Mrefseg}->{tmp_file_wH};
					delete $Host_Refseg_Href->{$self_Mrefseg}->{tmp_file};
					stout_and_sterr `date`."[INFO]\tstop writing, and reload tmp file of host refseg: $self_Mrefseg, and then delete it ok.\n";
				}
			}
			elsif($self_Mrefseg eq '*'){ # unmapped PE-reads
				stout_and_sterr `date`."[INFO]\tdealing paired unmapped reads\n";
			}
			# update last_Mrefseg
			$last_Mrefseg = $self_Mrefseg;
		}

		#--------------------------------------#
		# comprehensive judgment of cand_reads #
		#--------------------------------------#
		# prepare mapQ of mated end
		my $Mate_MapQ = &get_MateMapQ_from_bam_line($bam_line);
		# multiple judgments
		if(    # case 2, anyone end is lower than min_mapQ. it requires bam must have the MQ tag.
			   (  ($Mate_MapQ ne 'na') && ( ($Mate_MapQ < $V_Href->{min_mapQ}) || ($MapQ < $V_Href->{min_mapQ}) ) )
			   # case 3, both ends are unmap
			|| (  ($flag & 0x4) &&  ($flag & 0x8) && ($mate_Mrefseg eq '*') )
			   # case 3, this end is unmap, and the mate one is not unmap, so mate Mrefseg must be '='
			|| (  ($flag & 0x4) && !($flag & 0x8) && ($mate_Mrefseg eq '=') )
			   # case 3, this end is not unmap, and the mate one is unmap
			|| ( !($flag & 0x4) &&  ($flag & 0x8) )
			   # case 1, this end maps to non-host refseg
			|| ( !exists($Host_Refseg_Href->{$self_Mrefseg}) )
			   # case 1, the mate end maps to non-host refseg
			|| ( ($mate_Mrefseg ne '=') && ($mate_Mrefseg ne '*') && !exists($Host_Refseg_Href->{$mate_Mrefseg}) )
		  ){
		  	# RG info
		  	my $RG_ID = &get_RG_ID_from_bam_line($bam_line);
		  	# convert the non-host alignment to unmaps
		  	if(    # case 1, this end maps to non-host refseg
				   ( !exists($Host_Refseg_Href->{$self_Mrefseg}) )
				   # case 1, the mate end maps to non-host refseg
				|| ( ($mate_Mrefseg ne '=') && ($mate_Mrefseg ne '*') && !exists($Host_Refseg_Href->{$mate_Mrefseg}) )
		  	  ){
			  	# convert
			  	&convert_nonhost_alignment_to_unmap(
			  			{
			  				self_Mrefseg => $self_Mrefseg,
			  				mate_Mrefseg => $mate_Mrefseg,
			  				bamline_Sref => \$bam_line,
			  				bamline_Aref => \@bam_line
				  		}
				  	);
		  	}
		  	# output to specific file handle
			print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $bam_line;
			# update max_rlen
			my $end_NO = ($flag & 0x40) ? 0 : 1;
			$RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO] = max( $RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO], $rlen );
		}
		elsif( $Mate_MapQ eq 'na' ){ # bam lacks MQ tag
			if( $MapQ < $V_Href->{min_mapQ} ){ # case 2, mapQ lower than min_mapQ

				# RG info
				my $RG_ID = &get_RG_ID_from_bam_line($bam_line);
			  	# output to specific file handle
				print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $bam_line;
				# update max_rlen
				my $end_NO = ($flag & 0x40) ? 0 : 1;
				$RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO] = max( $RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO], $rlen );

				# store to tell the mate_end to output
				if(!exists($pid_mapQ{$pid})){
					if($mate_Mrefseg eq '='){ # same refseg, so store in memory for quick checking
						$pid_mapQ{$pid} = '';
					}
					elsif(exists($Host_Refseg_Href->{$mate_Mrefseg}) && exists($Host_Refseg_Href->{$mate_Mrefseg}->{tmp_file_wH})){ # this is should be store in tmp file of mate_Mrefseg for later checking
						print { $Host_Refseg_Href->{$mate_Mrefseg}->{tmp_file_wH} } "$pid\n"; # only record the pid
					}
					else{ # the mate_Mrefseg is too small to consider store by tmp file
						$pid_mapQ{$pid} = '';
					}
				}
				else{
					if(length($pid_mapQ{$pid}) != 0){ # not empty, means hasn't been output (high_mapQ).
						# RG info
						my $RG_ID = &get_RG_ID_from_bam_line($pid_mapQ{$pid});
					  	# output to specific file handle
						print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $pid_mapQ{$pid}; # SAM info of mate_end
						# update max_rlen
						my ($mate_rseq) = (split /\t+/,$pid_mapQ{$pid})[9]; # read seq is at column NO.10 in SAM format
						my $mate_rlen = length($mate_rseq);
						my $end_NO = ($flag & 0x40) ? 0 : 1;
						$RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO] = max( $RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO], $mate_rlen );
					}
					delete $pid_mapQ{$pid};
				}
			}
			else{ # mapQ >= min_mapQ: we cannot know the mapQ of the other end now, so need to store.
				if(!exists($pid_mapQ{$pid})){ # must record for next end
					if($mate_Mrefseg eq '='){ # same refseg, so store in memory for quick checking
						$pid_mapQ{$pid} = $bam_line;
					}
					elsif(exists($Host_Refseg_Href->{$mate_Mrefseg}) && exists($Host_Refseg_Href->{$mate_Mrefseg}->{tmp_file_wH})){ # this is should be store in tmp file of mate_Mrefseg for later checking
						print { $Host_Refseg_Href->{$mate_Mrefseg}->{tmp_file_wH} } $bam_line; # SAM info, not chomp
					}
					else{ # the mate_Mrefseg is too small to consider store by tmp file
						$pid_mapQ{$pid} = $bam_line;
					}
				}
				else{
					if(length($pid_mapQ{$pid}) == 0){ # empty, means has been output (low_mapQ already).
						# RG info
						my $RG_ID = &get_RG_ID_from_bam_line($bam_line);
					  	# output to specific file handle
						print { $RG_Info_Href->{$RG_ID}->{CandBam_wH} } $RG_Info_Href->{$RG_ID}->{RG_rID_pref} . $bam_line;
						# update max_rlen
						my $end_NO = ($flag & 0x40) ? 0 : 1;
						$RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO] = max( $RG_Info_Href->{$RG_ID}->{max_rlen}->[$end_NO], $rlen );
					}
					delete $pid_mapQ{$pid};
				}
			}
		}
	}
	close BAM;

	# close file-handle(s) for RG(s)
	for my $RG_ID (sort keys %$RG_Info_Href){
		close $RG_Info_Href->{$RG_ID}->{CandBam_wH};
		delete $RG_Info_Href->{$RG_ID}->{CandBam_wH};
	}

	# inform
	stout_and_sterr `date`."[INFO]\textract candidate reads from bam ok.\n";

	# close file handles of tmp files, and sweep them
	if($refseg_NO){
		for my $host_refseg (keys %$Host_Refseg_Href){
			if(exists($Host_Refseg_Href->{$host_refseg}->{tmp_file_wH})){
				close $Host_Refseg_Href->{$host_refseg}->{tmp_file_wH};
				delete $Host_Refseg_Href->{$host_refseg}->{tmp_file_wH};
				delete $Host_Refseg_Href->{$host_refseg}->{tmp_file};
			}
		}
		# sweep
		`rm -rf $V_Href->{CandFq_TmpFolder}`;
	}

	#-----------------------------#
	# to smart PE reads gzip file #
	#-----------------------------#
	for my $RG_ID (sort keys %$RG_Info_Href){

		my $RG_NO = $RG_Info_Href->{$RG_ID}->{RG_NO};
		my $candidate_bam = $RG_Info_Href->{$RG_ID}->{CandBam};
		my $candidate_bam_Nsorted = $RG_Info_Href->{$RG_ID}->{Cand_nSortBam};
		my $candidate_SmartFq_gz = $RG_Info_Href->{$RG_ID}->{CandFqGz};

		# sort reads in single fastq file for BWA mem smart pairing
		my $sort_tmp_folder = GetPath( filekey => 'Cand_nSortBam_TmpFolder', CandFq_RG_ID_NO => $RG_NO );
		`rm -rf $sort_tmp_folder`;
		`mkdir -p $sort_tmp_folder`;
		# sort bam by read-id
		stout_and_sterr `date`."[INFO]\tstarts sorting candidate bam for RG $RG_ID.\n";
		my $sort_command = "$V_Href->{samtools} sort --threads 4 -n -m 1G -T $sort_tmp_folder/sort_cand_bam -o $candidate_bam_Nsorted $candidate_bam";
		trible_run_for_success($sort_command, 'nSort_cand_BAM', {esdo_Nvb=>1});
		# sweep
		`rm -rf $sort_tmp_folder`;
		`rm -rf $candidate_bam`;
		delete $RG_Info_Href->{$RG_ID}->{CandBam};
		stout_and_sterr `date`."[INFO]\tsort candidate bam for RG $RG_ID ok.\n";

		# samtools to convert bam to fastq
		stout_and_sterr `date`."[INFO]\tstarts extracting fastq from nSorted candidate bam for RG $RG_ID.\n";
		&extract_fastq_from_nSortBam(
			{
				nSortBam => $candidate_bam_Nsorted,
				fastq_gz => $candidate_SmartFq_gz
			}
		);
		stout_and_sterr `date`."[INFO]\textract fastq from nSorted candidate bam for RG $RG_ID ok.\n";
	}

	# inform
	stout_and_sterr `date`."[INFO]\textract candidate PE reads to smart PE gz_file ok. from bam:\n\t $V_Href->{bam}\n";

	#-----------------------#
	# output RG_info_report #
	#-----------------------#
	open (RGRP,Try_GZ_Write($V_Href->{RG_info_report})) || die "fail $V_Href->{RG_info_report}: $!\n";
	print RGRP '#' . join("\t", ('RG_NO', 'RG_ID', 'RG_rID_pref', 'Max_rlen', 'RG_nSortBam', 'RG_smartFQ') ) . "\n"; # theme
	for my $RG_ID (sort { $RG_Info_Href->{$a}->{RG_NO} <=> $RG_Info_Href->{$b}->{RG_NO} } keys %$RG_Info_Href){
		my @max_rlen = @{ $RG_Info_Href->{$RG_ID}->{max_rlen} };
		if(min(@max_rlen) == 0){
			stout_and_sterr `date`."<WARN>\tRG $RG_ID contributes no candidate reads for virus analysis. Skip it.\n";
			delete $RG_Info_Href->{$RG_ID};
			next;
		}

		print RGRP      $RG_Info_Href->{$RG_ID}->{RG_NO};
		print RGRP "\t".$RG_ID;
		print RGRP "\t".$RG_Info_Href->{$RG_ID}->{RG_rID_pref};
		print RGRP "\t".join(',',@max_rlen);
		print RGRP "\t".basename($RG_Info_Href->{$RG_ID}->{Cand_nSortBam});
		print RGRP "\t".basename($RG_Info_Href->{$RG_ID}->{CandFqGz})."\n";
		# update whole max_rlen
		$V_Href->{max_rlen} = max( $V_Href->{max_rlen}, max(@max_rlen) );
	}
	close RGRP;
	stout_and_sterr `date`."[INFO]\twrite RG info report ok.\n";

	# run mode check, if not go on, just stop here.
	if( $V_Href->{run_mode} == 1 ){ # finish at this point in mode 1
		exit(0);
	}
}

#--- convert the non-host alignment to unmaps ---
sub convert_nonhost_alignment_to_unmap{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $self_Mrefseg = $Option_Href->{self_Mrefseg};
	my $mate_Mrefseg = $Option_Href->{mate_Mrefseg};
	my $bamline_Sref = $Option_Href->{bamline_Sref};
	my $bamline_Aref = $Option_Href->{bamline_Aref};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Host_Refseg_Href = $V_Href->{Host_Refseg_Href};

	#----------------#
	# common changes #
	#----------------#
	# discard the last '\n'
	$bamline_Aref->[-1] =~ s/\s$//;
	# modify flag, common changes
	$bamline_Aref->[1] -= 0x2 if( $bamline_Aref->[1] & 0x2 ); # minus 0x2 (P), means set 0x2 as 0
	# TLEN
	$bamline_Aref->[8] = 0; # TLEN

	# case 1, this end maps to non-host refseg
	if( !exists($Host_Refseg_Href->{$self_Mrefseg}) ){
		# reversed mapped to original info
		if( $bamline_Aref->[1] & 0x10 ){
			( $bamline_Aref->[9] = reverse uc($bamline_Aref->[9]) ) =~ tr/ACGT/TGCA/;
			$bamline_Aref->[10]  = reverse $bamline_Aref->[10];
		}
		# modify flag, common changes
		$bamline_Aref->[1] += 0x4 unless( $bamline_Aref->[1] & 0x4  ); # add   0x4 (u),  means set 0x4  as 1
		$bamline_Aref->[1] -= 0x10    if( $bamline_Aref->[1] & 0x10 ); # minus 0x10 (r), means set 0x10 as 0
		# this end becomes to unmapped
		$bamline_Aref->[4] = 0; # mapQ
		$bamline_Aref->[5] = '*'; # CIGAR
		# mate from host-refseg?
		if( exists($Host_Refseg_Href->{$mate_Mrefseg}) ){ # mate maps to host-refseg, keep it
			# modify the map info
			$bamline_Aref->[2] = $mate_Mrefseg; # mapped host-refseg of paired end
			$bamline_Aref->[3] = $bamline_Aref->[7]; # mapped position of paired end
			# modify the paired one
			$bamline_Aref->[6] = '='; # same host-refseg of paired end
			# Discard Tags
			for my $i ( 11 .. scalar(@$bamline_Aref)-1 ){
				$bamline_Aref->[$i] = '_DISCARD_' unless( $bamline_Aref->[$i] =~ /^(RG|XC):[iZ]:/ );
			}
		}
		else{ # mate is not host-refseg, so covert both end to unmapped
			# modify flag
			$bamline_Aref->[1] += 0x8 unless( $bamline_Aref->[1] & 0x8  ); # add   0x8  (U), means set 0x8 as 1
			$bamline_Aref->[1] -= 0x20    if( $bamline_Aref->[1] & 0x20 ); # minus 0x20 (R), means set 0x20 as 0
			# modify the map info
			$bamline_Aref->[2] = '*'; # map refseg
			$bamline_Aref->[3] = 0; # map position
			# modify the paired one
			$bamline_Aref->[6] = '*'; # map refseg of paired end
			$bamline_Aref->[7] = 0; # map position of paired end
			# Discard Tags
			for my $i ( 11 .. scalar(@$bamline_Aref)-1 ){
				$bamline_Aref->[$i] = '_DISCARD_' unless( $bamline_Aref->[$i] =~ /^RG:Z:/ );
			}
		}
	}
	# case 1, the mate end maps to non-host refseg
	elsif( ($mate_Mrefseg ne '=') && ($mate_Mrefseg ne '*') && !exists($Host_Refseg_Href->{$mate_Mrefseg}) ){
		# modify flag
		$bamline_Aref->[1] += 0x8 unless( $bamline_Aref->[1] & 0x8  ); # add   0x8  (U), means set 0x8 as 1
		$bamline_Aref->[1] -= 0x20    if( $bamline_Aref->[1] & 0x20 ); # minus 0x20 (R), means set 0x20 as 0
		# modify the paired one
		$bamline_Aref->[6] = $self_Mrefseg; # map refseg of paired end
		$bamline_Aref->[7] = $bamline_Aref->[3]; # map position of paired end
		# Discard Tags
		for my $i ( 11 .. scalar(@$bamline_Aref)-1 ){
			$bamline_Aref->[$i] = '_DISCARD_' if( $bamline_Aref->[$i] =~ /^(MQ|MC|Q2|R2):[iZ]:/ );
		}
	}

	# update the bam_line
	$$bamline_Sref = join("\t", grep !/^_DISCARD_$/, @$bamline_Aref) . "\n";
}

#--- extract smart fastq from nSortBam ---
## read length selection
## barcode concatenation
sub extract_fastq_from_nSortBam{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $nSortBam  = $Option_Href->{nSortBam};
	my $fastq_gz  = $Option_Href->{fastq_gz};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# check
	if( !-e $nSortBam ){
		warn_and_exit `date`."<ERROR>\tCannot find the nSortBam file.\n"
							."\t$nSortBam\n";
	}

	# samtools to convert bam to fastq
	## some reversed-compl mapped reads are not converted successfully, 2016-05-10
	# my $fastq_command = "$V_Href->{samtools} fastq $nSortBam | gzip -c > $fastq_gz";
	# trible_run_for_success($fastq_command, 'extract_fq', {esdo_Nvb=>1});

	my $action = $V_Href->{discard_DupReads} ? 'discard' : 'keep';
	stout_and_sterr `date`."[INFO]:\tWill $action duplicated paired-reads.\n";

	my %Paired_End;
	#   0x800: to avoid supplementary alignment
	# + 0x100: to avoid secondary alignment
	# = 0x900
	open (NSBAM,"$V_Href->{samtools} view -F 0x900 $nSortBam |") || die "fail samtools view $nSortBam: $!\n";
	open (SMFQ,Try_GZ_Write($fastq_gz)) || die "fail $fastq_gz: $!\n"; # smart fastq
	while (<NSBAM>){
		my ($pid, $flag, $rseq, $qual) = (split)[0,1,9,10];
		# duplication filter
		if(      $V_Href->{discard_DupReads}
			&&   $flag & 0x400
			&& !($flag & 0x8)
			&& !($flag & 0x4)
		){
			next;
		}
		# reversed-compl ?
		if( $flag & 0x10 ){
			( $rseq = reverse uc($rseq) ) =~ tr/ACGT/TGCA/;
			$qual = reverse $qual;
		}
		# store anyway
		my $endNO = ($flag & 0x40) ? 1 : 2;
		$Paired_End{$pid}{$endNO} = [ $rseq, $qual, ( length($rseq) >= $V_Href->{min_CandRead_length} ) ];
		# output? end '2' comes later.
		if( exists($Paired_End{$pid}{2}) && exists($Paired_End{$pid}{1}) ){
			# both meet minimum read-length criterion
			if ( $Paired_End{$pid}{2}->[2] && $Paired_End{$pid}{1}->[2] ){
				# barcode concatenation if possible
				my $barc_postfix = '';
				if(    defined($V_Href->{barcode_SamTag})
					&& length($V_Href->{barcode_SamTag}) != 0
					&& /\s$V_Href->{barcode_SamTag}(\S+)\s/
				  ){
					$barc_postfix = $V_Href->{barcode_pref} . $1;
				}
				# output fastq format
				for my $NO (1, 2){
					print SMFQ "\@$pid$barc_postfix/$NO\n";
					print SMFQ $Paired_End{$pid}{$NO}->[0] . "\n";
					print SMFQ '+' . "\n";
					print SMFQ $Paired_End{$pid}{$NO}->[1] . "\n";
				}
			}
			# sweep
			delete $Paired_End{$pid};
		}
	}
	close NSBAM;
	close SMFQ;
}

#--- load RG info report ---
sub load_RG_info_report{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $RG_Info_Href = $V_Href->{RG_Info_Href};

	open (RGRP,Try_GZ_Read($V_Href->{RG_info_report})) || die "fail read $V_Href->{RG_info_report}: $!\n";
	# theme
	(my $theme_line = lc(<RGRP>)) =~ s/^#//;
	my @theme_tag = split /\s+/, $theme_line;
	while(<RGRP>){
		next if(/^\#/);
		my @info = split;
		my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
		my @RG_max_rlen = split /,/, $rgOb{max_rlen};
		$RG_Info_Href->{ $rgOb{rg_id} } -> {RG_NO}    = $rgOb{rg_no};
		$RG_Info_Href->{ $rgOb{rg_id} } -> {RG_rID_pref} = $rgOb{rg_rid_pref};
		$RG_Info_Href->{ $rgOb{rg_id} } -> {max_rlen} = \@RG_max_rlen;
		$RG_Info_Href->{ $rgOb{rg_id} } -> {Cand_nSortBam} = GetPath( filekey => 'Cand_nSortBam', CandFq_RG_ID_NO => $rgOb{rg_no} );
		$RG_Info_Href->{ $rgOb{rg_id} } -> {CandFqGz}      = GetPath( filekey => 'CandFqGz', CandFq_RG_ID_NO => $rgOb{rg_no} );
		# update whole max_rlen
		$V_Href->{max_rlen} = max( $V_Href->{max_rlen}, max(@RG_max_rlen) );
		# check file existence
		if( !-e $RG_Info_Href->{ $rgOb{rg_id} }->{Cand_nSortBam} ){
			# warn_and_exit "<ERROR>\tCannot find nSortBam for ReadGroup (RG): $rgOb{rg_id}\n".
			# 			  "\tplease run command 'cand_reads'.\n";
		}
		if(    !-e $RG_Info_Href->{ $rgOb{rg_id} }->{CandFqGz}
			&& (   $V_Href->{run_mode} == 2   # PrepareInitialVirusRef
				|| $V_Href->{run_mode} == 3 ) # ReconstrcutVirusVariants
		){
			stout_and_sterr "<WARN>\tCannot find SmartFq for ReadGroup (RG): $rgOb{rg_id}. 'cand_virus' will construct it.\n";
		}
	}
	close RGRP;
	stout_and_sterr `date`."[INFO]\tload RG info report ok.\n";
}

#--- read host fai ---
sub read_host_fai{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Host_Refseg_Href = $V_Href->{Host_Refseg_Href};

	open (HFAI,$V_Href->{host_fai}) || die"fail read $V_Href->{host_fai}: $!\n";
	while(<HFAI>){
		my ($host_refseg) = (split)[0];
		$Host_Refseg_Href->{$host_refseg} = {};
	}
	close HFAI;
}

#--- get RG_ID from bam line ---
sub get_RG_ID_from_bam_line{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	if($_[0] =~ /RG:Z:(\S+)/){
		return $1;
	}
	else{
		warn_and_exit "<ERROR>\tcannot find RG_ID with RG:Z: tag from bam line:\n\t $_[0]\nfrom: $V_Href->{bam}\n";
	}
}

#--- get map quality of mate segment from bam line ---
sub get_MateMapQ_from_bam_line{
	if($_[0] =~ /MQ:i:(\d+)/){
		return $1;
	}
	else{
		return 'na';
	}
}

#--- get specfic prefix for read groups ---
sub get_RG_ID_prefix{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	return '_'.$V_Href->{Tool_Tag}.'rgNO'.$_[0].$V_Href->{Tool_Tag}.'_';
}

#--- 
1; ## tell the perl script the successful access of this module.
