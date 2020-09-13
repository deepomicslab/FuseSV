package FuseSV::Virus_Integ::LocalHaplotypeDetect::Allele_OB;

use strict;
use warnings;
use List::Util qw/ min max sum first /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr /;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::LocalMap_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Allele_OB';
#----- version --------
$VERSION = "0.12";
$DATE = '2018-05-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						get_NO
						is_valid
						get_UcycIDlistAref
						get_UcycOB
						get_LM_SimuSeqReads_Hash
						is_SimuLRbiEndSeq
						get_SimuSeqRlen
						LocalMapNO2OB
						get_LocalMapOBlist
						get_CN
						resolve_localmap
						prepare_LM_SimuSeqReads
						display
					 /;

#--- structure of object
# allele_OB -> NO = $NO
# allele_OB -> valid = 1/0, initilize as 1
# allele_OB -> CN = $CN
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB module
# allele_OB -> groupOB = $Group_OB
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Solution_OB module
# allele_OB -> solution_OB = $Solution_OB
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB module
# allele_OB -> UCYC -> { $UcycID } = { UcycOB => $UnitCycle_OB, count => $count }
# allele_OB -> FVGM_UCYC -> { $UcycID } = { UcycOB => $UnitCycle_OB, count => $count }
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::LocalMap_OB module
# allele_OB -> LocalMap -> { $LocalMap_NO => $LocalMap_OB }
# allele_OB -> LMsegStrMd5 -> { $segStrMd5 => $LM_NO }
# allele_OB -> att_Str -> $att_Str : general attributes string for output
## for simulated reads on local map
# allele_OB -> SimuSeqRead -> { $read_NO -> {stpos->$stpos, edpos->$edpos} }
# allele_OB -> SimuLRbiEndSeq = TURE/FALSE
# allele_OB -> SimuSeqRlen = $readLen;

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $NO = $parm{NO};
	my $CN = $parm{CN};
	my $UCYC_countStr = $parm{UCYC_countStr};
	my $FVGM_UCYC_countStr = $parm{FVGM_UCYC_countStr} || '';
	my $groupOB = $parm{groupOB};
	my $solution_OB = $parm{solution_OB};
	my $LineContent = $parm{LineContent}; # for output

	my $allele_OB = {};

	$allele_OB->{NO} = $NO;
	$allele_OB->{valid} = 1; # initilize as valid
	$allele_OB->{CN} = $CN;
	$allele_OB->{groupOB} = $groupOB;
	$allele_OB->{solution_OB} = $solution_OB;
	$allele_OB->{LineContent} = $LineContent;
	$allele_OB->{att_Str} = "SOLUTION_NO=".$solution_OB->get_NO.';'.
							"GROUP=".$groupOB->get_ID.';'.
							"ALLELE_NO=".$NO.';'.
							"ALLELE_CN=".$CN.';';

	# UCYC and FVGM_UCYC
	for my $Comb_Aref ( [      'UCYC',      $UCYC_countStr ],
						[ 'FVGM_UCYC', $FVGM_UCYC_countStr ] ){
		my ($key, $countStr) = @$Comb_Aref;
		$allele_OB->{$key} = {}; # initial
		for my $UcycInfo (split /[,;]/, $countStr){
			my ($UcycID, $count) = ($UcycInfo =~ /^(.+):(\d+)$/);
			my $UcycOB = $solution_OB->get_UcycOB( UcycID => $UcycID );
			# check
			$key =~ /FVGM/ && $UcycOB->get_fvgmID;
			# record
			$allele_OB->{$key}->{$UcycID} = { UcycOB => $UcycOB, count => $count };
		}
	}

	bless($allele_OB);
	return $allele_OB;
}

#--- return NO. of this allele
sub get_NO{
	my $allele_OB = shift;
	return $allele_OB->{NO};
}

#--- return valid status
sub is_valid{
	my $allele_OB = shift;
	return $allele_OB->{valid};
}

#--- return ref of list array of UCYC ID, see which key ('UCYC' or 'FVGM_UCYC')
sub get_UcycIDlistAref{
	my $allele_OB = shift;
	my %parm = @_;
	my $key = $parm{key} || 'UCYC';

	return [ sort keys %{$allele_OB->{$key}} ];
}

#--- return UCYC OB by UcycID, see which key ('UCYC' or 'FVGM_UCYC')
sub get_UcycOB{
	my $allele_OB = shift;
	my %parm = @_;
	my $key = $parm{key} || 'UCYC';
	my $UcycID = $parm{UcycID};

	return ( $allele_OB->{$key}->{$UcycID}->{UcycOB},
			 $allele_OB->{$key}->{$UcycID}->{count}  );
}

#--- return hash ref of simulated sequencing reads
sub get_LM_SimuSeqReads_Hash{
	my $allele_OB = shift;
	return $allele_OB->{SimuSeqRead};
}

#--- only get bi-end short reads ?
sub is_SimuLRbiEndSeq{
	my $allele_OB = shift;
	return $allele_OB->{SimuLRbiEndSeq};
}

#--- return read length of simulates reads
sub get_SimuSeqRlen{
	my $allele_OB = shift;
	return $allele_OB->{SimuSeqRlen};
}

#--- return Local Map OB according given Local Map NO
sub LocalMapNO2OB{
	my $allele_OB = shift;
	my %parm = @_;
	my $localmap_NO = $parm{localmap_NO};
	return $allele_OB->{LocalMap}->{$localmap_NO};
}

#--- return Local Map OB list array ref
sub get_LocalMapOBlist{
	my $allele_OB = shift;
	my @LMOBlist = map{ $allele_OB->{LocalMap}->{$_} } sort {$a<=>$b} keys %{$allele_OB->{LocalMap}};
	return \@LMOBlist;
}

#--- return copy number of this allele
sub get_CN{
	my $allele_OB = shift;
	return $allele_OB->{CN};
}

#--- construct localmap
sub resolve_localmap{
	my $allele_OB = shift;
	my %parm = @_;
	# belowed parms are for random mode
	my $randTime = $parm{randTime} || 0;
	my $segstrHref = $parm{segstrHref} || undef; # segstr(s) must have
	my $useFVGM  = $parm{useFVGM}  || 0;
	my $unitCnt  = $parm{unitCnt}  || 10;
	my $prbAdjt  = $parm{prbAdjt}  || 0.01;
	my $combMaxCnt = $parm{combMaxCnt} || 1; # min is 1
	my $combPrbGapT = $parm{combPrbGapT} || 2;

	# determine the local map amount
	my $localmap_amount = $randTime || 1;
	my $UCYCnumToMerge = scalar @{$allele_OB->get_UcycIDlistAref};
	# if only has one UCYC, just one local map
	if( $UCYCnumToMerge == 1 ){
		$localmap_amount = 1;
		stout_and_sterr "<WARN>\tJust run LocalMap process ONLY one time for allele NO.$allele_OB->{NO} as it has only one UCYC. (".
						"Solution NO.".$allele_OB->{solution_OB}->get_NO."; ".
						"Group ".$allele_OB->{groupOB}->get_ID.")\n";
	}

	# assign local-map amount to randProb layers
	my %randProbSetting = ( 1 => { prob=>[0], count=>1 } );
	if( $randTime != 0 ){
		## stratified random possibility
		my $unitPrb = sprintf "%.2f", 1.0 / $unitCnt;
		my @prbList = map { $unitPrb*$_-$prbAdjt } ( 1 .. $unitCnt ); # it's ascending order
		## combinations of uniq-prob
		my $combPrbGap = $combPrbGapT * $unitPrb;
		my $combIdx = scalar(keys %randProbSetting) + 1;
		my $probSum = 0; # count in prob=0
		for my $combCnt ( 1 .. $combMaxCnt ){
			if( $combCnt == 1 ){
				$randProbSetting{$combIdx++} = { prob=>[$_], count=>undef } for @prbList;
				$probSum += scalar(@prbList);
			}
			else{
				my $L1CombCnt = $combCnt - 1;
				for my $L1CombHref ( grep scalar(@{$_->{prob}}) == $L1CombCnt, values %randProbSetting ){
					next if( $L1CombHref->{prob}->[0] == 0 );
					my $prob_base = $L1CombHref->{prob}->[-1];
					for my $prob (@prbList){
						next if( $prob < $prob_base+$combPrbGap );
						$randProbSetting{$combIdx++} = { prob=>[@{$L1CombHref->{prob}},$prob], count=>undef };
						$probSum += $combCnt;
						$prob_base = $prob; # update
					}
				}
			}
		}
		## assign counts for each prob-combination
		for my $combHref (values %randProbSetting){
			my $combProbSum = scalar( @{$combHref->{prob}} );
			$combHref->{count} = $combHref->{prob}->[0] == 0 ? 1 : ( int( $localmap_amount * ($combProbSum / $probSum) ) || 1 );
		}
	}

	# form localmap in different randProbComb with different randTime
	my $total_time = 0;
	my $uniqLMcount = 0;
	my @cand_LM_NO = ( 1 .. $localmap_amount );
	while(    $uniqLMcount < $localmap_amount
		   && $total_time < $localmap_amount * 3
		   && !( $total_time != 0 && $randTime == 0 )
		   # && $total_time < $localmap_amount * 0.8 # debug for loop once
	){
		for my $combIdx (sort {$a<=>$b} keys %randProbSetting){
			# srand(time);
			my $combProbAref = $randProbSetting{$combIdx}->{prob};
			my $combDealTime = $randProbSetting{$combIdx}->{count};
			# update candidate localmap_NO
			if( $combDealTime > scalar(@cand_LM_NO) ){
				@cand_LM_NO = grep !exists $allele_OB->{LocalMap}->{$_}, ( 1 .. $localmap_amount );
			}
			my @localmap_NO = $combDealTime < scalar(@cand_LM_NO) ? @cand_LM_NO[ 0 .. $combDealTime-1 ] : @cand_LM_NO;
			@cand_LM_NO = @cand_LM_NO[ scalar(@localmap_NO) .. $#cand_LM_NO ];
			for my $localmap_NO ( @localmap_NO ){
				$total_time ++;
				# create local map object
				my @parm = ( NO => $localmap_NO, allele_OB => $allele_OB, groupOB => $allele_OB->{groupOB} );
				$allele_OB->{LocalMap}->{$localmap_NO} = FuseSV::Virus_Integ::LocalHaplotypeDetect::LocalMap_OB->new( @parm );
				# merge unit-cycle, construct the local map
				stout_and_sterr "[INFO]\t# $allele_OB->{att_Str}LOCALMAP_NO=$localmap_NO;RANDPROB=".join(',',@$combProbAref).";\n";
				@parm = ( source => $allele_OB->{groupOB}->get_sourceSegID, sink => $allele_OB->{groupOB}->get_sinkSegID, randProbAref => $combProbAref, useFVGM => $useFVGM );
				$allele_OB->{LocalMap}->{$localmap_NO}->merge_UnitCycle( @parm );
				# note that during merging, the local map might become non-valid
				if( ! $allele_OB->{LocalMap}->{$localmap_NO}->is_valid ){
					delete $allele_OB->{LocalMap}->{$localmap_NO};
				}
				else{ # then de-duplication based on md5 of segment structure
					# show local map segment structure
					stout_and_sterr "[INFO]\tSegStr:\n\t".$allele_OB->{LocalMap}->{$localmap_NO}->show_segStructure."\n";
					# check segstr must have
					if( defined $segstrHref ){
						my $all_match_bool = 1;
						my $segstrToWarn;
						for my $infoHref (values %$segstrHref){
							for my $regexType ( 'regex_fw', 'regex_rv' ){
								$all_match_bool &&= $allele_OB->{LocalMap}->{$localmap_NO}->hasSegStr(regex=>$infoHref->{$regexType});
								last if( $all_match_bool ); # just one of 'fw' and 'rv' matches is OK.
							}
							unless( $all_match_bool ){
								$segstrToWarn = $infoHref->{segstr};
								last;
							}
						}
						unless( $all_match_bool ){
							delete $allele_OB->{LocalMap}->{$localmap_NO};
							# inform
							stout_and_sterr "<WARN>\tdelete local-map NO.$localmap_NO due to not matching segstr ($segstrToWarn).\n";
						}
					}
					# skip process below
					next unless( exists $allele_OB->{LocalMap}->{$localmap_NO} );
					# md5 de-dup
					my $segStrMd5 = $allele_OB->{LocalMap}->{$localmap_NO}->get_segStrMd5;
					if( exists $allele_OB->{LMsegStrMd5}->{$segStrMd5} ){
						my $cmp_LM_NO = $allele_OB->{LMsegStrMd5}->{$segStrMd5};
						my $cmp_LM_OB = $allele_OB->{LocalMap}->{$cmp_LM_NO};
						if( $cmp_LM_OB->get_struHeight > $allele_OB->{LocalMap}->{$localmap_NO}->get_struHeight ){
							delete $allele_OB->{LocalMap}->{$cmp_LM_NO};
							$allele_OB->{LMsegStrMd5}->{$segStrMd5} = $localmap_NO;
						}
						else{
							delete $allele_OB->{LocalMap}->{$localmap_NO};
						}
					}
					else{
						$allele_OB->{LMsegStrMd5}->{$segStrMd5} = $localmap_NO;
						$uniqLMcount ++;
					}
				}
			}
			last if( $uniqLMcount >= $localmap_amount );
		}
	}
	# sweep
	delete $allele_OB->{LMsegStrMd5};
	# if no local-map found, mark non-valid
	$allele_OB->{valid} = 0 if( $uniqLMcount == 0 );
}

#--- simulate sequencing reads for calculating pair-wise links
## note that all local-maps of this allele has same $LM_SumLen,
## so record one set of reads' locations and use it on all local-maps
# for long-range PE sequencing, 10x or Hi-C
# note that, 1) MinIon or PacBio could use 10x mode.
#            2) Hi-C mode might need much more memory.
# define length of local-map is l.
# 1st, get a long-range length (L)
#      10x has prob.dist L of length; Hi-C just random in l.
# 2nd, get 5p/3p-edge position along local-map randomly
#      10x links all anchors in L; Hi-C just get bilateral reads.
# 3rd, based on the location of L, calculate the mean-depth of local-map
# do iteratively, till reaches the mena depth required
sub prepare_LM_SimuSeqReads{
	my $allele_OB = shift;
	my %parm = @_;
	my $localmap_NO = $parm{localmap_NO};
	my $seqType = $parm{seqType} || '10x';
	my $aimDepth = $parm{aimDepth} || 50;
	my $readLen = $parm{readLen} || 150; # for Hi-C, PE

	# take the first one local-map
	$localmap_NO ||= first {/^\d+$/} sort {$a<=>$b} keys %{$allele_OB->{LocalMap}};
	my $localmap_OB = $allele_OB->{LocalMap}->{$localmap_NO};
	if( !defined $localmap_OB ){
		warn_and_exit "<ERROR>\tCannot find the NO.$localmap_NO local map of allele NO.$allele_OB->{NO}.\n";
	}

	# only get bi-end short reads?
	my $LRbiEndSeq = 0;
	my $LRneedCovR = 0;
	if(    $seqType =~ /10x/i
		|| $seqType =~ /MinIon/i
		|| $seqType =~ /PacBio/i
	){
		$LRbiEndSeq = 0;
		$LRneedCovR = 1 if( $seqType =~ /10x/i );
	}
	elsif( $seqType =~ /Hi-C/i ){
		$LRbiEndSeq = 1;
	}
	else{
		warn_and_exit "<ERROR>\tseqType $seqType is not valid in localmap_OB->prepare_LM_SimuSeqReads func.\n";
	}

	# length of local map of this allele
	my $LMlen = $localmap_OB->get_length;
	# simulate sequencing
	my %SimuSeqRead;
	my $meanDepth = 0;
	my $read_NO = 0;
	while(   $meanDepth < $aimDepth
		  && $read_NO <= 5E6
	){
		$read_NO++;
		my $SimuSeqLen = $localmap_OB->{groupOB}->get_randLRLen( cand_StrLen => $LMlen );
		my $SimuSeqCovR = $LRneedCovR ? $localmap_OB->{groupOB}->get_randLRcovR( LRlen => $SimuSeqLen ) : 1;
		# print "read_NO:$read_NO, SimuSeqLen:$SimuSeqLen, SimuSeqCovR:$SimuSeqCovR\n";# debug
		# based on Seq length
		my $coverLMlen;
		if( $SimuSeqLen >= $LMlen ){
			my $randSelc = rand(1);
			if( $randSelc <= 0.33 ){ # full cover
				$SimuSeqRead{$read_NO} = {stpos=>1, edpos=>$LMlen};
				$coverLMlen = $LMlen;
			}
			else{
				my $randShift = int(rand($LMlen/2))+1;
				if( $randSelc >= 0.67 ){ # rand stpos
					$SimuSeqRead{$read_NO} = {stpos=>$randShift, edpos=>$LMlen};
					$coverLMlen = $LMlen - $randShift + 1;
				}
				else{ # rand edpos
					$SimuSeqRead{$read_NO} = {stpos=>1, edpos=>$LMlen-$randShift+1};
					$coverLMlen = $LMlen - $randShift + 1;
				}
			}
		}
		else{
			my $randShift = int(rand($LMlen-$SimuSeqLen))+1;
			$SimuSeqRead{$read_NO} = {stpos=>$randShift, edpos=>$randShift+$SimuSeqLen-1};
			$coverLMlen = $SimuSeqLen;
		}
		# record cover ratio
		$SimuSeqRead{$read_NO}{covR} = $SimuSeqCovR;
		# update
		if( $LRbiEndSeq ){
			$meanDepth += $readLen * 2 / $LMlen;
		}
		else{
			$meanDepth += $SimuSeqCovR * $coverLMlen/ $LMlen;
		}
	}

	# record SimuSeqRead
	$allele_OB->{SimuSeqRead} = \%SimuSeqRead;
	$allele_OB->{SimuLRbiEndSeq} = $LRbiEndSeq;
	$allele_OB->{SimuSeqRlen} = $readLen;

	# inform
	my $group_ID = $localmap_OB->{groupOB}->get_ID;
	stout_and_sterr "[INFO]\tGroup $group_ID Allele NO.$allele_OB->{NO}\n".
						  "\tGet simulated $read_NO $seqType long-range HMW/reads OK.\n";

	return $allele_OB->{SimuSeqRead}; # sometimes, necessary
}

#--- display allele's info and local-map(s)
sub display{
	my $allele_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID};

	my $return_Str = $allele_OB->{LineContent} . "\n"; # original contents from config
	# each local map
	for my $localmap_NO (sort {$a<=>$b} keys %{$allele_OB->{LocalMap}}){
		my $localmap_OB = $allele_OB->{LocalMap}->{$localmap_NO};
		# attributes to output
		$return_Str .=  "LOCALMAP\t".
						"NO=$localmap_NO;".
						"RANDPROB=".join(',',@{$localmap_OB->get_randProbAref}).";".
						$allele_OB->{att_Str}.
						"FVGMUSED=".$localmap_OB->get_fvgmInfoStr.";".
						"SEGSTRMD5=".$localmap_OB->get_segStrMd5.";".
						"SORTSEGSTRMD5=".$localmap_OB->get_nosSegStrMd5.";".
						"STRUDEPTH=".$localmap_OB->get_structureDepth.";\n";
		# human-reading structure
		$return_Str .= $localmap_OB->display( prefix => "STRUCTURE:\t", sayID => $sayID );
	}

	return $return_Str;
}

1; ## tell the perl script the successful access of this module.
