package FuseSV::Virus_Integ::LocalHaplotypeDetect::Solution_OB;

use strict;
use warnings;
use File::Temp;
use IO::File;
use Parallel::ForkManager;
use List::Util qw/ min max sum /;
use File::Basename;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Stat::PearsonCorr qw/ get_pearson_corr /;
use BioFuse::Stat::SpearmanCorr qw/ get_spearman_corr /;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Allele_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Solution_OB';
#----- version --------
$VERSION = "0.07";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						load_allele
						get_NO
						is_valid
						get_UcycOB
						resolve_localmap_of_allele
						display
						pwanchor_report
					 /;

#--- structure of object
# solution_OB -> NO = $NO
# solution_OB -> valid = 1/0, initilize as 1
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB module
# solution_OB -> groupOB = $Group_OB
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB module
# solution_OB -> UCYC -> { $UcycID } = { UcycOB => $UnitCycle_OB, count => $count }
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Allele_OB module
# solution_OB -> allele -> = { $allele_NO => $Allele_OB }
## local-map combinations
# solution_OB -> LocalMapComb -> { $LM_Comb_ID => {ID=>$ID, Corr=>$Corr, pwanchor=>{$pwanchorID=>{anchor_a=>,anchor_b=>,simuLC=>,userLC=>}}}, .. }
## corr
# solution_OB -> CorrMethod = 'Pearson' or 'Spearman'

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $NO = $parm{NO};
	my $UCYC_countStr = $parm{UCYC_countStr};
	my $groupOB = $parm{groupOB};
	my $UCYC_Href = $parm{UCYC_Href};

	my $solution_OB = {};

	$solution_OB->{NO} = $NO;
	$solution_OB->{valid} = 1; # initilize as valid
	$solution_OB->{groupOB} = $groupOB;
	# UCYC
	for my $UcycInfo (split /[,;]/, $UCYC_countStr){
		my ($UcycID, $count) = ($UcycInfo =~ /^(.+):(\d+)$/);
		# check existence
		if( !exists $UCYC_Href->{$UcycID} ){
			warn_and_exit "<ERROR>\tCannot recognize unit-cycle $UcycID from solution $NO of group " . $groupOB->get_ID . ":\n"
								."\t$UCYC_countStr\n";
		}
		# record
		$solution_OB->{UCYC}->{$UcycID} = { UcycOB => $UCYC_Href->{$UcycID}, count => $count };
	}

	bless($solution_OB);
	return $solution_OB;
}

#--- load an allele for this solution
sub load_allele{
	my $solution_OB = shift;
	my %parm = @_;
	my $allele_NO = $parm{allele_NO};
	my $allele_cn = $parm{allele_cn};
	my $allele_UCYC_countStr = $parm{allele_UCYC_countStr};
	my $allele_FVGM_UCYC_countStr = $parm{allele_FVGM_UCYC_countStr} || '';
	my $LineContent = $parm{LineContent}; # for output

	my @parm = ( NO => $allele_NO, UCYC_countStr => $allele_UCYC_countStr,
				 CN => $allele_cn, FVGM_UCYC_countStr => $allele_FVGM_UCYC_countStr,
				 groupOB => $solution_OB->{groupOB}, solution_OB => $solution_OB     );
	$solution_OB->{allele}->{$allele_NO} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Allele_OB->new( @parm, LineContent => $LineContent );
}

#--- return solution NO
sub get_NO{
	my $solution_OB = shift;
	return $solution_OB->{NO};
}

#--- return valid status
sub is_valid{
	my $solution_OB = shift;
	return $solution_OB->{valid};
}

#--- get the UCYC OB by UcycID
sub get_UcycOB{
	my $solution_OB = shift;
	my %parm = @_;
	my $UcycID = $parm{UcycID};

	if( !exists $solution_OB->{UCYC}->{$UcycID} ){
		warn_and_exit "<ERROR>\tCannot find unit-cycle $UcycID from solution $solution_OB->{NO} of group " . $solution_OB->{groupOB}->get_ID . "!\n";
	}

	return $solution_OB->{UCYC}->{$UcycID}->{UcycOB};
}

#--- construct localmap(s) of each allele
sub resolve_localmap_of_allele{
	my $solution_OB = shift;
	my %parm = @_;
	my $randTime = $parm{randTime} || 0; # random mode
	my $combMaxCnt = $parm{combMaxCnt} || 1; # random mode, min is 1
	my $segstrHref = $parm{segstrHref} || undef; # segstr(s) must have
	my $useFVGM  = $parm{useFVGM}  || 0; # random mode
	my $corrMethod = $parm{corrMethod} || 'p'; # 'p' for pearson, 's' for spearman
	# for long-range simu.seq work
	my $seqType = $parm{seqType} || '10x';
	my $aimDepth = $parm{aimDepth} || 50;
	my $readLen = $parm{readLen} || 150; # for Hi-C, PE
	# for fork
	my $forkNum = $parm{forkNum} || 1;
	   $forkNum = 1 if( $randTime == 0 ); # reset
	my $tmpPref = $parm{tmpPref};

	my @allele_NO = sort {$a<=>$b} keys %{$solution_OB->{allele}};
	# construct local map
	for my $i (0 .. $#allele_NO){
		my $allele_NO = $allele_NO[$i];
		my $allele_OB = $solution_OB->{allele}->{$allele_NO};
		# construct local map
		my @parm = ( randTime => $randTime, combMaxCnt => $combMaxCnt, useFVGM => $useFVGM );
		# add required segstr
		if(    defined $segstrHref # supplied by user
			&& scalar @{$allele_OB->get_UcycIDlistAref} != 1 # normally, sinlge UCYC is ref-allele
		){
			push @parm, ( segstrHref => $segstrHref );
		}
		# construct local map
		$allele_OB->resolve_localmap( @parm );
		# check allele valid
		unless( $allele_OB->is_valid ){
			# inform
			stout_and_sterr "<WARN>\tas allele NO.$allele_NO is non-valid, set solution NO.$solution_OB->{NO} non-valid.\n";
			delete $solution_OB->{allele}->{$allele_NO};
			$solution_OB->{valid} = 0;
			last;
		}
	}
	# not valid!
	return unless( $solution_OB->is_valid );

	# long-range simu.seq work ?
	my $groupOB = $solution_OB->{groupOB};
	return unless( $groupOB->do_LM_Comb_LR_Corr );
	my $groupID = $groupOB->get_ID;

	my $LM_SumLen = 0;
	my $LM_Comb = {};
	for my $i (0 .. $#allele_NO){
		my $allele_NO = $allele_NO[$i];
		my $allele_OB = $solution_OB->{allele}->{$allele_NO};
		my $LMOBlist_Aref = $allele_OB->get_LocalMapOBlist;
		# LM length sum
		$LM_SumLen += $LMOBlist_Aref->[0]->get_length * $allele_OB->get_CN;
		# local map combination
		if( $i == 0 ){
			push @{$LM_Comb->{$allele_NO}}, [$_] for @$LMOBlist_Aref;
		}
		else{
			for my $lastCombAref (@{$LM_Comb->{$allele_NO[$i-1]}}){
				push @{$LM_Comb->{$allele_NO}}, [@$lastCombAref, $_] for @$LMOBlist_Aref;
			}
			delete $LM_Comb->{$allele_NO[$i-1]};
		}
		# simu.seq to calculate pwanchor link counts
		my @parm = ( seqType => $seqType, aimDepth => $aimDepth, readLen => $readLen );
		# simulate reads for this allele
		## note that all local-maps of this allele has same $LM_SumLen,
		## so record one set of reads' locations and use it on all local-maps
		$allele_OB->prepare_LM_SimuSeqReads( @parm );
		# split local map OB based on fork number
		my $LMOB_Count = scalar(@$LMOBlist_Aref);
		my $unitForkLMnum = max( int( $LMOB_Count / $forkNum ), 1 );
		my %ForkIdxRange;
		for (my $stI = 0; $stI <= $LMOB_Count-1; ){
			$ForkIdxRange{$stI} = min( $stI + $unitForkLMnum - 1, $LMOB_Count - 1 );
			$stI += $unitForkLMnum;
		}
		# for fork temp file
		my $tmpTEMPLATE = basename($tmpPref).".AL$allele_NO.SimuPW_";
		my $tmpDIR = dirname($tmpPref);
		`rm -f $tmpDIR/$tmpTEMPLATE*`; # sweep
		# inform
		stout_and_sterr `date`."[INFO]\tcalculating pair-wise anchors simulated links. allele NO.$allele_NO; solution NO.$solution_OB->{NO}; group $groupID.\n";
		# fork manager
		my $pm;
		if($forkNum > 1){ $pm = new Parallel::ForkManager($forkNum) }
		for my $stI (sort {$a<=>$b} keys %ForkIdxRange){
			my $edI = $ForkIdxRange{$stI};
			# unit fork job
			if($forkNum > 1){ $pm->start and next; }
			# prepare temp file
			my $tmp_fh = File::Temp->new( TEMPLATE => $tmpTEMPLATE.'XXXX', UNLINK => 0, DIR => $tmpDIR );
			# simu.pwanchor
			for my $idx ($stI .. $edI){
				my $localmap_OB = $LMOBlist_Aref->[$idx];
				my $localmap_NO = $localmap_OB->get_NO;
				$localmap_OB->calculate_SimuPWanchorLink; # now SCCC work, 2018
				# record simu.pwanchor in temp file
				my $LM_SimuLinkHref = $localmap_OB->get_PWanchorSimuLink;
				print $tmp_fh "$localmap_NO\t$_\t$LM_SimuLinkHref->{$_}\n" for keys %$LM_SimuLinkHref;
			}
			close $tmp_fh;
			if($forkNum > 1){ $pm->finish; }
		}
		if($forkNum > 1){ $pm->wait_all_children; }
		# load all temp file
		if($forkNum > 1){
			stout_and_sterr `date`."[INFO]\tload tmpFiles to fetch pair-wise anchor simulated links.\n";
			my @tmpFile = glob "$tmpDIR/$tmpTEMPLATE*";
			for my $tmpFile (@tmpFile){
				open(TMP, Try_GZ_Read($tmpFile)) || die "fail to read simu.pwanchor tmp file: $!\n";
				while(<TMP>){
					my ($localmap_NO, $pwanchorID, $pwanchorLinkC) = (split)[0,1,2];
					my $localmap_OB = $allele_OB->LocalMapNO2OB( localmap_NO => $localmap_NO );
					$localmap_OB->load_PWanchorSimuLink( pwanchorID => $pwanchorID, pwanchorLinkC => $pwanchorLinkC );
				}
				close TMP;
				# inform
				stout_and_sterr "\tload $tmpFile OK.\n";
			}
		}
		# sweep
		`rm -f $tmpDIR/$tmpTEMPLATE*`;
	}

	$solution_OB->{CorrMethod} = ( $corrMethod =~ /^p$/i ? 'Pearson' : 'Spearman' );
	# inform
	stout_and_sterr `date`."[INFO]\tcalculating $solution_OB->{CorrMethod} Corr of pair-wise anchors links. solution NO.$solution_OB->{NO}; group $groupID.\n";
	# gather simu.link of pair-wise anchor in each local map combination
	for my $LM_Comb_Aref (@{$LM_Comb->{$allele_NO[-1]}}){
		my %PWanchorSimuLink;
		my @LM_Comb_ID;
		for my $i (0 .. $#allele_NO){
			my $allele_NO = $allele_NO[$i];
			my $allele_OB = $solution_OB->{allele}->{$allele_NO};
			my $LM_OB = $LM_Comb_Aref->[$i];
			push @LM_Comb_ID, "AL$allele_NO.LM".$LM_OB->get_NO;
			# summarize all pwanchor link counts
			my $ratio = ($LM_OB->get_length * $allele_OB->get_CN) / $LM_SumLen;
			my $LM_SimuLinkHref = $LM_OB->get_PWanchorSimuLink;
			for my $pwanchorID (keys %$LM_SimuLinkHref){
				$PWanchorSimuLink{$pwanchorID} += sprintf "%.2f", $LM_SimuLinkHref->{$pwanchorID} * $ratio;
			}
		}
		my $LM_Comb_ID = join(',', @LM_Comb_ID);
		# record
		$solution_OB->{LocalMapComb}->{$LM_Comb_ID} = { ID => $LM_Comb_ID };
		my $LM_Comb_Href = $solution_OB->{LocalMapComb}->{$LM_Comb_ID};
		my %LinkCount;
		for my $pwanchorID (sort keys %PWanchorSimuLink){
			my $simuLinkCount = $PWanchorSimuLink{$pwanchorID};
			push @{$LinkCount{simu}}, $simuLinkCount;
			my $pwanchorOB = $groupOB->PWanchorID2OB( pwanchorID => $pwanchorID );
			my $userLinkCount = $pwanchorOB->get_linkCount;
			push @{$LinkCount{user}}, $userLinkCount;
			$LM_Comb_Href->{pwanchor}->{$pwanchorID} = { anchor_a => $pwanchorOB->get_sideInfo(side=>'a'),
														 anchor_b => $pwanchorOB->get_sideInfo(side=>'b'),
														 simuLC => $simuLinkCount, userLC => $userLinkCount };
		}
		# Pearson/Spearman correlation coefficient
		my @corr_parm = ( list_a_Aref => $LinkCount{user}, list_b_Aref => $LinkCount{simu}, error_getNA => 1 );
		$LM_Comb_Href->{Corr} = ( $corrMethod =~ /^p$/i ? get_pearson_corr( @corr_parm ) : get_spearman_corr( @corr_parm ) );
	}
}

#--- display alleles' info and local-map(s)
sub display{
	my $solution_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID};

	my $return_Str = '';
	for my $allele_NO (sort keys %{$solution_OB->{allele}}){
		my $allele_OB = $solution_OB->{allele}->{$allele_NO};
		$return_Str .= $allele_OB->display( sayID => $sayID );
	}

	return $return_Str;
}

#--- generate pair-wise anchor report
sub pwanchor_report{
	my $solution_OB = shift;

	my $groupOB = $solution_OB->{groupOB};
	my $groupID = $groupOB->get_ID;
	my $LocalMapComb_Href = $solution_OB->{LocalMapComb};

	my @LM_Comb_ID = grep $LocalMapComb_Href->{$_}->{Corr} ne 'N/A', keys %$LocalMapComb_Href;
	   @LM_Comb_ID = sort { $LocalMapComb_Href->{$b}->{Corr} <=> $LocalMapComb_Href->{$a}->{Corr} } @LM_Comb_ID;
	push @LM_Comb_ID, $_ for grep $LocalMapComb_Href->{$_}->{Corr} eq 'N/A', keys %$LocalMapComb_Href;

	my $return_Str = '';
	my @theme = qw/ anchor_a anchor_b userLC simuLC /;
	for my $LM_Comb_ID ( @LM_Comb_ID ){
		my $LM_Comb_Href = $solution_OB->{LocalMapComb}->{$LM_Comb_ID};
		$return_Str .= '#'."Group:$groupID\tSolution:$solution_OB->{NO}\t$LM_Comb_ID\t$solution_OB->{CorrMethod}Corr=$LM_Comb_Href->{Corr}\n";
		$return_Str .= '#'.join("\t", @theme)."\n";
		for my $pwanchorID (sort keys %{$LM_Comb_Href->{pwanchor}}){
			my @info = map{ $LM_Comb_Href->{pwanchor}->{$pwanchorID}->{$_} } @theme;
			$return_Str .= join("\t", @info)."\n";
		}
	}

	return $return_Str;
}

1; ## tell the perl script the successful access of this module.
