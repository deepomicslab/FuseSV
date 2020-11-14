package FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Solution_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Group_OB';
#----- version --------
$VERSION = "0.09";
$DATE = '2018-05-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						load_SegOB
						load_SourceSegOB
						load_SinkSegOB
						load_JuncOB
						load_solution_allele
						load_PWanchorOB
						load_LRlenProbDist
						ProbToProbGrid
						get_ID
						is_valid
						get_sourceSegID
						get_sinkSegID
						get_PWanchorOB
						PWanchorID2OB
						segID2segOB
						do_LM_Comb_LR_Corr
						get_prob_distance
						get_randLRLen
						get_randLRcovR
						resolve_localmap_of_solution
						display
						pwanchor_report
					 /;

#--- structure of object
# group_OB -> ID = $ID
# group_OB -> valid = 1/0, initilize as 1
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB module
# group_OB -> segOB = { $segID_1 -> $segOB_1, ..  }
# group_OB -> sourceSegOB = $sourceSegOB
# group_OB -> sinkSegOB = $sinkSegOB
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB module
# group_OB -> juncOB = { $juncID_1 -> $juncOB_1, .. }
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Solution_OB module
# group_OB -> solution = { 1 => $Solution_OB, .. }
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB module
# group_OB -> pwanchorOB = { $pwanchorID_1 -> $pwanchorOB_1, .. }
# group_OB -> LRlenToAccuProb = { $LRlen -> $accuProb, .. }
# group_OB -> ProbGridToLRlen = { $ProbGrid -> [ $LRlen_1, $LRlen_2, .. ], .. }

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $ID = $parm{ID};

	my $group_OB = {};

	$group_OB->{ID} = $ID;
	$group_OB->{valid} = 1; # initilize as valid
	$group_OB->{segOB} = {};
	$group_OB->{sourceSegOB} = undef;
	$group_OB->{sinkSegOB} = undef;
	$group_OB->{juncOB} = {};
	$group_OB->{solution} = {};
	$group_OB->{pwanchorOB} = {};
	$group_OB->{LRlenToAccuProb} = undef;
	$group_OB->{ProbGridToLRlen} = undef;

	bless($group_OB);
	return $group_OB;
}

#--- load segment OB
sub load_SegOB{
	my $group_OB = shift;
	my %parm = @_;
	my $seg_OB = $parm{seg_OB};
	$group_OB->{segOB}->{$seg_OB->get_ID} = $seg_OB;
}

#--- load source segment OB
sub load_SourceSegOB{
	my $group_OB = shift;
	my %parm = @_;
	$group_OB->{sourceSegOB} = $parm{seg_OB};
}

#--- load sink segment OB
sub load_SinkSegOB{
	my $group_OB = shift;
	my %parm = @_;
	$group_OB->{sinkSegOB} = $parm{seg_OB};
}

#--- load junction OB
sub load_JuncOB{
	my $group_OB = shift;
	my %parm = @_;
	my $junc_OB = $parm{junc_OB};
	$group_OB->{juncOB}->{$junc_OB->get_ID} = $junc_OB;
}

#--- load allele in solution
sub load_solution_allele{
	my $group_OB = shift;
	my %parm = @_;
	my $solution_NO = $parm{solution_NO};
	my $solution_UCYC_countStr = $parm{solution_UCYC_countStr};
	my $allele_NO = $parm{allele_NO};
	my $allele_cn = $parm{allele_cn};
	my $allele_UCYC_countStr = $parm{allele_UCYC_countStr};
	my $allele_FVGM_UCYC_countStr = $parm{allele_FVGM_UCYC_countStr} || '';
	my $UCYC_Href = $parm{UCYC_Href};
	my $LineContent = $parm{LineContent}; # for output

	# create solution object if possible
	if( !exists $group_OB->{solution}->{$solution_NO} ){
		my @parm = ( NO => $solution_NO, UCYC_countStr => $solution_UCYC_countStr, groupOB => $group_OB );
		$group_OB->{solution}->{$solution_NO} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Solution_OB->new( @parm, UCYC_Href => $UCYC_Href );
	}
	# load alleles to solution
	my @parm = ( allele_NO => $allele_NO, allele_UCYC_countStr => $allele_UCYC_countStr,
				 allele_cn => $allele_cn, allele_FVGM_UCYC_countStr => $allele_FVGM_UCYC_countStr );
	$group_OB->{solution}->{$solution_NO}->load_allele( @parm, LineContent => $LineContent );
}

#--- load pair-wise anchor OB
sub load_PWanchorOB{
	my $group_OB = shift;
	my %parm = @_;
	my $pwanchorOB = $parm{pwanchorOB};
	$group_OB->{pwanchorOB}->{$pwanchorOB->get_ID} = $pwanchorOB;
}

#--- load long-range length prob distribution
sub load_LRlenProbDist{
	my $group_OB = shift;
	my %parm = @_;

	# LRlen to accuProb
	$group_OB->{LRlenToAccuProb} = $parm{LRlenToAccuProb_Href};
	# Prob in Grid to Rlen array, for quick query
	my ($maxProbGrid, $minProbGrid) = (-1, 1E6);
	for my $LRlen (keys %{$group_OB->{LRlenToAccuProb}}){
		my $accuProb = $group_OB->{LRlenToAccuProb}->{$LRlen};
		my $accuProbGrid = $group_OB->ProbToProbGrid( Prob => $accuProb );
		push @{$group_OB->{ProbGridToLRlen}->{$accuProbGrid}}, $LRlen;
		$maxProbGrid = max($maxProbGrid, $accuProbGrid);
		$minProbGrid = min($minProbGrid, $accuProbGrid);
	}
	# make up blanks
	## start blanks
	$group_OB->{ProbGridToLRlen}->{$_} = $group_OB->{ProbGridToLRlen}->{$minProbGrid} for ( 0 .. $minProbGrid-1 );
	## end blanks
	$group_OB->{ProbGridToLRlen}->{$_} = $group_OB->{ProbGridToLRlen}->{$maxProbGrid} for ( $maxProbGrid+1 .. $group_OB->ProbToProbGrid( Prob => 1 ) );
	## inner blanks
	my $last_LRlen_Aref;
	for my $ProbGrid ( $minProbGrid .. $maxProbGrid ){
		if( !exists $group_OB->{ProbGridToLRlen}->{$ProbGrid} ){
			$group_OB->{ProbGridToLRlen}->{$ProbGrid} = $last_LRlen_Aref;
		}
		else{
			$last_LRlen_Aref = $group_OB->{ProbGridToLRlen}->{$ProbGrid};
		}
	}

	# LRlen to cover ratio
	$group_OB->{LRlen2ratioA} = $parm{LRlen2ratioA_Href};
}

#--- convert Prob to ProbGrid
sub ProbToProbGrid{
	my $group_OB = shift;
	my %parm = @_;
	my $Prob = $parm{Prob};

	return int($Prob * 10000);
}

#--- return group ID
sub get_ID{
	my $group_OB = shift;
	return $group_OB->{ID};
}

#--- return valid status
sub is_valid{
	my $group_OB = shift;
	return $group_OB->{valid};
}

#--- return source segment ID
sub get_sourceSegID{
	my $group_OB = shift;
	return $group_OB->{sourceSegOB}->get_ID;
}

#--- return sink segment ID
sub get_sinkSegID{
	my $group_OB = shift;
	return $group_OB->{sinkSegOB}->get_ID;
}

#--- return ref of hash contains pwanchorOB
sub get_PWanchorOB{
	my $group_OB = shift;
	return $group_OB->{pwanchorOB};
}

#--- return pwanchorOB based on given pwanchorID
sub PWanchorID2OB{
	my $group_OB = shift;
	my %parm = @_;
	my $pwanchorID = $parm{pwanchorID};
	return $group_OB->{pwanchorOB}->{$pwanchorID};
}

#--- return segOB based on given segID
sub segID2segOB{
	my $group_OB = shift;
	my %parm = @_;
	my $segID = $parm{segID};
	return $group_OB->{segOB}->{$segID};
}

#--- check whether to do local-map combination to test corr with long-range data
## test number of pwanchorOB
## 0: no; 1: could calculate the simulated distance; >=2: corr work
sub do_LM_Comb_LR_Corr{
	my $group_OB = shift;
	return (scalar(keys %{$group_OB->{pwanchorOB}}) != 0);
}

#--- calculate prob sum of given interval
sub get_prob_distance{
	my $group_OB = shift;
	my %parm = @_;
	my $from = $parm{from};
	my $to = $parm{to};

	if( defined $group_OB->{LRlenToAccuProb} ){
		my $from_accuProb = 0;
		my $to_accuProb = 1;
		for my $LRlen (sort {$a<=>$b} keys %{$group_OB->{LRlenToAccuProb}}){
			if( $LRlen <= $from ){
				$from_accuProb = $group_OB->{LRlenToAccuProb}->{$LRlen};
			}
			if( $LRlen >= $to ){
				$to_accuProb = $group_OB->{LRlenToAccuProb}->{$LRlen};
				last;
			}
		}

		return $to_accuProb - $from_accuProb; # enclosed prob area
	}
	else{
		return 1000/$from; # if no ProbDist, use reciprocal 'from' distance
	}
}

#--- return rand long-range length based on prob.dist
# min: 500; max: 750000
sub get_randLRLen{
	my $group_OB = shift;
	my %parm = @_;
	my $cand_StrLen = $parm{cand_StrLen} || 1E4;
	my $min_Len = $parm{min_Len} || 500;
	my $max_Len = $parm{max_Len} || 750000;

	my $Prob = rand(1);
	my $randLRlen;
	if( exists $group_OB->{ProbGridToLRlen} ){
		my $ProbGrid = $group_OB->ProbToProbGrid( Prob => $Prob );
		my $LRlen_Aref = $group_OB->{ProbGridToLRlen}->{$ProbGrid};
		my $LRlenCandC = scalar(@$LRlen_Aref);
		if( $LRlenCandC == 1 ){
			$randLRlen = $LRlen_Aref->[0];
		}
		else{
			$randLRlen = $LRlen_Aref->[int(rand($LRlenCandC))];
		}
	}
	else{
		$randLRlen = int($cand_StrLen * $Prob);
	}

	return min( max( $min_Len, $randLRlen ), $max_Len );
}

#--- return rand cover ratio of long-range DNA based on prob.dist
# ONLY use for 10x, so far
sub get_randLRcovR{
	my $group_OB = shift;
	my %parm = @_;
	my $LRlen = $parm{LRlen} || 1E4;

	my $match_LRlen;
	if( exists $group_OB->{LRlen2ratioA}->{$LRlen} ){
		$match_LRlen = $LRlen;
	}
	else{
		my $last_SmallLRlen = undef;
		for my $cand_LRlen (sort {$a<=>$b} keys %{$group_OB->{LRlen2ratioA}}){
			if($cand_LRlen < $LRlen){
				$last_SmallLRlen = $cand_LRlen;
			}
			else{
				# more near
				$match_LRlen = abs($cand_LRlen-$LRlen) < abs($LRlen-$last_SmallLRlen) ? $cand_LRlen : $last_SmallLRlen;
				last;
			}
		}
		# if no LagerLRlen
		$match_LRlen ||= $last_SmallLRlen;
	}

	my $match_ratioAref = $group_OB->{LRlen2ratioA}->{$match_LRlen};
	return $match_ratioAref->[int(rand(scalar(@$match_ratioAref)))];
}

#--- construct localmap(s) of each solution
sub resolve_localmap_of_solution{
	my $group_OB = shift;
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
	my $tmpPref = $parm{tmpPref};

	# construct local-map for each allele
	for my $solution_NO (sort {$a<=>$b} keys %{$group_OB->{solution}}){
		my $solution_OB = $group_OB->{solution}->{$solution_NO};
		my @parm = ( randTime => $randTime, combMaxCnt => $combMaxCnt, segstrHref => $segstrHref,
					 useFVGM => $useFVGM, seqType => $seqType, corrMethod => $corrMethod,
					 aimDepth => $aimDepth, readLen => $readLen, forkNum => $forkNum, tmpPref => $tmpPref );
		$solution_OB->resolve_localmap_of_allele( @parm );
		# check valid
		unless( $solution_OB->is_valid ){
			# inform
			stout_and_sterr "<WARN>\tdelete non-valid solution NO.$solution_OB->{NO} from group NO.$group_OB->{ID}.\n";
			delete $group_OB->{solution}->{$solution_NO};
		}
	}
	# mark non-valid if has no valid solution
	if( scalar(keys %{$group_OB->{solution}}) == 0 ){
		# inform
		stout_and_sterr "<WARN>\tset group NO.$group_OB->{ID} non-valid, as it has no valid solution.\n";
		$group_OB->{valid} = 0;
	}
}

#--- display solution + allele's localmap(s)
sub display{
	my $group_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID};

	my $return_Str = '';
	for my $solution_NO (sort keys %{$group_OB->{solution}}){
		my $solution_OB = $group_OB->{solution}->{$solution_NO};
		$return_Str .= $solution_OB->display( sayID => $sayID );
	}

	return $return_Str;
}

#--- generate pair-wise anchor report
sub pwanchor_report{
	my $group_OB = shift;

	my $return_Str = '';
	for my $solution_NO (sort keys %{$group_OB->{solution}}){
		my $solution_OB = $group_OB->{solution}->{$solution_NO};
		$return_Str .= $solution_OB->pwanchor_report;
	}

	return $return_Str;
}

1; ## tell the perl script the successful access of this module.
