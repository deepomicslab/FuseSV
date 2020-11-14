package FuseSV::Virus_Integ::LocalHaplotypeDetect::LocalMap_OB;

use strict;
use warnings;
use List::Util qw/ min max sum shuffle first /;
use Digest::MD5 qw/ md5_hex /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Random qw/ GetRandBool /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'LocalMap_OB';
#----- version --------
$VERSION = "0.18";
$DATE = '2018-06-08';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						set_randProb
						is_valid
						get_free_UnitCycle_count
						get_structureDepth
						get_NO
						get_randProbAref
						get_fvgmInfoStr
						get_segStrMd5
						get_nosSegStrMd5
						get_struHeight
						get_length
						get_PWanchorSimuLink
						hasSegStr
						load_PWanchorSimuLink
						merge_UnitCycle
						find_basicUnitCycle
						reverse_ucyc
						reverse_innerSegments
						anchorSegToBasicUCYC
						retionalizeAnchorSeg
						insertNonbsUcycToBasicUcyc
						absorbRemainUCYCtoBasicUCYC
						release_ucyc
						simplify
						add_fvgm
						calculate_md5
						calculate_struHeight
						calculate_SimuPWanchorLink
						display
						show_segStructure
						show_noOritSortSegStr
					 /;

#--- structure of object
# localmap_OB -> NO = $NO
# localmap_OB -> valid = 1/0, initilize as 1
# localmap_OB -> groupOB = $groupOB
# localmap_OB -> allele_OB = $allele_OB
# localmap_OB -> randProbAref = $randProbAref
# localmap_OB -> UnitCycle -> $UcycID = { UcycOB=>$ucycOB, count=>$count, UcycToAvoidMerge=>{} }
## check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB module
# localmap_OB -> basicUcycOB = $UnitCycle_OB
## FVGM used in this localmap
# localmap_OB -> FVGM = { $FVGM_UcycID => $usedCount, ... };
## md5 of segments strucure
# localmap_OB -> segStrMd5 = $segStrMd5
## md5 of sorted-segments after oritation-ignorance
## nos: Non-Oritation Sorted
# localmap_OB -> nosSegStrMd5 = $nosSegStrMd5
## the height of the structure got from dispay func
# localmap_OB -> struHeight = $struHeight
## simulate links of pair-wise anchors
# localmap_OB -> PWanchorSimuLink = { $pwanchorID -> $simuLink, .. }

# deprecated!!!
## Pearson correlation coefficient based on distance-matrix
# localmap_OB -> score = $score
## prob distance of pair-wise anchors
## used in func 'calculate_prob_distance'
# localmap_OB -> PWanchorProbDist = { $pwanchorID -> $prob_dist, .. }

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $NO = $parm{NO} || 1;
	my $allele_OB = $parm{allele_OB};
	my $groupOB = $parm{groupOB};

	my $localmap_OB = {};

	$localmap_OB->{NO} = $NO;
	$localmap_OB->{valid} = 1; # initilize as valid
	$localmap_OB->{groupOB} = $groupOB;
	$localmap_OB->{allele_OB} = $allele_OB;
	$localmap_OB->{randProbAref} = [];
	$localmap_OB->{segStrMd5} = 'N/A';
	$localmap_OB->{nosSegStrMd5} = 'N/A';
	$localmap_OB->{struHeight} = 1;
	# load UCYC
	for my $UcycID ( @{$allele_OB->get_UcycIDlistAref} ){
		my ($UcycOB, $UcycCount) = $allele_OB->get_UcycOB( UcycID => $UcycID );
		my $cpUcycOB = $UcycOB->get_new_cpUcycOB;
		my $OrigID = $cpUcycOB->get_origID; # should be same as $UcycID
		$localmap_OB->{UnitCycle}->{$OrigID} = {
													UcycOB => $cpUcycOB,
													count => $UcycCount,
													UcycToAvoidMerge => { $OrigID => 1 } # never itself, but obsoleted
												};
	}

	bless($localmap_OB);
	return $localmap_OB;
}

#--- set rand probability
sub set_randProb{
	my $localmap_OB = shift;
	my %parm = @_;
	my $randProbAref = $parm{randProbAref} || [0]; # random mode

	# srand(time) if( $randProbAref->[0] != 0 );
	$localmap_OB->{randProbAref} = $randProbAref;
}

#--- return valid status
sub is_valid{
	my $localmap_OB = shift;
	return $localmap_OB->{valid};
}

#--- return count of remaining free UCYC(s)
sub get_free_UnitCycle_count{
	my $localmap_OB = shift;
	return scalar(keys %{$localmap_OB->{UnitCycle}});
}

#--- return structure depth
sub get_structureDepth{
	my $localmap_OB = shift;

	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	return $basicUcycOB->get_structureDepth;
}

#--- return the NO of this local map
sub get_NO{
	my $localmap_OB = shift;
	return $localmap_OB->{NO};
}

#--- return randProbAref of this local map
sub get_randProbAref{
	my $localmap_OB = shift;
	return $localmap_OB->{randProbAref};
}

#--- return fvgm info-str used in this local map
sub get_fvgmInfoStr{
	my $localmap_OB = shift;

	my @infoStr = map { $_ .':'. $localmap_OB->{FVGM}->{$_} } sort keys %{$localmap_OB->{FVGM}};
	return join(',',@infoStr) || 'N/A';
}

#--- return md5 of segments structure
sub get_segStrMd5{
	my $localmap_OB = shift;
	return $localmap_OB->{segStrMd5};
}

#--- return md5 of sorted-segments after oritation-ignorance
sub get_nosSegStrMd5{
	my $localmap_OB = shift;
	return $localmap_OB->{nosSegStrMd5};
}

#--- return the height of the structure got from dispay func
sub get_struHeight{
	my $localmap_OB = shift;
	return $localmap_OB->{struHeight};
}

#--- return local map length
sub get_length{
	my $localmap_OB = shift;

	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	return $basicUcycOB->get_length;
}

#--- return ref of hash contains simu.link of pw-anchor
sub get_PWanchorSimuLink{
	my $localmap_OB = shift;
	return $localmap_OB->{PWanchorSimuLink};
}

#--- test whether have segment string provided
sub hasSegStr{
	my $localmap_OB = shift;
	my %parm = @_;
	my $regex = $parm{regex};


	if( !defined $regex ){
		my $func_id = (caller(0))[3];
		warn_and_exit "<ERROR>\tplease set 'regex' for $func_id.\n";
	}

	my $segStr = $localmap_OB->show_segStructure;
	if( $segStr =~ /$regex/ ){
		return 1;
	}
	else{
		return 0;
	}
}

#--- load link count of given pwanchor
sub load_PWanchorSimuLink{
	my $localmap_OB = shift;
	my %parm = @_;
	my $pwanchorID = $parm{pwanchorID};
	my $pwanchorLinkC = $parm{pwanchorLinkC};
	$localmap_OB->{PWanchorSimuLink}->{$pwanchorID} = $pwanchorLinkC;
}

#--- merge Unit Cycle
sub merge_UnitCycle{
	my $localmap_OB = shift;
	my %parm = @_;
	my $source = $parm{source};
	my $sink = $parm{sink};
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	my $useFVGM  = $parm{useFVGM} || 0; # random mode

	# record randProbAref
	$localmap_OB->set_randProb( randProbAref => $randProbAref );

	# get the basic UCYC of the local haplotype
	$localmap_OB->find_basicUnitCycle( source => $source, sink => $sink );
	my $basicUcycOB = $localmap_OB->{basicUcycOB};

	# reverse whole-ucyc in random mode
	$localmap_OB->reverse_ucyc;

	# try to change all non-basis unit-cycle to has anchorSeg along basic unit-cycle (5p->3p)
	$localmap_OB->anchorSegToBasicUCYC;

	# make sure all UCYC except basic UCYC have their own anchorSeg in other UCYC(s)
	$localmap_OB->retionalizeAnchorSeg;

	# reverse innersegments in random mode
	$localmap_OB->reverse_innerSegments;

	# try to insert nonbsUCYC into basic UCYC
	$localmap_OB->insertNonbsUcycToBasicUcyc;

	# if still remains UCYC not inserted,
	# just add the count of itself that already in basic UCYC
	$localmap_OB->absorbRemainUCYCtoBasicUCYC;

	# add free virus genome randomly
	$localmap_OB->add_fvgm( useFVGM => $useFVGM );

	# make local map simple for clear display
	$localmap_OB->simplify;

	# calculate md5 of two types of segments structure
	$localmap_OB->calculate_md5;

	# calculate the height of the structure got from dispay func
	$localmap_OB->calculate_struHeight;
}

#--- find basic Unit Cycle
## 1) start from source, end at sink, 2) the longest, and 3) has only one copy.
sub find_basicUnitCycle{
	my $localmap_OB = shift;
	my %parm = @_;
	my $source = $parm{source};
	my $sink = $parm{sink};

	my $UCYC_Href= $localmap_OB->{UnitCycle};

	# find unit-cycle starts from source, ends at sink.
	my %cand = map { ($_, 1) }
			   grep    $UCYC_Href->{$_}->{UcycOB}->get_anchorSegID eq $source
					&& $UCYC_Href->{$_}->{UcycOB}->get_lastSegID   eq $sink
					&& $UCYC_Href->{$_}->{count} == 1, # copy of source->sink should be one in SINGLE COPY of local map
			   keys %$UCYC_Href;

	# sort by length
	my @cand = sort {    $UCYC_Href->{$b}->{UcycOB}->get_length
					 <=> $UCYC_Href->{$a}->{UcycOB}->get_length }
			   keys %cand;

	# if empty, alert
	if( scalar(@cand) == 0 ){
		warn_and_exit "<ERROR>\tCannot find source($source)-to-sink($sink) unit cyle in local map\n";
	}

	# record
	$localmap_OB->{basicUcycOB} = $UCYC_Href->{$cand[0]}->{UcycOB}->get_new_cpUcycOB;

	# inform
	stout_and_sterr "[INFO]\tFind UCYC ".$localmap_OB->{basicUcycOB}->get_origID." as basic_ucyc.\n";
}

#--- reverse whole-ucyc in random mode
sub reverse_ucyc{
	my $localmap_OB = shift;

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	return unless( $randProbAref->[0] );

	for my $nonbsUcycID ( shuffle grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href ){
		# replace with the complemented and reversed version
		if( GetRandBool(verbose=>0,Aref=>$randProbAref) ){
			my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
			$UCYC_Href->{$nonbsUcycID}->{UcycOB} = $nonbsUcycOB->get_new_crUcycOB;
			# inform
			stout_and_sterr "[INFO]\treverse ".$nonbsUcycOB->get_origID."\n";
		}
	}
}

#--- reverse innersegments of ucyc in random mode
sub reverse_innerSegments{
	my $localmap_OB = shift;

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	return unless( $randProbAref->[0] );

	for my $nonbsUcycID ( shuffle grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href ){
		# reverse valid inner-segments
		if( GetRandBool(verbose=>0,Aref=>$randProbAref) ){
			my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
			$nonbsUcycOB->reverseInnerSegs( randProbAref => $randProbAref );
		}
	}
}

#--- change non-basis unit-cycle to has anchorSeg along basic unit-cycle
# simple mode: 5p->3p; random mode: shuffle
sub anchorSegToBasicUCYC{
	my $localmap_OB = shift;

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	my @segCombHref = @{$basicUcycOB->{segComb}};
	   @segCombHref = shuffle @segCombHref if( $randProbAref->[0] );
	my @nonbsUcycID = grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href;
	   @nonbsUcycID = shuffle @nonbsUcycID if( $randProbAref->[0] );
	my %doneNobUcyc;
	for my $segCombHref ( @segCombHref ){
		next if( $segCombHref->{OBtype} ne 'SEG' ); # must be Segment
		next if( $segCombHref->{OB}->get_SpecTag ne 'H' ); # must be host genome
		next if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
		my $candAnchorSegID = $segCombHref->{OB}->get_ID;
		for my $nonbsUcycID ( grep !exists $doneNobUcyc{$_}, @nonbsUcycID ){
			my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
			next if( !$randProbAref->[0] && $nonbsUcycOB->is_anchorSegForBasis );
			$nonbsUcycOB->change_to_aimAnchorSeg( aimAnchorSeg => $candAnchorSegID, aimAnchorSegIsForBasis => 1, informPref => 'ToBasicUCYC', randProbAref => $randProbAref );
			# if successfully changes
			if( $nonbsUcycOB->get_anchorSegID eq $candAnchorSegID ){
				$doneNobUcyc{$nonbsUcycID} = 1; # marker for next skip
			}
			last if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
		}
	}
}

#--- make sure all UCYC except basic UCYC have their own anchorSeg in other UCYC(s)
sub retionalizeAnchorSeg{
	my $localmap_OB = shift;

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	while(1){
		my @nonbsUcycID = grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href;
		my $all_retional_bool = 1;
		for my $nonbsUcycID (@nonbsUcycID){
			my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
			my $nonbsUcycAnchorSeg = $nonbsUcycOB->get_anchorSegID;
			my @checkUcycID = grep $_ ne $nonbsUcycID, keys %$UCYC_Href;
			my $this_retional_bool = 0;
			for my $checkUcycID ( @checkUcycID ){
				my $checkUcycOB = $UCYC_Href->{$checkUcycID}->{UcycOB};
				my $locNO_Aref = $checkUcycOB->has_this_segID( segID => $nonbsUcycAnchorSeg, seekMode => 'SEG' );
				if( $locNO_Aref->[0] != 0 ){ # found
					$this_retional_bool = 1;
					last; # just found once is ok.
				}
			}
			# if not retional, change nonbsUcyc's AnchorSeg
			if(    $this_retional_bool == 0
				|| $randProbAref->[0] # must have at least one time 'anchor-change' in random mode
			){
				# update global retional bool
				$all_retional_bool &&= $this_retional_bool;
				# find one candidate anchor
				my $aimCombEleHref;
				for my $SpecTag ( 'H', 'V' ){ # host has priority
					$aimCombEleHref = first { defined $_ }
									  shuffle
									  grep    $_->{OB}->get_SpecTag eq $SpecTag
										   && $_->{OB}->get_ID ne $nonbsUcycAnchorSeg,
									  @{$nonbsUcycOB->{segComb}};
					last if( defined $aimCombEleHref );
				}
				# if not found! never happen!
				if( !defined $aimCombEleHref ){
					warn_and_exit "<ERROR>\t{RetionalAnchorSeg}, UCYC $nonbsUcycID: cannot find any proper Seg as aimAnchorSeg!\n";
				}
				my $aimAnchorSeg = $aimCombEleHref->{OB}->get_ID; # previous: $nonbsUcycOB->get_lastSegID
				$nonbsUcycOB->change_to_aimAnchorSeg( aimAnchorSeg => $aimAnchorSeg, informPref => 'Retionalize', randProbAref => $randProbAref );
				# inform
				stout_and_sterr "[INFO]\t{RetionalAnchorSeg}, UCYC $nonbsUcycID: change from $nonbsUcycAnchorSeg to $aimAnchorSeg.\n";
				# stop to check more nonbsUcyc, need to restart all loop
				last unless( $all_retional_bool );
			}
		}
		# if all nonbsUcyc's anchorSeg are retional, finish!
		last if( $all_retional_bool );
	}
}

#--- insert non-basis UCYCs to Basic UCYC
sub insertNonbsUcycToBasicUcyc{
	my $localmap_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!
	my $given_nonbsUcycID = $parm{given_nonbsUcycID}; # just deal with this given_nonbsUcycID if set

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	# basic UCYC tries to absorb non-basis UCYC 
	my $try_times = scalar(keys %$UCYC_Href) ** 2;
	my $while_time = 0;
	while( $while_time++ < $try_times ){
		# prepare list of non-basis UCYC
		my @nonbsUcycID = sort { $UCYC_Href->{$a}->{UcycOB}->get_length * $UCYC_Href->{$a}->{count}
								 <=>
								 $UCYC_Href->{$b}->{UcycOB}->get_length * $UCYC_Href->{$b}->{count}
								} grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href;
		# only keep the given nonbsUcycID if set
		if( defined $given_nonbsUcycID ){
			@nonbsUcycID = grep { $_ eq $given_nonbsUcycID } @nonbsUcycID;
		}
		# stop if no more nonbsUCYC to absorb
		last if( scalar(@nonbsUcycID) == 0 );
		# random mode?
		if( $randProbAref->[0] ){
			@nonbsUcycID = shuffle( @nonbsUcycID );
			# inform
			stout_and_sterr "[INFO]\ttry_NO.$while_time, Shuffle nonbsUcycID before insertion to basic_ucyc " . $basicUcycOB->get_origID . "\n";
		}
		# try one by one
		for my $nonbsUcycID ( @nonbsUcycID ){
			my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
			my $nonbsUcycOrigCount = $UCYC_Href->{$nonbsUcycID}->{count};
			my $nonbsUcycAvoidUcycHref = $UCYC_Href->{$nonbsUcycID}->{UcycToAvoidMerge};
			# let the basicUcyc to absorb this nonbsUcyc
			my $absorbCount = $basicUcycOB->absorbUcyc(
														thisUcycOrigCount => 1,
														absbUcycOB => $nonbsUcycOB,
														absbUcycOrigCount => $nonbsUcycOrigCount,
														avoidUcycHref => $nonbsUcycAvoidUcycHref,
														randProbAref => $randProbAref,
														sayID => $sayID
													  );
			# if absorb
			if( $absorbCount != 0 ){
				# release these absorbed
				$localmap_OB->release_ucyc( ucycOB => $nonbsUcycOB, delCount => $absorbCount );
				# start new loop after each absorb process
				last;
			}
		}
	}
}

#--- after insertion, may have remained UCYCs
## just add the remained count on itself found in Basic UCYC
sub absorbRemainUCYCtoBasicUCYC{
	my $localmap_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $UCYC_Href = $localmap_OB->{UnitCycle};

	while(1){
		my @nonbsUcycID = grep $_ ne $basicUcycOB->get_origID, keys %$UCYC_Href;
		last if( scalar(@nonbsUcycID) == 0 );
		# shuffle in random mode
		@nonbsUcycID = shuffle( @nonbsUcycID ) if( $randProbAref->[0] );
		# deal with the first one
		my $nonbsUcycID = $nonbsUcycID[0];
		my $nonbsUcycOB = $UCYC_Href->{$nonbsUcycID}->{UcycOB};
		my $nonbsUcycCount = $UCYC_Href->{$nonbsUcycID}->{count};
		# try to absorb
		my $asInfoAref = $basicUcycOB->absorbSameUcycItHas( testUcycOB => $nonbsUcycOB, testUcycCount => $nonbsUcycCount, 
															randProbAref => $randProbAref, sayID => $sayID );
		my ($addBool, $addCount) = @{$asInfoAref};
		if( $addBool ){ # successful, and release the used counts
			$localmap_OB->release_ucyc( ucycOB => $nonbsUcycOB, delCount => $addCount );
		}
		else{
			# try to change this nonbsUcyc's anchorSeg
			for my $try_time ( 1 .. scalar(@{$nonbsUcycOB->{segComb}}) ){
				my $idx = first {   $nonbsUcycOB->{segComb}->[$_]->{OBtype} eq 'SEG'
								 && $nonbsUcycOB->{segComb}->[$_]->{OB}->get_ID ne $nonbsUcycOB->get_anchorSegID
								} ( 1 .. scalar(@{$nonbsUcycOB->{segComb}})-1 );
				my $aimAnchorSeg = $nonbsUcycOB->{segComb}->[$idx]->{OB}->get_ID;
				# change this nonbsUcyc's anchorSeg
				$nonbsUcycOB->change_to_aimAnchorSeg( aimAnchorSeg => $aimAnchorSeg, informPref => 'absbRemUCYC', randProbAref => $randProbAref );
				# try to insert again
				$localmap_OB->insertNonbsUcycToBasicUcyc( given_nonbsUcycID => $nonbsUcycID, sayID => $sayID );
				# stop when successfully all counts inserted
				last if( !exists $UCYC_Href->{$nonbsUcycID} );
			}
			# this nonbsUcyc still has some count(s) left!
			if( exists $UCYC_Href->{$nonbsUcycID} ){
				warn "<WARN>\tIn local map NO.$localmap_OB->{NO} of allele NO.".$localmap_OB->{allele_OB}->get_NO."\n";
				warn "\tThese unit-cycle still remain, cannot merge with basic unit-cycle ($basicUcycOB->{$sayID}).\n"
					."\t@nonbsUcycID\n"
					."\tFailed at $nonbsUcycID (remain count=$UCYC_Href->{$nonbsUcycID}->{count}):\n"
					."\t".$nonbsUcycOB->show_segStructure."\n";
				warn "Current structure of this local map is:\n";
				warn $localmap_OB->display( prefix => "LOCALMAP:\t" );
				# not valid
				$localmap_OB->{valid} = 0;
				last;
			}
		}
	}
}

#--- release unit-cycle as its centain count has been used
sub release_ucyc{
	my $localmap_OB = shift;
	my %parm = @_;
	my $ucycOB = $parm{ucycOB};
	my $delCount = $parm{delCount};

	my $UCYC_Href = $localmap_OB->{UnitCycle};
	my $ucycOrigID = $ucycOB->get_origID; # $UCYC_Href contains initial ID, i.e., the origID

	# if given a count to subtract
	if( defined $delCount && $delCount >= 0 ){
		# subtract these inserted
		$UCYC_Href->{$ucycOrigID}->{count} -= $delCount;
	}
	# delete if none left
	if( $UCYC_Href->{$ucycOrigID}->{count} == 0 ){
		delete $UCYC_Href->{$ucycOrigID};
	}
}

#--- make local map simple for clear display
sub simplify{
	my $localmap_OB = shift;

	my $basicUcycOB = $localmap_OB->{basicUcycOB};

	# if possible, merge inner UCYC with flanking Element's Str (same as inner UCYC)
	$basicUcycOB->merge_innerUcycWithFlankEleStr;

	# if possible, let single copy inner UCYC to throw sequential host-Seg from tail
	$basicUcycOB->throwTailSegFromSingleInnerUCYC;

	# if possible, let single copy inner UCYC to absorb flanking sequential Seg
	$basicUcycOB->absorbSegToSingleInnerUCYC;
}

#--- add free vrius genome into local map in random mode
sub add_fvgm{
	my $localmap_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!
	my $useFVGM = $parm{useFVGM} || 0; # random mode

	my $randProbAref = $localmap_OB->{randProbAref}; # random mode
	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	my $allele_OB = $localmap_OB->{allele_OB};

	# add FVGM only in random mode
	return unless( $randProbAref->[0] );

	unless( $useFVGM ){
		# inform
		stout_and_sterr "<WARN>\tMerging FVGM unit-cycle to local map is disabled.\n";
		return;
	}

	my $avoidUcycHref = {};
	for my $FVGM_UcycID ( @{ $allele_OB->get_UcycIDlistAref( key => 'FVGM_UCYC' ) } ){
		my ($FVGM_UcycOB, $FVGM_UcycCount) = $allele_OB->get_UcycOB( key => 'FVGM_UCYC', UcycID => $FVGM_UcycID );
		# to avoid itself
		# note: $FVGM_UcycID is origID
		$avoidUcycHref->{ $FVGM_UcycID } = 1;
		# let the basicUcyc to absorb this FVGM_UcycOB
		while( $FVGM_UcycCount > 0 ){
			# random stop
			last if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
			# try to merge
			my $absorbCount = $basicUcycOB->absorbUcyc(
														thisUcycOrigCount => 1,
														absbUcycOB => $FVGM_UcycOB,
														absbUcycOrigCount => $FVGM_UcycCount,
														avoidUcycHref => $avoidUcycHref, # avoid this FVGM and previous FVGMs
														randProbAref => $randProbAref,
														sayID => $sayID
													  );
			# merge successfully?
			if( $absorbCount > 0 ){
				# inform
				stout_and_sterr "[INFO]\tMerge $absorbCount count of FVGM unit-cycle $FVGM_UcycOB->{$sayID} to local map.\n";
				# update
				$FVGM_UcycCount -= $absorbCount;
				# record
				$localmap_OB->{FVGM}->{$FVGM_UcycID} += $absorbCount;
			}
			else{
				last;
			}
		}
	}
}

#--- calculate md5 of two types of segments structure
sub calculate_md5{
	my $localmap_OB = shift;

	$localmap_OB->{segStrMd5} = 'N/A';
	$localmap_OB->{nosSegStrMd5} = 'N/A';

	# md5 of segments structure
	my $segStr = $localmap_OB->show_segStructure;
	$localmap_OB->{segStrMd5} = md5_hex( $segStr );
	# md5 of sorted-segments after oritation-ignorance
	my $segStr_NoOritSort = $localmap_OB->show_noOritSortSegStr;
	$localmap_OB->{nosSegStrMd5} = md5_hex( $segStr_NoOritSort );
}

#--- calculate the height of the structure got from dispay func
sub calculate_struHeight{
	my $localmap_OB = shift;
	$localmap_OB->{struHeight} = scalar( split /\n+/, $localmap_OB->display );
}

#--- simulate sequencing to calculate pair-wise links
sub calculate_SimuPWanchorLink{
	my $localmap_OB = shift;

	my $groupOB = $localmap_OB->{groupOB};
	my $allele_OB = $localmap_OB->{allele_OB};
	# load simulated sequencing reads done before
	my $SimuSeqRead_Href = $allele_OB->get_LM_SimuSeqReads_Hash;
	my $LRbiEndSeq_BOOL = $allele_OB->is_SimuLRbiEndSeq;
	my $SimuSeqRead_rLen = $allele_OB->get_SimuSeqRlen;
	# segment structure
	my @segStr = map{
					/^([HV]\d+)([\+\-])$/;
					({info=>$_, segID=>$1, segOrit=>$2});
				 } split /;/, $localmap_OB->show_segStructure;
	my $seg_LastIdx = $#segStr;
	# calculate each pair-wise anchor linkage
	for my $pwanchorOB ( values %{$groupOB->get_PWanchorOB} ){
		my ($anchorLoc_Aref, $valid) = $pwanchorOB->find_anchor_idx( segStr_Aref => \@segStr, seg_LastIdx => $seg_LastIdx );
		# valid? i.e., found both anchors
		next unless( $valid );
		# initialize link value
		my $pwanchorID = $pwanchorOB->get_ID;
		$localmap_OB->{PWanchorSimuLink}->{$pwanchorID} = 0;
		my $anchorLoc_count = scalar(@$anchorLoc_Aref);
		my $anchorFirst_loc = $anchorLoc_Aref->[0]->{accuLen};
		my $anchorLast_loc = $anchorLoc_Aref->[-1]->{accuLen};
		# check each simulated sequencing HMW/read
		# for my $rNO (sort {$a<=>$b} keys %$SimuSeqRead_Href){
		for my $rHref (sort {$a->{stpos} <=> $b->{stpos}} values %$SimuSeqRead_Href){
			my $stpos = $rHref->{stpos};
			last if( $stpos > $anchorLast_loc );
			my $edpos = $rHref->{edpos};
			next if( $edpos < $anchorFirst_loc );
			# deal with each pair-wise anchor a-b or b-a
			for my $i ( 0 .. $anchorLoc_count-2 ){
				my $i_anchor_accuLen = $anchorLoc_Aref->[$i]->{accuLen};
				next if( $i_anchor_accuLen < $stpos );
				last if( $i_anchor_accuLen > $edpos );
				next if( $LRbiEndSeq_BOOL && abs($i_anchor_accuLen-$stpos) > $SimuSeqRead_rLen );
				for my $j ( $i .. $anchorLoc_count-1 ){
					next if( $anchorLoc_Aref->[$j]->{side} eq $anchorLoc_Aref->[$i]->{side} );
					my $j_anchor_accuLen = $anchorLoc_Aref->[$j]->{accuLen};
					last if( $j_anchor_accuLen > $edpos );
					next if( $LRbiEndSeq_BOOL && abs($j_anchor_accuLen-$edpos) > $SimuSeqRead_rLen );
					# record this link
					$localmap_OB->{PWanchorSimuLink}->{$pwanchorID} += $rHref->{covR};
				}
			}
		}
	}

	# inform
	my $group_ID = $groupOB->get_ID;
	my $allele_NO = $localmap_OB->{allele_OB}->get_NO;
	my $localmap_NO = $localmap_OB->get_NO;
	stout_and_sterr "[INFO]\tGroup $group_ID Allele NO.$allele_NO LocalMap NO.$localmap_NO, get simulated PWanchor linkages OK.\n";
}

#--- calculate probabilistic distance
# !!! deprecated
sub calculate_prob_distance{
	# my $localmap_OB = shift;

	# # segment structure
	# my $segStr = $localmap_OB->show_segStructure;
	# my @segStr = split /;/, $segStr;
	# # process each pair-wise anchor
	# my $groupOB = $localmap_OB->{groupOB};
	# for my $pwanchorOB ( values %{$groupOB->get_PWanchorOB} ){
	# 	my $anchorIdx_Aref = $pwanchorOB->find_anchor_idx( segStr => $segStr );
	# 	# check both anchors found
	# 	my $pwanchor_a_count = grep $_->{side} eq 'a', @$anchorIdx_Aref;
	# 	my $pwanchor_b_count = grep $_->{side} eq 'b', @$anchorIdx_Aref;
	# 	next if( $pwanchor_a_count == 0 || $pwanchor_b_count == 0 );
	# 	# deal with each pair-wise anchor a-b or b-a
	# 	my $anchorIdx_count = scalar(@$anchorIdx_Aref);
	# 	for my $i ( 0 .. $anchorIdx_count-2 ){
	# 		my $this_anchor_Href = $anchorIdx_Aref->[$i];
	# 		my $next_diff_side_i = first { $anchorIdx_count->[$_]->{side} ne $this_anchor_Href->{side} } ( $i .. $anchorIdx_count-1 );
	# 		last if( !defined $next_diff_side_i );
	# 		my $diff_anchor_Href = $anchorIdx_Aref->[$next_diff_side_i];
	# 		# get basic distance
	# 		my $basic_dist = $pwanchorOB->get_anchor_distance( anchor_a_Href => $this_anchor_Href, anchor_b_Href => $diff_anchor_Href, segStr => $segStr );
	# 		# get bilater extend distance
	# 		my $up_anchor_Href = ( $i != 0 ? $anchorIdx_Aref->[$i-1] : {idx=>-1, type=>'JUNC'} ); # pretend to be JUNC from the start
	# 		my $up_extend_dist = $pwanchorOB->get_anchor_distance( anchor_a_Href => $up_anchor_Href, anchor_b_Href => $this_anchor_Href, segStr => $segStr );
	# 		my $dw_anchor_Href = ( $next_diff_side_i != $anchorIdx_count-1 ? $anchorIdx_Aref->[$next_diff_side_i+1] : {idx=>$#segStr, type=>'JUNC'} ); # pretend to be JUNC at the END
	# 		my $dw_extend_dist = $pwanchorOB->get_anchor_distance( anchor_a_Href => $diff_anchor_Href, anchor_b_Href => $dw_anchor_Href, segStr => $segStr );
	# 		# get ceiling distance
	# 		my $ceiling_dist = $basic_dist + max($up_extend_dist, $dw_extend_dist);
	# 		# get prob.dist from basix distance to ceiling distance
	# 		my $prob_dist = $groupOB->get_prob_distance( from => $basic_dist, to => $ceiling_dist );
	# 		# record prob.dist
	# 		$localmap_OB->{PWanchorProbDist}->{$pwanchorOB->get_ID} += $prob_dist;
	# 	}
	# }
}

#--- display the local map
# i.e., the basic unit-cycle
sub display{
	my $localmap_OB = shift;
	my %parm = @_;
	my $prefix = $parm{prefix} || '';
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	return $basicUcycOB->display( prefix=>$prefix, count=>1, sayID=>$sayID );
}

#--- show segment structure of this localmap
sub show_segStructure{
	my $localmap_OB = shift;

	my $basicUcycOB = $localmap_OB->{basicUcycOB};
	return $basicUcycOB->show_segStructure;
}

#--- show sorted-segments after oritation-ignorance
sub show_noOritSortSegStr{
	my $localmap_OB = shift;

	( my $segStr_NoOritSort = join(';', sort split ';', $localmap_OB->show_segStructure) ) =~ s/[\+\-]//g;
	return $segStr_NoOritSort;
}

1; ## tell the perl script the successful access of this module.
