package FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB;

use strict;
use warnings;
use List::Util qw/ min max sum shuffle /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr getStrUnitRepeatTime GetRandBool /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'UnitCycle_OB';
#----- version --------
$VERSION = "0.11";
$DATE = '2018-05-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						update_info
						get_ID
						get_origID
						get_length
						is_FVGM
						get_fvgmID
						get_anchorSegID
						get_lastSegID
						get_structureDepth
						change_to_aimAnchorSeg
						is_anchorSegForBasis
						breakSingleEleUCYC
						show_segStructure
						is_identical_UcycOB
						has_this_segID
						has_innerUCYC
						display
						insertUcycToReplaceGivenSeg
						get_new_cpUcycOB
						get_new_crUcycOB
						reverseInnerSegs
						absorbUcyc
						allot_innerUcyc
						absorbSameUcycItHas
						merge_innerUcycWithFlankEleStr
						throwTailSegFromSingleInnerUCYC
						absorbSegToSingleInnerUCYC
					 /;

#--- structure of object
# ucyc_OB -> ID = $ID
# ucyc_OB -> segm -> $segID -> $orit = $count
# ucyc_OB -> length = $length
# ucyc_OB -> segComb -> $A_ref = [ OBtype=>'SEG/UCYC', OB=>$segOB/$ucycOB, orit=>$orit, count=>$count ];
# ucyc_OB -> fvgm = [ $is_fvgm_Bool, $fvgmID ]
# ucyc_OB -> anchorSegID = $anchorSegID
# ucyc_OB -> lastSegID = $lastSegID
# ucyc_OB -> structureDepth = $structureDepth
# ucyc_OB -> anchorSegForBasis = 0/1

#--- shortname postfix of transformation of UCYC object
# '.cr' means (C)omplemented-(R)eversed
# '.cpN' means this is the NO.N time to make copy
# '.iN[xxx]' means insert UCYC xxx at loc NO.N
# '.a2' means the 2nd alloted group
# '.asNXXX' means (A)dd (S)ame UCYC (ID:XXX) at NO.N local
# '.tN' means (T)hrow N seg(s)
# '.b5' means absorb seg from 5-prime
# '.b3' means absorb seg from 3-prime

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $ID = $parm{ID}; # ID after tranformation
	my $origID = $parm{origID} || $ID; # original ID from initial input
	my $segComb = $parm{segComb};
	my $fvgmID = $parm{fvgmID};
	my $segHref = $parm{segHref};
	my $strOpmz = $parm{strOpmz} || 0;

	my $ucyc_OB = {};

	$ucyc_OB->{ID} = $ID;
	$ucyc_OB->{origID} = $origID;
	$ucyc_OB->{cpNO} = 0; # to control copy time at each transformation
	$ucyc_OB->{fvgm} = (defined $fvgmID) ? [ 1, $fvgmID ] : [ 0, 'N/A' ];

	if( defined $segComb ){
		# split
		my @segComb = split /[,;]/, $segComb;
		pop @segComb if( $segComb[0] eq $segComb[-1] ); # last one is duplicated
		# optimize unit-cycle structure
		## set anchorSeg as host segment
		## set the cycle-bkage at the real breakage
		# only optimize structure with host segments
		my $do_times = ( $strOpmz && $segComb =~ /\bH\d+\b/ ) ? 0 : -1;
		# process of optimize structure
		while(++$do_times){
			my ($firstSegType, $firstSegNum) = ( $segComb[0]  =~ /^(\D+)(\d+)/ );
			my ($lastSegType,  $lastSegNum)  = ( $segComb[-1] =~ /^(\D+)(\d+)/ );
			if(   $firstSegType ne 'H' # anchorSeg is not host segment
				  || (    $lastSegType eq 'H'
				       && $firstSegNum > $lastSegNum # adjust the 5-prime host segment at 5-prime of unit cycle
				     )
			){
				unshift @segComb, pop @segComb;
			}
			else{
				last;
			}
			# forever LOOP: whose unit-cycle has 10w segments? ^_^
			if( $do_times == 100000 ){
				warn_and_exit "<ERROR>\tCannot finish optimizing structure of unit-cycle $ID.\n";
			}
		}
		# records
		for my $i (0 .. $#segComb){
			my ($segID, $orit) = ($segComb[$i] =~ /^(.+)([\+\-])$/);
			if( !exists $segHref->{$segID} ){
				warn_and_exit "<ERROR>\tCannot recognize segment $segID from unit-cycle $ID\n";
			}
			$segComb[$i] = { OBtype=>'SEG', OB=>$segHref->{$segID}, orit=>$orit, count=>1 };
		}
		$ucyc_OB->{segComb} = \@segComb;
		$ucyc_OB->{anchorSegForBasis} = 0;
	}

	bless($ucyc_OB);
	return $ucyc_OB;
}

#--- load in the sub region details
sub update_info{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID';

	# check single segComb
	$ucyc_OB->breakSingleEleUCYC;

	# anchor Segment
	$ucyc_OB->{anchorSegID} =   $ucyc_OB->{segComb}->[0]->{OBtype} eq 'SEG'
							  ? $ucyc_OB->{segComb}->[0]->{OB}->get_ID
							  : $ucyc_OB->{segComb}->[0]->{OB}->get_anchorSegID;
	# last Segment
	$ucyc_OB->{lastSegID}   =   $ucyc_OB->{segComb}->[-1]->{OBtype} eq 'SEG'
							  ? $ucyc_OB->{segComb}->[-1]->{OB}->get_ID
							  : $ucyc_OB->{segComb}->[-1]->{OB}->get_lastSegID;

	# reset attributes
	$ucyc_OB->{segm} = {};
	$ucyc_OB->{length} = 0;
	$ucyc_OB->{structureDepth} = 1;

	my $innerUcyc_maxDepth = 0;
	# count segments, sum length, and innerUcyc max structure Depth
	for my $CombEleHref ( @{$ucyc_OB->{segComb}} ){
		my $CombEleOB      = $CombEleHref->{OB};
		my $CombEleOBtype  = $CombEleHref->{OBtype};
		my $CombEleOBorit  = $CombEleHref->{orit};
		my $CombEleOBcount = $CombEleHref->{count};
		## it's a segment
		if( $CombEleOBtype eq 'SEG' ){
			my $segID = $CombEleOB->get_ID;
			$ucyc_OB->{segm}->{$segID}->{$CombEleOBorit} += $CombEleOBcount;
		}
		## unit-cycle
		elsif( $CombEleOBtype eq 'UCYC' ){
			# check, all inner UCYC are required to have '+' orit
			if( $CombEleOBorit ne '+' ){
				warn_and_exit "<ERROR>\tInner unit-cycle ($CombEleOB->{$sayID}) has '-' orit in upper unit-cycle ($ucyc_OB->{$sayID})!\n";
			}
			# so, no need to convert the orit
			for my $innerSegID (keys %{$CombEleOB->{segm}}){
				for my $innerSegOrit (keys %{$CombEleOB->{segm}->{$innerSegID}}){
					my $innerSegCount = $CombEleOB->{segm}->{$innerSegID}->{$innerSegOrit};
					$ucyc_OB->{segm}->{$innerSegID}->{$innerSegOrit} += $CombEleOBcount * $innerSegCount;
				}
			}
			# innerUcyc structure depth
			$innerUcyc_maxDepth = max( $innerUcyc_maxDepth, $CombEleOB->get_structureDepth );
		}
		# add length
		$ucyc_OB->{length} += $CombEleOB->get_length * $CombEleOBcount;
	}
	# update structure depth
	$ucyc_OB->{structureDepth} = $innerUcyc_maxDepth + 1;
}

#--- return unit-cycle ID
sub get_ID{
	my $ucyc_OB = shift;
	return $ucyc_OB->{ID};
}

#--- return unit-cycle original ID
sub get_origID{
	my $ucyc_OB = shift;
	return $ucyc_OB->{origID};
}

#--- return unit-cycle length
sub get_length{
	my $ucyc_OB = shift;
	return $ucyc_OB->{length};
}

#--- bool value of FVGM
sub is_FVGM{
	my $ucyc_OB = shift;
	return $ucyc_OB->{fvgm}->[0];
}

#--- return FVGM ID
sub get_fvgmID{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID';

	# check
	unless( $ucyc_OB->is_FVGM ){
		warn_and_exit "<ERROR>\tUCYC $ucyc_OB->{$sayID} is not FVGM!\n";
	}
	return $ucyc_OB->{fvgm}->[1];
}

#--- return anchor SegID
sub get_anchorSegID{
	my $ucyc_OB = shift;
	return $ucyc_OB->{anchorSegID};
}

#--- return last SegID
sub get_lastSegID{
	my $ucyc_OB = shift;
	return $ucyc_OB->{lastSegID};
}

#--- return unit-cycle structure depth
sub get_structureDepth{
	my $ucyc_OB = shift;
	return $ucyc_OB->{structureDepth};
}

#--- show segment structure of this unit-cycle
sub show_segStructure{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID';

	my $SegStr = '';
	for my $CombEleHref ( @{$ucyc_OB->{segComb}} ){
		my $CombEleOB      = $CombEleHref->{OB};
		my $CombEleOBtype  = $CombEleHref->{OBtype};
		my $CombEleOBcount = $CombEleHref->{count};
		my $CombEleOBorit  = $CombEleHref->{orit};
		# segment
		if( $CombEleOBtype eq 'SEG' ){
			my $segID = $CombEleOB->get_ID;
			$SegStr .= "$segID$CombEleOBorit;" for ( 1 .. $CombEleOBcount ); # segment count must be one
		}
		# unit-cycle
		elsif( $CombEleOBtype eq 'UCYC' ){
			# check, all inner UCYC are required to have '+' orit
			if( $CombEleOBorit ne '+' ){
				warn_and_exit "<ERROR>\tInner UCYC ($CombEleOB->{$sayID}) has '-' orit in UCYC ($ucyc_OB->{$sayID})!\n";
			}
			# so, no need to convert the orit
			my $innerUcycSegStr = $CombEleOB->show_segStructure;
			$SegStr .= $innerUcycSegStr for ( 1 .. $CombEleOBcount );
		}
	}

	return $SegStr;
}

#--- display one unit-cycle
sub display{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $prefix = $parm{prefix} || '';
	my $count = $parm{count} || 1;
	my $recur_Count = $parm{recur_Count} || 0;
	my $subUcycIndtLen = $parm{subUcycIndtLen} || 4;
	my $subUcycIndt = ' ' x $subUcycIndtLen;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	# check recursion depth
	if( $recur_Count > 200 ){
		warn_and_exit "deep recursion ($recur_Count) in UCYC_OB->display function!\n";
	}

	my @segComb = @{$ucyc_OB->{segComb}};

	my $display_contect = '';
	my $lastCombEleOBtype = '';
	$display_contect .= $prefix.'(';
	for my $idx ( 0 .. $#segComb ){
		my $CombEleHref = $segComb[$idx];
		my $CombEleOBtype = $CombEleHref->{OBtype};
		my $CombEleOB = $CombEleHref->{OB};
		my $CombEleOrit  = $CombEleHref->{orit};
		my $CombEleCount = $CombEleHref->{count};
		# segment
		if( $CombEleOBtype eq 'SEG' ){
			$display_contect .= $prefix."|" if( $lastCombEleOBtype eq 'UCYC' );
			my $segID = $CombEleOB->get_ID;
			my $postBlank = (length($segID) <= 2 || $idx == $#segComb) ? ' ':'';
			$display_contect .= ' '.$segID.$CombEleOrit.$postBlank for ( 1 .. $CombEleCount );
		}
		# unit-cycle
		elsif( $CombEleOBtype eq 'UCYC' ){
			$display_contect .= "\n" if( $lastCombEleOBtype ne 'UCYC' );
			$display_contect .=
				$CombEleOB->display(
									prefix => $prefix.'|'.$subUcycIndt,
									count => $CombEleCount,
									subUcycIndtLen => $subUcycIndtLen,
									recur_Count => $recur_Count + 1,
									sayID => $sayID
								   );
		}
		# update
		$lastCombEleOBtype = $CombEleOBtype;
	}
	$display_contect .= $prefix if( $lastCombEleOBtype eq 'UCYC' );
	$display_contect .= '):' . $ucyc_OB->{$sayID} . ":$count\n";

	return $display_contect;
}

#--- check whether use anchorSeg for Basis
sub is_anchorSegForBasis{
	my $ucyc_OB = shift;
	return $ucyc_OB->{anchorSegForBasis};
}

#--- give two UCYC OB, to see whether they are identical
## identical: have same structure, use show_segStructure func
## two cases: Orignal or Complemented-Reversed
sub is_identical_UcycOB{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $test_UcycOB = $parm{test_UcycOB};

	my $this_ucyc_SegStr = $ucyc_OB->show_segStructure;

	# original tested one
	my $test_ucyc_SegStr = $test_UcycOB->show_segStructure;
	if( $this_ucyc_SegStr eq $test_ucyc_SegStr ){
		return 1;
	}

	# complemented and reversed version of tested one
	my $test_crUcycOB = $test_UcycOB->get_new_crUcycOB;
	my $test_crUcyc_SegStr = $test_crUcycOB->show_segStructure;
	if( $this_ucyc_SegStr eq $test_crUcyc_SegStr ){
		return -1;
	}

	# not identical
	return 0;
}

#--- find where given segID locates
sub has_this_segID{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $seekSegID = $parm{segID};
	my $seekMode = $parm{seekMode} || 'SEG';
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	# for 'UCYC' mode
	my $count = $parm{count}; # if set, the inner UCYC's count must be equal
	my $avoidUcycHref = $parm{avoidUcycHref}; # if set, must avoid

	# shuffle locNO in random mode
	my @locNO = ( 1 .. scalar(@{$ucyc_OB->{segComb}}) );
	@locNO = shuffle @locNO if( $randProbAref->[0] );
	# search
	for my $locNO ( @locNO ){
		my $CombEleHref = $ucyc_OB->{segComb}->[$locNO-1];
		next if( $CombEleHref->{OBtype} ne $seekMode );
		if( $seekMode eq 'SEG' ){
			my $segOB = $CombEleHref->{OB};
			my $segID = $segOB->get_ID;
			if( $segID eq $seekSegID ){
				return [ $locNO ];
			}
		}
		elsif( $seekMode eq 'UCYC' ){
			my $testUcycOB = $CombEleHref->{OB};
			my $subLocNO_Aref = $testUcycOB->has_this_segID(segID=>$seekSegID, seekMode=>'SEG', randProbAref=>$randProbAref);
			if(    $subLocNO_Aref->[0] != 0
				&& (   !defined $count
					|| $CombEleHref->{count} == $count
				   )
				&& (   !defined $avoidUcycHref
					|| !exists  $avoidUcycHref->{ $testUcycOB->get_origID }
				   )
			){
				return [ $locNO, $subLocNO_Aref->[0] ];
			}
		}
	}

	return [ 0 ];
}

#--- test whether has inner unit-cycle
## when given one inner UCYC OB, must find it
sub has_innerUCYC{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $innerUcycOB = $parm{innerUcycOB}; # optional
	my $randProbAref = $parm{randProbAref} || [0]; # random mode

	my @eleNO = ( 1 .. scalar(@{$ucyc_OB->{segComb}}) );
	@eleNO = shuffle( @eleNO ) if( $randProbAref->[0] );
	for my $eleNO (@eleNO) {
		my $CombEleHref = $ucyc_OB->{segComb}->[$eleNO-1];
		if(    $CombEleHref->{OBtype} eq 'UCYC' # must be UCYC
			&& (   !defined $innerUcycOB
				|| $CombEleHref->{OB}->is_identical_UcycOB( test_UcycOB => $innerUcycOB ) != 0
			   )
		){
			return $eleNO;
		}
	}

	return 0;
}

#--- make a copy of this UCYC OB
## use ID postfix '.cpN', means (C)o(P)y at the NO.N time
sub get_new_cpUcycOB{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $newID = $parm{newID};

	# copy NO control
	$ucyc_OB->{cpNO} ++;

	# prepare a new cp_UCYC object
	my $cp_UcycID = $newID || ( $ucyc_OB->get_ID . '.cp'. $ucyc_OB->{cpNO} );
	my $cp_UcycOB = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( ID=>$cp_UcycID, origID=>$ucyc_OB->get_origID );

	my @cp_SegComb = ();
	# iteratively make copy
	for my $CombEleHref (@{$ucyc_OB->{segComb}}){
		my $CombEleOBtype = $CombEleHref->{OBtype};
		my $CombEleOB     = $CombEleHref->{OB};
		my $CombEleOrit   = $CombEleHref->{orit};
		my $CombEleCount  = $CombEleHref->{count};
		# SEG or UCYC
		if( $CombEleOBtype eq 'SEG' ){
			# just copy the same
			push @cp_SegComb, { OBtype=>$CombEleOBtype, OB=>$CombEleOB, orit=>$CombEleOrit, count=>$CombEleCount };
		}
		elsif( $CombEleOBtype eq 'UCYC' ){
			# create cp_InnerUcycOB: copied
			my $cp_InnerUcycOB = $CombEleOB->get_new_cpUcycOB;
			push @cp_SegComb, { OBtype=>$CombEleOBtype, OB=>$cp_InnerUcycOB, orit=>$CombEleOrit, count=>$CombEleCount };
		}
	}

	# load info for cp_UcycOB
	$cp_UcycOB->{segComb} = \@cp_SegComb;
	$cp_UcycOB->update_info;

	return $cp_UcycOB;
}

#--- get new UCYC as complemented and reversed version of this UCYC
## use ID postfix '.cr', means (C)omplemented-(R)eversed
sub get_new_crUcycOB{
	my $ucyc_OB = shift;

	# prepare a new cr_UCYC object
	my $cr_UcycID = $ucyc_OB->get_ID . '.cr';
	my $cr_UcycOB = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( ID=>$cr_UcycID, origID=>$ucyc_OB->get_origID );

	my @cr_SegComb = ();
	# iteratively do compl_reversed
	for my $CombEleHref (@{$ucyc_OB->{segComb}}){
		my $CombEleOBtype = $CombEleHref->{OBtype};
		my $CombEleOB     = $CombEleHref->{OB};
		my $CombEleOrit   = $CombEleHref->{orit};
		my $CombEleCount  = $CombEleHref->{count};
		# SEG or UCYC
		if( $CombEleOBtype eq 'SEG' ){
			# reverse the orit
			my $cr_SegOrit = ($CombEleOrit eq '+' ? '-' : '+');
			unshift @cr_SegComb, { OBtype=>$CombEleOBtype, OB=>$CombEleOB, orit=>$cr_SegOrit, count=>$CombEleCount };
		}
		elsif( $CombEleOBtype eq 'UCYC' ){
			# create cr_InnerUcycOB: complemented and reversed
			my $cr_InnerUcycOB = $CombEleOB->get_new_crUcycOB;
			# here, I prefer to keep all inner UCYC '+' orit: let themselves to create new '.cr'
			unshift @cr_SegComb, { OBtype=>$CombEleOBtype, OB=>$cr_InnerUcycOB, orit=>'+', count=>$CombEleCount };
		}
	}

	# load info for cr_UcycOB
	$cr_UcycOB->{segComb} = \@cr_SegComb;
	$cr_UcycOB->update_info;

	return $cr_UcycOB;
}

#--- reverse valid inner-segments
sub reverseInnerSegs{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	# only work on full-seg UCYC
	next if( $ucyc_OB->has_innerUCYC );

	my @CombEleHref = @{$ucyc_OB->{segComb}};
	my $edit_bool = 0;
	my $min_idx_5p = 0;
	for my $idx_5p ( shuffle 1 .. $#CombEleHref-1 ){
		next unless( $idx_5p >= $min_idx_5p );
		my $CombEleHref_5pU = $CombEleHref[$idx_5p-1];
		my $SegID_5pU = $CombEleHref_5pU->{OB}->get_ID;
		my $Orit_5pU  = $CombEleHref_5pU->{orit};
		for my $idx_3p ( shuffle $idx_5p+1 .. $#CombEleHref-1 ){
			my $CombEleHref_3pD = $CombEleHref[$idx_3p+1];
			my $SegID_3pD = $CombEleHref_3pD->{OB}->get_ID;
			my $Orit_3pD = $CombEleHref_3pD->{orit};
			# must enclosed by same Seg but different orit ('+-/-+')
			# NOTE: we do not deal with same orit. but it might have
			if(    $SegID_5pU eq $SegID_3pD
				&&  $Orit_5pU ne  $Orit_3pD
			){
				next if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
				my @tempCombEleHref = reverse @CombEleHref[$idx_5p .. $idx_3p];
				# reverse inner segments
				for my $i (0 .. $#tempCombEleHref){
					my $i_CombEleHref = $tempCombEleHref[$i];
					my $r_Orit = ($i_CombEleHref->{orit} eq '+' ? '-' : '+');
					$CombEleHref[$idx_5p+$i] = { OBtype=>'SEG', OB=>$i_CombEleHref->{OB}, orit=>$r_Orit, count=>$i_CombEleHref->{count} };
				}
				# inform
				stout_and_sterr "[INFO]\t$ucyc_OB->{$sayID} InnerSeg reversed: idx.($idx_5p,$idx_3p).\n";
				# update
				$edit_bool = 1;
				# $min_idx_5p = $idx_3p + 1 if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
				$min_idx_5p = $idx_3p + 1;
				last;
			}
		}
	}

	if( $edit_bool ){
		# inform
		stout_and_sterr "[INFO]\tBEFORE: ".$ucyc_OB->show_segStructure."\n";
		# update
		$ucyc_OB->{segComb} = \@CombEleHref;
		$ucyc_OB->update_info;
		# inform
		stout_and_sterr "[INFO]\tAFTER:  ".$ucyc_OB->show_segStructure."\n";
	}
}

#--- change the unit-cycle structure to meet given anchor seg 
sub change_to_aimAnchorSeg{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $aimAnchorSeg = $parm{aimAnchorSeg};
	my $aimAnchorSegIsForBasis = $parm{aimAnchorSegIsForBasis} || 0;
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	my $informPref = $parm{informPref};
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	my $locNO_Aref = $ucyc_OB->has_this_segID( segID=>$aimAnchorSeg, seekMode=>'SEG', randProbAref=>$randProbAref );
	return if( $locNO_Aref->[0] == 0 );

	if( $ucyc_OB->get_anchorSegID ne $aimAnchorSeg ){
		push @{$ucyc_OB->{segComb}}, shift @{$ucyc_OB->{segComb}} for ( 1 .. $locNO_Aref->[0]-1 );
		# inform
		stout_and_sterr "[INFO]\t$informPref: Change anchorSeg of $ucyc_OB->{$sayID} to $aimAnchorSeg\n";
	}
	# update
	$ucyc_OB->{anchorSegForBasis} = ( $aimAnchorSegIsForBasis ? 1 : 0 );
	$ucyc_OB->update_info;
}

#--- break inner UCYC containing only one element
sub breakSingleEleUCYC{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID';

	return unless( $ucyc_OB->has_innerUCYC );

	for my $CombEleHref ( @{$ucyc_OB->{segComb}} ){
		next if( $CombEleHref->{OBtype} ne 'UCYC' ); # must be innerUcyc
		my $innerUcycOB = $CombEleHref->{OB};
		# do same to this innerUcyc
		$innerUcycOB->breakSingleEleUCYC( sayID => $sayID );
		# test itself
		my @innerUcycSegComb = @{$innerUcycOB->{segComb}};
		if( scalar(@innerUcycSegComb) == 1 ){
			# update this CombEleHref
			$CombEleHref->{OBtype} = $innerUcycSegComb[0]->{OBtype};
			$CombEleHref->{OB}     = $innerUcycSegComb[0]->{OB};
			$CombEleHref->{orit}   = ( $CombEleHref->{orit} eq $innerUcycSegComb[0]->{orit} ? '+' : '-' );
			$CombEleHref->{count} *= $innerUcycSegComb[0]->{count};
			# inform
			stout_and_sterr "[INFO]\tBreak SingleEle inner unit-cycle ($innerUcycOB->{$sayID}) in upper unit-cycle ($ucyc_OB->{$sayID})\n";
		}
	}
}

#--- absorb unit-cycle
## return the totally used count of the absbUcyc
sub absorbUcyc{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $thisUcycOrigCount = $parm{thisUcycOrigCount} || 1;
	my $absbUcycOB = $parm{absbUcycOB};
	my $absbUcycOrigCount = $parm{absbUcycOrigCount};
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	my $avoidUcycHref = $parm{avoidUcycHref} || {};
	my $sayID = $parm{sayID} || 'origID';

	# prepare the count of the absbUcyc could use in single copy of this ucyc_OB
	my $absbUcycCountToUse = int( $absbUcycOrigCount / $thisUcycOrigCount );
	return 0 if( $absbUcycCountToUse == 0 );

	# anchorSegID to find
	my $absbUcycAnchorSegID = $absbUcycOB->get_anchorSegID;

	# try to absort the absbUcyc at this ucyc_OB
	## 'SEG' mode: find matched 'SEG' in this ucyc_OB
	## 'UCYC' mode: find matched 'SEG' in innerUCYC of this ucyc_OB
	my @seekMode = ( 'SEG', 'UCYC' );
	if( $randProbAref->[0] ){
		# 'deepUCYC' means just go to innerUCYC of this UCYC
		@seekMode = GetRandBool(verbose=>0,Aref=>$randProbAref) ? ('deepUCYC') : shuffle( @seekMode );
		# inform
		stout_and_sterr "[INFO]\tShuffle/Reset seekMode to " . join(',',@seekMode) . "\n";
	}
	for my $seekMode ( @seekMode ){
		my @has_this_segID_parm = ( segID => $absbUcycAnchorSegID, seekMode => $seekMode, avoidUcycHref => $avoidUcycHref, randProbAref => $randProbAref );
		# SEG mode
		if( $seekMode eq 'SEG' ){
			# find the anchor seg in this UCYC
			my $locNO_Aref = $ucyc_OB->has_this_segID( @has_this_segID_parm );
			# can not find, try another seekMode
			next if( $locNO_Aref->[0] == 0 );
			# insert the absbUcyc with certain count(s)
			my $innerSegLocNO = $locNO_Aref->[0];
			my $absorbCount = ( GetRandBool(verbose=>0,Aref=>$randProbAref) ? (int(rand($absbUcycCountToUse))+1) : $absbUcycCountToUse );
			$ucyc_OB->insertUcycToReplaceGivenSeg( sayID => $sayID, insLocNO => $innerSegLocNO,
												   insUcycOB => $absbUcycOB, insCount => $absorbCount );
			# inform
			stout_and_sterr "[INFO]\tMerge $absbUcycOB->{$sayID} (anchor:$absbUcycAnchorSegID) to upper unit-cycle ($ucyc_OB->{$sayID})"
								  ." at NO.$innerSegLocNO loc with $absorbCount*$thisUcycOrigCount count.\n";
			# return totally used count of the absbUcyc
			return $absorbCount * $thisUcycOrigCount;
		}
		# UCYC mode
		elsif( $seekMode eq 'UCYC' ){
			# find the anchor seg in inner UCYC of basis
			my $locNO_Aref = [0];
			## try #1: equal count (only in non-random mode); then,
			## try #2: any count
			for my $try ( 1, 2 ){
				if( $try == 1 && !$randProbAref->[0] ){
					$locNO_Aref = $ucyc_OB->has_this_segID( @has_this_segID_parm, count => $absbUcycCountToUse );
				}
				else{
					$locNO_Aref = $ucyc_OB->has_this_segID( @has_this_segID_parm );
				}
				# stop when found
				last if( $locNO_Aref->[0] != 0 );
			}
			# can not find, try another seekMode
			next if( $locNO_Aref->[0] == 0 );
			# insert the absbUcyc to the innerUCYC
			my ($innerUcycLocNO, $innerUcycInsLocNO) = @$locNO_Aref;
			my $innerUcycOB    = $ucyc_OB->{segComb}->[$innerUcycLocNO-1]->{OB};
			my $innerUcycCount = $ucyc_OB->{segComb}->[$innerUcycLocNO-1]->{count};
			## based on the Counts
			if( $innerUcycCount <= $absbUcycCountToUse ){
				# inserted the absbUcyc, embed (as repeat) in innerUCYC
				my $absorbCountSingleTime = int( $absbUcycCountToUse / $innerUcycCount );
				   $absorbCountSingleTime = int( rand($absorbCountSingleTime) ) + 1 if( GetRandBool(verbose=>0,Aref=>$randProbAref) );
				$innerUcycOB->insertUcycToReplaceGivenSeg( sayID => $sayID, insLocNO => $innerUcycInsLocNO,
														   insUcycOB => $absbUcycOB, insCount => $absorbCountSingleTime );
				# update upper ucyc
				$ucyc_OB->update_info;
				# inform
				stout_and_sterr "[INFO]\tMerge $absbUcycOB->{$sayID} (anchor:$absbUcycAnchorSegID) to upper unit-cycle ($ucyc_OB->{$sayID})"
									  ." at NO.$innerUcycLocNO($innerUcycOB->{$sayID}) NO.$innerUcycInsLocNO loc with $absorbCountSingleTime*$innerUcycCount*$thisUcycOrigCount count.\n";
				# return totally used count of the absbUcyc
				my $absorbCount = $absorbCountSingleTime * $innerUcycCount;
				return $absorbCount * $thisUcycOrigCount;
			}
			else{
				# the innerUCYC should be alloted to two parts
				## note that the first  part (still be the same OB) is for absbUcyc to insert;
				##    while, the second part is a new copied innerUCYC.
				my $absorbCount = ( GetRandBool(verbose=>0,Aref=>$randProbAref) ? (int(rand($absbUcycCountToUse))+1) : $absbUcycCountToUse );
				$ucyc_OB->allot_innerUcyc( dvdLocNO => $innerUcycLocNO, firstPartCount => $absorbCount );
				# inform
				stout_and_sterr "[INFO]\tAllot inner unit-cycle ($innerUcycOB->{$sayID};locNO.$innerUcycLocNO) of upper unit-cycle ($ucyc_OB->{$sayID})"
									  ." to two inner unit-cycles (firstPartCount=$absorbCount).\n";
				# insert the absbUcyc to the first alloted part
				$innerUcycOB->insertUcycToReplaceGivenSeg( sayID => $sayID, insLocNO => $innerUcycInsLocNO,
														   insUcycOB => $absbUcycOB, insCount => 1			);
				# update upper ucyc
				$ucyc_OB->update_info;
				# inform
				stout_and_sterr "[INFO]\tMerge $absbUcycOB->{$sayID} (anchor:$absbUcycAnchorSegID) to upper unit-cycle ($ucyc_OB->{$sayID})"
									  ." at NO.$innerUcycLocNO($innerUcycOB->{$sayID}) NO.$innerUcycInsLocNO loc with 1*$absorbCount*$thisUcycOrigCount count.\n";
				# return totally used count of the absbUcyc
				return $absorbCount * $thisUcycOrigCount;
			}
		}
	}

	# if reach this, means cannot absorb the absbUcyc at this ucyc_OB
	# so, try to apply absorbUcyc func of innerUcyc(s) of this ucyc_OB
	for my $CombEleHref ( @{$ucyc_OB->{segComb}} ){
		next if( $CombEleHref->{OBtype} ne 'UCYC' ); # must be an innerUCYC
		my $innerUcycOB = $CombEleHref->{OB};
		my $innerUcycOrigCount = $CombEleHref->{count} * $thisUcycOrigCount;
		my $absorbTotalCount = $innerUcycOB->absorbUcyc( thisUcycOrigCount => $innerUcycOrigCount,
														 randProbAref => $randProbAref, sayID => $sayID,
														 absbUcycOB => $absbUcycOB,
														 absbUcycOrigCount => $absbUcycOrigCount,
														 avoidUcycHref => $avoidUcycHref           );
		return $absorbTotalCount if( $absorbTotalCount != 0 );
	}

	# cannot insert the absorbUcyc 
	return 0;
}

#--- insert unit-cycle at given loc to replace SEG
sub insertUcycToReplaceGivenSeg{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $insUcycOB = $parm{insUcycOB};
	my $insLocNO = $parm{insLocNO}; # count from 1, not 0
	my $insCount = $parm{insCount};
	my $sayID = $parm{sayID} || 'origID';

	# must be 'SEG'
	if( $ucyc_OB->{segComb}->[$insLocNO-1]->{OBtype} ne 'SEG' ){
		warn_and_exit "<ERROR>\tThe element that inserted unit-cycle ($insUcycOB->{$sayID}) to replace is not SEG type"
							  ." at NO.$insLocNO in upper unit-cycle $ucyc_OB->{$sayID}.\n";
	}

	my @newSegComb = ();
	my $insUcycID;
	# copy elements before insLocNO
	push @newSegComb, $ucyc_OB->{segComb}->[$_] for (0 .. $insLocNO-2);
	# next push the inserted UCYC based on the orientation
	my $insUcycAnchorSegOrit = $insUcycOB->{segComb}->[0]->{orit};
	my $recUcycInsLocSegOrit = $ucyc_OB->{segComb}->[$insLocNO-1]->{orit};
	# diiferent orit?
	if( $insUcycAnchorSegOrit ne $recUcycInsLocSegOrit ){ # diff
		# create crUcycOB: complemented and reversed
		my $cr_UcycOB = $insUcycOB->get_new_crUcycOB;
		$insUcycID = $cr_UcycOB->get_ID;
		# put the ins-ucyc behind the insLocNO
		push @newSegComb, $ucyc_OB->{segComb}->[$insLocNO-1];
		push @newSegComb, { OBtype=>'UCYC', OB=>$cr_UcycOB, orit=>'+', count=>$insCount };
		push @newSegComb, $ucyc_OB->{segComb}->[$_] for ($insLocNO .. scalar(@{$ucyc_OB->{segComb}})-1);
	}
	else{ # same
		# create cpUcycOB: copy of insUcycOB
		my $cp_UcycOB = $insUcycOB->get_new_cpUcycOB;
		$insUcycID = $cp_UcycOB->get_ID;
		# put the ins-ucyc before the insLocNO
		push @newSegComb, { OBtype=>'UCYC', OB=>$cp_UcycOB, orit=>'+', count=>$insCount };
		push @newSegComb, $ucyc_OB->{segComb}->[$_] for ($insLocNO-1 .. scalar(@{$ucyc_OB->{segComb}})-1);
	}

	# update
	$ucyc_OB->{ID} .= ".i$insLocNO\[$insUcycID\]";
	$ucyc_OB->{segComb} = \@newSegComb;
	$ucyc_OB->update_info;
}

#--- allot the inner unit-cycle at given loc
## use ID postfix '.a2', means the 2nd alloted group
## the first group keep to be the original one
sub allot_innerUcyc{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $dvdLocNO = $parm{dvdLocNO};
	my $firstPartCount = $parm{firstPartCount};

	if( $ucyc_OB->{segComb}->[$dvdLocNO-1]->{OBtype} ne 'UCYC' ){
		warn_and_exit "<ERROR>\tNot unit-cycle at NO.$dvdLocNO in unit cyle $ucyc_OB->{ID}.";
	}

	my @newSegComb = ();
	# copy elements before dvdLocNO
	push @newSegComb, $ucyc_OB->{segComb}->[$_] for (0 .. $dvdLocNO-2);

	# allot
	my $ToDvdUcycHref = $ucyc_OB->{segComb}->[$dvdLocNO-1];
	my $ToDvdUcycOrit  = $ToDvdUcycHref->{orit};
	my $ToDvdUcycOB    = $ToDvdUcycHref->{OB};
	my $ToDvdUcycCount = $ToDvdUcycHref->{count};

	## first part is still $ToDvdUcycOB but with less counts
	my $firstPartUcycOB = $ToDvdUcycOB;
	push @newSegComb, { OBtype=>'UCYC', OB=>$firstPartUcycOB, orit=>$ToDvdUcycOrit, count=>$firstPartCount };

	## second part is one NEW unit-cycle OB, just copy, and set new ID
	my $secondPartUcycID = $ToDvdUcycOB->get_ID . '.a2'; # '.a2' means the 2nd alloted group.
	my $secondPartUcycOB = $ToDvdUcycOB->get_new_cpUcycOB( newID => $secondPartUcycID );
	my $secondPartCount  = $ToDvdUcycCount - $firstPartCount;
	push @newSegComb, { OBtype=>'UCYC', OB=>$secondPartUcycOB, orit=>$ToDvdUcycOrit, count=>$secondPartCount };

	# copy elements after dvdLocNO
	push @newSegComb, $ucyc_OB->{segComb}->[$_] for ($dvdLocNO .. scalar(@{$ucyc_OB->{segComb}})-1);

	# update
	$ucyc_OB->{segComb} = \@newSegComb;
	$ucyc_OB->update_info;
}

#--- try to find whether it has the given UCYC
## if yes, absorb certain counts
## use ID prefix: '.asNXXX', means (A)dd (S)ame UCYC (ID:XXX) at NO.N local
sub absorbSameUcycItHas{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $testUcycOB = $parm{testUcycOB};
	my $testUcycCount = $parm{testUcycCount};
	my $randProbAref = $parm{randProbAref} || [0]; # random mode
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	# exceptions
	return [ 0 ] if( $testUcycCount == 0 );

	my $addBool = 0;
	my $addCount = 0;
	# check each elements of the basic UCYC
	my @eleNO = ( 1 .. scalar(@{$ucyc_OB->{segComb}}) );
	@eleNO = shuffle( @eleNO ) if( $randProbAref->[0] );
	for my $eleNO (@eleNO) {
		my $CombEleHref = $ucyc_OB->{segComb}->[$eleNO-1];
		next if( $CombEleHref->{OBtype} ne 'UCYC' ); # muse be an innerUCYC
		my $innerUcycOB = $CombEleHref->{OB};
		my $innerUcycCount = $CombEleHref->{count};
		my $identTestSign = $testUcycOB->is_identical_UcycOB( test_UcycOB => $innerUcycOB );
		# identical?
		if( $identTestSign != 0 ){ # yes, but no matter what orit
			$addCount = ( GetRandBool(verbose=>0,Aref=>$randProbAref) ? (int(rand($testUcycCount))+1) : $testUcycCount );
			$CombEleHref->{count} += $addCount;
			$ucyc_OB->{ID} .= ".as$eleNO".$testUcycOB->{ID}; # '.as' mean (A)dd (S)ame UCYC
			$addBool = 1;
		}
		else{ # not identical, then check whether this UCYC contains this testUcyc
			my $addCountSingleTime = int( $testUcycCount / $innerUcycCount ); # even 0 is ok, will check.
			my $asInfoAref = $innerUcycOB->absorbSameUcycItHas( testUcycOB => $testUcycOB, testUcycCount => $addCountSingleTime,
																randProbAref => $randProbAref, sayID => $sayID );
			next if( $asInfoAref->[0] == 0 );
			$addCount = $asInfoAref->[1] * $innerUcycCount;
			$addBool = 1;
		}
		# after add, should update_info, and update remain counts
		if( $addBool == 1 ){
			# inform
			stout_and_sterr "[INFO]\tAdd $testUcycOB->{$sayID} in unit-cycle ($ucyc_OB->{$sayID})"
								  ." at NO.$eleNO loc with $addCount count.\n";
			# update info
			$ucyc_OB->update_info;
			last; # not need to go on searching
		}
	}

	return [ $addBool, $addCount ];
}

#--- merge inner UCYC with flanking Element's Str (same as inner UCYC)
sub merge_innerUcycWithFlankEleStr{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	return unless( $ucyc_OB->has_innerUCYC );

	UCYC_FLANK_ELE:{
		my @segComb = @{$ucyc_OB->{segComb}};
		for my $idx ( 0 .. $#segComb ){
			my $CombEleHref = $segComb[$idx];
			next if( $CombEleHref->{OBtype} ne 'UCYC' );
			# inner UCYC should always be '+'
			next if( $CombEleHref->{orit} ne  '+' );
			# first check this innerUcyc
			my $innerUcycOB = $CombEleHref->{OB};
			$innerUcycOB->merge_innerUcycWithFlankEleStr( sayID => $sayID );
			# get the SegStr of the innerUcyc
			my $innerUcycSegStr = $innerUcycOB->show_segStructure;
			# forward (F) check
			my @F_SegStr = ();
			for my $F_idx ( reverse (0 .. $idx-1) ){
				my $F_CombEleOB      = $segComb[$F_idx]->{OB};
				my $F_CombEleOBtype  = $segComb[$F_idx]->{OBtype};
				my $F_CombEleOBcount = $segComb[$F_idx]->{count};
				if( $F_CombEleOBtype eq 'SEG' ){
					unshift @F_SegStr, ($F_CombEleOB->get_ID . $segComb[$F_idx]->{orit} . ';') x $F_CombEleOBcount;
				}
				else{
					unshift @F_SegStr, $F_CombEleOB->show_segStructure x $F_CombEleOBcount;
				}
				# test match, F_SegStr should match tail of the innerUcycSegStr or whole body (N times)
				my $F_repTime = getStrUnitRepeatTime( StrTest => join('', @F_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'end' );
				if( $F_repTime == 0 ){
					# discard one by one from last introduced one
					while( scalar(@F_SegStr) != 0 ){
						shift @F_SegStr;
						$F_repTime = getStrUnitRepeatTime( StrTest => join('', @F_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'end', MustInt => 1 );
						last if( $F_repTime != 0 ); # find integer-repeat(s)
					}
					# stop
					last;
				}
			}
			## F match?
			my $F_repTime = getStrUnitRepeatTime( StrTest => join('', @F_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'end', MustInt => 1 );
			if( $F_repTime != 0 ){
				my @newSegComb;
				my $F_st_idx = $idx - scalar(@F_SegStr);
				push @newSegComb, $segComb[$_] for (0 .. $F_st_idx-1);
				push @newSegComb, $segComb[$idx];
				$newSegComb[-1]->{count} += $F_repTime;
				push @newSegComb, $segComb[$_] for ($idx+1 .. $#segComb);
				# inform
				stout_and_sterr "[INFO]\tMerge SegStr ($innerUcycSegStr)*$F_repTime with its back-nearby inner unit-cycle ($innerUcycOB->{$sayID}) in upper unit-cycle ($ucyc_OB->{$sayID})\n";
				# update
				$ucyc_OB->{segComb} = \@newSegComb;
				$ucyc_OB->update_info;
				# redo
				redo UCYC_FLANK_ELE;
			}
			# backword (B) check
			my @B_SegStr = ();
			for my $B_idx ( $idx+1 .. $#segComb ){
				my $B_CombEleOB      = $segComb[$B_idx]->{OB};
				my $B_CombEleOBtype  = $segComb[$B_idx]->{OBtype};
				my $B_CombEleOBcount = $segComb[$B_idx]->{count};
				if( $B_CombEleOBtype eq 'SEG' ){
					push @B_SegStr, ($B_CombEleOB->get_ID . $segComb[$B_idx]->{orit} . ';') x $B_CombEleOBcount;
				}
				else{
					push @B_SegStr, $B_CombEleOB->show_segStructure x $B_CombEleOBcount;
				}
				# test match, B_SegStr should match head of the innerUcycSegStr or whole body (N times)
				my $B_repTime = getStrUnitRepeatTime( StrTest => join('', @B_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'start' );
				if( $B_repTime == 0 ){
					while( scalar(@B_SegStr) != 0 ){
						# discard one by one from last introduced one
						pop @B_SegStr;
						$B_repTime = getStrUnitRepeatTime( StrTest => join('', @B_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'start', MustInt => 1 );
						last if( $B_repTime != 0 ); # find integer-repeat(s)
					}
					# stop
					last;
				}
			}
			## B match?
			my $B_repTime = getStrUnitRepeatTime( StrTest => join('', @B_SegStr), StrUnit => $innerUcycSegStr, TestEnd => 'start', MustInt => 1 );
			if( $B_repTime != 0 ){
				my @newSegComb;
				push @newSegComb, $segComb[$_] for (0 .. $idx-1);
				push @newSegComb, $segComb[$idx];
				$newSegComb[-1]->{count} += $B_repTime;
				my $B_ed_idx = $idx + scalar(@B_SegStr);
				push @newSegComb, $segComb[$_] for ($B_ed_idx+1 .. $#segComb);
				# inform
				stout_and_sterr "[INFO]\tMerge SegStr ($innerUcycSegStr)*$B_repTime with its forw-nearby inner unit-cycle ($innerUcycOB->{$sayID}) in upper unit-cycle ($ucyc_OB->{$sayID})\n";
				# update
				$ucyc_OB->{segComb} = \@newSegComb;
				$ucyc_OB->update_info;
				# redo
				redo UCYC_FLANK_ELE;
			}
		}
	}
}

#--- let single copy inner UCYC to throw out sequential host-Seg from tail --
## use ID postfix '.tN', means (T)hrow N segments
sub throwTailSegFromSingleInnerUCYC{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	return unless( $ucyc_OB->has_innerUCYC );

	THROW_TAIL_HOST_SEG: {
		my @segComb = @{$ucyc_OB->{segComb}};
		for my $idx ( 1 .. $#segComb-1 ){
			my $CombEleHref = $segComb[$idx];
			next if( $CombEleHref->{OBtype} ne 'UCYC' ); # must be inner UCYC
			my $innerUcycOB = $CombEleHref->{OB};
			# firstly, do this func to this inner UCYC
			$innerUcycOB->throwTailSegFromSingleInnerUCYC( sayID => $sayID );
			# then go on do this func for this UCYC
			my $innerUcycOrit = $CombEleHref->{orit};
			next if( $innerUcycOrit ne  '+' ); # inner UCYC should be all '+'
			my $innerUcycCount = $CombEleHref->{count};
			next if( $innerUcycCount != 1 ); # only count=1 UCYC could throw 1-copy seg
			my @innerUcycSegComb = @{$innerUcycOB->{segComb}};
			# from tail to find sequential host-Seg following non-sequential host-Seg
			my $throw_idx; # throw the host-seg from this idx
			my $throw_bool = 0; # to throw or not
			for my $innerIdx (reverse 1 .. $#innerUcycSegComb ){ # reverse: start from tail
				# this innerCombEle
				my $this_innerCombEleHref = $innerUcycSegComb[$innerIdx];
				last if($this_innerCombEleHref->{OBtype} ne 'SEG'); # must be Seg
				last if($this_innerCombEleHref->{orit} ne '+'); # must forward orit
				my $this_innerSegOB = $this_innerCombEleHref->{OB};
				last if($this_innerSegOB->get_SpecTag ne 'H'); # must be host
				# last innerCombEle
				my $last_innerCombEleHref = $innerUcycSegComb[$innerIdx-1];
				my $last_innerSegOB = $last_innerCombEleHref->{OB}; # must be Seg, but check later
				## check
				if(    $last_innerCombEleHref->{OBtype} ne 'SEG' # must be Seg
					|| (   $last_innerSegOB->get_SpecTag eq 'H' # be host-seg
						 && (    $last_innerCombEleHref->{orit} ne '+' # but not forward orit
							  || ! $this_innerSegOB->test_sequential_seg( testSegOB=>$last_innerSegOB, testSegPrime=>'5p' ) # not sequential
						    )
					   )
				){
					$throw_idx = $innerIdx;
					$throw_bool = 1;
					last;
				}
				elsif( $last_innerSegOB->get_SpecTag ne 'H' ){ # not host, should not throw
					last;
				}
			}
			# can throw out sequential host-Seg
			if( $throw_bool ){
				my @throwSegComb = @innerUcycSegComb[ $throw_idx .. $#innerUcycSegComb ];
				my $throwSegCount = scalar(@throwSegComb);
				# CombEle before this inner UCYC
				my @newSegComb = @segComb[ 0 .. $idx-1 ];
				# create new inner UCYC
				my $newInnerUcycID = $innerUcycOB->get_ID . ".t$throwSegCount"; # '.tN' means (T)hrow N segments
				my $newInnerUcycOB = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( ID=>$newInnerUcycID, origID=>$innerUcycOB->get_origID );
				my @newInnerUcycSegComb = @innerUcycSegComb[ 0 .. $throw_idx-1 ];
				$newInnerUcycOB->{segComb} = \@newInnerUcycSegComb;
				$newInnerUcycOB->update_info;
				push @newSegComb, { OBtype=>'UCYC', OB=>$newInnerUcycOB, orit=>$innerUcycOrit, count=>$innerUcycCount };
				# host-seg throw
				push @newSegComb, @throwSegComb;
				# CombEle after this inner UCYC
				push @newSegComb, @segComb[ $idx+1 .. $#segComb ];
				# update
				$ucyc_OB->{segComb} = \@newSegComb;
				$ucyc_OB->update_info;
				# inform
				stout_and_sterr "[INFO]\tThrow $throwSegCount sequential host-Seg from tail of inner unit-cycle ($innerUcycOB->{$sayID};idx.$idx) in upper unit-cycle ($ucyc_OB->{$sayID}).\n";
				# run again
				redo THROW_TAIL_HOST_SEG;
			}
		}
	}
}

#--- let single copy inner UCYC to absorb flanking sequential Seg
## use ID postfix '.b5', means absorb seg from 5-prime
## use ID postfix '.b3', means absorb seg from 3-prime
sub absorbSegToSingleInnerUCYC{
	my $ucyc_OB = shift;
	my %parm = @_;
	my $sayID = $parm{sayID} || 'origID'; # or 'ID', longer!

	return unless( $ucyc_OB->has_innerUCYC );

	ABSORB_SEG:{
		my @segComb = @{$ucyc_OB->{segComb}};
		for my $idx ( 1 .. $#segComb-1 ){
			my $CombEleHref = $segComb[$idx];
			next if( $CombEleHref->{OBtype} ne 'UCYC' );
			my $innerUcycOB = $CombEleHref->{OB};
			# firstly, do this func to this inner UCYC
			$innerUcycOB->absorbSegToSingleInnerUCYC( sayID => $sayID );
			# then go on do this func for this UCYC
			my $innerUcycOrit  = $CombEleHref->{orit};
			next if( $innerUcycOrit ne  '+' ); # inner UCYC should be all '+'
			my $innerUcycCount = $CombEleHref->{count};
			next if( $innerUcycCount != 1 ); # only count=1 UCYC could absorb 1-copy seg
			# 'A'nchor
			my $innerUcycAsegOB = $innerUcycOB->{segComb}->[0]->{OB};
			my $innerUcycAsegOBtype = $innerUcycOB->{segComb}->[0]->{OBtype};
			my $innerUcycAsegOrit = $innerUcycOB->{segComb}->[0]->{orit};
			# 'L'ast
			my $innerUcycLsegOB = $innerUcycOB->{segComb}->[-1]->{OB};
			my $innerUcycLsegOBtype = $innerUcycOB->{segComb}->[-1]->{OBtype};
			my $innerUcycLsegOrit = $innerUcycOB->{segComb}->[-1]->{orit};
			# check 5-prime SEG
			my $idx_5p_max = ( $innerUcycAsegOBtype eq 'SEG' ? $idx-1 : -10 );
			my $idx_5p = undef;
			my @absorbSeg5p;
			my $flk_3p_segOB = $innerUcycAsegOB;
			for my $idx_5p_cand (reverse 1 .. $idx_5p_max){ # 1: avoid the anchor seg
				my $CombEleHref_5p = $segComb[$idx_5p_cand];
				last unless( $CombEleHref_5p->{OBtype} eq 'SEG' ); # must be SEG
				my $seg5pOB = $CombEleHref_5p->{OB};
				last unless( $seg5pOB->get_SpecTag eq $flk_3p_segOB->get_SpecTag ); # same species
				my $seg5pOrit = $CombEleHref_5p->{orit};
				last unless( $seg5pOrit eq $innerUcycAsegOrit ); # same orit with Aseg (always)
				if(    ( $seg5pOrit eq '+' && $flk_3p_segOB->test_stream_seg(testSegOB=>$seg5pOB, testSegPrime=>'5p') )
					|| ( $seg5pOrit eq '-' && $flk_3p_segOB->test_stream_seg(testSegOB=>$seg5pOB, testSegPrime=>'3p') )
				){
					$idx_5p = $idx_5p_cand;
					$flk_3p_segOB = $seg5pOB;
					unshift @absorbSeg5p, $seg5pOB->get_ID.$seg5pOrit;
				}
				else{
					last;
				}
			}
			if( defined $idx_5p ){
				my @newSegComb = @segComb[ 0 .. $idx_5p-1 ];
				# after absorbing, the innerUCYC should be replaced by novel UCYC
				my $newInnerUcycID = $innerUcycOB->get_ID . '.b5'; # '.b5' means absorb seg from 5-prime
				my $newInnerUcycOB = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( ID=>$newInnerUcycID, origID=>$innerUcycOB->get_origID );
				my @temp = @{$innerUcycOB->{segComb}};
				unshift @temp, { OBtype=>'SEG', OB=>$segComb[$_]->{OB}, orit=>$segComb[$_]->{orit}, count=>1 } for reverse $idx_5p .. $idx_5p_max;
				$newInnerUcycOB->{segComb} = \@temp;
				$newInnerUcycOB->update_info;
				push @newSegComb, { OBtype=>'UCYC', OB=>$newInnerUcycOB, orit=>$innerUcycOrit, count=>$innerUcycCount };
				# take followings
				push @newSegComb, @segComb[ $idx+1 .. $#segComb ];
				# update
				$ucyc_OB->{segComb} = \@newSegComb;
				$ucyc_OB->update_info;
				# inform
				my $absorbSeg5p = join(';', @absorbSeg5p);
				stout_and_sterr "[INFO]\tInner unit-cycle ($innerUcycOB->{$sayID};idx.$idx) absorbs its 5-prime Seg(s) ($absorbSeg5p) from upper unit-cycle ($ucyc_OB->{$sayID}).\n";
				# run again
				redo ABSORB_SEG;
			}
			# check 3-prime SEG
			my $idx_3p_min = ( $innerUcycLsegOBtype eq 'SEG' ? $idx+1 : 1E6 ); # 1E6: is enough?! a-ha~
			my $idx_3p = undef;
			my @absorbSeg3p;
			my $flk_5p_segOB = $innerUcycLsegOB;
			for my $idx_3p_cand ( $idx_3p_min .. $#segComb-1 ){ # $#segComb-1: avoid the last seg
				my $CombEleHref_3p = $segComb[$idx_3p_cand];
				last unless( $CombEleHref_3p->{OBtype} eq 'SEG' ); # must be SEG
				my $seg3pOB = $CombEleHref_3p->{OB};
				last unless( $seg3pOB->get_SpecTag eq $flk_5p_segOB->get_SpecTag ); # same species
				my $seg3pOrit = $CombEleHref_3p->{orit};
				last unless( $seg3pOrit eq $innerUcycLsegOrit ); # same orit with Lseg (always)
				if(    ( $seg3pOrit eq '+' && $flk_5p_segOB->test_stream_seg(testSegOB=>$seg3pOB, testSegPrime=>'3p') )
					|| ( $seg3pOrit eq '-' && $flk_5p_segOB->test_stream_seg(testSegOB=>$seg3pOB, testSegPrime=>'5p') )
				){
					$idx_3p = $idx_3p_cand;
					$flk_5p_segOB = $seg3pOB;
					push @absorbSeg3p, $seg3pOB->get_ID.$seg3pOrit;
				}
				else{
					last;
				}
			}
			if( defined $idx_3p ){
				my @newSegComb = @segComb[ 0 .. $idx-1 ];
				# after absorbing, the innerUCYC should be replaced by novel UCYC
				my $newInnerUcycID = $innerUcycOB->get_ID . '.b3'; # '.b3' means absorb seg from 3-prime
				my $newInnerUcycOB = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( ID=>$newInnerUcycID, origID=>$innerUcycOB->get_origID );
				my @temp = @{$innerUcycOB->{segComb}};
				push @temp, { OBtype=>'SEG', OB=>$segComb[$_]->{OB}, orit=>$segComb[$_]->{orit}, count=>1 } for $idx_3p_min .. $idx_3p;
				$newInnerUcycOB->{segComb} = \@temp;
				$newInnerUcycOB->update_info;
				push @newSegComb, { OBtype=>'UCYC', OB=>$newInnerUcycOB, orit=>$innerUcycOrit, count=>$innerUcycCount };
				# take followings
				push @newSegComb, @segComb[ $idx_3p+1 .. $#segComb ];
				# update
				$ucyc_OB->{segComb} = \@newSegComb;
				$ucyc_OB->update_info;
				# inform
				my $absorbSeg3p = join(';', @absorbSeg3p);
				stout_and_sterr "[INFO]\tInner unit-cycle ($innerUcycOB->{$sayID};idx.$idx) absorbs its 3-prime Seg(s) ($absorbSeg3p) from upper unit-cycle ($ucyc_OB->{$sayID}).\n";
				# run again
				redo ABSORB_SEG;
			}
		}
	}
}

1; ## tell the perl script the successful access of this module.
