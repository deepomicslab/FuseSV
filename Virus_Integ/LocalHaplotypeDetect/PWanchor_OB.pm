package FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'PWanchor_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-04-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						get_ID
						get_groupOB
						get_linkCount
						get_sideInfo
						has_rel_segID
						find_anchor_idx
						get_anchor_distance
					 /;

#--- structure of object
# pwanchor_OB -> ID = $ID
# pwanchor_OB -> linkCount = $linkCount
# pwanchor_OB -> info = { a -> $info_a, b -> $info_b }
# pwanchor_OB -> type = { a -> $type_a, b -> $type_a } : JUNC/POS
# pwanchor_OB -> obj = { a -> $obj_a, b -> $obj_b }
# pwanchor_OB -> groupOB = $groupOB
# pwanchor_OB -> rel_segID = { $segID -> {$side->1, .. }, .. }

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $ID = $parm{ID};
	my $linkCount = $parm{linkCount};
	my $info_a = $parm{info_a};
	my $info_b = $parm{info_b};
	my $type_a = $parm{type_a};
	my $type_b = $parm{type_b};
	my $obj_a  = $parm{obj_a};
	my $obj_b  = $parm{obj_b};

	my $pwanchor_OB = {};

	$pwanchor_OB->{ID} = $ID;
	$pwanchor_OB->{linkCount} = $linkCount;
	$pwanchor_OB->{info}->{a} = $info_a;
	$pwanchor_OB->{info}->{b} = $info_b;
	$pwanchor_OB->{type}->{a} = $type_a;
	$pwanchor_OB->{type}->{b} = $type_b;
	$pwanchor_OB->{obj}->{a}  = $obj_a;
	$pwanchor_OB->{obj}->{b}  = $obj_b;

	# check
	if(    $type_a eq $type_b
		&& (   ( $type_a eq 'JUNC' && $obj_a->get_ID eq $obj_b->get_ID )
		    || ( $type_a eq 'POS'  && $obj_a->get_segID eq $obj_b->get_segID )
		   )
	){
		warn_and_exit "<ERROR>\tpair-wise anchors ($ID) is not valid.\n";
	}

	# record related segment ID
	for my $side (qw/ a b /){
		my $type = $pwanchor_OB->{type}->{$side};
		my $obj  = $pwanchor_OB->{obj}->{$side};
		if( $type eq 'POS' ){
			$pwanchor_OB->{rel_segID}->{$obj->get_segID}->{$side} = 1;
		}
		elsif( $type eq 'JUNC' ){
			$pwanchor_OB->{rel_segID}->{$_}->{$side} = 1 for @{$obj->get_rel_segID_Aref};
		}
	}

	# get groupOB from objects
	my $groupOB_a = $obj_a->get_groupOB;
	my $groupOB_b = $obj_b->get_groupOB;
	if( defined $groupOB_a xor defined $groupOB_b ){
		$pwanchor_OB->{groupOB} = $groupOB_a || $groupOB_b;
	}
	elsif( defined $groupOB_a && defined $groupOB_b ){
		my $groupID_a = $groupOB_a->get_ID;
		my $groupID_b = $groupOB_b->get_ID;
		if( $groupID_a ne $groupID_b ){
			warn_and_exit "<ERROR>\tGroupOB of pair-wise anchors ($ID) is not one ($groupID_a, $groupID_b).\n";
		}
		else{
			$pwanchor_OB->{groupOB} = $groupOB_a;
		}
	}
	else{
		warn_and_exit "<ERROR>\tGroupOB of pair-wise anchors ($ID) doesn't exist ($info_a, $info_b).\n".
							 "\tNote that, at least one of the pair-wise anchors should be relevant to host segment (H).\n";
	}

	bless($pwanchor_OB);
	return $pwanchor_OB;
}

#--- return pwanchor ID
sub get_ID{
	my $pwanchor_OB = shift;
	return $pwanchor_OB->{ID};
}

#--- return group OB
sub get_groupOB{
	my $pwanchor_OB = shift;
	return $pwanchor_OB->{groupOB};
}

#--- return pwanchor linkCount
sub get_linkCount{
	my $pwanchor_OB = shift;
	return $pwanchor_OB->{linkCount};
}

#--- return Info of given side
sub get_sideInfo{
	my $pwanchor_OB = shift;
	my %parm = @_;

	if( $parm{side} eq 'a' ){
		return $pwanchor_OB->{info}->{a};
	}
	elsif( $parm{side} eq 'b' ){
		return $pwanchor_OB->{info}->{b};
	}
	else{
		warn_and_exit "<ERROR>\tWrong side parm ($parm{side}) in pwanchor_OB->get_sideInfo func.\n";
	}
}

#--- test whether relates with given segID
sub has_rel_segID{
	my $pwanchor_OB = shift;
	my %parm = @_;
	my $segID = $parm{segID};
	return exists $pwanchor_OB->{rel_segID}->{$segID};
}

#--- return idx of sides 'a' and 'b' in given segStr_array
sub find_anchor_idx{
	my $pwanchor_OB = shift;
	my %parm = @_;
	my $segStr_Aref = $parm{segStr_Aref};
	my $seg_LastIdx = $parm{seg_LastIdx};

	my $groupOB = $pwanchor_OB->get_groupOB;

	my $accuLen = 0;
	my @anchorLoc;
	my %Ct = ( a => 0, b => 0 );
	# get locations of 'a' and 'b' anchor
	for my $idx ( 0 .. $seg_LastIdx ){
		# my ($segID, $segOrit) = ($segStr_Aref->[$idx] =~ /^([HV]\d+)([\+\-])$/);
		my $segID = $segStr_Aref->[$idx]->{segID};
		my $segLen = $groupOB->segID2segOB(segID=>$segID)->get_length;
		# if has this seg, note that might not be 'a' or 'b'
		if( exists $pwanchor_OB->{rel_segID}->{$segID} ){
			for my $side ( keys %{$pwanchor_OB->{rel_segID}->{$segID}} ){
				my $OB = $pwanchor_OB->{obj}->{$side};
				if(    $pwanchor_OB->{type}->{$side} eq 'POS'
					&& $OB->get_segID eq $segID
				){
					my $edgeDist = $OB->get_edgeDistance( p => ( $segStr_Aref->[$idx]->{segOrit} eq '+' ? 5 : 3 ) );
					# push @anchorLoc, { idx => $idx, side => $side, type => 'POS', obj => $OB, accuLen => $accuLen+$edgeDist };
					push @anchorLoc, { side => $side, accuLen => $accuLen+$edgeDist };
					$Ct{$side}++;
				}
				elsif(   $pwanchor_OB->{type}->{$side} eq 'JUNC'
					  && $OB->has_rel_segID( segID => $segID )
					  && $idx != $seg_LastIdx # JUNC till the last 2nd
					  # && $OB->is_same( SegLink => $segStr_Aref->[$idx].','.$segStr_Aref->[$idx+1] )
					  && $OB->is_same( SegLink => $segStr_Aref->[$idx]->{info}.','.$segStr_Aref->[$idx+1]->{info} )
				){
					# push @anchorLoc, { idx => $idx, side => $side, type => 'JUNC', obj => $OB, accuLen => $accuLen+$segLen };
					push @anchorLoc, { side => $side, accuLen => $accuLen+$segLen };
					$Ct{$side}++;
				}
			}
		}
		# update
		$accuLen += $segLen;
	}

	# return anchors' LocInfo and valid sign (i.e., found both anchors)
	return ( \@anchorLoc, $Ct{a}*$Ct{b} );
}

#--- calculate distance between two given anchors on given segStr
# note that segStr START and END could pretend as JUNC from -1 or $#segStr, respectively
sub get_anchor_distance{
	my $pwanchor_OB = shift;
	my %parm = @_;
	my $segStr = $parm{segStr};
	my $anchor_a_Href = $parm{anchor_a_Href};
	my $anchor_b_Href = $parm{anchor_b_Href};

	my $groupOB = $pwanchor_OB->get_groupOB;
	my @segStr = split /;/, $segStr;
	my @encloseSegLen;
	# sort anchor by idx
	my ($fw_anchor_Href, $bh_anchor_HRef) = sort {$a->{idx} <=> $b->{idx}} ($anchor_a_Href, $anchor_b_Href);
	my $fw_idx = $fw_anchor_Href->{idx};
	my $bh_idx = $bh_anchor_HRef->{idx};
	# forward POS anchor
	if( $fw_anchor_Href->{type} eq 'POS' ){
		my $posOB = $fw_anchor_Href->{obj};
		push @encloseSegLen, $posOB->get_edgeDistance( p => ($segStr[$fw_idx] =~ /\+/ ? 3 : 5) );
	}
	# inner segment, and behind JUNC anchor
	my @in_idx = ( $fw_idx+1 .. $bh_idx-1 );
	push @in_idx, $bh_idx if( $bh_anchor_HRef->{type} eq 'JUNC' && $fw_idx != $bh_idx );
	for my $in_idx ( @in_idx ){
		my ($segID) = ($segStr[$in_idx] =~ /^([HV]\d+)[\+\-]/);
		my $segOB = $groupOB->segID2segOB( segID => $segID );
		push @encloseSegLen, $segOB->get_length;
	}
	# behind POS anchor
	if( $bh_anchor_HRef->{type} eq 'POS' ){
		my $posOB = $bh_anchor_HRef->{obj};
		# $fw_idx != $bh_idx, otherwise, never happen
		push @encloseSegLen, $posOB->get_edgeDistance( p => ($segStr[$fw_idx] =~ /\+/ ? 5 : 3) );
	}

	return sum(@encloseSegLen);
}

1; ## tell the perl script the successful access of this module.
