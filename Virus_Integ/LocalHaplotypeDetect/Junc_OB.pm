package FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB;

use strict;
use warnings;
use SOAPfuse::General_Operation qw/ warn_and_exit Get_Two_Seg_Olen /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Junc_OB';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-05-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						load_segOB
						get_ID
						get_segID
						get_segOB
						get_groupOB
						get_rel_segID_Aref
						has_rel_segID
						is_same
					 /;

#--- structure of object
# Junc_OB -> ID = $ID
# Junc_OB -> SegLink = $SegLink : e.g., V1+,H4-
# Junc_OB -> SegLink_rvcp = $SegLink_rvcp : H4+,V1-
# Junc_OB -> p5_segID = [HV]\d+
# Junc_OB -> p5_segOrit = [+-]
# Junc_OB -> p3_segID = [HV]\d+
# Junc_OB -> p3_segOrit = [+-]
# Junc_OB -> p5_segOB = $p5_segOB/free(virus)
# Junc_OB -> p3_segOB = $p3_segOB/free(virus)
# Junc_OB -> rel_segID = { $segID -> 1, .. }
# Junc_OB -> groupOB = $groupOB : may be undef, when it is V-vs-V

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $ID = $parm{ID};
	my $SegLink = $parm{SegLink};
	my $p5_segOB = $parm{p5_segOB};
	my $p3_segOB = $parm{p3_segOB};
	my $p5_segOrit = $parm{p5_segOrit};
	my $p3_segOrit = $parm{p3_segOrit};

	# check ID
	if( $ID !~ /^JUNC\d+$/ ){
		warn_and_exit "<ERROR>\tWrong Segment ID ($ID), cannot create segment object.\n"
							."\tJunc ID must has prefix 'JUNC' followed by digits.\n";
	}

	my $junc_OB = {};

	$junc_OB->{ID} = $ID;
	$junc_OB->{SegLink} = $SegLink;
	$junc_OB->{p5_segOB} = $p5_segOB;
	$junc_OB->{p3_segOB} = $p3_segOB;
	$junc_OB->{p5_segID} = $p5_segOB->get_ID;
	$junc_OB->{p3_segID} = $p3_segOB->get_ID;
	$junc_OB->{p5_segOrit} = $p5_segOrit;
	$junc_OB->{p3_segOrit} = $p3_segOrit;
	$junc_OB->{rel_segID}->{$junc_OB->{$_}} = 1 for qw/ p5_segID p3_segID /;

	# reversed complementary junc
	(my $rvcp_p5_segOrit = $p3_segOrit) =~ tr/\+\-/\-\+/;
	(my $rvcp_p3_segOrit = $p5_segOrit) =~ tr/\+\-/\-\+/;
	$junc_OB->{SegLink_rvcp} = $junc_OB->{p3_segID}.$rvcp_p5_segOrit.','.$junc_OB->{p5_segID}.$rvcp_p3_segOrit;

	# assign groupOB according to segOB
	my $p5_seg_groupOB = $p5_segOB->get_groupOB;
	my $p3_seg_groupOB = $p3_segOB->get_groupOB;
	if( defined $p5_seg_groupOB xor defined $p3_seg_groupOB ){
		$junc_OB->{groupOB} = $p5_seg_groupOB || $p3_seg_groupOB;
	}
	elsif( defined $p5_seg_groupOB && defined $p3_seg_groupOB ){
		my $p5_seg_groupID = $p5_seg_groupOB->get_ID;
		my $p3_seg_groupID = $p3_seg_groupOB->get_ID;
		if( $p5_seg_groupID ne $p3_seg_groupID ){
			warn_and_exit "<ERROR>\tjunc_OB $junc_OB->{ID} p5_seg_groupID ($p5_seg_groupID) diffs p3_seg_groupID ($p3_seg_groupID).\n";
		}
		else{
			$junc_OB->{groupOB} = $p5_seg_groupOB;
		}
	}

	bless($junc_OB);
	return $junc_OB;
}

#--- load segOB of given prime
sub load_segOB{
	my $junc_OB = shift;
	my %parm = @_;

	if( $parm{p} =~ /^5/ ){
		$junc_OB->{p5_segOB} = $parm{seg_OB};
	}
	elsif( $parm{p} =~ /^3/ ){
		$junc_OB->{p3_segOB} = $parm{seg_OB};
	}
	else{
		warn_and_exit "<ERROR>\tWrong p parm ($parm{p}) in junc_OB->load_segOB func.\n";
	}
}

#--- return segment ID
sub get_ID{
	my $junc_OB = shift;
	return $junc_OB->{ID};
}

#--- return segID of given prime
sub get_segID{
	my $junc_OB = shift;
	my %parm = @_;

	if( $parm{p} =~ /^5/ ){
		return $junc_OB->{p5_segID};
	}
	elsif( $parm{p} =~ /^3/ ){
		return $junc_OB->{p3_segID};
	}
	else{
		warn_and_exit "<ERROR>\tWrong p parm ($parm{p}) in junc_OB->get_segID func.\n";
	}
}

#--- return segOB of given prime
sub get_segOB{
	my $junc_OB = shift;
	my %parm = @_;

	if( $parm{p} =~ /^5/ ){
		return $junc_OB->{p5_segOB};
	}
	elsif( $parm{p} =~ /^3/ ){
		return $junc_OB->{p3_segOB};
	}
	else{
		warn_and_exit "<ERROR>\tWrong p parm ($parm{p}) in junc_OB->get_segOB func.\n";
	}
}

#--- return groupOB
sub get_groupOB{
	my $junc_OB = shift;
	return $junc_OB->{groupOB};
}

#--- return rel_segID array ref
sub get_rel_segID_Aref{
	my $junc_OB = shift;
	return [ keys %{$junc_OB->{rel_segID}} ];
}

#--- test whether relates with given segID
sub has_rel_segID{
	my $junc_OB = shift;
	my %parm = @_;
	my $segID = $parm{segID};
	return exists $junc_OB->{rel_segID}->{$segID};
}

#--- test whether given junc is this juncOB
sub is_same{
	my $junc_OB = shift;
	my %parm = @_;
	my $SegLink = $parm{SegLink}; # e.g., V1+,H4- or H4+,V1-

	if(    $SegLink eq $junc_OB->{SegLink}
		|| $SegLink eq $junc_OB->{SegLink_rvcp}
	){
		return 1;
	}
	else{
		return 0;
	}
}

1; ## tell the perl script the successful access of this module.
