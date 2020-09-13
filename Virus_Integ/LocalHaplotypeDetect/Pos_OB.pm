package FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB;

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

$MODULE_NAME = 'Pos_OB';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-05-05';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						get_groupOB
						get_edgeDistance
						get_segID
						is_valid
					 /;

#--- structure of object
# pos_OB -> chr = $chr
# pos_OB -> pos = $pos
# pos_OB -> segOB = $segOB
# pos_OB -> segID = $segID
# pos_OB -> edgeD = { 5->$DistToSegP5, 3->$DistToSegP3 }
# pos_OB -> valid = 1/0

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $chr = $parm{chr};
	my $pos = $parm{pos};
	my $seg_OB = $parm{seg_OB};
	my $minEdgeD = $parm{minEdgeD} || 100;
	my $keepError = $parm{keepError} || 0;

	# check chr
	if( $chr ne $seg_OB->get_chr ){
		warn_and_exit "<ERROR>\tPos_OB chr ($chr) diff from its seg_OB.\n";
	}
	# check pos
	if( !$seg_OB->test_pos_in( testPos => $pos ) ){
		warn_and_exit "<ERROR>\tPos_OB pos ($pos) is out range of its seg_OB.\n";
	}

	my $pos_OB = {};

	$pos_OB->{valid} = 1;
	$pos_OB->{chr} = $chr;
	$pos_OB->{pos} = $pos;
	$pos_OB->{segOB} = $seg_OB;
	$pos_OB->{segID} = $seg_OB->get_ID;
	$pos_OB->{edgeD}->{5} = $pos - $seg_OB->get_stpos;
	$pos_OB->{edgeD}->{3} = $seg_OB->get_edpos - $pos;

	# check edge distance
	if(    $pos_OB->{edgeD}->{5} < $minEdgeD
		|| $pos_OB->{edgeD}->{3} < $minEdgeD
	){
		if( !$keepError ){
			warn_and_exit "<ERROR>\tPos_OB pos ($pos) is too close to edges of its segment (".$seg_OB->get_ID.").\n".
								 "\tminimum distance is $minEdgeD; this pos: 5p is $pos_OB->{edgeD}->{5}, 3p is $pos_OB->{edgeD}->{3}.\n";
		}
		else{
			$pos_OB->{valid} = 0;
		}
	}

	bless($pos_OB);
	return $pos_OB;
}

#--- return group OB
sub get_groupOB{
	my $pos_OB = shift;
	return $pos_OB->{segOB}->get_groupOB;
}

#--- return distance to segOB's given prime edge
sub get_edgeDistance{
	my $pos_OB = shift;
	my %parm = @_;

	if( $parm{p} =~ /^5/ ){
		return $pos_OB->{edgeD}->{5};
	}
	elsif( $parm{p} =~ /^3/ ){
		return $pos_OB->{edgeD}->{3};
	}
	else{
		warn_and_exit "<ERROR>\tWrong p parm ($parm{p}) in Pos_OB->get_edgeDistance func.\n";
	}
}

#--- return segment ID it locotes
sub get_segID{
	my $pos_OB = shift;
	return $pos_OB->{segID};
}

#--- return valid sign
sub is_valid{
	my $pos_OB = shift;
	return $pos_OB->{valid};
}

1; ## tell the perl script the successful access of this module.
