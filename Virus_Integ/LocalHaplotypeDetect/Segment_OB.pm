package FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB;

use strict;
use warnings;
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

$MODULE_NAME = 'Segment_OB';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-04-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						get_chr
						get_stpos
						get_edpos
						get_length
						get_ID
						get_SpecTag
						get_segNO
						get_groupOB
						test_pos_in
						test_sequential_seg
						test_stream_seg
						test_overlap_seg
					 /;

#--- structure of object
# segment_OB -> ID = $ID
# segment_OB -> SpecTag = $SpecTag # 'H' or 'V'
# segment_OB -> SegNO = $SegNO
# segment_OB -> region = $region
# segment_OB -> chr = $chr
# segment_OB -> stpos = $stpos
# segment_OB -> edpos = $edpos
# segment_OB -> length = $length
# segment_OB -> groupOB = $groupOB
# segment_OB -> orig_cn = $orig_cn
# segment_OB -> cn_ref = $cn_ref, this is the selection of major/minor as normal-allele
# segment_OB -> cn_alt = $cn_alt
# segment_OB -> cn_alt_LP = $cn_alt_LP
# segment_OB -> LP_weight = $LP_weight
# segment_OB -> LP_offset = $LP_offset

#--- construction of object
sub new{
	my $type = shift;
	my %parm = @_;
	my $ID = $parm{ID};
	my $groupOB = $parm{groupOB};
	my $region = $parm{region};
	my $orig_cn = $parm{orig_cn};
	my $cn_ref = $parm{cn_ref}; # for host-segment
	my $cn_alt = $parm{cn_alt};
	my $cn_alt_LP = $parm{cn_alt_LP};
	my $LP_weight = $parm{LP_weight};
	my $LP_offset = $parm{LP_offset};

	# check ID
	if( $ID !~ /^[HV]\d+$/ ){
		warn_and_exit "<ERROR>\tWrong Segment ID ($ID), cannot create segment object.\n"
							."\tSegment ID must has prefix 'H' for host OR 'V' for virus, followed by digits.\n";
	}

	my $segment_OB = {};

	$segment_OB->{ID} = $ID;
	($segment_OB->{SpecTag}, $segment_OB->{SegNO}) = ($ID =~ /^(.*\D)(\d+)$/);
	$segment_OB->{region} = $region;
	$segment_OB->{groupOB} = $groupOB if( defined $groupOB );
	$segment_OB->{orig_cn} = $orig_cn if( defined $orig_cn );
	$segment_OB->{cn_ref} = $cn_ref if( defined $cn_ref );
	$segment_OB->{cn_alt} = $cn_alt;
	$segment_OB->{cn_alt_LP} = $cn_alt_LP if( defined $cn_alt_LP );
	$segment_OB->{LP_weight} = $LP_weight if( defined $LP_weight );
	$segment_OB->{LP_offset} = $LP_offset if( defined $LP_offset );

	if( $region =~ /^(.+):(\d+)-(\d+)$/ ){
		$segment_OB->{chr} = $1;
		$segment_OB->{stpos} = $2;
		$segment_OB->{edpos} = $3;
		$segment_OB->{length} = $3 - $2 + 1;
	}
	else{
		warn_and_exit "<ERROR>\tCannot analyze segment region: $region\n";
	}

	bless($segment_OB);
	return $segment_OB;
}

#--- return segment chr
sub get_chr{
	my $segment_OB = shift;
	return $segment_OB->{chr};
}

#--- return segment stpos
sub get_stpos{
	my $segment_OB = shift;
	return $segment_OB->{stpos};
}

#--- return segment edpos
sub get_edpos{
	my $segment_OB = shift;
	return $segment_OB->{edpos};
}

#--- return segment length
sub get_length{
	my $segment_OB = shift;
	return $segment_OB->{length};
}

#--- return segment ID
sub get_ID{
	my $segment_OB = shift;
	return $segment_OB->{ID};
}

#--- return species tag
sub get_SpecTag{
	my $segment_OB = shift;
	return $segment_OB->{SpecTag};
}

#--- return segment NO
sub get_segNO{
	my $segment_OB = shift;
	return $segment_OB->{SegNO};
}

#--- return group OB
sub get_groupOB{
	my $segment_OB = shift;
	return $segment_OB->{groupOB};
}

#--- test whether given pos is in this segment
sub test_pos_in{
	my $segment_OB = shift;
	my %parm = @_;
	my $testPos = $parm{testPos};

	if(    $testPos >= $segment_OB->{stpos}
		&& $testPos <= $segment_OB->{edpos}
	){
		return 1;
	}
	else{
		return 0;
	}
}

#--- test whether be sequential with given segment
sub test_sequential_seg{
	my $segment_OB = shift;
	my %parm = @_;
	my $testSegOB = $parm{testSegOB};
	my $testSegPrime = $parm{testSegPrime} || '3p';

	return 0 if( $segment_OB->{chr} ne $testSegOB->{chr} );

	# 3p: 5p -> this|test -> 3p
	# 5p: 5p -> test|this -> 3p
	if(    ( $testSegPrime eq '3p' && $segment_OB->{edpos} == $testSegOB->{stpos} - 1 )
		|| ( $testSegPrime eq '5p' && $segment_OB->{stpos} == $testSegOB->{edpos} + 1 )
	){
		return 1;
	}
	else{
		return 0;
	}
}

#--- test whether given segment locates given stream
sub test_stream_seg{
	my $segment_OB = shift;
	my %parm = @_;
	my $testSegOB = $parm{testSegOB};
	my $testSegPrime = $parm{testSegPrime} || '3p';

	return 0 if( $segment_OB->{chr} ne $testSegOB->{chr} );

	# 3p: 5p -> this|---|test -> 3p
	# 5p: 5p -> test|---|this -> 3p
	if(    ( $testSegPrime eq '3p' && $segment_OB->{edpos} <= $testSegOB->{stpos} - 1 )
		|| ( $testSegPrime eq '5p' && $segment_OB->{stpos} >= $testSegOB->{edpos} + 1 )
	){
		return 1;
	}
	else{
		return 0;
	}
}

#--- test whether overlaps with given segment ---
sub test_overlap_seg{
	my $segment_OB = shift;
	my %parm = @_;
	my $testSegOB = $parm{testSegOB};

	return Get_Two_Seg_Olen( $segment_OB->{stpos}, $segment_OB->{edpos},
							  $testSegOB->{stpos},  $testSegOB->{edpos}  );
}

1; ## tell the perl script the successful access of this module.
