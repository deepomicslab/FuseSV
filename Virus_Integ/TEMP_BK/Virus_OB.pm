package FuseSV::Virus_Integ::Virus_OB;

use strict;
use warnings;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr/;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::Reads_OB;
use FuseSV::Virus_Integ::ReadsGroup_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Virus_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-01-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
					 /;

#--- construction of object ---
sub new{
	my $type = shift;
	my %parm = @_;

	my $virus_OB = {};
	$virus_OB->{SegName} = $parm{SegName};
	$virus_OB->{IdxFa}   = $parm{IdxFa};
	$virus_OB->{subType} = $parm{subType};
	$virus_OB->{extLen}  = $parm{extLen};

	$virus_OB->{seq} = '';

	bless($virus_OB);
	return $virus_OB;
}

1; ## tell the perl script the successful access of this module.
