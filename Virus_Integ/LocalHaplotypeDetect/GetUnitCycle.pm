package FuseSV::Virus_Integ::LocalHaplotypeDetect::GetUnitCycle;

use strict;
use warnings;
use Getopt::Long;
use File::Spec qw/ catfile /;
use FindBin qw/ $RealBin /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  get_unit_cycle
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'GetUnitCycle';
#----- version --------
$VERSION = "0.02";
$DATE = '2020-09-13';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						get_unit_cycle
					 /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} get_UCYC <[Options]>
	 
	 Options:

	    -lhs [s]  path of the 'LocalHaplotypeSolver' folder. <required>
	              Download from $V_Href->{weblink}

	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
 # return "<ERROR>\tcannot find external script LocalHaplotypeSolver/LocalHap.py in $MODULE_NAME module.\n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$V_Href->{ $_->[0] } = $_->[1] for
		map {
			if( !exists $V_Href->{$_->[0]} ){
				( $_ );
			}
			else{
				warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
			}
		}
		(
			#--- input/output
			[ LHS_folder => 'null' ], # LocalHaplotypeSolver folder
			# external script to run
			[ script => undef ],
			[ weblink => 'https://github.com/deepomicslab/LocalHaplotypeSolver' ],
		);
}

#--- get options from command line ---
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	# GetOptions(
	# 	# input/output
	# 	"-lhs:s"	=> \$V_Href->{LHS_folder},
	# );
	for (my $i=0; $i<@{$V_Href->{argv_Aref}}; $i++){
		if($V_Href->{argv_Aref}->[$i] eq '-lhs'){
			$V_Href->{LHS_folder} = $V_Href->{argv_Aref}->[$i+1];
			splice @{$V_Href->{argv_Aref}}, $i, 2;
			last;
		}
	}
}

#--- test para and alert ---
sub para_alert{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	$V_Href->{script} = File::Spec->catfile($V_Href->{LHS_folder},'LocalHap.py');
	return  (   !file_exist( filePath => $V_Href->{LHS_folder} )
			 || !file_exist( filePath => $V_Href->{script} )
			);
}

#--- get unit-cycle ---
sub get_unit_cycle{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# run external script
	stout_and_sterr "[INFO]\tthis func encapsulates 'LocalHap.py', just use its options for '$V_Href->{command}'.\n"
						 ."\te.g., perl $V_Href->{MainName} $V_Href->{command} -h\n";
	system("python $V_Href->{script} @{$V_Href->{argv_Aref}}");
}

#--- 
1; ## tell the perl script the successful access of this module.
