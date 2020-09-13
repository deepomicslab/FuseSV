package FuseSV::Virus_Type::InitialMapSmartFQonVirus;

use strict;
use warnings;
use File::Spec qw/catfile/;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use FuseSV::LoadOn;
use FuseSV::Virus_Type::SmartAlignment;
use Parallel::ForkManager;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  initial_align_smfq_on_all_vref
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'InitialMapSmartFQonVirus';
#----- version --------
$VERSION = "0.04";
$DATE = '2017-05-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						initial_align_smfq_on_all_vref
					 /;

#--- initial alignment of smartFQ to virus ref ---
sub initial_align_smfq_on_all_vref{

	# variants
	my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
	my $v_initial_align_folder = $V_Href->{Init_align_folder};
	my $Fork_Number = $V_Href->{fork_number};

	# initial alignment folder
	`rm -rf $v_initial_align_folder`;
	`mkdir -p $v_initial_align_folder`;

	# potential multiple fork operations, initialization
	my $multiple_fork;
	if($Fork_Number > 1){
		$multiple_fork = new Parallel::ForkManager($Fork_Number);
		# inform
		stout_and_sterr `date`."[INFO]:\tOperate in multiple forks mode. Fork Number: $Fork_Number.\n";
	}

	# loop each virus
	my $OutBam_FilePref = "SmartFq.MapToVrefNO";
	my @sorted_v_refseg = sort {$Virus_Info_Href->{$a}->{vref_no} <=> $Virus_Info_Href->{$b}->{vref_no}} keys %$Virus_Info_Href;
	for my $v_refseg ( @sorted_v_refseg ){
		# potential multiple fork operations, child fork starts
		if($Fork_Number > 1){
			$multiple_fork->start and next;
		}

		# real operation
		my $v_refseg_ref_fa = $Virus_Info_Href->{$v_refseg}->{vref_fa};
		my $v_ref_NO = $Virus_Info_Href->{$v_refseg}->{vref_no};
		my $OutBam_pref = File::Spec->catfile($v_initial_align_folder, $OutBam_FilePref.$v_ref_NO);



		# align smart fq of all RGs against this v_ref
		my $init_mapbam = FuseSV::Virus_Type::SmartAlignment->SmartFq2MergeBam_against_given_ref(
			  	{
					RefForBwa	=>	$v_refseg_ref_fa,
					OutBamPref	=>	$OutBam_pref
				}
			);



		# if empty, warn it
		if( $init_mapbam eq 'EMPTY' ){
			stout_and_sterr `date`."[WARN]:\tNo read maps to virus reference NO.$v_ref_NO. Note as 'EMPTY'\n";
		}
		# store path on the disk
		`echo $init_mapbam > $OutBam_pref.mkdup.bam.list`;

		# potential multiple fork operations, child fork finishs
		if($Fork_Number > 1){
			$multiple_fork->finish;
		}
	}

	# potential multiple fork operations, wait all childern forks
	if($Fork_Number > 1){
		$multiple_fork->wait_all_children;
	}

	# retrieve the OutBam path
	for my $v_refseg ( @sorted_v_refseg ){
		my $v_ref_NO = $Virus_Info_Href->{$v_refseg}->{vref_no};
		# read from disk
		my $OutBam_list = File::Spec->catfile($v_initial_align_folder, $OutBam_FilePref.$v_ref_NO.'.mkdup.bam.list');
		chomp($Virus_Info_Href->{$v_refseg}->{init_mapbam} = `cat $OutBam_list`);
		# sweep
		`rm $OutBam_list`;
	}

	# inform
	stout_and_sterr `date`."[INFO]:\tInitially map the candidate smartPE reads to all virus reference.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
