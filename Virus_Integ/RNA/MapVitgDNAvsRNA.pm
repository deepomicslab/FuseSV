package FuseSV::Virus_Integ::RNA::MapVitgDNAvsRNA;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr trible_run_for_success Get_Two_Seg_Olen /;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  map_virus_integration_DNA_vs_RNA
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'MapVitgDNAvsRNA';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-02-19';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						map_virus_integration_DNA_vs_RNA
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} vitg_DvsR <[Options]>
	 
	 Options:

	   # Inputs and Outputs #


	    -h|help      Display this help info.

	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
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
			# input/output
			[ dHG_bam => {} ],
			## dHG_draw_region -> {}, later will change it to hash
			## dHG_draw_region -> refseg = $refseg
			## dHG_draw_region -> st_pos = $st_pos
			## dHG_draw_region -> ed_pos = $ed_pos
			[ dHG_draw_region => undef ],
			[ dHG_reg_resol_file => undef ],
			[ dG_figure => '' ],
			# database
			[ trans_psl_idx_bgz => undef ],
			[ rmsk_idx_bgz => undef ],
			[ dHG_rmsk_type => [] ], # later change to Hash
			[ dHG_show_gene => [] ], # later change to Hash
			[ dG_user_gene => [] ],
			# software
			[ samtools => '' ],
			[ tabix => '' ],
			# sample info
			[ dHG_meandepth => {case=>0, ctrl=>0} ],
			[ dHG_casePurity => 1 ],
			[ dHG_casePloidy => 2 ],
			[ dHG_ctrlPloidy => 2 ],
			[ dHG_break_pos => [] ], # later change to Hash
			[ dHG_sampleIDtoShow => '' ],
			# bam -> 'samtools depth'
			[ min_baseQ => 5 ],
			[ min_mapQ => 10 ],
			[ no_SoftClip => 0 ],
			# depth options
			[ dG_DepTrimRatio => 0 ],
			[ dG_DepEng_3SD => 0 ],
			[ dG_DepSmooth_flkWd => 2 ],
			[ dHG_ctrl_depthFluctR => 0 ],
			[ dHG_ctrl_prior_copyR => [] ],
			# dHG_CNbaseOnReg -> tissue -> $tissue
			# dHG_CNbaseOnReg -> refseg -> $refseg
			# dHG_CNbaseOnReg -> st_pos -> $st_pos
			# dHG_CNbaseOnReg -> ed_pos -> $ed_pos
			# dHG_CNbaseOnReg -> st_winNO -> $st_winNO
			# dHG_CNbaseOnReg -> ed_winNO -> $ed_winNO
			# dHG_CNbaseOnReg -> cn -> $cn
			# dHG_CNbaseOnReg -> mean_adjdepth -> $mean_adjdepth 
			# dHG_CNbaseOnReg -> depth2cn_ratio -> $depth2cn_ratio
			[ dHG_CNbaseOnReg => undef ],
			[ dHG_CNbaseOnReg_minDepthLimit => 5], # empirical
			[ dHG_TissueToShowCN => 'case' ],
			[ dHG_CNlineColor => 'blue' ],
			[ dHG_CNlineWidth => 1 ],
			# draw options
			[ dG_autoAdjustResol => 0 ],
			[ dG_DepthSpecColor => { case=>'red', ctrl=>'dodgerblue'} ],
			[ dG_extWidth => 150 ],
			[ dG_extHeight => 150 ],
			[ dG_axisZX => 0 ],
			[ dG_axisZY => 0 ],
			[ dG_x_axisLen => 0 ], # auto assign
			[ dG_x_LabBaseNum => 1000 ],
			[ dG_x_LabStep => 200 ],
			[ dG_x_LabFtsz => 12 ],
			[ dG_x_RegResolColor => undef ],
			[ dG_y_axisLen => 80 ],
			[ dG_y_resol => 0 ], # auto adjust
			[ dG_y_LabStep => 0 ], # auto adjust
			[ dG_y_LabFtsz => 12 ],
			[ dHG_minEleWidth => 2 ],
			[ dHG_x_RegResolColor => undef ],
			[ dHG_x_RegResolColorMinOpacity => 0.3 ],
			[ dHG_SegmentMode => undef ],
			[ dHG_SegmFillCol => 'lightblue' ],
			#--- intermediate variants
			[ DepthToDraw => undef ],
			[ dG_CaseVsCtrl_mode => undef ],
			[ case_Tcell_GMratio => undef ],
			[ case_Ncell_GMratio => undef ],
			# [ output_pref => undef ],
			[ CtrlPos2Depth => {} ],
			# check regResol_Href structure in FuseSV::Visual::ArrangeObj module
			[ dG_regResol => {} ],
			# check winDepth_Href structure in FuseSV::Visual::DepthWork module
			[ dG_winDepth => {} ],
			# dG_bkWindow -> $window_NO = $bkpos
			[ dG_bkWindow => {} ],
			# check TransInfo_Href structure in FuseSV::Visual::GeneWork module
			[ dG_TransInfo => {} ],
			# check RmskInfo_Href structure in FuseSV::Visual::RMSKwork module
			[ dG_RmskInfo => {} ],
			[ dG_maxDepthToDraw => 0 ],
			[ dG_maxWinNOtoDraw => 0 ],
			[ dG_bgStrokeCol => 'black' ],
			[ dG_DepthSpecBtom_Y => 0 ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['dG_figure'],
								  ['dHG_bam', 'case', 'filepath'],
								  ['dHG_bam', 'ctrl', 'filepath'],
								  ['dHG_reg_resol_file'],
								  ['trans_psl_idx_bgz'],
								  ['rmsk_idx_bgz'],
								  ['MH_seq_file'],
								  ['tabix'],
								  ['samtools']                    ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-o:s"		=> \$V_Href->{dG_figure},
		# help
		"-h|help"	=> \$V_Href->{HELP},
		# for debug
		"-debug"	=> \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 # || $V_Href->{dG_figure} !~ /\.([sp][vnd][gf])$/i
			);
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub map_virus_integration_DNA_vs_RNA{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	#--------------#
	# prepare work #
	#--------------#
	# &prepare;
	warn_and_exit "real_map_virus_integration_DNA_vs_RNA is in development.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
