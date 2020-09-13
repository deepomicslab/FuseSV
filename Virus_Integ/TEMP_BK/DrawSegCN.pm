package FuseSV::Virus_Integ::DrawSegCN;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr trible_run_for_success Get_Two_Seg_Olen file_exist /;
use SOAPfuse::SVG_Orthogonal_System_Elements qw/ draw_a_parallelogram /;
use SOAPfuse::SVG_Font qw/ show_text_in_line get_size_of_text_to_show /;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::Bam_OB;
use FuseSV::Visual::ArrangeObj qw/ load_regResolForWindowlize get_pos_windowNO allocate_object_bilater allocate_object_vertical /;
use FuseSV::Visual::DepthWork qw/ get_windowSmoDepth deal_window_depth_info get_ctrl_copyR draw_depth_spectrum get_givenWinItvalMeanDepth /;
use FuseSV::Visual::GeneWork qw/ load_transInfo_from_tpsl add_user_region_to_gene_info allocate_genes_and_draw /;
use FuseSV::Visual::RMSKwork qw/ load_rmskInfo_from_idxBgz allocate_RMSKele_and_draw /;
use FuseSV::Visual::AxisWork qw/ set_Y_axis_Resol_LabStep show_Y_axis show_X_axis /;
use FuseSV::Visual::RegionWork qw/ arrange_region_with_clip /;
use FuseSV::Virus_Integ::ObjCountTransfMixToPureT qw/ get_Tcell_GMpart get_Tcell_GMratio get_Ncell_GMpart get_Ncell_GMratio get_ObjCountOfPureTumorCell get_ObjSingleCNdepthInMixed /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  real_draw_segmemt_CN
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'DrawSegCN';
#----- version --------
$VERSION = "5.33";
$DATE = '2018-04-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						real_draw_segmemt_CN
						prepare
						deal_with_showCNopt
						load_window_adjdepth
						get_window_showdepth
						load_trans_info
						load_rmsk_info
						set_x_axis
						set_y_axis
						draw_rmsk
						draw_y_axis
						draw_depth
						draw_genes
						draw_x_axis
						show_breakpoints
						show_figure_information
						generate_figure_output
						get_HostSegmentID
						read_depthStat_to_window_adjdepth
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} draw_segCN <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -t_bam  [s]  indexed bam file of tumor (case) sample. <[optional]>
	    -t_mdp  [f]  average sequencing depth of tumor (case) sample. <required with '-t_bam'>
	    -n_bam  [s]  indexed bam file of normal (control) sample. <[optional]>
	    -n_mdp  [f]  average sequencing depth of normal (control) sample. <required with '-n_bam'>
	    -reg    [s]  genomic region to draw, e.g., 'chr12:100000-120000'. <required>
	    -rgsol  [s]  resolution config file of each region. <required>
	    -o      [s]  output figure file, allowed format: svg/png/pdf. <required>

	   # Database and Metadata #
	    -tpsl   [s]  transcript PSL annotation file, should be indexed bgz file. [optional]
	    -rmsk   [s]  repeat masker annotation file, should be indexed bgz file. [optional]
	    -eltp   [s]  RMSK element types to show. [all available]
	    -gene   [s]  gene to show. [all available]
	    -nreg   [s]  specific region(s) in '-reg' to display. [none]
	                  instance: -nreg 'TEST1:110000-110200' -nreg 'TEST2:110400-110700,111400-111700'

	   # Software Required #
	    -stl    [s]  SamTools, minimum version: 1.3. <required>
	    -tbi    [s]  Tabix path. <required>

	   # Sample Attributes #
	    -id     [s]  sample ID to show in the figure. [optional]
	    -t_puty [f]  purity of tumor (case) sample. [1]
	    -t_pldy [f]  ploidy of tumor (case) sample. [2]
	    -n_pldy [f]  ploidy of normal (control) sample. [2]
	    -junc   [s]  breakage junctions info file. [optional]
	                 denotes the virus integration position(s) and rearrangement break-point(s) if has.
	    -vseg   [s]  viral segments info. [optional]
	                 to calculate CN of virus segments for local genomic map construction.

	   # Options about Depth #
	    -baseQ  [i]  minimum base quality for 'SamTools depth'. [5]
	    -mapQ   [i]  minimum mapping quality for 'SamTools depth'. [10]
	    -no_sfcp     ignore the Softclip reads for depth. [disabled]
	    -tmrt   [f]  ratio to trim for even depth of one region. [0]
	                 opsitive value: bilateral operation; negative: only tail.
	    -eng3sd      apply engineer_3times_SD_filter for average depth. [disabled]
	                 Once enabled, '-tmrt' will be ignored.
	    -fwsm   [i]  how many bilateral flanking windows to make curve smooth. [2]
	                  note: value '2' means 2(fore)+1(this)+2(after); set '0' to disable.
	    -n_dptr [f]  depth ratio selected from normal tissue for regional depth adjustment. [0]
	                 defaults is zero, is disable this operation.
	    -nprcr  [s]  the abnormal copy ratio in normal (control) sample, multiple inputs ok.
	                 format: region_start_pos-region_end_pos:copyratio

	   # Options about Display #
	   ## colour
	    -T_col  [s]  the colour of depth spectrum from tumor (case) sample. [red]
	    -N_col  [s]  the colour of depth spectrum from normal (control) sample. [dodgerblue]
	    -bgbcol [s]  the colour of background boundary. [black]
	   ## axis
	    -x_lnb  [i]  label number basement shown on x-axis. [1000]
	    -x_lbs  [i]  how many screen pixels to show one position label on x-axis. [200]
	    -x_lbfz [i]  font size of label on x-axis. [12]
	    -x_rrc  [s]  color to show gradient region resolution. [none]
	    -y_len  [i]  length of y-axis, will be overwritten when '-y_res' is set. [80]
	    -y_res  [f]  how many depth (X) one spot in figure represents on y-axis. [auto]
	    -y_lbs  [i]  how many screen pixels to show one depth marker on y-axis, effective with '-y_res'. [auto]
	    -y_lbfz [i]  font size of label on y-axis. [12]
	   ## segmentation * effective with available '-bkpos' input *
	    -segm   [s]  show segmentation of this region. [u]
	                 note: 1) use 'u' for Upper-case (A-Z) naming; 'l' for Lower-case (a-z);
	                       2) other inputs will act as prefix postfixed by digital.
	                       3) first index can determined as inputs like 'u:C', 'l:d', or 'blabla:3'.
	                       4) blank (regex detection: \\s) is not allowed.
	    -segc   [s]  segmentation bg-color. ['lightblue']
	    -cn_bs  [s]  calculate copy number (CN) of each segments according to this region. [optional]
	                 format: 'T:chr12:105000-110000:2'
	                 note: 1) otherwise, use 'HapChrDepth mode' to calculate CN when tumor sample is provided.
	                       2) set this option if you want to calculate CN when only normal sample is provided.
	                       3) 'T' is tissue type, you could also use 'N' when '-n_bam' is given.
	                       4) the last number is copy number statement of this region.
	                       5) show CN of 'case' once '-t_bam' is set, otherwise show CN of 'control'.
	                       6) CN information will auto output to the file name 'prefix(para(-o)).SegCN_info.tsv'.
	   ## others
	    -auto_r [i]  auto adjusted resolution to make genomic region shown in such size of pixels. [disabled]
	    -miew   [f]  the minimum width of elements (exon/CDS/RMSKele) in the figure. [1]
	    -extwd  [i]  pixels to extend the width  for enough space to display. [minimum:150]
	    -extht  [i]  pixels to extend the height for enough space to display. [minimum:150]

	    -h|help      Display this help info.

	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}      #-bkpos  [i]  the virus integration position(s) OR rearrangement break-point(s). [none]

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
			[ dHG_InputForLocalHapAlgorithm => undef ],
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
			# [ dHG_break_pos => [] ], # later change to Hash
			[ dHG_juncInfoFile => undef ],
			[ dHG_vsegInfoFile => undef ],
			[ dHG_sampleIDtoShow => undef ],
			# bam -> 'samtools depth'
			[ min_baseQ => 5 ],
			[ min_mapQ => 10 ],
			[ no_SoftClip => 0 ],
			# depth options
			[ dG_DepTrimRatio => 0 ],
			[ dG_DepEng_3SD => 0 ],
			[ dG_DepSmooth_flkWd => 2 ],
			[ dHG_ctrl_depthFluctR => 0 ],
			[ dHG_ctrl_depthFluctR_mode => 0 ],
			[ dHG_ctrl_depthFluctR_flor => undef ],
			[ dHG_ctrl_depthFluctR_ceil => undef ],
			[ dHG_ctrl_depthFluctR_minEffectRegResol => 10 ],
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
			# dHG_CNbaseOnReg -> output_file -> "prefix(para(-o)).SegCN_info.tsv"
			[ dHG_CNbaseOnReg => undef ],
			[ dHG_CNbaseOnReg_minDepthLimit => 5], # empirical
			[ dHG_TissueToShowCN => 'case' ],
			[ dHG_CNlineColor => 'blue' ],
			[ dHG_CNlineWidth => 1 ],
			# draw options
			[ dG_autoAdjustResolFigSize => 0 ],
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
			[ dHG_SegmentMode => 'u' ],
			[ dHG_SegmFillCol => 'lightblue' ],
			#--- intermediate variants
			[ DepthToDraw => undef ],
			[ dG_CaseVsCtrl_mode => undef ],
			[ dG_HasT_mode => undef ],
			[ dG_HasN_mode => undef ],
			[ case_Tcell_GMratio => undef ],
			[ case_Ncell_GMratio => undef ],
			# dHG_JuncInfo -> refseg_5p -> $refseg_5p
			# dHG_JuncInfo -> bkpos_5p  -> $bkpos_5p
			# dHG_JuncInfo -> strd_5p   -> $strd_5p
			# dHG_JuncInfo -> segID_5p  -> $segID_5p
			# dHG_JuncInfo -> refseg_3p -> $refseg_3p
			# dHG_JuncInfo -> bkpos_3p  -> $bkpos_3p
			# dHG_JuncInfo -> strd_3p   -> $strd_3p
			# dHG_JuncInfo -> segID_3p  -> $segID_3p
			# dHG_JuncInfo -> raw_JRnum -> $JR_count
			# dHG_JuncInfo -> puT_JRnum -> $puT_JRnum
			# dHG_JuncInfo -> puT_copyNumber -> $puT_copyNumber
			# dHG_JuncInfo -> N_copyratio -> $CopyRatio_NormalCell
			[ dHG_JuncInfo => {} ],
			# dHG_HsegInfo -> Segm_ID -> $Segm_ID
			# dHG_HsegInfo -> st_pos -> $st_pos
			# dHG_HsegInfo -> ed_pos -> $ed_pos
			# dHG_HsegInfo -> st_winNO -> $st_winNO
			# dHG_HsegInfo -> ed_winNO -> $ed_winNO
			# dHG_HsegInfo -> raw_meanDepth -> $raw_meanDepth
			# dHG_HsegInfo -> puT_meanDepth -> $puT_meanDepth
			# dHG_HsegInfo -> puT_copyNumber -> $puT_copyNumber
			[ dHG_HsegInfo => {} ],
			# dHG_VsegInfo -> Segm_ID -> $vseg_ID
			# dHG_VsegInfo -> vrefseg -> $vrefseg
			# dHG_VsegInfo -> st_pos -> $st_pos
			# dHG_VsegInfo -> ed_pos -> $ed_pos
			# dHG_VsegInfo -> raw_meanDepth -> $raw_meanDepth (just read from file)
			# dHG_VsegInfo -> puT_meanDepth -> $raw_meanDepth (just read from file)
			# dHG_VsegInfo -> puT_copyNumber -> $puT_copyNumber
			[ dHG_VsegInfo => {} ],
			# [ output_pref => undef ],
			[ CtrlPos2Depth => {} ],
			# check regResol_Href structure in FuseSV::Visual::ArrangeObj module
			[ dG_regResol => {} ],
			[ dG_PosToWinNOSets => [] ],
			# check winDepth_Href structure in FuseSV::Visual::DepthWork module
			[ dG_winDepth => {} ],
			# dG_bkWindow -> $window_NO = $bkpos
			[ dG_bkWindow => {} ],
			# dG_bkReason -> $window_NO = 'intra-SV' or 'inter-SV/vitg'
			[ dG_bkReason => {} ],
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
								  ['dHG_juncInfoFile'],
								  ['dHG_vsegInfoFile'],
								  ['trans_psl_idx_bgz'],
								  ['rmsk_idx_bgz'],
								  ['tabix'],
								  ['samtools']                    ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# objects constructions, as these objects belong to the input options
	$V_Href->{dHG_bam}->{case} = FuseSV::Virus_Integ::Bam_OB->new( tissue => 'case', filepath => undef );
	$V_Href->{dHG_bam}->{ctrl} = FuseSV::Virus_Integ::Bam_OB->new( tissue => 'ctrl', filepath => undef );

	# get options
	GetOptions(
		# input/output
		"-o:s"		=> \$V_Href->{dG_figure},
		"-t_bam:s"	=> \$V_Href->{dHG_bam}->{case}->{filepath},
		"-n_bam:s"	=> \$V_Href->{dHG_bam}->{ctrl}->{filepath},
		"-reg:s"	=> \$V_Href->{dHG_draw_region},
		"-rgsol:s"	=> \$V_Href->{dHG_reg_resol_file},
		"-auto_r:i"	=> \$V_Href->{dG_autoAdjustResolFigSize},
		# database
		"-tpsl:s"	=> \$V_Href->{trans_psl_idx_bgz},
		"-rmsk:s"	=> \$V_Href->{rmsk_idx_bgz},
		"-eltp:s"	=> \@{$V_Href->{dHG_rmsk_type}},
		"-gene:s"	=> \@{$V_Href->{dHG_show_gene}},
		"-nreg:s"	=> \@{$V_Href->{dG_user_gene}},
		# software
		"-stl:s"	=> \$V_Href->{samtools},
		"-tbi:s"	=> \$V_Href->{tabix},
		# sample info
		"-t_mdp:f"	=> \$V_Href->{dHG_meandepth}->{case},
		"-n_mdp:f"	=> \$V_Href->{dHG_meandepth}->{ctrl},
		"-t_puty:f"	=> \$V_Href->{dHG_casePurity},
		"-t_pldy:f"	=> \$V_Href->{dHG_casePloidy},
		"-n_pldy:f"	=> \$V_Href->{dHG_ctrlPloidy},
		"-id:s"		=> \$V_Href->{dHG_sampleIDtoShow},
		"-junc:s"	=> \$V_Href->{dHG_juncInfoFile},
		"-vseg:s"	=> \$V_Href->{dHG_vsegInfoFile},
		# bam 'depth'
		"-baseQ:i"	=> \$V_Href->{min_baseQ},
		"-mapQ:i"	=> \$V_Href->{min_mapQ},
		"-no_sfcp"	=> \$V_Href->{no_SoftClip},
		# depth options
		"-tmrt:f"	=> \$V_Href->{dG_DepTrimRatio},
		"-eng3sd"	=> \$V_Href->{dG_DepEng_3SD},
		"-fwsm:i"	=> \$V_Href->{dG_DepSmooth_flkWd},
		"-n_dptr:f"	=> \$V_Href->{dHG_ctrl_depthFluctR},
		"-n_dprsl:f"=> \$V_Href->{dHG_ctrl_depthFluctR_minEffectRegResol}, # advanced
		"-nprcr:s"	=> \@{$V_Href->{dHG_ctrl_prior_copyR}},
		"-cn_bs:s"	=> \$V_Href->{dHG_CNbaseOnReg},
		# draw options
		"-T_col:s"	=> \$V_Href->{dG_DepthSpecColor}->{case},
		"-N_col:s"	=> \$V_Href->{dG_DepthSpecColor}->{ctrl},
		"-bgbcol:s"	=> \$V_Href->{dG_bgStrokeCol},
		"-x_lnb:i"	=> \$V_Href->{dG_x_LabBaseNum},
		"-x_lbs:i"	=> \$V_Href->{dG_x_LabStep},
		"-x_lbfz:i"	=> \$V_Href->{dG_x_LabFtsz},
		"-x_rrc:s"	=> \$V_Href->{dG_x_RegResolColor},
		"-y_len:i"	=> \$V_Href->{dG_y_axisLen},
		"-y_res:f"	=> \$V_Href->{dG_y_resol},
		"-y_lbs:i"	=> \$V_Href->{dG_y_LabStep},
		"-y_lbfz:i"	=> \$V_Href->{dG_y_LabFtsz},
		"-extwd:i"	=> \$V_Href->{dG_extWidth},
		"-extht:i"	=> \$V_Href->{dG_extHeight},
		"-miew:f"	=> \$V_Href->{dHG_minEleWidth},
		"-segm:s"	=> \$V_Href->{dHG_SegmentMode},
		"-segc:s"	=> \$V_Href->{dHG_SegmFillCol},
		# help
		"-h|help"	=> \$V_Href->{HELP},
		# for debug
		"-debug"	=> \$V_Href->{in_debug} # hidden option
	);
}	# "-bkpos:i"	=> \@{$V_Href->{dHG_break_pos}},

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || $V_Href->{dG_figure} !~ /\.([sp][vnd][gf])$/i
		     || ! file_exist( filePath => $V_Href->{samtools} )
		     || ! file_exist( filePath => $V_Href->{tabix} )
		     || !(    file_exist( filePath => $V_Href->{dHG_bam}->{case}->{filepath} ) 
		     	   || file_exist( filePath => $V_Href->{dHG_bam}->{ctrl}->{filepath} )
		     	 )
		     || ! file_exist( filePath => $V_Href->{dHG_reg_resol_file} )
		     || !( $V_Href->{dHG_meandepth}->{case} || $V_Href->{dHG_meandepth}->{ctrl} )
		     || !( $V_Href->{dHG_draw_region} && $V_Href->{dHG_draw_region} =~ /^.+:\d+\-\d+$/ )
		     || (   ( $V_Href->{dHG_bam}->{case}->{filepath} && $V_Href->{dHG_bam}->{ctrl}->{filepath} )
		     	 && ( $V_Href->{dHG_meandepth}->{case} <= 0 || $V_Href->{dHG_meandepth}->{ctrl} <= 0 )
				)
			 || $V_Href->{dG_DepSmooth_flkWd} < 0
			 || $V_Href->{dG_x_LabStep} <= 0
			 || $V_Href->{dG_y_resol} < 0
			 || $V_Href->{dG_y_LabStep} < 0
			 || $V_Href->{dHG_SegmentMode} =~ /\s/
			);
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub real_draw_segmemt_CN{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	#--------------#
	# prepare work #
	#--------------#
	&prepare;

	#---------------------------------#
	# calculate window depth from bam #
	#---------------------------------#
	&load_window_adjdepth;

	#--------------------------#
	# get show depth of window #
	#--------------------------#
	&get_window_showdepth;

	#-----------------------#
	# load genomic database #
	#-----------------------#
	# load transcript info from idx psl.bgz
	&load_trans_info;
	# load rmsk info from idx rmsk.bgz
	&load_rmsk_info if( defined $V_Href->{rmsk_idx_bgz} );

	#-------------#
	# draw figure #
	#-------------#
	# set attributes of x axis #
	&set_x_axis;
	# set attributes of y axis #
	&set_y_axis;
	# draw rmsk elements in this region #
	# create SVG object, maybe
	&draw_rmsk if( defined $V_Href->{rmsk_idx_bgz} );
	# draw y axis #
	# create SVG object, maybe
	&draw_y_axis;
	# draw depth spectrum #
	&draw_depth;
	# draw gene in this region #
	&draw_genes;
	# draw x axis #
	&draw_x_axis;
	# show breakpoints
	&show_breakpoints;
	# show basic information and contact
	&show_figure_information;
	# output SVG #
	&generate_figure_output;
}

#--- prepare work ---
sub prepare{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	#--- verify bam and bai
	for my $tissue (qw/ case ctrl /){
		my $bam_OB = $V_Href->{dHG_bam}->{$tissue};
		if( defined $bam_OB->{filepath} ){
			$bam_OB->verify_bam;
			$bam_OB->verify_index;
		}
	}

	#--- depth to draw
	$V_Href->{DepthToDraw} = $V_Href->{dHG_meandepth}->{case} || $V_Href->{dHG_meandepth}->{ctrl};

	#--- change array input to hash
	# for my $key (qw/ dHG_break_pos dHG_rmsk_type dHG_show_gene /){
	for my $key (qw/ dHG_rmsk_type dHG_show_gene /){
		my %temp = map{($_,1)} @{$V_Href->{$key}};
		$V_Href->{$key} = \%temp;
	}

	$V_Href->{dG_HasT_mode} = defined $V_Href->{dHG_bam}->{case}->{filepath};
	$V_Href->{dG_HasN_mode} = defined $V_Href->{dHG_bam}->{ctrl}->{filepath};
	$V_Href->{dG_CaseVsCtrl_mode} = ( $V_Href->{dG_HasT_mode} && $V_Href->{dG_HasN_mode} );
	#--- genome material ratio in case sample: mixed tumor cells and normal cells
	if( $V_Href->{dG_HasT_mode} ){
		my $case_Tcell_GMpart = get_Tcell_GMpart( TumorCellPurity=>$V_Href->{dHG_casePurity},  TumorCellPloidy=>$V_Href->{dHG_casePloidy} );
		my $case_Ncell_GMpart = get_Ncell_GMpart( TumorCellPurity=>$V_Href->{dHG_casePurity}, NormalCellPloidy=>$V_Href->{dHG_ctrlPloidy} );
		$V_Href->{case_Tcell_GMratio} = get_Tcell_GMratio( Tcell_GMpart=>$case_Tcell_GMpart, Ncell_GMpart=>$case_Ncell_GMpart );
		$V_Href->{case_Ncell_GMratio} = get_Ncell_GMratio( Tcell_GMpart=>$case_Tcell_GMpart, Ncell_GMpart=>$case_Ncell_GMpart );
		my $normalize_meandepth = sprintf "%.2f", $V_Href->{DepthToDraw} * $V_Href->{case_Tcell_GMratio};
		# inform
		stout_and_sterr `date`."[INFO]\tIn tumor sample, genome material ratio of Tumor  cell is $V_Href->{case_Tcell_GMratio}.\n"
							  ."[INFO]\tIn tumor sample, genome material ratio of Normal cell is $V_Href->{case_Ncell_GMratio}.\n"
							  # note that,
							  # for Normal-sample-depth, just scales to $normalize_meandepth.
							  # for  Tumor-sample-depth, subtracts depth belongs to the mixed normal-cells.
							  ."[INFO]\tFor comparison, all depth value will normalized to mean depth $normalize_meandepth X.\n";
	}
	else{ # sololy one sample, no comparison, no adjustment.
		stout_and_sterr `date`."[INFO]\tDraw in mode of mean depth: $V_Href->{DepthToDraw} X.\n";
	}

	#--- genomic region to draw
	my ($reg_Seg, $reg_stP, $reg_edP) = ( $V_Href->{dHG_draw_region} =~ /^(.+):(\d+)\-(\d+)$/ );
	if( $reg_stP >= $reg_edP ){
		warn_and_exit "<ERRO>\tError found in Region ($V_Href->{dHG_draw_region}) to draw.\n";
	}
	$V_Href->{dHG_draw_region} = { refseg=>$reg_Seg, st_pos=>$reg_stP, ed_pos=>$reg_edP };

	#--- load region resolution
	load_regResolForWindowlize(
								reg_resol_file => $V_Href->{dHG_reg_resol_file},
								regResol_Href => $V_Href->{dG_regResol},
								destReg_Href => $V_Href->{dHG_draw_region},
								destRegAutoFigureSize => $V_Href->{dG_autoAdjustResolFigSize}
							  );

	#--- prepare for get pos windowNO
	$V_Href->{dG_PosToWinNOSets} = [
									refseg => $V_Href->{dHG_draw_region}->{refseg},
									regResol_Href => $V_Href->{dG_regResol},
									stpos_1stwin => $V_Href->{dHG_draw_region}->{st_pos}
								   ];

	#--- prepare for OK-pos for windepth based on {dHG_ctrl_depthFluctR}
	if(    $V_Href->{dG_HasN_mode}
		&& $V_Href->{dHG_ctrl_depthFluctR} > 0
	){
		$V_Href->{dHG_ctrl_depthFluctR_mode} = 1;
		$V_Href->{dHG_ctrl_depthFluctR_flor} = $V_Href->{dHG_meandepth}->{ctrl} * ( 1 - $V_Href->{dHG_ctrl_depthFluctR} );
		$V_Href->{dHG_ctrl_depthFluctR_ceil} = $V_Href->{dHG_meandepth}->{ctrl} * ( 1 + $V_Href->{dHG_ctrl_depthFluctR} );
		stout_and_sterr "[INFO]\tAllowed depth region of Control position: $V_Href->{dHG_ctrl_depthFluctR_flor}X ~ $V_Href->{dHG_ctrl_depthFluctR_ceil}X\n";
	}

	#--- read junction info file ---
	&load_juncInfo;

	#--- read virus segment info file ---
	&load_vsegInfo;

	#--- map junction to segments ---
	&map_junction_to_segments;

	#--- priori known copy ratio exists in Control sample
	for my $i ( 1 .. scalar(@{$V_Href->{dHG_ctrl_prior_copyR}}) ){
		my $ctrl_prior_copyR_input = shift @{$V_Href->{dHG_ctrl_prior_copyR}};
		my @ele = split /[\-\:]/, $ctrl_prior_copyR_input;
		next if( $ele[0] ne $V_Href->{dHG_draw_region}->{refseg} );
		push @{$V_Href->{dHG_ctrl_prior_copyR}}, { refseg=>$ele[0], st_pos=>$ele[1], ed_pos=>$ele[2], copyR=>$ele[3] };
	}

	#--- for calculate CN of each segment
	&deal_with_showCNopt;

	#--- prepare the input file to haplotype algorithm
	my ($output_prefix) = ($V_Href->{dG_figure} =~ /^(.+)\.[sp][vnd][gf]$/i);
	$V_Href->{dHG_InputForLocalHapAlgorithm} = "$output_prefix.InputForLocalHapAlgorithm.txt";
}

#--- load juntion info ---
# including virus integration position(s) and rearrangement break-point(s)
# load up $V_Href->{dHG_JuncInfo}
# load up $V_Href->{dHG_HsegInfo}
# load up $V_Href->{dG_bkWindow}
sub load_juncInfo{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $JuncInfo_Href = $V_Href->{dHG_JuncInfo};
	my $HsegInfo_Href = $V_Href->{dHG_HsegInfo};
	my $BKwinInfo_Href= $V_Href->{dG_bkWindow};
	my $BKreason_Href = $V_Href->{dG_bkReason};

	return if( !defined $V_Href->{dHG_juncInfoFile} );

	my $junc_NO = 1;
	my $BKposForSegm_Href = {};
	open (JCIF, Try_GZ_Read($V_Href->{dHG_juncInfoFile})) || die "fail read juncInfo file: $!\n";
	while (<JCIF>){
		next if(/^#/);
		my ($refseg_5p, $bkpos_5p, $strd_5p, $refseg_3p, $bkpos_3p, $strd_3p, $JR_count, $CopyRatio_NormalCell) = (split)[0..7];
		$JuncInfo_Href->{$junc_NO} = {
										refseg_5p => $refseg_5p,
										bkpos_5p  => $bkpos_5p,
										strd_5p   => $strd_5p,
										segID_5p  => undef,
										refseg_3p => $refseg_3p,
										bkpos_3p  => $bkpos_3p,
										strd_3p   => $strd_3p,
										segID_3p  => undef,
										raw_JRnum => $JR_count,
										puT_JRnum => undef,
										puT_copyNumber => undef,
										N_copyratio => ($CopyRatio_NormalCell || 0) # default is somatic
									 };
		# update
		$junc_NO ++;

		# record break-point
		my $inDrawReg_5p = &is_in_draw_region( refseg=>$refseg_5p, pos_1=>$bkpos_5p );
		my $inDrawReg_3p = &is_in_draw_region( refseg=>$refseg_3p, pos_1=>$bkpos_3p );
		for my $Aref ( [ '5', $refseg_5p, $bkpos_5p, $strd_5p ],
					   [ '3', $refseg_3p, $bkpos_3p, $strd_3p ]
		){
			my ($prime, $refseg, $bkpos, $strd) = @$Aref;
			next unless &is_in_draw_region( refseg=>$refseg, pos_1=>$bkpos );
			# record for segmentation
			$BKposForSegm_Href->{$bkpos}->{"$prime,$strd"} = 1;
			# record for breakage annotation
			my $winNO = get_pos_windowNO( pos => $bkpos, @{$V_Href->{dG_PosToWinNOSets}} );
			$BKwinInfo_Href->{$winNO} = $bkpos;
			 $BKreason_Href->{$winNO} = ( $inDrawReg_5p && $inDrawReg_3p ) ? 'intra-SV' : 'inter-SV/vitg';
		}
	}
	close JCIF;

	# segmentations
	## load draw region
	$BKposForSegm_Href->{ $V_Href->{dHG_draw_region}->{st_pos} }->{'3,+'} = 1;
	$BKposForSegm_Href->{ $V_Href->{dHG_draw_region}->{ed_pos} }->{'5,+'} = 1;
	## gather bkpos-orientation
	my @BKposForSegm = ();
	for my $pos ( sort {$a<=>$b} keys %$BKposForSegm_Href ){
		# do first!
		if(    exists $BKposForSegm_Href->{$pos}->{'3,+'}
			|| exists $BKposForSegm_Href->{$pos}->{'5,-'}
		){
			push @BKposForSegm, [ $pos, -1 ];
		}
		# then do this!
		if(    exists $BKposForSegm_Href->{$pos}->{'5,+'}
			|| exists $BKposForSegm_Href->{$pos}->{'3,-'}
		){
			push @BKposForSegm, [ $pos, 0 ];
		}
	}
	## construct segments info
	my $Hseg_NO = 1;
	my $LastHsegEdpos = $BKposForSegm[0]->[0] - 1;
	shift @BKposForSegm; # first one must be $V_Href->{dHG_draw_region}->{st_pos} with '3,+'
	for my $Aref ( @BKposForSegm ){
		my ($pos, $shift) = @$Aref;
		# segment ends
		my $segm_stpos = $LastHsegEdpos + 1;
		my $segm_edpos = $pos + $shift;
		# filter, important!
		next if( $segm_edpos < $segm_stpos );
		# ends' window NO
		my $st_winNO = get_pos_windowNO( pos => $segm_stpos, @{$V_Href->{dG_PosToWinNOSets}} );
		my $ed_winNO = get_pos_windowNO( pos => $segm_edpos, @{$V_Href->{dG_PosToWinNOSets}} );
		# load host segment
		$HsegInfo_Href->{$Hseg_NO} = {
										Segm_ID => &get_HostSegmentID( Segm_NO => $Hseg_NO ),
										st_pos => $segm_stpos,
										ed_pos => $segm_edpos,
										st_winNO => $st_winNO,
										ed_winNO => $ed_winNO,
										raw_meanDepth => undef,
										puT_meanDepth => undef,
										puT_copyNumber => undef
									 };
		# update
		$Hseg_NO ++;
		$LastHsegEdpos = $segm_edpos;
	}

	# inform
	stout_and_sterr `date`."[INFO]\tload breakage junction information OK.\n";
}

#--- load virus segment info ---
# load up $V_Href->{dHG_VsegInfo}
sub load_vsegInfo{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $VsegInfo_Href = $V_Href->{dHG_VsegInfo};

	return if( !defined $V_Href->{dHG_vsegInfoFile} );

	my $Vseg_NO = 1;
	open (VSEG, Try_GZ_Read($V_Href->{dHG_vsegInfoFile})) || die "fail read vsegInfo file: $!\n";
	while (<VSEG>){
		next if(/^#/);
		my ($vseg_ID, $vrefseg, $st_pos, $ed_pos, $mean_depth) = (split)[0..4];
		$VsegInfo_Href->{$Vseg_NO} = {
										Segm_ID => $vseg_ID,
										vrefseg => $vrefseg,
										st_pos => $st_pos,
										ed_pos => $ed_pos,
										raw_meanDepth => $mean_depth,
										puT_meanDepth => $mean_depth, # consider virus segment is more than we need
										puT_copyNumber => undef
									 };
		# update
		$Vseg_NO ++;
	}
	close VSEG;

	# inform
	stout_and_sterr `date`."[INFO]\tload virus segment information OK.\n";
}

#--- map junctions to host/virus segments ---
sub map_junction_to_segments{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $JuncInfo_Href = $V_Href->{dHG_JuncInfo};
	my $HsegInfo_Href = $V_Href->{dHG_HsegInfo};
	my $VsegInfo_Href = $V_Href->{dHG_VsegInfo};

	return if( !defined $V_Href->{dHG_juncInfoFile} );

	for my $junc_NO ( sort {$a<=>$b} keys %$JuncInfo_Href ){
		my $this_junc_Href = $JuncInfo_Href->{$junc_NO};
		for my $primer ('5p', '3p'){
			my $refseg = $this_junc_Href->{"refseg_$primer"};
			my $bkpos  = $this_junc_Href->{"bkpos_$primer"};
			my $strd   = $this_junc_Href->{"strd_$primer"};
			my $BKtype = "$primer;$strd";
			my $posLoc = ( $BKtype eq "5p;+" || $BKtype eq "3p;-" ) ? 'ed_pos' : 'st_pos';
			my $segm_Href = &is_in_draw_region( refseg=>$refseg, pos_1=>$bkpos ) ? $HsegInfo_Href : $VsegInfo_Href;
			my $segmID = &get_SegmID_by_PosLoc( segm_Href=>$segm_Href, posLoc=>$posLoc, pos=>$bkpos );
			$this_junc_Href->{"segID_$primer"} = $segmID;
		}
	}

	# inform
	stout_and_sterr `date`."[INFO]\tmap junctions to host/virus segments OK.\n";
}

#--- get segment ID based on Position Location ---
sub get_SegmID_by_PosLoc{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $segm_Href = $parm{segm_Href};
	my $posLoc = $parm{posLoc};
	my $pos = $parm{pos};

	for my $segm_NO ( sort {$a<=>$b} keys %$segm_Href ){
		if( $segm_Href->{$segm_NO}->{$posLoc} == $pos ){
			return $segm_Href->{$segm_NO}->{Segm_ID};
		}
	}
	# cannot find
	warn_and_exit `date`."<ERROR>\tCannot find segment ID for position $pos as $posLoc.\n";
}

#--- test whether in draw region ---
sub is_in_draw_region{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $refseg = $parm{refseg};
	my $pos_1  = $parm{pos_1};
	my $pos_2  = $parm{pos_2};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	if( $refseg eq $V_Href->{dHG_draw_region}->{refseg} ){
		if( !defined $pos_2 ){
			if(    $pos_1  >= $V_Href->{dHG_draw_region}->{st_pos}
				&& $pos_1  <= $V_Href->{dHG_draw_region}->{ed_pos}
			){
				return 1;
			}
			else{
				return 0;
			}
		}
		else{
			my $ovp_len = Get_Two_Seg_Olen( $pos_1, $pos_2, $V_Href->{dHG_draw_region}->{st_pos}, $V_Href->{dHG_draw_region}->{ed_pos} );
			if( $ovp_len != 0 ){
				return 1;
			}
			else{
				return 0;
			}
		}
	}
	else{
		return 0;
	}
}

#--- for calculate CN of each segment ---
sub deal_with_showCNopt{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# set the tissue whose CN to show, 'case' has priori
	$V_Href->{dHG_TissueToShowCN} = ( $V_Href->{dG_HasT_mode} ? 'case' : 'ctrl' );

	return if( !defined $V_Href->{dHG_CNbaseOnReg} );

	my @CNbaseOnReg_ele = split /:/, $V_Href->{dHG_CNbaseOnReg};
	my ($tissue, $refseg, $interval, $cn);
	my ($st_pos, $ed_pos);
	if( $CNbaseOnReg_ele[0] =~ /^[TN]$/i ){
		($tissue, $refseg, $interval, $cn) = @CNbaseOnReg_ele;
		$tissue = ( $tissue =~ /t/i ? 'case' : 'ctrl' );
	}
	else{
		# 'control' has priori
		$tissue = ( $V_Href->{dG_HasN_mode} ? 'ctrl' : 'case' );
		($refseg, $interval, $cn) = @CNbaseOnReg_ele;
	}
	# default cn
	if(    !defined $cn
		|| $cn !~ /^\d+(\.\d+)?$/
		|| $cn <= 0.5
	){
		$cn = 2;
		stout_and_sterr "<WARN>\tReset CN to 2 as the cn in input '-cn_bs' is not available.\n";
	}
	# check refseg
	if( $refseg ne $V_Href->{dHG_draw_region}->{refseg} ){
		warn_and_exit "<ERROR>\tInput '-cn_bs' is from refseg ($refseg) different from your region ($V_Href->{dHG_draw_region}->{refseg}).\n";
	}
	# check interval
	if( $interval =~ /^(\d+)\-(\d+)$/ ){
		# check overlap with view region
		my @clip_interval = 
			arrange_region_with_clip(
										lftPos => $1,
										rgtpos => $2,
										min_Pos => $V_Href->{dHG_draw_region}->{st_pos},
										max_Pos => $V_Href->{dHG_draw_region}->{ed_pos}
									);
		# check rationality
		if(    @clip_interval == 0
			|| $clip_interval[0]->[0] >= $clip_interval[0]->[1]
		){
			warn_and_exit "<ERROR>\tRegion of '-cn_bs' input is not available.\n $clip_interval[0]->[0] ; $clip_interval[0]->[1]\n";
		}
		else{
			($st_pos, $ed_pos) = @{ $clip_interval[0] };
		}
	}
	else{
		warn_and_exit "<ERROR>\tCannot recognize region of '-cn_bs' input.\n";
	}
	# check breakpoint
	for my $bk_pos ( sort {$a<=>$b} values %{$V_Href->{dG_bkWindow}} ){
		if(    $bk_pos > $st_pos
			&& $bk_pos < $ed_pos
		){
			warn_and_exit "<ERROR>\tDo not allow any breakpoit ($bk_pos) exists in the region of '-cn_bs' input.\n";
		}
	}
	# load up
	$V_Href->{dHG_CNbaseOnReg} = {
									orig_input => $V_Href->{dHG_CNbaseOnReg},
									tissue => $tissue,
									refseg => $refseg,
									st_pos => $st_pos,
									ed_pos => $ed_pos,
									cn => $cn,
									depth2cn_ratio => undef
								};
	# get window NO
	$V_Href->{dHG_CNbaseOnReg}->{"${_}winNO"} = 
		get_pos_windowNO(
			pos => $V_Href->{dHG_CNbaseOnReg}->{"${_}pos"},
			@{$V_Href->{dG_PosToWinNOSets}}
		) for qw / st_ ed_ /;
}

#--- load depth of given region from bam ---
#--- transform to adjusted depth ---
sub load_window_adjdepth{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# use samtools depth
	# first is ctrl to get ok pos based on {dHG_ctrl_depthFluctR}
	# then  is case to discard once {dHG_ctrl_depthFluctR} works
	for my $tissue ('ctrl', 'case'){
		my $bam_OB = $V_Href->{dHG_bam}->{$tissue};
		next unless( $bam_OB->{filepath} );
		# depth stat of genomic positions in draw region
		my $FilePrefix = $V_Href->{dG_figure}.".tempfile.$tissue";
		$bam_OB->get_region_depth(
									refseg => $V_Href->{dHG_draw_region}->{refseg},
									st_pos => $V_Href->{dHG_draw_region}->{st_pos},
									ed_pos => $V_Href->{dHG_draw_region}->{ed_pos},
									min_baseQ => $V_Href->{min_baseQ},
									min_mapQ => $V_Href->{min_mapQ},
									no_SoftClip => $V_Href->{no_SoftClip},
									out_prefix => $FilePrefix,
								);
		# get depth of windows, tranform to adjdepth
		&read_depthStat_to_window_adjdepth( tissue => $tissue, depthStatFile => $bam_OB->{regionDepthFile} );
		# delete depth stat file
		$bam_OB->delete_regionDepthFile;
	}
	# release memory
	$V_Href->{CtrlPos2Depth} = {};
}

#--- get show depth (smooth) ---
sub get_window_showdepth{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# calculate show-depth of windows
	for my $tissue ('case', 'ctrl'){
		my $bam_OB = $V_Href->{dHG_bam}->{$tissue};
		next unless( $bam_OB->{filepath} );
		# smooth depth
		get_windowSmoDepth(
							winDepth_Href => $V_Href->{dG_winDepth},
							tissue => $tissue,
							DepSmooth_flkWd => $V_Href->{dG_DepSmooth_flkWd},
							smoDepthKey => 'showdepth',
							bkWindow_Href => $V_Href->{dG_bkWindow},
							maxDepth_Sref => \$V_Href->{dG_maxDepthToDraw},
							maxWinNO_Sref => \$V_Href->{dG_maxWinNOtoDraw}
						  );
		# inform
		if( $V_Href->{dG_DepSmooth_flkWd} != 0 ){
			stout_and_sterr `date`."[INFO]\tSmooth window depth of tissue $tissue OK.\n";
		}
	}

	# inform
	stout_and_sterr `date`."[INFO]\tMaximum depth to dispaly is $V_Href->{dG_maxDepthToDraw}X.\n";
}

#--- load transcript info ---
sub load_trans_info{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	if( defined $V_Href->{trans_psl_idx_bgz} ){
		# for idx bgz reading
		my $idx_readReg =    $V_Href->{dHG_draw_region}->{refseg}.
						 ':'.$V_Href->{dHG_draw_region}->{st_pos}.
						 '-'.$V_Href->{dHG_draw_region}->{ed_pos};
		my $g2tid_Href = {};
		# load trans info from idx tpsl.bgz file
		load_transInfo_from_tpsl(
									tpsl_file => $V_Href->{trans_psl_idx_bgz},
									TransInfo_Href => $V_Href->{dG_TransInfo},
									idx_readReg => $idx_readReg,
									tabix => $V_Href->{tabix},
									clip_min_Pos => $V_Href->{dHG_draw_region}->{st_pos},
									clip_max_Pos => $V_Href->{dHG_draw_region}->{ed_pos},
									onlyGene_Href => $V_Href->{dHG_show_gene},
									g2tid_Href => $g2tid_Href
								);
		# select trans
		for my $gene_name (keys %$g2tid_Href){
			my @trans_id = sort keys %{ $g2tid_Href->{$gene_name}};
			for my $keep_tid ( '001', '201', $trans_id[0] ){
				if( exists $g2tid_Href->{$gene_name}->{$keep_tid} ){
					delete $V_Href->{dG_TransInfo}->{"$gene_name-$_"} for grep $_ ne $keep_tid, @trans_id;
					last;
				}
			}
		}
		# inform
		stout_and_sterr `date`."[INFO]\tLoad transcript PSL file OK.\n"
									 ."\t$V_Href->{trans_psl_idx_bgz}\n";
	}

	# user defined region
	if( scalar( @{$V_Href->{dG_user_gene}} ) != 0 ){
		add_user_region_to_gene_info(
										user_gene_input_Aref => $V_Href->{dG_user_gene},
										TransInfo_Href => $V_Href->{dG_TransInfo},
										refseg => $V_Href->{dHG_draw_region}->{refseg},
										clip_min_Pos => $V_Href->{dHG_draw_region}->{st_pos},
										clip_max_Pos => $V_Href->{dHG_draw_region}->{ed_pos},
										input_optKey_forAlert => '-nreg'
									);
		# inform
		stout_and_sterr `date`."[INFO]:\tLoad user defined region OK.\n";
	}
}

#--- load rmsk info ---
sub load_rmsk_info{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	my $idx_readReg =    $V_Href->{dHG_draw_region}->{refseg}.
					 ':'.$V_Href->{dHG_draw_region}->{st_pos}.
					 '-'.$V_Href->{dHG_draw_region}->{ed_pos};
	# load rmsk info from idx rmsk.bgz file
	load_rmskInfo_from_idxBgz(
								rmsk_idx_bgz => $V_Href->{rmsk_idx_bgz},
								RmskInfoPool_Href => $V_Href->{dG_RmskInfo},
								idx_readReg => $idx_readReg,
								tabix => $V_Href->{tabix},
								clip_min_Pos => $V_Href->{dHG_draw_region}->{st_pos},
								clip_max_Pos => $V_Href->{dHG_draw_region}->{ed_pos},
								onlyEleType_Href => $V_Href->{dHG_rmsk_type}
							);

	# inform
	stout_and_sterr `date`."[INFO]\tLoad rmsk file OK.\n"
								 ."\t$V_Href->{rmsk_idx_bgz}\n";
}

#--- set attributes of x axis ---
sub set_x_axis{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$V_Href->{dG_axisZX} = $V_Href->{dG_extWidth};
	$V_Href->{dG_x_axisLen} = $V_Href->{dG_maxWinNOtoDraw};

	# inform
	stout_and_sterr `date`."[INFO]\tSet X axis OK.\n";
}

#--- set attributes of y axis ---
sub set_y_axis{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# set the resol and lab step of Y axis
	set_Y_axis_Resol_LabStep(
								yResol_Sref => \$V_Href->{dG_y_resol},
								yMaxValue_Sref => \$V_Href->{dG_maxDepthToDraw},
								yAxisLen_Sref => \$V_Href->{dG_y_axisLen},
								yLabStep_Sref => \$V_Href->{dG_y_LabStep},
								disMaxValueUpdate => 0
							);

	# inform
	stout_and_sterr `date`."[INFO]\tSet Y axis OK.\n";
}

#--- draw rmsk ele ---
## create SVG object, maybe
sub draw_rmsk{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# return if no rmsk ele to show
	if( scalar(keys %{$V_Href->{dG_RmskInfo}}) == 0 ){
		stout_and_sterr `date`."<WARN>:\tNone rmsk elements to show!\n";
		return;
	}

	allocate_RMSKele_and_draw(
								SVG_obj_Oref => \$V_Href->{general_SVG_obj},
								RmskInfoPool_Href => $V_Href->{dG_RmskInfo},
								pos2wNO_para_Aref => $V_Href->{dG_PosToWinNOSets},
								axisZX => $V_Href->{dG_axisZX},
								draw_baseY_Sref => \$V_Href->{general_BaseY},
								eleMinWidth => $V_Href->{dHG_minEleWidth},
								labelFtfm => $V_Href->{general_FontFam},
								labelFtsz => $V_Href->{general_FontSize},
								# for initialize SVG obj
								xAxisLen => $V_Href->{dG_x_axisLen},
								yAxisLen => $V_Href->{dG_y_axisLen},
								bi_extWidth => $V_Href->{dG_extWidth},
								bi_extHeight => $V_Href->{dG_extHeight},
								bgStrokeCol => $V_Href->{dG_bgStrokeCol}
							 );

	# inform
	stout_and_sterr `date`."[INFO]\tShow rmsk element OK.\n";
}

#--- draw y axis ---
## create SVG object, maybe
sub draw_y_axis{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# depth spectrum color legend info near Y axis
	my @yLegd_info = ();
	for my $tissue ( 'case', 'ctrl' ){
		next unless( $V_Href->{dHG_bam}->{$tissue}->{filepath} );
		push @yLegd_info, {
							shape => 'p',
							fill_col => $V_Href->{dG_DepthSpecColor}->{$tissue},
							text => "$tissue depth",
							font_family => $V_Href->{general_FontFam_Bold},
							font_size => $V_Href->{general_FontSize}
						  };
	}
	my $yLegd_Href = {loc=>'top', info_Aref=>\@yLegd_info};

	# draw Y axis
	show_Y_axis(
				 SVG_obj_Oref => \$V_Href->{general_SVG_obj},
				 # attributes
				 yResol => $V_Href->{dG_y_resol},
				 yMaxValue => $V_Href->{dG_maxDepthToDraw},
				 yAxisLen => $V_Href->{dG_y_axisLen},
				 yLabStep => $V_Href->{dG_y_LabStep},
				 # location
				 draw_baseY => $V_Href->{general_BaseY},
				 axisZX => $V_Href->{dG_axisZX},
				 axisZY_Sref => \$V_Href->{dG_axisZY},
				 bi_extWidth => $V_Href->{dG_extWidth},
				 bi_extHeight => $V_Href->{dG_extHeight},
				 # initialize SVG obj
				 xAxisLen => $V_Href->{dG_x_axisLen},
				 bgStrokeCol => $V_Href->{dG_bgStrokeCol},
				 # appearance
				 yLabelUnit => 'X',
				 yLabelFtsz => $V_Href->{dG_y_LabFtsz},
				 yLabelFtfm => $V_Href->{general_FontFam},
				 yLegd_Href => $yLegd_Href
			);

	# inform
	stout_and_sterr `date`."[INFO]\tDraw Y axis OK.\n";
}

#--- draw viral depth spectrum of each tissue ---
sub draw_depth{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# draw depth spectrum of each tissue
	draw_depth_spectrum(
							SVG_obj_Oref => \$V_Href->{general_SVG_obj},
							winDepth_Href => $V_Href->{dG_winDepth},
							smoDepthKey => 'showdepth',
							maxWinNOtoDraw => $V_Href->{dG_maxWinNOtoDraw},
							axisZX => $V_Href->{dG_axisZX},
							axisZY => $V_Href->{dG_axisZY},
							yResol => $V_Href->{dG_y_resol},
							SpmCol_Href => $V_Href->{dG_DepthSpecColor}
						);

	# record the bottom Y of depth spectrum field
	$V_Href->{dG_DepthSpecBtom_Y} = $V_Href->{dG_axisZY};

	# to sperate gene from depth spectrum, boundary
	if( scalar(keys %{$V_Href->{dG_TransInfo}}) != 0 ){
		$V_Href->{general_SVG_obj}
			   ->line(
						x1 => $V_Href->{dG_axisZX},
						y1 => $V_Href->{dG_axisZY},
						x2 => $V_Href->{dG_axisZX} + $V_Href->{dG_x_axisLen},
						y2 => $V_Href->{dG_axisZY},
						stroke => 'black',
						'stroke-width' => 1,
						'stroke-linecap' => 'round'
					 );
	}

	# inform
	stout_and_sterr `date`."[INFO]\tDraw depth spectrum OK.\n";
}

#--- draw genes ---
sub draw_genes{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# return if no rmsk ele to show
	if( scalar(keys %{$V_Href->{dG_TransInfo}}) == 0 ){
		stout_and_sterr `date`."<WARN>:\tNone genes to show!\n";
		return;
	}

	allocate_genes_and_draw(
							SVG_obj_Oref => \$V_Href->{general_SVG_obj},
							TransInfo_Href => $V_Href->{dG_TransInfo},
							pos2wNO_para_Aref => $V_Href->{dG_PosToWinNOSets},
							axisZX_Sref => \$V_Href->{dG_axisZX},
							axisZY_Sref => \$V_Href->{dG_axisZY},
							tnameFtfm => $V_Href->{general_FontFam},
							tnameFtsz => $V_Href->{general_FontSize},
							ExonMinWidth => $V_Href->{dHG_minEleWidth}
						   );

	# inform
	stout_and_sterr `date`."[INFO]\tShow genes OK.\n";
}

#--- draw x axis ---
sub draw_x_axis{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# x-axis Label Text
	my $xLabelText = $V_Href->{dHG_draw_region}->{refseg};
	$xLabelText .= "\n".$V_Href->{dHG_sampleIDtoShow} if( defined $V_Href->{dHG_sampleIDtoShow} );
	# show X axis
	show_X_axis(
				SVG_obj_Oref => \$V_Href->{general_SVG_obj},
				# location
				axisZX => $V_Href->{dG_axisZX},
				axisZY_Sref => \$V_Href->{dG_axisZY},
				draw_baseY_Sref => \$V_Href->{general_BaseY}, # for below display
				# attributes
				xAxisLen => $V_Href->{dG_x_axisLen},
				xAxisWidth => 1,
				xLabelSpanPix => $V_Href->{dG_x_LabStep},
				xLabelNumBase => $V_Href->{dG_x_LabBaseNum},
				xLabelFtsz => $V_Href->{dG_x_LabFtsz},
				xLabelFtfm => $V_Href->{general_FontFam},
				xLabelText => $xLabelText,
				# for region in different resolution
				xRegResolColor => $V_Href->{dG_x_RegResolColor},
				xResol_Href => $V_Href->{dG_regResol},
				regSegName => $V_Href->{dHG_draw_region}->{refseg},
				xRegMinPos => $V_Href->{dHG_draw_region}->{st_pos},
				xRegMaxPos => $V_Href->{dHG_draw_region}->{ed_pos}
			);

	# inform
	stout_and_sterr `date`."[INFO]\tDraw X axis OK.\n";
}

#--- the copy number informtion ---
# of host-segments, junctions, virus-segments
sub calculate_CopyNumber{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# calculate the depth2CN ratio
	if(    $V_Href->{dG_HasT_mode} # must have case sample
		&& !defined $V_Href->{dHG_CNbaseOnReg}
	){ # HapChrDepth mode
		$V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio} = 1.0 / $V_Href->{dHG_casePurity} /
			get_ObjSingleCNdepthInMixed(
										MixedCellMeanDepth => $V_Href->{dHG_meandepth}->{case},
										TumorCellPurity => $V_Href->{dHG_casePurity},
										TumorCellPloidy => $V_Href->{dHG_casePloidy},
										NormalCellPloidy => $V_Href->{dHG_ctrlPloidy}
									);
	}
	elsif( defined $V_Href->{dHG_CNbaseOnReg} ){ # user defined marker
		# first calculate the base-region mean_adjdepth
		$V_Href->{dHG_CNbaseOnReg}->{mean_adjdepth} = 
			get_givenWinItvalMeanDepth(
										winDepth_Href => $V_Href->{dG_winDepth},
										tissue => $V_Href->{dHG_CNbaseOnReg}->{tissue},
										depth_key => 'mean_adjdepth',
										st_winNO => $V_Href->{dHG_CNbaseOnReg}->{st_winNO},
										ed_winNO => $V_Href->{dHG_CNbaseOnReg}->{ed_winNO}
									);
		# check
		if( $V_Href->{dHG_CNbaseOnReg}->{mean_adjdepth} <= $V_Href->{dHG_CNbaseOnReg_minDepthLimit} ){
			stout_and_sterr "<WARN>\tmean depth ($V_Href->{dHG_CNbaseOnReg}->{mean_adjdepth}X) of CN-basement region is too low to be avaiable.\n"
								 ."\ttissue:$V_Href->{dHG_CNbaseOnReg}->{tissue}, $V_Href->{dHG_CNbaseOnReg}->{orig_input}\n";
			$V_Href->{dHG_CNbaseOnReg} = undef; # bye-bye
		}
		else{ # get the depth-To-cn transforming ratio
			$V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio} = $V_Href->{dHG_CNbaseOnReg}->{cn} / $V_Href->{dHG_CNbaseOnReg}->{mean_adjdepth};
		}
	}

	# calculate the copy number of host-segments, and output
	my $SOURCE_SEGID = undef;
	my $SINK_SEGID = undef;
	my $outstr_HOSTSEG = '';
	for my $Hseg_NO ( sort {$a<=>$b} keys %{$V_Href->{dHG_HsegInfo}} ){
		my $this_Segm_Href = $V_Href->{dHG_HsegInfo}->{$Hseg_NO};
		# get the mean depth of this segemnt in pure Tumor cell
		$this_Segm_Href->{puT_meanDepth} = sprintf "%.2f",
			get_givenWinItvalMeanDepth(
										winDepth_Href => $V_Href->{dG_winDepth},
										tissue => $V_Href->{dHG_TissueToShowCN},
										depth_key => 'mean_adjdepth',
										st_winNO => $this_Segm_Href->{st_winNO},
										ed_winNO => $this_Segm_Href->{ed_winNO}
									);
		# get the copy number of this host-segment in pure Tumor cell
		if( defined $V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio} ){
			$this_Segm_Href->{puT_copyNumber} = sprintf "%.2f", $this_Segm_Href->{puT_meanDepth} * $V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio};
		}
		# inform
		# SEG     ID=H1;CN_ORIG=3.25;W=1;INTERVAL=chr8:128255000-128277177;PTDP=59.17
		# PTDP: Pure Tumor-cell DePth
		my $Segm_ID_toShow = ($this_Segm_Href->{Segm_ID} =~ /^(.*[^\d\:])(\d+)$/) ? "$1$2" : $this_Segm_Href->{Segm_ID};
		my @output = (
						'SEG',
						join( ';', (
									 'ID='.$Segm_ID_toShow,
									 'CN_ORIG='.($this_Segm_Href->{puT_copyNumber} || 'N/A'),
									 'W=1',
									 'INTERVAL='.$V_Href->{dHG_draw_region}->{refseg}.':'.$this_Segm_Href->{st_pos}.'-'.$this_Segm_Href->{ed_pos},
									 'PTDP='.$this_Segm_Href->{puT_meanDepth}
							  		)
							)
					 );
		$outstr_HOSTSEG .= join("\t", @output)."\n";
		# source or sink
		$SOURCE_SEGID = $Segm_ID_toShow if( !defined $SOURCE_SEGID );
		$SINK_SEGID = $Segm_ID_toShow; # always update to last one
	}

	# only in somatic mode, deal with virus-segments and junctions
	my @outstr_FVGM = ();
	my $outstr_VSEG = '';
	my $outstr_JUNC = '';
	if( $V_Href->{dG_CaseVsCtrl_mode} ){
		# calculate the copy number of virus-segments, and output
		for my $Vseg_NO ( sort {$a<=>$b} keys %{$V_Href->{dHG_VsegInfo}} ){
			my $this_vSegm_Href = $V_Href->{dHG_VsegInfo}->{$Vseg_NO};
			# get the copy number of this virus-segment in pure Tumor cell
			if( defined $V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio} ){
				$this_vSegm_Href->{puT_copyNumber} = sprintf "%.2f", $this_vSegm_Href->{puT_meanDepth} * $V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio};
			}
			# inform
			# SEG     ID=V1;CN_ORIG=329.75;W=1;INTERVAL=HPV18:1-406;RAWDP=6000
			my @output = (
							'SEG',
							join( ';', ( 
										 'ID='.$this_vSegm_Href->{Segm_ID},
										 'CN_ORIG='.($this_vSegm_Href->{puT_copyNumber} || 'N/A'),
										 'W=1',
										 'INTERVAL='.$this_vSegm_Href->{vrefseg}.':'.$this_vSegm_Href->{st_pos}.'-'.$this_vSegm_Href->{ed_pos},
										 'RAWDP='.$this_vSegm_Href->{puT_meanDepth}
								  		)
								)
						 );
			$outstr_VSEG .= join("\t", @output)."\n";
			# whole virus genome string
			push @outstr_FVGM, $this_vSegm_Href->{Segm_ID}.'+';
		}

		# calculate the copy number of junctions, and output
		for my $junc_NO ( sort {$a<=>$b} keys %{$V_Href->{dHG_JuncInfo}} ){
			my $this_junc_Href = $V_Href->{dHG_JuncInfo}->{$junc_NO};
			# transform the raw junction reads to pure Tumor cell
			$this_junc_Href->{puT_JRnum} = int
				get_ObjCountOfPureTumorCell(
												ObjCount => $this_junc_Href->{raw_JRnum},
												MixedCellMeanDepth => $V_Href->{dHG_meandepth}->{case},
												Ncell_GMratio => $V_Href->{case_Ncell_GMratio},
												ObjCopyRatioInNcell => $this_junc_Href->{N_copyratio}
											);
			# get the copy number of this junction in pure Tumor cell
			$this_junc_Href->{puT_copyNumber} = sprintf "%.2f", $this_junc_Href->{puT_JRnum} * $V_Href->{dHG_CNbaseOnReg}->{depth2cn_ratio};
			# inform
			# JUNC    ID=JUNC1;SEGLINK=V1+,H4-;CN_ORIG=2.47;W=1;REFSEG_5P=HPV18;BKPOS_5P=406;STRD_5P=+;REFSEG_3P=chr8;BKPOS_3P=128283649;STRD_3P=-;NCPR=0;RAWJRC=45;PTJRC=45
			my $segID_5p_toShow = ($this_junc_Href->{segID_5p} =~ /^(.*[^\d\:])(\d+)$/) ? "$1$2" : $this_junc_Href->{segID_5p};
			my $segID_3p_toShow = ($this_junc_Href->{segID_3p} =~ /^(.*[^\d\:])(\d+)$/) ? "$1$2" : $this_junc_Href->{segID_3p};
			my @output = (
							'JUNC',
							join( ';', ( 
										 'ID=JUNC'.$junc_NO,
										 'SEGLINK='.$segID_5p_toShow.$this_junc_Href->{strd_5p}.','.$segID_3p_toShow.$this_junc_Href->{strd_3p},
										 'CN_ORIG='.($this_junc_Href->{puT_copyNumber} || 'N/A'),
										 'W=1',
										 'REFSEG_5P='.$this_junc_Href->{refseg_5p},
										 'BKPOS_5P='.$this_junc_Href->{bkpos_5p},
										 'STRD_5P='.$this_junc_Href->{strd_5p},
										 'REFSEG_3P='.$this_junc_Href->{refseg_3p},
										 'BKPOS_3P='.$this_junc_Href->{bkpos_3p},
										 'STRD_3P='.$this_junc_Href->{strd_3p},
										 'NCPR='.$this_junc_Href->{N_copyratio},
										 'RAWJRC='.$this_junc_Href->{raw_JRnum},
										 'PTJRC='.$this_junc_Href->{puT_JRnum}
								  		)
								)
						 );
			$outstr_JUNC .= join("\t", @output)."\n";
		}
	}

	# prepare the input file to haplotype algorithm
	open (FORHAP, Try_GZ_Write($V_Href->{dHG_InputForLocalHapAlgorithm})) || die "fail write input-file to haplotype algorithm: $!\n";
	print FORHAP "##VERSION=v1.1\n";
	print FORHAP "##_please_find_HEADER_in_online_demo_\n";
	# sample ID
	print FORHAP join( "\t", 'SAMPLE', $V_Href->{dHG_sampleIDtoShow} || '_please_fix_me_' )."\n";
	# ploidy from patchwork for major/minor Local Map
	print FORHAP join( "\t", 'PLOIDY', '_please_fix_me_' )."\n";
	# source segment
	print FORHAP join( "\t", 'SOURCE', $SOURCE_SEGID || '_please_fix_me_' )."\n";
	# sink segment
	print FORHAP join( "\t", 'SINK', $SINK_SEGID || '_please_fix_me_' )."\n";
	# free virus genome
	print FORHAP join( "\t", 'FVGM', 'ID=FVGM1;SEGSTR=' . (join(',',@outstr_FVGM) || '_please_fix_me_') )."\n";
	# host segments
	print FORHAP $outstr_HOSTSEG;
	# virus segments
	print FORHAP $outstr_VSEG;
	# junctions
	print FORHAP $outstr_JUNC;
	close FORHAP;
}

#--- show breakpoints ---
#-- try to show segments
#-- try to show CN of segment
sub show_breakpoints{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# has any avaiable breakpoint?
	return if( scalar(keys %{$V_Href->{dG_bkWindow}}) == 0 );

	# prepare for the copy number informtion of host-segments
	&calculate_CopyNumber;

	# segmentation
	# attributes
	my $SegmAxisGap = 10;
	my $SegmHeight = 10;
	my $Segm_Y = $V_Href->{general_BaseY} + $SegmAxisGap + $SegmHeight / 2;
	# get all Segments's edge
	my %SegmEdgePos = map { ($V_Href->{dHG_draw_region}->{$_}, 1) } qw/ st_pos ed_pos /;
	$SegmEdgePos{$_} = 1 for values %{$V_Href->{dG_bkWindow}};
	my @SegmEdgePos = sort {$a<=>$b} keys %SegmEdgePos;
	# show segmentations
	my $last_SegmEdgeWinNO = 1;
	for my $Hseg_NO ( sort {$a<=>$b} keys %{$V_Href->{dHG_HsegInfo}} ){
		my $this_Segm_Href = $V_Href->{dHG_HsegInfo}->{$Hseg_NO};
		my $st_winNO = $this_Segm_Href->{st_winNO};
		my $ed_winNO = $this_Segm_Href->{ed_winNO};
		my $Segm_ID  = $this_Segm_Href->{Segm_ID};
		my $puT_avgDP= $this_Segm_Href->{puT_meanDepth};
		my $puT_CN   = $this_Segm_Href->{puT_copyNumber};
		# show this segment body
		my $this_SegmSpanWinNum = $ed_winNO - $st_winNO + 1;
		my $this_Segm_X = $V_Href->{dG_axisZX} + ($st_winNO + $ed_winNO) / 2;
		draw_a_parallelogram(
								\$V_Href->{general_SVG_obj},
								x => $this_Segm_X,
								y => $Segm_Y,
								fill_color => $V_Href->{dHG_SegmFillCol},
								head_bottom_side_len => $this_SegmSpanWinNum,
								left_right_side_len => $SegmHeight,
								boundary_color => 'black',
								text_fill => $Segm_ID,
								text_font_fam => $V_Href->{general_FontFam_Bold},
								text_font_size => $V_Href->{general_FontSize},
								font_size_auto_adjust => 1,
								font_size_auto_adjust_wRatio => 0.8,
								font_size_auto_adjust_hRatio => 0.8
							);
		# show segment's copynumber
		if( defined $puT_CN ){
			## show CN line
			my $cnLine_X1 = $V_Href->{dG_axisZX} + $st_winNO;
			my $cnLine_X2 = $V_Href->{dG_axisZX} + $ed_winNO;
			my $cnLine_Y  = $V_Href->{dG_DepthSpecBtom_Y} - $puT_avgDP / $V_Href->{dG_y_resol};
			$V_Href->{general_SVG_obj}
				   ->line(
							x1 => $cnLine_X1,
							y1 => $cnLine_Y,
							x2 => $cnLine_X2,
							y2 => $cnLine_Y,
							stroke => $V_Href->{dHG_CNlineColor},
							'stroke-width' => $V_Href->{dHG_CNlineWidth},
							'stroke-linecap' => 'round'
						 );
			## show CN text
			my $cnText_X = ( $cnLine_X1 + $cnLine_X2 ) / 2;
			my $cnText_Y = $cnLine_Y - $V_Href->{dHG_CNlineWidth} * 2;
			show_text_in_line(
								\$V_Href->{general_SVG_obj},
								text_x => $cnText_X,
								text_y => $cnText_Y,
								text => "CN~$puT_CN",
								font_family => $V_Href->{general_FontFam},
								font_size => $V_Href->{general_FontSize},
								text_col => $V_Href->{dHG_CNlineColor},
								text_anchor => 'middle',
								height_adjust => 0,
								width_limit => $this_SegmSpanWinNum * 0.9,
								width_limit_fontsizeTorotate => -6
							);
		}
	}
	# update
	$V_Href->{general_BaseY} = $Segm_Y + $SegmHeight / 2;
	# inform
	stout_and_sterr `date`."[INFO]\tDraw segments according to breakpoints OK.\n";

	# show breakpoints
	my $bkPosTextAxisGap = 10;
	my $bkPosTextGap = 20;
	my $bkPosTextHeight = 0;
	my $bkPosTextLayerGap = 8;
	my $bkPosTextFtfm = $V_Href->{general_FontFam_Bold};
	my $bkPosTextFtsz = $V_Href->{general_FontSize};
	my $bkPosTextLineWidth = 1;
	my $bkPosTextLine_topY = $V_Href->{dG_DepthSpecBtom_Y};
	my $bkPosTextCol = 'red';
	# get bkPos Text Span Interval
	#- bkPosText_Href structure:
	#-- bkPosText_Href -> $bkPos -> 'body_drawItval' -> [[a,b], [c,d]]
	#-- bkPosText_Href -> $bkPos -> 'show_layer'
	my $bkPosText_Href = {};
	for my $bkPosWNO (sort {$a<=>$b} keys %{$V_Href->{dG_bkWindow}}){
		my $bkPos = $V_Href->{dG_bkWindow}->{$bkPosWNO};
		# text size
		my ($text_height, $text_width) = 
			get_size_of_text_to_show(
										font_family => $bkPosTextFtfm,
										font_size => $bkPosTextFtsz,
										text => $bkPos,
										return_HeightWid => 1
									);
		$bkPosText_Href->{$bkPos}->{text_size} = [ $text_height, $text_width ];
		# update
		$bkPosTextHeight = max( $bkPosTextHeight, $text_height );
		# recode bkpos reason
		$bkPosText_Href->{$bkPos}->{bk_reason} = $V_Href->{dG_bkReason}->{$bkPosWNO};
		# record Text Span Interval
		my $bkPosText_lftPix = $V_Href->{dG_axisZX} + $bkPosWNO - $text_width / 2;
		my $bkPosText_rgtPix = $V_Href->{dG_axisZX} + $bkPosWNO + $text_width / 2;
		$bkPosText_Href->{$bkPos}->{body_drawItval} = [ [ $bkPosText_lftPix, $bkPosText_rgtPix ] ];
	}
	# arrange bkPosText
	my $max_show_layerNO = 0;
	allocate_object_vertical(
								obj_Href => $bkPosText_Href,
								obj_itval_key => 'body_drawItval',
								obj_layer_key => 'show_layer',
								obj_gap => $bkPosTextGap,
								max_layer_Sref => \$max_show_layerNO
							);
	# display breakpoints
	for my $bkPos (sort {$bkPosText_Href->{$b}->{show_layer}<=>$bkPosText_Href->{$a}->{show_layer}} keys %$bkPosText_Href){
		my $bkPosTextLayerNO = $bkPosText_Href->{$bkPos}->{show_layer};
		my $bkPosText_X = sum( @{ $bkPosText_Href->{$bkPos}->{body_drawItval}->[0] } ) / 2;
		my $bkPosText_Y = $V_Href->{general_BaseY} + $bkPosTextAxisGap + ($bkPosTextHeight + $bkPosTextLayerGap) * ($bkPosTextLayerNO - 1);
		# show mark line
		my @line_para = (stroke => $bkPosTextCol, 'stroke-width' => $bkPosTextLineWidth, 'stroke-linecap' => 'round');
		## dash-line for SV whose both break-points locate in this draw region
		push @line_para, ('stroke-dasharray' => '5,2') if( $bkPosText_Href->{$bkPos}->{bk_reason} =~ /intra/ );
		$V_Href->{general_SVG_obj}
			   ->line(
						x1 => $bkPosText_X,
						y1 => $bkPosTextLine_topY,
						x2 => $bkPosText_X,
						y2 => $bkPosText_Y,
						@line_para
					 );
		# breakpoint text
		my ($text_height, $text_width) = @{ $bkPosText_Href->{$bkPos}->{text_size} };
		## show text square
		draw_a_parallelogram(
								\$V_Href->{general_SVG_obj},
								x => $bkPosText_X,
								y => $bkPosText_Y + $text_height / 2,
								fill_color => 'white',
								opacity => 0.618,
								head_bottom_side_len => $text_width * 0.85,
								left_right_side_len => $text_height,
								boundary_color => $bkPosTextCol,
								text_fill => $bkPos,
								text_col => $bkPosTextCol,
								text_font_fam => $bkPosTextFtfm,
								text_font_size => $bkPosTextFtsz,
								font_size_auto_adjust => 1,
								font_size_auto_adjust_wRatio => 1,
								font_size_auto_adjust_hRatio => 0.8
							);
	}

	# inform
	stout_and_sterr `date`."[INFO]\tDisplay breakpoints OK.\n";
}

#--- show basic figure information ---
# use a Golden rectangle
# figure name, version, author and institute
sub show_figure_information{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# basic values
	my $golden_ratio = ( sqrt(5) + 1 ) / 2; # (sqrt(5)+1)/2 = 1.6180339887498948
	my $botrig_y = min($V_Href->{general_BaseY} / 2, 80);
	my $botrig_x = $botrig_y * $golden_ratio;

	# background
	draw_a_parallelogram(
							\$V_Href->{general_SVG_obj},
							x => $botrig_x/2,
							y => $botrig_y/2,
							head_bottom_side_len => $botrig_x,
							left_right_side_len => $botrig_y,
							boud_dasharray => '5,5',
							fill_color => 'aliceblue',
							usage_print => 0
						);
	# contents
	my @contents = (
					["Host Local Genome", 'Times Bold', 1],
					["[------- v$VERSION -------]", 'Times', 2, {font_family=>{'Times Italic'=>[10 .. 9+length($VERSION)+1]}}],
					["author: $AUTHOR", 'Times', 2],
					["CityU of Hong Kong", 'Times', 2]
				   );
	my $text_left_x = $botrig_x * 0.1;
	my $text_height = $botrig_y / scalar(@contents);
	my $text_left_y = $text_height / 2 * 1.1;
	for my $content_Aref (@contents){
		show_text_in_line(
							\$V_Href->{general_SVG_obj},
							text_x => $text_left_x,
							text_y => $text_left_y,
							text => $content_Aref->[0],
							font_family => $content_Aref->[1],
							font_size => 12,
							text_col => 'black',
							text_anchor => 'start',
							height_adjust => 1,
							width_limit => $botrig_x - $text_left_x * $content_Aref->[2],
							height_limit => $text_height,
							features => ($content_Aref->[3] || {})
						);
		$text_left_y += $text_height;
	}
}

#--- generate figure file ---
sub generate_figure_output{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# file type (based on postfix)
	my ($figure_postfix) = ($V_Href->{dG_figure} =~ /\.([sp][vnd][gf])$/i);
	# temp svg
	my $out_svg = "$V_Href->{dG_figure}.svg";
	open (SVGO,">$out_svg") || die"fail write temp SVG figure: $!\n";
	print SVGO $V_Href->{general_SVG_obj}->xmlify;
	close SVGO;
	# user required file
	if( $figure_postfix =~ /svg/i ){
		if( $out_svg ne $V_Href->{dG_figure} ){
			`mv -f $out_svg $V_Href->{dG_figure}`;
		}
	}
	else{
		# https://stackoverflow.com/questions/9853325/how-to-convert-a-svg-to-a-png-with-image-magick
		# https://www.imagemagick.org/script/command-line-options.php
		`convert -density 1000 -resize 1500x1500 $out_svg $V_Href->{dG_figure}`;
		`rm -rf $out_svg`;
	}

	# inform
	stout_and_sterr `date`."[INFO]\tGenerate final figure OK.\n";
}

#--- get segment ID based on naming mode ---
sub get_HostSegmentID{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $Segm_NO = $parm{Segm_NO};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	if( $V_Href->{dHG_SegmentMode} =~ /^u(:[A-Z]{0,1})?$/ ){
		(my $stIdx = defined $1 ? $1 : '') =~ s/^://;
		$stIdx = ($stIdx =~ /^[A-Z]$/) ? $stIdx : 'A';
		my $segID = chr( ord($stIdx) + $Segm_NO - 1 );
		if( $segID !~ /^[A-Z]$/ ){
			warn_and_exit `date`."<ERROR>\tName of segment NO.$Segm_NO is not from [A-Z].\n";
		}
		return $segID;
	}
	elsif( $V_Href->{dHG_SegmentMode} =~ /^l(:[a-z]{0,1})?$/ ){
		(my $stIdx = defined $1 ? $1 : '') =~ s/^://;
		$stIdx = ($stIdx =~ /^[a-z]$/) ? $stIdx : 'a';
		my $segID = chr( ord($stIdx) + $Segm_NO - 1 );
		if( $segID !~ /^[a-z]$/ ){
			warn_and_exit `date`."<ERROR>\tName of segment NO.$Segm_NO is not from [a-z].\n";
		}
		return $segID;
	}
	elsif( $V_Href->{dHG_SegmentMode} =~ /^(\w+?)(:\d*)?$/ ){
		my $prefix = $1;
		(my $stIdx = defined $2 ? $2 : '') =~ s/^://;
		$stIdx = ($stIdx =~ /^\d+$/) ? $stIdx : 1;
		return $prefix . ( $stIdx + $Segm_NO - 1 );
	}
	else{
		warn_and_exit `date`."<ERROR>\tCannot recognize the Segment-Name mode ($V_Href->{dHG_SegmentMode}).\n";
	}
}

#--- get depth of window, tranform to adjdepth ---
sub read_depthStat_to_window_adjdepth{
	
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $depthStatFile = $parm{depthStatFile};
	my $tissue = $parm{tissue};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# for ctrl, OK-pos stat if set {dHG_ctrl_depthFluctR}
	my $ctrlSkipPosCount = 0;
	my $ctrlKeepPosCount = 0;
	my ($ctrl_copyratio, $ctrl_copyR_last_edp) = (1, 0); # for known CNV in ctrl

	my $last_window_NO = 0;
	my $last_window_last_pos = 0;
	my $last_window_resol = 0;
	open (DEPTH, Try_GZ_Read($depthStatFile)) || die "fail load depthStatFile: $!\n";
	while (<DEPTH>){
		my @info = split;
		my $refseg = shift @info;
		next if( $refseg ne $V_Href->{dHG_draw_region}->{refseg} );
		my $pos = shift @info;
		my $window_NO = $last_window_NO; # initialize as last window
		if( $pos > $last_window_last_pos ){
			my @winInfo = get_pos_windowNO(
											pos => $pos,
											@{$V_Href->{dG_PosToWinNOSets}},
											give_winLastPos => 1,
											winResol_Sref => \$last_window_resol
										  );
			my $sameWinLastPos = $winInfo[0];
			$window_NO = $winInfo[1];
			# deal with last window
			if( $last_window_NO != $window_NO ){
				if( exists $V_Href->{dG_winDepth}->{$tissue}->{$window_NO} ){
					warn_and_exit "<ERRO>\ttissue $tissue window NO.$window_NO ($pos) shows again.\n";
				}
				my $selected_okPos_Href = (    $V_Href->{dHG_ctrl_depthFluctR_mode}
											&& $last_window_resol >= $V_Href->{dHG_ctrl_depthFluctR_minEffectRegResol}
										  ) ? $V_Href->{CtrlPos2Depth} : undef;
				deal_window_depth_info(
										winDepth_Href => $V_Href->{dG_winDepth},
										window_NO => $last_window_NO,
										tissue => $tissue,
										role => 'host',
										selected_okPos_Href => $selected_okPos_Href,
										Eng_3SD => $V_Href->{dG_DepEng_3SD},
										TmmRatio => $V_Href->{dG_DepTrimRatio},
										HasT_mode => $V_Href->{dG_HasT_mode},
										h_DepthToDraw => $V_Href->{DepthToDraw},
										h_meanDepth_Href => $V_Href->{dHG_meandepth},
										h_InCase_Tcell_GMratio => $V_Href->{case_Tcell_GMratio},
										h_ctrl_prior_copyR_Aref => $V_Href->{dHG_ctrl_prior_copyR},
										h_InCase_Ncell_GMratio => $V_Href->{case_Ncell_GMratio}
									);
			}
			# update
			$last_window_NO = $window_NO;
			$last_window_last_pos = $sameWinLastPos;
		}
		# store depth info
		my $this_pos_depth = sum(@info);
		$V_Href->{dG_winDepth}->{$tissue}->{$window_NO}->{pos2depth}->{$pos} = $this_pos_depth;
		# control sample, for position selection based on depth fluctuation
		# basically, ctrl sample depth load first for this function, then ctrl and case could both use
		if(    $V_Href->{dHG_ctrl_depthFluctR_mode}
			&& $tissue eq 'ctrl'
		){
			# get known CNV in ctrl
			if( $pos > $ctrl_copyR_last_edp ){
				($ctrl_copyratio, $ctrl_copyR_last_edp) = get_ctrl_copyR( ctrl_prior_copyR_Aref=>$V_Href->{dHG_ctrl_prior_copyR}, pos=>$pos);
			}
			# test depth of this pos
			if(    $this_pos_depth < $V_Href->{dHG_ctrl_depthFluctR_flor} * $ctrl_copyratio
				|| $this_pos_depth > $V_Href->{dHG_ctrl_depthFluctR_ceil} * $ctrl_copyratio
			){ # out of OK-ctrl-depth-range, skip
				$ctrlSkipPosCount ++;
			}
			else{ # keep it
				$V_Href->{CtrlPos2Depth}->{$pos} = $this_pos_depth;
				$ctrlKeepPosCount ++;
			}
		}
	}
	close DEPTH;
	# inform 'ctrl pos selection' if set
	if(    $V_Href->{dHG_ctrl_depthFluctR_mode}
		&& $tissue eq 'ctrl'
	){
		stout_and_sterr `date`."[INFO]\tTotally, based on ctrl_depthFluctR ($V_Href->{dHG_ctrl_depthFluctR}), in ctrl sample:\n"
									."\t$ctrlSkipPosCount pos were discarded, $ctrlKeepPosCount pos were kept for depth-stat in both case and ctrl.\n";
	}
	# very last one window
	my $selected_okPos_Href = (    $V_Href->{dHG_ctrl_depthFluctR_mode}
								&& $last_window_resol >= $V_Href->{dHG_ctrl_depthFluctR_minEffectRegResol}
							  ) ? $V_Href->{CtrlPos2Depth} : undef;
	deal_window_depth_info(
							winDepth_Href => $V_Href->{dG_winDepth},
							window_NO => $last_window_NO,
							tissue => $tissue,
							role => 'host',
							selected_okPos_Href => $selected_okPos_Href,
							Eng_3SD => $V_Href->{dG_DepEng_3SD},
							TmmRatio => $V_Href->{dG_DepTrimRatio},
							HasT_mode => $V_Href->{dG_HasT_mode},
							h_DepthToDraw => $V_Href->{DepthToDraw},
							h_meanDepth_Href => $V_Href->{dHG_meandepth},
							h_InCase_Tcell_GMratio => $V_Href->{case_Tcell_GMratio},
							h_ctrl_prior_copyR_Aref => $V_Href->{dHG_ctrl_prior_copyR},
							h_InCase_Ncell_GMratio => $V_Href->{case_Ncell_GMratio}
						);
}

#--- 
1; ## tell the perl script the successful access of this module.
