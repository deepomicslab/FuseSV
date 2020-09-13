package FuseSV::Virus_Type::DrawViralGenome;

use strict;
use warnings;
use File::Basename qw/ basename dirname /;
use Cwd qw/ abs_path /;
use List::Util qw/ min max sum /;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::SVG_Radian_System_Elements qw/ %COLOR_DB draw_circle_seg /;
use SOAPfuse::SVG_Font qw/ show_text_in_line get_size_of_text_to_show /;
use SOAPfuse::SVG_Orthogonal_System_Elements qw/ draw_a_parallelogram draw_a_triangle /;
use FuseSV::LoadOn;
use FuseSV::GetPath;
use FuseSV::Virus_Type::VirusTypeMain;
use FuseSV::Visual::ArrangeObj qw/ get_pos_windowNO allocate_object_bilater /;
use FuseSV::Visual::DepthWork qw/ get_windowSmoDepth deal_window_depth_info draw_depth_spectrum /;
use FuseSV::Visual::GeneWork qw/ load_transInfo_from_tpsl add_user_region_to_gene_info allocate_genes_and_draw /;
use FuseSV::Visual::AxisWork qw/ set_Y_axis_Resol_LabStep show_Y_axis show_X_axis /;
use SVG;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  viral_genome_draw
			  show_figure_information
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'DrawViralGenome';
#----- version --------
$VERSION = "0.15";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						draw_viral_genome
						load_virus_depth
						load_virus_breakpoint
						load_viral_mutation
						obtain_mut_span_width
						arrange_mut_icon_location
						smooth_window_depth
						load_viral_trans_info
						set_y_axis
						draw_y_axis
						draw_mutations
						draw_genes
						draw_depth
						draw_x_axis
						generate_figure_output
						show_figure_information
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} draw_virus <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -vpf    [s]  virus database folder. <required>
	    -vtype  [s]  virus type. <required>
	                  it should be found at the column 'subtype' in file 'virus.info' in virus database folder.
	    -depth  [s]  the depth stat file of virus variant, provided by 'virus_depth' func. <required>
	    -o      [s]  output figure file, allowed format: svg/png/pdf. <required>

	   # Annotations #
	    -id     [s]  sample ID to show in the figure. [Optional]
	    -mutf   [s]  the modifications list of viral variant, provided by 'virus_type' func. [Optional]
	    -vbkl   [s]  break-points list of variant, such as virus integration sites, rearrangements. [Optional]
	                  TAB-delimited format:  vRef  VBkPos  Label
	                  Note: 'vRef' must match that in the depth stat file.
	    -glbi   [s]  use which block as gene label. [g]
	                  allowed options: 't': trans_name, 'g': gene_name, 'l': locus, 'n': note, 'p': protein_product.
	                  Note: 'locus_tag', 'note', 'protein_product' are from NCBI GenBank webpage text in virus database.
	    -gskip  [s]  skip virus gene (column NO.10 in trans.psl), allow use in multiple times. [none]
	    -no_gene     do not display viral genes. [Disabled]
	    -nreg   [s]  specific region(s) on '-vtype' to display. [NULL]
	                  instance: -nreg 'TEST1:100-700' -nreg 'TEST2:300-600,900-1500'

	   # Options about Depth #
	    -tmrt   [f]  ratio to trim for even depth of one region. [0]
	                  positive value: bilateral operation; negative: only tail.
	    -eng3sd      apply engineer_3times_SD_filter for average depth. [Disabled]
	                  Once enabled, '-tmrt' will be ignored.
	    -fwsm   [i]  how many bilateral flanking windows to make curve smooth. [2]
	                  note: value '2' means 2(fore)+1(this)+2(after); set '0' to disable.
	    -src    [i]  sequenced reads count of this sample for equivalent compare. [NULL]
	    -nsrc   [i]  sequenced reads count that normalize to for equivalent compare. [NULL]

	   # Options about Display #
	   ## colour
	    -sp_col [s]  the colour of depth spectrum. [red]
	    -bgbcol [s]  the colour of background boundary. [black]
	   ## region
	    -stpos  [i]  the minimum viral genomic position to show. [auto]
	    -edpos  [i]  the maximum viral genomic position to show. [auto]
	   ## axis
	    -x_res  [i]  how many bases (nt) one spot in figure represents on x-axis. [10]
	    -x_lnb  [i]  label number basement shown on x-axis. [100]
	    -x_lbs  [i]  how many screen pixels to show one position label on x-axis. [200]
	    -x_lbfz [i]  font size of label on x-axis. [12]
	    -y_len  [i]  length of y-axis, will be overwritten when '-y_res' is set. [80]
	    -y_res  [i]  how many depth (X) one spot in figure represents on y-axis. [auto]
	    -y_lbs  [i]  how many screen pixels to show one depth marker on y-axis, effective with '-y_res'. [auto]
	    -y_lbfz [i]  font size of label on y-axis. [12]
	   ## others
	    -sp_leg [s]  the text in legend of depth spectrum. [Optional]
	    -extwd  [i]  pixels to extend the width  for enough space to display. [minimum:150]
	    -extht  [i]  pixels to extend the height for enough space to display. [minimum:150]

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
			[ dVG_viral_subtype => '' ],
			[ dVG_DepFile => { case=>'', ctrl=>'' } ],
			[ dG_figure => '' ],
			# annotations
			[ dVG_viral_MutFile => undef ],
			[ dVG_viral_BkList => undef ],
			[ dVG_show_no_gene => 0 ],
			[ dVG_gene_label_idx => 'g' ],
			[ dVG_gene_to_skip => [] ],
			[ dG_user_gene => [] ],
			[ dVG_sampleIDtoShow => '' ],
			# depth options
			[ dG_DepTrimRatio => 0 ],
			[ dG_DepEng_3SD => 0 ],
			[ dG_DepSmooth_flkWd => 2 ],
			[ dVG_SeqReadsCount => { case=>0, ctrl=>0 } ],
			[ dVG_case2norm_seqRatio => 1 ], # default to no change
			# draw options
			[ dG_DepthSpecColor => { case=>'red', ctrl=>'dodgerblue' } ],
			[ dVG_spectrum_icon_text => undef ],
			[ dG_extWidth => 150 ],
			[ dG_extHeight => 150 ],
			[ dVG_WdSize => 10 ],
			[ dG_axisZX => 0 ],
			[ dG_axisZY => 0 ],
			[ dG_x_axisLen => 0 ], # auto assign
			[ dG_x_LabBaseNum => 100 ],
			[ dG_x_LabStep => 200 ],
			[ dG_x_LabFtsz => 12 ],
			[ dG_y_axisLen => 80 ],
			[ dG_y_resol => 0 ], # auto adjust
			[ dG_y_LabStep => 0 ], # auto adjust
			[ dG_y_LabFtsz => 12 ],
			[ dVG_pos_start => undef ],
			[ dVG_pos_end => undef ],
			[ dVG_mutLabelFontSize => 10 ],
			[ dVG_mutIconEdgeLen => 6 ],
			[ dVG_mutIconGap => undef ],
			# virus
			[ dVG_viral_orig_tpsl => '' ], # will auto detect
			[ dVG_viral_origLen => 0 ],
			# colour
			[ dVG_viral_mut_icon => {	 'snv' => { 'col' => 'green',  'icon' => 'circle'       },
										 'snp' => { 'col' => 'green',  'icon' => 'circle'       },
								 		 'ins' => { 'col' => 'yellow', 'icon' => 'rev_triangle' },
								 		 'del' => { 'col' => 'yellow', 'icon' => 'triangle'     },
									'integ_sv' => { 'col' => 'red',    'icon' => 'parallelogram'}   }  ],
			# intermediate variants
			[ dVG_vRef_to_draw => '' ],
			[ dG_maxDepthToDraw => 0 ],
			[ dG_maxWinNOtoDraw => 0 ],
			# check regResol_Href structure in FuseSV::Visual::ArrangeObj module
			[ dG_regResol => {} ],
			# check winDepth_Href structure in FuseSV::Visual::DepthWork module
			[ dG_winDepth => {} ],
			# dG_bkWindow -> $window_NO = $bkpos
			[ dG_bkWindow => {} ],
			# dVG_ViralMut -> $mut_count -> 'pos'
			# dVG_ViralMut -> $mut_count -> 'wNO'
			# dVG_ViralMut -> $mut_count -> 'show_wNO'
			# dVG_ViralMut -> $mut_count -> 'mut_type'
			# dVG_ViralMut -> $mut_count -> 'ref_allele'
			# dVG_ViralMut -> $mut_count -> 'mut_allele'
			# dVG_ViralMut -> $mut_count -> 'mut_label'
			# dVG_ViralMut -> $mut_count -> 'label_height'
			# dVG_ViralMut -> $mut_count -> 'mut_width'
			[ dVG_ViralMut => {} ],
			# check TransInfo_Href structure in FuseSV::Visual::GeneWork module
			[ dG_TransInfo => {} ],
			[ dVG_x_maxPos => -1 ],
			[ dVG_x_minPos => 1E10 ],
			[ dG_bgStrokeCol => 'black' ],
			[ dG_DepthSpecBtom_Y => 0 ],

            #--- containers ---#
            # id(subtype) -> virusOB
            # see FuseSV::Objects::Virus_OB
            [ virusRefSeg => {} ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['virus_db_folder'],
								  ['dVG_DepFile', 'case'],
								  ['dVG_DepFile', 'ctrl'],
								  ['dG_figure'],
								  ['dVG_viral_MutFile'],
								  ['dVG_viral_BkList']     ] ]
		);

	# load containers
	FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool( only_load_containers => 1 );
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-o:s"		=> \$V_Href->{dG_figure},
		"-depth:s"	=> \$V_Href->{dVG_DepFile}->{case},
		# "-N_df:s"	=> \$V_Href->{dVG_DepFile}->{ctrl},
		"-vpf:s"	=> \$V_Href->{virus_db_folder},
		"-vtype:s"	=> \$V_Href->{dVG_viral_subtype},
		# annotation
		"-id:s"		=> \$V_Href->{dVG_sampleIDtoShow},
		"-mutf:s"	=> \$V_Href->{dVG_viral_MutFile},
		"-vbkl:s"	=> \$V_Href->{dVG_viral_BkList},
		"-glbi:s"	=> \$V_Href->{dVG_gene_label_idx},
		"-gskip:s"	=> \@{ $V_Href->{dVG_gene_to_skip} },
		"-no_gene"	=> \$V_Href->{dVG_show_no_gene},
		"-nreg:s"	=> \@{$V_Href->{dG_user_gene}},
		# depth operation
		"-tmrt:f"	=> \$V_Href->{dG_DepTrimRatio},
		"-eng3sd"	=> \$V_Href->{dG_DepEng_3SD},
		"-fwsm:i"	=> \$V_Href->{dG_DepSmooth_flkWd},
		"-src:i"	=> \$V_Href->{dVG_SeqReadsCount}->{case},
		"-nsrc:i"	=> \$V_Href->{dVG_SeqReadsCount}->{ctrl},
		"-stpos:i"	=> \$V_Href->{dVG_pos_start},
		"-edpos:i"	=> \$V_Href->{dVG_pos_end},
		# draw options
		"-sp_col:s"	=> \$V_Href->{dG_DepthSpecColor}->{case},
		"-sp_leg:s"	=> \$V_Href->{dVG_spectrum_icon_text},
		"-bgbcol:s"	=> \$V_Href->{dG_bgStrokeCol},
		# "-N_col:s"	=> \$V_Href->{dG_DepthSpecColor}->{ctrl},
		"-x_res:i"	=> \$V_Href->{dVG_WdSize},
		"-x_lnb:i"	=> \$V_Href->{dG_x_LabBaseNum},
		"-x_lbs:i"	=> \$V_Href->{dG_x_LabStep},
		"-x_lbfz:i"	=> \$V_Href->{dG_x_LabFtsz},
		"-y_len:i"	=> \$V_Href->{dG_y_axisLen},
		"-y_res:f"	=> \$V_Href->{dG_y_resol},
		"-y_lbs:f"	=> \$V_Href->{dG_y_LabStep},
		"-y_lbfz:i"	=> \$V_Href->{dG_y_LabFtsz},
		"-extwd:i"	=> \$V_Href->{dG_extWidth},
		"-extht:i"	=> \$V_Href->{dG_extHeight},
		"-icgap:f"	=> \$V_Href->{dVG_mutIconGap}, # hidden option
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

	return  (  $V_Href->{HELP}
		     || !$V_Href->{dVG_viral_subtype}
			 || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder})
			 || (!$V_Href->{dVG_DepFile}->{case} || !-e $V_Href->{dVG_DepFile}->{case})
			 || $V_Href->{dG_figure} !~ /\.([sp][vnd][gf])$/i
			 || $V_Href->{dG_DepSmooth_flkWd} < 0
			 || $V_Href->{dVG_SeqReadsCount}->{case} < 0
			 || $V_Href->{dVG_SeqReadsCount}->{ctrl} < 0
			 || $V_Href->{dVG_WdSize} <= 0
			 || $V_Href->{dG_x_LabStep} <= 0
			 || $V_Href->{dG_y_resol} < 0
			 || $V_Href->{dG_y_LabStep} < 0
			 || $V_Href->{dVG_gene_label_idx} !~ /^[tglnp]$/i
			 # sub-region to draw
			 || (    defined($V_Href->{dVG_pos_start}) # input from user
			 	  && $V_Href->{dVG_pos_start} <= 0 # wrong input
			 	)
			 || (    defined($V_Href->{dVG_pos_end}) # input from user
			 	  && $V_Href->{dVG_pos_end} <= 0 # wrong input
			 	)
			 || (    defined($V_Href->{dVG_pos_start}) # input from user
			 	  && defined($V_Href->{dVG_pos_end}) # input from user
			 	  && ($V_Href->{dVG_pos_end} < $V_Href->{dVG_pos_start}) # wrong region
			 	)
			);
}

#--- depth stat of viral genome ---
sub draw_viral_genome{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	#---------------------#
	# validate virus type #
	#---------------------#
	&validate_virus_type;

	#------------#
	# load depth #
	#------------#
	my %WinDepth;
	&load_virus_depth;

	#-----------------------#
	# load viral breakpoint #
	# virus Integ, virusSV  #
	#-----------------------#
	&load_virus_breakpoint if( defined $V_Href->{dVG_viral_BkList} );

	#--------------------#
	# load mutation file #
	# snp, indel         #
	#--------------------#
	&load_viral_mutation if( defined $V_Href->{dVG_viral_MutFile} );

	#------------------------------------#
	# obtain span-width of each mutation #
	# based on text size of mut_labels   #
	#------------------------------------#
	&obtain_mut_span_width;

	#---------------------------#
	# arrange mut icon location #
	#---------------------------#
	&arrange_mut_icon_location;

	#---------------------#
	# smooth window depth #
	#---------------------#
	&smooth_window_depth;

	#----------------#
	# load trans psl #
	# it is gene     #
	#----------------#
	&load_viral_trans_info unless( $V_Href->{dVG_show_no_gene} );

	#-------------#
	# draw figure #
	#-------------#
	# set attributes of x axis #
	&set_x_axis;
	# set attributes of y axis #
	&set_y_axis;
	# draw y axis #
	# create SVG object
	&draw_y_axis;
	# draw mutations #
	&draw_mutations;
	# draw depth spectrum #
	&draw_depth;
	# draw gene in this region #
	&draw_genes unless( $V_Href->{dVG_show_no_gene} );
	# # draw x axis #
	&draw_x_axis;
	# show basic information and contact
	&show_figure_information;
	# output SVG #
	&generate_figure_output;
}

#--- validate virus type based on virus database ---
## assign related virus trans psl
sub validate_virus_type{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $virusRefSeg_Href = $V_Href->{virusRefSeg};

	# load 'virus.info'
	#FuseSV::Virus_Type::VirusInitAlign->load_virus_basic_information;
	FuseSV::Virus_Type::PrepareExtendVirusRef->load_virusBasicInfo;

	# check input virus type
	unless( exists $virusRefSeg_Href->{ $V_Href->{dVG_viral_subtype} } ){
		warn_and_exit `date`."<ERRO>\tCannot find input virus type ($V_Href->{dVG_viral_subtype}) in 'virus.info' file.\n"
									."\t$V_Href->{virus_info_doc}\n";
	}

	# assign related info
	$V_Href->{dVG_viral_origLen} = $virusRefSeg_Href->{ $V_Href->{dVG_viral_subtype} }->length;
	$V_Href->{dVG_viral_orig_tpsl} = FuseSV::GetPath->GetPath( filekey => 'single_viral_orig_tpsl', virus_subtype => $V_Href->{dVG_viral_subtype} );
	unless( -e $V_Href->{dVG_viral_orig_tpsl} ){
		warn_and_exit `date`."<ERRO>\tCannot find original trans.psl file of virus type ($V_Href->{dVG_viral_subtype}).\n"
									."\t$V_Href->{dVG_viral_orig_tpsl}\n";
	}
}

#--- load depth file ---
sub load_virus_depth{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# mark
	# $V_Href->{dG_CaseVsCtrl_mode} = ( $V_Href->{dVG_DepFile}->{case} && $V_Href->{dVG_DepFile}->{ctrl} );
	if( $V_Href->{dVG_SeqReadsCount}->{case} && $V_Href->{dVG_SeqReadsCount}->{ctrl} ){
		$V_Href->{dVG_case2norm_seqRatio} = $V_Href->{dVG_SeqReadsCount}->{ctrl} / $V_Href->{dVG_SeqReadsCount}->{case};
		# inform
		stout_and_sterr `date`."[INFO]:\tDepth will normlized with ratio $V_Href->{dVG_case2norm_seqRatio}\n";
	}

	# load vPos depth
	for my $read_mode (qw/ 1 2 /){
		for my $tissue ( qw/ case ctrl / ){
			my $DepthFile = $V_Href->{dVG_DepFile}->{$tissue};
			next if(!-e $DepthFile);
			# read file
			my $last_window_NO = 0;
			my $last_window_last_pos = 0;
			open (DEPTH, Try_GZ_Read($DepthFile)) || die "Cannot read DepthFile: $!\n";
			while(<DEPTH>){
				my @info = split;
				my $vRef = shift @info;
				my $vPos = shift @info;
				# check vRef
				next if( $vRef ne ($V_Href->{dVG_vRef_to_draw} ||= $vRef) );
				# draw sub-region, if user input
				next if( defined($V_Href->{dVG_pos_start}) && $vPos < $V_Href->{dVG_pos_start} );
				next if( defined($V_Href->{dVG_pos_end})   && $vPos > $V_Href->{dVG_pos_end} );
				# update max/min vPos
				if( $read_mode == 1 ){
					$V_Href->{dVG_x_maxPos} = max( $V_Href->{dVG_x_maxPos}, $vPos );
					$V_Href->{dVG_x_minPos} = min( $V_Href->{dVG_x_minPos}, $vPos );
					next;
				}
				# record depth, read-mode is 2
				my $window_NO = $last_window_NO; # initialize as last window
				if( $vPos > $last_window_last_pos ){
					my @winInfo = get_pos_windowNO(
													refseg => $V_Href->{dVG_vRef_to_draw},
													pos => $vPos,
													regResol_Href => $V_Href->{dG_regResol},
													stpos_1stwin => $V_Href->{dVG_x_minPos},
													give_winLastPos => 1
												);
					my $sameWinLastPos = $winInfo[0];
					$window_NO = $winInfo[1];
					# deal with last window
					if( $last_window_NO != $window_NO ){
						if( exists $V_Href->{dG_winDepth}->{$tissue}->{$window_NO} ){
							warn_and_exit "<ERRO>\ttissue $tissue window NO.$window_NO ($vPos) shows again.\n";
						}
						deal_window_depth_info(
												winDepth_Href => $V_Href->{dG_winDepth},
												window_NO => $last_window_NO,
												tissue => $tissue,
												role => 'virus',
												Eng_3SD => $V_Href->{dG_DepEng_3SD},
												TmmRatio => $V_Href->{dG_DepTrimRatio},
												CaseVsCtrl_mode => $V_Href->{dG_CaseVsCtrl_mode},
												v_case2norm_seqRatio => $V_Href->{dVG_case2norm_seqRatio}
											);
					}
					# update
					$last_window_NO = $window_NO;
					$last_window_last_pos = $sameWinLastPos;
				}
				# store depth info
				$V_Href->{dG_winDepth}->{$tissue}->{$window_NO}->{pos2depth}->{$vPos} = sum(@info);
			}
			close DEPTH;
			# very last one window
			if( $read_mode == 2 ){
				deal_window_depth_info(
										winDepth_Href => $V_Href->{dG_winDepth},
										window_NO => $last_window_NO,
										tissue => $tissue,
										role => 'virus',
										Eng_3SD => $V_Href->{dG_DepEng_3SD},
										TmmRatio => $V_Href->{dG_DepTrimRatio},
										CaseVsCtrl_mode => $V_Href->{dG_CaseVsCtrl_mode},
										v_case2norm_seqRatio => $V_Href->{dVG_case2norm_seqRatio}
									);
				# inform
				stout_and_sterr `date`."[INFO]:\tLoad depth of virus genome from file:\n"
											 ."\t$DepthFile\n";
			}
		}
		# create draw resolution hash
		if( $read_mode == 1 ){
			$V_Href->{dG_regResol}->{$V_Href->{dVG_vRef_to_draw}}->{$V_Href->{dVG_x_minPos}} = { st_pos => $V_Href->{dVG_x_minPos}, ed_pos=>$V_Href->{dVG_x_maxPos}, resol=>$V_Href->{dVG_WdSize} };
		}
	}
}

#--- read viral breakpoint and related labels ---
#--- default label is 'integ'.
sub load_virus_breakpoint{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $ViralBkpHref = $V_Href->{dG_bkWindow};
	my $ViralMutHref = $V_Href->{dVG_ViralMut};

	# load BkList
	my $mut_count = scalar( keys %$ViralMutHref );
	open(BKL, Try_GZ_Read($V_Href->{dVG_viral_BkList})) || die "fail read virus BkList: $!\n";
	while(<BKL>){
		next if(/^\#/);
		my ($vRef, $vBkPos, $label) = (split)[0,1,2];
		# check available
		if(    $vRef ne $V_Href->{dVG_vRef_to_draw} # wrong vRef
			|| $vBkPos < $V_Href->{dVG_x_minPos} # out of range
			|| $vBkPos > $V_Href->{dVG_x_maxPos} # out of range
		){
			next;
		}
		# count in
		$mut_count ++;
		# get window NO
		my $window_NO = get_pos_windowNO( refseg=>$V_Href->{dVG_vRef_to_draw}, pos=>$vBkPos, regResol_Href=>$V_Href->{dG_regResol}, stpos_1stwin=>$V_Href->{dVG_x_minPos} );
		$ViralBkpHref->{$window_NO} = $vBkPos;
		# record
		$ViralMutHref->{$mut_count}->{pos} = $vBkPos;
		# mutation line wNO
		$ViralMutHref->{$mut_count}->{wNO} = $window_NO;
		# mutation label wNO
		$ViralMutHref->{$mut_count}->{show_wNO} = $window_NO;
		# mutation
		$ViralMutHref->{$mut_count}->{mut_type} = 'integ_sv';
		$ViralMutHref->{$mut_count}->{ref_allele} = '';
		$ViralMutHref->{$mut_count}->{mut_allele} = '';
		$ViralMutHref->{$mut_count}->{mut_label} = "$vBkPos:".($label || 'integ');
	}
	close BKL;

	# inform
	stout_and_sterr `date`."[INFO]:\tLoad breakpoint(s) of virus genome OK.\n"
								 ."\t$V_Href->{dVG_viral_BkList}\n";
}

#--- load viral mutations ---
#--- types: snp, indel
sub load_viral_mutation{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $ViralMutHref = $V_Href->{dVG_ViralMut};

	# load mutations
	my $mut_count = scalar( keys %$ViralMutHref );
	open (MUT, Try_GZ_Read($V_Href->{dVG_viral_MutFile})) || die "fail read viral mutation file: $!\n";
	# theme
	(my $theme_line = lc(<MUT>)) =~ s/^#//;
	my @theme_tag = split /\s+/, $theme_line;
	while(<MUT>){
		next if(/^\#/);
		my @info = split;
		my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
		# check viral variant ref
		if( $rgOb{virusvariant} ne $V_Href->{dVG_vRef_to_draw} ){
			warn_and_exit `date`."[INFO]:\tEncounter different virus ref-seg name from mutation file.\n"
									   ."\t$rgOb{virusvariant} ne $V_Href->{dVG_vRef_to_draw}\n";
		}
		# skip minor mutation
		if(    $rgOb{majoralt_type} eq '-' # skip minor mutation
			|| $rgOb{originalpos} < $V_Href->{dVG_x_minPos} # out of range
			|| $rgOb{originalpos} > $V_Href->{dVG_x_maxPos} # out of range
		){
			next;
		}
		# count in
		$mut_count ++;
		# start record
		$ViralMutHref->{$mut_count}->{pos} = $rgOb{originalpos};
		# mutation line wNO
		$ViralMutHref->{$mut_count}->{wNO} = get_pos_windowNO( refseg=>$V_Href->{dVG_vRef_to_draw}, pos=>$rgOb{originalpos}, regResol_Href=>$V_Href->{dG_regResol}, stpos_1stwin=>$V_Href->{dVG_x_minPos} );
		# mutation label wNO
		$ViralMutHref->{$mut_count}->{show_wNO} = $ViralMutHref->{$mut_count}->{wNO};
		# mutation
		my $ref_allele = uc( $rgOb{ref_allele} );
		my $mut_allele = uc( $rgOb{majoralt_allele} );
		$ViralMutHref->{$mut_count}->{mut_type} = $rgOb{majoralt_type};
		$ViralMutHref->{$mut_count}->{ref_allele} = $ref_allele;
		$ViralMutHref->{$mut_count}->{mut_allele} = $mut_allele;
		# prepare mutation label
		if( $rgOb{majoralt_type} =~ /SN[PV]/i ){
			$ViralMutHref->{$mut_count}->{mut_label} = "$rgOb{originalpos}\:$ref_allele\>$mut_allele";
		}
		else{
			$ViralMutHref->{$mut_count}->{mut_label} = "$rgOb{originalpos}\:$rgOb{majoralt_type}\_$mut_allele";
		}
	}
	close MUT;

	# inform
	stout_and_sterr `date`."[INFO]:\tLoad mutations of virus genome OK.\n"
								 ."\t$V_Href->{dVG_viral_MutFile}\n";
}

#--- obtain icon span-width ---
sub obtain_mut_span_width{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $ViralMutHref = $V_Href->{dVG_ViralMut};

	for my $mut_NO (keys %$ViralMutHref){
		# get mut_label text size
		my ($label_height, $label_width) = 
			get_size_of_text_to_show(
										font_family => $V_Href->{general_FontFam},
										font_size => $V_Href->{dVG_mutLabelFontSize},
										text => $ViralMutHref->{$mut_NO}->{mut_label},
										return_HeightWid => 1
									);
		# pretend to rotate 90, so the text height is the span-x-range
		$ViralMutHref->{$mut_NO}->{label_height} = $label_height;
		$ViralMutHref->{$mut_NO}->{mut_width} = max( $V_Href->{dVG_mutIconEdgeLen}, $label_height );
	}
}

#--- arrange mut icon location iteratively ---
sub arrange_mut_icon_location{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $ViralMutHref = $V_Href->{dVG_ViralMut};

	allocate_object_bilater(
							obj_Href => $ViralMutHref,
							objIdx_key => 'pos',
							newLoc_key => 'show_wNO',
							objWid_key => 'mut_width',
							obj_gap => ($V_Href->{dVG_mutIconGap} || (scalar(keys %$ViralMutHref) > 50 ? 6 : 8)),
							adj_ratio => 0.25,
							timeToRpt => 1000,
							maxTryTime => 3E5
						   );

	# inform
	stout_and_sterr `date`."[INFO]:\tArrange location of mutation icon OK.\n";
}

#--- smooth window depth with bilateral flanking window ---
sub smooth_window_depth{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	for my $tissue ('case', 'ctrl'){
		my $DepthFile = $V_Href->{dVG_DepFile}->{$tissue};
		next if(!-e $DepthFile);
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
			stout_and_sterr `date`."[INFO]:\tSmooth window depth of virus genome of tissue $tissue OK.\n";
		}
	}

	# inform
	stout_and_sterr `date`."[INFO]:\tMaximum depth to dispaly is $V_Href->{dG_maxDepthToDraw}X.\n";
}

#--- load viral trans psl (it is gene) ---
sub load_viral_trans_info{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# change array input to hash
	for my $key (qw/ dVG_gene_to_skip /){
		my %temp = map{($_,1)} @{$V_Href->{$key}};
		$V_Href->{$key} = \%temp;
	}

	# database gene
	load_transInfo_from_tpsl(
								tpsl_file => $V_Href->{dVG_viral_orig_tpsl},
								TransInfo_Href => $V_Href->{dG_TransInfo},
								skipGene_Href => $V_Href->{dVG_gene_to_skip},
								clip_min_Pos => $V_Href->{dVG_x_minPos},
								clip_max_Pos => $V_Href->{dVG_x_maxPos},
								cirl_oLen => $V_Href->{dVG_viral_origLen},
								labelIdx => $V_Href->{dVG_gene_label_idx},
								fusesv_vmode => 1
							);
	# inform
	stout_and_sterr `date`."[INFO]:\tLoad transcript (viral gene) PSL file OK.\n"
								 ."\t$V_Href->{dVG_viral_orig_tpsl}\n";

	# user defined region
	if( scalar( @{$V_Href->{dG_user_gene}} ) != 0 ){
		add_user_region_to_gene_info(
										user_gene_input_Aref => $V_Href->{dG_user_gene},
										TransInfo_Href => $V_Href->{dG_TransInfo},
										refseg => $V_Href->{dVG_vRef_to_draw},
										clip_min_Pos => $V_Href->{dVG_x_minPos},
										clip_max_Pos => $V_Href->{dVG_x_maxPos},
										cirl_oLen => $V_Href->{dVG_viral_origLen},
										input_optKey_forAlert => '-nreg'
									);
		# inform
		stout_and_sterr `date`."[INFO]:\tLoad user defined region OK.\n";
	}
}

#--- set attributes of x axis ---
sub set_x_axis{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$V_Href->{dG_axisZX} = $V_Href->{dG_extWidth};
	$V_Href->{dG_x_axisLen} = $V_Href->{dG_maxWinNOtoDraw};

	# inform
	stout_and_sterr `date`."[INFO]:\tSet X axis OK.\n";
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
	stout_and_sterr `date`."[INFO]:\tSet Y axis OK.\n";
}

#--- draw y axis ---
## create SVG object
sub draw_y_axis{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# depth spectrum color legend info near Y axis
	my $yLegd_Href = undef;
	if( defined $V_Href->{dVG_spectrum_icon_text} ){
		my @yLegd_info = ();
		for my $tissue ( 'case', 'ctrl' ){
			next unless( $V_Href->{dVG_DepFile}->{$tissue} );
			push @yLegd_info, {
								shape => 'p',
								fill_col => $V_Href->{dG_DepthSpecColor}->{$tissue},
								text => $V_Href->{dVG_spectrum_icon_text}, # "$tissue depth"
								font_family => $V_Href->{general_FontFam_Bold},
								font_size => $V_Href->{general_FontSize}
							  };
		}
		$yLegd_Href = {loc=>'top', info_Aref=>\@yLegd_info};
	}

	# show Y axis
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
	stout_and_sterr `date`."[INFO]:\tDraw Y axis OK.\n";
}

#--- draw mutations icon and mut_label ---
sub draw_mutations{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $ViralMutHref = $V_Href->{dVG_ViralMut};

	# no mutation to draw
	my $mut_count = scalar(keys %$ViralMutHref);
	return if($mut_count == 0);

	# basic para
	my $base_layer_st_y = $V_Href->{dG_axisZY} - $V_Href->{dG_y_axisLen} - 10;
	my $base_layer_height = 4;
	my $mid_layer_height = ($mut_count > 50) ? 12 : 8;
	my $top_layer_height = 8;
	my $text_line_gap = 3;
	my $line_width = 1;
	# draw each mutation
	for my $mut_NO (sort {$a<=>$b} keys %$ViralMutHref){
		my $pos = $ViralMutHref->{$mut_NO}->{pos};
		my $wNO = $ViralMutHref->{$mut_NO}->{wNO};
		my $show_wNO = $ViralMutHref->{$mut_NO}->{show_wNO};
		my $mut_type = lc($ViralMutHref->{$mut_NO}->{mut_type});
		# path points
		my ($x1, $y1) = ($V_Href->{dG_axisZX}+$wNO, $V_Href->{dG_axisZY});
		my ($x2, $y2) = ($V_Href->{dG_axisZX}+$wNO, $base_layer_st_y-$base_layer_height);
		my ($x3, $y3) = ($V_Href->{dG_axisZX}+$show_wNO, $y2-$mid_layer_height);
		my ($x4, $y4) = ($V_Href->{dG_axisZX}+$show_wNO, $y3-$top_layer_height);
		# draw mutation mark lines
		$V_Href->{general_SVG_obj}->path( d=>"M$x1,$y1,L$x2,$y2", stroke=>'black', 'stroke-width'=>$line_width, 'stroke-linecap'=>'round' );
		$V_Href->{general_SVG_obj}->path( d=>"M$x2,$y2,L$x3,$y3", stroke=>'black', 'stroke-width'=>$line_width, 'stroke-linecap'=>'round' );
		$V_Href->{general_SVG_obj}->path( d=>"M$x3,$y3,L$x4,$y4", stroke=>'black', 'stroke-width'=>$line_width, 'stroke-linecap'=>'round' );
		# icon info
		my $icon_type = $V_Href->{dVG_viral_mut_icon}->{$mut_type}->{icon};
		if(!$icon_type){
			warn_and_exit `date`."<ERRO>\tCannot find icon for mutation type: $mut_type.\n";
		}
		# icon location
		my $icon_mx = $x3;
		my $icon_my = ($y3 + $y4) / 2;
		# draw icon
		if( $icon_type eq 'circle'){
			$V_Href->{general_SVG_obj}->circle(
												cx => $icon_mx,
												cy => $icon_my,
												r => $V_Href->{dVG_mutIconEdgeLen} / 2,
												stroke => 'black',
												'stroke-width' => 0.5,
												fill => $V_Href->{dVG_viral_mut_icon}->{$mut_type}->{col}
											);
		}
		elsif($icon_type =~ /triangle/){
			draw_a_triangle(
								\$V_Href->{general_SVG_obj},
								x => $icon_mx,
								y => $icon_my,
								fill_color => $V_Href->{dVG_viral_mut_icon}->{$mut_type}->{col},
								bottom_side_len => $V_Href->{dVG_mutIconEdgeLen},
								boundary_width => 0.5,
								rotate_degree => ($icon_type =~ /^rev_/) ? 180 : 0
							);
		}
		elsif($icon_type eq 'parallelogram'){
			draw_a_parallelogram(
									\$V_Href->{general_SVG_obj},
									x => $icon_mx,
									y => $icon_my,
									fill_color => $V_Href->{dVG_viral_mut_icon}->{$mut_type}->{col},
									head_bottom_side_len => $V_Href->{dVG_mutIconEdgeLen},
									boundary_width => 0.5
								);
		}
		# show mut-label text
		my $text_height = $ViralMutHref->{$mut_NO}->{label_height};
		show_text_in_line(
							\$V_Href->{general_SVG_obj},
							text_x => $icon_mx,
							text_y => $y4-$text_line_gap+$text_height/2,
							text => $ViralMutHref->{$mut_NO}->{mut_label},
							font_family => $V_Href->{general_FontFam},
							font_size => $V_Href->{general_FontSize}-2,
							text_col => 'black',
							text_anchor => 'start',
							height_adjust => 0,
							height_limit => 0,
							width_limit => 0,
							rotate_center_x => $icon_mx,
							rotate_center_y => $y4-$text_line_gap,
							rotate_degree => -45
						);
	}

	# inform
	stout_and_sterr `date`."[INFO]:\tDraw viral mutations OK.\n";
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
	unless( $V_Href->{dVG_show_no_gene} ){
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
	stout_and_sterr `date`."[INFO]:\tDraw viral depth spectrum OK.\n";
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

	my @get_pos_windowNO_setting = (
									refseg => $V_Href->{dVG_vRef_to_draw},
									regResol_Href => $V_Href->{dG_regResol},
									stpos_1stwin => $V_Href->{dVG_x_minPos}
								   );

	allocate_genes_and_draw(
							SVG_obj_Oref => \$V_Href->{general_SVG_obj},
							TransInfo_Href => $V_Href->{dG_TransInfo},
							pos2wNO_para_Aref => \@get_pos_windowNO_setting,
							axisZX_Sref => \$V_Href->{dG_axisZX},
							axisZY_Sref => \$V_Href->{dG_axisZY},
							tnameFtfm => $V_Href->{general_FontFam},
							tnameFtsz => $V_Href->{general_FontSize}
						   );

	# inform
	stout_and_sterr `date`."[INFO]:\tShow genes OK.\n";
}

#--- draw x axis ---
sub draw_x_axis{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# x-axis Label Text
	my $xLabelText = $V_Href->{dVG_viral_subtype};
	$xLabelText .= "\n".$V_Href->{dVG_sampleIDtoShow} if( defined $V_Href->{dVG_sampleIDtoShow} );
	# show x-axis
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
				# xRegResolColor => 'blue', # no-use
				xResol_Href => $V_Href->{dG_regResol},
				regSegName => $V_Href->{dVG_vRef_to_draw},
				xRegMinPos => $V_Href->{dVG_x_minPos},
				xRegMaxPos => $V_Href->{dVG_x_maxPos}
			);

	# inform
	stout_and_sterr `date`."[INFO]:\tDraw X axis OK.\n";
}

#--- generate figure file ---
sub generate_figure_output{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# file type (based on postfix)
	my ($figure_postfix) = ($V_Href->{dG_figure} =~ /\.([sp][vnd][gf])$/i);
	# temp svg
	my $out_svg = "$V_Href->{dG_figure}.svg";
	open (SVGO,">$out_svg")||die"fail write temp SVG figure: $!\n";
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
	stout_and_sterr `date`."[INFO]:\tGenerate final figure OK.\n";
}

#--- show basic figure information ---
# use a Golden rectangle
# figure name, version, author and institute
sub show_figure_information{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $fNameAref = $parm{fNameAref} || [ 'Viral Genome Map', 'after Reconstruction' ];
	my $fVersion = $parm{fVersion} || $VERSION;
	my $fAuthor = $parm{fAuthor} || $AUTHOR;
	my $fAffili = $parm{fAffili} || 'CityU of Hong Kong';

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
					# ["Viral Genome Map", 'Times Bold', 1],
					# ["after Reconstruction", 'Times Bold', 1],
					["[------- v$fVersion -------]", 'Times', 2, {font_family=>{'Times Italic'=>[10 .. 9+length($fVersion)+1]}}],
					["author: $fAuthor", 'Times', 2],
					[$fAffili, 'Times', 2]
				   );
	unshift @contents, [ $_, 'Times Bold', 1 ] for reverse @{$fNameAref};
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

#--- 
1; ## tell the perl script the successful access of this module.

