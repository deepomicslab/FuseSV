package FuseSV::Virus_Integ::LocalHaplotypeCircos::DrawLocalMapCircos;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/file_exist/;
use BioFuse::Visual::SVG_Util::RadSys qw/ $PI $deg2rad /;
use BioFuse::Visual::SVG_Util::Color qw/%COLOR_DB/;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Features qw/ load_features_source_type /;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal qw/ load_sample_seg_info /;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype qw/ load_local_haplotype /;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::SVG_draw qw / svg_wrok_pipe /;
# use Storable qw/dclone/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
				draw_local_map_circos
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
								 OTHER   => [qw()]);

$MODULE_NAME = 'DrawLocalMapCircos';
#----- version --------
$VERSION = "1.13";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
							return_HELP_INFO
							Load_moduleVar_to_pubVarPool
							Get_Cmd_Options
							para_alert
							draw_local_map_circos
						/;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
	 Usage:   perl $V_Href->{MainName} LMcircos <[Options]>

	 Options:

		 # Inputs and Outputs #
			-psg    [s]  List file of Segments information of samples. <required>
			-hap    [s]  The local haplotype of samples. <required>
			-fts    [s]  List file of Feature on host genome breakpoints. <required>
			-fet    [s]  Features of breakpoints of samples. <required>
			-rra    [s]  List file of Samples' RNA-seq reads amount for normalization. <required>
			-rpd    [s]  RNA depth of sorted positions for each samples. <required>
			-tpsl   [s]  Transcript PSL file provied by SOAPfuse. <required>
			-out    [s]  The output SVG file. <required>

		 # Options #
			-vnm    [s]  The virus name. ['your_virus']
			-vline       Use this to state the virus genome is linear. [false]
			              Default is false, meaning circular virus.
			-vp     [i]  The percentage that the virus segment will occupied in the whole circle. [25]
			-vmlp   [f]  The perctenage of the circle that the middle line of the virus segment
			             rotates from the twelve clock line. [50]
			              Note: 50 means let the virus at the bottom of the circle.
			-msp    [i]  The minimum perctenage that a host segment will occpy the whole circle. [7]
			-ct     [s]  The coordinate of centre point of the circle in the out SVG figure. ['300,300']
			              Note: In SVG, up-left corner is axis base point.
			                    X axis is heading right; Y axis is heading down.
			-hrds   [i]  The radius of the host segments on the circle. [200]
			-vrds   [i]  The radius of the virus segments on the circle. [150]
			-srds   [s]  To specify specific radius for some host region(s). [NULL]
			              Could be applied multiple times.
			              Instance: -srds chr8:128255000-128755000:100 -srds chr3:94795000-94850000:120
			-sts    [s]  Sorting String for segments to show in clock-wise. [1-22,X,Y]
			              Note: To optimize apparence, you may need to adjust the arrangement of segments,
			                    input someing like, '12,1,9,8,X,10'.
			-dgch   [s]  Signle PNC character to represent degenerated path arc.
			              Could be applied multiple times, e.g., -dgch '#' -dgch '='
			              Default is: @{$V_Href->{PNC_deg_char}}
			-arcreg      Just show the region on arc (path-seg), not the path number. [false]
			-hhf         Show and count in the features of host-to-host stucture variation.
			              In defaults, only the breakpoints between host and virus will be counted in.
			-osuft       Only to show the feature(s) used in the virus-integration haplotype(s). [disabled]
			-gene   [s]  genes only allowed to show on segments.
			              In defaults, all genes are allowed. Input as, 'TP53,VHL,ROS1,MYC'.
			-cnlg   [i]  Log Scale to display DNA copy number of segments. [0, usual linear]
			-eplg   [i]  Log Scale to display RNA expression of segments. [0, usual linear]

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
			#-------------------------#
			#--- inputs and output ---#
			#-------------------------#
			[ Feature_source => undef ],
			[ Sample_SegInfo_file => undef ],
			[ Sample_LocalHaplotype_file => undef ],
			[ Sample_Break_Feature => undef ],
			[ Sample_RNA_Reads_Amount_file => undef ],
			[ Sample_RNA_Pos_Depth_file => undef ],
			[ Trans_PSL => undef ],
			[ output_SVG_file => undef ],

			#-------------------------#
			#--- global containers ---#
			#-------------------------#
			## HostSite_Feature (theme_id)
			## HostSite_Feature -> 'type' -> {type} -> {}
			##                 {type} -> {icon=>'', colour=>''}
			## HostSite_Feature -> 'feature' -> {feature_id} -> { type=>'', show_text=>'', short_name=>'', sn_index=>'',
			##                            show_bool=>0/1, counts=>{'host'=>xx, 'virus'=>xx}}
			[ Feature_Basic_Href => { type => {
													cause =>  {icon=>'square',   colour=>'rgb(255,70,70)'},
													effect  =>  {icon=>'triangle', colour=>'yellow'}
												 },
									feature => {}
									} ],
			## Segment_Href -> $seg_NO = $seg_OB
			## NOTE: seg_NO = 0, virus; seg_NO > 0, host;
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB module
			[ Segment_Href => {} ],
			## Sample_Href -> $sam_id = $sam_OB
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeCircos::Sample_OB module
			[ Sample_Href => {} ],
			## Public_Container_Records for re-draw
			[ Public_Container_Records => {} ],

			#-------------#
			#--- basic ---#
			#-------------#

			#=== apparence settings ===#
			## virus name
			[ virus_name => 'your_virus' ],
			## virus part in the circle, x%
			[ virus_perct => 25 ],
			## The perctenage of the circle that the middle line of the virus
			## segment rotates from the twelve clock line. [50]
			[ virus_middleline_loc_perct => 50 ],
			## 1, this is a linear virus; 0, this is a circular virus
			[ viral_linear => 0 ],
			## minimum perctage each sub-seg occupies
			[ min_seg_perct => 3 ],
			## coordinates of circle center
			[ center => '300,300' ],
			## host section radius, this is for the most inner edge
			## unless sections given by specific radius
			[ host_radius => 200 ],
			## virus section radius
			[ virus_radius => 150 ],
			## radius specification
			[ radius_specific => [] ],
			## segment sorting
			## human being
			[ segment_sort_str => '1-22,X,Y' ],
			## hided parameter, only overlapped over such ratio.
			[ min_spec_radius_overlap_ratio => 0.3 ],
			## show arc region (path-seg), not the path number.
			## for adjustment on resolutions
			[ show_arc_reg => 0 ],
			## includes features of human-human connection points
			[ show_hhf => 0 ],
			## for dna copy-number scale
			[ dna_cn_log => 0 ], # use linear
			## for exp depth scale
			[ exp_log => 0 ], # use linear
			## required genes list
			[ allow_gene_str => 0 ],
			[ allow_gene_Href => {} ],
			[ genename_anchor_gTail => 0 ],
			[ genename_anchor_gHead => 0 ],
			[ genename_anchor_Href => {} ],
			## only show applied feature
			[ only_show_used_feature => 0 ],

			#----------------#
			#--- advanced ---#
			#----------------#

			[ circos_OB => '' ],
			#=== circos stucture ===#
			[ cx => 100 ],
			[ cy => 100 ],
			[ backend_width => 800 ],
			[ backend_height => 600 ],
			## gap size from svg-board-bound to figure elements
			[ figure_gap_from_svg_bound => 30 ],

			#=== segments apparence ===#
			## segment width, ratio to multiply with segment_radius
			[ segment_width_ratio => 0.06 ],
			[ segment_boundary_width => 0.5 ],
			[ segment_name_font_size => 8 ],
			[ segment_bandpos_font_size => 4.5 ],
			[ segment_bandpos_axisline_length => 2 ],
			[ segment_name_band_width => 8 ],
			## gap's rad between segments
			[ spacerad => 1.2 * $deg2rad ],
			## colour
			[ virus_segcol => 'antiquewhite' ],
			[ host_segcol => 'aliceblue' ],
			## segment axis
			[ axis_top_rad => [0.56 * $PI, 1.44 * $PI] ],

			#=== haplotype path ===#
			[ PNC_font_size => 8 ],
			[ path_extend => 0.5 ], # the extended part for path to go different Segment layer
			[ path_width => 0.5 ],
			[ path_line_opacity => 0.5 ],
			[ line_dash => '1,2' ],
			## about PNC
			[ PNC_radius => 4 ],
			[ PNC_deg_char => [ '+', '*', '#', '=', '%', '$', '/', '~', '&', 'T', 'x', 's',
								# '[', ']', '!', ':', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'
								] ],

			#=== text in figure ===#
			[ Font_family => "Times New Roman" ], # ArialNarrow, Arial
			[ Font_family_SegName => 'Arial Bold' ],

			#=== feature icon on sites ===#
			[ icon_edge_len => 8 ],
			[ icon_radial_gap => 2 ],
			[ icon_tangent_mingap => 2 ],
			[ icon_font_size => 6 ],
			[ icon_line_width => 0.5 ],
			[ icon_line_len => [2,2,2] ], # from inner to outer, icon line length
			[ icon_loc_adjust_rad_step => 0.002 ],
			[ icon_loc_max_arrange_time => 100000 ],

			#=== genes on segments ===#
			[ gene_col_plus => 'royalblue' ],
			[ gene_col_minus => 'hotpink' ],
			[ gene_col_region => 'orange' ],
			[ gene_name_fontfam => 'Arial Italic' ],
			[ gene_intron_width => 0.5 ],
			[ gene_exon_width => 4 ],
			[ gene_arrow_bottom_len => 3 ],
			[ min_exon_length => 0.5 ],
			[ gene_CDS_width => 5 ],
			[ gene_name_fontsize => 6 ],
			[ gene_name_height_limit => 6 ],
			[ gene_radial_gap => 2.5 ],
			[ gene_tangent_gap => 2 ],
			[ gene_name_gap => 2 ],

			#=== DNA copy number on segments ===#
			[ CN_bg_color => 'lightblue' ],
			[ CN_color_rgb => '255,0,0' ], # not the color-symbol, for gradiant changing of RGB
			[ CN_min => 1 ], # acturely, 1 is the minimum; change to 1E5, when want to get real minimum
			[ CN_max => 1 ],
			[ dna_CN_band_width => 3 ],

			[ CN_EXP_sample_note_rad_ratio => 0.7 ],

			#=== expression on segments ===#
			[ EXP_bg_color => 'lightblue' ],
			[ EXP_color_rgb => '0,129,0' ], # not the color-symbol ('green'), for gradiant changing of RGB
			[ EXP_min => 1E10 ],
			[ EXP_max => 0 ],
			[ rna_EXP_band_width => 3 ],
			[ EXP_win_arc_len => 1 ], # to automatically determine window_rad for each segment, fixed by sample at the innner layer
			[ RRA_max => 0 ], # record the rna reads number among all samples
			[ EXP_depth_buffer => 1500 ],
			[ EXP_Depth_Sample_Index_Href => {} ],

			#=== feature pan ===#
			[ ft_cx => 720 ],
			[ ft_cy => 220 ],
			[ ft_pan_gap_from_ceiling => 30 ],
			[ ft_pan_icon_width => 20 ],
			[ ft_pan_gap_from_segments => 30 ],
			[ ft_centre_v_cicle_radius => 8 ],
			[ ft_pan_radius_inner => 50 ],
			[ ft_pan_radius_outer => 60 ],
			[ ft_pan_circos_link_width => 3 ],

			#=== figure legend ===#
			[ legend_text_fontfam => 'Times' ],
			[ legend_SamId_fontsize => 7 ],
			[ legend_Samline_width => 2 ],
			[ legend_Samline_len => 10 ],
			[ legend_MainTheme_fontfam => 'Times Bold' ],
			[ legend_SubTheme_fontfam => 'Times Italic' ],
			[ legend_MainTheme_ft_Aref => [8, 30] ],
			[ legend_SubTheme_ft_Aref => [7, 29] ],
			[ legend_text_ft_Aref => [6, 28] ],
			[ legend_text_ft_size => 6 ],
			[ legend_icon_min_width => 8 ],
			[ legend_Theme_underline_width => 2 ],
			[ legend_text_underline_width => 1 ],
			[ legend_Theme_underline_gap => 4 ],
			[ legend_themeline_row_gap => 4 ],
			[ legend_theme_row_gap => 8 ],
			[ legend_icon_text_gap => 5 ],
			[ legend_column_gap => 15 ],
			[ legend_near_column_min_ratio => 0.5 ],
			#--- gradiant color for CN and expression
			[ legend_CnExp_grad_band_width => 30 ],
			[ legend_grad_band_axis_width => 0.5 ],
			[ legend_grad_band_axis_height => 1 ],
			[ legend_grad_band_axis_gap => 2 ],
			[ legend_grad_band_axis_step_NO => 4 ],
			[ legend_grad_band_axis_fontsize => 5 ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['Feature_source'],
									['Sample_SegInfo_file'],
									['Sample_LocalHaplotype_file'],
									['Sample_Break_Feature'],
									['Sample_RNA_Reads_Amount_file'],
									['Sample_RNA_Pos_Depth_file'],
									['Trans_PSL'],
									['output_SVG_file']       ] ]
		);

	#---- self adjustment, reload variants ---#
	#=== haplotype path ===#
	## gap between haplotype path layers
	$V_Href->{seg_path_space} = $V_Href->{PNC_radius} * 2;
	## minimum arc length
	$V_Href->{min_arc_len} = $V_Href->{PNC_radius} * 2.732 + 2;
	#=== feature icon on sites ===#
	$V_Href->{ft_line_tiny_head} = $V_Href->{icon_space};
	$V_Href->{ft_icon_remap_bandwidth} = $V_Href->{icon_edge_len};
	#=== feature pan ===#
	# $V_Href->{ft_cx} = $V_Href->{cx};
	# $V_Href->{ft_cy} = $V_Href->{ft_pan_gap_from_ceiling} + $V_Href->{ft_pan_radius_outer};
	$V_Href->{ncl_radius} = $V_Href->{ft_centre_v_cicle_radius} / 1.5;
	$V_Href->{ft_pan_icon_font_size} = $V_Href->{ft_pan_icon_width};

	# delete at last, just for normally test
	$V_Href->{Virus_seg_Href} = $V_Href->{Segment_Href}->{0};
}

#--- get options from command line ---
sub Get_Cmd_Options{

	# variants

	# get options
	GetOptions(
		# input files
		"-fts=s"  => \$V_Href->{Feature_source},
		"-psg=s"  => \$V_Href->{Sample_SegInfo_file},
		"-hap=s"  => \$V_Href->{Sample_LocalHaplotype_file},
		"-fet=s"  => \$V_Href->{Sample_Break_Feature},
		"-rra=s"  => \$V_Href->{Sample_RNA_Reads_Amount_file},
		"-rpd=s"  => \$V_Href->{Sample_RNA_Pos_Depth_file},
		"-tpsl=s" => \$V_Href->{Trans_PSL},
		# output file
		"-out=s"  => \$V_Href->{output_SVG_file},
		# apparence
		"-vnm=s"  => \$V_Href->{virus_name},
		"-vline"  => \$V_Href->{viral_linear},
		"-vp=i"   => \$V_Href->{virus_perct},
		"-vmlp=f" => \$V_Href->{virus_middleline_loc_perct},
		"-msp=i"  => \$V_Href->{min_seg_perct},
		"-ct=s"   => \$V_Href->{center},
		"-hrds=i" => \$V_Href->{host_radius},
		"-vrds=i" => \$V_Href->{virus_radius},
		"-srds=s" => \@{$V_Href->{radius_specific}},
		"-sts=s"  => \$V_Href->{segment_sort_str},
		"-dgch=s" => \@{$V_Href->{PNC_deg_char}},
		"-arcreg" => \$V_Href->{show_arc_reg},
		"-hhf"    => \$V_Href->{show_hhf},
		"-osuft"  => \$V_Href->{only_show_used_feature},
		"-cnlg=i" => \$V_Href->{dna_cn_log},
		"-eplg=i" => \$V_Href->{exp_log},
		# gene
		"-gene=s" => \$V_Href->{allow_gene_str},
		"-atgn=s" => \$V_Href->{genename_anchor_gTail},
		"-ahgn=s" => \$V_Href->{genename_anchor_gHead},
		# help
		"-h|help" => \$V_Href->{HELP},
		# for debug
		"-debug"  => \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants

	return  (   $V_Href->{HELP}
			 || !file_exist( filePath=>$V_Href->{Sample_SegInfo_file} )
			 || !file_exist( filePath=>$V_Href->{Sample_LocalHaplotype_file} )
			 || (    file_exist( filePath=>$V_Href->{Sample_Break_Feature} )
					 && !file_exist( filePath=>$V_Href->{Feature_source} )
					)
			 || (    file_exist( filePath=>$V_Href->{Sample_RNA_Pos_Depth_file} )
					 && !file_exist( filePath=>$V_Href->{Sample_RNA_Reads_Amount_file} )
					)
			 || !file_exist( filePath=>$V_Href->{Trans_PSL} )      
			);
}

#--- draw local map circos ---
sub draw_local_map_circos{

	# variants

	# features
	load_features_source_type;

	# sample information
	## both host segments and viral segments
	load_sample_seg_info;

	# load the haplotype of each sample
	load_local_haplotype;

	# pretend to draw
	svg_wrok_pipe( draw_bool => 1 );
}

#--- 
1; ## tell the perl script the successful access of this module.
