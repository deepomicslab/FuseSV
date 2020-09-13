package FuseSV::Virus_Integ::LocalHaplotypeCircos::SVG_draw;

use strict;
use warnings;
use SVG;
use Math::Trig;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use BioFuse::Util::Sys qw/file_exist/;
use BioFuse::Visual::SVG_Util::RadSys qw/$PI $deg2rad get_coordinate_on_circle normalize_radian/;
use BioFuse::Visual::SVG_Util::RadSysEle qw/draw_circle_seg/;
use BioFuse::Visual::SVG_Util::Color qw/%COLOR_DB/;
use BioFuse::Visual::SVG_Util::RectSysEle qw/draw_a_parallelogram draw_a_triangle draw_a_ellipse/;
use BioFuse::Visual::SVG_Util::Font qw /show_text_in_line show_text_on_arc get_size_of_text_to_show/;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Features;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
				svg_wrok_pipe
				draw_rna_exp_on_segments
				draw_copy_number_on_segments
				draw_genes_on_segments
				draw_pos_features
				draw_feature_stat_pan
				draw_segment_bands
				draw_haplotype_path
				makeup_PNC_ArcDot
				makeup_segment_sector
				show_figure_legends
				cn_exp_convert
				show_segment_name
				start_svg
				output_svg_file
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
								 OTHER   => [qw()]);

$MODULE_NAME = 'SVG_draw';
#----- version --------
$VERSION = "0.05";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
							svg_wrok_pipe
							draw_rna_exp_on_segments
							draw_copy_number_on_segments
							draw_genes_on_segments
							draw_pos_features
							draw_feature_stat_pan
							draw_segment_bands
							draw_haplotype_path
							makeup_PNC_ArcDot
							makeup_segment_sector
							show_figure_legends
							cn_exp_convert
							show_segment_name
							start_svg
							output_svg_file
						/;

#--- all works on SVG drawing
sub svg_wrok_pipe{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# create svg object
	&start_svg;

	# draw all segments, including host and virus
	&draw_segment_bands( draw_bool => $draw_bool );

	# draw haplotype path (update segment's outer_radius)
	&draw_haplotype_path( draw_bool => $draw_bool );

	# load features and show on radian segments
	&draw_pos_features( draw_bool => $draw_bool ) if file_exist( filePath=>$V_Href->{Sample_Break_Feature} );

	# show feature stat pan
	&draw_feature_stat_pan if( $draw_bool );

	# show genes
	&draw_genes_on_segments( draw_bool => $draw_bool );

	# show copy number of each segments of each sample
	&draw_copy_number_on_segments( draw_bool => $draw_bool );

	# show RNA expression
	&draw_rna_exp_on_segments( draw_bool => $draw_bool ) if file_exist( filePath=>$V_Href->{Sample_RNA_Pos_Depth_file} );

	# show segments' name
	&show_segment_name( draw_bool => $draw_bool );

	# show sample colour path
	# and, show degenerated PNCs
	# duplicated contigs
	# show features icon
	# show genes
	# show dna copy number, and rna expression
	&show_figure_legends( x => 60, y => 120 ) if( $draw_bool );

	&output_svg_file;
}

#--- draw segments expression of each sample
sub draw_rna_exp_on_segments{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# draw graded colored rna expression
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->load_sample_rna_reads_amount;

	# draw graded colored rna expression
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->draw_rna_exp_band( draw_bool => $draw_bool );
}

#--- show copy number of each segments of each sample
sub draw_copy_number_on_segments{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# draw graded colored CN band
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->draw_copy_number_band( draw_bool => $draw_bool );
}

#--- load trans.psl and draw genes on segments
sub draw_genes_on_segments{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# load genes to segment Objects from transcript psl file
	# select proper transcript to represent genes
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->load_gene;

	# draw genes on each segment
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->draw_gene( draw_bool => $draw_bool );

	# make up segment's feature icons sector
	&makeup_segment_sector( draw_bool => $draw_bool );
}

#--- load features and show on radian segments 
sub draw_pos_features{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# load features to sample Objects
	# and distribute these pos with features to segments
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Features->load_features_list;

	# draw feature icons
	FuseSV::Virus_Integ::LocalHaplotypeCircos::Features->draw_feature_icon_on_seg( draw_bool => $draw_bool );

	# make up segment's feature icons sector
	&makeup_segment_sector( draw_bool => $draw_bool );
}

#--- show feature stat pan
sub draw_feature_stat_pan{

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};

	# draw background host seg circle
	for my $tmp_Aref (  [$V_Href->{ft_pan_radius_outer}, $V_Href->{host_segcol}],
						[$V_Href->{ft_pan_radius_inner}, 'white']
					 ){
		draw_a_ellipse(
						svg_obj => $V_Href->{circos_OB},
						cx => $V_Href->{ft_cx},
						cy => $V_Href->{ft_cy},
						radius => $tmp_Aref->[0],
						boundary_width => 1,
						boundary_color => 'gray',
						boud_dasharray => $V_Href->{line_dash},
						fill_color => $tmp_Aref->[1]
					);
	}

	# allocate the feature squares
	my %features_to_show;
	my $feature_count = 0;
	for my $fet_id (keys %{$Feature_Basic_Href->{feature}}){
		my $fet_id_Href = $Feature_Basic_Href->{feature}->{$fet_id};
		next unless($fet_id_Href->{show_bool});
		my $fet_type = $fet_id_Href->{type};
		# if only show virus-partner
		if( !$V_Href->{show_hhf} && $fet_id_Href->{counts}->{virus}==0){
			next;
		}
		my $fet_shortname = $fet_id_Href->{short_name};
		$features_to_show{$fet_type}->{$fet_shortname}->{$_} += $fet_id_Href->{counts}->{$_} for keys %{$fet_id_Href->{counts}};
		$feature_count ++;
	}

	# linkage amount and feature icon
	if( $feature_count != 0 ){
		my $fet_pan_radius = ($V_Href->{ft_pan_radius_outer} + $V_Href->{ft_pan_radius_inner}) / 2;
		my $fet_rad_unit = 2 * $PI / $feature_count;
		my $index_to_show = 0;
		for my $fet_type (sort keys %features_to_show){
			for my $fet_shortname (sort keys %{$features_to_show{$fet_type}}){
				$index_to_show ++;
				my $fet_rad = $fet_rad_unit * ($index_to_show - 1);
				my ($ft_x,$ft_y) = get_coordinate_on_circle(cx=>$V_Href->{ft_cx}, cy=>$V_Href->{ft_cy}, rad=>$fet_rad, radius=>$fet_pan_radius);
				my $icon_shape = $Feature_Basic_Href->{type}->{$fet_type}->{icon};
				my $icon_fill_col = $Feature_Basic_Href->{type}->{$fet_type}->{colour};
				#----------------#
				# linkage amount #
				#----------------#
				# linkage amount circle for virus partner
				my $linkage_amount = $features_to_show{$fet_type}->{$fet_shortname}->{virus};
				next if $linkage_amount == 0;
				# radian path
				my $path_radius = $linkage_amount / 4 + ($fet_pan_radius ** 2) / (2 * $linkage_amount);
				my $add_rad = atan( $fet_pan_radius / $linkage_amount ) - atan( $linkage_amount / $fet_pan_radius );
				my $path = "M$ft_x,$ft_y,A$path_radius,$path_radius,0,0,1,$V_Href->{ft_cx},$V_Href->{ft_cy},A$path_radius,$path_radius,0,0,1,$ft_x,$ft_y";
				$V_Href->{circos_OB}->path(
											d=>$path,
											stroke=>'black',
											'stroke-width'=>0.5,
											'stroke-linecap'=>'round',
											'stroke-dasharray'=>' ',
											fill=>$icon_fill_col,
											opacity=>1
										 );
				# draw link number cicle (ncl)
				my ($ncl_x, $ncl_y) = ( ($ft_x+$V_Href->{ft_cx}) / 2, ($ft_y+$V_Href->{ft_cy}) / 2);
				my $ncl_radius = max( $V_Href->{ncl_radius}, log($linkage_amount)/log(3) * 3);
				# best font size
				my $ncl_font_size = min( int(0.8 * 4 * $ncl_radius / length($linkage_amount)) , int(0.8 * 8 * $ncl_radius / 3) );
				# circle and text
				draw_a_ellipse(
								svg_obj => $V_Href->{circos_OB},
								cx => $ncl_x,
								cy => $ncl_y,
								radius => $ncl_radius,
								boundary_width => 0.5,
								fill_color => $icon_fill_col,
								text_fill => $linkage_amount,
								text_font_size => $ncl_font_size,
								font_size_auto_adjust => 0
							);
				#--------------#
				# feature icon #
				#--------------#
				# draw this feature icon
				FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->draw_a_feature_icon(
																						x => $ft_x,
																						y => $ft_y,
																						icon_shape => $icon_shape,
																						icon_EdgeLen => $V_Href->{ft_pan_radius_outer} - $V_Href->{ft_pan_radius_inner},
																						fill_color => $icon_fill_col,
																						feature_text => $fet_shortname,
																						rotate_degree => 0
																					);
			}
		}
	}

	# at last, draw the virus circle and text
	draw_a_ellipse(
					svg_obj => $V_Href->{circos_OB},
					cx => $V_Href->{ft_cx},
					cy => $V_Href->{ft_cy},
					radius => $V_Href->{ft_centre_v_cicle_radius},
					boundary_width => 0.5,
					fill_color => $V_Href->{virus_segcol},
					text_fill => 'V',
					text_font_size => $V_Href->{ft_pan_icon_font_size},
					font_size_auto_adjust => 1
				);
}

#--- draw host and virus segments
sub draw_segment_bands{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};

	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		# draw radian seg
		my $seg_col = ($seg_NO>0) ? $V_Href->{host_segcol} : $V_Href->{virus_segcol};
		draw_circle_seg(
							svg_obj => $V_Href->{circos_OB},
							cx => $V_Href->{cx},
							cy => $V_Href->{cy},
							start_rad => $seg_OB->{start_rad},
							rad_size => $seg_OB->{rad_size}, 
							inner_radius => $seg_OB->{inner_radius},
							outer_radius => $seg_OB->{inner_radius} + $seg_OB->{segment_width},
							arrow_num => 4,
							seg_fill_color => $seg_col,
							seg_boundary_width => $V_Href->{segment_boundary_width},
							draw_bool => $draw_bool
						 );
		# show the start pos and end pos
		my $degree_limit = $seg_OB->{rad_size} / 2 / $deg2rad;
		my @PosInfo;
		push @PosInfo, [$seg_OB->{lf_pos}, $seg_OB->{start_rad}, 'start'];
		push @PosInfo, [$seg_OB->{rt_pos}, $seg_OB->{start_rad}+$seg_OB->{rad_size}, 'end'];
		my $last_text_span_rad = 0;
		for my $info_Aref ( @PosInfo ){
			my $position = $info_Aref->[0];
			normalize_radian(radian_Sref=>\$info_Aref->[1]);
			my $anchor_degree = $info_Aref->[1] / $deg2rad;
			my $text_anchor = $info_Aref->[2];
			my $text_to_center = ($info_Aref->[1] >= $V_Href->{axis_top_rad}->[0] && $info_Aref->[1] <= $V_Href->{axis_top_rad}->[1]) ? 'top':'toe';
			# get text size
			my ($text_height, $text_width) = get_size_of_text_to_show(
												font_family => $V_Href->{Font_family_SegName},
												font_size => $V_Href->{segment_bandpos_font_size},
												text => $position,
												return_HeightWid => 1
											);
			my $text_span_deg = $text_width / $seg_OB->{inner_radius} * 1.2 / $deg2rad;
			my $height_shift_ratio = ($text_span_deg > $degree_limit && $text_anchor eq 'end') ? 2.2 : 1;
			my $top_toe_shift = ($text_to_center eq 'top') ? 1 : 0;
			# the axis line
			my ($x_1, $y_1) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$info_Aref->[1], radius=>$seg_OB->{inner_radius} );
			my ($x_2, $y_2) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$info_Aref->[1], radius=>$seg_OB->{inner_radius}-$V_Href->{segment_bandpos_axisline_length}-($height_shift_ratio-1)*$text_height );
			$V_Href->{circos_OB}->line(
										x1 => $x_1,
										y1 => $y_1,
										x2 => $x_2,
										y2 => $y_2,
										stroke => 'black',
										'stroke-width' => 0.5,
										'stroke-linecap' => 'round',
										opacity => 1
										);
			# show position
			my $text_radius = $seg_OB->{inner_radius} - $V_Href->{segment_bandpos_axisline_length} - 0.2 * $text_height - $text_height * ($height_shift_ratio-$top_toe_shift);
			show_text_on_arc(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								radius => $text_radius,
								text => $position,
								font_family => $V_Href->{Font_family_SegName},
								font_size => $V_Href->{segment_bandpos_font_size},
								text_col => 'black',
								anchor_degree => $anchor_degree,
								text_anchor => $text_anchor,
								text_to_center => $text_to_center,
								height_adjust => 0,
								# height_limit=>$height_limit,
								# degree_limit=> $seg_OB->{rad_size} / $deg2rad,
								draw_bool=>$draw_bool
								);
		}
	}
}

#--- draw haplotype path
sub draw_haplotype_path{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};
	my $Segment_Href = $V_Href->{Segment_Href};
	my $extend = $V_Href->{path_extend};
	my $cx = $V_Href->{cx};
	my $cy = $V_Href->{cy};

	# draw all sub-regions' arc path on each segments
	# and store features of PNC(s) and arc_end_dots for next drawing (to be the top-svg-layer)
	for my $sam_id (sort keys %$Sample_Href){
		my $sam_OB = $Sample_Href->{$sam_id};
		my $PNC_count = 0;
		for my $hap_NO (sort {$a<=>$b} keys %{$sam_OB->{haplotype}}){
			my $color = $sam_OB->{pathcolour}->{$hap_NO};
			my $hap_Href = $sam_OB->{haplotype}->{$hap_NO};
			my %last_reg = (seg_NO => -1, st_rad => 0, ed_rad => 0, ori => 0, layer_NO => 0);
			for my $contig_NO (sort {$a<=>$b} keys %$hap_Href){
				my $contig_Aref = $hap_Href->{$contig_NO};
				my $region_Aref = $contig_Aref->[1];
				for my $reginfo_Aref (@$region_Aref){
					my ($seg_NO, $refseg, $st_pos, $ed_pos, $ori) = @$reginfo_Aref;
					my $reg_name = "$sam_id:$refseg:$st_pos-$ed_pos";
					my $seg_OB = $Segment_Href->{$seg_NO};
					# get the rad
					my $st_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $st_pos );
					my $ed_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $ed_pos );
					# new arc, PNC add by one
					$PNC_count ++;
					push @$reginfo_Aref, $PNC_count; # record this PNC_NO for this seg in this contig
					# layer to arrange this_seg (ts)
					my $ts_layer_NO = 0;
					# information
					my $this_reg_Href = {seg_NO=>$seg_NO, st_rad=>$st_rad, ed_rad=>$ed_rad, ori=>$ori, layer_NO=>$ts_layer_NO};
					## last seg comes from other segemnts, no matter host or virus
					if( $seg_NO != $last_reg{seg_NO} ){
						#----------------------#
						# arrange this seg arc #
						#----------------------#
						# as the last one is from other seg, so next to draw two radian lines and one bezier curve
						my $draw_lastarc_radian_bool = 1; # initial to draw radian line of last_arc
						my $draw_thisarc_radian_bool = 1; # initial to draw radian line of this_arc
						my $draw_bezier_bool = 1; # initial to draw bezier line linking last_arc and this_arc
						# try to find the previous same arc of this sample
						if( !$V_Href->{show_arc_reg} ){
							$ts_layer_NO = $sam_OB->find_previous_arc(
																		hap_NO => $hap_NO,
																		last_arc => \%last_reg,
																		this_arc => $this_reg_Href,
																		PNC_NO => $PNC_count,
																		last_radline => \$draw_lastarc_radian_bool,
																		this_radline => \$draw_thisarc_radian_bool,
																		bezier => \$draw_bezier_bool,
																		record => 1
																	);
						}
						# find the proper layer to draw
						if( $ts_layer_NO == 0 ){ # have not find any previous layer
							$ts_layer_NO = $seg_OB->which_layer_to_draw( 
																			st_rad => $st_rad,
																			ed_rad => $ed_rad,
																			min_layer => 0,
																			extend => $V_Href->{path_extend},
																			record => 1
																		);
							if($draw_bool){
								# draw the path arc
								# and PNC circle on it
								$seg_OB->draw_path_arc(
														layer_NO => $ts_layer_NO,
														st_rad => $st_rad,
														ed_rad => $ed_rad,
														ori => $ori, # 1 for forward, -1 for reversed
														sam_id => $sam_id,
														hap_NO => $hap_NO,
														arc_NO => $PNC_count,
														arc_name => $reg_name
													);
							}
						}
						## draw two radian lines and one bezier curve
						if( $last_reg{seg_NO} != -1){ # this is not the first segment
							if($draw_bool){
								# draw radian line for this segment's st_pos
								if( $draw_thisarc_radian_bool ){
									$seg_OB->draw_radian_line(
																rad => $st_rad,
																ed_layer_NO => $ts_layer_NO,
																color => $color
															);
								}
								# draw radian line for last segment's ed_pos
								my $last_seg_OB = $Segment_Href->{ $last_reg{seg_NO} };
								if( $draw_lastarc_radian_bool ){
									$last_seg_OB->draw_radian_line(
																	rad => $last_reg{ed_rad},
																	ed_layer_NO => $last_reg{layer_NO},
																	color => $color
																);
								}
								# draw bezier curve
								if( $draw_bezier_bool ){
									FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype->bezier_between_segments(
										seg1_OB => $last_seg_OB,
										seg1_rad => $last_reg{ed_rad},
										seg2_OB => $seg_OB,
										seg2_rad => $st_rad,
										color => $color
									);
								}
							}
						}
					}
					else{ # last seg is at the same seg with this seg, needs to deal with different cases
						# try to find the previous same arc of this sample
						if( !$V_Href->{show_arc_reg} ){
							$ts_layer_NO = $sam_OB->find_previous_arc(
																		hap_NO => $hap_NO,
																		last_arc => \%last_reg,
																		this_arc => $this_reg_Href,
																		PNC_NO => $PNC_count,
																		record => 1
																	);
						}
						# find the proper layer to draw
						if( $ts_layer_NO == 0){ # have not find any previous layer
							my ($last_ed_rad, $last_ori, $last_layer_NO) = @last_reg{('ed_rad', 'ori', 'layer_NO')};
							my ($ll_st_rad, $ll_ed_rad) = ($last_ed_rad, $st_rad);
							my $link_ori = ($ll_ed_rad > $ll_st_rad) ? 1 : -1;
							if(   ( $link_ori == $last_ori && $ori == $last_ori ) ## case NO.5, NO.7
								 || ( $link_ori != $last_ori && $ori == $link_ori ) ## case NO.2, NO.4
							){
								## link_line is with the this_seg
								my $MergeSeg_st_rad = $ll_st_rad;
								my $MergeSeg_ed_rad = $ed_rad;
								# try to find proper layer for merged segment
								for my $cand_layer_NO ( 1 .. 1000 ){ # could not touch 1000
									my $path_extend_ratio = 1;
									if(    ( $cand_layer_NO == $last_layer_NO )
											&& ( $link_ori == $last_ori && $ori == $last_ori ) ## case NO.5, NO.7, try to draw with last_seg
										){
											$path_extend_ratio = -1;
									}
									# ok at this cand_layer?
									if(     $seg_OB->ok_at_given_layer(   # merge seg OK?
																			layer_NO => $cand_layer_NO,
																			st_rad => $MergeSeg_st_rad,
																			ed_rad => $MergeSeg_ed_rad,
																			extend => $V_Href->{path_extend} * $path_extend_ratio
																	)
											&&  $seg_OB->ok_at_given_layer(   # this_seg ok?
																			layer_NO => $cand_layer_NO,
																			st_rad => $st_rad,
																			ed_rad => $ed_rad,
																			extend => $V_Href->{path_extend}
																	)
										){ # yes
										$ts_layer_NO = $cand_layer_NO;
										# record link_line at last_seg's layer
										$seg_OB->record_layer_occu(
																	layer_NO => $cand_layer_NO,
																	st_rad => $MergeSeg_st_rad,
																	ed_rad => $MergeSeg_ed_rad
																);
										# finish
										last;
									}
									elsif( $cand_layer_NO eq 1000 ){
										warn_and_exit "<ERROR>\tcan not find proper layer for segment $refseg:$st_pos-$ed_pos\n"
															."\tMergeSeg_st_rad:$MergeSeg_st_rad, MergeSeg_ed_rad:$MergeSeg_ed_rad, st_rad:$st_rad, ed_rad:$ed_rad\n";
									}
								}
								# draw link_line, and this seg
								if($draw_bool){
									# draw arc for link_line and this_seg
									#                     layer_NO    st_rad      ed_rad       ori      dash PNC_number
									my $l_arc_Aref = [$ts_layer_NO,  $ll_st_rad, $ll_ed_rad, $link_ori,  1,  0]; ## link_line
									my $t_arc_Aref = [$ts_layer_NO,  $st_rad,    $ed_rad,    $ori,       0,  $PNC_count]; ## this_seg
									for my $arc_Aref ( $l_arc_Aref, $t_arc_Aref ){
										$seg_OB->draw_path_arc(
																layer_NO => $arc_Aref->[0],
																st_rad => $arc_Aref->[1],
																ed_rad => $arc_Aref->[2],
																ori => $arc_Aref->[3], # 1 for forward, -1 for reversed
																sam_id => $sam_id,
																hap_NO => $hap_NO,
																dash_bool => $arc_Aref->[4],
																arc_NO => $arc_Aref->[5],
																arc_name => ($arc_Aref->[5])?$reg_name:''
															);
									}
									# if needed, radian line
									if( $ts_layer_NO != $last_layer_NO ){
										$seg_OB->draw_radian_line(
																	rad => $ll_st_rad,
																	st_layer_NO => $last_layer_NO,
																	ed_layer_NO => $ts_layer_NO,
																	color => $color
																);
									}
								}
							}
							elsif(   ( $link_ori == $last_ori && $ori != $last_ori ) ## case NO.6, NO.8
									|| ( $link_ori != $last_ori && $ori == $last_ori ) ## case NO.1, NO.3
							){
								# get the perc of link_line pos
								my $ll_st_perc = $seg_OB->rad2perc( rad => $ll_st_rad );
								my $ll_ed_perc = $seg_OB->rad2perc( rad => $ll_ed_rad );
								# specially for circular virus loop junction
								my $virus_loop =    ( $seg_NO == 0 && !$V_Href->{viral_linear} ) # this is seg from circular virus
												 && (   ($ll_st_perc <= 1  && $ll_ed_perc >= 99)    # virus loop junc, head->tail
													 || ($ll_st_perc >= 99 && $ll_ed_perc <= 1 ) ); # virus loop junc, tail->head
								# in all used layers, and the next new layer, find flanking layers to allocate
								my $ll_layer_NO = 0;
								for my $cand_layer_NO ( 1 .. 1000 ){ # could not touch 1000
									my $ll_path_extend_ratio = ($cand_layer_NO == $last_layer_NO)? -1 : 1;
									## layer ok_bool
									my $ll_ok_bool = 0; ## link_line(ll)
									my $ts_ok_bool = 0; ## this_seg(ts)
									## only when not virus_loop case, check ll
									## is link_line(ll) at this layer ok?
									unless( $virus_loop || $ll_st_rad == $ll_ed_rad ){
										$ll_ok_bool = $seg_OB->ok_at_given_layer(
																					layer_NO => $cand_layer_NO,
																					st_rad => $ll_st_rad,
																					ed_rad => $ll_ed_rad,
																					extend => $ll_path_extend_ratio * $V_Href->{path_extend},
																					arc_len_free => 1
																				);
									}
									## is this_seg(ts) at flanking layer ok[other_cases], or at this layer ok[virus_loop]?
									my @shift = ( $virus_loop || $ll_st_rad == $ll_ed_rad ) ? (0) : (-1, 1);
									my $layer_NO_for_ts;
									for my $shift (@shift){
										$layer_NO_for_ts = $cand_layer_NO + $shift;
										next if( $layer_NO_for_ts < 1 );
										$ts_ok_bool = $seg_OB->ok_at_given_layer(
																					layer_NO => $layer_NO_for_ts,
																					st_rad => $st_rad,
																					ed_rad => $ed_rad,
																					extend => $V_Href->{path_extend}
																				);
										# print "$shift\t$ts_ok_bool\n" if($shift == -1);
										last if( $ts_ok_bool );
									}
									# harvest
									if(    $ts_ok_bool
										&& ( ($virus_loop || $ll_st_rad == $ll_ed_rad) || $ll_ok_bool)
										){
										$ll_layer_NO = $cand_layer_NO;
										$ts_layer_NO = $layer_NO_for_ts;
										last;
									}
								}
								# record layers' usage
								#                     layer_NO       st_rad      ed_rad
								my $ll_layer_Aref = [$ll_layer_NO, $ll_st_rad, $ll_ed_rad]; ## link_line's layer
								my $ts_layer_Aref = [$ts_layer_NO, $st_rad,    $ed_rad]; ## this_seg's layer
								my @layer_Aref = ($ll_layer_Aref, $ts_layer_Aref);
								shift @layer_Aref if( $virus_loop || $ll_st_rad == $ll_ed_rad );
								for my $layer_Aref ( @layer_Aref ){
									$seg_OB->record_layer_occu(
																layer_NO => $layer_Aref->[0],
																st_rad => $layer_Aref->[1],
																ed_rad => $layer_Aref->[2]
															);
								}
								# draw link_line, and this seg
								if($draw_bool){
									# draw arc for link_line and this_seg
									#                    layer_NO      st_rad     ed_rad       ori     dash PNC_number
									my $ll_arc_Aref = [$ll_layer_NO, $ll_st_rad, $ll_ed_rad, $link_ori,  1,  0]; ## link_line(ll)
									my $ts_arc_Aref = [$ts_layer_NO, $st_rad,    $ed_rad,    $ori,       0,  $PNC_count]; ## this_seg(ts)
									my @arc_Arefs = ($ll_arc_Aref, $ts_arc_Aref);
									shift @arc_Arefs if( $virus_loop || $ll_st_rad == $ll_ed_rad );
									for my $arc_Aref ( @arc_Arefs ){
										$seg_OB->draw_path_arc(
																layer_NO => $arc_Aref->[0],
																st_rad => $arc_Aref->[1],
																ed_rad => $arc_Aref->[2],
																ori => $arc_Aref->[3], # 1 for forward, -1 for reversed
																sam_id => $sam_id,
																hap_NO => $hap_NO,
																dash_bool => $arc_Aref->[4],
																arc_NO => $arc_Aref->[5],
																arc_name => ($arc_Aref->[5])?$reg_name:''
															);
									}
									# radian line between link_line and seg(s)
									#                         rad      st_layer_NO     ed_layer_NO
									my $lt_radian_Aref = [$ll_ed_rad, $ll_layer_NO,   $ts_layer_NO]; # between link_line(l) and this_seg(s)
									my $Ll_radian_Aref = [$ll_st_rad, $last_layer_NO, $ll_layer_NO]; # between last_seg(L) and link_line(l)
									my @radian_Arefs = ($virus_loop) ? () : ($lt_radian_Aref, $Ll_radian_Aref);
									for my $radian_Aref ( @radian_Arefs ){
										next if( $radian_Aref->[1] == $radian_Aref->[2] );
										$seg_OB->draw_radian_line(
																	rad => $radian_Aref->[0],
																	st_layer_NO => $radian_Aref->[1],
																	ed_layer_NO => $radian_Aref->[2],
																	color => $color
																);
									}
								}
							}
						}
					}
					# update
					$this_reg_Href->{layer_NO} = $ts_layer_NO;
					# record arc details in sample object
					if( !$V_Href->{show_arc_reg} ){
						$sam_OB->load_path_arc_info(
														hap_NO => $hap_NO,
														last_arc => \%last_reg,
														this_arc => $this_reg_Href,
														PNC_text => $PNC_count
												);
					}
					# update last_seg hash
					%last_reg = (seg_NO=>$seg_NO, st_rad=>$st_rad, ed_rad=>$ed_rad, ori=>$ori, layer_NO=>$ts_layer_NO);
					stout_and_sterr "[INFO]\tshow seg_NO.$seg_NO ($refseg:$st_pos-$ed_pos) at layer_NO.$ts_layer_NO for hap_NO.$hap_NO of sample $sam_id\n";
				}
			}
			# warn
			stout_and_sterr "[INFO]\tDraw $sam_id NO.$hap_NO Local Haplotype Path ok.\n";
		}
	}

	# draw PNC(s) and arc_end_dots
	if( $draw_bool ){
		&makeup_PNC_ArcDot;
	}

	# update segment's outer_radius
	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		# update
		$seg_OB->{outer_radius} = $seg_OB->layer2radius( layer_NO => $seg_OB->{nowlayer} ) + 0.75 * $V_Href->{seg_path_space};
	}

	# make up segment's haplotype path sector
	&makeup_segment_sector( draw_bool => $draw_bool );
}

#--- draw PNC(s) and arc_end_dots, to be the top-svg-layer
sub makeup_PNC_ArcDot{

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};
	my $Sample_Href = $V_Href->{Sample_Href};

	# deal each segment
	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		# Arc end dots
		$seg_OB->draw_arc_end_dot(%{$_->{draw_opt}}) for @{$seg_OB->{svg_ele}->{ArcDot}};
		delete $seg_OB->{svg_ele}->{ArcDot};
		# PNC
		for my $PNC_info_Href ( @{$seg_OB->{svg_ele}->{PNC}} ){
			if( !$V_Href->{show_arc_reg} ){
				my $sam_id = $PNC_info_Href->{sam_id};
				my $hap_NO = $PNC_info_Href->{hap_NO};
				my $sam_OB = $Sample_Href->{$sam_id};
				# find characters for group PNCs (degenerate), or return the original PNC_NO
				$PNC_info_Href->{draw_opt}->{text_number}
							= $sam_OB->find_PNC_deg_char(
															hap_NO=>$hap_NO,
															PNC_NO=>$PNC_info_Href->{draw_opt}->{text_number}
														);
			}
			# draw
			$seg_OB->draw_PNC_on_arc( %{$PNC_info_Href->{draw_opt}} );
		}
		# sweep
		delete $seg_OB->{svg_ele}->{PNC};
	}
}

#--- make up segment's sector
sub makeup_segment_sector{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};

	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		# draw sector
		$seg_OB->makeup_sector( draw_bool => $draw_bool );
	}
}

#--- show sample legend
sub show_figure_legends{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $x = $parm{x};
	my $y = $parm{y};

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};

	#------------------#
	# sample path line #
	#------------------#
	for my $sam_id (sort keys %$Sample_Href){
		my $sam_OB = $Sample_Href->{$sam_id};
		# show sample_id
		my $text_size_Aref
			 = &show_text_in_line(
								svg_obj => $V_Href->{circos_OB},
								text_x => $x + $V_Href->{legend_Samline_len} * 1.5,
								text_y => $y,
								text => $sam_id,
								font_family => $V_Href->{legend_text_fontfam},
								font_size => $V_Href->{legend_SamId_fontsize},
								text_col => 'black',
								text_anchor => 'start',
								height_adjust => 0
							);
		# show path line
		$V_Href->{circos_OB}->line(
									x1 => $x,
									y1 => $y - $text_size_Aref->[0] * 0.5,
									x2 => $x + $V_Href->{legend_Samline_len},
									y2 => $y - $text_size_Aref->[0] * 0.5,
									stroke => $sam_OB->{pathcolour}->{1}, # defaults as first haplo-path
									'stroke-width' => $V_Href->{legend_Samline_width},
									'stroke-linecap' => 'round',
									opacity => 1
								 );
		# change y locus
		$y += $text_size_Aref->[0] * 1.5;
	}

	#------------------#
	# degernated PNCs  #
	#------------------#
	# drawing parameters
	my $y_span = 2.5 * $V_Href->{PNC_radius};
	# for degernated PNCs
	for my $sam_id (sort keys %$Sample_Href){
		my $sam_OB = $Sample_Href->{$sam_id};
		for my $hap_NO (sort {$a<=>$b} keys %{$sam_OB->{PNC_deg_char}}){
			my $circle_col = $sam_OB->{pathcolour}->{$hap_NO};
			my $text_col = ($circle_col =~ /\b$COLOR_DB{dark_col}/)?'white':'black';
			my $degchar_hap_Href = $sam_OB->{PNC_deg_char}->{$hap_NO};
			for my $PNC_deg_char ( @{$V_Href->{PNC_deg_char}} ){
				next unless( exists($degchar_hap_Href->{$PNC_deg_char}) );
				my $PNC_x = $x;
				my @PNC_text = ($PNC_deg_char);
				push @PNC_text, $_ for @{$degchar_hap_Href->{$PNC_deg_char}};
				# draw
				for my $i (0 .. $#PNC_text){
					my $PNC_text = $PNC_text[$i];
					my $text_font_size = FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype->PNC_text_size( PNC_text => $PNC_text );
					draw_a_ellipse(
									svg_obj => $V_Href->{circos_OB},
									cx => $PNC_x,
									cy => $y,
									radius => $V_Href->{PNC_radius},
									fill_color => $circle_col,
									text_fill => $PNC_text,
									text_font_size => $text_font_size,
									text_col => $text_col,
									font_size_auto_adjust => 0
								);
					# change x locus
					$PNC_x += 2.5 * $V_Href->{PNC_radius};
				}
				# change y locus
				$y += $y_span;
			}
		}
	}

	#--------------------#
	# duplicated contigs #
	# in local haplotype #
	#--------------------#
	$y += $y_span;
	# drawing parameters
	my $x_gap = $y_span * 0.5;
	my $strline_biend = $V_Href->{PNC_radius};
	# for duplicated contigs
	for my $sam_id (sort keys %$Sample_Href){
		my $sam_OB = $Sample_Href->{$sam_id};
		for my $hap_NO (sort {$a<=>$b} keys %{$sam_OB->{haplotype}}){
			my $circle_col = $sam_OB->{pathcolour}->{$hap_NO};
			my $text_col = ($circle_col =~ /\b$COLOR_DB{dark_col}/)?'white':'black';
			my $hap_Href = $sam_OB->{haplotype}->{$hap_NO};
			for my $ctg_NO (sort {$a<=>$b} keys %$hap_Href){
				my $repeat_time = $hap_Href->{$ctg_NO}->[0];
				next if($repeat_time <= 1);
				my $region_Aref = $hap_Href->{$ctg_NO}->[1];
				# show repeat time text
				my $rep_text = "${repeat_time}X:";
				my $text_size_Aref
					 = &show_text_in_line(
										svg_obj => $V_Href->{circos_OB},
										text_x => $x,
										text_y => $y,
										text => $rep_text,
										font_family => $V_Href->{legend_text_fontfam},
										font_size => $V_Href->{legend_SamId_fontsize},
										text_col => 'black',
										text_anchor => 'end',
										height_adjust => 0
									);
				# string middle line
				my $line_y = $y - $text_size_Aref->[0] * 0.5,
				my $line_x = $x + $x_gap;
				my $line_len = 2 * $strline_biend + scalar(@$region_Aref) * ($V_Href->{PNC_radius} * 3);
				$V_Href->{circos_OB}->line(
											x1 => $line_x,
											y1 => $line_y,
											x2 => $line_x + $line_len,
											y2 => $line_y,
											stroke => $circle_col,
											'stroke-width' => $V_Href->{legend_Samline_width},
											'stroke-linecap' => 'round',
											opacity => 1
										 );
				# draw PNC
				my $PNC_x = $line_x + $strline_biend + $V_Href->{PNC_radius};
				for my $reginfo_Aref (@$region_Aref){
					my $PNC_NO = $reginfo_Aref->[-1];
					## draw arrow
					my $arrow_cx = $PNC_x + sqrt(3) * $V_Href->{PNC_radius} * 0.5;
					draw_a_triangle(
										svg_obj => $V_Href->{circos_OB},
										x => $arrow_cx,
										y => $line_y,
										bottom_side_len => $V_Href->{PNC_radius} * 1.8,
										fill_color => 'white',
										rotate_degree => 90
									 );
					## draw the circle
					my $text_font_size = FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype->PNC_text_size( PNC_text => $PNC_NO );
					draw_a_ellipse(
									svg_obj => $V_Href->{circos_OB},
									cx => $PNC_x,
									cy => $line_y,
									radius => $V_Href->{PNC_radius},
									fill_color => $circle_col,
									text_fill => $PNC_NO,
									text_font_size => $text_font_size,
									text_col => $text_col,
									font_size_auto_adjust => 0
								);
					## change $PNC_x
					$PNC_x += $V_Href->{PNC_radius} * 3;
				}
				# change y locus
				$y += $y_span;
			}
		}
	}

	#---------------#
	# feature icons #
	#---------------#
	$y += $y_span;
	my $fet_basic_Href = $Feature_Basic_Href->{feature};
	for my $fet_type ('cause', 'effect'){
		my $icon_shape = $Feature_Basic_Href->{type}->{$fet_type}->{icon};
		my $fill_color = $Feature_Basic_Href->{type}->{$fet_type}->{colour};
		for my $fet_id (sort keys %$fet_basic_Href){
			next if(    $fet_basic_Href->{$fet_id}->{type} ne $fet_type
						 || $fet_basic_Href->{$fet_id}->{show_bool} == 0    );
			#--- show feature notes with underlined letter
			my $text_x = $x;
			my $text_y = $y;
			my $icon_fulltext = $fet_basic_Href->{$fet_id}->{show_text};
			my $icon_shortname = $fet_basic_Href->{$fet_id}->{short_name};
			my $icon_shortname_index = $fet_basic_Href->{$fet_id}->{sn_index};
			## get the text index
			my $index = 0;
			my $text_features_Href;
			for my $i ( 0 .. length($icon_fulltext)-1 ){
				if( substr($icon_fulltext, $i, 1) =~ /$icon_shortname/i){
					$index++;
					if( $index == $icon_shortname_index ){
						$text_features_Href = {
												font_family => {"$V_Href->{legend_text_fontfam} Bold"=>[$i+1]},
												underline => {$V_Href->{legend_text_underline_width}=>[$i+1]}
												};
						last;
					}
				}
			}
			## show feature legend text
			my $legend_size_Aref = 
				show_text_in_line(
									svg_obj => $V_Href->{circos_OB},
									text_x => $text_x,
									text_y => $text_y,
									text => $icon_fulltext,
									font_family => $V_Href->{legend_text_fontfam},
									font_size => $V_Href->{legend_SamId_fontsize},
									text_col => 'black',
									text_anchor => 'start',
									features => $text_features_Href
								);
			#--- show featture icon
			my $icon_x = $text_x - $V_Href->{legend_icon_text_gap} - $V_Href->{icon_edge_len} * 0.5;
			my $icon_y = $text_y - $legend_size_Aref->[0] * 0.5;
			FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->draw_a_feature_icon(
																					x => $icon_x,
																					y => $icon_y,
																					icon_shape => $icon_shape,
																					fill_color => $fill_color,
																					feature_text => $icon_shortname,
																					rotate_degree => 0
																				);
			# change y locus
			$y += $y_span * 1.25;
		}
	}

	#-------------#
	# gene strand #
	#-------------#
	$y += $y_span;
	my $gene_legend_Href = { plus=>$V_Href->{gene_col_plus}, minus=>$V_Href->{gene_col_minus} };
	for my $strd ('plus', 'minus'){
		my $text_x = $x;
		my $text_y = $y;
		# show strands
		my $legend_size_Aref = 
			show_text_in_line(
								svg_obj => $V_Href->{circos_OB},
								text_x => $text_x,
								text_y => $text_y,
								text => $strd,
								font_family => $V_Href->{legend_text_fontfam},
								font_size => $V_Href->{legend_SamId_fontsize},
								text_col => 'black',
								text_anchor => 'end'
							);
		# fill col
		my $strd_col = $gene_legend_Href->{$strd};
		my $icon_x = $text_x + $V_Href->{legend_icon_text_gap} + $V_Href->{icon_edge_len} * 0.5;
		my $icon_y = $text_y - $legend_size_Aref->[0] * 0.5;
		FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->draw_a_feature_icon(
																				x => $icon_x,
																				y => $icon_y,
																				icon_shape => 'square',
																				fill_color => $strd_col,
																				feature_text => '',
																				rotate_degree => 0
																			);
		# change y locus
		$y += $y_span * 1.25;
	}

	#-------------#
	# copy number #
	#-------------#
	for my $cn_or_exp ('CN', 'EXP'){
		$y += $y_span;
		my ($linear_text, $gradient_text);
		my ($log_value);
		my ($bg_color, $fill_color);
		my ($min_value, $max_value);
		my ($svg_grad_pattern_id);
		if( $cn_or_exp eq 'CN' ){
			$linear_text = 'CN';
			$log_value = $V_Href->{dna_cn_log};
			$gradient_text = "log[$log_value](CN)";
			$bg_color = $V_Href->{CN_bg_color};
			$fill_color = "RGB($V_Href->{CN_color_rgb})";
			$min_value = $V_Href->{CN_min};
			$max_value = $V_Href->{CN_max};
			$svg_grad_pattern_id = 'cn-grad';
		}
		else{
			$linear_text = 'EXP';
			$log_value = $V_Href->{exp_log};
			$gradient_text = "log[$log_value](EXP)";
			$bg_color = $V_Href->{EXP_bg_color};
			$fill_color = "RGB($V_Href->{EXP_color_rgb})";
			$min_value = $V_Href->{EXP_min};
			$max_value = $V_Href->{EXP_max};
			$svg_grad_pattern_id = 'exp-grad';
		}
		# text
		my $text_x = $x;
		my $text_y = $y;
		my $text = ($log_value) ? $gradient_text : $linear_text;
		my $legend_size_Aref = 
			show_text_in_line(
								svg_obj => $V_Href->{circos_OB},
								text_x => $text_x,
								text_y => $text_y,
								text => $text,
								font_family => $V_Href->{legend_text_fontfam},
								font_size => $V_Href->{legend_SamId_fontsize},
								text_col => 'black',
								text_anchor => 'end'
							);
		# non-overlap
		my $icon_x = $text_x + $V_Href->{legend_icon_text_gap} + $V_Href->{icon_edge_len} * 0.5;
		my $icon_y = $text_y - $legend_size_Aref->[0] * 0.5;
		FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->draw_a_feature_icon(
																				x => $icon_x,
																				y => $icon_y,
																				icon_shape => 'square',
																				fill_color => $bg_color,
																				feature_text => '',
																				rotate_degree => 0
																			);
		# color gradient
		my $rect_x = $icon_x + $V_Href->{legend_icon_text_gap};
		my $rect_y = $icon_y - $V_Href->{icon_edge_len} * 0.5;
		my $CN_patt = $V_Href->{circos_OB}->pattern();
		my $CN_grad = $CN_patt->gradient( id=>$svg_grad_pattern_id,
											x1=>"0%", y1=>"0%", x2=>"100%", y2=>"0%",
											spreadMethod=>"pad", gradientUnits=>"objectBoundingBox"
										);
		$CN_grad->stop(offset=>"0%", style=>{'stop-color'=>'white','stop-opacity'=>1});
		$CN_grad->stop(offset=>"100%", style=>{'stop-color'=>$fill_color,'stop-opacity'=>1});
		$V_Href->{circos_OB}->rect( x=>$rect_x, y=>$rect_y,
									width=>$V_Href->{legend_CnExp_grad_band_width}, height=>$V_Href->{icon_edge_len},
									'stroke-width'=>0.5, stroke=>'black', fill=>"url(#$svg_grad_pattern_id)"
									);
		# axis
		my $axis_y = $rect_y + $V_Href->{icon_edge_len} + $V_Href->{legend_grad_band_axis_gap}; # gap
		## draw axis line
		$V_Href->{circos_OB}->line( x1=>$rect_x, y1=>$axis_y,
									x2=>$rect_x+$V_Href->{legend_CnExp_grad_band_width}, y2=>$axis_y,
									stroke=>'black', 'stroke-linecap'=>'round',
									'stroke-width'=>$V_Href->{legend_grad_band_axis_width}
									);
		## prepare axis notes
		my $min_cord = sprintf("%.1f", &cn_exp_convert( value => $min_value, log => $log_value )) * 1.0;
		my $max_cord = sprintf("%.1f", &cn_exp_convert( value => $max_value, log => $log_value )) * 1.0;
		my %axis_note = ($min_cord=>1);
		if( $log_value == 0){
			my $step = int( ($max_cord - $min_cord) / $V_Href->{legend_grad_band_axis_step_NO} );
			$axis_note{int($min_cord)+$_*$step} = 1 for (1 .. 3);
		}
		else{
			$axis_note{0} = 1; # (x)^^0 = 1;
			my $step = sprintf("%.2f", $max_cord / $V_Href->{legend_grad_band_axis_step_NO});
			$axis_note{int($_*$step)} = 1 for (1 .. 3);
		}
		$axis_note{$max_cord} = 1;
		## draw axis pos
		for my $pos_note (sort {$a<=>$b} keys %axis_note){
			my $note_x = $rect_x + ($pos_note-$min_cord) / ($max_cord-$min_cord) * $V_Href->{legend_CnExp_grad_band_width};
			my $note_y = $axis_y+$V_Href->{legend_grad_band_axis_height};
			## short vertical line
			$V_Href->{circos_OB}->line( x1=>$note_x, y1=>$axis_y, x2=>$note_x, y2=>$note_y,
										stroke=>'black', 'stroke-linecap'=>'round',
										'stroke-width'=>$V_Href->{legend_grad_band_axis_width}
										);
			## draw text
			show_text_in_line(
								svg_obj => $V_Href->{circos_OB},
								text_x => $note_x,
								text_y => $note_y + $V_Href->{legend_grad_band_axis_fontsize},
								text => $pos_note,
								font_family => $V_Href->{legend_text_fontfam},
								font_size => $V_Href->{legend_grad_band_axis_fontsize},
								text_col => 'black',
								text_anchor => 'middle'
							);
		}

		# change y locus
		$y += $y_span * 1.25;
	}
}

#--- cn or exp value convert
sub cn_exp_convert{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $value = $parm{value};
	my $log_value = $parm{log};

	if( $log_value != 0 ){
		$value ||= 0.1;
		return log($value)/log($log_value);
	}
	else{
		return $value;
	}
}

#--- show segments' names
sub show_segment_name{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};

	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		# features
		my $seg_mid_rad = $seg_OB->{start_rad} + $seg_OB->{rad_size} * 0.5;
		&normalize_radian(radian_Sref=>\$seg_mid_rad);
		my $text_to_center = ($seg_mid_rad >= $V_Href->{axis_top_rad}->[0] && $seg_mid_rad <= $V_Href->{axis_top_rad}->[1]) ? 'top':'toe';
		my ($seg_show_name) = ($seg_OB->{name} =~ /^([^:]+):{0,1}/);
		my $seg_col = ($seg_NO>0) ? $V_Href->{host_segcol} : $V_Href->{virus_segcol};
		# draw band for seg_name
		draw_circle_seg(
							svg_obj => $V_Href->{circos_OB},
							cx => $V_Href->{cx},
							cy => $V_Href->{cy},
							start_rad => $seg_OB->{start_rad},
							rad_size => $seg_OB->{rad_size},
							inner_radius => $seg_OB->{outer_radius},
							outer_radius => $seg_OB->{outer_radius} + $V_Href->{segment_name_band_width},
							seg_fill_color => $seg_col,
							seg_boundary_width => $V_Href->{segment_boundary_width},
							text => $seg_show_name,
							font_family => $V_Href->{Font_family_SegName},
							font_size => $V_Href->{segment_name_font_size},
							text_col => 'black',
							text_to_center => $text_to_center,
							draw_bool => $draw_bool
						 );
		# update
		$seg_OB->{outer_radius} += $seg_OB->{segment_width};
	}
}

#--- start svg
sub start_svg{
	# variants
	my $svg_width = $V_Href->{backend_width};
	my $svg_height = $V_Href->{backend_height};

	$V_Href->{circos_OB} = SVG->new( width => $svg_width,
						height => $svg_height,
						author => $AUTHOR,
						'author-mail' => $EMAIL );
						# author => $V_Href->{Auther},
						# 'author-mail' => $V_Href->{Email} );

	# background
	draw_a_parallelogram(
							svg_obj => $V_Href->{circos_OB},
							x => $svg_width/2,
							y => $svg_height/2,
							head_bottom_side_len => $svg_width,
							left_right_side_len => $svg_height,
							fill_color => 'white',
							boundary_color => 'black',
							usage_print => 0
						);
}

#--- svg file
sub output_svg_file{

	# variants
	my $svgfile = $V_Href->{output_SVG_file};

	open OUT,">$svgfile" || die $!;
	print OUT $V_Href->{circos_OB}->xmlify;
}

#--- 
1; ## tell the perl script the successful access of this module.
