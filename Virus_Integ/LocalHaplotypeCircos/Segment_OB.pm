package FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB;

use strict;
use warnings;
use List::Util qw/min max sum/;
use Math::Trig;
use SVG;
use BioFuse::Util::GZfile qw/Try_GZ_Read/;
use BioFuse::Visual::SVG_Util::Font qw/show_text_on_arc/;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use BioFuse::Util::Interval qw/Get_Two_Seg_Olen/;
use BioFuse::Visual::SVG_Util::RadSys qw/$PI $deg2rad get_coordinate_on_circle normalize_radian draw_an_arc/;
use BioFuse::Visual::SVG_Util::RadSysEle qw/draw_circle_seg/;
use BioFuse::Visual::SVG_Util::Color qw/%COLOR_DB/;
use BioFuse::Visual::SVG_Util::RectSysEle qw/draw_a_parallelogram draw_a_triangle draw_a_ellipse/;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype;
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
$VERSION = "0.05";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						pos2rad
						which_layer_to_draw
						ok_at_given_layer
						record_layer_occu
						rad2perc
						layer2radius
						enable_draw
						draw_path_arc
						draw_PNC_on_arc
						draw_arc_end_dot
						draw_radian_line
						all_layer_NOs
						load_pos_features
						arrange_features_icon
						draw_features_icon
						draw_icon_stick
						draw_signle_feature_icon
						draw_a_feature_icon
						load_one_transcript
						draw_genes
						gene_at_which_layer
						draw_CN_band
						get_graded_RGB_color
						makeup_sector
						draw_EXP_band
						gather_exp_depth_rad_win
						get_given_sample_rna_depth
						deal_rna_depth_buffer
						calculate_mean_depth_of_rad_win
						pos2radwin
					 /;

#--- structure of object
# seg_OB -> NO = seg_NO
# seg_OB -> sample = []
# seg_OB -> ref_seg
# seg_OB -> lf_pos
# seg_OB -> rt_pos
# seg_OB -> resol -> {$st_pos} = [ed_pos, this_region_radsize]
# seg_OB -> rad_size
# seg_OB -> start_rad
# seg_OB -> nowlayer
# seg_OB -> layeroccu
# seg_OB -> draw_bool
# seg_OB -> name
# seg_OB -> inner_radius
# seg_OB -> outer_radius
# seg_OB -> segment_width
# seg_OB -> svg_ele -> {PNC / ArcDot}
# seg_OB -> pos_feature -> {pos} -> { sam_OB=>'', rad_ori=>'', rad_now=>'', features=>[ feature_1, feature_2, ... ] }
# seg_OB -> c -> {$gname} -> {$tname} -> { refseg=>'', strand=>'', gname=>'', tname=>'', lf_pos=>'', rt_pos=>'', ex_len=>'', ex_lfp=>''  }
# seg_OB -> exp_rad_win -> {$sam_id} -> {$rad_win_NO} -> {mean_depth=>, pos_depth=>{pos=>depth}} ## will del pos_depth

#--- construction of object
sub new{
	shift;
	my %parm = @_;

	my $seg_OB = {};
	$seg_OB->{NO} = $parm{NO};
	$seg_OB->{ref_seg} = $parm{refseg};
	$seg_OB->{lf_pos} = $parm{lfpos};
	$seg_OB->{rt_pos} = $parm{rtpos};
	$seg_OB->{rad_size} = $parm{rad_size};
	$seg_OB->{start_rad} = $parm{start_rad};
	$seg_OB->{name} = $parm{name};
	$seg_OB->{inner_radius} = $parm{inner_radius};
	$seg_OB->{segment_width} = $parm{inner_radius} * $parm{InRadToSegWidRatio};
	$seg_OB->{outer_radius} = $seg_OB->{inner_radius} + $seg_OB->{segment_width}; # already count in the segment_band
	$seg_OB->{nowlayer} = 0; # initially, no layer has been used
	$seg_OB->{layeroccu} = {};
	$seg_OB->{draw_bool} = 0; # initially, no haplotype linkage. temporarily, not used.
	# relevant sample
	@{ $seg_OB->{sample} } = @{ $parm{sampleid_Aref} };
	# sub-region with different resolutions
	my $linear_pixel = $parm{linear_pixel};
	for my $Aref (@{$parm{resol_AAref}}){
		my $st_pos = $Aref->[0];
		my $ed_pos = $Aref->[1];
		my $resol = $Aref->[2];
		my $sub_radsize = ( ($ed_pos-$st_pos+1) / $resol / $linear_pixel ) * $parm{rad_size};
		$seg_OB->{resol}->{$st_pos} = [ $ed_pos, $sub_radsize ];
	}

	bless($seg_OB);
	return $seg_OB;
}

#--- give a position, return its radian angle
sub pos2rad{
	my $seg_OB = shift;
	my %parm = @_;
	my $refseg = $parm{refseg};
	my $pos = $parm{pos};

	# examine
	if(    $refseg ne $seg_OB->{ref_seg}
		|| $pos < $seg_OB->{lf_pos}
		|| $pos > $seg_OB->{rt_pos}
		){
			warn_and_exit "<ERROR>\tPos ($refseg:$pos) is out of Segment NO.$seg_OB->{NO}.";
	}

	# find the sub_region where this pos locates
	my $pos_rad = $seg_OB->{start_rad};
	for my $st_pos (sort {$a<=>$b} keys %{$seg_OB->{resol}}){
		my ($ed_pos, $sub_rad) = @{ $seg_OB->{resol}->{$st_pos} };
		if( $pos >= $st_pos && $pos <= $ed_pos){
			$pos_rad += $sub_rad * ($pos-$st_pos+1) / ($ed_pos-$st_pos+1);
			return $pos_rad;
		}
		else{
			$pos_rad += $sub_rad;
		}
	}

	# warn but never meets
	warn_and_exit "<ERROR>\tCannot find radius for Pos ($refseg:$pos) on Segment NO.$seg_OB->{NO}.";
}

#--- determine which layer to draw the radius
sub which_layer_to_draw{
	my $seg_OB = shift;
	my %parm = @_;
	my $st_rad = $parm{st_rad};
	my $ed_rad = $parm{ed_rad};
	my $min_layer = $parm{min_layer} || 1;
	my $extend = $parm{extend};
	my $record = $parm{record};

	# find the layer proper for this reg
	my $ok_layer_NO = 0;
	for my $layer_NO ( 1 .. 1000 ){ # never touch 1000
		# only wants outer layers
		next if($layer_NO < $min_layer);
		# test this layer subregions
		if( $seg_OB->ok_at_given_layer( layer_NO=>$layer_NO, st_rad=>$st_rad, ed_rad=>$ed_rad, extend=>$extend ) ){
			# find this layer is porper
			$ok_layer_NO = $layer_NO;
			last;
		}
	}

	if( $ok_layer_NO == 0){ # not found
		# use the new layer
		$ok_layer_NO = $seg_OB->{nowlayer} + 1;
	}
	if($record){ # record the usage of this ok_layer of local haplotype path
		$seg_OB->record_layer_occu( layer_NO=>$ok_layer_NO, st_rad=>$st_rad, ed_rad=>$ed_rad ) ;
	}

	return $ok_layer_NO;
}

#--- test given layer ok to given sub-region
sub ok_at_given_layer{
	my $seg_OB = shift;
	my %parm = @_;
	my $layer_NO = $parm{layer_NO};
	my $st_rad = $parm{st_rad};
	my $ed_rad = $parm{ed_rad};
	my $extend = $parm{extend};
	# my $debug = $parm{debug} || 0;
	my $arc_len_free = $parm{arc_len_free} || 0;

	# variants
	my $min_arc_len = $V_Href->{min_arc_len};

	# check layer_NO
	if( $layer_NO < 1 ){
		warn_and_exit "Not allowed to judge ok on $layer_NO layer.\n";
	}

	# check arc length
	my $layer_radius = $seg_OB->layer2radius( layer_NO => $layer_NO );
	my $arc_len = $layer_radius * abs($ed_rad-$st_rad);
	if( !$arc_len_free && $arc_len < $min_arc_len){
		return 0;
	}


	# convert rad to percentage
	my $st_prec = $seg_OB->rad2perc( rad => $st_rad );
	my $ed_prec = $seg_OB->rad2perc( rad => $ed_rad );
	($st_prec, $ed_prec) = sort {$a<=>$b} ($st_prec, $ed_prec);

	for my $Aref (@{$seg_OB->{layeroccu}->{$layer_NO}}){
		my ($st_prec2, $ed_prec2) = @$Aref;
		# test overlap
		if(Get_Two_Seg_Olen($st_prec-$extend, $ed_prec+$extend, $st_prec2-$extend, $ed_prec2+$extend, 0.0001)){
			return 0;
		}
	}

	# not overlap, ok to this region
	return 1;
}

#--- to record the occupy of some radius
sub record_layer_occu{
	my $seg_OB = shift;
	my %parm = @_;
	my $layer_NO = $parm{layer_NO};
	my $st_rad = $parm{st_rad};
	my $ed_rad = $parm{ed_rad};

	# convert rad to percentage
	my $st_prec = $seg_OB->rad2perc( rad => $st_rad );
	my $ed_prec = $seg_OB->rad2perc( rad => $ed_rad );
	($st_prec, $ed_prec) = sort {$a<=>$b} ($st_prec, $ed_prec);

	# record
	push @{ $seg_OB->{layeroccu}->{$layer_NO} }, [$st_prec, $ed_prec];
	# and also update the most outer layer
	$seg_OB->{nowlayer} = max( $seg_OB->{nowlayer}, $layer_NO );
}

#--- to convert radius to perctage
sub rad2perc{
	my $seg_OB = shift;
	my %parm = @_;

	return int( ( $parm{rad} - $seg_OB->{start_rad} ) / $seg_OB->{rad_size} * 10000 ) / 100;
}

#--- get radius of given layer
sub layer2radius{
	my $seg_OB = shift;
	my %parm = @_;

	return (  $seg_OB->{inner_radius}
				+ $seg_OB->{segment_width}
				+ $parm{layer_NO} * ($V_Href->{path_width} + $V_Href->{seg_path_space}) );
}

#--- just set draw bool as 1
sub enable_draw{
	my $seg_OB = shift;
	$seg_OB->{draw_bool} = 1;
}

#--- draw a path on given layer of seg
sub draw_path_arc{
	my $seg_OB = shift;
	my %parm = @_;
	my $st_rad = $parm{st_rad};
	my $ed_rad = $parm{ed_rad};
	my $layer_NO = $parm{layer_NO};
	my $ori = $parm{ori}; # 1 for forward, -1 for reversed
	my $sam_id = $parm{sam_id};
	my $hap_NO = $parm{hap_NO};
	my $dash_bool = $parm{dash_bool} || 0;
	my $arc_NO = $parm{arc_NO} || 0;
	my $dot_free = $parm{dot_free} || 0;
	my $arc_name = $parm{arc_name};

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};
	my $show_arc_reg = $V_Href->{show_arc_reg};
	my $color = $Sample_Href->{$sam_id}->{pathcolour}->{$hap_NO};

	# draw this arc
	my $layer_radius = $seg_OB->layer2radius( layer_NO => $layer_NO );
	my ($x_st, $y_st) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$st_rad, radius=>$layer_radius );
	my ($x_ed, $y_ed) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$ed_rad, radius=>$layer_radius );
	my $flag_1 = (abs($st_rad-$ed_rad)>$PI)?1:0; # large arc(1) or small arc(0)
	my $flag_2 = ($ori == 1) ? 1 : 0; # clockwise(1) or anit-clockwise(0)
	my $path = "M$x_st,$y_st,A$layer_radius,$layer_radius,0,$flag_1,$flag_2,$x_ed,$y_ed";
	$V_Href->{circos_OB}->path(
								d => $path,
								stroke => $color,
								'stroke-width' => $V_Href->{path_width},
								'stroke-linecap' => 'round',
								'stroke-dasharray' => (($dash_bool) ? $V_Href->{line_dash} : ' '),
								fill => 'none',
								opacity => $V_Href->{path_line_opacity}
								);

	# draw two ends dot
	if( !$dot_free ){
		for my $dot_Aref ( [$x_st, $y_st], [$x_ed, $y_ed] ){
			my $draw_option_Href = { cx => $dot_Aref->[0], cy => $dot_Aref->[1] };
			# after all arc(s) are drawn, then draw arc_end_dot(s)
			push @{$seg_OB->{svg_ele}->{ArcDot}}, {sam_id=>$sam_id, hap_NO=>$hap_NO, draw_opt=>$draw_option_Href};
		}
	}

	# path_NO circle on arc
	# show arc_reg for resolution adjustment
	if( $arc_NO || $show_arc_reg ){
		my $draw_option_Href = {
								st_rad => $st_rad,
								ed_rad => $ed_rad,
								layer_NO => $layer_NO,
								circle_col => $color,
								ori => $ori,
								text_number => ($show_arc_reg)?$arc_name:$arc_NO
								 };
		# after all arc(s) are drawn, then draw PNC(s)
		push @{$seg_OB->{svg_ele}->{PNC}}, {sam_id=>$sam_id, hap_NO=>$hap_NO, draw_opt=>$draw_option_Href};
	}
}

#--- draw path_NO circle (PNC) on arc
sub draw_PNC_on_arc{
	my $seg_OB = shift;
	my %parm = @_;
	my $st_rad = $parm{st_rad};
	my $ed_rad = $parm{ed_rad};
	my $layer_NO = $parm{layer_NO};
	my $circle_col = $parm{circle_col};
	my $ori = $parm{ori}; # 1 for forward, -1 for reversed
	my $text_number = $parm{text_number};

	# for PNC
	my $PNC_rad = $st_rad + ($ed_rad - $st_rad) * (4 / 10);
	my $layer_radius = $seg_OB->layer2radius( layer_NO => $layer_NO );
	my ($cx_PNC, $cy_PNC) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$PNC_rad, radius=>$layer_radius);
	my $text_col = ($circle_col =~ /\b$COLOR_DB{dark_col}/)?'white':'black';
		 $text_col = 'black' if( $V_Href->{show_arc_reg} );
	my $text_font_size = FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype->PNC_text_size( PNC_text => $text_number );
	# for path direction arrow
	my $multiple_factor = $ori;
	my $arrow_half_span_rad = asin( $V_Href->{PNC_radius} / $layer_radius / 2) * 2;
	my $arrow_center_rad = $PNC_rad + $multiple_factor * $arrow_half_span_rad;
	my ($arrow_cx, $arrow_cy) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$arrow_center_rad, radius=>$layer_radius);
	my $arrow_bottom_sidelen = $V_Href->{PNC_radius} * sqrt(3);
	my $arrow_rotate_rad = $PNC_rad + $multiple_factor * ($arrow_half_span_rad + $PI / 2);
	my $arrow_rotate_degree = ($arrow_rotate_rad / $deg2rad * 100) / 100;

	# draw PNC and show text
	if( length($text_number) != 0){
		# draw arrow
		draw_a_triangle(
							svg_obj => $V_Href->{circos_OB},
							x => $arrow_cx,
							y => $arrow_cy,
							bottom_side_len => $arrow_bottom_sidelen,
							fill_color => 'white',
							rotate_degree => $arrow_rotate_degree
						 );
		## draw the circle
		draw_a_ellipse(
							svg_obj => $V_Href->{circos_OB},
							cx => $cx_PNC,
							cy => $cy_PNC,
							radius => $V_Href->{PNC_radius},
							fill_color => $circle_col,
							text_fill => $text_number,
							text_font_size => $text_font_size,
							text_col => $text_col,
							font_size_auto_adjust => 0
						);
	}
}

#--- draw path_NO circle (PNC) on arc
sub draw_arc_end_dot{
	my $seg_OB = shift;
	my %parm = @_;

	# draw the dot
	draw_a_ellipse(
					svg_obj => $V_Href->{circos_OB},
					cx => $parm{cx},
					cy => $parm{cy},
					radius => $V_Href->{path_width} * 0.5,
					fill_color => 'black'
					);
}

#--- draw radius line (dash)
sub draw_radian_line{
	my $seg_OB = shift;
	my %parm = @_;
	my $rad = $parm{rad};
	my $st_layer_NO = $parm{st_layer_NO} || 0;
	my $ed_layer_NO = $parm{ed_layer_NO};
	my $color = $parm{color};

	my $func_name = 'draw_radian_line';

	# check
	if( !defined($ed_layer_NO) || $ed_layer_NO == 0 ){
		warn_and_exit "<ERROR>\t$st_layer_NO\t$ed_layer_NO\t$seg_OB->{sample}->[0]\n"
							."\t[f($func_name)]: You need to assign ed_layer_NO as integer number lager than 0.\n";
	}

	# points coordinates
	my $st_layer_radius = ($st_layer_NO) ? $seg_OB->layer2radius( layer_NO => $st_layer_NO ) : $seg_OB->{inner_radius};
	my ($x_st, $y_st) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad, radius=>$st_layer_radius );
	my $ed_layer_radius = $seg_OB->layer2radius( layer_NO => $ed_layer_NO );
	my ($x_ed, $y_ed) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad, radius=>$ed_layer_radius );
	# draw dash line
	$V_Href->{circos_OB}->line(
								x1=>$x_st,
								y1=>$y_st,
								x2=>$x_ed,
								y2=>$y_ed,
								stroke=>$color,
								'stroke-width'=>$V_Href->{path_width},
								'stroke-dasharray'=>$V_Href->{line_dash},
								'stroke-linecap'=>'round',
								opacity=>$V_Href->{path_line_opacity}
								);
}

#--- return the seg_OB's all used layer(s)
## return reference of sorted array
sub all_layer_NOs{
	my $seg_OB = shift;
	my @sorted_layer_NOs = sort {$a<=>$b} keys %{$seg_OB->{layeroccu}};
	return \@sorted_layer_NOs;
}

#--- load position features
sub load_pos_features{
	my $seg_OB = shift;
	my %parm = @_;
	my $sam_OB = $parm{sam_OB};
	my $refseg = $parm{refseg};
	my $pos = $parm{pos};
	my $feature_Aref = $parm{feature_Aref};
	my $partner = $parm{partner};

	## seg_NO -> {pos_feature} -> {pos} -> { rad_ori=>'', rad_now=>'', features=>[ feature_1, feature_2, feature_3, ... ] }
	my $rad_ori = $seg_OB->pos2rad( refseg => $refseg, pos => $pos );
	$seg_OB->{pos_feature}->{$pos} = {
										sam_OB => $sam_OB,
										rad_ori => $rad_ori,
										rad_now => $rad_ori
									 };
	push @{$seg_OB->{pos_feature}->{$pos}->{features}}, [$partner, $feature_Aref];
}

#--- adjust the icon's radian to avoid overlap of feature icons
sub arrange_features_icon{
	my $seg_OB = shift;

	# variants
	my $icon_EdgeLen = $V_Href->{icon_edge_len};
	my $icon_tangent_mingap = $V_Href->{icon_tangent_mingap};
	my $icon_loc_adjust_rad_step = $V_Href->{icon_loc_adjust_rad_step};
	my $max_arrange_time = $V_Href->{icon_loc_max_arrange_time};

	my @sort_fet_pos = sort {$a<=>$b} keys %{$seg_OB->{pos_feature}};
	my $start_i = 0;
	my $arrange_time = 0;
	my $comp_radius = $seg_OB->{outer_radius} + sum(@{$V_Href->{icon_line_len}});
	my $icon_min_gaprad = $icon_tangent_mingap / $comp_radius;
	my $one_icon_rad = $icon_EdgeLen / $comp_radius;
	my $min_rad = $seg_OB->{start_rad} + $icon_min_gaprad + $one_icon_rad * 0.5;
	my $max_rad = $min_rad + $seg_OB->{rad_size} - $icon_min_gaprad - $one_icon_rad * 0.5;
	my $icon_line_min_gaprad = $one_icon_rad + $icon_min_gaprad;
	# check whether the seg_rad is enough for features
	my $fet_icon_occupy_rad = scalar(@sort_fet_pos) * $icon_line_min_gaprad + $icon_min_gaprad;
	if( $seg_OB->{rad_size} < $fet_icon_occupy_rad ){
		warn_and_exit "Seg NO.$seg_OB->{NO} rad_size $seg_OB->{rad_size} is less than feature icons need ($fet_icon_occupy_rad).\n";
	}
	ARRANGE: {
		$arrange_time ++;
		warn_and_exit " has done for seg NO.$seg_OB->{NO} reaches the maximum arrangement times for feature_icon ($arrange_time).\n"
			if( $arrange_time > $max_arrange_time );
		# warn "This is the $arrange_time time arrangement operation for seg NO.$seg_OB->{NO}. Start from index $start_i\n"
		#   if($arrange_time % 1000 == 0);
		for my $i ($start_i  .. $#sort_fet_pos){
			my $ts_pos = $sort_fet_pos[$i];
			my $ts_pos_Href = $seg_OB->{pos_feature}->{$ts_pos};
			my $ts_pos_rad = $ts_pos_Href->{rad_now};
			# check the limit region
			if( $ts_pos_rad < $min_rad ){
				$ts_pos_Href->{rad_now} = $min_rad;
				redo ARRANGE;
			}
			if( $ts_pos_rad > $max_rad ){
				$ts_pos_Href->{rad_now} = $max_rad;
				$start_i = $i - 1; # as rad_now is decreased, see whether it overlaps with last one (left one).
				redo ARRANGE;
			}
			# check the overlap with left pos (face from inner to outer layer).
			if( $i > 0 ){
				my $lf_pos = $sort_fet_pos[$i-1];
				my $lf_pos_Href = $seg_OB->{pos_feature}->{$lf_pos};
				my $lf_pos_rad = $lf_pos_Href->{rad_now};
				if( $ts_pos_rad - $lf_pos_rad < $icon_line_min_gaprad ){
					$lf_pos_Href->{rad_now} -= $icon_loc_adjust_rad_step;
					$ts_pos_Href->{rad_now} += $icon_loc_adjust_rad_step;
					$start_i = $i - 1;
					redo ARRANGE;
				}
			}
			# check the overlap with right pos (face from inner to outer layer).
			if( $i < $#sort_fet_pos ){
				my $rt_pos = $sort_fet_pos[$i+1];
				my $rt_pos_Href = $seg_OB->{pos_feature}->{$rt_pos};
				my $rt_pos_rad = $rt_pos_Href->{rad_now};
				if( $rt_pos_rad - $ts_pos_rad < $icon_line_min_gaprad ){
					$ts_pos_Href->{rad_now} -= $icon_loc_adjust_rad_step;
					$rt_pos_Href->{rad_now} += $icon_loc_adjust_rad_step;
					redo ARRANGE;
				}
			}
		}
	}
	stout_and_sterr "[INFO]\tSeg NO.$seg_OB->{NO} arranged feature icons in $arrange_time times.\n";
}

#--- draw feature icons
sub draw_features_icon{
	my $seg_OB = shift;
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};
	my $icon_EdgeLen = $V_Href->{icon_edge_len};
	my $icon_radial_gap = $V_Href->{icon_radial_gap};

	# draw feature icons for each position
	my @next_outer_radius;
	for my $pos (sort {$a<=>$b} keys %{$seg_OB->{pos_feature}}){
		my $pos_Href = $seg_OB->{pos_feature}->{$pos};
		my $pos_rad_ori = $pos_Href->{rad_ori};
		my $pos_rad_adj = $pos_Href->{rad_now};
		my $features = $pos_Href->{features};
		my $icon_count = sum( map {scalar(@{$_->[1]})} @$features );
		# draw the line
		my $pos_outer_radius = 
			$seg_OB->draw_icon_stick( rad_ori => $pos_rad_ori,
									  rad_adj => $pos_rad_adj,
									  icon_count => $icon_count,
									  colour => $pos_Href->{sam_OB}->{pathcolour}->{1}, # defaults as 1st haplo-path
									  draw_bool => $draw_bool
									);
		push @next_outer_radius, $pos_outer_radius;
		# draw icons
		my $inner_edge_radius = $seg_OB->{outer_radius} + sum(@{$V_Href->{icon_line_len}});
		for my $partner_features (@$features){
			my $partner = $partner_features->[0];
			my $feature_Aref = $partner_features->[1];
			for my $i ( 0 .. scalar(@$feature_Aref)-1 ){
				my $fet_id = $feature_Aref->[$i];
				my $fet_Href = $Feature_Basic_Href->{feature}->{$fet_id};
				my $fet_type = $fet_Href->{type};
				my $short_name = $fet_Href->{short_name};
				$seg_OB->draw_signle_feature_icon(  rad => $pos_rad_adj,
													feature_type => $fet_type,
													feature_text => $short_name,
													inner_edge_radius => $inner_edge_radius,
													partner => $partner
												 ) if( $draw_bool );
				$inner_edge_radius += $icon_EdgeLen + $icon_radial_gap;
			}
		}
	}

	# update this seg_OB's outer_radius
	$seg_OB->{outer_radius} = max( @next_outer_radius );
}

#--- draw the icon line
sub draw_icon_stick{
	my $seg_OB = shift;
	my %parm = @_;
	my $rad_ori = $parm{rad_ori};
	my $rad_adj = $parm{rad_adj};
	my $icon_count = $parm{icon_count};
	my $line_col = $parm{colour};
	my $draw_bool = $parm{draw_bool};

	# variants
	my $icon_EdgeLen = $V_Href->{icon_edge_len};
	my $icon_radial_gap = $V_Href->{icon_radial_gap};
	my $icon_line_width = $V_Href->{icon_line_width};
	my ($line_stlen, $line_inlen, $line_edlen) = @{$V_Href->{icon_line_len}};

	my $icon_occupy_len = $icon_count * ($icon_EdgeLen + $icon_radial_gap) - $icon_radial_gap;
	# prepare dot on line
	my ($x1, $y1) = get_coordinate_on_circle(cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad_ori, radius=>$seg_OB->{outer_radius});
	my @line_dot_xy = ($x1, $y1);
	if( $rad_ori == $rad_adj ){
		my $line_exd_len = $icon_occupy_len + sum(@{$V_Href->{icon_line_len}}) + $line_edlen;
		my ($x2, $y2) = get_coordinate_on_circle(cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad_ori, radius=>$seg_OB->{outer_radius}+$line_exd_len);
		push @line_dot_xy , $x2, $y2;
	}
	else{
		my $line_exd_len = $line_stlen;
		my ($x2, $y2) = get_coordinate_on_circle(cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad_ori, radius=>$seg_OB->{outer_radius}+$line_exd_len);
		push @line_dot_xy , $x2, $y2;
		$line_exd_len += $line_inlen;
		my ($x3, $y3) = get_coordinate_on_circle(cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad_adj, radius=>$seg_OB->{outer_radius}+$line_exd_len);
		push @line_dot_xy , $x3, $y3;
		$line_exd_len += $icon_occupy_len + 2 * $line_edlen;
		my ($x4, $y4) = get_coordinate_on_circle(cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad_adj, radius=>$seg_OB->{outer_radius}+$line_exd_len);
		push @line_dot_xy , $x4, $y4;
	}
	# draw line
	for my $j (1 .. scalar(@line_dot_xy)/2-1){
		$V_Href->{circos_OB}->line(
										x1 => $line_dot_xy[($j-1)*2],
										y1 => $line_dot_xy[($j-1)*2+1],
										x2 => $line_dot_xy[$j*2],
										y2 => $line_dot_xy[$j*2+1],
										stroke => $line_col,
										'stroke-width' => $icon_line_width,
										'stroke-linecap' => 'round',
										opacity => 1
									) if($draw_bool);
	}
	# return this pos stick outer radius
	return $seg_OB->{outer_radius} + $icon_occupy_len + sum(@{$V_Href->{icon_line_len}}) + $line_edlen;
}

#--- draw a feature icon at given location on segment
sub draw_signle_feature_icon{
	my $seg_OB = shift;
	my %parm = @_;
	my $rad = $parm{rad};
	my $feature_type = $parm{feature_type};
	my $feature_text = $parm{feature_text};
	my $inner_edge_radius = $parm{inner_edge_radius};
	my $partner = $parm{partner};

	# variants
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};
	my $icon_EdgeLen = $V_Href->{icon_edge_len};

	# prepare for icon
	my $icon_shape = $Feature_Basic_Href->{type}->{$feature_type}->{icon};
	my $fill_color = $Feature_Basic_Href->{type}->{$feature_type}->{colour};
	my $rotate_degree = $rad / $deg2rad;
	# draw the feature icon
	my ($x,$y);
	if($icon_shape eq 'square'){
		($x,$y) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad, radius=>$inner_edge_radius + $icon_EdgeLen * 0.5 );
	}
	elsif($icon_shape eq 'triangle'){
		($x,$y) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$rad, radius=>$inner_edge_radius + $icon_EdgeLen * sqrt(3) * 0.25 );
	}
	# partner defines the boundary
	my $boud_dasharray = $partner eq 'host' ? '1, 1' : ' ';
	# draw this feature icon
	$seg_OB->draw_a_feature_icon(
							x => $x,
							y => $y,
							icon_shape => $icon_shape,
							fill_color => $fill_color,
							feature_text => $feature_text,
							rotate_degree => $rotate_degree,
							boud_dasharray => $boud_dasharray
						);
}

#--- just draw a feature icon at the give position (x,y)
sub draw_a_feature_icon{
	my $seg_OB = shift;
	my %parm = @_;
	my $x = $parm{x};
	my $y = $parm{y};
	my $icon_shape = $parm{icon_shape};
	my $feature_text = $parm{feature_text};
	my $fill_color = $parm{fill_color};
	my $rotate_degree = $parm{rotate_degree};
	my $boud_dasharray = $parm{boud_dasharray};

	# variants
	my $icon_EdgeLen = $parm{icon_EdgeLen} || $V_Href->{icon_edge_len};
	my $icon_font_size = $parm{font_size} || $icon_EdgeLen;

	if($icon_shape eq 'square'){
		draw_a_parallelogram(
								svg_obj => $V_Href->{circos_OB},
								x => $x,
								y => $y,
								head_bottom_side_len => $icon_EdgeLen,
								fill_color => $fill_color,
								text_fill => $feature_text,
								rotate_degree => $rotate_degree,
								text_font_size => $icon_font_size + 1,
								font_size_auto_adjust => 0,
								boud_dasharray => $boud_dasharray
							);
	}
	elsif($icon_shape eq 'triangle'){
		draw_a_triangle(
							svg_obj => $V_Href->{circos_OB},
							x => $x,
							y => $y,
							bottom_side_len => $icon_EdgeLen + 2,
							fill_color => $fill_color,
							text_fill => $feature_text,
							rotate_degree => $rotate_degree,
							text_font_size => $icon_font_size,
							font_size_auto_adjust => 0,
							boud_dasharray => $boud_dasharray
						 );
	}
}

#--- load information of one transcript under a gene's name
sub load_one_transcript{
	my $seg_OB = shift;
	my %parm = @_;
	my $this_trans_Href = $parm{trans_Href};

	my $gene_name = $this_trans_Href->{gname};
	my $trans_name = $this_trans_Href->{tname};
	# load transcript information into seg_OB
	$seg_OB->{gene}->{$gene_name}->{$trans_name}->{$_} = $this_trans_Href->{$_} for keys %$this_trans_Href;
}

#--- draw genes
sub draw_genes{
	my $seg_OB = shift;
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	my $max_layer_radius = 0;
	my %genes_layer = (1=>[]);
	for my $gname (sort keys %{$seg_OB->{gene}}){
		stout_and_sterr "[INFO]\tdraw gene $gname.\n";
		my @tname = sort keys %{$seg_OB->{gene}->{$gname}};
		my $tname = $tname[0]; # select the first one
		my $t_Href = $seg_OB->{gene}->{$gname}->{$tname};
		# transcript information
		my $refseg = $t_Href->{refseg};
		my $strand = $t_Href->{strand};
		my $ex_len_Aref = $t_Href->{ex_len};
		my $ex_lfp_Aref = $t_Href->{ex_lfp};

		#-------------#
		# which layer #
		#-------------#
		# each exon
		my @exon;
		for my $i ( 0 .. scalar(@$ex_len_Aref)-1 ){
			my $lfp = $ex_lfp_Aref->[$i] + 1;
			my $rtp = $ex_lfp_Aref->[$i] + $ex_len_Aref->[$i];
			# check whether outside of this segments
			$lfp = max( $lfp, $seg_OB->{lf_pos} );
			$rtp = min( $rtp, $seg_OB->{rt_pos} );
			# pos -> rad
			my ($lfp_rad, $rtp_rad);
			if( $rtp < $seg_OB->{lf_pos} ){ # totally outside at upstream
				$rtp = $seg_OB->{lf_pos}; # for following judgment
				$lfp_rad = $seg_OB->{start_rad};
				$rtp_rad = $seg_OB->{start_rad};
			}
			elsif( $lfp > $seg_OB->{rt_pos} ){ # totally outside at downstream
				$lfp = $seg_OB->{rt_pos}; # for following judgment
				$lfp_rad = $seg_OB->{start_rad} + $seg_OB->{rad_size};
				$rtp_rad = $seg_OB->{start_rad} + $seg_OB->{rad_size};
			}
			else{
				$lfp_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $lfp );
				$rtp_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $rtp );
			}
			# for outside exons, only store one time
			if(    scalar(@exon) != 0
				&& (   ( $rtp_rad == $exon[-1]->{rtp_rad} && $rtp == $seg_OB->{lf_pos} )
						|| ( $lfp_rad == $exon[-1]->{lfp_rad} && $lfp == $seg_OB->{rt_pos} )
					 )
				){
				next;
			}
			# record
			push @exon, { lfp_rad=>$lfp_rad, rtp_rad=>$rtp_rad };
		}
		# get layer info
		my $layer_info_Href = $seg_OB->gene_at_which_layer( layer_Href => \%genes_layer, exon_Aref => \@exon, gene_name => $gname );
		my $layer_NO = $layer_info_Href->{layer_NO};
		#------------#
		# draw exons #
		#------------#
		# my $gene_col = $COLOR_DB{ (int(rand(100)) % scalar(grep /^\d+$/, keys %COLOR_DB) + 1) };
		# my $strd_col = ($gene_col =~ /$COLOR_DB{dark_col}/)?'white':'black';
		my $gene_col = ( $strand eq '-' ) ? $V_Href->{gene_col_minus} :
								 ( ( $strand eq '+' ) ? $V_Href->{gene_col_plus} : $V_Href->{gene_col_region});
		my $last_rt_rad;
		my $min_exon_rad = $V_Href->{min_exon_length} / $layer_info_Href->{inner_radius};
		for my $i ( 0 .. $#exon ){
			my $lfp_rad = $exon[$i]->{lfp_rad};
			my $rtp_rad = max( $exon[$i]->{rtp_rad}, $lfp_rad+$min_exon_rad );
			# draw exon
			unless(   $lfp_rad == ($seg_OB->{start_rad}+$seg_OB->{rad_size})
						 || $rtp_rad == $seg_OB->{start_rad}
			){ # not outside of segment
				draw_circle_seg(
									svg_obj => $V_Href->{circos_OB},
									cx => $V_Href->{cx},
									cy => $V_Href->{cy},
									start_rad => $lfp_rad,
									rad_size => $rtp_rad - $lfp_rad,
									inner_radius => $layer_info_Href->{inner_radius},
									outer_radius => $layer_info_Href->{inner_radius} + $V_Href->{gene_exon_width},
									seg_fill_color => $gene_col,
									seg_boundary_color => 'NONE',
									seg_boundary_width => 0,
									draw_bool => $draw_bool
								 );
				# update the layer region
				push @{$genes_layer{$layer_NO}}, [$lfp_rad, $rtp_rad];
			}
			# draw intron line
			if( defined($last_rt_rad) && $last_rt_rad < $lfp_rad ){
				draw_an_arc(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $last_rt_rad,
								rad_size => $lfp_rad - $last_rt_rad,
								radius => $layer_info_Href->{inner_radius} + $V_Href->{gene_exon_width} * 0.5,
								boundary_color => $gene_col,
								boundary_width => $V_Href->{gene_intron_width},
								draw_bool => $draw_bool
							 );
				# update the layer region
				push @{$genes_layer{$layer_NO}}, [$last_rt_rad, $lfp_rad];
			}
			# for intron
			$last_rt_rad = $rtp_rad;
		}
		#------------#
		# draw arrow #
		#------------#
		if($exon[0]->{lfp_rad} < $exon[-1]->{rtp_rad}){
			my $arrow_center_rad = ($exon[0]->{lfp_rad}+$exon[-1]->{rtp_rad}) / 2;
			my $arrow_layer_radius = $layer_info_Href->{inner_radius} + $V_Href->{gene_exon_width} / 2;
			my ($arrow_cx, $arrow_cy) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$arrow_center_rad, radius=>$arrow_layer_radius);
			my $arrow_rotate_rad = $arrow_center_rad + ($strand eq '+' ? 1 : -1) * ($PI / 2);
			my $arrow_rotate_degree = ($arrow_rotate_rad / $deg2rad * 100) / 100;
			# draw arrow
			draw_a_triangle(
								svg_obj => $V_Href->{circos_OB},
								x => $arrow_cx,
								y => $arrow_cy,
								bottom_side_len => $V_Href->{gene_arrow_bottom_len},
								fill_color => 'white',
								rotate_degree => $arrow_rotate_degree
							 );
		}
		#-----------#
		# gene name #
		#-----------#
		my $TG_gap_rad = $V_Href->{gene_name_gap} / $layer_info_Href->{inner_radius};
		my $text_anchor;
		my $anchor_rad;
		my $text_arc_radius = $layer_info_Href->{inner_radius};
		if( $layer_info_Href->{gname_anchor} == 0 ){
			$text_anchor = 'end';
			$anchor_rad = $exon[0]->{lfp_rad} - $TG_gap_rad;
			# update the layer region
			push @{$genes_layer{$layer_NO}}, [$anchor_rad-$layer_info_Href->{gname_arc_rad}, $exon[0]->{lfp_rad}];
		}
		elsif( $layer_info_Href->{gname_anchor} == -1 ){
			$text_anchor = 'start';
			$anchor_rad = $exon[-1]->{rtp_rad} + $TG_gap_rad;
			# update the layer region
			push @{$genes_layer{$layer_NO}}, [$exon[-1]->{rtp_rad}, $anchor_rad+$layer_info_Href->{gname_arc_rad}];
		}
		else{ # 1E10, draw gname at upper layer
			$text_anchor = 'start';
			$anchor_rad = $exon[0]->{lfp_rad};
			# to the next layer
			$text_arc_radius += $V_Href->{gene_radial_gap} + $V_Href->{gene_exon_width};
			# update the layer region
			push @{$genes_layer{$layer_NO+1}}, [$anchor_rad, $anchor_rad+$layer_info_Href->{gname_arc_rad}];
		}
		normalize_radian(radian_Sref=>\$anchor_rad);
		my $gname_height_limit = $V_Href->{gene_name_height_limit};
		my $text_to_center = ($anchor_rad >= $V_Href->{axis_top_rad}->[0] && $anchor_rad <= $V_Href->{axis_top_rad}->[1]) ? 'top':'toe';
		my $radius_shift_ratio = ($text_to_center eq 'toe') ? (0.45) : (0.3);
		show_text_on_arc(
							svg_obj => $V_Href->{circos_OB},
							cx => $V_Href->{cx},
							cy => $V_Href->{cy},
							radius => $text_arc_radius + $layer_info_Href->{gname_height} * $radius_shift_ratio,
							anchor_degree => $anchor_rad / $deg2rad,
							font_family => $V_Href->{gene_name_fontfam},
							text_anchor => $text_anchor,
							text_to_center => $text_to_center,
							text => $gname,
							text_col => 'black',
							height_adjust => 1,
							height_limit => $gname_height_limit,
							draw_bool => $draw_bool
						);
		#--- get maximum layer radius ---#
		$max_layer_radius = max( $max_layer_radius, $layer_info_Href->{inner_radius}, $text_arc_radius );
	}

	# update the outer radius
	if( $max_layer_radius != 0 ){ # draw some gene(s)
		$seg_OB->{outer_radius} = $max_layer_radius + $V_Href->{gene_radial_gap} + $V_Href->{gene_exon_width};
	}
}

#--- determine which layer where the given gene could be ploted
sub gene_at_which_layer{
	my $seg_OB = shift;
	my %parm = @_;
	my $layer_Href = $parm{layer_Href};
	my $exon_Aref = $parm{exon_Aref};
	my $gene_name = $parm{gene_name};

	# variants
	my $gname_height_limit = $V_Href->{gene_name_height_limit};
	my $layer_height = $V_Href->{gene_radial_gap} + $V_Href->{gene_exon_width};

	my $maxlayerNOtoTest = 1000;
	my $next_layer = 1;
	my $radius;
	my ($gname_arc_rad, $gname_height);
	for my $layer_NO ( 1 .. $maxlayerNOtoTest ){ # sort {$a<=>$b} keys %$layer_Href
		$radius = $seg_OB->{outer_radius} + $layer_NO * $layer_height - $V_Href->{gene_exon_width};
		$next_layer = $layer_NO + 1;
		# text size
		my $gname_size_Aref = 
			show_text_on_arc(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx}, cy => $V_Href->{cy}, radius => $radius,
								text => $gene_name, font_family => $V_Href->{gene_name_fontfam},
								height_adjust => 1, height_limit => $gname_height_limit,
								draw_bool => 0
								);
		$gname_height = $gname_size_Aref->[0];
		$gname_arc_rad = $gname_size_Aref->[1];
		# text-gene gap -> rad
		my $TG_gap_rad = $V_Href->{gene_name_gap} / $radius;
		# gene-gene gap -> rad
		my $GG_gap_rad = $V_Href->{gene_tangent_gap} / $radius;
		# test overlap
		## 0, draw gene name at left end; -1, draw gene name at right end; 1E10, draw gname at upper layer
		my @j = (0 , -1);
		## genename anchor reset
		pop   @j if exists $V_Href->{genename_anchor_Href}->{gHead}->{$gene_name};
		shift @j if exists $V_Href->{genename_anchor_Href}->{gTail}->{$gene_name};
		for my $j ( @j, 1E10 ){
			my $overlap_bool = 0;
			my $last_rt_rad;
			for my $i ( 0 .. scalar(@$exon_Aref)-1 ){
				my $lfp_rad = $exon_Aref->[$i]->{lfp_rad};
				my $rtp_rad = $exon_Aref->[$i]->{rtp_rad};
				# add text part
				if( $j == $i || $i == ($j+scalar(@$exon_Aref)) ){
					if( $j == 0 ){ $lfp_rad -= ($gname_arc_rad+$TG_gap_rad) } # left end
					else         { $rtp_rad += ($gname_arc_rad+$TG_gap_rad) } # right end
				}
				# check whether extend outside of segment
				if(    $lfp_rad < $seg_OB->{start_rad}
					|| $rtp_rad > ($seg_OB->{start_rad}+$seg_OB->{rad_size})
				){
					$overlap_bool = 1;
					last;
				}
				# this layer's all regions
				for my $reg_Aref (@{$layer_Href->{$layer_NO}}){
					my ($reg_lf, $reg_rt) = @$reg_Aref;
					if(    Get_Two_Seg_Olen($lfp_rad, $rtp_rad, $reg_lf-$GG_gap_rad, $reg_rt+$GG_gap_rad, 1E-6) # for exon
						|| (   defined($last_rt_rad) && $lfp_rad > $last_rt_rad # for intron
							&& Get_Two_Seg_Olen($last_rt_rad, $lfp_rad, $reg_lf-$GG_gap_rad, $reg_rt+$GG_gap_rad, 1E-6)
							 )
						){
						$overlap_bool = 1;
						last;
					}
				}
				# stop properly
				last if($overlap_bool);
				# update
				$last_rt_rad = $rtp_rad;
			}
			# specific for genename at upper layer
			if( $overlap_bool == 0 && $j == 1E10 ){
				my $gene_lfp_rad = $exon_Aref->[0]->{lfp_rad};
				my $gene_rtp_rad = $exon_Aref->[-1]->{rtp_rad};
				my $gene_radsize = $gene_rtp_rad - $gene_lfp_rad;
				if( $gene_radsize > $gname_arc_rad ){
					# next layer's all regions
					for my $reg_Aref (@{$layer_Href->{$next_layer}}){
						my ($reg_lf, $reg_rt) = @$reg_Aref;
						if( Get_Two_Seg_Olen($gene_lfp_rad, $gene_lfp_rad+$gname_arc_rad, $reg_lf-$GG_gap_rad, $reg_rt+$GG_gap_rad, 1E-6) ){
							$overlap_bool = 1;
							last;
						}
					}
				}
				else{
					$overlap_bool = 1;
				}
			}
			# is this layer ok?
			if( $overlap_bool == 0){ # yes, this layer , this gname-anchor, has no overlaps
				return { layer_NO=>$layer_NO, gname_anchor=>$j, inner_radius => $radius,
						 gname_arc_rad => $gname_arc_rad, gname_height => $gname_height };
			}
		}
	}

	# warn
	warn_and_exit "Cannot find proper layer (1-$maxlayerNOtoTest) for gene $gene_name on segment NO.$seg_OB->{NO}.\n";
}

#--- draw sample copy number band
sub draw_CN_band{

	my $seg_OB = shift;
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};

	# find its sample(s)
	my $draw_band_bool = 0;
	for my $sam_id (sort @{$seg_OB->{sample}}){
		my $sam_OB = $Sample_Href->{$sam_id};
		# check each region
		my $most_lf_pos = $seg_OB->{rt_pos};
		my $most_rt_pos = $seg_OB->{lf_pos};
		my $draw_this_sample_bool = 0;
		for my $reg_id (sort keys %{$sam_OB->{host}}){
			my $refseg = $sam_OB->{host}->{$reg_id}->{refseg};
			my $lf_pos = $sam_OB->{host}->{$reg_id}->{lfpos};
			my $rt_pos = $sam_OB->{host}->{$reg_id}->{rtpos};
			my $copy_number = $sam_OB->{host}->{$reg_id}->{copy_number};
			# on this segment?
			unless( $refseg eq $seg_OB->{ref_seg} && $lf_pos >= $seg_OB->{lf_pos} && $rt_pos <= $seg_OB->{rt_pos} ){
				next;
			}
			# record the minimum and maximum region boundary
			$most_lf_pos = min( $most_lf_pos, $lf_pos );
			$most_rt_pos = max( $most_rt_pos, $rt_pos );
			# here comes to draw, but firstly draw the background
			if( $draw_this_sample_bool == 0 ){
			}
			# draw this reg's copy number
			## reg CN color opacity
			my $fill_rgb = $seg_OB->get_graded_RGB_color(  value => $copy_number, log => $V_Href->{dna_cn_log},
															min => $V_Href->{CN_min}, max => $V_Href->{CN_max},
															orig_col_RGB => $V_Href->{CN_color_rgb}
														);
			## draw reg
			my $lf_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $lf_pos );
			my $rt_rad = $seg_OB->pos2rad( refseg => $refseg, pos => $rt_pos );
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $lf_rad,
								rad_size => $rt_rad - $lf_rad,
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{dna_CN_band_width},
								seg_fill_color => "RGB($fill_rgb)",
								seg_boundary_color => 'none',
								seg_boundary_width => 0,
								draw_bool => $draw_bool
							 );
			# yes, draw
			$draw_this_sample_bool = 1;
		}
		# does this sample has the non-covered region on this segment?
		if( $most_lf_pos != $seg_OB->{lf_pos} ){
			my $most_lf_rad = $seg_OB->pos2rad( refseg => $seg_OB->{ref_seg}, pos => $most_lf_pos );
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $seg_OB->{start_rad},
								rad_size => $most_lf_rad - $seg_OB->{start_rad},
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{dna_CN_band_width},
								seg_fill_color => $V_Href->{CN_bg_color},
								seg_boundary_color => 'none',
								seg_boundary_width => 0,
								seg_opacity => 1,
								draw_bool => $draw_bool
							 );
		}
		if( $most_rt_pos != $seg_OB->{rt_pos} ){
			my $most_rt_rad = $seg_OB->pos2rad( refseg => $seg_OB->{ref_seg}, pos => $most_rt_pos );
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $most_rt_rad,
								rad_size => $seg_OB->{rad_size} - ($most_rt_rad - $seg_OB->{start_rad}),
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{dna_CN_band_width},
								seg_fill_color => $V_Href->{CN_bg_color},
								seg_boundary_color => 'none',
								seg_boundary_width => 0,
								seg_opacity => 1,
								draw_bool => $draw_bool
							 );
		}
		# note sample
		if( $draw_this_sample_bool ){
			$draw_band_bool = 1;
			# sample color note bilaterally
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $_,
								rad_size => $V_Href->{spacerad} * $V_Href->{CN_EXP_sample_note_rad_ratio},
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{dna_CN_band_width},
								seg_fill_color => $sam_OB->{pathcolour}->{1},
								seg_boundary_color => 'black',
								seg_boundary_width => 0.3,
								seg_opacity => 1,
								draw_bool => $draw_bool
							 ) for ( $seg_OB->{start_rad} - $V_Href->{spacerad} * $V_Href->{CN_EXP_sample_note_rad_ratio},
										 $seg_OB->{start_rad} + $seg_OB->{rad_size}
									 );
			# sample arc seg
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $seg_OB->{start_rad},
								rad_size => $seg_OB->{rad_size},
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{dna_CN_band_width},
								seg_boundary_color => $sam_OB->{pathcolour}->{1},
								seg_boundary_width => 0.3,
								seg_opacity => 1,
								draw_bool => $draw_bool
							 );
			# update outer_radius
			$seg_OB->{outer_radius} += $V_Href->{dna_CN_band_width};
			stout_and_sterr "[INFO]\tShow Copy Number of sample $sam_id on semgent NO.$seg_OB->{NO}.\n";
		}
	}

	# make up segments
	if( $draw_band_bool ){
		$seg_OB->makeup_sector( draw_bool => $draw_bool );
	}
}

#--- get color opacity
sub get_graded_RGB_color{
	my $seg_OB = shift;
	my %parm = @_;
	my $value = $parm{value} || 0.1;
	my $min_value = $parm{min} || 0.1;
	my $max_value = $parm{max} || 0.1;
	my $log_value = $parm{log};
	my $orig_col_RGB = $parm{orig_col_RGB}; # e.g., '255,0,0'

	if( $log_value != 0 ){
		$_ = log($_)/log($log_value) for ($value, $min_value, $max_value);
	}

	my $ratio = ($value - $min_value) / (($max_value - $min_value) || 1);

	my $most_white = 255; # the lightest end
	my @return_rgb;
	for my $each_col_w (split /,/, $orig_col_RGB){
		push @return_rgb, min(int($each_col_w + (1.0 - $ratio) * ($most_white - $each_col_w)), $most_white);
	}

	return join(',', @return_rgb);
}

#--- make up segment's sector
sub makeup_sector{
	my $seg_OB = shift;
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $seg_boundary_color = $parm{seg_boundary_color} || 'black';
	my $seg_boundary_width = $parm{boundary_width} || $V_Href->{segment_boundary_width};

	draw_circle_seg(
						svg_obj => $V_Href->{circos_OB},
						cx => $V_Href->{cx},
						cy => $V_Href->{cy},
						start_rad => $seg_OB->{start_rad},
						rad_size => $seg_OB->{rad_size}, 
						inner_radius => $seg_OB->{inner_radius},
						outer_radius => $seg_OB->{outer_radius},
						seg_boundary_color => $seg_boundary_color,
						seg_boundary_width => $seg_boundary_width,
						# seg_side_bold=>(($seg_NO==0)?'0':'1'), # source side
						# seg_side_bold_opacity=>0.6,
						draw_bool => $draw_bool
					 );
}

#--- draw rna expression band
sub draw_EXP_band{
	my $seg_OB = shift;
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};

	# find its sample(s)
	my $win_rad = $V_Href->{EXP_win_arc_len} / $seg_OB->{outer_radius};
	my $draw_band_bool = 0;
	for my $sam_id (sort @{$seg_OB->{sample}}){
		my $sam_OB = $Sample_Href->{$sam_id};
		next unless(exists($seg_OB->{exp_rad_win}->{$sam_id}));
		my $win_depth_Href = $seg_OB->{exp_rad_win}->{$sam_id};
		# each rad window
		for my $rad_win (sort {$a<=>$b} keys %$win_depth_Href){
			my $mean_depth = $win_depth_Href->{$rad_win}->{mean_depth};
			## window expression color opacity
			my $fill_rgb = $seg_OB->get_graded_RGB_color(  value => $mean_depth, log => $V_Href->{exp_log},
															min => $V_Href->{EXP_min}, max => $V_Href->{EXP_max},
															orig_col_RGB => $V_Href->{EXP_color_rgb}
														);
				## draw rad window
			my $lf_rad = $seg_OB->{start_rad} + $win_rad * ($rad_win - 1);
			my $rt_rad = $seg_OB->{start_rad} + min( $win_rad * $rad_win, $seg_OB->{rad_size} );
			draw_circle_seg(
								svg_obj => $V_Href->{circos_OB},
								cx => $V_Href->{cx},
								cy => $V_Href->{cy},
								start_rad => $lf_rad,
								rad_size => $rt_rad - $lf_rad,
								inner_radius => $seg_OB->{outer_radius},
								outer_radius => $seg_OB->{outer_radius} + $V_Href->{rna_EXP_band_width},
								seg_fill_color => "RGB($fill_rgb)",
								seg_boundary_color => 'none',
								seg_boundary_width => 0,
								draw_bool => $draw_bool
							 );
		}
		## sample color note bilaterally
		draw_circle_seg(
							svg_obj => $V_Href->{circos_OB},
							cx => $V_Href->{cx},
							cy => $V_Href->{cy},
							start_rad => $_,
							rad_size => $V_Href->{spacerad} * $V_Href->{CN_EXP_sample_note_rad_ratio},
							inner_radius => $seg_OB->{outer_radius},
							outer_radius => $seg_OB->{outer_radius} + $V_Href->{rna_EXP_band_width},
							seg_fill_color => $sam_OB->{pathcolour}->{1},
							seg_boundary_color => 'black',
							seg_boundary_width => 0.3,
							seg_opacity => 1,
							draw_bool => $draw_bool
						 ) for ( $seg_OB->{start_rad} - $V_Href->{spacerad} * $V_Href->{CN_EXP_sample_note_rad_ratio},
								 $seg_OB->{start_rad} + $seg_OB->{rad_size}
								);
		# sample arc seg
		draw_circle_seg(
							svg_obj => $V_Href->{circos_OB},
							cx => $V_Href->{cx},
							cy => $V_Href->{cy},
							start_rad => $seg_OB->{start_rad},
							rad_size => $seg_OB->{rad_size},
							inner_radius => $seg_OB->{outer_radius},
							outer_radius => $seg_OB->{outer_radius} + $V_Href->{rna_EXP_band_width},
							seg_boundary_color => $sam_OB->{pathcolour}->{1},
							seg_boundary_width => 0.3,
							seg_opacity => 1,
							draw_bool => $draw_bool
						 );
		# update outer_radius
		$seg_OB->{outer_radius} += $V_Href->{rna_EXP_band_width};
		# record draw
		$draw_band_bool = 1;
		stout_and_sterr "[INFO]\tShow expression of sample $sam_id on semgent NO.$seg_OB->{NO}.\n";
	}

	# make up segments
	if( $draw_band_bool ){
		$seg_OB->makeup_sector( draw_bool => $draw_bool );
	}
}

#--- draw sample expression band
sub gather_exp_depth_rad_win{
	my $seg_OB = shift;

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};

	# find its sample(s)
	my $win_rad = $V_Href->{EXP_win_arc_len} / $seg_OB->{outer_radius};
	for my $sam_id (sort @{$seg_OB->{sample}}){
		my $sam_OB = $Sample_Href->{$sam_id};
		next unless($sam_OB->{rnaReads_Amount});
		# load the rna depth
		$seg_OB->{exp_rad_win}->{$sam_id} = {}; # {1}->{mean_depth=>, pos_depth=>{pos=>depth}}
		my $win_depth_Href = $seg_OB->{exp_rad_win}->{$sam_id};
		$seg_OB->get_given_sample_rna_depth(
											sam_id => $sam_id,
											win_depth_Href => $win_depth_Href,
											win_rad => $win_rad
										);
		# record minimum, maximum rna depth
		for my $rad_win (sort {$a<=>$b} keys %$win_depth_Href){
			# print "$sam_id\t$rad_win\t$win_depth_Href->{$rad_win}->{mean_depth}\n";
			$V_Href->{EXP_min} = min( $win_depth_Href->{$rad_win}->{mean_depth}, $V_Href->{EXP_min} );
			$V_Href->{EXP_max} = max( $win_depth_Href->{$rad_win}->{mean_depth}, $V_Href->{EXP_max} );
		}
		stout_and_sterr "[INFO]\tLoad RNA expression of sample $sam_id on semgent NO.$seg_OB->{NO}.\n";
	}
}

#--- read rna depth file for one sample at this segment
sub get_given_sample_rna_depth{
	my $seg_OB = shift;
	my %parm = @_;
	my $sam_id = $parm{sam_id};
	my $win_depth_Href = $parm{win_depth_Href};
	my $win_rad = $parm{win_rad};

	# variants
	my $Sample_RNA_Pos_Depth_file = $V_Href->{Sample_RNA_Pos_Depth_file};
	my $Sample_Href = $V_Href->{Sample_Href};
	my $Sample_Line_Index = $V_Href->{EXP_Depth_Sample_Index_Href};

	my $sam_OB = $Sample_Href->{$sam_id};
	warn_and_exit "Cannot find RNA reads amount of sample $sam_id.\n" unless($sam_OB->{rnaReads_Amount});
	my $norm_ratio = $V_Href->{RRA_max} / $sam_OB->{rnaReads_Amount};

	my @depth_buffer;
	my $last_rad_win = 1;
	my $sample_deal_bool = 0;
	open (RPD, Try_GZ_Read($Sample_RNA_Pos_Depth_file)) || die "fail read Sample_RNA_Pos_Depth_file: $!\n";
	(my $theme_line = lc(<RPD>)) =~ s/^#//;
	my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
	# quick read via index line number
	if( exists($Sample_Line_Index->{$sam_id}) ){
		# print "find $sam_id in Sample_Line_Index.\n";
		<RPD> for ( 2 .. $Sample_Line_Index->{$sam_id}-1 );
	}
	# pick data.
	while(<RPD>){
		next if(/^\#/);
		my @info = split;
		my %segif_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
		# for sample index
		unless( exists($Sample_Line_Index->{$segif_Ob{sampleid}}) ){
			$Sample_Line_Index->{$segif_Ob{sampleid}} = $.;
			# print "load $segif_Ob{sampleid} to Sample_Line_Index.\n";
		}
		# quick quit
		if( $segif_Ob{sampleid} ne $sam_id ){
			next if($sample_deal_bool == 0);
			last if($sample_deal_bool);
		}
		$sample_deal_bool = 1;
		my $refseg = $segif_Ob{refseg};
		my $pos = $segif_Ob{pos};
		my $rna_depth = $segif_Ob{rna_depth};
		next unless( $refseg eq $seg_OB->{ref_seg} && $pos >= $seg_OB->{lf_pos} && $pos <= $seg_OB->{rt_pos});
		push @depth_buffer, [ $pos, $rna_depth * $norm_ratio ];
		# distribute the pos depth to rad window
		if( $#depth_buffer > $V_Href->{EXP_depth_buffer} ){
			$seg_OB->deal_rna_depth_buffer(
											buffer_Aref => \@depth_buffer,
											win_depth_Href => $win_depth_Href,
											win_rad => $win_rad,
											last_rad_win_Sref => \$last_rad_win
										);
		}
	}
	close RPD;
	# deal with the remains
	$seg_OB->deal_rna_depth_buffer(
									buffer_Aref => \@depth_buffer,
									win_depth_Href => $win_depth_Href,
									win_rad => $win_rad,
									last_rad_win_Sref => \$last_rad_win
								);
	# deal last rad window
	$seg_OB->calculate_mean_depth_of_rad_win( rad_win => $last_rad_win, win_depth_Href => $win_depth_Href );
}

#--- extract depth buffer
sub deal_rna_depth_buffer{
	my $seg_OB = shift;
	my %parm = @_;
	my $buffer_Aref = $parm{buffer_Aref};
	my $win_depth_Href = $parm{win_depth_Href};
	my $win_rad = $parm{win_rad};
	my $last_rad_win_Sref = $parm{last_rad_win_Sref};

	while( scalar(@$buffer_Aref) != 0 ){ # till to buffer array empty
		my $st_rad_win = $seg_OB->pos2radwin( win_rad => $win_rad, pos => $buffer_Aref->[0]->[0] );
		# still last rad_win?
		if( $st_rad_win != $$last_rad_win_Sref ){
			# check sorted
			if( $st_rad_win < $$last_rad_win_Sref ){
				warn_and_exit "Position not sorted in RNA depth file.\n";
			}
			# deal passed rad windows in win_depth_Href Hash
			for my $rad_win ( $$last_rad_win_Sref .. $st_rad_win-1 ){
				$seg_OB->calculate_mean_depth_of_rad_win( rad_win => $rad_win, win_depth_Href => $win_depth_Href );
			}
			# update last_rad_win
			$$last_rad_win_Sref = $st_rad_win;
		}
		# transfer buffer to win_depth_Href Hash
		## find the continuous elements belong to st_rad_win
		my $i = scalar(@$buffer_Aref) - 1;
		my $i_rad_win = $seg_OB->pos2radwin( win_rad => $win_rad, pos => $buffer_Aref->[$i]->[0] );
		while( $i_rad_win != $st_rad_win ){
			$i = int($i * 0.5); # dichotomy
			last if( $i == 0 ); # quick
			$i_rad_win = $seg_OB->pos2radwin( win_rad => $win_rad, pos => $buffer_Aref->[$i]->[0] );
		}
		## now transfer
		$win_depth_Href->{$st_rad_win}->{pos_depth}->{$buffer_Aref->[$_]->[0]} = $buffer_Aref->[$_]->[1] for (0 .. $i);
		## discard these elements
		if( $i == scalar(@$buffer_Aref) - 1 ){ # quick
			@$buffer_Aref = ();
			last;
		}
		else{
			shift @$buffer_Aref for (0 .. $i);
		}
	}
}

#--- calculate mean depth of one rad window
sub calculate_mean_depth_of_rad_win{
	my $seg_OB = shift;
	my %parm = @_;
	my $win_depth_Href = $parm{win_depth_Href};
	my $rad_win = $parm{rad_win};

	# {1}->{mean_depth=>, pos_depth=>{pos=>depth}}
	# never record?
	if( !exists($win_depth_Href->{$rad_win}->{pos_depth}) ){
		$win_depth_Href->{$rad_win}->{mean_depth} = 0;
	}
	else{
		my $pos_depth_Href = $win_depth_Href->{$rad_win}->{pos_depth};
		my $sum_depth = sum(values %$pos_depth_Href);
		my @sorted_pos = sort {$a<=>$b} keys %$pos_depth_Href;
		$win_depth_Href->{$rad_win}->{mean_depth} = $sum_depth / ($sorted_pos[-1] - $sorted_pos[0] + 1);
		delete $win_depth_Href->{$rad_win}->{pos_depth};
	}
}

#--- get rad window NO of given pos
sub pos2radwin{
	my $seg_OB = shift;
	my %parm = @_;
	my $pos = $parm{pos};
	my $win_rad = $parm{win_rad};

	my $rad = $seg_OB->pos2rad( refseg => $seg_OB->{ref_seg}, pos => $pos );
	return int( ($rad - $seg_OB->{start_rad}) / $win_rad ) + 1;
}

1; ## tell the perl script the successful access of this module.
