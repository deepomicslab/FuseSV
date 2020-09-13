package FuseSV::Visual::RMSKwork;

use strict;
use warnings;
use List::Util qw/ max min sum /;
use SOAPfuse::General_Operation qw/ stout_and_sterr /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read /;
use FuseSV::Visual::RegionWork qw/ arrange_region_with_clip get_ObjSpanPixItval /;
use SOAPfuse::SVG_Radian_System_Elements qw/ %COLOR_DB draw_circle_seg /;
use SOAPfuse::SVG_Orthogonal_System_Elements qw/ draw_a_parallelogram /;
use SOAPfuse::SVG_Font qw/ show_text_in_line /;
use FuseSV::Visual::ArrangeObj qw/ get_pos_windowNO allocate_object_vertical /;
use FuseSV::Visual::AxisWork qw/ initialize_SVG_obj /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  load_rmskInfo_from_idxBgz
			  allocate_RMSKele_and_draw
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'RMSKwork';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-02-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						load_rmskInfo_from_idxBgz
						allocate_RMSKele_and_draw
						show_rmsk_ele
					 /;

# assume RmskInfoPool_Href has such structure:
# RmskInfoPool_Href = {}
# RmskInfoPool_Href -> $ele_NO -> 'type'
# RmskInfoPool_Href -> $ele_NO -> 'detail'
# RmskInfoPool_Href -> $ele_NO -> 'strand'
# RmskInfoPool_Href -> $ele_NO -> 'refseg'
# RmskInfoPool_Href -> $ele_NO -> 'st_pos'
# RmskInfoPool_Href -> $ele_NO -> 'ed_pos'
# RmskInfoPool_Href -> $ele_NO -> 'body_drawItval' -> [[a,b], [c,d]]
# RmskInfoPool_Href -> $ele_NO -> 'body_color'
# RmskInfoPool_Href -> $ele_NO -> 'show_layer'

#--- load rmsk info ---
sub load_rmskInfo_from_idxBgz{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $rmsk_idx_bgz = $parm{rmsk_idx_bgz};
	my $RmskInfoPool_Href = $parm{RmskInfoPool_Href};
	my $idx_readReg = $parm{idx_readReg};
	my $tabix = $parm{tabix};
	my $skipEleType_Href = $parm{skipEleType_Href};
	my $onlyEleType_Href = $parm{onlyEleType_Href};
	my $clip_min_Pos = $parm{clip_min_Pos}; # effective with $clip_max_Pos
	my $clip_max_Pos = $parm{clip_max_Pos};
	my $cirl_oLen = $parm{cirl_oLen}; # effective with $clip_max_Pos

	# read viral gene
	my $ele_NO = 0;
	my @openFH_para = ( $rmsk_idx_bgz );
	push @openFH_para, (idx_readReg=>$idx_readReg, tabix=>$tabix) if( defined $idx_readReg && defined $tabix );
	open (RMSK, Try_GZ_Read(@openFH_para)) || die "fail load trans_psl: $!\n";
	while(<RMSK>){
		my @info = split;
		my $refseg    = $info[5];
		my $st_pos    = $info[6];
		my $ed_pos    = $info[7];
		my $strand    = $info[9];
		my $eleDetail = $info[10];
		my $eleType   = $info[11];
		# gene/trans filtration
		next if( defined $skipEleType_Href && (exists $skipEleType_Href->{$eleType} || exists $skipEleType_Href->{$eleDetail}) );
		next if( defined $onlyEleType_Href && scalar(keys %$onlyEleType_Href) != 0 && !exists $onlyEleType_Href->{$eleType} && !exists $onlyEleType_Href->{$eleDetail} );
		# deal with circular extended part, if possible
		my @clip_opt = ( min_Pos => $clip_min_Pos, max_Pos => $clip_max_Pos, cirl_oLen => $cirl_oLen );
		my @eleRegion = defined $clip_max_Pos
						? arrange_region_with_clip( lftPos => $st_pos, rgtpos => $ed_pos, @clip_opt )
						: ( [ $st_pos , $ed_pos ] );
		# filter out-of-range
		next if( scalar(@eleRegion) == 0 );
		$ele_NO ++;
		# record info
		$RmskInfoPool_Href->{$ele_NO} = {
										type => $eleType,
										detail => $eleDetail,
										strand => $strand,
										refseg => $refseg,
										label_text => $eleType,
										eleBody_Aref => \@eleRegion
									};
	}
	close RMSK;
}

#--- allocate rmsk ele in vertical layer and draw
sub allocate_RMSKele_and_draw{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $RmskInfoPool_Href = $parm{RmskInfoPool_Href};
	my $pos2wNO_para_Aref = $parm{pos2wNO_para_Aref};
	my $axisZX = $parm{axisZX};
	my $draw_baseY_Sref = $parm{draw_baseY_Sref};
	  $$draw_baseY_Sref = 100 if( !defined $$draw_baseY_Sref );
	my $eleMinWidth = $parm{eleMinWidth} || 1;
	# create SVG object if possible
	my $xAxisLen = $parm{xAxisLen} || 1000;
	my $yAxisLen = $parm{yAxisLen} || 100;
	my $bi_extWidth = $parm{bi_extWidth};
	   $bi_extWidth = 50 if( !defined $bi_extWidth );
	my $bi_extHeight = $parm{bi_extHeight};
	   $bi_extHeight = 50 if( !defined $bi_extHeight );
	my $bgStrokeCol = $parm{bgStrokeCol} || 'none';
	# label
	my $labelFtfm = $parm{labelFtfm} || 'Arial';
	my $labelFtsz = $parm{labelFtsz} || 12;

	# return if no rmsk ele to show
	if( scalar(keys %$RmskInfoPool_Href) == 0 ){
		stout_and_sterr "<WARN>:\tNone rmsk elements to show!\n";
		return;
	}

	# basic para
	my $eleLayer_height = 5;
	my $ele_height = $eleLayer_height - 2;
	my $ele_strd_width = 0.5;
	my $ele_gap = 1;

	#--- calculate the span interval of each rmsk ele
	get_ObjSpanPixItval(
							ObjInfoPool_Href => $RmskInfoPool_Href,
							pos2xAxisPos_para_Aref => $pos2wNO_para_Aref,
							axisZX => $axisZX,
							bodyRegAref_key => 'eleBody_Aref',
							bodyPixItval_key => 'body_drawItval',
							bodyCol_key => 'body_color',
							bodyColRef_key => 'type'
						);

	#--- arrange the Y location of rmsk ele
	my $max_show_layerNO = 0;
	allocate_object_vertical(
								obj_Href => $RmskInfoPool_Href,
								obj_itval_key => 'body_drawItval',
								obj_layer_key => 'show_layer',
								obj_gap => $ele_gap,
								max_layer_Sref => \$max_show_layerNO
							);
	# update the baseY to give space of rmsk
	my $rmsk_space_Y_size = $max_show_layerNO * $eleLayer_height;
	$$draw_baseY_Sref = max( $bi_extHeight, $rmsk_space_Y_size + 20 );

	#--- initialize SVG object
	if( !defined $$SVG_obj_Oref ){
		initialize_SVG_obj(
							SVG_obj_Oref => $SVG_obj_Oref,
							bg_width => $axisZX + $xAxisLen + $bi_extWidth,
							bg_height => $$draw_baseY_Sref + $yAxisLen + $bi_extHeight,
							bg_stroke_col => $bgStrokeCol,
							AUTHOR => $AUTHOR,
							EMAIL => $EMAIL
						);
	}

	my $ele_st_Y = $$draw_baseY_Sref;
	# update the base Y for y-axis
	$$draw_baseY_Sref += 3;
	#--- draw rmsk ele
	&show_rmsk_ele(
					SVG_obj_Oref => $SVG_obj_Oref,
					RmskInfoPool_Href => $RmskInfoPool_Href,
					pos2wNO_para_Aref => $pos2wNO_para_Aref,
					axisZX => $axisZX,
					ele_st_Y => $ele_st_Y,
					eleLayer_height => $eleLayer_height,
					ele_height => $ele_height,
					eleStrdWidth => $ele_strd_width,
					eleMinWidth => $eleMinWidth
				 );
	#--- show 'RMSK' line+label
	$$SVG_obj_Oref->line(
							x1 => $axisZX - 1,
							y1 => $ele_st_Y,
							x2 => $axisZX - 1,
							y2 => $ele_st_Y - $rmsk_space_Y_size,
							stroke => 'black',
							"stroke-width" => 1,
							'stroke-linecap' => 'round'
						);
	my @opt = $max_show_layerNO <= 3
			  ? ('end',    1, $max_show_layerNO * $eleLayer_height - 2, 0, 0  )
			  : ('middle', 0, 0, $max_show_layerNO * $eleLayer_height - 2, -90);
	show_text_in_line(
						$SVG_obj_Oref,
						text_x => $axisZX - 5,
						text_y => $ele_st_Y - $rmsk_space_Y_size / 2,
						text => 'RMSK',
						font_family => $labelFtfm,
						font_size => $labelFtsz,
						text_col => 'black',
						text_anchor => $opt[0],
						height_adjust => $opt[1],
						height_limit => $opt[2],
						width_limit => $opt[3],
						rotate_degree => $opt[4]
					 );
	#--- show legend
	my $rmskLabelGap = $eleLayer_height;
	my $rmskLabelSize = $ele_height + 2;
	my $rmskLabelSpace = max($rmskLabelGap * 10, 30);
	my $rmskLabel_X = $axisZX + $xAxisLen + $rmskLabelSpace;
	my $rmskLabel_Y = $ele_st_Y - $rmsk_space_Y_size;
	# used colors
	my %rmskType2Color = map { ( $RmskInfoPool_Href->{$_}->{type},
								 $RmskInfoPool_Href->{$_}->{body_color} )
							 } keys %$RmskInfoPool_Href;
	for my $rmskType (sort keys %rmskType2Color){
		my $rmskBodyCol = $rmskType2Color{$rmskType};
		# text
		show_text_in_line(
							$SVG_obj_Oref,
							text_x => $rmskLabel_X,
							text_y => $rmskLabel_Y,
							text => $rmskType,
							font_family => $labelFtfm,
							font_size => $labelFtsz,
							text_col => 'black',
							text_anchor => 'end',
							height_adjust => 1,
							height_limit => $rmskLabelSize,
							width_limit => 0
						 );
		# two strands
		my %strd_opt = ( '+'=>1, '-'=>2 );
		for my $strand ( '+', '-' ){
			# minus strand ?
			my @strd_para;
			push @strd_para, (
								inner_line_ori => 'h',
								inner_line_wid => $ele_strd_width,
								inner_line_col => ($rmskBodyCol =~ /$COLOR_DB{dark_col}/ ? 'white' : 'black')
							 ) if( $strand eq '-' );
			# draw rmsk element
			draw_a_parallelogram(
									$SVG_obj_Oref,
									x => $rmskLabel_X - $rmskLabelSize / 2 + ($rmskLabelGap + $rmskLabelSize) * $strd_opt{$strand},
									y => $rmskLabel_Y,
									fill_color => $rmskBodyCol,
									head_bottom_side_len => $rmskLabelSize,
									left_right_side_len => $rmskLabelSize,
									boundary_color => 'none',
									@strd_para
								);
		}
		# update
		$rmskLabel_Y += $rmskLabelSize + 1;
	}
}

#--- draw rmsk elements
sub show_rmsk_ele{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $RmskInfoPool_Href = $parm{RmskInfoPool_Href};
	my $pos2wNO_para_Aref = $parm{pos2wNO_para_Aref};
	my $axisZX = $parm{axisZX};
	   $axisZX = 100 if( !defined $axisZX );
	my $ele_st_Y = $parm{ele_st_Y} || 50;
	my $eleLayer_height = $parm{eleLayer_height} || 7;
	my $ele_height = $parm{ele_height} || 5;
	my $eleStrdWidth = $parm{eleStrdWidth} || ($ele_height / 5);
	my $eleMinWidth = $parm{eleMinWidth} || 1;

	for my $ele_NO ( sort keys %$RmskInfoPool_Href ){
		my $rmskInfo_Href = $RmskInfoPool_Href->{$ele_NO};
		my $ele_y = $ele_st_Y - $eleLayer_height * ( $rmskInfo_Href->{show_layer} - 0.5 );
		#--- intron line, gene name
		for my $reg_Aref ( @{$rmskInfo_Href->{body_drawItval}} ){
			my $body_lftPix = $reg_Aref->[0];
			my $body_rgtPix = $reg_Aref->[1];
			# minus strand ?
			my @strd_para;
			push @strd_para, (
								inner_line_ori => 'h',
								inner_line_wid => $eleStrdWidth,
								inner_line_col => ($rmskInfo_Href->{body_color} =~ /$COLOR_DB{dark_col}/ ? 'white' : 'black')
							 ) if( $rmskInfo_Href->{strand} eq '-' );
			# draw rmsk element
			draw_a_parallelogram(
									$SVG_obj_Oref,
									x => ($body_lftPix+$body_rgtPix)/2,
									y => $ele_y,
									fill_color => $rmskInfo_Href->{body_color},
									head_bottom_side_len => max( abs($body_rgtPix-$body_lftPix)+1, $eleMinWidth ),
									left_right_side_len => $ele_height,
									boundary_color => 'none',
									@strd_para
								);
		}
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
