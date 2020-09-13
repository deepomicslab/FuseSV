package FuseSV::Visual::AxisWork;

use strict;
use warnings;
use List::Util qw/max min/;
use SOAPfuse::General_Operation qw/stout_and_sterr/;
use SOAPfuse::SVG_Font qw/show_text_in_line/;
use SOAPfuse::SVG_Orthogonal_System_Elements qw/draw_a_parallelogram draw_a_triangle/;
use SVG;
use FuseSV::Visual::ArrangeObj qw/ get_pos_windowNO windowNO_to_posInterval test_overlap_at_one_layer /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  initialize_SVG_obj
			  show_Y_axis
			  set_Y_axis_Resol_LabStep
			  show_X_axis
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'AxisWork';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-04-12';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						initialize_SVG_obj
						show_Y_axis
						set_Y_axis_Resol_LabStep
						show_X_axis
						CalculateRegResolColorOpacity
					 /;

#--- initilize SVG object ---
sub initialize_SVG_obj{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $bg_width  = $parm{bg_width};
	my $bg_height = $parm{bg_height};
	my $bg_col = $parm{bg_col} || 'white';
	my $bg_stroke_col = $parm{bg_stroke_col} || 'none';
	my $AUTHOR = $parm{AUTHOR};
	my $EMAIL = $parm{EMAIL};

	$$SVG_obj_Oref = SVG->new( width=>$bg_width, height=>$bg_height, author=>$AUTHOR, 'author-mail'=>$EMAIL);
	$$SVG_obj_Oref->rect(x=>0, y=>0, width=>$bg_width, height=>$bg_height, fill=>$bg_col, stroke=>$bg_stroke_col);
}

#--- whole work to draw Y axis ---
#--- default: rightVertical; set rightVerti as 1 for rightVertical
sub show_Y_axis{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	# create SVG object if possible
	my $xAxisLen = $parm{xAxisLen} || 1000;
	my $bgStrokeCol = $parm{bgStrokeCol} || 'none';
	# attributes
	my $yResol = $parm{yResol};
	my $yMaxValue = $parm{yMaxValue};
	my $yAxisLen = $parm{yAxisLen};
	my $yLabStep = $parm{yLabStep};
	my $disMaxValueUpdate = $parm{disMaxValueUpdate} || 0;
	# location
	my $draw_baseY = $parm{draw_baseY};
	   $draw_baseY = 100 if( !defined $draw_baseY );
	my $axisZX = $parm{axisZX};
	my $axisZY_Sref = $parm{axisZY_Sref};
	my $bi_extWidth = $parm{bi_extWidth};
	   $bi_extWidth = 50 if( !defined $bi_extWidth );
	my $bi_extHeight = $parm{bi_extHeight};
	   $bi_extHeight = 50 if( !defined $bi_extHeight );
	# appearance
	my $yStrokeCol = $parm{yStrokeCol} || 'black';
	my $yStrokeWid = $parm{yStrokeWid} || 1;
	my $yLabelLlen = $parm{yLabelLlen} || 3;
	my $yLabelUnit = $parm{yLabelUnit};
	my $yLabelRodg = $parm{yLabelRodg} || 0;
	my $yLabelFtsz = $parm{yLabelFtsz} || 12;
	my $yLabelFtfm = $parm{yLabelFtfm} || 'Arial';
	my $yLegd_Href = $parm{yLegd_Href};
	my $yLegIcSize = $parm{yLegIcSize} || 15;
	# right Vertical Y axis ?
	my $rightVerti = $parm{rightVerti} || 0;
	my $rgtvtXshift = $rightVerti ? $xAxisLen : 0;
	my $rgtvtRatio = $rightVerti ? 1 : -1;

	# location of axis Zero point
	$$axisZY_Sref = $draw_baseY + $yAxisLen;

	# initialize SVG object
	if( !defined $$SVG_obj_Oref ){
		&initialize_SVG_obj(
							SVG_obj_Oref => $SVG_obj_Oref,
							bg_width => $axisZX + $xAxisLen + $bi_extWidth,
							bg_height => $$axisZY_Sref + $bi_extHeight,
							bg_stroke_col => $bgStrokeCol,
							AUTHOR => $AUTHOR,
							EMAIL => $EMAIL
						);
	}


	# draw y axis
	$$SVG_obj_Oref->line(
							x1 => $axisZX + $rgtvtXshift + 1 * $rgtvtRatio,
							y1 => $$axisZY_Sref - $yAxisLen,
							x2 => $axisZX + $rgtvtXshift + 1 * $rgtvtRatio,
							y2 => $$axisZY_Sref,
							stroke => $yStrokeCol,
							'stroke-width' => $yStrokeWid,
							'stroke-linecap' => 'round'
						);

	# y axis text
	my $yLabValue = 0;
	my $yLabelMaxWid = 0;
	while( $yLabValue <= $yMaxValue ){
		my $y_loc = $$axisZY_Sref - $yLabValue / $yResol;
		# short Label line
		$$SVG_obj_Oref->line(
								x1 => $axisZX + $rgtvtXshift + (1 + $yLabelLlen) * $rgtvtRatio,
								y1 => $y_loc,
								x2 => $axisZX + $rgtvtXshift + 1 * $rgtvtRatio,
								y2 => $y_loc,
								stroke => $yStrokeCol,
								'stroke-width' => $yStrokeWid,
								'stroke-linecap' => 'round'
							);
		# label text
		my $show_text = ( !defined $yLabelUnit || $yLabValue != $yMaxValue ) ? $yLabValue : "($yLabelUnit) $yLabValue";
		my $text_size_Aref = 
			show_text_in_line(
								$SVG_obj_Oref,
							  	text_x => $axisZX + $rgtvtXshift + (1 + $yLabelLlen + 3) * $rgtvtRatio,
								text_y => $y_loc + $yLabelFtsz/3,
								text => $show_text,
								font_family => $yLabelFtfm,
								font_size => $yLabelFtsz,
								text_col => $yStrokeCol,
								text_anchor => ($rightVerti ? 'start' : 'end'),
								rotate_degree => $yLabelRodg
							  );
		# update
		$yLabValue += $yLabStep;
		$yLabelMaxWid = max( $yLabelMaxWid, $text_size_Aref->[1] );
	}

	# show legend near Y axis
	if( defined $yLegd_Href){
		# show legend-icon near Y axis
		my $LegTextIconGap = 3;
		my $LegIconWidth = $yLegIcSize;
		my $LegIconHeight = $LegIconWidth;
		my $LegIconStrokeWid = 0.5;
		my $legYshift = $yLegd_Href->{loc} eq 'top' ? $yAxisLen : 0;
		my $legYratio = $yLegd_Href->{loc} eq 'top' ? 1 : -1;
		my $LegVerOuterEdgeY = $$axisZY_Sref - $legYshift;
		my $LegHorInnerEdgeX = $axisZX + $rgtvtXshift + ($yLabelMaxWid + $LegTextIconGap * 5) * $rgtvtRatio;
		my $LegIcon_X = $LegHorInnerEdgeX + $LegIconWidth / 2 * $rgtvtRatio;
		# show legend
		for my $info_Href ( @{$yLegd_Href->{info_Aref}} ){
			# show icon
			my $LegIcon_Y = $LegVerOuterEdgeY + $LegIconHeight / 2 * $legYratio;
			if( $info_Href->{shape} eq 'p' ){
				draw_a_parallelogram(
										$SVG_obj_Oref,
										x => $LegIcon_X,
										y => $LegIcon_Y,
										fill_color => $info_Href->{fill_col},
										head_bottom_side_len => $LegIconWidth,
										left_right_side_len => $LegIconHeight,
										boundary_width => $LegIconStrokeWid
									);
			}
			elsif( $info_Href->{shape} eq 't' ){
				draw_a_triangle(
									$SVG_obj_Oref,
									x => $LegIcon_X,
									y => $LegIcon_Y,
									fill_color => $info_Href->{fill_col},
									bottom_side_len => $LegIconWidth,
									boundary_width => $LegIconStrokeWid
								);
			}
			elsif( $info_Href->{shape} eq 'c' ){
				$$SVG_obj_Oref->circle(
										cx => $LegIcon_X,
										cy => $LegIcon_Y,
										r => $LegIconWidth / 2,
										stroke => 'black',
										'stroke-width'=> $LegIconStrokeWid,
										fill => $info_Href->{fill_col}
									  );
			}
			# show text
			my $LegText_X = $LegHorInnerEdgeX + ($LegIconWidth + $LegTextIconGap) * $rgtvtRatio;
			my $LegText_Y = $LegIcon_Y;
			show_text_in_line(
								$SVG_obj_Oref,
							  	text_x => $LegText_X,
								text_y => $LegText_Y,
								text => $info_Href->{text},
								font_family => $info_Href->{font_family} || $yLabelFtfm,
								font_size => $info_Href->{font_size} || $yLabelFtsz,
								text_col => $info_Href->{text_col} || 'black',
								text_anchor =>  $rightVerti ? 'start' : 'end',
								height_adjust => 1,
								height_limit => $LegIconHeight
							);
			# update
			$LegVerOuterEdgeY += ( $LegIconHeight + $LegTextIconGap ) * $legYratio;
		}
	}
}

#--- automatically get the resol and lab step of Y axis ---
sub set_Y_axis_Resol_LabStep{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $yResol_Sref = $parm{yResol_Sref};
	my $yMaxValue_Sref = $parm{yMaxValue_Sref};
	my $yAxisLen_Sref = $parm{yAxisLen_Sref};
	my $yLabStep_Sref = $parm{yLabStep_Sref};
	my $disMaxValueUpdate = $parm{disMaxValueUpdate} || 0;

	# when Y axis is un-assigned
	$$yAxisLen_Sref = 80 if($$yAxisLen_Sref <= 0);
	# when no value to show
	$$yMaxValue_Sref = 10 if( $$yMaxValue_Sref == 0 );

	# get proper lab step of Y axis when need
	if(    $$yResol_Sref == 0
		|| $$yLabStep_Sref == 0
	){
		$$yResol_Sref ||= $$yMaxValue_Sref / $$yAxisLen_Sref;
		my $cand_yLabStep = $$yResol_Sref * 20;
		my @cand_yLabStep = ( 1E0, 2E0, 3E0, 4E0, 5E0, 6E0, 8E0,
							  1E1, 2E1, 3E1, 4E1, 5E1, 6E1, 8E1,
							  1E2, 2E2, 3E2, 4E2, 5E2, 6E2, 8E2,
							  1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 8E3,
							  1E4, 2E4, 3E4, 4E4, 5E4, 6E4, 8E4  );
		my @cand_yLabStep_idx = grep $cand_yLabStep[$_] <= $cand_yLabStep && $cand_yLabStep[$_+1] >= $cand_yLabStep, (0 .. $#cand_yLabStep-1);
		$$yLabStep_Sref = ( @cand_yLabStep_idx == 0 ? ( $cand_yLabStep > 10 ? 8E4 : 1 ) : $cand_yLabStep[$cand_yLabStep_idx[0]] );
	}

	# update max value of Y axis based on Y lab step
	unless($disMaxValueUpdate){
		$$yMaxValue_Sref = POSIX::ceil( $$yMaxValue_Sref / $$yLabStep_Sref ) * $$yLabStep_Sref;
	}

	# update length of Y axis
	$$yAxisLen_Sref = POSIX::ceil( $$yMaxValue_Sref / $$yResol_Sref );

	# inform
	stout_and_sterr `date`."[INFO]:\tY axis resolution: $$yResol_Sref\n"
								 ."\tY axis label_step: $$yLabStep_Sref\n";
}

#--- whole work to draw X axis ---
sub show_X_axis{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	# attributes
	my $axisZX = $parm{axisZX};
	my $axisZY_Sref = $parm{axisZY_Sref};
	my $draw_baseY_Sref = $parm{draw_baseY_Sref};
	my $xAxisLen = $parm{xAxisLen} || 1000;
	my $xAxisWidth = $parm{xAxisWidth} || 1;
	my $xLabelSpanPix = $parm{xLabelSpanPix} || 200;
	my $xLabelNumBase = $parm{xLabelNumBase}; # optional
	   $xLabelNumBase = int($xLabelNumBase) if( defined $xLabelNumBase);
	   $xLabelNumBase = 100 if( defined $xLabelNumBase && $xLabelNumBase <= 0 );
	my $xLabelLineHeight = $parm{xLabelLineHeight} || ($xAxisWidth * 3);
	my $xLabelFtsz = $parm{xLabelFtsz} || 12;
	my $xLabelFtfm = $parm{xLabelFtfm} || 'Arial';
	# for region in different resolution
	my $xResol_Href = $parm{xResol_Href};
	my $regSegName = $parm{regSegName};
	my $xRegMinPos = $parm{xRegMinPos};
	my $xRegMaxPos = $parm{xRegMaxPos};
	my $xRegResolColor = $parm{xRegResolColor};
	my $xRegResolColorMinOpacity = $parm{xRegResolColorMinOpacity};
	   $xRegResolColorMinOpacity = 0.3 if( !defined $xRegResolColorMinOpacity );
	my $xRegResolWidth = $xAxisWidth * 3;
	my $xLabelText = $parm{xLabelText};
	   $xLabelText = $regSegName if( !defined $xLabelText );

	#--- do resolution-color for region?
	my $MinRegResol = 1E10;
	my $MaxRegResol = 0;
	my $colorReg_bool = defined( $xRegResolColor ) &&
		&CalculateRegResolColorOpacity(
										regResol_Href => $xResol_Href,
										regSegName => $regSegName,
										MinPos => $xRegMinPos,
										MaxPos => $xRegMaxPos,
										MinOpacity => $xRegResolColorMinOpacity,
										MinResol_Sref => \$MinRegResol,
										MaxResol_Sref => \$MaxRegResol
									  );
	#--- show resolution-color line
	# note: $$axisZY_Sref will update.
	if( $colorReg_bool ){
		#--- gradual color-opacity
		$$axisZY_Sref += $xRegResolWidth / 2 + 1;
		my $line_Y = $$axisZY_Sref;
		for my $reg_st (sort {$a<=>$b} keys %{$xResol_Href->{$regSegName}} ){
			my $reg_ed = $xResol_Href->{$regSegName}->{$reg_st}->{ed_pos};
			next if($reg_ed < $xRegMinPos);
			last if($reg_st > $xRegMaxPos);
			my $resol  = $xResol_Href->{$regSegName}->{$reg_st}->{resol};
			my $opacity= $xResol_Href->{$regSegName}->{$reg_st}->{color_opacity};
			# region pos => win NO
			my %pos2WinNO = ();
			$pos2WinNO{$_} =
				get_pos_windowNO(
									refseg => $regSegName,
									pos => $xResol_Href->{$regSegName}->{$reg_st}->{$_},
									regResol_Href => $xResol_Href,
									stpos_1stwin => $xRegMinPos
								)
				for qw/ st_pos ed_pos /;
			# draw color line
			my $width = $pos2WinNO{ed_pos} - $pos2WinNO{st_pos} + 1;
			draw_a_parallelogram(
									$SVG_obj_Oref,
									x => $axisZX + $pos2WinNO{st_pos} - 1 + $width / 2,
									y => $line_Y,
									head_bottom_side_len => $width,
									left_right_side_len => $xRegResolWidth,
									fill_color => $xRegResolColor,
									opacity => $opacity,
									boundary_color => 'none'
								);
		}
		#--- show legend of gradual color-opacity at right hand
		my $axisLegendGap = 40;
		my $colorRegLegHeight = $xRegResolWidth * 3;
		my $colorRegLegWidth  = $xRegResolWidth * 8;
		my $colorRegLeg_X = $axisZX + $xAxisLen + $axisLegendGap + $colorRegLegWidth / 2;
		my $colorRegLeg_Y = $line_Y;
		draw_a_parallelogram(
								$SVG_obj_Oref,
								x => $colorRegLeg_X,
								y => $colorRegLeg_Y,
								colorGradOrit => 'h',
								colorGradStColor => $xRegResolColor,
								colorGradStOpacity => $xRegResolColorMinOpacity,
								colorGradEdColor => $xRegResolColor,
								colorGradEdOpacity => 1,
								head_bottom_side_len => $colorRegLegWidth,
								left_right_side_len => $colorRegLegHeight,
								boundary_color => 'none'
							);
		show_text_in_line(
							$SVG_obj_Oref,
						  	text_x => $colorRegLeg_X,
							text_y => $colorRegLeg_Y - $colorRegLegHeight / 2 - 1,
							text => 'resolution',
							font_family => $xLabelFtfm,
							font_size => $xLabelFtsz,
							text_col => 'black',
							text_anchor => 'middle',
							height_adjust => 0,
							width_limit => $colorRegLegWidth + 10
						);
		#--- show legend resol label
		show_text_in_line(
							$SVG_obj_Oref,
						  	text_x => $colorRegLeg_X + $colorRegLegWidth / 2 * ($_ eq 'Min' ? -1 : 1),
							text_y => $colorRegLeg_Y + $colorRegLegHeight / 2,
							text => ($_ eq 'Min' ? $MinRegResol : $MaxRegResol),
							font_family => $xLabelFtfm,
							font_size => $xLabelFtsz,
							text_col => 'black',
							text_anchor => 'middle',
							height_adjust => 2.2,
							width_limit => $colorRegLegWidth / 3 * length($_ eq 'Min' ? $MinRegResol : $MaxRegResol) / 2
						) for qw/ Min Max /;
		#--- update
		$$axisZY_Sref += $xAxisWidth + $xRegResolWidth / 2;
	}

	#--- draw X axis basic line
	$$SVG_obj_Oref->line(
							x1 => $axisZX,
							y1 => $$axisZY_Sref,
							x2 => $axisZX + $xAxisLen,
							y2 => $$axisZY_Sref,
							stroke => 'black',
							'stroke-width' => $xAxisWidth,
							'stroke-linecap' => 'round'
						);
	#--- show label-number
	my $xLabelAxisLineGap = 3;
	my $xLabel_Y = $$axisZY_Sref + $xLabelLineHeight + $xLabelAxisLineGap;
	my @xLabelOccuReg = ();
	my @xLabelLoc = (1, $xAxisLen);
	for my $t ( 1 .. int($xAxisLen/$xLabelSpanPix) ){
		if( $xLabelSpanPix * $t < $xAxisLen ){
			push @xLabelLoc, $xLabelSpanPix * $t;
		}
	}
	for my $idx ( 0 .. $#xLabelLoc ){
		my $xLabelLoc = $xLabelLoc[$idx];
		# get the position of this x-axis loc
		my @xPosItval = windowNO_to_posInterval(
													regResol_Href => $xResol_Href,
													refseg => $regSegName,
													windowNO => $xLabelLoc,
													stpos_1stwin => $xRegMinPos
												);
		# bilateral ends and no NumberBase required
		my $xLabelNum = ( $idx == 1 ? $xRegMaxPos : $xPosItval[0] );
		# find the nearest number in NumberBase, if set
		if(    defined $xLabelNumBase
			&& $xLabelNum % $xLabelNumBase != 0
			&& $idx != 0
			&& $idx != 1
		){
			# nearest Number in NumberBase
			my $remainder = $xLabelNum % $xLabelNumBase;
			$xLabelNum = $xLabelNumBase * 
						 (  ($remainder >= $xLabelNumBase / 2)
						  ? POSIX::ceil( $xLabelNum / $xLabelNumBase )
						  :         int( $xLabelNum / $xLabelNumBase ) );
			# update x-asix location
			$xLabelLoc = get_pos_windowNO(
											refseg => $regSegName,
											pos => $xLabelNum,
											regResol_Href => $xResol_Href,
											stpos_1stwin => $xRegMinPos
										 );
		}
		# show the x-axis label-number
		my $xLabel_X = $axisZX + $xLabelLoc;
		my $draw_bool = ($idx == 0 && $idx == 1) ? 1 : 0;
		for my $run_time (0, 1){
			# short line
			if( $draw_bool ){
				$$SVG_obj_Oref->line(
										x1 => $xLabel_X,
										y1 => $$axisZY_Sref,
										x2 => $xLabel_X,
										y2 => $xLabel_Y - $xLabelAxisLineGap,
										stroke => 'black',
										'stroke-width' => $xAxisWidth,
										'stroke-linecap' => 'round'
									);
			}
			# position
			my $text_size_Aref = show_text_in_line(
													$SVG_obj_Oref,
													text_x => $xLabel_X,
													text_y => $xLabel_Y,
													text => $xLabelNum,
													font_family => $xLabelFtfm,
													font_size => $xLabelFtsz,
													text_col => 'black',
													text_anchor => 'middle',
													height_adjust => 2,
													height_limit => 0,
													width_limit => 0,
													draw_bool => $draw_bool
												);
			my ($xLabelHeight, $xLabelWidth) = @$text_size_Aref;
			my $xLabelOccuReg_Aref = [ $xLabel_X - $xLabelWidth / 2, $xLabel_X + $xLabelWidth / 2 ];
			# update $$draw_baseY_Sref for below displaying
			$$draw_baseY_Sref = max( $$draw_baseY_Sref, $xLabel_Y + $xLabelHeight );
			# update xLabelOccuReg, and whether to show inner-label
			if( $draw_bool ){
				push @xLabelOccuReg, $xLabelOccuReg_Aref;
				last;
			}
			else{
				my $ovp_bool = test_overlap_at_one_layer(
															Sbjct_Reg_Aref => \@xLabelOccuReg,
															Query_Reg_Aref => [ $xLabelOccuReg_Aref ],
															obj_gap => 10
														 );
				if( $ovp_bool ){ # donot show this Label
					last;
				}
				else{ # allow show this Label
					$draw_bool = 1;
				}
			}
		}
	}
	#--- show lable-text at left-hand
	@xLabelOccuReg = sort { $a->[0] <=> $b->[0] } @xLabelOccuReg;
	my $xLabelText_X = $xLabelOccuReg[0]->[0] - 10;
	my @xLabelText = split /\n+/, $xLabelText;
	my $xLabel_yShift = 0;
	for my $i ( 0 .. $#xLabelText ){
		my $text_size_Aref = show_text_in_line(
												$SVG_obj_Oref,
												text_x => $xLabelText_X,
												text_y => $xLabel_Y + $xLabel_yShift,
												text => $xLabelText[$i],
												font_family => "$xLabelFtfm Bold",
												font_size => $xLabelFtsz,
												text_col => 'black',
												text_anchor => 'end',
												height_adjust => 2,
												height_limit => 0,
												width_limit => 0
											);
		# update, move upward
		$xLabel_yShift -= $text_size_Aref->[0] + 2;
	}
}

#--- to show reg resol in different opacity
sub CalculateRegResolColorOpacity{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $regResol_Href = $parm{regResol_Href};
	my $regSegName = $parm{regSegName};
	my $MinPos = $parm{MinPos};
	my $MaxPos = $parm{MaxPos};
	my $ratio_digit = $parm{ratio_digit} || 2;
	my $MinOpacity = $parm{MinOpacity} || 0.3;
	my $MinResol_Sref = $parm{MinResol_Sref};
	my $MaxResol_Sref = $parm{MaxResol_Sref};

	# get max resolution
	for my $reg_st (sort {$a<=>$b} keys %{$regResol_Href->{$regSegName}} ){
		my $reg_ed = $regResol_Href->{$regSegName}->{$reg_st}->{ed_pos};
		my $resol  = $regResol_Href->{$regSegName}->{$reg_st}->{resol};
		next if($reg_ed < $MinPos);
		last if($reg_st > $MaxPos);
		# record
		$$MinResol_Sref = min( $$MinResol_Sref, $resol );
		$$MaxResol_Sref = max( $$MaxResol_Sref, $resol );
	}
	# no, skip reg resol color work
	return 0 if( $$MinResol_Sref == $$MaxResol_Sref );

	# assign color opacity based on resolution
	$regResol_Href->{$regSegName}->{$_}->{color_opacity} = 
		sprintf "%.${ratio_digit}f",
		$MinOpacity + ($regResol_Href->{$regSegName}->{$_}->{resol}-$$MinResol_Sref) / ($$MaxResol_Sref-$$MinResol_Sref+1) * (1 - $MinOpacity)
		for keys %{$regResol_Href->{$regSegName}};
	# yes, do reg resol color work
	return 1;
}

#--- 
1; ## tell the perl script the successful access of this module.
