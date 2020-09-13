package FuseSV::Visual::DepthWork;

use strict;
use warnings;
use List::Util qw/ max min sum /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Dist::DistStat qw/ get_value_mean /;
use FuseSV::Virus_Integ::DrawSegCN::ObjCountTransfMixToPureT qw/ get_ObjCountOfPureTumorCell /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  get_windowSmoDepth
			  deal_window_depth_info
			  get_ctrl_copyR
			  draw_depth_spectrum
			  get_givenWinItvalMeanDepth
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'DepthWork';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-04-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						get_windowSmoDepth
						deal_window_depth_info
						get_ctrl_copyR
						draw_depth_spectrum
						get_givenWinItvalMeanDepth
					 /;

# assume winDepth_Href has such structure:
# winDepth_Href = {}
# winDepth_Href -> $tissue -> $window_NO -> 'pos2depth' -> pos = orig_depth # soon delete
# winDepth_Href -> $tissue -> $window_NO -> 'pos2adjdepth' -> pos = adjust_depth
# winDepth_Href -> $tissue -> $window_NO -> 'mean_adjdepth' = adjust_depth_mean
# winDepth_Href -> $tissue -> $window_NO -> 'showdepth' = smooth_depth_to_show

#--- get depth to show of each window ---
#--- note: may smooth depth with bilateral flanking window
sub get_windowSmoDepth{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $winDepth_Href = $parm{winDepth_Href};
	my $tissue = $parm{tissue};
	my $DepSmooth_flkWd = $parm{DepSmooth_flkWd};
	my $smoDepthKey = $parm{smoDepthKey} || 'showdepth';
	my $bkWindow_Href = $parm{bkWindow_Href};
	my $maxDepth_Sref = $parm{maxDepth_Sref};
	my $maxWinNO_Sref = $parm{maxWinNO_Sref};
	# default
	$DepSmooth_flkWd = 2 unless( defined $DepSmooth_flkWd );

	my $tissueWinDepth_Href = $winDepth_Href->{$tissue};
	for my $window_NO (sort {$a<=>$b} keys %$tissueWinDepth_Href){
		if( $DepSmooth_flkWd != 0 ){
			my %DepthForEven;
			my %BK_sideBool;
			for my $i ( -1 * $DepSmooth_flkWd .. $DepSmooth_flkWd ){
				my $this_window_NO = $window_NO + $i;
				if( exists $tissueWinDepth_Href->{$this_window_NO} ){
					$DepthForEven{$i} = $tissueWinDepth_Href->{$this_window_NO}->{mean_adjdepth};
				}
				# has bread-point
				if(    defined $bkWindow_Href
					&& exists $bkWindow_Href->{$this_window_NO}
				){
					my $side = ( $i == 0 ? 0 : abs($i)/$i );
					$BK_sideBool{$side} = 1;
				}
			}
			# specific for break-point
			if( exists $BK_sideBool{0} ){
				%DepthForEven = ( 0 => $DepthForEven{0} );
			}
			else{
				if( exists $BK_sideBool{-1} ){ delete $DepthForEven{$_} for (-1 * $DepSmooth_flkWd .. -1) }
				if( exists $BK_sideBool{1}  ){ delete $DepthForEven{$_} for (1 .. $DepSmooth_flkWd) }				
			}
			# even depth
			$tissueWinDepth_Href->{$window_NO}->{$smoDepthKey} = sprintf "%.2f", sum(values %DepthForEven) / scalar(keys %DepthForEven);
		}
		else{
			$tissueWinDepth_Href->{$window_NO}->{$smoDepthKey} = $tissueWinDepth_Href->{$window_NO}->{mean_adjdepth};
		}
		# for further svg display
		if( defined $maxDepth_Sref ){
			$$maxDepth_Sref = max( $$maxDepth_Sref, $tissueWinDepth_Href->{$window_NO}->{$smoDepthKey} );
		}
		if( defined $maxWinNO_Sref ){
			$$maxWinNO_Sref = max( $$maxWinNO_Sref, $window_NO );
		}
	}
}

#--- deal depth info of given window ---
#--- apllied in DrawSegCN.pm and DrawViralGenome.pm ---
sub deal_window_depth_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $winDepth_Href = $parm{winDepth_Href}; # $V_Href->{dG_winDepth}
	my $window_NO = $parm{window_NO};
	my $tissue = $parm{tissue};
	my $role = $parm{role} || 'host';
	# for already-selected OK-pos
	my $selected_okPos_Href = $parm{selected_okPos_Href}; # $V_Href->{CtrlPos2Depth}
	# how to get mean value
	my $Eng_3SD = $parm{Eng_3SD} || 0; # $V_Href->{dG_DepEng_3SD}
	my $TmmRatio = $parm{TmmRatio} || 0; # $V_Href->{dG_DepTrimRatio}
	# has-T mode
	my $HasT_mode = $parm{HasT_mode} || 0; # $V_Href->{dG_HasT_mode}
	# has-T mode, for host role
	my $h_DepthToDraw = $parm{h_DepthToDraw} || 50; # $V_Href->{DepthToDraw}
	my $h_meanDepth_Href = $parm{h_meanDepth_Href}; # $V_Href->{dHG_meandepth}
	# has-T mode, for host role, and for ctrl
	my $h_InCase_Tcell_GMratio = $parm{h_InCase_Tcell_GMratio}; # $V_Href->{case_Tcell_GMratio}
	# has-T mode, for host role, and for case
	my $h_ctrl_prior_copyR_Aref = $parm{h_ctrl_prior_copyR_Aref}; # $V_Href->{dHG_ctrl_prior_copyR}
	my $h_InCase_Ncell_GMratio = $parm{h_InCase_Ncell_GMratio}; # $V_Href->{case_Ncell_GMratio}
	# has-T mode, for virus role, and for case
	my $v_case2norm_seqRatio = $parm{v_case2norm_seqRatio}; # $V_Href->{dVG_case2norm_seqRatio}

	# variants
	my $this_winDepth_Href = $winDepth_Href->{$tissue}->{$window_NO};

	my ($ctrl_copyratio, $ctrl_copyR_last_edp) = (1, 0);
	my @pos_depth;
	for my $pos (sort {$a<=>$b} keys %{$this_winDepth_Href->{pos2depth}}){
		# is OK-pos if set?
		if(    defined $selected_okPos_Href
			&& !exists $selected_okPos_Href->{$pos}
		){
			next;
		}
		my $depth = $this_winDepth_Href->{pos2depth}->{$pos};
		# specific work for has-T mode
		if( $HasT_mode && $depth != 0 ){
			if( $role eq 'host' ){ # such as human genome
				if( $tissue eq 'ctrl' ){
					# normalize to draw_depth (case_meandepth)
					$depth *= $h_DepthToDraw / $h_meanDepth_Href->{$tissue};
					# shrink to level of pure-case-cell
					$depth *= $h_InCase_Tcell_GMratio;
				}
				else{
					# update the germline CNV
					if( $pos > $ctrl_copyR_last_edp ){
						($ctrl_copyratio, $ctrl_copyR_last_edp) = &get_ctrl_copyR( ctrl_prior_copyR_Aref=>$h_ctrl_prior_copyR_Aref, pos=>$pos );
					}
					# subtracts depth belongs to the mixed normal-cells.
					$depth = get_ObjCountOfPureTumorCell(
															ObjCount => $depth,
															MixedCellMeanDepth => $h_meanDepth_Href->{$tissue},
															Ncell_GMratio => $h_InCase_Ncell_GMratio,
															ObjCopyRatioInNcell => $ctrl_copyratio
														);
				}
				# reset
				$depth = 0 if($depth < 0);
			}
			else{ # role == 'virus', for DrawViralGenome.pm
				if( $tissue eq 'case' ){
					# normalize to case sequence amount level
					$depth *= $v_case2norm_seqRatio;
				}
			}
		}
		# for short
		$depth = sprintf "%.3f", $depth;
		# restore the pos depth
		$this_winDepth_Href->{pos2adjdepth}->{$pos} = $depth;
		# for mean depth
		push @pos_depth, $depth;
	}
	# sweep
	delete $this_winDepth_Href->{pos2depth};
	# get mean adjdepth
	if(    scalar(@pos_depth) == 0
		||    sum(@pos_depth) == 0
	){
		$this_winDepth_Href->{mean_adjdepth} = 0;
	}
	else{
		$this_winDepth_Href->{mean_adjdepth} = 
			get_value_mean(
							value_Aref => \@pos_depth,
							Eng_SD_Ntimes => ( $Eng_3SD ? 3 : 0 ),
							TrimRatio => $TmmRatio
			);
	}
}

#--- get copy ratio of control sample at given position ---
sub get_ctrl_copyR{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $pos = $parm{pos};
	my $ctrl_prior_copyR_Aref = $parm{ctrl_prior_copyR_Aref};

	for my $reg_Href (@$ctrl_prior_copyR_Aref){
		if( $pos >= $reg_Href->{st_pos} && $pos <= $reg_Href->{ed_pos} ){
			return ($reg_Href->{copyR}, $reg_Href->{ed_pos});
		}
	}

	return (1, 0);
}

#--- draw depth spectrum of each tissue ---
sub draw_depth_spectrum{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $winDepth_Href = $parm{winDepth_Href}; # next key is tissue
	my $smoDepthKey = $parm{smoDepthKey} || 'showdepth';
	my $maxWinNOtoDraw = $parm{maxWinNOtoDraw};
	my $axisZX = $parm{axisZX};
	my $axisZY = $parm{axisZY};
	my $yResol = $parm{yResol};
	my $SpmCol_Href = $parm{SpmCol_Href}; # next key is tissue
	my $maxDraw1st = $parm{maxDraw1st} || 0;

	# draw spectrum with large one at the backend
	for my $window_NO ( 1 .. $maxWinNOtoDraw ){
		# first to draw larger one
		my @sort_tis =  $maxDraw1st
					   # max comes first
					   ? ( sort { $winDepth_Href->{$b}->{$window_NO}->{$smoDepthKey}
								  <=>
								  $winDepth_Href->{$a}->{$window_NO}->{$smoDepthKey}
								} keys %$winDepth_Href )
					   # case -> ctrl
					   : ( sort keys %$winDepth_Href );
		# draw spectrum
		for my $tis_idx ( 0 .. $#sort_tis ){
			my $tissue = $sort_tis[$tis_idx];
			my $opacity = ($tis_idx == 0) ? 1 : 0.8;
			my $smoDepth = $winDepth_Href->{$tissue}->{$window_NO}->{$smoDepthKey};
			next if( !defined $smoDepth || $smoDepth == 0 );
			$$SVG_obj_Oref->line(
									x1 => $axisZX + $window_NO,
									y1 => $axisZY,
									x2 => $axisZX + $window_NO,
									y2 => $axisZY - $smoDepth / $yResol,
									stroke => $SpmCol_Href->{$tissue},
									'stroke-width' => 1,
									'stroke-linecap' => 'round',
									opacity => $opacity
								);
		}
	}
}

#--- get mean depth of given window interval ---
sub get_givenWinItvalMeanDepth{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $winDepth_Href = $parm{winDepth_Href};
	my $tissue = $parm{tissue};
	my $depth_key = $parm{depth_key};
	my $st_winNO = $parm{st_winNO};
	my $ed_winNO = $parm{ed_winNO};

	my $tissueWinDepth_Href = $winDepth_Href->{$tissue};
	my @winDepth_pool = map $tissueWinDepth_Href->{$_}->{$depth_key},
						grep   exists $tissueWinDepth_Href->{$_}
							&& exists $tissueWinDepth_Href->{$_}->{$depth_key},
						( $st_winNO .. $ed_winNO );

	if( scalar(@winDepth_pool) == 0 ){
		warn_and_exit `date`."<ERROR>\tNo available window found in window interval [$st_winNO, $ed_winNO].\n";
	}
	else{
		return 
			get_value_mean(
							value_Aref => \@winDepth_pool,
							Eng_SD_Ntimes => 3
			);
	}
}

#---
1; ## tell the perl script the successful access of this module.
