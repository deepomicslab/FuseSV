package FuseSV::Visual::ArrangeObj;

use strict;
use warnings;
use POSIX qw/ceil/;
use List::Util qw/max min sum/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr Get_Two_Seg_Olen/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  load_regResolForWindowlize
			  get_pos_windowNO
			  windowNO_to_posInterval
			  allocate_object_bilater
			  allocate_object_vertical
			  test_overlap_at_one_layer
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ArrangeObj';
#----- version --------
$VERSION = "0.08";
$DATE = '2018-04-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						load_regResolForWindowlize
						get_pos_windowNO
						windowNO_to_posInterval
						allocate_object_bilater
						allocate_object_vertical
						test_overlap_at_one_layer
					 /;

# assume regResol_Href has such structure:
# regResol_Href = {}
# regResol_Href -> $refseg -> $st_pos -> 'st_pos' = $st_pos
# regResol_Href -> $refseg -> $st_pos -> 'ed_pos' = $ed_pos
# regResol_Href -> $refseg -> $st_pos -> 'resol'  = $resol

#--- load region resolution for further windowlize ---
sub load_regResolForWindowlize{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $reg_resol_file = $parm{reg_resol_file};
	my $regResol_Href = $parm{regResol_Href}; # to store resolution of region
	# destReg_Href -> {}
	# destReg_Href -> refseg = $refseg
	# destReg_Href -> st_pos = $st_pos
	# destReg_Href -> ed_pos = $ed_pos
	my $destReg_Href = $parm{destReg_Href}; # for check completeness
	my $destRegAutoFigureSize = $parm{destRegAutoFigureSize}; # || (defined $destReg_Href ? 500 : 0);
	my $autoAdjustResol = $parm{autoAdjustResol} || $destRegAutoFigureSize;

	my %drawReg_edgeShowBool = (st_pos=>0, ed_pos=>0);

	open (RES, Try_GZ_Read($reg_resol_file)) || die "fail read resolution file: $!\n";
	while (<RES>){
		next if(/^#/);
		my ($reg_seg, $reg_st, $reg_ed, $resolution) = (split)[0..3];
		# required segment of region
		if(    defined $destReg_Href
			&& $reg_seg ne $destReg_Href->{refseg}
		){
			next;
		}
		# simple check
		if(    $reg_st > $reg_ed
			|| $resolution <= 0
		){
			warn_and_exit "<ERRO>\tWrong data in resolution file.\n"
								 ."\t$_"
								 ."\treg_resol_file\n";
		}
		# clip region with destReg
		$reg_st = max( $reg_st, $destReg_Href->{st_pos} );
		$reg_ed = min( $reg_ed, $destReg_Href->{ed_pos} );
		next if($reg_st >= $reg_ed);
		# adjust resolution
		my $orig_resolution = $resolution;
		if( $autoAdjustResol ){
			my $adjust_resol = $resolution;
			my $regSize = $reg_ed - $reg_st + 1;
			my $doTime = 1;
			my $ratio = 3;
			# find proper resol
			if( $resolution >= 1 ){
				$adjust_resol = $regSize;
				$adjust_resol = int( sqrt( $adjust_resol ) * $ratio ) for ( 1 .. $doTime );
			}
			# minimum regShowSize
			my $MinRegShowSize = 15;
			$adjust_resol = min( $adjust_resol, (sprintf "%.2f", $regSize / $MinRegShowSize) );
			# minimum resol is at least 1.0
			$adjust_resol = max( $adjust_resol, 1 );
			# update
			if( $resolution != $adjust_resol ){
				$resolution = $adjust_resol;
			}
		}
		# record
		if( defined $destReg_Href ){
			# only record for destReg
			if( Get_Two_Seg_Olen( $reg_st, $reg_ed, $destReg_Href->{st_pos}, $destReg_Href->{ed_pos} ) ){
				$regResol_Href->{$reg_seg}->{$reg_st} = { st_pos => $reg_st, ed_pos=>$reg_ed, resol=>$resolution, orig_resol=>$orig_resolution };
			}
			else{
				warn "<WARN>\tIgnore $reg_seg:$reg_st-$reg_ed from resolution file.\n";
			}
			# mark destReg edge found or not
			for my $edge_key (qw/ st_pos ed_pos /){
				if(    $reg_st  <= $destReg_Href->{$edge_key}
					&& $reg_ed  >= $destReg_Href->{$edge_key}
				){
					$drawReg_edgeShowBool{$edge_key} = 1;
				}
			}
		}
		else{
			$regResol_Href->{$reg_seg}->{$reg_st} = { st_pos => $reg_st, ed_pos=>$reg_ed, resol=>$resolution, orig_resol=>$orig_resolution };
		}
	}
	close RES;

	# destReg work
	if( defined $destReg_Href ){
		# check destReg is whole covered
		if( $drawReg_edgeShowBool{st_pos} * $drawReg_edgeShowBool{ed_pos} == 0 ){
			warn_and_exit "<ERROR>\tCannot find draw region edge-positions in resolution file.\n";
		}
		# resolution auto adjustment for whole destReg
		if( $autoAdjustResol ){
			# displayed pixel interval size
			my $destRegEdPos_winNO = get_pos_windowNO( 
														pos => $destReg_Href->{ed_pos},
														refseg => $destReg_Href->{refseg},
														regResol_Href => $regResol_Href,
														stpos_1stwin => $destReg_Href->{st_pos}
													);
			# adjust resol in proportion
			my $ratio = $destRegEdPos_winNO / $destRegAutoFigureSize;
			for my $reg_seg (sort keys %$regResol_Href){
				$regResol_Href->{$reg_seg}->{$_}->{resol} = max( 1, sprintf "%.2f", $regResol_Href->{$reg_seg}->{$_}->{resol} * $ratio ) for keys %{$regResol_Href->{$reg_seg}};
			}
		}
	}

	# check region linkage, region should connected with each other precisely
	for my $reg_seg (sort keys %$regResol_Href){
		my @reg_st = sort {$a<=>$b} keys %{$regResol_Href->{$reg_seg}};
		for my $i (0 .. $#reg_st){
			# check connection with previous region
			my ($st_now, $ed_now) = ( $reg_st[$i],   $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{ed_pos} );
			if( $i != 0 ){
				my ($st_pre, $ed_pre) = ( $reg_st[$i-1], $regResol_Href->{$reg_seg}->{$reg_st[$i-1]}->{ed_pos} );
				if( $st_now != $ed_pre + 1 ){
					warn_and_exit "<ERROR>\tfind non-linked regions on $reg_seg: [$st_pre, $ed_pre] and [$st_now, $ed_now]\n"
										 ."\t'pre_end_Pos($ed_pre)' must equal to 'now_start_Pos($st_now) - 1'.\n";
				}
			}
			# report resolution if changes.
			my $used_resol = $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{resol};
			my $orig_resol = $regResol_Href->{$reg_seg}->{$reg_st[$i]}->{orig_resol};
			if( $used_resol != $orig_resol ){
				stout_and_sterr "[INFO]\tAdjust resolution of region ($reg_seg:$st_now-$ed_now) from $orig_resol to $used_resol\n";
			}
		}
	}
}

#--- get window NO of given position ---
sub get_pos_windowNO{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $refseg = $parm{refseg};
	my $pos = $parm{pos};
	my $regResol_Href = $parm{regResol_Href};
	my $stpos_1stwin = $parm{stpos_1stwin};
	my $max_winNO = $parm{maxWinNOtoDraw};
	my $give_winLastPos = $parm{give_winLastPos}; # will give the last pos of the window
	my $winResol_Sref = $parm{winResol_Sref};

	if( $pos < $stpos_1stwin ){
		# warn "The pos to determine location is smaller than start position.\n$pos < $stpos_1stwin\n";
		return 1;
	}

	my $window_NO = 0;
	my $sameWinLastPos = 0;
	for my $reg_st (sort {$a<=>$b} keys %{$regResol_Href->{$refseg}} ){
		my $reg_ed = $regResol_Href->{$refseg}->{$reg_st}->{ed_pos};
		my $resol  = $regResol_Href->{$refseg}->{$reg_st}->{resol};
		next if($reg_ed < $stpos_1stwin);
		my $cal_st = max($reg_st, $stpos_1stwin);
		my $cal_ed = min($reg_ed, $pos);
		my $window_count = POSIX::ceil( ($cal_ed - $cal_st + 1) / $resol );
		$window_NO += $window_count;
		# found the point
		if($cal_ed == $pos){
			$sameWinLastPos = min( $reg_ed, int($cal_st+$window_count*$resol-1) );
			$$winResol_Sref = $resol if(defined $winResol_Sref);
			last;
		}
	}

	# reset to the first window
	$window_NO = 1 if($window_NO <= 0);
	# max window NO control
	$window_NO = $max_winNO if(defined($max_winNO) && $window_NO > $max_winNO);

	return $give_winLastPos
		   ? ($sameWinLastPos, $window_NO)
		   : ($window_NO);
}

#--- get position range represented by given window NO ---
sub windowNO_to_posInterval{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $refseg = $parm{refseg};
	my $windowNO = $parm{windowNO};
	my $regResol_Href = $parm{regResol_Href};
	my $stpos_1stwin = $parm{stpos_1stwin};

	if( $windowNO <= 0 ){
		warn_and_exit "Negative window_NO is not allowed in 'windowNO_to_posInterval' func.\n";
	}

	my $window_sum = 0;
	for my $reg_st (sort {$a<=>$b} keys %{$regResol_Href->{$refseg}} ){
		my $reg_ed = $regResol_Href->{$refseg}->{$reg_st}->{ed_pos};
		my $resol  = $regResol_Href->{$refseg}->{$reg_st}->{resol};
		next if($reg_ed < $stpos_1stwin);
		my $cal_st = max($reg_st, $stpos_1stwin);
		my $cal_ed = $reg_ed;
		my $window_count = POSIX::ceil( ($cal_ed - $cal_st + 1) / $resol );
		if( $window_sum + $window_count >= $windowNO ){
			my $itval_stPos = max( int( $reg_st + ($windowNO - $window_sum - 1) * $resol ), $stpos_1stwin );
			my $itval_edPos = max( int( $itval_stPos + $resol - 1 ), $itval_stPos );
			return ( $itval_stPos, $itval_edPos );
		}
		else{
			$window_sum += $window_count;
		}
	}

	warn_and_exit "<ERROR>\tCannot find pos-interval of given window_NO ($windowNO) via region-resol table in 'windowNO_to_posInterval' func.\n"
						."\treached sum window: $window_sum\n";
}

#--- arrange position of object bilaterally to avoid overlap ---
sub allocate_object_bilater{

	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $obj_Href = $parm{obj_Href};
	my $objIdx_key = $parm{objIdx_key};
	my $newLoc_key = $parm{newLoc_key};
	my $objWid_key = $parm{objWid_key};
	my $obj_gap  = $parm{obj_gap} || 5;
	my $adj_ratio = $parm{adj_ratio} || 0.25;
	my $timeToRpt = $parm{timeToRpt} || 1000;
	my $maxTryTime = $parm{maxTryTime} || 3E5;
	my $disable_L = $parm{disable_L} || 0;
	my $disable_R = ($disable_L) ? 0 : ($parm{disable_R} || 0);

	my $adjust_step = $obj_gap * $adj_ratio;
	my @sorted_objIdx_by_loc = sort {$obj_Href->{$a}->{$objIdx_key}<=>$obj_Href->{$b}->{$objIdx_key}} keys %$obj_Href;
	my $start_idx = 0;
	my $arrange_time = 0;
	ARRANGE: {
		if(++$arrange_time % $timeToRpt == 1){
			stout_and_sterr "$arrange_time time arrangement operation pass. Start from index $start_idx\n";
		}
		if($arrange_time > $maxTryTime){
			stout_and_sterr "$arrange_time time exceeds maximum ($maxTryTime), stop allocation operation now.\n";
			return;
		}
		for my $idx ($start_idx  .. $#sorted_objIdx_by_loc){
			my $objIdx = $sorted_objIdx_by_loc[$idx];
			my $newLoc = $obj_Href->{$objIdx}->{$newLoc_key};
			my $objWid = $obj_Href->{$objIdx}->{$objWid_key};
			my ($L_edge, $R_edge) = ($newLoc - $objWid / 2, $newLoc + $objWid / 2);
			# check the overlap with left object
			if(    !$disable_L
				&& $idx - 1 >= 0
			){
				my $left_objIdx = $sorted_objIdx_by_loc[$idx-1];
				my $left_newLoc = $obj_Href->{$left_objIdx}->{$newLoc_key};
				my $left_locWid = $obj_Href->{$left_objIdx}->{$objWid_key};
				my ($left_L_edge, $left_R_edge) = ($left_newLoc - $left_locWid / 2, $left_newLoc + $left_locWid / 2);
				if($L_edge - $left_R_edge < $obj_gap){
					$obj_Href->{$objIdx}->{$newLoc_key} += $adjust_step;
					$obj_Href->{$left_objIdx}->{$newLoc_key} -= $adjust_step;
					$start_idx = $idx - 1;
					redo ARRANGE;
				}
			}
			# check the overlap with right object
			if(    !$disable_R
				&& $idx + 1 <= $#sorted_objIdx_by_loc
			){
				my $right_objIdx = $sorted_objIdx_by_loc[$idx+1];
				my $right_newLoc = $obj_Href->{$right_objIdx}->{$newLoc_key};
				my $right_locWid = $obj_Href->{$right_objIdx}->{$objWid_key};
				my ($right_L_edge, $right_R_edge) = ($right_newLoc - $right_locWid / 2, $right_newLoc + $right_locWid / 2);
				if($right_L_edge - $R_edge < $obj_gap){
					$obj_Href->{$objIdx}->{$newLoc_key} -= $adjust_step;
					$obj_Href->{$right_objIdx}->{$newLoc_key} += $adjust_step;
					$start_idx = $idx + 1;
					redo ARRANGE;
				}
			}
		}
	}
}

#--- arrange layer of object vertically to avoid overlap ---
sub allocate_object_vertical{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $layer_href = $parm{layer_href};
	   $layer_href = {1=>[]} if( !defined $layer_href );
	my $obj_Href = $parm{obj_Href};
	my $obj_itval_key = $parm{obj_itval_key};
	my $obj_layer_key = $parm{obj_layer_key};
	my $obj_gap = $parm{obj_gap} || 5;
	my $max_layer_Sref = $parm{max_layer_Sref};

	# sort by left position
	my @sorted_obj = sort { 
							$obj_Href->{$a}->{$obj_itval_key}->[0]->[0]
							<=>
							$obj_Href->{$b}->{$obj_itval_key}->[0]->[0]
						  } keys %$obj_Href;
	# no objects
	return if( scalar(@sorted_obj) == 0 );
	# allocation
	for my $obj_name ( @sorted_obj ){
		my $obj_info_Href = $obj_Href->{$obj_name};

		#--- find the proper layer
		my $ok_layer_NO = 0;
		for my $layer_NO (sort {$a<=>$b} keys %$layer_href){
			my $ok_bool = ! &test_overlap_at_one_layer(
														Sbjct_Reg_Aref => $layer_href->{$layer_NO},
														Query_Reg_Aref => $obj_info_Href->{$obj_itval_key},
														obj_gap => $obj_gap
													  );
			if( $ok_bool ){
				$ok_layer_NO = $layer_NO;
				last;
			}
			else{
				$ok_layer_NO = $layer_NO + 1;
			}
		}
		# recode the layer of this obj
		$obj_info_Href->{$obj_layer_key} = $ok_layer_NO;
		# recode usage at this layer
		push @{ $layer_href->{$ok_layer_NO} }, @{ $obj_info_Href->{$obj_itval_key} };
	}
	# update
	$$max_layer_Sref = max( keys %$layer_href );
}

#--- test object overlap at one layer ---
sub test_overlap_at_one_layer{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $Sbjct_Reg_Aref = $parm{Sbjct_Reg_Aref};
	my $Query_Reg_Aref = $parm{Query_Reg_Aref};
	my $obj_gap = $parm{obj_gap} || 5;

	my $no_ovp_bool = 1;
	for my $occu_reg_Aref ( @$Sbjct_Reg_Aref ){
		my $occu_reg_lftEdge = $occu_reg_Aref->[0] - $obj_gap;
		my $occu_reg_rgtEdge = $occu_reg_Aref->[1] + $obj_gap;
		for my $obj_reg_Aref ( @$Query_Reg_Aref ){
			my $ovp_len = Get_Two_Seg_Olen( $occu_reg_lftEdge, $occu_reg_rgtEdge, $obj_reg_Aref->[0], $obj_reg_Aref->[1] );
			if( $ovp_len != 0 ){ # find overlap at this layer
				$no_ovp_bool = 0;
				last;
			}
		}
		# stop
		last if( $no_ovp_bool == 0 );
	}

	return !$no_ovp_bool;
}

#--- 
1; ## tell the perl script the successful access of this module.
