package FuseSV::Visual::RegionWork;

use strict;
use warnings;
use List::Util qw/max min sum/;
use SOAPfuse::OpenFile qw/Try_GZ_Read/;
use SOAPfuse::SVG_Radian_System_Elements qw/ %COLOR_DB /;
use SOAPfuse::SVG_Font qw/ show_text_in_line /;
use FuseSV::Visual::ArrangeObj qw/ get_pos_windowNO  /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  arrange_region_with_clip
			  get_ObjSpanPixItval
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'RegionWork';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-02-12';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						arrange_region_with_clip
						deal_circular_extended_part
						clip_region_edge
						get_ObjSpanPixItval
					 /;

#--- arrange one region in expected region ---
sub arrange_region_with_clip{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $lftPos = $parm{lftPos};
	my $rgtpos = $parm{rgtpos};
	my $min_Pos = $parm{min_Pos};
	my $max_Pos = $parm{max_Pos};
	my $cirl_oLen = $parm{cirl_oLen};

	# deal with circular extended part, if possible
	my @candReg_Aref = ( defined $cirl_oLen )
					  ? &deal_circular_extended_part( lftPos => $lftPos, rgtpos => $rgtpos, origLen => $cirl_oLen )
					  : ( [ $lftPos, $rgtpos ] );
	# clip region
	my @Region;
	for my $Aref ( @candReg_Aref ){
		$Aref->[0] = &clip_region_edge( LR_Edge => 'L', edgePos => $Aref->[0], min_Pos => $min_Pos, max_Pos => $max_Pos );
		$Aref->[1] = &clip_region_edge( LR_Edge => 'R', edgePos => $Aref->[1], min_Pos => $min_Pos, max_Pos => $max_Pos );
		# record
		push @Region, $Aref if( $Aref->[0] < $Aref->[1] );
	}

	return @Region;
}

#--- deal with the extended region ---
#--- such as for the circular virus
sub deal_circular_extended_part{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $lftPos = $parm{lftPos};
	my $rgtpos = $parm{rgtpos};
	my $origLen = $parm{origLen};

	if( $lftPos > $origLen ){
		return [ $lftPos - $origLen, $rgtpos - $origLen ];
	}
	elsif(   $lftPos < $origLen
		  && $rgtpos > $origLen
	){
		return ( [ $lftPos,          $origLen ],
				 [      1, $rgtpos - $origLen ]  );
	}
	else{
		return [ $lftPos, $rgtpos ];
	}
}

#--- clip region edge based on allowed min/max position ---
sub clip_region_edge{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $edgePos = $parm{edgePos};
	my $LR_Edge = $parm{LR_Edge};
	my $min_Pos = $parm{min_Pos};
	my $max_Pos = $parm{max_Pos};

	if( $LR_Edge =~ /^L/i ){
		return min( max( $edgePos, $min_Pos ), $max_Pos );
	}
	else{
		return max( min( $edgePos, $max_Pos ), $min_Pos );
	}
}

#--- calculate the span interval of each object ---
sub get_ObjSpanPixItval{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $ObjInfoPool_Href = $parm{ObjInfoPool_Href};
	my $pos2xAxisPos_para_Aref = $parm{pos2xAxisPos_para_Aref};
	my $axisZX = $parm{axisZX};
	   $axisZX = 100 if( !defined $axisZX );
	my $labelFtsz = $parm{labelFtsz} || 12;
	my $labelFtfm = $parm{labelFtfm} || 'Arial';
	my $labelHeight = $parm{labelHeight} || 15;
	my $bodyLabelGap = $parm{bodyLabelGap} || 5;
	my $bodyRegAref_key = $parm{bodyRegAref_key};
	my $bodyPixItval_key = $parm{bodyPixItval_key};
	my $bodyWithLabelPixItval_key = $parm{bodyWithLabelPixItval_key};
	my $labelText_key = $parm{labelText_key};
	my $bodyCol_key = $parm{bodyCol_key};
	my $bodyColRef_key = $parm{bodyColRef_key}; # for assign color if set, or else auto

	my %col_usage;
	for my $obj_name ( sort keys %$ObjInfoPool_Href ){
		my $objInfo_Href = $ObjInfoPool_Href->{$obj_name};
		#--- get gene draw X region
		for my $bodyRegAref (@{$objInfo_Href->{$bodyRegAref_key}}){
			# draw interval of object body
			my $body_lftPix = $axisZX + get_pos_windowNO( pos => $bodyRegAref->[0], @$pos2xAxisPos_para_Aref );
			my $body_rgtPix = $axisZX + get_pos_windowNO( pos => $bodyRegAref->[1], @$pos2xAxisPos_para_Aref );
			push @{$objInfo_Href->{$bodyPixItval_key}}, [ $body_lftPix, $body_rgtPix ];
			# label interval
			if(    defined $bodyWithLabelPixItval_key
				&& defined $labelText_key
			){
				my $temp_SVG_obj = 'temp';
				my $text_size_Aref = show_text_in_line(
														\$temp_SVG_obj,
														text_x => $body_lftPix - $bodyLabelGap,
														text_y => 100, # temp
														text => $objInfo_Href->{$labelText_key},
														font_family => $labelFtfm,
														font_size => $labelFtsz,
														text_col => 'black',
														text_anchor => 'end',
														height_adjust => 1,
														height_limit => $labelHeight,
														width_limit => 0,
														rotate_degree => 0,
														draw_bool => 0, # do not draw
													);
				# full width of this object info
				my $bodyLabel_lftPix = $body_lftPix - $bodyLabelGap - $text_size_Aref->[1];
				my $bodyLabel_rgtPix = $body_rgtPix;
				push @{$objInfo_Href->{$bodyWithLabelPixItval_key}}, [ $bodyLabel_lftPix, $bodyLabel_rgtPix ];
			}
		}
		#--- set object color
		if( defined $bodyCol_key ){
			my $objBody_col;
			if ( defined $bodyColRef_key ){
				my $bodyColRef = $objInfo_Href->{$bodyColRef_key};
				if( exists $col_usage{$bodyColRef} ){
					$objBody_col = $col_usage{$bodyColRef};
				}
				else{
					my %bodyColUsed  = map { ($_,1) } values %col_usage;
					my @COLOR_DB_idx = map { $COLOR_DB{$_} } sort {$a<=>$b} grep /^\d+$/ && !exists($bodyColUsed{$COLOR_DB{$_}}), keys %COLOR_DB;
					if( @COLOR_DB_idx == 0 ){
						@COLOR_DB_idx = map { $COLOR_DB{$_} } sort {$a<=>$b} grep /^\d+$/, keys %COLOR_DB;
						$objBody_col = $COLOR_DB_idx[int(rand(scalar(@COLOR_DB_idx)))];
					}
					else{
						$objBody_col = $COLOR_DB_idx[0];
					}
					# record
					$col_usage{$bodyColRef} = $objBody_col;
				}
			}
			else{
				my @COLOR_DB_idx = grep /^\d+$/, sort keys %COLOR_DB;
				my $count = sum( values %col_usage ) || 0;
				$objBody_col = $COLOR_DB{ $COLOR_DB_idx[ $count % (scalar(@COLOR_DB_idx) || 1) ] };
				# record
				$col_usage{$objBody_col} ++;
			}
			# assign color
			$objInfo_Href->{$bodyCol_key} = $objBody_col;
		}
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
