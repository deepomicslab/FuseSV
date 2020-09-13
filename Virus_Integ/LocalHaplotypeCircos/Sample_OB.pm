package FuseSV::Virus_Integ::LocalHaplotypeCircos::Sample_OB;

use strict;
use warnings;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Sample_OB';
#----- version --------
$VERSION = "0.03";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
							new
							load_subregion_info
							get_segid_details
							load_haplotype_contig
							load_path_arc_info
							find_previous_arc
							find_PNC_deg_char
							loadup_pos_features
							loadup_rna_reads_amount
						/;

#--- structure of object
# seg_OB -> NO = seg_NO
# sam_OB -> {sam_id}
# sam_OB -> {host} -> {reg_id (H1,H2)} -> {refseg, lfpos, rtpos, copy_number}
# sam_OB -> {virus} -> {reg_id (V1,V2)} -> {refseg, lfpos, rtpos}
# sam_OB -> {pathcolour} -> {$hap_NO}
# sam_OB -> {haplotype} -> {$hap_NO} -> {contig_NO} -> [ repeat_time,
#                         [ [ SegNO,refseg,reg_st,reg_ed,+/-1,PNC_NO ],
#                         ... ]
#                       ]
# sam_OB -> {path_arc} -> {$hap_NO} -> [ { last_arc=>info_Href, this_arc=>info_Href,
#                 PNC_text=>'', PNC_NO=>[]
#                 }, ... ]
#        :arc_info_Href -> { seg_NO=>$seg_NO, st_rad=>$st_rad, ed_rad=>$ed_rad, ori=>$ori, layer_NO=>$layer_NO }
# sam_OB -> {PNC_deg_char} -> {$hap_NO} -> {char1=>[PNC_NO1,PNC_NO2,..], char2=>[PNC_NO3,PNC_NO4,..], ...}
# sam_OB -> {feature} -> {partner} -> {chr:pos} -> [ feature_1, feature_2, feature_3, ... ];
# sam_OB -> {rnaReads_Amount}

#--- construction of object
sub new{
	shift;
	my %parm = @_;
	my $sam_id = $parm{sam_id};

	my $sam_OB = {};
	$sam_OB->{sam_id} = $sam_id;
	# host sub-region details
	$sam_OB->{host} = {};
	# virus sub-region details
	$sam_OB->{virus} = {};
	# haplotype details
	$sam_OB->{haplotype} = {};
	# path arc details
	$sam_OB->{path_arc} = {};

	bless($sam_OB);
	return $sam_OB;
}

#--- load in the sub region details
sub load_subregion_info{
	my $sam_OB = shift;
	my %parm = @_;

	# check valid?
	my $segtype = $parm{segtype};
	unless( exists($sam_OB->{$segtype}) ){
		warn_and_exit "Cannot recognise the segtype $segtype. 'host' or 'virus'?\n";
	}

	# check existence
	my $segid = $parm{segid};
	if ( exists($sam_OB->{$segtype}->{$segid}) ){
		stout_and_sterr "<WARN>\tEncounter $segtype region $segid of sample $sam_OB->{sam_id} the sceond time! Merge!\n";
		# check
		if( $parm{refseg} ne $sam_OB->{$segtype}->{$segid}->{refseg} ){
			warn_and_exit "Different refseg!\n";
		}
		if( $parm{copy_number} ne $sam_OB->{$segtype}->{$segid}->{copy_number} ){
			warn_and_exit "Different Copy Number!\n";
		}
		# update
		$sam_OB->{$segtype}->{$segid}->{lfpos} = min($sam_OB->{$segtype}->{$segid}->{lfpos}, $parm{lfpos});
		$sam_OB->{$segtype}->{$segid}->{rtpos} = max($sam_OB->{$segtype}->{$segid}->{rtpos}, $parm{rtpos});
	}
	else{ # load in details to new
		$sam_OB->{$segtype}->{$segid} = {
											refseg => $parm{refseg},
											lfpos => $parm{lfpos},
											rtpos => $parm{rtpos},
											copy_number => $parm{copy_number}
										};
	}
}

#--- return info of given segtype's segid
sub get_segid_details{
	my $sam_OB = shift;
	my %parm = @_;
	my $segtype = $parm{segtype};
	my $segid  = $parm{segid};

	# check existence
	unless( exists($sam_OB->{$segtype}) ){
		warn_and_exit "Cannot recognise the segtype $segtype. 'host' or 'virus'?\n";
	}
	unless( exists($sam_OB->{$segtype}->{$segid}) ){
		warn_and_exit "Cannot find information for $segtype $segid region of sample $sam_OB->{sam_id}.\n";
	}

	return ( $sam_OB->{$segtype}->{$segid}->{refseg},
			 $sam_OB->{$segtype}->{$segid}->{lfpos},
			 $sam_OB->{$segtype}->{$segid}->{rtpos}   );
}

#--- load in the contig details in given haplotype path
sub load_haplotype_contig{
	my $sam_OB = shift;
	my %parm = @_;
	my $hap_NO = $parm{hap_NO};
	my $hap_col = $parm{colour};
	my $ctg_NO = $parm{ctg_NO};
	my $repeat = $parm{repeat};
	my $reg_Aref = $parm{reg_Aref};

	# record
	$sam_OB->{haplotype}->{$hap_NO}->{$ctg_NO} = [ $repeat, $reg_Aref ];
	# this haplo-path colour
	if( exists($sam_OB->{pathcolour}->{$hap_NO}) ){
		if( $sam_OB->{pathcolour}->{$hap_NO} ne $hap_col ){
			warn_and_exit "Applied different colour for NO.$hap_NO path of sample $sam_OB->{sam_id}.\n";
		}
	}
	else{
		$sam_OB->{pathcolour}->{$hap_NO} = $hap_col;
	}
}

#--- load in path arc details
sub load_path_arc_info{
	my $sam_OB = shift;
	my %parm = @_;
	my $hap_NO = $parm{hap_NO};
	my $last_arc_Href = $parm{last_arc};
	my $this_arc_Href = $parm{this_arc};
	my $PNC_text = $parm{PNC_text};

	my $arc_Href = {};
	$arc_Href->{last_arc}->{$_} = $last_arc_Href->{$_} for keys %$last_arc_Href;
	$arc_Href->{this_arc}->{$_} = $this_arc_Href->{$_} for keys %$this_arc_Href;
	$arc_Href->{PNC_text} = $PNC_text; # later, this may change to the PNC_deg_char
	$arc_Href->{PNC_NO} = [$PNC_text]; # now, this only contains the PNC_NO
	# record
	push @{$sam_OB->{path_arc}->{$hap_NO}}, $arc_Href;
}

#--- find the previous arranged arc represents as same as the given arc
sub find_previous_arc{
	my $sam_OB = shift;
	my %parm = @_;
	my $hap_NO = $parm{hap_NO};
	my $last_arc_Href = $parm{last_arc};
	my $this_arc_Href = $parm{this_arc};
	my $PNC_NO = $parm{PNC_NO};
	my $record = $parm{record};
	# specific options
	my $draw_lastarc_radian_bool_Sref = $parm{last_radline} || 'NA';
	my $draw_thisarc_radian_bool_Sref = $parm{this_radline} || 'NA';
	my $draw_bezier_bool_Sref = $parm{bezier} || 'NA';

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};

	my $layer_NO = 0;
	for my $arc_Href ( @{$sam_OB->{path_arc}->{$hap_NO}} ){
		if(   (    $this_arc_Href->{seg_NO} == $arc_Href->{this_arc}->{seg_NO} # this_arc must be at same segments
						&& $this_arc_Href->{st_rad} == $arc_Href->{this_arc}->{st_rad} # this_arc must start at the same pos
						&& $this_arc_Href->{ed_rad} == $arc_Href->{this_arc}->{ed_rad} # this_arc must end at the same pos
					)
			 && (   $last_arc_Href->{seg_NO} != $this_arc_Href->{seg_NO} # last_arc's seg is different from this_seg's
					 || (    $arc_Href->{this_arc}->{seg_NO} == $arc_Href->{last_arc}->{seg_NO}
						&& $last_arc_Href->{ed_rad} == $arc_Href->{last_arc}->{ed_rad}   # last_arc must end at the same pos
					&& (   $last_arc_Href->{layer_NO} == $arc_Href->{last_arc}->{layer_NO} # last_arc must be at the same layer
										|| (    $this_arc_Href->{seg_NO} == 0 # special for virus loop junction
											 && (    (    $Segment_Href->{0}->rad2perc(rad => $last_arc_Href->{ed_rad} ) >= 99
															 && $Segment_Href->{0}->rad2perc(rad => $this_arc_Href->{st_rad} ) <= 1
														 )
													||
														 (    $Segment_Href->{0}->rad2perc( rad => $last_arc_Href->{ed_rad} ) <= 1
															 && $Segment_Href->{0}->rad2perc( rad => $this_arc_Href->{st_rad} ) >= 99
														 )
												)
											 )
									 )
							)
					)
			){
			# find the proper arc
			$layer_NO = $arc_Href->{this_arc}->{layer_NO};
			# record?
			if( $record ){
				push @{$arc_Href->{PNC_NO}}, $PNC_NO;
				# try to give a PNC degenerated char
				my $PNC_text = $arc_Href->{PNC_text};
				if( $PNC_text eq $arc_Href->{PNC_NO}->[0] ){ # still is PNC NO, not deg_char
					my $assign_bool = 0;
					for my $PNC_deg_char ( @{$V_Href->{PNC_deg_char}} ){
						if( !exists($sam_OB->{PNC_deg_char}->{$hap_NO}->{$PNC_deg_char}) ){ # available PNC_deg_char
							$arc_Href->{PNC_text} = $PNC_deg_char; # use this one
							$assign_bool = 1;
							# record this PNC_deg_char has been used
							$sam_OB->{PNC_deg_char}->{$hap_NO}->{$PNC_deg_char} = [@{$arc_Href->{PNC_NO}}];
							last;
						}
					}
					# check
					unless( $assign_bool ){
						warn_and_exit "<ERROR>\tNo new PNC_deg_char for sample $sam_OB->{sam_id} NO.$hap_NO haplotype NO.$PNC_NO path.\n";
					}
				}
				else{
					push @{$sam_OB->{PNC_deg_char}->{$hap_NO}->{$PNC_text}}, $PNC_NO;
				}
				# alert
				stout_and_sterr "[INFO]\tFind previous arc NO.$arc_Href->{PNC_NO}->[0] (char:$arc_Href->{PNC_text}) for sample $sam_OB->{sam_id} NO.$hap_NO haplotype NO.$PNC_NO path.\n";
			}
			# what to return
			if( $last_arc_Href->{seg_NO} == $this_arc_Href->{seg_NO} ){ # last_arc's seg is identical to this_arc's
				return $layer_NO;
			}
			else{ # last_arc's seg is different from this_arc's
				if( $arc_Href->{last_arc}->{seg_NO} != $arc_Href->{this_arc}->{seg_NO} ){ # have drawn this_arc's radian line
					$$draw_thisarc_radian_bool_Sref = 0; # no needs to draw this_arc's radian line
				}
				if( $last_arc_Href->{seg_NO} == $arc_Href->{last_arc}->{seg_NO} ){ # same last_seg
					if( $last_arc_Href->{ed_rad} == $arc_Href->{last_arc}->{ed_rad} ){ # end at same pos
						$$draw_bezier_bool_Sref = 0; # no needs to draw bezier line
						if( $last_arc_Href->{layer_NO} == $arc_Href->{last_arc}->{layer_NO} ){ # locate at same layer
							$$draw_lastarc_radian_bool_Sref = 0; # no need to draw last_arc's radian line
						}
					}
				}
				return $layer_NO;
			}
		}
	}

	# not found
	return $layer_NO;
}

#--- find characters for group PNCs (degenerate)
sub find_PNC_deg_char{
	my $sam_OB = shift;
	my %parm = @_;
	my $hap_NO = $parm{hap_NO};
	my $PNC_NO = $parm{PNC_NO};

	for my $arc_Href ( @{$sam_OB->{path_arc}->{$hap_NO}} ){
		if( $arc_Href->{PNC_NO}->[0] == $PNC_NO ){
			return $arc_Href->{PNC_text};
		}
	}

	warn_and_exit "Cannot find path_arc NO.$PNC_NO in haplotype NO.$hap_NO of sample $sam_OB->{sam_id}.\n";
}

#--- loadup position features
sub loadup_pos_features{
	my $sam_OB = shift;
	my %parm = @_;
	my $partner = $parm{partner};
	my $integ_case = $parm{integ_case};
	my $feature_id = $parm{feature_id};

	push @{$sam_OB->{feature}->{$partner}->{$integ_case}}, $feature_id;
}

#--- loadup rna sequenced reads amount
sub loadup_rna_reads_amount{
	my $sam_OB = shift;
	my %parm = @_;
	$sam_OB->{rnaReads_Amount} = $parm{reads_amount};
}

1; ## tell the perl script the successful access of this module.
