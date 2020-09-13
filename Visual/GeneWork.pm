package FuseSV::Visual::GeneWork;

use strict;
use warnings;
use List::Util qw/ max min sum /;
use SOAPfuse::General_Operation qw/ stout_and_sterr warn_and_exit /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read /;
use FuseSV::Visual::RegionWork qw/ arrange_region_with_clip get_ObjSpanPixItval /;
use SOAPfuse::SVG_Radian_System_Elements qw/ %COLOR_DB draw_circle_seg /;
use SOAPfuse::SVG_Orthogonal_System_Elements qw/ draw_a_parallelogram /;
use SOAPfuse::SVG_Font qw/ show_text_in_line /;
use FuseSV::Visual::ArrangeObj qw/ get_pos_windowNO allocate_object_vertical /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  load_transInfo_from_tpsl
			  add_user_region_to_gene_info
			  allocate_genes_and_draw
			  show_trans
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'GeneWork';
#----- version --------
$VERSION = "0.05";
$DATE = '2018-06-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						load_transInfo_from_tpsl
						add_user_region_to_gene_info
						allocate_genes_and_draw
						show_trans
					 /;

# assume TransInfo_Href has such structure:
# TransInfo_Href = {}
# TransInfo_Href -> $trans_name -> 'gene_name'
# TransInfo_Href -> $trans_name -> 'strand'
# TransInfo_Href -> $trans_name -> 'refseg'
# TransInfo_Href -> $trans_name -> 'left_pos'
# TransInfo_Href -> $trans_name -> 'right_pos'
# TransInfo_Href -> $trans_name -> 'exon_info_Aref'
# TransInfo_Href -> $trans_name -> 'cds_info_Aref'
# TransInfo_Href -> $trans_name -> 'body_drawItval' -> [[a,b], [c,d]]
# TransInfo_Href -> $trans_name -> 'bodyWithLabel_drawItval' -> [[a,b], [c,d]]
# TransInfo_Href -> $trans_name -> 'body_color'
# TransInfo_Href -> $trans_name -> 'show_layer'

#--- load transcript info ---
sub load_transInfo_from_tpsl{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $tpsl_file = $parm{tpsl_file};
	my $TransInfo_Href = $parm{TransInfo_Href};
	my $idx_readReg = $parm{idx_readReg};
	my $tabix = $parm{tabix};
	my $skipGene_Href = $parm{skipGene_Href};
	my $onlyGene_Href = $parm{onlyGene_Href};
	my $clip_min_Pos = $parm{clip_min_Pos}; # effective with $clip_max_Pos
	my $clip_max_Pos = $parm{clip_max_Pos};
	my $cirl_oLen = $parm{cirl_oLen}; # effective with $clip_max_Pos
	my $labelIdx = $parm{labelIdx} || 'g';
	my $fusesv_vmode = $parm{fusesv_vmode} || 0;
	my $g2tid_Href = $parm{g2tid_Href};

	my @clip_opt = ( min_Pos => $clip_min_Pos, max_Pos => $clip_max_Pos, cirl_oLen => $cirl_oLen );

	# read viral gene
	my @openFH_para = ( $tpsl_file );
	push @openFH_para, (idx_readReg=>$idx_readReg, tabix=>$tabix) if( defined $idx_readReg && defined $tabix );
	open (TPSL, Try_GZ_Read(@openFH_para)) || die "fail load trans_psl: $!\n";
	while(<TPSL>){
		my @info = split;
		my $locus_tag    = $info[2];
		my $note         = $info[3];
		my $product      = $info[7];
		my $strand       = $info[8];
		my $trans_name   = $info[9];
		my $refseg       = $info[13];
		my $left_pos     = $info[15] + 1;
		my $right_pos    = $info[16];
		my $exon_len     = $info[18];
		my $cds_info     = $info[19];
		my $exon_leftpos = $info[20];
		my $gene_name    = $info[-1];
		# gene/trans filtration
		next if( defined $skipGene_Href && (exists $skipGene_Href->{$gene_name} || exists $skipGene_Href->{$trans_name}) );
		next if( defined $onlyGene_Href && scalar(keys %$onlyGene_Href) != 0 && !exists $onlyGene_Href->{$gene_name} && !exists $onlyGene_Href->{$trans_name} );
		# meta-info
		$_ = ($_ eq 'N/A')?0:$_ for ($locus_tag, $product, $note);
		# deal with circular extended part, if possible
		my @geneRegion = ( defined $clip_max_Pos )
						 ? arrange_region_with_clip( lftPos => $left_pos, rgtpos => $right_pos, @clip_opt )
						 : ( [ $left_pos , $right_pos ] );
		# filter out-of-range
		next if( scalar(@geneRegion) == 0 );
		## exon info
		my @exon_len = split /,/, $exon_len;
		my @exon_leftpos = split /,/, $exon_leftpos;
		my @exon_info = map {
							  my $exonLFpos = $exon_leftpos[$_] + 1;
							  my $exonRTpos = $exon_leftpos[$_] + $exon_len[$_];
							  ( defined $clip_max_Pos )
							  ? arrange_region_with_clip( lftPos => $exonLFpos, rgtpos => $exonRTpos, @clip_opt )
							  : ( [ $exonLFpos, $exonRTpos ] );
						} (0 .. $#exon_len);
		## CDS info
		my @cds_info = map {
							if(/^(\d+)\((\d+)\)$/){
							  my $cdsLFpos = $1 + 1;
							  my $cdsRTpos = $1 + $2;
							  ( defined $clip_max_Pos )
							  ? arrange_region_with_clip( lftPos => $cdsLFpos, rgtpos => $cdsRTpos, @clip_opt )
							  : ( [ $cdsLFpos, $cdsRTpos ] );
							}
						} split /,/, $cds_info;
		# prepare label text
		if( $fusesv_vmode ){ $trans_name =~ s/^$gene_name\_// }
		my %label_text = ( t=>$trans_name, g=>$gene_name, p=>$product, n=>$note, l=>$locus_tag );
		if( $fusesv_vmode && $labelIdx =~ /^g$/i ){ ($label_text{g}=$trans_name) =~ s/\-\d+$// }
		# record info
		$TransInfo_Href->{$trans_name} = {
											gene_name => $gene_name,
											strand => $strand,
											refseg => $refseg,
											label_text => $label_text{lc($labelIdx)} || $trans_name,
											transBody_Aref => \@geneRegion,
											exon_info_Aref => \@exon_info,
											cds_info_Aref => \@cds_info
										};
		# map trans-id to gene
		if( defined $g2tid_Href ){
			$gene_name =~ s/_ENSG\d{10,}$//;
			(my $trans_id = $trans_name) =~ s/^$gene_name[\_\-]//;
			$g2tid_Href->{$gene_name}->{$trans_id} = 1;
		}
	}
	close TPSL;
}

#--- add the region provide by user ---
sub add_user_region_to_gene_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $user_gene_input_Aref = $parm{user_gene_input_Aref};
	my $TransInfo_Href = $parm{TransInfo_Href};
	my $refseg = $parm{refseg};
	my $clip_min_Pos = $parm{clip_min_Pos}; # effective with $clip_max_Pos
	my $clip_max_Pos = $parm{clip_max_Pos};
	my $cirl_oLen = $parm{cirl_oLen}; # effective with $clip_max_Pos
	my $input_optKey_forAlert = $parm{input_optKey_forAlert};

	my @clip_opt = ( min_Pos => $clip_min_Pos, max_Pos => $clip_max_Pos, cirl_oLen => $cirl_oLen );

	for my $user_gene_input ( @$user_gene_input_Aref ){
		if( $user_gene_input =~ /^(.+?):([\d\-,\+\:]+)$/ ){
			my ($RegName, $RegDetail) = ($1, $2);
			# check name
			warn_and_exit "<ERROR>\tName ($RegName) of '$input_optKey_forAlert' input ($user_gene_input) has been load before.\n"
				for grep $TransInfo_Href->{$_}->{label_text} eq $RegName, keys %$TransInfo_Href;
			# extract strand
			my $strand;
			if( $RegDetail =~ /^([\d\-,]+):([\-\+])$/ ){
				$RegDetail = $1;
				$strand = $2;
			}
			# check region
			my @this_user_gene_region;
			my @this_user_exon_region;
			for my $interval (split /,/, $RegDetail){
				my ($st_pos, $ed_pos) = ($interval =~ /^(\d+)-(\d+)$/);
				if(    !defined $st_pos || !defined $ed_pos
					|| $st_pos > $ed_pos
				){
					warn_and_exit "<ERROR>\tWrong region fromat of '$input_optKey_forAlert' input ('$user_gene_input'; [$interval]).\n";
				}
				# clip region
				my @this_interval_exon =
					 defined $clip_max_Pos
					 ? arrange_region_with_clip( lftPos => $st_pos, rgtpos => $ed_pos, @clip_opt )
					 : ( [ $st_pos , $ed_pos ] );
				# record as exon
				push @this_user_exon_region, @this_interval_exon;
				# record in gene body
				if( scalar(@this_user_gene_region) == 0 ){
					push @this_user_gene_region, @this_interval_exon;
				}
				else{
					for my $exon_Aref (@this_interval_exon){
						if(    defined $cirl_oLen
							&& $this_user_gene_region[-1]->[1] == $cirl_oLen
						){
							push @this_user_gene_region, $exon_Aref;
						}
						else{
							# re-assign Aref of new array!
							$this_user_gene_region[-1] = [ $this_user_gene_region[-1]->[0], $exon_Aref->[1] ];
						}
					}
				}
			}
			# filter out-of-range
			next if( scalar(@this_user_exon_region) == 0 );
			# load user-gene info
			$TransInfo_Href->{$RegName} = {
											gene_name => $RegName,
											strand => $strand || 'N/A',
											refseg => $refseg,
											label_text => $RegName,
											transBody_Aref => \@this_user_gene_region,
											exon_info_Aref => \@this_user_exon_region,
											cds_info_Aref => []
										  };
		}
		else{
			warn_and_exit "<ERROR>\tWrong fromat of '$input_optKey_forAlert' input ($user_gene_input).\n";
		}
	}
}

#--- allocate genes in vertical layer and draw
sub allocate_genes_and_draw{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $TransInfo_Href = $parm{TransInfo_Href};
	my $pos2wNO_para_Aref = $parm{pos2wNO_para_Aref};
	my $axisZX_Sref = $parm{axisZX_Sref};
	  $$axisZX_Sref = 100 if( !defined $$axisZX_Sref );
	my $axisZY_Sref = $parm{axisZY_Sref};
	  $$axisZY_Sref = 100 if( !defined $$axisZY_Sref );
	my $tnameFtfm = $parm{tnameFtfm} || 'Arial';
	my $tnameFtsz = $parm{tnameFtsz} || 12;
	my $ExonMinWidth = $parm{ExonMinWidth} || 1;

	# basic para
	my $transLayer_height = 10;
	my $intron_height = 1;
	my $exon_height = 4;
	my $CDS_height = 8;
	my $gene_gap = 10;
	my $gene_label_gap = 5;
	my $trans_st_y = $$axisZY_Sref + $transLayer_height / 2;

	#--- calculate the span interval of each transcript
	get_ObjSpanPixItval(
							ObjInfoPool_Href => $TransInfo_Href,
							pos2xAxisPos_para_Aref => $pos2wNO_para_Aref,
							axisZX => $$axisZX_Sref,
							labelFtfm => $tnameFtfm,
							labelFtsz => $tnameFtsz,
							labelHeight => $CDS_height,
							bodyLabelGap => $gene_label_gap,
							bodyRegAref_key => 'transBody_Aref',
							bodyPixItval_key => 'body_drawItval',
							bodyWithLabelPixItval_key => 'bodyWithLabel_drawItval',
							labelText_key => 'label_text',
							bodyCol_key => 'body_color'
						);

	#--- arrange the Y location of transcript
	my $max_show_layerNO = 0;
	allocate_object_vertical(
								obj_Href => $TransInfo_Href,
								obj_itval_key => 'bodyWithLabel_drawItval',
								obj_layer_key => 'show_layer',
								obj_gap => $gene_gap,
								max_layer_Sref => \$max_show_layerNO
							);
	#--- draw gene
	&show_trans(
				SVG_obj_Oref => $SVG_obj_Oref,
				TransInfo_Href => $TransInfo_Href,
				pos2wNO_para_Aref => $pos2wNO_para_Aref,
				axisZX => $$axisZX_Sref,
				trans_st_y => $trans_st_y,
				tnameFtfm => $tnameFtfm,
				tnameFtsz => $tnameFtsz,
				tnameHeight => $CDS_height,
				bodyLabelGap => $gene_label_gap,
				transLayer_height => $transLayer_height,
				intron_height => $intron_height,
				exon_height => $exon_height,
				CDS_height => $CDS_height,
				ExonMinWidth => $ExonMinWidth
			   );

	# update the axis Y location
	$$axisZY_Sref = $trans_st_y + $max_show_layerNO * $transLayer_height + 3;
}

#--- draw trans
sub show_trans{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $SVG_obj_Oref = $parm{SVG_obj_Oref};
	my $TransInfo_Href = $parm{TransInfo_Href};
	my $pos2wNO_para_Aref = $parm{pos2wNO_para_Aref};
	my $axisZX = $parm{axisZX};
	   $axisZX = 100 if( !defined $axisZX );
	my $trans_st_y = $parm{trans_st_y} || 50;
	my $tnameFtfm = $parm{tnameFtfm} || 'Arial';
	my $tnameFtsz = $parm{tnameFtsz} || 12;
	my $bodyLabelGap = $parm{bodyLabelGap} || 5;
	my $transLayer_height = $parm{transLayer_height} || 10;
	my $intron_height = $parm{intron_height} || 1;
	my $exon_height = $parm{exon_height} || 4;
	my $CDS_height = $parm{CDS_height} || 8;
	my $tnameHeight = $CDS_height;
	my $ExonMinWidth = $parm{ExonMinWidth} || 1;

	for my $trans_name ( sort keys %$TransInfo_Href ){
		my $transInfo_Href = $TransInfo_Href->{$trans_name};
		my $trans_y = $trans_st_y + $transLayer_height * ( $transInfo_Href->{show_layer} - 0.5 );
		#--- intron line, gene name
		for my $reg_Aref ( @{$transInfo_Href->{body_drawItval}} ){
			my $body_lftPix = $reg_Aref->[0];
			my $body_rgtPix = $reg_Aref->[1];
			my $body_PixSum = $body_rgtPix - $body_lftPix + 1;
			my $body_MidLoc = ( $body_lftPix + $body_rgtPix ) / 2;
			# draw intron
			draw_a_parallelogram(
									$SVG_obj_Oref,
									x => $body_MidLoc,
									y => $trans_y,
									fill_color => $transInfo_Href->{body_color},
									head_bottom_side_len => $body_PixSum,
									left_right_side_len => $intron_height,
									boundary_color => 'none'
								);
			# show gene name
			show_text_in_line(
								$SVG_obj_Oref,
								text_x => $body_lftPix - $bodyLabelGap,
								text_y => $trans_y,
								text => $transInfo_Href->{label_text},
								font_family => $tnameFtfm,
								font_size => $tnameFtsz,
								text_col => 'black',
								text_anchor => 'end',
								height_adjust => 1,
								height_limit => $tnameHeight,
								width_limit => 0,
								rotate_degree => 0,
								draw_bool => 1
							);
		}

		#--- show CDS, if available
		for my $cds_Aref (@{$transInfo_Href->{cds_info_Aref}}){
			my $cds_lftPix = $axisZX + get_pos_windowNO( pos => $cds_Aref->[0], @$pos2wNO_para_Aref );
			my $cds_rgtPix = $axisZX + get_pos_windowNO( pos => $cds_Aref->[1], @$pos2wNO_para_Aref );
			my $cds_spanPix = max( abs($cds_rgtPix-$cds_lftPix)+1, $ExonMinWidth );
			draw_a_parallelogram(
									$SVG_obj_Oref,
									x => ($cds_lftPix+$cds_rgtPix)/2,
									y => $trans_y,
									fill_color => $transInfo_Href->{body_color},
									head_bottom_side_len => $cds_spanPix,
									left_right_side_len => $CDS_height,
									boundary_color => 'none'
								);
		}

		#--- show exons
		my $trans_strand = $transInfo_Href->{strand};
		my $minSpaceForArrow = 10;
		my $strand_shown = undef;
		my $strand_arrow_col = ($transInfo_Href->{body_color} =~ /$COLOR_DB{dark_col}/) ? 'white' : 'black';
		my $strand_arrow_ori = ($trans_strand eq '+') ? 'clockwise' : 'anti-clockwise';
		my $intronSizeToLoc_Href = {}; # to store intron region for possible strand-arrow display
		my $last_exon_rgtPix = undef;
		for my $exon_Aref (@{$transInfo_Href->{exon_info_Aref}}){
			my $exon_lftPix = $axisZX + get_pos_windowNO( pos=>$exon_Aref->[0], @$pos2wNO_para_Aref );
			my $exon_rgtPix = $axisZX + get_pos_windowNO( pos=>$exon_Aref->[1], @$pos2wNO_para_Aref );
			my $exon_spanPix = max( abs($exon_rgtPix-$exon_lftPix)+1, $ExonMinWidth );
			my $mid_x = ($exon_lftPix + $exon_rgtPix) / 2;
			my $mid_y = $trans_y;
			my $radius = 1E6;
			my $rad_size = $exon_spanPix / $radius;
			# strand ?
			my @strand_arrow_parm;
			if( $trans_strand =~ /^[\-\+]$/ ){
				my $arrow_num = ($exon_spanPix >= 20) ? 2 : ($exon_spanPix >= $minSpaceForArrow ? 1 : 0);
				$strand_shown ||= $arrow_num;
				@strand_arrow_parm = (
										arrow_num => $arrow_num,
										arrow_color => $strand_arrow_col,
										arrow_ori => $strand_arrow_ori,
										arrow_span_perct => max( 4, 6 / $exon_spanPix * 100 ),
										arrow_maxArcLen => 10
									);
			}
			# draw
			draw_circle_seg(
								$SVG_obj_Oref,
								cx => $mid_x,
								cy => $mid_y + $radius,
								start_rad => $rad_size / 2 * -1,
								rad_size => $rad_size,
								inner_radius => $radius - $exon_height / 2,
								outer_radius => $radius + $exon_height / 2,
								seg_fill_color => $transInfo_Href->{body_color},
								seg_boundary_color => 'none',
								@strand_arrow_parm
							);
			# record of intron
			if( defined $last_exon_rgtPix ){
				my $intron_PixSum = $exon_lftPix - $last_exon_rgtPix + 1;
				my $intron_MidLoc = ( $exon_lftPix + $last_exon_rgtPix ) / 2;
				push @{$intronSizeToLoc_Href->{int($intron_PixSum)}}, $intron_MidLoc;
			}
			$last_exon_rgtPix = $exon_rgtPix;
		}

		#--- make up, when no strand-arrow is shown
		if( defined $strand_shown && $strand_shown == 0 ){
			my $shownArrowNum = 0;
			my $radius = 1E6;
			my $rad_size = 8 / $radius;
			for my $intron_PixSum (sort {$b<=>$a} keys %$intronSizeToLoc_Href){
				next if( $intron_PixSum < $minSpaceForArrow );
				for my $intronArrowLoc ( @{$intronSizeToLoc_Href->{$intron_PixSum}} ){
					draw_circle_seg(
										$SVG_obj_Oref,
										cx => $intronArrowLoc,
										cy => $trans_y + $radius,
										start_rad => $rad_size / 2 * -1,
										rad_size => $rad_size,
										inner_radius => $radius - $exon_height / 2,
										outer_radius => $radius + $exon_height / 2,
										seg_fill_color => 'none',
										seg_boundary_color => 'none',
										# arrow setting
										arrow_num => 1,
										arrow_color => $transInfo_Href->{body_color},
										arrow_ori => $strand_arrow_ori,
										arrow_span_perct => 100,
										arrow_maxArcLen => 8,
										arrow_draw_only => 1
									);
					$shownArrowNum ++;
					last if( $shownArrowNum == 3 );
				}
				last if( $shownArrowNum == 3 );
			}
		}

		# inform
		stout_and_sterr `date`."[INFO]:\tDraw viral gene $trans_name at $transInfo_Href->{show_layer} layer.\n";
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
