package FuseSV::Virus_Integ::LocalHaplotypeCircos::Haplotype;

use strict;
use warnings;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use BioFuse::Visual::SVG_Util::RadSys qw/get_coordinate_on_circle/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
        load_local_haplotype
      /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'Haplotype';
#----- version --------
$VERSION = "0.04";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
            load_local_haplotype
            transform_region_string
            regionid2region
            bezier_between_segments
            PNC_text_size
           /;

#--- load local haplotype of sample(s) ---
sub load_local_haplotype{

  # variants
  my $Sample_Href = $V_Href->{Sample_Href};
  my $Sample_LocalHaplotype_file = $V_Href->{Sample_LocalHaplotype_file};

  my %applied_path_col;
  open (LCHP,$Sample_LocalHaplotype_file) || die "fail open $Sample_LocalHaplotype_file: $!\n";
  # theme
  (my $theme_line = lc(<LCHP>)) =~ s/^#//;
  my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
  while (<LCHP>){
    next if(/^\#/);
    chomp;
    my @info = split /\t+/; # use '\t+' to split.
    my %hap_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
    my $sam_id = $hap_Ob{sampleid};
    my $hap_NO = $hap_Ob{haplotype_no};
    my $hap_col = $hap_Ob{colour};
    my $ctg_NO = $hap_Ob{contig_no};
    my $repeat = $hap_Ob{repeat_time};
    my $regidstr = $hap_Ob{regid_string};
    # check sample
    unless( exists($Sample_Href->{$sam_id}) ){
      next;
    }
    # check haplotype colour, overlap?
    my $sam_hap = "$sam_id,$hap_NO";
    if( exists($applied_path_col{$hap_col}) ){
      if( $applied_path_col{$hap_col} ne $sam_hap ){
        warn_and_exit "colour '$hap_col' has been applied for $applied_path_col{$hap_col} path.\n";
      }
    }
    else{
      $applied_path_col{$hap_col} = $sam_hap;
    }
    # transform region_id string to location string
    my $reg_Aref = &transform_region_string( sam_id => $sam_id, regidstr => $regidstr );
    # load haplotype info
    my $sam_OB = $Sample_Href->{$sam_id};
    $sam_OB->load_haplotype_contig(
                    hap_NO => $hap_NO,
                    colour => $hap_col,
                    ctg_NO => $ctg_NO,
                    repeat => $repeat,
                    reg_Aref => $reg_Aref
                  );
  }
  close LCHP;
}

#--- deal the region string of contig in haplotype
sub transform_region_string{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $sam_id = $parm{sam_id};
  my $regidstr = $parm{regidstr};

  # region id => region
  my @sub_region;
  &regionid2region( sam_id => $sam_id, regidstr => $regidstr, subreg_Aref => \@sub_region );

  # merge flanking sub region
  # instance: push @$subregion_Aref, [$hostseg_NO, $refseg, $st_pos, $ed_pos, $orientation];
  my @return_subreg;
  my @last_subregion = @{ $sub_region[0] };
  for my $i (1 .. $#sub_region){
    my $this_subreg_Aref = $sub_region[$i];
    if(    $last_subregion[0] == $this_subreg_Aref->[0] # at same radius segment
      && $last_subregion[4] == $this_subreg_Aref->[4] # same orientation
      && $last_subregion[3] == $this_subreg_Aref->[2] - ( $last_subregion[4] == 1 ? 1 : -1 ) # connected
      ){
        ## these two regions could be merged, extends subregion
        ## other information are same
      stout_and_sterr "[INFO]\tconnect $last_subregion[1] pos.$last_subregion[3] and pos.$this_subreg_Aref->[2] in sample $sam_id\n";
      $last_subregion[3] = $this_subreg_Aref->[3];
    }
    else{
      push @return_subreg, [ @last_subregion ];
      @last_subregion = @{ $this_subreg_Aref };
    }
  }
  # deal with the last one
  push @return_subreg, [ @last_subregion ];

  return \@return_subreg;
}

# convert region id to region
sub regionid2region{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $sam_id = $parm{sam_id};
  my $regidstr = $parm{regidstr};
  my $subregion_Aref = $parm{subreg_Aref};

  # variants
  my $Sample_Href = $V_Href->{Sample_Href};
  my $Segment_Href = $V_Href->{Segment_Href};

  # sample object
  my $sam_OB = $Sample_Href->{$sam_id};

  # region id => region
  for my $reg_id (split /,/, $regidstr){
    # region orientation in local haplotype
    my $orientation = 1; # forward
    if($reg_id =~ /^r_/){ # reverse signal
      $reg_id =~ s/^r_//;
      $orientation = -1;
    }
    # region type, host or virus?
    # new name: e.g., H1, V1
    my $regtype = ($reg_id =~ /^H/)?'host':'virus'; # OLD: /^[A-Z]$/
    my ($refseg, $lf_pos, $rt_pos) = $sam_OB->get_segid_details( segtype => $regtype, segid => $reg_id );
    my ($st_pos,$ed_pos) = ($orientation==1)?($lf_pos,$rt_pos):($rt_pos,$lf_pos);
    if($regtype =~ /host/){
      # judge which host segment
      my $hostseg_NO = FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->which_host_segment(
                refseg=>$refseg, lf_pos=>$lf_pos, rt_pos=>$rt_pos
               );
      # check
      if($hostseg_NO == -1){
        warn_and_exit "Cannot determine which host segment the $sam_id $reg_id belongs to.\n";
      }
      # store
      push @$subregion_Aref, [$hostseg_NO, $refseg, $st_pos, $ed_pos, $orientation];
      ## allow the host seg to draw
      my $seg_OB = $Segment_Href->{$hostseg_NO};
      $seg_OB->enable_draw();
    }
    else{ # use 0 as virus seg_NO
      push @$subregion_Aref, [0, $refseg, $st_pos, $ed_pos, $orientation];
    }
  }
}

#--- draw bezier line between two segments
sub bezier_between_segments{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $seg1_OB = $parm{seg1_OB};
  my $seg1_rad = $parm{seg1_rad};
  my $seg2_OB = $parm{seg2_OB};
  my $seg2_rad = $parm{seg2_rad};
  my $color = $parm{color};

  # variants
  my $circos = $V_Href->{circos_OB};

  # Q line
  my ($x1, $y1) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$seg1_rad, radius=>$seg1_OB->{inner_radius} );
  my ($x2, $y2) = get_coordinate_on_circle( cx=>$V_Href->{cx}, cy=>$V_Href->{cy}, rad=>$seg2_rad, radius=>$seg2_OB->{inner_radius} );
  my $path = "M$x1,$y1,Q$V_Href->{cx},$V_Href->{cy},$x2,$y2";
  # bezier line
  $circos->path(
          d=>$path,
          stroke=>$color,
          'stroke-width'=>$V_Href->{path_width},
          'stroke-linecap'=>'round',
          'stroke-dasharray'=>' ',
          fill=>"none",
          opacity=>$V_Href->{path_line_opacity}
         );
}

#--- adjust font_size of text in PNC
sub PNC_text_size{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $text_number = $parm{PNC_text};

  # variants

  my $text_font_size;
  if( !$V_Href->{show_arc_reg} ){ # path number, old_exp: $text_number =~ /^\d+$/
    $text_font_size = (length($text_number)==1) ? $V_Href->{PNC_font_size} : ($V_Href->{PNC_font_size}-1.5);
    $text_font_size -= 1 if($text_number =~ /^[\%]$/);
    $text_font_size += 1 if($text_number =~ /^[\*]$/);
  }
  else{ # should be showing the region for adjustment
    $text_font_size = 3;
  }

  return $text_font_size;
}

#--- 
1; ## tell the perl script the successful access of this module.
