package FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal;

use strict;
use warnings;
use List::Util qw/min max sum/;
use BioFuse::Visual::SVG_Util::RadSys qw/$PI normalize_radian/;
use BioFuse::Visual::SVG_Util::Color qw/%COLOR_DB/;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use BioFuse::Util::Interval qw/Get_Two_Seg_Olen/;
use BioFuse::Util::GZfile qw/Try_GZ_Read/;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Sample_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
        load_sample_seg_info
      /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'Segment_Deal';
#----- version --------
$VERSION = "0.03";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                       load_sample_seg_info
                       intial_segment
                       merge_segment
                       assign_resol_of_subreg
                       create_seg_objects
                       find_radius_for_hostseg
                       which_host_segment
                       load_gene
                       draw_gene
                       load_trans_psl
                       draw_copy_number_band
                       draw_rna_exp_band
                       load_sample_rna_reads_amount
                      /;

#--- load information of samples' segments
## both host segments and virus segments
sub load_sample_seg_info{

  # variants
  my $Sample_SegInfo_file = $V_Href->{Sample_SegInfo_file};
  my $Sample_Href = $V_Href->{Sample_Href};
  my $virus_name = $V_Href->{virus_name};

  my %Seg_TEMP;
  my $Seg_Count = 0;

  my %sam_segid2segNO;
  my ($last_samid, $last_segtype, $last_refseg, $last_RTpos);
  open (SEGIF, $Sample_SegInfo_file) || die "fail open $Sample_SegInfo_file: $!\n";
  # #SampleID Seg_type  Seg_ID  Ref_seg left_pos  right_pos size  res copy_number
  # T001  host  H1  chr12 72579000  72584510  5511  50  2.19
  (my $theme_line = lc(<SEGIF>)) =~ s/^#//;
  my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
  while(<SEGIF>){
    next if(/^\#/);
    my @info = split;
    my %segif_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
    my $samid  = $segif_Ob{sampleid};
    my $segtype = $segif_Ob{seg_type};
    my $ori_segid  = $segif_Ob{seg_id};
    my $refseg = $segif_Ob{ref_seg};
    my $LFpos  = $segif_Ob{left_pos};
    my $RTpos  = $segif_Ob{right_pos};
    my $resol  = $segif_Ob{res};
    my $copy_number = $segif_Ob{copy_number};
    # check valid
    if(    defined $last_samid
      &&    $samid eq $last_samid
         && $segtype eq $last_segtype
         && $refseg eq $last_refseg
         && $LFpos != $last_RTpos + 1
    ){
      warn_and_exit "<ERROR>\tplease give consecutive segments on $refseg of sample $samid:\n"
                ."\tlast_RTpos: $last_RTpos, now_LFpos is $LFpos, not last_RTpos+1\n";
    }
    else{ # update
      $last_samid = $samid;
      $last_segtype = $segtype;
      $last_refseg = $refseg;
      $last_RTpos = $RTpos;
    }
    # special for $segid
    ## sometimes, we may devide one region to sub-regions with different resolutions
    ## e.g., for region H11 harbouring MYC (upstream HPV-loci are far from MYC), sub-region H11_5 and H11_3 use diff-resolutions
    # new name: e.g., H1, V1
    my ($segid) = ($ori_segid =~ /^([HV]\d+)/); # OLD: /^([A-Za-z])/
    # store the segment info in sample object
    if( !exists($Sample_Href->{$samid}) ){
      # create sample object
      $Sample_Href->{$samid} = FuseSV::Virus_Integ::LocalHaplotypeCircos::Sample_OB->new( sam_id => $samid );
    }
    # sample object
    my $sam_OB = $Sample_Href->{$samid};
    # load in the sub-region details
    ## if meet one segid again, merge with old region
    $sam_OB->load_subregion_info(
                    segtype => $segtype,
                    segid => $segid,
                    refseg => $refseg,
                    lfpos => $LFpos,
                    rtpos => $RTpos,
                    copy_number => $copy_number
                  );
    ## record the minimum and maximum copy_number
    if( $segtype eq 'host' ){
      $V_Href->{CN_min} = min( $V_Href->{CN_min}, max($copy_number,1) );
      $V_Href->{CN_max} = max( $V_Href->{CN_max}, $copy_number );
    }
    # load into Seg_TEMP
    if( !exists($sam_segid2segNO{$samid}{$segid}) ){
      $Seg_Count++;
      $Seg_TEMP{$refseg}{$Seg_Count} = [ [$samid], $LFpos, $RTpos, [ [$LFpos, $RTpos, $resol] ] ];
      # record this sample's segid to this segNO
      $sam_segid2segNO{$samid}{$segid} = $Seg_Count;
    }
    else{ # sometimes, we may devide one region to sub-regions with different resolutions
      my $related_Seg_Count = $sam_segid2segNO{$samid}{$segid};
      my $seginfo_Aref = $Seg_TEMP{$refseg}{$related_Seg_Count};
      # update region
      $seginfo_Aref->[1] = min($seginfo_Aref->[1], $LFpos);
      $seginfo_Aref->[2] = max($seginfo_Aref->[2], $RTpos);
      # new resol sub-region
      push @{$seginfo_Aref->[3]}, [$LFpos, $RTpos, $resol];
    }
  }
  close SEGIF;

  # inform
  stout_and_sterr "[INFO]\tload samples' segments info OK.\n";

  # merge segments
  &intial_segment( SegTEMP_Href => \%Seg_TEMP, Seg_Count => $Seg_Count );
}

#--- merge the host segments and deal the resolution of sub-region
sub intial_segment{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $SegTEMP_Href = $parm{SegTEMP_Href};
  my $Seg_Count = $parm{Seg_Count};

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};

  # just merge segments, but leave the resol
  &merge_segment( SegTEMP_Href => $SegTEMP_Href, Seg_Count => $Seg_Count );

  # assign the different resols to different regions in one segment
  &assign_resol_of_subreg( SegTEMP_Href => $SegTEMP_Href );

  # create real segment object for both virus(0) and host(>=1) segments.
  &create_seg_objects( SegTEMP_Href => $SegTEMP_Href );
}

#--- just merge the segments in temp hash
sub merge_segment{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $SegTEMP_Href = $parm{SegTEMP_Href};
  my $Seg_Count = $parm{Seg_Count};

  for my $refseg (sort keys %$SegTEMP_Href){
    my $SegTEMPrefseg_Href = $SegTEMP_Href->{$refseg};
    MERGE_SEG: {
      for my $seg1_NO (sort {$a<=>$b} keys %$SegTEMPrefseg_Href){
        my $samid1_Aref = $SegTEMPrefseg_Href->{$seg1_NO}->[0];
        my $seg1_LFpos = $SegTEMPrefseg_Href->{$seg1_NO}->[1];
        my $seg1_RTpos = $SegTEMPrefseg_Href->{$seg1_NO}->[2];
        my $seg1_resol_AAref = $SegTEMPrefseg_Href->{$seg1_NO}->[3];
        for my $seg2_NO (grep $_ != $seg1_NO , sort {$a<=>$b} keys %$SegTEMPrefseg_Href){
          my $samid2_Aref = $SegTEMPrefseg_Href->{$seg2_NO}->[0];
          my $seg2_LFpos = $SegTEMPrefseg_Href->{$seg2_NO}->[1];
          my $seg2_RTpos = $SegTEMPrefseg_Href->{$seg2_NO}->[2];
          my $seg2_resol_AAref = $SegTEMPrefseg_Href->{$seg2_NO}->[3];
          my $merge_bool =    Get_Two_Seg_Olen($seg1_LFpos, $seg1_RTpos, $seg2_LFpos, $seg2_RTpos)
                   || $seg1_RTpos == $seg2_LFpos - 1
                   || $seg2_RTpos == $seg1_LFpos - 1;
          # print "$refseg:$seg1_LFpos-$seg1_RTpos with $refseg:$seg2_LFpos-$seg2_RTpos\t$merge_bool\n";
          # if( Get_Two_Seg_Olen($seg1_LFpos, $seg1_RTpos, $seg2_LFpos, $seg2_RTpos) ){
          if( $merge_bool ){
            my $new_seg_LFpos = min($seg1_LFpos, $seg1_RTpos, $seg2_LFpos, $seg2_RTpos);
            my $new_seg_RTpos = max($seg1_LFpos, $seg1_RTpos, $seg2_LFpos, $seg2_RTpos);
            # load new hostseg info
            my $new_seg_NO = ++$Seg_Count;
            $SegTEMPrefseg_Href->{$new_seg_NO} = [[], $new_seg_LFpos, $new_seg_RTpos, []];
            ## sample id
            my %hash;
            @{$SegTEMPrefseg_Href->{$new_seg_NO}->[0]} = sort grep {++$hash{$_}<2} (@$samid1_Aref, @$samid2_Aref);
            ## resolutions
            push @{$SegTEMPrefseg_Href->{$new_seg_NO}->[3]}, @$seg1_resol_AAref, @$seg2_resol_AAref;
            # sweep
            delete $SegTEMPrefseg_Href->{$seg1_NO};
            delete $SegTEMPrefseg_Href->{$seg2_NO};
            # inform
            stout_and_sterr "[INFO]\tmerge $refseg $seg1_LFpos-$seg1_RTpos with $seg2_LFpos-$seg2_RTpos as $new_seg_LFpos-$new_seg_RTpos.\n";
            # repeat again
            redo MERGE_SEG;
          }
        }
      }
    }
  }
}

#--- assign the resolutions to sub-regions
sub assign_resol_of_subreg{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $SegTEMP_Href = $parm{SegTEMP_Href};

  for my $refseg (sort keys %$SegTEMP_Href){
    my $SegTEMPrefseg_Href = $SegTEMP_Href->{$refseg};
    for my $seg_NO (sort {$a<=>$b} keys %$SegTEMPrefseg_Href){
      my $seg_resol_AAref = pop @{ $SegTEMPrefseg_Href->{$seg_NO} };
      my %segpos;
      $segpos{$_->[0]} = 1 for @$seg_resol_AAref;
      $segpos{$_->[1]} = 1 for @$seg_resol_AAref;
      my @segpos = sort {$a<=>$b} keys %segpos;
      # add pos to make cling regions
      my @new_segpos = ($segpos[0], $segpos[1], $segpos[1] + 1);
      my $next_min_pos = $new_segpos[-1] + 1; # lf_pos
      for my $i ( 0 .. $#segpos ){
        next if( $segpos[$i] < $next_min_pos );
        push @new_segpos, $segpos[$i]; # rt_pos
         if($i != $#segpos){
          push @new_segpos, $segpos[$i] + 1; # next lf_pos
          $next_min_pos = $new_segpos[-1] + 1; # update
         }
      }
      # check
      if( scalar(@new_segpos) % 2 != 0 ){
        push @new_segpos, $new_segpos[-1];
        stout_and_sterr "<WARN>\tCannot make cling regions on $refseg.\n"
                   ."\torig:\t@segpos\n"
                   ."\tnew:\t@new_segpos\n";
      }
      # replace
      @segpos = @new_segpos;
      # assign resol
      my $curr_min_pos = -1;
      for my $i (0 .. $#segpos){
        my $LFpos = $segpos[$i];
        next if( $LFpos < $curr_min_pos );
        my $RTpos = ( $i == $#segpos ? $segpos[$i] : $segpos[$i+1] );
        $curr_min_pos = $RTpos + 1;
        my @resol;
        for my $Aref (@$seg_resol_AAref){
          my $LFpos2 = $Aref->[0];
          my $RTpos2 = $Aref->[1];
          if($LFpos2 <= $LFpos && $RTpos2 >= $RTpos){
            push @resol, $Aref->[2];
          }
        }
        next if(scalar(@resol) == 0);
        # load the new region resol, always perfer to small resolution
        my $min_resol = min(@resol);
        push @{$SegTEMPrefseg_Href->{$seg_NO}->[3]}, [$LFpos, $RTpos, $min_resol];
        stout_and_sterr "[INFO]\tset segment $refseg:$LFpos-$RTpos resol as $min_resol\n";
      }
    }
  }
}

#--- create real segment objects
sub create_seg_objects{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $SegTEMP_Href = $parm{SegTEMP_Href};

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};
  my $seg_gap_rad = $V_Href->{spacerad};
  my $min_HostSeg_perct = $V_Href->{min_seg_perct};

  # get the linear pixel consumption of each segments
  my %seg_linear_pixel;
  my $full_pixel = 0;
  my $hostseg_count = 0;
  for my $refseg (sort keys %$SegTEMP_Href){
    my $SegTEMPrefseg_Href = $SegTEMP_Href->{$refseg};
    for my $seg_NO (sort {$a<=>$b} keys %$SegTEMPrefseg_Href){
      $hostseg_count++ if( $refseg ne $V_Href->{virus_name} ); # only host segments
      my $seg_resol_AAref = $SegTEMPrefseg_Href->{$seg_NO}->[3];
      for my $Aref (@$seg_resol_AAref){
        my $pixel = ($Aref->[1] - $Aref->[0] + 1) / $Aref->[2];
        $seg_linear_pixel{$refseg}{$seg_NO} += $pixel;
        $full_pixel += $pixel if( $refseg ne $V_Href->{virus_name} ); # only host segments
      }
    }
  }

  # 0 is virus segment, >=1 is host segment.
  my $real_NO = 0;

  #----------------#
  # virus segments #
  #----------------#
  my $virus_seg_rad = 2 * $PI * $V_Href->{virus_perct} / 100;
  my $virus_start_rad = (2 * $PI) * ($V_Href->{virus_middleline_loc_perct} - $V_Href->{virus_perct} / 2) / 100;
  normalize_radian(radian_Sref=>\$virus_start_rad);
  # check virus segment's existence
  unless( exists($SegTEMP_Href->{ $V_Href->{virus_name} }) ){
    warn_and_exit "Cannot find the segment data for virus $V_Href->{virus_name}.\n";
  }
  # virus genome segment should be only entire one
  my @virus_seg_NO = keys %{$seg_linear_pixel{ $V_Href->{virus_name} }};
  if( scalar(@virus_seg_NO) != 1){
    warn_and_exit "Segment count of virus $V_Href->{virus_name} exceeds one.\n";
  }
  # assign rad size for virus segment
  ## create virus segment objects
  $Segment_Href->{$real_NO} = FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->new(
                  NO => $real_NO,
                  sampleid_Aref => $SegTEMP_Href->{$V_Href->{virus_name}}->{$virus_seg_NO[0]}->[0],
                  refseg => $V_Href->{virus_name},
                  lfpos => $SegTEMP_Href->{$V_Href->{virus_name}}->{$virus_seg_NO[0]}->[1],
                  rtpos => $SegTEMP_Href->{$V_Href->{virus_name}}->{$virus_seg_NO[0]}->[2],
                  start_rad => $virus_start_rad,
                  rad_size => $virus_seg_rad,
                  name => $V_Href->{virus_name},
                  inner_radius => $V_Href->{virus_radius},
                  InRadToSegWidRatio => $V_Href->{segment_width_ratio},
                  resol_AAref => $SegTEMP_Href->{$V_Href->{virus_name}}->{$virus_seg_NO[0]}->[3],
                  linear_pixel => $seg_linear_pixel{$V_Href->{virus_name}}->{$virus_seg_NO[0]}
                );

  #---------------#
  # host segments #
  #---------------#
  # rad for host segments to split
  my $FullRad_HostSeg = 2 * $PI - ($hostseg_count + 1) * $seg_gap_rad - $virus_seg_rad;
  # assign rad size for each host segments
  ## create host segment objects
  my $accum_rad = 0;
  for my $chr_icon ( @{$V_Href->{segment_sort_Aref}} ){ # for human being, defaults as 1 .. 22, 'X', "Y"
    my $refseg = "chr$chr_icon";
    next unless(exists($seg_linear_pixel{$refseg}));
    for my $hostseg_NO (sort {$a<=>$b} keys %{$seg_linear_pixel{$refseg}}){
      my $sampleid_Aref = $SegTEMP_Href->{$refseg}->{$hostseg_NO}->[0];
      my $hostseg_LFpos = $SegTEMP_Href->{$refseg}->{$hostseg_NO}->[1];
      my $hostseg_RTpos = $SegTEMP_Href->{$refseg}->{$hostseg_NO}->[2];
      my $resol_AAref   = $SegTEMP_Href->{$refseg}->{$hostseg_NO}->[3];
      my $seg_linear_pixel = $seg_linear_pixel{$refseg}{$hostseg_NO};
      my $ratio = $seg_linear_pixel / $full_pixel;
      my $rad_size = $FullRad_HostSeg * $ratio;

      # examine the minimum segments percentage
      my $rad_perct = sprintf("%.2f", 100 * $rad_size / (2 * $PI));
      if( $rad_perct <  $min_HostSeg_perct){
        warn_and_exit "Host Segment ($refseg:$hostseg_LFpos-$hostseg_RTpos) section ($rad_perct) is less than minimum ($min_HostSeg_perct).\n";
      }
      else{
        $real_NO++;
        # the start rad
        my $start_rad =   $Segment_Href->{0}->{start_rad}
                + $Segment_Href->{0}->{rad_size} # virus part
                + $real_NO * $seg_gap_rad + $accum_rad;
        # update
        $accum_rad += $rad_size;
        # segment's radius
        my $this_hostseg_radius = &find_radius_for_hostseg( refseg => $refseg, lf_pos => $hostseg_LFpos, rt_pos => $hostseg_RTpos );
        # create Host Segment objects
        $Segment_Href->{$real_NO} = FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_OB->new(
                        NO => $real_NO,
                        sampleid_Aref => $sampleid_Aref,
                        refseg => $refseg,
                        lfpos => $hostseg_LFpos,
                        rtpos => $hostseg_RTpos,
                        start_rad => $start_rad,
                        rad_size => $rad_size,
                        name => "$refseg:$hostseg_LFpos-$hostseg_RTpos",
                        inner_radius => $this_hostseg_radius,
                        InRadToSegWidRatio => $V_Href->{segment_width_ratio},
                        resol_AAref => $resol_AAref,
                        linear_pixel => $seg_linear_pixel
                      );
      }
    }
  }
}

#--- find specific radius for given region
sub find_radius_for_hostseg{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $refseg = $parm{refseg};
  my $lf_pos = $parm{lf_pos};
  my $rt_pos = $parm{rt_pos};

  # variants
  my $spec_radius_AAref = $V_Href->{radius_specific};
  my $min_OP_ratio = $V_Href->{min_spec_radius_overlap_ratio};

  # test overlap with specific region.
  for my $Aref (@$spec_radius_AAref){
    my ($refseg2, $lf_pos2, $rt_pos2, $radius2) = @{$Aref};
    next if( $refseg ne $refseg2 );
    if( Get_Two_Seg_Olen($lf_pos, $rt_pos, $lf_pos2, $rt_pos2) >= $min_OP_ratio * abs($rt_pos - $lf_pos) ){
      stout_and_sterr "[INFO]\tSet radius of host segment ($refseg:$lf_pos-$rt_pos) as $radius2.\n";
      return $radius2;
    }
  }

  return $V_Href->{host_radius};
}

#--- judge which host segment given region belongs to
sub which_host_segment{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $refseg = $parm{refseg};
  my $lf_pos = $parm{lf_pos};
  my $rt_pos = $parm{rt_pos} || $lf_pos; # sometimes, we only want to check one single position

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};

  for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){ # could encounter virus segment
    my $seg_OB = $Segment_Href->{$seg_NO};
    if(    $seg_OB->{ref_seg} eq $refseg
      && $seg_OB->{lf_pos}  <= $lf_pos
      && $seg_OB->{rt_pos}  >= $rt_pos
      ){
      return $seg_NO;
    }
  }
  # invalid return
  return -1;
}

#--- load trans.psl and find genes on segment
sub load_gene{

  # variants
  my $Trans_PSL = $V_Href->{Trans_PSL};
  my $Segment_Href = $V_Href->{Segment_Href};

  # load the transript psl
  my %transcript;
  &load_trans_psl( Trans_PSL => $Trans_PSL,
           Trans_Href => \%transcript,
           Only_001_201 => 1,
           Filter_Href => ( ($V_Href->{allow_gene_str}) ? $V_Href->{allow_gene_Href} : 0 )
          );

  # assign genes to segments
  for my $seg_NO (sort keys %$Segment_Href){
    my $seg_OB = $Segment_Href->{$seg_NO};
    my $seg_refseg = $seg_OB->{ref_seg};
    next unless( exists($transcript{$seg_refseg}) );
    my $seg_lfpos = $seg_OB->{lf_pos};
    my $seg_rtpos = $seg_OB->{rt_pos};
    for my $tname (sort keys %{$transcript{$seg_refseg}}){
      my $trans_Href = $transcript{$seg_refseg}{$tname};
      # check whether overlap with this segment
      unless(   (   $trans_Href->{lf_pos} < $trans_Href->{rt_pos} # usual gene, linear
             && Get_Two_Seg_Olen($trans_Href->{lf_pos}, $trans_Href->{rt_pos}, $seg_lfpos, $seg_rtpos)
            )
           || (   $trans_Href->{lf_pos} > $trans_Href->{rt_pos} # circular gene
             && (   Get_Two_Seg_Olen($trans_Href->{lf_pos}, $seg_rtpos, $seg_lfpos, $seg_rtpos)
                 || Get_Two_Seg_Olen(1, $trans_Href->{rt_pos}, $seg_lfpos, $seg_rtpos)
                )
            )
      ){
        next;
      }
      # store this gene in this seg_OB
      $seg_OB->load_one_transcript( trans_Href => $trans_Href );
    }
  }

  # inform
  stout_and_sterr "[INFO]\tload transcript info from trans.PSL file OK.\n";
}

#--- draw genes on segments
sub draw_gene{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $draw_bool = $parm{draw_bool};

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};

  for my $seg_NO (sort keys %$Segment_Href){
    my $seg_OB = $Segment_Href->{$seg_NO};
    $seg_OB->draw_genes( draw_bool => $draw_bool );
  }
}

#--- load transcript psl to a given hash
sub load_trans_psl{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $Trans_PSL  = $parm{Trans_PSL};
  my $Trans_Href = $parm{Trans_Href};
  my $Filter_Href = $parm{Filter_Href} || 0;
  my $Only_001_201 = $parm{Only_001_201} || 0;

  open (PSL, Try_GZ_Read($Trans_PSL)) || die "fail read Trans_PSL: $!\n";
  while (<PSL>){
    my @ele = (split);
    my $strand = $ele[8];
    my $tname  = $ele[11];
    my $refseg = $ele[13];
    my $lf_pos = $ele[15] + 1;
    my $rt_pos = $ele[16];
    my $ex_len = $ele[18];
    my $ex_lfp = $ele[20];
    my $gname  = $ele[-1];
    # filter genes allowed
    if(    $Filter_Href
      && !exists($Filter_Href->{uc($gname)}) ){
      next;
    }
    # filter transcipt isforms
    if(    $Only_001_201
      && !($tname =~ /-001$/ || $tname =~ /-201$/) # arbitrary, -_-!!!
      ){
      next;
    }
    # loadup this trans-info
    my @ex_len = split /,/, $ex_len;
    my @ex_lfp = split /,/, $ex_lfp;
    $Trans_Href->{$refseg}->{$tname} = {
      refseg => $refseg, strand => $strand,
      gname => $gname, tname => $tname,
      lf_pos => $lf_pos,  rt_pos => $rt_pos,
      ex_len => \@ex_len, ex_lfp => \@ex_lfp
    };
  }
  close PSL;
}

#--- draw graded colored CN band
sub draw_copy_number_band{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $draw_bool = $parm{draw_bool};

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};

  for my $seg_NO (sort keys %$Segment_Href){
    my $seg_OB = $Segment_Href->{$seg_NO};
    $seg_OB->draw_CN_band( draw_bool => $draw_bool );
  }
}

#--- draw graded colored rna exp band
sub draw_rna_exp_band{
  # options
  shift if ($_[0] =~ /::$MODULE_NAME/);
  my %parm = @_;
  my $draw_bool = $parm{draw_bool};

  # variants
  my $Segment_Href = $V_Href->{Segment_Href};

  # gather data, get min_exp, max_exp
  for my $seg_NO (sort keys %$Segment_Href){
    next if( $seg_NO == 0 ); # virus segment
    my $seg_OB = $Segment_Href->{$seg_NO};
    $seg_OB->gather_exp_depth_rad_win;
  }

  # draw rna expression band
  for my $seg_NO (sort keys %$Segment_Href){
    next if( $seg_NO == 0 ); # virus segment
    my $seg_OB = $Segment_Href->{$seg_NO};
    $seg_OB->draw_EXP_band( draw_bool => $draw_bool );
  }
}

#--- load reads amount of each sample
sub load_sample_rna_reads_amount{

  # variants
  my $Sample_Href = $V_Href->{Sample_Href};
  my $Sample_RNA_Reads_Amount_file = $V_Href->{Sample_RNA_Reads_Amount_file};

  open (RRA, Try_GZ_Read($Sample_RNA_Reads_Amount_file)) || die "fail read Sample_RNA_Reads_Amount_file: $!\n";
  (my $theme_line = lc(<RRA>)) =~ s/^#//;
  my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
  while(<RRA>){
    next if(/^\#/);
    my @info = split;
    my %segif_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
    my $samid  = $segif_Ob{sampleid};
    my $reads_amount = $segif_Ob{reads_amount};
    next unless(exists($Sample_Href->{$samid}));
    my $sam_OB = $Sample_Href->{$samid};
    # record to sample Object
    $sam_OB->loadup_rna_reads_amount( reads_amount => $reads_amount );
    # maximum rna reads amount for following normalization
    $V_Href->{RRA_max} = max( $reads_amount, $V_Href->{RRA_max} );
  }
  close RRA;
}

#--- 
1; ## tell the perl script the successful access of this module.
