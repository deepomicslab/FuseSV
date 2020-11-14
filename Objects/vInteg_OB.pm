package FuseSV::Objects::vInteg_OB;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit /;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::MicroHomology::SeqFindMH qw/ Find_vtMH Find_sfMH /;
use BioFuse::Visual::SVG_Util::Font qw / show_text_in_line /;
use BioFuse::Visual::SVG_Util::RectSysEle qw/ draw_a_parallelogram draw_a_ellipse /;
use BioFuse::Visual::ArrangeObj qw/ allocate_object_vertical /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'FuseSV::Objects::vInteg_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-08-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_refseg
                        get_flank_region
                        get_vInteg_info
                        get_vInteg_info_Aref
                        get_MH_info
                        get_MH_iShift
                        load_flank_seq
                        has_complete_flkseq
                        has_inner_ins
                        find_flkseq_MH
                        get_MH_landscape
                        draw_MH
                     /;

#--- structure of object
# vInteg_OB -> refseg = { p5=>$refseg_p5, p3=>$refseg_p3 }
# vInteg_OB -> bkpos  = { p5=>$bkpos_p5,  p3=>$bkpos_p3 }
# vInteg_OB -> strand = { p5=>$strand_p5, p3=>$strand_p3 }
# vInteg_OB -> inner_ins = $inner_ins_seq / 'NONE'
# vInteg_OB -> vtMH -> $idx => { len=>, seq=>, extSign=>'-1/0/1' }
# vInteg_OB -> sfMH -> $idx => { len=>, seq=>, sft=>, opl=>, extSign=>'-5/-3/0/3/5' }
# vInteg_OB -> flkseq -> { p5=>{junc=> mkup=>}, p3=>{} }
# vInteg_OB -> flklen -> { p5=>{junc=> mkup=>}, p3=>{} }
# vInteg_OB -> {sr_count} = {T=>, N=>}
# vInteg_OB -> {jr_count} = {T=>, N=>}

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;

    my $vInteg_OB = {};
    $vInteg_OB->{refseg}->{p5} = $parm{refseg_p5};
    $vInteg_OB->{bkpos} ->{p5} = $parm{bkpos_p5};
    $vInteg_OB->{strand}->{p5} = $parm{strand_p5};
    $vInteg_OB->{refseg}->{p3} = $parm{refseg_p3};
    $vInteg_OB->{bkpos} ->{p3} = $parm{bkpos_p3};
    $vInteg_OB->{strand}->{p3} = $parm{strand_p3};
    $vInteg_OB->{inner_ins} = $parm{inner_ins} || 'NONE';
    # supporting reads, load from virus integratin result
    $vInteg_OB->{sr_count}->{T} = $parm{t_uniq_sr} if defined $parm{t_uniq_sr};
    $vInteg_OB->{sr_count}->{N} = $parm{n_uniq_sr} if defined $parm{n_uniq_sr};
    $vInteg_OB->{jr_count}->{T} = $parm{t_uniq_jr} if defined $parm{t_uniq_jr};
    $vInteg_OB->{jr_count}->{N} = $parm{n_uniq_jr} if defined $parm{n_uniq_jr};

    bless($vInteg_OB);
    return $vInteg_OB;
}

#--- return refseg of given prime ---
sub get_refseg{
    my $vInteg_OB = shift;
    my %parm = @_;
    return $vInteg_OB->{refseg}->{ $parm{p} };
}

#--- return host bk ---
sub get_host_bk{
    my $vInteg_OB = shift;
    my %parm = @_;
    my $hostRefHref = $parm{hostRefHref};

    for my $prime (qw/ p5 p3 /){
        next unless exists $hostRefHref->{ $vInteg_OB->{refseg}->{$prime} };
        return [ $vInteg_OB->{refseg}->{$prime},
                 $vInteg_OB->{bkpos}->{$prime},
                 $vInteg_OB->{strand}->{$prime}  ];
    }
    # alert
    warn_and_exit "<ERROR>\tcannot find host bk of virusInteg OB.\n".Dumper($vInteg_OB);
}

#--- return flank region ---
sub get_flank_region{
    my $vInteg_OB = shift;
    my %parm = @_;
    my $prime = $parm{p};
    my $flank = $parm{flank};
    my $type  = $parm{type};

    my $Flk_Href = {};

    # interval
    my $bkpos = $vInteg_OB->{bkpos}->{$prime};
    my $strd  = $vInteg_OB->{strand}->{$prime};
    if( $prime eq 'p5' ){
        if( $type eq 'junc' ){
            $Flk_Href->{stp} = ( $strd eq '+' ? ($bkpos-$flank+1) : $bkpos );
            $Flk_Href->{edp} = ( $strd eq '+' ? $bkpos : ($bkpos+$flank-1) );
        }
        else{
            $Flk_Href->{stp} = ( $strd eq '+' ? ($bkpos+1) : ($bkpos-$flank) );
            $Flk_Href->{edp} = ( $strd eq '+' ? ($bkpos+$flank) : ($bkpos-1) );
        }
    }
    elsif( $prime eq 'p3' ){
        if( $type eq 'junc' ){
            $Flk_Href->{stp} = ( $strd eq '+' ? $bkpos : ($bkpos-$flank+1) );
            $Flk_Href->{edp} = ( $strd eq '+' ? ($bkpos+$flank-1) : $bkpos );
        }
        else{
            $Flk_Href->{stp} = ( $strd eq '+' ? ($bkpos-$flank) : ($bkpos+1) );
            $Flk_Href->{edp} = ( $strd eq '+' ? ($bkpos-1) : ($bkpos+$flank) );
        }
    }
    # reverse-compl
    $Flk_Href->{rvp} = ( $strd eq '+' ? 0 : 1 );

    return $Flk_Href;
}

#--- return vInteg basic info ---
sub get_vInteg_info{
    my $vInteg_OB = shift;
    return join( ';', @{ $vInteg_OB->get_vInteg_info_Aref } );
}

#--- return ref of array storing vInteg basic info  ---
sub get_vInteg_info_Aref{
    my $vInteg_OB = shift;
    return  [   'REFSEG_P5='.$vInteg_OB->{refseg}->{p5},
                'BKPOS_P5=' .$vInteg_OB->{bkpos} ->{p5},
                'STRAND_P5='.$vInteg_OB->{strand}->{p5},
                'REFSEG_P3='.$vInteg_OB->{refseg}->{p3},
                'BKPOS_P3=' .$vInteg_OB->{bkpos} ->{p3},
                'STRAND_P3='.$vInteg_OB->{strand}->{p3},
                'INNER_INS='.$vInteg_OB->{inner_ins}
            ];
}

#--- return info of given type MH ---
sub get_MH_info{
    my $vInteg_OB = shift;
    my %parm = @_;
    my $MHtype = $parm{type};

    my @MH_info;
    my $MH_info_Href = $vInteg_OB->{$MHtype};
    my $iShift = $vInteg_OB->get_MH_iShift;
    for my $i ( sort {$a<=>$b} keys %$MH_info_Href ){
        if( $MHtype eq 'vtMH' ){
            push @MH_info, join( ';',
                                      'IDX='.($i-$iShift),
                                      'LEN='.$MH_info_Href->{$i}->{len},
                                      'SEQ='.$MH_info_Href->{$i}->{seq},
                                 'EXT_SIGN='.$MH_info_Href->{$i}->{extSign}
                                );
        }
        elsif( $MHtype eq 'sfMH' ){
            push @MH_info, join( ';',
                                   'IDX_P5='.($i-$iShift),
                                 'SHIFT_P3='.$MH_info_Href->{$i}->{sft},
                                      'LEN='.$MH_info_Href->{$i}->{len},
                                      'SEQ='.$MH_info_Href->{$i}->{seq},
                                    'OPLEN='.$MH_info_Href->{$i}->{opl},
                                 'EXT_SIGN='.$MH_info_Href->{$i}->{extSign}
                                );
        }
    }

    return \@MH_info;
}

#--- return MH i shift relative to zero ---
# p5 junc pos is zero idx
sub get_MH_iShift{
    my $vInteg_OB = shift;
    return length($vInteg_OB->{flkseq}->{p5}->{junc}) - 1;
}

#--- load flank sequence ---
sub load_flank_seq{
    my $vInteg_OB = shift;
    my %parm = @_;
    $vInteg_OB->{flkseq}->{ $parm{p} }->{ $parm{type} } = $parm{seq};
    $vInteg_OB->{flklen}->{ $parm{p} }->{ $parm{type} } = length($parm{seq});
}

#--- test whether has completed flkseq ---
sub has_complete_flkseq{
    my $vInteg_OB = shift;
    for my $prime ( 'p5', 'p3' ){
        for my $flkType ( 'junc', 'mkup' ){
            unless(exists $vInteg_OB->{flkseq}->{$prime}->{$flkType}){
                return 0; # lack
            }
        }
    }
    # completed
    return 1;
}

#--- test whether has inner-ins ---
sub has_inner_ins{
    my $vInteg_OB = shift;
    return $vInteg_OB->{inner_ins} !~ /NONE/i;
}

#--- find MH from flk seg ---
sub find_flkseq_MH{
    my $vInteg_OB = shift;
    my %parm = @_;
    my $MH_type = $parm{MH_type};
    my $min_MHlen = $parm{min_MHlen} || 2;
    my $max_MHlen = $parm{max_MHlen} || 1000;
    my $min_sfMHoplen = $parm{min_sfMHoplen} || 0;
    my $max_sfMHshift = $parm{max_sfMHshift} || 1000;

    my $flkseq_p5 = $vInteg_OB->{flkseq}->{p5}->{junc} . $vInteg_OB->{flkseq}->{p5}->{mkup};
    my $flkseq_p3 = $vInteg_OB->{flkseq}->{p3}->{mkup} . $vInteg_OB->{flkseq}->{p3}->{junc};

    my @parm = ( seq_a=>$flkseq_p5, seq_b=>$flkseq_p3, min_MHlen=>$min_MHlen, max_MHlen=>$max_MHlen );
    if( $MH_type eq 'vtMH' ){ # vertical MH
        $vInteg_OB->{vtMH} = Find_vtMH( @parm );
    }
    else{ # shift MH
        $vInteg_OB->{sfMH} = Find_sfMH(
                                @parm,
                                vtMH_Href => $vInteg_OB->{vtMH},
                                max_shift => $max_sfMHshift,
                                min_MHopl => $min_sfMHoplen
                            );
    }
}

#--- return MH landscape text ---
sub get_MH_landscape{
    my $vInteg_OB = shift;

    my $iShift = $vInteg_OB->get_MH_iShift;
    my @text;

    # index zero and masked axis ('x')
    push @text, "#IDX0:\t".(' ' x $iShift).'$~';
    # prepare inner part, lines and inner-ins
    my $inner_part = '|';
    my $inner_part_blk = '|';
    if( $vInteg_OB->has_inner_ins ){
        $inner_part .= $vInteg_OB->{inner_ins} . '|';
        $inner_part_blk .= ( ' ' x length($vInteg_OB->{inner_ins}) ) . '|';
        $text[0] .= ( '~' x length($vInteg_OB->{inner_ins}) ) . '~';
    }
    # sfMH
    my %sfMH;
    for my $p ( qw/ p5 p3 / ){
        @{$sfMH{$p}} = map {'.'} (1 .. $vInteg_OB->{flklen}->{$p}->{junc} + $vInteg_OB->{flklen}->{$p}->{mkup});
    }
    for my $i (sort {$a<=>$b} keys %{$vInteg_OB->{sfMH}}){
        my $len = $vInteg_OB->{sfMH}->{$i}->{len};
        my $sft = $vInteg_OB->{sfMH}->{$i}->{sft};
        $sfMH{p5}->[$_] = '*' for ( $i .. $i+$len-1 );
        $sfMH{p3}->[$_] = '*' for ( $i+$sft .. $i+$sft+$len-1 );
    }
    # sfMH p5
    push @text, "#sfP5:\t".
                join( '', @{$sfMH{p5}}[ 0 .. $iShift ] ).
                $inner_part_blk.
                join( '', @{$sfMH{p5}}[ $iShift+1 .. scalar(@{$sfMH{p5}})-1 ] );
    # Seq p5
    push @text, "#SQP5:\t".
                uc( $vInteg_OB->{flkseq}->{p5}->{junc} ).
                $inner_part.
                lc( $vInteg_OB->{flkseq}->{p5}->{mkup} );
    # vtMH
    my @vtMH = map {'.'} (1 .. $vInteg_OB->{flklen}->{p5}->{junc} + $vInteg_OB->{flklen}->{p5}->{mkup});
    for my $i (sort {$a<=>$b} keys %{$vInteg_OB->{vtMH}}){
        $vtMH[$_] = '+' for ( $i .. $i+$vInteg_OB->{vtMH}->{$i}->{len}-1 );
    }
    push @text, "#vtMH:\t".
                join( '', @vtMH[ 0 .. $iShift ] ).
                $inner_part_blk.
                join( '', @vtMH[ $iShift+1 .. $#vtMH ] );
    # Seq p3
    push @text, "#SQP3:\t".
                lc( $vInteg_OB->{flkseq}->{p3}->{mkup} ).
                $inner_part.
                uc( $vInteg_OB->{flkseq}->{p3}->{junc} );
    # sfMH p3
    push @text, "#sfP3:\t".
                join( '', @{$sfMH{p3}}[ 0 .. $iShift ] ).
                $inner_part_blk.
                join( '', @{$sfMH{p3}}[ $iShift+1 .. scalar(@{$sfMH{p3}})-1 ] );

    return \@text;
}

#--- show MH SVG ---
sub draw_MH{
    my $vInteg_OB = shift;
    my @font_parm = (
        font_family => $V_Href->{general_FontFam},
        font_size => $V_Href->{general_FontSize}
    );
    # vInteg info
    for my $info ( @{ $vInteg_OB->get_vInteg_info_Aref } ){
        my $sizeA = show_text_in_line(
                        \$V_Href->{general_SVG_obj},
                        text_x => $V_Href->{general_BaseX},
                        text_y => $V_Href->{general_BaseY},
                        text => $info,
                        text_col => 'black',
                        text_anchor => 'start',
                        @font_parm
                    );
        $V_Href->{general_BaseY} += $sizeA->[0] + $V_Href->{SVG_raw_gap};
    }
    # vInteg Junc Seq and MH
    my %SQ = ( P5=>{}, JC=>{}, P3=>{} );
    my %SQno = ( P5=>1, JC=>2, P3=>3 );
    my $iShift = $vInteg_OB->get_MH_iShift;
    my $innerIns_len = $vInteg_OB->has_inner_ins ? length($vInteg_OB->{inner_ins}) : 0;
    my %comb = ( P5=>'p5,junc,p5,mkup', JC=>'p5,junc,p3,junc', P3=>'p3,mkup,p3,junc' );
    for my $type ( keys %SQ ){
        my @comb = split /,/, $comb{$type};
        my @base = split //, $vInteg_OB->{flkseq}->{$comb[0]}->{$comb[1]} . $vInteg_OB->{flkseq}->{$comb[2]}->{$comb[3]};
        for my $i ( 0 .. $#base ){
            my $x_t = $i > $iShift ? ($i+$innerIns_len+1+$V_Href->{SVG_JuncSiteGapRatio}) : ($i+1);
            my $bgcol = $V_Href->{SVG_bgCol}->{ ( $i > $iShift ? $comb[2] : $comb[0] ) };
            $SQ{$type}{$i} = { base=>$base[$i], bgcol=>$bgcol, X_t=>$x_t };
        }
    }
    ## prepare sfMH
    my %sfMHloc = ( P5=>{}, P3=>{} );
    my $sfMHcount = 0;
    for my $i (sort {$a<=>$b} keys %{$vInteg_OB->{sfMH}}){
        $sfMHcount++;
        my $len = $vInteg_OB->{sfMH}->{$i}->{len};
        my $sft = $vInteg_OB->{sfMH}->{$i}->{sft};
        for my $type ( qw/ P5 P3 / ){
            my $sft = $type eq 'P5' ? 0 : $sft;
            my $i_X_t = $SQ{$type}{$i+$sft}{X_t};
            my $j_X_t = $SQ{$type}{$i+$sft+$len-1}{X_t};
            $sfMHloc{$type}{$sfMHcount} = { layer=>undef, X_t_itv=>[[$i_X_t*$V_Href->{SVG_X_t_toInt}, $j_X_t*$V_Href->{SVG_X_t_toInt}]] };
        }
    }
    ## arrange sfMH line
    my %maxLayer = ( P5=>undef, P3=>undef );
    if( $sfMHcount != 0 ){
        for my $type ( qw/ P5 P3 / ){
            allocate_object_vertical(
                obj_Href => $sfMHloc{$type},
                obj_itval_key => 'X_t_itv',
                obj_layer_key => 'layer',
                obj_gap => 0.2*$V_Href->{SVG_X_t_toInt},
                max_layer_Sref => \$maxLayer{$type}
            );
        }
        ## sfMH (P5)
        $vInteg_OB->show_sfMH( sfMHloc_Href => \%sfMHloc, p => 'P5', max_layer => $maxLayer{P5} );
    }
    ## vtMH
    for my $i (sort {$a<=>$b} keys %{$vInteg_OB->{vtMH}}){
        my $len = $vInteg_OB->{vtMH}->{$i}->{len};
        my $col = $V_Href->{SVG_bgCol}->{ ( $len == 1 ? 'sMH' : 'lMH' ) };
        for my $j ( $i .. $i+$len-1 ){
            $_->{$j}->{bgcol} = $col for values %SQ;
        }
    }
    ## sequences
    my $SeqRawHeight = $V_Href->{SVG_SeqBaseRectHeight} + $V_Href->{SVG_raw_gap};
    $V_Href->{general_BaseY} -= $SeqRawHeight;
    for my $type ( keys %SQ ){
        my $NO = $SQno{$type};
        my $Y_shift = $NO * $SeqRawHeight;
        # seq notes
        my $note = $type eq 'JC' ? 'JUNC:' : $vInteg_OB->{refseg}->{lc($type)}.':';
        show_text_in_line(
            \$V_Href->{general_SVG_obj},
            text_x => $V_Href->{general_BaseX},
            text_y => $V_Href->{general_BaseY} + $Y_shift,
            text => $note,
            text_col => 'black',
            text_anchor => 'end',
            height_adjust => 0.5,
            @font_parm
        );
        # shared parm
        my @shared_parm = (
            head_bottom_side_len => $V_Href->{SVG_SeqBaseRectWidth},
            left_right_side_len => $V_Href->{SVG_SeqBaseRectHeight},
            @font_parm
        );
        # reference bases
        for my $i ( sort {$a<=>$b} keys %{$SQ{$type}} ){
            draw_a_parallelogram(
                \$V_Href->{general_SVG_obj},
                x => $V_Href->{general_BaseX} + $SQ{$type}{$i}{X_t} * $V_Href->{SVG_SeqBaseRectWidth},
                y => $V_Href->{general_BaseY} + $Y_shift,
                fill_color => $SQ{$type}{$i}{bgcol},
                text_fill => $SQ{$type}{$i}{base},
                boundary_color => 'none',
                @shared_parm
            );
        }
        # inner ins
        if( $vInteg_OB->has_inner_ins ){
            my @base = split //, ( $type eq 'JC' ? lc($vInteg_OB->{inner_ins}) : ('*' x $innerIns_len) );
            my $baseXmoveRatio = $V_Href->{SVG_JuncSiteGapRatio} / 2;
            for my $i ( 0 .. $#base ){
                draw_a_parallelogram(
                    \$V_Href->{general_SVG_obj},
                    x => $V_Href->{general_BaseX} + ($iShift+1+$i+1+$baseXmoveRatio) * $V_Href->{SVG_SeqBaseRectWidth},
                    y => $V_Href->{general_BaseY} + $Y_shift,
                    fill_color => 'none',
                    text_fill => $base[$i],
                    boundary_color => ( $type eq 'JC' ? 'black' : 'none' ),
                    @shared_parm
                );
            }
        }
    }
    ## inner ins line
    my $InsLineWidth = $V_Href->{SVG_SeqBaseRectWidth} * $V_Href->{SVG_JuncSiteGapRatio} / 2;
    my @InsLineXmove = ( ($SQ{JC}{$iShift}{X_t}+0.5) * $V_Href->{SVG_SeqBaseRectWidth} + $InsLineWidth );
    if( $innerIns_len != 0 ){ # right-hand side line, if has
        push @InsLineXmove, ($SQ{JC}{$iShift+1}{X_t}-0.5) * $V_Href->{SVG_SeqBaseRectWidth} - $InsLineWidth;
    }
    for my $InsLineXmove ( @InsLineXmove ){
        $V_Href->{general_SVG_obj}->line(
            x1 => $V_Href->{general_BaseX} + $InsLineXmove,
            y1 => $V_Href->{general_BaseY} + 1 * $SeqRawHeight - $V_Href->{SVG_SeqBaseRectHeight} / 2,
            x2 => $V_Href->{general_BaseX} + $InsLineXmove,
            y2 => $V_Href->{general_BaseY} + 3 * $SeqRawHeight + $V_Href->{SVG_SeqBaseRectHeight} / 2,
            stroke => 'black',
            "stroke-width" => $InsLineWidth,
            opacity => 1,
            'stroke-linecap' => 'round',
            'stroke-dasharray' => '8,5'
        );
    }
    # ## sfMH black frame
    # for my $i (sort {$a<=>$b} keys %{$vInteg_OB->{sfMH}}){
    #   my $len = $vInteg_OB->{sfMH}->{$i}->{len};
    #   my $sft = $vInteg_OB->{sfMH}->{$i}->{sft};
    #   for my $type ( qw/ P5 P3 / ){
    #       my $NO = $SQno{$type};
    #       my $Y_shift = $NO * $SeqRawHeight;
    #       my $sft = $type eq 'P5' ? 0 : $sft;
    #       my $i_X_t = $SQ{$type}{$i+$sft}{X_t} - 0.5; # left edge
    #       my $j_X_t = $SQ{$type}{$i+$sft+$len-1}{X_t} + 0.5; # right edge
    #       my $X_span = $V_Href->{SVG_SeqBaseRectWidth} * ($j_X_t-$i_X_t);
    #       draw_a_parallelogram(
    #           \$V_Href->{general_SVG_obj},
    #           x => $V_Href->{general_BaseX} + $i_X_t * $V_Href->{SVG_SeqBaseRectWidth} + $X_span / 2,
    #           y => $V_Href->{general_BaseY} + $Y_shift,
    #           fill_color => 'none',
    #           boundary_color => 'black',
    #           head_bottom_side_len => $X_span,
    #           left_right_side_len => $V_Href->{SVG_SeqBaseRectHeight}
    #       );
    #   }
    # }
    # Y moves
    $V_Href->{general_BaseY} += 3 * $SeqRawHeight;
    ## sfMH (P3)
    if( $sfMHcount != 0 ){
        $vInteg_OB->show_sfMH( sfMHloc_Href => \%sfMHloc, p => 'P3', max_layer => $maxLayer{P3} );
    }
}

#--- display SVG sfMH ---
sub show_sfMH{
    my $vInteg_OB = shift;
    my %parm = @_;
    my $sfMHloc = $parm{sfMHloc_Href};
    my $prime = $parm{p};
    my $max_layer = $parm{max_layer};

    my $sfMH_height = $V_Href->{SVG_sfMH_cirNum} ? $V_Href->{SVG_sfMH_radius}*2+1 : $V_Href->{SVG_sfMH_dotR}*2+3;
    # Y moves
    $V_Href->{general_BaseY} += $sfMH_height * ($prime eq 'P5' ? $max_layer : 0.6);
    # sfMH
    my $ratio = $prime eq 'P5' ? -1 : 1;
    for my $sfMH_NO (sort { $prime eq 'P5'
                            ? $sfMHloc->{$prime}{$b}{layer}<=>$sfMHloc->{$prime}{$a}{layer}
                            : $sfMHloc->{$prime}{$a}{layer}<=>$sfMHloc->{$prime}{$b}{layer}
                          }
                     keys %{$sfMHloc->{$prime}})
    {
        my $Y = $V_Href->{general_BaseY} + $ratio * $sfMHloc->{$prime}{$sfMH_NO}{layer} * $sfMH_height;
        my $X_t_itv_Aref = $sfMHloc->{$prime}{$sfMH_NO}{X_t_itv}[0];
        my $st_X = $V_Href->{general_BaseX} + $X_t_itv_Aref->[0]/$V_Href->{SVG_X_t_toInt} * $V_Href->{SVG_SeqBaseRectWidth};
        my $ed_X = $V_Href->{general_BaseX} + $X_t_itv_Aref->[1]/$V_Href->{SVG_X_t_toInt} * $V_Href->{SVG_SeqBaseRectWidth};
        # draw span line
        $V_Href->{general_SVG_obj}->line(
            x1 => $st_X,
            y1 => $Y,
            x2 => $ed_X,
            y2 => $Y,
            stroke => 'black',
            "stroke-width" => 1,
            opacity => 1,
            'stroke-linecap' => 'round'
        );
        # draw circular number
        if( $V_Href->{SVG_sfMH_cirNum} ){
            draw_a_ellipse(
                \$V_Href->{general_SVG_obj},
                cx => $st_X,
                cy => $Y,
                radius => $V_Href->{SVG_sfMH_radius},
                fill_color => 'white',
                text_fill => $sfMH_NO,
                font_size_auto_adjust => 1,
                font_size => 3, # test
                font_family => $V_Href->{general_FontFam}
            );
        }
        else{
            draw_a_ellipse(
                \$V_Href->{general_SVG_obj},
                cx => $st_X,
                cy => $Y,
                radius => $V_Href->{SVG_sfMH_dotR},
                fill_color => 'black'
            );
        }
        # draw end dot
        draw_a_ellipse(
            \$V_Href->{general_SVG_obj},
            cx => $ed_X,
            cy => $Y,
            radius => $V_Href->{SVG_sfMH_dotR},
            fill_color => 'black'
        );
    }
    # Y moves
    $V_Href->{general_BaseY} += $sfMH_height * ($prime eq 'P5' ? 0.6 : $max_layer);
}

1; ## tell the perl script the successful access of this module.
