package FuseSV::Objects::PairEnd_OB;

use strict;
use warnings;
use List::Util qw/ min /;
use Data::Dumper;
use SOAPfuse::General_Operation qw/ warn_and_exit Get_Two_Seg_Olen /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'PairEnd_OB';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-10-17';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        load_reads_OB
                        get_pid
                        get_peIdx
                        get_reads_OB
                        get_sorted_reads_OB
                        discardAbnormalSP
                        printSAM
                     /;

#--- structure of object
# pe_OB -> {reads_OB} -> {1/2} = [ $reads_OB(mapped-1), $reads_OB(mapped-2) ... ]
# pe_OB -> {pe_Idx} = $pe_Idx

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $pe_Idx = $parm{pe_Idx} || 0; # this is for short-name of the pe-id

    my $pe_OB = {};
    $pe_OB->{pe_Idx} = $pe_Idx;

    bless($pe_OB);
    return $pe_OB;
}

#--- load reads OB ---
sub load_reads_OB{

    my $pe_OB = shift;
    my %parm = @_;
    my $reads_OB = $parm{reads_OB};

    my $reads_end = $reads_OB->get_endNO;
    push @{$pe_OB->{reads_OB}->{$reads_end}}, $reads_OB;
}

#--- return pid ---
sub get_pid{
    my $pe_OB = shift;
    return $pe_OB->{reads_OB}->{1}->[0]->get_pid;
}

#--- return pe-idx ---
sub get_peIdx{
    my $pe_OB = shift;
    return $pe_OB->{pe_Idx};
}

#--- return array-ref of given reads-end ---
sub get_reads_OB{
    my $pe_OB = shift;
    my %parm = @_;
    my $reads_end = $parm{reads_end};

    return $pe_OB->{reads_OB}->{$reads_end};
}

#--- return sorted reads OB array-ref ---
## ascending sort chr,pos
sub get_sorted_reads_OB{
    my $pe_OB = shift;
    my %parm = @_;
    my $rEndAref = $parm{rEndAref}; # [1,2] or [1] or [2]

    my @reads_OB = map { ( @{ $pe_OB->get_reads_OB( reads_end => $_ ) } ) } @$rEndAref;
    # chr,pos ascending sort
    @reads_OB = sort
                {   $a->get_mseg cmp $b->get_mseg
                 || $a->get_mpos <=> $b->get_mpos
                } @reads_OB;
    # return Aref
    return \@reads_OB;
}

#--- discard abnormal supplementary alignment ---
sub discardAbnormalSP{
    my $pe_OB = shift;

    my $removeSPbool = 0;
    for my $rEnd (1,2){
        my $rOB_Aref = $pe_OB->get_reads_OB( reads_end => $rEnd );
        # only deal with reads-end has two alignments
        next if( scalar(@$rOB_Aref) != 2 );
        # which is SP (supplementary) and FA (first-align)
        my $SP_i = ( $rOB_Aref->[0]->is_suppmap ? 0 : 1 );
        my $FA_i = ( $SP_i + 1 ) % 2;
        my $SP_rOB = $rOB_Aref->[$SP_i];
        my $FA_rOB = $rOB_Aref->[$FA_i];
        if( $FA_rOB->is_suppmap ){
            warn_and_exit "<ERROR>\tWrong PE-r$rEnd having two supplementary alignments.\n"
                                ."\t".Dumper($pe_OB)."\n";
        }
        #------------------------#
        # filter sum of mReadLen #
        #------------------------#
        my $mReadLenSum = $SP_rOB->get_mReadLen + $FA_rOB->get_mReadLen;
        my $origReadlen = $FA_rOB->get_rlen;
        if(    $mReadLenSum < $origReadlen * 0.8
            || $mReadLenSum > $origReadlen * 1.2
        ){
            splice( @$rOB_Aref, $SP_i, 1 );
            $removeSPbool = 1;
            next;
        }
        #-------------------------------#
        # filter overlap of mapped part #
        #-------------------------------#
        my $SP_fClipL = $SP_rOB->get_foreClipLen( clipType => 'SH' );
        my $SP_hClipL = $SP_rOB->get_hindClipLen( clipType => 'SH' );
        my $FA_fClipL = $FA_rOB->get_foreClipLen( clipType => 'SH' );
        my $FA_hClipL = $FA_rOB->get_hindClipLen( clipType => 'SH' );
        my @SP_mRange = ( $SP_fClipL + 1, $origReadlen - $SP_hClipL );
        my @FA_mRange = ( $FA_fClipL + 1, $origReadlen - $FA_hClipL );
        # alignment orientation diff?
        if(    ( $FA_rOB->is_fw_map && $SP_rOB->is_rv_map )
            || ( $FA_rOB->is_rv_map && $SP_rOB->is_fw_map )
        ){ # to reverse either one
            @SP_mRange = ( $SP_hClipL + 1, $origReadlen - $SP_fClipL );
        }
        # test overlp
        my $overlap_len = Get_Two_Seg_Olen( @FA_mRange[0,1] , @SP_mRange[0,1] );
        my $minMpartLen = min( $SP_rOB->get_mReadLen, $FA_rOB->get_mReadLen );
        if( $overlap_len > $minMpartLen * 0.33 ){
            splice( @$rOB_Aref, $SP_i, 1 );
            $removeSPbool = 1;
            next;
        }
    }

    return $removeSPbool;
}

#--- print PE in SAM format ---
sub printSAM{

    my $pe_OB = shift;
    my %parm = @_;
    my $skip_2ndmap = $parm{skip_2ndmap};
    my $skip_suppmap= $parm{skip_suppmap};
    my $skip_unmap  = $parm{skip_unmap};
    my $skip_dup    = $parm{skip_dup};
    my $skip_mltmap = $parm{skip_mltmap};
    my $skip_SHclip = $parm{skip_SHclip};

    my @SAM;
    for my $reads_end ( sort {$a<=>$b} keys %{$pe_OB->{reads_OB}} ){
        for my $reads_OB ( @{$pe_OB->{reads_OB}->{$reads_end}} ){
            next if( $skip_2ndmap  && $reads_OB->is_2ndmap );
            next if( $skip_suppmap && $reads_OB->is_suppmap );
            next if( $skip_unmap   && $reads_OB->is_unmap );
            next if( $skip_dup     && $reads_OB->is_dup );
            next if( $skip_mltmap  && $reads_OB->is_mltmap );
            next if( $skip_SHclip  && $reads_OB->is_clip );
            push @SAM, $reads_OB->printSAM;
        }
    }

    return \@SAM;
}

1; ## tell the perl script the successful access of this module.


























