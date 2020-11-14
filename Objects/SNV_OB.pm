package FuseSV::Objects::SNV_OB;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'SNV_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-04-21';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_imb
                        ini_hash_Href
                        has_HRC_sample
                     /;

#--- structure of object
# snv_OB -> chr = $chr
# snv_OB -> pos = $pos
# snv_OB -> refA = $refA
# snv_OB -> altA = $altA
# snv_OB -> refC = $refC
# snv_OB -> altC = $altC
# snv_OB -> imb = $imb
# snv_OB -> majorA = $majorA
# snv_OB -> minorA = $minorA
# snv_OB -> HRC = {} # phase info from HRC project
# snv_OB -> 10x = {} # linkage info from 10x-barcode
# snv_OB -> HiC = {} # linkage info from Hi-C reads

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $snv_OB = {};
    $snv_OB->{chr} = $parm{chr};
    $snv_OB->{pos} = $parm{pos};
    $snv_OB->{refA} = $parm{refA}; # ref allele
    $snv_OB->{altA} = $parm{altA}; # alt allele
    $snv_OB->{refC} = $parm{refC}; # ref reads count
    $snv_OB->{altC} = $parm{altC}; # alt reads count
    my $min_imb = $parm{min_imb} || 0.2; # allele imbalance, refC ~ altC
    my $max_imb = $parm{max_imb} || 0.9999; # allele imbalance, refC or altC ~ 0

    # check
    for my $key (qw/ chr pos refA altA refC altC /){
        if( !defined $snv_OB->{$key} ){
            warn_and_exit "<ERROR>\tplease give available variables (chr,pos,refA,altA,refC,altC) to SNV_OB.\n"
                                 ."Your Inputs for $MODULE_NAME->new func:\n"
                                 .Dumper(\%parm)."\n";
        }
    }

    # allele imbalance
    $snv_OB->{imb} = abs($snv_OB->{refC}-$snv_OB->{altC}) / max($snv_OB->{refC},$snv_OB->{altC});
    # major / minor
    if(    $snv_OB->{imb} >= $min_imb
        && $snv_OB->{imb} <= $max_imb
    ){
        ($snv_OB->{majorA}, $snv_OB->{minorA}) =    $snv_OB->{refC} > $snv_OB->{altC}
                                                 ? ($snv_OB->{refA} , $snv_OB->{altA})
                                                 : ($snv_OB->{altA} , $snv_OB->{refA})
    }
    else{
        ($snv_OB->{majorA}, $snv_OB->{minorA}) = (undef, undef);
    }

    bless($snv_OB);
    return $snv_OB;
}


#--- return chr ---
sub get_chr{
    my $snv_OB = shift;
    return $snv_OB->{chr};
}

#--- return pos ---
sub get_pos{
    my $snv_OB = shift;
    return $snv_OB->{pos};
}

#--- return ref allele ---
sub get_refAllele{
    my $snv_OB = shift;
    return $snv_OB->{refA};
}

#--- return alt allele ---
sub get_altAllele{
    my $snv_OB = shift;
    return $snv_OB->{altA};
}

#--- return minor allele ---
sub get_minorAllele{
    my $snv_OB = shift;
    return $snv_OB->{minorA};
}

#--- return major allele ---
sub get_majorAllele{
    my $snv_OB = shift;
    return $snv_OB->{majorA};
}

#--- return imb of this SNV ---
sub get_imb{
    my $snv_OB = shift;
    return $snv_OB->{imb};
}

#--- get chr-pos string ---
sub get_chr_pos{
    my $snv_OB = shift;
    my %parm = @_;
    my $delimit = $parm{delimit};
    my $region = $parm{region} || 0;

    my $str = $snv_OB->{chr}
             .(defined $delimit ? $delimit : ':')
             .$snv_OB->{pos}
             .(  $region
               ? (defined $delimit ? $delimit : '-').$snv_OB->{pos}
               : ''
              );
    return $str;
}

#--- return the L/R phased GT ---
sub get_HRC_sample_Genotp{
    my $snv_OB = shift;
    my %parm = @_;
    my $delimit = $parm{delimit};
    my $sampleID = $parm{sampleID};

    $snv_OB->has_HRC_sample( sampleID => $sampleID, give_err => 1 );

    my @GT =  ( $snv_OB->{HRC}->{$sampleID}->{GT_L},
                $snv_OB->{HRC}->{$sampleID}->{GT_R}  );

    if( !defined $delimit ){
        return @GT;
    }
    else{
        join($delimit, @GT);
    }
}

#--- return the L/R phased allele ---
sub get_HRC_sample_allele{
    my $snv_OB = shift;
    my %parm = @_;
    my $delimit = $parm{delimit};
    my $sampleID = $parm{sampleID};

    $snv_OB->has_HRC_sample( sampleID => $sampleID, give_err => 1 );

    my @alleles = ( $snv_OB->{HRC}->{$sampleID}->{Allele_L},
                    $snv_OB->{HRC}->{$sampleID}->{Allele_R}  );

    if( !defined $delimit ){
        return @alleles;
    }
    else{
        join($delimit, @alleles);
    }
}

#--- initialize linkage hash container ---
sub ini_hash_Href{
    my $snv_OB = shift;
    my %parm = @_;
    my $key_Aref = $parm{key_Aref};

    if(    !defined $key_Aref
        || scalar(@$key_Aref) == 0
    ){
        $key_Aref = ['test'];
    }
    # load
    my $temp_Href = $snv_OB;
    $temp_Href->{$key_Aref->[$_]} = {} for ( 0 .. scalar(@$key_Aref)-1 );
}

#--- check HRC sample existence ---
sub has_HRC_sample{
    my $snv_OB = shift;
    my %parm = @_;
    my $sampleID = $parm{sampleID};
    my $give_err = $parm{give_err} || 0;

    if( exists $snv_OB->{HRC}->{$sampleID} ){
        return 1;
    }
    elsif( $give_err ){
        warn_and_exit "<ERROR>\tCannot find sample ($sampleID) in HRC info of SNV_OB.\n".
                             "\t".$snv_OB->get_chr_pos."\n";
    }
    else{
        return 0;
    }
}

#--- check is het-SNV in HRC sample ---
sub is_het_HRC_sample{
    my $snv_OB = shift;
    my %parm = @_;
    my $sampleID = $parm{sampleID};

    $snv_OB->has_HRC_sample( sampleID => $sampleID, give_err => 1 );

    return $snv_OB->{HRC}->{$sampleID}->{Allele_L} ne $snv_OB->{HRC}->{$sampleID}->{Allele_R};
}

#--- load HRC sample phase info ---
sub load_HRC_sample_info{
    my $snv_OB = shift;
    my %parm = @_;
    my $sampleID = $parm{sampleID};

    for my $key ( qw/ GT_L Allele_L GT_R Allele_R / ){
        # check
        if( !defined $parm{$key} ){
            warn_and_exit "<ERROR>\tplease use available HRC sample variables (GT_L,Allele_L,GT_R,Allele_R).\n"
                                 ."Your Inputs for $MODULE_NAME->load_HRC_sample_info func:\n"
                                 .Dumper(\%parm)."\n";
        }
        # record
        $snv_OB->{HRC}->{$sampleID}->{$key} = $parm{$key};
    }
}

#--- load long-range link info ---
sub load_long_link_info{
    my $snv_OB = shift;
    my %parm = @_;
    my $link_type = $parm{link_type};

    $snv_OB->{$link_type} = $parm{allele2linkmarker};
}

#--- compare link-marker sharing of two SNV_OB ---
sub cmp_otherSNV_link{
    my $snv_OB = shift;
    my %parm = @_;
    my $oth_SNVob = $parm{oth_SNVob};
    my $link_type = $parm{link_type};
    my $cmp_mode = $parm{cmp_mode};

    # attributes
    my $slf_majorA = $snv_OB->get_majorAllele;
    my $slf_minorA = $snv_OB->get_minorAllele;
    my $oth_majorA = $oth_SNVob->get_majorAllele;
    my $oth_minorA = $oth_SNVob->get_minorAllele;

    # compare
    my $shared_linkmarker_count = 0;
    if( $cmp_mode eq 'mSame' ){
        $shared_linkmarker_count += grep exists $oth_SNVob->{$link_type}->{$oth_majorA}->{$_},
                                    keys %{$snv_OB->{$link_type}->{$slf_majorA}};
        $shared_linkmarker_count += grep exists $oth_SNVob->{$link_type}->{$oth_minorA}->{$_},
                                    keys %{$snv_OB->{$link_type}->{$slf_minorA}};
    }
    else{ # 'mDiff'
        $shared_linkmarker_count += grep exists $oth_SNVob->{$link_type}->{$oth_majorA}->{$_},
                                    keys %{$snv_OB->{$link_type}->{$slf_minorA}};
        $shared_linkmarker_count += grep exists $oth_SNVob->{$link_type}->{$oth_minorA}->{$_},
                                    keys %{$snv_OB->{$link_type}->{$slf_majorA}};
    }

    return $shared_linkmarker_count;
}

1; ## tell the perl script the successful access of this module.
