package FuseSV::Objects::PhasedMut_OB;

use strict;
use warnings;
use SOAPfuse::General_Operation qw/ warn_and_exit Get_Two_Seg_Olen /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'PhasedMut_OB';
#----- version --------
$VERSION = "0.10";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        load_posIdx
                        load_ref_allele
                        load_haplotype_allele
                        load_judgeReadEdgeDist
                        get_infoSummary
                        get_chr
                        get_pos
                        get_refAllele
                        get_alleleTypeToHapAref
                        get_hapInfo
                        get_phMutNO
                        get_judgeReadEdgeDist
                        is_hom
                        has_alleleType
                        has_ref
                        has_snv
                        has_del
                        has_ins
                        has_indel
                        allele2haploID
                        release_memory
                     /;

#--- structure of object
# phasedMut_OB -> chr = $chr
# phasedMut_OB -> pos = $pos
# phasedMut_OB -> posIdx = $posIdx
# phasedMut_OB -> phMutNO = $phMutNO
# phasedMut_OB -> refAllele = $refAllele
# phasedMut_OB -> alleleType = {'ref'=>[h1], 'snv'=>[h2], 'ins'=>[..], 'del'=>[..] }
# phasedMut_OB -> phaAllele = { h1=>{ alleleType=>, alleleSeq=>, delSize(might)=> }, h2=>.. }
# phasedMut_OB -> readEdgeDist = $readEdgeDist (currently, only for InDel alleles)

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $phMutNO = $parm{phMutNO} || undef;

    my $phasedMut_OB = {};

    $phasedMut_OB->{chr} = $chr;
    $phasedMut_OB->{pos} = $pos;
    $phasedMut_OB->{phMutNO} = $phMutNO;
    # allele info
    $phasedMut_OB->{refAllele} = undef;
    $phasedMut_OB->{alleleType} = {};
    $phasedMut_OB->{phaAllele} = {};

    bless($phasedMut_OB);
    return $phasedMut_OB;
}

#--- load pos idx ---
sub load_posIdx{
    my $phasedMut_OB = shift;
    my %parm = @_;
    $phasedMut_OB->{posIdx} = $parm{posIdx};
}

#--- load ref allele ---
sub load_ref_allele{
    my $phasedMut_OB = shift;
    my %parm = @_;
    $phasedMut_OB->{refAllele} = $parm{alleleSeq};
}

#--- load allele to given haplotype ---
sub load_haplotype_allele{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $haplo_ID = $parm{haplo_ID};
    my $alleleType = $parm{alleleType};
    my $alleleSeq = $parm{alleleSeq};

    $phasedMut_OB->{phaAllele}->{$haplo_ID} = { alleleType => $alleleType, alleleSeq => $alleleSeq };
    $phasedMut_OB->{phaAllele}->{$haplo_ID}->{delSize} = length($alleleSeq) if($alleleType eq 'del');
    push @{$phasedMut_OB->{alleleType}->{$alleleType}}, $haplo_ID;

    if( $alleleType eq 'ref' ){
        $phasedMut_OB->{refAllele} = $alleleSeq;
    }
}

#--- load read-edge distance to accept this allele ---
sub load_judgeReadEdgeDist{
    my $phasedMut_OB = shift;
    my %parm = @_;
    $phasedMut_OB->{readEdgeDist} = $parm{readEdgeDist};
}

#--- return info summary ---
sub get_infoSummary{
    my $phasedMut_OB = shift;
    return [
             'NO.'.$phasedMut_OB->get_phMutNO,
             $phasedMut_OB->get_chr,
             $phasedMut_OB->get_pos,
             ($phasedMut_OB->is_hom ? 'Hom' : 'Het'),
             map {
                $_ . '(' .
                join( ',', @{$phasedMut_OB->get_hapInfo(haplo_ID=>$_)},
                         , 'dist:'.( $phasedMut_OB->{readEdgeDist} || 'NULL' )
                ) . ')'
             }
             sort keys %{$phasedMut_OB->{phaAllele}}
           ]
}

#--- return chr ---
sub get_chr{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->{chr};
}

#--- return pos ---
sub get_pos{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->{pos};
}

#--- return refAllele ---
sub get_refAllele{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->{refAllele};
}

#--- return ref of hap-array for given allele type ---
sub get_alleleTypeToHapAref{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $type = $parm{type} || 'ref';

    return $phasedMut_OB->{alleleType}->{$type} || [];
}

#--- return haplotype allele info of given hapID---
sub get_hapInfo{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $haplo_ID = $parm{haplo_ID};

    if( !exists $phasedMut_OB->{phaAllele}->{$haplo_ID} ){
        warn_and_exit "Cannot find allele info of haplo_ID ($haplo_ID).\n";
    }

    return [ $phasedMut_OB->{phaAllele}->{$haplo_ID}->{alleleType},
             $phasedMut_OB->{phaAllele}->{$haplo_ID}->{alleleSeq}  ];
}

#--- return phMutNO ---
sub get_phMutNO{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->{phMutNO};
}

#--- return read-edge distance to accept this allele ---
sub get_judgeReadEdgeDist{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->{readEdgeDist} || 0;
}

#--- test whether is hom ---
sub is_hom{
    my $phasedMut_OB = shift;
    my %hapInfo = map {
                        ( join(',', @{$phasedMut_OB->get_hapInfo(haplo_ID=>$_)}), 1 )
                      }
                  keys %{$phasedMut_OB->{phaAllele}};
    return ( scalar(keys %hapInfo) == 1 );
}

#--- test whether has give alleleType ---
sub has_alleleType{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $typeAref = $parm{typeAref};
    my $oneIsOK = $parm{oneIsOK} || 0;

    for my $type (@$typeAref){
        if( $oneIsOK ){
            return 1 if(  exists $phasedMut_OB->{alleleType}->{$type} );
        }
        else{
            return 0 if( !exists $phasedMut_OB->{alleleType}->{$type} );
        }
    }

    return ($oneIsOK ? 0 : 1);
}

#--- test whether has ref ---
sub has_ref{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->has_alleleType( typeAref=>['ref'] );
}

#--- test whether has snv ---
sub has_snv{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->has_alleleType( typeAref=>['snv'] );
}

#--- test whether has del ---
sub has_del{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->has_alleleType( typeAref=>['del'] );
}

#--- test whether has ins ---
sub has_ins{
    my $phasedMut_OB = shift;
    return $phasedMut_OB->has_alleleType( typeAref=>['ins'] );
}

#--- test whether has ins and/or del ---
sub has_indel{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $hasBoth = $parm{hasBoth} || 0;
    return $phasedMut_OB->has_alleleType( typeAref=>['ins','del'], oneIsOK=>!$hasBoth );
}

#--- find haplo_ID matching given allele ---
sub allele2haploID{
    my $phasedMut_OB = shift;
    my %parm = @_;
    my $allele_OB = $parm{allele_OB};

    my $alleleType = $allele_OB->get_type;
    my $alleleSeq = $allele_OB->get_alleleSeq || '';

    for my $haplo_ID ( keys %{$phasedMut_OB->{phaAllele}} ){
        my $alleleInfo_Href = $phasedMut_OB->{phaAllele}->{$haplo_ID};
        if(    $alleleType =~ /$alleleInfo_Href->{alleleType}/i # note 'refORsnv' match 'ref', 'snv'
            && (   uc($alleleInfo_Href->{alleleSeq}) eq uc($alleleSeq)
                || (   $alleleType eq 'del'
                    && $allele_OB->get_delSize == $alleleInfo_Href->{delSize} # del-size match
                    && $allele_OB->get_deloffset == -1 # the next pos is the 1st-del-pos
                   )
               )
        ){
            return $haplo_ID;
        }
    }

    return undef; # not found
}

#--- delete several attributes to save memory ---
sub release_memory{
    my $phasedMut_OB = shift;
    for my $attrKey (qw/ chr posIdx phMutNO refAllele alleleType phaAllele readEdgeDist /){
        delete $phasedMut_OB->{attrKey};
    }
}

1; ## tell the perl script the successful access of this module.
