package FuseSV::Objects::Virus_OB;

use BioFuse::BioInfo::Objects::Segment::RefSeg_OB; # inheritance

use strict;
use warnings;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter BioFuse::BioInfo::Objects::Segment::RefSeg_OB);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'FuseSV::Objects::Virus_OB';
#----- version --------
$VERSION = "0.02";
$DATE = '2020-08-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        set_genebank
                        genebank
                        set_ncbi_url
                        ncbi_url
                        set_ncbi_txt
                        ncbi_txt
                        set_circular
                        is_circular
                        set_origLen
                        origLen
                        set_extLen
                        extLen
                        set_virus_no
                        no
                        set_ana_ref
                        ana_ref
                        set_init_bam
                        init_bam
                     /;

#--- structure of object
# basis, BioFuse::BioInfo::Objects::Segment::RefSeg_OB
# virus_OB -> circular = 0(no)/1(yes)
# virus_OB -> length = $origLen **basis**
# virus_OB -> extLen = $extLen
# virus_OB -> no = $NO
# virus_OB -> ncbi = {genebank=>, url=>, html=>}
# virus_OB -> ana_ref = $ana_ref
# virus_OB -> init_bam = $init_bam, BioFuse::BioInfo::Objects::SeqData::Bam_OB
# virus_OB -> variant = {}

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;
    my $virus_OB = BioFuse::BioInfo::Objects::Segment::RefSeg_OB->new(id=>$parm{id});
    bless($virus_OB);
    return $virus_OB;
}

#--- set genebank id ---
sub set_genebank{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{ncbi}->{genebank}  = $parm{genebank};
}

#--- return ncbi url ---
sub genebank{
    my $virus_OB = shift;
    return $virus_OB->{ncbi}->{genebank};
}

#--- set ncbi url ---
sub set_ncbi_url{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{ncbi}->{url}  = $parm{url};
}

#--- return ncbi url ---
sub ncbi_url{
    my $virus_OB = shift;
    return $virus_OB->{ncbi}->{url};
}

#--- set txt of ncbi url ---
sub set_ncbi_txt{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{ncbi}->{txt}  = $parm{txt};
}

#--- return txt of ncbi url ---
sub ncbi_txt{
    my $virus_OB = shift;
    return $virus_OB->{ncbi}->{txt};
}

#--- set circular attribute ---
sub set_circular{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{circular}  = $parm{circular};
}

#--- check whether circular or not ---
sub is_circular{
    my $virus_OB = shift;
    return ($virus_OB->{circular} eq 'yes' || $virus_OB->{circular});
}

#--- set extended length specific for circular virus ---
sub set_extLen{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{extLen}  = $parm{extLen};
}

#--- return extended length ---
sub extLen{
    my $virus_OB = shift;
    return $virus_OB->{extLen};
}

#--- set virus NO ---
sub set_virus_no{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{no}  = $parm{no};
}

#--- return virus NO ---
sub no{
    my $virus_OB = shift;
    return $virus_OB->{no};
}

#--- set analysis fasta ref ---
sub set_ana_ref{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{ana_ref}  = $parm{fa};
}

#--- return analysis fasta ref ---
sub ana_ref{
    my $virus_OB = shift;
    return $virus_OB->{ana_ref};
}

#--- set initial aligned bam ---
sub set_init_bam{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{init_bam}  = $parm{bam};
}

#--- return initial aligned bam ---
sub init_bam{
    my $virus_OB = shift;
    return $virus_OB->{init_bam};
}

1; ## tell the perl script the successful access of this module.
