package FuseSV::Virus_Type::PrepareExtendVirusRef;

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/ dirname /;
use File::Spec qw/ abs2rel catfile rel2abs /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ reset_folder file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file write_fasta_file BWA_index_fasta Faidx_Dict_fasta /;
use FuseSV::LoadOn;
use FuseSV::GetPath;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              PrepareExtendVirusRef
              load_virusBasicInfo
              make_virusAlignRef
              make_virusRefReport
              load_virusRefReport
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Type::PrepareExtendVirusRef';
#----- version --------
$VERSION = "0.06";
$DATE = '2019-06-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        PrepareExtendVirusRef
                        load_virusBasicInfo
                        make_virusAlignRef
                        make_singleVirusRef
                        make_virusRefReport
                        load_virusRefReport
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} extend_vref <[Options]>
     
     Options:

       # Inputs and Outputs #
        -ex     [i]  the length you want to extend the virus genome. <required>
                      You should set it as the read-length in the current analysis.
                      e.g., for WGS with PE90, set this option as 90.
        -vpf    [s]  virus database folder. <required>
                      You can make it via function 'make_virusdb'.
                      OR Download these widely used from FuseSV official website, e.g., HPV.
                      NOTE: the \033[4m \bWRITE \b\033[0m authority is required.

       # Software Required #
        -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
        -stl    [s]  SamTools, minimum version: 1.3. <required>

        -h|help   Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                cluck_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ vref_extLen => 0 ],
            # software
            [ bwa => '' ],
            [ samtools => '' ],

            #--- containers ---#
            # id(subtype) -> virusOB
            # see FuseSV::Objects::Virus_OB
            [ virusRefSeg => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['bwa'],
                                  ['samtools'],
                                  ['virus_db_folder'] ] ]
        );
}

#--- get options from command line
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-ex:i"     => \$V_Href->{vref_extLen},
        "-vpf:s"    => \$V_Href->{virus_db_folder},
        # software
        "-bwa:s"    => \$V_Href->{bwa},
        "-stl:s"    => \$V_Href->{samtools},
        # multiple fork operations
        "-fm:i"     => \$V_Href->{fork_number},
        # help
        "-h|help"   => \$V_Href->{HELP},
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || ($V_Href->{vref_extLen} <= 0)
             || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder} || !-w $V_Href->{virus_db_folder})
             || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
             || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
            );
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub PrepareExtendVirusRef{

    # load virus info from database
    &load_virusBasicInfo;

    # construct virus initial reference genome
    &make_virusAlignRef(
        virus_ref_dir => FuseSV::GetPath->GetPath(filekey => 'extend_vref_db_folder', vref_extLen => $V_Href->{vref_extLen}),
        extend_length => $V_Href->{vref_extLen}
    );
}

#--- load virus information for further operations ---
sub load_virusBasicInfo{
    open (VIF,Try_GZ_Read($V_Href->{virus_info_doc})) || die "fail read $V_Href->{virus_info_doc}: $!\n";
    # theme
    #subtype GeneBank_ID Length Circular Extend_Length
    (my $theme_line = lc(<VIF>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    my $vRefSegNO = 0;
    while(<VIF>){
        next if(/^\#/);
        my @info = split;
        my %vOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        my $virus_OB = FuseSV::Objects::Virus_OB->new(id => $vOb{subtype});
        $V_Href->{virusRefSeg}->{$vOb{subtype}} = $virus_OB;
        $virus_OB->set_genebank(genebank => $vOb{genebank_id});
        $virus_OB->set_circular(circular => $vOb{circular});
        $virus_OB->set_length(length => $vOb{length});
        $virus_OB->set_virus_no(no => ++$vRefSegNO);
    }
    close VIF;
    stout_and_sterr "[INFO]\tload virus basic information from database.\n";
}

#--- make indexed reference genome of all virus in database ---
sub make_virusAlignRef{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $virus_ref_dir = $parm{virus_ref_dir};
    my $extend_length = $parm{extend_length};

    read_fasta_file(
        FaFile => $V_Href->{virus_orig_genome_fa},
        subrtRef => \&make_singleVirusRef,
        subrtParmAref => \@_,
    );

    # inform
    stout_and_sterr "[INFO]\tmake extended ($extend_length nt once if) virus genome reference.\n";
    my @extend_virus = map {$_->id} grep $_->is_circular, sort {$a->no<=>$b->no} values %{$V_Href->{virusRefSeg}};
    stout_and_sterr "[INFO]\tthese circular virus genomes are extended by ${extend_length} nt.\n"
                         ."\t@extend_virus\n" if scalar @extend_virus;

    # generate virus analysis ref report
    &make_virusRefReport(
        report => File::Spec->catfile($virus_ref_dir, 'virus.ref.info'),
        hasmap => 0
    );
}

#--- make reference for single virus ---
sub make_singleVirusRef{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $vType = $parm{segName};
    my $vSeqSf = $parm{segSeq_Sref};
    my $virus_ref_dir = $parm{virus_ref_dir};
    my $extend_length = $parm{extend_length};

    # virus object
    unless(exists $V_Href->{virusRefSeg}->{$vType}){
        stout_and_sterr "[WARN]\tThe virus_refseg ($vType) is not in the virus_info file. Skip it.\n";
        return;
    }
    my $virus_OB = $V_Href->{virusRefSeg}->{$vType};
    # fasta
    my $name = 'VrefNO'.$virus_OB->no;
    my $singleVirusRefFolder = File::Spec->catfile($virus_ref_dir, $name);
    reset_folder(folder => $singleVirusRefFolder);
    my $singleVirusRefFa = File::Spec->catfile($singleVirusRefFolder, "$name.fa");
    $virus_OB->set_ana_ref(fa => $singleVirusRefFa);
    # extend?
    $virus_OB->set_extLen(extLen => $virus_OB->is_circular ? $extend_length : 0);
    # create the virus fasta reference file
    write_fasta_file(
        SeqSref => $vSeqSf,
        FaFile  => $singleVirusRefFa,
        SegName => $vType,
        LineBase => $V_Href->{line_base},
        CircleExtLen => $virus_OB->extLen
    );
    # create BWA reference for alignment
    BWA_index_fasta(FaFile => $singleVirusRefFa, bwa => $V_Href->{bwa});
    Faidx_Dict_fasta(FaFile => $singleVirusRefFa, samtools => $V_Href->{samtools});
}

#--- write virus ref report ---
sub make_virusRefReport{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $report = $parm{report};
    my $hasmap = $parm{hasmap} || 0;

    open (VRFRP, Try_GZ_Write($report)) || die "fail read $report: $!\n";
    my @theme_item = ('Vref_NO', 'Virus_type', 'Orig_length', 'Ext_length', 'Circle_or_not', 'Ref_fa');
    push @theme_item, 'Init_alignBam' if $hasmap;
    print VRFRP '#' . join("\t", @theme_item ) . "\n"; # theme
    for my $virus_OB (sort {$a->no<=>$b->no} values %{$V_Href->{virusRefSeg}}){
        my @content = ( $virus_OB->no,
                        $virus_OB->id,
                        $virus_OB->length,
                        $virus_OB->extLen,
                        $virus_OB->is_circular ? 'yes' : 'no',
                        File::Spec->abs2rel($virus_OB->ana_ref, dirname($report))
                      );
        if($hasmap){
            if(-e $virus_OB->init_bam->filepath){
                push @content, File::Spec->abs2rel($virus_OB->init_bam->filepath, dirname($report));
            }
            else{
                push @content, 'NA';
            }
        }
        print VRFRP join("\t", @content)."\n";
    }
    close VRFRP;
}

#--- load virus initial alignment report ---
# just make up certain attributes
# MUST run load_virusBasicInfo before to load most attributes
sub load_virusRefReport{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $report = $parm{report};
    my $hasmap = $parm{hasmap} || 0;

    open (VRFRP,Try_GZ_Read($report)) || die "fail read $report: $!\n";
    # theme
    (my $theme_line = lc(<VRFRP>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<VRFRP>){
        next if(/^\#/);
        my @info = split;
        my %vOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        if(!exists $V_Href->{virusRefSeg}->{$vOb{virus_type}}){
            cluck_and_exit "<ERROR>\tcannot find $vOb{virus_type} object.\n";
        }
        my $virus_OB = $V_Href->{virusRefSeg}->{$vOb{virus_type}};
        $virus_OB->set_extLen(extLen => $vOb{ext_length});
        $virus_OB->set_ana_ref(fa => File::Spec->rel2abs($vOb{ref_fa}, dirname($report)));
        unless(file_exist(filePath => $virus_OB->ana_ref)){
            cluck_and_exit "<ERROR>\tcannot find genome reference of virus $vOb{virus_type}.\n"
                                 ."\t".$virus_OB->ana_ref."\n";
        }
        if($hasmap){
            my $init_bamPath = File::Spec->rel2abs($vOb{init_alignbam}, dirname($report));
            if(file_exist(filePath => $init_bamPath)){
                my $init_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $init_bamPath, tag => "initBam");
                $init_bam->addTool(samtools => $V_Href->{samtools});
                $virus_OB->set_init_bam(bam => $init_bam);
            }
            else{
                cluck_and_exit "<ERROR>\tcannot find initial alignment bam of virus $vOb{virus_type}.\n"
                                     ."\t$vOb{init_mapbam}\n";
            }
        }
    }
    close VRFRP;
}

#--- 
1; ## tell the perl script the successful access of this module.
