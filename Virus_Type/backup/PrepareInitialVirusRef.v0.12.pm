package FuseSV::Virus_Type::PrepareInitialVirusRef;

use strict;
use warnings;
use File::Basename qw/ basename dirname /;
use File::Spec qw/ catfile /;
use List::Util qw/ min max sum /;
use Getopt::Long;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success reset_folder /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ write_fasta_file BWA_index_fasta Faidx_Dict_fasta /;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Virus_Type::InitialMapSmartFQonVirus;
use FuseSV::Virus_Type::VirusTypeMain;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              prepare_initial_virus_genome_ref
              load_virus_basic_information
              create_all_virus_IndexRef
              write_virus_analysis_report
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'PrepareInitialVirusRef';
#----- version --------
$VERSION = "0.12";
$DATE = '2019-05-02';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        prepare_initial_virus_genome_ref
                        create_all_virus_IndexRef
                        write_virus_analysis_report
                        load_virus_ref_for_analysis_report
                        load_virus_basic_information
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} cand_virus <[Options]>
     
     Options:

       # Inputs and Outputs #
       # Note: options ('-id' and '-odir') should be as \033[4m \bsame \b\033[0m as previous steps of 'virus_type' #
        -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
                      only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
        -odir   [s]  directory to store outputs, must exist already. <required>
        -vpf    [s]  virus database folder. <required>
                      You can make it via function 'make_virusdb'.
                      OR Download these widely used from FuseSV official website, e.g., HPV.

       # Software Required #
        -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
        -stl    [s]  SamTools, minimum version: 1.3. <required>
        -pic    [s]  Picard jar file, minimum version: 2.1.0. <required>
        -java   [s]  Java compiler, minimum version 1.8. ['java']

        -h|help   Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool;
}

#--- get options from command line
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-id:s"     => \$V_Href->{sample_ID},
        "-odir:s"   => \$V_Href->{virusAnaFolder},
        "-vpf:s"    => \$V_Href->{virus_package_folder},
        # software
        "-bwa:s"    => \$V_Href->{bwa},
        "-stl:s"    => \$V_Href->{samtools},
        "-pic:s"    => \$V_Href->{picard_jar},
        "-java:s"   => \$V_Href->{java},
        # multiple fork operations
        "-fm:i"     => \$V_Href->{fork_number}, # hidden option
        "-jvmem:s"  => \$V_Href->{mem_for_jvm}, # hidden option
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || (!$V_Href->{sample_ID} || $V_Href->{sample_ID} =~ /[^\w\-]/)
             || (!$V_Href->{virusAnaFolder} || !-d $V_Href->{virusAnaFolder})
             || (!$V_Href->{virus_package_folder} || !-d $V_Href->{virus_package_folder})
             || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
             || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
             || (!$V_Href->{picard_jar} || !-e $V_Href->{picard_jar})
            );
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub prepare_initial_virus_genome_ref{
    # file path
    $V_Href->{virus_ref_dir} = GetPath( filekey => 'virusRefFolderForInitAlign' );
    $V_Href->{virus_for_analysis_ref_report} = GetPath( filekey => 'virusRefInitAlign_report' );
    $V_Href->{Init_align_folder} = GetPath( filekey => 'virusRefInitAlign_folder' );

    #---------------------------#
    # mode 2. construct virus   #
    #  initial reference and    #
    #  align SmartFQ to them    #
    #---------------------------#

    # run mode check
    unless(    $V_Href->{run_mode} == 0 # run_whole mode must run this step
            || $V_Href->{run_mode} == 2 # just this step
          ){
        if(!-e $V_Href->{virus_for_analysis_ref_report}){
            cluck_and_exit "<ERROR>\tCannot find virus reference for analysis report:\n".
                          "\t $V_Href->{virus_for_analysis_ref_report}\n".
                          "\tplease run command 'cand_virus'.\n";
        }
        else{
            &load_virus_ref_for_analysis_report(
                    { CheckInitBam => 1 }
                );
        }
        return;
    }

    #------------------#
    #  load virus info #
    #  from database   #
    #------------------#

    &load_virus_basic_information;

    #---------------------#
    #  construct virus    #
    #  initial reference  #
    #---------------------#

    # prepare folder
    reset_folder(folder => $V_Href->{virus_ref_dir});
    # create the Virus specific genome
    &create_all_virus_IndexRef;

    #------------------------#
    #  initial alignment of  #
    #  smartFQ on each virus #
    #------------------------#

    # initial alignment of smartFQ to virus ref
    FuseSV::Virus_Type::InitialMapSmartFQonVirus->initial_align_smfq_on_all_vref;

    #-------------------------#
    # generate report of each #
    # virus and alignment     #
    #------------------------#

    # create virus analysis ref report
    &write_virus_analysis_report;

    # run mode check, if not go on, just stop here.
    if( $V_Href->{run_mode} == 2 ){ # finish at this point in mode 2
        exit(0);
    }
}

#--- create fasta file for all virus(es) ---
sub create_all_virus_IndexRef{
    # options
    shift if ($_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $virus_ref_dir = $parm{virus_ref_dir} || $V_Href->{virus_ref_dir};
    my $extend_length = $parm{extend_length} || $V_Href->{max_rlen};

    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    # create the Virus specific genome
    my $Virus_Ref_NO = 0;
    my @extend_virus;
    open (VOGF,Try_GZ_Read($V_Href->{virus_orig_genome_fa})) || die "fail read $V_Href->{virus_orig_genome_fa}: $!\n";
    $/ = '>'; <VOGF>; $/ = "\n";
    while(<VOGF>){
        chomp(my $refseg = $_);
        $/ = '>';
        chomp(my $segseq = <VOGF>);
        $/ = "\n";
        $segseq =~ s/\s//g;
        # check
        if(!exists($Virus_Info_Href->{$refseg})){
            stout_and_sterr `date`."[WARN]\tThe virus_refseg ($refseg) is not in the virus_info file. Skip it.\n";
            next;
        }

        # create this virus reference for analysis
        $Virus_Ref_NO++;
        my $Vrefname = "VrefNO$Virus_Ref_NO";
        my $single_virus_ref_folder = File::Spec->catfile($virus_ref_dir,"$Vrefname");
        my $single_virus_ref_fa = File::Spec->catfile($single_virus_ref_folder,"$Vrefname.fa");
        # record information
        $Virus_Info_Href->{$refseg}->{vref_no} = $Virus_Ref_NO;
        $Virus_Info_Href->{$refseg}->{vref_fa} = $single_virus_ref_fa;
        # prepare folder
        `rm -rf $single_virus_ref_folder`;

        # circular or not
        if($Virus_Info_Href->{$refseg}->{circle} =~ /yes/i){
            push @extend_virus, $refseg;
            $Virus_Info_Href->{$refseg}->{extend_len} = $extend_length;
        }
        else{
            $Virus_Info_Href->{$refseg}->{extend_len} = 0;  
        }

        # test whether the extended single virus indexed ref in database
        my $single_virus_extend_ref_db_dir = GetPath( filekey => 'single_extend_vref_dir', vref_extLen => $extend_length, virus_subtype => $Vrefname );
        if( -d $single_virus_extend_ref_db_dir ){ # just link from the virus database package
            `ln -sf $single_virus_extend_ref_db_dir $single_virus_ref_folder`;
            stout_and_sterr `date`."[INFO]\tLink prepared Indexed Ref (extend ${extend_length}nt) for $Vrefname.\n";
        }
        else{ # create in local workspace
            `mkdir -p $single_virus_ref_folder`;
            # create the virus fasta reference file
            write_fasta_file(
                SeqSref => \$segseq,
                FaFile  => $single_virus_ref_fa,
                SegName => $refseg,
                LineBase => $V_Href->{line_base},
                CircleExtLen => $Virus_Info_Href->{$refseg}->{extend_len}
            );
            # create BWA reference for alignment
            BWA_index_fasta(FaFile => $single_virus_ref_fa, bwa => $V_Href->{bwa});
            Faidx_Dict_fasta(FaFile => $single_virus_ref_fa, samtools => $V_Href->{samtools});
        }
    }
    close VOGF;

    # if extend
    if(scalar(@extend_virus)){
        stout_and_sterr `date`."[INFO]\tThese virus(es) are circular, and extend their genomes ${extend_length} nt.\n".
                               "\t@extend_virus\n";
    }

    # inform
    stout_and_sterr `date`."[INFO]\tcreate virus genome reference for analysis and BwaIndex-Fai-Dict ok.\n";
}

#--- write virus analysis ref report ---
sub write_virus_analysis_report{
    # options
    shift if ($_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $report = $parm{report} || $V_Href->{virus_for_analysis_ref_report};
    my $nonAlign = $parm{nonAlign} || 0;

    my $virus_ref_dir = dirname($report);
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

    open (VRFRP,Try_GZ_Write($report)) || die "fail $report: $!\n";
    my @theme_item = ('Vref_NO', 'Vref_name', 'Orig_length', 'Extd_length', 'circle_or_not', 'fa_ref', 'init_mapbam');
    pop @theme_item if($nonAlign); # skip the alignment results
    print VRFRP '#' . join("\t", @theme_item ) . "\n"; # theme
    for my $refseg (sort {$Virus_Info_Href->{$a}->{vref_no} <=> $Virus_Info_Href->{$b}->{vref_no}} keys %$Virus_Info_Href){
        (my $vref_fa = $Virus_Info_Href->{$refseg}->{vref_fa}) =~ s/^$virus_ref_dir\///; # this is for linux
        print VRFRP        $Virus_Info_Href->{$refseg}->{vref_no};
        print VRFRP "\t" . $refseg;
        print VRFRP "\t" . $Virus_Info_Href->{$refseg}->{length};
        print VRFRP "\t" . $Virus_Info_Href->{$refseg}->{extend_len};
        print VRFRP "\t" . $Virus_Info_Href->{$refseg}->{circle};
        print VRFRP "\t" . $vref_fa;
        # alignment results
        if( !$nonAlign ){
            if( exists($Virus_Info_Href->{$refseg}->{init_mapbam}) ){
                print VRFRP "\t" . basename($Virus_Info_Href->{$refseg}->{init_mapbam});
            }
            else{
                print VRFRP "\t" . 'NA';
            }
        }
        # new line
        print VRFRP "\n";
    }
    close VRFRP;

    # inform
    stout_and_sterr `date`."[INFO]\twrite virus genome reference report ok.\n"
                                ."\t$report\n";
}

#--- load virus ref for analysis report ---
sub load_virus_ref_for_analysis_report{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    open (VRFRP,Try_GZ_Read($V_Href->{virus_for_analysis_ref_report})) || die "fail read $V_Href->{virus_for_analysis_ref_report}: $!\n";
    # theme
    (my $theme_line = lc(<VRFRP>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<VRFRP>){
        next if(/^\#/);
        my @info = split;
        my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {vref_no} = $rgOb{vref_no};
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {length}  = $rgOb{orig_length};
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {extend_len}  = $rgOb{extd_length};
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {circle}  = $rgOb{circle_or_not};
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {vref_fa} = File::Spec->catfile($V_Href->{virus_ref_dir}, $rgOb{fa_ref});
        $Virus_Info_Href->{ $rgOb{vref_name} } -> {init_mapbam}  = $rgOb{init_mapbam};
        # check the existence of initial alignment bam file
        if($Option_Href->{CheckInitBam} && $rgOb{init_mapbam} !~ /^(EMPTY|SWEEP|NA)$/){
            $rgOb{init_mapbam} = File::Spec->catfile($V_Href->{Init_align_folder}, basename($rgOb{init_mapbam}));
            if(    !-e $rgOb{init_mapbam}
                && $V_Href->{run_mode} == 3 # ReconstrcutVirusVariants
            ){
                cluck_and_exit "<ERROR>:\tCannot find initial alignment bam file at line $. in virus_for_analysis_ref_report.\n".
                              "\t $rgOb{init_mapbam}\n";
            }
            else{
                $Virus_Info_Href->{ $rgOb{vref_name} } -> {init_mapbam}  = $rgOb{init_mapbam};
            }
        }
    }
    close VRFRP;

    # inform
    stout_and_sterr `date`."[INFO]\tload virus ref for analysis report ok.\n";
}

#--- load virus information for further operations ---
sub load_virus_basic_information{
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    open (VIF,Try_GZ_Read($V_Href->{virus_info_doc})) || die "fail read $V_Href->{virus_info_doc}: $!\n";
    # theme
    #subtype        GeneBank_ID     Length  Circular
    (my $theme_line = lc(<VIF>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<VIF>){
        next if(/^\#/);
        my @info = split;
        my %vOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        $Virus_Info_Href->{ $vOb{subtype} } -> {length} = $vOb{length};
        $Virus_Info_Href->{ $vOb{subtype} } -> {circle} = $vOb{circular};
    }
    close VIF;
    stout_and_sterr `date`."[INFO]\tload virus information from database packaege ok.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
