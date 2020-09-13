package FuseSV::Virus_Type::ReconstrcutVirusVariants;

use strict;
use warnings;
use File::Spec qw/ catfile abs2rel /;
use File::Basename qw/ dirname basename /;
use List::Util qw/ min max sum /;
use Getopt::Long;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ BWA_index_fasta Faidx_Dict_fasta /;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Virus_Type::SmartAlignment;
use FuseSV::Virus_Type::VCFwork;
use FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo;
use FuseSV::Virus_Type::BlastM0Mutation;
use FuseSV::Virus_Type::VirusTypeMain;
use Parallel::ForkManager;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              reconstruct_virus_variants
              loadon_virus_variant_report
              merge_bam_and_Nsort
              filter_Nsort_bam_and_to_next_variant_bam
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ReconstrcutVirusVariants';
#----- version --------
$VERSION = "0.48";
$DATE = '2018-01-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        reconstruct_virus_variants
                        rearrange_reads_to_proper_virus_variant
                        assign_pid_to_variant_NO
                        arrange_batch_reads_to_vNO
                        gather_virus_variant_info
                        resort_virus_variant_by_depth
                        compare_this_variant_with_V1
                        write_virus_variant_report
                        loadon_virus_variant_report
                        get_next_candidate_virus
                        get_depth_coverage_given_bam
                        filter_Nsort_bam_and_to_next_variant_bam
                        compare_prev_current_alignment_to_output
                        merge_bam_and_Nsort
                        count_virus_reads_in_bam
                        virus_variant_reconstruct_pipe
                        get_variant_finalBam_from_mkdupBam
                        summarize_virus_variant
                        summarize_all_mutations
                        extended_del_to_followed_reg
                        final_VCF_of_virus_variant
                        Construct_Mut2Orig_PosMapping
                        inform_modifications_in_one_cycle
                        output_modification
                        filter_mutation_from_VCF
                        blast_find_back_mutations
                        extract_reads_from_mut_region
                        extract_InDel_SoftClip_Sites_from_BAM
                        blast_new_to_orig_ref_get_mutations
                        check_whether_mutations_completed
                        merge_continuousSNV_to_MNP
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
     Usage:   perl $V_Href->{MainName} get_variant <[Options]>
     
     Options:

       # Inputs and Outputs #
       # Note: options ('-id' and '-odir') should be as \033[4m \bsame \b\033[0m as previous steps of 'virus_type' #
        -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
                      only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
        -odir   [s]  directory to store outputs, must exist already. <required>

       # Software Required #
        -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
        -stl    [s]  SamTools, minimum version: 1.3. <required>
        -btl    [s]  bcftools, minimum version: 1.3. <required>
        -gatk   [s]  GATK jar file, minimum version: 3.5. <required>
        -pic    [s]  Picard jar file, minimum version: 2.1.0. <required>
        -java   [s]  Java compiler, minimum version 1.8. ['java']
        -blta   [s]  blastall of blast. <required>
        -fordb  [s]  formatdb of blast. <required>

       # Reconstruction Options #
        -maxvn  [i]  the maximum amount of virus variants allowed to reconstruct. [3]
        -mincv  [f]  minimum coverage to consider one virus as a candidate variant. [0.5]
        -unprop      count the un-proper mapped PE reads into reads-counting to find next virus to deal. [Disabled]
        -nonsup      count the supplementary mapped reads into reads-counting. [Disabled]
        -cminmq [i]  count the reads whose mapQ NOT lower than INT into reads-counting. [30]
        -cmaxsp [f]  count the reads whose percent of SoftClip part lower than FLOAT into reads-counting. [0]
        -mindp  [i]  the minimum depth to accept one position on virus genome, else use 'N'. [1]
        -bstrd  [i]  set non-zero to require mutations supported by at least INT reads both strands, repectively. [1]
                      a) set as 0 means supported by reads from only one map-strand will also be good mutation.
                      b) note that, generally, we assign non-zero value; but for virus with linear genome (e.g., HPV59),
                         we'd better set it zero, because specially for the mutations located at the bilateral ends of
                         viral genome, which only provide either of reads-mapping orientations (e.g., 'F' or 'R').
        -minmr  [i]  the minimum amount of reads supporting to accept one mutation, sum of the both strands. [2]
        -use_idv     use 'IDV' value especially for sites have single INDEL mutation. [Disabled]
        -mindrt [f]  if NO.(>=2) virus variant(s) exist, the average_depth must reach FLOAT of that of NO.1. [0, <=1]
                      default value is 0, which means no requirement for the NO.(>=2) virus variant(s).
        -disdup      identical function to this in 'cand_reads' when need to re-generate candidate reads. [Disabled]

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

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # get options
    GetOptions(
        # input/output
        "-id:s"     => \$V_Href->{sample_ID},
        "-odir:s"   => \$V_Href->{virusAnaFolder},
        # software
        "-bwa:s"    => \$V_Href->{bwa},
        "-stl:s"    => \$V_Href->{samtools},
        "-btl:s"    => \$V_Href->{bcftools},
        "-java:s"   => \$V_Href->{java},
        "-gatk:s"   => \$V_Href->{gatk_jar},
        "-pic:s"    => \$V_Href->{picard_jar},
        "-blta:s"   => \$V_Href->{blastall},
        "-fordb:s"  => \$V_Href->{formatdb},
        # options
        ## multiple fork operations
        "-fm:i"     => \$V_Href->{fork_number}, # hidden option
        "-jvmem:s"  => \$V_Href->{mem_for_jvm}, # hidden option
        ## get_variant of virus
        "-maxvn:i"  => \$V_Href->{max_vNO},
        "-mincv:f"  => \$V_Href->{min_coverage},
        "-unprop"   => \$V_Href->{count_unproper_map},
        "-nonsup"   => \$V_Href->{count_non_supple},
        "-cminmq:i" => \$V_Href->{count_min_MAQ},
        "-cmaxsp:f" => \$V_Href->{count_max_Sratio},
        "-mindp:i"  => \$V_Href->{min_depth},
        "-bstrd:i"  => \$V_Href->{both_strd},
        "-minmr:i"  => \$V_Href->{min_mut_reads},
        "-use_idv"  => \$V_Href->{consider_IDV_bool},
        "-maxct:i"  => \$V_Href->{max_cycle_time}, # hidden option
        "-mindrt:f" => \$V_Href->{min_depth_ratio_of_V1},
        ## re-generate candidate reads
        "-disdup"   => \$V_Href->{discard_DupReads},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    return  (   $V_Href->{HELP}
             || (!$V_Href->{sample_ID} || $V_Href->{sample_ID} =~ /[^\w\-]/)
             || (!$V_Href->{virusAnaFolder} || !-d $V_Href->{virusAnaFolder})
             || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
             || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
             || (!$V_Href->{bcftools} || !-e $V_Href->{bcftools})
             || (!$V_Href->{gatk_jar} || !-e $V_Href->{gatk_jar})
             || (!$V_Href->{picard_jar} || !-e $V_Href->{picard_jar})
             || (!$V_Href->{blastall} || !-e $V_Href->{blastall})
             || (!$V_Href->{formatdb} || !-e $V_Href->{formatdb})
             || ($V_Href->{java} ne 'java' && !-e $V_Href->{java})
             || ($V_Href->{min_depth_ratio_of_V1} < 0 || $V_Href->{min_depth_ratio_of_V1} > 1)
             || ($V_Href->{max_vNO} <= 0 || $V_Href->{max_vNO} =~ /\D/)
            );
}

#--- reconstruct virus variants based on BWA+GATK+sample+VCF ---
sub reconstruct_virus_variants{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    # file path
    $V_Href->{Virus_variant_folder} = GetPath( filekey => 'Virus_variant_folder' );
    $V_Href->{Virus_variant_report} = GetPath( filekey => 'Virus_variant_report' );

    #---------------------------#
    # mode 3. reconstruct virus #
    #  variants iteratively     #
    #---------------------------#

    # run mode check
    unless(    $V_Href->{run_mode} == 0 # run_whole mode must run this step
            || $V_Href->{run_mode} == 3 # just this step
          ){
        # just load virus variant report
        &loadon_virus_variant_report;
        return;
    }

    #++++++++++++++++++++++++++++++++++++++++++#
    # stepwised construction of virus variants #
    #++++++++++++++++++++++++++++++++++++++++++#
    my $virus_variant_NO = 0;
    while(1){
        # get virus with maximum reads count
        my  ( $bestVirus_name,
              $bestVirus_coverage,
              $bestVirus_readcount,
              $bestVirus_alignbam,
              $PrepWorkfolder
            )
          = &get_next_candidate_virus(
                {
                    next_variant_NO => $virus_variant_NO + 1
                }
            );
        # if the first virus variant is null, just warn and exit
        if($bestVirus_readcount == 0){
            # sweep
            unless( $V_Href->{in_debug} ){
                `rm -rf $PrepWorkfolder`;
            }
            # jump here
            if($virus_variant_NO == 0){
                stout_and_sterr `date`."[WARN]:\tCannot find any credible virus variant.\n";
                exit(0);
            }
            else{
                stout_and_sterr `date`."[INFO]:\tNo more virus variant can be found.\n";
                last;
            }
        }
        else{
            $virus_variant_NO++;
            # inform
            stout_and_sterr `date`."[INFO]:\tFind $bestVirus_name as NO.$virus_variant_NO virus variant. Reconstruct now.\n"
                                  ."\tCoverage:" . ($bestVirus_coverage * 100) . "%; "
                                  ."AlignReadsCount:$bestVirus_readcount; "
                                  ."SelectionScore:" . ($bestVirus_coverage * $bestVirus_readcount) . ";\n";
            # reconstruct virus variant
            &virus_variant_reconstruct_pipe(
                    {
                        variant_NO  => $virus_variant_NO,
                        var_name    => $bestVirus_name,
                        var_OrigBam => $bestVirus_alignbam
                    }
                );
            # sweep
            unless( $V_Href->{in_debug} ){
                `rm -rf $PrepWorkfolder`;
            }
            # if finish normally
            if( exists($Variant_Info_Href->{$virus_variant_NO}) ){
                # inform
                stout_and_sterr `date`."[INFO]:\tNO.$virus_variant_NO virus variant. Reconstruct done.\n";
                # check the depth comparison of this variant and the first variant
                if( $virus_variant_NO > 1 ){
                    my $discard_this_variant_bool =
                        &compare_this_variant_with_V1(
                                {
                                    variant_NO => $virus_variant_NO,
                                }
                            );
                    # virus variant construction stops here
                    if( $discard_this_variant_bool ){
                        stout_and_sterr "[INFO]:\tSo virus variant construction work stops here.\n";
                        last;
                    }
                }
            }
            else{
                if( $virus_variant_NO == 1 ){
                    stout_and_sterr `date`."[WARN]:\tCannot find any credible virus variant.\n";
                    exit(0);
                }
                else{
                    last;
                }
            }
        }
        # if reach maximum virus variant NO, stop here.
        if( $virus_variant_NO == $V_Href->{max_vNO} ){
            stout_and_sterr `date`."[WARN]:\tReach the maximum amount of virus variant(s) allowed to reconstruct.\n"
                                  ."\tso, virus reconstruction work stops here.\n";
            last;
        }
    }

    #+++++++++++++++++++++++++++++++++++++++++++++++++#
    # try rearrange each read to proper virus variant #
    #+++++++++++++++++++++++++++++++++++++++++++++++++#
    if($virus_variant_NO > 1){
        &rearrange_reads_to_proper_virus_variant;
    }

    #+++++++++++++++++++++++++++++++++++++++++++++++++#
    # try rearrange each read to proper virus variant #
    #+++++++++++++++++++++++++++++++++++++++++++++++++#
    if($virus_variant_NO != 0){
        &gather_virus_variant_info;
    }

    #++++++++++++++++++++++++++++++#
    # resort the variants by depth #
    #++++++++++++++++++++++++++++++#
    # have at least two variants to sort
    if($virus_variant_NO > 1){
        &resort_virus_variant_by_depth;
    }

    #++++++++++++++++++++++++++#
    # output the modifications #
    #++++++++++++++++++++++++++#
    if($virus_variant_NO != 0){
        &write_virus_variant_report;
    }

    # run mode check, if not go on, just stop here.
    if( $V_Href->{run_mode} == 3 ){ # finish at this point in mode 3
        exit(0);
    }
}

#--- check the depth comparison of this variant and the first variant ---
sub compare_this_variant_with_V1{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $virus_variant_NO = $Option_Href->{variant_NO};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
    my $V1_variant_Href = $Variant_Info_Href->{1};

    # check whether to keep this variant
    my $discard_this_variant = 0;
    ## check the average_depth
    if( $this_variant_Href->{AverageDepth} < $V_Href->{min_depth} ){
        $discard_this_variant = 1;
    }
    ## compare the ratio with V1
    if(      $V_Href->{min_depth_ratio_of_V1} > 0
        && ( $this_variant_Href->{AverageDepth} < $V1_variant_Href->{AverageDepth} * $V_Href->{min_depth_ratio_of_V1} )
      ){
        $discard_this_variant = 2;
    }
    ## check the final coverage
    my $coverage = $this_variant_Href->{NonN_length} / $this_variant_Href->{Orig_length};
    if( $coverage < $V_Href->{min_coverage} ){
        $discard_this_variant = 3;
    }

    # operation to discard
    if( $discard_this_variant ){
        ## deletion
        unless( $V_Href->{in_debug} ){
            `rm -rf $this_variant_Href->{workfolder}`;
        }
        delete $Variant_Info_Href->{$virus_variant_NO};
        # inform
        if($discard_this_variant == 1){
            stout_and_sterr `date`."[WARN]:\tDiscard NO.$virus_variant_NO, as its "
                                  ."average_depth ($this_variant_Href->{AverageDepth}) is less than "
                                  ."min_depth ($V_Href->{min_depth}).\n";
        }
        elsif($discard_this_variant == 2){
            stout_and_sterr `date`."[WARN]:\tDiscard NO.$virus_variant_NO, as its "
                                  ."average_depth ($this_variant_Href->{AverageDepth}) is less than "
                                  ."min_depth_ratio ($V_Href->{min_depth_ratio_of_V1}) of "
                                  ."average_depth ($V1_variant_Href->{AverageDepth}) of NO.1 virus variant.\n";
        }
        elsif($discard_this_variant == 3){
            stout_and_sterr `date`."[WARN]:\tDiscard NO.$virus_variant_NO, as its "
                                  ."coverage ($coverage) is less than "
                                  ."min_coverage ($V_Href->{min_coverage}).\n";
        }
    }

    return $discard_this_variant;
}

#--- read virus variant report ---
sub loadon_virus_variant_report{
    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    if(!-e $V_Href->{Virus_variant_report}){
        cluck_and_exit "<ERROR>\tCannot find virus variants report:\n".
                      "\t $V_Href->{Virus_variant_report}\n".
                      "\tplease run command 'get_variant' or 'virus_type'.\n";
    }

    open (VVRP,Try_GZ_Read($V_Href->{Virus_variant_report})) || die "fail read $V_Href->{Virus_variant_report}: $!\n";
    # theme
    (my $theme_line = lc(<VVRP>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<VVRP>){
        next if(/^\#/);
        my @info = split;
        my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        $Variant_Info_Href->{ $rgOb{no} } -> {CycleTime}   = $rgOb{dealcycletime};
        $Variant_Info_Href->{ $rgOb{no} } -> {OrigRefSeg}  = $rgOb{originalvirus};
        $Variant_Info_Href->{ $rgOb{no} } -> {RefSeg}      = $rgOb{variantvirus};
        $Variant_Info_Href->{ $rgOb{no} } -> {Orig_length} = $rgOb{origlength};
        $Variant_Info_Href->{ $rgOb{no} } -> {NonN_length} = $rgOb{nonn_length};
        $Variant_Info_Href->{ $rgOb{no} } -> {Cricle}      = $rgOb{cricleornot};
        $Variant_Info_Href->{ $rgOb{no} } -> {Extend_length} = $rgOb{toextlength};
        #!!! this in-turn column is Coverage, skipped
        $Variant_Info_Href->{ $rgOb{no} } -> {AverageDepth} = $rgOb{averagedepth};
        #!!! this in-turn column is Modifications, skipped
        $Variant_Info_Href->{ $rgOb{no} } -> {ModCompBool} = ($rgOb{modcomp} eq 'Y');
        $Variant_Info_Href->{ $rgOb{no} } -> {Fa}     = File::Spec->catfile( $V_Href->{Virus_variant_folder}, $rgOb{reffasta} );
        $Variant_Info_Href->{ $rgOb{no} } -> {mkupFa} = File::Spec->catfile( $V_Href->{Virus_variant_folder}, $rgOb{makeupfasta} );
        $Variant_Info_Href->{ $rgOb{no} } -> {finalBam} = File::Spec->catfile( $V_Href->{Virus_variant_folder}, $rgOb{alignbam} );
        $Variant_Info_Href->{ $rgOb{no} } -> {Mods}   = File::Spec->catfile( $V_Href->{Virus_variant_folder}, $rgOb{moddetails} );

        # cover ratio
        $Variant_Info_Href->{ $rgOb{no} } -> {CoverRatio} = $rgOb{nonn_length} / $rgOb{origlength};

        # deal 'Coverage'
        for my $coverage_term (split /;/, $rgOb{coverage}){
            my ($depth, $coverage) = ($coverage_term =~ /(\d+)X:([\d\.]+)\%/i);
            $Variant_Info_Href->{ $rgOb{no} } -> {Coverage} -> {$depth} = $coverage;
        }

        # deal 'Modifications'
        ## however, we cannot load back the position-wise mutation;
        ## here, we only just load back the count of each mutation type
        for my $mods_term (split /;/, $rgOb{modifications}){
            my ($type, $count) = ($mods_term =~ /([a-z]+):(\d+)/i);
            $Variant_Info_Href->{ $rgOb{no} } -> {Mutation} -> {$type} = $count;
        }

        # check the existence of initial alignment bam file
        if( !-e $Variant_Info_Href->{ $rgOb{no} }->{mkupFa} ){
            cluck_and_exit "<ERROR>:\tCannot find virus variant makeup fasta file at line $. in virus_variant_report.\n".
                          "\t $Variant_Info_Href->{ $rgOb{no} }->{mkupFa}\n";
        }
        elsif( !-e $Variant_Info_Href->{ $rgOb{no} }->{finalBam} ){
            cluck_and_exit "<ERROR>:\tCannot find virus variant alignment bam file at line $. in virus_variant_report.\n".
                          "\t $Variant_Info_Href->{ $rgOb{no} }->{finalBam}\n";
        }
    }
    close VVRP;

    # inform
    stout_and_sterr `date`."[INFO]\tload virus variant report ok.\n";
}

#--- get the virus with most aligned segments as the next virus variant ---
sub get_next_candidate_virus{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $next_variant_NO = $Option_Href->{next_variant_NO};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Fork_Number = $V_Href->{fork_number};

    # folder
    my $PrepWorkfolder = File::Spec->catfile( $V_Href->{Virus_variant_folder}, "Prepare_for_V$next_variant_NO" );
    `rm -rf $PrepWorkfolder`;
    `mkdir -p $PrepWorkfolder`;

    # prepare information of the former variant(s)
    my @prev_vbam;
    my %prev_v_refseg;
    my %prev_v_origrefseg;
    for my $prev_vNO (sort {$a<=>$b} keys %$Variant_Info_Href){
        my $bam = $Variant_Info_Href->{$prev_vNO}->{finalBam};
        my $refseg = $Variant_Info_Href->{$prev_vNO}->{RefSeg};
        my $origrefseg = $Variant_Info_Href->{$prev_vNO}->{OrigRefSeg};
        push @prev_vbam, $bam;
        $prev_v_refseg{$refseg} = $prev_vNO;
        $prev_v_origrefseg{$origrefseg} = 1;
    }

    # count qulitified reads amount
    my $bestVirus_refseg = '';
    my $bestVirus_coverage = 0;
    my $bestVirus_readcount = 0;
    my $bestVirus_alignbam = '';

    # pick out the available virus ref segments
    my @available_v_refseg = grep    !exists($prev_v_origrefseg{$_})
                                  && $Virus_Info_Href->{$_}->{init_mapbam} !~ /^(EMPTY|SWEEP|NA)$/,
                             sort {$Virus_Info_Href->{$a}->{vref_no} <=> $Virus_Info_Href->{$b}->{vref_no}} keys %$Virus_Info_Href;

    # potential multiple fork operations, initialization
    my $multiple_fork;
    if($Fork_Number > 1){
        $multiple_fork = new Parallel::ForkManager($Fork_Number);
        # inform
        stout_and_sterr `date`."[INFO]:\tOperate in multiple forks mode. Fork Number: $Fork_Number.\n";
    }

    # check virus one by one
    for my $v_refseg ( @available_v_refseg ){

        # potential multiple fork operations, child fork starts
        if($Fork_Number > 1){
            $multiple_fork->start and next;
        }

        # virus info
        my $vref_no = $Virus_Info_Href->{$v_refseg}->{vref_no};
        my $vref_BwaRef = $Virus_Info_Href->{$v_refseg}->{vref_fa};
        my $source_bam = $Virus_Info_Href->{$v_refseg}->{init_mapbam};
        my $out_pref = File::Spec->catfile( $PrepWorkfolder, "NO.$vref_no" );
        my $bam_for_reads_count = $out_pref.".for_reads_count.sort.bam";

        #++++++++++++++++++++++++++++++#
        # merge with previous virus    #
        # variant before, select reads #
        #++++++++++++++++++++++++++++++#
        if($next_variant_NO == 1){
            `ln -sf $Virus_Info_Href->{$v_refseg}->{init_mapbam} $bam_for_reads_count`;
            `ln -sf $Virus_Info_Href->{$v_refseg}->{init_mapbam}.bai $bam_for_reads_count.bai`;
        }
        else{
            # merge bam and sort by read name
            my $merge_nSortBam =
            &merge_bam_and_Nsort({
                bam_pref => $out_pref,
                source_bam_Aref => [ $source_bam ],
                prev_bam_Aref => \@prev_vbam
            });
            # filter the Nsort bam and output the remained bam
            &filter_Nsort_bam_and_to_next_variant_bam({
                nSortBam => $merge_nSortBam,
                OutBam => $bam_for_reads_count,
                BwaRef => $vref_BwaRef,
                Prev_vref_Href => \%prev_v_refseg,
                keep_SQ => $v_refseg
            });
            # sweep
            `rm -rf $merge_nSortBam`;
        }

        #++++++++++++++++++++++++++#
        # count depth and coverage #
        #++++++++++++++++++++++++++#
        my $coverage =
        &get_depth_coverage_given_bam({
            Bam => $bam_for_reads_count,
            check_seg => $v_refseg,
            seg_orig_length => $Virus_Info_Href->{$v_refseg}->{length},
            seg_ext_length => $Virus_Info_Href->{$v_refseg}->{extend_len}
        });
        if( $coverage < $V_Href->{min_coverage}){
            # inform
            stout_and_sterr `date`."[INFO]:\tcoverage of $v_refseg is " . ($coverage * 100) . "%, bad for candidate virus variant NO.$next_variant_NO selection.\n";
        }
        else{
            #+++++++++++++#
            # reads count #
            #+++++++++++++#
            # count!
            my $reads_count =
            &count_virus_reads_in_bam({
                bam_for_reads_count => $bam_for_reads_count,
                proper_map => (! $V_Href->{count_unproper_map}),
                non_supple => (! $V_Href->{count_non_supple}),
                min_MAQ => $V_Href->{count_min_MAQ},
                max_Sratio => $V_Href->{count_max_Sratio}
            });
            # inform
            stout_and_sterr `date`."[INFO]:\tcoverage of $v_refseg is " . ($coverage * 100) . "%, ok for candidate virus variant NO.$next_variant_NO selection,\n"
                                  ."\tand there are $reads_count available reads aligned on $v_refseg.\n"
                                  ."\t[only_proper_map:" . ($V_Href->{count_unproper_map} ^ 1) . ';'
                                  .  "no_supple_map:"    . ($V_Href->{count_non_supple} ^ 1)   . ';'
                                  .  "minimum_mapQ:"     . $V_Href->{count_min_MAQ}            . ';'
                                  .  "maximum_Sratio:"   . $V_Href->{count_max_Sratio}         . "]\n"
                                  ."\tSelection_Score is "  . ($coverage * $reads_count) . "\n";
            # record the counts information
            `echo $coverage > $out_pref.Coverage.txt`;
            `echo $reads_count > $out_pref.ReadsCount.txt`;
            `echo $bam_for_reads_count > $out_pref.BamForCount.list`;
        }

        # potential multiple fork operations, child fork finishs
        if($Fork_Number > 1){
            $multiple_fork->finish;
        }
    }

    # potential multiple fork operations, wait all childern forks
    if($Fork_Number > 1){
        $multiple_fork->wait_all_children;
    }

    # find the maximum count
    for my $v_refseg ( @available_v_refseg ){
        my $vref_no = $Virus_Info_Href->{$v_refseg}->{vref_no};
        # read from disk
        my $out_pref = File::Spec->catfile( $PrepWorkfolder, "NO.$vref_no" );
        # get coverage
        my $Coverage_txt = "$out_pref.Coverage.txt";
        next if(!-e $Coverage_txt);
        chomp(my $coverage = `cat $Coverage_txt`);
        # get aligned-reads count
        my $ReadsCount_txt = "$out_pref.ReadsCount.txt";
        if(!-e $ReadsCount_txt){
            cluck_and_exit `data`."<ERROR>:\tThe aligned-reads count txt of virus $v_refseg does not exists:\n"
                                ."\t$ReadsCount_txt\n";
        }
        chomp(my $reads_count = `cat $ReadsCount_txt`);
        # Selection Score: coverage * reads_count
        if( $reads_count * $coverage > $bestVirus_coverage * $bestVirus_readcount ){
            $bestVirus_coverage = $coverage;
            $bestVirus_readcount = $reads_count;
            $bestVirus_refseg = $v_refseg;
            chomp($bestVirus_alignbam = `cat $out_pref.BamForCount.list`);
        }
        # sweep
        `rm $Coverage_txt $ReadsCount_txt $out_pref.BamForCount.list`;
    }

    return ( $bestVirus_refseg,
             $bestVirus_coverage,
             $bestVirus_readcount,
             $bestVirus_alignbam,
             $PrepWorkfolder
            );
}

#--- samtools depth to get the coverage ---
sub get_depth_coverage_given_bam{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $bam = $Option_Href->{Bam};
    my $check_seg = $Option_Href->{check_seg};
    my $seg_orig_length = $Option_Href->{seg_orig_length};
    my $seg_nonN_length = $Option_Href->{seg_nonN_length};
    my $seg_ext_length = $Option_Href->{seg_ext_length} || 0;
    my $need_details_bool = $Option_Href->{need_details} || 0;
    # check
    if($need_details_bool && !exists($Option_Href->{coverage_Href})){
        cluck_and_exit "<ERROR>:\tRequires 'coverage_Href' with 'need_details'.\n";
    }

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # read depth
    my %pos2depth;
    open (DEPTH, "$V_Href->{samtools} depth -a -r $check_seg $bam |") || die "fail samtools depth $bam: $!\n";
    while (<DEPTH>){
        my ($pos, $depth) = (split)[1,2];
        $pos2depth{$pos} = $depth;
    }
    close DEPTH;

    # sometimes, all alignments in the bam are duplicated; so no depth info obtained
    if(scalar(keys %pos2depth) == 0){
        return 0;
    }

    # deal circular virus
    if($seg_ext_length != 0){
        for my $i (1 .. $seg_ext_length){
            $pos2depth{$i} += $pos2depth{ $seg_orig_length + $i };
            delete $pos2depth{ $seg_orig_length + $i };
        }
    }

    # return what?
    if($need_details_bool){
        # empty coverage count
        $Option_Href->{coverage_Href}->{$_} = 0 for keys %{$Option_Href->{coverage_Href}};
        # stat on each pos
        my $okPos_sum_depth = 0;
        my $okPos_amount = 0;
        for my $pos (sort {$a<=>$b} keys %pos2depth){
            my $depth = $pos2depth{$pos};
            # coverage 
            for my $aim_depth (sort {$a<=>$b} keys %{$Option_Href->{coverage_Href}}){
                if($depth >= $aim_depth){
                    $Option_Href->{coverage_Href}->{$aim_depth} ++;
                }
                else{
                    last;
                }
            }
            # count available depth
            next if($depth < $V_Href->{min_depth});
            $okPos_amount ++;
            $okPos_sum_depth += $depth;
        }
        # gradient coverage
        for my $aim_depth (sort {$a<=>$b} keys %{$Option_Href->{coverage_Href}}){
            # standard transform, sometimes, 'N' position may have some depth
            $Option_Href->{coverage_Href}->{$aim_depth} = 
                sprintf "%0.4f", ( min($Option_Href->{coverage_Href}->{$aim_depth}, $seg_nonN_length) / $seg_nonN_length);
        }
        # average depth
        if(exists($Option_Href->{averagedepth_Sref})){
            ${$Option_Href->{averagedepth_Sref}} = sprintf "%0.2f", ($okPos_sum_depth / ($okPos_amount || 1) );
        }
    }
    else{
        # coverage judgement
        my $okPos_amount = grep $pos2depth{$_} >= $V_Href->{min_depth}, sort {$a<=>$b} keys %pos2depth;
        my $coverage = sprintf "%0.4f", ($okPos_amount / $seg_orig_length);
        return $coverage;
    }
}

#--- filter Nsort bam and extract the part for next candidate virus variant ---
sub filter_Nsort_bam_and_to_next_variant_bam{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $Nsort_bam = $Option_Href->{nSortBam};
    my $OutBam = $Option_Href->{OutBam};
    my $BwaRef = $Option_Href->{BwaRef} || '__NULL__'; # discard on 2016-05-04
    my $prev_v_refseg_Href = $Option_Href->{Prev_vref_Href}; # it is from $Variant_Info_Href
    my $keep_SQ = $Option_Href->{keep_SQ} || '__NA__';
    my $keep_prev_AlignInfo_bool = $Option_Href->{keep_prev_AlignInfo} || 0;
    my $discard_supple_align_bool = $Option_Href->{dis_SuppleAlign} || 0;
    my $exchange_AlignEnds_bool = $Option_Href->{exchange_AlignEnds} || 0;
    my $mark_nonPrev_SamLine_bool = $Option_Href->{mark_nonPrev_SamLine} || 0;
    my $mark_Prev_SamLine_bool = $Option_Href->{mark_Prev_SamLine} || 0;
    # my $RefSeg_Extend_Href = $Option_Href->{RefSeg_Extend_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $max_rlen = $V_Href->{max_rlen} || 0;

    # -T $BwaRef, discard on 2016-05-04
    open (my $OBAM_FH, "| $V_Href->{samtools} view -b -S -o $OutBam.unsort.bam -") || die "fail write $OutBam.unsort.bam: $!\n";

    # filter bam
    #   0x800, to avoid supplementary alignment
    # + 0x100: to avoid secondary alignment
    # = 0x900
    my $fileter_SuppleAlign = $discard_supple_align_bool ? '-F 0x900':'';
    my $last_pid = 'initial-value'; # this is initial value, will be deleted soon
    my $last_pid_new_vref_bool = 0;
    my $last_pid_prev_vref_bool = 0;
    my %pid_InfoPool = ( $last_pid => {} ); # use the initial value to load the initial hash
    my %only_has_new_vref_pid;
    my $amount_NeedToDeal = 0;
    open (NSBAM, "$V_Href->{samtools} view $fileter_SuppleAlign -h $Nsort_bam |") || die "fail read $Nsort_bam: $!\n";
    while (<NSBAM>){
        # header
        if(/^\@/){
            if( !/^\@SQ/ ){
                print $OBAM_FH $_;
            }
            elsif( $keep_SQ eq '__NA__' || /\sSN:$keep_SQ\s/ ){
                my ($Refseg, $Reflength) = (/\sSN:(\S+)\sLN:(\d+)/);
                # to extend the length?
                if(exists( $Option_Href->{RefSeg_Extend_Href}->{$Refseg} )){
                    $Reflength += $Option_Href->{RefSeg_Extend_Href}->{$Refseg}->{extend_len} - $max_rlen;
                    print $OBAM_FH "\@SQ\tSN:$Refseg\tLN:$Reflength\n";
                }
                else{
                    print $OBAM_FH $_;
                }
            }
            next;
        }
        # line info
        my @lineInfo = split;
        # correct the mapQ, especially for Q255 induced by alignment tool STAR/HISAT2.
        $lineInfo[4] = min( $lineInfo[4], 60 );
        # each attribute
        my ($pid, $flag, $refseg, $map_pos, $mapQ, $Cigar, $Tlen) = @lineInfo[0,1,2,3,4,5,8];
        # new pid to deal, count it
        if( $pid ne $last_pid ){
            # judge the last pid
            if($last_pid_new_vref_bool){ # have the new virus reference, keep to next compare.
                if($last_pid_prev_vref_bool){ # need to compare the prev and new in the following operation
                    # count this pid in, it has both new ref and prev ref
                    $amount_NeedToDeal ++;
                }
                else{ # only has new virus reference, record it for directly output in the following operation
                    $only_has_new_vref_pid{$last_pid} = 1;
                    # count this pid in, it only has new ref
                    $amount_NeedToDeal ++;
                }
            }
            else{ # only have the previous virus reference, OR all are unmaps
                # keep the non-new ref ?
                if( $keep_prev_AlignInfo_bool ){
                    for my $end_NO (sort {$a<=>$b} keys %{$pid_InfoPool{$last_pid}}){
                        for my $refseg (sort keys %{$pid_InfoPool{$last_pid}->{$end_NO}}){
                            for my $alignInfo_Href (@{$pid_InfoPool{$last_pid}->{$end_NO}->{$refseg}}){
                                print $OBAM_FH $alignInfo_Href->{LineInfo};
                                last if($refseg eq '*'); # this must come from both-unmapped paired-end, so just output one time
                            }
                        }
                    }
                }
                # sweep, discard the last pid
                delete $pid_InfoPool{$last_pid};
            }
            # reset
            $last_pid_new_vref_bool = 0;
            $last_pid_prev_vref_bool = 0;
            # deal groups of new-pids
            if( $amount_NeedToDeal >= $V_Href->{Nsort_MaxStore}){ # meet the store threshold
                &compare_prev_current_alignment_to_output(
                        {
                            file_handle => $OBAM_FH,
                            pid_InfoPool_Href => \%pid_InfoPool,
                            only_has_new_pid_Href => \%only_has_new_vref_pid,
                            Prev_vref_Href => $prev_v_refseg_Href,
                            keep_prev_AlignInfo => $keep_prev_AlignInfo_bool,
                            exchange_AlignEnds => $exchange_AlignEnds_bool,
                            mark_nonPrev_SamLine => $mark_nonPrev_SamLine_bool,
                            mark_Prev_SamLine => $mark_Prev_SamLine_bool
                        }
                    );
                # reset
                $amount_NeedToDeal = 0;
            }
        }
        # record
        # my ($edit_dist) = (/NM:i:(\d+)/);
        # my ($pairMapQ) = (/MQ:i:(\d+)/);
        my $end_NO = ($flag & 0x40) ? 1 : 2; # $flag & 0x80 + 1
        my ($align_score) = (/\bAS:i:(\d+)\b/); # even the read is unmap(0x4), this line still contain AS:i:0
        my $alignInfo_Href = {
                LineInfo => join("\t",@lineInfo)."\n",
                Refseg => $refseg,
                MapPos => $map_pos,
                Proper => $flag & 0x2,
                Unmapped => $flag & 0x4,
                RevcmlpMap => $flag & 0x10,
                MapQ => $mapQ, # MapQ of unmapped read is 0
                Cigar => $Cigar, # cigar of unmapped read is '*'
                Tlen => $Tlen, # Tlen of unmapped read is 0
                AlignScore => ($align_score || 0) # AlignScore of unmapped read is 0, unmap may lack AS:i:
                # PairMapQ => $pairMapQ,
                # EditDist => $edit_dist
            };
        # store in the InfoPool
        # note that, if one is mapped, another is unmapped, the unmapped end's refseg must be same as that of the mapped one.
        push @{ $pid_InfoPool{$pid}->{$end_NO}->{$refseg} }, $alignInfo_Href;
        # become the last pid
        $last_pid = $pid;
        # set the bool
        if( exists( $prev_v_refseg_Href->{$refseg} ) ){ # has the previous ref
            $last_pid_prev_vref_bool = 1;
        }
        elsif( !($flag & 0x4) || $refseg ne '*' ){ # [Double check], it is not unmap, So must be the new ref
            $last_pid_new_vref_bool = 1;
        }
    }
    close NSBAM;
    #--------------------------------------#
    # deal the last pid and the last group #
    #--------------------------------------#
    ## last pid
    if( $last_pid_new_vref_bool ){ # only has new virus reference
        if( !$last_pid_prev_vref_bool ){
            $only_has_new_vref_pid{$last_pid} = 1;
        }
    }
    else{ # only have the previous virus reference, OR all are unmaps
        # keep the non-new ref ?
        if( $keep_prev_AlignInfo_bool ){
            for my $end_NO (sort {$a<=>$b} keys %{$pid_InfoPool{$last_pid}}){
                for my $refseg (sort keys %{$pid_InfoPool{$last_pid}->{$end_NO}}){
                    for my $alignInfo_Href (@{$pid_InfoPool{$last_pid}->{$end_NO}->{$refseg}}){
                        print $OBAM_FH $alignInfo_Href->{LineInfo};
                        last if($refseg eq '*'); # this must come from both-unmapped paired-end, so just output one time
                    }
                }
            }
        }
        # sweep, discard the last pid
        delete $pid_InfoPool{$last_pid};
    }
    ## last group
    &compare_prev_current_alignment_to_output(
                    {
                        file_handle => $OBAM_FH,
                        pid_InfoPool_Href => \%pid_InfoPool,
                        only_has_new_pid_Href => \%only_has_new_vref_pid,
                        Prev_vref_Href => $prev_v_refseg_Href,
                        keep_prev_AlignInfo => $keep_prev_AlignInfo_bool,
                        exchange_AlignEnds => $exchange_AlignEnds_bool,
                        mark_nonPrev_SamLine => $mark_nonPrev_SamLine_bool,
                        mark_Prev_SamLine => $mark_Prev_SamLine_bool
                    }
                );

    close $OBAM_FH;

    # sort bam
    my $SortBam_command = "$V_Href->{samtools} sort --threads 3 -m 1G -T $OutBam.sort_tmp -o $OutBam $OutBam.unsort.bam && $V_Href->{samtools} index $OutBam";
    trible_run_for_success($SortBam_command, 'SortBam', {cmd_Nvb=>1, esdo_Nvb=>1});
    `rm -rf $OutBam.unsort.bam`;
}

#--- compare the previous variant and candidate virus variant ---
sub compare_prev_current_alignment_to_output{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $OBAM_FH = $Option_Href->{file_handle};
    my $pid_InfoPool_Href = $Option_Href->{pid_InfoPool_Href};
    my $prev_v_refseg_Href = $Option_Href->{Prev_vref_Href};
    my $only_has_new_pid_Href = $Option_Href->{only_has_new_pid_Href};
    my $keep_prev_AlignInfo_bool = $Option_Href->{keep_prev_AlignInfo} || 0;
    my $exchange_AlignEnds_bool = $Option_Href->{exchange_AlignEnds} || 0;
    my $mark_nonPrev_SamLine_bool = $Option_Href->{mark_nonPrev_SamLine} || 0;
    my $mark_Prev_SamLine_bool = $Option_Href->{mark_Prev_SamLine} || 0;

    for my $pid (keys %$pid_InfoPool_Href){
        # store the output SAM line info
        my $output_Href = {};
        # criterion to compare
        my %unmap_bool;
        my %align_score;
        my %mapQ;
        my %Cigar;
        # next will output the alignment info via this tag
        my $tag_to_output = '__NULL__';
        # load the criterion and prepare outputs
        for my $end_NO (sort {$a<=>$b} keys %{$pid_InfoPool_Href->{$pid}}){
            for my $refseg (keys %{$pid_InfoPool_Href->{$pid}->{$end_NO}}){
                # '*' means both-unmapped paired-end, useless here
                # because we have at least one aligned info on non-previous refseg
                next if($refseg eq '*');
                # which source? tag it.
                my $tag = exists($prev_v_refseg_Href->{$refseg}) ? 'Prev' : 'nonPrev';
                # infor comes
                for my $alignInfo_Href (@{$pid_InfoPool_Href->{$pid}->{$end_NO}->{$refseg}}){
                    # store for output
                    if(    $tag eq 'nonPrev' # prepare for output of non-previous
                        || $keep_prev_AlignInfo_bool # prepare for output of previous
                       ){
                        push @{$output_Href->{$tag}->{$end_NO}}, $alignInfo_Href;
                    }
                    # if only has the new virus reference, no need to store the following data
                    if( exists( $only_has_new_pid_Href->{$pid} ) ){
                        # mark the tag
                        $tag_to_output = 'nonPrev';
                        next;
                    }
                    # unmapped sign
                    $unmap_bool{$tag}->[$end_NO-1] = ($unmap_bool{$tag}->[$end_NO-1] || 0) || $alignInfo_Href->{Unmapped};
                    # Mapping Quality
                    $mapQ{$tag}->[$end_NO-1] = max( ($mapQ{$tag}->[$end_NO-1] || 0), $alignInfo_Href->{MapQ} );
                    # Cigar
                    push @{$Cigar{$tag}->[$end_NO-1]} , $alignInfo_Href->{Cigar};
                    # Align Score
                    $align_score{$tag}->[$end_NO-1] = max( ($align_score{$tag}->[$end_NO-1] || 0), $alignInfo_Href->{AlignScore} );
                }
            }
        }

        # exchange ends?
        ## one end is unmap in Prev, but map in nonPrev; while the other end is in same situation WITH v-refs exchanged.
        if(     $tag_to_output eq '__NULL__'  # not been assigned yet
             && $exchange_AlignEnds_bool
            ){
            if(    $unmap_bool{Prev}->[0]    xor $unmap_bool{Prev}->[1] # Previous align has one unmap
                && $unmap_bool{nonPrev}->[0] xor $unmap_bool{nonPrev}->[1] # non-Previous has one unmap
                && $unmap_bool{nonPrev}->[0] !=  $unmap_bool{Prev}->[0] # the unmapped end is not the same endNO
               ){
                # mark the tag
                $tag_to_output = 'Exchange';
                # prepare the sources
                my $first_end_source  = ($unmap_bool{Prev}->[0]) ? 'nonPrev' : 'Prev';
                my $second_end_source = ($unmap_bool{Prev}->[1]) ? 'nonPrev' : 'Prev';
                # assign, in my opinion, the array must contain only one element.
                $output_Href->{Exchange}->{1} = [ $output_Href->{ $first_end_source  }->{1}->[0] ];
                $output_Href->{Exchange}->{2} = [ $output_Href->{ $second_end_source }->{2}->[0] ];
                # update the information
                for my $endNO (1,2){
                    my $p_endNO = ($endNO==1) ? 2 : 1;
                    my $end_info_Href = $output_Href->{Exchange}->{$endNO}->[0];
                    my $p_end_info_Href = $output_Href->{Exchange}->{$p_endNO}->[0];
                    # info comes
                    my @end_info = split /\t+/, $end_info_Href->{LineInfo};
                    # update
                    ## mated end mapped refseg and position
                    $end_info[6] = $p_end_info_Href->{Refseg};
                    $end_info[7] = $p_end_info_Href->{MapPos};
                    ## set Tlen as zero
                    $end_info[8] = 0;
                    ## flags
                    $end_info[1] += 0x1 unless( $end_info[1] & 0x1 );
                    $end_info[1] -= 0x2     if( $end_info[1] & 0x2 );
                    $end_info[1] -= 0x4     if( $end_info[1] & 0x4 );
                    $end_info[1] -= 0x8     if( $end_info[1] & 0x8 );
                    $end_info[1] += 0x20    if( $p_end_info_Href->{RevcmlpMap}  && !($end_info[1] & 0x20) );
                    $end_info[1] -= 0x20    if( !$p_end_info_Href->{RevcmlpMap} && $end_info[1] & 0x20 );
                    $end_info[1] -= 0x100   if( $end_info[1] & 0x100 );
                    $end_info[1] -= 0x200   if( $end_info[1] & 0x200 );
                    $end_info[1] -= 0x400   if( $end_info[1] & 0x400 );
                    $end_info[1] -= 0x800   if( $end_info[1] & 0x800 );
                    ## MQ
                    my $new_LineInfo = join("\t", @end_info); # last element has "\n"
                    my $p_MQ = $p_end_info_Href->{MapQ};
                    $new_LineInfo =~ s/\bMQ:i:\d+\b/MQ:i:$p_MQ/;
                    # reassign LineInfo
                    $end_info_Href->{LineInfo} = $new_LineInfo;
                }
            }
        }

        # then judge to output whole non-previous or whole previous
        if(       $tag_to_output eq '__NULL__'  # not been assigned yet
               && sum(@{$mapQ{nonPrev}})        >= sum(@{$mapQ{Prev}})
               && sum(@{$align_score{nonPrev}}) >= sum(@{$align_score{Prev}})
           ){
            # calculate the map-length
            my %map_length;
            for my $tag (sort keys %Cigar){
                my $Cigar_tag_Aref = $Cigar{$tag};
                for my $i ( 0 .. scalar(@$Cigar_tag_Aref)-1 ){
                    my $max_mlen = 0;
                    my $cigar_Aref = $Cigar_tag_Aref->[$i];
                    for my $cigar (@$cigar_Aref){
                        next if($cigar eq '*'); # skip unmapped reads
                        (my $cigar_tmp = $cigar) =~ s/\d+[INSH]//g; # replace all non-placeholder
                        my $mlen = sum( split /\D+/, $cigar_tmp );
                        $max_mlen = max($max_mlen, $mlen);
                    }
                    # record
                    $map_length{$tag}->[$i] = $max_mlen;
                }
            }
            # mlen judgement, and at least one criterion is better on new virus reference
            if(     sum(@{$map_length{nonPrev}}) >= sum(@{$map_length{Prev}}) * 1.03
                 && ! (    sum(@{$mapQ{nonPrev}})        == sum(@{$mapQ{Prev}}) 
                        && sum(@{$align_score{nonPrev}}) == sum(@{$align_score{Prev}})
                        && sum(@{$map_length{nonPrev}})  == sum(@{$map_length{Prev}}) * 1.03
                      )
               ){
                # mark the tag
                $tag_to_output = 'nonPrev';
            }
        }

        # not Exchange and not non-Previous, if keep Previous, run this part
        if(     $tag_to_output eq '__NULL__'  # not been assigned yet
             && $keep_prev_AlignInfo_bool
           ){
            # mark the tag
            $tag_to_output = 'Prev';
        }

        # output related alignment data
        if( $tag_to_output ne '__NULL__' ){
            for my $end_NO (sort {$a<=>$b} %{$output_Href->{$tag_to_output}}){
                # print $OBAM_FH $_->{LineInfo} for @{$output_Href->{$tag_to_output}->{$end_NO}}; # old, commented on 20170719
                for my $end_info_Href (@{$output_Href->{$tag_to_output}->{$end_NO}}) {
                    if(    $mark_nonPrev_SamLine_bool
                        && $tag_to_output ne 'Prev' # nonPrev or exchanged
                    ){
                        chomp($end_info_Href->{LineInfo});
                        print $OBAM_FH $end_info_Href->{LineInfo} . "\tZS:i:1\n";
                    }
                    elsif(    $mark_Prev_SamLine_bool
                           && $tag_to_output ne 'nonPrev' # Prev or exchanged
                    ){
                        chomp($end_info_Href->{LineInfo});
                        print $OBAM_FH $end_info_Href->{LineInfo} . "\tZS:i:0\n";
                    }
                    else{
                        print $OBAM_FH $end_info_Href->{LineInfo};
                    }
                }
            }
        }

        # sweep
        delete $pid_InfoPool_Href->{$pid};
    }

    # reset
    %$pid_InfoPool_Href = ();
    %$only_has_new_pid_Href = ();
}

#--- merge bam and sort by readid ---
sub merge_bam_and_Nsort{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $bam_pref = $Option_Href->{bam_pref};
    my $source_bam_Aref = $Option_Href->{source_bam_Aref};
    my $prev_bam_Aref = $Option_Href->{prev_bam_Aref};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # merge
    # -h $source_bam, discard on 2016-05-03
    # -h FILE    Copy the header in FILE to <out.bam> [in1.bam], although mannual is this, but it combines all headers' SQ
    my $merge_command = "$V_Href->{samtools} merge -f -c -p $bam_pref.merged.bam @$source_bam_Aref @$prev_bam_Aref";
    trible_run_for_success($merge_command, 'MergeBam', {cmd_Nvb=>1, esdo_Nvb=>1});

    my $merge_nSortBam = $bam_pref.".merged.sortN.bam";
    my $nsort_command = "$V_Href->{samtools} sort --threads 4 -n -m 1G -T $merge_nSortBam.sort_tmp -o $merge_nSortBam $bam_pref.merged.bam";
    trible_run_for_success($nsort_command, 'nSortBam', {cmd_Nvb=>1, esdo_Nvb=>1});
    `rm -rf $bam_pref.merged.bam`;

    return $merge_nSortBam;
}

#--- count mapping reads based on the given v_refseg Href ---
sub count_virus_reads_in_bam{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $bam_for_reads_count = $Option_Href->{bam_for_reads_count};
    ## properly mapped, -f 0x2
    my $proper_map_bool = $Option_Href->{proper_map};
    ## not supplement alignment, -F 0x800 + 0x100
    my $non_supple_map_bool = $Option_Href->{non_supple};
    ## mapping quality >= 30, -q 30
    my $min_MAQ = $Option_Href->{min_MAQ};
    ## Softclip bases ratio, <= 0.00
    my $max_Sratio = $Option_Href->{max_Sratio};

    # empty bam
    if($bam_for_reads_count =~ /^(EMPTY|SWEEP|NA)$/){
        return 0;
    }

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # firstly count the standard reads
    my $reads_count = 0;
    my $proper_map_cmd = ($proper_map_bool) ? '-f 0x2' : '';
    my $non_supple_cmd = ($non_supple_map_bool) ? '-F 0x900' : ''; # 0x800 + 0x100
    open (BAM, "$V_Href->{samtools} view $proper_map_cmd $non_supple_cmd -q $min_MAQ $bam_for_reads_count |") || die "fail samtools view $bam_for_reads_count: $!\n";
    while(<BAM>){
        my ($cigar, $rseq) = (split)[5,9];
        my $rlen = length($rseq);
        my $S_sum = 0;
        while($cigar =~ /(\d+)S/){
            $S_sum += $1;
            $cigar =~ s/$1S//;
        }
        if($S_sum/$rlen <= $max_Sratio){
            $reads_count++;
        }
    }
    close BAM;

    return $reads_count;
}

#--- details about virus variant reconstruction work ---
sub virus_variant_reconstruct_pipe{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $virus_variant_NO = $Option_Href->{variant_NO};
    my $virus_refseg = $Option_Href->{var_name};
    my $MapOrigVrefBam = $Option_Href->{var_OrigBam};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $SampleID = $V_Href->{sample_ID};

    # info about this variant v_ref
    my $this_variant_Href = \%{ $Variant_Info_Href->{$virus_variant_NO} };

    # folder
    $this_variant_Href->{workfolder} = GetPath( filekey => 'single_virus_variant_folder', virus_variant_NO => $virus_variant_NO );
    `rm -rf $this_variant_Href->{workfolder}`;
    `mkdir -p $this_variant_Href->{workfolder}`;

    # aimmed files
    my $this_variant_vFa = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.fa" );
    my $this_variant_mkup_vFa = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.makeup.fa" );
    my $this_variant_mkdupBam = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.mkdup.bam" );
    my $this_variant_finalBam = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.bam" );
    my $this_variant_mods_file = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.modifications.list" );
    my $this_variant_New2Orig_PosMapping_file = File::Spec->catfile( $this_variant_Href->{workfolder}, "V$virus_variant_NO.regmap.list" );

    # record this variant info
    $this_variant_Href->{variant_NO} = $virus_variant_NO;
    $this_variant_Href->{OrigRefSeg} = $virus_refseg;
    ## for re-sort virus variants NO. based on DEPTH, we avoid variant_NO in refseg
    ## start to use new mut-Vrefseg name from 2015-05-04.
    # $this_variant_Href->{RefSeg} = "$SampleID-V$virus_variant_NO";
    $this_variant_Href->{RefSeg} = "$SampleID-$virus_refseg-RC"; # 'RC' means ReConstructed
    $this_variant_Href->{Orig_length} = $Virus_Info_Href->{$virus_refseg}->{length}; # updates later when confirm the re-construction
    $this_variant_Href->{NonN_length} = $Virus_Info_Href->{$virus_refseg}->{length}; # updates later when confirm the re-construction
    $this_variant_Href->{Extend_length} = $Virus_Info_Href->{$virus_refseg}->{extend_len};
    $this_variant_Href->{Cricle} = $Virus_Info_Href->{$virus_refseg}->{circle};
    $this_variant_Href->{Fa} = $this_variant_vFa;
    $this_variant_Href->{Fa_md5} = {};
    $this_variant_Href->{mkupFa} = $this_variant_mkup_vFa;
    $this_variant_Href->{mkdupBam} = $this_variant_mkdupBam;
    $this_variant_Href->{finalBam} = $this_variant_finalBam;
    $this_variant_Href->{VCF} = $this_variant_finalBam.'.vcf';
    $this_variant_Href->{CycleTime} = 0;
    $this_variant_Href->{Coverage} = {1=>0, 4=>0, 10=>0, 30=>0, 50=>0, 100=>0, $V_Href->{min_depth}=>0};
    $this_variant_Href->{AverageDepth} = 0;
    $this_variant_Href->{Mutation} = {};
    $this_variant_Href->{Mods} = $this_variant_mods_file;
    $this_variant_Href->{ModCompBool} = 0;
    $this_variant_Href->{PosMapList} = $this_variant_New2Orig_PosMapping_file;

    # to store the cycle work
    my %ConsCycle;

    #++++++++++++++++++++++++++++++++++++++++#
    # pretent the initial bam as the Cycle 0 #
    #++++++++++++++++++++++++++++++++++++++++#
    $ConsCycle{0}->{workFolder} = File::Spec->catfile( $this_variant_Href->{workfolder}, 'C0-From_origVref' );
    `rm -rf $ConsCycle{0}->{workFolder}`;
    `mkdir -p $ConsCycle{0}->{workFolder}`;
    ## the Cycle 0 v_ref for modification
    ## need to copy to prepare real file because if V1 is no-modification, we need to move these files
    ## note that the FA-header is not that related to this virus-variant.NO
    $ConsCycle{0}->{Orig_length} = $Virus_Info_Href->{$virus_refseg}->{length};
    $ConsCycle{0}->{NonN_length} = $Virus_Info_Href->{$virus_refseg}->{length};
    $ConsCycle{0}->{Fa} = File::Spec->catfile( $ConsCycle{0}->{workFolder}, 'C0.reconstruct.virus.fa' );
    $ConsCycle{0}->{makeupFa} = File::Spec->catfile( $ConsCycle{0}->{workFolder}, 'C0.reconstruct.virus.makeup.fa' );
    `cp -fL $Virus_Info_Href->{$virus_refseg}->{vref_fa} $ConsCycle{0}->{Fa}`;
    `cp -fL $Virus_Info_Href->{$virus_refseg}->{vref_fa} $ConsCycle{0}->{makeupFa}`;
    ## the Cycle 0 v_ref for Bwa alignment
    ## note that the FA-header is not that related to this virus-variant.NO
    $ConsCycle{0}->{vBwaRef} = $Virus_Info_Href->{$virus_refseg}->{vref_fa};
    ## prefix
    $ConsCycle{0}->{BamPref} = File::Spec->catfile( $ConsCycle{0}->{workFolder}, "V$virus_variant_NO.align_origVref" );
    ## link the original aligment bam to this Vxx work folder
    ## do not worry, for V1, this is just sofr-link created in function 'get_next_candidate_virus'
    ## Or, for other variant_NO, this is real file created in function 'get_next_candidate_virus'
    ## note that the SAM-header is not that related to this virus-variant.NO
    `mv -f $MapOrigVrefBam $ConsCycle{0}->{BamPref}.mkdup.bam`;
    `mv -f $MapOrigVrefBam.bai $ConsCycle{0}->{BamPref}.mkdup.bam.bai`;
    ## initialize the 'has_mut_sign' as 1 to avoid the first cycle use 'N' for low depth position
    $ConsCycle{0}->{HasMutSign} = 1;

    #++++++++++++++++++++++++++++++++++++++#
    # try to construct vref_fa iteratively #
    #++++++++++++++++++++++++++++++++++++++#
    # thoughts: initially do not use 'N' from the Cycle NO.1, and till
    # meets the cycle has not mutations, then allow 'N' in next cycle.
    my $Use_N_for_lowdepth_sign = 0;
    my $Cycle_time = 0;
    while ( ++$Cycle_time ){
        #------------------------------#
        # deal with the last cycle bam #
        #------------------------------#
        # last cycle time
        my $last_Cycle_time = $Cycle_time - 1;

        # GATK-realn-Samtools-calmd pipeline
        my $lastCycle_bam_for_vcf =
        FuseSV::Virus_Type::SmartAlignment->GATK_relan_Samtools_calmd({
            OutPref     => $ConsCycle{$last_Cycle_time}->{BamPref}, # auto add 'mkdup.bam'
            RefForBwa   => $ConsCycle{$last_Cycle_time}->{vBwaRef},
            SkipCalmd   => 1 # Calmd will introduce loss of MNP
        });

        # call vcf via SamTools mpileup
        my $lastCycle_vcf =
        FuseSV::Virus_Type::VCFwork->vcf_from_bam({
            Bam => $lastCycle_bam_for_vcf,
            OutPref => $ConsCycle{$last_Cycle_time}->{BamPref}, # auto add '.vcf'
            RefForBwa   => $ConsCycle{$last_Cycle_time}->{vBwaRef}
        });

        #-----------------------------#
        # this cycle work begins here #
        #-----------------------------#
        # whether use 'N' to replace the low-depth site in this cycle
        $Use_N_for_lowdepth_sign ||= (! $ConsCycle{$last_Cycle_time}->{HasMutSign});
        # inform
        stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO NO.$Cycle_time cycle of construction starts.\n"
                              ."\tIn this cycle, " . ($Use_N_for_lowdepth_sign?'will':'do not') . " use 'N' to replace the low-depth sites.\n";
        # workspace folder
        $ConsCycle{$Cycle_time}->{workFolder} = File::Spec->catfile( $this_variant_Href->{workfolder}, "C$Cycle_time-reconstruct" );
        `rm -rf $ConsCycle{$Cycle_time}->{workFolder}`;
        `mkdir -p $ConsCycle{$Cycle_time}->{workFolder}`;

        #------------------------------------------------#
        # read credible mutations from VCF of last cycle #
        #------------------------------------------------#
        # inform progress
        stout_and_sterr `date`."[INFO]:\tDeal with mutation from Bwa-SamTools-VCF cascade.\n";
        # filter mutations from vcf of last cycle
        ## hash stores each refpos information
        my $RefPos_Href = {};
        &filter_mutation_from_VCF({
            RefPos_Href => $RefPos_Href,
            VCF => $lastCycle_vcf,
            vBwaRef_Fa => $ConsCycle{$last_Cycle_time}->{vBwaRef},
            vMakeupFa => $ConsCycle{$last_Cycle_time}->{makeupFa},
            V_refseg => $this_variant_Href->{OrigRefSeg} #$virus_refseg,
        });

        #------------------------------------------------#
        # use blast to find back mutations missed by VCF #
        # just add these mutation(s) to %{$RefPos_Href}  #
        #  > especially for hyper-mutated blocks         #
        #------------------------------------------------#
        # inform progress
        stout_and_sterr `date`."[INFO]:\tTry find back missed mutations by blast cascade.\n";
        # prepare work_folder
        my $blast_works_folder = File::Spec->catfile( $ConsCycle{$Cycle_time}->{workFolder}, 'blast_MissMuSc' );
        &blast_find_back_mutations({
            vRefName => ($Cycle_time == 1) ? $this_variant_Href->{OrigRefSeg} : $this_variant_Href->{RefSeg},
            RefPos_Href => $RefPos_Href,
            orig_vExt => $Virus_Info_Href->{$virus_refseg}->{extend_len}, # to subtract the extend part, if possible
            lastCycle_bam => $lastCycle_bam_for_vcf,
            lastCycle_vMakeupFa => $ConsCycle{$last_Cycle_time}->{makeupFa},
            works_folder => $blast_works_folder
        });

        # sweep last cycle bam/bai files
        # unless( $V_Href->{in_debug} ){
        #   `rm -rf $lastCycle_bam_for_vcf $lastCycle_bam_for_vcf.bai`;
        # }

        #-------------------------------------------------#
        # modify vref_fa for this cycle base on mutations #
        #-------------------------------------------------#
        my $mody_VirusFa = File::Spec->catfile( $ConsCycle{$Cycle_time}->{workFolder}, "C$Cycle_time.reconstruct.virus.fa" );
        my $mody_makeup_VirusFa = File::Spec->catfile( $ConsCycle{$Cycle_time}->{workFolder}, "C$Cycle_time.reconstruct.virus.makeup.fa" );
        # INSTANCE: $modify_sign_Aref = [$ever_modify_sign, $has_mutation_sign]
        # my ($modify_sign_Aref, $alt_info_Href, $thisCycle_vRef_origLen, $thisCycle_vRef_nonNLen) = 
        my $reconstruct_details_Href = 
            FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->reconstruct_virus_genome_BaseOn_mutation({
                RefPos_Href => $RefPos_Href,
                Mod_V_G => $mody_VirusFa,
                Mod_V_G_refname => $this_variant_Href->{RefSeg},
                Mod_makeup_V_G => $mody_makeup_VirusFa,
                Use_N_for_lowdepth => $Use_N_for_lowdepth_sign
            });

        #----------------------------------------#
        # inform the modifications in this cycle #
        #----------------------------------------#
        &inform_modifications_in_one_cycle({
            modify_inform_str => '(BOOL:Mod[' . $reconstruct_details_Href->{ever_modify_sign} . ']' .
                                 ',Mut['      . $reconstruct_details_Href->{has_mutation_sign} . ']' .
                                 ',LowD['     . $reconstruct_details_Href->{has_low_depth_pos} . '])',
            # modify_sign_Aref => $modify_sign_Aref,
            alt_info_Href => $reconstruct_details_Href->{alt_info_Href},
            Cycle_time => $Cycle_time,
            virus_variant_NO => $virus_variant_NO
        });

        #---------------------------#
        # this cycle work ends here #
        stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO NO.$Cycle_time cycle of construction finishes.\n";
        #---------------------------#

        #----------------------------------------#
        # IF no modification in this cycle, stop #
        # ELSE, just prepare for next cycle ...  #
        #----------------------------------------#
        my $new_vRef_sign = 1;
        # if has modification, try whether comes before via MD5 check
        if( $reconstruct_details_Href->{ever_modify_sign} ){
            # calculate this cycle's variant genome md5 value
            chomp(my $this_cycle_vFa_md5 = `md5sum $mody_VirusFa | awk '{print \$1}'`);
            # find the same MD5 in pre-cycle(s)
            if( exists( $this_variant_Href->{Fa_md5}->{$this_cycle_vFa_md5} ) ){
                # inform
                my $last_md5_cycleTime = $this_variant_Href->{Fa_md5}->{$this_cycle_vFa_md5};
                stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO NO.$Cycle_time cycle [MD5] is identical to NO.$last_md5_cycleTime cycle.\n";
                # has low_depth pos to deal with 'N'
                if( $reconstruct_details_Href->{has_low_depth_pos} ){
                    # start to assign 'N' to low_depth pos
                    # and, donot need to record this md5 value
                    $Use_N_for_lowdepth_sign = 1;
                    # inform
                    stout_and_sterr "\tBut certain low-depth positions exist in this cycle, so allow 'N' for low-depth-pos NOW and Go On to next cycle.\n";
                }
                # donot have 'N' to deal with
                else{
                    # this is not new reference of this virus variant
                    $new_vRef_sign = 0;
                    # record the cycle time of this virus variant as the cycle_time where this MD5 firstly appears
                    $this_variant_Href->{CycleTime} = $last_md5_cycleTime;
                    # inform
                    stout_and_sterr "\tAnd no low-depth to modify any more, so use NO.$last_md5_cycleTime results.\n";
                }
            }
            # new MD5 value, record it with this cycle
            else{
                $this_variant_Href->{Fa_md5}->{$this_cycle_vFa_md5} = $Cycle_time;
            }
        }
        # no modification in this Cycle, stop here, and use results of last Cycle
        else{
            $new_vRef_sign = 0;
            # record the cycle time, this time no-change, so record the last time
            $this_variant_Href->{CycleTime} = $last_Cycle_time;
            # inform
            stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO NO.$Cycle_time cycle is identical to NO.$last_Cycle_time cycle.\n";
        }
        # no more new modification to virus ref, stop here, and arrange files
        if( !$new_vRef_sign ){
            # $this_variant_Href->{CycleTime} is just assigned above
            my $accept_Cycle_time = $this_variant_Href->{CycleTime};
            # update the refseg length of this virus variant
            $this_variant_Href->{Orig_length} = $ConsCycle{$accept_Cycle_time}->{Orig_length};
            $this_variant_Href->{NonN_length} = $ConsCycle{$accept_Cycle_time}->{NonN_length};
            # use files of accepted cycle
            ## not the very start (cycle 0), just move and copy
            if( $accept_Cycle_time != 0 ){
                `mv -f  $ConsCycle{$accept_Cycle_time}->{Fa} $this_variant_vFa`;
                `mv -f  $ConsCycle{$accept_Cycle_time}->{makeupFa} $this_variant_mkup_vFa`;
                `cp -fL $ConsCycle{$accept_Cycle_time}->{BamPref}.mkdup.bam $this_variant_mkdupBam`;
                `cp -fL $ConsCycle{$accept_Cycle_time}->{BamPref}.mkdup.bam.bai $this_variant_mkdupBam.bai`;
            }
            ## it is the very start (cycle 0), the FA-header and SAM-header need to update for this virus-variant.NO
            else{
                # reload the virus fa and subtract the extended part
                for my $FAs_Aref ( [ $ConsCycle{0}->{Fa}, $this_variant_vFa ],
                                   [ $ConsCycle{0}->{makeupFa}, $this_variant_mkup_vFa ]
                ){
                    FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
                        Ref => $FAs_Aref->[0],
                        OutFaFile => $FAs_Aref->[1],
                        SegName => $this_variant_Href->{RefSeg},
                        LineBase => $V_Href->{line_base},
                        # negative value means subtract such length seq from tail
                        CircleExtLen => -1 * $this_variant_Href->{Extend_length}
                    });
                }
                # new SAM-header for reheader
                `$V_Href->{samtools} view -H $ConsCycle{0}->{BamPref}.mkdup.bam | sed 's/\\bSN:$virus_refseg\\b/SN:$this_variant_Href->{RefSeg}/' > $ConsCycle{0}->{BamPref}.mkdup.bam.new.header.sam`;
                # reheader the finalBam
                `$V_Href->{samtools} reheader $ConsCycle{0}->{BamPref}.mkdup.bam.new.header.sam $ConsCycle{0}->{BamPref}.mkdup.bam > $this_variant_mkdupBam`;
                `$V_Href->{samtools} index $this_variant_mkdupBam`;
            }
            # empty all cycles workspace of this virus variant
            unless( $V_Href->{in_debug} ){
                `rm -rf $ConsCycle{$_}->{workFolder}` for ( 0 .. $Cycle_time );
            }
            last;
        }
        # donot sweep last_cycle workfolder, empty them one by one at last.
        # unless( $V_Href->{in_debug} ){
        #   `rm -rf $ConsCycle{$last_Cycle_time}->{workFolder}`;
        # }

        #----------------------------------------------------------------#
        # construct indexed reference for NEW virus genome of this cycle #
        #----------------------------------------------------------------#
        # prepare fa file for indexing
        my $Indexed_VrefFa = File::Spec->catfile( $ConsCycle{$Cycle_time}->{workFolder}, "C$Cycle_time.reconstruct.virus.for.index.fa" );
        # specially extend for circular virus
        if( $this_variant_Href->{Cricle} =~ /yes/i ){
            # reload the virus fa and extend
            FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
                Ref => $mody_VirusFa,
                OutFaFile => $Indexed_VrefFa,
                SegName => $this_variant_Href->{RefSeg},
                LineBase => $V_Href->{line_base},
                CircleExtLen => $this_variant_Href->{Extend_length}
            });
        }
        else{ # just soft link
            `ln -sf $mody_VirusFa $Indexed_VrefFa`;
        }
        # create the indexed reference of new constructed vref_fa
        BWA_index_fasta(FaFile => $Indexed_VrefFa, bwa => $V_Href->{bwa});
        Faidx_Dict_fasta(FaFile => $Indexed_VrefFa, samtools => $V_Href->{samtools});

        #-------------------------------------------------------#
        # align orginal reads to new virus genome of this cycle #
        #-------------------------------------------------------#
        # use orginal smart fq to align this new reference
        my $mody_Virus_BamPref = File::Spec->catfile( $ConsCycle{$Cycle_time}->{workFolder}, "C$Cycle_time.reconstruct.SmartFqMap" );
        # here, deal with alignment against single virus reference, multiple forks could be applied here.
        my $mody_Virus_mkdup_bam = 
        FuseSV::Virus_Type::SmartAlignment->SmartFq2MergeBam_against_given_ref({
            RefForBwa   =>  $Indexed_VrefFa,
            OutBamPref  =>  $mody_Virus_BamPref,
            fork_number =>  $V_Href->{fork_number}
        });
        # check the mkdup_bam
        # if this virus variant has no alignment at this cycle.
        if( $mody_Virus_mkdup_bam eq 'EMPTY' ){
            # delete all data of this virus variant, note it is not one cycle, it is whole data of this virus variant
            unless( $V_Href->{in_debug} ){
                `rm -rf $this_variant_Href->{workfolder}`;
            }
            delete $Variant_Info_Href->{$virus_variant_NO};
            # inform
            stout_and_sterr `date`."[WARN]:\tvariant V$virus_variant_NO lacks alignment at Cycle NO.$Cycle_time. "
                                  ."\tSo stops here.\n";
            return;
        }

        #-----------------------------------------------------------------#
        # select best alignment of variants for reads from 2nd virus type #
        #-----------------------------------------------------------------#
        # merge with previous detected virus variant and extract reads for this virus variant
        if($virus_variant_NO != 1){
            ## prepare information of the former variant(s)
            my @prev_vbam;
            my %prev_v_refseg;
            for my $prev_vNO (sort {$a<=>$b} keys %$Variant_Info_Href){
                next if($prev_vNO == $virus_variant_NO); # this variant_NO recorded just above
                my $bam = $Variant_Info_Href->{$prev_vNO}->{finalBam};
                my $refseg = $Variant_Info_Href->{$prev_vNO}->{RefSeg};
                push @prev_vbam, $bam;
                $prev_v_refseg{$refseg} = $prev_vNO;
            }
            ## merge bam and sort by read name
            my $merge_nSortBam = 
            &merge_bam_and_Nsort({
                bam_pref => $mody_Virus_BamPref,
                source_bam_Aref => [ $mody_Virus_BamPref.".mkdup.bam" ],
                prev_bam_Aref => \@prev_vbam
            });
            ## filter the Nsort bam and output the remained bam
            &filter_Nsort_bam_and_to_next_variant_bam({
                nSortBam => $merge_nSortBam,
                OutBam => $mody_Virus_BamPref.".mkdup.bam",
                BwaRef => $Indexed_VrefFa,
                Prev_vref_Href => \%prev_v_refseg,
                keep_SQ => $this_variant_Href->{RefSeg}
            });
            # sweep
            `rm -rf $merge_nSortBam`;
            # inform
            stout_and_sterr `date`."[INFO]:\tmerge variant V$virus_variant_NO Cycle NO.$Cycle_time with pre-confirmed virus variants.\n";
        }

        #----------------------------------------#
        # IF it's the maximum cycle time, stop   #
        #----------------------------------------#
        # once the cycle time reaches the maximum, stop
        if($Cycle_time == $V_Href->{max_cycle_time}){
            # record the cycle time
            $this_variant_Href->{CycleTime} = $Cycle_time;
            # update the refseg length of this virus variant
            $this_variant_Href->{Orig_length} = $reconstruct_details_Href->{orig_length};
            $this_variant_Href->{NonN_length} = $reconstruct_details_Href->{nonN_length};
            # inform
            stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO reaches the maximum construction cycle time ($V_Href->{max_cycle_time}), so stop here.\n"
                                  ."<WARN>:\tGenerally, reaching the maximum time is not common. You'd better check the result, OR contact with the author.\n";
            # use this_cycle files
            `mv -f $mody_VirusFa $this_variant_vFa`;
            `mv -f $mody_makeup_VirusFa $this_variant_mkup_vFa`;
            `mv -f $mody_Virus_BamPref.mkdup.bam $this_variant_mkdupBam`;
            `mv -f $mody_Virus_BamPref.mkdup.bam.bai $this_variant_mkdupBam.bai`;
            # empty all cycles workspace of this virus variant
            unless( $V_Href->{in_debug} ){
                `rm -rf $ConsCycle{$_}->{workFolder}` for ( 0 .. $Cycle_time );
            }
            last;
        }

        # record remained information
        # {workFolder} has been assigned at the very start of this cycle work
        $ConsCycle{$Cycle_time}->{HasMutSign} = $reconstruct_details_Href->{has_mutation_sign};
        $ConsCycle{$Cycle_time}->{Fa} = $mody_VirusFa;
        $ConsCycle{$Cycle_time}->{makeupFa} = $mody_makeup_VirusFa;
        $ConsCycle{$Cycle_time}->{Orig_length} = $reconstruct_details_Href->{orig_length};
        $ConsCycle{$Cycle_time}->{NonN_length} = $reconstruct_details_Href->{nonN_length};
        $ConsCycle{$Cycle_time}->{vBwaRef} = $Indexed_VrefFa;
        $ConsCycle{$Cycle_time}->{BamPref} = $mody_Virus_BamPref;
    }

    # at last, generate finalBam from mkdupBam
    ## of course, should make indexed virus-ref of the last mod-vref
    &get_variant_finalBam_from_mkdupBam({
        variant_Href => $this_variant_Href
    });

    # stat depth and coverage on this variant alignment bam
    ## here is for the judgment after this variant to go on re-construction or not
    ## here does not need details of coverage
    ## note that it will operate again after 're-arrange reads of each virus variant' to get the REAL information
    &get_depth_coverage_given_bam({
        Bam => $this_variant_Href->{finalBam},
        check_seg => $this_variant_Href->{RefSeg},
        seg_orig_length => $this_variant_Href->{Orig_length},
        seg_nonN_length => $this_variant_Href->{NonN_length},
        seg_ext_length => $this_variant_Href->{Extend_length},
        need_details => 1,
        coverage_Href => $this_variant_Href->{Coverage}, # used when need_details is 1
        averagedepth_Sref => \$this_variant_Href->{AverageDepth} # used when need_details is 1
    });
}

#--- filter VCF and get modified virus fa ---
sub filter_mutation_from_VCF{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $RefPos_Href = $Option_Href->{RefPos_Href};
    my $vcf = $Option_Href->{VCF};
    my $vBwaRef_fa = $Option_Href->{vBwaRef_Fa};
    my $makeup_vFa = $Option_Href->{vMakeupFa} || undef;
    my $v_refseg = $Option_Href->{V_refseg};
    my $JustLoad_NoFilter = $Option_Href->{NoFilter} || 0;

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

    # virus attributes
    my $v_circle = $Virus_Info_Href->{$v_refseg}->{circle};
    my $v_ExtLen = $Virus_Info_Href->{$v_refseg}->{extend_len};

    # load the virus reference
    FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
        Ref => $vBwaRef_fa,
        makeup_Ref => $makeup_vFa, # specific because makeup ref seq requires non-N
        RefPos_Href => $RefPos_Href
    });

    # read VCF and load each base info
    FuseSV::Virus_Type::VCFwork->load_VCF_and_refpos_info({
        RefPos_Href => $RefPos_Href,
        VCF => $vcf
    });

    # merge bases in the extended part
    if($v_circle eq 'yes'){
        FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->merge_circular_virus_extend_part({
            RefPos_Href => $RefPos_Href,
            Ext_Len => $v_ExtLen
        });
    }

    # skip the mutation filtration
    return if($JustLoad_NoFilter);

    # filter refpos and mutations from VCF
    FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->filter_refpos_mutation_from_VCF({
        RefPos_Href => $RefPos_Href
    });
}

#--- use blast to find back the mutations missed by VCF ---
#--- specially, the hyper-mutated blocks ---
sub blast_find_back_mutations{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $Vrefname = $Option_Href->{vRefName};
    my $RefPos_Href = $Option_Href->{RefPos_Href};
    my $lastCycle_bam = $Option_Href->{lastCycle_bam};
    my $lastCycle_vMakeupFa = $Option_Href->{lastCycle_vMakeupFa};
    my $orig_vExt = $Option_Href->{orig_vExt};
    my $works_folder = $Option_Href->{works_folder};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    my $CIGAR_Blast_Mut_Min_Dist = 3;

    # read bam to get the INDEL and SoftClip sites
    my %Mut_SoftClip_Sites; # stores these mutations
    my %SoftClip_Sites_Sup; # just store soft-clip info
    &extract_InDel_SoftClip_Sites_from_BAM(
        {
            RefPos_Href => $RefPos_Href,
            lastCycle_bam => $lastCycle_bam,
            lastCycle_vMakeupFa => $lastCycle_vMakeupFa,
            orig_vExt => $orig_vExt,
            max_vPos => max(keys %$RefPos_Href),
            Vrefname => $Vrefname,
            MuSc_Sites_Href => \%Mut_SoftClip_Sites,
            SoftClip_Href => \%SoftClip_Sites_Sup,
            # SoftClip_Min_Size => 10,
            MuSc_Min_Support => $V_Href->{min_mut_reads}
        }
    );
    # check mutations
    if( scalar(keys %Mut_SoftClip_Sites) == 0 ){
        stout_and_sterr `date`."[INFO]:\tNo mutated sites picked for blast find-back operation.\n";
        return;
    }
    # inform
    stout_and_sterr `date`."[INFO]:\tPick mutated sites for blast find-back operation OK.\n";

    # extract the reads for blastall
    `rm -rf $works_folder`;
    `mkdir -p $works_folder`;
    my $read_fa_file = File::Spec->catfile( $works_folder, 'MissMuSc.reads.fa' );
    &extract_reads_from_mut_region(
        {
            Vrefname => $Vrefname,
            read_fa_file => $read_fa_file,
            lastCycle_bam => $lastCycle_bam,
            MuSc_Sites_Href => \%Mut_SoftClip_Sites
        }
    );
    # inform
    stout_and_sterr `date`."[INFO]:\tExtract reads from mutated sites for blast find-back operation OK.\n";

    # formatdb and blastall
    ## construct blast database(reference)
    my $blastall_ref = File::Spec->catfile( $works_folder, 'lastVref.for_MissMuSc_blast.fa' );
    `ln -sf $lastCycle_vMakeupFa $blastall_ref`;
    FuseSV::Virus_Type::BlastM0Mutation->format_database( { ref_fa_file => $blastall_ref } );
    ## blastall query
    my $blast_m0 = File::Spec->catfile( $works_folder, 'MissMuSc.reads.blastall.m0' );
    FuseSV::Virus_Type::BlastM0Mutation->blastall_query_m0(
        {
            ref_fa_file => $blastall_ref,
            query_fa_file => $read_fa_file,
            blast_m0 => $blast_m0
        }
    );
    # inform
    stout_and_sterr `date`."[INFO]:\tBlast find-back mutations blastall OK.\n";

    # upload mutations(modifications)
    ## in this function, merge the extended part with the original part.
    my %FindBack_Mutation;
    FuseSV::Virus_Type::BlastM0Mutation->read_blast_m0_for_mutation(
        {
            blast_m0 => $blast_m0,
            # query_seq_fa => $read_fa_file, # skip as reads maybe too much to handle memory
            orig_fa_ref => $blastall_ref,
            orig_ext => $orig_vExt,
            mutation_Href => \%FindBack_Mutation,
            min_depth => $V_Href->{min_depth},
            both_strd => $V_Href->{both_strd},
            min_mut_reads => $V_Href->{min_mut_reads}
        }
    );
    # inform
    stout_and_sterr `date`."[INFO]:\tRead m0 result for blast find-back operation OK.\n";

    # pick 'Blast' find-back mutation(s) back to mutation pool
    ## filter this site if it already has mutation
    for my $mut_pos (sort {$a<=>$b} keys %FindBack_Mutation){
        next if( exists($RefPos_Href->{$mut_pos}) && exists($RefPos_Href->{$mut_pos}->{mutation}));
        # my ($mut_type, $mut_seq, $sum_sup, $plus_sup, $minus_sup) = @{ $FindBack_Mutation{$mut_pos} };
        my $mut_type  = $FindBack_Mutation{$mut_pos}->{mut_type};
        my $mut_seq   = $FindBack_Mutation{$mut_pos}->{mut_seq};
        my $sum_sup   = $FindBack_Mutation{$mut_pos}->{mut_sum_sup};
        my $plus_sup  = $FindBack_Mutation{$mut_pos}->{mut_plus_sup};
        my $minus_sup = $FindBack_Mutation{$mut_pos}->{mut_minus_sup};
        # load 'Blast' find-back mutation into the credible mutation pool
        $RefPos_Href->{$mut_pos}->{mutation} = {
                                                 Sum_Sup => $sum_sup,
                                                 Plus_Sup => $plus_sup,
                                                 Minus_Sup => $minus_sup,
                                                 Mut_Type => $mut_type,
                                                 Mut_Seq => $mut_seq,
                                                 Method => 'Blast'
                                               };
        # inform
        stout_and_sterr "[INFO]:\taccept Blast mutation: $mut_pos ".
                        join(',', ($mut_type, $mut_seq, $sum_sup, $plus_sup, $minus_sup)) . "\n";
        # remove ref_allele counts given by 'insertion' specially
        if( $mut_type eq 'ins' ){
            $RefPos_Href->{$mut_pos}->{refallel_depth} -= $sum_sup;
            $RefPos_Href->{$mut_pos}->{refallel_depth} = max( $RefPos_Href->{$mut_pos}->{refallel_depth}, 0 );
        }
        # delete the mutations obtained by 'CIGAR-analysis' nearby blast-mutations
        my $cigar_MutPos_exclude_Reg = [ $mut_pos-$CIGAR_Blast_Mut_Min_Dist+1, $mut_pos+$CIGAR_Blast_Mut_Min_Dist-1 ];
        for my $exclu_pos ( $cigar_MutPos_exclude_Reg->[0] .. $cigar_MutPos_exclude_Reg->[1] ){
            if( exists($Mut_SoftClip_Sites{$exclu_pos}) ){
                if( $V_Href->{in_debug} && exists($Mut_SoftClip_Sites{$exclu_pos}->{mutation}) ){
                    warn "discard CIGAR mutation: $exclu_pos:" . join(',',@{$Mut_SoftClip_Sites{$exclu_pos}->{mutation}}) . "\n";
                }
                delete $Mut_SoftClip_Sites{$exclu_pos};
            }
        }
    }

    # pick 'CIGAR' find-back mutation(s) back to mutation pool
    for my $mut_pos (sort {$a<=>$b} keys %Mut_SoftClip_Sites){
        if(    ! exists($RefPos_Href->{$mut_pos}) # ignore extended part, if possible
            || ! exists($Mut_SoftClip_Sites{$mut_pos}->{mutation}) # no CIGAR mutation on this site
            ||   exists($RefPos_Href->{$mut_pos}->{mutation}) # already has mutation record on this site
          ){
            if($V_Href->{in_debug}){
                warn "Skip CIGAR mutation site: $mut_pos, basic filtration.\n";
            }
            next;
        }
        my ($mut_type, $mut_seq, $sum_sup, $plus_sup, $minus_sup) = @{ $Mut_SoftClip_Sites{$mut_pos}->{mutation} };
        # compare with refallele-depth
        my $refallel_depth_for_Comp = $RefPos_Href->{$mut_pos}->{refallel_depth};
        # remove ref_allele counts given by 'insertion' specially
        if( $mut_type eq 'ins' ){
            $refallel_depth_for_Comp -= $sum_sup;
            if( exists($SoftClip_Sites_Sup{$mut_pos}->{softclip3}) ){
                my $softclip_3p_depth = $SoftClip_Sites_Sup{$mut_pos}->{softclip3};
                $refallel_depth_for_Comp -= $softclip_3p_depth;
                if($V_Href->{in_debug}){
                    warn " subtract softclip induced depth ($softclip_3p_depth) from ref_allele_depth. [arbitrary action]\n";
                }
            }
        }
        if( $sum_sup < $refallel_depth_for_Comp ){
            if($V_Href->{in_debug}){
                warn "discard CIGAR mutation site: $mut_pos, Depth ($sum_sup) is lower than "
                    ."refallel_depth ($RefPos_Href->{$mut_pos}->{refallel_depth}, $refallel_depth_for_Comp).\n";
            }
            next;
        }
        # load mutation find back into the credible mutation pool
        $RefPos_Href->{$mut_pos}->{mutation} = {
                                                 Sum_Sup => $sum_sup,
                                                 Plus_Sup => $plus_sup,
                                                 Minus_Sup => $minus_sup,
                                                 Mut_Type => $mut_type,
                                                 Mut_Seq => $mut_seq,
                                                 Method => 'CIGAR'
                                               };
        # update refallel_depth, 'ins'(-$sum_sup), 'del'(no_change)
        $RefPos_Href->{$mut_pos}->{refallel_depth} = max( $refallel_depth_for_Comp, 0 );
        # inform
        stout_and_sterr "[INFO]:\taccept CIGAR mutation: $mut_pos " . join(',', @{ $Mut_SoftClip_Sites{$mut_pos}->{mutation} }) . "\n";
    }
}

#--- extract undetected InDel SoftClip sites from Bam ---
sub extract_InDel_SoftClip_Sites_from_BAM{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $RefPos_Href = $Option_Href->{RefPos_Href};
    my $lastCycle_bam = $Option_Href->{lastCycle_bam};
    my $lastCycle_vMakeupFa = $Option_Href->{lastCycle_vMakeupFa};
    my $Vrefname = $Option_Href->{Vrefname};
    my $orig_vExt = $Option_Href->{orig_vExt};
    my $max_vPos = $Option_Href->{max_vPos};
    my $MuSc_Sites_Href = $Option_Href->{MuSc_Sites_Href};
    my $SoftClip_Href = $Option_Href->{SoftClip_Href};
    # my $SoftClip_Min_Size = $Option_Href->{SoftClip_Min_Size};
    my $MuSc_Min_Support = $Option_Href->{MuSc_Min_Support}; # this is for mutation sites selection for blast re-alignment

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $both_strd = $V_Href->{both_strd};
    # for mutations from CIGAR-mutation analysis
    my $min_mut_reads = $V_Href->{min_mut_reads};
    my $min_mut_tailD = $V_Href->{I16_min_mean_tailDist};
    my $RepeatMutFilterFlankSize = $V_Href->{cigar_RepeatMut_FilterFlankSize};
    my $RepeatMutFilterAlertCount = $V_Href->{cigar_RepeatMut_FilterAlertCount};

    # seek InDel and Softclip sites
    open (BAM, "$V_Href->{samtools} view -F 0x400 $lastCycle_bam |") || die "cannot read $lastCycle_bam: $!\n";
    while (<BAM>){
        chomp;
        my ($FLAG, $mPos, $CIGAR, $ReadSeq) = (split /\t+/)[1, 3, 5, 9];
        # my $read_len = length($ReadSeq);
        # filter
        if(    $FLAG & 0x400 # duplication
            || $FLAG & 0x100 # secondary alignment
            || $FLAG & 0x800 # supplementary alignment
            || $FLAG & 0x4   # unmapped
            || $CIGAR !~ /[IDS]/ # have InDel or SoftClip
            || $CIGAR =~ /[^IDSM\d]/ # have other cigar operators
          ){
            next;
        }
        # deal mutations
        my $i = 0;
        while(length($CIGAR)!=0){
            my ($mut_pos, $mut_type, $mut_seq);
            my ($len, $type) = ($CIGAR =~ /^(\d+)(\D)/);
            $CIGAR =~ s/^\d+\D//;
            if($type eq 'M'){
                $i += $len;
                $ReadSeq = substr($ReadSeq, $len);
                next;
            }
            elsif($type eq 'I'){
                $mut_pos = $mPos + $i - 1; # assign insertion to the forefront-pos
                $mut_type = 'ins';
                $mut_seq = uc(substr($ReadSeq, 0, $len));
                $ReadSeq = substr($ReadSeq, $len);
                next if( $i < $min_mut_tailD || length($ReadSeq) < $min_mut_tailD); # located at tail of this read
            }
            elsif($type eq 'D'){
                $mut_pos = $mPos + $i; # assign deletion to the first pos of the deleted part
                $i += $len;
                $mut_type = 'del';
                $mut_seq = 'N'; # assign the seq later
                next if( $i < $min_mut_tailD || length($ReadSeq) < $min_mut_tailD); # located at tail of this read
            }
            elsif($type eq 'S'){
                if( $i == 0 ){ # left prime of read
                    $mut_pos = $mPos;
                    $mut_type = 'softclip5';
                }
                else{ # right prime of read
                    $mut_pos = $mPos + $i - 1;
                    $mut_type = 'softclip3';
                }
                $i += $len;
                $ReadSeq = substr($ReadSeq, $len);
            }
            # filter this site if it already has mutation
            next if( exists($RefPos_Href->{$mut_pos}) && exists($RefPos_Href->{$mut_pos}->{mutation}));
            # # softclip length limitation
            # next if( $type eq 'S' && $len < $SoftClip_Min_Size);
            # record the mutated pos with supporting
            $MuSc_Sites_Href->{$mut_pos}->{mut_whole_depth} ++;
            # store supporting of soft-clip for later subtraction from ref_allele_depth
            if( $type eq 'S' ){
                $SoftClip_Href->{$mut_pos}->{$mut_type}++ if( $mut_type eq 'softclip3' ); # only primer3 impact depth
                next; # do not record softclip as mutation
            }
            # store insertion and deletion info
            my $mut_id = "$mut_type,$len,$mut_seq";
            unless( exists($MuSc_Sites_Href->{$mut_pos}->{mutation}->{$mut_id}) ){
                $MuSc_Sites_Href->{$mut_pos}->{mutation}->{$mut_id} = [0, 0, 0];
            }
            my $mut_Aref = $MuSc_Sites_Href->{$mut_pos}->{mutation}->{$mut_id};
            $mut_Aref->[0] ++; # all support
            if( $FLAG & 0x10 ){ # reverse complemented support
                $mut_Aref->[2] ++;
            }
            else{ # forward support
                $mut_Aref->[1] ++;
            }
        }
    }
    close BAM;

    # specially merge info in extended region for circular virus
    for my $rel_pos (1 .. $orig_vExt){
        my $sor_RefPos = $max_vPos + $rel_pos;
        # merge once the extended pos has CIGAR mutation
        if( exists($MuSc_Sites_Href->{$sor_RefPos}) ){
            unless( exists($MuSc_Sites_Href->{$rel_pos}) ){ # initialize if needed
                $MuSc_Sites_Href->{$rel_pos}->{mut_whole_depth} = 0;
                $MuSc_Sites_Href->{$rel_pos}->{mutation} = {};
            }
            my $aim_RefPos_Href = $MuSc_Sites_Href->{$rel_pos};
            my $sor_RefPos_Href = $MuSc_Sites_Href->{$sor_RefPos};
            # merge mutation whole depth
            $aim_RefPos_Href->{mut_whole_depth} += $sor_RefPos_Href->{mut_whole_depth};
            # merge each mutation
            for my $mut_id (sort keys %{$sor_RefPos_Href->{mutation}}){
                if(!exists($aim_RefPos_Href->{mutation}->{$mut_id})){
                    # FORMAT (Sum, Plus, Minus)
                    $aim_RefPos_Href->{mutation}->{$mut_id} = [ 0, 0, 0 ];
                }
                # merge
                for my $i ( 0 .. scalar(@{$sor_RefPos_Href->{mutation}->{$mut_id}})-1 ){
                    $aim_RefPos_Href->{mutation}->{$mut_id}->[$i] += $sor_RefPos_Href->{mutation}->{$mut_id}->[$i];
                }
            }
            # sweep the extended pos
            delete $MuSc_Sites_Href->{$sor_RefPos};
        }
    }

    # find and discard tandem-appeared mutations
    ## induced by indel at repeat region
    my %PosRepMutToDel;
    for my $mut_pos (sort {$a<=>$b} keys %$MuSc_Sites_Href){
        my $mut_Href = $MuSc_Sites_Href->{$mut_pos}->{mutation};
        for my $mut_id (keys %$mut_Href){
            # both sides to count same as this mutation
            my @UpHasSameMutPos = grep    exists($MuSc_Sites_Href->{$_})
                                       && exists($MuSc_Sites_Href->{$_}->{mutation}->{$mut_id}), ( $mut_pos-$RepeatMutFilterFlankSize .. $mut_pos );
            my @DwHasSameMutPos = grep    exists($MuSc_Sites_Href->{$_})
                                       && exists($MuSc_Sites_Href->{$_}->{mutation}->{$mut_id}), ( $mut_pos .. $mut_pos+$RepeatMutFilterFlankSize );
            # if reach alert amount, record to delete
            for my $Pos_Aref (\@UpHasSameMutPos, \@DwHasSameMutPos){
                if( scalar(@$Pos_Aref) >= $RepeatMutFilterAlertCount ){
                    $PosRepMutToDel{$_}->{$mut_id} = 1 for @$Pos_Aref;
                }
            }
        }
    }
    # delete crazy tandem-appeared mutations
    for my $mut_pos (sort {$a<=>$b} keys %PosRepMutToDel){
        for my $mut_id (keys %{$PosRepMutToDel{$mut_pos}}){
            # discard
            delete $MuSc_Sites_Href->{$mut_pos}->{mutation}->{$mut_id};
            # warn
            if( $V_Href->{in_debug} ){
                warn "Discard CIGAR-mutation at $Vrefname:$mut_pos ($mut_id) as tandem-appeared.\n";
            }
        }
    }

    # filter mutation supports
    ## assign mut_seq for the deltetion
    for my $mut_pos (sort {$a<=>$b} keys %$MuSc_Sites_Href){
        # depth filter
        if( $MuSc_Sites_Href->{$mut_pos}->{mut_whole_depth} < $MuSc_Min_Support ){
            delete $MuSc_Sites_Href->{$mut_pos};
            next;
        }
        # filter mutation
        ## assign mut_seq for the deletion
        my @mut_id;
        my $mut_Href = $MuSc_Sites_Href->{$mut_pos}->{mutation};
        for my $mut_id (keys %$mut_Href){
            my $mut_Aref = $mut_Href->{$mut_id};
            # FORMAT (Sum, Plus, Minus)
            my ($sum_sup, $plus_sup, $minus_sup) = @$mut_Aref;
            # filter reads supports
            if(    ( $sum_sup < $min_mut_reads )
                || ( $both_strd && ($plus_sup < $both_strd || $minus_sup < $both_strd) ) # if required, both strand supported
              ){
                delete $mut_Href->{$mut_id};
                next;
            }
            # update deletions
            my ($mut_type, $len, $mut_seq) = split /,/, $mut_id;
            if($mut_type eq 'del'){
                my $faidx_reg = "$Vrefname:$mut_pos-".($mut_pos+$len-1);
                ($mut_seq = uc(`$V_Href->{samtools} faidx $lastCycle_vMakeupFa $faidx_reg 2>/dev/null | sed -n '2,\$p'`)) =~ s/\s//g;
                my $new_mut_id = "$mut_type,$len,$mut_seq";
                # update mutation information
                delete $mut_Href->{$mut_id};
                $mut_Href->{$new_mut_id} = [ $sum_sup, $plus_sup, $minus_sup ];
                # record new_mut_id for best mutation selection
                push @mut_id, $new_mut_id;
            }
            else{
                # record mut_id for best mutation selection
                push @mut_id, $mut_id;
            }
        }
        # no mutations
        if(scalar(@mut_id) == 0){
            # deprive mutation container of this site(position)
            ## cannot delete this site, because later provide it to reads pick-up for blast re-alignment
            delete $MuSc_Sites_Href->{$mut_pos}->{mutation};
        }
        else{ # has mutations
            # sort by the sum_sup
            @mut_id = sort { $mut_Href->{$b}->[0] <=> $mut_Href->{$a}->[0] } @mut_id;
            my $best_mut_id = $mut_id[0];
            my ($mut_type, $len, $mut_seq) = split /,/, $best_mut_id; # first one
            my ($sum_sup, $plus_sup, $minus_sup) = @{ $mut_Href->{$best_mut_id} };
            # only record best mutation
            $MuSc_Sites_Href->{$mut_pos}->{mutation} = [ $mut_type, $mut_seq, $sum_sup, $plus_sup, $minus_sup ];
        }
    }

    # debug
    if( $V_Href->{in_debug} ){
        for my $mut_pos (sort {$a<=>$b} keys %$MuSc_Sites_Href){
            warn "CIGAR $mut_pos: " .
                 ((exists($MuSc_Sites_Href->{$mut_pos}->{mutation})) ? join(',', @{$MuSc_Sites_Href->{$mut_pos}->{mutation}}) : '_Filtered_') .
                 "\n";
        }
    }
}

#--- extract reads related to mut-sites from bam ---
sub extract_reads_from_mut_region{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $Vrefname = $Option_Href->{Vrefname};
    my $read_fa_file = $Option_Href->{read_fa_file};
    my $lastCycle_bam = $Option_Href->{lastCycle_bam};
    my $MuSc_Sites_Href = $Option_Href->{MuSc_Sites_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # group site into regions
    my @check_sites;
    my $site_gap = 5;
    my @sort_sites = sort {$a<=>$b} keys %$MuSc_Sites_Href;
    my $left_i = 0;
    for (my $i = 1; $i < scalar(@sort_sites); $i++){
        if( $sort_sites[$i] > $sort_sites[$i-1] + $site_gap){
            if($V_Href->{in_debug}){
                warn "Group these sites into one region: ".$sort_sites[$left_i]."-".$sort_sites[$i-1]."\n".
                                "\t" . join(',', @sort_sites[$left_i .. $i-1]) . "\n";
            }
            push @check_sites, "$Vrefname:".$sort_sites[$left_i]."-".$sort_sites[$i-1];
            $left_i = $i;
        }
    }
    push @check_sites, "$Vrefname:".$sort_sites[$left_i]."-".$sort_sites[-1];

    # read region bam part by part
    my $max_check_sites_count = 10;
    my %Reads;
    for (my $i = 0; $i <= $#check_sites; $i += $max_check_sites_count){
        my $last_site_i = min($i+$max_check_sites_count-1, $#check_sites);
        my $region_para = join(' ', @check_sites[$i .. $last_site_i]);
        if($V_Href->{in_debug}){
            warn "blast-bamview-region_para: $region_para\n";
        }
        # get ALL reads from InDel and Softclip sites
        open (BAM, "$V_Href->{samtools} view $lastCycle_bam $region_para |") || die "cannot read $lastCycle_bam: $!\n";
        while (<BAM>){
            chomp;
            my ($PairID, $FLAG, $ReadSeq) = (split /\t+/)[0, 1, 9];
            # filter
            if(    $FLAG & 0x400 # duplication
                || $FLAG & 0x100 # secondary alignment
                || $FLAG & 0x800 # supplementary alignment
                || $FLAG & 0x4   # unmapped
              ){
                next;
            }
            my $end_NO = ($FLAG & 0x40) ? 1 : 2;
            my $ReadID = $PairID.'/'.$end_NO;
            if( exists($Reads{$ReadID}) ){ # already load this read
                next;
            }
            # record this read
            if( $FLAG & 0x10 ){ # reverse complemented aligned
                ( $ReadSeq = reverse uc($ReadSeq) ) =~ tr/ACGT/TGCA/;
            }
            $Reads{$ReadID} = $ReadSeq;
        }
        close BAM;
    }

    # output the mutated reads
    open (FA, ">$read_fa_file") || die "fail write $read_fa_file: $!\n";
    print FA ">$_\n$Reads{$_}\n" for sort keys %Reads;
    close FA;
}

#--- just inform the modifications in this cycle ---
sub inform_modifications_in_one_cycle{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $modify_inform_str = $Option_Href->{modify_inform_str};
    my $alt_info_Href = $Option_Href->{alt_info_Href};
    my $Cycle_time = $Option_Href->{Cycle_time};
    my $virus_variant_NO = $Option_Href->{virus_variant_NO};

    # inform the modifications counts
    my $sum_alt_count = sum(values %{$alt_info_Href->{sum}});
    my $alt_count_str = "sum:$sum_alt_count";
    $alt_count_str .= '; '.$_.':'.$alt_info_Href->{sum}->{$_} for ('snp', 'ins', 'del');
    $alt_count_str .= '; '.$_.':'.$alt_info_Href->{N_count} for ('N_count');
    stout_and_sterr `date`."[INFO]:\tvariant V$virus_variant_NO NO.$Cycle_time cycle modifications"
                          ."$modify_inform_str" . ':'
                          ."$alt_count_str\n";
    # inform the modifications details, ONLY mutations
    if( $sum_alt_count ){
        # theme
        my @theme = qw/ Position Ref_Allel Mut_Type Mut_Seq Ref_Allel_Sup Mut_Sum_Sup Mut_Plus_Sup Mut_Minus_Sup Method /;
        my @tmlen = map length($_), @theme;
        # meet the widest
        for my $refpos (sort {$a<=>$b} keys %{$alt_info_Href->{altpos}}){
            my $mut_info_Href = $alt_info_Href->{altpos}->{$refpos};
            $tmlen[$_] = max( $tmlen[$_], length($mut_info_Href->{$theme[$_]}) ) for (0 .. $#theme);
        }
        # extend BLANK
        $_ += 3 for @tmlen;
        # print Theme
        stout_and_sterr "[INFO]:\t" . join( '', ( map sprintf("%-$tmlen[$_]s",$theme[$_]), (0 .. $#theme) ) ) . "\n";
        # print the mutations
        for my $refpos (sort {$a<=>$b} keys %{$alt_info_Href->{altpos}}){
            my $mut_info_Href = $alt_info_Href->{altpos}->{$refpos};
            stout_and_sterr "[INFO]:\t" . join( '', ( map sprintf("%-$tmlen[$_]s",$mut_info_Href->{$theme[$_]}), (0 .. $#theme) ) ) . "\n";
        }
    }
}

#--- generate finalBam from mkdupBam ---
## including making indexed virus-ref of the last mod-vref as the first step
## this indexed virus-ref will also be applied in function 'summarize_virus_variant'
## and it will be deleted in that function
sub get_variant_finalBam_from_mkdupBam{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    my $variant_ref_folder = File::Spec->catfile( $variant_Href->{workfolder}, 'indexed_variant_ref' );
    `rm -rf $variant_ref_folder`;
    `mkdir -p $variant_ref_folder`;

    # construct indexed reference for NEW virus genome
    my $Indexed_VrefFa = File::Spec->catfile( $variant_ref_folder, 'variant_genome.fa' );
    # specially extend for circular virus
    if( $variant_Href->{Cricle} =~ /yes/i ){
        # reload the virus fa and extend
        FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
            Ref => $variant_Href->{Fa},
            OutFaFile => $Indexed_VrefFa,
            SegName => $variant_Href->{RefSeg},
            LineBase => $V_Href->{line_base},
            CircleExtLen => $variant_Href->{Extend_length}
        });
    }
    else{ # just soft link
        `ln -sf $variant_Href->{Fa} $Indexed_VrefFa`;
    }
    # create the indexed ref
    BWA_index_fasta(FaFile => $Indexed_VrefFa, bwa => $V_Href->{bwa});
    Faidx_Dict_fasta(FaFile => $Indexed_VrefFa, samtools => $V_Href->{samtools});
    # record
    $variant_Href->{FinalIndexVref} = $Indexed_VrefFa;

    # GATK-realn-Samtools-calmd pipeline
    my $temp_bam =
    FuseSV::Virus_Type::SmartAlignment->GATK_relan_Samtools_calmd({
        mkdup_bam   => $variant_Href->{mkdupBam},
        OutPref     => $Indexed_VrefFa, # auto add 'mkdup.bam'
        RefForBwa   => $Indexed_VrefFa,
        SkipCalmd   => 1 # Calmd will introduce loss of MNP
    });
    # store this temp_bam as final-bam
    `mv -f $temp_bam $variant_Href->{finalBam}`;
    `mv -f $temp_bam.bai $variant_Href->{finalBam}.bai`;
    `rm -rf $variant_Href->{mkdupBam} $variant_Href->{mkdupBam}.bai`;

    # infrom
    stout_and_sterr `date`."[INFO]:\tGenerate finalBam for V$variant_Href->{variant_NO} virus variant OK.\n"
                          ."\tBut please note that maybe re-arrange reads later.\n";
}

#--- arrange each reads to its best aligned virus variant ---
## as comparsion has been done during variant reconstruction pipe
## so reads follow the priority from later vNO to former vNO
sub rearrange_reads_to_proper_virus_variant{
    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    # workspace
    my $reads_pool_folder = File::Spec->catfile( $V_Href->{Virus_variant_folder}, 'RearrangeReadsTempSpace' );
    `rm -rf $reads_pool_folder`;
    `mkdir -p $reads_pool_folder`;

    # extract reads from each virus variant's alignment result
    my $readsInfo_list = File::Spec->catfile( $reads_pool_folder, 'AllMappedReadsListForSort.gz' );
    open (RILT, Try_GZ_Write($readsInfo_list)) || die "fail write $readsInfo_list: $!\n";
    for my $virus_variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
        my $this_variant_Bam = $this_variant_Href->{finalBam};
        #   0x800: to avoid supplementary alignment
        # + 0x100: to avoid secondary alignment
        # = 0x900
        open (VBAM,"$V_Href->{samtools} view -F 0x900 $this_variant_Bam |") || die "fail samtools view $this_variant_Bam: $!\n";
        while(<VBAM>){
            my ($pid, $flag) = (split)[0,1];
            my $endNO = ($flag & 0x40) ? 1 : 2;
            my $mBoll = ($flag & 0x4) ? 0 : 1;
            print RILT "$pid\t$endNO\t$mBoll\t$virus_variant_NO\n";
        }
        close VBAM;
    }
    close RILT;

    # sort read-id
    my $sort_Tdir = File::Spec->catfile( $reads_pool_folder, 'AllMappedReadsListSort_Tdir' );
    my $readsInfo_sortlist = File::Spec->catfile( $reads_pool_folder, 'AllMappedReadsList.sorted.gz' );
    `rm -rf $sort_Tdir`;
    `mkdir -p $sort_Tdir`;
    `zcat $readsInfo_list | sort -S 600000 -T $sort_Tdir | gzip -c > $readsInfo_sortlist`;
    `rm -rf $sort_Tdir`;

    my $buffer_maxSize = 1E6;
    # read sorted readslist to make judgement of each pid
    my $readSeletion_list = File::Spec->catfile( $reads_pool_folder, 'Reads.Selection.list.gz' );
    open (my $RSEL_fh, Try_GZ_Write($readSeletion_list)) || die "fail write $readSeletion_list: $!\n";
    my %ReadsInfo;
    my $last_pid = '__NA__';
    my $buffer_count = 0;
    open (SRILT, Try_GZ_Read($readsInfo_sortlist)) || die "fail read $readsInfo_sortlist: $!\n";
    while(<SRILT>){
        my ($pid, $endNO, $mBoll, $vNO) = (split)[0,1,2,3];
        # meet new pid
        if(    $pid ne $last_pid
            && ++$buffer_count > $buffer_maxSize
        ){
            # output pid selection
            &assign_pid_to_variant_NO({
                Out_FH_Sref => \$RSEL_fh,
                ReadsInfo_Href => \%ReadsInfo
            });
            # empty buffer
            %ReadsInfo = ();
            $buffer_count = 1;
        }
        # record readinfo
        $ReadsInfo{$pid}->{$vNO}->{$endNO} = $mBoll;
        # update
        $last_pid = $pid;
    }
    close SRILT;
    # deal with last part
    &assign_pid_to_variant_NO({
        Out_FH_Sref => \$RSEL_fh,
        ReadsInfo_Href => \%ReadsInfo
    });
    close $RSEL_fh;
    # empty buffer
    %ReadsInfo = ();
    $buffer_count = 0;

    my %Reads2Vno;
    my $BamUpdateTime = 0;
    # arrange each pid in virus variant bam
    open (RSEL, Try_GZ_Read($readSeletion_list)) || die "fail read $readSeletion_list: $!\n";
    while(<RSEL>){
        my ($pid, $vNO) = (split);
        $Reads2Vno{$pid} = $vNO;
        if( ++$buffer_count > $buffer_maxSize ){
            &arrange_batch_reads_to_vNO({
                Reads2Vno_Href => \%Reads2Vno,
                BamUpdateTime => ++$BamUpdateTime
            });
            # empty buffer
            %Reads2Vno = ();
            $buffer_count = 0;
        }
    }
    close RSEL;
    # deal with last part
    &arrange_batch_reads_to_vNO({
        Reads2Vno_Href => \%Reads2Vno,
        BamUpdateTime => ++$BamUpdateTime
    });
    # empty buffer
    %Reads2Vno = ();
    $buffer_count = 0;

    # sweep
    unless( $V_Href->{in_debug} ){
        `rm -rf $reads_pool_folder`;
    }

    # inform
    stout_and_sterr `date`."[INFO]:\trearrange related-reads of virus variants ok.\n";
}

#--- deal with %ReadsInfo in function 'rearrange_reads_to_proper_virus_variant'
sub assign_pid_to_variant_NO{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $Out_FH_Sref = $Option_Href->{Out_FH_Sref};
    my $ReadsInfo_Href = $Option_Href->{ReadsInfo_Href};

    for my $pid (sort keys %$ReadsInfo_Href){
        my @vNO = sort {$b<=>$a} keys %{$ReadsInfo_Href->{$pid}};
        # only one vNO
        next if( scalar(@vNO) == 1 );
        # later vNO has best priority
        my $latest_vNO = $vNO[0];
        print { $$Out_FH_Sref } "$pid\t$latest_vNO\n";
        # if one end is unmap, consider former vNO
        my $mBool_1 = $ReadsInfo_Href->{$pid}->{$latest_vNO}->{1};
        my $mBool_2 = $ReadsInfo_Href->{$pid}->{$latest_vNO}->{2};
        unless( $mBool_1 && $mBool_2 ){
            my $check_Mend = ($mBool_1) ? 2 : 1;
            for my $i ( 1 .. $#vNO ){
                my $former_vNO = $vNO[$i];
                if( $ReadsInfo_Href->{$pid}->{$former_vNO}->{$check_Mend} ){
                    print { $$Out_FH_Sref } "$pid\t$former_vNO\n";
                    last; # find OK former vNO, so stop
                }
            }
        }
    }
}

#--- arrange part of reads to its related vNO ---
## because we use buffer to save the memory, so deal part by part
sub arrange_batch_reads_to_vNO{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $Reads2Vno_Href = $Option_Href->{Reads2Vno_Href};
    my $BamUpdateTime = $Option_Href->{BamUpdateTime};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    my @virus_variant_NO = sort {$b<=>$a} keys %$Variant_Info_Href;
    # ignore the last vNO
    for my $virus_variant_NO ( @virus_variant_NO[ 1 .. $#virus_variant_NO ] ){
        my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
        my $this_variant_Bam = $this_variant_Href->{finalBam};
        # prepare output filtered bam
        my $this_variant_fBam = "$this_variant_Bam.ArrangeReads.bam";
        open (FBAM,"| $V_Href->{samtools} view -b -S -o $this_variant_fBam") || die"fail samtools write $this_variant_fBam: $!\n";
        # read original final-bam
        #   0x800: to avoid supplementary alignment
        # + 0x100: to avoid secondary alignment
        # = 0x900
        open (VBAM,"$V_Href->{samtools} view -h -F 0x900 $this_variant_Bam |") || die "fail samtools view $this_variant_Bam: $!\n";
        while(<VBAM>){
            # SAM header
            if(/^\@/){
                print FBAM;
                next;
            }
            my ($pid) = (split)[0];
            if(    !exists($Reads2Vno_Href->{$pid}) # donot consider this pid
                || $Reads2Vno_Href->{$pid} == $virus_variant_NO # this pid belongs to this vNO
            ){
                print FBAM;
            }
        }
        close VBAM;
        close FBAM;
        # change bam name
        ## delete original final-bam and index bai file
        `rm -rf $this_variant_Bam $this_variant_Bam.bai`;
        ## rename new filtered bam
        `mv -f $this_variant_fBam $this_variant_Bam`;
        ## index new filtered bam
        `$V_Href->{samtools} index $this_variant_Bam`;
        # inform
        stout_and_sterr `date`."[INFO]:\trearrange batch of reads of virus variants V$virus_variant_NO (time:$BamUpdateTime) ok.\n";
    }
}

#--- gather each virus variant information ---
## 1) depth, coverage
## 2) mutations on original v-ref
sub gather_virus_variant_info{
    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    for my $virus_variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
        # summarise info of this variant
        &summarize_virus_variant({
            variant_Href => $this_variant_Href
        });
    }

    # inform
    stout_and_sterr `date`."[INFO]:\tgather depth/coverage/mutations info of each virus variant ok.\n";
}

#--- summarize mutations from virus variant determined ---
sub summarize_virus_variant{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # get the cover ratio
    ## $Virus_Info_Href->{ $variant_Href->{OrigRefSeg} }->{length}
    $variant_Href->{CoverRatio} = $variant_Href->{NonN_length} / $variant_Href->{Orig_length};

    # call VCF from finalBam for mutations
    # generate $variant_Href->{VCF}
    &final_VCF_of_virus_variant({
        variant_Href => $variant_Href,
    });

    # stat depth and coverage on this variant alignment bam
    ## here is deal with the bam after 're-arrange reads of each virus variant' to get the REAL information
    &get_depth_coverage_given_bam({
        Bam => $variant_Href->{finalBam},
        check_seg => $variant_Href->{RefSeg},
        seg_orig_length => $variant_Href->{Orig_length},
        seg_nonN_length => $variant_Href->{NonN_length},
        seg_ext_length => $variant_Href->{Extend_length},
        need_details => 1,
        coverage_Href => $variant_Href->{Coverage}, # used when need_details is 1
        averagedepth_Sref => \$variant_Href->{AverageDepth} # used when need_details is 1
    });

    # use blast to map the original ref get major mutations
    &blast_new_to_orig_ref_get_mutations({
        variant_Href => $variant_Href
    });

    # record mutation supportings
    # record other mutations, must be minor mutations
    if( $variant_Href->{ModCompBool} ){

        # mapping table of positions
        my $Mut2Orig_PosMapping_Href = {};
        my $MajorMut_OrigRefAle_Href = {};
        &Construct_Mut2Orig_PosMapping({
            variant_Href => $variant_Href,
            Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
            MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href
        });

        # read mutations from vcf
        my $RefPos_Href = {};
        &filter_mutation_from_VCF({
            RefPos_Href => $RefPos_Href,
            VCF => $variant_Href->{VCF},
            vBwaRef_Fa => $variant_Href->{FinalIndexVref},
            V_refseg => $variant_Href->{OrigRefSeg},
            NoFilter => 1
        });

        # debug
        if( $V_Href->{in_debug} ){
            for my $pos (sort {$a<=>$b} keys %$RefPos_Href){
                my $pos_Href = $RefPos_Href->{$pos};
                my $mut_Href = $pos_Href->{mutation};
                for my $mut_id (sort keys %$mut_Href){
                    my ($sum_sup, $plus_sup, $minus_sup, $refSup_toCmp) = @{ $pos_Href->{mutation}->{$mut_id} };
                    warn "<finalVCFMut> pos $pos,$mut_id: ".
                         "sum_sup:$sum_sup, plus_sup:$plus_sup, minus_sup:$minus_sup, refSup_toCmp:$refSup_toCmp\n";
                }
            }
        }

        # mutation summary
        # get supports of major-mut and their ref-allel
        # get any other credible minor-mut discarded before
        &summarize_all_mutations({
            variant_Href => $variant_Href,
            RefPos_Href => $RefPos_Href,
            Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
            MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href
        });
    }

    # sweep
    my $Indexed_VrefFa_folder = dirname( $variant_Href->{FinalIndexVref} );
    unless( $V_Href->{in_debug} ){
        `rm -rf $Indexed_VrefFa_folder`;
    }

    #++++++++++++++++++++++++++#
    # output the modifications #
    #++++++++++++++++++++++++++#
    &output_modification({
        variant_Href=>$variant_Href
    });
}

#--- call VCF for virus variant ---
sub final_VCF_of_virus_variant{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};

    # call vcf via SamTools mpileup
    my $vcf =
    FuseSV::Virus_Type::VCFwork->vcf_from_bam({
        Bam => $variant_Href->{finalBam},
        OutPref => $variant_Href->{FinalIndexVref}, # auto add '.vcf'
        RefForBwa => $variant_Href->{FinalIndexVref}
    });

    # store VCF to right place
    `mv -f $vcf $variant_Href->{VCF}`;

    # inform
    stout_and_sterr `date`."[INFO]:\tGenerate final VCF for V$variant_Href->{variant_NO} virus variant OK.\n";
}

#--- use blast to get all modifications on the original virus ref ---
# choose to use 'MakeUpMutVref' to get the mutations by BLAST alignment, 2017-05-04
# previously used 'MutVref', has some uncertain situation when has 'lots of N'.
sub blast_new_to_orig_ref_get_mutations{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

    # prepare workspace
    my $work_folder = File::Spec->catfile( $variant_Href->{workfolder}, 'blast_for_mutations' );
    `rm -rf $work_folder`;
    `mkdir -p $work_folder`;
    # files needed
    my $orig_refname = $variant_Href->{OrigRefSeg};
    my $mut_ref = $variant_Href->{mkupFa}; # previous is {Fa}.
    my $mut_ext = $variant_Href->{Extend_length};
    my $vari_mut_Href = $variant_Href->{Mutation};
    my $virus_variant_NO = $variant_Href->{variant_NO};
    my $orig_vRef = $Virus_Info_Href->{$orig_refname}->{vref_fa};
    my $orig_vExt = $Virus_Info_Href->{$orig_refname}->{extend_len};

    # generate the alignment ref for blastall
    my $blastall_ref = File::Spec->catfile( $work_folder, 'blast_ref.fa' );
    ## use the original fa, already extended if possible
    FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
        Ref => $orig_vRef,
        OutFaFile => $blastall_ref,
        SegName => $orig_refname,
        LineBase => $V_Href->{line_base}
        # CircleExtLen => ( -1 * $orig_vExt ) # to trim the extended part if possible
    });
    ## construct blast reference
    FuseSV::Virus_Type::BlastM0Mutation->format_database({
        ref_fa_file => $blastall_ref
    });

    # prepare the query sequence
    my $query_fafile = File::Spec->catfile( $work_folder, 'blast_query.fa' );
    # `rm -rf $query_fafile`;
    # use mutated virus ref splited by 'N' if possible
    ## as makeup-mut-vref has none 'N', so ignore 'split_N' Ture or False.
    FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa({
        Ref => $mut_ref,
        OutFaFile => $query_fafile,
        SegName => 'query_MakeUpMutVref', # 'query_SpN'
        LineBase => $V_Href->{line_base},
        CircleExtLen => $mut_ext, # for good-alignment of circular virus
        split_N => 1
    });
    # # use whole mutated virus ref
    # `cat $mut_ref >> $query_fafile`;
    # blastall query
    my $blast_m0 = File::Spec->catfile( $work_folder, 'blast_query.m0' );
    FuseSV::Virus_Type::BlastM0Mutation->blastall_query_m0({
        ref_fa_file => $blastall_ref,
        query_fa_file => $query_fafile,
        blast_m0 => $blast_m0
    });

    # upload mutations(modifications)
    ## in this function, merge the extended part with the original part.
    my $RefPos_Href =
    FuseSV::Virus_Type::BlastM0Mutation->read_blast_m0_for_mutation({
        blast_m0 => $blast_m0,
        query_seq_fa => $query_fafile,
        use_all_map => 1, # query_seq may break down to several parts
        orig_fa_ref => $blastall_ref,
        orig_ext => $orig_vExt,
        ext_merge_Consider_Region => $orig_vExt,
        mutation_Href => $vari_mut_Href,
        min_depth => 1,
        both_strd => 0,
        min_mut_reads => 1,
        need_RefPos_Hash => 1,
        without_supports => 1
    });

    # inform
    stout_and_sterr `date`."[INFO]:\tblast virus mod-fa (V$virus_variant_NO) against initial database ref-fa ($orig_refname) OK.\n";

    # sweep
    unless( $V_Href->{in_debug} ){
        `rm -rf $work_folder`;
    }

    # check whether all mutations are listed
    my $mutation_complete_sign =
    &check_whether_mutations_completed({
        mut_Ref => $mut_ref,
        mutation_Href => $vari_mut_Href,
        RefPos_Href => $RefPos_Href
    });
    # warn
    if( $mutation_complete_sign == 0 ){
        stout_and_sterr `date`."<WARN>:\tMutations obtained via blast-map are not completed to construct mutated v-ref (V$virus_variant_NO).\n";
    }
    else{
        stout_and_sterr `date`."[INFO]:\tMutations obtained via blast-map are completed to construct mutated v-ref (V$virus_variant_NO).\n";
    }
    # record the complete sign
    $variant_Href->{ModCompBool} = $mutation_complete_sign;

    # merge the SNV to MNP
    # &merge_continuousSNV_to_MNP({
    #   mutation_Href => $vari_mut_Href
    # });
}

#--- check whether mutations are completed ---
sub check_whether_mutations_completed{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $mut_Href = $Option_Href->{mutation_Href};
    my $RefPos_Href = $Option_Href->{RefPos_Href};
    my $mut_Ref = $Option_Href->{mut_Ref};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # concatenate ref and mutated bases
    my $cat_seq = '';
    my $pos_filter = 0;
    for my $pos (sort {$a<=>$b} keys %$RefPos_Href){
        next if( $pos < $pos_filter);
        my $ref_allel = $RefPos_Href->{$pos}->{ref_allel};
        if( exists($mut_Href->{$pos}) ){
            my $mut_type = $mut_Href->{$pos}->{mut_type};
            my $mut_seq = $mut_Href->{$pos}->{mut_seq};
            if( $mut_type eq 'snp' ){
                $cat_seq .= $mut_seq;
            }
            elsif( $mut_type eq 'ins' ){
                $cat_seq .= $ref_allel.$mut_seq;
            }
            elsif( $mut_type eq 'del' ){
                $pos_filter = $pos + length($mut_seq);
            }
            # debug
            if( $V_Href->{in_debug} ){
                warn "BlastCheckModCompleteMod:\t$pos\t$mut_type\t$mut_seq\n";
            }
        }
        else{
            $cat_seq .= $ref_allel;
        }
    }
    # debug
    if( $V_Href->{in_debug} ){
        warn `date`."BlastCheckModComplete_ModFormed_Vseq:\n$cat_seq\n";
    }
    # mut ref fasta
    my $mut_ref_seq = `cat $mut_Ref | grep -v '^>'`;
    $mut_ref_seq =~ s/\s//g;
    # compare
    my $mutation_complete_sign = 1;
    if( length($cat_seq) == length($mut_ref_seq) ){
        my @cat_seq_base = split //, uc($cat_seq);
        my @mut_ref_base = split //, uc($mut_ref_seq);
        for my $i (0 .. $#cat_seq_base){
            if(    $cat_seq_base[$i] ne $mut_ref_base[$i]
                && $cat_seq_base[$i] ne 'N'
                && $mut_ref_base[$i] ne 'N'
              ){
                $mutation_complete_sign = 0;
                # inform
                stout_and_sterr `date`."ModsNotComplete_Reason:\tPos_index [$i] is different from mut_vref: $cat_seq_base[$i] ne $mut_ref_base[$i].\n";
                last;
            }
        }
    }
    else{
        $mutation_complete_sign = 0;
        # inform
        stout_and_sterr `date`."ModsNotComplete_Reason:\tLength is different.\n";
    }

    return $mutation_complete_sign;
}

#--- merge continuous SNV to MNP
sub merge_continuousSNV_to_MNP{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $mut_Href = $Option_Href->{mutation_Href};

    my @SNV_pos = grep $mut_Href->{$_}->[0] =~ /snp/i, (sort {$a<=>$b} keys %$mut_Href);
    for (my $i=0; $i<$#SNV_pos; $i++){
        my $now_index = $i;
        while( $SNV_pos[$i+1] == $SNV_pos[$i]+1 ){ # yes, next pos is still a SNP
            $i++;
            last if($i >= $#SNV_pos);
        }
        if( $i != $now_index ){ # index has updated to following ones
            $mut_Href->{$SNV_pos[$now_index]}->[0] = 'mnp';
            for my $flow_i ($now_index+1 .. $i){
                $mut_Href->{$SNV_pos[$now_index]}->[1] .= $mut_Href->{$SNV_pos[$flow_i]}->[1];
                delete $mut_Href->{$SNV_pos[$flow_i]};
            }
        }
    }
}

#--- position mapping table from mut-ref to orig-ref
# the mapping table is position <-> position
sub Construct_Mut2Orig_PosMapping{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};
    my $Mut2Orig_PosMapping_Href = $Option_Href->{Mut2Orig_PosMapping_Href};
    my $MajorMut_OrigRefAle_Href = $Option_Href->{MajorMut_OrigRefAle_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

    my $orig_refname = $variant_Href->{OrigRefSeg};
    my $vmod_refname = $variant_Href->{RefSeg};
    my $vari_mut_Href = $variant_Href->{Mutation};
    my $orig_vLen = $Virus_Info_Href->{$orig_refname}->{length};
    my $orig_vExt = $Virus_Info_Href->{$orig_refname}->{extend_len};

    # read throughout the major mutations
    my $last_MutPos  = 0;
    my $last_OrigPos = 0;
    for my $OrigPos (sort {$a<=>$b} keys %$vari_mut_Href){
        # deal with positions in last interval
        my $interval_size = $OrigPos - $last_OrigPos - 1;
        for my $i (1 .. $interval_size){
            $Mut2Orig_PosMapping_Href->{ $last_MutPos+$i } = $last_OrigPos+$i;
        }
        # deal with this major mutation
        my $ref_allel = $vari_mut_Href->{$OrigPos}->{ref_allel};
        my $mut_type = $vari_mut_Href->{$OrigPos}->{mut_type};
        my $mut_seq = $vari_mut_Href->{$OrigPos}->{mut_seq};
        # different scenarios
        ## majorMut-snp is orig_Vref-snp
        if( $mut_type eq 'snp' ){
            my $Pos_shift = $OrigPos - $last_OrigPos;
            my $MutPos = $last_MutPos + $Pos_shift;
            $Mut2Orig_PosMapping_Href->{$MutPos} = $OrigPos;
            $MajorMut_OrigRefAle_Href->{$MutPos} = { mut_type=>'snp', mut_seq=>$ref_allel, mutHref=>$vari_mut_Href->{$OrigPos} };
            $last_MutPos = $MutPos;
            $last_OrigPos = $OrigPos;
        }
        ## majorMut-del is orig_Vref-ins
        elsif( $mut_type eq 'del' ){
            $last_OrigPos = $OrigPos + length($mut_seq) - 1;
            $last_MutPos += $interval_size;
            # assign 'insertion' to the position before it
            $MajorMut_OrigRefAle_Href->{$last_MutPos} = { mut_type=>'ins', mut_seq=>$mut_seq, mutHref=>$vari_mut_Href->{$OrigPos} };
        }
        ## majorMut-ins is orig_Vref-del
        elsif( $mut_type eq 'ins' ){
            $last_OrigPos += $interval_size + 1;
            $last_MutPos  += $interval_size + 1;
            $Mut2Orig_PosMapping_Href->{$last_MutPos} = $last_OrigPos;
            # assign 'deletion' to the first position of deleted part, so add 1
            $MajorMut_OrigRefAle_Href->{$last_MutPos+1} = { mut_type=>'del', mut_seq=>$mut_seq, mutHref=>$vari_mut_Href->{$OrigPos} };
            $last_MutPos  += length($mut_seq);
        }
    }
    # remained region
    for my $i (1 .. $orig_vLen-$last_OrigPos){
        $Mut2Orig_PosMapping_Href->{ $last_MutPos+$i } = $last_OrigPos+$i;
    }

    # output PosMapping Table
    my @theme = qw/ VirusVariant
                    ModRegion_start
                    ModRegion_end
                    OriginalVirus
                    OrigRegion_start
                    OrigRegion_end
                  /;
    open (PML, Try_GZ_Write($variant_Href->{PosMapList})) || die "fail to write PosMapping List: $!\n";
    # theme line
    print PML '#' . join("\t", @theme) . "\n";
    # region records
    ($last_MutPos, $last_OrigPos) = (0, 0); # reuse
    my ($st_MutPos, $st_OrigPos) = (0, 0);
    for my $MutPos (sort {$a<=>$b} keys %$Mut2Orig_PosMapping_Href){
        my $OrigPos = $Mut2Orig_PosMapping_Href->{$MutPos};
        # first get
        if( $st_MutPos == 0 ){
            $st_MutPos = $MutPos;
            $st_OrigPos = $OrigPos;
        }
        # region segmentation, pos not continuous
        if( abs($MutPos-$st_MutPos) != abs($OrigPos-$st_OrigPos) ){
            print PML join("\t", $vmod_refname, $st_MutPos, $last_MutPos, $orig_refname, $st_OrigPos, $last_OrigPos) . "\n";
            # update
            $st_MutPos = $MutPos;
            $st_OrigPos = $OrigPos;
        }
        # update
        $last_MutPos = $MutPos;
        $last_OrigPos = $OrigPos;
    }
    # last part
    print PML join("\t", $vmod_refname, $st_MutPos, $last_MutPos, $orig_refname, $st_OrigPos, $last_OrigPos) . "\n";
    close PML;

    if( $V_Href->{in_debug} ){
        for my $MutPos (sort {$a<=>$b} keys %$Mut2Orig_PosMapping_Href){
            warn "<PosMapping>:\tM $MutPos : O $Mut2Orig_PosMapping_Href->{$MutPos}\n";
        }
        for my $MutPos (sort {$a<=>$b} keys %$MajorMut_OrigRefAle_Href){
            my $majorMut_Href = $MajorMut_OrigRefAle_Href->{$MutPos};
            warn "<MajMut_RefA>:\tPos $MutPos\t$majorMut_Href->{mut_type}  $majorMut_Href->{mut_seq}\n";
        }
    }

    # inform
    stout_and_sterr `date`."[INFO]:\tObtain MutRef-to-OrigRef Pos-mapping table for V$variant_Href->{variant_NO} virus variant OK.\n"
                                 ."\tOutput MutRef-to-OrigRef Pos-mapping list OK.\n";
}

#--- distinguish and record Major/Minor mutations for virus variant ---
sub summarize_all_mutations{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};
    my $RefPos_Href = $Option_Href->{RefPos_Href};
    my $Mut2Orig_PosMapping_Href = $Option_Href->{Mut2Orig_PosMapping_Href};
    my $MajorMut_OrigRefAle_Href = $Option_Href->{MajorMut_OrigRefAle_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};
    my $min_depth = $V_Href->{min_depth};
    my $both_strd = $V_Href->{both_strd};
    my $min_mut_reads = $V_Href->{min_mut_reads};

    # load original virus variant reference sequence for followed possible substr for origMut-del
    # use it in very essential scenario
    my $OrigVrefFa = $Virus_Info_Href->{ $variant_Href->{OrigRefSeg} }->{vref_fa};
    chomp(my $OrigVrefSeqStr = `cat $OrigVrefFa | grep -v '^>' | perl -ne 'chomp;print;'`);
    # debug
    if( $V_Href->{in_debug} ){
        warn "OrigVrefSeqStr:\n$OrigVrefSeqStr\n";
    }

    my $VariantMut_Href = $variant_Href->{Mutation};
    # filter mutations
    for my $pos (sort {$a<=>$b} keys %$RefPos_Href){
        my $pos_Href = $RefPos_Href->{$pos};
        # count the ref-allel depth to major-mut in virus variant's mod-Vref
        if( exists($MajorMut_OrigRefAle_Href->{$pos}) ){
            ## this is major-mutation's supports
            my $majorMut_Href = $MajorMut_OrigRefAle_Href->{$pos}->{mutHref};
            ## NOTE: modVrefpos is tranformed for orig_Vref-allele action
            ## 1) majorMut-ins is recorded as orig_Vref-del
            ##    so, majorMut-ins's modVrefpos is recorded as the first base of majorMut-ins_seq for orig_Vref-del
            ## 2) majorMut-del is recorded as orig_Vref-ins
            ##    so, majorMut-del's modVrefpos is recorded as the base before majorMut-del_seq for orig_Vref-ins
            ## 3) majorMut-snp is never minded.
            $majorMut_Href->{modVrefpos} = $pos;
            $majorMut_Href->{mut_sum_sup}   = $pos_Href->{refallel_depth};
            $majorMut_Href->{mut_plus_sup}  = $pos_Href->{refallel_D_Fw};
            $majorMut_Href->{mut_minus_sup} = $pos_Href->{refallel_D_Rv};
        }
        # record mutations
        if(   $pos_Href->{depth} < $min_depth   # depth filtration
           || $pos_Href->{ref_allel} =~ /^N$/i  # do not allow mutation at 'N' locus
          ){
            delete $pos_Href->{mutation}; # delete it anyway
        }
        # has mutation(s), so check whether it is good enough
        elsif( exists $pos_Href->{mutation} ){
            # deal with each mutation
            for my $mut_id (sort keys %{$pos_Href->{mutation}}){
                my ($mut_type, $mut_seq) = split /,/, $mut_id;
                my ($sum_sup, $plus_sup, $minus_sup, $refSup_toCmp) = @{ $pos_Href->{mutation}->{$mut_id} };
                # filter: read-support
                if(    ( $sum_sup < $min_mut_reads ) # at least so many supported reads
                    || ( $both_strd && ($plus_sup < $both_strd || $minus_sup < $both_strd) ) # if required, both strand supported
                  ){
                    next;
                }

                #---------------------------------------------------#
                # 1) whether is related to major-mut? Y or N
                # 2) this mut-pos has mapping original pos? Y or N
                #---------------------------------------------------#

                if(    exists($MajorMut_OrigRefAle_Href->{$pos}) # 1) Y
                    && exists($Mut2Orig_PosMapping_Href->{$pos}) # 2) Y, majorMut type must be 'snp' or 'del'
                ){
                    # mapping orig-pos
                    my $orig_pos = $Mut2Orig_PosMapping_Href->{$pos};
                    # this is orig-Vref allele's action
                    my $majorMut_OrigRefAle_mutType = $MajorMut_OrigRefAle_Href->{$pos}->{mut_type};
                    # this is orig-Vref allele's action seq
                    my $majorMut_OrigRefAle_mutSeq  = $MajorMut_OrigRefAle_Href->{$pos}->{mut_seq};
                    # this is ref of majotMut info hash
                    my $majorMut_Href = $MajorMut_OrigRefAle_Href->{$pos}->{mutHref};
                    # it is the original allele appearance at the major-mut site
                    # so, this is original-ref-allel's supports
                    if(       $mut_type  eq    $majorMut_OrigRefAle_mutType
                        && uc($mut_seq)  eq uc($majorMut_OrigRefAle_mutSeq)
                      ){
                        $majorMut_Href->{ref_sum_sup}   = $sum_sup;
                        $majorMut_Href->{ref_plus_sup}  = $plus_sup;
                        $majorMut_Href->{ref_minus_sup} = $minus_sup;
                    }
                    # although at the major-mut site, but is other minor mutation
                    # note that this mutation info is against mod-Vref
                    else{
                        # reocrd the original modification
                        my $orig_mod;
                        # the majorMut is 'snp'
                        if( $majorMut_Href->{mut_type} eq 'snp' ){
                            ## for 'snp' or 'ins', just record.
                            if(    $mut_type eq 'snp'
                                || $mut_type eq 'ins'
                              ){
                                # no contradiction with majotMut
                                $orig_mod = "$orig_pos,$mut_type," . uc( $mut_seq );
                            }
                            ## for 'del', try to check overlap with other followed majorMut
                            elsif( $mut_type eq 'del' ){
                                my @del_seq = split //, $mut_seq;
                                # deleted part should start wtih the orig-base
                                my $orig_del_seq = uc($majorMut_OrigRefAle_mutSeq);
                                # already base on this pos, discard this base (may-not-orig-allele)
                                shift @del_seq;
                                # check overlap with other followed majorMut
                                my $judge_pos = $pos + 1;
                                my $Poss_OrigMod_ToPhase = '';
                                &extended_del_to_followed_reg({
                                    del_seq_Aref => \@del_seq,
                                    MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href,
                                    Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
                                    judge_pos => $judge_pos,
                                    orig_del_seq_Sref => \$orig_del_seq,
                                    orig_vref_str_Sref => \$OrigVrefSeqStr,
                                    Poss_OrigMod_ToPhase_Sref => \$Poss_OrigMod_ToPhase
                                });
                                # orig modification
                                # $Poss_OrigMod_ToPhase may be empty
                                $orig_mod = "$orig_pos,del,$orig_del_seq" . $Poss_OrigMod_ToPhase;
                            }
                        }
                        # the majorMut is 'del'
                        elsif( $majorMut_Href->{mut_type} eq 'del' ){
                            my $majorMut_del_seq = $majorMut_Href->{mut_seq};
                            my $majorMut_del_len = length($majorMut_del_seq);
                            ## for 'snp', just record
                            if( $mut_type eq 'snp' ){
                                # no contradiction with majotMut
                                $orig_mod = "$orig_pos,snp," . uc( $mut_seq );
                            }
                            ## for 'ins'
                            ## 1) compare with del_seq;
                            ## 2) and try to see the relative location of ins_seq and del_seq
                            elsif( $mut_type eq 'ins' ){
                                my $ins_len = length($mut_seq);
                                # ins-seq is shorter than del-seq
                                # AND, ins-seq is from del-seq
                                if(    $ins_len < $majorMut_del_len
                                    && $majorMut_del_seq =~ /$mut_seq/i
                                ){
                                    # 'del_seq harbors ins_seq'
                                    # ins_seq makeup the del-tail
                                    if( $majorMut_del_seq =~ /$mut_seq$/i ){
                                        my $orig_del_seq = uc( substr($majorMut_del_seq, 0, $majorMut_del_len-$ins_len) );
                                        # deleted_1st_pos is next base of this recorded orig_pos
                                        my $orig_del_st_pos = $orig_pos + 1;
                                        $orig_mod = "$orig_del_st_pos,del,$orig_del_seq";
                                    }
                                    # 'del_seq harbors ins_seq'
                                    # ins_seq makeup the del-head
                                    elsif( $majorMut_del_seq =~ /^$mut_seq/i ){
                                        my $orig_del_seq = uc( substr($majorMut_del_seq, $ins_len) );
                                        my $orig_del_st_pos = $orig_pos + $ins_len + 1;
                                        $orig_mod = "$orig_del_st_pos,del,$orig_del_seq";
                                    }
                                    # 'del_seq harbors ins_seq'
                                    # ins_seq locates in the del_seq
                                    # this devides the majorMut-del into two small phased deletions
                                    elsif( $majorMut_del_seq =~ /[ACGTN]+$mut_seq[ACGTN]+/i ){
                                        my $ins_idx = index( uc($majorMut_del_seq), uc($mut_seq) );
                                        # first small deletion, 5-primer
                                        my $orig_del_seq = uc( substr($majorMut_del_seq, 0, $ins_idx) );
                                        # deleted_1st_pos is next base of this recorded orig_pos
                                        my $orig_del_st_pos = $orig_pos + 1;
                                        $orig_mod = "$orig_del_st_pos,del,$orig_del_seq";
                                        # second small deletion, 3-primer
                                        my $next_del_seq = uc( substr($majorMut_del_seq, $ins_idx+$ins_len) );
                                        my $next_del_st_pos = $orig_del_st_pos + $ins_idx + $ins_len;
                                        $orig_mod .= '|' . "$next_del_st_pos,del,$next_del_seq";
                                    }
                                }
                                # ins-seq is longer than del-seq
                                # AND, ins-seq makes up del-seq, and extends
                                elsif(    $ins_len > $majorMut_del_len
                                       && $mut_seq =~ /$majorMut_del_seq/i
                                ){
                                    # 'ins_seq harbors del_seq'
                                    # del_seq is ins-tail
                                    if( $mut_seq =~ /$majorMut_del_seq$/i ){
                                        my $orig_ins_seq = uc( substr($mut_seq, 0, $ins_len-$majorMut_del_len) );
                                        $orig_mod = "$orig_pos,ins,$orig_ins_seq";
                                    }
                                    # 'ins_seq harbors del_seq'
                                    # del_seq is ins-head
                                    elsif( $mut_seq =~ /^$majorMut_del_seq/i ){
                                        my $orig_ins_seq = uc( substr($mut_seq, $majorMut_del_len) );
                                        my $orig_ins_st_pos = $orig_pos + $majorMut_del_len;
                                        $orig_mod = "$orig_ins_st_pos,ins,$orig_ins_seq";
                                    }
                                    # 'ins_seq harbors del_seq'
                                    # del_seq locates in ins_seq
                                    # this devides the ins into two small phased insertions
                                    elsif( $mut_seq =~ /[ACGTN]+$majorMut_del_seq[ACGTN]+/i ){
                                        my $del_idx = index( uc($mut_seq), uc($majorMut_del_seq) );
                                        # first small insertion, 5-primer
                                        my $orig_ins_seq = uc( substr($mut_seq, 0, $del_idx) );
                                        $orig_mod = "$orig_pos,ins,$orig_ins_seq";
                                        # second small insertion, 3-primer
                                        my $next_ins_seq = uc( substr($mut_seq, $del_idx+$majorMut_del_len) );
                                        my $next_ins_st_pos = $orig_pos + $majorMut_del_len;
                                        $orig_mod .= '|' . "$next_ins_st_pos,ins,$next_ins_seq";
                                    }
                                }
                                # 'del_seq doesn't harbor ins_seq'
                                # never mind the sizes of this minorMut-ins and majorMut-del
                                elsif(    $majorMut_del_seq !~ /$mut_seq/i
                                       && $mut_seq !~ /$majorMut_del_seq/i
                                ){
                                    # first record this minorMut-ins
                                    $orig_mod = "$orig_pos,ins," . uc( $mut_seq );
                                    # then, must record the majorMut-del to phase with this minorMut ins
                                    # deleted_1st_pos is next base of this recorded orig_pos
                                    my $majorDel_orig_st_pos = $orig_pos + 1;
                                    $orig_mod .= '|' . "$majorDel_orig_st_pos,del," . uc( $majorMut_del_seq );
                                }
                            }
                            ## for 'del', try to see whether overlap with other followed majorMut
                            elsif( $mut_type eq 'del' ){
                                my @del_seq = split //, $mut_seq;
                                # already base on this pos
                                # use orig-allele from orig-vref-seq
                                my $orig_del_seq = uc( substr($OrigVrefSeqStr, $orig_pos-1, 1) );
                                # already base on this pos, discard this base (may-not-orig-allele)
                                shift @del_seq;
                                # concatenate majorMut deletion
                                # $majorMut_del_seq is already orig-allele(s)
                                $orig_del_seq .= uc( $majorMut_del_seq );
                                # check overlap with other followed majorMut
                                my $judge_pos = $pos + 1;
                                my $Poss_OrigMod_ToPhase = '';
                                &extended_del_to_followed_reg({
                                    del_seq_Aref => \@del_seq,
                                    MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href,
                                    Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
                                    judge_pos => $judge_pos,
                                    orig_del_seq_Sref => \$orig_del_seq,
                                    orig_vref_str_Sref => \$OrigVrefSeqStr,
                                    Poss_OrigMod_ToPhase_Sref => \$Poss_OrigMod_ToPhase
                                });
                                # orig modification
                                # $Poss_OrigMod_ToPhase may be empty
                                $orig_mod = "$orig_pos,del,$orig_del_seq" . $Poss_OrigMod_ToPhase;
                            }
                        }
                        # record this other deletion and its related original modification
                        my $mut_id_ToRecord = "$pos,$mut_id";
                        $majorMut_Href->{other_mutation}->{$mut_id_ToRecord}
                                                    = {
                                                       mut_pos => $pos,
                                                       mut_type => $mut_type,
                                                       mut_seq => uc($mut_seq),
                                                       mut_id => $mut_id,
                                                       mut_sum_sup => $sum_sup,
                                                       mut_plus_sup => $plus_sup,
                                                       mut_minus_sup => $minus_sup,
                                                       orig_mod => $orig_mod
                                                      };
                    }
                }
                elsif(      exists($MajorMut_OrigRefAle_Href->{$pos}) # 1) Y
                       && ! exists($Mut2Orig_PosMapping_Href->{$pos}) # 2) N, majorMut must be 'ins'
                ){
                    # nearby 5p base be must has mapping orig_pos
                    my $last5p_pos = $pos - 1;
                    unless( exists($Mut2Orig_PosMapping_Href->{$last5p_pos}) ){
                        cluck_and_exit `date`."<ERROR>\tmajorMut ins unmappable pos ($pos) 5p nearest major-mut ($last5p_pos) "
                                            ."still be unmappable in virus variant V$variant_Href->{variant_NO}.\n";
                    }
                    # mapping orig-pos
                    my $orig_pos = $Mut2Orig_PosMapping_Href->{$last5p_pos};
                    # this is orig-Vref allele's action
                    my $majorMut_OrigRefAle_mutType = $MajorMut_OrigRefAle_Href->{$pos}->{mut_type};
                    # this is orig-Vref allele's action seq
                    my $majorMut_OrigRefAle_mutSeq  = $MajorMut_OrigRefAle_Href->{$pos}->{mut_seq};
                    # this is ref of majotMut info hash
                    my $majorMut_Href = $MajorMut_OrigRefAle_Href->{$pos}->{mutHref};
                    # it is the original allele appearance at the major-mut site
                    # so, this is original-ref-allel's supports
                    if(       $mut_type  eq    $majorMut_OrigRefAle_mutType
                        && uc($mut_seq)  eq uc($majorMut_OrigRefAle_mutSeq)
                      ){
                        $majorMut_Href->{ref_sum_sup}   = $sum_sup;
                        $majorMut_Href->{ref_plus_sup}  = $plus_sup;
                        $majorMut_Href->{ref_minus_sup} = $minus_sup;
                    }
                    # although at the major-mut site, but is other minor mutation
                    # note that this mutation info is against mod-Vref
                    else{
                        # reocrd the original modification
                        my $orig_mod;
                        # the majorMut is 'ins'
                        if( $majorMut_Href->{mut_type} eq 'ins' ){
                            my $majorMut_ins_seq = $majorMut_Href->{mut_seq};
                            my $majorMut_ins_len = length($majorMut_ins_seq);
                            ## for 'snp', replace first base
                            if( $mut_type eq 'snp' ){
                                my $orig_ins_seq = uc(  $mut_seq
                                                      . substr($majorMut_ins_seq, 1) );
                                $orig_mod = "$orig_pos,ins,$orig_ins_seq";
                            }
                            ## for 'ins', extend iinsertion
                            elsif( $mut_type eq 'ins' ){
                                my $orig_ins_seq = uc(   substr($majorMut_ins_seq, 0, 1)     # first base of majorMut-ins 
                                                       . $mut_seq                            # new ins_seq
                                                       . substr($majorMut_ins_seq, 1)     ); # remains of majorMut-ins
                                $orig_mod = "$orig_pos,ins,$orig_ins_seq";
                            }
                            ## for 'del'
                            ## 1) compare with ins_seq;
                            ## 2) and try to see whether overlap with other followed majorMut
                            elsif( $mut_type eq 'del' ){
                                my $ins_len_diff = $majorMut_ins_len - length($mut_seq);
                                # still keep 3-prime part of majorMut-ins
                                if( $ins_len_diff > 0 ){
                                    my $orig_ins_seq = substr($majorMut_ins_seq, length($mut_seq));
                                    $orig_mod = "$orig_pos,ins,$orig_ins_seq";
                                }
                                # delete more than majorMut-ins
                                elsif( $ins_len_diff < 0 ){
                                    # try to check overlap with other majorMut
                                    my @del_seq = split //, substr($mut_seq, $majorMut_ins_len);
                                    # already delete the base after majorMut-ins
                                    my $del_1st_pos = $pos + $majorMut_ins_len;
                                    my $del_1st_orig_pos = $Mut2Orig_PosMapping_Href->{$del_1st_pos};
                                    my $orig_del_seq = uc( substr($OrigVrefSeqStr, $del_1st_orig_pos-1, 1) );
                                    # already base on this pos, discard this base (may-not-orig-allele)
                                    shift @del_seq;
                                    # check overlap with other followed majorMut
                                    my $judge_pos = $del_1st_pos + 1;
                                    my $Poss_OrigMod_ToPhase = '';
                                    &extended_del_to_followed_reg({
                                        del_seq_Aref => \@del_seq,
                                        MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href,
                                        Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
                                        judge_pos => $judge_pos,
                                        orig_del_seq_Sref => \$orig_del_seq,
                                        orig_vref_str_Sref => \$OrigVrefSeqStr,
                                        Poss_OrigMod_ToPhase_Sref => \$Poss_OrigMod_ToPhase
                                    });
                                    # orig modification
                                    # $Poss_OrigMod_ToPhase may be empty
                                    $orig_mod = "$del_1st_orig_pos,del,$orig_del_seq" . $Poss_OrigMod_ToPhase;
                                }
                            }
                        }
                        else{
                            cluck_and_exit `date`."<ERROR>\tmajorMut ins unmappable pos ($pos) 5p nearest major-mut ($last5p_pos) "
                                                ."still be unmappable in virus variant V$variant_Href->{variant_NO}.\n";
                        }
                        # record this other deletion and its related original modification
                        my $mut_id_ToRecord = "$pos,$mut_id";
                        $majorMut_Href->{other_mutation}->{$mut_id_ToRecord}
                                                    = {
                                                       mut_pos => $pos,
                                                       mut_type => $mut_type,
                                                       mut_seq => uc($mut_seq),
                                                       mut_id => $mut_id,
                                                       mut_sum_sup => $sum_sup,
                                                       mut_plus_sup => $plus_sup,
                                                       mut_minus_sup => $minus_sup,
                                                       orig_mod => $orig_mod
                                                      };
                    }
                }
                # it is not at major-mut site, so it must be minor-mut (not included in mod-Vref), never known before
                # here, create this KEY of this pos in virus variant mutation hash to recode this minor-mut
                # it is in the pos-mapping, so must be the pos from orig-Vref
                # this pos may be the one very before majorMut-ins, let them compare with each other to get the best *at last*.
                elsif(    ! exists($MajorMut_OrigRefAle_Href->{$pos}) # 1) N
                       &&   exists($Mut2Orig_PosMapping_Href->{$pos}) # 2) Y
                ){
                    my $orig_pos = $Mut2Orig_PosMapping_Href->{$pos};
                    # reocrd the original modification
                    my $orig_mod;
                    ## for 'snp' or 'ins', just record
                    if( $mut_type eq 'snp' ){
                        $orig_mod = "$orig_pos,$mut_type," . uc( $mut_seq );
                    }
                    elsif( $mut_type eq 'ins' ){
                        # check next pos, is it recorded as a majorMut-Ins, only Ins
                        # if it is, need to connect this the majorMut seq
                        my $next_pos = $pos + 1;
                        if(    exists($MajorMut_OrigRefAle_Href->{$next_pos})
                            && $MajorMut_OrigRefAle_Href->{$next_pos}->{mutHref}->{mut_type} eq 'ins'
                        ){
                            my $majorMut_InsSeq = $MajorMut_OrigRefAle_Href->{$next_pos}->{mut_seq};
                            $orig_mod = "$orig_pos,$mut_type," . uc( $mut_seq.$majorMut_InsSeq );
                        }
                        else{
                            $orig_mod = "$orig_pos,$mut_type," . uc( $mut_seq );
                        }
                    }
                    ## for 'del', try to see whether overlap with other followed majorMut
                    elsif( $mut_type eq 'del' ){
                        # try to check overlap with other majorMut
                        my @del_seq = split //, $mut_seq;
                        # already base on this pos
                        ## as this is not majorMut site
                        ## so, at this base, mod-Vref is identical with orig-Vref
                        my $orig_del_seq = uc( shift @del_seq );
                        # check overlap with other followed majorMut
                        my $judge_pos = $pos + 1;
                        my $Poss_OrigMod_ToPhase = '';
                        &extended_del_to_followed_reg({
                            del_seq_Aref => \@del_seq,
                            MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href,
                            Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
                            judge_pos => $judge_pos,
                            orig_del_seq_Sref => \$orig_del_seq,
                            orig_vref_str_Sref => \$OrigVrefSeqStr,
                            Poss_OrigMod_ToPhase_Sref => \$Poss_OrigMod_ToPhase
                        });
                        # orig modification
                        # $Poss_OrigMod_ToPhase may be empty
                        $orig_mod = "$orig_pos,del,$orig_del_seq" . $Poss_OrigMod_ToPhase;
                    }
                    # initialize this minor mutation that not included in mod-Vref
                    if( ! exists($VariantMut_Href->{$orig_pos}) ){
                        $VariantMut_Href->{$orig_pos} = {
                                                           modVrefpos => $pos,
                                                           ref_allel => $pos_Href->{ref_allel},
                                                           ref_sum_sup => $pos_Href->{refallel_depth},
                                                           ref_plus_sup => $pos_Href->{refallel_D_Fw},
                                                           ref_minus_sup => $pos_Href->{refallel_D_Rv},
                                                           major_mut_bool => 0,
                                                           other_mutation => {}
                                                        };
                    }
                    # alway record as the other mutation (minorMut)
                    my $mut_id_ToRecord = "$pos,$mut_id";
                    $VariantMut_Href->{$orig_pos}->{other_mutation}->{$mut_id_ToRecord}
                                                             = {
                                                                   mut_pos => $pos,
                                                                   mut_type => $mut_type,
                                                                   mut_seq => uc($mut_seq),
                                                                   mut_id => $mut_id,
                                                                   mut_sum_sup => $sum_sup,
                                                                   mut_plus_sup => $plus_sup,
                                                                   mut_minus_sup => $minus_sup,
                                                                   orig_mod => $orig_mod
                                                               };
                }
                # it is not in the pos-mapping, so must be out of the pos from orig-Vref
                # so, it must locate in the insertion introduced by major-mut, mut_type must be 'ins'
                # if snp/ins, just update the inserted string, and record as 'other_mutation' at major-mut site.
                # if del, must check whether it reached out the inserted string?
                #    if not, deal it similar to 'snp/ins'.
                #    if yes, besides action similar to 'snp/ins', adds anther new minor-mut (deletion)
                elsif(    ! exists($MajorMut_OrigRefAle_Href->{$pos}) # 1) N
                       && ! exists($Mut2Orig_PosMapping_Href->{$pos}) # 2) N
                ){
                    # find the adjacent mappable pos to 5' end
                    my $next5p_mappable_Pos = $pos - 1;
                    while( !exists($Mut2Orig_PosMapping_Href->{$next5p_mappable_Pos}) ) {
                        $next5p_mappable_Pos--;
                    }
                    my $major_ins_orig_pos = $Mut2Orig_PosMapping_Href->{$next5p_mappable_Pos};
                    # next pos must be major-mut insertion
                    my $major_ins_pos = $next5p_mappable_Pos + 1;
                    # check major-mut existence
                    unless( exists($MajorMut_OrigRefAle_Href->{$major_ins_pos}) ){
                        cluck_and_exit `date`."<ERROR>\tunmappable pos ($pos) 5p nearest major-mut ($major_ins_pos) is not major-mut site in virus variant V$variant_Href->{variant_NO}.\n";
                    }
                    # this is ref of majotMut info hash
                    my $majorMut_Href = $MajorMut_OrigRefAle_Href->{$major_ins_pos}->{mutHref};
                    # check mut_type
                    unless( $majorMut_Href->{mut_type} eq 'ins' ){
                        cluck_and_exit `date`."<ERROR>\tunmappable pos ($pos) 5p nearest major-mut ($major_ins_pos) is not ins, "
                                            ."but $majorMut_Href->{mut_type} in virus variant V$variant_Href->{variant_NO}.\n";
                    }
                    # check distance
                    my $majorMut_ins_seq = $majorMut_Href->{mut_seq};
                    my $majorMut_ins_len = length($majorMut_ins_seq); 
                    unless( $pos - $major_ins_pos + 1 <= $majorMut_ins_len ){
                        cluck_and_exit `date`."<ERROR>\tunmappable pos ($pos) 5p nearest major-mut ($major_ins_pos) is ins, "
                                            ."but distance out-range of insertion '$majorMut_ins_seq' (len=$majorMut_ins_len) in virus variant V$variant_Href->{variant_NO}.\n";
                    }
                    # reocrd the original modification
                    my $orig_mod;
                    # snp? ins? or del?
                    my $update_ins_seq;
                    ## for 'snp', replace the base in the majorMut-ins-seq
                    if( $mut_type eq 'snp' ){
                        $update_ins_seq = uc(   substr($majorMut_ins_seq, 0, $pos-$major_ins_pos)     # 5p-part of majorMut-ins
                                              . $mut_seq                                              # the snv
                                              . substr($majorMut_ins_seq, $pos-$major_ins_pos+1)   ); # 3p-part of majorMut-ins
                        $orig_mod = "$major_ins_orig_pos,ins,$update_ins_seq";
                    }
                    ## for 'ins', add the ins-seq in the majorMut-ins-seq
                    elsif( $mut_type eq 'ins' ){
                        $update_ins_seq = uc(   substr($majorMut_ins_seq, 0, $pos-$major_ins_pos+1)    # 5p-part of majorMut-ins
                                              . $mut_seq                                               # new ins_seq
                                              . substr($majorMut_ins_seq, $pos-$major_ins_pos+1)    ); # 3p-part of majorMut-ins
                        $orig_mod = "$major_ins_orig_pos,ins,$update_ins_seq";
                    }
                    ## for 'del'
                    ## 1) compare length of minorMut-del and majorMut-ins
                    elsif( $mut_type eq 'del' ){
                        # the part of majorMut-ins-seq before minorMut-del-pos
                        $update_ins_seq = uc( substr($majorMut_ins_seq, 0, $pos-$major_ins_pos) );
                        # compare del-seq and remained ins-seq
                        my $del_len = length($mut_seq);
                        my $ins_rm_len = $majorMut_ins_len - ( $pos - $major_ins_pos );
                        # majorMut-ins-tail remains
                        if( $ins_rm_len > $del_len ){
                            $update_ins_seq .= uc( substr($majorMut_ins_seq, $pos-$major_ins_pos+$del_len) );
                            # insertion
                            $orig_mod = "$major_ins_orig_pos,ins,$update_ins_seq";
                        }
                        # just delete majorMut-ins-tail
                        elsif( $ins_rm_len == $del_len ){
                            # insertion, next del follow
                            $orig_mod = "$major_ins_orig_pos,ins,$update_ins_seq";
                        }
                        # minorMut-del covers all majorMut-ins, and also has extendings to 3-prime.
                        elsif( $ins_rm_len < $del_len ){
                            # insertion, next del follow
                            $orig_mod = "$major_ins_orig_pos,ins,$update_ins_seq";
                            # here add new minor-mut;
                            my $del_1st_pos = $major_ins_pos + $majorMut_ins_len;
                            # check mapped orig_pos in orig-Vref
                            unless( exists($Mut2Orig_PosMapping_Href->{$del_1st_pos}) ){
                                cluck_and_exit `date`."<ERROR>\textended deletion pos ($del_1st_pos) followed major-mut (ins:$major_ins_pos) can not find mapping original pos "
                                                    ."in virus variant V$variant_Href->{variant_NO}.\n";
                            }
                            my $del_1st_orig_pos = $Mut2Orig_PosMapping_Href->{$del_1st_pos};
                            # try to see whether overlap with other followed majorMut
                            my @del_seq = split //, substr($mut_seq, $ins_rm_len);
                            # already delete the base after majorMut-ins
                            my $orig_del_seq = uc( substr($OrigVrefSeqStr, $del_1st_orig_pos-1, 1) );
                            # already base on this pos, discard this base (may-not-orig-allele)
                            shift @del_seq;
                            # check overlap with other followed majorMut
                            my $judge_pos = $del_1st_pos + 1;
                            my $Poss_OrigMod_ToPhase = '';
                            &extended_del_to_followed_reg({
                                del_seq_Aref => \@del_seq,
                                MajorMut_OrigRefAle_Href => $MajorMut_OrigRefAle_Href,
                                Mut2Orig_PosMapping_Href => $Mut2Orig_PosMapping_Href,
                                judge_pos => $judge_pos,
                                orig_del_seq_Sref => \$orig_del_seq,
                                orig_vref_str_Sref => \$OrigVrefSeqStr,
                                Poss_OrigMod_ToPhase_Sref => \$Poss_OrigMod_ToPhase
                            });
                            # orig modification updates
                            # $Poss_OrigMod_ToPhase may be empty
                            $orig_mod .= '|' . "$del_1st_orig_pos,del,$orig_del_seq" . $Poss_OrigMod_ToPhase;
                        }
                    }
                    # add the update insertion into the 'other_mutation' at majorMut-ins site.
                    my $mut_id_ToRecord = "$pos,$mut_id";
                    $majorMut_Href->{other_mutation}->{$mut_id_ToRecord}
                                                = {
                                                   mut_pos => $pos,
                                                   mut_type => $mut_type,
                                                   mut_seq => uc($mut_seq),
                                                   mut_id => $mut_id,
                                                   mut_sum_sup => $sum_sup,
                                                   mut_plus_sup => $plus_sup,
                                                   mut_minus_sup => $minus_sup,
                                                   orig_mod => $orig_mod
                                                  };
                }
            }
        }
    }

    # inform
    stout_and_sterr `date`."[INFO]:\tSummarize all credible mutations from final VCF for V$variant_Href->{variant_NO} virus variant OK.\n";
}

#--- for deletion, check overlap with other followed majorMut
sub extended_del_to_followed_reg{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $del_seq_Aref = $Option_Href->{del_seq_Aref};
    my $MajorMut_OrigRefAle_Href = $Option_Href->{MajorMut_OrigRefAle_Href};
    my $Mut2Orig_PosMapping_Href = $Option_Href->{Mut2Orig_PosMapping_Href};
    my $judge_pos = $Option_Href->{judge_pos};
    my $orig_del_seq_Sref = $Option_Href->{orig_del_seq_Sref};
    my $orig_vref_str_Sref = $Option_Href->{orig_vref_str_Sref};
    my $Poss_OrigMod_ToPhase_Sref = $Option_Href->{Poss_OrigMod_ToPhase_Sref};

    while( scalar(@$del_seq_Aref) != 0 ){
        # this judge_pos has majorMut
        if( exists($MajorMut_OrigRefAle_Href->{$judge_pos}) ){
            my $judge_majorMut_orig_seq = $MajorMut_OrigRefAle_Href->{$judge_pos}->{mut_seq};
            my $judge_majorMut_mut_type = $MajorMut_OrigRefAle_Href->{$judge_pos}->{mutHref}->{mut_type};
            # which majorMut-type?
            if( $judge_majorMut_mut_type eq 'ins' ){
                my $majorInsLen = length($judge_majorMut_orig_seq);
                my $remainDelLen = scalar(@$del_seq_Aref);
                # majorMut-ins is in this deletion
                if( $majorInsLen <= $remainDelLen ){
                    # discard the full insertion
                    shift @$del_seq_Aref for ( 1 .. $majorInsLen );
                    # update judge_pos(mod-Vref) to skip the inserted part
                    $judge_pos += $majorInsLen;
                }
                # only part of insertion is in this deletion
                else{
                    # empty remained del_seq
                    # AND, no need to update $judge_pos
                    @$del_seq_Aref = ();
                    # here, the majorIns has been modified to smaller insertion.
                    # i.e., this is a minorIns to record.
                    my $this_minorInsSeq = uc( substr($judge_majorMut_orig_seq, $remainDelLen) );
                    # the pos before this majorIns is what we need to record this minorIns
                    my $this_majorInsStPos = $judge_pos - 1;
                    my $this_minorIns_OrigPos = $Mut2Orig_PosMapping_Href->{$this_majorInsStPos};
                    # tell back this minorIns, must have '|' meant 'to phase'
                    $$Poss_OrigMod_ToPhase_Sref = '|' . "$this_minorIns_OrigPos,ins,$this_minorInsSeq";
                }
            }
            elsif( $judge_majorMut_mut_type eq 'snp' ){
                # reocrd the orig-base at this pos
                $$orig_del_seq_Sref .= uc( $judge_majorMut_orig_seq );
                # discard this pos from del_seq_Aref
                shift @$del_seq_Aref;
                # to next judge_pos(mod-Vref)
                $judge_pos ++;
            }
            elsif( $judge_majorMut_mut_type eq 'del' ){
                # firstly, record this base before the del-part
                # use orig-allele from orig-vref-seq
                my $judge_origpos = $Mut2Orig_PosMapping_Href->{$judge_pos};
                $$orig_del_seq_Sref .= uc( substr($$orig_vref_str_Sref, $judge_origpos-1, 1) );
                # already base on this pos, discard this base (may-not-orig-allele)
                shift @$del_seq_Aref;
                # then, concatenate the majorMut_del_seq
                $$orig_del_seq_Sref .= uc( $judge_majorMut_orig_seq );
                # to next judge_pos(mod-Vref)
                $judge_pos ++;
            }
        }
        # this is just pos between majorMuts
        else{
            # just record the orig-allele (already in del_seq_Aref)
            $$orig_del_seq_Sref .= uc( shift @$del_seq_Aref );
            # to next judge_pos(mod-Vref)
            $judge_pos ++;
        }
    }
}

#--- output modification list ---
sub output_modification{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_Href = $Option_Href->{variant_Href};

    my @theme = qw/ VirusVariant
                    OriginalVirus
                    OriginalPos
                    Ref_Allele
                    Ref_Supp
                    MajorAlt_Type
                    MajorAlt_Allele
                    MajorAlt_Supp
                    In_ModVirus
                    MutVrefPos
                    MinorAltModVref(OVref)
                  /;
    open (MODF,Try_GZ_Write($variant_Href->{Mods})) || die "fail write $variant_Href->{Mods}: $!\n";
    # theme
    print MODF '#' . join("\t", @theme ) . "\n";
    for my $orig_pos (sort {$a<=>$b} keys %{ $variant_Href->{Mutation} }){
        # VirusVariant
        print MODF        $variant_Href->{RefSeg};
        # OriginalVirus
        print MODF "\t" . $variant_Href->{OrigRefSeg};
        # OriginalPos
        print MODF "\t" . $orig_pos;

        my $Mut_Href = $variant_Href->{Mutation}->{$orig_pos};
        # OrigV Ref_Allele
        print MODF "\t" . uc($Mut_Href->{ref_allel});
        # OrigV Ref_Supp
        if( $variant_Href->{ModCompBool} ){
            print MODF "\t" . join(',', $Mut_Href->{ref_sum_sup},
                                        $Mut_Href->{ref_plus_sup},
                                        $Mut_Href->{ref_minus_sup} );
        }
        else{
            print MODF "\t" . '-';
        }

        # majorMut
        if( $Mut_Href->{major_mut_bool} ){
            # MajorAlt_Type
            print MODF "\t" . uc($Mut_Href->{mut_type});
            # MajorAlt_Allele
            print MODF "\t" . uc($Mut_Href->{mut_seq});
            # MajorAlt_Supp
            if( $variant_Href->{ModCompBool} ){
                print MODF "\t" . join(',', $Mut_Href->{mut_sum_sup},
                                            $Mut_Href->{mut_plus_sup},
                                            $Mut_Href->{mut_minus_sup} );
            }
            else{
                print MODF "\t" . '-';
            }
            # is it MajorAlt?
            print MODF "\t" . 'Y';
        }
        # minorMut(s) are all in key 'other_mutation'
        else{
            # MajorAlt_Type
            print MODF "\t" . '-';
            # MajorAlt_Allele
            print MODF "\t" . '-';
            # MajorAlt_Supp
            print MODF "\t" . '-';
            # is it MajorAlt?
            print MODF "\t" . '-';
        }

        # Pos on V-ref modified
        # even no majorMut exists at this pos, the related mod-vref pos should show here if possible
        if( $variant_Href->{ModCompBool} ){
            print MODF "\t" . ($Mut_Href->{modVrefpos} || '-');
        }
        else{
            print MODF "\t" . '-';
        }

        # OtherAlt_Info
        my $other_mut_Href = $Mut_Href->{other_mutation};
        my @other_mut_id = sort { $other_mut_Href->{$b}->{mut_sum_sup}
                                  <=>
                                  $other_mut_Href->{$a}->{mut_sum_sup}
                                } keys %$other_mut_Href;
        if(    $variant_Href->{ModCompBool}
            && scalar(@other_mut_id)!=0
          ){
            my @other_mut_info;
            for my $mut_id (@other_mut_id){
                my $mut_Href = $other_mut_Href->{$mut_id};
                push @other_mut_info, join(',', $mut_Href->{mut_pos},
                                                $mut_Href->{mut_type},
                                                $mut_Href->{mut_seq},
                                                $mut_Href->{mut_sum_sup},
                                                $mut_Href->{mut_plus_sup},
                                                $mut_Href->{mut_minus_sup},
                                                "($mut_Href->{orig_mod})"  );
            }
            print MODF "\t" . join(';',@other_mut_info);
        }
        else{
            print MODF "\t" . '-'; # MajorAlt_Allele
        }
        print MODF "\n";
    }
    close MODF;
}

#--- resort the virus variants by average-depth
sub resort_virus_variant_by_depth{
    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    my %depth2vNO;
    # gather virus_variants' average-depth
    for my $virus_variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
        my $average_depth = $this_variant_Href->{AverageDepth};
        push @{ $depth2vNO{$average_depth} }, $virus_variant_NO;
    }

    my @resorted_vNO;
    # resort virus_variant NO. by average_depth
    for my $average_depth (sort {$b<=>$a} keys %depth2vNO){
        push @resorted_vNO, sort {$a<=>$b} @{ $depth2vNO{$average_depth} };
    }

    my $exchange_bool = 0;
    # check whether needs to exchange files of virus variants
    $exchange_bool ||= ( $resorted_vNO[$_] != ($_+1) ) for ( 0 .. $#resorted_vNO );
    unless( $exchange_bool ){
        return;
    }

    # exchange files!
    my $exchange_temp_folder = File::Spec->catfile( $V_Href->{Virus_variant_folder}, 'ExchangeTempSpace' );
    `rm -rf $exchange_temp_folder`;
    `mkdir -p $exchange_temp_folder`;
    ## move files to temp space and rename!
    for my $i ( 0 .. $#resorted_vNO ){
        my $aim_vNO = $i + 1;
        my $sor_vNO = $resorted_vNO[$i];
        # no need to exchange
        next if($aim_vNO == $sor_vNO);
        # change results' name
        my $sor_variant_Href = $Variant_Info_Href->{$sor_vNO};
        my $sor_variant_folder = $sor_variant_Href->{workfolder};
        ## change file name
        for my $filekey (qw/ Fa mkupFa Mods VCF finalBam PosMapList /){
            my $sor_file = $sor_variant_Href->{$filekey};
            # change name
            ( my $aim_file_basename = basename($sor_file) ) =~ s/V$sor_vNO/V$aim_vNO/;
            my $aim_file = File::Spec->catfile( $sor_variant_folder, $aim_file_basename );
            `mv $sor_file $aim_file` if( -e $sor_file );
            ## bai file
            `mv $sor_file.bai $aim_file.bai` if($filekey eq 'finalBam');
        }
        ## 1) move workfolder to temp space and 2) change folder name
        ( my $aim_folder_basename = basename($sor_variant_folder) ) =~ s/V$sor_vNO/V$aim_vNO/;
        my $aim_folder = File::Spec->catfile( $exchange_temp_folder, $aim_folder_basename );
        `mv $sor_variant_folder $aim_folder`;
        # record the exchange
        my $aim_variant_Href = $Variant_Info_Href->{$aim_vNO};
        $aim_variant_Href->{ExchangeSorVno} = $sor_vNO;
    }
    ## move all workfolder out of temp space
    `mv $exchange_temp_folder/V* $V_Href->{Virus_variant_folder}/`;
    `rm -rf $exchange_temp_folder`;

    # inform
    stout_and_sterr `date`."[INFO]:\tresort virus variant by average-depth ok.\n";
}

#--- write virus variant report ---
sub write_virus_variant_report{
    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

    my @theme = (
                    'NO',
                    'DealCycleTime',
                    'OriginalVirus',
                    'VariantVirus',
                    'OrigLength',
                    'nonN_Length',
                    'CricleOrNot',
                    'ToExtLength',
                    'Coverage',
                    'AverageDepth',
                    'Modifications',
                    'ModComp',
                    'RefFasta',
                    'MakeupFasta',
                    'AlignBam',
                    'ModDetails'
                );

    open (VVRP,Try_GZ_Write($V_Href->{Virus_variant_report})) || die "fail write $V_Href->{Virus_variant_report}: $!\n";
    print VVRP '#' . join("\t", @theme) . "\n"; # theme
    # each virus variant
    for my $virus_variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        my $this_variant_Href = $Variant_Info_Href->{$virus_variant_NO};
        my $real_variant_Href = (exists $this_variant_Href->{ExchangeSorVno}) ?
                                      $Variant_Info_Href->{ $this_variant_Href->{ExchangeSorVno} }
                                    : $this_variant_Href;
        # basic information
        print VVRP        $virus_variant_NO;
        print VVRP "\t" . $real_variant_Href->{CycleTime};
        print VVRP "\t" . $real_variant_Href->{OrigRefSeg};
        print VVRP "\t" . $real_variant_Href->{RefSeg};
        print VVRP "\t" . $Virus_Info_Href->{ $real_variant_Href->{OrigRefSeg} }->{length};
        print VVRP "\t" . $real_variant_Href->{NonN_length};
        print VVRP "\t" . $real_variant_Href->{Cricle};
        print VVRP "\t" . $real_variant_Href->{Extend_length};
        # coverage
        print VVRP "\t";
        for my $depth (sort {$a<=>$b} keys %{ $real_variant_Href->{Coverage} }){
            print VVRP $depth . 'X:' . ($real_variant_Href->{Coverage}->{$depth} * 100) . '%;';
        }
        # average depth
        print VVRP "\t" . $real_variant_Href->{AverageDepth};
        # mutations
        my %type2count = (sum=>0, snp=>0, mnp=>0, ins=>0, del=>0);
        for my $pos (sort {$a<=>$b} keys %{ $real_variant_Href->{Mutation} }){
            next unless($real_variant_Href->{Mutation}->{$pos}->{major_mut_bool});
            my $type = $real_variant_Href->{Mutation}->{$pos}->{mut_type};
            $type2count{$type} ++;
            $type2count{sum} ++;
        }
        print VVRP "\t";
        print VVRP $_ . ':' . $type2count{$_} . ';' for ('sum', 'snp', 'ins', 'del'); # 'mnp', 
        # whether modifications are completed for construct v-ref
        print VVRP "\t" . ($real_variant_Href->{ModCompBool} ? 'Y' : 'N');
        # Fa
        print VVRP "\t" . File::Spec->abs2rel( $this_variant_Href->{Fa}, $V_Href->{Virus_variant_folder} );
        # makeup Fa
        print VVRP "\t" . File::Spec->abs2rel( $this_variant_Href->{mkupFa}, $V_Href->{Virus_variant_folder} );
        # Bam
        print VVRP "\t" . File::Spec->abs2rel( $this_variant_Href->{finalBam}, $V_Href->{Virus_variant_folder} );
        # Modification details
        print VVRP "\t" . File::Spec->abs2rel( $this_variant_Href->{Mods}, $V_Href->{Virus_variant_folder} ) . "\n";
    }
    close VVRP;

    # inform
    stout_and_sterr `date`."[INFO]:\tgenerate report of virus variant ok.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
