package FuseSV::Virus_Type::VirusTypeMain;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::BioInfo::Objects::Sample::Sample_OB;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Virus_Type::CandReads qw/ CandReads /;
use FuseSV::Virus_Type::VirusInitAlign qw/ VirusInitAlign /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
                confirmViursTypeVariants
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Type::VirusTypeMain';
#----- version --------
$VERSION = "0.03";
$DATE = '2019-06-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        confirmViursTypeVariants
                        prepare
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} virus_type <[Options]>
     
     Options:

       # Inputs and Outputs #
        -bam    [s]  the indexed-bam alignment file. <required>
                      Assume to be in Sanger-quality format (column NO.11).
        -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
                      only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
        -odir   [s]  directory to store outputs, must exist already. <required>
        -host   [s]  the .fai file of only host genome reference. <required>
                      such as pure human reference fai file, created by samtools faidx.
        -vpf    [s]  virus database folder. <required>
                      You can make it via function 'make_virusdb'.
                      OR Download these widely used from FuseSV official website, e.g., HPV.

       # Software Required #
        -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
        -stl    [s]  SamTools, minimum version: 1.3. <required>
        -btl    [s]  bcftools, minimum version: 1.3. <required>
        -gatk   [s]  GATK jar file, minimum version: 3.5. <required>
        -pic    [s]  Picard jar file, minimum version: 2.1.0. <required>
        -java   [s]  Java compiler, minimum version 1.8. ['java']
        -blta   [s]  blastall of blast. <required>
        -fordb  [s]  formatdb of blast. <required>

       # General Options #
        -st_step [i]  start wrokflow from certain step. [1]
        -ed_step [i]  stop workflow at certain step. [$V_Href->{totalStepNO}]
                       Note: step-NO. list:  1: cand_reads; 2: cand_virus; 3: get_variant; 4: prep_integ;

       # Candidate Reads Option #
        -minmq  [i]  PE reads whose either mapQ lower than INT will be extract for virus analysis. [30]
                      The properly mapping PE-reads (i.e., 0x2 in flag) will be all ignored.
        -minsc  [i]  PE reads whose either softclip length larger than INT will be extracted. [10]
        -minlen [i]  minimum candidate read length. [50]
        -disdup      discard duplicated reads with SAM-flag 0x400. [Disabled]
        -bartag [s]  TAG in bam to extract barcode info, such as 'BX:Z:' of 10X technology. [none]
                      The extracted barcode will be concatenated with paired-end ID.

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

        # Prepare Integration Options #
        -skipitg     Skip the module 'PrepareForIntegDetection' (when just want 'VirusVariants' results). [Disabled]
        -simple      to get simplified bam for quicker but maybe more false-negative integration detection. [Disabled]
        -vno    [s]  to select virus variants (its NO.) for virus integration detection. [1]
                      default value is 1, which means only use the FIRST virus variant;
                      set as 'all' to use ALL virus variants;
                      set as '1,2' to use FIRST and SECOND variants (if it has, automatically judgement).
         ** Note: options below just work on virus variants (besides NO.1) selected by option '-vno'.
        -midp   [i]  minimum average depth of viral genome. [5]
        -micr   [f]  minimum covered ratio of viral genome. [0.5]

        -h|help   Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
        # -tstag  [s]  add ONE uppercase prefix to read-id to distinguish reads from different tissues. [Disabled]
        # -mxnc   [i]  maximum 'N-string' count allowed in each virus variant genome. [5]
        # -minl   [i]  minimum 'N-string' length to count in each virus variant genome. [5]
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{

    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $only_load_containers = $parm{only_load_containers} || 0;

    $V_Href->{ $_->[0] } = $_->[1] for
        (
            #------------#
            # containers #
            #------------#
            # id -> refsegOB
            # see BioFuse::BioInfo::Objects::Segment::RefSeg_OB
            [ hostRefSeg => {} ],
            ## refseg
            ## refseg -> genebank_id
            ## refseg -> vref_no
            ## refseg -> length
            ## refseg -> extend_len
            ## refseg -> circle
            ## refseg -> vref_fa
            ## refseg -> init_mapbam
            ## refseg -> ncbi_url
            # [ Virus_Info_Href => {} ],
            # id(subtype) -> virusOB
            # see FuseSV::Objects::Virus_OB
            [ virusRefSeg => {} ],
            # id -> RG_OB
            # see BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB
            [ RG => {} ],
            ## virus_variant_NO (vNO)
            ## vNO -> variant_NO
            ## vNO -> workfolder
            ## vNO -> OrigRefSeg
            ## vNO -> RefSeg
            ## vNO -> Orig_length
            ## vNO -> NonN_length
            ## vNO -> Extend_length
            ## vNO -> Cricle
            ## vNO -> Fa
            ## vNO -> Fa_md5 -> { md5 => cycle_time }
            ## vNO -> mkupFa
            ## vNO -> mkdupBam
            ## vNO -> finalBam
            ## vNO -> FinalIndexVref, this is intermediate file
            ## vNO -> VCF
            ## vNO -> CycleTime
            ## vNO -> Coverage -> {>=1X, >=4X, >=10X, >=30X, >=50X, >=min_depth}
            ## vNO -> AverageDepth
            ## vNO -> Mutation -> {pos} -> {  modVrefpos,
            ##                                ref_allel, ref_sum_sup, ref_plus_sup, ref_minus_sup,
            ##                                mut_type, mut_seq, mut_id,
            ##                                mut_sum_sup, mut_plus_sup, mut_minus_sup,
            ##                                major_mut_bool=>0/1, other_mutation=>{} }
            ## vNO -> Mods
            ## vNO -> ModCompBool
            ## vNO -> ExchangeAimVno
            [ Variant_Info_Href => {} ],
            ## GenePSL of virus
            ## TransPSL of virus
            ## Variant_Extend_Href
            ## GeneList
            ## MakeupRef
            ## BamForInteg
            ## ParaTxt
            ## AllowVNO
            [ Prepare_For_Integ => {} ]
        );

    return if $only_load_containers;

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
            [ virusAnaFolder => undef ],
            [ bam => undef ],
            [ sample_ID => undef ],
            # database
            [ host_fai => undef ], # such as pure human reference fai file, created by samtools faidx
            # software
            [ bwa => undef ], # bwa 0.7.9
            [ samtools => undef ], # samtools 1.3
            [ bcftools => undef ], # bcftools 1.3
            [ java => 'java' ], # java 1.8
            [ gatk_jar => undef ], # gatk 3.5
            [ picard_jar => undef ], # picard 2.1.0
            [ blastall => undef ],
            [ formatdb => undef ],
            # other options
            [ mem_for_jvm => '2g' ],
            # 1: cand_reads; 2: cand_virus; 3: get_variant; 4: prep_integ;
            [ stepToStart => 1 ],
            [ stepToStop => 4 ],
            [ totalStepNO => 4 ],
            # extract candidate fastq #
            [ min_mapQ => 30 ],
            [ min_Slen => 10 ],
            [ min_Rlen => 50 ],
            [ discard_DupReads => 0 ], # discard dup-reads when extract smart-reads from smart-bam(alway-has-dup-reads)
            # [ tissue_tag => undef ], # for distinguish reads from different tissues
            [ barcSAMtag => undef ], # such as BX:Z:'' of 10x reads
            [ barcPidPrf => ':__BARC__:' ], # barcode prefix in pid
            # determine virus variant, construct virus genome #
            [ max_vNO => 3 ], # allow such number of virus variant to reconstruct.
            [ min_coverage => 0.5 ], # minimum coverage of one virus to count into the candidate virus(es) group for reads-counting
            [ count_unproper_map => 0 ], # count the initial mapped reads on each virus sub-type
            [ count_non_supple => 0 ],
            [ count_min_MAQ => 30 ],
            [ count_max_Sratio => 0 ],
            [ Nsort_MaxStore => 2000 ],
            # use_refallele => 0, # discarded
            [ min_depth => 1 ],
            [ both_strd => 1 ], # generally, we require this non-zero; but for virus with linear genome, should set as 0 for their bilateral mutations.
            [ min_mut_reads => 2 ],
            [ min_mut_freq => 0.1 ],
            [ I16_min_mean_baseQ => 20 ],
            [ I16_min_mean_mapQ => 20 ],
            [ I16_min_mean_tailDist => 5 ],
            [ adjacent_indel_min_dist => 5 ], # for adjacent indel exclusion
            [ consider_IDV_bool => 0 ], # for INDEL, IDV provides more precise supports.
            [ max_cycle_time => 10 ], # maximum cycle time of reconstruction work for each virus variant
            [ min_depth_ratio_of_V1 => 0 ], # from NO.2 variant, min depth ratio of V1's average depth
            [ cigar_RepeatMut_FilterFlankSize => 10 ], # cigar mutation control, this is for one side
            [ cigar_RepeatMut_FilterAlertCount => 3 ], # cigar mutation control, this is for one side
            [ Virus_variant_report => '' ],
            # prepare virus integration #
            [ Skip_VirusItegPrep => 0 ],
            [ Get_simple_BamForInteg => 0 ],
            [ VariantNO_For_Analysis => 1 ],
            [ MinDepth_For_IntegAna => 5 ],
            [ MinCoverRatio_For_IntegAna => 0.5 ],
            [ MaxAllow_Nstr_Count => 5 ],
            [ Nstr_MinLength => -1 ],

            # containers
            [ SampleOB => undef ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['bam'],
                                  ['virusAnaFolder'],
                                  ['host_fai'],
                                  ['virus_db_folder'],
                                  ['bwa'],
                                  ['samtools'],
                                  ['bcftools'],
                                  ['java'],
                                  ['gatk_jar'],
                                  ['picard_jar'],
                                  ['blastall'],
                                  ['formatdb']             ] ]
        );
}

#--- get options from command line
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-bam:s"    => \$V_Href->{bam},
        "-id:s"     => \$V_Href->{sample_ID},
        "-odir:s"   => \$V_Href->{virusAnaFolder},
        "-host:s"   => \$V_Href->{host_fai},
        "-vpf:s"    => \$V_Href->{virus_db_folder},
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
        ## steps
        "-st_step:i"=> \$V_Href->{stepToStart},
        "-ed_step:i"=> \$V_Href->{stepToStop},
        ## multiple fork operations
        "-fm:i"     => \$V_Href->{fork_number}, # hidden option
        "-jvmem:s"  => \$V_Href->{mem_for_jvm}, # hidden option
        ## cand_reads for alignment on virus
        # "-tstag:s"  => \$V_Href->{tissue_tag},
        "-minmq:i"  => \$V_Href->{min_mapQ},
        "-minsc:i"  => \$V_Href->{min_Slen},
        "-minlen:i" => \$V_Href->{min_Rlen},
        "-disdup"   => \$V_Href->{discard_DupReads},
        "-bartag:s" => \$V_Href->{barcSAMtag},
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
        ## prepare for integration
        "-skipitg"  => \$V_Href->{Skip_VirusItegPrep},
        "-simple"   => \$V_Href->{Get_simple_BamForInteg},
        "-vno:s"    => \$V_Href->{VariantNO_For_Analysis},
        "-midp:i"   => \$V_Href->{MinDepth_For_IntegAna},
        "-micr:f"   => \$V_Href->{MinCoverRatio_For_IntegAna},
        "-mxnc:i"   => \$V_Href->{MaxAllow_Nstr_Count}, # hidden option
        "-minl:i"   => \$V_Href->{Nstr_MinLength}, # hidden option
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{bam})
             || !file_exist(filePath=>$V_Href->{virusAnaFolder})
             || !file_exist(filePath=>$V_Href->{host_fai})
             || !file_exist(filePath=>$V_Href->{virus_db_folder})
             || !file_exist(filePath=>$V_Href->{bwa})
             || !file_exist(filePath=>$V_Href->{samtools})
             || !file_exist(filePath=>$V_Href->{bcftools})
             || !file_exist(filePath=>$V_Href->{gatk_jar})
             || !file_exist(filePath=>$V_Href->{picard_jar})
             || !file_exist(filePath=>$V_Href->{blastall})
             || !file_exist(filePath=>$V_Href->{formatdb})
             || ($V_Href->{java} ne 'java' && !-e $V_Href->{java})
             || $V_Href->{stepToStart} < 1 || $V_Href->{stepToStart} > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < 1 || $V_Href->{stepToStop}  > $V_Href->{totalStepNO}
             || $V_Href->{stepToStop}  < $V_Href->{stepToStart}
             || (!$V_Href->{sample_ID} || $V_Href->{sample_ID} =~ /[^\w\-]/)
             # || (defined $V_Href->{tissue_tag} && $V_Href->{tissue_tag} !~ /^[A-Z]$/)
             || ($V_Href->{min_depth_ratio_of_V1} < 0 || $V_Href->{min_depth_ratio_of_V1} > 1)
             || ($V_Href->{max_vNO} <= 0 || $V_Href->{max_vNO} =~ /\D/)
             || $V_Href->{MinCoverRatio_For_IntegAna} < 0
             || $V_Href->{MinCoverRatio_For_IntegAna} > 1
             || $V_Href->{MinDepth_For_IntegAna} < 0
             # || $V_Href->{Nstr_MinLength} < 0
             || $V_Href->{MaxAllow_Nstr_Count} < 0
            );
}

#--- confirm virus type and variants ---
sub confirmViursTypeVariants{
    # preparation
    &prepare;

    # cand_reads
    CandReads;

    # initial alignment
    VirusInitAlign;
}

#--- prepare objects ---
sub prepare{
    # sample object
    $V_Href->{SampleOB} = BioFuse::BioInfo::Objects::Sample::Sample_OB->new(id=>$V_Href->{sample_ID});
    # orig bam object
    my $bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath=>$V_Href->{bam}, tag=>'origBam');
    $bam->addTool(samtools=>$V_Href->{samtools});
    $V_Href->{SampleOB}->addBam(bam=>$bam);
    # step 1, cand_reads
    $V_Href->{CandReadsFolder} = GetPath(filekey => 'CandReadsFolder');
    $V_Href->{Cand_RgInfoReport} = GetPath(filekey => 'Cand_RgInfoReport');
    # step 2, cand_virus
    $V_Href->{virusInitAlignFolder} = GetPath(filekey => 'virusInitAlignFolder');
    $V_Href->{virusInitAlignReport} = GetPath(filekey => 'virusInitAlignReport');
    $V_Href->{virusInitAlignBamFolder} = GetPath(filekey => 'virusInitAlignBamFolder');
    $V_Href->{virusInitAlignRefFolder} = GetPath(filekey => 'virusInitAlignRefFolder');

}

#--- 
1; ## tell the perl script the successful access of this module.
