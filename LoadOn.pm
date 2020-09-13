package FuseSV::LoadOn;

use strict;
use warnings;
use Parallel::ForkManager;
use BioFuse::Util::Log qw/ cluck_and_exit /;
use FindBin qw/ $RealBin /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              $V_Href
              forkSetting
            /;
              # load_variants_dict
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::LoadOn';
#----- version --------
$VERSION = "0.49";
$DATE = '2018-07-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#------ variants -------
our $V_Href = {

    #++++++++++++++#
    # Main Version #
    #++++++++++++++#

    MainName => 'FuseSV.pl',
    Version => '0.28',
    Date => '2020-09-12',
    AUTHOR => 'Wenlong Jia',
    EMAIL => 'wenlongkxm@gmail.com',

    # functions
    RealBin => $RealBin,
    func => {},
    command => undef,
    argv_Aref => undef,
    run_mode => 0,

    ## manual
    HELP => 0,
    HELP_INFO => {},

    # for debug to keep some key intermediate folders/files
    in_debug => 0,

    # software pre-defined
    Tool_Tag => 'FS', # for readid-prefix distinguish

    #+++++++++++++++++++++++++++++#
    # intermediate file/folder(s) #
    #+++++++++++++++++++++++++++++#

    ## virus database files
    virus_db_folder => undef,
    virus_orig_genome_fa => '',
    virus_info_doc => '',
    virus_orig_gene_psl => '',
    virus_orig_trans_psl => '',
    ## folder and files
    CandReadsFolder => '', # SampleID.CandReads
    Cand_RgInfoReport => '', # RG.info

    virusInitAlignFolder => '', # VirusInitAlign
    virusInitAlignReport => '', # virus.initAlign.report
    virusInitAlignBamFolder => '', # InitAlignBAM, under virusInitAlignFolder
    virusInitAlignRefFolder => '', # InitAlignREF, under virusInitAlignFolder

    Virus_variant_folder => '', # VirusVariants

    Prepare_Integ_folder => 'Prepare_DetectInteg',

    ForVitgPrefix => 'For_Detect_Integ', # file prefix under Prepare_Integ_folder

    vIteg_workFolder => 'Virus_Integration',

    #+++++++++++++++++#
    # general options #
    #+++++++++++++++++#

    ## maximum read length from bam
    max_rlen => 0,
    ## multiple fork operations
    fork_number => 0,
    ## fasta
    line_base => 50,
    ## for visual
    general_SVG_obj => undef,
    general_FontFam => 'Arial',
    general_FontFam_Bold => 'Arial Bold',
    general_FontSize => 12,
    general_BaseX => 150,
    general_BaseY => 150
};

#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_variants_dict
                        forkSetting
                     /;

#--- fork manager setting ---
sub forkSetting{
    if($V_Href->{fork_number} > 1){
        my $pm = new Parallel::ForkManager($V_Href->{fork_number});
        $pm->run_on_finish(
            sub {
                my ($pid, $exit_code, $tag, $signal) = @_;
                cluck_and_exit "<ERROR>\t$tag child-process($pid) failed.\n" if $exit_code || $signal;
            }
        );
        return ($pm, 1);
    }
    else{
        return (undef, 0);
    }
}

#--- load variant dict Href ---
sub load_variants_dict{
    return $V_Href;
}

#--- 
1; ## tell the perl script the successful access of this module.
