package FuseSV::Virus_Integ::Enrichment::IntegEnrichCheck;

use strict;
use warnings;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::Enrichment::hostBKlocEnrich qw/ test_locHost_para host_bkLoc_enrich /;
use FuseSV::Virus_Integ::Enrichment::virusBKlocEnrich qw/ test_locVirus_para virus_bkLoc_enrich /;
use FuseSV::Virus_Integ::Enrichment::BKcharMHEnrich qw/ test_charMH_para bkCharMH_enrich /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              check_vinteg_enrichment
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Integ::Enrichment::IntegEnrichCheck';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-23';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        check_vinteg_enrichment
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} integ_enrich <[Options]>
     
     Options:

       # Inputs and Outputs #
        -vitg    [s]  virus integration list obtain from 'virus_integ'. <required>
        -odir    [s]  folder to store results and temporary files. <required>

       # DataBase #
        -href    [s]  host reference fasta file with fai file. <required>
        -nonN    [s]  bed file, nonN region of host reference. <required>
                       Note: required by 'loc:host' and 'char:MH' mode.
        -bed     [s]  bed file, concerned region.
                       Note: 1) required by 'loc:host' and 'loc:virus' mode.
                             2) format: chr3  1000  2000  [nameOfRegion]
        -vrefl   [s]  list of virus reference.
                       Note: required by 'loc:virus' and 'char:MH' mode.

       # Options #
        -mode    [s]  mode of enrichment. [loc:host]
                       Note: use one of 'loc:host', 'loc:virus' 'char:MH'.
        -simur   [s]  round of simulation for expected. [1E5]
        -simuc   [s]  count of case to simulate in each round. [same as observed]
        -dist    [s]  distance to virus integration. [0]

        -h|help      Display this help info.

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
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ sampleToVitg_list => undef ],
            [ nonN_bed => undef ],
            [ outdir => undef ],
            [ enrichRegion_bed => undef ],
            [ hostRef_fa => undef ],
            [ virusRef_list => undef ],

            # options
            [ enrich_mode => 'loc:host' ],
            [ vItg_distance => [] ], # soon reset to Hash
            [ simulation_round => 100000 ],
            [ simulation_count => undef ],
            [ positiveRateUnit => 0.05 ],

            # setting
            [ modeToFunc => { 'loc:host' => \&host_bkLoc_enrich,
                              'loc:virus' => \&virus_bkLoc_enrich,
                              'char:MH' => \&bkCharMH_enrich
                            } ],

            # intermediate variants
            [ hostRefseg => {} ],
            [ hostNonNRegion => [] ],
            [ ConcernRegion => {} ],
            [ winSize => 10000 ],
            [ SampleToVitg => {} ],
            [ realVitgCount => 0 ],
            [ vItgLocCount => {} ],
            [ positiveRate => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['sampleToVitg_list'],
                                  ['nonN_bed'],
                                  ['outdir'],
                                  ['enrichRegion_bed'],
                                  ['hostRef_fa'],
                                  ['virusRef_list']   ] ]
        );
}

#--- get options from command line
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-vitg:s"   => \$V_Href->{sampleToVitg_list},
        "-odir:s"   => \$V_Href->{outdir},
        "-nonN:s"   => \$V_Href->{nonN_bed},
        "-bed:s"    => \$V_Href->{enrichRegion_bed},
        "-href:s"   => \$V_Href->{hostRef_fa},
        "-vrefl:s"  => \$V_Href->{virusRef_list},
        # options
        "-mode:s"   => \$V_Href->{enrich_mode},
        "-simur:s"  => \$V_Href->{simulation_round},
        "-simuc:s"  => \$V_Href->{simulation_count},
        "-dist:s"   => \@{$V_Href->{vItg_distance}}, # soon reset to Hash
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{sampleToVitg_list})
             || !file_exist(filePath=>$V_Href->{hostRef_fa})
             || !defined $V_Href->{outdir}   || !-d $V_Href->{outdir}
             || !exists $V_Href->{modeToFunc}->{$V_Href->{enrich_mode}}
             || test_locHost_para
             || test_locVirus_para
             || test_charMH_para
            );
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub check_vinteg_enrichment{
    # preparation
    &prepare;
    # specific enrichment work
    &{ $V_Href->{modeToFunc}->{$V_Href->{enrich_mode}} };
}

#--- prepare work ---
sub prepare{
    # hash virus distance
    my @regionDist = map {if(/(.+):(\d+)$/){($1,$2)}} grep /:\d+$/, @{$V_Href->{vItg_distance}};
    my @globalDist = sort {$b<=>$a} grep {/^\d+$/ && $_>0} @{$V_Href->{vItg_distance}};
    $V_Href->{vItg_distance} =  scalar @regionDist
                              ? {@regionDist}
                              : (  scalar @globalDist
                                 ? {'_ALL_' => $globalDist[0]}
                                 : {}
                                );
    if(scalar @regionDist || scalar @globalDist){
        my $tmpHf = $V_Href->{vItg_distance};
        stout_and_sterr "[INFO]\tBelow region(s) are set additional distance to virus integration.\n"
                             .join('', map{"\t$_\t$tmpHf->{$_}\n"} sort keys %$tmpHf);
    }
}

#--- 
1; ## tell the perl script the successful access of this module.




































