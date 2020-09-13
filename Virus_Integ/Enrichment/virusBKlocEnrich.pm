package FuseSV::Virus_Integ::Enrichment::virusBKlocEnrich;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use BioFuse::Util::Sys qw/ file_exist /;
use FuseSV::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              virus_bkLoc_enrich
              test_locVirus_para
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Integ::Enrichment::virusBKlocEnrich';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        virus_bkLoc_enrich
                        test_locVirus_para
                     /;

#--- do enrichment test on host break-point distribution ---
sub virus_bkLoc_enrich{
}

#--- test on required options ---
sub test_locVirus_para{
    return  (   $V_Href->{enrich_mode} eq 'loc:virus'
             && (   !file_exist(filePath=>$V_Href->{enrichRegion_bed})
                 || !file_exist(filePath=>$V_Href->{virusRef_list})
                )
            );
}

#--- 
1; ## tell the perl script the successful access of this module.

























