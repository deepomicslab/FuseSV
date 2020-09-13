package FuseSV::Virus_Integ::DrawSegCN::ObjCountTransfMixToPureT;

use strict;
use warnings;
use List::Util qw/ min max sum /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              get_Tcell_GMpart
              get_Tcell_GMratio
              get_Ncell_GMpart
              get_Ncell_GMratio
              get_ObjCountOfPureTumorCell
              get_ObjSingleCNdepthInMixed
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ObjCountTransfMixToPureT';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-10-19';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        get_Tcell_GMpart
                        get_Tcell_GMratio
                        get_Ncell_GMpart
                        get_Ncell_GMratio
                        get_ObjCountOfPureTumorCell
                        get_ObjSingleCNdepthInMixed
                     /;

#--- get genoome meterial part of Pure Tumor Cell in Mixed cells ---
sub get_Tcell_GMpart{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $TumorCellPurity = $parm{TumorCellPurity} || 1; # default is full-tumor-cell
    my $TumorCellPloidy = $parm{TumorCellPloidy} || 2; # default is diploid

    return $TumorCellPurity * $TumorCellPloidy;
}

#--- get genoome meterial ratio of Pure Tumor Cell in Mixed cells ---
sub get_Tcell_GMratio{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $Tcell_GMpart = $parm{Tcell_GMpart};
    my $Ncell_GMpart = $parm{Ncell_GMpart};

    return $Tcell_GMpart / ($Tcell_GMpart + $Ncell_GMpart);
}

#--- get genoome meterial part of Pure Normal Cell in Mixed cells ---
sub get_Ncell_GMpart{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $TumorCellPurity  = $parm{TumorCellPurity} || 1; # default is full-tumor-cell
    my $NormalCellPloidy = $parm{NormalCellPloidy} || 2; # default is diploid

    return ( 1 - $TumorCellPurity ) * $NormalCellPloidy;
}

#--- get genoome meterial ratio of Pure Normal Cell in Mixed cells ---
sub get_Ncell_GMratio{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $Tcell_GMpart = $parm{Tcell_GMpart};
    my $Ncell_GMpart = $parm{Ncell_GMpart};

    return $Ncell_GMpart / ($Tcell_GMpart + $Ncell_GMpart);
}

#--- subtracts object's count of normal-cells from mixed object's count ---
sub get_ObjCountOfPureTumorCell{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $ObjCount = $parm{ObjCount};
    my $MixedCellMeanDepth = $parm{MixedCellMeanDepth};
    my $Ncell_GMratio = $parm{Ncell_GMratio} || 0; # default is full-tumor-cell
    my $ObjCopyRatioInNcell = $parm{ObjCopyRatioInNcell} || 0; # default is somatic

    my $ObjCountOfPureTumorCell = $ObjCount - $MixedCellMeanDepth * $Ncell_GMratio * $ObjCopyRatioInNcell;

    return $ObjCountOfPureTumorCell;
}

#--- get depth one haplotype in mixed situation ---
sub get_ObjSingleCNdepthInMixed{
    # options
    shift if (@_ && $_[0] =~ /::$MODULE_NAME/);
    my %parm = @_;
    my $MixedCellMeanDepth = $parm{MixedCellMeanDepth};
    my $TumorCellPurity = $parm{TumorCellPurity} || 1; # default is full-tumor-cell
    my $TumorCellPloidy = $parm{TumorCellPloidy} || 2; # default is diploid
    my $NormalCellPloidy = $parm{NormalCellPloidy} || 2; # default is diploid

    my $SingleCNdepth = $MixedCellMeanDepth / ( $TumorCellPurity * $TumorCellPloidy + ( 1 - $TumorCellPurity ) * $NormalCellPloidy );

    return $SingleCNdepth;
}

#--- 
1; ## tell the perl script the successful access of this module.
