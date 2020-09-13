package FuseSV::RunFunc;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use JSON;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read /;
use FuseSV::LoadOn;
use FuseSV::Check;
use FuseSV::GetPath;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              options_alert_and_run
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::RunFunc';
#----- version --------
$VERSION = "0.49";
$DATE = '2018-07-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        options_alert_and_run
                        func_run
                        load_functions
                        return_HELP_INFO
                        HelpInfo_ExtractOpt
                        para_alert
                     /;

#--- get options ---
sub options_alert_and_run{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $argv_Aref = $parm{argv_Aref};

    if(scalar(@$argv_Aref) == 0){
        $V_Href->{command} = 'none';
    }
    else{
        $V_Href->{command} = shift @$argv_Aref;
        $V_Href->{argv_Aref} = $argv_Aref;
    }

    &load_functions;
    &HelpInfo_ExtractOpt;
    &para_alert;
    &func_run;
}

#--- run given func ---
sub func_run{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $func_name = $parm{func_name} || $V_Href->{command};

    if( !exists $V_Href->{func}->{$func_name} ){
        warn_and_exit "<ERROR>\tcannot recognize func $func_name.\n";
    }

    # pre-run functions
    if(    exists $V_Href->{func}->{$func_name}->{pre_run}
        && $V_Href->{run_mode} == $V_Href->{func}->{$func_name}->{mode_NO}
    ){
        my $preRun_func = $V_Href->{func}->{$func_name}->{pre_run};
        &func_run( func_name => $_ ) for @$preRun_func;
    }

    # main function
    my $main_func = $V_Href->{func}->{$func_name}->{main_func};
    if( ref($main_func) && ref($main_func) eq 'ARRAY' ){ # workflow
        &func_run( func_name => $_ ) for @$main_func;
    }
    else{
        my $module = $V_Href->{func}->{$func_name}->{module};
        "$module"->$main_func;
    }
}

#--- load function information from json config ---
sub load_functions{
    my $func_json = FuseSV::GetPath->GetPath( filekey => 'func_json' );
    # read json
    my $json_txt;
    open (FCJSON, Try_GZ_Read($func_json)) || die "cannot read functions.json";
    $json_txt .= $_ while(<FCJSON>);
    close FCJSON;
    # load func
    $V_Href->{func} = from_json( $json_txt, { utf8  => 1 } );
    # validate
    my %prev;
    for my $InfoHref (values %{$V_Href->{func}}){
        for my $key (qw/ mode_NO module main_func /){
            my $value = $InfoHref->{$key};
            if( exists $prev{$key}{$value} ){
                warn_and_exit "<ERROR>\t$key($value) is duplicated in function.json.\n";
            }
            $prev{$key}{$value} = 1;
        }
    }
}

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} <command>
     
     Commands:

     -- Local Genomic Map
        draw_segCN     $V_Href->{func}->{draw_segCN}->{descp}
        get_UCYC       $V_Href->{func}->{get_UCYC}->{descp}
        UCYCtoLM       $V_Href->{func}->{UCYCtoLM}->{descp}

     Version:
        V$V_Href->{Version} at $V_Href->{Date}

     Author:
        $V_Href->{AUTHOR} ($V_Href->{EMAIL})
 \n";
        # LMcircos       $V_Href->{func}->{LMcircos}->{descp}
        # find_MH        $V_Href->{func}->{find_MH}->{descp}
}

#--- extract options for command ---
sub HelpInfo_ExtractOpt{

    my $second_cmd  = $V_Href->{command};

    if($second_cmd eq 'none'){
        # help info
        $V_Href->{HELP_INFO}->{$second_cmd} = &return_HELP_INFO;
    }
    elsif( exists $V_Href->{func}->{$second_cmd} ){
        # module
        my $module = $V_Href->{func}->{$second_cmd}->{module};
        # set run_mode
        $V_Href->{run_mode} = $V_Href->{func}->{$second_cmd}->{mode_NO};
        # load module variants
        "$module"->Load_moduleVar_to_pubVarPool;
        # help info
        $V_Href->{HELP_INFO}->{$second_cmd} = "$module"->return_HELP_INFO;
        # get options
        "$module"->Get_Cmd_Options;
    }
    else{
        # help info
        $V_Href->{HELP_INFO}->{none} = &return_HELP_INFO;
    }
}

#--- test para and alert
sub para_alert{

    my $second_cmd = $V_Href->{command};

    # alert help if possible
    if( !exists($V_Href->{HELP_INFO}->{$second_cmd}) ){
        if(    length($second_cmd) != 0
            && $second_cmd !~ /^\-{1,2}help$/
            && $second_cmd !~ /^\-{1,2}h$/
          ){
            warn  "\n <ERROR>: Cannot recognize command '$second_cmd'.\n"
                 ." Check Usage:\n";
        }
        warn_and_exit $V_Href->{HELP_INFO}->{none};
    }
    elsif($second_cmd eq 'none'){
        warn_and_exit $V_Href->{HELP_INFO}->{$second_cmd};
    }

    # abs_path, firstly
    for my $abs_aim_Aref ( @{ $V_Href->{ToAbsPath_Aref} } ){
        my $temp_Href = $V_Href;
        # hierarchical keys
        $temp_Href = $temp_Href->{$abs_aim_Aref->[$_]} for ( 0 .. scalar(@$abs_aim_Aref)-2 );
        if(    defined $temp_Href->{$abs_aim_Aref->[-1]}
            && $temp_Href->{$abs_aim_Aref->[-1]} ne ''
            && $temp_Href->{$abs_aim_Aref->[-1]} ne 'java' # avoids OS java
          ){
            $temp_Href->{$abs_aim_Aref->[-1]} = abs_path( $temp_Href->{$abs_aim_Aref->[-1]} );
        }
    }

    # then, module build-in para_check
    my $module = $V_Href->{func}->{$second_cmd}->{module};
    if( "$module"->para_alert ){
        warn_and_exit $V_Href->{HELP_INFO}->{$second_cmd};
    }

    # then, more specific checks
    FuseSV::Check->check;
}

#--- 
1; ## tell the perl script the successful access of this module.
