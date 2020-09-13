package FuseSV::Virus_Integ::Enrichment::hostBKlocEnrich;

use strict;
use warnings;
use List::Util qw/ max min sum first /;
use Data::Dumper;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::FileHeader qw/ getHeaderTag lineInfoToHash /;
use BioFuse::BioInfo::BED qw/ read_bed_file /;
use BioFuse::Util::Index qw/ IndexRegion FindOverlapIdxRegion /;
use BioFuse::BioInfo::Position qw/ load_region_for_randPos get_random_pos /;
use FuseSV::Objects::vInteg_OB;
use BioFuse::Stat::ChiSquareTest::FourFoldTable;
use BioFuse::Visual::Objects::Axis;
use BioFuse::Visual::Objects::BiAxis;
use BioFuse::Visual::SVG_Util::SVGWork qw/ initialize_SVG_obj output_SVG_file /;
use BioFuse::Visual::Objects::Histogram;
use FuseSV::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              host_bkLoc_enrich
              test_locHost_para
              load_allSampleVitg
              load_oneSampleVitg
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Integ::Enrichment::hostBKlocEnrich';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        host_bkLoc_enrich
                        test_locHost_para
                        load_host_refseg
                        load_concern_region
                        load_allSampleVitg
                        load_oneSampleVitg
                        host_realBkLoc_analysis
                        judge_BK_loc
                        load_host_nonN
                        host_simuBkLoc_analysis
                        draw_enrich_spectrum
                     /;

#--- do enrichment test on host break-point distribution ---
sub host_bkLoc_enrich{
    # load host refseg list
    &load_host_refseg;
    # load region bed
    &load_concern_region;
    # load virus integration list of samples
    &load_allSampleVitg;
    # analysis on host-bk of real virus integration
    &host_realBkLoc_analysis;
    # load nonN region of host reference
    &load_host_nonN;
    # analysis on host-bk of simu virus integration
    &host_simuBkLoc_analysis;
    # draw enrich spectrum
    &draw_enrich_spectrum;
}

#--- test on required options ---
sub test_locHost_para{
    return  (   $V_Href->{enrich_mode} eq 'loc:host'
             && (   !file_exist(filePath=>$V_Href->{enrichRegion_bed})
                 || !file_exist(filePath=>$V_Href->{nonN_bed})
                )
            );
}

#--- load refseg list of host ---
sub load_host_refseg{
    # read host ref fai
    open (FAI, Try_GZ_Read($V_Href->{hostRef_fa}.'.fai')) || die "cannot read ref_fai: $!\n";
    while(<FAI>){
        my ($refseg, $reflen) = (split)[0,1];
        $V_Href->{hostRefseg}->{$refseg} = $reflen;
    }
    close FAI;
    # inform
    stout_and_sterr "[INFO]\tload host reference segment list OK.\n";
}

#--- load concern region from bed file ---
sub load_concern_region{
    # read bed file
    $V_Href->{ConcernRegion} = read_bed_file(bedFile => $V_Href->{enrichRegion_bed}, loadAsBED => 0);
    # reset to regions' index
    for my $regNameHf (values %{$V_Href->{ConcernRegion}}){
        for my $refseg (keys %$regNameHf){
            $regNameHf->{$refseg} = IndexRegion(ItvAref => $regNameHf->{$refseg}, winSize => $V_Href->{winSize});
        }
    }
    # inform
    stout_and_sterr "[INFO]\tload concern region bed file OK.\n";
}

#--- load virus integration list of samples ---
## record vitg to sample
sub load_allSampleVitg{
    open (SAMVI, Try_GZ_Read($V_Href->{sampleToVitg_list})) || die "fail to read sampleToVitg list: $!\n";
    while(<SAMVI>){
        next if /^#/;
        my ($sample, $vItg_list) = (split)[0,1];
        $V_Href->{SampleToVitg}->{$sample} = [];
        &load_oneSampleVitg(vItgList => $vItg_list, vItgAref => $V_Href->{SampleToVitg}->{$sample});
    }
    close SAMVI;
    # inform
    stout_and_sterr "[INFO]\tload virus integrations of all samples OK.\n";
}

#--- load virus integration of one sample ---
sub load_oneSampleVitg{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $vItgList = $parm{vItgList};
    my $vItgAref = $parm{vItgAref};

    my @objPara = qw/ refseg_p5 bkpos_p5 strand_p5 refseg_p3 bkpos_p3 strand_p3 inner_ins t_uniq_sr t_uniq_jr n_uniq_sr n_uniq_jr /;
    open (VITG, Try_GZ_Read($vItgList)) || die "fail to read vItg list: $!\n";
    my $headTagAf = getHeaderTag(header => ($_=<VITG>)); # get header
    while(<VITG>){
        next if /^\#/;
        s/\(\d+\)//; # alternative virus bkpos
        my $infoHf = lineInfoToHash(headTagAf => $headTagAf, lineInfo => $_);
        my @para = map{ $_ => $infoHf->{$_} } @objPara;
        my $vInteg_OB = FuseSV::Objects::vInteg_OB->new(@para);
        # create vItg object
        push @$vItgAref, $vInteg_OB;
    }
    close VITG;
}

#--- analysis host bk-pos reigon of real vItg ---
sub host_realBkLoc_analysis{
    my $realVitgReport = catfile($V_Href->{outdir}, 'real.vitg.region.report');
    open (RVIR, Try_GZ_Write($realVitgReport)) || die "fail to write real vItg report: $!\n";
    print RVIR '#'.join("\t", qw/ sample vitg_info region /)."\n";
    for my $sample (sort keys %{$V_Href->{SampleToVitg}}){
        $V_Href->{realVitgCount} += scalar @{$V_Href->{SampleToVitg}->{$sample}};
        for my $vInteg_OB (@{$V_Href->{SampleToVitg}->{$sample}}){
            my $hostBkAf = $vInteg_OB->get_host_bk(hostRefHref => $V_Href->{hostRefseg});
            my $hostLcHf = &judge_BK_loc(bkseg => $hostBkAf->[0], bkpos => $hostBkAf->[1]);
            # output
            my @output = ($sample);
            push @output, $vInteg_OB->get_vInteg_info;
            push @output, join(';', map {"$_:$hostLcHf->{$_}"} sort keys %$hostLcHf);
            print RVIR join("\t", @output)."\n";
            # record
            $V_Href->{vItgLocCount}->{real}->{$_}->{($hostLcHf->{$_}?1:0)}++ for keys %$hostLcHf;
        }
    }
    close RVIR;
    # inform
    stout_and_sterr "[INFO]\tanalyze host breakpoint of real virus integations OK.\n";
}

#--- judegment on location of given bkpos ---
sub judge_BK_loc{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $bkseg = $parm{bkseg};
    my $bkpos = $parm{bkpos};

    my $BkLocHf = {};
    # global distance?
    my $global_dist = $V_Href->{vItg_distance}->{_ALL_} || 0;
    for my $regName (sort keys %{$V_Href->{ConcernRegion}}){
        if(!exists $V_Href->{ConcernRegion}->{$regName}->{$bkseg}){
            $BkLocHf->{$regName} = 0;
        }
        else{
            # specific distance
            my $reg_dist = $V_Href->{vItg_distance}->{$regName} || 0;
            my $distance = $global_dist + $reg_dist;
            my $ovpRegAf = FindOverlapIdxRegion( IdxItvHref => $V_Href->{ConcernRegion}->{$regName}->{$bkseg},
                                                 regionAref => [$bkpos-$distance, $bkpos+$distance],
                                                 winSize => $V_Href->{winSize},
                                                 queryMode => 'itv'
                                               );
            $BkLocHf->{$regName} = scalar @$ovpRegAf;
        }
    }
    return $BkLocHf;
}

#--- load host nonN region ---
sub load_host_nonN{
    load_region_for_randPos( ItvFile => $V_Href->{nonN_bed},
                             ItvAref => $V_Href->{hostNonNRegion}
                           );
    # inform
    stout_and_sterr "[INFO]\tload host nonN region OK.\n";
}

#--- analysis host bk-pos reigon of simu vItg ---
sub host_simuBkLoc_analysis{
    # simu vItg reports
    my $simuVitgReport = catfile($V_Href->{outdir}, 'simu.vitg.region.report');
    open (SVIR, Try_GZ_Write($simuVitgReport)) || die "fail to write simu vItg report: $!\n";
    print SVIR '#'.join("\t", qw/ simu_R simu_C vitg_info region /)."\n";
    # enrichment test report
    my $enrichmentReport = catfile($V_Href->{outdir}, 'enrichment.report');
    open (ERCR, Try_GZ_Write($enrichmentReport)) || die "fail to write vItg enrichment report: $!\n";
    print ERCR '#'.join("\t", qw/ simu_R region real_Y real_N simu_Y simu_N real_r simu_r method ChiSquare Pvalue oddLog /)."\n";
    # simulation
    $V_Href->{simulation_count} = $V_Href->{realVitgCount} unless defined $V_Href->{simulation_count};
    for my $simu_R (1 .. $V_Href->{simulation_round}){
        # get randoom bk-pos
        my $simuBKAf = get_random_pos( ItvAref => $V_Href->{hostNonNRegion},
                                       randCnt => $V_Href->{simulation_count}
                                     );
        # analysis location
        $V_Href->{vItgLocCount}->{simu} = {}; # reset
        for my $i (0 .. scalar(@$simuBKAf)-1){
            my $hostLcHf = &judge_BK_loc(bkseg => $simuBKAf->[$i]->[0], bkpos => $simuBKAf->[$i]->[1]);
            # output
            my @output = ($simu_R, $i+1);
            push @output, join(';', @{$simuBKAf->[$i]});
            push @output, join(';', map {"$_:$hostLcHf->{$_}"} sort keys %$hostLcHf);
            print SVIR join("\t", @output)."\n";
            # record
            $V_Href->{vItgLocCount}->{simu}->{$_}->{($hostLcHf->{$_}?1:0)}++ for keys %$hostLcHf;
        }
        # real-vs-simu compare
        for my $regName (sort keys %{$V_Href->{vItgLocCount}->{simu}}){
            my $real_Y = $V_Href->{vItgLocCount}->{real}->{$regName}->{1} || 0;
            my $real_N = $V_Href->{vItgLocCount}->{real}->{$regName}->{0} || 0;
            my $simu_Y = $V_Href->{vItgLocCount}->{simu}->{$regName}->{1} || 0;
            my $simu_N = $V_Href->{vItgLocCount}->{simu}->{$regName}->{0} || 0;
            my $ffTable = BioFuse::Stat::ChiSquareTest::FourFoldTable->new(Value_Href => {r1c1=>$real_Y, r1c2=>$real_N, r2c1=>$simu_Y, r2c2=>$simu_N});
            my $method = $ffTable->get_method;
            my $Pvalue = $ffTable->get_P_value;
            my $ChiSquare = $ffTable->get_chi_square;
            my $real_r = sprintf "%0.3f",($real_Y/($real_Y+$real_N));
            my $simu_r = sprintf "%0.3f",($simu_Y/($simu_Y+$simu_N));
            my $ORhref = $ffTable->get_odds_ratio;
            my $oddLog = sprintf "%0.3f", log($ORhref->{OR});
            print ERCR join("\t", $simu_R, $regName, $real_Y, $real_N, $simu_Y, $simu_N, $real_r, $simu_r, $method, $ChiSquare, $Pvalue, $oddLog)."\n";
            # record rate
            $real_r = int($real_r / $V_Href->{positiveRateUnit} + 1) * $V_Href->{positiveRateUnit};
            $V_Href->{positiveRate}->{$regName}->{real} = $real_r;
            $simu_r = int($simu_r / $V_Href->{positiveRateUnit} + 1) * $V_Href->{positiveRateUnit};
            $V_Href->{positiveRate}->{$regName}->{simu}->{$simu_r} ++;
        }
        # inform
        stout_and_sterr "[INFO]\tfinish $simu_R round simulations of virus integration and analyze OK.\n" if $simu_R % 1000 == 0;
    }
    close SVIR;
    close ERCR;
    # inform
    stout_and_sterr "[INFO]\tTotally, $V_Href->{simulation_round} round simulations of virus integration and analyze OK.\n";
}

#--- draw enrichment histogram figure ---
sub draw_enrich_spectrum{
    for my $regName (sort keys %{$V_Href->{positiveRate}}){
        my $spec_svg_file = catfile($V_Href->{outdir}, "$regName.enrich.histogram.svg");
        my $simu_r2C_Hf = $V_Href->{positiveRate}->{$regName}->{simu};
        # histogram in bi-axis
        my $x_axis = BioFuse::Visual::Objects::Axis->new(origP_X=>100, origP_Y=>300, axisLen=>300, headAng=>90);
        my $y_axis = BioFuse::Visual::Objects::Axis->new(origP_X=>100, origP_Y=>300, axisLen=>200, headAng=>0);
        my $biAxis = BioFuse::Visual::Objects::BiAxis->new(axis_1=>$x_axis, axis_2=>$y_axis);
        my $histogram = BioFuse::Visual::Objects::Histogram->new(bi_axis => $biAxis, pillarWid => 15);
        # axis setting
        $x_axis->add_resol(minValue => 0, maxValue => 1);
        $y_axis->add_resol(minValue => 0, maxValue => 100);
        $x_axis->set_tic(len => 3, width => 1, clockwise => 1);
        $y_axis->set_tic(len => 3, width => 1, clockwise => 0);
        $x_axis->set_stub(stubStep => 0.1);
        $y_axis->set_stub(stubStep => 10, verticalToBL => 1, stubText => "($V_Href->{simulation_round},%)");
        # histogram data
        for my $rate (keys %$simu_r2C_Hf){
            my $perc = sprintf "%.2f", 100 * $simu_r2C_Hf->{$rate} / $V_Href->{simulation_round};
            $histogram->load_data(vA => $rate, vB=>[$perc]);
        }
        # draw
        my $svg = initialize_SVG_obj;
        $histogram->draw(svg_obj => $svg);
        output_SVG_file(svg_obj => $svg, svg_file => $spec_svg_file);
        # inform
        stout_and_sterr "[INFO]\toutput $regName enrichment histogram figure OK.\n";
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
