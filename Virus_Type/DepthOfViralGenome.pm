package FuseSV::Virus_Type::DepthOfViralGenome;

use strict;
use warnings;
use File::Basename qw/basename dirname/;
use File::Spec qw/catfile/;
use Cwd qw/abs_path/;
use List::Util qw[min max sum];
use Getopt::Long;
use SOAPfuse::General_Operation qw/warn_and_exit trible_run_for_success stout_and_sterr/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use FuseSV::LoadOn;
use FuseSV::Virus_Type::VirusTypeMain;
use FuseSV::GetPath qw/ GetPath /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              viral_genome_depth_stat
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'DepthOfViralGenome';
#----- version --------
$VERSION = "0.03";
$DATE = '2017-06-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        viral_genome_depth_stat
                        regmap_to_orig_pos
                        find_virus_vNO_for_depth_stat
                        prepare_temp_space_for_depth_work
                        get_depth_of_viral_region
                        filter_softclip_length
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
     Usage:   perl $V_Href->{MainName} virus_depth <[Options]>
     
     Options:

       # Inputs and Outputs #
        -id     [s]  symbol name that used as '-id' in step 'virus_type'. <required>
                      only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
        -vdir   [s]  directory that used as '-odir' in step 'virus_type'. <required>
        -o      [s]  file to store pos depth list. <required>
                      Could output gzip/bgzip/plain_text according to the postfix of file-name.

       # Software Required #
        -stl    [s]  SamTools, minimum version: 1.3. <required>
        -bgz    [s]  bgzip tool. <required>
        -tbi    [s]  tabix tool. <required>

       # Depth Options #
        -r      [s]  region to stat position depth. <required>
                      Format: 'vseg' or 'vseg:from-to'
                      Note that the 'vseg' is the column '\033[4m \bVariantVirus \b\033[0m ' in virus_variant report.
        -m      [i]  maximum coverage depth. [100000]
        -bqual  [i]  base quality threshold. [0]
        -mqual  [i]  mapping quality threshold. [0]
        -minml  [i]  allowed minimum match part length for one read. [30]
        -sfcpr  [f]  allowed maximum softclip part ratio for one read. [0.7]
        -orig        output depth of original viral position. [Disabled]
        -mean        use this option to get mean depth of give region. [Disabled]

        -h|help   Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
       # # Other Options #
       #  -daddag [s]  Specially deal with additional alignments (0x100/0x800). [Disabled]
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    $V_Href->{ $_->[0] } = $_->[1] for
        (
            # input/output
            [ sample_ID => '' ],
            [ outdir => '' ],
            [ viral_depth_output => '' ],
            # software
            [ samtools => '' ],
            [ bgzip => '' ],
            [ tabix => '' ],
            # intermediate file
            [ Virus_variant_report => '' ],
            # depth
            [ viral_depth_region => '' ],
            [ viral_depth_vseg => '' ],
            [ viral_depth_vNO => 0 ],
            [ viral_depth_maxdepth => 100000 ],
            [ viral_depth_minbqual => 0 ],
            [ viral_depth_minmqual => 0 ],
            [ viral_depth_MinAllowMatchLen => 30 ],
            [ viral_depth_MaxAllowSoftClipPart => 0.7 ],
            [ viral_depth_getorigpos => 0 ],
            [ viral_depth_getmean => 0 ],
            [ viral_depth_tempdir => '' ],
            [ viral_depth_dealadd => 0 ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['outdir'],
                                  ['viral_depth_output'],
                                  ['samtools'],
                                  ['bgzip'],
                                  ['tabix']               ] ]
        );

    # load containers
    FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool( only_load_containers => 1 );
}

#--- get options from command line
sub Get_Cmd_Options{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # get options
    GetOptions(
        # input/output
        "-id:s"     => \$V_Href->{sample_ID},
        "-vdir:s"   => \$V_Href->{outdir},
        "-o:s"      => \$V_Href->{viral_depth_output},
        # software
        "-stl:s"    => \$V_Href->{samtools},
        "-bgz:s"    => \$V_Href->{bgzip},
        "-tbi:s"    => \$V_Href->{tabix},
        # depth
        "-r:s"      => \$V_Href->{viral_depth_region},
        "-m:i"      => \$V_Href->{viral_depth_maxdepth},
        "-bqual:i"  => \$V_Href->{viral_depth_minbqual},
        "-mqual:i"  => \$V_Href->{viral_depth_minmqual},
        "-minml:i"  => \$V_Href->{viral_depth_MinAllowMatchLen},
        "-sfcpr:f"  => \$V_Href->{viral_depth_MaxAllowSoftClipPart},
        "-orig"     => \$V_Href->{viral_depth_getorigpos},
        "-mean"     => \$V_Href->{viral_depth_getmean},
        "-daddag"   => \$V_Href->{viral_depth_dealadd}, # hidden option
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
             || (!$V_Href->{outdir} || !-d $V_Href->{outdir})
             || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
             || (!$V_Href->{bgzip} || !-e $V_Href->{bgzip})
             || (!$V_Href->{tabix} || !-e $V_Href->{tabix})
             || !$V_Href->{viral_depth_region}
             || !$V_Href->{viral_depth_output}
            );
}

#--- depth stat of viral genome ---
sub viral_genome_depth_stat{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # file path
    $V_Href->{Virus_variant_folder} = GetPath(filekey => 'Virus_variant_folder');
    $V_Href->{Virus_variant_report} = GetPath(filekey => 'Virus_variant_report');

    # read virus variant report
    FuseSV::Virus_Type::ReconstrcutVirusVariants->loadon_virus_variant_report;

    # find the virus seg for depth stat
    &find_virus_vNO_for_depth_stat;

    # prepare temp workspace
    &prepare_temp_space_for_depth_work;

    # get depth of required viral region
    &get_depth_of_viral_region;

    # sweep temp workspace
    `rm -rf $V_Href->{viral_depth_tempdir}`;
}

#--- according to the regmap file, map ModPos to OrigPos, update the temp depbgz ---
sub regmap_to_orig_pos{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $regmap = $Option_Href->{regmap};
    my $depbgz = $Option_Href->{depbgz};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # read mapping region
    my @mapreg;
    my $shift_bool = 0;
    open (MR, Try_GZ_Read($regmap)) || die "fail read regmap: $!\n";
    while(<MR>){
        next if(/^#/);
        my ($vRefName, $newRegSt, $newRegEd, $origRegSt, $origRegEd) = (split)[0,1,2,4,5];
        # check
        if( $vRefName ne $V_Href->{viral_depth_vseg} ){
            warn_and_exit `date`."<ERROR>:\tWrong virus variant ref-name read from regmap file.\n"
                                ."\tName: $vRefName ne $V_Href->{viral_depth_vseg}\n"
                                ."\tFile: $regmap\n";
        }
        # record
        push @mapreg, {newSt=>$newRegSt, newEd=>$newRegEd, origSt=>$origRegSt, origEd=>$origRegEd};
        # need to translate?
        $shift_bool ||= ( $newRegSt != $origRegSt || $newRegEd != $origRegEd );
    }
    close MR;

    unless( $shift_bool ){
        # inform
        stout_and_sterr `date`."[INFO]:\tNo need to map ModPos to OrigPos as no shift mapping viral region.\n";
        return;
    }

    # read depbgz and output
    open (DEPIN, "$V_Href->{tabix} $depbgz $V_Href->{viral_depth_region} | ") || die "fail read depbgz: $!\n";
    open (DEPOP, "| $V_Href->{bgzip} -c > $depbgz.origpos.bgz") || die "fail write depbgz.origpos.bgz: $!\n";
    my $RegIndex = 0;
    my $last_OrigPos = 1;
    while(<DEPIN>){
        my ($vRefName, $newPos, $depth) = (split)[0,1,2];
        # change to other index?
        while(    $newPos > $mapreg[$RegIndex]->{newEd} # out of range
               && $RegIndex < $#mapreg # not overflow
        ){
            $RegIndex ++;
        }
        # check end pos, error case
        if( $newPos > $mapreg[$RegIndex]->{newEd} ){
            warn_and_exit `date`."<ERROR>:\tMeet ModPos ($newPos) out of mapping region.\n";
        }
        # insertion
        if( $newPos < $mapreg[$RegIndex]->{newSt} ){
            next;
        }
        # find the mapping original pos
        my $origPos = $mapreg[$RegIndex]->{origSt} + ( $newPos - $mapreg[$RegIndex]->{newSt} );
        # output
        ## any gap? should be deletion
        if( $origPos > $last_OrigPos+1 ){
            print DEPOP join("\t", $vRefName, $_, 0) . "\n" for ($last_OrigPos+1 .. $origPos-1);
        }
        ## this pos
        print DEPOP join("\t", $vRefName, $origPos, $depth) . "\n";
        ## update
        $last_OrigPos = $origPos;
    }
    close DEPOP;
    close DEPIN;

    # update depbgz
    `rm -f $depbgz $depbgz.tbi`;
    `mv -f $depbgz.origpos.bgz $depbgz`;
    `$V_Href->{tabix} -s 1 -b 2 -e 2 $depbgz`;

    # inform
    stout_and_sterr `date`."[INFO]:\tMap ModPos to OrigPos OK.\n";
}

#--- find the virus vNO that will calculate depth ---
sub find_virus_vNO_for_depth_stat{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

    # match the virus variant segment name
    for my $variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        if( $Variant_Info_Href->{$variant_NO}->{RefSeg} eq $V_Href->{viral_depth_vseg} ){
            $V_Href->{viral_depth_vNO} = $variant_NO;
            last;
        }
    }

    # check and warn
    if( $V_Href->{viral_depth_vNO} == 0 ){
        warn_and_exit `date`."<ERROR>:\tCannot find the virus variant match your input.\n"
                            ."\t$V_Href->{viral_depth_vseg}.\n";
    }
}

#--- prepare temp workspace ---
sub prepare_temp_space_for_depth_work{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $TdirName = "$V_Href->{sample_ID}.DepthStat.temp_dir";

    if( $V_Href->{viral_depth_output} ){
        my $output_dir = abs_path( dirname($V_Href->{viral_depth_output}) );
        if( !-w $output_dir ){
            warn_and_exit `date`."<ERROR>:\tCannot create temp workspace. E1\n";
        }
        else{
            $V_Href->{viral_depth_tempdir} = File::Spec->catfile($output_dir, $TdirName);
        }
    }
    # elsif( -w $V_Href->{Virus_variant_folder} ){
    #   $V_Href->{viral_depth_tempdir} = File::Spec->catfile($V_Href->{Virus_variant_folder}, $TdirName);
    # }
    else{
        warn_and_exit `date`."<ERROR>:\tCannot create temp workspace. E2\n";
    }

    # create folder
    `rm -rf $V_Href->{viral_depth_tempdir}`;
    `mkdir -p $V_Href->{viral_depth_tempdir}`;
}

#--- get depth of required viral region ---
sub get_depth_of_viral_region{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $viral_ext_len = $Variant_Info_Href->{$V_Href->{viral_depth_vNO}}->{Extend_length};

    my $bam_for_stat = $Variant_Info_Href->{$V_Href->{viral_depth_vNO}}->{finalBam};

    #------------------------------------#
    # operation on additional alignments #
    #------------------------------------#
    if( $V_Href->{viral_depth_dealadd} ){
        my $M900_mkdup_bam = File::Spec->catfile($V_Href->{viral_depth_tempdir}, "$V_Href->{viral_depth_vseg}.for_depth.M900_mkdup.bam");
        FuseSV::Virus_Type::SmartAlignment->mark_additional_alignments(
            {
                source_bam => $bam_for_stat,
                output_bam => $M900_mkdup_bam
            }
        );
        # index
        `$V_Href->{samtools} index $M900_mkdup_bam`;
        # update
        $bam_for_stat = $M900_mkdup_bam;
        # inform
        stout_and_sterr `date`."[INFO]:\tMark duplication of 0x100/0x800 reads OK.\n";
    }

    #----------------------------------#
    # operation on softclip alignments #
    #----------------------------------#
    my $filter_SoftClip_bam = File::Spec->catfile($V_Href->{viral_depth_tempdir}, "$V_Href->{viral_depth_vseg}.for_depth.filter_SoftClip.bam");
    &filter_softclip_length(
        {
            source_bam => $bam_for_stat,
            output_bam => $filter_SoftClip_bam,
            MaxAllowSoftClipPart => $V_Href->{viral_depth_MaxAllowSoftClipPart},
            MinAllowMatchLen => $V_Href->{viral_depth_MinAllowMatchLen}
        }
    );
    # index
    `$V_Href->{samtools} index $filter_SoftClip_bam`;
    # sweep the internal M900_mkdup_bam
    if( $V_Href->{viral_depth_dealadd} ){
        `rm -rf $bam_for_stat $bam_for_stat.bai`;
    }
    # update
    $bam_for_stat = $filter_SoftClip_bam;
    # inform
    stout_and_sterr `date`."[INFO]:\tFilter softclip part and match length OK.\n";

    #--------------------------------#
    # samtools depth to get raw info #
    #--------------------------------#
    my $depth_PosList_bgz = File::Spec->catfile($V_Href->{viral_depth_tempdir}, "$V_Href->{viral_depth_vseg}.all_pos_depth.raw.bgz");
    # my $r_para = ($viral_ext_len == 0) ? $V_Href->{viral_depth_region} : $V_Href->{viral_depth_vseg};
    my $depth_cmd =  "$V_Href->{samtools} depth -a -d $V_Href->{viral_depth_maxdepth} -q $V_Href->{viral_depth_minbqual} -Q $V_Href->{viral_depth_minmqual} -r $V_Href->{viral_depth_vseg} $bam_for_stat "
                    ."| $V_Href->{bgzip} -c > $depth_PosList_bgz";
    trible_run_for_success($depth_cmd, 'Depth', {esdo_Nvb=>1});
    # inform
    stout_and_sterr `date`."[INFO]:\tSamTools depth to get pos_depth raw info OK.\n";

    #--------------------------------#
    # merge depth for circular virus #
    #--------------------------------#
    my $depth_PosMerged_bgz = File::Spec->catfile($V_Href->{viral_depth_tempdir}, "$V_Href->{viral_depth_vseg}.all_pos_depth.merged.bgz");
    if( $viral_ext_len != 0 ){
        chomp(my $viral_genome_size = `$V_Href->{samtools} view -H $bam_for_stat | grep '^\@SQ' | grep '$V_Href->{viral_depth_vseg}' | xargs -n1 | grep '^LN' | sed 's/LN://'`);
        my $viral_noext_len = $viral_genome_size - $viral_ext_len;
        my $merge_cmd =  "( gzip -cd $depth_PosList_bgz | awk '\$2>$viral_noext_len' | cut -f3 > $depth_PosMerged_bgz.ext_part ) && "
                        ."( gzip -cd $depth_PosList_bgz | awk '\$2<=$viral_noext_len' "
                        ."| paste - $depth_PosMerged_bgz.ext_part "
                        .'| perl -ne \'my @ele=split; my $sum_depth=0; $sum_depth+=$_ for @ele[2..$#ele]; print join("\t", @ele[0,1], $sum_depth)."\n";\' '
                        ."| $V_Href->{bgzip} -c > $depth_PosMerged_bgz ) && "
                        ."( rm -rf $depth_PosMerged_bgz.ext_part )";
        trible_run_for_success($merge_cmd, 'MegerDepth', {esdo_Nvb=>1});
        # inform
        stout_and_sterr `date`."[INFO]:\tMerge extended part to get pos_depth merged info OK.\n";
    }
    else{
        `gzip -cd $depth_PosList_bgz | $V_Href->{bgzip} -c > $depth_PosMerged_bgz`;
    }
    # tabix index
    `$V_Href->{tabix} -s 1 -b 2 -e 2 $depth_PosMerged_bgz`;
    # inform
    stout_and_sterr `date`."[INFO]:\tIndex pos_depth list OK.\n";

    #---------------------------------#
    # translate to original viral pos #
    #---------------------------------#
    if(    $V_Href->{viral_depth_region} eq $V_Href->{viral_depth_vseg} # just whole viral segment
        || $V_Href->{viral_depth_getorigpos} # required to convert
    ){
        # check MutCompletedBool
        if( $Variant_Info_Href->{$V_Href->{viral_depth_vNO}}->{ModCompBool} ){
            my $variant_outdir = dirname( $Variant_Info_Href->{$V_Href->{viral_depth_vNO}}->{finalBam} );
            my $variant_regmap = File::Spec->catfile($variant_outdir, "V$V_Href->{viral_depth_vNO}.regmap.list");
            # update depth_PosMerged_bgz
            &regmap_to_orig_pos(
                {
                    regmap => $variant_regmap,
                    depbgz => $depth_PosMerged_bgz
                }
            );
        }
        else{
            # inform
            stout_and_sterr `date`."[WARN]:\tCannot output the depth info of region from original vrial genome.\n"
                                         ."\tAs selected virus variant's mutations is not completed mapped.\n";
        }
    }
    else{
        # inform
        stout_and_sterr `date`."[INFO]:\tOutput the depth info of region from reconstrcut vrial genome.\n";
    }

    #---------------------------------#
    # report results based on options #
    #---------------------------------#
    my $fifo_cmd;
    if   ( $V_Href->{viral_depth_output} =~ /\.gz$/ ) { $fifo_cmd = '| gzip -c ' }
    elsif( $V_Href->{viral_depth_output} =~ /\.bgz$/ ){ $fifo_cmd = "| $V_Href->{bgzip} -c " }
    else { $fifo_cmd = '' }
    # extract depth list
    `$V_Href->{tabix} $depth_PosMerged_bgz $V_Href->{viral_depth_region} $fifo_cmd > $V_Href->{viral_depth_output}`;
    # inform
    stout_and_sterr `date`."[INFO]:\tGenerate pos_depth list file OK.\n"
                          ."\t$V_Href->{viral_depth_output}\n";
    # get mean depth?
    if( $V_Href->{viral_depth_getmean} ){
        chomp(my $pos_count = `$V_Href->{tabix} $depth_PosMerged_bgz $V_Href->{viral_depth_region} | wc -l`);
        chomp(my $depth_sum = `$V_Href->{tabix} $depth_PosMerged_bgz $V_Href->{viral_depth_region} | awk 'BEGIN{sum=0} {sum+=\$3} END{print sum}'`);
        my $mean_depth = sprintf ("%.2f", $depth_sum / $pos_count);
        # inform
        stout_and_sterr `date`."[INFO]:\tMean depth of region ($V_Href->{viral_depth_region}) is $mean_depth\n";
    }
}

#--- filter softclip ratio and match length ---
## just mark them as duplication.
sub filter_softclip_length{
    # options
    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $source_bam = $Option_Href->{source_bam};
    my $output_bam = $Option_Href->{output_bam};
    my $MaxAllowSoftClipPart = $Option_Href->{MaxAllowSoftClipPart};
    my $MinAllowMatchLen = $Option_Href->{MinAllowMatchLen};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    open (my $obam_fh,"| $V_Href->{samtools} view -b -S -o $output_bam")|| die"fail samtools write $output_bam: $!\n";
    open (SAM,"$V_Href->{samtools} view -h $source_bam |") || die"fail samtools view header of $source_bam: $!\n";
    while(<SAM>){
        # SAM header
        if(/^@/){
            print {$obam_fh} $_;
            next;
        }
        # Read alignment info
        my @ele = split;
        my ($cigar, $rseq) = @ele[5,9];
        # deal with cigar
        my $M_sum = 0;
        my $S_sum = 0;
        while($cigar =~ /^(\d+)(\D)/){
            if( $2 eq 'M' ){
                $M_sum += $1;
            }
            elsif( $2 eq 'S' ){
                $S_sum += $1;
            }
            $cigar =~ s/^\d+\D//;
        }
        # judge to output
        if(    $M_sum < $MinAllowMatchLen
            || $S_sum / length($rseq) > $MaxAllowSoftClipPart
        ){
            # mark duplication (0x400)
            $ele[1] += 0x400 unless( $ele[1] & 0x400 );
            print $obam_fh join("\t", @ele) . "\n";
        }
        else{
            print $obam_fh $_;
        }
    }
    close SAM;
    close $obam_fh;
}

#--- 
1; ## tell the perl script the successful access of this module.
