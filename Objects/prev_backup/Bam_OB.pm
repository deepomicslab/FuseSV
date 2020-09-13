package FuseSV::Objects::Bam_OB;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use BioFuse::Util::Sys qw/ trible_run_for_success file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Objects::Reads_OB;
use FuseSV::Objects::ReadsGroup_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Bam_OB';
#----- version --------
$VERSION = "0.08";
$DATE = '2018-05-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        verify_bam
                        verify_index
                        get_filepath
                        load_reads_for_ReadsGroup
                        rg_count_need_reads_ForIns
                        extract_ReadsGroup_OB
                        add_ReadsGroup_OBs
                        get_region_depth
                        delete_regionDepthFile
                        get_region_alt_vcf_gz
                        get_pos_marker_stat
                        get_allele_marker_stat
                     /;

#--- structure of object
# bam_OB -> tissue = $tissue
# bam_OB -> filepath = $filepath
# bam_OB -> rgOB = { rgID->$rgOB }, check FuseSV::Objects::ReadsGroup_OB
# bam_OB -> regionDepthFile = $regionDepthFile

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $bam_OB = {};
    $bam_OB->{tissue}   = $parm{tissue};
    $bam_OB->{filepath} = $parm{filepath};
    $bam_OB->{rgOB} = {};

    if( defined $bam_OB->{filepath} ){
        $bam_OB->{filepath} = abs_path( $bam_OB->{filepath} );
    }

    bless($bam_OB);
    return $bam_OB;
}

#--- verify bam file ---
sub verify_bam{

    my $bam_OB = shift;

    $bam_OB->{filepath} = abs_path( $bam_OB->{filepath} );

    if( ! file_exist(filePath=>$bam_OB->{filepath}) ){
        warn_and_exit "<ERROR>\tCannot find bam:\n"
                            ."\t$bam_OB->{filepath}\n";
    }
}

#--- verify index (bai file) ---
sub verify_index{

    my $bam_OB = shift;

    my $bai_filepath1 = $bam_OB->{filepath} . '.bai';
    (my $bai_filepath2 = $bam_OB->{filepath}) =~ s/bam$/bai/;
    if(    ! file_exist(filePath=>$bai_filepath1)
        && ! file_exist(filePath=>$bai_filepath2)
    ){
        warn_and_exit "<ERROR>\tCannot find index bai file of bam:\n"
                            ."\t$bam_OB->{filepath}\n";
    }
}

#--- return file path ---
sub get_filepath{
    my $bam_OB = shift;
    return $bam_OB->{filepath};
}

#--- load reads and INS-stat for all ReadsGroup ---
sub load_reads_for_ReadsGroup{

    my $bam_OB = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};
    my $only_SoftClip  = $parm{only_SoftClip};
    my $pr_DropProb    = $parm{pr_DropProb} || 0;
    my $pr_AimCount    = $parm{pr_AimCount} || 10000;
    my $Tool_Tag       = $parm{Tool_Tag} || '__NULL__';
    my $samtools       = $parm{samtools};

    # still has reads group to load reads?
    my $rgNeedReads = $bam_OB->rg_count_need_reads_ForIns;
    return if( $rgNeedReads == 0 );

    # to store the pid
    my $tissue = $bam_OB->{tissue};
    # FLAG: -f 0x2 (P), -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -f 0x2 -F 0xD00 $bam_OB->{filepath} |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 0, _rc_rdstr => 0 );
        # check RG_ID
        my $rg_OB = $reads_OB->get_available_rgOB( bam_OB => $bam_OB, rgid2rgOB_Href => $rgid2rgOB_Href );
        ## enough reads for ins evalue
        next if( $rg_OB->{EnoughReadsBool} );
        # extract read-id FS prefix to rg_OB, and prepare for update
        $reads_OB->extract_FS_readid_prefix( rg_OB => $rg_OB, Tool_Tag => $Tool_Tag ) unless( defined $rg_OB->{rID_prefix} );
        # count in this paired reads for RG ins evalue
        $rg_OB->load_reads_for_ins_evalue( reads_OB => $reads_OB, only_SoftClip => $only_SoftClip, pr_AimCount => $pr_AimCount, pr_DropProb => $pr_DropProb );
        # stop when all RG(s) are full
        $rgNeedReads-- if( $rg_OB->{EnoughReadsBool} );
        last if( $rgNeedReads == 0 );
    }
    close BAM;

    # inform
    stout_and_sterr "[INFO]\tLoad paired-end info for reads group (only_SoftClip:$only_SoftClip) from bam file.\n"
                         ."\t$bam_OB->{filepath}\n";
}

#--- whether all RG(s) of this bam get enough reads for ins evaluation ---
sub rg_count_need_reads_ForIns{

    my $bam_OB = shift;

    my $count = 0;
    for my $rg_OB ( values %{$bam_OB->{rgOB}} ){
        if( $rg_OB->{EnoughReadsBool} == 0 ){
            $count++;
        }
    }

    return $count;
}

#--- extract reads group and create related objects ---
sub extract_ReadsGroup_OB{

    my $bam_OB = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};
    my $samtools       = $parm{samtools};
    # my $Allow_Lack_RG  = $parm{Allow_Lack_RG};

    # read bam header
    open (BAMHEADER,"$samtools view -H $bam_OB->{filepath} |") || die"fail read bam: $!\n";
    while(<BAMHEADER>){
        if(/^\@RG/){ # the read group (RG)
            my ($RG_ID, $LB_ID); # add LB_ID for supporting-reads de-duplication
            if(/\sID\:(\S+).*\sLB\:(\S+)/){
                ($RG_ID, $LB_ID) = ($1, $2);
            }
            else{
                warn_and_exit "<ERROR>:\tCannot get RG_ID or LB_ID from header line $. of  bam file:\n".
                                      "\t$bam_OB->{filepath}\n";
            }
            # check recurrence
            if( exists $rgid2rgOB_Href->{$RG_ID} ){
                warn_and_exit "<ERROR>:\tThe RG_ID $RG_ID from header line $. of $bam_OB->{tissue} bam file is recurrent.\n".
                                      "\t$bam_OB->{filepath}\n";
            }
            # create reads group object (rg_OB)
            $rgid2rgOB_Href->{$RG_ID} = FuseSV::Objects::ReadsGroup_OB->new( bam_OB => $bam_OB, RG_ID => $RG_ID, LB_ID => $LB_ID );
            # link the rg_OB with the bam_OB
            $bam_OB->add_ReadsGroup_OBs( rgOB_Aref => [ $rgid2rgOB_Href->{$RG_ID} ] );
        }
    }
    close BAMHEADER;

    # once the bam lacks reads group
    if( scalar(keys %{$bam_OB->{rgOB}}) == 0 ){
        if( 0 ){ # $Allow_Lack_RG
            # future updates 
        }
        else{
            warn_and_exit "<ERROR>:\tThe $bam_OB->{tissue} bam file lacks reads group info.\n".
                                  "\t$bam_OB->{filepath}\n";
        }
    }

    # inform
    stout_and_sterr "[INFO]\tCreate reads group objects from $bam_OB->{tissue} bam file.\n"
                         ."\t$bam_OB->{filepath}\n";
}

#--- add reads group object(s) to this bam_OB ---
sub add_ReadsGroup_OBs{

    my $bam_OB = shift;
    my %parm = @_;
    my $rgOB_Aref = $parm{rgOB_Aref};

    for my $rgOB (@$rgOB_Aref){
        $bam_OB->{rgOB}->{$rgOB->{RG_ID}} = $rgOB;
    }
}

#--- get depth of certain region ---
sub get_region_depth{

    my $bam_OB = shift;
    my %parm = @_;
    my $region_refseg = $parm{refseg};
    my $region_st_pos = $parm{st_pos};
    my $region_ed_pos = $parm{ed_pos};
    my $min_baseQ = exists($parm{min_baseQ}) ? $parm{min_baseQ} : 5;
    my $min_mapQ  = exists($parm{min_mapQ})  ? $parm{min_mapQ}  : 10;
    my $no_SoftClip = $parm{no_SoftClip} || 0;
    my $outFilePrefix = $parm{out_prefix};
    my $samtools = $parm{samtools};

    unless( $region_refseg && $region_st_pos && $region_ed_pos && $outFilePrefix){
        warn_and_exit "bam_OB->get_region_depth( refseg=>?, st_pos=>?, ed_pos=>?, out_prefix=>?)\n";
    }

    unless( defined($samtools) && -e $samtools ){
        warn_and_exit "<ERROR>:\tCannot locate samtools in function bam_OB->get_region_depth().\n";
    }

    my $bamForDepth = $bam_OB->{filepath};
    my $viewRegParm = "$region_refseg:$region_st_pos-$region_ed_pos";

    # filter soft-clip
    if( $no_SoftClip ){
        my $tmp_bam = "$outFilePrefix.noSoftClip.bam";
        my $tmpbam_cmd = "( $samtools view $bam_OB->{filepath} $viewRegParm | awk '\$6 !~ /S/' | $samtools view -S -b -o $tmp_bam - ) && ( $samtools index $tmp_bam )";
        trible_run_for_success($tmpbam_cmd, 'tmpNoSoftClipBam', {esdo_Nvb=>1});
        $bamForDepth = $tmp_bam;
    }

    # get depth
    $bam_OB->{regionDepthFile} = "$outFilePrefix.depth.gz";
    my $depth_cmd = "( $samtools depth -a -r $viewRegParm -q $min_baseQ -Q $min_mapQ $bamForDepth | gzip -c > $bam_OB->{regionDepthFile} )";
    trible_run_for_success($depth_cmd, 'bamDepth', {esdo_Nvb=>1});

    # sweep
    if( $bamForDepth ne $bam_OB->{filepath} ){
        `rm -rf $bamForDepth`;
    }
}

#--- sweep region depth stat file ---
sub delete_regionDepthFile{
    my $bam_OB = shift;
    `rm -rf $bam_OB->{regionDepthFile}` if(-e $bam_OB->{regionDepthFile});
}

#--- get mutation vcf.gz in certain region ---
sub get_region_alt_vcf_gz{

    my $bam_OB = shift;
    my %parm = @_;
    my $vcfgz = $parm{vcfgz};
    my $cmd_name = $parm{cmd_name} || 'get_alt_vcf_gz';
    my $samtools = $parm{samtools};
    my $bcftools = $parm{bcftools};
    my $pos_list = $parm{pos_list};
    my $idx_ref = $parm{idx_ref};
    my $ploidy_set = $parm{ploidy_set} || 'GRCh37';
    my $no_indel = $parm{no_indel} || 0;
    my $snpGap = $parm{snpGap} || 5;
    my $IndelGap = $parm{IndelGap} || 5;
    my $min_altQual = $parm{min_altQual} || 20;
    my $min_mapQual = $parm{min_mapQual} || 30;
    my $min_posDepth = $parm{min_posDepth} || 10;
    my $min_alleleDepth = $parm{min_alleleDepth} || 5;
    my $min_strdDepth = $parm{min_strdDepth} || 2;
    my $only_hetAlt = $parm{only_hetAlt} || 0;

    # commands
    my $hetSNP_cmd = "($samtools mpileup -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR'".
                     (   defined $pos_list
                       ? " -l $pos_list "
                       : ""
                     ).
                     "-uvf $idx_ref $bam_OB->{filepath} | ".
                     "$bcftools call --ploidy $ploidy_set -cv".
                     (   $no_indel
                       ? " -V indels | "
                       : " | "
                     ).
                     "$bcftools filter -g $snpGap -G $IndelGap ".
                     "-e 'QUAL<$min_altQual ".
                     "|| DP<$min_posDepth ".
                     (   $only_hetAlt
                       ? "|| AD[0]<$min_alleleDepth ".
                         "|| DP4[0]<$min_strdDepth ".
                         "|| DP4[1]<$min_strdDepth "
                       : ""
                     ).
                     "|| AD[1]<$min_alleleDepth ".
                     "|| DP4[2]<$min_strdDepth ".
                     "|| DP4[3]<$min_strdDepth ".
                     "|| MQ<$min_mapQual' ".
                     "-O z -o $vcfgz) && ".
                     "($bcftools index $vcfgz)";
    # run
    trible_run_for_success( $hetSNP_cmd, $cmd_name );
}

#--- do stat on marker at given position ---
## maker: barcode (10x) or PE-id (Hi-C, MinIon and PacBio)
sub get_pos_marker_stat{

    my $bam_OB = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $marker = $parm{marker}; # '10x', 'Hi-C', 'MinIon', 'PacBio'
    my $samtools = $parm{samtools};

    my $PmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam_OB->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        next if( $reads_OB->is_unmap );
        # get marker of pos
        my $marker_str = ( $marker =~ /10x/i ? $reads_OB->get_10x_barc : $reads_OB->get_pid );
        $PmakerStat_Href->{$marker_str} ++ if( $marker_str ne 'NA' );
    }
    close BAM;

    return $PmakerStat_Href;
}

#--- do stat on marker of alleles at given position ---
## maker: barcode (10x) or PE-id (Hi-C, MinIon and PacBio)
sub get_allele_marker_stat{

    my $bam_OB = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $marker = $parm{marker}; # '10x', 'Hi-C', 'MinIon', 'PacBio'
    my $samtools = $parm{samtools};

    my $AmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam_OB->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        next if( $reads_OB->is_unmap );
        # get marker of pos allele
        my $alleleInfo_Href = $reads_OB->get_pos_allele( chr => $chr, pos => $pos );
        # record with marker
        my $alleleType = $alleleInfo_Href->{type};
        # if( $allele =~ /^[ACGT]$/ ){ # former version, only for refORsnv
        if( $alleleType ne 'miss' ){
            my $allele;
            if( $alleleType eq 'refORsnv' ){
                $allele = $alleleInfo_Href->{allele};
            }
            elsif( $alleleType eq 'ins' ){
                $allele = 'I,'.$alleleInfo_Href->{allele};
            }
            elsif( $alleleType eq 'del' ){
                $allele = 'D,'.$alleleInfo_Href->{delSize}.','.$alleleInfo_Href->{offset};
            }
            my $marker_str = ( $marker =~ /10x/i ? $reads_OB->get_10x_barc : $reads_OB->get_pid );
            $AmakerStat_Href->{$allele}->{$marker_str} ++ if( $marker_str ne 'NA' );
        }
    }
    close BAM;

    return $AmakerStat_Href;
}

1; ## tell the perl script the successful access of this module.
