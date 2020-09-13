package FuseSV::Virus_Type::VirusInitAlign;

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/ basename dirname /;
use File::Spec qw/ catfile /;
use List::Util qw/ min max sum /;
use POSIX qw/ ceil /;
use Getopt::Long;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ reset_folder /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::Alignment::BWA qw/ PEfqToSortBam /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Objects::Virus_OB;
use FuseSV::Virus_Type::PrepareExtendVirusRef qw/ load_virusBasicInfo make_virusAlignRef load_virusRefReport make_virusRefReport /;

use FuseSV::Virus_Type::InitialMapSmartFQonVirus;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              VirusInitAlign
              map_CandFqToGivenVirus
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Type::VirusInitAlign';
#----- version --------
$VERSION = "0.13";
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
                        VirusInitAlign
                        prep_virusAlignRef
                        map_CandFqToVirus
                        map_CandFqToGivenVirus
                        load_virusAlignReport
                        make_virusAlignReport
                     /;

#--- based on the virus info, prepare virus ref for first alignment ---
sub VirusInitAlign{
    # load virus info from database
    load_virusBasicInfo;

    # start from step after current step
    if($V_Href->{stepToStart} > 2){
        # load virus initial alignment report
        &load_virusAlignReport;
        return;
    }
    else{ # workspace
        reset_folder(folder => $V_Href->{virusInitAlignFolder});
    }

    # make/link the virus extend? ref genome
    &prep_virusAlignRef;

    # initial alignment of smartFQ to virus ref
    # FuseSV::Virus_Type::InitialMapSmartFQonVirus->initial_align_smfq_on_all_vref;
    &map_CandFqToVirus;

    # make virus initial alignment report
    &make_virusAlignReport;

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 2;
}

#--- prepare virus extend? ref genome ---
sub prep_virusAlignRef{
    my $extend_length = ceil($V_Href->{max_rlen} / 100) * 100;
    my $extend_vref_db_folder = FuseSV::GetPath->GetPath(filekey => 'extend_vref_db_folder', vref_extLen => $extend_length);
    if(-e $extend_vref_db_folder){ # if database has
        `ln -sf $extend_vref_db_folder $V_Href->{virusInitAlignRefFolder}`;
        # load virusRefReport to make up certain attributes
        load_virusRefReport(
            report => File::Spec->catfile($V_Href->{virusInitAlignRefFolder}, 'virus.ref.info'),
            hasmap => 0
        );
    }
    else{
        make_virusAlignRef(
            virus_ref_dir => $V_Href->{virusInitAlignRefFolder},
            extend_length => $extend_length
        );
    }
}

#--- map smartFQ of all reads group to all virus ref ---
sub map_CandFqToVirus{
    reset_folder(folder => $V_Href->{virusInitAlignBamFolder});
    # fork manager
    my ($pm, $fork_DO) = &forkSetting; # FuseSV::LoadOn
    # read paired bam files and split
    for my $virus_OB (sort {$a->no<=>$b->no} values %{$V_Href->{virusRefSeg}}){
        # fork job starts
        my $fork_tag = $virus_OB->id . ' init_align';
        my $filePref = File::Spec->catfile($V_Href->{virusInitAlignBamFolder}, 'CandFqMapToVrefNO'.$virus_OB->no);
        my $tmpList = "$filePref.initBam.tmp.list";
        if($fork_DO){ $pm->start($fork_tag) and next }
        eval{
            my $init_bam = &map_CandFqToGivenVirus(filePref => $filePref, virusRef => $virus_OB->ana_ref);
            `echo "$_" > $tmpList` for $init_bam->filepath;
        };
        if($@){
            if($fork_DO){ warn $@; $pm->finish(1) }
            else{ cluck_and_exit $@; }
        }
        # fork job finishes
        if($fork_DO){ $pm->finish(0) }
        # collect alignBam
        chomp(my $initBamPath = `cat $tmpList`);
        my $init_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $initBamPath, tag => "initBam");
        $init_bam->addTool(samtools => $V_Href->{samtools});
        $virus_OB->set_init_bam(bam => $init_bam);
        `rm -f $tmpList`;
    }
    # collect fork jobs
    if($fork_DO){ $pm->wait_all_children }
}

#--- align candidate fastq to given virus reference ---
sub map_CandFqToGivenVirus{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $filePref = $parm{filePref};
    my $virusRef = $parm{virusRef};

    # each rg
    my @cSortBam;
    for my $rgOB (sort {$a->RG_NO<=>$b->RG_NO} values %{$V_Href->{RG}}){
        my $bamPath = "$filePref.RG_NO".$rgOB->RG_NO.'.cSort.bam';
        push @cSortBam, BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $bamPath, tag => $rgOB->RG_ID.' cSortBam');
        $cSortBam[-1]->addTool(samtools => $V_Href->{samtools});
        my $RG_para = $rgOB->RG_para(sampleID => $V_Href->{sample_ID});
        PEfqToSortBam(
            ref => $virusRef,
            fqAf => [$rgOB->affix(key=>'Cand_smartFqGz')],
            bam => $cSortBam[-1],
            bwa => $V_Href->{bwa},
            func => 'mem',
            RG_para => "-R '$RG_para'",
            samtools => $V_Href->{samtools},
            sort_mode => 'c',
            discBuPair => 1
        );
    }
    # merge
    my $mergeBamPath = $filePref.'.merge.bam';
    my $mergeBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $mergeBamPath, tag => 'mergeBam');
    $mergeBam->addTool(samtools => $V_Href->{samtools});
    if(@cSortBam == 1){
        my $cSortBamPath = $cSortBam[0]->filepath;
        `mv -f $cSortBamPath $mergeBamPath`;
    }
    else{
        $mergeBam->merge(bamAf => \@cSortBam);
    }
    $mergeBam->toIndex;
    $_->delete_file for @cSortBam;
    # mark dup
    my $mkdupBamPath = $filePref.'.mkdup.bam';
    my $mkdupBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $mkdupBamPath, tag => 'mkdupBam');
    $mkdupBam->addTool(samtools => $V_Href->{samtools});
    $mergeBam->toMarkdup(markdupBam => $mkdupBam);
    $mkdupBam->toIndex;
    $mergeBam->delete_file;

    return $mkdupBam;
}

#--- load virus initial alignment report ---
# format is similar to virusRefReport
sub load_virusAlignReport{
    load_virusRefReport(
        report => $V_Href->{virusInitAlignReport},
        hasmap => 1
    );
    # inform
    stout_and_sterr "[INFO]\tload virus initial alignment report.\n";
}

#--- make virus initial alignment report ---
# format is similar to virusRefReport
sub make_virusAlignReport{
    make_virusRefReport(
        report => $V_Href->{virusInitAlignReport},
        hasmap => 1
    );
    # inform
    stout_and_sterr "[INFO]\tmake virus initial alignment report.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
