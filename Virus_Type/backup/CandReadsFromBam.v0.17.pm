package FuseSV::Virus_Type::CandReadsFromBam;

use strict;
use warnings;
use File::Basename qw/ basename dirname /;
use List::Util qw/ min max sum /;
use BioFuse::Util::Log qw/ stout_and_sterr cluck_and_exit /;
use BioFuse::Util::Sys qw/ trible_run_for_success reset_folder /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ read_Fai /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use BioFuse::BioInfo::Objects::SeqData::FastQ_OB;
use BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
                CandReads
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Type::CandReadsFromBam';
#----- version --------
$VERSION = "0.17";
$DATE = '2019-05-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        setReadsGroup
                        CandReads
                        getCandReadsFromBam
                        PEobToCand
                        makeRGreport
                        loadRGreport
                     /;


#--- work of extracting candidate reads ---
sub CandReads{
    # start from step after current step
    if($V_Href->{stepToStart} > 1){
        &loadRGreport;
        return;
    }
    else{ # workspace
        reset_folder(folder => $V_Href->{CandFqFolder});
    }

    # load host-refseg
    read_Fai(fai => $V_Href->{host_fai}, cHf => $V_Href->{hostRefSeg});

    # reads group
    &setReadsGroup;

    # get candidate reads of each RG
    &getCandReadsFromBam;

    # write reads group report
    &makeRGreport;

    # stop at current step
    exit(0) if $V_Href->{stepToStop} == 1;
}

#--- set reads group from bam ---
sub setReadsGroup{
    # extract RG from bam
    for my $bam (@{$V_Href->{SampleOB}->bam}){
        $bam->pick_rgOB(rgid2rgOB_Href=>$V_Href->{RG});
    }
    # set RG info
    my $RG_NO = 0;
    for my $RGid (sort keys %{$V_Href->{RG}}){
        $RG_NO++;
        my $rgOB = $V_Href->{RG}->{$RGid};
        $rgOB->set_RG_NO(RG_NO => $RG_NO);
        # candidate reads: nSort.bam
        my $Cand_nSortBamPath = GetPath(filekey => 'Cand_nSortBam', CandFq_RG_ID_NO => $RG_NO);
        my $Cand_nSortBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $Cand_nSortBamPath, tag => 'Cand_nSortBam');
        $Cand_nSortBam->addTool(samtools=>$V_Href->{samtools});
        $rgOB->addAffix(Cand_nSortBam => $Cand_nSortBam);
        # candidate reads: smart fastq.gz
        my $CandFqGzPath = GetPath(filekey => 'CandFqGz', CandFq_RG_ID_NO => $RG_NO);
        my $CandFqGz = BioFuse::BioInfo::Objects::SeqData::FastQ_OB->new(filepath => $CandFqGzPath, tag => 'CandFqGzPath');
        $rgOB->addAffix(CandFqGz => $CandFqGz);
    }
}

#--- get candidate reads from each bam ---
sub getCandReadsFromBam{
    for my $bam (@{$V_Href->{SampleOB}->bam}){
        my $isNsort = $bam->isNsort;
        my $nSortBam = $bam; # initilize
        # to nSortBam if not N-sort
        unless($isNsort){
            my $nSortBamPath = GetPath(filekey => 'Tmp_nSortBam', basename => basename($bam->filepath));
            $nSortBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $nSortBamPath, tag => 'Tmp_nSortBam');
            $nSortBam->addTool(samtools=>$V_Href->{samtools});
            $bam->toNsort(nSortBam => $nSortBam);
        }
        # start rgOB Cand_nSortBam and CandFqGz
        for my $rgOB (values %{$bam->rgOB_Hf}){
            my $RGid = $rgOB->RG_ID;
            # Cand_nSortBam
            $rgOB->affix(key=>'Cand_nSortBam')->start_write;
            $rgOB->affix(key=>'Cand_nSortBam')->write(content => $_) for grep !/^\@RG/ || /ID:$RGid\s/, @{$nSortBam->header_Af};
            # CandFqGz
            $rgOB->affix(key=>'CandFqGz')->start_write;
        }
        # extract candidate reads
        ## 0x800 + 0x100 = 0x900
        $nSortBam->smartBam_PEread(viewOpt => '-F 0x900', deal_peOB_pool => 1, subrtRef => \&PEobToCand);
        # stop writing Cand_nSortBam and CandFqGz
        $_->affix(key=>'Cand_nSortBam')->stop_write for values %{$bam->rgOB_Hf};
        $_->affix(key=>'CandFqGz')->stop_write for values %{$bam->rgOB_Hf};
        # sweep temp nSort origBam
        unless($isNsort){
            my $nSortBamPath = $nSortBam->filepath;
            `rm -f $nSortBamPath`;
        }
    }
}

#--- write reads from chrPairBam to getHapBam ---
# 1. low_mapQ<$min_mapQ (host); 2. long-softclip; 3. unmapped reads; 4. non-host ref;
sub PEobToCand{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};

    for my $pe_OB (@$pe_OB_poolAf){
        my $r1_OB = $pe_OB->rOB_Af(reads_end=>1)->[0];
        my $r2_OB = $pe_OB->rOB_Af(reads_end=>2)->[0];
        # filter
        if(  (   $V_Href->{discard_DupReads} # if discard dup
              && (   $r1_OB->is_dup
                  || $r2_OB->is_dup
                 )
             )
          || (   $r1_OB->rlen < $V_Href->{min_CandRead_length} # min rlen
              || $r2_OB->rlen < $V_Href->{min_CandRead_length}
             )
        ){
            next;
        }
        # cases of candidate reads
        if(  (   $r1_OB->mapQ < $V_Href->{min_mapQ} # low mapQ
              || $r2_OB->mapQ < $V_Href->{min_mapQ}
             )
          || (   $r1_OB->biClipLen > $V_Href->{min_Slen} # long-softclip
              || $r2_OB->biClipLen > $V_Href->{min_Slen}
             )
          || (   $r1_OB->is_unmap # unmap
              || $r2_OB->is_unmap
             )
          || (   !exists $V_Href->{hostRefSeg}->{$r1_OB->mseg}
              || !exists $V_Href->{hostRefSeg}->{$r2_OB->mseg}
             )
        ){
            my $rgOB = $r1_OB->find_rgOB(rgid2rgOB_Href => $V_Href->{RG});
            for my $rOB ($r1_OB, $r2_OB){
                $rgOB->update_maxRlen(endNO => $rOB->endNO, rLen => $rOB->rlen);
                $rgOB->affix(key=>'Cand_nSortBam')->write(content => $rOB->printSAM."\n");
                $rgOB->affix(key=>'CandFqGz')->write(content => $rOB->printFQ."\n");
            }
        }
    }
}

#--- write reads group report ---
sub makeRGreport{
    open (RGRP,Try_GZ_Write($V_Href->{RG_info_report})) || die "fail $V_Href->{RG_info_report}: $!\n";
    # theme
    print RGRP '#' . join("\t", ('RG_NO', 'RG_ID', 'LB_ID', 'MAX_rlen', 'RG_nSortBam', 'RG_smartFQ') ) . "\n";
    for my $rgOB (sort {$a->RG_NO<=>$b->RG_NO} values %{$V_Href->{RG}}){
        my $RG_ID = $rgOB->RG_ID;
        my $rg_maxRlen_Hf = $rgOB->maxRlen_Hf;
        # no reads?
        if(min(values %$rg_maxRlen_Hf) == 0){
            stout_and_sterr "<WARN>\treads group $RG_ID contributes no candidate reads for virus analysis. Skip it.\n";
            delete $V_Href->{RG}->{$RG_ID};
            next;
        }
        # content
        print RGRP join ("\t",  $rgOB->RG_NO,
                                $rgOB->RG_ID,
                                $rgOB->LB_ID,
                                join(',', map {"$_:$rg_maxRlen_Hf->{$_}"} sort {$a<=>$b} keys %$rg_maxRlen_Hf),
                                basename($rgOB->affix(key=>'Cand_nSortBam')->filepath),
                                basename($rgOB->affix(key=>'CandFqGz')->filepath)
                        )."\n";
        # update global max_rlen
        $V_Href->{max_rlen} = max($V_Href->{max_rlen}, values %$rg_maxRlen_Hf);
    }
    close RGRP;
    # inform
    stout_and_sterr "[INFO]\twrite RG info report ok.\n";
}

#--- load RG info report ---
sub loadRGreport{
    open (RGRP,Try_GZ_Read($V_Href->{RG_info_report})) || die "fail read $V_Href->{RG_info_report}: $!\n";
    # theme
    (my $theme_line = lc(<RGRP>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<RGRP>){
        next if(/^\#/);
        my @info = split;
        my %k2v = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        # reads group object
        my $rgOB = BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB->new(RG_ID => $k2v{rg_id}, LB_ID => $k2v{lb_id});
        $rgOB->set_RG_NO(RG_NO => $k2v{rg_no});
        $V_Href->{RG}->{$k2v{rg_id}} = $rgOB;
        # candidate reads: nSort.bam
        my $Cand_nSortBamPath = GetPath(filekey => 'Cand_nSortBam', CandFq_RG_ID_NO => $k2v{rg_no});
        my $Cand_nSortBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $Cand_nSortBamPath, tag => 'Cand_nSortBam');
        $Cand_nSortBam->addTool(samtools=>$V_Href->{samtools});
        $rgOB->addAffix(Cand_nSortBam => $Cand_nSortBam);
        # candidate reads: smart fastq.gz
        my $CandFqGzPath = GetPath(filekey => 'CandFqGz', CandFq_RG_ID_NO => $k2v{rg_no});
        my $CandFqGz = BioFuse::BioInfo::Objects::SeqData::FastQ_OB->new(filepath => $CandFqGzPath, tag => 'CandFqGzPath');
        $rgOB->addAffix(CandFqGz => $CandFqGz);
        # max rlen
        my %maxRlen = map {(split /:/)[0,1]} split /,/, $k2v{max_rlen};
        $rgOB->update_maxRlen(endNO => $_, rLen => $maxRlen{$_}) for keys %maxRlen;
        # update global max_rlen
        $V_Href->{max_rlen} = max($V_Href->{max_rlen}, values %maxRlen);
        # check file existence
        if(    !-e $CandFqGzPath
            && (   $V_Href->{stepToStart} >= 2
                || $V_Href->{stepToStop} <= 3
               )
        ){
            stout_and_sterr "<WARN>\tCannot find SmartFq for ReadGroup (RG): $k2v{rg_id}. 'cand_virus' will construct it.\n";
        }
    }
    close RGRP;
    # inform
    stout_and_sterr "[INFO]\tload RG info report ok.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
