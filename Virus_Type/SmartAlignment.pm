package FuseSV::Virus_Type::SmartAlignment;

use strict;
use warnings;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation qw/warn_and_exit trible_run_for_success stout_and_sterr/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use FuseSV::LoadOn;
use FuseSV::Virus_Type::CandReads;
use Parallel::ForkManager;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              SmartFq2MergeBam_against_given_ref
              smartfq_singleRG_bwa
              merge_bam_mkdup
              GATK_relan_Samtools_calmd
              mark_additional_alignments
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'SmartAlignment';
#----- version --------
$VERSION = "0.08";
$DATE = '2017-06-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        SmartFq2MergeBam_against_given_ref
                        smartfq_singleRG_bwa
                        merge_bam_mkdup
                        mark_additional_alignments
                        Judge_Dup_for_SuppleM
                        GATK_relan_Samtools_calmd
                     /;

#--- align smartfq of all RGs to one given ref ---
sub SmartFq2MergeBam_against_given_ref{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
    my $OutBam_pref    = $Option_Href->{OutBamPref};
    my $Fork_Number    = $Option_Href->{fork_number} || 1;

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $RG_Info_Href = $V_Href->{RG_Info_Href};

    # potential multiple fork operations, initialization
    my $multiple_fork;
    if($Fork_Number > 1){
        $multiple_fork = new Parallel::ForkManager($Fork_Number);
        # inform
        stout_and_sterr `date`."[INFO]:\tOperate in multiple forks mode. Fork Number: $Fork_Number.\n";
    }

    # deal each RG
    for my $RG_ID (sort keys %$RG_Info_Href){
        # potential multiple fork operations, child fork starts
        if($Fork_Number > 1){
            $multiple_fork->start and next;
        }

        my $RG_NO = $RG_Info_Href->{$RG_ID}->{RG_NO};
        my $SmartFqGz = $RG_Info_Href->{$RG_ID}->{CandFqGz};
        # check the smartfq existence
        if(!-e $SmartFqGz){
            my $nSortBam = $RG_Info_Href->{$RG_ID}->{Cand_nSortBam};
            # samtools to convert bam to fastq
            stout_and_sterr `date`."[INFO]\tstarts extracting fastq from nSorted candidate bam for RG $RG_ID.\n";
            FuseSV::Virus_Type::CandReads->extract_fastq_from_nSortBam(
                    {
                        nSortBam => $nSortBam,
                        fastq_gz => $SmartFqGz
                    }
                );
            stout_and_sterr `date`."[INFO]\textract fastq from nSorted candidate bam for RG $RG_ID ok.\n";
        }
        # align the smartfq to ref
        my $RG_OutBam_pref = "$OutBam_pref.RG_NO$RG_NO";
        my $RG_bam = 
            &smartfq_singleRG_bwa(
                    {
                        CandFqGz    =>  $SmartFqGz,
                        RefForBwa   =>  $Ref_IndexByBwa,
                        OutBamPref  =>  $RG_OutBam_pref,
                        RG_ID       =>  $RG_ID
                    }
                );
        # bam list of this RG
        my $RG_bam_list = "$RG_OutBam_pref.bam.list";
        `rm -rf $RG_bam_list` if(-e $RG_bam_list);
        if($RG_bam ne 'EMPTY'){
            `echo $RG_bam > $RG_bam_list`;
        }

        # potential multiple fork operations, child fork finishs
        if($Fork_Number > 1){
            $multiple_fork->finish;
        }
    }

    # potential multiple fork operations, wait all childern forks
    if($Fork_Number > 1){
        $multiple_fork->wait_all_children;
    }

    # gather all RG bam(s)
    my @BamList = glob("$OutBam_pref.RG_NO*.bam.list");
    # all RG to merged mkdup bam
    if( scalar(@BamList) != 0 ){
        chomp(my @RG_Bam = `cat @BamList`);
        my $mkdup_bam = 
            &merge_bam_mkdup(
                    {
                        SortBamAref =>  \@RG_Bam,
                        RefForBwa   =>  $Ref_IndexByBwa,
                        OutBamPref  =>  $OutBam_pref
                    }
                );
        return $mkdup_bam;
    }
    else{
        return 'EMPTY';
    }
}

#--- do BWA+samtools pipeline on given smartFq and ref ---
sub smartfq_singleRG_bwa{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
    my $OutBam_pref    = $Option_Href->{OutBamPref};
    my $RG_ID          = $Option_Href->{RG_ID};
    my $SmartFqGz      = $Option_Href->{CandFqGz} || 'NULL';
    # additional functions
    my $Map_func = $Option_Href->{MapFunc} || 'mem'; # generally, >70bp uses mem, <=70bp uses aln
    my $FqGz_1 = $Option_Href->{FqGz1} || 'NULL';
    my $FqGz_2 = $Option_Href->{FqGz2} || 'NULL';
    my $Max_ins = $Option_Href->{MaxIns} || 500; # just for 'aln' function, '-a' parameter
    my $Save0Bam = $Option_Href->{Save0Bam} || 0;
    my $F_filter_flag = exists($Option_Href->{F_filter}) ? "-F $Option_Href->{F_filter}" : ''; # use for samtools' '-F' parameter filtration

    # function check
    ## lack one
    if(      !exists($Option_Href->{CandFqGz})
        && ( !exists($Option_Href->{FqGz1}) || !exists($Option_Href->{FqGz2}) )
       ){
        warn_and_exit `date`."<ERROR>:\tYou must input single 'CandFqGz' OR 'FqGz1' && 'FqGz2'.\n";
    }
    ## duplicate
    if(      exists($Option_Href->{CandFqGz})
        && ( exists($Option_Href->{FqGz1}) || exists($Option_Href->{FqGz2}) )
       ){
        warn_and_exit `date`."<ERROR>:\tYou cannot input 'CandFqGz' together with 'FqGz1' && 'FqGz2'.\n";
    }
    ## aln function needs more parameters
    if( $Map_func eq 'aln' && !exists($Option_Href->{MaxIns}) ){
        stout_and_sterr `date`."<WARN>:\tYou select bwa aln for mapping, but lacks maximum insertsize.\n"
                              ."\tSo automatically set it as 500.\n";
    }
    ## aln requires two fq-files
    if(      $Map_func eq 'aln'
        && ( !exists($Option_Href->{FqGz1}) || !exists($Option_Href->{FqGz2}) )
       ){
        warn_and_exit `date`."<ERROR>:\tYou must use 'aln' together with 'FqGz1' && 'FqGz2'.\n";
    }

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # bwa mem to orig.bam
    # use sample_ID as Lib_ID for de-duplication in whole data
    # Platform must be recognized one, or GATK cannot run
    # -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
    # -U INT        penalty for an unpaired read pair [17]
    # -Y, soft-clip for supplementary alignments
    # -p, smart pairing
    # use awk to filter the bwa-mem standout sam format pipe.
    #    awk '
    #            $0~/^\@/                                    # include SAM header
    #         || (    !(and($2,0x4) && and($2,0x8))          # Cannot be u and U simultaneously
    #              && xor(and($2,0x40), and($2,0x80))        # must have and only have one of end_NO 1 and 2
    #              && and($2,0x1)                            # must be paired sequencing, mean 'p'
    #            )
    #        '
    my $orig_bam = "$OutBam_pref.orig.bam";
    my $OrigBam_command;
    my $RG_para = "'\@RG\\tID:$RG_ID\\tLB:$V_Href->{sample_ID}\\tPL:ILLUMINA\\tCN:NULL\\tPU:$V_Href->{sample_ID}\\tSM:$V_Href->{sample_ID}'";
    # my $OrigBam_command = "$V_Href->{bwa} mem -t 8 -Y -p -R '\@RG\\tID:$RG_ID\\tLB:$V_Href->{sample_ID}\\tPL:ILLUMINA\\tCN:NULL\\tPU:$V_Href->{sample_ID}\\tSM:$V_Href->{sample_ID}' $Ref_IndexByBwa $SmartFqGz | awk '\$0~/^\@/ || (!(and(\$2,0x4) && and(\$2,0x8)) && xor(and(\$2,0x40), and(\$2,0x80)) && and(\$2,0x1))' | $V_Href->{samtools} view -b -S -T $Ref_IndexByBwa - > $orig_bam";
    if( $Map_func eq 'mem' ){ # bwa mem
        $OrigBam_command =    "$V_Href->{bwa} mem -t 8 -Y "
                            . "-R $RG_para "
                            . ( ($SmartFqGz ne 'NULL') ? '-p ' : '' )
                            . "$Ref_IndexByBwa "
                            . ( ($SmartFqGz ne 'NULL') ? "$SmartFqGz " : "$FqGz_1 $FqGz_2 " )
                            . "2>/dev/null "
                            . "| awk '\$0~/^\@/ || (!(and(\$2,0x4) && and(\$2,0x8)) && xor(and(\$2,0x40), and(\$2,0x80)) && and(\$2,0x1))' "
                            . "| $V_Href->{samtools} view -b -S $F_filter_flag -T $Ref_IndexByBwa -o $orig_bam - 2>/dev/null";
    }
    else{ # bwa aln
        my $sai_1 = "$OutBam_pref.orig.end_1.sai";
        my $sai_2 = "$OutBam_pref.orig.end_2.sai";
        $OrigBam_command =    "($V_Href->{bwa} aln -L -m 800000 -l 27 -i 10 -k 3 -t 8 -e 10 -f $sai_1 $Ref_IndexByBwa $FqGz_1 2>/dev/null) && "
                            . "($V_Href->{bwa} aln -L -m 800000 -l 27 -i 10 -k 3 -t 8 -e 10 -f $sai_2 $Ref_IndexByBwa $FqGz_2 2>/dev/null) && "
                            . "("
                            . "$V_Href->{bwa} sampe -a $Max_ins -r $RG_para $Ref_IndexByBwa $sai_1 $sai_2 $FqGz_1 $FqGz_2 2>/dev/null "
                            . "| awk '\$0~/^\@/ || (!(and(\$2,0x4) && and(\$2,0x8)) && xor(and(\$2,0x40), and(\$2,0x80)) && and(\$2,0x1))' "
                            . "| $V_Href->{samtools} view -b -S $F_filter_flag -T $Ref_IndexByBwa -o $orig_bam - 2>/dev/null"
                            . ") && "
                            . "(rm -rf $sai_1 $sai_2)";
    }
    trible_run_for_success($OrigBam_command, 'OrigBam', {esdo_Nvb=>1});
    # inform
    stout_and_sterr `date`."[INFO]\tBWA $Map_func finished.\n"
                          ."\tRef: $Ref_IndexByBwa\n"
                          ."\tFq:  $_\n" for ( ($SmartFqGz ne 'NULL') ? ($SmartFqGz) : ($FqGz_1, $FqGz_2) );

    # check the bam, not empty
    my $CheckBam_command = "$V_Href->{samtools} view $orig_bam | head | wc -l";
    chomp(my $LineNumber = `$CheckBam_command`);
    if($LineNumber == 0){
        if($Save0Bam){
            return "EMPTY:$orig_bam";
        }
        else{
            `rm -rf $orig_bam`;
            return 'EMPTY';
        }
    }

    # sort by read-name and fixmate
    my $fixmate_bam = "$OutBam_pref.fixmate.bam";
    my $fixmate_command = "$V_Href->{samtools} sort --threads 4 -m 1G -n -T $fixmate_bam.pre_sort_tmp $orig_bam | $V_Href->{samtools} fixmate -O bam - $fixmate_bam";
    trible_run_for_success($fixmate_command, 'FixmateBam', {esdo_Nvb=>1});
    `rm -rf $orig_bam`;
    # inform
    stout_and_sterr `date`."[INFO]\tSamTools fixmate finished.\n";

    # coordinate-sort and index
    my $sort_bam = "$OutBam_pref.sort.bam";
    my $SortBam_command = "$V_Href->{samtools} sort --threads 4 -m 1G -T $sort_bam.sort_tmp -o $sort_bam $fixmate_bam && $V_Href->{samtools} index $sort_bam";
    trible_run_for_success($SortBam_command, 'SortBam', {esdo_Nvb=>1});
    `rm -rf $fixmate_bam`;
    # inform
    stout_and_sterr `date`."[INFO]\tSamTools sort bam finished.\n";

    return $sort_bam;
}

#--- do Merge-Dup pipeline on given sort.bam and ref ---
sub merge_bam_mkdup{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $Sort_bam_Aref  = $Option_Href->{SortBamAref};
    my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
    my $OutBam_pref    = $Option_Href->{OutBamPref};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # merge RGs' bams
    my $merge_bam = "$OutBam_pref.merged.bam";
    my $MergeBam_command = "$V_Href->{samtools} merge -c -p --threads 4 $merge_bam @$Sort_bam_Aref";
    trible_run_for_success($MergeBam_command, 'MergeBam', {esdo_Nvb=>1});
    `rm -rf $_*` for @$Sort_bam_Aref;
    # inform
    stout_and_sterr `date`."[INFO]\tSamTools merge bam finished.\n";

    # mark dup and index, picard
    my $mkdup_bam = "$OutBam_pref.mkdup.bam";
    my $metrics = "$OutBam_pref.mkdup.bam.mat";
    my $mkdup_tmpdir = "$OutBam_pref.mkdup.bam.tmpdir";
    `rm -rf $mkdup_tmpdir`;
    `mkdir -p $mkdup_tmpdir`;
    my $MkdupBam_command = "($V_Href->{java} -Xmx$V_Href->{mem_for_jvm} -jar $V_Href->{picard_jar} MarkDuplicates REMOVE_DUPLICATES=false I=$merge_bam O=$mkdup_bam METRICS_FILE=$metrics TMP_DIR=$mkdup_tmpdir 2>/dev/null) && $V_Href->{samtools} index $mkdup_bam";
    trible_run_for_success($MkdupBam_command, 'MkdupBam', {esdo_Nvb=>1});
    `rm -rf $merge_bam`;
    `rm -rf $metrics $mkdup_tmpdir`;
    # inform
    stout_and_sterr `date`."[INFO]\tPicard MarkDuplicates bam finished.\n";

    # specially deal with additional alignments, 0x100, 0x800
    my $M900_mkdup_bam = "$OutBam_pref.mkdup.M900.bam";
    &mark_additional_alignments(
        {
            source_bam => $mkdup_bam,
            output_bam => $M900_mkdup_bam
        }
    );
    `rm -rf $mkdup_bam $mkdup_bam.bai`;
    `mv -f $M900_mkdup_bam $mkdup_bam`; # rename
    `$V_Href->{samtools} index $mkdup_bam`;
    # inform
    stout_and_sterr `date`."[INFO]\tAdditional alignments MarkDuplicates bam finished.\n";

    return $mkdup_bam;
}

#--- make up for picard, mkdup for 0x100 and 0x800
sub mark_additional_alignments{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $source_bam = $Option_Href->{source_bam};
    my $output_bam = $Option_Href->{output_bam};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    open (my $obam_fh,"| $V_Href->{samtools} view -b -S -o $output_bam")|| die"fail samtools write $output_bam: $!\n";
    my %SuppleM_Info;
    my $last_mpos = 0;
    my $record_NO = 0;
    open (SAM,"$V_Href->{samtools} view -h $source_bam |") || die"fail samtools view header of $source_bam: $!\n";
    while(<SAM>){
        # SAM header
        if(/^@/){
            print {$obam_fh} $_;
            next;
        }
        # Read alignment info
        my @line_ele = split;
        my $pid   = $line_ele[0];
        my $flag  = $line_ele[1];
        my $mseg  = $line_ele[2];
        my $mpos  = $line_ele[3];
        my $mqual = $line_ele[4];
        my $cigar = $line_ele[5];
        my $p_mseg= $line_ele[6];
        my $p_mpos= $line_ele[7];
        # deal with reads group of last mpos
        if(    $last_mpos != 0
            && $last_mpos != $mpos
            && $record_NO != 0
        ){
            &Judge_Dup_for_SuppleM({out_fh=>$obam_fh, minfo_Href=>\%SuppleM_Info});
            # sweep
            %SuppleM_Info = ();
            $record_NO = 0;
        }
        # already judged as duplication
        if( $flag & 0x400 ){
            print {$obam_fh} $_;
        }
        # additional alignment
        elsif(    $flag & 0x100 # secondary alignment
               || $flag & 0x800 # supplementary alignment
        ){
            # NOTE that the alignInfo does not include the mapping-orientation!!!
            my $alignInfo = join(',', $mseg, $mpos, $cigar, $p_mseg, $p_mpos);
            $SuppleM_Info{$alignInfo}->{$pid} = {
                line_Aref => \@line_ele,
                mqual => $mqual
            };
            # count in
            $record_NO ++;
            # update
            $last_mpos = $mpos;
        }
        # others has been dealt by previous tool(e.g., picard)
        else{
            print {$obam_fh} $_;
        }
    }
    close SAM;
    # last group if possible
    &Judge_Dup_for_SuppleM({out_fh=>$obam_fh, minfo_Href=>\%SuppleM_Info}) if($record_NO != 0);
    close $obam_fh;
}

#--- judge duplication from additional alignments
sub Judge_Dup_for_SuppleM{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $out_fh = $Option_Href->{out_fh};
    my $minfo_Href = $Option_Href->{minfo_Href};

    for my $alignInfo (sort keys %$minfo_Href){
        my $Info_Href = $minfo_Href->{$alignInfo};
        my @pid = sort keys %$Info_Href;
        # find the best pid
        my $best_mqual = -1; # initial
        my $best_pid;
        for my $pid (@pid){
            # always keep the pid with best mapping-quality
            if( $Info_Href->{$pid}->{mqual} > $best_mqual ){
                $best_mqual = $Info_Href->{$pid}->{mqual};
                $best_pid = $pid;
            }
        }
        # output each pid
        for my $pid (@pid){
            my $line_Aref = $Info_Href->{$pid}->{line_Aref};
            if( $pid eq $best_pid ){
                # release from dup status when necessary
                $line_Aref->[1] -= 0x400 if( $line_Aref->[1] & 0x400 );
            }
            else{
                # mark as dup anyway
                $line_Aref->[1] += 0x400 unless( $line_Aref->[1] & 0x400 );
            }
            print {$out_fh} join("\t",@$line_Aref)."\n";
        }
    }
}

#--- do GATK-realn pipeline on given mkdup bam and ref ---
sub GATK_relan_Samtools_calmd{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

    my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
    my $OutBam_pref    = $Option_Href->{OutPref};
    my $SkipCalmd_bool = $Option_Href->{SkipCalmd} || 0;

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    my $mkdup_bam = $Option_Href->{mkdup_bam} || "$OutBam_pref.mkdup.bam";

    # GATK realn
    my $realn_bam = "$OutBam_pref.realn.bam";
    (my $realn_bam_bai = $realn_bam) =~ s/bam$/bai/;
    my $intervals = "$OutBam_pref.realn.bam.intervals";
    my $RealnBam_command = "($V_Href->{java} -Xmx$V_Href->{mem_for_jvm} -jar $V_Href->{gatk_jar} -T RealignerTargetCreator -R $Ref_IndexByBwa -I $mkdup_bam -o $intervals 2>/dev/null) && ($V_Href->{java} -Xmx$V_Href->{mem_for_jvm} -jar $V_Href->{gatk_jar} -T IndelRealigner -R $Ref_IndexByBwa -I $mkdup_bam -o $realn_bam -targetIntervals $intervals -maxInMemory 300000 -maxReads 100000 -l INFO 2>/dev/null)";
    trible_run_for_success($RealnBam_command, 'RealnBam', {esdo_Nvb=>0});
    # `rm -rf $mkdup_bam*`;
    `rm -rf $intervals`;
    # inform
    stout_and_sterr `date`."[INFO]\tGATK realn bam finished.\n";

    if($SkipCalmd_bool){
        # return GATK realn bam
        `mv -f $realn_bam_bai $realn_bam.bai`;
        return $realn_bam;
    }
    else{
        # samtools calmd
        my $calmd_bam = "$OutBam_pref.calmd.bam";
        my $CalmdBam_command = "($V_Href->{samtools} calmd -b -r $realn_bam $Ref_IndexByBwa >$calmd_bam 2>/dev/null) && ($V_Href->{samtools} index $calmd_bam)";
        trible_run_for_success($CalmdBam_command, 'CalmdBam', {esdo_Nvb=>1});
        `rm -rf $realn_bam $realn_bam_bai`;
        # inform
        stout_and_sterr `date`."[INFO]\tSamTools calmd bam finished.\n";
        # return SamTools calmd bam
        return $calmd_bam;
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
