package FuseSV::Virus_Type::ConstructVirusDatabase;

use strict;
use warnings;
use FindBin qw/ $RealBin /;
use Getopt::Long;
use List::Util qw/ max min /;
use POSIX qw/ ceil /;
use BioFuse::Util::Web qw/ urlToHtmlText /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success reset_folder /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ write_fasta_file BWA_index_fasta Faidx_Dict_fasta /;
use BioFuse::BioInfo::Objects::GeneAnno::GTF_OB;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Objects::Virus_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              ConstructVirusDatabase
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Type::ConstructVirusDatabase';
#----- version --------
$VERSION = "0.11";
$DATE = '2020-08-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        ConstructVirusDatabase
                        load_virus_ncbi_list
                        get_virus_web_data
                        analyze_virus_ncbi_web_data
                        query_virus_html
                        write_virus_info
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} make_virusdb <[Options]>
     
     Options:

       # Inputs and Outputs #
        -vlist  [s]  virus list for online NCBI query. <required>
                      This list contains a header prefixed by '\033[1m \b# \b\033[0m'.
                      TAB delimited Header: #subtype \\t GeneBank_ID
                      'subtype' only allows letter (A-Z and a-z), number (0-9) and underline ('_').
        -vpf    [s]  virus database folder. <required>
                      Folder to store virus database, PLS make this folder manually.
                      NOTE: the \033[4m \bWRITE \b\033[0m authority is required.

       # Software Required #
        -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
        -stl    [s]  SAMtools, minimum version: 1.3. <required>

       # Options #
        -maxlen [i]  Maximum length of single virus allowed. [20000]

        -h|help   Display this help info.
 
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
                cluck_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ virus_db_list => '' ],
            # software
            [ bwa => '' ],
            [ samtools => '' ],
            # options
            [ db_virus_max_length => 20000 ],
            # basic
            [ ncbi_url_pf => 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=' ],
            [ Start_codon => ['ATG'] ],

            #--- containers ---#
            # id(subtype) -> virusOB
            # see FuseSV::Objects::Virus_OB
            [ virusRefSeg => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['bwa'],
                                  ['samtools'],
                                  ['virus_db_list'],
                                  ['virus_db_folder'] ] ]
        );
}

#--- get options from command line
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-vlist:s"  => \$V_Href->{virus_db_list},
        "-vpf:s"    => \$V_Href->{virus_db_folder},
        # software
        "-bwa:s"    => \$V_Href->{bwa},
        "-stl:s"    => \$V_Href->{samtools},
        # options
        "-maxlen:i" => \$V_Href->{db_virus_max_length},
        # help
        "-h|help"   => \$V_Href->{HELP},
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || (!$V_Href->{virus_db_list} || !-e $V_Href->{virus_db_list})
             || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder} || !-w $V_Href->{virus_db_folder})
             || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
             || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
            );
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub ConstructVirusDatabase{
    # read virus subtype list
    &load_virus_ncbi_list;

    # obtain information from web of all virus
    &get_virus_web_data;

    # generate virus info list
    &write_virus_info;
}

#--- load virus list for database ---
sub load_virus_ncbi_list{
    open (VLIST,Try_GZ_Read($V_Href->{virus_db_list})) || die "fail read $V_Href->{virus_db_list}: $!\n";
    # theme
    #subtype        GeneBank_ID
    (my $theme_line = lc(<VLIST>)) =~ s/^#//;
    if( $theme_line !~ /subtype|genebank_id/i){
        cluck_and_exit "<ERROR>\tvirus list lacks header.\n";
    }
    my @theme_tag = split /\s+/, $theme_line;
    while(<VLIST>){
        next if(/^\#/);
        my @info = split;
        my %vOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        if( $vOb{subtype} =~ /\W/ ){
            cluck_and_exit "<ERROR>\tsubtype $vOb{subtype} uses unallowed char.\n";
        }
        my $virus_OB = FuseSV::Objects::Virus_OB->new(id => $vOb{subtype});
        $V_Href->{virusRefSeg}->{$vOb{subtype}} = $virus_OB;
        $virus_OB->set_genebank(genebank => $vOb{genebank_id});
        $virus_OB->set_ncbi_url(url => $V_Href->{ncbi_url_pf}.$vOb{genebank_id});
    }
    close VLIST;
    # inform
    stout_and_sterr "[INFO]\tload virus list for database.\n";
}

#--- get each virus webdata ---
sub get_virus_web_data{
    my $NCBI_HtmlText_folder = GetPath(filekey => 'NCBI_HtmlText_folder');
    my $Single_virus_folder = GetPath(filekey => 'single_virus_data_folder');
    reset_folder(folder => $_) for $NCBI_HtmlText_folder, $Single_virus_folder;

    for my $virus_subtype (sort keys %{$V_Href->{virusRefSeg}}){
        # web query virus info
        &query_virus_html(virus_subtype => $virus_subtype);
        # analysis
        &analyze_virus_ncbi_web_data(virus_subtype => $virus_subtype);
        # check virus availability
        my $virus_length = $V_Href->{virusRefSeg}->{$virus_subtype}->length;
        if( $virus_length > $V_Href->{db_virus_max_length} ){
            cluck_and_exit "<WARN>\tthe original virus genome length ($virus_length nt) of $virus_subtype exceeds allowed maximum value ($V_Href->{db_virus_max_length} nt).\n";
        }
        # wait for while, avoid robot detection
        sleep 3;
    }

    # merge virus annotation
    my $virus_anno_folder = GetPath(filekey => 'virus_gene_anno_folder');
    reset_folder(folder => $virus_anno_folder);
    my $merged_gpsl = GetPath(filekey => 'virus_orig_gene_psl');
    `cat $Single_virus_folder/*.gene.psl > $merged_gpsl`;
    my $merged_tpsl = GetPath(filekey => 'virus_orig_trans_psl');
    `cat $Single_virus_folder/*.trans.psl > $merged_tpsl`;

    # merge virus reference
    my $virus_orig_ref_folder = GetPath(filekey => 'virus_orig_ref_folder');
    reset_folder(folder => $virus_orig_ref_folder);
    my $merged_vfa = GetPath(filekey => 'virus_orig_genome_fa');
    `cat $Single_virus_folder/*.complete.genome.fa > $merged_vfa`;
    # create virus BWA reference for alignment
    BWA_index_fasta(FaFile => $merged_vfa, bwa => $V_Href->{bwa});
    Faidx_Dict_fasta(FaFile => $merged_vfa, samtools => $V_Href->{samtools});

    # inform
    stout_and_sterr "[INFO]\tmerge virus PSL and fasta (index,faidx,dict) OK.\n";
}

#--- query virus html ---
sub query_virus_html{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $virus_subtype = $parm{virus_subtype};

    # web query
    my $ncbi_url = $V_Href->{virusRefSeg}->{$virus_subtype}->ncbi_url;
    my $html_text = urlToHtmlText( url => $ncbi_url );
    ## check
    while( $html_text !~ /200 OK/ ){
        cluck_and_exit "[ERROR]\tfail to connect virus $virus_subtype [ $ncbi_url ].\n";
    }
    ## extract the info
    $html_text =~ s/<pre class="genbank">//g;
    $html_text =~ s/<\/pre>//g;
    $html_text =~ s/<a [^>]+>//g;
    $html_text =~ s/<\/a>//g;
    # ## record
    # $V_Href->{virusRefSeg}->{$virus_subtype}->{html_text} = $html_text;

    # output text
    my $ncbi_txt = GetPath(filekey => 'virus_NCBI_page_txt', virus_subtype => $virus_subtype);
    $V_Href->{virusRefSeg}->{$virus_subtype}->set_ncbi_txt(txt => $ncbi_txt);
    open (NCBIPAGE, Try_GZ_Write($ncbi_txt)) || die "cannot write $ncbi_txt: $!\n";
    print NCBIPAGE $html_text;
    close NCBIPAGE;

    # inform
    stout_and_sterr "[INFO]\tget virus $virus_subtype NCBI data.\n";
}

#--- analyze the ncbi html text ---
sub analyze_virus_ncbi_web_data{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $VIRUS_TYPE = $parm{virus_subtype};
    # files
    my $single_virus_folder = GetPath(filekey => 'single_virus_data_folder');
    my $COMP_fa  = File::Spec->catfile($single_virus_folder, "$VIRUS_TYPE.complete.genome.fa");
    my $ANAY_fa  = File::Spec->catfile($single_virus_folder, "$VIRUS_TYPE.extended.for.analysis.fa");
    my $gtfPath  = File::Spec->catfile($single_virus_folder, "$VIRUS_TYPE.extended.gtf"); # based on extended.for.analysis.fa
    my $gPSLpath = File::Spec->catfile($single_virus_folder, "$VIRUS_TYPE.extended.gene.psl");
    my $tPSLpath = File::Spec->catfile($single_virus_folder, "$VIRUS_TYPE.extended.trans.psl");
    # genome data
    my $VERSION = 'NULL';
    my $CIRCLE = 0;
    my $LENGTH = -1;
    my $Genome_Seq = '';
    my $Extend_len = 1000; # for the general longest insertsize (800+3*SD) of WGS
    my $Ref_Name_for_Analysis;
    my $Ref_Name_for_Complete = $VIRUS_TYPE;
    # gene data
    my %gene_record;
    my ($gene_NO,$trans_NO,$exon_wNO,$protein_NO) = (0,1,1,1);
    my $gene_biotype = 'protein_coding'; # only read CDS, however maybe 3'/5'-partial
    my @GTF_OutPut = ();
    # analyze ncbi webpae data
    my $NCBI_TXT = $V_Href->{virusRefSeg}->{$VIRUS_TYPE}->ncbi_txt;
    open (NCBI, Try_GZ_Read($NCBI_TXT)) || die"fail read $NCBI_TXT: $!\n";
    while(<NCBI>){
        if(/^LOCUS/){
            ($LENGTH) = (/(\d+)\s+bp/);
            $CIRCLE = (/circular/)?1:0;
            $Ref_Name_for_Analysis = $VIRUS_TYPE . (($CIRCLE)?'_Extended':'');
        }
        elsif(/^VERSION/){
          ($VERSION) = (/^VERSION\s+(\S+)/);
          # $Ref_Name_for_Analysis = "virus_${VIRUS_TYPE}_${VERSION}_OL$LENGTH".(($CIRCLE)?'_Extended':'');
          # $Ref_Name_for_Complete = "virus_${VIRUS_TYPE}_${VERSION}_OL$LENGTH";
        }
        elsif(/^\s{5}CDS\s/){ # the prefix blanks count is limited
            # get the region string
            (my ($region_string) = (/CDS\s+(.+)\s+$/)) =~ s/\s//g;
            # region written in mutilple lines
            while($region_string =~ /,$/){
                (my $next_region_string = <NCBI>) =~ s/\s//g;
                $region_string .= $next_region_string;
            }
            # strand
            my $strand;
            if($region_string =~ /complement/){
                $strand = '-';
                $region_string =~ s/complement\(//;
                $region_string =~ s/\)$//;
            }
            else{
                $strand = '+';
            }
            # get the sequential region
            $region_string =~ s/join\(//;
            $region_string =~ s/[\(\)]//g; # no bracket again
            $region_string =~ s/\&gt;/>/g; # translate the html '>'
            $region_string =~ s/\&lt;/</g; # translate the html '<'
            # has '<' or '>', 3'/5'-partial
            if($region_string =~ /[\>\<]/){
                $gene_biotype = 'partial_transcript'; # defined by me
            }
            my @CDS_region = split /,/,$region_string;
            my @Comp_CDS_region;
            my ($Comp_CDS_st,$Comp_CDS_ed) = (-1,-1);
            for my $CDS_region (@CDS_region){
                my ($CDS_st,$CDS_ed) = ($CDS_region =~ /[\>\<]?(\d+)[\>\<]?\.+[\>\<]?(\d+)[\>\<]?/);
                if(!defined $CDS_st){
                    cluck_and_exit "$CDS_region\n";
                }
                if($CDS_st < $Comp_CDS_ed){ 
                    if(   $CIRCLE
                       && abs($CDS_st-$Comp_CDS_ed) > $LENGTH * 0.5
                    ){# span the virus circular junction point
                        $CDS_st += $LENGTH;
                        $CDS_ed += $LENGTH;
                        my $ex_len = ceil(($CDS_ed - $LENGTH) / 100) * 100;
                        $Extend_len = max($Extend_len, $ex_len); # keep structure of longest circular gene
                        $Extend_len = min($Extend_len, $LENGTH); # whole length control, max: two folds of original length
                    }
                    else{
                        warn "<WARN>\tCDS back overlap\n"
                                  ."\t$CDS_region\n"
                                  ."\t@CDS_region\n";
                    }
                }
                if($CDS_st != $Comp_CDS_ed + 1){ # not sequential
                    if($Comp_CDS_st != -1){ # not initial
                        push @Comp_CDS_region , [$Comp_CDS_st , $Comp_CDS_ed];
                    }
                    $Comp_CDS_st = $CDS_st;
                }
                $Comp_CDS_ed = $CDS_ed;
            }
            push @Comp_CDS_region , [$Comp_CDS_st , $Comp_CDS_ed]; # remained
            # get the tag information
            my ($locus_tag, $gene, $product, $protein_id, $note);
            while(<NCBI>){
                if(/\/gene=\"(.+)\"/){
                    ($gene = $1) =~ s/\s+/-/g;
                }
                elsif(/\/locus_tag=\"(.+)\"/){
                    ($locus_tag = $1) =~ s/\s+/-/g;
                }
                elsif(/\/product=\"(.+)\"/){
                    ($product = $1) =~ s/\s+/-/g;
                }
                elsif(/\/protein_id=\"(.+)\"/){
                    ($protein_id = $1) =~ s/\s+/-/g;
                }
                elsif(/\/note=\"(.+)\"/){
                    ($note = $1) =~ s/\s+/-/g;
                }
                elsif(/^\s+\/translation=\"/){
                    last;
                }
            }
            # record gene as any one with priority from left to right
            my $gene_record = $gene || $locus_tag || $note || $product || $protein_id;
            $gene_record =~ s/[\/\-\^]/_/g;
            $gene_record = $VIRUS_TYPE.'_'.$gene_record if $gene_record !~ /^$VIRUS_TYPE/;
            $gene_NO++ unless(exists $gene_record{$gene_record});
            # gene info
            my $gene_Info =  "gene_id \"${VIRUS_TYPE}G$gene_NO\"; "
                            ."gene_name \"$gene_record\";";
            # GTF, gene line
            unless(exists $gene_record{$gene_record}){
                push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
                push @GTF_OutPut , "$gene_biotype\t";
                push @GTF_OutPut , "gene\t$Comp_CDS_region[0][0]\t$Comp_CDS_region[-1][-1]\t.\t$strand\t.\t";
                push @GTF_OutPut , "$gene_Info gene_biotype \"$gene_biotype\";\n"; # gene info has gene_biotype (copy from Ensembl)
            }
            # trans
            my $trans_num = ++$gene_record{$gene_record};
            my $transcript_name = $gene_record.'-'.('0' x (3-length($trans_num))).$trans_num;
            # trans info
            my $trans_Info =  "transcript_id \"${VIRUS_TYPE}T$trans_NO\"; "
                             ."gene_biotype \"$gene_biotype\"; "
                             ."transcript_name \"$transcript_name\";";
            # GTF, trans line
            push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
            push @GTF_OutPut , "$gene_biotype\t";
            push @GTF_OutPut , "transcript\t$Comp_CDS_region[0][0]\t$Comp_CDS_region[-1][-1]\t.\t$strand\t.\t";
            push @GTF_OutPut , "$gene_Info $trans_Info\n";
            # protein info
            my $protein_info =  "protein_id \"${VIRUS_TYPE}P$protein_NO\"; "
                               ."protein_code \"$protein_id\"";
            # meta data, put it at the tail
            $_ ||= 'N/A' for ($locus_tag, $gene, $product, $protein_id, $note);
            my $meta_Info = "locus_tag \"$locus_tag\"; "
                            ."note \"$note\"; "
                            ."product \"$product\";";
            # output gtf info
            for my $exon_NO (1 .. scalar(@Comp_CDS_region)){
                my $i = ($strand eq '-')?(scalar(@Comp_CDS_region)-$exon_NO):($exon_NO-1);
                my ($CDS_st,$CDS_ed) = @{$Comp_CDS_region[$i]};
                # exon info
                my $exon_Info =  "exon_number \"$exon_NO\"; "
                                ."exon_id \"${VIRUS_TYPE}E$exon_wNO\";";
                # GTF, exon line
                push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
                push @GTF_OutPut , "$gene_biotype\t";
                push @GTF_OutPut , "exon\t$CDS_st\t$CDS_ed\t.\t$strand\t.\t";
                push @GTF_OutPut , "$gene_Info $trans_Info $exon_Info $meta_Info\n";
                # GTF, CDS line
                push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
                push @GTF_OutPut , "$gene_biotype\t";
                push @GTF_OutPut , "CDS\t$CDS_st\t$CDS_ed\t.\t$strand\t.\t";
                push @GTF_OutPut , "$gene_Info $trans_Info $exon_Info $protein_info $meta_Info\n";
                # GTF, start codon line
                if($exon_NO == 1){
                    my $codon_st = ($strand eq '+') ? $CDS_st : ($CDS_ed-2);
                    my $codon_ed = ($strand eq '+') ? ($CDS_st+2) : $CDS_ed;
                    push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
                    push @GTF_OutPut , "$gene_biotype\t";
                    push @GTF_OutPut , "start_codon\t$codon_st\t$codon_ed\t.\t$strand\t.\t";
                    push @GTF_OutPut , "$gene_Info $trans_Info $exon_Info $meta_Info\n";
                }
                # GTF, stop codon line
                if($exon_NO == scalar(@Comp_CDS_region)){
                    my $codon_st = ($strand eq '-') ? $CDS_st : ($CDS_ed-2);
                    my $codon_ed = ($strand eq '-') ? ($CDS_st+2) : $CDS_ed;
                    push @GTF_OutPut , "$Ref_Name_for_Analysis\t";
                    push @GTF_OutPut , "$gene_biotype\t";
                    push @GTF_OutPut , "stop_codon\t$codon_st\t$codon_ed\t.\t$strand\t.\t";
                    push @GTF_OutPut , "$gene_Info $trans_Info $exon_Info $meta_Info\n";
                }
                # global count
                $exon_wNO ++;
            }
            $trans_NO++;
            $protein_NO++;
        }
        elsif(/^ORIGIN/){
            while(<NCBI>){
                last if(!/^\s+\d+/);
                s/[\d\s]//g;
                $Genome_Seq .= $_;
            }
        }
    }
    close NCBI;

    # check
    cluck_and_exit "<ERROR>\t$VIRUS_TYPE genome length is not $LENGTH.\n" unless length($Genome_Seq) == $LENGTH;

    $V_Href->{virusRefSeg}->{$VIRUS_TYPE}->set_circular(circular => $CIRCLE);
    $V_Href->{virusRefSeg}->{$VIRUS_TYPE}->set_extLen(extLen => $Extend_len * $CIRCLE);
    $V_Href->{virusRefSeg}->{$VIRUS_TYPE}->set_length(length => $LENGTH); # original virus genome length
    # modify the $Ref_Name_for_Analysis
    $Ref_Name_for_Analysis .= ($CIRCLE) ? $Extend_len : '';

    # Complete Genome FA
    write_fasta_file(
        SeqSref => \$Genome_Seq,
        FaFile  => $COMP_fa,
        SegName => $Ref_Name_for_Complete,
        LineBase => $V_Href->{line_base},
        CircleExtLen => 0
    );

    # extended genome fa for analysis
    ## extended for circular virus
    write_fasta_file(
        SeqSref => \$Genome_Seq,
        FaFile  => $ANAY_fa,
        SegName => $Ref_Name_for_Analysis,
        LineBase => $V_Href->{line_base},
        CircleExtLen => $Extend_len * $CIRCLE
    );

    # change the Ref_Name_for_Analysis name, as it has been changed
    for my $i (1 .. scalar(@GTF_OutPut)){
        if($i % 4 == 1){
            $GTF_OutPut[$i-1] = "$Ref_Name_for_Analysis\t";
        }
    }

    # GTF file
    open (GTF, Try_GZ_Write($gtfPath)) || die"fail write gtf: $!\n";
    print GTF "#!genome-build-accession NCBI:$VERSION\n";
    print GTF join('',@GTF_OutPut);
    close GTF;

    # generate PSL
    ## gtf object
    my $gtf = BioFuse::BioInfo::Objects::GeneAnno::GTF_OB->new(filePath => $gtfPath);
    ## load GTF
    $gtf->load_GTF;
    ## gene PSL
    $gtf->create_gene_PSL(gPSLpath => $gPSLpath);
    ## trans PSL
    $gtf->create_trans_PSL(tPSLpath => $tPSLpath, Start_codon => $V_Href->{Start_codon}, whole_genome => $ANAY_fa);

    # inform
    stout_and_sterr "[INFO]\tanalysis virus $VIRUS_TYPE NCBI data, get orig/extend fasta, and gene GTF/PSL.\n";
}

#--- write virus virus_info file ---
sub write_virus_info{
    my $virus_info = GetPath(filekey => 'virus_info_doc');
    open (VRFRP,Try_GZ_Write($virus_info)) || die "fail $virus_info: $!\n";
    my @theme_item = ('Subtype', 'GeneBank_ID', 'Length', 'Circular', 'Extend_Length');
    print VRFRP '#' . join("\t", @theme_item ) . "\n"; # theme
    for my $virus_subtype (sort keys %{$V_Href->{virusRefSeg}}){
        my $virus_OB = $V_Href->{virusRefSeg}->{$virus_subtype};
        print VRFRP join("\t", $virus_subtype,
                               $virus_OB->genebank,
                               $virus_OB->length,
                               $virus_OB->is_circular ? 'yes' : 'no',
                               $virus_OB->extLen
                        )."\n";
    }
    close VRFRP;
    # attach 'readme'
    my $readme = GetPath(filekey => 'virus_db_readme');
    my $source = File::Spec->catfile($RealBin, 'Virus_Type', (split /::/,$MODULE_NAME)[-1].".readme.txt");
    `cp -f $source $readme`;
    # inform
    stout_and_sterr "[INFO]\twrite virus database info.\n"
                         ."\t$virus_info\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
