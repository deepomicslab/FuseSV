package FuseSV::Virus_Type::VirusMutLinkage;

use strict;
use warnings;
use File::Spec qw/catfile/;
use File::Basename qw/basename dirname/;
use Getopt::Long;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use FuseSV::LoadOn;
use FuseSV::Virus_Type::VirusTypeMain;
use FuseSV::GetPath qw/ GetPath /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              find_virus_mutation_linkage
            /;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'VirusMutLinkage';
#----- version --------
$VERSION = "0.09";
$DATE = '2018-01-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        find_virus_mutation_linkage
                        load_mutation_of_virus_variant
                        load_variants_modifications
                        get_read_mut_chain_of_virus_variant
                        load_variants_ReadMutChain_finalBam
                     /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
     Usage:   perl $V_Href->{MainName} mut_link <[Options]>

     Options:

       # Inputs and Outputs #
       # Note: options ('-id' and '-odir') should be as \033[4m \bsame \b\033[0m as you used in 'virus_type' #
        -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
                      only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
        -odir   [s]  directory to store outputs, must exist already. <required>

       # Software Required #
        -stl    [s]  SamTools, minimum version: 1.3. <required>

        # Mutation Linkage Options #
        -vno    [s]  to select virus variants (its NO.) for virus mutation linkage analysis. [1]
                      default value is 1, which means only use the FIRST virus variant;
                      set as 'all' to use ALL virus variants;
                      set as '1,2' to use FIRST and SECOND variants (if it has, automatically judgement).
        -minsr  [i]  the minimum amount of supporting reads to accept one allele, sum of the both strands. [2]
        -bstrd  [i]  set non-zero to require allele supported by at least INT reads both strands, repectively. [1]
                      a) set as 0 means supported by reads from only one map-strand will also be good for analysis.
                      b) note that, generally, we assign non-zero value; but for virus with linear genome (e.g., HPV59),
                         we'd better set it zero, because specially for the mutations located at the bilateral ends of
                         viral genome, which only provide either of reads-mapping orientations (e.g., 'F' or 'R').

        -h|help   Display this help info.

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            [ outdir => '' ],
            [ sample_ID => '' ],
            # software
            [ samtools => '' ],
            # intermediate file
            [ Virus_variant_report => '' ],
            # draw options
            [ VariantNO_For_Analysis => 1 ],
            [ both_strd => 1 ],
            [ min_mut_reads => 2 ],
            # container
            ## -> AllowVNO
            ## -> mutations
            [ Variant_Mut_Linkage => {} ],

            # list to abs-path
            [ ToAbsPath_Aref => [ ['outdir'],
                                  ['samtools'] ] ]
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
        "-odir:s"   => \$V_Href->{outdir},
        # software
        "-stl:s"    => \$V_Href->{samtools},
        # options
        "-vno:s"    => \$V_Href->{VariantNO_For_Analysis},
        "-bstrd:i"  => \$V_Href->{both_strd},
        "-minsr:i"  => \$V_Href->{min_mut_reads},
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
            );
}

#--- find the mutation linkage chain of detected virus variants ---
sub find_virus_mutation_linkage{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # file path
    $V_Href->{Virus_variant_folder} = GetPath(filekey => 'Virus_variant_folder');
    $V_Href->{Virus_variant_report} = GetPath(filekey => 'Virus_variant_report');

    # read virus variant report
    FuseSV::Virus_Type::ReconstrcutVirusVariants->loadon_virus_variant_report;

    # read modification list of virus variants
    &load_mutation_of_virus_variant;

    # read final bam of virus variants
    &get_read_mut_chain_of_virus_variant;
}

#--- read modification list ---
sub load_mutation_of_virus_variant{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Variant_Mut_Linkage_Href = $V_Href->{Variant_Mut_Linkage};

    # deal with each virus variant
    for my $variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
        # only selected virus variant
        unless( exists($Variant_Mut_Linkage_Href->{AllowVNO}->{$variant_NO}) ){
            next;
        }
        # are mutations completed?
        unless( $Variant_Info_Href->{$variant_NO}->{ModCompBool} ){
            stout_and_sterr `date`."<WARN>:\tskip virus variant V$variant_NO as its listed mutations are not completed.\n";
            next;
        }
        # load modification list for mutations
        &load_variants_modifications({
            variant_NO => $variant_NO
        });
        # skip if no mutation for mut_link analysis
        unless( $Variant_Info_Href->{$variant_NO}->{HasMutToLink} ){
            stout_and_sterr `date`."<WARN>:\tskip virus variant V$variant_NO as it has no good mutation for mut_link analysis.\n";
            next;
        }
        # check modifications
        if($V_Href->{in_debug}){
            my $Var_Mut_Href = $Variant_Mut_Linkage_Href->{mutations}->{$variant_NO};
            for my $modvpos (sort {$a<=>$b} keys %$Var_Mut_Href){
                for my $mut_id (sort keys %{$Var_Mut_Href->{$modvpos}}){
                    my $mutInfo_Href = $Var_Mut_Href->{$modvpos}->{$mut_id};
                    my $out_str = "$variant_NO\t$modvpos\t$mut_id";
                    $out_str .= "\tsupport:" . $mutInfo_Href->{support};
                    $out_str .= "\tallele_attr:" . $mutInfo_Href->{allele_attr};
                    $out_str .= "\treal_mut_info:" . $mutInfo_Href->{real_mut_info};
                    $out_str .= "\tsnp_map_extlen:" . $mutInfo_Href->{snp_map_extlen};
                    $out_str .= "\tfindreads:" . $mutInfo_Href->{findreads};
                    warn "$out_str\n";
                }
            }
        }
    }
}

#--- load modification list for mutations of given virus variant ---
sub load_variants_modifications{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_NO = $Option_Href->{variant_NO};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Variant_Mut_Linkage_Href = $V_Href->{Variant_Mut_Linkage};
    my $V2O_PosMap_Href = \%{ $Variant_Mut_Linkage_Href->{v2o_posmap}->{$variant_NO} };
    my $Var_Mut_Href = \%{ $Variant_Mut_Linkage_Href->{mutations}->{$variant_NO} };

    my $this_variant_Href = $Variant_Info_Href->{$variant_NO};
    open (MOD,Try_GZ_Read($this_variant_Href->{Mods})) || die "fail read $this_variant_Href->{Mods}: $!\n";
    # theme
    (my $theme_line = lc(<MOD>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<MOD>){
        next if(/^\#/);
        my @info = split;
        my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);

        # record this position mappings
        $V2O_PosMap_Href->{ $rgOb{mutvrefpos} } = $rgOb{originalpos};

        # is orig-vref-allele good enough?
        my $orig_allele_is_good = &whether_allele_enough_supported( { support_count => $rgOb{ref_supp} } );

        # record all alleles
        # different scenarios
        # Notes of 'allele_attr':
        #  1) 'R' means orig_Vref-allele;
        #  2) 'M' means major_allele;
        #  3) 'm' means minor_allele;
        #     major-Allele-attr is action refer to [orig-vref]
        #  4) 'S' means major-Allele 'Snp';
        #  5) 'I' means major-Allele 'Insertion';
        #  6) 'D' means major-Allele 'Deletion'.
        #     minor-Allele-attr is action refer to [mod-vref]
        #  7) 's' means minor-Allele 'sNP';
        #  8) 'i' means minor-Allele 'iNSERTION';
        #  9) 'd' means minor-Allele 'dELETION'.
        # 10) 'A' means 'Artificial' minor-Allele concatenated by two alleles.
        ## this position in mod-Vref is modified with major_alt
        if( $rgOb{in_modvirus} eq 'Y' ){
            # is mod-vref-allele good enough?
            my $modv_allele_is_good = &whether_allele_enough_supported( { support_count => $rgOb{majoralt_supp} } );
            # as long as one allele is good
            if( $orig_allele_is_good || $modv_allele_is_good ){
                # different mut_type, different records
                my ($orig_Vref_allele_mut_id, $orig_Vref_allele_attr);
                my ($mod_Vref_allele_mut_id,  $mod_Vref_allele_attr);
                my ($mod_Vref_snp_map_extlen);
                ## majorMut is snp
                if( lc($rgOb{majoralt_type}) eq 'snp' ){
                    ## record original Vref allele
                    $orig_Vref_allele_mut_id = 'snp,' . uc($rgOb{ref_allele});
                    $orig_Vref_allele_attr   = 'Rms'; 
                    ## record reconstructed Vref allele
                    $mod_Vref_allele_mut_id = 'snp,' . uc($rgOb{majoralt_allele});
                    $mod_Vref_allele_attr   = 'MS';
                    ## no need to check next base
                    $mod_Vref_snp_map_extlen = 0;
                }
                ## majorMut is del
                elsif( lc($rgOb{majoralt_type}) eq 'del' ){
                    ## record original Vref allele
                    $orig_Vref_allele_mut_id = 'ins,' . uc($rgOb{majoralt_allele});
                    $orig_Vref_allele_attr   = 'Rmi'; 
                    ## record reconstructed Vref allele
                    ## pretend as 'snp' with additional extended mapping
                    $mod_Vref_allele_mut_id = 'snp,' . '-'; # '-' means no need to check this base allele
                    $mod_Vref_allele_attr   = 'MD';
                    ## need to check more 1bp for spanning the del-gap
                    $mod_Vref_snp_map_extlen = 1;
                }
                ## majorMut is ins
                elsif( lc($rgOb{majoralt_type}) eq 'ins' ){
                    ## record original Vref allele
                    $orig_Vref_allele_mut_id = 'del,' . uc($rgOb{majoralt_allele});
                    $orig_Vref_allele_attr   = 'Rmd'; 
                    ## record reconstructed Vref allele
                    ## pretend as 'snp' with additional extended mapping
                    $mod_Vref_allele_mut_id = 'snp,' . '-'; # '-' means no need to check this base allele
                    $mod_Vref_allele_attr   = 'MI';
                    ## spanning the whole ins-seq
                    $mod_Vref_snp_map_extlen = length($rgOb{majoralt_allele});
                }
                # to record alleles' info
                my $MutPos_Href = \%{ $Var_Mut_Href->{ $rgOb{mutvrefpos} } };
                # record orig-Vref allele
                if( $orig_allele_is_good ){
                    $MutPos_Href->{$orig_Vref_allele_mut_id}
                                            = {
                                                support => $rgOb{ref_supp},
                                                allele_attr => $orig_Vref_allele_attr,
                                                real_mut_info => join(',', $rgOb{originalpos}, 'ref', uc($rgOb{ref_allele})),
                                                snp_map_extlen => 0,
                                                findreads => 0
                                              };
                }
                # record mod-Vref allele
                if( $modv_allele_is_good ){
                    $MutPos_Href->{$mod_Vref_allele_mut_id}
                                            = {
                                                support => $rgOb{majoralt_supp},
                                                allele_attr => $mod_Vref_allele_attr,
                                                real_mut_info => join(',', $rgOb{originalpos}, lc($rgOb{majoralt_type}), uc($rgOb{majoralt_allele})),
                                                snp_map_extlen => $mod_Vref_snp_map_extlen,
                                                findreads => 0
                                              };
                }
                # record has_mutation sign
                $this_variant_Href->{HasMutToLink} = 1;
            }
        }
        ## this position in mod-Vref is identical with orig-Vref
        elsif( $orig_allele_is_good ){
            # to record alleles' info
            my $MutPos_Href = \%{ $Var_Mut_Href->{ $rgOb{mutvrefpos} } };
            ## record original Vref allele
            ## pretend as 'SNP'
            my $mut_id = 'snp,' . uc($rgOb{ref_allele});
            $MutPos_Href->{$mut_id} = {
                                        support => $rgOb{ref_supp},
                                        allele_attr => 'RMS',
                                        real_mut_info => join(',', $rgOb{originalpos}, 'ref', uc($rgOb{ref_allele})),
                                        snp_map_extlen => 0, # no need to check next base
                                        findreads => 0
                                      };
            # record has_mutation sign
            $this_variant_Href->{HasMutToLink} = 1;
        }

        # record other alteration allele (minorMut)
        if( $rgOb{'minoraltmodvref(ovref)'} ne '-' ){
            my @otheralt_info = split /;/, $rgOb{'minoraltmodvref(ovref)'};
            for my $alt_info (@otheralt_info){
                my @alt_info = split /,/, $alt_info;
                my $support_count = join(',', @alt_info[3,4,5]);
                # this minorMut is good enough
                if( &whether_allele_enough_supported( { support_count => $support_count } ) ){
                    my ($real_mut_info) = ($alt_info =~ /\(([^\(\)]+)\)$/);
                    my $mut_modVref_pos = $alt_info[0];
                    my $mut_id = lc($alt_info[1]) . ',' . uc($alt_info[2]);
                    # to record alleles' info
                    # must link the mut_id to its right modVref-pos
                    my $MutPos_Href = \%{ $Var_Mut_Href->{$mut_modVref_pos} };
                    # record minorMut info
                    $MutPos_Href->{$mut_id} = {
                                                support => $support_count,
                                                allele_attr => 'm'.lc(substr($alt_info[1],0,1)),
                                                real_mut_info => $real_mut_info,
                                                snp_map_extlen => 0,
                                                findreads => 0
                                              };
                    # record has_mutation sign
                    $this_variant_Href->{HasMutToLink} = 1;
                }
            }
        }
    }
    close MOD;

    # inform
    stout_and_sterr `date`."[INFO]\tload all mutations of virus variant V$variant_NO OK.\n";
}

#--- check whether allele supported by enough reads ---
sub whether_allele_enough_supported{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $support_count = $Option_Href->{support_count};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    my ($sum_supp, $forward_supp, $reverse_supp) = (split /,/, $support_count);
    # defined?
    if ( !defined($sum_supp) || !defined($forward_supp) || !defined($reverse_supp) ){
        warn_and_exit `date`."<ERROR>:\tNot available allele's reads support info: $support_count\n";
    }
    # reads requirements
    if(    $sum_supp >= $V_Href->{min_mut_reads}
        && (   $V_Href->{both_strd} == 0
            || min($forward_supp, $reverse_supp) >= $V_Href->{both_strd}
           )
    ){
        return 1; # good
    }
    else{
        return 0; # bad
    }
}

#--- read final bam of virus variants for read mut chain ---
sub get_read_mut_chain_of_virus_variant{

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Mut_Href = $V_Href->{Variant_Mut_Linkage}->{mutations};

    # deal with each virus variant
    for my $variant_NO (sort {$a<=>$b} keys %$Mut_Href){
        # skip if no mutation for mut_link analysis
        unless( $Variant_Info_Href->{$variant_NO}->{HasMutToLink} ){
            next;
        }
        # load bam for reads linkage of mutations
        &load_variants_ReadMutChain_finalBam({
            variant_NO => $variant_NO
        });
    }
}

#---load bam for reads linkage of mutations ---
sub load_variants_ReadMutChain_finalBam{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $variant_NO = $Option_Href->{variant_NO};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();
    my $Variant_Info_Href = $V_Href->{Variant_Info_Href};
    my $Variant_Mut_Linkage_Href = $V_Href->{Variant_Mut_Linkage};
    my $Var_Mut_Href = $Variant_Mut_Linkage_Href->{mutations}->{$variant_NO};
    my $Var_MutLink_Href = \%{ $Variant_Mut_Linkage_Href->{read_mut_link}->{$variant_NO} };

    my $this_variant_Href = $Variant_Info_Href->{$variant_NO};

    # get mod-Vref un-extended length for circular extended part
    my $Vref_extlen = $this_variant_Href->{Extend_length};
    chomp(my $modVref_length = `sed -n '2,\$p' $this_variant_Href->{Fa} | perl -ne 'chomp; print;' | wc -L`);

    # prepare allele read linkage work folder
    $this_variant_Href->{MutLink_work} = File::Spec->catfile($V_Href->{Virus_variant_folder}, "V$variant_NO\_reconstruction", "V$variant_NO.Allele_Linkage");
    `rm -rf $this_variant_Href->{MutLink_work}`;
    `mkdir -p $this_variant_Href->{MutLink_work}`;

    # prepare allele read linkage output file
    my $AlleleReadsLink_gz = File::Spec->catfile($this_variant_Href->{MutLink_work}, 'Allele_ReadsLink.gz');
    open (ARLK, Try_GZ_Write($AlleleReadsLink_gz)) || die "cannot write $AlleleReadsLink_gz: $!\n";
    # '#' prefixed theme line
    my @theme = qw/PID ModVref NO.1-mModVrefPos NO.1-mOrit NO.1-Allele NO.2-mModVrefPos NO.2-mOrit NO.2-Allele/;
    print ARLK '#' . join("\t", @theme) . "\n";

    # seek InDel and Softclip sites
    my $deal_pe_sum = 0;
    my $finalBam = $Variant_Info_Href->{$variant_NO}->{finalBam};
    open (BAM, "$V_Href->{samtools} view -F 0x400 $finalBam |") || die "cannot read $finalBam: $!\n";
    while (<BAM>){
        chomp;
        my ($PID, $FLAG, $mSeq, $mPos, $CIGAR, $ReadSeq) = (split /\t+/)[0, 1, 2, 3, 5, 9];
        # filter
        if(    $FLAG & 0x400 # duplication
            || $FLAG & 0x100 # secondary alignment
            || $FLAG & 0x800 # supplementary alignment
            || $FLAG & 0x4   # unmapped
            || $CIGAR =~ /[^IDSM\d]/ # have other cigar operators
          ){
            next;
        }
        # specically for pos in extended part of circular virus
        if(    $this_variant_Href->{Cricle} eq 'yes'
            && $mPos > $modVref_length
        ){
            $mPos -= $modVref_length;
        }
        # store basic info
        my $end_NO = ($FLAG & 0x40) ? 1 : 2;
        my $m_Orit = ($FLAG & 0x10) ? '-' : '+';
        $Var_MutLink_Href->{$PID}->{$end_NO}->{mPos} = $mPos;
        $Var_MutLink_Href->{$PID}->{$end_NO}->{mOrit} = $m_Orit;
        # deal mutations from CIGAR
        &analysis_CIGAR_for_readMut_Linkage({
            PID => $PID,
            end_NO => $end_NO,
            mPos => $mPos,
            ReadSeq => $ReadSeq,
            CIGAR => $CIGAR,
            Var_Mut_Href => $Var_Mut_Href,
            Var_MutLink_Href => $Var_MutLink_Href
        });
        # output this pair-end allele supporting situation
        my $next_end_NO = ($FLAG & 0x40) ? 2 : 1;
        if(    $FLAG & 0x8 # next fragment is unmapped
            || exists($Var_MutLink_Href->{$PID}->{$next_end_NO})
        ){
            my @output_info = ($PID, $mSeq);
            my $output_sign = 0;
            for my $end (1, 2){
                # mapped read
                if( exists($Var_MutLink_Href->{$PID}->{$end}) ){
                    my $end_Href = $Var_MutLink_Href->{$PID}->{$end};
                    # has support Allele at Mut site
                    if( exists($end_Href->{supp_Mut}) ){
                        push @output_info, $end_Href->{mPos};
                        push @output_info, $end_Href->{mOrit};
                        push @output_info, join( ';', @{ $end_Href->{supp_Mut} } );
                        $output_sign = 1;
                    }
                    # no support Allele at Mut site
                    else{
                        push @output_info, (('N/A') x 3);
                    }
                }
                # unmapped
                else{
                    push @output_info, (('UMP') x 3);
                }
            }
            # output?
            if( $output_sign ){
                print ARLK join("\t", @output_info) . "\n";
            }
            # sweep
            delete $Var_MutLink_Href->{$PID};
            if( ++$deal_pe_sum % $V_Href->{Handle_MutLinkage_PE_alert_unit} == 0){
                stout_and_sterr `date`."[INFO]\thandled $deal_pe_sum pair-ends from finalBam of virus variant V$variant_NO.\n";
            }
        }
    }
    close BAM;

    close ARLK;

    # inform
    stout_and_sterr `date`."[INFO]\tTotally handled $deal_pe_sum pair-ends from finalBam of virus variant V$variant_NO.\n";
    stout_and_sterr `date`."[INFO]\tload finalBam for allele reads' linkage of virus variant V$variant_NO OK.\n";
    stout_and_sterr `date`."[INFO]\toutput allele reads' linkage of virus variant V$variant_NO OK.\n"
                          ."\t$AlleleReadsLink_gz\n";
}

#--- link pid/[12] to certain mapped orgi_mut ---
sub analysis_CIGAR_for_readMut_Linkage{

    my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
    my $PID = $Option_Href->{PID};
    my $end_NO = $Option_Href->{end_NO};
    my $mPos = $Option_Href->{mPos};
    my $ReadSeq = $Option_Href->{ReadSeq};
    my $CIGAR = $Option_Href->{CIGAR};
    my $Var_Mut_Href = $Option_Href->{Var_Mut_Href};
    my $Var_MutLink_Href = $Option_Href->{Var_MutLink_Href};

    # variants
    my $V_Href = FuseSV::LoadOn->load_variants_dict();

    # check CIGAR 'S'
    if(    $CIGAR =~ /S/
        && (   $CIGAR !~ /^\d+S/ # 5-prime soft-clip
            && $CIGAR !~ /S$/    # 3-prime soft-clip
           )
    ){
        warn_and_exit `date`."<ERROR>:\tencounter abnormal SoftClip CIGAR ($CIGAR) of end $end_NO from PID:\n"
                            ."\t$PID\n";
    }

    # deal mutations from CIGAR
    my $i = 0;
    while(length($CIGAR)!=0){
        my ($len, $type) = ($CIGAR =~ /^(\d+)(\D)/);
        $CIGAR =~ s/^\d+\D//;
        # mapped with possible mismatch ONLY
        if($type eq 'M'){
            my $deal_seq = substr($ReadSeq, 0, $len);
            my @deal_seq = split //, $deal_seq;
            for my $shift (0 .. $len-1){
                my $judge_pos = $mPos + $i + $shift;
                # skip this pos if no mutation at the site
                next unless( exists( $Var_Mut_Href->{$judge_pos} ) );
                my $allele = uc($deal_seq[$shift]);
                # just single base
                my @cand_mut_id = ( "snp,$allele", 'snp,-' );
                for my $mut_id ( @cand_mut_id ){
                    # skip it is not allele we found at this pos before
                    next unless( exists( $Var_Mut_Href->{$judge_pos}->{$mut_id} ) );
                    # check whether needs to extend
                    my $found_mut_id_Href = $Var_Mut_Href->{$judge_pos}->{$mut_id};
                    my $snp_map_extlen = $found_mut_id_Href->{snp_map_extlen};
                    if( $snp_map_extlen != 0 ){
                        my $remain_Mlen = $len - 1 - $shift;
                        if( $remain_Mlen < $snp_map_extlen ){
                            next;
                        }
                    }
                    # count in the mutation read supports
                    $found_mut_id_Href->{findreads} ++;
                    # link this read with the mutation
                    push @{ $Var_MutLink_Href->{$PID}->{$end_NO}->{supp_Mut} }, $found_mut_id_Href->{real_mut_info};
                }
            }
            # update for next block
            $i += $len;
            $ReadSeq = substr($ReadSeq, $len);
        }
        elsif($type eq 'I'){
            my $ins_seq = substr($ReadSeq, 0, $len);
            # assign insertion to the forefront-pos
            my $mut_pos = $mPos + $i - 1;
            my $mut_id = 'ins,' . uc($ins_seq);
            # this site has allele we recorded before
            if(    exists( $Var_Mut_Href->{$mut_pos} )
                && exists( $Var_Mut_Href->{$mut_pos}->{$mut_id} )
            ){
                my $found_mut_id_Href = $Var_Mut_Href->{$mut_pos}->{$mut_id};
                # count in the mutation read supports
                $found_mut_id_Href->{findreads} ++;
                # link this read with the mutation
                push @{ $Var_MutLink_Href->{$PID}->{$end_NO}->{supp_Mut} }, $found_mut_id_Href->{real_mut_info};
            }
            # update for next block
            $ReadSeq = substr($ReadSeq, $len);
        }
        elsif($type eq 'D'){
            my $del_len = $len;
            # assign deletion to the first pos of the deleted part
            my $mut_pos = $mPos + $i;
            # this site has allele we recorded before
            if( exists( $Var_Mut_Href->{$mut_pos} ) ){
                # as del_seq is unknown without mod-Vref, so traverse mut_id at this pos
                # only del, see whether del_len is identical
                for my $mut_id (grep /^del,/, keys %{$Var_Mut_Href->{$mut_pos}}){
                    my ($del_seq) = (split /,/, $mut_id)[1];
                    # skip, if deleted length is different
                    next if( length($del_seq) != $del_len );
                    # here, this deletion is recorded minorMut-del
                    my $found_mut_id_Href = $Var_Mut_Href->{$mut_pos}->{$mut_id};
                    # next, we need to check whether its front pos is recorded with majorMut-del
                    # if it is, connect them two together as one DEL mutaton
                    my $found_majorDel_sign = 0;
                    my $last_pos = $mut_pos - 1;
                    if( exists($Var_Mut_Href->{$last_pos}) ){
                        for my $lastPos_mut_id (sort keys %{$Var_Mut_Href->{$last_pos}}){
                            my $lastPos_mut_id_Href = $Var_Mut_Href->{$last_pos}->{$lastPos_mut_id};
                            my $allele_attr = $lastPos_mut_id_Href->{allele_attr};
                            # it is majorMut-del
                            if( $allele_attr =~ /M/ && $allele_attr =~ /D/ ){
                                my @majorDel_info = split /,/, $lastPos_mut_id_Href->{real_mut_info};
                                my @minorDel_info = split /,/, $found_mut_id_Href->{real_mut_info};
                                my $real_mut_info = join( ',', $majorDel_info[0], 'del', uc($majorDel_info[2].$minorDel_info[2]) );
                                # here, we cannot create this artificial mutation Href, neither the 'findreads' accumulation count.
                                # link this read with the mutation
                                push @{ $Var_MutLink_Href->{$PID}->{$end_NO}->{supp_Mut} }, $real_mut_info;
                                # update connection status
                                $found_majorDel_sign = 1;
                                # inform
                                if( $V_Href->{in_debug} ){
                                    warn `date` . "[INFO]:\tConcatenate MinorMut-Del ($found_mut_id_Href->{real_mut_info}) with MajorMut-Del ($lastPos_mut_id_Href->{real_mut_info}) as one Del ($real_mut_info) on Read:\n"
                                                . "\t$PID/$end_NO\n";
                                }
                                # no more work do next
                                last;
                            }
                        }
                    }
                    # the front pos is not recorded with majorMut-del
                    if( $found_majorDel_sign == 0 ){
                        # count in the mutation read supports
                        $found_mut_id_Href->{findreads} ++;
                        # link this read with the mutation
                        push @{ $Var_MutLink_Href->{$PID}->{$end_NO}->{supp_Mut} }, $found_mut_id_Href->{real_mut_info};
                    }
                }
            }
            # update for next block
            $i += $len;
        }
        elsif($type eq 'S'){
            # just skip S part, update for next block
            $ReadSeq = substr($ReadSeq, $len);
        }
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
