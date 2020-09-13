package FuseSV::Check;

use strict;
use warnings;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              check
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'Check';
#----- version --------
$VERSION = "0.51";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        check
                        check_Find_MH_setting
                        check_LMcircos_setting
                        check_draw_genome_para_format
                        check_virus_depth_para_format
                        inform_Allow_vNO
                        check_java_version
                        check_virus_db
                     /;

#--- check para and database
sub check{

    my $second_cmd = $V_Href->{command};

    # version check, this part consumes lots of vmem: use '-Xmx1m', ~1.5G; OR >9G.
    if($second_cmd eq 'virus_type'){
        &check_java_version;
    }

    # virus database check
    if(    $second_cmd eq 'virus_type'
        || $second_cmd eq 'extend_vref'
        || $second_cmd eq 'mut_virus'
        || $second_cmd eq 'virus_tree'
        || $second_cmd eq 'draw_virus'
      ){
        &check_virus_db;
    }

    if(    ( $second_cmd eq 'virus_type'
        &&   $V_Href->{stepToStop} >= 4 )
        ||   $second_cmd eq 'mut_link'
      ){
        &inform_Allow_vNO;
    }

    # simple bam for integration
    if(    $second_cmd eq 'virus_type'
        && $V_Href->{stepToStop} >= 4
        && $V_Href->{Get_simple_BamForInteg}
      ){
        # inform
        stout_and_sterr "[INFO]\twill generate simple Bam for integration.\n";
    }

    # maximum amount of virus variants to reconstruct
    if(    ( $second_cmd eq 'virus_type'
        &&   $V_Href->{stepToStop} >= 3 )
      ){
        # inform
        stout_and_sterr "[INFO]\twill reconstruct up to $V_Href->{max_vNO} virus variant(s) if possible.\n";
    }

    # # modification of the tissue_tag
    # if(    $second_cmd eq 'virus_type'
    #     || $second_cmd eq 'cand_reads'
    #   ){
    #     &reset_tissue_tag;
    # }

    # virus depth
    if( $second_cmd eq 'virus_depth' ){
        &check_virus_depth_para_format;
    }

    # draw virus
    if(    $second_cmd eq 'draw_virus'
        || $second_cmd eq 'draw_segCN'
    ){
        &check_draw_genome_para_format;
    }

    # local map
    if( $second_cmd eq 'LMcircos' ){
        &check_LMcircos_setting;
    }

    # search MH
    if( $second_cmd eq 'find_MH' ){
        &check_Find_MH_setting;
    }
}

#--- check setting of MH searching ---
sub check_Find_MH_setting{
    # check sfMH's min_overlap and max_len
    if( $V_Href->{min_sfMHoplen} > $V_Href->{max_sfMH} ){
        warn_and_exit "<ERROR>\t'-minsfo' is larger than '-max_sfMH', can not find any shifting MH.\n";
    }
    # check sfMH's min_overlap and min_len
    if( $V_Href->{min_sfMHoplen} > $V_Href->{min_sfMH} ){
        warn "<WARN>\t'-minsfo' is larger than '-minsfm', reset '-minsfm' as '-minsfo'.\n";
    }
}

#--- check settings of local map circos and alert ---
sub check_LMcircos_setting{
    # check the centres
    ( $V_Href->{cx}, $V_Href->{cy} ) = ( $V_Href->{center} =~ /(\d*\.?\d*),(\d*\.?\d*)/ );
    if( !$V_Href->{cx} || !$V_Href->{cy} ){
        warn_and_exit "<ERROR>\twrong input centre position.\n";
    }

    # check the specific radius
    my @spec_radius;
    for ( 1 .. scalar(@{$V_Href->{radius_specific}}) ){
        my $spec_radius_str = shift @{$V_Href->{radius_specific}};
        my ($refseg, $st_pos, $ed_pos, $spec_radius) = ($spec_radius_str =~ /^([^:]+):(\d+)-(\d+):(\d+)$/);
        unless( $refseg && $st_pos && $ed_pos && $spec_radius){
            warn_and_exit "<ERROR>\tyour specific radius inputs is abnormal: $spec_radius_str\n";
        }
        # quick
        next if( $spec_radius == $V_Href->{host_radius} );
        # for next value-check
        push @spec_radius, $spec_radius;
        # restore
        push @{$V_Href->{radius_specific}}, [ $refseg, $st_pos, $ed_pos, $spec_radius ];
    }

    # check the radius availability
    if(    $V_Href->{cx} < 1.25 * max($V_Href->{host_radius}, $V_Href->{virus_radius}, @spec_radius)
        || $V_Href->{cy} < 1.25 * max($V_Href->{host_radius}, $V_Href->{virus_radius}, @spec_radius)
      ){
        my $aim_radius = int( min( $V_Href->{cx}, $V_Href->{cy} ) / 1.25 );
        warn_and_exit "<ERROR>\tplease reset radius (-hrds, -vrds and -srds) less than $aim_radius.\n";
    }

    # check the PNC for degenerated path
    for my $PNC_deg_char ( @{$V_Href->{PNC_deg_char}} ){
        if( length($PNC_deg_char) != 1 ){
            warn_and_exit "<ERROR>\tplease use single char for PNC for degenerated path. (-dgch)\n";
        }
    }

    # check the gene list allowed to show
    if( $V_Href->{allow_gene_str} ){
        $V_Href->{allow_gene_Href}->{uc($_)} = 1 for split /,/,$V_Href->{allow_gene_str};
    }

    # check the gene name anchor position
    if( $V_Href->{genename_anchor_gTail} || $V_Href->{genename_anchor_gHead} ){
        $V_Href->{genename_anchor_Href}->{gTail}->{uc($_)} = 1 for split /,/,$V_Href->{genename_anchor_gTail};
        $V_Href->{genename_anchor_Href}->{gHead}->{uc($_)} = 1 for split /,/,$V_Href->{genename_anchor_gHead};
    }

    # segments' sorting
    my %show_seg;
    for my $seg_reg ( split /,/, $V_Href->{segment_sort_str} ){
        if( $seg_reg =~ /^(\d+)-(\d+)$/ ){
            for ($1 .. $2){
                warn_and_exit "<ERROR>\tsegment id $_ has been stated before (-sts).\n" if(exists($show_seg{$_}));
                push @{$V_Href->{segment_sort_Aref}}, $_;
                $show_seg{$_} = 1;
            }
        }
        elsif( $seg_reg =~ /^(\d+|[XY])$/ ){
            warn_and_exit "<ERROR>\tsegment id $1 has been stated before (-sts).\n" if(exists($show_seg{$1}));
            push @{$V_Href->{segment_sort_Aref}}, $1;
            $show_seg{$1} = 1;
        }
        else{
            warn_and_exit "<ERROR>\tcannot recognize segment-id ($seg_reg) from '-sts' inputs.\n";
        }
    }

    # log scale
    if( $V_Href->{dna_cn_log} < 0 || $V_Href->{exp_log} < 0 ){
        warn_and_exit "<ERROR>\tlog scale of dna_copy_number or expression is negative.\n";
    }
}

#--- check the para format of draw virus cmd ---
sub check_draw_genome_para_format{
    #--- check normalize depth of Normal tissue
    if( $V_Href->{dVG_case_SeqReadsCount} && $V_Href->{dVG_ctrl_SeqReadsCount} ){
        # ratio to normalize depth
        $V_Href->{dVG_ctrl2case_depthRatio} = $V_Href->{dVG_case_SeqReadsCount} / $V_Href->{dVG_ctrl_SeqReadsCount};
        # inform
        stout_and_sterr "[INFO]\twill normalize depth of Normal tissue with ratio " . sprintf("%.3f",$V_Href->{dVG_ctrl2case_depthRatio}) . "\n";
    }
    elsif( $V_Href->{dVG_case_DepFile} && $V_Href->{dVG_ctrl_DepFile} ){
        # warn
        stout_and_sterr "<WARN>\tno normalization will operate on depth of Normal tissue as Options '-tsrc' or '-nsrc' is unset.\n"
                                       ."\tYou may get un-equivalent depth comparsion of Tumor and Normal.\n";
    }
    #--- check color duplication
    if(    $V_Href->{dVG_case_DepFile} && $V_Href->{dVG_ctrl_DepFile}
        && $V_Href->{dG_DepthSpecColor}->{case} eq $V_Href->{dG_DepthSpecColor}->{ctrl}
    ){
        warn_and_exit "<ERROR>\tcannot use same color ($V_Href->{dG_DepthSpecColor}->{case}) for Tumor and Normal depth spectrum.\n";
    }
    #--- draw options
    $V_Href->{dG_extWidth}  = max($V_Href->{dG_extWidth}, 150);
    $V_Href->{dG_extHeight} = max($V_Href->{dG_extHeight}, 150);
}

#--- reset tissue_tag if possible ---
# sub reset_tissue_tag{
#     if( defined $V_Href->{tissue_tag} ){
#         $V_Href->{tissue_tag} = '_'.$V_Href->{Tool_Tag}.'ts'.$V_Href->{tissue_tag}.$V_Href->{Tool_Tag};
#     }
#     else{
#         $V_Href->{tissue_tag} = '';
#     }
# }

#--- check the para format of virus depth cmd ---
sub check_virus_depth_para_format{
    #--- check region Format for viral depth
    if(    $V_Href->{viral_depth_region} !~ /^([^:]+)$/ # chr
        && $V_Href->{viral_depth_region} !~ /^([^:]+):(\d+)-(\d+)$/ # chr:st-ed
    ){
        # warn
        warn_and_exit "<ERROR>\tthe viral region to stat depth is in wrong Format.\n"
                                     . "\tYour input: $V_Href->{viral_depth_region}\n"
                                     . "\tStd Format: 'vseg' or 'vseg:from-to'\n";
    }
    else{
        # get the viral seg name
        ($V_Href->{viral_depth_vseg}) = ($V_Href->{viral_depth_region} =~ /^([^:]+)/);
    }
    #--- check the minimum match length of read
    if( $V_Href->{viral_depth_MinAllowMatchLen} < 30 ){
        $V_Href->{viral_depth_MinAllowMatchLen} = 30;
        # warn
        stout_and_sterr "<WARN>\tReset minimum match length of read to 30 nt.\n";
    }
}

#--- inform vNO allowed for analysis ---
sub inform_Allow_vNO{
    my $second_cmd = $V_Href->{command};
    #--- different context
    my ($container, $info_phrase);
    if($second_cmd eq 'virus_type'){
        $container = 'Prepare_For_Integ';
        $info_phrase = 'preparation for integration detection';
    }
    elsif( $second_cmd eq 'mut_link' ){
        $container = 'Variant_Mut_Linkage';
        $info_phrase = 'virus mutation linkage analysis';
    }
    #--- selection of virus variants, if set
    ## if use 'all', no 'AllowVNO' will assign to allow all virus
    my $inform_vNOs = 'ALL';
    if( $V_Href->{VariantNO_For_Analysis} !~ /\ball\b/i ){
        $V_Href->{$container}->{AllowVNO}->{$_} = 1 for grep /^\d+$/ && $_ > 0, (split /,/, $V_Href->{VariantNO_For_Analysis});
        $inform_vNOs = 'NO.' . join(', NO.', sort {$a<=>$b} keys %{$V_Href->{$container}->{AllowVNO}});
    }
    else{ # pretend to accept all: No sample could reach 100 virus variants
        $V_Href->{$container}->{AllowVNO}->{$_} = 1 for ( 1 .. 100 );
    }
    # inform
    stout_and_sterr "[INFO]\twill select $inform_vNOs virus variants in $info_phrase.\n";
}

#--- check java version ---
sub check_java_version{
    my ($java_version) = (`$V_Href->{java} -Xmx1m -version 2>&1` =~ /\s+version\D+(\d+\.\d+)/);
    if($java_version < 1.8){
        warn_and_exit "<ERROR>\tthe version of java must be at least 1.8, but yours is $java_version\n";
    }
}

#--- test para and alert ---
sub check_virus_db{
    # virus genomes
    $V_Href->{virus_orig_genome_fa} = GetPath( filekey => 'virus_orig_genome_fa' );
    $V_Href->{virus_info_doc} = GetPath( filekey => 'virus_info_doc' );
    my $virus_orig_genome_fai = GetPath( filekey => 'virus_orig_genome_fai' );
    my $virus_orig_genome_bwt = GetPath( filekey => 'virus_orig_genome_bwt' );
    if( !-e $V_Href->{virus_orig_genome_fa} || !-e $V_Href->{virus_info_doc} ){
        warn_and_exit "<ERROR>\tthe virus database folder must contain virus.genome.fa and virus.info under 'virus_original_reference'.\n";
    }
    if( !-e $virus_orig_genome_fai || !-e $virus_orig_genome_bwt ){
        warn_and_exit "<ERROR>\tvirus.genome.fa in virus database folder hasn't been indexed by BWA.\n";
    }

    # check virus annotation files
    $V_Href->{virus_orig_gene_psl}  = GetPath( filekey => 'virus_orig_gene_psl' );
    $V_Href->{virus_orig_trans_psl} = GetPath( filekey => 'virus_orig_trans_psl' );
    if( !-e $V_Href->{virus_orig_gene_psl} || !-e $V_Href->{virus_orig_trans_psl} ){
        warn_and_exit "<ERROR>\tthe virus database folder lacks PSL file under 'virus_gene_annotation' folder.\n";
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
