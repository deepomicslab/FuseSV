package FuseSV::GetPath;

use strict;
use warnings;
use File::Spec qw/ catfile /;
use BioFuse::Util::Log qw/ stout_and_sterr cluck_and_exit /;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GetPath
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::GetPath';
#----- version --------
$VERSION = "0.05";
$DATE = '2019-06-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--- get path ---
sub GetPath{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $filekey = $parm{filekey};
    my $virus_subtype = $parm{virus_subtype} || undef; # virus name, such as HPV18
    my $CandFq_RG_ID_NO = $parm{CandFq_RG_ID_NO} || undef;
    my $refseg_NO = $parm{refseg_NO} || undef;
    my $basename = $parm{basename} || undef;
    my $vref_extLen = $parm{vref_extLen} || undef;
    my $virus_variant_NO = $parm{virus_variant_NO} || undef;

    # func json
    if( $filekey eq 'func_json' ){
        return File::Spec->catfile($V_Href->{RealBin},'functions.json');
    }

    # virus db
    elsif( $filekey eq 'virus_db_readme' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'readme.txt');
    }
    elsif( $filekey eq 'NCBI_HtmlText_folder' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'NCBI_page');
    }
    elsif( $filekey eq 'virus_NCBI_page_txt' ){
        return File::Spec->catfile(GetPath(filekey => 'NCBI_HtmlText_folder'), $virus_subtype.'.ncbi_page.txt.gz');
    }
    elsif( $filekey eq 'virus_orig_ref_folder' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'virus_original_reference');
    }
    elsif( $filekey eq 'virus_orig_genome_fa' ){
        return File::Spec->catfile(GetPath(filekey => 'virus_orig_ref_folder'), 'virus.genome.fa');
    }
    elsif( $filekey eq 'virus_orig_genome_fai' ){
        return $V_Href->{virus_orig_genome_fa}.'.fai';
    }
    elsif( $filekey eq 'virus_orig_genome_bwt' ){
        return $V_Href->{virus_orig_genome_fa}.'.bwt';
    }
    elsif( $filekey eq 'virus_info_doc' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'virus.info');
    }
    elsif( $filekey eq 'virus_gene_anno_folder' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'virus_gene_annotation');
    }
    elsif( $filekey eq 'virus_orig_gene_psl' ){
        return File::Spec->catfile(GetPath(filekey => 'virus_gene_anno_folder'), 'virus.extended.gene.psl');
    }
    elsif( $filekey eq 'virus_orig_trans_psl' ){
        return File::Spec->catfile(GetPath(filekey => 'virus_gene_anno_folder'), 'virus.extended.trans.psl');
    }
    elsif( $filekey eq 'single_virus_data_folder' ){
        return File::Spec->catfile($V_Href->{virus_db_folder}, 'single_virus_data');
    }
    elsif( $filekey eq 'single_viral_orig_tpsl' ){
        return File::Spec->catfile(GetPath(filekey => 'single_virus_data_folder'), $virus_subtype.'.extended.trans.psl');
    }

    # prepare extended virus db
    elsif( $filekey eq 'extend_vref_db_folder' ){
        return File::Spec->catfile( $V_Href->{virus_db_folder}, 'virus_divided_IndexRef', 'ExtendLen_'.$vref_extLen );
    }

    # simulate mutated virus
    elsif( $filekey eq 'simu_virus_mutation_list' ){
        return File::Spec->catfile( $V_Href->{mut_virus_outdir}, $virus_subtype.'.mutation.list' );
    }
    elsif( $filekey eq 'simu_virus_mutated_fa' ){
        return File::Spec->catfile( $V_Href->{mut_virus_outdir}, $virus_subtype.'.mutated.fa' );
    }

    # get candidate reads
    elsif( $filekey eq 'CandReadsFolder' ){
        return File::Spec->catfile( $V_Href->{virusAnaFolder}, $V_Href->{sample_ID}.'.CandReads' );
    }
    elsif( $filekey eq 'Cand_RgInfoReport' ){
        return File::Spec->catfile( $V_Href->{CandReadsFolder}, 'RG.info' );
    }
    elsif( $filekey eq 'Cand_smartFqGz' ){
        return File::Spec->catfile( $V_Href->{CandReadsFolder}, 'RG_NO'.$CandFq_RG_ID_NO.'.candidate.smart.fastq.gz' );
    }
    elsif( $filekey eq 'Cand_nSortBam' ){
        return File::Spec->catfile( $V_Href->{CandReadsFolder}, 'RG_NO'.$CandFq_RG_ID_NO.'.candidate.nSort.bam' );
    }
    elsif( $filekey eq 'Tmp_nSortBam' ){
        return File::Spec->catfile( $V_Href->{CandReadsFolder}, '_'.sprintf("%0*d",5,int(rand(10000))).'_'.$basename );
    }

    # initial alignment on each virus ref
    elsif( $filekey eq 'virusInitAlignFolder' ){
        return File::Spec->catfile( $V_Href->{virusAnaFolder}, $V_Href->{sample_ID}.'.VirusInitAlign' );
    }
    elsif( $filekey eq 'virusInitAlignReport' ){
        return File::Spec->catfile( $V_Href->{virusInitAlignFolder}, 'virus.initAlign.report' );
    }
    elsif( $filekey eq 'virusInitAlignBamFolder' ){
        return File::Spec->catfile( $V_Href->{virusInitAlignFolder}, 'InitAlignBAM' );
    }
    elsif( $filekey eq 'virusInitAlignRefFolder' ){
        return File::Spec->catfile( $V_Href->{virusInitAlignFolder}, 'InitAlignREF' );
    }

    # virus variant
    elsif( $filekey eq 'Virus_variant_folder' ){
        return File::Spec->catfile( $V_Href->{virusAnaFolder}, $V_Href->{sample_ID}.'.VirusVariants' );
    }
    elsif( $filekey eq 'Virus_variant_report' ){
        return File::Spec->catfile( $V_Href->{Virus_variant_folder}, $V_Href->{sample_ID}.'.VirusVariants.report' );
    }
    elsif( $filekey eq 'single_virus_variant_folder' ){
        return File::Spec->catfile( $V_Href->{Virus_variant_folder}, 'V'.$virus_variant_NO.'_reconstruction' );
    }

    # prepare for virus integration
    elsif( $filekey eq 'Prepare_Integ_folder' ){
        return File::Spec->catfile( $V_Href->{virusAnaFolder}, $V_Href->{sample_ID}.'.Prepare_DetectInteg' );
    }

    # find MH of virus integrations
    elsif( $filekey eq 'FindMH_report' ){
        return File::Spec->catfile( $V_Href->{output_dir}, $V_Href->{sample_ID}.'.FlankMH_report.txt' );
    }
    elsif( $filekey eq 'FindMH_SVG_figure' ){
        return File::Spec->catfile( $V_Href->{output_dir}, $V_Href->{sample_ID}.'.FlankMH.svg' );
    }

    # none
    else{
        cluck_and_exit "<ERROR>\tunrecognized filekey $filekey in FuseSV::GetPath.\n";
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
