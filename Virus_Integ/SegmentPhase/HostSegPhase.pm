package FuseSV::Virus_Integ::SegmentPhase::HostSegPhase;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use Data::Dumper;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success file_exist /;
use BioFuse::Dist::DistStat qw/ get_value_mean /;
use FuseSV::LoadOn;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use FuseSV::Objects::SNV_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  phase_host_segment
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'HostSegPhase';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						phase_host_segment
						prepare
						read_segment_list
						extract_segmVCF_from_SNP_source
						extract_segmVCF_from_gmBAM
						extract_segmVCF_from_HRC
						select_ImbHet_sourceSNP_in_csaeBam
						load_HRC_sample_list
						load_hetImbSNP_from_HRC
						extract_het_sourceSNP_from_caseBAM
						calculate_segm_mVSm_HRCphaseScore
						calculate_segm_mVsm_LongRangeLinkScore
						output_report
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} SegPhase <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -o      [s]  output report. <required>
	    -cs_bam [s]  indexed bam file of case sample. <required>
	    -segm   [s]  segment region list. <required>
	                 table-delimited format: seg_id  chr  start_pos  end_pos

	   # Database and Metadata #
	    -ref    [s]  reference genome indexed fasta file for samtools mpileup func. <required>
	    -hrc    [s]  vcf file from Haplotype Reference Consortium, should be indexed. [optional]
	                 e.g., HRC.r1-1.EGA.GRCh37.chr10.haplotypes.vcf.gz from 1000G project.
	    -hrc_sl [s]  HRC sample list to extract from HRC vcf. <required with '-hrc'>
	                 note: one column one sample-ID, such as, NA18525.
	    -gm_bam [s]  indexed bam file of control (case's germline) sample. [optional]
	                 note: this option conflicts with '-hrc', use either of them.

	   # Software Required #
	    -bfl    [s]  bcfTools, minimum version: 1.3. <required>
	    -stl    [s]  SamTools, minimum version: 1.3. <required>

	   # Options #
	    -ploidy [s]  plodiy setting for bcftools call func. [GRCh37]
	    -minimb [f]  minimum imbalance ratio for Het SNP. [0.2]
	                 allow values from 0.1 to 0.9.
	    -minpc  [i]  minimum count of het SNP with enough imb in one segment. [1]
	    -minvq  [i]  minimum QUAL of one het SNP selected from your bam. [20]
	    -minmq  [i]  minimum MQ of one het SNP selected from your bam. [30]
	    -minpdp [i]  minimum depth of one het SNP position. [20]
	    -minadp [i]  minimum depth of one allele of het SNP. [5]
	    -minsdp [i]  minimum depth on one strand of one allele of het SNP. [2]
	    -bar10x      calculate 10x-barcode linkage. [optional]
	    -hic         calculate Hi-C reads linkage. [optional]

	    -h|help      Display this help info.

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
			[ segPhase_report => '' ],
			[ hsp_bam => {} ], # (h)ost (s)eg (p)hsae
			[ segment_list => undef ],
			# database
			[ idx_ref_fa => undef ],
			[ HRC_VCF => undef ],
			[ HRC_SampleList => undef ],
			[ gml_bam => {} ], # germline sample of (h)ost (s)eg (p)hsae sample
			# software
			[ bcftools => '' ],
			[ samtools => '' ],
			# options
			[ ploidy_set => 'GRCh37' ],
			[ HetSNP_minImb => 0.2 ],
			[ segm_HetSourceSNP_minCount => 1 ],
			[ minimum_SNP_quality => 20 ],
			[ minimum_MAP_quality => 30 ],
			[ minimum_Pos_depth => 20 ],
			[ minimum_Allele_depth => 5 ],
			[ minimum_Strand_depth => 2 ],
			[ use_10x_barc => 0 ],
			[ use_hic_read => 0 ],
			#--- intermediate variants
			[ temp_workspace => undef ],
			[ SNP_source => undef ], # gmBAM / HRC
			[ link_type => undef ], # 10x / HiC
			[ SNP_source_viewRegList => undef ],
			[ SNPsource_segm_VCF_gz => undef ],
			[ SNPsource_segm_SNP_list => undef ],
			[ caseBAM_HetSourceSNP_vcf_gz => undef ],
			# report_content -> theme = $theme
			# report_content -> HRCscore = $HRCscoreInEachSample
			# report_content -> longlink = $10x_barcScore / $hic_readScore
			[ report_content => {} ],
			# SegmentInfo -> id = $segm_id,
			# SegmentInfo -> chr = $chr,
			# SegmentInfo -> stp = $stpos,
			# SegmentInfo -> edp = $edpos,
			# SegmentInfo -> full_cn = $full_cn,
			# SegmentInfo -> minor_cn = $minor_cn,
			# SegmentInfo -> ImbHetSNP -> $SNP_pos = SNV_OB, check FuseSV::Objects::SNV_OB
			# SegmentInfo -> hetImbSNPCount = $hetImbSNPCount
			[ SegmentInfo => {} ],
			# HRC_Samples -> NO = $.
			# HRC_Samples -> HRC_idx = $VCF_theme_idx
			[ HRC_Samples => {} ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['segPhase_report'],
								  ['hsp_bam', 'filepath'],
								  ['gml_bam', 'filepath'],
								  ['segment_list'],
								  ['HRC_VCF'],
								  # ['idx_ref_fa'],
								  ['HRC_SampleList'],
								  ['samtools'],
								  ['bcftools']                    ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# objects constructions, as these objects belong to the input options
	$V_Href->{hsp_bam} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'case',    filepath => undef );
	$V_Href->{gml_bam} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'control', filepath => undef );

	# get options
	GetOptions(
		# input/output
		"-o:s"		=> \$V_Href->{segPhase_report},
		"-cs_bam:s"	=> \$V_Href->{hsp_bam}->{filepath},
		"-segm:s"	=> \$V_Href->{segment_list},
		# database
		"-ref:s"	=> \$V_Href->{idx_ref_fa},
		"-hrc:s"	=> \$V_Href->{HRC_VCF},
		"-hrc_sl:s"	=> \$V_Href->{HRC_SampleList},
		"-gm_bam:s"	=> \$V_Href->{gml_bam}->{filepath},
		# software
		"-bfl:s"	=> \$V_Href->{bcftools},
		"-stl:s"	=> \$V_Href->{samtools},
		# options
		"-plodiy:s"	=> \$V_Href->{ploidy_set},
		"-minimb:f"	=> \$V_Href->{HetSNP_minImb},
		"-minpc:i"	=> \$V_Href->{segm_HetSourceSNP_minCount},
		"-minvq:i"	=> \$V_Href->{minimum_SNP_quality},
		"-minmq:i"	=> \$V_Href->{minimum_MAP_quality},
		"-minpdp:i"	=> \$V_Href->{minimum_Pos_depth},
		"-minadp:i"	=> \$V_Href->{minimum_Allele_depth},
		"-minsdp:i"	=> \$V_Href->{minimum_Strand_depth},
		"-bar10x"	=> \$V_Href->{use_10x_barc},
		"-hic"		=> \$V_Href->{use_hic_read},
		# help
		"-h|help"	=> \$V_Href->{HELP},
		# for debug
		"-debug"	=> \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || ! $V_Href->{segPhase_report}
		     || ! file_exist( filePath => $V_Href->{hsp_bam}->{filepath} )
		     || ! file_exist( filePath => $V_Href->{segment_list} )
		     || ! file_exist( filePath => $V_Href->{idx_ref_fa} )
		     || ! (     file_exist( filePath => $V_Href->{gml_bam}->{filepath} )
		     		|| (   file_exist( filePath => $V_Href->{HRC_VCF} )
		     			&& file_exist( filePath => $V_Href->{HRC_SampleList} )
		     		   )
		     	  )
		     || ! file_exist( filePath => $V_Href->{bcftools} )
		     || ! file_exist( filePath => $V_Href->{samtools} )
		     || $V_Href->{segm_HetSourceSNP_minCount} <= 0
		     || ( $V_Href->{HetSNP_minImb} < 0.1 || $V_Href->{HetSNP_minImb} > 0.9 )
			);
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub phase_host_segment{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	#--------------#
	# prepare work #
	#--------------#
	&prepare;

	#----------------------#
	# load info of segment #
	#----------------------#
	&read_segment_list;

	#---------------------------------#
	# extract segments from SNPsource #
	#---------------------------------#
	&extract_segmVCF_from_SNP_source;

	#-------------------------------------#
	# extract het-sourceSNP from case bam #
	#-------------------------------------#
	&extract_het_sourceSNP_from_caseBAM;

	#--------------------------------------#
	# select imb het-source-SNP in caseBam #
	#--------------------------------------#
	&select_ImbHet_sourceSNP_in_csaeBam;

	#-------------------------------#
	# calculate m-m HRC phase-score #
	#-------------------------------#
	&calculate_segm_mVSm_HRCphaseScore;

	#-----------------------------#
	# 10x-barc/hic m-m link-score #
	#-----------------------------#
	&calculate_segm_mVsm_LongRangeLinkScore;

	#-----------------------------#
	# generate phase score report #
	#-----------------------------#
	&output_report;
}

#--- prepare work ---
sub prepare{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# workspace name
	$V_Href->{temp_workspace} = $V_Href->{segPhase_report}.'-tempDir';

	# check HRC amd germline-bam
	if( file_exist( filePath => $V_Href->{gml_bam}->{filepath} ) ){
		if( file_exist( filePath => $V_Href->{HRC_VCF} ) ){
			warn_and_exit "<ERROR>\tboth of HRC vcf ('-hrc' parm) and germline bam ('-gm_bam' parm) are assigned.\n"
								."\tplease use either of '-hrc' and '-gm_bam' parameters.\n";
		}
		elsif( file_exist( filePath => $V_Href->{HRC_SampleList} )
		){
			warn_and_exit "<ERROR>\tHRC sample list ('-hrc_sl' parm) should be with HRC vcf ('-hrc' parm), not for germline bam ('-gm_bam' parm).\n"
								."\tdo not use '-gm_bam' and '-hrc_sl' parameters together.\n";
		}
		$V_Href->{SNP_source} = 'gmBAM';
	}
	else{
		$V_Href->{SNP_source} = 'HRC';
	}

	# long linkage
	if( $V_Href->{use_10x_barc} ){
		$V_Href->{link_type} = '10x';
	}
	elsif( $V_Href->{use_hic_read} ){
		$V_Href->{link_type} = 'HiC';
	}

	# gmBAM should be with linkage
	if(    file_exist( filePath => $V_Href->{gml_bam}->{filepath} )
		&& !defined $V_Href->{link_type}
	){
		warn_and_exit "<ERROR>\tgermline bam ('-gm_bam' parm) should couple with either of '-bar10x' and '-hic' parameters.\n";
	}
}

#--- load segment information ---
sub read_segment_list{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};

	my $last_chr = undef;
	my $last_edpos = -1;
	open (SEGM, Try_GZ_Read($V_Href->{segment_list})) || die "fail read segment list: $!\n";
	while(<SEGM>){
		next if(/^#/); # comments
		my ($segm_id, $chr, $stpos, $edpos, $fullCN) = (split);
		# test recurrence
		if( exists $Segment_Href->{$segm_id} ){
			warn_and_exit "<ERROR>\tsegment (id=$segm_id) is duplicated.\n";
		}
		# test position retional
		if( $stpos > $edpos ){
			warn_and_exit "<ERROR>\tsegment (id=$segm_id) has reversed interval, start_pos=$stpos, end_pos=$edpos.\n";
		}
		# test chr
		if( defined $last_chr && $chr ne $last_chr ){
			warn_and_exit "<ERROR>\tsegment (id=$segm_id) has chr ($chr), different from previous chr ($last_chr).\n";
		}
		# test segment partition, must be sorted
		if( $last_edpos >= $stpos ){
			warn_and_exit "<ERROR>\tsegment (id=$segm_id) has start_pos ($stpos) smaller than last end_pos ($last_edpos).\n";
		}
		# test fullCN, optional
		my ($full_cn, $minor_cn) = (0, 0);
		if( length($fullCN) != 0 ){
			($full_cn, $minor_cn) = ($fullCN =~ /^c?n?(\d+)m(\d+)$/);
			if(    !defined $full_cn
				|| !defined $minor_cn
				|| $minor_cn >= ($full_cn - $minor_cn) # M is M, m is m, and should imbalance (!=)
				|| $minor_cn == 0 # must het
			){
				warn_and_exit "<ERROR>\tfullCN ($fullCN) of segment (id=$segm_id) is not available.\n";
			}
		}
		# load up
		$Segment_Href->{$segm_id} = {
										id => $segm_id,
										chr => $chr,
										stp => $stpos,
										edp => $edpos,
										full_cn => $full_cn,
										minor_cn => $minor_cn,
										ImbHetSNP => {},
										hetImbSNPCount => 0
									};
		# update
		$last_chr = $chr;
		$last_edpos = $edpos;
	}
	close SEGM;

	# inform
	stout_and_sterr `date`."[INFO]\tLoad segment list OK.\n"
								 ."\t$V_Href->{segment_list}\n";
}

#--- extract segments from VCF of SNP_source ---
sub extract_segmVCF_from_SNP_source{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};

	# make workspace
	`rm -rf $V_Href->{temp_workspace}`;
	`mkdir -p $V_Href->{temp_workspace}`;

	# make region list
	$V_Href->{SNP_source_viewRegList} = File::Spec->catfile($V_Href->{temp_workspace}, 'SNPsource.'.$V_Href->{SNP_source}.'.view_region.list');
	open (VREG, Try_GZ_Write($V_Href->{SNP_source_viewRegList})) || die "fail to write view_region list: $!\n";
	for my $segm_id (sort {$Segment_Href->{$a}->{stp} <=> $Segment_Href->{$b}->{stp}} keys %$Segment_Href){
		my @content = map $Segment_Href->{$segm_id}->{$_}, qw/ chr stp edp /;
		print VREG join( "\t", @content )."\n";
	}
	close VREG;

	# SNP source segmVCF related files
	$V_Href->{SNPsource_segm_VCF_gz} = File::Spec->catfile($V_Href->{temp_workspace}, 'SNPsource.'.$V_Href->{SNP_source}.'.segments.vcf.gz');
	# select SNP source
	if( $V_Href->{SNP_source} eq 'HRC' ){
		# from HRC
		&extract_segmVCF_from_HRC;
	}
	elsif( $V_Href->{SNP_source} eq 'gmBAM' ){
		# from germline bam
		&extract_segmVCF_from_gmBAM;
	}

	$V_Href->{SNPsource_segm_SNP_list} = File::Spec->catfile($V_Href->{temp_workspace}, 'SNPsource.'.$V_Href->{SNP_source}.'.segments.SNP.pos.list');
	# obtain SNPsource SNP positions
	my $SNPsource_segm_SNPlist_cmd = "$V_Href->{bcftools} view -H $V_Href->{SNPsource_segm_VCF_gz} | cut -f1,2 > $V_Href->{SNPsource_segm_SNP_list}";
	trible_run_for_success($SNPsource_segm_SNPlist_cmd, 'get_SNPsource_SNP_pos_list');

	# inform
	stout_and_sterr `date`."[INFO]\tExtract segments' $V_Href->{SNP_source} VCF OK.\n"
								 ."\t$V_Href->{SNPsource_segm_VCF_gz}\n"
								 ."\t$V_Href->{SNPsource_segm_SNP_list}\n";
}

#--- extract segments VCF of gmBAM ---
sub extract_segmVCF_from_gmBAM{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get het SNP
	$V_Href->{gml_bam}->get_region_alt_vcf_gz(
		vcfgz => $V_Href->{SNPsource_segm_VCF_gz},
		cmd_name => 'germline_bam_het_SNP',
		samtools => $V_Href->{samtools},
		bcftools => $V_Href->{bcftools},
		pos_list => $V_Href->{SNP_source_viewRegList},
		idx_ref => $V_Href->{idx_ref_fa},
		ploidy_set => $V_Href->{ploidy_set},
		no_indel => 1,
		snpGap => 5,
		IndelGap => 5,
		min_altQual => $V_Href->{minimum_SNP_quality},
		min_mapQual => $V_Href->{minimum_MAP_quality},
		min_posDepth => $V_Href->{minimum_Pos_depth},
		min_alleleDepth => $V_Href->{minimum_Allele_depth},
		min_strdDepth => $V_Href->{minimum_Strand_depth},
		only_hetAlt => 1
	);
}

#--- extract segments VCF of HRC ---
sub extract_segmVCF_from_HRC{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# extract segment from HRC
	my $HRC_segm_VCF_cmd = "($V_Href->{bcftools} view -R $V_Href->{SNP_source_viewRegList} -O z -o $V_Href->{SNPsource_segm_VCF_gz} $V_Href->{HRC_VCF}) && ($V_Href->{bcftools} index $V_Href->{SNPsource_segm_VCF_gz})";
	trible_run_for_success($HRC_segm_VCF_cmd, 'HRCvcf_segm_view');
}

#--- select imbalance het-source-SNP ---
sub select_ImbHet_sourceSNP_in_csaeBam{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};

	for my $segm_id (sort keys %$Segment_Href){
		my $this_segm = $Segment_Href->{$segm_id};
		my $chr = $this_segm->{chr};
		my $stp = $this_segm->{stp};
		my $edp = $this_segm->{edp};
		# load up het source SNP with enough imb
		my $ImbHetSNP_Href = $this_segm->{ImbHetSNP};
		my $region = "$chr:$stp-$edp";
		open (SEGVCF, "$V_Href->{bcftools} view -r $region -H $V_Href->{caseBAM_HetSourceSNP_vcf_gz}|") || die "fail read het_Source_VCF: $!\n";
		while(<SEGVCF>){
			my ($chr, $pos, $refA, $altA) = (split)[0,1,3,4];
			my ($refC, $altC) = (/AD=(\d+),(\d+);/);
			# load snv object
			my $SNV_OB = FuseSV::Objects::SNV_OB->new(
				chr => $chr,
				pos => $pos,
				refA => $refA,
				altA => $altA,
				refC => $refC,
				altC => $altC,
				min_imb => $V_Href->{HetSNP_minImb}
			);
			# min-imb filter
			next if( $SNV_OB->get_imb < $V_Href->{HetSNP_minImb} );
			# initialize hash container
			$SNV_OB->ini_hash_Href( key_Aref => ['HRC'] ) if( $V_Href->{SNP_source} eq 'HRC' );
			# record
			$ImbHetSNP_Href->{$pos} = $SNV_OB;
			# count in
			$this_segm->{hetImbSNPCount} ++;
		}
		close SEGVCF;
		# check avaible SNP count
		if( $this_segm->{hetImbSNPCount} < $V_Href->{segm_HetSourceSNP_minCount} ){
			$this_segm->{goodForUse} = 0;
			# inform
			stout_and_sterr "<WARN>\tcount of het SNP with enough imb in segment ($segm_id) is less than required.\n"
								 ."\tits count: $this_segm->{hetImbSNPCount}, required count: $V_Href->{segm_HetSourceSNP_minCount}.\n";
		}
		else{
			$this_segm->{goodForUse} = 1;
			# inform
			stout_and_sterr `date`."[INFO]\tload het imb source SNP in segment $segm_id OK.\n";
		}
	}
}

#--- load HRC sample list ---
sub load_HRC_sample_list{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Samples_Href = $V_Href->{HRC_Samples};

	open (SL, Try_GZ_Read($V_Href->{HRC_SampleList})) || die "fail read HRC sample list: $!\n";
	while(<SL>){
		next if(/^#/);
		my ($sampleID) = (split)[0];
		if( exists $Samples_Href->{$sampleID} ){
			warn "<WARN>\t$sampleID is duplicated in HRC_sample.list.\n";
			next;
		}
		$Samples_Href->{$sampleID}->{NO} = $.;
		$Samples_Href->{$sampleID}->{HRC_idx} = 0;
	}
	close SL;

	# inform
	stout_and_sterr `date`."[INFO]\tRead HRC sample list OK.\n";
}

#--- load het-imb-SNP from HRC based on sample list ---
sub load_hetImbSNP_from_HRC{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};
	my $Samples_Href = $V_Href->{HRC_Samples};

	for my $segm_id (sort keys %$Segment_Href){
		my $this_segm = $Segment_Href->{$segm_id};
		next unless($this_segm->{goodForUse});
		my $this_segm_chr = $this_segm->{chr};
		my $ImbHetSNP_Href = $this_segm->{ImbHetSNP};

		# create region list for HRC vcf bcftools view
		my $SNP_region_list = File::Spec->catfile($V_Href->{temp_workspace}, "segm.$segm_id.hetImbSNP.region.list");
		open (SNPL, Try_GZ_Write($SNP_region_list)) || die "fail create $segm_id.hetImbSNP.region.list: $!\n";
		print SNPL $ImbHetSNP_Href->{$_}->get_chr_pos(delimit=>"\t",region=>1)."\n" for sort {$a<=>$b} keys %$ImbHetSNP_Href;
		close SNPL;

		# read HRC vcf of these SNP pos in samples one by one
		open (HRC, "$V_Href->{bcftools} view -R $SNP_region_list $V_Href->{SNPsource_segm_VCF_gz} |") || die "fail read SNPsource_segm_VCF_gz: $!\n";
		while(<HRC>){
			next if(/^##/); # double ##
			if(/^#[^#]/){ # VCF theme
				my @theme = split;
				# find samples' idx
				for my $idx ( 9 .. $#theme ){
					my $sampleID = $theme[$idx];
					if( exists $Samples_Href->{$sampleID} ){
						$Samples_Href->{$sampleID}->{HRC_idx} = $idx;
					}
				}
				# delete samples not found
				for my $sampleID (sort keys %$Samples_Href){
					if( $Samples_Href->{$sampleID}->{HRC_idx} == 0 ){
						stout_and_sterr "<WARN>\tsample $sampleID does not exist in HRC VCF, and skip it.\n";
						delete $Samples_Href->{$sampleID};
					}
				}
			}
			else{ # content
				my @info = split;
				my ($pos, $refA, $altA) = @info[1,3,4];
				next unless( exists $ImbHetSNP_Href->{$pos} );
				my $SNV_OB = $ImbHetSNP_Href->{$pos};
				# test ref alt equal to your_bam
				if(    $refA ne $SNV_OB->get_refAllele
					|| $altA ne $SNV_OB->get_altAllele
				){
					stout_and_sterr "<WARN>\tHRC alleles at pos $pos are different from your_bam. skip it. ".
								 			"HRC (refA=$refA, altA=$altA); YourBam (refA=$SNV_OB->{refA}, altA=$SNV_OB->{altA}).\n";
					delete $ImbHetSNP_Href->{$pos};
					next;
				}
				# each sample
				for my $sampleID (sort keys %$Samples_Href){
					my $idx = $Samples_Href->{$sampleID}->{HRC_idx};
					my $phGT = $info[$idx];
					# phase info
					my ($GT_L, $GT_R) = (split /\|/, $phGT)[0,1];
					$SNV_OB->load_HRC_sample_info(
						sampleID => $sampleID,
						GT_L => $GT_L,
						Allele_L => ( $GT_L == 0 ? $refA : $altA ),
						GT_R => $GT_R,
						Allele_R => ( $GT_R == 0 ? $refA : $altA )
					);
				}
			}
		}
		close HRC;
		# inform
		stout_and_sterr "[INFO]\tload HRC sample phase info of hetImbSNP in segment $segm_id OK.\n";

		# output the details of your_bam Het-HRC-SNP in all HRC samples
		my $segm_SNP_details = File::Spec->catfile($V_Href->{temp_workspace}, "segm.$segm_id.hetImbSNP.HRC_phased_info.txt");
		open (SNPDT, Try_GZ_Write($segm_SNP_details)) || die "fail create $segm_id.hetImbSNP.HRC_phased_info.txt: $!\n";
		my @theme = qw/ segm chr pos refA altA Major minor /;
		my @samples = sort {$Samples_Href->{$a}->{NO} <=> $Samples_Href->{$b}->{NO}} keys %$Samples_Href;
		push @theme, @samples;
		print SNPDT '#'.join("\t", @theme)."\n";
		for my $pos (sort {$a<=>$b} keys %$ImbHetSNP_Href){
			my $SNV_OB = $ImbHetSNP_Href->{$pos};
			my @contents = ($segm_id, $this_segm_chr, $pos);
			push @contents, $SNV_OB->get_refAllele;
			push @contents, $SNV_OB->get_altAllele;
			push @contents, $SNV_OB->get_majorAllele;
			push @contents, $SNV_OB->get_minorAllele;
			for my $sampleID (@samples){
				if( $SNV_OB->has_HRC_sample(sampleID=>$sampleID) ){
					push @contents, $SNV_OB->get_HRC_sample_allele( sampleID=>$sampleID, delimit=>'|').','.
									$SNV_OB->get_HRC_sample_Genotp( sampleID=>$sampleID, delimit=>'|');
				}
				else{
					push @contents, '-|-,-|-';
				}
			}
			print SNPDT join("\t", @contents)."\n";
		}
		close SNPDT;
	}
}

#--- extract het-sourceSNP from case bam ---
sub extract_het_sourceSNP_from_caseBAM{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$V_Href->{caseBAM_HetSourceSNP_vcf_gz} = File::Spec->catfile($V_Href->{temp_workspace}, 'caseBam.source_'.$V_Href->{SNP_source}.'.Het_SNP.vcf.gz');
	# get het Source SNP
	$V_Href->{hsp_bam}->get_region_alt_vcf_gz(
		vcfgz => $V_Href->{caseBAM_HetSourceSNP_vcf_gz},
		cmd_name => 'caseBam_het_source_SNP',
		samtools => $V_Href->{samtools},
		bcftools => $V_Href->{bcftools},
		pos_list => $V_Href->{SNPsource_segm_SNP_list},
		idx_ref => $V_Href->{idx_ref_fa},
		ploidy_set => $V_Href->{ploidy_set},
		no_indel => 1,
		snpGap => 5,
		IndelGap => 5,
		min_altQual => $V_Href->{minimum_SNP_quality},
		min_mapQual => $V_Href->{minimum_MAP_quality},
		min_posDepth => $V_Href->{minimum_Pos_depth},
		min_alleleDepth => $V_Href->{minimum_Allele_depth},
		min_strdDepth => $V_Href->{minimum_Strand_depth},
		only_hetAlt => 1
	);

	# inform
	stout_and_sterr `date`."[INFO]\tExtract het sourceSNP from case bam OK.\n"
								 ."\t$V_Href->{caseBAM_HetSourceSNP_vcf_gz}\n";
}

#--- calculate m-m phase-score in each sample ---
sub calculate_segm_mVSm_HRCphaseScore{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};
	my $Samples_Href = $V_Href->{HRC_Samples};

	# not use 
	return if( $V_Href->{SNP_source} ne 'HRC' );

	#----------------------------#
	# load HRC check sample list #
	#----------------------------#
	&load_HRC_sample_list;

	#---------------------------#
	# load het-imb-SNP from HRC #
	#---------------------------#
	&load_hetImbSNP_from_HRC;

	# sort the segments from smaller pos to larger pos
	my @segm = sort {$Segment_Href->{$a}->{stp} <=> $Segment_Href->{$b}->{stp}} keys %$Segment_Href;
	## m-m phase score across samples
	$V_Href->{report_content}->{HRCscore} = '';
	my %segm_pair_NPC = map { ( $_, { N=>[], P=>[] } ) } ( 0 .. $#segm-1 );
	for my $sampleID (sort {$Samples_Href->{$a}->{NO} <=> $Samples_Href->{$b}->{NO}} keys %$Samples_Href){
		my @contents = ($sampleID);
		for my $idx ( 0 .. $#segm-1 ){
			my $this_segm_id = $segm[$idx];
			my $next_segm_id = $segm[$idx+1];
			# prepare HRC het imb SNP for phase score in this sample
			my $this_segm_SNP_Href = $Segment_Href->{$this_segm_id}->{ImbHetSNP};
			my @this_segm_Pos = grep $this_segm_SNP_Href->{$_}->is_het_HRC_sample( sampleID=>$sampleID ),
								keys %$this_segm_SNP_Href;
			my $this_segm_PosCount = scalar(@this_segm_Pos);
			my $next_segm_SNP_Href = $Segment_Href->{$next_segm_id}->{ImbHetSNP};
			my @next_segm_Pos = grep $next_segm_SNP_Href->{$_}->is_het_HRC_sample( sampleID=>$sampleID ),
								keys %$next_segm_SNP_Href;
			my $next_segm_PosCount = scalar(@next_segm_Pos);
			my $PosCount_Str = $this_segm_PosCount.','.$next_segm_PosCount;
			# not enough SNP positions
			if(    $this_segm_PosCount < $V_Href->{segm_HetSourceSNP_minCount}
				|| $next_segm_PosCount < $V_Href->{segm_HetSourceSNP_minCount}
			){
				push @contents, $PosCount_Str.',NA';
				next;
			}
			# calculate phase score
			my $phase_score = 0;
			for my $pos1 (@this_segm_Pos){
				my $pos1_minor_A = $this_segm_SNP_Href->{$pos1}->get_minorAllele;
				my ($pos1_lft_A, $pos1_rgt_A) = $this_segm_SNP_Href->{$pos1}->get_HRC_sample_allele( sampleID=>$sampleID );
				for my $pos2 (@next_segm_Pos){
					my $pos2_minor_A = $next_segm_SNP_Href->{$pos2}->get_minorAllele;
					my ($pos2_lft_A, $pos2_rgt_A) = $next_segm_SNP_Href->{$pos2}->get_HRC_sample_allele( sampleID=>$sampleID );
					# score
					if(    ( $pos1_minor_A eq $pos1_lft_A && $pos2_minor_A eq $pos2_lft_A )
						|| ( $pos1_minor_A eq $pos1_rgt_A && $pos2_minor_A eq $pos2_rgt_A )
					){
						$phase_score ++;
					}
					else{
						$phase_score --;
					}
				}
			}
			$phase_score = sprintf "%.3f", $phase_score / ( $this_segm_PosCount * $next_segm_PosCount );
			push @contents, $PosCount_Str.','.$phase_score;
			# store for all-calculation
			my $sign = ( $phase_score < 0 ? 'N' : 'P' );
			push @{$segm_pair_NPC{$idx}->{$sign}}, $phase_score;
		}
		$V_Href->{report_content}->{HRCscore} .= join("\t", @contents)."\n";
	}
	# calculate all samples' average
	my $content_of_all = 'All_Sample';
	for my $idx (sort {$a<=>$b} keys %segm_pair_NPC){
		my $segm_pair_NPC_Href = $segm_pair_NPC{$idx};
		my @segm_pair_content = ();
		for my $sign (qw / N P /){
			my $count = scalar( @{$segm_pair_NPC_Href->{$sign}} );
			my $TMM_mean =   $count == 0
						   ? 'NA'
						   : get_value_mean( value_Aref => $segm_pair_NPC_Href->{$sign}, TrimRatio => 0.1, ratio_digit => 3 );
			push @segm_pair_content, "$sign:$count:$TMM_mean";
		}
		$content_of_all .= "\t".join(',',@segm_pair_content);
	}
	$V_Href->{report_content}->{HRCscore} = $content_of_all."\n".$V_Href->{report_content}->{HRCscore};

	# inform
	stout_and_sterr `date`."[INFO]\tcalculate segment phase score in HRC samples OK.\n";
}

#--- calculate linkage score of m-m from long-range DNA info ---
## based on 10x-barcode or Hi-C pair-end-id
sub calculate_segm_mVsm_LongRangeLinkScore{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};

	return if( !defined $V_Href->{link_type} );

	my $link_type = $V_Href->{link_type};
	my @longlink_report_content = ("$link_type-Link");

	my @segm = sort {$Segment_Href->{$a}->{stp} <=> $Segment_Href->{$b}->{stp}} keys %$Segment_Href;
	# record marker stat of each segm, and compare with last segm
	my $last_seg_id = undef;
	for my $segm_id (@segm){
		my $segm_SNP_Href = $Segment_Href->{$segm_id}->{ImbHetSNP};
		# record allele marker stat
		for my $pos (sort {$a<=>$b} keys %$segm_SNP_Href){
			my $SNV_OB = $segm_SNP_Href->{$pos};
			$SNV_OB->ini_hash_Href( key_Aref => [$link_type] ); # initialize
			# stat allele maker
			my $AmakerStat_Href = $V_Href->{hsp_bam}->get_allele_marker_stat(
				chr => $SNV_OB->get_chr,
				pos => $pos,
				marker => $link_type,
				samtools => $V_Href->{samtools}
			);
			# load link stat
			$SNV_OB->load_long_link_info( allele2linkmarker => $AmakerStat_Href, link_type => $link_type );
		}
		# inform
		stout_and_sterr "[INFO]\tload segment $segm_id $link_type linkage OK.\n";
		# compare with last segm
		if( defined $last_seg_id ){
			# details file
			my $segmlinkCmpDetails_file = File::Spec->catfile($V_Href->{temp_workspace}, "linkCmp.$link_type.segm.$last_seg_id.vs.$segm_id.list");
			open (CMP, Try_GZ_Write($segmlinkCmpDetails_file)) || die "fail write segmlinkCmpDetails: $!\n";
			print CMP '#'.join("\t",qw/ chr segm_1 pos_1 segm_2 pos_2 mSameC mDiffC mSameImb /)."\n";
			my ($mSame_Nc, $mSame_Pc) = (0,0);
			my $last_segm_SNP_Href = $Segment_Href->{$last_seg_id}->{ImbHetSNP};
			for my $last_segm_pos (sort {$a<=>$b} keys %$last_segm_SNP_Href){
				my $ltsg_SNV_OB = $last_segm_SNP_Href->{$last_segm_pos};
				for my $this_segm_pos (sort {$a<=>$b} keys %$segm_SNP_Href){
					my $thsg_SNV_OB = $segm_SNP_Href->{$this_segm_pos};
					# count allele linkage
					my $mSame_link_count = $thsg_SNV_OB->cmp_otherSNV_link(oth_SNVob=>$ltsg_SNV_OB, link_type=>$link_type, cmp_mode=>'mSame');
					my $mDiff_link_count = $thsg_SNV_OB->cmp_otherSNV_link(oth_SNVob=>$ltsg_SNV_OB, link_type=>$link_type, cmp_mode=>'mDiff');
					# mSame linkage imb
					my $mSame_imb;
					if(    $mSame_link_count == 0
						&& $mDiff_link_count == 0
					){
						$mSame_imb = 'NA';
					}
					else{
						$mSame_imb = sprintf "%.3f", ($mSame_link_count-$mDiff_link_count) / max($mSame_link_count,$mDiff_link_count);
						$mSame_Nc ++ if( $mSame_imb <= -0.2 );
						$mSame_Pc ++ if( $mSame_imb >=  0.2 );
					}
					# output details
					print CMP join("\t", $ltsg_SNV_OB->get_chr,
										 $last_seg_id,
										 $ltsg_SNV_OB->get_pos,
										 $segm_id,
										 $thsg_SNV_OB->get_pos,
										 $mSame_link_count,
										 $mDiff_link_count,
										 $mSame_imb
								  )."\n";
				}
			}
			close CMP;
			# prepare outputs
			my ($mSame_Nr, $mSame_Pr);
			if(    $mSame_Pc == 0
				&& $mSame_Nc == 0
			){
				($mSame_Nr, $mSame_Pr) = ('NA', 'NA');
			}
			else{
				$mSame_Nr = sprintf "%.3f", $mSame_Nc / ($mSame_Pc+$mSame_Nc) * -1;
				$mSame_Pr = sprintf "%.3f", $mSame_Pc / ($mSame_Pc+$mSame_Nc);
			}
			push @longlink_report_content, "N:$mSame_Nc:$mSame_Nr,P:$mSame_Pc:$mSame_Pr";
			# inform
			stout_and_sterr "[INFO]\t$link_type SharedLinkCmp of segments $last_seg_id.vs.$segm_id OK.\n";
		}
		# update
		$last_seg_id = $segm_id;
	}

	$V_Href->{report_content}->{longlink} = join("\t", @longlink_report_content)."\n";
}

#--- write phase/linkage report ---
sub output_report{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{SegmentInfo};

	# report theme
	my @theme;
	my @theme_i0;
	push @theme_i0, 'HRC_idvl' if( $V_Href->{SNP_source} eq 'HRC' );
	push @theme_i0, $V_Href->{link_type}.'_link' if( defined $V_Href->{link_type} );
	push @theme, join(';', @theme_i0);
	## sort the segments from smaller pos to larger pos
	my @segm = sort {$Segment_Href->{$a}->{stp} <=> $Segment_Href->{$b}->{stp}} keys %$Segment_Href;
	push @theme, $segm[$_]  .'_m('.$Segment_Href->{$segm[$_]}->{hetImbSNPCount}  .'),'.
				 $segm[$_+1].'_m('.$Segment_Href->{$segm[$_+1]}->{hetImbSNPCount}.')' for ( 0 .. $#segm-1 );
	$V_Href->{report_content}->{theme} = '#'.join("\t", @theme)."\n";

	# full report
	open (REPORT, Try_GZ_Write($V_Href->{segPhase_report})) || die "fail to write segPhase_report: $!.\n";
	for my $info_key ( qw/ theme longlink HRCscore / ){
		if( exists $V_Href->{report_content}->{$info_key} ){
			print REPORT $V_Href->{report_content}->{$info_key};
		}
	}
	close REPORT;

	# inform
	stout_and_sterr `date`."[INFO]\toutput segment phase report OK.\n"
								 ."\t$V_Href->{segPhase_report}\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
