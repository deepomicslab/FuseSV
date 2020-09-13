package FuseSV::Virus_Type::PrepareForIntegDetection;

use strict;
use warnings;
use File::Spec qw/ catfile /;
use List::Util qw/ min max sum /;
use Getopt::Long;
use BioFuse::Util::Log qw/ stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ BWA_index_fasta Faidx_Dict_fasta /;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Virus_Type::VirusTypeMain;
# use FuseSV::Virus_Type::CandReads;
use FuseSV::Virus_Type::ReconstrcutVirusVariants;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  prepare_for_virus_integration_detection
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'PrepareForIntegDetection';
#----- version --------
$VERSION = "0.13";
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
						prepare_for_virus_integration_detection
						get_parameter_for_integration_detection
						get_makeupRef_for_integration_detection
						get_PSLfile_Genelist_for_integration_detection
						get_merge_bam_for_integration_detection
						select_reads_alignments_for_integration_bam
						choose_vv_source_alignment_to_output
						simplify_bam_for_integration
						is_this_virus_variants_been_selected
						is_this_virus_variants_Nstr_OK
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} prep_integ <[Options]>

	 Note: If you want to use FuseSV to detect viral integraion, this operation is needed.

	 Options:

	   # Inputs and Outputs #
	   # Note: these options below should be as \033[4m \bsame \b\033[0m as previous steps of 'virus_type' #
	    -id     [s]  assign a symbol name for this process, e.g., sample-id. <required>
	                  only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
	    -odir   [s]  directory to store outputs, must exist already. <required>
	    -host   [s]  the .fai file of only host genome reference. <required>
	                  such as pure human reference fai file, created by samtools faidx.
	    -vpf    [s]  virus database folder. <required>
	                  You can make it via function '\033[3m \bmake_virusdb \b\033[0m'.
	                  OR Download these widely used from FuseSV official website, e.g., HPV.

	   # Software Required #
	    -bwa    [s]  BWA tools, minimum version: 0.7.9. <required>
	    -stl    [s]  SamTools, minimum version: 1.3. <required>
	    -pic    [s]  Picard jar file, minimum version: 2.1.0. <required>
	    -java   [s]  Java compiler, minimum version 1.8. ['java']

	    # Prepare Integration Options #
	    -simple      to get simplified bam for quicker but maybe more false-negative integration detection. [Disabled]
	    -vno    [s]  to select virus variants (its NO.) for virus integration detection. [1]
	                  default value is 1, which means only use the FIRST virus variant;
	                  set as 'all' to use ALL virus variants;
	                  set as '1,2' to use FIRST and SECOND variants (if it has, automatically judgement).
	     ** Note: options below just work on virus variants (besides NO.1) selected by option '-vno'.
	    -midp   [i]  minimum average depth of viral genome. [5]
	    -micr   [f]  minimum covered ratio of viral genome. [0.5]

	    -h|help   Display this help info.

 	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
	    # -mxnc   [i]  maximum 'N-string' count allowed in each virus variant genome. [5]
	    # -minl   [i]  minimum 'N-string' length to count in each virus variant genome. [5]
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
	FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool;
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-id:s"		=> \$V_Href->{sample_ID},
		"-odir:s"	=> \$V_Href->{virusAnaFolder},
		"-host:s"	=> \$V_Href->{host_fai},
		"-vpf:s"	=> \$V_Href->{virus_db_folder},
		# software
		"-bwa:s"	=> \$V_Href->{bwa},
		"-stl:s"	=> \$V_Href->{samtools},
		"-java:s"	=> \$V_Href->{java},
		"-pic:s"	=> \$V_Href->{picard_jar},
		# options
		## multiple fork operations
		"-jvmem:s"	=> \$V_Href->{mem_for_jvm}, # hidden option
		## select virus variants
		"-simple"	=> \$V_Href->{Get_simple_BamForInteg},
		"-vno:s"	=> \$V_Href->{VariantNO_For_Analysis},
		"-midp:i"	=> \$V_Href->{MinDepth_For_IntegAna},
		"-micr:f"	=> \$V_Href->{MinCoverRatio_For_IntegAna},
		"-mxnc:i"	=> \$V_Href->{MaxAllow_Nstr_Count}, # hidden option
		"-minl:i"	=> \$V_Href->{Nstr_MinLength}, # hidden option
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
			 || (!$V_Href->{sample_ID} || $V_Href->{sample_ID} =~ /[^\w\-]/)
			 || (!$V_Href->{virusAnaFolder} || !-d $V_Href->{virusAnaFolder})
			 || (!$V_Href->{host_fai} || !-e $V_Href->{host_fai})
			 || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder})
			 || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
			 || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
			 || $V_Href->{MinCoverRatio_For_IntegAna} < 0
			 || $V_Href->{MinCoverRatio_For_IntegAna} > 1
			 || $V_Href->{MinDepth_For_IntegAna} < 0
			 # || $V_Href->{Nstr_MinLength} < 0
			 || $V_Href->{MaxAllow_Nstr_Count} < 0
			);
}

#--- gather the virus variants file and information for next virus integration detection ---
sub prepare_for_virus_integration_detection{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Prepare_For_Integ_Href = $V_Href->{Prepare_For_Integ};

	# file path
	$V_Href->{Prepare_Integ_folder} = GetPath( filekey => 'Prepare_Integ_folder' );

	#---------------------------#
	# mode 4. prepare files for #
	# virus integration detect  #
	#---------------------------#

	# run mode check
	unless(    ( $V_Href->{run_mode} == 0 && ! $V_Href->{Skip_VirusItegPrep} ) # run_whole and donot skip this step
			||   $V_Href->{run_mode} == 4 # just this step
		  ){
		return;
	}
	# if just start from this step, the basic info need to reload.
	if( $V_Href->{run_mode} == 4 ){
		# load the host segment as the previous seg
		# FuseSV::Virus_Type::CandReads->read_host_fai;
	}

	# prepare folder
	`rm -rf $V_Href->{Prepare_Integ_folder}`;
	`mkdir -p $V_Href->{Prepare_Integ_folder}`;

	#-------------------------------------------------#
	# prepare the virus PSL file, gene.psl; trans.psl #
	# and gene list for FuseSV, and get ext-length    #
	#-------------------------------------------------#
	$Prepare_For_Integ_Href->{GenePSL}  = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.virus.gene.psl');
	$Prepare_For_Integ_Href->{TransPSL} = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.virus.trans.psl');
	$Prepare_For_Integ_Href->{Variant_Extend_Href} = {};
	$Prepare_For_Integ_Href->{GeneList} = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.virus.gene.list');
	&get_PSLfile_Genelist_for_integration_detection(
			{
				GenePSL => $Prepare_For_Integ_Href->{GenePSL},
				TransPSL => $Prepare_For_Integ_Href->{TransPSL},
				Variant_Extend_Href => $Prepare_For_Integ_Href->{Variant_Extend_Href},
				GeneList => $Prepare_For_Integ_Href->{GeneList}
			}
		);

	#--------------------------------------------------------#
	# prepare the virus makeup reference, extended if needed #
	#--------------------------------------------------------#
	my $makeup_ref_folder = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.virus_makeup_ref');
	`rm -rf $makeup_ref_folder`;
	`mkdir -p $makeup_ref_folder`;
	$Prepare_For_Integ_Href->{MakeupRef} = File::Spec->catfile($makeup_ref_folder, $V_Href->{ForVitgPrefix}.'.virus.makeup.fa');
	&get_makeupRef_for_integration_detection(
			{
				MakeupRef_Fa => $Prepare_For_Integ_Href->{MakeupRef},
				Variant_Extend_Href => $Prepare_For_Integ_Href->{Variant_Extend_Href}
			}
		);

	#--------------------------------#
	# merge the Nsort source bam and #
	# the virus-variant aligned bam, #
	# select best alignment for each #
	# paired-end reads.              #
	#--------------------------------#
	$Prepare_For_Integ_Href->{BamForInteg} = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.bam');
	&get_merge_bam_for_integration_detection(
			{
				BamForInteg => $Prepare_For_Integ_Href->{BamForInteg},
				Variant_Extend_Href => $Prepare_For_Integ_Href->{Variant_Extend_Href}
			}
		);

	#--------------------------------------#
	# prepare the parameter doc for FuseSV #
	#--------------------------------------#
	$Prepare_For_Integ_Href->{ParaTxt} = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, $V_Href->{ForVitgPrefix}.'.parameters.txt');
	&get_parameter_for_integration_detection(
			{
				ParaTxt => $Prepare_For_Integ_Href->{ParaTxt},
				Variant_Extend_Href => $Prepare_For_Integ_Href->{Variant_Extend_Href}
			}
		);

	# run mode check, if not go on, just stop here.
	if( $V_Href->{run_mode} == 4 ){ # finish at this point in mode 4
		exit(0);
	}
}

#--- get gene psl and trans psl file for integration detection ---
sub get_PSLfile_Genelist_for_integration_detection{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $GenePSLforInteg  = $Option_Href->{GenePSL};
	my $TransPSLforInteg = $Option_Href->{TransPSL};
	my $Variant_Extend_Href = $Option_Href->{Variant_Extend_Href};
	my $GeneListForInteg = $Option_Href->{GeneList};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

	# get all virus variants related segname
	my %Variant_SegName = map { ( $Variant_Info_Href->{$_}->{OrigRefSeg},
									{
										variant_NO => $_,
										RefSeg => $Variant_Info_Href->{$_}->{RefSeg}
									} 
								) } keys %$Variant_Info_Href;

	#----------#
	# Gene PSL #
	#----------#
	my %yes_psl_refseg;
	open (S_GPSL,Try_GZ_Read($V_Href->{virus_orig_gene_psl})) || die "fail read $V_Href->{virus_orig_gene_psl}: $!\n";
	open (O_GPSL,Try_GZ_Write($GenePSLforInteg)) || die "fail write $GenePSLforInteg: $!\n";
	open (GLIST,Try_GZ_Write($GeneListForInteg)) || die "fail write $GeneListForInteg: $!\n";
	while (<S_GPSL>){
		my @psl_info = split /\t+/;
		my $virus_refseg_Sref = \ $psl_info[13];
		my ($virus_refseg, $extend_len);
		# extended or not
		if( $$virus_refseg_Sref =~ /^(.+)_Extended(\d+)$/i ){ # extended virus
			($virus_refseg, $extend_len) = ($1, $2);
		}
		else{ # not-extended
			($virus_refseg, $extend_len) = ($$virus_refseg_Sref, 0);
		}
		# is virus variant?
		if( exists($Variant_SegName{$virus_refseg}) ){
			my $mod_virus_refseg = $Variant_SegName{$virus_refseg}->{RefSeg};
			my $mod_virus_vNO = $Variant_SegName{$virus_refseg}->{variant_NO};
			# virus variants selection filtration
			next unless( &is_this_virus_variants_been_selected({vNO=>$mod_virus_vNO}) );
			# record for trans psl
			$yes_psl_refseg{$$virus_refseg_Sref} = $mod_virus_refseg;
			# update the psl file refseg as the variant segname
			$$virus_refseg_Sref = $mod_virus_refseg;
			# record variant needs to extend
			if( $extend_len != 0 ){
				$Variant_Extend_Href->{$mod_virus_refseg} = {
					extend_len => $extend_len,
					variant_NO => $mod_virus_vNO
				};
			}
			# ouput to the Gene PSL for integration detection
			print O_GPSL join("\t", @psl_info); # last one has "\n"
			# output the gene list specially for FuseSV
			print GLIST join("\t", (@psl_info[9,8,13], 'XXX', 'XXX', 'XXX', 'XXX', 'XXX', 'XXX', 'XXX') ) . "\n";
			print GLIST join("\t", ('XXX', 'XXX', 'XXX', 'XXX', 'XXX', @psl_info[9,8,13], 'XXX', 'XXX') ) . "\n";
		}
	}
	close S_GPSL;
	close O_GPSL;
	close GLIST;

	#-----------#
	# Trans PSL #
	#-----------#
	open (S_TPSL,Try_GZ_Read($V_Href->{virus_orig_trans_psl})) || die "fail read $V_Href->{virus_orig_trans_psl}: $!\n";
	open (O_TPSL,Try_GZ_Write($TransPSLforInteg)) || die "fail write $TransPSLforInteg: $!\n";
	while (<S_TPSL>){
		my @psl_info = split /\t+/;
		my $virus_refseg_Sref = \ $psl_info[13];
		# check whether this virus refseg has been acceptted in Gene PSL
		if( exists($yes_psl_refseg{$$virus_refseg_Sref}) ){
			# update the psl file refseg as the variant segname
			$$virus_refseg_Sref = $yes_psl_refseg{$$virus_refseg_Sref};
			# ouput to the Gene PSL for integration detection
			print O_TPSL join("\t", @psl_info); # last one has "\n"
		}
	}
	close S_TPSL;
	close O_TPSL;

	# inform
	stout_and_sterr `date`."[INFO]\tgenerate gene-psl, trans-psl, and gene-list for integration detection ok.\n";
}

#--- get makeup virus variant reference for integration detection ---
sub get_makeupRef_for_integration_detection{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $MakeupRef_Fa = $Option_Href->{MakeupRef_Fa};
	my $Variant_Extend_Href = $Option_Href->{Variant_Extend_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

	# deal each virus variant, later merge them together
	my @virus_variant_makeupFa;
	for my $variant_NO (sort {$a<=>$b} keys %$Variant_Info_Href){
		# virus variants selection filtration
		next unless( &is_this_virus_variants_been_selected({vNO=>$variant_NO}) );
		# basic information
		my $this_variant_makeupFa = $Variant_Info_Href->{$variant_NO}->{mkupFa};
		my $tmp_aim_fa = $MakeupRef_Fa.".vNO$variant_NO.fa";
		my $this_variant_ModSegname = $Variant_Info_Href->{$variant_NO}->{RefSeg};
		# extend ?
		my $extend_len = 0;
		if( exists($Variant_Extend_Href->{$this_variant_ModSegname}) ){
			$extend_len = $Variant_Extend_Href->{$this_variant_ModSegname}->{extend_len};
		}
		# generate tmp fa file
		FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->load_refpos_allel_from_ref_v_fa(
				{
					Ref => $this_variant_makeupFa,
					OutFaFile => $tmp_aim_fa,
					SegName	=> $this_variant_ModSegname,
					LineBase => $V_Href->{line_base},
					CircleExtLen => $extend_len
				}
			);
		# record
		push @virus_variant_makeupFa, $tmp_aim_fa;
	}

	# merge all tmp makeup fasta
	my $mergeFa_command = "cat @virus_variant_makeupFa > $MakeupRef_Fa";
	trible_run_for_success($mergeFa_command, 'MergeMakeupFA', {esdo_Nvb=>1, cmd_Nvb=>1});
	`rm -rf @virus_variant_makeupFa`;

	# create the indexed reference of merge virus variant makeup fasta
    BWA_index_fasta(FaFile => $MakeupRef_Fa, bwa => $V_Href->{bwa});
    Faidx_Dict_fasta(FaFile => $MakeupRef_Fa, samtools => $V_Href->{samtools});

	# inform
	stout_and_sterr `date`."[INFO]\tgenerate merged virus variant makeup fasta reference for integration detection ok.\n";
}

#--- prepare bam for virus integration detection ---
sub get_merge_bam_for_integration_detection{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $BamForInteg = $Option_Href->{BamForInteg};
	my $Variant_Extend_Href = $Option_Href->{Variant_Extend_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

	#------------------------------------#
	# merge all virus variants' bam and  #
	# select the best aligment for reads #
	# according to the selection, if set #
	#------------------------------------#
	my $VirusVariant_MergeBam = File::Spec->catfile($V_Href->{Prepare_Integ_folder},'VirusVariants.merge.bam');
	my @virus_variants_NO = grep &is_this_virus_variants_been_selected({vNO=>$_}), # only the selected virus variants
								 sort {$a<=>$b} keys %$Variant_Info_Href;
	my %previous_variants;
	my $last_bam_pref = 'null'; # initialize
	for my $i (0 .. $#virus_variants_NO){
		my $this_vNO = $virus_variants_NO[$i];
		my $this_bam_pref  = File::Spec->catfile($V_Href->{Prepare_Integ_folder},"VirusVariant.to_V$this_vNO.merge");
		if($i == 0){ # first one, just link
			`ln -sf $Variant_Info_Href->{$this_vNO}->{finalBam} $this_bam_pref.bam`;
		}
		else{
			## merge bams and sort by read name
			my $merge_nSortBam = FuseSV::Virus_Type::ReconstrcutVirusVariants->merge_bam_and_Nsort(
					{
						bam_pref => $this_bam_pref,
						source_bam_Aref => [ $Variant_Info_Href->{$this_vNO}->{finalBam} ],
						prev_bam_Aref => [ $last_bam_pref.'.bam' ]
					}
				);
			# select reads
			FuseSV::Virus_Type::ReconstrcutVirusVariants->filter_Nsort_bam_and_to_next_variant_bam(
					{
						nSortBam => $merge_nSortBam,
						OutBam => $this_bam_pref.'.bam',
						Prev_vref_Href => \%previous_variants,
						keep_prev_AlignInfo => 1
					}
				);
			# sweep
			unless( $V_Href->{in_debug} ){
				`rm -rf $last_bam_pref.bam`;
				`rm -rf $last_bam_pref.bam.bai`;
				`rm -rf $merge_nSortBam`;
			}
		}
		# load this virus variant into previous hash
		$previous_variants{ $Variant_Info_Href->{$this_vNO}->{RefSeg} } = 1;
		# change to the standard result
		if( $i == $#virus_variants_NO ){
			`mv -f $this_bam_pref.bam $VirusVariant_MergeBam`;
			`rm -rf $this_bam_pref.bam.bai`;
		}
		# update 'last' info
		$last_bam_pref = $this_bam_pref;
	}

	#--------------------------------#
	# merge the Nsort source bam and #
	# the virus-variant aligned bam  #
	#--------------------------------#
	my @RGs_Nsort_bam = map $V_Href->{RG_Info_Href}->{$_}->{Cand_nSortBam}, sort keys %{$V_Href->{RG_Info_Href}};
	my $MergeBam_Pref = File::Spec->catfile($V_Href->{Prepare_Integ_folder}, 'Source_VirusVar');
	## merge bam and sort by read name
	my $merge_nSortBam = FuseSV::Virus_Type::ReconstrcutVirusVariants->merge_bam_and_Nsort(
			{
				bam_pref => $MergeBam_Pref,
				source_bam_Aref => \@RGs_Nsort_bam, # here, to use host ref seg in front in BAM header
				prev_bam_Aref => [ $VirusVariant_MergeBam ]
			}
		);
	# sweep
	unless( $V_Href->{in_debug} ){
		`rm -rf $VirusVariant_MergeBam`;
	}

	#---------------------------------------#
	# select proper alignments of PE-reads  #
	# from source bam and virus-variant bam #
	#---------------------------------------#
	# select reads
	&select_reads_alignments_for_integration_bam(
		{
			nSortBam => $merge_nSortBam,
			OutBam => "$BamForInteg.for_mkdup.bam",
			vref_Href => \%previous_variants,
			RefSeg_Extend_Href => $Variant_Extend_Href,
		}
	);
	# sweep
	unless( $V_Href->{in_debug} ){
		`rm -rf $merge_nSortBam`;
	}

	# #--------------------------------------------#
	# # select best alignment for each reads that  #
	# # shared by source bam and virus-variant bam #
	# #--------------------------------------------#
	# # select reads
	# FuseSV::Virus_Type::ReconstrcutVirusVariants->filter_Nsort_bam_and_to_next_variant_bam(
	# 		{
	# 			nSortBam => $merge_nSortBam,
	# 			OutBam => "$BamForInteg.for_mkdup.bam",
	# 			Prev_vref_Href => $V_Href->{Host_Refseg_Href},
	# 			keep_prev_AlignInfo => 1,
	# 			dis_SuppleAlign => 1,
	# 			RefSeg_Extend_Href => $Variant_Extend_Href,
	# 			exchange_AlignEnds => 1,
	# 			mark_Prev_SamLine => 1 # this time, host is the 'Prev', so should mark any out of 'nonPrev'.
	# 		}
	# 	);
	# # sweep
	# unless( $V_Href->{in_debug} ){
	# 	`rm -rf $merge_nSortBam`;
	# }

	#---------------------------------------------#
	# mark duplications in merged BamForInteg bam #
	#---------------------------------------------#
	my $metrics = "$BamForInteg.mat";
	my $mkdup_tmpdir = "$BamForInteg.mkdup_tmpdir";
	`rm -rf $mkdup_tmpdir`;
	`mkdir -p $mkdup_tmpdir`;
	my $MkdupBam_command = "($V_Href->{java} -Xmx$V_Href->{mem_for_jvm} -jar $V_Href->{picard_jar} MarkDuplicates REMOVE_DUPLICATES=false I=$BamForInteg.for_mkdup.bam O=$BamForInteg METRICS_FILE=$metrics TMP_DIR=$mkdup_tmpdir 2>/dev/null) && ($V_Href->{samtools} index $BamForInteg)";
	trible_run_for_success($MkdupBam_command, 'BamForInteg_MkdupBam', {esdo_Nvb=>1});
	# sweep
	unless( $V_Href->{in_debug} ){
		`rm -rf $metrics $mkdup_tmpdir`;
		`rm -rf $BamForInteg.for_mkdup.bam $BamForInteg.for_mkdup.bam.bai`;
	}

	#------------------------------#
	# simplify the BamForInteg bam #
	# discard duplicated PE-reads  #
	#------------------------------#
	if( $V_Href->{Get_simple_BamForInteg} ){
		my $unsimple_bam = "$BamForInteg.unsimple.bam";
		`mv $BamForInteg $unsimple_bam`;
		# discard duplicated reads
		my $DisPEdup_cmd = "($V_Href->{samtools} view -h -F 0x900 $unsimple_bam | awk '/^\@/ || !and(\$2,0x400) || and(\$2,0xC)' | $V_Href->{samtools} view -b -S -o $BamForInteg) && ($V_Href->{samtools} index $BamForInteg)";
		trible_run_for_success($DisPEdup_cmd, 'DisPEdupForSimpleBam', {esdo_Nvb=>1});
		# &simplify_bam_for_integration(
		# 	{
		# 		input_bam => $unsimple_bam,
		# 		output_bam => $BamForInteg,
		# 		vNO_Aref => \@virus_variants_NO
		# 	}
		# );
		# sweep
		unless( $V_Href->{in_debug} ){
			`rm -rf $unsimple_bam`;
		}
	}

	# inform
	stout_and_sterr `date`."[INFO]\tgenerate merged bam for integration detection ok.\n";
}

#--- select proper alignments of PE-reads for integration detection ---
sub select_reads_alignments_for_integration_bam{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Nsort_bam = $Option_Href->{nSortBam};
	my $OutBam = $Option_Href->{OutBam};
	my $v_refseg_Href = $Option_Href->{vref_Href};
	my $RefSeg_Extend_Href = $Option_Href->{RefSeg_Extend_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $max_rlen = $V_Href->{max_rlen} || 0;

	open (my $OBAM_FH, "| $V_Href->{samtools} view -b -S -o $OutBam.unsort.bam -") || die "fail write $OutBam.unsort.bam: $!\n";

	# filter bam
	my $last_pid = 'initial-value'; # initial to start
	my %pid_InfoPool;
	my $amount_NeedToDeal = 0;
	#   0x800, to avoid supplementary alignment
	# + 0x100: to avoid secondary alignment
	# = 0x900
	open (NSBAM, "$V_Href->{samtools} view -F 0x900 -h $Nsort_bam |") || die "fail read $Nsort_bam: $!\n";
	while (<NSBAM>){
		# header
		if(/^\@/){
			if( !/^\@SQ/ ){
				print $OBAM_FH $_;
			}
			else{
				my ($Refseg, $Reflength) = (/\sSN:(\S+)\sLN:(\d+)/);
				# to extend the length?
				if(exists( $RefSeg_Extend_Href->{$Refseg} )){
					$Reflength += $RefSeg_Extend_Href->{$Refseg}->{extend_len} - $max_rlen;
					print $OBAM_FH "\@SQ\tSN:$Refseg\tLN:$Reflength\n";
				}
				else{
					print $OBAM_FH $_;
				}
			}
			next;
		}
		# line info
		my @lineInfo = split;
		# correct the mapQ, especially for Q255 induced by alignment tool STAR/HISAT2.
		$lineInfo[4] = min( $lineInfo[4], 60 );
		# each attribute
		my ($pid, $flag, $refseg, $map_pos, $mapQ, $Cigar, $Tlen) = @lineInfo[0,1,2,3,4,5,8];
		# new pid to deal, count it
		if( $pid ne $last_pid ){
			# meet the store threshold
			# deal groups of new-pids
			if( ++$amount_NeedToDeal >= $V_Href->{Nsort_MaxStore}){
				&choose_vv_source_alignment_to_output(
					{
						file_handle => $OBAM_FH,
						pid_InfoPool_Href => \%pid_InfoPool
					}
				);
				# reset
				$amount_NeedToDeal = 0;
			}
			# become the last pid
			$last_pid = $pid;
		}
		# record
		my $end_NO = ($flag & 0x40) ? 1 : 2; # $flag & 0x80 + 1
		my ($align_score) = (/\bAS:i:(\d+)\b/); # even the read is unmap(0x4), this line still contain AS:i:0
		my $alignInfo_Href = {
				LineInfo => join("\t",@lineInfo)."\n",
				MapQ => $mapQ # MapQ of unmapped read is 0
			};
		# store in the InfoPool
		# note that, if one is mapped, another is unmapped, the unmapped end's refseg must be same as that of the mapped one.
		my $info_tag = (exists $v_refseg_Href->{$refseg}) ? 'V' : 'H';
		$info_tag .= ($flag & 0x4) ? 'U' : 'M';
		$pid_InfoPool{$pid}->{$end_NO}->{$info_tag} = $alignInfo_Href;
	}
	close NSBAM;
	# deal the last pid and the last group
	&choose_vv_source_alignment_to_output(
		{
			file_handle => $OBAM_FH,
			pid_InfoPool_Href => \%pid_InfoPool
		}
	);
	close $OBAM_FH;

	# sort bam
	my $SortBam_command = "$V_Href->{samtools} sort --threads 3 -m 1G -T $OutBam.sort_tmp -o $OutBam $OutBam.unsort.bam && $V_Href->{samtools} index $OutBam";
	trible_run_for_success($SortBam_command, 'SortBam', {cmd_Nvb=>1, esdo_Nvb=>1});
	`rm -rf $OutBam.unsort.bam`;
}

#--- choose alignments from virus variants and source ---
# In_VvBam     Action
# _1:M;_2:M    exchange at most one end.
# _1:M;_2:U    only exchange U end.
# _1:U;_2:M    only exchange U end.
# _1:U;_2:U    not exists, just output source bam (non-Simple-mode).
sub choose_vv_source_alignment_to_output{
	# options
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $OBAM_FH = $Option_Href->{file_handle};
	my $pid_InfoPool_Href = $Option_Href->{pid_InfoPool_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	for my $pid (keys %$pid_InfoPool_Href){
		my $end_1_Href = $pid_InfoPool_Href->{$pid}->{1};
		my $end_2_Href = $pid_InfoPool_Href->{$pid}->{2};
		my ($end_1_tag, $end_2_tag);
		# _1:VM;_2:VM
		if( exists $end_1_Href->{VM} && exists $end_2_Href->{VM} ){
			my $end_1_HmapQ_offset = (exists $end_1_Href->{HM}) ? ($end_1_Href->{HM}->{MapQ} - $end_1_Href->{VM}->{MapQ}) : -10000;
			my $end_2_HmapQ_offset = (exists $end_2_Href->{HM}) ? ($end_2_Href->{HM}->{MapQ} - $end_2_Href->{VM}->{MapQ}) : -10000;
			# pick the larger optimized
			if( $end_1_HmapQ_offset >= $end_2_HmapQ_offset){
				$end_1_tag = (exists $end_1_Href->{HM}) ? 'HM' : 'VM';
				$end_2_tag = 'VM';
			}
			else{
				$end_1_tag = 'VM';
				$end_2_tag = (exists $end_2_Href->{HM}) ? 'HM' : 'VM';
			}
		}
		# _1:VM;_2:VU
		elsif( exists $end_1_Href->{VM} && exists $end_2_Href->{VU} ){
			$end_1_tag = 'VM';
			$end_2_tag = (exists $end_2_Href->{HM}) ? 'HM' : 'VU';
		}
		# _1:VU;_2:VM
		elsif( exists $end_1_Href->{VU} && exists $end_2_Href->{VM} ){
			$end_1_tag = (exists $end_1_Href->{HM}) ? 'HM' : 'VU';
			$end_2_tag = 'VM';
		}
		# _1:VU;_2:VU, should only have 'H' prefixed tag.
		else{
			unless( $V_Href->{Get_simple_BamForInteg} ){
				$end_1_tag = (exists $end_1_Href->{HM}) ? 'HM' : 'HU';
				$end_2_tag = (exists $end_2_Href->{HM}) ? 'HM' : 'HU';
			}
		}
		# output
		if( $end_1_tag && $end_2_tag ){
			my $Info_Href = { 1=>$end_1_Href->{$end_1_tag}, 2=>$end_2_Href->{$end_2_tag} };
			# different segment, so should update line info
			if( substr($end_1_tag,0,1) ne substr($end_2_tag,0,1) ){
				for my $Aref ( [1,2], [2,1] ){
					my $thisEnd_Href = $Info_Href->{$Aref->[0]};
					my $pairEnd_Href = $Info_Href->{$Aref->[1]};
					my @thisEnd_info = split /\t+/, $thisEnd_Href->{LineInfo};
					my @pairEnd_info = split /\t+/, $pairEnd_Href->{LineInfo};
					# update
					## mated end mapped refseg and position
					$thisEnd_info[6] = $pairEnd_info[2];
					$thisEnd_info[7] = $pairEnd_info[3];
					## set Tlen as zero
					$thisEnd_info[8] = 0;
					## flags
					$thisEnd_info[1] += 0x1 unless( $thisEnd_info[1] & 0x1 );
					$thisEnd_info[1] -= 0x2     if( $thisEnd_info[1] & 0x2 );
					$thisEnd_info[1] -= 0x4     if( $thisEnd_info[1] & 0x4 );
					$thisEnd_info[1] -= 0x8     if( $thisEnd_info[1] & 0x8 );
					$thisEnd_info[1] += 0x20    if(   $pairEnd_info[1] & 0x10  && !($thisEnd_info[1] & 0x20) );
					$thisEnd_info[1] -= 0x20    if( !($pairEnd_info[1] & 0x10) &&   $thisEnd_info[1] & 0x20  );
					$thisEnd_info[1] -= 0x100   if( $thisEnd_info[1] & 0x100 );
					$thisEnd_info[1] -= 0x200   if( $thisEnd_info[1] & 0x200 );
					$thisEnd_info[1] -= 0x400   if( $thisEnd_info[1] & 0x400 );
					$thisEnd_info[1] -= 0x800   if( $thisEnd_info[1] & 0x800 );
					## mapQ
					my $thisEnd_info_updated = join("\t", @thisEnd_info); # last element has "\n"
					my $pair_mapQ = $pairEnd_Href->{MapQ};
					$thisEnd_info_updated =~ s/\bMQ:i:\d+\b/MQ:i:$pair_mapQ/;
					# reassign LineInfo of thisEnd
					$thisEnd_Href->{LineInfo} = $thisEnd_info_updated;
				}
			}
			# output
			print $OBAM_FH $Info_Href->{$_}->{LineInfo} for (1, 2);
		}
		# sweep
		delete $pid_InfoPool_Href->{$pid};
	}

	# reset
	%$pid_InfoPool_Href = ();
}

#--- simplify the bam, filter reads ---
## duplicated PE-reads 0x400 but not 0x8
sub simplify_bam_for_integration{

	# my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	# my $input_bam  = $Option_Href->{input_bam};
	# my $output_bam = $Option_Href->{output_bam};
	# my $vNO_Aref   = $Option_Href->{vNO_Aref};

	# # variants
	# my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

	# my %vSegName = map { ( $Variant_Info_Href->{$_}->{RefSeg}, 1 ) } @$vNO_Aref;

	# open (SMBAM, "| $V_Href->{samtools} view -b -S -o $output_bam") || die "fail write simple bam: $!\n";
	# open (USMBAM, "$V_Href->{samtools} view -h -F 0x900 $input_bam |") || die "fail read unsimple bam: $!\n";
	# while(<USMBAM>){
	# 	if(/^\@/){ # SAM header
	# 		print SMBAM;
	# 		next;
	# 	}
	# 	# read info
	# 	my ($flag, $self_Mrefseg, $Cigar, $mate_Mrefseg) = (split)[1,2,5,6];
	# 	# filter reads
	# 	if(	$flag & 0x400 && !($flag & 0x8) ){
	# 		next;
	# 	}
	# 	else{
	# 		print SMBAM;
	# 	}
	# }
	# close USMBAM;
	# close SMBAM;

	# # command
	# my $smbam_index_command = "$V_Href->{samtools} index $output_bam";
	# trible_run_for_success($smbam_index_command, 'SimpleBamIndex', {esdo_Nvb=>1});
}

#--- get parameter txt doc for integration detection ---
sub get_parameter_for_integration_detection{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $ParaTxt = $Option_Href->{ParaTxt};
	my $Variant_Extend_Href = $Option_Href->{Variant_Extend_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Prepare_For_Integ_Href = $V_Href->{Prepare_For_Integ};

	# generate
	open (PARA,Try_GZ_Write($ParaTxt)) || die "fail write $ParaTxt: $!\n";
	print PARA '-fgl '  . $Prepare_For_Integ_Href->{GeneList} . ' ';
	print PARA '-tbm '  . $V_Href->{max_rlen} . ',' . $Prepare_For_Integ_Href->{BamForInteg} . ' ';
	print PARA '-vfa '  . $Prepare_For_Integ_Href->{MakeupRef} . ' ';
	print PARA '-gpsl ' . $Prepare_For_Integ_Href->{GenePSL} . ' ';
	print PARA '-tpsl ' . $Prepare_For_Integ_Href->{TransPSL} . ' ';
	print PARA ' -cvel ' . "'$_=$Variant_Extend_Href->{$_}->{extend_len}'"
			for sort { $Variant_Extend_Href->{$a}->{variant_NO} <=> $Variant_Extend_Href->{$b}->{variant_NO} } 
				keys %$Variant_Extend_Href;
	print PARA "\n";
	close PARA;

	# inform
	stout_and_sterr `date`."[INFO]\tgenerate suggested FuseSV parameter for integration detection ok.\n";
}

#--- virus variants selection filtration ---
sub is_this_virus_variants_been_selected{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $virus_vNO = $Option_Href->{vNO};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Prepare_For_Integ_Href = $V_Href->{Prepare_For_Integ};
	my $Variant_Info_Href = $V_Href->{Variant_Info_Href};

	# check whether this vNO is been selected
	# if non-exist 'AllowVNO', means all are selected.
	if( exists($Prepare_For_Integ_Href->{AllowVNO}) ){
		if( exists($Prepare_For_Integ_Href->{AllowVNO}->{$virus_vNO}) ){
			if(    $virus_vNO == 1
				|| (    $Variant_Info_Href->{$virus_vNO}->{AverageDepth} >= $V_Href->{MinDepth_For_IntegAna}
					 && &is_this_virus_variants_Nstr_OK( { vfa=>$Variant_Info_Href->{$virus_vNO}->{Fa} } )
					 && $Variant_Info_Href->{$virus_vNO}->{CoverRatio} >= $V_Href->{MinCoverRatio_For_IntegAna}  )
			){
				return 1;
			}
			else{
				# sweep from candidate vNO for no-coming again
				delete $Prepare_For_Integ_Href->{AllowVNO}->{$virus_vNO};
				# inform
				stout_and_sterr `date`."<WARN>:\tDiscard virus variant NO.$virus_vNO as low-depth or long-Nstr.\n";
				return 0;
			}
		}
		else{
			return 0;
		}
	}
	else{
		return 1;
	}
}

#--- virus variants Nstr count ---
sub is_this_virus_variants_Nstr_OK{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $vfa = $Option_Href->{vfa};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# shut down this operation
	if( $V_Href->{Nstr_MinLength} < 0 ){
		return 1;
	}

	# count the Nstr exceeds 'Nstr_MinLength'
	(my $fa_str = `cat $vfa | grep -v '^>'`) =~ s/\s//g;
	my $Nstr_Count = grep { length($_) >= $V_Href->{Nstr_MinLength} } (split /[^Nn]+/, $fa_str);
	# return
	if( $Nstr_Count > $V_Href->{MaxAllow_Nstr_Count} ){
		return 0;
	}
	else{
		return 1;
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
