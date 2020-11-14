package FuseSV::Virus_Integ::DetectInteg;

use strict;
use warnings;
use Getopt::Long;
use File::Spec qw/catfile/;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use FuseSV::LoadOn;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use FuseSV::Objects::Virus_OB;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  detect_virus_integration
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'DetectInteg';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-05-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						detect_virus_integration
						evaluate_ReadGroup_Ins
						prepare_vInteg_workspace
						write_viteg_work_tag
						get_sampleID_from_virusDir
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} virus_integ <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -o_pref [s]  file name prefix of all output results. <required>
	                  only allows letter (A-Z and a-z), number (0-9) and under/short-line ('-' and '_').
	    -casevd [s]  out directory used by \033[4m \bvirus_type \b\033[0m function for \033[4m \bCase \b\033[0m. e.g., tumor sample. <required>
	                  all results of virus integration will be stored in sub-folder of \033[4m \bCase \b\033[0m's directory.
	    -ctrlvd [s]  out directory used by \033[4m \bvirus_type \b\033[0m function for \033[4m \bControl \b\033[0m. e.g., normal sample.
	    -casebm [s]  the original source indexed-bam file for \033[4m \bCase \b\033[0m. i.e., the '-bam' option. <required>
	    -ctrlbm [s]  the original source indexed-bam file for \033[4m \bControl \b\033[0m. effective with '-ctrlvd'.
	    -maxrdl [i]  the maximum read length. <required>

	   # Database #
	    -hostg  [s]  the host genome fasta file, indexed by BWA. <required>

	   # Software Required #
	    -bwa    [s]  BWA tool, minimum version: 0.7.9. <required>
	    -stl    [s]  SamTools, minimum version: 1.3. <required>
	    -blat   [s]  Blat tool. <required>

	   # Detection Options #
	    -flank  [i]  the flanking bases to extend chromosome region of anchor reads. [1000]
	    -mincsp [i]  the minimum deduplicated span-reads allowed for one candidate integration cluster. [3]
	    -dpmode [i]  mode of de-deuplication of supporting reads from same library, set as 1 or 2. [2]
	                  set 1 for position-based; set 2 for sequence-based.
	    -minjr  [i]  the minimum deduplicated junc-reads allowed for one viral integration. [3]
	    -minjsl [i]  the minimum length spanned the integration break-point by one junc-reads. [5]
	    -minjbl [i]  the minimum bilateral length covered by all junc-reads for one viral integration. [10]
	    -no_tdchk    Disable the check of tandem repeats in the region covered by junc-reads. [enabled]

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
			[ vItg_OutPrefix => '' ],
			[ virusAna_dir => { case => '', ctrl => '' }  ],
			[ SourceBam => {} ], # tissue => bam_OB (source)
			[ MaxReadLength => 0 ],
			# database
			[ host_genome_Idx_fa => '' ],
			[ dVG_viral_BkList => '' ],
			[ dVG_show_no_gene => 0 ],
			# software
			[ bwa => '' ],
			[ samtools => '' ],
			[ blat => '' ],
			# filter
			[ MinSpanPairForCand => 3 ],
			[ MaxSpanPairForCand => 1E10 ],
			[ SupReadsDeDupMode => 2 ],
			[ MinJrCount => 3 ],
			[ MinLenJrSpanBreakage => 5 ],
			[ MinLenAllJrSpanBreakageBilateral => 10 ],
			[ DisableTandemCheckOfJrCoverSeg => 0 ],
			# other
			[ MaxMemoryForForkNumber => '4G' ],
			[ SkipInsEvalue => 0 ],
			[ StandNumPairForInsEvalue => 1E5 ],
			[ MinPairForInsEvalue => 1000 ],
			[ PRdropProbForInsEvalue => 0 ],
			[ RG3pOverlapPerctCutoff => 0.2 ],
			[ Max_InsMeanAddNsd => 0 ], # ins_mean + {Max_InsSD_Ntimes} * ins_sd
			[ Max_InsSD_Ntimes => 3 ],
			[ FlankSizeExtendHostRegion => 1000 ],
			# annotation
			[ gene_psl => '' ],
			[ trans_psl => '' ],
			[ cytoBand_file => '' ],
			[ ClusterSizeForMajorMinorCase => 2E6 ], # 2M
			[ FlankSizeForMHcheck => 100 ],
			[ MinUnitOfMH => 2 ],
			[ BiSegExtractAtMH => 20 ],
			# intermediate variants
			[ vIteg_SampleID => { case => '', ctrl => '' } ],
			[ vIteg_tagFile => '' ],
			[ VirusBam => {} ], # tissue => bam_OB (virus)
			[ vIteg_rgid2rgOB => {} ], # RG_ID => rg_OB
			[ vIteg_rIDpref2rgOB => {} ], # rid-prefix => rg_OB
			[ vIteg_Virus => {} ], # virusSegName => virus_OB

			# list to abs-path
			[ ToAbsPath_Aref => [ ['bwa'],
								  ['samtools'],
								  ['blat'],
								  ['host_genome_Idx_fa'],
						  		  ['virusAna_dir', 'case'],
						  		  ['virusAna_dir', 'ctrl'],
								  ['SourceBam', 'case', 'filepath'],
								  ['SourceBam', 'ctrl', 'filepath']  ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# objects constructions, as these objects belong to the input options
	$V_Href->{SourceBam}->{case} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'case', filepath => undef );
	$V_Href->{SourceBam}->{ctrl} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'ctrl', filepath => undef );

	# get options
	GetOptions(
		# input/output
		"-o_pref:s"	=> \$V_Href->{vItg_OutPrefix},
		"-casevd:s"	=> \$V_Href->{virusAna_dir}->{case},
		"-ctrlvd:s"	=> \$V_Href->{virusAna_dir}->{ctrl},
		"-casebm:s"	=> \$V_Href->{SourceBam}->{case}->{filepath},
		"-ctrlbm:s"	=> \$V_Href->{SourceBam}->{ctrl}->{filepath},
		"-maxrdl:i"	=> \$V_Href->{MaxReadLength},
		# database
		"-hostg=s"	=> \$V_Href->{host_genome_Idx_fa},
		# "-virusg=s"	=> \$V_Href->{virus_genome_Idx_fa},
		# "-cir_vel=s"=> \%{$V_Href->{CircleVirusExtLen}},
		# software
		"-bwa:s"	=> \$V_Href->{bwa},
		"-stl:s"	=> \$V_Href->{samtools},
		"-blat:s"	=> \$V_Href->{blat},
		# options
		## multiple fork operations
		"-maxmem:i"	=> \$V_Href->{MaxMemoryForForkNumber}, # hidden option
		## insert size evaluation
		"-skip_eins"=> \$V_Href->{SkipInsEvalue}, # hidden option
		"-mp_eins:i"=> \$V_Href->{MinPairForInsEvalue}, # hidden option
		"-pd_eins:f"=> \$V_Href->{PRdropProbForInsEvalue}, # hidden option
		## host genomic region selection
		"-flank:i"	=> \$V_Href->{FlankSizeExtendHostRegion},
		## filteration
		"-mincsp:i"	=> \$V_Href->{MinSpanPairForCand},
		"-maxcsp:i"	=> \$V_Href->{MaxSpanPairForCand}, # hidden option
		"-dpmode=i"	=> \$V_Href->{SupReadsDeDupMode},
		"-minjr:i"	=> \$V_Href->{MinJrCount},
		"-minjsl:i"	=> \$V_Href->{MinLenJrSpanBreakage},
		"-minjbl:i"	=> \$V_Href->{MinLenAllJrSpanBreakageBilateral},
		"-no_tdchk"	=> \$V_Href->{DisableTandemCheckOfJrCoverSeg},
		## annotation
		# "-gpsl=s"	=> \$V_Href->{gene_psl},
		# "-tpsl=s"	=> \$V_Href->{trans_psl},
		# "-cbd=s"	=> \$V_Href->{cytoBand_file},
		# "-mjclut=i"	=> \$V_Href->{ClusterSizeForMajorMinorCase},
		# "-mh_flk=i"	=> \$V_Href->{FlankSizeForMHcheck},
		# "-mh_minu=i"=> \$V_Href->{MinUnitOfMH},
		# "-mh_blen=i"=> \$V_Href->{BiSegExtractAtMH},
		# help
		"-h|help"	=> \$V_Href->{HELP},
		# for debug
		"-debug=f"	=> \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || (!$V_Href->{vItg_OutPrefix} || $V_Href->{vItg_OutPrefix} =~ /[^\w\-]/)
			 || (!$V_Href->{virusAna_dir}->{case} || !-d $V_Href->{virusAna_dir}->{case})
			 || (!-e $V_Href->{SourceBam}->{case}->{filepath})
			 || (-d $V_Href->{virusAna_dir}->{ctrl} xor -e $V_Href->{SourceBam}->{ctrl}->{filepath})
			 || $V_Href->{MaxReadLength} < 50
			 || (!$V_Href->{host_genome_Idx_fa} || !-e $V_Href->{host_genome_Idx_fa})
			 || (!$V_Href->{bwa} || !-e $V_Href->{bwa})
			 || (!$V_Href->{samtools} || !-e $V_Href->{samtools})
			 || (!$V_Href->{blat} || !-e $V_Href->{blat})
			 || $V_Href->{MinSpanPairForCand} <= 0
			 || $V_Href->{MaxSpanPairForCand} <= 0
			 || $V_Href->{SupReadsDeDupMode} !~ /^[12]$/
			 || $V_Href->{MinLenJrSpanBreakage} <= 0
			 || $V_Href->{MinLenAllJrSpanBreakageBilateral} <= 0
			);
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub detect_virus_integration{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# prepare virus integration workspace
	&prepare_vInteg_workspace;

	# evaluate insertsize of reads groups
	&evaluate_ReadGroup_Ins;

	# load Case sample's virus
	&load_virus_for_integration;
}

#--- load virus of Case sample for viral integrations ---
sub load_virus_for_integration{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $virusIdxFa = File::Spec->catfile( $V_Href->{virusAna_dir}->{case},
										  "$V_Href->{vIteg_SampleID}->{case}.$V_Href->{Prepare_Integ_folder}",
										  "$V_Href->{ForVitgPrefix}.virus_makeup_ref",
										  "$V_Href->{ForVitgPrefix}.virus.makeup.fa" );
	my $virusList = File::Spec->catfile( $V_Href->{virusAna_dir}->{case}, "$V_Href->{vIteg_SampleID}->{case}.$V_Href->{Prepare_Integ_folder}", "$V_Href->{ForVitgPrefix}.virus.list");
	open (VL, $virusList) || die "fail read virus list: $!\n";
	while(<VL>){
		my ($virusSegName, $virusExtLen) = (split);
		# instance: T0001-HPV16-RC
		my ($SampleID, $virusType) = ($virusSegName =~ /(.+)-([^\-]+)-RC/);
		if( $SampleID ne $V_Href->{vIteg_SampleID}->{case} ){
			warn_and_exit "<ERROR>\tDifferent case sample-ID detected from virusSegName $virusSegName\n";
		}
		# create virus_OB
		$V_Href->{vIteg_Virus}->{$virusSegName} = FuseSV::Objects::Virus_OB->new( SegName => $virusSegName, IdxFa => $virusIdxFa, subType => $virusType, extLen => $virusExtLen );
	}
	close VL;
}

#--- evaluate insertsize ---
sub evaluate_ReadGroup_Ins{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $rgid2rgOB_Href    = $V_Href->{vIteg_rgid2rgOB};
	my $rIDpref2rgOB_Href = $V_Href->{vIteg_rIDpref2rgOB};
	my $pr_AimCount = $V_Href->{StandNumPairForInsEvalue};
	my $Tool_Tag = $V_Href->{Tool_Tag};
	my $samtools = $V_Href->{samtools};

	# work-folder
	my $RGinsStoreDir = File::Spec->catfile($V_Href->{vIteg_workFolder}, "$V_Href->{vItg_OutPrefix}.evaluate_ins");
	`rm -rf $RGinsStoreDir && mkdir -p $RGinsStoreDir` unless( $V_Href->{SkipInsEvalue} );

	for my $tissue ('case', 'ctrl'){
		next unless( defined $V_Href->{vIteg_SampleID}->{$tissue} );
		# loadup For_Detect_Integ bam objects
		# $V_Href->{Prepare_Integ_folder} is initialized as only the postfix.
		my $VirusBamPath = File::Spec->catfile( $V_Href->{virusAna_dir}->{$tissue}, "$V_Href->{vIteg_SampleID}->{$tissue}.$V_Href->{Prepare_Integ_folder}", "$V_Href->{ForVitgPrefix}.bam");
		$V_Href->{VirusBam}->{$tissue} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => $tissue, filepath => $VirusBamPath );
		# load reads group from virus bam, as RG info of 'virus bam' is officially designed
		$V_Href->{VirusBam}->{$tissue}->pick_rgOB( rgid2rgOB_Href => $rgid2rgOB_Href, samtools => $samtools ); # $V_Href->{Allow_Lack_RG}
		# load previous or calculate
		if( $V_Href->{SkipInsEvalue} ){
			# find and load rg_report
			my $rg_report_glob = File::Spec->catfile($RGinsStoreDir, "$tissue.*.report");
			for my $rg_report ( sort glob($rg_report_glob) ){
				chomp(my $RG_ID = `cat $rg_report | grep '^RG_ID' | cut -f2`);
				next unless( exists $rgid2rgOB_Href->{$RG_ID} );
				my $rg_OB = $rgid2rgOB_Href->{$RG_ID};
				$rg_OB->load_report( rg_report => $rg_report );
				# inform
				stout_and_sterr `date`."[INFO]\tLoad report of $tissue reads group $RG_ID.\n";
			}
			# check
			for my $RG_ID ( sort keys %$rgid2rgOB_Href ){
				my $rg_OB = $rgid2rgOB_Href->{$RG_ID};
				next if( $rg_OB->{tissue} ne $tissue);
				unless( defined $rg_OB->{rID_prefix} ){
					warn_and_exit "<ERROR>\tMiss report of $tissue reads group $RG_ID.\n";
				}
			}
		}
		else{
			# link the RG of virus bam to source bam
			$V_Href->{SourceBam}->{$tissue}->{rgOB} = $V_Href->{VirusBam}->{$tissue}->{rgOB};
			# load reads for reads group
			for my $opt_Aref (  [ 'VirusBam',  0, 0 ],
								[ 'SourceBam', 0, $V_Href->{PRdropProbForInsEvalue} ],
								[ 'VirusBam',  1, 0 ],
								[ 'SourceBam', 1, $V_Href->{PRdropProbForInsEvalue} ]
			){
				my $bam_OB = $V_Href->{$opt_Aref->[0]}->{$tissue};
				$bam_OB->load_reads_for_ReadsGroup( rgid2rgOB_Href => $rgid2rgOB_Href, only_SoftClip => $opt_Aref->[1], pr_DropProb => $opt_Aref->[2],
													pr_AimCount => $pr_AimCount, Tool_Tag => $Tool_Tag, samtools => $samtools );
			}
			# ins evalue of RG
			for my $RG_ID ( sort keys %$rgid2rgOB_Href ){
				my $rg_OB = $rgid2rgOB_Href->{$RG_ID};
				next if( $rg_OB->{tissue} ne $tissue);
				# ins evaluation
				$rg_OB->evalue_ins( MinPairForInsEvalue => $V_Href->{MinPairForInsEvalue} );
				# 3 prime overlap: PCW
				$rg_OB->test_3p_overlap( perct_cutoff => $V_Href->{RG3pOverlapPerctCutoff} );
				# output stat files
				$rg_OB->write_report( folder => $RGinsStoreDir );
				$rg_OB->write_insDistLog( folder => $RGinsStoreDir );
				# inform
				stout_and_sterr `date`."[INFO]\tOutput info of $tissue reads group $RG_ID.\n";
			}
		}
	}

	# others
	for my $RG_ID ( sort keys %$rgid2rgOB_Href ){
		my $rg_OB = $rgid2rgOB_Href->{$RG_ID};
		# map readid-prefix to rg_OB
		$rIDpref2rgOB_Href->{ $rg_OB->{rID_prefix} } = $rg_OB;
		# max ( ins + 3 * sd )
		$V_Href->{Max_InsMeanAddNsd} = max( $V_Href->{Max_InsMeanAddNsd}, $rg_OB->{Ins_mean} + $V_Href->{Max_InsSD_Ntimes} * $rg_OB->{Ins_sd});
	}
}

#--- prepare workspace for virus integration ---
sub prepare_vInteg_workspace{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get ID of Case sample
	my $CaseID_Aref = &get_sampleID_from_virusDir( virusDir => $V_Href->{virusAna_dir}->{case} );
	if( $CaseID_Aref->[1] == 0 ){
		warn_and_exit "<ERROR>\tCannot find $V_Href->{Prepare_Integ_folder} folder under Case sample's virus analysis directory.\n";
	}
	else{
		$V_Href->{vIteg_SampleID}->{case} = $CaseID_Aref->[0];
	}

	# get ID of Control sample, if set
	if( -d $V_Href->{virusAna_dir}->{ctrl} ){
		my $ControlID_Aref = &get_sampleID_from_virusDir( virusDir => $V_Href->{virusAna_dir}->{ctrl} );
		if( $ControlID_Aref->[1] == 0 ){
			warn "<WARN>\tCannot find $V_Href->{Prepare_Integ_folder} folder under Control sample's virus analysis directory.\n".
					   "\tNote that FuseSV will ignore all Control data.\n";
		}
		$V_Href->{vIteg_SampleID}->{ctrl} = $ControlID_Aref->[0];
	}
	else{
		$V_Href->{vIteg_SampleID}->{ctrl} = undef;
	}

	# create folder for virus integration
	$V_Href->{vIteg_workFolder} = File::Spec->catfile($V_Href->{virusAna_dir}->{case}, "$V_Href->{vIteg_SampleID}->{case}.$V_Href->{vIteg_workFolder}");
	`mkdir -p $V_Href->{vIteg_workFolder}`;
	# inform
	stout_and_sterr `date`."[INFO]\tWork Folder for virus integration has been created:\n".
								 "\t$V_Href->{vIteg_workFolder}\n";

	# write work tag file
	&write_viteg_work_tag;
}

# write work tag file
sub write_viteg_work_tag{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	$V_Href->{vIteg_tagFile} = File::Spec->catfile($V_Href->{vIteg_workFolder}, "$V_Href->{vItg_OutPrefix}.virus_integ.tag");

	# structure of report
	my @Structure = (
						[ 'Time',     chomp(my $date = `date`)                              ],
						[ 'O_pref',   $V_Href->{vItg_OutPrefix}                             ],
						[ 'CaseVD',   $V_Href->{virusAna_dir}->{case}                      ],
						[ 'CaseID',   $V_Href->{vIteg_SampleID}->{case}                     ],
						[ 'CaseSBM',  $V_Href->{SourceBam}->{case}->{filepath}              ],
						[ 'HasCtrl',  -d $V_Href->{virusAna_dir}->{ctrl} ? 'Y' : 'N'       ],
						[ 'CtrlVD',   $V_Href->{virusAna_dir}->{ctrl} || 'N/A'             ],
						[ 'CtrlID',   $V_Href->{vIteg_SampleID}->{ctrl} || 'N/A'            ],
						[ 'CtrlPVI',  defined $V_Href->{vIteg_SampleID}->{ctrl} ? 'Y' : 'N' ],
						[ 'CtrlSBM',  $V_Href->{SourceBam}->{ctrl}->{filepath} || 'N/A'     ],
						[ 'MaxRlen',  $V_Href->{MaxReadLength}                              ],
						[ 'EvalIns',  $V_Href->{SkipInsEvalue} ? 'N' : 'Y'                  ],
						[ 'PrdpIns',  $V_Href->{PRdropProbForInsEvalue}                     ],
						[ 'HostRef',  $V_Href->{host_genome_Idx_fa}                         ],
						[ 'GRflank',  $V_Href->{FlankSizeExtendHostRegion}                  ],
						[ 'MinCSrP',  $V_Href->{MinSpanPairForCand}                         ],
						[ 'DPmode',   $V_Href->{SupReadsDeDupMode}                          ],
						[ 'MinJr',    $V_Href->{MinJrCount}                                 ],
						[ 'MinJSL',   $V_Href->{MinLenJrSpanBreakage}                       ],
						[ 'MinJBL',   $V_Href->{MinLenAllJrSpanBreakageBilateral}           ],
						[ 'TDRchk',   $V_Href->{DisableTandemCheckOfJrCoverSeg} ? 'N' : 'Y' ]
					);

	open (TAG, Try_GZ_Write($V_Href->{vIteg_tagFile})) || die "fail create virus integration tag file: $!.\n";
	print TAG join("\t", @$_)."\n" for @Structure;
	close TAG;
	# inform
	stout_and_sterr `date`."[INFO]\tTag file generated.\n";
}

#--- get sample ID from $V_Href->{Prepare_Integ_folder} dirname
sub get_sampleID_from_virusDir{
	
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $virusDir = $parm{virusDir};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	opendir (VDIR, $virusDir) || die "fail to read virusDir: $!\n";
	for my $sub_file (readdir VDIR){
		if( $sub_file =~ /^(.+)\.$V_Href->{Prepare_Integ_folder}/ ){
			closedir VDIR;
			return [$1, 1];
		}
	}
	closedir VDIR;
	return [undef, 0];
}









































#--- 
1; ## tell the perl script the successful access of this module.
