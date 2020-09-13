package FuseSV::Extension::GetCtrlVCFfromCase;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum first /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr trible_run_for_success Get_Two_Seg_Olen /;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  get_CtrlVCF_from_CaseVCF
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'GetCtrlVCFfromCase';
#----- version --------
$VERSION = "0.05";
$DATE = '2018-03-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						get_CtrlVCF_from_CaseVCF
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} TtoN_vcf <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -tvcf   [s]  VCF file of case sample. <required>
	    -o      [s]  VCF file of control sample to output. <required>
	    -cnv    [s]  CNV file of pathchwork based on output of mode #1.
	                  note: a) it is required input for mode #2.
	                        b) regions of Xm0 will apply mode #1 specially.

	   # Database and Tools #
	    -dbSNP  [s]  dbSNP file, gzip and indexed. <required>
	    -btl    [s]  bcftools, minimum version: 1.3. <required>

	   # Options #
	    -mode   [i]  run mode. [1]
	                 #1(oneRef): always have at least one ref-allele.
	                 #2(origCs): follow the original alleles' read ratio in 'tvcf'.
	    -normdp [i]  depth normalize to. [30]
	    -region [s]  deal with this region, chr:from-to OR chr. [none]
	    -mtype  [s]  allowed mutation types: 'snp', 'indel', 'all'. [snp]
	    -minMTQ [i]  minimum mutation quality, i.e., the QUAL info. [20]
	    -minDP  [i]  minimum depth of position, i.e., the DP tag. [10]
	    -minMPQ [i]  minimum mapping quality, i.e., the MQ tag. [30]
	    -kpxmx  [s]  more Xmx except Xm0 to keep from '-cnv' file. [none]
	                  format: -kpxmx 4m1,5m1,5m2

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
			[ case_vcf => undef ],
			[ ctrl_vcf => undef ],
			[ pw_csv => undef ],
			# database
			[ dbSNP_idx_bgz => undef ],
			# software
			[ bcftools => '' ],
			# options
			[ depth_mode => 1 ],
			[ region_str => undef ],
			[ region => undef ],
			[ depthToNorm => 30 ],
			[ mut_type => 'snp' ],
			[ min_mut_quality => 20 ],
			[ min_pos_depth => 10 ],
			[ min_map_quality => 30 ],
			[ moreXmxToKeep => '' ],
			# setting
			[ regBase => 1E5 ],
			# intermediate variants
			[ Xm0RegInM1 => {} ],
			[ Xm0RegInM1_ForOutput => [] ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['case_vcf'],
								  ['ctrl_vcf'],
								  ['pw_csv'],
								  ['dbSNP_idx_bgz'],
								  ['bcftools']       ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-tvcf:s"	=> \$V_Href->{case_vcf},
		"-o:s"		=> \$V_Href->{ctrl_vcf},
		"-cnv:s"	=> \$V_Href->{pw_csv},
		# database
		"-dbSNP:s"	=> \$V_Href->{dbSNP_idx_bgz},
		# tool
		"-btl:s"	=> \$V_Href->{bcftools},
		# options
		"-mode:i"	=> \$V_Href->{depth_mode},
		"-normdp:i"	=> \$V_Href->{depthToNorm},
		"-region:s"	=> \$V_Href->{region_str},
		"-mtype:s"	=> \$V_Href->{mut_type},
		"-minMTQ:i"	=> \$V_Href->{min_mut_quality},
		"-minDP:i"	=> \$V_Href->{min_pos_depth},
		"-minMPQ:i"	=> \$V_Href->{min_map_quality},
		"-kpxmx:s"	=> \$V_Href->{moreXmxToKeep},
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
			 || ( !defined $V_Href->{case_vcf} || !-e $V_Href->{case_vcf} )
			 || ( !defined $V_Href->{dbSNP_idx_bgz} || !-e $V_Href->{dbSNP_idx_bgz} )
			 || ( !defined $V_Href->{bcftools} || !-e $V_Href->{bcftools} )
			 || ( $V_Href->{mut_type} !~ /^(snp|indel|all)$/ )
			 || ( $V_Href->{depth_mode} !~ /^[12]$/)
			 || (    $V_Href->{depth_mode} == 2
			 	  && ( !defined $V_Href->{pw_csv} || !-e $V_Href->{pw_csv} )
			 	)
			);
}

#--- get VCF of control sample from case sample's ---
sub get_CtrlVCF_from_CaseVCF{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	&prepare;

	# quality filter
	my $filterCaseVCF = "$V_Href->{ctrl_vcf}.filterCase.vcf.gz";
	my $exclude_parm = "QUAL<$V_Href->{min_mut_quality} || DP<$V_Href->{min_pos_depth} || MQ<$V_Href->{min_map_quality}";
	$exclude_parm .= " || TYPE!=\"$V_Href->{mut_type}\"" if( $V_Href->{mut_type} ne 'all' );
	my $region_parm = ( defined $V_Href->{region_str} ) ? "-r $V_Href->{region_str}" : '';
	my $filter_cmd = "( $V_Href->{bcftools} filter -e '$exclude_parm' $region_parm -O z -o $filterCaseVCF $V_Href->{case_vcf} ) && ( $V_Href->{bcftools} index $filterCaseVCF )";
	trible_run_for_success( $filter_cmd, 'filterCaseVCF' );

	# pick dbSNP sites
	my $dbSNP_CaseVCF = "$V_Href->{ctrl_vcf}.filterCase.dbSNP.vcf.gz";
	my $dbSNP_cmd = "$V_Href->{bcftools} filter -R $V_Href->{dbSNP_idx_bgz} $filterCaseVCF | uniq | $V_Href->{bcftools} view -O z -o $dbSNP_CaseVCF";
	trible_run_for_success( $dbSNP_cmd, 'Pick_dbSNP_sites' );

	# load Xm0 region from patch-work
	&load_Xm0_from_m1pwCNV_for_m2;

	# normalize to required depth
	# reset as 2m0 or 2m1
	&normalize_depth(
					  INVCF => $dbSNP_CaseVCF,
					  OUTVCF => $V_Href->{ctrl_vcf},
					  depthToNorm => $V_Href->{depthToNorm},
					  depth_mode => $V_Href->{depth_mode}
					);

	# sweep
	`rm -f $filterCaseVCF* $dbSNP_CaseVCF`;
}

#--- prepare work ---
sub prepare{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# moreXmxToKeep
	my %temp = map {($_,1)} split /,+/, $V_Href->{moreXmxToKeep};
	$V_Href->{moreXmxToKeep} = \%temp;

	# deal region
	if( defined $V_Href->{region_str} ){
		if( $V_Href->{region_str} =~ /^(.+):(\d+)-(\d+)$/ ){
			$V_Href->{region} = { chr => $1, st_pos=>$2, ed_pos=>$3 };
			# inform
			print "[INFO]\tDeal region is on chromosome $1, position interval is [ $2, $3 ]\n";
		}
		elsif( $V_Href->{region_str} =~ /^(.+):(\d+)\-?$/ ){
			warn_and_exit "<ERROR>\tDeal region $V_Href->{region_str} should follow format: chr:from-to OR chr\n";
		}
		else{
			$V_Href->{region} = { chr => $V_Href->{region} };
			# inform
			print "[INFO]\tDeal region is whole chromosome $V_Href->{region}->{chr}\n";
		}
	}
}

#--- load Xm0 region from patch-work based on mode #1
sub load_Xm0_from_m1pwCNV_for_m2{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return if( !defined $V_Href->{pw_csv} );

	# write selected region
	my $cnvtsv = "$V_Href->{ctrl_vcf}.TtoN_vcf.M$V_Href->{depth_mode}.select_CNVreg.tsv";
	open (REGTSV, Try_GZ_Write($cnvtsv)) || die "cannot write select_CNVreg tsv: $!\n";
	print REGTSV join("\t", qw/ Chr From To Xmx /)."\n";

	# read cnv.csv
	open (PWCSV, Try_GZ_Read($V_Href->{pw_csv})) || die "fail read pw_cnv: $!\n";
	while(<PWCSV>){
		next unless(/^"\d+"/);
		s/"//g;
		my ($Xmx) = (/,(cn\d+m\d+),/);
		next if( $Xmx !~ /m0$/ && !exists $V_Href->{moreXmxToKeep}->{$Xmx} );
		my ($chr, $stp, $edp) = (split /,/)[1,2,3];
		# filter by deal region
		if(    defined $V_Href->{region}
			&& (   $chr ne $V_Href->{region}->{chr}
				|| (   exists $V_Href->{region}->{st_pos}
					&& exists $V_Href->{region}->{ed_pos}
					&& Get_Two_Seg_Olen( $stp, $edp, $V_Href->{region}->{st_pos}, $V_Href->{region}->{ed_pos} ) == 0
				   )
			   )
		){
			next;
		}
		# record for query
		$V_Href->{Xm0RegInM1}->{$chr}->{$_}->{$stp} = $edp
			for ( int( $stp / $V_Href->{regBase} ) .. int( $edp / $V_Href->{regBase} ) );
		# output selected region
		print REGTSV join("\t", ($chr, $stp, $edp, $Xmx))."\n";
	}
	close PWCSV;
	close REGTSV;

	# inform
	print "[INFO]\tload Xm0 and keep_Xmx(-kpxmx) region from patch-work based on mode #1 for mode #2 OK.\n";

	# inform
	print "[INFO]\tcheck selected CNV region:\n".
				"\t$cnvtsv\n";
}

#--- check whether given pos is in Xm0 region from patch-work based on mode #1
sub is_from_Xm0RegInM1{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $chr = $parm{chr};
	my $pos = $parm{pos};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# check chr
	if( !exists $V_Href->{Xm0RegInM1}->{$chr} ){
		return 0;
	}
	# check region Base index
	my $regBaseIdx = int( $pos / $V_Href->{regBase} );
	if( !exists $V_Href->{Xm0RegInM1}->{$chr}->{$regBaseIdx} ){
		return 0;
	}
	# check locate-in-region
	my $regBase_Href = $V_Href->{Xm0RegInM1}->{$chr}->{$regBaseIdx};
	my $reg_stp = first { $_<=$pos && $regBase_Href->{$_}>=$pos } sort {$b<=>$a} keys %$regBase_Href;
	return ( defined $reg_stp ? 1 : 0 );
}

#--- normalize all positions' depth to required depth ---
# reset all to 2m0 or 2m1
sub normalize_depth{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $INVCF = $parm{INVCF};
	my $OUTVCF = $parm{OUTVCF};
	my $depthToNorm = $parm{depthToNorm};
	my $depth_mode = $parm{depth_mode};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	my %chrHetHomCount;
	my @chr;
	my $MARK_INFO_bool = 0;
	open (OUTVCF, "|$V_Href->{bcftools} view -O z -o $OUTVCF") || die "cannot write out_vcf: $!\n";
	open (INVCF, Try_GZ_Read($INVCF)) || die "cannot read in_vcf: $!\n";
	while(<INVCF>){
		if(/^#/){
			if( /^##INFO/ && !$MARK_INFO_bool ){
				print OUTVCF '##INFO=<ID=MARK,Number=1,Type=String,Description="the source of this control VCF">'."\n";
				$MARK_INFO_bool = 1;
			}
			print OUTVCF;
			next;
		}
		my ($chr, $pos) = (split)[0,1];
		my ($ref_f,$ref_r,$alt_f,$alt_r) = (/DP4=(\d+),(\d+),(\d+),(\d+)/);
		my $ref_sum = $ref_f + $ref_r;
		my $alt_sum = $alt_f + $alt_r;
		my $AtXm0RegInM1 = ( $depth_mode == 2 && &is_from_Xm0RegInM1( chr=>$chr, pos=>$pos ) );
		my $ref_bool = (   $depth_mode == 1 # mode #1 (oneRef)
						|| (    $AtXm0RegInM1 # mode #2 (origCs), Xm0 region use mode #1 (oneRef)
							 || $ref_sum >= 5 # mode #2 (origCs)
						   )
					   ) ? 1 : 0;
		my $alt_bool = ( $alt_sum >= 5 ) ? 1 : 0;
		if($ref_bool || $alt_bool){
			# DP4=
			my $allele_sum = ($ref_bool && $alt_bool) ? int( $depthToNorm / 2 ) : $depthToNorm;
			my $half_sum = int( $allele_sum / 2 );
			## ref
			my $ref_f_mayCt = int( $allele_sum * $ref_f / ( $ref_sum || 1 ) );
			($ref_f,$ref_r) = (!$ref_bool)
							  ? ( 0, 0 ) 
							  : ( $ref_sum >= 5
							  	  ? ( $ref_f_mayCt, $allele_sum - $ref_f_mayCt )
							  	  : (    $half_sum, $allele_sum - $half_sum )
							  	);
			## alt
			my $alt_f_mayCt = int( $allele_sum * $alt_f / ( $alt_sum || 1 ) );
			($alt_f,$alt_r) = (!$alt_bool)
							  ? ( 0, 0 )
							  : ( $alt_f_mayCt, $allele_sum - $alt_f_mayCt );
			## replace
			s/DP4=\d+,\d+,\d+,\d+/DP4=$ref_f,$ref_r,$alt_f,$alt_r/;
			# DP=
			my $DP_sum = sum($ref_f,$ref_r,$alt_f,$alt_r);
			## replace, and add mark noting the source of this control VCF
			my $mark = "NfT_M$depth_mode";
			$mark .= '_Xm0_M1' if( $AtXm0RegInM1 );
			s/DP=\d+;/DP=$DP_sum;MARK=$mark;/;
			# output
			print OUTVCF;
			# record chr orders
			if( !exists $chrHetHomCount{$chr} ){
				push @chr, $chr;
			}
			# het/hom count
			if( $ref_bool && $alt_bool ){
				$chrHetHomCount{$chr}->{'2m1'} ++;
			}
			elsif( $ref_bool ){
				$chrHetHomCount{$chr}->{'2m0_ref'} ++;
			}
			elsif( $alt_bool ){
				$chrHetHomCount{$chr}->{'2m0_alt'} ++;
			}
		}
	}
	close INVCF;
	close OUTVCF;

	# inform
	print "[INFO]\tnormalize (mode $depth_mode) allele-depth to ${depthToNorm}X level as 2m0 or 2m1.\n";

	# report stat
	my $report = "$OUTVCF.TtoN_vcf.M$depth_mode.report";
	open (REPORT, Try_GZ_Write($report)) || die "cannot write het/hom pos count report: $!\n";
	print REPORT join("\t", qw/ Chr Het_2m1 Hom_2m0_ref Hom_2m0_alt /)."\n";
	for my $chr (@chr){
		my $statHref = $chrHetHomCount{$chr};
		print REPORT join("\t", ( $chr, $statHref->{'2m1'}||0, $statHref->{'2m0_ref'}||0, $statHref->{'2m0_alt'}||0 ) )."\n";
	}
	close REPORT;
	# inform
	print "[INFO]\tcheck het/hom pos count:\n".
				"\t$report\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
