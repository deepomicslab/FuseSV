package FuseSV::Virus_Type::VCFwork;

use strict;
use warnings;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation qw/warn_and_exit trible_run_for_success/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  vcf_from_bam
			  mpileup_from_bam
			  load_VCF_and_refpos_info
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'VCFwork';
#----- version --------
$VERSION = "0.08";
$DATE = '2017-04-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						vcf_from_bam
						mpileup_from_bam
						load_VCF_and_refpos_info
						Is_I16_Alt_OK
					 /;

#--- samtools mpileup call mutations in vcf ---
sub vcf_from_bam{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
	my $Out_pref    = $Option_Href->{OutPref};
	my $bam_for_call   = $Option_Href->{Bam};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# mpileup, initially call mutation vcf4.2 for BQSR
	# disable BAQ (per-Base Alignment Quality)
	my $vcf = "$Out_pref.vcf";
	my $call_vcf_command = "$V_Href->{samtools} mpileup -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -A -B -d 1000000 -v -O -s -m 3 -L 1000000 -f $Ref_IndexByBwa -o $vcf $bam_for_call";
	trible_run_for_success($call_vcf_command, 'CallVCF', {esdo_Nvb=>1});

	return $vcf;
}

#--- samtools mpileup call mutations ---
sub mpileup_from_bam{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $Ref_IndexByBwa = $Option_Href->{RefForBwa};
	my $Out_pref       = $Option_Href->{OutPref};
	my $bam_for_call   = $Option_Href->{Bam};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# mpileup
	# disable BAQ (per-Base Alignment Quality)
	my $mpileup = "$Out_pref.mpileup";
	my $call_vcf_command = "$V_Href->{samtools} mpileup -A -B -d 1000000 -O -s -m 3 -L 1000000 -f $Ref_IndexByBwa -o $mpileup $bam_for_call";
	trible_run_for_success($call_vcf_command, 'Mpileup', {esdo_Nvb=>1});

	return $mpileup;
}

#--- determine whether alt is ok based on I16 information ---
#--- this provides merged situation of the non-ref Q13 bases ---
sub Is_I16_Alt_OK{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $I16_Aref = $Option_Href->{I16_Aref};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $I16_min_mean_baseQ = $V_Href->{I16_min_mean_baseQ};
	my $I16_min_mean_mapQ = $V_Href->{I16_min_mean_mapQ};
	my $I16_min_mean_tailDist = $Option_Href->{min_mean_tailDist} || $V_Href->{I16_min_mean_tailDist};

	# these are all based on Q13
	my $alt_reads_count = sum($I16_Aref->[2], $I16_Aref->[3]);
	my $nonRef_baseQsum = $I16_Aref->[6];
	my $nonRef_mapQsum  = $I16_Aref->[10];
	my $nonRef_tailDsum = $I16_Aref->[14];

	# only judge the read tail distance
	if(    $Option_Href->{only_judge_ReadTail}
		&& $nonRef_tailDsum / $alt_reads_count < $I16_min_mean_tailDist
	){
		return 0;
	}

	# baseQ && mapQ && tail-distance
	if(    $nonRef_baseQsum / $alt_reads_count < $I16_min_mean_baseQ
		|| $nonRef_mapQsum  / $alt_reads_count < $I16_min_mean_mapQ
		|| $nonRef_tailDsum / $alt_reads_count < $I16_min_mean_tailDist
	){
		return 0;
	}
	else{
		return 1;
	}
}

#--- read VCF and load each base info ---
sub load_VCF_and_refpos_info{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $vcf = $Option_Href->{VCF};
	my $RefPos_Href = $Option_Href->{RefPos_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	my @theme_tag;
	open (VCF, "$V_Href->{bcftools} view $vcf |") || die "Cannot bcftools read $vcf: $!\n";
	while(<VCF>){
		if(/^\#\#/){ # description
			next;
		}
		elsif(/^\#/){ # theme
			(my $theme_line = lc($_)) =~ s/^#//;
			@theme_tag = split /\s+/, $theme_line;
		}
		else{
			my @info = split;
			my %rgOb = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);

			# information
			my $refpos = $rgOb{pos};
			(my $alt = $rgOb{alt}) =~ s/,?\<\*\>$//; # v4.2; has the annoying <*>
			my $info = $rgOb{info};
			my $is_indel = ($info =~ /INDEL;/) ? 1 : 0;

			# whole depth
			# my ($DP) = ($info =~ /DP=(\d+);/);
			# $RefPos_Href->{$refpos}->{depth} ||= $DP; # one refpos may be in two lines, DP is duplicated.
			# allels
			my (@Allel_D, @Allel_D_Fw, @Allel_D_Rv, @I16_info);
			if($info =~ /AD=([^;]+);/) { @Allel_D    = split /,/, $1; }
			if($info =~ /ADF=([^;]+);/){ @Allel_D_Fw = split /,/, $1; }
			if($info =~ /ADR=([^;]+);/){ @Allel_D_Rv = split /,/, $1; }
			if($info =~ /I16=([^;]+);/){ @I16_info   = split /,/, $1; }
			## whole depth, donot use DP tag, its raw reads depth
			if(!$is_indel){ # only consider the non-indel line
				$RefPos_Href->{$refpos}->{depth} += sum(@Allel_D);
			}
			## ref-allel
			$RefPos_Href->{$refpos}->{refallel_depth} ||= $Allel_D[0];
			my $ref_allel_depth = shift @Allel_D;
			$RefPos_Href->{$refpos}->{refallel_D_Fw} = shift @Allel_D_Fw;
			$RefPos_Href->{$refpos}->{refallel_D_Rv} = shift @Allel_D_Rv;
			## mutations
			my @alt = split /,/, $alt;
			if(    scalar(@alt) != 0
				&& (  # consider IDV especially for single INDEL
					  # still filter the distance from nearest read-end
					  (    $is_indel
						&& $V_Href->{consider_IDV_bool}
						&& scalar(@alt) == 1
						&& &Is_I16_Alt_OK( { I16_Aref=>\@I16_info, only_judge_ReadTail=>1 } )
					  )
					  # full filtration
					|| &Is_I16_Alt_OK( { I16_Aref=>\@I16_info } )
				   )
			   ){
				if(!$is_indel){ # SNP, i donot want to use unless-else; once it isn't noted as INDEL, it must be one letter.
					for my $i (0 .. $#alt){
						my $alt_allel = $alt[$i];
						# record
						my $mutation_id = 'snp,' . uc($alt_allel);
						# assign the SNV to current position
						my $aimPos_Href = $RefPos_Href->{$refpos};
						$aimPos_Href->{mutation}->{$mutation_id} = [ $Allel_D[$i], $Allel_D_Fw[$i], $Allel_D_Rv[$i], $ref_allel_depth ];
					}
				}
				else{ # indels
					# IDV
					my $IDV_value = 0;
					my $use_IDV_bool = ( $V_Href->{consider_IDV_bool} && scalar(@alt) == 1 );
					if( $use_IDV_bool ){ ($IDV_value) = (/IDV=(\d+);/) }
					# this is the string around this pos
					my $ref_allel = $rgOb{ref};
					for my $i (0 .. $#alt){
						my $alt_allel = $alt[$i];
						if(length($ref_allel) < length($alt_allel)){ # insertion
							my @ref_bases = split //, $ref_allel;
							my @alt_bases = split //, $alt_allel;
							# check postfixed identical base
							my $postfix_ibase = 0;
							while( scalar(@ref_bases) != 0 && uc($ref_bases[-1]) eq uc($alt_bases[-1]) ){
								pop @ref_bases;
								pop @alt_bases;
								$postfix_ibase ++;
							}
							# check prefixed identical base
							my $prefix_ibase = 0;
							while( scalar(@ref_bases) != 0 && uc($ref_bases[0]) eq uc($alt_bases[0]) ){
								shift @ref_bases;
								shift @alt_bases;
								$prefix_ibase ++;
							}
							# get the inserted seq
							my $insert_seq = join('', @alt_bases);
							# check
							if( scalar(@ref_bases) != 0 || scalar(@alt_bases) == 0 ){
								warn_and_exit "<ERROR>:\tWrong VCF insertion info.\nref_allel:$ref_allel; alt_allel:$alt_allel\n";
							}
							# record
							my $mutation_id = 'ins,' . uc($insert_seq);
							# assign the insertion to the position before its first base
							my $aimPos_Href = $RefPos_Href->{ $refpos + $prefix_ibase - 1 };
							if( $use_IDV_bool ){ # use IDV to pretend
								$aimPos_Href->{mutation}->{$mutation_id} = [ $IDV_value, int($IDV_value/2), int($IDV_value/2), $ref_allel_depth ];
							}
							else{
								$aimPos_Href->{mutation}->{$mutation_id} = [ $Allel_D[$i], $Allel_D_Fw[$i], $Allel_D_Rv[$i], $ref_allel_depth ];
							}
						}
						elsif(length($ref_allel) > length($alt_allel)){ # deletion
							my @ref_bases = split //, $ref_allel;
							my @alt_bases = split //, $alt_allel;
							# check postfixed identical base
							my $postfix_ibase = 0;
							while( scalar(@alt_bases) != 0 && uc($ref_bases[-1]) eq uc($alt_bases[-1]) ){
								pop @ref_bases;
								pop @alt_bases;
								$postfix_ibase ++;
							}
							# check prefixed identical base
							my $prefix_ibase = 0;
							while( scalar(@alt_bases) != 0 && uc($ref_bases[0]) eq uc($alt_bases[0]) ){
								shift @ref_bases;
								shift @alt_bases;
								$prefix_ibase ++;
							}
							# get the deleted seq
							my $delete_seq = join('', @ref_bases);
							# check
							if( scalar(@alt_bases) != 0 || scalar(@ref_bases) == 0 ){
								warn_and_exit "<ERROR>:\tWrong VCF deletion info.\nref_allel:$ref_allel; alt_allel:$alt_allel\n";
							}
							# record
							my $mutation_id = 'del,' . uc($delete_seq);
							# assign the deletion to the first base of itself
							my $aimPos_Href = $RefPos_Href->{ $refpos + $prefix_ibase };
							$aimPos_Href->{depth} += $Allel_D[$i]; # count the deleted part into depth also
							# last element is for the aimPos, not this one
							if( $use_IDV_bool ){ # use IDV to pretend
								$aimPos_Href->{mutation}->{$mutation_id} = [ $IDV_value, int($IDV_value/2), int($IDV_value/2), 0 ];
							}
							else{
								$aimPos_Href->{mutation}->{$mutation_id} = [ $Allel_D[$i], $Allel_D_Fw[$i], $Allel_D_Rv[$i], 0 ];
							}
						}
					}
				}
			}
		}
	}
	close VCF;
}

#--- 
1; ## tell the perl script the successful access of this module.
