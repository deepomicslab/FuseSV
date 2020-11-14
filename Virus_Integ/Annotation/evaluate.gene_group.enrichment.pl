#!/usr/bin/perl -w

use strict;
use FindBin qw/$RealBin $RealScript/;
use File::Basename qw/basename dirname/;
use File::Spec qw/catfile abs2rel/;
use Getopt::Long;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;

# my ($gene_group_list);
my ($vinteg_list);
my ($vInteg_geneDist_list);
my ($min_vitg_dist)=(500000);
my ($gene_psl);
my ($simu_region_list);
my ($simu_reg_step)=(10000);
my ($simulation_time)=(10000);
my ($max_genedist)=(750000);
my ($max_JRstNO)=(5);
my ($Help);

GetOptions(
	# "-ggl:s"	=> \$gene_group_list,
	"-vil:s"	=> \$vinteg_list,
	"-vgdl:s"	=> \$vInteg_geneDist_list,
	"-gpsl:s"	=> \$gene_psl,
	"-reg:s"	=> \$simu_region_list,
	"-vidist:i"	=> \$min_vitg_dist,
	"-gdist:i"	=> \$max_genedist,
	"-jrstno:i"	=> \$max_JRstNO,
	"-simut:i"	=> \$simulation_time,
	"-help|h"	=> \$Help
);

#-------------- show help ----------------
if(    $Help
	# || !$gene_group_list || !-e $gene_group_list
	|| !$vinteg_list || !-e $vinteg_list
	|| !$vInteg_geneDist_list || !-e $vInteg_geneDist_list
	|| !$gene_psl || !-e $gene_psl
	|| !$simu_region_list || !-e $simu_region_list
){
	warn "perl $RealScript <-vil> <-vgdl> <-gpsl> <-reg> [-vidist/-gdist/-jrstno/-simut]\n"; #<-ggl> 
	exit(1);
}

# load vinteg info
my %Sample2vInteg;
open (VIL, $vinteg_list) || die "fail read vinteg_list: $!\n";
while(<VIL>){
	next if(/^#/);
	chomp;
	my ($sample, $itg_chr, $itg_cyb, $itg_pos, $JRstNO) = (split)[0,1,2,3,8];
	next if($JRstNO > $max_JRstNO);
	$Sample2vInteg{$sample}{$itg_chr}{$itg_pos} = 1;
}
close VIL;
warn "read viurs integion list ok.\n";

# discard close vInteg pos
if($min_vitg_dist){
	for my $sample (sort keys %Sample2vInteg){
		for my $itg_chr (sort keys %{$Sample2vInteg{$sample}}){
			FILTER:{
				my @itg_pos = sort {$a<=>$b} keys %{$Sample2vInteg{$sample}{$itg_chr}};
				for my $i (1..$#itg_pos){
					if( $itg_pos[$i] - $itg_pos[$i-1] < $min_vitg_dist){
						delete $Sample2vInteg{$sample}{$itg_chr}{$itg_pos[$i]};
						redo FILTER;
					}
				}
			}
		}
	}
	warn "filter closed vInteg.\n";
}

# get vInteg count of each sample
my %Sample_vInteg_Count;
for my $sample (sort keys %Sample2vInteg){
	$Sample_vInteg_Count{$sample} = 0;
	for my $itg_chr (sort keys %{$Sample2vInteg{$sample}}){
			$Sample_vInteg_Count{$sample} += scalar(keys %{$Sample2vInteg{$sample}{$itg_chr}});
	}
	warn "Sample2vItgCount: $sample\t$Sample_vInteg_Count{$sample}\n";
}
warn "Get vInteg count of samples\n";

# load vinteg gene dist info
my %Obv_Group_Posiive_Sample;
my %GeneGroup;
my %Gene;
open (VDL, $vInteg_geneDist_list) || die "fail read vInteg_geneDist_list: $!\n";
while(<VDL>){
	next if(/^#/);
	chomp;
	my ($sample, $chr, $cytoband, $JRstNO, $gene_group, $gene_name, $gene_dist) = (split)[0,1,2,8,9,10,11];
	next if($JRstNO > $max_JRstNO || $gene_dist > $max_genedist);
	# mark as positive
	$Obv_Group_Posiive_Sample{$gene_group}{$sample} = 1;
	# load gene group info
	$GeneGroup{$gene_group}->{gene_list}->{$gene_name} = 1;
	$Gene{$gene_name}->{chr} = $chr;
	$Gene{$gene_name}->{cyb}->{$cytoband} = 1;
}
close VDL;
warn "get observed position samples of groups.\n";

warn "$_\t".join(',', sort keys %{$GeneGroup{$_}->{gene_list}})."\n" for sort keys %GeneGroup;

warn "$_\t".scalar(keys %{$Obv_Group_Posiive_Sample{$_}})."\n" for sort keys %Obv_Group_Posiive_Sample;

# load gene group info
# my %GeneGroup;
# my %Gene;
# open (GGL, $gene_group_list) || die "fail read gene_group_list: $!\n";
# while(<GGL>){
# 	next if(/^#/);
# 	my ($gene_group, $gene_name, $chr, $cytoband) = (split)[0,1,2,3];
# 	$GeneGroup{$gene_group}->{gene_list}->{$gene_name} = 1;
# 	$Gene{$gene_name}->{chr} = $chr;
# 	$Gene{$gene_name}->{cyb}->{$cytoband} = 1;
# }
# close GGL;
# warn "load gene group ok.\n";

# load gene info
open (GPSL, Try_GZ_Read($gene_psl)) || die "fail read gene_psl: $!\n";
while(<GPSL>){
	my ($strd, $gene_name, $chr, $cytoband, $st_Gpos, $ed_Gpos) = (split)[8,9,13,14,15,16];
	next unless( exists($Gene{$gene_name}) );
	$Gene{$gene_name}->{strd} = $strd;
	$Gene{$gene_name}->{stGp} = $st_Gpos;
	$Gene{$gene_name}->{edGp} = $ed_Gpos;
}
close GPSL;
warn "load gene information ok.\n";

warn "GENEs:\t".join(',', sort keys %Gene)."\n";

# load region for simulation
my (@simu_region);
open (SIMUR,$simu_region_list) || die "fail read simu_region_list: $!\n";
while(<SIMUR>){
	my ($refseg,$st,$ed) = (split);
	next if($refseg !~ /^chr[\d+X]+$/i);
	for (my $i=$st; $i<=$ed; $i+=$simu_reg_step) {
		my $load_st = $i;
		my $load_ed = $i + $simu_reg_step - 1;
		$load_ed = $ed if($load_ed > $ed);
		push @simu_region, [$refseg, $load_st, $load_ed];
	}
}
close SIMUR;
warn "load simulated region.\n";

# simulation work
my %GeneGroup_Positive_time;
for my $simu_time (1 .. $simulation_time){
	warn "simu_time: $simu_time\n" if($simu_time % 100 == 0);
	my %GeneGroup_Simu_positive_sam;
	for my $sample (sort keys %Sample_vInteg_Count){
		my $vInteg_Count = $Sample_vInteg_Count{$sample};
		while($vInteg_Count--){
			## position
			my $simu_reg_Aref = $simu_region[int(rand(scalar(@simu_region)))];
			my ($itg_chr,$simu_reg_st,$simu_reg_ed) = @$simu_reg_Aref;
			my $itg_pos = $simu_reg_st + int(rand($simu_reg_ed-$simu_reg_st+1));
			## test position
			for my $gene_group (sort keys %GeneGroup){
				next if( exists $GeneGroup_Simu_positive_sam{$gene_group}{$sample} );
				for my $gene_name (sort keys %{$GeneGroup{$gene_group}->{gene_list}}){
					my $chr = $Gene{$gene_name}->{chr};
					next if($chr ne $itg_chr);
					my $gene_stGp = $Gene{$gene_name}->{stGp};
					my $gene_edGp = $Gene{$gene_name}->{edGp};
					my $dist;
					if( $itg_pos < $gene_stGp ){
						$dist = $gene_stGp - $itg_pos;
					}
					elsif( $itg_pos > $gene_edGp ){
						$dist = $itg_pos - $gene_edGp;
					}
					else{
						$dist = 0;
					}
					if($dist <= $max_genedist){
						$GeneGroup_Simu_positive_sam{$gene_group}{$sample} = 1;
						last;
					}
				}
			}
		}
	}
	# calculate positive case count of each gene group
	for my $gene_group (sort keys %GeneGroup){
		my $obv_sample_count = scalar(keys %{$Obv_Group_Posiive_Sample{$gene_group}});
		my $sim_sample_count = scalar(keys %{$GeneGroup_Simu_positive_sam{$gene_group}});
		if($sim_sample_count >= $obv_sample_count){
			$GeneGroup_Positive_time{$gene_group}++;
		}
	}
}

# output the positive ratio
for my $gene_group (sort keys %GeneGroup){
	my $positive_time = $GeneGroup_Positive_time{$gene_group} || 0;
	my $prob = $positive_time / $simulation_time;
	print "$gene_group\t$simulation_time\t$positive_time\t$prob\n";
}
