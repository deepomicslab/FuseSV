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

my ($gene_group_list);
my ($vinteg_list);
my ($output);
my ($gene_psl);
my ($Help);

GetOptions(
	"-ggl:s"	=> \$gene_group_list,
	"-vil:s"	=> \$vinteg_list,
	"-gpsl:s"	=> \$gene_psl,
	"-output:s"	=> \$output,
	"-help|h"	=> \$Help
);

#-------------- show help ----------------
if(    $Help
	|| !$gene_group_list || !-e $gene_group_list
	|| !$vinteg_list || !-e $vinteg_list
	|| !$gene_psl || !-e $gene_psl
){
	warn "perl $RealScript <-ggl> <-vil> <-gpsl> <-output>\n";
	exit(1);
}

# load gene group info
my %GeneGroup;
my %Gene;
open (GGL, $gene_group_list) || die "fail read gene_group_list: $!\n";
while(<GGL>){
	next if(/^#/);
	my ($gene_group, $gene_name, $chr, $cytoband) = (split)[0,1,2,3];
	$GeneGroup{$gene_group}->{gene_list}->{$gene_name} = 1;
	$Gene{$gene_name}->{chr} = $chr;
	$Gene{$gene_name}->{cyb}->{$cytoband} = 1;
}
close GGL;

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

# load vinteg info
my %vItg;
open (VIL, $vinteg_list) || die "fail read vinteg_list: $!\n";
while(<VIL>){
	next if(/^#/);
	chomp;
	my ($itg_chr, $itg_cyb, $itg_pos) = (split)[1,2,3];
	push @{$vItg{$itg_chr}{$itg_cyb}{$itg_pos}}, $_;
}
close VIL;

# compare
open (OUT, Try_GZ_Write($output)) || die "fail write output: $!\n";
for my $gene_group (sort keys %GeneGroup){
	for my $gene_name (sort keys %{$GeneGroup{$gene_group}->{gene_list}}){
		my $chr = $Gene{$gene_name}->{chr};
		next unless(exists($vItg{$chr}));
		for my $cyb (sort keys %{$Gene{$gene_name}->{cyb}}){
			next unless(exists($vItg{$chr}{$cyb}));
			my $gene_stGp = $Gene{$gene_name}->{stGp};
			my $gene_edGp = $Gene{$gene_name}->{edGp};
			if(!defined($gene_stGp) || !defined($gene_edGp)){
				die "$gene_name\n";
			}
			for my $itg_pos (sort {$a<=>$b} keys %{$vItg{$chr}{$cyb}}){
				my $dist;
				# die "$gene_name\t$gene_stGp\t$gene_edGp\t$itg_pos\n";
				if( $itg_pos < $gene_stGp ){
					$dist = $gene_stGp - $itg_pos;
				}
				elsif( $itg_pos > $gene_edGp ){
					$dist = $itg_pos - $gene_edGp;
				}
				else{
					$dist = 0;
				}
				print OUT $_."\t$gene_group\t$gene_name\t$dist\n" for @{$vItg{$chr}{$cyb}{$itg_pos}};
			}
		}
	}
}
close OUT;
