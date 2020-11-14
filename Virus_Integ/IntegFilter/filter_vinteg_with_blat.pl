#!/usr/bin/perl -w

use strict;
use FindBin qw/$RealBin $RealScript/;
use File::Basename qw/basename dirname/;
use File::Spec qw/catfile abs2rel/;
use Getopt::Long;
use List::Util qw/min max sum/;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;

my ($seq_for_PCR_file);
my ($blat_bin);
my ($host_seq_m8);
my ($out_dir);
my ($host_genome);
my ($Help);

GetOptions(
	"-sqfp:s"	=> \$seq_for_PCR_file,
	"-blat:s"	=> \$blat_bin,
	"-hostg:s"	=> \$host_genome,
	"-hostm8:s"	=> \$host_seq_m8,
	"-odir:s"	=> \$out_dir,
	"-help|h"	=> \$Help
);

#-------------- show help ----------------
if(    $Help
	|| !$seq_for_PCR_file || !-e $seq_for_PCR_file
	|| (!$blat_bin && !$host_seq_m8)
	|| !$host_genome || !-e $host_genome
	|| !$out_dir || !-d $out_dir
){
	warn "perl $RealScript <-sqfp> <-hostg> [-hostm8/-blat] <-odir>\n";
	exit(1);
}

my $host_seq_fa = File::Spec->catfile( $out_dir, 'host_seq_for_blat.fa' );

my %vInteg;
open (HSEQ,">$host_seq_fa") || die "fail write host_seq_fa: $!\n";
open (SEQP,$seq_for_PCR_file) || die "fail read seq_for_PCR_file: $!\n";
while(<SEQP>){
	my ($case, $itg_seq) = split;
	my @ele = split /\//,$case;
	my $refseg_5 = $ele[1];
	my @itg_seq = split /\|/, $itg_seq;
	my $host_seq = ($refseg_5 =~ /chr/)?$itg_seq[0]:$itg_seq[-1];
	$host_seq =~ s/^[a-z]+//;
	$host_seq =~ s/[a-z]+$//;
	my $case_id = join('/',@ele[1,2,4,5,6]);
	print HSEQ ">$case_id\n";
	print HSEQ "$host_seq\n";
	$vInteg{$case_id}->{host_chr} = ($refseg_5 =~ /chr/)?$ele[1]:$ele[4];
	$vInteg{$case_id}->{host_pos} = ($refseg_5 =~ /chr/)?$ele[2]:$ele[5];
	$vInteg{$case_id}->{host_loc} = ($refseg_5 =~ /chr/)?'Up':'Dw';
	$vInteg{$case_id}->{host_seq} = $host_seq;
	$vInteg{$case_id}->{hseq_len} = length($host_seq);
}
close SEQP;
close HSEQ;

stout_and_sterr "Host Seq for Blat ok.\n";


if( !$host_seq_m8 || !-e $host_seq_m8 ){
	$host_seq_m8 = File::Spec->catfile( $out_dir, 'host_seq.m8' );
	my $blat_cmd = "$blat_bin -t=dna -q=dna -out=blast8 $host_genome $host_seq_fa $host_seq_m8";
	trible_run_for_success($blat_cmd, 'blat_hseq', {cmd_Nvb=>1, esdo_Nvb=>1});
}


open (M8,$host_seq_m8) || die "fail read host_seq m8: $!\n";
while(<M8>){
	my ($case_id, $mseg, $mlen, $host_st, $host_ed, $score) = (split)[0,1,3,8,9,11];
	push @{$vInteg{$case_id}->{blat}}, {mseg=>$mseg, mlen=>$mlen, host_st=>$host_st, host_ed=>$host_ed, score=>$score};
}
close M8;

stout_and_sterr "Read Host seq blat m8 ok.\n";

my $marked_vInteg = File::Spec->catfile( $out_dir, 'vInteg_marked.list' );

open (MITG,">$marked_vInteg") || die "fail write marked_vInteg: $!\n";
for my $case_id (sort keys %vInteg){
	unless( exists $vInteg{$case_id}->{blat} ){
		print MITG "$case_id\t$vInteg{$case_id}->{host_loc}-TooShortToBlat\n";
		next;
	}
	my @mHref_sortByScore = sort {$b->{score} <=> $a->{score}} @{$vInteg{$case_id}->{blat}};
	my $MoreBlatMapBool = (scalar(@mHref_sortByScore) > 1);
	my $score_diff_abs = $MoreBlatMapBool ? abs($mHref_sortByScore[0]->{score} - $mHref_sortByScore[1]->{score}) : 0;
	my $mlen_diff_abs = $MoreBlatMapBool ? abs($mHref_sortByScore[0]->{mlen} - $mHref_sortByScore[1]->{mlen}) : 0;
	if(    $mHref_sortByScore[0]->{mseg} ne $vInteg{$case_id}->{host_chr}
		|| (   abs($mHref_sortByScore[0]->{host_st} - $vInteg{$case_id}->{host_pos}) > 3
			&& abs($mHref_sortByScore[0]->{host_ed} - $vInteg{$case_id}->{host_pos}) > 3
		   )
		|| (   $MoreBlatMapBool
			&& $mlen_diff_abs <= 3
			&& $score_diff_abs <= 5
		   )
	){
		print MITG "$case_id\t$vInteg{$case_id}->{host_loc}-MultipleMappedEqually\n";
	}
	elsif(   $MoreBlatMapBool
		  && $score_diff_abs <= 10
	){
		print MITG "$case_id\t$vInteg{$case_id}->{host_loc}-MultipleMappedButBest\n";
	}
	elsif($vInteg{$case_id}->{hseq_len} < 100){
		print MITG "$case_id\t$vInteg{$case_id}->{host_loc}-Short-But-SeqMap-OK\n";
	}
	else{
		print MITG "$case_id\t$vInteg{$case_id}->{host_loc}-SeqMap-OK\n";
	}
}
close MITG;

stout_and_sterr "Filter vInteg with Blat m8 ok.\n";
