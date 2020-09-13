#!/bin/usr/perl -w
use strict;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use Getopt::Long;
use List::Util qw[min max sum];
use File::Spec qw/catfile/;
use FindBin qw/$RealBin $RealScript/;
use FuseSV::Virus_Type::CreateMutatedVirusGenome;

my ($VirusDatabaseSorDir);
my ($MutatedVirusOutDir);
my ($snv_freq) = (0.01);
my ($indel_freq) = (0.001);
my ($MNP_freq) = (0.001);
my ($indel_localDiff) = (1);
my ($max_InDel_Length) = (15);
my ($max_MNP_Length) = (3);
my (@select_virus) = ();
my ($help);

GetOptions(
	# input
	"-vpf=s"	=>	\$VirusDatabaseSorDir,
	"-odir=s"	=>	\$MutatedVirusOutDir,
	# mutations
	"-virus=s"	=>	\@select_virus,
	"-SNVf=f"	=>	\$snv_freq,
	"-INDELf=f"	=>	\$indel_freq,
	"-MNPf=f"	=>	\$MNP_freq,
	"-mxIDL=i"	=>	\$max_InDel_Length,
	"-mxMNPL=i"	=>	\$max_MNP_Length,
	"-idlocd=i"	=>	\$indel_localDiff,
	# help
	"-h|help"	=>	\$help
);

if(!$VirusDatabaseSorDir || !$MutatedVirusOutDir || $help){
	warn "
	 Usage:   perl $RealScript <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -vpf    [s]  virus package folder. <required>
	                  Download it from FuseSV official website, such as HBV, HPV package.
	    -odir   [s]  directory to store outputs, must exist already. <required>

	   # Mutation Option #
	    -virus  [s]  Select certain virus. [all]
	    -SNVf   [f]  SNV frequency. [0.01]
	    -INDELf [f]  INDEL frequency. [0.001]
	    -MNPf   [f]  INDEL frequency. [0.001]
	    -mxIDL  [i]  Maximum length of random InDel. [15]
	    -mxMNPL [i]  Maximum length of random MNP. [3]
	    -idlocd [i]  Require certain local different bases to avoid shift of InDel. [1]

	    -h|help   Display this help info.

	 Version: 1.01
	 Date:    2017-02-12
	 Author:  Wenlong Jia, wenlongkxm\@gmail.com
	\n";
	exit(1);
}

#-----------------#
# load virus info #
#-----------------#
my $virus_info_file = File::Spec->catfile($VirusDatabaseSorDir, 'virus.info');
my %select_virus = map{($_,1)} @select_virus;
my %virusInfo;
FuseSV::Virus_Type::CreateMutatedVirusGenome->load_virus_info(
		{
		 virusInfo_Href=>\%virusInfo,
		 virusInfo_file=>$virus_info_file,
		 needvirus_Href=>(scalar(@select_virus)==0)?undef:\%select_virus
		}
	);

#---------------------#
# load each virus seq #
#---------------------#
my $virus_orig_genome_fa = File::Spec->catfile($VirusDatabaseSorDir, 'virus_original_reference', 'virus.genome.fa');
FuseSV::Virus_Type::CreateMutatedVirusGenome->load_virus_genome(
		{
		 virusInfo_Href=>\%virusInfo,
		 virus_orig_fa=>$virus_orig_genome_fa
		}
	);

#--------------------#
# generate mutations #
#--------------------#
FuseSV::Virus_Type::CreateMutatedVirusGenome->mutate_virus(
		{
		 virusInfo_Href=>\%virusInfo,
		 snv_freq=>$snv_freq,
		 indel_freq=>$indel_freq,
		 MNP_freq=>$MNP_freq,
		 indel_localDiff=>$indel_localDiff,
		 max_InDel_Length=>$max_InDel_Length,
		 max_MNP_Length=>$max_MNP_Length
		}
	);

#----------------------#
# output mutated virus #
# and mutations on it  #
#----------------------#
FuseSV::Virus_Type::CreateMutatedVirusGenome->output_mutated_virus(
		{
		 virusInfo_Href=>\%virusInfo,
		 MutatedVirusOutDir=>$MutatedVirusOutDir
		}
	);
