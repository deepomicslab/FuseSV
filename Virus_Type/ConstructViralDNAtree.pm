package FuseSV::Virus_Type::ConstructViralDNAtree;

use strict;
use warnings;
use File::Spec qw/catfile/;
use File::Basename qw/basename dirname/;
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ trible_run_for_success /;
use FuseSV::LoadOn;
use FuseSV::Virus_Type::VirusTypeMain;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  construct_viral_DNA_tree
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ConstructViralDNAtree';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-01-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Get_Cmd_Options
						construct_viral_DNA_tree
						merge_fa_for_tree
						multi_align_by_MAFFT
						conserved_blocks_by_Gblocks
						check_available_fasta
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} virus_tree <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -vfa    [s]  viral genome(s) file, fasta format. <required>
	                  This file contains all viral variants' genomes in your tree.
	    -vpf    [s]  virus database folder. <required>
	                  You can make it via function 'make_virusdb'.
	                  OR Download these widely used from FuseSV official website, e.g., HPV.
	    -odir   [s]  directory to store outputs, must exist already. <required>

	   # Software Required #
	    -MAFFT    [s]  tool MAFFT, minimum version: 7.221. <required>
	    -Gblocks  [s]  tool Gblocks, minimum version: 0.91b. <required>
	    -raxmlHPC [s]  tool raxmlHPC, minimum version: 8.2.0. <required>

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
		(
			# input/output
			[ tree_outdir => '' ],
			[ virus_fa_for_tree => '' ],
			# software
			[ MAFFT => '' ],
			[ Gblocks => '' ],
			[ raxmlHPC => '' ],
			# basic
			[ virus_tree_prefix => 'viral_genome.for.tree' ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['virus_fa_for_tree'],
								  ['virus_db_folder'],
								  ['tree_outdir'],
								  ['MAFFT'],
								  ['Gblocks'],
								  ['raxmlHPC']			    ] ]
		);

	# load containers
	FuseSV::Virus_Type::VirusTypeMain->Load_moduleVar_to_pubVarPool( only_load_containers => 1 );
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-vfa:s"	=> \$V_Href->{virus_fa_for_tree},
		"-vpf:s"	=> \$V_Href->{virus_db_folder},
		"-odir:s"	=> \$V_Href->{tree_outdir},
		# software
		"-MAFFT:s"	=> \$V_Href->{MAFFT},
		"-Gblocks:s"=> \$V_Href->{Gblocks},
		"-raxmlHPC:s"=> \$V_Href->{raxmlHPC},
		# help
		"-h|help"	=> \$V_Href->{HELP},
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder})
			 || (!$V_Href->{tree_outdir} || !-d $V_Href->{tree_outdir} || !-w $V_Href->{tree_outdir})
			 || (!$V_Href->{virus_fa_for_tree} || !-e $V_Href->{virus_fa_for_tree})
			 || (!$V_Href->{MAFFT} || !-e $V_Href->{MAFFT})
			 || (!$V_Href->{Gblocks} || !-e $V_Href->{Gblocks})
			 || (!$V_Href->{raxmlHPC} || !-e $V_Href->{raxmlHPC})
			);
}

#--- based on the virus info, prepare virus ref for first alignment ---
sub construct_viral_DNA_tree{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Virus_Info_Href = $V_Href->{Virus_Info_Href};

	# prepare all viral genome fasta sequences in one single file
	my $merged_viral_fa = &merge_fa_for_tree;

	# MAFFT multi-Alignment
	my $multiAligned = &multi_align_by_MAFFT( { merged_fa => $merged_viral_fa } );

	# Gblocks Selection of conserved blocks
	my $blocks_file = &conserved_blocks_by_Gblocks( { multiAligned => $multiAligned } );

	# raxmlHPC to get the newick tree file
	&get_tree_newick_by_raxmlHPC( { blocks_file => $blocks_file, merged_fa => $merged_viral_fa } );
}

#--- merge virus variant fasta file with vrial database package ---
sub merge_fa_for_tree{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# check virus variant fasta file
	&check_available_fasta( { fasta_file => $V_Href->{virus_fa_for_tree} } );

	# merge
	my $merged_viral_fa = File::Spec->catfile($V_Href->{tree_outdir},"$V_Href->{virus_tree_prefix}.fa");
	my $mergeFA_command = "cat $V_Href->{virus_fa_for_tree} $V_Href->{virus_orig_genome_fa} > $merged_viral_fa";
	trible_run_for_success($mergeFA_command, 'mergeFA', {cmd_Nvb=>1,esdo_Nvb=>1});

	return $merged_viral_fa;
}

#--- MAFFT multi-Alignment ---
# http://mafft.cbrc.jp/alignment/software/
sub multi_align_by_MAFFT{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $merged_viral_fa = $Option_Href->{merged_fa};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	stout_and_sterr `date`."[INFO]\tMAFFT multi-Alignment starts.\n";

	my $MAFFT_command = "$V_Href->{MAFFT} --quiet --auto $merged_viral_fa > $merged_viral_fa.multiAligned";
	trible_run_for_success($MAFFT_command, 'MAFFT_Multi_Align', {cmd_Nvb=>1,esdo_Nvb=>1});

	stout_and_sterr `date`."[INFO]\tMAFFT multi-Alignment finishes.\n";

	return "$merged_viral_fa.multiAligned";
}

#--- Gblocks Selection of conserved blocks ---
# http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html
sub conserved_blocks_by_Gblocks{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $multiAligned = $Option_Href->{multiAligned};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	stout_and_sterr `date`."[INFO]\tGblocks Selection of conserved blocks starts.\n";

	my $Gblocks_command = "$V_Href->{Gblocks} $multiAligned -t=d -v=100";
	# Gblocks always returns '0'.
	# trible_run_for_success($Gblocks_command, 'Gblocks', {cmd_Nvb=>1,esdo_Nvb=>1});
	`$Gblocks_command`;

	stout_and_sterr `date`."[INFO]\tGblocks Selection of conserved blocks finishes.\n";

	return "$multiAligned-gb";
}

#--- raxmlHPC to get the newick tree file ---
# http://evolution.genetics.washington.edu/phylip/newicktree.html
# http://sco.h-its.org/exelixis/web/software/raxml/hands_on.html
sub get_tree_newick_by_raxmlHPC{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $blocks_file = $Option_Href->{blocks_file};
	my $merged_viral_fa = $Option_Href->{merged_fa};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	stout_and_sterr `date`."[INFO]\traxmlHPC for the newick tree file starts.\n";

	my $raxmlHPC_command = "( cd $V_Href->{tree_outdir}; $V_Href->{raxmlHPC} -f a -m GTRCAT -x 12345 -p 12345 -N 20 -d -s $blocks_file -n $V_Href->{virus_tree_prefix}.Tree )";
	trible_run_for_success($raxmlHPC_command, 'raxmlHPC_tree', {cmd_Nvb=>1,esdo_Nvb=>1});

	# rename
	my @RAxML_results = glob( File::Spec->catfile($V_Href->{tree_outdir},'RAxML_*.Tree') );
	for my $filepath (@RAxML_results){
		my $filename = basename($filepath);
		(my $new_postfix = $filename) =~ s/\.$V_Href->{virus_tree_prefix}.Tree//;
		my $new_filepath = File::Spec->catfile($V_Href->{tree_outdir}, "$V_Href->{virus_tree_prefix}.Tree.$new_postfix");
		`mv $filepath $new_filepath`;
	}

	stout_and_sterr `date`."[INFO]\traxmlHPC for the newick tree file finishes.\n"
						  ."\tYou can import software 'Dendroscope' with newick format results, such as,\n"
						  ."\t $V_Href->{virus_tree_prefix}.Tree.RAxML_bestTree\n";
}

#--- check fasta file availability ---
sub check_available_fasta{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $fasta_file = $Option_Href->{fasta_file};

	# line number must >= 2
	my $line_number = `cat $fasta_file | head -10 | wc -l`;
	# check '>' must >= 1
	my $FAheader_line_number = `cat $fasta_file | head -10 | grep '^>' | wc -l`;
	if(    $line_number < 2
		|| $FAheader_line_number < 1
	  ){
		warn_and_exit `date`."<ERROR>\tincomplete fasta_file.\n"
							."\t$fasta_file.\n";
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
