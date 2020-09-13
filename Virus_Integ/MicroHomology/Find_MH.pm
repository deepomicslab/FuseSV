package FuseSV::Virus_Integ::MicroHomology::Find_MH;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use POSIX qw/ strftime /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::FileHeader qw/ getHeaderTag lineInfoToHash /;
use FuseSV::LoadOn;
use FuseSV::GetPath qw/ GetPath /;
use FuseSV::Objects::vInteg_OB;
use FuseSV::Virus_Type::DrawViralGenome qw/ show_figure_information /;
use SVG;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  find_vInteg_MH
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSV::Virus_Integ::MicroHomology::Find_MH';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-11-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						find_vInteg_MH
						load_vIteg_list
						load_seq
						search_MH
						output_MH_report
						draw_MH_figure
					 /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
	 Usage:   perl $V_Href->{MainName} find_MH <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -vitg   [s]  virus integration list. <required>
	    -odir   [s]  folder to store results. <required>
	    -samid  [s]  sample-ID to be the output prefix. <required>

	   # DataBase #
	    -href   [s]  host reference fasta file with fai file. <[optinal]>
	    -vref   [s]  virus reference for integration. <[optinal]>

	   # Tool #
	    -samt   [s]  SAMtools. <required>

	   # Options #
	    -flank  [i]  size of flank region to check micro-homology. [20]
	    -minvtm [i]  minimum size of vertial micro-homology. [2]
	    -maxvtm [i]  maximum size of vertial micro-homology. [10]
	    -minsfm [i]  minimum size of shifting micro-homology. [2]
	    -maxsfm [i]  maximum size of shifting micro-homology. [10]
	    -minsfo [i]  minimum overlap size of shifting micro-homology. [0]
	    -maxsft [i]  maximum shift distance of shifting micro-homology. [5]
	    -sfcnb       enable showing circular number of shifting MH. [disabled]

	    -h|help      Display this help info.

	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
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
			#--- input/output
			[ VirusIntegList => '' ],
			[ output_dir => '' ],
			[ sample_ID => undef ],
			#--- database
			[ host_ref => undef ],
			[ vref_forItg => undef ],
			#--- tools
			[ samtools => '' ],
			#--- options
			[ FlankSize => 20 ],
			[ min_vtMH => 2 ],
			[ max_vtMH => 10 ],
			[ min_sfMH => 2 ],
			[ max_sfMH => 10 ],
			[ min_sfMHoplen => 0 ],
			[ max_sfMHshift => 5 ],
			#--- intermediate variants
			[ vIntegOBs_Aref => [] ],
			[ refsegToVitgOB_Href => {} ],
			[ validCaseForMH => 0 ],

			#--- SVG parm
			[ SVG_raw_gap => 3 ],
			[ SVG_vItg_gap => 40 ],
			[ SVG_bgCol => { p5=>'lime', p3=>'yellow', sMH=>'pink', lMH=>'red' } ],
			[ SVG_SeqBaseRectWidth => 9 ],
			[ SVG_SeqBaseRectHeight => 14 ],
			[ SVG_JuncSiteGapRatio => 0.5 ], # to 'SVG_SeqBaseRectWidth'
			[ SVG_sfMH_cirNum => 0 ],
			[ SVG_sfMH_radius => 4 ],
			[ SVG_sfMH_dotR => 1.5 ],
			[ SVG_X_t_toInt => 10 ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['VirusIntegList'],
								  ['output_dir'],
								  ['host_ref'],
								  ['vref_forItg'],
								  ['samtools']  ] ]
		);
}

#--- get options from command line ---
sub Get_Cmd_Options{
	# get options
	GetOptions(
		# input/output
		"-vitg:s"	=> \$V_Href->{VirusIntegList},
		"-odir:s"	=> \$V_Href->{output_dir},
		"-samid:s"	=> \$V_Href->{sample_ID},
		# tools
		"-samt:s"	=> \$V_Href->{samtools},
		# database
		"-href:s"	=> \$V_Href->{host_ref},
		"-vref:s"	=> \$V_Href->{vref_forItg},
		# options
		"-flank:i"	=> \$V_Href->{FlankSize},
		"-minvtm:i"	=> \$V_Href->{min_vtMH},
		"-maxvtm:i"	=> \$V_Href->{max_vtMH},
		"-minsfm:i"	=> \$V_Href->{min_sfMH},
		"-maxsfm:i"	=> \$V_Href->{max_sfMH},
		"-minsfo:i"	=> \$V_Href->{min_sfMHoplen},
		"-maxsft:i"	=> \$V_Href->{max_sfMHshift},
		"-sfcnb"	=> \$V_Href->{SVG_sfMH_cirNum},
		# help
		"-h|help"	=> \$V_Href->{HELP},
	);
}

#--- test para and alert ---
sub para_alert{
	return  (   $V_Href->{HELP}
			 || !file_exist( filePath=>$V_Href->{VirusIntegList} )
			 || !-d $V_Href->{output_dir}
			 || !defined $V_Href->{sample_ID}
			 || !file_exist( filePath=>$V_Href->{samtools} )
			 # || !file_exist( filePath=>$V_Href->{host_ref} )
			 # || !file_exist( filePath=>$V_Href->{vref_forItg} )
			 || $V_Href->{FlankSize} < 10
			 || $V_Href->{min_vtMH} < 1
			 || $V_Href->{max_vtMH} < 1
			 || $V_Href->{min_sfMH} < 1
			 || $V_Href->{max_sfMH} < 1
			 || $V_Href->{min_sfMHoplen} < 0
			 || $V_Href->{max_sfMHshift} < 1
			);
}

#--- find micro-homology around virus integrations ---
sub find_vInteg_MH{

	&load_vIteg_list;

	&load_seq( ref_type => 'host' );

	&load_seq( ref_type => 'virus' );

	&search_MH;

	&output_MH_report;

	&draw_MH_figure;
}

#--- load virus Integtions
sub load_vIteg_list{
	# read virus integration list
	my @objPara = qw/ refseg_p5 bkpos_p5 strand_p5 refseg_p3 bkpos_p3 strand_p3 inner_ins /;
	open (VITG, Try_GZ_Read($V_Href->{VirusIntegList})) || die "fail read VirusIntegList: $!\n";
    my $headTagAf = getHeaderTag(header => ($_=<VITG>)); # get header
	while(<VITG>){
		next if /^#/;
        s/\(\d+\)//; # alternative virus bkpos
		my @ele = split;
        my $infoHf = lineInfoToHash(headTagAf => $headTagAf, lineInfo => $_);
        my @para = map{ $_ => $infoHf->{$_} } @objPara;
		my $vInteg_OB = FuseSV::Objects::vInteg_OB->new( @para );
		push @{$V_Href->{vIntegOBs_Aref}}, $vInteg_OB;
		push @{$V_Href->{refsegToVitgOB_Href}->{$vInteg_OB->get_refseg(p=>$_)}->{$_}}, $vInteg_OB for ('p5','p3');
	}
	close VITG;
	# inform
	stout_and_sterr "[INFO]\tread virus integration list OK.\n";
}

#--- load sequence from host ref
sub load_seq{
	# options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	my %parm = @_;
	my $ref_type = $parm{ref_type};
	# which ref
	my $ref_fa = ( $ref_type eq 'host' ? $V_Href->{host_ref} : $V_Href->{vref_forItg} );
	return unless( defined $ref_fa );
	# read ref.fai
	open (FAI, Try_GZ_Read($ref_fa.'.fai')) || die "cannot read ref_fai: $!\n";
	while(<FAI>){
		my ($refseg, $reflen) = (split)[0,1];
		next unless( exists $V_Href->{refsegToVitgOB_Href}->{$refseg} );
		# extract seq
		my $refseg_Href = $V_Href->{refsegToVitgOB_Href}->{$refseg};
		for my $prime (sort keys %$refseg_Href){
			for my $vInteg_OB (@{$refseg_Href->{$prime}}){
				for my $flkType (qw/ junc mkup /){
					my $FlkReg_Href = $vInteg_OB->get_flank_region( p=>$prime, flank=>$V_Href->{FlankSize}, type=>$flkType );
					my $stp = $FlkReg_Href->{stp};
					my $edp = $FlkReg_Href->{edp};
					if( $stp < 1 || $edp > $reflen ){
						warn "<WARN>\t".$vInteg_OB->get_vInteg_info." $flkType flank interval is out of range.\n"
								  ."\t $refseg: [$stp, $edp]\n";
						last;
					}
					(my $seq = uc(`$V_Href->{samtools} faidx $ref_fa $refseg:$stp-$edp | grep -v '>'`)) =~ s/\s//g;
					($seq = reverse $seq) =~ tr/ACGT/TGCA/ if( $FlkReg_Href->{rvp} );
					$vInteg_OB->load_flank_seq( p=>$prime, type=>$flkType, seq=>$seq );
				}
			}
		}
	}
	close FAI;
	# inform
	stout_and_sterr "[INFO]\tload $ref_type flank sequences of vIntegs OK.\n";
}

#--- search MH from flank seq
sub search_MH{
	# find MH for each vInteg
	for my $vInteg_OB (@{$V_Href->{vIntegOBs_Aref}}){
		unless( $vInteg_OB->has_complete_flkseq ){
			warn "<WARN>\t".$vInteg_OB->get_vInteg_info." lacks completed flank seg.\n";
			next;
		}
		# vertial MH
		$vInteg_OB->find_flkseq_MH(
			MH_type => 'vtMH',
			min_MHlen => $V_Href->{min_vtMH},
			max_MHlen => $V_Href->{max_vtMH}
		);
		# shift MH
		$vInteg_OB->find_flkseq_MH(
			MH_type=>'sfMH',
			min_MHlen => $V_Href->{min_sfMH},
			max_MHlen => $V_Href->{max_sfMH},
			min_sfMHoplen => $V_Href->{min_sfMHoplen},
			max_sfMHshift => $V_Href->{max_sfMHshift}
		);
		# count valid cases
		$V_Href->{validCaseForMH} ++;
	}
	# inform
	stout_and_sterr "[INFO]\tsearch MH in flank region of vInteg OK.\n";
}

#--- output MH report ---
sub output_MH_report{
	my $MHreport = GetPath( filekey => 'FindMH_report' );
	open (MHRT, Try_GZ_Write($MHreport)) || die "fail write mh_report: $!\n";
	print MHRT join("\t", '#DOCS', 'FORMAT=MHRTv1.0;DATE='.strftime("%Y%m%d",localtime())."\n");
	for my $vInteg_OB (@{$V_Href->{vIntegOBs_Aref}}){
		next unless( $vInteg_OB->has_complete_flkseq );
		print MHRT join("\t", '#INFO', $vInteg_OB->get_vInteg_info).";FLANK=$V_Href->{FlankSize}\n";
		print MHRT join("\t", '#vtMH', $_)."\n" for @{ $vInteg_OB->get_MH_info(type=>'vtMH') };
		print MHRT join("\t", '#sfMH', $_)."\n" for @{ $vInteg_OB->get_MH_info(type=>'sfMH') };
		print MHRT "$_\n" for @{ $vInteg_OB->get_MH_landscape };
	}
	close MHRT;
	# inform
	stout_and_sterr "[INFO]\tgenerate MH report OK.\n";
}

#--- draw MH SVG ---
sub draw_MH_figure{
	# figure file
	my $MH_SVGfigure = GetPath( filekey => 'FindMH_SVG_figure' );
	`rm $MH_SVGfigure` if( -e $MH_SVGfigure );
	# no valid case
	if( $V_Href->{validCaseForMH} == 0 ){
		stout_and_sterr "<WARN>\tno valid vInteg case to show MH figure!\n";
		return;
	}
	# draw
	$V_Href->{general_BaseX} += 30;
	my $FigureWidth = $V_Href->{general_BaseX} + $V_Href->{FlankSize} * 3 * $V_Href->{SVG_SeqBaseRectWidth};
	my $oneCaseHeight = 10 * $V_Href->{SVG_SeqBaseRectHeight} + 10 * ($V_Href->{SVG_sfMH_cirNum}?$V_Href->{SVG_sfMH_radius}:$V_Href->{SVG_sfMH_dotR});
	my $FigureHeight= $V_Href->{general_BaseY} + $V_Href->{validCaseForMH} * ($oneCaseHeight+$V_Href->{SVG_vItg_gap});
	$V_Href->{general_SVG_obj} = SVG->new( width=>$FigureWidth, height=>$FigureHeight );
	for my $vInteg_OB (@{$V_Href->{vIntegOBs_Aref}}){
		next unless( $vInteg_OB->has_complete_flkseq );
		$vInteg_OB->draw_MH;
		$V_Href->{general_BaseY} += $V_Href->{SVG_vItg_gap};
	}
	# add figure information
	show_figure_information(
		fNameAref => [ 'Micro-Homology around', "Virus-Integ Flank($V_Href->{FlankSize})" ],
		fVersion => $VERSION,
		fAuthor => $AUTHOR,
		fAffili => 'CityU of Hong Kong'
	);
	# generate figure file
	open (MHF,">$MH_SVGfigure") || die"fail generate SVG: $!\n";
	print MHF $V_Href->{general_SVG_obj}->xmlify;
	close MHF;
	# inform
	stout_and_sterr "[INFO]\tgenerate MH figure OK.\n";
}

#---
1; ## tell the perl script the successful access of this module.
