package FuseSV::Virus_Integ::LocalHaplotypeDetect::PairWiseAnchorLinkStat;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr file_exist /;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycleToBioContig qw/ deal_SEG_info deal_JUNC_info check_segments /;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  pairwise_anchors_link_stat
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'PairWiseAnchorLinkStat';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-05-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						pairwise_anchors_link_stat
						read_config
						read_Junc_marker_list
						collect_marker_HostSegPosAnchor
						calculate_PWanchor_link
						generate_output
					 /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Get links (i.e., shared markers) of pair-wise anchors.
	 Markers could be barcodes (10x) or PE-id (Hi-C, MinIon and PacBio).

	 Usage:   perl $V_Href->{MainName} PWanStat <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -config [s]  unit-cycle config from local haplotype algorithm. <required>
	    -pwlink [s]  output file recodes pair-wise link counts. <required>

	   # Software Required #
	    -stl    [s]  SamTools, minimum version: 1.3. [optional]
	                  Note: effective with '-hbam'.

	   # Options #
	    -hbam   [s]  indexed host-ref bam file. [optional]
	    -stype  [s]  sequencing type of long-range DNA. ['10x']
	                  Note: valid types are '10x', 'Hi-C', MinIon', and 'PacBio'.
	    -jcmkl  [s]  marker list of junction. [optional]
	                  Note: two column format, JUNC_ID and marker_ID.
	    -segpnm [i]  number of pos-type anchor in one segment. [1]
	                  Note: 1) use this number to divide the segment equally.
	                        2) so, default value 1 means the middle-point.
	    -minpds [i]  min distance of one pos-type anchor from edge of segment. [100]
	    -matrix      turn on matrix mode of output. [disabeld]
	                  Note: non-matrix mode is for func 'UCYCtoLM'.

	    -h|help      Display this help info.

	 Version:
	    $VERSION at $DATE

	 Author:
	    $AUTHOR ($EMAIL)
 \n";
}
	    # -vbam   [s]  indexed virus-ref bam file. [optional]
	    # -vref   [s]  name of virus-ref in 'vbam' file. [undef]
	    #               Note: must provided with 'vbam'.

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
			#--- input/output
			[ config => undef ],
			[ output => undef ],
			# software
			[ samtools => '' ],
			#--- options
			[ hostBam => {} ],
			[ virusBam => {} ],
			[ virusRef => undef ],
			[ seqType => '10x' ],
			[ junc_marker_list => undef ],
			[ seg_posAnchor_num => 1 ],
			[ min_edgeDist_posAnchor => 100 ],
			[ output_matrix => 0 ],
			#--- intermediate variants
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB module
			## Group_OB contains Solution_OB -> Allele_OB -> LocalMap_OB
			[ Group => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB module
			[ Segment => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB module
			[ Junc => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB module
			[ Pos => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB module
			[ PWanchor => {} ],
			# store anchors
			[ Anchor => {} ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['config'],
								  ['samtools'],
								  ['junc_marker_list'],
								  ['output']  ] ]
		);
}

#--- get options from command line ---
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# objects constructions, as these objects belong to the input options
	$V_Href->{hostBam}  = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'case', filepath => undef );
	$V_Href->{virusBam} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new( tissue => 'case', filepath => undef );

	# get options
	GetOptions(
		# input/output
		"-config:s"	=> \$V_Href->{config},
		"-pwlink:s"	=> \$V_Href->{output},
		# software
		"-stl:s"	=> \$V_Href->{samtools},
		# options
		"-hbam:s"	=> \$V_Href->{hostBam}->{filepath},
		"-vbam:s"	=> \$V_Href->{virusBam}->{filepath},
		"-vref:s"	=> \$V_Href->{virusRef},
		"-stype:s"	=> \$V_Href->{seqType},
		"-jcmkl:s"	=> \$V_Href->{junc_marker_list},
		"-segpnm:i"	=> \$V_Href->{seg_posAnchor_num},
		"-minpds:i"	=> \$V_Href->{min_edgeDist_posAnchor},
		"-matrix"	=> \$V_Href->{output_matrix},
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
			 || !file_exist( filePath=>$V_Href->{config} )
			 || !defined $V_Href->{output}
			 || (   defined $V_Href->{hostBam}->{filepath}
			 	 && !file_exist( filePath=>$V_Href->{hostBam}->{filepath} )
			 	)
			 || (   defined $V_Href->{virusBam}->{filepath}
			 	 && !file_exist( filePath=>$V_Href->{virusBam}->{filepath} )
			 	)
			 || (   defined $V_Href->{junc_marker_list}
			 	 && !file_exist( filePath=>$V_Href->{junc_marker_list} )
			 	)
			 || (   !file_exist( filePath=>$V_Href->{samtools} )
			 	 && (   defined $V_Href->{hostBam}->{filepath}
			 		 || defined $V_Href->{virusBam}->{filepath}
			 	 	)
			 	)
			 || (    defined $V_Href->{virusBam}->{filepath}
			 	 && !defined $V_Href->{virusRef}
			 	)
			 || (   !defined $V_Href->{hostBam}->{filepath}
			 	 && !defined $V_Href->{virusBam}->{filepath}
			 	 && !defined $V_Href->{junc_marker_list}
			 	)
			 || $V_Href->{seg_posAnchor_num} < 1
			 || $V_Href->{min_edgeDist_posAnchor} < 1
			);
}

#--- calculate links of pair-wise anchors ---
sub pairwise_anchors_link_stat{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# read config
	&read_config;

	# check the availability of all segments
	check_segments;

	# load JUNC marker file
	&read_Junc_marker_list;

	# collect marker of pos-anchor in host segments
	&collect_marker_HostSegPosAnchor;

	# calculate pair-wise anchors' link
	&calculate_PWanchor_link;

	# write output
	&generate_output;
}

#--- read unit-cycle config from local haplotype algorithm ---
sub read_config{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	open (CONFIG, Try_GZ_Read($V_Href->{config})) || die "fail read UnitCycle_configFile: $!\n";
	while(<CONFIG>){
		chomp;
		next if(/^\#/); # should be comments
		# diff type content
		if(/^\@DOC\s/){ # skip document information
			next;
		}
		elsif(/^\@[HV]SEG\s/){ # segment
			deal_SEG_info( content => $_ );
		}
		elsif(/^\@JUNC\s/){ # junction
			deal_JUNC_info( content => $_ );
		}
		elsif(/^\@UCYC\s/){ # skip unit-cycle
			next;
		}
		elsif(/^SOLUTION\s/){ # skip solution
			next;
		}
	}
	close CONFIG;

	# inform
	stout_and_sterr `date`."[FILE]\tRead config file OK.\n";
}

#--- load markers of Junc ---
sub read_Junc_marker_list{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Anchor_Href = $V_Href->{Anchor};

	return if( !defined $V_Href->{junc_marker_list} );

	open (JCML, Try_GZ_Read($V_Href->{junc_marker_list})) || die "fail read junc_marker_list: $!\n";
	while(<JCML>){
		chomp;
		next if(/^\#/); # should be comments
		my ($juncID, $markers) = (split)[0,1];
		$Anchor_Href->{$juncID}->{marker}->{$_} ++ for split /[,;]+/, $markers;
	}
	close JCML;

	# examine juncID
	for my $juncID (keys %$Anchor_Href){
		if( !exists $V_Href->{Junc}->{$juncID} ){
			warn_and_exit "<ERROR>\tJuncID ($juncID) is not valid in junc_marker_list.\n";
		}
		else{
			# use juncID as anchorID
			$Anchor_Href->{$juncID}->{type} = 'JUNC';
			$Anchor_Href->{$juncID}->{obj} = $V_Href->{Junc}->{$juncID};
		}
	}

	# inform
	stout_and_sterr "[FILE]\tRead junc_marker list OK.\n";
}

#--- collect marker of pos-anchor in host segments ---
sub collect_marker_HostSegPosAnchor{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{Segment};
	my $Anchor_Href = $V_Href->{Anchor};

	return unless(defined $V_Href->{hostBam}->{filepath});

	my %RemainSegOB;
	my @UnitLength;
	# inner segments
	for my $segID (sort keys %$Segment_Href){
		my $segOB = $Segment_Href->{$segID};
		next if($segOB->get_SpecTag ne 'H');
		my $groupOB = $segOB->get_groupOB;
		if($segID eq $groupOB->get_sourceSegID){
			push @{$RemainSegOB{Source}}, $segOB;
		}
		elsif($segID eq $groupOB->get_sinkSegID){
			push @{$RemainSegOB{Sink}}, $segOB;
		}
		else{
			my $unitLength = int($segOB->get_length / ($V_Href->{seg_posAnchor_num}+1));
			my $seg_chr = $segOB->get_chr;
			my $seg_stpos = $segOB->get_stpos;
			for my $i (1 .. $V_Href->{seg_posAnchor_num}){
				my $pos = $seg_stpos + $i * $unitLength - 1;
				my $posID = "$seg_chr:$pos";
				my @parm = ( chr => $seg_chr, pos => $pos, seg_OB => $segOB, minEdgeD => $V_Href->{min_edgeDist_posAnchor} );
				$V_Href->{Pos}->{$posID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB->new( @parm, keepError => 1 );
				if( ! $V_Href->{Pos}->{$posID}->is_valid ){
					stout_and_sterr "<WARN>\tPass $posID in host segment ($segID) as min_edgeDist filter.\n";
					delete $V_Href->{Pos}->{$posID};
					next;
				}
				else{
					push @UnitLength, $unitLength;
					my $anchorID = "$segID,$posID";
					$Anchor_Href->{$anchorID}->{type} = 'POS';
					$Anchor_Href->{$anchorID}->{obj} = $V_Href->{Pos}->{$posID};
					my @parm = ( chr => $seg_chr, pos => $pos, marker => $V_Href->{seqType}, samtools => $V_Href->{samtools} );
					$Anchor_Href->{$anchorID}->{marker} = $V_Href->{hostBam}->get_pos_marker_stat( @parm );
					stout_and_sterr "[INFO]\tGet markers of $posID in host segment ($segID).\n";
				}
			}
		}
	}

	# check pos anchors
	if( scalar(@UnitLength) == 0 ){
		stout_and_sterr "<WARN>\tNot valid Pos anchor found in inner segments of Local Map.\n";
		return;
	}

	# source and sink segments
	my $avg_UnitLen = max(int(sum(@UnitLength) / scalar(@UnitLength)), $V_Href->{min_edgeDist_posAnchor});
	for my $type (keys %RemainSegOB){
		for my $segOB (@{$RemainSegOB{$type}}){
			my $segID = $segOB->get_ID;
			my $seg_chr = $segOB->get_chr;
			my $posShift = min($avg_UnitLen, int($segOB->get_length/2));
			my $pos = ( $type eq 'Source' ? ($segOB->get_edpos - $posShift) : ($segOB->get_stpos + $posShift) );
			my $posID = "$seg_chr:$pos";
			my @parm = ( chr => $seg_chr, pos => $pos, seg_OB => $segOB, minEdgeD => $V_Href->{min_edgeDist_posAnchor} );
			$V_Href->{Pos}->{$posID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB->new( @parm, keepError => 1 );
			if( ! $V_Href->{Pos}->{$posID}->is_valid ){
				stout_and_sterr "<WARN>\tPass $posID in host segment ($segID) as min_edgeDist filter.\n";
				delete $V_Href->{Pos}->{$posID};
				next;
			}
			else{
				my $anchorID = "$segID,$posID";
				$Anchor_Href->{$anchorID}->{type} = 'POS';
				$Anchor_Href->{$anchorID}->{obj} = $V_Href->{Pos}->{$posID};
				my @parm = ( chr => $seg_chr, pos => $pos, marker => $V_Href->{seqType}, samtools => $V_Href->{samtools} );
				$Anchor_Href->{$anchorID}->{marker} = $V_Href->{hostBam}->get_pos_marker_stat( @parm );
				stout_and_sterr "[INFO]\tGet markers of $seg_chr:$pos in host segment ($segID).\n";
			}
		}
	}

	# inform
	stout_and_sterr "[INFO]\tGet markers of pos anchors OK.\n";
}

#--- calculate pair-wise anchors' link ---
sub calculate_PWanchor_link{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Anchor_Href = $V_Href->{Anchor};

	my @anchorID = sort keys %$Anchor_Href;
	return if( scalar(@anchorID) == 0 );

	# pair-wise anchor
	for my $i (0 .. $#anchorID-1){
		my $i_anchorID = $anchorID[$i];
		my $i_anchorType = $Anchor_Href->{$i_anchorID}->{type};
		my $i_anchorOB = $Anchor_Href->{$i_anchorID}->{obj};
		my $i_anchorMKHref = $Anchor_Href->{$i_anchorID}->{marker};
		for my $j ($i+1 .. $#anchorID){
			my $j_anchorID = $anchorID[$j];
			my $j_anchorType = $Anchor_Href->{$j_anchorID}->{type};
			my $j_anchorOB = $Anchor_Href->{$j_anchorID}->{obj};
			my $j_anchorMKHref = $Anchor_Href->{$j_anchorID}->{marker};
			# for Pos Anchor cannot from same segment
			if(    $i_anchorType eq 'POS' && $j_anchorType eq 'POS'
				&& $i_anchorOB->get_segID eq $j_anchorOB->get_segID
			){
				next;
			}
			# shared marker
			my $LinkCount = grep exists $i_anchorMKHref->{$_}, keys %$j_anchorMKHref;
			# record
			$Anchor_Href->{$i_anchorID}->{links}->{$j_anchorID} = $LinkCount;
			# $Anchor_Href->{$j_anchorID}->{links}->{$i_anchorID} = $LinkCount;
		}
	}

	# inform
	stout_and_sterr "[INFO]\tCount shared markers of pair-wise anchors as link count OK.\n";
}

#--- output pair-wise anchor links ---
sub generate_output{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Anchor_Href = $V_Href->{Anchor};

	my @anchorID = sort keys %$Anchor_Href;
	return if( scalar(@anchorID) == 0 );

	# marker links
	open (PWLINK, Try_GZ_Write($V_Href->{output})) || die "fail write pair-wise link file: $!\n";
	if( $V_Href->{output_matrix} ){
		print PWLINK '#'.join("\t", 'anchor_i/j', @anchorID)."\n";
		for my $i_anchorID (@anchorID){
			my @content;
			push @content, $i_anchorID;
			for my $j_anchorID (@anchorID){
				if( $i_anchorID eq $j_anchorID ){
					push @content, '-';
				}
				else{
					my $link =    $Anchor_Href->{$i_anchorID}->{links}->{$j_anchorID}
							   || $Anchor_Href->{$j_anchorID}->{links}->{$i_anchorID};
					push @content, ( defined $link ? $link : '-');
				}
			}
			print PWLINK join("\t", @content)."\n";
		}
	}
	else{
		print PWLINK '#'.join("\t", qw/ anchor_i anchor_j link /)."\n";
		for my $i_anchorID (@anchorID){
			next unless(exists $Anchor_Href->{$i_anchorID}->{links});
			for my $j_anchorID (sort keys %{$Anchor_Href->{$i_anchorID}->{links}}){
				print PWLINK join("\t", $i_anchorID, $j_anchorID, $Anchor_Href->{$i_anchorID}->{links}->{$j_anchorID})."\n";
			}
		}
	}
	close PWLINK;

	# marker list
	my $AnchorMarkerList = $V_Href->{output}.'marker.list';
	open (MKLIST, Try_GZ_Write($AnchorMarkerList)) || die "fail write anchor marker list: $!\n";
	print MKLIST '#'.join("\t", qw/ anchorID markerID /)."\n";
	for my $anchorID (@anchorID){
		next unless(exists $Anchor_Href->{$anchorID}->{marker});
		print MKLIST join("\t",$anchorID,$_)."\n" for sort keys %{$Anchor_Href->{$anchorID}->{marker}};
	}
	close MKLIST;

	# inform
	stout_and_sterr `date`."[FILE]\tGenerate pair-wise link file and anchor marker list OK.\n"
								."\t$V_Href->{output}\n"
								."\t$AnchorMarkerList\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
