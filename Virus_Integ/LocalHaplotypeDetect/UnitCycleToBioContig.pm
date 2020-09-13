package FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycleToBioContig;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max sum /;
use SOAPfuse::OpenFile qw/ Try_GZ_Read Try_GZ_Write /;
use SOAPfuse::General_Operation qw/ warn_and_exit stout_and_sterr file_exist /;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB;
use FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  merge_unitcycle_to_biocontig
			  deal_SEG_info
			  deal_JUNC_info
			  check_segments
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'UnitCycleToBioContig';
#----- version --------
$VERSION = "0.13";
$DATE = '2018-05-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						merge_unitcycle_to_biocontig
						read_config
						deal_DOC_info
						deal_SEG_info
						deal_JUNC_info
						deal_UCYC_info
						deal_SOLUTION_info
						link_Vseg_allGroup
						check_segments
						load_long_range_info
						load_longRangeLenDist
						Load_PWanchorLinkCount
						resolve_localmap
						generate_output
					 /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} UCYCtoLM <[Options]>
	 
	 Options:

	   # Inputs and Outputs #
	    -config [s]  unit-cycle config from local haplotype algorithm. <required>
	    -output [s]  output local map file. <required>

	   # Options #
	    -randtm [i]  how many times to try construct local-map(s) in random mode. [0]
	                  Note: 0 means simplest local-map mode.
	    -segstr [s]  local-map must have such segments string(s) in random mode.
	                  Note: 1) use as, -segstr 'H2+;V3-;H3+;H5+'
	                        2) allows multiple times: -segstr 'xx' -segstr 'xx' ...
	                        3) allows omit seg-orit (+/-), e.g., 'H2;V3-;H3;H5+'
	    -solu_y [s]  select these SOLUTION(s) to analysis. [disabled]
	                  Note: 1) use like, -solu_y 1,3,6
	                        2) find 'SOLUTION NO.' in config file.
	    -solu_n [s]  avoid these SOLUTION(s) from analysis. [disabled]
	                  Note: 1) use like, -solu_n 2,7,10
	                        2) it has higher priority than '-solu_y'.
	    -usefvgm     use FVGM (Free Viral GenoMe, if has) in random mode. [disabled]
	    -pwlink [s]  the file recodes pair-wise anchor link counts. [optional]
	                  Note: please check FuseSV online instructions to know about.
	    -lrdist [s]  file recodes distribution of long-range DNA length. [optional]
	                  Note: 1) e.g., the barcode cover length in 10x-seq data.
	                        2) effective with '-pwlink'.
	                        3) check FuseSV online instructions to know about.
	    -corrmd [s]  method to calculate pw-anchor links' correlations. [p]
	                  Note: 'p' for pearson, 's' for spearman.
	    -stype  [s]  sequencing type of long-range DNA. ['10x']
	                  Note: valid types are '10x', 'Hi-C', MinIon', and 'PacBio'.
	    -depth  [i]  simualted sequencing depth of the long-range DNA. [50]
	    -perlen [i]  read length of paired-end (PE) sequencing in simulation. [150]
	                  Note: 1) effective with '-stype' as 'Hi-C'.
	                        2) use it for both ends of PE-reads.
	    -fork   [i]  to run simualted sequencing with N forks in parallel. [1]

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
			#--- input/output
			[ config => undef ],
			[ output => undef ],
			#--- options
			[ randTime => 0 ], # random-times
			[ SegStrLMhas => [] ],
			[ probCombCount => 3 ], # hidden
			[ useFVGM => 0 ],  # add FVGM during random mode
			[ SolutionNO_select => '' ], # to Hash next
			[ SolutionNO_giveup => '' ], # to Hash next
			# long-range simulation work
			[ pairWiseLinkFile => undef ],
			[ LongRangeLenDistFile => undef ],
			[ CorrMethod => 'p' ],
			[ seqType => '10x' ],
			[ aimDepth => 50 ],
			[ readLen => 150 ],
			[ forkNum => 1 ],
			#--- intermediate variants
			[ ConfigHeader => [] ], # for output
			[ GroupOrder => {} ], # for output
			[ SampleID => undef ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB module
			## Group_OB contains Solution_OB -> Allele_OB -> LocalMap_OB
			[ Group => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB module
			[ Segment => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB module
			[ Junc => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB module
			[ UnitCycle => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB module
			[ Pos => {} ],
			# check structure in FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB module
			[ PWanchor => {} ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['config'],
								  ['pairWiseLinkFile'],
								  ['LongRangeLenDistFile'],
								  ['output']  ] ]
		);
}

#--- get options from command line ---
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-config:s"	=> \$V_Href->{config},
		"-output:s"	=> \$V_Href->{output},
		# options
		"-randtm:i"	=> \$V_Href->{randTime},
		"-segstr:s"	=> \@{$V_Href->{SegStrLMhas}},
		"-pbcmbc:i"	=> \$V_Href->{probCombCount}, # hidden
		"-usefvgm"	=> \$V_Href->{useFVGM},
		"-pwlink:s"	=> \$V_Href->{pairWiseLinkFile},
		"-lrdist:s"	=> \$V_Href->{LongRangeLenDistFile},
		"-corrmd:s"	=> \$V_Href->{CorrMethod},
		"-stype:s"	=> \$V_Href->{seqType},
		"-depth:i"	=> \$V_Href->{aimDepth},
		"-perlen:i"	=> \$V_Href->{readLen},
		"-fork:i"	=> \$V_Href->{forkNum},
		"-solu_y:s"	=> \$V_Href->{SolutionNO_select},
		"-solu_n:s"	=> \$V_Href->{SolutionNO_giveup},
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
			 || $V_Href->{randTime} < 0
			 || $V_Href->{probCombCount} < 1
			 || (   defined $V_Href->{pairWiseLinkFile}
			 	 && !file_exist( filePath=>$V_Href->{pairWiseLinkFile} )
			 	)
			 || (   defined $V_Href->{LongRangeLenDistFile}
			 	 && (   !file_exist( filePath=>$V_Href->{LongRangeLenDistFile} )
			 	 	 || !file_exist( filePath=>$V_Href->{pairWiseLinkFile} )
			 	 	 || $V_Href->{CorrMethod} !~ /^[ps]$/i
			 	 	)
			 	)
			);
}

#--- merge unit-cycle to bio-contig ---
sub merge_unitcycle_to_biocontig{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# read config
	&read_config;

	# check the availability of all segments
	&check_segments;

	# link all virus segment to all groups
	&link_Vseg_allGroup;

	# read long-range info
	&load_long_range_info;

	# construct local-map(s)
	&resolve_localmap;

	# write output
	&generate_output;
}

#--- read unit-cycle config from local haplotype algorithm ---
sub read_config{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# solution selection
	$V_Href->{$_} = { map {($_,1)} grep /^\d+$/, (split /,/, $V_Href->{$_}) } for qw/ SolutionNO_select SolutionNO_giveup /;
	$V_Href->{SolutionNO_select} = undef if(scalar(keys %{$V_Href->{SolutionNO_select}}) == 0);

	open (CONFIG, Try_GZ_Read($V_Href->{config})) || die "fail read UnitCycle_configFile: $!\n";
	while(<CONFIG>){
		chomp;
		next if(/^\#/); # should be comments
		# store header for output
		if(/^\@/){
			push @{$V_Href->{ConfigHeader}}, $_;
		}
		# diff type content
		if(/^\@DOC\s/){ # document information
			&deal_DOC_info( content => $_ );
		}
		elsif(/^\@[HV]SEG\s/){ # segment
			&deal_SEG_info( content => $_ );
		}
		elsif(/^\@JUNC\s/){ # junction
			&deal_JUNC_info( content => $_ );
		}
		elsif(/^\@UCYC\s/){ # unit-cycle
			&deal_UCYC_info( content => $_ );
		}
		elsif(/^SOLUTION\s/){ # solution
			&deal_SOLUTION_info( content => $_ );
		}
	}
	close CONFIG;

	# inform
	stout_and_sterr `date`."[FILE]\tRead config file OK.\n";
}

#--- deal DOC line content ---
sub deal_DOC_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $content = $parm{content};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$_ = $content;
	#@DOC	VERSION=1.0.0;SAMPLE_ID=T014
	($V_Href->{SampleID}) = (/SAMPLE_ID=([^;]+)/);
	if( !defined $V_Href->{SampleID} ){
		$$V_Href->{SampleID} = 'NULL';
		stout_and_sterr "<WARN>\treset sample-id as 'NULL'.\n";
	}
}

#--- deal SEG line content ---
sub deal_SEG_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $content = $parm{content};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$_ = $content;
	# @HSEG	ID=H1;GROUP=DEF;INTERVAL=chr8:128300000-128364420;CN_REF=xx;CN_ORIG=3.00;CN_ALT=xx;CN_ALT_LP=3.00;LP_WEIGHT=1.00;LP_OFFSET%=17.19%;NOTES=SOURCE
	# @VSEG	ID=V1;INTERVAL=HPV18:1-782;CN_ORIG=23.00;CN_ALT=xx;CN_ALT_LP=23.00;LP_WEIGHT=1.00;LP_OFFSET%=0.00%
	my ($segID)     = /\bID=([^;]+)/;
	my ($groupID)   = /\bGROUP=(\d+)/;
	my ($region)    = /\bINTERVAL=([^;]+)/;
	my ($orig_cn)   = /\bCN_ORIG=([^;]+)/;
	my ($cn_ref)    = /\bCN_REF=([^;]+)/;
	my ($cn_alt)    = /\bCN_ALT=([^;]+)/;
	my ($cn_alt_LP) = /\bCN_ALT_LP=([^;]+)/;
	my ($LP_weight) = /\bLP_WEIGHT=([^;]+)/;
	my ($LP_offset) = /\bLP_OFFSET\%=([^;]+)/;
	# parmeters for creating segment OB
	my @parm = ( ID => $segID, region => $region, orig_cn => $orig_cn, cn_alt => $cn_alt, cn_alt_LP => $cn_alt_LP, LP_weight => $LP_weight, LP_offset => $LP_offset );
	push @parm, ( cn_ref => $cn_ref ) if( defined $cn_ref );
	# for group
	if( defined $groupID ){
		# create group object
		if( !exists $V_Href->{Group}->{$groupID} ){
			$V_Href->{Group}->{$groupID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Group_OB->new( ID => $groupID );
		}
		push @parm, ( groupOB => $V_Href->{Group}->{$groupID} );
	}
	elsif( /^\@HSEG\s/ ){ # host-seg must have GROUP tag
		warn_and_exit "<ERROR>\tHost-Segment lacks digital 'GROUP' tag (e.g., GROUP=1;).\n"
					 ."Content:$_\n";
	}
	# create segment OB
	$V_Href->{Segment}->{$segID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Segment_OB->new( @parm );
	# link segment OB with group OB
	if( defined $groupID ){
		# link this segment with this group
		$V_Href->{Group}->{$groupID}->load_SegOB( seg_OB => $V_Href->{Segment}->{$segID} );
		# source / sink
		if(/\bNOTES=SOURCE\b/){
			$V_Href->{Group}->{$groupID}->load_SourceSegOB( seg_OB => $V_Href->{Segment}->{$segID} );
		}
		elsif(/\bNOTES=SINK\b/){
			$V_Href->{Group}->{$groupID}->load_SinkSegOB( seg_OB => $V_Href->{Segment}->{$segID} );
		}
	}
}

#--- deal JUNC line content ---
sub deal_JUNC_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $content = $parm{content};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$_ = $content;
	#@JUNC   ID=JUNC1;SEGLINK=V1+,H4-;CN_ORIG=3.32;CN_ALT=3.32;CN_ALT_LP=2.00;LP_WEIGHT=0.20;LP_OFFSET%=-79.52%;W=1;REFSEG_5P=HPV18;BKPOS_5P=406;STRD_5P=+;REFSEG_3P=chr8;BKPOS_3P=128283649;STRD_3P=-;NCPR=0;RAWJRC=45;PTJRC=45
	my ($juncID) = /\bID=([^;]+)/;
	my ($SegLink) = /\bSEGLINK=([^;]+)/;
	my (%segID, %segOrit, %segOB);
	($segID{5}, $segOrit{5}, $segID{3}, $segOrit{3}) = ($SegLink =~ /^([HV]\d+)([\+\-]),([HV]\d+)([\+\-])$/);
	for my $prime (qw/ 5 3 /){
		my $segID = $segID{$prime};
		if( !defined $segID ){
			warn_and_exit "<ERROR>\tCannot analyze SegLink ($SegLink) of Junc $juncID.\n";
		}
		if( !exists $V_Href->{Segment}->{$segID} ){
			warn_and_exit "<ERROR>\tCannot find the $prime-prime SegOB of JUNC $juncID: $SegLink.\n";
		}
		$segOB{$prime} = $V_Href->{Segment}->{$segID};
	}
	# create new junc_OB
	my @parm = ( ID => $juncID, SegLink => $SegLink,
				 p5_segOB => $segOB{5}, p5_segOrit => $segOrit{5},
				 p3_segOB => $segOB{3}, p3_segOrit => $segOrit{3}  );
	$V_Href->{Junc}->{$juncID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Junc_OB->new( @parm );
	# link JUNC with group
	my $groupOB = $V_Href->{Junc}->{$juncID}->get_groupOB;
	if( defined $groupOB ){ # this junc must have 'H' seg
		$groupOB->load_JuncOB( junc_OB => $V_Href->{Junc}->{$juncID} );
	}
}

#--- deal UCYC line content ---
sub deal_UCYC_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $content = $parm{content};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$_ = $content;
	# @UCYC	ID=UCYC1;SEGSTR=H3+,H4+,V2-,V1-,V4-,H3+
	# @UCYC	ID=UCYC5;SEGSTR=V1+,V2+,V3+,V4+,V1+;FVGM_ID=1
	my ($UcycID)  = /\bID=([^;]+)/;
	my ($segComb) = /\bSEGSTR=([^;]+)/;
	my ($fvgmID)  = /\bFVGM_ID=([^;]+)/;
	if( !defined $UcycID || !defined $segComb ){
		warn_and_exit "<ERROR>\tUnable to recognize UCYC information.\n"
					 ."Content:$_\n";
	}
	if( exists $V_Href->{UnitCycle}->{$UcycID} ){
		warn_and_exit "<ERROR>\tUnit-cycle $UcycID is duplicated. Not allowed.\n"
					 ."Content:$_\n";
	}
	my @parm = ( ID => $UcycID, segComb => $segComb, segHref => $V_Href->{Segment}, strOpmz => ($V_Href->{randTime}>0?0:1) );
	push @parm, ( fvgmID => $fvgmID ) if( defined $fvgmID );
	$V_Href->{UnitCycle}->{$UcycID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::UnitCycle_OB->new( @parm );
	$V_Href->{UnitCycle}->{$UcycID}->update_info;
}

#--- deal SOLUTION line content ---
sub deal_SOLUTION_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $content = $parm{content};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	$_ = $content;
	# SOLUTION	NO=1;GROUP=DEF;SOLUTION_UCYC=UCYC1:1,UCYC2:26,UCYC3:3,UCYC4:22,UCYC5:8;UCYC6:5;ALLELE_NO=1;ALLELE_CN=2;ALLELE_UCYC=UCYC2:13,UCYC3:1,UCYC4:11;ALLELE_FVGM_UCYC=UCYC5:2;UCYC6:1
	my ($NO) = /\bNO=(\d+)/;
	# filter
	if(    exists $V_Href->{SolutionNO_giveup}->{$NO}
		|| (   defined $V_Href->{SolutionNO_select}
		    && !exists $V_Href->{SolutionNO_select}->{$NO}
		   )
	){
		# inform
		stout_and_sterr "<WARN>\tSkip solution $NO according to options ('-solu_y', '-solu_n').\n";
		return;
	}
	# more info
	my ($groupID) = /\bGROUP=([^;]+)/;
	my ($solution_UCYC) = /\bSOLUTION_UCYC=([^;]+)/;
	my ($allele_NO) = /\bALLELE_NO=(\d+)/;
	my ($allele_cn) = /\bALLELE_CN=(\d+)/; # \d+
	my ($allele_UCYC) = /\bALLELE_UCYC=([^;]+)/;
	my ($allele_FVGM_UCYC) = /\bALLELE_FVGM_UCYC=([^;]+)/;
	# must have GROUP tag
	if(    !defined $groupID
		|| !exists $V_Href->{Group}->{$groupID}
	){
		warn_and_exit "<ERROR>\tSolution lacks valid 'GROUP' tag.\n"
					 ."Content:$_\n";
	}
	# record
	my $group_OB = $V_Href->{Group}->{$groupID};
	my $order_NO = scalar(keys %{$V_Href->{GroupOrder}}) + 1;
	if( !exists $V_Href->{GroupOrder}->{$groupID} ){
		$V_Href->{GroupOrder}->{$groupID} = { order_NO => $order_NO, group_OB => $group_OB };
	}
	# load solution+allele information
	my @parm = ( solution_NO => $NO, solution_UCYC_countStr => $solution_UCYC,
				 allele_NO => $allele_NO, allele_UCYC_countStr => $allele_UCYC, allele_cn => $allele_cn );
	push @parm, ( allele_FVGM_UCYC_countStr => $allele_FVGM_UCYC ) if( defined $allele_FVGM_UCYC && $allele_FVGM_UCYC ne 'N/A' );
	$group_OB->load_solution_allele( @parm, UCYC_Href => $V_Href->{UnitCycle}, LineContent => $_ );
}

#--- link all virus segments to all groups
sub link_Vseg_allGroup{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	my $inform_bool = 0;
	for my $segOB (values %{$V_Href->{Segment}}){
		next if($segOB->get_SpecTag ne 'V');
		for my $groupOB (values %{$V_Href->{Group}}){
			$groupOB->load_SegOB( seg_OB => $segOB );
		}
		$inform_bool = 1;
	}

	# inform
	if( $inform_bool ){
		stout_and_sterr "[INFO]\tLink all virus segments to all groups OK.\n";
	}
}

#--- check all segments, make sure the availability ---
sub check_segments{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Segment_Href = $V_Href->{Segment};

	# semgents could not have overlap with each other
	my @segID = sort keys %$Segment_Href;
	for my $idx_a ( 0 .. $#segID-1 ){
		my $segID_a = $segID[$idx_a];
		my $segOB_a = $Segment_Href->{$segID_a};
		for my $idx_b ( $idx_a+1 .. $#segID ){
			my $segID_b = $segID[$idx_b];
			my $segOB_b = $Segment_Href->{$segID_b};
			# test overlap of two segments
			if( $segOB_a->test_overlap_seg(testSegOB=>$segOB_b) ){
				warn_and_exit "<ERROR>\tSegment $segID_a overlaps with segment $segID_b.\n";
			}
		}
	}

	# check segment string required for local-maps
	if( scalar(@{$V_Href->{SegStrLMhas}}) != 0 ){
		my $segstrHref = {};
		my $segstrNO = 0;
		for my $segstr (@{$V_Href->{SegStrLMhas}}){
			$segstrNO ++;
			my (@regex_fw, @regex_rv);
			# check valid and load regex_fw
			for my $seg (split /;/, $segstr){
				my ($segID, $orit) = ($seg =~ /^(\D\d+)([\+\-]?)$/);
				if(    !defined $segID
					|| !exists $Segment_Href->{$segID}
				){
					warn_and_exit "<ERROR>\tcannot recognize segID from $seg ('-segstr' parm).\n";
				}
				push @regex_fw, $segID . ( $orit ? "\\$orit" : '[\+\-]' );
			}
			# regex_rv
			@regex_rv = reverse @regex_fw;
			tr/+-/-+/ for @regex_rv;
			# record
			$segstrHref->{$segstrNO} = { segstr => $segstr,
										 regex_fw => join(';', @regex_fw),
										 regex_rv => join(';', @regex_rv)
										};
		}
		# all valid, load up!
		$V_Href->{SegStrLMhas} = $segstrHref;
	}
	else{
		$V_Href->{SegStrLMhas} = undef;
	}

	# inform
	stout_and_sterr "[INFO]\tCheck segments OK.\n";
}

#--- read pwlink and lrdist file ---
sub load_long_range_info{

	# read pwanchor-link file
	&Load_PWanchorLinkCount;

	# read long-range DNA length distribution
	&load_longRangeLenDist;
}

#--- read pwanchor link file ---
sub load_longRangeLenDist{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return unless( defined $V_Href->{LongRangeLenDistFile} );

	my %LRlen2count;
	my %LRlen2ratio;
	open (LRLD, Try_GZ_Read($V_Href->{LongRangeLenDistFile})) || die "fail read long-range length dist: $!\n";
	while(<LRLD>){
		next if(/^#/); # comments
		my ($LRlen, $counts, $ratio2C) = (split)[0,1,2];
		$LRlen2count{$LRlen} = $counts;
		$ratio2C ||= '1(1)';
		for my $r2C (split /,/, $ratio2C){
			my ($ratio, $count) = ($r2C =~ /^([\d\.]+)\((\d+)\)$/);
			$ratio = min(1, sprintf "%.2f", $ratio);
			$LRlen2ratio{$LRlen}{$ratio} += $count;
			# check
			if( $ratio > 1 ){
				warn_and_exit "<ERROR>\tencounter conver ratio larger than 1.0 in LongRangeLenDistFile.\n"
									."\t$V_Href->{LongRangeLenDistFile}\n";
			}
		}
	}
	close LRLD;

	# get prob distribution
	my $all_count = sum( values %LRlen2count );
	my $accu_count = 0;
	my %LRlen2accuProb = map{ $accu_count += $LRlen2count{$_}; ( $_, sprintf "%.6f", $accu_count/$all_count ) } sort {$a<=>$b}  keys %LRlen2count;

	# get ratio count distribution
	my %LRlen2ratioA;
	for my $LRlen (sort {$a<=>$b}  keys %LRlen2ratio){
		my $sumC = 0;
		$sumC += $LRlen2ratio{$LRlen}{$_} for keys %{$LRlen2ratio{$LRlen}};
		my $normSumC = min($sumC, 100);
		for my $ratio (sort {$a<=>$b}  keys %{$LRlen2ratio{$LRlen}}){
			my $normC = int($LRlen2ratio{$LRlen}{$ratio} * $normSumC / $sumC);
			push @{$LRlen2ratioA{$LRlen}}, $ratio for (1 .. $normC);
		}
	}

	# assign lRlen.dist to all groups
	for my $groupID (sort keys %{$V_Href->{Group}}){
		my $group_OB = $V_Href->{Group}->{$groupID};
		$group_OB->load_LRlenProbDist( LRlenToAccuProb_Href => \%LRlen2accuProb, LRlen2ratioA_Href => \%LRlen2ratioA );
	}

	# inform
	stout_and_sterr "[FILE]\tRead long-range length distribution OK.\n";
}

#--- read pwanchor link file ---
sub Load_PWanchorLinkCount{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return unless( defined $V_Href->{pairWiseLinkFile} );

	open (PWLC, Try_GZ_Read($V_Href->{pairWiseLinkFile})) || die "fail read pair-wise link counts: $!\n";
	while(<PWLC>){
		next if(/^#/); # comments
		my (%anchor, %type, %obj, $linkCount);
		($anchor{a}, $anchor{b}, $linkCount) = (split)[0,1,2];
		for my $i (qw/ a b /){
			if( $anchor{$i} =~ /^JUNC\d+$/ ){ # JUNC
				my $juncID = $anchor{$i};
				if( !exists $V_Href->{Junc}->{$juncID} ){
					warn_and_exit "<ERROR>\tCannot find junc_OB for junc_ID ($juncID) of anchor_$i.\n";
				}
				$type{$i} = 'JUNC';
				$obj{$i} = $V_Href->{Junc}->{$juncID};
			}
			elsif( $anchor{$i} =~ /^([HV]\d+),(.+):(\d+)$/ ){
				my ($segID, $chr, $pos) = ($1, $2, $3);
				if( !exists $V_Href->{Segment}->{$segID} ){
					warn_and_exit "<ERROR>\tCannot find seg_OB for seg_ID ($segID) of anchor_$i.\n";
				}
				# create pos_OB if possible
				my $posID = "$chr:$pos";
				if( !exists $V_Href->{Pos}->{$posID} ){
					my @parm = ( chr => $chr, pos => $pos, seg_OB => $V_Href->{Segment}->{$segID} );
					$V_Href->{Pos}->{$posID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::Pos_OB->new( @parm );
				}
				$type{$i} = 'POS';
				$obj{$i} = $V_Href->{Pos}->{$posID};
			}
			else{
				warn_and_exit "<ERROR>\tCannot analyze anchor_$i ($anchor{$i}).\n";
			}
		}
		# get new PWanchor_OB
		my $PWanchorID = join(';', sort ($anchor{a}, $anchor{b}));
		if( exists $V_Href->{PWanchor}->{$PWanchorID} ){
			warn_and_exit "<ERROR>\tPair-wise anchors ($PWanchorID) is duplicated.\n";
		}
		my @parm = ( ID => $PWanchorID, linkCount => $linkCount,
					 info_a => $anchor{a}, type_a => $type{a}, obj_a => $obj{a},
					 info_b => $anchor{b}, type_b => $type{b}, obj_b => $obj{b}  );
		$V_Href->{PWanchor}->{$PWanchorID} = FuseSV::Virus_Integ::LocalHaplotypeDetect::PWanchor_OB->new( @parm );
		# link PWanchor to Group
		my $groupOB = $V_Href->{PWanchor}->{$PWanchorID}->get_groupOB;
		$groupOB->load_PWanchorOB( pwanchorOB => $V_Href->{PWanchor}->{$PWanchorID} );
	}
	close PWLC;

	# inform
	stout_and_sterr "[FILE]\tRead pair-wise link file OK.\n";
}

#--- merge unit-cycle to local map ---
sub resolve_localmap{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	for my $groupID (sort keys %{$V_Href->{Group}}){
		my $group_OB = $V_Href->{Group}->{$groupID};
		$group_OB->resolve_localmap_of_solution( randTime => $V_Href->{randTime},
												 combMaxCnt => $V_Href->{probCombCount},
												 segstrHref => $V_Href->{SegStrLMhas},
												 useFVGM => $V_Href->{useFVGM},
												 seqType => $V_Href->{seqType},
												 aimDepth => $V_Href->{aimDepth},
												 readLen => $V_Href->{readLen},
												 forkNum => $V_Href->{forkNum},
												 corrMethod => $V_Href->{CorrMethod},
												 tmpPref => $V_Href->{output}
												);
	}
}

#--- display local map (merged UnitCycle) ---
## and pair-wise anchor report if set
sub generate_output{
	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	my $pwanchor_report_bool = 0;
	my $pwanchor_report_info = '';

	my @group_OB = map  { $_->{group_OB} }
				   sort { $a->{order_NO} <=> $b->{order_NO} }
				   grep $_->{group_OB}->is_valid, # only valid group
				   values %{$V_Href->{GroupOrder}};

	# local map output
	if( scalar(@group_OB) != 0 ){
		open (LOCALMAP, Try_GZ_Write($V_Href->{output})) || die "fail write local map file: $!\n";
		# config header, just copy
		print LOCALMAP "$_\n" for @{$V_Href->{ConfigHeader}};
		# each group
		for my $group_OB ( @group_OB ){
			# localmap
			print LOCALMAP $group_OB->display( sayID => 'origID' );
			# pair-wise anchors and pearson corr
			if( $group_OB->do_LM_Comb_LR_Corr ){
				$pwanchor_report_bool = 1;
				$pwanchor_report_info .= $group_OB->pwanchor_report;
			}
		}
		close LOCALMAP;
		# inform
		stout_and_sterr `date`."[FILE]\tGenerate local map file OK.\n".
									 "\t$V_Href->{output}\n";
	}
	else{
		# inform
		stout_and_sterr `date`."<WARN>\tNo valid group to output.\n";
	}

	# report of pair-wise anchor
	if( $pwanchor_report_bool ){
		my $pwanchor_report = $V_Href->{output}.'.pwanchor.report';
		open (PWAR, Try_GZ_Write($pwanchor_report)) || die "fail write pwanchor report: $!\n";
		print PWAR $pwanchor_report_info;
		close PWAR;
		# inform
		stout_and_sterr "[FILE]\tGenerate pair-wise anchor report OK.\n".
							  "\t$pwanchor_report\n";
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
