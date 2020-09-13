package FuseSV::Virus_Integ::GeneratePrimers;

use strict;
use warnings;
use File::Spec qw/catfile/;
use Getopt::Long;
use SOAPfuse::General_Operation qw/write_fasta_file warn_and_exit stout_and_sterr trible_run_for_success/;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  real_generate_primers
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'GeneratePrimers';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-01-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						real_generate_primers
						Prepare_Global_Tags
						Generate_Primer_For_FuseSV
						Generate_Primer_For_SOAPfuse_Fusion
						Get_Primer
						Extract_opt_primer
						FuseSV_Primer_Output
						SOAPfuse_Fusion_Primer_Output
						Get_Exclude_Region_From_Fseq
						Extract_Seq_from_FastaFile
						Get_Seq_from_Fasta
						Generate_Primer_For_General_Mode
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} primers <[Options]>

	 Options:

	   # Inputs and Outputs #
	    -odir     [s]   Directory to store the primers results. <required>
	    -prefix   [s]   Id will be used as the prefix of result files in outdir. <required>

	   # Software Required #
	    -primer3  [s]   Filepath of primer3_score bin. <required>

	   # For prepared (RT-)PCR files #
	    -SOAPfuse [s]   The Seq-for-RT-PCR file created by SOAPfuse for your Fusion analysis.
	    -FuseSV   [s]   The Seq-for-PCR file created by FuseSV for your FuseSV analysis.
	    -mindist  [i]   Minimum bases between the primer-3-end and the BreakPoint. [50]

	   # General primer construction #
	    -sor_fa   [s]   Fasta format file to be the sequence soruce for you PCR target.
	    -region   [s]   File contains the region you want to PCR out, the region is from the '-sor_fa' file.
	                     TAB delimited format:   Seq-ID \\t RefSeg-ID \\t start_pos \\t end_pos
	    -target   [s]   Interval list of target region. ['200,200']
	                     This is refer to the '-sor_fa' OR selected regions ('-region') from '-sor_fa' file.
	                     Format: <start>,<length>
	                     <start> is the index of the first base of a Target, and <length> is its length.

	   # Gereral options for primer3 #
	   ## product of PCR
	    -product  [s]   Size range of PCR product. ['100-300']
	   ## melting temperature (TM)
	    -opttm    [f]   Optimum TM(Celsius) for a primer oligo. [57]
	    -mintm    [f]   Minimum acceptable TM(Celsius) for a primer oligo. [53]
	    -maxtm    [f]   Maximum acceptable TM(Celsius) for a primer oligo. [61]
	    -maxdiff  [f]   Maximum acceptable (unsigned) difference between the TMs of the left and right primers. [5.0]
	   ## GC content
	    -maxgc    [i]   Maximum allowable percentage of Gs and Cs in any primer. [60]
	    -mingc    [i]   Minimum allowable percentage of Gs and Cs in any primer. [40]
	   ## other options
	    -OtherTag [s]   Other global tags of primer3, this para can use multiple times, once one option.
	                     You can not use this para to reset the options aboved. You can input like this:
	                      ... -OtherTag 'PRIMER_MIN_SIZE=16' -OtherTag 'PRIMER_MAX_SIZE=25' ...
	                     a) http://www.molquest.com/help/2.3/programs/Primer3/parameters.html
	                     b) http://primer3.sourceforge.net/primer3_manual.htm

	    -h|help   Display this help info.

 	 Version:
	    $FuseSV::Virus_Integ::GeneratePrimers::VERSION at $FuseSV::Virus_Integ::GeneratePrimers::DATE

	 Author:
	    $FuseSV::Virus_Integ::GeneratePrimers::AUTHOR ($FuseSV::Virus_Integ::GeneratePrimers::EMAIL)
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
			[ primer_outdir => '' ],
			[ primer_out_prefix => '' ],
			# software
			[ primer3_bin => '' ],
			# for FuseSV/SOAPfuse
			[ Seq_file_For_PCR_FuseSV => '' ], # FuseSV output
			[ Seq_file_For_PCR_SOAPfuse => '' ], # SOAPfuse output
			[ primer_minDist_from_breakpoint => 50 ],
			# for general fasta
			[ For_primer_Source_fasta => '' ], # general fasta
			[ target_interval_for_primer => '200,200' ],
			[ Region_file_for_primer => '' ], # region to selection from genome file
			# options
			[ MIN_Region_For_Primer => 30 ],
			# primer3 options
			[ PRIMER_General_Opt => { 	PRIMER_PRODUCT_SIZE_RANGE => '100-300',
										PRIMER_OPT_TM => 57,
										PRIMER_MIN_TM => 53,
										PRIMER_MAX_TM => 61,
										PRIMER_MAX_DIFF_TM => 5.0,
										PRIMER_MAX_GC => 60,
										PRIMER_MIN_GC => 40,
										PRIMER_FIRST_BASE_INDEX => 1            } ],
			[ PRIMER_other_opt => [] ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['primer_outdir'],
								  ['primer3_bin'],
								  ['Seq_file_For_PCR_FuseSV'],
								  ['Seq_file_For_PCR_SOAPfuse'],
								  ['For_primer_Source_fasta'],
								  ['Region_file_for_primer']     ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-odir:s"		=> \$V_Href->{primer_outdir},
		"-prefix:s"		=> \$V_Href->{primer_out_prefix},
		# software
		"-primer3:s"	=> \$V_Href->{primer3_bin},
		# For SOAPfuse(fusion) or FuseSV(SV), which knows the breakpoint of (RT-)PCR
		"-FuseSV:s"		=> \$V_Href->{Seq_file_For_PCR_FuseSV},
		"-SOAPfuse:s"	=> \$V_Href->{Seq_file_For_PCR_SOAPfuse},
		"-mindist:i"	=> \$V_Href->{primer_minDist_from_breakpoint},
		# For General primer construction from Genome
		"-sor_fa:s"		=> \$V_Href->{For_primer_Source_fasta},
		"-target:s"		=> \$V_Href->{target_interval_for_primer}, # with aimfa
		"-region:s"		=> \$V_Href->{Region_file_for_primer}, # with genome
		# Gereral Para for Primer3
		## target
		"-product:s"	=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_PRODUCT_SIZE_RANGE},
		## melting temperature
		"-opttm:f"		=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_OPT_TM},
		"-mintm:f"		=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_MIN_TM},
		"-maxtm:f"		=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_MAX_TM},
		"-maxdiff:f"	=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_MAX_DIFF_TM},
		## GC content
		"-maxgc:i"		=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_MAX_GC},
		"-mingc:i"		=> \$V_Href->{PRIMER_General_Opt}->{PRIMER_MIN_GC},
		## other general options of primer3
		"-OtherTag:s"	=>  $V_Href->{PRIMER_other_opt},
		# help
		"-h|help"		=> \$V_Href->{HELP},
		# for debug
		"-debug"		=> \$V_Href->{in_debug} # hidden option
	);
}

#--- test para and alert ---
sub para_alert{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	return  (   $V_Href->{HELP}
			 || !$V_Href->{primer_outdir}
			 || !$V_Href->{primer_out_prefix}
			 || (!$V_Href->{primer3_bin} || !-e $V_Href->{primer3_bin})
			 || (    (!$V_Href->{Seq_file_For_PCR_FuseSV}   || !-e $V_Href->{Seq_file_For_PCR_FuseSV})
			 	  && (!$V_Href->{Seq_file_For_PCR_SOAPfuse} || !-e $V_Href->{Seq_file_For_PCR_SOAPfuse})
			 	  && (!$V_Href->{For_primer_Source_fasta}      || !-e $V_Href->{For_primer_Source_fasta})
			 	  && (!$V_Href->{Region_file_for_primer}    || !-e $V_Href->{Region_file_for_primer})
			 	)
			 || (-e $V_Href->{Region_file_for_primer} && !-e $V_Href->{For_primer_Source_fasta})
			);
}


#--- introduc mutated into virus genome
sub real_generate_primers{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# prepare the global tags
	&Prepare_Global_Tags;

	# prepare workspace
	## details dir
	$V_Href->{primer_detail_dir} = File::Spec->catfile($V_Href->{primer_outdir}, "$V_Href->{primer_out_prefix}-primer-details");
	`rm -rf $V_Href->{primer_detail_dir}`;
	`mkdir -p $V_Href->{primer_detail_dir}`;
	## output
	$V_Href->{primer_out} = File::Spec->catfile($V_Href->{primer_outdir}, "$V_Href->{primer_out_prefix}.primers.result");
	$V_Href->{primer_fail} = File::Spec->catfile($V_Href->{primer_outdir}, "$V_Href->{primer_out_prefix}.primers.fail");
	`rm -f $V_Href->{primer_out}`;
	`rm -f $V_Href->{primer_fail}`;

	#---------------------------------------------#
	#---------|      Main-function     |----------#
	#---------------------------------------------#

	if( $V_Href->{Seq_file_For_PCR_FuseSV} ){ # FuseSV
		&Generate_Primer_For_FuseSV;
	}
	elsif( $V_Href->{Seq_file_For_PCR_SOAPfuse} ){ # SOAPfuse
		&Generate_Primer_For_SOAPfuse_Fusion;
	}
	elsif( $V_Href->{For_primer_Source_fasta} ){ # General mode
		my %Seq;
		# prepare seq for primer, maybe selected region from source fasta
		&Extract_Seq_from_FastaFile( { Seq_Href => \%Seq } );
		# generate primers for each seq
		&Generate_Primer_For_General_Mode( { Seq_Href => \%Seq } );
	}

	stout_and_sterr `date`."[INFO]\tGenerate Primers of $V_Href->{primer_out_prefix} successfully.\n";
}

#--- check and prepare the global tags
sub Prepare_Global_Tags{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $General_Opt_Href = $V_Href->{PRIMER_General_Opt};
	my $Other_Opt_Aref = $V_Href->{PRIMER_other_opt};

	# deal the other general options
	for my $other_opt ( @$Other_Opt_Aref ){
		my ($tag, $value) = ($other_opt =~ /^(.+)=(.+)$/);
		if( exists($General_Opt_Href->{$tag}) && $General_Opt_Href->{$tag} ne '' ){
			warn_and_exit `date`."<Error>\tTag $tag has been set as $General_Opt_Href->{$tag} before.\n";
		}
		if($tag eq 'TARGET'){
			warn_and_exit `date`."<Error>\tPlease set Tag $tag via '-target'.\n"
								."\tFor SOAPfuse or FuseSV, ignore this para, as it will be set and adjusted automatically.\n";
		}
		$General_Opt_Href->{$tag} = $value;
	}
	# update
	$V_Href->{PRIMER_General_Opt_String} = join("\n", map $_.'='.$General_Opt_Href->{$_} , sort keys %$General_Opt_Href);
}

#--- this is for FuseSV primers
sub Generate_Primer_For_FuseSV{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $MIN_Region_For_Primer = $V_Href->{MIN_Region_For_Primer};
	my $primer_minDist_from_breakpoint = $V_Href->{primer_minDist_from_breakpoint};

	my %FuseSV_Seq;
	# read the Fuse-SV seq
	open (SEQ, $V_Href->{Seq_file_For_PCR_FuseSV}) || die "fail $V_Href->{Seq_file_For_PCR_FuseSV}: $!\n";
	while(<SEQ>){
		my ($sv_case, $sv_seq) = (split)[0,1];
		$FuseSV_Seq{$sv_case} = $sv_seq;
	}
	close SEQ;

	my %Fail;
	my %Primer;
	my $primer_time = 1;
	for my $sv_case (sort keys %FuseSV_Seq){
		my $sv_seq = $FuseSV_Seq{$sv_case};
		my ($up_SVseq) = ($sv_seq =~ /^([^\|]+)\|/); # for the case that have inner insertion
		my ($dw_SVseq) = ($sv_seq =~ /\|([^\|]+)$/);
		# check whether too short
		if(length($up_SVseq) < ($MIN_Region_For_Primer + $primer_minDist_from_breakpoint)){
			$Fail{$sv_case} = 'Upstream_SVseq_is_Too_Short';
			next;
		}
		if(length($dw_SVseq) < ($MIN_Region_For_Primer + $primer_minDist_from_breakpoint)){
			$Fail{$sv_case} = 'Downstream_SVseq_is_Too_Short';
			next;
		}
		# inner insertion bases
		my ($inner_ins) = ($sv_case =~ /\/(\w+)$/);
		my $inner_ins_len = ($inner_ins eq 'NONE') ? 0 : length($inner_ins);
		# target region
		my $target_start = length($up_SVseq) - $primer_minDist_from_breakpoint;
		my $target_length = $inner_ins_len + $primer_minDist_from_breakpoint * 2;
		my $target_para = "TARGET=$target_start,$target_length";
		# get the exclude region
		(my $SVseq = $sv_seq) =~ s/\|//g;
		my $Exclude_Region_Para = 'EXCLUDED_REGION='.&Get_Exclude_Region_From_Fseq( { Seq => $SVseq } );
		# generate primers
		my ($opt_primer,$opt_info) = &Get_Primer(
											{
												primer_time => $primer_time,
												Seq_ID => $sv_case,
												Seq => $SVseq,
												target_para => $target_para,
												Exclude_Region_Para => $Exclude_Region_Para
											}
									 );
		# if null, no exclude region again
		if($opt_primer eq 'null'){
			($opt_primer,$opt_info) = &Get_Primer(
											{
												primer_time => $primer_time,
												Seq_ID => $sv_case,
												Seq => $SVseq,
												target_para => $target_para
											}
									  );
		}
		# still is null
		if($opt_primer eq 'null'){
			$Fail{$sv_case} = "NO.$primer_time-No_available_Primer";
		}
		else{
			my ($left_TM, $right_TM, $left_pos_len, $right_pos_len, $product_size) = (split /;/,$opt_info)[0..4];
			my ($left_pos, $left_len) = ($left_pos_len =~ /^(\d+),(\d+)$/);
			my ($right_pos, $right_len) = ($right_pos_len =~ /^(\d+),(\d+)$/);
			# store the primer for next output
			my $primer_id = "$opt_primer;$left_TM;$right_TM;$left_len;$right_len";
			push @{$Primer{$primer_id}->{$product_size}} , [ $primer_time, "$left_pos,$right_pos", $sv_case ];
		}
		$primer_time++;
	}

	# Output the primer result and fails
	&FuseSV_Primer_Output(
		{
			Primer_Href => \%Primer,
			Fail_Href => \%Fail
		}
	);
}

#--- this is for SOAPfuse fusion primers
sub Generate_Primer_For_SOAPfuse_Fusion{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $MIN_Region_For_Primer = $V_Href->{MIN_Region_For_Primer};
	my $primer_minDist_from_breakpoint = $V_Href->{primer_minDist_from_breakpoint};

	my %Fusion_Seq;
	# read the fusion seq
	open (SEQ, $V_Href->{Seq_file_For_PCR_SOAPfuse}) || die "fail $V_Href->{Seq_file_For_PCR_SOAPfuse}: $!\n";
	while(<SEQ>){
		my ($gene_case,$trans_case,$fusion_seq) = (split)[0,1,2];
		push @{$Fusion_Seq{$gene_case}{$fusion_seq}} , $trans_case;
	}
	close SEQ;

	my %Fail;
	my %Primer;
	my $primer_time = 1;
	for my $gene_case (sort keys %Fusion_Seq){
		for my $fusion_seq (sort keys %{$Fusion_Seq{$gene_case}}){
			my ($up_fseq, $dw_fseq) = (split /\|/,$fusion_seq)[0,1];
			# SEQ ID
			my $trans_case_str = join(';',sort @{$Fusion_Seq{$gene_case}{$fusion_seq}});
			my $Seq_ID = $gene_case.':'.$trans_case_str;
			# check whether too short
			if(length($up_fseq) < ($MIN_Region_For_Primer + $primer_minDist_from_breakpoint)){
				$Fail{$Seq_ID} = 'Upstream_Fseq_is_Too_Short';
				next;
			}
			if(length($dw_fseq) < ($MIN_Region_For_Primer + $primer_minDist_from_breakpoint)){
				$Fail{$Seq_ID} = 'Downstream_Fseq_is_Too_Short';
				next;
			}
			# target region
			my $target_start = length($up_fseq) - $primer_minDist_from_breakpoint;
			my $target_length = $primer_minDist_from_breakpoint * 2;
			my $target_para = "TARGET=$target_start,$target_length";
			# get the exclude region
			my $Seq = $up_fseq.$dw_fseq;
			my $Exclude_Region_Para = 'EXCLUDED_REGION=' . &Get_Exclude_Region_From_Fseq( { Seq => $Seq } );
			# generate primers
			my ($opt_primer, $opt_info) = &Get_Primer(
											{
												primer_time => $primer_time,
												Seq_ID => $Seq_ID,
												Seq => $Seq,
												target_para => $target_para,
												Exclude_Region_Para => $Exclude_Region_Para
											}
										  );
			# if null, no exclude region again
			if($opt_primer eq 'null'){
				($opt_primer, $opt_info) = &Get_Primer(
											{
												primer_time => $primer_time,
												Seq_ID => $Seq_ID,
												Seq => $Seq,
												target_para => $target_para
											}
										   );
			}
			# still is null
			if($opt_primer eq 'null'){
				$Fail{$Seq_ID} = "NO.$primer_time-No_available_Primer";
			}
			else{
				my ($left_TM, $right_TM, $left_pos_len, $right_pos_len, $product_size) = (split /;/,$opt_info)[0..4];
				my ($left_pos, $left_len) = ($left_pos_len =~ /^(\d+),(\d+)$/);
				my ($right_pos, $right_len) = ($right_pos_len =~ /^(\d+),(\d+)$/);
				# store the primer for next output
				my $primer_id = "$opt_primer;$left_TM;$right_TM;$left_len;$right_len";
				push @{$Primer{$primer_id}->{$product_size}->{$gene_case}} , [ "$left_pos,$right_pos", $trans_case_str, $primer_time ];
			}
			$primer_time++;
		}
	}

	# Output the primer result and fails
	&SOAPfuse_Fusion_Primer_Output(
		{
			Primer_Href => \%Primer,
			Fail_Href => \%Fail
		}
	);
}

#--- generate primers
sub Get_Primer{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $primer_time = $Option_Href->{primer_time};
	my $Seq_ID = $Option_Href->{Seq_ID};
	my $Seq = $Option_Href->{Seq};
	my $target_para = $Option_Href->{target_para};
	my $Exclude_Region_Para = $Option_Href->{Exclude_Region_Para} || '';

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $primer_detail_dir = $V_Href->{primer_detail_dir};
	my $PRIMER_General_Opt_String = $V_Href->{PRIMER_General_Opt_String};
	my $primer3_bin = $V_Href->{primer3_bin};

	# prepare the input for primer3
	my $input_for_primers = File::Spec->catfile($primer_detail_dir, "Case$primer_time.input.for.primer3");
	open (INP,">$input_for_primers") || die "fail $input_for_primers: $!\n";
	print INP "PRIMER_SEQUENCE_ID=$Seq_ID\n";
	print INP "SEQUENCE=$Seq\n";
	print INP "$target_para\n";
	print INP "$Exclude_Region_Para\n" if($Exclude_Region_Para);
	print INP "$PRIMER_General_Opt_String\n";
	print INP "=\n";
	close INP;

	# primer3_core
	my $format_out = File::Spec->catfile($primer_detail_dir, "Case$primer_time.primer3.format.out.txt");
	my $format_command = "$primer3_bin -format_output <$input_for_primers >$format_out";
	trible_run_for_success($format_command, 'format_output_primer', {cmd_Nvb=>1,esdo_Nvb=>1});
	my $strict_out = File::Spec->catfile($primer_detail_dir, "Case$primer_time.primer3.strict.out.txt");
	my $strict_command = "$primer3_bin -strict_tags <$input_for_primers >$strict_out";
	trible_run_for_success($strict_command, 'strict_output_primer', {cmd_Nvb=>1,esdo_Nvb=>1});
	# opt primer?
	return &Extract_opt_primer( { strict_out => $strict_out } );
}

#--- extract the optimized primer from primer3 strict out
sub Extract_opt_primer{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $strict_out = $Option_Href->{strict_out};

	my %INFO;
	open (ST,$strict_out) || die "fail $strict_out: $!\n";
	while(<ST>){
		chomp;
		my ($ele,$value) = (split /\=/)[0,1];
		if(/^(\S+)=(\S+)$/){
			$INFO{$1} = $2;
		}
	}
	close ST;

	# return
	if(exists($INFO{PRIMER_LEFT_SEQUENCE})){ # have primer
		my $opt_primer = join(';',@INFO{('PRIMER_LEFT_SEQUENCE','PRIMER_RIGHT_SEQUENCE')});
		my $opt_info = join(';',@INFO{('PRIMER_LEFT_TM','PRIMER_RIGHT_TM','PRIMER_LEFT','PRIMER_RIGHT','PRIMER_PRODUCT_SIZE')});
		return ( $opt_primer, $opt_info );
	}
	else{
		return ( 'null', 'null' );
	}
}

#--- Output the primer result of SOAPfuse fusion
##-- also prepare for general mode
sub FuseSV_Primer_Output{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Primer_Href = $Option_Href->{Primer_Href};
	my $Fail_Href = $Option_Href->{Fail_Href};
	my $General_mode = $Option_Href->{General_mode} || 0;

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $primer_out = $V_Href->{primer_out};
	my $primer_fail = $V_Href->{primer_fail};

	# output the primers
	open (POUT,">$primer_out") || die "fail $primer_out: $!\n";
	my @theme = ( '#left_primer',
				  'right_primer',
				  'left_primer_length',
				  'right_primer_length',
				  'left_TM',
				  'right_TM',
				  'product_size',
				  ($General_mode) ? 'involved_aim_Seq-ID_list' : 'involved_sv_case_list',
				  # (primers_detail_NO:left_primer_position,right_primer_position:aim_Seq-ID)
				  # (primers_detail_NO:left_primer_position,right_primer_position:sv_case)
				  ($General_mode) ? 'detailed_info_of_aim_Seq' : 'detailed_info_of_sv_case'
				);
	print POUT join("\t",@theme) . "\n";
	for my $primer_info (sort keys %$Primer_Href){
		my ($left_primer,$right_primer,$left_TM,$right_TM,$left_len,$right_len) = (split /;/,$primer_info);
		for my $product_size (sort keys %{$Primer_Href->{$primer_info}}){
			my @sv_info;
			my %sv_case;
			for my $Aref (@{$Primer_Href->{$primer_info}->{$product_size}}){
				my ($primer_time,$primer_pos,$sv_case) = @$Aref;
				push @sv_info , "NO.$primer_time:$primer_pos:$sv_case";
				$sv_case{$sv_case} = 1;
			}
			my @output_line = ( $left_primer,
								$right_primer,
								$left_len,
								$right_len,
								$left_TM,
								$right_TM,
								$product_size,
								join(';',sort keys %sv_case),
								join('|',@sv_info)
							  );
			print POUT join("\t", @output_line) . "\n";
		}
	}
	close POUT;

	# output the fails
	return if(scalar(keys %$Fail_Href) == 0); # none
	open (FAIL,">$primer_fail") || die "fail $primer_fail: $!\n";
	print FAIL '#Cases in this file may need manual primers construction'."\n";
	print FAIL join("\t",(($General_mode)?'#aim_Seq-ID':'#sv_case','fail_reason'))."\n";
	for my $sv_case (sort keys %$Fail_Href){
		my $reason = $Fail_Href->{$sv_case};
		print FAIL join("\t", ($sv_case, $reason) )."\n";
	}
	close FAIL;
}

#--- Output the primer result of SOAPfuse fusion
sub SOAPfuse_Fusion_Primer_Output{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Primer_Href = $Option_Href->{Primer_Href};
	my $Fail_Href = $Option_Href->{Fail_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $primer_out = $V_Href->{primer_out};
	my $primer_fail = $V_Href->{primer_fail};

	# output the primers
	open (POUT,">$primer_out") || die "fail $primer_out: $!\n";
	my @theme = ( '#left_primer',
				  'right_primer',
				  'left_primer_length',
				  'right_primer_length',
				  'left_TM',
				  'right_TM',
				  'product_size',
				  'involved_fusion_gene_case_list',
				  'detailed_fusion_case_primer_info' # （primers_detail_NO:left_primer_position,right_primer_position:[fusion_gene_case:fusion_trans_case_list]）
				  );
	print POUT join("\t", @theme). "\n";
	for my $primer_info (sort keys %$Primer_Href){
		my ($left_primer, $right_primer, $left_TM, $right_TM, $left_len, $right_len) = (split /;/,$primer_info);
		for my $product_size (sort keys %{$Primer_Href->{$primer_info}}){
			my $fusion_case = '';
			for my $gene_case (sort keys %{$Primer_Href->{$primer_info}->{$product_size}}){
				my @info;
				for my $Aref (@{$Primer_Href->{$primer_info}->{$product_size}->{$gene_case}}){
					my ($primer_pos, $trans_case_str, $primer_time) = @$Aref;
					push @info , "NO.$primer_time:$primer_pos:[$gene_case:$trans_case_str]";
				}
				$fusion_case .= '(' . join('|',@info) . ')';
			}
			my $gene_case_list = join(';',sort keys %{$Primer_Href->{$primer_info}->{$product_size}});
			my @output_line = ( $left_primer,
								$right_primer,
								$left_len,
								$right_len,
								$left_TM,
								$right_TM,
								$product_size,
								$gene_case_list,
								$fusion_case
							  );
			print POUT join("\t", @output_line) . "\n";
		}
	}
	close POUT;

	# output the fails
	return if(scalar(keys %$Fail_Href) == 0); # none
	open (FAIL,">$primer_fail") || die "fail $primer_fail: $!\n";
	print FAIL '#Cases in this file may need manual primers construction'."\n";
	print FAIL join("\t",('#fusion_gene_case','fail_reason','fusion_trans_case'))."\n";
	for my $Seq_ID (sort keys %$Fail_Href){
		my ($gene_case, $trans_case) = (split /\:/,$Seq_ID)[0,1];
		my $reason = $Fail_Href->{$Seq_ID};
		print FAIL join("\t", ($gene_case, $reason, $trans_case) ) . "\n";
	}
	close FAIL;
}

#--- get the lowerbase-or-N region to exclude primers
sub Get_Exclude_Region_From_Fseq{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $fseq = $Option_Href->{Seq};

	my @fseq = split //,$fseq;
	# get the exclude region end loc
	my @exclude_region = ();
	for my $i (0 .. $#fseq){
		# first base
		if($i == 0 && $fseq[$i] =~ /[acgtnN]/){
			push @exclude_region , $i+1;
		}
		# last base
		elsif($i == $#fseq && $fseq[$i] =~ /[acgtnN]/){
			push @exclude_region , $i+1;
			# if last one is just upper -> lower|N, store again
			if($fseq[$i-1] =~ /[ACGT]/ && $fseq[$i] =~ /[acgtnN]/){
				push @exclude_region , $i+1;
			}
		}
		# inner base
		elsif($i != 0){ # not the first base
			if($fseq[$i-1] =~ /[acgtnN]/ && $fseq[$i] =~ /[ACGT]/){ # lower|N -> upper
				push @exclude_region , $i;
			}
			elsif($fseq[$i-1] =~ /[ACGT]/ && $fseq[$i] =~ /[acgtnN]/){ # upper -> lower|N
				push @exclude_region , $i+1;
			}
		}
	}
	warn_and_exit `date`."<ERROR>\tWrong exclude region for fseq:\n$fseq\n@exclude_region\n" if( scalar(@exclude_region) % 2 != 0);
	# to interval
	my @exclude_para = map $exclude_region[2*$_].','.($exclude_region[2*$_+1]-$exclude_region[2*$_]+1) , (0 .. scalar(@exclude_region)/2-1);
	return join(' ', @exclude_para);
}

#--- extract seq from genome based on region file
sub Extract_Seq_from_FastaFile{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Seq_Href = $Option_Href->{Seq_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $region_file = $V_Href->{Region_file_for_primer};
	my $For_primer_Source_fasta = $V_Href->{For_primer_Source_fasta};

	# read region file
	my %region;
	if( -e $region_file ){
		open (RG, $region_file) || die "fail $region_file: $!\n";
		while(<RG>){
			my ($Seq_name, $chr, $st, $ed) = (split)[0..3];
			$region{$chr}->{$Seq_name} = [ $st, $ed ];
		}
		close RG;
	}

	# read selected region from genome file
	&Get_Seq_from_Fasta(
		{
			Seq_Href => $Seq_Href,
			source_fa => $For_primer_Source_fasta,
			region_Href => ((-e $region_file) ? \%region : 0)
		}
	);
}

#--- read seq from fasta file
sub Get_Seq_from_Fasta{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Seq_Href = $Option_Href->{Seq_Href};
	my $source_fa = $Option_Href->{source_fa};
	my $region_Href = $Option_Href->{region_Href} || 0;

	# read fasta file supplied by user
	open (FA, $source_fa) || die "fail $source_fa: $!\n";
	$/=">"; <FA>; $/="\n";
	while(<FA>){
		chomp(my $Seq_name = $_); ## remove the last "\n"
		$/=">";
		chomp(my $Seq = <FA>); ## remove the last '>'
		$/="\n";
		if( $region_Href ){
			next unless( exists($region_Href->{$Seq_name}) );
			$Seq =~ s/\s+//g;
			foreach my $SubSeq_name (keys %{$region_Href->{$Seq_name}}) {
				my ($st, $ed) = @{$region_Href->{$Seq_name}->{$SubSeq_name}};
				my $SubSeq = substr($Seq, $st-1, $ed-$st+1);
				$Seq_Href->{$SubSeq_name} = uc($SubSeq);
			}
		}
		else{
			$Seq =~ s/\s+//g;
			$Seq_Href->{$Seq_name} = uc($Seq);
		}
	}
	close FA;
}

#--- generate the primers
sub Generate_Primer_For_General_Mode{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Seq_Href = $Option_Href->{Seq_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $target_region = $V_Href->{target_interval_for_primer};

	my $target_para = 'TARGET='.$target_region;
	my %Fail;
	my %Primer;
	my $primer_time = 1;
	for my $Seq_name (sort keys %$Seq_Href){
		my $Seq = $Seq_Href->{$Seq_name};
		my ($opt_primer, $opt_info) = &Get_Primer(
											{
												primer_time => $primer_time,
												Seq_ID => $Seq_name,
												Seq => $Seq,
												target_para => $target_para
											}
									  );
		if($opt_primer eq 'null'){
			$Fail{$Seq_name} = "NO.$primer_time-No_available_Primer";
		}
		else{
			my ($left_TM, $right_TM, $left_pos_len, $right_pos_len, $product_size) = (split /;/,$opt_info)[0..4];
			my ($left_pos, $left_len) = ($left_pos_len =~ /^(\d+),(\d+)$/);
			my ($right_pos, $right_len) = ($right_pos_len =~ /^(\d+),(\d+)$/);
			# store the primer for next output
			my $primer_id = "$opt_primer;$left_TM;$right_TM;$left_len;$right_len";
			push @{$Primer{$primer_id}->{$product_size}} , [ $primer_time, "$left_pos,$right_pos", $Seq_name ];
		}
		$primer_time++;
	}

	# Output the primer result and fails
	&FuseSV_Primer_Output(
		{
			Primer_Href => \%Primer,
			Fail_Href => \%Fail,
			General_mode => 1
		}
	);
}

#--- 
1; ## tell the perl script the successful access of this module.
