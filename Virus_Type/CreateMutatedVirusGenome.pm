package FuseSV::Virus_Type::CreateMutatedVirusGenome;

use strict;
use warnings;
use File::Basename qw/basename dirname/;
use List::Util qw[min max sum];
use Getopt::Long;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::BioInfo::FASTA qw/ write_fasta_file /;
use FuseSV::LoadOn;
use FuseSV::GetPath;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  create_mutated_virus_genome
			  rand_mutation
			  rand_snv
			  rand_ins
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'CreateMutatedVirusGenome';
#----- version --------
$VERSION = "0.09";
$DATE = '2019-05-02';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						return_HELP_INFO
						Load_moduleVar_to_pubVarPool
						Get_Cmd_Options
						para_alert
						create_mutated_virus_genome
						output_mutated_virus
						load_virus_info
						load_virus_genome
						mutate_virus
						rand_mutation
						rand_snv
						rand_ins
					 /;

#--- return HELP_INFO
sub return_HELP_INFO{
 my $V_Href = FuseSV::LoadOn->load_variants_dict();
 return "
	 Usage:   perl $V_Href->{MainName} mut_virus <[Options]>

	 Options:

	   # Inputs and Outputs #
	    -odir   [s]  directory to store outputs, must exist already. <required>
	    -vpf    [s]  virus package folder. <required>
	                  You can make it via function 'make_virusdb'.
	                  OR Download these widely used from FuseSV official website, e.g., HPV.

	   # Mutation Option #
	    -virus  [s]  Select certain virus. [all]
	    -SNVf   [f]  SNV frequency. [0.01]
	    -INDELf [f]  INDEL frequency. [0.001]
	    -MNPf   [f]  MNP frequency. [0.001]
	    -mxIDL  [i]  Maximum length of random InDel. [15]
	    -mxMNPL [i]  Maximum length of random MNP. [3]
	    -idlocd [i]  Require certain local different bases to avoid shift of InDel. [1]

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
			[ mut_virus_outdir => '' ],
			# mutate options
			[ mut_select_virus => [] ], # array Aref
			[ mut_snv_freq => 0.01 ],
			[ mut_indel_freq => 0.001 ],
			[ mut_MNP_freq => 0.001 ],
			[ mut_max_InDel_Length => 15 ],
			[ mut_max_MNP_Length => 3 ],
			[ mut_indel_localDiff => 1 ],

			# list to abs-path
			[ ToAbsPath_Aref => [ ['virus_db_folder'],
								  ['mut_virus_outdir']      ] ]
		);
}

#--- get options from command line
sub Get_Cmd_Options{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# get options
	GetOptions(
		# input/output
		"-vpf:s"	=> \$V_Href->{virus_db_folder},
		"-odir:s"	=> \$V_Href->{mut_virus_outdir},
		# mutate options
		"-virus=s"	=> \@{$V_Href->{mut_select_virus}},
		"-SNVf=f"	=> \$V_Href->{mut_snv_freq},
		"-INDELf=f"	=> \$V_Href->{mut_indel_freq},
		"-MNPf=f"	=> \$V_Href->{mut_MNP_freq},
		"-mxIDL=i"	=> \$V_Href->{mut_max_InDel_Length},
		"-mxMNPL=i"	=> \$V_Href->{mut_max_MNP_Length},
		"-idlocd=i"	=> \$V_Href->{mut_indel_localDiff},
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
			 || (!$V_Href->{virus_db_folder} || !-d $V_Href->{virus_db_folder})
			 || (!$V_Href->{mut_virus_outdir} || !-d $V_Href->{mut_virus_outdir} || !-w $V_Href->{mut_virus_outdir})
			);
}

#--- introduc mutated into virus genome
sub create_mutated_virus_genome{

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my @select_virus = @{$V_Href->{mut_select_virus}};

	#-----------------#
	# load virus info #
	#-----------------#
	my %select_virus = map{($_,1)} @select_virus;
	my %virusInfo;
	&load_virus_info(
		virusInfo_Href => \%virusInfo,
		virusInfo_file => $V_Href->{virus_info_doc},
		needvirus_Href => (scalar(@select_virus)==0)?undef:\%select_virus
	);

	#---------------------#
	# load each virus seq #
	#---------------------#
	&load_virus_genome(
		virusInfo_Href => \%virusInfo,
		virus_orig_fa => $V_Href->{virus_orig_genome_fa}
	);

	#--------------------#
	# generate mutations #
	#--------------------#
	&mutate_virus(
		virusInfo_Href => \%virusInfo,
		snv_freq => $V_Href->{mut_snv_freq},
		indel_freq => $V_Href->{mut_indel_freq},
		MNP_freq => $V_Href->{mut_MNP_freq},
		indel_localDiff => $V_Href->{mut_indel_localDiff},
		max_InDel_Length => $V_Href->{mut_max_InDel_Length},
		max_MNP_Length => $V_Href->{mut_max_MNP_Length}
	);

	#----------------------#
	# output mutated virus #
	# and mutations on it  #
	#----------------------#
	&output_mutated_virus(
		virusInfo_Href => \%virusInfo
	);
}

#--- output mutated virus genome and list
sub output_mutated_virus{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $virusInfo_Href = $parm{virusInfo_Href};

	for my $vname (keys %$virusInfo_Href){
		my $vlength = $virusInfo_Href->{$vname}->{length};
		my $mut_Href = $virusInfo_Href->{$vname}->{mutations};
		my @ref_base = split //, $virusInfo_Href->{$vname}->{orig_seq};
		my $mut_vseq = '';
		# mutations list
		my $mutation_list = FuseSV::GetPath->GetPath( filekey => 'simu_virus_mutation_list', virus_subtype => $vname );
		open (MLT, Try_GZ_Write($mutation_list)) || die "fail write $mutation_list: $!\n";
		for my $pos (1 .. $vlength){
			if(!exists($mut_Href->{$pos})){
				$mut_vseq .= $ref_base[$pos-1];
			}
			else{
				my $mutType = $mut_Href->{$pos}->{type};
				my $mutBase = $mut_Href->{$pos}->{mutbase};
				next if($mutBase eq '');
				if( $mutType ne 'del' ){
					$mut_vseq .= $mutBase;
				}
				if( $mutType eq 'ins' ){
					$mutBase = substr($mutBase, 1);
				}
				print MLT "$vname\t$pos\t$mutType\t$mutBase\n";
			}
		}
		close MLT;
		# create the mutated virus fasta reference file
		my $mutated_virus_fa = FuseSV::GetPath->GetPath( filekey => 'simu_virus_mutated_fa', virus_subtype => $vname );
		write_fasta_file(
			SeqSref	=> \$mut_vseq,
			FaFile	=> $mutated_virus_fa,
			SegName	=> "$vname-mut",
			LineBase => 50,
			CircleExtLen => ($virusInfo_Href->{$vname}->{circle} =~ /yes/) ? 500 : 0
		);
	}
}

#--- load virus information
sub load_virus_info{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $virusInfo_file = $parm{virusInfo_file};
	my $virusInfo_Href = $parm{virusInfo_Href};
	my $needvirus_Href = $parm{needvirus_Href};

	open (VINFO, Try_GZ_Read($virusInfo_file)) || die "fail to read $virusInfo_file: $!\n";
	while(<VINFO>){
		next if(/^#/);
		my @ele = split;
		my $vname = $ele[0];
		next if(defined($needvirus_Href) && !exists($needvirus_Href->{$vname}));
		$virusInfo_Href->{$vname} = { length=>$ele[2], circle=>$ele[3], exdlen=>$ele[4],
									  orig_seq=>'', mutations=>{} 
									};
	}
	close VINFO;
}

#--- read virus genome fasta
sub load_virus_genome{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $virus_orig_genome_fa = $parm{virus_orig_fa};
	my $virusInfo_Href = $parm{virusInfo_Href};

	open (VOGF, Try_GZ_Read($virus_orig_genome_fa)) || die "fail read $virus_orig_genome_fa: $!\n";
	$/ = '>'; <VOGF>; $/ = "\n";
	while(<VOGF>){
		chomp(my $refseg = $_);
		$/ = '>';
		chomp(my $segseq = <VOGF>);
		$/ = "\n";
		$segseq =~ s/\s//g;
		# check
		if(!exists($virusInfo_Href->{$refseg})){
			# warn "The virus ($refseg) is not in the virus_info file. Skip it.\n";
			next;
		}
		# load
		$virusInfo_Href->{$refseg}->{orig_seq} = lc($segseq);
	}
	close VOGF;
}

#--- generate mutations for virus
sub mutate_virus{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $snv_freq = $parm{snv_freq} || 0.01;
	my $indel_freq = $parm{indel_freq} || 0.001;
	my $MNP_freq = $parm{MNP_freq} || 0.001;
	my $indel_localDiff = $parm{indel_localDiff} || 1; # at least 1
	my $max_InDel_Length = $parm{max_InDel_Length} || 15;
	my $max_MNP_Length = $parm{max_MNP_Length} || 3;
	my $virusInfo_Href = $parm{virusInfo_Href};

	for my $vname (keys %$virusInfo_Href){
		my $vlength = $virusInfo_Href->{$vname}->{length};
		my $mut_Href = $virusInfo_Href->{$vname}->{mutations};
		my @ref_base = split //, $virusInfo_Href->{$vname}->{orig_seq};
		#****************
		# mutation: snv
		#****************
		warn "generate SNV randomly for virus $vname.\n";
		my $snv_amount = int($vlength * $snv_freq);
		while( $snv_amount != 0 ){
			my ($rand_pos, $mutBase) = &rand_mutation(
										ref_range=>$vlength,
										mut_type=>'snv',
										refBase_Aref=>\@ref_base
									   );
			# new alt pos, no neighboors
			my $mut_bool = 0;
			$mut_bool ||= exists($mut_Href->{$rand_pos+$_}) for (-1 .. 1);
			next if($mut_bool);
			$snv_amount--;
			$mut_Href->{$rand_pos} = {type=>'snv', mutbase=>$mutBase};
		}
		#****************
		# mutation: INDEL
		#****************
		warn "generate INDEL randomly for virus $vname.\n";
		my $indel_amount = int($vlength * $indel_freq);
		while( $indel_amount != 0 ){
			my $type = (rand(2)>1) ? 'ins':'del';
			if($type eq 'ins'){
				my ($rand_pos, $insBase) = &rand_mutation(
											ref_range=>$vlength,
											mut_type=>$type,
											refBase_Aref=>\@ref_base,
											max_InDel_Length=>$max_InDel_Length
										   );
				# new alt pos, no neighboors
				my $mut_bool = 0;
				$mut_bool ||= exists($mut_Href->{$rand_pos+$_}) for (-1 .. 1);
				next if($mut_bool);
				# avoid indel shift because local identical bases
				my $localSame_bool = 0;
				## head of ins seq cmp with bases behind the pos
				# $localSame_bool ||= (uc($ref_base[$rand_pos+$_-1]) eq uc(substr($insBase,1*$_-1,1))) for (1 .. min($indel_localDiff,length($insBase)));
				## tail of ins seq cmp with bases before the pos (including)
				# $localSame_bool ||= (uc($ref_base[$rand_pos-$_+0]) eq uc(substr($insBase,-1*$_+0,1))) for (1 .. min($indel_localDiff,length($insBase)));
				for my $end ('head', 'tail'){
					my ($factor, $offset) = ($end eq 'head')?(1,-1):(-1,0);
					for (1 .. min($indel_localDiff,length($insBase))){
						my $locBaseForCmp_i = $rand_pos + $factor * $_ + $offset;
						my $insBaseForCmp_i = $factor * $_ + $offset;
						# index out of range? not allowed
						if( $locBaseForCmp_i < 0 || $locBaseForCmp_i >= $vlength ){
							$localSame_bool = 1;
							last;
						}
						$localSame_bool ||= (uc($ref_base[$locBaseForCmp_i]) eq uc(substr($insBase, $insBaseForCmp_i, 1)));
						last if($localSame_bool);
					}
					last if($localSame_bool);
				}
				# bilateral identical bases, not allowed
				next if($localSame_bool);
				## accept this insertion
				$indel_amount--;
				$mut_Href->{$rand_pos} = {type=>$type, mutbase=>$ref_base[$rand_pos-1].$insBase};
			}
			else{
				my ($rand_pos, $del_Len) = &rand_mutation(
											ref_range=>$vlength,
											mut_type=>$type,
											refBase_Aref=>\@ref_base,
											max_InDel_Length=>$max_InDel_Length
										   );
				# new alt pos, no neighboors
				my $mut_bool = 0;
				$mut_bool ||= exists($mut_Href->{$rand_pos+$_}) for (-1 .. $del_Len);
				next if($mut_bool);
				# avoid indel shift because local identical bases
				my $localSame_bool = 0;
				## head of del seq cmp with bases behind the last del-pos
				# $localSame_bool ||= (uc($ref_base[$rand_pos+$_-2]) eq uc($ref_base[$rand_pos+$del_Len+$_-2])) for (1 .. min($indel_localDiff,$del_Len));
				## tail of del seq cmp with bases before the first del-pos
				# $localSame_bool ||= (uc($ref_base[$rand_pos-$_-1]) eq uc($ref_base[$rand_pos+$del_Len-$_-1])) for (1 .. min($indel_localDiff,$del_Len));
				for my $end ('head', 'tail'){
					my ($factor, $offset) = ($end eq 'head')?(1,-2):(-1,-1);
					for (1 .. min($indel_localDiff,$del_Len)){
						my $locBaseForCmp_i = $rand_pos + $factor * $_ + $offset;
						my $delBaseForCmp_i = $rand_pos + $factor * $_ + $offset + $del_Len;
						# index out of range? not allowed
						if(    $locBaseForCmp_i < 0 || $locBaseForCmp_i >= $vlength
							|| $delBaseForCmp_i < 0 || $delBaseForCmp_i >= $vlength
						  ){
							$localSame_bool = 1;
							last;
						}
						$localSame_bool ||= (uc($ref_base[$locBaseForCmp_i]) eq uc($ref_base[$delBaseForCmp_i]));
						last if($localSame_bool);
					}
					last if($localSame_bool);
				}
				## bilateral identical bases, not allowed
				next if($localSame_bool);
				# accept this deletion
				$indel_amount--;
				$mut_Href->{$rand_pos} = { type=>$type, mutbase=>join('', @ref_base[$rand_pos-1 .. $rand_pos+$del_Len-2]) };
				$mut_Href->{$rand_pos+$_} = {type=>$type, mutbase=>''} for (1 .. $del_Len-1);
			}
		}
		#****************
		# mutation: MNP
		#****************
		warn "generate MNP randomly for virus $vname.\n";
		my $mnp_amount = int($vlength * $MNP_freq);
		while( $mnp_amount != 0 ){
			my ($rand_pos, $mutBase) = &rand_mutation(
										ref_range=>$vlength,
										mut_type=>'mnp',
										refBase_Aref=>\@ref_base,
										max_MNP_Length=>$max_MNP_Length
									   );
			# new alt pos, no neighboors
			my $mnp_len = length($mutBase);
			my $mut_bool = 0;
			$mut_bool ||= exists($mut_Href->{$rand_pos+$_}) for (-1 .. $mnp_len);
			next if($mut_bool);
			$mnp_amount--;
			$mut_Href->{$rand_pos} = {type=>'mnp', mutbase=>$mutBase};
			$mut_Href->{$rand_pos+$_} = {type=>'mnp', mutbase=>''} for (1 .. $mnp_len-1);
		}
	}
}

#--- randomly arrange mutations based on given frequency
sub rand_mutation{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $ref_range = $parm{ref_range};
	my $mut_type  = $parm{mut_type};
	my $refBase_Aref  = $parm{refBase_Aref};
	my $max_MNP_Length = $parm{max_MNP_Length} || 3;
	my $max_InDel_Length = $parm{max_InDel_Length} || 15;

	my $rand_pos = int(rand($ref_range))+1;

	if( $mut_type =~ /snv/i ){
		my $mutBase = &rand_snv( refbase=>$refBase_Aref->[$rand_pos-1] );
		return ($rand_pos, $mutBase);
	}
	elsif( $mut_type =~ /ins/i ){
		my $insBase = &rand_ins( max_InDel_Length=>$max_InDel_Length );
		return ($rand_pos, $insBase);
	}
	elsif( $mut_type =~ /del/i ){
		my $del_Len = int(rand($max_InDel_Length))+1;
		return ($rand_pos, $del_Len);
	}
	elsif( $mut_type =~ /mnp/i ){
		my $mnpLen = int(rand($max_MNP_Length-1))+2;
		my $mutBase = '';
		$mutBase .= &rand_snv( refbase=>$refBase_Aref->[$rand_pos+$_-1] ) for (0..$mnpLen-1);
		return ($rand_pos, $mutBase);
	}
}

#--- randomly generate a snv
sub rand_snv{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $refbase = $parm{refbase};

	(my $candbase = 'ACGT') =~ s/$refbase//i;
	my $mutbase = (split //,$candbase)[int(rand(3))];
	return $mutbase;
}

#--- randomly generate inserted sequence
sub rand_ins{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $max_InDel_Length = $parm{max_InDel_Length} || 15;

	my $ins_len = int(rand($max_InDel_Length))+1;
	my $mutbase = '';
	$mutbase .= (split //,'ACGT')[int(rand(4))] for (1..$ins_len);
	return $mutbase;
}

#--- 
1; ## tell the perl script the successful access of this module.
