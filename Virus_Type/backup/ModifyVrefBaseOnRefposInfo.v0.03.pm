package FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo;

use strict;
use warnings;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  filter_refpos_and_reconstruct_virus_genome
			  merge_circular_virus_extend_part
			  load_refpos_allel_from_ref_v_fa
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ModifyVrefBaseOnRefposInfo';
#----- version --------
$VERSION = "0.03";
$DATE = '2016-04-27';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						filter_refpos_and_reconstruct_virus_genome
						merge_circular_virus_extend_part
						load_refpos_allel_from_ref_v_fa
					 /;

#--- filter refpos and mutations, reconstruct virus genome ---
sub filter_refpos_and_reconstruct_virus_genome{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{RefPos_Href};
	my $modify_vfa = $Option_Href->{Mod_V_G};
	my $modify_vfa_refname = $Option_Href->{Mod_V_G_refname};
	my $modify_makeup_vfa = $Option_Href->{Mod_makeup_V_G};
	my $use_N_for_lowdepth_sign = $Option_Href->{Use_N_for_lowdepth};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $min_depth = $V_Href->{min_depth};
	my $both_strd = $V_Href->{both_strd};
	my $min_mut_reads = $V_Href->{min_mut_reads};
	my $line_base = $V_Href->{line_base};

	# each refpos and its mutations
	my $ever_modify_sign = 0;
	my $has_mutation_sign = 0;
	my $modified_ref_seq = '';
	my $modified_ref_seq_for_makeup = '';
	my $sbjct_filter_pos = 1;
	my %alt_info = ( sum=>{snp=>0, ins=>0, del=>0}, altpos=>{}, N_count=>0 );
	for my $sbjct_pos (sort {$a<=>$b} keys %$RefPos_Href){
		my $ref_allel = $RefPos_Href->{$sbjct_pos}->{ref_allel};
		my $makeup_allel = $RefPos_Href->{$sbjct_pos}->{makeup_allel};
		my $depth = $RefPos_Href->{$sbjct_pos}->{depth};
		my $refallel_depth = $RefPos_Href->{$sbjct_pos}->{refallel_depth};
		# filter the sbjct_pos
		if($sbjct_pos < $sbjct_filter_pos){
			next;
		}
		# output
		if($depth < $min_depth){ # depth filtration
			$modified_ref_seq .= ( $use_N_for_lowdepth_sign ? 'N' : $ref_allel );
			$modified_ref_seq_for_makeup .= $makeup_allel;
			# the original ref-allele is not 'N', so must be changed now.
			if($ref_allel !~ /^N$/i){
				$ever_modify_sign = 1; # record it should be modified
				$alt_info{N_count} ++ if( $use_N_for_lowdepth_sign );
			}
			$sbjct_filter_pos = $sbjct_pos + 1;
		}
		else{
			if(    $ref_allel !~ /^N$/i                           # do not allow mutation at 'N' locus
				&& exists $RefPos_Href->{$sbjct_pos}->{mutation}  # has mutation(s), so check whether it is good enough
			  ){
				my $mutation_Href = $RefPos_Href->{$sbjct_pos}->{mutation};
				# select the best available mutation
				my $find_good_mut_bool = 0;
				my $best_mut_info_Href = {
					Sum_Sup => 0,
					Plus_Sup => 0,
					Minus_Sup => 0,
					Mut_Type => '',
					Mut_Seq => ''
				};
				# filter all mutations
				for my $mutation_id (sort keys %$mutation_Href){
					my ($sum_sup, $plus_sup, $minus_sup) = @{ $mutation_Href->{$mutation_id} };
					# filter the mutation
					if(    ( $sum_sup < $min_mut_reads ) # at least so many supported reads
					    || ( $both_strd && ($plus_sup < $both_strd || $minus_sup < $both_strd) ) # if required, both strand supported
					    || ( $sum_sup <= $refallel_depth ) # must be more than ref-allel
					  ){
						next;
					}
					elsif($sum_sup > $best_mut_info_Href->{Sum_Sup}){
						# update bool
						$find_good_mut_bool = 1;
						# information update
						my ($mut_type, $mut_seq) = split /,/, $mutation_id;
						$best_mut_info_Href->{Sum_Sup}   = $sum_sup;
						$best_mut_info_Href->{Plus_Sup}  = $plus_sup;
						$best_mut_info_Href->{Minus_Sup} = $minus_sup;
						$best_mut_info_Href->{Mut_Type}  = $mut_type;
						$best_mut_info_Href->{Mut_Seq}   = $mut_seq;
					}
				}
				# has good mutations?
				if( !$find_good_mut_bool ){ # no good mutations
					$modified_ref_seq .= $ref_allel;
					$modified_ref_seq_for_makeup .= $makeup_allel;
					$sbjct_filter_pos = $sbjct_pos + 1;
				}
				else{ # use the best mutations to modify the genome
					# here, update modify the ref sign
					$ever_modify_sign = 1;
					$has_mutation_sign = 1;
					# the mutation information
					my $mut_type = $best_mut_info_Href->{Mut_Type};
					my $mut_seq  = $best_mut_info_Href->{Mut_Seq};
					# different types, different modification methods
					if($mut_type eq 'snp'){
						$modified_ref_seq .= $mut_seq;
						$modified_ref_seq_for_makeup .= $mut_seq;
						$sbjct_filter_pos = $sbjct_pos + 1;
					}
					elsif($mut_type eq 'ins'){
						$modified_ref_seq .= $ref_allel.$mut_seq;
						$modified_ref_seq_for_makeup .= $makeup_allel.$mut_seq;
						$sbjct_filter_pos = $sbjct_pos + 1;
					}
					elsif($mut_type eq 'del'){
						#$modified_ref_seq .= ''; # deletion, no output
						#$modified_ref_seq_for_makeup .= '';
						$sbjct_filter_pos = $sbjct_pos + length($mut_seq);
					}
					# record alt count
					$alt_info{sum}->{$mut_type}++;
					$alt_info{altpos}->{$sbjct_pos} = {
						Position      => $sbjct_pos,
						Ref_Allel     => $ref_allel,
						Ref_Allel_Sup => $refallel_depth,
						Mut_Sum_Sup   => $best_mut_info_Href->{Sum_Sup},
						Mut_Plus_Sup  => $best_mut_info_Href->{Plus_Sup},
						Mut_Minus_Sup => $best_mut_info_Href->{Minus_Sup},
						Mut_Type      => $best_mut_info_Href->{Mut_Type},
						Mut_Seq       => $best_mut_info_Href->{Mut_Seq}
					};
				}
			}
			else{ # no mutations at this refpos
				$modified_ref_seq .= $ref_allel;
				$modified_ref_seq_for_makeup .= $makeup_allel;
				$sbjct_filter_pos = $sbjct_pos + 1;
			}
		}
	}

	# create modified virus genome fasta file no matter it is changed or same
	SOAPfuse::General_Operation->write_fasta_file(
			{
				SeqSref	=> \$modified_ref_seq,
				FaFile	=> $modify_vfa,
				SegName	=> $modify_vfa_refname,
				LineBase => $V_Href->{line_base},
				CircleExtLen => 0 # no matter circle or not, this part should be original, donot extend here
			}
		);

	# create modified virus genome fasta file, the makeup version
	SOAPfuse::General_Operation->write_fasta_file(
			{
				SeqSref	=> \$modified_ref_seq_for_makeup,
				FaFile	=> $modify_makeup_vfa,
				SegName	=> $modify_vfa_refname,
				LineBase => $V_Href->{line_base},
				CircleExtLen => 0 # no matter circle or not, this part should be original, donot extend here
			}
		);

	# return changed or not
	my $orig_length = length($modified_ref_seq);
	$modified_ref_seq =~ s/[Nn]//g;
	my $nonN_length = length($modified_ref_seq);
	return ( [$ever_modify_sign, $has_mutation_sign], \%alt_info, $orig_length, $nonN_length );
}

#--- specially merge info in extended region for circular virus ---
sub merge_circular_virus_extend_part{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{RefPos_Href};
	my $v_ExtLen = $Option_Href->{Ext_Len};

	# get the maximum position in reference
	my $MaxPosInRef = max(keys %$RefPos_Href);

	for my $rel_pos (1 .. $v_ExtLen){
		my $aim_RefPos_Href = $RefPos_Href->{$rel_pos};
		my $sor_RefPos = $MaxPosInRef - $v_ExtLen + $rel_pos;
		my $sor_RefPos_Href = $RefPos_Href->{$sor_RefPos};

		# check, debug
		if(!exists($aim_RefPos_Href->{depth})){
			warn "aim_RefPos_Href\t$rel_pos\n";
			exit(0);
		}
		if(!exists($sor_RefPos_Href->{depth})){
			warn "sor_RefPos_Href\t$sor_RefPos\n";
			exit(0);
		}

		# merge
		$aim_RefPos_Href->{depth} += $sor_RefPos_Href->{depth};
		$aim_RefPos_Href->{refallel_depth} += $sor_RefPos_Href->{refallel_depth};
		for my $mutation_id (sort keys %{$sor_RefPos_Href->{mutation}}){
			if(!exists($aim_RefPos_Href->{mutation}->{$mutation_id})){
				$aim_RefPos_Href->{mutation}->{$mutation_id} = [0,0,0];
			}
			# merge
			for my $i (0, 1, 2){
				$aim_RefPos_Href->{mutation}->{$mutation_id}->[$i] += $sor_RefPos_Href->{mutation}->{$mutation_id}->[$i];
			}
		}
		# sweep the extended pos
		delete $RefPos_Href->{$sor_RefPos};
	}
}

#--- load ref-allel ---
## two functions
sub load_refpos_allel_from_ref_v_fa{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $origi_fa_ref = $Option_Href->{Ref};

	open (REF, $origi_fa_ref) || die "fail read $origi_fa_ref: $!\n";
	$/ = '>'; <REF>; $/ = "\n";
	while(<REF>){
		chomp(my $refseg = $_);
		$/ = '>';
		chomp(my $segseq = <REF>);
		$/ = "\n";
		$segseq =~ s/\s//g;
		# do as the function required
		if( exists($Option_Href->{RefPos_Href}) ){ # function 1st
			my @segseq = split //, $segseq;
			if( !exists($Option_Href->{makeup_sign}) ){ # load the original fasta
				$Option_Href->{RefPos_Href}->{$_} = {
														ref_allel => $segseq[$_-1],
														depth => 0,
														refallel_depth => 0,
														mutation => {}
													}   for ( 1 .. scalar(@segseq) );
			}
			else{ # load makeup fasta
				$Option_Href->{RefPos_Href}->{$_}->{makeup_allel} = $segseq[$_-1] for ( 1 .. scalar(@segseq) );
			}
		}
		elsif( exists($Option_Href->{OutFaFile}) ){ # function 2nd
			SOAPfuse::General_Operation->write_fasta_file(
					{
						SeqSref	=> \$segseq,
						FaFile => $Option_Href->{OutFaFile},
						SegName	=> $Option_Href->{SegName},
						LineBase => $Option_Href->{LineBase},
						CircleExtLen => ($Option_Href->{CircleExtLen} || 0),
						split_N => ($Option_Href->{split_N} || 0)
					}
				);
		}
	}
	close REF;

	# for the makeup ref
	if( exists($Option_Href->{makeup_Ref}) ){
		&load_refpos_allel_from_ref_v_fa(
			{
				Ref => $Option_Href->{makeup_Ref},
				makeup_sign => 1,
				RefPos_Href => $Option_Href->{RefPos_Href}
			}
		);
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
