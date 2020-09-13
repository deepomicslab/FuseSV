package FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use BioFuse::BioInfo::FASTA qw/ write_fasta_file /;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  filter_refpos_mutation_from_VCF
			  reconstruct_virus_genome_BaseOn_mutation
			  merge_circular_virus_extend_part
			  load_refpos_allel_from_ref_v_fa
			/;
push @EXPORT , qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ModifyVrefBaseOnRefposInfo';
#----- version --------
$VERSION = "0.11";
$DATE = '2019-05-02';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						filter_refpos_mutation_from_VCF
						reconstruct_virus_genome_BaseOn_mutation
						merge_circular_virus_extend_part
						load_refpos_allel_from_ref_v_fa
					 /;

#--- filter refpos and mutations from VCF ---
sub filter_refpos_mutation_from_VCF{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{RefPos_Href};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $min_depth = $V_Href->{min_depth};
	my $both_strd = $V_Href->{both_strd};
	my $min_mut_reads = $V_Href->{min_mut_reads};
	my $adjacent_indel_min_dist = $V_Href->{adjacent_indel_min_dist};

	# filter mutations
	# basic criterion
	# adjacent indel exclusion
	my $last_indel_pos;
	for my $sbjct_pos (sort {$a<=>$b} keys %$RefPos_Href){
		my $sbjct_pos_Href = $RefPos_Href->{$sbjct_pos};
		my $ref_allel = $sbjct_pos_Href->{ref_allel};
		# my $makeup_allel = $sbjct_pos_Href->{makeup_allel};
		my $depth = $sbjct_pos_Href->{depth};
		my $refallel_depth = $sbjct_pos_Href->{refallel_depth};
		# sweep secondary mutations
		if(   $depth < $min_depth   # depth filtration
		   || $ref_allel =~ /^N$/i  # do not allow mutation at 'N' locus
		  ){
			delete $sbjct_pos_Href->{mutation}; # delete it anyway
		}
		elsif( exists $sbjct_pos_Href->{mutation} ){ # has mutation(s), so check whether it is good enough
			my $mutation_Href = $sbjct_pos_Href->{mutation};
			# select the best available mutation
			my $find_good_mut_bool = 0;
			my $best_mut_info_Href = { Sum_Sup => 0, Plus_Sup => 0, Minus_Sup => 0, Mut_Type => '', Mut_Seq => '' };
			# filter all mutations
			for my $mutation_id (sort keys %$mutation_Href){
				my ($mut_type, $mut_seq) = split /,/, $mutation_id;
				my ($sum_sup, $plus_sup, $minus_sup, $refSup_toCmp) = @{ $mutation_Href->{$mutation_id} };
				# filter the mutation
				if(    ( $sum_sup < $min_mut_reads ) # at least so many supported reads
				    || ( $both_strd && ($plus_sup < $both_strd || $minus_sup < $both_strd) ) # if required, both strand supported
				    || ( $mut_type !~ /del/ && $sum_sup <= $refSup_toCmp )   # must be dominant against ref-allel, for non-del
				    || ( $mut_type =~ /del/ && $sum_sup <= $refallel_depth ) # must be dominant against ref-allel, for del specifically
				  ){
					next;
				}
				elsif($sum_sup > $best_mut_info_Href->{Sum_Sup}){
					# update bool
					$find_good_mut_bool = 1;
					# information update
					$best_mut_info_Href = { Sum_Sup=>$sum_sup, Plus_Sup=>$plus_sup, Minus_Sup=>$minus_sup, Mut_Type=>$mut_type, Mut_Seq=>$mut_seq, Method=>'Bwa-SamTools-VCF' };
				}
			}
			# has good mutations?
			if( !$find_good_mut_bool ){ # no
				delete $sbjct_pos_Href->{mutation};
			}
			else{ # yes, update to best mutation
				$sbjct_pos_Href->{mutation} = $best_mut_info_Href; # here, change the mutation_Href structure
				# compare adjacent indel supports
				if( $sbjct_pos_Href->{mutation}->{Mut_Type} =~ /ins|del/ ){
					if(    defined($last_indel_pos)
						&& $sbjct_pos - $last_indel_pos < $adjacent_indel_min_dist
					  ){
						my $last_indel_Sum_Sup = $RefPos_Href->{$last_indel_pos}->{mutation}->{Sum_Sup};
						my $this_indel_Sum_Sup = $sbjct_pos_Href->{mutation}->{Sum_Sup};
						if( $last_indel_Sum_Sup <= $this_indel_Sum_Sup ){
							delete $RefPos_Href->{$last_indel_pos}->{mutation};
						}
						else{
							delete $sbjct_pos_Href->{mutation};
						}
					}
					else{
						$last_indel_pos = $sbjct_pos;
					}
				}
			}
		}
	}
}

#--- reconstruct virus genome base on mutation list ---
sub reconstruct_virus_genome_BaseOn_mutation{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{RefPos_Href};
	my $modify_vfa = $Option_Href->{Mod_V_G};
	my $modify_vfa_refname = $Option_Href->{Mod_V_G_refname};
	my $modify_makeup_vfa = $Option_Href->{Mod_makeup_V_G};
	my $use_N_for_lowdepth_sign = $Option_Href->{Use_N_for_lowdepth};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $min_depth = $V_Href->{min_depth};
	my $line_base = $V_Href->{line_base};

	# each refpos and its mutations
	# my $ever_modify_sign = 0;
	my $has_low_depth_pos = 0;
	my $has_mutation_sign = 0;
	my $modified_ref_seq = '';
	my $modified_ref_seq_for_makeup = '';
	my $sbjct_filter_pos = 1;
	my %alt_info = ( sum=>{snp=>0, ins=>0, del=>0}, altpos=>{}, N_count=>0 );
	for my $sbjct_pos (sort {$a<=>$b} keys %$RefPos_Href){
		my $sbjct_pos_Href = $RefPos_Href->{$sbjct_pos};
		# info
		my $ref_allel = $sbjct_pos_Href->{ref_allel};
		my $makeup_allel = $sbjct_pos_Href->{makeup_allel};
		my $depth = $sbjct_pos_Href->{depth};
		my $refallel_depth = $sbjct_pos_Href->{refallel_depth};
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
				# $ever_modify_sign = 1; # record it should be modified
				$has_low_depth_pos = 1;
				$alt_info{N_count} ++ if( $use_N_for_lowdepth_sign );
			}
			$sbjct_filter_pos = $sbjct_pos + 1;
		}
		else{
			if(    $ref_allel !~ /^N$/i                # do not allow mutation at 'N' locus
				&& exists $sbjct_pos_Href->{mutation}  # has the best mutation
			  ){
				my $mutation_Href = $sbjct_pos_Href->{mutation};
				# use the best mutations to modify the genome
				# here, update modify the ref sign
				# $ever_modify_sign = 1;
				$has_mutation_sign = 1;
				# the mutation information
				my $mut_type = $mutation_Href->{Mut_Type};
				my $mut_seq  = $mutation_Href->{Mut_Seq};
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
					Mut_Sum_Sup   => $mutation_Href->{Sum_Sup},
					Mut_Plus_Sup  => $mutation_Href->{Plus_Sup},
					Mut_Minus_Sup => $mutation_Href->{Minus_Sup},
					Mut_Type      => $mutation_Href->{Mut_Type},
					Mut_Seq       => $mutation_Href->{Mut_Seq},
					Method        => $mutation_Href->{Method}
				};
			}
			else{ # no mutations at this refpos
				$modified_ref_seq .= $ref_allel;
				$modified_ref_seq_for_makeup .= $makeup_allel;
				$sbjct_filter_pos = $sbjct_pos + 1;
			}
		}
	}

	# create modified virus genome fasta file no matter it is changed or same
	write_fasta_file(
		SeqSref	=> \$modified_ref_seq,
		FaFile	=> $modify_vfa,
		SegName	=> $modify_vfa_refname,
		LineBase => $V_Href->{line_base},
		CircleExtLen => 0 # no matter circle or not, this part should be original, donot extend here
	);

	# create modified virus genome fasta file, the makeup version
	write_fasta_file(
		SeqSref	=> \$modified_ref_seq_for_makeup,
		FaFile	=> $modify_makeup_vfa,
		SegName	=> $modify_vfa_refname,
		LineBase => $V_Href->{line_base},
		CircleExtLen => 0 # no matter circle or not, this part should be original, donot extend here
	);

	# return changed or not
	my $orig_length = length($modified_ref_seq);
	$modified_ref_seq =~ s/[Nn]//g;
	my $nonN_length = length($modified_ref_seq);
	# return ( [$ever_modify_sign, $has_mutation_sign], \%alt_info, $orig_length, $nonN_length );
	return {
		ever_modify_sign  => ($has_low_depth_pos || $has_mutation_sign),
		has_low_depth_pos => $has_low_depth_pos,
		has_mutation_sign => $has_mutation_sign,
		alt_info_Href => \%alt_info,
		orig_length => $orig_length,
		nonN_length => $nonN_length
	};
}

#--- specially merge info in extended region for circular virus ---
sub merge_circular_virus_extend_part{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{RefPos_Href};
	my $v_ExtLen = $Option_Href->{Ext_Len};
	my $consider_region = $Option_Href->{Consider_Region} || $v_ExtLen;

	# get the maximum position in reference
	my $MaxPosInRef = max(keys %$RefPos_Href);

	for my $rel_pos (1 .. $v_ExtLen){
		my $aim_RefPos_Href = $RefPos_Href->{$rel_pos};
		my $sor_RefPos = $MaxPosInRef - $v_ExtLen + $rel_pos;

		# merge
		if(    $rel_pos <= $consider_region
			&& exists( $RefPos_Href->{$sor_RefPos} ) # must has the information of the pos in the extended part
		){
			my $sor_RefPos_Href = $RefPos_Href->{$sor_RefPos};
			$aim_RefPos_Href->{depth} += $sor_RefPos_Href->{depth};
			$aim_RefPos_Href->{refallel_depth} += $sor_RefPos_Href->{refallel_depth};
			$aim_RefPos_Href->{refallel_D_Fw} += $sor_RefPos_Href->{refallel_D_Fw};
			$aim_RefPos_Href->{refallel_D_Rv} += $sor_RefPos_Href->{refallel_D_Rv};
			for my $mutation_id (sort keys %{$sor_RefPos_Href->{mutation}}){
				if(!exists($aim_RefPos_Href->{mutation}->{$mutation_id})){
					$aim_RefPos_Href->{mutation}->{$mutation_id} = [ 0, 0, 0, 0 ];
				}
				# merge
				for my $i (0 .. 3){
					$aim_RefPos_Href->{mutation}->{$mutation_id}->[$i] += $sor_RefPos_Href->{mutation}->{$mutation_id}->[$i];
				}
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

	my $orig_ref_fa = $Option_Href->{Ref};

	open (REF, $orig_ref_fa) || die "fail read $orig_ref_fa: $!\n";
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
														refallel_D_Fw => 0,
														refallel_D_Rv => 0,
														mutation => {}
													}   for ( 1 .. scalar(@segseq) );
			}
			else{
				# load makeup fasta
				# specific because makeup ref seq requires non-N
				$Option_Href->{RefPos_Href}->{$_}->{makeup_allel} = $segseq[$_-1] for ( 1 .. scalar(@segseq) );
			}
		}
		elsif( exists($Option_Href->{OutFaFile}) ){ # function 2nd
			write_fasta_file(
				SeqSref	=> \$segseq,
				FaFile => $Option_Href->{OutFaFile},
				SegName	=> $Option_Href->{SegName},
				LineBase => $Option_Href->{LineBase},
				CircleExtLen => ($Option_Href->{CircleExtLen} || 0),
				split_N => ($Option_Href->{split_N} || 0)
			);
		}
	}
	close REF;

	# for the makeup ref
	# specific because makeup ref seq requires non-N
	if(    exists($Option_Href->{makeup_Ref})
		&& defined($Option_Href->{makeup_Ref})
	    && -e $Option_Href->{makeup_Ref}
	  ){
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
