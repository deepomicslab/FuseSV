package FuseSV::Virus_Type::BlastM0Mutation;

use strict;
use warnings;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr trible_run_for_success/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  format_database
			  blastall_query_m0
			  read_blast_m0_for_mutation
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BlastM0Mutation';
#----- version --------
$VERSION = "0.15";
$DATE = '2017-06-20';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						format_database
						blastall_query_m0
						read_blast_m0_for_mutation
						upload_mutation_list
						summary_alignment_result
						load_refpos_allel_from_ref_fa
						read_query_mapinfo_from_blast_m0
						adjust_ins_after_mismatch
						concatenate_GapInducedByLargeInDel
					 /;

#--- construct blastall database ---
sub format_database{
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $ref_fa_file = $Option_Href->{ref_fa_file};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# construct blast reference
	my $formatdb_cmd = "$V_Href->{formatdb} -i $ref_fa_file -p F -l $ref_fa_file.log";
	trible_run_for_success($formatdb_cmd, 'formatdb_VirusRef', {esdo_Nvb=>1});
}

#--- blastall queries database ---
sub blastall_query_m0{
	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $ref_fa_file = $Option_Href->{ref_fa_file};
	my $query_fafile = $Option_Href->{query_fa_file};
	my $blast_m0 = $Option_Href->{blast_m0};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# construct blast reference
	my $blastall_cmd = "$V_Href->{blastall} -m 0 -a 8 -p blastn -F F -e 1e-06 -i $query_fafile -d $ref_fa_file -o $blast_m0";
	trible_run_for_success($blastall_cmd, 'blastall_VirusRef', {esdo_Nvb=>1});
}

#--- read blast m0 for mutation ---
sub read_blast_m0_for_mutation{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $blast_m0 = $Option_Href->{blast_m0};
	my $query_seq_fa = $Option_Href->{query_seq_fa} || '';
	my $queryid_list = $Option_Href->{queryid_list} || '';
	my $use_all_map = $Option_Href->{use_all_map} || 0;
	my $orig_fa_ref = $Option_Href->{orig_fa_ref};
	my $orig_ext = $Option_Href->{orig_ext} || 0;
	my $extend_merge_Consider_Region = $Option_Href->{ext_merge_Consider_Region} || 20;
	my $min_depth = $Option_Href->{min_depth} || 1;
	my $both_strd = $Option_Href->{both_strd} || 0;
	my $min_mut_reads = $Option_Href->{min_mut_reads} || 2;
	my $mutation_Href = $Option_Href->{mutation_Href};
	my $need_RefPos_Hash = $Option_Href->{need_RefPos_Hash} || 0;
	my $without_supports = $Option_Href->{without_supports} || 0;

	# read blast m0 result
	my %Query_MapInfo;
	&read_query_mapinfo_from_blast_m0({
		blast_m0 => $blast_m0,
		query_seq_fa => $query_seq_fa, # maybe empty string
		sbjct_seq_fa => $orig_fa_ref,
		queryid_list => $queryid_list, # maybe empty string
		query_mapinfo_Href => \%Query_MapInfo
	});

	# construct the basic info hash of all ref pos
	my %RefPos;
	&load_refpos_allel_from_ref_fa({
		orig_fa_ref => $orig_fa_ref,
		refpos_Href => \%RefPos
	});

	# load mutations
	&summary_alignment_result({
		query_mapinfo_Href => \%Query_MapInfo,
		refpos_Href => \%RefPos,
		only_first_align => !$use_all_map
	});

	# merge bases in the extended part
	if( $orig_ext != 0 ){
		FuseSV::Virus_Type::ModifyVrefBaseOnRefposInfo->merge_circular_virus_extend_part({
			RefPos_Href => \%RefPos,
			# pos in this part will be deleted
			Ext_Len => $orig_ext,
			# only such long region forefront of circular virus genome to merge
			Consider_Region => $extend_merge_Consider_Region
		});
	}

	# feedback all qualified major mutations
	&upload_mutation_list({
		mutation_Href => $mutation_Href,
		refpos_Href => \%RefPos,
		min_depth => $min_depth,
		both_strd => $both_strd,
		min_mut_reads => $min_mut_reads,
		without_supports => $without_supports
	});

	if( $need_RefPos_Hash ){
		return \%RefPos;
	}
}

#--- define mutation or not for each base, and output ---
sub upload_mutation_list{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{refpos_Href};
	my $major_mutation_Href = $Option_Href->{mutation_Href};
	my $min_depth = $Option_Href->{min_depth};
	my $both_strd = $Option_Href->{both_strd};
	my $min_mut_reads = $Option_Href->{min_mut_reads};
	my $without_supports = $Option_Href->{without_supports} || 0;

	my $sbjct_filter_pos = 1;

	# mutations
	for my $sbjct_pos (sort {$a<=>$b} keys %$RefPos_Href){
		my $ref_allel = $RefPos_Href->{$sbjct_pos}->{ref_allel};
		my $depth = $RefPos_Href->{$sbjct_pos}->{depth} || 0;
		my $refallel_depth = $RefPos_Href->{$sbjct_pos}->{refallel_depth} || 0;
		# filter the sbjct_pos
		if($sbjct_pos < $sbjct_filter_pos){
			next;
		}
		# output
		if($depth < $min_depth){
			$sbjct_filter_pos = $sbjct_pos + 1;
		}
		else{
			if(exists $RefPos_Href->{$sbjct_pos}->{mutation}){
				my $mutation_Href = $RefPos_Href->{$sbjct_pos}->{mutation};
				# select the mutations
				my @mutation_id;
				for my $mutation_id (keys %$mutation_Href){
					# FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
					my ($sum_sup, $plus_sup, $minus_sup) = @{ $mutation_Href->{$mutation_id} };
					# filter the mutation
					if(    ( $sum_sup < $refallel_depth ) # must be dominant against ref-allel
				  	    || ( $sum_sup < $min_mut_reads )
				  	    || ( $both_strd && ($plus_sup < $both_strd || $minus_sup < $both_strd) ) # if required, both strand supported
				  	  )
					{
						next;
					}
					else{
						push @mutation_id, $mutation_id;
					}
				}
				# no mutations
				if(scalar(@mutation_id) == 0){
					$sbjct_filter_pos = $sbjct_pos + 1;
				}
				else{ # has mutations
					# sort by the sum_sup
					@mutation_id = sort { $mutation_Href->{$b}->[0] <=> $mutation_Href->{$a}->[0] } @mutation_id;
					my $best_mutation_id = $mutation_id[0];
					my ($mut_type, $mut_seq) = split /,/, $best_mutation_id; # first one
					my ($sum_sup, $plus_sup, $minus_sup) = @{ $mutation_Href->{$best_mutation_id} };
					# record mutation
					$major_mutation_Href->{$sbjct_pos} = {
														   ref_allel => $ref_allel,
														   ref_sum_sup => 0,
														   ref_plus_sup => 0,
														   ref_minus_sup => 0,
														   mut_type => $mut_type,
														   mut_seq => $mut_seq,
														   mut_id => $best_mutation_id,
														   mut_sum_sup => ($without_supports?0:$sum_sup),
														   mut_plus_sup => ($without_supports?0:$plus_sup),
														   mut_minus_sup => ($without_supports?0:$minus_sup),
														   major_mut_bool => 1,
														   other_mutation => {}
														 };
					# update filter pos
					if($mut_type eq 'snp'){
						$sbjct_filter_pos = $sbjct_pos + 1;
					}
					elsif($mut_type eq 'ins'){
						$sbjct_filter_pos = $sbjct_pos + 1;
					}
					elsif($mut_type eq 'del'){
						$sbjct_filter_pos = $sbjct_pos + length($mut_seq);
					}
				}
			}
			else{
				$sbjct_filter_pos = $sbjct_pos + 1;
			}
		}
	}
}

#--- stat the depth and mutations ---
sub summary_alignment_result{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $Query_MapInfo_real_Href = $Option_Href->{query_mapinfo_Href};
	my $RefPos_Href = $Option_Href->{refpos_Href};
	my $only_first_align = $Option_Href->{only_first_align} || 0;

	for my $query_id (sort keys %$Query_MapInfo_real_Href){
		for my $map_idx ( 0 .. scalar(@{$Query_MapInfo_real_Href->{$query_id}})-1 ){
			# only first alignment?
			last if($only_first_align && $map_idx != 0);
			my $Query_MapInfo_Href = $Query_MapInfo_real_Href->{$query_id}->[$map_idx];
			# information
			my $sbjct_strd =    $Query_MapInfo_Href->{sbjct_strd};
			my $align_strd =    $Query_MapInfo_Href->{align_strd};
			my $align_info =    $Query_MapInfo_Href->{map_info};
			my @query_info = @{ $Query_MapInfo_Href->{query_info} };
			my @sbjct_info = @{ $Query_MapInfo_Href->{sbjct_info} };
			# make the subject plus orientated
			if($sbjct_strd == -1){ # Minus
				@sbjct_info = reverse @sbjct_info;
				($sbjct_info[1] = reverse $sbjct_info[1]) =~ tr/acgtACGT/tgcaTGCA/;
				@query_info = reverse @query_info;
				($query_info[1] = reverse $query_info[1]) =~ tr/acgtACGT/tgcaTGCA/;
				$align_info = reverse $align_info;
				# check
				if($sbjct_info[2] < $sbjct_info[0]){
					warn_and_exit("wrong sbjct info:$sbjct_info[2] < $sbjct_info[0] for query_id $query_id\n");
				}
			}
			# details
			my @align_info = split //, $align_info;
			my @query_seq  = split //, $query_info[1];
			my @sbjct_seq  = split //, $sbjct_info[1];
			my $sbjct_st = $sbjct_info[0];
			my $sbjct_ed = $sbjct_info[2];
			my $last_sbjct_pos = $sbjct_st - 1;
			# judge each base, and operate according to the mutation type
			for (my $i=0; $i<=$#align_info; $i++){
				my $map_sign = $align_info[$i];
				my $query_allel = $query_seq[$i];
				my $sbjct_allel = $sbjct_seq[$i];
				# mutation types
				if($map_sign eq '|'){ # same
					my $sbjct_pos = $last_sbjct_pos + 1;
					# accumulate depth
					$RefPos_Href->{$sbjct_pos}->{depth} ++;
					$RefPos_Href->{$sbjct_pos}->{refallel_depth} ++;
					# update last sbjct_pos
					$last_sbjct_pos = $sbjct_pos;
				}
				elsif($query_allel ne '-' && $sbjct_allel ne '-'){ # mismatch
					my $sbjct_pos = $last_sbjct_pos + 1;
					#check
					if($query_allel eq $sbjct_allel){
						warn_and_exit( "Should be a Mismatch: $sbjct_pos\n"
									  ."Q: $query_allel\n"
									  ."S: $sbjct_allel\n");
					}
					# accumulate depth
					$RefPos_Href->{$sbjct_pos}->{depth} ++;
					# record mutation
					my $mutation_id = 'snp,' . uc($query_allel);
					# FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0,0];
					my $mutation_info_Aref = $RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id};
					$mutation_info_Aref->[0] ++;
					if($align_strd > 0){
						$mutation_info_Aref->[1] ++;
					}
					else{
						$mutation_info_Aref->[2] ++;
					}
					# update last sbjct_pos
					$last_sbjct_pos = $sbjct_pos;
				}
				elsif($sbjct_allel eq '-'){ # insertion
					my $sbjct_pos = $last_sbjct_pos; # insertion is assigned to the base before it.
					# find the inserted part from query sequence
					my $j = 1;
					while($sbjct_seq[$i+$j] eq '-'){ # insertion extend
						$j ++;
					}
					$j --; # adjust
					my $insert_seq = join('',@query_seq[$i .. $i+$j]);
					$i += $j; # jump to the last base of the insertion
					# do not need to accumulate depth for this pos 'again', because just added before
					#$RefPos_Href->{$sbjct_pos}->{depth} ++;
					# decrease one from 'refallel_depth' of this pos, because just added before
					$RefPos_Href->{$sbjct_pos}->{refallel_depth} --;
					# record mutation
					my $mutation_id = 'ins,' . uc($insert_seq);
					# FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0,0];
					my $mutation_info_Aref = $RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id};
					$mutation_info_Aref->[0] ++;
					if($align_strd > 0){
						$mutation_info_Aref->[1] ++;
					}
					else{
						$mutation_info_Aref->[2] ++;
					}
					# update last sbjct_pos
					$last_sbjct_pos = $sbjct_pos;
				}
				elsif($query_allel eq '-'){ # deletion
					my $sbjct_pos = $last_sbjct_pos + 1; # deletion is assigned to the first base of itself.
					my $j = 1;
					while($query_seq[$i+$j] eq '-'){ # insertion extend
						$j ++;
					}
					$j --; # adjust
					my $delete_seq = join('',@sbjct_seq[$i .. $i+$j]);
					$i += $j; # jump to the last base of the deletion
					# accumulate depth, it is just used to measure the reads' supportings of this pos
					$RefPos_Href->{$_}->{depth} ++ for ($sbjct_pos .. $sbjct_pos+$j);
					# record mutation
					my $mutation_id = 'del,' . uc($delete_seq);
					# FORMAT (Sum, Plus, Minus, RefAllelDepth[deprecated])
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0,0];
					my $mutation_info_Aref = $RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id};
					$mutation_info_Aref->[0] ++;
					if($align_strd > 0){
						$mutation_info_Aref->[1] ++;
					}
					else{
						$mutation_info_Aref->[2] ++;
					}
					# update last sbjct_pos
					$last_sbjct_pos = $sbjct_pos + $j;
				}
				else{
					warn_and_exit("Wrong align info:\n$query_info[1]\n$align_info\n$sbjct_info[1]\n");
				}
			}
			# last check
			if($last_sbjct_pos != $sbjct_ed){
				warn_and_exit("The last_sbjct_pos:$last_sbjct_pos is not the sbjct_ed:$sbjct_ed\n");
			}
		}
	}
}

#--- load ref-allel ---
## two functions, 2017-04-28, by Wenlong Jia
## provide 'refpos_Href' or 'segseq_Href'
sub load_refpos_allel_from_ref_fa{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $orig_fa_ref = $Option_Href->{orig_fa_ref};

	open (REF, $orig_fa_ref) || die "fail $orig_fa_ref: $!\n";
	$/ = '>'; <REF>; $/ = "\n";
	while(<REF>){
		chomp(my $refseg = $_);
		$/ = '>';
		chomp(my $segseq = <REF>);
		$/ = "\n";
		# filter non-needed refseg if required
		if(     exists( $Option_Href->{needed_refseg_Href} )
			&& !exists( $Option_Href->{needed_refseg_Href}->{$refseg} )
		){
			next;
		}
		# eliminate global '\n'
		$segseq =~ s/\s//g;
		# function 1st, load each pos as key
		if( exists($Option_Href->{refpos_Href}) ){
			my $RefPos_Href = $Option_Href->{refpos_Href};
			my @segseq = split //, $segseq;
			$RefPos_Href->{$_+1} = {
									ref_allel => $segseq[$_],
									depth => 0,
									refallel_depth => 0,
									refallel_D_Fw => 0,
									refallel_D_Rv => 0
								   } for (0 .. $#segseq);
			# just load the first segment
			# it's ok, always deal with single virus-ref
			last;
		}
		# function 2nd, load each segment as key
		elsif( exists ($Option_Href->{segseq_Href}) ){
			my $SegSeq_Href = $Option_Href->{segseq_Href};
			$SegSeq_Href->{$refseg} = $segseq;
		}
	}
	close REF;
}

#--- read query mapping info from blast m0 result ---
sub read_query_mapinfo_from_blast_m0{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $blast_m0 = $Option_Href->{blast_m0};
	my $query_seq_fa = $Option_Href->{query_seq_fa} || '';
	my $sbjct_seq_fa = $Option_Href->{sbjct_seq_fa};
	my $queryid_list = $Option_Href->{queryid_list} || '';
	my $Query_MapInfo_real_Href = $Option_Href->{query_mapinfo_Href};

	# my $IUPAC_code = 'acgtryswkmndhvn'; # \., -, use later

	my %Required_Query_ID;
	if(length($queryid_list) != 0){
		open (QID,Try_GZ_Read($queryid_list)) || die "fail $queryid_list: $!\n";
		while(<QID>){
			my ($query_id) = (split)[0];
			$Required_Query_ID{$query_id} = 1;
		}
		close QID;
	}

	# read m0 file
	open (M0,Try_GZ_Read($blast_m0)) || die "fail $blast_m0: $!\n";
	while (<M0>){
		if(!/^Query=/){
			next;
		}
		else{
			# result(s) of one query_id starts here
			my ($query_id) = (/Query=\s(\S+)/);
			if($queryid_list && !exists($Required_Query_ID{$query_id})){ # set required query_id
				next;
			}
			# subject mapped by seq of this query_id
			my ($sbjct_id);
			# SIGNAL for control
			## 1, initial, 'wait result'
			## 2, meet Score
			my $signal = 1;
			while($signal){
				my $line_info = <M0>;
				if($line_info !~ /^\s+Score\s+=/){
					# no hits, just skip
					if( $line_info =~ /No\s+hits\s+found/ ){
						last;
					}
					# mapped subject id
					elsif( $line_info =~ /^>(.+)\s$/ ){
						$sbjct_id = $1;
						# reset the signal as 'wait result'.
						$signal = 1;
						next;
					}
					# search next line for certain query result
					elsif( $signal == 1 ){
						next;
					}
					# empty lines
					elsif( $line_info =~ /^\s+$/ ){
						next;
					}
					# other not-allowed info, means finish of this query_id's result
					else{
						last;
					}
				}
				# meet Score:
				$signal = 2;
				# to store the map info
				my $Query_MapInfo_Href = {};
				# instance: Score =  149 bits (75), Expect = 6e-40
				my ($score, $expect_E) = ($line_info =~ /Score\s+=\s+([\d\.]+)\s+bits.*Expect\s+=\s+(\S+)/);
				#die "$line_info\n$score, $expect_E\n"; ## debug
				# instance: Identities = 84/87 (96%)
				my $Identities_line_info = <M0>;
				my ($mapped_length, $identities) = ($Identities_line_info =~ /Identities\s+=\s+\d+\/(\d+)\s+\(([^\(\)]+)\)/);
				#die "$Identities_line_info\n$mapped_length, $identities\n"; ## debug
				($identities = '0.' . $identities) =~ s/\%$//;
				# instance: Strand = Plus / Plus
				my $Strand_line_info = <M0>;
				my ($query_strd, $sbjct_strd) = ($Strand_line_info =~ /Strand\s+=\s+(\S+)\s+\/\s+(\S+)/);
				#die "$Strand_line_info\n$query_strd, $sbjct_strd\n"; ## debug
				$query_strd = ($query_strd =~ /Plus/) ? 1 : -1;
				$sbjct_strd = ($sbjct_strd =~ /Plus/) ? 1 : -1;
				# read query mapping details repeatedly
				my $cat_map_info;
				my @cat_query_info;
				my @cat_sbjct_info;
				while (<M0>){
					if(/^Query:/){ # find Query
						chomp(my $query_info = $_);
						chomp(my $map_info = <M0>);
						chomp(my $sbjct_info = <M0>);
						# split info for details
						my ($prefix, $query_st, $query_seq, $query_ed) = ($query_info =~ /(Query:\s+(\d+)\s+)(\S+)\s+(\d+)/i); # /(Query:\s+(\d+)\s+)([acgt\-n]+)\s+(\d+)/i
						my $pf_len = length($prefix);
						$map_info =~ s/^\s{$pf_len,$pf_len}//;
						my ($sbjct_st, $sbjct_seq, $sbjct_ed) = ($sbjct_info =~ /(\d+)\s+(\S+)\s+(\d+)/i); # /(\d+)\s+([$IUPAC_code\.\-]+)\s+(\d+)/i
						# check
						if(length($query_seq) != length($sbjct_seq) || length($query_seq) != length($map_info)){
							warn_and_exit  "query_id: $query_id\nquery seq length is not equal to sbjct seq length:\n"
										  ."query: $query_seq\n"
										  ."align: $map_info\n"
										  ."sbjct: $sbjct_seq\n"
										  ."$query_info\n"
										  ."$sbjct_info\n";
						}
						# try to concatnate the mapping details
						if(!defined($cat_map_info)){ # first meets
							# just load on
							$cat_map_info = $map_info;
							@cat_query_info = ($query_st, $query_seq, $query_ed);
							@cat_sbjct_info = ($sbjct_st, $sbjct_seq, $sbjct_ed);
						}
						else{
							# smart check the continuous attribute
							if( ($query_st-$cat_query_info[2]) != $query_strd || ($sbjct_st-$cat_sbjct_info[2]) != $sbjct_strd ){
								warn_and_exit  "the mapping is not continuous, Query_ID: $query_id\n"
											  ."query: prev:$cat_query_info[0]-$cat_query_info[2], now:$query_st-$query_ed\n"
									 		  ."sbjct: prev:$cat_sbjct_info[0]-$cat_sbjct_info[2], now:$sbjct_st-$sbjct_ed\n";
							}
							# concatnate
							$cat_map_info      .= $map_info;
							$cat_query_info[1] .= $query_seq;
							$cat_query_info[2]  = $query_ed; # update the end
							$cat_sbjct_info[1] .= $sbjct_seq;
							$cat_sbjct_info[2]  = $sbjct_ed; # update the end
						}
						# check the map_info length, and trigger breaking
						if(length($cat_map_info) == $mapped_length){
							last;
						}
					}
				}
				# specifically, adjust the map-info for 'insertion' after 'mismatch' position
				# only change 'sbjct_seq'
				&adjust_ins_after_mismatch(
					{
						query_seq_Sref => \$cat_query_info[1],
						map_info_Sref  => \$cat_map_info,
						sbjct_seq_Sref => \$cat_sbjct_info[1]
					}
				);
				# store this mapping for this Query
				$Query_MapInfo_Href->{score}      = $score;
				$Query_MapInfo_Href->{expect_E}   = $expect_E;
				$Query_MapInfo_Href->{ident}      = $identities;
				$Query_MapInfo_Href->{query_strd} = $query_strd;
				$Query_MapInfo_Href->{sbjct_strd} = $sbjct_strd;
				$Query_MapInfo_Href->{align_strd} = $sbjct_strd * $query_strd;
				$Query_MapInfo_Href->{sbjct_id}   = $sbjct_id;
				$Query_MapInfo_Href->{map_info}   = $cat_map_info;
				$Query_MapInfo_Href->{query_info} = \@cat_query_info;
				$Query_MapInfo_Href->{sbjct_info} = \@cat_sbjct_info;
				# store
				push @{$Query_MapInfo_real_Href->{$query_id}}, $Query_MapInfo_Href;
				#print "$cat_query_info[1]\n$cat_map_info\n$cat_sbjct_info[1]\n"; # debug
			}
		}
	}
	close M0;

	# Try to concatenate devided seq by large-DEL for one same query_id
	if( $query_seq_fa || -e $query_seq_fa ){
		&concatenate_GapInducedByLargeInDel(
			{
				query_mapinfo_Href => $Query_MapInfo_real_Href,
				query_seq_fa => $query_seq_fa,
				sbjct_seq_fa => $sbjct_seq_fa
			}
		);
	}
}

# specifically, adjust the map-info for 'insertion' after 'mismatch' position
# only change 'sbjct_seq'
sub adjust_ins_after_mismatch{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $querySeq_Sref = $Option_Href->{query_seq_Sref};
	my $mapInfo_Sref  = $Option_Href->{map_info_Sref};
	my $sbjctSeq_Sref = $Option_Href->{sbjct_seq_Sref};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();

	# find the insertion pos
	my $update_sign = 0;
	my @sbjctSeq = split //, $$sbjctSeq_Sref;
	my @querySeq = split //, $$querySeq_Sref;
	for ( my $i = 0; $i <= $#sbjctSeq; $i++ ){
		# insertion follows mismatch
		if(    $sbjctSeq[$i] eq '-' # insertion start here
			&& $sbjctSeq[$i-1] ne '-' # last pos is not insertion, {never be}.
			&& $querySeq[$i-1] ne '-' # last pos is not deletion.
			&& lc($sbjctSeq[$i-1]) ne lc($querySeq[$i-1]) # last pos is mismatch
		){
			# first, check the mimatch extending to 5-prime
			my $mism_i = $i-1;
			while(    $mism_i > 0 # cannot reach the 5-prime end
				   && $sbjctSeq[$mism_i-1] ne '-' # last pos is not insertion, {never be}.
				   && $querySeq[$mism_i-1] ne '-' # last pos is not deletion.
				   && lc($sbjctSeq[$mism_i-1]) ne lc($querySeq[$mism_i-1]) # last one is still mismatch
			){
				$mism_i --;
			}
			# second, check the insertion extending to 3-prime
			my $ins_len = 1; # already is one-bp long
			while(    $i+$ins_len <= $#sbjctSeq
				   && $sbjctSeq[$i+$ins_len] eq '-'
			){
				$ins_len ++;
			}
			# third, use map_info to validate
			my $substr_idx = $mism_i;
			my $substr_len = $i - $mism_i + $ins_len;
			my $thisReg_map_Info = substr($$mapInfo_Sref,  $substr_idx, $substr_len);
			my $thisReg_querySeq = substr($$querySeq_Sref, $substr_idx, $substr_len);
			my $thisReg_sbjctSeq = substr($$sbjctSeq_Sref, $substr_idx, $substr_len);
			## map info of this region should be all blanks
			if( $thisReg_map_Info =~ /\S/ ){
				warn_and_exit `date`."<ERROR>:\tmeet wrong map-info when adjust insertion [idx:$i] after mismatch from Blast result:\n"
									."\tquerySeq: >$thisReg_querySeq<\n"
									."\tmap_info: >$thisReg_map_Info<\n"
									."\tsbjctSeq: >$thisReg_sbjctSeq<\n";
			}
			# fourth, exchange location of insertion and mismatch
			## only need to change 'sbjctSeq'
			my $thisReg_new_sbjctSeq = substr($thisReg_sbjctSeq, $i-$mism_i) . substr($thisReg_sbjctSeq, 0, $i-$mism_i);
			for my $j ( 1 .. $substr_len ){
				$sbjctSeq[$mism_i+$j-1] = substr($thisReg_new_sbjctSeq, $j-1, 1);
			}
			## record update sign
			$update_sign = 1;
			## infrom this exchange in debug mode
			if( $V_Href->{in_debug} ){
				warn  `date`."Adjust insertion [idx:$i] after mismatch from Blast result.\n"
							." From:\n"
							."\tquerySeq: >$thisReg_querySeq<\n"
							."\tmap_info: >$thisReg_map_Info<\n"
							."\tsbjctSeq: >$thisReg_sbjctSeq<\n"
							." To:\n"
							."\tquerySeq: >$thisReg_querySeq<\n"
							."\tmap_info: >$thisReg_map_Info<\n"
							."\tsbjctSeq: >$thisReg_new_sbjctSeq<\n";
			}
			# update the $i
			$i += $ins_len;
		}
	}

	# update the 'sbjct_seq'
	if( $update_sign ){
		$$sbjctSeq_Sref = join('', @sbjctSeq);
	}
}

#--- try to fill the gap induced by large insertion or deletion ---
sub concatenate_GapInducedByLargeInDel{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $Query_MapInfo_real_Href = $Option_Href->{query_mapinfo_Href};
	my $query_seq_fa = $Option_Href->{query_seq_fa};
	my $sbjct_seq_fa = $Option_Href->{sbjct_seq_fa};

	# load query_seq of each query_id
	my %query_seq;
	&load_refpos_allel_from_ref_fa( { orig_fa_ref=>$query_seq_fa, segseq_Href=>\%query_seq } );
	# load sbjct_seq of each refseg
	my %sbjct_seq;
	&load_refpos_allel_from_ref_fa( { orig_fa_ref=>$sbjct_seq_fa, segseq_Href=>\%sbjct_seq } );

	# search for candidate segments to concatenate
	for my $query_id (sort keys %$Query_MapInfo_real_Href){
		my $MapInfo_Aref = $Query_MapInfo_real_Href->{$query_id};
		my %discard_idx;
		# find the connection iteratively until no more found
		while(1){
			my $map_count = scalar(@$MapInfo_Aref);
			# no more mappings to deal with
			if( $map_count - scalar(keys %discard_idx) == 1 ){
				last;
			}
			# update sign
			my $update_sign = 0;
			# loop to select fore mapping
			for my $fore_map_idx ( 0 .. $map_count-1 ){
				# already discard, old mapping
				next if( exists($discard_idx{$fore_map_idx}) );
				# fore mapping information
				my $fore_mapInfo_Href = $MapInfo_Aref->[$fore_map_idx];
				my $fore_map_sbjct = $fore_mapInfo_Href->{sbjct_id};
				my $fore_map_Qstrd = $fore_mapInfo_Href->{query_strd};
				my $fore_map_Sstrd = $fore_mapInfo_Href->{sbjct_strd};
				my $fore_map_QedP  = $fore_mapInfo_Href->{query_info}->[2];
				my $fore_map_SedP  = $fore_mapInfo_Href->{sbjct_info}->[2];
				my $fore_map_SpanL = length $fore_mapInfo_Href->{map_info};
				# loop to test each other mappings as follower mapping
				for my $cand_map_idx ( 0 .. $map_count-1 ){
					# already discard, old mapping
					next if( exists($discard_idx{$cand_map_idx}) );
					# cand mapping information
					my $cand_mapInfo_Href = $MapInfo_Aref->[$cand_map_idx];
					my $cand_map_sbjct = $cand_mapInfo_Href->{sbjct_id};
					my $cand_map_Qstrd = $cand_mapInfo_Href->{query_strd};
					my $cand_map_Sstrd = $cand_mapInfo_Href->{sbjct_strd};
					my $cand_map_QstP  = $cand_mapInfo_Href->{query_info}->[0];
					my $cand_map_SstP  = $cand_mapInfo_Href->{sbjct_info}->[0];
					my $cand_map_SpanL = length $cand_mapInfo_Href->{map_info};
					# filters
					if(    $cand_map_idx == $fore_map_idx # not itself
						|| $fore_map_sbjct ne $cand_map_sbjct # must map to same subject
						|| $fore_map_Qstrd ne $cand_map_Qstrd # query mapping strand
						|| $fore_map_Sstrd ne $cand_map_Sstrd # sbjct mapping strand
					){
						next;
					}
					# the two segment must be continuous at position
					## large DEL: 'none' OR 'slight' gap on query_seq, and 'large' gap on sbjct_seq
					## large INS: 'none' OR 'slight' gap on sbjct_seq, and 'large' gap on query_seq
					## large INS is from 'unkown source', Here
					if(    $fore_map_Qstrd * $fore_map_QedP < $cand_map_Qstrd * $cand_map_QstP
						&& $fore_map_Sstrd * $fore_map_SedP < $cand_map_Sstrd * $cand_map_SstP
						&& (  # query_seq has large deletion
							  (   abs($cand_map_QstP - $fore_map_QedP) <= 10 # empirical distance in practise
							   && abs($cand_map_SstP - $fore_map_SedP) >= 10
							  )
							|| # query_seq has large insertion (unkown source)
							  (   abs($cand_map_QstP - $fore_map_QedP) >= 10
							   && abs($cand_map_SstP - $fore_map_SedP) <= 10
							  )
						   )
					){
						# prepare makeup seq for query_seq gap
						my $makeup_Qidx = min($fore_map_QedP, $cand_map_QstP);
						my $makeup_Qlen = abs($cand_map_QstP - $fore_map_QedP) - 1;
						my $makeup_Qseq = uc( substr($query_seq{$query_id}, $makeup_Qidx, $makeup_Qlen) );
						if( $fore_map_Qstrd == -1 && $makeup_Qlen > 0 ){
							($makeup_Qseq = reverse uc($makeup_Qseq)) =~ tr/ACGT/TGCA/;
						}
						# prepare makeup seq for sbjct_seq gap
						my $makeup_Sidx = min($fore_map_SedP, $cand_map_SstP);
						my $makeup_Slen = abs($cand_map_SstP - $fore_map_SedP) - 1;
						my $makeup_Sseq = uc( substr($sbjct_seq{$fore_map_sbjct}, $makeup_Sidx, $makeup_Slen) );
						if( $fore_map_Sstrd == -1 && $makeup_Slen > 0 ){
							($makeup_Sseq = reverse uc($makeup_Sseq)) =~ tr/ACGT/TGCA/;
						}
						# compare and makeup map-info
						my $diff_len = $makeup_Slen - $makeup_Qlen;
						my $makeup_MapI = join('', (' ' x $diff_len));
						# gap part
						## query_seq has large deletion
						if( $diff_len > 0 ){
							$makeup_Qseq = join('', ('-' x $diff_len)) . $makeup_Qseq;
						}
						## query_seq has large insertion
						else{
							$makeup_Sseq = join('', ('-' x $diff_len)) . $makeup_Sseq;
						}
						# remained part, maybe mismatch(s)
						for my $i ( reverse ( 1 .. min($makeup_Qlen, $makeup_Slen) ) ){
							my $Q_allele = uc( substr($makeup_Qseq, -1*$i, 1) );
							my $S_allele = uc( substr($makeup_Sseq, -1*$i, 1) );
							$makeup_MapI .= (($Q_allele eq $S_allele)?'|':' ');
						}
						# concatenate 'makeup'+'cand' to 'fore'
						## still record in 'fore_mapInfo_Href'
						$fore_mapInfo_Href->{score} += $cand_mapInfo_Href->{score};
						$fore_mapInfo_Href->{ident} = ($fore_mapInfo_Href->{ident}+$cand_mapInfo_Href->{ident}) / 2;
						$fore_mapInfo_Href->{map_info} .= ($makeup_MapI . $cand_mapInfo_Href->{map_info});
						$fore_mapInfo_Href->{query_info}->[1] .= lc($makeup_Qseq . $cand_mapInfo_Href->{query_info}->[1]);
						$fore_mapInfo_Href->{query_info}->[2] = $cand_mapInfo_Href->{query_info}->[2];
						$fore_mapInfo_Href->{sbjct_info}->[1] .= lc($makeup_Sseq . $cand_mapInfo_Href->{sbjct_info}->[1]);
						$fore_mapInfo_Href->{sbjct_info}->[2] = $cand_mapInfo_Href->{sbjct_info}->[2];
						# inform
						stout_and_sterr `date`."[INFO]:\tConnect two blast mapping blocks (large DEL/uks-INS) of query_id($query_id):\n"
											  ."\tfore_mapping interval, idx=$fore_map_idx:\n"
											  ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_map_QedP;"
											  .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_map_SedP\n"
											  ."\ttail_mapping interval, idx=$cand_map_idx:\n"
											  ."\t Qst:$cand_map_QstP; Qed:$cand_mapInfo_Href->{query_info}->[2];"
											  .  " Sst:$cand_map_SstP; Sed:$cand_mapInfo_Href->{sbjct_info}->[2]\n"
											  ."\tnew_mapping interval, idx=$fore_map_idx:\n"
											  ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_mapInfo_Href->{query_info}->[2];"
											  .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_mapInfo_Href->{sbjct_info}->[2]\n"
											  # debug, show the new connected mapping block
											  # ."\t '$fore_mapInfo_Href->{query_info}->[1]'\n"
											  # ."\t '$fore_mapInfo_Href->{map_info}'\n"
											  # ."\t '$fore_mapInfo_Href->{sbjct_info}->[1]'\n"
											  ."\tmakeup_mapInfo:\n"
											  ."\t Q: '$makeup_Qseq'\n"
											  ."\t M: '$makeup_MapI'\n"
											  ."\t S: '$makeup_Sseq'\n";
						# record cand_map_idx as old_discarded
						$discard_idx{$cand_map_idx} = 1;
						# reocrd update sign
						$update_sign = 1;
						# stop to start the whole new detection loop
						last;
					}
					## large INS: 'none' OR 'slight' overlap on query_seq, and 'large' overlap on sbjct_seq
					## large INS is tandem duplicated, Here
					elsif(    $fore_map_Qstrd * $fore_map_QedP >= $cand_map_Qstrd * $cand_map_QstP - 1
						   && $fore_map_Sstrd * $fore_map_SedP >  $cand_map_Sstrd * $cand_map_SstP
						   &&  # query_seq has large insertion (tandem duplicated)
							  (   abs($cand_map_QstP - $fore_map_QedP) <= 10 # empirical distance in practise
							   && abs($cand_map_SstP - $fore_map_SedP) >= 10
							   && abs($cand_map_SstP - $fore_map_SedP) <= 50 # it is the max ins-size detected from traditional PE seq
							   && abs($cand_map_SstP - $fore_map_SedP) <= $fore_map_SpanL / 5 # fore maps long enough, empirical value
							   && abs($cand_map_SstP - $fore_map_SedP) <= $cand_map_SpanL / 5 # cand maps long enough, empirical value
							  )
					){
						# 1st, discard possible extended end of cand query_seq to make it continuous
						my $QextLen = ( $fore_map_QedP - $cand_map_QstP ) * $fore_map_Qstrd + 1;
						# modify the candidate info
						my @cand_Qseq = split //, $cand_mapInfo_Href->{query_info}->[1];
						my $cand_Qidx = -1;
						while( $QextLen > 0 ){
							$cand_Qidx ++;
							if( $cand_Qseq[$cand_Qidx] ne '-' ){
								$QextLen --;
							}
						}
						if( $cand_Qidx >= 0 ){
							my $DiscLen = $cand_Qidx + 1;
							# map info
							$cand_mapInfo_Href->{map_info} = substr( $cand_mapInfo_Href->{map_info}, $DiscLen );
							# query_seq
							$cand_mapInfo_Href->{query_info}->[1] = substr( $cand_mapInfo_Href->{query_info}->[1], $DiscLen );
							# sbjct_seq
							my $Slen = grep !/-/, ( split //, substr( $cand_mapInfo_Href->{sbjct_info}->[1], 0, $DiscLen ) );
							$cand_mapInfo_Href->{sbjct_info}->[1] = substr( $cand_mapInfo_Href->{sbjct_info}->[1], $DiscLen );
							$cand_map_SstP += $Slen * $fore_map_Sstrd;
						}
						# 2nd, introduce insertion to make subject_seq continuous
						my $SextLen = ( $fore_map_SedP - $cand_map_SstP ) * $fore_map_Sstrd + 1;
						# modify the candidate info
						my @cand_Sseq = split //, $cand_mapInfo_Href->{sbjct_info}->[1];
						my $cand_Sidx = -1;
						while( $SextLen > 0 ){
							$cand_Sidx ++;
							if( $cand_Sseq[$cand_Sidx] ne '-' ){
								$SextLen --;
							}
						}
						if( $cand_Sidx >= 0 ){ # must be
							my $ReplcLen = $cand_Sidx + 1;
							# map info
							$cand_mapInfo_Href->{map_info}        = join('', (' ' x $ReplcLen)) . substr($cand_mapInfo_Href->{map_info},        $ReplcLen);
							# sbjct_seq
							$cand_mapInfo_Href->{sbjct_info}->[1] = join('', ('-' x $ReplcLen)) . substr($cand_mapInfo_Href->{sbjct_info}->[1], $ReplcLen);
						}
						# concatenate 'makeup'+'cand' to 'fore'
						## still record in 'fore_mapInfo_Href'
						$fore_mapInfo_Href->{score} += $cand_mapInfo_Href->{score};
						$fore_mapInfo_Href->{ident} = ($fore_mapInfo_Href->{ident}+$cand_mapInfo_Href->{ident}) / 2;
						$fore_mapInfo_Href->{map_info} .= $cand_mapInfo_Href->{map_info};
						$fore_mapInfo_Href->{query_info}->[1] .= $cand_mapInfo_Href->{query_info}->[1];
						$fore_mapInfo_Href->{query_info}->[2]  = $cand_mapInfo_Href->{query_info}->[2];
						$fore_mapInfo_Href->{sbjct_info}->[1] .= $cand_mapInfo_Href->{sbjct_info}->[1];
						$fore_mapInfo_Href->{sbjct_info}->[2]  = $cand_mapInfo_Href->{sbjct_info}->[2];
						# inform
						stout_and_sterr `date`."[INFO]:\tConnect two blast mapping blocks (large tdp-INS) of query_id($query_id):\n"
											  ."\tfore_mapping interval, idx=$fore_map_idx:\n"
											  ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_map_QedP;"
											  .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_map_SedP\n"
											  ."\ttail_mapping interval, idx=$cand_map_idx:\n"
											  ."\t Qst:$cand_mapInfo_Href->{query_info}->[0]; Qed:$cand_mapInfo_Href->{query_info}->[2];"
											  .  " Sst:$cand_mapInfo_Href->{sbjct_info}->[0]; Sed:$cand_mapInfo_Href->{sbjct_info}->[2]\n"
											  ."\tnew_mapping interval, idx=$fore_map_idx:\n"
											  ."\t Qst:$fore_mapInfo_Href->{query_info}->[0]; Qed:$fore_mapInfo_Href->{query_info}->[2];"
											  .  " Sst:$fore_mapInfo_Href->{sbjct_info}->[0]; Sed:$fore_mapInfo_Href->{sbjct_info}->[2]\n";
											  # debug, show the new connected mapping block
											  # ."\t '$fore_mapInfo_Href->{query_info}->[1]'\n"
											  # ."\t '$fore_mapInfo_Href->{map_info}'\n"
											  # ."\t '$fore_mapInfo_Href->{sbjct_info}->[1]'\n"
						# record cand_map_idx as old_discarded
						$discard_idx{$cand_map_idx} = 1;
						# reocrd update sign
						$update_sign = 1;
						# stop to start the whole new detection loop
						last;
					}
				}
				# stop to start the whole new detection loop
				if( $update_sign ){
					last;
				}
			}
			# stop when no update operation is done
			unless( $update_sign ){
				last;
			}
		}
		# discard any recorded idx
		if( scalar(keys %discard_idx) != 0 ){
			my @keep_idx;
			for my $idx ( 0 .. scalar(@$MapInfo_Aref)-1 ){
				unless( exists($discard_idx{$idx}) ){
					push @keep_idx, $MapInfo_Aref->[$idx];
				}
			}
			# re-assign new mapping result array-ref to this query_id
			$Query_MapInfo_real_Href->{$query_id} = \@keep_idx;
		}
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
