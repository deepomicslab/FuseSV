package FuseSV::Virus_Type::BlastM0Mutation;

use strict;
use warnings;
use List::Util qw[min max sum];
use SOAPfuse::General_Operation qw/warn_and_exit/;
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  read_blast_m0_for_mutation
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BlastM0Mutation';
#----- version --------
$VERSION = "0.02";
$DATE = '2016-04-12';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						read_blast_m0_for_mutation
						upload_mutation_list
						summary_alignment_result
						load_refpos_allel_from_ref_fa
						read_query_mapinfo_from_blast_m0
					 /;

#--- read blast m0 for mutation ---
sub read_blast_m0_for_mutation{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $blast_m0 = $Option_Href->{blast_m0};
	my $queryid_list = exists($Option_Href->{queryid_list}) ? $Option_Href->{queryid_list} : '';
	my $orig_fa_ref = $Option_Href->{orig_fa_ref};
	my $min_depth = $Option_Href->{min_depth} || 1;
	my $both_strd = $Option_Href->{both_strd} || 0;
	my $min_mut_reads = $Option_Href->{min_mut_reads} || 2;
	my $mutation_Href = $Option_Href->{mutation_Href};

	my %Query_MapInfo;
	&read_query_mapinfo_from_blast_m0(
			{
				blast_m0 => $blast_m0,
				queryid_list => $queryid_list,
				query_mapinfo_Href => \%Query_MapInfo
			}
		);

	my %RefPos;
	&load_refpos_allel_from_ref_fa(
			{
				orig_fa_ref => $orig_fa_ref,
				refpos_Href => \%RefPos
			}
		);

	&summary_alignment_result(
			{
				query_mapinfo_Href => \%Query_MapInfo,
				refpos_Href => \%RefPos
			}
		);

	&upload_mutation_list(
			{
				mutation_Href => $mutation_Href,
				refpos_Href => \%RefPos,
				min_depth => $min_depth,
				both_strd => $both_strd,
				min_mut_reads => $min_mut_reads
			}
		);
}
#--- define mutation or not for each base, and output ---
sub upload_mutation_list{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $RefPos_Href = $Option_Href->{refpos_Href};
	my $major_mutation_Href = $Option_Href->{mutation_Href};
	my $min_depth = $Option_Href->{min_depth};
	my $both_strd = $Option_Href->{both_strd};
	my $min_mut_reads = $Option_Href->{min_mut_reads};

	my $sbjct_filter_pos = 1;

	# mutations
	for my $sbjct_pos (sort {$a<=>$b} keys %$RefPos_Href){
		my $ref_allel = $RefPos_Href->{$sbjct_pos}->{ref_allel};
		my $depth = $RefPos_Href->{$sbjct_pos}->{depth} || 0;
		my $allel_depth = $RefPos_Href->{$sbjct_pos}->{refallel_depth} || 0;
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
				##warn "Part:3\n"; # debug
				my $mutation_Href = $RefPos_Href->{$sbjct_pos}->{mutation};
				# select the mutations
				my @mutation_id;
				for my $mutation_id (keys %$mutation_Href){
					##warn "Part:3-1\n"; # debug
					my ($sum_sup, $plus_sup, $minus_sup) = @{ $mutation_Href->{$mutation_id} };
					# filter the mutation
					unless(    $sum_sup < $min_mut_reads
					  	    || ( $both_strd && (!$plus_sup || !$minus_sup) )
					  	  )
					{
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
					my ($mut_type,$mut_seq) = split /,/, $mutation_id[0]; # first one
					# record mutation
					$major_mutation_Href->{$sbjct_pos} = [$mut_type, $mut_seq];
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

	for my $query_id (keys %$Query_MapInfo_real_Href){
		for my $Query_MapInfo_Href (@{ $Query_MapInfo_real_Href->{$query_id} }){
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
						warn_and_exit("Should be a Mismatch: $sbjct_pos\nQ: $query_allel\n    \nS: $sbjct_allel\n");
					}
					# accumulate depth
					$RefPos_Href->{$sbjct_pos}->{depth} ++;
					# record mutation
					my $mutation_id = "snp,$query_allel";
					# format (Sum, Plus, Minus)
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0];
					my $mutation_info_Aref = $RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id};
					$mutation_info_Aref->[0] ++;
					if($align_strd > 0){
						$mutation_info_Aref->[1] ++;
					}
					else{
						$mutation_info_Aref->[2] ++;
					}
					##warn "@{$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id}}\n"; # debug
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
					# do not need to accumulate depth for this pos
					#$RefPos_Href->{$sbjct_pos}->{depth} ++;
					# record mutation
					my $mutation_id = "ins,$insert_seq";
					# format (Sum, Plus, Minus)
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0];
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
					# accumulate depth
					$RefPos_Href->{$_}->{depth} ++ for ($sbjct_pos .. $sbjct_pos+$j);
					# record mutation
					my $mutation_id = "del,$delete_seq";
					# format (Sum, Plus, Minus)
					$RefPos_Href->{$sbjct_pos}->{mutation}->{$mutation_id} ||= [0,0,0];
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
sub load_refpos_allel_from_ref_fa{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $orig_fa_ref = $Option_Href->{orig_fa_ref};
	my $RefPos_Href = $Option_Href->{refpos_Href};

	open (REF, $orig_fa_ref) || die "fail $orig_fa_ref: $!\n";
	$/ = '>'; <REF>; $/ = "\n";
	while(<REF>){
		chomp(my $refseg = $_);
		$/ = '>';
		chomp(my $segseq = <REF>);
		$/ = "\n";
		$segseq =~ s/\s//g;
		my @segseq = split //, $segseq;
		$RefPos_Href->{$_+1}->{ref_allel} = $segseq[$_] for (0 .. $#segseq);
	}
	close REF;
}

#--- read query mapping info from blast m0 result ---
sub read_query_mapinfo_from_blast_m0{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];

	my $blast_m0 = $Option_Href->{blast_m0};
	my $queryid_list = exists($Option_Href->{queryid_list}) ? $Option_Href->{queryid_list} : '';
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
			my ($query_id) = (/Query=\s(\S+)/);
			if($queryid_list && !exists($Required_Query_ID{$query_id})){ # set required query_id
				next;
			}
			my $signal = 1;
			while($signal){
				my $line_info = <M0>;
				if($line_info !~ /^\s+Score\s+=/){
					if($line_info =~ /No\s+hits\s+found/){last;} # no hits, just skip
					elsif($signal == 1) {next;}
					elsif($line_info =~ /^\s+$/){next;}
					else{last;}
				}
				# meet Score:
				$signal = 2;
				# to store the map info
				my $Query_MapInfo_Href = {};
				# instance: Score =  149 bits (75), Expect = 6e-40
				my ($score, $expect_E) = ($line_info =~ /Score\s+=\s+(\d+)\s+bits.*Expect\s+=\s+(\S+)/);
				#die "$line_info\n$score, $expect_E\n"; ## debug
				# instance: Identities = 84/87 (96%)
				my $Identities_line_info = <M0>;
				my ($mapped_length, $identities) = ($Identities_line_info =~ /Identities\s+=\s+\d+\/(\d+)\s+\(([^\(\)]+)\)/);
				#die "$Identities_line_info\n$mapped_length, $identities\n"; ## debug
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
							warn "query_id: $query_id\nquery seq length is not equal to sbjct seq length:\n" .
								 "query: $query_seq\n" .
								 "align: $map_info\n" .
								 "sbjct: $sbjct_seq\n" .
								 "$query_info\n" .
								 "$sbjct_info\n";
							exit(1);
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
								warn "the mapping is not continuous, Query_ID: $query_id\n" .
									 "query: prev:$cat_query_info[0]-$cat_query_info[2], now:$query_st-$query_ed\n" .
									 "sbjct: prev:$cat_sbjct_info[0]-$cat_sbjct_info[2], now:$sbjct_st-$sbjct_ed\n";
								exit(1);
							}
							# concatnate
							$cat_map_info .= $map_info;
							$cat_query_info[1] .= $query_seq;
							$cat_query_info[2]  = $query_ed;
							$cat_sbjct_info[1] .= $sbjct_seq;
							$cat_sbjct_info[2]  = $sbjct_ed;
						}
						# check the map_info length, and trigger breaking
						if(length($cat_map_info) == $mapped_length){
							last;
						}
					}
				}
				# store this mapping for this Query
				$Query_MapInfo_Href->{score}      = $score;
				$Query_MapInfo_Href->{expect_E}   = $expect_E;
				$Query_MapInfo_Href->{ident}      = $identities;
				$Query_MapInfo_Href->{query_strd} = $query_strd;
				$Query_MapInfo_Href->{sbjct_strd} = $sbjct_strd;
				$Query_MapInfo_Href->{align_strd} = $sbjct_strd * $query_strd;
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
}

#--- 
1; ## tell the perl script the successful access of this module.
