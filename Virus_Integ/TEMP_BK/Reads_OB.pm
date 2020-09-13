package FuseSV::Virus_Integ::Reads_OB;

use strict;
use warnings;
use List::Util qw/max/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'Reads_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2017-12-21';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						new
						get_available_rgOB
						update_rgOB_maxRlen
						is_fw_map
						is_rv_map
						is_unmap
						is_good_cigar
						is_softclip
						fore_softclip_length
						hind_softclip_length
						tlen_FixTlen_with_S
						extract_FS_readid_prefix
					 /;

#--- construction of object ---
sub new{
	my $type = shift;
	my %parm = @_;

	chomp($parm{ReadsLineText});
	my @term = split /\t+/, $parm{ReadsLineText};

	my $reads_OB = {};
	# basic infos
	$reads_OB->{pid}   = $term[0];
	$reads_OB->{flag}  = $term[1];
	$reads_OB->{mseg}  = $term[2];
	$reads_OB->{mpos}  = $term[3];
	$reads_OB->{mapQ}  = $term[4];
	$reads_OB->{cigar} = $term[5];
	# paired-end info
	$reads_OB->{p_mseg} = $term[6];
	$reads_OB->{p_mpos} = $term[7];
	$reads_OB->{tlen}   = $term[8];
	# read string
	$reads_OB->{rlen}   = length($term[9]); # if not Hard-clip
	if( $parm{_rc_rdstr} ){
		$reads_OB->{rseq}  = $term[9];
		$reads_OB->{baseQ} = $term[10];
	}
	# optional fields
	if( $parm{_rc_optfd} ){
		$reads_OB->{optfd} = join("\t", @term[11 .. $#term]);
	}
	# reads group
	$reads_OB->{rg_id} = ($parm{ReadsLineText} =~ /RG\:Z\:(\S+)/) ? $1 : undef;
	# reads end NO.
	$reads_OB->{endNO} = ($reads_OB->{flag} & 0x40) ? 1 : 2; # 0x80 is '2'

	bless($reads_OB);
	return $reads_OB;
}

#--- find rg_OB based on rg_id ---
sub get_available_rgOB{

	my $reads_OB = shift;
	my %parm = @_;
	my $bam_OB = $parm{bam_OB};
	my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};

	if(    !defined $reads_OB->{rg_id}
		|| !exists $rgid2rgOB_Href->{$reads_OB->{rg_id}}
	  ){
		warn_and_exit "<ERROR>\tCan not get available RG_ID from reads $reads_OB->{pid}/$reads_OB->{endNO} of $bam_OB->{tissue} bam:\n"
							."\t$bam_OB->{filepath}\n";
	}
	else{
		$reads_OB->{rg_OB} = $rgid2rgOB_Href->{$reads_OB->{rg_id}};
		return $rgid2rgOB_Href->{$reads_OB->{rg_id}};
	}
}

#--- update the maximum read length of its rg_OB ---
sub update_rgOB_maxRlen{
	my $reads_OB = shift;
	my $rg_OB = $reads_OB->{rg_OB};
	my $endNO = $reads_OB->{endNO};
	$rg_OB->{read_Len}->{$endNO} = max( $rg_OB->{read_Len}->{$endNO}, $reads_OB->{rlen} );
}

#--- test whether it is forward-mapped ---
sub is_fw_map{
	my $reads_OB = shift;
	return !($reads_OB->{flag} & 0x4) && !($reads_OB->{flag} & 0x10);
}

#--- test whether it is reversed-mapped ---
sub is_rv_map{
	my $reads_OB = shift;
	return !($reads_OB->{flag} & 0x4) && ($reads_OB->{flag} & 0x10);
}

#--- test whether it is un-mapped ---
sub is_unmap{
	my $reads_OB = shift;
	return $reads_OB->{flag} & 0x4;
}

#--- test whether it has good cigar ---
# good: lacks 'N', 'H' tag
# 'N' may affect RNAseq-alignment
sub is_good_cigar{
	my $reads_OB = shift;
	return $reads_OB->{cigar} !~ /[NH]/;
}

#--- test whether it is soft-clip ---
sub is_softclip{
	my $reads_OB = shift;
	return $reads_OB->{cigar} =~ /S/;
}

#--- get fore-softclip length if possible ---
sub fore_softclip_length{
	my $reads_OB = shift;
	return ($reads_OB->{cigar} =~ /^(\d+)S/) ? $1 : 0;
}

#--- get hind-softclip length if possible ---
sub hind_softclip_length{
	my $reads_OB = shift;
	return ($reads_OB->{cigar} =~ /(\d+)S$/) ? $1 : 0;
}

#--- get the length to fix Tlen based on soft-clip Edge part ---
sub tlen_FixTlen_with_S{
	my $reads_OB = shift;

	if( $reads_OB->is_fw_map ){
		return $reads_OB->fore_softclip_length;
	}
	elsif( $reads_OB->is_rv_map ){
		return $reads_OB->hind_softclip_length;
	}
}

#--- extract and record new read-id FS prefix ---
sub extract_FS_readid_prefix{

	my $reads_OB = shift;
	my %parm = @_;
	my $rg_OB = $parm{rg_OB};

	# variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	my $Tool_Tag = $V_Href->{Tool_Tag};

	if( $reads_OB->{pid} =~ /^(.*_${Tool_Tag}rgNO(\d+)${Tool_Tag}_)(.+)$/ ){
		$rg_OB->{rID_prefix_prev} = $1;
		$rg_OB->{RG_NO} = $2;
		# add tissue tag
		my $curr_tissue_tag = ucfirst( $rg_OB->{tissue} );
		my $prev_tissue_tag = '';
		my $rg_info_prefix = '';
		if( $rg_OB->{rID_prefix_prev} =~ /_${Tool_Tag}ts(.+)${Tool_Tag}(.+)$/ ){
			$prev_tissue_tag = $1;
			$rg_info_prefix = $2;
		}
		else{ # no tissue tag added before
			$rg_info_prefix = $rg_OB->{rID_prefix_prev};
		}
		# construct current read-id prefix
		$rg_OB->{rID_prefix} = "_${Tool_Tag}ts$curr_tissue_tag$prev_tissue_tag${Tool_Tag}$rg_info_prefix";
	}
	else{
		warn_and_exit "<ERROR>\tUnrecognized read-id prefix ($reads_OB->{pid}) from $rg_OB->{tissue} bam file:\n"
							."\t$rg_OB->{bam_OB}->{filepath}\n";
	}

	# inform
	stout_and_sterr "[INFO]\tExtract and update reads-id prefix of reads group $rg_OB->{RG_ID}:\n"
						 ."\tfrom $rg_OB->{rID_prefix_prev} to $rg_OB->{rID_prefix}\n";
}

1; ## tell the perl script the successful access of this module.
