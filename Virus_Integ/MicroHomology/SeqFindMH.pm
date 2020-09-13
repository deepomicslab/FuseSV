package FuseSV::Virus_Integ::MicroHomology::SeqFindMH;

use strict;
use warnings;
use List::Util qw/ min max sum /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
	          Find_vtMH
	          Find_sfMH
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'SeqFindMH';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-08-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						Find_vtMH
						Find_sfMH
					 /;

#--- find vertical MH for given seq pair ---
sub Find_vtMH{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $seq_a = uc( $parm{seq_a} );
	my $seq_b = uc( $parm{seq_b} );
	my $min_MHlen = $parm{min_MHlen} || 2;
	my $max_MHlen = $parm{max_MHlen} || 1000;

	my @seq_a = split //, $seq_a;
	my @seq_b = split //, $seq_b;

	my %vt_MH;
	my $max_i = min($#seq_a, $#seq_b);
	for (my $i = 0; $i <= $max_i; ){
		if( $seq_a[$i] eq $seq_b[$i] ){
			my $j = $i + 1;
			while(   $j <= $#seq_a
				  && $seq_a[$j] eq $seq_b[$j]
			){
				$j++;
			}
			my $len = $j - $i;
			if(    $len >= $min_MHlen
				&& $len <= $max_MHlen
			){
				my $extSign = 0;
				if( $i == 0 ){
					$extSign = -1; # might extend to up-stream
				}
				elsif( $j-1 == $#seq_a ){
					$extSign = 1; # might extend to down-stream
				}
				# record
				$vt_MH{$i} = { len=>$len, seq=>substr($seq_a,$i,$len), extSign=>$extSign };
			}
			$i = $j+1;
		}
		else{
			$i++;
		}
	}

	return \%vt_MH;
}

#--- find shift MH for given seq pair ---
sub Find_sfMH{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $seq_a = uc( $parm{seq_a} );
	my $seq_b = uc( $parm{seq_b} );
	my $min_sfMHlen = $parm{min_MHlen} || 2;
	my $max_sfMHlen = $parm{max_MHlen} || 1000;
	my $min_sfMHopl = $parm{min_MHopl} || 0;
	my $max_shift = $parm{max_shift} || 1000;
	my $vtMH_Href = $parm{vtMH_Href} || {};

	my @seq_a = split //, $seq_a;
	my @seq_b = split //, $seq_b;

	my %sf_MH;
	my $last_vtMHlen = 0;
	my $last_vtMH_ed = -1;
	my $last_sfMH_ed_a = -1;
	my $last_sfMH_st_b = -1;
	my $last_sfMH_ed_b = -1;
	for (my $i = 0; $i <= $#seq_a; ){
		my $sum_vtMHlen = ($i > $last_vtMH_ed ? 0 : $last_vtMHlen);
		my $i_add = 1;
		for (my $j = $i; $j <= $#seq_a; $j++){
			if( exists $vtMH_Href->{$j} ){
				$sum_vtMHlen += $vtMH_Href->{$j}->{len};
				if( $j == $i ){ # new vtMH, update!
					$last_vtMHlen = $vtMH_Href->{$j}->{len};
					$last_vtMH_ed = $j + $vtMH_Href->{$j}->{len} - 1;
				}
			}
			my $len = $j - $i + 1;
			next if( $len <= $sum_vtMHlen );
			next if( $len < $min_sfMHlen );
			next if( $len <= $min_sfMHopl );
			last if( $len > $max_sfMHlen );
			# MH seq to scan
			my $itv_ex = $len - $min_sfMHopl;
			my $itv_st = max(0, $i - $itv_ex);
			my $itv_ed = min($#seq_b, $j + $itv_ex);
			my $q_seq = join('',@seq_a[$i..$j]);
			my $s_seq = join('',@seq_b[$itv_st..$itv_ed]);
			# find sfMH with minimum abs(shift)
			my $mShift = $max_shift + 1;
			my $cutLen = 0;
			while( $s_seq =~ /$q_seq/ ){
				my $pre_len = length($`);
				my $shift = $itv_st + $pre_len + $cutLen - $i;
				if(    $shift != 0
					&& abs($shift) <= $max_shift
				){
					$mShift = abs($shift) < abs($mShift) ? $shift : $mShift;
				}
				# destroy the mapped part: just cut the first mapped base
				$s_seq = substr($s_seq, $pre_len+1);
				$cutLen += $pre_len + 1;
			}
			if(    $mShift != $max_shift+1 # found sfMH
				&& (   $last_sfMH_ed_a < $j
					|| $last_sfMH_st_b > $i + $mShift
					|| $last_sfMH_ed_b < $i + $mShift + $len - 1
				   )
			){
				my $extSign = 0;
				if( $i == 0 ){
					$extSign = -5; # might extend to P5 up-stream
				}
				elsif( $i + $mShift == 0 ){
					$extSign = -3; # might extend to P3 up-stream
				}
				elsif( $j == $#seq_a ){
					$extSign = 5; # might extend to P5 down-stream
				}
				elsif( $i + $mShift + $len - 1 == $#seq_b ){
					$extSign = 3; # might extend to P3 down-stream
				}
				# record
				$sf_MH{$i} = { len=>$len, sft=>$mShift, seq=>substr($seq_a, $i, $len), opl=>$len-abs($mShift), extSign=>$extSign };
				# update
				$last_sfMH_ed_a = $j;
				$last_sfMH_st_b = $i + $mShift;
				$last_sfMH_ed_b = $i + $mShift + $len - 1;
				# $i_add = $len;
			}
		}
		$i += $i_add;
	}

	return \%sf_MH;
}

#---
1; ## tell the perl script the successful access of this module.
