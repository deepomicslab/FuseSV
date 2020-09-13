package FuseSV::Virus_Integ::FuseSVFigure::FuseSVFigure;

use strict;
use warnings;
use File::Spec qw/catfile/;
use List::Util qw[min max sum];
use SOAPfuse::OpenFile qw/Try_GZ_Read Try_GZ_Write/;
use SOAPfuse::General_Operation qw/warn_and_exit stout_and_sterr trible_run_for_success merge_genome_region/;
use FuseSV::LoadOn;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  stat_depth_from_bam_files
			  read_samtools_depth_out
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'FuseSVFigure';
#----- version --------
$VERSION = "0.02";
$DATE = '2016-07-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						stat_depth_from_bam_files
						read_samtools_depth_out
					 /;

#-------- stat the depth of each gene pos from samtools-depth command output ---
sub read_samtools_depth_out{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $Depth_Stat_Href = $Option_Href->{Depth_Stat_Href};
	my $Pos_Depth_Href = $Option_Href->{Pos_Depth_Href};
	my $region_chr = $Option_Href->{region_chr};
	my $region_st = $Option_Href->{region_st};
	my $region_ed = $Option_Href->{region_ed};
	# about normalize the depth
	my $std_depth = $Option_Href->{std_depth};
	my $T_Average_Depth = $Option_Href->{T_Average_Depth};
	my $N_Average_Depth = $Option_Href->{N_Average_Depth};
	# this chr may be circular virus
	# for viral loop genome depth
	my $VLoop_ExtendedLen = $Option_Href->{VLoop_ExtendedLen};
	my $region_chr_EndPos = $Option_Href->{region_chr_EndPos};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	# only for normal and tumor compare mode
	my %depth_amp;
	$depth_amp{T} = ($std_depth)?($std_depth/$T_Average_Depth):1;
	$depth_amp{N} = ($std_depth)?($std_depth/$N_Average_Depth):1;

	# read the depth stat
	for my $tissue (keys %$Depth_Stat_Href){
		# initialize depth of each pos
		$Pos_Depth_Href->{$_}->{$tissue} ||= 0 for ($region_st .. $region_ed);
		# read depth file
		my $depth_stat_out_file = $Depth_Stat_Href->{$tissue};
		open (DEPTH,Try_GZ_Read($depth_stat_out_file)) || die "fail $depth_stat_out_file: $!\n";
		while(<DEPTH>){
			my @ele = split;
			my $segname = shift @ele;
			next if($segname ne $region_chr);
			my $pos = shift @ele;
			next unless($pos>=$region_st && $pos<=$region_ed);
			# may be from several bam files, so use 'sum'
			my $sum_depth = sum(@ele) * $depth_amp{$tissue};
			# may be several lines to one position, so use '+='
			$Pos_Depth_Href->{$pos}->{$tissue} += $sum_depth;
		}
		close DEPTH;
		# for viral loop genome depth, maybe
		if(    $VLoop_ExtendedLen
			&& (   $region_st <= $VLoop_ExtendedLen
				|| $region_ed > $region_chr_EndPos
				)
			){
			for my $ori_pos (1 .. $VLoop_ExtendedLen){
				my $ext_pos = $ori_pos + $region_chr_EndPos;
				my ($ori_pos_depth, $ext_pos_depth) = (0, 0);
				if( exists($Pos_Depth_Href->{$ori_pos}) ){
					$ori_pos_depth = $Pos_Depth_Href->{$ori_pos}->{$tissue};
				}
				if( exists($Pos_Depth_Href->{$ext_pos}) ){
					$ext_pos_depth = $Pos_Depth_Href->{$ext_pos}->{$tissue};
				}
				# original pos depth
				if( exists($Pos_Depth_Href->{$ori_pos}) ){
					$Pos_Depth_Href->{$ori_pos}->{$tissue} += $ext_pos_depth;
				}
				# extended pos depth
				if( exists($Pos_Depth_Href->{$ext_pos}) ){
					$Pos_Depth_Href->{$ext_pos}->{$tissue} += $ori_pos_depth;
				}
			}
		}
	}
}

#--- read bam for insertsize of properly mapped PE reads ---
sub stat_depth_from_bam_files{

	my $Option_Href = ($_[0]=~/::$MODULE_NAME/) ? $_[1] : $_[0];
	my $bam_Href = $Option_Href->{bam_Href};
	my $source_bam_Href = $Option_Href->{source_bam_Href};
	my $bed_file = $Option_Href->{bed_file};
	my $samtools = $Option_Href->{samtools};
	my $Depth_Stat_Href = $Option_Href->{Depth_Stat_Href};

	# system variants
	my $V_Href = FuseSV::LoadOn->load_variants_dict();
	# my $xxx = $V_Href->{xxx};

	# view region para
	my $view_region_para = '';
	chomp(my $region_count = `cat $bed_file | wc -l`);
	if($region_count <= 20){
		my %View_Region;
		my $Region_NO = 1;
		open (REG, Try_GZ_Read($bed_file)) || die "fail to read $bed_file: $!\n";
		while(<REG>){
			my ($view_chr, $reg_st, $reg_ed) = (split)[0, 1, 2];
			$reg_st ++;
			$View_Region{$view_chr}{$Region_NO++} = [$reg_st, $reg_ed];
		}
		close REG;

		SOAPfuse::General_Operation->merge_genome_region(
				{
					Region_Href => \%View_Region,
					Next_RegNO => $Region_NO
				}
			);

		for my $ref_seg (sort keys %View_Region){
			$view_region_para .= "'$ref_seg:".
					  $View_Region{$ref_seg}{$_}->[0]."-".
					  $View_Region{$ref_seg}{$_}->[1]."' "
					    for sort {$a<=>$b} keys %{ $View_Region{$ref_seg} };
		}
	}
	else{
		$view_region_para = "-L $bed_file";
	}

	# gather depth info from bam files
	my @bam_Href_array = ($bam_Href, $source_bam_Href);
	for my $tissue (sort keys %$bam_Href){
		my $depth_stat_out_file = $Depth_Stat_Href->{$tissue};
		unlink $depth_stat_out_file; # delete files
		for my $i (0, 1){
			my $useful_bam_Href = $bam_Href_array[$i];
			next unless( exists($useful_bam_Href->{$tissue}) );
			my @bam = sort keys %{$useful_bam_Href->{$tissue}};
			# samtools depth
			for my $bam (@bam){
				my $tmp_bam = "$depth_stat_out_file.$tissue.tmp.bam";
				my $command;
				$command = "($samtools view -h -F 0xD00 ";
				if($view_region_para =~ /^\-L/){
					$command .= "$view_region_para $bam | ";
				}
				else{
					$command .= "$bam $view_region_para | ";
				}
				$command .= 'perl -ne \'if(/^@/){print;} else{@ele=split; $rlen=length($ele[9]); $cigar=$ele[5]; $Mlen=0; while($cigar=~/(\d+)M/){$Mlen+=$1; $cigar=~s/$1M//;} if($Mlen/$rlen>=0.33){print;}}\' | '; # if($i == 0);
				$command .= 'awk \'/^@/ || $3~/-V[0-9]$/\' | ' if($i == 0); # only stat for the virus segments.
				$command .= "$samtools view -b -S -o $tmp_bam -) && ";
				$command .= "($samtools sort -T $tmp_bam.sort_tmp -o $tmp_bam.sort.bam $tmp_bam) && ";
				$command .= "($samtools index $tmp_bam.sort.bam) && ";
				$command .= "($samtools depth -b $bed_file -d 100000 $tmp_bam.sort.bam | gzip -c >> $depth_stat_out_file)";
				trible_run_for_success($command, 'depth_from_bam', {esdo_Nvb=>1});
				unlink $tmp_bam, "$tmp_bam.sort.bam", "$tmp_bam.sort.bam.bai";
			}
		}
		# inform
		stout_and_sterr `date`."Stat gene-region depth of $tissue Bams ok.\n";
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
