package FuseSV::Virus_Integ::LocalHaplotypeCircos::Features;

use strict;
use warnings;
use BioFuse::Util::Log qw/warn_and_exit stout_and_sterr/;
use FuseSV::LoadOn;
use FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
				load_features_source_type
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
								 OTHER   => [qw()]);

$MODULE_NAME = 'Features';
#----- version --------
$VERSION = "0.03";
$DATE = '2020-07-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						 load_features_source_type
						 load_features_list
						 draw_feature_icon_on_seg
						/;

#--- load features source types information ---
sub load_features_source_type{

	# variants
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};
	my $Feature_source = $V_Href->{Feature_source};

	# read feature list setting
	open (FTS,$Feature_source) || die "fail open $Feature_source: $!\n";
	# theme
	# theme_ID \t type \t show_text \t short_name \t sn_index
	(my $theme_line = lc(<FTS>)) =~ s/^#//;
	my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
	while (<FTS>){
		next if(/^\#/);
		chomp;
		my @info = split /\t+/; # use '\t+' to split.
		my %ft_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
		# check the existence and duplications
		## feature_id
		if( exists($Feature_Basic_Href->{feature}->{ $ft_Ob{feature_id} } ) ){
			stout_and_sterr "<WARN>:\tfeature_id $ft_Ob{feature_id} is duplicated.\n$_\n";
			next;
		}
		## feature short_name
		for my $ft_id (keys %{$Feature_Basic_Href->{feature}}){
			if( $Feature_Basic_Href->{feature}->{$ft_id}->{short_name} eq $ft_Ob{short_name} ){
				warn_and_exit "<WARN>:\tfeature short_name $ft_Ob{short_name} is occupied by $ft_id.\n$_\n";
			}
		}
		## feature type
		unless( exists($Feature_Basic_Href->{type}->{ $ft_Ob{type} }) ){
			warn_and_exit "<WARN>:\tfeature type $ft_Ob{type} is unknown.\n$_\n";
		}
		# load feature attributes
		$Feature_Basic_Href->{feature}->{ $ft_Ob{feature_id} } = {};
		my $fet_id_Href = $Feature_Basic_Href->{feature}->{ $ft_Ob{feature_id} };
		$fet_id_Href->{type}        = $ft_Ob{type};
		$fet_id_Href->{show_text}   = $ft_Ob{show_text};
		$fet_id_Href->{short_name}  = $ft_Ob{short_name};
		$fet_id_Href->{sn_index}    = $ft_Ob{sn_index};
		$fet_id_Href->{show_bool}   = 0;
		$fet_id_Href->{counts}      = { host=>0, virus=>0 };
	}
	close FTS;

	# inform
	stout_and_sterr "[INFO]\tload features source type OK.\n";
}

#--- load features source types information ---
sub load_features_list{

	# variants
	my $Sample_Href = $V_Href->{Sample_Href};
	my $Segment_Href = $V_Href->{Segment_Href};
	my $Feature_Basic_Href = $V_Href->{Feature_Basic_Href};
	my $Feature_list = $V_Href->{Sample_Break_Feature};

	# read feature list setting
	open (FL,$Feature_list) || die "fail open $Feature_list: $!\n";
	# theme
	# Sample \t chr \t Position\t Partner \t Feature_1 \t Feature_2
	(my $theme_line = <FL>) =~ s/^#//;
	my @theme_tag = split /\s+/, $theme_line; # use '\s+' to split.
	my $basic_info_count = 4;
	$theme_tag[$_] = lc($theme_tag[$_]) for (0 .. $basic_info_count-1); # convert fixed-info to lower cases.
	# check Feature
	for my $i ($basic_info_count .. $#theme_tag){
		warn_and_exit "Cannot find feature with id '$theme_tag[$i]'.\n"
			unless exists($Feature_Basic_Href->{feature}->{ $theme_tag[$i] });
	}
	# load in
	while (<FL>){
		next if(/^\#/);
		chomp;
		my @info = split /\t+/; # use '\t+' to split.
		my %ft_Ob = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
		# check the existence and duplications
		my $sample_id = $ft_Ob{sample};
		next unless(exists($Sample_Href->{$sample_id}));
		my $sam_OB = $Sample_Href->{$sample_id};
		my $refseg = $ft_Ob{chr};
		my $position = $ft_Ob{position};
		my $partner = $ft_Ob{partner};
		my $integ_case = "$refseg:$position";
		# check partner
		if( $partner ne 'host' && $partner ne 'virus' ){
			warn_and_exit "Only allow 'host' and 'virus' partners.\n$Feature_list\n";
		}
		# find the seg_OB
		my $seg_NO = FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->which_host_segment( refseg=>$refseg, lf_pos=>$position );
		next if( $seg_NO == -1 ); # cannot find segments
		# load the features to samples
		for my $i ($basic_info_count .. $#theme_tag){
			my $feature_id = $theme_tag[$i];
			if($ft_Ob{$feature_id} eq 'Y'){ # has this feature.
				$sam_OB->loadup_pos_features( partner=>$partner, integ_case=>$integ_case, feature_id=>$feature_id );
				# note this feature id show in legend later.
				if( $partner eq 'virus' || $V_Href->{show_hhf} ){
					$Feature_Basic_Href->{feature}->{$feature_id}->{show_bool} = 1;
				}
				# count this feature
				$Feature_Basic_Href->{feature}->{$feature_id}->{counts}->{$partner} ++;
			}
		}
	}
	close FL;

	# distribute positions with features to segments
	for my $sample_id (sort keys %$Sample_Href){
		my $sam_OB = $Sample_Href->{$sample_id};
		for my $partner (keys %{$sam_OB->{feature}}){
			next if( !$V_Href->{show_hhf} && $partner ne 'virus');
			for my $integcase (sort keys %{$sam_OB->{feature}->{$partner}}){
				my ($refseg, $pos) = (split /:/, $integcase)[0,1];
				# find the seg_OB
				my $seg_NO = FuseSV::Virus_Integ::LocalHaplotypeCircos::Segment_Deal->which_host_segment( refseg=>$refseg, lf_pos=>$pos );
				next if( $seg_NO == -1 ); # cannot find segments, filter before, again~
				my $seg_OB = $Segment_Href->{$seg_NO};
				# load this pos-features into this seg_OB
				$seg_OB->load_pos_features( sam_OB=>$sam_OB, refseg=>$refseg, pos=>$pos,
											feature_Aref=>$sam_OB->{feature}->{$partner}->{$integcase},
											partner=>$partner
										  );
			}
		}
	}

	# to avoid overlap of feature icons, adjust the icon's radian for each segment
	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		next unless(exists($seg_OB->{pos_feature}));
		$seg_OB->arrange_features_icon;
	}
}

#--- draw feature icons on each segment
sub draw_feature_icon_on_seg{
	# options
	shift if ($_[0] =~ /::$MODULE_NAME/);
	my %parm = @_;
	my $draw_bool = $parm{draw_bool};

	# variants
	my $Segment_Href = $V_Href->{Segment_Href};

	for my $seg_NO (sort {$a<=>$b} keys %$Segment_Href){
		my $seg_OB = $Segment_Href->{$seg_NO};
		next unless(exists($seg_OB->{pos_feature}));
		$seg_OB->draw_features_icon( draw_bool => $draw_bool );
	}
}

#--- 
1; ## tell the perl script the successful access of this module.
