#verion=v1.1
#Wenlong JIA, 2021-01-11

sample=null;
st_step=1;
ed_step=3;
randtm=0;
# Parsing arguments
usage="Usage: `basename $0` -s <sample> -p <st_step[1]> -e <ed_step[3]> -t <randtm[0]>"
while getopts 'p:e:t:s:h' OPT; do
    case $OPT in
        h)
            echo $usage; exit;;
        p)
            st_step=$OPTARG;;
        e)
            ed_step=$OPTARG;;
        s)
            sample=$OPTARG;;
        t)
            randtm=$OPTARG;;
        ?)
            echo $usage; exit;;
    esac
done

if [ $sample = "null" ]; then
  echo $usage; exit 1;
fi

# load variables
source ./setting.sh;

pm_dir=`dirname $FuseSV`;
pm_dir=$pm_dir/..;
if [ ! -e $pm_dir ] || [ ! -e "$pm_dir/FuseSV" ]; then
  echo "cannot find perl module folder: $pm_dir";
  exit 1;
fi

sample_dir=$PWD/$sample;
if [ ! -e $sample_dir ]; then
  echo "cannot find sample folder: $sample_dir";
  exit 1;
fi
# add module
PERL5LIB=$pm_dir:$PERL5LIB; export PERL5LIB;

cd $sample_dir/find_lgm;

# path check
for file in FuseSV LHS_DIR SAMtools tabix tpsl_bgz; do
  if [ ! -e $(eval echo \$$file) ]; then
    echo "please use correct absolute path for '$file' in setting.sh";
    exit;
  fi;
done

#--------------------------------#
# Step 1, create config for LGM  #
# and draw SVG of local genome.  #
#--------------------------------#
if [ $st_step -le 1 ] && [ $ed_step -ge 1 ]; then
  echo -e "\n#-------- step 01 get_segCN --------#\n";
  # get step01 config
  if [ ! -e ./step01_config.txt ]; then
    echo "cannot find step01_config.txt.";
    exit 1;
  fi
  # load variables
  source ./step01_config.txt;
  # draw_segCN
  perl $FuseSV draw_segCN \
    -id $caseID \
    -reg $chr:$st_pos-$ed_pos \
    -t_bam $caseBam \
    -t_mdp $t_mean_dp \
    -t_puty $t_purity \
    -t_pldy $t_ploidy \
    -stl $SAMtools -tbi $tabix \
    -tpsl $tpsl_bgz \
    $reg_resol_para $juncFile_para $vsegFile_para \
    $x_axis_para $y_axis_para $segm_para $bg_para \
    -o $out_svg;
  sleep 1;
  # add ploidy
  step01_dir=`dirname $out_svg`;
  cat $s01_out_pref.InputForLocalHapAlgorithm.txt | sed "s/_please_fix_me_/$ploidy/" > $step01_dir/$patientID.InputForLocalHapAlgorithm.txt;
fi

#-------------------------#
# Step 2, get unit-cycles #
#-------------------------#
if [ $st_step -le 2 ] && [ $ed_step -ge 2 ]; then
  echo -e "\n#-------- step 02 get_ucyc --------#\n";
  # para variable
  para="-s 4 -SL 1 -m 100 --cap_sv_weight 0.5 --cap_jun_weight 0.8 --auto_weight --lp_with_original_ploidy";
  # prepare input
  input=step01_segCN/$sample.InputForLocalHapAlgorithm.txt;
  if [ ! -e $input ]; then
    echo "cannot find input from step01_segCN.";
    exit 1;
  fi
  # major and minor allele
  ploidy=`grep '^PLOIDY' $input | cut -f2`;
  allCN=${ploidy:0:1};
  minorCN=${ploidy:2:1};
  ((majorCN=$allCN-$minorCN));
  # get_UCYC
  mkdir -p step02_ucyc;
  # different allele
  for ref_allele in Major minor; do
    # cannot drop Major when LOH
    if [ $minorCN = 0 ] && [ $ref_allele = "Major" ]; then continue; fi
    # skip drop Major when M=m
    if [ $minorCN = $majorCN ] && [ $ref_allele = "Major" ]; then continue; fi
    # tag
    ref_allele_para=${ref_allele:0:1};
    tag="vItg_on_minor"; # intial, virus integrated at the mior allele
    if [ $ref_allele = "minor" ]; then
      if [ $minorCN = $majorCN ]; then
        tag="vItg"; # balance
      else
        tag="vItg_on_Major"; # virus integrated at the Major allele
      fi
    fi
    if [ $minorCN = 0 ] && [ $ref_allele = "minor" ]; then
      tag="LOH"; # LOH
    fi
    # run!
    output=step02_ucyc/$sample.$tag.ucyc.txt;
    if [ -e $output ]; then rm $output; fi # sweep first
    perl $FuseSV get_UCYC \
      -lhs $LHS_DIR \
      -o $output \
      $para --equally_assign_cycles_to_alleles \
      --drop_ref_allele $ref_allele_para \
      $input;
  done
  sleep 1;
fi

#----------------------------------#
# Step 3, merge unit-cycles to LGM #
#----------------------------------#
if [ $st_step -le 3 ] && [ $ed_step -ge 3 ]; then
  echo -e "\n#-------- step 03 get_LGM --------#\n";
  # settings
  ## randtm=0 just get simplest local haplotype
  ## randtm=N>0 could enable random process in N times
  ls ./step02_ucyc/$sample*.ucyc.txt | while read input; do
    pref=`basename $input | sed 's/.ucyc.txt//'`;
    # input
    ucyc_config=./step02_ucyc/$pref.ucyc.txt;
    # output
    output_LGM=./step03_LGM/$pref.LGM.randtm_$randtm.txt;
    # UCYCtoLM
    mkdir -p step03_LGM;
    perl $FuseSV UCYCtoLM \
      -randtm $randtm \
      -config $ucyc_config \
      -output $output_LGM;
  done
  sleep 1;
fi

