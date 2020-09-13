This demo case is the 'HPV16_LRP1B_Complex' simulated case in our submitted paper.
More cases could be found in the paper and the published repository.

The host genome reference is GRCh37.

< Content >

1. bam folder
   a) The alignemnt of sequencing reads of the simulated HPV16-integrated local genomic map.
   b) alignment reference is GRCh37.

2. LGM folder
   a) step01_segCN folder stores inputs and outputs of the step01 (FuseSV draw_segCN).
     - input
       > junctions.info.txt
          This case is simulated to have two HPV16 integrations and two host SVs at gene LRP1B locus.
       > region.resol.txt
          The local host genomic region is divided into four segments by the junction breakpoints.
       > virusSeg.info.txt
          The HPV16 genome (NC_001526.4) is divided into three segments by the junction breakpoints.
     - output
       > HPV16_LRP1B_Complex.chr2_141740000_141806500.svg
       > HPV16_LRP1B_Complex.InputForLocalHapAlgorithm.txt
          You may find two 'InputForLocalHapAlgorithm.txt' files. The one with
          longer name is the initial output of step01_segCN. And after updating the
          '_please_fix_me_' field with results from Patchwork, we get the one with
          a shorter name, which is the input of the second step.
   b) step02_ucyc folder stores outputs of the step02 (FuseSV get_UCYC).
      > HPV16_LRP1B_Complex.LOH.ucyc.txt
         It contains the basic contigs that form the LGM. We call them unit-cycle (UCYC).
         This file is the input of the third step.
   c) step03_LGM folder stores outputs of the step03 (FuseSV UCYCtoLM).
      > HPV16_LRP1B_Complex.LOH.LGM.randtm_0.txt
         This is the final result showing LGM of the virus integrations. The
         result corresponds to the 'Simplest LGM' report in our paper.
     b) If non-zero randtm is set, e.g., 100, the file will be named with 'randtm_100'.

< Operation >

1. set variables in the 'setting.sh' file
2. run the 'run.sh' file in this folder
    bash run.sh -s HPV16_LRP1B_Complex


Wenlong Jia
wenlongkxm@gmail.com
City University of Hong Kong
2020-09-12
