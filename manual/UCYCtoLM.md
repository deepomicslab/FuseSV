# UCYCtoLM
`UCYCtoLM` merges unit-cycles to local genomic map.

## Usage

```
Usage:   perl FuseSV.pl UCYCtoLM <[Options]>

Options:

  # Inputs and Outputs #
   -config [s]  unit-cycle config from local haplotype algorithm. <required>
   -output [s]  output local map file. <required>

  # Options #
   -randtm [i]  how many times to try construct local-map(s) in random mode. [0]
                 Note: 0 means simplest local-map mode.
   -segstr [s]  local-map must have such segments string(s) in random mode.
                 Note: 1) use as, -segstr 'H2+;V3-;H3+;H5+'
                       2) allows multiple times: -segstr 'xx' -segstr 'xx' ...
                       3) allows omit seg-orit (+/-), e.g., 'H2;V3-;H3;H5+'
   -solu_y [s]  select these SOLUTION(s) to analysis. [disabled]
                 Note: 1) use like, -solu_y 1,3,6
                       2) find 'SOLUTION NO.' in config file.
   -solu_n [s]  avoid these SOLUTION(s) from analysis. [disabled]
                 Note: 1) use like, -solu_n 2,7,10
                       2) it has higher priority than '-solu_y'.
   -usefvgm     use FVGM (Free Viral GenoMe, if has) in random mode. [disabled]
   -pwlink [s]  the file recodes pair-wise anchor link counts. [optional]
                 Note: please check FuseSV online instructions to know about.
   -lrdist [s]  file recodes distribution of long-range DNA length. [optional]
                 Note: 1) e.g., the barcode cover length in 10x-seq data.
                       2) effective with '-pwlink'.
                       3) check FuseSV online instructions to know about.
   -corrmd [s]  method to calculate pw-anchor links' correlations. [p]
                 Note: 'p' for pearson, 's' for spearman.
   -stype  [s]  sequencing type of long-range DNA. ['10x']
                 Note: valid types are '10x', 'Hi-C', MinIon', and 'PacBio'.
   -depth  [i]  simualted sequencing depth of the long-range DNA. [50]
   -perlen [i]  read length of paired-end (PE) sequencing in simulation. [150]
                 Note: 1) effective with '-stype' as 'Hi-C'.
                       2) use it for both ends of PE-reads.
   -fork   [i]  to run simualted sequencing with N forks in parallel. [1]

   -h|help      Display this help info.

Version:
   0.13 at 2018-05-31

Author:
   Wenlong Jia (wenlongkxm@gmail.com)
```
