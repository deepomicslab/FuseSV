# draw_segCN
`draw_segCN` visualizes host local genome segmentations, gets the input of [`get_UCYC`](./get_UCYC.md) .

## Usage

```
Usage:   perl FuseSV.pl draw_segCN <[Options]>

Options:

 # Inputs and Outputs #
  -t_bam  [s]  indexed bam file of tumor (case) sample. <[optional]>
  -t_mdp  [f]  average sequencing depth of tumor (case) sample. <required with '-t_bam'>
  -n_bam  [s]  indexed bam file of normal (control) sample. <[optional]>
  -n_mdp  [f]  average sequencing depth of normal (control) sample. <required with '-n_bam'>
  -reg    [s]  genomic region to draw, e.g., 'chr12:100000-120000'. <required>
  -rgsol  [s]  resolution config file of each region. <required>
  -o      [s]  output figure file, allowed format: svg/png/pdf. <required>

 # Database and Metadata #
  -tpsl   [s]  transcript PSL annotation file, should be indexed bgz file. [optional]
  -rmsk   [s]  repeat masker annotation file, should be indexed bgz file. [optional]
  -eltp   [s]  RMSK element types to show. [all available]
  -gene   [s]  gene to show. [all available]
  -nreg   [s]  specific region(s) in '-reg' to display. [none]
                instance: -nreg 'TEST1:110000-110200' -nreg 'TEST2:110400-110700,111400-111700'

 # Software Required #
  -stl    [s]  SamTools, minimum version: 1.3. <required>
  -tbi    [s]  Tabix path. <required>

 # Sample Attributes #
  -id     [s]  sample ID to show in the figure. [optional]
  -t_puty [f]  purity of tumor (case) sample. [1]
  -t_pldy [f]  ploidy of tumor (case) sample. [2]
  -n_pldy [f]  ploidy of normal (control) sample. [2]
  -junc   [s]  breakage junctions info file. [optional]
               denotes the virus integration position(s) and rearrangement break-point(s) if has.
  -vseg   [s]  viral segments info. [optional]
               to calculate CN of virus segments for local genomic map construction.

 # Options about Depth #
  -baseQ  [i]  minimum base quality for 'SamTools depth'. [5]
  -mapQ   [i]  minimum mapping quality for 'SamTools depth'. [10]
  -no_sfcp     ignore the Softclip reads for depth. [disabled]
  -tmrt   [f]  ratio to trim for even depth of one region. [0]
               opsitive value: bilateral operation; negative: only tail.
  -eng3sd      apply engineer_3times_SD_filter for average depth. [disabled]
               Once enabled, '-tmrt' will be ignored.
  -fwsm   [i]  how many bilateral flanking windows to make curve smooth. [2]
                note: value '2' means 2(fore)+1(this)+2(after); set '0' to disable.
  -n_dptr [f]  depth ratio selected from normal tissue for regional depth adjustment. [0]
               defaults is zero, is disable this operation.
  -nprcr  [s]  the abnormal copy ratio in normal (control) sample, multiple inputs ok.
               format: region_start_pos-region_end_pos:copyratio

 # Options about Display #
 ## colour
  -T_col  [s]  the colour of depth spectrum from tumor (case) sample. [red]
  -N_col  [s]  the colour of depth spectrum from normal (control) sample. [dodgerblue]
  -bgbcol [s]  the colour of background boundary. [black]
 ## axis
  -x_lnb  [i]  label number basement shown on x-axis. [1000]
  -x_lbs  [i]  how many screen pixels to show one position label on x-axis. [200]
  -x_lbfz [i]  font size of label on x-axis. [12]
  -x_rrc  [s]  color to show gradient region resolution. [none]
  -y_len  [i]  length of y-axis, will be overwritten when '-y_res' is set. [80]
  -y_res  [f]  how many depth (X) one spot in figure represents on y-axis. [auto]
  -y_lbs  [i]  how many screen pixels to show one depth marker on y-axis, effective with '-y_res'. [auto]
  -y_lbfz [i]  font size of label on y-axis. [12]
 ## segmentation * effective with available '-bkpos' input *
  -segm   [s]  show segmentation of this region. [u]
               note: 1) use 'u' for Upper-case (A-Z) naming; 'l' for Lower-case (a-z);
                     2) other inputs will act as prefix postfixed by digital.
                     3) first index can determined as inputs like 'u:C', 'l:d', or 'blabla:3'.
                     4) blank (regex detection: \s) is not allowed.
  -segc   [s]  segmentation bg-color. ['lightblue']
  -cn_bs  [s]  calculate copy number (CN) of each segments according to this region. [optional]
               format: 'T:chr12:105000-110000:2'
               note: 1) otherwise, use 'HapChrDepth mode' to calculate CN when tumor sample is provided.
                     2) set this option if you want to calculate CN when only normal sample is provided.
                     3) 'T' is tissue type, you could also use 'N' when '-n_bam' is given.
                     4) the last number is copy number statement of this region.
                     5) show CN of 'case' once '-t_bam' is set, otherwise show CN of 'control'.
                     6) CN information will auto output to the file name 'prefix(para(-o)).SegCN_info.tsv'.
 ## others
  -auto_r [i]  auto adjusted resolution to make genomic region shown in such size of pixels. [disabled]
  -miew   [f]  the minimum width of elements (exon/CDS/RMSKele) in the figure. [1]
  -extwd  [i]  pixels to extend the width  for enough space to display. [minimum:150]
  -extht  [i]  pixels to extend the height for enough space to display. [minimum:150]

  -h|help      Display this help info.

Version:
  5.40 at 2021-01-11

Author:
  Wenlong Jia (wenlongkxm@gmail.com)

```
