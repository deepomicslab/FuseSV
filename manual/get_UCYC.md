# get_UCYC
`get_UCYC` gets unit-cycles of local genomic map, the input of `[UCYCtoLM](./UCYCtoLM.md)`.

## Usage

```
Usage:   perl FuseSV.pl get_UCYC <[Options]>

Options:

  -lhs [s]  path of the 'LocalHaplotypeSolver' folder. <required>
            Download from https://github.com/deepomicslab/LocalHaplotypeSolver

Version:
  0.02 at 2020-09-13

Author:
  Wenlong Jia (wenlongkxm@gmail.com)
```

### LocalHaplotypeSolver

`get_UCYC` packages `[LocalHaplotypeSolver](https://github.com/deepomicslab/LocalHaplotypeSolver)`.
Please check the detailed helps info in its repository.

```
usage: LocalHap.py [-h] [-o OUTPUT_FILE] [--host_seg_weight_factor]
                   [--virus_seg_weight_factor] [--junction_weight_factor]
                   [--cap_sv_weight] [--cap_integ_weight]
                   [--inferred_jun_weight] [--cap_jun_weight] [--auto_weight]
                   [--weight_by_length] [--weight_by_cn] [--no_infer_normal]
                   [--keep_ref_allele] [--lp_with_original_ploidy]
                   [--equally_assign_cycles_to_alleles] [--drop_ref_allele]
                   [-e] [--linear_virus] [-s] [-v] [-m] [-S] [-L] [--write_lp]
                   input_file
```
