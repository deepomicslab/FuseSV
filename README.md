# FuseSV

Tools for oncovirus analysis from NGS data.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Version
0.29

## Installation

FuseSV is written in PERL. All of its functions are packaged into a standalone PERL module. Besides, it also requires other additional modules.

To install FuseSV, you need to download FuseSV and add the current directory to the `PERL5LIB` path.
```bash
git clone https://github.com/deepomicslab/FuseSV.git
PERL5LIB=$PERL5LIB:$PWD; export PERL5LIB
```
List of additional PERL modules required:
- [JSON](https://metacpan.org/pod/JSON)
- [List::Util](https://metacpan.org/pod/List::Util)
- [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager)
- [Math::Trig](https://metacpan.org/pod/Math::Trig)
- [POSIX](https://metacpan.org/pod/distribution/perl/ext/POSIX/lib/POSIX.pod)
- [SVG](https://metacpan.org/pod/SVG)
- [BioFuse](https://github.com/Nobel-Justin/BioFuse)

If you encounter problems, please open an issue at the [project on Github](https://github.com/deepomicslab/FuseSV).

## Commands

Currently, FuseSV provides commands.

- [draw_segCN](./manual/draw_segCN.md)

  visualize host local genome segmentations, get input of `get_UCYC`.

- [get_UCYC](./manual/get_UCYC.md)
  
  get unit-cycles of local genomic map, input of `UCYCtoLM`.

  It packages [`LocalHaplotypeSolver`](https://github.com/deepomicslab/LocalHaplotypeSolver), which requires Python2.7 and PuLP-1.6.8 package for linear programming.

  ```bash
  python -m pip install pulp==1.6.8
  ```

- [UCYCtoLM](./manual/UCYCtoLM.md)

  merge unit-cycles to local genomic map.

## Demo case

### Local Genomic Map (LGM)

one demo case is provided for the LGM analysis pipeline. Please check [demo](./demo/LGM).

PS. More cases could be found in the submitted paper and the published repository.
