# FuseSV

Tools for oncovirus analysis from NGS data.

- Author: Wenlong Jia
- Email:  wenlongkxm@gmail.com

## Version
0.28

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
- [SOAPfuse](https://github.com/Nobel-Justin/SOAPfuse)
- [BioFuse](https://github.com/Nobel-Justin/BioFuse)

If you encounter problems, please open an issue at the [project on Github](https://github.com/deepomicslab/FuseSV).

## Commands

Currently, FuseSV provides commands.

- draw_segCN

  visualize host local genome segmentations, get 'get_UCYC' config.

- get_UCYC
  
  get unit-cycles of local map.

  It packages '[LocalHaplotypeSolver](https://github.com/deepomicslab/LocalHaplotypeSolver)', which requires Python2.7 and PuLP package for linear programming.

  ```bash
  python -m pip install pulp
  ```

- UCYCtoLM

  merge unit-cycle to bio-contig in local map of virus integration.
