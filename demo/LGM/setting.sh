# TOOLS
# FuseSV
# https://github.com/deepomicslab/FuseSV
FuseSV=$PWD/../../FuseSV.pl;
# LocalHaplotypeSolver
# https://github.com/deepomicslab/LocalHaplotypeSolver
LHS_DIR=__THE_LOCALHAPLOTYPESOLVER_FOLDER_PATH__;
# SAMtools, version: 1.3.1
SAMtools=__THE_SAMTOOLS_PATH__;
# tabix, version: 0.2.5
tabix=__THE_TABIX_PATH__;

# DATABASE
# transcript psl bgzip file with tabix index
# see BioFuse, https://github.com/Nobel-Justin/BioFuse
# see SOAPfuse blog, https://sourceforge.net/p/soapfuse/blog
#     'New PSL file format applied by SOAPfuse'
# or, you can use the shared file (hg19, ensembl v75):
#     https://www.dropbox.com/sh/qsjsvfgeq5lcgbs/AADgFgv7B8_CyoNco-T4_Bh-a?dl=0
tpsl_bgz=__THE_TPSL_BGZ_PATH__;
