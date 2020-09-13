#!/usr/bin/perl -w
use strict;
use FuseSV::FuseSV;
use FuseSV::RunFunc;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "0.13";
$DATE = '2018-07-24';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

FuseSV::RunFunc->options_alert_and_run( argv_Aref => \@ARGV );
