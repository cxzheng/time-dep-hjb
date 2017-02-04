#! /usr/bin/perl

use strict;
use warnings;

use Term::ANSIColor;

my $DATA_DIR = 'test-07';
my $MIN_DT = 0.000138;
my $REP_NUM = 7;
my $LERROR = '../script/Lerror.py';

sub print_cmd {
    my $cmd = "$_[0]";
    print colored['bold yellow'], 'Run cmd: ';
    print $cmd."\n";
}

my $cmd;
my $ret;
my $MIN_2048_DT = $MIN_DT / 2.;
my $dt=$MIN_DT;
open(ERR_FD, ">${DATA_DIR}/exp_err.txt") or die "Cannot open file\n";
for(my $res = 1024;$res >= 128; $res /= 2) {
    chomp($ret = `$LERROR ${DATA_DIR}/exp/u2048_${MIN_2048_DT}_4.npy ${DATA_DIR}/exp/u${res}_${dt}_4.npy`);
    print $dt.' '.$ret."\n";
    print ERR_FD $dt.' '.$ret."\n";
    $dt *= 2;
}
close ERR_FD;
print "===================\n";

$dt=$MIN_DT;
for(my $res=1024;$res >= 256;$res /= 2) {
    my $cdt = $dt;
    open(ERR_FD, ">${DATA_DIR}/imp_mix_err_${res}.txt") or die "Cannot open file\n";
    for(my $i=0;$i<$REP_NUM;++$i) {
        chomp($ret = `$LERROR ${DATA_DIR}/exp/u2048_${MIN_2048_DT}_4.npy ${DATA_DIR}/imp/u${res}_${cdt}_4.npy`);
        print ERR_FD $cdt.' '.$ret;
        chomp($ret = `$LERROR ${DATA_DIR}/exp/u2048_${MIN_2048_DT}_4.npy ${DATA_DIR}/mix/u${res}_${cdt}_4.npy`);
        print ERR_FD ' '.$ret."\n";
        #print ERR_FD "\n";
        print $ret."\n";
        $cdt *= 2;
    }
    close ERR_FD;
    $dt *= 2;
}
