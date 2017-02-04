#! /usr/bin/perl

use strict;
use warnings;

use Benchmark;
use Term::ANSIColor;
use File::Path qw(mkpath);

my $DATA_DIR = 'test-07';
my $MIN_DT = 0.000138;
my $REP_NUM = 7;

sub print_cmd {
    my $cmd = "$_[0]";
    print colored['bold yellow'], 'Run cmd: ';
    print $cmd."\n";
}
## Check the existence of dir. if not, create it.
sub make_dir {
    my $dir = $_[0];
    if ( -e "$dir" ) {
        die colored ['bold red'], $dir.'exists but not a directory'."\n" if ! -d "$dir";
    } else {
        mkpath $dir;
    }
}

make_dir("$DATA_DIR/exp");
make_dir("$DATA_DIR/imp");
make_dir("$DATA_DIR/mix");

# --------------------------------------------------------------------------------------

my $cmd;
###############################################################################################
my $dt=$MIN_DT/3.;

$cmd = "src/ExplicitTDHJ2D -X 3072 -Y 3072 -T 4 -d $dt -o ${DATA_DIR}/exp/u3072_${dt}_4.npy";
print_cmd $cmd;
system $cmd;
###############################################################################################

$dt=$MIN_DT;
open(EXPT_FD, '>exp_time.txt') or die "Cannot open file\n";
for(my $res=1024;$res >= 128;$res /= 2) {
    $cmd = "src/ExplicitTDHJ2D -X $res -Y $res -T 4 -d $dt -o ${DATA_DIR}/exp/u${res}_${dt}_4.npy";
    print_cmd $cmd;

    #### run the command 
    my $t0 = Benchmark->new;
    system $cmd;
    my $t1 = Benchmark->new;
    my $diffs = timediff($t1, $t0);
    printf STDERR "\nreal %.2f\nuser %.2f\nsys  %.2f\n", @$diffs[0,3,4];
    print EXPT_FD @$diffs[3]."\n";
    $dt *= 2;
}
close EXPT_FD;

$dt=$MIN_DT;
open(IMPT_FD, '>imp_time.txt') or die "Cannot open file\n";
for(my $res=1024;$res >= 256;$res /= 2) {
    my $cdt = $dt;
    for(my $i=0;$i<$REP_NUM;++$i) {
        $cmd = "src/ImplicitTDHJ2D -X $res -Y $res -T 4 -d $cdt -o ${DATA_DIR}/imp/u${res}_${cdt}_4.npy";
        print_cmd $cmd;

        #### run the command 
        my $t0 = Benchmark->new;
        system $cmd;
        my $t1 = Benchmark->new;
        my $diffs = timediff($t1, $t0);
        printf STDERR "\nreal %.2f\nuser %.2f\nsys  %.2f\n", @$diffs[0,3,4];
        print IMPT_FD @$diffs[3].' ';

        $cdt *= 2;
    }
    print IMPT_FD "\n";
    $dt *= 2;
}
close IMPT_FD;

$dt=$MIN_DT;
open(MIXT_FD, '>mix_time.txt') or die "Cannot open file\n";
for(my $res=1024;$res >= 256;$res /= 2) {
    my $cdt = $dt;
    for(my $i=0;$i<$REP_NUM;++$i) {
        $cmd = "src/MixedTDHJ2D -X $res -Y $res -T 4 -d $cdt -o ${DATA_DIR}/mix/u${res}_${cdt}_4.npy";
        print_cmd $cmd;

        #### run the command 
        my $t0 = Benchmark->new;
        system $cmd;
        my $t1 = Benchmark->new;
        my $diffs = timediff($t1, $t0);
        printf STDERR "\nreal %.2f\nuser %.2f\nsys  %.2f\n", @$diffs[0,3,4];
        print MIXT_FD @$diffs[3].' ';

        $cdt *= 2;
    }
    print MIXT_FD "\n";
    $dt *= 2;
}
close MIXT_FD;
