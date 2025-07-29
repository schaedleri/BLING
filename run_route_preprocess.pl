#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;
use Parallel::ForkManager;

# === ベースディレクトリ指定 ===
my $base_dir = $ARGV[0] // '.';
my $done_dir = "$base_dir/done";
make_path($done_dir) unless -d $done_dir;

# === ステップ定義 ===
my @steps = (
    { name => 'assembly',       script => 'scripts/microbiome_assembly.pl',       depends_on => [] },
    { name => 'get_datasets',   script => 'scripts/microbiome_get_datasets.pl',   depends_on => ['assembly'] },
    { name => 'extract',        script => 'scripts/microbiome_extract.pl',        depends_on => ['get_datasets'] },

);

# === 実行 ===
my %done;
foreach my $step (@steps) {
    my $flag = "$done_dir/$step->{name}.done";
    if (-e $flag) {
        print "? Skip $step->{name} (done)\n";
        $done{$step->{name}} = 1;
        next;
    }

    # 依存チェック
    for my $dep (@{$step->{depends_on}}) {
        unless ($done{$dep}) {
            die "? Cannot run $step->{name} before $dep is completed\n";
        }
    }

    # 実行
    print "? Running $step->{name}\n";
    my $cmd = "perl $base_dir/$step->{script}";
    my $status = system($cmd);

    if ($status != 0) {
        die "? Step $step->{name} failed (exit: $status)\n";
    }

    open my $fh, '>', $flag or die "? Cannot write $flag\n";
    print $fh "done\n";
    close $fh;
    $done{$step->{name}} = 1;
    print "? Completed $step->{name}\n";
}

print "\n? All steps completed successfully\n";
