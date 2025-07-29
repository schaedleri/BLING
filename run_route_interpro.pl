#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename;

# === オプション処理 ===
my %opt = (
    base_dir => '.',
    cpu      => 8,
    appl     => 'TIGRFAM,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,SMART,CDD,PRINTS,PIRSR,ProSitePatterns,Pfam',
);
GetOptions(
    'base_dir=s' => \$opt{base_dir},
    'cpu=i'      => \$opt{cpu},
    'appl=s'     => \$opt{appl},
) or die "Usage: $0 --base_dir DIR [--cpu N] [--appl APP1,APP2,...]\n";

my $done_dir = "$opt{base_dir}/done";
make_path($done_dir) unless -d $done_dir;

# === ステップ定義 ===
my @steps = (
    { name => 'interpro',     script => 'scripts/microbiome_run_interpro.pl', depends_on => ['extract'], pass_opts => { base_dir => 1, cpu => 1, appl => 1 } },
    { name => 'interval',     script => 'scripts/microbiome_interval.pl',     depends_on => ['interpro'], pass_opts => { base_dir => 1 } },
    { name => 'merge_GO',     script => 'scripts/microbiome_merge_GO.pl',     depends_on => ['interval'], pass_opts => { base_dir => 1 } },
);

# === 実行管理 ===
my %done;
foreach my $step (@steps) {
    my $flag = "$done_dir/$step->{name}.done";

    if (-e $flag) {
        print "? Skip $step->{name} (done)\n";
        $done{$step->{name}} = 1;
        next;
    }

    for my $dep (@{ $step->{depends_on} }) {
        die "? Cannot run $step->{name} before $dep is completed\n" unless $done{$dep};
    }

    print "\n=== STEP: $step->{name} ===\n";

    my @cmd = ('perl', "$opt{base_dir}/$step->{script}");
    if ($step->{pass_opts}) {
        push @cmd, ('--base-dir', $opt{base_dir}) if $step->{pass_opts}{base_dir};
        push @cmd, ('--cpu', $opt{cpu})           if $step->{pass_opts}{cpu};
        push @cmd, ('--appl', $opt{appl})         if $step->{pass_opts}{appl};
    }

    print "? Running: @cmd\n";
    my $status = system(@cmd);
    if ($status != 0) {
        die "? Step $step->{name} failed (exit: $status)\n";
    }

    open my $fh, '>', $flag or die "? Cannot write $flag\n";
    print $fh "done\n";
    close $fh;

    $done{$step->{name}} = 1;
    print "? Completed $step->{name}\n";
}

print "\n? All steps completed successfully.\n";
