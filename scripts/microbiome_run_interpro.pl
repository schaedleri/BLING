#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use File::Glob ':glob';
use Cwd qw(abs_path);

# === ÉIÉvÉVÉáÉìéÊìæ ===
my $base_dir = '.';
my $interpro_dir_opt;
my $cpu = 8;
my $appl = 'TIGRFAM,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,SMART,CDD,PRINTS,PIRSR,ProSitePatterns,Pfam';

GetOptions(
    'base_dir=s'          => \$base_dir,
    'interproscan-dir=s'  => \$interpro_dir_opt,
    'cpu=i'               => \$cpu,
    'appl=s'              => \$appl,
) or die "Usage: $0 [--base_dir DIR] [--interproscan-dir DIR] [--cpu N] [--appl APPL1,APPL2,...]\n";

$base_dir = abs_path($base_dir);

my $interpro_base_dir = "$base_dir/tools";
my $interpro_dir;

if (defined $interpro_dir_opt) {
    $interpro_dir = $interpro_dir_opt;
} else {
    my @candidates = bsd_glob("$interpro_base_dir/interproscan-*");
    die "Error: Could not find interproscan directory under $interpro_base_dir\n" unless @candidates;
    $interpro_dir = $candidates[0];
}

my $interproscan = "$interpro_dir/interproscan.sh";
die "Error: $interproscan does not exist or is not executable\n" unless -x $interproscan;

my @levels = qw(Order Family Genus Species);
my %done_basename;

for my $level (@levels) {
    my @faa_files = glob("$base_dir/microbiome/${level}_*/faa/*.faa");
    for my $input_path (@faa_files) {
        my ($class_dir, $file) = $input_path =~ m{^($base_dir/microbiome/(?:Genus|Family|Order|Species)_[^/]+)/faa/([^/]+\.faa)$};
        next unless $class_dir && $file;

        my ($basename) = fileparse($file, qr/\.[^.]*/);
        my $interpro_out_dir = "$class_dir/interpro";
        make_path($interpro_out_dir) unless -d $interpro_out_dir;
        my $output_path = "$interpro_out_dir/$basename.tsv";

        if (-e $output_path) {
            print "[?] Found existing InterProScan output for $basename ? skipping.\n";
            $done_basename{$basename} = $output_path;
            next;
        }

        if (!$done_basename{$basename}) {
            print "[?] Running InterProScan for $input_path\n";
            my $cmd = "$interproscan -i $input_path -f tsv -o $output_path --goterms -cpu $cpu -appl $appl";
            system($cmd) == 0 or warn "[?] Failed to run InterProScan for $file\n";
            $done_basename{$basename} = $output_path;
        } else {
            my $from = $done_basename{$basename};
            print "[Å®] Copying result from $from to $output_path\n";
            unless (my_copy($from, $output_path)) {
                warn "[?] Failed to copy $from to $output_path\n";
            }
        }
    }
}

sub my_copy {
    my ($from, $to) = @_;
    open my $in,  '<', $from or return 0;
    open my $out, '>', $to   or return 0;
    binmode $in;
    binmode $out;
    while (my $buf = <$in>) {
        print $out $buf or return 0;
    }
    close $in;
    close $out;
    return 1;
}