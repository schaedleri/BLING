#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Glob ':glob';

# === オプション ===
my $base_dir = '.';
my $obo_file = '';
GetOptions(
    'base-dir=s' => \$base_dir,
    'go-obo=s'   => \$obo_file,
) or die "Usage: $0 --base-dir <path> [--go-obo <path>]\n";

$obo_file ||= "$base_dir/data/go.obo";
die "Cannot find GO OBO file at $obo_file\n" unless -e $obo_file;

# ========== 入力ファイルパターン ==========
my @patterns = (
    "$base_dir/microbiome/*/go_annotations/*.tsv",
    "$base_dir/microbiome/*/interval_test/*/IP_GO_pairs.tsv",
    "$base_dir/microbiome/*/interval_test/*/nonGO_GO_pairs*.tsv"
);

# ========== STEP 1: go.oboの親子関係を構築 ==========
my %is_a;
{
    open(my $fh, "<", $obo_file) or die "Cannot open $obo_file: $!";
    my ($current_go, @parents);
    while (<$fh>) {
        chomp;
        if (/^\[Term\]/) {
            $current_go = undef;
            @parents    = ();
        } elsif (/^id:\s+(GO:\d+)/) {
            $current_go = $1;
        } elsif (/^is_a:\s+(GO:\d+)/) {
            push @parents, $1;
        } elsif (/^(part_of|regulates|enabled_by|contributes_to):\s+(GO:\d+)/) {
            push @parents, $2;
        } elsif (/^\s*\$/ && $current_go) {
            $is_a{$current_go} = [@parents];
        }
    }
    close $fh;
}

# ========== STEP 2: GO-Geneペアを収集 ==========
my %go_pair;
my %ipr_pair;

sub get_all_ancestors {
    my ($go_id, $visited) = @_;
    $visited ||= {};
    return () if $visited->{$go_id}++;
    my @parents = @{ $is_a{$go_id} || [] };
    return (@parents, map { get_all_ancestors($_, $visited) } @parents);
}

foreach my $pattern (@patterns) {
    foreach my $file (glob($pattern)) {
        open my $fh, '<', $file or do { warn "Cannot open $file: $!"; next; };
        <$fh>;
        while (<$fh>) {
            chomp;
            my @cols = split /\t/;
            my ($go, $gene, $ipr_field);

            if ($file =~ m{/go_annotations/}) {
                next unless @cols >= 4;
                $gene = $cols[0];
                $go   = $cols[2] if $cols[2] =~ /^GO:/;
                $ipr_field = $cols[3];
            } elsif ($file =~ m{IP_GO_pairs\.tsv\$}) {
                next unless @cols >= 2;
                $gene = $cols[0];
                $go   = $cols[1] if $cols[1] =~ /^GO:/;
            } elsif ($file =~ m{nonGO_GO_pairs}) {
                next unless @cols >= 2;
                $gene = $cols[0];
                $go   = $cols[1] if $cols[1] =~ /^GO:/;
            } else {
                next;
            }

            if ($go) {
                my %all_gos = map { $_ => 1 } ($go, get_all_ancestors($go));
                $go_pair{"$_\t$gene"} = 1 for keys %all_gos;
            }

            if ($ipr_field && $ipr_field ne '-' && $ipr_field ne '') {
                my @iprs = grep { /^IPR\d+/ } split(/[;, ]+/, $ipr_field);
                for my $ipr (@iprs) {
                    $ipr_pair{"$ipr\t$gene"} = 1;
                }

            }
        }
        close $fh;
    }
}

my @interpro_files = glob("$base_dir/microbiome/*/interpro/*.tsv");
foreach my $interpro_file (@interpro_files) {
    open my $ipfh, '<', $interpro_file or do { warn "Cannot open $interpro_file: $!"; next; };
    while (<$ipfh>) {
        chomp;
        my @cols = split /\t/;
        next unless @cols > 11;
        my $gene = $cols[0];
        my $ipr_field = $cols[11];
        next unless $ipr_field && $ipr_field ne '-' && $ipr_field ne '';
        my @iprs = grep { /^IPR\d+/ } split(/[;, ]+/, $ipr_field);
        for my $ipr (@iprs) {
            $ipr_pair{"$ipr\t$gene"} = 1;
        }

    }
    close $ipfh;
}

open my $go_out, '>', 'go_gene.tsv' or die "Cannot write go_gene.tsv: $!";
print $go_out "$_\n" for sort keys %go_pair;
close $go_out;

open my $ipr_out, '>', 'ipr_gene.tsv' or die "Cannot write ipr_gene.tsv: $!";
print $ipr_out "$_\n" for sort keys %ipr_pair;
close $ipr_out;

my %gene;
for my $line (keys %go_pair, keys %ipr_pair) {
    my ($id, $gene) = split /\t/, $line;
    $gene{$gene} = 1 if defined $gene;
}
open my $bg, '>', 'background_list.txt' or die "Cannot write background_list.txt: $!";
print $bg "$_\n" for sort keys %gene;
close $bg;
