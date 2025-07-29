#!/usr/bin/perl
use strict;
use warnings;
use File::Find;
use File::Basename;
use Getopt::Long;
use File::Spec;

my $base_dir;
GetOptions('base-dir=s' => \$base_dir) or die "Usage: $0 --base-dir <path> [統合ディレクトリ名]\n";
die "Error: --base-dir is required\n" unless defined $base_dir;

my @target_dirs = @ARGV;
die "Error: No target directories specified\n" unless @target_dirs == 1;

my $merged_dir = "$base_dir/microbiome/$target_dirs[0]";
my $faa_root   = "$merged_dir/faa";

# --- サブディレクトリごとにcd-hit ---
my @all_fasta_paths;
opendir(my $dh, $faa_root) or die "Cannot open $faa_root: $!";
my @subdirs = grep { -d "$faa_root/$_" && !/^\./ } readdir($dh);
closedir($dh);

for my $sub (@subdirs) {
    my $fasta_dir = "$faa_root/$sub";
    my @fasta_files;
    find(
        sub {
            push @fasta_files, $File::Find::name if /\.faa$/i && -f $_;
        },
        $fasta_dir
    );
    unless (@fasta_files) {
        warn "?? No faa files found under $fasta_dir\n";
        next;
    }
    my $output_fasta = "$fasta_dir/${sub}_all_sequences.fasta";
    my $cdhit_output = "$fasta_dir/${sub}_all_sequences_cdhit";
    open my $out, '>', $output_fasta or die "? Cannot open $output_fasta: $!";
    for my $file (@fasta_files) {
        open my $in, '<', $file or die "? Cannot open $file: $!";
        print $out $_ while (<$in>);
        close $in;
    }
    close $out;
    my $cmd = "cd-hit -i $output_fasta -p 1 -s 0.95 -aS 0.95 -o $cdhit_output -c 0.9";
    system($cmd) == 0
        or warn "? Failed to run cd-hit for $sub\n";
    print "? $sub のcd-hit出力: $cdhit_output\n";
    push @all_fasta_paths, $output_fasta;
}

# --- すべてのall_sequences.fastaをまとめて再cd-hit ---
my $merged_fasta = "$merged_dir/all_groups_merged.fasta";
my $merged_cdhit = "$merged_dir/all_groups_merged_cdhit";
if (@all_fasta_paths > 1) {
    open my $mout, '>', $merged_fasta or die "? Cannot open $merged_fasta: $!";
    for my $fasta (@all_fasta_paths) {
        open my $min, '<', $fasta or die "? Cannot open $fasta: $!";
        print $mout $_ while (<$min>);
        close $min;
    }
    close $mout;

    my $merge_cmd = "cd-hit -i $merged_fasta -p 1 -s 0.95 -aS 0.95 -o $merged_cdhit -c 0.9";
    system($merge_cmd) == 0
        or warn "? Failed to run merged cd-hit\n";

    print "? すべてのグループをまとめたcd-hit出力: $merged_cdhit\n";
}

# --- 各サブディレクトリのclstrファイルから代表配列IDを集計 ---
my %id_in_dirs; # 配列ID => { サブディレクトリ名 => 1, ... }
for my $sub (@subdirs) {
    my $clstr = "$faa_root/$sub/${sub}_all_sequences_cdhit.clstr";
    next unless -e $clstr;
    open my $cfh, '<', $clstr or die "Cannot open $clstr: $!";
    while (<$cfh>) {
        if (/^\d+\s+\d+aa, >(\S+)\.\.\.\s*\*$/) {
            $id_in_dirs{$1}{$sub} = 1;
        }
    }
    close $cfh;
}

# --- merged.fastaの配列をIDごとに記録 ---
my %merged_seq;
open my $fh, '<', $merged_fasta or die "Cannot open $merged_fasta: $!";
my ($id, $seq) = ('', '');
while (<$fh>) {
    if (/^>(\S+)/) {
        if ($id) { $merged_seq{$id} = $seq; }
        $id = $1;
        $seq = $_;
    } else {
        $seq .= $_;
    }
}
if ($id) { $merged_seq{$id} = $seq; }
close $fh;

# --- merged_cdhit.clstrをパースし、クラスターごとに代表配列と構成配列を集める ---
my $clstr_file = "$merged_dir/all_groups_merged_cdhit.clstr";
open my $cfh, '<', $clstr_file or die "Cannot open $clstr_file: $!";
my %cluster_to_ids;      # クラスター番号 => [配列ID, ...]
my %cluster_to_rep;      # クラスター番号 => 代表配列ID
my $cluster_num = -1;
while (<$cfh>) {
    if (/^>Cluster\s+(\d+)/) {
        $cluster_num = $1;
        next;
    }
    if (/^\d+\s+\d+aa, >(\S+)\.\.\.\s*(\*?)/) {
        my $seqid = $1;
        push @{ $cluster_to_ids{$cluster_num} }, $seqid;
        $cluster_to_rep{$cluster_num} = $seqid if $2; # *がついていれば代表
    }
}
close $cfh;

# --- クラスターごとに構成配列の由来サブディレクトリを調べて分類 ---
my %comb_to_repseqs; # 組み合わせ => [代表配列ID, ...]
for my $clnum (keys %cluster_to_ids) {
    my %dirs;
    for my $id (@{ $cluster_to_ids{$clnum} }) {
        for my $dir (keys %{ $id_in_dirs{$id} // {} }) {
            $dirs{$dir} = 1;
        }
    }
    my @comb = sort keys %dirs;
    my $comb = @comb ? join("_", @comb) : "none";
    my $rep_id = $cluster_to_rep{$clnum};
    push @{ $comb_to_repseqs{$comb} }, $rep_id if $rep_id;
}

# --- 組み合わせごとに代表配列をFASTA出力 ---
for my $comb (keys %comb_to_repseqs) {
    my $file = "$merged_dir/merged_cluster_representative_by_dir_${comb}.fasta";
    open my $out, '>', $file or die "Cannot open $file: $!";
    my %done;
    for my $id (@{ $comb_to_repseqs{$comb} }) {
        next if $done{$id}++;
        print $out $merged_seq{$id} if exists $merged_seq{$id};
    }
    close $out;
    print "? $file を出力しました\n";
}


