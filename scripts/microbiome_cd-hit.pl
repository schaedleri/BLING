#!/usr/bin/perl
use strict;
use warnings;
use File::Find;
use File::Basename;
use Getopt::Long;

my $base_dir;
GetOptions('base-dir=s' => \$base_dir) or die "Usage: $0 --base-dir <path> [Genus_Xxx Species_Yyy ...]\n";
die "Error: --base-dir is required\n" unless defined $base_dir;

# === 引数チェック ===
my @target_dirs = @ARGV;
die "Error: No target directories specified (e.g., Genus_Xxx)\n" unless @target_dirs;

# === 実行 ===
my $microbiome_dir = "$base_dir/microbiome";

for my $group (@target_dirs) {
    my $group_dir  = "$microbiome_dir/$group";
    my $fasta_dir  = "$group_dir/fasta";

    unless (-d $fasta_dir) {
        warn "?? No fasta directory: $fasta_dir\n";
        next;
    }

    # FASTA ファイル探索
    my @fasta_files;
    find(
        sub {
            push @fasta_files, $File::Find::name if /\.fasta$/i && -f $_;
        },
        $fasta_dir
    );

    unless (@fasta_files) {
        warn "?? No fasta files found under $fasta_dir\n";
        next;
    }

    # 出力ファイル名
    my $output_fasta = "$fasta_dir/${group}_all_sequences.fasta";
    my $cdhit_output = "$fasta_dir/${group}_all_sequences_cdhit";

    # 結合
    open my $out, '>', $output_fasta or die "? Cannot open $output_fasta: $!";
    for my $file (@fasta_files) {
        open my $in, '<', $file or die "? Cannot open $file: $!";
        print $out $_ while (<$in>);
        close $in;
    }
    close $out;

    # cd-hit 実行
    my $cmd = "cd-hit -i $output_fasta -p 1 -s 0.95 -aS 0.95 -o $cdhit_output -c 0.9";
    system($cmd) == 0
        or warn "? Failed to run cd-hit for $group\n";

    print "? $group のcd-hit出力: $cdhit_output\n";
}
