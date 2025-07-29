#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use File::Path qw(make_path);
use File::Basename;
use Parallel::ForkManager;
use Cwd qw(abs_path);

# === ベースディレクトリと出力先 ===
my $base_dir = abs_path($ARGV[0] // '.');
my $data_dir = "$base_dir/data";
my $microbiome_dir = "$base_dir/microbiome";

# rawgbffだけは data/ に配置
my $rawgbff_dir = "$data_dir/rawgbff";
make_path($rawgbff_dir) unless -d $rawgbff_dir;

# === 失敗ログファイルは microbiome 配下に
my $failed_file = "$microbiome_dir/failed_microbes_bacteria_taxid.txt";
make_path($microbiome_dir) unless -d $microbiome_dir;
open my $failed_fh, '>>', $failed_file or die "Cannot open failed file: $!";

my $pm = Parallel::ForkManager->new(1);  # 並列数調整可能

# === 入力TSVの探索 ===
my @tsv_files;
find(
    sub {
        return unless -f $_;
        return unless $_ =~ /^(Genus|Family|Order|Species)_.+\.tsv$/;
        push @tsv_files, $File::Find::name;
    },
    $microbiome_dir
);

my %seen;
foreach my $tsv_file (@tsv_files) {
    open my $fh, '<', $tsv_file or die "Cannot open $tsv_file: $!";
    my $header = <$fh>;  # skip header

    while (<$fh>) {
        chomp;
        next if /^\s*$/;
        my @cols = split /\t/;
        my $acc = $cols[0];
        next unless $acc;
        next if $seen{$acc}++;
        my $org = $cols[1] // 'unknown';
        $org =~ s/[^\w]/_/g;

        my $zipfile = "$rawgbff_dir/${org}_protein.zip";
        my $unzip_dir = "$rawgbff_dir/${org}_protein";

        # すでに解凍済みならスキップ
        if (-d $unzip_dir) {
            print "Skip existing: $unzip_dir\n";
            next;
        }

        $pm->start and next;

        # ダウンロード
        my $cmd = qq(datasets download genome accession $acc --assembly-source 'RefSeq' --reference --filename "$zipfile" --include gbff);
        my $result = system($cmd);

        my $attempt = 0;
        my $max_retries = 3;
        while ($result != 0 && $attempt < $max_retries) {
            $attempt++;
            warn "Retry $attempt: $acc ($org)\n";
            $result = system($cmd);
        }

        if ($result != 0) {
            my $exit_code = $? >> 8;
            print $failed_fh "$org (Accession: $acc) - Failed with exit code $exit_code\n";
            $pm->finish;
            next;
        }

        # 解凍
        make_path($unzip_dir);
        my $unzip_cmd = qq(unzip -q "$zipfile" -d "$unzip_dir");
        my $unzip_result = system($unzip_cmd);
        if ($unzip_result != 0) {
            warn "Failed to unzip: $zipfile\n";
        } else {
            print "Unzipped: $zipfile -> $unzip_dir\n";
        }

        $pm->finish;
    }
    close $fh;
}
$pm->wait_all_children;
close $failed_fh;

print "Download & unzip complete. Files in $rawgbff_dir\n";
print "Failed list: $failed_file\n";
