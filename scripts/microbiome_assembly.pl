#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use File::Basename;

# === 入力パス設定 ===
my $base_dir = abs_path($ARGV[0] // '.');
my $data_dir = "$base_dir/data";
my $input_file = "$data_dir/level-7.tsv";

# 出力ディレクトリ
my $output_root = "$base_dir/microbiome";
make_path($output_root) unless -d $output_root;
chdir $output_root or die "chdir failed: $!";

# コマンド出力ファイル
open my $in,     '<', $input_file or die "Cannot open $input_file: $!";
open my $cmdfh,  '>', 'commands.txt' or die $!;

# ヘッダー処理
my $header_line = <$in>;
chomp $header_line;
my @headers = split /\t/, $header_line;

while (<$in>) {
    chomp;
    my @cols = split /\t/;
    my $taxid = $cols[11];
    $taxid =~ s/^\s+|\s+$//g;
    next unless $taxid && $taxid =~ /^\d+$/;

    # 最も下位の分類名を取得
    my ($name, $header) = ('', '');
    for my $i (reverse 1..8) {
        if (defined $cols[$i] && $cols[$i] =~ /\S/) {
            $name   = $cols[$i];
            $header = $headers[$i];
            last;
        }
    }
    next unless $name;

    # 使えない文字を置換
    $name   =~ s/[\\\/:*?"<>|]/_/g;
    $name   =~ s/\s+/_/g;
    $header =~ s/[\\\/:*?"<>|]/_/g;
    $header =~ s/\s+/_/g;

    my $dirname = "${header}_$name";
    make_path($dirname) unless -d $dirname;
    my $output = "$dirname/${header}_$name.tsv";

    my $cmd = "datasets summary genome taxon $taxid --as-json-lines --reference | dataformat tsv genome --fields accession,organism-name,assminfo-level,organism-tax-id > $output";
    print $cmdfh "$cmd\n";
}
close $cmdfh;
close $in;

# コマンドの実行
open my $cmdin, '<', 'commands.txt' or die "Can't open commands.txt: $!";
while (my $cmd = <$cmdin>) {
    chomp $cmd;
    next unless $cmd =~ /\S/;
    print "Do: $cmd\n";
    system($cmd) == 0 or warn "failed: $cmd\n";
}
close $cmdin;
