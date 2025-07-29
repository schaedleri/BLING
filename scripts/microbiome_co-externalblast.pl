#!/usr/bin/env perl
use strict;
use warnings;
use File::Glob ':bsd_glob';
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use Fcntl qw(:flock);
use Parallel::ForkManager;
use Getopt::Long;


# --- ベースディレクトリの指定 ---
my $base_dir;
GetOptions(
    'base-dir=s' => \$base_dir,
) or die "Usage: $0 --base-dir BASEDIR SELECTED_GENUS...\n";
die "Usage: $0 --base-dir BASEDIR SELECTED_GENUS...\n" unless $base_dir && @ARGV;

foreach my $selected_genus (@ARGV) {
    my $input_dir = File::Spec->catdir($base_dir, "microbiome", $selected_genus);
    my $output_dir;
    if ($selected_genus =~ /\+/) {
        $output_dir = File::Spec->catdir($input_dir, "BLAST_specify");
        unless (-d $output_dir) {
            make_path($output_dir) or die "Cannot create $output_dir: $!";
        }
    } else {
        $output_dir = undef;
    }
    my $cpu = 10;
    my $taxid_str;
    my $pm = Parallel::ForkManager->new($cpu);
    my $project_root = $base_dir;

    my @dbs = (
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB1/BSDB1", number => "DB1"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB2/BSDB2", number => "DB2"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB3/BSDB3", number => "DB3"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB4/BSDB4", number => "DB4"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB5/BSDB5", number => "DB5"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB6/BSDB6", number => "DB6"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB7/BSDB7", number => "DB7"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB8/BSDB8", number => "DB8"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB9/BSDB9", number => "DB9"  },
        { db => "$project_root/DB/bacteria_strain_taxid_DB/DB0/BSDB0", number => "DB10" },
    );

    # --- Co-occurence.txtからマルチFASTAを切り出し ---
    my $cooccur_fasta = File::Spec->catfile($input_dir, "Co-occurence.txt");
    my @seqs;
    if (-e $cooccur_fasta) {
        open my $fh, '<', $cooccur_fasta or die "Cannot open $cooccur_fasta: $!";
        my ($header, $seq);
        while (my $line = <$fh>) {
            chomp $line;
            if ($line =~ /^>/) {
                if (defined $header && defined $seq) {
                    push @seqs, [$header, $seq];
                }
                $header = $line;
                $seq = '';
            } else {
                $seq .= $line;
            }
        }
        push @seqs, [$header, $seq] if defined $header && defined $seq;
        close $fh;
    }

    # --- taxid情報の読み込みとtaxid_str生成 ---
    my @all_dirs = grep { -d $_ && File::Basename::basename($_) !~ /\+/ }
        glob(File::Spec->catdir($base_dir, "microbiome", "*"));
    my @tsv_files = map { glob(File::Spec->catfile($_, "*.tsv")) } @all_dirs;
    my %all_taxids;
    for my $tsv_file (@tsv_files) {
        open my $tsv_fh, "<", $tsv_file or die "Cannot open $tsv_file: $!";
        <$tsv_fh>; # ヘッダー行をスキップ
        while (<$tsv_fh>) {
            chomp;
            my @cols = map { trim($_) } split /\t/;
            my $taxid = $cols[3];
            next unless $taxid;
            $all_taxids{$taxid} = 1;
        }
        close $tsv_fh;
    }

    # 指定ディレクトリ内の全tsvからtaxid取得
    my %selected_taxids;
    my @selected_tsvs = glob(File::Spec->catfile($input_dir, "*.tsv"));
    for my $selected_tsv (@selected_tsvs) {
        open my $mt_fh, "<", $selected_tsv or die "Cannot open $selected_tsv: $!";
        <$mt_fh>; # ヘッダー行をスキップ
        while (<$mt_fh>) {
            chomp;
            my @cols = map { trim($_) } split /\t/;
            my $taxid = $cols[3];
            next unless $taxid;
            $selected_taxids{$taxid} = 1;
        }
        close $mt_fh;
    }

    # ①-②をtaxid_strに
    my @taxid_diff = grep { !$selected_taxids{$_} } keys %all_taxids;
    $taxid_str = join(",", sort @taxid_diff);

    # taxidリストを標準出力に表示
    print "BLASTに渡すtaxidリスト:\n$taxid_str\n";

    # --- 各配列ごとにBLAST ---
    foreach my $seqinfo (@seqs) {
        next unless $output_dir;
        my ($header, $seq) = @$seqinfo;
        my ($id) = $header =~ /^>(\S+)/;
        next unless $id;
        # 一時FASTAファイル作成
        my $tmp_fasta = File::Spec->catfile($output_dir, "${id}.fasta");
        open my $tfh, '>', $tmp_fasta or die "Cannot open $tmp_fasta: $!";
        print $tfh "$header\n$seq\n";
        close $tfh;

        $pm->start and next;
        foreach my $db_info (@dbs) {
            my $db_path   = $db_info->{db};
            my $db_number = $db_info->{number};
            my $out = File::Spec->catfile($output_dir, "${id}_$db_number.tsv");
            my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';
            my $cmd = "blastp -query \"$tmp_fasta\" -db \"$db_path\" -max_target_seqs 10 -evalue 1e-10 -taxids $taxid_str -outfmt \"$outfmt\" -out \"$out\"";
            system($cmd) == 0 or warn "[ERROR] blastp failed for $id DB$db_number\n";
        }
        $pm->finish(0);
    }
    $pm->wait_all_children;

    # --- BLAST結果の集計（BLAST_specify直下のみ） ---
    my %result_cooccurence;
    for my $file (bsd_glob(File::Spec->catfile($output_dir, "*_DB1.tsv"))) {
        if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
            my $id = $1;
            my $all_empty = 1;
            for my $i (1..10) {
                my $f = File::Spec->catfile($output_dir, "${id}_DB${i}.tsv");
                if (-e $f && -s $f) {
                    $all_empty = 0;
                    last;
                }
            }
            if ($all_empty) {
                $result_cooccurence{$id} = 1;
            }
        }
    }
    my $result_cooccurence_path = File::Spec->catfile($output_dir, "result_co-occurence.txt");
    open my $out_cooccurence, '>', $result_cooccurence_path or die "Cannot open $result_cooccurence_path: $!";
    print $out_cooccurence "$_\n" for sort keys %result_cooccurence;
    close $out_cooccurence;
}




sub trim {
    my $s = shift // '';
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    return $s;
}

sub load_ids {
    my ($file) = @_;
    my %hash;
    if (-e $file) {
        open my $fh, '<', $file or die "Cannot open $file: $!";
        while (<$fh>) {
            chomp;
            my $id = trim($_);
            $hash{$id} = 1 if $id;
        }
        close $fh;
    }
    return %hash;
}

sub load_done {
    my ($file) = @_;
    my %hash;
    if (-e $file) {
        open my $fh, '<', $file or die "Cannot open $file: $!";
        while (<$fh>) {
            chomp;
            my $id = trim($_);
            $hash{$id} = 1 if $id;
        }
        close $fh;
    }
    return %hash;
}

sub save_ids {
    my ($hashref, $file) = @_;
    open my $fh, '>', $file or die "Cannot open $file: $!";
    for my $id (keys %$hashref) {
        print $fh trim($id), "\n";
    }
    close $fh;
}

