#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Fcntl qw(:flock);
use Parallel::ForkManager;
use Getopt::Long;
use File::Spec;
use File::Glob ':glob';

my $project_root;

my $cpu = 10;

GetOptions(
    'base-dir=s' => \$project_root,
    'cpu=i'      => \$cpu,
) or die "Usage: $0 --base-dir <BLING root path> [--cpu <n>] <Genus_or_Species_dir> ...\n";

die "Error: --base-dir is required\n" unless defined $project_root;
die "Usage: $0 --base-dir <BLING root path> [--cpu <n>] <Genus_or_Species_dir> ...\n" unless @ARGV;

my $base_dir = "$project_root/microbiome";

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

# 複数Genusを+で連結したディレクトリ名を作成
my $genus_dir_name = join("+", @ARGV);
my $target_dir = "$base_dir/$genus_dir_name";


# ARGVが連結名の場合は分割
my @genus_list = map { split /\+/, $_ } @ARGV;

# genusごとにディレクトリとtaxid情報を保持
my %genus_info;
foreach my $genus (@genus_list) {
    my $genus_path = "$base_dir/$genus";
    my $tsv_file = "$genus_path/${genus}.tsv";
    my @taxids;
    if (-e $tsv_file) {
        open my $tfh, '<', $tsv_file or next;
        my $header = <$tfh>;
        if (defined $header) {
            while (<$tfh>) {
                chomp;
                my @cols = split /\t/;
                push @taxids, $cols[3] if defined $cols[3] && $cols[3] ne '';
            }
        }
        close $tfh;
    }
    # 自分自身のtaxidを除去
    @taxids = grep { $_ ne $genus } @taxids;
    $genus_info{$genus} = {
        path   => $genus_path,
        taxids => \@taxids,
    };
}

foreach my $genus (@genus_list) {
    my $genus_path = $genus_info{$genus}{path};
    my @taxids    = @{ $genus_info{$genus}{taxids} };

    unless (@taxids) {
        warn "[INFO] $genus のtsvにtaxidが無いためBLAST処理をスキップします。\n";
        next;
    }

    # ★ すべてのfastaファイルを取得
    my @fasta_files = bsd_glob(File::Spec->catfile($target_dir, "merged_cluster_representative_by_dir_*.fasta"));
    unless (@fasta_files) {
        warn "[INFO] No merged_cluster_representative_by_dir_*.fasta found in $target_dir\n";
        next;
    }

    my $pm = Parallel::ForkManager->new($cpu);

    foreach my $fasta_file (@fasta_files) {
        my $basename = File::Basename::basename($fasta_file, '.fasta');
        my $result_dir = File::Spec->catdir($target_dir, $basename);
        unless (-d $result_dir) {
            make_path($result_dir) or die "Cannot create result dir: $result_dir";
        }
        my @file_genus = $basename =~ /Genus_[A-Za-z0-9]+/g;
        my %file_genus_hash = map { $_ => 1 } @file_genus;

        # 追加：ファイル由来genusが選択genusすべてと一致ならBLAST処理をスキップ
        my $all_included = 1;
        for my $g (@genus_list) {
            $all_included = 0 unless exists $file_genus_hash{$g};
        }
        $all_included = 0 unless scalar(@file_genus) == scalar(@genus_list);

        if ($all_included) {
            warn "[INFO] $basename のfastaは全genus由来のためBLAST処理をスキップします。\n";
            next;
        }

        open my $in, '<', $fasta_file or die "Cannot open $fasta_file: $!";
        my ($seq_id, $seq);
        while (my $line = <$in>) {
            chomp $line;
            if ($line =~ /^>(WP_\d+\.\d+)/) {
                if (defined $seq_id && defined $seq && $seq ne '') {
                    my $pid = $pm->start;
                    if ($pid == 0) { # 子プロセス
                        my $tmp_fasta = "$target_dir/tmp_${seq_id}.fasta";
                        open my $tf, '>', $tmp_fasta or die $!;
                        print $tf ">$seq_id\n$seq\n";
                        close $tf;
                        for my $db (@dbs) {
                            my $db_path   = $db->{db};
                            my $blast_task = 'blastp-fast';
                            foreach my $target_genus (@genus_list) {
                                next if exists $file_genus_hash{$target_genus};
                                # ここでtaxidをgenusごとに取得
                                my @target_taxids = @{ $genus_info{$target_genus}{taxids} };
                                next unless @target_taxids;
                                my $result_file = File::Spec->catfile($result_dir, "result_${basename}_${seq_id}_vs_${target_genus}_$db->{number}.tsv");
                                my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';
                                my @cmd = (
                                    'blastp',
                                    '-query', $tmp_fasta,
                                    '-db', $db_path,
                                    '-task', $blast_task,
                                    '-max_target_seqs', 10,
                                    '-evalue', '1e-10',
                                    '-taxids', join(",", @target_taxids),
                                    '-outfmt', $outfmt,
                                    '-out', $result_file
                                );
                                system(@cmd) == 0 or warn "[ERROR] blastp failed for $basename $seq_id vs $target_genus $db->{number}\n";
                            }
                        }
                        unlink $tmp_fasta;
                        $pm->finish;
                        exit;
                    }
                }
                $seq_id = $1;
                $seq = '';
            } elsif (defined $seq_id) {
                $seq .= $line;
            }
        }
        # 最後の配列も同様に
        if (defined $seq_id && defined $seq && $seq ne '') {
            my $pid = $pm->start;
            if ($pid == 0) { # 子プロセス
                my $tmp_fasta = "$target_dir/tmp_${seq_id}.fasta";
                open my $tf, '>', $tmp_fasta or die $!;
                print $tf ">$seq_id\n$seq\n";
                close $tf;
                for my $db (@dbs) {
                    my $db_path   = $db->{db};
                    my $blast_task = 'blastp-fast';
                    foreach my $target_genus (@genus_list) {
                        next if exists $file_genus_hash{$target_genus};
                        my @target_taxids = @{ $genus_info{$target_genus}{taxids} };
                        next unless @target_taxids;
                        my $result_file = File::Spec->catfile($result_dir, "result_${basename}_${seq_id}_vs_${target_genus}_$db->{number}.tsv");
                        my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';
                        my @cmd = (
                            'blastp',
                            '-query', $tmp_fasta,
                            '-db', $db_path,
                            '-task', $blast_task,
                            '-max_target_seqs', 10,
                            '-evalue', '1e-10',
                            '-taxids', join(",", @target_taxids),
                            '-outfmt', $outfmt,
                            '-out', $result_file
                        );
                        system(@cmd) == 0 or warn "[ERROR] blastp failed for $basename $seq_id vs $target_genus $db->{number}\n";
                    }
                }
                unlink $tmp_fasta;
                $pm->finish;
                exit;
            }
        }
        close $in;
    }
    $pm->wait_all_children;

    # --- 統合処理：Co-occurence.txt作成 ---
    my $cooccur_file = File::Spec->catfile($target_dir, "Co-occurence.txt");

    # すべてのfastaファイルからIDと配列を取得
    my %seqs;
    foreach my $fasta_file (@fasta_files) {
        open my $fh, '<', $fasta_file or next;
        my ($id, $seq);
        while (my $line = <$fh>) {
            chomp $line;
            if ($line =~ /^>(WP_\d+\.\d+)/) {
                $id = $1;
                $seq = '';
            } elsif (defined $id) {
                $seq .= $line;
            }
            if ((($line =~ /^>/) && defined $id && $seq ne '') || (eof($fh) && defined $id && $seq ne '')) {
                $seqs{$id} = $seq;
            }
        }
        close $fh;
    }

    # 1. DBごとにIDを集計（vs_GenusごとにDB1~10のIDを統合）
    my %vs_id_hit; # vs_Genus => { ID => 1 }
    foreach my $vs_genus (@genus_list) {
        my %id_union;
        for my $db (@dbs) {
            my @result_files = bsd_glob(File::Spec->catfile(
                $target_dir, "*", "result_*_vs_${vs_genus}_$db->{number}.tsv"
            ));
            foreach my $rf (@result_files) {
                open my $rfh, '<', $rf or next;
                while (my $line = <$rfh>) {
                    next if $line =~ /^\s*$/;
                    if ($line =~ /^(WP_\d+\.\d+)/) {
                        $id_union{$1} = 1;
                    }
                }
                close $rfh;
            }
        }
        $vs_id_hit{$vs_genus} = \%id_union;
    }

    # 2. すべてのvs属性でヒットしているIDのみ抽出
    my @vs_genus_list = sort keys %vs_id_hit;
    my %common_ids = %{ $vs_id_hit{ $vs_genus_list[0] } || {} };
    for my $vs (@vs_genus_list[1..$#vs_genus_list]) {
        my %next_ids = %{ $vs_id_hit{$vs} || {} };
        for my $id (keys %common_ids) {
            delete $common_ids{$id} unless exists $next_ids{$id};
        }
    }

    # 3. FASTA形式で出力
    open my $cofh, '>', $cooccur_file or die "Cannot open $cooccur_file: $!";
    for my $id (sort keys %common_ids) {
        if (exists $seqs{$id}) {
            print $cofh ">$id\n$seqs{$id}\n";
        }
    }
    close $cofh;
}
# === サブルーチン ===
sub load_ids {
    my $file = shift;
    my %h;
    return %h unless -e $file;
    open my $fh, '<', $file;
    while (<$fh>) {
        chomp;
        if (/WP_\d+\.\d+/) {
            $h{$&} = 1;
        }
    }
    close $fh;
    return %h;
}

sub save_ids {
    my ($h, $file) = @_;
    open my $fh, '>', $file;
    print $fh "$_\n" for sort keys %$h;
    close $fh;
}

sub write_hit_query_pairs {
    my ($cdhit_dir) = @_;
    my @tsv_files = bsd_glob(File::Spec->catfile($cdhit_dir, '*', '*.tsv'));
    my @pairs;
    foreach my $file (@tsv_files) {
        open my $fh, '<', $file or next;
        while (<$fh>) {
            chomp;
            my @cols = split /\t/;
            next unless @cols >= 6;
            my $query = $cols[0];
            my $subject_full = $cols[5];
            my ($subject) = $subject_full =~ /(WP_\d+\.\d+)/;
            next unless $subject;
            push @pairs, [$subject, $query];
        }
        close $fh;
    }
    my $out_file = File::Spec->catfile($cdhit_dir, 'hit_query_pairs.tsv');
    open my $out, '>', $out_file or die "Cannot open $out_file: $!";
    print $out join("\t", @$_), "\n" for @pairs;
    close $out;
}