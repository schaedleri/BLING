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
foreach my $genus (@ARGV) {
    my $done_flag = File::Spec->catfile($project_root, "done", "blastscreen.$genus.done");
    if (-e $done_flag) {
        print "[INFO] $done_flag exists. Skipping BLAST, only genus.txt/species.txt regeneration.\n";
        next;
    }

    # --- speciesディレクトリの場合はID抽出のみでBLASTは行わない ---
    if ($genus =~ /^Species_/i) {
        my $species_path = "$base_dir/$genus";
        my $multi_fasta = "$species_path/fasta/${genus}_all_sequences_cdhit";
        my $species_txt = "$species_path/CDhit/species.txt";
        make_path("$species_path/CDhit");

        print "[DEBUG] species.txt path: $species_txt\n";

        open my $in,  '<', $multi_fasta or die "Cannot open $multi_fasta: $!";
        open my $out, '>', $species_txt or die "Cannot open $species_txt: $!";

        my $count = 0;
        while (my $line = <$in>) {
            chomp $line;
            $line =~ s/\r//g;
            if ($line =~ /^>(WP_\d+\.\d+)/) {
                print "[DEBUG] found: $1\n";
                print $out "$1\n";
                $count++;
            }
        }
        close $in;
        close $out;
        print "[INFO] Wrote $count species IDs to $species_txt\n";
        next;
    }

    # --- genusディレクトリは従来通りBLAST等を実行 ---
    my $genus_path = "$base_dir/$genus";
    my $fasta_file = "$genus_path/fasta/${genus}_all_sequences_cdhit";
    my $tsv_file   = "$genus_path/${genus}.tsv";
    my $output_dir = "$genus_path/CDhit";

    # ディレクトリがなければ作成
    make_path($output_dir) unless -d $output_dir;

    my $used_hits_file    = "$output_dir/used_hits.txt";
    my $used_queries_file = "$output_dir/used_queries.txt";
    # ↓↓↓ done_files.txt関連の変数・処理を削除 ↓↓↓
    # my $done_file         = "$output_dir/done_files.txt";
    # for my $file ($used_hits_file, $used_queries_file, $done_file) {
    for my $file ($used_hits_file, $used_queries_file) {
        open my $fh, ">>", $file or die "Cannot open $file: $!";
        close $fh;
    }
    my %used_hits    = load_ids($used_hits_file);
    my %used_queries = load_ids($used_queries_file);
    # my %done         = load_done($done_file);
    my @hit_query_pairs;

    my $pm = Parallel::ForkManager->new($cpu);

    $pm->run_on_finish(sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
        return unless $data;
        $used_hits{$_}    = 1 for @{ $data->{hits} };
        $used_queries{$_} = 1 for @{ $data->{queries} };
        # $done{$_}         = 1 for @{ $data->{done_keys} };
    });

    $SIG{INT} = sub {
        save_ids(\%used_hits,    $used_hits_file);
        save_ids(\%used_queries, $used_queries_file);
        # save_ids(\%done,         $done_file);
        exit 1;
    };

        # taxid情報のハッシュを初期化・セット
    my (%species2taxid, %genus_taxids);
    open my $tsv_fh, "<", $tsv_file or die "Cannot open $tsv_file: $!";
    <$tsv_fh>; # ヘッダー行をスキップ
    while (<$tsv_fh>) {
        chomp;
        my @cols = split /\t/;
        my $org_name = $cols[1];
        my $taxid    = $cols[3];
        next unless $org_name && $taxid;
        my ($genus_name, $species) = (split(/\s+/, $org_name))[0,1];
        next unless $genus_name && $species;
        my $species_name = "${genus_name}_$species";
        $species2taxid{$species_name} = $taxid;
        push @{ $genus_taxids{$genus_name} }, $taxid;
    }
    close $tsv_fh;

    # blastpは常に -task blastp-fast 指定（speciesはスキップ済みのため）
    my $blast_task = 'blastp-fast';

    open my $mf, '<', $fasta_file or die $!;
    local $/ = "\n>";
    while (my $record = <$mf>) {
        $record =~ s/^>//;
        my ($header, @seq_lines) = split /\n/, $record;
        my $seq = join("", @seq_lines);
        $seq =~ s/[^A-Za-z]//g;

        my ($id) = $header =~ /(WP_\d+\.\d+)/ or next;
        my ($org_label) = $header =~ /\[([^\]]+)\]/;
        $org_label //= 'other';
        $org_label =~ s/\s+/_/g;

        my $org_out = "$output_dir/$org_label";
        make_path($org_out);

        next if $used_hits{$id};

        my $tmp_fasta = "$org_out/$id.fasta";
        open my $tf, '>', $tmp_fasta or die $!;
        print $tf ">$header\n$seq\n";
        close $tf;

        $pm->start and next;

        my (@new_hits, @new_queries);

        foreach my $db_info (@dbs) {
            my $db_path   = $db_info->{db};
            my $db_number = $db_info->{number};
            # my $done_key  = "$org_label:$id:$db_number";

            # next if $done{$done_key};

            my $out = "$org_out/${id}_$db_number.tsv";

            my ($genus_name, $species) = split /_/, $org_label;
            my $self_taxid = $species2taxid{"${genus_name}_$species"};
            my @taxids = grep { $_ ne $self_taxid } @{ $genus_taxids{$genus_name} // [] };
            my $taxid_str = join(",", @taxids);

            # taxid_strが空ならBLASTをスキップし、ID抽出のみ(genus.txtに追記)
            if ($taxid_str eq '') {
                warn "[INFO] taxidが空のためBLASTをスキップ: $org_label $id\n";
                my $genus_txt = "$output_dir/genus.txt";
                open my $gfh, '>>', $genus_txt or die "Cannot open $genus_txt: $!";
                flock($gfh, LOCK_EX);
                print $gfh "$id\n";
                flock($gfh, LOCK_UN);
                close $gfh;
                next;
            }

            my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';

            my @cmd = (
                'blastp',
                '-query', $tmp_fasta,
                '-db', $db_path,
                '-task', $blast_task,
                '-max_target_seqs', 10,
                '-evalue', '1e-10',
                '-taxids', $taxid_str,
                '-outfmt', $outfmt,
                '-out', $out
            );

            system(@cmd) == 0 or warn "[ERROR] blastp failed for $id DB$db_number\n";

            # 結果解析
            if (-s $out) {
                open my $fh, '<', $out;
                while (<$fh>) {
                    chomp;
                    my @c = split /\t/;
                    next unless @c >= 16;
                    my $is_good = 0;
                    if ($c[15] >= 95 && $c[14] >= 95) {
                        my ($qstart, $sstart) = ($c[4], $c[9]);
                        if (
                            defined $qstart && defined $sstart &&
                            $qstart =~ /^\d+$/ && $sstart =~ /^\d+$/ &&
                            $sstart > 0
                        ) {
                            my $ratio = $qstart / $sstart;
                            $is_good = 1 if $ratio >= 0.95 && $ratio <= 1.05;
                        }
                    }
                    if ($is_good) {
                        my ($s_clean) = $c[5] =~ /(WP_\d+\.\d+)/;
                        next unless $s_clean;
                        push @new_hits,    $s_clean;
                        push @new_queries, $c[0];
                        push @hit_query_pairs, [$s_clean, $c[0]];
                    }
                }
                close $fh;
                # push @new_done, $done_key;
            }
        }

        # 排他制御付きでusedファイルに追記
        for my $file_data (
            [$used_hits_file,    \@new_hits],
            [$used_queries_file, \@new_queries]
            # [$done_file,         \@new_done]
        ) {
            my ($file, $arr_ref) = @$file_data;
            next unless @$arr_ref;
            open my $fh, '>>', $file or die "Cannot open $file: $!";
            flock($fh, LOCK_EX);
            print $fh "$_\n" for @$arr_ref;
            flock($fh, LOCK_UN);
            close $fh;
        }

        unlink $tmp_fasta;

        $pm->finish(0, {
            hits      => \@new_hits,
            queries   => \@new_queries,
            # done_keys => \@new_done,
        });
    }
    close $mf;

    $pm->wait_all_children;

    # 各ハッシュの内容を表示
    print "=== used_hits ===\n";
    print "$_ => $used_hits{$_}\n" for sort keys %used_hits;

    print "=== used_queries ===\n";
    print "$_ => $used_queries{$_}\n" for sort keys %used_queries;

    # print "=== done ===\n";
    # print "$_ => $done{$_}\n" for sort keys %done;

    # ここでhit_query_pairs.tsvを出力
    write_hit_query_pairs($output_dir);

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
# === genus.txtの自動生成 ===
 foreach my $selected_genus (@ARGV) {
     my $cdhit_root = File::Spec->catdir($base_dir, $selected_genus, "CDhit");
     next unless -d $cdhit_root;
 
     # 既存genus.txtのIDを読み込む
     my %already;
     my $genus_txt_path = File::Spec->catfile($cdhit_root, "genus.txt");
     if (-e $genus_txt_path) {
         open my $gfhr, '<', $genus_txt_path or die "Cannot open $genus_txt_path: $!";
         while (<$gfhr>) {
             chomp;
             $already{$_} = 1 if $_;
         }
         close $gfhr;
     }
 
     my @ids;
     for my $subdir (glob(File::Spec->catdir($cdhit_root, "*"))) {
         next unless -d $subdir;
         for my $file (bsd_glob(File::Spec->catfile($subdir, "*_DB1.tsv"))) {
             if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                 push @ids, { id => $1, dir => $subdir };
             }
         }
    }
 
     open my $out_fh2, '>', $genus_txt_path or die "Cannot open genus.txt: $!";
 
     # 既存IDをまず書き戻す
     print $out_fh2 "$_\n" for sort keys %already;
 
     foreach my $rec (@ids) {
         my $id  = $rec->{id};
         my $dir = $rec->{dir};
         # 既にgenus.txtにあるIDはスキップ
         next if $already{$id};
         my $has_nonzero = 0;
         for my $i (1..10) {
             my $f = File::Spec->catfile($dir, "${id}_DB${i}.tsv");
             if (-e $f && -s $f) {
                $has_nonzero = 1;
                 last;
             }
         }
         print $out_fh2 "$id\n" if $has_nonzero;
     }
 
     close $out_fh2;
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