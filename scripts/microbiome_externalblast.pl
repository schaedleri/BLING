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
    # ここに今までの処理（$selected_genusを使う部分）をすべて入れる
    # --- Genus/Species自動判別 ---
    my $mode = "";
    if    ($selected_genus =~ /^Genus_/)   { $mode = "genus"; }
    elsif ($selected_genus =~ /^Species_/) { $mode = "species"; }
    else {
        die "ディレクトリ名はGenus_またはSpecies_で始まる必要があります\n";
    }

    # --- 設定 ---
    my $input_dir         = File::Spec->catdir($base_dir, "microbiome", $selected_genus);
    my $fasta_root        = File::Spec->catdir($input_dir, "fasta");
    my $output_dir        = File::Spec->catdir($input_dir, "BLAST_specify");
    my $used_hits_file    = File::Spec->catfile($output_dir, "used_hits.txt");
    my $used_queries_file = File::Spec->catfile($output_dir, "used_queries.txt");
    my $done_file         = File::Spec->catfile($output_dir, "done_files.txt");
    my $cpu               = 10;
    my $taxid_str;

    # --- taxid情報の読み込み ---
    my @tsv_files = glob(File::Spec->catfile($base_dir, "microbiome", "*", "*.tsv"));
    my (%species2taxid, %genus_taxids);

    for my $tsv_file (@tsv_files) {
        open my $tsv_fh, "<", $tsv_file or die "Cannot open $tsv_file: $!";
        <$tsv_fh>; # ヘッダー行をスキップ
        while (<$tsv_fh>) {
            chomp;
            my @cols = map { trim($_) } split /\t/;
            my $org_name = $cols[2];
            my $taxid    = $cols[3];
            next unless $org_name && $taxid;
            my ($genus, $species) = map { trim($_) } (split(/\s+/, $org_name))[0,1];
            next unless $genus && $species;
            my $species_name = "${genus}_$species";
            $species2taxid{$species_name} = $taxid;
            push @{ $genus_taxids{$genus} }, $taxid;
        }
        close $tsv_fh;
    }

    # --- すべてのtsvからtaxid取得（①） ---
    my @all_dirs = grep { -d $_ && File::Basename::basename($_) !~ /\+/ }
        glob(File::Spec->catdir($base_dir, "microbiome", "*"));
    @tsv_files = map { glob(File::Spec->catfile($_, "*.tsv")) } @all_dirs;
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

    # --- 指定Genusのtsvからtaxid取得（②） ---
    my %selected_taxids;
    my $selected_tsv = File::Spec->catfile($input_dir, "$selected_genus.tsv");
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

    # --- ①-②をtaxid_strに ---
    my @taxid_diff = grep { !$selected_taxids{$_} } keys %all_taxids;
    $taxid_str = join(",", sort @taxid_diff);

    # taxid_strの中身を確認
    #print STDERR "[DEBUG] taxid_str = $taxid_str\n";

    # --- species/genus情報を選択Genusのtsvから取得 ---
    open my $tsv_fh, "<", $selected_tsv or die "Cannot open $selected_tsv: $!";
    <$tsv_fh>; # skip header

    my %species_of_dir;
    my %genus_of_dir;
    while (<$tsv_fh>) {
        chomp;
        my @cols = map { trim($_) } split /\t/;
        my $org_dir = $cols[1]; # ディレクトリ名
        my $org_name = $cols[2]; # 種名
        next unless $org_dir && $org_name;
        my ($genus, $species) = map { trim($_) } (split(/\s+/, $org_name))[0,1];
        next unless $genus && $species;
        $species_of_dir{$org_dir} = "${genus}_$species";
        $genus_of_dir{$org_dir}   = $genus;
    }
    close $tsv_fh;

    my $screening_dir = $output_dir; # 出力先をBLAST_specify直下に

    # --- 静的DBリスト ---
    my @dbs = (
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB1", "BSDB1"),  number => "DB1"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB2", "BSDB2"),  number => "DB2"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB3", "BSDB3"),  number => "DB3"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB4", "BSDB4"),  number => "DB4"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB5", "BSDB5"),  number => "DB5"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB6", "BSDB6"),  number => "DB6"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB7", "BSDB7"),  number => "DB7"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB8", "BSDB8"),  number => "DB8"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB9", "BSDB9"),  number => "DB9"  },
        { db => File::Spec->catfile($base_dir, "DB", "bacteria_strain_taxid_DB", "DB0", "BSDB0"),  number => "DB10" },
    );

    # --- 初期化 ---
    make_path($output_dir);

    # ↓↓↓ done_files.txt関連のロジックを削除 ↓↓↓
    # my %done = load_done($done_file);
    my $pm = Parallel::ForkManager->new($cpu);

    $pm->run_on_finish(
        sub {
            # my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
            # return unless $data;
            # $done{$_} = 1 for @{ $data->{done_keys} };
        }
    );

    $SIG{INT} = sub {
        # save_ids(\%done, $done_file);
        exit 1;
    };

    # --- メイン処理 ---
    # ここからfasta探索方法を変更
    for my $fasta_file (bsd_glob(File::Spec->catfile($fasta_root, "*", "WP_*.fasta"))) {
        my ($org_dir, $id) = (File::Basename::dirname($fasta_file), undef);
        if ($fasta_file =~ /([^\/\\]+)\.fasta$/) {
            $id = $1;
        } else {
            next;
        }
        my $species_dir = File::Basename::basename($org_dir);
        my $base_out = File::Spec->catdir($output_dir, $species_dir);
        my $genus_out   = File::Spec->catdir($base_out, "genus");
        my $species_out = File::Spec->catdir($base_out, "species");
        make_path($genus_out);
        make_path($species_out);

        # --- CDhit直下のIDリストを取得 ---
        my $cdhit_root = File::Spec->catdir($base_dir, "microbiome", $selected_genus, "CDhit");
        my $cdhit_species_txt = File::Spec->catfile($cdhit_root, "species.txt");
        my $cdhit_genus_txt   = File::Spec->catfile($cdhit_root, "genus.txt");

        my %cdhit_species_ids;
        my %cdhit_genus_ids;
        if (-e $cdhit_species_txt) {
            open my $sfh, '<', $cdhit_species_txt or die "Cannot open $cdhit_species_txt: $!";
            while (<$sfh>) {
                chomp;
                $cdhit_species_ids{$_} = 1 if $_;
            }
            close $sfh;
        }
        if (-e $cdhit_genus_txt) {
            open my $gfh, '<', $cdhit_genus_txt or die "Cannot open $cdhit_genus_txt: $!";
            while (<$gfh>) {
                chomp;
                $cdhit_genus_ids{$_} = 1 if $_;
            }
            close $gfh;
        }

        # --- modeによる分岐 ---
        if ($mode eq "genus" && $cdhit_genus_ids{$id}) {
            next unless -s $fasta_file;
            $pm->start and next;
            # my @new_done;
            foreach my $db_info (@dbs) {
                my $db_path   = $db_info->{db};
                my $db_number = $db_info->{number};
                # my $done_key  = "$species_dir:$id:$db_number:genus";
                # next if $done{$done_key};
                my $out = File::Spec->catfile($genus_out, "${id}_$db_number.tsv");
                my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';
                my $cmd = "blastp -query \"$fasta_file\" -db \"$db_path\" -max_target_seqs 10 -evalue 1e-10 -taxids $taxid_str -outfmt \"$outfmt\" -out \"$out\"";
                system($cmd) == 0 or warn "[ERROR] blastp failed for $id DB$db_number\n";
                # push @new_done, $done_key;
            }
            $pm->finish(0); # , { done_keys => \@new_done }
        }
        if ($mode eq "species" && $cdhit_species_ids{$id}) {
            next unless -s $fasta_file;
            $pm->start and next;
            # my @new_done;
            foreach my $db_info (@dbs) {
                my $db_path   = $db_info->{db};
                my $db_number = $db_info->{number};
                # my $done_key  = "$species_dir:$id:$db_number:species";
                # next if $done{$done_key};
                my $out = File::Spec->catfile($species_out, "${id}_$db_number.tsv");
                my $outfmt = '6 qseqid qacc qstart qend qlen sseqid sacc sstart send slen staxids salltitles evalue bitscore qcovs pident ppos qseq sseq';
                my $cmd = "blastp -query \"$fasta_file\" -db \"$db_path\" -max_target_seqs 10 -evalue 1e-10 -taxids $taxid_str -outfmt \"$outfmt\" -out \"$out\"";
                system($cmd) == 0 or warn "[ERROR] blastp failed for $id DB$db_number\n";
                # push @new_done, $done_key;
            }
            $pm->finish(0); # , { done_keys => \@new_done }
        }
    }
    $pm->wait_all_children;
    # save_ids(\%done, $done_file);

    # --- BLAST結果の集計 ---
    my %result_genus;
    my %result_species;

    for my $dir (bsd_glob(File::Spec->catdir($output_dir, "*/"))) {
        next unless -d $dir;
        my $species_dir = basename($dir);

        # genus/speciesディレクトリ内のIDリスト
        my @genus_ids;
        for my $file (bsd_glob(File::Spec->catfile($dir, "genus", "*_DB1.tsv"))) {
            if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                push @genus_ids, $1;
            }
        }
        my @species_ids;
        for my $file (bsd_glob(File::Spec->catfile($dir, "species", "*_DB1.tsv"))) {
            if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                push @species_ids, $1;
            }
        }

        # genus: 10DBすべて空ならresult_genus.txt用に記録
        foreach my $id (@genus_ids) {
            my $all_empty = 1;
            for my $i (1..10) {
                my $f = File::Spec->catfile($dir, "genus", "${id}_DB${i}.tsv");
                if (-e $f && -s $f) {
                    $all_empty = 0;
                    last;
                }
            }
            if ($all_empty) {
                $result_genus{$id} = 1;
            }
        }

        # species: 10DBすべて空ならresult_species.txt用に記録
        foreach my $id (@species_ids) {
            my $all_empty = 1;
            for my $i (1..10) {
                my $f = File::Spec->catfile($dir, "species", "${id}_DB${i}.tsv");
                if (-e $f && -s $f) {
                    $all_empty = 0;
                    last;
                }
            }
            if ($all_empty) {
                $result_species{$id} = 1;
            }
        }
    }

    # --- クラスターファイルから*付きIDのクラスターメンバーを追記 ---
    my ($cluster_file) = glob(File::Spec->catfile($fasta_root, "*.clstr"));
    if ($cluster_file && -e $cluster_file) {
        open my $cfh, '<', $cluster_file or die "Cannot open $cluster_file: $!";
        my @cluster_ids;
        my $star_id;
        my @clusters;
        while (<$cfh>) {
            chomp;
            if (/^>/) {
                push @clusters, { star_id => $star_id, ids => [@cluster_ids] } if @cluster_ids;
                @cluster_ids = ();
                $star_id = undef;
                next;
            }
            if (/\>\s*(WP_\S+)\.\.\.\s*\*$/) {
                $star_id = $1;
                push @cluster_ids, $star_id;
            } elsif (/\>\s*(WP_\S+)\.\.\./) {
                push @cluster_ids, $1;
            }
        }
        push @clusters, { star_id => $star_id, ids => [@cluster_ids] } if @cluster_ids;
        close $cfh;

        # genus
        if (%result_genus) {
            my %add_ids;
            for my $cluster (@clusters) {
                my $star = $cluster->{star_id};
                next unless $star && $result_genus{$star};
                $add_ids{$_} = 1 for @{ $cluster->{ids} };
            }
            $result_genus{$_} = 1 for keys %add_ids;
        }
        # species
        if (%result_species) {
            my %add_ids;
            for my $cluster (@clusters) {
                my $star = $cluster->{star_id};
                next unless $star && $result_species{$star};
                $add_ids{$_} = 1 for @{ $cluster->{ids} };
            }
            $result_species{$_} = 1 for keys %add_ids;
        }
    }

    # --- BLAST_specify直下に1つだけファイルを出力 ---
    my $result_genus_path   = File::Spec->catfile($output_dir, "result_genus.txt");
    my $result_species_path = File::Spec->catfile($output_dir, "result_species.txt");

    open my $out_genus,   '>', $result_genus_path   or die "Cannot open $result_genus_path: $!";
    open my $out_species, '>', $result_species_path or die "Cannot open $result_species_path: $!";

    print $out_genus   "$_\n" for sort keys %result_genus;
    print $out_species "$_\n" for sort keys %result_species;

    close $out_genus;
    close $out_species;
}

# --- ここから下を追加してください ---
# --- 選択したBLAST_specify直下のresult_genus.txt/result_species.txtを統合してgene_list.txtを作成 ---
my %gene_ids;
foreach my $selected_genus (@ARGV) {
    my $input_dir  = File::Spec->catdir($base_dir, "microbiome", $selected_genus);
    my $output_dir = File::Spec->catdir($input_dir, "BLAST_specify");
    my $mode = "";
    if    ($selected_genus =~ /^Genus_/)   { $mode = "genus"; }
    elsif ($selected_genus =~ /^Species_/) { $mode = "species"; }
    else { next; }

    my $result_file = ($mode eq "genus")
        ? File::Spec->catfile($output_dir, "result_genus.txt")
        : File::Spec->catfile($output_dir, "result_species.txt");

    if (-e $result_file) {
        open my $fh, '<', $result_file or die "Cannot open $result_file: $!";
        while (<$fh>) {
            chomp;
            $gene_ids{$_} = 1 if $_;
        }
        close $fh;
    }
}

my $gene_list_path = File::Spec->catfile($base_dir, "gene_list.txt");
open my $gene_fh, '>', $gene_list_path or die "Cannot open $gene_list_path: $!";
print $gene_fh "$_\n" for sort keys %gene_ids;
close $gene_fh;


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

