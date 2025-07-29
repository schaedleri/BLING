#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use File::Spec;
use File::Glob ':glob';

# === オプション ===
my $base_dir = '.';
my $dir_input = '';
my $dir_names = '';
my $all_flag = 0;  # 追加
GetOptions(
    'base-dir=s'   => \$base_dir,
    'dir=s'        => \$dir_input,    # 例: 1,2,5
    'dir-names=s'  => \$dir_names,    # 例: Genus_Adlercreutzia,Species_Ecoli
    'all'          => \$all_flag      # 追加: --all で全ディレクトリ対象
) or die "Usage: $0 [--base-dir DIR] [--dir 1,2,3] [--dir-names name1,name2] [--all]\n";

# === パス設定 ===
$base_dir = abs_path($base_dir);
my $script_dir     = "$base_dir/scripts";
my $microbiome_dir = "$base_dir/microbiome";
my $done_root      = "$base_dir/done";
make_path($done_root) unless -d $done_root;

# === ディレクトリ列挙 ===
opendir(my $dh, $microbiome_dir) or die "Cannot open $microbiome_dir: $!";
my @subdirs = sort grep { /^Genus_|^Species_/ && -d "$microbiome_dir/$_" } readdir($dh);
closedir($dh);
die "No target directories found.\n" unless @subdirs;

# === 対象ディレクトリの選択 ===
my @targets;
if ($all_flag) {
    @targets = @subdirs;
}
elsif ($dir_names) {
    my @names = split /,\s*/, $dir_names;
    for my $name (@names) {
        if (grep { $_ eq $name } @subdirs) {
            push @targets, $name;
        } else {
            warn "Warning: '$name' not found in $microbiome_dir\n";
        }
    }
}
elsif ($dir_input) {
    my @indices = split /,\s*/, $dir_input;
    for my $i (@indices) {
        if ($i =~ /^\d+$/ && $i >= 1 && $i <= @subdirs) {
            push @targets, $subdirs[$i - 1];
        } else {
            warn "Invalid index: $i\n";
        }
    }
}
else {
    # 対話モード
    print "Select directories to process (comma-separated index numbers):\n";
    for my $i (0 .. $#subdirs) {
        printf "  [%2d] %s\n", $i + 1, $subdirs[$i];
    }
    print "\nYour choice: ";
    chomp(my $input = <STDIN>);
    my @indices = split /,\s*/, $input;
    for my $i (@indices) {
        if ($i =~ /^\d+$/ && $i >= 1 && $i <= @subdirs) {
            push @targets, $subdirs[$i - 1];
        } else {
            warn "Invalid index: $i\n";
        }
    }
}
die "No valid directories selected.\n" unless @targets;

# === ステップ定義 ===
my @steps = (
    {
        name       => 'cdhit',
        script     => "$script_dir/microbiome_cd-hit.pl",
        depends_on => [],
        mode       => 'single',
    },
    {
        name       => 'blastscreen',
        script     => "$script_dir/microbiome_cdhit_blastscreening.pl",
        depends_on => ['cdhit'],
        mode       => 'single',
    },
    {
        name       => 'externalblast',
        script     => "$script_dir/microbiome_externalblast.pl",
        depends_on => ['blastscreen'],
        mode       => 'single',
    },
);

# === blastscreen.*.doneが存在するディレクトリのみgenus.txt/species.txtを再生成 ===
for my $dir (@targets) {
    my $blast_done = "$done_root/blastscreen.$dir.done";
    next unless -e $blast_done;

    my $cdhit_dir = "$microbiome_dir/$dir/CDhit";
    next unless -d $cdhit_dir;

    if ($dir =~ /^Genus_/) {
        # Genusディレクトリの場合のみgenus.txtを再生成
        my $genus_txt = "$cdhit_dir/genus.txt";
        open my $gf, '>', $genus_txt or die "Cannot write $genus_txt: $!";

        my %unique_ids;
        for my $subdir (glob(File::Spec->catdir($cdhit_dir, "*"))) {
            next unless -d $subdir;
            for my $file (glob(File::Spec->catfile($subdir, "*_DB1.tsv"))) {
                if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                    my $id = $1;
                    my $has_nonzero = 0;
                    for my $i (1..10) {
                        my $f = File::Spec->catfile($subdir, "${id}_DB${i}.tsv");
                        if (-e $f && -s $f) {
                            $has_nonzero = 1;
                            last;
                        }
                    }
                    $unique_ids{$id} = 1 if $has_nonzero;
                }
            }
        }
        print $gf "$_\n" for sort keys %unique_ids;
        close $gf;
    }
    elsif ($dir =~ /^Species_/) {
        # Speciesディレクトリの場合のみspecies.txtを再生成
        my $species_txt = "$cdhit_dir/species.txt";
        open my $sf, '>', $species_txt or die "Cannot write $species_txt: $!";

        my %unique_ids;
        for my $subdir (glob(File::Spec->catdir($cdhit_dir, "*"))) {
            next unless -d $subdir;
            for my $file (glob(File::Spec->catfile($subdir, "*_DB1.tsv"))) {
                if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                    my $id = $1;
                    my $has_nonzero = 0;
                    for my $i (1..10) {
                        my $f = File::Spec->catfile($subdir, "${id}_DB${i}.tsv");
                        if (-e $f && -s $f) {
                            $has_nonzero = 1;
                            last;
                        }
                    }
                    $unique_ids{$id} = 1 if $has_nonzero;
                }
            }
        }
        print $sf "$_\n" for sort keys %unique_ids;
        close $sf;
    }
}

# === ステップごとに実行 ===
for my $step (@steps) {
    print "\n=== STEP: $step->{name} ===\n";

    my @need_run;
    for my $dir (@targets) {
        my $flag = "$done_root/$step->{name}.$dir.done";
        if (-e $flag) {
            print "  ? Skip $step->{name} for $dir (already done)\n";
        } else {
            for my $dep (@{ $step->{depends_on} }) {
                my $dep_flag = "$done_root/$dep.$dir.done";
                die "  ? Cannot run $step->{name} before $dep is done for $dir\n" unless -e $dep_flag;
            }
            push @need_run, $dir;
        }
    }
    next unless @need_run;

    if ($step->{mode} eq 'single') {
        print "  ? Running $step->{name} for: @need_run\n";
        for my $dir (@need_run) {
            my $status = system('perl', $step->{script}, '--base-dir', $base_dir, $dir);
            if ($status != 0) {
                die "  ? Step $step->{name} failed for $dir (exit $status)\n";
            }
            open my $fh, '>', "$done_root/$step->{name}.$dir.done" or die "Cannot write done file\n";
            print $fh "done\n";
            close $fh;
            print "  ? Completed $step->{name} for $dir\n";
        }
    } else {
        die "Unsupported mode: $step->{mode}\n";
    }

    print "  ? Completed $step->{name}\n";
}

print "\n? All selected directories processed successfully.\n";

# === externalblast.*.doneが存在するディレクトリのみresult_genus.txt/result_species.txtおよびgene_list.txtを再生成 ===
for my $dir (@targets) {
    my $external_done = "$done_root/externalblast.$dir.done";
    next unless -e $external_done;

    my $input_dir  = File::Spec->catdir($microbiome_dir, $dir);
    my $output_dir = File::Spec->catdir($input_dir, "BLAST_specify");
    my $fasta_root = File::Spec->catdir($input_dir, "fasta");

    my %result_genus;
    my %result_species;

    # サブディレクトリごとの集計
    for my $subdir (bsd_glob(File::Spec->catdir($output_dir, "*/"))) {
        next unless -d $subdir;
        my @genus_ids;
        for my $file (bsd_glob(File::Spec->catfile($subdir, "genus", "*_DB1.tsv"))) {
            if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                push @genus_ids, $1;
            }
        }
        my @species_ids;
        for my $file (bsd_glob(File::Spec->catfile($subdir, "species", "*_DB1.tsv"))) {
            if ($file =~ /([^\/\\]+)_DB1\.tsv$/) {
                push @species_ids, $1;
            }
        }
        foreach my $id (@genus_ids) {
            my $all_empty = 1;
            for my $i (1..10) {
                my $f = File::Spec->catfile($subdir, "genus", "${id}_DB${i}.tsv");
                if (-e $f && -s $f) {
                    $all_empty = 0;
                    last;
                }
            }
            $result_genus{$id} = 1 if $all_empty;
        }
        foreach my $id (@species_ids) {
            my $all_empty = 1;
            for my $i (1..10) {
                my $f = File::Spec->catfile($subdir, "species", "${id}_DB${i}.tsv");
                if (-e $f && -s $f) {
                    $all_empty = 0;
                    last;
                }
            }
            $result_species{$id} = 1 if $all_empty;
        }
    }

    # クラスターファイルから*付きIDのクラスターメンバーを追記
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

        if (%result_genus) {
            my %add_ids;
            for my $cluster (@clusters) {
                my $star = $cluster->{star_id};
                next unless $star && $result_genus{$star};
                $add_ids{$_} = 1 for @{ $cluster->{ids} };
            }
            $result_genus{$_} = 1 for keys %add_ids;
        }
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

    # BLAST_specify直下に1つだけファイルを出力
    my $result_genus_path   = File::Spec->catfile($output_dir, "result_genus.txt");
    my $result_species_path = File::Spec->catfile($output_dir, "result_species.txt");

    open my $out_genus,   '>', $result_genus_path   or die "Cannot open $result_genus_path: $!";
    open my $out_species, '>', $result_species_path or die "Cannot open $result_species_path: $!";

    print $out_genus   "$_\n" for sort keys %result_genus;
    print $out_species "$_\n" for sort keys %result_species;

    close $out_genus;
    close $out_species;
}

# --- gene_list.txtを作成 ---
my %gene_ids;
for my $dir (@targets) {
    my $input_dir  = File::Spec->catdir($microbiome_dir, $dir);
    my $output_dir = File::Spec->catdir($input_dir, "BLAST_specify");
    my $mode = "";
    if    ($dir =~ /^Genus_/)   { $mode = "genus"; }
    elsif ($dir =~ /^Species_/) { $mode = "species"; }
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