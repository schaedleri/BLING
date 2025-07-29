#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Getopt::Long;
use Cwd qw(abs_path);

# === オプション処理 ===
my $base_dir  = '.';
my $input_dir = '';
GetOptions(
    'base-dir=s'  => \$base_dir,
    'input-dir=s' => \$input_dir,
) or die "Usage: $0 [--base-dir DIR] [--input-dir DIR]\n";

$base_dir  = abs_path($base_dir);
$input_dir = $input_dir ? abs_path($input_dir) : "$base_dir/data/rawgbff";
my $tsv_root = "$base_dir/microbiome";

# --- 出力用サブディレクトリ名 ---
my %subdirs = (
    faa                 => 'faa',
    go_annotations     => 'go_annotations',
    genomic_annotations => 'genomic_annotations',
    fasta              => 'fasta',
);

# --- tsvファイル探索 ---
my @tsv_files;
find(
    sub {
        return unless -f $_;
        return unless $_ =~ /^(Genus|Family|Order|Species)_.+\.tsv$/;
        push @tsv_files, $File::Find::name;
    },
    $tsv_root,
);

die "No TSV files found under $tsv_root\n" unless @tsv_files;

foreach my $tsv_file (@tsv_files) {
    my $dir = dirname($tsv_file);

    # 出力ディレクトリ作成
    foreach my $subdir (values %subdirs) {
        make_path("$dir/$subdir") unless -d "$dir/$subdir";
    }

    open my $fh, '<', $tsv_file or die "Could not open $tsv_file: $!";
    <$fh>; # ヘッダー行スキップ

    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @cols = split /\t/, $line;
        my $gcf = $cols[0];
        next unless $gcf;

        $gcf =~ s/^\s+|\s+$//g;
        my $organism = $cols[1] // 'unknown';
        $organism =~ s/^\s+|\s+$//g;
        $organism =~ s/\s+/_/g;
        $organism = 'unknown' if $organism eq '';

        my @gbff_files = glob("$input_dir/**/ncbi_dataset/data/$gcf/genomic.gbff");
        foreach my $gbff_file (@gbff_files) {
            my $faa_output_file   = "$dir/$subdirs{faa}/$gcf.faa";
            my $go_output_file    = "$dir/$subdirs{go_annotations}/${gcf}_go_annotations.tsv";
            my $tsv_output_file   = "$dir/$subdirs{genomic_annotations}/${gcf}_genomic.tsv";
            my $gcf_fasta_dir     = "$dir/$subdirs{fasta}/$gcf";
            make_path($gcf_fasta_dir) unless -d $gcf_fasta_dir;

            open my $in,       "<", $gbff_file         or die "Cannot open $gbff_file: $!";
            open my $faa_out,  ">", $faa_output_file   or die "Cannot open $faa_output_file: $!";
            open my $go_out,   ">", $go_output_file    or die "Cannot open $go_output_file: $!";
            open my $tsv_out,  ">", $tsv_output_file   or die "Cannot open $tsv_output_file: $!";

            print $go_out join("\t", "Accession", "Organism", "GO_ID", "GO_Type"), "\n";
            print $tsv_out join("\t", qw(
                Accession Organism Begin End Strand Product Gene Locus_Tag Protein_Length Protein_ID
            )), "\n";

            local $/ = "//\n";  # GenBankエントリの区切り

            while (my $entry = <$in>) {
                my ($accession) = $entry =~ /\sACCESSION\s+(\S+)/;
                my ($org)       = $entry =~ /  ORGANISM +([^\n]+)/;
                next unless $accession && $org;
                my @annotations;

                while ($entry =~ /     CDS             ([^\n]+)\n(.*?)\n(?= {5}\S|\Z)/sg) {
                    my ($location, $cds_block) = ($1, $2);
                    my ($begin, $end, $strand) = parse_location_info($location);
                    my ($protein_id)  = $cds_block =~ /\/protein_id="([^"]+)"/;
                    next unless $protein_id;

                    my ($product_raw) = $cds_block =~ /\/product="((?:[^"]|"\n)+?)"/;
                    my ($gene)        = $cds_block =~ /\/gene="([^"]+)"/;
                    my ($locus_tag)   = $cds_block =~ /\/locus_tag="([^"]+)"/;
                    my ($translation) = $cds_block =~ /\/translation="((?:[^"]|"\n)+?)"/;

                    $product_raw =~ s/"//g if $product_raw;
                    $product_raw =~ s/\n/ /g if $product_raw;
                    $product_raw =~ s/^\s+|\s+$//g if $product_raw;
                    $product_raw =~ s/\s{2,}/ /g if $product_raw;
                    my $product = $product_raw // 'unknown';
                    $gene ||= "unknown";
                    $locus_tag ||= "unknown";
                    my $protein_length = $translation ? length($translation =~ s/\s+//gr) : 0;

                    push @annotations, {
                        accession      => $accession,
                        organism       => $org,
                        begin          => $begin,
                        end            => $end,
                        strand         => $strand,
                        product        => $product,
                        gene           => $gene,
                        locus_tag      => $locus_tag,
                        protein_length => $protein_length,
                        protein_id     => $protein_id,
                    };

                    my $has_go = 0;
                    while ($cds_block =~ /\/GO_(function|component|process)\s*=\s*"([^"]+)"/g) {
                        my ($go_type, $go_value) = ($1, $2);
                        if ($go_value =~ /(GO:\d+)/) {
                            print $go_out join("\t", $protein_id, $org, $1, $go_type), "\n";
                            $has_go = 1 if $go_type =~ /^(function|process)$/;
                        }
                    }

                    # GOアノテーション有無に関わらずfaa出力
                    if ($translation) {
                        $translation =~ s/\s+//g;
                        print $faa_out ">$protein_id [$org] $product\n$translation\n";
                    }

                    if ($translation) {
                        $translation =~ s/\s+//g;
                        my $fasta_file = "$gcf_fasta_dir/$protein_id.fasta";
                        open my $fasta_out, ">", $fasta_file or die "Cannot open $fasta_file: $!";
                        print $fasta_out ">$protein_id [$org] $product\n$translation\n";
                        close $fasta_out;
                    }
                }

                foreach my $ann (@annotations) {
                    print $tsv_out join("\t",
                        map { $_ // 'unknown' } (
                            $ann->{accession},
                            $ann->{organism},
                            $ann->{begin},
                            $ann->{end},
                            $ann->{strand},
                            $ann->{product},
                            $ann->{gene},
                            $ann->{locus_tag},
                            $ann->{protein_length},
                            $ann->{protein_id}
                        )
                    ), "\n";
                }
            }

            close $in;
            close $faa_out;
            close $go_out;
            close $tsv_out;

            print "Processed: $gbff_file -> $faa_output_file, $go_output_file, $tsv_output_file\n";
        }
    }

    close $fh;
}

print "All files processed.\n";

# --- location情報パース ---
sub parse_location_info {
    my ($location) = @_;
    my ($strand, $begin, $end) = ('plus', '', '');

    if ($location =~ /^complement\((.+)\)$/) {
        $strand = 'minus';
        $location = $1;
    }

    if ($location =~ /^join\((.+)\)$/) {
        my @coords = $location =~ /<?(\d+)\.\.>?\d+/g;
        my @ends   = $location =~ /<?\d+\.\.>?(\d+)/g;
        $begin = $coords[0] if @coords;
        $end   = $ends[-1]  if @ends;
    } elsif ($location =~ /<?(\d+)\.\.>?(\d+)/) {
        ($begin, $end) = ($1, $2);
    }

    return ($begin, $end, $strand);
}
