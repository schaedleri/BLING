#!/usr/bin/perl
use strict;
use warnings;
use Carp;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use Cwd qw(abs_path);

# コマンドライン引数からベースディレクトリとgo.oboパスを受け取る
my $base_dir = '.';
GetOptions(
    'base-dir=s' => \$base_dir,
) or die "Usage: $0 --base-dir <path>\n";

$base_dir = abs_path($base_dir);
my $data_dir="$base_dir/microbiome";

# go.oboのパスは固定
my $go_obo_file = "$base_dir/data/go.obo";

die "Error: go.obo file not found at $go_obo_file\n" unless -e $go_obo_file;

# go.oboの読み込み
my %go_classifications = load_go_classifications($go_obo_file);

# Family_Christensenellaceaeなどのサブディレクトリを取得
my @family_dirs = grep { -d $_ } glob("$data_dir/{Family,Genus,Order,Species}_*");

foreach my $family_dir (@family_dirs) {
    next unless -d $family_dir;
    my $genomic_dir = "$family_dir/genomic_annotations";
    next unless -d $genomic_dir;

    my @genomic_files = glob("$genomic_dir/*.tsv");

    foreach my $genomic_file (@genomic_files) {
        my ($basename) = $genomic_file =~ m{/([^/]+)_genomic\.tsv$};
        next unless $basename;

        print "Processing $basename in $family_dir...\n";

        my $output_dir = "$family_dir/interval_test/$basename";
        make_path($output_dir) unless -d $output_dir;

        my $go_file = "$family_dir/go_annotations/${basename}_go_annotations.tsv";
        my $ip_file = "$family_dir/interpro/${basename}.tsv";

        unless (-e $go_file) {
            warn "File doesn't exist: $go_file\n";
        }
        unless (-e $ip_file) {
            warn "File doesn't exist: $ip_file\n";
        }

        my %go_data = read_go_data($go_file);

        my %ip_data;
        my %ip_go_data;
        if (-e $ip_file) {
            %ip_data = read_ip_data($ip_file, \%go_classifications);
            %ip_go_data = read_ip_go_data($ip_file, \%go_classifications);
        }

        my %genomic_data = read_genomic_data($genomic_file);

        my @clusters = generate_clusters(\%genomic_data);

        open my $fh_all, '>', "$output_dir/interval_afterpsi100.tsv" or die "Can't open file: $!";
        open my $fh_filtered, '>', "$output_dir/filtered_interval_afterpsi100.tsv" or die "Can't open file: $!";
        open my $fh_nearest, '>', "$output_dir/nearest_results100.tsv" or die "Can't open file: $!";
        open my $fh_pairs, '>', "$output_dir/nonGO_GO_pairs100.tsv" or die "Can't open file: $!";
        open my $fh_ipgo, '>', "$output_dir/IP_GO_pairs.tsv" or die "Can't open file: $!";
        print $fh_ipgo "IP_Protein_ID\tGO_ID\tGO_Type\n";

        if (-e $ip_file) {
            output_IP_GO_pairs_from_interpro($ip_file, \%go_classifications, $fh_ipgo);
        }

        print $fh_all join("\t", qw(Accession Organism Begin End Strand Product Gene Locus_Tag Protein_Length Protein_ID Status)) . "\n";
        print $fh_filtered join("\t", qw(Accession Organism Begin End Strand Product Gene Locus_Tag Protein_Length Protein_ID Status)) . "\n";
        print $fh_nearest join("\t", qw(Cluster_ID NonGO_Protein_ID Neighbor_Status Distance Neighbor_Protein_ID)) . "\n";
        print $fh_pairs join("\t", qw(NonGO_Protein_ID GO_Protein_ID GO_Type)) . "\n";

        my $cluster_id = 1;
        foreach my $cluster (@clusters) {
            annotate_cluster($cluster, \%go_data, \%ip_data);
            print_cluster($fh_all, $cluster);
            print_cluster($fh_filtered, $cluster) if is_mixed_cluster($cluster);
            find_nearest($cluster, $cluster_id, $fh_nearest);
            $cluster_id++;
        }

        foreach my $cluster (@clusters) {
            output_nonGO_GO_pairs($cluster, $fh_pairs, \%go_data, \%ip_go_data);
            output_IP_self_GO_pairs($cluster, $fh_ipgo, \%go_data);
        }

        close $fh_all;
        close $fh_filtered;
        close $fh_nearest;
        close $fh_pairs;
        close $fh_ipgo;

        print "Finished processing $basename in $family_dir.\n";
    }
}

# --- サブルーチン ---

sub load_go_classifications {
    my ($file) = @_;
    my %classifications;
    open my $fh, "<", $file or die "Cannot open $file: $!";
    local $/ = "\n\n";  # [Term]単位で読み込み
    while (my $block = <$fh>) {
        if ($block =~ /id: (GO:\d+)/) {
            my $goid = $1;
            my $ns = $block =~ /namespace: (\S+)/ ? $1 : '';
            $classifications{$goid} = $ns;
        }
    }
    close $fh;
    return %classifications;
}

sub read_go_data {
    my ($file) = @_;
    my %go_data;
    return %go_data unless -e $file; # ファイルなければ空hash返す
    open my $fh, '<', $file or die "Can't open $file: $!";
    <$fh>; # ヘッダー読み飛ばし
    while (<$fh>) {
        chomp;
        my ($acc, $org, $go_id, $go_type) = split /\t/;
        if ($go_type eq 'function' || $go_type eq 'process' || $go_type eq 'cellular_component') {
            push @{ $go_data{$acc} }, { id => $go_id, type => $go_type };
        }
    }
    close $fh;
    return %go_data;
}

sub read_ip_data {
    my ($file, $go_class_ref) = @_;
    my %ip_data;
    open my $fh, '<', $file or die "Can't open $file: $!";
    while (<$fh>) {
        chomp;
        my @cols = split /\t/;
        my $prot_id = $cols[0];
        my $n_col = $cols[12];
        my $go_field = $cols[13];
        my @go_terms = $go_field =~ /(GO:\d+)/g;
        my %ns_count;
        foreach my $go_id (@go_terms) {
            my $ns = $go_class_ref->{$go_id} // "unknown";
            $ns_count{$ns}++;
        }
        next if keys(%ns_count) == 1 && exists $ns_count{'cellular_component'};
        if ($n_col ne '-') {
            $ip_data{$prot_id} = 1;
        }
    }
    close $fh;
    return %ip_data;
}

sub read_ip_go_data {
    my ($file, $go_class_ref) = @_;
    my %ip_go_data;
    open my $fh, '<', $file or die "Can't open $file: $!";
    while (<$fh>) {
        chomp;
        my @cols = split /\t/;
        my $prot_id = $cols[0];
        my $go_field = $cols[13];
        my @go_terms = $go_field =~ /(GO:\d+)/g;
        foreach my $go_id (@go_terms) {
            my $ns = $go_class_ref->{$go_id} // "";
            next unless $ns eq 'molecular_function' || $ns eq 'biological_process';
            push @{ $ip_go_data{$prot_id} }, { id => $go_id, type => ($ns eq 'molecular_function' ? 'function' : 'process') };
        }
    }
    close $fh;
    return %ip_go_data;
}

sub read_genomic_data {
    my ($file) = @_;
    my %data;
    open my $fh, '<', $file or die "Can't open $file: $!";
    <$fh>; # ヘッダー読み飛ばし
    while (<$fh>) {
        chomp;
        my ($acc, $org, $begin, $end, $strand, $prod, $gene, $locus, $plen, $pid) = split /\t/;
        push @{ $data{$acc} }, {
            accession => $acc,
            organism => $org,
            begin => $begin,
            end => $end,
            strand => $strand,
            product => $prod,
            gene => $gene,
            locus_tag => $locus,
            protein_length => $plen,
            protein_id => $pid,
        };
    }
    close $fh;
    return %data;
}

sub generate_clusters {
    my ($genomic_ref) = @_;
    my @clusters;
    my $interval = 100;
    for my $acc (keys %$genomic_ref) {
        my @genes = @{ $genomic_ref->{$acc} };
        my @current_cluster = ($genes[0]);
        for my $i (1..$#genes) {
            my $gap = $genes[$i]{begin} - $genes[$i-1]{end};
            my $same_strand = $genes[$i]{strand} eq $genes[$i-1]{strand};
            if ($gap <= $interval && $same_strand) {
                push @current_cluster, $genes[$i];
            } else {
                push @clusters, [@current_cluster] if @current_cluster > 2;
                @current_cluster = ($genes[$i]);
            }
        }
        push @clusters, [@current_cluster] if @current_cluster > 2;
    }
    return @clusters;
}

sub annotate_cluster {
    my ($cluster, $go_data_ref, $ip_data_ref) = @_;
    foreach my $gene (@$cluster) {
        if ($go_data_ref->{$gene->{protein_id}}) {
            $gene->{status} = 'GO';
        } elsif ($ip_data_ref->{$gene->{protein_id}}) {
            $gene->{status} = 'IP';
        } else {
            $gene->{status} = 'nonGO';
        }
    }
}

sub is_mixed_cluster {
    my ($cluster) = @_;
    my ($has_nonGO, $has_other) = (0,0);
    foreach my $gene (@$cluster) {
        if ($gene->{status} eq 'nonGO') {
            $has_nonGO = 1;
        } else {
            $has_other = 1;
        }
        return 1 if $has_nonGO && $has_other;
    }
    return 0;
}

sub print_cluster {
    my ($fh, $cluster) = @_;
    foreach my $gene (@$cluster) {
        print $fh join("\t",
            $gene->{accession},
            $gene->{organism} // '',
            $gene->{begin} // '',
            $gene->{end} // '',
            $gene->{strand} // '',
            $gene->{product} // '',
            $gene->{gene} // '',
            $gene->{locus_tag} // '',
            $gene->{protein_length} // '',
            $gene->{protein_id},
            $gene->{status} // ''
        ), "\n";
    }
    print $fh "\n";
}

sub find_nearest {
    my ($cluster, $cluster_id, $fh) = @_;
    foreach my $gene (@$cluster) {
        next unless $gene->{status} eq 'nonGO';
        foreach my $other (@$cluster) {
            next if $gene == $other;
            next unless $other->{status} eq 'GO' || $other->{status} eq 'IP';
            my $dist;
            if ($gene->{begin} > $other->{end}) {
                $dist = $gene->{begin} - $other->{end};
            } elsif ($gene->{end} < $other->{begin}) {
                $dist = $other->{begin} - $gene->{end};
            } else {
                $dist = 0;
            }
            next if $dist > 300;
            print $fh join("\t",
                $cluster_id,
                $gene->{protein_id},
                "neighbor_" . $other->{status},
                $dist,
                $other->{protein_id}
            ), "\n";
        }
    }
}

sub output_nonGO_GO_pairs {
    my ($cluster, $fh, $go_data_ref, $ip_go_data_ref) = @_;
    foreach my $gene (@$cluster) {
        next unless $gene->{status} eq 'nonGO';
        foreach my $other (@$cluster) {
            next if $gene == $other;
            next unless $other->{status} eq 'GO' || $other->{status} eq 'IP';
            my $dist;
            if ($gene->{begin} > $other->{end}) {
                $dist = $gene->{begin} - $other->{end};
            } elsif ($gene->{end} < $other->{begin}) {
                $dist = $other->{begin} - $gene->{end};
            } else {
                $dist = 0;
            }
            next if $dist > 300;

            if ($other->{status} eq 'GO' && $go_data_ref->{$other->{protein_id}}) {
                for my $go (@{ $go_data_ref->{$other->{protein_id}} }) {
                    next unless $go->{type} eq 'function' || $go->{type} eq 'process';
                    print $fh join("\t",
                        $gene->{protein_id},
                        $go->{id},
                        $go->{type}
                    ), "\n";
                }
            }
            if ($other->{status} eq 'IP' && $ip_go_data_ref->{$other->{protein_id}}) {
                for my $go (@{ $ip_go_data_ref->{$other->{protein_id}} }) {
                    print $fh join("\t",
                        $gene->{protein_id},
                        $go->{id},
                        $go->{type}
                    ), "\n";
                }
            }
        }
    }
}

sub output_IP_self_GO_pairs {
    my ($cluster, $fh, $go_data_ref) = @_;
    foreach my $gene (@$cluster) {
        next unless $gene->{status} eq 'IP';
        if ($go_data_ref->{$gene->{protein_id}}) {
            for my $go (@{ $go_data_ref->{$gene->{protein_id}} }) {
                print $fh join("\t",
                    $gene->{protein_id},
                    $go->{id},
                    $go->{type}
                ), "\n";
            }
        }
    }
}

sub output_IP_GO_pairs_from_interpro {
    my ($ip_file, $go_class_ref, $fh) = @_;
    open my $fh_in, '<', $ip_file or die "Can't open $ip_file: $!";
    while (<$fh_in>) {
        chomp;
        my @cols = split /\t/;
        my $prot_id = $cols[0];
        my $go_field = $cols[13];
        my @go_terms = $go_field =~ /(GO:\d+)/g;
        foreach my $go_id (@go_terms) {
            my $ns = $go_class_ref->{$go_id} // "unknown";
            print $fh join("\t",
                $prot_id,
                $go_id,
                $ns
            ), "\n";
        }
    }
    close $fh_in;
}
