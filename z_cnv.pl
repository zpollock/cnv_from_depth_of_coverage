#!/usr/bin/perl
use Data::Dumper;

if ($#ARGV < 1 || !($#ARGV % 2)) {
    $setupMsg = "
    Make sure to include one or more pairs of GATK depth of coverage interval summaries as command line parameters. 
    Each pair should consist of the normal sample as the first entry and the tumor sample as the second entry. 
    For example: \"perl ./z_cnv patient_1_N.interval_summary patient_1_T.interval_summary\".";    
    print $setupMsg;
    exit;
}

my @samples = (@ARGV);
my %meta;
my ($annotations, $ann_totals) = loadAnnotations("./EXOME.interval_list.annotated");
my $specificity = 1;

while (my ($normal, $tumor) = splice(@samples, 0, 2)) {
    my @n = loadFile($normal);
    my @t = loadFile($tumor);
    my $normal_avg = $meta{$normal}[2] / $meta{$normal}[3];
    my $tumor_avg = $meta{$tumor}[2] / $meta{$tumor}[3];
    my $normalizer = ($normal_avg / $tumor_avg);
    my $metaName = $meta{$normal}[0] . '-' . $meta{$tumor}[0];
    my $fileOut = $metaName . ".totals";
    
    my @ratio;
    my $sum = 0;
    for my $i (0 .. $#n) {
        if ($n[$i][2] == 0) {
            $ratio[$i] = 1;
        } else {
            $ratio[$i] = ($t[$i][2] * $normalizer) / $n[$i][2];
        }
        $sum += $ratio[$i];
    }

    my $avg = $sum / scalar @ratio;
    my $std = stdev(\@ratio, $avg);

    print "Writing output to $fileOut\n";
    print "Depth of Coverage STD=$std, AVG=$avg\n";

    open(my $fh, '>', $fileOut);
    open(my $fh2, '>', $fileOut . ".cnv");
    print $fh "Samples\tLocation\tGene\tNormalDepth\tTumorNormalizedDepth\tRatio\tZScore\n";
    my $last_gene = "";
    my $gene_count = 0;
    my $amp = 1;
    my $del = 1;
    for my $i (0 .. $#ratio) {
        my $z = ($ratio[$i] - $avg) / $std;
        my $gene = $annotations->{$n[$i][0]};
        print $fh "$metaName\t$n[$i][0]\t$gene\t$n[$i][2]\t" . ($t[$i][2] * $normalizer) . "\t$ratio[$i]\t$z\n";
        if ($last_gene ne $gene) {
            $last_gene = $gene;
            $gene_count = 1;
        } else {
            $gene_count++;
        }
        if ($z >= $specificity) {
            $del = 0;
            $amp &= 1;
        } elsif ($z <= -($specificity)) {
            $del &= 1;
            $amp = 0;
        }
        if ($gene_count == $ann_totals->{$gene}) {
            if ($amp) {
                print $fh2 "$metaName\t$n[$i][0]\t$gene\t$n[$i][2]\t" . ($t[$i][2] * $normalizer) . "\t$ratio[$i]\t$z\tamp\n";
            } elsif ($del) {
                print $fh2 "$metaName\t$n[$i][0]\t$gene\t$n[$i][2]\t" . ($t[$i][2] * $normalizer) . "\t$ratio[$i]\t$z\tdel\n";
            }
        }
    }
    close $fh;
    close $fh2;
}

sub loadFile {
    my $sample = $_[0];
    print "Loading " . $sample . "\n";
    open(my $in, '<', $sample);
    my $header = <$in>;
    my @array;
    my $sum = 0;
    while (<$in>) {
        chomp;
        my @cols = split(/\t/);
        $sum += $cols[2];
        push @array, \@cols;
    }
    close $in;
    my @file_parts = split(/\./, $sample);
    $meta{$sample} = [$file_parts[0], $header, $sum, scalar(@array)];
    return @array;
}

sub loadAnnotations {
    my $file = $_[0];
    my %annotations;
    my %ann_totals;
    open(my $in, '<', $file);
    while (<$in>) {
        chomp;
        my @cols = split(/\t/);
        $annotations{$cols[0]} = $cols[1];
        if (exists $ann_totals{$cols[1]}) {
            $ann_totals{$cols[1]} += 1;
        } else {
            $ann_totals{$cols[1]} = 1;
        }
    }
    close $in;
    return (\%annotations, \%ann_totals);
}

sub stdev {
    my ($data, $average) = @_;
    my $sqtotal = 0;
    foreach (@$data) {
        $sqtotal += ($average - $_) ** 2;
    }
    my $std = ($sqtotal / (@$data - 1)) ** 0.5;
    return $std;
}
