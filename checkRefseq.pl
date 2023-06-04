use strict;
use warnings;

my $referenceSeq;

while (my $filename = <STDIN>) {
    chomp $filename;

    if ($filename =~ /\.gz$/) {
        open(my $fh, '-|', 'gzip', '-dc', $filename) or die "Cannot open file '$filename': $!";
        my $gatk_command_line_exists = 0;

        while (my $line = <$fh>) {
            chomp $line;

            if ($line =~ /^##GATKCommandLine/) {
                $gatk_command_line_exists = 1;
                if ($line =~ /--reference\s+(\S+)/) {
                    $referenceSeq = $1;
                }
                last;
            }
            last if $line !~ /^##/;
        }

        close($fh);
        my $reference_output = $gatk_command_line_exists ? $referenceSeq : '-';
        print "$filename\t$gatk_command_line_exists\t$reference_output\n";
    }
    else {
        print "$filename\t-\t-\n";
    }
}
