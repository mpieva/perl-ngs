package track;

# In-memory representation of a BED file, both reasonably compact and
# reasonably idiomatic.

require Exporter ;
@ISA = qw( Exporter ) ;
@EXPORT_OK = qw( read_track contains ) ;

use strict;

# Reads a BED file containing regions.  Being a BED file, coordinates
# are zero-based half-open intervals.
sub read_track {
    my ($targetfile, $one_based) = @_;
    my ($counter, $interval, %target);
    open TARGETS, $targetfile or die "could not read $targetfile: $?!\n";
    print STDERR "[track.pm] Reading track\n";
    while (my $line = <TARGETS>) {
        next if $line =~ /^#/;
        next if $line =~ /^\s+$/;
        $counter++; $interval++; if ($interval == 100_000) {print STDERR "\r$counter"; $interval = 0}
        chomp $line;
        my ($chromosome, $start, $end) = split /\t/, $line;
        $start++ unless defined $one_based; #make 0-based bed file 1-based
        die "target end must be after target start: $chromosome $start $end\nForgot -one?" unless $end >= $start;
        $target{$chromosome} = [] unless $target{$chromosome} ;
        push @{$target{$chromosome}}, (pack "LL", $start, $end);
    }
    while (my ($k,$v) = each %target) {
        my @arr = @$v ;
        undef $target{$k} ;
        @arr = sort {(unpack "L", $a) <=> (unpack "L", $b)} @arr ;
        $target{$k} = \@arr ;
    }
    print STDERR "\n";
    return \%target ;
}

sub contains {
    my ($posns, $pos) = @_;
    return 0 unless defined $posns ;

    my $l1 = 0 ;
    my $r1 = scalar @{$posns} ;
    while( $l1 != $r1 ) {
        my $m = ($l1 + $r1) >> 1; 
        my ($start,$end) = unpack "LL", $posns->[$m] ;
        if( $pos > $end ) { $l1 = $m+1 ; }
        else { $r1 = $m ; }
    }
    return 0 if $l1 == @{$posns} ;

    my ($start,$end) = unpack "LL", $posns->[$l1] ;
    return ($start <= $pos && $pos < $end)+0 ;
}

1;
