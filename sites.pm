package sites;

# Reads a file with informative positions, keeps them in a reasonably
# compact, but still fairly idiomatic format.

require Exporter ;
@ISA = qw( Exporter ) ;
@EXPORT_OK = qw( read_sites target_positions ) ;

use strict;

# Reads a file listing informative sites.  Sites must actual sites
# (start == end), with one-based(!) coordinates.
sub read_sites {
    my ($targetfile, $one_based) = @_;
    my ($counter, $interval, %target, @labels, %nlabels, %clabels, $total);
    open TARGETS, $targetfile or die "could not read $targetfile: $?!\n";
    print STDERR "[sites.pm] Reading sites\n";
    while (my $line = <TARGETS>) {
        next if $line =~ /^#/;
        next if $line =~ /^\s+$/;
        $counter++; $interval++; if ($interval == 100_000) {print STDERR "\r$counter"; $interval = 0}
        chomp $line;
        my ($chromosome, $start, $end, $anc_state, $der_state, $lineage) = split /\t/, $line;
        $start++ unless defined $one_based; #make 0-based bed file 1-based
        die "start and end should be equal for _sites_: $chromosome $start $end\n" unless $end == $start;
        $target{$chromosome} = [] unless $target{$chromosome} ;
        unless( exists $nlabels{$lineage} ) {
            $nlabels{$lineage} = scalar @labels ;
            push @labels, $lineage ;
        }
        push @{$target{$chromosome}}, (pack "La1a1C", $start, uc $anc_state, uc $der_state, $nlabels{$lineage});
        $clabels{$lineage}++ ;
        $total++ ;
    }
    while (my ($k,$v) = each %target) {
        my @arr = @$v ;
        undef $target{$k} ;
        @arr = sort {(unpack "L", $a) <=> (unpack "L", $b)} @arr ;
        $target{$k} = \@arr ;
    }
    print STDERR "\n";
    return( \%target, \@labels, \%clabels, $total ) ;
}

# Extracts a hash table of positions overlapping an interval.  
sub target_positions {
    my ($posns, $start, $end) = @_;
    return () unless defined $posns ;

    my $l1 = 0 ;
    my $r1 = scalar @{$posns} ;
    while( $l1 != $r1 ) {
        my $m = ($l1 + $r1) >> 1; 
        my $posm = unpack "L", $posns->[$m] ;
        if( $start > $posm ) { $l1 = $m+1 ; }
        else { $r1 = $m ; }
    }
    my $l2 = 0 ;
    my $r2 = scalar @{$posns} ;
    while( $l2 != $r2 ) {
        my $m = ($l2 + $r2) >> 1;
        my $posm = unpack "L", $posns->[$m] ;
        if( $end > $posm ) { $l2 = $m+1 ; }
        else { $r2 = $m ; }
    }
    my %result;
    for (@{$posns}[$l1 .. $r2-1]) {
        my ($start,@rest) = unpack "La1a1C", $_ ;
        $result{$start} = \@rest;
    }
    return %result;
}

1;
