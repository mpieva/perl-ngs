package uarray ;

# "Unboxed array".  This is an array of constant size, 8 byte long
# strings, stored in a string.  The stored strings are occasionally
# unpacked into integers, see below for details.
# 
# Two use cases are supported:  a BED file can be read with read_track()
# and queried for overlap using contains(), or a list of sites can be
# read with read_sites() and queried using target_positions().  Because
# of this, this module can serve as a more compact drop-in replacement
# for both sites.pm and track.pm

require Exporter ;
@ISA = qw(Exporter) ;
@EXPORT_OK = qw( read_track read_sites contains target_positions ) ;

use strict;

sub median3 {
    my ($a,$b,$c) = @_ ;
    if( $a < $b ) {
        if( $b < $c ) {
            return $b ;
        } elsif ( $a < $c ) {
            return $c ;
        } else {
            return $a ;
        }
    } else {
        if( $a < $c ) {
            return $a ;
        } elsif( $b < $c )  {
            return $c ;
        } else {
            return $b ;
        }
    }
}

sub do_usort {
    my ($x, $s, $e) = @_ ;

    while( $e - $s >= 8 ) {
        # Quicksort:  select median of left/middle/right as pivot
        my $l = $s ;
        my $r = $e ;
        my $m = (($s + $e) >> 4) << 3 ;
        my $p = median3 ( (unpack "L", substr($$x,$l,8)) 
                        , (unpack "L", substr($$x,$m,8)) 
                        , (unpack "L", substr($$x,$r,8)) ) ;
        while( $l < $r ) {
            $l += 8 while ( (unpack "L", substr($$x,$l,8)) < $p && $l < $r ) ;
            $r -= 8 while ( (unpack "L", substr($$x,$r,8)) > $p && $l < $r ) ;
            if( $l < $r ) {
                my $t = substr($$x,$l,8) ;
                substr($$x,$l,8, substr($$x,$r,8) ) ;
                substr($$x,$r,8,$t) ;
                $l += 8 ;
                $r -= 8 ;
            }
        }
        if( $l-8-$s < $e-$l ) {
            do_usort( $x, $s, $l-8 ) ;
            $s = $l+8 ;
        } else {
            do_usort( $x, $l+8, $e ) ;
            $e = $l-8 ;
        }
    }
}


sub usort {
    my ($x) = @_ ;
    for( my $i = 0 ; $i+8 < length $$x ; $i+=8 ) {
        if( (unpack "L", substr($$x, $i, 8)) > (unpack "L", substr($$x, $i+8, 8)) ) {
            do_usort( $x, 0, length($$x) - 8 ) ;
            return ;
        }
    }
    return ;
}

# Reads a file listing informative sites.  Sites must actual sites
# (start == end), with one-based(!) coordinates.
sub read_sites($) {
    my ($targetfile, $one_based) = @_;
    my ($counter, $interval, %target, @labels, %nlabels, %clabels, $total);
    open TARGETS, $targetfile or die "could not read $targetfile: $?!\n";
    print STDERR "[uarray.pm] Reading sites\n";
    while (my $line = <TARGETS>) {
        next if $line =~ /^#/;
        next if $line =~ /^\s+$/;
        $counter++; $interval++; if ($interval == 100_000) {print STDERR "\r$counter"; $interval = 0}
        chomp $line;
        my ($chromosome, $start, $end, $anc_state, $der_state, $lineage) = split /\t/, $line;
        die "start and end should be equal for _sites_: $chromosome $start $end\n" unless $end == $start;
        $target{$chromosome} = "" unless $target{$chromosome} ;
        unless( exists $nlabels{$lineage} ) {
            $nlabels{$lineage} = scalar @labels ;
            push @labels, $lineage ;
        }
        # redundant 0 to pad to exactly 8 bytes
        $target{$chromosome} .= pack "La1a1CC", $start-1, uc $anc_state, uc $der_state, $nlabels{$lineage}, 0 ;
        $clabels{$lineage}++ ;
        $total++ ;
    }
    usort( \$target{$_} ) for (keys %target) ;
    print STDERR "\n";
    return( \%target, \@labels, \%clabels, $total ) ;
}

# Extracts a hash table of positions overlapping an interval.  
sub target_positions($$$) {
    my ($posns, $start, $end) = @_;
    return () unless $posns ;

    my $l1 = 0 ;
    my $r1 = length $posns ;
    while( $l1 < $r1 ) {
        my $m = (($l1 + $r1) >> 4) << 3; 
        my $posm = unpack "L", substr($posns,$m,8) ;
        if( $start > $posm ) { $l1 = $m+8 ; }
        else { $r1 = $m ; }
    }
    my $l2 = 0 ;
    my $r2 = length $posns ;
    while( $l2 < $r2 ) {
        my $m = (($l2 + $r2) >> 4) << 3; 
        my $posm = unpack "L", substr($posns,$m,8) ;
        if( $end > $posm ) { $l2 = $m+8 ; }
        else { $r2 = $m ; }
    }
    my %result;
    for( my $i = $l1; $i < $r2; $i+=8 ) {
        my ($start,@rest) = unpack "La1a1C", substr($posns,$i,8) ;
        $result{$start} = \@rest;
    }
    return %result;
}

# Reads a BED file containing regions.  Being a BED file, coordinates
# are zero-based half-open intervals.
sub read_track($) {
    my ($targetfile) = @_;
    my ($counter, $interval, %target);
    open TARGETS, $targetfile or die "could not read $targetfile: $?!\n";
    print STDERR "[uarray.pm] Reading track\n";
    while (my $line = <TARGETS>) {
        next if $line =~ /^#/;
        next if $line =~ /^\s+$/;
        $counter++; $interval++; if ($interval == 100_000) {print STDERR "\r$counter"; $interval = 0}
        chomp $line;
        my ($chromosome, $start, $end) = split /\t/, $line;
        die "target end must be after target start: $chromosome $start $end\nForgot -one?" unless $end >= $start;
        $target{$chromosome} = "" unless $target{$chromosome} ;
        $target{$chromosome} .= pack "LL", $start, $end ;
    }
    usort( \$target{$_} ) for (keys %target) ;
    print STDERR "\n";
    return \%target ;
}

sub contains($$) {
    my ($posns, $pos) = @_;
    return 0 unless $posns ;

    my $l1 = 0 ;
    my $r1 = length $posns ;
    while( $l1 < $r1 ) {
        my $m = (($l1 + $r1) >> 4) << 3; 
        my ($start,$end) = unpack "LL", substr($posns,$m,8) ;
        if( $pos > $end ) { $l1 = $m+8 ; }
        else { $r1 = $m ; }
    }
    return 0 if $l1 == length $posns ;

    my ($start,$end) = unpack "LL", substr($posns,$l1,8) ;
    return ($start <= $pos && $pos < $end)+0 ;
}

1;
