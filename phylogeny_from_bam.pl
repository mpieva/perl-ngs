#!/usr/bin/perl -w
#
# A restricted parsimony method to place an ancient sample into a known
# phylogeny.

use strict;
use Getopt::Long;
use File::Basename;

use lib dirname($0);
use uarray qw( read_sites read_track contains target_positions ) ;
use Bam qw( cigarOp cigarNum ) ;

my ($help, $informative_sites, $sample_sites, $exclude_damage, $track, $bq_filter, $relaxed, $minlength, $maxlength);
my $mapqual = 0;
GetOptions( 
    "informative=s" => \$informative_sites,
    "track=s" => \$track,
    "exclude" => \$exclude_damage,
    "relaxed" => \$relaxed,
    "basequal=i" => \$bq_filter,
    "mapqual=i" => \$mapqual,
    "minlength=i" => \$minlength,
    "maxlength=i" => \$maxlength,
    "help" => \$help);

my $bamfile = shift @ARGV or &help;
&help unless $informative_sites;

my ($sites, $labels, $lineages, $total) = read_sites $informative_sites ;

if ($track) {
    print STDERR "reading track $track\n";
    $track = read_track $track ;
}

my %stats = ( "informative_sites" => $total ) ;

my $bam = Bam::In->new( $bamfile ) ;
print STDERR "reading bam file...\n";
my ($counter, $interval) = (0, 0);
while (my $rec = $bam->fetch()) {
    $interval++; if ($interval == 100_000) {print STDERR "\r$rec->{rname}\t$rec->{pos}\t\t"; $interval = 0;}
    next if defined $mapqual && $rec->{mapq} < $mapqual ;
    next if defined $minlength && length $rec->{seq} < $minlength ; 
    next if defined $maxlength && length $rec->{seq} > $maxlength ;

    die "no MD field seen in\n $rec->{qname}\n" unless $rec->{opt_fields}->{MD};
    my $end = $rec->{pos} + ($rec->{cigar} ? $rec->{alnlength} : length $rec->{seq}) ;
    my %target_positions = target_positions($sites->{$rec->{rname}},$rec->{pos},$end);
    next unless %target_positions ;
    my $genomic_position = $rec->{pos}-1;   # so we can increment first, then process
    for (@{$rec->{cigar}}) {
        my $number = cigarNum $_ ;
        my $type = cigarOp $_ ;

        next if $number == 0;
        if ($type eq 'M') {
            foreach (1..$number) {
                $genomic_position++;
                $rec->{seq} =~ s/^(.)// or die "parsing error in cigar field\n";
                my $base = $1;
                my $basequal = ord($rec->{qual}) ;
                $rec->{qual} = substr($rec->{qual},1);

                if (exists $target_positions{$genomic_position}) {
                    $stats{"match"}++;
                    next if $track && !contains($track->{$rec->{rname}}, $genomic_position) ;
                    $stats{"track"}++;
                    my ($ancestral_base, $derived_base, $lineage_sym) = @{$target_positions{$genomic_position}} ;
                    my $lineage = $labels->[$lineage_sym] ;

                    next unless ($base eq $ancestral_base || $base eq $derived_base) ;
                    $stats{"state_passed"}++;
                    next if $bq_filter && $basequal < $bq_filter;
                    $stats{"BQpassed"}++;
                    if ($relaxed) {
                        unless (isReversed $rec) {
                            next if $base =~ /T/;
                        } else {
                            next if $base =~ /A/;
                        }
                    } elsif ($exclude_damage) {
                        unless (isReversed $rec) {
                            next if $ancestral_base =~ /C/;
                            next if $derived_base =~ /C/;
                        } else {
                            next if $ancestral_base =~ /G/;
                            next if $derived_base =~ /G/;
                        }
                    }
                    $stats{"damage_passed"}++;
                    if ($base eq $ancestral_base) {
                        $stats{"counts"}{$lineage}{"anc"}++;
                    } elsif ($base eq $derived_base) {
                        $stats{"counts"}{$lineage}{"der"}++;
                    }
                    $stats{"final"}++;
                }
            }
        } elsif ($type eq 'I') {
            foreach (1..$number) {
                $rec->{seq} =~ s/^.// or die "parsing error2 in cigar field\n";
                $rec->{qual} = substr($rec->{qual},1);
            }
        } elsif ($type eq 'D') {
            $genomic_position += $number ;
        }
    }
}

print "#lineage\tsample_anc\tsample_der\tsample_total\tfraction_der\n";
foreach my $lineage (sort @{$labels}) {
    print "$lineage\t";
    my $anc = $stats{"counts"}{$lineage}{"anc"} || 0;
    my $der = $stats{"counts"}{$lineage}{"der"} || 0;
    my $total = $anc + $der;
    print "$anc\t$der\t$total\t";
    my $percent_der;
    if ($total) {
        printf "%.2f", $der / $total;
        print "\n";
    } else {
        print "NA\n";
    }
}

print "#informative sites\t", $stats{"informative_sites"} || 0, "\n";
print "#overlapping sequences\t", $stats{"match"} || 0, "\n";
print "#in track\t", $stats{"track"} || 0, "\n";
print "#ancestral/derived state\t", $stats{"state_passed"} || 0, "\n";
print "#base quality passed\t", $stats{"BQpassed"} || 0, "\n";
print "#damage passed\t", $stats{"damage_passed"} || 0, "\n";
print "#final\t", $stats{"final"} || 0, "\n";



sub help {
    print "
    This script reads a table with phylogenetically informative positions (tab-delimited: chromosome, start, end, ancestral_state, derived_state, lineage) and a bam file and outputs the fraction of derived alleles in support of each lineage. If a position is covered more than once, it is counted more than once. 

    [usage] 
    phylogenetic_position.pl [-options] in.bam

    [options]
    -informative      table with informative sites (tab delimited, 1-based) [required]
    -track            specify mappability track (chr start end; 0-based) and disregard all positions outside
    -exclude          ignore sequence on plus strand of the informative alleles is C (and on minus strand if G) 
    -relaxed          ignore sequence with T on plus strand and A on minus strand
    -basequal         ignore sequence with basequal < X (phred scale) at informative site [default 0]
    -mapqual          require mapqual [default 0]
    -minlength        ignore sequences shorter than X
    -maxlength        ignore sequences longer than X\n";
    exit;
}
