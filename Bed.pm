package Bed ;

use strict ;
use warnings ;

# 100% pure Perl reader/writer for Bed
# Works on the same structures the Bam reader/writer use.

sub new {
    my $class = shift ;
    my $input = shift ;
    my $bed = { fh => $input, } ;
    return bless( $bed, $class ) ;
}

sub show {
    my( $pipe, $rec ) = @_ ;
    unless( $rec->{flag} & 0x4 ) {
	print {$pipe} join( "\t", (
		$rec->{rname} // '',
		$rec->{pos} // 0,
		($rec->{pos} // 0) + ($rec->{alnlength} // 0),
		$rec->{qname} // '', 
		$rec->{score} // 0, (($rec->{flag} // 0) & 0x10 ? '-' : '+') )
	    ), "\n" ;
    }
}

sub fetch {
    my $bed = shift->{fh} ;
    if( defined( my $l = <$bed> ) ) {
	chomp $l ;
	my @el = split "\t", $l ;
	return {
	    rname => $el[0],
	    pos => $el[1],
	    alnlength => $el[2] - $el[1],
	    qname => $el[3],
	    score => $el[4],
	    flag => ($el[5] eq '-' ? 16 : 0)
	} ;
    }
    return undef ;
}

sub dequeue { return fetch( @_ ) ; }

1 ;
