package TwobitFile ;

use strict ;
use warnings ;
use constant MAGIC => 0x1A412743 ;

use Fcntl qw ( SEEK_SET ) ;

my %bitMap  =  ('00' => 'T', '01' => 'C', '10' => 'A', '11' => 'G');
our @byteMap = map{ 
    join '', map $bitMap{ $_ }, unpack '(A2)4', unpack 'B8', chr 
} 0 .. 255 ;

sub new {
    my ($class, $file, %argv) = @_ ;
    open my $fd, '<', $file or die "open $file: $!" ;

    my $buf ;
    read $fd, $buf, 4*4 ;
    my( $unpack, $magic, $version, $seqcount ) = ('V', unpack "VVV", $buf) ;
    ( $unpack, $magic, $version, $seqcount ) = ('N', unpack "NNN", $buf) unless $magic == MAGIC ;
    die "$file: not a 2bit file" unless $magic == MAGIC ;
    die "$file: wrong version $version" unless $version == 0 ;

    my %seqs ;
    for my $seqnum (1 .. $seqcount) {
	read $fd, $buf, 1 ;

	my $name ;
	read $fd, $name, (unpack "C", $buf) ;
	read $fd, $buf, 4 ;
	$seqs{$name} = [ unpack $unpack, $buf ] ;
    }

    bless { -file => $file, -fd => $fd, -unp => $unpack, -seqs => \%seqs }, $class ;
}
	
sub _get_seq_record($$) {
    my( $self, $id ) = @_ ;
    my $arr = $self->{-seqs}->{$id} or die "$id: not found in $self->{-file}" ;
    
    if( @{$arr} == 1 ) {
	my $buf ;
	my $fd = $self->{-fd} ;
	my $unp = $self->{-unp} ;

	my $offset = $arr->[0] ;
	seek $fd, $offset, SEEK_SET ;
	read $fd, $buf, 2*4 ;
	my ($dnasize, $nblockcount) = unpack $unp x 2, $buf ;

	read $fd, $buf, 4*(1+2*$nblockcount) ;
	my @nblocks = unpack "$unp*", $buf ;
	my $mblockcount = pop @nblocks ;

	read $fd, $buf, 4*(1+2*$mblockcount) ;
	my @mblocks = unpack "$unp*", $buf ;
	my $reserved = pop @mblocks ;

	my $dnaoffs = tell $fd ;
	@{$arr} = (
	    $dnaoffs, $dnasize,
	    [ @nblocks[0 .. $nblockcount-1] ],
	    [ @nblocks[$nblockcount..$#nblocks] ],
	    [ @mblocks[0 .. $mblockcount-1] ],
	    [ @mblocks[$mblockcount..$#mblocks] ]
	) ;
    }
    return $arr ;
}

sub _for_overlaps {
    my( $start, $len, $starts, $sizes, $op ) = @_ ;
    return unless @{$starts} ;
    _for_overlaps_r( 0, $#$starts, @_ ) ;
}

sub _for_overlaps_r {
    my $l = shift ;
    my $r = shift ;
    my( $start, $len, $starts, $sizes, $op ) = @_ ;

    # leftmost interval too far to the right?
    return if $starts->[$l] >= $start+$len ;

    # rightmost interval too far to the left?
    return if $starts->[$r] + $sizes->[$r] <= $start ;

    # more than one left?  (split and recurse)
    if( $l != $r ) {
	my $m = ($r+$l) >> 1 ;
	_for_overlaps_r( $l,   $m, @_ ) ;
	_for_overlaps_r( $m+1, $r, @_ ) ;
    } else {
	# one interval left, and it does overlap
	# clamp left...
	my $istart = $starts->[$l] ;
	my $ilen   = $sizes->[$l] ;
	if( $istart < $start ) {
	    $ilen -= $start - $istart ;
	    $istart = $start ;
	}
	if( $istart + $ilen > $start + $len ) {
	    $ilen = $start + $len - $istart ;
	}
	&{$op}( $istart, $ilen ) ;
    }
}

sub length {
    my( $self, $name ) = @_ ;
    my $rec = $self->_get_seq_record( $name ) ;
    return $rec->[1] ;
}

sub seq {
    my( $self, $name, $raw_start, $raw_end, $mask ) = @_ ;
    my $rec = $self->_get_seq_record( $name ) ;

    unless( defined $raw_start ) { $raw_start = 1 ; }
    unless( defined $raw_end ) { $raw_end = $rec->[1] ; }

    my( $start, $len, $rc ) = 
    	$raw_start <= $raw_end
	? ( $raw_start-1, $raw_end-$raw_start+1, 0 )
	: ( $raw_end-1, $raw_start-$raw_end+1, 1 ) ;

    my $buf ;
    seek $self->{-fd}, $rec->[0] + ($start >> 2) , SEEK_SET ;
    read $self->{-fd}, $buf, ($len+($start%4)+3) >> 2 ;

    my $DNA = join '', map{ $byteMap[ $_ ] } unpack 'C*', $buf ;
    $DNA = substr( $DNA, $start%4, $len ) ;

    _for_overlaps( $start, $len, $rec->[2], $rec->[3], sub {
	    my($s,$l) = @_ ; substr( $DNA, $s - $start, $l ) = 'N' x $l ;
	} ) ;
    if( $mask && $mask eq 'hardmask' ) {
	_for_overlaps( $start, $len, $rec->[4], $rec->[5], sub {
		my($s,$l) = @_ ; substr( $DNA, $s - $start, $l ) = 'N' x $l ; } ) ;
    }
    elsif( ! $mask || $mask ne 'nomask' ) {
	_for_overlaps( $start, $len, $rec->[4], $rec->[5], sub {
		my($s,$l) = @_ ;
		substr( $DNA, $s - $start, $l ) = lc substr( $DNA, $s - $start, $l ) ; } ) ;
    }
    
    if( $rc ) {
	$DNA = reverse $DNA ;
	$DNA =~ y/ACGTacgt/TGCAtgca/ ;
    }
    $DNA ;
}

sub ids($) { return keys %{$_[0]->{-seqs}} ; }

1;

__END__

=head1 NAME

Twobit -- Fast indexed access to a 2bit file

=head1 SYNOPSIS

  use Twobit;

  my $db      = TwobitFile->new('/path/to/2bit/file');

  my $seq      = $db->seq('CHROMOSOME_I',4_000_000 => 4_100_000);
  my $revseq   = $db->seq('CHROMOSOME_I',4_100_000 => 4_000_000);
  my @ids      = $db->ids;
  my $length   = $db->length('CHROMOSOME_I');


=head1 DESCRIPTION

Twobit provides indexed access to 2bit (DNA) files.  It
provides random access to each sequence entry, and to subsequences
within each entry, allowing you to retrieve portions of very large
sequences without bringing the entire sequence into memory.

When a 2bit file is opened, the module reads the global header and
memoizes the list of sequence entries.  When a sequence is first
accessed, the module reads the entry's header and memoizes it.  The
actual sequence is then read from the file.

The interface is losely modelled after Bio::DB:Fasta, the coordinate
system follows the same conventions.

=head1 OBJECT METHODS

=over 10

=item $db = TwobitFile-E<gt>new($twobit_path)

Create a new TwobitFile object from a 2bit file.  ndicated by
$fasta_path.  If successful, new() will return the database accessor
object, else it dies with an appropriate message.

=item $raw_seq = $db-E<gt>seq($id [,$start [,$stop [,$mask]]])

Return the raw sequence (a string) given an ID and optionally a start
and stop position in the sequence.  If $id doesn't exist, seq() dies.
If $stop is less than $start, then the reverse complement of the
sequence is returned.  If $stop is omitted, the length of sequence $id
is used.  If $start is omitted, it is assumed as 1.  

If $mask is "nomask", the returned sequence is rendered in all capital
letters. If $mask is "hardmask", masked parts of the returned sequence
are replaced by Ns.  otherwise, masked parts are return as lowercase
letters.


=item $length = $db-E<gt>length($id)

Return the length of the indicated sequence.

=item $ids = $db-E<gt>ids()

Return the list of sequence IDs present in the file.

=back

=head1 SEE ALSO

http://genome.ucsc.edu/FAQ/FAQformat.html#format7

=head1 AUTHOR

Udo Stenzel E<lt>udo_stenzel@eva.mpg.deE<gt>.  

Copyright (c) 2011 MPI EVAN.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

