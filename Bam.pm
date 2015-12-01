# 100% pure Perl reader & writer for BAM
# (Why?  Because I can, and because a dependency on libbam isn't really
# helpful here.  Also because reading the samtools source hurt.)
#
# TODO:  - test support of B-type fields
#          There is support for optional fields of array type now, but
#          it's untested.  Testing is a bit hard right now, since there
#          are no tools that actually produce such fields.
#        - access the index (BAI or CSI) for selective reading of
#          regions

package Bam ;

require Exporter ;
@ISA = qw( Exporter ) ;
@EXPORT_OK = qw( qual_to_fastq qual_from_fastq cigarOp cigarNum mkCigar) ;

sub qual_to_fastq($) {
    my $s = shift ;
    for( my $i = 0 ; $i != length $s ; ++$i )
    {
        substr( $s, $i, 1 ) = chr( ord( substr $s, $i, 1 ) + 33 ) ;
    }
    return $s ;
}

sub qual_from_fastq($) {
    my $s = shift ;
    for( my $i = 0 ; $i != length $s ; ++$i )
    {
        substr( $s, $i, 1 ) = chr( ord( substr $s, $i, 1 ) - 33 ) ;
    }
    return $s ;
}

sub cigarOp($)  { substr( 'MIDNSHPM', $_[0] & 0x7, 1 ) ; }
sub cigarNum($) { $_[0] >> 4 ; }

sub mkCigar($$) { 
    my ($num,$op) = @_ ;
    ( $op == 'M' ? 0 : $op == 'I' ? 1 :
      $op == 'D' ? 2 : $op == 'N' ? 3 :
      $op == 'S' ? 4 : $op == 'H' ? 5 :
      $op == 'P' ? 6 : 7 ) | ($num << 4) ;
}
      
1;

package Bam::In ;

use strict ;
use warnings ;
use IO::Uncompress::Gunzip ;

sub unpack_array {
    my ($c,$n,$s) = @_ ;
    my @dreck = unpack "${c}${n}a*", $s ;
    my $rest = pop @dreck ;
    return (\@dreck, $rest) ;
}

# Map from codes to procedures to extract optional fields.  Given the
# raw binary code, should return the field value and the remainder.
# NOTE:  - float (code 'f') may not work, since it relies on the
#          machine's representation of floating point numbers  
#        - double (code 'd') is supported when reading , but will not be
#          distinguishable from float
#        - code 'H' is ill-specified, apparently it's the same as Z with a hint
#          to print it in hex in SAM.  we treat it just like Z.

my %xxtype = ( 
    'c' => sub { unpack_array( "c", @_ ) },
    'C' => sub { unpack_array( "C", @_ ) },
    'S' => sub { unpack_array( "v", @_ ) },
    'I' => sub { unpack_array( "V", @_ ) },
    'f' => sub { unpack_array( "f", @_ ) },            # inherently non-portable
    'd' => sub { unpack_array( "d", @_ ) },            # inherently non-portable
    's' => sub { my ($x,$y) = unpack_array( "v", @_ ) ;    # quite portable
                 return ( [ map { $_ >= 0x8000 ? $_-0x10000 : $_ } @$x ], $y ) ; },
    'i' => sub { my ($x,$y) = unpack_array( "v", 2*$_[0], $_[1] ) ;    # quite portable
                 return ( [ map { ($x->[$_+1] >= 0x8000 ? $x->[$_+1]-0x10000 : $x->[$_+1]) << 16 | $x->[$_] }
                                0 .. $_[0]-1 ], $y ) ; } ) ;

my %xtype = ( 
    'c' => sub { unpack "ca*", $_[0] },
    'C' => sub { unpack "Ca*", $_[0] },
    'S' => sub { unpack "va*", $_[0] },
    'I' => sub { unpack "Va*", $_[0] },
    'A' => sub { unpack "aa*", $_[0] },
    'f' => sub { unpack "fa*", $_[0] },            # inherently non-portable
    'd' => sub { unpack "da*", $_[0] },            # inherently non-portable
    'Z' => sub { unpack "Z*a*", $_[0] },           # unclear if this is right
    'H' => sub { unpack "Z*a*", $_[0] },           # same encoding as Z
    's' => sub { my ($x,$y) = unpack "va*", $_[0] ;        # quite portable
                 return ( $x >= 0x8000 ? $x-0x10000 : $x, $y ) ; },
    'i' => sub { my ($x,$y,$z) = unpack "vva*", $_[0] ;    # quite portable
                 return ( ($y >= 0x8000 ? $y-0x10000 : $y) << 16 | $x, $z ) ; },
    'B' => sub { my($c,$n,$rest) = unpack "aVa*", $_[0] ; 
                 if( exists $xxtype{$c} ) { $xxtype{$c}($n,$rest) ; }
                 else { die "unknown array type $c" } } ) ;

sub all { $_ || return 0 for @_; 1 }
sub notall { $_ || return 1 for @_; 0 }

sub get_code {
    if( notall( map { $_ + 0 eq $_ } @_ ) )
    {
        # Not all numbers.  Single string is good, else it's an error.
        @_ == 1 or die "unsuitable array in opt. field" ;
        return ('A','A') if length $_[0] == 1 ;
        return ('Z','Z') ;
    }
    return ('f','f') if( notall( map { $_ == round($_) } @_ ) ) ;            # not all integers, so floats
    return ('c','c') if( all( map {   -0x80 <= $_ && $_ <=   0x7f } @_ ) ) ;
    return ('C','C') if( all( map {       0 <= $_ && $_ <=   0xff } @_ ) ) ;
    return ('s','s') if( all( map { -0x8000 <= $_ && $_ <= 0x7fff } @_ ) ) ; # is this portable?
    return ('v','S') if( all( map {       0 <= $_ && $_ <= 0xffff } @_ ) ) ;
    return ('V','I') if( all( map {       0 <= $_ } @_ ) ) ;
    return ('l','i') ;                                                       # is this portable?
}

sub ensure($$) {
    my ($bam, $len) = @_ ;
    while(length $bam->{buf} < $len) {
        my $buf ;
        if( !read( $bam->{fh}, $buf, 4096 ) ) {
            # no more data, @ret cannot change anymore
            return 0 ;
        }
        $bam->{buf} .= $buf ;
    }
    return 1;
}
	
sub get_int_32($) {
    my ($bam, $len) = @_ ;
    ensure( $bam, 4 ) ;
    my( $x, $s ) = unpack "Va*", $bam->{buf} ;
    $bam->{buf} = $s ;
    return $x ;
}
sub get_string($$) {
    my ($bam, $len) = @_ ;
    ensure( $bam, $len ) ;
    my $s = substr( $bam->{buf}, 0, $len ) ;
    $bam->{buf} = substr( $bam->{buf}, $len ) ;
    return $s ;
}
sub get_string_nul($$) {
    my ($bam, $len) = @_ ;
    ensure( $bam, $len ) ;
    my $s = substr( $bam->{buf}, 0, $len-1 ) ;
    $bam->{buf} = substr( $bam->{buf}, $len ) ;
    return $s ;
}

sub get_aln_length {
    my %good = ( 0 => 1, 2 => 1, 3 => 1 ) ;
    my $len = 0 ;
    for( @_ ) {
        $len += $_ >> 4 if $good{$_ & 0xf} ;
    }
    return $len ;
}

sub new {
    my $class = shift ;
    my $input = shift ;
    my $bam = {
        fh => IO::Uncompress::Gunzip->new( $input, MultiStream => 1 ),
        buf => '',
        refs => []
    } ;

    die "not enough magic" if get_string( $bam, 4 ) ne "BAM\1" ;
    $bam->{header} = get_string( $bam, get_int_32( $bam ) ) ;

    my $nrefs = get_int_32( $bam ) ;
    for (1..$nrefs)
    {
        push @{$bam->{refs}}, [
            get_string_nul( $bam, get_int_32( $bam ) ),
            get_int_32( $bam ) ] ;
    }
    return bless( $bam, $class ) ;
}

sub to_base($) {
    return substr( "=AC.G...T......N", $_[0], 1 ) ;
}

sub fetch {
    my $bam = shift ;
    return undef if( !($bam->{buf}) && (eof $bam->{fh}) ) ;

    my( $block ) = get_string( $bam, get_int_32( $bam ) ) ;
    return undef unless defined $block ;

    my( $rid, $pos, $bin_mq_nl, $flag_nc, $read_len, $mate_rID,
	$mate_pos, $ins_size, $rest ) = unpack "VVVVVVVVa*", $block ;

    my $read_name_len = ($bin_mq_nl & 0xff) - 1 ;
    my $cigar_len = 4 * ($flag_nc & 0xffff) ;
    my $seq_len = int(( $read_len+1 )/2) ;

    my( $read_name, $cigar, $seq, $qual, $tags ) = unpack(
        "a[$read_name_len] x a[$cigar_len] a[$seq_len] ".
        "a[$read_len] a*", $rest ) ;

    $seq = substr( join( '', 
            map( {  to_base( ($_ >> 4) & 0xf ) . 
                to_base( $_ & 0xf ) } 
            unpack( 'C*', $seq ) ) ), 
        0, $read_len ) ;

    my %opt_fields :shared ;
    while( $tags ) {
        my( $tag, $type, $more ) = unpack "a2aa*", $tags ;
        if( exists $xtype{$type} ) {
            ( $opt_fields{$tag}, $tags ) = $xtype{$type}( $more ) ;
        } else {
            die "unknown optional field type '$type.'" ;
        }
    }

    # NOTE: 
    # - pos is 0 based
    # - qual is 0 based (i.e. Phred scale, but no FastQ encoding!)
    # - cigar is 'native' (as in BAM, a list of integers)
    
    my %aln :shared = (
        qname => $read_name,
        flag => $flag_nc >> 16 ,
        seq => $seq,
        qual => $qual,
        opt_fields => \%opt_fields
    ) ;

    unless( $flag_nc & 0x40000 ) {
        my @cig :shared = unpack "V*", $cigar ;
        %aln = (
            %aln,
            mapq => ($bin_mq_nl >> 8) & 0xff,
            cigar => \@cig,
            pos => $pos,
            rname => $bam->{refs}->[$rid]->[0],
            alnlength => get_aln_length( @cig ) 
        ) ;
    }

    $aln{mrnm} = $bam->{refs}->[$mate_rID]->[0] if $mate_rID != 0xffffffff ;
    $aln{mpos} = $mate_pos if $mate_pos != 0xffffffff ;
    $aln{isize} = $ins_size if $ins_size ;

    if( $opt_fields{UQ} ) {
        $aln{score} = $opt_fields{UQ} ;
    }
    elsif( $opt_fields{AS} ) {
        $aln{score} = $opt_fields{AS} ;
    }

    return \%aln ;
}

sub dequeue { return fetch( @_ ) ; }

sub isPaired($)         { $_[0]->{flag} &    1 ; }
sub isProperlyPaired($) { $_[0]->{flag} &    2 ; }
sub isUnmapped($)       { $_[0]->{flag} &    4 ; }
sub isMateUnmapped($)   { $_[0]->{flag} &    8 ; }
sub isReversed($)       { $_[0]->{flag} &   16 ; }
sub isMateReversed($)   { $_[0]->{flag} &   32 ; }
sub isFirstMate($)      { $_[0]->{flag} &   64 ; }
sub isSecondMate($)     { $_[0]->{flag} &  128 ; }
sub isSecondary($)      { $_[0]->{flag} &  256 ; }
sub isLowQuality($)     { $_[0]->{flag} &  512 ; }
sub isDuplicate($)      { $_[0]->{flag} & 1024 ; }
sub isSupplementary($)  { $_[0]->{flag} & 2048 ; }

sub setPaired($)         { $_[0]->{flag} |=    1 ; }
sub setProperlyPaired($) { $_[0]->{flag} |=    2 ; }
sub setUnmapped($)       { $_[0]->{flag} |=    4 ; }
sub setMateUnmapped($)   { $_[0]->{flag} |=    8 ; }
sub setReversed($)       { $_[0]->{flag} |=   16 ; }
sub setMateReversed($)   { $_[0]->{flag} |=   32 ; }
sub setFirstMate($)      { $_[0]->{flag} |=   64 ; }
sub setSecondMate($)     { $_[0]->{flag} |=  128 ; }
sub setSecondary($)      { $_[0]->{flag} |=  256 ; }
sub setLowQuality($)     { $_[0]->{flag} |=  512 ; }
sub setDuplicate($)      { $_[0]->{flag} |= 1024 ; }
sub setSupplementary($)  { $_[0]->{flag} |= 2048 ; }

sub clearPaired($)         { $_[0]->{flag} &=    ~1 ; }
sub clearProperlyPaired($) { $_[0]->{flag} &=    ~2 ; }
sub clearUnmapped($)       { $_[0]->{flag} &=    ~4 ; }
sub clearMateUnmapped($)   { $_[0]->{flag} &=    ~8 ; }
sub clearReversed($)       { $_[0]->{flag} &=   ~16 ; }
sub clearMateReversed($)   { $_[0]->{flag} &=   ~32 ; }
sub clearFirstMate($)      { $_[0]->{flag} &=   ~64 ; }
sub clearSecondMate($)     { $_[0]->{flag} &=  ~128 ; }
sub clearSecondary($)      { $_[0]->{flag} &=  ~256 ; }
sub clearLowQuality($)     { $_[0]->{flag} &=  ~512 ; }
sub clearDuplicate($)      { $_[0]->{flag} &= ~1024 ; }
sub clearSupplementary($)  { $_[0]->{flag} &= ~2048 ; }

1;

package Bam::Out ;

use strict ;
use warnings ;
use IO::Compress::Gzip qw( gzip ) ;

sub new {
    my( $class, $fh, $refs, $header ) = @_ ;

    my $self = {
        fh => $fh,
        refs => {},
        buf => join '',
                pack( "a4 V/a* V", "BAM\1", $header, scalar @$refs ),
                map { pack "V/a* V", $_->[0]."\0", $_->[1] } @$refs 
    } ;
    for (0..$#$refs) {
        $self->{refs}->{$refs->[$_]->[0]} = $_ ;
    }
    return bless( $self, $class ) ;
}

sub put {
    my( $self, $rec ) = @_ ;

    my %b2c = ( A => 1, C => 2, G => 4, T => 8, N => 15 ) ;

    my $ciglen = scalar @{$rec->{cigar} // []} ;
    my $seq ;
    for( my $i = 0 ; $i < length $rec->{seq} ; $i+=2 ) 
    {
        $seq .= chr( ( $b2c{ substr( $rec->{seq}, $i, 1 ) } << 4) |
            ( $i+1 == length $rec->{seq} ? 0 : $b2c{ substr( $rec->{seq}, $i+1, 1 ) } ) );
    }

    my $bin = 0 ;
    if( $rec->{pos} ) {
	my $beg = $rec->{pos} ;
	my $end = $rec->{pos} + $rec->{alnlength} -1 ;
	$bin = ($beg >> 14) == ($end >> 14) ? ((1<<15)-1)/7 + ($beg >> 14) :
	       ($beg >> 17) == ($end >> 17) ? ((1<<12)-1)/7 + ($beg >> 17) :
	       ($beg >> 20) == ($end >> 20) ? ((1<< 9)-1)/7 + ($beg >> 20) :
	       ($beg >> 23) == ($end >> 23) ? ((1<< 6)-1)/7 + ($beg >> 23) :
	       ($beg >> 26) == ($end >> 26) ? ((1<< 3)-1)/7 + ($beg >> 26) : 0 ;
    }

    my $buf = pack "VVVVVVVV a* V[$ciglen] a* a*",
    	exists $rec->{rname} ? $self->{refs}->{ $rec->{rname} } // 0xffffffff : 0xffffffff,
        $rec->{pos} // 0xffffffff,
        ($bin << 16) | (($rec->{mapq} // 0) << 8) | (1 + length ($rec->{qname} // '*')),
        ($rec->{flag} << 16) | $ciglen,
        length $rec->{seq},
        exists $rec->{mate_rID} ? $self->{refs}->{ $rec->{mate_rID} } // 0xffffffff : 0xffffffff,
        $rec->{mate_pos} // 0xffffffff,
        $rec->{ins_size} // 0,
        ($rec->{qname} // '*') . "\0",
        @{$rec->{cigar} // []},
        $seq,
        $rec->{qual} // '\xff' x length $rec->{seq} ;

        while( my( $tag, $val ) = each %{$rec->{opt_fields}} ) {
            if( ref($val) eq 'ARRAY' )
            {
                my( $pc, $bc ) = get_code( @$val ) ;
                my $n = scalar @$val ;
                $buf .= pack "a2 a a V $pc$n", $tag, 'B', $bc, $n, @$val ;
            }
            else
            {
                my( $pc, $bc ) = get_code( $val ) ;
                $buf .= pack "a2 a $pc", $tag, $bc, $val ;
            }
        }

    $self->{buf} .= pack "V", length $buf ;
    $self->{buf} .= $buf ;
    $self->flush() ;
}

sub flush {
    my( $self, $force ) = @_ ;

    while( $force || length $self->{buf} > 65000 ) 
    {
	my $chunk = substr( $self->{buf}, 0, 65000 ) ;
	$self->{buf} = length $self->{buf} > 65000 ? substr( $self->{buf}, 65000 ) : '' ;

	my $compdata ;
	gzip \$chunk => \$compdata ;
	my $l = length $compdata ;
	substr( $compdata, 3, 1 ) = chr( ord( substr( $compdata, 3, 1 )) | 4 ) ;
	substr( $compdata, 10, 0 ) = pack( "S< a2 S< S<", 6, "BC", 2, 7 + $l ) ;
	print {$self->{fh}} $compdata ;
	$force = 0 ;
    }
}
	
sub DESTROY {
    my $self = shift ;
    $self->flush( 1 ) ; # flush any remaining stuff
}

1 ;
__END__

=head1 NAME

Bam -- Reading and writing of BAM files

=head1 SYNOPSIS

    use Bam;

    my $input = Bam::In::new($input);
    while( my $rec = $input->fetch() ) { ... }

    my $output = Bam::Out::new ($fh, $refs);
    $output->put ($rec);

    $fastq_string = Bam::qual_to_fastq($rec->{qual});
    $rec->{qual} = Bam::qual_from_fastq($fastq_string);

=head1 DESCRIPTION

Bam provides simple streaming access to BAM (alignment) files, both
reading and writing.

=head1 OBJECT METHODS

=over 10

=item $input = Bam::In->new($input)

Opens a BAM file for reading.  The $input can be a file name or a file
handle.  Dies when anything goes wrong.  Else, $input-E<gt>{refs}
contains a list of reference sequences, where each sequence is
represented as a list containing the sequence name and its length.

=item $rec = $input->fetch()

Reads the next record from $input.  Returns undef if no more record are
available.  Records have named fields for the mandatory fields according
to the BAM specification and a single hash of all the optional fields
found.

=item $output = Bam::Out->new ($fh, $refs)

Opens a BAM file for writing.  The $fh must be file handle and $refs
must be a table of reference sequences, which could be obtained from the
refs entry of an input BAM file.

=item $output->put($rec)

Writes a record to an output file, in the same format as obtained from
$input->fetch().

=back

=head1 HELPER FUNCTIONS

=over 10

=item $fastq_string = Bam::qual_to_fastq($rec->{qual})

Convert the raw quality values in a BAM record into a string suitable
for dumping into a FASTQ file.

=item $rec->{qual} = Bam::qual_from_fastq($fastq_string)

Convert a quality string obtained from a FASTQ file into the raw
quality values needed in a BAM record.

=back

=head1 NOTES

The 'd' format specifier, apparently a samtools extension, is supported
when reading, but not when writing BAM files.  This may mean loss of
precision, but the output definitely follows the specification.

=back

=head1 SEE ALSO

http://samtools.sf.net

=head1 AUTHOR

Udo Stenzel E<lt>udo_stenzel@eva.mpg.deE<gt>.  

Copyright (c) 2011 MPI EVAN.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

