perl-ngs
========

This is just a collection of stuff that's generally useful, but wouldn't
fit well anywhere else, generally related to the processing of Next
Generation Sequencing data in Perl.  So far:

* **Twobit.pm**:   Self-contained Perl module to access 2bit files.  The
                   interface is loosely modelled after Bio::DB::Fasta, but it
                   doesn't use BioPerl.

* **Bam.pm**:      Perl module to read and write Bam files.  Not quite
                   complete, but useable.  If you have threaded Perl 2.8 or
                   later, some data structures are declared "shared".  To make
                   use of that, you must "use threads;" before loading the Bam
                   module.

* **Bed.pm**:      Very small Bed parser and printer, interfaces easily to
                   Bam.pm.  Probably only useful when translating between these
                   formats.

phylogeny\_from\_bam
--------------------

A tool that implements a restricted form of parsimony analysis to place
an ancient sample into an otherwise fixed phylogeny.  It uses
**sites.pm** and **track.pm**, or in a more compact but less readable
alternative, **uarray.pm**.  It also shows how to work with **Bam.pm**.
