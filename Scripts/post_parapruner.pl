#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;

our $DELIM = '\|';

MAIN: {
    my $min_num_sp = $ARGV[0] or die "usage: $0 MIN_NUM_SP FA [FA FA...]\n";    
    my @fas = @ARGV;
    shift @fas; # $min_num_sp
    foreach my $fa (@fas) {
        my $fp = JFR::Fasta->new($fa);
        my %sp = ();
        while (my $rec = $fp->get_record()) {
            my @f = split /$DELIM/, $rec->{'def'};
            $sp{$f[0]}++;
        }
        print "$fa\n" if (scalar(keys(%sp)) >= $min_num_sp);
    }
}

