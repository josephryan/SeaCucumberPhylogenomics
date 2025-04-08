#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;
use File::Copy;

#created: Wed Sep 12 09:54:41 EDT 2018
our $VERSION = 0.01;

our $DELIM = '\|';
our $PARA_DIR = '/bwdata1/jessie/01-SEACUCUMBER/06-PARAPRUNE/02-ALN';
our $MAX_PARA = 5;
our $MIN_SEQ_LEN  = 50;

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
        my $para_count = scalar(keys(%sp));
        next unless ($para_count >= $min_num_sp);
        $fa =~ m/(OG\d+)/ or die "unexpected";
        my $fp2 = JFR::Fasta->new("$PARA_DIR/$1.mafft-gb.res");
        my %sp2 = ();
        my $longest = 0; # longest non gap
        while (my $rec = $fp2->get_record()) {
            $rec->{'seq'} =~ s/-//;
            $longest = length($rec->{'seq'}) unless ($longest > length($rec->{'seq'}));
            my @f = split /\./, $rec->{'def'};
            $sp2{$f[0]}++;
        }
        next unless ($longest >= $MIN_SEQ_LEN);
        my $orig_count = scalar(keys(%sp2));
        print "$fa\n" if (($orig_count - $para_count) <= $MAX_PARA);
    }
}

    


