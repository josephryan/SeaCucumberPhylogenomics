#!/usr/bin/perl

use lib qw(/home/s9/jfryan/lib);
use strict;
use warnings;
use Getopt::Long;
use JFR::Fasta;
use Pod::Usage;

our $VERSION = 0.06;
# version 0.04 allows for Trinity2 deflines
# version 0.05 allows for Trinity2.4.0 deflines
# version 0.06 corrects a bug in 0.05 which supposedly allowed for Trinity2.4.0
our $AUTHOR  = 'Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>';

MAIN: {
    my $opt_version = 0;
    my $opt_help = 0;
    my $opt_trinity = 0;
    my $opt_trinity2 = 0;
    my $opt_pad = '';
    my $opt_fasta = '';
    my $opt_prefix = '';

    my $opt_results = Getopt::Long::GetOptions(  "version" => \$opt_version,
                                   "trinity" => \$opt_trinity,
                                   "trinity2" => \$opt_trinity2,
                                     "pad=i" => \$opt_pad,
                                   "fasta=s" => \$opt_fasta,
                                  "prefix=s" => \$opt_prefix,
                                      "help" => \$opt_help);

    $opt_version && version();
    pod2usage({-exitval => 0, -verbose => 2}) if $opt_help;
    $opt_fasta || usage();
    $opt_prefix || usage();
    my $count = 0;
    my $fp = JFR::Fasta->new($opt_fasta);
    while (my $rec = $fp->get_record()) {
        if ($opt_trinity) {
            if ($rec->{'def'} =~ m/^>comp\d+_c\d+_seq(\d+)/) {
                my $iso = $1;
                $count++ if ($iso == 1);
                $count = sprintf("%0${opt_pad}d", $count) if ($opt_pad);
                print ">$opt_prefix.$count.$iso\n$rec->{'seq'}\n";
            } elsif ($rec->{'def'} =~ m/^>TRINITY\S+_i(\d+) /) {
                my $iso = $1;
                $count++ if ($iso == 1);
                $count = sprintf("%0${opt_pad}d", $count) if ($opt_pad);
                print ">$opt_prefix.$count.$iso\n$rec->{'seq'}\n";
            } else {
                die "expected trinity formatted defline: $rec->{'def'}";
            }
        } elsif ($opt_trinity2) {
            if ($rec->{'def'} =~ m/^>TR(\d+)\|(\S+)/) {
                my $iso   = $2;
                $count++ if ($iso == 1);
                $count = sprintf("%0${opt_pad}d", $count) if ($opt_pad);
                print ">$opt_prefix.$count.$iso\n$rec->{'seq'}\n";
            } elsif ($rec->{'def'} =~ m/^>TRINITY_[^_]+_c\d+_g\d+_i(\d+)_/) {
                my $iso   = $1;
                $count++ if ($iso == 1);
                $count = sprintf("%0${opt_pad}d", $count) if ($opt_pad);
                print ">$opt_prefix.$count.$iso\n$rec->{'seq'}\n";
            } else {
                die "expected trinity ver. 2 formatted defline: $rec->{'def'}";
            }
        } else {
            $count++;
            $count = sprintf("%0${opt_pad}d", $count) if ($opt_pad);
            print ">$opt_prefix.$count\n$rec->{'seq'}\n";
        }
    }
}

sub usage {
    die "$0 [--version] [--help] [--trinity] [--trinity2] [--pad=INT] --fasta=FASTA --prefix=DEFLINE_PREFIX\n";
}

sub version {
    die "replace_deflines.pl $VERSION\n";
}

__END__

=head1 NAME

B<replace_deflines.pl>

=head1 AUTHOR

Joseph F. Ryan <josephryan@yahoo.com>

=head1 SYNOPSIS

 replace_deflines.pl [--version] [--help] [--trinity] [--trinity2] [--pad=INT] --fasta=FASTA --prefix=DEFLINE_PREFIX\n";

=head1 OPTIONS

=item B<--trinity>

add this flag if you have trinity deflines and you would like to preserve isoform information

=item B<--pad=INT>

instead of >blah.1 --pad=3 would give >blah.001

=item B<--fasta=FASTA>

FASTA to have it's deflines replaced

=item B<--version>

print version and exit

=item B<--help>

Print this manual


=head1 BUGS

Please report them to <josephryan@yahoo.com>

=head1 COPYRIGHT

Copyright (C) 2012,2013 Joseph F. Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
