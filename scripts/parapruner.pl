#!/usr/bin/perl

# PIPES do not work as delimiter need to fix

# need to roll in /bwdata1/ariane/09-PARAPHYLY_PRUNING/process.pl

# we might want to check that $frac is less than 0.2 (or some other cutoff)
# print number of taxa retained in output

use strict;
use warnings;
use JFR::Fasta;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

our $VERSION = 0.12;

our $PY_GET_NODES = 'parapruner_get_nodes.py';
our $PY_IS_MONO   = 'parapruner_is_mono.py';
our $PY_PRUNE     = 'parapruner_prune.py';

MAIN: {
    my $rh_o = process_options();
#print "DELIM: $rh_o->{'delim'}\n";
    make_outdir($rh_o->{'outdir'});

#print STDERR "TREE: $rh_o->{'tree'}\n";
#print STDERR "ALN : $rh_o->{'aln'}\n";

    # get nodes from tree using $PY_GET_NODES
    my $ra_nodes = get_nodes($rh_o->{'tree'});

    # build has with key = defline value = sequence
    my $rh_aln = get_aln($rh_o->{'aln'});

    # check that all nodes in alignment are in tree and vice versa
    check_nodes_match_aln($ra_nodes,$rh_aln);

    # builds hash of arrays keys being species and values being identifiers
    my $rh_sp_nodes = get_sp_nodes($ra_nodes,$rh_o->{'tree'},$rh_o->{'delim'});

    # identify a species with only one identifier; die if none exist
    # the "singlet" will be used to root the tree 
    # a future version could do a 2-part rooting if no singlets existed
    my $root = get_single_sp($rh_sp_nodes);

    # use $PY_IS_MONO to test if multiple identifiers from a single species
    #    are monophyletic. return the non-monophyletic subset of $rh_sp_nodes 
    my $rh_to_prune = get_to_prune($root,$rh_o->{'tree'},$rh_sp_nodes); 

    # compare $rh_sp_nodes (complete) to $rh_to_prune (subset)
    #   and print the numbers that need to be pruned
    print_fraction_pruned($rh_o->{'tree'},$rh_sp_nodes,$rh_to_prune);

    # prune trees of species with non-monophyletic sequences
    prune($rh_o->{'tree'},$rh_to_prune,$rh_sp_nodes,$rh_o->{'outdir'});

    # prune alignments of species with non-monophyletic sequences
    prune_aln($rh_to_prune,$rh_aln,$rh_o->{'aln'},$rh_o->{'outdir'},$rh_o->{'delim'});
}

sub process_options {
    my $rh_o = {};
    my $opt_results = Getopt::Long::GetOptions(
                                "aln=s" => \$rh_o->{'aln'},
                               "tree=s" => \$rh_o->{'tree'},
                             "outdir=s" => \$rh_o->{'outdir'},
                              "delim=s" => \$rh_o->{'delim'},
                              "version" => \$rh_o->{'version'},
                                 "help" => \$rh_o->{'help'});
    die "$VERSION\n" if ($rh_o->{'version'});
    usage() if $rh_o->{'help'};
#    pod2usage({-exitval => 0, -verbose => 2}) if($rh_o->{'help'});
    print "  --aln is a required opt\n" unless ($rh_o->{'aln'});
    print "  --tree is a required opt\n" unless ($rh_o->{'tree'});
    print "  --delim is a required opt\n" unless ($rh_o->{'delim'});
    print "  --outdir is a required opt\n" unless ($rh_o->{'outdir'});
    if ($rh_o->{'tree'}) {
        die "treefile ($rh_o->{'tree'}) doesnt exist\n" unless (-e $rh_o->{'tree'});
    } 
    if ($rh_o->{'aln'}) {
        die "alignment ($rh_o->{'aln'}) doesnt exist\n" unless (-e $rh_o->{'aln'});
    }
    usage() unless ($rh_o->{'aln'} && $rh_o->{'tree'} && $rh_o->{'delim'} && $rh_o->{'outdir'});
    die "cannot use pipe as delimiter\n" if ($rh_o->{'delim'} eq '|' || $rh_o->{'delim'} eq '\|');
#    $rh_o->{'qm_delim'} = quotemeta($rh_o->{'delim'}); # escape special chars
    return $rh_o;
}

sub make_outdir {
    my $dir = shift;
#    die "$dir exists" if (-d $dir);
    mkdir $dir or die "cannot mkdir $dir:$!";
}

sub usage {
    print "parapruner.pl --aln FASTA_ALN --tree NEWICK_TREE --delim = DELIMITER --outdir OUTDIR [--help] [--version]\n";
    exit;
}

sub prune_aln {
    my $rh_tp  = shift;
    my $rh_aln = shift;
    my $aln    = shift;
    my $outdir = shift;
    my $delim = shift;
    my $filename = File::Basename::fileparse($aln);
    open OUT, ">$outdir/$filename.parapruned" or die "cannot open >$outdir/$filename.parapruned:$!";
    foreach my $id (keys %{$rh_aln}) {
        $id =~ m/^([^$delim]+)$delim/ or die "id $id does not have delimiter $delim";
        my $sp = $1; 
        next if ($rh_tp->{$sp});
        print OUT ">$id\n$rh_aln->{$id}\n";
    }
    close OUT;
}

sub check_nodes_match_aln {
    my $ra_n = shift;
    my $rh_a = shift;
    my %nodes = ();
    my $missing = '';
    foreach my $id (@{$ra_n}) {
        $nodes{$id}++;
        unless ($rh_a->{$id}) {
            $missing .= "$id\n";
        }
    }
    my $missing_from_tree = '';
    foreach my $id (keys %{$rh_a}) {
        unless ($nodes{$id}) {
            $missing_from_tree .= "$id\n";
        }
    }
    if ($missing || $missing_from_tree) {
        print "Missing from aln\n\n$missing" if ($missing);
        print "Missing from tree\n\n$missing_from_tree" if ($missing_from_tree);
        die;
    }
}

sub get_aln {
    my $fa  = shift;
    my %aln = ();
    my $fp  = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $aln{$id} = $rec->{'seq'};
    }
    return \%aln;
}

sub print_fraction_pruned {
    my $tree = shift;
    my $rh_sp_nodes = shift;
    my $rh_to_prune = shift;

    my $val1 = scalar(keys(%{$rh_sp_nodes}));
    my $val2 = scalar(keys(%{$rh_to_prune}));
    my $frac = $val2/$val1;
    print "$tree,$frac,$val2,$val1\n";
}

sub prune {
    my $tree = shift;
    my $rh_to_prune = shift;
    my $rh_sp_nodes = shift;
    my $outdir = shift;
    my @keepers = ();
    foreach my $sp (keys %{$rh_sp_nodes}) {
        next if ($rh_to_prune->{$sp});
        foreach my $leaf (@{$rh_sp_nodes->{$sp}}) {
            push @keepers, $leaf;
        }
    }
    my $name = File::Basename::fileparse($tree);
    my $csv = join ',', @keepers;
#    print "$csv\n";
    system "$PY_PRUNE -t $tree -p '$csv' -o $outdir/$name.parapruned > /dev/null";
}

sub get_to_prune {
    my $root = shift;
    my $tree = shift;
    my $rh_sp = shift;
    my %to_prune = ();
    foreach my $sp (keys %{$rh_sp}) {
        next if (scalar(@{$rh_sp->{$sp}}) == 1);
        my $csv = join ',', @{$rh_sp->{$sp}};
        my $out = `$PY_IS_MONO -t $tree -l '$csv' -r '$root'`;
        if ($out =~ m/^\(False/) {
            $to_prune{$sp}++;
        }
    }
    return \%to_prune;
}

sub get_single_sp {
    my $rh_sp = shift;
    foreach my $sp (keys %{$rh_sp}) {
        return $rh_sp->{$sp}->[0] if (scalar(@{$rh_sp->{$sp}}) == 1);
    }
    die "no species with only one sequence";
}

sub get_sp_nodes {
    my $ra_n  = shift;
    my $tree  = shift;
    my $delim = shift;
    my %sp   = ();
    foreach my $node (@{$ra_n}) {
        my $sp_ = '';
        if ($node =~ m/^([^$delim]+)$delim/) {
            $sp_ = $1;
        } else {
            die "unexpected: $node has no delimiter ($delim) in $tree";
        }
        push @{$sp{$sp_}}, $node;
    }
    return \%sp;
}

sub get_nodes {
    my $tree = shift;
    my @nodes = `$PY_GET_NODES -t $tree`;
    chomp @nodes;
    return \@nodes;
}

__END__

=head1 NAME

B<parapruner.pl> - prune paraphyletic taxa from tree and alignment

=head1 AUTHOR

Joseph Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

parapruner.pl --aln FASTA_ALN --tree NEWICK_TREE --delim DELIMITER --outdir OUTDIR [--help] [--version]

=head1 DESCRIPTION

--delim should divide species and unique id. for example pipe (--delim=\|) SPECIES|UNIQ_ID or period --delim=. 

=head1 BUGS

Please report them to the author.

=head1 COPYRIGHT

Copyright (C) 2018, Joseph Ryan

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
