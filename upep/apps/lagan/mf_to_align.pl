#!/usr/bin/env perl

# Copyright (c) 2011, The University of Queensland
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Queensland nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY OF QUEENSLAND BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(min max);

BEGIN {
    select STDERR; $| = 1;
    select STDOUT; $| = 1;
}

my $options = check_params();

# Check if BioPerl is installed
eval "use Bio::SeqIO";

# Use regex to split FASTA file if no BioPerl, use Bio::SeqIO otherwise
my $seq1;
my $seq2;
if ($@) {
    open(my $fh, $options->{'f'});
    my $fasta_string = '';
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            $fasta_string .= ">";
        } else {
            $fasta_string .= $line;
        }
    } 
    close($fh);
    my $seq1;
    my $seq2; 
    if ($fasta_string =~ /^>([A-Za-z\-]*)>([A-Za-z\-]*)/) {
        $seq1 = $1;
        $seq2 = $2;
    }
} else {
    my $fasta_stream = Bio::SeqIO->new( "-format" => "fasta",
                                    "-file" => $options->{'f'});
    $seq1 = $fasta_stream->next_seq()->seq();
    $seq2 = $fasta_stream->next_seq()->seq();
}

if ($options->{'out'}) {
    open(my $out_fh, ">", $options->{'out'});
    print {$out_fh} $seq1."\n";
    print {$out_fh} get_match_string($seq1,$seq2)."\n";
    print {$out_fh} $seq2."\n";
    close($out_fh);
} else {
    print $seq1."\n";
    print get_match_string($seq1,$seq2)."\n";
    print $seq2."\n";
}

sub get_match_string {
    my $seq1 = shift;
    my $seq2 = shift;
    my @seq1array = split(//, $seq1);
    my @seq2array = split(//, $seq2);
    my $matchstring;
    for (my $i = 0; $i <= $#seq1array; $i++) {
        if ($seq1array[$i] eq $seq2array[$i]) {
            $matchstring .= '|'
        } else {
            $matchstring .= ' ';
        }
    }
    return $matchstring;
}

sub check_params {

    my @standard_options = ( "help+", "man+");
    my %options;

    GetOptions( \%options, @standard_options, "f:s", "out:s" );
    exec("pod2usage $0") if $options{'help'};
    exec("perldoc $0")   if $options{'man'};
    exec("pod2usage $0")
        if ( !( $options{'f'} && $options{'out'}));

    return \%options;
}

__DATA__

=head1 NAME

    mf_to_align.pl
   
=head1 DESCRIPTION

    Converts a multifasta alignment (between two sequences) to the
    alignment format required by uPEPeroni.

=head1 SYNOPSIS

    mf_to_align.pl -f <multifasta.fa> [-out <outfile>]

        -f    Multifasta file containing alignment of two sequences.
        -out  Output file (default: stdout) 

=cut



