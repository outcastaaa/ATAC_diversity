#!/usr/bin/perl

use strict;
use warnings;


my %file1_hash;
open my $file1, '<', 'list1.txt' or die "Cannot open list1.txt: $!";
while (<$file1>) {
    chomp;
    my ($col1, $col2) = split /\t/;
    $file1_hash{$col2} = $col1;
}
close $file1;

open my $file2, '<', 'list2.txt' or die "Cannot open list2.txt: $!";
while (<$file2>) {
    chomp;
    my ($col1, $col2) = split /\t/;
    if (exists $file1_hash{$col2} && $file1_hash{$col2} ne $col1) {
        print "$col2\n";
    }
}
close $file2;
