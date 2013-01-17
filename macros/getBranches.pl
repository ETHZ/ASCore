#!/usr/bin/perl -s
#
# getBranches.pl
#
# Perl script to extract the list of ETH tree branches used in a list of files
# The output is the list of used branches. It should be redirected in a text file.
# This file can then be fed into makeTreeClassBase.py.
# 
# Usage: getBranches.pl <file1.cc> [file2.cc ...]
#

if ( @ARGV==0 || $h || $help ) {
    die "Usage: ".$0." <file1.cc> [file2.cc ...]\n";
}

my %branches ;
foreach $file ( @ARGV ) {

    open(FIN,$file) or die "Couldn't open $file: $!";
    while (<FIN>) {
        chomp();
        my $line = $_;
        next if ( $line =~ /^\s*\/\// );
        next if ( $line !~ /fTR->/);
        #print "-----\n$line\n";
        @list = split(/fTR/,$line);
        foreach $entry ( @list ) {
            next if ( $entry !~ /^->/);
            $entry =~ s/^->([\w]+).*/\1/;
        #    print "$entry\n";
            $branches{$entry} = 1;
        }
    }
    close(FIN);
}

foreach $branch (keys %branches) {
    print "$branch\n";
}

