#!/usr/bin/env perl

use strict;
use warnings;

package StringHelper;

=head1 StringHelper

Provides some string operations.

=head cutSuffix 

=cut

sub cutSuffix {
    my $string = shift;
    $string =~ s/\.([^\.]+$)//;
    return $string;
}

=head2 replaceOverhangChar

Usage:
    
    my $string = replaceOverhangChar($string, $char, $replace_char);

Function:

    Replaces $char on left and right side of $string with $replace_char.

=cut

sub replaceOverhangChar {
    my $string       = shift;
    my $char         = shift;
    my $replace_char = shift;
    if ( $string =~ /^($char+)/ ) {
        substr $string, 0, length($1), '';
        $string = $replace_char x length($1) . $string;
    }
    if ( $string =~ /($char+)$/ ) {
        substr $string, length($string) - length($1), length($1), '';
        $string = $string . $replace_char x length($1);
    }
    return $string;
}

=head2 getOverhang

Usage: 

    my ($left, $right) = getOverhang($string, $char))

Example:

    my ($left, $right) = getOverhang("XXXXABCXX", "X")

Function: 
    
    Returns number of characters on left and right side of $string. 

Returns:
    
    List with both integers.

=cut

sub getOverhang {
    my $string = shift;
    my $char   = shift;

    my ( $left, $right ) = ( 0, 0 );
    if ( $string =~ /^($char+)/ ) {
        $left = length($1);
    }

    if ( $string =~ /($char+)$/ ) {
        $right = length($1);
    }

    return ( $left, $right );
}

=head2 getRange
                  123456789012
    my $string = "NNNNAAAAANNN";
    my ($start, end) = StringHelper::getRange("$string", "N");
    start = 5, end = 9

    my $string = "NNNNNNNNN";
    start = 1, end = 1
=cut

sub getRange {
    my $string = shift;
    my $char   = shift;

    my ( $left, $right ) = ( 0, 0 );
    if ( $string =~ /^($char+)/ ) {
        $left = length($1);
    }

    if ( $string =~ /($char+)$/ ) {
        $right = length($1);
    }

    my $end = length($string) - $right;
    $end = 1 if ($end == 0); 
        

    return ( $left + 1, $end );
}


=head2 stringSplitter

Usage:

    my @result = @{stringSplitter($string, $width)};

Function:

    Returns splitted $string by $width as array reference.

=cut

sub stringSplitter {
    my $string = shift;
    my $width  = shift;

    my @result;
    my $prev_i = 0;
    my $i      = $width;
    my $l      = 0;
    for ( ; $i <= length($string) ; $i += $width ) {
        push @result, substr $string, $prev_i, $width;
        $prev_i = $i;
    }
    if ( $prev_i < length($string) ) {
        push @result, substr $string, $prev_i;
    }

    return \@result;
}

1;

__END__
