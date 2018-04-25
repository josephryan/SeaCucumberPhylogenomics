#We used this script to split the blastp commands into .sh files of 5 commands or less. This allowed us to run multiple blastp searches at the same time using the appropriate amount of server space.

#!/usr/bin/env perl
use strict;
use warnings;

my $DIR = "/bwdata1/jessie/01-SEACUCUMBER/05-ORTHOFINDER";
my $file = "of.STD.out";
my $OUTDIR = "/bwdata1/jessie/01-SEACUCUMBER/05-ORTHOFINDER/files";


open(FILE, "<", "$DIR/$file") or die("Can't open file");
my @line_array = <FILE>;
my $round_count = 0;
my $line_count = 0;
my @new_file = ();
my $file_number_name = 1;
my $limit = 5;
my $five_count = 1;
my $count = 0;
while($count<@line_array){
        open(my $out_file, '>>', "$OUTDIR/file$file_number_name") or die;
        print $out_file $line_array[$count];
     if($five_count eq 5){
         $file_number_name ++;
         $five_count = 0;
     }


$five_count++;
$count++;
}
