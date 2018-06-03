#!/usr/bin/perl

# Rong Chen 6/24/2002
# modified Brian Pierce 11/22/2010
# modified by Fayyaz 5/20/2013 (Added support for reading and writing files to directories other than the current one)
# This program creates the structural predictions from a zdock output file

#Usage: create.pl zdock_output_file num_preds (opt) outfilename (opt) pdb_dir (opt) no_ms (opt)
#	zdock_output_file: path of output file from ZDOCK
#	num_preds: Number of predicted structures to produce (0: for all, positive: for number of structures, when negative it will produce only that structure)

use strict;

my $outfile = $ARGV[0];

if ($outfile eq "")
{
    die("usage: create.pl zdock_output_file num_preds (opt) outfilename (opt) pdb_dir (opt) no_ms (opt)\n");
}
my $num_preds = $ARGV[1];
if (($num_preds eq "") || ($num_preds == 0)) { $num_preds = 100000; }
## Added by Fayyaz
my $sgl=0;
if ($num_preds < 0) {$num_preds=-$num_preds; $sgl=1; }
my $ofilename = $ARGV[2];
if ($ofilename eq "") {$ofilename = "complex"}
my $bdir = $ARGV[3];
#if ($bdir eq "") {$bdir = ""}
if (!(-e "./create_lig")) { die("error: need to have create_lig executable linked or copied to current directory\n"); }
my $noms = $ARGV[4];
if (($noms eq "") || ($noms <= 0)) { $noms = 0; }
# open the zdock output file
open (ZDOUT, "$outfile") || die "Unable to open file $outfile!\n";
my @zdout_lines = <ZDOUT>;
chomp(@zdout_lines);
close(ZDOUT);

# parse the header of the zdock output file
(my $n, my $spacing, my $switch_num)=split(" ", $zdout_lines[0]);
my $line_num = 1;
my $rec_rand1 = 0.0, my $rec_rand2 = 0.0, my $rec_rand3 = 0.0;
if ($switch_num ne "")
{
    ($rec_rand1, $rec_rand2, $rec_rand3)= split(" ", $zdout_lines[$line_num++]);
}
(my $lig_rand1, my $lig_rand2, my $lig_rand3)=split(" ", $zdout_lines[$line_num++]);
(my $rec, my $r1, my $r2, my $r3) = split (" ", $zdout_lines[$line_num++]);
(my $lig, my $l1, my $l2, my $l3) = split (" ", $zdout_lines[$line_num++]);

if ($switch_num eq "1")
{
    my $temp_name = $rec;
    $rec = $lig;
    $lig = $temp_name;
}
##Added by Fayyaz
$rec = $bdir . $rec;
$lig = $bdir . $lig;
if ($noms) #Remove .ms
{
$rec = substr($rec,0,-3);
$lig = substr($lig,0,-3);
}
# generate the predictions
my $pred_num = 1;
for (my $i = $line_num; ($i < @zdout_lines) && ($pred_num <= $num_preds); $i++)
{
	if (($sgl==0) || ($pred_num==$num_preds))
	{
		(my $angl_x, my $angl_y, my $angl_z, my $tran_x, my $tran_y, my $tran_z, my $score) = split ( " ", $zdout_lines[$i] );
		my $newligfile = $outfile . "." . $pred_num;
		my $create_cmd = "./create_lig $newligfile $lig $lig_rand1 $lig_rand2 $lig_rand3 $r1 $r2 $r3 $l1 $l2 $l3 $angl_x $angl_y $angl_z $tran_x $tran_y $tran_z $n $spacing\n";
		if ($switch_num ne "")
		{
		$create_cmd = "./create_lig $newligfile $lig $switch_num $rec_rand1 $rec_rand2 $rec_rand3 $lig_rand1 $lig_rand2 $lig_rand3 $r1 $r2 $r3 $l1 $l2 $l3 $angl_x $angl_y $angl_z $tran_x $tran_y $tran_z $n $spacing\n";
		}

		#print "executing: $create_cmd\n";
		system($create_cmd);
		system "cat $rec $newligfile > $ofilename." . "$pred_num.pdb";
		system "rm $newligfile\n";
	}
    $pred_num++;
}