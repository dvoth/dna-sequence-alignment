use strict;
use warnings;

# The original global alignment program is modified significantly by Gunes to make an Affine alignment

# matrices @Hmatrix and @Vmatrix added, and original @matrix renamed @Dmatrix
# $gapextend and $gapstart scores replace the single $gapscore
my($seq1, $len1, $seq2, $len2, $data, $i, $j, $x, $y, $val1, $val2);
my(@Dmatrix);
my(@Hmatrix);
my(@Vmatrix);
my($val3, $pathrow, $pathcol, $seq1loc, $seq2loc, $gapstart, $gapextend, $matchscore, $mismatchscore);
my($neginf); # added to simulate the "negative infinity" scores of some initializations
my($hval1, $hval2, $vval1, $vval2);
my($lastD, $lastH, $lastV);
my($layer); # this will help with finding the path back to origin
my ($s); # substitution score
my($forwardpath); # path in forwards direction

# more declarations
my($filename1,$filename2);
my($path);
my($print1);
my($pathloc);
my($alignseq1);
my($alignseq2);
my($count);


# Modified by Gunes to prompt user for filename:
# print "\nEnter filename for first sequence: ";
$filename1 = "affineseq1.txt";
# print "\nEnter filename for second sequence: ";
$filename2 = "affineseq2.txt";
print "\n";

# Also modified to initialize seq1 and seq2:
$seq1 = "";
$seq2 = "";

if (!open(infile1,$filename1)) {
    print "error opening input file 1\n";
    exit;
}
$data = <infile1>;   #ignore FASTA comment
while ($data = <infile1>){
   #chomp $data; MODIFIED OUT
   $seq1 = $seq1 . $data;
}

if (!open(infile2,$filename2)){
    print "error opening input file 2\n";
    exit;
}
$data = <infile2>;   #ignore FASTA comment
while ($data = <infile2>){
   #chomp $data; MODIFIED OUT
   $seq2 = $seq2 . $data;
}

# Following added by Gunes to replace the individual "chomps"
# It is the alternative way of replacing every newline,carriage return, and space with nothing.
$seq1 =~ s/\n//g; $seq1 =~ s/\r//g; $seq1 =~ s/" "//g;
$seq2 =~ s/\n//g; $seq2 =~ s/\r//g; $seq2 =~ s/" "//g;

# Following also added by Gunes to convert everything to UPPERCASE
$seq1 =~ tr/a-z/A-Z/;
$seq2 =~ tr/a-z/A-Z/;


# adding extra characters so sequences align with matrix
# saves some calculations later on
$seq1 = " " . $seq1;
$seq2 = " " . $seq2;
$len1 = length($seq1);
$len2 = length($seq2);

print "Enter the match score: ";
$matchscore=<STDIN>;
chomp $matchscore;

print "Enter the mismatch score: ";
$mismatchscore=<STDIN>;
chomp $mismatchscore;

print "Enter the initial gap score: ";
$gapstart=<STDIN>;
chomp $gapstart;

print "Enter the gap extension score: ";
$gapextend=<STDIN>;
chomp $gapextend;

print "Print scoring matrix and path to screen? (Y/N): ";
$print1 = <STDIN>;
chomp $print1;

# declare all two dimensional arrays and initialize to spaces
# arrays must contain one extra row and one extra column
@Dmatrix = (); @Hmatrix = (); @Vmatrix = ();
# We need to simulate "negative infinity" using the worst possible alignment score of gapstart*(len1+len2)
$neginf = ($len1+$len2)*$gapstart;
for($i = 0; $i < $len1; $i++){
   for($j = 0; $j < $len2; $j++){
      $Dmatrix[$i][$j] = int($neginf);
      $Hmatrix[$i][$j] = int($neginf);
      $Vmatrix[$i][$j] = int($neginf);
#      print "\nDEBUG: Dmatrix[$i][$j] is $DMatrix[$i][$j]\n";
#      print "\nDEBUG: Hmatrix[$i][$j] is $HMatrix[$i][$j]\n";
#      print "\nDEBUG: Vmatrix[$i][$j] is $VMatrix[$i][$j]\n";
   }
}

# now correctly initialize Origin of OPT, first Column of Vmatrix and first Row of Hmatrix
$Dmatrix[0][0] = 0;
$Hmatrix[0][1] = $gapstart;
$Vmatrix[1][0] = $gapstart;
#print "\nDEBUG: Vmatrix[1][0] is $VMatrix[1][0]\n";
#print "\nDEBUG: Hmatrix[0][1] is $HMatrix[0][1]\n";

# consider 0th column only
for ($i = 2; $i < $len1; $i ++){
    #$Vmatrix[$i][0] = $Vmatrix[$i-1][0] + $gapextend;
    $Vmatrix[$i][0] = $gapstart + ($i-1)*$gapextend;
#    print "\nDEBUG: Vmatrix[$i][0] is $VMatrix[$i][0]\n";
}


# consider 0th row only
for ($i = 2; $i < $len2; $i ++){
    #$Hmatrix[0][$i] = $Hmatrix[0][$i-1] + $gapextend;
    $Hmatrix[0][$i] = $gapstart + ($i-1)*$gapextend;
#    print "\nDEBUG: Hmatrix[0][$i] is $HMatrix[0][$i]\n";
}


# Fill in rest of matrices using the following rules:
# 
# determine three possible values for Dmatrix[x][y]
# value 1 = Hmatrix[x-1][y-1]+sub-or-match score
# value 2 = Vmatrix[x-1][y-1]+sub-or-match score
# value 3 = Dmatrix[x-1][y-1]+sub-or-match score
#
# determine two possible values for Hmatrix[x][y]
# hvalue 1 = Hmatrix[x][y-1] + $gapextend
# hvalue 2 = Dmatrix[x][y-1] + $gapstart
#
# determine two possible values for Vmatrix[x][y]
# vvalue 1 = Vmatrix[x-1][y] + $gapextend
# vvalue 2 = Dmatrix[x-1][y] + $gapstart
#
# place the largest of the respective values into the respective matrix
#
# Best alignment score appears in matrix[$len1][$len2].

for($x = 1; $x < $len1; $x++){
   for($y = 1; $y < $len2; $y++){
       # first calculate the possible substitution score
       if (substr($seq1, $x, 1) eq substr($seq2, $y, 1)){
	   $s = $matchscore;
       }
       else { $s = $mismatchscore; }
	   
       # calculate possible horizontal value
       if ($Hmatrix[$x][$y-1] == $neginf) { $hval1 = $neginf; }
       else { $hval1 = $Hmatrix[$x][$y-1] + $gapextend; }
       
       if ($Dmatrix[$x][$y-1] == $neginf) { $hval2 = $neginf; }
       else { $hval2 = $Dmatrix[$x][$y-1] + $gapstart; }

       # calculate possible vertical value
       if ($Vmatrix[$x-1][$y] == $neginf) { $vval1 = $neginf; }
       else { $vval1 = $Vmatrix[$x-1][$y] + $gapextend; }
       
       if ($Dmatrix[$x-1][$y] == $neginf) { $vval2 = $neginf; }
       else { $vval2 = $Dmatrix[$x-1][$y] + $gapstart; }
       
       # calculate possible middle value
       if ($Hmatrix[$x-1][$y-1] == $neginf) { $val1 = $neginf; }
       else { $val1 = $Hmatrix[$x-1][$y-1] + $s; }

       if ($Vmatrix[$x-1][$y-1] == $neginf) { $val2 = $neginf; }
       else { $val2 = $Vmatrix[$x-1][$y-1] + $s; }
       
       if ($Dmatrix[$x-1][$y-1] == $neginf) {
	   $val3 = $neginf;
       }
       else {
	   $val3 = $Dmatrix[$x-1][$y-1] + $s;
       }
       
       # now to determine max for OPT matrix
       if (($val1 >= $val2) && ($val1 >= $val3)){
	   $Dmatrix[$x][$y] = $val1;
       }
       elsif (($val2 >= $val1) && ($val2 >= $val3)){
	   $Dmatrix[$x][$y] = $val2;
       }
       else{
	   $Dmatrix[$x][$y] = $val3;
       }

       # now to determine max for Hmatrix
       if ($hval1 >= $hval2) {	$Hmatrix[$x][$y] = $hval1; }
       else { $Hmatrix[$x][$y] = $hval2; }
       
       # and finally determine max for Vmatrix
       if ($vval1 >= $vval2) { $Vmatrix[$x][$y] = $vval1; }
       else { $Vmatrix[$x][$y] = $vval2; }
       
   }
}


if (($print1 eq 'Y') || ($print1 eq 'y')){
    # Display scoring matrices
    print "Middle MATRIX:\n";
    for($x = 0; $x < $len1; $x++){
	for($y = 0; $y < $len2; $y++){
	    $lastD = $Dmatrix[$x][$y];
	    print "$Dmatrix[$x][$y] ";
	}
	print "\n";
    }
    print "\n";
    print "Horizontal MATRIX:\n";
    for($x = 0; $x < $len1; $x++){
	for($y = 0; $y < $len2; $y++){
	    $lastH = $Hmatrix[$x][$y];
	    print "$Hmatrix[$x][$y] ";
	}
	print "\n";
    }
    print "\n";
    print "Vertical MATRIX:\n";
    for($x = 0; $x < $len1; $x++){
	for($y = 0; $y < $len2; $y++){
	    $lastV = $Vmatrix[$x][$y];
	    print "$Vmatrix[$x][$y] ";
	}
	print "\n";
    }
    print "\n";
}

# Now to find the alignment path
# Begin at the MAX valued last matrix entry over the three matrices!
# This is a modification: That value is NOT necessarily at the DMatrix!
# Find a path from max_{M=D,H,V} (Mmatrix[n][m]) entry BACK TO the origin [0,0] 
# Build string to hold path pattern by concatenating either 
# 'H' (current cell produced by cell horizontally to left), 
# 'D' (current cell produced by cell on diagonal), 
# 'V' (current cell produced by cell vertically above)
# It is possible for cell to have been produced by more
# than one of these possible choices, if multiple optimal 
# alignments exist. For now, choose only one.
# 
$pathrow = int($len1-1);
$pathcol = int($len2-1);

#$layer = 'D'; # this helps keep track of matrices, and is initialized to the middle layer BUT should be initialized to whatever layer has maximum last entry!
# MODIFICATION:
if (($lastV >= $lastH) && ($lastV >= $lastD)) {$layer = 'V';}
if (($lastD >= $lastH) && ($lastD >= $lastV)) {$layer = 'D';}
if (($lastH >= $lastV) && ($lastH >= $lastD)) {$layer = 'H';}
print "Initialized backtrack layer: $layer\n"; 
###


while (($pathrow != 0) || ($pathcol != 0)){
    # first calculate $s the substitution score again
    if (substr($seq1, $pathrow, 1) eq substr($seq2, $pathcol, 1)){
	$s = $matchscore;
    }
    else { $s = $mismatchscore; }
    
    if ($pathrow == 0){
       # When you're on the 0th row:
       # must be via a horizontal gap
       $path = $path . 'H';
       $pathcol = $pathcol - 1;
       $layer = 'H';
    }
    elsif ($pathcol == 0){
       # When you're on the 0th column:
       # must be via a vertical gap
       $path = $path . 'V';
       $pathrow = $pathrow - 1;
       $layer = 'V';
    }

    # Otherwise more cases to consider:
    elsif ($layer eq 'D') { # if we were in the middle layer
	# must be because of a substitution move
	# now need to possibly update the layer if previous move was from H or V
	if ($Dmatrix[$pathrow][$pathcol] == $Hmatrix[$pathrow-1][$pathcol-1] + $s) {
	    $layer = 'H'; # need to move to the horizontal matrix
	}
	elsif ($Dmatrix[$pathrow][$pathcol] == $Vmatrix[$pathrow-1][$pathcol-1] + $s) {
	    $layer = 'V'; # need to move to the vertical matrix
	} # else don't change the layer anyway
	# do the following updates in any case:
	$path = $path . 'D';
	$pathrow = $pathrow - 1;
	$pathcol = $pathcol - 1;
    }

    elsif ($layer eq 'H') { # if we were in the horizontal layer
	# but need to switch to the middle layer if this was a gap start move
	if ($Hmatrix[$pathrow][$pathcol] == $Dmatrix[$pathrow][$pathcol-1] + $gapstart) {
	    $layer = 'D';
	}
	# do the following updates in any case
	$path = $path . 'H'; # current move is a gap
	$pathcol = $pathcol - 1;
    }

    elsif ($layer eq 'V') { # only remaining situation is that we're in the vertical layer
	# but need to switch to the middle layer if this was a gap start move
	if ($Vmatrix[$pathrow][$pathcol] == $Dmatrix[$pathrow-1][$pathcol] + $gapstart) {
	    $layer = 'D';
	}
	# do the following updates in any case
	$path = $path . 'V'; # current move is a gap
	$pathrow = $pathrow - 1;
    }  
} # this ends the while	



# following added to display forward path:
$forwardpath = reverse($path);
if (($print1 eq 'Y') || ($print1 eq 'y')){
    print "Backtrack path is $path\n";
    print "Forwards path is $forwardpath\n";
}

# STEP 3:
# Determine alignment pattern by reading path string 
# created in step 2.
# Create two string variables ($alignseq1 and $alignseq2) to hold 
# the sequences for alignment.
# Reading backwards from $seq1, $seq2 and path string, 
#   if string character is 'D', then 
#      concatenate to front of $alignseq1, last char in $seq1 
#      and to the front of $alignseq2, last char in $seq2
#   if string character is 'V', then 
#      concatenate to front of $alignseq1, last char in $seq1 
#      and to the front of $alignseq2, the gap char
#   if string character is 'H', then 
#      concatenate to front of $alignseq1 the gap char
#      and to the front of $alignseq2, last char in $seq2
# Continue process until path string has been traversed.
# Result appears in string $alignseq1 and $seq2
#

$seq1loc = $len1-1;
$seq2loc = $len2-1;
$pathloc = 0;

# initializing aligned strings
$alignseq1 = "";
$alignseq2 = "";

while ($pathloc < length($path)){
   if (substr($path, $pathloc, 1) eq 'D'){
      $alignseq1 = substr($seq1, $seq1loc, 1) . $alignseq1;
      $alignseq2 = substr($seq2, $seq2loc, 1) . $alignseq2;
      $seq1loc--;
      $seq2loc--;
   }
   elsif (substr($path, $pathloc, 1) eq 'V'){
      $alignseq1 = substr($seq1, $seq1loc, 1) . $alignseq1;
      $alignseq2 = '-' . $alignseq2;
      $seq1loc--;
   }  
   elsif  (substr($path, $pathloc, 1) eq 'H') {  # must be an H
      $alignseq1 = '-' . $alignseq1;
      $alignseq2 = substr($seq2, $seq2loc, 1) . $alignseq2;
      $seq2loc--;
   }  
   $pathloc++;
}

print "\nAligned Sequences:\n";
for ($i=0; $i+50<length($alignseq1); $i=$i+50){
	print substr($alignseq1, $i, 50) . "\n";
	print substr($alignseq2, $i, 50) . "\n\n";
}
# print whatever is left
print substr($alignseq1, $i) . "\n";
print substr($alignseq2, $i) . "\n";

$count = 0;
for ($i=0; $i < length($alignseq1); $i++){
   if (substr($alignseq1, $i, 1) eq substr($alignseq2, $i, 1)){
      $count++;
   }
}
print "Total number of aligned characters: $count\n\n";
print "Percentage match: ", $count / length($alignseq1), "\n\n";
# statement may be needed to hold output screen
print "Press any key to exit program\n";
$x = <STDIN>;	
