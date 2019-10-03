#Exploring Bioinformatics: A Project-Based Approach
#Caroline St.Clair and Jonathan Visick

#Chapter 3 - Skills 1-3 Solution  

my($seq1, $len1, $seq2, $len2, $data, @matrix, $i, $j, $x, $y, $val1, $val2);
my($val3, $pathrow, $pathcol, $seq1loc, $seq2loc, $gapscore, $matchscore, $mismatchscore);
# modified by Gunes to display path from origin
my($forwardpath);

if(!open(infile1, 'fig16seq1.txt')){
    print "error opening input file 1\n";
    exit;
}
$data = <infile1>;   #ignore FASTA comment
while ($data = <infile1>){
   chomp $data;
   $seq1 = $seq1 . $data;
}

if(!open(infile2, 'fig16seq2.txt')){
    print "error opening input file 2\n";
    exit;
}
$data = <infile2>;   #ignore FASTA comment
while ($data = <infile2>){
   chomp $data;
   $seq2 = $seq2 . $data;
}

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

print "Enter the gap score: ";
$gapscore=<STDIN>;
chomp $gapscore;

print "Print scoring matrix and path to screen? (Y/N): ";
$print = <STDIN>;
chomp $print;

# declare a two dimensional array and initialize to spaces
# array must contain one extra row and one extra column
@matrix = ();
for($i = 0; $i < $len1; $i++){
   for($j = 0; $j < $len2; $j++){
      $matrix[$i][$j] = ' ';
   }
}

# initialize 1st row and 1st column of matrix
$matrix[0][0] = 0;
for ($i = 1; $i < $len1; $i ++){
    $matrix[$i][0] = $matrix[$i-1][0] + $gapscore;
}
for ($i = 1; $i < $len2; $i ++){
    $matrix[0][$i] = $matrix[0][$i-1] + $gapscore;
}


# STEP 1:
# Fill in rest of matrix using the following rules:
# determine three possible values for matrix[x][y]
# value 1 = add gap score to matrix[x][y-1]
# value 2 = add gap score to matrix[x-1][y]
# value 3 = add match score or mismatch score to 
#           matrix[x-1][y-1] depending on nucleotide 
#           match for position x of $seq1 and position y
#           of seq2
# place the largest of the three values in matrix[x][y]
#
# Best alignment score appears in matrix[$len1][$len2].

for($x = 1; $x < $len1; $x++){
   for($y = 1; $y < $len2; $y++){
	$val1 = $matrix[$x][$y-1] + $gapscore;
	$val2 = $matrix[$x-1][$y] + $gapscore;
	if (substr($seq1, $x, 1) eq substr($seq2, $y, 1)){
           $val3 = $matrix[$x-1][$y-1] + $matchscore;
	}
	else{
	   $val3 = $matrix[$x-1][$y-1] + $mismatchscore;
	}
	if (($val1 >= $val2) && ($val1 >= $val3)){
	   $matrix[$x][$y] = $val1;
	}
	elsif (($val2 >= $val1) && ($val2 >= $val3)){
	   $matrix[$x][$y] = $val2;
	}
	else{
	   $matrix[$x][$y] = $val3;
	}
   }
}


if (($print eq 'Y') || ($print eq 'y')){
   # Display scoring matrix
   print "MATRIX:\n";
   for($x = 0; $x < $len1; $x++){
      for($y = 0; $y < $len2; $y++){
	print "$matrix[$x][$y] ";
      }
      print "\n";
   }
   print "\n";
}
# STEP 2:
# Begin at matrix[$len1][$len2] and find a path to 
# matrix[0][0].
# Build string to hold path pattern by concatenating either 
# 'H' (current cell produced by cell horizontally to left), 
# 'D' (current cell produced by cell on diagonal), 
# 'V' (current cell produced by cell vertically above)
# It is possible for cell to have been produced by more
# than one of these possible choices, if multiple optimal 
# alignments exist. For now, choose only one.
# 
$pathrow = $len1-1;
$pathcol = $len2-1;
while (($pathrow != 0) || ($pathcol != 0)){
    if ($pathrow == 0){
       # must be from cell to left
       $path = $path . 'H';
       $pathcol = $pathcol - 1;
    }
    elsif ($pathcol == 0){
       # must be from cell above
       $path = $path . 'V';
       $pathrow = $pathrow - 1;
    }
    # could be from any direction
    elsif  (($matrix[$pathrow][$pathcol-1] + $gapscore) == $matrix[$pathrow][$pathcol]){
       # from left
       $path = $path . 'H';
       $pathcol = $pathcol - 1;
    }
    elsif (($matrix[$pathrow-1][$pathcol] + $gapscore) == $matrix[$pathrow][$pathcol]){
       # from above
       $path = $path . 'V';
       $pathrow = $pathrow - 1;
    } 
    else{
       # must be from diagonal
       $path = $path . 'D';
       $pathrow = $pathrow - 1;
       $pathcol = $pathcol - 1;
    }
}	

# following added to display forward path:
$forwardpath = reverse($path);

if (($print eq 'Y') || ($print eq 'y')){
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
   else{  # must be an H
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
