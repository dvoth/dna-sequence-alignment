/*
 * File: Main.java
 * Name: Dalton Voth
 * Date: 9-19-2019
 * Course: CS-490 Bioinformatics
 * Desc: Takes 2 input files in DNA FASTA and performs sequence alignment
 */

package SequenceAlignment;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
/**
 *
 * @author dvoth
 */
public class Main {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        char choice = 'b';


        while(Character.toUpperCase(choice) != 'A' && Character.toUpperCase(choice) != 'G' && Character.toUpperCase(choice) != 'L') {
          System.out.print("Choose (G)lobal, (L)ocal, or (A)ffinity alignment: ");
          choice = scanner.next().charAt(0);
        }

        if (choice == 'L') {
          localAlignment();
        } else if (choice == 'A') {
          affinityAlignment();
        } else {
          globalAlignment();
        }
    }
    
    public static String getFileNameFromUserInput() {
        Scanner scanner = new Scanner(System.in); 
        String inputFileName = "";
        boolean fileExists = false;

        while (!fileExists) {
            inputFileName = scanner.nextLine().trim();
            File input = new File("Ch3/" + inputFileName); 
            if (!input.exists()) {
                System.out.print("File " + input.getName() + " not found, try again: ");
            } else {
                fileExists = true;
            }
        }
        
        
        return inputFileName;
    }
    
    public static String getSequenceFromFile(File file) {
        String sequence = "", fileHeader = "";
        
        try {
            Scanner scanner = new Scanner(file);
            // Get the file header
            fileHeader = scanner.nextLine();
            
            // Read the rest of the file as a single string, remove whitespace, to uppercase
            sequence = scanner.useDelimiter("\\Z").next().replaceAll("\\s+","").toUpperCase();
            
            // Output header and sequence separately
            // System.out.println("Header: " + fileHeader + "\n");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(0);
        }
        
        return sequence;
    }

    public static void localAlignment() {
        String seq1FileName, seq2FileName;
        File seq1File, seq2File;
        String seq1, seq2, path = "", alignseq1 = "", alignseq2 = "", forwardpath = "";
        Scanner scanner = new Scanner(System.in); 
        int gapscore, matchscore, mismatchscore, len1, len2, val1, val2, val3, pathrow, pathcol, seq1loc, seq2loc, pathloc;
        char print;
        int[][] matrix;
        
        System.out.print("Sequence 1 filename: ");
        seq1FileName = getFileNameFromUserInput();
        System.out.print("Sequence 2 filename: ");
        seq2FileName = getFileNameFromUserInput();

        System.out.print("Enter the match score: ");
        matchscore = scanner.nextInt();

        System.out.print("Enter the mismatch score: ");
        mismatchscore = scanner.nextInt();

        System.out.print("Enter the gap score: ");
        gapscore = scanner.nextInt();

        System.out.print("Print scoring matrix and path to screen? (Y/N): ");
        print = scanner.next().charAt(0);
        
        seq1File = new File("Ch3/" + seq1FileName);
        seq2File = new File("Ch3/" + seq2FileName);
        seq1 = getSequenceFromFile(seq1File);
        seq2 = getSequenceFromFile(seq2File);

        // adding extra characters so sequences align with matrix
        // saves some calculations later on
        seq1 = " " + seq1;
        seq2 = " " + seq2;

        len1 = seq1.length();
        len2 = seq2.length();

        matrix = new int[seq1.length() + 1][seq2.length() + 1];

        // initialize 1st row and 1st column of matrix
        matrix[0][0] = 0;
        for (int i = 1; i < seq1.length(); i ++){
            matrix[i][0] = 0;
        }
        for (int i = 1; i < seq2.length(); i ++){
            matrix[0][i] = 0;
        }

        // STEP 1:
        // Fill in rest of matrix using the following rules:
        // determine three possible values for matrix[x][y]
        // value 1 = add gap score to matrix[x][y-1]
        // value 2 = add gap score to matrix[x-1][y]
        // value 3 = add match score or mismatch score to 
        //           matrix[x-1][y-1] depending on nucleotide 
        //           match for position x of seq1 and position y
        //           of seq2
        // place the largest of the three values in matrix[x][y]
        //
        // Best alignment score appears in matrix[len1][len2].
        for(int x = 1; x < len1; x++){
           for(int y = 1; y < len2; y++){
            val1 = matrix[x][y-1] + gapscore;
            val2 = matrix[x-1][y] + gapscore;
            if (seq1.charAt(x) == (seq2.charAt(y))) {
                   val3 = matrix[x-1][y-1] + matchscore;
            }
            else{
               val3 = matrix[x-1][y-1] + mismatchscore;
            }
            if ((val1 >= val2) && (val1 >= val3)){
               matrix[x][y] = val1;
            }
            else if ((val2 >= val1) && (val2 >= val3)){
               matrix[x][y] = val2;
            }
            else{
               matrix[x][y] = val3;
            }

            if (matrix[x][y] < 0) {
              matrix[x][y] = 0;
            }
           }
        }

        if ((print == 'Y') || (print == 'y')){
           // Display scoring matrix
           System.out.println("MATRIX:");
           for(int x = 0; x < len1; x++){
              for(int y = 0; y < len2; y++){
                System.out.print(matrix[x][y] + " ");
              }
              System.out.print("\n");
           }
           System.out.println("\n");
        }
        // STEP 2:
        // Begin at matrix[len1][len2] and find a path to 
        // matrix[0][0].
        // Build string to hold path pattern by concatenating either 
        // 'H' (current cell produced by cell horizontally to left), 
        // 'D' (current cell produced by cell on diagonal), 
        // 'V' (current cell produced by cell vertically above)
        // It is possible for cell to have been produced by more
        // than one of these possible choices, if multiple optimal 
        // alignments exist. For now, choose only one.
        // 
        pathrow = len1-1;
        pathcol = len2-1;
        while ((pathrow != 0) || (pathcol != 0)){
            if (pathrow == 0){
               // must be from cell to left
               path = path + 'H';
               pathcol = pathcol - 1;
            }
            else if (pathcol == 0){
               // must be from cell above
               path = path + 'V';
               pathrow = pathrow - 1;
            }
            // could be from any direction
            else if  ((matrix[pathrow][pathcol-1] + gapscore) == matrix[pathrow][pathcol]){
               // from left
               path = path + 'H';
               pathcol = pathcol - 1;
            }
            else if ((matrix[pathrow-1][pathcol] + gapscore) == matrix[pathrow][pathcol]){
               // from above
               path = path + 'V';
               pathrow = pathrow - 1;
            } 
            else{
               // must be from diagonal
               path = path + 'D';
               pathrow = pathrow - 1;
               pathcol = pathcol - 1;
            }
        }   

        // following added to display forward path:
        forwardpath = new StringBuilder(path).reverse().toString();

        if ((print == 'Y') || (print == 'y')){
            System.out.println("Backtrack path is " + path);
            System.out.println("Forwards path is " + forwardpath);
        }

        // STEP 3:
        // Determine alignment pattern by reading path string 
        // created in step 2.
        // Create two string variables (alignseq1 and alignseq2) to hold 
        // the sequences for alignment.
        // Reading backwards from seq1, seq2 and path string, 
        //   if string character is 'D', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, last char in seq2
        //   if string character is 'V', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, the gap char
        //   if string character is 'H', then 
        //      concatenate to front of alignseq1 the gap char
        //      and to the front of alignseq2, last char in seq2
        // Continue process until path string has been traversed.
        // Result appears in string alignseq1 and seq2
        //

        seq1loc = len1-1;
        seq2loc = len2-1;
        pathloc = 0;
        while (pathloc < path.length()){
           if (path.charAt(pathloc) == ('D')){
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq1loc--;
              seq2loc--;
           }
           else if (path.charAt(pathloc) == ('V')) {
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = '-' + alignseq2;
              seq1loc--;
           }  
           else{  // must be an H
              alignseq1 = '-' + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq2loc--;
           }  
           pathloc++;
        }

        System.out.println("\nAligned Sequences:\n");
        for (int i=0; i<alignseq1.length(); i++){
            System.out.print(alignseq1.charAt(i));
        }

        System.out.println("");

        for (int i=0; i<alignseq1.length(); i++){
            System.out.print(alignseq2.charAt(i));
        }

        System.out.println("\n");

        // print whatever is left
        // System.out.println(alignseq1.charAt(i));
        // System.out.println(alignseq2.charAt(i));

        int count = 0;
        for (int i=0; i < alignseq1.length(); i++){
           if (alignseq1.charAt(i) == (alignseq2.charAt(i))){
              count++;
           }
        }
        System.out.println("Total number of aligned characters: " + count);
        System.out.println("Percentage match: " + ((count * 1.0) / alignseq1.length()) + "\n\n");
    }

    public static void affinityAlignment() {
        String seq1FileName, seq2FileName;
        File seq1File, seq2File;
        String seq1, seq2, path = "", alignseq1 = "", alignseq2 = "", forwardpath = "";
        Scanner scanner = new Scanner(System.in); 
        int gapstart, gapextend, matchscore, mismatchscore, len1, len2, val1, val2, val3, pathrow, pathcol, seq1loc, seq2loc, pathloc;
        int hval1, hval2, vval1, vval2, vval, neginf, lastD = 0, lastH = 0, lastV = 0, s, x, y, i, j, count;
        char print1, layer=' ';
        int[][] Dmatrix, Hmatrix, Vmatrix;

        
        System.out.print("Sequence 1 filename: ");
        seq1FileName = getFileNameFromUserInput();
        System.out.print("Sequence 2 filename: ");
        seq2FileName = getFileNameFromUserInput();

        System.out.print("Enter the match score: ");
        matchscore = scanner.nextInt();

        System.out.print("Enter the mismatch score: ");
        mismatchscore = scanner.nextInt();

        System.out.print("Enter the gap start: ");
        gapstart = scanner.nextInt();

        System.out.print("Enter the gap extension: ");
        gapextend = scanner.nextInt();

        System.out.print("Print scoring matrix and path to screen? (Y/N): ");
        print1 = scanner.next().charAt(0);
        
        seq1File = new File("Ch3/" + seq1FileName);
        seq2File = new File("Ch3/" + seq2FileName);

        seq1 = getSequenceFromFile(seq1File);
        seq2 = getSequenceFromFile(seq2File);

        // adding extra characters so sequences align with matrix
        // saves some calculations later on
        seq1 = " " + seq1;
        seq2 = " " + seq2;

        len1 = seq1.length();
        len2 = seq2.length();

        Dmatrix = new int[seq1.length()][seq2.length()];
        Hmatrix = new int[seq1.length()][seq2.length()];
        Vmatrix = new int[seq1.length()][seq2.length()];

        neginf = (len1+len2)*gapstart;

        // declare all two dimensional arrays and initialize to spaces
        // arrays must contain one extra row and one extra column
        for(i = 0; i < len1; i++){
           for(j = 0; j < len2; j++){
              Dmatrix[i][j] = neginf;
              Hmatrix[i][j] = neginf;
              Vmatrix[i][j] = neginf;
           }
        }

        // now correctly initialize Origin of OPT, first Column of Vmatrix and first Row of Hmatrix
        Dmatrix[0][0] = 0;
        Hmatrix[0][1] = gapstart;
        Vmatrix[1][0] = gapstart;

        // consider 0th column only
        for (i = 2; i < len1; i ++){
            //Vmatrix[i][0] = Vmatrix[i-1][0] + gapextend;
            Vmatrix[i][0] = gapstart + (i-1)*gapextend;
        //    print "\nDEBUG: Vmatrix[i][0] is VMatrix[i][0]\n";
        }


        // consider 0th row only
        for (i = 2; i < len2; i ++){
            //Hmatrix[0][i] = Hmatrix[0][i-1] + gapextend;
            Hmatrix[0][i] = gapstart + (i-1)*gapextend;
        //    print "\nDEBUG: Hmatrix[0][i] is HMatrix[0][i]\n";
        }


        // Fill in rest of matrices using the following rules:
        // 
        // determine three possible values for Dmatrix[x][y]
        // value 1 = Hmatrix[x-1][y-1]+sub-or-match score
        // value 2 = Vmatrix[x-1][y-1]+sub-or-match score
        // value 3 = Dmatrix[x-1][y-1]+sub-or-match score
        //
        // determine two possible values for Hmatrix[x][y]
        // hvalue 1 = Hmatrix[x][y-1] + gapextend
        // hvalue 2 = Dmatrix[x][y-1] + gapstart
        //
        // determine two possible values for Vmatrix[x][y]
        // vvalue 1 = Vmatrix[x-1][y] + gapextend
        // vvalue 2 = Dmatrix[x-1][y] + gapstart
        //
        // place the largest of the respective values into the respective matrix
        //
        // Best alignment score appears in matrix[len1][len2].

        for(x = 1; x < len1; x++){
           for(y = 1; y < len2; y++){
               // first calculate the possible substitution score
               if (seq1.charAt(x) == seq2.charAt(y)){
                s = matchscore;
               }
               else { s = mismatchscore; }
               
               // calculate possible horizontal value
               if (Hmatrix[x][y-1] == neginf) { hval1 = neginf; }
               else { hval1 = Hmatrix[x][y-1] + gapextend; }
               
               if (Dmatrix[x][y-1] == neginf) { hval2 = neginf; }
               else { hval2 = Dmatrix[x][y-1] + gapstart; }

               // calculate possible vertical value
               if (Vmatrix[x-1][y] == neginf) { vval1 = neginf; }
               else { vval1 = Vmatrix[x-1][y] + gapextend; }
               
               if (Dmatrix[x-1][y] == neginf) { vval2 = neginf; }
               else { vval2 = Dmatrix[x-1][y] + gapstart; }
               
               // calculate possible middle value
               if (Hmatrix[x-1][y-1] == neginf) { val1 = neginf; }
               else { val1 = Hmatrix[x-1][y-1] + s; }

               if (Vmatrix[x-1][y-1] == neginf) { val2 = neginf; }
               else { val2 = Vmatrix[x-1][y-1] + s; }
               
               if (Dmatrix[x-1][y-1] == neginf) {
               val3 = neginf;
               }
               else {
               val3 = Dmatrix[x-1][y-1] + s;
               }
               
               // now to determine max for OPT matrix
               if ((val1 >= val2) && (val1 >= val3)){
               Dmatrix[x][y] = val1;
               }
               else if ((val2 >= val1) && (val2 >= val3)){
               Dmatrix[x][y] = val2;
               }
               else{
               Dmatrix[x][y] = val3;
               }

               // now to determine max for Hmatrix
               if (hval1 >= hval2) {  Hmatrix[x][y] = hval1; }
               else { Hmatrix[x][y] = hval2; }
               
               // and finally determine max for Vmatrix
               if (vval1 >= vval2) { Vmatrix[x][y] = vval1; }
               else { Vmatrix[x][y] = vval2; }
               
           }
        }


        if ((print1 == 'Y') || (print1 == 'y')){
            // Display scoring matrices
            System.out.println("Middle MATRIX:");
            for(x = 0; x < len1; x++){
                for(y = 0; y < len2; y++){
                    lastD = Dmatrix[x][y];
                    System.out.print(Dmatrix[x][y] + " ");
                }
                System.out.println("");
            }
            System.out.println("\n");
            System.out.println("Horizontal MATRIX:");
            for(x = 0; x < len1; x++){
            for(y = 0; y < len2; y++){
                lastH = Hmatrix[x][y];
                System.out.print(Hmatrix[x][y] + " ");
            }
            System.out.println("");
            }
            System.out.println("");
            System.out.println("Vertical MATRIX:");
            for(x = 0; x < len1; x++){
            for(y = 0; y < len2; y++){
                lastV = Vmatrix[x][y];
                System.out.print(Vmatrix[x][y] + " ");
            }
            System.out.println("");
            }
            System.out.println("");
        }

        // Now to find the alignment path
        // Begin at the MAX valued last matrix entry over the three matrices!
        // This is a modification: That value is NOT necessarily at the DMatrix!
        // Find a path from max_{M=D,H,V} (Mmatrix[n][m]) entry BACK TO the origin [0,0] 
        // Build string to hold path pattern by concatenating either 
        // 'H' (current cell produced by cell horizontally to left), 
        // 'D' (current cell produced by cell on diagonal), 
        // 'V' (current cell produced by cell vertically above)
        // It is possible for cell to have been produced by more
        // than one of these possible choices, if multiple optimal 
        // alignments exist. For now, choose only one.
        // 
        pathrow = len1-1;
        pathcol = len2-1;

        //layer = 'D'; // this helps keep track of matrices, and is initialized to the middle layer BUT should be initialized to whatever layer has maximum last entry!
        // MODIFICATION:
        if ((lastV >= lastH) && (lastV >= lastD)) {layer = 'V';}
        if ((lastD >= lastH) && (lastD >= lastV)) {layer = 'D';}
        if ((lastH >= lastV) && (lastH >= lastD)) {layer = 'H';}
        System.out.println("Initialized backtrack layer " + layer); 
        //////


        while ((pathrow != 0) || (pathcol != 0)){
            // first calculate s the substitution score again
            if (seq1.charAt(pathrow) == seq2.charAt(pathcol)){
            s = matchscore;
            }
            else { s = mismatchscore; }
            
            if (pathrow == 0){
               // When you're on the 0th row:
               // must be via a horizontal gap
               path = path + 'H';
               pathcol = pathcol - 1;
               layer = 'H';
            }
            else if (pathcol == 0){
               // When you're on the 0th column:
               // must be via a vertical gap
               path = path + 'V';
               pathrow = pathrow - 1;
               layer = 'V';
            }

            // Otherwise more cases to consider:
            else if (layer == 'D') { // if we were in the middle layer
            // must be because of a substitution move
            // now need to possibly update the layer if previous move was from H or V
            if (Dmatrix[pathrow][pathcol] == Hmatrix[pathrow-1][pathcol-1] + s) {
                layer = 'H'; // need to move to the horizontal matrix
            }
            else if (Dmatrix[pathrow][pathcol] == Vmatrix[pathrow-1][pathcol-1] + s) {
                layer = 'V'; // need to move to the vertical matrix
            } // else don't change the layer anyway
            // do the following updates in any case:
            path = path + 'D';
            pathrow = pathrow - 1;
            pathcol = pathcol - 1;
            }

            else if (layer == 'H') { // if we were in the horizontal layer
            // but need to switch to the middle layer if this was a gap start move
            if (Hmatrix[pathrow][pathcol] == Dmatrix[pathrow][pathcol-1] + gapstart) {
                layer = 'D';
            }
            // do the following updates in any case
            path = path + 'H'; // current move is a gap
            pathcol = pathcol - 1;
            }

            else if (layer == 'V') { // only remaining situation is that we're in the vertical layer
            // but need to switch to the middle layer if this was a gap start move
            if (Vmatrix[pathrow][pathcol] == Dmatrix[pathrow-1][pathcol] + gapstart) {
                layer = 'D';
            }
            // do the following updates in any case
            path = path + 'V'; // current move is a gap
            pathrow = pathrow - 1;
            }  
        } // this ends the while 



        // following added to display forward path:
        forwardpath = new StringBuilder(path).reverse().toString();
        if ((print1 == 'Y') || (print1 == 'y')){
            System.out.println("Backtrack path is " + path);
            System.out.println("Forwards path is " + forwardpath);
        }

        // STEP 3:
        // Determine alignment pattern by reading path string 
        // created in step 2.
        // Create two string variables (alignseq1 and alignseq2) to hold 
        // the sequences for alignment.
        // Reading backwards from seq1, seq2 and path string, 
        //   if string character is 'D', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, last char in seq2
        //   if string character is 'V', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, the gap char
        //   if string character is 'H', then 
        //      concatenate to front of alignseq1 the gap char
        //      and to the front of alignseq2, last char in seq2
        // Continue process until path string has been traversed.
        // Result appears in string alignseq1 and seq2
        //

        seq1loc = len1-1;
        seq2loc = len2-1;
        pathloc = 0;

        // initializing aligned strings
        alignseq1 = "";
        alignseq2 = "";

        while (pathloc < path.length()){
           if (path.charAt(pathloc) == 'D'){
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq1loc--;
              seq2loc--;
           }
           else if (path.charAt(pathloc) == 'V'){
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = '-' + alignseq2;
              seq1loc--;
           }  
           else if(path.charAt(pathloc) == 'H') {  // must be an H
              alignseq1 = '-' + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq2loc--;
           }  
           pathloc++;
        }

        System.out.println("\nAligned Sequences:");
        for (i=0; i<alignseq1.length(); i++){
            System.out.print(alignseq1.charAt(i));
        }

        System.out.println("");

        for (i=0; i<alignseq2.length(); i++){
            System.out.print(alignseq2.charAt(i));
        }
        System.out.println("\n");

        count = 0;
        for (i=0; i < alignseq1.length(); i++){
           if (alignseq1.charAt(i) == alignseq2.charAt(i)){
              count++;
           }
        }
        
        System.out.println("Total number of aligned characters: " + count);
        System.out.println("Percentage match: " + ((count * 1.0) / alignseq1.length()));
    }

    public static void globalAlignment() {
        String seq1FileName, seq2FileName;
        File seq1File, seq2File;
        String seq1, seq2, path = "", alignseq1 = "", alignseq2 = "", forwardpath = "";
        Scanner scanner = new Scanner(System.in); 
        int gapscore, matchscore, mismatchscore, len1, len2, val1, val2, val3, pathrow, pathcol, seq1loc, seq2loc, pathloc;
        char print;
        int[][] matrix;
        
        System.out.print("Sequence 1 filename: ");
        seq1FileName = getFileNameFromUserInput();
        System.out.print("Sequence 2 filename: ");
        seq2FileName = getFileNameFromUserInput();

        System.out.print("Enter the match score: ");
        matchscore = scanner.nextInt();

        System.out.print("Enter the mismatch score: ");
        mismatchscore = scanner.nextInt();

        System.out.print("Enter the gap score: ");
        gapscore = scanner.nextInt();

        System.out.print("Print scoring matrix and path to screen? (Y/N): ");
        print = scanner.next().charAt(0);
        
        seq1File = new File("Ch3/" + seq1FileName);
        seq2File = new File("Ch3/" + seq2FileName);

        seq1 = getSequenceFromFile(seq1File);
        seq2 = getSequenceFromFile(seq2File);

        // adding extra characters so sequences align with matrix
        // saves some calculations later on
        seq1 = " " + seq1;
        seq2 = " " + seq2;

        len1 = seq1.length();
        len2 = seq2.length();

        matrix = new int[seq1.length() + 1][seq2.length() + 1];

        // initialize 1st row and 1st column of matrix
        matrix[0][0] = 0;
        for (int i = 1; i < seq1.length(); i ++){
            matrix[i][0] = matrix[i-1][0] + gapscore;
        }
        for (int i = 1; i < seq2.length(); i ++){
            matrix[0][i] = matrix[0][i-1] + gapscore;
        }

        // STEP 1:
        // Fill in rest of matrix using the following rules:
        // determine three possible values for matrix[x][y]
        // value 1 = add gap score to matrix[x][y-1]
        // value 2 = add gap score to matrix[x-1][y]
        // value 3 = add match score or mismatch score to 
        //           matrix[x-1][y-1] depending on nucleotide 
        //           match for position x of seq1 and position y
        //           of seq2
        // place the largest of the three values in matrix[x][y]
        //
        // Best alignment score appears in matrix[len1][len2].
        for(int x = 1; x < len1; x++){
           for(int y = 1; y < len2; y++){
            val1 = matrix[x][y-1] + gapscore;
            val2 = matrix[x-1][y] + gapscore;
            if (seq1.charAt(x) == (seq2.charAt(y))) {
                   val3 = matrix[x-1][y-1] + matchscore;
            }
            else{
               val3 = matrix[x-1][y-1] + mismatchscore;
            }
            if ((val1 >= val2) && (val1 >= val3)){
               matrix[x][y] = val1;
            }
            else if ((val2 >= val1) && (val2 >= val3)){
               matrix[x][y] = val2;
            }
            else{
               matrix[x][y] = val3;
            }
           }
        }

        if ((print == 'Y') || (print == 'y')){
           // Display scoring matrix
           System.out.println("MATRIX:");
           for(int x = 0; x < len1; x++){
              for(int y = 0; y < len2; y++){
                System.out.print(matrix[x][y] + " ");
              }
              System.out.print("\n");
           }
           System.out.println("\n");
        }
        // STEP 2:
        // Begin at matrix[len1][len2] and find a path to 
        // matrix[0][0].
        // Build string to hold path pattern by concatenating either 
        // 'H' (current cell produced by cell horizontally to left), 
        // 'D' (current cell produced by cell on diagonal), 
        // 'V' (current cell produced by cell vertically above)
        // It is possible for cell to have been produced by more
        // than one of these possible choices, if multiple optimal 
        // alignments exist. For now, choose only one.
        // 
        pathrow = len1-1;
        pathcol = len2-1;
        while ((pathrow != 0) || (pathcol != 0)){
            if (pathrow == 0){
               // must be from cell to left
               path = path + 'H';
               pathcol = pathcol - 1;
            }
            else if (pathcol == 0){
               // must be from cell above
               path = path + 'V';
               pathrow = pathrow - 1;
            }
            // could be from any direction
            else if  ((matrix[pathrow][pathcol-1] + gapscore) == matrix[pathrow][pathcol]){
               // from left
               path = path + 'H';
               pathcol = pathcol - 1;
            }
            else if ((matrix[pathrow-1][pathcol] + gapscore) == matrix[pathrow][pathcol]){
               // from above
               path = path + 'V';
               pathrow = pathrow - 1;
            } 
            else{
               // must be from diagonal
               path = path + 'D';
               pathrow = pathrow - 1;
               pathcol = pathcol - 1;
            }
        }   

        // following added to display forward path:
        forwardpath = new StringBuilder(path).reverse().toString();

        if ((print == 'Y') || (print == 'y')){
            System.out.println("Backtrack path is " + path);
            System.out.println("Forwards path is " + forwardpath);
        }

        // STEP 3:
        // Determine alignment pattern by reading path string 
        // created in step 2.
        // Create two string variables (alignseq1 and alignseq2) to hold 
        // the sequences for alignment.
        // Reading backwards from seq1, seq2 and path string, 
        //   if string character is 'D', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, last char in seq2
        //   if string character is 'V', then 
        //      concatenate to front of alignseq1, last char in seq1 
        //      and to the front of alignseq2, the gap char
        //   if string character is 'H', then 
        //      concatenate to front of alignseq1 the gap char
        //      and to the front of alignseq2, last char in seq2
        // Continue process until path string has been traversed.
        // Result appears in string alignseq1 and seq2
        //

        seq1loc = len1-1;
        seq2loc = len2-1;
        pathloc = 0;
        while (pathloc < path.length()){
           if (path.charAt(pathloc) == ('D')){
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq1loc--;
              seq2loc--;
           }
           else if (path.charAt(pathloc) == ('V')) {
              alignseq1 = seq1.charAt(seq1loc) + alignseq1;
              alignseq2 = '-' + alignseq2;
              seq1loc--;
           }  
           else{  // must be an H
              alignseq1 = '-' + alignseq1;
              alignseq2 = seq2.charAt(seq2loc) + alignseq2;
              seq2loc--;
           }  
           pathloc++;
        }

        System.out.println("\nAligned Sequences:\n");
        for (int i=0; i<alignseq1.length(); i++){
            System.out.print(alignseq1.charAt(i));
        }

        System.out.println("");

        for (int i=0; i<alignseq1.length(); i++){
            System.out.print(alignseq2.charAt(i));
        }

        System.out.println("\n");

        // print whatever is left
        // System.out.println(alignseq1.charAt(i));
        // System.out.println(alignseq2.charAt(i));

        int count = 0;
        for (int i=0; i < alignseq1.length(); i++){
           if (alignseq1.charAt(i) == (alignseq2.charAt(i))){
              count++;
           }
        }
        System.out.println("Total number of aligned characters: " + count);
        System.out.println("Percentage match: " + ((count * 1.0) / alignseq1.length()) + "\n\n");
    }
}
