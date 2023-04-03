# semi-global-alignment-calculator
Dependency used: 
Bio 

How to run the programs:
Part1:
python3 semi_global_alignment_calc.py input.fasta output.txt
The first argument is the input fasta file and the second argument is the file to write the output to.

Part2: 
python3 best_align.py query.fasta library.fasta output1.txt
The first argument is the input file for the query sequence and the second argument is the input file for the library sequences and the third argument is the file to write output to. 

Scoring metrics: 
We use the following score scheme and linear gap penalty: match = 1, mismatch = -1, indel=-1

Since we are trying to find a substring T’ of T such that the global sequence alignment score between S and T’ is maximized. We can achieve this by performing some modifications to the pairwise sequence alignment matrix. Since we are free to take the substring of T but not S, we should not penalize any gaps at the beginning of T. Hence, for our matrix, the rows will be sequence S and the columns will be T. The first rows of the matrix will all be zeros and the first column will be decrementing by one each row from 0. 

The recurrence relationship will be the same as the global alignment. 

To get the optimal score, we look at the last row and find the max value in the last row if the length of sequence T is longer than S. Otherwise, we find the last score in the len(T)-th row.  Then, we perform the same backtrack method from that optimal score and the termination condition is when we reach back at the first row

