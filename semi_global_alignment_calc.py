import numpy as np
from cell import Cell
import sys
from Bio import SeqIO

# Score scheme contants 
gap_score = -1
match_score = 1
mismatch_score = -1

# Get the score of a cell in the matrix
def get_score(mat, i, j):
    return mat[i][j].score

# Recursion rule for evaluating a new cell score
def eval_score(mat, i, j, seq1, seq2):
  S_ij = 1 if seq1[i] == seq2[j] else -1

  down_val = get_score(mat, i-1, j) + gap_score
  right_val = get_score(mat, i, j-1) + gap_score
  diag_val = get_score(mat, i-1, j-1) + S_ij

  val_list = [down_val, right_val, diag_val] 
  array = np.array(val_list)
  max_index = array.argmax()
  return max_index, val_list[max_index]

# Print the score matrix for debug
def print_score_mat(mat, seq1, seq2):
  print('  ', end = ' ')
  for i in range(len(seq2)):
    print(seq2[i], end = '  ')
  print("")

  for i in range(len(seq1)):
    print(seq1[i], end = ' ')
    for j in range(len(seq2)):
      print("{:02d}".format(mat[i][j].score), end = ' ')
    print("")
    
# Print the direction matrix for debug
def print_direction_mat(mat, seq1, seq2):
  print(' ', end = ' ')
  for i in range(len(seq2)):
    print(seq2[i], end = ' ')
  print("")

  for i in range(len(seq1)):
    print(seq1[i], end = ' ')
    for j in range(len(seq2)):
      print(mat[i][j].direction, end = ' ')
    print("")

# Find the optimal score for the alignment
def find_opt_score(mat, seq1, seq2): 
  last_row = [cell.score for cell in mat[len(seq1)-1]]
  arr = np.array(last_row)
  max_index = arr.argmax()
  return (len(seq1)-1 ,max_index), last_row[max_index]

# compute where the last cell is on the route 
def step_trace_back(mat, i, j):
  if mat[i][j].direction == 0:
    return i-1, j
  elif mat[i][j].direction == 1:
    return i, j-1
  elif mat[i][j].direction == 2:
    return i-1, j-1

# Compute the route for the alignment 
def trace_route(mat, seq1, seq2) :
  cur_x = len(seq1)-1
  cur_y = find_opt_score(mat, seq1, seq2)[0][1]

  route = []

  while cur_x != 0:
    route.insert(0, (cur_x, cur_y))
    cur_x, cur_y = step_trace_back(mat, cur_x, cur_y)
     
  return route

# Print the sub sequences fpr 
def print_sub_sequence(mat, route, seq1, seq2):
  seq1_align = []
  seq2_align = []
  for i in range(len(route)):
    if mat[route[i][0]][route[i][1]].direction == 0:
      seq1_align.append(seq1[route[i][0]])
      seq2_align.append("-")
    elif mat[route[i][0]][route[i][1]].direction == 1:
      seq1_align.append("-")
      seq2_align.append(seq2[route[i][1]])
    elif mat[route[i][0]][route[i][1]].direction == 2:
      seq1_align.append(seq1[route[i][0]])
      seq2_align.append(seq2[route[i][1]])

  return seq1_align, seq2_align

def fit_align(sequence1, sequence2):
  seq1 = []
  seq2 = []
  seq1[:0] = "-" + sequence1
  seq2[:0] = "-" + sequence2

  # Creating the alignment matrix 
  mat = []

  for i in range(len(seq1)):
    row = []
    for j in range(len(seq2)):
      row.append(Cell(99, None))
    mat.append(row)

  # initialize first row to be 0 and right direction
  for i in range(len(seq2)):
    mat[0][i].score = 0
    mat[0][i].direction = 1

  # initialize first column to decrement and down direction
  for i in range(len(seq1)):
    mat[i][0].score = -i
    mat[i][0].direction = 0


  for i in range(1, len(seq1)):
    for j in range(1, len(seq2)):
      mat[i][j].score = eval_score(mat, i, j, seq1, seq2)[1]
      mat[i][j].direction = eval_score(mat, i, j, seq1, seq2)[0]



  max_score = find_opt_score(mat, seq1, seq2)[1]

  best_sub_seq1 = print_sub_sequence(mat, trace_route(mat, seq1, seq2), seq1, seq2)[0]
  best_sub_seq2 = print_sub_sequence(mat, trace_route(mat, seq1, seq2), seq1, seq2)[1]

  print_score_mat(mat, seq1, seq2)

  #print_direction_mat(mat, seq1, seq2)
  
  #print(trace_route(mat, seq1, seq2))

  return max_score, best_sub_seq1, best_sub_seq2


def main():
  # parse the input 
  records = SeqIO.parse(sys.argv[1], "fasta")
  sequences = [r.seq for r in records]

  # Write to the output file 
  f = open(sys.argv[2], "w")
  max_score, best_sub_seq1, best_sub_seq2 = fit_align(str(sequences[0]), str(sequences[1]))
  f.write(str(max_score))
  f.write("\n")
  f.write(''.join(map(str, best_sub_seq1)))
  f.write("\n")
  f.write(''.join(map(str, best_sub_seq2)))
  f.write("\n")

if __name__ == "__main__":
    main()






