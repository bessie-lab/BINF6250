# BINF6250 - Project 06

# Introduction


# Pseudocode
def read_fasta(filename: str) -> Dict[str, str]:
    sequences = empty dictionary (name to sequence)
    current_name = nothing yet

    Open the file and read it line by line:
        strip spaces and the newline off the end of the line
        IF the line is empty:
            #skip it
        IF the line starts with ">":
            #current_name = last word of the line (i.e. "DM1")
            #sequences[current_name] = empty string
        ELSE (it's part of a sequence):
            add the line to the end of sequences[current_name]

    Return sequences

Changes from Project 05 smith_waterman:
    1. cal_score is kept as-is; it still calculates each cell's score. The traceback arrow it records is not used here.
    2. Drop traceback, along with the traceback_matrix lines and the max_pos line, since nothing uses them anymore.
    3. Keep the original table-filling loop and the best_score line in a helper function, and remove both print statements.
    4. smith_waterman uses the Project 06 name and inputs, calls the helper, and returns the normalized score:
       normalized = raw(seq1, seq2) / smaller of [raw(seq1, seq1), raw(seq2, seq2)]

def cal_score(matrix, seq1, seq2, i, j, match, mismatch, gap):  #Original calc_score from Project05 untouched.  Needed to run Smith_Waterman function from Project05

def project5_smith_waterman_score(seq1, seq2, match, mismatch, gap):  #Original smith_waterman function from Project05, turned into a helper function for the new smith_waterman function in Project06

def smith_waterman(
    sequence1: str,
    sequence2: str,
    match_score: float = 2.0,
    mismatch_penalty: float = -1.0,
    gap_penalty: float = -1.0
) -> float:
    raw_score = project5_smith_waterman_score(sequence1 vs sequence2, using match_score, mismatch_penalty, gap_penalty)
    self_score1 = project5_smith_waterman_score(sequence1 vs sequence1, same scores)
    self_score2 = project5_smith_waterman_score(sequence2 vs sequence2, same scores)

    smaller_self_score = smaller of self_score1 and self_score2

    IF smaller_self_score is 0 (only happens if a sequence is empty):
        return 0.0 (avoids dividing by zero)
    
    Return (raw_score / smaller_self_score) #Note: 1.0 = identical, closer to 0 = less similar

def build_distance_matrix(sequences: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
    ids = list of the sequence names, in file order (i.e. ["DM1", "DM2", ...])
    n = number of sequences
    distance_matrix = n x n grid of zeros (diagonal stays 0: each sequence vs itself)

    For each row i from 0 to n-1:
        FOR each column j from i+1 to n-1 (only pairs not yet done; skip the diagonal):
            similarity = smith_waterman (sequence for ids[i], sequence for ids[j])
            distance = 1 - similarity ***my suggestion for the group because can't produce undefined values and for those highly similar sequences it gives nearly the same distances as -log(similarity) 
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  (Table is symmetric, so fill both halves)

    Return distance_matrix, ids

def neighbor_joining(
    distance_matrix: np.ndarray,
    labels: List[str]
) -> str:
    PROPOSAL FOR DECISIONS IN PSEUDOCODE/CODE
        Final branch between the last two nodes is split evenly between them
        Branch lengths are written to 4 decimals places in the Newick text (calculations use full precision)
        All 20 sequences are used.  The docstring expects 19, note it in README
        Ties in the D* score: the first pair found is kept
        Negative branch lengths: set to 0, differences added to the other branch
    
    SETUP
    D = copy of distance_matrix (so original table isn't changed)
    current_labels = copy of labels (each starts as just a name, i.e. "DM1")

    WHILE more than 2 nodes remain:
        n = number of rows in D

        Step 1. Total Distance:  totals = sum of each for row of D

        Step 2. D* score for every pair, keep the lowest:  lowest_score = a very large number
        FOR each row i from 0 to n-1:
            FOR each column j from i+1 to n-1:
                score = (n-2) * D[i, j] - totals[i] - totals[j]
                IF score < lowest_score:
                    lowest_score = score
                    best_i = i
                    best_j = j

        Step 3. Branch lengths for best_i and best_j
        branch_i = (D[best_i, best_j] + (totals[best_i] - totals[best_j]) / (n-2)) / 2
        branch_j = D[best_i, best_j] - branch_i
        IF either branch is negative:
            set it to 0 and add the difference to the other branch
        
        Step 4. Distances from the new node to every other node k
        new_distances = empty list
        FOR each k that is not best_i or best_j:
            new_distances.append((D[best_i, k] + D[best_j, k] - D[best_i, best_j]) / 2)

        New_label for the joined pair:
        new_label = "(" + current_labels[best_i] + ":" + f"{branch_i:.4f}" + "," + current_labels[best_j] + ":" + f"{branch_j:.4f}" + ")"

        Step 5. Shrink the table
        keep = every index except best_i and best_j
        new_D = (n-1) x (n-1) grid of zeros
        copy D's row/columns for the kept indices into the top left of new_D
        put new_distances in the last row and last column of new_D (new node goes last)
        current_labels = labels at the kept indices, then new_label added at the end
        D = new_D

    Step 6. Exactly 2 nodes left: connect them  (The final branch is split evenly between two sides so the tree can be drawn)
    final_distance = D[0, 1] / 2
    tree = "(" + current_labels[0] + ":" + f"{final_distance:.4f}" + "," + current_labels[1] + ":" + f"{final_distance:.4f}" + ");"

    Return tree    

# Results Output

# Successes

# Struggles

# Reflection and Adaptive Learning

# Personal Reflections
## Group leader 

## Other member

## Other member

# Generative AI Appendix
